import functools
import io
import sys
import threading
import time
from collections.abc import Callable, Collection, Generator, Iterable, Sized
from collections.abc import Sequence as PySeq
from typing import Any, ParamSpec, Self, TypeVar

from tqdm import notebook, tqdm

from cogent3.util import parallel as PAR
from cogent3.util.misc import in_jupyter

P = ParamSpec("P")
R = TypeVar("R")
T = TypeVar("T")


class LogFileOutput:
    """A fake progress bar for when progress bars are impossible"""

    def __init__(self, **kw: Any) -> None:
        self.n: float = 0
        self.message = ""
        self.t0 = time.time()
        self.lpad = ""
        self.output = sys.stdout  # sys.stderr

    def set_description(self, desc: str = "", refresh: bool = False) -> None:
        self.message = desc

    def close(self) -> None:
        pass

    def refresh(self) -> None:
        if self.message:
            delta = f"+{int(time.time() - self.t0)}"
            progress = int(100 * self.n + 0.5)
            print(
                "%s %5s %3i%% %s" % (self.lpad, delta, progress, str(self.message)),
                file=self.output,
            )

    @classmethod
    def write(cls, *args: Any, **kwargs: Any) -> None:
        pass


class ProgressContext:
    def __init__(
        self,
        progress_bar_type: type[tqdm | notebook.tqdm | LogFileOutput] | None = None,
        depth: int = -1,
        message: str | None = None,
        mininterval: float = 1.0,
    ) -> None:
        self.progress_bar_type = progress_bar_type
        self.progress_bar: tqdm[Any] | notebook.tqdm[Any] | LogFileOutput | None = None
        self.progress: float = 0
        self.depth = depth
        self.message = message
        self.mininterval = mininterval

    def set_new_progress_bar(self) -> None:
        if self.progress_bar_type:
            self.progress_bar = self.progress_bar_type(
                total=1,
                position=self.depth,
                leave=True,
                bar_format="{desc} {percentage:3.0f}%|{bar}|{elapsed}<{remaining}",
                mininterval=self.mininterval,
                dynamic_ncols=True,
            )

    def subcontext(self, *args: Any, **kw: Any):
        return ProgressContext(
            progress_bar_type=self.progress_bar_type,
            depth=self.depth + 1,
            message=self.message,
            mininterval=self.mininterval,
        )

    def display(self, msg: str | None = None, progress: float | None = None) -> None:
        if not self.progress_bar:
            self.set_new_progress_bar()
        updated = False
        assert self.progress_bar is not None
        if progress is not None:
            self.progress = min(progress, 1.0)
            self.progress_bar.n = self.progress
            updated = True
        else:
            self.progress_bar.n = 1
        if msg is not None and msg != self.message:
            self.message = msg
            self.progress_bar.set_description(self.message, refresh=False)
            updated = True
        if updated:
            self.progress_bar.refresh()

    def done(self) -> None:
        if self.progress_bar:
            self.progress_bar.close()
            self.progress_bar = None

    def series(
        self,
        items: Iterable[T],
        noun: str = "",
        labels: PySeq[str] | None = None,
        start: float | None = None,
        end: float = 1.0,
        count: int | None = None,
    ) -> Generator[T]:
        """Wrap a looped-over list with a progress bar"""
        # TODO optimise label creation
        if count is None:
            if not isinstance(items, Sized):
                items = list(items)
            count = len(items)
        if count == 0:
            # nothing to do
            return

        if start is None:
            start = 0.0
        step = (end - start) / count
        if labels:
            assert len(labels) == count
        elif count == 1:
            labels = [""]
        else:
            if noun:
                noun += " "
            template = f"{noun}%{len(str(count))}d/{count}"
            labels = [template % (i + 1) for i in range(count)]
        for i, item in enumerate(items):
            self.display(msg=labels[i], progress=start + step * i)
            yield item
        self.display(progress=end)

    def write(self, *args: Any, **kw: Any) -> None:
        if self.progress_bar_type and len(kw) < 3 and not in_jupyter():
            self.progress_bar_type.write(*args, **kw)
        else:
            pass

    def imap(
        self,
        f: Callable[[T], R],
        s: Collection[T],
        mininterval: float = 1.0,
        parallel: bool = False,
        par_kw: dict[str, Any] | None = None,
        **kw: Any,
    ) -> Generator[R]:
        self.mininterval = mininterval
        if parallel:
            # TODO document parallel.map arguments
            par_kw = par_kw or {}
            results: Iterable[R] = PAR.imap(f, s, **par_kw)
        else:
            results = map(f, s)
        yield from self.series(results, count=len(s), **kw)

    def map(self, f: Callable[[T], R], s: Collection[T], **kw: Any):
        return list(self.imap(f, s, **kw))


class NullContext(ProgressContext):
    """A UI context which discards all output.  Useful on secondary MPI cpus,
    and other situations where all output is suppressed"""

    def subcontext(self, *args: Any, **kw: Any) -> Self:
        return self

    def display(self, *args: Any, **kw: Any) -> None:
        pass

    def done(self) -> None:
        pass


NULL_CONTEXT = NullContext()
CURRENT = threading.local()
CURRENT.context = None


def display_wrap(slow_function: Callable[P, R]) -> Callable[P, R]:
    """Decorator which give the function its own UI context.
    The function will receive an extra argument, 'ui',
    which is used to report progress etc."""

    @functools.wraps(slow_function)
    def f(*args: P.args, **kw: P.kwargs) -> R:
        if getattr(CURRENT, "context", None) is None:
            klass: type[tqdm | notebook.tqdm | LogFileOutput] | None
            if sys.stdout.isatty():
                klass = tqdm
            elif in_jupyter():
                klass = notebook.tqdm
            elif isinstance(sys.stdout, io.FileIO):
                klass = LogFileOutput
            else:
                klass = None

            if klass is None:
                CURRENT.context = NULL_CONTEXT
            else:
                CURRENT.context = ProgressContext(klass)
        parent = CURRENT.context
        show_progress = kw.pop("show_progress", None)
        subcontext = NULL_CONTEXT if show_progress is False else parent.subcontext()
        kw["ui"] = CURRENT.context = subcontext
        try:
            result = slow_function(*args, **kw)
        finally:
            CURRENT.context = parent
            subcontext.done()
        return result

    return f


@display_wrap
def subdemo(ui: ProgressContext) -> None:
    for _j in ui.series(list(range(10))):
        time.sleep(0.1)


@display_wrap
def demo(ui: ProgressContext) -> None:
    ui.write("non-linebuffered output, tricky but look:")
    for i in ui.series(list(range(10))):
        time.sleep(0.6)
        for i in imap(fun, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], ui):
            ui.write(str(i))
        ui.write(str(i) + ".")
    ui.write("done")


@display_wrap
def imap(f: Callable[[T], R], s: Collection[T], ui: ProgressContext) -> Generator[R]:
    yield from ui.imap(f, s)


def fun(inp: T) -> T:
    time.sleep(0.1)
    return inp


if __name__ == "__main__":
    demo(ProgressContext(tqdm))
