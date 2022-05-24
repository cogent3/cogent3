import functools
import io
import sys
import threading
import time

from tqdm import notebook, tqdm

from cogent3.util import parallel as PAR


__author__ = "Sheng Han Moses Koh"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Sheng Han Moses Koh"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"

from cogent3.util.misc import in_jupyter


class LogFileOutput:
    """A fake progress bar for when progress bars are impossible"""

    def __init__(self, **kw):
        self.n = 0
        self.message = ""
        self.t0 = time.time()
        self.lpad = ""
        self.output = sys.stdout  # sys.stderr

    def set_description(self, desc="", refresh=False):
        self.message = desc

    def close(self):
        pass

    def refresh(self):
        if self.message:
            delta = f"+{int(time.time() - self.t0)}"
            progress = int(100 * self.n + 0.5)
            print(
                "%s %5s %3i%% %s" % (self.lpad, delta, progress, str(self.message)),
                file=self.output,
            )


class ProgressContext:
    def __init__(self, progress_bar_type=None, depth=-1, message=None, mininterval=1.0):
        self.progress_bar_type = progress_bar_type
        self.progress_bar = None
        self.progress = 0
        self.depth = depth
        self.message = message
        self.mininterval = mininterval

    def set_new_progress_bar(self):
        if self.progress_bar_type:
            self.progress_bar = self.progress_bar_type(
                total=1,
                position=self.depth,
                leave=True,
                bar_format="{desc} {percentage:3.0f}%|{bar}|{elapsed}<{remaining}",
                mininterval=self.mininterval,
                dynamic_ncols=True,
            )

    def subcontext(self, *args, **kw):
        return ProgressContext(
            progress_bar_type=self.progress_bar_type,
            depth=self.depth + 1,
            message=self.message,
            mininterval=self.mininterval,
        )

    def display(self, msg=None, progress=None):
        if not self.progress_bar:
            self.set_new_progress_bar()
        updated = False
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

    def done(self):
        if self.progress_bar:
            self.progress_bar.close()
            self.progress_bar = None

    def series(self, items, noun="", labels=None, start=None, end=1.0, count=None):
        """Wrap a looped-over list with a progress bar"""
        # todo optimise label creation
        if count is None:
            if not hasattr(items, "__len__"):
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
            labels = [template % (i + 1) for i in range(0, count)]
        for (i, item) in enumerate(items):
            self.display(msg=labels[i], progress=start + step * i)
            yield item
        self.display(progress=end)

    def write(self, *args, **kw):
        if self.progress_bar_type and len(kw) < 3 and not in_jupyter():
            self.progress_bar_type.write(*args, **kw)
        else:
            print(*args, **kw)

    def imap(self, f, s, mininterval=1.0, parallel=False, par_kw=None, **kw):
        self.mininterval = mininterval
        if parallel:
            # todo document parallel.map arguments
            par_kw = par_kw or {}
            results = PAR.imap(f, s, **par_kw)
        else:
            results = map(f, s)
        for result in self.series(results, count=len(s), **kw):
            yield result

    def map(self, f, s, **kw):
        return list(self.imap(f, s, **kw))


class NullContext(ProgressContext):
    """A UI context which discards all output.  Useful on secondary MPI cpus,
    and other situations where all output is suppressed"""

    def subcontext(self, *args, **kw):
        return self

    def display(self, *args, **kw):
        pass

    def done(self):
        pass


NULL_CONTEXT = NullContext()
CURRENT = threading.local()
CURRENT.context = None


def display_wrap(slow_function):
    """Decorator which give the function its own UI context.
    The function will receive an extra argument, 'ui',
    which is used to report progress etc."""

    @functools.wraps(slow_function)
    def f(*args, **kw):
        if getattr(CURRENT, "context", None) is None:
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
        if show_progress is False:
            subcontext = NULL_CONTEXT
        else:
            subcontext = parent.subcontext()
        kw["ui"] = CURRENT.context = subcontext
        try:
            result = slow_function(*args, **kw)
        finally:
            CURRENT.context = parent
            subcontext.done()
        return result

    return f


@display_wrap
def subdemo(ui):
    for j in ui.series(list(range(10))):
        time.sleep(0.1)
    return


@display_wrap
def demo(ui):
    ui.write("non-linebuffered output, tricky but look:")
    for i in ui.series(list(range(10))):
        time.sleep(0.6)
        for i in imap(fun, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]):
            ui.write(str(i))
        ui.write(str(i) + ".")
    ui.write("done")


@display_wrap
def imap(f, s, ui):
    for result in ui.imap(f, s):
        yield result


def fun(inp):
    time.sleep(0.1)
    return inp


if __name__ == "__main__":
    demo()
