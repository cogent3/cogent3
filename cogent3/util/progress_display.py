import time
import functools
import threading
import sys
import io
import os
from tqdm import tqdm, tqdm_notebook
from cogent3.util import parallel

__author__ = "Sheng Han Moses Koh"
__copyright__ = ""
__credits__ = ["Peter Maxwell", "Sheng Han Moses Koh"]
__license__ = "GPL"
__version__ = ""


class LogFileOutput:
    """A fake progress bar for when progress bars are impossible"""

    def __init__(self, total=1, depth=0, leave=False, bar_format=None):
        self.n = 0
        self.message = ''
        self.t0 = time.time()
        self.lpad = ''
        self.output = sys.stdout  # sys.stderr

    def set_description(self, desc='', refresh=False):
        self.message = desc

    def close(self):
        pass

    def refresh(self):
        if self.message:
            delta = '+%s' % int(time.time() - self.t0)
            progress = int(100 * self.n + 0.5)
            print("%s %5s %3i%% %s" % (
                self.lpad, delta, progress,
                str(self.message)), file=self.output)


class ProgressContext:
    def __init__(self, progress_bar_type=None, depth=-1,
                                    message=None, rate=1.0):
        self.progress_bar_type = progress_bar_type
        self.progress_bar = None
        self.progress = 0
        self.depth = depth
        self.message = message
        self.rate = rate

    def set_new_progress_bar(self):
        if self.progress_bar_type:
            self.progress_bar = self.progress_bar_type(total=1,
                            position=self.depth,
                            leave=True,
                            bar_format='{desc} {percentage:3.0f}%|{bar}| ')

    def subcontext(self, *args, **kw):
        return ProgressContext(
            progress_bar_type=self.progress_bar_type,
            depth=self.depth+1,
            message=self.message,
            rate=self.rate)

    def display(self, msg=None, progress=None, current=0.0):
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
            self.progress_bar.set_description(self.message,
                                                refresh=False)
            updated = True
        if updated:
            self.progress_bar.refresh()

    def done(self):
        if self.progress_bar:
            self.progress_bar.close()
            self.progress_bar = None

    def series(self, items, noun='', labels=None, start=None, end=1.0,
               count=None):
        """Wrap a looped-over list with a progress bar"""
        if count is None:
            if not hasattr(items, '__len__'):
                items = list(items)
            count = len(items)
        if start is None:
            start = 0.0
        step = (end - start) / count
        if labels:
            assert len(labels) == count
        elif count == 1:
            labels = ['']
        else:
            if noun:
                noun += ' '
            template = '%s%%%sd/%s' % (noun, len(str(count)), count)
            labels = [template % i for i in range(0, count)]
        for (i, item) in enumerate(items):
            self.display(msg=labels[i], progress=start +
                                                 step * i)
            yield item
        self.display(progress=end)

    def write(self, *args, **kw):
        if self.progress_bar_type and len(kw) < 3 and not using_notebook():
            self.progress_bar_type.write(*args, **kw)
        else:
            print(*args, **kw)

    def imap(self, f, s, workers=None, pure=True, **kw):
        if pure:
            results = parallel.map(f, s, workers)
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


def using_notebook():
    try:
        get_ipython()
        return True
    except NameError:
        return False


def display_wrap(slow_function):
    """Decorator which give the function its own UI context.
    The function will receive an extra argument, 'ui',
    which is used to report progress etc."""
    depth = 0

    @functools.wraps(slow_function)
    def f(*args, **kw):
        nonlocal depth
        if getattr(CURRENT, 'context', None) is None:
            rate = None
            if sys.stdout.isatty():
                klass = tqdm
            elif using_notebook():
                klass = tqdm_notebook
            elif isinstance(sys.stdout, io.FileIO):
                klass = LogFileOutput
                if rate is None:
                    rate = 5.0
            else:
                klass = None

            if klass is None:
                CURRENT.context = NULL_CONTEXT
            else:
                if rate is None:
                    rate = 0.1
                CURRENT.context = ProgressContext(klass, rate=rate)
        parent = CURRENT.context
        show_progress = kw.pop('show_progress', None)
        if show_progress is False:
            # PendingDeprecationWarning?
            subcontext = NULL_CONTEXT
        else:
            subcontext = parent.subcontext()
        kw['ui'] = CURRENT.context = subcontext
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
        time.sleep(.6)
        if i == 5:
            ui.write('halfway through, a new line: ')
        if i % 2:
            subdemo()
        ui.write(str(i)+".")
    ui.write("done")

if __name__ == '__main__':
    demo()
