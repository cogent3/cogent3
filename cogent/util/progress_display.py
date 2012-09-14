# -*- coding: utf-8 -*-

"""The hook for new user interfaces to take control of progress bar and status
display is pass setupRootUiContext an instance of their UI context class, with
the same methods as the _Context class defined here.

Long-running functions can be decorated with @display_wrap, and will then be
given the extra argument 'ui'.  'ui' is a ProgressContext instance with methods
.series(), .imap(), .map() and .display(), any one of which will cause a 
progress-bar to be displayed.

@display_wrap
def long_running_function(..., ui)
    ui.display(msg, progress)  # progress is between 0.0 and 1.0
  or
    for item in ui.map(items, function)
  or
    for item in ui.imap(items, function)
  or
    for item in ui.series(items)
"""

from __future__ import with_statement, division
import sys, time, contextlib, functools, warnings
import os, atexit
import threading
import itertools
from cogent.util import parallel, terminal

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.5.3"

try:
    curses_terminal = terminal.CursesOutput()
except terminal.TerminalUnavailableError:
    curses_terminal = None
else:
    CODES = curses_terminal.getCodes()
    if not (CODES['BOL'] and CODES['UP'] and CODES['CLEAR_EOL']):
        curses_terminal = None  # Terminal too primitive
    else:
        BOL = CODES['BOL']
        CLEAR = CODES['UP'] + BOL + CODES['CLEAR_EOL']
        if CODES['GREEN']:
            bar_template = CODES['GREEN'] + '%s' + CODES['NORMAL'] + '%s'
            def terminal_progress_bar(dots, width):
                return bar_template % ('█' * dots, '█' * (width-dots))
        else:
            def terminal_progress_bar(dots, width):
                return '.' * dots

class TextBuffer(object):
    """A file-like object which accumulates written text.  Specialised for 
    output to a curses terminal in that it uses CLEAR and re-writing to extend
    incomplete lines instead of just outputting or buffering them.  That
    allows the output to always end at a newline, ready for a progress bar
    to be shown, without postponing output of any incomplete last line."""
    
    def __init__(self):
        self.chunks = []
        self.pending_eol = False
        
    def write(self, text):
        self.chunks.append(text)
    
    # multiprocessing calls these
    def flush(self):
        pass
    def isatty(self):
        return False
        
    def regurgitate(self, out):
        if self.chunks:
            text = ''.join(self.chunks)
            if self.pending_eol:
                out.write(CLEAR)
            #out.write(CODES['YELLOW'])
            out.write(text)
            if text.endswith('\n'):
                self.pending_eol = False
                self.chunks = []
            else:
                self.pending_eol = True
                self.chunks = [text.split('\n')[-1]]
                out.write('\n')
            #out.write(CODES['NORMAL'])


class ProgressContext(object):
    """The interface by which cogent algorithms can report progress to the
    user interface.  Calls self.progress_bar.set(progress, message)"""
    
    def __init__(self, progress_bar=None, prefix=None, base=0.0, segment=1.0, 
            parent=None, rate=1.0):
        self.progress_bar = progress_bar
        self.desc = ''
        self.base = base
        self.segment = segment
        self.progress = 0
        self.current = 1
        if parent is None:
            self.depth = 0
            self.parent = None
            self.t_last = 0
        else:
            assert progress_bar is parent.progress_bar
            self.depth = parent.depth + 1
            self.parent = parent
            self.t_last = parent.t_last
        self.msg = ''
        self.prefix = prefix or []
        self.message = self.prefix + [self.msg]
        self._max_text_len = 0
        self.max_depth = 2
        self.rate = rate
        
    def subcontext(self):
        """For any sub-task which may want to report its own progress, but should
        not get its own progress bar."""
        if self.depth == self.max_depth:
            return NullContext()
        return ProgressContext(
            progress_bar = self.progress_bar, 
            prefix = self.message,
            base = self.base+self.progress*self.segment,
            segment = self.current*self.segment,
            parent = self,
            rate = self.rate) 
    
    def display(self, msg=None, progress=None, current=0.0):
        """Inform the UI that we are are at 'progress' of the way through and 
        will be doing 'msg' until we reach and report at progress+current.
        """
        if self.depth > 0:
            msg = None
            
        updated = False
        if progress is not None:
            self.progress = min(progress, 1.0)
            updated = True

        if current is not None:
            self.current = current
            updated = True

        if msg is not None and msg != self.msg:
            self.msg = self.message[-1] = msg
            updated = True

        if updated and (
                (self.depth==0 and self.progress in [0.0, 1.0]) or 
                time.time() > self.t_last + self.rate):
            self.render()

    def render(self):
        self.progress_bar.set(self.base+self.progress*self.segment, self.message[0])
        self.t_last = time.time()

    def done(self):
        if self.depth == 0:
            self.progress_bar.done()
    
    # Not much point while cogent is still full of print statements, but
    # .info() (and maybe other logging analogues such as .warning()) would
    # avoid the need to capture stdout:
    
    #def info(self, text):
    #    """Display some information which may be more than fleetingly useful, 
    #    such as a summary of intermediate statistics or a very mild warning.  
    #    A GUI should make this information retrievable but not intrusive.
    #    For terminal UIs this is equivalent to printing"""
    #    raise NotImplementedError
        
    def series(self, items, noun='', labels=None, start=None, end=1.0, count=None):
        """Wrap a looped-over list with a progress bar"""
        if count is None:
            if not hasattr(items, '__len__'):
                items = list(items)
            count = len(items)
        if start is None:
            start = 0.0
        step = (end-start) / count
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
            self.display(msg=labels[i], progress=start+step*i, current=step)
            yield item
        self.display(progress=end, current=0)
        
    def imap(self, f, s, pure=True, **kw):
        """Like itertools.imap() but with a progress bar"""
        results = (parallel if pure else itertools).imap(f, s)
        for result in self.series(results, count=len(s), **kw):
            yield result
    
    def eager_map(self, f, s, **kw):
        """Like regular Python2 map() but with a progress bar"""
        return list(self.imap(f,s, **kw))
    
    def map(self, f, s, **kw):
        """Synonym for eager_map, unlike in Python3"""
        return self.eager_map(f, s, **kw)


class NullContext(ProgressContext):
    """A UI context which discards all output.  Useful on secondary MPI cpus, 
    and other situations where all output is suppressed"""
    def subcontext(self, *args, **kw):
        return self
        
    def display(self, *args, **kw):
        pass

    def done(self):
        pass


class LogFileOutput(object):
    """A fake progress bar for when progress bars are impossible"""
    def __init__(self):
        self.t0 = time.time()
        self.lpad = ''
        self.output = sys.stdout # sys.stderr
    
    def done(self):
        pass
    
    def set(self, progress, message):        
        if message:
            delta = '+%s' % int(time.time() - self.t0)
            progress = int(100*progress+0.5)
            print >>self.output, "%s %5s %3i%% %s" % (
                    self.lpad, delta, progress,
                    str(message.encode('utf8')))
            

class CursesTerminalProgressBar(object):
    """Wraps stdout and stderr, displaying a progress bar via simple 
    ascii/curses art and scheduling other output around its redraws."""
    def __init__(self):
        global curses_terminal
        assert curses_terminal is not None
        self.curses_terminal = curses_terminal
        self.stdout = sys.stdout
        self.stderr = sys.stderr
        self.stdout_log = TextBuffer()
        self.stderr_log = TextBuffer()
        self.lines = []
        self.chunks = []
        self.pending_eol = False
        self.line_count = 0
        (sys.stdout, sys.stderr, self._stdout, self._stderr) = (
                self.stdout_log, self.stderr_log, sys.stdout, sys.stderr)
    
    def done(self):
        self.set(None, None)
        (sys.stdout, sys.stderr) = (self._stdout, self._stderr)
    
    def set(self, progress, message):
        """Clear the existing progress bar, write out any accumulated
        stdout and stderr, then draw the updated progress bar."""
        cols = self.curses_terminal.getColumns()
        width = cols - 1
        if progress is not None:
            assert 0.0 <= progress <= 1.0, progress
            dots = int(progress * width)
            bar = terminal_progress_bar(dots, width)

        if self.line_count:
            self.stderr.write(CLEAR * (self.line_count))
        else:
            self.stderr.write(BOL)
        self.stdout_log.regurgitate(self.stdout)
        self.stderr_log.regurgitate(self.stderr)
        
        if progress is not None:
            self.stderr.writelines([bar, '\n'])
        
        if message is not None:
            self.stderr.writelines([str(message[:width].encode('utf8')), u'\n'])
        self.line_count = (progress is not None) + (message is not None)


NULL_CONTEXT = NullContext()
CURRENT = threading.local()
CURRENT.context = None

class RootProgressContext(object):
    """The context between long running jobs, when there is no progress bar"""
    
    def __init__(self, pbar_constructor, rate):
        self.pbar_constructor = pbar_constructor
        self.rate = rate
        
    def subcontext(self):
        pbar = self.pbar_constructor()
        return ProgressContext(pbar, rate=self.rate)


def setupRootUiContext(progressBarConstructor=None, rate=None):
    """Select a UI Context type depending on system environment"""
    if parallel.getCommunicator().Get_rank() != 0:
        klass = None
    elif progressBarConstructor is not None:
        klass = progressBarConstructor
    elif curses_terminal and sys.stdout.isatty():
        klass = CursesTerminalProgressBar
    elif isinstance(sys.stdout, file):
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
        CURRENT.context = RootProgressContext(klass, rate)


def display_wrap(slow_function):
    """Decorator which give the function its own UI context.
    The function will receive an extra argument, 'ui', 
    which is used to report progress etc."""
    @functools.wraps(slow_function)
    def f(*args, **kw):
        if getattr(CURRENT, 'context', None) is None:
            setupRootUiContext()
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
    for j in ui.series(range(10)):
        time.sleep(0.1)
    return
    
@display_wrap
def demo(ui):
    print "non-linebuffered output, tricky but look:",
    for i in ui.series(range(10)):
        time.sleep(.6)
        if i == 5:
            print '\nhalfway through, a new line: ',
        if i % 2:
            subdemo()
        print i, ".", 
    print "done"

if __name__ == '__main__':
    #setupRootUiContext(rate=0.2)
    demo()
    
# This messes up interactive shells a bit:
#CURRENT.start()
#atexit.register(CURRENT.done)



