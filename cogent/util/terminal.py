# Copyright: 2008 Nadia Alramli
# License: BSD

import sys
try:
    import curses
except ImportError:
    curses = None

COLORS = "BLUE GREEN CYAN RED MAGENTA YELLOW WHITE BLACK".split()

# List of terminal controls, you can add more to the list.
_CONTROLS = {
    'BOL':'cr', 'UP':'cuu1', 'DOWN':'cud1', 'LEFT':'cub1', 'RIGHT':'cuf1',
    'CLEAR_SCREEN':'clear', 'CLEAR_EOL':'el', 'CLEAR_BOL':'el1',
    'CLEAR_EOS':'ed', 'BOLD':'bold', 'BLINK':'blink', 'DIM':'dim',
    'REVERSE':'rev', 'UNDERLINE':'smul', 'NORMAL':'sgr0',
    'HIDE_CURSOR':'cinvis', 'SHOW_CURSOR':'cnorm'
        }

class TerminalUnavailableError(RuntimeError):
    pass
    
class CursesOutput(object):
    def __init__(self):
        if curses is None or not hasattr(sys.stdout, 'fileno'):
            raise TerminalUnavailableError("No curses modules")
        try:
            curses.setupterm()
        except curses.error, detail:
            raise TerminalUnavailableError(detail)
    
    def getColumns(self):
        return curses.tigetnum('cols')
    
    def getLines(self):
        return curses.tigetnum('lines')
    
    def getCodes(self):
        # Get the color escape sequence template or '' if not supported
        # setab and setaf are for ANSI escape sequences
        bgColorSeq = curses.tigetstr('setab') or curses.tigetstr('setb') or ''
        fgColorSeq = curses.tigetstr('setaf') or curses.tigetstr('setf') or ''
        codes = {}

        for color in COLORS:
            # Get the color index from curses
            colorIndex = getattr(curses, 'COLOR_%s' % color)
            # Set the color escape sequence after filling the template with index
            codes[color] = curses.tparm(fgColorSeq, colorIndex)
            # Set background escape sequence
            codes['BG_%s' % color] = curses.tparm(bgColorSeq, colorIndex)
        for control in _CONTROLS:
            # Set the control escape sequence
            codes[control] = curses.tigetstr(_CONTROLS[control]) or ''
        
        return codes
