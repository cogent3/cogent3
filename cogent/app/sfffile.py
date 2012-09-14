#!/usr/bin/env python

from cogent.app.parameters import FlagParameter, ValuedParameter
from cogent.app.util import CommandLineApplication, ResultPath


"""Application controller for sfffile"""


__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"


class Sfffile(CommandLineApplication):
    """Simple sfffile application controller.
    """
    _options = {
        # output filepath
        '-o': ValuedParameter('-', 'o', Delimiter=' '),
        # file of accession numbers to be included
        '-i': ValuedParameter('-', 'i', Delimiter=' '),
        # file of accession numbers to be excluded
        '-e': ValuedParameter('-', 'e', Delimiter=' '),
        # file of custom trim points
        '-t': ValuedParameter('-', 't', Delimiter=' '),
        # number of cycles in output sff
        '-c': ValuedParameter('-', 'c', Delimiter=' '),
        # shortcut for -c 42
        '-gs20': FlagParameter('-', 'gs20'),
        # shortcut for -c 100
        '-gsflx': FlagParameter('-', 'gsflx'),
        # split multiplexed reads
        '-s': ValuedParameter('-', 's', Delimiter=' '),
        # custom MID configuration file
        '-mcf': ValuedParameter('-', 'mcf', Delimiter=' '),        
        # prevent propagation of sff index
        '-nmft': FlagParameter('-', 'nmft'),
        }
    _parameters = {}
    _parameters.update(_options)
    _input_handler = '_input_as_path'
    _command = 'sfffile'

    def _get_result_paths(self, data):
        """Collect the resultant SFF file in the results.

        Because cogent.app.util.CommandLineAppResult opens output
        files in text mode, this method may not be portable for
        Windows users.  A more portable solution would be to not use
        the app controller results, but instead specify the output SFF
        filepath manually via the '-o' parameter.
        """
        if self.Parameters['-o'].isOn():
            sff_path = self.Parameters['-o'].Value
        else:
            sff_path = '454Reads.sff'
        return {'sff': ResultPath(sff_path)}

    def _accept_exit_status(self, exit_status):
        """Accept an exit status of 0 for the sfffile program.
        """
        return exit_status == 0
