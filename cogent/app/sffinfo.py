#!/usr/bin/env python

from cogent.app.parameters import FlagParameter, ValuedParameter
from cogent.app.util import CommandLineApplication, ResultPath, FilePath,\
        ApplicationError


"""Application controller for sffinfo"""


__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"


class ManyValuedParameter(ValuedParameter):
    """Like a ValuedParameter, but stores a list of values."""

    def __init__(self, Prefix, Name, Values=None, Delimiter=None, Quote=None,
                 IsPath=False, ValueDelimiter=None):
        """Initialize a ManyValuedParameter object.

        Supports all keyword arguments of ValuedParameter, plus an
        additional keyword, ValueDelimiter.  This specifies the
        delimiter between consecutive parameter values.  The Delimiter
        keyword, as before, specifies the delimiter between the flag
        and the values.
        """
        self.ValueDelimiter = ValueDelimiter
        super(ManyValuedParameter, self).__init__(
            Prefix, Name, Value=Values, Delimiter=Delimiter, Quote=Quote,
            IsPath=IsPath)

    def on(self, val):
        """Alias for append().
        """
        return self.append(val)

    def append(self, val):
        """Appends val to the current list of values for this parameter.
        """
        if self.Value is None:
            self.Value = []
        if self.IsPath:
            val = FilePath(val)
        self.Value.append(val)

    def __str__(self):
        """Return the parameter as a string
        
        When turned on: string representation of the parameter
        When turned off: empty string
        """
        if self.isOff():
            return ''

        quote = self.Quote or ''
        delim = self.ValueDelimiter or ''
        value_str = delim.join(
            ['%s%s%s' % (quote, v, quote) for v in self.Value])

        parts = [self.Prefix, self.Name, self.Delimiter, value_str]
        return ''.join(map(str, filter(None, parts)))


class Sffinfo(CommandLineApplication):
    """Simple sffinfo application controller.
    """
    _options = {
        '-a': FlagParameter('-', Name='a'),
        '-s': FlagParameter('-', Name='s'),
        '-q': FlagParameter('-', Name='q'),
        '-f': FlagParameter('-', Name='f'),
        '-t': FlagParameter('-', Name='t'),
        '-n': FlagParameter('-', Name='n'),
        '-m': FlagParameter('-', Name='m'),
        'accno': ManyValuedParameter(None, None, ValueDelimiter=' ')
        }
    _parameters = {}
    _parameters.update(_options)
    _trailing_parameters = set(['accno'])
    _input_handler = '_input_as_path'
    _command = 'sffinfo'

    ## We have hacked both _get_base_command and InputHandler to
    ## support a positional command line parameter that must be
    ## specified after the input filename.
    ##
    ## This approach was inspired by the Gtcmpca app controller; see
    ## gctmpca.py for the official solution to this problem.

    def _get_base_command(self):
        """Returns the full command string 

        input_arg: the argument to the command which represents the
            input to the program, this will be a string, either
            representing input or a filename to get input from
        """
        #
        # Overridden from parent class, changing only 2 lines in
        # definition of parameters variable.
        #

        command_parts = []
        # Append a change directory to the beginning of the command to change 
        # to self.WorkingDir before running the command
        # WorkingDir should be in quotes -- filenames might contain spaces
        cd_command = ''.join(['cd ',str(self.WorkingDir),';'])
        if self._command is None:
            raise ApplicationError, '_command has not been set.'
        command = self._command
        parameters = dict([
            (name, param) for (name, param) in self.Parameters.items()
            if name not in self._trailing_parameters])
        
        command_parts.append(cd_command)
        command_parts.append(command)
        command_parts.append(self._command_delimiter.join(filter(\
            None,(map(str,parameters.values())))))
      
        return self._command_delimiter.join(command_parts).strip()

    BaseCommand = property(_get_base_command)

    def _set_input_handler(self, method_name):
        """Stores the selected input handler in a private attribute.
        """
        self.__InputHandler = method_name

    def _get_input_handler(self):
        return '_input_handler_decorator'

    InputHandler = property(_get_input_handler, _set_input_handler)

    def _input_handler_decorator(self, data):
        """Appends trailing parameters to selected input_handler's results.
        """
        # begin with output from real input handler...
        f = getattr(self, self.__InputHandler)
        input_parts = [str(f(data))]
        
        # ...and append trailing parameters
        for name in self._trailing_parameters:
            param_str = str(self.Parameters[name])
            if param_str:
                input_parts.append(param_str)

        return self._command_delimiter.join(input_parts).strip()

    def _accept_exit_status(self, exit_status):
        """Accept an exit status of 0 for the sfffile program.
        """
        return exit_status == 0


def sffinfo_from_file(sff_filepath):
    """Run sffinfo on sff_filepath and return output as file object.
    """
    a = Sffinfo()
    results = a(sff_filepath)
    return results['StdOut']
