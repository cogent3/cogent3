#!/usr/bin/env python

"""Application controller for NUPACK v1.2 application

IMPORTANT!!!!
The file; dataS_G.rna, must be present in exec directory!!!
"""

from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    MixedParameter,Parameters, _find_synonym, is_not_None

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"  

class Nupack(CommandLineApplication):
    """Application controller for Nupack_1.2 application 

    Predicts RNA secondary structure
    All pseudoknot-free secondary structures are allowed, as well as simple 
    pseudoknots.

    IMPORTANT!!!!
    The file; dataS_G.rna, must be present in exec directory!!!
    """

    _command = 'Fold.out'
    _input_handler = '_input_as_string'

    def _input_as_lines(self,data):
        """ Write a seq of lines to a temp file and return the filename string
        
            data: a sequence to be written to a file, each element of the 
                sequence will compose a line in the file

            Note: '\n' will be stripped off the end of each sequence element
                before writing to a file in order to avoid multiple new lines
                accidentally be written to a file
        """
        filename = self._input_filename = self.getTmpFilename(self.WorkingDir)
        data_file = open(filename,'w')
        data_to_file = '\n'.join([str(d).strip('\n') for d in data])
        data_file.write(data_to_file)
        data_file.write('\n') #needs a new line att the end of input
        data_file.close()
        return filename
    
    def _get_result_paths(self,data):
        """Return a dict of ResultPath objects representing all possible output

        This dictionary will have keys based
        on the name that you'd like to access the file by in the 
        CommandLineAppResult object that will be created, and the values
        which are ResultPath objects.
        """
        result = {}
        try:
            f = open((self.WorkingDir+'out.pair'))
            f.close()
            result['pair'] =\
                ResultPath(Path=(self.WorkingDir+'out.pair'))
        except IOError:
            pass
        
        try:
            f = open((self.WorkingDir+'out.ene'))
            f.close()
            result['ene'] =\
                ResultPath(Path=(self.WorkingDir+'out.ene'))
        except IOError:
            pass
        
        return result
    
