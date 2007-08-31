#!/usr/bin/env python

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

class Consan(CommandLineApplication):
    """Application controller for CONSAN v1.1"""

    _parameters = {
        '-m':ValuedParameter(Prefix='-',Name='m'),
        '-M':ValuedParameter(Prefix='-',Name='M'), 
        '-C':ValuedParameter(Prefix='-',Name='C'), 
        '-P':ValuedParameter(Prefix='-',Name='P'), 
        '-V':FlagParameter(Prefix='-',Name='V'), 
        '-f':FlagParameter(Prefix='-',Name='f'),
        '-x':FlagParameter(Prefix='-',Name='x'),
        '-t':FlagParameter(Prefix='-',Name='t')}
    _command = 'sfold'
    _input_handler='_input_as_string'


    def _input_as_string(self,data):
        """
        Takes two files in a list as input
        eg. data = [path1,path2]
        """
        inputFiles = ' '.join(data)
        return inputFiles

    def _input_as_lines(self,data):
        """ 
        Writes to first sequences(fasta) in a list to two temp files
       
        data: a sequence to be written to a file, each element of the 
        sequence will compose a line in the file

        Data should be in the following format:
        data = ['>tag1','sequence1','>tag2','sequence2']

        Note: '\n' will be stripped off the end of each sequence element
        before writing to a file in order to avoid multiple new lines
        accidentally be written to a file
        """
        inputFiles = ''
        for i in range(2):
            filename = self._input_filename = self.getTmpFilename(self.WorkingDir)
            
            data_file = open(filename,'w')
            if i == 0:
                data_to_file = '\n'.join(data[:2])
                tmp1 = filename
            else:
                data_to_file = '\n'.join(data[2:])
                tmp2 = filename
            data_file.write(data_to_file)
            data_file.close()
        inputFiles = ' '.join([tmp1,tmp2])
        return inputFiles
