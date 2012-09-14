#!/usr/bin/env python

"""Application controller for NUPACK v1.2 application

"""

import shutil
from os import remove, system, environ
from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath, FilePath, ApplicationError
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    MixedParameter,Parameters, _find_synonym, is_not_None

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"  

if 'NUPACK_DATA' not in environ:
    raise RuntimeError, \
        'NUPACK app controller requires the NUPACK_DATA environment variable'
nupack_data_dir = environ['NUPACK_DATA']
nupack_data_dna = 'dataS_G.dna'
nupack_data_rna = 'dataS_G.rna'

class Nupack(CommandLineApplication):
    """Application controller for Nupack_1.2 application 

    Predicts RNA secondary structure
    All pseudoknot-free secondary structures are allowed, as well as simple 
    pseudoknots.
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
    def __call__(self,data=None, remove_tmp=True):
        """Run the application with the specified kwargs on data
        
            data: anything that can be cast into a string or written out to
                a file. Usually either a list of things or a single string or 
                number. input_handler will be called on this data before it 
                is passed as part of the command-line argument, so by creating
                your own input handlers you can customize what kind of data
                you want your application to accept

            remove_tmp: if True, removes tmp files
        """
        input_handler = self.InputHandler
        suppress_stdout = self.SuppressStdout
        suppress_stderr = self.SuppressStderr
        if suppress_stdout:
            outfile = FilePath('/dev/null')
        else:
            outfile = self.getTmpFilename(self.TmpDir)
        if suppress_stderr:
            errfile = FilePath('/dev/null')
        else:
            errfile = FilePath(self.getTmpFilename(self.TmpDir))
        if data is None:
            input_arg = ''
        else:
            input_arg = getattr(self,input_handler)(data)

        # Build up the command, consisting of a BaseCommand followed by
        # input and output (file) specifications
        command = self._command_delimiter.join(filter(None,\
            [self.BaseCommand,str(input_arg),'>',str(outfile),'2>',\
                str(errfile)]))
        if self.HaltExec:
            raise AssertionError, "Halted exec with command:\n" + command

        # copy over data files
        nupack_data_dna_src = '/'.join([nupack_data_dir, nupack_data_dna])
        nupack_data_rna_src = '/'.join([nupack_data_dir, nupack_data_rna])
        shutil.copy(nupack_data_dna_src, self.WorkingDir)
        shutil.copy(nupack_data_rna_src, self.WorkingDir)

        # The return value of system is a 16-bit number containing the signal 
        # number that killed the process, and then the exit status. 
        # We only want to keep the exit status so do a right bitwise shift to 
        # get rid of the signal number byte
        # NOTE: we copy the data files to the working directory first
        exit_status = system(command) >> 8

        # remove data files
        nupack_data_dna_dst = ''.join([self.WorkingDir, nupack_data_dna])
        nupack_data_rna_dst = ''.join([self.WorkingDir, nupack_data_rna])
        remove(nupack_data_dna_dst)
        remove(nupack_data_rna_dst)

        # Determine if error should be raised due to exit status of 
        # appliciation
        if not self._accept_exit_status(exit_status):
            raise ApplicationError, \
             'Unacceptable application exit status: %s, command: %s'\
                % (str(exit_status),command)

        # open the stdout and stderr if not being suppressed
        out = None
        if not suppress_stdout:
            out = open(outfile,"r")
        err = None
        if not suppress_stderr:
            err = open(errfile,"r")

        result =  CommandLineAppResult(out,err,exit_status,\
            result_paths=self._get_result_paths(data))

        # Clean up the input file if one was created
        if remove_tmp:
            if self._input_filename:
                remove(self._input_filename)
                self._input_filename = None

        return result
    
