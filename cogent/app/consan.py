#!/usr/bin/env python

from os import remove, system, rmdir, mkdir
from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath, FilePath, ApplicationError
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    MixedParameter,Parameters, _find_synonym, is_not_None

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
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
        self._input_filename = data
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
        self._input_filename = []    
        for i in range(2):
            filename = self.getTmpFilename(self.WorkingDir)
            self._input_filename.append(filename)
            
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

    # need to override __call__ to remove files properly
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
        # The return value of system is a 16-bit number containing the signal 
        # number that killed the process, and then the exit status. 
        # We only want to keep the exit status so do a right bitwise shift to 
        # get rid of the signal number byte
        tmp_dir = ''.join([self.WorkingDir, 'tmp'])
        mkdir(tmp_dir)
        exit_status = system(command) >> 8
        rmdir(tmp_dir)

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
                for f in self._input_filename:
                    remove(f)
                self._input_filename = None

        return result

