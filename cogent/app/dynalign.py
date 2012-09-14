#!/usr/bin/env python

from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath, ApplicationError
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    MixedParameter,Parameters

from sys import platform
from os  import remove,system,mkdir,getcwd,close

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class Dynalign(CommandLineApplication):
    """Application controller for Dynalign

    Input:         Sequences for input in Dynalign are assumed to be in the
                   following format: ;(first line of file) Comments must start
                   with a ;semicolon ;There can be any number of comments A
                   single line title must immediately follow:
                   AAA GCGG UUTGTT UTCUTaaTCTXXXXUCAGG1  where the terminal 1
                   is required at the end.
    
    pknotsRG is a tool for thermodynamic folding of RNA secondary
    structures, including the class of canonical simple recursive
    pseudoknots.

    -alignment          The alignment file that will be created
    -M		        Maximum seperation parameter (larger then diff in length
                        of seq)
    -gap_cost	        0.4 good value
    -max_precent_diff	maximum % diff in free energy in suboptimal struct.
                        20 good starting point
    -bp_window	        specifies how different the suboptimal structures must
                        be from each other. 2 good starting point
    -align_window	specifies how different alignments must be, 1 good starting
                        point.
    -single_bp_inserts	specifies whether single base pair inserts are allowed;
                        0 = no; 1 = yes.
    -[savefile]	        optional

    """
    _parameters = {
        '-alignment':ValuedParameter(Prefix='',Name='',Value='align',Delimiter=' '),
        '-M':ValuedParameter(Prefix='',Name='',Value=10,Delimiter=' '),
        '-gap_cost':ValuedParameter(Prefix='',Name='',Value=0.4,Delimiter=' '),
        '-max_structures':ValuedParameter(Prefix='',Name='',Value=8,Delimiter=' '),
        '-max_percent_diff':ValuedParameter(Prefix='',Name='',Value=20,\
                                            Delimiter=' '),
        '-bp_window':ValuedParameter(Prefix='',Name='',Value=2,Delimiter=' '),
        '-align_window':ValuedParameter(Prefix='',Name='',Value=1,Delimiter=' '),
        '-single_bp_inserts':ValuedParameter(Prefix='',Name='',Value=1,\
                                             Delimiter=' '),}

    _command = 'dynalign'
    _input_handler = '_input_as_string'

    def _input_as_string(self,data):
        """Return data as a string

        Data = list with two file paths
        ex: data = ['path1,'path2'']"""

        inputFiles = ''
        self._input_filename = []
        for i in data:
            self._input_filename.append(i)
            inputFiles  = ' '.join([inputFiles,i])
        outputFile = ' '.join(['ct1', 'ct2'])
        inputFiles = ' '.join([inputFiles,outputFile])
        return inputFiles

    def _input_as_lines(self,data):
        """ Write a seq of lines to a temp file and return the filename string
        
        data: a sequence to be written to a file, each element of the 
        sequence will compose a line in the file. 
        ex. data = [file1,file2], file1 and 2 is lists of lines -> 
        [[List],[List]]
        
        Note: '\n' will be stripped off the end of each sequence element
        before writing to a file in order to avoid multiple new lines
        accidentally be written to a file
        """
        inputFiles = ''
        self._input_filename = []
        for el in data:
            filename = self.getTmpFilename(self.WorkingDir)
            self._input_filename.append(filename)
            data_file = open(filename,'w')
            data_to_file = '\n'.join([str(d).strip('\n') for d in el])
            data_file.write(data_to_file)
            data_file.close()
            inputFiles =  ' '.join([inputFiles,filename])
        outputFiles = ' '.join(['ct1', 'ct2'])
        inputFiles = ' '.join([inputFiles,outputFiles])
        return inputFiles
        

    
    def _get_result_paths(self,data):
        """
        data: the data the instance of the application is called on
        """
        result = {}
        result['seq_1_ct'] =\
            ResultPath(Path=(self.WorkingDir+'ct1'))
        result['seq_2_ct'] =\
            ResultPath(Path=(self.WorkingDir+'ct2'))
        result['alignment'] =\
            ResultPath(Path=(self.WorkingDir+'align'))

        return result

#Below from Cogent/app/util.py modified to accommodate dynalign
    
    def __call__(self,data=None):
        """Run the application with the specified kwargs on data

        Overides the __call__ function in util.py becasue of the special
        circumstance surrounding the command line input.
        
        data: anything that can be cast into a string or written out to
        a file. Usually either a list of things or a single string or 
        number. input_handler will be called on this data before it 
        is passed as part of the command-line argument, so by creating
        your own input handlers you can customize what kind of data
        you want you application to accept
        """
        input_handler = self.InputHandler
        suppress_stdout = self.SuppressStdout
        suppress_stderr = self.SuppressStderr
        if suppress_stdout:
            outfile = '/dev/null'
        else:
            outfile = self.getTmpFilename(self.WorkingDir)
        if suppress_stderr:
            errfile = '/dev/null'
        else:
            errfile = self.getTmpFilename(self.WorkingDir)
        if data is None:
            input_arg = ''
        else:
            input_arg = getattr(self,input_handler)(data)

        # Build up the command, consisting of a BaseCommand followed by
        # input and output (file) specifications
        first,second=self.BaseCommand
        command = self._command_delimiter.join(filter(None,\
            [first,input_arg,second,'>',outfile,'2>',errfile]))
        if self.HaltExec: 
            raise AssertionError, "Halted exec with command:\n" + command
        # The return value of system is a 16-bit number containing the signal 
        # number that killed the process, and then the exit status. 
        # We only want to keep the exit status so do a right bitwise shift to 
        # get rid of the signal number byte
        exit_status = system(command) >> 8
      
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
        if self._input_filename:
            for f in self._input_filename:
                remove(f)
            self._input_filename = None

        return result


    def _get_base_command(self):
        """ Returns the full command string

            Overides the __call__ function in util.py becasue of the special
            circumstance surrounding the command line input.

            input_arg: the argument to the command which represents the input 
                to the program, this will be a string, either 
                representing input or a filename to get input from
        """
        command_part1 = []
        command_part2 = []
        # Append a change directory to the beginning of the command to change 
        # to self.WorkingDir before running the command
        cd_command = ''.join(['cd ',self.WorkingDir,';'])
        if self._command is None:
            raise ApplicationError, '_command has not been set.'
        command = self._command
        
        command_part1.append(cd_command)
        command_part1.append(command)

        lista = [self.Parameters['-alignment'],\
                  self.Parameters['-M'],\
                  self.Parameters['-gap_cost'],\
                  self.Parameters['-max_structures'],\
                  self.Parameters['-max_percent_diff'],\
                  self.Parameters['-bp_window'],\
                  self.Parameters['-align_window'],\
                  self.Parameters['-single_bp_inserts']]
        
        command_part2.append(self._command_delimiter.join(filter(\
            None,(map(str,lista)))))
      
        return self._command_delimiter.join(command_part1).strip(),\
               self._command_delimiter.join(command_part2).strip()
    
    BaseCommand = property(_get_base_command)
