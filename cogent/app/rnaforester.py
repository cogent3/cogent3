#!/usr/bin/env python

from os import mkdir,getcwd,system,remove,close
from random import choice
from cogent.util.misc import app_path
from cogent.app.util import CommandLineApplication, CommandLineAppResult,\
    ResultPath, ApplicationNotFoundError, ApplicationError
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    MixedParameter,Parameters, _find_synonym, is_not_None

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman","Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class RNAforester(CommandLineApplication):
    """Application controller for RNAforester application"""
     #Not all parameters added!!
    _parameters = {
        '-d':FlagParameter(Prefix='-',Name='d',Value=False),
        '-r':FlagParameter(Prefix='-',Name='r',Value=False),
        '-m':FlagParameter(Prefix='-',Name='m',Value=True),
        '-p':FlagParameter(Prefix='-',Name='p',Value=False)}

    _command1 = 'RNAshapes -C -c 20 -f'
    _command2 = 'RNAforester'
    _input_handler = '_input_as_string'

    def _input_as_lines(self,data):
        """
        """
        data = ' '.join([super(RNAforester,self)._input_as_lines(data),'-o f','|'])
        return data

    def _input_as_string(self,data):
        """Makes data the value of a specific parameter
        
        This method returns the empty string. The parameter will be printed
        automatically once set.
        """
        data = ' '.join([data,'-o f','|'])
        return data

    def _error_on_missing_application(self,params):
        """ Raise an ApplicationNotFoundError if the app is not accessible
        """
        if not app_path('RNAforester'):
            raise ApplicationNotFoundError,\
             "Cannot find RNAforester. Is it installed? Is it in your path?"
        if not app_path('RNAshapes'):
            raise ApplicationNotFoundError,\
             "Cannot find RNAshapes. Is it installed? Is it in your path?"

#Override these functions to biuld up the command  
    def __call__(self,data=None, remove_tmp=True):
        """Run the application with the specified args on data
        
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

        first,second = self.BaseCommand
        command = self._command_delimiter.join(filter(None,\
            [first,input_arg,second,'>',outfile,'2>',errfile]))

        #print 'COMMAND',command

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
        if remove_tmp:
            if self._input_filename:
                remove(self._input_filename)
                self._input_filename = None

        remove(''.join([self.WorkingDir, 'cluster.dot']))
        remove(''.join([self.WorkingDir, 'test.out']))
        remove(''.join([self.WorkingDir, 'ShapesStderr']))

        return result

    def _get_base_command(self):
        """ Returns the full command string 

            input_arg: the argument to the command which represents the input 
                to the program, this will be a string, either 
                representing input or a filename to get input from
        """
        command_part1 = []
        command_part2 = []
        # Append a change directory to the beginning of the command to change 
        # to self.WorkingDir before running the command
        cd_command = ''.join(['cd ',self.WorkingDir,';'])
        if self._command1 is None:
            raise ApplicationError, '_command has not been set.'

        parameters = self.Parameters
        command1 = self._command1
        command2 = self._command2

        command_part1.append(cd_command)
        command_part1.append(command1)
        command_part1.append(''.join(['2> ', self.WorkingDir, 'ShapesStderr']))

        command_part2.append(command2)

        command_part2.append(self._command_delimiter.join(filter(\
            None,(map(str,parameters.values())))))
      
        return self._command_delimiter.join(command_part1).strip(),\
               self._command_delimiter.join(command_part2).strip()
    
    BaseCommand = property(_get_base_command)
    
