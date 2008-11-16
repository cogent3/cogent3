#!/usr/bin/env python
import commands
from sys import platform
from os import remove,system,mkdir,getcwd,close,sep
from os.path import isabs
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    MixedParameter,Parameters, _find_synonym, is_not_None, FilePath
from cogent.util.misc import if_
from random import choice
from numpy import zeros, array, nonzero, max

__author__ = "Sandra Smit and Greg Caporaso"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Greg Caporaso", "Sandra Smit", "Micah Hamady",
                    "Jeremy Widmann", "Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"
 
#the following are used to create temp file names       
_chars = "abcdefghigklmnopqrstuvwxyz"
_all_chars = _chars + _chars.upper() + "0123456790"


class ApplicationError(OSError):
    pass
   
class ResultPath(object):
    """ Hold a file path a boolean value specifying whether file was written
    """
    def __init__(self,Path,IsWritten=True):
        """ Initialize the ResultPath object

            Path: a string representing the absolute or relative path where
                the file can be found
            IsWritten: a boolean specifying whether the file has been written,
                default = True
        """
        self.Path = FilePath(Path)
        self.IsWritten = IsWritten

class CommandLineAppResult(dict):
    """ Class for holding the result of a CommandLineApplication run """

    def __init__(self,out,err,exit_status,result_paths):
        """Initialization of CommandLineAppResult

        out: a file handler to the file containing the stdout
        err: a file handler to the file containing the stderr
        exit_status: the exit status of the program, 0 if run ok, 1 else.
        result_paths: dictionary containing ResultPath objects for each 
            output file that could be written
        """
        
        self['StdOut'] = out
        self['StdErr'] = err
        self['ExitStatus'] = exit_status
       
        self.file_keys = result_paths.keys()
        for key,value in result_paths.items():
            if value.IsWritten:
                try:
                    self[key] = open(value.Path)
                except IOError:
                    raise ApplicationError, 'Could not open %s' %value.Path
            else:
                self[key] = None

    def cleanUp(self):
        """ Delete files that are written by CommandLineApplication from disk
            
            WARNING: after cleanUp() you may still have access to part of 
                your result data, but you should be aware that if the file
                size exceeds the size of the buffer you will only have part 
                of the file. To be safe, you should not use cleanUp() until 
                you are done with the file or have copied it to a different 
                location.
        """
        file_keys = self.file_keys
        for item in file_keys:
            if self[item] is not None:
                self[item].close()
                remove(self[item].name)

        # remove input handler temp files
        if hasattr(self, "_input_filename"): 
            remove(self._input_filename)

    def __del__(self):
        """ Delete temporary files created by the CommandLineApplication 
        """
        if self['StdOut'] is not None:
            remove(self['StdOut'].name)
        if self['StdErr'] is not None:
            remove(self['StdErr'].name)

class Application(object):
    """ Generic Class for controlling an application """

    _command = None
    _command_delimiter = ' '
    _parameters = {}
    _synonyms = {}
    
    def __init__(self,params=None):
        """
            params: a dict of parameters which should be turned on where the 
                key is either the parameter id or a synonym for the parameter
                and the value is either the value for the parameter or None
        """
        self.Parameters = Parameters(self._parameters, self._synonyms)
        if params:
            for key,v in params.items():
                try:
                    self.Parameters[key].on(v)
                except TypeError:
                    self.Parameters[key].on()

class CommandLineApplication(Application):
    """ Generic class for controlling command line applications 
    """

    _input_handler = '_input_as_string'
    _suppress_stderr = False
    _suppress_stdout = False
    _working_dir = None

    def __init__(self,params=None,InputHandler=None,SuppressStderr=None,\
        SuppressStdout=None,WorkingDir=None,TmpDir='/tmp', \
        TmpNameLen=20, HALT_EXEC=False):
        """ Initialize the CommandLineApplication object
        
            params: a dictionary mapping the Parameter id or synonym to its
                value (or None for FlagParameters or MixedParameters in flag
                mode) for Parameters that should be turned on
            InputHandler: this is the method to be run on data when it is
                passed into call. This should be a string containing the
                method name. The default is _input_as_string which casts data
                to a string before appending it to the command line argument
            SuppressStderr: if set to True, will route standard error to
                /dev/null, False by default
            SuppressStdout: if set to True, will route standard out to
                /dev/null, False by default
            WorkingDir: the directory where you want the application to run,
                default is the current working directory, but is useful to 
                change in cases where the program being run creates output
                to its current working directory and you either don't want
                it to end up where you are running the program, or the user 
                running the script doesn't have write access to the current 
                working directory
                WARNING: WorkingDir MUST be an absolute path!
            TmpDir: the directory where temp files will be created, /tmp
                by default 
            TmpNameLen: the length of the temp file name
            HALT_EXEC: if True, raises exception w/ command output just
            before execution, doesn't clean up temp files. Default False.
        """
        # set attributes to parameter that was passed in or class default
        if InputHandler is not None:
            self.InputHandler = InputHandler
        else:
            self.InputHandler = self._input_handler
        if SuppressStderr is not None:
            self.SuppressStderr = SuppressStderr
        else:
            self.SuppressStderr = self._suppress_stderr
        if SuppressStdout is not None:
            self.SuppressStdout = SuppressStdout
        else:
            self.SuppressStdout = self._suppress_stdout
        if WorkingDir is not None:
            working_dir = WorkingDir
        else:
            working_dir = self._working_dir or getcwd()
        self.WorkingDir = FilePath(working_dir)
        self.TmpDir = FilePath(TmpDir)
        self.TmpNameLen = TmpNameLen
        self.HaltExec = HALT_EXEC
        #===========================
        #try: 
        #    mkdir(self.WorkingDir)
        #except OSError:
            # Directory already exists
        #    pass 
        #===========================
        # create a variable to hold the name of the file being used as
        # input to the application. this is important especially when 
        # you are using an input handler which creates a temporary file
        # and the output filenames are based on the input filenames
        self._input_filename = None
        
        super(CommandLineApplication,self).__init__(params=params)
    
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

        return result
   
    def _input_as_string(self,data):
        """ Return data as a string """
        return str(data)

    def _input_as_multiline_string(self, data):
        """Write a multiline string to a temp file and return the filename.

            data: a multiline string to be written to a file.

           * Note: the result will be the filename as a FilePath object 
            (which is a string subclass).

        """
        filename = self._input_filename = \
            FilePath(self.getTmpFilename(self.TmpDir))
        data_file = open(filename,'w')
        data_file.write(data)
        data_file.close()
        return filename
   
    def _input_as_lines(self,data):
        """ Write a seq of lines to a temp file and return the filename string
        
            data: a sequence to be written to a file, each element of the 
                sequence will compose a line in the file
           * Note: the result will be the filename as a FilePath object 
            (which is a string subclass).

           * Note: '\n' will be stripped off the end of each sequence element
                before writing to a file in order to avoid multiple new lines
                accidentally be written to a file
        """
        filename = self._input_filename = \
            FilePath(self.getTmpFilename(self.TmpDir))
        filename = FilePath(filename)
        data_file = open(filename,'w')
        data_to_file = '\n'.join([str(d).strip('\n') for d in data])
        data_file.write(data_to_file)
        data_file.close()
        return filename

    def _input_as_path(self,data):
        """ Return data as string with the path wrapped in quotes
            
            data: path or filename, most likely as a string
           
            * Note: the result will be the filename as a FilePath object 
            (which is a string subclass).

        """
        return FilePath(data)

    def _input_as_paths(self,data):
        """ Return data as a space delimited string with each path quoted
            
            data: paths or filenames, most likely as a list of 
             strings

        """
        return self._command_delimiter.join(\
            map(str,map(self._input_as_path,data)))

    def _absolute(self,path):
        """ Convert a filename to an absolute path """
        path = FilePath(path)
        if isabs(path):
            return path
        else:
            # these are both Path objects, so joining with + is acceptable
            return self.WorkingDir + path
 
    def _get_base_command(self):
        """ Returns the full command string 

            input_arg: the argument to the command which represents the input 
                to the program, this will be a string, either 
                representing input or a filename to get input from
         tI"""
        command_parts = []
        # Append a change directory to the beginning of the command to change 
        # to self.WorkingDir before running the command
        # WorkingDir should be in quotes -- filenames might contain spaces
        cd_command = ''.join(['cd ',str(self.WorkingDir),';'])
        if self._command is None:
            raise ApplicationError, '_command has not been set.'
        command = self._command
        parameters = self.Parameters
        synonyms = self._synonyms
        
        command_parts.append(cd_command)
        command_parts.append(command)
        command_parts.append(self._command_delimiter.join(filter(\
            None,(map(str,parameters.values())))))
      
        return self._command_delimiter.join(command_parts).strip()
    
    BaseCommand = property(_get_base_command)
    
    def _get_WorkingDir(self):
        """Gets the working directory"""
        return self._curr_working_dir

    def _set_WorkingDir(self,path):
        """Sets the working directory
        
        Appends a slash to the end of path
        The reasoning behind this is that the user may or may not pass
        in a path with a '/' at the end. Since having multiple
        '/' at the end doesn't hurt anything, it's convienient to 
        be able to rely on it, and not have to check for it
        """
        self._curr_working_dir = FilePath(path) + '/'
        try:
            mkdir(self.WorkingDir)
        except OSError:
            # Directory already exists
            pass 

    WorkingDir = property(_get_WorkingDir,_set_WorkingDir)

    
    def _accept_exit_status(self,exit_status):
        """ Return False to raise an error due to exit_status of applciation
        
            This method should be overwritten if you'd like to raise an error
            based on certain exit statuses of the application that was run. The
            default is that no value of exit_status will raise an error.
        """
        return True

    def _get_result_paths(self,data):
        """ Return a dict of ResultPath objects representing all possible output
            
            This method should be overwritten if the application creates output
            other than stdout and stderr.  This dictionary will have keys based
            on the name that you'd like to access the file by in the 
            CommandLineAppResult object that will be created, and the values
            which are ResultPath objects. For an example of how this should be
            written see the rnaview or vienna_package classes.
            WARNING: be sure that the path that you give a file is accurate
                from any directory where the program could be running. For that
                reason, absolute paths are very good. Relative paths can also be
                used as long as you are careful. For cases where the application
                leaves files in the current working directory, you should append
                self.WorkingDir to the beginning of the file name.
                It would be a very bad idea to just use a file name as the path,
                in some cases that you might not be testing for.
        """
        return {}


    def getTmpFilename(self, tmp_dir="/tmp",prefix='tmp',suffix='.txt',\
        include_class_id=False):
        """ Return a temp filename

            tmp_dir: path for temp file
            prefix: text to append to start of file name
            suffix: text to append to end of file name
            include_class_id: if True, will append a class identifier (built
             from the class name) to the filename following prefix. This is 
             False by default b/c there is some string processing overhead
             in getting the class name. This will probably be most useful for
             testing: if temp files are being left behind by tests, you can
             turn this on in here (temporarily) to find out which tests are
             leaving the temp files.
        """

        # check not none
        if not tmp_dir:
            tmp_dir = self.TmpDir
        # if not current directory, append "/" if not already on path
        elif not tmp_dir.endswith("/"):
            tmp_dir += "/"

        if include_class_id:
            # Append the classname to the prefix from the class name 
            # so any problematic temp files can be associated with 
            # the class that created them. This should be especially 
            # useful for testing, but is turned off by default to
            # avoid the string-parsing overhead.
            class_id = str(self.__class__())
            prefix = ''.join([prefix,\
             class_id[class_id.rindex('.')+1:class_id.index(' ')]])
        
        try:
            mkdir(tmp_dir)
        except OSError:
            # Directory already exists
            pass
        # note: it is OK to join FilePath objects with + 
        return FilePath(tmp_dir) + FilePath(prefix) + \
            FilePath(''.join([choice(_all_chars) \
             for i in range(self.TmpNameLen)])) +\
            FilePath(suffix)

class ParameterEnumerator:
    """Enumerates a set of parameters for a specific application"""
    def __init__(self, Application, Parameters, AlwaysOn=None):
        """Initialize the ParameterEnumerator

        Application : A CommandLineApplication subclass
        Parameters  : A dict keyed by the application paramter, value by
                      the range of parameters to enumerate over. For 
                      FlagParameters, unless specified in AlwaysOn, the value
                      will cycle between True/False (on/off). For 
                      MixedParameters, include [None] specifically to utilize
                      flag functionality.
        AlwaysOn    : List of parameters that will always be on

        Parameters is checked against the applications known parameters, but
        only performed superficially: only keys are validated. AlwaysOn
        values must have entries within Parameters.

        NOTE: If the parameter is not specified in AlwaysOn, a False value
        is appended so that the parameter can be turned off. Multiple False
        states for a parameter will result if False is specified without
        adding the parameter to AlwaysOn. If a parameter has a default value,
        then that parameter is implicitly always on.
        """
        self._app_params = Application._parameters
       
        # Validate Parameters
        for k in Parameters:
            if k not in self._app_params:
                raise ValueError, "Parameter %s not present in app" % k
            values = Parameters[k]

            # Make sure the parameters are specified in a list
            if not isinstance(values, list):
                Parameters[k] = [values]
        _my_params = Parameters
        
        # Validate AlwaysOn
        for k in AlwaysOn:
            if k not in _my_params:
                raise ValueError, "AlwaysOn %s not passed with Parameters" % k

        # Append "off states" to relevant parameters
        for k in set(_my_params.keys()) - set(AlwaysOn):
            _my_params[k] = _my_params[k] + [False]

        # Create seperate key/value lists preserving index relation
        self._keys, self._values = zip(*sorted(_my_params.items()))

        # Construct generator
        self._generator = self._get_combinations()

    def _get_combinations(self):
        """Enumerates all possible combinations of parameters"""
        num_items = [len(i) - 1 for i in self._values]
        state = zeros(len(self._values), dtype=int)
        finished = array(num_items, dtype=int)

        yield self._make_app_params(state)

        while True:
            if state[-1] != num_items[-1]:
                state[-1] += 1
                yield self._make_app_params(state)
            else:
                incrementable = nonzero(state != finished)[0]
                if not len(incrementable):
                    raise StopIteration
                rightmost = max(incrementable)
                state[rightmost] += 1
                state[rightmost+1:] = 0
                yield self._make_app_params(state)

    def _make_app_params(self, state):
        """Returns an app's param dict with values set as described by state"""
        app_params = self._app_params.copy()
        for key, values, state_idx in zip(self._keys, self._values, state):
            val = values[state_idx]
            if val is False:
                app_params[key].off()
            elif val is True:
                app_params[key].on()
            else:
                app_params[key].on(val)
        return app_params
   
    def __iter__(self):
        return self

    def next(self):
        return self._generator.next()

def get_tmp_filename(tmp_dir="/tmp", prefix="tmp", suffix=".txt"):
    """
    Generate temp filename
    """
    # check not none
    if not tmp_dir:
        tmp_dir = ""
    # if not current directory, append "/" if not already on path
    elif not tmp_dir.endswith("/"):
        tmp_dir += "/"

    chars = "abcdefghigklmnopqrstuvwxyz"
    picks = chars + chars.upper() + "0123456790"
    return FilePath(tmp_dir) + FilePath(prefix) +\
        FilePath("%s%s" % \
        (''.join([choice(picks) for i in range(20)]),suffix))

def guess_input_handler(seqs, add_seq_names=False):
    """Returns the name of the input handler for seqs."""
    if isinstance(seqs, str):
        if '\n' in seqs:    #can't be a filename...
            ih = '_input_as_multiline_string'
        else:               #assume it was a filename
            ih = '_input_as_string'
            #Uncommenting the next line causes errors in muscle tests - Micah?
            #ih = '_input_as_path'
    elif add_seq_names:
        ih = '_input_as_seqs'
    else:
        ih = '_input_as_lines'
    return ih


