#!/usr/bin/env python
"""Application controller for rdp_classifier-2.0
"""

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.4.0.dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"


from os import remove
from cogent.app.parameters import Parameter, ValuedParameter, Parameters
from cogent.app.util import CommandLineApplication, CommandLineAppResult, \
    FilePath, ResultPath, get_tmp_filename, guess_input_handler, system

class RdpClassifier(CommandLineApplication):
    """RDP Classifier application controller

    The RDP Classifier program is distributed as a java archive (.jar)
    file, so an option '-jar' is provided to optionally specify the
    full path to the archive.  By default, the file
    'rdp_classifier-2.0.jar' is called in the current directory.

    The RDP Classifier often requires memory in excess of Java's
    default 64M. To correct this situation, the authors recommend
    increasing the maximum heap size for the java virtual machine.  An
    option '-Xmx' (default 1000M) is provided for this purpose.
    Details on this option may be found at
    http://java.sun.com/j2se/1.5.0/docs/tooldocs/solaris/java.html

    The classifier may optionally use a custom training set.  This
    ability is not yet implemented.
    """
    _input_handler = '_input_as_multiline_string'
    _command = "rdp_classifier-2.0.jar"
    _options ={}
    _java_vm_parameters = {
        # JAVA VM OPTIONS 
        #
        # These are extracted from the parameters dict to construct
        # the command for the Java virtual machine.
        #
        # Maximum heap size for JVM.
        '-Xmx': ValuedParameter('-', Name='Xmx', Delimiter='', Value='1000m'),
        #
        # Location of the java archive file.
        '-jar': ValuedParameter('-', Name='jar', Delimiter=' ', Value=_command, IsPath=True),
        }
    _parameters = {}
    _parameters.update(_options)
    _parameters.update(_java_vm_parameters)

    def getHelp(self):
        """Returns documentation string"""
        # Summary paragraph copied from rdp_classifier-2.0, which is
        # licensed under the GPL 2.0 and Copyright 2008 Michigan State
        # University Board of Trustees
        help_str =\
        """
        Ribosomal Database Project - Classifier
        http://rdp.cme.msu.edu/classifier/

        The RDP Classifier is a naive Bayesian classifier which was
        developed to provide rapid taxonomic placement based on rRNA
        sequence data. The RDP Classifier can rapidly and accurately
        classify bacterial 16s rRNA sequences into the new
        higher-order taxonomy proposed by Bergey's Trust. It provides
        taxonomic assignments from domain to genus, with confidence
        estimates for each assignment. The RDP Classifier is not
        limited to using the bacterial taxonomy proposed by the
        Bergey's editors. It worked equally well when trained on the
        NCBI taxonomy. The RDP Classifier likely can be adapted to
        additional phylogenetically coherent bacterial taxonomies.

        The following paper should be cited if this resource is used:

        Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007.
        Naive Bayesian Classifier for Rapid Assignment of rRNA
        Sequences into the New Bacterial Taxonomy.  Appl Environ
        Microbiol. 73(16):5261-7.
        """
        return help_str

    # The __call__ method has been copied directly from superclass,
    # for the sake of a single minor change in the command list.  It
    # may be worth factoring out a _get_full_command method from
    # __call__ to avoid this.
    def __call__(self, data=None, remove_tmp=True):
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
            [self.BaseCommand, str(input_arg), str(outfile),'2>',\
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

    # Overridden to pull out JVM-specific command-line arguments.
    def _get_base_command(self):
        """Returns the base command plus command-line options.

        Does not include input file, output file, and training set.
        """
        # Repeated pattern; may be useful in superclass
        def commandline_join(words):
            return self._command_delimiter.join(map(str, words)).strip()

        # Copy method does not work for self.Parameters, so substitute
        # with simple key-by-key copy.
        def copy_by_key(parameters):
            result = {}
            for key in parameters:
                result[key] = parameters[key]
            return Parameters(result)

        def extract(dict, keys):
            result = {}
            for key in keys:
                if key in dict:
                    result[key] = dict.pop(key)
            return result

        # Divy up parameters using local copy of the Parameters class
        # attribute
        parameters = copy_by_key(self.Parameters)

        # Append a change directory to the beginning of the command to change 
        # to self.WorkingDir before running the command
        # WorkingDir should be in quotes -- filenames might contain spaces
        cd_command = ''.join(['cd ',str(self.WorkingDir),';'])

        jvm_command = "java"
        jvm_parameters = extract(parameters, ['-Xmx'])
        jvm_arguments = commandline_join(jvm_parameters.values())

        if self._command is None:
            raise ApplicationError, '_command has not been set.'

        if '-jar' in parameters:
            jar_command = parameters.pop('-jar')
        else:
            jar_command = "-jar %s" % self._command

        jar_arguments = commandline_join(parameters.values())

        result = commandline_join([cd_command, jvm_command, jvm_arguments,
                                   jar_command, jar_arguments])
        return result
    
    BaseCommand = property(_get_base_command)

