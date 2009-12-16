#!/usr/bin/env python
"""Provides an application controller for the commandline version of
mothur Version 1.6.0
"""


from __future__ import with_statement
from os import path, getcwd, mkdir, remove
from shutil import copyfile
from subprocess import Popen
from cogent.app.parameters import ValuedParameter
from cogent.app.util import CommandLineApplication, ResultPath, \
    CommandLineAppResult, ApplicationError
from cogent.parse.mothur import parse_otu_list


__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"


class Mothur(CommandLineApplication):
    """Mothur application controller
    """
    _options = {
        # Clustering algorithm.  Choices are furthest, nearest, and
        # average
        'method': ValuedParameter(
            Name='method', Value='furthest', Delimiter='=', Prefix=''),
        # Cutoff distance for the distance matrix
        'cutoff': ValuedParameter(
            Name='cutoff', Value=None, Delimiter='=', Prefix=''),
        # Minimum pairwise distance to consider for clustering
        'precision': ValuedParameter(
            Name='precision', Value=None, Delimiter='=', Prefix=''),
        }
    _parameters = {}
    _parameters.update(_options)
    _input_handler = '_input_as_multiline_string'
    _command = 'mothur'
    
    def __init__(self, params=None, InputHandler=None, SuppressStderr=None,
                 SuppressStdout=None, WorkingDir=None, TmpDir='/tmp',
                 TmpNameLen=20, HALT_EXEC=False):
        """Initialize a Mothur application controller
        
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
        super(Mothur, self).__init__(
            params=params, InputHandler=InputHandler, 
            SuppressStderr=SuppressStderr, SuppressStdout=SuppressStdout,
            WorkingDir='', TmpDir='', TmpNameLen=TmpNameLen, 
            HALT_EXEC=HALT_EXEC)
        # Prevent self.WorkingDir from being explicitly cast as a
        # FilePath object.  This behavior does not seem necessary in
        # the parent's __init__() method, since the casting is
        # repeated in _set_WorkingDir().
        if WorkingDir is not None:
            working_dir = WorkingDir
        else:
            working_dir = self._working_dir or getcwd()
        self.WorkingDir = working_dir
        self.TmpDir = TmpDir

    @staticmethod
    def getHelp():
        """Returns link to online manual"""
        help = (
            'See manual, available on the MOTHUR wiki:\n'
            'http://schloss.micro.umass.edu/mothur/'
            )
        return help

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
        # Process the input data.  Input filepath is stored in
        # self._input_filename
        getattr(self, self.InputHandler)(data)

        if self.SuppressStdout:
            outfile = None
        else:
            outfile = open(self.getTmpFilename(self.TmpDir), 'w')
        if self.SuppressStderr:
            errfile = None
        else:
            errfile = open(self.getTmpFilename(self.TmpDir), 'w')

        args = [self._command, self._compile_mothur_script()]
        process = Popen(
            args, stdout=outfile, stderr=errfile, cwd=self.WorkingDir)
        exit_status = process.wait()
        if not self._accept_exit_status(exit_status):
            raise ApplicationError(
                'Unacceptable application exit status: %s, command: %s' % \
                    (exit_status, args))

        if outfile is not None:
            outfile.seek(0)
        if errfile is not None:
            errfile.seek(0)
        result = CommandLineAppResult(
            outfile, errfile, exit_status, result_paths=self._get_result_paths())

        # Clean up the input file if one was created
        if remove_tmp:
            if self._input_filename:
                remove(self._input_filename)
                self._input_filename = None

        return result

    def _accept_exit_status(self, status):
        return int(status) == 0

    def _compile_mothur_script(self):
        """Returns a Mothur batch script as a string"""
        def format_opts(*opts):
            """Formats a series of options for a Mothur script"""
            return ', '.join(filter(None, map(str, opts)))
        vars = {
            'in': self._input_filename,
            'unique': self.__derive_unique_path(),
            'dist': self.__derive_dist_path(),
            'names': self.__derive_names_path(),
            'cluster_opts' : format_opts(
                self.Parameters['method'],
                self.Parameters['cutoff'],
                self.Parameters['precision'],
                ),
            }
        script = (
            '#'
            'unique.seqs(fasta=%(in)s); '
            'dist.seqs(fasta=%(unique)s); '
            'read.dist(column=%(dist)s, name=%(names)s); '
            'cluster(%(cluster_opts)s)' % vars
            )
        return script

    def _get_result_paths(self):
        paths = {
            'distance matrix': self.__derive_dist_path(),
            'otu list': self.__derive_list_path(),
            'rank abundance': self.__derive_rank_abundance_path(),
            'species abundance': self.__derive_species_abundance_path(),
            'unique names': self.__derive_names_path(),
            'unique seqs': self.__derive_unique_path(),
            'log': self.__derive_log_path(),
            }
        return dict([(k, ResultPath(v)) for (k,v) in paths.items()])

    # Methods to derive/guess output pathnames produced by MOTHUR.
    # TODO: test for input files that do not have a filetype extension

    def __derive_log_path(self):
        """Guess logfile path produced by Mothur"""
        dir = path.dirname(self._input_filename)
        return path.join(dir, 'mothur.logFile')

    def __derive_unique_path(self):
        """Guess unique sequences path produced by Mothur"""
        base, ext = path.splitext(self._input_filename)
        return '%s.unique%s' % (base, ext)

    def __derive_dist_path(self):
        """Guess distance matrix path produced by Mothur"""
        base, ext = path.splitext(self._input_filename)
        return '%s.unique.dist' % base

    def __derive_names_path(self):
        """Guess unique names file path produced by Mothur"""
        base, ext = path.splitext(self._input_filename)
        return '%s.names' % base

    def __get_method_abbrev(self):
        """Abbreviated form of clustering method parameter.

        Used to guess output filenames for MOTHUR.
        """
        abbrevs = {
            'furthest': 'fn',
            'nearest':  'nn',
            'average':  'an',
        }
        if self.Parameters['method'].isOn():
            method = self.Parameters['method'].Value
        else:
            method = self.Parameters['method'].Default
        return abbrevs[method]

    def __derive_list_path(self):
        """Guess otu list file path produced by Mothur"""
        base, ext = path.splitext(self._input_filename)
        return '%s.unique.%s.list' % (base, self.__get_method_abbrev())
        
    def __derive_rank_abundance_path(self):
        """Guess rank abundance file path produced by Mothur"""
        base, ext = path.splitext(self._input_filename)
        return '%s.unique.%s.rabund' % (base, self.__get_method_abbrev())
        
    def __derive_species_abundance_path(self):
        """Guess species abundance file path produced by Mothur"""
        base, ext = path.splitext(self._input_filename)
        return '%s.unique.%s.sabund' % (base, self.__get_method_abbrev())

    def getTmpFilename(self, tmp_dir='/tmp', prefix='tmp', suffix='.txt'):
        """Returns a temporary filename

        Similar interface to tempfile.mktmp()
        """
        # Override to change default constructor to str(). FilePath
        # objects muck up the Mothur script.
        return super(Mothur, self).getTmpFilename(
            tmp_dir=tmp_dir, prefix=prefix, suffix=suffix, 
            result_constructor=str)

    # Temporary input file needs to be in the working directory, so we
    # override all input handlers.

    def _input_as_multiline_string(self, data):
        """Write multiline string to temp file, return filename

        data: a multiline string to be written to a file.
        """
        self._input_filename = self.getTmpFilename(
            self.WorkingDir, suffix='.fasta')
        with open(self._input_filename, 'w') as f:
            f.write(data)
        return self._input_filename

    def _input_as_lines(self, data):
        """Write sequence of lines to temp file, return filename
        
        data: a sequence to be written to a file, each element of the
            sequence will compose a line in the file

        * Note: '\n' will be stripped off the end of each sequence
            element before writing to a file in order to avoid
            multiple new lines accidentally be written to a file
        """
        self._input_filename = self.getTmpFilename(
            self.WorkingDir, suffix='.fasta')
        with open(self._input_filename, 'w') as f:
            # Use lazy iteration instead of list comprehension to
            # prevent reading entire file into memory
            for line in data:
                f.write(str(line).strip('\n'))
                f.write('\n')
        return self._input_filename

    def _input_as_path(self, data):
        """Copys the provided file to WorkingDir and returns the new filename

        data: path or filename
        """
        self._input_filename = self.getTmpFilename(
            self.WorkingDir, suffix='.fasta')
        copyfile(data, self._input_filename)
        return self._input_filename

    def _input_as_paths(self, data):
        raise NotImplementedError('Not applicable for MOTHUR controller.')

    def _input_as_string(self, data):
        raise NotImplementedError('Not applicable for MOTHUR controller.')

    # FilePath objects muck up the Mothur script, so we override the
    # property methods for self.WorkingDir

    def _get_WorkingDir(self):
        """Gets the working directory"""
        return self._curr_working_dir

    def _set_WorkingDir(self, path):
        """Sets the working directory
        """
        self._curr_working_dir = path
        try:
            mkdir(self.WorkingDir)
        except OSError:
            # Directory already exists
            pass

    WorkingDir = property(_get_WorkingDir, _set_WorkingDir)

def mothur_from_file(file):
    app = Mothur(InputHandler='_input_as_lines')
    result = app(file)
    # Force evaluation, so we can safely clean up files
    otus = list(parse_otu_list(result['otu list']))
    result.cleanUp()
    return otus

