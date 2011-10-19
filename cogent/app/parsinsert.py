#!/usr/bin/env python
"""Application controller for ParsInsert

designed for ParsInsert v1.03 """

from cogent.app.parameters import ValuedParameter, FlagParameter, \
       MixedParameter
from cogent.app.util import CommandLineApplication, FilePath, system, \
       CommandLineAppResult, ResultPath, remove, ApplicationError
from cogent.core.tree import PhyloNode
from cogent.parse.tree import DndParser
from cogent.core.moltype import DNA, RNA, PROTEIN
from cogent.core.alignment import SequenceCollection
from os.path import splitext, join
__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"

class ParsInsert(CommandLineApplication):
    """ParsInsert application Controller"""

    _command = 'ParsInsert'
    _input_handler = '_input_as_multiline_string'
    _parameters = {
                    # read mask from this file
                    '-m':ValuedParameter('-',Name='m',Delimiter=' '),
                    
                    # read core tree sequences from this file
                    '-s':ValuedParameter('-',Name='s',Delimiter=' '),
                    
                    # read core tree from this file
                    '-t':ValuedParameter('-',Name='t',Delimiter=' '),

                    # read core tree taxomony from this file 
                    '-x':ValuedParameter('-',Name='x',Delimiter=' '),

                    # output taxonomy for each insert sequence to this file
                    '-o':ValuedParameter('-',Name='o',Delimiter=' '),

                    # create log file
                    '-l':ValuedParameter('-',Name='l',Delimiter=' '),
                    
                    # number of best matches to display
                    '-n':ValuedParameter('-',Name='n',Delimiter=' '),

                    #percent threshold cutoff
                    '-c':ValuedParameter('-',Name='c',Delimiter=' '),
                   }

    def __call__(self,data=None, remove_tmp=True):
        """Run the application with the specified kwargs on data
        
            data: anything that can be cast into a string or written out to
                a file. Usually either a list of things or a single string or 
                number. input_handler will be called on this data before it 
                is passed as part of the command-line argument, so by creating
                your own input handlers you can customize what kind of data
                you want your application to accept

            remove_tmp: if True, removes tmp files

            NOTE: Override of the base class to handle redirected output
        """
        input_handler = self.InputHandler
        suppress_stdout = self.SuppressStdout
        suppress_stderr = self.SuppressStderr

        outfile = self.getTmpFilename(self.TmpDir)
        self._outfile = outfile
        
        if suppress_stderr:
            errfile = FilePath('/dev/null')
        else:
            errfile = FilePath(self.getTmpFilename(self.TmpDir))
        
        if data is None:
            input_arg = ''
        else:
            input_arg = getattr(self,input_handler)(data)
        
        self._tree_fname=input_arg
        
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

        out = open(outfile,"r")
        err = None
        if not suppress_stderr:
            err = open(errfile,"r")
        
        # capture the command-call error, since the normal error only says file
        # not found
        try:
            result =  CommandLineAppResult(out,err,exit_status,\
                result_paths=self._get_result_paths())
        except ApplicationError:
            raise ApplicationError, 'Pplacer failed to produce an output file due to the following error: \n\n%s ' % open(errfile).read()

        # Clean up the input file if one was created
        if remove_tmp:
            if self._input_filename:
                remove(self._input_filename)
                self._input_filename = None
        
        return result

    def _get_result_paths(self):
        result = {}
        result['Tree'] = ResultPath(Path=splitext(self._tree_fname)[0]+'.tree')
        return result

def test_build_tree_from_alignment_using_params(aln, moltype, params={}):
    """Returns a tree from placement of sequences
    """

    # verify seqs are DNA
    if moltype != DNA:
        raise ValueError, \
                "ParsInsert does not support moltype: %s" % moltype.label

    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = aln.getIntMap()
    
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)

    app = ParsInsert(params=params)
    result = app(aln.toFasta())
    
    # parse tree
    tree = DndParser(result['Tree'].read(), constructor=PhyloNode)
    
    # cleanup files
    result.cleanUp()
    
    return tree
    
