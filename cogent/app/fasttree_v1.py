#!/usr/bin/env python
"""Application controller for FastTree v1.0"""

from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, FilePath, system, \
       CommandLineAppResult, ResultPath, remove, ApplicationError
from cogent.core.tree import PhyloNode
from cogent.parse.tree import DndParser
from cogent.core.moltype import DNA, RNA, PROTEIN

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class FastTree(CommandLineApplication):
    """FastTree application Controller"""

    _command = 'FastTree'
    _input_handler = '_input_as_multiline_string'
    _parameters = {
            '-quiet':FlagParameter('-',Name='quiet'),
            '-boot':ValuedParameter('-',Delimiter=' ',Name='boot'),
            '-seed':ValuedParameter('-',Delimiter=' ',Name='seed'),
            '-nni':ValuedParameter('-',Delimiter=' ',Name='nni'),
            '-slow':FlagParameter('-',Name='slow'),
            '-fastest':FlagParameter('-',Name='fastest'),
            '-top':FlagParameter('-',Name='top'),
            '-notop':FlagParameter('-',Name='notop'),
            '-topm':ValuedParameter('-',Delimiter=' ',Name='topm'),
            '-close':ValuedParameter('-',Delimiter=' ',Name='close'),
            '-refresh':ValuedParameter('-',Delimiter=' ',Name='refresh'),
            '-matrix':ValuedParameter('-',Delimiter=' ',Name='matrix'),
            '-nomatrix':FlagParameter('-',Name='nomatrix'),
            '-nj':FlagParameter('-',Name='nj'),
            '-bionj':FlagParameter('-',Name='bionj'),
            '-nt':FlagParameter('-',Name='nt'),
            '-n':ValuedParameter('-',Delimiter=' ',Name='n')}

    #FastTree [-quiet] [-boot 1000] [-seed 1253] [-nni 10] [-slow | -fastest]
    #      [-top | -notop] [-topm 1.0 [-close 0.75] [-refresh 0.8]]
    #      [-matrix Matrix | -nomatrix] [-nj | -bionj]
    #      [-nt] [-n 100] [alignment] > newick_tree

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

        result =  CommandLineAppResult(out,err,exit_status,\
            result_paths=self._get_result_paths(data))

        # Clean up the input file if one was created
        if remove_tmp:
            if self._input_filename:
                remove(self._input_filename)
                self._input_filename = None

        return result

    def _get_result_paths(self, data):
        result = {}
        result['Tree'] = ResultPath(Path=self._outfile)
        return result
 
def build_tree_from_alignment(aln, moltype, best_tree=False, params=None):
    """Returns a tree from alignment
    
    Will check MolType of aln object
    """
    if params is None:
        params = {}

    if moltype == DNA or moltype == RNA:
        params['-nt'] = True
    elif moltype == PROTEIN:
        params['-nt'] = False
    else:
        raise ValueError, \
                "FastTree does not support moltype: %s" % moltype.label

    app = FastTree(params=params)
    
    if best_tree:
        raise NotImplementedError, "best_tree not implemented yet"
    result = app(aln.toFasta())
    tree = DndParser(result['Tree'].read(), constructor=PhyloNode)
    return tree
   
