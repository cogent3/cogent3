#!/usr/bin/env python
"""Application controller for guppy 1.1"""

from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, FilePath, system, \
       CommandLineAppResult, ResultPath, remove, ApplicationError
from cogent.core.alignment import Alignment
from os.path import splitext,split,join
from os import listdir
from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode
       
__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Prototype"

class Guppy(CommandLineApplication):
    """guppy Application Controller
    """

    _command = 'guppy'
    _input_handler = '_input_as_multiline_string'
    _parameters = {
    #visualizations
        # makes trees with edges fattened in proportion to the number of reads
        'fat': FlagParameter('', Name='fat'),
        
        # maps an an arbitrary vector of the correct length to the tree
        'heat': FlagParameter('', Name='heat'),
        
        # writes a taxonomically annotated reference tree and an induced 
        # taxonomic tree
        'ref_tree': FlagParameter('', Name='ref_tree'),
        
        # makes one tree for each query sequence, showing uncertainty
        'sing': FlagParameter('', Name='sing'),
        
        # makes a tree with each of the reads represented as a pendant edge
        'tog': FlagParameter('', Name='tog'),
      
      #statistical comparison
        # draws the barycenter of a placement collection on the reference tree
        'bary': FlagParameter('', Name='bary'),
        
        # makes a phyloXML tree showing the bootstrap values
        'bootviz': FlagParameter('', Name='bootviz'),
        
        # calculates the EDPL uncertainty values for a collection of pqueries
        'edpl': FlagParameter('', Name='edpl'),
        
        # calculates the Kantorovich-Rubinstein distance and corresponding 
        # p-values
        'kr': FlagParameter('', Name='kr'),
        
        # makes a heat tree
        'kr_heat': FlagParameter('', Name='kr_heat'),
        
        # performs edge principal components
        'pca': FlagParameter('', Name='pca'),
        
        # writes out differences of masses for the splits of the tree
        'splitify': FlagParameter('', Name='splitify'),
        
        # performs squash clustering
        'squash': FlagParameter('', Name='squash'),
      
      #classification
        # outputs classification information in a tabular or SQLite format
        'classify': FlagParameter('', Name='classify'),
      
      #utilities
        # check a reference package
        'check_refpkg': FlagParameter('', Name='check_refpkg'),
        
        # splits apart placements with multiplicity, undoing a round procedure
        'demulti': FlagParameter('', Name='demulti'),
        
        # prints out a pairwise distance matrix between the edges
        'distmat': FlagParameter('', Name='distmat'),
        
        # filters one or more placefiles by placement name
        'filter': FlagParameter('', Name='filter'),
        
        # writes the number of leaves of the reference tree and the number of 
        # pqueries
        'info': FlagParameter('', Name='info'),
        
        # merges placefiles together
        'merge': FlagParameter('', Name='merge'),
        
        # restores duplicates to deduped placefiles
        'redup': FlagParameter('', Name='redup'),
        
        # clusters the placements by rounding branch lengths
        'round': FlagParameter('', Name='round'),
        
        # makes SQL enabling taxonomic querying of placement results
        'taxtable': FlagParameter('', Name='taxtable'),
        
        # converts old-style .place files to .json placement files
        'to_json': FlagParameter('', Name='to_json'),
        
        # Run the provided batch file of guppy commands
        'batch': FlagParameter('--', Name='batch'),
        
        # Print version and exit
        'version': FlagParameter('--', Name='version'),
        
        # Print a list of the available commands.
        'cmds': FlagParameter('--', Name='cmds'),
        
        # Display this list of options
        '--help': FlagParameter('--', Name='help'),
        
        # Display this list of options
        '-help': FlagParameter('-', Name='help'),
    }
 
    def __call__(self,fname=None, remove_tmp=True):
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
        
        if fname is None:
            self._input_filename = ''
        else:
            self._input_filename = fname
        
        # Build up the command, consisting of a BaseCommand followed by
        # input and output (file) specifications
        command = self._command_delimiter.join(filter(None,\
            [self.BaseCommand,str(self._input_filename),'>',str(outfile),'2>',\
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
            result_paths=self._get_result_paths())
        
        return result

    def _get_result_paths(self):
        basepath,basename=split(splitext(self._input_filename)[0])
        outfile_list=listdir(split(self._input_filename)[0])
        result = {}
        for i in outfile_list:
            if i.startswith(basename) and not i.endswith('.json'):
                result['result'] = ResultPath(Path=join(basepath,i))
                
        return result
    
def build_tree_from_json_using_params(fname,output_dir='/tmp/',params={}):
    """Returns a tree from Alignment object aln.

    aln: an xxx.Alignment object, or data that can be used to build one.

    moltype: cogent.core.moltype.MolType object

    params: dict of parameters to pass in to the RAxML app controller.

    The result will be an xxx.Alignment object, or None if tree fails.
    """

    # convert aln to fasta in case it is not already a fasta file
    
    ih = '_input_as_multiline_string'    

    guppy_app = Guppy(params=params,
                      InputHandler=ih,
                      WorkingDir=output_dir,
                      SuppressStderr=True,
                      SuppressStdout=True)
                      
    guppy_result = guppy_app(fname)
    
    new_tree=guppy_result['result'].read()

    tree = DndParser(new_tree, constructor=PhyloNode)
    
    guppy_result.cleanUp()

    return tree



    
