#!/usr/bin/env python
"""Provides an application controller for the commandline version of:
Clearcut v1.0.8
"""
from cogent.app.parameters import FlagParameter, ValuedParameter, \
    MixedParameter
from cogent.app.util import CommandLineApplication, ResultPath, get_tmp_filename
from cogent.core.alignment import SequenceCollection, Alignment
from cogent.core.moltype import DNA, RNA, PROTEIN
from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

MOLTYPE_MAP = {'DNA':'-D',
               'RNA':'-D',
               'PROTEIN':'-P',
               }

class Clearcut(CommandLineApplication):
    """ clearcut application controller 
   
    The parameters are organized by function to give some idea of how the 
    program works. However, no restrictions are put on any combinations 
    of parameters. Misuse of parameters can lead to errors or otherwise
    strange results.
    """
    #General options.
    _general = {\
        # --verbose.  More Output. (Default:OFF)
        '-v':FlagParameter('-',Name='v'),
        # --quiet.  Silent operation. (Default: ON)
        '-q':FlagParameter('-',Name='q',Value=True),
        # --seed=<seed>.  Explicitly set the PRNG seed to a specific value.
        '-s':ValuedParameter('-',Name='s',Delimiter='='),
        # --norandom.  Attempt joins deterministically.  (Default: OFF)
        '-r':FlagParameter('-',Name='r'),
        # --shuffle.  Randomly shuffle the distance matrix.  (Default: OFF)
        '-S':FlagParameter('-',Name='S'),
        #--neighbor.  Use traditional Neighbor-Joining algorithm. (Default: OFF)
        '-N':FlagParameter('-',Name='N'),
        
        }
         

    # Input file is distance matrix or alignment.  Default expects distance
    # matrix.  Output file is tree created by clearcut.
    _input = {\
        # --in=<infilename>.  Input file
        '--in':ValuedParameter('--',Name='in',Delimiter='=',IsPath=True),
        # --stdin.  Read input from STDIN. 
        '-I':FlagParameter('-',Name='I'),
        # --distance.  Input file is a distance matrix. (Default: ON)
        '-d':FlagParameter('-',Name='d',Value=True),
        # --alignment.  Input file is a set of aligned sequences.
        #     (Default: OFF)
        '-a':FlagParameter('-',Name='a'),
        # --DNA.  Input alignment are DNA sequences.
        '-D':FlagParameter('-',Name='D'),
        # --protein.  Input alignment are protein sequences.
        '-P':FlagParameter('-',Name='P'),
        }
  
  
    #Correction model for computing distance matrix (Default: NO Correction):
    _correction={\
        # --jukes.  Use Jukes-Cantor correction for computing distance matrix.
        '-j':FlagParameter('-',Name='j'),
        # --kimura.  Use Kimura correction for distance matrix.
        '-k':FlagParameter('-',Name='k'),
        
        }
    
    _output={\
        # --out=<outfilename>.  Output file
        '--out':ValuedParameter('--',Name='out',Delimiter='=',IsPath=True),
        # --stdout.  Output tree to STDOUT.
        '-O':FlagParameter('-',Name='O'),
        # --matrixout=<file> Output distance matrix to specified file.
        '-m':ValuedParameter('-',Name='m',Delimiter='='),
        # --ntrees=<n>.  Output n trees.  (Default: 1)
        '-n':ValuedParameter('-',Name='n',Delimiter='='),
        # --expblen.  Exponential notation for branch lengths. (Default: OFF)
        '-e':FlagParameter('-',Name='e'),
        # --expdist.  Exponential notation in distance output. (Default: OFF)
        '-E':FlagParameter('-',Name='E'),
        
        }

    
        #NOT SUPPORTED
        #'-h':FlagParameter('-','h'),       #Help
        #'-V':FlagParameter('-','V'),       #Version


    _parameters = {}
    _parameters.update(_general)
    _parameters.update(_input)
    _parameters.update(_correction)
    _parameters.update(_output)
 
    _command = 'clearcut'
   
    def getHelp(self):
        """Method that points to the Clearcut documentation."""
        help_str =\
        """
        See Clearcut homepage at:
        http://bioinformatics.hungry.com/clearcut/
        """
        return help_str
   
    def _input_as_multiline_string(self, data):
        """Writes data to tempfile and sets -infile parameter

        data -- list of lines
        """
        if data:
            self.Parameters['--in']\
                .on(super(Clearcut,self)._input_as_multiline_string(data))
        return ''

    def _input_as_lines(self,data):
        """Writes data to tempfile and sets -infile parameter

        data -- list of lines, ready to be written to file
        """
        if data:
            self.Parameters['--in']\
                .on(super(Clearcut,self)._input_as_lines(data))
        return ''

    def _input_as_seqs(self,data):
        """writes sequences to tempfile and sets -infile parameter

        data -- list of sequences

        Adds numbering to the sequences: >1, >2, etc.
        """
        lines = []
        for i,s in enumerate(data):
            #will number the sequences 1,2,3,etc.
            lines.append(''.join(['>',str(i+1)]))
            lines.append(s)
        return self._input_as_lines(lines)

    def _input_as_string(self,data):
        """Makes data the value of a specific parameter
    
        This method returns the empty string. The parameter will be printed
        automatically once set.
        """
        if data:
            self.Parameters['--in'].on(data)
        return ''
    
    def _tree_filename(self):
        """Return name of file containing the alignment
        
        prefix -- str, prefix of alignment file.
        """
        if self.Parameters['--out']:
            aln_filename = self._absolute(self.Parameters['--out'].Value)
        else:
            raise ValueError, "No tree output file specified."
        return aln_filename

    def _get_result_paths(self,data):
        """Return dict of {key: ResultPath}
        """
        result = {}
        if self.Parameters['--out'].isOn():
            out_name = self._tree_filename()
            result['Tree'] = ResultPath(Path=out_name,IsWritten=True)
        return result      


        
#SOME FUNCTIONS TO EXECUTE THE MOST COMMON TASKS


def align_unaligned_seqs(seqs, moltype, params=None):
    """Returns an Alignment object from seqs.

    seqs: SequenceCollection object, or data that can be used to build one.
    
    moltype: a MolType object.  DNA, RNA, or PROTEIN.

    params: dict of parameters to pass in to the Clearcut app controller.
    
    Result will be an Alignment object.
    """
    #Clearcut does not support alignment
    raise NotImplementedError, """Clearcut does not support alignment."""
    
def align_and_build_tree(seqs, moltype, best_tree=False, params={}):
    """Returns an alignment and a tree from Sequences object seqs.
    
    seqs: SequenceCollection object, or data that can be used to build one.
    
    best_tree: if True (default:False), uses a slower but more accurate
    algorithm to build the tree.

    params: dict of parameters to pass in to the Clearcut app controller.

    The result will be a tuple containing an Alignment object and a
    cogent.core.tree.PhyloNode object (or None for the alignment and/or tree
    if either fails).
    """
    #Clearcut does not support alignment
    raise NotImplementedError, """Clearcut does not support alignment."""
    
def build_tree_from_alignment(aln, moltype, best_tree=False, params={},\
    working_dir='/tmp'):
    """Returns a tree from Alignment object aln.

    aln: an cogent.core.alignment.Alignment object, or data that can be used
    to build one.
        -  Clearcut only accepts aligned sequences.  Alignment object used to
        handle unaligned sequences.
    
    moltype: a cogent.core.moltype object.
        - NOTE: If moltype = RNA, we must convert to DNA since Clearcut v1.0.8
        gives incorrect results if RNA is passed in.  'U' is treated as an 
        incorrect character and is excluded from distance calculations.

    best_tree: if True (default:False), uses a slower but more accurate
    algorithm to build the tree.

    params: dict of parameters to pass in to the Clearcut app controller.

    The result will be an cogent.core.tree.PhyloNode object, or None if tree
    fails.
    """
    params['--out'] = get_tmp_filename(working_dir)
    
    # Create instance of app controller, enable tree, disable alignment
    app = Clearcut(InputHandler='_input_as_multiline_string', params=params, \
                   WorkingDir=working_dir, SuppressStdout=True,\
                   SuppressStderr=True)
    #Input is an alignment
    app.Parameters['-a'].on()
    #Turn off input as distance matrix
    app.Parameters['-d'].off()
    
    #If moltype = RNA, we must convert to DNA.
    if moltype == RNA:
        moltype = DNA
    
    if best_tree:
        app.Parameters['-N'].on()
    
    #Turn on correct moltype
    moltype_string = moltype.label.upper()
    app.Parameters[MOLTYPE_MAP[moltype_string]].on()    

    # Setup mapping. Clearcut clips identifiers. We will need to remap them.
    # Clearcut only accepts aligned sequences.  Let Alignment object handle
    # unaligned sequences.
    seq_aln = Alignment(aln,MolType=moltype)
    #get int mapping
    int_map, int_keys = seq_aln.getIntMap()
    #create new Alignment object with int_map
    int_map = Alignment(int_map)

    # Collect result
    result = app(int_map.toFasta())
    
    # Build tree
    tree = DndParser(result['Tree'].read(), constructor=PhyloNode)
    for node in tree.tips():
        node.Name = int_keys[node.Name]

    # Clean up
    result.cleanUp()
    del(seq_aln, app, result, int_map, int_keys, params)

    return tree
    
def add_seqs_to_alignment(seqs, aln, params=None):
    """Returns an Alignment object from seqs and existing Alignment.

    seqs: an cogent.core.sequence.Sequence object, or data that can be used
    to build one.

    aln: an cogent.core.alignment.Alignment object, or data that can be used
    to build one

    params: dict of parameters to pass in to the Clearcut app controller.
    """
    #Clearcut does not support alignment
    raise NotImplementedError, """Clearcut does not support alignment."""

def align_two_alignments(aln1, aln2, params=None):
    """Returns an Alignment object from two existing Alignments.

    aln1, aln2: cogent.core.alignment.Alignment objects, or data that can be
    used to build them.

    params: dict of parameters to pass in to the Clearcut app controller.
    """
    #Clearcut does not support alignment
    raise NotImplementedError, """Clearcut does not support alignment."""

    