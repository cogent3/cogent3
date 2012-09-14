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
from cogent.util.dict2d import Dict2D
from cogent.format.table import phylipMatrix

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
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

    
def build_tree_from_distance_matrix(matrix, best_tree=False, params={},\
    working_dir='/tmp'):
    """Returns a tree from a distance matrix.

    matrix: a square Dict2D object (cogent.util.dict2d)
    
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
    #Turn off input as alignment
    app.Parameters['-a'].off()
    #Input is a distance matrix
    app.Parameters['-d'].on()
    
    if best_tree:
        app.Parameters['-N'].on()
    
    # Turn the dict2d object into the expected input format
    matrix_input, int_keys = _matrix_input_from_dict2d(matrix)

    # Collect result
    result = app(matrix_input)
    
    # Build tree
    tree = DndParser(result['Tree'].read(), constructor=PhyloNode)

    # reassign to original names
    for node in tree.tips():
        node.Name = int_keys[node.Name]

    # Clean up
    result.cleanUp()
    del(app, result, params)

    return tree

def _matrix_input_from_dict2d(matrix):
    """makes input for running clearcut on a matrix from a dict2D object"""
    #clearcut truncates names to 10 char- need to rename before and 
    #reassign after
    
    #make a dict of env_index:full name
    int_keys = dict([('env_' + str(i), k) for i,k in \
            enumerate(sorted(matrix.keys()))])
    #invert the dict
    int_map = {}
    for i in int_keys:
        int_map[int_keys[i]] = i

    #make a new dict2D object with the integer keys mapped to values instead of
    #the original names
    new_dists = []
    for env1 in matrix:
        for env2 in matrix[env1]:
            new_dists.append((int_map[env1], int_map[env2], matrix[env1][env2]))
    int_map_dists = Dict2D(new_dists)
    
    #names will be fed into the phylipTable function - it is the int map names
    names = sorted(int_map_dists.keys())
    rows = []
    #populated rows with values based on the order of names
    #the following code will work for a square matrix only
    for index, key1 in enumerate(names):
        row = []
        for key2 in names:
            row.append(str(int_map_dists[key1][key2]))
        rows.append(row)
    input_matrix = phylipMatrix(rows, names)
    #input needs a trailing whitespace or it will fail!
    input_matrix += '\n'
   
    return input_matrix, int_keys

