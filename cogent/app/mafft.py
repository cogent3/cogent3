#!/usr/bin/env python
"""
Provides an application controller for the commandline version of:
MAFFT v6.602
"""
from cogent.app.parameters import FlagParameter, ValuedParameter, FilePath
from cogent.app.util import CommandLineApplication, ResultPath, \
    get_tmp_filename
from random import choice
from cogent.parse.fasta import MinimalFastaParser
from cogent.core.moltype import DNA, RNA, PROTEIN
from cogent.core.alignment import SequenceCollection, Alignment
from cogent.core.tree import PhyloNode
from cogent.parse.tree import DndParser
from os import remove

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

MOLTYPE_MAP = {'DNA':'--nuc',\
               'RNA':'--nuc',\
               'PROTEIN':'--amino',\
               }

class Mafft(CommandLineApplication):
    """Mafft application controller"""
    
    
    _options ={
    # Algorithm
    
    # Automatically selects an appropriate strategy from L-INS-i, FFT-NS-i
    # and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
    '--auto':FlagParameter(Prefix='--',Name='auto'),\

    # Distance is calculated based on the number of shared 6mers. Default: on
    '--6merpair':FlagParameter(Prefix='--',Name='6merpair'),\

    # All pairwise alignments are computed with the Needleman-Wunsch algorithm.
    # More accurate but slower than --6merpair. Suitable for a set of globally
    # alignable sequences. Applicable to up to ~200 sequences. A combination
    # with --maxiterate 1000 is recommended (G-INS-i). Default: off 
    # (6mer distance is used)
    '--globalpair':FlagParameter(Prefix='--',Name='globalpair'),\

    # All pairwise alignments are computed with the Smith-Waterman algorithm.
    # More accurate but slower than --6merpair. Suitable for a set of locally
    # alignable sequences. Applicable to up to ~200 sequences. A combination
    # with --maxiterate 1000 is recommended (L-INS-i). Default: off
    # (6mer distance is used)
    '--localpair':FlagParameter(Prefix='--',Name='localpair'),\

    # All pairwise alignments are computed with a local algorithm with the
    # generalized affine gap cost (Altschul 1998). More accurate but slower than 
    # --6merpair. Suitable when large internal gaps are expected. Applicable to
    # up to ~200 sequences. A combination with --maxiterate 1000 is recommended
    # (E-INS-i). Default: off (6mer distance is used)
    '--genafpair':FlagParameter(Prefix='--',Name='genafpair'),\

    # All pairwise alignments are computed with FASTA (Pearson and Lipman 1988). 
    # FASTA is required. Default: off (6mer distance is used)
    '--fastapair':FlagParameter(Prefix='--',Name='fastapair'),\

    # Weighting factor for the consistency term calculated from pairwise
    # alignments. Valid when either of --blobalpair, --localpair, --genafpair,
    # --fastapair or --blastpair is selected. Default: 2.7
    '--weighti':ValuedParameter(Prefix='--',Name='weighti',Delimiter=' '),\

    # Guide tree is built number times in the progressive stage. Valid with 6mer 
    # distance. Default: 2
    '--retree':ValuedParameter(Prefix='--',Name='retree',Delimiter=' '),\

    # number cycles of iterative refinement are performed. Default: 0
    '--maxiterate':ValuedParameter(Prefix='--',Name='maxiterate',\
        Delimiter=' '),\
  
    # Use FFT approximation in group-to-group alignment. Default: on
    '--fft':FlagParameter(Prefix='--',Name='fft'),\
    
    # Do not use FFT approximation in group-to-group alignment. Default: off
    '--nofft':FlagParameter(Prefix='--',Name='nofft'),\

    #Alignment score is not checked in the iterative refinement stage. Default:
    # off (score is checked)
    '--noscore':FlagParameter(Prefix='--',Name='noscore'),\

    # Use the Myers-Miller (1988) algorithm. Default: automatically turned on 
    # when the alignment length exceeds 10,000 (aa/nt).
    '--memsave':FlagParameter(Prefix='--',Name='memsave'),\

    # Use a fast tree-building method (PartTree, Katoh and Toh 2007) with the
    # 6mer distance. Recommended for a large number (> ~10,000) of sequences are 
    # input. Default: off
    '--parttree':FlagParameter(Prefix='--',Name='parttree'),\

    # The PartTree algorithm is used with distances based on DP. Slightly more
    # accurate and slower than --parttree. Recommended for a large number
    # (> ~10,000) of sequences are input. Default: off
    '--dpparttree':FlagParameter(Prefix='--',Name='dpparttree'),\

    # The PartTree algorithm is used with distances based on FASTA. Slightly
    # more accurate and slower than --parttree. Recommended for a large number
    # (> ~10,000) of sequences are input. FASTA is required. Default: off
    '--fastaparttree':FlagParameter(Prefix='--',Name='fastaparttree'),\

    # The number of partitions in the PartTree algorithm. Default: 50
    '--partsize':ValuedParameter(Prefix='--',Name='partsize',Delimiter=' '),\

    # Do not make alignment larger than number sequences. Valid only with the
    # --*parttree options. Default: the number of input sequences
    '--groupsize':ValuedParameter(Prefix='--',Name='groupsize',Delimiter=' '),\
 
    # Parameter

    # Gap opening penalty at group-to-group alignment. Default: 1.53
    '--op':ValuedParameter(Prefix='--',Name='op',Delimiter=' '),\

    # Offset value, which works like gap extension penalty, for group-to-group
    # alignment. Deafult: 0.123
    '--ep':ValuedParameter(Prefix='--',Name='ep',Delimiter=' '),\

    # Gap opening penalty at local pairwise alignment. Valid when the
    # --localpair or --genafpair option is selected. Default: -2.00
    '--lop':ValuedParameter(Prefix='--',Name='lop',Delimiter=' '),\

    # Offset value at local pairwise alignment. Valid when the --localpair or 
    # --genafpair option is selected. Default: 0.1
    '--lep':ValuedParameter(Prefix='--',Name='lep',Delimiter=' '),\

    # Gap extension penalty at local pairwise alignment. Valid when the
    # --localpair or --genafpair option is selected. Default: -0.1
    '--lexp':ValuedParameter(Prefix='--',Name='lexp',Delimiter=' '),\

    # Gap opening penalty to skip the alignment. Valid when the --genafpair
    # option is selected. Default: -6.00
    '--LOP':ValuedParameter(Prefix='--',Name='LOP',Delimiter=' '),\

    # Gap extension penalty to skip the alignment. Valid when the --genafpair
    # option is selected. Default: 0.00
    '--LEXP':ValuedParameter(Prefix='--',Name='LEXP',Delimiter=' '),\

    # BLOSUM number matrix (Henikoff and Henikoff 1992) is used. number=30, 45,
    # 62 or 80. Default: 62
    '--bl':ValuedParameter(Prefix='--',Name='bl',Delimiter=' '),\

    # JTT PAM number (Jones et al. 1992) matrix is used. number>0.
    # Default: BLOSUM62
    '--jtt':ValuedParameter(Prefix='--',Name='jtt',Delimiter=' '),\

    # Transmembrane PAM number (Jones et al. 1994) matrix is used. number>0.
    # Default: BLOSUM62
    '--tm':ValuedParameter(Prefix='--',Name='tm',Delimiter=' '),\

    # Use a user-defined AA scoring matrix. The format of matrixfile is the same 
    # to that of BLAST. Ignored when nucleotide sequences are input.
    # Default: BLOSUM62
    '--aamatrix':ValuedParameter(Prefix='--',Name='aamatrix',Delimiter=' '),\

    # Incorporate the AA/nuc composition information into the scoring matrix.
    # Deafult: off
    '--fmodel':FlagParameter(Prefix='--',Name='fmodel'),\
    
    # Output

    # Output format: clustal format. Default: off (fasta format)
    '--clustalout':FlagParameter(Prefix='--',Name='clustalout'),\

    # Output order: same as input. Default: on
    '--inputorder':FlagParameter(Prefix='--',Name='inputorder'),\

    # Output order: aligned. Default: off (inputorder)
    '--reorder':FlagParameter(Prefix='--',Name='reorder'),\

    # Guide tree is output to the input.tree file. Default: off
    '--treeout':FlagParameter(Prefix='--',Name='treeout'),\

    # Do not report progress. Default: off
    '--quiet':FlagParameter(Prefix='--',Name='quiet'),\

# Input

    # Assume the sequences are nucleotide. Deafult: auto
    '--nuc':FlagParameter(Prefix='--',Name='nuc'),\

    # Assume the sequences are amino acid. Deafult: auto
    '--amino':FlagParameter(Prefix='--',Name='amino'),\

    # Seed alignments given in alignment_n (fasta format) are aligned with
    # sequences in input. The alignment within every seed is preserved.
    '--seed':ValuedParameter(Prefix='--',Name='seed',Delimiter=' '),\
    }
    
    _parameters = {}
    _parameters.update(_options)
    _command = "mafft"
    _suppress_stderr=True
    
    def _input_as_seqs(self,data):
        lines = []
        for i,s in enumerate(data):
            #will number the sequences 1,2,3,etc.
            lines.append(''.join(['>',str(i+1)]))
            lines.append(s)
        return self._input_as_lines(lines)
    
    def _tree_out_filename(self):
        if self.Parameters['--treeout'].isOn():
            tree_filename = self._absolute(str(self._input_filename))+'.tree'
        else:
            raise ValueError, "No tree output file specified."
        return tree_filename
    
    def _tempfile_as_multiline_string(self, data):
        """Write a multiline string to a temp file and return the filename.

            data: a multiline string to be written to a file.

           * Note: the result will be the filename as a FilePath object 
            (which is a string subclass).

        """
        filename = FilePath(self.getTmpFilename(self.TmpDir))
        data_file = open(filename,'w')
        data_file.write(data)
        data_file.close()
        return filename
        
    def getHelp(self):
        """Method that points to the Mafft documentation."""
        
        help_str = \
        """
        See Mafft documentation at:
        http://align.bmr.kyushu-u.ac.jp/mafft/software/manual/manual.html
        """
        return help_str
    
    def _get_result_paths(self,data):
        result = {}
        if self.Parameters['--treeout'].isOn():
            out_name = self._tree_out_filename()
            result['Tree'] = ResultPath(Path=out_name,IsWritten=True)
        return result

def align_unaligned_seqs(seqs,moltype,params=None,accurate=False):
    """Aligns unaligned sequences

    seqs: either list of sequence objects or list of strings
    add_seq_names: boolean. if True, sequence names are inserted in the list
        of sequences. if False, it assumes seqs is a list of lines of some
        proper format that the program can handle
    """
    #create SequenceCollection object from seqs
    seq_collection = SequenceCollection(seqs,MolType=moltype)
    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seq_collection.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)
    #Create Mafft app.
    app = Mafft(InputHandler='_input_as_multiline_string',params=params)
    
    #Turn on correct moltype
    moltype_string = moltype.label.upper()
    app.Parameters[MOLTYPE_MAP[moltype_string]].on()
    
    #Do not report progress
    app.Parameters['--quiet'].on()
    
    #More accurate alignment, sacrificing performance.
    if accurate:
        app.Parameters['--globalpair'].on()
        app.Parameters['--maxiterate'].Value=1000
    
    #Get results using int_map as input to app
    res = app(int_map.toFasta())
    #Get alignment as dict out of results
    alignment = dict(MinimalFastaParser(res['StdOut'].readlines()))
    #Make new dict mapping original IDs
    new_alignment = {}
    for k,v in alignment.items():
        new_alignment[int_keys[k]]=v
    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment,MolType=moltype)
    #Clean up
    res.cleanUp()
    del(seq_collection,int_map,int_keys,app,res,alignment)

    return new_alignment


def align_and_build_tree(seqs, moltype, best_tree=False, params={}):
    """Returns an alignment and a tree from Sequences object seqs.
    
    seqs: SequenceCollection object, or data that can be used to build one.
    
    best_tree: if True (default:False), uses a slower but more accurate
    algorithm to build the tree.

    params: dict of parameters to pass in to the Mafft app controller.

    The result will be a tuple containing an Alignment object and a
    cogent.core.tree.PhyloNode object (or None for the alignment and/or tree
    if either fails).
    """
    #Current version of Mafft does not support tree building.
    raise NotImplementedError, """Current version of Mafft does not support tree building."""
    
def build_tree_from_alignment(aln, moltype, best_tree=False, params={},\
    working_dir='/tmp'):
    """Returns a tree from Alignment object aln.

    aln: a cogent.core.alignment.Alignment object, or data that can be used
    to build one.

    best_tree: if True (default:False), uses a slower but more accurate
    algorithm to build the tree.
        NOTE: Mafft does not necessarily support best_tree option.
        Will only return guide tree used to align sequences.  Passing 
        best_tree = True will construct the guide tree 100 times instead
        of default 2 times.
        
        ***Mafft does allow you to get the guide tree back, but the IDs in the
        output guide tree do not match the original IDs in the fasta file
        and are impossible to map.  Sent bug report to Mafft authors; possibly
        expect this option in future version.***

    params: dict of parameters to pass in to the Mafft app controller.

    The result will be an cogent.core.tree.PhyloNode object, or None if tree
    fails.
    """
    #Current version of Mafft does not support tree building.
    raise NotImplementedError, """Current version of Mafft does not support tree building."""
    
def add_seqs_to_alignment(seqs, aln, moltype, params=None, accurate=False):
    """Returns an Alignment object from seqs and existing Alignment.

    seqs: a cogent.core.sequence.Sequence object, or data that can be used
    to build one.

    aln: an cogent.core.alignment.Alignment object, or data that can be used
    to build one

    params: dict of parameters to pass in to the Mafft app controller.
    """
    #create SequenceCollection object from seqs
    seq_collection = SequenceCollection(seqs,MolType=moltype)
    #Create mapping between abbreviated IDs and full IDs
    seq_int_map, seq_int_keys = seq_collection.getIntMap()
    #Create SequenceCollection from int_map.
    seq_int_map = SequenceCollection(seq_int_map,MolType=moltype)
    
    #create Alignment object from aln
    aln = Alignment(aln,MolType=moltype)
    #Create mapping between abbreviated IDs and full IDs
    aln_int_map, aln_int_keys = aln.getIntMap(prefix='seqn_')
    #Create SequenceCollection from int_map.
    aln_int_map = Alignment(aln_int_map,MolType=moltype)
    
    #Update seq_int_keys with aln_int_keys
    seq_int_keys.update(aln_int_keys)
    
    #Create Mafft app.
    app = Mafft(InputHandler='_input_as_multiline_string',\
        params=params,
        SuppressStderr=True)
    
    #Turn on correct moltype
    moltype_string = moltype.label.upper()
    app.Parameters[MOLTYPE_MAP[moltype_string]].on()
    
    #Do not report progress
    app.Parameters['--quiet'].on()
    
    #Add aln_int_map as seed alignment
    app.Parameters['--seed'].on(\
        app._tempfile_as_multiline_string(aln_int_map.toFasta()))
        
    #More accurate alignment, sacrificing performance.
    if accurate:
        app.Parameters['--globalpair'].on()
        app.Parameters['--maxiterate'].Value=1000
    
    #Get results using int_map as input to app
    res = app(seq_int_map.toFasta())
    #Get alignment as dict out of results
    alignment = dict(MinimalFastaParser(res['StdOut'].readlines()))
    
    #Make new dict mapping original IDs
    new_alignment = {}
    for k,v in alignment.items():
        key = k.replace('_seed_','')
        new_alignment[seq_int_keys[key]]=v
    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment,MolType=moltype)
    #Clean up
    res.cleanUp()
    remove(app.Parameters['--seed'].Value)
    del(seq_collection,seq_int_map,seq_int_keys,\
        aln,aln_int_map,aln_int_keys,app,res,alignment)

    return new_alignment

def align_two_alignments(aln1, aln2, moltype, params=None):
    """Returns an Alignment object from two existing Alignments.

    aln1, aln2: cogent.core.alignment.Alignment objects, or data that can be
    used to build them.
        - Mafft profile alignment only works with aligned sequences. Alignment
        object used to handle unaligned sequences.

    params: dict of parameters to pass in to the Mafft app controller.
    """
    #create SequenceCollection object from seqs
    aln1 = Alignment(aln1,MolType=moltype)
    #Create mapping between abbreviated IDs and full IDs
    aln1_int_map, aln1_int_keys = aln1.getIntMap()
    #Create SequenceCollection from int_map.
    aln1_int_map = Alignment(aln1_int_map,MolType=moltype)
    
    #create Alignment object from aln
    aln2 = Alignment(aln2,MolType=moltype)
    #Create mapping between abbreviated IDs and full IDs
    aln2_int_map, aln2_int_keys = aln2.getIntMap(prefix='seqn_')
    #Create SequenceCollection from int_map.
    aln2_int_map = Alignment(aln2_int_map,MolType=moltype)
    
    #Update aln1_int_keys with aln2_int_keys
    aln1_int_keys.update(aln2_int_keys)
    
    #Create Mafft app.
    app = Mafft(InputHandler='_input_as_paths',\
        params=params,
        SuppressStderr=False)
    app._command = 'mafft-profile'
    
    aln1_path = app._tempfile_as_multiline_string(aln1_int_map.toFasta())
    aln2_path = app._tempfile_as_multiline_string(aln2_int_map.toFasta())
    filepaths = [aln1_path,aln2_path]
    
    #Get results using int_map as input to app
    res = app(filepaths)

    #Get alignment as dict out of results
    alignment = dict(MinimalFastaParser(res['StdOut'].readlines()))
    
    #Make new dict mapping original IDs
    new_alignment = {}
    for k,v in alignment.items():
        key = k.replace('_seed_','')
        new_alignment[aln1_int_keys[key]]=v
    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment,MolType=moltype)
    #Clean up
    res.cleanUp()
    remove(aln1_path)
    remove(aln2_path)
    remove('pre')
    remove('trace')
    del(aln1,aln1_int_map,aln1_int_keys,\
        aln2,aln2_int_map,aln2_int_keys,app,res,alignment)

    return new_alignment

    
