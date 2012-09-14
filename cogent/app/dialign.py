#!/usr/bin/env python
"""
Application controller for dialign2-2
"""

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


from cogent.app.parameters import FlagParameter, ValuedParameter
from cogent.app.util import CommandLineApplication, ResultPath, \
    get_tmp_filename, guess_input_handler
from random import choice
from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode
from cogent.parse.fasta import MinimalFastaParser

class Dialign(CommandLineApplication):
    """Dialign application controller"""
    
    _options ={
        # -afc            Creates additional output file "*.afc" containing data of
        #                 all fragments considered for alignment
        #                 WARNING: this file can be HUGE !
        '-afc':FlagParameter(Prefix='-',Name='afc'),
        # -afc_v          like "-afc" but verbose: fragments are explicitly printed
        #                 WARNING: this file can be EVEN BIGGER !
        '-afc_v':FlagParameter(Prefix='-',Name='afc_v'),
        # -anc            Anchored alignment. Requires a file <seq_file>.anc
        #                 containing anchor points.
        '-anc':FlagParameter(Prefix='-',Name='anc'),
        # -cs             if segments are translated, not only the `Watson strand'
        #                 but also the `Crick strand' is looked at.
        '-cs':FlagParameter(Prefix='-',Name='cs'),
        # -cw             additional output file in CLUSTAL W format.
        '-cw':FlagParameter(Prefix='-',Name='cw'),
        # -ds             `dna alignment speed up' - non-translated nucleic acid
        #                 fragments are taken into account only if they start with
        #                 at least two matches. Speeds up DNA alignment at the expense
        #                 of sensitivity.
        '-ds':FlagParameter(Prefix='-',Name='ds'),
        # -fa             additional output file in FASTA format.
        '-fa':FlagParameter(Prefix='-',Name='fa'),
        # -ff             Creates file *.frg containing information about all
        #                 fragments that are part of the respective optimal pairwise
        #                 alignmnets plus information about consistency in the multiple
        #                 alignment
        '-ff':FlagParameter(Prefix='-',Name='ff'),
        # -fn <out_file>  output files are named <out_file>.<extension> .
        '-fn':ValuedParameter('-',Name='fn',Delimiter=' ', IsPath=True),
        #
        #
        # -fop            Creates file *.fop containing coordinates of all fragments
        #                 that are part of the respective pairwise alignments.
        '-fop':FlagParameter(Prefix='-',Name='fop'),
        # -fsm            Creates file *.fsm containing coordinates of all fragments
        #                 that are part of the final alignment
        '-fsm':FlagParameter(Prefix='-',Name='fsm'),
        # -iw             overlap weights switched off (by default, overlap weights are
        #                 used if up to 35 sequences are aligned). This option
        #                 speeds up the alignment but may lead to reduced alignment
        #                 quality.
        '-iw':FlagParameter(Prefix='-',Name='iw'),
        # -lgs            `long genomic sequences' - combines the following options:
        #                 -ma, -thr 2, -lmax 30, -smin 8, -nta, -ff,
        #                 -fop, -ff, -cs, -ds, -pst
        '-lgs':FlagParameter(Prefix='-',Name='lgs'),
        # -lgs_t          Like "-lgs" but with all segment pairs assessed at the
        #                 peptide level (rather than 'mixed alignments' as with the
        #                 "-lgs" option). Therefore faster than -lgs but not very
        #                 sensitive for non-coding regions.
        '-lgs_t':FlagParameter(Prefix='-',Name='lgs_t'),
        # -lmax <x>       maximum fragment length = x  (default: x = 40 or x = 120
        #                 for `translated' fragments). Shorter x speeds up the program
        #                 but may affect alignment quality.
        '-lmax':ValuedParameter('-',Name='lmax',Delimiter=' '),
        # -lo             (Long Output) Additional file *.log with information abut
        #                 fragments selected for pairwise alignment and about
        #                 consistency in multi-alignment proceedure
        '-lo':FlagParameter(Prefix='-',Name='lo'),
        # -ma             `mixed alignments' consisting of P-fragments and N-fragments
        #                 if nucleic acid sequences are aligned.
        '-ma':FlagParameter(Prefix='-',Name='ma'),
        # -mask           residues not belonging to selected fragments are replaced
        #                 by `*' characters in output alignment (rather than being
        #                 printed in lower-case characters)
        '-mask':FlagParameter(Prefix='-',Name='mask'),
        # -mat            Creates file *mat with substitution counts derived from the
        #                 fragments that have been selected for alignment
        '-mat':FlagParameter(Prefix='-',Name='mat'),
        # -mat_thr <t>    Like "-mat" but only fragments with weight score > t
        #                 are considered
        '-mat_thr':ValuedParameter('-',Name='mat_thr',Delimiter=' '),
        # -max_link       "maximum linkage" clustering used to construct sequence tree
        #                 (instead of UPGMA).
        '-max_link':FlagParameter(Prefix='-',Name='max_link'),
        # -min_link       "minimum linkage" clustering used.
        '-min_link':FlagParameter(Prefix='-',Name='min_link'),
        #
        # -mot            "motif" option.
        '-mot':FlagParameter(Prefix='-',Name='mot'),
        # -msf            separate output file in MSF format.
        '-msf':FlagParameter(Prefix='-',Name='msf'),
        # -n              input sequences are nucleic acid sequences. No translation
        #                 of fragments.
        '-n':FlagParameter(Prefix='-',Name='n'),
        # -nt             input sequences are nucleic acid sequences and `nucleic acid
        #                 segments' are translated to `peptide segments'.
        '-nt':FlagParameter(Prefix='-',Name='nt'),
        # -nta            `no textual alignment' - textual alignment suppressed. This
        #                 option makes sense if other output files are of intrest --
        #                 e.g. the fragment files created with -ff, -fop, -fsm or -lo
        '-nta':FlagParameter(Prefix='-',Name='nta'),
        # -o              fast version, resulting alignments may be slightly different.
        '-o':FlagParameter(Prefix='-',Name='o'),
        #
        # -ow             overlap weights enforced (By default, overlap weights are
        #                 used only if up to 35 sequences are aligned since calculating
        #                 overlap weights is time consuming). Warning: overlap weights
        #                 generally improve alignment quality but the running time
        #                 increases in the order O(n^4) with the number of sequences.
        #                 This is why, by default, overlap weights are used only for
        #                 sequence sets with < 35 sequences.
        '-ow':FlagParameter(Prefix='-',Name='ow'),
        # -pst            "print status". Creates and updates a file *.sta with
        #                 information about the current status of the program run.
        #                 This option is recommended if large data sets are aligned
        #                 since it allows the user to estimate the remaining running
        #                 time.
        '-pst':FlagParameter(Prefix='-',Name='pst'),
        # -smin <x>       minimum similarity value for first residue pair (or codon
        #                 pair) in fragments. Speeds up protein alignment or alignment
        #                 of translated DNA fragments at the expense of sensitivity.
        '-smin':ValuedParameter('-',Name='smin',Delimiter=' '),
        # -stars <x>      maximum number of `*' characters indicating degree of
        #                 local similarity among sequences. By default, no stars
        #                 are used but numbers between 0 and 9, instead.
        '-stars':ValuedParameter('-',Name='stars',Delimiter=' '),
        # -stdo           Results written to standard output.
        '-stdo':FlagParameter(Prefix='-',Name='stdo'),
        # -ta             standard textual alignment printed (overrides suppression
        #                 of textual alignments in special options, e.g. -lgs)
        '-ta':FlagParameter(Prefix='-',Name='ta'),
        # -thr <x>        Threshold T = x.
        '-thr':ValuedParameter('-',Name='thr',Delimiter=' '),
        # -xfr            "exclude fragments" - list of fragments can be specified
        #                 that are NOT considered for pairwise alignment
        '-xfr':FlagParameter(Prefix='-',Name='xfr'),
    
    }
    
    _parameters = {}
    _parameters.update(_options)
    _command = "dialign2-2"
    
    def _input_as_seqs(self,data):
        lines = []
        for i,s in enumerate(data):
            #will number the sequences 1,2,3,etc.
            lines.append(''.join(['>',str(i+1)]))
            lines.append(s)
        return self._input_as_lines(lines)
    
    def _align_out_filename(self):
        
        if self.Parameters['-fn'].isOn():
            aln_filename = self._absolute(str(self.Parameters['-fn'].Value))
        else:
            raise ValueError, "No output file specified."
        return aln_filename
    
    def _get_result_paths(self,data):
        
        result = {}
        if self.Parameters['-fn'].isOn():
            out_name = self._align_out_filename()
            result['Align'] = ResultPath(Path=out_name,IsWritten=True)
        return result
    
    def getHelp(self):
        """Dialign help"""
        
        help_str = """
"""
        return help_str
    
