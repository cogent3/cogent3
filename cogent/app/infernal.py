#!/usr/bin/env python
"""
Provides an application controller for the commandline version of:
Infernal 1.0 and 1.0.2 only.
"""
from cogent.app.parameters import FlagParameter, ValuedParameter, FilePath
from cogent.app.util import CommandLineApplication, ResultPath, get_tmp_filename
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.rfam import MinimalRfamParser, ChangedSequence, \
    ChangedRnaSequence, ChangedDnaSequence
from cogent.parse.infernal import CmsearchParser
from cogent.core.moltype import DNA, RNA
from cogent.core.alignment import SequenceCollection, Alignment, DataError
from cogent.format.stockholm import stockholm_from_alignment
from cogent.struct.rna2d import ViennaStructure, wuss_to_vienna
from os import remove

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

MOLTYPE_MAP = {'DNA':'--dna',\
                DNA:'--dna',\
               'RNA':'--rna',\
                RNA:'--rna',\
               }
               
SEQ_CONSTRUCTOR_MAP = {'DNA':ChangedDnaSequence,\
                        DNA:ChangedDnaSequence,\
                       'RNA':ChangedRnaSequence,\
                        RNA:ChangedRnaSequence,\
                       }

class Cmalign(CommandLineApplication):
    """cmalign application controller."""
    _options = {
 
    # -o <f> Save the alignment in Stockholm format to a file <f>. The default
    #   is to write it to standard output.
    '-o':ValuedParameter(Prefix='-',Name='o',Delimiter=' '),\
    
    # -l Turn on the local alignment algorithm. Default is global.
    '-l':FlagParameter(Prefix='-',Name='l'),\
    
    # -p Annotate the alignment with posterior probabilities calculated using
    #   the Inside and Outside algorithms.
    '-p':FlagParameter(Prefix='-',Name='p'),\
    
    # -q Quiet; suppress the verbose banner, and only print the resulting
    #   alignment to stdout.
    '-q':FlagParameter(Prefix='-',Name='q'),\
    
    # --informat <s> Assert that the input seqfile is in format <s>. Do not run
    #   Babelfish format autodection. Acceptable formats are: FASTA, EMBL,
    #   UNIPROT, GENBANK, and DDBJ. <s> is case-insensitive.
    '--informat':ValuedParameter(Prefix='--',Name='informat',Delimiter=' '),\
    
    # --mpi Run as an MPI parallel program.  (see User's Guide for details).
    '--mpi':FlagParameter(Prefix='--',Name='mpi'),\
    
    # Expert Options 
    
    # --optacc Align sequences using the Durbin/Holmes optimal accuracy 
    #   algorithm. This is default behavior, so this option is probably useless. 
    '--optacc':FlagParameter(Prefix='--',Name='optacc'),\
    
    # --cyk Do not use the Durbin/Holmes optimal accuracy alignment to align the 
    #   sequences, instead use the CYK algorithm which determines the optimally
    #   scoring alignment of the sequence to the model. 
    '--cyk':FlagParameter(Prefix='--',Name='cyk'),\
    
    # --sample Sample an alignment from the posterior distribution of
    #   alignments.
    '--sample':FlagParameter(Prefix='--',Name='sample'),\
    
    # -s <n> Set the random number generator seed to <n>, where <n> is a 
    #   positive integer. This option can only be used in combination with 
    #   --sample. The default is to use time() to generate a different seed for
    #   each run, which means that two different runs of cmalign --sample on the
    #   same alignment will give slightly different results. You can use this
    #   option to generate reproducible results.
    '-s':ValuedParameter(Prefix='-',Name='s',Delimiter=' '),\
    
    # --viterbi Do not use the CM to align the sequences, instead use the HMM
    #   Viterbi algorithm to align with a CM Plan 9 HMM.
    '--viterbi':FlagParameter(Prefix='--',Name='viterbi'),\
    
    # --sub Turn on the sub model construction and alignment procedure.
    '--sub':FlagParameter(Prefix='--',Name='sub'),\
    
    # --small Use the divide and conquer CYK alignment algorithm described in 
    #   SR Eddy, BMC Bioinformatics 3:18, 2002.
    '--small':FlagParameter(Prefix='--',Name='small'),\
    
    # --hbanded This option is turned on by default. Accelerate alignment by
    #   pruning away regions of the CM DP matrix that are deemed negligible by
    #   an HMM.
    '--hbanded':FlagParameter(Prefix='--',Name='hbanded'),\
    
    # --nonbanded Turns off HMM banding.
    '--nonbanded':FlagParameter(Prefix='--',Name='nonbanded'),\
    
    # --tau <x> Set the tail loss probability used during HMM band calculation
    #   to <x>.
    '--tau':ValuedParameter(Prefix='--',Name='tau',Delimiter=' '),\
    
    # --mxsize <x> Set the maximum allowable DP matrix size to <x> megabytes.
    '--mxsize':ValuedParameter(Prefix='--',Name='mxsize',Delimiter=' '),\

    # --rna Output the alignments as RNA sequence alignments. This is true by
    #   default.
    '--rna':FlagParameter(Prefix='--',Name='rna'),\
    
    # --dna Output the alignments as DNA sequence alignments.
    '--dna':FlagParameter(Prefix='--',Name='dna'),\
    
    # --matchonly Only include match columns in the output alignment, do not
    #   include any insertions relative to the consensus model.
    '--matchonly':FlagParameter(Prefix='--',Name='matchonly'),\
    
    # --resonly Only include match columns in the output alignment that have at
    #   least 1 residue (non-gap character) in them.
    '--resonly':FlagParameter(Prefix='--',Name='resonly'),\
    
    # --fins Change the behavior of how insert emissions are placed in the 
    #   alignment.
    '--fins':FlagParameter(Prefix='--',Name='fins'),\
    
    # --onepost Modifies behavior of the -p option. Use only one character
    #   instead of two to annotate the posterior probability of each aligned
    #   residue.
    '--onepost':FlagParameter(Prefix='--',Name='onepost'),\
    
    # --withali <f> Reads an alignment from file <f> and aligns it as a single
    #   object to the CM; e.g. the alignment in <f> is held fixed.
    '--withali':ValuedParameter(Prefix='--',Name='withali',Delimiter=' '),\
    
    # --withpknots Must be used in combination with --withali <f>. Propogate
    #   structural information for any pseudoknots that exist in <f> to the
    #   output alignment.
    '--withpknots':FlagParameter(Prefix='--',Name='withpknots'),\
    
    # --rf Must be used in combination with --withali <f>. Specify that the
    #   alignment in <f> has the same "#=GC RF" annotation as the alignment file
    #   the CM was built from using cmbuild and further that the --rf option was 
    #   supplied to cmbuild when the CM was constructed.
    '--rf':FlagParameter(Prefix='--',Name='rf'),\
    
    # --gapthresh <x> Must be used in combination with --withali <f>. Specify
    #   that the --gapthresh <x> option was supplied to cmbuild when the CM was
    #   constructed from the alignment file <f>.
    '--gapthresh':ValuedParameter(Prefix='--',Name='gapthresh',Delimiter=' '),\
    
    # --tfile <f> Dump tabular sequence tracebacks for each individual sequence
    #   to a file <f>. Primarily useful for debugging.
    '--tfile':ValuedParameter(Prefix='--',Name='tfile',Delimiter=' '),\
    
    
    }
    _parameters = {}
    _parameters.update(_options)
    _command = "cmalign"
    _suppress_stderr=True

    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str
    
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

    def _alignment_out_filename(self):
        
        if self.Parameters['-o'].isOn():
            refined_filename = self._absolute(str(\
                self.Parameters['-o'].Value))
        else:
            raise ValueError, 'No alignment output file specified.'
        return refined_filename

    def _get_result_paths(self,data):
        result = {}
        if self.Parameters['-o'].isOn():
            out_name = self._alignment_out_filename()
            result['Alignment'] = ResultPath(Path=out_name,IsWritten=True)
        
        return result

class Cmbuild(CommandLineApplication):
    """cmbuild application controller."""
    _options = {
    
    # -n <s> Name the covariance model <s>. (Does not work if alifile contains
    #   more than one alignment).
    '-n':ValuedParameter(Prefix='-',Name='n',Delimiter=' '),\
    
    # -A Append the CM to cmfile, if cmfile already exists.
    '-A':FlagParameter(Prefix='-',Name='A'),\
    
    # -F Allow cmfile to be overwritten. Normally, if cmfile already exists,
    #   cmbuild exits with an error unless the -A or -F option is set.
    '-F':FlagParameter(Prefix='-',Name='F'),\
    
    # -v Run in verbose output mode instead of using the default single line
    #   tabular format. This output format is similar to that used by older
    #   versions of Infernal.
    '-v':FlagParameter(Prefix='-',Name='v'),\
    
    # --iins Allow informative insert emissions for the CM. By default, all CM
    #   insert emission scores are set to 0.0 bits.
    '--iins':FlagParameter(Prefix='--',Name='iins'),\
    
    # --Wbeta<x> Set the beta tail loss probability for query-dependent banding
    #   (QDB) to <x> The QDB algorithm is used to determine the maximium length
    #   of a hit to the model. For more information on QDB see (Nawrocki and
    #   Eddy, PLoS Computational Biology 3(3): e56).
    '--Wbeta':ValuedParameter(Prefix='--',Name='Wbeta',Delimiter=' '),\
    
    # Expert Options
    
    # --rsearch <f> Parameterize emission scores a la RSEARCH, using the
    #   RIBOSUM matrix in file <f>. For more information see the RSEARCH 
    #   publication (Klein and Eddy, BMC Bioinformatics 4:44, 2003). Actually,
    #   the emission scores will not exactly With --rsearch enabled, all
    #   alignments in alifile must contain exactly one sequence or the --call
    #   option must also be enabled.
    '--rsearch':ValuedParameter(Prefix='--',Name='rsearch',Delimiter=' '),\
     
    # --binary Save the model in a compact binary format. The default is a more
    #   readable ASCII text format.
    '--binary':FlagParameter(Prefix='--',Name='binary'),\
    
    # --rf Use reference coordinate annotation (#=GC RF line, in Stockholm) to
    #   determine which columns are consensus, and which are inserts.
    '--rf':FlagParameter(Prefix='--',Name='rf'),\
    
    # --gapthresh <x> Set the gap threshold (used for determining which columns
    #   are insertions versus consensus; see --rf above) to <x>. The default is
    #   0.5.
    '--gapthresh':ValuedParameter(Prefix='--',Name='gapthresh',Delimiter=' '),\
    
    # --ignorant Strip all base pair secondary structure information from all
    #   input alignments in alifile before building the CM(s).
    '--ignorant':FlagParameter(Prefix='--',Name='ignorant'),\
    
    # --wgsc Use the Gerstein/Sonnhammer/Chothia (GSC) weighting algorithm.
    #   This is the default unless the number of sequences in the alignment
    #   exceeds a cutoff (see --pbswitch), in which case the default becomes
    #   the faster Henikoff position-based weighting scheme.
    '--wgsc':FlagParameter(Prefix='--',Name='wgsc'),\
    
    # --wblosum Use the BLOSUM filtering algorithm to weight the sequences,
    #   instead of the default GSC weighting.
    '--wblosum':FlagParameter(Prefix='--',Name='wblosum'),\
    
    # --wpb Use the Henikoff position-based weighting scheme. This weighting
    #   scheme is automatically used (overriding --wgsc and --wblosum) if the
    #   number of sequences in the alignment exceeds a cutoff (see --pbswitch).
    '--wpb':FlagParameter(Prefix='--',Name='wpb'),\
    
    # --wnone Turn sequence weighting off; e.g. explicitly set all sequence
    #   weights to 1.0.
    '--wnone':FlagParameter(Prefix='--',Name='wnone'),\
    
    # --wgiven Use sequence weights as given in annotation in the input
    #   alignment file. If no weights were given, assume they are all 1.0.
    #   The default is to determine new sequence weights by the Gerstein/
    #   Sonnhammer/Chothia algorithm, ignoring any annotated weights.
    '--wgiven':FlagParameter(Prefix='--',Name='wgiven'),\
    
    # --pbswitch <n> Set the cutoff for automatically switching the weighting
    #   method to the Henikoff position-based weighting scheme to <n>. If the
    #   number of sequences in the alignment exceeds <n> Henikoff weighting is
    #   used. By default <n> is 5000.
    '--pbswitch':ValuedParameter(Prefix='--',Name='pbswitch',Delimiter=' '),\
    
    # --wid <x> Controls the behavior of the --wblosum weighting option by
    #   setting the percent identity for clustering the alignment to <x>.
    '--wid':ValuedParameter(Prefix='--',Name='wid',Delimiter=' '),\
    
    # --eent Use the entropy weighting strategy to determine the effective
    #   sequence number that gives a target mean match state relative entropy.
    '--wgiven':FlagParameter(Prefix='--',Name='wgiven'),\
    
    # --enone Turn off the entropy weighting strategy. The effective sequence
    #   number is just the number of sequences in the alignment.
    '--wgiven':FlagParameter(Prefix='--',Name='wgiven'),\
    
    # --ere <x> Set the target mean match state entropy as <x>. By default the
    #   target entropy 1.46 bits.
    '--ere':ValuedParameter(Prefix='--',Name='ere',Delimiter=' '),\
    
    # --null <f> Read a null model from <f>. The null model defines the
    #   probability of each RNA nucleotide in background sequence, the default
    #   is to use 0.25 for each nucleotide.
    '--null':ValuedParameter(Prefix='--',Name='null',Delimiter=' '),\
    
    # --prior <f> Read a Dirichlet prior from <f>, replacing the default mixture 
    #   Dirichlet.
    '--prior':ValuedParameter(Prefix='--',Name='prior',Delimiter=' '),\
    
    # --ctarget <n> Cluster each alignment in alifile by percent identity.
    #   find a cutoff percent id threshold that gives exactly <n> clusters and
    #   build a separate CM from each cluster. If <n> is greater than the number 
    #   of sequences in the alignment the program will not complain, and each
    #   sequence in the alignment will be its own cluster. Each CM will have a
    #   positive integer appended to its name indicating the order in which it
    #   was built.
    '--ctarget':ValuedParameter(Prefix='--',Name='ctarget',Delimiter=' '),\
    
    # --cmaxid <x> Cluster each sequence alignment in alifile by percent
    #   identity. Define clusters at the cutoff fractional id similarity of <x>
    #   and build a separate CM from each cluster.
    '--cmaxid':ValuedParameter(Prefix='--',Name='cmaxid',Delimiter=' '),\
    
    # --call Build a separate CM from each sequence in each alignment in
    #   alifile. Naming of CMs takes place as described above for --ctarget.
    '--call':FlagParameter(Prefix='--',Name='call'),\
    
    # --corig After building multiple CMs using --ctarget, --cmindiff or --call
    #   as described above, build a final CM using the complete original
    #   alignment from alifile.
    '--corig':FlagParameter(Prefix='--',Name='corig'),\
    
    # --cdump<f> Dump the multiple alignments of each cluster to <f> in
    #   Stockholm format. This option only works in combination with --ctarget,
    #   --cmindiff or --call.
    '--cdump':ValuedParameter(Prefix='--',Name='cdump',Delimiter=' '),\
    
    # --refine <f> Attempt to refine the alignment before building the CM using
    #   expectation-maximization (EM). The final alignment (the alignment used
    #   to build the CM that gets written to cmfile) is written to <f>.
    '--refine':ValuedParameter(Prefix='--',Name='refine',Delimiter=' '),\
    
    # --gibbs Modifies the behavior of --refine so Gibbs sampling is used
    #   instead of EM.
    '--gibbs':FlagParameter(Prefix='--',Name='gibbs'),\
    
    # -s <n> Set the random seed to <n>, where <n> is a positive integer.
    #   This option can only be used in combination with --gibbs. The default is 
    #   to use time() to generate a different seed for each run, which means
    #   that two different runs of cmbuild --refine <f> --gibbs on the same
    #   alignment will give slightly different results. You can use this option
    #   to generate reproducible results.
    '-s':ValuedParameter(Prefix='-',Name='s',Delimiter=' '),\
    
    # -l With --refine, turn on the local alignment algorithm, which allows the
    #   alignment to span two or more subsequences if necessary (e.g. if the
    #   structures of the query model and target sequence are only partially
    #   shared), allowing certain large insertions and deletions in the
    #   structure to be penalized differently than normal indels. The default is 
    #   to globally align the query model to the target sequences.
    '-l':ValuedParameter(Prefix='-',Name='l',Delimiter=' '),\
    
    # -a With --refine, print the scores of each individual sequence alignment.
    '-a':ValuedParameter(Prefix='-',Name='a',Delimiter=' '),\
    
    # --cyk With --refine, align with the CYK algorithm.
    '--cyk':FlagParameter(Prefix='--',Name='cyk'),\
    
    # --sub With --refine, turn on the sub model construction and alignment
    #   procedure.
    '--sub':FlagParameter(Prefix='--',Name='sub'),\
    
    # --nonbanded With --refine, do not use HMM bands to accelerate alignment.
    #   Use the full CYK algorithm which is guaranteed to give the optimal
    #   alignment. This will slow down the run significantly, especially for
    #   large models.
    '--nonbanded':FlagParameter(Prefix='--',Name='nonbanded'),\
    
    # --tau <x> With --refine, set the tail loss probability used during HMM
    #   band calculation to <f>. This is the amount of probability mass within
    #   the HMM posterior probabilities that is considered negligible. The
    #   default value is 1E-7. In general, higher values will result in greater
    #   acceleration, but increase the chance of missing the optimal alignment
    #   due to the HMM bands.
    '--tau':ValuedParameter(Prefix='--',Name='tau',Delimiter=' '),\
    
    # --fins With --refine, change the behavior of how insert emissions are
    #   placed in the alignment.
    '--fins':FlagParameter(Prefix='--',Name='fins'),\
    
    # --mxsize <x> With --refine, set the maximum allowable matrix size for
    #   alignment to <x> megabytes.
    '--mxsize':ValuedParameter(Prefix='--',Name='mxsize',Delimiter=' '),\
    
    # --rdump<x> With --refine, output the intermediate alignments at each
    #   iteration of the refinement procedure (as described above for --refine )
    #   to file <f>.
    '--rdump':ValuedParameter(Prefix='--',Name='rdump',Delimiter=' '),\
    
    }
    _parameters = {}
    _parameters.update(_options)
    _command = "cmbuild"
    _suppress_stderr=True
    
    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str
    
    def _refine_out_filename(self):
        
        if self.Parameters['--refine'].isOn():
            refined_filename = self._absolute(str(\
                self.Parameters['--refine'].Value))
        else:
            raise ValueError, 'No refine output file specified.'
        return refined_filename
    
    def _cm_out_filename(self):
        
        if self.Parameters['-n'].isOn():
            refined_filename = self._absolute(str(\
                self.Parameters['-n'].Value))
        else:
            raise ValueError, 'No cm output file specified.'
        return refined_filename
    
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
    
    def _get_result_paths(self,data):
        result = {}
        if self.Parameters['--refine'].isOn():
            out_name = self._refine_out_filename()
            result['Refined'] = ResultPath(Path=out_name,IsWritten=True)
        if self.Parameters['-n'].isOn():
            cm_name = self._cm_out_filename()
            result['CmFile'] = ResultPath(Path=cm_name,IsWritten=True)
        
        return result

    
class Cmcalibrate(CommandLineApplication):
    """cmcalibrate application controller."""
    _options = {
    
    # -s <n> Set the random number generator seed to <n>, where <n> is a
    #   positive integer. The default is to use time() to generate a different
    #   seed for each run, which means that two different runs of cmcalibrate on 
    #   the same CM will give slightly different E-value and HMM filter
    #   threshold parameters. You can use this option to generate reproducible
    #   results.
    '-s':ValuedParameter(Prefix='-',Name='s',Delimiter=' '),\
    
    # --forecast <n> Predict the running time of the calibration for cmfile and
    #   provided options and exit, DO NOT perform the calibration.
    '--forecast':ValuedParameter(Prefix='--',Name='forecast',Delimiter=' '),\
    
    # --mpi Run as an MPI parallel program.
    '--mpi':FlagParameter(Prefix='--',Name='mpi'),\
    
    # Expert Options
    
    # --exp-cmL-glc <x> Set the length of random sequence to search for the CM
    #   glocal exponential tail fits to <x> megabases (Mb).
    '--exp-cmL-glc':ValuedParameter(Prefix='--',Name='exp-cmL-glc',\
        Delimiter=' '),\
    
    # --exp-cmL-loc <x> Set the length of random sequence to search for the CM
    #   local exponential tail fits to <x> megabases (Mb).
    '--exp-cmL-loc':ValuedParameter(Prefix='--',Name='exp-cmL-loc',\
        Delimiter=' '),\
    
    # --exp-hmmLn-glc <x> Set the minimum random sequence length to search for
    #   the HMM glocal exponential tail fits to <x> megabases (Mb).
    '--exp-hmmLn-glc':ValuedParameter(Prefix='--',Name='exp-hmmLn-glc',\
        Delimiter=' '),\
    
    # --exp-hmmLn-loc <x> Set the minimum random sequence length to search for
    #   the HMM local exponential tail fits to <x> megabases (Mb).
    '--exp-hmmLn-loc':ValuedParameter(Prefix='--',Name='exp-hmmLn-loc',\
        Delimiter=' '),\
    
    # --exp-hmmLx <x> Set the maximum random sequence length to search when
    #   determining HMM E-values to <x> megabases (Mb).
    '--exp-hmmLx':ValuedParameter(Prefix='--',Name='exp-hmmLx',Delimiter=' '),\
    
    # --exp-fract <x> Set the HMM/CM fraction of dynamic programming
    #   calculations to <x>.
    '--exp-fract':ValuedParameter(Prefix='--',Name='exp-fract',Delimiter=' '),\
    
    # --exp-tailn-cglc <x> During E-value calibration of glocal CM search modes
    #   fit the exponential tail to the high scores in the histogram tail that
    #   includes <x> hits per Mb searched.
    '--exp-tailn-cglc':ValuedParameter(Prefix='--',Name='exp-tailn-cglc',\
        Delimiter=' '),\
    
    # --exp-tailn-cloc <x> During E-value calibration of local CM search modes
    #   fit the exponential tail to the high scores in the histogram tail that
    #   includes <x> hits per Mb searched.
    '--exp-tailn-cloc':ValuedParameter(Prefix='--',Name='exp-tailn-cloc',\
        Delimiter=' '),\
    
    # --exp-tailn-hglc <x> During E-value calibration of glocal HMM search modes
    #   fit the exponential tail to the high scores in the histogram tail that
    #   includes <x> hits per Mb searched.
    '--exp-tailn-hglc':ValuedParameter(Prefix='--',Name='exp-tailn-hglc',\
        Delimiter=' '),\
    
    # --exp-tailn-hloc <x> During E-value calibration of local HMM search modes
    #   fit the exponential tail to the high scores in the histogram tail that
    #   includes <x> hits per Mb searched.
    '--exp-tailn-hloc':ValuedParameter(Prefix='--',Name='exp-tailn-hloc',\
        Delimiter=' '),\
    
    # --exp-tailp <x> Ignore the --exp-tailn prefixed options and fit the <x>
    #   fraction right tail of the histogram to exponential tails, for all
    #   search modes.
    '--exp-tailp':ValuedParameter(Prefix='--',Name='exp-tailp',Delimiter=' '),\
    
    # --exp-tailxn <n> With --exp-tailp enforce that the maximum number of hits
    #   in the tail that is fit is <n>.
    '--exp-tailxn':ValuedParameter(Prefix='--',Name='exp-tailxn',\
        Delimiter=' '),\
    
    # --exp-beta <x> During E-value calibration, by default query-dependent
    #   banding (QDB) is used to accelerate the CM search algorithms with a beta
    #   tail loss probability of 1E-15.
    '--exp-beta':ValuedParameter(Prefix='--',Name='exp-beta',Delimiter=' '),\
    
    # --exp-no-qdb Turn of QDB during E-value calibration. This will slow down
    #   calibration, and is not recommended unless you plan on using --no-qdb in
    #   cmsearch.
    '--exp-no-qdb':FlagParameter(Prefix='--',Name='exp-no-qdb'),\
    
    # --exp-hfile <f> Save the histograms fit for the E-value calibration to
    #   file <f>. The format of this file is two tab delimited columns.
    '--exp-hfile':ValuedParameter(Prefix='--',Name='exp-hfile',Delimiter=' '),\
    
    # --exp-sfile <f> Save a survival plot for the E-value calibration to file
    #   <f>. The format of this file is two tab delimited columns.
    '--exp-sfile':ValuedParameter(Prefix='--',Name='exp-sfile',Delimiter=' '),\
    
    # --exp-qqfile <f> Save a quantile-quantile plot for the E-value calibration
    #   to file <f>. The format of this file is two tab delimited columns.
    '--exp-qqfile':ValuedParameter(Prefix='--',Name='exp-qqfile',\
        Delimiter=' '),\
    
    # --exp-ffile <f> Save statistics on the exponential tail statistics to file
    #   <f>. The file will contain the lambda and mu values for exponential
    #   tails fit to tails of different sizes.
    '--exp-ffile':ValuedParameter(Prefix='--',Name='exp-ffile',Delimiter=' '),\
    
    # --fil-N <n> Set the number of sequences sampled and searched for the HMM
    #   filter threshold calibration to <n>. By default, <n> is 10,000.
    '--fil-N':ValuedParameter(Prefix='--',Name='fil-N',Delimiter=' '),\
    
    # --fil-F <x> Set the fraction of sample sequences the HMM filter must be
    #   able to recognize, and allow to survive, to <x>, where <x> is a positive
    #   real number less than or equal to 1.0. By default, <x> is 0.995.
    '--fil-F':ValuedParameter(Prefix='--',Name='fil-F',Delimiter=' '),\
    
    # --fil-xhmm <x> Set the target number of dynamic programming calculations
    #   for a HMM filtered CM QDB search with beta = 1E-7 to <x> times the
    #   number of calculations required to do an HMM search. By default, <x> is
    #   2.0.
    '--fil-xhmm':ValuedParameter(Prefix='--',Name='fil-xhmm',Delimiter=' '),\
    
    # --fil-tau <x> Set the tail loss probability during HMM band calculation
    #   for HMM filter threshold calibration to <x>.
    '--fil-tau':ValuedParameter(Prefix='--',Name='fil-tau',Delimiter=' '),\
    
    # --fil-gemit During HMM filter calibration, always sample sequences from a
    #   globally configured CM, even when calibrating local modes.
    '--fil-gemit':FlagParameter(Prefix='--',Name='fil-gemit'),\
    
    # --fil-dfile <f> Save statistics on filter threshold calibration, including
    #   HMM and CM scores for all sampled sequences, to file <f>.
    '--fil-dfile':ValuedParameter(Prefix='--',Name='fil-dfile',Delimiter=' '),\
    
    # --mxsize <x> Set the maximum allowable DP matrix size to <x> megabytes.
    '--mxsize':ValuedParameter(Prefix='--',Name='mxsize',Delimiter=' '),\

    }
    
    _parameters = {}
    _parameters.update(_options)
    _command = "cmcalibrate"
    _suppress_stderr=True
    
    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str

class Cmemit(CommandLineApplication):
    """cmemit application controller."""
    _options = {
    
    # -o <f> Save the synthetic sequences to file <f> rather than writing them
    #   to stdout.
    '-o':ValuedParameter(Prefix='-',Name='o',Delimiter=' '),\
    
    # -n <n> Generate <n> sequences. Default is 10.
    '-n':ValuedParameter(Prefix='-',Name='n',Delimiter=' '),\
    
    # -u Write the generated sequences in unaligned format (FASTA). This is the
    # default, so this option is probably useless.
    '-u':FlagParameter(Prefix='-',Name='u'),\
    
    # -a Write the generated sequences in an aligned format (STOCKHOLM) with
    #   consensus structure annotation rather than FASTA.
    '-a':FlagParameter(Prefix='-',Name='a'),\
    
    # -c Predict a single majority-rule consensus sequence instead of sampling
    #   sequences from the CM's probability distribution.
    '-c':FlagParameter(Prefix='-',Name='c'),\
    
    # -l Configure the CMs into local mode before emitting sequences. See the
    #   User's Guide for more information on locally configured CMs.
    '-l':FlagParameter(Prefix='-',Name='l'),\
    
    # -s <n> Set the random seed to <n>, where <n> is a positive integer. The
    #   default is to use time() to generate a different seed for each run,
    #   which means that two different runs of cmemit on the same CM will give
    #   different results. You can use this option to generate reproducible
    #   results.
    '-s':ValuedParameter(Prefix='-',Name='s',Delimiter=' '),\
    
    # --rna Specify that the emitted sequences be output as RNA sequences. This
    #   is true by default.
    '--rna':FlagParameter(Prefix='--',Name='rna'),\
    
    # --dna Specify that the emitted sequences be output as DNA sequences. By
    #   default, the output alphabet is RNA.
    '--dna':FlagParameter(Prefix='--',Name='dna'),\
    
    # --tfile <f> Dump tabular sequence parsetrees (tracebacks) for each emitted
    #   sequence to file <f>. Primarily useful for debugging.
    '--tfile':ValuedParameter(Prefix='--',Name='tfile',Delimiter=' '),\
    
    # --exp <x> Exponentiate the emission and transition probabilities of the CM
    #   by <x> and then renormalize those distributions before emitting
    #   sequences.
    '--exp':ValuedParameter(Prefix='--',Name='exp',Delimiter=' '),\
    
    # --begin <n> Truncate the resulting alignment by removing all residues
    #   before consensus column <n>, where <n> is a positive integer no greater
    #   than the consensus length of the CM. Must be used in combination with
    #   --end and either -a or --shmm (a developer option).
    '--begin':ValuedParameter(Prefix='--',Name='begin',Delimiter=' '),\
    
    # --end <n> Truncate the resulting alignment by removing all residues after
    #   consensus column <n>, where <n> is a positive integer no greater than
    #   the consensus length of the CM. Must be used in combination with --begin
    #   and either -a or --shmm (a developer option).
    '--end':ValuedParameter(Prefix='--',Name='end',Delimiter=' '),\
    
    }
    _parameters = {}
    _parameters.update(_options)
    _command = "cmemit"
    _suppress_stderr=True
    
    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str
    
class Cmscore(CommandLineApplication):
    """cmscore application controller."""
    _options = {
    
    # -n <n> Set the number of sequences to generate and align to <n>. This
    #   option is incompatible with the --infile option.
    '-n':ValuedParameter(Prefix='-',Name='n',Delimiter=' '),\
    
    # -l Turn on the local alignment algorithm, which allows the alignment to
    #   span two or more subsequences if necessary (e.g. if the structures of
    #   the query model and target sequence are only partially shared), allowing
    #   certain large insertions and deletions in the structure to be penalized
    #   differently than normal indels. The default is to globally align the
    #   query model to the target sequences.
    '-l':FlagParameter(Prefix='-',Name='l'),\
    
    # -s <n> Set the random seed to <n>, where <n> is a positive integer. The
    #   default is to use time() to generate a different seed for each run,
    #   which means that two different runs of cmscore on the same CM will give
    #   different results. You can use this option to generate reproducible
    #   results. The random number generator is used to generate sequences to
    #   score, so -s is incompatible with the --infile option which supplies
    #   the sequences to score in an input file.
    '-s':ValuedParameter(Prefix='-',Name='s',Delimiter=' '),\
    
    # -a Print individual timings and score comparisons for each sequence in
    #   seqfile. By default only summary statistics are printed.
    '-a':FlagParameter(Prefix='-',Name='a'),\
    
    # --sub Turn on the sub model construction and alignment procedure.
    '--sub':FlagParameter(Prefix='--',Name='sub'),\
    
    # --mxsize <x> Set the maximum allowable DP matrix size to <x> megabytes.
    '--mxsize':ValuedParameter(Prefix='--',Name='mxsize',Delimiter=' '),\
    
    # --mpi Run as an MPI parallel program.
    '--mpi':FlagParameter(Prefix='--',Name='mpi'),\
    
    # Expert Options
    
    # --emit Generate sequences to score by sampling from the CM.
    '--emit':FlagParameter(Prefix='--',Name='emit'),\
    
    # --random Generate sequences to score by sampling from the CMs null
    #   distribution. This option turns the --emit option off.
    '--random':FlagParameter(Prefix='--',Name='random'),\
    
    # --infile <f> Sequences to score are read from the file <f>. All the
    #   sequences from <f> are read and scored, the -n and -s options are
    #   incompatible with --infile.
    '--infile':ValuedParameter(Prefix='--',Name='infile',Delimiter=' '),\
    
    # --outfile <f> Save generated sequences that are scored to the file <f> in
    #   FASTA format. This option is incompatible with the --infile option.
    '--outfile':ValuedParameter(Prefix='--',Name='outfile',Delimiter=' '),\
    
    # --Lmin <n1> Must be used in combination with --random and --Lmax <n2>.
    '--Lmin':ValuedParameter(Prefix='--',Name='Lmin',Delimiter=' '),\
    
    # --pad Must be used in combination with --emit and --search. Add <n> cm->W
    #   (max hit length) minus L (sequence <x> length) residues to the 5' and 3'
    #   end of each emitted sequence <x>.
    '--pad':FlagParameter(Prefix='--',Name='pad'),\
    
    # --hbanded Specify that the second stage alignment algorithm be HMM banded
    #   CYK. This option is on by default.
    '--hbanded':FlagParameter(Prefix='--',Name='hbanded'),\
    
    # --tau <x> For stage 2 alignment, set the tail loss probability used during
    #   HMM band calculation to <x>.
    '--tau':ValuedParameter(Prefix='--',Name='tau',Delimiter=' '),\
    
    # --aln2bands With --search, when calculating HMM bands, use an HMM
    #   alignment algorithm instead of an HMM search algorithm.
    '--aln2bands':FlagParameter(Prefix='--',Name='aln2bands'),\
    
    # --hsafe For stage 2 HMM banded alignment, realign any sequences with a
    #   negative alignment score using non-banded CYK to guarantee finding the
    #   optimal alignment.
    '--hsafe':FlagParameter(Prefix='--',Name='hsafe'),\
    
    # --nonbanded Specify that the second stage alignment algorithm be standard,
    #   non-banded, non-D&C CYK. When --nonbanded is enabled, the program fails
    #   with a non-zero exit code and prints an error message if the parsetree
    #   score for any sequence from stage 1 D&C alignment and stage 2 alignment
    #   differs by more than 0.01 bits. In theory, this should never happen as
    #   both algorithms are guaranteed to determine the optimal parsetree. For
    #   larger RNAs (more than 300 residues) if memory is limiting, --nonbanded
    #   should be used in combination with --scoreonly.
    '--nonbanded':FlagParameter(Prefix='--',Name='nonbanded'),\
    
    # --scoreonly With --nonbanded during the second stage standard non-banded
    #   CYK alignment, use the "score only" variant of the algorithm to save
    #   memory, and don't recover a parse tree.
    '--scoreonly':FlagParameter(Prefix='--',Name='scoreonly'),\
    
    # --viterbi Specify that the second stage alignment algorithm be Viterbi to
    #   a CM Plan 9 HMM.
    '--viterbi':FlagParameter(Prefix='--',Name='viterbi'),\
    
    # --search Run all algorithms in scanning mode, not alignment mode.
    '--search':FlagParameter(Prefix='--',Name='search'),\
    
    # --inside With --search Compare the non-banded scanning Inside algorithm to
    #   the HMM banded scanning Inside algorith, instead of using CYK versions.
    '--inside':FlagParameter(Prefix='--',Name='inside'),\
    
    # --forward With --search Compare the scanning Forward scoring algorithm
    #   against CYK.
    '--forward':FlagParameter(Prefix='--',Name='forward'),\
    
    # --taus <n> Specify the first alignment algorithm as non-banded D&C CYK,
    #   and multiple stages of HMM banded CYK alignment. The first HMM banded
    #   alignment will use tau=1E-<x>, which will be the highest value of tau
    #   used. Must be used in combination with --taue.
    '--taus':ValuedParameter(Prefix='--',Name='taus',Delimiter=' '),\
    
    # --taue <n> Specify the first alignment algorithm as non-banded D&C CYK,
    #   and multiple stages of HMM banded CYK alignment. The final HMM banded
    #   alignment will use tau=1E-<x>, which will be the lowest value of tau
    #   used. Must be used in combination with --taus.
    '--taue':ValuedParameter(Prefix='--',Name='taue',Delimiter=' '),\
    
    # --tfile <f> Print the parsetrees for each alignment of each sequence to
    #   file <f>.
    '--tfile':ValuedParameter(Prefix='--',Name='tfile',Delimiter=' '),\

    }
    _parameters = {}
    _parameters.update(_options)
    _command = "cmscore"
    _suppress_stderr=True
    
    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str

class Cmsearch(CommandLineApplication):
    """cmsearch application controller."""
    _options = {
    
    # -o <f> Save the high-scoring alignments of hits to a file <f>. The default
    #   is to write them to standard output.
    '-o':ValuedParameter(Prefix='-',Name='o',Delimiter=' '),\
    
    # -g <f> Turn on the 'glocal' alignment algorithm, local with respect to the
    #   target database, and global with respect to the model. By default, the
    #   local alignment algorithm is used which is local with respect to both
    #   the target sequence and the model.
    '-g':ValuedParameter(Prefix='-',Name='g',Delimiter=' '),\
    
    # -p Append posterior probabilities to alignments of hits.
    '-p':FlagParameter(Prefix='-',Name='p'),\
    
    # -x Annotate non-compensatory basepairs and basepairs that include a gap in
    #   the left and/or right half of the pair with x's in the alignments of
    #   hits.
    '-x':FlagParameter(Prefix='-',Name='x'),\
    
    # -Z <x> Calculate E-values as if the target database size was <x> megabases
    #   (Mb). Ignore the actual size of the database. This option is only valid
    #   if the CM file has been calibrated. Warning: the predictions for timings
    #   and survival fractions will be calculated as if the database was of size
    #   <x> Mb, which means they will be inaccurate.
    '-Z':ValuedParameter(Prefix='-',Name='Z',Delimiter=' '),\
    
    # --toponly Only search the top (Watson) strand of the sequences in seqfile.
    #   By default, both strands are searched.
    '--toponly':FlagParameter(Prefix='--',Name='toponly'),\
    
    # --bottomonly Only search the bottom (Crick) strand of the sequences in
    #   seqfile. By default, both strands are searched.
    '--bottomonly':FlagParameter(Prefix='--',Name='bottomonly'),\
    
    # --forecast <n> Predict the running time of the search with provided files
    #   and options and exit, DO NOT perform the search. This option is only
    #   available with calibrated CM files.
    '--forecast':ValuedParameter(Prefix='--',Name='forecast',Delimiter=' '),\
    
    # --informat <s> Assert that the input seqfile is in format <s>. Do not run
    #   Babelfish format autodection. This increases the reliability of the
    #   program somewhat, because the Babelfish can make mistakes; particularly
    #   recommended for unattended, high-throughput runs of @PACKAGE@. <s> is
    #   case-insensitive. Acceptable formats are: FASTA, EMBL, UNIPROT, GENBANK,
    #   and DDBJ. <s> is case-insensitive.
    '--informat':ValuedParameter(Prefix='--',Name='informat',Delimiter=' '),\
    
    # --mxsize <x> Set the maximum allowable DP matrix size to <x> megabytes.
    '--mxsize':ValuedParameter(Prefix='--',Name='mxsize',Delimiter=' '),\
    
    # --mpi Run as an MPI parallel program.
    '--mpi':FlagParameter(Prefix='--',Name='mpi'),\
    
    # Expert Options
    
    # --inside Use the Inside algorithm for the final round of searching. This
    #   is true by default.
    '--inside':FlagParameter(Prefix='--',Name='inside'),\
    
    # --cyk Use the CYK algorithm for the final round of searching.
    '--cyk':FlagParameter(Prefix='--',Name='cyk'),\
    
    # --viterbi Search only with an HMM. This is much faster but less sensitive
    #   than a CM search. Use the Viterbi algorithm for the HMM search.
    '--viterbi':FlagParameter(Prefix='--',Name='viterbi'),\
    
    # --forward Search only with an HMM. This is much faster but less sensitive
    #   than a CM search. Use the Forward algorithm for the HMM search.
    '--forward':FlagParameter(Prefix='--',Name='forward'),\
    
    # -E <x> Set the E-value cutoff for the per-sequence/strand ranked hit list
    #   to <x>, where <x> is a positive real number.
    '-E':ValuedParameter(Prefix='-',Name='E',Delimiter=' '),\
    
    # -T <x> Set the bit score cutoff for the per-sequence ranked hit list to
    #   <x>, where <x> is a positive real number.
    '-T':ValuedParameter(Prefix='-',Name='T',Delimiter=' '),\
    
    # --nc Set the bit score cutoff as the NC cutoff value used by Rfam curators
    #   as the noise cutoff score.
    '--nc':FlagParameter(Prefix='--',Name='nc'),\
    
    # --ga Set the bit score cutoff as the GA cutoff value used by Rfam curators
    #   as the gathering threshold.
    '--ga':FlagParameter(Prefix='--',Name='ga'),\
    
    # --tc Set the bit score cutoff as the TC cutoff value used by Rfam curators
    #   as the trusted cutoff.
    '--tc':FlagParameter(Prefix='--',Name='tc'),\
    
    # --no-qdb Do not use query-dependent banding (QDB) for the final round of
    #   search.
    '--no-qdb':FlagParameter(Prefix='--',Name='no-qdb'),\
    
    # --beta " <x>" For query-dependent banding (QDB) during the final round of
    #   search, set the beta parameter to <x> where <x> is any positive real
    #   number less than 1.0.
    '--beta':ValuedParameter(Prefix='--',Name='beta',Delimiter=' '),\
    
    # --hbanded Use HMM bands to accelerate the final round of search.
    #   Constraints for the CM search are derived from posterior probabilities
    #   from an HMM. This is an experimental option and it is not recommended
    #   for use unless you know exactly what you're doing.
    '--hbanded':FlagParameter(Prefix='--',Name='hbanded'),\
    
    # --tau <x> Set the tail loss probability during HMM band calculation to
    #   <x>.
    '--tau':ValuedParameter(Prefix='--',Name='tau',Delimiter=' '),\
    
    # --fil-no-hmm Turn the HMM filter off.
    '--fil-no-hmm':FlagParameter(Prefix='--',Name='fil-no-hmm'),\
    
    # --fil-no-qdb Turn the QDB filter off.
    '--fil-no-qdb':FlagParameter(Prefix='--',Name='fil-no-qdb'),\
    
    # --fil-beta For the QDB filter, set the beta parameter to <x> where <x> is
    #   any positive real number less than 1.0.
    '--fil-beta':FlagParameter(Prefix='--',Name='fil-beta'),\
    
    # --fil-T-qdb <x> Set the bit score cutoff for the QDB filter round to <x>,
    #   where <x> is a positive real number.
    '--fil-T-qdb':ValuedParameter(Prefix='--',Name='fil-T-qdb',Delimiter=' '),\
    
    # --fil-T-hmm <x> Set the bit score cutoff for the HMM filter round to <x>,
    #   where <x> is a positive real number.
    '--fil-T-hmm':ValuedParameter(Prefix='--',Name='fil-T-hmm',Delimiter=' '),\
    
    # --fil-E-qdb <x> Set the E-value cutoff for the QDB filter round. <x>,
    #   where <x> is a positive real number. Hits with E-values better than
    #   (less than) or equal to this threshold will survive and be passed to the
    #   final round. This option is only available if the CM file has been
    #   calibrated.
    '--fil-E-qdb':ValuedParameter(Prefix='--',Name='fil-E-qdb',Delimiter=' '),\
    
    # --fil-E-hmm <x> Set the E-value cutoff for the HMM filter round. <x>,
    #   where <x> is a positive real number. Hits with E-values better than
    #   (less than) or equal to this threshold will survive and be passed to the
    #   next round, either a QDB filter round, or if the QDB filter is disable,
    #   to the final round of search. This option is only available if the CM
    #   file has been calibrated.
    '--fil-E-hmm':ValuedParameter(Prefix='--',Name='fil-E-hmm',Delimiter=' '),\
    
    # --fil-Smax-hmm <x> Set the maximum predicted survival fraction for an HMM
    #   filter as <x>, where <x> is a positive real number less than 1.0.
    '--fil-Smax-hmm':ValuedParameter(Prefix='--',Name='fil-Smax-hmm',\
        Delimiter=' '),\
    
    # --noalign Do not calculate and print alignments of each hit, only print
    #   locations and scores.
    '--noalign':FlagParameter(Prefix='--',Name='noalign'),\
    
    # --aln-hbanded Use HMM bands to accelerate alignment during the hit
    #   alignment stage.
    '--aln-hbanded':FlagParameter(Prefix='--',Name='aln-hbanded'),\
    
    # --aln-optacc Calculate alignments of hits from final round of search using
    #   the optimal accuracy algorithm which computes the alignment that
    #   maximizes the summed posterior probability of all aligned residues given
    #   the model, which can be different from the highest scoring one.
    '--aln-optacc':FlagParameter(Prefix='--',Name='aln-optacc'),\
    
    # --tabfile <f> Create a new output file <f> and print tabular results to
    #   it.
    '--tabfile':ValuedParameter(Prefix='--',Name='tabfile',Delimiter=' '),\
    
    # --gcfile <f> Create a new output file <f> and print statistics of the GC
    #   content of the sequences in seqfile to it.
    '--gcfile':ValuedParameter(Prefix='--',Name='gcfile',Delimiter=' '),\
    
    # --rna Output the hit alignments as RNA sequences alignments. This is true
    #   by default.
    '--rna':FlagParameter(Prefix='--',Name='rna'),\
    
    # --dna Output the hit alignments as DNA sequence alignments.
    '--dna':FlagParameter(Prefix='--',Name='dna'),\

    }
    _parameters = {}
    _parameters.update(_options)
    _command = "cmsearch"
    _suppress_stderr=True
    
    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str
    
    def _tabfile_out_filename(self):
        
        if self.Parameters['--tabfile'].isOn():
            tabfile_filename = self._absolute(str(\
                self.Parameters['--tabfile'].Value))
        else:
            raise ValueError, 'No tabfile output file specified.'
        return tabfile_filename
    
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
    
    def _get_result_paths(self,data):
        result = {}
        if self.Parameters['--tabfile'].isOn():
            out_name = self._tabfile_out_filename()
            result['SearchResults'] = ResultPath(Path=out_name,IsWritten=True)
        
        return result

class Cmstat(CommandLineApplication):
    """cmstat application controller."""
    _options = {
    
    # -g Turn on the 'glocal' alignment algorithm, local with respect to the
    #   target database, and global with respect to the model. By default, the
    #   model is configured for local alignment which is local with respect to
    #   both the target sequence and the model.
    '-g':FlagParameter(Prefix='-',Name='g'),\
    
    # -m print general statistics on the models in cmfile and the alignment it
    #   was built from.
    '-m':FlagParameter(Prefix='-',Name='m'),\
    
    # -Z <x> Calculate E-values as if the target database size was <x> megabases
    #   (Mb). Ignore the actual size of the database. This option is only valid
    #   if the CM file has been calibrated.
    '-Z':ValuedParameter(Prefix='-',Name='Z',Delimiter=' '),\
    
    # --all print all available statistics
    '--all':FlagParameter(Prefix='--',Name='all'),\
    
    # --le print local E-value statistics. This option only works if cmfile has
    #   been calibrated with cmcalibrate.
    '--le':FlagParameter(Prefix='--',Name='le'),\
    
    # --ge print glocal E-value statistics. This option only works if cmfile has
    #   been calibrated with cmcalibrate.
    '--ge':FlagParameter(Prefix='--',Name='ge'),\
    
    # --beta <x> With the --search option set the beta parameter for the query-
    #   dependent banding algorithm stages to <x> Beta is the probability mass
    #   considered negligible during band calculation. The default is 1E-7.
    '--beta':ValuedParameter(Prefix='--',Name='beta',Delimiter=' '),\
    
    # --qdbfile <f> Save the query-dependent bands (QDBs) for each state to file
    #   <f>
    '--qdbfile':ValuedParameter(Prefix='--',Name='qdbfile',Delimiter=' '),\
    
    # Expert Options
    
    # --lfi Print the HMM filter thresholds for the range of relevant CM bit
    #   score cutoffs for searches with locally configured models using the
    #   Inside algorithm.
    '--lfi':FlagParameter(Prefix='--',Name='lfi'),\
    
    # --gfi Print the HMM filter thresholds for the range of relevant CM bit
    #   score cutoffs for searches with globally configured models using the
    #   Inside algorithm.
    '--gfi':FlagParameter(Prefix='--',Name='gfi'),\
    
    # --lfc Print the HMM filter thresholds for the range of relevant CM bit
    #   score cutoffs for searches with locally configured models using the CYK
    #   algorithm.
    '--lfc':FlagParameter(Prefix='--',Name='lfc'),\
    
    # --gfc Print the HMM filter thresholds for the range of relevant CM bit
    #   score cutoffs for searches with globally configured models using the CYK
    #   algorithm.
    '--gfc':FlagParameter(Prefix='--',Name='gfc'),\
    
    # -E <x> Print filter threshold statistics for an HMM filter if a final CM
    #   E-value cutoff of <x> were to be used for a run of cmsearch on 1 MB of
    #   sequence.
    '-E':ValuedParameter(Prefix='-',Name='E',Delimiter=' '),\
    
    # -T <x> Print filter threshold statistics for an HMM filter if a final CM
    #   bit score cutoff of <x> were to be used for a run of cmsearch.
    '-T':ValuedParameter(Prefix='-',Name='T',Delimiter=' '),\
    
    # --nc Print filter threshold statistics for an HMM filter if a CM bit score
    #   cutoff equal to the Rfam NC cutoff were to be used for a run of
    #   cmsearch.
    '--nc':FlagParameter(Prefix='--',Name='nc'),\
    
    # --ga Print filter threshold statistics for an HMM filter if a CM bit score
    #   cutoff of Rfam GA cutoff value were to be used for a run of cmsearch.
    '--ga':FlagParameter(Prefix='--',Name='ga'),\
    
    # --tc Print filter threshold statistics for an HMM filter if a CM bit score
    #   cutoff equal to the Rfam TC cutoff value were to be used for a run of
    #   cmsearch.
    '--tc':FlagParameter(Prefix='--',Name='tc'),\
    
    # --seqfile <x> With the -E option, use the database size of the database in
    #   <x> instead of the default database size of 1 MB.
    '--seqfile':ValuedParameter(Prefix='--',Name='seqfile',Delimiter=' '),\
    
    # --toponly In combination with --seqfile <x> option, only consider the top
    #   strand of the database in <x> instead of both strands. --search perform
    #   an experiment to determine how fast the CM(s) can search with different
    #   search algorithms.
    '--toponly':FlagParameter(Prefix='--',Name='toponly'),\
    
    # --cmL <n> With the --search option set the length of sequence to search
    #   with CM algorithms as <n> residues. By default, <n> is 1000.
    '--cmL':ValuedParameter(Prefix='--',Name='cmL',Delimiter=' '),\
    
    # --hmmL <n> With the --search option set the length of sequence to search
    #   with HMM algorithms as <n> residues. By default, <n> is 100,000.
    '--hmmL':ValuedParameter(Prefix='--',Name='hmmL',Delimiter=' '),\
    
    # --efile <f> Save a plot of cmsearch HMM filter E value cutoffs versus CM
    #   E-value cutoffs in xmgrace format to file <f>. This option must be used
    #   in combination with --lfi, --gfi, --lfc or --gfc.
    '--efile':ValuedParameter(Prefix='--',Name='efile',Delimiter=' '),\
    
    # --bfile <f> Save a plot of cmsearch HMM bit score cutoffs versus CM bit
    #   score cutoffs in xmgrace format to file <f>. This option must be used in
    #   combination with --lfi, --gfi, --lfc or --gfc.
    '--bfile':ValuedParameter(Prefix='--',Name='bfile',Delimiter=' '),\
    
    # --sfile <f> Save a plot of cmsearch predicted survival fraction from the
    #   HMM filter versus CM E value cutoff in xmgrace format to file <f>. This
    #   option must be used in combination with --lfi, --gfi, --lfc or --gfc.
    '--sfile':ValuedParameter(Prefix='--',Name='sfile',Delimiter=' '),\
    
    # --xfile <f> Save a plot of 'xhmm' versus CM E value cutoff in xmgrace
    #   format to file <f> 'xhmm' is the ratio of the number of dynamic
    #   programming calculations predicted to be required for the HMM filter and
    #   the CM search of the filter survivors versus the number of dynamic
    #   programming calculations for the filter alone. This option must be
    #   used in combination with --lfi, --gfi, --lfc or --gfc.
    '--xfile':ValuedParameter(Prefix='--',Name='xfile',Delimiter=' '),\
    
    # --afile <f> Save a plot of the predicted acceleration for an HMM filtered
    #   search versus CM E value cutoff in xmgrace format to file <f>. This
    #   option must be used in combination with --lfi, --gfi, --lfc or --gfc.
    '--afile':ValuedParameter(Prefix='--',Name='afile',Delimiter=' '),\
    
    # --bits With --efile, --sfile, --xfile, and --afile use CM bit score
    #   cutoffs instead of CM E value cutoffs for the x-axis values of the plot.
    '--bits':FlagParameter(Prefix='--',Name='bits'),\

    }
    _parameters = {}
    _parameters.update(_options)
    _command = "cmstat"
    _suppress_stderr=True
    
    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str

def cmbuild_from_alignment(aln, structure_string, refine=False, \
    return_alignment=False,params=None):
    """Uses cmbuild to build a CM file given an alignment and structure string.
    
        - aln: an Alignment object or something that can be used to construct
            one.  All sequences must be the same length.
        - structure_string: vienna structure string representing the consensus
            stucture for the sequences in aln.  Must be the same length as the
            alignment.
        - refine: refine the alignment and realign before building the cm.
            (Default=False)
        - return_alignment: Return (in Stockholm format) alignment file used to
            construct the CM file.  This will either be the original alignment
            and structure string passed in, or the refined alignment if --refine 
            was used. (Default=False)
            - Note.  This will be a string that can either be written to a file
                or parsed.
    """
    aln = Alignment(aln)
    if len(structure_string) != aln.SeqLen:
        raise ValueError, """Structure string is not same length as alignment.  Structure string is %s long. Alignment is %s long."""%(len(structure_string),\
        aln.SeqLen)
    else:
        struct_dict = {'SS_cons':structure_string}
    #Make new Cmbuild app instance.
    app = Cmbuild(InputHandler='_input_as_paths',WorkingDir='/tmp',\
        params=params)
    
    #turn on refine flag if True.
    if refine:
        app.Parameters['--refine'].on(get_tmp_filename(app.WorkingDir))
        
    #Get alignment in Stockholm format
    aln_file_string = stockholm_from_alignment(aln,GC_annotation=struct_dict)
    
    #get path to alignment filename
    aln_path = app._input_as_multiline_string(aln_file_string)
    cm_path = aln_path.split('.txt')[0]+'.cm'
    app.Parameters['-n'].on(cm_path)
    
    filepaths = [cm_path,aln_path]
    
    res = app(filepaths)
    
    cm_file = res['CmFile'].read()
    
    if return_alignment:
        #If alignment was refined, return refined alignment and structure,
        # otherwise return original alignment and structure.
        if refine:
            aln_file_string = res['Refined'].read()
        res.cleanUp()
        return cm_file, aln_file_string
    #Just return cm_file
    else:
        res.cleanUp()
        return cm_file


def cmbuild_from_file(stockholm_file_path, refine=False,return_alignment=False,\
    params=None):
    """Uses cmbuild to build a CM file given a stockholm file.
    
        - stockholm_file_path: a path to a stockholm file.  This file should
            contain a multiple sequence alignment formated in Stockholm format. 
            This must contain a sequence structure line:
                #=GC SS_cons <structure string>
        - refine: refine the alignment and realign before building the cm.
            (Default=False)
        - return_alignment: Return alignment and structure string used to
            construct the CM file.  This will either be the original alignment
            and structure string passed in, or the refined alignment if
            --refine was used. (Default=False)
    """
    #get alignment and structure string from stockholm file.
    info, aln, structure_string = \
        list(MinimalRfamParser(open(stockholm_file_path,'U'),\
            seq_constructor=ChangedSequence))[0]
    
    #call cmbuild_from_alignment.
    res = cmbuild_from_alignment(aln, structure_string, refine=refine, \
        return_alignment=return_alignment,params=params)
    return res

def cmalign_from_alignment(aln, structure_string, seqs, moltype,\
    include_aln=True,refine=False, return_stdout=False,params=None,\
    cmbuild_params=None):
    """Uses cmbuild to build a CM file, then cmalign to build an alignment.
    
        - aln: an Alignment object or something that can be used to construct
            one.  All sequences must be the same length.
        - structure_string: vienna structure string representing the consensus
            stucture for the sequences in aln.  Must be the same length as the
            alignment.
        - seqs: SequenceCollection object or something that can be used to
            construct one, containing unaligned sequences that are to be aligned 
            to the aligned sequences in aln.
        - moltype: Cogent moltype object.  Must be RNA or DNA.
        - include_aln: Boolean to include sequences in aln in final alignment.
            (Default=True)
        - refine: refine the alignment and realign before building the cm.
            (Default=False)
        - return_stdout: Boolean to return standard output from infernal.  This
            includes alignment and structure bit scores and average
            probabilities for each sequence. (Default=False)
    """
    #NOTE: Must degap seqs or Infernal well seg fault!
    seqs = SequenceCollection(seqs,MolType=moltype).degap()
    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seqs.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)
    
    cm_file, aln_file_string = cmbuild_from_alignment(aln, structure_string,\
        refine=refine,return_alignment=True,params=cmbuild_params)
    
    if params is None:
        params = {}    
    params.update({MOLTYPE_MAP[moltype]:True})
    
    app = Cmalign(InputHandler='_input_as_paths',WorkingDir='/tmp',\
        params=params)
    app.Parameters['--informat'].on('FASTA')
    
    #files to remove that aren't cleaned up by ResultPath object
    to_remove = []    
    #turn on --withali flag if True.
    if include_aln:
        app.Parameters['--withali'].on(\
            app._tempfile_as_multiline_string(aln_file_string))
        #remove this file at end
        to_remove.append(app.Parameters['--withali'].Value)
    
    seqs_path = app._input_as_multiline_string(int_map.toFasta())
    cm_path = app._tempfile_as_multiline_string(cm_file)
    
    #add cm_path to to_remove
    to_remove.append(cm_path)
    paths = [cm_path,seqs_path]

    app.Parameters['-o'].on(get_tmp_filename(app.WorkingDir))
    
    res = app(paths)
    
    info, aligned, struct_string = \
        list(MinimalRfamParser(res['Alignment'].readlines(),\
            seq_constructor=SEQ_CONSTRUCTOR_MAP[moltype]))[0]
    
    #Make new dict mapping original IDs
    new_alignment={}
    for k,v in aligned.NamedSeqs.items():
        new_alignment[int_keys.get(k,k)]=v
    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment,MolType=moltype)
    
    std_out = res['StdOut'].read()
    #clean up files
    res.cleanUp()
    for f in to_remove: remove(f)
    
    if return_stdout:
        return new_alignment, struct_string, std_out
    else:
        return new_alignment, struct_string
    

def cmalign_from_file(cm_file_path, seqs, moltype, alignment_file_path=None,\
    include_aln=False,return_stdout=False,params=None):
    """Uses cmalign to align seqs to alignment in cm_file_path.
        
        - cm_file_path: path to the file created by cmbuild, containing aligned
            sequences. This will be used to align sequences in seqs.
        - seqs: unaligned sequendes that are to be aligned to the sequences in
            cm_file.
        - moltype: cogent.core.moltype object.  Must be DNA or RNA
        - alignment_file_path: path to stockholm alignment file used to create
            cm_file.
            __IMPORTANT__: This MUST be the same file used by cmbuild
            originally.  Only need to pass in this file if include_aln=True.
            This helper function will NOT check if the alignment file is correct
            so you must use it correctly.
        - include_aln: Boolean to include sequences in aln_file in final
            alignment. (Default=False)
        - return_stdout: Boolean to return standard output from infernal.  This
            includes alignment and structure bit scores and average
            probabilities for each sequence. (Default=False)
    """
    #NOTE: Must degap seqs or Infernal well seg fault!
    seqs = SequenceCollection(seqs,MolType=moltype).degap()
    
    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seqs.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)
    
    if params is None:
        params = {}
    params.update({MOLTYPE_MAP[moltype]:True})
    
    app = Cmalign(InputHandler='_input_as_paths',WorkingDir='/tmp',\
        params=params)
    app.Parameters['--informat'].on('FASTA')
        
    #turn on --withali flag if True.
    if include_aln:
        if alignment_file_path is None:
            raise DataError, """Must have path to alignment file used to build CM if include_aln=True."""
        else:
            app.Parameters['--withali'].on(alignment_file_path)
                
    seqs_path = app._input_as_multiline_string(int_map.toFasta())
    paths = [cm_file_path,seqs_path]
    
    app.Parameters['-o'].on(get_tmp_filename(app.WorkingDir))
    res = app(paths)
    
    info, aligned, struct_string = \
        list(MinimalRfamParser(res['Alignment'].readlines(),\
            seq_constructor=SEQ_CONSTRUCTOR_MAP[moltype]))[0]
    
    
    #Make new dict mapping original IDs
    new_alignment={}
    for k,v in aligned.items():
        new_alignment[int_keys.get(k,k)]=v
    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment,MolType=moltype)
    std_out = res['StdOut'].read()
    res.cleanUp()
    if return_stdout:
        return new_alignment, struct_string, std_out
    else:
        return new_alignment, struct_string
    
def cmsearch_from_alignment(aln, structure_string, seqs, moltype, cutoff=0.0,\
    refine=False,params=None):
    """Uses cmbuild to build a CM file, then cmsearch to find homologs.
    
        - aln: an Alignment object or something that can be used to construct
            one.  All sequences must be the same length.
        - structure_string: vienna structure string representing the consensus
            stucture for the sequences in aln.  Must be the same length as the
            alignment.
        - seqs: SequenceCollection object or something that can be used to
            construct one, containing unaligned sequences that are to be
            searched.
        - moltype: cogent.core.moltype object.  Must be DNA or RNA
        - cutoff: bitscore cutoff.  No sequences < cutoff will be kept in
            search results. (Default=0.0).  Infernal documentation suggests
            a cutoff of log2(number nucleotides searching) will give most
            likely true homologs.
        - refine: refine the alignment and realign before building the cm.
            (Default=False)
    """
    #NOTE: Must degap seqs or Infernal well seg fault!
    seqs = SequenceCollection(seqs,MolType=moltype).degap()
    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seqs.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)
    
    cm_file, aln_file_string = cmbuild_from_alignment(aln, structure_string,\
        refine=refine,return_alignment=True)
    
    app = Cmsearch(InputHandler='_input_as_paths',WorkingDir='/tmp',\
        params=params)
    app.Parameters['--informat'].on('FASTA')
    app.Parameters['-T'].on(cutoff)
    
    to_remove = []
    
    seqs_path = app._input_as_multiline_string(int_map.toFasta())
    cm_path = app._tempfile_as_multiline_string(cm_file)
    paths = [cm_path,seqs_path]
    to_remove.append(cm_path)
    
    app.Parameters['--tabfile'].on(get_tmp_filename(app.WorkingDir))
    res = app(paths)
    
    search_results = list(CmsearchParser(res['SearchResults'].readlines()))
    if search_results:
        for i,line in enumerate(search_results):
            label = line[1]
            search_results[i][1]=int_keys.get(label,label)
    
    res.cleanUp()
    for f in to_remove:remove(f)
    
    return search_results

def cmsearch_from_file(cm_file_path, seqs, moltype, cutoff=0.0, params=None):
    """Uses cmbuild to build a CM file, then cmsearch to find homologs.
    
        - cm_file_path: path to the file created by cmbuild, containing aligned
            sequences. This will be used to search sequences in seqs.
        - seqs: SequenceCollection object or something that can be used to
            construct one, containing unaligned sequences that are to be
            searched.
        - moltype: cogent.core.moltype object.  Must be DNA or RNA
        - cutoff: bitscore cutoff.  No sequences < cutoff will be kept in
            search results. (Default=0.0).  Infernal documentation suggests
            a cutoff of log2(number nucleotides searching) will give most
            likely true homologs.
    """
    #NOTE: Must degap seqs or Infernal well seg fault!
    seqs = SequenceCollection(seqs,MolType=moltype).degap()
    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seqs.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)
    
    app = Cmsearch(InputHandler='_input_as_paths',WorkingDir='/tmp',\
        params=params)
    app.Parameters['--informat'].on('FASTA')
    app.Parameters['-T'].on(cutoff)
    
    seqs_path = app._input_as_multiline_string(int_map.toFasta())

    paths = [cm_file_path,seqs_path]
    
    app.Parameters['--tabfile'].on(get_tmp_filename(app.WorkingDir))
    res = app(paths)
    
    search_results = list(CmsearchParser(res['SearchResults'].readlines()))
    
    if search_results:    
        for i,line in enumerate(search_results):
            label = line[1]
            search_results[i][1]=int_keys.get(label,label)
    
    res.cleanUp()

    return search_results

