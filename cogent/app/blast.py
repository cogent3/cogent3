#!/usr/bin/env python
"""Application controllers for blast family
"""
from string import strip
from os import remove, access, F_OK, environ, path
from cogent.app.parameters import FlagParameter, ValuedParameter, MixedParameter
from cogent.app.util import CommandLineApplication, ResultPath, \
    get_tmp_filename, guess_input_handler, ApplicationNotFoundError
from cogent.parse.fasta import FastaFinder, LabeledRecordFinder, is_fasta_label
from cogent.parse.blast import LastProteinIds9, QMEBlast9, QMEPsiBlast9, BlastResult
from cogent.util.misc import app_path
from random import choice
from copy import copy

__author__ = "Micah Hamady"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Zongzhi Liu", "Micah Hamady", "Jeremy Widmann",
                    "Catherine Lozupone", "Rob Knight","Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Micah Hamady"
__email__ = "hamady@colorado.edu"
__status__ = "Prototype"

class Blast(CommandLineApplication):
    """BLAST generic application controller"""

    _common_options ={
        # defaults to non-redundant database
        #WARNING: This will only work if BLASTDB environment variable is set
        '-d':ValuedParameter('-',Name='d',Delimiter=' ', Value="nr"),
        
        # query file
        '-i':ValuedParameter('-',Name='i',Delimiter=' '),

        # Multiple Hits window size [Integer]
        '-A':ValuedParameter('-',Name='A',Delimiter=' '),

        # Threshold for extending hits [Integer] 
        '-f':ValuedParameter('-',Name='f',Delimiter=' '),

        # Expectation value (E) [Real] 
        '-e':ValuedParameter('-',Name='e',Delimiter=' ', Value="10.0"),

        # alignment view options: 
        # 0 = pairwise,
        # 1 = query-anchored showing identities,
        # 2 = query-anchored no identities,
        # 3 = flat query-anchored, show identities,
        # 4 = flat query-anchored, no identities,
        # 5 = query-anchored no identities and blunt ends,
        # 6 = flat query-anchored, no identities and blunt ends,
        # 7 = XML Blast output,
        # 8 = Tabular output,
        # 9 = Tabular output with comments
        # 10 = ASN, text
        # 11 = ASN, binary [Integer]
        '-m':ValuedParameter('-',Name='m',Delimiter=' ', Value="9"),

        # Output File for Alignment [File Out]  Optional
        '-o':ValuedParameter('-',Name='o',Delimiter=' '),

         # Filter query sequence with SEG [String] 
        '-F':ValuedParameter('-',Name='F',Delimiter=' '),

         #  Cost to open a gap [Integer] 
        '-G':ValuedParameter('-',Name='G',Delimiter=' '),

        # Cost to extend a gap [Integer] 
        '-E':ValuedParameter('-',Name='E',Delimiter=' '),

        # X dropoff value for gapped alignment (in bits) [Integer]
        # blastn 30, megablast 20, tblastx 0, all others 15 [Integer]
        '-X':ValuedParameter('-',Name='X',Delimiter=' '),

         # Show GI's in deflines [T/F] 
        '-I':ValuedParameter('-',Name='I',Delimiter=' '),
   
        # Number of database seqs to show one-line descriptionss for [Integer] 
        '-v':ValuedParameter('-',Name='v',Delimiter=' '),

        # Number of database sequence to show alignments for (B) [Integer] 
        '-b':ValuedParameter('-',Name='b',Delimiter=' '),

        # Perform gapped alignment (not available with tblastx) [T/F] 
        '-g':ValuedParameter('-',Name='g',Delimiter=' '),

        # Number of processors to use [Integer] 
        '-a':ValuedParameter('-',Name='a',Delimiter=' ', Value="1"),

        # Believe the query defline [T/F] 
        '-J':ValuedParameter('-',Name='J',Delimiter=' '),

        # SeqAlign file ('Believe the query defline' must be TRUE) [File Out]
        # Optional 
        '-O':ValuedParameter('-',Name='O',Delimiter=' '),

        # Matrix [String] 
        '-M':ValuedParameter('-',Name='M',Delimiter=' ', Value="BLOSUM62"),

        # Word size [Integer]  (blastn 11, megablast 28, all others 3)
        '-W':ValuedParameter('-',Name='W',Delimiter=' '),

        # Effective length of the database (use zero for the real size) [Real] 
        '-z':ValuedParameter('-',Name='z',Delimiter=' '),

        # Number of best hits from a region to keep [Integer] 
        '-K':ValuedParameter('-',Name='K',Delimiter=' '),

        # 0 for multiple hit, 1 for single hit [Integer] 
        '-P':ValuedParameter('-',Name='P',Delimiter=' '),

        # Effective length of the search space (use zero for real size) [Real] 
        '-Y':ValuedParameter('-',Name='Y',Delimiter=' '),

        # Produce HTML output [T/F] 
        '-T':ValuedParameter('-',Name='T',Delimiter=' ', Value="F"),

        # Restrict search of database to list of GI's [String]  Optional 
        '-l':ValuedParameter('-',Name='l',Delimiter=' '),

        # Use lower case filtering of FASTA sequence [T/F] Optional  
        '-U':ValuedParameter('-',Name='U',Delimiter=' '),

        # Dropoff (X) for blast extensions in bits (default if zero) [Real] 
        # blastn 20, megablast 10, all others 7
        '-y':ValuedParameter('-',Name='y',Delimiter=' '),

        # X dropoff value for final gapped alignment (in bits) [Integer] 
        # blastn/megablast 50, tblastx 0, all others 25
        '-Z':ValuedParameter('-',Name='Z',Delimiter=' '),

        # Input File for PSI-BLAST Restart [File In]  Optional 
        '-R':ValuedParameter('-',Name='R',Delimiter=' '),

    }

    _executable = 'blastall'

    _parameters = {}
    _parameters.update(_common_options)

    def __init__(self, cur_options, command, blast_mat_root=None,
                 extra_env="",
                 params=None,InputHandler=None,
                 SuppressStderr=None, SuppressStdout=None,WorkingDir=None,\
                 HALT_EXEC=False):
        """ Initialize blast """
        # update options
        self._parameters.update(cur_options)

        # check if need to set env variable (for cgi calls)
        if blast_mat_root:
            self._command = "export BLASTMAT=%s;%s%s" % (blast_mat_root, 
                                                    extra_env, command)
        else:
            # Determine if blast is installed and raise an ApplicationError 
            # if not -- this is done here so the user will get the most 
            # informative error message available.  
            self._error_on_missing_application(params)
                 
            # Otherwise raise error about $BLASTMAT not being set
            if not ('BLASTMAT' in environ or \
                    access(path.expanduser("~/.ncbirc"), F_OK) or \
                    access(".ncbirc", F_OK)):
                ## SHOULD THIS BE CHANGED TO RAISE AN ApplicationError?
                raise RuntimeError, blastmat_error_message
            self._command = command

        super(Blast, self).__init__(params=params,
                    InputHandler=InputHandler,SuppressStderr=SuppressStderr,
                    SuppressStdout=SuppressStdout,WorkingDir=WorkingDir,\
                    HALT_EXEC=HALT_EXEC)

    def _error_on_missing_application(self,params):
        """ Raise an ApplicationNotFoundError if the app is not accessible
        """
        if not app_path('blastall'):
            raise ApplicationNotFoundError,\
             "Cannot find blastall. Is it installed? Is it in your path?"

    def _input_as_seqs(self,data):
        lines = []
        for i,s in enumerate(data):
            #will number the sequences 1,2,3,etc.
            lines.append(''.join(['>',str(i+1)]))
            lines.append(s)
        return self._input_as_lines(lines)
        
    def _input_as_seq_id_seq_pairs(self,data):
        lines = []
        for seq_id,seq in data:
            lines.append(''.join(['>',str(seq_id)]))
            lines.append(seq)
        return self._input_as_lines(lines)

    def _input_as_lines(self,data):
        if data:
            self.Parameters['-i']\
                .on(super(Blast,self)._input_as_lines(data))

        return ''

    def _input_as_string(self,data):
        """Makes data the value of a specific parameter
    
        This method returns the empty string. The parameter will be printed
        automatically once set.
        """
        if data:
            self.Parameters['-i'].on(str(data))
        return ''

    def _input_as_multiline_string(self, data):
        if data:
            self.Parameters['-i']\
                .on(super(Blast,self)._input_as_multiline_string(data))
        return ''

    def _align_out_filename(self):

        if self.Parameters['-o'].isOn():
            aln_filename = self._absolute(str(self.Parameters['-o'].Value))
        else:
            raise ValueError, "No output file specified." 
        return aln_filename

    def _get_result_paths(self,data):
        
        result = {}
        if self.Parameters['-o'].isOn():
            out_name = self._align_out_filename()
            result['BlastOut'] = ResultPath(Path=out_name,IsWritten=True)
        return result

blastmat_error_message =\
"""BLAST cannot run if the BLASTMAT environment variable is not set.

Usually, the BLASTMAT environment variable points to the NCBI data directory,
which contains matrices like PAM30 and PAM70, etc.

Alternatively, you may create a .ncbirc file to define these variables.

From help file:

2) Create a .ncbirc file. In order for Standalone BLAST to operate, you
have will need to have a .ncbirc file that contains the following lines:

[NCBI] 
Data="path/data/"

Where "path/data/" is the path to the location of the Standalone BLAST
"data" subdirectory. For Example: 

Data=/root/blast/data

The data subdirectory should automatically appear in the directory where
the downloaded file was extracted. Please note that in many cases it may
be necessary to delimit the entire path including the machine name and
or the net work you are located on. Your systems administrator can help
you if you do not know the entire path to the data subdirectory.

Make sure that your .ncbirc file is either in the directory that you
call the Standalone BLAST program from or in your root directory.
"""

class PsiBlast(Blast):
    """PSI-BLAST application controller - Prototype"""
    _options ={

        # ASN.1 Scoremat input of checkpoint data: 
        # 0: no scoremat input
        # 1: Restart is from ASCII scoremat checkpoint file,
        # 2: Restart is from binary scoremat checkpoint file [Integer]  Optional
        '-q':ValuedParameter('-',Name='q',Delimiter=' '),

        # Output File for PSI-BLAST Matrix in ASCII [File Out] Optional 
        '-Q':ValuedParameter('-',Name='Q',Delimiter=' '),

        # Start of required region in query [Integer] 
        '-S':ValuedParameter('-',Name='S',Delimiter=' ', Value="1"),

        # ASN.1 Scoremat output of checkpoint data: 
        # 0: no scoremat output
        # 1: Output is ASCII scoremat checkpoint file (requires -J),
        # 2: Output is binary scoremat checkpoint file (requires -J) Optional
        '-u':ValuedParameter('-',Name='u',Delimiter=' '),

        # Cost to decline alignment (disabled when 0) [Integer] 
        '-L':ValuedParameter('-',Name='L',Delimiter=' ', Value="0"),

        # program option for PHI-BLAST [String] 
        '-p':ValuedParameter('-',Name='p',Delimiter=' ', Value="blastpgp"),

        # Use composition based statistics [T/F] 
        '-t':ValuedParameter('-',Name='t',Delimiter=' ', Value="T"),

        # Input Alignment File for PSI-BLAST Restart [File In] Optional 
        '-B':ValuedParameter('-',Name='B',Delimiter=' '),

        # Number of bits to trigger gapping [Real] 
        '-N':ValuedParameter('-',Name='N',Delimiter=' ', Value="22.0"),

        # End of required region in query (-1 indicates end of query) [Integer] 
        '-H':ValuedParameter('-',Name='H',Delimiter=' ', Value="-1"),

        # e-value threshold for inclusion in multipass model [Real] 
        '-h':ValuedParameter('-',Name='h',Delimiter=' ', Value="0.001"),

        # Constant in pseudocounts for multipass version [Integer] 
        '-c':ValuedParameter('-',Name='c',Delimiter=' ', Value="9"),

        # Maximum number of passes to use in  multipass version [Integer] 
        '-j':ValuedParameter('-',Name='j',Delimiter=' ', Value="1"),

        # Output File for PSI-BLAST Checkpointing [File Out]  Optional 
        '-C':ValuedParameter('-',Name='C',Delimiter=' '),

        # Compute locally optimal Smith-Waterman alignments [T/F] 
        '-s':ValuedParameter('-',Name='s',Delimiter=' ', Value="F"),

        # Hit File for PHI-BLAST [File In] 
        '-k':ValuedParameter('-',Name='k',Delimiter=' '),

    }

    def __init__(self, blast_mat_root=None, params=None,
                 extra_env="",
                 InputHandler=None,SuppressStderr=None,
                 SuppressStdout=None,WorkingDir=None,
                 HALT_EXEC=False):
        """ Initialize the Psi-Blast"""
        super(PsiBlast, self).__init__(self._options,
                    "blastpgp",
                    extra_env=extra_env,
                    blast_mat_root=blast_mat_root,
                    params=params,
                    InputHandler=InputHandler,SuppressStderr=SuppressStderr,
                    SuppressStdout=SuppressStdout,WorkingDir=WorkingDir,
                    HALT_EXEC=HALT_EXEC)


# should probably go into blastall superclass. it's late, works for now
BLASTALL_OPTIONS ={
        # Use lower case filtering of FASTA sequence [T/F] Optional  
        '-U':ValuedParameter('-',Name='U',Delimiter=' '),

        # Penalty for a nucleotide mismatch (blastn only) [Integer] 
        # default = -3 
        '-q':ValuedParameter('-',Name='q',Delimiter=' '),

        # Reward for a nucleotide match (blastn only) [Integer] 
        '-r':ValuedParameter('-',Name='r',Delimiter=' '),

        # Query Genetic code to use [Integer] default = 1
        '-Q':ValuedParameter('-',Name='Q',Delimiter=' '),

        # DB Genetic code (for tblast[nx] only) [Integer] 
        '-D':ValuedParameter('-',Name='D',Delimiter=' '),

        # Query strands to search against database (for blast[nx], and tblastx)
        # 3 is both, 1 is top, 2 is bottom [Integer] 
        '-S':ValuedParameter('-',Name='S',Delimiter=' '),

        # Program Name 
        '-p':ValuedParameter('-',Name='p',Delimiter=' '),

        # MegaBlast search [T/F] 
        '-n':ValuedParameter('-',Name='n',Delimiter=' '),

        # Location on query sequence [String]  Option 
        '-L':ValuedParameter('-',Name='L',Delimiter=' '),

        # Frame shift penalty (OOF algorithm for blastx) [Integer] 
        '-w':ValuedParameter('-',Name='w',Delimiter=' '),

        # Length of the largest intron allowed in tblastn for linking HSPs 
        #(0 disables linking) [Integer] 
        '-t':ValuedParameter('-',Name='t',Delimiter=' '),

        # Number of concatenated queries, for blastn and tblastn [Integer]
        '-B':ValuedParameter('-',Name='B',Delimiter=' '),
    }


class Blastall(Blast):
    """blastall application controller - Prototype """

    def __init__(self, blast_mat_root=None, params=None,
                 extra_env="",
                 InputHandler=None,SuppressStderr=None,
                 SuppressStdout=None,WorkingDir=None,
                 HALT_EXEC=False):
        """ Initialize the blastall"""
        super(Blastall, self).__init__(BLASTALL_OPTIONS,
                    "blastall", 
                    blast_mat_root=blast_mat_root,
                    extra_env=extra_env,
                    params=params,
                    InputHandler=InputHandler,SuppressStderr=SuppressStderr,
                    SuppressStdout=SuppressStdout,WorkingDir=WorkingDir,
                    HALT_EXEC=HALT_EXEC)
class MpiBlast(Blast):
    """mpblast application controller - Prototype """

    _mpi_options ={
        # Produces verbose debugging output for each node, optionally logs the 
        # output to a file
        '--debug':ValuedParameter('-',Name='--debug',Delimiter='='),

        # Set the scheduler process' MPI Rank (default is 1). Because the 
        # scheduler uses very little CPU it can be useful to force the 
        # scheduler to run on the same physical machine as the writer (rank 0).
        '--scheduler-rank':ValuedParameter('-',Name='--scheduler-rank',
                                           Delimiter='='),

        # Print the Altschul. et. al. 1997 paper reference instead of the 
        # mpiBLAST paper reference. With this option mpiblast output is nearly 
        # identical to NCBI-BLAST output.
        '--altschul-reference':FlagParameter(Prefix='--',
                                             Name='altschul-reference'),

        #Removes the local copy of the database from each node before 
        # terminating execution
        '--removedb':FlagParameter(Prefix='--', Name='removedb'),

        # Sets the method of copying files that each worker will use. 
        #  Default = "cp"
        # * cp : use standard file system "cp" command. 
        #        Additional option is --concurrent.
        # * rcp : use rsh "rcp" command. Additonal option is --concurrent.
        # * scp : use ssh "scp" command. Additional option is --concurrent.
        # * mpi : use MPI_Send/MPI_Recv to copy files. 
        #         Additional option is --mpi-size.
        # * none : do not copy files,instead use shared storage as local storage
        '--copy-via':ValuedParameter('-',Name='--copy-via', Delimiter='='),


        # set the number of concurrent accesses to shared storage. Default = 1
        '--concurrent':ValuedParameter('-',Name='--concurrent', Delimiter='='),

    
        # in bytes, set the maximum buffer size that MPI will use to send data 
        # when transferring files. Default = 65536
        '--mpi-size':ValuedParameter('-',Name='--mpi-size', Delimiter='='),


        # set whether file locking should be used to manage local fragment 
        # lists. Defaults to off. When --concurrency > 1 defaults to on
        # [on|off]
        '--lock':ValuedParameter('-',Name='--lock', Delimiter='='),

        # When set, the writer will use the database on shared storage for 
        # sequence lookup. Can drastically reduce overhead for some blastn 
        # searches.
        '--disable-mpi-db':FlagParameter(Prefix='--', Name='disable-mpi-db'),

        # Under unix, sets the nice value for each mpiblast process.
        '--nice':ValuedParameter('-',Name='--nice', Delimiter='='),

        # Under unix, sets the nice value for each mpiblast process.
        '--config-file':ValuedParameter('--',Name='config-file', Delimiter='='),


        # Experimental. When set, mpiblast will read the output file and 
        # attempt to continue a previously aborted run where it left off
        '--resume-run':FlagParameter(Prefix='--', Name='resume-run'),

        # print the mpiBLAST version
        '--version':FlagParameter(Prefix='--', Name='version'),
    }

    _mpi_options.update(BLASTALL_OPTIONS)

    def __init__(self, blast_mat_root=None, params=None,
                 mpiblast_root="/usr/local/bin/",
                 local_root="/var/scratch/mpiblastdata/",
                 shared_root="/quicksand/hamady/data/blast/mpidb/",
                 config_file="/quicksand2/downloads2/mpiblast/mpiblast.conf",
                 num_db_frags=40,
                 InputHandler=None,SuppressStderr=None,
                 SuppressStdout=None,WorkingDir=None,
                 HALT_EXEC=False):
        """ Initialize mpiblast"""
        if config_file:
            params["--config-file"] = config_file
        super(MpiBlast, self).__init__(self._mpi_options,
                    "mpirun -np %d %smpiblast" % ((num_db_frags + 2),  
                                                    mpiblast_root),
                    blast_mat_root=blast_mat_root,
                    extra_env="export Local=%s; export Shared=%s;" %(local_root,
                        shared_root),
                    params=params,
                    InputHandler=InputHandler,SuppressStderr=SuppressStderr,
                    SuppressStdout=SuppressStdout,WorkingDir=WorkingDir,
                    HALT_EXEC=HALT_EXEC)

class FastaCmd(CommandLineApplication):
    """FastaCmd application controller - Prototype"""

    _options ={
        # Database [String]  Optional 
        '-d':ValuedParameter('-',Name='d',Delimiter=' '),

        # Type of file 
        # G - guess mode (look for protein, then nucleotide)
        # T - protein
        # F - nucleotide [String]  Optional
        '-p':ValuedParameter('-',Name='p',Delimiter=' ', Value="G"),

        # Search str: GIs, accessions and loci may be used delimited by comma 
        '-s':ValuedParameter('-',Name='s',Delimiter=' '),

        # Input file wilth GIs/accessions/loci for batch retrieval Optional 
        '-i':ValuedParameter('-',Name='i',Delimiter=' '),

        # Retrieve duplicate accessions [T/F]  Optional 
        '-a':ValuedParameter('-',Name='a',Delimiter=' ', Value='F'),

        # Line length for sequence [Integer]  Optional 
        '-l':ValuedParameter('-',Name='l',Delimiter=' '),

        # Definition line should contain target gi only [T/F]  Optional 
        '-t':ValuedParameter('-',Name='t',Delimiter=' '),

        # Output file [File Out]  Optional
        '-o':ValuedParameter('-',Name='o',Delimiter=' '),

        # Use Ctrl-A's as non-redundant defline separator [T/F]  Optional 
        '-c':ValuedParameter('-',Name='c',Delimiter=' '),

        # Dump the entire database in fasta format [T/F]  Optional 
        '-D':ValuedParameter('-',Name='D',Delimiter=' '),

        # Range of sequence to extract (Format: start,stop) 
        # 0 in 'start' refers to the beginning of the sequence
        # 0 in 'stop' refers to the end of the sequence [String]  Optional
        '-L':ValuedParameter('-',Name='L',Delimiter=' '),

        # Strand on subsequence (nucleotide only): 1 is top, 2 is bottom [Int] 
        '-S':ValuedParameter('-',Name='S',Delimiter=' '),

        # Print taxonomic information for requested sequence(s) [T/F] 
        '-T':ValuedParameter('-',Name='T',Delimiter=' '),

        # Print database information only (overrides all other options) [T/F]
        '-I':ValuedParameter('-',Name='I',Delimiter=' '),

        #  Retrieve sequences with this PIG [Integer]  Optional 
        '-P':ValuedParameter('-',Name='P',Delimiter=' '),
        
    }
    _parameters = {}
    _parameters.update(_options)
    _command = 'fastacmd'
  
    def _input_as_lines(self,data):
        if data:
            self.Parameters['-i']\
                .on(super(FastaCmd,self)._input_as_lines(data))
        return ''

    def _input_as_seqs(self,data):
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
            self.Parameters['-s'].on(data)
        return ''

    def _out_filename(self):

        if self.Parameters['-o'].isOn():
            aln_filename = self._absolute(str(self.Parameters['-o'].Value))
        else:
            raise ValueError, "No output file specified." 
        return aln_filename

    def _get_result_paths(self,data):
        
        result = {}
        if self.Parameters['-o'].isOn():
            out_name = self._out_filename()
            result['FastaOut'] = ResultPath(Path=out_name,IsWritten=True)
        return result

def seqs_to_stream(seqs, ih):
    """Converts seqs into stream of FASTA records, depending on input handler.

    Each FASTA record will be a list of lines.
    """
    if ih == '_input_as_multiline_string':
        recs = FastaFinder(seqs.split('\n'))
    elif ih == '_input_as_string':
        recs = FastaFinder(open(seqs))
    elif ih == '_input_as_seqs':
        recs = [['>'+str(i), s] for i, s in enumerate(seqs)]
    elif ih == '_input_as_lines':
        recs = FastaFinder(seqs)
    else:
        raise TypeError, "Unknown input handler %s" % ih
    return recs
 
#SOME FUNCTIONS TO EXECUTE THE MOST COMMON TASKS
def blast_seqs(seqs,
                 blast_constructor,
                 blast_db=None,
                 blast_mat_root=None,
                 params={},
                 add_seq_names=True,
                 out_filename=None,
                 WorkingDir=None,
                 SuppressStderr=None,
                 SuppressStdout=None,
                 input_handler=None,
                 HALT_EXEC=False
                 ):
    """Blast list of sequences.

    seqs: either file name or list of sequence objects or list of strings or
    single multiline string containing sequences.
    
    WARNING: DECISION RULES FOR INPUT HANDLING HAVE CHANGED. Decision rules 
    for data are as follows. If it's s list, treat as lines, unless 
    add_seq_names is true (in which case treat as list of seqs). If it's a
    string, test whether it has newlines. If it doesn't have newlines, assume
    it's a filename. If it does have newlines, it can't be a filename, so
    assume it's a multiline string containing sequences.

    If you want to skip the detection and force a specific type of input
    handler, use input_handler='your_favorite_handler'.
   
    add_seq_names: boolean. if True, sequence names are inserted in the list
        of sequences. if False, it assumes seqs is a list of lines of some
        proper format that the program can handle
    """

    # set num keep
    
    if blast_db:
        params["-d"] = blast_db

    if out_filename:
        params["-o"] = out_filename

    ih = input_handler or guess_input_handler(seqs, add_seq_names)
           
    blast_app = blast_constructor(
                   params=params,
                   blast_mat_root=blast_mat_root,
                   InputHandler=ih,
                   WorkingDir=WorkingDir,
                   SuppressStderr=SuppressStderr,
                   SuppressStdout=SuppressStdout,
                   HALT_EXEC=HALT_EXEC)

    return blast_app(seqs)


def fasta_cmd_get_seqs(acc_list,
                 blast_db=None,
                 is_protein=None,
                 out_filename=None,
                 params={},
                 WorkingDir="/tmp",
                 SuppressStderr=None,
                 SuppressStdout=None):
    """Retrieve sequences for list of accessions """

    if is_protein is None:
        params["-p"] = 'G' 
    elif is_protein:
        params["-p"] = 'T' 
    else:
        params["-p"] = 'F' 

    if blast_db:
        params["-d"] = blast_db

    if out_filename:
        params["-o"] = out_filename

    # turn off duplicate accessions
    params["-a"] = "F"

    # create Psi-BLAST
    fasta_cmd = FastaCmd(params=params,
                       InputHandler='_input_as_string',
                       WorkingDir=WorkingDir,
                       SuppressStderr=SuppressStderr,
                       SuppressStdout=SuppressStdout)

    # return results
    return fasta_cmd("\"%s\"" % ','.join(acc_list))

def fastacmd_is_crap(line):
    """Handles missing ids..."""
    return (not line) or line.isspace() or line.startswith('[')

FastaCmdFinder = LabeledRecordFinder(is_fasta_label, ignore=fastacmd_is_crap)

def seqs_from_fastacmd(acc_list, blast_db,is_protein=True):
    """Get dict of description:seq from fastacmd."""
    fasta_cmd_res = fasta_cmd_get_seqs(acc_list, blast_db=blast_db, \
        is_protein=is_protein)
    recs = FastaCmdFinder(fasta_cmd_res['StdOut'])
    result = {}
    for rec in recs:
        try:
            result[rec[0][1:].strip()] = ''.join(map(strip, rec[1:]))
        except IndexError:  #maybe we didn't get a sequence?
            pass
    fasta_cmd_res.cleanUp()
    return result

def psiblast_n_neighbors(seqs,
                 n=100,
                 blast_db=None,
                 core_threshold=1e-50,
                 extra_threshold=1e-10,
                 lower_threshold=1e-6,
                 step=100,
                 method="two-step",
                 blast_mat_root=None,
                 params={},
                 add_seq_names=False,
                 WorkingDir=None,
                 SuppressStderr=None,
                 SuppressStdout=None,
                 input_handler=None,
                 scorer=3,   #shotgun with 3 hits needed to keep
                 second_db=None
                 ):
    """PsiBlasts sequences, stopping when n neighbors are reached.

    core_threshold: threshold for the core profile (default: 1e-50)
    extra_threshold: threshold for pulling in additional seqs (default:1e-10)
    lower_threshold: threshold for seqs in final round (default:1e-6)

    seqs: either file name or list of sequence objects or list of strings or
    single multiline string containing sequences.
    If you want to skip the detection and force a specific type of input
    handler, use input_handler='your_favorite_handler'.
   
    add_seq_names: boolean. if True, sequence names are inserted in the list
        of sequences. if False, it assumes seqs is a list of lines of some
        proper format that the program can handle
    """
    if blast_db:
        params["-d"] = blast_db
    
    ih = input_handler or guess_input_handler(seqs, add_seq_names)
    recs = seqs_to_stream(seqs, ih) #checkpointing can only handle one seq...
    
    #set up the parameters for the core and additional runs
    max_iterations = params['-j']
    params['-j'] = 2    #won't checkpoint with single iteration
    
    app = PsiBlast(params=params,
                   blast_mat_root=blast_mat_root,
                   InputHandler='_input_as_lines',
                   WorkingDir=WorkingDir,
                   SuppressStderr=SuppressStderr,
                   SuppressStdout=SuppressStdout,
                   )
    result = {}
    for seq in recs:
        query_id = seq[0][1:].split(None,1)[0]
        if method == "two-step":
            result[query_id] = ids_from_seq_two_step(seq, n, max_iterations, \
                app, core_threshold, extra_threshold, lower_threshold, second_db)
        elif method == "lower_threshold":
            result[query_id] = ids_from_seq_lower_threshold(seq, n, \
                max_iterations, app, core_threshold, lower_threshold, step)
        elif method == "iterative":
            result[query_id] = ids_from_seqs_iterative(seq, app, \
               QMEPsiBlast9, scorer, params['-j'], n)
        else:
            raise TypeError, "Got unknown method %s" % method
            
    params['-j'] = max_iterations
    return result

def ids_from_seq_two_step(seq, n, max_iterations, app, core_threshold, \
    extra_threshold, lower_threshold, second_db=None):
    """Returns ids that match a seq, using a 2-tiered strategy.
    
    Optionally uses a second database for the second search.
    """
    #first time through: reset 'h' and 'e' to core
    #-h is the e-value threshold for including seqs in the score matrix model
    app.Parameters['-h'].on(core_threshold)
    #-e is the e-value threshold for the final blast
    app.Parameters['-e'].on(core_threshold)
    checkpoints = []
    ids = []
    last_num_ids = None
    for i in range(max_iterations):
        if checkpoints:
            app.Parameters['-R'].on(checkpoints[-1])
        curr_check = 'checkpoint_%s.chk' % i
        app.Parameters['-C'].on(curr_check)

        output = app(seq)
        #if we didn't write a checkpoint, bail out
        if not access(curr_check, F_OK):
            break
        #if we got here, we wrote a checkpoint file
        checkpoints.append(curr_check)
        result = list(output.get('BlastOut', output['StdOut']))
        output.cleanUp()
        if result:
            ids = LastProteinIds9(result,keep_values=True,filter_identity=False)
        num_ids = len(ids)
        if num_ids >= n:
            break
        if num_ids == last_num_ids:
            break
        last_num_ids = num_ids

    #if we didn't write any checkpoints, second run won't work, so return ids
    if not checkpoints:
        return ids

    #if we got too many ids and don't have a second database, return the ids we got
    if (not second_db) and num_ids >= n:
        return ids
    
    #second time through: reset 'h' and 'e' to get extra hits, and switch the
    #database if appropriate
    app.Parameters['-h'].on(extra_threshold)
    app.Parameters['-e'].on(lower_threshold)
    if second_db:
        app.Parameters['-d'].on(second_db)
    for i in range(max_iterations): #will always have last_check if we get here
        app.Parameters['-R'].on(checkpoints[-1])
        curr_check = 'checkpoint_b_%s.chk' % i
        app.Parameters['-C'].on(curr_check)
        output = app(seq)
        #bail out if we couldn't write a checkpoint
        if not access(curr_check, F_OK):
            break
        #if we got here, the checkpoint worked
        checkpoints.append(curr_check)
        result = list(output.get('BlastOut', output['StdOut']))
        if result:
            ids = LastProteinIds9(result,keep_values=True,filter_identity=False)
        num_ids = len(ids)
        if num_ids >= n:
            break
        if num_ids == last_num_ids:
            break
        last_num_ids = num_ids
    #return the ids we got. may not be as many as we wanted.
    for c in checkpoints:
        remove(c)
    return ids
    
class ThresholdFound(Exception): pass

def ids_from_seq_lower_threshold(seq, n, max_iterations, app, core_threshold, \
    lower_threshold, step=100):
    """Returns ids that match a seq, decreasing the sensitivity."""
    last_num_ids = None
    checkpoints = []
    cp_name_base = make_unique_str()

    # cache ides for each iteration
    # store { iteration_num:(core_threshold, [list of matching ids]) }
    all_ids = {}
    try:
        i=0
        while 1:
            #-h is the e-value threshold for inclusion in the score matrix model
            app.Parameters['-h'].on(core_threshold)
            app.Parameters['-e'].on(core_threshold)
            if core_threshold > lower_threshold:
                raise ThresholdFound
            if checkpoints:
                #-R restarts from a previously stored file
                app.Parameters['-R'].on(checkpoints[-1])
            #store the score model from this iteration
            curr_check = 'checkpoint_' + cp_name_base + '_' + str(i) + \
                    '.chk'
            app.Parameters['-C'].on(curr_check)
            output = app(seq)
            result = list(output.get('BlastOut', output['StdOut']))
            #sometimes fails on first try -- don't know why, but this seems
            #to fix problem
            while not result:
                output = app(seq)
                result = list(output.get('BlastOut', output['StdOut']))

            ids = LastProteinIds9(result,keep_values=True,filter_identity=False)
            output.cleanUp()
            all_ids[i + 1] = (core_threshold, copy(ids))
            if not access(curr_check, F_OK):
                raise ThresholdFound
            checkpoints.append(curr_check)
            num_ids = len(ids)
            if num_ids >= n:
                raise ThresholdFound
            last_num_ids = num_ids
            core_threshold *= step
            if i >= max_iterations - 1: #because max_iterations is 1-based
                raise ThresholdFound
            i += 1
    except ThresholdFound:
        for c in checkpoints:
            remove(c)
        #turn app.Parameters['-R'] off so that for the next file it does not
        #try and read in a checkpoint file that is not there
        app.Parameters['-R'].off()
        return ids, i + 1, all_ids

def make_unique_str(num_chars=20):
    """make a random string of characters for a temp filename"""
    chars = 'abcdefghigklmnopqrstuvwxyz'
    all_chars = chars + chars.upper() + '01234567890'
    picks = list(all_chars)
    return ''.join([choice(picks) for i in range(num_chars)])

def make_subject_match_scorer(count):
    def subject_match_scorer(checked_ids):
        """From {subject:{query:score}} returns subject ids w/ >= count hits.
        
        Useful for elminating subjects with few homologs.
        """
        return [key for key, val in checked_ids.items() if len(val) >= count]
    return subject_match_scorer

def make_shotgun_scorer(count):
    def shotgun_scorer(checked_ids):
        """From {subject:{query:score}} returns any ids w/ >= count hits.
        
        A hit counts towards a sequence's score if it was either the subject
        or the query, but we don't double-count (subject, query) pairs, i.e.
        if A hits B and B hits A, only one (A,B) hit will be counted, although
        it will be counted as both (A,B) and (B,A) (i.e. it will help preserve
        both A and B).
        """
        result = {}
        for subject, val in checked_ids.items():
            for query in val.keys():
                if subject not in result:
                    result[subject] = {}
                result[subject][query] = True
                if query not in result:
                    result[query] = {}
                result[query][subject] = True
        return [key for key, val in result.items() if len(val) >= count]
    return shotgun_scorer

def keep_everything_scorer(checked_ids):
    """Returns every query and every match in checked_ids, with best score."""
    result = checked_ids.keys()
    for i in checked_ids.values():
        result.extend(i.keys())
    return dict.fromkeys(result).keys()

def ids_from_seqs_iterative(seqs, app, query_parser, \
    scorer=keep_everything_scorer, max_iterations=None, blast_db=None,\
    max_seqs=None, ):
    """Gets the ids from each seq, then does each additional id until all done.

    If scorer is passed in as an int, uses shotgun scorer with that # hits.
    """
    if isinstance(scorer, int):
        scorer = make_shotgun_scorer(scorer)
    seqs_to_check = list(seqs)
    checked_ids = {}
    curr_iteration = 0
    while seqs_to_check:
        unchecked_ids = {}
        #pass seqs to command
        all_output = app(seqs_to_check)
        output = all_output.get('BlastOut', all_output['StdOut'])

        for query_id, match_id, match_score in query_parser(output):
            if query_id not in checked_ids:
                checked_ids[query_id] = {}
            checked_ids[query_id][match_id] = match_score
            if match_id not in checked_ids:
                unchecked_ids[match_id] = True
        all_output.cleanUp()
        if unchecked_ids:
            seq_file = fasta_cmd_get_seqs(unchecked_ids.keys(),
                app.Parameters['-d'].Value)['StdOut']
            seqs_to_check = []
            for s in FastaCmdFinder(fasta_cmd_get_seqs(\
                unchecked_ids.keys(), app.Parameters['-d'].Value)['StdOut']):
                seqs_to_check.extend(s)
        else:
            seqs_to_check = []
        #bail out if max iterations or max seqs was defined and we've reached it
        curr_iteration += 1
        if max_iterations and (curr_iteration >= max_iterations):
            break
        if max_seqs:
            curr = scorer(checked_ids)
            if len(curr) >= max_seqs:
                return curr
    return scorer(checked_ids)  #scorer should return list of good ids


def blastp(seqs, blast_db="nr", e_value="1e-20", max_hits=200, 
           working_dir="/tmp", blast_mat_root=None, extra_params={}):
    """
    Returns BlastResult from input seqs, using blastp.
    
    Need to add doc string   
    """

    # set up params to use with blastp
    params = {
        # matrix
        "-M":"BLOSUM62",

        # max procs
        "-a":"1",

        # expectation
        "-e":e_value,

        # max seqs to show
        "-b":max_hits,

        # max one line descriptions
        "-v":max_hits,

        # program
        "-p":"blastp"
    }
    params.update(extra_params)
 
    # blast
    blast_res =  blast_seqs(seqs, 
        Blastall,
        blast_mat_root=blast_mat_root,
        blast_db=blast_db,
        params=params, 
        add_seq_names=False,
        WorkingDir=working_dir
        )

    # get prot id map
    if blast_res['StdOut']:
        lines = [x for x in blast_res['StdOut']]
        return BlastResult(lines)

    return None 

def blastn(seqs, blast_db="nt", e_value="1e-20", max_hits=200, 
           working_dir="/tmp", blast_mat_root=None, extra_params={}):
    """
    Returns BlastResult from input seqs, using blastn.
    
    Need to add doc string   
    """

    # set up params to use with blastp
    params = {
        # matrix
        "-M":"BLOSUM62",

        # max procs
        "-a":"1",

        # expectation
        "-e":e_value,

        # max seqs to show
        "-b":max_hits,

        # max one line descriptions
        "-v":max_hits,

        # program
        "-p":"blastn"
    }
    params.update(extra_params)
 
    # blast
    blast_res =  blast_seqs(seqs, 
        Blastall,
        blast_mat_root=blast_mat_root,
        blast_db=blast_db,
        params=params, 
        add_seq_names=False,
        WorkingDir=working_dir
        )

    # get prot id map
    if blast_res['StdOut']:
        lines = [x for x in blast_res['StdOut']]
        return BlastResult(lines)

    return None 



def blastx(seqs, params=None):
    """Returns BlastResults from input seqs, using blastx."""
    raise NotImplementedError

def tblastx(seqs, params=None):
    """Returns BlastResults from input seqs, using tblastx."""
    raise NotImplementedError

def psiblast(seqs, params=None):
    """Returns BlastResults from input seqs, using psiblast."""
    raise NotImplementedError

def reciprocal_best_blast_hit(query_id, db_1, db_2, exclude_self_hits=True,\
    params=None):
    """Returns best hit in db_2 that maps back to query_id in db_1, or None.
    
    exclude_self_hits: if True (the default), returns the best hit that 
    doesn't have the same id. Otherwise, will return the same id if it is in
    both databases (assuming it's the same sequence in both).
    """
    raise NotImplementedError

    #make with factory functions for the blast hits


if __name__ == "__main__":

    print "Debug. examples of how i've been using."
    
    print "Example of straightforward BLAST"

# WARNING: I changed a bunch of stuff to make testing easier, since nr doesn't
# fit in memory on my laptop. I created a database 'eco' using formatdb on the
# E. coli K12 fasta file from this URL:
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K12/NC_000913.faa
# Because we're blasting an archaeal sequence against one bacterial genome, I
# relaxed the inclusion thresholds substantially. DO NOT USE THESE AGAINST NR!

    in_filename = "test_seq.fasta" 
    out_filename = "test.out" 
    # if blast env variable set, can just say 'nr'
    #BLAST_DB = "/home/hamady/quicksand/data/blast/db/nr"  
    BLAST_DB = 'nr' #'nr'
    BLAST_MAT_ROOT="/home/hamady/apps/blast-2.2.9/data"
    #BLAST_MAT_ROOT='/Users/rob/ncbi/data'
    # set up params to use with iterative 

    #print seqs_from_fastacmd(['16766313'], 'nr', True)
    #raise ValueError, "dbug"
    params = {

        # matrix
        "-M":"PAM70",
        # max procs 
        "-a":2,
         # expect 
        "-e":1e-15,
 
# blastall  
#        # program   
#        "-p":"blastp",
 
# psi-blast
        # max iterations 
        "-j":2,
       
        # max seqs to show 
        "-b":50,
         # inclusion 
        "-h":1e-2,
    }
   
    in_seqs = """>stm:STMabcdef  thrA; aspartokinase I 
    MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTIGGQDA
    LPNISDAERIFSDLLAGLASAQPGFPLARLKMVVEQEFAQIKHVLHGISLLGQCPDSINA
    ALICRGEKMSIAIMAGLLEARGHRVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASQIP
    ADHMILMAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQV
    PDARLLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASDS
    DDNLPVKGISNLNNMAMFSVSGPGMKGMIGMAARVFAAMSRAGISVVLITQSSSEYSISF
    CVPQSDCARARRAMQDEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAAL
    ARANINIVAIAQGSSERSISVVVNNDDATTGVRVTHQMLFNTDQVIEVFVIGVGGVGGAL"""

# The following should now give the same output:
#
#    in_seqs = 'tiny.faa'    #tiny.faa in cwd contains the sequence above
#
#    in_seqs = """>gi|2501594|sp|Q57997|Y577_METJA PROTEIN MJ0577
#MSVMYKKILYPTDFSETAEIALKHVKAFKTLKAEEVILLHVIDEREIKKRDIFSLLLGVAGLNKSVEEFE
#NELKNKLTEEAKNKMENIKKELEDVGFKVKDIIVVGIPHEEIVKIAEDEGVDIIIMGSHGKTNLKEILLG
#SVTENVIKKSNKPVLVVKRKNS""".split()  #lines instead of multiline string
#   
    blast_res =  blast_seqs(in_seqs, Blastall,
                        blast_mat_root=BLAST_MAT_ROOT,
                        add_seq_names=False,
                        blast_db=BLAST_DB, 
                        params={'-p': 'blastp','-e': '1','-m': 9},
                        out_filename=out_filename)
                            
    print [x for x in blast_res['StdOut']]
    print [x for x in blast_res['StdErr']]
    print blast_res
    #for x in blast_res['BlastOut']:
    #    print x.rstrip()
    blast_res.cleanUp()
    #print '\n\n'
    #print "Example of psiblast_n_neighbors"
    #print "Method 1: two-step with high- and low-confidence matches"
    #print psiblast_n_neighbors(in_seqs, n=10, blast_db=BLAST_DB, \
    #    method="two-step", blast_mat_root=BLAST_MAT_ROOT,params=params,\
    #    core_threshold=1e-5, extra_threshold=1e-2, lower_threshold=1e-1)
    #print
    #print "Method 2: keep lowering threshold"
    #print psiblast_n_neighbors(in_seqs, n=10, blast_db=BLAST_DB, \
    #    method="lower_threshold", blast_mat_root=BLAST_MAT_ROOT,params=params,
    #    core_threshold=1e-6, lower_threshold=1e-2)
    #print
    #print "Method 3: psi-blast shotgun"
    #print psiblast_n_neighbors(in_seqs, n=10, blast_db=BLAST_DB, \
    #    method="iterative", blast_mat_root=BLAST_MAT_ROOT,params=params,
    #    core_threshold=1e-5, lower_threshold=1e-2)
    #print 
    #print "Method 4: two-step with high- and low-confidence matches, diff dbs"
    #print psiblast_n_neighbors(in_seqs, n=10, blast_db=BLAST_DB, \
    #    method="two-step", blast_mat_root=BLAST_MAT_ROOT,params=params,\
    #    core_threshold=1e-5, extra_threshold=1e-2, lower_threshold=1e-1, second_db='stm')
    #print

