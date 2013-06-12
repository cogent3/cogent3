#! /usr/bin/env python
# file: flash.py

# Application controller for FLASh v1.2.6 
# Fast Length Adjustment of Short reads:
# http://ccb.jhu.edu/software/FLASH/

from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, ResultPath, \
    ApplicationError
from os import mkdir

__author__ = "Michael Robeson"
__copyright__ = "Copyright 2007-2013, The Cogent Project"
__credits__ = ["Michael Robeson"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Michael Robeson"
__email__ = "robesonms@ornl.gov"
__status__ = "Development"

class Flash(CommandLineApplication):
    """FLASh (v1.2.6) application controller for paired-end illumina data"""
    _command = 'flash'
    _parameters = {
    # Descriptions of parameters copied directly from 'flash -h'
    # and pasted below. NOTE: FLASh does not have flags for infiles.
    # These will be handled in separate convenience functions below.

    # -m, --min-overlap
    # The minimum required overlap length between two
    # reads to provide a confident overlap.  Default:
    # 10bp. 
    '-m':ValuedParameter(Prefix='-', Delimiter=' ', Name='m', Value='10'),
        
    # -M, --max-overlap
    # Maximum overlap length expected in approximately
    # 90% of read pairs.  It is by default set to 70bp,
    # which works well for 100bp reads generated from a
    # 180bp library, assuming a normal distribution of
    # fragment lengths.  Overlaps longer than the maximum
    # overlap parameter are still considered as good
    # overlaps, but the mismatch density (explained below)
    # is calculated over the first max_overlap bases in
    # the overlapped region rather than the entire
    # overlap.  Default: 70bp, or calculated from the
    # specified read length, fragment length, and fragment
    # length standard deviation.
    '-M':ValuedParameter(Prefix='-', Delimiter=' ', Name='M', Value='70'),

    # -x, --max-mismatch-density
    # Maximum allowed ratio between the number of
    # mismatched base pairs and the overlap length.
    # Two reads will not be combined with a given overlap
    # if that overlap results in a mismatched base density
    # higher than this value.  Note: Any occurence of an
    # 'N' in either read is ignored and not counted
    # towards the mismatches or overlap length.  Our
    # experimental results suggest that higher values of
    # the maximum mismatch density yield larger
    # numbers of correctly merged read pairs but at
    # the expense of higher numbers of incorrectly
    # merged read pairs.  Default: 0.25.
    '-x':ValuedParameter(Prefix='-', Delimiter=' ', Name='x', Value='0.25'),

    # -p, --phred-offset=OFFSET
    # The smallest ASCII value of the characters used to
    # represent quality values of bases in FASTQ files.
    # It should be set to either 33, which corresponds
    # to the later Illumina platforms and Sanger
    # platforms, or 64, which corresponds to the
    # earlier Illumina platforms.  Default: 33.
    '-p':ValuedParameter(Prefix='-', Delimiter=' ', Name='p', Value='33'),

    # -r, --read-len=LEN
    # -f, --fragment-len=LEN
    # -s, --fragment-len-stddev=LEN
    # Average read length, fragment length, and fragment
    # standard deviation.  These are convenience parameters
    # only, as they are only used for calculating the
    # maximum overlap (--max-overlap) parameter.
    # The maximum overlap is calculated as the overlap of
    # average-length reads from an average-size fragment
    # plus 2.5 times the fragment length standard
    # deviation.  The default values are -r 100, -f 180,
    # and -s 18, so this works out to a maximum overlap of
    # 70 bp.  If --max-overlap is specified, then the
    # specified value overrides the calculated value.
    #  If you do not know the standard deviation of the
    #  fragment library, you can probably assume that the
    #  standard deviation is 10% of the average fragment
    #  length.
    '-r':ValuedParameter(Prefix='-', Delimiter=' ', Name='r', Value='100'), 
    '-f':ValuedParameter(Prefix='-', Delimiter=' ', Name='f', Value='180'),
    '-s':ValuedParameter(Prefix='-', Delimiter=' ', Name='s', Value='18'),

    # --interleaved-input     
    # Instead of requiring files MATES_1.FASTQ and
    # MATES_2.FASTQ, allow a single file MATES.FASTQ that
    # has the paired-end reads interleaved.  Specify "-"
    # to read from standard input.
    '--interleaved-input':FlagParameter(Prefix='--', Name='interleaved-input'),    

    # --interleaved-output
    # Write the uncombined pairs in interleaved format.
    '--interleaved-output':FlagParameter(Prefix='--', Name='interleaved-output'),    
        
    # -I, --interleaved       
    # Equivalent to specifying both --interleaved-input
    # and --interleaved-output.
    '-I':FlagParameter(Prefix='-', Name='I'),
        
    # -o, --output-prefix=PREFIX
    #  Prefix of output files.  Default: "out".
    '-o':ValuedParamter(Prefix='-', Delimiter=' ', Name='o', Value='out'),
        
    # -d, --output-directory=DIR
    #  Path to directory for output files.  Default:
    #  current working directory.
    '-d':ValuedParameter(Prefix='-', Delimmiter=' ', Name='d'), #Value='./'),
        
    # -c, --to-stdout
    # Write the combined reads to standard output; do not
    # write uncombined reads to anywhere.
    '-c':FlagParameter(Prefix='-', Name='c'),
        
    #  -z, --compress
    #  Compress the FASTQ output files directly with zlib.
    #  Similar to specifying --compress-prog=gzip and
    #  --suffix=gz, but may be slightly faster.
    '-z':FlagParameter(Prefix='-', Name='z'),
        
    # --compress-prog=PROG    
    # Pipe the output through the compression program
    # PROG, which will be called as `PROG -c -',
    # plus any arguments specified by --compress-prog-args.
    # PROG must read uncompressed data from standard input
    # and write compressed data to standard output.
    # Examples: gzip, bzip2, xz, pigz.
    '--compress-prog':FlagProgram(Prefix='--', Name='compress-prog'),
        
    # --compress-prog-args=ARGS
    # A string of arguments that will be passed to the
    # compression program if one is specified with
    # --compress-prog.  Note: the argument -c is already
    # assumed.
    '--compress-prog-args':ValueParameter(Prefix='--', Delimiter=' ', Name='compress-prog-args'),
        
    # --suffix=SUFFIX, --output-suffix=SUFFIX
    # Use SUFFIX as the suffix of the output files
    # after ".fastq".  A dot before the suffix is assumed,
    # unless an empty suffix is provided.  Default:
    # nothing; or 'gz' if -z is specified; or PROG if
    # --compress-prog is specified.
    '--suffix':ValuedParameter(Prefix='--', Delimiter=' ', Name='suffix'),
    '--output-suffix':ValuedParameter(Prefix='--', Delimiter=' ', Name='output-suffix'),
    
    # -t, --threads=NTHREADS  
    # Set the number of worker threads.  This is in
    # addition to the I/O threads.  Default: number of
    # processors.  Note: if you need FLASH's output to
    # appear deterministically or in the same order as
    # the original reads, you must specify -t 1
    # (--threads=1).
    '-t':ValuedParameter(Prefix='-', Delimiter=' ', Name='t', Value='1'),
        
    # -q, --quiet
    # Do not print informational messages.  (Implied with
    # --to-stdout.)
    '-q':FlagParameter(Prefix='-', Name='q'),
        
    # -h, --help
    # Display this help and exit.
    '-h':FlagParameter(Prefix='-', Name='h'),
  
    # -v, --version
    # Display version.
    '-v':ValuedParameter(Prefix='-', Name='v')}

    _synonyms = {
    '--min-overlap':'-m',
    '--max-overlap':'-M',
    '--max-mismatch-density':'-x',
    '--phred-offset':'-p',
    '--read-len':'-r',
    '--fragment-len':'-f',
    '--fragment-len-stddev':'-s',
    '--interleaved':'-I',
    '--output-prefix':'-o',
    '--output-directory':'-d',
    '--to-stdout':'-c',
    '--compress':'-z',
    '--threads':'-t',
    '--quiet':'q',
    '--help':'-h',
    '--version':'-v'}


    _input_handler = '_input_as_paths'

    def _get_WorkingDir(self):
        """Gets the working directory"""
        return self._curr_working_dir

    def _set_WorkingDir(self, data):
       """Sets the working diretory"""
        self._curr_working_dir = self.Parameters['-d'] + '/'
        try: 
            mkdir(self.WorkingDir)
        except OSError: # Directory already exists
            pass
    
    WorkingDir = property(_get_WorkingDir,_set_WorkingDir)

    _suppress_stdout = True

    _suppress_stderr = True

    def _get_result_paths(self, data):
        """Captures FLASh output paths.
            FLASh defaults writing output to 3 files:
            - the assembled reads stored as *.extendedFrags.fastq
            - reads1 that failed to assemble as *.notCombined_1.fastq'
            - reads2 that failed to assemble as *.notCombined_2.fastq'
        """
        # NOTE: Need to add infile names to ResultPath.
        # However, there are no infile parameters. They need to be read in as
        # _input_as_paths within the convenience function below. I'm not 
        # sure how to capture outfile paths here as they are based on the
        # infile names.
        # I'll change these to be similar in syntax to 'muscle.py'
        result = {}
        result['assebled_reads'] = ResultPath(Path=self.WorkingDir+'.extendedFrags.fastq',\
            IsWritten=True)
        result['unassembled_reads1'] = ResultPath(Path=self.WorkingDir+'.notCombined_1.fastq',\
            IsWritten=True)
        result['unassembled_reads2'] = ResultPath(Path=self.WorkingDir+'.notCombined_2.fastq',\
            IsWritten=True)

    def getHelp(self):
        """FLASh (v1.2.6) description and help."""
        help_str =\
        """
        For basic help, type the following at the command line:
        'flash -h'
        
        Website:
            http://ccb.jhu.edu/software/FLASH/
        
        For questions / comments send e-mail to:
             flash.comment@gmail.com
        """
        return help_str



###################################################
# SOME FUNCTIONS TO EXECUTE THE MOST COMMON TASKS #
###################################################

def default_assemble(\  
    reads1_infile_path,
    reads2_infile_path,
    outfile_dir,
    output_prefix,
    read_length='100',
    frag_length='180',
    frag_std_dev='18',
    mis_match_density='0.25',
    min_overlap='10',
    max_overlap=None
    params={},
    WorkingDir=None,      # Do I need to set this given 'outfile_dir' above?
    SuppressStderr=None,  #   ''
    SuppressStdout=None): #   ''
    """Runs FLASh, with HISEQ default parameters to assemble paired-end reads.

        -reads1_infile_path : reads1.fastq infile path
        -reads2_infile_path : reads2.fastq infile path
        -outfile_dir : directory to write output
        -output_prefix : prefix to append to output files
        -read_length : average length of individual reads
        -frag_length : average length of assembled reads
        -frag_std_dev : fragment length standard deviation, ~ 10% of frag_length
        -mis_match_identity : max allowable ratio of mismatched bases and 
            overlap length. Reads above this value will not be assembled.
        -min_overlap : minimum allowable overlab to assemble reads
        -max_overlap : if set this will override the settings specified by 
            '-r','-s', and '-f'. These three parameters are used to dynamically
            calculate max_overlap when max_overlap is not provided.
        
        For HISEQ a good default 'max_overlap' would be '200'.
        For MISEQ use these parameters:
            read_length='250' frag_length='340' frag_std_dev='34'
            or: max_overlap='500'
    """
    
    # There are no input options for fastq infiles. So, we check if they exist
    # and store them as a list for later input via '_input_as_paths'
    # for the default '_input_handler'.
    if not exists(reads1_infile_path):
        raise IOError, 'File not found at: %s' % reads1_infile_path
    if not exists(reads2_infile_path):
        raise IOError, 'File not found at: %s' % reads2_infile_path

    infile_paths = [reads1_infile_path, reads2_infile_path]

    # required params
    params['-d'] = outfile_dir # set to absolut path!
    params['-o'] = output_prefix # string to append to all outfiles
    params['-x'] = mis_match_density
    params['-m'] = min_overlap

    # optional params
    if max_overlap == None:
        params['-r'] = read_length
        params['-f'] = frag_length
        params['-s'] = frag_std_dev
    else:
        params['-M'] = max_overlap

    # run assembler
    flash_app = Flash(\
        params=params,
        WorkingDir=WorkingDir,
        SuppressStderr=SuppressStderr,
        SuppressStdout=SuppressStdout)

    app_result = flash_app(infile_paths) # use default '_input_as_paths'
    #result = app_result
    #app_result.cleanUp()
    #return result




#if __name__ == "__main__":
 

  
