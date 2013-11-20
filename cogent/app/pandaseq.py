#!/usr/bin/env python
# file: pandaseq.py

# Application controller for pandaseq (v2.4) 
# https://github.com/neufeld/pandaseq
# 

from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, ResultPath, \
    ApplicationError
from os.path import exists, abspath, join
import tempfile

__author__ = "Michael Robeson"
__copyright__ = "Copyright 2007-2013, The Cogent Project"
__credits__ = ["Michael Robeson"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Michael Robeson"
__email__ = "robesonms@ornl.gov"
__status__ = "Development"


class PandaSeq(CommandLineApplication):
    """pandaseq application controller for joining paired-end reads """
    _command = 'pandaseq'
    _parameters = {
    # pandaseq 2.4 <andre@masella.name>
    # Usage: pandaseq -f forward.fastq -r reverse.fastq [-6] [-a] [-B] 
    #    [-C module1 -C module2 ...] [-d flags] [-F] [-j] [-L maxlen] 
    #    [-l minlen] [-N] [-o minoverlap] [-p forwardprimer] 
    #    [-q reverseprimer] [-T threads] [-t threshold] > assembled.fastq

   # -6  Use PHRED+64 (CASAVA 1.3-1.7) instead of PHRED+33 (CASAVA 1.8+).
    '-6':FlagParameter(Prefix='-', Name='6'),

    # -a  Strip the primers after assembly, rather than before.
    '-a':FlagParameter(Prefix='-', Name='a'),

    # -B  Allow unbarcoded sequences (try this for BADID errors).
    '-B':FlagParameter(Prefix='-', Name='B'),

    # -C module   Load a sequence validation module.
    '-C':FlagParameter(Prefix='-', Name='C'),

    # -d flags    Control the logging messages. Capital to enable; small to disable.
    #    (R)econstruction detail.
    #    Sequence (b)uilding information.
    #    (F)ile processing.
    #    (k)-mer table construction.
    #    Show every (m)ismatch.
    '-d':ValuedParameter(Prefix='-', Delimiter=' ', Name='d'),

    #    Optional (s)tatistics.
    # -f  Input FASTQ file containing forward reads.
    '-f':ValuedParameter(Prefix='-', Delimiter=' ', Name='f'),

    # -F  Output FASTQ instead of FASTA.
    '-F':FlagParameter(Prefix='-', Name='F'),

    # -j  Input files are bzipped.
    '-j':FlagParameter(Prefix='-', Name='j'),

    # -k kmers    The number of k-mers in the table.
    '-k':ValuedParameter(Prefix='-', Delimiter=' ', Name='k'),

    # -L maxlen   Maximum length for a sequence
    '-L':ValuedParameter(Prefix='-', Delimiter=' ', Name='L'),

    # -l minlen   Minimum length for a sequence
    '-l':ValuedParameter(Prefix='-', Delimiter=' ', Name='l'),

    # -N  Eliminate all sequences with unknown nucleotides in the output.
    '-N':FlagParameter(Prefix='-', Name='N'),

    # -o minoverlap   Minimum overlap between forward and reverse reads (default = 1)
    '-o':ValuedParameter(Prefix='-', Delimiter=' ', Name='o'),

    # -p  Forward primer sequence or number of bases to be removed.
    '-p':ValuedParameter(Prefix='-', Delimiter=' ', Name='p'),

    # -q  Reverse primer sequence or number of bases to be removed.
    '-q':ValuedParameter(Prefix='-', Delimiter=' ', Name='q'),

    # -r  Input FASTQ file containing reverse reads.
    '-r':ValuedParameter(Prefix='-', Delimiter=' ', Name='r'),

    # -T thread   Run with a number of parallel threads.
    '-T':ValuedParameter(Prefix='-', Delimiter=' ', Name='T'),

    # -t  The minimum probability that a sequence must have to match a primer.
    #     (default = 6.000000e-01)
    '-t':ValuedParameter(Prefix='-', Delimiter=' ', Name='t'),
    }


    # No _get_results_path needed as all results (the merged paired-ends)
    # are sent to STDOUT.

    def getHelp(self):
        """pandaseq help"""
        help_Str = """
        For basic help, type the following at the command line:
            'pandaseq' or 'pandaseq -h'

        Website:
            https://github.com/neufeld/pandaseq
        """


def run_pandaseq(
    reads1_infile_path,
    reads2_infile_path,
    outfile_label='pandaseq',
    phred_64=False,
    fastq=True,
    threads=1,
    params={},
    working_dir=tempfile.gettempdir(),
    SuppressStderr=True,
    SuppressStdout=False,
    HALT_EXEC=False):
    """ Runs pandaseq with default parameters to assemble paired-end reads.
        -reads1_infile_path : reads1.fastq infile path
        -reads2_infile_path : reads2.fastq infile path
        -outfile_label : base name for outfile
        -fastq : output assembly as fastq (True)
                 or Fasta (False)
        -threads : number of cpus to use. Default 1.
        -phred_64 : if you are using phred 64 scores instead of
                    phred 33
        -params : other optional pandaseq parameters
    """

    abs_r1_path = abspath(reads1_infile_path)
    abs_r2_path = abspath(reads2_infile_path)
    
    infile_paths = [abs_r1_path, abs_r2_path]
    # check / make absolute infile paths
    for p in infile_paths:
        if not exists(p):
            raise IOError, 'Infile not found at: %s' % p
    
    # required by pandaseq to assemble
    params['-f'] = abs_r1_path
    params['-r'] = abs_r2_path
    params['-T'] = threads
  
    # set up controller
    pandaseq_app = PandaSeq(
        params=params,
        WorkingDir=working_dir,
        SuppressStderr=SuppressStderr,
        SuppressStdout=SuppressStdout,
        HALT_EXEC=HALT_EXEC)
    
    # Fastq?
    if fastq:
        pandaseq_app.Parameters['-F'].on()

    # if using phred 64:
    if phred_64:
        pandaseq_app.Parameters['-6'].on()


    # run assembler
    result = pandaseq_app()
    
    # write STDOUT (assembly) to file 
    # NOTE: res['StdOut'] will be empty after this.
    # We do this so that the actual output remains saved to
    # disk and can be accessed outside python.
    #assembled_output_file_name = get_tmp_filename(prefix='assembled_pandaseq_')
    assembled_output_file_name = join(working_dir, outfile_label) + '_joined.fastq'
    pandaseq_outfile_handle = open(assembled_output_file_name,'w')

    for line in result['StdOut']:
        pandaseq_outfile_handle.write(line)
    pandaseq_outfile_handle.close()

    result.cleanUp()

    # We store ouput file path within a dictionary in order to kepp the 
    # expected output between various paired-end assembly wrappers
    # identical.
    path_dict = {}
    path_dict['Assembled'] = assembled_output_file_name

    # sanity check that files actually exist in path lcoations
    for path in path_dict.values():
        if not exists(path):
            raise IOError, 'Output file not found at: %s' % path

    return path_dict





