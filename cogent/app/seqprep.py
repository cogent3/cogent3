#! /usr/bin/env python
# file: seqprep.py

# Application controller for SeqPrep 
# https://github.com/jstjohn/SeqPrep 
# 

from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, ResultPath, \
    ApplicationError
from os.path import isabs,exists 

__author__ = "Michael Robeson"
__copyright__ = "Copyright 2007-2013, The Cogent Project"
__credits__ = ["Michael Robeson"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Michael Robeson"
__email__ = "robesonms@ornl.gov"
__status__ = "Development"

# SeqPrep help:
# Usage:
# SeqPrep [Required Args] [Options]
# NOTE 1: The output is always gziped compressed.
# NOTE 2: If the quality strings in the output contain characters less than 
# ascii 33 on an ascii table (they look like lines from a binary file), try 
# running again with or without the -6 option.
#

class SeqPrep(CommandLineApplication):
    """SeqPrep application controller for joining paired-end reads"""
    _command = 'SeqPrep'
    _parameters = {
    # Required Arguments
    # -f <first read input fastq filename>
    # -r <second read input fastq filename>
    # -1 <first read output fastq filename>
    # -2 <second read output fastq filename>
    '-f':ValuedParameter(Prefix='-', Delimiter=' ', Name='f'),
    '-r':ValuedParameter(Prefix='-', Delimiter=' ', Name='r'),
    '-1':ValuedParameter(Prefix='-', Delimiter=' ', Name='1'),
    '-2':ValuedParameter(Prefix='-', Delimiter=' ', Name='2'),
    
    # General Arguments (Optional):
    # -3 <first read discarded fastq filename>
    # -4 <second read discarded fastq filename>
    # -h Display this help message and exit (also works with no args) 
    # -6 Input sequence is in phred+64 rather than phred+33 format, the output will still be phred+33 
    # -q <Quality score cutoff for mismatches to be counted in overlap; default = 13>
    # -L <Minimum length of a trimmed or merged read to print it; default = 30>
    '-3':ValuedParameter(Prefix='-', Delimiter=' ', Name='3'),
    '-4':ValuedParameter(Prefix='-', Delimiter=' ', Name='4'),
    '-h':FlagParameter(Prefix='-', Name='h'),
    '-6':FlagParameter(Prefix='-', Name='6'),
    '-q':ValuedParameter(Prefix='-', Delimiter=' ', Name='q'),
    '-L':ValuedParameter(Prefix='-', Delimiter=' ', Name='L'),
  
    # Arguments for Adapter/Primer Trimming (Optional):
    # -A <forward read primer/adapter sequence to trim as it would appear at the 
    #   end of a read (recommend about 20bp of this)
    #	(should validate by grepping a file); 
    #   default (genomic non-multiplexed adapter1) = AGATCGGAAGAGCGGTTCAG>
    # -B <reverse read primer/adapter sequence to trim as it would appear at the 
    #   end of a read (recommend about 20bp of this)
    #	(should validate by grepping a file); 
    #   default (genomic non-multiplexed adapter2) = AGATCGGAAGAGCGTCGTGT>
    # -O <minimum overall base pair overlap with adapter sequence to trim; 
    #   default = 10>
    # -M <maximum fraction of good quality mismatching bases for primer/adapter
    #    overlap; default = 0.020000>
    # -N <minimum fraction of matching bases for primer/adapter overlap; 
    #   default = 0.870000>
    # -b <adapter alignment band-width; default = 50>
    # -Q <adapter alignment gap-open; default = 8>
    # -t <adapter alignment gap-extension; default = 2>
    # -e <adapter alignment gap-end; default = 2>
    # -Z <adapter alignment minimum local alignment score cutoff 
    #   [roughly (2*num_hits) - (num_gaps*gap_open) - (num_gaps*gap_close) - 
    #   (gap_len*gap_extend) - (2*num_mismatches)]; default = 26>
    # -w <read alignment band-width; default = 50>
    # -W <read alignment gap-open; default = 26>
    # -p <read alignment gap-extension; default = 9>
    # -P <read alignment gap-end; default = 5>
    # -X <read alignment maximum fraction gap cutoff; default = 0.125000>
    '-A':ValuedParameter(Prefix='-', Delimiter=' ', Name='A'),
    '-B':ValuedParameter(Prefix='-', Delimiter=' ', Name='B'),
    '-O':ValuedParameter(Prefix='-', Delimiter=' ', Name='O'),
    '-M':ValuedParameter(Prefix='-', Delimiter=' ', Name='M'),
    '-N':ValuedParameter(Prefix='-', Delimiter=' ', Name='N'),
    '-b':ValuedParameter(Prefix='-', Delimiter=' ', Name='b'),
    '-Q':ValuedParameter(Prefix='-', Delimiter=' ', Name='Q'),
    '-t':ValuedParameter(Prefix='-', Delimiter=' ', Name='t'),
    '-e':ValuedParameter(Prefix='-', Delimiter=' ', Name='e'),
    '-Z':ValuedParameter(Prefix='-', Delimiter=' ', Name='Z'),
    '-w':ValuedParameter(Prefix='-', Delimiter=' ', Name='w'),
    '-W':ValuedParameter(Prefix='-', Delimiter=' ', Name='W'),
    '-p':ValuedParameter(Prefix='-', Delimiter=' ', Name='p'),
    '-P':ValuedParameter(Prefix='-', Delimiter=' ', Name='P'),
    '-X':ValuedParameter(Prefix='-', Delimiter=' ', Name='X'),

    # Optional Arguments for Merging:
    # -y <maximum quality score in output ((phred 33) default = ']' )>
    # -g <print overhang when adapters are present and stripped (use this if 
    #   reads are different length)>
    # -s <perform merging and output the merged reads to this file>
    # -E <write pretty alignments to this file for visual Examination>
    # -x <max number of pretty alignments to write (if -E provided);
    #   default = 10000>
    # -o <minimum overall base pair overlap to merge two reads; default = 15>
    # -m <maximum fraction of good quality mismatching bases to overlap reads;
    #   default = 0.020000>
    # -n <minimum fraction of matching bases to overlap reads;
    #   default = 0.900000>
    '-y':ValuedParameter(Prefix='-', Delimiter=' ', Name='y'),
    '-g':FlagParameter(Prefix='-', Name='y'),
    '-s':ValuedParameter(Prefix='-', Delimiter=' ', Name='s'),
    '-E':ValuedParameter(Prefix='-', Delimiter=' ', Name='E'),
    '-x':ValuedParameter(Prefix='-', Delimiter=' ', Name='x'),
    '-o':ValuedParameter(Prefix='-', Delimiter=' ', Name='o'),
    '-m':ValuedParameter(Prefix='-', Delimiter=' ', Name='m'),
    '-n':ValuedParameter(Prefix='-', Delimiter=' ', Name='n')}


    def _get_result_paths(self, data):
        """Captures SeqPrep output.
        
        """
        result = {}
        
        # required for assembly
        result['Reads1Out'] = ResultPath(Path = , IsWritten=True)
        result['Reads2Out'] = ResultPath(Path = , IsWritten=True)
        result['Assembled'] = ResultPath(Path = , IsWritten=True)

        # optional
        result['Reads1Discarded'] = ResultPath(Path = , IsWritten=True)
        result['Reads2Discarded'] = ResultPath(Path = , IsWritten=True)
        result['PrettyAlignments'] ResultPath(Path = , IsWritten=True)
        
        return result


    def getHelp(self):
        """seqprep help"""
        help_str =\
        """
        For basic help, type the following at the command line:
            'SeqPrep -h'

        Website:
            https://github.com/jstjohn/SeqPrep
        """
        return help_str









    
