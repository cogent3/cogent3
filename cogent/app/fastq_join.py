#! /usr/bin/env python
# file: ea_utils.py

# Application controller for ea-utils v1.1.2-537 
# fastq processing utilities
# http://code.google.com/p/ea-utils/
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

class FastqJoin(CommandLineApplication):
    """fastq-join (v1.1.2) application controller for joining paired-end reads."""
    
    _command = 'fastq-join'
    
    _parameters = {
    # Description coopied from 'fastq-join'
    # Usage: fastq-join [options] <read1.fq> <read2.fq> [mate.fq] -o <read.%.fq>
    
    # Output: 
    # You can supply 3 -o arguments, for un1, un2, join files, or one 
    # argument as a file name template.  The suffix 'un1, un2, or join' is 
    # appended to the file, or they replace a %-character if present.
    # If a 'mate' input file is present (barcode read), then the files
    # 'un3' and 'join2' are also created.
    
    # -o FIL:  See 'Output' above
    '-o':ValuedParameter(Prefix='-', Delimiter=' ', Name='o'),

    # -v C:  Verifies that the 2 files probe id's match up to char C
    # use ' ' (space) for Illumina reads
    '-v':ValuedParameter(Prefix='-', Delimiter=' ', Name='v'),

    # -p N:  N-percent maximum difference (8)
    '-p':ValuedParameter(Prefix='-', Delimiter=' ', Name='p'),
    
    # -m N:  N-minimum overlap (6)
    '-m':ValuedParameter(Prefix='-', Delimiter=' ', Name='m'),
   
    # -r FIL:  Verbose stitch length report
    '-r':ValuedParameter(Prefix='-', Delimiter=' ', Name='m')}

    _input_handler = '_input_as_paths'

    def getHelp(self):
    """fastq-join (v1.1.2) help"""
    help_str =\
    """
    For issues with the actual program 'fastq-join', see the following:
    
    For basic help, type the following at the command line:
        'fastq-join'

    Website:
        http://code.google.com/p/ea-utils/

    For questions / comments subit an issue to:
    http://code.google.com/p/ea-utils/issues/list
    """
    return help_str

class FastqClipper(CommandLineApplication):
    """fastq-clipper (v1.1.2) application controller for joining paired-end reads."""
    # usage: fastq-clipper [options] <fastq-file> <adapters>

    # Removes one or more adapter sequences from the fastq file.
    # Adapter sequences are colon-delimited.
    # Stats go to stderr, unless -o is specified.

    #Options:
    # -h  This help
    # -o FIL  Output file (stats to stdout)
    # -p N    Maximum difference percentage (10)
    # -m N    Minimum clip length (1)
    # -l N    Minimum remaining sequence length (15)
    # -x [N]  Extra match length past adapter length, 
    #     N =-1 : search all
    #     N = 0 : search only up to adapter length
    # -e  End-of-line (default)
    # -b  Beginning-of-line (not supported yet)
