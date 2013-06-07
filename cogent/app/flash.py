#! /usr/bin/env python
# file: flash.py

# Application controller for FLASh v1.2.6 
# Fast Length Adjustment of Short reads:
# http://ccb.jhu.edu/software/FLASH/

from cogent.app.parameters import ValuedParameter, FlagParameter, \
    MixedParameter
from cogent.app.util import CommandLineApplication, ResultPath, \
    ApplicationError


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
        # and pasted below.

        # --min-overlap
        # The minimum required overlap length between two
        # reads to provide a confident overlap.  Default:
        # 10bp. 
        '-m':ValuedParameter('-', Delimiter=' ', Name='m', Value='10'),
        
        # --max-overlap
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
        '-M':ValuedParameter('-', Delimiter=' ', Name='M', Value='70'),
                 }

        # --max-mismatch-density
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
        '-x':ValuedParameter('-', Delimiter=' ', Name='x', Value='0.25'),

        # --phred-offset=OFFSET
        # The smallest ASCII value of the characters used to
        # represent quality values of bases in FASTQ files.
        # It should be set to either 33, which corresponds
        # to the later Illumina platforms and Sanger
        # platforms, or 64, which corresponds to the
        # earlier Illumina platforms.  Default: 33.
        '-p':ValuedParameter('-',Delimiter=' ', Delimiter, Name='p', Value='33')

 


    #TODO
    # make convenience functions with parameters for typical HISEQ vs MISEQ
    # that is default FLASh == HISEQ. For MISEQ use something similar to:
    # -r 250 -f 340 -s 34 -M 500   
