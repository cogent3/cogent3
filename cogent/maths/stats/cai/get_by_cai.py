#!/usr/bin/env python
"""Produces data that matches the CAI for a gene against its P3."""
from cogent.parse.cutg import CutgParser
from cogent.parse.fasta import MinimalFastaParser
from cogent.core.usage import UnsafeCodonUsage as CodonUsage, \
    UnsafeCodonsFromString
from cogent.maths.stats.cai.util import cais
from cogent.maths.stats.cai.adaptor import consolidate, read_nt

__author__ = "Stephanie Wilson"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Stephanie Wilson"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class CaiFilter(object):
    """Returns filter that checks objects for CAI range. Abstract."""
    def __init__(self, training_set, min_cai, max_cai=1, cai_type='2g'):
        self._cai_type = cai_type
        self._training_set = training_set
        self._min_cai = min_cai
        self._max_cai = max_cai
        
    def __call__(self, *args, **kwargs):
        raise NotImplemented

class CaiSeqFilter(CaiFilter):
    """Returns filter that checks seqs for CAI range."""
    def __call__(self, seq):
        """Returns True if within CAI threshold."""
        u = UnsafeCodonsFromString(seq.upper().replace('T','U'))
        return self._min_cai < cais[self._cai_type]([self._training_set], u) \
            < self._max_cai

class CaiUsageFilter(CaiFilter):
    """Returns filter that checks usages for CAI range."""
    def __call__(self, u):
        """Returns True if within CAI threshold."""
        return self._min_cai < cais[self._cai_type]([self._training_set], u) \
            < self._max_cai

def filter_seqs(infile, f):
    """Iterates (seq, label) from Fasta-format infile"""
    for label, seq in MinimalFastaParser(infile):
        if f(seq):
            yield label, seq

#run from command-line
if __name__ == '__main__':
    from sys import argv
    infile = open(argv[1])
    training = read_nt(open(argv[2]))
    training_freqs = consolidate(training)
    min_cai = float(argv[3])
    try:
        max_cai = float(argv[4])
    except:
        max_cai = 1
    f = CaiSeqFilter(training_freqs, min_cai, max_cai)
    for label, seq in filter_seqs(infile, f):
        print '>'+label+'\n'+seq+'\n'
