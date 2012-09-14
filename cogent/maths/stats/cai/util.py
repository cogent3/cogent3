#!/usr/bin/env python
""" Provides codon usage bias measurements, the Codon Adaptation Index (CAI).

Codon bias patterns for a gene have been shown to coorelate with translation
levels and the accuracy of translation.  A CAI is derived from codon preference
statistics.  This CAI evaluation measures the extent to which codon usage in a
particular gene matches the usage pattern of a highly expressed gene set. Three
models have been implemented here, each as a function of either an arithmetic
or geometric mean to test multiplicative versus additive fitness:
1) Assumes each gene selects codons and amino acids to maximize translation
   rate.  Gene length may also be selected.
2) Assumes each gene selects codons to maximize translation rate, but not amino
   amino acids.  Assumes the absolute translation rate, independent of amino
   acid frequency, is what's important (Sharp and Li, 1987).
3) Same as (2) above, but assumes the per codon translation rate is what's
   important (Eyre-Walker, 1996).

Note:
-   Frequencies for single codon families (Met and Trp) are omitted from the
    calculations since they can not exhibit bias.
-   Stop codon frequencies are not omitted.
-   The Bulmer (1988) correction is used to prevent small normalized
    frequencies from creating a bias.
"""
from __future__ import division
from numpy import log , exp
from cogent.core.genetic_code import GeneticCodes

__author__ = "Michael Eaton"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Michael Eaton", "Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

# Default cutoff for small codon frequencies to avoid bias
FREQUENCY_MIN = 1.0e-3
# Default Standard Genetic Code NCBI number
SGC = 1
# Default set of codons -- note: need to be RNA for compatibility with
#the standard CodonUsage object.
bases = 'UGAC'
cu = dict.fromkeys([i+j+k for i in bases for j in bases for k in bases], 0.0)

def as_rna(s): 
    """Converts string to uppercase RNA"""
    return s.upper().replace('T','U')

def synonyms_to_rna(syn):
    """Converts values in synonyms dict to uppercase RNA"""
    result = {}
    for k, v in syn.items():
        result[k] = map(as_rna, v)
    return result

def get_synonyms(genetic_code=SGC, singles_removed=True, stops_removed=True):
    """Gets synonymous codon blocks as dict, keyed by amino acid."""
    GC = GeneticCodes[genetic_code]
    synonyms = GC.Synonyms.copy()   #don't modify original
    if stops_removed:
        if '*' in synonyms:
            del synonyms['*']
    if singles_removed:
        for aa, family in synonyms.items():
            if len(family) < 2: #delete any empty ones as well
                del synonyms[aa]
    return synonyms_to_rna(synonyms)

def sum_codon_freqs(freqs):
    """Sums a set of individual codon freqs, assuming dicts.

    Omits invalid keys.
    """
    result = cu.copy()
    for f in freqs:
        for k, v in f.items():
            if k in cu:
                cu[k] += v
    return cu

def norm_to_max(vals):
    """Normalizes items in vals relative to max val: returns copy."""
    best = max(vals)
    return [i/best for i in vals]

def arithmetic_mean(vals, freqs=None):
    """Returns arithmetic mean of vals."""
    if freqs is None:
        return sum(vals)/float(len(vals))
    else:
        return sum([v*i for v, i in zip(vals, freqs)])/sum(freqs)

def geometric_mean(vals, freqs=None):
    """Returns geometric mean of vals."""
    if freqs is None:
        return exp(log(vals).sum()/float(len(vals)))
    else:
        return exp(sum([v*i for v,i in zip(log(vals), freqs)])/sum(freqs))

def codon_adaptiveness_all(freqs):
    """Calculates relative codon adaptiveness, using all codons."""
    k = freqs.keys()
    v = freqs.values()
    n = norm_to_max(freqs.values())
    return dict(zip(freqs.keys(), norm_to_max(freqs.values())))

def codon_adaptiveness_blocks(freqs, blocks):
    """Calculates relative codon adaptiveness, using codon blocks."""
    result = freqs.copy()
    for b, codons in blocks.items():
        codon_vals = norm_to_max([result[c] for c in codons])
        for c, v in zip(codons, codon_vals):
            result[c] = v
    return result

def set_min(freqs, threshold):
    """Sets all values in freqs below min to specified threshold, in-place."""
    for k, v in freqs.items():
        if v < threshold:
            freqs[k] = threshold

def valid_codons(blocks):
    """Gets all valid codons from blocks"""
    result = []
    for b in blocks.values():
        result.extend(b)
    return result

def make_cai_1(genetic_code=SGC):
    """Returns function that calculates CAI model 1 using specified gen code."""
    blocks = get_synonyms(genetic_code)
    codons = frozenset(valid_codons(blocks))

    def cai_1(ref_freqs, gene_freqs, average=arithmetic_mean, \
        threshold=FREQUENCY_MIN):
        """ Assumes codon and amino acid selection to maximize translation rate.

        Gene length may also be under selection.

        ref_freqs: dict mapping reference set of codons for
            highly expressed genes to their frequencies of occurrence.
        gene_freqs: dictionary mapping codons to their usage in gene of
            interest.
        average: function for normalizing and averaging gene of interest codon
            frequencies.
        threshold: cutoff for small normalized codon frequencies to avoid bias.
        """
        r = ref_freqs.copy()
        set_min(r, threshold)
        adaptiveness_values = codon_adaptiveness_all(r)
        curr_codons = [k for k in gene_freqs if k in codons]
        return average([adaptiveness_values[i] for i in curr_codons],\
            [gene_freqs[i] for i in curr_codons])
    return cai_1

cai_1 = make_cai_1()

def make_cai_2(genetic_code=SGC):
    """Returns function that calculates CAI model 2 using specified gen code."""
    blocks = get_synonyms(genetic_code)
    codons = frozenset(valid_codons(blocks))

    def cai_2(ref_freqs, gene_freqs, average=arithmetic_mean, \
        threshold=FREQUENCY_MIN):
        """ Assumes codon, but not amino acid, selection to maximize translation
            rate (using geometric mean - Sharp and Li, 1987).

        ref_freqs: dict mapping reference set of codons for
            highly expressed genes to their frequencies of occurrence.
        gene_freqs: dictionary mapping codons to their usage in the gene of
            interest.
        average: function for normalizing and averaging the gene of interest codon
            frequencies.
        threshold: cutoff for small normalized codon frequencies to avoid bias.
        """
        r = ref_freqs.copy()
        set_min(r, threshold)
        adaptiveness_values = codon_adaptiveness_blocks(r, blocks)
        curr_codons = [k for k in gene_freqs if k in codons]
        return average([adaptiveness_values[i] for i in curr_codons],\
            [gene_freqs[i] for i in curr_codons])
    return cai_2

cai_2 = make_cai_2()

def make_cai_3(genetic_code=SGC):
    """Returns function that calculates CAI model 3 using specified gen code."""
    blocks = get_synonyms(genetic_code)

    def cai_3(ref_freqs, gene_freqs, average=arithmetic_mean, \
        threshold=FREQUENCY_MIN):
        """ Assumes codon, but not amino acid, selection to maximize translation
            rate, and per codon translation rate is paramount (using geometric 
            mean -- Eyre-Walker, 1996).

        ref_freqs: dict mapping reference set of codons for
            highly expressed genes to their frequencies of occurrence.

        gene_freqs: dictionary mapping codons to their usage in the gene of
            interest.

        average: function for normalizing and averaging the gene of interest 
            codon frequencies. if average is the string 'eyre-walker',
            will calculate CAI exactly according to the formula in his 1996 
            paper:

            CAI = exp(sum_i(sum_j(X_ij*ln_w_ij)/sum(X_ij))/m)
            ...where i is the codon block, j is the codon, X_ij is the freq,
            and w_ij is the relative adaptiveness of each codon.

            In other words, this version averages the log of the geometric
            means for each codon family, then exponentiates the result. Note
            that this actually produces exactly the same results as performing
            the geometric mean at both steps, so in general you'll just want to
            pass in the geometric mean ass the averaging function if you want to
            do this.

            Otherwise, if average is a function of the weights and frequencies,
            the same averaging function will be applied within and between
            families.

        threshold: cutoff for small normalized codon frequencies to avoid bias.
        """
        if average == 'eyre_walker':
            eyre_walker = True
            average = geometric_mean
        else:
            eyre_walker = False
        r = ref_freqs.copy()
        set_min(r, threshold)
        adaptiveness_values = codon_adaptiveness_blocks(r, blocks)
        block_results = []
        for b, codons in blocks.items():
            vals = [adaptiveness_values[i] for i in codons if i in gene_freqs]
            freqs = [gene_freqs[i] for i in codons if i in gene_freqs]
            #will skip if freqs missing
            if sum(freqs):
                block_results.append(average(vals, freqs))
        if eyre_walker:
            return exp(arithmetic_mean(log(block_results)))
        else:
            return average(block_results)
    return cai_3

cai_3 = make_cai_3()

#Define dict containing standard CAI variants
cais = {
    '1a' : lambda ref, gene, threshold=FREQUENCY_MIN: \
        cai_1(ref, gene, arithmetic_mean, threshold),
    '1g' : lambda ref, gene, threshold=FREQUENCY_MIN: \
        cai_1(ref, gene, geometric_mean, threshold),

    '2a' : lambda ref, gene, threshold=FREQUENCY_MIN: \
        cai_2(ref, gene, arithmetic_meani, threshold),
    '2g' : lambda ref, gene, threshold=FREQUENCY_MIN: \
        cai_2(ref, gene, geometric_mean, threshold),

    '3a' : lambda ref, gene, threshold=FREQUENCY_MIN: \
        cai_3(ref,gene, arithmetic_mean, threshold),
    '3g' : lambda ref, gene, threshold=FREQUENCY_MIN: \
        cai_3(ref,gene, geometric_mean, threshold),
}
