#/usr/bin/env python
"""Adaptors to fit data read from CUTG or fasta-format files into CAI graphs.
"""
from cogent.core.usage import UnsafeCodonUsage as CodonUsage
from cogent.core.info import Info
from cogent.parse.cutg import CutgParser
from cogent.parse.fasta import MinimalFastaParser
from cogent.core.genetic_code import GeneticCode
from numpy import array, arange, searchsorted, sort, sqrt, zeros, \
    concatenate, transpose
from sys import argv
from string import split
from cogent.maths.stats.cai.util import cais

__author__ = "Stephanie Wilson"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Stephanie Wilson"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

def data_from_file(lines):
    """reads lines returns array
     """
    return array([map(float, i) for i in map(split, lines)])

empty_codons = dict.fromkeys([i+j+k for i in 'TCAG' for j in 'TCAG' for k in 'TCAG'], 0.0)

def read_cutg(lines):
    """Returns list of CUTG objects from file-like object lines.
    
    Warning: reads whole file into memory as objects.
    """
    return list(CutgParser(lines))

def get_ribosomal_proteins(usages):
    """Returns list of ribosomal proteins."""
    keys = ['rpl', 'rps']
    result=[]
    for u in usages:
        for k in keys:
            if str(u.Gene).lower().startswith(k):
                result.append(u)
                break
    return result

def consolidate(usages):
    """Sums frequencies of a list of usages into one usage."""
    result = CodonUsage()
    for u in usages:
        result += u
    result.normalize()
    return result

def make_output(training_freqs, usages, funcs=cais):
    """Makes results as table."""
    result = []
    header = funcs.keys()
    header.sort()
    result.append(['gene', 'P3'] + header)

    for u in usages:
        u.normalize()
        curr_line = []
        curr_line.append(u.Gene)
        thirdpos = u.positionalBases(purge_unwanted=True).Third
        thirdpos.normalize()
        curr_line.append(thirdpos['G'] + thirdpos['C'])
        for h in header:
            curr_line.append(funcs[h](training_freqs, u))
        result.append(curr_line)
    return result
        
def read_nt(infile):
    """Returns list of usage objects from Fasta-format infile"""
    result = []
    for label, seq in MinimalFastaParser(infile):
        u = UnsafeCodonsFromString(seq.upper().replace('T','U'))
        u.Gene = label.split()[1]
        result.append(u)
    return result

def print_output(table):
    """Prints table as tab-delimited text."""
    for line in table:
        print '\t'.join(map(str, line))


def seq_to_codon_dict(seq):
    """Converts sequence into codon dict."""
    leftover = len(seq) % 3
    if leftover:
        seq += 'A' * (3-leftover)
    result = empty_codons.copy()
    for i in range(0, len(seq), 3):
        curr = seq[i:i+3]
        if curr in result:  #ignore others
            result[curr] += 1
    return result

def kegg_fasta_to_codon_list(lines):
    """Reads list of CodonUsage objects from KEGG-format FASTA file."""
    result = []
    for label, seq in MinimalFastaParser(lines):
        seq = seq.upper()
        curr_info = {}
        fields = label.split()
        curr_info['SpeciesAbbreviation'], curr_info['GeneId'] = \
            fields[0].split(':')
        if len(fields) > 1: #additional annotation
            first_word = fields[1]
            if first_word.endswith(';'):    #gene label
                curr_info['Gene'] = first_word[:-1]
                curr_info['Description'] = ' '.join(fields[2:])
            else:
                curr_info['Description'] = ' '.join(fields[1:])
        curr_codon_usage = CodonUsage(seq_to_codon_dict(seq), Info=curr_info)
        curr_codon_usage.__dict__.update(curr_info)
        result.append(curr_codon_usage)
    return result
        
def group_codon_usages_by_ids(codon_usages, ids, field='GeneId'):
    """Sorts codon usages into 2 lists: non-matching and matching ids."""
    result = [[],[]]
    for c in codon_usages:
        result[getattr(c, field) in ids].append(c)
    return result

def file_to_codon_list(infilename):
    """converts a file from the cutg parser
    to a list of codon usages
    """
    return list(CutgParser(open(infilename), constructor=CodonUsage))

def adapt_fingerprint(codon_usages, which_blocks='quartets', \
    include_mean=True, normalize=True):
    """takes a sequence of CodonUsage objects
    and returns an array for a fingerprint plot with:
    x: the g3/(g3+c3)
    y: the a3/(a3+u3)
    frequency: total of the base/total of all

    in the order:
    alanine, arginine4, glycine, leucine4,
    proline, serine4, threonine, valine (if quartets_only is True).

    codon_usages:   list of CodonUsage objects
    quartets_only:  return only the quartets that all code for the same aa(True)
    quartets_only set to false yeilds a 16 fingerprint
    include_mean:   include a point for the mean in the result (True)
    normalize:      ensure the frequencies returned sum to 1 (True)
    """
    tot_codon_usage = CodonUsage()
    for idx, c in enumerate(codon_usages):
        tot_codon_usage += c
    return tot_codon_usage.fingerprint(which_blocks=which_blocks, \
        include_mean=include_mean, normalize=normalize)
    
def make_bin(lowerbound, upperbound, binwidth):
    """Returns range suitable for use in searchsorted(), incl. upper bound.
    
    takes: lowerbound, upperbound and bin width of
    number to be sorted
    outputs: a bin that will correspond to array indexes
    """
    return arange((lowerbound+binwidth),(upperbound-binwidth+0.0001),binwidth)

def bin_by_p3(codon_usages, bin_lowbound=0.0, bin_upbound=1.0, bin_width=0.1):
    """takes a list of one or more codon usage and bin range
    returns list of lists of codon usages, split into specified lists
    """
    p3_bin = make_bin(bin_lowbound, bin_upbound, bin_width)
    p3_list = [[] for i in p3_bin]
    p3_list.append([])
    for cu in codon_usages:
        third_usage=cu.positionalBases().Third
        #calculating the P3 value (overall)
        third_CG=(third_usage['C']+third_usage['G'])
        third_AT=(third_usage['A']+third_usage['U'])    #stored as RNA
        P3=array((third_CG/float(third_CG+third_AT)))#test if value works
        cur_p3=searchsorted(p3_bin,P3)
        p3_list[cur_p3].append(cu)
    return p3_list

def adapt_pr2_bias(codon_usages, block='GC', bin_lowbound=0.0, bin_upbound=1.0,\
    binwidth=0.1):
    """Returns the bin midpoint and the PR2 biases for each bin of GC3."""
    result = []
    for i, bin in enumerate(bin_by_p3(codon_usages, bin_lowbound, bin_upbound, \
        binwidth)):
        if not bin:
            continue
        try:
            tot_usage = CodonUsage()
            for c in bin:
                tot_usage += c
            curr_pr2 = tot_usage.pr2bias(block)
            midbin = bin_lowbound + (i+0.5)*binwidth
            result.append([midbin]+list(curr_pr2))
        except (ZeroDivisionError, FloatingPointError):
            pass
    return array(result)
        
def adapt_p12(codon_usages, purge_unwanted=True):
    """From list of codon usages, returns [P3, (P1+P2)/2] for each usage.

    purge_unwanted: get rid of singleton codons and stop codons (True).

    P3 is the resulting x axis; P12 is the resulting y axis.
    """
    data = array([c.positionalGC(purge_unwanted) for c in codon_usages])
    return [data[:,3],(data[:,1]+data[:,2])/2]

def adapt_p123gc(codon_usages, purge_unwanted=True):
    """From list of codon usages, returns [P1,P2,P3,GC] for each usage.

    purge_unwanted: get rid of singleton codons and stop codons (True).
    """
    return(transpose(array(\
        [c.positionalGC(purge_unwanted) for c in codon_usages])))

def cu_gene(obj):
    """Extracts the gene name from a codon usage object in std. format"""
    return str(obj.Gene).lower()

def bin_codon_usage_by_patterns(codon_usages, patterns, extract_gene_f=cu_gene):
    """Returns two lists of codon usages, matching vs not matching the strings.

    no_match is result[0], match is result[1].
    """
    result = [[],[]]
    for u in codon_usages:
        matched = 0
        norm_str = extract_gene_f(u)
        for p in patterns:
            if norm_str.startswith(p):
                matched = 1
                break
        result[matched].append(u)
    return result
    
def adapt_p3_histogram(codon_usages, purge_unwanted=True):
    """Returns P3 from each set of codon usage for feeding to hist()."""
    return [array([c.positionalGC(purge_unwanted=True)[3] for c in curr])\
        for curr in codon_usages]

def adapt_cai_histogram(codon_usages, cai_model='2g', patterns=['rpl','rps'],\
    purge_unwanted=True):
    """Returns two arrays of CAIs for non-training and training genes.

    Result is suitable for feeding to hist().

    codon_usages:   list of codon usages to examine.
    cai_model:      1a, 1g, 2a, 2g, 3a, 3g depending on model and 
                    arithmetic mean. See CAI.py for documentation.
    patterns:       list of patterns used to pick the training set. Default is
                    ['rpl', 'rps'] for ribosomal proteins.
    
    Returns list of results for each set. Each result for a set is an array of 
    CAI for each gene.
    """
    normal, training = bin_codon_usage_by_patterns(codon_usages, patterns)
    cai_f = cais[cai_model]
    ref = consolidate(training)
    result = []
    for series in [normal, training]:
        curr_result = []
        for c in series:
            c = c.copy()
            c.normalize()
            curr_result.append(cai_f(ref,c))
        result.append(array(curr_result))
    return result
    
def adapt_cai_p3(codon_usages,cai_model='2g', patterns=['rpl','rps'], \
    purge_unwanted=True, both_series=False):
    """Returns array of [P3,CAI] for each gene, based on training set."""
    normal, training = bin_codon_usage_by_patterns(codon_usages, patterns)
    cai_f = cais[cai_model]
    ref = consolidate(training)
    result = []
    if both_series:
        for series in [normal, training]:
            curr_result = []
            for c in series:
                c = c.copy()
                c.normalize()
                curr_result.append([c.positionalGC()[3], cai_f(ref,c)])
            result.extend(transpose(array(curr_result)))
        return result
    #if only one series was given, add all of them to the list.
    for c in codon_usages:
        c = c.copy()
        c.normalize()
        result.append([c.positionalGC()[3], cai_f(ref, c)])
    return transpose(array(result))
