#!/usr/bin/env python
# pairs_util.py
"""Provides functions related to a Pairs object

Functions to adjust Pairs in several ways (e.g. from gapped to ungapped
or from ungapped to gapped. Works on strings or Sequence objects,
on list of tuples or Pairs objects.

The module also contains several function for measuring the distance
(or similarity) between structures. 
The metrics from Gardner and Giegerich 2004 are provided.
"""

from __future__ import division
from string import strip
from numpy import array, sqrt, searchsorted, flatnonzero, take, sum
from cogent.struct.rna2d import Pairs
from cogent.parse.fasta import MinimalFastaParser

__author__ = "Sandra Smit and Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Sandra Smit", "Shandy Wikman", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"


class PairsAdjustmentError(Exception):
    pass

# ==================================================================
# Adjustment functions for Pairs objects
# ==================================================================

def adjust_base(pairs, offset):
    """Returns new Pairs with values shifted by offset

    pairs: Pairs object or list of tuples
    offset: integer

    Adjusts the base of a pairs object or a list of pairs according to
        the given offset.
    There's no validation in here! It is possible negative values are
        returned -> user responsibility.
    This method treats all pairs as equal. It'll return a pairs object
        of exactly the same length as the input, including pairs containing
        None, and duplicates.

    Example: adjust_base(Pairs([(2,8),(4,None)]), 2) --> [(4,10),(6,None)]
    """
    if not isinstance(offset, int):
        raise PairsAdjustmentError("adjust_base: offset should be integer")
    result = Pairs()
    for x, y in pairs:
        if x is not None:
            new_x = x + offset
        else:
            new_x = x
        if y is not None:
            new_y = y + offset
        else:
            new_y = y
        result.append((new_x, new_y))
    assert len(result) == len(pairs)
    return result

def adjust_base_structures(structures, offset):
    """Adjusts the base of all structures by offset

    structures: list of Pairs objects
    offset: integer
    """
    result = []
    for struct in structures:
        result.append(adjust_base(struct, offset))
    return result

def adjust_pairs_from_mapping(pairs, mapping):
    """Returns new Pairs object with numbers adjusted according to map

    pairs: list of tuples or Pairs object
    mapping: dictionary containing mapping of positions from
        one state to the other (e.g. ungapped to gapped)
         For example:
        {0: 0, 1: 1, 2: 3, 3: 4, 4: 6, 5: 7, 6: 9, 7: 10, 8: 12}

    When the Pairs object corresponds to an ungapped sequence and
        you want to insert gaps, use a mapping from ungapped to gapped.
    When the Pairs object corresponds to a gapped sequence and you
        want to degap it, use a mapping from gapped to ungapped.
    """
    result = Pairs()
    for x,y in pairs:
        if x is None:
            new_x = None
        elif x not in mapping:
            continue
        else:
            new_x = mapping[x]
        if y is None:
            new_y = None
        elif y not in mapping:
            continue
        else:
            new_y = mapping[y]
        result.append((new_x, new_y))

    return result

def delete_gaps_from_pairs(pairs, gap_list):
    """Returns Pairs object with pairs adjusted to gap_list

    pairs: list of tuples or Pairs object
    gap_list: list or array of gapped positions that should be removed
        from the pairs object

    Base pairs of which one of the partners or both of them are in 
        the gap list are removed. If both of them are not in the gap_list, the
        numbering is adjusted according to the gap_list.
    When at least one of the two pair members is in the gap_list, the
        pair will be removed. The rest of the structure will be left
        intact. Pairs containing None, duplicates, pseudoknots, and 
        conflicts will be maintained and adjusted according to the gap_list.
    """
    if not gap_list:
        result = Pairs()
        result.extend(pairs)
        return result

    g = array(gap_list)
    result = Pairs()
    for up, down in pairs:
        if up in g or down in g:
            continue
        else:
            if up is not None:
                new_up = up - g.searchsorted(up)
            else:
                new_up = up
            if down is not None:
                new_down = down - g.searchsorted(down)
            else:
                new_down = down
            result.append((new_up, new_down))
    return result

def insert_gaps_in_pairs(pairs, gap_list):
    """Adjusts numbering in pairs according to the gap list.

    pairs: Pairs object
    gap_list: list of integers, gap positions in a sequence

    The main assumptionis that all positions in pairs correspond to
    ungapped positions. If this is not true, the result will be meaningless.
    """
    if not gap_list:
        new = Pairs()
        new.extend(pairs)
        return new

    ungapped = []
    for idx in range(max(gap_list)+2):
        if idx not in gap_list:
            ungapped.append(idx)
    new = Pairs()
    for x,y in pairs:
        if x is not None:
            try:
                new_x = ungapped[x]
            except IndexError:
                new_x = ungapped[-1] + (x-len(ungapped)+1)
        else:
            new_x = x
        if y is not None:
            try:
                new_y = ungapped[y]
            except IndexError:
                new_y = ungapped[-1] + (y-len(ungapped)+1)
        else:
            new_y = y
        new.append((new_x, new_y))
    return new


def get_gap_symbol(seq):
    """Return gap symbol.

    seq: Sequence object or plain string.

    Should be able to handle cogent.core Sequence and ModelSequence object.
    If the input sequence doesn't have a MolType, '-' will be returned
    as default.
    """
    try:
        gap = seq.MolType.Alphabet.Gap
    except AttributeError:
        gap = '-'
    return gap

def get_gap_list(gapped_seq, gap_symbol=None):
    """Return list of gapped positions.

    gapped_seq: string of sequence object. Should be able to handle
        old_cogent.base.sequence object, cogent.core Sequence and
        ModelSequence object, or plain strings.
    gap_symbol: gap symbol. Will be used for plain strings.
    """
    try: 
        gap_list = gapped_seq.gapList() #should work for RnaSequence
    except AttributeError:
        try:
            gap_list = flatnonzero(gapped_seq.gaps())
        except AttributeError:
            gap_list = flatnonzero(array(gapped_seq,'c') == gap_symbol)
    try: # if gap_list is array, convert it to list
        gap_list = gap_list.tolist()
    except AttributeError: #already a list
        pass
    return gap_list

def degap_model_seq(seq):
    """Returns ungapped copy of self, not changing alphabet.
    
    This function should actually be a method of ModelSequence. Right
        now the ungapped method is broken, so this is a temporary 
        replacement.
    """
    if seq.Alphabet.Gap is None:
        return seq.copy()
    d = take(seq._data, flatnonzero(seq.nongaps()))
    return seq.__class__(d, Alphabet=seq.Alphabet, Name=seq.Name, \
        Info=seq.Info)

def degap_seq(gapped_seq, gap_symbol=None):
    """Return ungapped copy of sequence.

    Should be able to handle
        old_cogent.base.sequence object, cogent.core Sequence and
        ModelSequence object, or plain strings.
    """
    # degap the sequence
    try: #should work for old and new RnaSequence
        ungapped_seq = gapped_seq.degap()
    except AttributeError:
        try:
            ungapped_seq = degap_model_seq(gapped_seq)
        except AttributeError:
            ungapped_symbols = take(array(list(gapped_seq)),\
                flatnonzero((array(list(gapped_seq)) != gap_symbol)))
            ungapped_seq = ''.join(ungapped_symbols)
    return ungapped_seq

def gapped_to_ungapped(gapped_seq, gapped_pairs):
    """Returns ungapped sequence and corresponding Pairs object
    
    gapped_seq: string of characters (can handle Sequence, ModelSequence,
        str, or old_cogent Sequence objects).
    gapped_pairs: Pairs object, e.g. [(3,7),(4,6)]. The Pairs object should
        correspond to the gapped sequence version.

    The gap_symbol will be extracted from the sequence object. In case
    the gapped_seq is a simple str, a '-' will be used as default.
    """
    gap_symbol = get_gap_symbol(gapped_seq)
    gap_list = get_gap_list(gapped_seq, gap_symbol)
    ungapped_seq = degap_seq(gapped_seq, gap_symbol) 
    ungapped_pairs = delete_gaps_from_pairs(gapped_pairs, gap_list)
    return ungapped_seq, ungapped_pairs

def ungapped_to_gapped(gapped_seq, ungapped_pairs):
    """Returns gapped sequence (same obj) and corresponding Pairs object
    
    gapped_seq: string of characters (can handle Sequence, ModelSequence,
        str, or old_cogent Sequence objects).
    ungapped_pairs: Pairs object, e.g. [(3,7),(4,6)]. The Pairs object should
        correspond to the ungapped sequence version.

    The gap_symbol will be extracted from the sequence object. In case
    the gapped_seq is a simple str, a '-' will be used as default.
    """
    gap_symbol = get_gap_symbol(gapped_seq)
    gap_list = get_gap_list(gapped_seq, gap_symbol)
    gapped_pairs = insert_gaps_in_pairs(ungapped_pairs, gap_list)
    return gapped_seq, gapped_pairs


# ==================================================================
# Distance/similarity measures and logical operations
# Pairs comparisons
# ==================================================================

def pairs_intersection(one, other):
    """Returns Pairs object with pairs common to one and other

    one: list of tuples or Pairs object
    other: list of tuples or Pairs object

    one and other should map onto a sequence of the same length.
    """
    pairs1 = frozenset(Pairs(one).directed()) #removes duplicates
    pairs2 = frozenset(Pairs(other).directed())
    return Pairs(pairs1&pairs2)

def pairs_union(one, other):
    """Returns the intersection of one and other

    one: list of tuples or Pairs object
    other: list of tuples or Pairs object

    one and other should map onto a sequence of the same length.
    """
    pairs1 = frozenset(Pairs(one).directed()) #removes duplicates
    pairs2 = frozenset(Pairs(other).directed())
    return Pairs(pairs1 | pairs2)

def compare_pairs(one, other):
    """Returns size of intersection divided by size of union between two Pairs

    Use as a similiraty measure for comparing secondary structures.
    Returns the number of base pairs common to both structures divided by
    the number of base pairs that is in one or the other structure:
    (A AND B)/(A OR B) (intersection/union)

    one: list of tuples or Pairs object
    other: list of tuples or Pairs object
    """
    if one.hasConflicts() or other.hasConflicts():
        raise ValueError("Can't handle conflicts in the structure""")
    if not one and not other:
        return 1.0
    pairs1 = frozenset(Pairs(one).directed()) #removes duplicates
    pairs2 = frozenset(Pairs(other).directed())
    return len(pairs1 & pairs2)/len(pairs1|pairs2)

def compare_random_to_correct(one, other):
    """Returns fraction of bp in one that is in other (correct)
    
    one: list of tuples or Pairs object
    other: list of tuples or Pairs object

    Note: the second structure is the one compared against (the correct
        structure)
    """
    if not one and not other:
        return 1.0
    if not one or not other:
        return 0.0
    pairs1 = frozenset(Pairs(one).directed()) #removes duplicates
    pairs2 = frozenset(Pairs(other).directed())
    return len(pairs1 & pairs2)/len(pairs1)

def compare_pairs_mapping(one, other, one_to_other):
    """Returns intersection/union given a mapping from the first pairs to second

    Use in case the numbering of the two Pairs object don't correspond.
    Sort of aligning two ungapped sequences and comparing their Pairs
    object via a mapping. 
    
    one: list of tuples or Pairs object
    other: list of tuples or Pairs object
    one_to_other: mapping of positions in first pairs object to positions
        in second pairs object.
    For example:
    # pos in first seq, base, pos in second seq
    #1 U 0
    #2 C 1
    #3 G 2
    #4 A 3
    #  A 4
    #5 C 5
    #6 C 6
    #7 U 
    #8 G 7

    mapping = {1:0, 2:1, 3:2, 4:3, 5:5, 6:6, 7:None, 8:7}
    """
    if not one and not other:
        return 1.0
    just_in_first = 0
    just_in_second = 0
    in_both = 0
    pairs1 = Pairs(one).directed() #removes duplicates
    pairs2 = Pairs(other).directed()
    for x,y in pairs1:
        other_match = (one_to_other[x],one_to_other[y])
        if other_match in pairs2:
            in_both += 1
            pairs2.remove(other_match)
        else:
            just_in_first += 1
    just_in_second += len(pairs2)

    return in_both/(just_in_first + in_both + just_in_second)

# ===========================================================
# Gardner & Giegerich 2004 metrics
# ===========================================================
ACCEPTED = dict.fromkeys(map(tuple,["GC","CG","AU","UA","GU","UG"]))

def check_structures(ref, predicted):
    """Raise ValueError if one of the two structures contains conflicts"""
    if ref.hasConflicts():
        raise ValueError("Reference structure contains conflicts")
    if predicted.hasConflicts():
        raise ValueError("Predicted structure contains conflicts")

def get_all_pairs(sequences, min_dist=4):
    """Return number of possible base pairs in the sequece
    
    sequences: list of Sequence objects or strings
    min_dist: integer, minimum distance between two members of a 
        base pair. Default is 4 (i.e. minimum of 3 unpaired bases in a loop)

    The number of pairs is defined as all possible GC, AU, and GU pairs, 
        respecting the minimum distance between the two members of a
        base pair.
    This method returns the average number of possible base pairs over
        all provided sequences.
    """
    if min_dist < 1:
        raise ValueError("Minimum distance should be >= 1")
    if not sequences:
        return 0.0
    tn_counts = []
    for seq in sequences:
        seq_str = str(seq).upper()
        seq_count = 0
        #print 'xrange', range(len(seq)-min_dist)
        for x in range(len(seq)-min_dist):
            for y in range(x+min_dist,len(seq)):
                if (seq_str[x],seq_str[y]) in ACCEPTED:
                    #print x,y, seq_str[x], seq_str[y], 'Y'
                    seq_count += 1
                else:
                    pass
                    #print x,y, seq_str[x], seq_str[y], 'N'
        tn_counts.append(seq_count)
    return sum(tn_counts)/len(tn_counts)

def get_counts(ref, predicted, split_fp=False, sequences=None, min_dist=4):
    """Return TP, TN, FPcont, FPconf FPcomp, FN counts"""

    result = dict.fromkeys(['TP','TN','FN','FP',\
        'FP_INCONS','FP_CONTRA','FP_COMP'],0)

    ref_set = frozenset(Pairs(ref).directed())
    pred_set = frozenset(Pairs(predicted).directed())

    ref_dict = dict(ref.symmetric())
    pred_dict = dict(predicted.symmetric())

    tp_pairs = ref_set.intersection(pred_set)
    fn_pairs = ref_set.difference(pred_set)
    fp_pairs = pred_set.difference(ref_set)
    result['TP'] = len(tp_pairs)
    result['FN'] = len(fn_pairs)
    result['FP'] = len(fp_pairs)
    if split_fp:
        fp_incons = []
        fp_contra = []
        fp_comp = []
        for x,y in fp_pairs:
            if x in ref_dict or y in ref_dict:
                #print "Conflicting: %d - %d"%(x,y)
                fp_incons.append((x,y))
            else:
                five_prime = x
                three_prime = y
                contr_found = False
                for idx in range(x,y+1):
                    if idx in ref_dict and\
                        (ref_dict[idx] < five_prime or\
                            ref_dict[idx] > three_prime):
                        #print "Contradicting: %d - %d"%(x,y)
                        contr_found = True
                        fp_contra.append((x,y))
                        break
                if not contr_found:
                    #print "Comatible: %d - %d"%(x,y)
                    fp_comp.append((x,y))


        result['FP_INCONS'] = len(fp_incons)
        result['FP_CONTRA'] = len(fp_contra)
        result['FP_COMP'] = len(fp_comp)
        assert result['FP_INCONS'] + result['FP_CONTRA'] + result['FP_COMP'] ==\
            result['FP']
    if sequences:
        num_possible_pairs = get_all_pairs(sequences, min_dist)
        result['TN'] = num_possible_pairs - result['TP'] -\
            result['FP_INCONS'] - result['FP_CONTRA']

    return result

def extract_seqs(seqs):
    """Return list of sequences as strings.

    seqs could either be:
    -- a long string in fasta format:
    ">seq1\nACGUAGC\n>seq2\nGGUAGCG"
    -- a list of lines in fasta format:
    [">seq1","ACGUAGC",">seq2","GGUAGCG"]
    -- a list of sequences (strings or objects):
    ['ACGUAGC','GGUAGCG']
    """
    if isinstance(seqs, str): #assume fasta string
        result = [v for (l,v) in list(MinimalFastaParser(seqs.split('\n')))]
    elif isinstance(seqs, list):
        seq_strings = map(strip,map(str, seqs))
        if seq_strings[0].startswith('>'): #list of fasta lines
            result = [v for l,v in list(MinimalFastaParser(seq_strings))]
        else:
            result = seq_strings
    else:
        raise Exception
    result = [s.replace('T','U') for s in result]
    return result

def sensitivity_formula(counts):
    """Return sensitivity
    
    counts: dict of counts, containing at least TP and FN
    """
    tp = counts['TP']
    fn = counts['FN']

    if not tp and not fn:
        return 0.0
    
    sensitivity = tp/(tp + fn)
    return sensitivity

def selectivity_formula(counts):
    """Return selectivity
    
    counts: dict of counts, containing at least TP, FP, and FP_COMP
    """
    tp = counts['TP']
    fp = counts['FP']
    fp_comp = counts['FP_COMP']
    
    if not tp and fp==fp_comp:
        return 0.0

    selectivity = tp/(tp + (fp - fp_comp))
    return selectivity

def ac_formula(counts):
    """Return approximate correlation
    
    counts: dict of counts, containing at least TP, FP, and FP_COMP
    """
    sens = sensitivity_formula(counts)
    sel = selectivity_formula(counts)
    return (sens+sel)/2

def cc_formula(counts):
    """Return correlation coefficient
    
    counts: dict of counts, containing at least TP, TN, FN, FP, and FP_COMP
    """
    tp = counts['TP']
    tn = counts['TN']
    fp = counts['FP']
    fn = counts['FN']
    comp = counts['FP_COMP']
    
    sens = sensitivity_formula(counts)
    sel = selectivity_formula(counts)

    N = tp+ (fp-comp) + fn + tn
    cc = 0.0
    if tp >0:
        cc = (N*sens*sel-tp)/sqrt((N*sens-tp)*(N*sel-tp))
    return cc

def mcc_formula(counts):
    """Return correlation coefficient
    
    counts: dict of counts, containing at least TP, TN, FN, FP, and FP_COMP
    """
    tp = counts['TP']
    tn = counts['TN']
    fp = counts['FP']
    fn = counts['FN']
    comp = counts['FP_COMP']
    
    mcc_quotient = (tp+fp-comp)*(tp+fn)*(tn+fp-comp)*(tn+fn)
    if mcc_quotient > 0:
        mcc = (tp*tn-(fp-comp)*fn)/sqrt(mcc_quotient)
    else:
        raise ValueError("mcc_quotient <= 0: %.2f"%(mcc_quotient))

    return mcc

def sensitivity(ref, predicted):
    """Return sensitivity of the predicted structure

    ref: Pairs object -> reference structure (true structure)
    predicted: Pairs object -> predicted structure

    Formula: sensitivity = tp/(tp + fn)
    tp = True positives
    fn = False negatives
    """
    check_structures(ref, predicted)
    if not ref and not predicted:
        return 1.0
    elif not predicted:
        return 0.0

    counts = get_counts(ref, predicted)
    return sensitivity_formula(counts)

def selectivity(ref,predicted):
    """Return selectivity of the predicted structure
    
    ref: Pairs object -> reference structure (true structure)
    predicted: Pairs object -> predicted structure

    Formula: selectivity = tp/(tp+fp-fp_comp)
    tp = True positives
    fp = False positives
    fp_comp = compatible fp pairs
    """
    check_structures(ref, predicted)
    if not ref and not predicted:
        return 1.0
    elif not predicted:
        return 0.0

    counts = get_counts(ref, predicted, split_fp=True)
    return selectivity_formula(counts)
    
def selectivity_simple(ref, predicted):
    """Return selectivity without subtracting compatible false positives
    
    ref: Pairs object -> reference structure (true structure)
    predicted: Pairs object -> predicted structure

    Formula: selectivity = tp/(tp+fp)
    tp = True positives
    fp = False positives

    Not considering compatible false positives.
    As implemented in Dowell 2004
    """
    check_structures(ref, predicted)
    if not ref and not predicted:
        return 1.0
    elif not predicted:
        return 0.0

    counts = get_counts(ref, predicted)
    tp = counts['TP']
    fp = counts['FP']
    
    if not tp:  #and fp==fp_comp:
        return 0.0

    selectivity = tp/(tp + fp)
    return selectivity

def approximate_correlation(ref, predicted, seqs):
    """Return the approximate correlation between sensitivity and selectivity

    ref: Pairs object -> reference structure (true structure)
    predicted: Pairs object -> predicted structure

    For the specific case of RNA structure comparisons, Matthews
        correlation coefficient can be approximated by the arithmetic-mean
        or geometrix-mean of sensitivity and selectivity
    
    Formula: ac = (sensitivity+selectivity)/2
    """
    check_structures(ref, predicted)
    counts = get_counts(ref, predicted, split_fp=True)
    return ac_formula(counts)

def correlation_coefficient(ref, predicted, seqs, min_dist=4):
    """Return correlation coefficient to relate sensitivity and selectivity

    Implementation copied from compare_ct.pm
    Always same result as MCC?
    """
    check_structures(ref, predicted)
    sequences = extract_seqs(seqs)
    counts = get_counts(ref, predicted, sequences=sequences, split_fp=True,\
        min_dist=min_dist)
    return cc_formula(counts)
    
def mcc(ref, predicted, seqs, min_dist=4):
    """Return the Matthews correlation coefficient

    ref: Pairs object -> reference structure (true structure)
    predicted: Pairs object -> predicted structure
    seqs: list of sequences, necessary to compute the number of true
        negatives. See documentation of extract_seqs function for 
        accepted formats.
    min_dist: minimum distance required between two members of a base pair.
        Needed to calculate the number of true negatives.
    """
    check_structures(ref, predicted)
    if not ref and not predicted:
        return 1.0
    elif not predicted:
        return 0.0
    elif not seqs:
        raise ValueError, 'No sequence provided!'

    sequences = extract_seqs(seqs)
    counts = get_counts(ref, predicted, sequences=sequences, split_fp=True,\
        min_dist=min_dist)
    return mcc_formula(counts)

def all_metrics(ref, predicted, seqs, min_dist=4):
    """Return dictionary containing the values of five metrics

    ref: Pairs object -> reference structure (true structure)
    predicted: Pairs object -> predicted structure
    seqs: list of sequences, necessary to compute the number of true
        negatives. See documentation of extract_seqs function for 
        accepted formats.
    min_dist: minimum distance required between two members of a base pair.
        Needed to calculate the number of true negatives.

    the metrics returned are:
        sensitivity, selectivity, approximate correlation,
        correlation coefficient, and Matthews correlation coefficient
    """
    check_structures(ref, predicted)
    
    result = {}
    if not ref and not predicted: # set all to 1.0
        for i in ['SENSITIVITY','SELECTIVITY','AC','CC','MCC']:
            result[i] = 1.0
        return result
    elif not predicted: # set all to 0.0
        for i in ['SENSITIVITY','SELECTIVITY','AC','CC','MCC']:
            result[i] = 0.0
        return result
    elif not seqs:
        raise ValueError, 'No sequence provided!'

    sequences = extract_seqs(seqs)
    counts = get_counts(ref, predicted, sequences=sequences, split_fp=True,\
        min_dist=min_dist)

    result['SENSITIVITY'] = sensitivity_formula(counts)
    result['SELECTIVITY'] = selectivity_formula(counts)
    result['AC'] = ac_formula(counts)
    result['CC'] = cc_formula(counts)
    result['MCC'] = mcc_formula(counts)
    return result


if __name__ == "__main__":
    pass
