""" Description
File created on 03 May 2007.

Functions to perform coevolutionary analyses on
pre-aligned biological sequences. Coevolutionary analyses detect
correlated substitutions between alignment positions. Analyses
can be performed to look for covariation between a pair of
alignment positions, in which case a single 'coevolve score' is
returned. (The nature of this coevolve score is determined by the
method used to detect coevolution.) Alternatively, coevolution can
be calculated between one position and all other positions in an
alignment, in which case a vector of coevolve scores is returned.
Finally, coevolution can be calculated over all pairs of positions
in an alignment, in which case a matrix (usually, but not necessarily,
symmetric) is returned.

The functions providing the core functionality here are:

coevolve_pair: coevolution between a pair of positions (float returned)
coevolve_position: coevolution between a position and all other
    positions in the alignment (vector returned)
coevolve_alignment: coevolution between all pairs of positions in an
    alignment (matrix returned)

Each of these functions takes a coevolution calculator, an alignment, and
 any additional keyword arguments that should be passed to the coevolution
 calculator. More information on these functions and how they should be used
 is available as executable documentation in coevolution.rst.

The methods provided for calculating coevolution are:
Mutual Information (Shannon 19xx)
Normalized Mutual Information (Martin 2005)
Statistical Coupling Analysis (Suel 2003)
*Ancestral states (Tuffery 2000 -- might not be the best ref,
 a better might be Shindyalov, Kolchannow, and Sander 1994, but so far I
 haven't been able to get my hands on that one).
* These methods require a phylogenetic tree, in addition to an alignment.
 Trees are calculated on-the-fly, by neighbor-joining, if not provided.

This file can be applied as a script to calculate a coevolution matrix given
an alignment. For information, run python coevolution.py -h from the command
line.

"""

from os.path import basename
from pickle import Pickler, Unpickler
from random import shuffle

from numpy import (
    array,
    e,
    float,
    greater_equal,
    isnan,
    less_equal,
    log,
    nan,
    nonzero,
    ones,
    ravel,
    zeros,
)
from numpy.linalg import norm

from cogent3 import PROTEIN, make_aligned_seqs
from cogent3.core.alignment import ArrayAlignment
from cogent3.core.alphabet import CharAlphabet
from cogent3.core.moltype import IUPAC_gap, IUPAC_missing
from cogent3.core.sequence import Sequence
from cogent3.evolve.substitution_model import (
    EmpiricalProteinMatrix,
    Parametric,
)
from cogent3.maths.stats.distribution import binomial_exact
from cogent3.maths.stats.number import CategoryCounter, CategoryFreqs
from cogent3.maths.stats.special import ROUND_ERROR


__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Greg Caporaso",
    "Gavin Huttley",
    "Brett Easton",
    "Sandra Smit",
    "Rob Knight",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Beta"

DEFAULT_EXCLUDES = "".join([IUPAC_gap, IUPAC_missing])
DEFAULT_NULL_VALUE = nan


def build_rate_matrix(count_matrix, freqs, aa_order="ACDEFGHIKLMNPQRSTVWY"):

    epm = EmpiricalProteinMatrix(count_matrix, freqs)
    word_probs = array([freqs[aa] for aa in aa_order])
    num = word_probs.shape[0]
    mprobs_matrix = ones((num, num), float) * word_probs

    return epm.calcQ(word_probs, mprobs_matrix)


# Mutual Information Analysis
# Mutual Information Calculators


def mi(h1, h2, joint_h):
    """Calc Mutual Information given two entropies and their joint entropy"""
    return h1 + h2 - joint_h


def normalized_mi(h1, h2, joint_h):
    """MI normalized by joint entropy, as described in Martin 2005"""
    return mi(h1, h2, joint_h) / joint_h


nmi = normalized_mi

# Other functions used in MI calculations


def join_positions(pos1, pos2):
    """Merge two positions and return as a list of strings

    pos1: iterable object containing the first positions data
    pos2: iterable object containing the second positions data

    Example:
        >>> join_positions('ABCD','1234')
            ['A1', 'B2', 'C3', 'D4']
    """
    return ["".join([r1, r2]) for r1, r2 in zip(pos1, pos2)]


def joint_entropy(pos1, pos2):
    """Calculate the joint entroy of a pair of positions"""
    return CategoryCounter(join_positions(pos1, pos2)).entropy


# Exclude handlers (functions for processing position strings with exclude
# characters)


def ignore_excludes(pos, excludes=DEFAULT_EXCLUDES):
    """Return position data as-is (results in excludes treated as other chars)"""
    return pos


# Functions for scoring coevolution on the basis of Mutual Information


def mi_pair(
    alignment,
    pos1,
    pos2,
    h1=None,
    h2=None,
    mi_calculator=mi,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
):
    """Calculate mutual information of a pair of alignment positions

    alignment: the full alignment object
    pos1: index of 1st position in alignment to be compared
     (zero-based, not one-based)
    pos2: index of 2nd position in alignment to be compared
     (zero-based, not one-based)
    h1: entropy of pos1, if already calculated (to avoid time to recalc)
    h2: entropy of pos2, if already calculated (to avoid time to recalc)
    mi_calculator: a function which calculated MI from two entropies and
     their joint entropy -- see mi and normalized_mi for examples
    null_value: the value to be returned if mi cannot be calculated (e.g.,
     if mi_calculator == normalized_mi and joint_h = 0.0)
    excludes: iterable objects containing characters that require special
     handling -- by default, if a position contains an exclude, null_value
     will be returned. For non-default handling, pass an exclude_handler
    exclude_handler: a function which takes position data and returns it
     with exclude characters processed in someway. Position data should be
     an iterable object containing the characters present at each position.
     f(position_data,excludes=gDefaultExcludes) -> position_data

    """
    col1 = list(alignment[pos1].positions)[0]
    col2 = list(alignment[pos2].positions)[0]
    # Detect and process exclude characters.
    # This bit of code is slow, and not necessary if
    # exclude_hanlder == ignore_excludes, so I explicitly
    # check, and bypass this block if possible.
    if exclude_handler != ignore_excludes:
        for col in (col1, col2):
            states = set(col)
            for exclude in excludes:
                if exclude in states:
                    try:
                        _ = exclude_handler(col, excludes)
                        break
                    except TypeError:
                        return null_value

    # Calculate entropy of pos1 & pos2, if they weren't passed in.
    if not h1:
        h1 = CategoryCounter(col1).entropy
    if not h2:
        h2 = CategoryCounter(col2).entropy
    # Calculate the joint entropy of pos1 & pos2
    joint_h = joint_entropy(col1, col2)
    # Calculate MI using the specified method -- return null_value when
    # the specified MI cannot be calculated
    # (e.g., mi_calculator=nmi and joint_h=0.0)
    try:
        result = mi_calculator(h1, h2, joint_h)
        if result <= ROUND_ERROR:
            result = 0.0
    except ZeroDivisionError:
        result = null_value
    return result


def mi_position(
    alignment,
    position,
    positional_entropies=None,
    mi_calculator=mi,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
):
    """Calc mi b/w position and all other positions in an alignment

    alignment: the full alignment object
    position: the position number of interest -- NOTE: this is the
     position index, not the sequenece position (so zero-indexed, not
    one-indexed)
    positional_entropies: a list containing the entropy of each position in
     the alignment -- these can be passed in to avoid recalculating if
     calling this function over more than one position (e.g., in
     mi_alignment)
    mi_calculator: a function which calculated MI from two entropies and
     their joint entropy -- see mi and normalized_mi for examples
    null_value: the value to be returned if mi cannot be calculated (e.g.,
     if mi_calculator == normalized_mi and joint_h = 0.0)
    excludes: iterable objects containing characters that require special
     handling -- by default, if a position contains an exclude, null_value
     will be returned. For non-default handling, pass an exclude_handler
    exclude_handler: a function which takes a position and returns it
     with exclude characters processed in someway.

    """
    aln_length = len(alignment)
    # Create result vector
    result = zeros(aln_length, float)

    # compile positional entropies if not passed in
    if positional_entropies is None:
        positional_entropies = [CategoryCounter(p).entropy for p in alignment.positions]

    # Will want to make a change here so that we don't need to recalculate
    # all values when calling from mi_alignment
    for i in range(aln_length):
        result[i] = mi_pair(
            alignment,
            pos1=position,
            pos2=i,
            h1=positional_entropies[position],
            h2=positional_entropies[i],
            mi_calculator=mi_calculator,
            null_value=null_value,
            excludes=excludes,
            exclude_handler=exclude_handler,
        )
    return result


def mi_alignment(
    alignment,
    mi_calculator=mi,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
):
    """Calc mi over all position pairs in an alignment

    alignment: the full alignment object
    mi_calculator: a function which calculated MI from two entropies and
     their joint entropy -- see mi and normalized_mi for examples
    null_value: the value to be returned if mi cannot be calculated (e.g.,
     if mi_calculator == normalized_mi and joint_h = 0.0)
    excludes: iterable objects containing characters that require special
     handling -- by default, if a position contains an exclude, null_value
     will be returned. For non-default handling, pass an exclude_handler
    exclude_handler: a function which takes a position and returns it
     with exclude characters processed in someway.

    """
    aln_length = len(alignment)
    # Create result matrix
    result = zeros((aln_length, aln_length), float)

    # Compile postional entropies for each position in the alignment
    # I believe I started using this rather than alignment.uncertainties
    # b/c the latter relies on converting a ArrayAlignment to an Alignment --
    # need to check into this.
    positional_entropies = alignment.entropy_per_pos()

    # Calculate pairwise MI between position_number and all alignment
    # positions, and return the results in a vector.
    for i in range(aln_length):
        for j in range(i + 1):
            result[i, j] = mi_pair(
                alignment,
                pos1=i,
                pos2=j,
                h1=positional_entropies[i],
                h2=positional_entropies[j],
                mi_calculator=mi_calculator,
                null_value=null_value,
                excludes=excludes,
                exclude_handler=exclude_handler,
            )
    # copy the lower triangle to the upper triangle to make
    # the matrix symmetric
    ltm_to_symmetric(result)
    return result


# End Mutual Information Analysis

# Start Normalized Mutual Information Analysis (Martin 2005)


def normalized_mi_pair(
    alignment,
    pos1,
    pos2,
    h1=None,
    h2=None,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
):
    """Calc normalized mutual information of a pair of alignment positions

    alignment: the full alignment object
    pos1: index of 1st position in alignment to be compared
     (zero-based, not one-based)
    pos2: index of 2nd position in alignment to be compared
     (zero-based, not one-based)
    h1: entropy of pos1, if already calculated (to avoid time to recalc)
    h2: entropy of pos2, if already calculated (to avoid time to recalc)
    null_value: the value to be returned if mi cannot be calculated (e.g.,
     if mi_calculator == normalized_mi and joint_h = 0.0)
    excludes: iterable objects containing characters that require special
     handling -- by default, if a position contains an exclude, null_value
     will be returned. For non-default handling, pass an exclude_handler
    exclude_handler: a function which takes a position and returns it
     with exclude characters processed in someway.

    """
    return mi_pair(
        alignment,
        pos1,
        pos2,
        h1=h1,
        h2=h2,
        mi_calculator=nmi,
        null_value=null_value,
        excludes=excludes,
        exclude_handler=exclude_handler,
    )


nmi_pair = normalized_mi_pair


def normalized_mi_position(
    alignment,
    position,
    positional_entropies=None,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
):
    """Calc normalized mi b/w position and all other positions in an alignment

    alignment: the full alignment object
    position: the position number of interest -- NOTE: this is the
     position index, not the sequenece position (so zero-indexed, not
    one-indexed)
    positional_entropies: a list containing the entropy of each position in
     the alignment -- these can be passed in to avoid recalculating if
     calling this function over more than one position (e.g., in
     mi_alignment)
    null_value: the value to be returned if mi cannot be calculated (e.g.,
     if mi_calculator == normalized_mi and joint_h = 0.0)
    excludes: iterable objects containing characters that require special
     handling -- by default, if a position contains an exclude, null_value
     will be returned. For non-default handling, pass an exclude_handler
    exclude_handler: a function which takes a position and returns it
     with exclude characters processed in someway.

    """
    return mi_position(
        alignment,
        position,
        positional_entropies=positional_entropies,
        mi_calculator=nmi,
        null_value=null_value,
        excludes=excludes,
        exclude_handler=exclude_handler,
    )


nmi_position = normalized_mi_position


def normalized_mi_alignment(
    alignment,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
):
    """Calc normalized mi over all position pairs in an alignment

    alignment: the full alignment object
    null_value: the value to be returned if mi cannot be calculated (e.g.,
     if mi_calculator == normalized_mi and joint_h = 0.0)
    excludes: iterable objects containing characters that require special
     handling -- by default, if a position contains an exclude, null_value
     will be returned. For non-default handling, pass an exclude_handler
    exclude_handler: a function which takes a position and returns it
     with exclude characters processed in someway.
    """
    return mi_alignment(
        alignment=alignment,
        mi_calculator=normalized_mi,
        null_value=null_value,
        excludes=excludes,
        exclude_handler=exclude_handler,
    )


nmi_alignment = normalized_mi_alignment
# End Normalized Mutual Information Analysis


# Start Statistical coupling analysis (SCA) (Suel 2003)
class SCAError(Exception):
    pass


# PROTEIN's alphabet contains U, so redefining the alphabet for now
# rather than use PROTEIN.alphabet. May want to revist this decision...
AAGapless = CharAlphabet("ACDEFGHIKLMNPQRSTVWY")
default_sca_alphabet = AAGapless
# AAGapless = PROTEIN.alphabet

# Dictionary of mean AA-frequencies in all natural proteins
# Compiled by Rama Ranganathan from 36,498 unique eukaryotic proteins
# from the Swiss-Prot database
protein_dict = {
    "A": 0.072658,
    "C": 0.024692,
    "D": 0.050007,
    "E": 0.061087,
    "F": 0.041774,
    "G": 0.071589,
    "H": 0.023392,
    "I": 0.052691,
    "K": 0.063923,
    "L": 0.089093,
    "M": 0.02315,
    "N": 0.042931,
    "P": 0.052228,
    "Q": 0.039871,
    "R": 0.052012,
    "S": 0.073087,
    "T": 0.055606,
    "V": 0.063321,
    "W": 0.01272,
    "Y": 0.032955,
}
default_sca_freqs = protein_dict


def freqs_to_array(f, alphabet):
    """Takes data in freqs object and turns it into array.

    f = dict or CategoryCounter object
    alphabet = Alphabet object or just a list that specifies the order
        of things to appear in the resulting array
    """
    return array([f.get(i, 0) for i in alphabet])


def get_allowed_perturbations(counts, cutoff, alphabet, num_seqs=100):
    """Returns list of allowed perturbations as characters

    count: Profile object of raw character counts at each position
    num_seqs: number of sequences in the alignment
    cutoff: minimum number of sequences in the subalignment (as fraction
    of the total number of seqs in the alignment.

    A perturbation is allowed if the subalignment of sequences that
    contain the specified char at the specified position is larger
    that the cutoff value * the total number of sequences in the alignment.

    """
    result = []
    abs_cutoff = cutoff * num_seqs

    for char, count in zip(alphabet, counts):
        if count >= abs_cutoff:
            result.append(char)
    return result


def probs_from_dict(d, alphabet):
    """Convert dict of alphabet char probabilities to list in alphabet's order

    d: probabilities of observing each character in alphabet (dict indexed
     by char)
    alphabet: the characters in the alphabet -- provided for list order.
     Must iterate over the ordered characters in the alphabet (e.g., a list
     of characters or an Alphabet object)

    """
    return array([d[c] for c in alphabet])


def freqs_from_aln(aln, alphabet, scaled_aln_size=100):
    """Return the frequencies in aln of chars in alphabet's order

    aln: the alignment object
    alphabet: the characters in the alphabet -- provided for list order.
     Must iterate over the ordered characters in the alphabet (e.g., a list
     of characters or an Alphabet object)
    scaled_aln_size: the scaled number of sequences in the alignment. The
     original SCA implementation treats all alignments as if they contained
     100 sequences when calculating frequencies and probabilities. 100 is
     therefore the default value.

    *Warning: characters in aln that are not in alphabet are silently
        ignored. Is this the desired behavior?

    Need to combine this function with get_position_frequences (and renamed
     that one to be more generic) since they're doing the same thing now.

    """
    alphabet_as_indices = array([aln.alphabet.to_indices(alphabet)]).transpose()
    aln_data = ravel(aln.array_positions)
    return (alphabet_as_indices == aln_data).sum(1) * (scaled_aln_size / len(aln_data))


def get_positional_frequencies(aln, position_number, alphabet, scaled_aln_size=100):
    """Return the freqs in aln[position_number] of chars in alphabet's order

    aln: the alignment object
    position_number: the index of the position of interest in aln
     (note: zero-based alignment indexing)
    alphabet: the characters in the alphabet -- provided for list order.
     Must iterate over the ordered characters in the alphabet (e.g., a list
     of characters or an Alphabet object)
    scaled_aln_size: the scaled number of sequences in the alignment. The
     original SCA implementation treats all alignments as if they contained
     100 sequences when calculating frequencies and probabilities. 100 is
     therefore the default value.

    *Warning: characters in aln that are not in alphabet are silently
        ignored. Is this the desired behavior?

    """
    alphabet_as_indices = array([aln.alphabet.to_indices(alphabet)]).transpose()
    position_data = aln.array_positions[position_number]
    return (alphabet_as_indices == position_data).sum(1) * (
        scaled_aln_size / len(position_data)
    )


def get_positional_probabilities(pos_freqs, natural_probs, scaled_aln_size=100):
    """Get probs of observering the freq of each char given it's natural freq
    In Suel 2003 supplementary material, this step is defined as:
     "... each element is the binomial probability of observing each
      amino acid residue at position j given its mean frequency in
      all natural proteins."
    This function performs the calculate for a single position.

    pos_freqs: the frequencies of each char in the alphabet at a
     position-of-interest in the alignment (list of floats, typically
     output of get_positional_frequencies)
    natural_probs: the natural probabilities of observing each char
     in the alphabet (list of floats: typically output of probs_from_dict)
    scaled_aln_size: the scaled number of sequences in the alignment. The
     original SCA implementation treats all alignments as if they contained
     100 sequences when calculating frequencies and probabilities. 100 is
     therefore the default value.

    Note: It is critical that the values in pos_freqs and natural_probs are
     in the same order, which should be the order of chars in the alphabet.

    """
    results = []
    for pos_freq, natural_prob in zip(pos_freqs, natural_probs):
        try:
            results.append(binomial_exact(pos_freq, scaled_aln_size, natural_prob))
        # Because of the scaling of alignments to scaled_aln_size, pos_freq is
        # a float rather than an int. So, if a position is perfectly conserved,
        # pos_freq as a float could be greater than scaled_aln_size.
        # In this case I cast it to an int. I don't like this alignment
        # scaling stuff though.
        except ValueError:
            results.append(binomial_exact(int(pos_freq), scaled_aln_size, natural_prob))
    return array(results)


def get_subalignments(aln, position, selections):
    """returns subalns w/ seq[pos] == selection for each in selections
    aln: an alignment object
    position: int in alignment to be checked for each perturbation
    selections: characters which must be present at seq[pos] for
        seq to be in subalignment

    Note: This method returns a list of subalignments corresponding
        to the list of selections. So, if you specify selections as
        ['A','G'], you would get two subalignments back -- the first
        containing sequences with 'A' at position, and the second
        containing sequences with 'G' at position. If you want all
        sequences containing either 'A' or 'G', merge the resulting
        subalignments.

    """
    result = []
    for s in aln.alphabet.to_indices(selections):
        seqs_to_keep = nonzero(aln.array_seqs[:, position] == s)[0]
        result.append(aln.get_sub_alignment(seqs=seqs_to_keep))
    return result


def get_dg(position_probs, aln_probs):
    """Return delta_g vector

    position_probs: the prob of observing each alphabet chars frequency in
     the alignment position-of-interest, given it's background frequency
     in all proteins (list of floats, typically the output of
     get_positional_probabilities)
    aln_probs: the prob of observing each alphabet chars frequency in the
     full alignment, given it's background frequency (list of floats)

    """
    results = []
    for position_prob, aln_prob in zip(position_probs, aln_probs):
        results.append(log(position_prob / aln_prob))
    return array(results)


def get_dgg(all_dgs, subaln_dgs, scaled_aln_size=100):
    """Return delta_delta_g value

    all_dgs: the dg vector for a position-of-interest in the alignment
     (list of floats, typically the output of get_dg)
    subaln_dgs: the dg vector for a sub-alignment of the position-of-
     interest in the alignment (list of floats, typically the output
     of get_dg applied to a sub-alignment)
    scaled_aln_size: the scaled number of sequences in the alignment. The
     original SCA implementation treats all alignments as if they contained
     100 sequences when calculating frequencies and probabilities. 100 is
     therefore the default value.

    * There are two weird issues in this function with respect to the
    desciption of the algorithm in the Suel 2003 supplementary material.
    In order to get the values presented in their GPCR paper, we need to
    (1) divide the euclidian norm by the scaled_aln_size, and then (2)
    multiply the result by e.
    ** IT IS CRITICAL TO UNDERSTAND WHY
    WE NEED TO APPLY THESE STEPS BEFORE PUBLISHING ANYTHING THAT USES
    THIS CODE.**

    * A possible reason for the mysterious e scaling is that we are
    misinterpreting what they mean when they say ddg is 'the magnitude of
    this difference vector.' We are assuming they are referring to the
    Euclidian norm, but until I see their code, I can't be sure about
    this.
    """
    return norm(all_dgs - subaln_dgs) / scaled_aln_size * e


def sca_pair(
    alignment,
    pos1,
    pos2,
    cutoff,
    position_freqs=None,
    position_probs=None,
    dgs=None,
    perturbations=None,
    scaled_aln_size=100,
    null_value=DEFAULT_NULL_VALUE,
    return_all=False,
    alphabet=default_sca_alphabet,
    background_freqs=default_sca_freqs,
):
    """ Calculate statistical coupling b/w a pair of alignment columns

        alignment: full alignment object
        pos1: the first position used to probe for statistical coupling
         (subalignments will be generated based on allowed perturbations
         at this position) -- int, zero-based indexing into alignment
        pos2: the second position used to probe for statistical coupling
         -- int, zero-based indexing into alignment
        cutoff: the percentage of sequences that must contain a specific
         char at a specific pos1 to result in an allowed sub-alignment.
         (According to the Ranganathan papers, this should be the value
         determined by their 3rd criteria.)
        position_freqs: if precalculated, a matrix containing the output
         of get_positional_frequencies for each position in the alignment.
         This will typically be used only when sca_pair is called from
         sca_position, and these values are therefore pre-calculated.
        position_probs: if precalculated, a matrix containing the output
         of get_positional_probabilities for each position in the alignment.
         This will typically be used only when sca_pair is called from
         sca_position, and these values are therefore pre-calculated.
        dgs: if precalculated, a matrix containing the output
         of get_dg for each position in the alignment.
         This will typically be used only when sca_pair is called from
         sca_position, and these values are therefore pre-calculated.
        perturbations: if precalculated, a matrix containing the output
         of get_allowed_perturbations for each position in the alignment.
         This will typically be used only when sca_pair is called from
         sca_position, and these values are therefore pre-calculated.
        scaled_aln_size: the scaled number of sequences in the alignment. The
         original SCA implementation treats all alignments as if they contained
         100 sequences when calculating frequencies and probabilities. 100 is
         therefore the default value.
        null_value: the value which should be returned if SCA cannot or
         should not be calculated (e.g., no allowed perturbations or
         pos1==pos2, respectively).
        return_all: if cutoff <= 0.50, it is possible that there will be more
         than one allowed_perturbation per position. In these cases, either all
         of the values could be returned (return_all=True) or the max of the
         values can be returned (return_all=False, default). If you'd like one
         value, but not the max, wrap this function with return_all=True, and
         handle the return value as desired.
        alphabet: an ordered iterable object containing the characters in the
         alphabet. For example, this can be a CharAlphabet object, a list,
         or a string.

        **IMPORTANT NOTE: SCA, unlike (all?) other methods implemented here,
         requires the full alignment, even to calculate coupling between just
         a pair of positions. Because frequencies of characters in the full
         alignment are compared with frequencies at each position, you cannot
         simply pull out two columns of the alignment, and pass them to this
         function as a subalignment. Your results would differ from calculating
         coupling of the same positions with the full alignment. For example:
            sca_pair(aln,10,20,0.85) != \
            sca_pair(aln.take_positions([10,20]),0,1,0.85)
    """

    # Calculate frequency distributions
    natural_probs = probs_from_dict(background_freqs, alphabet)
    aln_freqs = freqs_from_aln(alignment, alphabet, scaled_aln_size)
    aln_probs = get_positional_probabilities(aln_freqs, natural_probs, scaled_aln_size)

    # get positional frequencies
    if position_freqs:
        pos1_freqs = position_freqs[pos1]
        pos2_freqs = position_freqs[pos2]
    else:
        pos1_freqs = get_positional_frequencies(
            alignment, pos1, alphabet, scaled_aln_size
        )
        pos2_freqs = get_positional_frequencies(
            alignment, pos2, alphabet, scaled_aln_size
        )
    # get positional probability vectors ("... each element is the binomial
    # probability of observing each amino acid residue at position j given its
    # mean frequency in all natural proteins." Suel 2003 supplementary
    # material)
    if position_probs:
        pos2_probs = position_probs[pos2]
    else:
        pos2_probs = get_positional_probabilities(
            pos2_freqs, natural_probs, scaled_aln_size
        )

    # get statistical energies for pos2 in full alignment
    if dgs:
        pos2_dg = dgs[pos2]
    else:
        pos2_dg = get_dg(pos2_probs, aln_probs)

    # determine allowed perturbations
    if perturbations:
        allowed_perturbations = perturbations[pos1]
    else:
        allowed_perturbations = get_allowed_perturbations(
            pos1_freqs, cutoff, alphabet, scaled_aln_size
        )
    # should we do something different here on return_all == True?
    if not allowed_perturbations:
        return null_value

    # generate the subalignments which contain each allowed
    # perturbation residue at pos1
    subalignments = get_subalignments(alignment, pos1, allowed_perturbations)

    # calculate ddg for each allowed perturbation
    ddg_values = []
    for subalignment in subalignments:
        # Calculate dg for the subalignment
        subaln_freqs = freqs_from_aln(subalignment, alphabet, scaled_aln_size)
        subaln_probs = get_positional_probabilities(
            subaln_freqs, natural_probs, scaled_aln_size
        )
        subaln_pos2_freqs = get_positional_frequencies(
            subalignment, pos2, alphabet, scaled_aln_size
        )
        subaln_pos2_probs = get_positional_probabilities(
            subaln_pos2_freqs, natural_probs, scaled_aln_size
        )
        subaln_dg = get_dg(subaln_pos2_probs, subaln_probs)
        ddg_values.append(get_dgg(pos2_dg, subaln_dg, scaled_aln_size))

    if return_all:
        return list(zip(allowed_perturbations, ddg_values))
    else:
        return max(ddg_values)


def sca_position(
    alignment,
    position,
    cutoff,
    position_freqs=None,
    position_probs=None,
    dgs=None,
    perturbations=None,
    scaled_aln_size=100,
    null_value=DEFAULT_NULL_VALUE,
    return_all=False,
    alphabet=default_sca_alphabet,
    background_freqs=default_sca_freqs,
):
    """Calculate statistical coupling b/w a column and all other columns

    alignment: full alignment object
    position: the position of interest to probe for statistical coupling
     (subalignments will be generated based on allowed perturbations
     at this position) -- int, zero-based indexing into alignment
    cutoff: the percentage of sequences that must contain a specific
     char at a specific pos1 to result in an allowed sub-alignment.
     (According to the Ranganathan papers, this should be the value
     determined by their 3rd criteria.)
    position_freqs: if precalculated, a matrix containing the output
     of get_positional_frequencies for each position in the alignment.
     This will typically be used only when sca_position is called from
     sca_alignment, and these values are therefore pre-calculated.
    position_probs: if precalculated, a matrix containing the output
     of get_positional_probabilities for each position in the alignment.
     This will typically be used only when sca_position is called from
     sca_alignment, and these values are therefore pre-calculated.
    dgs: if precalculated, a matrix containing the output
     of get_dg for each position in the alignment.
     This will typically be used only when sca_position is called from
     sca_alignment, and these values are therefore pre-calculated.
    perturbations: if precalculated, a matrix containing the output
     of get_allowed_perturbations for each position in the alignment.
     This will typically be used only when sca_position is called from
     sca_alignment, and these values are therefore pre-calculated.
    scaled_aln_size: the scaled number of sequences in the alignment. The
     original SCA implementation treats all alignments as if they contained
     100 sequences when calculating frequencies and probabilities. 100 is
     therefore the default value.
    null_value: the value which should be returned if SCA cannot or
     should not be calculated (e.g., no allowed perturbations or
    pos1==pos2, respectively).
    return_all: if cutoff <= 0.50, it is possible that there will be more
     than one allowed_perturbation per position. In these cases, either all
     of the values could be returned (return_all=True) or the max of the
     values can be returned (return_all=False, default). If you'd like one
     value, but not the max, wrap this function with return_all=True, and
     handle the return value as desired.
    alphabet: an ordered iterable object containing the characters in the
     alphabet. For example, this can be a CharAlphabet object, a list,
     or a string.

    """
    natural_probs = probs_from_dict(background_freqs, alphabet)
    aln_freqs = freqs_from_aln(alignment, alphabet, scaled_aln_size)
    aln_probs = get_positional_probabilities(aln_freqs, natural_probs, scaled_aln_size)
    if not position_freqs:
        position_freqs = []
        for i in range(len(alignment)):
            position_freqs.append(
                get_positional_frequencies(alignment, i, alphabet, scaled_aln_size)
            )

    if not position_probs:
        position_probs = []
        for i in range(len(alignment)):
            position_probs.append(
                get_positional_probabilities(
                    position_freqs[i], natural_probs, scaled_aln_size
                )
            )
    if not dgs:
        dgs = []
        for i in range(len(alignment)):
            dgs.append(get_dg(position_probs[i], aln_probs))

    if not perturbations:
        perturbations = []
        for i in range(len(alignment)):
            perturbations.append(
                get_allowed_perturbations(
                    position_freqs[i], cutoff, alphabet, scaled_aln_size
                )
            )

    result = []
    for i in range(len(alignment)):
        result.append(
            sca_pair(
                alignment,
                position,
                i,
                cutoff,
                position_freqs=position_freqs,
                position_probs=position_probs,
                dgs=dgs,
                perturbations=perturbations,
                scaled_aln_size=scaled_aln_size,
                null_value=null_value,
                return_all=return_all,
                alphabet=alphabet,
                background_freqs=background_freqs,
            )
        )
    return array(result)


def sca_alignment(
    alignment,
    cutoff,
    null_value=DEFAULT_NULL_VALUE,
    scaled_aln_size=100,
    return_all=False,
    alphabet=default_sca_alphabet,
    background_freqs=default_sca_freqs,
):
    """Calculate statistical coupling b/w all columns in alignment

    alignment: full alignment object
    cutoff: the percentage of sequences that must contain a specific
     char at a specific pos1 to result in an allowed sub-alignment.
     (According to the Ranganathan papers, this should be the value
     determined by their 3rd criteria.)
    scaled_aln_size: the scaled number of sequences in the alignment. The
     original SCA implementation treats all alignments as if they contained
     100 sequences when calculating frequencies and probabilities. 100 is
     therefore the default value.
    null_value: the value which should be returned if SCA cannot or
     should not be calculated (e.g., no allowed perturbations or
    pos1==pos2, respectively).
    return_all: if cutoff <= 0.50, it is possible that there will be more
     than one allowed_perturbation per position. In these cases, either all
     of the values could be returned (return_all=True) or the max of the
     values can be returned (return_all=False, default). If you'd like one
     value, but not the max, wrap this function with return_all=True, and
     handle the return value as desired.
    alphabet: an ordered iterable object containing the characters in the
     alphabet. For example, this can be a CharAlphabet object, a list,
     or a string.

    """
    natural_probs = probs_from_dict(background_freqs, alphabet)
    aln_freqs = freqs_from_aln(alignment, alphabet, scaled_aln_size)
    aln_probs = get_positional_probabilities(aln_freqs, natural_probs, scaled_aln_size)
    # get all positional frequencies
    position_freqs = []
    for i in range(len(alignment)):
        position_freqs.append(
            get_positional_frequencies(alignment, i, alphabet, scaled_aln_size)
        )

    # get all positional probabilities
    position_probs = []
    for i in range(len(alignment)):
        position_probs.append(
            get_positional_probabilities(
                position_freqs[i], natural_probs, scaled_aln_size
            )
        )

    # get all delta_g vectors
    dgs = []
    for i in range(len(alignment)):
        dgs.append(get_dg(position_probs[i], aln_probs))

    # get all allowed perturbations
    perturbations = []
    for i in range(len(alignment)):
        perturbations.append(
            get_allowed_perturbations(
                position_freqs[i], cutoff, alphabet, scaled_aln_size
            )
        )

    result = []
    for i in range(len(alignment)):
        result.append(
            sca_position(
                alignment,
                i,
                cutoff,
                position_freqs=position_freqs,
                position_probs=position_probs,
                dgs=dgs,
                perturbations=perturbations,
                scaled_aln_size=scaled_aln_size,
                null_value=null_value,
                return_all=return_all,
                alphabet=alphabet,
                background_freqs=background_freqs,
            )
        )
    return array(result)


# End statistical coupling analysis

# Start Resampled Mutual Information Analysis
# (developed by Hutley and Easton, and first published in
# Caporaso et al., 2008)


def make_weights(counts, n):
    """Return the weights for replacement states for each possible character.
    We compute the weight as the normalized frequency of the replacement state
    divided by 2*n."""
    char_prob = list(counts.to_freqs().items())
    weights = []
    for C, P in char_prob:
        alts = CategoryFreqs({c: p for c, p in char_prob if c != C})
        alts = alts.to_normalized()
        alts = CategoryCounter({c: w / (2 * n) for c, w in list(alts.items())})
        weights += [(C, alts)]
    return weights


def calc_pair_scale(seqs, obs1, obs2, weights1, weights2):
    """Return entropies and weights for comparable alignment.
    A comparable alignment is one in which, for each paired state ij, all
    alternate observable paired symbols are created. For instance, let the
    symbols {A,C} be observed at position i and {A,C} at position j. If we
    observe the paired types {AC, AA}. A comparable alignment would involve
    replacing an AC pair with a CC pair."""
    # scale is calculated as the product of mi from col1 with alternate
    # characters. This means the number of states is changed by swapping
    # between the original and selected alternate, calculating the new mi

    pair_freqs = CategoryCounter(seqs)
    weights1 = dict(weights1)
    weights2 = dict(weights2)
    scales = []
    for a, b in list(pair_freqs.keys()):
        weights = weights1[a]

        pr = a + b
        pair_freqs -= pr
        obs1 -= a

        # make comparable alignments by mods to col 1
        for c, w in list(weights.items()):
            new_pr = c + b
            pair_freqs += new_pr
            obs1 += c

            entropy = mi(obs1.entropy, obs2.entropy, pair_freqs.entropy)
            scales += [(pr, entropy, w)]

            pair_freqs -= new_pr
            obs1 -= c

        obs1 += a
        # make comparable alignments by mods to col 2
        weights = weights2[b]
        obs2 -= b
        for c, w in list(weights.items()):
            new_pr = a + c
            pair_freqs += new_pr
            obs2 += c

            entropy = mi(obs1.entropy, obs2.entropy, pair_freqs.entropy)
            scales += [(pr, entropy, w)]

            obs2 -= c
            pair_freqs -= new_pr

        obs2 += b

        pair_freqs += pr
    return scales


def resampled_mi_pair(
    alignment,
    pos1,
    pos2,
    weights=None,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
    null_value=DEFAULT_NULL_VALUE,
):
    """returns scaled mutual information for a pair.

    Parameters
    ----------
    alignment
        Alignment instance
    pos1, pos2
        alignment positions to be assessed
    weights
        Freq objects of weights for pos1, pos2
    excludes
        states to be excluded.

    """
    positions = list(alignment.positions)
    col1 = positions[pos1]
    col2 = positions[pos2]
    seqs = ["".join(p) for p in zip(col1, col2)]
    for col in (col1, col2):
        states = {}.fromkeys(col)
        for exclude in excludes:
            if exclude in states:
                try:
                    _ = exclude_handler(col, excludes)
                    break
                except TypeError:
                    return null_value

    excludes = excludes or []
    num = len(seqs)
    col1 = CategoryCounter(col1)
    col2 = CategoryCounter(col2)
    seq_freqs = CategoryCounter(seqs)
    if weights:
        weights1, weights2 = weights
    else:
        weights1 = make_weights(col1, num)
        weights2 = make_weights(col2, num)

    entropy = mi(col1.entropy, col2.entropy, seq_freqs.entropy)
    scales = calc_pair_scale(seqs, col1, col2, weights1, weights2)
    scaled_mi = 1 - sum([w * seq_freqs[pr] for pr, e, w in scales if entropy <= e])

    return scaled_mi


rmi_pair = resampled_mi_pair


def resampled_mi_position(
    alignment,
    position,
    positional_entropies=None,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
    null_value=DEFAULT_NULL_VALUE,
):
    aln_length = len(alignment)
    result = zeros(aln_length, float)

    if positional_entropies is None:
        positional_entropies = alignment.entropy_per_pos()

    for i in range(aln_length):
        result[i] = resampled_mi_pair(
            alignment,
            pos1=position,
            pos2=i,
            excludes=excludes,
            exclude_handler=exclude_handler,
            null_value=null_value,
        )
    return result


def resampled_mi_alignment(
    alignment,
    excludes=DEFAULT_EXCLUDES,
    exclude_handler=None,
    null_value=DEFAULT_NULL_VALUE,
):
    """returns scaled mutual information for all possible pairs."""
    aln_length = len(alignment)
    result = zeros((aln_length, aln_length), float)
    positional_entropies = alignment.entropy_per_pos()

    for i in range(aln_length):
        result[i] = resampled_mi_position(
            alignment=alignment,
            position=i,
            positional_entropies=positional_entropies,
            excludes=excludes,
            exclude_handler=exclude_handler,
            null_value=null_value,
        )
    return result


# End Resampled Mutual Information Analysis

# Begin ancestral_states analysis


def get_ancestral_seqs(aln, tree, sm=None, pseudocount=1e-6, optimise=True):
    """Calculates ancestral sequences by maximum likelihood

    Parameters
    ----------
    sm
        a Parametric instance. If not provided, one is
        constructed from the alignment alphabet
    pseudocount
        unobserved sequence states must not be zero, this value
        is assigned to sequence states not observed in the alignment.
    optimise
        whether to optimise the likelihood function.

        Note: for the sake of reduced alphabets, we calculate the
         substitution model from the alignment. This also appears
         to be what what described in Tuffery 2000, although they're
         not perfectly clear about it.
    """
    sm = sm or Parametric(aln.alphabet, recode_gaps=True)
    lf = sm.make_likelihood_function(tree, sm.motif_probs)
    lf.set_alignment(aln, motif_pseudocount=pseudocount)
    if optimise:
        lf.optimise(local=True)
    return ArrayAlignment(lf.likely_ancestral_seqs(), moltype=aln.moltype)


def ancestral_state_alignment(
    aln, tree, ancestral_seqs=None, null_value=DEFAULT_NULL_VALUE
):
    ancestral_seqs = ancestral_seqs or get_ancestral_seqs(aln, tree)
    result = []
    for i in range(len(aln)):
        row = [null_value] * len(aln)
        for j in range(i + 1):
            row[j] = ancestral_state_pair(aln, tree, i, j, ancestral_seqs, null_value)
        result.append(row)
    return ltm_to_symmetric(array(result))


def ancestral_state_position(
    aln, tree, position, ancestral_seqs=None, null_value=DEFAULT_NULL_VALUE
):

    ancestral_seqs = ancestral_seqs or get_ancestral_seqs(aln, tree)
    result = []
    for i in range(len(aln)):
        result.append(
            ancestral_state_pair(aln, tree, position, i, ancestral_seqs, null_value)
        )
    return array(result)


def ancestral_state_pair(
    aln, tree, pos1, pos2, ancestral_seqs=None, null_value=DEFAULT_NULL_VALUE
):
    """ """
    ancestral_seqs = ancestral_seqs or get_ancestral_seqs(aln, tree)
    ancestral_names_to_seqs = dict(
        list(zip(ancestral_seqs.names, ancestral_seqs.array_seqs))
    )
    distances = tree.get_distances()
    tips = tree.get_node_names(tipsonly=True)
    # map names to nodes (there has to be a built-in way to do this
    # -- what is it?)
    nodes = dict([(n, tree.get_node_matching_name(n)) for n in tips])
    # add tip branch lengths as distance b/w identical tips -- this is
    # necessary for my weighting step, where we want correlated changes
    # occuring on a single branch to be given the most weight
    distances.update(dict([((n, n), nodes[n].length) for n in nodes]))
    result = 0
    names_to_seqs = dict(list(zip(aln.names, aln.array_seqs)))
    for i in range(len(tips)):
        org1 = tips[i]
        seq1 = names_to_seqs[org1]
        for j in range(i, len(tips)):
            org2 = tips[j]
            seq2 = names_to_seqs[org2]
            ancestor = nodes[org1].last_common_ancestor(nodes[org2]).name
            if ancestor == org1 == org2:
                # we're looking for correlated change along a
                # single branch
                ancestral_seq = ancestral_names_to_seqs[nodes[org1].ancestors()[0].name]
            else:
                # we're looking for correlated change along different
                # branches (most cases)
                ancestral_seq = ancestral_names_to_seqs[ancestor]

            # get state of pos1 in org1, org2, and ancestor
            org1_p1 = seq1[pos1]
            org2_p1 = seq2[pos1]
            ancestor_p1 = ancestral_seq[pos1]

            # if pos1 has changed in both organisms since their lca,
            # this is a position of interest
            if org1_p1 != ancestor_p1 and org2_p1 != ancestor_p1:
                # get state of pos2 in org1, org2, and ancestor
                org1_p2 = seq1[pos2]
                org2_p2 = seq2[pos2]
                ancestor_p2 = ancestral_seq[pos2]
                # if pos2 has also changed in both organisms since their lca,
                # then we add a count for a correlated change
                if org1_p2 != ancestor_p2 and org2_p2 != ancestor_p2:
                    # There are a variety of ways to score. The simplest is
                    # to increment by one, which seems to be what was done
                    # in other papers.) This works well, but in a quick test
                    # (alpha helices/myoglobin with several generally
                    # high scoring alphabets) weighting works better. A more
                    # detailed analysis is in order.
                    # result += 1
                    # Now I weight based on distance so
                    # changes in shorter time are scored higher than
                    # in longer time. (More ancient changes
                    # are more likely to be random than more recent changes,
                    # b/c more time has passed for the changes to occur in.)
                    # This gives results
                    # that appear to be better under some circumstances,
                    # and at worst, about the same as simply incrementing
                    # by 1.
                    result += 1 / distances[(org1, org2)]
                    # Another one to try might involve discounting the score
                    # for a pair when one changes and the other doesn't.
    return result


# End ancestral_states analysis


# Methods for running coevolutionary analyses on sequence data.
method_abbrevs_to_names = {
    "mi": "Mutual Information",
    "nmi": "Normalized Mutual Information",
    "sca": "Statistical Coupling Analysis",
    "an": "Ancestral States",
    "rmi": "Resampled Mutual Information",
}

# Method-specific error checking functions
# Some of the coevolution algorithms require method-specific input validation,
# but that code isn't included in the alrogithm-specific functions (e.g.
# sca_alignment,
# sca_pair) because those are sometimes run many times. For example,
# sca_alignment makes many calls to sca_pair, so we can't have sca_pair
# perform validation every time it's called. My solution is to have the
# coevolve_* functions perform the input validation, and recommend that
# users always perform analyses via these functions. So, in the above example,
# the user would access sca_alignment via coevolve_alignment('sca', ...). Since
# sca_alignment makes calls to sca_pair, not coevolve_pair, the input
# validation
# is only performed once by coevolve_alignment.


def sca_input_validation(alignment, **kwargs):
    """SCA specific validations steps"""

    # check that all required parameters are present in kwargs
    required_parameters = ["cutoff"]
    # users must provide background frequencies for MolTypes other
    # than PROTEIN -- by default, protein background frequencies are used.
    if alignment.moltype != PROTEIN:
        required_parameters.append("background_freqs")
    for rp in required_parameters:
        if rp not in kwargs:
            raise ValueError("Required parameter was not provided: " + rp)

    # check that the value provided for cutoff is valid (ie. between 0 and 1)
    if not 0.0 <= kwargs["cutoff"] <= 1.0:
        raise ValueError("Cutoff must be between zero and one.")

    # check that the set of chars in alphabet and background_freqs are
    # identical
    try:
        alphabet = kwargs["alphabet"]
    except KeyError:
        # We want to use the PROTEIN alphabet minus the U character for
        # proteins since we don't have a background frequency for U
        if alignment.moltype == PROTEIN:
            alphabet = AAGapless
        else:
            alphabet = alignment.moltype.alphabet
    try:
        background_freqs = kwargs["background_freqs"]
    except KeyError:
        background_freqs = default_sca_freqs
    validate_alphabet(alphabet, background_freqs)


def validate_alphabet(alphabet, freqs):
    """SCA validation: ValueError if set(alphabet) != set(freqs.keys())"""
    alphabet_chars = set(alphabet)
    freq_chars = set(freqs.keys())
    if alphabet_chars != freq_chars:
        raise ValueError(
            "Alphabet and background freqs must contain identical sets of chars."
        )


def ancestral_states_input_validation(alignment, **kwargs):
    """Ancestral States (AS) specific validations steps"""
    # check that all required parameters are present in kwargs
    required_parameters = ["tree"]
    for rp in required_parameters:
        if rp not in kwargs:
            raise ValueError("Required parameter was not provided: " + rp)

    # validate the tree
    validate_tree(alignment, kwargs["tree"])

    # if ancestral seqs are provided, validate them. (If calculated on the fly,
    # we trust them.)
    if "ancestral_seqs" in kwargs:
        validate_ancestral_seqs(alignment, kwargs["tree"], kwargs["ancestral_seqs"])


def validate_ancestral_seqs(alignment, tree, ancestral_seqs):
    """AS validation: ValueError if incompatible aln, tree, & ancestral seqs

    Incompatibility between the alignment and the ancestral_seqs is
        different sequence lengths. Incompatbility between the tree and
        the ancestral seqs is imperfect overlap between the names of the
        ancestors in the tree and the ancestral sequence names.
    """
    if len(alignment) != len(ancestral_seqs):
        raise ValueError("Alignment and ancestral seqs are different lengths.")
    # is there a better way to get all the ancestor names? why doesn't
    # tree.ancestors() do this?
    edges = set(tree.get_node_names()) - set(tree.get_tip_names())
    seqs = set(ancestral_seqs.names)
    if edges != seqs:
        raise ValueError(
            "Must be ancestral seqs for all edges and root in tree, and no more."
        )


def validate_tree(alignment, tree):
    """AS validation: ValueError if tip and seq names aren't same"""
    if set(tree.get_tip_names()) != set(alignment.names):
        raise ValueError("Tree tips and seqs must have perfectly overlapping names.")


# End method-specific error checking functions

# General (opposed to algorithm-specific) validation functions


def validate_position(alignment, position):
    """ValueError if position is outside the range of the alignment"""
    if not 0 <= position < len(alignment):
        raise ValueError(
            "Position is outside the range of the alignment: " + str(position)
        )


def validate_alignment(alignment):
    """ValueError on ambiguous alignment characters"""
    bad_seqs = []
    for name, ambiguous_pos in list(alignment.get_ambiguous_positions().items()):
        if ambiguous_pos:
            bad_seqs.append(name)
    if bad_seqs:
        raise ValueError(
            f"Ambiguous characters in sequences: {'; '.join(map(str, bad_seqs))}"
        )


def coevolve_alignments_validation(
    method, alignment1, alignment2, min_num_seqs, max_num_seqs, **kwargs
):
    """Validation steps required for intermolecular coevolution analyses"""
    valid_methods_for_different_moltypes = {}.fromkeys(
        [mi_alignment, nmi_alignment, resampled_mi_alignment]
    )
    if (
        alignment1.moltype != alignment2.moltype
    ) and method not in valid_methods_for_different_moltypes:
        raise AssertionError(
            "Different MolTypes only supported for %s"
            % " ".join(map(str, list(valid_methods_for_different_moltypes.keys())))
        )

    alignment1_names = set([n.split("+")[0].strip() for n in alignment1.names])
    alignment2_names = set([n.split("+")[0].strip() for n in alignment2.names])

    if "tree" in kwargs:
        tip_names = set(
            [n.split("+")[0].strip() for n in kwargs["tree"].get_tip_names()]
        )
        assert (
            alignment1_names == alignment2_names == tip_names
        ), "Alignment and tree sequence names must perfectly overlap"
    else:
        # no tree passed in
        assert (
            alignment1_names == alignment2_names
        ), "Alignment sequence names must perfectly overlap"

    # Determine if the alignments have enough sequences to proceed.
    if alignment1.num_seqs < min_num_seqs:
        raise ValueError(
            "Too few sequences in merged alignment: %d < %d"
            % (alignment1.num_seqs, min_num_seqs)
        )

    # Confirm that min_num_seqs <= max_num_seqs
    if max_num_seqs and min_num_seqs > max_num_seqs:
        raise ValueError(
            "min_num_seqs (%d) cannot be greater than max_num_seqs (%d)."
            % (min_num_seqs, max_num_seqs)
        )


# End general validation functions

# Start alignment-wide intramolecular coevolution analysis


# coevolve alignment functions: f(alignment,**kwargs) -> 2D array
coevolve_alignment_functions = {
    "mi": mi_alignment,
    "nmi": normalized_mi_alignment,
    "rmi": resampled_mi_alignment,
    "sca": sca_alignment,
    "an": ancestral_state_alignment,
}


def coevolve_alignment(method, alignment, **kwargs):
    """Apply coevolution method to alignment (for intramolecular coevolution)

    method: f(alignment,**kwargs) -> 2D array of coevolution scores
    alignment: alignment object for which coevolve scores should be
        calculated
    **kwargs: parameters to be passed to method()
    """
    # Perform method specific validation steps
    if method == sca_alignment:
        sca_input_validation(alignment, **kwargs)
    if method == ancestral_state_alignment:
        ancestral_states_input_validation(alignment, **kwargs)
    validate_alignment(alignment)
    return method(alignment, **kwargs)


# End alignment-wide intramolecular coevolution analysis

# Start intermolecular coevolution analysis


# Mapping between coevolve_alignment functions and coevolve_pair functions.
# These are used in coevolve_alignments, b/c under some circumstance the
# alignment function is used, and under other circumstance the pair
# function is used, but the user shouldn't have to know anything about
# that.
coevolve_alignment_to_coevolve_pair = {
    mi_alignment: mi_pair,
    normalized_mi_alignment: normalized_mi_pair,
    resampled_mi_alignment: resampled_mi_pair,
    sca_alignment: sca_pair,
    ancestral_state_alignment: ancestral_state_pair,
}


def merge_alignments(alignment1, alignment2):
    """Append alignment 2 to the end of alignment 1

    This function is used by coevolve_alignments to merge two alignments
     so they can be evaluated by coevolve_alignment.
    """
    result = {}
    # Created maps from the final seq ids (i.e., seq id before plus) to the
    # seq ids in the original alignments
    aln1_name_map = dict([(n.split("+")[0].strip(), n) for n in alignment1.names])
    aln2_name_map = dict([(n.split("+")[0].strip(), n) for n in alignment2.names])

    try:
        for merged_name, orig_name in list(aln1_name_map.items()):
            result[merged_name] = alignment1.get_gapped_seq(
                orig_name
            ) + alignment2.get_gapped_seq(aln2_name_map[merged_name])
    except ValueError:  # Differing MolTypes
        for merged_name, orig_name in list(aln1_name_map.items()):
            result[merged_name] = Sequence(
                alignment1.get_gapped_seq(orig_name)
            ) + Sequence(alignment2.get_gapped_seq(aln2_name_map[merged_name]))
    except KeyError:
        raise KeyError(
            "A sequence identifier is in alignment2 "
            + "but not alignment1 -- did you filter out sequences identifiers"
            + " not common to both alignments?"
        )
    return make_aligned_seqs(result, array_align=True)


def n_random_seqs(alignment, n):
    """Given alignment, return n random seqs in a new alignment object.

    This function is used by coevolve_alignments.

    """
    seq_names = alignment.names
    shuffle(seq_names)
    return alignment.take_seqs(seq_names[:n])


def coevolve_alignments(
    method,
    alignment1,
    alignment2,
    return_full=False,
    merged_aln_filepath=None,
    min_num_seqs=2,
    max_num_seqs=None,
    sequence_filter=n_random_seqs,
    **kwargs,
):
    """Apply method to a pair of alignments (for intermolecular coevolution)

    method: the *_alignment function to be applied
    alignment1: alignment of first molecule (ArrayAlignment)
    alignment2: alignment of second molecule (ArrayAlignment)
    return_full: if True, returns intra- and inter-molecular
     coevolution data in a square matrix (default: False)
    merged_aln_filepath: if provided, will write the merged
     alignment to file (useful for running post-processing filters)
    min_num_seqs: the minimum number of sequences that should be
     present in the merged alignment to perform the analysis
     (default: 2)
    max_num_seqs: the maximum number of sequences to include
     in an analysis - if the number of sequences exceeds
     max_num_seqs, a random selection of max_num_seqs will be
     used. This is a time-saving step as too many sequences can
     slow things down a lot. (default: None, any number of
     sequences is allowed)
    sequence_filter: function which takes an alignment and an int
     and returns the int number of sequences from the alignment in
     a new alignment object (defualt: util.n_random_seqs(alignment,n))
     if None, a ValueError will be raised if there are more than
     max_num_seqs

    This function allows for calculation of coevolve scores between
     pairs of alignments. The results are returned in a rectangular
     len(alignment1) x len(alignment2) matrix.

    There are some complications involved in preparing alignments for
     this function, because it needs to be obvious how to associate the
     putative interacting sequences. For example, if looking for
     interactions between mammalian proteins A and B, sequences are
     required from the same sets of species, and it must be apparant how
     to match the sequences that are most likely to be involved in
     biologically meaningful interactions. This typically means matching
     the sequences of proteins A&B that come from the same species. In
     other words, interaction of T. aculeatus proteinA and
     H. sapien proteinB likely don't form a biologically relevant
     interaction, because the species are so diverged.

     Matching of sequences is performed via the identifiers, but it is
     the responsibility of the user to correctly construct the sequence
     identifiers before passing the alignments (and tree, if applicable)
     to this function. To faciliate matching sequence identifiers, but not
     having to discard the important information already present in a
     sequence identifier obtained from a database such as KEGG or RefSeq,
     sequence identifiers may contain a plus symbol (+). The characters
     before the + are used to match sequences between the alignments and
     tree. The characters after the + are ignored by this function. So, a
     good strategy is to make the text before the '+' a taxonomic
     identifier and leave the text after the '+' as the original sequence
     identifier. For example, your sequence/tip names could look like:

     alignment1: 'H. sapien+gi|123', 'T. aculeatus+gi|456'
     alignment2: 'T. aculeatus+gi|999', 'H. sapien+gi|424'
     tree: 'T. aculeatus+gi|456', 'H. sapien'

     If there is no plus, the full sequence identifier will be used for the
     matching (see H. sapien in tree).  The order of sequences in the
     alignments is not important. Also note that we can't split on a colon,
     as would be convenient for pulling sequences from KEGG, because colons
     are special characters in newick.

     A WORD OF WARNING ON SEQUENCE IDENTIFIER CONSTRUCTION:
     A further complication is that in some cases, an organism will have
     multiple copies of proteins involved in a complex, but proteinA from
     locus 1 will not form a functional comples with proteinB from locus 2.
     An example of this is the three T6SSs in P. aeuroginosa. Make sure
     this is handled correctly when building your sequence identifiers!
     Sequence identifiers are used to match the sequences which are
     suspected to form a functional complex, which may not simply mean
     sequences from the same species.

    """
    # Perform general validation step
    coevolve_alignments_validation(
        method, alignment1, alignment2, min_num_seqs, max_num_seqs, **kwargs
    )
    # Append alignment 2 to the end of alignment 1 in a new alignment object
    merged_alignment = merge_alignments(alignment1, alignment2)
    validate_alignment(merged_alignment)

    if max_num_seqs and merged_alignment.num_seqs > max_num_seqs:
        try:
            merged_alignment = sequence_filter(merged_alignment, max_num_seqs)
        except TypeError:
            raise ValueError("Too many sequences for covariation analysis.")

    # If the user provided a filepath for the merged alignment, write it to
    # disk. This is sometimes useful for post-processing steps.
    if merged_aln_filepath:
        merged_aln_file = open(merged_aln_filepath, "w")
        merged_aln_file.write(merged_alignment.to_fasta())
        merged_aln_file.close()

    if return_full:
        # If the user requests the full result matrix (inter and intra
        # molecular coevolution data), call coevolve_alignment on the
        # merged alignment. Calling coevolve_alignment ensures that
        # the correct validations are performed, rather than directly
        # calling method.
        result = coevolve_alignment(method, merged_alignment, **kwargs)
        return result

    # Note: we only get here if the above if statement comes back False,
    # i.e., if we only want the intermolecular coevolution and don't care
    # about the intramolecular coevolution.

    # Get the appropriate method (need the pair method,
    # not the alignment method)
    try:
        method = coevolve_alignment_to_coevolve_pair[method]
    except KeyError:
        # may have passed in the coevolve_pair function, so just
        # continue -- will fail (loudly) soon enough if not.
        pass

    # Cache the alignment lengths b/c we use them quite a bit, and build
    # the result object to be filled in.
    len_alignment1 = len(alignment1)
    len_alignment2 = len(alignment2)
    result = array([[DEFAULT_NULL_VALUE] * len_alignment1] * len_alignment2)

    # Some of the methods run much faster if relevant data is computed once,
    # and passed in -- that is done here, but there is a lot of repeated code.
    # I'm interested in suggestions for how to make this block of code more
    # compact (e.g., can I be making better use of kwargs?).
    if method == mi_pair or method == nmi_pair or method == normalized_mi_pair:
        positional_entropies = [
            CategoryCounter(p).entropy for p in merged_alignment.positions
        ]
        for i in range(len_alignment1):
            for j in range(len_alignment2):
                result[j, i] = method(
                    merged_alignment,
                    j + len_alignment1,
                    i,
                    h1=positional_entropies[j + len_alignment1],
                    h2=positional_entropies[i],
                    **kwargs,
                )
    elif method == ancestral_state_pair:
        # Perform method-specific validations so we can safely work
        # directly with method rather than the coevolve_pair wrapper,
        # and thereby avoid validation steps on each call to method.
        ancestral_states_input_validation(merged_alignment, **kwargs)
        ancestral_seqs = get_ancestral_seqs(merged_alignment, kwargs["tree"])
        for i in range(len_alignment1):
            for j in range(len_alignment2):
                result[j, i] = method(
                    aln=merged_alignment,
                    pos1=j + len_alignment1,
                    pos2=i,
                    ancestral_seqs=ancestral_seqs,
                    **kwargs,
                )
    else:
        # Perform method-specific validations so we can safely work
        # directly with method rather than the coevolve_pair wrapper,
        # and thereby avoid validation steps on each call to method.
        if method == sca_pair:
            sca_input_validation(merged_alignment, **kwargs)
        for i in range(len_alignment1):
            for j in range(len_alignment2):
                result[j, i] = method(merged_alignment, j + len_alignment1, i, **kwargs)
    return result


# End intermolecular coevolution analysis

# Start positional coevolution analysis
# coevolve position functions: f(alignment,position,**kwargs) -> 1D array
coevolve_position_functions = {
    "mi": mi_position,
    "nmi": normalized_mi_position,
    "rmi": resampled_mi_position,
    "sca": sca_position,
    "an": ancestral_state_position,
}


def coevolve_position(method, alignment, position, **kwargs):
    """Apply provided coevolution method to a column in alignment

    method: f(alignment,position,**kwargs) -> array of coevolution scores
    alignment: alignment object for which coevolve scores should be
        calculated (ArrayAlignment)
    position: position of interest for coevolution analysis (int)
    **kwargs: parameters to be passed to method()
    """
    # Perform method-specific validation steps
    if method == sca_position:
        sca_input_validation(alignment, **kwargs)
    if method == ancestral_state_position:
        ancestral_states_input_validation(alignment, **kwargs)
    # Perform general validation steps
    validate_position(alignment, position)
    validate_alignment(alignment)
    # Perform the analysis and return the result vector
    return method(alignment, position=position, **kwargs)


# End positional coevolution analysis


# Start pairwise coevolution analysis
# coevolve pair functions: f(alignment,pos1,pos2,**kwargs) -> float
coevolve_pair_functions = {
    "mi": mi_pair,
    "nmi": normalized_mi_pair,
    "rmi": resampled_mi_pair,
    "sca": sca_pair,
    "an": ancestral_state_pair,
}


def coevolve_pair(method, alignment, pos1, pos2, **kwargs):
    """Apply provided coevolution method to columns pos1 & pos2 of alignment

    method: f(alignment,pos1,pos2,**kwargs) -> coevolution score
    alignment: alignment object for which coevolve score should be
        calculated (ArrayAlignment)
    pos1, pos2: positions to evaluate coevolution between (int)
    **kwargs: parameters to be passed to method()

    """
    # Perform method-specific validation steps
    if method == sca_pair:
        sca_input_validation(alignment, **kwargs)
    if method == ancestral_state_pair:
        ancestral_states_input_validation(alignment, **kwargs)
    # Perform general validation steps
    validate_position(alignment, pos1)
    validate_position(alignment, pos2)
    validate_alignment(alignment)
    # Perform the analysis and return the result score
    return method(alignment, pos1=pos1, pos2=pos2, **kwargs)


# End pairwise coevolution analysis
# End methods for running coevolutionary analyses on sequence data


# Coevolution matrix filters: the following functions are used as
# post-processing filters for coevolution result matrices.


def filter_threshold_based_multiple_interdependency(
    aln,
    coevolution_matrix,
    threshold=0.95,
    max_cmp_threshold=1,
    cmp_function=greater_equal,
    intermolecular_data_only=False,
):
    """Filters positions with more than max_cmp_threshold scores >= threshold

    This post-processing filter is based on the idea described in:
     "Using multiple interdependency to separate functional from
      phylogenetic correlations in protein alignments"
      Tillier and Lui, 2003

    The idea is that when a position achieved a high covariation score
     with many other positions, the covariation is more likely to arise
     from the phylogeny than from coevolution. They illustrate that this
     works in their paper, and I plan to test it with my alpha-helix-based
     analysis. Note that you can change cmp_function to change whether
     you're looking for high values to indicate covarying positions
     (cmp_function=greater_equal, used for most coevolution algorithms) or
     low values to indicate covarying positions (cmp_function=less_equal,
     used, e.g., for p-value matrices).

    aln: alignment used to generate the coevolution matrix -- this
     isn't actually used, but is required to maintain the same interface
     as other post-processing filters. Pass None if that's more convenient.
    coevolution_matrix: the 2D numpy array to be filtered. This should
     be a rectangular matrix for intermoelcular coevolution data (in which
     case intermolecular_data_only must be set to True) or a symmetric
     square matrix (when intermolecular_data_only=False)
    threshold: the threshold coevolution score that other scores should be
     compared to
    max_cmp_threshold: the max number of scores that are allowed to be
     True with respect to cmp_function and threshold (e.g., the max number
     of positions that may be greater than the threhsold) before setting
     all values associated that position to gDefaultNullValue (default: 1)
    cmp_function: the function that compares each score in
     coevolution_matrix to threshold (default: ge (greater than)) -
     function should return True if the score is one that your looking
     (e.g. score >= threshold) or False otherwise
    intermolecular_data_only: True if coevolution_matrix is a rectangular
     matrix representing an intermolecular coevolution study, and False
     if the matrix is a symmetric square matrix

    NOTE: IF intermolecular_data_only == True, coevolution_matrix MUST BE
     SYMMETRIC, NOT LOWER TRIANGULAR OR OTHERWISE NON-SYMMETRIC!!
    """
    # Determine which rows need to be filtered (but don't filter them
    # right away or subsequent counts could be off)
    filtered_rows = []
    for row_n in range(coevolution_matrix.shape[0]):
        count_cmp_threshold = 0
        for v in coevolution_matrix[row_n, :]:
            if v != DEFAULT_NULL_VALUE and cmp_function(v, threshold):
                count_cmp_threshold += 1
                if count_cmp_threshold > max_cmp_threshold:
                    filtered_rows.append(row_n)
                    break

    # if the matrix is not symmetric, determine which cols need to be filtered
    if intermolecular_data_only:
        filtered_cols = []
        for col_n in range(coevolution_matrix.shape[1]):
            count_cmp_threshold = 0
            for v in coevolution_matrix[:, col_n]:
                if v != DEFAULT_NULL_VALUE and cmp_function(v, threshold):
                    count_cmp_threshold += 1
                    if count_cmp_threshold > max_cmp_threshold:
                        filtered_cols.append(col_n)
                        break
        # filter the rows and cols in a non-symmetric matrix
        for row_n in filtered_rows:
            coevolution_matrix[row_n, :] = DEFAULT_NULL_VALUE
        for col_n in filtered_cols:
            coevolution_matrix[:, col_n] = DEFAULT_NULL_VALUE
    else:
        # filter the rows and cols in a symmetric matrix
        for row_n in filtered_rows:
            coevolution_matrix[row_n, :] = coevolution_matrix[
                :, row_n
            ] = DEFAULT_NULL_VALUE

    # return the result
    return coevolution_matrix


def is_parsimony_informative(
    column_freqs,
    minimum_count=2,
    minimum_differences=2,
    ignored=DEFAULT_EXCLUDES,
    strict=False,
):
    """Return True is aln_position is parsimony informative

    column_freqs: dict of characters at alignmnet position mapped
     to their counts -- this is the output of call alignment.column_freqs()
    minimum_count: the minimum number of times a character must show up
     for it to be acceptable (default: 2)
    minimum_differences: the minimum number of different characters
     that must show up at the alignment position (default: 2)
    ignored: characters that should not be counted toward
     minimum_differences (default are exclude characters)
    strict: if True, requires that all amino acids showing up at least
     once at the alignment position show up at least minimum_counts
     times, rather than only requiring that minimum_differences
     amino acids show up minimum_counts times. (default: False)

    The term parsimony informative comes from Codoner, O'Dea,
     and Fares 2008, Reducing the false positive rate in the non-
     parametric analysis of molecular coevolution. In the paper
     they find that if positions which don't contain at least two
     different amino acids, and where each different amino acid doesnt
     show up at least twice each are ignored (i.e., treated as though
     there is not enough information) the positive predictive value
     (PPV) and sensitivity (SN) increase on simulated alignments. They
     term this quality parsimony informative.
     I implemented this as a filter, but include some generalization.
     To determine if a column in an alignment is parsimony informative
     in the exact manner described in Codoner et al., the following
     parameter settings are required:
      minimum_count = 2 (default)
      minimum_differences = 2 (default)
      strict = True (default is False)
     To generalize this function, minimum_count and minimum_differences
     can be passed in so at least minimum_differences different amino
     acids must show up, and each amino acid must show up at least
     minimum_count times.
     In additional variation, strict=False can be passed requiring
     that only minimum_differences number of amino acids show up at least
     minimum_counts times (opposed to requiring that ALL amino acids show
     up minimum_counts times). This is the default behavior.
     By default, the default exclude characters (- and ?) don't count.

    """
    try:
        column_freqs = column_freqs.to_dict()
    except AttributeError:
        pass
    ignored = None if not ignored else list(set(ignored) & set(column_freqs.keys()))
    if ignored:
        for e_ in ignored:
            try:
                del column_freqs[e_]
            except KeyError:
                pass

    if len(column_freqs) < minimum_differences:
        return False
    count_gte_minimum = 0
    for count in list(column_freqs.values()):
        # if not strict, only minimum_differences of the counts
        # must be greater than or equal to minimum_count, so
        # count those occurrences (this is different than the
        # exact technique presented in Codoner et al.)
        if count >= minimum_count:
            count_gte_minimum += 1
        # if strict, all counts must be greater than minimum_count
        # so return False here if we find one that isn't. This is how
        # the algorithm is described in Codoner et al.
        elif strict:
            return False
    return count_gte_minimum >= minimum_differences


def filter_non_parsimony_informative(
    aln,
    coevolution_matrix,
    null_value=DEFAULT_NULL_VALUE,
    minimum_count=2,
    minimum_differences=2,
    ignored=DEFAULT_EXCLUDES,
    intermolecular_data_only=False,
    strict=False,
):
    """Replaces scores in coevolution_matrix with null_value for positions
    which are not parsimony informative.

    See is_parsimony_informative doc string for definition of
     parsimony informative.

    aln: the input alignment used to generate the coevolution matrix;
     if the alignment was recoded, this should be the recoded alignment.
    coevolution_matrix: the result matrix
    null_value: the value to place in positions which are not
     parsimony informative
    """
    if intermolecular_data_only:
        len_aln1 = coevolution_matrix.shape[1]
    column_frequencies = aln.counts_per_pos()
    for i in range(len(column_frequencies)):
        if not is_parsimony_informative(
            column_frequencies[i], minimum_count, minimum_differences, ignored, strict
        ):
            if not intermolecular_data_only:
                coevolution_matrix[i, :] = coevolution_matrix[:, i] = null_value
            else:
                try:
                    coevolution_matrix[:, i] = null_value
                except IndexError:
                    coevolution_matrix[i - len_aln1, :] = null_value


def make_positional_exclude_percentage_function(excludes, max_exclude_percent):
    """return function to identify aln positions with > max_exclude_percent"""
    excludes = {}.fromkeys(excludes)

    def f(col):
        exclude_count = 0
        for c in col:
            if c in excludes:
                exclude_count += 1
        return exclude_count / len(col) > max_exclude_percent

    return f


def filter_exclude_positions(
    aln,
    coevolution_matrix,
    max_exclude_percent=0.1,
    null_value=DEFAULT_NULL_VALUE,
    excludes=DEFAULT_EXCLUDES,
    intermolecular_data_only=False,
):
    """Assign null_value to positions with > max_exclude_percent excludes

    aln: the ArrayAlignment object
    coevolution_matrix: the 2D numpy array -- this will be modified
    max_exclude_percent: the maximimu percent of characters that
     may be exclude characters in any alignment position (column).
     if the percent of exclude characters is greater than this value,
     values in this position will be replaced with null_value
     (default = 0.10)
    null_value: the value to be used as null (default: gDefaultNullValue)
    excludes: the exclude characters (default: gDefaultExcludes)
    intermolecular_data_only: True if the coevolution result
     matrix contains only intermolecular data (default: False)

    """
    # construct the function to be passed to aln.get_position_indices
    f = make_positional_exclude_percentage_function(excludes, max_exclude_percent)
    # identify the positions containing too many exclude characters
    exclude_positions = aln.get_position_indices(f)

    # replace values from exclude_positions with null_value
    if not intermolecular_data_only:
        # if working with intramolecular data (or inter + intra molecular data)
        # this is easy
        for p in exclude_positions:
            coevolution_matrix[p, :] = coevolution_matrix[:, p] = null_value
    else:
        # if working with intermolecular data only, this is more complicated --
        # must convert from alignment positions to matrix positions
        len_aln1 = coevolution_matrix.shape[1]
        for p in exclude_positions:
            try:
                coevolution_matrix[:, p] = null_value
            except IndexError:
                coevolution_matrix[p - len_aln1, :] = null_value


# Functions for archiving/retrieiving coevolve results
# These functions are extremely general -- should they go
# somewhere else, or should I be using pre-existing code?


def pickle_coevolution_result(coevolve_result, out_filepath="output.pkl"):
    """Pickle coevolve_result and store it at output_filepath

    coevolve_result: result from a coevolve_* function (above); this can
     be a float, an array, or a 2D array (most likely it will be one of the
     latter two, as it will usually be fast enough to compute a single
     coevolve value on-the-fly.
    out_filepath: path where the pickled result should be stored
    """
    try:
        infile = open(out_filepath, "wb")
        p = Pickler(infile)
    except IOError:
        err = "Can't access filepath. Do you have write access? " + out_filepath
        raise IOError(err)
    p.dump(coevolve_result)
    infile.close()


def unpickle_coevolution_result(in_filepath):
    """Read in coevolve_result from a pickled file

    in_filepath: filepath to unpickle
    """
    try:
        infile = open(in_filepath, "rb")
        u = Unpickler(infile)
    except IOError:
        err = (
            "Can't access filepath. Does it exist? Do you have read access? "
            + in_filepath
        )
        raise IOError(err)
    r = u.load()
    infile.close()
    return r


def coevolution_matrix_to_csv(coevolve_matrix, out_filepath="output.csv"):
    """Write coevolve_matrix as csv file at output_filepath

    coevolve_result: result from a coevolve_alignment function (above);
     this should be a 2D numpy array
    out_filepath: path where the csv result should be stored
    """
    try:
        f = open(out_filepath, "w")
    except IOError:
        err = "Can't access filepath. Do you have write access? " + out_filepath
        raise IOError(err)
    f.write("\n".join([",".join([str(v) for v in row]) for row in coevolve_matrix]))
    f.close()


def csv_to_coevolution_matrix(in_filepath):
    """Read a coevolution matrix from a csv file

    in_filepath: input filepath
    """
    try:
        f = open(in_filepath)
    except IOError:
        err = (
            "Can't access filepath. Does it exist? Do you have read access? "
            + in_filepath
        )
        raise IOError(err)
    result = []
    for line in f:
        values = line.strip().split(",")
        result.append(list(map(float, values)))
    f.close()
    return array(result)


# End functions for archiving/retrieiving coevolve results

# Start functions for analyzing the results of a coevolution run.


def identify_aln_positions_above_threshold(
    coevolution_matrix, threshold, aln_position, null_value=DEFAULT_NULL_VALUE
):
    """Returns the list of alignment positions which achieve a
    score >= threshold with aln_position.
    Coevolution matrix should be symmetrical or you
    may get weird results -- scores are pulled from the row describing
    aln_position.
    """
    coevolution_scores = coevolution_matrix[aln_position]
    results = []
    for i in range(len(coevolution_scores)):
        s = coevolution_scores[i]
        if s != null_value and s >= threshold:
            results.append(i)
    return results


def aln_position_pairs_cmp_threshold(
    coevolution_matrix,
    threshold,
    cmp_function,
    null_value=DEFAULT_NULL_VALUE,
    intermolecular_data_only=False,
):
    """Returns list of position pairs with score >= threshold

    coevolution_matrix: 2D numpy array
    threshold: value to compare matrix positions against
    cmp_function: function which takes a value and theshold
     and returns a boolean (e.g., ge(), le())
    null_value: value representing null scores -- these are
     ignored
    intermolecular_data_only: True if the coevolution result
     matrix contains only intermolecular data (default: False)
    """
    if not intermolecular_data_only:
        assert (
            coevolution_matrix.shape[0] == coevolution_matrix.shape[1]
        ), "Non-square matrices only supported for intermolecular-only data."
    results = []
    # compile the matrix positions with cmp(value,threshold) == True
    for i, row in enumerate(coevolution_matrix):
        for j, value in enumerate(row):
            if value != null_value and cmp_function(value, threshold):
                results.append((i, j))

    # if working with intermolecular data only, need to convert
    # matrix positions to alignment positions
    if intermolecular_data_only:
        # convert matrix positions to alignment positions
        adjustment = coevolution_matrix.shape[1]
        results = [(j, i + adjustment) for i, j in results]
    return results


def aln_position_pairs_ge_threshold(
    coevolution_matrix,
    threshold,
    null_value=DEFAULT_NULL_VALUE,
    intermolecular_data_only=False,
):
    """wrapper function for aln_position_pairs_cmp_threshold"""
    return aln_position_pairs_cmp_threshold(
        coevolution_matrix,
        threshold,
        greater_equal,
        null_value,
        intermolecular_data_only,
    )


def aln_position_pairs_le_threshold(
    coevolution_matrix,
    threshold,
    null_value=DEFAULT_NULL_VALUE,
    intermolecular_data_only=False,
):
    """wrapper function for aln_position_pairs_cmp_threshold"""
    return aln_position_pairs_cmp_threshold(
        coevolution_matrix, threshold, less_equal, null_value, intermolecular_data_only
    )


def count_cmp_threshold(
    m,
    threshold,
    cmp_function,
    null_value=DEFAULT_NULL_VALUE,
    symmetric=False,
    ignore_diagonal=False,
):
    """Returns a count of the values in m >= threshold, ignoring nulls.

    m: coevolution matrix (numpy array)
    thresold: value to compare against scores in matrix (float)
    cmp_function: function used to compare value to threshold
     (e.g., greater_equal, less_equal)
    """

    total_non_null = 0
    total_hits = 0
    if not symmetric:
        if ignore_diagonal:
            values = [
                m[i, j] for i in range(m.shape[0]) for j in range(m.shape[1]) if i != j
            ]
        else:
            values = m.flat
    else:
        if ignore_diagonal:
            # has to be a better way to do this... tril doesn't work b/c it
            # sets the upper triangle to zero -- if i could get it to set
            # that to null_value, and then apply flat, that'd be fine.
            # values = tril(m,-1)
            values = [m[i, j] for i in range(len(m)) for j in range(i)]
        else:
            # values = tril(m)
            values = [m[i, j] for i in range(len(m)) for j in range(i + 1)]

    if isnan(null_value):

        def is_not_null_value(v):
            return not isnan(v)

    else:

        def is_not_null_value(v):
            return isnan(v) or v != null_value

    for value in values:
        if is_not_null_value(value):
            total_non_null += 1
            if cmp_function(value, threshold):
                total_hits += 1
    return total_hits, total_non_null


def count_ge_threshold(
    m, threshold, null_value=DEFAULT_NULL_VALUE, symmetric=False, ignore_diagonal=False
):
    """wrapper function for count_cmp_threshold"""
    return count_cmp_threshold(
        m, threshold, greater_equal, null_value, symmetric, ignore_diagonal
    )


def count_le_threshold(
    m, threshold, null_value=DEFAULT_NULL_VALUE, symmetric=False, ignore_diagonal=False
):
    """wrapper function for count_cmp_threshold"""
    return count_cmp_threshold(
        m, threshold, less_equal, null_value, symmetric, ignore_diagonal
    )


def ltm_to_symmetric(m):
    """Copies values from lower triangle to upper triangle"""
    assert (
        m.shape[0] == m.shape[1]
    ), "Making matrices symmetric only supported for square matrices"

    for i in range(len(m)):
        for j in range(i):
            m[j, i] = m[i, j]
    return m


# End functions for analyzing the results of a coevolution run


# Script functionality
def build_coevolution_matrix_filepath(
    input_filepath, output_dir="./", method=None, alphabet=None, parameter=None
):
    """ Build filepath from input filename, output dir, and list of suffixes

        input_filepath: filepath to be used for generating the output
            filepath. The path and the final suffix will be stripped to
            get the 'base' filename.
        output_dir: the path to append to the beginning of the base filename
        method: string indicating method that should be appended to filename
        alphabet: string indicating an alphabet recoding which should be
            appended to filename, or None
        parameter: parameter that should be appended to the filename,
            or None (ignored if method doesn't require parameter)

        Examples:
         >>> build_coevolution_matrix_filepath(\
          './p53.fasta','/output/path','mi','charge')
         /output/path/p53.charge.mi
         >>> build_coevolution_matrix_filepath(\
          './p53.new.fasta','/output/path','mi','charge')
         /output/path/p53.new.charge.mi
         >>> build_coevolution_matrix_filepath(\
          './p53.fasta','/output/path','sca','charge',0.75)
         /output/path/p53.charge.sca_75

    """
    if method == "sca":
        try:
            cutoff_str = str(parameter)
            point_index = cutoff_str.rindex(".")
            method = "_".join([method, cutoff_str[point_index + 1 : point_index + 4]])
        except ValueError:
            raise ValueError("Cutoff must be provided when method == 'sca'")

    suffixes = [_f for _f in [alphabet, method] if _f]

    # strip path
    try:
        result = input_filepath[input_filepath.rindex("/") + 1 :]
    except ValueError:
        result = input_filepath
    # strip final suffix
    try:
        result = result[: result.rindex(".")]
    except ValueError:
        pass
    # append output path
    if output_dir.endswith("/"):
        result = "".join([output_dir, result])
    else:
        result = "".join([output_dir, "/", result])
    # append output suffixes
    result = ".".join([_f for _f in [result] + suffixes if _f])

    return result


def parse_coevolution_matrix_filepath(filepath):
    """Parses a coevolution matrix filepath into constituent parts.

    Format is very specific. Will only work on filenames such as:
     path/alignment_identifier.alphabet_id.method.pkl
     path/alignment_identifier.alphabet_id.method.csv

    This format is the recommended naming convention for coevolution
     matrices. To ensure filepaths compatible with this function, use
     cogent3.evolve.coevolution.build_coevolution_matrix_filepath to build
     the filepaths for your coevolution matrices.


     Examples:
     parse_coevolution_matrix_filepath('pkls/myosin_995.a1_4.nmi.pkl')
        => ('myosin_995', 'a1_4', 'nmi')
     parse_coevolution_matrix_filepath('p53.orig.mi.csv')
        => ('p53','orig','mi')
    """
    filename = basename(filepath)
    fields = filename.split(".")
    try:
        alignment_id = fields[0]
        alphabet_id = fields[1]
        method_id = fields[2]
        _ = fields[3]  # extension
    except IndexError:
        raise ValueError(
            "output filepath not in parsable format: %s. See doc string for format definition."
            % filepath
        )

    return (alignment_id, alphabet_id, method_id)
