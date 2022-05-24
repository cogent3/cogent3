"""This file contains functions for recoding alignment objects with
 reduced-state alphabets, and also defines some reduced-state alphabets.

Example alphabet definitions:
charge_his_3: Three states for +/-/no charge, with histidine counted as
 a charged residue.
size_2: Two states for large/small. Splits were manually derived by
    ordering the Mr of the residues and finding a natural split in the
    middle of the data set (residues with Mr <= 133 -> small;
    Mr >= 146 -> large).
    Ordered Mr: [75,89,105,115,117,119,121,131,131,132,133,146,146,147,
                 149,155,165,174,181,204]
    Differences between neighboring Mr values: [14,16,10,2,2,2,10,0,1,1,13,
                                                0,1,2,6,10,9,7,23]
    The difference between 133 and 146 (or D and K/Q) seems like the most
        natural split.

Alphabet naming convention: standard alphabets (those defined in this file)
 are named based on an identifier (e.g., charge_his) followed by an underscore,
 followed by the number of states in the reduced alphabet (e.g., 3).

Alphabet definition convention: When defining reduced alphabets, a reduced
 state should be represented by the first char listed in that state. This
 allows for alignments to not have to support new characeters (by changing
 the MolType of the alignment) and it allows for easier interpretation of the
 alignment states.

Many of the alphabets presented here are discussed in:
Detecting Coevolution by Disregarding Evolution? Tree-Ignorant Metrics of
Coevolution Perform As Well As Tree-Aware Metrics; J. Gregory Caporaso,
Sandra Smit, Brett C. Easton, Lawrence Hunter, Gavin A. Huttley, and
Rob Knight. BMC Evolutionary Biology, 2008.


"""


from numpy import array, take, zeros

from cogent3.core.alignment import ArrayAlignment
from cogent3.evolve.models import DSO78_freqs, DSO78_matrix


__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Greg Caporaso"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Beta"


class RecodeError(Exception):
    """A generic error to be raised when errors occur in recoding"""

    pass


alphabets = {
    # Note: 3-STATE CHARGE ALPHAS ASSIGN B/Z TO UNCHARGED -- I FIGURE THAT'S
    # SAFER THAT ASSIGNING THEM TO EITHER +/-, BUT MAYBE THAT'S NOT ACCURATE.
    # AN ALTERNATIVE WOULD BE TO ASSIGN THEM TO WHATEVER THEIR MORE COMMON
    # VALUE WAS IN (e.g.) ALL PROTEINS
    "charge_2": [("K", "KRDEBZ"), ("A", "ACFGHILMNPQSTVWYX")],
    "charge_3": [("K", "KR"), ("D", "DE"), ("A", "ACFGHILMNPQSTVWYBXZ")],
    "charge_his_2": [("K", "KRDEHBZ"), ("A", "ACFGILMNPQSTVWYX")],
    "charge_his_3": [("K", "KRH"), ("D", "DE"), ("A", "ACFGILMNPQSTVWYBXZ")],
    "size_2": [("G", "GAVLISPTCNDXB"), ("M", "MFYWQKHREZ")],
    "hydropathy_3": [("R", "RKDENQHBZ"), ("Y", "YWSTG"), ("P", "PAMCFLVIX")],
    # B/Z are assigned to the acidic residues, since D/E are more common --
    # need to figure out if this is the best way to handle them
    "polarity_his_4": [
        ("D", "DEBZ"),
        ("R", "RHK"),
        ("A", "AILMFPWV"),
        ("G", "GSTCYNQ"),
    ],
    # This is a modified A1_4 alphabet to capture natural breaks in the metric
    "a1_4m": [("C", "CVILF"), ("M", "MWA"), ("T", "GSTPYH"), ("N", "QNDERK")],
    "a1_2": [("C", "CVILFMWAGS"), ("T", "TPYHQNDERK")],
    "a1_3": [("C", "CVILFMW"), ("A", "AGSTPY"), ("H", "HQNDERK")],
    "a1_4": [("C", "CVILF"), ("M", "MWAGS"), ("T", "TPYHQ"), ("N", "NDERK")],
    "a1_5": [("C", "CVIL"), ("F", "FMWA"), ("G", "GSTP"), ("Y", "YHQN"), ("D", "DERK")],
    "a1_6": [
        ("C", "CVI"),
        ("L", "LFMW"),
        ("A", "AGS"),
        ("T", "TPY"),
        ("H", "HQND"),
        ("E", "ERK"),
    ],
    "a1_7": [
        ("C", "CVI"),
        ("L", "LFM"),
        ("W", "WAG"),
        ("S", "ST"),
        ("P", "PYH"),
        ("Q", "QND"),
        ("E", "ERK"),
    ],
    "a1_8": [
        ("C", "CVI"),
        ("L", "LF"),
        ("M", "MWA"),
        ("G", "GS"),
        ("T", "TPY"),
        ("H", "HQ"),
        ("N", "NDE"),
        ("R", "RK"),
    ],
    "a1_9": [
        ("C", "CV"),
        ("I", "IL"),
        ("F", "FMW"),
        ("A", "AG"),
        ("S", "ST"),
        ("P", "PY"),
        ("H", "HQN"),
        ("D", "DE"),
        ("R", "RK"),
    ],
    "a1_10": [
        ("C", "CV"),
        ("I", "IL"),
        ("F", "FM"),
        ("W", "WA"),
        ("G", "GS"),
        ("T", "TP"),
        ("Y", "YH"),
        ("Q", "QN"),
        ("D", "DE"),
        ("R", "RK"),
    ],
    "a2_2": [("M", "MEALFKIHVQ"), ("R", "RWDTCNYSGP")],
    "a2_3": [("M", "MEALFKI"), ("H", "HVQRWD"), ("T", "TCNYSGP")],
    "a2_4": [("M", "MEALF"), ("K", "KIHVQ"), ("R", "RWDTC"), ("N", "NYSGP")],
    "a2_5": [("M", "MEAL"), ("F", "FKIH"), ("V", "VQRW"), ("D", "DTCN"), ("Y", "YSGP")],
    "a2_6": [
        ("M", "MEA"),
        ("L", "LFKI"),
        ("H", "HVQ"),
        ("R", "RWD"),
        ("T", "TCNY"),
        ("S", "SGP"),
    ],
    "a2_7": [
        ("M", "MEA"),
        ("L", "LFK"),
        ("I", "IHV"),
        ("Q", "QR"),
        ("W", "WDT"),
        ("C", "CNY"),
        ("S", "SGP"),
    ],
    "a2_8": [
        ("M", "MEA"),
        ("L", "LF"),
        ("K", "KIH"),
        ("V", "VQ"),
        ("R", "RWD"),
        ("T", "TC"),
        ("N", "NYS"),
        ("G", "GP"),
    ],
    "a2_9": [
        ("M", "ME"),
        ("A", "AL"),
        ("F", "FKI"),
        ("H", "HV"),
        ("Q", "QR"),
        ("W", "WD"),
        ("T", "TCN"),
        ("Y", "YS"),
        ("G", "GP"),
    ],
    "a2_10": [
        ("M", "ME"),
        ("A", "AL"),
        ("F", "FK"),
        ("I", "IH"),
        ("V", "VQ"),
        ("R", "RW"),
        ("D", "DT"),
        ("C", "CN"),
        ("Y", "YS"),
        ("G", "GP"),
    ],
    "a3_2": [("S", "SDQHPLCAVK"), ("W", "WNGERFITMY")],
    "a3_3": [("S", "SDQHPLC"), ("A", "AVKWNG"), ("E", "ERFITMY")],
    "a3_4": [("S", "SDQHP"), ("L", "LCAVK"), ("W", "WNGER"), ("F", "FITMY")],
    "a3_5": [("S", "SDQH"), ("P", "PLCA"), ("V", "VKWN"), ("G", "GERF"), ("I", "ITMY")],
    "a3_6": [
        ("S", "SDQ"),
        ("H", "HPLC"),
        ("A", "AVK"),
        ("W", "WNG"),
        ("E", "ERFI"),
        ("T", "TMY"),
    ],
    "a3_7": [
        ("S", "SDQ"),
        ("H", "HPL"),
        ("C", "CAV"),
        ("K", "KW"),
        ("N", "NGE"),
        ("R", "RFI"),
        ("T", "TMY"),
    ],
    "a3_8": [
        ("S", "SDQ"),
        ("H", "HP"),
        ("L", "LCA"),
        ("V", "VK"),
        ("W", "WNG"),
        ("E", "ER"),
        ("F", "FIT"),
        ("M", "MY"),
    ],
    "a3_9": [
        ("S", "SD"),
        ("Q", "QH"),
        ("P", "PLC"),
        ("A", "AV"),
        ("K", "KW"),
        ("N", "NG"),
        ("E", "ERF"),
        ("I", "IT"),
        ("M", "MY"),
    ],
    "a3_10": [
        ("S", "SD"),
        ("Q", "QH"),
        ("P", "PL"),
        ("C", "CA"),
        ("V", "VK"),
        ("W", "WN"),
        ("G", "GE"),
        ("R", "RF"),
        ("I", "IT"),
        ("M", "MY"),
    ],
    "a4_2": [("W", "WHCMYQFKDN"), ("E", "EIPRSTGVLA")],
    "a4_3": [("W", "WHCMYQF"), ("K", "KDNEIP"), ("R", "RSTGVLA")],
    "a4_4": [("W", "WHCMY"), ("Q", "QFKDN"), ("E", "EIPRS"), ("T", "TGVLA")],
    "a4_5": [("W", "WHCM"), ("Y", "YQFK"), ("D", "DNEI"), ("P", "PRST"), ("G", "GVLA")],
    "a4_6": [
        ("W", "WHC"),
        ("M", "MYQF"),
        ("K", "KDN"),
        ("E", "EIP"),
        ("R", "RSTG"),
        ("V", "VLA"),
    ],
    "a4_7": [
        ("W", "WHC"),
        ("M", "MYQ"),
        ("F", "FKD"),
        ("N", "NE"),
        ("I", "IPR"),
        ("S", "STG"),
        ("V", "VLA"),
    ],
    "a4_8": [
        ("W", "WHC"),
        ("M", "MY"),
        ("Q", "QFK"),
        ("D", "DN"),
        ("E", "EIP"),
        ("R", "RS"),
        ("T", "TGV"),
        ("L", "LA"),
    ],
    "a4_9": [
        ("W", "WH"),
        ("C", "CM"),
        ("Y", "YQF"),
        ("K", "KD"),
        ("N", "NE"),
        ("I", "IP"),
        ("R", "RST"),
        ("G", "GV"),
        ("L", "LA"),
    ],
    "a4_10": [
        ("W", "WH"),
        ("C", "CM"),
        ("Y", "YQ"),
        ("F", "FK"),
        ("D", "DN"),
        ("E", "EI"),
        ("P", "PR"),
        ("S", "ST"),
        ("G", "GV"),
        ("L", "LA"),
    ],
    "a5_2": [("D", "DSQPVLECWA"), ("H", "HFINMTYKGR")],
    "a5_3": [("D", "DSQPVLE"), ("C", "CWAHFI"), ("N", "NMTYKGR")],
    "a5_4": [("D", "DSQPV"), ("L", "LECWA"), ("H", "HFINM"), ("T", "TYKGR")],
    "a5_5": [("D", "DSQP"), ("V", "VLEC"), ("W", "WAHF"), ("I", "INMT"), ("Y", "YKGR")],
    "a5_6": [
        ("D", "DSQ"),
        ("P", "PVLE"),
        ("C", "CWA"),
        ("H", "HFI"),
        ("N", "NMTY"),
        ("K", "KGR"),
    ],
    "a5_7": [
        ("D", "DSQ"),
        ("P", "PVL"),
        ("E", "ECW"),
        ("A", "AH"),
        ("F", "FIN"),
        ("M", "MTY"),
        ("K", "KGR"),
    ],
    "a5_8": [
        ("D", "DSQ"),
        ("P", "PV"),
        ("L", "LEC"),
        ("W", "WA"),
        ("H", "HFI"),
        ("N", "NM"),
        ("T", "TYK"),
        ("G", "GR"),
    ],
    "a5_9": [
        ("D", "DS"),
        ("Q", "QP"),
        ("V", "VLE"),
        ("C", "CW"),
        ("A", "AH"),
        ("F", "FI"),
        ("N", "NMT"),
        ("Y", "YK"),
        ("G", "GR"),
    ],
    "a5_10": [
        ("D", "DS"),
        ("Q", "QP"),
        ("V", "VL"),
        ("E", "EC"),
        ("W", "WA"),
        ("H", "HF"),
        ("I", "IN"),
        ("M", "MT"),
        ("Y", "YK"),
        ("G", "GR"),
    ],
    # orig does no recoding, but is provided for convenience so if you want to
    # iterate over all reduced alphabets and the full alphabet, you can do that
    # without having specify the original alphabet differently.
    "orig": list(zip("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY")),
}


def build_alphabet_map(alphabet_id=None, alphabet_def=None):
    """return dict mapping old alphabet chars to new alphabet chars

    alphabet_id: string identifying an alphabet in
        cogent3.util.recode_alignment.alphabets.
        (See cogent3.util.recode_alignment.alphabets.keys()
        for valid alphabet_ids.)
    alphabet_def: list of two-element tuples where first element is
        the new alphabet character and the second elements is an iterable
        object containing the old alphabet chars which should be mapped to
        the new char.
        e.g., [('A','CVILFMWAGSTPYH'),('B','QNDERKBZ')]
        (See cogent3.util.recode_alignment.alphabets.values()
        for more examples.)

    NOTE: Only one of the two parameters should be provided -- you either
        provide the alphabet, or it is looked up. If you do provide both,
        the alphabet_id is ignored.

    """
    try:
        alphabet_def = alphabet_def or alphabets[alphabet_id]
    except KeyError:
        if not alphabet_id:
            raise ValueError("Must provide an alphabet_id or alphabet definiton.")
        raise ValueError("Invalid alphabet id.")

    result = {}
    for new, old in alphabet_def:
        for old_c in old:
            result[old_c] = new

    return result


def recode_dense_alignment(aln, alphabet_id=None, alphabet_def=None):
    """Return new ArrayAlignment recoded in the provided reduced-state alphabet

    aln: the ArrayAlignment object to be recoded
    alphabet_id: string identifying an alphabet in
        cogent3.util.recode_alignment.alphabets.
        (See cogent3.util.recode_alignment.alphabets.keys()
        for valid alphabet_ids.)
    alphabet_def: list of two-element tuples where first element is
        the new alphabet character and the second elements is an iterable
        object containing the old alphabet chars which should be mapped to
        the new char.
        e.g., [('A','CVILFMWAGSTPYH'),('B','QNDERKBZ')]
        (See cogent3.util.recode_alignment.alphabets.values()
        for more examples.)

    Note: either alphabet_id OR alphabet_def must be passed. Either
        provide the alphabet, or have it is looked up. If both are provided
        the alphabet_id is ignored.

    """

    # Construct a dict mapping from UInt8s in alignment to their
    # associated characters. This dict is then used for looking
    # up chars in the new and old alphabets.
    byte_map = dict(list(zip(aln.alphabet, list(range(len(aln.alphabet))))))

    # Construct a dict mapping old characters to new characters.
    alphabet_map = build_alphabet_map(
        alphabet_id=alphabet_id, alphabet_def=alphabet_def
    )

    # Create the recoded version of seqs.alphabet
    new_indices = list(range(len(aln.alphabet)))
    for old, new in list(alphabet_map.items()):
        new_indices[byte_map[old]] = byte_map[new]

    # Map the old alphabet onto the new alphabet. Note: characters that
    # that are not mapped are ignored. Returns a new ArrayAlignment.
    return ArrayAlignment(
        take(new_indices, aln.array_seqs).transpose(), aln.names[:], moltype=aln.moltype
    )


def recode_alignment(aln, alphabet_id=None, alphabet_def=None):
    raise NotImplementedError


def recode_freq_vector(alphabet_def, freqs, ignores="BXZ"):
    """recode the bg_freqs to reflect the recoding defined in alphabet_def

    alphabet_def: list of tuples where new char is first tuple element
        and sequence of old chars is second tuple element. (For examples,
        see cogent3.util.recode_alignment.alphabets.values())
    freqs: dict mapping chars to their frequencies
    ignores: the degenerate characters -- we don't want to include these
        in the new freqs, b/c they'll be counted multiple times. Also,
        isn't clear what should be done if an alphabet were to split them
        apart.

    Note: there is no error-checking here, so you need to be sure that
     the alphabet and the frequencies are compatible (i.e., freqs and the
     old characters must overlap perfectly, with the exception of the
     degenerate characters, which are ignored by default).
    """
    result = {}
    for new, olds in alphabet_def:
        for old in olds:
            if old in ignores:
                continue
            try:
                result[new] += freqs[old]
            except KeyError:
                result[new] = freqs[old]
    return result


# The following code is for recoding substitution matrices


def square_matrix_to_dict(matrix, key_order="ACDEFGHIKLMNPQRSTVWY"):
    result = {}
    for c, row in zip(key_order, matrix):
        result[c] = dict(list(zip(key_order, row)))
    return result


def recode_count_matrix(alphabet, count_matrix, aa_order):
    """Recodes a subsitution count matrix

    alphabet: the alphabet to be used for recoding the matrix
     (see cogent3.util.recode_alignment.alphabets.values()) for
     examples
    count_matrix: matrix to be recoded (e.g.,
     cogent3.evolve.models.DSO78_matrix)
    aa_order: the order of the rows/cols in the matrix as a string
     (for cogent3.evolve.models.DSO78_matrix this would be
     'ACDEFGHIKLMNPQRSTVWY')

    """
    m = square_matrix_to_dict(count_matrix, aa_order)
    result = zeros(len(aa_order) ** 2).reshape(len(aa_order), len(aa_order))
    result = square_matrix_to_dict(result, aa_order)
    for row_new, row_olds in alphabet:
        for col_new, col_olds in alphabet:
            if row_new not in col_olds:
                new_count = 0.0
                for row_old in row_olds:
                    for col_old in col_olds:
                        try:
                            new_count += m[row_old][col_old]
                        except KeyError:
                            # hit a char that's not in the sub matrix --
                            # probablyan ambiguous residue (i.e., B, X, or Z)
                            pass
                result[row_new][col_new] = new_count
    cm = []
    for row_c in aa_order:
        r = []
        for col_c in aa_order:
            r.append(result[row_c][col_c])
        cm.append(r)
    return array(cm)


def recode_counts_and_freqs(
    alphabet,
    count_matrix=DSO78_matrix,
    freqs=DSO78_freqs,
    aa_order="ACDEFGHIKLMNPQRSTVWY",
):
    """recode a substituion count matrix and a vector of character freqs"""

    recoded_freqs = recode_freq_vector(alphabet, freqs)
    for aa in aa_order:
        if aa not in recoded_freqs:
            recoded_freqs[aa] = 0.0
    recoded_counts = recode_count_matrix(alphabet, count_matrix, aa_order)
    return recoded_counts, recoded_freqs
