"""
moltype.py

MolType provides services for resolving ambiguities, or providing the
correct ambiguity for recoding. It also maintains the mappings between
different kinds of alphabets, sequences and alignments.

One issue with MolTypes is that they need to know about Sequence, Alphabet,
and other objects, but, at the same time, those objects need to know about
the MolType. It is thus essential that the connection between these other
types and the MolType can be made after the objects are created.
"""

__author__ = "Peter Maxwell, Gavin Huttley and Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight", "Daniel McDonald"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


import json
import re

from collections import defaultdict
from copy import deepcopy
from random import choice
from string import ascii_letters as letters

import numpy

from cogent3.core.alignment import (
    Alignment,
    ArrayAlignment,
    SequenceCollection,
)
from cogent3.core.alphabet import (
    Alphabet,
    AlphabetError,
    CharAlphabet,
    Enumeration,
    _make_complement_array,
)
from cogent3.core.genetic_code import get_code
from cogent3.core.sequence import (
    ABSequence,
    ArrayDnaCodonSequence,
    ArrayDnaSequence,
    ArrayProteinSequence,
    ArrayProteinWithStopSequence,
    ArrayRnaCodonSequence,
    ArrayRnaSequence,
    ArraySequence,
    ByteSequence,
    DnaSequence,
    NucleicAcidSequence,
    ProteinSequence,
    ProteinWithStopSequence,
    RnaSequence,
)
from cogent3.core.sequence import Sequence as DefaultSequence
from cogent3.data.molecular_weight import DnaMW, ProteinMW, RnaMW
from cogent3.util.misc import (
    FunctionWrapper,
    add_lowercase,
    get_object_provenance,
    iterable,
)
from cogent3.util.transform import KeepChars, first_index_in_set


Float = numpy.core.numerictypes.sctype2char(float)
Int = numpy.core.numerictypes.sctype2char(int)

maketrans = str.maketrans
translate = str.translate

IUPAC_gap = "-"

IUPAC_missing = "?"

IUPAC_DNA_chars = ["T", "C", "A", "G"]
IUPAC_DNA_ambiguities = {
    "N": ("A", "C", "T", "G"),
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "W": ("A", "T"),
    "S": ("C", "G"),
    "K": ("T", "G"),
    "M": ("C", "A"),
    "B": ("C", "T", "G"),
    "D": ("A", "T", "G"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
}
IUPAC_DNA_ambiguities_complements = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "-": "-",
    "M": "K",
    "K": "M",
    "N": "N",
    "R": "Y",
    "Y": "R",
    "W": "W",
    "S": "S",
    "X": "X",  # not technically an IUPAC ambiguity, but used by repeatmasker
    "V": "B",
    "B": "V",
    "H": "D",
    "D": "H",
}

IUPAC_DNA_complements = {"A": "T", "C": "G", "G": "C", "T": "A", "-": "-"}

# note change in standard order from DNA
IUPAC_RNA_chars = ["U", "C", "A", "G"]
IUPAC_RNA_ambiguities = {
    "N": ("A", "C", "U", "G"),
    "R": ("A", "G"),
    "Y": ("C", "U"),
    "W": ("A", "U"),
    "S": ("C", "G"),
    "K": ("U", "G"),
    "M": ("C", "A"),
    "B": ("C", "U", "G"),
    "D": ("A", "U", "G"),
    "H": ("A", "C", "U"),
    "V": ("A", "C", "G"),
}

IUPAC_RNA_ambiguities_complements = {
    "A": "U",
    "C": "G",
    "G": "C",
    "U": "A",
    "-": "-",
    "M": "K",
    "K": "M",
    "N": "N",
    "R": "Y",
    "Y": "R",
    "W": "W",
    "S": "S",
    "X": "X",  # not technically an IUPAC ambiguity, but used by repeatmasker
    "V": "B",
    "B": "V",
    "H": "D",
    "D": "H",
}

IUPAC_RNA_complements = {"A": "U", "C": "G", "G": "C", "U": "A", "-": "-"}


# Standard RNA pairing: GU pairs count as 'weak' pairs
RnaStandardPairs = {
    ("A", "U"): True,  # True vs False for 'always' vs 'sometimes' pairing
    ("C", "G"): True,
    ("G", "C"): True,
    ("U", "A"): True,
    ("G", "U"): False,
    ("U", "G"): False,
}

# Watson-Crick RNA pairing only: GU pairs don't count as pairs
RnaWCPairs = {("A", "U"): True, ("C", "G"): True, ("G", "C"): True, ("U", "A"): True}

# RNA pairing with GU counted as standard pairs
RnaGUPairs = {
    ("A", "U"): True,
    ("C", "G"): True,
    ("G", "C"): True,
    ("U", "A"): True,
    ("G", "U"): True,
    ("U", "G"): True,
}

# RNA pairing with GU, AA, GA, CA and UU mismatches allowed as weak pairs
RnaExtendedPairs = {
    ("A", "U"): True,
    ("C", "G"): True,
    ("G", "C"): True,
    ("U", "A"): True,
    ("G", "U"): False,
    ("U", "G"): False,
    ("A", "A"): False,
    ("G", "A"): False,
    ("A", "G"): False,
    ("C", "A"): False,
    ("A", "C"): False,
    ("U", "U"): False,
}
# Standard DNA pairing: only Watson-Crick pairs count as pairs
DnaStandardPairs = {
    ("A", "T"): True,
    ("C", "G"): True,
    ("G", "C"): True,
    ("T", "A"): True,
}


# protein letters & ambiguity codes
IUPAC_PROTEIN_code_aa = {
    "A": "Alanine",
    "C": "Cysteine",
    "D": "Aspartic Acid",
    "E": "Glutamic Acid",
    "F": "Phenylalanine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "K": "Lysine",
    "L": "Leucine",
    "M": "Methionine",
    "N": "Asparagine",
    "P": "Proline",
    "Q": "Glutamine",
    "R": "Arginine",
    "S": "Serine",
    "T": "Threonine",
    "V": "Valine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "*": "STOP",
}

IUPAC_PROTEIN_chars = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "U",
    "V",
    "W",
    "Y",
]

PROTEIN_WITH_STOP_chars = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "U",
    "V",
    "W",
    "Y",
    "*",
]

IUPAC_PROTEIN_ambiguities = {"B": ["N", "D"], "X": IUPAC_PROTEIN_chars, "Z": ["Q", "E"]}

PROTEIN_WITH_STOP_ambiguities = {
    "B": ["N", "D"],
    "X": PROTEIN_WITH_STOP_chars,
    "Z": ["Q", "E"],
}


class FoundMatch(Exception):
    """Raised when a match is found in a deep loop to skip many levels"""

    pass


def make_matches(monomers=None, gaps=None, degenerates=None):
    """Makes a dict of symbol pairs (i,j) -> strictness.

    Strictness is True if i and j always match and False if they sometimes
    match (e.g. A always matches A, but W sometimes matches R).
    """
    result = {}
    # allow defaults to be left blank without problems
    monomers = monomers or {}
    gaps = gaps or {}
    degenerates = degenerates or {}
    # all monomers always match themselves and no other monomers
    for i in monomers:
        result[(i, i)] = True
    # all gaps always match all other gaps
    for i in gaps:
        for j in gaps:
            result[(i, j)] = True
    # monomers sometimes match degenerates that contain them
    for i in monomers:
        for j in degenerates:
            if i in degenerates[j]:
                result[(i, j)] = False
                result[(j, i)] = False
    # degenerates sometimes match degenerates that contain at least one of
    # the same monomers
    for i in degenerates:
        for j in degenerates:
            try:
                for i_symbol in degenerates[i]:
                    if i_symbol in degenerates[j]:
                        result[(i, j)] = False
                        raise FoundMatch
            except FoundMatch:
                pass  # flow control: break out of doubly nested loop
    return result


def make_pairs(pairs=None, monomers=None, gaps=None, degenerates=None):
    """Makes a dict of symbol pairs (i,j) -> strictness.

    Expands pairs into all possible pairs using degen symbols.
    Strictness is True if i and j always pair, and False if they 'weakly' pair
    (e.g. GU pairs or if it is possible that they pair).

    If you want to make GU pairs count as 'always matching', pass in pairs
    that have (G,U) and (U, G) mapped to True rather than False.
    """
    result = {}
    # allow defaults to be left blank without problems
    pairs = pairs or {}
    monomers = monomers or {}
    gaps = gaps or {}
    degenerates = degenerates or {}
    # add in the original pairs: should be complete monomer pairs
    result.update(pairs)
    # all gaps 'weakly' pair with each other
    for i in gaps:
        for j in gaps:
            result[(i, j)] = False
    # monomers sometimes pair with degenerates if the monomer's complement
    # is in the degenerate symbol
    for i in monomers:
        for j in degenerates:
            found = False
            try:
                for curr_j in degenerates[j]:
                    # check if (i,curr_j) and/or (curr_j,i) is a valid pair:
                    # not mutually required if pairs are not all commutative!
                    if (i, curr_j) in pairs:
                        result[(i, j)] = False
                        found = True
                    if (curr_j, i) in pairs:
                        result[(j, i)] = False
                        found = True
                    if found:
                        raise FoundMatch
            except FoundMatch:
                pass  # flow control: break out of nested loop
    # degenerates sometimes pair with each other if the first degenerate
    # contains the complement of one of the bases in the second degenerate
    for i in degenerates:
        for j in degenerates:
            try:
                for curr_i in degenerates[i]:
                    for curr_j in degenerates[j]:
                        if (curr_i, curr_j) in pairs:
                            result[(i, j)] = False
                            raise FoundMatch
            except FoundMatch:
                pass  # just using for flow control
    # don't forget the return value!
    return result


# RnaPairingRules is a dict of {name:(base_pairs,degen_pairs)} where base_pairs
# is a dict with the non-degenerate pairing rules and degen_pairs is a dict with
# both the degenerate and non-degenerate pairing rules.
# NOTE: uses make_pairs to augment the initial dict after construction.
RnaPairingRules = {
    "Standard": RnaStandardPairs,
    "WC": RnaWCPairs,
    "GU": RnaGUPairs,
    "Extended": RnaExtendedPairs,
}

for k, v in list(RnaPairingRules.items()):
    RnaPairingRules[k] = (v, make_pairs(v))


class CoreObjectGroup(object):
    """Container relating gapped, ungapped, degen, and non-degen objects."""

    _types = ["base", "degen", "gap", "degen_gap"]

    def __init__(self, base, degen=None, gapped=None, degen_gapped=None):
        """Returns new CoreObjectGroup. Only base is required"""
        self.base = base
        self.degen = degen
        self.gapped = gapped
        self.degen_gapped = degen_gapped
        self._items = [base, degen, gapped, degen_gapped]
        self._set_relationships()

    def _set_relationships(self):
        """Sets relationships between the different "flavors"."""
        self.base.gapped = self.gapped
        self.base.ungapped = self.base
        self.base.degen = self.degen
        self.base.non_degen = self.base

        statements = [
            "self.degen.gapped = self.degen_gapped",
            "self.degen.ungapped = self.degen",
            "self.degen.degen = self.degen",
            "self.degen.non_degen = self.base",
            "self.gapped.gapped = self.gapped",
            "self.gapped.ungapped = self.base",
            "self.gapped.degen = self.degen_gapped",
            "self.gapped.non_degen = self.gapped",
            "self.degen_gapped.gapped = self.degen_gapped",
            "self.degen_gapped.ungapped = self.degen",
            "self.degen_gapped.degen = self.degen_gapped",
            "self.degen_gapped.non_degen = self.gapped",
        ]
        for s in statements:
            try:
                exec(s)
            except AttributeError:
                pass

    def __getitem__(self, i):
        """Allows container to be indexed into, by type of object (e.g. gap)."""
        return self.__dict__[i]

    def which_type(self, a):
        """Returns the type of an alphabet in self, or None if not present."""
        return self._types[self._items.find(a)]


class AlphabetGroup(CoreObjectGroup):
    """Container relating gapped, ungapped, degen, and non-degen alphabets."""

    def __init__(
        self,
        chars,
        degens,
        gap=IUPAC_gap,
        missing=IUPAC_missing,
        moltype=None,
        constructor=None,
    ):
        """Returns new AlphabetGroup."""
        if constructor is None:
            if max(list(map(len, chars))) == 1:
                constructor = CharAlphabet
                chars = "".join(chars)
                degens = "".join(degens)
            else:
                constructor = Alphabet  # assume multi-char

        super(AlphabetGroup, self).__init__(
            base=constructor(chars, moltype=moltype),
            degen=constructor(chars + degens, moltype=moltype),
            gapped=constructor(chars + gap, gap, moltype=moltype),
            degen_gapped=constructor(
                chars + gap + degens + missing, gap, moltype=moltype
            ),
        )

        self._items = [self.base, self.degen, self.gapped, self.degen_gapped]
        self._set_relationships()
        # set complements if MolType was specified
        if moltype is not None:
            comps = moltype.complements
            for i in self._items:
                i._complement_array = _make_complement_array(i, comps)


# colours for HTML representation


def _expand_colors(base, colors):
    base = base.copy()
    base.update({ch: clr for chars, clr in colors.items() for ch in chars})
    return base


class _DefaultValue:
    def __init__(self, value):
        self.value = value

    def __call__(self):
        return self.value


_gray = _DefaultValue("gray")
_base_colors = defaultdict(_gray)

NT_COLORS = _expand_colors(
    _base_colors, {"A": "#FF0102", "C": "black", "G": "green", "T": "blue", "U": "blue"}
)

AA_COLORS = _expand_colors(
    _base_colors,
    {
        "GAVLI": "#009999",
        "FYW": "#ff6600",
        "CM": "orange",
        "ST": "#009900",
        "KRH": "#FF0102",
        "DE": "blue",
        "NQ": "#993300",
        "P": "#cc0099",
    },
)


class MolType(object):
    """MolType: Handles operations that depend on the sequence type (e.g. DNA).

    The MolType knows how to connect alphabets, sequences, alignments, and so
    forth, and how to disambiguate ambiguous symbols and perform base
    pairing (where appropriate).

    WARNING: Objects passed to a MolType become associated with that MolType,
    i.e. if you pass ProteinSequence to a new MolType you make up, all
    ProteinSequences will now be associated with the new MolType. This may
    not be what you expect. Use preserve_existing_moltypes=True if you
    don't want to reset the moltype.
    """

    def __init__(
        self,
        motifset,
        gap=IUPAC_gap,
        missing=IUPAC_missing,
        gaps=None,
        seq_constructor=None,
        ambiguities=None,
        label=None,
        complements=None,
        pairs=None,
        mw_calculator=None,
        add_lower=False,
        preserve_existing_moltypes=False,
        make_alphabet_group=False,
        array_seq_constructor=None,
        colors=None,
    ):
        """Returns a new MolType object. Note that the parameters are in flux.

        Parameters
        ----------
        motifset
            Alphabet or sequence of items in the default
            alphabet. Does not include degenerates.
        gap
            default gap symbol
        missing
            symbol for missing data
        gaps
            any other symbols that should be treated as gaps (doesn't have
            to include gap or missing; they will be silently added)
        seq_constructor
            Class for constructing sequences.
        ambiguities
            dict of char:tuple, doesn't include gaps (these are
            hard-coded as - and ?, and added later.
        label
            text label, don't know what this is used for. Unnecessary?
        complements
            dict of symbol:symbol showing how the non-degenerate
            single characters complement each other. Used for constructing
            on the fly the complement table, incl. support for must_pair and
            can_pair.
        pairs
            dict in which keys are pairs of symbols that can pair
            with each other, values are True (must pair) or False (might
            pair). Currently, the meaning of GU pairs as 'weak' is conflated
            with the meaning of degenerate symbol pairs (which might pair
            with each other but don't necessarily, depending on how the
            symbol is resolved). This should be refactored.
        mw_calculator
            f(seq) -> molecular weight.
        add_lower
            if True (default: False) adds the lowercase versions of
            everything into the alphabet. Slated for deletion.
        preserve_existing_moltypes
            if True (default: False), does not
            set the MolType of the things added in **kwargs to self.
        make_alphabet_group
            if True, makes an AlphabetGroup relating
            the various alphabets to one another.
        array_seq_constructor
            sequence type for array sequence
        colors
            dict mapping moltype characters to colors for display

        Note on "degenerates" versus "ambiguities": self.degenerates contains
        _only_ mappings for degenerate symbols, whereas self.ambiguities
        contains mappings for both degenerate and non-degenerate symbols.
        Sometimes you want one, sometimes the other, so both are provided.
        """
        self._serialisable = {k: v for k, v in locals().items() if k != "self"}
        self.gap = gap
        self.missing = missing
        self.gaps = frozenset([gap, missing])
        if gaps:
            self.gaps = self.gaps.union(frozenset(gaps))
        self.label = label
        # set the sequence constructor
        if seq_constructor is None:
            seq_constructor = "".join  # safe default string constructor
        elif not preserve_existing_moltypes:
            seq_constructor.moltype = self
        self._make_seq = seq_constructor

        # set the ambiguities
        ambigs = {self.missing: tuple(motifset) + (self.gap,), self.gap: (self.gap,)}
        if ambiguities:
            ambigs.update(ambiguities)
        for c in motifset:
            ambigs[c] = (c,)
        self.ambiguities = ambigs

        # set complements -- must set before we make the alphabet group
        self.complements = complements or {}

        if make_alphabet_group:  # note: must use _original_ ambiguities here
            self.alphabets = AlphabetGroup(motifset, ambiguities, moltype=self)
            self.alphabet = self.alphabets.base
        else:
            if isinstance(motifset, Enumeration):
                self.alphabet = motifset
            elif max(len(motif) for motif in motifset) == 1:
                self.alphabet = CharAlphabet(motifset, moltype=self)
            else:
                self.alphabet = Alphabet(motifset, moltype=self)
        # set the other properties
        self.degenerates = ambiguities and ambiguities.copy() or {}
        self.degenerates[self.missing] = "".join(motifset) + self.gap
        self.matches = make_matches(motifset, self.gaps, self.degenerates)
        self.pairs = pairs and pairs.copy() or {}
        self.pairs.update(make_pairs(pairs, motifset, self.gaps, self.degenerates))
        self.mw_calculator = mw_calculator

        # add lowercase characters, if we're doing that
        if add_lower:
            self._add_lowercase()
        # cache various other data that make the calculations faster
        self._make_all()
        self._make_comp_table()
        # a gap can be a true gap char or a degenerate character, typically '?'
        # we therefore want to ensure consistent treatment across the definition
        # of characters as either gap or degenerate
        self.gap_string = "".join(self.gaps)
        strict_gap = "".join(set(self.gap_string) - set(self.degenerates))
        self.strip_degenerate = FunctionWrapper(
            KeepChars(strict_gap + "".join(self.alphabet))
        )
        self.strip_bad = FunctionWrapper(KeepChars("".join(self.All)))
        to_keep = set(self.alphabet) ^ set(self.degenerates) - set(self.gaps)
        self.strip_bad_and_gaps = FunctionWrapper(KeepChars("".join(to_keep)))

        # make inverse degenerates from degenerates
        # ensure that lowercase versions also exist if appropriate
        inv_degens = {}
        for key, val in list(self.degenerates.items()):
            inv_degens[frozenset(val)] = key.upper()
            if add_lower:
                inv_degens[frozenset("".join(val).lower())] = key.lower()
        for m in self.alphabet:
            inv_degens[frozenset(m)] = m
            if add_lower:
                inv_degens[frozenset("".join(m).lower())] = m.lower()
        for m in self.gaps:
            inv_degens[frozenset(m)] = m
        self.inverse_degenerates = inv_degens

        # set array type for modeling alphabets
        try:
            self.array_type = self.alphabet.array_type
        except AttributeError:
            self.array_type = None

        # set modeling sequence
        self._make_array_seq = array_seq_constructor

        self._colors = colors or defaultdict(_DefaultValue("black"))

    def __repr__(self):
        """String representation of MolType.

        WARNING: This doesn't allow you to reconstruct the object in its present
        incarnation.
        """
        return f"MolType({self.alphabet})"

    def __getnewargs_ex__(self, *args, **kw):
        data = self.to_rich_dict(for_pickle=True)
        return (), data

    def to_rich_dict(self, for_pickle=False):
        data = deepcopy(self._serialisable)
        if not for_pickle:  # we rely on reconstruction from label
            data = dict(type=get_object_provenance(self), moltype=self.label)
            data["version"] = __version__
        return data

    def to_json(self):
        """returns result of json formatted string"""
        data = self.to_rich_dict(for_pickle=False)
        return json.dumps(data)

    def to_regex(self, seq):
        """returns a regex pattern with ambiguities expanded to a character set"""
        if not self.is_valid(seq):
            raise ValueError(f"'{seq}' is invalid for this moltype")

        degen_indices = self.get_degenerate_positions(sequence=seq, include_gap=False)
        seq = list(seq)  # seq can now be modified
        for index in degen_indices:
            expanded = self.ambiguities[seq[index]]
            seq[index] = f"[{''.join(expanded)}]"
        return "".join(seq)

    def gettype(self):
        """Returns type, e.g. 'dna', 'rna', 'protein'. Delete?"""
        return self.label

    def make_seq(self, seq, name=None, **kwargs):
        """Returns sequence of correct type."""
        return self._make_seq(seq, name, **kwargs)

    def make_array_seq(self, seq, name=None, **kwargs):
        """
        creates an array sequence

        Parameters
        ----------
        seq
            characters or array
        name : str
        kwargs
            keyword arguments for the ArraySequence constructor.

        Returns
        -------
        ArraySequence
        """
        alphabet = kwargs.pop("alphabet", None)
        if alphabet is None and hasattr(self, "alphabets"):
            alphabet = self.alphabets.degen_gapped
        elif alphabet is None:
            alphabet = self.alphabet
        return self._make_array_seq(seq, alphabet=alphabet, name=name, **kwargs)

    def verify_sequence(self, seq, gaps_allowed=True, wildcards_allowed=True):
        """Checks whether sequence is valid on the default alphabet.

        Has special-case handling for gaps and wild-cards. This mechanism is
        probably useful to have in parallel with the validation routines that
        check specifically whether the sequence has gaps, degenerate symbols,
        etc., or that explicitly take an alphabet as input.
        """
        alpha = frozenset(self.ambiguities)
        if gaps_allowed:
            alpha = alpha.union(self.gaps)
        if wildcards_allowed:
            alpha = alpha.union(self.missing)
        try:
            nonalpha = re.compile(f"[^{re.escape(''.join(alpha))}]")
            badchar = nonalpha.search(seq)
            if badchar:
                motif = badchar.group()
                raise AlphabetError(motif)
        except TypeError:  # not alphabetic sequence: try slow method
            for motif in seq:
                if motif not in alpha:
                    raise AlphabetError(motif)

    def is_ambiguity(self, querymotif):
        """Return True if querymotif is an amibiguity character in alphabet.

        Parameters
        ----------
        querymotif
            the motif being queried.

        """

        return len(self.ambiguities[querymotif]) > 1

    def _what_ambiguity(self, motifs):
        """The code that represents all of 'motifs', and minimal others.

        Does this duplicate DegenerateFromSequence directly?
        """
        most_specific = len(self.alphabet) + 1
        result = self.missing
        for (code, motifs2) in list(self.ambiguities.items()):
            for c in motifs:
                if c not in motifs2:
                    break
            else:
                if len(motifs2) < most_specific:
                    most_specific = len(motifs2)
                    result = code
        return result

    def what_ambiguity(self, motifs):
        """The code that represents all of 'motifs', and minimal others.

        Does this duplicate DegenerateFromSequence directly?
        """
        if not hasattr(self, "_reverse_ambiguities"):
            self._reverse_ambiguities = {}
        motifs = frozenset(motifs)
        if motifs not in self._reverse_ambiguities:
            self._reverse_ambiguities[motifs] = self._what_ambiguity(motifs)
        return self._reverse_ambiguities[motifs]

    def _add_lowercase(self):
        """Adds lowercase versions of keys and vals to each internal dict."""
        for name in [
            "alphabet",
            "degenerates",
            "gaps",
            "complements",
            "pairs",
            "matches",
        ]:
            curr = getattr(self, name)
            # temp hack to get around re-ordering
            if isinstance(curr, Alphabet):
                curr = tuple(curr)
            new = add_lowercase(curr)
            setattr(self, name, new)

    def _make_all(self):
        """Sets self.All, which contains all the symbols self knows about.

        Note that the value of items in self.All will be the string containing
        the possibly degenerate set of symbols that the items expand to.
        """
        all = {}
        for i in self.alphabet:
            all[i] = i
        for key, val in list(self.degenerates.items()):
            all[key] = val
        for i in self.gaps:
            all[i] = i
        self.All = all

    def _make_comp_table(self):
        """Sets self.ComplementTable, which maps items onto their complements.

        Note: self.ComplementTable is only set if self.complements exists.
        """
        if self.complements:
            self.ComplementTable = maketrans(
                "".join(list(self.complements.keys())),
                "".join(list(self.complements.values())),
            )

    def complement(self, item):
        """Returns complement of item, using data from self.complements.

        Always tries to return same type as item: if item looks like a dict,
        will return list of keys.
        """
        if not self.complements:
            raise TypeError(
                "Tried to complement sequence using alphabet without complements."
            )
        try:
            return item.translate(self.ComplementTable)
        except (AttributeError, TypeError):
            item = iterable(item)
            get = self.complements.get
            return item.__class__([get(i, i) for i in item])

    def rc(self, item):
        """Returns reverse complement of item w/ data from self.complements.

        Always returns same type as input.
        """
        comp = list(self.complement(item))
        comp.reverse()
        if isinstance(item, str):
            return item.__class__("".join(comp))
        else:
            return item.__class__(comp)

    def strand_symmetric_motifs(self, motif_length=1):
        """returns ordered pairs of strand complementary motifs"""
        if not self.pairs:
            raise TypeError("moltype must be DNA or RNA")

        motif_set = self.alphabet.get_word_alphabet(word_length=motif_length)
        motif_pairs = []
        for m in motif_set:
            pair = tuple(sorted([m, self.complement(m)]))
            motif_pairs.append(pair)

        motif_pairs = set(motif_pairs)
        return motif_pairs

    def __contains__(self, item):
        """A MolType contains every character it knows about."""
        return item in self.All

    def __iter__(self):
        """A MolType iterates only over the characters in its Alphabet.."""
        return iter(self.alphabet)

    def is_gap(self, char):
        """Returns True if char is a gap."""
        return char in self.gaps

    def is_gapped(self, sequence):
        """Returns True if sequence contains gaps."""
        return self.first_gap(sequence) is not None

    def is_degenerate(self, sequence):
        """Returns True if sequence contains degenerate characters."""
        return self.first_degenerate(sequence) is not None

    def is_valid(self, sequence):
        """Returns True if sequence contains no items that are not in self."""
        try:
            return self.first_invalid(sequence) is None
        except:
            return False

    def is_strict(self, sequence):
        """Returns True if sequence contains only items in self.alphabet."""
        try:
            return (len(sequence) == 0) or (self.first_non_strict(sequence) is None)
        except:
            return False

    def valid_on_alphabet(self, sequence, alphabet=None):
        """Returns True if sequence contains only items in alphabet.

        alphabet can actually be anything that implements __contains__.
        Defaults to self.alphabet if not supplied.
        """
        if alphabet is None:
            alphabet = self.alphabet
        return first_index_in_set(sequence, alphabet) is not None

    def first_not_in_alphabet(self, sequence, alphabet=None):
        """Returns index of first item not in alphabet, or None.

        Defaults to self.alphabet if alphabet not supplied.
        """
        if alphabet is None:
            alphabet = self.alphabet
        return first_index_in_set(sequence, alphabet)

    def first_gap(self, sequence):
        """Returns the index of the first gap in the sequence, or None."""
        gap = self.gaps
        for i, s in enumerate(sequence):
            if s in gap:
                return i
        return None

    def first_degenerate(self, sequence):
        """Returns the index of first degenerate symbol in sequence, or None."""
        degen = self.degenerates
        for i, s in enumerate(sequence):
            if s in degen:
                return i
        return None

    def first_invalid(self, sequence):
        """Returns the index of first invalid symbol in sequence, or None."""
        all = self.All
        for i, s in enumerate(sequence):
            if s not in all:
                return i
        return None

    def first_non_strict(self, sequence):
        """Returns the index of first non-strict symbol in sequence, or None."""
        monomers = self.alphabet
        for i, s in enumerate(sequence):
            if s not in monomers:
                return i
        return None

    def disambiguate(self, sequence, method="strip"):
        """Returns a non-degenerate sequence from a degenerate one.

        method can be 'strip' (deletes any characters not in monomers or gaps)
        or 'random'(assigns the possibilities at random, using equal
        frequencies).
        """
        if method == "strip":
            try:
                return sequence.__class__(self.strip_degenerate(sequence))
            except:
                ambi = self.degenerates

                def not_ambiguous(x):
                    return x not in ambi

                return sequence.__class__(list(filter(not_ambiguous, sequence)))

        elif method == "random":
            degen = self.degenerates
            result = []
            for i in sequence:
                if i in degen:
                    result.append(choice(degen[i]))
                else:
                    result.append(i)
            if isinstance(sequence, str):
                return sequence.__class__("".join(result))
            else:
                return sequence.__class__(result)
        else:
            raise NotImplementedError(f"Got unknown method {method}")

    def degap(self, sequence):
        """Deletes all gap characters from sequence."""
        try:
            trans = dict([(i, None) for i in map(ord, self.gaps)])
            return sequence.__class__(sequence.translate(trans))
        except AttributeError:
            gap = self.gaps

            def not_gap(x):
                return x not in gap

            return sequence.__class__(list(filter(not_gap, sequence)))

    def gap_indices(self, sequence):
        """Returns list of indices of all gaps in the sequence, or []."""
        gaps = self.gaps
        return [i for i, s in enumerate(sequence) if s in gaps]

    def gap_vector(self, sequence):
        """Returns list of bool indicating gap or non-gap in sequence."""
        return list(map(self.is_gap, sequence))

    def gap_maps(self, sequence):
        """Returns tuple containing dicts mapping between gapped and ungapped.

        First element is a dict such that d[ungapped_coord] = gapped_coord.
        Second element is a dict such that d[gapped_coord] = ungapped_coord.

        Note that the dicts will be invalid if the sequence changes after the
        dicts are made.

        The gaps themselves are not in the dictionary, so use d.get() or test
        'if pos in d' to avoid KeyErrors if looking up all elements in a gapped
        sequence.
        """
        ungapped = {}
        gapped = {}
        num_gaps = 0
        for i, is_gap in enumerate(self.gap_vector(sequence)):
            if is_gap:
                num_gaps += 1
            else:
                ungapped[i] = i - num_gaps
                gapped[i - num_gaps] = i
        return gapped, ungapped

    def count_gaps(self, sequence):
        """Counts the gaps in the specified sequence."""
        gaps = self.gaps
        gap_count = sum(1 for s in sequence if s in gaps)
        return gap_count

    def get_degenerate_positions(self, sequence, include_gap=True):
        """returns indices matching degenerate characters"""
        degen = list(self.degenerates)
        if include_gap:
            degen.append(self.gap)

        return [i for i, c in enumerate(sequence) if c in degen]

    def count_degenerate(self, sequence):
        """Counts the degenerate bases in the specified sequence."""
        degen = self.degenerates
        degen_count = 0
        for s in sequence:
            if s in degen:
                degen_count += 1
        return degen_count

    def possibilities(self, sequence):
        """Counts number of possible sequences matching the sequence.

        Uses self.degenerates to decide how many possibilites there are at
        each position in the sequence.
        """
        degen = self.degenerates
        count = 1
        for s in sequence:
            if s in degen:
                count *= len(degen[s])
        return count

    def mw(self, sequence, method="random", delta=None):
        """Returns the molecular weight of the sequence.

        If the sequence is ambiguous, uses method (random or strip) to
        disambiguate the sequence.

        if delta is present, uses it instead of the standard weight adjustment.
        """
        if not sequence:
            return 0
        try:
            return self.mw_calculator(sequence, delta)
        except KeyError:  # assume sequence was ambiguous
            return self.mw_calculator(self.disambiguate(sequence, method), delta)

    def can_match(self, first, second):
        """Returns True if every pos in 1st could match same pos in 2nd.

        Truncates at length of shorter sequence.
        gaps are only allowed to match other gaps.
        """
        m = self.matches
        for pair in zip(first, second):
            if pair not in m:
                return False
        return True

    def can_mismatch(self, first, second):
        """Returns True if any position in 1st could cause a mismatch with 2nd.

        Truncates at length of shorter sequence.
        gaps are always counted as matches.
        """
        m = self.matches
        if not first or not second:
            return False

        for pair in zip(first, second):
            if not m.get(pair, None):
                return True
        return False

    def must_match(self, first, second):
        """Returns True if all positions in 1st must match positions in second."""
        return not self.can_mismatch(first, second)

    def can_pair(self, first, second):
        """Returns True if first and second could pair.

        Pairing occurs in reverse order, i.e. last position of second with
        first position of first, etc.

        Truncates at length of shorter sequence.
        gaps are only allowed to pair with other gaps, and are counted as 'weak'
        (same category as GU and degenerate pairs).

        NOTE: second must be able to be reverse
        """
        p = self.pairs
        sec = list(second)
        sec.reverse()
        for pair in zip(first, sec):
            if pair not in p:
                return False
        return True

    def can_mispair(self, first, second):
        """Returns True if any position in 1st could mispair with 2nd.

        Pairing occurs in reverse order, i.e. last position of second with
        first position of first, etc.

        Truncates at length of shorter sequence.
        gaps are always counted as possible mispairs, as are weak pairs like GU.
        """
        p = self.pairs
        if not first or not second:
            return False

        sec = list(second)
        sec.reverse()
        for pair in zip(first, sec):
            if not p.get(pair, None):
                return True
        return False

    def must_pair(self, first, second):
        """Returns True if all positions in 1st must pair with second.

        Pairing occurs in reverse order, i.e. last position of second with
        first position of first, etc.
        """
        return not self.can_mispair(first, second)

    def degenerate_from_seq(self, sequence):
        """Returns least degenerate symbol corresponding to chars in sequence.

        First tries to look up in self.inverse_degenerates. Then disambiguates
        and tries to look up in self.inverse_degenerates. Then tries converting
        the case (tries uppercase before lowercase). Raises TypeError if
        conversion fails.
        """
        symbols = frozenset(sequence)
        # check if symbols are already known
        inv_degens = self.inverse_degenerates
        result = inv_degens.get(symbols, None)
        if result:
            return result
        # then, try converting the symbols
        degens = self.All
        converted = set()
        for sym in symbols:
            for char in degens[sym]:
                converted.add(char)
        symbols = frozenset(converted)
        result = inv_degens.get(symbols, None)
        if result:
            return result
        # then, try converting case
        symbols = frozenset([s.upper() for s in symbols])
        result = inv_degens.get(symbols, None)
        if result:
            return result
        symbols = frozenset([s.lower() for s in symbols])
        result = inv_degens.get(symbols, None)
        if result:
            return result
        # finally, try to find the minimal subset containing the symbols
        symbols = frozenset([s.upper() for s in symbols])
        lengths = {}
        for i in inv_degens:
            if symbols.issubset(i):
                lengths[len(i)] = i
        if lengths:  # found at least some matches
            sorted = list(lengths.keys())
            sorted.sort()
            return inv_degens[lengths[sorted[0]]]

        # if we got here, nothing worked
        raise TypeError(f"Cannot find degenerate char for symbols: {symbols}")

    def get_css_style(self, colors=None, font_size=12, font_family="Lucida Console"):
        """returns string of CSS classes and {character: <CSS class name>, ...}

        Parameters
        ----------
        colors
            {char
        font_size
            in points
        font_family
            name of a monospace font

        """
        colors = colors or self._colors
        # !important required to stop some browsers over-riding the style sheet ...!!
        template = (
            '.%s_%s{font-family: "%s",monospace !important; '
            "font-size: %dpt !important; color: %s; }"
        )
        label = self.label or ""
        styles = _style_defaults[label].copy()
        styles.update(
            {c: "_".join([c, label]) for c in list(self.alphabet) + ["terminal_ambig"]}
        )

        css = [
            template % (char, label, font_family, font_size, colors[char])
            for char in list(styles) + ["ambig"]
        ]

        return css, styles


ASCII = MolType(
    # A default type for text read from a file etc. when we don't
    # want to prematurely assume DNA or Protein.
    seq_constructor=DefaultSequence,
    motifset=letters,
    ambiguities={},
    label="text",
    array_seq_constructor=ArraySequence,
)

DNA = MolType(
    seq_constructor=DnaSequence,
    motifset=IUPAC_DNA_chars,
    ambiguities=IUPAC_DNA_ambiguities,
    label="dna",
    mw_calculator=DnaMW,
    complements=IUPAC_DNA_ambiguities_complements,
    pairs=DnaStandardPairs,
    make_alphabet_group=True,
    array_seq_constructor=ArrayDnaSequence,
    colors=NT_COLORS,
)

RNA = MolType(
    seq_constructor=RnaSequence,
    motifset=IUPAC_RNA_chars,
    ambiguities=IUPAC_RNA_ambiguities,
    label="rna",
    mw_calculator=RnaMW,
    complements=IUPAC_RNA_ambiguities_complements,
    pairs=RnaStandardPairs,
    make_alphabet_group=True,
    array_seq_constructor=ArrayRnaSequence,
    colors=NT_COLORS,
)

PROTEIN = MolType(
    seq_constructor=ProteinSequence,
    motifset=IUPAC_PROTEIN_chars,
    ambiguities=IUPAC_PROTEIN_ambiguities,
    mw_calculator=ProteinMW,
    make_alphabet_group=True,
    array_seq_constructor=ArrayProteinSequence,
    label="protein",
    colors=AA_COLORS,
)

PROTEIN_WITH_STOP = MolType(
    seq_constructor=ProteinWithStopSequence,
    motifset=PROTEIN_WITH_STOP_chars,
    ambiguities=PROTEIN_WITH_STOP_ambiguities,
    mw_calculator=ProteinMW,
    make_alphabet_group=True,
    array_seq_constructor=ArrayProteinWithStopSequence,
    label="protein_with_stop",
    colors=AA_COLORS,
)

BYTES = MolType(
    # A default type for arbitrary chars read from a file etc. when we don't
    # want to prematurely assume _anything_ about the data.
    seq_constructor=ByteSequence,
    motifset=list(map(chr, list(range(256)))),
    ambiguities={},
    array_seq_constructor=ArraySequence,
    label="bytes",
)

# the None value catches cases where a moltype has no label attribute
_style_defaults = {
    getattr(mt, "label", ""): defaultdict(
        _DefaultValue(f"ambig_{getattr(mt, 'label', '')}")
    )
    for mt in (ASCII, BYTES, DNA, RNA, PROTEIN, PROTEIN_WITH_STOP, None)
}

# following is a two-state MolType useful for testing
AB = MolType(
    seq_constructor=ABSequence,
    motifset="ab",
    ambiguities={},
    array_seq_constructor=ArraySequence,
    label="ab",
)


class _CodonAlphabet(Alphabet):
    """Codon alphabets are DNA TupleAlphabets with a genetic code attribute and some codon-specific methods"""

    def _with(self, motifs):
        a = Alphabet._with(self, motifs)
        a.__class__ = type(self)
        a._gc = self._gc
        return a

    def is_sense_codon(self, codon):
        return not self._gc.is_stop(codon)

    def is_stop_codon(self, codon):
        return self._gc.is_stop(codon)

    def get_genetic_code(self):
        return self._gc


def CodonAlphabet(gc=1, include_stop_codons=False):
    if isinstance(gc, (int, str)):
        gc = get_code(gc)
    if include_stop_codons:
        motifset = list(gc.codons)
    else:
        motifset = list(gc.sense_codons)
    motifset = [codon.upper().replace("U", "T") for codon in motifset]
    a = _CodonAlphabet(motifset, moltype=DNA)
    a._gc = gc
    return a


def _method_codon_alphabet(ignore, *args, **kwargs):
    """If CodonAlphabet is set as a property, it gets self as extra 1st arg."""
    return CodonAlphabet(*args, **kwargs)


STANDARD_CODON = CodonAlphabet()

# Modify NucleicAcidSequence to avoid circular import
NucleicAcidSequence.codon_alphabet = _method_codon_alphabet
NucleicAcidSequence.protein = PROTEIN
ArrayRnaSequence.moltype = RNA
ArrayRnaSequence.alphabet = RNA.alphabets.degen_gapped

ArrayDnaSequence.moltype = DNA
ArrayDnaSequence.alphabet = DNA.alphabets.degen_gapped

ArrayProteinSequence.moltype = PROTEIN
ArrayProteinSequence.alphabet = PROTEIN.alphabets.degen_gapped

ArrayProteinWithStopSequence.moltype = PROTEIN_WITH_STOP
ArrayProteinWithStopSequence.alphabet = PROTEIN_WITH_STOP.alphabets.degen_gapped

ArraySequence.alphabet = BYTES.alphabet

ArrayAlignment.alphabet = BYTES.alphabet
ArrayAlignment.moltype = BYTES

ArrayDnaCodonSequence.alphabet = DNA.alphabets.base.Triples
ArrayRnaCodonSequence.alphabet = RNA.alphabets.base.Triples

# Modify Alignment to avoid circular import
Alignment.moltype = ASCII
SequenceCollection.moltype = BYTES


def _make_moltype_dict():
    env = globals()
    moltypes = {}
    for key in env:
        obj = env[key]
        if not isinstance(obj, MolType):
            continue
        if obj.label is not None:
            moltypes[obj.label] = obj

    return moltypes


moltypes = _make_moltype_dict()


def get_moltype(name):
    """returns the moltype with the matching name attribute"""
    if isinstance(name, MolType):
        return name
    name = name.lower()
    if name not in moltypes:
        raise ValueError(f"unknown moltype {name!r}")
    return moltypes[name]


def available_moltypes():
    """returns Table listing the available moltypes"""
    from cogent3.util.table import Table

    rows = []
    for n, m in moltypes.items():
        v = str(m)
        num = len(list(m))
        if num > 10:
            v = f"{v[:39]}..."
        rows.append([n, num, v])
    header = ["Abbreviation", "Number of states", "Moltype"]
    title = "Specify a moltype by the Abbreviation (case insensitive)."

    result = Table(header=header, data=rows, title=title, index_name="Abbreviation")
    result = result.sorted(columns=["Number of states", "Abbreviation"])
    result.format_column("Abbreviation", repr)
    return result
