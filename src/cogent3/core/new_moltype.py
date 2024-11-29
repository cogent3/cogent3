import dataclasses
import functools
import itertools
import typing
from collections import defaultdict
from string import ascii_letters

import numpy

from cogent3.core import new_alphabet, new_sequence
from cogent3.data.molecular_weight import DnaMW, ProteinMW, RnaMW, WeightCalculator

OptStr = typing.Optional[str]
OptFloat = typing.Optional[float]
OptCallable = typing.Optional[typing.Callable]
SeqStrType = typing.Union[list[str], tuple[str, ...]]
StrORBytes = typing.Union[str, bytes]
StrORArray = typing.Union[str, numpy.ndarray]
StrORBytesORArray = typing.Union[str, bytes, numpy.ndarray]
SeqStrBytesType = typing.Union[list[StrORBytes], tuple[StrORBytes, ...]]

IUPAC_gap = "-"

IUPAC_missing = "?"

IUPAC_DNA_chars = "T", "C", "A", "G"
IUPAC_DNA_ambiguities = {
    "N": frozenset(("A", "C", "T", "G")),
    "R": frozenset(("A", "G")),
    "Y": frozenset(("C", "T")),
    "W": frozenset(("A", "T")),
    "S": frozenset(("C", "G")),
    "K": frozenset(("T", "G")),
    "M": frozenset(("C", "A")),
    "B": frozenset(("C", "T", "G")),
    "D": frozenset(("A", "T", "G")),
    "H": frozenset(("A", "C", "T")),
    "V": frozenset(("A", "C", "G")),
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
    "?": "?",
}

IUPAC_DNA_complements = {"A": "T", "C": "G", "G": "C", "T": "A", "-": "-"}
# Standard DNA pairing: only Watson-Crick pairs count as pairs
DNA_STANDARD_PAIRS = {
    frozenset(("A", "T")): True,
    frozenset(("C", "G")): True,
}

# note change in standard order from DNA
IUPAC_RNA_chars = "U", "C", "A", "G"
IUPAC_RNA_ambiguities = {
    "N": frozenset(("A", "C", "U", "G")),
    "R": frozenset(("A", "G")),
    "Y": frozenset(("C", "U")),
    "W": frozenset(("A", "U")),
    "S": frozenset(("C", "G")),
    "K": frozenset(("U", "G")),
    "M": frozenset(("C", "A")),
    "B": frozenset(("C", "U", "G")),
    "D": frozenset(("A", "U", "G")),
    "H": frozenset(("A", "C", "U")),
    "V": frozenset(("A", "C", "G")),
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
    "?": "?",
}

IUPAC_RNA_complements = {"A": "U", "C": "G", "G": "C", "U": "A", "-": "-"}

# Standard RNA pairing: GU pairs count as 'weak' pairs
RNA_STANDARD_PAIRS = {
    frozenset(("A", "U")): True,  # True vs False for 'always' vs 'sometimes' pairing
    frozenset(("C", "G")): True,
    frozenset(("G", "U")): False,
}

# Watson-Crick RNA pairing only: GU pairs don't count as pairs
RNA_W_C_PAIRS = {
    frozenset(("A", "U")): True,
    frozenset(("C", "G")): True,
    frozenset(("U", "A")): True,
}

# RNA pairing with GU counted as standard pairs
RNA_G_U_PAIRS = {
    frozenset(("A", "U")): True,
    frozenset(("C", "G")): True,
    frozenset(("G", "C")): True,
    frozenset(("U", "A")): True,
    frozenset(("G", "U")): True,
    frozenset(("U", "G")): True,
}

# RNA pairing with GU, AA, GA, CA and UU mismatches allowed as weak pairs
RNA_EXTENDED_PAIRS = {
    frozenset({"A", "U"}): True,
    frozenset({"C", "G"}): True,
    frozenset({"G", "U"}): False,
    frozenset({"A"}): False,
    frozenset({"A", "G"}): False,
    frozenset({"A", "C"}): False,
    frozenset({"U"}): False,
}


class MolTypeError(TypeError): ...


def make_pairs(
    *,
    pairs: dict[tuple[str, str], bool] = None,
    monomers: tuple[str, ...] = None,
    gaps: str = None,
    degenerates: dict[str, set[str]] = None,
) -> dict[frozenset[str], bool]:
    """Makes a dict of symbol pairs (i,j) -> strictness.

    Expands pairs into all possible pairs using degen symbols.
    Strictness is True if i and j always pair, and False if they 'weakly' pair
    (e.g. GU pairs or if it is possible that they pair).

    If you want to make GU pairs count as 'always matching', pass in pairs
    that have (G,U) and (U, G) mapped to True rather than False.
    """
    result = {}
    pairs = pairs or {}
    monomers = monomers or ()
    gaps = gaps or ()
    degenerates = degenerates or {}

    result |= pairs
    result |= {frozenset((i, j)): False for i in gaps for j in gaps}
    for b, d in itertools.product(monomers, degenerates):
        if any(frozenset((b, e)) in pairs for e in degenerates[d]):
            result[frozenset((b, d))] = False

    for d1, d2 in itertools.combinations_with_replacement(degenerates, 2):
        if any(
            frozenset((e1, e2)) in pairs
            for e1 in degenerates[d1]
            for e2 in degenerates[d2]
        ):
            result[frozenset((d1, d2))] = False

    return result


def make_matches(
    monomers: tuple[str, ...] = None,
    gaps: str = None,
    degenerates: dict[str, set[str]] = None,
) -> dict[frozenset[str], bool]:
    """Makes a dict of symbol pairs (i,j) -> strictness.

    Strictness is True if i and j always match and False if they sometimes
    match (e.g. A always matches A, but W sometimes matches R).
    """

    monomers = monomers or {}
    gaps = gaps or {}
    degenerates = degenerates or {}
    result = {}

    # all monomers always match themselves and no other monomers
    for i in monomers:
        result[(i, i)] = True

    # all gaps always match all other gaps
    for i, j in itertools.combinations_with_replacement(gaps, 2):
        result[(i, j)], result[(j, i)] = True, True

    # monomers sometimes match degenerates that contain them
    for i in monomers:
        for j in degenerates:
            if i in degenerates[j]:
                result[(i, j)], result[(j, i)] = False, False

    # degenerates sometimes match degenerates if the sets intersect
    for i in degenerates:
        for j in degenerates:
            if degenerates[i] & degenerates[j]:
                result[(i, j)], result[(j, i)] = False, False
    return result


# RNA_PAIRING_RULES is a dict of {name:(base_pairs,degen_pairs)} where base_pairs
# is a dict with the non-degenerate pairing rules and degen_pairs is a dict with
# both the degenerate and non-degenerate pairing rules.
# NOTE: uses make_pairs to augment the initial dict after construction.
def _build_pairing_rules() -> dict[frozenset[str], bool]:
    pairing_rules = {
        "Standard": RNA_STANDARD_PAIRS,
        "WC": RNA_W_C_PAIRS,
        "GU": RNA_G_U_PAIRS,
        "Extended": RNA_EXTENDED_PAIRS,
    }
    for k, v in list(pairing_rules.items()):
        pairing_rules[k] = (v, make_pairs(pairs=v))
    return pairing_rules


RNA_PAIRING_RULES = _build_pairing_rules()

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

IUPAC_PROTEIN_chars = (
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
)

PROTEIN_WITH_STOP_chars = (
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
)

IUPAC_PROTEIN_ambiguities = {
    "B": frozenset(["N", "D"]),
    "X": frozenset(IUPAC_PROTEIN_chars),
    "Z": frozenset(["Q", "E"]),
}

PROTEIN_WITH_STOP_ambiguities = {
    "B": frozenset(["N", "D"]),
    "X": frozenset(PROTEIN_WITH_STOP_chars),
    "Z": frozenset(["Q", "E"]),
}

# styling for moltype display


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


def _strictly_upper(monomers: tuple[StrORBytes]):
    """whether all elements correspond to upper case characters"""
    cast = str if isinstance(monomers[0], str) else chr
    return all(cast(c).isupper() for c in monomers)


@dataclasses.dataclass
class MolType:
    """MolType handles operations that depend on the sequence type.

    Notes
    -----
    The only way to create sequences is via a MolType instance.
    The instance defines different alphabets that are used for data
    conversions.
    Create a moltype using the ``get_moltype()`` function.
    """

    name: str
    monomers: dataclasses.InitVar[StrORBytes]
    make_seq: dataclasses.InitVar[typing.Type]
    gap: OptStr = IUPAC_gap
    missing: OptStr = IUPAC_missing
    complements: dataclasses.InitVar[typing.Optional[dict[str, str]]] = None
    ambiguities: typing.Optional[dict[str, tuple[str, ...]]] = None
    colors: dataclasses.InitVar[typing.Optional[dict[str, str]]] = None
    pairing_rules: typing.Optional[dict[str, dict[frozenset[str], bool]]] = None
    mw_calculator: typing.Optional[WeightCalculator] = None

    # private attributes to be delivered via properties
    _monomers: new_alphabet.CharAlphabet = dataclasses.field(init=False)
    _gapped: new_alphabet.CharAlphabet = dataclasses.field(init=False)
    _gapped_missing: new_alphabet.CharAlphabet = dataclasses.field(init=False)
    _degen: new_alphabet.CharAlphabet = dataclasses.field(init=False)
    _degen_gapped: new_alphabet.CharAlphabet = dataclasses.field(init=False)
    _colors: dict[str, str] = dataclasses.field(init=False)

    # how to connect this to the sequence constructor and avoid
    # circular imports
    _make_seq: typing.Callable = dataclasses.field(init=False)
    _complement: OptCallable = dataclasses.field(init=False, default=None)

    def __post_init__(
        self,
        monomers: StrORBytes,
        make_seq: typing.Type,
        complements: typing.Optional[dict[str, str]],
        colors: typing.Optional[dict[str, str]],
    ):
        self._colors = colors or defaultdict(_DefaultValue("black"))
        self._make_seq = make_seq
        gap = new_alphabet._coerce_to_type(monomers, self.gap or "")
        missing = new_alphabet._coerce_to_type(monomers, self.missing or "")
        ambigs = new_alphabet._coerce_to_type(monomers, "".join(self.ambiguities or ""))

        self._monomers = new_alphabet.make_alphabet(
            chars=monomers, gap=None, missing=None, moltype=self
        )
        self._degen = (
            new_alphabet.make_alphabet(
                chars=monomers + ambigs + missing,
                gap=None,
                missing=missing or None,
                moltype=self,
            )
            if ambigs
            else None
        )
        self._gapped = (
            new_alphabet.make_alphabet(
                chars=monomers + gap, gap=self.gap, missing=None, moltype=self
            )
            if gap
            else None
        )
        self._gapped_missing = (
            new_alphabet.make_alphabet(
                chars=monomers + gap + missing,
                gap=self.gap,
                missing=missing or None,
                moltype=self,
            )
            if missing and gap
            else None
        )
        self._degen_gapped = (
            new_alphabet.make_alphabet(
                chars=monomers + gap + ambigs + missing,
                gap=self.gap,
                missing=missing or None,
                moltype=self,
            )
            if ambigs and gap
            else None
        )

        if complements:
            # assume we have a nucleic acid moltype
            dest = "".join(complements[c] for c in self.degen_gapped_alphabet)
            self._complement = new_alphabet.convert_alphabet(
                self.degen_gapped_alphabet.as_bytes(),
                dest.encode("utf8"),
            )

    def __repr__(self):
        name = self.__class__.__name__
        return f"{name}({self.alphabet})"

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return id(self) == id(other)

    def __len__(self) -> int:
        return len(self._monomers)

    def __iter__(self):
        yield from self._monomers

    @property
    def label(self):
        """synonym for name"""
        return self.name

    @property
    def alphabet(self):
        """monomers"""
        return self._monomers

    @property
    def degen_alphabet(self):
        """monomers + ambiguous characters"""
        return self._degen

    @property
    def gapped_alphabet(self):
        """monomers + gap"""
        return self._gapped

    @property
    def degen_gapped_alphabet(self):
        """monomers + gap + ambiguous characters"""
        return self._degen_gapped

    @property
    def gaps(self) -> frozenset:  # refactor: docstring
        gaps = [char for char in (self.gap, self.missing) if char is not None]
        return frozenset(gaps)

    @property
    def matching_rules(self) -> dict[frozenset[str], bool]:
        return make_matches(
            monomers=self._monomers, gaps=self.gaps, degenerates=self.ambiguities
        )

    def is_valid(self, seq: StrORArray) -> bool:
        """checks against most degenerate alphabet"""
        alpha = self.most_degen_alphabet()
        return alpha.is_valid(seq)

    def iter_alphabets(self):
        """yield alphabets in order of most to least degenerate"""
        alphas = (
            self._degen_gapped,
            self._degen,
            self._gapped_missing,
            self._gapped,
            self._monomers,
        )
        yield from (a for a in alphas if a)

    def is_compatible_alphabet(
        self, alphabet: new_alphabet.CharAlphabet, strict: bool = True
    ) -> bool:
        """checks that characters in alphabet are equal to a bound alphabet

        Parameters
        ----------
        alphabet
            an Alphabet instance
        strict
            the order of elements must match
        """
        if not strict:
            query = set(alphabet)
            return any(set(alpha) == query for alpha in self.iter_alphabets())

        return any(alpha == alphabet for alpha in self.iter_alphabets())

    def make_seq(
        self,
        *,
        seq: typing.Union[str, "SeqViewABC"],
        name: OptStr = None,
        check_seq: bool = True,
        **kwargs,
    ) -> new_sequence.Sequence:
        """creates a Sequence object corresponding to the molecular type of
        this instance.

        Parameters
        ----------
        seq
            the raw sequence data
        name
            the name of the sequence
        check_seq
            whether to validate the sequence data against the molecular type
            if True, performs validation checks; set to False if the sequence
            data is already known to be valid.
        **kwargs
            additional keyword arguments that may be required for creating the
            Sequence object
        """
        if hasattr(seq, "moltype"):
            return seq if seq.moltype is self else seq.to_moltype(self)

        if check_seq:
            assert self.is_valid(
                seq
            ), f"{seq[:4]!r} not valid for moltype {self.name!r}"
        seq = "" if seq is None else seq
        return self._make_seq(moltype=self, seq=seq, name=name, **kwargs)

    @functools.singledispatchmethod
    def complement(self, seq: StrORBytesORArray, validate: bool = True) -> str:
        """converts a string or bytes into it's nucleic acid complement

        Parameters
        ----------
        seq
            sequence to be complemented
        validate
            if True, checks the sequence is validated against the most
            degenerate alphabet
        """
        raise TypeError(f"{type(seq)} not supported")

    @complement.register
    def _(self, seq: str, validate: bool = True) -> str:
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )

        return self.complement(seq.encode("utf8"), validate=False).decode("utf8")

    @complement.register
    def _(self, seq: bytes, validate: bool = True) -> str:
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )
        if not self.is_nucleic:
            raise MolTypeError(f"{self.name!r} cannot complement")
        return self._complement(seq)

    @complement.register
    def _(self, seq: numpy.ndarray, validate: bool = True) -> numpy.ndarray[int]:
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )

        return self.degen_gapped_alphabet.to_indices(
            self.complement(
                self.degen_gapped_alphabet.array_to_bytes(seq), validate=False
            )
        )

    @property
    def is_nucleic(self) -> bool:
        """is a nucleic acid moltype

        Notes
        -----
        nucleic moltypes can be used for complementing and translating
        into amino acids.
        """
        return callable(self._complement)

    def rc(self, seq: str, validate: bool = True) -> str:
        """reverse reverse complement of a sequence"""
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )

        return self.complement(seq)[::-1]

    @functools.singledispatchmethod
    def is_degenerate(self, seq: StrORBytesORArray, validate: bool = True) -> bool:
        """checks if a sequence contains degenerate characters"""
        raise TypeError(f"{type(seq)} not supported")

    @is_degenerate.register
    def _(self, seq: bytes, validate: bool = True) -> bool:  # refactor: docstring
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )
        return self.is_degenerate(
            self.most_degen_alphabet().to_indices(seq), validate=False
        )

    @is_degenerate.register
    def _(self, seq: str, validate: bool = True) -> bool:
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )
        return self.is_degenerate(
            self.most_degen_alphabet().to_indices(seq), validate=False
        )

    @is_degenerate.register
    def _(self, seq: numpy.ndarray, validate: bool = True) -> bool:
        if self.degen_alphabet is None:
            return False

        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )

        first_degen = self.most_degen_alphabet().gap_index + 1

        return (seq >= first_degen).any()

    @functools.singledispatchmethod
    def is_gapped(self, seq, validate: bool = True) -> bool:
        """checks if a sequence contains gaps"""
        raise TypeError(f"{type(seq)} not supported")

    @is_gapped.register
    def _(self, seq: str, validate: bool = True) -> bool:
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )
        return any(gap in seq for gap in self.gaps)

    @is_gapped.register
    def _(self, seq: bytes, validate: bool = True) -> bool:
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )
        return self.is_gapped(seq.decode("utf8"), validate=False)

    @is_gapped.register
    def _(self, seq: bytes, validate: bool = True) -> bool:
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )
        return self.is_gapped(seq.decode("utf8"), validate=False)

    @is_gapped.register
    def _(self, seq: numpy.ndarray, validate: bool = True) -> bool:
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )
        alpha = self.most_degen_alphabet()
        gaps_indices = [alpha.index(gap) for gap in self.gaps]
        return numpy.isin(seq, gaps_indices).any()

    def get_degenerate_positions(
        self, seq: StrORBytesORArray, include_gap: bool = True, validate: bool = True
    ) -> list[int]:  # refactor: docstring
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )

        if self.degen_alphabet is None:
            return []

        alpha = self.most_degen_alphabet()
        seq = alpha.to_indices(seq)
        index = len(self) if include_gap else len(self) + 1
        degens = seq >= index
        return numpy.where(degens)[0].tolist()

    @functools.singledispatchmethod
    def degap(self, seq, validate: bool = True) -> StrORBytesORArray:
        """removes all gap and missing characters from a sequence"""
        # refactor: design
        # cache translation callables (using alpabet.convert_alphabet)
        # previous implementation had support for tuples -- is this necessary?
        raise TypeError(f"{type(seq)} not supported")

    @degap.register
    def _(self, seq: bytes, validate: bool = True) -> bytes:
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )

        if self.gap is None:
            raise TypeError(f"no gap character defined for {self.name!r}")

        return seq.translate(None, delete="".join(self.gaps).encode("utf8"))

    @degap.register
    def _(self, seq: str, validate: bool = True) -> str:
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )
        degapped = self.degap(seq.encode("utf8"), validate=False)
        return degapped.decode("utf8")

    @degap.register
    def _(self, seq: numpy.ndarray, validate: bool = True) -> numpy.ndarray:
        if validate and not self.is_valid(seq):
            raise new_alphabet.AlphabetError(
                f"{seq[:4]!r} not valid for moltype {self.name!r}"
            )

        degen = self.most_degen_alphabet()
        if not degen.gap_char:
            raise TypeError(f"no gap character defined for {self.name!r}")

        degapped = degen.array_to_bytes(seq)
        return degen.to_indices(self.degap(degapped, validate=False))

    def is_ambiguity(self, query_motif: str, validate: bool = True) -> bool:
        """Return True if querymotif is an amibiguity character in alphabet.

        Parameters
        ----------
        query_motif
            the motif being queried.

        """
        if validate and not self.is_valid(query_motif):
            raise new_alphabet.AlphabetError(
                f"{query_motif[:4]!r} not valid for moltype {self.name!r}"
            )

        ambigs_missing = {
            self.missing,
            *frozenset(self.ambiguities.keys()),
        }
        return query_motif in ambigs_missing

    def resolve_ambiguity(
        self,
        ambig_motif: str,
        alphabet: typing.Optional[new_alphabet.CharAlphabet] = None,
        allow_gap: bool = False,
        validate: bool = True,
    ) -> typing.Tuple[str]:
        """Returns tuple of all possible canonical characters corresponding
        to ambig_motif

        Parameters
        ----------
        ambig_motif
            the string to be expanded
        alphabet
            optional, disambiguated motifs not present in alphabet will be
            excluded. This could be a codon alphabet where stop codons are
            not present.
        allow_gap
            whether the gap character is allowed in output. Only
            applied when alphabet is None.

        Notes
        -----
        If ambig_motif is > 1 character long and alphabet is None, we construct
        a word alphabet with the same length.
        """
        if validate and not self.is_valid(ambig_motif):
            raise new_alphabet.AlphabetError(
                f"{ambig_motif[:4]!r} not valid for moltype {self.name!r}"
            )

        # refactor: simplify
        ambiguities = {
            **self.ambiguities,
            self.gap: frozenset(self.gap),
            self.missing: set([*self.alphabet, self.gap]),
        }
        for m in self.alphabet:
            ambiguities[m] = frozenset(m)
        if alphabet is None:
            word_alpha = self.alphabet.get_kmer_alphabet(
                k=len(ambig_motif), include_gap=allow_gap
            )
            # refactor: design
            # with_gap_motif() is not implemented on KmerAlphabet, check whether
            # allow_gap makes sense in this case
            alphabet = word_alpha.with_gap_motif() if allow_gap else word_alpha
            if not allow_gap:
                ambiguities["?"] = tuple(c for c in ambiguities["?"] if c != self.gap)

        if ambig_motif in alphabet:
            return (ambig_motif,)

        try:
            resolved = [ambiguities[c] for c in ambig_motif]
        except KeyError as e:
            raise new_alphabet.AlphabetError(ambig_motif) from e

        result = tuple("".join(e) for e in itertools.product(*resolved))
        if alphabet:
            result = tuple(e for e in result if e in alphabet)

        if not result:
            raise new_alphabet.AlphabetError(ambig_motif)

        return result

    @functools.singledispatchmethod
    def strip_degenerate(self, seq: StrORBytesORArray) -> StrORBytesORArray:
        """removes degenerate characters"""
        raise TypeError(f"{type(seq)} not supported")

    @strip_degenerate.register
    def _(self, seq: bytes) -> bytes:
        monomers = self.alphabet
        ambigs = set(self.degen_alphabet) - set(monomers)
        func = new_alphabet.convert_alphabet(
            src=monomers.as_bytes(),
            dest=monomers.as_bytes(),
            delete="".join(ambigs).encode("utf-8"),
        )
        return func(seq)

    @strip_degenerate.register
    def _(self, seq: str) -> str:
        return self.strip_degenerate(seq.encode("utf-8")).decode("utf-8")

    @strip_degenerate.register
    def _(self, seq: numpy.ndarray) -> numpy.ndarray:
        return self.degen_gapped_alphabet.to_indices(
            self.strip_degenerate(self.degen_gapped_alphabet.array_to_bytes(seq))
        )

    @functools.singledispatchmethod
    def strip_bad(self, seq: StrORBytes) -> StrORBytes:
        """Removes any symbols not in the alphabet."""
        raise TypeError(f"{type(seq)} not supported")

    def _strip_bad(self, seq: numpy.ndarray) -> numpy.ndarray:
        return seq[seq < len(self.degen_gapped_alphabet)]

    @strip_bad.register
    def _(self, seq: bytes) -> bytes:
        return self.degen_gapped_alphabet.array_to_bytes(
            self._strip_bad(self.degen_gapped_alphabet.to_indices(seq))
        )

    @strip_bad.register
    def _(self, seq: str) -> str:
        return self.strip_bad(seq.encode("utf8")).decode("utf8")

    @functools.singledispatchmethod
    def strip_bad_and_gaps(self, seq: StrORBytes) -> StrORBytes:
        """Removes any symbols not in the alphabet, and any gaps.
        Since missing could be a gap, it is also removed."""
        raise TypeError(f"{type(seq)} not supported")

    def _strip_bad_and_gaps(self, seq: numpy.ndarray) -> numpy.ndarray:
        return seq[
            numpy.logical_and(
                numpy.logical_and(
                    seq < len(self.degen_gapped_alphabet),
                    seq != self.degen_gapped_alphabet.gap_index,
                ),
                seq != self.degen_gapped_alphabet.missing_index,
            )
        ]

    @strip_bad_and_gaps.register
    def _(self, seq: bytes) -> bytes:
        return self.degen_gapped_alphabet.array_to_bytes(
            self._strip_bad_and_gaps(self.degen_gapped_alphabet.to_indices(seq))
        )

    @strip_bad_and_gaps.register
    def _(self, seq: str) -> str:
        return self.strip_bad_and_gaps(seq.encode("utf8")).decode("utf8")

    def disambiguate(
        self, seq: StrORBytesORArray, method: str = "strip"
    ) -> StrORBytesORArray:
        """Returns a non-degenerate sequence from a degenerate one.

        Parameters
        ----------
        seq
            the sequence to be disambiguated
        method
            how to disambiguate the sequence, one of "strip", "random"
            strip: removes degenerate characters
            random: randomly selects a non-degenerate character
        """
        if method == "strip":
            return self.strip_degenerate(seq)

        elif method == "random":
            return self.random_disambiguate(seq)

        else:
            raise NotImplementedError(f"method={method} not implemented")

    @functools.singledispatchmethod
    def random_disambiguate(self, seq: StrORBytesORArray) -> StrORBytesORArray:
        """disambiguates a sequence by randomly selecting a non-degenerate character"""
        raise TypeError(f"{type(seq)} not supported")

    @random_disambiguate.register
    def _(self, seq: str) -> str:
        return "".join(
            numpy.random.choice(list(self.ambiguities.get(n, n))) for n in seq
        )

    @random_disambiguate.register
    def _(self, seq: bytes) -> bytes:
        return self.random_disambiguate(seq.decode("utf8")).encode("utf8")

    @random_disambiguate.register
    def _(self, seq: numpy.ndarray) -> numpy.ndarray:
        return self.degen_gapped_alphabet.to_indices(
            self.random_disambiguate(self.degen_gapped_alphabet.array_to_bytes(seq))
        )

    @functools.singledispatchmethod
    def count_gaps(self, seq: StrORBytes) -> int:
        """returns the number of gap characters in a sequence"""
        raise TypeError(f"{type(seq)} not supported")

    def _count_gaps(self, seq: numpy.ndarray) -> int:
        return numpy.sum(seq == self.degen_gapped_alphabet.gap_index)

    @count_gaps.register
    def _(self, seq: bytes) -> int:
        return self._count_gaps(self.degen_gapped_alphabet.to_indices(seq))

    @count_gaps.register
    def _(self, seq: str) -> int:
        return self.count_gaps(seq.encode("utf8"))

    @functools.singledispatchmethod
    def count_degenerate(self, seq: StrORBytes) -> int:
        """returns the number of degenerate characters in a sequence"""
        raise TypeError(f"{type(seq)} not supported")

    def _count_degenerate(self, seq: numpy.ndarray) -> int:
        return numpy.sum(seq >= self.degen_gapped_alphabet.gap_index)

    @count_degenerate.register
    def _(self, seq: bytes) -> int:
        return self._count_degenerate(self.degen_gapped_alphabet.to_indices(seq))

    @count_degenerate.register
    def _(self, seq: str) -> int:
        return self.count_degenerate(seq.encode("utf8"))

    @functools.singledispatchmethod
    def count_variants(self, seq: StrORBytes) -> int:
        """Counts number of possible sequences matching the sequence, given
        any ambiguous characters in the sequence.

        Notes
        -----
        Uses self.ambiguitues to decide how many possibilities there are at
        each position in the sequence and calculates the permutations.
        """
        raise TypeError(f"{type(seq)} not supported")

    @count_variants.register
    def _(self, seq: bytes) -> int:
        return self.count_variants(seq.decode("utf8"))

    @count_variants.register
    def _(self, seq: str) -> int:
        return numpy.prod([len(self.ambiguities.get(c, c)) for c in seq])

    def degenerate_from_seq(self, seq: str) -> str:
        """Returns least degenerate symbol that encompasses a set of characters"""
        symbols = frozenset(seq)

        # if the length of the set is 1, return the single character
        if len(symbols) == 1:
            return next(iter(symbols))

        # all degenerate symbols should be added to the sets that they encompass
        degens = self.ambiguities.copy()
        for degen1 in degens:
            for degen2, char_set in self.ambiguities.items():
                if self.ambiguities[degen1] <= char_set:
                    degens[degen2] = degens[degen2].union(frozenset(degen1))

        # reverse the mapping between degenerate symbols and their encompassing sets
        inv_degens = {frozenset(val): key for key, val in list(degens.items())}

        # add gap and missing characters to the mapping
        inv_degens[frozenset(self.gap)] = self.gap
        inv_degens[frozenset(self.degen_gapped_alphabet)] = self.missing

        # if we exactly match a set of symbols, return the corresponding degenerate
        if result := inv_degens.get(symbols):
            return result

        # if we don't exactly match, sort all encompassing sets and return the
        # least degenerate one
        encompassing = [chars for chars in inv_degens if chars.issuperset(symbols)]
        encompassing = sorted(encompassing, key=len)
        return inv_degens[encompassing[0]]

    def strand_symmetric_motifs(
        self, motif_length: int = 1
    ) -> set[tuple[str, str]]:  # refactor: docstring
        """returns ordered pairs of strand complementary motifs"""
        if len(self.alphabet) != 4:
            raise TypeError("moltype must be DNA or RNA")

        motif_set = self.alphabet.get_kmer_alphabet(k=motif_length)
        return {tuple(sorted([m, self.complement(m)])) for m in motif_set}

    def mw(self, seq: str, method: str = "random", delta: OptFloat = None) -> float:
        """Returns the molecular weight of the sequence. If the sequence is
        ambiguous, uses method to disambiguate the sequence.

        Parameters
        ----------
        seq
            the sequence whose molecular weight is to be calculated.
        method
            the method provided to .disambiguate() to disambiguate the sequence.
            either "random" or "first".
        delta
            if delta is present, uses it instead of the standard weight adjustment.

        """
        if not seq:
            return 0
        try:
            return self.mw_calculator(seq, delta)
        except KeyError:  # assume sequence was ambiguous
            return self.mw_calculator(self.disambiguate(seq, method), delta)

    def can_match(self, first: str, second: str) -> bool:
        """Returns True if every pos in 1st could match same pos in 2nd.

        Notes
        -----
        Truncates at length of shorter sequence. gaps are only allowed to match
        other gaps.
        """
        # refactor: design
        # is this method necessary?
        m = self.matching_rules
        return all(pair in m for pair in zip(first, second))

    def can_pair(self, first: str, second: str) -> bool:
        pairs = make_pairs(
            pairs=self.pairing_rules,
            monomers=self._monomers,
            gaps=self.gaps,
            degenerates=self.ambiguities,
        )
        sec = list(second)
        sec.reverse()
        for pair in zip(first, sec):
            if frozenset(pair) not in pairs:
                return False
        return True

    def can_mispair(self, first: str, second: str) -> bool:
        """Returns True if any position in self could mispair with other.

        Notes
        -----
        Pairing occurs in reverse order, i.e. last position of other with
        first position of self, etc.

        Truncates at length of shorter sequence.

        Gaps are always counted as possible mispairs, as are weak pairs like GU.
        """
        p = self.pairing_rules
        if not first or not second:
            return False

        for pair in zip(first, second[::-1]):
            if not p.get(frozenset(pair), None):
                return True
        return False

    def get_css_style(
        self,
        colors: typing.Optional[dict[str, str]] = None,
        font_size: int = 12,
        font_family="Lucida Console",
    ):
        """returns string of CSS classes and {character: <CSS class name>, ...}

        Parameters
        ----------
        colors
             A dictionary mapping characters to CSS color values.
        font_size
            Font size in points.
        font_family
            Name of a monospace font.

        """
        colors = colors or self._colors
        # !important required to stop some browsers over-riding the style sheet ...!!
        template = (
            '.%s_%s{font-family: "%s",monospace !important; '
            "font-size: %dpt !important; color: %s; }"
        )
        label = self.label or ""
        styles = _STYLE_DEFAULTS[label].copy()
        styles.update(
            {c: "_".join([c, label]) for c in list(self.alphabet) + ["terminal_ambig"]}
        )

        css = [
            template % (char, label, font_family, font_size, colors[char])
            for char in list(styles) + ["ambig"]
        ]

        return css, styles

    @functools.cache
    def most_degen_alphabet(self):
        """returns the most degenerate alphabet for this instance"""
        return next(self.iter_alphabets())

    def to_regex(self, seq: str) -> str:
        """returns a regex pattern with ambiguities expanded to a character set"""
        if not self.is_valid(seq):
            raise ValueError(f"'{seq}' is invalid for this moltype")

        degen_indices = self.get_degenerate_positions(seq=seq, include_gap=False)
        seq = list(seq)  # seq can now be modified
        for index in degen_indices:
            expanded = self.ambiguities[seq[index]]
            seq[index] = f"[{''.join(sorted(expanded))}]"
        return "".join(seq)


def _make_moltype_dict() -> dict[str, MolType]:
    """make a dictionary of local name space molecular types"""
    env = globals()
    moltypes = {}
    for obj in env.values():
        if not isinstance(obj, MolType):
            continue
        moltypes[obj.name] = obj

    return moltypes


def get_moltype(name: typing.Union[str, MolType]) -> MolType:
    """returns the moltype with the matching name attribute"""
    if isinstance(name, MolType):
        return name
    # refactor: simplify
    # intoduced to check whether we have a old-style moltype object,
    # remove when support for the old-style moltype objects is removed
    if hasattr(name, "label"):
        name = name.label
    name = name.lower()
    if name not in _moltypes:
        raise ValueError(f"unknown moltype {name!r}")
    return _moltypes[name]


def available_moltypes():
    """returns Table listing the available moltypes"""
    from cogent3.util.table import Table

    rows = []
    for n, m in _moltypes.items():
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


# constant instances of the core molecular types
ASCII = MolType(
    # A default type for text read from a file etc. when we don't
    # want to prematurely assume DNA or Protein. We therefore need to include
    # characters that could be present in any of those files.
    monomers="".join(ascii_letters),
    name="text",
    make_seq=new_sequence.Sequence,
    missing=IUPAC_missing,
)

DNA = MolType(
    monomers="".join(IUPAC_DNA_chars),
    ambiguities=IUPAC_DNA_ambiguities,
    name="dna",
    complements=IUPAC_DNA_ambiguities_complements,
    colors=NT_COLORS,
    make_seq=new_sequence.DnaSequence,
    pairing_rules=DNA_STANDARD_PAIRS,
    mw_calculator=DnaMW,
)

RNA = MolType(
    monomers="".join(IUPAC_RNA_chars),
    ambiguities=IUPAC_RNA_ambiguities,
    name="rna",
    complements=IUPAC_RNA_ambiguities_complements,
    colors=NT_COLORS,
    make_seq=new_sequence.RnaSequence,
    pairing_rules=RNA_STANDARD_PAIRS,
    mw_calculator=RnaMW,
)
#
PROTEIN = MolType(
    monomers="".join(IUPAC_PROTEIN_chars),
    ambiguities=IUPAC_PROTEIN_ambiguities,
    name="protein",
    colors=AA_COLORS,
    make_seq=new_sequence.ProteinSequence,
    mw_calculator=ProteinMW,
)

PROTEIN_WITH_STOP = MolType(
    monomers="".join(PROTEIN_WITH_STOP_chars),
    ambiguities=PROTEIN_WITH_STOP_ambiguities,
    name="protein_with_stop",
    colors=AA_COLORS,
    make_seq=new_sequence.ProteinWithStopSequence,
    mw_calculator=ProteinMW,
)
BYTES = MolType(
    # A default type for arbitrary chars read from a file etc. when we don't
    # want to prematurely assume _anything_ about the data.
    monomers=bytes(bytearray(range(2**8))),
    name="bytes",
    gap=None,
    missing=None,
    make_seq=new_sequence.ByteSequence,
)

# the None value catches cases where a moltype has no label attribute
_STYLE_DEFAULTS = {
    getattr(mt, "label", ""): defaultdict(
        _DefaultValue(f"ambig_{getattr(mt, 'label', '')}")
    )
    for mt in (ASCII, BYTES, DNA, RNA, PROTEIN, PROTEIN_WITH_STOP, None)
}

# build this at end of file
_moltypes = _make_moltype_dict()
