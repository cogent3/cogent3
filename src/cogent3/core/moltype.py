import dataclasses
import functools
import itertools
import json
import warnings
from collections import defaultdict
from collections.abc import Callable, Generator, Mapping
from string import ascii_letters
from typing import (
    TYPE_CHECKING,
    Any,
    Generic,
    Literal,
    TypeVar,
    cast,
    overload,
)

import numpy
import numpy.typing as npt

from cogent3.core import alphabet as c3_alphabet
from cogent3.core import sequence as c3_sequence
from cogent3.data.molecular_weight import DnaMW, ProteinMW, RnaMW, WeightCalculator
from cogent3.util import warning as c3warn
from cogent3.util.deserialise import register_deserialiser
from cogent3.util.misc import get_object_provenance

if TYPE_CHECKING:  # pragma: no cover
    from cogent3.core.table import Table

NumpyIntArrayType = npt.NDArray[numpy.integer]
StrOrBytes = TypeVar("StrOrBytes", str, bytes)

IUPAC_gap = "-"

IUPAC_missing = "?"

IUPAC_DNA_chars = "T", "C", "A", "G"
IUPAC_DNA_ambiguities: dict[str, frozenset[str]] = {
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
DNA_STANDARD_PAIRS: dict[frozenset[str], bool] = {
    frozenset(("A", "T")): True,
    frozenset(("C", "G")): True,
}

# note change in standard order from DNA
IUPAC_RNA_chars = "U", "C", "A", "G"
IUPAC_RNA_ambiguities: dict[str, frozenset[str]] = {
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
RNA_STANDARD_PAIRS: dict[frozenset[str], bool] = {
    frozenset(("A", "U")): True,  # True vs False for 'always' vs 'sometimes' pairing
    frozenset(("C", "G")): True,
    frozenset(("G", "U")): False,
}

# Watson-Crick RNA pairing only: GU pairs don't count as pairs
RNA_W_C_PAIRS: dict[frozenset[str], bool] = {
    frozenset(("A", "U")): True,
    frozenset(("C", "G")): True,
    frozenset(("U", "A")): True,
}

# RNA pairing with GU counted as standard pairs
RNA_G_U_PAIRS: dict[frozenset[str], bool] = {
    frozenset(("A", "U")): True,
    frozenset(("C", "G")): True,
    frozenset(("G", "C")): True,
    frozenset(("U", "A")): True,
    frozenset(("G", "U")): True,
    frozenset(("U", "G")): True,
}

# RNA pairing with GU, AA, GA, CA and UU mismatches allowed as weak pairs
RNA_EXTENDED_PAIRS: dict[frozenset[str], bool] = {
    frozenset({"A", "U"}): True,
    frozenset({"C", "G"}): True,
    frozenset({"G", "U"}): False,
    frozenset({"A"}): False,
    frozenset({"A", "G"}): False,
    frozenset({"A", "C"}): False,
    frozenset({"U"}): False,
}


def ambigs_and_missing(
    ambiguities: dict[str, frozenset[str]],
    gap_char: str | None,
) -> dict[str, frozenset[str]]:
    """returns new ambiguity dict with IUPAC missing as a key to most general code

    Notes
    -----
    Includes the gap character if gap is not None
    """
    val = max(ambiguities.values(), key=len)
    if gap_char:
        val = val | frozenset((gap_char,))
    return {**ambiguities, IUPAC_missing: val}


class MolTypeError(TypeError): ...


def make_pairs(
    *,
    pairs: dict[frozenset[str], bool] | None = None,
    monomers: tuple[str, ...] | None = None,
    gaps: str | frozenset[str] | None = None,
    degenerates: Mapping[str, set[str] | frozenset[str]] | None = None,
) -> dict[frozenset[str], bool]:
    """Makes a dict of symbol pairs (i,j) -> strictness.

    Expands pairs into all possible pairs using degen symbols.
    Strictness is True if i and j always pair, and False if they 'weakly' pair
    (e.g. GU pairs or if it is possible that they pair).

    If you want to make GU pairs count as 'always matching', pass in pairs
    that have (G,U) and (U, G) mapped to True rather than False.
    """
    pairs = pairs or {}
    monomers = monomers or ()
    gaps = gaps or ""
    degenerates = degenerates or {}

    result = pairs.copy()
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
    monomers: tuple[StrOrBytes, ...] | None = None,
    gaps: frozenset[StrOrBytes] | None = None,
    degenerates: dict[StrOrBytes, frozenset[StrOrBytes]] | None = None,
) -> dict[tuple[StrOrBytes, StrOrBytes], bool]:
    """Makes a dict of symbol pairs (i,j) -> strictness.

    Strictness is True if i and j always match and False if they sometimes
    match (e.g. A always matches A, but W sometimes matches R).
    """

    monomers = monomers or ()
    gaps = gaps or frozenset()
    degenerates = degenerates or {}
    result = {(i, i): True for i in monomers}
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
def _build_pairing_rules() -> dict[
    str, tuple[dict[frozenset[str], bool], dict[frozenset[str], bool]]
]:
    pairing_rules = {
        "Standard": RNA_STANDARD_PAIRS,
        "WC": RNA_W_C_PAIRS,
        "GU": RNA_G_U_PAIRS,
        "Extended": RNA_EXTENDED_PAIRS,
    }
    return {k: (v, make_pairs(pairs=v)) for k, v in pairing_rules.items()}


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

IUPAC_PROTEIN_ambiguities: dict[str, frozenset[str]] = {
    "B": frozenset(["N", "D"]),
    "X": frozenset(IUPAC_PROTEIN_chars),
    "Z": frozenset(["Q", "E"]),
}

PROTEIN_WITH_STOP_ambiguities: dict[str, frozenset[str]] = {
    "B": frozenset(["N", "D"]),
    "X": frozenset(PROTEIN_WITH_STOP_chars),
    "Z": frozenset(["Q", "E"]),
}


T = TypeVar("T")


class _DefaultValue(Generic[T]):
    def __init__(self, value: T) -> None:
        self.value = value

    def __call__(self) -> T:
        return self.value


# styling for moltype display


def _expand_colors(base: dict[str, str], colors: dict[str, str]) -> dict[str, str]:
    base = base.copy()
    base.update({ch: clr for chars, clr in colors.items() for ch in chars})
    return base


_gray = _DefaultValue("gray")
_base_colors: defaultdict[str, str] = defaultdict(_gray)

NT_COLORS = _expand_colors(
    _base_colors,
    {"A": "#FF0102", "C": "black", "G": "green", "T": "blue", "U": "blue"},
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


def _combined_chars(*parts: StrOrBytes) -> StrOrBytes:
    concat = parts[0]
    for part in parts[1:]:
        if not part or part in concat:
            continue
        concat = concat + part
    return concat


coerce_to_rna = c3_alphabet.convert_alphabet(
    b"tTucag-nrywskmbdhv?", b"UUUCAG-NRYWSKMBDHV?"
)
coerce_to_dna = c3_alphabet.convert_alphabet(
    b"uUtcag-nrywskmbdhv?", b"TTTCAG-NRYWSKMBDHV?"
)
coerce_to_protein = c3_alphabet.convert_alphabet(
    b"acdefghiklmnpqrstuvwy*-bxz?", b"ACDEFGHIKLMNPQRSTUVWY*-BXZ?"
)


@dataclasses.dataclass
class MolType(Generic[StrOrBytes]):
    """MolType handles operations that depend on the sequence type.

    Notes
    -----
    The only way to create sequences is via a MolType instance.
    The instance defines different alphabets that are used for data
    conversions.
    Create a moltype using the ``get_moltype()`` function.
    """

    name: str
    monomers: dataclasses.InitVar[StrOrBytes]
    make_seq: dataclasses.InitVar[type[c3_sequence.Sequence]]
    gap: str | None = IUPAC_gap
    missing: str | None = IUPAC_missing
    complements: dataclasses.InitVar[dict[str, str] | None] = None
    ambiguities: dict[str, frozenset[str]] | None = None
    colors: dataclasses.InitVar[dict[str, str] | None] = None
    pairing_rules: dict[frozenset[str], bool] | None = None
    mw_calculator: WeightCalculator | None = None
    coerce_to: Callable[[bytes], bytes] | None = None

    # private attributes to be delivered via properties
    _monomers: c3_alphabet.CharAlphabet[StrOrBytes] = dataclasses.field(init=False)
    _gapped: c3_alphabet.CharAlphabet[StrOrBytes] | None = dataclasses.field(init=False)
    _gapped_missing: c3_alphabet.CharAlphabet[StrOrBytes] | None = dataclasses.field(
        init=False
    )
    _degen: c3_alphabet.CharAlphabet[StrOrBytes] | None = dataclasses.field(init=False)
    _degen_gapped: c3_alphabet.CharAlphabet[StrOrBytes] | None = dataclasses.field(
        init=False
    )
    _colors: dict[str, str] = dataclasses.field(init=False)

    # how to connect this to the sequence constructor and avoid
    # circular imports
    _make_seq: Callable[..., c3_sequence.Sequence] = dataclasses.field(init=False)
    _complement: Callable[[bytes], bytes] | None = dataclasses.field(
        init=False, default=None
    )

    def __post_init__(
        self,
        monomers: StrOrBytes,
        make_seq: type,
        complements: dict[str, str] | None,
        colors: dict[str, str] | None,
    ) -> None:
        self._colors = colors or defaultdict(_DefaultValue("black"))
        self._make_seq = make_seq
        gap = c3_alphabet._coerce_to_type(monomers, self.gap or "")
        missing = c3_alphabet._coerce_to_type(monomers, self.missing or "")
        ambigs = c3_alphabet._coerce_to_type(monomers, "".join(self.ambiguities or ""))

        self._monomers = c3_alphabet.make_alphabet(
            chars=monomers,
            gap=None,
            missing=None,
            moltype=self,
        )
        self._degen = (
            c3_alphabet.make_alphabet(
                chars=_combined_chars(monomers, ambigs, missing),
                gap=None,
                missing=missing or None,
                moltype=self,
            )
            if ambigs
            else None
        )
        self._gapped = (
            c3_alphabet.make_alphabet(
                chars=_combined_chars(monomers, gap),
                gap=gap,
                missing=None,
                moltype=self,
            )
            if gap
            else None
        )
        self._gapped_missing = (
            c3_alphabet.make_alphabet(
                chars=_combined_chars(monomers, gap, missing),
                gap=gap,
                missing=missing,
                moltype=self,
            )
            if missing and gap
            else None
        )
        self._degen_gapped = (
            c3_alphabet.make_alphabet(
                chars=_combined_chars(monomers, gap, ambigs, missing),
                gap=gap,
                missing=missing,
                moltype=self,
            )
            if ambigs and gap
            else None
        )

        if self.missing and self.ambiguities:
            # this MUST occur after the creation of the different alphabets otherwise
            # missing data is added twice into an alphabet
            self.ambiguities = ambigs_and_missing(self.ambiguities, self.gap)

        if complements:
            # assume we have a nucleic acid moltype
            # any character not in defined complement set is defined as its
            # own complement
            if self.degen_gapped_alphabet is None:
                msg = (
                    "Degen gapped alphabet is None but needed for handling complements."
                )
                raise ValueError(msg)
            dest = "".join(
                complements.get(c, c)
                for c in cast(
                    "c3_alphabet.CharAlphabet[str]", self.degen_gapped_alphabet
                )
            )
            self._complement = c3_alphabet.convert_alphabet(
                self.degen_gapped_alphabet.as_bytes(),
                dest.encode("utf8"),
            )

    def __repr__(self) -> str:
        name = self.__class__.__name__
        return f"{name}({self.alphabet})"

    def __hash__(self) -> int:
        return id(self)

    def __eq__(self, other: object) -> bool:
        return id(self) == id(other)

    def __len__(self) -> int:
        return len(self._monomers)

    def __iter__(self) -> Generator[StrOrBytes]:
        yield from self._monomers

    @property
    def label(self) -> str:
        """synonym for name"""
        return self.name

    @property
    def alphabet(self) -> c3_alphabet.CharAlphabet[StrOrBytes]:
        """monomers"""
        return self._monomers

    @property
    def degen_alphabet(self) -> c3_alphabet.CharAlphabet[StrOrBytes] | None:
        """monomers + ambiguous characters"""
        return self._degen

    @property
    def gapped_alphabet(self) -> c3_alphabet.CharAlphabet[StrOrBytes] | None:
        """monomers + gap"""
        return self._gapped

    @property
    def degen_gapped_alphabet(self) -> c3_alphabet.CharAlphabet[StrOrBytes] | None:
        """monomers + gap + ambiguous characters"""
        return self._degen_gapped

    @property
    def gaps(self) -> frozenset[str]:  # refactor: docstring
        gaps = [char for char in (self.gap, self.missing) if char is not None]
        return frozenset(gaps)

    @property
    def matching_rules(self) -> dict[tuple[StrOrBytes, StrOrBytes], bool]:
        # Assumes no gaps or degnerates when monomers are bytes
        return make_matches(
            monomers=self._monomers,
            gaps=cast("frozenset[StrOrBytes]", self.gaps),
            degenerates=cast(
                "dict[StrOrBytes, frozenset[StrOrBytes]]", self.ambiguities
            ),
        )

    def is_valid(self, seq: str | bytes | NumpyIntArrayType) -> bool:
        """checks against most degenerate alphabet"""
        alpha = self.most_degen_alphabet()
        return alpha.is_valid(seq)

    def iter_alphabets(self) -> Generator[c3_alphabet.CharAlphabet[StrOrBytes]]:
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
        self,
        alphabet: c3_alphabet.CharAlphabet[Any],
        strict: bool = True,
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

    # This overshadows the make_seq attribute
    def make_seq(  # type: ignore[no-redef]
        self,
        *,
        seq: "str | bytes | c3_sequence.Sequence",
        name: str | None = None,
        check_seq: bool = True,
        **kwargs: Any,
    ) -> c3_sequence.Sequence:
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

        Notes
        -----
        If seq is a string, and the moltype has a coerce_to attribute, the
        string will be converted via that callable into a character set
        compatible with the moltype. Only applies to nucleic acid moltypes.
        """
        if hasattr(seq, "moltype"):
            seq = cast("c3_sequence.Sequence", seq)
            return seq if seq.moltype is self else seq.to_moltype(self)

        seq = cast("str | bytes", seq)
        if isinstance(seq, str):
            seq = self.coerce_to(seq.encode("utf8")) if self.coerce_to else seq

        if check_seq and not self.is_valid(seq):
            alpha = self.most_degen_alphabet()
            s = alpha.from_indices(seq)
            values = tuple(set(s) - set(map(str, alpha)))
            msg = f"{values} not valid for moltype {self.name!r} alphabet {alpha}"
            raise c3_alphabet.AlphabetError(msg)

        return self._make_seq(moltype=self, seq=seq, name=name, **kwargs)

    @c3warn.deprecated_callable(
        "2025.9", "no longer has an effect", is_discontinued=True
    )
    def make_array_seq(
        self,
        *args: Any,
        **kwargs: Any,
    ) -> c3_sequence.Sequence:
        """discontinued, use make_seq instead"""
        if args:
            kwargs["seq"] = args[0]

        return self.make_seq(**kwargs)  # type: ignore[attr-defined]

    @overload
    def complement(self, seq: str, validate: bool = True) -> str: ...

    @overload
    def complement(self, seq: bytes, validate: bool = True) -> bytes: ...

    @overload
    def complement(
        self, seq: NumpyIntArrayType, validate: bool = True
    ) -> NumpyIntArrayType: ...

    def complement(
        self, seq: str | bytes | NumpyIntArrayType, validate: bool = True
    ) -> str | bytes | NumpyIntArrayType:
        """converts a string or bytes into it's nucleic acid complement

        Parameters
        ----------
        seq
            sequence to be complemented
        validate
            if True, checks the sequence is validated against the most
            degenerate alphabet
        """
        if not isinstance(seq, (str, bytes, numpy.ndarray)):
            msg = f"{type(seq)} not supported"
            raise TypeError(msg)

        if validate and not self.is_valid(seq):
            msg = f"{seq[:4]!r} not valid for moltype {self.name!r}"
            raise c3_alphabet.AlphabetError(
                msg,
            )

        if isinstance(seq, str):
            return self.complement(seq.encode("utf8"), validate=False).decode("utf8")

        if isinstance(seq, bytes):
            if not self.is_nucleic:
                msg = f"{self.name!r} cannot complement"
                raise MolTypeError(msg)
            assert self._complement is not None
            return self._complement(seq)

        # Is a numpy array
        assert self.degen_gapped_alphabet is not None
        return self.degen_gapped_alphabet.to_indices(
            self.complement(
                self.degen_gapped_alphabet.array_to_bytes(seq),
                validate=False,
            ),
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

    def has_ambiguity(self, seq: str | bytes | NumpyIntArrayType) -> bool:
        """whether sequence has an ambiguity character"""
        if isinstance(seq, str):
            return any(self.is_ambiguity(c) for c in seq) if self.ambiguities else False

        if isinstance(seq, bytes):
            return self.has_ambiguity(seq.decode("utf8"))

        if isinstance(seq, numpy.ndarray):
            if not self.ambiguities:
                return False

            min_canonical = len(self.alphabet) - 1
            assert self.degen_alphabet is not None
            max_degen = len(self.degen_alphabet)
            result = ((seq > min_canonical) & (seq < max_degen)).any()
            return bool(result)

        msg = f"{type(seq)} not supported"
        raise TypeError(msg)

    def rc(self, seq: str, validate: bool = True) -> str:
        """reverse reverse complement of a sequence"""
        if validate and not self.is_valid(seq):
            msg = f"{seq[:4]!r} not valid for moltype {self.name!r}"
            raise c3_alphabet.AlphabetError(
                msg,
            )

        return self.complement(seq)[::-1]

    def is_degenerate(
        self, seq: str | bytes | NumpyIntArrayType, validate: bool = True
    ) -> bool:
        """checks if a sequence contains degenerate characters"""
        if not isinstance(seq, (str, bytes, numpy.ndarray)):
            msg = f"{type(seq)} not supported"
            raise TypeError(msg)

        if validate and not self.is_valid(seq):
            msg = f"{seq[:4]!r} not valid for moltype {self.name!r}"
            raise c3_alphabet.AlphabetError(
                msg,
            )

        if isinstance(seq, (str, bytes)):
            seq = self.most_degen_alphabet().to_indices(seq)

        if self.degen_alphabet is None:
            return False

        first_degen = cast("int", self.most_degen_alphabet().gap_index) + 1

        return bool((seq >= first_degen).any())

    def is_gapped(
        self, seq: str | bytes | NumpyIntArrayType, validate: bool = True
    ) -> bool:
        """checks if a sequence contains gaps"""
        if not isinstance(seq, (str, bytes, numpy.ndarray)):
            msg = f"{type(seq)} not supported"
            raise TypeError(msg)

        if validate and not self.is_valid(seq):
            msg = f"{seq[:4]!r} not valid for moltype {self.name!r}"
            raise c3_alphabet.AlphabetError(
                msg,
            )

        if isinstance(seq, bytes):
            seq = seq.decode("utf8")

        if isinstance(seq, str):
            return any(gap in seq for gap in self.gaps)

        # Is a numpy array
        alpha = self.most_degen_alphabet()
        gaps_indices = [alpha.index(gap) for gap in self.gaps]
        return bool(numpy.isin(seq, gaps_indices).any())

    def get_degenerate_positions(
        self,
        seq: str | bytes | NumpyIntArrayType,
        include_gap: bool = True,
        validate: bool = True,
    ) -> list[int]:
        """Return List of position indexs of degenerate characters in the sequence.

        Parameters
        ----------
        seq : str | bytes | NumpyIntArrayType
            the sequence to be used for getting degenerate positions
        include_gap : bool, optional
            if True, then the gap state together with ’canonical’ sates (A,C,G,T for DNA) will be considered non-ambiguous.
        validate : bool, optional
            if True, checks the sequence is validated for the alphabet
        """
        if validate and not self.is_valid(seq):
            msg = f"{seq[:4]!r} not valid for moltype {self.name!r}"
            raise c3_alphabet.AlphabetError(
                msg,
            )

        if self.degen_alphabet is None:
            return []

        alpha = self.most_degen_alphabet()
        seq_array = alpha.to_indices(seq)
        index = len(self) if include_gap else len(self) + 1
        degens = seq_array >= index
        return numpy.where(degens)[0].tolist()

    @overload
    def degap(self, seq: str, validate: bool = True) -> str: ...

    @overload
    def degap(self, seq: bytes, validate: bool = True) -> bytes: ...

    @overload
    def degap(
        self, seq: NumpyIntArrayType, validate: bool = True
    ) -> NumpyIntArrayType: ...

    def degap(
        self, seq: str | bytes | NumpyIntArrayType, validate: bool = True
    ) -> str | bytes | NumpyIntArrayType:
        """removes all gap and missing characters from a sequence"""
        # refactor: design
        # cache translation callables (using alpabet.convert_alphabet)
        # previous implementation had support for tuples -- is this necessary?
        if not isinstance(seq, (str, bytes, numpy.ndarray)):
            msg = f"{type(seq)} not supported"
            raise TypeError(msg)

        if validate and not self.is_valid(seq):
            msg = f"{seq[:4]!r} not valid for moltype {self.name!r}"
            raise c3_alphabet.AlphabetError(
                msg,
            )

        if isinstance(seq, bytes):
            if self.gap is None:
                msg = f"no gap character defined for {self.name!r}"
                raise TypeError(msg)

            return seq.translate(None, delete="".join(self.gaps).encode("utf8"))

        if isinstance(seq, str):
            degapped = self.degap(seq.encode("utf8"), validate=False)
            return degapped.decode("utf8")

        # Is a numpy array
        degen = self.most_degen_alphabet()
        if not degen.gap_char:
            msg = f"no gap character defined for {self.name!r}"
            raise TypeError(msg)

        degapped = degen.array_to_bytes(seq)
        return degen.to_indices(self.degap(degapped, validate=False))

    def is_ambiguity(self, query_motif: str, validate: bool = True) -> bool:
        """Return True if querymotif is an amibiguity character in alphabet.

        Parameters
        ----------
        query_motif
            the motif being queried.
        validate
            if True, checks the sequence is validated against the most
            degenerate alphabet
        """
        if not self.ambiguities:
            return False

        if validate and not self.is_valid(query_motif):
            msg = f"{query_motif[:4]!r} not valid for moltype {self.name!r}"
            raise c3_alphabet.AlphabetError(
                msg,
            )

        return query_motif in self.ambiguities or query_motif == self.missing

    def resolve_ambiguity(
        self,
        ambig_motif: str,
        alphabet: c3_alphabet.CharAlphabet[str]
        | c3_alphabet.KmerAlphabet[str]
        | None = None,
        allow_gap: bool = False,
        validate: bool = True,
    ) -> tuple[str, ...]:
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
            msg = f"{ambig_motif[:4]!r} not valid for moltype {self.name!r}"
            raise c3_alphabet.AlphabetError(
                msg,
            )

        # refactor: simplify
        ambigs = cast(
            "dict[str | None, frozenset[str] | set[str] | tuple[str, ...]]",
            self.ambiguities or {},
        )

        self_alphabet = cast("c3_alphabet.CharAlphabet[str]", self.alphabet)
        assert self.gap is not None
        ambiguities: dict[str | None, frozenset[str] | set[str] | tuple[str, ...]] = {
            **ambigs,
            self.gap: frozenset(self.gap),
            self.missing: {
                *cast("c3_alphabet.CharAlphabet[str]", self.alphabet),
                self.gap,
            },
        }
        for m in self_alphabet:
            ambiguities[m] = frozenset(m)
        if alphabet is None:
            word_alpha = self_alphabet.get_kmer_alphabet(
                k=len(ambig_motif),
                include_gap=allow_gap,
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
            raise c3_alphabet.AlphabetError(ambig_motif) from e

        result = tuple("".join(e) for e in itertools.product(*resolved))
        if alphabet:
            result = tuple(e for e in result if e in alphabet)

        if not result:
            raise c3_alphabet.AlphabetError(ambig_motif)

        return result

    @overload
    def strip_degenerate(self, seq: str) -> str: ...

    @overload
    def strip_degenerate(self, seq: bytes) -> bytes: ...

    @overload
    def strip_degenerate(self, seq: NumpyIntArrayType) -> NumpyIntArrayType: ...

    def strip_degenerate(
        self, seq: str | bytes | NumpyIntArrayType
    ) -> str | bytes | NumpyIntArrayType:
        """removes degenerate characters"""
        if isinstance(seq, bytes):
            monomers = self.alphabet

            assert self.degen_alphabet is not None
            ambigs = cast("set[str]", set(self.degen_alphabet) - set(monomers))
            func = c3_alphabet.convert_alphabet(
                src=monomers.as_bytes(),
                dest=monomers.as_bytes(),
                delete="".join(ambigs).encode("utf-8"),
            )
            return func(seq)

        if isinstance(seq, str):
            return self.strip_degenerate(seq.encode("utf-8")).decode("utf-8")

        if isinstance(seq, numpy.ndarray):
            assert self.degen_gapped_alphabet is not None
            return self.degen_gapped_alphabet.to_indices(
                self.strip_degenerate(self.degen_gapped_alphabet.array_to_bytes(seq)),
            )
        msg = f"{type(seq)} not supported"
        raise TypeError(msg)

    @overload
    def strip_bad(self, seq: str) -> str: ...
    @overload
    def strip_bad(self, seq: bytes) -> bytes: ...

    def strip_bad(self, seq: str | bytes) -> str | bytes:
        """Removes any symbols not in the alphabet."""
        if isinstance(seq, bytes):
            assert self.degen_gapped_alphabet is not None
            return self.degen_gapped_alphabet.array_to_bytes(
                self._strip_bad(self.degen_gapped_alphabet.to_indices(seq)),
            )

        if isinstance(seq, str):
            return self.strip_bad(seq.encode("utf8")).decode("utf8")

        msg = f"{type(seq)} not supported"
        raise TypeError(msg)

    def _strip_bad(self, seq: NumpyIntArrayType) -> NumpyIntArrayType:
        return seq[
            seq < len(cast("c3_alphabet.CharAlphabet[Any]", self.degen_gapped_alphabet))
        ]

    @overload
    def strip_bad_and_gaps(self, seq: str) -> str: ...
    @overload
    def strip_bad_and_gaps(self, seq: bytes) -> bytes: ...

    def strip_bad_and_gaps(self, seq: str | bytes) -> str | bytes:
        """Removes any symbols not in the alphabet, and any gaps.
        Since missing could be a gap, it is also removed."""

        if isinstance(seq, bytes):
            assert self.degen_gapped_alphabet is not None
            return self.degen_gapped_alphabet.array_to_bytes(
                self._strip_bad_and_gaps(self.degen_gapped_alphabet.to_indices(seq)),
            )

        if isinstance(seq, str):
            return self.strip_bad_and_gaps(seq.encode("utf8")).decode("utf8")

        msg = f"{type(seq)} not supported"
        raise TypeError(msg)

    def _strip_bad_and_gaps(self, seq: NumpyIntArrayType) -> NumpyIntArrayType:
        assert self.degen_gapped_alphabet is not None
        assert self.degen_gapped_alphabet.gap_index is not None
        assert self.degen_gapped_alphabet.missing_index is not None

        return seq[
            numpy.logical_and(
                numpy.logical_and(
                    seq < len(self.degen_gapped_alphabet),
                    seq != self.degen_gapped_alphabet.gap_index,
                ),
                seq != self.degen_gapped_alphabet.missing_index,
            )
        ]

    @overload
    def disambiguate(self, seq: str, method: str) -> str: ...

    @overload
    def disambiguate(self, seq: bytes, method: str) -> bytes: ...

    @overload
    def disambiguate(
        self, seq: NumpyIntArrayType, method: str
    ) -> NumpyIntArrayType: ...

    def disambiguate(
        self,
        seq: str | bytes | NumpyIntArrayType,
        method: str = "strip",
    ) -> str | bytes | NumpyIntArrayType:
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

        if method == "random":
            return self.random_disambiguate(seq)

        msg = f"method={method} not implemented"
        raise NotImplementedError(msg)

    @overload
    def random_disambiguate(self, seq: str) -> str: ...
    @overload
    def random_disambiguate(self, seq: bytes) -> bytes: ...
    @overload
    def random_disambiguate(self, seq: NumpyIntArrayType) -> NumpyIntArrayType: ...

    def random_disambiguate(
        self, seq: str | bytes | NumpyIntArrayType
    ) -> str | bytes | NumpyIntArrayType:
        """disambiguates a sequence by randomly selecting a non-degenerate character"""

        if isinstance(seq, str):
            assert self.ambiguities is not None
            rng = numpy.random.default_rng()
            return "".join(rng.choice(list(self.ambiguities.get(n, n))) for n in seq)

        if isinstance(seq, bytes):
            return self.random_disambiguate(seq.decode("utf8")).encode("utf8")

        if isinstance(seq, numpy.ndarray):
            assert self.degen_gapped_alphabet is not None
            return self.degen_gapped_alphabet.to_indices(
                self.random_disambiguate(
                    self.degen_gapped_alphabet.array_to_bytes(seq)
                ),
            )

        msg = f"{type(seq)} not supported"
        raise TypeError(msg)

    def count_gaps(self, seq: str | bytes) -> int:
        """returns the number of gap characters in a sequence"""
        if isinstance(seq, str):
            seq = seq.encode("utf8")

        if isinstance(seq, bytes):
            assert self.degen_gapped_alphabet is not None
            return self._count_gaps(self.degen_gapped_alphabet.to_indices(seq))

        msg = f"{type(seq)} not supported"
        raise TypeError(msg)

    def _count_gaps(self, seq: NumpyIntArrayType) -> int:
        assert self.degen_gapped_alphabet is not None
        return numpy.sum(seq == self.degen_gapped_alphabet.gap_index)

    def count_degenerate(self, seq: str | bytes) -> int:
        """returns the number of degenerate characters in a sequence"""
        if isinstance(seq, str):
            seq = seq.encode("utf8")

        if isinstance(seq, bytes):
            assert self.degen_gapped_alphabet is not None
            return self._count_degenerate(self.degen_gapped_alphabet.to_indices(seq))

        msg = f"{type(seq)} not supported"
        raise TypeError(msg)

    def _count_degenerate(self, seq: NumpyIntArrayType) -> int:
        assert self.degen_gapped_alphabet is not None
        assert self.degen_gapped_alphabet.gap_index is not None
        return int(numpy.sum(seq >= self.degen_gapped_alphabet.gap_index))

    def count_variants(self, seq: str | bytes) -> int:
        """Counts number of possible sequences matching the sequence, given
        any ambiguous characters in the sequence.

        Notes
        -----
        Uses self.ambiguitues to decide how many possibilities there are at
        each position in the sequence and calculates the permutations.
        """
        if isinstance(seq, bytes):
            seq = seq.decode("utf8")

        if isinstance(seq, str):
            assert self.ambiguities is not None
            return int(numpy.prod([len(self.ambiguities.get(c, c)) for c in seq]))

        msg = f"{type(seq)} not supported"
        raise TypeError(msg)

    def degenerate_from_seq(self, seq: str) -> str:
        """Returns least degenerate symbol that encompasses a set of characters"""
        if not self.ambiguities:
            assert (
                self.missing is not None
            )  # Is this always valid? or should None be added to the return type annotation?
            return self.missing

        symbols = frozenset(seq)

        # if the length of the set is 1, return the single character
        if len(symbols) == 1:
            return next(iter(symbols))

        # all degenerate symbols should be added to the sets that they encompass
        ambigs = self.ambiguities or {}
        degens = ambigs.copy()
        for degen1 in degens:
            for degen2, char_set in self.ambiguities.items():
                if self.ambiguities[degen1] <= char_set:
                    degens[degen2] = degens[degen2].union(frozenset(degen1))

        # reverse the mapping between degenerate symbols and their encompassing sets
        inv_degens = {frozenset(val): key for key, val in list(degens.items())}

        # add gap and missing characters to the mapping
        assert self.gap is not None
        assert self.degen_gapped_alphabet is not None
        assert self.missing is not None

        inv_degens[frozenset(self.gap)] = self.gap
        inv_degens[
            frozenset(cast("c3_alphabet.CharAlphabet[str]", self.degen_gapped_alphabet))
        ] = self.missing

        # if we exactly match a set of symbols, return the corresponding degenerate
        if result := inv_degens.get(symbols):
            return result

        # if we don't exactly match, sort all encompassing sets and return the
        # least degenerate one
        encompassing = [chars for chars in inv_degens if chars.issuperset(symbols)]
        encompassing = sorted(encompassing, key=len)
        return inv_degens[encompassing[0]]

    def strand_symmetric_motifs(
        self,
        motif_length: int = 1,
    ) -> set[tuple[str, str]]:  # refactor: docstring
        """returns ordered pairs of strand complementary motifs"""
        if len(self.alphabet) != 4:
            msg = "moltype must be DNA or RNA"
            raise TypeError(msg)

        motif_set = self.alphabet.get_kmer_alphabet(k=motif_length)
        return {
            cast("tuple[str, str]", tuple(sorted([m, self.complement(m)])))
            for m in motif_set
        }

    def mw(self, seq: str, method: str = "random", delta: float | None = None) -> float:
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

        if self.mw_calculator is None:
            msg = "moltype has no molecurlar weight calculator specified"
            raise ValueError(msg)
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
        return all(pair in m for pair in zip(first, second, strict=False))

    def can_pair(self, first: str, second: str) -> bool:
        pairs = make_pairs(
            pairs=self.pairing_rules,
            monomers=cast("c3_alphabet.CharAlphabet[str]", self._monomers),
            gaps=self.gaps,
            degenerates=self.ambiguities,
        )
        sec = list(second)
        sec.reverse()
        return all(frozenset(pair) in pairs for pair in zip(first, sec, strict=False))

    def can_mispair(self, first: str, second: str) -> bool:
        """Returns True if any position in self could mispair with other.

        Notes
        -----
        Pairing occurs in reverse order, i.e. last position of other with
        first position of self, etc.

        Truncates at length of shorter sequence.

        Gaps are always counted as possible mispairs, as are weak pairs like GU.
        """
        p = cast("dict[frozenset[str], bool]", self.pairing_rules)
        if not first or not second:
            return False

        for pair in zip(first, second[::-1], strict=False):
            if not p.get(frozenset(pair), None):
                return True
        return False

    def get_css_style(
        self,
        colors: dict[str, str] | None = None,
        font_size: int = 12,
        font_family: str = "Lucida Console",
    ) -> tuple[list[str], defaultdict[str | bytes, str]]:
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
        css_colors = cast("dict[str | bytes, str]", colors or self._colors)
        # !important required to stop some browsers over-riding the style sheet ...!!
        template = (
            '.%s_%s{font-family: "%s",monospace !important; '
            "font-size: %dpt !important; color: %s; }"
        )
        label = self.name if self.name in _STYLE_DEFAULTS else ""
        styles = _STYLE_DEFAULTS[label].copy()
        styles.update(
            {
                c: f"{c}_{label}"
                for c in cast(
                    "list[str]", [*list(self.alphabet), "terminal_ambig"]
                )  # the type of the colors dict is a lie as it could be a bytes alphabet
            },
        )

        style_chars = list(styles)
        style_chars.append("ambig")

        css = [
            template % (char, label, font_family, font_size, css_colors[char])
            for char in style_chars
        ]

        return css, styles

    def most_degen_alphabet(self) -> c3_alphabet.CharAlphabet[StrOrBytes]:
        """returns the most degenerate alphabet for this instance"""
        return _most_degen_alphabet(self)

    def to_regex(self, seq: str) -> str:
        """returns a regex pattern with ambiguities expanded to a character set"""
        if not self.is_valid(seq):
            msg = f"'{seq}' is invalid for this moltype"
            raise ValueError(msg)

        degen_indices = self.get_degenerate_positions(seq=seq, include_gap=False)
        seq_list = list(seq)  # seq can now be modified

        assert len(degen_indices) == 0 or self.ambiguities is not None
        for index in degen_indices:
            expanded = cast("dict[str, frozenset[str]]", self.ambiguities)[
                seq_list[index]
            ]
            seq_list[index] = f"[{''.join(sorted(expanded))}]"
        return "".join(seq_list)

    def to_rich_dict(self, **kwargs: Any) -> dict[str, str]:
        """returns dict suitable for serialisation"""
        return {"type": get_object_provenance(self), "moltype": self.label}

    def to_json(self) -> str:
        """returns result of json formatted string"""
        data = self.to_rich_dict()
        return json.dumps(data)


@functools.cache
def _most_degen_alphabet(
    mt: MolType[StrOrBytes],
) -> c3_alphabet.CharAlphabet[StrOrBytes]:
    """returns the most degenerate alphabet for this instance"""
    return next(mt.iter_alphabets())


@register_deserialiser(
    get_object_provenance(MolType),
    "cogent3.core.moltype.MolType",
    "cogent3.core.c3_moltype.MolType",
)
def deserialise_c3_moltype(data: dict[str, str]) -> MolType[Any]:
    return get_moltype(data["moltype"])


def _make_moltype_dict() -> dict[str, MolType[Any]]:
    """make a dictionary of local name space molecular types"""
    env = globals()
    moltypes: dict[str, MolType[Any]] = {}
    for obj in env.values():
        if not isinstance(obj, MolType):
            continue
        moltypes[obj.name] = obj

    return moltypes


MolTypeLiteral = Literal["dna", "rna", "protein", "protein_with_stop", "text", "bytes"]


@c3warn.deprecated_args("2025.9", "no longer has an effect", discontinued="new_type")
def get_moltype(name: MolTypeLiteral | MolType[Any] | None) -> MolType[Any]:
    """returns the moltype with the matching name attribute"""
    if name is None:
        msg = "no moltype specified, defaulting to ASCII"
        warnings.warn(msg, UserWarning, stacklevel=3)
        return ASCII

    if isinstance(name, MolType):
        return name

    if name.lower() not in _moltypes:
        msg = f"unknown moltype {name!r}"
        raise ValueError(msg)

    return _moltypes[name.lower()]


def available_moltypes() -> "Table":
    """returns Table listing the available moltypes"""
    from cogent3.core.table import Table

    rows: list[list[str | int]] = []
    for n, m in _moltypes.items():
        v = str(m)
        num = len(list(m))
        if num > 10:
            v = f"{v[:39]}..."
        rows.append([n, num, v])

    header = ["Abbreviation", "Number of states", "Moltype"]
    title = "Specify a moltype by the Abbreviation (case insensitive)."

    result = Table(header=header, data=rows, title=title, index_name="Abbreviation")
    result.format_column("Abbreviation", lambda x: repr(str(x)))
    return result.sorted(columns=["Number of states", "Abbreviation"])


# constant instances of the core molecular types
ASCII = MolType(
    # A default type for text read from a file etc. when we don't
    # want to prematurely assume DNA or Protein. We therefore need to include
    # characters that could be present in any of those files.
    monomers="".join(ascii_letters),
    name="text",
    make_seq=c3_sequence.Sequence,
)

DNA = MolType(
    monomers="".join(IUPAC_DNA_chars),
    ambiguities=IUPAC_DNA_ambiguities,
    name="dna",
    complements=IUPAC_DNA_ambiguities_complements,
    colors=NT_COLORS,
    make_seq=c3_sequence.DnaSequence,
    pairing_rules=DNA_STANDARD_PAIRS,
    mw_calculator=DnaMW,
    coerce_to=coerce_to_dna,
)

RNA = MolType(
    monomers="".join(IUPAC_RNA_chars),
    ambiguities=IUPAC_RNA_ambiguities,
    name="rna",
    complements=IUPAC_RNA_ambiguities_complements,
    colors=NT_COLORS,
    make_seq=c3_sequence.RnaSequence,
    pairing_rules=RNA_STANDARD_PAIRS,
    mw_calculator=RnaMW,
    coerce_to=coerce_to_rna,
)
PROTEIN = MolType(
    monomers="".join(IUPAC_PROTEIN_chars),
    ambiguities=IUPAC_PROTEIN_ambiguities,
    name="protein",
    colors=AA_COLORS,
    make_seq=c3_sequence.ProteinSequence,
    mw_calculator=ProteinMW,
    coerce_to=coerce_to_protein,
)

PROTEIN_WITH_STOP = MolType(
    monomers="".join(PROTEIN_WITH_STOP_chars),
    ambiguities=PROTEIN_WITH_STOP_ambiguities,
    name="protein_with_stop",
    colors=AA_COLORS,
    make_seq=c3_sequence.ProteinWithStopSequence,
    mw_calculator=ProteinMW,
    coerce_to=coerce_to_protein,
)
BYTES = MolType(
    # A default type for arbitrary chars read from a file etc. when we don't
    # want to prematurely assume _anything_ about the data.
    monomers=bytes(bytearray(range(2**8))),
    name="bytes",
    make_seq=c3_sequence.ByteSequence,
    gap=None,
    missing=None,
)

# the None value catches cases where a moltype has no label attribute
_STYLE_DEFAULTS: dict[str, defaultdict[str | bytes, str]] = {
    getattr(mt, "name", ""): defaultdict(
        _DefaultValue(f"ambig_{getattr(mt, 'name', '')}"),
    )
    for mt in (ASCII, BYTES, DNA, RNA, PROTEIN, PROTEIN_WITH_STOP, None)
}

# build this at end of file
_moltypes = _make_moltype_dict()
