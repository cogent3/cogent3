import dataclasses
import functools
import itertools
import typing

from collections import defaultdict
from string import ascii_letters

import numpy

from cogent3.core import new_alphabet, new_sequence


OptStr = typing.Optional[str]
OptCallable = typing.Optional[typing.Callable]
SeqStrType = typing.Union[list[str], tuple[str, ...]]
StrORBytes = typing.Union[str, bytes]
StrORBytesORArray = typing.Union[str, bytes, numpy.ndarray]
SeqStrBytesType = typing.Union[list[StrORBytes], tuple[StrORBytes, ...]]
StrORArray = typing.Union[str, numpy.ndarray]

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
IUPAC_RNA_chars = ["U", "C", "A", "G"]
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

IUPAC_PROTEIN_ambiguities = {"B": ["N", "D"], "X": IUPAC_PROTEIN_chars, "Z": ["Q", "E"]}

PROTEIN_WITH_STOP_ambiguities = {
    "B": ["N", "D"],
    "X": PROTEIN_WITH_STOP_chars,
    "Z": ["Q", "E"],
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


@dataclasses.dataclass
class MolType:
    name: str
    monomers: dataclasses.InitVar[StrORBytes]
    make_seq: dataclasses.InitVar[typing.Type]
    gap: OptStr = IUPAC_gap
    missing: OptStr = IUPAC_missing
    complements: dataclasses.InitVar[typing.Optional[dict[str, str]]] = None
    ambiguities: typing.Optional[dict[str, tuple[str, ...]]] = None
    colors: dataclasses.InitVar[typing.Optional[dict[str, str]]] = None
    pairing_rules: typing.Optional[dict[str, dict[frozenset[str], bool]]] = None

    # private attributes to be delivered via properties
    _monomers: new_alphabet.CharAlphabet = dataclasses.field(init=False)
    _gapped: new_alphabet.CharAlphabet = dataclasses.field(init=False)
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

    def is_valid(self, seq: StrORArray) -> bool:
        """checks against most degenerate alphabet"""
        alpha = next(
            alpha
            for alpha in (self._degen_gapped, self._degen, self._gapped, self._monomers)
            if alpha
        )
        return alpha.is_valid(seq)

    def iter_alphabets(self):
        """yield the different defined alphabets"""
        alphas = (self._monomers, self._gapped, self._degen, self._degen_gapped)
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
        self, *, seq: str, name: OptStr = None, check_seq=True, **kwargs
    ):  # refactor: docstring
        if check_seq:
            assert self.is_valid(
                seq
            ), f"{seq[:4]!r} not valid for moltype {self.name!r}"
        return self._make_seq(moltype=self, seq=seq or "", name=name, **kwargs)

    @functools.singledispatchmethod
    def complement(self, seq: StrORBytesORArray) -> str:
        """converts a string or bytes into it's nucleic acid complement"""
        raise TypeError(f"{type(seq)} not supported")

    @complement.register
    def _(self, seq: str) -> str:
        return self.complement(seq.encode("utf8"))

    @complement.register
    def _(self, seq: bytes) -> str:
        return self._complement(seq).decode("utf8")

    @complement.register
    def _(self, seq: numpy.ndarray) -> str:
        return self.complement(self.degen_gapped_alphabet.array_to_bytes(seq))

    def rc(self, seq: str) -> str:
        """reverse reverse complement of a sequence"""
        return self.complement(seq)[::-1]

    @functools.singledispatchmethod
    def is_degenerate(self, seq: StrORBytesORArray) -> bool:
        """checks if a sequence contains degenerate characters"""
        raise TypeError(f"{type(seq)} not supported")

    @is_degenerate.register
    def _(self, seq: bytes) -> bool:  # refactor: docstring
        return self.is_degenerate(self.degen_gapped_alphabet.to_indices(seq))

    @is_degenerate.register
    def _(self, seq: str) -> bool:
        return self.is_degenerate(self.degen_gapped_alphabet.to_indices(seq))

    @is_degenerate.register
    def _(self, seq: numpy.ndarray) -> bool:
        # what index is the first degenerate character
        for index, val in enumerate(self.degen_gapped_alphabet):
            if val in self.ambiguities:
                break
        else:
            return False
        return (seq >= index).any()

    @functools.singledispatchmethod
    def is_gapped(self, seq) -> bool:
        """checks if a sequence contains gaps"""
        raise TypeError(f"{type(seq)} not supported")

    @is_gapped.register
    def _(self, seq: str) -> bool:
        return any(gap in seq for gap in self.gaps)

    @is_gapped.register
    def _(self, seq: bytes) -> bool:
        return self.is_gapped(seq.decode("utf8"))

    @is_gapped.register
    def _(self, seq: bytes) -> bool:
        return self.is_gapped(seq.decode("utf8"))

    @is_gapped.register
    def _(self, seq: numpy.ndarray) -> bool:
        gaps_indices = [self.degen_gapped_alphabet.index(gap) for gap in self.gaps]
        return numpy.isin(seq, gaps_indices).any()

    def get_degenerate_positions(
        self, seq: StrORBytesORArray, include_gap: bool = True
    ) -> list[int]:  # refactor: docstring
        seq = self.degen_gapped_alphabet.to_indices(seq)
        for index, val in enumerate(self.degen_gapped_alphabet):
            if include_gap and val in self.gap or val in self.ambiguities:
                break
        degens = seq >= index
        return numpy.where(degens)[0].tolist()

    @functools.singledispatchmethod
    def degap(self, seq) -> StrORBytesORArray:
        """removes all gap and missing characters from a sequence"""
        # refactor: design
        # cache translation callables (using alpabet.convert_alphabet)
        # previous implementation had support for tuples -- is this necessary?
        raise TypeError(f"{type(seq)} not supported")

    @degap.register
    def _(self, seq: bytes) -> bytes:
        return seq.translate(None, delete="".join(self.gaps).encode("utf8"))

    @degap.register
    def _(self, seq: str) -> str:
        degapped = self.degap(seq.encode("utf8"))
        return degapped.decode("utf8")

    @degap.register
    def _(self, seq: numpy.ndarray) -> numpy.ndarray:
        degapped = self.degen_gapped_alphabet.array_to_bytes(seq)
        return self.degen_gapped_alphabet.to_indices(self.degap(degapped))

    def is_ambiguity(self, query_motif: str) -> bool:
        """Return True if querymotif is an amibiguity character in alphabet.

        Parameters
        ----------
        query_motif
            the motif being queried.

        """
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

    def strand_symmetric_motifs(
        self, motif_length: int = 1
    ) -> set[tuple[str, str]]:  # refactor: docstring
        """returns ordered pairs of strand complementary motifs"""
        if len(self.alphabet) != 4:
            raise TypeError("moltype must be DNA or RNA")

        motif_set = self.alphabet.get_kmer_alphabet(k=motif_length)
        return {tuple(sorted([m, self.complement(m)])) for m in motif_set}

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
)

RNA = MolType(
    monomers="".join(IUPAC_RNA_chars),
    ambiguities=IUPAC_RNA_ambiguities,
    name="rna",
    complements=IUPAC_RNA_ambiguities_complements,
    colors=NT_COLORS,
    make_seq=new_sequence.RnaSequence,
    pairing_rules=RNA_STANDARD_PAIRS,
)
#
PROTEIN = MolType(
    monomers="".join(IUPAC_PROTEIN_chars),
    ambiguities=IUPAC_PROTEIN_ambiguities,
    name="protein",
    colors=AA_COLORS,
    make_seq=new_sequence.ProteinSequence,
)

PROTEIN_WITH_STOP = MolType(
    monomers="".join(PROTEIN_WITH_STOP_chars),
    ambiguities=PROTEIN_WITH_STOP_ambiguities,
    name="protein_with_stop",
    colors=AA_COLORS,
    make_seq=new_sequence.ProteinWithStopSequence,
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
