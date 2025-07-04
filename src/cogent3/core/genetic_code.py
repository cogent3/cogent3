"""Translates RNA or DNA string to amino acid sequence.

Notes
 ----

* is used to denote termination (as per NCBI standard).

Although the genetic code objects convert DNA to RNA and vice
versa, lists of codons that they produce will be provided in DNA format.
"""

import collections
import contextlib
import dataclasses
import functools
import itertools
import typing

import numpy

from cogent3.core import alphabet as c3_alphabet
from cogent3.core import moltype as c3_moltype
from cogent3.core.table import Table
from cogent3.util import warning as c3warn

if typing.TYPE_CHECKING:  # pragma: no cover
    from cogent3.core.sequence import Sequence

OptStr = typing.Optional[str]
SetStr = set[str]
ConverterType = typing.Callable[[bytes, bytes], bytes]
StrORInt = typing.Union[str, int]
StrORBytesORArray = typing.Union[str, numpy.ndarray]
GeneticCodeChoiceType = typing.Union[int, str, "GeneticCode"]


class GeneticCodeError(Exception):
    pass


class GeneticCodeInitError(ValueError, GeneticCodeError):
    pass


class InvalidCodonError(KeyError, GeneticCodeError):
    pass


def _make_mappings(
    codons: c3_alphabet.KmerAlphabet,
    code_sequence: str,
) -> tuple[dict[str, str], dict[str, SetStr], tuple[str, ...]]:
    """makes amino acid / codon mappings and stop codon group

    Parameters
    ----------
    codons
        alphabet of codons
    code_sequence
        64-character string containing NCBI representation of the genetic code.

    Returns
    -------
    codon to amino acid mapping, the reverse mapping, the set of stop codons
    """
    stops = []
    codon_to_aa = {}
    aa_to_codon = collections.defaultdict(set)
    for codon, aa in zip(codons, code_sequence, strict=False):
        if aa == "*":
            stops.append(codon)
        codon_to_aa[codon] = aa
        aa_to_codon[aa].add(codon)
    return codon_to_aa, aa_to_codon, tuple(stops)


def _get_start_codon_indices(start_codon_map: str) -> tuple[int, ...]:
    return tuple(i for i, start in enumerate(start_codon_map) if start == "M")


def _make_converter(
    kmer_alpha: c3_alphabet.KmerAlphabet,
    codons: tuple[str, ...],
    code_sequence: str,
) -> typing.Callable[[bytes, bytes], bytes]:
    """returns a converter of codon indices into amino acid indices

    Parameters
    ----------
    kmer_alpha
        a complete trinucleotide alphabet
    codons
        the order of codons corresponding to code_sequence
    code_sequence
        NCBI genetic code sequence plus gap and missing states

    Returns
    -------
    callable for converting codon sequence as bytes to amino acid sequence
    as bytes
    """
    # we get the index of the codon in the kmer alphabet
    kmers = numpy.array(
        [kmer_alpha.to_index(codon) for codon in codons],
        dtype=numpy.uint8,
    )
    return c3_alphabet.convert_alphabet(kmers.tobytes(), code_sequence.encode("utf8"))


@dataclasses.dataclass
class GeneticCode:
    """Holds codon to amino acid mapping, and vice versa.

    Notes
    -----
    We add additional states to the genetic code to represent gapped codons
    and missing data.
    """

    ID: int
    name: str
    ncbi_code_sequence: dataclasses.InitVar[str]
    ncbi_start_codon_map: dataclasses.InitVar[str]
    moltype: c3_moltype.MolType = c3_moltype.DNA
    _codon_to_aa: dict[str, str] = dataclasses.field(init=False, default_factory=dict)
    _aa_to_codon: dict[str, SetStr] = dataclasses.field(
        init=False,
        default_factory=dict,
    )
    _sense_codons: tuple[str, ...] = dataclasses.field(init=False)
    _stop_codons: tuple[str, ...] = dataclasses.field(init=False)
    _start_codons: tuple[str, ...] = dataclasses.field(init=False)
    codons: c3_alphabet.KmerAlphabet = dataclasses.field(init=False)
    anticodons: tuple[str, ...] = dataclasses.field(init=False)
    # callables for translating on the plus strand, or the minus strand
    _translate_plus: ConverterType = dataclasses.field(init=False)
    _translate_minus: ConverterType = dataclasses.field(init=False)

    def __post_init__(self, ncbi_code_sequence: str, ncbi_start_codon_map: str) -> None:
        alpha = self.moltype.alphabet.with_gap_motif(include_missing=True)
        trinuc_alpha = alpha.get_kmer_alphabet(
            k=3,
            include_gap=True,
        )
        code_seq = f"{ncbi_code_sequence}-X"
        start_map = f"{ncbi_start_codon_map}--"
        self._codon_to_aa, self._aa_to_codon, self._stop_codons = _make_mappings(
            trinuc_alpha,
            code_seq,
        )
        self.codons = trinuc_alpha
        self._start_codons = tuple(
            self.codons[i] for i in _get_start_codon_indices(start_map)
        )
        self._sense_codons = tuple(
            c for c in self.codons if self._codon_to_aa[c] not in "*-X"
        )
        self.anticodons = tuple(self.moltype.rc(codon) for codon in self.codons)
        self._translate_plus = _make_converter(trinuc_alpha, self.codons, code_seq)
        self._translate_minus = _make_converter(trinuc_alpha, self.anticodons, code_seq)

    def __hash__(self) -> int:
        return hash(self.name)

    def __eq__(self, other: object) -> bool:
        """Allows two GeneticCode objects to be compared to each other."""
        return self.name == getattr(other, "name", None)

    @property
    def stop_codons(self) -> tuple[str, ...]:
        return self._stop_codons

    @property
    def start_codons(self) -> tuple[str, ...]:
        return self._start_codons

    @property
    def sense_codons(self) -> tuple[str, ...]:
        return self._sense_codons

    def to_table(self) -> Table:
        """returns aa to codon mapping as a cogent3 Table"""
        rows = []
        headers = ["aa", "IUPAC code", "codons"]
        for code, aa in c3_moltype.IUPAC_PROTEIN_code_aa.items():
            codons = ",".join(self[code])
            row = [aa, code, codons]
            rows.append(row)
        return Table(header=headers, data=rows, title=self.name)

    def __repr__(self) -> str:
        display = self.to_table()
        return str(display)

    def _repr_html_(self) -> str:
        """Returns the html representation of GeneticCode."""
        display = self.to_table()
        display.set_repr_policy(show_shape=False)
        return display._repr_html_()

    def __getitem__(self, item) -> str | SetStr:
        """Returns amino acid corresponding to codon, or codons for an aa.

        Returns [] for empty list of codons, 'X' for unknown amino acid.
        """
        item = str(item)
        if len(item) == 1:  # amino acid
            return self._aa_to_codon.get(item, set())

        if len(item) != 3:
            msg = f"Codon or aa {item!r} has wrong length"
            raise InvalidCodonError(msg)

        key = item.upper()
        key = key.replace("U", "T")
        return self._codon_to_aa.get(key, "X")

    def translate(
        self,
        dna: StrORBytesORArray,
        start: int = 0,
        rc: bool = False,
        incomplete_ok: bool = True,
    ) -> str:
        """Translates DNA to protein.

        Parameters
        ----------
        dna
            a string of nucleotides
        start
            position to begin translation (used to implement frames)
        rc
            if True, returns the translation of the reverse complement sequence
        incomplete_ok
            if True, translates codons that are a mix of gaps and bases
            as a gap. If False, raises an AlphabetError on those incomplete
            cases.

        Notes
        -----
        Sequences are truncated to be a multiple of 3.
        Codons containing ambiguous nucleotides are translated as 'X',
        codons containing a gap character are translated as '-' unless
        incomplete_ok is False. Codons with a mix of ambiguous nucleotides
        are translated as 'X'.

        Returns
        -------
        The amino acid sequence as a string.
        """
        # refactor: design
        # previous implementation converted mix of ambigs and nt to '?', potentially
        # we should change the behaviour to match this.
        if start:
            dna = dna[start:]

        if diff := len(dna) % 3:
            dna = dna[:-diff]

        # convert to indices and then bytes
        if incomplete_ok:
            seq = self.codons.to_indices(dna)
            missing = seq >= len(self.codons)
        else:
            codons = self.get_alphabet(
                include_stop=True,
                include_gap=True,
                include_missing=True,
            )
            seq = codons.to_indices(dna)
            missing = seq >= len(codons)
        # any codon indices include missing data, they will have an index ==
        # length of the alphabet employed, so we need to modify that index
        # to represent the missing state which is X
        if missing.any():
            seq[missing] = len(self._codon_to_aa) - 1

        if rc:
            return self._translate_minus(seq.tobytes()).decode("utf8")[::-1]

        return self._translate_plus(seq.tobytes()).decode("utf8")

    @functools.singledispatchmethod
    def sixframes(self, seq) -> typing.Iterable[tuple[str, int, str]]:
        """Returns the six reading frames of the genetic code.

        Returns
        -------
        A dictionary with keys (strand, start) where strand is "+"/"-"
        """
        return self.sixframes(str(seq))

    @sixframes.register
    def _(self, seq: str) -> typing.Iterable[tuple[str, int, str]]:
        """Returns the six reading frames of the genetic code.

        Returns
        -------
        A dictionary with keys (strand, start) where strand is "+"/"-"
        """

        for strand, start in itertools.product(("+", "-"), range(3)):
            yield strand, start, self.translate(seq, start, rc=strand == "-")

    def is_stop(self, codon: str) -> bool:
        """Returns True if codon is a stop codon, False otherwise."""
        return self[codon] == "*"

    def is_start(self, codon: str) -> bool:
        """Returns True if codon is a start codon, False otherwise."""
        fixed_codon = codon.upper().replace("U", "T")
        return fixed_codon in self.start_codons

    def get_alphabet(
        self,
        *,
        include_gap: bool = False,
        include_missing: bool = False,
        include_stop: bool = False,
    ) -> c3_alphabet.SenseCodonAlphabet:
        """returns a codon alphabet

        Parameters
        ----------
        include_gap
            alphabet includes the gap motif
        include_missing
            alphabet includes the missing state as 3 * IUPAC_missing
        include_stop
            if True, this is just a kmer alphabet

        Notes
        -----
        If include_stop, the returned alphabet includes all codons,
        and is thus a KmerAlphabet. Otherwise its a SenseCodonAlphabet.
        """
        if include_stop:
            words = tuple(self.moltype.alphabet.get_kmer_alphabet(k=3))
        else:
            words = tuple(self.sense_codons)
        return _get_code_alphabet(
            words,
            self.moltype.gapped_alphabet.gap_char,
            self.moltype.degen_gapped_alphabet.missing_char,
            include_gap,
            include_missing,
            self.moltype.degen_gapped_alphabet,
        )

    def to_regex(self, seq: typing.Union[str, "Sequence"]) -> str:
        """returns a regex pattern with an amino acid expanded to its codon set

        Parameters
        ----------
        seq
            a Sequence or string of amino acids
        """
        from .moltype import PROTEIN_WITH_STOP_ambiguities as ambigs

        seq = list(str(seq))
        mappings = []
        for aa in seq:
            aa = ambigs.get(aa, [aa])
            codons = []
            for a in aa:
                codons.extend(self[a])

            # we create a regex non-capturing group for each amino acid
            mappings.append(f"(?:{'|'.join(codons)})")

        return "".join(mappings)

    @property
    def synonyms(self):
        return self._aa_to_codon


@functools.cache
def _get_code_alphabet(
    words: tuple[str, ...],
    gap_char: str | None,
    missing_char: str | None,
    include_gap: bool,
    include_missing: bool,
    monomers: c3_alphabet.CharAlphabet,
) -> c3_alphabet.SenseCodonAlphabet:
    gap = gap_char * 3 if include_gap else None
    missing = missing_char * 3 if include_missing else None
    if include_gap:
        words += (gap,)

    if include_missing:
        words += (missing,)

    return c3_alphabet.SenseCodonAlphabet(
        words=words,
        monomers=monomers,
        gap=gap,
    )


_mapping_cols = "ncbi_code_sequence", "ID", "name", "ncbi_start_codon_map"
# code mappings are based on the product of bases in order TCAG
code_mapping = (
    (
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        1,
        "Standard",
        "---M---------------M---------------M----------------------------",
    ),
    (
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        2,
        "Vertebrate Mitochondrial",
        "--------------------------------MMMM---------------M------------",
    ),
    (
        "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        3,
        "Yeast Mitochondrial",
        "----------------------------------MM---------------M------------",
    ),
    (
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        4,
        "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate "
        "Mitochondrial; Mycoplasma; Spiroplasma",
        "--MM---------------M------------MMMM---------------M------------",
    ),
    (
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        5,
        "Invertebrate Mitochondrial",
        "---M----------------------------MMMM---------------M------------",
    ),
    (
        "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        6,
        "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
        "-----------------------------------M----------------------------",
    ),
    (
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        9,
        "Echinoderm Mitochondrial; Flatworm Mitochondrial",
        "-----------------------------------M---------------M------------",
    ),
    (
        "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        10,
        "Euplotid Nuclear",
        "-----------------------------------M----------------------------",
    ),
    (
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        11,
        "Bacterial, Archaeal and Plant Plastid",
        "---M---------------M------------MMMM---------------M------------",
    ),
    (
        "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        12,
        "Alternative Yeast Nuclear",
        "-------------------M---------------M----------------------------",
    ),
    (
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        13,
        "Ascidian Mitochondrial",
        "---M------------------------------MM---------------M------------",
    ),
    (
        "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        14,
        "Alternative Flatworm Mitochondrial",
        "-----------------------------------M----------------------------",
    ),
    (
        "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        15,
        "Blepharisma Macronuclear",
        "-----------------------------------M----------------------------",
    ),
    (
        "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        16,
        "Chlorophycean Mitochondrial",
        "-----------------------------------M----------------------------",
    ),
    (
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        21,
        "Trematode Mitochondrial",
        "-----------------------------------M---------------M------------",
    ),
    (
        "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        22,
        "Scenedesmus obliquus Mitochondrial",
        "-----------------------------------M----------------------------",
    ),
    (
        "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        23,
        "Thraustochytrium Mitochondrial",
        "--------------------------------M--M---------------M------------",
    ),
    (
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        24,
        "Rhabdopleuridae Mitochondrial",
        "---M---------------M---------------M---------------M------------",
    ),
    (
        "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        25,
        "Candidate Division SR1 and Gracilibacteria",
        "---M-------------------------------M---------------M------------",
    ),
    (
        "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        26,
        "Pachysolen tannophilus Nuclear",
        "-------------------M---------------M----------------------------",
    ),
    (
        "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        27,
        "Karyorelict Nuclear",
        "-----------------------------------M----------------------------",
    ),
    (
        "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        28,
        "Condylostoma Nuclear",
        "-----------------------------------M----------------------------",
    ),
    (
        "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        29,
        "Mesodinium Nuclear",
        "-----------------------------------M----------------------------",
    ),
    (
        "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        30,
        "Peritrich Nuclear",
        "-----------------------------------M----------------------------",
    ),
    (
        "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        31,
        "Blastocrithidia Nuclear",
        "-----------------------------------M----------------------------",
    ),
    (
        "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        32,
        "Balanophoraceae Plastid",
        "---M---------------M------------MMMM---------------M------------",
    ),
    (
        "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        33,
        "Cephalodiscidae Mitochondrial",
        "---M---------------M---------------M---------------M------------",
    ),
)
_CODES = {}
for mapping in code_mapping:
    code = GeneticCode(**dict(zip(_mapping_cols, mapping, strict=False)))
    _CODES[code.ID] = code
    _CODES[code.name] = code
    _CODES[code] = code


DEFAULT = _CODES[1]


@c3warn.deprecated_args("2025.9", "no longer has an effect", discontinued="new_type")
def get_code(code_id: StrORInt = 1) -> GeneticCode:
    """returns the genetic code

    Parameters
    ----------
    code_id
        genetic code identifier, name, number or string(number), defaults to
        standard genetic code
    """
    # refactor: simplify
    # added for compatibility with previous API, should be removed when
    # support for old style is dropped
    code_id = code_id or 1

    with contextlib.suppress((ValueError, TypeError)):
        code_id = int(code_id)

    if code_id not in _CODES:
        msg = f"Unknown genetic code {code_id!r}"
        raise GeneticCodeError(msg)
    return _CODES[code_id]


def available_codes():
    """returns Table listing the available genetic codes"""
    rows = [(k, code.name) for k, code in _CODES.items() if isinstance(k, int)]
    header = ["Code ID", "Name"]
    return Table(
        header=header,
        data=rows,
        index_name="Code ID",
        title="Specify a genetic code using either 'Name' or "
        "Code ID (as an integer or string)",
    )
