#!/usr/bin/env python
"""Translates RNA or DNA string to amino acid sequence.

NOTE: * is used to denote termination (as per NCBI standard).
NOTE: Although the genetic code objects convert DNA to RNA and vice
versa, lists of codons that they produce will be provided in DNA format.
"""
import re

from itertools import product

from cogent3.util.table import Table


__author__ = "Greg Caporaso and Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Greg Caporaso", "Rob Knight", "Peter Maxwell", "Thomas La"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Greg Caporaso"
__email__ = "caporaso@colorado.edu"
__status__ = "Production"

maketrans = str.maketrans


class GeneticCodeError(Exception):
    pass


class GeneticCodeInitError(ValueError, GeneticCodeError):
    pass


class InvalidCodonError(KeyError, GeneticCodeError):
    pass


_dna_trans = maketrans("TCAG", "AGTC")


def _simple_rc(seq):
    """simple reverse-complement: works only on unambiguous uppercase DNA"""
    return seq.translate(_dna_trans)[::-1]


_bases = "TCAG"


class GeneticCode:
    """Holds codon to amino acid mapping, and vice versa.

    Use the `get_code()` function to get one of the included code instances. These are created as follows.

    >>> code_sequence = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    >>> gc = GeneticCode(code_sequence)
    >>> sgc['UUU'] == 'F'
    >>> sgc['TTT'] == 'F'
    >>> sgc['F'] == ['TTT', 'TTC']          #in arbitrary order
    >>> sgc['*'] == ['TAA', 'TAG', 'TGA']   #in arbitrary order

    code_sequence : 64 character string containing NCBI genetic code translation

    GeneticCode is immutable once created.
    """

    # class data: need the bases, the list of codons in UUU -> GGG order, and
    # a mapping from positions in the list back to codons. These should be the
    # same for all GeneticCode instances, and are immutable (therefore
    # private).
    _nt = _bases
    _codons = tuple(map("".join, product(_bases, _bases, _bases)))

    def __init__(self, code_sequence, ID=None, name=None, start_codon_sequence=None):
        """Returns new GeneticCode object.

        code_sequence :

        Parameters
        ----------
        code_sequence : str
            64-character string containing NCBI representation of the genetic code.
        ID
            Identifier
        name
            name of the Genetic code
        start_codon_sequence
            64-character string where the '-' character indicates the corresponding
            position of code_sequence **is not** a start codon
        """
        if len(code_sequence) != 64:
            raise GeneticCodeInitError(
                "code_sequence: %s has length %d, but expected 64"
                % (code_sequence, len(code_sequence))
            )

        self.code_sequence = code_sequence
        self.ID = ID
        self.name = name
        self.start_codon_sequence = start_codon_sequence
        start_codons = {}
        if start_codon_sequence:
            for codon, aa in zip(self._codons, start_codon_sequence):
                if aa != "-":
                    start_codons[codon] = aa
        self.start_codons = start_codons
        codon_lookup = dict(list(zip(self._codons, code_sequence)))
        self.codons = codon_lookup
        # create synonyms for each aa
        aa_lookup = {}
        for codon in self._codons:
            aa = codon_lookup[codon]
            if aa not in aa_lookup:
                aa_lookup[aa] = [codon]
            else:
                aa_lookup[aa].append(codon)
        self.synonyms = aa_lookup
        sense_codons = codon_lookup.copy()
        # create sense codons
        stop_codons = self["*"]
        for c in stop_codons:
            del sense_codons[c]
        self.sense_codons = sense_codons
        # create anticodons
        ac = {}
        for aa, codons in list(self.synonyms.items()):
            ac[aa] = list(map(_simple_rc, codons))
        self.anticodons = ac

    def _analyze_quartet(self, codons, aa):
        """Analyzes a quartet of codons and amino acids: returns list of lists.

        Each list contains one block, splitting at R/Y if necessary.

        codons should be a list of 4 codons.
        aa should be a list of 4 amino acid symbols.

        Possible states:
            - All amino acids are the same: returns list of one quartet.
            - Two groups of 2 aa: returns list of two doublets.
            - One group of 2 and 2 groups of 1: list of one doublet, 2 singles.
            - 4 groups of 1: four singles.

        Note: codon blocks like Ile in the standard code (AUU, AUC, AUA) will
        be split when they cross the R/Y boundary, so [[AUU, AUC], [AUA]]. This
        would also apply to a block like AUC AUA AUG -> [[AUC],[AUA,AUG]],
        although this latter pattern is not observed in the standard code.
        """
        if aa[0] == aa[1]:
            first_doublet = True
        else:
            first_doublet = False
        if aa[2] == aa[3]:
            second_doublet = True
        else:
            second_doublet = False
        if first_doublet and second_doublet and aa[1] == aa[2]:
            return [codons]
        else:
            blocks = []
            if first_doublet:
                blocks.append(codons[:2])
            else:
                blocks.extend([[codons[0]], [codons[1]]])
            if second_doublet:
                blocks.append(codons[2:])
            else:
                blocks.extend([[codons[2]], [codons[3]]])
            return blocks

    def _get_blocks(self):
        """Returns list of lists of codon blocks in the genetic code.

        A codon block can be:
            - a quartet, if all 4 XYn codons have the same amino acid.
            - a doublet, if XYt and XYc or XYa and XYg have the same aa.
            - a singlet, otherwise.

        Returns a list of the quartets, doublets, and singlets in the order
        UUU -> GGG.

        Note that a doublet cannot span the purine/pyrimidine boundary, and
        a quartet cannot span the boundary between two codon blocks whose first
        two bases differ.
        """
        if hasattr(self, "_blocks"):
            return self._blocks
        else:
            blocks = []
            curr_codons = []
            curr_aa = []
            for index, codon, aa in zip(
                list(range(64)), self._codons, self.code_sequence
            ):
                # we're in a new block if it's a new quartet or a different aa
                new_quartet = not index % 4
                if new_quartet and curr_codons:
                    blocks.extend(self._analyze_quartet(curr_codons, curr_aa))
                    curr_codons = []
                    curr_aa = []
                curr_codons.append(codon)
                curr_aa.append(aa)
            # don't forget to append last block
            if curr_codons:
                blocks.extend(self._analyze_quartet(curr_codons, curr_aa))
            self._blocks = blocks
            return self._blocks

    blocks = property(_get_blocks)

    def to_table(self):
        """returns aa to codon mapping as a cogent3 Table"""
        from cogent3.core.moltype import IUPAC_PROTEIN_code_aa

        rows = []
        headers = ["aa", "IUPAC code", "codons"]
        for code, aa in IUPAC_PROTEIN_code_aa.items():
            codons = ",".join(self[code])
            row = [aa, code, codons]
            rows.append(row)
        return Table(header=headers, data=rows, title=self.name)

    def __str__(self):
        """Returns code_sequence that constructs the GeneticCode."""
        return self.code_sequence

    def __repr__(self):
        display = self.to_table()
        return str(display)

    def _repr_html_(self):
        """Returns the html representation of GeneticCode."""
        display = self.to_table()
        display.set_repr_policy(show_shape=False)
        return display._repr_html_()

    def __eq__(self, other):
        """Allows two GeneticCode objects to be compared to each other.
        Two GeneticCode objects are equal if they have equal CodeSequences.
        """
        return str(self) == str(other)

    def __getitem__(self, item):
        """Returns amino acid corresponding to codon, or codons for an aa.

        Returns [] for empty list of codons, 'X' for unknown amino acid.
        """
        item = str(item)
        if len(item) == 1:  # amino acid
            return self.synonyms.get(item, [])
        elif len(item) == 3:  # codon
            key = item.upper()
            key = key.replace("U", "T")
            return self.codons.get(key, "X")
        else:
            raise InvalidCodonError(f"Codon or aa {item} has wrong length")

    def translate(self, dna, start=0):
        """Translates DNA to protein with current GeneticCode.

        Parameters
        ----------
        dna: str
            a string of nucleotides
        start: int
            position to begin translation (used to implement frames)

        Returns
        -------
        String containing amino acid sequence. Translates the entire sequence.
        It is the caller's responsibility to find open reading frames.
        """
        if not dna:
            return ""
        if start + 1 > len(dna):
            raise ValueError("Translation starts after end of RNA")
        return "".join([self[dna[i : i + 3]] for i in range(start, len(dna) - 2, 3)])

    def get_stop_indices(self, dna, start=0):
        """returns indexes for stop codons in the specified frame"""
        stops = self["*"]
        stop_pattern = f"({'|'.join(stops)})"
        stop_pattern = re.compile(stop_pattern)
        seq = str(dna)
        found = [hit.start() for hit in stop_pattern.finditer(seq)]
        found = [index for index in found if index % 3 == start]
        return found

    def sixframes(self, dna):
        """Returns six-frame translation as dict containing {frame:translation}"""
        reverse = dna.rc()
        return [self.translate(dna, start) for start in range(3)] + [
            self.translate(reverse, start) for start in range(3)
        ]

    def is_start(self, codon):
        """Returns True if codon is a start codon, False otherwise."""
        fixed_codon = codon.upper().replace("U", "T")
        return fixed_codon in self.start_codons

    def is_stop(self, codon):
        """Returns True if codon is a stop codon, False otherwise."""
        return self[codon] == "*"

    def changes(self, other):
        """Returns dict of {codon:'XY'} for codons that differ.

        X is the string representation of the amino acid in self, Y is the
        string representation of the amino acid in other. Always returns a
        2-character string.
        """
        changes = {}
        try:
            other_code = other.code_sequence
        except AttributeError:  # try using other directly as sequence
            other_code = other
        for codon, old, new in zip(self._codons, self.code_sequence, other_code):
            if old != new:
                changes[codon] = old + new
        return changes

    def to_regex(self, seq):
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
            if aa in ambigs:
                aa = ambigs[aa]
            else:
                aa = [aa]

            codons = []
            for a in aa:
                codons.extend(self[a])

            # we create a regex non-capturing group for each amino acid
            mappings.append(f"(?:{'|'.join(codons)})")

        return "".join(mappings)


NcbiGeneticCodeData = [
    GeneticCode(*data)
    for data in [
        [
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            1,
            "Standard Nuclear",
            "---M---------------M---------------M----------------------------",
        ],
        [
            "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
            2,
            "Vertebrate Mitochondrial",
            "--------------------------------MMMM---------------M------------",
        ],
        [
            "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            3,
            "Yeast Mitochondrial",
            "----------------------------------MM----------------------------",
        ],
        [
            "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            4,
            "Mold, Protozoan, and Coelenterate Mitochondrial, and Mycoplasma/Spiroplasma Nuclear",
            "--MM---------------M------------MMMM---------------M------------",
        ],
        [
            "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
            5,
            "Invertebrate Mitochondrial",
            "---M----------------------------MMMM---------------M------------",
        ],
        [
            "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            6,
            "Ciliate, Dasycladacean and Hexamita Nuclear",
            "-----------------------------------M----------------------------",
        ],
        [
            "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
            9,
            "Echinoderm and Flatworm Mitochondrial",
            "-----------------------------------M---------------M------------",
        ],
        [
            "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            10,
            "Euplotid Nuclear",
            "-----------------------------------M----------------------------",
        ],
        [
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            11,
            "Bacterial Nuclear and Plant Plastid",
            "---M---------------M------------MMMM---------------M------------",
        ],
        [
            "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            12,
            "Alternative Yeast Nuclear",
            "-------------------M---------------M----------------------------",
        ],
        [
            "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
            13,
            "Ascidian Mitochondrial",
            "-----------------------------------M----------------------------",
        ],
        [
            "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
            14,
            "Alternative Flatworm Mitochondrial",
            "-----------------------------------M----------------------------",
        ],
        [
            "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            15,
            "Blepharisma Nuclear",
            "-----------------------------------M----------------------------",
        ],
        [
            "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            16,
            "Chlorophycean Mitochondrial",
            "-----------------------------------M----------------------------",
        ],
        [
            "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
            20,
            "Trematode Mitochondrial",
            "-----------------------------------M---------------M------------",
        ],
        [
            "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            22,
            "Scenedesmus obliquus Mitochondrial",
            "-----------------------------------M----------------------------",
        ],
        [
            "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            23,
            "Thraustochytrium Mitochondrial",
        ],
    ]
]

# build dict of GeneticCodes keyed by ID (as int, not str)
GeneticCodes = dict([(i.ID, i) for i in NcbiGeneticCodeData])
# add str versions for convenience
for key, value in list(GeneticCodes.items()):
    GeneticCodes[str(key)] = value

DEFAULT = GeneticCodes[1]


def get_code(code_id=1):
    """returns the genetic code

    Parameters
    ----------
    code_id
        genetic code identifier, name, number or string(number), defaults to
        standard genetic code
    """
    code_id = code_id or 1
    if isinstance(code_id, GeneticCode):
        return code_id

    code = None
    if str(code_id).isdigit():
        code = GeneticCodes[code_id]
    else:
        for gc in GeneticCodes.values():
            if gc.name == code_id:
                code = gc

    if code is None:
        raise ValueError(f'No genetic code matching "{code_id}"')

    return code


def available_codes():
    """returns Table listing the available genetic codes"""
    from cogent3.util.table import Table

    all_keys = sorted({int(k) for k in GeneticCodes if str(k).isdigit()})
    rows = [(k, GeneticCodes[k].name) for k in all_keys]
    header = ["Code ID", "Name"]
    return Table(
        header=header,
        data=rows,
        index_name="Code ID",
        title="Specify a genetic code using either 'Name' or "
        "Code ID (as an integer or string)",
    )
