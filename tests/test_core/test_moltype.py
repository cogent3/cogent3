#!/usr/bin/env python
import pickle

from unittest import TestCase, main

from cogent3.core import sequence
from cogent3.core.moltype import (
    DNA,
    PROTEIN,
    RNA,
    STANDARD_CODON,
    AlphabetError,
    AlphabetGroup,
    CodonAlphabet,
    CoreObjectGroup,
    DnaStandardPairs,
    IUPAC_DNA_ambiguities,
    IUPAC_DNA_ambiguities_complements,
    IUPAC_DNA_chars,
    IUPAC_RNA_ambiguities,
    IUPAC_RNA_ambiguities_complements,
    IUPAC_RNA_chars,
    MolType,
    RnaStandardPairs,
    available_moltypes,
    get_moltype,
    make_matches,
    make_pairs,
)
from cogent3.data.molecular_weight import DnaMW, RnaMW


__author__ = "Gavin Huttley, Peter Maxwell, and Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley", "Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

# ind some of the standard alphabets to reduce typing
from numpy.testing import assert_allclose


RnaBases = RNA.alphabets.base
DnaBases = DNA.alphabets.base
AminoAcids = PROTEIN.alphabets.base

# the following classes are to preserve compatibility for older test code
# that assumes mixed-case is OK.
RnaMolType = MolType(
    seq_constructor=sequence.RnaSequence,
    motifset=IUPAC_RNA_chars,
    ambiguities=IUPAC_RNA_ambiguities,
    label="rna_with_lowercase",
    mw_calculator=RnaMW,
    complements=IUPAC_RNA_ambiguities_complements,
    pairs=RnaStandardPairs,
    add_lower=True,
    preserve_existing_moltypes=True,
    make_alphabet_group=True,
)
DnaMolType = MolType(
    seq_constructor=sequence.DnaSequence,
    motifset=IUPAC_DNA_chars,
    ambiguities=IUPAC_DNA_ambiguities,
    label="dna_with_lowercase",
    mw_calculator=DnaMW,
    complements=IUPAC_DNA_ambiguities_complements,
    pairs=DnaStandardPairs,
    add_lower=True,
    preserve_existing_moltypes=True,
    make_alphabet_group=True,
)
ProteinMolType = PROTEIN


class make_matches_tests(TestCase):
    """Tests of the make_matches top-level function"""

    def test_init_empty(self):
        """make_matches should init ok with no parameters"""
        self.assertEqual(make_matches(), {})

    def test_init_monomers(self):
        """make_matches with only monomers should produce {(i,i):True}"""
        m = make_matches("")
        self.assertEqual(m, {})
        m = make_matches("qaz")
        self.assertEqual(m, {("q", "q"): True, ("a", "a"): True, ("z", "z"): True})

    def test_init_gaps(self):
        """make_matches with only gaps should match all gaps to each other"""
        m = make_matches("", "~!")
        self.assertEqual(
            m, {("~", "~"): True, ("!", "!"): True, ("!", "~"): True, ("~", "!"): True}
        )

    def test_init_degen(self):
        """make_matches with only degen should work as expected"""
        m = make_matches(None, None, {"x": "ab", "y": "bc", "z": "cd", "n": "bcd"})
        self.assertEqual(
            m,
            {
                ("x", "x"): False,
                ("x", "y"): False,
                ("x", "n"): False,
                ("y", "x"): False,
                ("y", "y"): False,
                ("y", "z"): False,
                ("y", "n"): False,
                ("z", "y"): False,
                ("z", "z"): False,
                ("z", "n"): False,
                ("n", "x"): False,
                ("n", "y"): False,
                ("n", "z"): False,
                ("n", "n"): False,
            },
        )
        self.assertNotIn(("x", "z"), m)

    def test_init_all(self):
        """make_matches with everything should produce correct dict"""
        m = make_matches(
            "ABC", ("-", "~"), {"X": "AB", "Y": ("B", "C"), "N": list("ABC")}
        )
        exp = {
            ("-", "-"): True,
            ("~", "~"): True,
            ("-", "~"): True,
            ("~", "-"): True,
            ("A", "A"): True,
            ("B", "B"): True,
            ("C", "C"): True,
            ("A", "X"): False,
            ("X", "A"): False,
            ("B", "X"): False,
            ("X", "B"): False,
            ("B", "Y"): False,
            ("Y", "B"): False,
            ("C", "Y"): False,
            ("Y", "C"): False,
            ("A", "N"): False,
            ("N", "A"): False,
            ("B", "N"): False,
            ("N", "B"): False,
            ("C", "N"): False,
            ("N", "C"): False,
            ("X", "X"): False,
            ("Y", "Y"): False,
            ("N", "N"): False,
            ("X", "Y"): False,
            ("Y", "X"): False,
            ("X", "N"): False,
            ("N", "X"): False,
            ("Y", "N"): False,
            ("N", "Y"): False,
        }
        self.assertEqual(m, exp)


class make_pairs_tests(TestCase):
    """Tests of the top-level make_pairs factory function."""

    def setUp(self):
        """Define some standard pairs and other data"""
        self.pairs = {("U", "A"): True, ("A", "U"): True, ("G", "U"): False}

    def test_init_empty(self):
        """make_pairs should init ok with no parameters"""
        self.assertEqual(make_pairs(), {})

    def test_init_pairs(self):
        """make_pairs with just pairs should equal the original"""
        self.assertEqual(make_pairs(self.pairs), self.pairs)
        self.assertIsNot(make_pairs(self.pairs), self.pairs)

    def test_init_monomers(self):
        """make_pairs with pairs and monomers should equal just the pairs"""
        self.assertEqual(make_pairs(self.pairs, "ABCDEFG"), self.pairs)
        self.assertIsNot(make_pairs(self.pairs, "ABCDEFG"), self.pairs)

    def test_init_gaps(self):
        """make_pairs should add all combinations of gaps as weak pairs"""
        p = make_pairs(self.pairs, None, "-~")
        self.assertNotEqual(p, self.pairs)
        self.pairs.update(
            {("~", "~"): False, ("-", "~"): False, ("-", "-"): False, ("~", "-"): False}
        )
        self.assertEqual(p, self.pairs)

    def test_init_degen(self):
        """make_pairs should add in degenerate combinations as weak pairs"""
        p = make_pairs(self.pairs, "AUG", "-", {"R": "AG", "Y": "CU", "W": "AU"})
        self.assertNotEqual(p, self.pairs)
        self.pairs.update(
            {
                ("-", "-"): False,
                ("A", "Y"): False,
                ("Y", "A"): False,
                ("A", "W"): False,
                ("W", "A"): False,
                ("U", "R"): False,
                ("R", "U"): False,
                ("U", "W"): False,
                ("W", "U"): False,
                ("G", "Y"): False,
                ("G", "W"): False,
                ("R", "Y"): False,
                ("R", "W"): False,
                ("Y", "R"): False,
                ("Y", "W"): False,
                ("W", "R"): False,
                ("W", "Y"): False,
                ("W", "W"): False,
            }
        )
        self.assertEqual(p, self.pairs)


class CoreObjectGroupTests(TestCase):
    """Tests of the CoreObjectGroup class."""

    def test_init(self):
        """CoreObjectGroup should init with basic list of objects."""

        class o(object):
            def __init__(self, s):
                self.s = s

        base = o("base")
        c = CoreObjectGroup(base)
        self.assertIs(c.base, base)
        self.assertIs(c.degen, None)
        self.assertIs(c.base.degen, None)

        base, degen, gap, degengap = list(map(o, ["base", "degen", "gap", "degengap"]))
        c = CoreObjectGroup(base, degen, gap, degengap)
        self.assertIs(c.base, base)
        self.assertIs(c.base.degen, degen)
        self.assertIs(c.degen.gapped, degengap)


class AlphabetGroupTests(TestCase):
    """Tests of the AlphabetGroup class."""

    def test_init(self):
        """AlphabetGroup should initialize successfully"""
        chars = "AB"
        degens = {"C": "AB"}
        g = AlphabetGroup(chars, degens)
        self.assertEqual("".join(g.base), "AB")
        self.assertEqual("".join(g.degen), "ABC")
        self.assertEqual("".join(g.gapped), "AB-")
        self.assertEqual("".join(g.degen_gapped), "AB-C?")


class MolTypeTests(TestCase):
    """Tests of the MolType class. Should support same API as old Alphabet."""

    def test_to_regex(self):
        """returns a valid regex"""
        seq = "ACYGR"
        regular_expression = DNA.to_regex(seq=seq)
        self.assertEqual(regular_expression, "AC[CT]G[AG]")
        # raises an exception if a string is already a regex, or invalid
        with self.assertRaises(ValueError):
            DNA.to_regex("(?:GAT|GAC)(?:GGT|GGC|GGA|GGG)(?:GAT|GAC)(?:CAA|CAG)")

    def test_pickling(self):
        pkl = pickle.dumps(DNA)
        got = pickle.loads(pkl)
        self.assertIsInstance(got, type(DNA))
        self.assertEqual(DNA.ambiguities, got.ambiguities)

    def test_get_moltype(self):
        """correctly return a moltype by name"""
        for label in ("dna", "rna", "protein", "protein_with_stop"):
            mt = get_moltype(label)
            self.assertEqual(mt.label, label)
            mt = get_moltype(label.upper())
            self.assertEqual(mt.label, label)

        mt = get_moltype(DNA)
        self.assertEqual(mt.label, "dna")
        with self.assertRaises(ValueError):
            _ = get_moltype("blah")

    def test_available_moltypes(self):
        """available_moltypes returns table with correct columns and rows"""
        available = available_moltypes()
        self.assertEqual(available.shape, (7, 3))
        self.assertEqual(available[1, "Number of states"], 4)
        self.assertEqual(available["dna", "Number of states"], 4)
        txt = repr(available)
        self.assertIn("'dna'", txt)

    def test_init_minimal(self):
        """MolType should init OK with just monomers"""
        a = MolType("Abc")
        self.assertIn("A", a.alphabet)
        self.assertNotIn("a", a.alphabet)  # case-sensitive
        self.assertIn("b", a.alphabet)
        self.assertNotIn("B", a.alphabet)
        self.assertNotIn("x", a.alphabet)

    def test_init_everything(self):
        """MolType should init OK with all parameters set"""
        k = dict.fromkeys
        a = MolType(
            k("Abc"),
            ambiguities={"d": "bc"},
            gaps=k("~"),
            complements={"b": "c", "c": "b"},
            pairs={},
            add_lower=False,
        )
        for i in "Abcd~":
            self.assertIn(i, a)
        self.assertEqual(a.complement("b"), "c")
        self.assertEqual(a.complement("AbcAA"), "AcbAA")
        self.assertEqual(a.first_degenerate("AbcdA"), 3)
        self.assertEqual(a.first_gap("a~c"), 1)
        self.assertEqual(a.first_invalid("Abcx"), 3)

    def test_strip_degenerate(self):
        """MolType strip_degenerate should remove any degenerate bases"""
        s = RnaMolType.strip_degenerate
        self.assertEqual(s("UCAG-"), "UCAG-")
        self.assertEqual(s("NRYSW"), "")
        self.assertEqual(s("USNG"), "UG")

    def test_strip_bad(self):
        """MolType strip_bad should remove any non-base, non-gap chars"""
        s = RnaMolType.strip_bad
        self.assertEqual(s("UCAGwsnyrHBND-D"), "UCAGwsnyrHBND-D")
        self.assertEqual(s("@#^*($@!#&()!@QZX"), "")
        self.assertEqual(s("aaa ggg ---!ccc"), "aaaggg---ccc")

    def test_strip_bad_and_gaps(self):
        """MolType strip_bad_and_gaps should remove gaps and bad chars"""
        s = RnaMolType.strip_bad_and_gaps
        self.assertEqual(s("UCAGwsnyrHBND-D"), "UCAGwsnyrHBNDD")
        self.assertEqual(s("@#^*($@!#&()!@QZX"), "")
        self.assertEqual(s("aaa ggg ---!ccc"), "aaagggccc")

    def test_complement(self):
        """MolType complement should correctly complement sequence"""
        self.assertEqual(RnaMolType.complement("UauCG-NR"), "AuaGC-NY")
        self.assertEqual(DnaMolType.complement("TatCG-NR"), "AtaGC-NY")
        self.assertEqual(RnaMolType.complement(""), "")
        self.assertRaises(TypeError, ProteinMolType.complement, "ACD")
        # if it wasn't a string, result should be a list
        self.assertEqual(RnaMolType.complement(list("UauCG-NR")), list("AuaGC-NY"))
        self.assertEqual(RnaMolType.complement(("a", "c")), ("u", "g"))
        # constructor should fail for a dict
        self.assertRaises(ValueError, RnaMolType.complement, {"a": "c"})

    def test_rc(self):
        """MolType rc should correctly reverse-complement sequence"""
        self.assertEqual(RnaMolType.rc("U"), "A")
        self.assertEqual(RnaMolType.rc(""), "")
        self.assertEqual(RnaMolType.rc("R"), "Y")
        self.assertEqual(RnaMolType.rc("UauCG-NR"), "YN-CGauA")
        self.assertEqual(DnaMolType.rc("TatCG-NR"), "YN-CGatA")
        self.assertRaises(TypeError, ProteinMolType.rc, "ACD")
        # if it wasn't a string, result should be a list
        self.assertEqual(RnaMolType.rc(list("UauCG-NR")), list("YN-CGauA"))
        self.assertEqual(RnaMolType.rc(("a", "c")), ("g", "u"))
        # constructor should fail for a dict
        self.assertRaises(ValueError, RnaMolType.rc, {"a": "c"})

    def test_contains(self):
        """MolType contains should return correct result"""
        for i in "UCAGWSMKRYBDHVN-" + "UCAGWSMKRYBDHVN-".lower():
            self.assertIn(i, RnaMolType)
        for i in "x!@#$%^&ZzQq":
            self.assertNotIn(i, RnaMolType)

        a = MolType(dict.fromkeys("ABC"), add_lower=True)
        for i in "abcABC":
            self.assertIn(i, a)
        self.assertNotIn("x", a)
        b = MolType(dict.fromkeys("ABC"), add_lower=False)
        for i in "ABC":
            self.assertIn(i, b)
        for i in "abc":
            self.assertNotIn(i, b)

    def test_iter(self):
        """MolType iter should iterate over monomer order"""
        self.assertEqual(list(RnaMolType), ["U", "C", "A", "G", "u", "c", "a", "g"])
        a = MolType("ZXCV")
        self.assertEqual(list(a), ["Z", "X", "C", "V"])

    def test_is_gapped(self):
        """MolType is_gapped should return True if gaps in seq"""
        g = RnaMolType.is_gapped
        self.assertFalse(g(""))
        self.assertFalse(g("ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN"))
        self.assertTrue(g("-"))
        self.assertTrue(g("--"))
        self.assertTrue(g("CAGUCGUACGUCAGUACGUacucauacgac-caguACUG"))
        self.assertTrue(g("CA--CGUAUGCA-----g"))
        self.assertTrue(g("CAGU-"))

    def test_is_gap(self):
        """MolType is_gap should return True if char is a gap"""
        g = RnaMolType.is_gap
        # True for the empty string
        self.assertFalse(g(""))
        # True for all the standard and degenerate symbols
        s = "ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN"
        self.assertFalse(g(s))
        for i in s:
            self.assertFalse(g(i))
        # should be true for a single gap
        self.assertTrue(g("-"))
        # note that it _shouldn't_ be true for a run of gaps: use a.is_gapped()
        self.assertFalse(g("--"))

    def test_is_degenerate(self):
        """MolType is_degenerate should return True if degen symbol in seq"""
        d = RnaMolType.is_degenerate
        self.assertFalse(d(""))
        self.assertFalse(d("UACGCUACAUGuacgucaguGCUAGCUA---ACGUCAG"))
        self.assertTrue(d("N"))
        self.assertTrue(d("R"))
        self.assertTrue(d("y"))
        self.assertTrue(d("GCAUguagcucgUCAGUCAGUACgUgcasCUAG"))
        self.assertTrue(d("ACGYAUGCUGYEWEWNFMNfuwbybcwuybcjwbeiwfub"))

    def test_is_valid(self):
        """MolType is_valid should return True if any unknown symbol in seq"""
        v = RnaMolType.is_valid
        self.assertFalse(v(3))
        self.assertFalse(v(None))
        self.assertTrue(v("ACGUGCAUGUCAYCAYGUACGcaugacyugc----RYNCYRNC"))
        self.assertTrue(v(""))
        self.assertTrue(v("a"))
        self.assertFalse(v("ACIUBHFWUIXZKLNJUCIHBICNSOWMOINJ"))
        self.assertFalse(v("CAGUCAGUCACA---GACCAUG-_--cgau"))

    def test_is_strict(self):
        """MolType is_strict should return True if all symbols in Monomers"""
        s = RnaMolType.is_strict
        self.assertFalse(s(3))
        self.assertFalse(s(None))
        self.assertTrue(s(""))
        self.assertTrue(s("A"))
        self.assertTrue(s("UAGCACUgcaugcauGCAUGACuacguACAUG"))
        self.assertFalse(s("CAGUCGAUCA-cgaucagUCGAUGAC"))
        self.assertFalse(s("ACGUGCAUXCAGUCAG"))

    def test_first_gap(self):
        """MolType first_gap should return index of first gap symbol, or None"""
        g = RnaMolType.first_gap
        self.assertEqual(g(""), None)
        self.assertEqual(g("a"), None)
        self.assertEqual(g("uhacucHuhacUIUIhacan"), None)
        self.assertEqual(g("-abc"), 0)
        self.assertEqual(g("b-ac"), 1)
        self.assertEqual(g("abcd-"), 4)

    def test_first_degenerate(self):
        """MolType first_degenerate should return index of first degen symbol"""
        d = RnaMolType.first_degenerate
        self.assertEqual(d(""), None)
        self.assertEqual(d("a"), None)
        self.assertEqual(d("UCGACA--CU-gacucaguacgua"), None)
        self.assertEqual(d("nCAGU"), 0)
        self.assertEqual(d("CUGguagvAUG"), 7)
        self.assertEqual(d("ACUGCUAacgud"), 11)

    def test_first_invalid(self):
        """MolType first_invalid should return index of first invalid symbol"""
        i = RnaMolType.first_invalid
        self.assertEqual(i(""), None)
        self.assertEqual(i("A"), None)
        self.assertEqual(i("ACGUNVBuacg-wskmWSMKYRryNnN--"), None)
        self.assertEqual(i("x"), 0)
        self.assertEqual(i("rx"), 1)
        self.assertEqual(i("CAGUNacgunRYWSwx"), 15)

    def test_first_non_strict(self):
        """MolType first_non_strict should return index of first non-strict symbol"""
        s = RnaMolType.first_non_strict
        self.assertEqual(s(""), None)
        self.assertEqual(s("A"), None)
        self.assertEqual(s("ACGUACGUcgaucagu"), None)
        self.assertEqual(s("N"), 0)
        self.assertEqual(s("-"), 0)
        self.assertEqual(s("x"), 0)
        self.assertEqual(s("ACGUcgAUGUGCAUcaguX"), 18)
        self.assertEqual(s("ACGUcgAUGUGCAUcaguX-38243829"), 18)

    def test_disambiguate(self):
        """MolType disambiguate should remove degenerate bases"""
        d = RnaMolType.disambiguate
        self.assertEqual(d(""), "")
        self.assertEqual(d("AGCUGAUGUA--CAGU"), "AGCUGAUGUA--CAGU")
        self.assertEqual(d("AUn-yrs-wkmCGwmrNMWRKY", "strip"), "AU--CG")
        self.assertEqual(d(tuple("AUn-yrs-wkmCGwmrNMWRKY"), "strip"), tuple("AU--CG"))
        s = "AUn-yrs-wkmCGwmrNMWRKY"
        t = d(s, "random")
        u = d(s, "random")
        for i, j in zip(s, t):
            if i in RnaMolType.degenerates:
                self.assertIn(j, RnaMolType.degenerates[i])
            else:
                self.assertEqual(i, j)
        self.assertNotEqual(t, u)
        self.assertEqual(d(tuple("UCAG"), "random"), tuple("UCAG"))
        self.assertEqual(len(s), len(t))
        self.assertIs(RnaMolType.first_degenerate(t), None)
        # should raise exception on unknown disambiguation method
        self.assertRaises(NotImplementedError, d, s, "xyz")

    def test_degap(self):
        """MolType degap should remove all gaps from sequence"""
        g = RnaMolType.degap
        self.assertEqual(g(""), "")
        self.assertEqual(g("GUCAGUCgcaugcnvuincdks"), "GUCAGUCgcaugcnvuincdks")
        self.assertEqual(g("----------------"), "")
        self.assertEqual(g("gcuauacg-"), "gcuauacg")
        self.assertEqual(g("-CUAGUCA"), "CUAGUCA")
        self.assertEqual(g("---a---c---u----g---"), "acug")
        self.assertEqual(g(tuple("---a---c---u----g---")), tuple("acug"))

    def test_gap_indices(self):
        """MolType gap_indices should return correct gap positions"""
        g = RnaMolType.gap_indices
        self.assertEqual(g(""), [])
        self.assertEqual(g("ACUGUCAGUACGHFSDKJCUICDNINS"), [])
        self.assertEqual(g("GUACGUIACAKJDC-SDFHJDSFK"), [14])
        self.assertEqual(g("-DSHFUHDSF"), [0])
        self.assertEqual(g("UACHASJAIDS-"), [11])
        self.assertEqual(
            g("---CGAUgCAU---ACGHc---ACGUCAGU---"),
            [0, 1, 2, 11, 12, 13, 19, 20, 21, 30, 31, 32],
        )
        a = MolType({"A": 1}, gaps=dict.fromkeys("!@#$%"))
        g = a.gap_indices
        self.assertEqual(g(""), [])
        self.assertEqual(g("!!!"), [0, 1, 2])
        self.assertEqual(g("!@#$!@#$!@#$"), list(range(12)))
        self.assertEqual(g("cguua!cgcuagua@cguasguadc#"), [5, 14, 25])

    def test_gap_vector(self):
        """MolType gap_vector should return correct gap positions"""
        g = RnaMolType.gap_vector
        self.assertEqual(g(""), [])
        self.assertEqual(g("ACUGUCAGUACGHFSDKJCUICDNINS"), [False] * 27)
        self.assertEqual(
            g("GUACGUIACAKJDC-SDFHJDSFK"),
            list(map(bool, list(map(int, "000000000000001000000000")))),
        )
        self.assertEqual(g("-DSHFUHDSF"), list(map(bool, list(map(int, "1000000000")))))
        self.assertEqual(
            g("UACHASJAIDS-"), list(map(bool, list(map(int, "000000000001"))))
        )
        self.assertEqual(
            g("---CGAUgCAU---ACGHc---ACGUCAGU---"),
            list(map(bool, list(map(int, "111000000001110000011100000000111")))),
        )
        a = MolType({"A": 1}, gaps=dict.fromkeys("!@#$%"))
        g = a.gap_vector
        self.assertEqual(g(""), [])
        self.assertEqual(g("!!!"), list(map(bool, [1, 1, 1])))
        self.assertEqual(g("!@#$!@#$!@#$"), [True] * 12)
        self.assertEqual(
            g("cguua!cgcuagua@cguasguadc#"),
            list(map(bool, list(map(int, "00000100000000100000000001")))),
        )

    def test_gap_maps(self):
        """MolType gap_maps should return dicts mapping gapped/ungapped pos"""
        empty = ""
        no_gaps = "aaa"
        all_gaps = "---"
        start_gaps = "--abc"
        end_gaps = "ab---"
        mid_gaps = "--a--b-cd---"
        gm = RnaMolType.gap_maps
        self.assertEqual(gm(empty), ({}, {}))
        self.assertEqual(gm(no_gaps), ({0: 0, 1: 1, 2: 2}, {0: 0, 1: 1, 2: 2}))
        self.assertEqual(gm(all_gaps), ({}, {}))
        self.assertEqual(gm(start_gaps), ({0: 2, 1: 3, 2: 4}, {2: 0, 3: 1, 4: 2}))
        self.assertEqual(gm(end_gaps), ({0: 0, 1: 1}, {0: 0, 1: 1}))
        self.assertEqual(
            gm(mid_gaps), ({0: 2, 1: 5, 2: 7, 3: 8}, {2: 0, 5: 1, 7: 2, 8: 3})
        )

    def test_count_gaps(self):
        """MolType count_gaps should return correct gap count"""
        c = RnaMolType.count_gaps
        self.assertEqual(c(""), 0)
        self.assertEqual(c("ACUGUCAGUACGHFSDKJCUICDNINS"), 0)
        self.assertEqual(c("GUACGUIACAKJDC-SDFHJDSFK"), 1)
        self.assertEqual(c("-DSHFUHDSF"), 1)
        self.assertEqual(c("UACHASJAIDS-"), 1)
        self.assertEqual(c("---CGAUgCAU---ACGHc---ACGUCAGU---"), 12)
        a = MolType({"A": 1}, gaps=dict.fromkeys("!@#$%"))
        c = a.count_gaps
        self.assertEqual(c(""), 0)
        self.assertEqual(c("!!!"), 3)
        self.assertEqual(c("!@#$!@#$!@#$"), 12)
        self.assertEqual(c("cguua!cgcuagua@cguasguadc#"), 3)

    def test_count_degenerate(self):
        """MolType count_degenerate should return correct degen base count"""
        d = RnaMolType.count_degenerate
        self.assertEqual(d(""), 0)
        self.assertEqual(d("GACUGCAUGCAUCGUACGUCAGUACCGA"), 0)
        self.assertEqual(d("N"), 1)
        self.assertEqual(d("NRY"), 3)
        self.assertEqual(d("ACGUAVCUAGCAUNUCAGUCAGyUACGUCAGS"), 4)

    def test_possibilites(self):
        """MolType possibilities should return correct # possible sequences"""
        p = RnaMolType.possibilities
        self.assertEqual(p(""), 1)
        self.assertEqual(p("ACGUgcaucagUCGuGCAU"), 1)
        self.assertEqual(p("N"), 4)
        self.assertEqual(p("R"), 2)
        self.assertEqual(p("H"), 3)
        self.assertEqual(p("nRh"), 24)
        self.assertEqual(p("AUGCnGUCAg-aurGauc--gauhcgauacgws"), 96)

    def test_MW(self):
        """MolType MW should return correct molecular weight"""
        r = RnaMolType.mw
        p = ProteinMolType.mw
        self.assertEqual(p(""), 0)
        self.assertEqual(r(""), 0)
        assert_allclose(p("A"), 89.09)
        assert_allclose(r("A"), 375.17)
        assert_allclose(p("AAA"), 231.27)
        assert_allclose(r("AAA"), 1001.59)
        assert_allclose(r("AAACCCA"), 2182.37)

    def test_can_match(self):
        """MolType can_match should return True if all positions can match"""
        m = RnaMolType.can_match
        self.assertTrue(m("", ""))
        self.assertTrue(m("UCAG", "UCAG"))
        self.assertFalse(m("UCAG", "ucag"))
        self.assertTrue(m("UCAG", "NNNN"))
        self.assertTrue(m("NNNN", "UCAG"))
        self.assertTrue(m("NNNN", "NNNN"))
        self.assertFalse(m("N", "x"))
        self.assertFalse(m("N", "-"))
        self.assertTrue(m("UCAG", "YYRR"))
        self.assertTrue(m("UCAG", "KMWS"))

    def test_can_mismatch(self):
        """MolType can_mismatch should return True on any possible mismatch"""
        m = RnaMolType.can_mismatch
        self.assertFalse(m("", ""))
        self.assertTrue(m("N", "N"))
        self.assertTrue(m("R", "R"))
        self.assertTrue(m("N", "r"))
        self.assertTrue(m("CGUACGCAN", "CGUACGCAN"))
        self.assertTrue(m("U", "C"))
        self.assertTrue(m("UUU", "UUC"))
        self.assertTrue(m("UUU", "UUY"))
        self.assertFalse(m("UUU", "UUU"))
        self.assertFalse(m("UCAG", "UCAG"))
        self.assertFalse(m("U--", "U--"))

    def test_must_match(self):
        """MolType must_match should return True when no possible mismatches"""
        m = RnaMolType.must_match
        self.assertTrue(m("", ""))
        self.assertFalse(m("N", "N"))
        self.assertFalse(m("R", "R"))
        self.assertFalse(m("N", "r"))
        self.assertFalse(m("CGUACGCAN", "CGUACGCAN"))
        self.assertFalse(m("U", "C"))
        self.assertFalse(m("UUU", "UUC"))
        self.assertFalse(m("UUU", "UUY"))
        self.assertTrue(m("UU-", "UU-"))
        self.assertTrue(m("UCAG", "UCAG"))

    def test_can_pair(self):
        """MolType can_pair should return True if all positions can pair"""
        p = RnaMolType.can_pair
        self.assertTrue(p("", ""))
        self.assertFalse(p("UCAG", "UCAG"))
        self.assertTrue(p("UCAG", "CUGA"))
        self.assertFalse(p("UCAG", "cuga"))
        self.assertTrue(p("UCAG", "NNNN"))
        self.assertTrue(p("NNNN", "UCAG"))
        self.assertTrue(p("NNNN", "NNNN"))
        self.assertFalse(p("N", "x"))
        self.assertFalse(p("N", "-"))
        self.assertTrue(p("-", "-"))
        self.assertTrue(p("UCAGU", "KYYRR"))
        self.assertTrue(p("UCAG", "KKRS"))
        self.assertTrue(p("U", "G"))

        d = DnaMolType.can_pair
        self.assertFalse(d("T", "G"))

    def test_can_mispair(self):
        """MolType can_mispair should return True on any possible mispair"""
        m = RnaMolType.can_mispair
        self.assertFalse(m("", ""))
        self.assertTrue(m("N", "N"))
        self.assertTrue(m("R", "Y"))
        self.assertTrue(m("N", "r"))
        self.assertTrue(m("CGUACGCAN", "NUHCHUACH"))
        self.assertTrue(m("U", "C"))
        self.assertTrue(m("U", "R"))
        self.assertTrue(m("UUU", "AAR"))
        self.assertTrue(m("UUU", "GAG"))
        self.assertFalse(m("UUU", "AAA"))
        self.assertFalse(m("UCAG", "CUGA"))
        self.assertTrue(m("U--", "--U"))

        d = DnaMolType.can_pair
        self.assertTrue(d("TCCAAAGRYY", "RRYCTTTGGA"))

    def test_must_pair(self):
        """MolType must_pair should return True when no possible mispairs"""
        m = RnaMolType.must_pair
        self.assertTrue(m("", ""))
        self.assertFalse(m("N", "N"))
        self.assertFalse(m("R", "Y"))
        self.assertFalse(m("A", "A"))
        self.assertFalse(m("CGUACGCAN", "NUGCGUACG"))
        self.assertFalse(m("U", "C"))
        self.assertFalse(m("UUU", "AAR"))
        self.assertFalse(m("UUU", "RAA"))
        self.assertFalse(m("UU-", "-AA"))
        self.assertTrue(m("UCAG", "CUGA"))

        d = DnaMolType.must_pair
        self.assertTrue(d("TCCAGGG", "CCCTGGA"))
        self.assertTrue(d("tccaggg", "ccctgga"))
        self.assertFalse(d("TCCAGGG", "NCCTGGA"))


class RnaMolTypeTests(TestCase):
    """Spot-checks of alphabet functionality applied to RNA alphabet."""

    def test_contains(self):
        """RnaMolType should __contain__ the expected symbols."""
        keys = "ucagrymkwsbhvdn?-"
        for k in keys:
            self.assertIn(k, RnaMolType)
        for k in keys.upper():
            self.assertIn(k, RnaMolType)
        self.assertNotIn("X", RnaMolType)

    def test_degenerate_from_seq(self):
        """RnaMolType degenerate_from_seq should give correct results"""
        d = RnaMolType.degenerate_from_seq
        # check monomers
        self.assertEqual(d("a"), "a")
        self.assertEqual(d("C"), "C")
        # check seq of monomers
        self.assertEqual(d("aaaaa"), "a")
        # check some 2- to 4-way cases
        self.assertEqual(d("aaaaag"), "r")
        self.assertEqual(d("ccgcgcgcggcc"), "s")
        self.assertEqual(d("accgcgcgcggcc"), "v")
        self.assertEqual(d("aaaaagcuuu"), "n")
        # check some cases with gaps
        self.assertEqual(d("aa---aaagcuuu"), "?")
        self.assertEqual(d("aaaaaaaaaaaaaaa-"), "?")
        self.assertEqual(d("----------------"), "-")
        # check mixed case example
        self.assertEqual(d("AaAAaa"), "A")
        # check example with degenerate symbols in set
        self.assertEqual(d("RS"), "V")
        self.assertEqual(d("RN-"), "?")
        # check that it works for proteins as well
        p = ProteinMolType.degenerate_from_seq
        self.assertEqual(p("A"), "A")
        self.assertEqual(p("AAA"), "A")
        self.assertEqual(p("DN"), "B")
        self.assertEqual(p("---"), "-")
        self.assertEqual(p("ACD"), "X")
        self.assertEqual(p("ABD"), "X")
        self.assertEqual(p("ACX"), "X")
        self.assertEqual(p("AC-"), "?")

    def test_get_css_styles(self):
        """generates css style mapping and css"""
        css, styles = DNA.get_css_style()
        # styles should consist of core characters appended by DNA.label
        states = list(DNA) + ["ambig", "terminal_ambig"]
        got = {styles[b] for b in states}
        expect = {b + "_dna" for b in states}
        self.assertEqual(got, expect)
        for state in expect:
            self.assertTrue(state in "\n".join(css))

        # check subset of protein
        css, styles = PROTEIN.get_css_style()
        states = list(PROTEIN) + ["ambig", "terminal_ambig"]
        got = {styles[b] for b in states}
        expect = {b + "_protein" for b in states}
        self.assertEqual(got, expect)

    def test_get_css_no_label(self):
        """should not fail when moltype has no label"""
        dna = get_moltype("dna")
        orig_label, dna.label = dna.label, None
        _ = dna.get_css_style()
        dna.label = orig_label


class _AlphabetTestCase(TestCase):
    def assertEqualSeqs(self, a, b):
        """For when we don't care about the type, just the elements"""
        self.assertEqual(list(a), list(b))

    def assertEqualSets(self, a, b):
        self.assertEqual(set(a), set(b))


class DNAAlphabet(_AlphabetTestCase):
    def setUp(self):
        self.alpha = DNA.alphabet

    def test_exclude(self):
        """Nucleotide alphabet testing excluding gap motif"""
        self.assertEqualSeqs(self.alpha, ["T", "C", "A", "G"])

    def test_include(self):
        """Nucleotide alphabet testing including gap motif"""
        self.assertEqualSets(self.alpha.with_gap_motif(), ["A", "C", "G", "T", "-"])

    def test_usesubset(self):
        """testing using a subset of motifs."""
        self.assertEqualSets(self.alpha.with_gap_motif(), ["A", "C", "G", "T", "-"])
        alpha = self.alpha.get_subset(motif_subset=["A"])
        self.assertEqualSets(alpha, ["A"])
        # self.assertRaises(AlphabetError, self.alpha.get_subset, ['A','C'])
        alpha = DNA.alphabet
        self.assertEqualSets(alpha, ["T", "C", "A", "G"])
        alpha = alpha.get_subset(motif_subset=["A", "T", "G"])
        self.assertEqualSets(alpha, ["A", "G", "T"])


class DinucAlphabet(_AlphabetTestCase):
    def setUp(self):
        self.alpha = DNA.alphabet.with_gap_motif().get_word_alphabet(2)

    def test_exclude(self):
        """Dinucleotide alphabet testing excluding gap motif"""
        expected = [
            "-A",
            "-C",
            "-G",
            "-T",
            "A-",
            "AA",
            "AC",
            "AG",
            "AT",
            "C-",
            "CA",
            "CC",
            "CG",
            "CT",
            "G-",
            "GA",
            "GC",
            "GG",
            "GT",
            "T-",
            "TA",
            "TC",
            "TG",
            "TT",
        ]

        self.assertEqualSets(self.alpha.get_subset(["--"], excluded=True), expected)

    def test_include(self):
        """Dinucleotide alphabet testing including gap motif"""
        expected = [
            "--",
            "-A",
            "-C",
            "-G",
            "-T",
            "A-",
            "AA",
            "AC",
            "AG",
            "AT",
            "C-",
            "CA",
            "CC",
            "CG",
            "CT",
            "G-",
            "GA",
            "GC",
            "GG",
            "GT",
            "T-",
            "TA",
            "TC",
            "TG",
            "TT",
        ]
        self.assertEqualSets(self.alpha, expected)

    def test_usesubset(self):
        """testing using a subset of motifs."""
        alpha = self.alpha.get_subset(motif_subset=["AA", "CA", "GT"])
        self.assertEqualSeqs(alpha, ["AA", "CA", "GT"])

        self.assertRaises(
            AlphabetError, alpha.get_subset, motif_subset=["AA", "CA", "GT", "TT"]
        )

    def test_usesubsetbyfreq(self):
        """testing using a subset of motifs by using motif probs."""
        motif_freqs = {
            "--": 0,
            "-A": 0.0,
            "-C": 0,
            "-G": 0,
            "-T": 0,
            "A-": 0,
            "AA": 1.0,
            "AC": 0.0,
            "AG": 0,
            "AT": 0,
            "C-": 0,
            "CA": 1,
            "CC": 0,
            "CG": 0,
            "CT": 0,
            "G-": 0,
            "GA": 0,
            "GC": 0,
            "GG": 0,
            "GT": 1,
            "T-": 0,
            "TA": 0,
            "TC": 0,
            "TG": 0,
            "TT": 0,
        }

        alpha = self.alpha.get_subset(motif_freqs)
        self.assertEqualSets(alpha, ["AA", "CA", "GT"])

    def test_strand_symmetric_motifs(self):
        """construction of strand symmetric motif sets"""
        # fails for a moltype with no strand complement
        with self.assertRaises(TypeError):
            PROTEIN.strand_symmetric_motifs()

        got = DNA.strand_symmetric_motifs(motif_length=1)
        expect = set([("A", "T"), ("C", "G")])
        self.assertEqual(got, expect)
        got = RNA.strand_symmetric_motifs(motif_length=1)
        expect = set([("A", "U"), ("C", "G")])
        self.assertEqual(got, expect)
        got = DNA.strand_symmetric_motifs(motif_length=2)
        self.assertEqual(len(got), 8)
        got = DNA.strand_symmetric_motifs(motif_length=3)
        self.assertEqual(len(got), 32)


class TestCodonAlphabet(_AlphabetTestCase):
    def setUp(self):
        self.alpha = STANDARD_CODON

    def test_ambiguous_gaps(self):
        alpha = self.alpha.with_gap_motif()
        self.assertEqual(len(alpha.resolve_ambiguity("AT?")), 4)
        self.assertRaises(Exception, alpha.resolve_ambiguity, "at-")
        self.assertEqual(len(alpha.resolve_ambiguity("???")), 62)
        self.assertEqual(len(alpha.resolve_ambiguity("---")), 1)

        alpha = self.alpha
        self.assertEqual(len(alpha.resolve_ambiguity("AT?")), 4)
        self.assertRaises(Exception, alpha.resolve_ambiguity, "at-")
        self.assertEqual(len(alpha.resolve_ambiguity("???")), 61)
        self.assertRaises(Exception, alpha.resolve_ambiguity, "---")

    def test_exclude(self):
        """testing excluding gap motif"""
        alpha = self.alpha
        expected = [
            "AAA",
            "AAC",
            "AAG",
            "AAT",
            "ACA",
            "ACC",
            "ACG",
            "ACT",
            "AGA",
            "AGC",
            "AGG",
            "AGT",
            "ATA",
            "ATC",
            "ATG",
            "ATT",
            "CAA",
            "CAC",
            "CAG",
            "CAT",
            "CCA",
            "CCC",
            "CCG",
            "CCT",
            "CGA",
            "CGC",
            "CGG",
            "CGT",
            "CTA",
            "CTC",
            "CTG",
            "CTT",
            "GAA",
            "GAC",
            "GAG",
            "GAT",
            "GCA",
            "GCC",
            "GCG",
            "GCT",
            "GGA",
            "GGC",
            "GGG",
            "GGT",
            "GTA",
            "GTC",
            "GTG",
            "GTT",
            "TAC",
            "TAT",
            "TCA",
            "TCC",
            "TCG",
            "TCT",
            "TGC",
            "TGG",
            "TGT",
            "TTA",
            "TTC",
            "TTG",
            "TTT",
        ]
        self.assertEqualSets(alpha, expected)

    def test_include(self):
        """testing including gap motif"""
        alpha = self.alpha.with_gap_motif()
        expected = [
            "---",
            "AAA",
            "AAC",
            "AAG",
            "AAT",
            "ACA",
            "ACC",
            "ACG",
            "ACT",
            "AGA",
            "AGC",
            "AGG",
            "AGT",
            "ATA",
            "ATC",
            "ATG",
            "ATT",
            "CAA",
            "CAC",
            "CAG",
            "CAT",
            "CCA",
            "CCC",
            "CCG",
            "CCT",
            "CGA",
            "CGC",
            "CGG",
            "CGT",
            "CTA",
            "CTC",
            "CTG",
            "CTT",
            "GAA",
            "GAC",
            "GAG",
            "GAT",
            "GCA",
            "GCC",
            "GCG",
            "GCT",
            "GGA",
            "GGC",
            "GGG",
            "GGT",
            "GTA",
            "GTC",
            "GTG",
            "GTT",
            "TAC",
            "TAT",
            "TCA",
            "TCC",
            "TCG",
            "TCT",
            "TGC",
            "TGG",
            "TGT",
            "TTA",
            "TTC",
            "TTG",
            "TTT",
        ]
        self.assertEqualSets(alpha, expected)

    def test_constructing_from_func(self):
        """the CodonAlphabet function should support genetic code names as well"""
        alpha_int = CodonAlphabet(1)
        alpha_name = CodonAlphabet("Standard Nuclear")
        self.assertEqual(alpha_int, alpha_name)
        alpha_int = CodonAlphabet(2)
        alpha_name = CodonAlphabet("Vertebrate Mitochondrial")
        self.assertEqual(alpha_int, alpha_name)


if __name__ == "__main__":
    main()
