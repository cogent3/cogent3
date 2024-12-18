import pickle
from unittest import TestCase

import pytest

# ind some of the standard alphabets to reduce typing
from numpy.testing import assert_allclose

from cogent3.core import sequence
from cogent3.core.moltype import (
    DNA,
    PROTEIN,
    RNA,
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
        assert make_matches() == {}

    def test_init_monomers(self):
        """make_matches with only monomers should produce {(i,i):True}"""
        m = make_matches("")
        assert m == {}
        m = make_matches("qaz")
        assert m == {("q", "q"): True, ("a", "a"): True, ("z", "z"): True}

    def test_init_gaps(self):
        """make_matches with only gaps should match all gaps to each other"""
        m = make_matches("", "~!")
        assert m == {
            ("~", "~"): True,
            ("!", "!"): True,
            ("!", "~"): True,
            ("~", "!"): True,
        }

    def test_init_degen(self):
        """make_matches with only degen should work as expected"""
        m = make_matches(None, None, {"x": "ab", "y": "bc", "z": "cd", "n": "bcd"})
        assert m == {
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
        }
        assert ("x", "z") not in m

    def test_init_all(self):
        """make_matches with everything should produce correct dict"""
        m = make_matches(
            "ABC",
            ("-", "~"),
            {"X": "AB", "Y": ("B", "C"), "N": list("ABC")},
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
        assert m == exp


class make_pairs_tests(TestCase):
    """Tests of the top-level make_pairs factory function."""

    def setUp(self):
        """Define some standard pairs and other data"""
        self.pairs = {("U", "A"): True, ("A", "U"): True, ("G", "U"): False}

    def test_init_empty(self):
        """make_pairs should init ok with no parameters"""
        assert make_pairs() == {}

    def test_init_pairs(self):
        """make_pairs with just pairs should equal the original"""
        assert make_pairs(self.pairs) == self.pairs
        assert make_pairs(self.pairs) is not self.pairs

    def test_init_monomers(self):
        """make_pairs with pairs and monomers should equal just the pairs"""
        assert make_pairs(self.pairs, "ABCDEFG") == self.pairs
        assert make_pairs(self.pairs, "ABCDEFG") is not self.pairs

    def test_init_gaps(self):
        """make_pairs should add all combinations of gaps as weak pairs"""
        p = make_pairs(self.pairs, None, "-~")
        assert p != self.pairs
        self.pairs.update(
            {
                ("~", "~"): False,
                ("-", "~"): False,
                ("-", "-"): False,
                ("~", "-"): False,
            },
        )
        assert p == self.pairs

    def test_init_degen(self):
        """make_pairs should add in degenerate combinations as weak pairs"""
        p = make_pairs(self.pairs, "AUG", "-", {"R": "AG", "Y": "CU", "W": "AU"})
        assert p != self.pairs
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
            },
        )
        assert p == self.pairs


class CoreObjectGroupTests(TestCase):
    """Tests of the CoreObjectGroup class."""

    def test_init(self):
        """CoreObjectGroup should init with basic list of objects."""

        class o:
            def __init__(self, s):
                self.s = s

        base = o("base")
        c = CoreObjectGroup(base)
        assert c.base is base
        assert c.degen is None
        assert c.base.degen is None

        base, degen, gap, degengap = list(map(o, ["base", "degen", "gap", "degengap"]))
        c = CoreObjectGroup(base, degen, gap, degengap)
        assert c.base is base
        assert c.base.degen is degen
        assert c.degen.gapped is degengap


class AlphabetGroupTests(TestCase):
    """Tests of the AlphabetGroup class."""

    def test_init(self):
        """AlphabetGroup should initialize successfully"""
        chars = "AB"
        degens = {"C": "AB"}
        g = AlphabetGroup(chars, degens)
        assert "".join(g.base) == "AB"
        assert "".join(g.degen) == "ABC"
        assert "".join(g.gapped) == "AB-"
        assert "".join(g.degen_gapped) == "AB-C?"


class MolTypeTests(TestCase):
    """Tests of the MolType class. Should support same API as old Alphabet."""

    def test_to_regex(self):
        """returns a valid regex"""
        seq = "ACYGR"
        regular_expression = DNA.to_regex(seq=seq)
        assert regular_expression == "AC[CT]G[AG]"
        # raises an exception if a string is already a regex, or invalid
        with pytest.raises(ValueError):
            DNA.to_regex("(?:GAT|GAC)(?:GGT|GGC|GGA|GGG)(?:GAT|GAC)(?:CAA|CAG)")

    def test_pickling(self):
        pkl = pickle.dumps(DNA)
        got = pickle.loads(pkl)
        assert isinstance(got, type(DNA))
        assert DNA.ambiguities == got.ambiguities

    def test_get_moltype(self):
        """correctly return a moltype by name"""
        for label in ("dna", "rna", "protein", "protein_with_stop"):
            mt = get_moltype(label)
            assert mt.label == label
            mt = get_moltype(label.upper())
            assert mt.label == label

        mt = get_moltype(DNA)
        assert mt.label == "dna"
        with pytest.raises(ValueError):
            _ = get_moltype("blah")

    def test_available_moltypes(self):
        """available_moltypes returns table with correct columns and rows"""
        available = available_moltypes()
        assert available.shape == (7, 3)
        assert available[1, "Number of states"] == 4
        assert available["dna", "Number of states"] == 4
        txt = repr(available)
        assert "'dna'" in txt

    def test_init_minimal(self):
        """MolType should init OK with just monomers"""
        a = MolType("Abc")
        assert "A" in a.alphabet
        assert "a" not in a.alphabet  # case-sensitive
        assert "b" in a.alphabet
        assert "B" not in a.alphabet
        assert "x" not in a.alphabet

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
            assert i in a
        assert a.complement("b") == "c"
        assert a.complement("AbcAA") == "AcbAA"
        assert a.first_degenerate("AbcdA") == 3
        assert a.first_gap("a~c") == 1
        assert a.first_invalid("Abcx") == 3

    def test_strip_degenerate(self):
        """MolType strip_degenerate should remove any degenerate bases"""
        s = RnaMolType.strip_degenerate
        assert s("UCAG-") == "UCAG-"
        assert s("NRYSW") == ""
        assert s("USNG") == "UG"

    def test_strip_bad(self):
        """MolType strip_bad should remove any non-base, non-gap chars"""
        s = RnaMolType.strip_bad
        assert s("UCAGwsnyrHBND-D") == "UCAGwsnyrHBND-D"
        assert s("@#^*($@!#&()!@QZX") == ""
        assert s("aaa ggg ---!ccc") == "aaaggg---ccc"

    def test_strip_bad_and_gaps(self):
        """MolType strip_bad_and_gaps should remove gaps and bad chars"""
        s = RnaMolType.strip_bad_and_gaps
        assert s("UCAGwsnyrHBND-D") == "UCAGwsnyrHBNDD"
        assert s("@#^*($@!#&()!@QZX") == ""
        assert s("aaa ggg ---!ccc") == "aaagggccc"

    def test_complement(self):
        """MolType complement should correctly complement sequence"""
        assert RnaMolType.complement("UauCG-NR") == "AuaGC-NY"
        assert DnaMolType.complement("TatCG-NR") == "AtaGC-NY"
        assert RnaMolType.complement("") == ""
        self.assertRaises(TypeError, ProteinMolType.complement, "ACD")
        # if it wasn't a string, result should be a list
        assert RnaMolType.complement(list("UauCG-NR")) == list("AuaGC-NY")
        assert RnaMolType.complement(("a", "c")) == ("u", "g")
        # constructor should fail for a dict
        self.assertRaises(ValueError, RnaMolType.complement, {"a": "c"})

    def test_rc(self):
        """MolType rc should correctly reverse-complement sequence"""
        assert RnaMolType.rc("U") == "A"
        assert RnaMolType.rc("") == ""
        assert RnaMolType.rc("R") == "Y"
        assert RnaMolType.rc("UauCG-NR") == "YN-CGauA"
        assert DnaMolType.rc("TatCG-NR") == "YN-CGatA"
        self.assertRaises(TypeError, ProteinMolType.rc, "ACD")
        # if it wasn't a string, result should be a list
        assert RnaMolType.rc(list("UauCG-NR")) == list("YN-CGauA")
        assert RnaMolType.rc(("a", "c")) == ("g", "u")
        # constructor should fail for a dict
        self.assertRaises(ValueError, RnaMolType.rc, {"a": "c"})

    def test_contains(self):
        """MolType contains should return correct result"""
        for i in "UCAGWSMKRYBDHVN-" + "UCAGWSMKRYBDHVN-".lower():
            assert i in RnaMolType
        for i in "x!@#$%^&ZzQq":
            assert i not in RnaMolType

        a = MolType(dict.fromkeys("ABC"), add_lower=True)
        for i in "abcABC":
            assert i in a
        assert "x" not in a
        b = MolType(dict.fromkeys("ABC"), add_lower=False)
        for i in "ABC":
            assert i in b
        for i in "abc":
            assert i not in b

    def test_iter(self):
        """MolType iter should iterate over monomer order"""
        assert list(RnaMolType) == ["U", "C", "A", "G", "u", "c", "a", "g"]
        a = MolType("ZXCV")
        assert list(a) == ["Z", "X", "C", "V"]

    def test_is_gapped(self):
        """MolType is_gapped should return True if gaps in seq"""
        g = RnaMolType.is_gapped
        assert not g("")
        assert not g("ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN")
        assert g("-")
        assert g("--")
        assert g("CAGUCGUACGUCAGUACGUacucauacgac-caguACUG")
        assert g("CA--CGUAUGCA-----g")
        assert g("CAGU-")

    def test_is_gap(self):
        """MolType is_gap should return True if char is a gap"""
        g = RnaMolType.is_gap
        # True for the empty string
        assert not g("")
        # True for all the standard and degenerate symbols
        s = "ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN"
        assert not g(s)
        for i in s:
            assert not g(i)
        # should be true for a single gap
        assert g("-")
        # note that it _shouldn't_ be true for a run of gaps: use a.is_gapped()
        assert not g("--")

    def test_is_degenerate(self):
        """MolType is_degenerate should return True if degen symbol in seq"""
        d = RnaMolType.is_degenerate
        assert not d("")
        assert not d("UACGCUACAUGuacgucaguGCUAGCUA---ACGUCAG")
        assert d("N")
        assert d("R")
        assert d("y")
        assert d("GCAUguagcucgUCAGUCAGUACgUgcasCUAG")
        assert d("ACGYAUGCUGYEWEWNFMNfuwbybcwuybcjwbeiwfub")

    def test_is_valid(self):
        """MolType is_valid should return True if any unknown symbol in seq"""
        v = RnaMolType.is_valid
        assert not v(3)
        assert not v(None)
        assert v("ACGUGCAUGUCAYCAYGUACGcaugacyugc----RYNCYRNC")
        assert v("")
        assert v("a")
        assert not v("ACIUBHFWUIXZKLNJUCIHBICNSOWMOINJ")
        assert not v("CAGUCAGUCACA---GACCAUG-_--cgau")

    def test_is_strict(self):
        """MolType is_strict should return True if all symbols in Monomers"""
        s = RnaMolType.is_strict
        assert not s(3)
        assert not s(None)
        assert s("")
        assert s("A")
        assert s("UAGCACUgcaugcauGCAUGACuacguACAUG")
        assert not s("CAGUCGAUCA-cgaucagUCGAUGAC")
        assert not s("ACGUGCAUXCAGUCAG")

    def test_first_gap(self):
        """MolType first_gap should return index of first gap symbol, or None"""
        g = RnaMolType.first_gap
        assert g("") is None
        assert g("a") is None
        assert g("uhacucHuhacUIUIhacan") is None
        assert g("-abc") == 0
        assert g("b-ac") == 1
        assert g("abcd-") == 4

    def test_first_degenerate(self):
        """MolType first_degenerate should return index of first degen symbol"""
        d = RnaMolType.first_degenerate
        assert d("") is None
        assert d("a") is None
        assert d("UCGACA--CU-gacucaguacgua") is None
        assert d("nCAGU") == 0
        assert d("CUGguagvAUG") == 7
        assert d("ACUGCUAacgud") == 11

    def test_first_invalid(self):
        """MolType first_invalid should return index of first invalid symbol"""
        i = RnaMolType.first_invalid
        assert i("") is None
        assert i("A") is None
        assert i("ACGUNVBuacg-wskmWSMKYRryNnN--") is None
        assert i("x") == 0
        assert i("rx") == 1
        assert i("CAGUNacgunRYWSwx") == 15

    def test_first_non_strict(self):
        """MolType first_non_strict should return index of first non-strict symbol"""
        s = RnaMolType.first_non_strict
        assert s("") is None
        assert s("A") is None
        assert s("ACGUACGUcgaucagu") is None
        assert s("N") == 0
        assert s("-") == 0
        assert s("x") == 0
        assert s("ACGUcgAUGUGCAUcaguX") == 18
        assert s("ACGUcgAUGUGCAUcaguX-38243829") == 18

    def test_disambiguate(self):  # ported
        """MolType disambiguate should remove degenerate bases"""
        d = RnaMolType.disambiguate
        assert d("") == ""
        assert d("AGCUGAUGUA--CAGU") == "AGCUGAUGUA--CAGU"
        assert d("AUn-yrs-wkmCGwmrNMWRKY", "strip") == "AU--CG"
        assert d(tuple("AUn-yrs-wkmCGwmrNMWRKY"), "strip") == tuple("AU--CG")
        s = "AUn-yrs-wkmCGwmrNMWRKY"
        t = d(s, "random")
        u = d(s, "random")
        for i, j in zip(s, t, strict=False):
            if i in RnaMolType.degenerates:
                assert j in RnaMolType.degenerates[i]
            else:
                assert i == j
        assert t != u
        assert d(tuple("UCAG"), "random") == tuple("UCAG")
        assert len(s) == len(t)
        assert RnaMolType.first_degenerate(t) is None
        # should raise exception on unknown disambiguation method
        self.assertRaises(NotImplementedError, d, s, "xyz")

    def test_degap(self):  # ported
        """MolType degap should remove all gaps from sequence"""
        g = RnaMolType.degap
        assert g("") == ""
        assert g("GUCAGUCgcaugcnvuincdks") == "GUCAGUCgcaugcnvuincdks"
        assert g("----------------") == ""
        assert g("gcuauacg-") == "gcuauacg"
        assert g("-CUAGUCA") == "CUAGUCA"
        assert g("---a---c---u----g---") == "acug"
        assert g(tuple("---a---c---u----g---")) == tuple("acug")

    def test_gap_indices(self):
        """MolType gap_indices should return correct gap positions"""
        g = RnaMolType.gap_indices
        assert g("") == []
        assert g("ACUGUCAGUACGHFSDKJCUICDNINS") == []
        assert g("GUACGUIACAKJDC-SDFHJDSFK") == [14]
        assert g("-DSHFUHDSF") == [0]
        assert g("UACHASJAIDS-") == [11]
        assert g("---CGAUgCAU---ACGHc---ACGUCAGU---") == [
            0,
            1,
            2,
            11,
            12,
            13,
            19,
            20,
            21,
            30,
            31,
            32,
        ]
        a = MolType({"A": 1}, gaps=dict.fromkeys("!@#$%"))
        g = a.gap_indices
        assert g("") == []
        assert g("!!!") == [0, 1, 2]
        assert g("!@#$!@#$!@#$") == list(range(12))
        assert g("cguua!cgcuagua@cguasguadc#") == [5, 14, 25]

    def test_gap_vector(self):
        """MolType gap_vector should return correct gap positions"""
        g = RnaMolType.gap_vector
        assert g("") == []
        assert g("ACUGUCAGUACGHFSDKJCUICDNINS") == [False] * 27
        assert g("GUACGUIACAKJDC-SDFHJDSFK") == list(
            map(bool, list(map(int, "000000000000001000000000"))),
        )
        assert g("-DSHFUHDSF") == list(map(bool, list(map(int, "1000000000"))))
        assert g("UACHASJAIDS-") == list(map(bool, list(map(int, "000000000001"))))
        assert g("---CGAUgCAU---ACGHc---ACGUCAGU---") == list(
            map(bool, list(map(int, "111000000001110000011100000000111"))),
        )
        a = MolType({"A": 1}, gaps=dict.fromkeys("!@#$%"))
        g = a.gap_vector
        assert g("") == []
        assert g("!!!") == list(map(bool, [1, 1, 1]))
        assert g("!@#$!@#$!@#$") == [True] * 12
        assert g("cguua!cgcuagua@cguasguadc#") == list(
            map(bool, list(map(int, "00000100000000100000000001"))),
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
        assert gm(empty) == ({}, {})
        assert gm(no_gaps) == ({0: 0, 1: 1, 2: 2}, {0: 0, 1: 1, 2: 2})
        assert gm(all_gaps) == ({}, {})
        assert gm(start_gaps) == ({0: 2, 1: 3, 2: 4}, {2: 0, 3: 1, 4: 2})
        assert gm(end_gaps) == ({0: 0, 1: 1}, {0: 0, 1: 1})
        assert gm(mid_gaps) == ({0: 2, 1: 5, 2: 7, 3: 8}, {2: 0, 5: 1, 7: 2, 8: 3})

    def test_count_gaps(self):
        """MolType count_gaps should return correct gap count"""
        c = RnaMolType.count_gaps
        assert c("") == 0
        assert c("ACUGUCAGUACGHFSDKJCUICDNINS") == 0
        assert c("GUACGUIACAKJDC-SDFHJDSFK") == 1
        assert c("-DSHFUHDSF") == 1
        assert c("UACHASJAIDS-") == 1
        assert c("---CGAUgCAU---ACGHc---ACGUCAGU---") == 12
        a = MolType({"A": 1}, gaps=dict.fromkeys("!@#$%"))
        c = a.count_gaps
        assert c("") == 0
        assert c("!!!") == 3
        assert c("!@#$!@#$!@#$") == 12
        assert c("cguua!cgcuagua@cguasguadc#") == 3

    def test_count_degenerate(self):  # ported
        """MolType count_degenerate should return correct degen base count"""
        d = RnaMolType.count_degenerate
        assert d("") == 0
        assert d("GACUGCAUGCAUCGUACGUCAGUACCGA") == 0
        assert d("N") == 1
        assert d("NRY") == 3
        assert d("ACGUAVCUAGCAUNUCAGUCAGyUACGUCAGS") == 4

    def test_possibilites(self):  # ported
        """MolType possibilities should return correct # possible sequences"""
        p = RnaMolType.possibilities
        assert p("") == 1
        assert p("ACGUgcaucagUCGuGCAU") == 1
        assert p("N") == 4
        assert p("R") == 2
        assert p("H") == 3
        assert p("nRh") == 24
        assert p("AUGCnGUCAg-aurGauc--gauhcgauacgws") == 96

    def test_MW(self):
        """MolType MW should return correct molecular weight"""
        r = RnaMolType.mw
        p = ProteinMolType.mw
        assert p("") == 0
        assert r("") == 0
        assert_allclose(p("A"), 89.09)
        assert_allclose(r("A"), 375.17)
        assert_allclose(p("AAA"), 231.27)
        assert_allclose(r("AAA"), 1001.59)
        assert_allclose(r("AAACCCA"), 2182.37)

    def test_can_match(self):
        """MolType can_match should return True if all positions can match"""
        m = RnaMolType.can_match
        assert m("", "")
        assert m("UCAG", "UCAG")
        assert not m("UCAG", "ucag")
        assert m("UCAG", "NNNN")
        assert m("NNNN", "UCAG")
        assert m("NNNN", "NNNN")
        assert not m("N", "x")
        assert not m("N", "-")
        assert m("UCAG", "YYRR")
        assert m("UCAG", "KMWS")

    def test_can_mismatch(self):
        """MolType can_mismatch should return True on any possible mismatch"""
        m = RnaMolType.can_mismatch
        assert not m("", "")
        assert m("N", "N")
        assert m("R", "R")
        assert m("N", "r")
        assert m("CGUACGCAN", "CGUACGCAN")
        assert m("U", "C")
        assert m("UUU", "UUC")
        assert m("UUU", "UUY")
        assert not m("UUU", "UUU")
        assert not m("UCAG", "UCAG")
        assert not m("U--", "U--")

    def test_must_match(self):
        """MolType must_match should return True when no possible mismatches"""
        m = RnaMolType.must_match
        assert m("", "")
        assert not m("N", "N")
        assert not m("R", "R")
        assert not m("N", "r")
        assert not m("CGUACGCAN", "CGUACGCAN")
        assert not m("U", "C")
        assert not m("UUU", "UUC")
        assert not m("UUU", "UUY")
        assert m("UU-", "UU-")
        assert m("UCAG", "UCAG")

    def test_can_pair(self):
        """MolType can_pair should return True if all positions can pair"""
        p = RnaMolType.can_pair
        assert p("", "")
        assert not p("UCAG", "UCAG")
        assert p("UCAG", "CUGA")
        assert not p("UCAG", "cuga")
        assert p("UCAG", "NNNN")
        assert p("NNNN", "UCAG")
        assert p("NNNN", "NNNN")
        assert not p("N", "x")
        assert not p("N", "-")
        assert p("-", "-")
        assert p("UCAGU", "KYYRR")
        assert p("UCAG", "KKRS")
        assert p("U", "G")

        d = DnaMolType.can_pair
        assert not d("T", "G")

    def test_can_mispair(self):
        """MolType can_mispair should return True on any possible mispair"""
        m = RnaMolType.can_mispair
        assert not m("", "")
        assert m("N", "N")
        assert m("R", "Y")
        assert m("N", "r")
        assert m("CGUACGCAN", "NUHCHUACH")
        assert m("U", "C")
        assert m("U", "R")
        assert m("UUU", "AAR")
        assert m("UUU", "GAG")
        assert not m("UUU", "AAA")
        assert not m("UCAG", "CUGA")
        assert m("U--", "--U")

        d = DnaMolType.can_pair
        assert d("TCCAAAGRYY", "RRYCTTTGGA")

    def test_must_pair(self):
        """MolType must_pair should return True when no possible mispairs"""
        m = RnaMolType.must_pair
        assert m("", "")
        assert not m("N", "N")
        assert not m("R", "Y")
        assert not m("A", "A")
        assert not m("CGUACGCAN", "NUGCGUACG")
        assert not m("U", "C")
        assert not m("UUU", "AAR")
        assert not m("UUU", "RAA")
        assert not m("UU-", "-AA")
        assert m("UCAG", "CUGA")

        d = DnaMolType.must_pair
        assert d("TCCAGGG", "CCCTGGA")
        assert d("tccaggg", "ccctgga")
        assert not d("TCCAGGG", "NCCTGGA")


class RnaMolTypeTests(TestCase):
    """Spot-checks of alphabet functionality applied to RNA alphabet."""

    def test_contains(self):
        """RnaMolType should __contain__ the expected symbols."""
        keys = "ucagrymkwsbhvdn?-"
        for k in keys:
            assert k in RnaMolType
        for k in keys.upper():
            assert k in RnaMolType
        assert "X" not in RnaMolType

    def test_degenerate_from_seq(self):
        """RnaMolType degenerate_from_seq should give correct results"""
        d = RnaMolType.degenerate_from_seq
        # check monomers
        assert d("a") == "a"
        assert d("C") == "C"
        # check seq of monomers
        assert d("aaaaa") == "a"
        # check some 2- to 4-way cases
        assert d("aaaaag") == "r"
        assert d("ccgcgcgcggcc") == "s"
        assert d("accgcgcgcggcc") == "v"
        assert d("aaaaagcuuu") == "n"
        # check some cases with gaps
        assert d("aa---aaagcuuu") == "?"
        assert d("aaaaaaaaaaaaaaa-") == "?"
        assert d("----------------") == "-"
        # check mixed case example
        assert d("AaAAaa") == "A"
        # check example with degenerate symbols in set
        assert d("RS") == "V"
        assert d("RN-") == "?"
        # check that it works for proteins as well
        p = ProteinMolType.degenerate_from_seq
        assert p("A") == "A"
        assert p("AAA") == "A"
        assert p("DN") == "B"
        assert p("---") == "-"
        assert p("ACD") == "X"
        assert p("ABD") == "X"
        assert p("ACX") == "X"
        assert p("AC-") == "?"

    def test_get_css_styles(self):
        """generates css style mapping and css"""
        css, styles = DNA.get_css_style()
        # styles should consist of core characters appended by DNA.label
        states = [*list(DNA), "ambig", "terminal_ambig"]
        got = {styles[b] for b in states}
        expect = {b + "_dna" for b in states}
        assert got == expect
        for state in expect:
            assert state in "\n".join(css)

        # check subset of protein
        css, styles = PROTEIN.get_css_style()
        states = [*list(PROTEIN), "ambig", "terminal_ambig"]
        got = {styles[b] for b in states}
        expect = {b + "_protein" for b in states}
        assert got == expect

    ("new MolType will not support setting label")

    def test_get_css_no_label(self):
        """should not fail when moltype has no label"""
        dna = get_moltype("dna")
        orig_label, dna.label = dna.label, None
        _ = dna.get_css_style()
        dna.label = orig_label


class _AlphabetTestCase(TestCase):
    def assertEqualSeqs(self, a, b):
        """For when we don't care about the type, just the elements"""
        assert list(a) == list(b)

    def assertEqualSets(self, a, b):
        assert set(a) == set(b)


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
            AlphabetError,
            alpha.get_subset,
            motif_subset=["AA", "CA", "GT", "TT"],
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

    def test_strand_symmetric_motifs(self):  # ported
        """construction of strand symmetric motif sets"""
        # fails for a moltype with no strand complement
        with pytest.raises(TypeError):
            PROTEIN.strand_symmetric_motifs()

        got = DNA.strand_symmetric_motifs(motif_length=1)
        expect = {("A", "T"), ("C", "G")}
        assert got == expect
        got = RNA.strand_symmetric_motifs(motif_length=1)
        expect = {("A", "U"), ("C", "G")}
        assert got == expect
        got = DNA.strand_symmetric_motifs(motif_length=2)
        assert len(got) == 8
        got = DNA.strand_symmetric_motifs(motif_length=3)
        assert len(got) == 32


class TestCodonAlphabet(_AlphabetTestCase):
    def setUp(self):
        from cogent3.core.moltype import STANDARD_CODON

        self.alpha = STANDARD_CODON

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
        alpha_name = CodonAlphabet("Standard")
        assert alpha_int == alpha_name
        alpha_int = CodonAlphabet(2)
        alpha_name = CodonAlphabet("Vertebrate Mitochondrial")
        assert alpha_int == alpha_name


def test_resolve_ambiguity_nucs():  # ported
    got = DNA.resolve_ambiguity("AT?", allow_gap=False)
    assert len(got) == 4
    assert len(got[0]) == 3


def test_resolve_ambiguity_codons():  # ported
    from cogent3 import get_code

    gc = get_code(1)
    codon_alpha = gc.get_alphabet(include_stop=False)
    codon_alpha_w_gap = codon_alpha.with_gap_motif()
    assert len(DNA.resolve_ambiguity("AT?", alphabet=codon_alpha_w_gap)) == 4
    assert len(DNA.resolve_ambiguity("???", alphabet=codon_alpha_w_gap)) == 62
    assert len(DNA.resolve_ambiguity("---", alphabet=codon_alpha_w_gap)) == 1

    assert len(DNA.resolve_ambiguity("AT?", alphabet=codon_alpha)) == 4
    assert len(DNA.resolve_ambiguity("???", alphabet=codon_alpha)) == 61

    with pytest.raises(AlphabetError):
        DNA.resolve_ambiguity("at-")
    with pytest.raises(AlphabetError):
        DNA.resolve_ambiguity("---", alphabet=codon_alpha)


def test_make_seq_on_seq():
    seq = DNA.make_seq(seq="ACGG")
    got = DNA.make_seq(seq=seq)
    assert got is seq


("Dropping support for coerce_str")


def test_make_seq_diff_moltype():
    seq = RNA.make_seq(seq="ACGG")
    seq.add_feature(biotype="gene", name="test", spans=[(0, 2)])
    got = DNA.make_seq(seq=seq)
    assert len(got.annotation_db) == 1


def test_is_compatible_alphabet():
    from cogent3.core.alphabet import CharAlphabet

    dna = get_moltype("dna")
    alpha = CharAlphabet("TCAG")
    assert dna.is_compatible_alphabet(alpha)
    rna = get_moltype("rna")
    assert not rna.is_compatible_alphabet(alpha)
    alpha = CharAlphabet(list(dna.ambiguities))
    prot = get_moltype("protein")
    assert not prot.is_compatible_alphabet(alpha)


def test_is_compatible_alphabet_strict():
    from cogent3.core.alphabet import CharAlphabet

    dna = get_moltype("dna")
    alpha1 = CharAlphabet("TCAG")
    assert dna.is_compatible_alphabet(alpha1, strict=True)
    # returns False if the order is not exactly the same
    alpha1 = CharAlphabet("CTAG")
    assert not dna.is_compatible_alphabet(alpha1, strict=True)
