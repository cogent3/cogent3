import json
import os
import re
import typing
from pickle import dumps

import numpy
import pytest

import cogent3
from cogent3._version import __version__
from cogent3.core import new_alphabet, new_genetic_code, new_moltype, new_sequence
from cogent3.util.deserialise import deserialise_object
from cogent3.util.misc import get_object_provenance


@pytest.fixture(scope="function")
def dna_alphabet():
    return new_moltype.DNA.degen_gapped_alphabet


@pytest.fixture(scope="function")
def ascii_alphabet():
    return new_moltype.ASCII.alphabet


@pytest.fixture
def bytes_alphabet():
    return new_moltype.BYTES.most_degen_alphabet()


@pytest.fixture(scope="function")
def integer_seq(bytes_alphabet):
    """Used for slicing tests"""
    return new_sequence.SeqView(parent="0123456789", alphabet=bytes_alphabet)


@pytest.mark.parametrize("name", ("dna", "rna", "protein", "protein_with_stop", "text"))
def test_moltype_make_seq(name):
    raw = "ACGGA"
    moltype = new_moltype.get_moltype(name)
    seq = moltype.make_seq(name="s1", seq=raw)
    assert seq.moltype.name == name
    assert str(seq) == raw


def test_moltype_make_bytes_seq():
    raw = "ACGGA"
    name = "bytes"
    moltype = new_moltype.get_moltype(name)
    seq = moltype.make_seq(name="s1", seq=raw)
    assert seq.moltype.name == name
    assert str(seq) == raw


# Tests from test_sequence.py


@pytest.mark.parametrize(
    "moltype", ("dna", "rna", "protein", "protein_with_stop", "text")
)
def test_sequence_copy(moltype):
    """correctly returns a copy version of self"""
    mt = new_moltype.get_moltype(moltype)
    s = mt.make_seq(seq="CCCCCCCCCCCCCAAAA", name="test_copy")
    annot1 = s.add_feature(biotype="exon", name="annot1", spans=[(0, 10)])
    annot2 = s.add_feature(biotype="exon", name="annot2", spans=[(10, 14)])
    got = s.copy()
    got_annot1 = list(got.get_features(biotype="exon", name="annot1"))[0]
    got_annot2 = list(got.get_features(biotype="exon", name="annot2"))[0]
    assert got is not s
    assert got_annot1 is not annot1
    assert got_annot2 is not annot2
    assert got.name == s.name
    assert got.info == s.info
    assert str(got) == str(s)
    assert got.moltype == s.moltype
    annot1_slice = str(annot1.get_slice())
    annot2_slice = str(annot2.get_slice())
    got1_slice = str(got_annot1.get_slice())
    got2_slice = str(got_annot2.get_slice())
    assert annot1_slice == got1_slice
    assert annot2_slice == got2_slice


@pytest.mark.parametrize(
    "moltype", ("dna", "rna", "protein", "protein_with_stop", "text")
)
@pytest.mark.parametrize("seq", ("ACG", "AC-G", "-A-C"))
def test_sequence_compare_to_string(moltype, seq):
    """Sequence should compare equal to same string."""
    mt = new_moltype.get_moltype(moltype)
    s = mt.make_seq(seq=seq)
    assert s == seq


def test_sequence_slice():
    """Sequence slicing should work as expected"""
    r = new_moltype.RNA.make_seq(seq="UCAGG")
    assert r[0] == "U"
    assert r[-1] == "G"
    assert r[1:3] == "CA"


def test_sequence_to_dna():
    """Returns copy of self as DNA."""
    r = new_moltype.RNA.make_seq(seq="UCA")
    assert str(r) == "UCA"
    assert str(r.to_dna()) == "TCA"


def test_sequence_to_rna():
    """Returns copy of self as RNA."""
    r = new_moltype.DNA.make_seq(seq="TCA")
    assert str(r) == "TCA"
    assert str(r.to_rna()) == "UCA"


def test_sequence_to_fasta():
    """Sequence.to_fasta() should return Fasta-formatted string"""
    even = "TCAGAT"
    odd = f"{even}AAA"
    even_dna = new_moltype.DNA.make_seq(seq=even, name="even")
    odd_dna = new_moltype.DNA.make_seq(seq=odd, name="odd")
    assert even_dna.to_fasta() == ">even\nTCAGAT\n"
    # set line wrap to small number so we can test that it works
    assert even_dna.to_fasta(block_size=2) == ">even\nTC\nAG\nAT\n"
    assert odd_dna.to_fasta(block_size=2) == ">odd\nTC\nAG\nAT\nAA\nA\n"
    # check that changing the linewrap again works
    assert even_dna.to_fasta(block_size=4) == ">even\nTCAG\nAT\n"


def test_sequence_serialize():
    """Sequence should be serializable"""
    r = new_moltype.RNA.make_seq(seq="UGAGG")
    assert dumps(r)


def test_sequence_to_moltype():
    """correctly convert to specified moltype"""
    s = new_moltype.ASCII.make_seq(seq="TTTTTTTTTTAAAA", name="test1")
    s.add_feature(biotype="exon", name="fred", spans=[(0, 10)])
    s.add_feature(biotype="exon", name="trev", spans=[(10, 14)])
    got = s.to_moltype("dna")
    fred = list(got.get_features(name="fred"))[0]
    assert str(got[fred]) == "TTTTTTTTTT"
    trev = list(got.get_features(name="trev"))[0]
    assert str(got[trev]) == "AAAA"

    # should raise exception if moltype not compatible with sequence data
    with pytest.raises(ValueError):
        s.to_moltype("rna")

    # calling with a null object should raise an exception
    with pytest.raises(ValueError):
        s.to_moltype(None)

    with pytest.raises(ValueError):
        s.to_moltype("")


@pytest.mark.parametrize(
    "seq, expect", (("UCAG-", "UCAG-"), ("NRYSW", ""), ("USNG", "UG"))
)
def test_sequence_strip_degenerate(seq, expect):
    """Sequence strip_degenerate should remove any degenerate bases"""
    seq = new_moltype.RNA.make_seq(seq=seq)
    got = seq.strip_degenerate()
    assert got == expect


def test_add():
    mt = new_moltype.DNA
    seq1 = mt.make_seq(seq="AAA")
    seq2 = mt.make_seq(seq="CCC")
    got = seq1 + seq2
    assert got == "AAACCC"
    assert got.moltype == mt


def test_add_bad():
    mt = new_moltype.DNA
    seq1 = mt.make_seq(seq="AAA")
    with pytest.raises(new_alphabet.AlphabetError):
        _ = seq1 + "s8d3j%31 s-']"


@pytest.mark.parametrize(
    "moltype,label",
    [
        (new_moltype.ASCII, "text"),
        (new_moltype.BYTES, "bytes"),
        (new_moltype.DNA, "dna"),
    ],
)
def test_get_type(moltype, label):
    seq = moltype.make_seq(seq="ARCGT")
    got = seq.get_type()
    assert got == label


@pytest.mark.parametrize(
    "s, expect",
    [
        ("ARC", [set(("A")), set(("A", "G")), set(("C",))]),
        ("AGC", [set(("A",)), set(("G",)), set(("C",))]),
    ],
)
def test_resolved_ambiguities(s, expect):
    mt = new_moltype.DNA
    seq = mt.make_seq(seq=s)
    got = seq.resolved_ambiguities()
    assert got == expect


@pytest.mark.parametrize(
    "seq, expect",
    (
        ("UCXXXAGWSNYRHBNZZZD-D", "UCAGWSNYRHBND-D"),
        ("@#^*($@!#&()!@QZX", ""),
        ("aaaxggg---!ccc", "---"),
    ),
)
def test_sequence_strip_bad(seq, expect):
    """Sequence strip_bad should remove any non-base, non-gap chars"""
    # have to turn off check to get bad data in
    seq = new_moltype.RNA.make_seq(seq=seq, check_seq=False)
    got = seq.strip_bad()
    assert str(got) == expect


@pytest.mark.parametrize(
    "seq, expect",
    (
        ("UXXCAGWSNYRHBNZ#!D-D", "UCAGWSNYRHBNDD"),
        ("@#^*($@!#&()!@QZX", ""),
        ("AAA GGG ---!CCC", "AAAGGGCCC"),
    ),
)
def test_sequence_strip_bad_and_gaps(seq, expect):
    """Sequence strip_bad_and_gaps should remove gaps and bad chars"""
    # have to turn off check to get bad data in; no longer preserves case
    seq = new_moltype.RNA.make_seq(seq=seq, check_seq=False)
    got = seq.strip_bad_and_gaps()
    assert str(got) == expect


def test_sequence_shuffle():
    """Sequence shuffle should return new random sequence w/ same monomers"""
    r = new_moltype.RNA.make_seq(seq="UUUUCCCCAAAAGGGG")
    s = r.shuffle()
    assert r != s
    # assert the number of counts of each monomer is the same
    assert r.counts() == s.counts()


def test_sequence_complement():
    """Sequence complement should correctly complement sequence"""
    got = new_moltype.RNA.make_seq(seq="UAUCG-NR").complement()
    assert got == "AUAGC-NY"
    got = new_moltype.DNA.make_seq(seq="TATCG-NR").complement()
    assert got == "ATAGC-NY"
    got = new_moltype.DNA.make_seq(seq="").complement()
    assert got == ""
    with pytest.raises(AttributeError):
        new_moltype.PROTEIN.make_seq(seq="ACD").complement()


def test_sequence_rc():
    """Sequence.rc() should correctly reverse-complement sequence"""
    # no longer preserves case!
    assert new_moltype.DNA.make_seq(seq="TATCG-NR").rc() == "YN-CGATA"
    assert new_moltype.RNA.make_seq(seq="").rc() == ""
    assert new_moltype.RNA.make_seq(seq="UAUCG-NR").rc() == "YN-CGAUA"
    assert new_moltype.RNA.make_seq(seq="A").rc() == "U"
    with pytest.raises(AttributeError):
        new_moltype.PROTEIN.make_seq(seq="ACD").rc()


def test_sequence_contains():
    """Sequence contains should return correct result"""
    r = new_moltype.RNA.make_seq(seq="UCA")
    assert "U" in r
    assert "CA" in r
    assert "X" not in r
    assert "G" not in r


def test_sequence_iter():
    """Sequence iter should iterate over sequence"""
    p = new_moltype.PROTEIN.make_seq(seq="QWE")
    assert list(p) == ["Q", "W", "E"]


def test_sequence_is_gapped():
    """Sequence is_gapped should return True if gaps in seq"""
    assert not new_moltype.RNA.make_seq(seq="").is_gapped()
    assert not new_moltype.RNA.make_seq(seq="ACGUCAGUACGUCAGNRCGAUYRNRYRN").is_gapped()
    assert new_moltype.RNA.make_seq(seq="-").is_gapped()
    assert new_moltype.PROTEIN.make_seq(seq="--").is_gapped()
    assert new_moltype.RNA.make_seq(seq="CAGUCGUACGUCAGUACGU-ACUG").is_gapped()
    assert new_moltype.RNA.make_seq(seq="CA--CGUAUGCA-----G").is_gapped()
    assert new_moltype.RNA.make_seq(seq="CAGU-").is_gapped()


def test_sequence_is_degenerate():
    """Sequence is_degenerate should return True if degen symbol in seq"""
    assert not new_moltype.RNA.make_seq(seq="").is_degenerate()
    assert not new_moltype.RNA.make_seq(
        seq="UACGCUACAUGGCUAGCUA---ACGUCAG"
    ).is_degenerate()
    assert new_moltype.RNA.make_seq(seq="N").is_degenerate()
    assert new_moltype.RNA.make_seq(seq="R").is_degenerate()
    assert new_moltype.RNA.make_seq(seq="Y").is_degenerate()
    assert new_moltype.RNA.make_seq(seq="GCSUA").is_degenerate()
    assert new_moltype.RNA.make_seq(seq="ACGYAUGCUGYWWNMN").is_degenerate()


def test_sequence_is_strict():
    """Sequence is_strict should return True if all symbols in Monomers"""
    assert new_moltype.RNA.make_seq(seq="").is_strict()
    assert new_moltype.PROTEIN.make_seq(seq="A").is_strict()
    assert new_moltype.RNA.make_seq(seq="UAGCACU").is_strict()
    assert not new_moltype.RNA.make_seq(seq="CAGUCGAUCA-").is_strict()


def test_sequence_disambiguate():
    """Sequence disambiguate should remove degenerate bases

    Notes
    -----
    This test relies on random generation not being the same twice!
    """
    assert new_moltype.RNA.make_seq(seq="").disambiguate() == ""
    assert (
        new_moltype.RNA.make_seq(seq="AGCUGAUGUA--CAGU").disambiguate()
        == "AGCUGAUGUA--CAGU"
    )
    assert new_moltype.RNA.make_seq(seq="AU-YRS-CG").disambiguate("strip") == "AU--CG"
    s = new_moltype.RNA.make_seq(seq="AUN-YRS-WKMCGWMRNMWRKY")
    t = s.disambiguate("random")
    u = s.disambiguate("random")
    for i, j in zip(str(s), str(t)):
        if i in s.moltype.ambiguities:
            assert j in s.moltype.ambiguities[i]
        else:
            assert i == j
    assert t != u
    assert len(s) == len(t)


def test_sequence_degap():
    """Sequence degap should remove all gaps from sequence"""
    # doesn't preserve case
    assert new_moltype.RNA.make_seq(seq="").degap() == ""
    assert new_moltype.RNA.make_seq(seq="GUCAGUC").degap() == "GUCAGUC"
    assert new_moltype.RNA.make_seq(seq="----------------").degap() == ""
    assert new_moltype.RNA.make_seq(seq="GCUAUACG-").degap() == "GCUAUACG"
    assert new_moltype.RNA.make_seq(seq="-CUAGUCA").degap() == "CUAGUCA"
    assert new_moltype.RNA.make_seq(seq="---A---C---U----G---").degap() == "ACUG"
    assert new_moltype.RNA.make_seq(seq="?A-").degap() == "A"


def test_sequence_gap_indices():
    """Sequence gap_indices should return correct gap positions"""
    got = new_moltype.RNA.make_seq(seq="").gap_indices()
    expect = numpy.array([])
    assert numpy.array_equal(got, expect)

    got = new_moltype.RNA.make_seq(seq="ACUGUCAGUACGHSDKCUCDNNS").gap_indices()
    expect = numpy.array([])
    assert numpy.array_equal(got, expect)

    got = new_moltype.RNA.make_seq(seq="GUACGUACAKDC-SDHDSK").gap_indices()
    expect = numpy.array([12])
    assert numpy.array_equal(got, expect)

    got = new_moltype.RNA.make_seq(seq="-DSHUHDS").gap_indices()
    expect = numpy.array([0])
    assert numpy.array_equal(got, expect)

    got = new_moltype.RNA.make_seq(seq="UACHASADS-").gap_indices()
    expect = numpy.array([9])
    assert numpy.array_equal(got, expect)

    got = new_moltype.RNA.make_seq(
        seq="---CGAUGCAU---ACGHC---ACGUCAGU---"
    ).gap_indices()
    expect = numpy.array([0, 1, 2, 11, 12, 13, 19, 20, 21, 30, 31, 32])
    assert numpy.array_equal(got, expect)


def test_sequence_gap_vector():
    """Sequence gap_vector should return correct gap positions"""

    def g(x):
        return new_moltype.RNA.make_seq(seq=x).gap_vector()

    assert g("") == []
    assert g("ACUGUCAGUACGHCSDKCCUCCDNCNS") == [False] * 27
    assert g("GUACGUAACAKADC-SDAHADSAK") == list(
        map(bool, list(map(int, "000000000000001000000000")))
    )
    assert g("-DSHSUHDSS") == list(map(bool, list(map(int, "1000000000"))))
    assert g("UACHASCAGDS-") == list(map(bool, list(map(int, "000000000001"))))
    assert g("---CGAUCAU---ACGH---ACGUCAGU--?") == list(
        map(bool, list(map(int, "1110000000111000011100000000111")))
    )


def test_count_gaps():
    """Sequence.count_gaps should return correct gap count"""
    assert new_moltype.RNA.make_seq(seq="").count_gaps() == 0
    assert new_moltype.RNA.make_seq(seq="ACUGUCAGUACGHSDKCUCDNNS").count_gaps() == 0
    assert new_moltype.RNA.make_seq(seq="GUACGUACAKDC-SDHDSK").count_gaps() == 1
    assert new_moltype.RNA.make_seq(seq="-DSHUHDS").count_gaps() == 1
    assert new_moltype.RNA.make_seq(seq="UACHASADS-").count_gaps() == 1
    assert (
        new_moltype.RNA.make_seq(seq="---CGAUGCAU---ACGHC---ACGUCAGU---").count_gaps()
        == 12
    )


def test_count_degenerate():
    """Sequence.count_degenerate should return correct degen base count"""
    assert new_moltype.RNA.make_seq(seq="").count_degenerate() == 0
    assert (
        new_moltype.RNA.make_seq(seq="GACUGCAUGCAUCGUACGUCAGUACCGA").count_degenerate()
        == 0
    )
    assert new_moltype.RNA.make_seq(seq="N").count_degenerate() == 1
    assert new_moltype.PROTEIN.make_seq(seq="N").count_degenerate() == 0
    assert new_moltype.RNA.make_seq(seq="NRY").count_degenerate() == 3
    assert (
        new_moltype.RNA.make_seq(
            seq="ACGUAVCUAGCAUNUCAGUCAGYUACGUCAGS"
        ).count_degenerate()
        == 4
    )


def test_sequence_count_variants():
    """Sequence count_variants should return correct # possible sequences"""
    assert new_moltype.RNA.make_seq(seq="").count_variants() == 1
    assert new_moltype.RNA.make_seq(seq="ACGUGCAU").count_variants() == 1
    assert new_moltype.RNA.make_seq(seq="N").count_variants() == 4
    assert new_moltype.RNA.make_seq(seq="R").count_variants() == 2
    assert new_moltype.RNA.make_seq(seq="H").count_variants() == 3
    assert new_moltype.RNA.make_seq(seq="NRH").count_variants() == 24
    assert (
        new_moltype.RNA.make_seq(
            seq="AUGCNGUCAG-AURGAUC--GAUHCGAUACGWS"
        ).count_variants()
        == 96
    )


def test_mw():
    """Sequence MW should return correct molecular weight"""
    assert new_moltype.PROTEIN.make_seq(seq="").mw() == 0
    assert new_moltype.RNA.make_seq(seq="").mw() == 0
    assert numpy.allclose(new_moltype.PROTEIN.make_seq(seq="A").mw(), 89.09)
    assert numpy.allclose(new_moltype.RNA.make_seq(seq="A").mw(), 375.17)
    assert numpy.allclose(new_moltype.PROTEIN.make_seq(seq="AAA").mw(), 231.27)
    assert numpy.allclose(new_moltype.RNA.make_seq(seq="AAA").mw(), 1001.59)
    assert numpy.allclose(new_moltype.RNA.make_seq(seq="AAACCCA").mw(), 2182.37)


def test_can_match():
    """Sequence can_match should return True if all positions can match"""
    assert new_moltype.RNA.make_seq(seq="").can_match("")
    assert new_moltype.RNA.make_seq(seq="UCAG").can_match("UCAG")
    assert new_moltype.RNA.make_seq(seq="UCAG").can_match("NNNN")
    assert new_moltype.RNA.make_seq(seq="NNNN").can_match("UCAG")
    assert new_moltype.RNA.make_seq(seq="NNNN").can_match("NNNN")
    assert not new_moltype.RNA.make_seq(seq="N").can_match("X")
    assert not new_moltype.RNA.make_seq(seq="N").can_match("-")
    assert new_moltype.RNA.make_seq(seq="UCAG").can_match("YYRR")
    assert new_moltype.RNA.make_seq(seq="UCAG").can_match("KMWS")


def test_can_pair():
    """Sequence can_pair should return True if all positions can pair"""
    assert new_moltype.RNA.make_seq(seq="").can_pair("")
    assert not new_moltype.RNA.make_seq(seq="UCAG").can_pair("UCAG")
    assert new_moltype.RNA.make_seq(seq="UCAG").can_pair("CUGA")
    assert not new_moltype.RNA.make_seq(seq="UCAG").can_pair("cuga")
    assert new_moltype.RNA.make_seq(seq="UCAG").can_pair("NNNN")
    assert new_moltype.RNA.make_seq(seq="NNNN").can_pair("UCAG")
    assert new_moltype.RNA.make_seq(seq="NNNN").can_pair("NNNN")
    assert not new_moltype.RNA.make_seq(seq="N").can_pair("x")
    assert not new_moltype.RNA.make_seq(seq="N").can_pair("-")
    assert new_moltype.RNA.make_seq(seq="-").can_pair("-")
    assert new_moltype.RNA.make_seq(seq="UCAGU").can_pair("KYYRR")
    assert new_moltype.RNA.make_seq(seq="UCAG").can_pair("KKRS")
    assert new_moltype.RNA.make_seq(seq="U").can_pair("G")

    assert not new_moltype.DNA.make_seq(seq="T").can_pair("G")


def test_can_mispair():
    """Sequence can_mispair should return True on any possible mispair"""
    assert not new_moltype.RNA.make_seq(seq="").can_mispair("")
    assert new_moltype.RNA.make_seq(seq="N").can_mispair("N")
    assert new_moltype.RNA.make_seq(seq="R").can_mispair("Y")
    assert new_moltype.RNA.make_seq(seq="N").can_mispair("r")
    assert new_moltype.RNA.make_seq(seq="CGUACGCAN").can_mispair("NUHCHUACH")
    assert new_moltype.RNA.make_seq(seq="U").can_mispair("C")
    assert new_moltype.RNA.make_seq(seq="U").can_mispair("R")
    assert new_moltype.RNA.make_seq(seq="UUU").can_mispair("AAR")
    assert new_moltype.RNA.make_seq(seq="UUU").can_mispair("GAG")
    assert not new_moltype.RNA.make_seq(seq="UUU").can_mispair("AAA")
    assert not new_moltype.RNA.make_seq(seq="UCAG").can_mispair("CUGA")
    assert new_moltype.RNA.make_seq(seq="U--").can_mispair("--U")
    assert new_moltype.DNA.make_seq(seq="TCCAAAGRYY").can_mispair("RRYCTTTGGA")


def test_must_pair():
    """Sequence must_pair should return True when no possible mispairs"""
    assert new_moltype.RNA.make_seq(seq="").must_pair("")
    assert not new_moltype.RNA.make_seq(seq="N").must_pair("N")
    assert not new_moltype.RNA.make_seq(seq="R").must_pair("Y")
    assert not new_moltype.RNA.make_seq(seq="A").must_pair("A")
    assert not new_moltype.RNA.make_seq(seq="CGUACGCAN").must_pair("NUGCGUACG")
    assert not new_moltype.RNA.make_seq(seq="U").must_pair("C")
    assert not new_moltype.RNA.make_seq(seq="UUU").must_pair("AAR")
    assert not new_moltype.RNA.make_seq(seq="UUU").must_pair("RAA")
    assert not new_moltype.RNA.make_seq(seq="UU-").must_pair("-AA")
    assert new_moltype.RNA.make_seq(seq="UCAG").must_pair("CUGA")

    assert new_moltype.DNA.make_seq(seq="TCCAGGG").must_pair("CCCTGGA")
    assert new_moltype.DNA.make_seq(seq="TCCAGGG").must_pair(
        new_moltype.DNA.make_seq(seq="CCCTGGA")
    )
    assert not new_moltype.DNA.make_seq(seq="TCCAGGG").must_pair("NCCTGGA")


def test_diff():
    """Sequence diff should count 1 for each difference between sequences"""
    assert new_moltype.RNA.make_seq(seq="UGCUGCUC").diff("") == 0
    assert new_moltype.RNA.make_seq(seq="UGCUGCUC").diff("U") == 0
    assert new_moltype.RNA.make_seq(seq="UGCUGCUC").diff("UCCCCCUC") == 3
    assert new_moltype.RNA.make_seq(seq="AAAAA").diff("CCCCC") == 5
    # raises TypeError if other not iterable
    with pytest.raises(TypeError):
        new_moltype.RNA.make_seq(seq="AAAAA").diff(5)


def test_distance():
    """Sequence distance should calculate correctly based on function"""

    def f(a, b):
        if a == b:
            return 0
        if (a in "UC" and b in "UC") or (a in "AG" and b in "AG"):
            return 1
        else:
            return 10

    # uses identity function by default
    assert new_moltype.RNA.make_seq(seq="UGCUGCUC").distance("") == 0
    assert new_moltype.RNA.make_seq(seq="UGCUGCUC").distance("U") == 0
    assert new_moltype.RNA.make_seq(seq="UGCUGCUC").distance("UCCCCCUC") == 3
    assert new_moltype.RNA.make_seq(seq="AAAAA").distance("CCCCC") == 5
    # should use function if supplied
    assert new_moltype.RNA.make_seq(seq="UGCUGCUC").distance("", f) == 0
    assert new_moltype.RNA.make_seq(seq="UGCUGCUC").distance("U", f) == 0
    assert new_moltype.RNA.make_seq(seq="UGCUGCUC").distance("C", f) == 1
    assert new_moltype.RNA.make_seq(seq="UGCUGCUC").distance("G", f) == 10
    assert new_moltype.RNA.make_seq(seq="UGCUGCUC").distance("UCCCCCUC", f) == 21
    assert new_moltype.RNA.make_seq(seq="AAAAA").distance("CCCCC", f) == 50


def test_matrix_distance():
    """Sequence matrix_distance should look up distances from a matrix"""
    # note that the score matrix must contain 'diagonal' elements m[i][i]
    # to avoid failure when the sequences match.
    m = {"U": {"U": 0, "C": 1, "A": 5}, "C": {"C": 0, "A": 2, "G": 4}}
    assert new_moltype.RNA.make_seq(seq="UUUCCC").matrix_distance("UCACGG", m) == 14
    assert new_moltype.RNA.make_seq(seq="UUUCCC").matrix_distance("", m) == 0
    assert new_moltype.RNA.make_seq(seq="UUU").matrix_distance("CAC", m) == 7
    with pytest.raises(KeyError):
        new_moltype.RNA.make_seq(seq="UUU").matrix_distance("CAG", m)


def test_frac_same():
    """Sequence frac_same should return similarity between sequences"""
    s1 = new_moltype.RNA.make_seq(seq="ACGU")
    s2 = new_moltype.RNA.make_seq(seq="AACG")
    s3 = new_moltype.RNA.make_seq(seq="GG")
    s4 = new_moltype.RNA.make_seq(seq="A")
    e = new_moltype.RNA.make_seq(seq="")
    assert s1.frac_same(e) == 0
    assert s1.frac_same(s2) == 0.25
    assert s1.frac_same(s3) == 0
    assert s1.frac_same(s4) == 1.0  # note truncation


def test_frac_diff():
    """Sequence frac_diff should return difference between sequences"""
    s1 = new_moltype.RNA.make_seq(seq="ACGU")
    s2 = new_moltype.RNA.make_seq(seq="AACG")
    s3 = new_moltype.RNA.make_seq(seq="GG")
    s4 = new_moltype.RNA.make_seq(seq="A")
    e = new_moltype.RNA.make_seq(seq="")
    assert s1.frac_diff(e) == 0
    assert s1.frac_diff(s2) == 0.75
    assert s1.frac_diff(s3) == 1
    assert s1.frac_diff(s4) == 0  # note truncation


def test_frac_same_gaps():
    """Sequence frac_same_gaps should return similarity in gap positions"""
    s1 = new_moltype.RNA.make_seq(seq="AAAA")
    s2 = new_moltype.RNA.make_seq(seq="GGGG")
    s3 = new_moltype.RNA.make_seq(seq="----")
    s4 = new_moltype.RNA.make_seq(seq="A-A-")
    s5 = new_moltype.RNA.make_seq(seq="-G-G")
    s6 = new_moltype.RNA.make_seq(seq="UU--")
    s7 = new_moltype.RNA.make_seq(seq="-")
    s8 = new_moltype.RNA.make_seq(seq="GGG")
    e = new_moltype.RNA.make_seq(seq="")
    assert s1.frac_same_gaps(s1) == 1
    assert s1.frac_same_gaps(s2) == 1
    assert s1.frac_same_gaps(s3) == 0
    assert s1.frac_same_gaps(s4) == 0.5
    assert s1.frac_same_gaps(s5) == 0.5
    assert s1.frac_same_gaps(s6) == 0.5
    assert s1.frac_same_gaps(s7) == 0
    assert s1.frac_same_gaps(e) == 0
    assert s3.frac_same_gaps(s3) == 1
    assert s3.frac_same_gaps(s4) == 0.5
    assert s3.frac_same_gaps(s7) == 1.0
    assert e.frac_same_gaps(e) == 0.0
    assert s4.frac_same_gaps(s5) == 0.0
    assert s4.frac_same_gaps(s6) == 0.5
    assert numpy.allclose(s6.frac_same_gaps(s8), 2 / 3.0)


def test_frac_diff_gaps():
    """Sequence frac_diff_gaps should return difference in gap positions"""
    s1 = new_moltype.RNA.make_seq(seq="AAAA")
    s2 = new_moltype.RNA.make_seq(seq="GGGG")
    s3 = new_moltype.RNA.make_seq(seq="----")
    s4 = new_moltype.RNA.make_seq(seq="A-A-")
    s5 = new_moltype.RNA.make_seq(seq="-G-G")
    s6 = new_moltype.RNA.make_seq(seq="UU--")
    s7 = new_moltype.RNA.make_seq(seq="-")
    s8 = new_moltype.RNA.make_seq(seq="GGG")
    e = new_moltype.RNA.make_seq(seq="")
    assert s1.frac_diff_gaps(s1) == 0
    assert s1.frac_diff_gaps(s2) == 0
    assert s1.frac_diff_gaps(s3) == 1
    assert s1.frac_diff_gaps(s4) == 0.5
    assert s1.frac_diff_gaps(s5) == 0.5
    assert s1.frac_diff_gaps(s6) == 0.5
    assert s1.frac_diff_gaps(s7) == 1
    assert s1.frac_diff_gaps(e) == 0
    assert s3.frac_diff_gaps(s3) == 0
    assert s3.frac_diff_gaps(s4) == 0.5
    assert s3.frac_diff_gaps(s7) == 0.0
    assert e.frac_diff_gaps(e) == 0.0
    assert s4.frac_diff_gaps(s5) == 1.0
    assert s4.frac_diff_gaps(s6) == 0.5
    assert numpy.allclose(s6.frac_diff_gaps(s8), 1 / 3.0)


def test_frac_same_non_gaps():
    """Sequence frac_same_non_gaps should return similarities at non-gaps"""
    s1 = new_moltype.RNA.make_seq(seq="AAAA")
    s2 = new_moltype.RNA.make_seq(seq="AGGG")
    s3 = new_moltype.RNA.make_seq(seq="GGGG")
    s4 = new_moltype.RNA.make_seq(seq="AG--GA-G")
    s5 = new_moltype.RNA.make_seq(seq="CU--CU-C")
    s6 = new_moltype.RNA.make_seq(seq="AC--GC-G")
    s7 = new_moltype.RNA.make_seq(seq="--------")
    s8 = new_moltype.RNA.make_seq(seq="AAAA----")
    s9 = new_moltype.RNA.make_seq(seq="A-GG-A-C")
    e = new_moltype.RNA.make_seq(seq="")

    def test(x, y, z):
        assert numpy.allclose(x.frac_same_non_gaps(y), z)

    test(s1, s2, 0.25)
    test(s1, s3, 0)
    test(s2, s3, 0.75)
    test(s1, s4, 0.5)
    test(s4, s5, 0)
    test(s4, s6, 0.6)
    test(s4, s7, 0)
    test(s4, s8, 0.5)
    test(s4, s9, 2 / 3.0)
    test(e, s4, 0)


def test_frac_diffNonGaps():
    """Sequence frac_diff_non_gaps should return differences at non-gaps"""
    s1 = new_moltype.RNA.make_seq(seq="AAAA")
    s2 = new_moltype.RNA.make_seq(seq="AGGG")
    s3 = new_moltype.RNA.make_seq(seq="GGGG")
    s4 = new_moltype.RNA.make_seq(seq="AG--GA-G")
    s5 = new_moltype.RNA.make_seq(seq="CU--CU-C")
    s6 = new_moltype.RNA.make_seq(seq="AC--GC-G")
    s7 = new_moltype.RNA.make_seq(seq="--------")
    s8 = new_moltype.RNA.make_seq(seq="AAAA----")
    s9 = new_moltype.RNA.make_seq(seq="A-GG-A-C")
    e = new_moltype.RNA.make_seq(seq="")

    def test(x, y, z):
        assert numpy.allclose(x.frac_diff_non_gaps(y), z)

    test(s1, s2, 0.75)
    test(s1, s3, 1)
    test(s2, s3, 0.25)
    test(s1, s4, 0.5)
    test(s4, s5, 1)
    test(s4, s6, 0.4)
    test(s4, s7, 0)
    test(s4, s8, 0.5)
    test(s4, s9, 1 / 3.0)
    test(e, s4, 0)


def test_frac_similar():
    """Sequence frac_similar should return the fraction similarity"""
    transitions = dict.fromkeys(
        [
            ("A", "A"),
            ("A", "G"),
            ("G", "A"),
            ("G", "G"),
            ("U", "U"),
            ("U", "C"),
            ("C", "U"),
            ("C", "C"),
        ]
    )

    s1 = new_moltype.RNA.make_seq(seq="UCAGGCAA")
    s2 = new_moltype.RNA.make_seq(seq="CCAAAUGC")
    s3 = new_moltype.RNA.make_seq(seq="GGGGGGGG")
    e = new_moltype.RNA.make_seq(seq="")

    def test(x, y, z):
        assert numpy.allclose(x.frac_similar(y, transitions), z)

    test(e, e, 0)
    test(s1, e, 0)
    test(s1, s1, 1)
    test(s1, s2, 7.0 / 8)
    test(s1, s3, 5.0 / 8)
    test(s2, s3, 4.0 / 8)


def test_with_termini_unknown():
    """with_termini_unknown should reset termini to unknown char"""
    s1 = new_moltype.RNA.make_seq(seq="-?--AC--?-")
    s2 = new_moltype.RNA.make_seq(seq="AC")
    s3 = new_moltype.RNA.make_seq(seq="-----")
    assert s1.with_termini_unknown() == "????AC????"
    assert s2.with_termini_unknown() == "AC"
    assert s3.with_termini_unknown() == "?????"


def test_consistent_gap_degen_handling():
    """gap degen character should be treated consistently"""
    # the degen character '?' can be a gap, so when we strip either gaps or
    # degen characters it should be gone too
    raw_seq = "---??-??TC-GGCG-GCA-G-GC-?-C-TAN-GCGC-CCTC-AGGA?-???-??--"
    raw_ungapped = re.sub("[-?]", "", raw_seq)
    raw_no_ambigs = re.sub("[N?]+", "", raw_seq)
    dna = new_moltype.DNA.make_seq(seq=raw_seq)
    assert dna.degap() == raw_ungapped
    assert dna.strip_degenerate() == raw_no_ambigs
    assert dna.strip_bad_and_gaps() == raw_ungapped


def test_counts():
    """count motifs of different sizes, +/- ambiguities"""
    # test DNA seq
    orig = "AACCGGTTAN-T"
    seq = new_moltype.DNA.make_seq(seq=orig)
    # no gaps, no ambiguities
    got = seq.counts()
    expect = dict(A=3, C=2, G=2, T=3)
    assert dict(got) == expect
    # gaps allowed
    got = seq.counts(allow_gap=True)
    expect = dict(A=3, C=2, G=2, T=3)
    expect.update({"-": 1})
    assert dict(got) == expect
    # ambig allowed
    got = seq.counts(include_ambiguity=True)
    expect = dict(A=3, C=2, G=2, T=3, N=1)
    assert dict(got) == expect
    # ambig and gap allowed
    got = seq.counts(include_ambiguity=True, allow_gap=True)
    expect = dict(A=3, C=2, G=2, T=3, N=1)
    expect.update({"-": 1})
    assert dict(got) == expect

    # test DNA seq motif length of 2
    got = seq.counts(motif_length=2)
    expect = dict(AA=1, CC=1, GG=1, TT=1)
    assert dict(got) == expect
    # gap allowed
    got = seq.counts(motif_length=2, allow_gap=True)
    expect = dict(AA=1, CC=1, GG=1, TT=1)
    expect.update({"-T": 1})
    # ambig allowed
    got = seq.counts(motif_length=2, include_ambiguity=True)
    expect = dict(AA=1, CC=1, GG=1, TT=1, AN=1)
    assert dict(got) == expect
    # ambig and gap allowed
    got = seq.counts(motif_length=2, include_ambiguity=True, allow_gap=True)
    expect = dict(AA=1, CC=1, GG=1, TT=1, AN=1)
    expect.update({"-T": 1})
    assert dict(got) == expect

    # handle '?'
    orig = "AACCGGTTAN-T?"
    seq = new_moltype.DNA.make_seq(seq=orig)
    got = seq.counts()
    expect = dict(A=3, C=2, G=2, T=3)
    assert dict(got) == expect
    got = seq.counts(allow_gap=True, include_ambiguity=True)
    expect.update({"-": 1, "N": 1, "?": 1})
    assert dict(got) == expect


def test_strand_symmetry():
    """correctly compute test of strand symmetry"""
    from cogent3 import get_moltype
    from cogent3.core.alignment import Aligned

    seq = new_moltype.DNA.make_seq(seq="ACGGCTGAAGCGCTCCGGGTTTAAAACG")
    ssym = seq.strand_symmetry(motif_length=1)
    assert numpy.allclose(ssym.observed.array, [[7, 5], [7, 9]])
    assert numpy.allclose(ssym.expected.array, [[6, 6], [8, 8]])

    # RNA too
    seq = seq.to_rna()
    ssym = seq.strand_symmetry(motif_length=1)
    assert numpy.allclose(ssym.observed.array, [[7, 5], [7, 9]])

    # Aligned
    seq = new_moltype.DNA.make_seq(seq="ACGGCTGAAGCGCTCCGGGTTTAAAACG")
    m, s = seq.parse_out_gaps()
    seq = Aligned(m, s)
    ssym = seq.strand_symmetry(motif_length=1)
    assert numpy.allclose(ssym.observed.array, [[7, 5], [7, 9]])

    with pytest.raises(TypeError):
        text = get_moltype("text")
        m, s = text.make_seq(seq="ACGGCTGAAGCGCTCCGGGTTTAAAACG").parse_out_gaps()
        s.strand_symmetry(motif_length=1)

    # with motif_length=2
    seq = new_moltype.DNA.make_seq(
        seq="AC GG CT GA AG CG CT CC GG GT TT AA AA CG".replace(" ", "")
    )
    ssym = seq.strand_symmetry(motif_length=2)
    assert len(ssym.observed.keys()) <= 8
    assert numpy.allclose(ssym.observed["AA"].to_array(), [2, 1])
    assert numpy.allclose(ssym.observed["CC"].to_array(), [1, 2])


def test_is_annotated():
    """is_annotated operates correctly"""
    s = new_moltype.DNA.make_seq(seq="ACGGCTGAAGCGCTCCGGGTTTAAAACG")
    if hasattr(s, "annotation_db"):
        assert not s.is_annotated()
        _ = s.add_feature(biotype="gene", name="blah", spans=[(0, 10)])
        assert s.is_annotated()
    else:
        with pytest.raises(AttributeError):
            s.is_annotated()


def test_to_html():
    """produce correct html formatted text"""
    seq = new_moltype.DNA.make_seq(seq="ACGGTGGGGGGGGG")
    got = seq.to_html(wrap=50)
    # ensure balanced tags are in the txt
    for tag in ["<style>", "</style>", "<div", "</div>", "<table>", "</table>"]:
        assert tag in got

    seq_row = (
        '<tr><td class="label">None</td>'
        '<td><span class="A_dna">A</span>'
        '<span class="C_dna">C</span>'
        '<span class="G_dna">G</span>'
        '<span class="G_dna">G</span>'
        '<span class="T_dna">T</span>'
        '<span class="G_dna">G</span>'
        '<span class="G_dna">G</span>'
        '<span class="G_dna">G</span>'
        '<span class="G_dna">G</span>'
        '<span class="G_dna">G</span>'
        '<span class="G_dna">G</span>'
        '<span class="G_dna">G</span>'
        '<span class="G_dna">G</span>'
        '<span class="G_dna">G</span></td></tr>'
    )

    assert seq_row in got


def test_repr_html():
    """correctly uses set_repr and the environment variable settings"""
    token = 'class="label"'
    moltype = new_moltype.get_moltype("text")
    seq = new_sequence.Sequence(seq="AAAAA", moltype=moltype)

    orig = [l for l in seq._repr_html_().splitlines() if token in l][0]
    orig_num = len(re.findall(r"\bA\b", orig))
    assert orig_num == 5

    # using environment variable
    env_name = "COGENT3_ALIGNMENT_REPR_POLICY"
    os.environ[env_name] = "num_pos=2"
    got = [l for l in seq._repr_html_().splitlines() if token in l][0]
    got_num = len(re.findall(r"\bA\b", got))
    assert got_num == 2
    os.environ.pop(env_name, None)


def _check_mix_add(s1, s2):
    s1s2 = s1 + s2
    s2s1 = s2 + s1
    assert s1s2.name is None
    assert s2s1.name is None
    assert s1s2 == str(s1) + str(s2)
    assert s2s1 == str(s2) + str(s1)


def test_sequence_add_propogates_name():
    """name property correctly handled in sequence add"""
    moltype = new_moltype.get_moltype("text")
    a1 = new_sequence.Sequence(moltype=moltype, seq="AAA", name="1")
    a2 = new_sequence.Sequence(moltype=moltype, seq="CC", name="1")
    a = a1 + a2
    assert a.name == "1"
    assert a == "AAACC"

    b = new_sequence.Sequence(moltype=moltype, seq="GGGG", name="2")
    _check_mix_add(a1, b)
    c = new_sequence.Sequence(moltype=moltype, seq="TT")
    _check_mix_add(a1, c)

    e = "AA"
    be = b + e
    assert be.name is None
    assert be == str(b) + e


@pytest.mark.parametrize(
    "size, expect",
    [
        (60, ">even\nTCAGAT\n"),  # default block_size
        (2, ">even\nTC\nAG\nAT\n"),
        (4, ">even\nTCAG\nAT\n"),
    ],
)
def test_to_fasta_even(size, expect):
    mt = new_moltype.DNA
    seq = mt.make_seq(seq="TCAGAT", name="even")
    got = seq.to_fasta(block_size=size)
    assert got == expect


@pytest.mark.parametrize(
    "size, expect",
    [
        (60, ">odd\nTCAGATAAA\n"),  # default block_size
        (2, ">odd\nTC\nAG\nAT\nAA\nA\n"),
        (4, ">odd\nTCAG\nATAA\nA\n"),
    ],
)
def test_to_fasta_odd(size, expect):
    mt = new_moltype.DNA
    seq = mt.make_seq(seq="TCAGATAAA", name="odd")
    got = seq.to_fasta(block_size=size)
    assert got == expect


def test_to_phylip():
    mt = new_moltype.DNA
    seq = mt.make_seq(seq="TCAGAT", name="xyz")
    got = seq.to_phylip()
    expect = f"xyz{' '*27}TCAGAT"
    assert got == expect


@pytest.mark.parametrize("start", (None, 0, 1, 10, -1, -10))
@pytest.mark.parametrize("stop", (None, 10, 8, 1, 0, -1, -11))
@pytest.mark.parametrize("step", (None, 1, 2, -1, -2))
def test_seqview_initialisation(start, stop, step, bytes_alphabet):
    """Initialising a SeqView should work with range of provided values"""
    seq_data = "0123456789"
    slice_record = new_sequence.SliceRecord(
        start=start, stop=stop, step=step, parent_len=len(seq_data)
    )
    got = new_sequence.SeqView(
        parent=seq_data, slice_record=slice_record, alphabet=bytes_alphabet
    )
    expected = seq_data[start:stop:step]
    assert got.str_value == expected


@pytest.mark.parametrize("index", (-10, -5, 0, 5, 9))  # -10 and 9 are boundary
def test_seqview_index(index, bytes_alphabet):
    """SeqView with default values can be sliced with a single index, when within the length of the sequence"""
    seq_data = "0123456789"
    sv = new_sequence.SeqView(parent=seq_data, alphabet=bytes_alphabet)
    got = sv[index]
    expected = seq_data[index]
    assert got.str_value == expected
    assert len(got) == 1


def test_seqview_index_null(ascii_alphabet):
    "Indexing a SeqView of length 0 should return an IndexError"
    sv = new_sequence.SeqView(parent="", alphabet=ascii_alphabet)
    with pytest.raises(IndexError):
        _ = sv[0]


def test_seqview_step_0(bytes_alphabet):
    "Initialising or slicing a SeqView with a step of 0 should return an IndexError"
    sv = new_sequence.SeqView(parent="0123456789", alphabet=bytes_alphabet)
    with pytest.raises(ValueError):
        _ = sv[::0]
    with pytest.raises(ValueError):
        sr = new_sequence.SliceRecord(step=0, parent_len=10)
        _ = new_sequence.SeqView(
            parent="0123456789", alphabet=bytes_alphabet, slice_record=sr
        )


@pytest.mark.parametrize("start", (0, 2, 4))
def test_seqview_invalid_index(start, bytes_alphabet):
    "indexing out of bounds with a forward step should raise an IndexError"
    seq = "0123456789"
    length = abs(start - len(seq))
    pos_boundary_index = length
    neg_boundary_index = -length - 1

    sr = new_sequence.SliceRecord(start=start, parent_len=len(seq))
    sv = new_sequence.SeqView(parent=seq, slice_record=sr, alphabet=bytes_alphabet)
    with pytest.raises(IndexError):
        _ = sv[pos_boundary_index]
    with pytest.raises(IndexError):
        _ = sv[neg_boundary_index]


@pytest.mark.parametrize("start", (0, 2, 4))
def test_seqview_invalid_index_positive_step_gt_1(start, bytes_alphabet):
    "boundary condition for indexing out of bounds with a forward step greater than 1"
    seq = "0123456789"
    step = 2
    length = abs((start - len(seq)) // step)
    neg_boundary_index = -length - 1
    pos_boundary_index = length

    sr = new_sequence.SliceRecord(start=start, step=step, parent_len=len(seq))
    sv = new_sequence.SeqView(parent=seq, slice_record=sr, alphabet=bytes_alphabet)
    with pytest.raises(IndexError):
        _ = sv[pos_boundary_index]
    with pytest.raises(IndexError):
        _ = sv[neg_boundary_index]


@pytest.mark.parametrize("stop", (0, 2, -11))
def test_seqview_invalid_index_reverse_step(stop, bytes_alphabet):
    "boundary condition for indexing out of bounds with a reverse step"
    seq = "0123456789"
    step = -1
    start = len(seq)
    length = abs((start - stop) // step)
    neg_boundary_index = -length - 1
    pos_boundary_index = length

    sr = new_sequence.SliceRecord(
        start=start, stop=stop, step=step, parent_len=len(seq)
    )
    sv = new_sequence.SeqView(parent=seq, slice_record=sr, alphabet=bytes_alphabet)
    with pytest.raises(IndexError):
        _ = sv[pos_boundary_index]
    with pytest.raises(IndexError):
        _ = sv[neg_boundary_index]


@pytest.mark.parametrize("stop", (0, 2, -6))
def test_seqview_invalid_index_reverse_step_gt_1(stop, bytes_alphabet):
    "boundary condition for indexing out of bounds with a reverse step less than -1"
    seq = "0123456789"
    step = -2
    start = len(seq)
    length = abs((start - stop) // step)
    neg_boundary_index = -length - 1
    pos_boundary_index = length

    sr = new_sequence.SliceRecord(
        start=start, stop=stop, step=step, parent_len=len(seq)
    )
    sv = new_sequence.SeqView(parent=seq, slice_record=sr, alphabet=bytes_alphabet)
    with pytest.raises(IndexError):
        _ = sv[pos_boundary_index]
    with pytest.raises(IndexError):
        _ = sv[neg_boundary_index]


def test_seqview_slice_null(ascii_alphabet):
    sv = new_sequence.SeqView(parent="", alphabet=ascii_alphabet)
    assert len(sv) == 0
    got = sv[2:]
    assert len(got) == 0


def test_seqview_start_out_of_bounds(bytes_alphabet):
    "boundary condition for start index out of bounds"
    seq = "0123456789"
    init_start, init_stop, init_step = 2, 10, 1
    boundary = abs((init_start - init_stop) // init_step)
    sr = new_sequence.SliceRecord(
        start=init_start, stop=init_stop, step=init_step, parent_len=len(seq)
    )
    sv = new_sequence.SeqView(
        parent=seq,
        slice_record=sr,
        alphabet=bytes_alphabet,
    )
    got = sv[boundary::].str_value
    assert got == ""


def test_seqview_start_out_of_bounds_step_gt_1(bytes_alphabet):
    "boundary condition for start index out of bounds with step greater than 1"
    seq = "0123456789"
    init_start, init_stop, init_step = 2, 10, 2
    boundary = abs((init_start - init_stop) // init_step)
    sr = new_sequence.SliceRecord(
        start=init_start, stop=init_stop, step=init_step, parent_len=len(seq)
    )
    sv = new_sequence.SeqView(
        parent=seq,
        slice_record=sr,
        alphabet=bytes_alphabet,
    )
    got = sv[boundary::].str_value
    assert got == ""


def test_seqview_start_out_of_bounds_reverse_step(bytes_alphabet):
    "boundary condition for start index out of bounds with reverse step"
    seq = "0123456789"
    init_start, init_stop, init_step = 2, 10, -2
    boundary_pos = abs((init_start - init_stop) // init_step)
    boundary_neg = -abs((init_start - init_stop) // init_step) - 1

    sr = new_sequence.SliceRecord(
        start=init_start, stop=init_stop, step=init_step, parent_len=len(seq)
    )
    sv = new_sequence.SeqView(
        parent=seq,
        slice_record=sr,
        alphabet=bytes_alphabet,
    )

    assert sv[boundary_pos::].str_value == ""
    assert sv[boundary_neg::].str_value == ""


@pytest.mark.parametrize(
    "simple_slices",
    (
        slice(None, None, 1),
        slice(None, 3, None),
        slice(1, None, None),
        slice(1, 3, None),
        slice(None, None, None),
    ),
)
def test_seqview_defaults(simple_slices, bytes_alphabet):
    """SeqView should accept slices with all combinations of default parameters"""
    seq = "0123456789"
    got = new_sequence.SeqView(parent=seq, alphabet=bytes_alphabet)[simple_slices]
    expected = seq[simple_slices]
    assert got.str_value == expected


@pytest.mark.parametrize("index", (-8, -5, 0, 5, 8))
@pytest.mark.parametrize(
    "simple_slices",
    (
        slice(None, None, 1),
        slice(None, 10, None),
        slice(1, None, None),
        slice(1, 10, None),
        slice(1, 10, 1),
        slice(None, None, None),
    ),
)
def test_seqview_sliced_index(index, simple_slices, bytes_alphabet):
    """SeqView that has been sliced with default parameters, can then be indexed"""
    seq = "0123456789"
    sv = new_sequence.SeqView(parent=seq, alphabet=bytes_alphabet)
    got = sv[simple_slices][index]
    expected = seq[simple_slices][index]
    assert got.str_value == expected


@pytest.mark.parametrize("first_step", (1, 2, -1, -2))
@pytest.mark.parametrize("second_step", (1, 2, -1, -2))
def test_seqview_reverse_slice(first_step, second_step, bytes_alphabet):
    """subsequent slices may reverse the previous slice"""
    seq = "0123456789"
    sr = new_sequence.SliceRecord(step=first_step, parent_len=len(seq))
    sv = new_sequence.SeqView(parent=seq, slice_record=sr, alphabet=bytes_alphabet)
    got = sv[::second_step]
    expected = seq[::first_step][::second_step]
    assert got.str_value == expected


@pytest.mark.parametrize("seq", ("0123456789", "01234567890"))
@pytest.mark.parametrize("index", (-10, -4, 0, 6, 10))
@pytest.mark.parametrize("start", (None, 10, -1, -10))
@pytest.mark.parametrize("stop", (None, 9, -10, -11))
@pytest.mark.parametrize("step", (-1, -2))
def test_seqview_rev_sliced_index(index, start, stop, step, seq, bytes_alphabet):
    """SeqView that has been reverse sliced, can then be sliced with a single index"""
    seq_data = seq
    try:  # if python slicing raises an index error, we expect SeqView to also throw error
        expected = seq_data[start:stop:step][index]
    except IndexError:
        with pytest.raises(IndexError):
            sr = new_sequence.SliceRecord(
                start=start, stop=stop, step=step, parent_len=len(seq_data)
            )
            _ = new_sequence.SeqView(
                parent=seq_data,
                slice_record=sr,
                alphabet=bytes_alphabet,
            )[index].str_value
    else:  # if no index error, SeqView should match python slicing
        sr = new_sequence.SliceRecord(
            start=start, stop=stop, step=step, parent_len=len(seq_data)
        )
        got = new_sequence.SeqView(
            parent=seq_data, slice_record=sr, alphabet=bytes_alphabet
        )[index].str_value
        assert got == expected


@pytest.mark.parametrize("seq", ("0123456789", "012345678"))
@pytest.mark.parametrize("start", (None, 0, 1, 9, -1, -10))
@pytest.mark.parametrize("stop", (None, 0, 10, -7, -11))
@pytest.mark.parametrize("step", (1, 2, -1, -2))
def test_seqview_init_with_negatives(seq, start, stop, step, bytes_alphabet):
    "SeqView initialisation should handle any combination of positive and negative slices"
    sr = new_sequence.SliceRecord(
        start=start, stop=stop, step=step, parent_len=len(seq)
    )
    got = new_sequence.SeqView(parent=seq, slice_record=sr, alphabet=bytes_alphabet)
    expected = seq[start:stop:step]
    assert got.str_value == expected


@pytest.mark.parametrize("seq", ("0123456789", "012345678"))
@pytest.mark.parametrize("start", (None, 0, 1, 9, -1, -10))
@pytest.mark.parametrize("stop", (None, 0, 10, -7, -11))
@pytest.mark.parametrize("step", (1, 2, -1, -2))
def test_seqview_slice_with_negatives(seq, start, stop, step, bytes_alphabet):
    """SeqView should handle any combination of positive and negative slices"""
    sv = new_sequence.SeqView(parent=seq, alphabet=bytes_alphabet)
    got = sv[start:stop:step]
    expected = seq[start:stop:step]
    assert got.str_value == expected


@pytest.mark.parametrize("start", (None, 0, 2))
@pytest.mark.parametrize("stop", (None, 5, 7, 10))
@pytest.mark.parametrize("step", (1, 2))
@pytest.mark.parametrize("start_2", (None, 0, 1, 2))
@pytest.mark.parametrize("stop_2", (None, 2, 4, 10))
@pytest.mark.parametrize("step_2", (1, 2))
def test_subsequent_slice_forward(
    start, stop, step, start_2, stop_2, step_2, bytes_alphabet
):
    """SeqView should handle subsequent forward slice"""
    seq = "0123456789"
    sv = new_sequence.SeqView(parent=seq, alphabet=bytes_alphabet)
    got = sv[start:stop:step][start_2:stop_2:step_2]
    expected = seq[start:stop:step][start_2:stop_2:step_2]
    assert got.str_value == expected
    assert len(got) == len(expected)


@pytest.mark.parametrize(
    "slice_1, slice_2",
    (
        # WITH DEFAULTS
        # first stop -ve
        (slice(None, -3, None), slice(None, None, None)),
        # second stop -ve
        (slice(None, None, None), slice(None, -1, None)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(None, -3, None), slice(None, -5, None)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(None, -5, None), slice(None, -3, None)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(None, -3, None), slice(None, -8, None)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(None, -8, None), slice(None, -3, None)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(None, -2, None), slice(None, 7, None)),
        # first stop -ve, second stop +ve, second slice OUTSIDE first
        (slice(None, -6, None), slice(None, 7, None)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(None, 6, None), slice(None, -2, None)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(None, 6, None), slice(None, -7, None)),
        # WITH FIRST STEP > 1
        # first stop -ve
        (slice(None, -3, 2), slice(None, None, None)),
        # second stop -ve
        (slice(None, None, 2), slice(None, -1, None)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(None, -1, 2), slice(None, -3, None)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(None, -3, 2), slice(None, -2, None)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(None, -3, 2), slice(None, -8, None)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(None, -8, 2), slice(None, -3, None)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(None, -2, 2), slice(None, 3, None)),
        # first stop -ve, second stop +ve, second slice OVERLAP first
        (slice(None, -6, 2), slice(None, 7, None)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(None, 6, 2), slice(None, -2, None)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(None, 6, 2), slice(None, -7, None)),
        # WITH SECOND STEP > 1
        # first stop -ve
        (slice(None, -3, None), slice(None, None, 3)),
        # second stop -ve
        (slice(None, None, None), slice(None, -1, 3)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(None, -2, None), slice(None, -4, 2)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(None, -4, None), slice(None, -3, 2)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(None, -3, None), slice(None, -8, 2)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(None, -8, None), slice(None, -3, 2)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(None, -2, None), slice(None, 7, 2)),
        # first stop -ve, second stop +ve, second slice OVERLAP first
        (slice(None, -6, None), slice(None, 7, 2)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(None, 9, None), slice(None, -2, 3)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(None, 6, None), slice(None, -7, 3)),
        # WITH BOTH STEP > 1
        # first stop -ve
        (slice(None, -3, 2), slice(None, None, 3)),
        # second stop -ve
        (slice(None, None, 2), slice(None, -1, 3)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(None, -1, 3), slice(None, -2, 2)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(None, -2, 2), slice(None, -1, 2)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(None, -3, 3), slice(None, -8, 2)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(None, -8, 2), slice(None, -3, 2)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(None, -2, 3), slice(None, 7, 2)),
        # first stop -ve, second stop +ve, second slice OVERLAP first
        (slice(None, -3, 3), slice(None, 7, 2)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(None, 9, 2), slice(None, -1, 3)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(None, 6, 2), slice(None, -7, 3)),
        # NON-ZERO START
        # first stop -ve
        (slice(1, -3, 2), slice(None, None, 3)),
        # second stop -ve
        (slice(1, None, 2), slice(None, -1, 3)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(1, -1, 3), slice(None, -2, 2)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(1, -2, 2), slice(None, -1, 2)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(1, -3, 3), slice(None, -8, 2)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(1, -8, 2), slice(None, -3, 2)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(1, -2, 3), slice(None, 7, 2)),
        # first stop -ve, second stop +ve, second slice OVERLAP first
        (slice(1, -3, 3), slice(None, 7, 2)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(1, 10, 2), slice(None, -1, 3)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(1, 6, 2), slice(None, -7, 3)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
    ),
)
def test_subsequent_slice_neg_stop(slice_1, slice_2, ascii_alphabet):
    """SeqView should handle subsequence slices with >=1 negative stop values,
    subsequent slices may overlap or be within previous slices
    """
    seq_data = "abcdefghijk"
    sv = new_sequence.SeqView(parent=seq_data, alphabet=ascii_alphabet)
    assert sv[slice_1][slice_2].str_value == seq_data[slice_1][slice_2]


@pytest.mark.parametrize(
    "slice_1, slice_2",
    (
        # WITH DEFAULTS
        # first start -ve
        (slice(-6, None, None), slice(None, None, None)),
        # second start -ve
        (slice(None, None, None), slice(-6, None, None)),
        # both start -ve, (first < second), second slice WITHIN first
        (slice(-6, None, None), slice(-4, None, None)),
        # both start -ve, (first > second), second slice OUTSIDE first
        (slice(-4, None, None), slice(-6, None, None)),
        # first start -ve, second start +ve, second slice WITHIN first
        (slice(-8, None, None), slice(2, None, None)),
        # first start -ve, second start +ve, second slice OUTSIDE first
        (slice(-6, None, None), slice(7, None, None)),
        # first start +ve, second start -ve, second slice WITHIN first
        (slice(2, None, None), slice(-7, None, None)),
        # first start +ve, second start -ve, second slice OUTSIDE first
        (slice(5, None, None), slice(-6, None, None)),
        # WITH FIRST STEP > 1
        # first start -ve
        (slice(-6, None, 2), slice(None, None, None)),
        # second start -ve
        (slice(None, None, 2), slice(-6, None, None)),
        # both start -ve, (first < second), second slice WITHIN first
        (slice(-8, None, 2), slice(-6, None, None)),
        # both start -ve, (first > second), second slice OUTSIDE first
        (slice(-7, None, 2), slice(-9, None, None)),
        # first start -ve, second start +ve, second slice WITHIN first
        (slice(-9, None, 2), slice(2, None, None)),
        # first start -ve, second start +ve, second slice OUTSIDE first
        (slice(-6, None, 2), slice(7, None, None)),
        # first start +ve, second start -ve, second slice WITHIN first
        (slice(2, None, 2), slice(-7, None, None)),
        # first start +ve, second start -ve, second slice OUTSIDE first
        (slice(3, None, 2), slice(-9, None, None)),
        # WITH SECOND STEP > 1
        # first start -ve
        (slice(-6, None, None), slice(None, None, 2)),
        # second start -ve
        (slice(None, None, None), slice(-6, None, 2)),
        # both start -ve, (first < second), second slice WITHIN first
        (slice(-8, None, None), slice(-6, None, 2)),
        # both start -ve, (first > second), second slice OUTSIDE first
        (slice(-7, None, None), slice(-9, None, 2)),
        # first start -ve, second start +ve, second slice WITHIN first
        (slice(-9, None, None), slice(2, None, 2)),
        # first start -ve, second start +ve, second slice OUTSIDE first
        (slice(-6, None, None), slice(7, None, 2)),
        # first start +ve, second start -ve, second slice WITHIN first
        (slice(2, None, None), slice(-7, None, 2)),
        # first start +ve, second start -ve, second slice OUTSIDE first
        (slice(3, None, None), slice(-9, None, 2)),
        # WITH BOTH STEP > 1
        # first start -ve
        (slice(-6, None, 3), slice(None, None, 2)),
        # second start -ve
        (slice(None, None, 3), slice(-6, None, 2)),
        # both start -ve, (first < second), second slice WITHIN first
        (slice(-9, None, 3), slice(-7, None, 2)),
        # both start -ve, (first > second), second slice OUTSIDE first
        (slice(-7, None, 3), slice(-9, None, 2)),
        # first start -ve, second start +ve, second slice WITHIN first
        (slice(-9, None, 3), slice(2, None, 2)),
        # first start -ve, second start +ve, second slice OUTSIDE first
        (slice(-6, None, 2), slice(7, None, 2)),
        # first start +ve, second start -ve, second slice WITHIN first
        (slice(2, None, 3), slice(-7, None, 2)),
        # first start +ve, second start -ve, second slice OUTSIDE first
        (slice(3, None, 3), slice(-9, None, 2)),
        (slice(-9, 7, 3), slice(-2, None, None)),
    ),
)
def test_subsequent_slice_neg_start(slice_1, slice_2, ascii_alphabet):
    """SeqView should handle subsequence slices with >=1 negative start values,
    subsequent slices may or may not overlap or be within previous slices
    """
    seq_data = "abcdefghijk"
    sv = new_sequence.SeqView(parent=seq_data, alphabet=ascii_alphabet)
    assert sv[slice_1][slice_2].str_value == seq_data[slice_1][slice_2]


@pytest.mark.parametrize(
    "slice_1, slice_2",
    (
        # WITH DEFAULTS
        # first step -ve
        (slice(None, None, -1), slice(None, None, None)),
        # second step -ve
        (slice(None, None, None), slice(None, None, -1)),
        # both step -ve, start/stop -ve, second slice WITHIN first
        (slice(-1, -11, -2), slice(-1, -5, -3)),
        # both step -ve, start/stop -ve, second slice OUTSIDE first
        (slice(-1, -11, -2), slice(-1, -11, -3)),
        # both step -ve, start/stop +ve, second slice WITHIN first
        (slice(10, 0, -2), slice(5, 0, -3)),
        # both step -ve, start/stop +ve, second slice OUTSIDE first
        (slice(10, 0, -2), slice(10, 0, -3)),
        # first step -ve, second step +ve, second slice WITHIN first
        (slice(10, 0, -2), slice(1, 5, 2)),
        # first step -ve, second step +ve, second slice OUTSIDE first
        (slice(10, 0, -2), slice(0, 10, 2)),
        # first step +ve, second step -ve, second slice WITHIN first
        (slice(0, 10, 2), slice(4, 0, -2)),
        # first step +ve, second step -ve, second slice OUTSIDE first
        (slice(0, 10, 3), slice(10, 0, -2)),
        # first step -ve, second step +ve, second start/stop +ve
        (slice(10, 1, -1), slice(-8, 11, 2)),
        # first step -ve, second step +ve, second start/stop +ve
        (slice(10, 1, -1), slice(-19, 0, -2)),
    ),
)
def test_subsequent_slice_neg_step(slice_1, slice_2, ascii_alphabet):
    """SeqView should handle subsequence slices with negative step values,
    subsequent slices may overlap or be within previous slices
    """
    seq_data = "0123456789"
    sv = new_sequence.SeqView(parent=seq_data, alphabet=ascii_alphabet)
    assert sv[slice_1][slice_2].str_value == seq_data[slice_1][slice_2]


@pytest.mark.parametrize(
    "sub_slices_triple",
    (
        (slice(None, None, None), slice(None, None, None), slice(None, None, None)),
        (slice(1, 9, 1), slice(2, 8, 1), slice(3, 7, 1)),
        (slice(1, 9, 1), slice(2, 8, 1), slice(3, 9, 1)),
        (slice(1, 9, 1), slice(2, 8, 2), slice(3, 7, -3)),
    ),
)
def test_subslice_3(sub_slices_triple, ascii_alphabet):
    """SeqView should handle three subsequent slices"""
    seq_data = "abcdefghijk"
    sv = new_sequence.SeqView(parent=seq_data, alphabet=ascii_alphabet)
    slice_1, slice_2, slice_3 = sub_slices_triple
    assert (
        sv[slice_1][slice_2][slice_3].str_value == seq_data[slice_1][slice_2][slice_3]
    )


@pytest.mark.parametrize("start", (0, 2, -1))
@pytest.mark.parametrize("stop", (7, 10, -11))
@pytest.mark.parametrize("step", (1, -2))
@pytest.mark.parametrize("start_2", (0, 2, -8))
@pytest.mark.parametrize("stop_2", (2, 4))
@pytest.mark.parametrize("step_2", (2, -1))
@pytest.mark.parametrize("start_3", (0, 1, -6))
@pytest.mark.parametrize("stop_3", (4, 10, -10))
@pytest.mark.parametrize("step_3", (2, -2))
def test_triple_slice(
    integer_seq, start, stop, step, start_2, stop_2, step_2, start_3, stop_3, step_3
):
    """SeqView should handle subsequent forward slice"""
    seq = integer_seq.parent
    got = integer_seq[start:stop:step][start_2:stop_2:step_2][start_3:stop_3:step_3]
    expected = seq[start:stop:step][start_2:stop_2:step_2][start_3:stop_3:step_3]

    assert got.str_value == expected
    assert len(got) == len(expected)


def test_seqview_repr():
    alpha = new_moltype.DNA.most_degen_alphabet()
    # Short sequence, defaults
    seq = "ACGT"
    view = new_sequence.SeqView(parent=seq, alphabet=alpha)
    expected = f"SeqView(seqid=None, parent='ACGT', slice_record={view.slice_record!r})"
    assert repr(view) == expected

    # Long sequence
    seq = "ACGT" * 10
    view = new_sequence.SeqView(parent=seq, alphabet=alpha)
    expected = f"SeqView(seqid=None, parent='ACGTACGTAC...TACGT', slice_record={view.slice_record!r})"
    assert repr(view) == expected

    # Non-zero slice record
    seq = "ACGT" * 10
    sr = new_sequence.SliceRecord(start=5, stop=35, step=2, parent_len=len(seq))
    view = new_sequence.SeqView(parent=seq, alphabet=alpha, slice_record=sr)
    expected = f"SeqView(seqid=None, parent='ACGTACGTAC...TACGT', slice_record={view.slice_record!r})"
    assert repr(view) == expected

    # slice record with an offset
    seq = "ACGT"
    sr = new_sequence.SliceRecord(offset=5, parent_len=len(seq))
    view = new_sequence.SeqView(parent=seq, alphabet=alpha, slice_record=sr)
    expected = f"SeqView(seqid=None, parent='ACGT', slice_record={view.slice_record!r})"
    assert repr(view) == expected

    # seqid
    seq = "ACGT"
    view = new_sequence.SeqView(parent=seq, seqid="seq1", alphabet=alpha)
    expected = (
        f"SeqView(seqid='seq1', parent='ACGT', slice_record={view.slice_record!r})"
    )
    assert repr(view) == expected


def test_slice_record_repr():
    # defaults
    sr = new_sequence.SliceRecord(parent_len=20)
    expected = "SliceRecord(start=0, stop=20, step=1, parent_len=20, offset=0)"
    assert repr(sr) == expected

    sr = new_sequence.SliceRecord(start=1, stop=10, step=2, parent_len=20, offset=5)
    expected = "SliceRecord(start=1, stop=10, step=2, parent_len=20, offset=5)"
    assert repr(sr) == expected


@pytest.mark.parametrize("k", range(1, 7))
def test_iter_kmers_returns_generator(k):
    orig = "TCAGGA"
    mt = new_moltype.DNA
    r = mt.make_seq(seq=orig)
    got = r.iter_kmers(k=k)
    assert isinstance(got, typing.Generator)


def test_iter_kmers_empty():
    orig = ""
    mt = new_moltype.DNA
    r = mt.make_seq(seq=orig)
    got = r.iter_kmers(k=1)
    assert not list(got)


@pytest.mark.parametrize("k", (0, -1, 1.1))
def test_iter_kmers_handles_invalid(k):
    orig = ""
    mt = new_moltype.DNA
    r = mt.make_seq(seq=orig)
    got = r.iter_kmers(k=k)
    with pytest.raises(ValueError):
        list(got)


@pytest.mark.parametrize(
    "k, strict, expect",
    [
        (1, True, ["T", "C", "A", "G", "G", "A"]),
        (2, True, ["TC", "CA", "AG", "GG", "GA"]),
        (3, True, ["TCA", "CAG", "AGG", "GGA"]),
        (4, True, ["TCAG", "CAGG", "AGGA"]),
        (5, True, ["TCAGG", "CAGGA"]),
        (6, True, ["TCAGGA"]),
        (7, True, []),
        (1, False, ["T", "C", "A", "G", "G", "A", "N"]),
        (2, False, ["TC", "CA", "AG", "GG", "GA", "AN"]),
        (3, False, ["TCA", "CAG", "AGG", "GGA", "GAN"]),
        (4, False, ["TCAG", "CAGG", "AGGA", "GGAN"]),
        (5, False, ["TCAGG", "CAGGA", "AGGAN"]),
        (6, False, ["TCAGGA", "CAGGAN"]),
        (7, False, ["TCAGGAN"]),
        (8, False, []),
    ],
)
def test_get_kmers_strict_dna(k, strict, expect):
    orig = "TCAGGAN"
    mt = new_moltype.DNA
    r = mt.make_seq(seq=orig)
    got = r.get_kmers(k, strict=strict)
    assert got == expect


@pytest.mark.parametrize(
    "k, strict, expect",
    [
        (1, True, ["U", "C", "A", "G", "G", "A"]),
        (2, True, ["UC", "CA", "AG", "GG", "GA"]),
        (3, True, ["UCA", "CAG", "AGG", "GGA"]),
        (4, True, ["UCAG", "CAGG", "AGGA"]),
        (5, True, ["UCAGG", "CAGGA"]),
        (6, True, ["UCAGGA"]),
        (7, True, []),
        (1, False, ["U", "C", "A", "G", "G", "A", "N"]),
        (2, False, ["UC", "CA", "AG", "GG", "GA", "AN"]),
        (3, False, ["UCA", "CAG", "AGG", "GGA", "GAN"]),
        (4, False, ["UCAG", "CAGG", "AGGA", "GGAN"]),
        (5, False, ["UCAGG", "CAGGA", "AGGAN"]),
        (6, False, ["UCAGGA", "CAGGAN"]),
        (7, False, ["UCAGGAN"]),
        (8, False, []),
    ],
)
def test_get_kmers_strict_rna(k, strict, expect):
    orig = "UCAGGAN"
    mt = new_moltype.RNA
    r = mt.make_seq(seq=orig)
    got = r.get_kmers(k, strict=strict)
    assert got == expect


@pytest.mark.parametrize(
    "k, strict, expect",
    [
        (1, True, ["C", "E", "F", "G", "M", "N"]),
        (2, True, ["CE", "EF", "FG", "GM", "MN"]),
        (3, True, ["CEF", "EFG", "FGM", "GMN"]),
        (4, True, ["CEFG", "EFGM", "FGMN"]),
        (5, True, ["CEFGM", "EFGMN"]),
        (6, True, ["CEFGMN"]),
        (1, False, ["C", "E", "F", "G", "M", "N", "X"]),
        (2, False, ["CE", "EF", "FG", "GM", "MN", "NX"]),
        (3, False, ["CEF", "EFG", "FGM", "GMN", "MNX"]),
        (4, False, ["CEFG", "EFGM", "FGMN", "GMNX"]),
        (5, False, ["CEFGM", "EFGMN", "FGMNX"]),
        (6, False, ["CEFGMN", "EFGMNX"]),
        (7, False, ["CEFGMNX"]),
    ],
)
def test_get_kmers_strict_protein(k, strict, expect):
    orig = "CEFGMNX"
    mt = new_moltype.PROTEIN
    r = mt.make_seq(seq=orig)
    got = r.get_kmers(k, strict=strict)
    assert got == expect


@pytest.mark.parametrize(
    "k, strict, expect",
    [
        (1, True, ["T", "C", "A", "G", "A", "T"]),
        (2, True, ["TC", "CA", "GA", "AT"]),
        (3, True, ["TCA", "GAT"]),
        (4, True, []),
        (1, False, ["T", "C", "A", "-", "G", "A", "T"]),
        (2, False, ["TC", "CA", "A-", "-G", "GA", "AT"]),
        (3, False, ["TCA", "CA-", "A-G", "-GA", "GAT"]),
        (4, False, ["TCA-", "CA-G", "A-GA", "-GAT"]),
        (5, False, ["TCA-G", "CA-GA", "A-GAT"]),
        (6, False, ["TCA-GA", "CA-GAT"]),
        (7, False, ["TCA-GAT"]),
        (8, False, []),
    ],
)
def test_get_kmers_strict_dna_gaps(k, strict, expect):
    orig = "TCA-GAT"
    mt = new_moltype.DNA
    r = mt.make_seq(seq=orig)
    got = r.get_kmers(k, strict=strict)
    assert got == expect


@pytest.mark.parametrize(
    "k, strict, expect",
    [
        (1, True, ["U", "C", "A", "G", "A", "U"]),
        (2, True, ["UC", "CA", "GA", "AU"]),
        (3, True, ["UCA", "GAU"]),
        (4, True, []),
        (1, False, ["U", "C", "A", "-", "G", "A", "U"]),
        (2, False, ["UC", "CA", "A-", "-G", "GA", "AU"]),
        (3, False, ["UCA", "CA-", "A-G", "-GA", "GAU"]),
        (4, False, ["UCA-", "CA-G", "A-GA", "-GAU"]),
        (5, False, ["UCA-G", "CA-GA", "A-GAU"]),
        (6, False, ["UCA-GA", "CA-GAU"]),
        (7, False, ["UCA-GAU"]),
        (8, False, []),
    ],
)
def test_get_kmers_strict_rna_gaps(k, strict, expect):
    orig = "UCA-GAU"
    mt = new_moltype.RNA
    r = mt.make_seq(seq=orig)
    got = r.get_kmers(k, strict=strict)
    assert got == expect


@pytest.mark.parametrize(
    "k, strict, expect",
    [
        (1, True, ["C", "E", "F", "G", "M", "N"]),
        (2, True, ["CE", "EF", "GM", "MN"]),
        (3, True, ["CEF", "GMN"]),
        (4, True, []),
        (1, False, ["C", "E", "F", "-", "G", "M", "N"]),
        (2, False, ["CE", "EF", "F-", "-G", "GM", "MN"]),
        (3, False, ["CEF", "EF-", "F-G", "-GM", "GMN"]),
        (4, False, ["CEF-", "EF-G", "F-GM", "-GMN"]),
        (5, False, ["CEF-G", "EF-GM", "F-GMN"]),
        (6, False, ["CEF-GM", "EF-GMN"]),
        (7, False, ["CEF-GMN"]),
        (8, False, []),
    ],
)
def test_get_kmers_strict_protein_gaps(k, strict, expect):
    orig = "CEF-GMN"
    mt = new_moltype.PROTEIN
    r = mt.make_seq(seq=orig)
    got = r.get_kmers(k, strict=strict)
    assert got == expect


@pytest.mark.parametrize(
    "moltype", (new_moltype.DNA, new_moltype.RNA, new_moltype.PROTEIN)
)
def test_get_kmers_allgap(moltype):
    orig = "-------"
    expect = ["-", "-", "-", "-", "-", "-", "-"]
    r = moltype.make_seq(seq=orig)
    got = r.get_kmers(1, strict=False)
    assert got == expect

    expect = []
    got = r.get_kmers(1, strict=True)
    assert got == expect


@pytest.mark.parametrize(
    "k, strict, expect",
    [
        (1, True, ["G", "A", "T", "A"]),
        (2, True, ["GA", "TA"]),
        (3, True, []),
        (1, False, ["N", "G", "A", "S", "T", "A", "H"]),
        (2, False, ["NG", "GA", "AS", "ST", "TA", "AH"]),
        (3, False, ["NGA", "GAS", "AST", "STA", "TAH"]),
        (4, False, ["NGAS", "GAST", "ASTA", "STAH"]),
        (5, False, ["NGAST", "GASTA", "ASTAH"]),
        (6, False, ["NGASTA", "GASTAH"]),
        (7, False, ["NGASTAH"]),
        (8, False, []),
    ],
)
def test_get_kmers_mixed_ambiguities_dna(k, strict, expect):
    mt = new_moltype.DNA
    r = mt.make_seq(seq="NGASTAH")
    got = r.get_kmers(k, strict=strict)
    assert got == expect


@pytest.mark.parametrize(
    "k, strict, expect",
    [
        (1, True, ["G", "A", "U", "A"]),
        (2, True, ["GA", "UA"]),
        (3, True, []),
        (1, False, ["R", "G", "A", "W", "U", "A", "D"]),
        (2, False, ["RG", "GA", "AW", "WU", "UA", "AD"]),
        (3, False, ["RGA", "GAW", "AWU", "WUA", "UAD"]),
        (4, False, ["RGAW", "GAWU", "AWUA", "WUAD"]),
        (5, False, ["RGAWU", "GAWUA", "AWUAD"]),
        (6, False, ["RGAWUA", "GAWUAD"]),
        (7, False, ["RGAWUAD"]),
        (8, False, []),
    ],
)
def test_get_kmers_mixed_ambiguities_rna(k, strict, expect):
    mt = new_moltype.RNA
    r = mt.make_seq(seq="RGAWUAD")
    got = r.get_kmers(k, strict=strict)
    assert got == expect


@pytest.mark.parametrize(
    "k, strict, expect",
    [
        (1, True, ["Q", "M", "N", "R"]),
        (2, True, ["QM", "NR"]),
        (3, True, []),
        (1, False, ["B", "Q", "M", "X", "N", "R", "Z"]),
        (2, False, ["BQ", "QM", "MX", "XN", "NR", "RZ"]),
        (3, False, ["BQM", "QMX", "MXN", "XNR", "NRZ"]),
        (4, False, ["BQMX", "QMXN", "MXNR", "XNRZ"]),
        (5, False, ["BQMXN", "QMXNR", "MXNRZ"]),
        (6, False, ["BQMXNR", "QMXNRZ"]),
        (7, False, ["BQMXNRZ"]),
        (8, False, []),
    ],
)
def test_get_kmers_mixed_ambiguities_protein(k, strict, expect):
    mt = new_moltype.PROTEIN
    r = mt.make_seq(seq="BQMXNRZ")
    got = r.get_kmers(k, strict=strict)
    assert got == expect


@pytest.fixture(scope="function")
def one_seq():
    return new_moltype.DNA.make_seq(seq="AACCTGGAACC")


def test_seq_repr(one_seq):
    pat = re.compile("[ACGT]+")
    expect = str(one_seq)
    seq = one_seq

    got = pat.findall(repr(seq))[0]
    assert expect.startswith(got), (expect, got)


def test_seq_repr_rc(one_seq):
    pat = re.compile("[ACGT]+")
    dna = one_seq.moltype
    expect = dna.rc(str(one_seq))
    seq = one_seq.rc()

    got = pat.findall(repr(seq))[0]
    assert expect.startswith(got), (expect, got)


def test_annotation_from_slice_with_stride():
    seq = new_moltype.DNA.make_seq(seq="AAACGCGCGAAAAAAA", name="s1")
    seq.add_feature(biotype="exon", name="ex1", spans=[(3, 9)])
    f = list(seq.get_features(name="ex1"))[0]
    assert str(f.get_slice()) == "CGCGCG"
    s1 = seq[1::2]
    f = list(s1.get_features(name="ex1"))[0]
    assert str(f.get_slice()) == "CCC"


def test_absolute_position_base_cases(one_seq):
    """with no offset or view, the absolute index should remain unchanged"""
    got = one_seq._seq.slice_record.absolute_position(5)
    assert got == 5

    # an index outside the range of the sequence should raise an IndexError
    with pytest.raises(IndexError):
        one_seq._seq.slice_record.absolute_position(20)

    with pytest.raises(IndexError):
        one_seq._seq.slice_record.absolute_position(-20)


def test_absolute_position_positive(one_seq):
    # with an offset, the abs index should be offset + index
    one_seq._seq = one_seq._seq.with_offset(2)
    view = one_seq._seq
    got = view.slice_record.absolute_position(2)
    assert got == 2 + 2

    # with an offset and start, the abs index should be offset + start + index
    seq = one_seq[2::]
    view = seq._seq  # todo: do we want the annotation_offset to be preserved when slicing? I think yes
    got = view.slice_record.absolute_position(2)
    assert got == 2 + 2 + 2

    # with an offset, start and step, the abs index should be offset + start + index * step
    seq = one_seq[2::2]
    got = seq._seq.slice_record.absolute_position(2)
    assert got == 2 + 2 + 2 * 2


def test_relative_position_base_cases(one_seq):
    """with no offset or view, the absolute index should remain unchanged"""
    got = one_seq._seq.slice_record.relative_position(5)
    assert got == 5

    # a -ve index  should raise an IndexError
    with pytest.raises(IndexError):
        one_seq._seq.slice_record.relative_position(-5)


def test_relative_position(integer_seq):
    """This test checks if the method returns the correct relative positions when
    the given index precedes or exceeds the range of the SeqView."""

    view = integer_seq[1:9:]
    # view = "12345678"
    got = view.slice_record.relative_position(0)
    # precedes the view, so should return -1
    assert got == -1
    # exceeds the view, but still returns a value
    got = view.slice_record.relative_position(10)
    assert got == 9


def test_relative_position_step_GT_one(integer_seq):
    """This test checks if the method returns the correct relative positions when
    the given index precedes or exceeds the range of the SeqView with a step greater than one.
    """

    # precedes the view, with step > 1
    view = integer_seq[2:7:2]
    # view = "246", precedes the view by 1 step
    got = view.slice_record.relative_position(0)
    assert got == -1
    # precedes the view by 0.5 step, default behaviour is to round up to 0
    got = view.slice_record.relative_position(1)
    assert got == 0
    # exceeds the view by two steps, len(view) + 2 = 4
    got = view.slice_record.relative_position(10)
    assert got == 4


@pytest.mark.parametrize("sliced", (False, True))
@pytest.mark.parametrize("rev", (False, True))
def test_seqview_copy(sliced, rev, integer_seq):
    raw_data = integer_seq.parent
    integer_seq = integer_seq[::-1] if rev else integer_seq
    raw_data = raw_data[::-1] if rev else raw_data

    slice_start = 2
    slice_end = 4
    sv = integer_seq[slice_start:slice_end]
    copied = sv.copy(sliced=sliced)

    assert copied.str_value == raw_data[slice_start:slice_end]
    assert copied.is_reversed == integer_seq.is_reversed
    assert (
        sliced and copied.parent is not sv.parent or copied.parent is integer_seq.parent
    )


def test_relative_position_with_remainder(integer_seq):
    """tests relative_position when the index given is excluded from the view as it falls on
    a position that is 'stepped over'"""
    view = integer_seq[1:9:2]
    # view = "1357"
    got = view.slice_record.relative_position(2)
    # 2 is stepped over in the view, so we return the index of 3 (which is 1)
    assert got == 1

    # setting the arg stop=True will adjust to the largest number, smaller than the given abs value, that is in the view
    got = view.slice_record.relative_position(8, stop=True)
    # 8 is excluded from the view, so we return the index of 7 (which is 3)
    assert got == 3


@pytest.mark.parametrize("value", (0, 3))
@pytest.mark.parametrize("offset", (None, 1, 2))
@pytest.mark.parametrize("start", (None, 1, 2))
@pytest.mark.parametrize("stop", (None, 10, 11))
@pytest.mark.parametrize("step", (None, 1, 2))
def test_absolute_relative_roundtrip(one_seq, value, offset, start, stop, step):
    # a round trip from relative to absolute then from absolute to relative, should return the same value we began with
    seq = one_seq[start:stop:step]
    seq._seq.with_offset(offset or 0)
    abs_val = seq._seq.slice_record.absolute_position(value)
    rel_val = seq._seq.slice_record.relative_position(abs_val)
    assert rel_val == value


@pytest.mark.parametrize("value", (0, 2))
@pytest.mark.parametrize("offset", (None, 1, 2))
@pytest.mark.parametrize("start", (None, -1, -2))
@pytest.mark.parametrize("stop", (None, -10))
@pytest.mark.parametrize("step", (-1, -2))
def test_absolute_relative_roundtrip_reverse(
    integer_seq, value, offset, start, stop, step
):
    # a round trip from relative to absolute then from absolute to relative, should return the same value we began with
    view = integer_seq[start:stop:step]
    view.slice_record.offset = offset or 0
    abs_val = view.slice_record.absolute_position(value)
    rel_val = view.slice_record.relative_position(abs_val)
    assert view.slice_record.offset == (offset or 0)
    assert (view[rel_val]).str_value == view[value].str_value


def test_annotate_gff_nested_features(DATA_DIR):
    """correctly annotate a sequence with nested features"""
    # the synthetic example
    #          1111111111222222222333333333334
    # 1234567890123456789012345678901234567890
    #  **** biological_region
    #                                     ** biological_region
    #                                       * biological_region
    #      *******************************  gene
    #         *********************   mRNA
    #            *********            exon
    #                       *****     exon
    # ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC...
    seq = new_moltype.DNA.make_seq(
        seq="ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC", name="22"
    )
    gff3_path = DATA_DIR / "ensembl_sample.gff3"
    # TODO: directly assign an annotation_db, annotate_from_gff to be discontinued
    seq.annotate_from_gff(gff3_path)
    # we have 8 records in the gff file
    assert seq.annotation_db.num_matches() == 8

    # get the gene and check it has a single annotation and that
    # its slice is correct
    ann = list(seq.get_features(biotype="gene"))
    assert len(ann) == 1
    ann_seq = ann[0].get_slice()
    assert str(ann_seq) == "GGAAAATTTTTTTTTAAGGGGGAAAAAAAAA"
    # the gene has 1 transcript
    gene = ann[0]
    mrna = list(gene.get_children(biotype="mRNA"))
    assert len(mrna) == 1
    mrna = mrna[0]
    ann_seq = mrna.get_slice()
    assert str(ann_seq) == "AAAATTTTTTTTTAAGGGGGAAA"

    # the transcript has 2 exons, from the parent feature
    exons = list(mrna.get_children(biotype="exon"))
    assert len(exons) == 2
    # or the sequence
    ann = list(seq.get_features(biotype="exon"))
    assert len(ann) == 2
    exon_seqs = ("TTTTTTTTT", "GGGGG")
    assert tuple(str(ex.get_slice()) for ex in exons) == exon_seqs


def test_to_moltype_dna():
    """to_moltype("dna") ensures conversion from T to U"""
    seq = new_moltype.DNA.make_seq(seq="AAAAGGGGTTT", name="seq1")
    rna = seq.to_moltype("rna")
    assert "T" not in rna


def test_to_moltype_rna():
    """to_moltype("rna") ensures conversion from U to T"""
    seq = new_moltype.RNA.make_seq(seq="AAAAGGGGUUU", name="seq1")
    rna = seq.to_moltype("dna")
    assert "U" not in rna


def test_to_rich_dict():
    """Sequence to_dict works"""
    dna = new_moltype.DNA
    r = dna.make_seq(seq="AAGGCC", name="seq1")
    got = r.to_rich_dict()
    seq = new_sequence.SeqView(
        parent="AAGGCC", seqid="seq1", alphabet=dna.most_degen_alphabet()
    ).to_rich_dict()

    expect = {
        "name": "seq1",
        "seq": seq,
        "moltype": r.moltype.label,
        "info": None,
        "type": get_object_provenance(r),
        "version": __version__,
        "annotation_offset": 0,
    }

    assert got == expect


def test_sequence_to_json():
    """to_json roundtrip recreates to_dict"""
    dna = new_moltype.DNA
    r = dna.make_seq(seq="AAGGCC", name="seq1")
    got = json.loads(r.to_json())
    seq = new_sequence.SeqView(
        parent="AAGGCC", seqid="seq1", alphabet=dna.most_degen_alphabet()
    ).to_rich_dict()

    expect = {
        "name": "seq1",
        "seq": seq,
        "moltype": dna.label,
        "info": None,
        "type": get_object_provenance(r),
        "version": __version__,
        "annotation_offset": 0,
    }

    assert got == expect


@pytest.mark.parametrize("coord", ("start", "stop"))
def test_seqview_to_rich_dict(coord, dna_alphabet):
    parent = "ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC"
    sv = new_sequence.SeqView(parent=parent, alphabet=dna_alphabet)
    plus = sv.to_rich_dict()
    minus = sv[::-1].to_rich_dict()
    plus = plus.pop("init_args")
    minus = minus.pop("init_args")
    assert plus.pop("parent") == minus.pop("parent")
    assert (
        plus["slice_record"]["init_args"]["step"]
        == -minus["slice_record"]["init_args"]["step"]
    )
    assert coord not in plus
    assert coord not in minus


@pytest.mark.parametrize("reverse", (False, True))
def test_sliced_seqview_rich_dict(reverse, dna_alphabet):
    parent = "ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC"
    sl = slice(2, 13)
    sv = new_sequence.SeqView(parent=parent, alphabet=dna_alphabet)[sl]
    sv = sv[::-1] if reverse else sv
    rd = sv.to_rich_dict()
    assert rd["init_args"]["parent"] == parent[sl]


@pytest.mark.parametrize(
    "sl",
    (
        slice(2, 5, 1),  # positive indices, positive step
        slice(-8, -5, 1),  # negative indices, positive step
        slice(4, 1, -1),  # positive indices, negative step
        slice(-6, -9, -1),  # negative indices, negative step
    ),
)
@pytest.mark.parametrize("offset", (4, 0))
def test_parent_start_stop(sl, offset, ascii_alphabet):
    data = "0123456789"
    # check our slice matches the expectation for rest of test
    expect = "234" if sl.step > 0 else "432"
    sv = new_sequence.SeqView(parent=data, alphabet=ascii_alphabet, offset=offset)
    sv = sv[sl]
    assert sv.offset == offset
    assert sv.str_value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.slice_record.parent_start, sv.slice_record.parent_stop) == (
        2 + offset,
        5 + offset,
    )


@pytest.mark.parametrize(
    "sl",
    (
        slice(None, None, 1),  # slice whole sequence plus strand
        slice(None, None, -1),  # slice whole sequence minus strand
    ),
)
def test_parent_start_stop_limits(sl, ascii_alphabet):
    data = "0123456789"
    # check our slice matches the expectation for rest of test
    expect = data[sl]
    sv = new_sequence.SeqView(parent=data, alphabet=ascii_alphabet)
    sv = sv[sl]
    assert sv.str_value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.slice_record.parent_start, sv.slice_record.parent_stop) == (0, 10)


@pytest.mark.parametrize("rev", (False, True))
def test_parent_start_stop_empty(rev, ascii_alphabet):
    data = "0123456789"
    # check our slice matches the expectation for rest of test
    expect = ""
    sv = new_sequence.SeqView(parent=data, alphabet=ascii_alphabet)
    sv = sv[0 : 0 : -1 if rev else 1]
    assert sv.str_value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.slice_record.parent_start, sv.slice_record.parent_stop) == (0, 0)


@pytest.mark.parametrize("rev", (False, True))
@pytest.mark.parametrize("index", range(9))
def test_parent_start_stop_singletons(index, rev, ascii_alphabet):
    data = "0123456789"
    start, stop = (-(10 - index), -(10 - index + 1)) if rev else (index, index + 1)
    sl = slice(start, stop, -1 if rev else 1)
    # check our slice matches the expectation for rest of test
    expect = data[sl]
    sv = new_sequence.SeqView(parent=data, alphabet=ascii_alphabet)
    sv = sv[sl]
    assert sv.str_value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.slice_record.parent_start, sv.slice_record.parent_stop) == (
        index,
        index + 1,
    )


def test_get_drawable(DATA_DIR):
    seq = cogent3.load_seq(DATA_DIR / "annotated_seq.gb")
    seq = seq[2000:4000]
    biotypes = "CDS", "gene", "mRNA"
    for feat in seq.get_features(biotype=biotypes, allow_partial=True):
        draw = feat.get_drawable()
        assert "(incomplete)" in draw.text

    full = seq.get_drawable(biotype=biotypes)
    # should only include elements that overlap the segment
    assert len(full.traces) == len(biotypes)
    # and their names should indicate they're incomplete
    for trace in full.traces:
        assert "(incomplete)" in trace.text


@pytest.mark.parametrize("gc,seq", ((1, "TCCTGA"), (1, "ACGTAA---"), (2, "TCCAGG")))
def test_has_terminal_stop_true(gc, seq):
    gc = new_genetic_code.get_code(gc)
    seq = new_moltype.DNA.make_seq(seq=seq)
    assert seq.has_terminal_stop(gc=gc)


@pytest.mark.parametrize(
    "gc,seq", ((1, "TCCAGG"), (2, "TCCAAA"), (1, "CCTGA"), (2, "CCAGG"))
)
def test_has_terminal_stop_false(gc, seq):
    gc = new_genetic_code.get_code(gc)
    seq = new_moltype.DNA.make_seq(seq=seq)
    assert not seq.has_terminal_stop(gc=gc)


def test_has_terminal_stop_strict():
    gc = new_genetic_code.get_code(1)
    seq = new_moltype.DNA.make_seq(seq="TCCAG")
    with pytest.raises(new_alphabet.AlphabetError):
        seq.has_terminal_stop(gc=gc, strict=True)


@pytest.mark.parametrize(
    "gc,seq",
    (
        (2, "TCCAGG"),
        (1, "TAATGA"),
        (1, "ACGTGA---"),
        (1, "--AT-CTGA"),
    ),
)
def test_trim_terminal_stop_true(gc, seq):
    gc = new_genetic_code.get_code(gc)
    expect = re.sub("(TGA|AGG)(?=[-]*$)", "---" if "-" in seq else "", seq)

    seq = new_moltype.DNA.make_seq(seq=seq)
    got = str(seq.trim_stop_codon(gc=gc))
    assert got == expect


@pytest.mark.parametrize("gc,seq", ((1, "T?CTGC"), (2, "TCCAAG")))
def test_trim_terminal_stop_nostop(gc, seq):
    gc = new_genetic_code.get_code(gc)
    seq = new_moltype.DNA.make_seq(seq=seq)
    got = seq.trim_stop_codon(gc=gc)
    assert str(got) == str(seq)
    # since there's no stop, we just return the same object
    assert got is seq


@pytest.mark.parametrize(
    "gc,seq", ((1, "TCCAGG"), (2, "TCCAAA"), (1, "CCTGA"), (2, "CCAGG"))
)
def test_trim_terminal_stop_false(gc, seq):
    gc = new_genetic_code.get_code(gc)
    seq = new_moltype.DNA.make_seq(seq=seq)
    assert str(seq.trim_stop_codon(gc=gc)) == str(seq)


def test_trim_terminal_stop_strict():
    gc = new_genetic_code.get_code(1)
    seq = new_moltype.DNA.make_seq(seq="TCCAG")
    with pytest.raises(new_alphabet.AlphabetError):
        seq.trim_stop_codon(gc=gc, strict=True)


@pytest.mark.parametrize("cast", (int, numpy.int32, numpy.int64, numpy.uint8))
def test_index_a_seq(cast):
    seq = new_moltype.DNA.make_seq(seq="TCCAG")
    got = seq[cast(1)]
    assert isinstance(got, new_sequence.Sequence)


@pytest.mark.parametrize("cast", (float, numpy.float32))
def test_index_a_seq_float_fail(cast):
    seq = new_moltype.DNA.make_seq(seq="TCCAG")
    index = cast(1)
    with pytest.raises(TypeError):
        seq[index]  # pylint: disable=W0104


@pytest.mark.parametrize("moltype", ("dna", "protein"))
def test_same_moltype(moltype):
    moltype = new_moltype.get_moltype(moltype)
    seq = moltype.make_seq(seq="TCCAG")
    got = seq.to_moltype(moltype)
    assert got is seq


def test_gapped_by_map_segment_iter():
    moltype = new_moltype.DNA
    m, seq = moltype.make_seq(seq="-TCC--AG").parse_out_gaps()
    g = list(seq.gapped_by_map_segment_iter(m, allow_gaps=True, recode_gaps=False))
    assert g == ["-", "TCC", "--", "AG"]


@pytest.mark.parametrize("rev", (False, True))
@pytest.mark.parametrize("sliced", (False, True))
@pytest.mark.parametrize("start_stop", ((None, None), (3, 7)))
def test_copied_parent_coordinates(sliced, rev, start_stop):
    orig_name = "orig"
    seq = new_moltype.DNA.make_seq(seq="ACGGTGGGAC", name=orig_name)
    start, stop = start_stop
    start = start or 0
    stop = stop or len(seq)
    sl = slice(start, stop)
    seq = seq[sl]
    sliced_name = "sliced"
    seq.name = sliced_name
    assert seq.name == sliced_name
    seq = seq.rc() if rev else seq
    copied = seq.copy(sliced=sliced)
    assert copied.name == sliced_name
    # matches original
    assert copied.parent_coordinates() == seq.parent_coordinates()
    # and expected -- the coordinate name always reflects the underlying sequence
    assert copied.parent_coordinates() == (orig_name, start, stop, -1 if rev else 1)


@pytest.mark.parametrize("rev", (False, True))
def test_parent_coordinates(one_seq, rev):
    seq = one_seq[1:1]
    seq = seq.rc() if rev else seq
    seq.name = "sliced"  # this assignment does not affect the
    # note that when a sequence has zero length, the parent seqid is None
    assert seq.parent_coordinates() == (None, 0, 0, 1)


@pytest.mark.parametrize("cls", (str, bytes))
def test_coerce_to_seqview_str_bytes(cls, dna_alphabet):
    seq = "AC--GGTGGGAC"
    seqid = "seq1"
    s = bytes(seq, "utf8") if cls == bytes else seq
    got = new_sequence._coerce_to_seqview(s, seqid, alphabet=dna_alphabet, offset=0)
    assert got.str_value == seq
    assert isinstance(got, new_sequence.SeqView)


def test_coerce_to_seqview_sequence(dna_alphabet):
    seq = "AC--GGTGGGAC"
    seqid = "seq1"
    got = new_sequence._coerce_to_seqview(
        new_moltype.DNA.make_seq(seq=seq), seqid, alphabet=dna_alphabet, offset=0
    )
    assert got.str_value == seq
    assert isinstance(got, new_sequence.SeqView)


def test_coerce_to_seqview_already_seqview(dna_alphabet):
    seq = "AC--GGTGGGAC"
    seqid = "seq1"
    got = new_sequence._coerce_to_seqview(
        new_sequence.SeqView(parent=seq, alphabet=dna_alphabet),
        seqid,
        alphabet=dna_alphabet,
        offset=0,
    )
    assert got.str_value == seq
    assert isinstance(got, new_sequence.SeqView)


def test_seqview_seqid(dna_alphabet):
    sv = new_sequence.SeqView(parent="ACGGTGGGAC", alphabet=dna_alphabet)
    assert sv.seqid is None

    sv = new_sequence.SeqView(parent="ACGGTGGGAC", seqid="seq1", alphabet=dna_alphabet)
    assert sv.seqid == "seq1"


def test_seqview_to_rich_dict_seqid(dna_alphabet):
    sv = new_sequence.SeqView(parent="ACGGTGGGAC", seqid="seq1", alphabet=dna_alphabet)
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seqid"] == "seq1"

    sv = new_sequence.SeqView(parent="ACGGTGGGAC", alphabet=dna_alphabet)
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seqid"] is None


def test_seqview_slice_propagates_seqid(dna_alphabet):
    sv = new_sequence.SeqView(parent="ACGGTGGGAC", seqid="seq1", alphabet=dna_alphabet)
    sliced_sv = sv[1:8:2]
    assert sliced_sv.seqid == "seq1"

    copied_sv = sliced_sv.copy(sliced=False)
    assert copied_sv.seqid == "seq1"

    copied_sliced_sv = sliced_sv.copy(sliced=True)
    assert copied_sliced_sv.seqid == "seq1"


def test_sequences_propogates_seqid():
    # creating a name Sequence propagates the seqid to the SeqView.
    seq = new_moltype.DNA.make_seq(seq="ACGGTGGGAC", name="seq1")
    assert seq._seq.seqid == "seq1"

    # renaming the Sequence doesnt change the seqid of the SeqView.
    seq.name = "seq2"
    assert seq.name == "seq2"
    assert seq._seq.seqid == "seq1"


def test_sequences_propogates_seqid_seqview(dna_alphabet):
    # creating a Sequence with a seqview does not change the seqid of the SeqView.
    seq = new_moltype.DNA.make_seq(
        seq=new_sequence.SeqView(
            parent="ACGGTGGGAC", seqid="parent_name", alphabet=dna_alphabet
        ),
        name="seq_name",
    )
    assert seq.name == "seq_name"
    assert seq._seq.seqid == "parent_name"

    # creating a Sequence with an unnamed seqview does not name the SeqView.
    seq = new_moltype.DNA.make_seq(
        seq=new_sequence.SeqView(parent="ACGGTGGGAC", alphabet=dna_alphabet),
        name="seq_name",
    )
    assert seq.name == "seq_name"
    assert seq._seq.seqid is None


def test_make_seq_assigns_to_seqview():
    seq = new_moltype.DNA.make_seq(seq="ACGT", name="s1")
    assert seq.name == seq._seq.seqid == "s1"


def test_empty_seqview_translate_position(dna_alphabet):
    sv = new_sequence.SeqView(parent="", alphabet=dna_alphabet)
    assert sv.slice_record.absolute_position(0) == 0
    assert sv.slice_record.relative_position(0) == 0


@pytest.mark.parametrize("start", (None, 0, 1, 10, -1, -10))
@pytest.mark.parametrize("stop", (None, 10, 8, 1, 0, -1, -11))
@pytest.mark.parametrize("step", (None, 1, 2, -1, -2))
@pytest.mark.parametrize("length", (1, 8, 999))
def test_seqview_seq_len_init(start, stop, step, length, dna_alphabet):
    # seq_len is length of seq when None
    seq_data = "A" * length
    sr = new_sequence.SliceRecord(start=start, stop=stop, step=step, parent_len=length)
    sv = new_sequence.SeqView(parent=seq_data, slice_record=sr, alphabet=dna_alphabet)
    expect = len(seq_data)
    # Check property and slot
    assert sv.parent_len == expect
    assert sv._parent_len == expect


def test_seqview_copy_propagates_seq_len(dna_alphabet):
    seq = "ACGGTGGGAC"
    sv = new_sequence.SeqView(parent=seq, alphabet=dna_alphabet)
    copied = sv.copy()
    assert copied.parent_len == len(seq)


def test_seqview_seq_len_modified_seq(dna_alphabet):
    seq = "ACGGTGGGAC"
    sv = new_sequence.SeqView(parent=seq, alphabet=dna_alphabet)

    sv.parent = "ATGC"  # this should not modify seq_len
    assert sv.parent_len == len(seq)


@pytest.mark.parametrize("offset", [0, 4])
def test_seqview_with_offset(offset, dna_alphabet):
    seq = "ACGGTGGGAC"
    sv = new_sequence.SeqView(parent=seq, alphabet=dna_alphabet)
    got = sv.with_offset(offset)
    assert got is not sv
    assert got.offset == offset


@pytest.mark.parametrize("offset", [0, 4])
def test_seqview_with_offset_fails(offset, dna_alphabet):
    seq = "ACGGTGGGAC"
    sv = new_sequence.SeqView(parent=seq, alphabet=dna_alphabet, offset=1)
    with pytest.raises(ValueError):
        _ = sv.with_offset(offset)


def test_sequence_str_bytes_array():
    data = "ACGGTGGGAC"
    seq = new_moltype.DNA.make_seq(seq=data)
    assert str(seq) == data
    assert bytes(seq) == data.encode("utf8")
    assert numpy.array_equal(
        numpy.array(seq), new_moltype.DNA.alphabet.to_indices(data)
    )


@pytest.mark.parametrize("seq,rc", (("ATGTTT", False), ("AAACAT", True)))
def test_translation(seq, rc):
    seq = new_moltype.DNA.make_seq(seq=seq)
    if rc:
        seq = seq.rc()
    get_str = str(seq)
    assert get_str == "ATGTTT"
    aa = seq.get_translation()
    assert str(aa) == "MF"


def test_get_translation_include_stop():
    s = new_moltype.DNA.make_seq(seq="ATTTAACTT", name="s1")
    aa = s.get_translation(include_stop=True)
    assert str(aa) == "I*L"


def test_get_translation_trim_stop():
    s = new_moltype.DNA.make_seq(seq="ATTTCCTGA", name="s1")
    aa = s.get_translation(trim_stop=True)
    assert str(aa) == "IS"
    # no effect on internal stops
    s = new_moltype.DNA.make_seq(seq="ATTTAACTT", name="s1")
    aa = s.get_translation(include_stop=True, trim_stop=True)
    assert str(aa) == "I*L"


@pytest.mark.parametrize(
    "moltype, data",
    (
        ("dna", "ACGT"),
        ("rna", "ACGU"),
        ("protein", "ACDE"),
        ("protein_with_stop", "ACDE*"),
        ("bytes", "ACDE"),
    ),
)
def test_sequence_serialisation_round_trip(moltype, data):
    moltype = new_moltype.get_moltype(moltype)
    seq = moltype.make_seq(seq=data, name="seq1")

    rd = seq.to_rich_dict()
    got = deserialise_object(rd)
    assert isinstance(got, type(seq))
    assert got.to_rich_dict() == seq.to_rich_dict()
