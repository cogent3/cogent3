import numpy
import pytest

from cogent3.core import moltype, new_alphabet, new_moltype, new_sequence


def test_make_pairs():
    orig = moltype.make_pairs(
        pairs=moltype.RnaStandardPairs,
        monomers=moltype.IUPAC_RNA_chars,
        gaps=moltype.IUPAC_gap,
        degenerates=moltype.IUPAC_RNA_ambiguities,
    )
    new = new_moltype.make_pairs(
        pairs=new_moltype.RNA_STANDARD_PAIRS,
        monomers=new_moltype.IUPAC_RNA_chars,
        gaps=new_moltype.IUPAC_gap,
        degenerates=new_moltype.IUPAC_RNA_ambiguities,
    )
    # convert the old object to a frozen set keyed one for comparison
    orig = {frozenset(k): v for k, v in orig.items()}
    assert new == orig


def test_is_compatible_alphabet():
    from cogent3.core.new_alphabet import CharAlphabet

    dna = new_moltype.get_moltype("dna")
    alpha = CharAlphabet("TCAG")
    assert dna.is_compatible_alphabet(alpha)
    rna = new_moltype.get_moltype("rna")
    assert not rna.is_compatible_alphabet(alpha)
    alpha = CharAlphabet("".join(dna.ambiguities))
    prot = new_moltype.get_moltype("protein")
    assert not prot.is_compatible_alphabet(alpha)


def test_is_compatible_alphabet_strict():
    from cogent3.core.alphabet import CharAlphabet

    dna = new_moltype.get_moltype("dna")
    alpha1 = CharAlphabet("TCAG")
    assert dna.is_compatible_alphabet(alpha1, strict=True)
    # returns False if the order is not exactly the same
    alpha1 = CharAlphabet("CTAG")
    assert not dna.is_compatible_alphabet(alpha1, strict=True)


@pytest.mark.parametrize(
    "name", ("dna", "rna", "protein", "protein_with_stop", "bytes", "text")
)
def test_get_moltype(name):
    """correctly return a moltype by name"""
    mt = new_moltype.get_moltype(name)
    assert mt.name == name
    mt = new_moltype.get_moltype(name.upper())
    assert mt.name == name
    got = new_moltype.get_moltype(mt)
    assert got is mt


def test_available_moltypes():
    t = new_moltype.available_moltypes()
    assert t.shape[0] == 6


def test_str_moltype():
    dna = new_moltype.get_moltype("dna")
    text = str(dna)
    assert isinstance(text, str)
    assert text == f"MolType({tuple('TCAG')})"


@pytest.mark.parametrize(
    "seq, data_type",
    (
        ("ACCCG", str),
        (b"ACCCG", bytes),
        (numpy.array([2, 1, 1, 1, 3], dtype=numpy.uint8), numpy.ndarray),
    ),
)
@pytest.mark.parametrize("name", ("dna", "rna"))
def test_complement(name, seq, data_type):
    moltype = new_moltype.get_moltype(name)
    expect = "TGGGC" if name == "dna" else "UGGGC"
    got = moltype.complement(seq)
    assert (
        got == make_typed(expect, data_type, moltype)
        if data_type is not numpy.ndarray
        else numpy.array_equal(got, make_typed(expect, data_type, moltype))
    )


def make_typed(seq, data_type, moltype):
    if data_type is numpy.ndarray:
        seq = moltype.most_degen_alphabet().to_indices(seq)
    elif data_type is bytes:
        seq = seq.encode("utf-8")
    return seq


@pytest.mark.parametrize("data_type", (str, bytes, numpy.ndarray))
@pytest.mark.parametrize(
    "seq",
    (
        "N",
        "R",
        "Y",
        "?",  # IUPAC missing is also a degenerate
        "GCAUGUAGCUCGUCAGUCAGUACGUGCASCUAG",
        "ACGYAUGCUGYEWEWNFMNFUWBYBCWUYBCJWBEIWFUB",
    ),
)
def test_is_degenerate(seq, data_type):
    seq = make_typed(seq, data_type, new_moltype.RNA)
    # note that the last sequence is NOT a valid RNA sequence, so
    # we need to turn off validation
    assert new_moltype.RNA.is_degenerate(seq, validate=False)


@pytest.mark.parametrize("data_type", (str, bytes, numpy.ndarray))
@pytest.mark.parametrize(
    "seq",
    (
        "",
        "A",
        "UACGCUACAUGUACGUCAGUGCUAGCUA",
    ),
)
def test_is_not_degenerate(seq, data_type):
    seq = make_typed(seq, data_type, new_moltype.RNA)
    assert not new_moltype.RNA.is_degenerate(seq)


def test_is_degenerate_invalid():
    with pytest.raises(TypeError):
        new_moltype.RNA.is_degenerate(list("GAG"))


@pytest.mark.parametrize("data_type", (str, bytes, numpy.ndarray))
@pytest.mark.parametrize(
    "seq",
    (
        "-",
        "Y-",
        "GC--A",
        "-ACGYA",
    ),
)
def test_is_gapped(seq, data_type):
    seq = make_typed(seq, data_type, new_moltype.RNA)
    assert new_moltype.RNA.is_gapped(seq)


@pytest.mark.parametrize("data_type", (str, bytes, numpy.ndarray))
@pytest.mark.parametrize(
    "seq",
    (
        "",
        "Y",
        "GCA",
        "ACGYA",
    ),
)
def test_not_is_gapped(seq, data_type):
    seq = make_typed(seq, data_type, new_moltype.RNA)
    assert not new_moltype.RNA.is_gapped(seq)


@pytest.mark.parametrize("moltype", (new_moltype.DNA, new_moltype.RNA))
def test_gap_index_constant(moltype):
    # make sure gap index is always the same
    assert moltype.gapped_alphabet.gap_index == moltype.degen_gapped_alphabet.gap_index


@pytest.mark.parametrize("data_type", (str, bytes, numpy.ndarray))
@pytest.mark.parametrize("moltype", (new_moltype.DNA, new_moltype.RNA))
def test_get_degenerate_positions(data_type, moltype):
    seq = make_typed("ASA", data_type, moltype)
    got = moltype.get_degenerate_positions(seq)
    expect = [1]
    assert got == expect

    seq = make_typed("A-SA", data_type, moltype)
    got = moltype.get_degenerate_positions(seq)
    expect = [1, 2]
    assert got == expect

    got = moltype.get_degenerate_positions(seq, include_gap=False)
    expect = [2]
    assert got == expect

    seq = make_typed("BAB", data_type, moltype)
    got = moltype.get_degenerate_positions(seq)
    expect = [0, 2]
    assert got == expect

    seq = make_typed("---", data_type, moltype)
    got = moltype.get_degenerate_positions(seq)
    expect = [0, 1, 2]
    assert got == expect

    seq = make_typed("", data_type, moltype)
    got = moltype.get_degenerate_positions(seq)
    expect = []
    assert got == expect


def test_resolve_ambiguity_nucs():
    got = new_moltype.DNA.resolve_ambiguity("AT?", allow_gap=False)
    assert len(got) == 4
    assert len(got[0]) == 3


def test_resolve_ambiguity_codons():
    # refactor to using new_genetic_code when codon alphabet class is implemented
    from cogent3 import get_code

    gc = get_code(1)
    codon_alpha = gc.get_alphabet(include_stop=False)
    codon_alpha_w_gap = codon_alpha.with_gap_motif()
    assert (
        len(new_moltype.DNA.resolve_ambiguity("AT?", alphabet=codon_alpha_w_gap)) == 4
    )
    assert (
        len(new_moltype.DNA.resolve_ambiguity("???", alphabet=codon_alpha_w_gap)) == 62
    )
    assert (
        len(new_moltype.DNA.resolve_ambiguity("---", alphabet=codon_alpha_w_gap)) == 1
    )

    assert len(new_moltype.DNA.resolve_ambiguity("AT?", alphabet=codon_alpha)) == 4
    assert len(new_moltype.DNA.resolve_ambiguity("???", alphabet=codon_alpha)) == 61

    with pytest.raises(new_alphabet.AlphabetError):
        new_moltype.DNA.resolve_ambiguity("at-")
    with pytest.raises(new_alphabet.AlphabetError):
        new_moltype.DNA.resolve_ambiguity("---", alphabet=codon_alpha)


@pytest.mark.parametrize("char", ("N", "R", "Y", "W", "S", "M", "?"))
def test_is_ambiguity_true(char):
    assert new_moltype.DNA.is_ambiguity(char)


@pytest.mark.parametrize("char", ("-", "A", "T", "C", "G"))
def test_is_ambiguity_false(char):
    assert not new_moltype.DNA.is_ambiguity(char)


@pytest.mark.parametrize("data_type", (str, bytes, numpy.ndarray))
@pytest.mark.parametrize("moltype", ("text", "rna"))
@pytest.mark.parametrize(
    "seq, expect",
    (
        ("", ""),
        ("GAUGgaug", "GAUGgaug"),
        ("----", ""),
        ("--GAUG--", "GAUG"),
        ("?gaug-", "gaug"),
        ("-gaug", "gaug"),
        ("gaug-", "gaug"),
        ("-g?a-u?g-", "gaug"),
    ),
)
def test_degap(seq, expect, data_type, moltype):
    """MolType degap should remove all gaps from sequence"""
    moltype = new_moltype.get_moltype(moltype)
    seq = seq.upper() if moltype.name == "rna" else seq
    expect = expect.upper() if moltype.name == "rna" else expect
    seq = make_typed(seq, data_type, moltype)
    expect = make_typed(expect, data_type, moltype)
    got = moltype.degap(seq)
    assert (
        numpy.array_equal(got, expect) if data_type == numpy.ndarray else got == expect
    )


def test_strand_symmetric_motifs():
    """construction of strand symmetric motif sets"""
    # fails for a moltype with no strand complement
    with pytest.raises(TypeError):
        new_moltype.PROTEIN.strand_symmetric_motifs()

    got = new_moltype.DNA.strand_symmetric_motifs(motif_length=1)
    expect = set([("A", "T"), ("C", "G")])
    assert got == expect
    got = new_moltype.RNA.strand_symmetric_motifs(motif_length=1)
    expect = set([("A", "U"), ("C", "G")])
    assert got == expect
    got = new_moltype.DNA.strand_symmetric_motifs(motif_length=2)
    assert len(got) == 8
    got = new_moltype.DNA.strand_symmetric_motifs(motif_length=3)
    assert len(got) == 32


@pytest.mark.parametrize(
    "moltype",
    (
        new_moltype.ASCII,
        new_moltype.DNA,
        new_moltype.RNA,
        new_moltype.PROTEIN,
        new_moltype.PROTEIN_WITH_STOP,
    ),
)
def test_gaps(moltype):
    got = moltype.gaps
    expect = frozenset({"-", "?"})
    assert got == expect


def test_gaps_bytes():
    got = new_moltype.get_moltype("bytes").gaps
    expect = frozenset()
    assert got == expect


def test_gaps_none():
    mt = new_moltype.MolType(
        "no_gap",
        monomers="".join(new_moltype.IUPAC_DNA_chars),
        make_seq=new_sequence.DnaSequence,
        gap=None,
    )

    got = mt.gaps
    expect = frozenset({"?"})
    assert got == expect

    mt = new_moltype.MolType(
        "no_missing",
        monomers="".join(new_moltype.IUPAC_DNA_chars),
        make_seq=new_sequence.DnaSequence,
        missing=None,
    )
    got = mt.gaps
    expect = frozenset({"-"})
    assert got == expect


@pytest.mark.parametrize(
    "label",
    (
        "bytes",
        "dna",
        "rna",
        "protein",
        "protein_with_stop",
        "text",
    ),
)
def test_most_degenerate_alphabet(label):
    moltype = new_moltype.get_moltype(label)
    got = moltype.most_degen_alphabet()
    # expected value for number of characters is
    # length of monomers + len(moltype gaps) + len(ambiguities)
    num_ambigs = len(moltype.ambiguities or [])
    expected = len(moltype.alphabet) + len(moltype.gaps or []) + num_ambigs
    assert len(got) == expected


def test_validate_seq():
    moltype = new_moltype.ASCII
    alpha = moltype.most_degen_alphabet()
    assert alpha.is_valid("?gau")


def test_degap_raises():
    with pytest.raises(new_alphabet.AlphabetError):
        new_moltype.RNA.degap("-g?a-u?g-", validate=True)


def test_is_invalid_rna():
    rna = new_moltype.RNA
    seq = "ACGYAUGCUGYEWEWNFMNFUWBYBCWUYBCJWBEIWFUB"
    assert not rna.is_valid(seq)
