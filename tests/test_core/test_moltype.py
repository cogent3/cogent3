import numpy
import pytest

from cogent3.core import (
    alphabet as c3_alphabet,
)
from cogent3.core import (
    genetic_code as c3_genetic_code,
)
from cogent3.core import (
    moltype as c3_moltype,
)
from cogent3.core import (
    sequence as c3_sequence,
)


def test_is_compatible_alphabet():
    from cogent3.core.alphabet import CharAlphabet

    dna = c3_moltype.get_moltype("dna")
    alpha = CharAlphabet("TCAG")
    assert dna.is_compatible_alphabet(alpha)
    rna = c3_moltype.get_moltype("rna")
    assert not rna.is_compatible_alphabet(alpha)
    alpha = CharAlphabet("".join(dna.ambiguities))
    prot = c3_moltype.get_moltype("protein")
    assert not prot.is_compatible_alphabet(alpha)


@pytest.mark.parametrize(
    "name",
    ["dna", "rna", "protein", "protein_with_stop", "bytes", "text"],
)
def test_get_moltype(name):
    """correctly return a moltype by name"""
    mt = c3_moltype.get_moltype(name)
    assert mt.name == name
    mt = c3_moltype.get_moltype(name.upper())
    assert mt.name == name
    got = c3_moltype.get_moltype(mt)
    assert got is mt


def test_available_moltypes():
    t = c3_moltype.available_moltypes()
    assert t.shape[0] == 6
    got = str(t)
    assert "np." not in got
    assert "'dna'" in got


def test_str_moltype():
    dna = c3_moltype.get_moltype("dna")
    text = str(dna)
    assert isinstance(text, str)
    assert text == f"MolType({tuple('TCAG')})"


@pytest.mark.parametrize(
    ("seq", "data_type"),
    [
        ("ACCCG", str),
        (b"ACCCG", bytes),
        (numpy.array([2, 1, 1, 1, 3], dtype=numpy.uint8), numpy.ndarray),
    ],
)
@pytest.mark.parametrize("name", ["dna", "rna"])
def test_complement(name, seq, data_type):
    moltype = c3_moltype.get_moltype(name)
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


@pytest.mark.parametrize("data_type", [str, bytes, numpy.ndarray])
@pytest.mark.parametrize(
    "seq",
    [
        "N",
        "R",
        "Y",
        "?",  # IUPAC missing is also a degenerate
        "GCAUGUAGCUCGUCAGUCAGUACGUGCASCUAG",
        "ACGYAUGCUGYEWEWNFMNFUWBYBCWUYBCJWBEIWFUB",
    ],
)
def test_is_degenerate(seq, data_type):
    seq = make_typed(seq, data_type, c3_moltype.RNA)
    # note that the last sequence is NOT a valid RNA sequence, so
    # we need to turn off validation
    assert c3_moltype.RNA.is_degenerate(seq, validate=False)


@pytest.mark.parametrize("data_type", [str, bytes, numpy.ndarray])
@pytest.mark.parametrize(
    "seq",
    [
        "",
        "A",
        "UACGCUACAUGUACGUCAGUGCUAGCUA",
    ],
)
def test_is_not_degenerate(seq, data_type):
    seq = make_typed(seq, data_type, c3_moltype.RNA)
    assert not c3_moltype.RNA.is_degenerate(seq)


def test_is_degenerate_invalid():
    with pytest.raises(TypeError):
        c3_moltype.RNA.is_degenerate(list("GAG"))


@pytest.mark.parametrize("data_type", [str, bytes, numpy.ndarray])
@pytest.mark.parametrize(
    "seq",
    [
        "",
        "QWERTYUIOPASDFGHJKLZXCVBNM",
    ],
)
def test_text_moltype_is_not_degenerate(seq, data_type):
    """Text moltype should not be degenerate"""
    seq = make_typed(seq, data_type, c3_moltype.ASCII)
    assert not c3_moltype.ASCII.is_degenerate(seq)


@pytest.mark.parametrize("data_type", [str, bytes, numpy.ndarray])
@pytest.mark.parametrize(
    "seq",
    [
        "-",
        "Y-",
        "GC--A",
        "-ACGYA",
    ],
)
def test_is_gapped(seq, data_type):
    seq = make_typed(seq, data_type, c3_moltype.RNA)
    assert c3_moltype.RNA.is_gapped(seq)


@pytest.mark.parametrize("data_type", [str, bytes, numpy.ndarray])
@pytest.mark.parametrize(
    "seq",
    [
        "",
        "Y",
        "GCA",
        "ACGYA",
    ],
)
def test_not_is_gapped(seq, data_type):
    seq = make_typed(seq, data_type, c3_moltype.RNA)
    assert not c3_moltype.RNA.is_gapped(seq)


@pytest.mark.parametrize("moltype", [c3_moltype.DNA, c3_moltype.RNA])
def test_gap_index_constant(moltype):
    # make sure gap index is always the same
    assert moltype.gapped_alphabet.gap_index == moltype.degen_gapped_alphabet.gap_index


@pytest.mark.parametrize("data_type", [str, bytes, numpy.ndarray])
@pytest.mark.parametrize("moltype", [c3_moltype.DNA, c3_moltype.RNA])
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
    got = c3_moltype.DNA.resolve_ambiguity("AT?", allow_gap=False)
    assert len(got) == 4
    assert len(got[0]) == 3


def test_resolve_ambiguity_codons():
    gc = c3_genetic_code.get_code(1)
    codon_alpha = gc.get_alphabet(include_stop=False)
    codon_alpha_w_gap = codon_alpha.with_gap_motif()
    assert len(c3_moltype.DNA.resolve_ambiguity("AT?", alphabet=codon_alpha_w_gap)) == 4
    assert (
        len(c3_moltype.DNA.resolve_ambiguity("???", alphabet=codon_alpha_w_gap)) == 62
    )
    assert len(c3_moltype.DNA.resolve_ambiguity("---", alphabet=codon_alpha_w_gap)) == 1

    assert len(c3_moltype.DNA.resolve_ambiguity("AT?", alphabet=codon_alpha)) == 4
    assert len(c3_moltype.DNA.resolve_ambiguity("???", alphabet=codon_alpha)) == 61

    with pytest.raises(c3_alphabet.AlphabetError):
        c3_moltype.DNA.resolve_ambiguity("at-")
    with pytest.raises(c3_alphabet.AlphabetError):
        c3_moltype.DNA.resolve_ambiguity("---", alphabet=codon_alpha)


@pytest.mark.parametrize("char", ["N", "R", "Y", "W", "S", "M", "?"])
def test_is_ambiguity_true(char):
    assert c3_moltype.DNA.is_ambiguity(char)


@pytest.mark.parametrize("char", ["-", "A", "T", "C", "G"])
def test_is_ambiguity_false(char):
    assert not c3_moltype.DNA.is_ambiguity(char)


@pytest.mark.parametrize("data_type", [str, bytes, numpy.ndarray])
@pytest.mark.parametrize("moltype", ["text", "rna"])
@pytest.mark.parametrize(
    ("seq", "expect"),
    [
        ("", ""),
        ("GAUGgaug", "GAUGgaug"),
        ("----", ""),
        ("--GAUG--", "GAUG"),
        ("?gaug-", "gaug"),
        ("-gaug", "gaug"),
        ("gaug-", "gaug"),
        ("-g?a-u?g-", "gaug"),
    ],
)
def test_degap(seq, expect, data_type, moltype):
    """MolType degap should remove all gaps from sequence"""
    moltype = c3_moltype.get_moltype(moltype)
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
        c3_moltype.PROTEIN.strand_symmetric_motifs()

    got = c3_moltype.DNA.strand_symmetric_motifs(motif_length=1)
    expect = {("A", "T"), ("C", "G")}
    assert got == expect
    got = c3_moltype.RNA.strand_symmetric_motifs(motif_length=1)
    expect = {("A", "U"), ("C", "G")}
    assert got == expect
    got = c3_moltype.DNA.strand_symmetric_motifs(motif_length=2)
    assert len(got) == 8
    got = c3_moltype.DNA.strand_symmetric_motifs(motif_length=3)
    assert len(got) == 32


@pytest.mark.parametrize(
    "moltype",
    [
        c3_moltype.ASCII,
        c3_moltype.DNA,
        c3_moltype.RNA,
        c3_moltype.PROTEIN,
        c3_moltype.PROTEIN_WITH_STOP,
    ],
)
def test_gaps(moltype):
    got = moltype.gaps
    expect = frozenset({"-", "?"})
    assert got == expect


def test_gaps_bytes():
    got = c3_moltype.get_moltype("bytes").gaps
    expect = frozenset()
    assert got == expect


def test_gaps_none():
    mt = c3_moltype.MolType(
        "no_gap",
        monomers="".join(c3_moltype.IUPAC_DNA_chars),
        make_seq=c3_sequence.DnaSequence,
        gap=None,
    )

    got = mt.gaps
    expect = frozenset({"?"})
    assert got == expect

    mt = c3_moltype.MolType(
        "no_missing",
        monomers="".join(c3_moltype.IUPAC_DNA_chars),
        make_seq=c3_sequence.DnaSequence,
        missing=None,
    )
    got = mt.gaps
    expect = frozenset({"-"})
    assert got == expect


@pytest.mark.parametrize(
    "label",
    [
        "bytes",
        "dna",
        "rna",
        "protein",
        "protein_with_stop",
        "text",
    ],
)
def test_most_degenerate_alphabet(label):
    moltype = c3_moltype.get_moltype(label)
    got = moltype.most_degen_alphabet()
    # expected value for number of characters is
    # length of monomers + len(moltype gaps) + len(ambiguities)
    ambigs = set(moltype.ambiguities or ())
    ambigs |= set(moltype.missing) if moltype.missing else set()
    num_ambigs = len(ambigs)
    expected = len(moltype.alphabet) + len(moltype.gap or "") + num_ambigs
    assert len(got) == expected


def test_validate_seq():
    moltype = c3_moltype.ASCII
    alpha = moltype.most_degen_alphabet()
    assert alpha.is_valid("?gau")


def test_degap_raises():
    with pytest.raises(c3_alphabet.AlphabetError):
        c3_moltype.RNA.degap("-g?a-u?g-", validate=True)


def test_is_invalid_rna():
    rna = c3_moltype.RNA
    seq = "ACGYAUGCUGYEWEWNFMNFUWBYBCWUYBCJWBEIWFUB"
    assert not rna.is_valid(seq)


@pytest.mark.parametrize("data_type", [str, bytes, numpy.ndarray])
def test_strip_degenerate(data_type):
    seq = make_typed("N-TCGA-N-TCNGA", data_type, c3_moltype.DNA)
    got = c3_moltype.DNA.strip_degenerate(seq)
    expect = make_typed("-TCGA--TCGA", data_type, c3_moltype.DNA)
    assert (
        numpy.array_equal(got, expect) if data_type is numpy.ndarray else got == expect
    )

    seq = make_typed("AST-B-ASZT-", data_type, c3_moltype.PROTEIN)
    got = c3_moltype.PROTEIN.strip_degenerate(seq)
    expect = make_typed("AST--AST-", data_type, c3_moltype.PROTEIN)
    assert (
        numpy.array_equal(got, expect) if data_type is numpy.ndarray else got == expect
    )

    seq = make_typed("AST-B*-ASZT-", data_type, c3_moltype.PROTEIN_WITH_STOP)
    got = c3_moltype.PROTEIN_WITH_STOP.strip_degenerate(seq)
    expect = make_typed("AST-*-AST-", data_type, c3_moltype.PROTEIN_WITH_STOP)
    assert (
        numpy.array_equal(got, expect) if data_type is numpy.ndarray else got == expect
    )


@pytest.mark.parametrize("data_type", [str, bytes])
def test_strip_bad(data_type):
    """removes characters that are not in the alphabet"""

    # DNA - Q not in DNA alphabet/ambiguities
    seq = make_typed("TCGA-Q?NRYWSKMBDHV", data_type, c3_moltype.DNA)
    got = c3_moltype.DNA.strip_bad(seq)
    expect = make_typed("TCGA-?NRYWSKMBDHV", data_type, c3_moltype.DNA)
    assert got == expect

    # Protein - J and * not in Protein alphabet/ambiguities
    seq = make_typed("ASTBJZ*-X", data_type, c3_moltype.PROTEIN)
    got = c3_moltype.PROTEIN.strip_bad(seq)
    expect = make_typed("ASTBZ-X", data_type, c3_moltype.PROTEIN)
    assert got == expect

    # Protein with stop - J  not in Protein alphabet/ambiguities
    seq = make_typed("ASTBJZ*-X", data_type, c3_moltype.PROTEIN_WITH_STOP)
    got = c3_moltype.PROTEIN_WITH_STOP.strip_bad(seq)
    expect = make_typed("ASTBZ*-X", data_type, c3_moltype.PROTEIN_WITH_STOP)
    assert got == expect


@pytest.mark.parametrize("data_type", [str, bytes])
def test_strip_bad_and_gaps(data_type):
    """removes characters that are not in the alphabet and gaps"""

    # DNA - Q not in DNA alphabet/ambiguities
    seq = make_typed("TCGA-QNRYWSKMBDHV", data_type, c3_moltype.DNA)
    got = c3_moltype.DNA.strip_bad_and_gaps(seq)
    expect = make_typed("TCGANRYWSKMBDHV", data_type, c3_moltype.DNA)
    assert got == expect

    # Protein - J and * not in Protein alphabet/ambiguities
    seq = make_typed("ASTBJZ*X", data_type, c3_moltype.PROTEIN)
    got = c3_moltype.PROTEIN.strip_bad(seq)
    expect = make_typed("ASTBZX", data_type, c3_moltype.PROTEIN)
    assert got == expect

    # Protein with stop - J  not in Protein alphabet/ambiguities
    seq = make_typed("ASTBJZ*X", data_type, c3_moltype.PROTEIN_WITH_STOP)
    got = c3_moltype.PROTEIN_WITH_STOP.strip_bad(seq)
    expect = make_typed("ASTBZ*X", data_type, c3_moltype.PROTEIN_WITH_STOP)
    assert got == expect


@pytest.mark.parametrize("data_type", [str, bytes, numpy.ndarray])
def test_disambiguate_empty(data_type):
    d = c3_moltype.RNA.disambiguate
    seq = make_typed("", data_type, c3_moltype.RNA)
    assert (
        numpy.array_equal(d(seq), seq) if data_type is numpy.ndarray else d(seq) == seq
    )


@pytest.mark.parametrize("data_type", [str, bytes, numpy.ndarray])
def test_disambiguate_strip(data_type):
    """MolType disambiguate should remove degenerate bases"""
    d = c3_moltype.RNA.disambiguate
    seq = make_typed("AUN-YRS-WKMCGWMRNMWRKY", data_type, c3_moltype.RNA)
    got = d(seq, "strip")
    expect = make_typed("AU--CG", data_type, c3_moltype.RNA)
    assert (
        numpy.array_equal(got, expect) if data_type is numpy.ndarray else got == expect
    )


def test_disambiguate_random_str():
    # TODO: add tests for bytes and numpy.array
    d = c3_moltype.RNA.disambiguate
    s = "AUN-YRS-WKMCGWMRNMWRKY"
    t = d(s, "random")

    assert len(s) == len(t)

    for i, j in zip(s, t, strict=False):
        if i in c3_moltype.RNA.ambiguities:
            assert j in c3_moltype.RNA.ambiguities[i]

    assert not c3_moltype.RNA.is_degenerate(t)
    with pytest.raises(NotImplementedError):
        d(s, method="xyz")


@pytest.mark.parametrize("data_type", [str, bytes])
@pytest.mark.parametrize(
    ("data", "expect"),
    [("---CGAUGCAU---ACGHC---ACGUCAGU---", 12), ("", 0), ("ACGU", 0), ("-", 1)],
)
def test_count_gaps(data, expect, data_type):
    seq = make_typed(data, data_type, c3_moltype.RNA)
    got = c3_moltype.RNA.count_gaps(seq)
    assert got == expect


def test_count_degenerate():
    """MolType count_degenerate should return correct degen base count"""
    d = c3_moltype.RNA.count_degenerate
    assert d("") == 0
    assert d("GACUGCAUGCAUCGUACGUCAGUACCGA") == 0
    assert d("N") == 1
    assert d("NRY") == 3
    assert d("ACGUAVCUAGCAUNUCAGUCAGyUACGUCAGS") == 4


def test_count_variants():
    """MolType count_variants should return correct # possible sequences"""
    p = c3_moltype.RNA.count_variants
    assert p("") == 1
    assert p("ACGUGCAUCAGUCGUGCAU") == 1
    assert p("N") == 4
    assert p("R") == 2
    assert p("H") == 3
    assert p("NRH") == 24
    assert p("AUGCNGUCAG-AURGAUC--GAUHCGAUACGWS") == 96


def test_degenerate_from_seq():
    """RnaMolType degenerate_from_seq should give correct results"""
    rna = c3_moltype.get_moltype("rna")
    d = rna.degenerate_from_seq
    # check monomers
    assert d("A") == "A"
    assert d("C") == "C"
    # check seq of monomers
    assert d("AAA") == "A"
    # check some 2- to 4-way cases
    assert d("AG") == "R"
    assert d("CG") == "S"
    assert d("ACG") == "V"
    assert d("AGCU") == "N"
    # check some cases with gaps
    assert d("A---ACU") == "?"
    assert d("AAA-") == "?"
    assert d("---") == "-"
    # check mixed case example
    assert d("AAAAAA") == "A"
    # check example with degenerate symbols in set
    assert d("RS") == "V"
    assert d("RN-") == "?"
    # check that it works for proteins as well
    protein = c3_moltype.get_moltype("protein")
    p = protein.degenerate_from_seq
    assert p("A") == "A"
    assert p("AAA") == "A"
    assert p("DN") == "B"
    assert p("---") == "-"
    assert p("ACD") == "X"
    assert p("ABD") == "X"
    assert p("ACX") == "X"
    assert p("AC-") == "?"


def test_degenerate_from_seq_general():
    dna = c3_moltype.get_moltype("dna")
    assert dna.degenerate_from_seq("ACGT--") == "?"


@pytest.mark.parametrize(
    ("moltype", "zero", "second", "last"),
    [
        ("DNA", "T", "A", "G"),
        ("RNA", "U", "A", "G"),
        ("PROTEIN", "A", "D", "Y"),
        ("PROTEIN_WITH_STOP", "A", "D", "*"),
    ],
)
def test_char_order(moltype, zero, second, last):
    prefix = "IUPAC_" if "STOP" not in moltype else ""
    attr = getattr(c3_moltype, f"{prefix}{moltype}_chars")
    assert attr[0] == zero
    assert attr[2] == second
    assert attr[-1] == last


@pytest.mark.parametrize(
    "seq",
    [
        mt.make_seq(seq="ACGG")
        for mt in (
            c3_moltype.DNA,
            c3_moltype.RNA,
            c3_moltype.PROTEIN,
            c3_moltype.ASCII,
        )
    ],
)
def test_make_seq_with_seq(seq):
    new_seq = c3_moltype.DNA.make_seq(seq=seq)
    assert new_seq.moltype is c3_moltype.DNA
    assert (new_seq is seq) if seq.moltype.name == "dna" else seq is not new_seq


def test_make_seq_with_seq_invalid_moltype():
    seq = c3_moltype.BYTES.make_seq(seq=b"".join([bytes([i]) for i in range(4)]))
    with pytest.raises(c3_moltype.MolTypeError):
        _ = c3_moltype.DNA.make_seq(seq=seq)


def test_to_regex():
    """returns a valid regex"""
    seq = "ACYGR"
    regular_expression = c3_moltype.DNA.to_regex(seq=seq)
    assert regular_expression == "AC[CT]G[AG]"
    # raises an exception if a string is already a regex, or invalid
    with pytest.raises(ValueError):
        c3_moltype.DNA.to_regex("(?:GAT|GAC)(?:GGT|GGC|GGA|GGG)(?:GAT|GAC)(?:CAA|CAG)")


@pytest.mark.parametrize("moltype", ["dna", "rna"])
def test_is_nucleic(moltype):
    mt = c3_moltype.get_moltype(moltype)
    assert mt.is_nucleic


@pytest.mark.parametrize("moltype", ["protein", "text", "bytes", "protein_with_stop"])
def test_not_is_nucleic(moltype):
    mt = c3_moltype.get_moltype(moltype)
    assert not mt.is_nucleic


def test_moltype_coerce_seqs():
    dna = c3_moltype.get_moltype("dna")
    rna_seq = "AUUG"
    dna_seq = "ATTG"
    rna = c3_moltype.get_moltype("rna")
    assert str(dna.make_seq(seq=rna_seq)) == dna_seq
    assert str(rna.make_seq(seq=dna_seq)) == rna_seq

    # no coercion in protein seq
    prot = c3_moltype.get_moltype("protein")
    assert str(prot.make_seq(seq=rna_seq)) == rna_seq


@pytest.mark.parametrize(
    "moltype",
    ["protein", "text", "bytes", "protein_with_stop", "dna", "rna"],
)
def test_json_roundtrip_moltype(moltype):
    from cogent3.util.deserialise import deserialise_object

    mt = c3_moltype.get_moltype(moltype)
    got = deserialise_object(mt.to_json())
    # because our moltype instances are singletons,
    # we expect the same object
    assert got is mt


@pytest.mark.parametrize(
    ("moltype", "seq", "expect"),
    [
        (c3_moltype.DNA, "ACGT", False),  # No ambiguity
        (c3_moltype.DNA, "ACGTN", True),  # N is ambiguous
        (c3_moltype.DNA, "ACGTY", True),  # Y is ambiguous
        (c3_moltype.RNA, "ACGU", False),  # No ambiguity
        (c3_moltype.RNA, "ACGUN", True),  # N is ambiguous
        (c3_moltype.RNA, "ACGUY", True),  # Y is ambiguous
        (c3_moltype.PROTEIN, "ACDEFGHIKLMNPQRSTVWY", False),  # No ambiguity
        (c3_moltype.PROTEIN, "ACDEFGHIKLMNPQRSTVWYX", True),  # X is ambiguous
        (
            c3_moltype.PROTEIN_WITH_STOP,
            "ACDEFGHIKLM*NPQRSTVWY",
            False,
        ),  # * is not ambiguous (stop)
    ],
)
@pytest.mark.parametrize("cast", [str, bytes, numpy.ndarray])
def test_has_ambiguity(moltype, seq, expect, cast):
    seq = make_typed(seq, cast, moltype)
    assert moltype.has_ambiguity(seq) == expect


@pytest.mark.parametrize(
    "seq",
    ["ACGT", "ACGTN"],
)
@pytest.mark.parametrize("moltype", ["text", "bytes"])
@pytest.mark.parametrize("cast", [str, bytes, numpy.ndarray])
def test_has_ambiguity_all_false(moltype, seq, cast):
    moltype = c3_moltype.get_moltype(moltype)
    # these moltypes don't have an ambuity code
    seq = make_typed(seq, cast, moltype)
    assert not moltype.has_ambiguity(seq)


def test_has_ambiguity_validation():
    """raises TypeError if input not a string"""
    with pytest.raises(TypeError):
        c3_moltype.DNA.has_ambiguity(None)

    with pytest.raises(TypeError):
        c3_moltype.DNA.has_ambiguity([])

    with pytest.raises(TypeError):
        c3_moltype.DNA.has_ambiguity(1)


def test_custom_moltype():
    mt = c3_moltype.MolType(
        name="dna-gapped",
        make_seq=c3_sequence.DnaSequence,
        monomers="".join(c3_moltype.IUPAC_DNA_chars),
        ambiguities=c3_moltype.IUPAC_DNA_ambiguities,
        complements=c3_moltype.IUPAC_DNA_ambiguities_complements,
        pairing_rules=c3_moltype.DNA_STANDARD_PAIRS,
        gap=".",
    )
    assert "." in mt.gaps
