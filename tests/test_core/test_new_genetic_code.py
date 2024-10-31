import pytest

from cogent3.core import new_alphabet, new_genetic_code


@pytest.mark.parametrize("gc", (1, "Standard", new_genetic_code.DEFAULT))
def test_get_code(gc):
    code = new_genetic_code.get_code(gc)
    assert code.name == "Standard"


def test_get_aa():
    code = new_genetic_code.DEFAULT
    assert code["ATG"] == "M"
    assert code["TAA"] == "*"


def test_stop_codons():
    code = new_genetic_code.DEFAULT
    assert code.stop_codons == {"TAA", "TAG", "TGA"}


@pytest.mark.parametrize("codons", ("TAA", "TAG", "TGA"))
def test_is_stop_true(codons):
    code = new_genetic_code.DEFAULT
    assert code.is_stop(codons)


def test_start_codons():
    code = new_genetic_code.DEFAULT
    dna = code.moltype.alphabet
    trinucs = code.codons
    # the starts are in the ncbi_start_codon_map
    expect = {
        dna.from_indices(trinucs.from_index(i))
        for i, e in enumerate(new_genetic_code.code_mapping[0][3])
        if e == "M"
    }
    assert code.start_codons == expect


def test_sense_codons():
    code = new_genetic_code.DEFAULT
    assert len(code.sense_codons) == 61


@pytest.mark.parametrize("invalid", ("AT", "TAAA", 23))
def test_invalid_index(invalid):
    with pytest.raises(new_genetic_code.InvalidCodonError):
        new_genetic_code.DEFAULT[invalid]


def test_translate():
    code = new_genetic_code.DEFAULT
    assert code.translate("ATGGGG") == code["ATG"] + code["GGG"]


def test_translate_rc():
    code = new_genetic_code.DEFAULT
    seq = "ATGGGG"
    rc = code.moltype.rc(seq)
    assert code.translate(rc, rc=True) == code["ATG"] + code["GGG"]


def test_translate_incomplete_ok():
    code = new_genetic_code.DEFAULT
    assert code.translate("ATGG--", incomplete_ok=True) == code["ATG"] + "-"


def test_translate_incomplete_not_ok():
    code = new_genetic_code.DEFAULT
    with pytest.raises(new_alphabet.AlphabetError):
        assert code.translate("ATGG--", incomplete_ok=False)


def test_to_table():
    code = new_genetic_code.DEFAULT
    table = code.to_table()
    assert table.shape[0] == 21


def test_available_codes():
    codes = new_genetic_code.available_codes()
    assert len(codes) == len(new_genetic_code.code_mapping)


def test_sizeframes():
    code = new_genetic_code.DEFAULT
    seq = "ATGGGGTAACAT"
    # following also validated against old and new genetic codes
    expect = {
        ("+", 0, "MG*H"),
        ("+", 1, "WGN"),
        ("+", 2, "GVT"),
        ("-", 0, "MLPH"),
        ("-", 1, "VTP"),
        ("-", 2, "CYP"),
    }
    got = set(code.sixframes(seq))
    assert got == expect


@pytest.mark.parametrize("code_id", range(1, 5))
def test_eq_true(code_id):
    gc0 = new_genetic_code.get_code(code_id)
    gc1 = new_genetic_code.get_code(code_id)
    assert gc0 == gc1


@pytest.mark.parametrize("code_id", range(1, 5))
def test_eq_false(code_id):
    gc0 = new_genetic_code.get_code(code_id)
    gc1 = new_genetic_code.get_code(code_id + 1)
    assert gc0 != gc1


@pytest.mark.parametrize(
    "seq, expect",
    [
        ("?GATCT", "XS"),
        ("GAT-T-", "D-"),
        ("GAT---", "D-"),
        ("GAT?CT", "DX"),
        ("GAT?C-", "DX"),
    ],
)
def test_translate_incomplete(seq, expect):
    gc = new_genetic_code.DEFAULT
    got = gc.translate(seq)
    assert got == expect


@pytest.mark.parametrize(
    "seq, expect",
    [
        ("AGATC?", "XS"),
        ("-A-ATC", "D-"),
        ("---ATC", "D-"),
        ("AG?ATC", "DX"),
        ("-G?ATC", "DX"),
    ],
)
def test_translate_incomplete_rc(seq, expect):
    gc = new_genetic_code.DEFAULT
    got = gc.translate(seq, rc=True)
    assert got == expect


def test_get_alphabet():
    gc = new_genetic_code.get_code(1)
    alpha = gc.get_alphabet()
    assert len(alpha) == 61
    assert alpha.to_index("TTT") == 0

    alpha_gap = gc.get_alphabet(include_gap=True)
    assert len(alpha_gap) == 62


def test_get_alphabet_with_stop():
    gc = new_genetic_code.get_code(1)
    alpha = gc.get_alphabet(include_stop=True)
    assert len(alpha) == 64

    alpha_gap = gc.get_alphabet(include_gap=True, include_stop=True)
    assert len(alpha_gap) == 65


@pytest.mark.parametrize("repr_method", ("__str__", "__repr__", "_repr_html_"))
def test_code_repr(repr_method):
    gc = new_genetic_code.get_code(1)
    got = getattr(gc, repr_method)()
    assert isinstance(got, str)


def test_to_regex():
    """creates a regex from aa seq to match a DNA sequence"""
    import re

    from cogent3 import make_seq

    dna = "ACCGAACAGGGC"
    aa = "TEQG"
    pattern = new_genetic_code.DEFAULT.to_regex(aa)
    assert "".join(re.findall(pattern, dna)) == dna
    # note that Z is Q or E
    aa = "TZQG"
    pattern = new_genetic_code.DEFAULT.to_regex(aa)
    assert "".join(re.findall(pattern, dna)) == dna
    aa = make_seq(aa, moltype="protein", new_type=True)
    pattern = new_genetic_code.DEFAULT.to_regex(aa)
    assert "".join(re.findall(pattern, dna)) == dna
