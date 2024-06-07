import pytest

from cogent3.core import new_genetic_code


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


def test_translate():
    code = new_genetic_code.DEFAULT
    assert code.translate("ATGGGG") == code["ATG"] + code["GGG"]


def test_translate_rc():
    code = new_genetic_code.DEFAULT
    seq = "ATGGGG"
    rc = code.moltype.rc(seq)
    assert code.translate(rc, rc=True) == code["ATG"] + code["GGG"]


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


@pytest.mark.parametrize("id", range(1, 5))
def test_compare_codes_true(id):
    gc = new_genetic_code.get_code(id)
    assert gc == gc


@pytest.mark.parametrize("id", range(1, 5))
def test_compare_codes_false(id):
    gc0 = new_genetic_code.get_code(id)
    gc1 = new_genetic_code.get_code(id + 1)
    assert gc0 != gc1
