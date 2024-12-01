import pathlib

import pytest

from cogent3.parse import sequence


def test_line_based_wrap_invalid():
    parser = sequence.LineBasedParser(lambda x: x)
    with pytest.raises(TypeError):
        parser({"a": "GAT"})


def test_line_based_wrap_invalid_path():
    parser = sequence.LineBasedParser(lambda x: x)
    path = pathlib.Path("0-0-swedfghjkl")
    with pytest.raises(FileNotFoundError):
        list(parser(path))


def test_line_based_wrap_invalid_types():
    parser = sequence.LineBasedParser(lambda x: x)
    path = pathlib.Path("0-0-swedfghjkl")
    with pytest.raises(FileNotFoundError):
        list(parser(path))


@pytest.fixture(params=(str, pathlib.Path))
def gde_path(DATA_DIR, request):
    path = DATA_DIR / "formattest.gde"
    return request.param(path)


@pytest.fixture(params=(list, tuple))
def gde_data(DATA_DIR, request):
    path = DATA_DIR / "formattest.gde"
    return request.param(path.read_text().splitlines())


def test_line_based_wrap_load_valid_path_types(gde_path):
    parser = sequence.get_parser("gde")
    got = dict(parser(gde_path))
    assert len(got) == 10


def test_line_based_wrap_load_valid_input_types(gde_data):
    parser = sequence.get_parser("gde")
    got = dict(parser(gde_data))
    assert len(got) == 10


def test_get_unknown_format():
    with pytest.raises(ValueError):
        sequence.get_parser("nope")


@pytest.mark.internet
def test_line_based_url(DATA_DIR):
    fname = "long_testseqs.fasta"
    url = f"https://raw.githubusercontent.com/cogent3/cogent3/develop/doc/data/{fname}"

    parser = sequence.LineBasedParser(lambda x: x)
    from_url = list(parser(url))
    from_local = list(parser(DATA_DIR / fname))
    assert from_url == from_local


@pytest.mark.parametrize("fmt", ("gb", "gbk", "gbff", "genbank"))
def test_is_genbank(fmt):
    assert sequence.is_genbank(fmt)


@pytest.mark.parametrize("fmt", ("blah", "fa", "xml", "nex", None))
def test_is_not_genbank(fmt):
    assert not sequence.is_genbank(fmt)


@pytest.fixture(params=("phylip", "phy"))
def phylip_file(DATA_DIR, tmp_path, request):
    with open(DATA_DIR / "interleaved.phylip") as f:
        data = f.read()
    outpath = tmp_path / f"interleaved.{request.param}"
    outpath.write_text(data)
    return outpath


def test_select_parser_phylip_suffixes(phylip_file):
    from cogent3 import load_aligned_seqs

    got = load_aligned_seqs(phylip_file)
    assert set(got.names) == {"human", "chimp", "mouse"}
