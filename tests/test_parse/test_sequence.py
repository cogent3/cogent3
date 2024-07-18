# the delegator module
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
    yield request.param(path)


@pytest.fixture(params=(list, tuple))
def gde_data(DATA_DIR, request):
    path = DATA_DIR / "formattest.gde"
    yield request.param(path.read_text().splitlines())


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
