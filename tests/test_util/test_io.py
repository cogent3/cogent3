import pathlib

import pytest
from scinexus.io_util import open_

from cogent3.util.io import (
    _path_relative_to_zip_parent,
    iter_line_blocks,
    iter_splitlines,
)


@pytest.fixture
def home_file(DATA_DIR, HOME_TMP_DIR) -> str:
    """makes a temporary directory with file"""
    fn = "sample.tsv"
    contents = (DATA_DIR / fn).read_text()
    (HOME_TMP_DIR / fn).expanduser().write_text(contents)
    return str(HOME_TMP_DIR / fn)


def test_path_relative_to_zip_parent():
    """correctly generates member paths for a zip archive"""
    zip_path = pathlib.Path("some/path/to/a/data.zip")
    for member in ("data/member.txt", "member.txt", "a/b/c/member.txt"):
        got = _path_relative_to_zip_parent(zip_path, pathlib.Path(member))
        assert got.parts[0] == "data"


def test_iter_splitlines_one(tmp_path):
    # file has a single line
    path = tmp_path / "one-line.txt"
    value = "We have text on one line."
    path.write_text(value)
    got = list(iter_splitlines(path))
    assert got == [value]


@pytest.mark.parametrize("newline", ["\n", "\r\n"])
def test_iter_splitlines_line_diff_newline(tmp_path, newline):
    path = tmp_path / "multi-line.txt"
    value = ["We have some", "text on different lines", "which load"]
    with open_(path, mode="w", newline=newline) as out:
        out.write("\n".join(value))
    # we use a very small chunk size
    got = list(iter_splitlines(path, chunk_size=5))
    assert got == value


@pytest.mark.parametrize("newline", ["\n", "\r\n"])
def test_iter_splitlines_file_endswith_newline(tmp_path, newline):
    path = tmp_path / "multi-line.txt"
    value = ["We have some", "text on different lines", "which load"]
    with open_(path, mode="w", newline=newline) as out:
        out.write("\n".join(value) + "\n")
    # we use a very small chunk size
    got = list(iter_splitlines(path, chunk_size=5))
    assert got == value


def test_iter_splitlines_chunk_size_exceeds_file_size(tmp_path):
    path = tmp_path / "multi-line.txt"
    value = ["We have some", "text on different lines", "which load"]
    path.write_text("\n".join(value))
    # we use a massive chunk size
    got = list(iter_splitlines(path, chunk_size=5_000_000))
    assert got == value


@pytest.mark.parametrize(
    "value",
    [
        # creates a one line block ending on newline
        "With text\nending on a\nended in newline.",
        # creates a two line block ending on newline
        "With text\nending\non a\nended in newline.",
    ],
)
def test_iter_splitlines_chunk_endswith_newline(tmp_path, value):
    path = tmp_path / "multi-line.txt"
    # character 22 is a newline
    value = value.splitlines()
    path.write_text("\n".join(value))
    # we use a chunk size that ends with a newline
    got = list(iter_splitlines(path, chunk_size=11))
    assert got == value


def test_iter_splitlines_chunk_empty_file(tmp_path):
    path = tmp_path / "zero.txt"
    path.write_text("")
    got = list(iter_splitlines(path))
    assert not got


@pytest.mark.parametrize("transform", [str, pathlib.Path])
def test_iter_splitlines_tilde(home_file, transform):
    expect = pathlib.Path(home_file).expanduser().read_text().splitlines()
    got = list(iter_splitlines(transform(home_file)))
    assert len(got) == len(expect)


def test_iter_line_blocks_correct_size(tmp_path):
    # correctly break up
    path = tmp_path / "multi-line.txt"
    value = ["We have some", "text on different lines", "which load"]
    path.write_text("\n".join(value))
    # we use a massive chunk size
    got = list(iter_line_blocks(path, num_lines=2, chunk_size=5))
    expect = [value[:2], value[-1:]]
    assert got == expect


def test_iter_line_blocks_empty(tmp_path):
    path = tmp_path / "zero.txt"
    path.write_text("")
    # we use a massive chunk size
    got = list(iter_line_blocks(path, num_lines=2))
    assert not got


def test_iter_line_blocks_one(tmp_path):
    # file has a single line
    path = tmp_path / "one-line.txt"
    value = "We have text on one line."
    path.write_text(value)
    got = list(iter_line_blocks(path, num_lines=2))
    assert got == [[value]]


def test_iter_line_blocks_none_num_lines(tmp_path):
    # correctly break up
    path = tmp_path / "multi-line.txt"
    value = ["We have some", "text on different lines", "which load"]
    path.write_text("\n".join(value))
    # we use a massive chunk size
    got = list(iter_line_blocks(path, num_lines=None))
    expect = [value]
    assert got == expect
