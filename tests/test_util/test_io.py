import bz2
import gzip
import pathlib
import tempfile
import zipfile
from urllib.parse import urlparse

import pytest

from cogent3.app.composable import NotCompleted
from cogent3.util.io import (
    _path_relative_to_zip_parent,
    atomic_write,
    get_format_suffixes,
    is_url,
    iter_line_blocks,
    iter_splitlines,
    open_,
    open_url,
    path_exists,
    remove_files,
)


@pytest.fixture
def tmp_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("test_io")


@pytest.fixture
def home_file(DATA_DIR) -> str:
    """makes a temporary directory with file"""
    import tempfile

    HOME = pathlib.Path("~")
    fn = "sample.tsv"
    contents = (DATA_DIR / fn).read_text()
    with tempfile.TemporaryDirectory(dir=HOME.expanduser()) as dn:
        outpath = HOME / pathlib.Path(dn).name / fn
        outpath.expanduser().write_text(contents)
        yield str(outpath)


@pytest.mark.parametrize("transform", (str, pathlib.Path))
def test_open_home(DATA_DIR, home_file, transform):
    """expands tilde for opening / writing to home"""
    data_path = DATA_DIR / "sample.tsv"
    expect = data_path.read_text()
    with open_(transform(home_file)) as infile:
        got = infile.read()
        assert got == expect


def test_does_not_write_if_exception(tmp_dir):
    """file does not exist if an exception raised before closing"""
    test_filepath = tmp_dir / "Atomic_write_test"
    with pytest.raises(AssertionError):
        with atomic_write(test_filepath, mode="w") as f:
            f.write("abc")
            raise AssertionError
    assert not test_filepath.exists()


@pytest.mark.parametrize("suffix", ("gz", "bz2", "zip", "lmza", "xz"))
def test_writes_compressed_formats(DATA_DIR, tmp_dir, suffix):
    """correctly writes / reads different compression formats"""
    fpath = DATA_DIR / "sample.tsv"
    expect = pathlib.Path(fpath).read_text()
    outpath = tmp_dir / f"{fpath.name}.{suffix}"
    with atomic_write(outpath, mode="wt") as f:
        f.write(expect)

    with open_(outpath) as infile:
        got = infile.read()

    assert got == expect, f"write failed for {suffix}"


def test_rename(tmp_dir):
    """Renames file as expected"""
    test_filepath = tmp_dir / "Atomic_write_test"
    # touch the filepath so it exists
    f = open(test_filepath, "w").close()
    assert test_filepath.exists()
    # file should overwrite file if file already exists
    with atomic_write(test_filepath, mode="w") as f:
        f.write("abc")


def test_atomic_write_noncontext(tmp_dir):
    """atomic write works as more regular file object"""
    path = tmp_dir / "foo.txt"
    zip_path = path.parent / f"{path.name}.zip"
    aw = atomic_write(path, in_zip=zip_path, mode="w")
    aw.write("some data")
    aw.close()
    with open_(zip_path) as ifile:
        got = ifile.read()
    assert got == "some data"


def test_open_handles_bom(tmp_dir):
    """handle files with a byte order mark"""
    text = "some text"

    # plain text
    textfile = tmp_dir / "sample.txt"
    textfile.write_text(text, encoding="utf-8-sig")

    # gzipped
    gzip_file = tmp_dir / "sample.txt.gz"
    with gzip.open(gzip_file, "wt", encoding="utf-8-sig") as outfile:
        outfile.write(text)

    # bzipped
    bzip_file = tmp_dir / "sample.txt.bz2"
    with bz2.open(bzip_file, "wt", encoding="utf-8-sig") as outfile:
        outfile.write(text)

    # zipped
    zip_file = tmp_dir / "sample.zip"
    with zipfile.ZipFile(zip_file, "w") as outfile:
        outfile.write(textfile, "sample.txt")

    for path in (bzip_file, gzip_file, textfile, zip_file):
        with open_(path) as infile:
            got = infile.read()
            assert got == text, f"failed reading {path}"


@pytest.mark.parametrize("non", (None, ""))
def test_open_empty_raises(non):
    with pytest.raises(ValueError):
        open_(non)


def test_aw_zip_from_path(tmp_dir):
    """supports inferring zip archive name from path"""
    path = tmp_dir / "foo.txt"
    zip_path = path.parent / f"{path.name}.zip"
    aw = atomic_write(zip_path, in_zip=True, mode="w")
    aw.write("some data")
    aw.close()
    with open_(zip_path) as ifile:
        got = ifile.read()
        assert got == "some data"

    path = tmp_dir / "foo2.txt"
    zip_path = path.parent / f"{path.name}.zip"
    aw = atomic_write(path, in_zip=zip_path, mode="w")
    aw.write("some data")
    aw.close()
    with open_(zip_path) as ifile:
        got = ifile.read()
        assert got == "some data"


def test_expanduser(tmp_dir):
    """expands user correctly"""
    # create temp file directory
    home = pathlib.Path("~").expanduser()
    test_filepath = tmp_dir / "Atomic_write_test"
    test_filepath = str(test_filepath).replace(str(home), "~")
    with atomic_write(test_filepath, mode="w") as f:
        f.write("abc")


def test_path_relative_to_zip_parent():
    """correctly generates member paths for a zip archive"""
    zip_path = pathlib.Path("some/path/to/a/data.zip")
    for member in ("data/member.txt", "member.txt", "a/b/c/member.txt"):
        got = _path_relative_to_zip_parent(zip_path, pathlib.Path(member))
        assert got.parts[0] == "data"


def test_remove_files():
    """Remove files functions as expected"""
    import os

    # create list of temp file paths
    test_filepaths = [
        tempfile.NamedTemporaryFile(prefix="remove_files_test").name for i in range(5)
    ]

    # try to remove them with remove_files and verify that an IOError is
    # raises
    pytest.raises(OSError, remove_files, test_filepaths)
    # now get no error when error_on_missing=False
    remove_files(test_filepaths, error_on_missing=False)

    # touch one of the filepaths so it exists
    open(test_filepaths[2], "w").close()
    # check that an error is raised on trying to remove the files...
    pytest.raises(OSError, remove_files, test_filepaths)
    # ... but that the existing file was still removed
    assert not os.path.exists(test_filepaths[2])

    # touch one of the filepaths so it exists
    open(test_filepaths[2], "w").close()
    # no error is raised on trying to remove the files
    # (although 4 don't exist)...
    remove_files(test_filepaths, error_on_missing=False)
    # ... and the existing file was removed
    assert not os.path.exists(test_filepaths[2])


@pytest.mark.parametrize(
    "name,expect",
    (
        ("suffixes.GZ", (None, "gz")),
        ("suffixes.ABCD", ("abcd", None)),
        ("suffixes.ABCD.BZ2", ("abcd", "bz2")),
        ("suffixes.abcd.BZ2", ("abcd", "bz2")),
        ("suffixes.ABCD.bz2", ("abcd", "bz2")),
    ),
)
def test_get_format_suffixes_returns_lower_case(name, expect):
    """should always return lower case"""
    suffixes = get_format_suffixes(name)
    assert suffixes == expect


@pytest.mark.parametrize(
    "name,expect",
    (
        ("no_suffixes", (None, None)),
        ("suffixes.gz", (None, "gz")),
        ("suffixes.abcd", ("abcd", None)),
        ("suffixes.abcd.bz2", ("abcd", "bz2")),
        ("suffixes.zip", (None, "zip")),
    ),
)
def test_get_format_suffixes(name, expect):
    """correctly return suffixes for compressed etc.. formats"""
    suffixes = get_format_suffixes(name)
    assert suffixes == expect


@pytest.mark.parametrize(
    "name,expect",
    (
        ("no_suffixes", (None, None)),
        ("suffixes.gz", (None, "gz")),
        ("suffixes.abcd", ("abcd", None)),
        ("suffixes.abcd.bz2", ("abcd", "bz2")),
        ("suffixes.zip", (None, "zip")),
    ),
)
def test_get_format_suffixes_pathlib(name, expect):
    """correctly return suffixes for compressed etc.. formats from pathlib"""
    suffixes = get_format_suffixes(pathlib.Path(name))
    assert suffixes == expect


@pytest.mark.parametrize(
    "val,expect",
    (
        ({}, False),
        ("not an existing path", False),
        ("(a,b,(c,d))", False),
        ("(a:0.1,b:0.1,(c:0.1,d:0.1):0.1)", False),
        (__file__, True),
        (pathlib.Path(__file__), True),
        (NotCompleted("FAIL", "test", message="none", source="unknown"), False),
    ),
)
def test_path_exists(val, expect):
    """robustly identifies whether an object is a valid path and exists"""
    assert path_exists(val) == expect


def test_open_reads_zip(tmp_dir):
    """correctly reads a zip compressed file"""
    text_path = tmp_dir / "foo.txt"
    with open(text_path, "w") as f:
        f.write("any str")

    zip_path = tmp_dir / "foo.zip"
    with zipfile.ZipFile(zip_path, "w") as zip:
        zip.write(text_path)

    with open_(zip_path) as got:
        assert got.readline() == "any str"


def test_open_writes_zip(tmp_dir):
    """correctly writes a zip compressed file"""
    zip_path = tmp_dir / "foo.txt.zip"

    with open_(zip_path, "w") as f:
        f.write("any str")

    with zipfile.ZipFile(zip_path, "r") as zip:
        name = zip.namelist()[0]
        got = zip.open(name).read()
        assert got == b"any str"


def test_open_zip_multi(tmp_dir):
    """zip with multiple records cannot be opened using open_"""
    text_path1 = tmp_dir / "foo.txt"
    with open(text_path1, "w") as f:
        f.write("any str")

    text_path2 = tmp_dir / "bar.txt"
    with open(text_path2, "w") as f:
        f.write("any str")

    zip_path = tmp_dir / "foo.zip"
    with zipfile.ZipFile(zip_path, "w") as zip:
        zip.write(text_path1)
        zip.write(text_path2)

    with pytest.raises(ValueError):
        open_(zip_path)


@pytest.mark.parametrize(
    "mode",
    ("r", "rb", "rt", None),
)
@pytest.mark.internet
def test_open_url(DATA_DIR, mode):
    """different open mode's all work"""
    # None value for open_url mode defaults to "r"
    file_name = "gff2_test.gff"
    remote_root = (
        "https://raw.githubusercontent.com/cogent3/cogent3/develop/tests/data/{}"
    )

    with open_(DATA_DIR / file_name, mode=mode) as infile:
        local_data = infile.read()

    with open_url(remote_root.format(file_name), mode=mode) as infile:
        remote_data = infile.read()

    assert remote_data.splitlines() == local_data.splitlines()

    # Test using a ParseResult for url
    with open_url(urlparse(remote_root.format(file_name)), mode=mode) as infile:
        remote_data = infile.read()
    assert remote_data.splitlines() == local_data.splitlines()


def test_open_url_local(DATA_DIR):
    """using file:///"""
    file_name = "gff2_test.gff"
    local_path = DATA_DIR / file_name
    with open_(local_path) as infile:
        local_data = infile.read()

    # make absolute path
    with open_url(local_path.absolute().as_uri()) as infile:
        remote_data = infile.read()

    assert remote_data.splitlines() == local_data.splitlines()


@pytest.mark.internet
def test_open_url_compressed(DATA_DIR):
    """comparing compressed file handling"""
    file_name = "formattest.fasta.gz"
    remote_root = (
        "https://raw.githubusercontent.com/cogent3/cogent3/develop/tests/data/{}"
    )

    with open_(DATA_DIR / file_name) as infile:
        local_data = infile.read()

    with open_url(remote_root.format(file_name)) as infile:
        remote_data = infile.read()

    assert remote_data.splitlines() == local_data.splitlines()


def test_open_url_write_exceptions():
    """Test 'w' mode (should raise Exception)"""
    with pytest.raises(Exception):
        _ = open_url(
            "https://raw.githubusercontent.com/cogent3/cogent3/develop/tests/data/gff2_test.gff",
            mode="w",
        )


def test_open_url_exceptions():
    """non-http(s) address for url (should raise Exception)"""
    with pytest.raises(Exception):
        _ = open_url(
            "ftp://raw.githubusercontent.com/cogent3/cogent3/develop/tests/data/gff2_test.gff",
        )


def test_iter_splitlines_one(tmp_path):
    # file has a single line
    path = tmp_path / "one-line.txt"
    value = "We have text on one line."
    path.write_text(value)
    got = list(iter_splitlines(path))
    assert got == [value]


@pytest.mark.parametrize("newline", ("\n", "\r\n"))
def test_iter_splitlines_line_diff_newline(tmp_path, newline):
    path = tmp_path / "multi-line.txt"
    value = ["We have some", "text on different lines", "which load"]
    with open_(path, mode="w", newline=newline) as out:
        out.write("\n".join(value))
    # we use a very small chunk size
    got = list(iter_splitlines(path, chunk_size=5))
    assert got == value


@pytest.mark.parametrize("newline", ("\n", "\r\n"))
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
    (
        # creates a one line block ending on newline
        "With text\nending on a\nended in newline.",
        # creates a two line block ending on newline
        "With text\nending\non a\nended in newline.",
    ),
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


@pytest.mark.parametrize("transform", (str, pathlib.Path))
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


@pytest.mark.parametrize(
    "url",
    ("https://example.com", pathlib.Path("https://example.com"), b"file://example.txt"),
)
def test_is_url(url):
    assert is_url(url)


@pytest.mark.parametrize(
    "url", ("example.txt", pathlib.Path("example.txt"), b"example.txt")
)
def test_not_is_url(url):
    assert not is_url(url)
