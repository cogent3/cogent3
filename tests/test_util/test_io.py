import bz2
import gzip
import pathlib
import tempfile
import zipfile

from urllib.parse import urlparse

import pytest

from cogent3.util.io import (
    _path_relative_to_zip_parent,
    atomic_write,
    get_format_suffixes,
    open_,
    open_url,
    path_exists,
    remove_files,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

DATADIR = pathlib.Path(__file__).parent.parent / "data"


@pytest.fixture
def tmp_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("test_io")


@pytest.fixture
def home_file() -> str:
    """makes a temporary directory with file"""
    import tempfile

    HOME = pathlib.Path("~")
    fn = "sample.tsv"
    contents = (DATADIR / fn).read_text()
    with tempfile.TemporaryDirectory(dir=HOME.expanduser()) as dn:
        outpath = HOME / pathlib.Path(dn).name / fn
        outpath.expanduser().write_text(contents)
        yield str(outpath)


@pytest.mark.parametrize("transform", (str, pathlib.Path))
def test_open_home(home_file, transform):
    """expands tilde for opening / writing to home"""
    data_path = DATADIR / "sample.tsv"
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


def test_writes_compressed_formats(tmp_dir):
    """correctly writes / reads different compression formats"""
    fpath = DATADIR / "sample.tsv"
    expect = pathlib.Path(fpath).read_text()
    for suffix in ["gz", "bz2", "zip"]:
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
def test_open_url(mode):
    """different open mode's all work"""
    # None value for open_url mode defaults to "r"
    file_name = "gff2_test.gff"
    remote_root = (
        "https://raw.githubusercontent.com/cogent3/cogent3/develop/tests/data/{}"
    )

    with open_(DATADIR / file_name, mode=mode) as infile:
        local_data = infile.read()

    with open_url(remote_root.format(file_name), mode=mode) as infile:
        remote_data = infile.read()

    assert remote_data.splitlines() == local_data.splitlines()

    # Test using a ParseResult for url
    with open_url(urlparse(remote_root.format(file_name)), mode=mode) as infile:
        remote_data = infile.read()
    assert remote_data.splitlines() == local_data.splitlines()


def test_open_url_local():
    """using file:///"""
    file_name = "gff2_test.gff"
    local_path = DATADIR / file_name
    with open_(local_path) as infile:
        local_data = infile.read()

    # make absolute path
    with open_url(local_path.absolute().as_uri()) as infile:
        remote_data = infile.read()

    assert remote_data.splitlines() == local_data.splitlines()


def test_open_url_compressed():
    """comparing compressed file handling"""
    file_name = "formattest.fasta.gz"
    remote_root = (
        "https://raw.githubusercontent.com/cogent3/cogent3/develop/tests/data/{}"
    )

    with open_(DATADIR / file_name) as infile:
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
