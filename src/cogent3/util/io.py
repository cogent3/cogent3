import contextlib
import functools
import shutil
import uuid
from bz2 import open as bzip_open
from gzip import open as gzip_open
from io import TextIOWrapper
from lzma import open as lzma_open
from os import PathLike, remove
from pathlib import Path, PurePath
from tempfile import mkdtemp
from typing import IO, Callable, Iterator, Optional, Tuple, Union
from urllib.parse import ParseResult, urlparse
from urllib.request import urlopen
from zipfile import ZipFile

from chardet import detect

from cogent3.util.misc import _wout_period

PathType = Union[str, PathLike, PurePath]


@functools.singledispatch
def is_url(path: Union[str, bytes, Path]) -> bool:
    """whether a path is a url"""
    return False


@is_url.register
def _(path: str) -> bool:
    return is_url(urlparse(path))


@is_url.register
def _(path: bytes) -> bool:
    return is_url(urlparse(path.decode("utf8")))


@is_url.register
def _(path: ParseResult) -> bool:
    return path.scheme in {"http", "https", "file"}


def _get_compression_open(
    path: Optional[PathType] = None, compression: Optional[str] = None
) -> Optional[Callable]:
    """returns function for opening compression formats

    Parameters
    ----------
    path
        file path or url
    compression
        file compression suffix

    Returns
    -------
    function for opening compressed files or None if unknown compression
    """
    assert path or compression
    if compression is None:
        _, compression = get_format_suffixes(path)
    return _compression_handlers.get(compression, None)


def open_zip(filename: PathType, mode: str = "r", **kwargs) -> IO:
    """open a single member zip-compressed file

    Note
    ----
    If mode="r". The function raises ValueError if zip has > 1 record.
    The returned object is wrapped by TextIOWrapper with latin encoding
    (so it's not a bytes string).

    If mode="w", returns an atomic_write() instance.
    """
    # import of standard library io module as some code quality tools
    # confuse this with a circular import
    mode = mode or "r"
    binary_mode = "b" in mode
    mode = mode[:1]

    encoding = kwargs.pop("encoding") if "encoding" in kwargs else "latin-1"
    if mode.startswith("w"):
        return atomic_write(filename, mode=mode, in_zip=True)

    mode = mode.strip("t")
    with ZipFile(filename) as zf:
        if len(zf.namelist()) != 1:
            raise ValueError("Archive is supposed to have only one record.")

        opened = zf.open(zf.namelist()[0], mode=mode, **kwargs)

        return opened if binary_mode else TextIOWrapper(opened, encoding=encoding)


_compression_handlers = {
    "gz": gzip_open,
    "bz2": bzip_open,
    "zip": open_zip,
    "xz": lzma_open,
    "lzma": lzma_open,
}


def open_(filename: PathType, mode="rt", **kwargs) -> IO:
    """open that handles different compression

    Parameters
    ----------
    filename
        path or url, if a url delegates processing to open_url
    mode
        standard file opening mode
    kwargs
        passed to open functions

    Returns
    -------
    an object compatible with the file protocol
    """
    if not filename:
        raise ValueError(f"{filename} not a valid file name or url")

    if is_url(filename):
        return open_url(filename, mode=mode, **kwargs)

    mode = mode or "rt"
    filename = Path(filename).expanduser()
    op = _get_compression_open(filename) or open

    encoding = kwargs.pop("encoding", None)
    need_encoding = mode.startswith("r") and "b" not in mode
    if need_encoding and "encoding" not in kwargs:
        with op(filename, mode="rb") as infile:
            data = infile.read(100)

        encoding = detect(data)
        encoding = encoding["encoding"]

    return op(filename, mode, encoding=encoding, **kwargs)


def open_url(url: Union[str, ParseResult], mode="r", **kwargs) -> IO:
    """open a url

    Parameters
    ----------
    url : urllib.parse.ParseResult or str
        A url of file in http or https web address
    mode : str
        mode of reading file, 'rb', 'rt', 'r'

    Raises
    ------
    If not http(s).

    Notes
    -----
    If mode='b' or 'rb' (binary read), the function returns file object to read
    else returns TextIOWrapper to read text with specified encoding in the URL
    """
    _, compression = get_format_suffixes(
        getattr(url, "path", url)
    )  # handling possibility of ParseResult
    mode = "rb" if compression else mode or "r"

    if not is_url(url):
        raise ValueError(
            f"URL scheme must be http, https or file, not {str(url)[:20]!r}"
        )

    url_parsed = url if isinstance(url, ParseResult) else urlparse(url)

    if "r" not in mode:
        raise ValueError("opening a url only allowed in read mode")

    response = urlopen(url_parsed.geturl(), timeout=10)
    encoding = response.headers.get_content_charset()
    if compression:
        response = _get_compression_open(compression=compression)(response)
        mode = "r"  # wrap it as text

    return response if "b" in mode else TextIOWrapper(response, encoding=encoding)


def _path_relative_to_zip_parent(zip_path, member_path):
    """returns member_path relative to zip_path

    Parameters
    ----------
    zip_path: Path
    member_path: Path

    Notes
    -----
    with zip_path = "parentdir/named.zip", then member_path="named/member.tsv"
    or path="member.tsv" will return "named/member.tsv"
    """
    zip_name = zip_path.name.replace(".zip", "")
    if zip_name not in member_path.parts:
        return Path(zip_name) / member_path

    return Path(*member_path.parts[member_path.parts.index(zip_name) :])


class atomic_write:
    """performs atomic write operations, cleans up if fails"""

    def __init__(
        self, path: PathType, tmpdir=None, in_zip=None, mode="w", encoding=None
    ):
        """

        Parameters
        ----------
        path
            path to file, or relative to directory specified by in_zip
        tmpdir
            directory where temporary file will be created
        in_zip
            path to the zip archive containing path,
            e.g. if in_zip="path/to/data.zip", then path="data/seqs.tsv"
            Decompressing the archive will produce the "data/seqs.tsv"
        mode
            file writing mode
        encoding
            text encoding
        """
        path = Path(path).expanduser()
        in_zip = Path(in_zip) if isinstance(in_zip, str) else in_zip
        _, cmp = get_format_suffixes(path)
        if in_zip and cmp == "zip":
            in_zip = path if isinstance(in_zip, bool) else in_zip
            path = Path(str(path)[: str(path).rfind(".zip")])

        if in_zip:
            path = _path_relative_to_zip_parent(in_zip, path)

        self._path = path
        self._cmp = cmp
        self._mode = mode
        self._file = None
        self._encoding = encoding
        self._in_zip = in_zip
        self._tmppath = self._make_tmppath(tmpdir)

        self.succeeded = None
        self._close_func = (
            self._close_rename_zip if in_zip else self._close_rename_standard
        )

    def _make_tmppath(self, tmpdir):
        """returns path of temporary file

        Parameters
        ----------
        tmpdir: Path
            to directory

        Returns
        -------
        full path to a temporary file

        Notes
        -----
        Uses a random uuid as the file name, adds suffixes from path
        """
        suffixes = (
            "".join(self._path.suffixes)
            if not self._in_zip
            else "".join(self._path.suffixes[:-1])
        )
        parent = self._in_zip.parent if self._in_zip else self._path.parent
        name = f"{uuid.uuid4()}{suffixes}"
        tmpdir = Path(mkdtemp(dir=parent)) if tmpdir is None else Path(tmpdir)

        if not tmpdir.exists():
            raise FileNotFoundError(f"{tmpdir} directory does not exist")

        tmp_path = tmpdir / name
        return tmp_path

    def _get_fileobj(self):
        """returns file to be written to"""
        if self._file is None:
            self._file = open_(self._tmppath, self._mode, encoding=self._encoding)

        return self._file

    def __enter__(self) -> IO:
        return self._get_fileobj()

    def _close_rename_standard(self, src):
        dest = Path(self._path)
        try:
            dest.unlink()
        except FileNotFoundError:
            pass
        finally:
            src.rename(dest)

        shutil.rmtree(src.parent)

    def _close_rename_zip(self, src):
        with ZipFile(self._in_zip, "a") as out:
            out.write(str(src), arcname=self._path)

        shutil.rmtree(src.parent)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._file.close()
        if exc_type is None:
            self._close_func(self._tmppath)
            self.succeeded = True
        else:
            self.succeeded = False
            shutil.rmtree(self._tmppath.parent)

    def write(self, text):
        """writes text to file"""
        fileobj = self._get_fileobj()
        fileobj.write(text)

    def close(self):
        """closes file"""
        self.__exit__(None, None, None)


T = Optional[str]


def get_format_suffixes(filename: PathType) -> Tuple[T, T]:
    """returns file, compression suffixes"""
    filename = Path(filename)
    if not filename.suffix:
        return None, None

    suffixes = [_wout_period.sub("", sfx).lower() for sfx in filename.suffixes[-2:]]
    if suffixes[-1] in _compression_handlers:
        cmp_suffix = suffixes[-1]
    else:
        cmp_suffix = None

    if len(suffixes) == 2 and cmp_suffix is not None:
        suffix = suffixes[0]
    elif cmp_suffix is None:
        suffix = suffixes[-1]
    else:
        suffix = None
    return suffix, cmp_suffix


def remove_files(list_of_filepaths, error_on_missing=True):
    """Remove list of filepaths, optionally raising an error if any are missing"""
    missing = []
    for fp in list_of_filepaths:
        try:
            remove(fp)
        except OSError:
            missing.append(fp)

    if error_on_missing and missing:
        raise OSError("Some filepaths were not accessible: %s" % "\t".join(missing))


def path_exists(path: PathType) -> bool:
    """whether path is a valid path and it exists"""
    with contextlib.suppress(Exception):
        return Path(path).exists()
    return False


def iter_splitlines(
    path: PathType, chunk_size: Optional[int] = 1_000_000
) -> Iterator[str]:
    """yields line from file

    Parameters
    ----------
    path
        data file
    chunk_size
        number of bytes to load in one go from path

    Notes
    -----
    Loads chunks of data from the file, yields one line at a time
    """
    if is_url(path):
        chunk_size = None
    else:
        path = Path(path).expanduser()
        if chunk_size and path.stat().st_size < chunk_size:
            # file is smaller than provided chunk_size, just
            # load it all
            chunk_size = None

    with open_(path) as infile:
        last = ""
        while True:
            data = infile.read(chunk_size)
            if not data:  # end of file
                break

            data = last + data
            end_is_newline = data.endswith("\n")
            lines = data.splitlines()
            last = lines.pop(-1)
            if end_is_newline:
                # even if text is from Windows and uses "\r\n", pythons
                # string splitlines() will respect \n
                last += "\n"

            if not len(lines):
                # we have not seen a newline
                continue

            yield from lines

        if last:
            yield from last.splitlines()


def iter_line_blocks(
    path: PathType,
    num_lines: Optional[int] = 1000,
    chunk_size: Optional[int] = 5_000_000,
) -> Iterator[list[str]]:
    """yields list with num_lines str from path

    Parameters
    ----------
    path
        data file
    num_lines
        number of lines per block. If None just returns all lines.
    chunk_size
        number of bytes to load in one go from path
    """
    lines = []
    for line in iter_splitlines(path, chunk_size=chunk_size):
        lines.append(line)
        if len(lines) == num_lines:
            yield lines
            lines = []

    if lines:
        yield lines
