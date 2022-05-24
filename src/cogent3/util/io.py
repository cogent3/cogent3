import shutil
import uuid

from bz2 import open as bzip_open
from gzip import open as gzip_open
from os import path as os_path
from os import remove
from pathlib import Path
from tempfile import mkdtemp
from typing import Union
from zipfile import ZipFile

from chardet import detect

from cogent3.util.misc import _wout_period


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


def open_zip(filename: Union[str, Path], mode="r", **kwargs):
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
    from io import TextIOWrapper

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

        if binary_mode:
            return opened

        return TextIOWrapper(opened, encoding=encoding)


def open_(filename: Union[str, Path], mode="rt", **kwargs):
    """open that handles different compression"""

    filename = Path(filename).expanduser().absolute()
    op = {".gz": gzip_open, ".bz2": bzip_open, ".zip": open_zip}.get(
        filename.suffix, open
    )

    encoding = kwargs.pop("encoding", None)
    need_encoding = mode.startswith("r") and "b" not in mode
    if need_encoding:
        if "encoding" not in kwargs:
            with op(filename, mode="rb") as infile:
                data = infile.read(100)

            encoding = detect(data)
            encoding = encoding["encoding"]

    return op(filename, mode, encoding=encoding, **kwargs)


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
        self, path: Union[str, Path], tmpdir=None, in_zip=None, mode="w", encoding=None
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

    def __enter__(self):
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


def get_format_suffixes(filename: Union[str, Path]):
    """returns file, compression suffixes"""
    filename = Path(filename)
    if not filename.suffix:
        return None, None

    compression_suffixes = ("bz2", "gz", "zip")
    suffixes = [_wout_period.sub("", sfx).lower() for sfx in filename.suffixes[-2:]]
    if suffixes[-1] in compression_suffixes:
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


def path_exists(path):
    """whether path is a valid path and it exists"""
    if not (isinstance(path, str) or isinstance(path, Path)):
        return False
    try:
        is_path = os_path.exists(str(path))
    except (ValueError, TypeError):
        is_path = False
    return is_path
