import contextlib
import json
import os
import pickle
import zipfile

from enum import Enum
from functools import singledispatch
from gzip import compress as gzip_compress
from gzip import decompress as gzip_decompress
from pathlib import Path
from typing import Any, Optional, Union

import numpy

from cogent3.core.alignment import ArrayAlignment, SequenceCollection
from cogent3.core.moltype import get_moltype
from cogent3.core.profile import (
    make_motif_counts_from_tabular,
    make_motif_freqs_from_tabular,
    make_pssm_from_tabular,
)
from cogent3.evolve.fast_distance import DistanceMatrix
from cogent3.format.alignment import FORMATTERS
from cogent3.parse.sequence import PARSERS
from cogent3.util.deserialise import deserialise_object
from cogent3.util.misc import extend_docstring_from
from cogent3.util.table import Table

from .composable import LOADER, WRITER, NotCompleted, define_app
from .data_store import (
    READONLY,
    DataStoreABC,
    DataStoreDirectory,
    Mode,
    ReadOnlyDataStoreZipped,
    get_unique_id,
    load_record_from_json,
    make_record_for_json,
)
from .sqlite_data_store import _MEMORY, DataStoreSqlite
from .typing import (
    AlignedSeqsType,
    IdentifierType,
    SeqsCollectionType,
    SerialisableType,
    TabularType,
    UnalignedSeqsType,
)


_datastore_reader_map = {}


class register_datastore_reader:
    """
    registration decorator for read only data store classes

    The registration key must be a string that of the file format suffix
    (more than one suffix can be registered at a time).

    Parameters
    ----------
    args: str or sequence of str
        must be unique, a preceding '.' will be added if not already present
    """

    def __init__(self, *args):
        args = list(args)
        for i, suffix in enumerate(args):
            if suffix is None:
                assert (
                    suffix not in _datastore_reader_map
                ), f"{suffix!r} already in {list(_datastore_reader_map)}"
                continue

            if not isinstance(suffix, str):
                raise TypeError(f"{suffix!r} is not a string")

            if suffix.strip() == suffix and not suffix:
                raise ValueError("cannot have white-space suffix")

            suffix = suffix.strip()
            if suffix:
                suffix = suffix if suffix[0] == "." else f".{suffix}"

            assert (
                suffix not in _datastore_reader_map
            ), f"{suffix!r} already in {list(_datastore_reader_map)}"
            args[i] = suffix

        self._type_str = tuple(args)

    def __call__(self, func):
        for type_str in self._type_str:
            _datastore_reader_map[type_str] = func
        return func


# register the main readers
register_datastore_reader("zip")(ReadOnlyDataStoreZipped)
register_datastore_reader(None)(DataStoreDirectory)
register_datastore_reader("sqlitedb")(DataStoreSqlite)


def open_data_store(
    base_path: Union[str, Path],
    suffix: Optional[str] = None,
    limit: Optional[int] = None,
    mode: Union[str, Mode] = READONLY,
    **kwargs,
) -> DataStoreABC:
    """returns DataStore instance of a type specified by the path suffix

    Parameters
    ----------
    base_path
        path to directory or db
    suffix
        suffix of filenames
    limit
        the number of matches to return
    mode
        opening mode, either r, w, a as per file opening modes
    """
    mode = Mode(mode)
    if not isinstance(suffix, (str, type(None))):
        raise ValueError(f"suffix {type(suffix)} not one of string or None")

    kwargs = {"limit": limit, "mode": mode, "suffix": suffix, **kwargs}
    base_path = Path(base_path)
    base_path = (
        base_path if base_path.name == ":memory:" else base_path.expanduser().absolute()
    )
    if base_path.is_dir():
        ds_suffix = None
    elif base_path.suffix == ".sqlitedb" or base_path.name == _MEMORY:
        ds_suffix = ".sqlitedb"
        kwargs.pop("suffix")
    elif zipfile.is_zipfile(base_path):
        ds_suffix = ".zip"
    elif base_path.suffix:
        ds_suffix = base_path.suffix
    else:
        # triggered when mode="w"
        ds_suffix = None

    if base_path.name == _MEMORY and mode is READONLY:
        raise NotImplementedError("in memory readonly sqlitedb")

    if ds_suffix is None and suffix is None:
        raise ValueError("a suffix is required if using a directory data store")

    klass = _datastore_reader_map[ds_suffix]

    return klass(base_path, **kwargs)


@define_app(skip_not_completed=False)
def pickle_it(data: SerialisableType) -> bytes:
    """Serialises data using pickle."""
    return pickle.dumps(data)


@define_app(skip_not_completed=False)
def unpickle_it(data: bytes) -> SerialisableType:
    "Deserialises pickle data."
    return pickle.loads(data)


@define_app(skip_not_completed=False)
class compress:
    """Compresses bytes data."""

    def __init__(self, compressor: callable = gzip_compress):
        """
        Parameters
        ----------
        compressor
            function for compressing bytes data, defaults to gzip
        """
        self.compressor = compressor

    def main(self, data: bytes) -> bytes:
        return self.compressor(data)


@define_app(skip_not_completed=False)
class decompress:
    """Decompresses data."""

    def __init__(self, decompressor: callable = gzip_decompress):
        """
        Parameters
        ----------
        decompressor
            a function for decompression, defaults to the gzip decompress
            function
        """
        self.decompressor = decompressor

    def main(self, data: bytes) -> bytes:
        return self.decompressor(data)


def as_dict(obj: Any) -> dict:
    """returns result of to_rich_dict method if it exists"""
    with contextlib.suppress(AttributeError):
        obj = obj.to_rich_dict()
    return obj


@define_app(skip_not_completed=False)
class to_primitive:
    """convert an object to primitive python types suitable for serialisation"""

    def __init__(self, convertor: callable = as_dict):
        self.convertor = convertor

    def main(self, data: SerialisableType) -> SerialisableType:
        """returns dict from a cogent3 object"""
        return self.convertor(data)


@define_app(skip_not_completed=False)
class from_primitive:
    """deserialises from primitive python types"""

    def __init__(self, deserialiser: callable = deserialise_object):
        self.deserialiser = deserialiser

    def main(self, data: SerialisableType) -> SerialisableType:
        """either json or a dict from a cogent3 object"""
        return self.deserialiser(data)


@define_app
def to_json(data: SerialisableType) -> str:
    """Convert primitive python types to json string."""
    return json.dumps(data)


@define_app
def from_json(data: str) -> SerialisableType:
    """Convert json string to primitive python types."""
    return json.loads(data)


class tabular(Enum):
    table = "table"
    distances = "distances"
    motif_counts = "motif_counts"
    motif_freqs = "motif_freqs"
    pssm = "pssm"


@singledispatch
def _read_it(path) -> str:
    try:
        data = path.read()
    except AttributeError:
        raise IOError(f"unexpected type {type(path)}")
    return data


@_read_it.register
def _(path: os.PathLike) -> str:
    path = path.expanduser().absolute()
    return path.read_text()


@_read_it.register
def _(path: str) -> os.PathLike:
    return _read_it(Path(path))


def _load_seqs(path, klass, parser, moltype):
    data = _read_it(path)
    data = data.splitlines()
    data = dict(iter(parser(data)))
    seqs = klass(data=data, moltype=moltype)
    seqs.info.source = str(path)
    return seqs


@define_app(app_type=LOADER)
class load_aligned:
    """Loads aligned sequences. Returns an Alignment object."""

    def __init__(self, moltype=None, format="fasta"):
        """
        Parameters
        ----------
        moltype
            molecular type, string or instance
        format : str
            sequence file format
        """
        self.moltype = moltype if moltype is None else get_moltype(moltype)
        self._parser = PARSERS[format.lower()]

    T = Union[SerialisableType, AlignedSeqsType]

    def main(self, path: IdentifierType) -> T:
        """returns alignment"""
        return _load_seqs(path, ArrayAlignment, self._parser, self.moltype)


@define_app(app_type=LOADER)
class load_unaligned:
    """Loads unaligned sequences. Returns a SequenceCollection."""

    def __init__(self, *, moltype=None, format="fasta"):
        """
        Parameters
        ----------
        moltype
            molecular type, string or instance
        format : str
            sequence file format
        """
        self.moltype = moltype if moltype is None else get_moltype(moltype)
        self._parser = PARSERS[format.lower()]

    T = Union[SerialisableType, UnalignedSeqsType]

    def main(self, path: IdentifierType) -> T:
        """returns sequence collection"""
        seqs = _load_seqs(path, SequenceCollection, self._parser, self.moltype)
        return seqs.degap()


@define_app(app_type=LOADER)
class load_tabular:
    """Loads delimited data. Returns a Table."""

    def __init__(
        self,
        with_title=False,
        with_header=True,
        limit=None,
        sep="\t",
        strict=True,
        as_type: tabular = "table",
    ):
        """
        Parameters
        ----------
        with_title
            files have a title
        with_header
            files have a header
        limit
            number of records to read
        sep
            field delimiter
        as_type
            cogent3 type of tabular data
        strict
            all rows MUST have the same number of records
        """
        self._sep = sep
        self._with_title = with_title
        self._with_header = with_header
        self._limit = limit
        self.strict = strict
        self.as_type = tabular(as_type)

    def _parse(self, data):
        """returns header, records, title"""
        title = header = None
        sep = self._sep
        strict = self.strict
        lines = _read_it(data).splitlines()
        if self._limit:
            lines = lines[: self._limit]
        if self._with_title:
            title = lines.pop(0).strip()
        if self._with_header:
            header = [e.strip() for e in lines.pop(0).strip().split(sep)]

        num_records = None if header is None else len(header)
        rows = []
        for line in lines:
            line = line.strip()
            line = [e.strip() for e in line.split(sep)]
            if num_records is None:
                num_records = len(line)
            if strict and len(line) != num_records:
                raise AssertionError(
                    f"Inconsistent number of fields: {len(line)} != {num_records}"
                )
            rows.append(line)

        records = []
        for record in zip(*rows):
            record = numpy.array(record, dtype="O")
            try:
                record = record.astype(int)
            except ValueError:
                try:
                    record = record.astype(float)
                except ValueError:
                    pass
            records.append(record)
        records = numpy.array(records, dtype="O").T
        return header, records, title

    def main(self, path: IdentifierType) -> TabularType:
        try:
            header, data, title = self._parse(path)
        except Exception as err:
            return NotCompleted("ERROR", self, err.args[0], source=str(path))

        if self.as_type is tabular.table:
            return Table(header=header, data=data, title=title)

        assert data.shape[1] == 3, "Invalid tabular data"

        if self.as_type is tabular.distances:
            # records is of the form [ [dim-1, dim-2, value] for entries in DistanceMatrix ]
            return DistanceMatrix({(e[0], e[1]): e[2] for e in data})

        func = {
            tabular.motif_counts: make_motif_counts_from_tabular,
            tabular.motif_freqs: make_motif_freqs_from_tabular,
            tabular.pssm: make_pssm_from_tabular,
        }

        return func[self.as_type](data)


@define_app(app_type=LOADER)
class load_json:
    """Loads json serialised cogent3 objects from a json file.
    Returns whatever object type was stored."""

    def main(self, path: IdentifierType) -> SerialisableType:
        """returns object deserialised from json at path"""
        data = _read_it(path)
        identifier, data, completed = load_record_from_json(data)

        result = deserialise_object(data)
        if hasattr(result, "info"):
            result.info["source"] = result.info.get("source", identifier)
        else:
            try:
                identifier = getattr(result, "source", identifier)
                setattr(result, "source", identifier)
            except AttributeError:
                pass
        return result


DEFAULT_DESERIALISER = unpickle_it() + from_primitive()


@define_app(app_type=LOADER)
class load_db:
    """Loads serialised cogent3 objects from a db.
    Returns whatever object type was stored."""

    def __init__(self, deserialiser: callable = DEFAULT_DESERIALISER):
        self.deserialiser = deserialiser

    def main(self, identifier: IdentifierType) -> SerialisableType:
        """returns deserialised object"""
        try:
            data = identifier.read()
        except AttributeError:
            raise AttributeError(
                f"{identifier} failed because its of type {type(identifier)}"
            )

        # do we need to inject identifier attribute?
        result = self.deserialiser(data)
        if hasattr(result, "info"):
            result.info["source"] = result.info.get("source", identifier)
        else:
            with contextlib.suppress(AttributeError):
                identifier = getattr(result, "source", identifier)
                setattr(result, "source", identifier)
        return result


@define_app(app_type=WRITER)
class write_json:
    """Writes data in json format."""

    def __init__(
        self,
        data_store: DataStoreABC,
        id_from_source: callable = get_unique_id,
    ):
        """
        Parameters
        ----------
        data_store
            A writeable data store.
        id_from_source
            A function for creating a unique identifier based on the data
            source. The default function removes path information and
            filename + compression suffixes.
        """
        if not isinstance(data_store, DataStoreABC):
            raise TypeError(f"invalid type {type(data_store)!r} for data_store")
        self.data_store = data_store
        self._format = "json"
        self._id_from_source = id_from_source

    def main(
        self, /, data: SerialisableType, *, identifier: Optional[str] = None
    ) -> IdentifierType:
        identifier = identifier or self._id_from_source(data)
        if isinstance(data, NotCompleted):
            return self.data_store.write_not_completed(
                unique_id=f"{identifier}.json", data=data.to_json()
            )

        out = make_record_for_json(identifier, data, True)
        data = json.dumps(out)
        return self.data_store.write(unique_id=identifier, data=data)


@define_app(app_type=WRITER)
class write_seqs:
    """Write sequences in standard formats."""

    @extend_docstring_from(write_json.__init__, pre=False)
    def __init__(
        self,
        data_store: DataStoreABC,
        id_from_source: callable = get_unique_id,
        format: str = "fasta",
    ):
        """
        format
            sequence format
        """
        if not isinstance(data_store, DataStoreABC):
            raise TypeError(f"invalid type {type(data_store)!r} for data_store")
        self.data_store = data_store
        self._formatter = FORMATTERS[format]
        self._id_from_source = id_from_source

    def main(
        self, /, data: SeqsCollectionType, *, identifier: Optional[str] = None
    ) -> IdentifierType:
        identifier = identifier or self._id_from_source(data)
        if isinstance(data, NotCompleted):
            return self.data_store.write_not_completed(
                unique_id=f"{identifier}.json", data=data.to_json()
            )

        data = self._formatter(data.to_dict())
        return self.data_store.write(unique_id=identifier, data=data)


@define_app(app_type=WRITER)
class write_tabular:
    """Writes tabular data in text format supported by the cogent3 Table object."""

    @extend_docstring_from(write_json.__init__, pre=False)
    def __init__(
        self,
        data_store: DataStoreABC,
        id_from_source: callable = get_unique_id,
        format: str = "tsv",
    ):
        """
        format
            tabular format, e.g. 'csv' or 'tsv'
        """
        if not isinstance(data_store, DataStoreABC):
            raise TypeError(f"invalid type {type(data_store)!r} for data_store")
        self.data_store = data_store
        self._id_from_source = id_from_source
        self._format = format

    def main(
        self, /, data: TabularType, *, identifier: Optional[str] = None
    ) -> IdentifierType:
        identifier = identifier or self._id_from_source(data)
        if isinstance(data, NotCompleted):
            return self.data_store.write_not_completed(
                unique_id=f"{identifier}.json", data=data.to_json()
            )

        output = data.to_string(format=self._format)
        return self.data_store.write(unique_id=identifier, data=output)


DEFAULT_SERIALISER = to_primitive() + pickle_it()


@define_app(app_type=WRITER, skip_not_completed=False)
class write_db:
    """Write serialised objects to a database instance."""

    @extend_docstring_from(write_json.__init__, pre=False)
    def __init__(
        self,
        data_store: DataStoreABC,
        id_from_source: callable = get_unique_id,
        serialiser: callable = DEFAULT_SERIALISER,
    ):
        """
        serialiser
            A callable for serialising input data. By default, it converts
            data into primitive python data types, pickling the result.
        """
        if not isinstance(data_store, DataStoreABC):
            raise TypeError(f"invalid type {type(data_store)!r} for data_store")
        self.data_store = data_store
        self._serialiser = serialiser
        self._id_from_source = id_from_source

    def main(
        self, /, data: SerialisableType, *, identifier: Optional[str] = None
    ) -> IdentifierType:
        identifier = identifier or self._id_from_source(data)
        blob = self._serialiser(data)

        if isinstance(data, NotCompleted):
            return self.data_store.write_not_completed(unique_id=identifier, data=blob)

        if self.data_store.record_type is None:
            self.data_store.record_type = data

        return self.data_store.write(unique_id=identifier, data=blob)
