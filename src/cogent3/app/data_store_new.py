from __future__ import annotations

import contextlib
import inspect
import json
import re
import reprlib

from abc import ABC, abstractmethod
from collections import defaultdict
from enum import Enum
from functools import singledispatch
from pathlib import Path
from typing import Optional, Union

from scitrack import get_text_hexdigest

from cogent3.app.typing import TabularType
from cogent3.core.alignment import (
    Alignment,
    ArrayAlignment,
    SequenceCollection,
)
from cogent3.util.deserialise import deserialise_object
from cogent3.util.io import get_format_suffixes, open_
from cogent3.util.parallel import is_master_process
from cogent3.util.table import Table


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

_NOT_COMPLETED_TABLE = "not_completed"
_LOG_TABLE = "logs"
_MD5_TABLE = "md5"

StrOrBytes = Union[str, bytes]
NoneType = type(None)


class Mode(Enum):
    r = "r"
    w = "w"
    a = "a"


APPEND = Mode.a
OVERWRITE = Mode.w
READONLY = Mode.r


class DataMemberABC(ABC):
    """Abstract base class for DataMember

    A data member is a handle to a record in a DataStore. It has a reference
    to its data store and a unique identifier.
    """

    @property
    @abstractmethod
    def data_store(self) -> DataStoreABC:
        ...

    @property
    @abstractmethod
    def unique_id(self):
        ...

    def __str__(self):
        return self.unique_id

    def __repr__(self):
        return f"{self.__class__.__name__}(data_store={self.data_store.source}, unique_id={self.unique_id})"

    def read(self) -> StrOrBytes:
        return self.data_store.read(self.unique_id)

    def __eq__(self, other):
        """to check equality of members and check existence of a
        member in a list of members"""
        return type(self) == type(other) and (self.data_store, self.unique_id) == (
            other.data_store,
            other.unique_id,
        )

    @property
    def md5(self):
        return self.data_store.md5(self.unique_id)


class DataStoreABC(ABC):
    """Abstract base class for DataStore"""

    def __new__(klass, *args, **kwargs):
        obj = object.__new__(klass)

        init_sig = inspect.signature(klass.__init__)
        bargs = init_sig.bind_partial(klass, *args, **kwargs)
        bargs.apply_defaults()
        init_vals = bargs.arguments
        init_vals.pop("self", None)

        obj._init_vals = init_vals
        obj._completed = []
        obj._not_completed = []
        return obj

    @property
    @abstractmethod
    def source(self) -> str | Path:
        """string that references connecting to data store, override in subclass constructor"""
        ...

    @property
    @abstractmethod
    def mode(self) -> Mode:
        """string that references datastore mode, override in subclass constructor"""
        ...

    @property
    @abstractmethod
    def limit(self):
        ...

    def __repr__(self):
        name = self.__class__.__name__
        construction = ", ".join(f"{k}={v}" for k, v in self._init_vals.items())
        return f"{name}({construction})"

    def __str__(self):
        num = len(self.members)
        name = self.__class__.__name__
        sample = f"{list(self[:2])}..." if num > 2 else list(self)
        return f"{num}x member {name}(source='{self.source}', members={sample})"

    def __getitem__(self, index):
        return self.members[index]

    def __len__(self):
        return len(self.members)

    def __contains__(self, identifier):
        """whether relative identifier has been stored"""
        return any(identifier == m.unique_id for m in self)

    @abstractmethod
    def read(self, unique_id: str) -> StrOrBytes:
        ...

    def _check_writable(self, unique_id: str):
        if self.mode is READONLY:
            raise IOError("datastore is readonly")
        elif unique_id in self and self.mode is APPEND:
            raise IOError("cannot overwrite existing record in append mode")

    @abstractmethod
    def write(self, *, unique_id: str, data: StrOrBytes) -> None:
        self._check_writable(unique_id)

    @abstractmethod
    def write_not_completed(self, *, unique_id: str, data: StrOrBytes) -> None:
        self._check_writable(unique_id)

    @abstractmethod
    def write_log(self, *, unique_id: str, data: StrOrBytes) -> None:
        self._check_writable(unique_id)

    @property
    def members(self) -> list[DataMemberABC]:
        return self.completed + self.not_completed

    def __iter__(self):
        yield from self.members

    @property
    @abstractmethod
    def logs(self) -> list[DataMemberABC]:
        ...

    @property
    @abstractmethod
    def completed(self) -> list[DataMemberABC]:
        ...

    @property
    @abstractmethod
    def not_completed(self) -> list[DataMemberABC]:
        ...

    @property
    def summary_logs(self) -> TabularType:
        """returns a table summarising log files"""
        rows = []
        for record in self.logs:
            data = record.read().splitlines()
            first = data.pop(0).split("\t")
            row = [first[0], record.unique_id]
            key = None
            mapped = {}
            for line in data:
                line = line.split("\t")[-1].split(" : ", maxsplit=1)
                if len(line) == 1:
                    mapped[key] += line[0]
                    continue

                key = line[0]
                mapped[key] = line[1]

            data = mapped
            row.extend(
                [
                    data["python"],
                    data["user"],
                    data["command_string"],
                    data.get("composable function", ""),
                ]
            )
            rows.append(row)
        return Table(
            header=["time", "name", "python version", "who", "command", "composable"],
            data=rows,
            title="summary of log files",
        )

    @property
    def summary_not_completed(self) -> TabularType:
        """returns a table summarising not completed results"""
        # detect last exception line
        return summary_not_completeds(self.not_completed)

    @property
    def describe(self) -> TabularType:
        title = "Directory datastore"
        num_not_completed = len(self.not_completed)
        num_completed = len(self.completed)
        num_logs = len(self.logs)
        return Table(
            header=["record type", "number"],
            data=[
                ["completed", num_completed],
                ["not_completed", num_not_completed],
                ["logs", num_logs],
            ],
            title=title,
        )

    @abstractmethod
    def drop_not_completed(self):
        ...

    def validate(self) -> TabularType:
        correct_md5 = len(self.members)
        missing_md5 = 0
        for m in self.members:
            data = m.read()
            md5 = self.md5(m.unique_id)
            if md5 is None:
                missing_md5 += 1
                correct_md5 -= 1
            elif md5 != get_text_hexdigest(data):
                correct_md5 -= 1

        incorrect_md5 = len(self.members) - correct_md5 - missing_md5

        return Table(
            header=["Condition", "Value"],
            data=[
                ["Num md5sum correct", correct_md5],
                ["Num md5sum incorrect", incorrect_md5],
                ["Num md5sum missing", missing_md5],
                ["Has log", len(self.logs) > 0],
            ],
            title="validate status",
            index_name="Condition",
        )

    @abstractmethod
    def md5(self, unique_id: str) -> Union[str, NoneType]:
        """
        Parameters
        ----------
        unique_id
            name of data store member

        Returns
        -------
        md5 checksum for the member, if available, None otherwise
        """


class DataMember(DataMemberABC):
    """Generic DataMember class, bound to a data store. All read operations
    delivered by the parent."""

    def __init__(self, *, data_store: DataStoreABC, unique_id: str):
        self._data_store = data_store
        self._unique_id = str(unique_id)

    @property
    def data_store(self) -> DataStoreABC:
        return self._data_store

    @property
    def unique_id(self):
        return self._unique_id


def summary_not_completeds(
    not_completed: list[DataMemberABC], deserialise: Optional[callable] = None
) -> Table:
    """
    Parameters
    ----------
    not_completed
        list of DataMember instances for notcompleted records
    deserialise
        a callable for converting not completed contents, the result of member.read() must be a json string
    """
    err_pat = re.compile(r"[A-Z][a-z]+[A-Z][a-z]+\:.+")
    types = defaultdict(list)
    indices = "type", "origin"
    num_bytes = 0
    for member in not_completed:
        record = member.read()
        if deserialise:
            record = deserialise(record)
        if isinstance(record, bytes):
            num_bytes += 1
            continue
        record = deserialise_object(record)
        key = tuple(getattr(record, k, None) for k in indices)
        match = err_pat.findall(record.message)
        types[key].append([match[-1] if match else record.message, record.source])
    header = list(indices) + ["message", "num", "source"]
    if num_bytes == len(not_completed):
        return Table(
            header=header,
            title="Cannot summarise not_completed as they are all bytes, "
            "use an appropriate reader",
        )

    rows = []
    maxtring = reprlib.aRepr.maxstring
    reprlib.aRepr.maxstring = 45
    for record in types:
        messages, sources = list(zip(*types[record]))
        messages = reprlib.repr(", ".join(m.splitlines()[-1] for m in set(messages)))
        sources = reprlib.repr(", ".join(s.splitlines()[-1] for s in sources))
        row = list(record) + [
            messages,
            len(types[record]),
            sources,
        ]
        rows.append(row)
    reprlib.aRepr.maxstring = maxtring  # restoring original val
    return Table(header=header, data=rows, title="not completed records")


class DataStoreDirectory(DataStoreABC):
    def __init__(
        self,
        source: Union[str, Path],
        mode: Mode = READONLY,
        suffix: Optional[str] = None,
        limit: int = None,
        verbose=False,
    ):
        self._mode = Mode(mode)
        suffix = suffix or ""
        if suffix != "*":  # wild card search for all
            suffix = re.sub(r"^[\s.*]+", "", suffix)  # tidy the suffix
        source = Path(source)
        self._source = source.expanduser()
        self.suffix = suffix
        self._verbose = verbose
        self._source_check_create(mode)
        self._limit = limit

    def __contains__(self, item: str):
        item = f"{item}.{self.suffix}" if self.suffix not in item else item
        return super().__contains__(item)

    def _source_check_create(self, mode: Mode) -> None:
        if not is_master_process():
            return

        sub_dirs = [_NOT_COMPLETED_TABLE, _LOG_TABLE, _MD5_TABLE]
        source = self.source
        if mode is READONLY and not source.exists():
            raise IOError(f"'{source}' does not exist")
        elif mode is READONLY:
            return

        if not source.exists():
            source.mkdir(parents=True, exist_ok=True)

        for sub_dir in sub_dirs:
            (source / sub_dir).mkdir(parents=True, exist_ok=True)

    @property
    def source(self) -> str | Path:
        """string that references connecting to data store, override in subclass constructor"""
        return self._source

    @property
    def mode(self) -> Mode:
        """string that references datastore mode, override in override in subclass constructor"""
        return self._mode

    @property
    def limit(self):
        return self._limit

    def read(self, unique_id: str) -> str:
        """reads data corresponding to identifier"""
        with open_(self.source / unique_id) as infile:
            data = infile.read()
        return data

    def drop_not_completed(self):
        nc_dir = self.source / _NOT_COMPLETED_TABLE
        md5_dir = self.source / _MD5_TABLE
        for file in nc_dir.glob("*.json"):
            file.unlink()
            md5_file = md5_dir / f"{file.stem}.txt"
            md5_file.unlink()
        Path(self.source / _NOT_COMPLETED_TABLE).rmdir()
        # reset _not_completed list to force not_completed function to make it again
        self._not_completed = []

    @property
    def logs(self) -> list[DataMember]:
        log_dir = self.source / _LOG_TABLE
        return (
            [
                DataMember(data_store=self, unique_id=Path(_LOG_TABLE) / m.name)
                for m in log_dir.glob("*")
            ]
            if log_dir.exists()
            else []
        )

    @property
    def completed(self) -> list[DataMember]:
        if not self._completed:
            self._completed = []
            for i, m in enumerate(self.source.glob(f"*.{self.suffix}")):
                if self.limit and i == self.limit:
                    break
                self._completed.append(DataMember(data_store=self, unique_id=m.name))
        return self._completed

    @property
    def not_completed(self) -> list[DataMember]:
        if not self._not_completed:
            self._not_completed = []
            for i, m in enumerate((self.source / _NOT_COMPLETED_TABLE).glob(f"*.json")):
                if self.limit and i == self.limit:
                    break
                self._not_completed.append(
                    DataMember(
                        data_store=self, unique_id=Path(_NOT_COMPLETED_TABLE) / m.name
                    )
                )
        return self._not_completed

    def _write(
        self, *, subdir: str, unique_id: str, suffix: str, data: str
    ) -> DataMember:
        super().write(unique_id=unique_id, data=data)
        assert suffix, "Must provide suffix"
        # check suffix compatible with this datastore
        sfx, cmp = get_format_suffixes(unique_id)
        if sfx != suffix:
            unique_id = f"{Path(unique_id).stem}.{suffix}"

        unique_id = (
            unique_id.replace(self.suffix, suffix)
            if self.suffix and self.suffix != suffix
            else unique_id
        )
        if suffix != "log" and unique_id in self:
            return None

        with open_(self.source / subdir / unique_id, mode="w") as out:
            out.write(data)

        if subdir == _LOG_TABLE:
            return None
        if subdir == _NOT_COMPLETED_TABLE:
            member = DataMember(
                data_store=self, unique_id=Path(_NOT_COMPLETED_TABLE) / unique_id
            )
        elif not subdir:
            member = DataMember(data_store=self, unique_id=unique_id)

        md5 = get_text_hexdigest(data)
        unique_id = unique_id.replace(suffix, "txt")
        unique_id = unique_id if cmp is None else unique_id.replace(f".{cmp}", "")
        with open_(self.source / _MD5_TABLE / unique_id, mode="w") as out:
            out.write(md5)

        return member

    def write(self, *, unique_id: str, data: str) -> DataMember:
        member = self._write(
            subdir="", unique_id=unique_id, suffix=self.suffix, data=data
        )
        if member is not None:
            self._completed.append(member)
        return member

    def write_not_completed(self, *, unique_id: str, data: str) -> DataMember:
        (self.source / _NOT_COMPLETED_TABLE).mkdir(parents=True, exist_ok=True)
        member = self._write(
            subdir=_NOT_COMPLETED_TABLE, unique_id=unique_id, suffix="json", data=data
        )
        if member is not None:
            self._not_completed.append(member)
        return member

    def write_log(self, *, unique_id: str, data: str) -> None:
        (self.source / _LOG_TABLE).mkdir(parents=True, exist_ok=True)
        _ = self._write(subdir=_LOG_TABLE, unique_id=unique_id, suffix="log", data=data)

    def md5(self, unique_id: str) -> Union[str, NoneType]:
        """
        Parameters
        ----------
        unique_id
            name of data store member

        Returns
        -------
        md5 checksum for the member, if available, None otherwise
        """
        unique_id = Path(unique_id)
        unique_id = re.sub(rf"[.]({self.suffix}|json)$", ".txt", unique_id.name)
        path = self.source / _MD5_TABLE / unique_id

        return path.read_text() if path.exists() else None


def get_unique_id(name: str) -> str:
    """strips any format suffixes from name"""
    name = get_data_source(name)
    suffixes = ".".join(sfx for sfx in get_format_suffixes(name) if sfx)
    return re.sub(fr"[.]{suffixes}$", "", name)


@singledispatch
def get_data_source(data) -> str:
    source = getattr(data, "source", None)
    if source is None:
        return None
    return get_data_source(source)


@get_data_source.register
def _(data: SequenceCollection):
    return get_data_source(data.info)


@get_data_source.register
def _(data: ArrayAlignment):
    return get_data_source(data.info)


@get_data_source.register
def _(data: Alignment):
    return get_data_source(data.info)


@get_data_source.register
def _(data: str):
    return get_data_source(Path(data))


@get_data_source.register
def _(data: Path):
    return str(data.name)


@get_data_source.register
def _(data: dict):
    try:
        source = data.get("info", {})["source"]
    except KeyError:
        source = data.get("source", None)
    return get_data_source(source)


@get_data_source.register
def _(data: DataMemberABC):
    return str(data.unique_id)


def convert_directory_datastore(
    inpath: Path, outpath: Path, suffix: Optional[str] = None
) -> DataStoreABC:
    out_dstore = DataStoreDirectory(source=outpath, mode=OVERWRITE, suffix=suffix)
    filenames = inpath.glob(f"*{suffix}")
    for fn in filenames:
        out_dstore.write(unique_id=fn.name, data=fn.read_text())
    return out_dstore


def convert_tinydb_to_sqlite(source: Path, dest: Optional[Path] = None) -> DataStoreABC:
    from datetime import datetime
    from fnmatch import translate

    from .composable import CachingLogger, _make_logfile_name
    from .data_store import load_record_from_json
    from .io_new import write_db
    from .sqlite_data_store import _LOG_TABLE, DataStoreSqlite

    try:
        from tinydb import Query, TinyDB
        from tinydb.middlewares import CachingMiddleware
        from tinydb.storages import JSONStorage
    except ImportError as e:
        raise ImportError(
            "You need to install tinydb to be able to migrate to new datastore."
        ) from e

    source = Path(source)
    storage = CachingMiddleware(JSONStorage)
    tinydb = TinyDB(str(source), storage=storage)
    pattern = translate("*")
    query = Query().identifier.matches(pattern)

    id_list = []
    data_list = []
    lock_id = None
    for record in tinydb.search(query):
        if record["identifier"] == "LOCK":
            lock_id = record["pid"]
            continue
        id_list.append(record["identifier"])
        data_list.append(load_record_from_json(record))

    dest = dest or Path(source.parent) / f"{source.stem}.sqlitedb"
    if dest.exists():
        raise IOError(
            f"Destination file {str(dest)} already exists. Delete or define new dest."
        )

    LOGGER = CachingLogger(create_dir=True)
    log_file_path = source.parent / _make_logfile_name("convert_tinydb_to_sqlite")
    LOGGER.log_file_path = log_file_path

    dstore = DataStoreSqlite(source=dest, mode=OVERWRITE)
    writer = write_db(data_store=dstore)
    for id, data, is_completed in data_list:
        if id.endswith(".log"):
            cmnd = f"UPDATE {_LOG_TABLE} SET data =?, log_name =?"
            values = (data, id)
            with contextlib.suppress(ValueError):
                date = datetime.strptime(
                    data.split("\t", maxsplit=1)[0], "%Y-%m-%d %H:%M:%S"
                )
                cmnd = f"{cmnd}, date=?"
                values += (date,)

            cmnd = f"{cmnd} WHERE log_id=?"
            values += (dstore._log_id,)

            dstore.db.execute(cmnd, values)
        else:
            writer.main(data, identifier=id)

    # add a new log, recording this conversion
    LOGGER.shutdown()
    dstore.close()
    dstore = DataStoreSqlite(source=dest, mode=APPEND)
    dstore.write_log(unique_id=log_file_path.name, data=log_file_path.read_text())
    log_file_path.unlink()
    if lock_id is not None or dstore._lock_id:
        cmnd = f"UPDATE state SET lock_pid =? WHERE state_id == 1"
        dstore.db.execute(cmnd, (lock_id,))

    return dstore


def make_record_for_json(identifier, data, completed):
    """returns a dict for storage as json"""
    try:
        data = data.to_rich_dict()
    except AttributeError:
        pass

    data = json.dumps(data)
    return dict(identifier=identifier, data=data, completed=completed)


def load_record_from_json(data):
    """returns identifier, data, completed status from json string"""
    if type(data) == str:
        data = json.loads(data)

    value = data["data"]
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except json.JSONDecodeError:
            pass

    return data["identifier"], value, data["completed"]
