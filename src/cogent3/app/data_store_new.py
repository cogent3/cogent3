from __future__ import annotations

import contextlib
import inspect
import re
import reprlib
import shutil
import traceback

from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from functools import singledispatch
from pathlib import Path
from typing import Optional, Union

from scitrack import get_text_hexdigest

from cogent3.app.data_store import _db_lockid
from cogent3.app.typing import TabularType
from cogent3.core.alignment import SequenceCollection
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


class Mode(Enum):
    r = "r"
    w = "w"
    a = "a"


APPEND = Mode.a
OVERWRITE = Mode.w
READONLY = Mode.r


@dataclass
class DataMemberABC(ABC):
    data_store: DataStoreABC
    unique_id: str

    #    @property
    #    @abstractmethod
    #    def data_store(self) -> DataStoreABC:
    #        ...

    #   @property
    #   @abstractmethod
    #   def unique_id(self):
    #       ...

    def __str__(self):
        return self.unique_id

    def __repr__(self):
        return f"{self.__class__.__name__}( data_store={self.data_store.source}, unique_id={self.unique_id})"

    def read(self) -> StrOrBytes:
        return self.data_store.read(self.unique_id)

    def write(self, data: StrOrBytes) -> None:
        self.data_store.write(unique_id=self.unique_id, data=data)

    def __eq__(self, other):
        """to check equality of members and check existence of a
        member in a list of members"""
        return self.unique_id == other.unique_id

    @property
    def md5(self):
        return self.data_store.md5(self.unique_id)


@dataclass
class DataStoreABC(ABC):
    source: Union[str, Path]
    mode: Mode
    limit: int

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

    #     @property
    #     @abstractmethod
    #     def source(self) -> str | Path:
    #         """string that references connecting to data store, override in subclass constructor"""
    #         ...

    # @source.setter
    # @abstractmethod
    # def source(self, value):
    #     ...

    #     @property
    #     @abstractmethod
    #     def mode(self) -> Mode:
    #         """string that references datastore mode, override in override in subclass constructor"""
    #         ...

    # @mode.setter
    # @abstractmethod
    # def mode(self, value):
    #     ...

    # @property
    # @abstractmethod
    # def limit(self):
    #     ...

    def __repr__(self):
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
        return any(Path(identifier) == Path(m.unique_id) for m in self)

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
    def write_not_completed(self, unique_id: str, data: StrOrBytes) -> None:
        self._check_writable(unique_id)

    @abstractmethod
    def write_log(self, unique_id: str, data: StrOrBytes) -> None:
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
                    data["composable function"],
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
        err_pat = re.compile(r"[A-Z][a-z]+[A-Z][a-z]+\:.+")
        types = defaultdict(list)
        indices = "type", "origin"
        for member in self.not_completed:
            record = member.read()
            record = deserialise_object(record)
            key = tuple(getattr(record, k, None) for k in indices)
            match = err_pat.findall(record.message)
            types[key].append([match[-1] if match else record.message, record.source])

        header = list(indices) + ["message", "num", "source"]
        rows = []
        maxtring = reprlib.aRepr.maxstring
        reprlib.aRepr.maxstring = 45

        for record in types:
            messages, sources = list(zip(*types[record]))
            messages = reprlib.repr(
                ", ".join(m.splitlines()[-1] for m in set(messages))
            )
            sources = reprlib.repr(", ".join(s.splitlines()[-1] for s in sources))
            row = list(record) + [
                messages,
                len(types[record]),
                sources,
            ]
            rows.append(row)

        reprlib.aRepr.maxstring = maxtring  # restoring original val
        return Table(header=header, data=rows, title="not completed records")

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
        for m in self.members:
            data = m.read()
            md5 = self.md5(m.unique_id)
            if md5 != get_text_hexdigest(data):
                correct_md5 -= 1
        incorrect_md5 = len(self.members) - correct_md5

        return Table(
            header=["Condition", "Value"],
            data=[
                ["Num md5sum correct", correct_md5],
                ["Num md5sum incorrect", incorrect_md5],
                ["Has log", len(self.logs) > 0],
            ],
            title="validate status",
            index_name="Condition",
        )

    def md5(self, unique_id: str) -> str:
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


class DataMember(DataMemberABC):
    def __init__(self, *, data_store: DataStoreABC, unique_id: str):
        self.data_store = data_store
        self.unique_id = str(unique_id)


class DataStoreDirectory(DataStoreABC):
    def __init__(
        self,
        source: Union[str, Path],
        mode: Mode = READONLY,
        suffix: Optional[str] = None,
        limit: int = None,
        verbose=False,
    ):
        self.mode = Mode(mode)
        suffix = suffix or ""
        if suffix != "*":  # wild card search for all
            suffix = re.sub(r"^[\s.*]+", "", suffix)  # tidy the suffix
        source = Path(source)
        self.source = source.expanduser()
        self.suffix = suffix
        self._verbose = verbose
        self._source_check_create(mode)
        self.limit = limit

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

    def write_not_completed(self, unique_id: str, data: str) -> DataMember:
        member = self._write(
            subdir=_NOT_COMPLETED_TABLE, unique_id=unique_id, suffix="json", data=data
        )
        if member is not None:
            self._not_completed.append(member)
        return member

    def write_log(self, unique_id: str, data: str) -> None:
        _ = self._write(subdir=_LOG_TABLE, unique_id=unique_id, suffix="log", data=data)


@singledispatch
def get_data_source(data) -> str:
    source = getattr(data, "source", None)
    if source is None:
        raise NotImplementedError(
            f"Cannot resolve a unique identifier from {type(data)}"
        )
    return source


@get_data_source.register
def _(data: SequenceCollection):
    return get_data_source(data.info.source)


@get_data_source.register
def _(data: str):
    return get_data_source(Path(data))


@get_data_source.register
def _(data: Path):
    return str(data.name)


@get_data_source.register
def _(data: dict):
    try:
        source = data["info"]["source"]
    except KeyError:
        source = data.get("source", None)
    return get_data_source(source)


@get_data_source.register
def _(data: DataMemberABC):
    return str(data.unique_id)


def upgrade_data_store(
    inpath: Path, outpath: Path, insuffix, suffix: Optional[str] = None
) -> DataStoreABC:
    from cogent3.app.io import get_data_store
    from cogent3.app.sqlite_data_store import DataStoreSqlite

    insuffix = insuffix.lower()
    suffix = suffix.lower()

    if suffix == ".tinydb":
        raise "cogent3 does not support tinydb. You can implement your own DataStore derived from DataStoreABC"
    in_dstore = get_data_store(base_path=inpath, suffix=insuffix)

    klass = ""
    if suffix == ".sqlitedb":
        klass = DataStoreSqlite
    elif suffix is None or suffix == ".fasta" or suffix == "":
        klass = DataStoreDirectory

    out_dstore = klass(source=outpath, mode=OVERWRITE, suffix=suffix)
    for member in in_dstore:
        out_dstore.write(unique_id=member.name, data=member.read())

    return out_dstore


def convert_tinydb_to_sqlite(source: Path, dest: Optional[Path] = None) -> DataStoreABC:
    try:
        from fnmatch import fnmatch, translate

        from tinydb import Query, TinyDB
        from tinydb.middlewares import CachingMiddleware
        from tinydb.storages import JSONStorage

        from cogent3.app.data_store import load_record_from_json
        from cogent3.app.sqlite_data_store import DataStoreSqlite
    except ImportError as e:
        raise ImportError(
            "You need to install tinydb to be able to migrate to new datastore."
        ) from e

    storage = CachingMiddleware(JSONStorage)
    storage.WRITE_CACHE_SIZE = 50  # todo support for user specifying
    tinydb = TinyDB(str(source), storage=storage)
    pattern = translate("*")
    query = Query().identifier.matches(pattern)

    id_list = []
    data_list = []
    for record in tinydb.search(query):
        id_list.append(record["identifier"])
        data_list.append(load_record_from_json(record))

    dest = dest or Path(source.parent) / f"{source.stem}.sqlitedb"
    if dest.exists():
        raise IOError(
            f"Destination file {str(dest)} already exists. Delete or define new dest."
        )
    dstore = DataStoreSqlite(source=dest, mode=OVERWRITE)

    for id, data, is_completed in data_list:
        if id == "LOCK":
            cmnd = f"UPDATE state SET lock_pid =? WHERE state_id == 1"
            dstore.db.execute(cmnd, (data,))
            # todo make sure _lockid of sqlitedb matches lockid from tinydb
            assert dstore._lock_id == _db_lockid(str(source))
        else:
            if is_completed:
                dstore.write(unique_id=id, data=data)
            else:
                dstore.write_not_completed(unique_id=id, data=data)

    return dstore
