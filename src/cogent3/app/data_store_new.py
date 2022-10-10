from __future__ import annotations

import contextlib
import re
import reprlib
import shutil

from abc import ABC, abstractmethod
from collections import defaultdict
from enum import Enum
from functools import singledispatch
from pathlib import Path

from scitrack import get_text_hexdigest

from cogent3.app.typing import TabularType
from cogent3.core.alignment import SequenceCollection
from cogent3.util.deserialise import deserialise_not_completed
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


class IfExist(Enum):
    SKIP = "skip"
    OVERWRITE = "overwrite"
    RAISE = "raise"
    READONLY = "readonly"


SKIP = IfExist.SKIP
OVERWRITE = IfExist.OVERWRITE
RAISE = IfExist.RAISE
READONLY = IfExist.READONLY


class DataMemberABC(ABC):
    def __init__(self, data_store: "DataStoreABC" = None):
        self.data_store = data_store

    def __str__(self):
        return self.unique_id

    @property
    @abstractmethod
    def unique_id(self):
        ...

    def read(self) -> str | bytes:
        return self.data_store.read(self.unique_id)

    def write(self, data: str) -> None:
        self.data_store.write(self.unique_id, data)


class DataStoreABC(ABC):
    def __init__(
        self,
        source: str | Path,
        if_dest_exists: IfExist = READONLY,
        if_member_exists: IfExist = RAISE,
    ):
        self._if_dest_exists = IfExist(if_dest_exists)
        self._if_member_exists = IfExist(if_member_exists)

    @abstractmethod
    def __repr__(self):
        ...

    def __getitem__(self, index):
        return self.members[index]

    def __len__(self):
        return len(self.members)

    def __contains__(self, identifier):
        """whether relative identifier has been stored"""
        return any(identifier == m.unique_id for m in self)

    @abstractmethod
    def read(self, unique_id: str) -> str | bytes:
        ...

    def _check_writable(self, unique_id: str):
        if self._if_dest_exists is READONLY:
            raise IOError("datastore is readonly")
        if self._if_member_exists in (READONLY, RAISE) and unique_id in self:
            raise IOError(
                f"DataStore member_exists is {self._if_member_exists!r} and {unique_id!r} exists"
            )

    @abstractmethod
    def write(self, unique_id: str, data: str | bytes) -> None:
        self._check_writable(unique_id)

    @abstractmethod
    def write_not_completed(self, unique_id: str, data: str | bytes) -> None:
        self._check_writable(unique_id)

    @abstractmethod
    def write_log(self, unique_id: str, data: str | bytes) -> None:
        self._check_writable(unique_id)

    @abstractmethod
    def describe(self) -> TabularType:
        ...

    @property
    @abstractmethod
    def members(self) -> list[DataMemberABC]:
        ...

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
    @abstractmethod
    def summary_logs(self) -> TabularType:
        ...

    @property
    @abstractmethod
    def summary_not_completed(self) -> TabularType:
        ...

    @abstractmethod
    def drop_not_completed(self):
        ...


class DataMember(DataMemberABC):
    def __init__(self, data_store: DataStoreABC = None, unique_id: str = None):
        super().__init__(data_store)
        self._unique_id = str(unique_id)

    def __repr__(self):
        return self._unique_id

    @property
    def unique_id(self):
        return self._unique_id

    @property
    def source(self):
        return self.unique_id

    @property
    def md5(self):
        return self.data_store.md5(self.unique_id)


class DataMember(DataMemberABC):
    def __init__(self, data_store: DataStoreABC = None, unique_id: str = None):
        super().__init__(data_store)
        self._unique_id = unique_id

    def __repr__(self):
        return self._unique_id

    @property
    def unique_id(self):
        return self._unique_id


class DataMember(DataMemberABC):
    def __init__(self, data_store: DataStoreABC = None, unique_id: str = None):
        super().__init__(data_store)
        self._unique_id = unique_id

    def __repr__(self):
        return self._unique_id

    @property
    def unique_id(self):
        return self._unique_id


class DataStoreDirectory(DataStoreABC):
    def __init__(
        self,
        source,
        if_dest_exists: IfExist = READONLY,
        if_member_exists: IfExist = RAISE,
        suffix=None,
        limit=None,
        verbose=False,
        md5=True,
    ):
        super().__init__(source, if_dest_exists, if_member_exists)
        suffix = suffix or ""
        if suffix != "*":  # wild card search for all
            suffix = re.sub(r"^[\s.*]+", "", suffix)  # tidy the suffix
        source = Path(source)
        self.source = source.expanduser()
        self.suffix = suffix
        self._completed = []
        self._not_completed = []
        self._limit = limit
        self._verbose = verbose
        self._md5 = md5
        self._source_check_create(if_dest_exists)

    def _source_check_create(self, if_dest_exists):
        if not is_master_process():
            return
        sub_dirs = [_NOT_COMPLETED_TABLE, _LOG_TABLE, _MD5_TABLE]
        source = self.source
        if if_dest_exists is READONLY and not source.exists():
            raise IOError(f"'{source}' does not exist")
        elif if_dest_exists is RAISE and source.exists():
            raise IOError(f"'{source}' exists")
        elif if_dest_exists is READONLY:
            return

        if source.exists() and if_dest_exists is OVERWRITE:
            for sub_dir in sub_dirs:
                with contextlib.suppress((FileNotFoundError, NotADirectoryError)):
                    shutil.rmtree(source / sub_dir)
        elif not source.exists():
            source.mkdir(parents=True, exist_ok=True)

        for sub_dir in sub_dirs:
            (source / sub_dir).mkdir(parents=True, exist_ok=True)

    def read(self, unique_id: str) -> str:
        """reads data corresponding to identifier"""
        with open_(self.source / unique_id) as infile:
            data = infile.read()

        return data

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
            [DataMember(self, Path(_LOG_TABLE) / m.name) for m in log_dir.glob("*")]
            if log_dir.exists()
            else []
        )

    @property
    def completed(self) -> list[DataMember]:
        if not self._completed:
            self._completed = []
            for i, m in enumerate(self.source.glob(f"*.{self.suffix}")):
                if self._limit and i == self._limit:
                    break
                self._completed.append(DataMember(self, m.name))
        return self._completed

    @property
    def not_completed(self) -> list[DataMember]:
        if not self._not_completed:
            self._not_completed = []
            for i, m in enumerate((self.source / _NOT_COMPLETED_TABLE).glob(f"*.json")):
                if self._limit and i == self._limit:
                    break
                self._not_completed.append(
                    DataMember(self, Path(_NOT_COMPLETED_TABLE) / m.name)
                )
        return self._not_completed

    @property
    def members(self) -> list[DataMember]:
        return self.completed + self.not_completed

    # need to update write methods to append to _completed and _not_completed

    def _write(
        self, subdir: str, unique_id: str, suffix: str, data: str, md5: bool
    ) -> DataMember:
        super().write(unique_id, data)
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
            return

        with open_(self.source / subdir / unique_id, mode="w") as out:
            out.write(data)

        if subdir == _NOT_COMPLETED_TABLE:
            member = DataMember(self, Path(_NOT_COMPLETED_TABLE) / unique_id)
        elif not subdir:
            member = DataMember(self, unique_id)
        else:
            assert subdir == _LOG_TABLE
            member = None

        if md5:
            md5 = get_text_hexdigest(data)
            unique_id = unique_id.replace(suffix, "txt")
            unique_id = unique_id if cmp is None else unique_id.replace(f".{cmp}", "")
            if not (self.source / _MD5_TABLE).exists():
                (self.source / _MD5_TABLE).mkdir(parents=True, exist_ok=True)
            with open_(self.source / _MD5_TABLE / unique_id, mode="w") as out:
                out.write(md5)

        return member

    def write(self, unique_id: str, data: str) -> DataMember:
        member = self._write("", unique_id, self.suffix, data, True)
        if member is not None:
            self._completed.append(member)
        return member

    def write_not_completed(self, unique_id: str, data: str) -> DataMember:
        member = self._write(_NOT_COMPLETED_TABLE, unique_id, "json", data, True)
        if member is not None:
            self._not_completed.append(member)
        return member

    def write_log(self, unique_id: str, data: str) -> None:
        _ = self._write(_LOG_TABLE, unique_id, "log", data, False)

    def describe(self) -> TabularType:
        num_notcompleted = len(self.not_completed)
        num_completed = len(self.completed)
        num_logs = len(self.logs)
        return Table(
            header=["record type", "number"],
            data=[
                ["completed", num_completed],
                ["not_completed", num_notcompleted],
                ["logs", num_logs],
            ],
        )

    @property
    def summary_logs(self) -> TabularType:
        """returns a table summarising log files"""
        rows = []
        for record in self.logs:
            data = record.read().splitlines()
            first = data.pop(0).split("\t")
            row = [first[0], record.name]
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
            record = deserialise_not_completed(record)
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

    def __repr__(self):
        num = len(self.members)
        name = self.__class__.__name__
        sample = f"{list(self[:2])}..." if num > 2 else list(self)
        return f"{num}x member {name}(source='{self.source}', members={sample})"


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
