from __future__ import annotations

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
        source,
        if_dest_exists: IfExist = READONLY,
        if_member_exists: IfExist = RAISE,
    ):
        self._if_dest_exists = if_dest_exists
        self._if_member_exists = if_member_exists

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

    @abstractmethod
    def write(self, unique_id: str, data: str | bytes) -> None:
        if self._if_dest_exists is READONLY:
            raise IOError("datastore is readonly")
        if self._if_member_exists in (READONLY, RAISE) and unique_id in self:
            raise IOError(
                f"DataStore member_exists is {self._if_member_exists!r} and {unique_id!r} exists"
            )

    @abstractmethod
    def write_not_completed(self, unique_id: str, data: str | bytes) -> None:
        if self.if_member_exists is RAISE and unique_id in self:
            raise IOError("member exists")

    @abstractmethod
    def write_log(self, unique_id: str, data: str | bytes) -> None:
        if self.if_member_exists is RAISE and unique_id in self:
            raise IOError("member exists")

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
        self._members = []
        self._limit = limit
        self._verbose = verbose
        self._checksums = {}
        self._md5 = md5
        self._source_check_create(if_dest_exists)

    def _source_check_create(self, if_dest_exists):
        if not is_master_process():
            return
        sub_dirs = [_NOT_COMPLETED_TABLE, _LOG_TABLE, _MD5_TABLE]
        if if_dest_exists is READONLY and not self.source.exists():
            raise IOError("must exist")
        elif if_dest_exists is RAISE and self.source.exists():
            raise IOError("must not exist")
        elif if_dest_exists is READONLY:
            return

        path = Path(self.source)
        if path.exists() and if_dest_exists == RAISE:
            raise FileExistsError(f"'{self.source}' exists")
        elif path.exists() and if_dest_exists == OVERWRITE:
            for sub_dir in sub_dirs:
                try:
                    shutil.rmtree(self.source / sub_dir)
                except:
                    pass
        elif not path.exists():
            path.mkdir(parents=True, exist_ok=True)

        for sub_dir in sub_dirs:
            (self.source / sub_dir).mkdir(parents=True, exist_ok=True)

    def read(self, unique_id: str) -> str:
        """reads data corresponding to identifier"""
        # before open should make absolute path
        with open_(self.source / unique_id) as input:
            data = input.read()
            if self._md5 and isinstance(data, str):
                self._checksums[unique_id] = get_text_hexdigest(data)

        return data

    def md5(self, unique_id: str, force=True) -> str:
        """
        Parameters
        ----------
        identifier
            name of data store member
        force : bool
            forces reading of data if not already done

        Returns
        -------
        md5 checksum for the member, if available, None otherwise
        """
        md5_setting = self._md5  # for restoring automatic md5 calc setting
        if force and unique_id not in self._checksums:
            self._md5 = True
            _ = self.read(unique_id)

        self._md5 = md5_setting
        return self._checksums.get(unique_id, None)

    def drop_not_completed(self):
        nc_dir = self.source / _NOT_COMPLETED_TABLE
        md5_dir = self.source / _MD5_TABLE
        for file in list(nc_dir.glob("*.json")):
            file.unlink()
            md5_file = md5_dir / f"{file.stem}.txt"
            md5_file.unlink()

    @property
    def members(self) -> list[DataMember]:
        if not self._members:  # members in completed and not_completed
            self._members = []
            for path in list(self.source.glob(f"*.{self.suffix}")) + list(
                (self.source / _NOT_COMPLETED_TABLE).glob("*.json")
            ):
                if self._limit and len(self._members) >= self._limit:
                    break
                self._members.append(DataMember(self, Path(path).name))

        return self._members

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
        path = Path(self.source)
        return (
            [DataMember(self, path / m.name) for m in path.glob(f"*.{self.suffix}")]
            if path.exists()
            else []
        )

    @property
    def not_completed(self) -> list[DataMember]:
        nc_dir = self.source / _NOT_COMPLETED_TABLE
        return (
            [
                DataMember(self, Path(_NOT_COMPLETED_TABLE) / m.name)
                for m in nc_dir.glob("*.json")
            ]
            if nc_dir.exists()
            else []
        )

    def _write(
        self, subdir: str, unique_id: str, suffix: str, data: str, md5: bool
    ) -> None:
        super().write(unique_id, data)
        assert suffix, "Must provide suffix"
        # check suffix compatible with this datastore
        sfx, cmp = get_format_suffixes(unique_id)
        if sfx != suffix:
            unique_id = f"{unique_id}.{suffix}"
        unique_id = (
            unique_id.replace(self.suffix, suffix)
            if self.suffix and self.suffix != suffix
            else unique_id
        )
        with open_(self.source / subdir / unique_id, mode="w") as out:
            out.write(data)
        if suffix != "log" and not self._limit or len(self) < self._limit:
            self._members.append(DataMember(self, unique_id))
        if md5:
            md5 = get_text_hexdigest(data)
            unique_id = unique_id.replace(suffix, "txt")
            unique_id = unique_id if cmp is None else unique_id.replace(f".{cmp}", "")
            with open_(self.source / _MD5_TABLE / unique_id, mode="w") as out:
                out.write(md5)

    def write(self, unique_id: str, data: str) -> None:
        self._write("", unique_id, self.suffix, data, True)

    def write_not_completed(self, unique_id: str, data: str) -> None:
        self._write(_NOT_COMPLETED_TABLE, unique_id, "json", data, True)

    def write_log(self, unique_id: str, data: str) -> None:
        self._write(_LOG_TABLE, unique_id, "log", data, False)

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
        num_completed = len(self)
        num_not_completed = len(self.not_completed)
        name = self.__class__.__name__
        sample = f"{list(self[:2])}..." if num_completed > 2 else list(self)
        return f"{num_completed + num_not_completed}x member {name}(source='{self.source}', members={sample})"


@singledispatch
def get_data_source(data) -> str:
    ...


@get_data_source.register
def get_str_source(data: str):
    return data


@get_data_source.register
def get_str_source(data: Path):
    return str(data)


@get_data_source.register
def get_str_source(data: dict):
    return str(data.get("source", None))


@get_data_source.register
def get_str_source(data: DataMemberABC):
    return str(data.unique_id)
