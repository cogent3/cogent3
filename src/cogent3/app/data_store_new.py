import glob
import re
import shutil

from abc import ABC, abstractmethod
from enum import Enum
from functools import singledispatch
from pathlib import Path

from scitrack import get_text_hexdigest

from cogent3.util.io import open_
from cogent3.util.parallel import is_master_process


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

_INCOMPLETE_TABLE = "incomplete"
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

    def read(self) -> str:
        return self.data_store.read(self.unique_id)

    def write(self, data: str) -> None:
        self.data_store.write(self.unique_id, data)


class DataMember(DataMemberABC):
    def __init__(self, data_store: "DataStoreDirectory" = None, unique_id: str = None):
        super().__init__(data_store)
        self._unique_id = unique_id

    def __repr__(self):
        return self._unique_id

    @property
    def unique_id(self):
        return self._unique_id


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
    def read(self, unique_id: str) -> str:
        ...

    @abstractmethod
    def write(self, unique_id: str, data: str) -> None:
        if self._if_dest_exists is READONLY:
            raise IOError("datastore is readonly")
        if self._if_member_exists is RAISE and unique_id in self:
            raise IOError("member exists")

    @abstractmethod
    def write_incomplete(self, unique_id: str, data):
        if self.if_member_exists is RAISE and unique_id in self:
            raise IOError("member exists")

    @abstractmethod
    def write_log(self, unique_id: str, data):
        if self.if_member_exists is RAISE and unique_id in self:
            raise IOError("member exists")

    @abstractmethod
    def describe(self):
        ...

    @property
    @abstractmethod
    def members(self):
        ...

    @abstractmethod
    def open(self, member):
        ...

    def __iter__(self):
        yield from self.members

    @property
    @abstractmethod
    def logs(self):
        ...

    @property
    @abstractmethod
    def incomplete(self):
        ...

    @property
    @abstractmethod
    def summary_logs(self):
        ...

    @property
    @abstractmethod
    def summary_incomplete(self):
        ...


class DataStoreDirectory(DataStoreABC):
    store_suffix = None

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
        self.suffix = suffix
        if self.store_suffix and not source.name.endswith(self.store_suffix):
            raise ValueError(f"{source.name} doesn't match {self.store_suffix}")
        self.source = Path(source).expanduser()
        self._members = []
        self._limit = limit
        self._verbose = verbose
        self._checksums = {}
        self._md5 = md5
        self._source_check_create(if_dest_exists)

    def _source_check_create(self, if_dest_exists):
        if not is_master_process():
            return

        sub_dirs = [_INCOMPLETE_TABLE, _LOG_TABLE, _MD5_TABLE]

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
            (self.source.parent / sub_dir).mkdir(parents=True, exist_ok=True)

    def read(self, unique_id: str) -> str:
        """reads data corresponding to identifier"""
        # before open should make absolute path
        with open_(self.source / unique_id) as input:
            data = input.read()
            if self._md5 and isinstance(data, str):
                self._checksums[unique_id] = get_text_hexdigest(data)

        return data

    def md5(self, unique_id: str, force=True):
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

    @property
    def members(self):
        if not self._members:  # members in completed and incomplete
            pattern = f"{self.source}/**/*.{self.suffix}"
            paths = list(glob.iglob(pattern, recursive=True))
            self._members = []
            for i, path in enumerate(paths):
                if self._limit and i >= self._limit:
                    break
                self._members.append(DataMember(self, Path(path).name))

        return self._members

    @property
    def logs(self):
        log_dir = self.source / _LOG_TABLE
        return (
            [DataMember(self, Path(_LOG_TABLE) / m.name) for m in log_dir.glob("*")]
            if log_dir.exists()
            else []
        )

    @property
    def incomplete(self):
        ic_dir = self.source / _INCOMPLETE_TABLE
        return (
            [
                DataMember(self, Path(_INCOMPLETE_TABLE) / m.name)
                for m in ic_dir.glob("*.json")
            ]
            if ic_dir.exists()
            else []
        )

    def write(self, unique_id: str, data: str) -> None:
        super().write(unique_id, data)
        with open_(self.source / unique_id, mode="w") as out:
            out.write(data)
            if self._md5 and isinstance(data, str):
                self._checksums[unique_id] = get_text_hexdigest(data)

    def write_incomplete(self, unique_id: str, data: str) -> None:
        super().write(unique_id, data)
        with open_(self.source / _INCOMPLETE_TABLE / unique_id, mode="w") as out:
            out.write(data)

    def write_log(self, unique_id: str, data: str) -> None:
        super().write(unique_id, data)
        with open_(self.source / _LOG_TABLE / unique_id, mode="w") as out:
            out.write(data)

    def describe(self):
        ...

    def open(self, member):
        return open_(self.source / member.unique_id, mode="r")

    @property
    def summary_logs(self):
        ...

    @property
    def summary_incomplete(self):
        ...

    def __repr__(self):
        incomplete_path = Path(self.source / _INCOMPLETE_TABLE)
        if not incomplete_path.exists():
            return "This isn't a results directory"
        num_completed += len(self)
        num_incomplete = len(list(incomplete_path.glob("*")))
        name = self.__class__.__name__
        sample = f"{str(list(self[:3]))[:-1]}..." if num_completed > 3 else list(self)
        return f"{num_completed+num_incomplete}x member {name}(source='{self.source}', members={sample})"


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
