from __future__ import annotations

import contextlib
import datetime
import os
import re
import sqlite3

from pathlib import Path
from typing import Optional, Union

from scitrack import get_text_hexdigest

from cogent3.app.data_store_new import (
    _LOG_TABLE,
    APPEND,
    OVERWRITE,
    READONLY,
    DataMember,
    DataMemberABC,
    DataStoreABC,
    Mode,
    StrOrBytes,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.10.31a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

_RESULT_TABLE = "results"
_MEMORY = ":memory:"
_mem_pattern = re.compile(r"^\s*[:]{0,1}memory[:]{0,1}\s*$")
NoneType = type(None)

# create db
def open_sqlite_db_rw(path: Union[str, Path]):
    """creates a new sqlitedb for read/write at path, can be an in-memory db

    Notes
    -----
    This function embeds the schema. There are three tables:

    - results: analysis objects, may be completed or not completed
    - logs: log-file contents
    - state: whether db is locked to a process

    Returns
    -------
    Handle to a sqlite3 session
    """
    db = sqlite3.connect(
        path,
        isolation_level=None,
        detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES,
    )
    db.row_factory = sqlite3.Row
    create_template = "CREATE TABLE IF NOT EXISTS {};"
    # note it is essential to use INTEGER for the autoincrement of primary key to work
    creates = [
        "state(state_id INTEGER PRIMARY KEY, record_type TEXT, lock_pid INTEGER)",
        f"{_LOG_TABLE}(log_id INTEGER PRIMARY KEY, log_name TEXT, date timestamp, data BLOB)",
        f"{_RESULT_TABLE}(record_id TEXT PRIMARY KEY, log_id INTEGER, md5 BLOB, is_completed INTEGER, data BLOB)",
    ]
    for table in creates:
        db.execute(create_template.format(table))
    return db


def open_sqlite_db_ro(path):
    """returns db opened as read only
    Returns
    -------
    Handle to a sqlite3 session
    """
    db = sqlite3.connect(
        f"file:{path}?mode=ro",
        isolation_level=None,
        detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES,
        uri=True,
    )
    db.row_factory = sqlite3.Row
    assert has_valid_schema(db)
    return db


def has_valid_schema(db):
    # todo: should be a full schema check
    query = "SELECT name FROM sqlite_master WHERE type='table'"
    result = db.execute(query).fetchall()
    table_names = {r["name"] for r in result}
    return table_names == {_RESULT_TABLE, _LOG_TABLE, "state"}


class DataStoreSqlite(DataStoreABC):
    store_suffix = "sqlitedb"

    def __init__(
        self,
        source,
        mode: Mode = READONLY,
        limit=None,
        verbose=False,
    ):
        if _mem_pattern.search(str(source)):
            self._source = _MEMORY
        else:
            source = Path(source).expanduser()
            self._source = (
                source
                if source.suffix[1:] == self.store_suffix  # sliced to remove "."
                else Path(f"{source}.{self.store_suffix}")
            )
        self._mode = Mode(mode)
        if mode is not READONLY and limit is not None:
            raise ValueError(
                "Using limit argument is only valid for readonly datastores"
            )
        self._limit = limit
        self._verbose = verbose
        self._db = None
        self._open = False
        self._log_id = None

    def __getstate__(self):
        return {**self._init_vals}

    def __setstate__(self, state):
        # this will reset connections to read only db's
        obj = self.__class__(**state)
        self.__dict__.update(obj.__dict__)

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

    @property
    def db(self):
        if self._db is None:
            db_func = open_sqlite_db_ro if self.mode is READONLY else open_sqlite_db_rw
            self._db = db_func(self.source)
            self._open = True
            self.lock()

        return self._db

    def _init_log(self):
        timestamp = datetime.datetime.now()
        self.db.execute(f"INSERT INTO {_LOG_TABLE}(date) VALUES (?)", (timestamp,))
        self._log_id = self._db.execute(
            f"SELECT log_id FROM {_LOG_TABLE} where date = ?", (timestamp,)
        ).fetchone()["log_id"]

    def close(self):
        with contextlib.suppress(sqlite3.ProgrammingError):
            self.db.close()
        self._open = False

    def read(self, identifier: str) -> StrOrBytes:
        """
        identifier string formed from Path(table_name) / identifier
        """
        identifier = Path(identifier)
        table_name = str(identifier.parent)
        if table_name not in (
            ".",
            _LOG_TABLE,
        ):
            raise ValueError(f"unknown table for {str(identifier)!r}")

        if table_name != _LOG_TABLE:
            cmnd = f"SELECT * FROM {_RESULT_TABLE} WHERE record_id = ?"
            result = self.db.execute(cmnd, (identifier.name,)).fetchone()
            return result["data"]

        cmnd = f"SELECT * FROM {_LOG_TABLE} WHERE log_name = ?"
        result = self.db.execute(cmnd, (identifier.name,)).fetchone()

        return result["data"]

    @property
    def completed(self):
        if not self._completed:
            self._completed = self._select_members(
                table_name=_RESULT_TABLE, is_completed=True
            )
        return self._completed

    @property
    def not_completed(self):
        """returns database records of type NotCompleted"""
        if not self._not_completed:
            self._not_completed = self._select_members(
                table_name=_RESULT_TABLE, is_completed=False
            )
        return self._not_completed

    def _select_members(
        self, *, table_name: str, is_completed: bool
    ) -> list[DataMemberABC]:
        limit = f"LIMIT {self.limit}" if self.limit else ""
        cmnd = self.db.execute(
            f"SELECT record_id FROM {table_name} WHERE is_completed=? {limit}",
            (is_completed,),
        )
        return [
            DataMember(data_store=self, unique_id=r["record_id"])
            for r in cmnd.fetchall()
        ]

    @property
    def logs(self):
        """returns all log records"""
        cmnd = self.db.execute(f"SELECT log_name FROM {_LOG_TABLE}")
        return [
            DataMember(data_store=self, unique_id=Path(_LOG_TABLE) / r["log_name"])
            for r in cmnd.fetchall()
            if r["log_name"]
        ]

    def _write(
        self, *, table_name: str, unique_id: str, data: str, is_completed: bool
    ) -> DataMemberABC | None:
        """
        Parameters
        ----------
        table_name: str
            name of table to save data. It must be _RESULT_TABLE or _LOG_TABLE.
        unique_id : str
            unique identifier that data will be saved under.
        data: str
            data to be saved.
        is_completed: bool
            flag to identify NotCompleted results

        Returns
        -------
        DataMember instance or None when writing to _LOG_TABLE
        """
        if self._log_id is None:
            self._init_log()

        if table_name == _LOG_TABLE:
            # todo how to evaluate whether writing a new log?
            cmnd = f"UPDATE {table_name} SET data =?, log_name =? WHERE log_id=?"
            self.db.execute(cmnd, (data, unique_id, self._log_id))
            return None

        md5 = get_text_hexdigest(data)

        if unique_id in self and self.mode is not APPEND:
            cmnd = f"UPDATE {table_name} SET data= ?, log_id=?, md5=? WHERE record_id=?"
            self.db.execute(cmnd, (data, self._log_id, md5, unique_id))
        else:
            cmnd = f"INSERT INTO {table_name} (record_id,data,log_id,md5,is_completed) VALUES (?,?,?,?,?)"
            self.db.execute(cmnd, (unique_id, data, self._log_id, md5, is_completed))

        return DataMember(data_store=self, unique_id=unique_id)

    def drop_not_completed(self) -> None:
        self.db.execute(f"DELETE FROM {_RESULT_TABLE} WHERE is_completed=0")
        self._not_completed = []

    @property
    def _lock_id(self) -> Optional[int]:
        """returns lock_pid"""
        result = self.db.execute("SELECT lock_pid FROM state").fetchone()
        return result[0] if result else result

    @property
    def locked(self) -> bool:
        """returns if lock_pid is NULL or doesn't exist."""
        return self._lock_id is not None

    def lock(self):
        """if writable, and not locked, locks the database to this pid"""
        # if mode=w and the data store exists AND has a lock_pid
        # value already, then we should fail. The user might
        # inadvertently overwrite something otherwise.
        # BUT if mode=a, as the user expects to be modifying the data
        # store then we have to update the value
        if self.mode is READONLY:
            return

        result = self._db.execute("SELECT state_id,lock_pid FROM state").fetchall()
        locked = result[0]["lock_pid"] if result else None
        if locked and self.mode is OVERWRITE:
            raise IOError(
                f"You are trying to OVERWRITE {str(self.source)!r} which is "
                "locked. Use APPEND mode or unlock."
            )

        if result:
            # we will update an existing
            state_id = result[0]["state_id"]
            cmnd = "UPDATE state SET lock_pid=? WHERE state_id=?"
            vals = [os.getpid(), state_id]
        else:
            cmnd = "INSERT INTO state(lock_pid) VALUES (?)"
            vals = [os.getpid()]
        self._db.execute(cmnd, tuple(vals))

    def unlock(self, force=False):
        """remove a lock if pid matches. If force, ignores pid. ignored if mode is READONLY"""
        if self.mode is READONLY:
            return

        lock_id = self._lock_id
        if lock_id is None:
            return

        if lock_id == os.getpid() or force:
            self.db.execute("UPDATE state SET lock_pid=NULL WHERE state_id=1")

        return

    def write(self, *, unique_id: str, data: StrOrBytes) -> DataMemberABC:
        if unique_id.startswith(_RESULT_TABLE):
            unique_id = Path(unique_id).name

        super().write(unique_id=unique_id, data=data)
        member = self._write(
            table_name=_RESULT_TABLE,
            unique_id=unique_id,
            data=data,
            is_completed=True,
        )
        if (
            member is not None and member not in self._completed
        ):  # new to check existence
            self._completed.append(member)
        return member

    def write_log(self, *, unique_id: str, data: StrOrBytes) -> None:
        if unique_id.startswith(_LOG_TABLE):
            unique_id = Path(unique_id).name

        super().write_log(unique_id=unique_id, data=data)
        _ = self._write(
            table_name=_LOG_TABLE, unique_id=unique_id, data=data, is_completed=False
        )

    def write_not_completed(self, *, unique_id: str, data: StrOrBytes) -> DataMemberABC:
        if unique_id.startswith(_RESULT_TABLE):
            unique_id = Path(unique_id).name

        super().write_not_completed(unique_id=unique_id, data=data)
        member = self._write(
            table_name=_RESULT_TABLE, unique_id=unique_id, data=data, is_completed=False
        )
        if member is not None:
            self._not_completed.append(member)
        return member

    def md5(self, unique_id: str) -> Union[str, NoneType]:  # we have it in base class
        """
        Parameters
        ----------
        unique_id
            name of data store member
        Returns
        -------
        md5 checksum for the member, if available, None otherwise
        """
        cmnd = f"SELECT * FROM {_RESULT_TABLE} WHERE record_id = ?"
        result = self.db.execute(cmnd, (unique_id,)).fetchone()

        return result["md5"] if result else None

    @property
    def describe(self):
        if self.locked and self._lock_id != os.getpid():
            title = f"Locked db store. Locked to pid={self._lock_id}, current pid={os.getpid()}."
        elif self.locked:
            title = "Locked to the current process."
        else:
            title = "Unlocked db store."
        table = super().describe
        table.title = title
        return table

    @property
    def record_type(self) -> str:
        """class name of completed results"""
        result = self.db.execute("SELECT record_type FROM state").fetchone()
        return result["record_type"]

    @record_type.setter
    def record_type(self, obj):
        from cogent3.util.misc import get_object_provenance

        rt = self.record_type
        if self.mode is OVERWRITE and rt:
            raise IOError(f"cannot overwrite existing record_type {rt}")

        n = get_object_provenance(obj)
        self.db.execute("UPDATE state SET record_type=? WHERE state_id=1", (n,))
