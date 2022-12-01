from __future__ import annotations

import contextlib
import datetime
import os
import re
import sqlite3

from pathlib import Path
from typing import Union

from scitrack import get_text_hexdigest

from cogent3.app.data_store_new import (
    _LOG_TABLE,
    APPEND,
    READONLY,
    DataMember,
    DataMemberABC,
    DataStoreABC,
    Mode,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

_RESULT_TABLE = "result_table"
_MEMORY = ":memory:"
_mem_pattern = re.compile(r"^\s*[:]{0,1}memory[:]{0,1}\s*$")

# create db
def open_sqlite_db_rw(path: Union[str, Path]):
    """creates a new sqlitedb for read/write at path, can be an in-memory db

    Notes
    -----
    This function embeds the schema. There are four tables:
    - complete: completed analysis objects
    - incomplete: NotCompleted objects
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
        f"{_RESULT_TABLE}(record_id TEXT PRIMARY KEY, log_id INTEGER, md5 BLOB, not_completed INTEGER, data BLOB)",
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
    return db


class SqliteDataStore(DataStoreABC):
    store_suffix = "sqlitedb"

    def __init__(
        self,
        source,
        mode: Mode,
        limit=None,
        verbose=False,
    ):
        if _mem_pattern.search(str(source)):
            self.source = _MEMORY
        else:
            source = Path(source).expanduser()
            self.source = (
                source
                if source.suffix[1:] == self.store_suffix  # sliced to remove "."
                else Path(f"{source}.{self.store_suffix}")
            )

        super().__init__(
            source=source,
            mode=mode,
        )
        if mode is not READONLY and limit is not None:
            raise ValueError(
                "Using limit argument is only valid for readonly datastores"
            )
        self._limit = limit
        self._verbose = verbose
        self._checksums = {}
        self._db = None
        self._members = []
        self._completed = []
        self._not_completed = []
        self._open = False
        self._log_id = None

    def __getstate__(self):
        return {**self._init_vals}

    def __setstate__(self, state):
        # this will reset connections to read only db's
        obj = self.__class__(**state)
        self.__dict__.update(obj.__dict__)

    def __contains__(self, unique_id: str):
        """whether relative identifier has been stored"""
        # todo add limit here
        query = f"SELECT count(*) as c FROM {_RESULT_TABLE} where record_id == ?"
        result = self.db.execute(query, (unique_id,)).fetchone()
        if result["c"] > 0:
            return True
        return False

    @property
    def db(self):
        if self._db is None:
            db_func = open_sqlite_db_ro if self._mode is READONLY else open_sqlite_db_rw
            self._db = db_func(self.source)
            self._open = True

            if self._mode is not READONLY:
                pid = os.getpid()
                self._db.execute("INSERT INTO state(lock_pid) VALUES (?)", (pid,))

        # todo: lock_id comes from process id, into state table  #new  check this on write

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

    @property
    def _lock_id(self):
        """returns lock_pid"""
        result = self.db.execute("SELECT lock_pid FROM state").fetchone()
        return result[0] if result else result

    def read(self, identifier: str) -> str | bytes:
        """
        identifier string formed from Path(table_name) / identifier
        """
        identifier = Path(identifier)
        table_name = str(identifier.parent)
        if table_name not in (
            _RESULT_TABLE,
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
                table_name=_RESULT_TABLE, not_completed=False
            )
        return self._completed

    @property
    def not_completed(self):
        """returns database records of type NotCompleted"""
        if not self._not_completed:
            self._not_completed = self._select_members(
                table_name=_RESULT_TABLE, not_completed=True
            )
        return self._not_completed

    def _select_members(
        self, table_name: str, not_completed: bool
    ) -> list[DataMemberABC]:
        limit = f"LIMIT {self._limit}" if self._limit else ""
        cmnd = self.db.execute(
            f"SELECT record_id FROM {table_name} WHERE not_completed=? {limit}",
            (not_completed,),
        )
        return [
            DataMember(data_store=self, unique_id=Path(table_name) / r["record_id"])
            for r in cmnd.fetchall()
        ]

    @property
    def logs(self):
        """returns all log records"""
        cmnd = self.db.execute(f"SELECT data, log_name FROM {_LOG_TABLE}")
        return [
            DataMember(data_store=self, unique_id=Path(_LOG_TABLE) / r["log_name"])
            for r in cmnd.fetchall()
        ]

    def _write(
        self, *, table_name: str, unique_id: str, data: str, not_completed: bool
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
        not_completed: bool
            flag to identify not_completed results

        Returns
        -------
        DataMember instance or None when writing to _LOG_TABLE
        """
        if self._log_id is None:
            self._init_log()
        if table_name == _LOG_TABLE:
            cmnd = f"UPDATE {table_name} SET data =?, log_name =? WHERE log_id=?"
            self.db.execute(cmnd, (data, unique_id, self._log_id))
            return None
        md5 = get_text_hexdigest(data)

        if unique_id in self and self._mode is not APPEND:
            cmnd = f"UPDATE {table_name} SET data= ?, log_id=?, md5=? WHERE record_id=?"
            self.db.execute(cmnd, (data, self._log_id, md5, unique_id))
        else:
            cmnd = f"INSERT INTO {table_name} (record_id,data,log_id,md5,not_completed) VALUES (?,?,?,?,?)"
            self.db.execute(cmnd, (unique_id, data, self._log_id, md5, not_completed))
        return DataMember(self, Path(table_name) / unique_id)

    def drop_not_completed(self) -> None:
        # self.db.execute(f"DELETE FROM {_MD5_TABLE} WHERE record_id IN (SELECT record_id FROM {_RESULT_TABLE} WHERE not_completed=1)")
        self.db.execute(f"DELETE FROM {_RESULT_TABLE} WHERE not_completed=1")
        self._not_completed = []

    def write(self, *, unique_id: str, data: str) -> DataMember:
        if unique_id.startswith(_RESULT_TABLE):
            unique_id = Path(unique_id).name

        super().write(unique_id=unique_id, data=data)
        member = self._write(
            table_name=_RESULT_TABLE,
            unique_id=unique_id,
            data=data,
            not_completed=False,
        )
        if member is not None:
            self._completed.append(member)
        return member

    def write_log(self, *, unique_id: str, data: str) -> None:
        if unique_id.startswith(_LOG_TABLE):
            unique_id = Path(unique_id).name

        super().write_log(unique_id=unique_id, data=data)
        _ = self._write(
            table_name=_LOG_TABLE, unique_id=unique_id, data=data, not_completed=False
        )

    def write_not_completed(self, *, unique_id: str, data: str) -> DataMember:
        if unique_id.startswith(_RESULT_TABLE):
            unique_id = Path(unique_id).name

        super().write_not_completed(unique_id=unique_id, data=data)
        member = self._write(
            table_name=_RESULT_TABLE, unique_id=unique_id, data=data, not_completed=True
        )
        if member is not None:
            self._not_completed.append(member)
        return member

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
        unique_id = re.sub(rf"[.]({self.store_suffix}|json)$", ".txt", unique_id.name)

        cmnd = f"SELECT * FROM {_RESULT_TABLE} WHERE record_id = ?"
        result = self.db.execute(cmnd, (unique_id,)).fetchone()

        return result["md5"]

    @property
    def describe(self):
        if self._lock_id:
            title = f"Locked db store. Locked to pid={self._lock_id}, current pid={os.getpid()}."
        else:
            title = "Unlocked db store."
        table = super().describe
        table.title = title
        return table