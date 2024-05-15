from __future__ import annotations

import collections
import copy
import inspect
import io
import json
import pathlib
import sqlite3
import sys
import typing
import warnings

import numpy
import typing_extensions

import cogent3.util.warning as c3warn

from cogent3._version import __version__
from cogent3.parse.gff import merged_gff_records
from cogent3.util.deserialise import register_deserialiser
from cogent3.util.io import PathType, iter_line_blocks
from cogent3.util.misc import extend_docstring_from, get_object_provenance
from cogent3.util.progress_display import display_wrap
from cogent3.util.table import Table


OptionalInt = typing.Optional[int]
OptionalStr = typing.Optional[str]
OptionalStrList = typing.Optional[typing.Union[str, typing.List[str]]]
OptionalStrContainer = typing.Optional[typing.Union[str, typing.Sequence]]
OptionalBool = typing.Optional[bool]
OptionalDbCursor = typing.Optional[sqlite3.Cursor]
ReturnType = typing.Tuple[str, tuple]  # the sql statement and corresponding values
# data type for sqlitedb constructors
T = typing.Optional[typing.Iterable[dict]]

# used for presence of sqlite feature
_is_ge_3_11 = (sys.version_info.major, sys.version_info.minor) >= (3, 11)


# Define custom types for storage in sqlite
# https://stackoverflow.com/questions/18621513/python-insert-numpy-array-into-sqlite3-database
def array_to_sqlite(data: numpy.ndarray) -> bytes:
    with io.BytesIO() as out:
        numpy.save(out, data)
        out.seek(0)
        output = out.read()
    return output


def sqlite_to_array(data: bytes) -> numpy.ndarray:
    with io.BytesIO(data) as out:
        out.seek(0)
        try:
            result = numpy.load(out)
        except ValueError:
            # array is not stored in the numpy.save format
            # attempt to read from the old format where the
            # array was saved using numpy.ndarray.tobytes
            result = numpy.frombuffer(data, dtype=int)
            dim = result.shape[0] // 2
            result = result.reshape((dim, 2))

            warnings.warn(
                "Old OS dependent database file format detected. "
                "Update the file format using cogent3.core.annotation_db.update_file_format() "
                "For reason see https://github.com/cogent3/cogent3/issues/1776.",
                UserWarning,
            )

    return result


def dict_to_sqlite_as_json(data: dict) -> str:
    return json.dumps(data)


def sqlite_json_to_dict(data):
    return json.loads(data)


# registering the conversion functions with sqlite
sqlite3.register_adapter(numpy.ndarray, array_to_sqlite)
sqlite3.register_converter("array", sqlite_to_array)
sqlite3.register_adapter(dict, dict_to_sqlite_as_json)
sqlite3.register_converter("json", sqlite_json_to_dict)


class FeatureDataType(typing.TypedDict):
    seqid: str  # name of the sseq
    biotype: str  # rename type attr of cogent3 Annotatables to match this?
    name: str  # rename to name to match cogent3 Annotatable.name?
    spans: list[tuple[int, int]]
    strand: str  # "-" if feature on reverse strand
    on_alignment: bool  # True if feature on an alignment


@typing.runtime_checkable
class SerialisableType(typing.Protocol):  # pragma: no cover
    def to_rich_dict(self) -> dict: ...

    def to_json(self) -> str: ...

    def from_dict(self, data): ...


@typing.runtime_checkable
class SupportsQueryFeatures(typing.Protocol):  # pragma: no cover
    # should be defined centrally
    @c3warn.deprecated_args(
        "2024.9",
        reason="replace usage of SQL keyword as column name",
        old_new=[("end", "stop")],
    )
    def get_features_matching(
        self,
        *,
        seqid: OptionalStr = None,
        biotype: OptionalStr = None,
        name: OptionalStr = None,
        start: OptionalInt = None,
        stop: OptionalInt = None,
        strand: OptionalStr = None,
        attributes: OptionalStr = None,
        on_alignment: OptionalBool = None,
    ) -> typing.Iterator[FeatureDataType]: ...

    @c3warn.deprecated_args(
        "2024.9",
        reason="replace usage of SQL keyword as column name",
        old_new=[("end", "stop")],
    )
    def get_feature_children(
        self,
        *,
        name: str,
        start: OptionalInt = None,
        stop: OptionalInt = None,
        **kwargs,
    ) -> typing.List[FeatureDataType]: ...

    @c3warn.deprecated_args(
        "2024.9",
        reason="replace usage of SQL keyword as column name",
        old_new=[("end", "stop")],
    )
    def get_feature_parent(
        self,
        *,
        name: str,
        start: OptionalInt = None,
        stop: OptionalInt = None,
        **kwargs,
    ) -> typing.List[FeatureDataType]: ...

    def num_matches(
        self,
        *,
        seqid: OptionalStr = None,
        biotype: OptionalStr = None,
        name: OptionalStr = None,
        strand: OptionalStr = None,
        attributes: OptionalStr = None,
        on_alignment: OptionalBool = None,
    ) -> int: ...

    @c3warn.deprecated_args(
        "2024.9",
        reason="replace usage of SQL keyword as column name",
        old_new=[("end", "stop")],
    )
    def subset(
        self,
        *,
        seqid: OptionalStr = None,
        biotype: OptionalStr = None,
        name: OptionalStr = None,
        start: OptionalInt = None,
        stop: OptionalInt = None,
        strand: OptionalStr = None,
        attributes: OptionalStr = None,
    ): ...


@typing.runtime_checkable
class SupportsWriteFeatures(typing.Protocol):  # pragma: no cover
    # should be defined centrally
    def add_feature(
        self,
        *,
        biotype: str,
        name: str,
        spans: typing.List[typing.Tuple[int, int]],
        seqid: OptionalStr = None,
        parent_id: OptionalStr = None,
        attributes: OptionalStr = None,
        strand: OptionalStr = None,
        on_alignment: bool = False,
    ) -> None: ...

    def add_records(
        self,
        *,
        records: typing.Sequence[dict],
        # seqid required for genbank
        seqid: OptionalStr = None,
    ) -> None: ...

    def update(self, annot_db, seqids: OptionalStrList = None) -> None:
        # update records with those from an instance of the same type
        ...

    def union(self, annot_db):
        # returns a new instance of the more complex class
        ...


@typing.runtime_checkable
class SupportsFeatures(
    SupportsQueryFeatures,
    SupportsWriteFeatures,
    SerialisableType,
    typing.Sized,
    typing.Protocol,
):  # pragma: no cover
    # should be defined centrally
    def __len__(self):
        # the number of records
        ...

    def __eq__(self, other):
        # equality based on class and identity of the bound db
        ...


def _make_table_sql(
    table_name: str,
    columns: dict,
) -> str:
    """makes the SQL for creating a table

    Parameters
    ----------
    table_name : str
        name of the table
    columns : dict
        {<column name>: <column SQL type>, ...}

    Returns
    -------
    str
    """
    columns_types = ", ".join([f"{name} {ctype}" for name, ctype in columns.items()])
    return f"CREATE TABLE IF NOT EXISTS {table_name} ({columns_types});"


# for making sqlite regex queries see https://stackoverflow.com/questions/67209494/python3-sqlite-and-regex-query


def _add_record_sql(
    table_name: str,
    data: dict,
) -> ReturnType:
    """creates SQL defining the table

    Parameters
    ----------
    table_name : str
        name of the table
    data : dict
        {<column name>: column type}

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """
    cols = ", ".join(k.lower() for k in data)
    pos = ", ".join("?" * len(data))
    sql = f"INSERT INTO {table_name} ({cols}) VALUES ({pos});"
    return sql, tuple(data.values())


def _matching_conditions(
    conditions: dict,
    allow_partial: bool = True,
):
    """creates WHERE clause

    Parameters
    ----------
    conditions : dict
        column name and values to be matched
    allow_partial : bool, optional
        if False, only records within start, stop are included. If True,
        all records that overlap the segment defined by start, stop are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """
    start = conditions.pop("start", None)
    stop = conditions.pop("stop", None)

    sql = []
    vals = ()
    if conditions:
        conds = []
        vals = []
        for col, val in conditions.items():
            # conditions are filtered for None before here, so we should add
            # an else where the op is assigned !=
            if isinstance(val, (tuple, set, list)):
                placeholders = ",".join(["?" for _ in val])
                op = "IN"
                conds.append(f"{col} {op} ({placeholders})")
                vals.extend(val)
            elif val is not None:
                op = "LIKE" if isinstance(val, str) and "%" in val else "="
                conds.append(f"{col} {op} ?")
                vals.append(val)
        sql.append(" AND ".join(conds))
        vals = tuple(vals)

    if start is not None and stop is not None:
        if allow_partial:
            # allow matches that overlap the segment
            cond = [
                f"(start >= {start} AND stop <= {stop})",  # lies within the segment
                f"(start <= {start} AND stop > {start})",  # straddles beginning of segment
                f"(start < {stop} AND stop >= {stop})",  # straddles stop of segment
                f"(start <= {start} AND stop >= {stop})",  # includes segment
            ]
            cond = " OR ".join(cond)
        else:
            # only matches within bounds
            cond = f"start >= {start} AND stop <= {stop}"
        sql.append(f"({cond})")
    elif start is not None:
        # if query has no stop, then any feature containing start
        cond = f"(start <= {start} AND {start} < stop)"
        sql.append(f"({cond})")
    elif stop is not None:
        # if query has no start, then any feature containing stop
        cond = f"(start <= {stop} AND {stop} < stop)"
        sql.append(f"({cond})")

    sql = f"{' AND '.join(sql)}"
    return sql, vals


def _del_records_sql(
    table_name: str,
    conditions: dict,
    start: OptionalInt = None,
    stop: OptionalInt = None,
    allow_partial=True,
) -> ReturnType:
    """creates the SQL and values for identifying records to be deleted

    Parameters
    ----------
    table_name : str
        table to have records deleted from
    conditions : dict
        column name and values to be matched
    start, stop : OptionalInt
        select records whose (start, stop) values lie between start and stop,
        or overlap them if (allow_partial is True)
    allow_partial : bool, optional
        if False, only records within start, stop are included. If True,
        all records that overlap the segment defined by start, stop are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """
    where, vals = _matching_conditions(
        conditions=conditions, start=start, stop=stop, allow_partial=allow_partial
    )
    sql = f"DELETE FROM {table_name}"
    if not where:
        return sql

    sql = f"{sql} WHERE {where};"
    return sql, vals


def _select_records_sql(
    table_name: str,
    conditions: dict,
    columns: typing.Optional[typing.List[str]] = None,
    start: OptionalInt = None,
    stop: OptionalInt = None,
    allow_partial=True,
) -> ReturnType:
    """create SQL select statement and values

    Parameters
    ----------
    table_name : str
        containing the data to be selected from
    columns : Tuple[str]
        values to select
    conditions : dict
        the WHERE conditions
    start, stop : OptionalInt
        select records whose (start, stop) values lie between start and stop,
        or overlap them if (allow_partial is True)
    allow_partial : bool, optional
        if False, only records within start, stop are included. If True,
        all records that overlap the segment defined by start, stop are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """

    where, vals = _matching_conditions(
        conditions=conditions, allow_partial=allow_partial
    )
    columns = f"{', '.join(columns)}" if columns else "*"
    sql = f"SELECT {columns} FROM {table_name}"
    if not where:
        return sql, None

    sql = f"{sql} WHERE {where};"
    return sql, vals


def _count_records_sql(
    table_name: str,
    conditions: dict,
    columns: typing.Optional[typing.List[str]] = None,
    start: OptionalInt = None,
    stop: OptionalInt = None,
    allow_partial=True,
) -> ReturnType:
    """create SQL count statement and values

    Parameters
    ----------
    table_name : str
        containing the data to be selected from
    columns : Tuple[str]
        values to select
    conditions : dict
        the WHERE conditions
    start, stop : OptionalInt
        select records whose (start, stop) values lie between start and stop,
        or overlap them if (allow_partial is True)
    allow_partial : bool, optional
        if False, only records within start, stop are included. If True,
        all records that overlap the segment defined by start, stop are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """

    where, vals = _matching_conditions(
        conditions=conditions, allow_partial=allow_partial
    )
    sql = f"SELECT COUNT(*) FROM {table_name}"
    if not where:
        return sql, None

    sql = f"{sql} WHERE {where};"
    return sql, vals


def _compatible_schema(db: sqlite3.Connection, schema: dict[str, set]) -> bool:
    """ensures the db instance is compatible with schema"""
    for table in db.execute(
        "SELECT name FROM sqlite_master WHERE type='table'"
    ).fetchall():
        table = table["name"]
        if table not in schema:
            return False

        db_schema = {
            (row["name"], row["type"])
            for row in db.execute(f"pragma table_info({table!r})").fetchall()
        }

        if len(db_schema & schema[table]) != len(db_schema):
            return False

    return True


def _rename_column_if_exists(
    db: sqlite3.Connection, table_name: str, old_column: str, new_column: str
) -> None:
    """Rename a column in a sqlite3 database only if it exists.

    Parameters
    ----------
    db : sqlite3.Connection
        the sqlite3 connection
    table_name : str
        the table name
    old_column : str
        the column to rename if it exists
    new_column : str
        the new name of old_column
    """
    cur = db.cursor()
    table_info = cur.execute(f"PRAGMA table_info({table_name});").fetchall()

    # Check if the column exists
    for column_info in table_info:
        if column_info["name"] == old_column:
            break
    else:
        return  # There is no column to rename

    sql = f'ALTER TABLE {table_name} RENAME COLUMN "{old_column}" TO {new_column}'

    cur.execute(sql)
    db.commit()


class SqliteAnnotationDbMixin:
    # table schema for user provided annotations
    _user_schema = {
        "biotype": "TEXT",
        "seqid": "TEXT",
        "name": "TEXT",
        "parent_id": "TEXT",
        "start": "INTEGER",
        "stop": "INTEGER",
        "strand": "TEXT",
        "spans": "array",
        "attributes": "TEXT",
        "on_alignment": "INT",
    }
    # args to exclude from serialisation init
    _exclude_init = "db", "data"

    def __new__(cls, *args, **kwargs):
        obj = object.__new__(cls)
        init_sig = inspect.signature(cls.__init__)
        bargs = init_sig.bind_partial(cls, *args, **kwargs)
        bargs.apply_defaults()
        init_vals = bargs.arguments
        for param in cls._exclude_init:
            init_vals.pop(param)
        init_vals.pop("self", None)

        obj._serialisable = init_vals
        return obj

    def __repr__(self):
        name = self.__class__.__name__
        total_records = len(self)
        args = ", ".join(
            f"{k}={repr(v) if isinstance(v, str) else v}"
            for k, v in self._serialisable.items()
            if k != "data"
        )
        return f"{name}({args}, total_records={total_records})"

    def __deepcopy__(self, memodict=None):
        memodict = memodict or {}
        if _is_ge_3_11:
            new = self.__class__(source=self.source)
            new.db.deserialize(self._db.serialize())
            return new

        # use rich dict
        rd = self.to_rich_dict()
        return type(self).from_dict(rd)

    def __getstate__(self):
        if _is_ge_3_11:
            return {"data": self._db.serialize(), "source": self.source}

        return self.to_rich_dict()

    def __setstate__(self, state):
        if _is_ge_3_11:
            new = self.__class__(source=state.pop("source", None))
            new._db.deserialize(state["data"])
            self.__dict__.update(new.__dict__)
            return self

        # from the rich dict method
        data = type(self).from_dict(state)
        self.__dict__.update(data.__dict__)
        return self

    def __len__(self):
        return self.num_matches()

    def __eq__(self, other):
        return isinstance(other, self.__class__) and other.db is self.db

    @property
    def table_names(self) -> tuple[str]:
        return self._table_names

    def _setup_db(self, db: typing.Union[SupportsFeatures, sqlite3.Connection]) -> None:
        """initialises the db, using the db passed to the constructor"""
        if isinstance(db, self.__class__):
            self._db = db.db
            return

        if isinstance(db, sqlite3.Connection):
            schema = {}
            for table_name in self.table_names:
                attr = getattr(self, f"_{table_name}_schema")
                schema[table_name] = set(attr.items())

            if not _compatible_schema(db, schema):
                raise TypeError("incompatible schema")

            self._db = db
            self._init_tables()
            return

        if db and not self.compatible(db):
            raise TypeError(f"cannot initialise annotation db from {type(db)}")

        self._init_tables()

        if db and len(db):
            # update self with data from other
            self.update(db)

    def _init_tables(self) -> None:
        # bit of magic
        # assumes schema attributes named as `_<table name>_schema`
        for table_name in self.table_names:
            attr = getattr(self, f"_{table_name}_schema")
            sql = _make_table_sql(table_name, attr)
            self._execute_sql(sql)

    @property
    def db(self) -> sqlite3.Connection:
        if self._db is None:
            # todo gah understand serialisation issue
            # the check_same_thread=False is required for multiprocess, even
            # when the db is empty (tests fail). This  appears unrelated to
            # our code, and does not affect pickling/unpickling on the same
            # thread
            self._db = sqlite3.connect(
                self.source,
                detect_types=sqlite3.PARSE_DECLTYPES,
                check_same_thread=False,
            )
            self._db.row_factory = sqlite3.Row

            for table_name in self.table_names:
                _rename_column_if_exists(self._db, table_name, "end", "stop")

        return self._db

    def _execute_sql(self, cmnd, values=None):
        with self.db:
            # context manager ensures safe transactions
            cursor = self.db.cursor()
            cursor.execute(cmnd, values or [])
            return cursor

    def add_feature(
        self,
        *,
        seqid: str,
        biotype: str,
        name: str,
        spans: typing.List[typing.Tuple[int, int]],
        parent_id: OptionalStr = None,
        strand: OptionalStr = None,
        attributes: OptionalStr = None,
        on_alignment: OptionalBool = False,
    ) -> None:
        """adds a record to user table

        Parameters
        ----------
        seqid : str
            name of the sequence feature resides on
        biotype : str
            biological type of the record
        name : str
            the name of a record, an identifier
        spans : typing.List[typing.Tuple[int, int]], optional
            this will be sorted
        strand : str, optional
            either +, -. Defaults to '+'
        attributes : str, optional
            additional attributes as a string
        on_alignment : bool, optional
            whether the annotation is an alignment annotation
        """
        spans = numpy.array(sorted(sorted(coords) for coords in spans), dtype=int)
        # the start, stop variables are added into the record via the loop
        # on local variables
        start = int(spans.min())
        stop = int(spans.max())
        # we define record as all defined variables from local name space,
        # excluding "self"
        record = {k: v for k, v in locals().items() if k != "self" and v is not None}
        sql, values = _add_record_sql("user", record)
        self._execute_sql(sql, values=values)

    def _get_records_matching(
        self, table_name: str, **kwargs
    ) -> typing.Iterator[sqlite3.Row]:
        """return all fields"""
        if kwargs.get("attributes", None) and "%%" not in kwargs["attributes"]:
            kwargs["attributes"] = f'%{kwargs["attributes"]}%'
        columns = kwargs.pop("columns", None)
        allow_partial = kwargs.pop("allow_partial", False)
        sql, vals = _select_records_sql(
            table_name=table_name,
            conditions=kwargs,
            columns=columns,
            allow_partial=allow_partial,
        )
        yield from self._execute_sql(sql, values=vals)

    def _get_feature_by_id(
        self,
        table_name: str,
        columns: typing.Optional[typing.Iterable[str]],
        column: str,
        name: str,
        start: OptionalInt = None,
        stop: OptionalInt = None,
        biotype: OptionalStr = None,
        allow_partial: bool = False,
    ) -> typing.List[FeatureDataType]:
        # we return the parent_id because `get_feature_parent()` requires it
        sql, vals = _select_records_sql(
            table_name=table_name,
            conditions={column: name, "biotype": biotype},
            columns=columns,
            start=start,
            stop=stop,
            allow_partial=allow_partial,
        )
        for result in self._execute_sql(sql, values=vals):
            result = dict(zip(result.keys(), result))
            result["on_alignment"] = result.get("on_alignment")
            result["spans"] = [tuple(c) for c in result["spans"]]
            yield result

    def get_feature_children(
        self, name: str, biotype: OptionalStr = None, **kwargs
    ) -> typing.Iterator[FeatureDataType]:
        """yields children of name"""
        # kwargs is used because other classes need start / stop
        # just uses search matching
        for table_name in self.table_names:
            columns = "seqid", "biotype", "spans", "strand", "name", "parent_id"
            if table_name == "user":
                columns += ("on_alignment",)
            for result in self._get_feature_by_id(
                table_name=table_name,
                columns=columns,
                column="parent_id",
                name=f"%{name}%",
                biotype=biotype,
                **kwargs,
            ):
                result.pop("parent_id")  # remove invalid field for the FeatureDataType
                yield result

    def get_feature_parent(
        self, name: str, **kwargs
    ) -> typing.Iterator[FeatureDataType]:
        """yields parents of name"""
        for table_name in self.table_names:
            columns = "seqid", "biotype", "spans", "strand", "name", "parent_id"
            if table_name == "user":
                columns += ("on_alignment",)
            for result in self._get_feature_by_id(
                table_name=table_name,
                columns=columns,
                column="name",
                name=f"%{name}%",
            ):
                # multiple values for parent means this is better expressed
                # as an OR clause
                # todo modify the conditional SQL generation
                if not result["parent_id"]:
                    return

                for name in result["parent_id"].replace(" ", "").split(","):
                    if parent := list(
                        self._get_feature_by_id(
                            table_name=table_name,
                            columns=columns,
                            column="name",
                            name=name,
                        )
                    ):
                        parent = parent[0]
                        parent.pop(
                            "parent_id"
                        )  # remove invalid field for the FeatureDataType
                        yield parent

    @c3warn.deprecated_args(
        "2024.9",
        reason="replace usage of SQL keyword as column name",
        old_new=[("end", "stop")],
    )
    def get_records_matching(
        self,
        *,
        biotype: str = None,
        seqid: str = None,
        name: str = None,
        start: int = None,
        stop: int = None,
        strand: OptionalStr = None,
        attributes: OptionalStr = None,
        on_alignment: bool = None,
        allow_partial: bool = False,
    ) -> typing.Iterator[dict]:
        """return all fields for matching records"""
        # a record is Everything, a Feature is a subset
        # we define query as all defined variables from local name space,
        # excluding "self" and kwargs at default values
        kwargs = {k: v for k, v in locals().items() if k != "self" and v is not None}
        # alignment features are created by the user specific
        table_names = ["user"] if on_alignment else self.table_names
        for table_name in table_names:
            for result in self._get_records_matching(table_name, **kwargs):
                yield dict(zip(result.keys(), result))

    @c3warn.deprecated_args(
        "2024.9",
        reason="replace usage of SQL keyword as column name",
        old_new=[("end", "stop")],
    )
    def get_features_matching(
        self,
        *,
        biotype: str = None,
        seqid: str = None,
        name: str = None,
        start: int = None,
        stop: int = None,
        strand: OptionalStr = None,
        attributes: OptionalStr = None,
        on_alignment: bool = None,
        allow_partial: bool = False,
    ) -> typing.Iterator[FeatureDataType]:
        # returns essential values to create a Feature
        # we define query as all defined variables from local name space,
        # excluding "self" and kwargs with a default value of None
        kwargs = {k: v for k, v in locals().items() if k != "self" and v is not None}
        # alignment features are created by the user specific
        table_names = ["user"] if on_alignment else self.table_names
        for table_name in table_names:
            columns = ("seqid", "biotype", "spans", "strand", "name")
            query_args = {**kwargs}

            if table_name == "user":
                columns += ("on_alignment",)
            else:
                query_args.pop("on_alignment", None)

            for result in self._get_records_matching(
                table_name=table_name, columns=columns, **query_args
            ):
                result = dict(zip(result.keys(), result))
                result["on_alignment"] = result.get("on_alignment")
                result["spans"] = [tuple(c) for c in result["spans"]]
                yield result

    def num_matches(
        self,
        *,
        seqid: OptionalStr = None,
        biotype: OptionalStr = None,
        name: OptionalStr = None,
        strand: OptionalStr = None,
        attributes: OptionalStr = None,
        on_alignment: OptionalBool = None,
    ) -> int:
        """return the number of records matching condition"""
        kwargs = {k: v for k, v in locals().items() if k != "self"}
        num = 0
        for table_name in self.table_names:
            sql, values = _count_records_sql(table_name, conditions=kwargs)
            result = list(self._execute_sql(sql, values=values).fetchone())[0]
            num += result
        return num

    StrOrBool = typing.Union[str, bool]

    def count_distinct(
        self,
        *,
        seqid: StrOrBool = False,
        biotype: StrOrBool = False,
        name: StrOrBool = False,
    ) -> typing.Optional[Table]:
        """return table of counts of distinct values

        Parameters
        ----------
        seqid, biotype, name
            if a string, selects the subset of rows matching the provided values
            and counts distinct values for the other fields whose value is True.

        Returns
        -------
        Table with columns corresponding to argument whose value was True

        Examples
        --------
        To compute copy number by gene name within each genome

        >>> counts_table = db.count_distinct(seqid=True, biotype='gene', name=True)
        """
        columns = {k for k, v in locals().items() if v is True}
        if not columns:
            return

        if constraints := {k: v for k, v in locals().items() if isinstance(v, str)}:
            where_clause, values = _matching_conditions(constraints)
            where_clause = f"WHERE {where_clause}"
        else:
            where_clause = ""
            values = ()

        header = list(columns)
        placeholder = "{}"
        sql_template = (
            f"SELECT {', '.join(header)}, COUNT(*) as count FROM {placeholder}"
            f" {where_clause} GROUP BY {', '.join(header)};"
        )
        data = []
        for table in self.table_names:
            sql = sql_template.format(table)
            if result := self._execute_sql(sql, values).fetchall():
                data.extend(tuple(r) for r in result)

        return Table(
            header=header + ["count"],
            data=data,
            column_templates=dict(count=lambda x: f"{x:,}"),
        )

    @property
    def describe(self) -> Table:
        """top level description of the annotation db"""
        sql_template = "SELECT {}, COUNT(*) FROM {} GROUP BY {};"
        data = {}
        for column in ("seqid", "biotype"):
            for table in self.table_names:
                sql = sql_template.format(column, table, column)
                if result := self._execute_sql(sql).fetchall():
                    counts = dict(tuple(r) for r in result)
                    for distinct, count in counts.items():
                        key = f"{column}({distinct!r})"
                        data[key] = data.get(key, 0) + count

        row_counts = []
        for table in self.table_names:
            result = self._execute_sql(
                f"SELECT COUNT(*) as count FROM {table}"
            ).fetchone()
            row_counts.append(result["count"])

        return Table(
            data={
                "": list(data.keys()) + [f"num_rows({t!r})" for t in self.table_names],
                "count": list(data.values()) + row_counts,
            },
            column_templates=dict(count=lambda x: f"{x:,}"),
        )

    def biotype_counts(self) -> dict:
        """return counts of biological types across all tables and seqids"""
        sql_template = "SELECT biotype FROM {}"
        counts = collections.Counter()
        for table in self.table_names:
            sql = sql_template.format(table)
            if result := self._execute_sql(sql).fetchall():
                counts.update(v for r in result if (v := r["biotype"]))
        return counts

    def to_rich_dict(self) -> dict:
        """returns a dict suitable for json serialisation"""
        result = {
            "type": get_object_provenance(self),
            "version": __version__,
            "tables": {},
            "init_args": {**self._serialisable},
        }
        tables = result["tables"]
        for table_name in self.table_names:
            table_data = []
            for record in self._get_records_matching(table_name):
                store = {k: v for k, v in zip(record.keys(), record) if v is not None}
                store["spans"] = store["spans"].tolist()
                table_data.append(store)
            tables[table_name] = table_data
        return result

    def compatible(self, other_db: SupportsFeatures, symmetric=True) -> bool:
        """checks whether table_names are compatible

        Parameters
        ----------
        other_db
            the other annotation db instance
        symmetric
            checks only that tables of other_db equal, or are a subset, of
            mine
        """
        if not isinstance(other_db, SupportsFeatures):
            raise TypeError(f"{type(other_db)} does not support features")
        mine = set(self.table_names)
        theirs = set(other_db.table_names)
        return mine <= theirs or mine > theirs if symmetric else mine >= theirs

    def update(
        self, annot_db: SupportsFeatures, seqids: OptionalStrList = None
    ) -> None:
        """update records with those from an instance of the same type"""
        if not isinstance(annot_db, SupportsFeatures):
            raise TypeError(f"{type(annot_db)} does not satisfy SupportsFeatures")
        elif not self.compatible(annot_db, symmetric=False):
            raise TypeError(f"{type(self)} cannot be updated from {type(annot_db)}")

        if not annot_db or not len(annot_db):
            return

        self._update_db_from_other_db(annot_db, seqids=seqids)

    def union(self, annot_db: SupportsFeatures) -> SupportsFeatures:
        """returns a new instance with merged records with other

        Parameters
        ----------
        annot_db
            an annotation db whose schema is either a subset, or superset of
            self

        Returns
        -------
        The class whose schema contains the other
        """
        if annot_db and not isinstance(annot_db, SupportsFeatures):
            raise TypeError(f"{type(annot_db)} does not satisfy SupportsFeatures")
        elif not annot_db:
            return copy.deepcopy(self)

        if self.compatible(annot_db, symmetric=False):
            cls = type(self)
        elif self.compatible(annot_db, symmetric=True):
            cls = type(annot_db)
        else:
            raise TypeError(
                f"cannot make a union between {type(self)} and {type(annot_db)}"
            )

        db = cls()
        db.update(self)
        db.update(annot_db)
        return db

    @classmethod
    def from_dict(cls, data: dict):
        # make an empty db
        init_args = data.pop("init_args")
        db = cls(**init_args)
        db._update_db_from_rich_dict(data)
        return db

    def _update_db_from_rich_dict(self, data: dict, seqids: OptionalStr = None):
        data.pop("type", None)
        data.pop("version")
        if isinstance(seqids, str):
            seqids = {seqids}
        elif seqids is not None:
            seqids = set(seqids) - {None}  # make sure None is not part of this!

        # todo gah prevent duplication of existing records
        for table_name, records in data["tables"].items():
            for record in records:
                if seqids and record["seqid"] not in seqids:
                    continue
                record["spans"] = numpy.array(record["spans"], dtype=int)
                sql, vals = _add_record_sql(
                    table_name, {k: v for k, v in record.items() if v is not None}
                )
                self._execute_sql(sql, vals)

    def _update_db_from_other_db(
        self, other_db: SupportsFeatures, seqids: OptionalStrContainer = None
    ):
        conditions = {"seqid": seqids} if seqids else {}
        table_names = other_db.table_names

        col_order = {
            tname: [
                row[1]
                for row in other_db.db.execute(f"PRAGMA table_info({tname})").fetchall()
            ]
            for tname in table_names
        }

        for tname in other_db.table_names:
            sql, vals = _select_records_sql(table_name=tname, conditions=conditions)
            data = other_db._execute_sql(sql, vals)
            val_placeholder = ", ".join("?" * len(col_order[tname]))
            sql = f"INSERT INTO {tname} ({', '.join(col_order[tname])}) VALUES ({val_placeholder})"

            self.db.executemany(sql, data)

    def to_json(self) -> str:
        return json.dumps(self.to_rich_dict())

    def write(self, path: PathType) -> None:
        """writes db as bytes to path"""
        backup = sqlite3.connect(path)
        with self.db:
            self.db.backup(backup)
        backup.close()

    @c3warn.deprecated_args(
        "2024.9",
        reason="replace usage of SQL keyword as column name",
        old_new=[("end", "stop")],
    )
    def subset(
        self,
        *,
        source: PathType = ":memory:",
        biotype: str = None,
        seqid: str = None,
        name: str = None,
        start: OptionalInt = None,
        stop: OptionalInt = None,
        strand: OptionalStr = None,
        attributes: OptionalStr = None,
        allow_partial: bool = False,
    ) -> typing_extensions.Self:
        """returns a new db instance with records matching the provided conditions"""
        # make sure python, not numpy, integers
        start = start if start is None else int(start)
        stop = stop if stop is None else int(stop)

        kwargs = {k: v for k, v in locals().items() if k not in {"self", "source"}}

        result = self.__class__(source=source)
        if not len(self):
            return result

        for table_name in self.table_names:
            cols = [
                r["name"]
                for r in self.db.execute(f"PRAGMA table_info({table_name})").fetchall()
            ]
            pos = ", ".join("?" * len(cols))
            sql = f"INSERT INTO {table_name} ({','.join(cols)}) VALUES ({pos});"
            records = []
            for record in self._get_records_matching(table_name=table_name, **kwargs):
                records.append(tuple(record[c] for c in cols))

            with result.db as cursor:
                cursor.executemany(sql, records)
                cursor.commit()

        return result

    def close(self):
        """closes the db"""
        self.db.close()

    def _make_index(self, *, table_name: str, col_names: tuple[str, ...]):
        """index columns for faster search"""
        sql = f"CREATE INDEX IF NOT EXISTS %s on {table_name}(%s)"
        for col in col_names:
            self._execute_sql(sql % (col, col))

    def make_indexes(self):
        """adds db indexes for core attributes"""
        for table_name in self.table_names:
            self._make_index(
                table_name=table_name,
                col_names=("biotype", "seqid", "name", "start", "stop", "parent_id"),
            )


class BasicAnnotationDb(SqliteAnnotationDbMixin):
    """Provides a user table for annotations. This can be merged with
    either the Gff or Genbank versions.

    Notes
    -----
    This is the default db on Sequence, SequenceCollection and Alignment
    """

    _table_names = ("user",)

    def __init__(
        self, *, data: T = None, db: SupportsFeatures = None, source=":memory:"
    ):
        """
        Parameters
        ----------
        data
            data for entry into the database
        db
            a compatible annotation db instance. If db is the same class, it's
            db will be bound to self and directly modified.
        source
            location to store the db, defaults to in memory only
        """
        data = data or []
        # note that data is destroyed
        self._num_fakeids = 0
        self.source = source
        self._db = None
        self._setup_db(db)

        self.add_records(data)

    def add_records(self, data: T) -> None:
        table_name = self.table_names[0]  # only one name for this class
        for record in data:
            record["spans"] = numpy.array(record["spans"], dtype=int)
            cmnd, vals = _add_record_sql(
                table_name,
                {
                    k: v
                    for k, v in record.items()
                    if isinstance(v, numpy.ndarray) or v not in (".", None)
                },
            )
            self._execute_sql(cmnd=cmnd, values=vals)


def _merge_spans(old: numpy.ndarray, new: list[list[int]]) -> numpy.ndarray:
    """returns sorted, merged old and new spans"""
    if len(new) == old.shape[0] and (old == new).all():
        return old

    new = numpy.array(sorted(new), dtype=old.dtype)
    return numpy.unique(numpy.concatenate([old, new]), axis=0)


class GffAnnotationDb(SqliteAnnotationDbMixin):
    """Support for annotations from gff files. Records that span multiple
    rows in the gff are merged into a single record."""

    _table_names = "gff", "user"
    # We are relying on an attribute name structured as _<table name>_schema
    _gff_schema = {
        "seqid": "TEXT",
        "source": "TEXT",
        "biotype": "TEXT",  # type in GFF
        "start": "INTEGER",
        "stop": "INTEGER",
        "score": "TEXT",  # check defn
        "strand": "TEXT",
        "phase": "TEXT",
        "attributes": "TEXT",
        "comments": "TEXT",
        "spans": "array",  # aggregation of coords across records
        "name": "TEXT",
        "parent_id": "TEXT",
    }

    @extend_docstring_from(BasicAnnotationDb.__init__)
    def __init__(
        self, *, data: T = None, db: SupportsFeatures = None, source=":memory:"
    ):
        data = data or []
        # note that data is destroyed
        self._num_fakeids = 0
        # if a db instance is passed in, we use that for source
        self.source = getattr(db, "source", source)
        # ensure serialisable state reflects this
        self._serialisable["source"] = self.source
        self._db = None
        self._setup_db(db)
        data, self._num_fakeids = merged_gff_records(data, self._num_fakeids)
        self.add_records(data)

    def add_records(self, reduced: dict) -> None:
        col_order = [
            r["name"] for r in self.db.execute("PRAGMA table_info(gff)").fetchall()
        ]
        val_placeholder = ", ".join("?" * len(col_order))
        sql = f"INSERT INTO gff ({', '.join(col_order)}) VALUES ({val_placeholder})"

        # Can we really trust text case of "ID" and "Parent" to be consistent
        # across sources of gff?? I doubt it, so more robust regex likely
        # required
        rows = []
        for record in reduced.values():
            # our Feature code assumes start always < stop,
            # we record direction using Strand
            spans = numpy.array(sorted(record["spans"]), dtype=int)  # sorts the rows
            spans.sort(axis=1)
            record["start"] = int(spans.min())
            record["stop"] = int(spans.max())
            record["spans"] = spans
            rows.append(tuple(record.get(c) for c in col_order))

        self.db.executemany(sql, rows)

        self.db.commit()
        del reduced

    def update_record_spans(self, *, name: str, spans: list[list[int]]) -> None:
        """updates spans attribute of a gff table record if present

        Notes
        -----
        Has no effect if name is not present.
        """
        if not len(spans):
            return

        result = self._execute_sql(
            cmnd="SELECT spans from gff WHERE name = ?", values=(name,)
        ).fetchone()

        if result is None:
            return

        old_spans = _merge_spans(result["spans"], spans)
        self._execute_sql(
            cmnd="UPDATE gff SET spans = ? WHERE name = ?", values=(old_spans, name)
        )


# The GenBank format is less clear on the relationship between identifiers,
# e.g. it can be gene -> SRP_RNA -> exon, gene -> CDS, etc... Without having
# the explicit type-hierarchy from NCBI (which I have not yet found) the only way
# to establish the Parent of a feature is by knowing what other
# features have the same ID and then which of those are candidates to contain
# a feature based on the hierarchy of relationships.
# In other words, we assume that children have the same ID as their parent BUT
# their span lies within the parents. Conversely, for a parent query, the parent
# has the same ID as their child but their span contains the childs.


class GenbankAnnotationDb(SqliteAnnotationDbMixin):
    """Support for annotations from Genbank files.

    Notes
    -----
    Extended attributes are stored as json in the gb, attributes column.
    """

    _table_names = "gb", "user"
    # We are relying on an attribute name structured as _<table name>_schema
    _gb_schema = {
        "seqid": "TEXT",
        "source": "TEXT",
        "biotype": "TEXT",  # type in GFF
        "start": "INTEGER",
        "stop": "INTEGER",
        "strand": "TEXT",
        "comments": "TEXT",
        "spans": "array",  # aggregation of coords across records
        "name": "TEXT",
        "symbol": "TEXT",
        "parent_id": "TEXT",
        "attributes": "json",
    }

    @extend_docstring_from(BasicAnnotationDb.__init__)
    def __init__(
        self,
        *,
        data: T = None,
        seqid: OptionalStr = None,
        db: SupportsFeatures = None,
        source=":memory:",
        namer: typing.Callable = None,
    ):
        """
        seqid
            name of the sequence data is associated with
        namer
            callable that takes a record dict and returns a
            [name]
        """
        data = data or []
        # note that data is destroyed
        self._num_fakeids = 0
        self.source = source
        self._db = None
        self._setup_db(db)
        if db:
            # guessing there's multiple seqid's
            self._serialisable["seqid"] = "<multiple seqids>"

        self._namer = namer if callable(namer) else self._default_namer
        self.add_records(data, seqid)

    def add_records(self, records, seqid):
        col_order = [
            r["name"] for r in self.db.execute("PRAGMA table_info(gb)").fetchall()
        ]
        val_placeholder = ", ".join("?" * len(col_order))
        sql = f"INSERT INTO gb ({', '.join(col_order)}) VALUES ({val_placeholder})"

        # need to capture genetic code from genbank
        rows = []
        exclude = {"translation", "location", "type"}  # location is grabbed directly
        for record in records:
            # our Feature code assumes start always < stop,
            store = {"seqid": seqid}
            # we create the location data directly
            if location := record.get("location", None):
                store["spans"] = numpy.array(location.get_coordinates(), dtype=int)
                if strand := location.strand:
                    store["strand"] = "-" if strand == -1 else "+"
                store["start"] = int(store["spans"].min())
                store["stop"] = int(store["spans"].max())

            attrs_keys = record.keys() - exclude
            store["biotype"] = record["type"]
            name = self._namer(record)
            if name is None:
                name = self._make_fake_id(record)
            store["attributes"] = {k: record[k] for k in attrs_keys}
            store["name"] = ",".join(name)
            rows.append(tuple(store.get(c) for c in col_order))

        self.db.executemany(sql, rows)
        self.db.commit()
        del records

    def _default_namer(self, record: dict) -> typing.Union[typing.List[str], None]:
        # we evaluate potential tokens in the genbank record in order of
        # preference for naming. If none of these are found, a fake name
        # will be generated
        for key in (
            "gene",
            "locus_tag",
            "strain",
            "rpt_unit_seq",
            "db_xref",
            "bound_moiety",
            "regulatory_class",
        ):
            if key in record:
                return record[key]

        if element := record.get("mobile_element_type"):
            return element[0].split(":")[-1:]

        return None

    def _make_fake_id(self, record: dict) -> typing.List[str]:
        name = [f"{record.get('type', 'fakeid')}-{self._num_fakeids}"]
        self._num_fakeids += 1
        return name

    @c3warn.deprecated_args(
        "2024.9",
        reason="replace usage of SQL keyword as column name",
        old_new=[("end", "stop")],
    )
    def get_feature_children(
        self,
        name: str,
        biotype: OptionalStr = None,
        exclude_biotype: OptionalStr = None,
        start: OptionalInt = None,
        stop: OptionalInt = None,
    ) -> typing.Iterator[FeatureDataType]:
        """yields children of name"""
        # children must lie within the parent coordinates
        if start is None or stop is None:
            name = self.__class__.__name__
            msg = f"coordinates required to query children for {name!r}"
            raise ValueError(msg)

        for table_name in self.table_names:
            columns = (
                "seqid",
                "biotype",
                "spans",
                "strand",
                "name",
                "parent_id",
                "start",
                "stop",
            )
            if table_name == "user":
                columns += ("on_alignment",)

            for feat in self._get_feature_by_id(
                table_name=table_name,
                columns=columns,
                column="name",
                name=name,
                biotype=biotype,
                start=start,
                stop=stop,
                allow_partial=False,
            ):
                if feat["biotype"] == exclude_biotype:
                    continue
                cstart, cstop = feat.pop("start"), feat.pop("stop")
                if not (start <= cstart < stop and start < cstop <= stop):
                    continue

                feat.pop("parent_id")  # remove invalid field for the FeatureDataType
                yield feat

    @c3warn.deprecated_args(
        "2024.9",
        reason="replace usage of SQL keyword as column name",
        old_new=[("end", "stop")],
    )
    def get_feature_parent(
        self,
        name: str,
        exclude_biotype: OptionalStr = None,
        start: OptionalInt = None,
        stop: OptionalInt = None,
    ) -> typing.Iterator[FeatureDataType]:
        """yields parents of name"""
        if start is None or stop is None:
            name = self.__class__.__name__
            msg = f"coordinates required to query parent for {name!r}"
            raise ValueError(msg)

        # parent must match or lie outside the coordinates
        for table_name in self.table_names:
            columns = (
                "seqid",
                "biotype",
                "spans",
                "strand",
                "name",
                "parent_id",
                "start",
                "stop",
            )
            if table_name == "user":
                columns += ("on_alignment",)
            for feat in self._get_feature_by_id(
                table_name=table_name,
                columns=columns,
                column="name",
                name=name,
                start=start,
                stop=stop,
                allow_partial=False,
            ):
                # add support for != operation to SQL where clause generation
                if feat["biotype"] == exclude_biotype:
                    continue

                cstart, cstop = feat.pop("start"), feat.pop("stop")
                if cstart > start or stop > cstop:
                    continue

                feat.pop("parent_id")  # remove invalid field for the FeatureDataType
                yield feat


@register_deserialiser(get_object_provenance(BasicAnnotationDb))
def deserialise_basic_db(data: dict):
    return BasicAnnotationDb.from_dict(data)


@register_deserialiser(get_object_provenance(GffAnnotationDb))
def deserialise_gff_db(data: dict):
    return GffAnnotationDb.from_dict(data)


@register_deserialiser(get_object_provenance(GenbankAnnotationDb))
def deserialise_gb_db(data: dict):
    return GenbankAnnotationDb.from_dict(data)


@register_deserialiser("annotation_to_annotation_db")
def convert_annotation_to_annotation_db(data: dict) -> SupportsFeatures:
    from cogent3.util.deserialise import deserialise_map_spans

    db = BasicAnnotationDb()

    seqid = data.pop("name", data.pop("seqid", None))
    anns = data.pop("data")
    for ann in anns:
        ann = ann.pop("annotation_construction")
        m = deserialise_map_spans(ann.pop("map"))
        spans = m.get_coordinates()
        strand = "-" if m.reverse else "+"
        biotype = ann.pop("type")
        name = ann.pop("name")
        db.add_feature(
            seqid=seqid, biotype=biotype, name=name, spans=spans, strand=strand
        )

    return db


@display_wrap
def _db_from_genbank(
    path: PathType, db: typing.Optional[SupportsFeatures], write_path, **kwargs
):
    from cogent3 import open_
    from cogent3.parse.genbank import MinimalGenbankParser

    paths = pathlib.Path(path)
    paths = list(paths.parent.glob(paths.name))

    ui = kwargs.pop("ui")
    one_valid_path = False
    for path in ui.series(paths):
        with open_(path) as infile:
            rec = list(MinimalGenbankParser(infile))[0]
            db = GenbankAnnotationDb(
                source=write_path,
                data=rec.pop("features", None),
                seqid=rec["locus"],
                db=db,
            )
            one_valid_path = True

    if not one_valid_path:
        raise IOError(f"{str(path)!r} not found")

    db.make_indexes()
    return db


def _leave_attributes(*attrs):
    return attrs[0]


@display_wrap
def _db_from_gff(
    path: PathType,
    seqids: OptionalStrContainer,
    db: typing.Optional[SupportsFeatures],
    write_path: PathType,
    num_lines: OptionalInt,
    **kwargs,
) -> SupportsFeatures:
    from cogent3.parse.gff import gff_parser, is_gff3

    paths = pathlib.Path(path)
    paths = list(paths.parent.glob(paths.name))

    ui = kwargs.pop("ui")
    one_valid_path = False
    seen_ids = set()
    for path in ui.series(paths):
        num_fake_ids = 0
        gff3 = is_gff3(path)
        db = GffAnnotationDb(source=write_path, db=db)
        for block in iter_line_blocks(path, num_lines=num_lines):
            data = list(
                gff_parser(
                    block, seqids=seqids, attribute_parser=_leave_attributes, gff3=gff3
                )
            )
            data, num_fake_ids = merged_gff_records(data, num_fake_ids)
            if already_seen := seen_ids & data.keys():
                for name in already_seen:
                    db.update_record_spans(name=name, spans=data[name].spans)

            seen_ids |= data.keys()
            db.add_records(data)
        one_valid_path = True
    if not one_valid_path:
        raise IOError(f"{str(path)!r} not found")

    db.make_indexes()
    return db


def load_annotations(
    *,
    path: PathType,
    seqids: OptionalStr = None,
    db: typing.Optional[SupportsFeatures] = None,
    write_path: PathType = ":memory:",
    lines_per_block: OptionalInt = 500_000,
    show_progress: bool = False,
) -> SupportsFeatures:
    """loads annotations from flatfile into a db

    Parameters
    ----------
    path
        path to a plain text file containing features
    seqids
        only features whose seqid matches a provided identifier are returned,
        the default is all features.
    db
        an existing feature db to add these records to. Must be of a
        compatible type.
    write_path
        where the constructed database should be written, defaults to
        memory only
    lines_per_block
        number of lines to insert into the db per iteration. This can help with
        memory usage. Only applies to gff files.
    show_progress
        applied only if loading features from multiple files

    Notes
    -----
    We DO NOT check if a provided db already contains records from a flatfile.
    """
    if seqids is not None:
        seqids = {seqids} if isinstance(seqids, str) else set(seqids)
    path = pathlib.Path(path)
    return (
        _db_from_genbank(
            path, db=db, write_path=write_path, show_progress=show_progress
        )
        if {".gb", ".gbk"} & set(path.suffixes)
        else _db_from_gff(
            path,
            seqids=seqids,
            db=db,
            write_path=write_path,
            show_progress=show_progress,
            num_lines=lines_per_block,
        )
    )


def _update_array_format(data: bytes) -> bytes:
    """Convert from the old to the new sqlite numpy array format.

    Converts the previous format saved with numpy.ndarray.tobytes
    to the .npy format generated with numpy.save.

    Ignores any entries saved in the new format

    Parameters
    ----------
    data : bytes
        The sqlite3 representation of the numpy array.

    Returns
    -------
    bytes
        The new sqlite3 representation of the numpy array.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        new_data = array_to_sqlite(sqlite_to_array(data))
    return new_data


def update_file_format(
    source_path: PathType,
    db_class: typing.Union[
        type[BasicAnnotationDb],
        type[GenbankAnnotationDb],
        type[GffAnnotationDb],
    ],
    backup: bool = True,
) -> None:
    """Update the database file to the latest format.

    Fixes an OS dependent issue under the previous representation.

    Ensure update_file_format is run on the same OS used to generate the
    database file.

    By default performs a backup of the database to another file before
    updating the given file to the new format. This is in case the function
    is run on a different OS than the one used to generate the file where
    the spans column may become corrupted. The backup file is saved as the
    original file path appended with ".bak".

    See https://github.com/cogent3/cogent3/issues/1776 for more details.

    Parameters
    ----------
    source_path : PathType
        The database file to reformat.
    db_class : typing.Union[ type[BasicAnnotationDb], type[GenbankAnnotationDb], type[GffAnnotationDb], ]
        The type of database the file is.
    backup : bool, optional
        If True (default), performs a backup of the database before updating.
        Otherwise does not perform a backup prior to update (not recommended).
    """
    source_path = pathlib.Path(source_path).expanduser()

    if not source_path.exists():
        raise OSError(f"File {source_path} does not exist.")

    anno_db = db_class(source=source_path)

    if backup:
        backup_path = source_path.parent / f"{source_path.name}.bak"
        if backup_path.exists():
            raise FileExistsError(
                f"Backup file already exists for {source_path}. "
                f"If there was a problem with the conversion process, "
                f"update_file_format should be run on the backed up file. "
                f"Ensure update_file_format is run on the same OS used to generate the file."
            )
        anno_db.write(backup_path)

    conn = anno_db.db
    conn.create_function("update_array_format", 1, _update_array_format)
    for table_name in anno_db.table_names:
        cursor = conn.cursor()
        array_columns = [
            r["name"]
            for r in cursor.execute(f"PRAGMA table_info({table_name});").fetchall()
            if r["type"] == "array"
        ]

        for column in array_columns:
            cursor.execute(
                f"UPDATE {table_name} SET {column}=update_array_format({column});"
            )
        conn.commit()
