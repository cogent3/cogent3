from __future__ import annotations

import collections
import inspect
import json
import os
import pathlib
import re
import sqlite3
import sys
import typing

import numpy

from cogent3._version import __version__
from cogent3.util.deserialise import register_deserialiser
from cogent3.util.misc import get_object_provenance
from cogent3.util.table import Table


OptionalInt = typing.Optional[int]
OptionalStr = typing.Optional[str]
OptionalStrList = typing.Optional[typing.Union[str, typing.List[str]]]
OptionalBool = typing.Optional[bool]
OptionalDbCursor = typing.Optional[sqlite3.Cursor]
ReturnType = typing.Tuple[str, tuple]  # the sql statement and corresponding values
# data type for sqlitedb constructors
T = typing.Optional[typing.Iterable[dict]]

# used for presence of sqlite feature
_is_ge_3_11 = (sys.version_info.major, sys.version_info.minor) >= (3, 11)


# Define custom types for storage in sqlite
# https://stackoverflow.com/questions/18621513/python-insert-numpy-array-into-sqlite3-database
def array_to_sqlite(data):
    return data.tobytes()


def sqlite_to_array(data):
    result = numpy.frombuffer(data, dtype=int)
    dim = result.shape[0] // 2
    return result.reshape((dim, 2))


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
    reversed: bool  # True if feature on reverse strand
    on_alignment: bool  # True if feature on an alignment


@typing.runtime_checkable
class SerialisableType(typing.Protocol):
    def to_rich_dict(self) -> dict:
        ...

    def to_json(self) -> str:
        ...

    def from_dict(self, data):
        ...


@typing.runtime_checkable
class SupportsQueryFeatures(typing.Protocol):  # should be defined centrally
    def get_features_matching(
        self,
        *,
        seqid: OptionalStr = None,
        biotype: OptionalStr = None,
        name: OptionalStr = None,
        start: OptionalInt = None,
        end: OptionalInt = None,
        strand: OptionalStr = None,
        attributes: OptionalStr = None,
        on_alignment: OptionalBool = None,
    ) -> typing.Iterator[FeatureDataType]:
        ...

    def get_feature_children(
        self, *, name: str, start: OptionalInt = None, end: OptionalInt = None
    ) -> typing.List[FeatureDataType]:
        ...

    def get_feature_parent(
        self, *, name: str, start: OptionalInt = None, end: OptionalInt = None
    ) -> typing.List[FeatureDataType]:
        ...

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
        ...


@typing.runtime_checkable
class SupportsWriteFeatures(typing.Protocol):  # should be defined centrally
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
    ) -> None:
        ...

    def add_records(
        self,
        *,
        records: typing.Sequence[dict],
        # seqid required for genbank
        seqid: OptionalStr = None,
    ) -> None:
        ...

    def update(self, annot_db, seqids: OptionalStrList = None) -> None:
        # update records with those from an instance of the same type
        ...


@typing.runtime_checkable
class SupportsFeatures(
    SupportsQueryFeatures, SupportsWriteFeatures, SerialisableType, typing.Protocol
):
    @property
    def db(self):
        # pointer to the actual db
        ...

    def __len__(self):
        # the number of records
        ...

    def __eq__(self):
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
        if False, only records within start, end are included. If True,
        all records that overlap the segment defined by start, end are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """
    # todo this needs to support OR operation on some conditions, e.g. if the value of
    # a condition is a tuple, do OR

    start = conditions.pop("start", None)
    end = conditions.pop("end", None)

    sql = []
    vals = ()
    if conditions:
        conds = []
        vals = []
        for col, val in conditions.items():
            # conditions are filtered for None before here, so we should add
            # an else where the op is assigned !=
            if val is not None:
                op = "LIKE" if isinstance(val, str) and "%" in val else "="
                conds.append(f"{col} {op} ?")
                vals.append(val)
        sql.append(" AND ".join(conds))
        vals = tuple(vals)

    if isinstance(start, int) and isinstance(end, int):
        if allow_partial:
            # allow matches that overlap the segment
            cond = [
                f"(start >= {start} AND end <= {end})",  # lies within the segment
                f"(start <= {start} AND end > {start})",  # straddles beginning of segment
                f"(start < {end} AND end >= {end})",  # straddles end of segment
                f"(start <= {start} AND end >= {end})",  # includes segment
            ]
            cond = " OR ".join(cond)
        else:
            # only matches within bounds
            cond = f"start >= {start} AND end <= {end}"
        sql.append(f"({cond})")
    elif isinstance(start, int):
        # if query has no end, then any feature containing start
        cond = f"(start <= {start} AND {start} < end)"
        sql.append(f"({cond})")
    elif isinstance(end, int):
        # if query has no start, then any feature containing end
        cond = f"(start <= {end} AND {end} < end)"
        sql.append(f"({cond})")

    sql = f"{' AND '.join(sql)}"
    return sql, vals


def _del_records_sql(
    table_name: str,
    conditions: dict,
    start: OptionalInt = None,
    end: OptionalInt = None,
    allow_partial=True,
) -> ReturnType:
    """creates the SQL and values for identifying records to be deleted

    Parameters
    ----------
    table_name : str
        table to have records deleted from
    conditions : dict
        column name and values to be matched
    start, end : OptionalInt
        select records whose (start, end) values lie between start and end,
        or overlap them if (allow_partial is True)
    allow_partial : bool, optional
        if False, only records within start, end are included. If True,
        all records that overlap the segment defined by start, end are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """
    where, vals = _matching_conditions(
        conditions=conditions, start=start, end=end, allow_partial=allow_partial
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
    end: OptionalInt = None,
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
    start, end : OptionalInt
        select records whose (start, end) values lie between start and end,
        or overlap them if (allow_partial is True)
    allow_partial : bool, optional
        if False, only records within start, end are included. If True,
        all records that overlap the segment defined by start, end are included.

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
    end: OptionalInt = None,
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
    start, end : OptionalInt
        select records whose (start, end) values lie between start and end,
        or overlap them if (allow_partial is True)
    allow_partial : bool, optional
        if False, only records within start, end are included. If True,
        all records that overlap the segment defined by start, end are included.

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


class SqliteAnnotationDbMixin:
    # table schema for user provided annotations
    _user_schema = {
        "biotype": "TEXT",
        "seqid": "TEXT",
        "name": "TEXT",
        "parent_id": "TEXT",
        "start": "INTEGER",
        "end": "INTEGER",
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

    def __deepcopy__(self, memodict=None):
        memodict = memodict or {}
        if _is_ge_3_11:
            new = self.__class__(source=self.source)
            new._db.deserialize(self._db.serialize())
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
        return type(self) == type(other) and other.db is self.db

    @property
    def table_names(self) -> tuple[str]:
        return self._table_names

    def _init_tables(self):
        # bit of magic
        # assumes schema attributes named as `_<table name>_schema`
        for table_name in self.table_names:
            attr = getattr(self, f"_{table_name}_schema")
            sql = _make_table_sql(table_name, attr)
            self._execute_sql(sql)

    @property
    def db(self) -> sqlite3.Connection:
        if self._db is None:
            self._db = sqlite3.connect(
                self.source, detect_types=sqlite3.PARSE_DECLTYPES
            )
            self._db.row_factory = sqlite3.Row

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
        # the start, end variables are added into the record via the loop
        # on local variables
        start = int(spans.min())
        end = int(spans.max())
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
        columns: list[str],
        column: str,
        name: str,
        start: OptionalInt = None,
        end: OptionalInt = None,
        biotype: OptionalStr = None,
        allow_partial: bool = False,
    ) -> typing.List[FeatureDataType]:
        # we return the parent_id because `get_feature_parent()` requires it
        sql, vals = _select_records_sql(
            table_name=table_name,
            conditions={column: name, "biotype": biotype},
            columns=columns,
            start=start,
            end=end,
        )
        for result in self._execute_sql(sql, values=vals):
            result = dict(zip(result.keys(), result))
            result["on_alignment"] = result.get("on_alignment", None)
            result["spans"] = [tuple(c) for c in result["spans"]]
            result["reversed"] = result.pop("strand", None) == "-"
            yield result

    def get_feature_children(
        self, name: str, biotype: OptionalStr = None, **kwargs
    ) -> typing.Iterator[FeatureDataType]:
        """yields children of name"""
        # kwargs is used because other classes need start / end
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
                table_name=table_name, columns=columns, column="name", name=f"%{name}%"
            ):
                # multiple values for parent means this is better expressed
                # as an OR clause
                # todo modify the conditional SQL generation
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

    def get_records_matching(
        self,
        biotype: str = None,
        seqid: str = None,
        name: str = None,
        start: int = None,
        end: int = None,
        strand: bool = None,
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

    def get_features_matching(
        self,
        *,
        biotype: str = None,
        seqid: str = None,
        name: str = None,
        start: int = None,
        end: int = None,
        strand: bool = None,
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
                result["on_alignment"] = result.get("on_alignment", None)
                result["spans"] = [tuple(c) for c in result["spans"]]
                result["reversed"] = result.pop("strand", None) == "-"
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

    def count_distinct(
        self,
        *,
        seqid: bool = False,
        biotype: bool = False,
        name: bool = False,
    ) -> typing.Optional[Table]:
        """return table of counts of distinct values"""
        conditions = {k for k, v in locals().items() if k != "self" and v}
        header = list(conditions)
        columns = ", ".join(conditions)
        if not columns:
            return
        placeholder = "{}"
        if len(conditions) == 1:
            sql_template = f"SELECT {columns}, COUNT(DISTINCT {columns}) FROM {placeholder} GROUP BY {columns};"
        else:
            sql_template = f"SELECT {columns}, COUNT(*) as count FROM {placeholder} GROUP BY {columns};"

        data = []
        for table in self.table_names:
            sql = sql_template.format(table)
            if result := self._execute_sql(sql).fetchall():
                data.extend(tuple(r) for r in result)
        return Table(header=header + ["count"], data=data)

    @property
    def describe(self) -> Table:
        """top level description of the annotation db"""
        from cogent3 import make_table

        sql_template = "SELECT {}, COUNT(DISTINCT {}) FROM {} GROUP BY {};"
        data = {}
        for column in ("seqid", "biotype"):
            for table in self.table_names:
                sql = sql_template.format(column, column, table, column)
                if result := self._execute_sql(sql).fetchall():
                    counts = dict(tuple(r) for r in result)
                    for distinct, count in counts.items():
                        key = f"{column}({distinct!r})"
                        data[key] = data.get(key, 0) + count

        row_counts = []
        for table in self.table_names:
            result = self._execute_sql(f"SELECT COUNT(*) FROM {table}").fetchone()
            row_counts.append(result["COUNT(*)"])

        table = make_table(
            data={
                "": list(data.keys()) + [f"num_rows({t!r})" for t in self.table_names],
                "count": [v for v in data.values()] + row_counts,
            }
        )
        return table

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

    def update(
        self, annot_db: SupportsFeatures, seqids: OptionalStrList = None
    ) -> None:
        """update records with those from an instance of the same type"""
        if not isinstance(annot_db, type(self)):
            raise TypeError(f"{type(annot_db)} != {type(self)}")

        self._update_db_from_rich_dict(annot_db.to_rich_dict(), seqids=seqids)

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

    def to_json(self) -> str:
        return json.dumps(self.to_rich_dict())


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
        "end": "INTEGER",
        "score": "TEXT",  # check defn
        "strand": "TEXT",
        "phase": "TEXT",
        "attributes": "TEXT",
        "comments": "TEXT",
        "spans": "array",  # aggregation of coords across records
        "name": "TEXT",
        "parent_id": "TEXT",
    }

    def __init__(
        self, *, data: T = None, db: OptionalDbCursor = None, source=":memory:"
    ):
        data = data or []
        # note that data is destroyed
        self._num_fakeids = 0
        self.source = source
        self._db = db
        if db is None:
            self._init_tables()
        data = self._merged_data(data)
        self.add_records(data)

    def add_records(self, reduced: dict) -> None:
        # Can we really trust text case of "ID" and "Parent" to be consistent
        # across sources of gff?? I doubt it, so more robust regex likely
        # required
        for record in reduced.values():
            # our Feature code assumes start always < end,
            # we record direction using Strand
            spans = numpy.array(sorted(record["spans"]), dtype=int)  # sorts the rows
            spans.sort(axis=1)
            record["Start"] = int(spans.min())
            record["End"] = int(spans.max())
            record["spans"] = spans
            cmnd, vals = _add_record_sql(
                self.table_names[0],
                {
                    k: v
                    for k, v in record.items()
                    if isinstance(v, numpy.ndarray) or v not in (".", None)
                },
            )
            self._execute_sql(cmnd=cmnd, values=vals)

        del reduced

    def _merged_data(self, records) -> dict:
        field_template = r"(?<={}=)[^;\s]+"
        name = re.compile(field_template.format("ID"))
        parent_id = re.compile(field_template.format("Parent"))

        reduced = collections.OrderedDict()
        # collapse records with ID's occurring in multiple rows into one
        # row, converting their coordinates
        # extract the name from ID and add this into the table
        # I am not convinced we can rely on gff files to be ordered,
        # if we could, we could do this as one pass over the data
        while records:
            record = records.pop(0)
            record["biotype"] = record.pop("Type")
            attrs = record["Attributes"] or ""
            if match := name.search(attrs):
                record_id = match.group()
            else:
                record_id = f"unknown-{self._num_fakeids}"
                self._num_fakeids += 1

            record["name"] = record_id
            if pid := parent_id.search(attrs):
                record["parent_id"] = pid.group()

            if record_id not in reduced:
                reduced[record_id] = record
                reduced[record_id]["spans"] = []

            # should this just be an append?
            reduced[record_id]["spans"].append((record["Start"], record["End"]))
        return reduced


# The GenBank format is less clear on the relationship between identifiers,
# e.g. it can be gene -> SRP_RNA -> exon, gene -> CDS, etc... Without having
# the explicit type-hierarchy from NCBI (which I have not yet found) the only way
# to establish the Parent of a feature is by knowing what other
# features have the same ID and then which of those are candidates to contain
# a feature based on the hierarchy of relationships.
# In other words, we assume that children have the same ID as their parent BUT but
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
        "end": "INTEGER",
        "strand": "TEXT",
        "comments": "TEXT",
        "spans": "array",  # aggregation of coords across records
        "name": "TEXT",
        "symbol": "TEXT",
        "parent_id": "TEXT",
        "attributes": "json",
    }

    def __init__(
        self,
        *,
        data: T = None,
        seqid: OptionalStr = None,
        db: OptionalDbCursor = None,
        source=":memory:",
    ):
        data = data or []
        # note that data is destroyed
        self._db = db
        self._num_fakeids = 0
        self.source = source
        self._init_tables()
        self.add_records(data, seqid)

    def add_records(self, records, seqid):
        # need to capture genetic code from genbank, but what a

        col_key_map = {"type": "biotype", "locus_tag": "name"}
        exclude = {"translation", "location"}  # location is grabbed directly
        for record in records:
            # our Feature code assumes start always < end,
            store = {"seqid": seqid}
            # we create the location data directly
            if location := record.get("location", None):
                store["spans"] = numpy.array(location.get_coordinates(), dtype=int)
                if strand := location.strand:
                    store["strand"] = "-" if strand == -1 else "+"
                store["start"] = int(store["spans"].min())
                store["end"] = int(store["spans"].max())

            record_keys = col_key_map.keys() & record.keys()
            attrs_keys = record.keys() - col_key_map.keys() - exclude
            store.update({col_key_map[k]: record[k] for k in record_keys})
            store["attributes"] = {k: record[k] for k in attrs_keys}
            if "name" not in store:
                store["name"] = [f"fakeid-{self._num_fakeids}"]
                self._num_fakeids += 1

            store["name"] = ",".join(store["name"])
            cmnd, vals = _add_record_sql(
                self.table_names[0],
                store,
            )
            self._execute_sql(cmnd=cmnd, values=vals)

    def get_feature_children(
        self,
        name: str,
        biotype: OptionalStr = None,
        exclude_biotype: OptionalStr = None,
        start: OptionalInt = None,
        end: OptionalInt = None,
    ) -> typing.Iterator[FeatureDataType]:
        """yields children of name"""
        # children must lie within the parent coordinates
        for table_name in self.table_names:
            columns = (
                "seqid",
                "biotype",
                "spans",
                "strand",
                "name",
                "parent_id",
                "start",
                "end",
            )
            if table_name == "user":
                columns += ("on_alignment",)

            for feat in self._get_feature_by_id(
                table_name=table_name,
                columns=columns,
                column="name",
                name=name,
                biotype=biotype,
                start=int(start),
                end=int(end),
                allow_partial=False,
            ):
                if feat["biotype"] == exclude_biotype:
                    continue
                cstart, cend = feat.pop("start"), feat.pop("end")
                if not (start <= cstart < end and start < cend <= end):
                    continue
                yield feat

    def get_feature_parent(
        self,
        name: str,
        exclude_biotype: OptionalStr = None,
        start: OptionalInt = None,
        end: OptionalInt = None,
    ) -> typing.Iterator[FeatureDataType]:
        """yields parents of name"""
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
                "end",
            )
            if table_name == "user":
                columns += ("on_alignment",)
            for feat in self._get_feature_by_id(
                table_name=table_name,
                columns=columns,
                column="name",
                name=name,
                start=int(start),
                end=int(end),
                allow_partial=False,
            ):
                # add support for != operation to SQL where clause generation
                if feat["biotype"] == exclude_biotype:
                    continue

                cstart, cend = feat.pop("start"), feat.pop("end")
                if cstart > start or end > cend:
                    continue
                yield feat


@register_deserialiser(get_object_provenance(GffAnnotationDb))
def deserialise_gff_db(data: dict):
    return GffAnnotationDb.from_dict(data)


@register_deserialiser(get_object_provenance(GenbankAnnotationDb))
def deserialise_gb_db(data: dict):
    return GenbankAnnotationDb.from_dict(data)


@register_deserialiser("annotation_to_annotation_db")
def convert_annotation_to_annotation_db(data: dict) -> dict:
    from cogent3.util.deserialise import deserialise_map_spans

    db = GffAnnotationDb()

    seqid = data.pop("name", None)
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


def _db_from_genbank(path, db):
    from cogent3 import open_
    from cogent3.parse.genbank import RichGenbankParser

    with open_(path) as infile:
        data = list(RichGenbankParser(infile, db=db))[0][1]

    return data.annotation_db


def _leave_attributes(*attrs):
    return attrs[0]


def _db_from_gff(path, seqids, db):
    from cogent3.parse.gff import gff_parser

    data = list(
        gff_parser(
            path,
            seqids=seqids,
            attribute_parser=_leave_attributes,
        )
    )
    return GffAnnotationDb(data=data, db=db)


def load_annotations(
    path: os.PathLike, seqids: OptionalStr = None, db: OptionalDbCursor = None
) -> SupportsFeatures:
    if seqids is not None:
        seqids = {seqids} if isinstance(seqids, str) else set(seqids)
    path = pathlib.Path(path)
    return (
        _db_from_genbank(path, db=db)
        if ".gb" in path.suffixes
        else _db_from_gff(path, seqids=seqids, db=db)
    )
