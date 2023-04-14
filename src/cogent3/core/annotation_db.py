import collections
import json
import os
import pathlib
import re
import sqlite3
import typing

import numpy

from cogent3.util.deserialise import register_deserialiser
from cogent3.util.misc import get_object_provenance
from cogent3.util.table import Table


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


OptionalInt = typing.Optional[int]
OptionalStr = typing.Optional[str]
OptionalBool = typing.Optional[bool]
ReturnType = typing.Tuple[str, tuple]  # the sql statement and corresponding values


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
    biotype: str  # rename type attr of cogent3 Annotatables to match this?
    name: str  # rename to name to match cogent3 Annotatable.name?
    spans: list[tuple[int, int]]
    reversed: bool  # True if feature on reverse strand


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
        biotype: OptionalStr = None,
        name: OptionalStr = None,
        start: OptionalInt = None,
        end: OptionalInt = None,
        strand: OptionalStr = None,
        on_alignment: OptionalBool = None,
    ) -> typing.Iterator[FeatureDataType]:
        ...

    def get_feature_children(
        self, name: str, start: OptionalInt = None, end: OptionalInt = None
    ) -> typing.List[FeatureDataType]:
        ...

    def get_feature_parent(
        self, name: str, start: OptionalInt = None, end: OptionalInt = None
    ) -> typing.List[FeatureDataType]:
        ...


@typing.runtime_checkable
class SupportsWriteFeatures(typing.Protocol):  # should be defined centrally
    def add_feature(
        self,
        seqid: str,  # is seqid and name referring to the same thing?
        biotype: str,
        name: str,
        spans: typing.List[typing.Tuple[int, int]],
        strand: str = None,
        on_alignment: bool = False,
    ) -> None:
        ...


@typing.runtime_checkable
class SupportsFeatures(SupportsQueryFeatures, SupportsWriteFeatures, typing.Protocol):
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
    partial: bool = True,
):
    """creates WHERE clause

    Parameters
    ----------
    conditions : dict
        column name and values to be matched
    partial : bool, optional
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
            if val:
                op = "LIKE" if "%" in val else "="
                conds.append(f"{col} {op} ?")
                vals.append(val)
        sql.append(" AND ".join(conds))
        vals = tuple(vals)

    if isinstance(start, int) and isinstance(end, int):
        # only matches within bounds
        if partial:
            cond = [
                f"(start <= {start} AND end >= {start})",  # straddles beginning of segment
                f"(start < {end} AND end >= {end})",  # straddles end of segment
                f"(start <= {start} AND end >= {end})",  # includes segment
            ]
            cond = " OR ".join(cond)
        else:
            cond = f"start >= {start} AND end <= {end}"
        sql.append(f"({cond})")
    elif isinstance(start, int):
        if partial:
            cond = " OR ".join(
                [f"start >= {start}", f"(start <= {start} AND end >= {start})"]
            )
        else:
            cond = f"start >= {start}"
        sql.append(f"({cond})")
    elif isinstance(end, int):
        if partial:
            cond = " OR ".join([f"end <= {end}", f"(start <= {end} AND end >= {end})"])
        else:
            cond = f"end <= {end}"
        sql.append(f"({cond})")

    sql = f"{' AND '.join(sql)}"
    return sql, vals


def _del_records_sql(
    table_name: str,
    conditions: dict,
    start: OptionalInt = None,
    end: OptionalInt = None,
    partial=True,
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
        or overlap them if (partial is True)
    partial : bool, optional
        if False, only records within start, end are included. If True,
        all records that overlap the segment defined by start, end are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """
    where, vals = _matching_conditions(
        conditions=conditions, start=start, end=end, partial=partial
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
    partial=True,
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
        or overlap them if (partial is True)
    partial : bool, optional
        if False, only records within start, end are included. If True,
        all records that overlap the segment defined by start, end are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """

    where, vals = _matching_conditions(conditions=conditions, partial=partial)
    columns = f"{', '.join(columns)}" if columns else "*"
    sql = f"SELECT {columns} FROM {table_name}"
    if not where:
        return sql

    sql = f"{sql} WHERE {where};"
    return sql, vals


# todo add support for querying for text within the additional feature,
# for example, add a "attrs_like" argument which users can provide text
# that will be treated as surrounded by wild-cards
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
        "on_alignment": "INT",
    }

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
        seqid: str,
        biotype: str,
        name: str,
        spans: typing.List[typing.Tuple[int, int]],
        strand: str = None,
        on_alignment: bool = False,
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
    ) -> typing.Iterator[FeatureDataType]:
        """return all fields"""
        columns = kwargs.pop("columns", None)
        sql, vals = _select_records_sql(
            table_name=table_name, conditions=kwargs, columns=columns
        )
        yield from self._execute_sql(sql, values=vals)

    def get_records_matching(
        self,
        biotype: str = None,
        seqid: str = None,
        name: str = None,
        start: int = None,
        end: int = None,
        strand: bool = None,
        on_alignment: bool = None,
    ) -> typing.Iterator[dict]:
        """return all fields for matching records"""
        # a record is Everything, a Feature is a subset
        # we define query as all defined variables from local name space,
        # excluding "self" and kwargs at default values
        kwargs = {k: v for k, v in locals().items() if k != "self" and v is not None}
        for table_name in self.table_names:
            for result in self._get_records_matching(table_name, **kwargs):
                yield {k: result[k] for k in result.keys()}

    def get_features_matching(
        self,
        biotype: str = None,
        seqid: str = None,
        name: str = None,
        start: int = None,
        end: int = None,
        strand: bool = None,
        on_alignment: bool = None,
    ) -> typing.Iterator[FeatureDataType]:
        # returns essential values to create a Feature
        # we define query as all defined variables from local name space,
        # excluding "self" and kwargs at default values
        kwargs = {k: v for k, v in locals().items() if k != "self" and v is not None}
        columns = ("biotype", "spans", "strand", "name")
        for table_name in self.table_names:
            for result in self._get_records_matching(
                table_name=table_name, columns=columns, **kwargs
            ):
                yield FeatureDataType(
                    biotype=result["biotype"],
                    name=result["name"],
                    spans=[tuple(c) for c in result["spans"]],
                    reversed=result["strand"] == "-",
                )

    @property
    def describe(self) -> Table:
        """top level description of the annotation db"""
        sql_template = "SELECT DISTINCT {} FROM {};"
        data = {}
        for column in ("seqid", "biotype"):
            data[column] = set()
            for table in self.table_names:
                sql = sql_template.format(column, table)
                if result := self._execute_sql(sql).fetchall():
                    data[column] |= {r[column] for r in result}

        row_counts = []
        for table in self.table_names:
            result = self._execute_sql(f"SELECT COUNT(*) FROM {table}").fetchone()
            row_counts.append(result["COUNT(*)"])

        from cogent3 import make_table

        table = make_table(
            data={
                "biotype": list(data.keys())
                + [f"num_rows({t!r})" for t in self.table_names],
                "count": [len(v) for v in data.values()] + row_counts,
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

    def to_rich_dict(self):
        # todo drop columns from records with no value
        # top level dict for each table_name
        records = {}
        for table in self.table_names:
            records[table] = []
            for record in self._execute_sql(f"SELECT * from {table};"):
                record = {k: record[k] for k in record.keys()}
                record["spans"] = record[
                    "spans"
                ].tolist()  # required for json serialisation
                records[table].append(
                    {k: v for k, v in record.items() if v is not None}
                )

        records["type"] = get_object_provenance(self)
        return records

    @classmethod
    def from_dict(cls, data: dict):
        data.pop("type")
        # make an empty db
        db = cls([])
        for table_name, records in data.items():
            for record in records:
                record["spans"] = numpy.array(record["spans"], dtype=int)
                sql, vals = _add_record_sql(
                    table_name, {k: v for k, v in record.items() if v is not None}
                )
                db._execute_sql(sql, vals)
        return db


class GffAnnotationDb(SqliteAnnotationDbMixin):
    _table_names = "gff", "user"
    # am relying on name structured as _<table name>_schema
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

    def __init__(self, data, db=None, source=":memory:"):
        # note that data is destroyed
        self._num_fakeids = 0
        self.source = source
        if db is None:
            self._db = sqlite3.connect(
                self.source, detect_types=sqlite3.PARSE_DECLTYPES
            )
            self._db.row_factory = sqlite3.Row
        else:
            self._db = db
        self._init_tables()
        data = self._merged_data(data)
        self._populate_from_records(data)

    def _populate_from_records(self, reduced):
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

    def _get_feature_by_id(
        self, table_name: str, column: str, name: str, biotype: OptionalStr = None
    ) -> typing.List[FeatureDataType]:
        # we return the parent_id because `get_feature_parent()` requires it
        sql, vals = _select_records_sql(
            table_name=table_name,
            conditions={column: name, "biotype": biotype},
            columns=["biotype", "spans", "strand", "name", "parent_id"],
        )
        for result in self._execute_sql(sql, values=vals):
            yield {
                "biotype": result["biotype"],
                "name": result["name"],
                "spans": [tuple(c) for c in result["spans"]],
                "reversed": result["strand"] == "-",
                "parent_id": result["parent_id"],
            }

    def get_feature_children(
        self, name: str, biotype: OptionalStr = None, **kwargs
    ) -> typing.Iterator[FeatureDataType]:
        """yields children of name"""
        # kwargs is used because other classes need start / end
        # just uses search matching
        for table_name in self.table_names:
            for result in self._get_feature_by_id(
                table_name=table_name,
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
            for result in self._get_feature_by_id(
                table_name=table_name, column="name", name=f"%{name}%"
            ):
                # multiple values for parent means this is better expressed
                # as an OR clause
                # todo modify the conditional SQL generation
                for name in result["parent_id"].replace(" ", "").split(","):
                    if parent := list(
                        self._get_feature_by_id(
                            table_name=table_name,
                            column="name",
                            name=name,
                        )
                    ):
                        parent = parent[0]
                        parent.pop(
                            "parent_id"
                        )  # remove invalid field for the FeatureDataType
                        yield parent


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
    _table_names = "gb", "user"
    # am relying on name structured as _<table name>_schema
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

    def __init__(self, data: typing.List[dict], seqid: str, source=":memory:", db=None):
        # note that data is destroyed
        self._db = db
        self._num_fakeids = 0
        self.source = source
        self._init_tables()
        self._populate_from_records(data, seqid)

    # todo: rename add_records, make public method - can be used to add genbank records to an existing GenbankAnnotationDb
    def _populate_from_records(self, records, seqid):
        # Can we really trust case to be consistent across sources of gff??
        # I doubt it, so more robust regex likely required
        # Can we even rely on these concepts being represented by the text
        # ID and Parent?

        # need to capture genetic code from genbank, but what a

        col_key_map = {"type": "biotype", "locus_tag": "name"}
        exclude = {"translation", "location"}  # location is grabbed directly
        for record in records:
            # our Feature code assumes start always < end,
            store = {"seqid": seqid}
            # we create the location data directly
            if location := record.get("location", None):
                store["spans"] = numpy.array(
                    sorted([sorted((s.first(), s.last())) for s in location]), dtype=int
                )
                if strand := location.strand:
                    store["strand"] = strand
                store["start"] = int(store["spans"].min())
                store["end"] = int(store["spans"].max())

            record_keys = col_key_map.keys() & record.keys()
            attrs_keys = record.keys() - col_key_map.keys() - exclude
            store |= {col_key_map[k]: record[k] for k in record_keys}
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

    def _get_feature_by_id(
        self,
        table_name: str,
        name: str,
        start: int,
        end: int,
        biotype: OptionalStr = None,
        partial: bool = False,
    ) -> typing.List[FeatureDataType]:
        # we return the parent_id because `get_feature_parent()` requires it
        sql, vals = _select_records_sql(
            table_name=table_name,
            conditions={"name": name, "biotype": biotype},
            start=start,
            end=end,
            columns=["biotype", "start", "end", "spans", "strand", "name", "parent_id"],
            partial=partial,
        )
        for result in self._execute_sql(sql, values=vals):
            yield {
                "biotype": result["biotype"],
                "name": name,
                "start": result["start"],
                "end": result["end"],
                "spans": [tuple(c) for c in result["spans"]],
                "reverse": result["strand"] == "-",
            }

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
        for table in self.table_names:
            for feat in self._get_feature_by_id(
                table_name=table,
                name=name,
                biotype=biotype,
                start=int(start),
                end=int(end),
                partial=False,
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
        for table in self.table_names:
            for feat in self._get_feature_by_id(
                table_name=table,
                name=name,
                start=int(start),
                end=int(end),
                partial=False,
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


def _db_from_genbank(path, db):
    from cogent3 import open_
    from cogent3.parse.genbank import MinimalGenbankParser

    with open_(path) as infile:
        data = list(MinimalGenbankParser(infile))

    return GenbankAnnotationDb(data[0]["features"], data[0]["locus"], db=db)


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


def load_annotations(path: os.PathLike, seqids=None, db=None) -> SupportsFeatures:
    if seqids is not None:
        seqids = {seqids} if isinstance(seqids, str) else set(seqids)
    path = pathlib.Path(path)
    return (
        _db_from_genbank(path, db=db)
        if ".gb" in path.suffixes
        else _db_from_gff(path, seqids=seqids, db=db)
    )
