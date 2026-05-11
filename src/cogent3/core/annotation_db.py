from __future__ import annotations

import abc
import collections
import contextlib
import copy
import inspect
import io
import json
import pathlib
import sqlite3
import warnings
from collections.abc import Callable, Iterable, Iterator
from collections.abc import Sequence as PySeq
from typing import (
    TYPE_CHECKING,
    Any,
    ClassVar,
    Protocol,
    Self,
    TypedDict,
    cast,
    runtime_checkable,
)

import numpy
import numpy.typing as npt
from scinexus import warning as snx_warn
from scinexus.deserialise import deserialise_object, register_deserialiser
from scinexus.io_util import get_format_suffixes, iter_line_blocks
from scinexus.misc import extend_docstring_from, get_object_provenance
from scinexus.progress import Progress, get_progress

from cogent3._version import __version__
from cogent3.core.location import Strand, deserialise_map_spans
from cogent3.parse.gff import GffRecordABC, merged_gff_records
from cogent3.util.io import PathType

if TYPE_CHECKING:  # pragma: no cover
    from cogent3.core.table import Table

NumpyIntArrayType = npt.NDArray[numpy.integer]

ANNDB_SUFFIX_GENBANK = "c3gbdb"
ANNDB_SUFFIX_GFF = "c3gffdb"
ANNDB_SUFFIX_BASIC = "c3andb"


# Define custom types for storage in sqlite
# https://stackoverflow.com/questions/18621513/python-insert-numpy-array-into-sqlite3-database
def array_to_sqlite(data: npt.NDArray[Any]) -> bytes:
    with io.BytesIO() as out:
        numpy.save(out, data)
        out.seek(0)
        return out.read()


def sqlite_to_array(data: bytes) -> npt.NDArray[Any]:
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
                stacklevel=2,
            )

    return result


def dict_to_sqlite_as_json(data: dict[str, Any]) -> str:
    return json.dumps(data)


def sqlite_json_to_dict(data: str | bytes) -> dict[str, Any]:
    return json.loads(data)


# registering the conversion functions with sqlite
sqlite3.register_adapter(numpy.ndarray, array_to_sqlite)
sqlite3.register_converter("array", sqlite_to_array)
sqlite3.register_adapter(dict, dict_to_sqlite_as_json)
sqlite3.register_converter("json", sqlite_json_to_dict)

# Schema version for tracking database format changes
# Version 1: Original schema with numpy array serialization for all spans
# Version 2: Revised schema with inline single-span, hierarchy table
SCHEMA_VERSION = 3

# SQL statements for creating schema tables
_LOOKUP_TABLES_SQL = {
    "seqids": """
        CREATE TABLE IF NOT EXISTS seqids (
            id INTEGER PRIMARY KEY,
            name TEXT UNIQUE NOT NULL
        )
    """,
    "biotypes": """
        CREATE TABLE IF NOT EXISTS biotypes (
            id INTEGER PRIMARY KEY,
            name TEXT UNIQUE NOT NULL
        )
    """,
    "sources": """
        CREATE TABLE IF NOT EXISTS sources (
            id INTEGER PRIMARY KEY,
            name TEXT UNIQUE NOT NULL
        )
    """,
}

_FEATURE_SPANS_SQL = """
    CREATE TABLE IF NOT EXISTS feature_spans (
        feature_id INTEGER NOT NULL,
        table_name TEXT NOT NULL,
        span_index INTEGER NOT NULL,
        start INTEGER NOT NULL,
        stop INTEGER NOT NULL,
        PRIMARY KEY (feature_id, table_name, span_index)
    )
"""

_FEATURE_HIERARCHY_SQL = """
    CREATE TABLE IF NOT EXISTS feature_hierarchy (
        child_id INTEGER NOT NULL,
        parent_id INTEGER NOT NULL,
        table_name TEXT NOT NULL,
        PRIMARY KEY (child_id, parent_id, table_name)
    )
"""

_SCHEMA_VERSION_SQL = """
    CREATE TABLE IF NOT EXISTS schema_version (
        version INTEGER NOT NULL,
        created_at TEXT DEFAULT CURRENT_TIMESTAMP
    )
"""


def _get_schema_version(db: sqlite3.Connection) -> int:
    """Get the schema version from a database, or 1 if not present."""
    try:
        cursor = db.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='schema_version'"
        )
        if not cursor.fetchone():
            return 1
        cursor = db.execute(
            "SELECT version FROM schema_version ORDER BY version DESC LIMIT 1"
        )
        row = cursor.fetchone()
        return row[0] if row else 1
    except sqlite3.Error:
        return 1


def _set_schema_version(db: sqlite3.Connection, version: int) -> None:
    """Set the schema version in the database."""
    db.execute(_SCHEMA_VERSION_SQL)
    db.execute("INSERT INTO schema_version (version) VALUES (?)", (version,))
    db.commit()


class LookupTableCache:
    """Manages caching of lookup table IDs for seqid, biotype, source normalization."""

    def __init__(self, db: sqlite3.Connection) -> None:
        self._db = db
        self._seqid_cache: dict[str, int] = {}
        self._biotype_cache: dict[str, int] = {}
        self._source_cache: dict[str, int] = {}
        self._seqid_reverse: dict[int, str] = {}
        self._biotype_reverse: dict[int, str] = {}
        self._source_reverse: dict[int, str] = {}

    def _get_or_create_id(
        self,
        table: str,
        name: str,
        cache: dict[str, int],
        reverse: dict[int, str],
    ) -> int:
        """Get ID for a name, creating it if necessary."""
        if name in cache:
            return cache[name]

        cursor = self._db.execute(f"SELECT id FROM {table} WHERE name = ?", (name,))
        row = cursor.fetchone()
        if row:
            cache[name] = row[0]
            reverse[row[0]] = name
            return row[0]

        cursor = self._db.execute(f"INSERT INTO {table} (name) VALUES (?)", (name,))
        self._db.commit()
        new_id = cursor.lastrowid
        assert new_id is not None
        cache[name] = new_id
        reverse[new_id] = name
        return new_id

    def _get_name_by_id(
        self,
        table: str,
        id_val: int,
        cache: dict[str, int],
        reverse: dict[int, str],
    ) -> str | None:
        """Get name for an ID."""
        if id_val in reverse:
            return reverse[id_val]

        cursor = self._db.execute(f"SELECT name FROM {table} WHERE id = ?", (id_val,))
        row = cursor.fetchone()
        if row:
            name = row[0]
            cache[name] = id_val
            reverse[id_val] = name
            return name
        return None

    def get_seqid_id(self, name: str) -> int:
        return self._get_or_create_id(
            "seqids", name, self._seqid_cache, self._seqid_reverse
        )

    def get_biotype_id(self, name: str) -> int:
        return self._get_or_create_id(
            "biotypes", name, self._biotype_cache, self._biotype_reverse
        )

    def get_source_id(self, name: str) -> int:
        return self._get_or_create_id(
            "sources", name, self._source_cache, self._source_reverse
        )

    def get_seqid_name(self, id_val: int) -> str | None:
        return self._get_name_by_id(
            "seqids", id_val, self._seqid_cache, self._seqid_reverse
        )

    def get_biotype_name(self, id_val: int) -> str | None:
        return self._get_name_by_id(
            "biotypes", id_val, self._biotype_cache, self._biotype_reverse
        )

    def get_source_name(self, id_val: int) -> str | None:
        return self._get_name_by_id(
            "sources", id_val, self._source_cache, self._source_reverse
        )

    def load_all(self) -> None:
        """Pre-load all lookup tables into cache."""
        for table, cache, reverse in [
            ("seqids", self._seqid_cache, self._seqid_reverse),
            ("biotypes", self._biotype_cache, self._biotype_reverse),
            ("sources", self._source_cache, self._source_reverse),
        ]:
            try:
                cursor = self._db.execute(f"SELECT id, name FROM {table}")
                for row in cursor:
                    cache[row[1]] = row[0]
                    reverse[row[0]] = row[1]
            except sqlite3.OperationalError:
                # Table doesn't exist (old schema)
                pass


class FeatureDataType(  # When does this have a start and stop key? same with parent_id
    TypedDict
):
    seqid: str  # name of the sseq
    biotype: str  # rename type attr of cogent3 Annotatables to match this?
    name: str  # rename to name to match cogent3 Annotatable.name?
    spans: list[tuple[int, int]] | NumpyIntArrayType
    strand: int | str  # will be transformed by Strand Enum
    on_alignment: bool | None  # True if feature on an alignment
    xattr: dict[str, Any]  # extra attributes


@runtime_checkable
class SerialisableType(Protocol):  # pragma: no cover
    def to_rich_dict(self) -> dict[str, Any]: ...

    def to_json(self) -> str: ...

    def from_dict(self, data: dict[str, Any]) -> None: ...


@runtime_checkable
class SupportsQueryFeatures(Protocol):  # pragma: no cover
    @snx_warn.deprecated_callable(
        version="2026.6", reason="use AnnotationDbABC instead", is_discontinued=True
    )
    def __init_subclass__(cls): ...

    # should be defined centrally
    def get_features_matching(
        self,
        *,
        biotype: str | tuple[str, ...] | list[str] | set[str] | None = None,
        seqid: str | None = None,
        name: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        strand: str | None = None,
        attributes: str | None = None,
        on_alignment: bool | None = None,
        allow_partial: bool = False,
    ) -> Iterator[FeatureDataType]: ...

    def get_feature_children(
        self,
        name: str,
        biotype: str | None = None,
        **kwargs: Any,
    ) -> Iterator[FeatureDataType]: ...

    def get_feature_parent(
        self,
        name: str,
        **kwargs: Any,
    ) -> Iterator[FeatureDataType]: ...

    def num_matches(
        self,
        *,
        seqid: str | None = None,
        biotype: str | tuple[str, ...] | list[str] | set[str] | None = None,
        name: str | None = None,
        strand: str | None = None,
        attributes: str | None = None,
        on_alignment: bool | None = None,
    ) -> int: ...

    def subset(
        self,
        *,
        source: PathType = ":memory:",
        biotype: str | None = None,
        seqid: str | None = None,
        name: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        strand: str | None = None,
        attributes: str | None = None,
        allow_partial: bool = False,
    ) -> Self: ...


@runtime_checkable
class SupportsWriteFeatures(Protocol):  # pragma: no cover
    @snx_warn.deprecated_callable(
        version="2026.6", reason="use AnnotationDbABC instead", is_discontinued=True
    )
    def __init_subclass__(cls): ...

    def add_feature(
        self,
        *,
        seqid: str,
        biotype: str,
        name: str,
        spans: list[tuple[int, int]] | NumpyIntArrayType,
        parent_id: str | None = None,
        strand: str | int | Strand | None = None,
        attributes: str | None = None,
        on_alignment: bool | None = False,
    ) -> None: ...

    def add_records(
        self,
        records: Iterable[dict[str, Any]],
        **kwargs: Any,
    ) -> None: ...

    def update(
        self,
        annot_db: AnnotationDbABC,
        seqids: str | PySeq[str] | None = None,
        **kwargs: Any,
    ) -> None:
        # update records with those from an instance of the same type
        ...

    def union(self, annot_db: AnnotationDbABC) -> AnnotationDbABC:
        # returns a new instance of the more complex class
        ...


class AnnotationDbLoaderBase(abc.ABC):
    """Base class for annotation database format loaders."""

    @property
    @abc.abstractmethod
    def name(self) -> str:
        """Name of the storage backend (e.g., 'sqlite', 'hdf5')"""
        ...

    @property
    @abc.abstractmethod
    def supported_suffixes(self) -> set[str]:
        """File suffixes this loader supports (e.g., {'.gff', '.gff3'})"""
        ...

    @abc.abstractmethod
    def load(
        self,
        path: PathType,
        seqids: str | Iterable[str] | None = None,
        db: AnnotationDbABC | None = None,
        write_path: PathType = ":memory:",
        format_name: str | None = None,
        **kwargs: Any,
    ) -> AnnotationDbABC:
        """Load annotations from file(s) into a database.

        Parameters
        ----------
        path
            Path to annotation file (may contain glob patterns)
        seqids
            Filter to specific sequence IDs
        db
            Existing database to add records to
        write_path
            Where to write the database
        format_name
            Explicit format override. If None, auto-detect from file suffix.
        **kwargs
            Additional loader-specific arguments (includes 'ui' for progress,
            'lines_per_block' for GFF chunked loading)

        Returns
        -------
        AnnotationDbABC instance with loaded annotations
        """
        ...


class AnnotationDbABC(abc.ABC):
    @abc.abstractmethod
    def __len__(self) -> int: ...

    @abc.abstractmethod
    def __eq__(self, other: object) -> bool: ...

    @abc.abstractmethod
    def get_features_matching(
        self,
        *,
        biotype: str | tuple[str, ...] | list[str] | set[str] | None = None,
        seqid: str | None = None,
        name: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        strand: str | None = None,
        attributes: str | None = None,
        on_alignment: bool | None = None,
        allow_partial: bool = False,
    ) -> Iterator[FeatureDataType]: ...

    @abc.abstractmethod
    def get_records_matching(
        self,
        *,
        biotype: str | None = None,
        seqid: str | None = None,
        name: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        strand: str | None = None,
        attributes: str | None = None,
        on_alignment: bool | None = None,
        allow_partial: bool = False,
    ) -> Iterator[dict[str, Any]]: ...

    @abc.abstractmethod
    def get_feature_children(
        self,
        name: str,
        biotype: str | None = None,
        **kwargs: Any,
    ) -> Iterator[FeatureDataType]: ...

    @abc.abstractmethod
    def get_feature_parent(
        self,
        name: str,
        **kwargs: Any,
    ) -> Iterator[FeatureDataType]: ...

    @abc.abstractmethod
    def num_matches(
        self,
        *,
        seqid: str | None = None,
        biotype: str | tuple[str, ...] | list[str] | set[str] | None = None,
        name: str | None = None,
        strand: str | None = None,
        attributes: str | None = None,
        on_alignment: bool | None = None,
    ) -> int: ...

    def subset(
        self,
        *,
        source: PathType = ":memory:",
        biotype: str | None = None,
        seqid: str | None = None,
        name: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        strand: str | None = None,
        attributes: str | None = None,
        allow_partial: bool = False,
    ) -> Self:
        # override in subclass
        msg = f"{self.__class__.__name__!r} does not support subset()"
        raise NotImplementedError(msg)

    def add_feature(
        self,
        *,
        seqid: str,
        biotype: str,
        name: str,
        spans: list[tuple[int, int]] | NumpyIntArrayType,
        parent_id: str | None = None,
        strand: str | int | Strand | None = None,
        attributes: str | None = None,
        on_alignment: bool | None = False,
    ) -> None:
        # override in subclass
        msg = f"{self.__class__.__name__!r} does not support add_feature()"
        raise NotImplementedError(msg)

    def add_records(
        self,
        records: Iterable[dict[str, Any]],
        **kwargs: Any,
    ) -> None:
        # override in subclass
        msg = f"{self.__class__.__name__!r} does not support add_records()"
        raise NotImplementedError(msg)

    @abc.abstractmethod
    def compatible(self, other_db: Self, symmetric: bool = True) -> bool: ...

    def update(
        self,
        annot_db: Self,
        seqids: str | PySeq[str] | None = None,
        **kwargs: Any,
    ) -> None:
        # override in subclass
        msg = f"{self.__class__.__name__!r} does not support update()"
        raise NotImplementedError(msg)

    def union(self, annot_db: Self) -> Self:
        # override in subclass
        msg = f"{self.__class__.__name__!r} does not support union"
        raise NotImplementedError(msg)

    def count_distinct(
        self,
        *,
        seqid: str | bool = False,
        biotype: str | bool = False,
        name: str | bool = False,
    ) -> Table | None:
        # override in subclass
        msg = f"{self.__class__.__name__!r} does not support count_distinct()"
        raise NotImplementedError(msg)

    @abc.abstractmethod
    def close(self) -> None: ...


def _make_table_sql(
    table_name: str,
    columns: dict[Any, Any],
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
    data: dict[str, Any],
) -> tuple[str, tuple[Any, ...]]:
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
    conditions: dict[str, Any],
    allow_partial: bool = True,
) -> tuple[str, tuple[Any, ...]]:
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

    sql: list[str] = []
    vals: tuple[Any, ...] = ()
    if conditions:
        conds: list[str] = []
        cond_vals: list[Any] = []
        for col, val in conditions.items():
            # conditions are filtered for None before here, so we should add
            # an else where the op is assigned !=
            if isinstance(val, tuple | set | list):
                v = cast("Iterable[Any]", val)
                placeholders = ",".join(["?" for _ in v])
                op = "IN"
                conds.append(f"{col} {op} ({placeholders})")
                cond_vals.extend(v)
            elif val is not None:
                op = "LIKE" if isinstance(val, str) and "%" in val else "="
                conds.append(f"{col} {op} ?")
                cond_vals.append(val)
        sql.append(" AND ".join(conds))
        vals = tuple(cond_vals)

    if start is not None and stop is not None:
        if allow_partial:
            # allow matches that overlap the segment
            cond = " OR ".join(
                (
                    f"(start >= {start} AND stop <= {stop})",  # lies within the segment
                    f"(start <= {start} AND stop > {start})",  # straddles beginning of segment
                    f"(start < {stop} AND stop >= {stop})",  # straddles stop of segment
                    f"(start <= {start} AND stop >= {stop})",  # includes segment
                )
            )
        else:
            # only matches within bounds
            cond = f"start >= {start} AND stop <= {stop}"
        sql.append(f"({cond})")
    elif start is not None:
        if allow_partial:
            # any feature that overlaps region [start, ∞)
            cond = f"stop > {start}"
        else:
            # features fully within [start, ∞)
            cond = f"start >= {start}"
        sql.append(f"({cond})")
    elif stop is not None:
        if allow_partial:
            # any feature that overlaps region [0, stop)
            cond = f"start < {stop}"
        else:
            # features fully within [0, stop)
            cond = f"stop <= {stop}"
        sql.append(f"({cond})")

    return f"{' AND '.join(sql)}", vals


def _select_records_sql(
    table_name: str,
    conditions: dict[str, Any],
    columns: Iterable[str] | None = None,
    allow_partial: bool = True,
) -> tuple[str, tuple[Any, ...] | None]:
    """create SQL select statement and values

    Parameters
    ----------
    table_name
        containing the data to be selected from
    columns
        values to select
    conditions
        the WHERE conditions
    allow_partial
        if False, only records within start, stop are included. If True,
        all records that overlap the segment defined by start, stop are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """

    where, vals = _matching_conditions(
        conditions=conditions,
        allow_partial=allow_partial,
    )
    columns_str = f"{', '.join(columns)}" if columns else "*"
    sql = f"SELECT {columns_str} FROM {table_name}"
    if not where:
        return sql, None

    sql = f"{sql} WHERE {where};"
    return sql, vals


def _count_records_sql(
    table_name: str,
    conditions: dict[str, Any],
    allow_partial: bool = True,
) -> tuple[str, tuple[Any, ...] | None]:
    """create SQL count statement and values

    Parameters
    ----------
    table_name
        containing the data to be selected from
    columns
        values to select
    conditions
        the WHERE conditions
    start, stop
        select records whose (start, stop) values lie between start and stop,
        or overlap them if (allow_partial is True)
    allow_partial
        if False, only records within start, stop are included. If True,
        all records that overlap the segment defined by start, stop are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """

    where, vals = _matching_conditions(
        conditions=conditions,
        allow_partial=allow_partial,
    )
    sql = f"SELECT COUNT(*) FROM {table_name}"
    if not where:
        return sql, None

    sql = f"{sql} WHERE {where};"
    return sql, vals


# Tables that are part of the schema and should be ignored in compatibility checks
_SCHEMA_TABLES = frozenset(
    {
        "seqids",
        "biotypes",
        "sources",
        "feature_spans",
        "feature_hierarchy",
        "schema_version",
    }
)


def _make_db_connection(
    source: PathType,
    table_names: tuple[str, ...] | None = None,
) -> sqlite3.Connection:
    """Create a sqlite3 connection with standard settings.

    Parameters
    ----------
    source
        Path to database file or ":memory:" for in-memory database.
    table_names
        Optional tuple of table names to check for legacy 'end' column rename.
        If provided, renames 'end' to 'stop' in these tables.

    Returns
    -------
    sqlite3.Connection configured with PARSE_DECLTYPES and Row factory.
    """
    conn = sqlite3.connect(
        source,
        detect_types=sqlite3.PARSE_DECLTYPES,
        check_same_thread=False,
    )
    conn.row_factory = sqlite3.Row

    # Handle legacy databases with 'end' column instead of 'stop'
    if table_names is not None:
        for table_name in table_names:
            _rename_column_if_exists(conn, table_name, "end", "stop")

    return conn


def _get_name_to_rowid_batch(
    conn: sqlite3.Connection,
    table_name: str,
    names: set[str],
) -> dict[str, int]:
    """Get rowids for multiple features by name in a single query.

    Parameters
    ----------
    conn
        Database connection
    table_name
        Name of table to query
    names
        Set of feature names to look up

    Returns
    -------
    Dict mapping feature name to rowid
    """
    if not names:
        return {}
    placeholders = ",".join("?" * len(names))
    names_list = list(names)
    sql = f"SELECT name, rowid FROM {table_name} WHERE name IN ({placeholders})"
    cursor = conn.execute(sql, tuple(names_list))
    return {row[0]: row[1] for row in cursor.fetchall()}


def _populate_hierarchy_table(
    conn: sqlite3.Connection,
    table_name: str,
    entries: list[tuple[int, str]],
) -> None:
    """Populate the feature_hierarchy table from child-parent entries.

    Parameters
    ----------
    conn
        Database connection with feature_hierarchy table
    table_name
        Name of feature table (gff, gb, user)
    entries
        List of (child_rowid, parent_id_string) tuples where parent_id_string
        may contain comma-separated parent names
    """
    if not entries:
        return

    # Collect all parent names for batch lookup
    parent_names: set[str] = set()
    for _, parent_id in entries:
        for p in parent_id.replace(" ", "").split(","):
            if p := p.strip():
                parent_names.add(p)

    # Batch lookup parent names to rowids
    name_to_rowid = _get_name_to_rowid_batch(conn, table_name, parent_names)

    # Build hierarchy rows
    hierarchy_rows = []
    for child_rowid, parent_id in entries:
        for p in parent_id.replace(" ", "").split(","):
            p = p.strip()
            if p and p in name_to_rowid:
                hierarchy_rows.append((child_rowid, name_to_rowid[p], table_name))

    if hierarchy_rows:
        conn.executemany(
            "INSERT OR IGNORE INTO feature_hierarchy (child_id, parent_id, table_name) VALUES (?, ?, ?)",
            hierarchy_rows,
        )
        conn.commit()


def _migrate_to_normalized_schema(
    conn: sqlite3.Connection,
    table_names: list[str],
    lookup_cache: LookupTableCache,
) -> None:
    """Migrate old TEXT columns (seqid, biotype, source) to new *_id INTEGER columns.

    Parameters
    ----------
    conn
        Database connection
    table_names
        List of table names to migrate (gff, gb, user)
    lookup_cache
        LookupTableCache instance for getting/creating ID mappings

    Notes
    -----
    This function handles the schema migration from old-format databases that used
    TEXT columns directly to the new normalized format using lookup tables.
    """
    # Mapping of old column names to new column names
    column_mappings = [
        ("seqid", "seqid_id", lookup_cache.get_seqid_id),
        ("biotype", "biotype_id", lookup_cache.get_biotype_id),
        ("source", "source_id", lookup_cache.get_source_id),
    ]

    for table_name in table_names:
        # Get current table columns
        table_info = conn.execute(f"PRAGMA table_info({table_name})").fetchall()
        existing_columns = {row["name"] for row in table_info}

        # Check which old columns need to be migrated
        columns_to_migrate = []
        for old_col, new_col, get_id_func in column_mappings:
            if old_col in existing_columns and new_col not in existing_columns:
                columns_to_migrate.append((old_col, new_col, get_id_func))

        if not columns_to_migrate:
            continue

        # Add new INTEGER columns
        for old_col, new_col, _ in columns_to_migrate:
            conn.execute(
                f"ALTER TABLE {table_name} ADD COLUMN {new_col} INTEGER DEFAULT 0"
            )
        conn.commit()

        # Get all unique values for each column to migrate and populate lookup tables
        # Then update all rows with the new IDs
        for old_col, new_col, get_id_func in columns_to_migrate:
            # Get all distinct values from the old column
            cursor = conn.execute(
                f"SELECT DISTINCT {old_col} FROM {table_name} WHERE {old_col} IS NOT NULL AND {old_col} != ''"
            )
            distinct_values = [row[old_col] for row in cursor.fetchall()]

            # Populate lookup table and get ID mappings
            value_to_id = {}
            for value in distinct_values:
                if value:
                    value_to_id[value] = get_id_func(value)

            # Update rows with new IDs
            for value, id_val in value_to_id.items():
                conn.execute(
                    f"UPDATE {table_name} SET {new_col} = ? WHERE {old_col} = ?",
                    (id_val, value),
                )

        conn.commit()

        # Now we need to drop the old columns
        # SQLite doesn't support DROP COLUMN directly in older versions,
        # so we need to recreate the table without those columns
        # Get the full schema to recreate the table
        table_info = conn.execute(f"PRAGMA table_info({table_name})").fetchall()
        old_column_names = {old_col for old_col, _, _ in columns_to_migrate}

        # Build column list for new table (excluding old TEXT columns)
        new_columns = []
        column_defs = []
        for col in table_info:
            if col["name"] not in old_column_names:
                new_columns.append(col["name"])
                col_type = col["type"]
                col_def = f'"{col["name"]}" {col_type}'
                if col["notnull"]:
                    col_def += " NOT NULL"
                if col["dflt_value"] is not None:
                    col_def += f" DEFAULT {col['dflt_value']}"
                column_defs.append(col_def)

        # Create temporary table with new schema
        temp_table = f"{table_name}_migration_temp"
        create_sql = f"CREATE TABLE {temp_table} ({', '.join(column_defs)})"
        conn.execute(create_sql)

        # Copy data to temporary table
        columns_str = ", ".join(f'"{c}"' for c in new_columns)
        conn.execute(
            f"INSERT INTO {temp_table} ({columns_str}) SELECT {columns_str} FROM {table_name}"
        )

        # Drop old table and rename temporary table
        conn.execute(f"DROP TABLE {table_name}")
        conn.execute(f"ALTER TABLE {temp_table} RENAME TO {table_name}")

        conn.commit()


def _compatible_schema(
    db: sqlite3.Connection, schema: dict[str, set[tuple[str, str]]]
) -> bool:
    """ensures the db instance is compatible with schema"""
    for table in db.execute(
        "SELECT name FROM sqlite_master WHERE type='table'",
    ).fetchall():
        table = table["name"]  # noqa: PLW2901
        # Skip auxiliary tables
        if table in _SCHEMA_TABLES:
            continue
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
    db: sqlite3.Connection,
    table_name: str,
    old_column: str,
    new_column: str,
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
    _table_names: ClassVar[tuple[str, ...]]
    _suffix: ClassVar[str]

    _user_schema: ClassVar = {
        "biotype_id": "INTEGER",
        "seqid_id": "INTEGER",
        "name": "TEXT",
        "parent_id": "TEXT",
        "start": "INTEGER",
        "stop": "INTEGER",
        "strand": "INTEGER",
        "spans": "array",
        "attributes": "TEXT",
        "on_alignment": "INT",
    }
    # args to exclude from serialisation init
    _exclude_init = "db", "data"

    def __new__(cls, *args: Any, **kwargs: Any) -> Self:
        obj = object.__new__(cls)
        init_sig = inspect.signature(cls.__init__)
        bargs = init_sig.bind_partial(cls, *args, **kwargs)
        bargs.apply_defaults()
        init_vals = bargs.arguments
        for param in cls._exclude_init:
            init_vals.pop(param)
        init_vals.pop("self", None)

        obj._serialisable = init_vals  # noqa: SLF001
        return obj

    def __init__(self, *args: Any, **kwargs: Any) -> None:  # pragma: no cover
        self._serialisable: dict[str, Any]
        self.source: PathType
        self._schema_version: int = 1
        self._lookup_cache: LookupTableCache | None = None

    def __repr__(self) -> str:
        name = self.__class__.__name__
        total_records = len(self)
        args = ", ".join(
            f"{k}={repr(v) if isinstance(v, str) else v}"
            for k, v in self._serialisable.items()
            if k != "data"
        )
        return f"{name}({args}, total_records={total_records})"

    def __deepcopy__(
        self, memodict: dict[int, SqliteAnnotationDbMixin] | None = None
    ) -> Self:
        memodict = memodict or {}
        try:
            new = self.__class__(source=self.source)
            new.db.deserialize(self._db.serialize())  # type: ignore[union-attr]
        except AttributeError:
            # if the db is not serialisable, we use the rich dict
            # representation to create a new instance
            rd = self.to_rich_dict()
            new = type(self).from_dict(rd)
        return new

    def __getstate__(self) -> dict[str, Any]:
        try:
            result: dict[str, Any] = {
                "data": self._db.serialize(),  # type: ignore[union-attr]
                "source": self.source,
            }
        except AttributeError:
            # if the db is not serialisable, we use the rich dict
            # representation to create a new instance
            result = self.to_rich_dict()

        return result

    def __setstate__(self, state: dict[str, Any]) -> Self:
        if "type" in state:
            data = type(self).from_dict(state)
            self.__dict__.update(data.__dict__)
        else:
            new = self.__class__(source=state.pop("source", None))
            new._db.deserialize(state["data"])  # type: ignore[union-attr]  # noqa: SLF001
            self.__dict__.update(new.__dict__)

        return self

    def __len__(self) -> int:
        return self.num_matches()

    def __eq__(self, other: object) -> bool:
        return isinstance(other, self.__class__) and other.db is self.db

    @property
    def table_names(self) -> tuple[str, ...]:
        return self._table_names

    def _setup_db(
        self, db: SqliteAnnotationDbMixin | sqlite3.Connection | None
    ) -> None:
        """initialises the db, using the db passed to the constructor"""
        if isinstance(db, self.__class__):
            self._db: sqlite3.Connection | None = db.db
            self._schema_version = db._schema_version
            self._lookup_cache = db._lookup_cache
            return

        if isinstance(db, sqlite3.Connection):
            schema: dict[str, set[tuple[str, str]]] = {}
            for table_name in self.table_names:
                attr = getattr(self, f"_{table_name}_schema")
                schema[table_name] = set(attr.items())

            if not _compatible_schema(db, schema):
                msg = "incompatible schema"
                raise TypeError(msg)

            self._db = db
            self._init_tables()
            return

        if db and not self.compatible(cast("AnnotationDbABC", db)):
            msg = f"cannot initialise annotation db from {type(db)}"
            raise TypeError(msg)

        self._init_tables()

        if db and len(db):
            # update self with data from other
            self.update(cast("AnnotationDbABC", db))

    def _init_tables(self) -> None:
        # bit of magic
        # assumes schema attributes named as `_<table name>_schema`
        for table_name in self.table_names:
            attr = getattr(self, f"_{table_name}_schema")
            sql = _make_table_sql(table_name, attr)
            self._execute_sql(sql)

        self._schema_version = _get_schema_version(self.db)

        if self._schema_version >= SCHEMA_VERSION:
            # Already at current version, just set up caches
            self._setup_caches()
            return

        # Database needs upgrade - create the schema v2 tables
        self._create_v2_tables()
        _set_schema_version(self.db, SCHEMA_VERSION)
        self._schema_version = SCHEMA_VERSION
        self._setup_caches()

    def _create_v2_tables(self) -> None:
        # Create lookup tables
        for sql in _LOOKUP_TABLES_SQL.values():
            self._execute_sql(sql)

        # Create feature_spans table for multi-span features
        self._execute_sql(_FEATURE_SPANS_SQL)

        # Create feature_hierarchy table
        self._execute_sql(_FEATURE_HIERARCHY_SQL)

    def _setup_caches(self) -> None:
        """Set up caches for schema lookups."""
        if self._schema_version >= SCHEMA_VERSION:
            self._lookup_cache = LookupTableCache(self.db)
            self._lookup_cache.load_all()
        else:
            self._lookup_cache = None

    @property
    def schema_version(self) -> int:
        """Return the schema version of this database."""
        return self._schema_version

    @property
    def db(self) -> sqlite3.Connection:
        if self._db is None:
            # TODO: gah understand serialisation issue
            # the check_same_thread=False is required for multiprocess, even
            # when the db is empty (tests fail). This  appears unrelated to
            # our code, and does not affect pickling/unpickling on the same
            # thread
            self._db = _make_db_connection(self.source, table_names=self.table_names)

            # Detect schema version for existing databases
            self._schema_version = _get_schema_version(self._db)
            self._setup_caches()

        return self._db

    def _execute_sql(
        self,
        cmnd: str,
        values: tuple[Any, ...] | None = None,
    ) -> sqlite3.Cursor:
        with self.db:
            # context manager ensures safe transactions
            cursor = self.db.cursor()
            cursor.execute(cmnd, values or [])
            return cursor

    def _add_to_hierarchy(
        self,
        table_name: str,
        child_id: int,
        parent_id: str | None,
    ) -> None:
        """Add parent-child relationships to the hierarchy table."""
        if not self._can_use_hierarchy() or parent_id is None:
            return

        parent_rowids = self._resolve_parent_rowids(table_name, parent_id)
        for parent_rowid in parent_rowids:
            self._insert_hierarchy_entry(child_id, parent_rowid, table_name)

    def _can_use_hierarchy(self) -> bool:
        """Check if hierarchy table operations are available."""
        return self._schema_version >= SCHEMA_VERSION

    def _get_lookup_ids(
        self, seqid: str | None, biotype: str | None, source: str | None = None
    ) -> tuple[int, int, int]:
        """Get lookup table IDs for seqid, biotype, and source."""
        if self._lookup_cache is None:
            return 0, 0, 0
        seqid_id = self._lookup_cache.get_seqid_id(seqid) if seqid else 0
        biotype_id = self._lookup_cache.get_biotype_id(biotype) if biotype else 0
        source_id = self._lookup_cache.get_source_id(source) if source else 0
        return seqid_id, biotype_id, source_id

    def _translate_record_to_ids(self, record: dict[str, Any]) -> dict[str, Any]:
        """Convert seqid/biotype/source strings to IDs for insertion."""
        if self._lookup_cache is None:
            return record

        result = dict(record)
        if "seqid" in result:
            seqid = result.pop("seqid")
            result["seqid_id"] = self._lookup_cache.get_seqid_id(seqid) if seqid else 0
        if "biotype" in result:
            biotype = result.pop("biotype")
            result["biotype_id"] = (
                self._lookup_cache.get_biotype_id(biotype) if biotype else 0
            )
        if "source" in result:
            source = result.pop("source")
            result["source_id"] = (
                self._lookup_cache.get_source_id(source) if source else 0
            )
        return result

    def _translate_record_from_ids(self, record: dict[str, Any]) -> dict[str, Any]:
        """Convert seqid_id/biotype_id/source_id to strings for query results."""
        if self._lookup_cache is None:
            return record

        result = dict(record)
        if "seqid_id" in result:
            seqid_id = result.pop("seqid_id")
            result["seqid"] = (
                self._lookup_cache.get_seqid_name(seqid_id) if seqid_id else ""
            )
        if "biotype_id" in result:
            biotype_id = result.pop("biotype_id")
            result["biotype"] = (
                self._lookup_cache.get_biotype_name(biotype_id) if biotype_id else ""
            )
        if "source_id" in result:
            source_id = result.pop("source_id")
            result["source"] = (
                self._lookup_cache.get_source_name(source_id) if source_id else ""
            )
        return result

    def _translate_query_conditions(self, conditions: dict[str, Any]) -> dict[str, Any]:
        """Convert seqid/biotype/source strings to IDs for query conditions."""
        if self._lookup_cache is None:
            return conditions

        result = dict(conditions)
        if "seqid" in result:
            seqid = result.pop("seqid")
            if seqid is not None:
                if isinstance(seqid, str):
                    result["seqid_id"] = self._lookup_cache.get_seqid_id(seqid)
                else:
                    # Multiple seqids (tuple, list, set)
                    result["seqid_id"] = tuple(
                        self._lookup_cache.get_seqid_id(s) for s in seqid
                    )
        if "biotype" in result:
            biotype = result.pop("biotype")
            if biotype is not None:
                if isinstance(biotype, str):
                    result["biotype_id"] = self._lookup_cache.get_biotype_id(biotype)
                else:
                    # Multiple biotypes (tuple, list, set)
                    result["biotype_id"] = tuple(
                        self._lookup_cache.get_biotype_id(b) for b in biotype
                    )
        if "source" in result:
            source = result.pop("source")
            if source is not None:
                result["source_id"] = self._lookup_cache.get_source_id(source)
        return result

    def _get_feature_rowid(
        self,
        table_name: str,
        name: str,
        biotype: str | None = None,
    ) -> int | None:
        """Get rowid for a feature by name (and optionally biotype)."""
        if biotype and self._lookup_cache is not None:
            biotype_id = self._lookup_cache.get_biotype_id(biotype)
            sql = f"SELECT rowid FROM {table_name} WHERE name = ? AND biotype_id = ?"
            cursor = self._execute_sql(sql, (name, biotype_id))
        else:
            sql = f"SELECT rowid FROM {table_name} WHERE name = ?"
            cursor = self._execute_sql(sql, (name,))
        row = cursor.fetchone()
        return row[0] if row else None

    def _get_feature_rowids_batch(
        self,
        table_name: str,
        names: set[str],
    ) -> dict[str, int]:
        """Get rowids for multiple features by name in a single query."""
        if not names:
            return {}
        # Use a single query with IN clause for efficiency
        placeholders = ",".join("?" * len(names))
        names_list = list(names)
        sql = f"SELECT name, rowid FROM {table_name} WHERE name IN ({placeholders})"
        cursor = self._execute_sql(sql, tuple(names_list))
        return {row[0]: row[1] for row in cursor.fetchall()}

    def _resolve_parent_rowids(
        self,
        table_name: str,
        parent_id: str,
        name_to_rowid: dict[str, int] | None = None,
    ) -> list[int]:
        """Parse parent_id string and resolve to list of rowids.

        If name_to_rowid is provided, use it as a cache to avoid SQL queries.
        """
        parent_names = [
            p.strip() for p in parent_id.replace(" ", "").split(",") if p.strip()
        ]
        rowids = []
        for parent_name in parent_names:
            if name_to_rowid is not None:
                rowid = name_to_rowid.get(parent_name)
            else:
                rowid = self._get_feature_rowid(table_name, parent_name)
            if rowid is not None:
                rowids.append(rowid)
        return rowids

    def _insert_hierarchy_entry(
        self,
        child_id: int,
        parent_id: int,
        table_name: str,
    ) -> None:
        """Insert a single entry into the hierarchy table."""
        try:
            self._execute_sql(
                "INSERT OR IGNORE INTO feature_hierarchy (child_id, parent_id, table_name) VALUES (?, ?, ?)",
                (child_id, parent_id, table_name),
            )
        except sqlite3.IntegrityError:
            pass  # Already exists

    def add_feature(
        self,
        *,
        seqid: str,
        biotype: str,
        name: str,
        spans: list[tuple[int, int]] | NumpyIntArrayType,
        parent_id: str | None = None,
        strand: str | int | Strand | None = None,
        attributes: str | None = None,
        on_alignment: bool | None = False,
    ) -> None:
        """adds a record to user table

        Parameters
        ----------
        seqid
            name of the sequence feature resides on
        biotype
            biological type of the record
        name
            the name of a record, an identifier
        spans
            this will be sorted
        strand
            either +, -. Defaults to '+'
        attributes
            additional attributes as a string
        on_alignment
            whether the annotation is an alignment annotation
        """
        spans = numpy.array(sorted(sorted(coords) for coords in spans), dtype=int)
        # the start, stop variables are added into the record via the loop
        # on local variables
        start = int(spans.min())
        stop = int(spans.max())
        # we define record as all defined variables from local name space,
        # excluding "self"
        local_vars = locals()
        record = {k: v for k, v in local_vars.items() if k != "self" and v is not None}
        if "strand" in record:
            record["strand"] = Strand.from_value(strand).value
        # Translate seqid/biotype to IDs
        record = self._translate_record_to_ids(record)
        sql, values = _add_record_sql("user", record)
        cursor = self._execute_sql(sql, values=values)
        feature_id = cursor.lastrowid

        # Add to hierarchy table if available
        if feature_id is not None:
            self._add_to_hierarchy("user", feature_id, parent_id)

    def _translate_columns_to_ids(
        self, columns: Iterable[str] | None
    ) -> tuple[str, ...] | None:
        """Translate column names from strings to ID columns for normalized schema."""
        if columns is None:
            return None
        mapping = {"seqid": "seqid_id", "biotype": "biotype_id", "source": "source_id"}
        return tuple(mapping.get(c, c) for c in columns)

    def _get_records_matching(
        self,
        table_name: str,
        **kwargs: Any,
    ) -> Iterator[sqlite3.Row]:
        """return all fields"""
        if kwargs.get("attributes") and "%%" not in kwargs["attributes"]:
            kwargs["attributes"] = f"%{kwargs['attributes']}%"
        columns = kwargs.pop("columns", None)
        allow_partial = kwargs.pop("allow_partial", False)

        # Translate query conditions and column names for normalized schema
        if self._lookup_cache is not None:
            kwargs = self._translate_query_conditions(kwargs)
            columns = self._translate_columns_to_ids(columns)

        sql, vals = _select_records_sql(
            table_name=table_name,
            conditions=kwargs,
            columns=columns,
            allow_partial=allow_partial,
        )
        with contextlib.suppress(sqlite3.ProgrammingError):
            # garbage collection issue
            yield from self._execute_sql(sql, values=vals)

    # Return type is instead larger than FeatureDataType
    def _get_feature_by_id(
        self,
        table_name: str,
        columns: Iterable[str] | None,
        column: str,
        name: str,
        biotype: str | None = None,
        allow_partial: bool = False,
        **kwargs: Any,  # noqa: ARG002
    ) -> Iterator[FeatureDataType]:
        # we return the parent_id because `get_feature_parent()` requires it
        # Translate conditions and columns for normalized schema
        conditions: dict[str, Any] = {column: name, "biotype": biotype}
        if self._lookup_cache is not None:
            conditions = self._translate_query_conditions(conditions)
            columns = self._translate_columns_to_ids(columns)
        sql, vals = _select_records_sql(
            table_name=table_name,
            conditions=conditions,
            columns=columns,
            allow_partial=allow_partial,
        )
        for result in self._execute_sql(sql, values=vals):
            result = cast(  # noqa: PLW2901
                "FeatureDataType", dict(zip(result.keys(), result, strict=False))
            )
            # Translate IDs back to strings
            result = self._translate_record_from_ids(result)  # type: ignore[arg-type]
            result["on_alignment"] = result.get("on_alignment")
            result["spans"] = [
                cast("tuple[int, int]", tuple(c)) for c in result["spans"]
            ]
            yield cast("FeatureDataType", result)

    def _has_hierarchy_table(self) -> bool:
        """Check if the feature_hierarchy table exists and has data."""
        try:
            cursor = self.db.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name='feature_hierarchy'"
            )
            return cursor.fetchone() is not None
        except sqlite3.Error:
            return False

    def _get_children_via_hierarchy(
        self,
        table_name: str,
        parent_name: str,
        columns: tuple[str, ...],
        biotype: str | None = None,
    ) -> Iterator[FeatureDataType]:
        # First, find the parent's rowid
        cursor = self._execute_sql(
            f"SELECT rowid FROM {table_name} WHERE name = ?",
            (parent_name,),
        )
        row = cursor.fetchone()
        if not row:
            return
        parent_rowid = row[0]

        # Translate column names to ID columns for normalized schema
        translated_columns = self._translate_columns_to_ids(columns) or columns
        columns_str = ", ".join(f"f.{c}" for c in translated_columns)
        sql = f"""
            SELECT {columns_str} FROM {table_name} f
            INNER JOIN feature_hierarchy h ON f.rowid = h.child_id
            WHERE h.parent_id = ? AND h.table_name = ?
        """
        vals: list[Any] = [parent_rowid, table_name]

        if biotype is not None and self._lookup_cache is not None:
            biotype_id = self._lookup_cache.get_biotype_id(biotype)
            sql += " AND f.biotype_id = ?"
            vals.append(biotype_id)

        for result in self._execute_sql(sql, values=tuple(vals)):
            res = dict(zip(result.keys(), result, strict=False))
            # Translate IDs back to strings
            res = self._translate_record_from_ids(res)  # type: ignore[arg-type]
            res["on_alignment"] = res.get("on_alignment")
            res["spans"] = [cast("tuple[int, int]", tuple(c)) for c in res["spans"]]
            # Remove parent_id if present
            res.pop("parent_id", None)  # type: ignore[typeddict-item]
            yield cast("FeatureDataType", res)

    def _get_parents_via_hierarchy(
        self,
        table_name: str,
        child_name: str,
        columns: tuple[str, ...],
    ) -> Iterator[FeatureDataType]:
        # First, find the child's rowid
        cursor = self._execute_sql(
            f"SELECT rowid FROM {table_name} WHERE name = ?",
            (child_name,),
        )
        row = cursor.fetchone()
        if not row:
            return
        child_rowid = row[0]

        # Translate column names to ID columns for normalized schema
        translated_columns = self._translate_columns_to_ids(columns) or columns
        columns_str = ", ".join(f"f.{c}" for c in translated_columns)
        sql = f"""
            SELECT {columns_str} FROM {table_name} f
            INNER JOIN feature_hierarchy h ON f.rowid = h.parent_id
            WHERE h.child_id = ? AND h.table_name = ?
        """
        vals = child_rowid, table_name

        for result in self._execute_sql(sql, values=vals):
            res = dict(zip(result.keys(), result, strict=False))
            # Translate IDs back to strings
            res = self._translate_record_from_ids(res)  # type: ignore[arg-type]
            res["on_alignment"] = res.get("on_alignment")
            res["spans"] = [cast("tuple[int, int]", tuple(c)) for c in res["spans"]]
            # Remove parent_id if present
            res.pop("parent_id", None)  # type: ignore[typeddict-item]
            yield cast("FeatureDataType", res)

    def _old_get_feature_children(
        self,
        name: str,
        biotype: str | None = None,
        **kwargs: Any,
    ) -> Iterator[FeatureDataType]:  # pragma: no cover
        for table_name in self.table_names:
            columns: tuple[str, ...] = (
                "seqid",
                "biotype",
                "spans",
                "strand",
                "name",
                "parent_id",
            )
            if table_name == "user":
                columns += ("on_alignment",)

            # Fall back to LIKE-based search
            for result in self._get_feature_by_id(
                table_name=table_name,
                columns=columns,
                column="parent_id",
                name=f"%{name}%",
                biotype=biotype,
                **kwargs,
            ):
                # remove invalid field for the FeatureDataType #
                result.pop("parent_id")  # type: ignore[typeddict-item]
                yield result

    def get_feature_children(
        self,
        name: str,
        biotype: str | None = None,
        **kwargs: Any,
    ) -> Iterator[FeatureDataType]:
        """yields children of name"""
        if not self._can_use_hierarchy():
            yield from self._old_get_feature_children(
                name=name, biotype=biotype, **kwargs
            )
            return

        for table_name in self.table_names:
            columns: tuple[str, ...] = (
                "seqid",
                "biotype",
                "spans",
                "strand",
                "name",
                "parent_id",
            )
            if table_name == "user":
                columns += ("on_alignment",)

            yield from self._get_children_via_hierarchy(
                table_name=table_name,
                parent_name=name,
                columns=columns,
                biotype=biotype,
            )

    def _old_get_feature_parent(
        self,
        name: str,
        **kwargs: Any,  # noqa: ARG002
    ) -> Iterator[FeatureDataType]:  # pragma: no cover
        """yields parents of name"""
        for table_name in self.table_names:
            columns: tuple[str, ...] = (
                "seqid",
                "biotype",
                "spans",
                "strand",
                "name",
                "parent_id",
            )
            if table_name == "user":
                columns += ("on_alignment",)

            # Older version uses LIKE based search
            for result in self._get_feature_by_id(
                table_name=table_name,
                columns=columns,
                column="name",
                name=f"%{name}%",
            ):
                # multiple values for parent means this is better expressed
                # as an OR clause
                # TODO modify the conditional SQL generation
                # When and when isn't parent_id in FeatureDataType? Should it be in the TypedDict? Or should _get_feature_by_id have a different return type?
                if not result["parent_id"]:  # type: ignore[typeddict-item]
                    return

                for name in (  # noqa: PLR1704
                    cast("str", result["parent_id"]).replace(" ", "").split(",")  # type: ignore[typeddict-item]
                ):
                    if parent := list(
                        self._get_feature_by_id(
                            table_name=table_name,
                            columns=columns,
                            column="name",
                            name=name,
                        ),
                    ):
                        par = parent[0]
                        par.pop(
                            "parent_id",  # type: ignore[typeddict-item]
                        )  # remove invalid field for the FeatureDataType
                        yield par

    def get_feature_parent(
        self,
        name: str,
        **kwargs: Any,
    ) -> Iterator[FeatureDataType]:
        """yields parents of name"""
        if not self._can_use_hierarchy():
            yield from self._old_get_feature_parent(name=name, **kwargs)
            return

        for table_name in self.table_names:
            columns: tuple[str, ...] = (
                "seqid",
                "biotype",
                "spans",
                "strand",
                "name",
                "parent_id",
            )
            if table_name == "user":
                columns += ("on_alignment",)

            yield from self._get_parents_via_hierarchy(
                table_name=table_name,
                child_name=name,
                columns=columns,
            )

    def get_records_matching(
        self,
        *,
        biotype: str | None = None,
        seqid: str | None = None,
        name: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        strand: str | None = None,
        attributes: str | None = None,
        on_alignment: bool | None = None,
        allow_partial: bool = False,
    ) -> Iterator[dict[str, Any]]:
        """return all fields for matching records"""
        # a record is Everything, a Feature is a subset
        # we define query as all defined variables from local name space,
        # excluding "self" and kwargs at default values
        local_vars = locals()
        kwargs = {k: v for k, v in local_vars.items() if k != "self" and v is not None}
        if "strand" in kwargs:
            kwargs["strand"] = Strand.from_value(strand).value

        # alignment features are created by the user specific
        table_names = ["user"] if on_alignment else self.table_names
        for table_name in table_names:
            for result in self._get_records_matching(table_name, **kwargs):
                res = dict(zip(result.keys(), result, strict=False))
                yield self._translate_record_from_ids(res)

    def get_features_matching(
        self,
        *,
        biotype: str | tuple[str, ...] | list[str] | set[str] | None = None,
        seqid: str | None = None,
        name: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        strand: str | None = None,
        attributes: str | None = None,
        on_alignment: bool | None = None,
        allow_partial: bool = False,
    ) -> Iterator[FeatureDataType]:
        # returns essential values to create a Feature
        # we define query as all defined variables from local name space,
        # excluding "self" and kwargs with a default value of None
        local_vars = locals()
        kwargs = {k: v for k, v in local_vars.items() if k != "self" and v is not None}
        if "strand" in kwargs:
            kwargs["strand"] = Strand.from_value(strand).value

        # alignment features are created by the user specific
        table_names = ["user"] if on_alignment else self.table_names
        for table_name in table_names:
            columns: tuple[str, ...] = ("seqid", "biotype", "spans", "strand", "name")
            query_args = {**kwargs}

            if table_name == "user":
                columns += ("on_alignment",)
            else:
                query_args.pop("on_alignment", None)

            for result in self._get_records_matching(
                table_name=table_name,
                columns=columns,
                **query_args,
            ):
                res = dict(zip(result.keys(), result, strict=False))
                # Translate IDs back to strings
                res = self._translate_record_from_ids(res)  # type: ignore[arg-type]
                res["on_alignment"] = res.get("on_alignment")
                res["spans"] = [cast("tuple[int, int]", tuple(c)) for c in res["spans"]]
                yield cast("FeatureDataType", res)

    def num_matches(
        self,
        *,
        seqid: str | None = None,
        biotype: str | tuple[str, ...] | list[str] | set[str] | None = None,
        name: str | None = None,
        strand: str | None = None,
        attributes: str | None = None,
        on_alignment: bool | None = None,
    ) -> int:
        """return the number of records matching condition"""
        local_vars = locals()
        kwargs = {k: v for k, v in local_vars.items() if k != "self" and v is not None}
        if "strand" in kwargs:
            kwargs["strand"] = Strand.from_value(kwargs["strand"]).value
        # Translate conditions for normalized schema
        if self._lookup_cache is not None:
            kwargs = self._translate_query_conditions(kwargs)
        num = 0
        for table_name in self.table_names:
            sql, values = _count_records_sql(table_name, conditions=kwargs)
            result = next(iter(self._execute_sql(sql, values=values).fetchone()))
            num += result
        return num

    StrOrBool = str | bool

    def count_distinct(
        self,
        *,
        seqid: StrOrBool = False,
        biotype: StrOrBool = False,
        name: StrOrBool = False,
    ) -> Table | None:
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

        >>> counts_table = db.count_distinct(seqid=True, biotype="gene", name=True)
        """
        from cogent3.core.table import Table

        columns = {k for k, v in locals().items() if v is True}
        if not columns:
            return None

        constraints = {k: v for k, v in locals().items() if isinstance(v, str)}

        header = list(columns)
        data: list[tuple[Any, ...]] = []

        for table in self.table_names:
            # Build SELECT and GROUP BY with proper joins for normalized columns
            if self._lookup_cache is not None:
                # Use joins for normalized columns
                select_parts = []
                group_by_parts = []
                joins = []
                for col in header:
                    if col == "seqid":
                        select_parts.append("s.name as seqid")
                        group_by_parts.append("t.seqid_id")
                        joins.append("LEFT JOIN seqids s ON t.seqid_id = s.id")
                    elif col == "biotype":
                        select_parts.append("b.name as biotype")
                        group_by_parts.append("t.biotype_id")
                        joins.append("LEFT JOIN biotypes b ON t.biotype_id = b.id")
                    else:
                        select_parts.append(f"t.{col}")
                        group_by_parts.append(f"t.{col}")

                if translated_constraints := self._translate_query_conditions(
                    constraints
                ):
                    # Prefix with t. for table alias
                    prefixed_constraints = {
                        f"t.{k}": v for k, v in translated_constraints.items()
                    }
                    where_clause, values = _matching_conditions(prefixed_constraints)
                    where_clause = f"WHERE {where_clause}"
                else:
                    where_clause = ""
                    values = ()

                join_clause = " ".join(dict.fromkeys(joins))  # Remove duplicates
                sql = (
                    f"SELECT {', '.join(select_parts)}, COUNT(*) as count "
                    f"FROM {table} t {join_clause} {where_clause} "
                    f"GROUP BY {', '.join(group_by_parts)};"
                )
            else:
                # No lookup cache, use original column names
                if constraints:
                    where_clause, values = _matching_conditions(constraints)
                    where_clause = f"WHERE {where_clause}"
                else:
                    where_clause = ""
                    values = ()
                sql = (
                    f"SELECT {', '.join(header)}, COUNT(*) as count FROM {table}"
                    f" {where_clause} GROUP BY {', '.join(header)};"
                )

            if result := self._execute_sql(sql, values).fetchall():
                data.extend(tuple(r) for r in result)

        return Table(
            header=[*header, "count"],
            data=data,
            column_templates={"count": lambda x: f"{x:,}"},
        )

    @property
    def describe(self) -> Table:
        """top level description of the annotation db"""
        from cogent3.core.table import Table

        data: dict[str, int] = {}
        for column in ("seqid", "biotype"):
            for table in self.table_names:
                if self._lookup_cache is not None:
                    # Use joins for normalized columns
                    if column == "seqid":
                        sql = f"SELECT s.name, COUNT(*) FROM {table} t LEFT JOIN seqids s ON t.seqid_id = s.id GROUP BY t.seqid_id;"
                    else:  # biotype
                        sql = f"SELECT b.name, COUNT(*) FROM {table} t LEFT JOIN biotypes b ON t.biotype_id = b.id GROUP BY t.biotype_id;"
                else:
                    sql = f"SELECT {column}, COUNT(*) FROM {table} GROUP BY {column};"
                if result := self._execute_sql(sql).fetchall():
                    counts = dict(tuple(r) for r in result)
                    for distinct, count in counts.items():
                        key = f"{column}({distinct!r})"
                        data[key] = data.get(key, 0) + count

        row_counts: list[int] = []
        for table in self.table_names:
            result = self._execute_sql(
                f"SELECT COUNT(*) as count FROM {table}",
            ).fetchone()
            row_counts.append(result["count"])

        return Table(
            data={
                "": list(data.keys()) + [f"num_rows({t!r})" for t in self.table_names],
                "count": list(data.values()) + row_counts,
            },
            column_templates={"count": lambda x: f"{x:,}"},
        )

    def biotype_counts(self) -> dict[str, int]:
        """return counts of biological types across all tables and seqids"""
        counts: dict[str, int] = collections.Counter()
        for table in self.table_names:
            if self._lookup_cache is not None:
                sql = f"SELECT b.name as biotype FROM {table} t LEFT JOIN biotypes b ON t.biotype_id = b.id"
            else:
                sql = f"SELECT biotype FROM {table}"
            if result := self._execute_sql(sql).fetchall():
                counts.update(v for r in result if (v := r["biotype"]))
        return counts

    def to_rich_dict(self) -> dict[str, Any]:
        """returns a dict suitable for json serialisation"""
        result: dict[str, Any] = {
            "type": get_object_provenance(self),
            "version": __version__,
            "tables": {},
            "init_args": {**self._serialisable},
        }
        tables = result["tables"]
        for table_name in self.table_names:
            table_data: list[dict[str, Any]] = []
            for record in self._get_records_matching(table_name):
                store = {
                    k: v
                    for k, v in zip(record.keys(), record, strict=False)
                    if v is not None
                }
                # Translate IDs back to strings for serialization
                store = self._translate_record_from_ids(store)
                store["spans"] = store["spans"].tolist()
                table_data.append(store)
            tables[table_name] = table_data
        return result

    def compatible(self, other_db: AnnotationDbABC, symmetric: bool = True) -> bool:
        """checks whether table_names are compatible

        Parameters
        ----------
        other_db
            the other annotation db instance
        symmetric
            checks only that tables of other_db equal, or are a subset, of
            mine
        """
        if not isinstance(other_db, AnnotationDbABC):
            msg = f"{type(other_db)} does not support features"
            raise TypeError(msg)
        if not isinstance(other_db, SqliteAnnotationDbMixin):
            return False

        mine = set(self.table_names)
        theirs = set(other_db.table_names)
        return mine <= theirs or mine > theirs if symmetric else mine >= theirs

    def update(
        self,
        annot_db: AnnotationDbABC,
        seqids: str | PySeq[str] | None = None,
        **kwargs: Any,
    ) -> None:
        """update records with those from an instance of the same type"""
        if not isinstance(annot_db, AnnotationDbABC):
            msg = f"{type(annot_db)} does not satisfy AnnotationDbABC"
            raise TypeError(msg)
        if not self.compatible(annot_db, symmetric=False):
            msg = f"{type(self)} cannot be updated from {type(annot_db)}"
            raise TypeError(msg)

        if not annot_db or not len(annot_db):
            return

        self._update_db_from_other_db(
            cast("SqliteAnnotationDbMixin", annot_db), seqids=seqids, **kwargs
        )

    def union(self, annot_db: AnnotationDbABC) -> Self:
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
        if annot_db and not isinstance(annot_db, AnnotationDbABC):
            msg = f"{type(annot_db)} does not satisfy AnnotationDbABC"
            raise TypeError(msg)
        if not annot_db:
            return copy.deepcopy(self)

        cls: type[Self]
        if self.compatible(annot_db, symmetric=False):
            cls = type(self)
        elif self.compatible(annot_db, symmetric=True):
            cls = type(annot_db)  # type: ignore[assignment]
        else:
            msg = f"cannot make a union between {type(self)} and {type(annot_db)}"
            raise TypeError(
                msg,
            )

        db = cls()
        db.update(self)  # type: ignore[arg-type]
        db.update(annot_db)
        return db

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        # make an empty db
        init_args = data.pop("init_args")
        db = cls(**init_args)
        db._update_db_from_rich_dict(data)
        return db

    def _update_db_from_rich_dict(
        self, data: dict[str, Any], seqids: str | Iterable[str] | None = None
    ) -> None:
        data.pop("type", None)
        data.pop("version")
        if isinstance(seqids, str):
            seqids = {seqids}
        elif seqids is not None:
            seqids = set(seqids) - {None}  # make sure None is not part of this!

        # TODO gah prevent duplication of existing records
        for table_name, records in data["tables"].items():
            for record in records:
                if seqids and record.get("seqid") not in seqids:
                    continue
                record["spans"] = numpy.array(record["spans"], dtype=int)
                # Translate string columns to IDs for normalized schema
                record = self._translate_record_to_ids(record)
                sql, vals = _add_record_sql(
                    table_name,
                    {k: v for k, v in record.items() if v is not None},
                )
                self._execute_sql(sql, vals)

    def _update_db_from_other_db(
        self,
        other_db: SqliteAnnotationDbMixin,
        seqids: str | PySeq[str] | None = None,
    ) -> None:
        if other_db == self:
            return

        if not isinstance(other_db, SqliteAnnotationDbMixin):
            msg = f"cannot update from {type(other_db)}"
            raise TypeError(msg)

        conditions: dict[str, Any] = {"seqid": seqids} if seqids else {}
        # Translate conditions for other_db's schema
        if other_db._lookup_cache is not None:
            conditions = other_db._translate_query_conditions(conditions)
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
            data = list(other_db._execute_sql(sql, vals))
            val_placeholder = ", ".join("?" * len(col_order[tname]))
            sql = f"INSERT INTO {tname} ({', '.join(col_order[tname])}) VALUES ({val_placeholder})"

            self.db.executemany(sql, data)
            self.db.commit()

    def to_json(self) -> str:
        return json.dumps(self.to_rich_dict())

    def write(self, path: PathType) -> None:
        """writes db as bytes to path

        Parameters
        ----------
        path
            Path to write the database. Must have suffix matching
            the class's _suffix attribute.

        Raises
        ------
        ValueError
            If the file suffix doesn't match the expected suffix for this class.
        """
        path = pathlib.Path(path).expanduser()
        expected_suffix = f".{self._suffix}"
        if not path.name.endswith(expected_suffix):
            msg = f"Expected file suffix {expected_suffix!r}, got {path.suffix!r}"
            raise ValueError(msg)
        backup = sqlite3.connect(path)
        with self.db:
            self.db.backup(backup)
        backup.close()

    @classmethod
    def from_file(cls, path: PathType) -> Self:
        """Load an annotation database from a file.

        Parameters
        ----------
        path
            Path to the saved database file. Must have suffix matching
            the class's _suffix attribute.

        Returns
        -------
        An instance of the annotation database class with the loaded data.

        Raises
        ------
        ValueError
            If the file suffix doesn't match the expected suffix for this class.
        OSError
            If the file doesn't exist.
        """
        path = pathlib.Path(path).expanduser()

        if not path.exists():
            msg = f"File {path} does not exist."
            raise OSError(msg)

        expected_suffix = f".{cls._suffix}"
        if not path.name.endswith(expected_suffix):
            msg = f"Expected file suffix {expected_suffix!r}, got {path.suffix!r}"
            raise ValueError(msg)

        conn = _make_db_connection(path, table_names=cls._table_names)
        return cls(source=path, db=conn)

    def subset(
        self,
        *,
        source: PathType = ":memory:",
        biotype: str | None = None,
        seqid: str | None = None,
        name: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        strand: str | None = None,
        attributes: str | None = None,
        allow_partial: bool = False,
    ) -> Self:
        """returns a new db instance with records matching the provided conditions"""
        # make sure python, not numpy, integers
        start = start if start is None else int(start)
        stop = stop if stop is None else int(stop)
        local_vars = locals()

        kwargs = {
            k: v
            for k, v in local_vars.items()
            if k not in {"self", "source"} and v is not None
        }
        if "strand" in kwargs:
            kwargs["strand"] = Strand.from_value(strand).value

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
            records = [
                tuple(record[c] for c in cols)
                for record in self._get_records_matching(
                    table_name=table_name,
                    **kwargs,
                )
            ]
            with result.db as cursor:
                cursor.executemany(sql, records)
                cursor.commit()

        return result

    def close(self) -> None:
        """closes the db"""
        if self._db is not None:
            self._db.close()
            self._db = None

    def _make_index(self, *, table_name: str, col_names: tuple[str, ...]) -> None:
        """index columns for faster search"""
        sql = f"CREATE INDEX IF NOT EXISTS %s on {table_name}(%s)"
        for col in col_names:
            self._execute_sql(sql % (col, col))

    def make_indexes(self) -> None:
        """adds db indexes for core attributes"""
        for table_name in self.table_names:
            self._make_index(
                table_name=table_name,
                col_names=(
                    "biotype_id",
                    "seqid_id",
                    "name",
                    "start",
                    "stop",
                    "parent_id",
                ),
            )


class BasicAnnotationDb(SqliteAnnotationDbMixin, AnnotationDbABC):
    """Provides a user table for annotations. This can be merged with
    either the Gff or Genbank versions.

    Notes
    -----
    This is the default db on Sequence, SequenceCollection and Alignment
    """

    _table_names: ClassVar = ("user",)
    _suffix = ANNDB_SUFFIX_BASIC

    def __init__(
        self,
        *,
        data: Iterable[dict[str, Any]] | None = None,
        db: SqliteAnnotationDbMixin | sqlite3.Connection | None = None,
        source: PathType = ":memory:",
    ) -> None:
        """
        Parameters
        ----------
        data
            data for entry into the database
        db
            a compatible annotation db instance. If db is the same class, it's
            db will be bound to self and directly modified.
        source
            location to store the db, defaults to in memory only. Must not be
            an existing file - use from_file() to load existing databases.
        """
        # Reject existing files - users should use from_file() instead
        # Only check when db is not provided (db=None means we're not binding to existing)
        if db is None and source != ":memory:":
            source_path = pathlib.Path(source).expanduser()
            if source_path.exists():
                msg = (
                    f"File {source_path} already exists. "
                    f"Use {self.__class__.__name__}.from_file() to load existing databases."
                )
                raise ValueError(msg)

        data = data or []
        # note that data is destroyed
        self._num_fakeids = 0
        self.source = source
        self._db = None
        # Initialise schema attributes
        self._schema_version = 1
        self._lookup_cache = None
        self._setup_db(db)

        self.add_records(data)

    def add_records(self, records: Iterable[dict[str, Any]], **kwargs: Any) -> None:
        table_name = self.table_names[0]  # only one name for this class
        for record in records:
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


def _merge_spans(
    old: NumpyIntArrayType, new: list[tuple[int, int]]
) -> NumpyIntArrayType:
    """returns sorted, merged old and new spans"""
    if len(new) == old.shape[0] and (old == new).all():
        return old

    new_array = numpy.array(sorted(new), dtype=old.dtype)
    return numpy.unique(numpy.concatenate([old, new_array]), axis=0)


class GffAnnotationDb(SqliteAnnotationDbMixin, AnnotationDbABC):
    """Support for annotations from gff files. Records that span multiple
    rows in the gff are merged into a single record."""

    _table_names: ClassVar = "gff", "user"
    # We are relying on an attribute name structured as _<table name>_schema
    _gff_schema: ClassVar = {
        "seqid_id": "INTEGER",
        "source_id": "INTEGER",
        "biotype_id": "INTEGER",  # type in GFF
        "start": "INTEGER",
        "stop": "INTEGER",
        "score": "TEXT",  # check defn
        "strand": "INTEGER",
        "phase": "TEXT",
        "attributes": "TEXT",
        "comments": "TEXT",
        "spans": "array",  # aggregation of coords across records
        "name": "TEXT",
        "parent_id": "TEXT",
    }
    _suffix = ANNDB_SUFFIX_GFF

    @extend_docstring_from(BasicAnnotationDb.__init__)
    def __init__(
        self,
        *,
        data: list[GffRecordABC] | None = None,
        db: SqliteAnnotationDbMixin | sqlite3.Connection | None = None,
        source: PathType = ":memory:",
    ) -> None:
        # Reject existing files - users should use from_file() instead
        # Only check when db is not provided (db=None means we're not binding to existing)
        if db is None and source != ":memory:":
            source_path = pathlib.Path(source).expanduser()
            if source_path.exists():
                msg = (
                    f"File {source_path} already exists. "
                    f"Use {self.__class__.__name__}.from_file() to load existing databases."
                )
                raise ValueError(msg)

        data = data or []
        # note that data is destroyed
        self._num_fakeids = 0
        # if a db instance is passed in, we use that for source
        self.source = getattr(db, "source", source)
        # ensure serialisable state reflects this
        self._serialisable["source"] = self.source
        self._db = None
        # Initialise schema attributes
        self._schema_version = 1
        self._lookup_cache = None
        self._setup_db(db)
        merged_data, self._num_fakeids = merged_gff_records(data, self._num_fakeids)
        self.add_records(merged_data)

    # This definition is incompatible with the base class.
    def add_records(self, reduced: dict[str, GffRecordABC], **kwargs: Any) -> None:  # type: ignore[override]
        col_order = [
            r["name"] for r in self.db.execute("PRAGMA table_info(gff)").fetchall()
        ]
        val_placeholder = ", ".join("?" * len(col_order))
        sql = f"INSERT INTO gff ({', '.join(col_order)}) VALUES ({val_placeholder})"

        # Can we really trust text case of "ID" and "Parent" to be consistent
        # across sources of gff?? I doubt it, so more robust regex likely
        # required
        rows: list[tuple[Any, ...]] = []
        hierarchy_entries: list[tuple[str, str | None]] = []

        for record in reduced.values():
            # our Feature code assumes start always < stop,
            # we record direction using Strand
            spans = numpy.array(sorted(record["spans"]), dtype=int)  # sorts the rows
            spans.sort(axis=1)
            start = int(spans.min())
            stop = int(spans.max())

            # Build a row dict with translated IDs (since GffRecord objects are read-only)
            row_data: dict[str, Any] = {
                "start": start,
                "stop": stop,
                "spans": spans,
                "strand": Strand.from_value(record.get("strand")).value,
                "name": record.get("name"),
                "parent_id": record.get("parent_id"),
                "score": record.get("score"),
                "phase": record.get("phase"),
                "attributes": record.get("attributes"),
                "comments": record.get("comments"),
            }

            # Translate seqid/biotype/source to IDs
            if self._lookup_cache is not None:
                row_data["seqid_id"] = (
                    self._lookup_cache.get_seqid_id(record.get("seqid"))
                    if record.get("seqid")
                    else 0
                )
                row_data["biotype_id"] = (
                    self._lookup_cache.get_biotype_id(record.get("biotype"))
                    if record.get("biotype")
                    else 0
                )
                row_data["source_id"] = (
                    self._lookup_cache.get_source_id(record.get("source"))
                    if record.get("source")
                    else 0
                )
            else:
                row_data["seqid_id"] = 0
                row_data["biotype_id"] = 0
                row_data["source_id"] = 0

            rows.append(tuple(row_data.get(c) for c in col_order))

            # Collect data for hierarchy
            name = record.get("name") or ""
            if parent_id := record.get("parent_id"):
                hierarchy_entries.append((name, parent_id))

        self.db.executemany(sql, rows)
        self.db.commit()

        # Populate hierarchy table
        self._populate_hierarchy_batch("gff", hierarchy_entries)

        del reduced

    def _populate_hierarchy_batch(
        self,
        table_name: str,
        entries: list[tuple[str, str | None]],
    ) -> None:
        """Populate hierarchy table for a batch of records.

        Parameters
        ----------
        table_name
            Name of feature table (gff, gb, user)
        entries
            List of (child_name, parent_id_string) tuples
        """
        if not self._can_use_hierarchy():
            return

        # Collect child names that have parent_id
        child_names: set[str] = {
            child_name for child_name, parent_id in entries if parent_id
        }
        if not child_names:
            return

        # Batch lookup child names to rowids
        name_to_rowid = self._get_feature_rowids_batch(table_name, child_names)

        # Convert to (child_rowid, parent_id) format for shared helper
        rowid_entries = [
            (name_to_rowid[child_name], parent_id)
            for child_name, parent_id in entries
            if parent_id and child_name in name_to_rowid
        ]

        _populate_hierarchy_table(self.db, table_name, rowid_entries)

    def update_record_spans(self, *, name: str, spans: list[tuple[int, int]]) -> None:
        """updates spans attribute of a gff table record if present

        Notes
        -----
        Has no effect if name is not present.
        """
        if not len(spans):
            return

        result = self._execute_sql(
            cmnd="SELECT spans from gff WHERE name = ?",
            values=(name,),
        ).fetchone()

        if result is None:
            return

        old_spans = _merge_spans(result["spans"], spans)
        self._execute_sql(
            cmnd="UPDATE gff SET spans = ? WHERE name = ?",
            values=(old_spans, name),
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


class GenbankAnnotationDb(SqliteAnnotationDbMixin, AnnotationDbABC):
    """Support for annotations from Genbank files.

    Notes
    -----
    Extended attributes are stored as json in the gb, attributes column.
    """

    _table_names: ClassVar = "gb", "user"
    # We are relying on an attribute name structured as _<table name>_schema
    _gb_schema: ClassVar = {
        "seqid_id": "INTEGER",
        "source_id": "INTEGER",
        "biotype_id": "INTEGER",  # type in GFF
        "start": "INTEGER",
        "stop": "INTEGER",
        "strand": "INTEGER",
        "comments": "TEXT",
        "spans": "array",  # aggregation of coords across records
        "name": "TEXT",
        "symbol": "TEXT",
        "parent_id": "TEXT",
        "attributes": "json",
    }
    _suffix = ANNDB_SUFFIX_GENBANK

    @extend_docstring_from(BasicAnnotationDb.__init__)
    def __init__(
        self,
        *,
        data: Iterable[dict[str, Any]] | None = None,
        seqid: str | None = None,
        db: SqliteAnnotationDbMixin | sqlite3.Connection | None = None,
        source: PathType = ":memory:",
        namer: Callable[[dict[str, Any]], list[str] | None] | None = None,
    ) -> None:
        """
        seqid
            name of the sequence data is associated with
        namer
            callable that takes a record dict and returns a
            [name]
        """
        # Reject existing files - users should use from_file() instead
        # Only check when db is not provided (db=None means we're not binding to existing)
        if db is None and source != ":memory:":
            source_path = pathlib.Path(source).expanduser()
            if source_path.exists():
                msg = (
                    f"File {source_path} already exists. "
                    f"Use {self.__class__.__name__}.from_file() to load existing databases."
                )
                raise ValueError(msg)

        data = data or []
        # note that data is destroyed
        self._num_fakeids = 0
        self.source = source
        self._db = None
        # Initialise schema attributes
        self._schema_version = 1
        self._lookup_cache = None
        self._setup_db(db)
        if db:
            # guessing there's multiple seqid's
            self._serialisable["seqid"] = "<multiple seqids>"

        self._namer = namer if callable(namer) else self._default_namer
        self.add_records(data, seqid)

    def add_records(
        self, records: Iterable[dict[str, Any]], seqid: str | None = None, **kwargs: Any
    ) -> None:
        col_order = [
            r["name"] for r in self.db.execute("PRAGMA table_info(gb)").fetchall()
        ]

        val_placeholder = ", ".join("?" * len(col_order))
        sql = f"INSERT INTO gb ({', '.join(col_order)}) VALUES ({val_placeholder})"

        # need to capture genetic code from genbank
        rows: list[tuple[Any, ...]] = []
        exclude = {"translation", "location", "type"}  # location is grabbed directly
        for record in records:
            # our Feature code assumes start always < stop,
            store: dict[str, Any] = {}
            # we create the location data directly
            if location := record.get("location", None):
                store["spans"] = numpy.array(location.get_coordinates(), dtype=int)
                if strand := location.strand:
                    store["strand"] = Strand.from_value(strand).value
                store["start"] = int(store["spans"].min())
                store["stop"] = int(store["spans"].max())

            attrs_keys = record.keys() - exclude
            biotype = record["type"]
            name = self._namer(record)
            if name is None:
                name = self._make_fake_id(record)
            store["attributes"] = {k: record[k] for k in attrs_keys}
            store["name"] = ",".join(name)

            # Translate seqid/biotype to IDs
            if self._lookup_cache is not None:
                store["seqid_id"] = (
                    self._lookup_cache.get_seqid_id(seqid) if seqid else 0
                )
                store["biotype_id"] = (
                    self._lookup_cache.get_biotype_id(biotype) if biotype else 0
                )
            else:
                store["seqid_id"] = 0
                store["biotype_id"] = 0

            rows.append(tuple(store.get(c) for c in col_order))

        self.db.executemany(sql, rows)
        self.db.commit()

        del records

    def _default_namer(self, record: dict[str, Any]) -> list[str] | None:
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

    def _make_fake_id(self, record: dict[str, Any]) -> list[str]:
        name = [f"{record.get('type', 'fakeid')}-{self._num_fakeids}"]
        self._num_fakeids += 1
        return name

    def get_feature_children(
        self,
        name: str,
        biotype: str | None = None,
        exclude_biotype: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        **kwargs: Any,
    ) -> Iterator[FeatureDataType]:
        """yields children of name"""
        # children must lie within the parent coordinates
        if start is None or stop is None:
            name = self.__class__.__name__
            msg = f"coordinates required to query children for {name!r}"
            raise ValueError(msg)

        for table_name in self.table_names:
            columns: tuple[str, ...] = (
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
                cstart: int = feat.pop("start")  # type: ignore[typeddict-item]
                cstop: int = feat.pop("stop")  # type: ignore[typeddict-item]
                if not (start <= cstart < stop and start < cstop <= stop):
                    continue

                # remove invalid field for the FeatureDataType
                feat.pop("parent_id")  # type: ignore[typeddict-item]
                yield feat

    def get_feature_parent(
        self,
        name: str,
        exclude_biotype: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        **kwargs: Any,
    ) -> Iterator[FeatureDataType]:
        """yields parents of name"""
        if start is None or stop is None:
            name = self.__class__.__name__
            msg = f"coordinates required to query parent for {name!r}"
            raise ValueError(msg)

        # parent must match or lie outside the coordinates
        for table_name in self.table_names:
            columns: tuple[str, ...] = (
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

                cstart: int = feat.pop("start")  # type: ignore[typeddict-item]
                cstop: int = feat.pop("stop")  # type: ignore[typeddict-item]
                if cstart > start or stop > cstop:
                    continue

                # remove invalid field for the FeatureDataType
                feat.pop("parent_id")  # type: ignore[typeddict-item]
                yield feat


@register_deserialiser(get_object_provenance(BasicAnnotationDb))
def deserialise_basic_db(data: dict[str, Any]) -> BasicAnnotationDb:
    return BasicAnnotationDb.from_dict(data)


@register_deserialiser(get_object_provenance(GffAnnotationDb))
def deserialise_gff_db(data: dict[str, Any]) -> GffAnnotationDb:
    return GffAnnotationDb.from_dict(data)


@register_deserialiser(get_object_provenance(GenbankAnnotationDb))
def deserialise_gb_db(data: dict[str, Any]) -> GenbankAnnotationDb:
    return GenbankAnnotationDb.from_dict(data)


@register_deserialiser("annotation_to_annotation_db")
def convert_annotation_to_annotation_db(data: dict[str, Any]) -> AnnotationDbABC:
    db = BasicAnnotationDb()

    minus = Strand.MINUS
    plus = Strand.PLUS

    seqid = data.pop("name", data.pop("seqid", None))
    anns = data.pop("data")
    for ann in anns:
        ann = ann.pop("annotation_construction")  # noqa: PLW2901
        m = deserialise_map_spans(ann.pop("map"))
        spans = m.get_coordinates()
        strand = minus if any(s.reverse for s in m.iter_non_lost_spans()) else plus
        biotype = ann.pop("type")
        name = ann.pop("name")
        db.add_feature(
            seqid=seqid,
            biotype=biotype,
            name=name,
            spans=spans,
            strand=strand,
        )

    return db


def _db_from_genbank(
    path: PathType,
    db: SqliteAnnotationDbMixin | None,
    write_path: PathType,
    **kwargs: Any,
) -> SqliteAnnotationDbMixin:
    from cogent3.parse.genbank import minimal_parser

    kwargs.pop("ui", None)

    rec = next(iter(minimal_parser(path)))
    data = cast("Iterable[dict[str, Any]] | None", rec.pop("features", None))
    locus = cast("str", rec["locus"])
    db = GenbankAnnotationDb(
        source=write_path,
        data=data,
        seqid=locus,
        db=db,
    )

    db.make_indexes()
    return db


def _leave_attributes(*attrs: str) -> str:
    return attrs[0]


def _db_from_gff(
    path: PathType,
    seqids: str | Iterable[str] | None,
    db: SqliteAnnotationDbMixin | None,
    write_path: PathType,
    num_lines: int | None,
    **kwargs: Any,
) -> SqliteAnnotationDbMixin:
    from cogent3.parse.gff import gff_parser, is_gff3

    kwargs.pop("ui", None)

    num_fake_ids = 0
    seen_ids: set[str] = set()
    gff3 = is_gff3(path)
    db = GffAnnotationDb(source=write_path, db=db)
    for block in iter_line_blocks(path, num_lines=num_lines):
        data = list(
            gff_parser(
                block,
                seqids=seqids,
                attribute_parser=_leave_attributes,
                gff3=gff3,
            ),
        )
        merged_data, num_fake_ids = merged_gff_records(data, num_fake_ids)
        if already_seen := seen_ids & merged_data.keys():
            for name in already_seen:
                db.update_record_spans(
                    name=name,
                    spans=cast("list[tuple[int, int]]", merged_data[name].spans),
                )

        seen_ids |= merged_data.keys()
        db.add_records(merged_data)

    assert db is not None
    db.make_indexes()
    return db


class SqliteAnnotationDbLoader(AnnotationDbLoaderBase):
    """Loader for annotation files using SQLite storage backend."""

    @property
    def name(self) -> str:
        return "c3anndb"

    @property
    def supported_suffixes(self) -> set[str]:
        return {
            ".gff",
            ".gff3",
            ".gb",
            ".gbk",
            ".gbff",
            ".json",
            f".{ANNDB_SUFFIX_GENBANK}",
            f".{ANNDB_SUFFIX_GFF}",
            f".{ANNDB_SUFFIX_BASIC}",
        }

    def load(
        self,
        path: PathType,
        seqids: str | Iterable[str] | None = None,
        db: AnnotationDbABC | None = None,
        write_path: PathType = ":memory:",
        format_name: str | None = None,
        **kwargs: Any,
    ) -> AnnotationDbABC:
        """Load annotations using SQLite backend.

        Supports GFF, GenBank, and JSON formats.

        Parameters
        ----------
        path
            Path to annotation file (may contain glob patterns)
        seqids
            Filter to specific sequence IDs
        db
            Existing database instance to add records to
        write_path
            Path for database file (":memory:" for in-memory)
        format_name
            Explicit format override ('gff', 'genbank', 'json').
            If None, auto-detect from file suffix.
        **kwargs
            Additional arguments (e.g., ui for progress display,
            lines_per_block for GFF chunked loading)
        """
        # was this an already saved db
        path = pathlib.Path(path)
        if path.name.endswith(f".{ANNDB_SUFFIX_GENBANK}"):
            return GenbankAnnotationDb.from_file(path)
        if path.name.endswith(f".{ANNDB_SUFFIX_GFF}"):
            return GffAnnotationDb.from_file(path)
        if path.name.endswith(f".{ANNDB_SUFFIX_BASIC}"):
            return BasicAnnotationDb.from_file(path)

        # Determine format from format_name parameter or file suffix
        if format_name:
            fmt = format_name.lower()
        else:
            suffix, _ = get_format_suffixes(path)
            if suffix is None:
                msg = f"Unknown annotation format for {path!r}"
                raise ValueError(msg)
            suffix = (suffix if suffix.startswith(".") else f".{suffix}").lower()
            if suffix == ".json":
                fmt = "json"
            elif suffix in {".gb", ".gbk", ".gbff"}:
                fmt = "genbank"
            elif suffix in {".gff", ".gff3"}:
                fmt = "gff"
            else:
                msg = f"Unknown annotation format for suffix {suffix!r}"
                raise ValueError(msg)

        # Check for glob pattern
        p = pathlib.Path(path)
        if fmt == "json":
            path_str = str(p)
            has_glob = "*" in path_str or "?" in path_str or "[" in path_str

            if has_glob:
                # JSON format doesn't support glob patterns
                msg = "Glob patterns not supported for JSON format"
                raise NotImplementedError(msg)
            return deserialise_object(path)

        # Select parsing function and prepare kwargs before loop
        parse_func: Callable[..., SqliteAnnotationDbMixin]
        if fmt == "genbank":
            parse_func = _db_from_genbank
            func_kwargs = {"write_path": write_path, **kwargs}
        elif fmt == "gff":
            parse_func = _db_from_gff
            func_kwargs = {
                "seqids": seqids,
                "write_path": write_path,
                "num_lines": 500_000,
                "lines_per_block": 500_000,
                **kwargs,
            }
        else:
            msg = f"Unknown format {format_name!r}"
            raise ValueError(msg)

        # Expand glob pattern
        paths = sorted(p.parent.glob(p.name))
        if not paths:
            msg = f"No files found matching pattern {path!r}"
            raise OSError(msg)

        # Get progress display
        show_progress = kwargs.pop("show_progress", False)
        kwargs.pop("ui", None)  # remove legacy ui kwarg
        progress = (
            get_progress(show_progress=True, **show_progress)
            if isinstance(show_progress, dict)
            else get_progress(show_progress)
        )
        series = progress(paths, msg="Loading annotations")

        # Process each file, passing db from one to the next
        db_result: SqliteAnnotationDbMixin | None = cast(
            "SqliteAnnotationDbMixin | None", db
        )
        for file_path in series:
            db_result = parse_func(path=file_path, db=db_result, **func_kwargs)

        if db_result is None:
            msg = f"Failed to load annotations into database {str(path)!r}"
            raise RuntimeError(msg)
        return db_result


def load_annotations(
    *,
    path: PathType,
    seqids: str | Iterable[str] | None = None,
    db: AnnotationDbABC | None = None,
    write_path: PathType = ":memory:",
    lines_per_block: int | None = 500_000,
    show_progress: bool | Progress | dict[str, Any] = False,  # type: ignore[type-arg]
    format_name: str | None = None,
    storage_backend: str | None = None,
) -> AnnotationDbABC:
    """loads annotations from flatfile into a db

    Parameters
    ----------
    path
        path to a plain text file containing annotations, or a json file
        of a serialised cogent3 annotation db object
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
    format_name
        explicitly specify annotation format ('gff', 'genbank', 'json').
        If not provided, format is auto-detected from file suffix.
    storage_backend
        storage backend to use for the annotation database (e.g., 'c3anndb').
        If not provided, selects first compatible third-party plugin, falling back
        to cogent3 built-in SQLite loaders.

    Notes
    -----
    We DO NOT check if a provided db already contains records from a flatfile.
    """
    from cogent3._plugin import get_annotation_loader_plugin

    if seqids is not None:
        seqids = {seqids} if isinstance(seqids, str) else set(seqids)

    path = pathlib.Path(path).expanduser()

    # Determine file suffix for format detection
    if format_name:
        # Explicit format specified - map to suffix
        suffix = format_name if format_name.startswith(".") else f".{format_name}"
    else:
        # Auto-detect from file suffix
        suffix, _ = get_format_suffixes(path)
        if suffix is None:
            msg = f"Cannot auto-detect format for {path!r}, use format_name parameter"
            raise ValueError(msg)
        suffix = (suffix if suffix.startswith(".") else f".{suffix}").lower()

    # Get the loader plugin
    loader = get_annotation_loader_plugin(
        storage_backend=storage_backend,
        file_suffix=suffix,
    )

    return loader.load(
        path=path,
        seqids=seqids,
        db=db,
        write_path=write_path,
        format_name=format_name,
        lines_per_block=lines_per_block,
        show_progress=show_progress,
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
        return array_to_sqlite(sqlite_to_array(data))


def upgrade_annotation_db(
    source_path: PathType,
    backup: bool = True,
) -> None:
    """Upgrade an existing annotation database to the schema.

    This function migrates old-format annotation databases to use the new
    schema, which includes:
    - Lookup tables for seqid, biotype, source normalization
    - Feature hierarchy table for fast parent/child lookups

    Parameters
    ----------
    source_path
        Path to the existing database file
    backup
        If True, creates backup at source_path.bak before migration

    Notes
    -----
    - Detects schema version and skips if already upgraded
    - Preserves all existing data
    - Creates lookup tables and hierarchy table
    - Running twice is safe

    Raises
    ------
    OSError
        If source_path does not exist
    FileExistsError
        If backup=True and backup file already exists

    Examples
    --------
    >>> from cogent3.core.annotation_db import upgrade_annotation_db
    >>> upgrade_annotation_db("my_annotations.gffdb")  # doctest: +SKIP
    """
    source_path = pathlib.Path(source_path).expanduser()

    if not source_path.exists():
        msg = f"File {source_path} does not exist."
        raise OSError(msg)

    # Open the database to check schema version
    conn = _make_db_connection(source_path)

    current_version = _get_schema_version(conn)
    if current_version >= SCHEMA_VERSION:
        conn.close()
        return  # Already at current version

    # Create backup if requested
    if backup:
        backup_path = source_path.parent / f"{source_path.name}.bak"
        if backup_path.exists():
            conn.close()
            msg = (
                f"Backup file already exists at {backup_path}. "
                "Remove or rename it before upgrading."
            )
            raise FileExistsError(msg)

        backup_conn = sqlite3.connect(backup_path)
        conn.backup(backup_conn)
        backup_conn.close()

    # Detect which table type this is (gff, gb, or just user)
    table_names: list[str] = []
    for table in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'"
    ).fetchall():
        name = table["name"]
        if name in ("gff", "gb", "user"):
            table_names.append(name)

    # Create schema tables
    for sql in _LOOKUP_TABLES_SQL.values():
        conn.execute(sql)
    conn.execute(_FEATURE_SPANS_SQL)
    conn.execute(_FEATURE_HIERARCHY_SQL)
    conn.commit()

    # Build lookup tables from existing data
    lookup_cache = LookupTableCache(conn)

    # Migrate old TEXT columns to new *_id INTEGER columns
    _migrate_to_normalized_schema(conn, table_names, lookup_cache)

    for table_name in table_names:
        # Get all records from this table
        cursor = conn.execute(f"SELECT rowid, * FROM {table_name}")
        records = cursor.fetchall()

        # Build hierarchy entries as (child_rowid, parent_id_string) tuples
        hierarchy_entries: list[tuple[int, str]] = []

        for row in records:
            row_dict = dict(zip(row.keys(), row, strict=False))
            rowid = row_dict.get("rowid")
            parent_id = row_dict.get("parent_id")

            # Collect hierarchy entries
            if parent_id and rowid is not None:
                hierarchy_entries.append((rowid, parent_id))

        # Populate hierarchy table using shared helper
        _populate_hierarchy_table(conn, table_name, hierarchy_entries)

    # Set schema version
    _set_schema_version(conn, SCHEMA_VERSION)

    conn.close()


DEFAULT_ANNOTATION_DB = BasicAnnotationDb

# This is a design note about the annotation_db property of SequenceCollection, Alignment,
#   Aligned and Sequence classes.
#   We use an empty list as the default value for two reasons. First, many applications
#   of these objects do not employ annotations and thus to reduce memory overhead we
#   create a db instance lazily. Second, members of a collection can be accessed directly
#   and have their annotation db accessed. Construction of those members is controlled
#   by the container class. To avoid tightly coupling the container and its enclosed
#   sequences, we provide a mutable data object (an empty list) as the starting value.
#   This list is provided to the constructed elements as their starting value. So an
#   Aligned instance, for example, that is created and then used to add a feature, has
#   created the db instance shared by all members of the collection and the collection
#   itself. The list can only have 0 or 1 element.


class AnnotatableMixin:
    """class handling an annotation database for a collection, aligned or sequence"""

    def __init__(self) -> None:  # pragma: no cover
        self._annotation_db: list[AnnotationDbABC] = []

    def _init_annot_db_value(
        self, value: AnnotationDbABC | list[AnnotationDbABC] | None
    ) -> list[AnnotationDbABC]:
        """returns the value for assignment to self._annotation_db given value

        Parameters
        ----------
        value
            the provided value to interpret for assignment to self._annotation_db.
            Ccan be None, an annotation db instance or a list containing an annotation
            db instance.
        """
        if isinstance(value, list):
            return value

        if value is None:
            return []

        # if annotation_db is not a list, assume it is a single
        # annotation database instance
        return [value]

    @property
    def annotation_db(self) -> AnnotationDbABC:
        """the annotation database for the collection"""
        if not self._annotation_db:
            # if no annotation db is set, use the default
            self._annotation_db.append(DEFAULT_ANNOTATION_DB())
        return self._annotation_db[0]

    @annotation_db.setter
    def annotation_db(self, value: AnnotationDbABC | None) -> None:
        # Without knowing the contents of the db we cannot
        # establish whether self.moltype is compatible, so
        # we rely on the user to get that correct
        # one approach to support validation might be to add
        # to the AnnotationDbABC protocol a is_nucleic flag,
        # for both DNA and RNA. But if a user trys get_slice()
        # on a '-' strand feature, they will get a TypeError.
        # I think that's enough.
        self.replace_annotation_db(value, check=False)

    def replace_annotation_db(
        self,
        value: AnnotationDbABC | list[AnnotationDbABC] | None,
        check: bool = True,
    ) -> None:
        """public interface to assigning the annotation_db

        Parameters
        ----------
        value
            the annotation db instance
        check
            whether to check value supports the feature interface

        Notes
        -----
        The check can be very expensive, so if you're confident set it to False
        """
        if not value and isinstance(value, list):
            self._annotation_db = value
            return

        if isinstance(value, list) and len(value) > 0:
            value = value[0]

        if value in self._annotation_db:
            return

        if value is None:
            self._annotation_db = self._init_annot_db_value(None)
            return

        if check and value and not isinstance(value, AnnotationDbABC):
            msg = f"{type(value)} does not satisfy AnnotationDbABC"
            raise TypeError(msg)

        value = cast("AnnotationDbABC", value)
        if not self._annotation_db:
            # if no annotation db is set, use the default
            self._annotation_db.append(value)
        elif self._annotation_db[0] is not value:
            # we close the current annotation db and replace it
            self._annotation_db[0].close()
            self._annotation_db[0] = value

        return
