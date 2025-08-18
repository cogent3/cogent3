from __future__ import annotations

import collections
import functools
import io
import pathlib
import re
from abc import ABC
from collections.abc import Callable, Iterable
from collections.abc import Sequence as PySeq
from typing import IO, Any, cast

from cogent3.util.io import PathType, open_


@functools.singledispatch
def is_gff3(f: PathType) -> bool:
    """True if gff-version is 3"""
    msg = f"unsupported type {type(f)}"
    raise TypeError(msg)


@is_gff3.register
def _(f: pathlib.Path) -> bool:
    with open_(f) as f:
        return is_gff3(f)


@is_gff3.register
def _(f: str) -> bool:
    with open_(f) as f:
        return is_gff3(f)


@is_gff3.register
def _(f: io.IOBase) -> bool:
    gff3_header = "gff-version 3"
    is_gff3 = gff3_header in f.readline()
    f.seek(0)
    return is_gff3


@is_gff3.register
def _(f: list) -> bool:
    return is_gff3(tuple(f))


@is_gff3.register
def _(f: tuple) -> bool:
    gff3_header = "gff-version 3"
    return f and gff3_header in f[0]


_gff_item_map = {
    "type": "biotype",
    "end": "stop",
    "attributes": "attrs",
}


class GffRecordABC(ABC):
    """abstract base class for gff records"""

    # attributes that need to be on the instance
    __slots__ = (
        "attrs",
        "biotype",
        "comments",
        "name",
        "parent_id",
        "phase",
        "score",
        "seqid",
        "source",
        "spans",
        "start",
        "stop",
        "strand",
    )

    def __init__(
        self,
    ) -> None:  # pragma: no cover
        self.seqid: str | None
        self.source: str | None
        self.biotype: str | None
        self.name: str | None
        self.parent_id: str | None
        self.start: int | None
        self.stop: int | None
        self.spans: list[tuple[int, int]] | None
        self.score: str | None
        self.strand: str | None
        self.phase: str | None
        self.comments: str | None
        self.attrs: str | dict[str, str] | Any | None

    def __repr__(self) -> str:
        name = self.__class__.__name__
        return (
            f"{name}(seqid={self.seqid!r}, name={self.name!r}, "
            f"biotype={self.biotype!r}, start={self.start}, stop={self.stop}, "
            f"strand={self.strand!r}, spans={self.spans})"
        )

    def __getitem__(self, item: str) -> Any:
        item = item.lower()
        return getattr(self, _gff_item_map.get(item, item))

    def __setitem__(self, key: str, value: Any) -> None:
        key = key.lower()
        setattr(self, _gff_item_map.get(key, key), value)

    def update(self, values: dict[str, Any]) -> None:
        for key, value in values.items():
            self[key] = value

    def get(self, item: str) -> Any:
        item = item.lower()
        return getattr(self, _gff_item_map.get(item, item), None)

    def to_dict(self) -> dict[str, Any]:
        return {k: getattr(self, k) for k in self.__slots__}


class GffRecord(GffRecordABC):
    """Container for single GFF record. Elements can be indexed as a dict
    or directly as attributes.

    Notes
    -----
    Some fields defined in the GFF spec are aliased, or truncated.
    """

    # should be a dataclass, but in py 3.9 they don't support slots
    # so until then ...
    __slots__ = (
        "attrs",
        "biotype",
        "comments",
        "name",
        "parent_id",
        "phase",
        "score",
        "seqid",
        "source",
        "spans",
        "start",
        "stop",
        "strand",
    )

    def __init__(
        self,
        seqid: str | None = None,
        source: str | None = None,
        biotype: str | None = None,
        name: str | None = None,
        parent_id: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        spans: list[tuple[int, int]] | None = None,
        score: str | None = None,
        strand: str | None = None,
        phase: str | None = None,
        comments: str | None = None,
        attrs: str | dict[str, str] | Any | None = None,
    ) -> None:
        self.seqid = seqid
        self.source = source
        self.biotype = biotype
        self.name = name
        self.parent_id = parent_id
        self.start = start
        self.stop = stop
        self.spans = spans
        self.score = score
        self.strand = strand
        self.phase = phase
        self.comments = comments
        self.attrs = attrs


@functools.singledispatch
def gff_parser(
    f: PathType | IO[str] | str | PySeq[str] | None,
    attribute_parser: Callable[[str], Any] | None = None,
    seqids: str | Iterable[str] | None = None,
    gff3: bool | None = None,
    make_record: type[GffRecord] | None = None,
) -> Iterable[GffRecordABC]:
    """parses a gff file

    Parameters
    -----------
    f
        accepts string path or pathlib.Path or file-like object (e.g. StringIO)
    attribute_parser
        callback function for custom handling of attributes field. Defaults to
        a function that converts Attributes field into a dict based on the gff
        version.
    seqids
        only gff records matching these sequence IDs are returned
    gff3
        True if the format is gff3
    make_record
        callable that constructs a GffRecordABC compatible record. Defaults to
        GffRecord. Over-rides attribute_parser if provided.

    Returns
    -------
    iterable of GffRecordABCs
    """
    gff3 = gff3 or is_gff3(f)
    yield from _gff_parser(
        f,
        attribute_parser=attribute_parser,
        seqids=seqids,
        gff3=gff3,
        make_record=make_record,
    )


@gff_parser.register
def _(
    f: pathlib.Path,
    attribute_parser: Callable[[str], Any] | None = None,
    seqids: str | Iterable[str] | None = None,
    gff3: bool | None = None,
    make_record: type[GffRecord] | None = None,
) -> Iterable[GffRecordABC]:
    with open_(f) as infile:
        yield from gff_parser(
            infile,
            attribute_parser=attribute_parser,
            seqids=seqids,
            gff3=gff3,
            make_record=make_record,
        )


@gff_parser.register
def _(
    f: str,
    attribute_parser: Callable[[str], Any] | None = None,
    seqids: str | Iterable[str] | None = None,
    gff3: bool | None = None,
    make_record: type[GffRecord] | None = None,
) -> Iterable[GffRecordABC]:
    with open_(f) as infile:
        yield from gff_parser(
            infile,
            attribute_parser=attribute_parser,
            seqids=seqids,
            gff3=gff3,
            make_record=make_record,
        )


@gff_parser.register
def _(
    f: io.IOBase,
    attribute_parser: Callable[[str], Any] | None = None,
    seqids: str | Iterable[str] | None = None,
    gff3: bool | None = None,
    make_record: type[GffRecord] | None = None,
) -> Iterable[GffRecordABC]:
    yield from _gff_parser(
        f,
        attribute_parser=attribute_parser,
        seqids=seqids,
        gff3=gff3,
        make_record=make_record,
    )


def _gff_parser(
    f: Iterable[str],
    attribute_parser: Callable[[str], Any] | None = None,
    seqids: str | Iterable[str] | None = None,
    gff3: bool | None = True,
    make_record: type[GffRecord] | None = None,
) -> Iterable[GffRecordABC]:
    """parses a gff file"""
    seqids = seqids or set()
    seqids = {seqids} if isinstance(seqids, str) else set(seqids)

    if attribute_parser is None and make_record is None:
        attribute_parser = parse_attributes_gff3 if gff3 else parse_attributes_gff2
    elif callable(make_record):
        attribute_parser = None

    if make_record is None:
        make_record = GffRecord

    for line in f:
        # comments and blank lines
        if "#" in line:
            (line, comments) = line.split("#", 1)
        else:
            comments = None
        line = line.strip()
        if not line:
            continue

        cols = [c.strip() for c in line.split("\t")]
        # the final column (attributes) may be empty
        if len(cols) == 8:
            cols.append("")
        assert len(cols) == 9, repr(cols)
        seqid, source, type_, start_, end_, score, strand, phase, attributes = cols

        if seqids and seqid not in seqids:
            continue

        # adjust for 0-based indexing
        start, end = int(start_) - 1, int(end_)
        # start is always meant to be less than end in GFF
        # features that extend beyond sequence have negative indices
        if start < 0 or end < 0:
            start, end = abs(start), abs(end)
        if start > end:
            start, end = end, start

        attributes = attribute_parser(attributes) if attribute_parser else attributes

        yield make_record(
            seqid=seqid,
            source=source,
            biotype=type_,
            start=start,
            stop=end,
            score=score,
            strand=strand,
            phase=phase,
            attrs=attributes,
            comments=comments,
        )


def parse_attributes_gff2(attributes: str) -> dict[str, str]:
    """Returns a dict with name and info keys"""
    name = attributes[attributes.find('"') + 1 :]
    if '"' in name:
        name = name[: name.find('"')]
    return {"ID": name, "Info": attributes}


def parse_attributes_gff3(attributes: str) -> dict[str, str | list[str]]:
    """Returns a dictionary containing all the attributes"""
    attributes = attributes.strip(";")
    attributes_split = attributes.split(";")
    attributes_dict: dict[str, str | list[str]] = cast(
        "dict[str, str | list[str]]",
        dict(t.split("=") for t in attributes_split) if attributes_split[0] else {},
    )

    if "Parent" in attributes_dict:
        # There may be multiple parents
        if "," in attributes_dict["Parent"]:
            attributes_dict["Parent"] = cast("str", attributes_dict["Parent"]).split(
                ","
            )
        else:
            attributes_dict["Parent"] = [cast("str", attributes_dict["Parent"])]
    attributes_dict["ID"] = attributes_dict.get("ID", "")
    return attributes_dict


def merged_gff_records(
    records: list[GffRecordABC],
    num_fake_ids: int,
) -> tuple[dict[str, GffRecordABC], int]:
    """merges GFF records that have the same ID

    Parameters
    ----------
    records
        list of GFF data rows as dict
    num_fake_ids
        starting number for new fakeids

    Returns
    -------
    a dictionary with IDs as keys and values being all attributes

    Notes
    -----
    Merging means the coordinates are combined into a single "spans"
    entry. Only records which have an ID field in the attributes get merged.
    """
    field_template = r"(?<={}=)[^;\s]+"
    name = re.compile(field_template.format("ID"))
    parent_id = re.compile(field_template.format("Parent"))

    reduced: dict[str, GffRecordABC] = collections.OrderedDict()
    # collapse records with ID's occurring in multiple rows into one
    # row, converting their coordinates
    # extract the name from ID and add this into the table
    # I am not convinced we can rely on gff files to be ordered,
    # if we could, we could do this as one pass over the data
    # all keys need to be lower case
    for record in records:
        attrs = cast("str", record.attrs or "")
        if match := name.search(attrs):
            record_id = match.group()
        else:
            record_id = f"unknown-{num_fake_ids}"
            num_fake_ids += 1

        record.name = record_id
        if pid := parent_id.search(attrs):
            record.parent_id = pid.group()

        if record_id not in reduced:
            reduced[record_id] = record
            reduced[record_id].spans = []

        record_spans_list = cast("list[tuple[int, int]]", reduced[record_id].spans)
        record_spans_list.append((record["start"], record["stop"]))

    del records

    return reduced, num_fake_ids
