from __future__ import annotations

import collections
import functools
import io
import re
import typing

from cogent3.util.io import PathType, open_


OptionalCallable = typing.Optional[typing.Callable]
OptionalStrContainer = typing.Optional[typing.Union[str, typing.Sequence[str]]]
OptionalIntList = typing.Optional[list[list[int]]]
OptionalStr = typing.Optional[str]
OptionalInt = typing.Optional[str]
OptionalStrDict = typing.Optional[typing.Union[str, dict[str, str]]]
OptionalBool = typing.Optional[bool]


@functools.singledispatch
def is_gff3(f) -> bool:
    """True if gff-version is 3"""
    raise TypeError(f"unsopported type type {type(f)}")


@is_gff3.register
def _(f: PathType) -> bool:
    gff3_header = "gff-version 3"
    with open_(f) as f:
        is_gff3 = gff3_header in f.readline()
        f.seek(0)
    return is_gff3


@is_gff3.register
def _(f: io.IOBase) -> bool:
    gff3_header = "gff-version 3"
    is_gff3 = gff3_header in f.readline()
    f.seek(0)
    return is_gff3


@is_gff3.register
def _(f: typing.Union[list, tuple]) -> bool:
    gff3_header = "gff-version 3"
    return f and gff3_header in f[0]


_gff_item_map = {
    "type": "biotype",
    "end": "stop",
    "attributes": "attrs",
}


class GffRecord:
    """Container for single GFF record. Elements can be indexed as a dict
    or directly as attributes.

    Notes
    -----
    Some fields defined in the GFF spec are aliased, or truncated.
    """

    # should be a dataclass, but in py 3.9 they don't support slots
    # so until then ...
    __slots__ = (
        "seqid",
        "source",
        "biotype",
        "name",
        "parent_id",
        "start",
        "stop",
        "spans",
        "strand",
        "score",
        "phase",
        "attrs",
        "comments",
    )

    def __init__(
        self,
        seqid: OptionalStr = None,
        source: OptionalStr = None,
        biotype: OptionalStr = None,
        name: OptionalStr = None,
        parent_id: OptionalStr = None,
        start: OptionalInt = None,
        stop: OptionalInt = None,
        spans: OptionalIntList = None,
        score: OptionalStr = None,
        strand: OptionalStr = None,
        phase: OptionalInt = None,
        comments: OptionalStr = None,
        attrs: OptionalStrDict = None,
    ):
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

    def __repr__(self):
        name = self.__class__.__name__
        return (
            f"{name}(seqid={self.seqid!r}, name={self.name!r}, "
            f"start={self.start}, stop={self.stop}, strand={self.strand!r}, "
            f"spans={self.spans})"
        )

    def __getitem__(self, item: str):
        item = item.lower()
        return getattr(self, _gff_item_map.get(item, item))

    def __setitem__(self, key: str, value: typing.Any):
        key = key.lower()
        setattr(self, _gff_item_map.get(key, key), value)

    def update(self, values: dict[str, typing.Any]):
        for key, value in values.items():
            self[key] = value

    def get(self, item: str) -> typing.Any:
        item = item.lower()
        return getattr(self, _gff_item_map.get(item, item), None)

    def to_dict(self) -> dict[str, typing.Any]:
        return {k: getattr(self, k) for k in self.__slots__}


@functools.singledispatch
def gff_parser(
    f: typing.Union[PathType, typing.IO, OptionalStrContainer],
    attribute_parser: OptionalCallable = None,
    seqids: OptionalStrContainer = None,
    gff3: OptionalBool = None,
) -> typing.Iterable[GffRecord]:
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

    Returns
    -------
    dict
        contains each of the 9 parameters specified by gff3, and comments.
    """
    gff3 = gff3 or is_gff3(f)
    yield from _gff_parser(
        f, attribute_parser=attribute_parser, seqids=seqids, gff3=gff3
    )


@gff_parser.register
def _(
    f: PathType,
    attribute_parser: OptionalCallable = None,
    seqids: OptionalStrContainer = None,
    gff3: OptionalBool = None,
) -> typing.Iterable[GffRecord]:

    with open_(f) as infile:
        yield from gff_parser(
            infile, attribute_parser=attribute_parser, seqids=seqids, gff3=gff3
        )


def _gff_parser(
    f: typing.Union[PathType, typing.IO, OptionalStrContainer],
    attribute_parser: OptionalCallable = None,
    seqids: OptionalStrContainer = None,
    gff3: bool = True,
) -> typing.Iterable[GffRecord]:
    """parses a gff file"""
    seqids = seqids or set()
    seqids = {seqids} if isinstance(seqids, str) else set(seqids)

    if attribute_parser is None:
        attribute_parser = parse_attributes_gff3 if gff3 else parse_attributes_gff2

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
        seqid, source, type_, start, end, score, strand, phase, attributes = cols

        if seqids and seqid not in seqids:
            continue

        # adjust for 0-based indexing
        start, end = int(start) - 1, int(end)
        # start is always meant to be less than end in GFF
        # features that extend beyond sequence have negative indices
        if start < 0 or end < 0:
            start, end = abs(start), abs(end)
        if start > end:
            start, end = end, start

        # all attributes have an "ID" but this may not be unique
        attributes = attribute_parser(attributes, (start, end))

        yield GffRecord(
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


def parse_attributes_gff2(attributes: str, span: typing.Tuple[int, int]) -> dict:
    """Returns a dict with name and info keys"""
    name = attributes[attributes.find('"') + 1 :]
    if '"' in name:
        name = name[: name.find('"')]
    return {"ID": name, "Info": attributes}


def parse_attributes_gff3(attributes: str, span: typing.Tuple[int, int]) -> dict:
    """Returns a dictionary containing all the attributes"""
    attributes = attributes.strip(";")
    attributes = attributes.split(";")
    attributes = dict(t.split("=") for t in attributes) if attributes[0] else {}
    if "Parent" in attributes:
        # There may be multiple parents
        if "," in attributes["Parent"]:
            attributes["Parent"] = attributes["Parent"].split(",")
        else:
            attributes["Parent"] = [attributes["Parent"]]
    attributes["ID"] = attributes.get("ID", "")
    return attributes


def merged_gff_records(records: list[GffRecord], num_fake_ids: int) -> tuple[dict, int]:
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

    reduced = collections.OrderedDict()
    # collapse records with ID's occurring in multiple rows into one
    # row, converting their coordinates
    # extract the name from ID and add this into the table
    # I am not convinced we can rely on gff files to be ordered,
    # if we could, we could do this as one pass over the data
    # all keys need to be lower case
    for record in records:
        attrs = record.attrs or ""
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

        # should this just be an append?
        reduced[record_id].spans.append((record["start"], record["stop"]))

    del records

    return reduced, num_fake_ids
