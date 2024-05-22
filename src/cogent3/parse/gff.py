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
def gff_parser(
    f: typing.Union[PathType, typing.IO, OptionalStrContainer],
    attribute_parser: OptionalCallable = None,
    seqids: OptionalStrContainer = None,
) -> typing.Iterable[dict]:
    """parses a gff file

    Parameters
    -----------
    f
        accepts string path or pathlib.Path or file-like object (e.g. StringIO)
    attribute_parser
        callback function for custom handling of attributes field. Defaults to
        a function that converts Attributes field into a dict based on the gff
        version.

    Returns
    -------
    dict
        contains each of the 9 parameters specified by gff3, and comments.
    """
    yield from _gff_parser(f, attribute_parser=attribute_parser, seqids=seqids)


@gff_parser.register
def _(
    f: str,
    attribute_parser: OptionalCallable = None,
    seqids: OptionalStrContainer = None,
) -> typing.Iterable[dict]:
    with open_(f) as infile:
        yield from gff_parser(infile, attribute_parser=attribute_parser, seqids=seqids)


@gff_parser.register
def _(
    f: PathType,
    attribute_parser: OptionalCallable = None,
    seqids: OptionalStrContainer = None,
) -> typing.Iterable[dict]:
    with open_(f) as infile:
        yield from gff_parser(infile, attribute_parser=attribute_parser, seqids=seqids)


def _gff_parser(
    f: typing.Union[PathType, typing.IO, OptionalStrContainer],
    attribute_parser: OptionalCallable = None,
    seqids: OptionalStrContainer = None,
    """parses a gff file"""
    seqids = seqids or set()
    seqids = {seqids} if isinstance(seqids, str) else set(seqids)

    gff3_header = "gff-version 3"
    if isinstance(f, list):
        gff3 = f and gff3_header in f[0]
    else:
        gff3 = gff3_header in f.readline()
        f.seek(0)

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
        assert len(cols) == 9, len(line)
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
        # reverse indices when the feature is on the opposite strand
        if strand == "-":
            (start, end) = (end, start)

        # all attributes have an "ID" but this may not be unique
        attributes = attribute_parser(attributes, (start, end))

        yield {
            "SeqID": seqid,
            "Source": source,
            "Type": type_,
            "Start": start,
            "End": end,
            "Score": score,
            "Strand": strand,
            "Phase": phase,
            "Attributes": attributes,
            "Comments": comments,
        }


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
