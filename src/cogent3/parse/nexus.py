#!/usr/bin/env python

"""
parses Nexus formatted tree files and Branchlength info in log files
"""

import re
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from functools import singledispatch
from pathlib import Path
from typing import TypedDict, cast

from scinexus.io_util import iter_splitlines, open_

from cogent3.parse.record import RecordError

strip = str.strip


@dataclass(frozen=True, slots=True)
class Partition:
    """a single Nexus charset, resolved for slicing an alignment

    Ranges use zero-based, half-open slices so they apply directly to an
    alignment. The class holds only the information needed to obtain and
    match an alignment, never an alignment instance itself.
    """

    name: str
    ranges: tuple[slice, ...]
    data_type: str | None = None
    model: str | None = None
    align_name: str | None = None
    tree_len: float | None = None


@dataclass(frozen=True, slots=True)
class PartitionSet:
    """an ordered set of Partitions defined by a Nexus charpartition"""

    name: str
    partitions: tuple[Partition, ...]

    def __len__(self) -> int:
        return len(self.partitions)

    def __iter__(self) -> Iterator[Partition]:
        return iter(self.partitions)

    def __getitem__(self, name: str) -> Partition:
        for partition in self.partitions:
            if partition.name == name:
                return partition
        raise KeyError(name)


def range_to_slice(text: str) -> slice:
    """converts a Nexus charset range token to a zero-based, half-open slice

    Nexus ranges are 1-based with an inclusive end. A trailing "\\k" gives a
    step, "." as the end means "to the end", and "*" means the whole file.
    """
    if text == "*":
        return slice(0, None)

    step = None
    if "\\" in text:
        text, step_s = text.split("\\", 1)
        step = int(step_s)

    if "-" not in text:
        start = int(text)
        return slice(start - 1, start)

    start_s, end_s = text.split("-", 1)
    start = int(start_s) - 1
    stop = None if end_s == "." else int(end_s)
    return slice(start, stop, step)


_type_prefix = re.compile(r"^\s*([A-Za-z][A-Za-z0-9]*)\s*,\s*(.*)$", re.DOTALL)


def parse_charset(text: str) -> Partition:
    """parses one Nexus charset statement into a Partition

    Expects the statement body with the "charset" keyword and trailing ";"
    already removed, e.g. "part1 = aln.phy:CODON, 1-900".
    """
    name, spec = text.split("=", 1)
    name = name.strip()
    spec = spec.strip()

    align_name = None
    if ":" in spec:
        file_part, spec = spec.split(":", 1)
        align_name = file_part.strip() or None

    data_type = None
    if match := _type_prefix.match(spec):
        data_type, spec = match.groups()

    tokens = re.split(r"[\s,]+", spec.strip())
    ranges = tuple(range_to_slice(token) for token in tokens if token)
    return Partition(
        name=name, ranges=ranges, data_type=data_type, align_name=align_name
    )


_data_type_kw = re.compile(r"^(DNA|AA|BIN|MORPH|CODON\d*|NT2AA\d*)$")


class _Slot(TypedDict):
    model: str | None
    data_type: str | None
    tree_len: float | None


def parse_charpartition(text: str) -> tuple[str, dict[str, _Slot]]:
    """parses a Nexus charpartition into a scheme name and per-charset slots

    Expects the statement body with the "charpartition" keyword and trailing
    ";" removed. Each slot maps a charset to a model (or, when the slot is a
    data-type keyword, a data type) and an optional {value} branch length.
    """
    scheme, body = text.split("=", 1)
    scheme = scheme.strip()

    mapping: dict[str, _Slot] = {}
    for entry in body.split(","):
        slot, charset = entry.split(":", 1)
        slot = slot.strip()
        charset = charset.strip()

        tree_len: float | None = None
        if "{" in slot:
            slot, value = slot.split("{", 1)
            tree_len = float(value.rstrip("}"))
            slot = slot.strip()

        if _data_type_kw.match(slot):
            mapping[charset] = {"model": None, "data_type": slot, "tree_len": tree_len}
        else:
            mapping[charset] = {"model": slot, "data_type": None, "tree_len": tree_len}
    return scheme, mapping


@singledispatch
def get_sets_block(data: str | Path | Iterable[str]) -> str | None:
    """returns the text of the Nexus 'sets' block with [ ... ] comments removed"""
    in_block = False
    collected = []
    for line in cast("Iterable[str]", data):
        stripped = line.lower().lstrip()
        if stripped.startswith("begin sets;"):
            in_block = True
            continue
        if in_block:
            if stripped.startswith(("end;", "endblock;")):
                break
            collected.append(line)

    if not collected:
        return None
    return re.sub(r"\[.*?\]", " ", "".join(collected), flags=re.DOTALL)


@get_sets_block.register(str)
@get_sets_block.register(Path)
def _(data: str | Path) -> str | None:
    return get_sets_block(iter_splitlines(data))


def parse_nexus_partitions(data: str | Path | Iterable[str]) -> PartitionSet:
    """parses the partition scheme from the 'sets' block of a Nexus file"""
    block = get_sets_block(data)
    if block is None:
        msg = "no sets block found"
        raise RecordError(msg)

    charsets: dict[str, Partition] = {}
    scheme_name = ""
    slots: dict[str, _Slot] = {}
    for raw in block.split(";"):
        statement = raw.strip()
        if not statement:
            continue
        keyword, rest = statement.split(None, 1)
        keyword = keyword.lower()
        if keyword == "charset":
            charset = parse_charset(rest)
            charsets[charset.name] = charset
        elif keyword == "charpartition":
            scheme_name, slots = parse_charpartition(rest)

    partitions = []
    for name, slot in slots.items():
        base = charsets[name]
        partitions.append(
            Partition(
                name=base.name,
                ranges=base.ranges,
                data_type=base.data_type or slot["data_type"],
                model=slot["model"],
                align_name=base.align_name,
                tree_len=slot["tree_len"],
            )
        )
    return PartitionSet(name=scheme_name, partitions=tuple(partitions))


def parse_nexus_tree(tree_f):
    """returns a dict mapping taxa # to name from the translation table,
    and a dict mapping tree name to dnd string;
    takes a handle for a Nexus formatted file as input"""
    trans_table = None
    tree_info = get_tree_info(tree_f)
    check_tree_info(tree_info)
    header_s, trans_table_s, dnd_s = split_tree_info(tree_info)
    if trans_table_s:
        trans_table = parse_trans_table(trans_table_s)
    dnd = parse_dnd(dnd_s)
    return trans_table, dnd


def get_tree_info(tree_f):
    """returns the trees section of a Nexus file:
    takes a handle for a Nexus formatted file as input:
    returns the section describing trees as a list of strings"""
    in_tree = False
    result = []
    for line in tree_f:
        # get lines from the 'Begin trees;' tag to the 'End;' tag
        line_lower = line.lower()
        if line_lower.startswith("begin trees;"):
            in_tree = True
        if in_tree:
            if line_lower.startswith(("end;", "endblock;")):
                return result
            result.append(line)
    return None


def check_tree_info(tree_info) -> None:
    """makes sure that there is a tree section in the file"""
    if tree_info:
        pass
    else:
        msg = "not a valid Nexus Tree File"
        raise RecordError(msg)


def split_tree_info(tree_info):
    """Returns header, table, and dnd info from tree section of Nexus file.:

    Expects to receive the output of get_tree_info"""
    header = []
    trans_table = []
    dnd = []
    state = "in_header"

    for line in tree_info:
        line_lower = line.lower()
        if state == "in_header":
            header.append(line)
            if line_lower.strip() == "translate":
                state = "in_trans"
            elif line_lower.startswith("tree"):
                state = "in_dnd"
                dnd.append(line)

        elif state == "in_trans":
            trans_table.append(line)
            if line.strip() == ";":
                state = "in_dnd"

        elif state == "in_dnd":
            dnd.append(line)
    return header, trans_table, dnd


def parse_trans_table(trans_table):
    """returns a dict with the taxa names indexed by number"""
    result = {}
    for line in trans_table:
        line = line.strip()
        if line != ";":
            label, name = line.split(None, 1)
            # take comma out of name if it is there
            name = name.removesuffix(",")
            # remove single quotes
            if name.startswith("'") and name.endswith("'"):
                name = name[1:-1]
            result[label] = name
    return result


def parse_dnd(dnd):  # get rooted info
    """returns a dict with dnd indexed by name"""
    dnd_dict = {}
    for line in dnd:
        line = line.strip()
        name, dnd_s = list(map(strip, line.split("=", 1)))
        # get dnd from dnd_s and populate
        dnd_index = dnd_s.find("(")
        data = dnd_s[dnd_index:]
        dnd_dict[name] = data
    return dnd_dict


def get_BL_table(branch_lengths):
    """returns the section of the log file with the BL table
    as a list of strings"""

    in_table = 0
    result = []
    beg_tag = re.compile(r"\s+Node\s+to node\s+length")
    end_tag = re.compile("Sum")
    for line in branch_lengths:
        if end_tag.match(line):
            in_table = 0
        if beg_tag.match(line):
            in_table = 1
        if in_table == 1 and (
            not line.startswith("---")
            and not beg_tag.match(line)
            and line.strip() != ""
        ):
            result.append(line)
    return result


def find_fields(line, field_order=None, field_delims=None):
    """takes line from BL table and returns dict with field names mapped to info

    field order is the order of field names to extract from the file and
    field_delims is a list of index numbers indicating where the field is split
    """

    field_order = field_order or ["taxa", "parent", "bl"]
    field_delims = field_delims or [0, 21, 36, 49]

    field_dict = {}
    for i, f in enumerate(field_order):
        start = field_delims[i]
        try:
            end = field_delims[i + 1]
        except IndexError:
            end = None
        field_dict[f] = line[start:end].strip()
    return field_dict


def parse_taxa(taxa_field):
    """gets taxa # from taxa field extracted with find_fields"""

    if not (term_match := re.search(r"\(\d+\)", taxa_field)):
        return taxa_field
    term = term_match[0]
    data_match = re.search(r"\d+", term)
    return data_match[0]


def parse_PAUP_log(branch_lengths):
    """gets branch length info from a PAUP log file
    returns a dictionary mapping the taxon number to the parent number
    and the branch length"""
    BL_table = get_BL_table(branch_lengths)
    BL_dict = {}
    for line in BL_table:
        info = find_fields(line)
        parent = info["parent"]
        bl = float(info["bl"])
        taxa = parse_taxa(info["taxa"])

        BL_dict[taxa] = (parent, bl)

    return BL_dict


def MinimalNexusAlignParser(align_path):
    """returns {label: seq, ...}"""
    infile = open_(align_path) if type(align_path) == str else align_path

    isblock = re.compile(r"begin\s+(data|characters)").search
    inblock = False
    try:
        line = infile.readline()
    except AttributeError:
        # guessing it's a list of strings from a nexus file
        line = next(iter(infile))

    if not line.lower().startswith("#nexus"):
        msg = "not a nexus file"
        raise ValueError(msg)

    block = []
    index = None
    for line in infile:
        if isblock(line.lower()):
            inblock = True
        elif inblock and line.lower().startswith("end;"):
            break
        elif inblock:
            line = line.strip()
            if line.lower().startswith("matrix"):
                index = len(block)
            elif not line.startswith(";"):
                block.append(line)

    if hasattr(infile, "close"):
        infile.close()

    if not block:
        msg = "not found DATA or CHARACTER block"
        raise ValueError(msg)
    elif index is None:
        msg = "malformed block, no 'matrix' line"
        raise RecordError(msg)

    block = block[index:]
    seqs = defaultdict(list)
    for line in block:
        if not line or (line.startswith("[") and line.endswith("]")):
            # blank or comment line
            continue

        line = line.split()
        seqs[line[0]].append("".join(line[1:]))

    for n, s in seqs.items():
        yield n, "".join(s)
