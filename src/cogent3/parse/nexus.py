#!/usr/bin/env python

"""
parses Nexus formatted tree files and Branchlength info in log files
"""

import re
from collections import defaultdict

from cogent3.parse.record import RecordError
from cogent3.util.io import open_

strip = str.strip


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
            if name.endswith(","):
                name = name[:-1]
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
