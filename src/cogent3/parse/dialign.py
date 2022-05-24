#!/usr/bin/env python

import re

from cogent3 import ASCII


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

_header = re.compile(r"^\s+[=]+")
_quality_scores = re.compile(r"^ +\d+[\s\d]*$")


def align_block_lines(lines):
    counter = 0
    for line in lines:
        if "Alignment (DIALIGN format):" in line:
            counter += 1
            continue
        elif counter == 1 and _header.findall(line):
            counter += 2
            continue
        elif not counter or not line:
            continue
        elif "Sequence tree:" in line:
            break
        yield line


def parse_data_line(line):
    if _quality_scores.findall(line):
        line = line.split()
        name = None
        seq = "".join(line)
    elif line[0].isspace():
        name, seq = None, None
    else:
        line = line.split()
        name = line[0]
        seq = "".join(line[2:])
    return name, seq


def DialignParser(lines, seq_maker=None, get_scores=False):
    """Yields label, sequence pairs.

    The alignment quality info is recorded in the sequence
    case and the score line. Font info can be handled by
    providing a custom seq_maker function. The quality
    scores are returned as the last value pair with
    name 'QualityScores' when get_scores is True."""

    if seq_maker is None:
        seq_maker = ASCII.make_seq
    seqs = {}
    quality_scores = []
    for line in align_block_lines(lines):
        name, seq = parse_data_line(line)
        if seq is None:
            continue
        elif name is None and seq:
            quality_scores.append(seq)
            continue
        if name in seqs:
            seqs[name].append(seq)
        else:
            seqs[name] = [seq]

    # concat sequence blocks
    for name, seq_segs in list(seqs.items()):
        seq = "".join(seq_segs)
        yield name, seq_maker(seq, name=name)

    if get_scores:
        yield "QualityScores", "".join(quality_scores)
