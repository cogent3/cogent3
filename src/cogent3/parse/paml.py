#!/usr/bin/env python
from io import TextIOWrapper


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def PamlParser(data):
    if isinstance(data, TextIOWrapper):
        data = data.read().splitlines()
    num_seqs, seq_len = [int(v) for v in data[0].split()]
    curr_seq = []
    curr_length = 0
    seqname = None
    n = 0
    for line in data[1:]:
        line = line.strip()
        if not line:
            continue

        if seqname is None:
            seqname = line
            continue

        curr_length += len(line)
        curr_seq.append(line)
        if curr_length == seq_len:
            yield seqname, "".join(curr_seq)

            seqname = None
            curr_seq = []
            curr_length = 0
            n += 1

    if n != num_seqs:
        raise ValueError(f"read {n} seqs, expected {num_seqs}")
