#!/usr/bin/env python
from io import TextIOWrapper


def PamlParser(data):
    if isinstance(data, TextIOWrapper):
        data = data.read().splitlines()

    data = list(data)

    num_seqs, seq_len = [int(v) for v in data.pop(0).split()]
    curr_seq = []
    curr_length = 0
    seqname = None
    n = 0
    for line in data:
        line = line.strip()
        if not line:
            continue

        if seqname is None:
            seqname = line
            continue

        curr_length += len(line)
        curr_seq.append(line)
        if curr_length == seq_len:
            yield seqname, "".join(curr_seq).upper()

            seqname = None
            curr_seq = []
            curr_length = 0
            n += 1

    if n != num_seqs:
        raise ValueError(f"read {n} seqs, expected {num_seqs}")
