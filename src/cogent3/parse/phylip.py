#!/usr/bin/env python
from cogent3.core.alignment import Alignment
from cogent3.parse.record import RecordError


__author__ = "Micah Hamady"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Micah Hamady", "Peter Maxwell", "Gavin Huttley", "Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Micah Hamady"
__email__ = "hamady@colorado.edu"
__status__ = "Prototype"


def is_blank(x):
    """Checks if x is blank."""
    return not x.strip()


def _get_header_info(line):
    """
    Get number of sequences and length of sequence
    """
    header_parts = line.split()
    num_seqs, length = list(map(int, header_parts[:2]))
    is_interleaved = len(header_parts) > 2
    return num_seqs, length, is_interleaved


def _split_line(line, id_offset):
    """
    First 10 chars must be blank or contain id info
    """
    if not line or not line.strip():
        return None, None

    # extract id and sequence
    curr_id = line[0:id_offset].strip()
    curr_seq = line[id_offset:].strip().replace(" ", "")

    return curr_id, curr_seq


def MinimalPhylipParser(data, id_map=None, interleaved=True):
    """Yields successive sequences from data as (label, seq) tuples.

    **Need to implement id map.

    **NOTE if using phylip interleaved format, will cache entire file in
        memory before returning sequences. If phylip file not interleaved
        then will yield each successive sequence.

    data: sequence of lines in phylip format (an open file, list, etc)
    id_map: optional id mapping from external ids to phylip labels - not sure
        if we're going to implement this


    returns (id, sequence) tuples
    """

    seq_cache = {}
    interleaved_id_map = {}
    id_offset = 10
    curr_ct = -1

    for line in data:
        if curr_ct == -1:
            # get header info
            num_seqs, seq_len, interleaved = _get_header_info(line)

            if not num_seqs or not seq_len:
                return
            curr_ct += 1
            continue

        curr_id, curr_seq = _split_line(line, id_offset)

        # skip blank lines
        if not curr_id and not curr_seq:
            continue

        if not interleaved:
            if curr_id:
                if seq_cache:
                    yield seq_cache[0], "".join(seq_cache[1:])
                seq_cache = [curr_id, curr_seq]
            else:
                seq_cache.append(curr_seq)
        else:
            curr_id_ix = curr_ct % num_seqs

            if (curr_ct + 1) % num_seqs == 0:
                id_offset = 0

            if curr_id_ix not in interleaved_id_map:
                interleaved_id_map[curr_id_ix] = curr_id
                seq_cache[curr_id_ix] = []

            seq_cache[curr_id_ix].append(curr_seq)
        curr_ct += 1

    # return joined sequences if interleaved
    if interleaved:
        for curr_id_ix, seq_parts in list(seq_cache.items()):
            join_seq = "".join(seq_parts)

            if len(join_seq) != seq_len:
                raise RecordError(
                    "Length of sequence '%s' is not the same as in header "
                    "Found %d, Expected %d"
                    % (interleaved_id_map[curr_id_ix], len(join_seq), seq_len)
                )

            yield interleaved_id_map[curr_id_ix], join_seq
    # return last seq if not interleaved
    else:
        if seq_cache:
            yield seq_cache[0], "".join(seq_cache[1:])


def get_align_for_phylip(data, id_map=None):
    """
    Convenience function to return aligment object from phylip data

    data: sequence of lines in phylip format (an open file, list, etc)
    id_map: optional id mapping from external ids to phylip labels - not sure
        if we're going to implement this

    returns Alignment object
    """

    mpp = MinimalPhylipParser(data, id_map)

    tuples = []
    for tup in mpp:
        tuples.append(tup)
    return Alignment(tuples)
