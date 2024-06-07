import pytest

from cogent3.core import new_moltype, new_sequence


@pytest.mark.parametrize("name", ("dna", "rna", "protein", "protein_with_stop", "text"))
def test_moltype_make_seq(name):
    raw = "ACGGA"
    moltype = new_moltype.get_moltype(name)
    seq = moltype.make_seq(name="s1", seq=raw)
    assert seq.moltype.name == name
    assert str(seq) == raw


def test_moltype_make_bytes_seq():
    raw = "ACGGA"
    name = "bytes"
    moltype = new_moltype.get_moltype(name)
    seq = moltype.make_seq(name="s1", seq=raw)
    assert seq.moltype.name == name
    assert str(seq) == raw


# From test_sequence.py
def test_seqview_seqid():
    sv = new_sequence.SeqView(seq="ACGGTGGGAC")
    assert sv.seqid is None

    sv = new_sequence.SeqView(seq="ACGGTGGGAC", seqid="seq1")
    assert sv.seqid == "seq1"


def test_seqview_rich_dict_round_trip_seqid():
    sv = new_sequence.SeqView(seq="ACGGTGGGAC", seqid="seq1")
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seqid"] == "seq1"

    got = new_sequence.SeqView.from_rich_dict(rd)
    assert got.seqid == "seq1"

    sv = new_sequence.SeqView(seq="ACGGTGGGAC")
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seqid"] is None

    got = new_sequence.SeqView.from_rich_dict(rd)
    assert got.seqid is None


def test_seqview_slice_propagates_seqid():
    sv = new_sequence.SeqView(seq="ACGGTGGGAC", seqid="seq1")
    sliced_sv = sv[1:8:2]
    assert sliced_sv.seqid == "seq1"

    copied_sv = sliced_sv.copy(sliced=False)
    assert copied_sv.seqid == "seq1"

    copied_sliced_sv = sliced_sv.copy(sliced=True)
    assert copied_sliced_sv.seqid == "seq1"


def test_sequences_propogates_seqid():
    # creating a name Sequence propagates the seqid to the SeqView.
    seq = new_moltype.DNA.make_seq(seq="ACGGTGGGAC", name="seq1")
    assert seq._seq.seqid == "seq1"

    # renaming the Sequence doesnt change the seqid of the SeqView.
    seq.name = "seq2"
    assert seq.name == "seq2"
    assert seq._seq.seqid == "seq1"


def test_make_seq_assigns_to_seqview():
    seq = new_moltype.DNA.make_seq(seq="ACGT", name="s1")
    assert seq.name == seq._seq.seqid == "s1"


def test_empty_seqview_translate_position():
    sv = new_sequence.SeqView(seq="")
    assert sv.absolute_position(0) == 0
    assert sv.relative_position(0) == 0


@pytest.mark.parametrize("start", (None, 0, 1, 10, -1, -10))
@pytest.mark.parametrize("stop", (None, 10, 8, 1, 0, -1, -11))
@pytest.mark.parametrize("step", (None, 1, 2, -1, -2))
@pytest.mark.parametrize("length", (1, 8, 999))
def test_seqview_seq_len_init(start, stop, step, length):
    # seq_len is length of seq when None
    seq_data = "A" * length
    sv = new_sequence.SeqView(seq=seq_data, start=start, stop=stop, step=step)
    expect = len(seq_data)
    # Check property and slot
    assert sv.seq_len == expect
    assert sv._seq_len == expect


@pytest.mark.parametrize("seq, seq_len", [("A", 0), ("", 1), ("A", 2)])
def test_seqview_seq_len_mismatch(seq, seq_len):
    # If provided, seq_len must match len(seq)
    with pytest.raises(AssertionError):
        new_sequence.SeqView(seq=seq, seq_len=seq_len)


def test_seqview_copy_propagates_seq_len():
    seq = "ACGGTGGGAC"
    sv = new_sequence.SeqView(seq=seq)
    copied = sv.copy()
    assert copied.seq_len == len(seq)


def test_seqview_seq_len_modified_seq():
    seq = "ACGGTGGGAC"
    sv = new_sequence.SeqView(seq=seq)

    sv.seq = "ATGC"  # this should not modify seq_len
    assert sv.seq_len == len(seq)
