import json
import re

import numpy
import pytest

import cogent3

from cogent3._version import __version__
from cogent3.core import (
    new_alphabet,
    new_genetic_code,
    new_moltype,
    new_sequence,
)
from cogent3.util.deserialise import deserialise_object
from cogent3.util.misc import get_object_provenance


@pytest.fixture(scope="function")
def dna_alphabet():
    return new_moltype.DNA.degen_gapped_alphabet


@pytest.fixture(scope="function")
def ascii_alphabet():
    return new_moltype.ASCII.alphabet


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


# Tests from test_sequence.py
@pytest.mark.xfail(
    reason="AttributeError: 'MolType' object has no attribute 'coerce_str'"
)
def test_to_moltype_dna():
    """to_moltype("dna") ensures conversion from T to U"""
    seq = new_moltype.DNA.make_seq(seq="AAAAGGGGTTT", name="seq1")
    rna = seq.to_moltype("rna")

    assert "T" not in rna


@pytest.mark.xfail(
    reason="AttributeError: 'MolType' object has no attribute 'coerce_str'"
)
def test_to_moltype_rna():
    """to_moltype("rna") ensures conversion from U to T"""
    seq = new_moltype.RNA.make_seq(seq="AAAAGGGGUUU", name="seq1")
    rna = seq.to_moltype("dna")

    assert "U" not in rna


def test_to_rich_dict():
    """Sequence to_dict works"""
    seq = "AAGGCC"
    r = new_moltype.DNA.make_seq(seq="AAGGCC", name="seq1")
    got = r.to_rich_dict()
    seq = new_sequence.SeqView(seq=seq, seqid="seq1").to_rich_dict()

    expect = {
        "name": "seq1",
        "seq": seq,
        "moltype": r.moltype.label,
        "info": None,
        "type": get_object_provenance(r),
        "version": __version__,
        "annotation_offset": 0,
    }

    assert got == expect


@pytest.mark.xfail(reason="refactor: how to serialise an alphabet")
def test_to_json():
    """to_json roundtrip recreates to_dict"""
    dna = new_moltype.DNA
    r = dna.make_seq(seq="AAGGCC", name="seq1")
    got = json.loads(r.to_json())
    seq = new_sequence.SeqView(
        seq="AAGGCC", seqid="seq1", alphabet=dna.most_degen_alphabet()
    ).to_rich_dict()

    expect = {
        "name": "seq1",
        "seq": seq,
        "moltype": dna.label,
        "info": None,
        "type": get_object_provenance(r),
        "version": __version__,
        "annotation_offset": 0,
    }

    assert got == expect


@pytest.mark.xfail(
    reason="NotImplementedError: deserialising 'cogent3.core.new_sequence.DnaSequence' from json"
)
def test_offset_with_multiple_slices(DATA_DIR):
    seq = new_moltype.DNA.make_seq(
        seq="ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC", name="22"
    )
    gff3_path = DATA_DIR / "ensembl_sample.gff3"
    # TODO: directly assign an annotation_db, annotate_from_gff to be discontinued
    seq.annotate_from_gff(gff3_path)
    rd = seq[2:].to_rich_dict()
    s1 = deserialise_object(rd)
    assert s1.annotation_offset == 2
    rd = s1[3:].to_rich_dict()
    s2 = deserialise_object(rd)
    assert s2.annotation_offset == 5
    expect = {(f.seqid, f.biotype, f.name) for f in seq.get_features(start=5)}
    got = {(f.seqid, f.biotype, f.name) for f in s2.get_features()}
    assert got == expect


@pytest.mark.parametrize("coord", ("start", "stop"))
def test_seqview_to_rich_dict(coord, dna_alphabet):
    parent = "ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC"
    sv = new_sequence.SeqView(seq=parent, alphabet=dna_alphabet)
    plus = sv.to_rich_dict()
    minus = sv[::-1].to_rich_dict()
    plus = plus.pop("init_args")
    minus = minus.pop("init_args")
    assert plus.pop("seq") == minus.pop("seq")
    assert plus["step"] == -minus["step"]
    assert coord not in plus
    assert coord not in minus


@pytest.mark.xfail(
    reason="NotImplementedError: deserialising 'cogent3.core.new_sequence.SeqView' from json"
)
@pytest.mark.parametrize("reverse", (False, True))
def test_seqview_round_trip(reverse, dna_alphabet):
    parent = "ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC"
    sv = new_sequence.SeqView(seq=parent, alphabet=dna_alphabet)
    sv = sv[::-1] if reverse else sv

    rd = sv.to_rich_dict()
    got = deserialise_object(rd)
    assert isinstance(got, new_sequence.SeqView)
    assert got.to_rich_dict() == sv.to_rich_dict()


@pytest.mark.parametrize("reverse", (False, True))
def test_sliced_seqview_rich_dict(reverse, dna_alphabet):
    parent = "ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC"
    sl = slice(2, 13)
    sv = new_sequence.SeqView(seq=parent, alphabet=dna_alphabet)[sl]
    sv = sv[::-1] if reverse else sv
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seq"] == parent[sl]
    assert rd["init_args"]["offset"] == 2


@pytest.mark.parametrize(
    "sl",
    (
        slice(2, 5, 1),  # positive indices, positive step
        slice(-8, -5, 1),  # negative indices, positive step
        slice(4, 1, -1),  # positive indices, negative step
        slice(-6, -9, -1),  # negative indices, negative step
    ),
)
@pytest.mark.parametrize("offset", (4, 0))
def test_parent_start_stop(sl, offset, ascii_alphabet):
    data = "0123456789"
    # check our slice matches the expectation for rest of test
    expect = "234" if sl.step > 0 else "432"
    sv = new_sequence.SeqView(seq=data, alphabet=ascii_alphabet)
    sv.offset = offset
    sv = sv[sl]
    assert sv.str_value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.parent_start, sv.parent_stop) == (2 + offset, 5 + offset)


@pytest.mark.parametrize(
    "sl",
    (
        slice(None, None, 1),  # slice whole sequence plus strand
        slice(None, None, -1),  # slice whole sequence minus strand
    ),
)
def test_parent_start_stop_limits(sl, ascii_alphabet):
    data = "0123456789"
    # check our slice matches the expectation for rest of test
    expect = data[sl]
    sv = new_sequence.SeqView(seq=data, alphabet=ascii_alphabet)
    sv = sv[sl]
    assert sv.str_value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.parent_start, sv.parent_stop) == (0, 10)


@pytest.mark.parametrize("rev", (False, True))
def test_parent_start_stop_empty(rev, ascii_alphabet):
    data = "0123456789"
    # check our slice matches the expectation for rest of test
    expect = ""
    sv = new_sequence.SeqView(seq=data, alphabet=ascii_alphabet)
    sv = sv[0 : 0 : -1 if rev else 1]
    assert sv.str_value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.parent_start, sv.parent_stop) == (0, 0)


@pytest.mark.parametrize("rev", (False, True))
@pytest.mark.parametrize("index", range(9))
def test_parent_start_stop_singletons(index, rev, ascii_alphabet):
    data = "0123456789"
    start, stop = (-(10 - index), -(10 - index + 1)) if rev else (index, index + 1)
    sl = slice(start, stop, -1 if rev else 1)
    # check our slice matches the expectation for rest of test
    expect = data[sl]
    sv = new_sequence.SeqView(seq=data, alphabet=ascii_alphabet)
    sv = sv[sl]
    assert sv.str_value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.parent_start, sv.parent_stop) == (index, index + 1)


def test_get_drawable(DATA_DIR):
    seq = cogent3.load_seq(DATA_DIR / "annotated_seq.gb")
    seq = seq[2000:4000]
    biotypes = "CDS", "gene", "mRNA"
    for feat in seq.get_features(biotype=biotypes, allow_partial=True):
        draw = feat.get_drawable()
        assert "(incomplete)" in draw.text

    full = seq.get_drawable(biotype=biotypes)
    # should only include elements that overlap the segment
    assert len(full.traces) == len(biotypes)
    # and their names should indicate they're incomplete
    for trace in full.traces:
        assert "(incomplete)" in trace.text


@pytest.mark.parametrize("gc,seq", ((1, "TCCTGA"), (1, "ACGTAA---"), (2, "TCCAGG")))
def test_has_terminal_stop_true(gc, seq):
    gc = new_genetic_code.get_code(gc)
    seq = new_moltype.DNA.make_seq(seq=seq)
    assert seq.has_terminal_stop(gc=gc)


@pytest.mark.parametrize(
    "gc,seq", ((1, "TCCAGG"), (2, "TCCAAA"), (1, "CCTGA"), (2, "CCAGG"))
)
def test_has_terminal_stop_false(gc, seq):
    gc = new_genetic_code.get_code(gc)
    seq = new_moltype.DNA.make_seq(seq=seq)
    assert not seq.has_terminal_stop(gc=gc)


def test_has_terminal_stop_strict():
    gc = new_genetic_code.get_code(1)
    seq = new_moltype.DNA.make_seq(seq="TCCAG")
    with pytest.raises(new_alphabet.AlphabetError):
        seq.has_terminal_stop(gc=gc, strict=True)


@pytest.mark.parametrize(
    "gc,seq",
    (
        (2, "TCCAGG"),
        (1, "TAATGA"),
        (1, "ACGTGA---"),
        (1, "--AT-CTGA"),
    ),
)
def test_trim_terminal_stop_true(gc, seq):
    gc = new_genetic_code.get_code(gc)
    expect = re.sub("(TGA|AGG)(?=[-]*$)", "---" if "-" in seq else "", seq)

    seq = new_moltype.DNA.make_seq(seq=seq)
    got = str(seq.trim_stop_codon(gc=gc))
    assert got == expect


@pytest.mark.parametrize("gc,seq", ((1, "T?CTGC"), (2, "TCCAAG")))
def test_trim_terminal_stop_nostop(gc, seq):
    gc = new_genetic_code.get_code(gc)
    seq = new_moltype.DNA.make_seq(seq=seq)
    got = seq.trim_stop_codon(gc=gc)
    assert str(got) == str(seq)
    # since there's no stop, we just return the same object
    assert got is seq


@pytest.mark.parametrize(
    "gc,seq", ((1, "TCCAGG"), (2, "TCCAAA"), (1, "CCTGA"), (2, "CCAGG"))
)
def test_trim_terminal_stop_false(gc, seq):
    gc = new_genetic_code.get_code(gc)
    seq = new_moltype.DNA.make_seq(seq=seq)
    assert str(seq.trim_stop_codon(gc=gc)) == str(seq)


def test_trim_terminal_stop_strict():
    gc = new_genetic_code.get_code(1)
    seq = new_moltype.DNA.make_seq(seq="TCCAG")
    with pytest.raises(new_alphabet.AlphabetError):
        seq.trim_stop_codon(gc=gc, strict=True)


@pytest.mark.parametrize("cast", (int, numpy.int32, numpy.int64, numpy.uint8))
def test_index_a_seq(cast):
    seq = new_moltype.DNA.make_seq(seq="TCCAG")
    got = seq[cast(1)]
    assert isinstance(got, new_sequence.Sequence)


@pytest.mark.parametrize("cast", (float, numpy.float32))
def test_index_a_seq_float_fail(cast):
    seq = new_moltype.DNA.make_seq(seq="TCCAG")
    index = cast(1)
    with pytest.raises(TypeError):
        seq[index]  # pylint: disable=W0104


@pytest.mark.parametrize("moltype", ("dna", "protein"))
def test_same_moltype(moltype):
    moltype = new_moltype.get_moltype(moltype)
    seq = moltype.make_seq(seq="TCCAG")
    got = seq.to_moltype(moltype)
    assert got is seq


def test_gapped_by_map_segment_iter():
    moltype = new_moltype.DNA
    m, seq = moltype.make_seq(seq="-TCC--AG").parse_out_gaps()
    g = list(seq.gapped_by_map_segment_iter(m, allow_gaps=True, recode_gaps=False))
    assert g == ["-", "TCC", "--", "AG"]


@pytest.mark.parametrize("rev", (False, True))
@pytest.mark.parametrize("sliced", (False, True))
@pytest.mark.parametrize("start_stop", ((None, None), (3, 7)))
def test_copied_parent_coordinates(sliced, rev, start_stop):
    orig_name = "orig"
    seq = new_moltype.DNA.make_seq(seq="ACGGTGGGAC", name=orig_name)
    start, stop = start_stop
    start = start or 0
    stop = stop or len(seq)
    sl = slice(start, stop)
    seq = seq[sl]
    sliced_name = "sliced"
    seq.name = sliced_name
    assert seq.name == sliced_name
    seq = seq.rc() if rev else seq
    copied = seq.copy(sliced=sliced)
    assert copied.name == sliced_name
    # matches original
    assert copied.parent_coordinates() == seq.parent_coordinates()
    # and expected -- the coordinate name always reflects the underlying sequence
    assert copied.parent_coordinates() == (orig_name, start, stop, -1 if rev else 1)


@pytest.mark.parametrize("rev", (False, True))
def test_parent_coordinates(rev):
    seq = new_moltype.DNA.make_seq(seq="ACGGTGGGAC")
    seq = seq[1:1]
    seq = seq.rc() if rev else seq
    seq.name = "sliced"  # this assignment does not affect the
    # note that when a sequence has zero length, the parent seqid is None
    assert seq.parent_coordinates() == (None, 0, 0, 1)


@pytest.mark.parametrize("cls", (str, bytes))
def test_coerce_to_seqview_str_bytes(cls, dna_alphabet):
    seq = "AC--GGTGGGAC"
    seqid = "seq1"
    s = bytes(seq, "utf8") if cls == bytes else seq
    got = new_sequence._coerce_to_seqview(s, seqid, alphabet=dna_alphabet)
    assert got.str_value == seq
    assert isinstance(got, new_sequence.SeqView)


def test_coerce_to_seqview_sequence(dna_alphabet):
    seq = "AC--GGTGGGAC"
    seqid = "seq1"
    got = new_sequence._coerce_to_seqview(
        new_moltype.DNA.make_seq(seq=seq), seqid, alphabet=dna_alphabet
    )
    assert got.str_value == seq
    assert isinstance(got, new_sequence.SeqView)


def test_coerce_to_seqview_already_seqview(dna_alphabet):
    seq = "AC--GGTGGGAC"
    seqid = "seq1"
    got = new_sequence._coerce_to_seqview(
        new_sequence.SeqView(seq=seq, alphabet=dna_alphabet),
        seqid,
        alphabet=dna_alphabet,
    )
    assert got.str_value == seq
    assert isinstance(got, new_sequence.SeqView)


def test_seqview_seqid(dna_alphabet):
    sv = new_sequence.SeqView(seq="ACGGTGGGAC", alphabet=dna_alphabet)
    assert sv.seqid is None

    sv = new_sequence.SeqView(seq="ACGGTGGGAC", seqid="seq1", alphabet=dna_alphabet)
    assert sv.seqid == "seq1"


def test_seqview_rich_dict_round_trip_seqid(dna_alphabet):
    sv = new_sequence.SeqView(seq="ACGGTGGGAC", seqid="seq1", alphabet=dna_alphabet)
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seqid"] == "seq1"

    got = new_sequence.SeqView.from_rich_dict(rd)
    assert got.seqid == "seq1"

    sv = new_sequence.SeqView(seq="ACGGTGGGAC", alphabet=dna_alphabet)
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seqid"] is None

    got = new_sequence.SeqView.from_rich_dict(rd)
    assert got.seqid is None


def test_seqview_slice_propagates_seqid(dna_alphabet):
    sv = new_sequence.SeqView(seq="ACGGTGGGAC", seqid="seq1", alphabet=dna_alphabet)
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


@pytest.mark.xfail(reason="no SeqView dispatch for new_alphabet.to_indices")
def test_sequences_propogates_seqid_seqview():
    # creating a Sequence with a seqview does not change the seqid of the SeqView.
    seq = new_moltype.DNA.make_seq(
        seq=new_sequence.SeqView(seq="ACGGTGGGAC", seqid="parent_name"), name="seq_name"
    )
    assert seq.name == "seq_name"
    assert seq._seq.seqid == "parent_name"

    # creating a Sequence with an unnamed seqview does not name the SeqView.
    seq = new_moltype.DNA.make_seq(
        new_sequence.SeqView(seq="ACGGTGGGAC"), name="seq_name"
    )
    assert seq.name == "seq_name"
    assert seq._seq.seqid is None


def test_make_seq_assigns_to_seqview():
    seq = new_moltype.DNA.make_seq(seq="ACGT", name="s1")
    assert seq.name == seq._seq.seqid == "s1"


def test_empty_seqview_translate_position(dna_alphabet):
    sv = new_sequence.SeqView(seq="", alphabet=dna_alphabet)
    assert sv.absolute_position(0) == 0
    assert sv.relative_position(0) == 0


@pytest.mark.parametrize("start", (None, 0, 1, 10, -1, -10))
@pytest.mark.parametrize("stop", (None, 10, 8, 1, 0, -1, -11))
@pytest.mark.parametrize("step", (None, 1, 2, -1, -2))
@pytest.mark.parametrize("length", (1, 8, 999))
def test_seqview_seq_len_init(start, stop, step, length, dna_alphabet):
    # seq_len is length of seq when None
    seq_data = "A" * length
    sv = new_sequence.SeqView(
        seq=seq_data, start=start, stop=stop, step=step, alphabet=dna_alphabet
    )
    expect = len(seq_data)
    # Check property and slot
    assert sv.seq_len == expect
    assert sv._seq_len == expect


@pytest.mark.parametrize("seq, seq_len", [("A", 0), ("", 1), ("A", 2)])
def test_seqview_seq_len_mismatch(seq, seq_len, dna_alphabet):
    # If provided, seq_len must match len(seq)
    with pytest.raises(AssertionError):
        new_sequence.SeqView(seq=seq, seq_len=seq_len, alphabet=dna_alphabet)


def test_seqview_copy_propagates_seq_len(dna_alphabet):
    seq = "ACGGTGGGAC"
    sv = new_sequence.SeqView(seq=seq, alphabet=dna_alphabet)
    copied = sv.copy()
    assert copied.seq_len == len(seq)


def test_seqview_seq_len_modified_seq(dna_alphabet):
    seq = "ACGGTGGGAC"
    sv = new_sequence.SeqView(seq=seq, alphabet=dna_alphabet)

    sv.seq = "ATGC"  # this should not modify seq_len
    assert sv.seq_len == len(seq)


def test_sequence_str_bytes_array():
    data = "ACGGTGGGAC"
    seq = new_moltype.DNA.make_seq(seq=data)
    print(type(seq))
    assert str(seq) == data
    assert bytes(seq) == data.encode("utf8")
    assert numpy.array_equal(
        numpy.array(seq), new_moltype.DNA.alphabet.to_indices(data)
    )
