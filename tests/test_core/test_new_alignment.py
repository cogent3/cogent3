import numpy as numpy
import pytest

from cogent3 import get_moltype
from cogent3._version import __version__
from cogent3.core.alphabet import CharAlphabet
from cogent3.core.new_alignment import (
    AlignedData,
    AlignedDataView,
    SeqData,
    SeqDataView,
    SequenceCollection,
    _names_list_tuple,
    make_unaligned_seqs,
    seq_index,
    seq_to_gap_coords,
    validate_names,
)


@pytest.fixture
def seq1():
    return "ACTG"


@pytest.fixture
def simple_dict():
    return dict(seq1="ACGT", seq2="GTTTGCA")


@pytest.fixture
def simple_dict_arr():
    return dict(seq1=numpy.array([2, 1, 3, 0]), seq2=numpy.array([3, 0, 0, 0, 3, 1, 2]))


@pytest.fixture
def dna_alpha():
    moltype = get_moltype("dna")
    return moltype.alphabets.degen_gapped


@pytest.fixture
def sd_demo(simple_dict: dict[str, str], dna_alpha):
    return SeqData(data=simple_dict, alphabet=dna_alpha)


@pytest.fixture
def int_arr():
    return numpy.arange(17, dtype=numpy.uint8)


@pytest.fixture
def sdv_s2(sd_demo: SeqData) -> SeqDataView:
    return sd_demo.get_seq_view(seqid="seq2")


@pytest.fixture
def simple_dna_seq_coll(sd_demo: SeqData) -> SequenceCollection:
    moltype = get_moltype("dna")
    return SequenceCollection(seq_data=sd_demo, moltype=moltype)


@pytest.fixture
def aligned_dict():
    return dict(seq1="ACG--T", seq2="-CGAAT")


@pytest.fixture
def ad_demo(aligned_dict: dict[str, str]):
    return AlignedData.from_gapped_seqs(aligned_dict)


def test_seqdata_default_attributes(sd_demo: SeqData):
    assert sd_demo._names == ("seq1", "seq2")
    assert isinstance(sd_demo._alphabet, CharAlphabet)


def test_seqdata_seq_if_str(seq1: str, dna_alpha):
    with pytest.raises(NotImplementedError):
        SeqData(seq1, alphabet=dna_alpha)
    # assert SeqData(data)._names != ("A", "C", "T", "G")


def test_validate_names_dict_happy(simple_dict):
    got = validate_names(simple_dict, None)
    assert got == ("seq1", "seq2")
    got = validate_names(simple_dict, names=("seq2", "seq1"))
    assert got == ("seq2", "seq1")


@pytest.mark.parametrize(
    "bad_names", [("bad"), ("bad",), ("bad2", "bad1"), ("seq1",), "seq1"]
)
def test_validate_names_dict_bad(simple_dict, bad_names):
    with pytest.raises(ValueError):
        validate_names(simple_dict, bad_names)


@pytest.mark.parametrize(
    "names", (["seq2"], ["seq2", "seq1"], ("seq2",), ("seq2", "seq1"))
)
def test_names_tuple_list_happy(names, sd_demo):
    correct_names = sd_demo._names
    got = validate_names(correct_names, names)
    assert got == tuple(names)


@pytest.mark.parametrize("bad_names", [("bad"), ("bad",), ("bad2", "bad1"), "seq1"])
def test_names_tuple_list_bad(bad_names, sd_demo):
    correct_names = sd_demo._names
    with pytest.raises(ValueError):
        validate_names(correct_names, bad_names)


def test_names_tuple_list_None(sd_demo):
    correct_names = sd_demo._names
    got = validate_names(correct_names, None)
    assert got == correct_names


@pytest.mark.parametrize(
    "bad_names", [("bad"), ("bad",), ("bad2", "bad1"), ("seq1",), "seq1"]
)
def test_names_init(simple_dict, dna_alpha, bad_names):
    sd = SeqData(simple_dict, alphabet=dna_alpha)
    assert sd._names == ("seq1", "seq2")
    sd = SeqData(simple_dict, alphabet=dna_alpha, names=("seq2", "seq1"))
    assert sd._names == ("seq2", "seq1")

    with pytest.raises(ValueError):
        SeqData(simple_dict, alphabet=dna_alpha, names=bad_names)


def test_seqdataview_zero_step_raises(sd_demo):
    with pytest.raises(ValueError):
        SeqDataView(seq=sd_demo, step=0)


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
def test_seqdataview_repr_default(sd_demo: SeqData, seqid: str):
    seq = sd_demo.get_seq_str(seqid=seqid)
    seq_len = len(seq)
    expect = f"SeqDataView(seq={seq}, start=0, stop={seq_len}, step=1, offset=0, seqid='{seqid}', seq_len={seq_len})"
    got = sd_demo.get_seq_view(seqid)
    assert repr(got) == expect


def test_seqdataview_repr_default_long(dna_alpha):
    longseq = "CGATCGTAGTACGTGTCAAGTCTGAC"
    seq_len = len(longseq)
    trunc = f"{longseq[:10]}...{longseq[-5:]}"
    expect = f"SeqDataView(seq={trunc}, start=0, stop={seq_len}, step=1, offset=0, seqid='long', seq_len={seq_len})"

    d = {"long": longseq}
    sd = SeqData(d, alphabet=dna_alpha)
    got = sd.get_seq_view(seqid="long")
    assert repr(got) == expect


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
@pytest.mark.parametrize("step", [1, -1])
def test_seqdataview_to_rich_dict(sd_demo: SeqData, seqid, step):
    seq = sd_demo.get_seq_str(seqid=seqid)
    sdv = sd_demo.get_seq_view(seqid)
    sdv = sdv[::step]
    expect = {
        "type": "cogent3.core.new_alignment.SeqDataView",
        "version": __version__,
        "init_args": {
            "seq": seq[::step],
            "seqid": seqid,
            "seq_len": sdv.seq_len,
            "step": step,
            "offset": 0,
        },
    }

    got = sdv.to_rich_dict()
    assert expect == got


def test_seqdataview_rich_dict_round_trip(sdv_s2):
    sdv_rd = sdv_s2.to_rich_dict()
    # creating rd twice because from_rich_dict modifies data in place
    sdv_rd2 = sdv_s2.to_rich_dict()
    got = SeqDataView.from_rich_dict(sdv_rd2)

    for init_arg, init_arg_val in sdv_rd["init_args"].items():
        assert getattr(got, init_arg) == init_arg_val


@pytest.mark.xfail(
    reason="copied SDV.seq is not a SD instance, but a string when sliced==True"
)
@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
@pytest.mark.parametrize("sliced", (False, True))
@pytest.mark.parametrize("step", (1, -1))
def test_seqdataview_copy(dna_alpha, seqid, sliced, step):
    seq1 = "ATGTTCTC"
    seq2 = "ATGAACTCATT"
    start, stop = 2, 6
    data = {"seq1": seq1, "seq2": seq2}

    sd = SeqData(data=data, alphabet=dna_alpha)
    sdv = sd.get_seq_view(seqid=seqid)
    sliced_sdv = sdv[start:stop:step]

    copied_sdv = sliced_sdv.copy(sliced=sliced)

    assert copied_sdv.str_value == data[seqid][start:stop:step]


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
def test_seqdata_get_seq_view(simple_dict, dna_alpha, seqid):
    sd = SeqData(simple_dict, alphabet=dna_alpha)
    seq = simple_dict[seqid]
    seq_len = len(seq)
    got = sd.get_seq_view(seqid)
    assert got.seq == sd
    assert got.stop == seq_len
    assert got.seqid == seqid
    assert got.seq_len == seq_len


@pytest.mark.parametrize("seq", ("seq1", "seq2"))
@pytest.mark.parametrize("start", (None, -1, 0, 1, 4))
@pytest.mark.parametrize("stop", (None, -1, 0, 1, 4))
def test_get_seq_str(simple_dict, dna_alpha, seq, start, stop):
    # slicing should be tested in test_get_seq_array
    expect = simple_dict[seq][start:stop]
    sd = SeqData(data=simple_dict, alphabet=dna_alpha)
    got = sd.get_seq_str(seqid=seq, start=start, stop=stop)
    assert expect == got


def test_get_seq_str_empty(sd_demo: SeqData):
    with pytest.raises(TypeError):
        sd_demo.get_seq_str()


def test_names_get(simple_dict, dna_alpha):
    expect = simple_dict.keys()
    sd = SeqData(simple_dict, alphabet=dna_alpha)
    got = sd.names
    # returns iterator
    assert list(got) == list(expect)


@pytest.mark.parametrize(
    "seq, moltype_name",
    [("TCAG-NRYWSKMBDHV?", "dna"), ("UCAG-NRYWSKMBDHV?", "rna")],
)
def test_seq_index_str(seq, moltype_name, int_arr):
    alpha = get_moltype(moltype_name).alphabets.degen_gapped
    got = seq_index(seq, alpha)
    assert numpy.array_equal(got, int_arr)


@pytest.mark.parametrize(
    "seq, moltype_name",
    [(b"TCAG-NRYWSKMBDHV?", "dna"), (b"UCAG-NRYWSKMBDHV?", "rna")],
)
def test_seq_index_bytes(seq, moltype_name, int_arr):
    alpha = get_moltype(moltype_name).alphabets.degen_gapped
    got = seq_index(seq, alpha)
    assert numpy.array_equal(got, int_arr)


@pytest.mark.parametrize("moltype_name", ("dna", "rna", "protein", "protein_with_stop"))
def test_seq_index_arr(moltype_name, int_arr):
    alpha = get_moltype(moltype_name).alphabets.degen_gapped
    got = seq_index(int_arr, alpha)
    assert numpy.array_equal(got, int_arr)
    assert got.dtype == int_arr.dtype


@pytest.mark.parametrize("moltype_name", ("dna", "rna", "protein"))
def test_seq_index_raises(moltype_name):
    bad_seq = None
    alpha = get_moltype(moltype_name).alphabets.degen_gapped

    with pytest.raises(NotImplementedError):
        seq_index(bad_seq, alpha)


def test_get_seq_array(simple_dict, dna_alpha):
    # TODO: slicing should be tested here, not get_seq_str
    # todo: kath, edit and parametrize for all moltypes
    # seq1
    expect = numpy.array([2, 1, 3, 0], dtype="uint8")
    sd = SeqData(data=simple_dict, alphabet=dna_alpha)
    got = sd.get_seq_array(seqid="seq1")
    assert numpy.array_equal(got, expect)

    # seq2
    expect = numpy.array([3, 0, 0, 0, 3, 1, 2], dtype="uint8")
    sd = SeqData(data=simple_dict, alphabet=dna_alpha)
    got = sd.get_seq_array(seqid="seq2")
    assert numpy.array_equal(got, expect)


def test_get_seq_bytes(sd_demo: SeqData):
    # getting seqid and slicing tested in test_get_seq_str
    got = sd_demo.get_seq_bytes(seqid="seq1")
    assert isinstance(got, bytes)


# SeqData set_seq_maker tests
def test_seqdataview_make_seq_default(sd_demo):
    assert sd_demo.make_seq is None


@pytest.mark.parametrize(
    "maker", (get_moltype("dna").make_seq, get_moltype("dna").make_array_seq)
)
def test_seqdataview_make_seq_setget(sd_demo, maker):
    sd_demo.make_seq = maker
    assert sd_demo.make_seq is maker


# SeqData __getitem__ tests
@pytest.mark.parametrize("seq", ("seq1", "seq2"))
def test_getitem_str_1(sd_demo, seq):
    got = sd_demo[seq]
    assert got.seq == sd_demo
    assert got.seqid == seq


@pytest.mark.xfail(reason="'SeqDataView' object has no attribute 'replace'")
@pytest.mark.parametrize("idx", (0, 1))
@pytest.mark.parametrize("make_seq", (True, False))
def test_getitem_int(simple_dict, dna_alpha, idx, make_seq):
    sd = SeqData(simple_dict, alphabet=dna_alpha)
    ms = get_moltype("dna").make_seq if make_seq else None
    sd.make_seq = ms
    got = sd[idx]
    assert got.seq == sd
    assert got.seqid == list(simple_dict)[idx]


@pytest.mark.xfail(reason="'SeqDataView' object has no attribute 'replace'")
@pytest.mark.parametrize("seqid", ("seq1", "seq2"))
@pytest.mark.parametrize("make_seq", (True, False))
def test_getitem_str(sd_demo, seqid, make_seq):
    ms = get_moltype("dna").make_seq if make_seq else None
    sd_demo.make_seq = ms
    got = sd_demo[seqid]
    assert got.seq == sd_demo
    assert got.seqid == seqid


@pytest.mark.parametrize("make_seq", (True, False))
def test_getitem_raises(sd_demo, make_seq):
    ms = get_moltype("dna").make_seq if make_seq else None
    sd_demo.make_seq = ms
    invalid_index = ["this", "shouldn't", "work"]
    with pytest.raises(NotImplementedError):
        sd_demo[invalid_index]


# SeqDataView tests for returning an instance of itself
def test_seqdataview_returns_self(sd_demo: SeqData):
    sdv = sd_demo.get_seq_view("seq1")
    assert isinstance(sdv, SeqDataView)


@pytest.mark.parametrize(
    "index",
    [
        slice(1, 3, None),
        slice(None, 2, None),
        slice(3, None, None),
        slice(0, 10, None),
        slice(None, None, 2),
        slice(None),
    ],
)
def test_seqdataview_slice_returns_self(seq1: str, index: slice):
    sdv = SeqDataView(seq=seq1, seqid="seq1", seq_len=len(seq1))
    got = sdv[index]
    assert isinstance(got, SeqDataView)


# SeqDataView tests for value properties
@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seqdataview_value(simple_dict: dict, dna_alpha, start, stop, step):
    seq = "seq2"
    expect = simple_dict[seq][start:stop:step]
    sd = SeqData(data=simple_dict, alphabet=dna_alpha)
    # Get SeqDataView on seq
    sdv = sd.get_seq_view(seqid=seq)
    sdv2 = sdv[start:stop:step]
    got = sdv2.str_value
    assert got == expect


@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seqdata_array_value(simple_dict_arr: dict, dna_alpha, start, stop, step):
    seq = "seq2"
    expect = simple_dict_arr[seq][start:stop:step]
    sd = SeqData(data=simple_dict_arr, alphabet=dna_alpha)
    # Get SeqDataView on seq
    sdv = sd.get_seq_view(seqid=seq)
    got = sdv.array_value[start:stop:step]
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seqdata_bytes_value(simple_dict: dict, dna_alpha, start, stop, step):
    seq = "seq2"
    expect = simple_dict[seq][start:stop:step]
    expect = expect.encode("utf8")
    sd = SeqData(data=simple_dict, alphabet=dna_alpha)
    # Get SeqDataView on seq
    sdv = sd.get_seq_view(seqid=seq)
    got = sdv.bytes_value[start:stop:step]
    assert expect == got


# SeqDataView tests for special methods that access "value" properties
def test_array(sdv_s2: SeqDataView):
    expect = sdv_s2.array_value
    got = numpy.array(sdv_s2)
    assert numpy.array_equal(expect, got)


def test_bytes(sdv_s2: SeqDataView):
    expect = sdv_s2.bytes_value
    got = bytes(sdv_s2)
    assert expect == got


# AlignedSeqData tests
def test_from_string_unequal_seqlens():
    data = dict(seq1="A-A", seq2="AAAAAAA--")
    with pytest.raises(ValueError):
        AlignedData.from_gapped_seqs(data=data)


def test_aligned_from_string_returns_self(aligned_dict):
    got = AlignedData.from_gapped_seqs(data=aligned_dict)
    assert isinstance(got, AlignedData)
    # assert gap lengths


@pytest.mark.parametrize("seqid", ("seq1", "seq2"))
def test_get_aligned_view(aligned_dict, seqid):
    ad = AlignedData.from_gapped_seqs(aligned_dict)
    got = ad.get_aligned_view(seqid)
    assert isinstance(got, AlignedDataView)
    assert got.seq == ad
    assert got.stop == ad.align_len
    assert got.seq_len == ad.align_len


# AlignedData seq to gaps
def test_seq_to_gap_coords_str_all_gaps():
    parent_seq = "-----"
    expect_gaplen = numpy.array([len(parent_seq)])
    got_ungap, got_map = seq_to_gap_coords(parent_seq, moltype=get_moltype("dna"))
    assert got_ungap == ""
    assert got_map.cum_gap_lengths == expect_gaplen


def test_seq_to_gap_coords_str_no_gaps():
    parent_seq = "ACTGC"
    got_ungap, got_empty_arr = seq_to_gap_coords(parent_seq, moltype=get_moltype("dna"))
    assert got_ungap == parent_seq
    assert got_empty_arr.size == 0


def test_seq_to_gap_coords_arr_all_gaps():
    alpha = get_moltype("dna").alphabets.degen_gapped
    parent_seq = seq_index("-----", alpha)
    got_ungap, got_map = seq_to_gap_coords(parent_seq, moltype=get_moltype("dna"))
    assert got_ungap.size == 0
    assert got_map.get_gap_coordinates() == [[0, 5]]


def test_seq_to_gap_coords_arr_no_gaps():
    alpha = get_moltype("dna").alphabets.degen_gapped
    parent_seq = seq_index("ACTGC", alpha)
    got_ungap, got_empty_arr = seq_to_gap_coords(parent_seq, moltype=get_moltype("dna"))
    assert numpy.array_equal(got_ungap, parent_seq)
    assert got_empty_arr.size == 0


@pytest.fixture
def gap_seqs():
    return [
        ("A---CTG-C", [[1, 3], [4, 1]]),
        ("-GTAC--", [[0, 1], [4, 2]]),
        ("---AGC--TGC--", [[0, 3], [3, 2], [6, 2]]),
    ]


@pytest.mark.parametrize("i", range(3))  # range(len(gap_seqs()))
def test_seq_to_gap_coords_str(gap_seqs, i):
    seq, gap_coords = gap_seqs[i]
    got_ungapped, got_map = seq_to_gap_coords(seq, moltype=get_moltype("dna"))
    assert got_ungapped == seq.replace("-", "")
    assert got_map.get_gap_coordinates() == gap_coords


@pytest.mark.parametrize("i", range(3))  # range(len(gap_seqs()))
def test_seq_to_gap_coords_arr(gap_seqs, i):
    seq, gap_coords = gap_seqs[i]
    alpha = get_moltype("dna").alphabets.degen_gapped
    seq = seq_index(seq, alpha)  # convert to array repr
    got_ungapped, got_map = seq_to_gap_coords(seq, moltype=get_moltype("dna"))
    assert numpy.array_equal(got_ungapped, seq[seq != 4])  # gap_char = 4
    assert got_map.get_gap_coordinates() == gap_coords


def test_seq_coll_init(simple_dna_seq_coll):
    assert isinstance(simple_dna_seq_coll.seqs, SeqData)
    assert simple_dna_seq_coll.seqs.make_seq is not None


@pytest.mark.parametrize("moltype", ("dna", "rna", "protein", "protein_with_stop"))
def test_make_unaligned_seqs_dict(moltype):
    """test SequenceCollection constructor utility function"""
    data = {"a": "AGGCCC", "b": "AGAAAA"}
    got = make_unaligned_seqs(data=data, moltype=moltype)
    assert isinstance(got, SequenceCollection)
    assert got._seq_data.make_seq == get_moltype(moltype).make_seq


@pytest.mark.parametrize("collection_type", (list, tuple))
def test_make_unaligned_seqs_list(collection_type):
    """test SequenceCollection constructor utility function"""
    data = [collection_type(["a", "AGGTT"]), collection_type(["b", "AG"])]
    got = make_unaligned_seqs(data=data, moltype="dna")
    assert isinstance(got, SequenceCollection)


def test_make_unaligned_seqs_no_names():
    """test SequenceCollection constructor utility function"""
    data = ["ATTGCGG", "ATTCCCGCC"]
    got = make_unaligned_seqs(data=data, moltype="dna")
    assert isinstance(got, SequenceCollection)


def test_make_unaligned_seqs_list_of_arrays():
    data = [numpy.array([0, 1, 2, 3, 2]), numpy.array([0, 1, 2, 3])]
    got = make_unaligned_seqs(data=data, moltype="dna")
    assert isinstance(got, SequenceCollection)


@pytest.mark.parametrize("moltype_name", ("dna", "rna", "protein", "protein_with_stop"))
def test_make_unaligned_seqs_array(moltype_name):
    data = numpy.array([[0, 1, 2, 3, 2], [0, 1, 2, 3, 3]])
    got = make_unaligned_seqs(data=data, moltype=moltype_name)
    assert isinstance(got, SequenceCollection)


def test_make_unaligned_seqs_raises():
    data = "AGTCCTGA"
    with pytest.raises(NotImplementedError):
        make_unaligned_seqs(data=data, moltype="dna")
