import numpy as numpy
import pytest

from cogent3 import Sequence, get_moltype
from cogent3.core.alphabet import CharAlphabet
from cogent3.core.new_alignment import (
    AlignedData,
    AlignedDataView,
    SeqData,
    SeqDataView,
    process_name_order,
    seq_index,
    seq_to_gap_coords,
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
def alpha():
    moltype = get_moltype("dna")
    return moltype.alphabets.degen_gapped


@pytest.fixture
def sd_demo(simple_dict: dict[str, str], alpha):
    return SeqData(data=simple_dict, alphabet=alpha)


@pytest.fixture
def int_arr():
    return numpy.arange(17, dtype=numpy.uint8)


@pytest.fixture
def sdv_s2(sd_demo: SeqData) -> SeqDataView:
    return sd_demo.get_seq_view(seqid="seq2")


@pytest.fixture
def aligned_dict():
    return dict(seq1="ACG--T", seq2="-CGAAT")


@pytest.fixture
def ad_demo(aligned_dict: dict[str, str]):
    return AlignedData.from_gapped_seqs(aligned_dict)


def test_seqdata_default_attributes(sd_demo: SeqData):
    assert sd_demo._name_order == ("seq1", "seq2")
    assert isinstance(sd_demo._alphabet, CharAlphabet)


def test_seqdata_seq_if_str(seq1: str, alpha):
    with pytest.raises(NotImplementedError):
        SeqData(seq1, alphabet=alpha)
    # assert SeqData(data)._name_order != ("A", "C", "T", "G")


def test_process_name_order_dict_happy(simple_dict):
    got = process_name_order(simple_dict, None)
    assert got == ("seq1", "seq2")
    got = process_name_order(simple_dict, name_order=("seq2", "seq1"))
    assert got == ("seq2", "seq1")


@pytest.mark.parametrize(
    "bad_names", [("bad"), ("bad",), ("bad2", "bad1"), ("seq1",), "seq1"]
)
def test_process_name_order_dict_bad(simple_dict, bad_names):
    with pytest.raises(ValueError):
        process_name_order(simple_dict, bad_names)


@pytest.mark.parametrize(
    "names", (["seq2"], ["seq2", "seq1"], ("seq2",), ("seq2", "seq1"))
)
def test_name_order_tuple_list_happy(names, sd_demo):
    correct_names = sd_demo._name_order
    got = process_name_order(correct_names, names)
    assert got == tuple(names)


@pytest.mark.parametrize("bad_names", [("bad"), ("bad",), ("bad2", "bad1"), "seq1"])
def test_name_order_tuple_list_bad(bad_names, sd_demo):
    correct_names = sd_demo._name_order
    with pytest.raises(ValueError):
        process_name_order(correct_names, bad_names)


@pytest.mark.parametrize(
    "bad_names", [("bad"), ("bad",), ("bad2", "bad1"), ("seq1",), "seq1"]
)
def test_name_order_init(simple_dict, alpha, bad_names):
    sd = SeqData(simple_dict, alphabet=alpha)
    assert sd._name_order == ("seq1", "seq2")
    sd = SeqData(simple_dict, alphabet=alpha, name_order=("seq2", "seq1"))
    assert sd._name_order == ("seq2", "seq1")

    with pytest.raises(ValueError):
        SeqData(simple_dict, alphabet=alpha, name_order=bad_names)


def test_seqdata_get_seq_view(sd_demo: SeqData):
    got = sd_demo.get_seq_view("seq1")
    expect = f"SeqDataView(seq={sd_demo}, start=0, stop=4, step=1, offset=0, seqid='seq1', seq_len=4)"
    assert repr(got) == expect

    got = sd_demo.get_seq_view("seq2")
    expect = f"SeqDataView(seq={sd_demo}, start=0, stop=7, step=1, offset=0, seqid='seq2', seq_len=7)"
    assert repr(got) == expect


@pytest.mark.parametrize("seq", ("seq1", "seq2"))
@pytest.mark.parametrize("start", (None, -1, 0, 1, 4))
@pytest.mark.parametrize("stop", (None, -1, 0, 1, 4))
def test_get_seq_str(simple_dict, alpha, seq, start, stop):
    # slicing should be tested in test_get_seq_array
    expect = simple_dict[seq][start:stop]
    sd = SeqData(data=simple_dict, alphabet=alpha)
    got = sd.get_seq_str(seqid=seq, start=start, stop=stop)
    assert expect == got


def test_get_seq_str_empty(sd_demo: SeqData):
    with pytest.raises(TypeError):
        sd_demo.get_seq_str()


def test_iter_names(sd_demo: SeqData):
    got = sd_demo.iter_names()
    assert tuple(got) == ("seq1", "seq2")

    got = sd_demo.iter_names(name_order=["seq2", "seq1"])
    assert tuple(got) == ("seq2", "seq1")

    got = sd_demo.iter_names(name_order=["seq2"])
    assert tuple(got) == ("seq2",)

    got = sd_demo.iter_names(name_order=("bad", "bad2"))
    with pytest.raises(ValueError):
        tuple(got)


@pytest.mark.parametrize("idx, seq", ([0, "seq1"], [1, "seq2"]))
def test_iter_seq_view(sd_demo: SeqData, idx, seq):
    got = tuple(sd_demo.iter_seq_view())
    assert len(got) == 2
    assert got[idx].seq == sd_demo
    assert got[idx].seqid == seq


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


@pytest.mark.parametrize("moltype_name", ("dna", "rna"))
def test_seq_index_arr(moltype_name, int_arr):
    alpha = get_moltype(moltype_name).alphabets.degen_gapped
    got = seq_index(int_arr, alpha)
    assert numpy.array_equal(got, int_arr)
    assert got.dtype == int_arr.dtype


def test_get_seq_array(simple_dict, alpha):
    # TODO: slicing should be tested here, not get_seq_str
    # seq1
    expect = numpy.array([2, 1, 3, 0], dtype="uint8")
    sd = SeqData(data=simple_dict, alphabet=alpha)
    got = sd.get_seq_array(seqid="seq1")
    assert numpy.array_equal(got, expect)

    # seq2
    expect = numpy.array([3, 0, 0, 0, 3, 1, 2], dtype="uint8")
    sd = SeqData(data=simple_dict, alphabet=alpha)
    got = sd.get_seq_array(seqid="seq2")
    assert numpy.array_equal(got, expect)


def test_get_seq_bytes(sd_demo: SeqData):
    # getting seqid and slicing tested in test_get_seq_str
    got = sd_demo.get_seq_bytes(seqid="seq1")
    assert isinstance(got, bytes)


# SeqData set_seq_maker tests
def test_set_set_seq_maker_default(sd_demo):
    assert sd_demo.set_seq_maker is None


@pytest.mark.parametrize(
    "maker", (get_moltype("dna").make_seq, get_moltype("dna").make_array_seq)
)
def test_set_seq_maker_setget(sd_demo, maker):
    alpha = sd_demo._alphabet
    sd_demo.set_seq_maker = (maker, alpha)
    assert sd_demo.set_seq_maker == maker


def test_set_seq_maker_wrong_alpha(sd_demo):
    mt = get_moltype("dna")
    maker = mt.make_seq
    alpha = mt.alphabet
    with pytest.raises(ValueError):
        sd_demo.set_seq_maker = maker, alpha


# SeqData __getitem__ tests
@pytest.mark.parametrize("seq", ("seq1", "seq2"))
def test_getitem_str(sd_demo, seq):
    got = sd_demo[seq]
    assert got.seq == sd_demo
    assert got.seqid == seq


@pytest.mark.parametrize("idx", (0, 1))
def test_getitem_int(simple_dict, alpha, idx):
    sd = SeqData(simple_dict, alphabet=alpha)
    got = sd[idx]
    assert got.seq == sd
    assert got.seqid == list(simple_dict)[idx]


# @pytest.mark.parametrize("seqid", ("seq1", "seq2"))
# def test_getitem_str_when_set_seq_maker(sd_demo, seqid):
#     mt = get_moltype("dna")
#     sd_demo.set_seq_maker = mt.make_seq, mt.alphabets.degen_gapped
#     got = sd_demo[seqid] # get item
#     assert got.seq == sd_demo
#     assert got.name == seqid


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
def test_seqdataview_value(simple_dict: dict, alpha, start, stop, step):
    seq = "seq2"
    expect = simple_dict[seq][start:stop:step]
    sd = SeqData(data=simple_dict, alphabet=alpha)
    # Get SeqDataView on seq
    sdv = sd.get_seq_view(seqid=seq)
    sdv2 = sdv[start:stop:step]
    got = sdv2.value
    assert got == expect


@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_array_value(simple_dict_arr: dict, alpha, start, stop, step):
    seq = "seq2"
    expect = simple_dict_arr[seq][start:stop:step]
    sd = SeqData(data=simple_dict_arr, alphabet=alpha)
    # Get SeqDataView on seq
    sdv = sd.get_seq_view(seqid=seq)
    got = sdv.array_value[start:stop:step]
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_bytes_value(simple_dict: dict, alpha, start, stop, step):
    seq = "seq2"
    expect = simple_dict[seq][start:stop:step]
    expect = expect.encode("utf8")
    sd = SeqData(data=simple_dict, alphabet=alpha)
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
