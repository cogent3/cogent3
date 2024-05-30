import numpy as numpy
import pytest

from cogent3 import get_moltype
from cogent3._version import __version__
from cogent3.core import new_alignment as new_aln
from cogent3.core.alphabet import CharAlphabet


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
    return new_aln.SeqsData(data=simple_dict, alphabet=dna_alpha)


@pytest.fixture
def int_arr():
    return numpy.arange(17, dtype=numpy.uint8)


@pytest.fixture
def sdv_s2(sd_demo: new_aln.SeqsData) -> new_aln.SeqDataView:
    return sd_demo.get_seq_view(seqid="seq2")


@pytest.fixture
def simple_dna_seq_coll(sd_demo: new_aln.SeqsData) -> new_aln.SequenceCollection:
    moltype = get_moltype("dna")
    return new_aln.SequenceCollection(seq_data=sd_demo, moltype=moltype)


@pytest.fixture
def aligned_dict():
    return dict(seq1="ACG--T", seq2="-CGAAT")


@pytest.fixture
def ad_demo(aligned_dict: dict[str, str]):
    return new_aln.AlignedData.from_gapped_seqs(aligned_dict)


def test_seqsdata_default_attributes(sd_demo: new_aln.SeqsData):
    assert sd_demo.names == ["seq1", "seq2"]
    assert isinstance(sd_demo._alphabet, CharAlphabet)


@pytest.mark.xfail(reason="SeqsData currently expects correctly formatted data")
def test_seqsdata_seq_if_str(seq1: str, dna_alpha):
    with pytest.raises(NotImplementedError):
        new_aln.SeqsData(seq1, alphabet=dna_alpha)


@pytest.mark.xfail(reason="not currently supporting setting names")
@pytest.mark.parametrize(
    "bad_names", [("bad"), ("bad",), ("bad2", "bad1"), ("seq1",), "seq1"]
)
def test_names_init(simple_dict, dna_alpha, bad_names):
    sd = new_aln.SeqsData(simple_dict, alphabet=dna_alpha)
    assert sd.names == ["seq1", "seq2"]
    sd = new_aln.SeqsData(simple_dict, alphabet=dna_alpha, names=("seq2", "seq1"))
    assert sd.names == ["seq2", "seq1"]

    with pytest.raises(ValueError):
        new_aln.SeqsData(simple_dict, alphabet=dna_alpha, names=bad_names)


def test_seqdataview_zero_step_raises(sd_demo):
    with pytest.raises(ValueError):
        new_aln.SeqDataView(seq=sd_demo, step=0)


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
def test_seqdataview_repr_default(sd_demo: new_aln.SeqsData, seqid: str):
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
    sd = new_aln.SeqsData(d, alphabet=dna_alpha)
    got = sd.get_seq_view(seqid="long")
    assert repr(got) == expect


@pytest.mark.xfail(reason="do we want to support copying and/or copying with slicing?")
@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
@pytest.mark.parametrize("sliced", (False, True))
@pytest.mark.parametrize("step", (1, -1))
def test_seqdataview_copy(dna_alpha, seqid, sliced, step):
    seq1 = "ATGTTCTC"
    seq2 = "ATGAACTCATT"
    start, stop = 2, 6
    data = {"seq1": seq1, "seq2": seq2}

    sd = new_aln.SeqsData(data=data, alphabet=dna_alpha)
    sdv = sd.get_seq_view(seqid=seqid)
    sliced_sdv = sdv[start:stop:step]

    copied_sdv = sliced_sdv.copy(sliced=sliced)

    assert copied_sdv.str_value == data[seqid][start:stop:step]


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
def test_seqsdata_get_seq_view(simple_dict, dna_alpha, seqid):
    sd = new_aln.SeqsData(simple_dict, alphabet=dna_alpha)
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
    sd = new_aln.SeqsData(data=simple_dict, alphabet=dna_alpha)
    got = sd.get_seq_str(seqid=seq, start=start, stop=stop)
    assert expect == got


def test_get_seq_str_empty(sd_demo: new_aln.SeqsData):
    with pytest.raises(TypeError):
        sd_demo.get_seq_str()


def test_names_get(simple_dict, dna_alpha):
    expect = simple_dict.keys()
    sd = new_aln.SeqsData(simple_dict, alphabet=dna_alpha)
    got = sd.names
    # returns iterator
    assert list(got) == list(expect)


@pytest.mark.parametrize(
    "seq, moltype_name",
    [("TCAG-NRYWSKMBDHV?", "dna"), ("UCAG-NRYWSKMBDHV?", "rna")],
)
def test_seq_index_str(seq, moltype_name, int_arr):
    alpha = get_moltype(moltype_name).alphabets.degen_gapped
    got = new_aln.seq_index(seq, alpha)
    assert numpy.array_equal(got, int_arr)


@pytest.mark.parametrize(
    "seq, moltype_name",
    [(b"TCAG-NRYWSKMBDHV?", "dna"), (b"UCAG-NRYWSKMBDHV?", "rna")],
)
def test_seq_index_bytes(seq, moltype_name, int_arr):
    alpha = get_moltype(moltype_name).alphabets.degen_gapped
    got = new_aln.seq_index(seq, alpha)
    assert numpy.array_equal(got, int_arr)


@pytest.mark.parametrize("moltype_name", ("dna", "rna", "protein", "protein_with_stop"))
def test_seq_index_arr(moltype_name, int_arr):
    alpha = get_moltype(moltype_name).alphabets.degen_gapped
    got = new_aln.seq_index(int_arr, alpha)
    assert numpy.array_equal(got, int_arr)
    assert got.dtype == int_arr.dtype


@pytest.mark.parametrize("moltype_name", ("dna", "rna", "protein"))
def test_seq_index_raises(moltype_name):
    bad_seq = None
    alpha = get_moltype(moltype_name).alphabets.degen_gapped

    with pytest.raises(NotImplementedError):
        new_aln.seq_index(bad_seq, alpha)


def test_get_seq_array(simple_dict, dna_alpha):
    # TODO: slicing should be tested here, not get_seq_str
    # todo: kath, edit and parametrize for all moltypes
    # seq1
    expect = numpy.array([2, 1, 3, 0], dtype="uint8")
    sd = new_aln.SeqsData(data=simple_dict, alphabet=dna_alpha)
    got = sd.get_seq_array(seqid="seq1")
    assert numpy.array_equal(got, expect)

    # seq2
    expect = numpy.array([3, 0, 0, 0, 3, 1, 2], dtype="uint8")
    sd = new_aln.SeqsData(data=simple_dict, alphabet=dna_alpha)
    got = sd.get_seq_array(seqid="seq2")
    assert numpy.array_equal(got, expect)


def test_get_seq_bytes(sd_demo: new_aln.SeqsData):
    # getting seqid and slicing tested in test_get_seq_str
    got = sd_demo.get_seq_bytes(seqid="seq1")
    assert isinstance(got, bytes)


def test_seqdataview_make_seq_default(sd_demo):
    assert sd_demo.make_seq is None


@pytest.mark.parametrize(
    "maker", (get_moltype("dna").make_seq, get_moltype("dna").make_array_seq)
)
def test_seqdataview_make_seq_setget(sd_demo, maker):
    sd_demo.make_seq = maker
    assert sd_demo.make_seq is maker


@pytest.mark.parametrize("seq", ("seq1", "seq2"))
def test_getitem_str_1(sd_demo, seq):
    got = sd_demo[seq]
    assert got.seq == sd_demo
    assert got.seqid == seq


@pytest.mark.xfail(reason="'SeqDataView' object has no attribute 'replace'")
@pytest.mark.parametrize("idx", (0, 1))
@pytest.mark.parametrize("make_seq", (True, False))
def test_getitem_int(simple_dict, dna_alpha, idx, make_seq):
    sd = new_aln.SeqsData(simple_dict, alphabet=dna_alpha)
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
def test_seqdataview_returns_self(sd_demo: new_aln.SeqsData):
    sdv = sd_demo.get_seq_view("seq1")
    assert isinstance(sdv, new_aln.SeqDataView)


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
    sdv = new_aln.SeqDataView(seq=seq1, seqid="seq1", seq_len=len(seq1))
    got = sdv[index]
    assert isinstance(got, new_aln.SeqDataView)


# SeqDataView tests for value properties
@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seqdataview_value(simple_dict: dict, dna_alpha, start, stop, step):
    seq = "seq2"
    expect = simple_dict[seq][start:stop:step]
    sd = new_aln.SeqsData(data=simple_dict, alphabet=dna_alpha)
    # Get SeqDataView on seq
    sdv = sd.get_seq_view(seqid=seq)
    sdv2 = sdv[start:stop:step]
    got = sdv2.str_value
    assert got == expect


@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seqsdata_array_value(simple_dict_arr: dict, dna_alpha, start, stop, step):
    seq = "seq2"
    expect = simple_dict_arr[seq][start:stop:step]
    sd = new_aln.SeqsData(data=simple_dict_arr, alphabet=dna_alpha)
    # Get SeqDataView on seq
    sdv = sd.get_seq_view(seqid=seq)
    got = sdv.array_value[start:stop:step]
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seqsdata_bytes_value(simple_dict: dict, dna_alpha, start, stop, step):
    seq = "seq2"
    expect = simple_dict[seq][start:stop:step]
    expect = expect.encode("utf8")
    sd = new_aln.SeqsData(data=simple_dict, alphabet=dna_alpha)
    # Get SeqDataView on seq
    sdv = sd.get_seq_view(seqid=seq)
    got = sdv.bytes_value[start:stop:step]
    assert expect == got


# SeqDataView tests for special methods that access "value" properties
def test_array(sdv_s2: new_aln.SeqDataView):
    expect = sdv_s2.array_value
    got = numpy.array(sdv_s2)
    assert numpy.array_equal(expect, got)


def test_bytes(sdv_s2: new_aln.SeqDataView):
    expect = sdv_s2.bytes_value
    got = bytes(sdv_s2)
    assert expect == got


# AlignedSeqData tests
def test_from_string_unequal_seqlens():
    data = dict(seq1="A-A", seq2="AAAAAAA--")
    with pytest.raises(ValueError):
        new_aln.AlignedData.from_gapped_seqs(data=data)


@pytest.mark.xfail(reason="AlignedData not yet functional")
def test_aligned_from_string_returns_self(aligned_dict):
    got = new_aln.AlignedData.from_gapped_seqs(data=aligned_dict)
    assert isinstance(got, new_aln.AlignedData)
    # assert gap lengths


@pytest.mark.xfail(reason="AlignedData not yet functional")
@pytest.mark.parametrize("seqid", ("seq1", "seq2"))
def test_get_aligned_view(aligned_dict, seqid):
    ad = new_aln.AlignedData.from_gapped_seqs(aligned_dict)
    got = ad.get_aligned_view(seqid)
    assert isinstance(got, new_aln.AlignedDataView)
    assert got.seq == ad
    assert got.stop == ad.align_len
    assert got.seq_len == ad.align_len


# AlignedData seq to gaps
def test_seq_to_gap_coords_str_all_gaps():
    parent_seq = "-----"
    expect_gaplen = numpy.array([len(parent_seq)])
    got_ungap, got_map = new_aln.seq_to_gap_coords(
        parent_seq, moltype=get_moltype("dna")
    )
    assert got_ungap == ""
    assert got_map.cum_gap_lengths == expect_gaplen


def test_seq_to_gap_coords_str_no_gaps():
    parent_seq = "ACTGC"
    got_ungap, got_empty_arr = new_aln.seq_to_gap_coords(
        parent_seq, moltype=get_moltype("dna")
    )
    assert got_ungap == parent_seq
    assert got_empty_arr.size == 0


def test_seq_to_gap_coords_arr_all_gaps():
    alpha = get_moltype("dna").alphabets.degen_gapped
    parent_seq = new_aln.seq_index("-----", alpha)
    got_ungap, got_map = new_aln.seq_to_gap_coords(
        parent_seq, moltype=get_moltype("dna")
    )
    assert got_ungap.size == 0
    assert got_map.get_gap_coordinates() == [[0, 5]]


def test_seq_to_gap_coords_arr_no_gaps():
    alpha = get_moltype("dna").alphabets.degen_gapped
    parent_seq = new_aln.seq_index("ACTGC", alpha)
    got_ungap, got_empty_arr = new_aln.seq_to_gap_coords(
        parent_seq, moltype=get_moltype("dna")
    )
    assert numpy.array_equal(got_ungap, parent_seq)
    assert got_empty_arr.size == 0


@pytest.fixture
def gap_seqs():
    return [
        ("A---CTG-C", [[1, 3], [4, 1]]),
        ("-GTAC--", [[0, 1], [4, 2]]),
        ("---AGC--TGC--", [[0, 3], [3, 2], [6, 2]]),
    ]


@pytest.mark.parametrize("i", range(3))
def test_seq_to_gap_coords_str(gap_seqs, i):
    seq, gap_coords = gap_seqs[i]
    got_ungapped, got_map = new_aln.seq_to_gap_coords(seq, moltype=get_moltype("dna"))
    assert got_ungapped == seq.replace("-", "")
    assert got_map.get_gap_coordinates() == gap_coords


@pytest.mark.parametrize("i", range(3))
def test_seq_to_gap_coords_arr(gap_seqs, i):
    seq, gap_coords = gap_seqs[i]
    alpha = get_moltype("dna").alphabets.degen_gapped
    seq = new_aln.seq_index(seq, alpha)  # convert to array repr
    got_ungapped, got_map = new_aln.seq_to_gap_coords(seq, moltype=get_moltype("dna"))
    assert numpy.array_equal(got_ungapped, seq[seq != 4])  # gap_char = 4
    assert got_map.get_gap_coordinates() == gap_coords


def test_seq_coll_init(simple_dna_seq_coll):
    assert isinstance(simple_dna_seq_coll.seqs, new_aln.SeqsData)
    assert simple_dna_seq_coll.seqs.make_seq is not None


@pytest.mark.parametrize("moltype", ("dna", "rna", "protein", "protein_with_stop"))
def test_make_unaligned_seqs_dict(moltype):
    """test SequenceCollection constructor utility function"""
    data = {"a": "AGGCCC", "b": "AGAAAA"}
    got = new_aln.make_unaligned_seqs(data=data, moltype=moltype)
    assert isinstance(got, new_aln.SequenceCollection)
    assert got._seq_data.make_seq == get_moltype(moltype).make_seq


def test_make_unaligned_seqs_raises():
    data = "AGTCCTGA"
    with pytest.raises(NotImplementedError):
        new_aln.make_unaligned_seqs(data=data, moltype="dna")


@pytest.mark.xfail(
    reason="todo: kath, AttributeError: 'SeqDataView' object has no attribute 'replace'"
)
def test_iter_seqs_ragged_padded():
    """SequenceCollection.iter_seqs() method should support reordering of seqs"""
    ragged_padded = new_aln.make_unaligned_seqs(
        data={"a": "AAAAAA", "b": "AAA---", "c": "AAAA--"}, moltype="dna"
    )
    seqs = list(ragged_padded.iter_seqs())
    assert seqs == ["AAAAAA", "AAA---", "AAAA--"]
    seqs = list(ragged_padded.iter_seqs(seq_order=["b", "a", "a"]))
    assert seqs == ["AAA---", "AAAAAA", "AAAAAA"]
    assert seqs[1] is seqs[2]
    assert seqs[0], ragged_padded.seqs["b"]


@pytest.mark.xfail(
    reason="todo: kath, AttributeError: 'SeqDataView' object has no attribute 'replace'"
)
def test_iter_seqs_ragged(self):
    """SequenceCollection iter_seqs() method should support reordering of seqs"""
    ragged = new_aln.make_unaligned_seqs(
        data={"a": "AAAAAA", "b": "AAA", "c": "AAAA"}, moltype="dna"
    )
    seqs = list(ragged.iter_seqs())
    assert seqs == ["AAAAAA", "AAA", "AAAA"]
    seqs = list(self.ragged.iter_seqs(seq_order=["b", "a", "a"]))
    assert seqs == ["AAA", "AAAAAA", "AAAAAA"]
    assert seqs[1] is seqs[2]
    assert seqs[0], ragged.seqs["b"]


@pytest.mark.xfail(reason="New sc repr not functional yet")
def test_sequence_collection_repr():
    data = {
        "ENSMUSG00000056468": "GCCAGGGGGAAAAGGGAGAA",
        "ENSMUSG00000039616": "GCCCTTCAAATTT",
    }
    seqs = new_aln.make_unaligned_seqs(data=data, moltype="dna")
    assert (
        repr(seqs)
        == "2x (ENSMUSG00000039616[GCCCTTCAAA...], ENSMUSG00000056468[GCCAGGGGGA...]) dna seqcollection"
    )

    data = {
        "ENSMUSG00000039616": "GCCCTTCAAATTT",
        "ENSMUSG00000056468": "GCCAGGGGGAAAAGGGAGAA",
    }
    seqs = new_aln.make_unaligned_seqs(data=data, moltype="dna")
    assert (
        repr(seqs)
        == "2x (ENSMUSG00000039616[GCCCTTCAAA...], ENSMUSG00000056468[GCCAGGGGGA...]) dna seqcollection"
    )

    data = {
        "a": "TCGAT",
    }
    seqs = new_aln.make_unaligned_seqs(data=data, moltype="dna")
    assert repr(seqs) == "1x (a[TCGAT]) dna seqcollection"

    data = {
        "a": "TCGAT" * 2,
    }
    seqs = new_aln.make_unaligned_seqs(data=data, moltype="dna")
    assert repr(seqs) == "1x (a[TCGATTCGAT]) dna seqcollection"

    data = {
        "a": "A" * 11,
        "b": "B" * 3,
        "c": "C" * 3,
        "d": "D" * 11,
        "e": "E" * 8,
    }
    seqs = new_aln.make_unaligned_seqs(data=data, moltype="ASCII")
    assert repr(seqs) == "5x (b[BBB], ..., d[DDDDDDDDDD...]) text seqcollection"
