import os

import numpy as numpy
import pytest

import cogent3.core.new_alignment as new_aln
import cogent3.core.new_alphabet as new_alpha
import cogent3.core.new_moltype as new_moltype
import cogent3.core.new_sequence as new_seq

from cogent3._version import __version__


@pytest.fixture
def seq1():
    return "ACTG"


@pytest.fixture
def str_seqs_dict(moltype="dna"):
    if moltype == "dna":
        return dict(seq1="ACGT", seq2="GTTTGCA", seq3="ACGTACGT")


@pytest.fixture
def arr_seqs_dict():
    return dict(
        seq1=numpy.array([2, 1, 3, 0]),
        seq2=numpy.array([3, 0, 0, 0, 3, 1, 2]),
        seq3=numpy.array([2, 1, 3, 0, 2, 1, 3, 0]),
    )


@pytest.fixture
def alpha(moltype="dna"):
    moltype = new_moltype.get_moltype(moltype)
    return moltype.degen_gapped_alphabet


@pytest.fixture
def dna_sd(str_seqs_dict: dict[str, str], alpha):
    return new_aln.SeqsData(data=str_seqs_dict, alphabet=alpha)


@pytest.fixture
def int_arr():
    return numpy.arange(17, dtype=numpy.uint8)


@pytest.fixture
def sdv_s2(dna_sd: new_aln.SeqsData) -> new_aln.SeqDataView:
    return dna_sd.get_seq_view(seqid="seq2")


@pytest.fixture
def seqs() -> new_aln.SequenceCollection:
    data = {"seq1": "AAAAAA", "seq2": "TTTT", "seq3": "ATTCCCC"}
    return new_aln.make_unaligned_seqs(data=data, moltype="dna")


@pytest.fixture
def ragged_padded():
    return new_aln.make_unaligned_seqs(
        data={"a": "AAAAAA", "b": "AAA---", "c": "AAAA--"}, moltype="dna"
    )


@pytest.fixture
def ragged():
    return new_aln.make_unaligned_seqs(
        data={"a": "AAAAAA", "b": "AAA", "c": "AAAA"}, moltype="dna"
    )


@pytest.fixture
def aligned_dict():
    return dict(seq1="ACG--T", seq2="-CGAAT")


@pytest.fixture
def ad_demo(aligned_dict: dict[str, str]):
    return new_aln.AlignedData.from_gapped_seqs(aligned_dict)


@pytest.fixture
def gap_seqs():
    return [
        ("A---CTG-C", [[1, 3], [4, 1]]),
        ("-GTAC--", [[0, 1], [4, 2]]),
        ("---AGC--TGC--", [[0, 3], [3, 2], [6, 2]]),
    ]


def test_seqs_data_default_attributes(dna_sd: new_aln.SeqsData):
    assert dna_sd.names == ["seq1", "seq2", "seq3"]
    assert isinstance(dna_sd._alphabet, new_alpha.CharAlphabet)


@pytest.mark.xfail(reason="SeqsData currently expects correctly formatted data")
def test_seqs_data_seq_if_str(seq1: str, alpha):
    with pytest.raises(NotImplementedError):
        new_aln.SeqsData(seq1, alphabet=alpha)


def test_seqs_data_view_zero_step_raises(dna_sd):
    with pytest.raises(ValueError):
        new_aln.SeqDataView(seqs=dna_sd, step=0, seq_len=4)


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
def test_seqs_data_view_repr_default(dna_sd: new_aln.SeqsData, seqid: str):
    seq = dna_sd.get_seq_str(seqid=seqid)
    seq_len = len(seq)
    expect = f"SeqDataView(seq={seq}, start=0, stop={seq_len}, step=1, offset=0, seqid='{seqid}', seq_len={seq_len})"
    got = dna_sd.get_seq_view(seqid)
    assert repr(got) == expect


def test_seqs_data_view_repr_default_long(alpha):
    longseq = "CGATCGTAGTACGTGTCAAGTCTGAC"
    seq_len = len(longseq)
    trunc = f"{longseq[:10]}...{longseq[-5:]}"
    expect = f"SeqDataView(seq={trunc}, start=0, stop={seq_len}, step=1, offset=0, seqid='long', seq_len={seq_len})"

    d = {"long": longseq}
    sd = new_aln.SeqsData(d, alphabet=alpha)
    got = sd.get_seq_view(seqid="long")
    assert repr(got) == expect


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
@pytest.mark.parametrize("sliced", (False, True))
@pytest.mark.parametrize("step", (1, -1))
def test_seqs_data_view_copy(alpha, seqid, sliced, step):
    seq1 = "ATGTTCTC"
    seq2 = "ATGAACTCATT"
    start, stop = 2, 6
    data = {"seq1": seq1, "seq2": seq2}

    sd = new_aln.SeqsData(data=data, alphabet=alpha)
    sdv = sd.get_seq_view(seqid=seqid)
    sliced_sdv = sdv[start:stop:step]

    copied_sdv = sliced_sdv.copy(sliced=sliced)

    assert copied_sdv.str_value == data[seqid][start:stop:step]


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
def test_seqs_data_get_seq_view(str_seqs_dict, alpha, seqid):
    sd = new_aln.SeqsData(str_seqs_dict, alphabet=alpha)
    seq = str_seqs_dict[seqid]
    seq_len = len(seq)
    got = sd.get_seq_view(seqid)
    assert got.seqs == sd
    assert got.stop == seq_len
    assert got.seqid == seqid
    assert got.seq_len == seq_len


@pytest.mark.parametrize("seq", ("seq1", "seq2"))
@pytest.mark.parametrize("start", (None, -1, 0, 1, 4))
@pytest.mark.parametrize("stop", (None, -1, 0, 1, 4))
def test_seqs_data_get_seq_str(str_seqs_dict, alpha, seq, start, stop):
    # slicing should be tested in test_get_seq_array
    expect = str_seqs_dict[seq][start:stop]
    sd = new_aln.SeqsData(data=str_seqs_dict, alphabet=alpha)
    got = sd.get_seq_str(seqid=seq, start=start, stop=stop)
    assert expect == got


def test_seqs_data_get_seq_str_empty(dna_sd: new_aln.SeqsData):
    with pytest.raises(TypeError):
        dna_sd.get_seq_str()


def test_seqs_data_names(str_seqs_dict, alpha):
    expect = str_seqs_dict.keys()
    sd = new_aln.SeqsData(str_seqs_dict, alphabet=alpha)
    got = sd.names
    # returns iterator
    assert list(got) == list(expect)


def test_seqs_data_get_seq_array(str_seqs_dict, alpha):
    # todo: kath, edit and parametrize for all moltypes
    # seq1
    expect = numpy.array([2, 1, 3, 0], dtype="uint8")
    sd = new_aln.SeqsData(data=str_seqs_dict, alphabet=alpha)
    got = sd.get_seq_array(seqid="seq1")
    assert numpy.array_equal(got, expect)

    # seq2
    expect = numpy.array([3, 0, 0, 0, 3, 1, 2], dtype="uint8")
    sd = new_aln.SeqsData(data=str_seqs_dict, alphabet=alpha)
    got = sd.get_seq_array(seqid="seq2")
    assert numpy.array_equal(got, expect)


def test_seqs_data_get_seq_bytes(dna_sd: new_aln.SeqsData):
    # getting seqid and slicing tested in test_get_seq_str
    got = dna_sd.get_seq_bytes(seqid="seq1")
    assert isinstance(got, bytes)


def test_seqs_data_make_seq_default(dna_sd):
    assert dna_sd.make_seq is None


def test_seqs_data_make_seq_setget(dna_sd):
    dna_sd.make_seq = new_moltype.get_moltype("dna").make_seq
    assert dna_sd.make_seq == new_moltype.get_moltype("dna").make_seq


@pytest.mark.parametrize("seq", ("seq1", "seq2"))
def test_seqs_data_getitem_str_1(dna_sd, seq):
    got = dna_sd[seq]
    assert got.seqs == dna_sd
    assert got.seqid == seq


@pytest.mark.parametrize("idx", (0, 1))
def test_seqs_data_getitem_int(str_seqs_dict, alpha, idx):
    sd = new_aln.SeqsData(str_seqs_dict, alphabet=alpha)
    got = sd[idx]
    assert got.seqs == sd
    assert got.seqid == list(str_seqs_dict)[idx]


@pytest.mark.parametrize("seqid", ("seq1", "seq2"))
def test_seqs_data_getitem_str(
    dna_sd,
    seqid,
):
    got = dna_sd[seqid]
    assert got.seqs == dna_sd
    assert got.seqid == seqid


@pytest.mark.parametrize("make_seq", (True, False))
def test_seqs_data_getitem_raises(dna_sd, make_seq):
    ms = new_moltype.get_moltype("dna").make_seq if make_seq else None
    dna_sd.make_seq = ms
    invalid_index = ["this", "shouldn't", "work"]
    with pytest.raises(NotImplementedError):
        _ = dna_sd[invalid_index]


def test_seqs_data_get_seq_view(dna_sd: new_aln.SeqsData):
    sdv = dna_sd.get_seq_view("seq1")
    assert isinstance(sdv, new_aln.SeqDataView)


def tests_seqs_data_subset(dna_sd):
    got = dna_sd.subset("seq1")
    assert isinstance(got, new_aln.SeqsData)
    assert got.names == ["seq1"]
    assert got.get_seq_str(seqid="seq1") == dna_sd.get_seq_str(seqid="seq1")

    got = dna_sd.subset(["seq1", "seq2"])
    assert isinstance(got, new_aln.SeqsData)
    assert got.names == ["seq1", "seq2"]
    assert got.get_seq_str(seqid="seq1") == dna_sd.get_seq_str(seqid="seq1")
    assert got.get_seq_str(seqid="seq2") == dna_sd.get_seq_str(seqid="seq2")

    got = dna_sd.subset(["seq2", "seq1"])
    assert isinstance(got, new_aln.SeqsData)
    assert got.names == ["seq2", "seq1"]
    assert got.get_seq_str(seqid="seq1") == dna_sd.get_seq_str(seqid="seq1")
    assert got.get_seq_str(seqid="seq2") == dna_sd.get_seq_str(seqid="seq2")

    got = dna_sd.subset(["seq3"])
    assert isinstance(got, new_aln.SeqsData)
    assert got.names == ["seq3"]
    assert got.get_seq_str(seqid="seq3") == dna_sd.get_seq_str(seqid="seq3")


@pytest.mark.xfail(
    reason="todo: kath, decide how to support equality testing between SeqsData"
)
def test_seqs_data_eq(alpha):
    data = dict(seq1="ACGT", seq2="GTTTGCA", seq3="ACGTACGT")
    sd_1 = new_aln.SeqsData(data=data, alphabet=alpha)
    sd_2 = new_aln.SeqsData(data=data, alphabet=alpha)

    assert sd_1 == sd_2
    assert sd_1 == data


@pytest.mark.xfail(
    reason="todo: kath, decide how to support equality testing between SeqsData"
)
def test_seqs_data_eq_diff_alphabets():
    rna_data = dict(seq1="ACGU", seq2="GUUUGCA", seq3="ACGUACGU")
    dna_data = dict(seq1="ACGT", seq2="GTTTGCA", seq3="ACGTACGT")

    dna_alpha = new_moltype.get_moltype("dna").degen_gapped_alphabet
    dna_sd = new_aln.SeqsData(data=dna_data, alphabet=dna_alpha)

    assert dna_data == dna_sd
    assert dna_sd == rna_data


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
def test_seq_data_view_slice_returns_self(seq1: str, index: slice):
    sdv = new_aln.SeqDataView(seqs=seq1, seqid="seq1", seq_len=len(seq1))
    got = sdv[index]
    assert isinstance(got, new_aln.SeqDataView)


# SeqDataView tests for value properties
@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seq_data_view_value(str_seqs_dict: dict, alpha, start, stop, step):
    seq = "seq2"
    expect = str_seqs_dict[seq][start:stop:step]
    sd = new_aln.SeqsData(data=str_seqs_dict, alphabet=alpha)
    # Get SeqDataView on seq
    sdv = sd.get_seq_view(seqid=seq)
    sdv2 = sdv[start:stop:step]
    got = sdv2.str_value
    assert got == expect


@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seqs_data_array_value(arr_seqs_dict: dict, alpha, start, stop, step):
    seq = "seq2"
    expect = arr_seqs_dict[seq][start:stop:step]
    sd = new_aln.SeqsData(data=arr_seqs_dict, alphabet=alpha)
    # Get SeqDataView on seq
    sdv = sd.get_seq_view(seqid=seq)
    got = sdv.array_value[start:stop:step]
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seqs_data_bytes_value(str_seqs_dict: dict, alpha, start, stop, step):
    seq = "seq2"
    expect = str_seqs_dict[seq][start:stop:step]
    expect = expect.encode("utf8")
    sd = new_aln.SeqsData(data=str_seqs_dict, alphabet=alpha)
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
    assert got.seqs == ad
    assert got.stop == ad.align_len
    assert got.seq_len == ad.align_len


@pytest.mark.xfail(reason="AlignedData not yet functional")
def test_seq_to_gap_coords_str_all_gaps():
    parent_seq = "-----"
    expect_gaplen = numpy.array([len(parent_seq)])
    got_ungap, got_map = new_aln.seq_to_gap_coords(
        parent_seq, moltype=new_moltype.get_moltype("dna")
    )
    assert got_ungap == ""
    assert got_map.cum_gap_lengths == expect_gaplen


@pytest.mark.xfail(reason="AlignedData not yet functional")
def test_seq_to_gap_coords_str_no_gaps():
    parent_seq = "ACTGC"
    got_ungap, got_empty_arr = new_aln.seq_to_gap_coords(
        parent_seq, moltype=new_moltype.get_moltype("dna")
    )
    assert got_ungap == parent_seq
    assert got_empty_arr.size == 0


@pytest.mark.xfail(reason="AlignedData not yet functional")
def test_seq_to_gap_coords_arr_all_gaps():
    alpha = new_moltype.get_moltype("dna").degen_gapped_alphabet
    parent_seq = alpha.to_indices("-----", alpha)
    got_ungap, got_map = new_aln.seq_to_gap_coords(
        parent_seq, moltype=new_moltype.get_moltype("dna")
    )
    assert got_ungap.size == 0
    assert got_map.get_gap_coordinates() == [[0, 5]]


@pytest.mark.xfail(reason="AlignedData not yet functional")
def test_seq_to_gap_coords_arr_no_gaps():
    alpha = new_moltype.get_moltype("dna").degen_gapped_alphabet
    parent_seq = alpha.to_indices("ACTGC", alpha)
    got_ungap, got_empty_arr = new_aln.seq_to_gap_coords(
        parent_seq, moltype=new_moltype.get_moltype("dna")
    )
    assert numpy.array_equal(got_ungap, parent_seq)
    assert got_empty_arr.size == 0


@pytest.mark.xfail(reason="AlignedData not yet functional")
@pytest.mark.parametrize("i", range(3))
def test_seq_to_gap_coords_str(gap_seqs, i):
    seq, gap_coords = gap_seqs[i]
    got_ungapped, got_map = new_aln.seq_to_gap_coords(
        seq, moltype=new_moltype.get_moltype("dna")
    )
    assert got_ungapped == seq.replace("-", "")
    assert got_map.get_gap_coordinates() == gap_coords


@pytest.mark.xfail(reason="AlignedData not yet functional")
@pytest.mark.parametrize("i", range(3))
def test_seq_to_gap_coords_arr(gap_seqs, i):
    seq, gap_coords = gap_seqs[i]
    alpha = new_moltype.get_moltype("dna").degen_gapped_alphabet
    seq = alpha.to_indices(seq, alpha)  # convert to array repr
    got_ungapped, got_map = new_aln.seq_to_gap_coords(
        seq, moltype=new_moltype.get_moltype("dna")
    )
    assert numpy.array_equal(got_ungapped, seq[seq != 4])  # gap_char = 4
    assert got_map.get_gap_coordinates() == gap_coords


@pytest.mark.parametrize("moltype", ("dna", "rna", "protein", "protein_with_stop"))
def test_make_unaligned_seqs_dict(moltype):
    """test SequenceCollection constructor utility function"""
    data = {"a": "AGGCCC", "b": "AGAAAA"}
    got = new_aln.make_unaligned_seqs(data=data, moltype=moltype)
    assert isinstance(got, new_aln.SequenceCollection)
    assert got._seqs_data.make_seq == new_moltype.get_moltype(moltype).make_seq

    # should also work if seqs are arrays
    data = {"a": numpy.array(["AGGCCC"]), "b": numpy.array(["AGAAAA"])}
    got = new_aln.make_unaligned_seqs(data=data, moltype=moltype)
    assert isinstance(got, new_aln.SequenceCollection)
    assert got._seqs_data.make_seq == new_moltype.get_moltype(moltype).make_seq


@pytest.mark.parametrize("moltype", ("dna", "rna", "protein", "protein_with_stop"))
def test_make_unaligned_seqs_list(moltype):
    """test SequenceCollection constructor utility function"""
    data = ["AGGCCC", "AGAAAA"]
    got = new_aln.make_unaligned_seqs(data=data, moltype=moltype)
    assert isinstance(got, new_aln.SequenceCollection)
    assert got._seqs_data.make_seq == new_moltype.get_moltype(moltype).make_seq

    # should also work if seqs are arrays
    data = [numpy.array(["AGGCCC"]), numpy.array(["AGAAAA"])]
    got = new_aln.make_unaligned_seqs(data=data, moltype=moltype)
    assert isinstance(got, new_aln.SequenceCollection)
    assert got._seqs_data.make_seq == new_moltype.get_moltype(moltype).make_seq


def test_make_unaligned_seqs_raises():
    data = "AGTCCTGA"
    with pytest.raises(NotImplementedError):
        new_aln.make_unaligned_seqs(data=data, moltype="dna")


def test_make_unaligned_seqs_incompatible_moltype(dna_sd):
    with pytest.raises(ValueError):
        _ = new_aln.make_unaligned_seqs(data=dna_sd, moltype="rna")


def test_make_unaligned_seqs_no_seqs():
    data = {}
    with pytest.raises(ValueError):
        new_aln.make_unaligned_seqs(data=data, moltype="dna")


def test_sequence_collection_init(seqs):
    assert isinstance(seqs.seqs, new_aln.SeqsData)
    assert seqs.seqs.make_seq is not None


@pytest.mark.xfail(reason="todo: kath, havent implemented label_to_name")
def test_sequence_collection_label_to_name():
    """SequenceCollection init should allow name mapping function"""
    data = {"a": "AAAAA", "b": "BBBBB"}

    def f(x):
        return x.upper()

    seqs = new_aln.make_unaligned_seqs(data=data, moltype="dna", label_to_name=f)
    assert list(seqs.names) == ["A", "B"]
    assert seqs.seqs.names == ["A", "B"]


def test_iter_seqs_ragged_padded(ragged_padded):
    """SequenceCollection.iter_seqs() method should support reordering of seqs"""
    seqs = list(ragged_padded.iter_seqs())
    assert seqs == ["AAAAAA", "AAA---", "AAAA--"]
    seqs = list(ragged_padded.iter_seqs(seq_order=["b", "a", "a"]))
    assert seqs == ["AAA---", "AAAAAA", "AAAAAA"]
    assert seqs[1] == seqs[2]
    assert seqs[0] == ragged_padded.seqs["b"]


def test_sequence_collection_iter_seqs_ragged(ragged):
    """SequenceCollection iter_seqs() method should support reordering of seqs"""
    seqs = list(ragged.iter_seqs())
    assert seqs == ["AAAAAA", "AAA", "AAAA"]
    seqs = list(ragged.iter_seqs(seq_order=["b", "a", "a"]))
    assert seqs == ["AAA", "AAAAAA", "AAAAAA"]
    assert seqs[1] == seqs[2]
    assert seqs[0] == ragged.seqs["b"]


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


def test_sequence_collection_set_wrap_affects_repr_html():
    """the wrap argument affects the number of columns"""
    # indirectly tested via counting number of occurrences of 'class="label"'
    seqs = new_aln.make_unaligned_seqs(data={"a": "AAAAA", "b": "AAA--"}, moltype="dna")
    orig = seqs._repr_html_()
    seqs.set_repr_policy(wrap=3)  # break alignment into 2
    got = seqs._repr_html_()
    token = 'class="label"'
    assert got.count(token) == 2 * orig.count(token)

    # using environment variable
    env_name = "COGENT3_ALIGNMENT_REPR_POLICY"
    os.environ[env_name] = "wrap=2"
    seqs = new_aln.make_unaligned_seqs(data={"a": "AAAAA", "b": "AAA--"}, moltype="dna")
    got = seqs._repr_html_()
    os.environ.pop(env_name, None)
    assert got.count(token) == 3 * orig.count(token)


@pytest.mark.xfail(reason="todo: kath, new SequenceCollection does not support slicing")
def test_sequence_collection_repr_html(seqs):
    """exercises method normally invoked in notebooks"""
    seqs.set_repr_policy(num_seqs=5, num_pos=40)
    assert seqs[:2]._repr_policy == seqs._repr_policy
    row_a = '<tr><td class="label">a</td>'
    row_b = '<tr><td class="label">b</td>'
    # default order is longest sequence at top
    got = seqs._repr_html_()
    assert got.find(row_a) < got.find(row_b)
    # change order, a should now be last
    seqs.set_repr_policy(num_seqs=5, num_pos=40, ref_name="seq2")
    got = seqs._repr_html_()
    assert got.find(row_a) > got.find(row_b)


def test_sequence_collection_repr_html_applies_policy(seqs):
    # tests repr policy has been successfully applied
    seqs.set_repr_policy(num_seqs=2)
    got = seqs._repr_html_()
    assert got.count("</tr>") == 3
    seqs.set_repr_policy(num_seqs=3)
    got = seqs._repr_html_()
    assert got.count("</tr>") == 4
    seqs.set_repr_policy(num_seqs=len(seqs.seqs))
    got = seqs._repr_html_()
    assert got.count("</tr>") == len(seqs.seqs) + 1


def test_sequence_collection_repr_html_correct_num_seqs(seqs):
    # tests _repr_html_ displays correct number of sequences
    got = seqs._repr_html_()
    seq_lens = numpy.array(
        [len(seqs.seqs.get_seq_str(seqid=name)) for name in seqs.names]
    )
    seq_stats = {
        "min": seq_lens.min(),
        "median": numpy.median(seq_lens),
        "max": seq_lens.max(),
    }
    assert (
        f"{seqs.num_seqs} x {{min={seq_lens.min()}, median={numpy.median(seq_lens)}, max={seq_lens.max()}}}"
        in got.splitlines()[-2]
    )


def test_sequence_collection_set_repr_policy_no_input(seqs):
    """repr_policy should remain unchanged"""
    seqs.set_repr_policy(num_seqs=None, num_pos=None)
    assert seqs._repr_policy == dict(
        num_seqs=10, num_pos=60, ref_name="longest", wrap=60
    )


def test_sequence_collection_set_repr_policy_invalid_input(seqs):
    """repr_policy should remain unchanged"""
    invalid_args = (
        dict(num_seqs="foo", err=TypeError),
        dict(num_pos=4.2, err=TypeError),
        dict(ref_name="blah", err=ValueError),
        dict(wrap=3.1, err=TypeError),
    )

    for arg in invalid_args:
        err = arg.pop("err")
        with pytest.raises(err):
            seqs.set_repr_policy(**arg)
        assert seqs._repr_policy == dict(
            num_seqs=10, num_pos=60, ref_name="longest", wrap=60
        )


def test_sequence_collection_set_repr_policy_valid_input(seqs):
    """repr_policy should be set to new values"""
    seqs.set_repr_policy(num_seqs=5, num_pos=40, ref_name="seq1", wrap=10)
    assert seqs._repr_policy == dict(num_seqs=5, num_pos=40, ref_name="seq1", wrap=10)


@pytest.mark.xfail(reason="todo: kath, need support for __eq__ betwen collections")
def test_sequence_collection_take_seqs(ragged_padded):
    """SequenceCollection take_seqs should return new SequenceCollection with selected seqs."""
    a = ragged_padded.take_seqs(list("bc"))
    assert isinstance(a, new_aln.SequenceCollection)
    assert a == {"b": "AAA---", "c": "AAAA--"}
    # should be able to negate
    a = ragged_padded.take_seqs(list("bc"), negate=True)
    assert a == {"a": "AAAAAA"}


@pytest.mark.xfail(reason="todo: kath, need support for __eq__ betwen collections")
def test_sequence_collection_take_seqs_str(ragged_padded):
    """string arg to SequenceCollection take_seqs should work."""
    a = ragged_padded.take_seqs("a")
    assert a == {"a": "AAAAAA"}

    # should be able to negate
    a = ragged_padded.take_seqs("a", negate=True)
    assert isinstance(a, new_aln.SequenceCollection)
    assert a == {"b": "AAA---", "c": "AAAA--"}


def test_sequence_collection_take_seqs_info():
    """take_seqs should preserve info attribute"""
    moltype = new_moltype.get_moltype("dna")
    orig = new_aln.make_unaligned_seqs(
        data={"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"},
        moltype=moltype,
        info={"key": "value"},
    )
    subset = orig.take_seqs(list("ab"))
    assert set(subset.info) == set(orig.info)


def test_sequence_collection_take_seqs_moltype():
    """take_seqs should preserve the MolType"""
    moltype = new_moltype.get_moltype("dna")
    orig = new_aln.make_unaligned_seqs(
        data={"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"}, moltype=moltype
    )
    subset = orig.take_seqs(list("ab"))
    assert set(subset.moltype) == set(orig.moltype)


def test_sequence_collection_take_seqs_empty_names():
    moltype = new_moltype.get_moltype("dna")
    orig = new_aln.make_unaligned_seqs(
        data={"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"}, moltype=moltype
    )
    subset = orig.take_seqs([])
    assert subset == {}


def test_sequence_collection_num_seqs():
    """SequenceCollection.num_seqs should count seqs."""

    seqs = new_aln.make_unaligned_seqs(
        data={"seq1": "ACGU", "seq2": "CGUA", "seq3": "CCGU"}, moltype="rna"
    )
    assert seqs.num_seqs == 3


@pytest.mark.parametrize("index", [(0, "seq1"), (1, "seq2"), (2, "seq3")])
def test_sequence_collection_getitem(seqs, index):
    got1 = seqs.seqs[index[0]]
    got2 = seqs.seqs[index[1]]

    assert isinstance(got1, new_seq.Sequence)
    assert isinstance(got2, new_seq.Sequence)
    assert got1 == got2


# Functions for testing SequenceCollection methods that accept a function as an argument
def is_long(x):
    return len(x) > 10


def is_med(x):
    return len(str(x).replace("-", "")) > 3


def is_any(x):
    return len(x) > 0


def test_get_seq_indices(ragged_padded):
    """SequenceCollection.get_seq_indices should return names of seqs where f(row) is True"""

    assert ragged_padded.get_seq_indices(is_long) == []
    assert ragged_padded.get_seq_indices(is_med) == ["a", "c"]
    # return order should reflect names when updated
    ragged_padded.names = ["b", "c", "a"]
    assert ragged_padded.get_seq_indices(is_med) == ["c", "a"]
    assert ragged_padded.get_seq_indices(is_med, negate=True) == ["b"]
    assert ragged_padded.get_seq_indices(is_any) == ["b", "c", "a"]
    assert ragged_padded.get_seq_indices(is_any, negate=True) == []


@pytest.mark.xfail(reason="todo: kath, need support for __eq__ betwen collections")
def test_take_seqs_if(ragged_padded):
    """SequenceCollection take_seqs_if should return seqs where f(seq) is True"""

    assert ragged_padded.take_seqs_if(is_long) == {}
    assert ragged_padded.take_seqs_if(is_any) == ragged_padded
    assert isinstance(ragged_padded.take_seqs_if(is_med), new_aln.SequenceCollection)
    assert ragged_padded.take_seqs_if(is_any, negate=True) == {}
