import os

import numpy
import pytest

import cogent3.core.new_alignment as new_aln
import cogent3.core.new_alphabet as new_alpha
import cogent3.core.new_moltype as new_moltype
import cogent3.core.new_sequence as new_seq

# todo: kath, update this to using new Sequence load_seq when it is ready
from cogent3 import load_seq
from cogent3.core.annotation_db import GffAnnotationDb, load_annotations


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
    return new_aln.make_unaligned_seqs(data, moltype="dna")


@pytest.fixture
def ragged_padded():
    return new_aln.make_unaligned_seqs(
        {"a": "AAAAAA", "b": "AAA---", "c": "AAAA--"}, moltype="dna"
    )


@pytest.fixture
def ragged():
    return new_aln.make_unaligned_seqs(
        {"a": "AAAAAA", "b": "AAA", "c": "AAAA"}, moltype="dna"
    )


@pytest.fixture
def unordered():
    return new_aln.make_unaligned_seqs({"a": "AAAAA", "c": "CCCCC"}, moltype="dna")


@pytest.fixture
def ordered1():
    seqs = new_aln.make_unaligned_seqs({"a": "AAAAA", "c": "CCCCC"}, moltype="dna")
    seqs.names = ["a", "c"]
    return seqs


@pytest.fixture
def ordered2():
    seqs = new_aln.make_unaligned_seqs({"a": "AAAAA", "c": "CCCCC"}, moltype="dna")
    seqs.names = ["c", "a"]
    return seqs


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


@pytest.fixture(scope="function")
def gb_db(DATA_DIR):
    return load_annotations(path=DATA_DIR / "annotated_seq.gb")


@pytest.fixture(scope="function")
def gff_db(DATA_DIR):
    return load_annotations(path=DATA_DIR / "simple.gff")


@pytest.fixture(scope="session")
def seqcoll_db():
    fasta_path = os.path.join("data/c_elegans_WS199_dna_shortened.fasta")
    gff3_path = os.path.join("data/c_elegans_WS199_shortened_gff.gff3")
    seq = load_seq(fasta_path, moltype="dna")
    seq_coll = new_aln.make_unaligned_seqs({seq.name: seq}, moltype="dna")
    seq_coll.annotate_from_gff(gff3_path)
    return seq_coll


def test_seqs_data_default_attributes(dna_sd: new_aln.SeqsData):
    assert dna_sd.names == ["seq1", "seq2", "seq3"]
    assert isinstance(dna_sd.alphabet, new_alpha.CharAlphabet)


@pytest.mark.xfail(reason="SeqsData currently expects correctly formatted data")
def test_seqs_data_seq_if_str(seq1: str, alpha):
    with pytest.raises(NotImplementedError):
        new_aln.SeqsData(data=seq1, alphabet=alpha)


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
    sd = new_aln.SeqsData(data=d, alphabet=alpha)
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
    sd = new_aln.SeqsData(data=str_seqs_dict, alphabet=alpha)
    seq = str_seqs_dict[seqid]
    seq_len = len(seq)
    got = sd.get_seq_view(seqid)
    assert isinstance(got, new_aln.SeqDataView)
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
    sd = new_aln.SeqsData(data=str_seqs_dict, alphabet=alpha)
    got = sd.names
    # returns iterator
    assert list(got) == list(expect)


def test_seqs_data_seq_lengths(str_seqs_dict, arr_seqs_dict, alpha):
    expect = {k: len(v) for k, v in str_seqs_dict.items()}
    sd = new_aln.SeqsData(data=str_seqs_dict, alphabet=alpha)
    got = sd.seq_lengths()
    assert got == expect

    expect = {k: len(v) for k, v in arr_seqs_dict.items()}
    sd = new_aln.SeqsData(data=arr_seqs_dict, alphabet=alpha)
    got = sd.seq_lengths()
    assert got == expect


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
    sd = new_aln.SeqsData(data=str_seqs_dict, alphabet=alpha)
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
    assert got.names == ["seq2", "seq1"]
    assert got.get_seq_str(seqid="seq1") == dna_sd.get_seq_str(seqid="seq1")
    assert got.get_seq_str(seqid="seq2") == dna_sd.get_seq_str(seqid="seq2")

    got = dna_sd.subset(["seq3"])
    assert got.names == ["seq3"]
    assert got.get_seq_str(seqid="seq3") == dna_sd.get_seq_str(seqid="seq3")


def test_seqs_data_subset_raises(dna_sd):
    # should handle missing seqids
    got = dna_sd.subset(["seq1", "seq4"])
    assert got.names == ["seq1"]

    # but will raise error if resulting in empty selection
    with pytest.raises(ValueError):
        _ = dna_sd.subset(["seq4"])


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
    got = new_aln.make_unaligned_seqs(data, moltype=moltype)
    assert isinstance(got, new_aln.SequenceCollection)
    assert got._seqs_data.make_seq == new_moltype.get_moltype(moltype).make_seq

    # should also work if seqs are arrays
    data = {"a": numpy.array(["AGGCCC"]), "b": numpy.array(["AGAAAA"])}
    got = new_aln.make_unaligned_seqs(data, moltype=moltype)
    assert isinstance(got, new_aln.SequenceCollection)
    assert got._seqs_data.make_seq == new_moltype.get_moltype(moltype).make_seq


@pytest.mark.parametrize("moltype", ("dna", "rna", "protein", "protein_with_stop"))
def test_make_unaligned_seqs_list(moltype):
    """test SequenceCollection constructor utility function"""
    data = ["AGGCCC", "AGAAAA"]
    got = new_aln.make_unaligned_seqs(data, moltype=moltype)
    assert isinstance(got, new_aln.SequenceCollection)
    assert got._seqs_data.make_seq == new_moltype.get_moltype(moltype).make_seq

    # should also work if seqs are arrays
    data = [numpy.array(["AGGCCC"]), numpy.array(["AGAAAA"])]
    got = new_aln.make_unaligned_seqs(data, moltype=moltype)
    assert isinstance(got, new_aln.SequenceCollection)
    assert got._seqs_data.make_seq == new_moltype.get_moltype(moltype).make_seq


def test_make_unaligned_seqs_label_to_name():
    """test SequenceCollection constructor utility function"""
    data = {"a": "AGGCCC", "b": "AGAAAA"}

    def f(x):
        return x.upper()

    got = new_aln.make_unaligned_seqs(data, moltype="dna", label_to_name=f)
    assert list(got.names) == ["A", "B"]
    assert got._seqs_data.make_seq == new_moltype.get_moltype("dna").make_seq


def test_make_unaligned_seqs_raises():
    data = "AGTCCTGA"
    with pytest.raises(NotImplementedError):
        new_aln.make_unaligned_seqs(data, moltype="dna")


def test_make_unaligned_seqs_incompatible_moltype(dna_sd):
    with pytest.raises(ValueError):
        _ = new_aln.make_unaligned_seqs(dna_sd, moltype="rna")


def test_make_unaligned_seqs_no_seqs():
    data = {}
    with pytest.raises(ValueError):
        new_aln.make_unaligned_seqs(data, moltype="dna")


def test_sequence_collection_init(seqs):
    assert isinstance(seqs.seqs, new_aln.SeqsData)
    assert seqs.seqs.make_seq is not None


def test_sequence_collection_names_is_list():
    """expected to be a list"""
    seqs = new_aln.make_unaligned_seqs({"a": b"AAAAA", "b": b"TTTTT"}, moltype="dna")
    assert isinstance(seqs.names, list)


def test_sequence_collection_names(seqs):
    assert seqs.names == ["seq1", "seq2", "seq3"]
    seqs.names = ["seq2", "seq3", "seq1"]
    assert seqs.names == ["seq2", "seq3", "seq1"]
    seqs.names = ["seq1", "seq2"]
    assert seqs.names == ["seq1", "seq2"]
    with pytest.raises(ValueError):
        seqs.names = ["seq1", "seq2", "seq3", "seq4"]


def test_init_seq():
    """SequenceCollection init from list of sequences should use indices as keys"""
    seqs = ["TTTTT", "CCCCC", "GGGGG"]
    a = new_aln.make_unaligned_seqs(seqs, moltype="dna")
    assert len(a.seqs) == 3
    assert a.seqs["seq_0"] == "TTTTT"
    assert a.seqs["seq_1"] == "CCCCC"
    assert a.seqs["seq_2"] == "GGGGG"
    assert a.names == ["seq_0", "seq_1", "seq_2"]


@pytest.mark.xfail(
    reason="todo: decide if we support init from list of (key,val) pairs"
)
def test_init_pairs():
    """SequenceCollection init from list of (key,val) pairs should work correctly"""
    seqs = [["a", "AAA"], ["t", "TTT"], ["c", "CCC"]]
    a = new_aln.make_unaligned_seqs(seqs, moltype="dna")
    assert len(a.seqs) == 3
    assert a.seqs["a"] == "AAA"
    assert a.seqs["t"] == "TTT"
    assert a.seqs["c"] == "CCC"
    assert a.names == ["a", "t", "c"]


def test_init_ordered(ordered1, ordered2):
    """SequenceCollection should iterate over seqs correctly even if ordered"""
    first = ordered1
    sec = ordered2

    assert first.names == ["a", "c"]
    assert sec.names == ["c", "a"]


def test_sequence_collection_init_ambig():
    """SequenceCollection should tolerate ambiguous chars"""
    _ = new_aln.make_unaligned_seqs(["AAA", "CCC"], moltype="dna")
    _ = new_aln.make_unaligned_seqs(["ANS", "CWC"], moltype="dna")
    _ = new_aln.make_unaligned_seqs(["A-A", "CC-"], moltype="dna")
    _ = new_aln.make_unaligned_seqs(["A?A", "CC-"], moltype="dna")


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
    seqs = new_aln.make_unaligned_seqs(data, moltype="dna")
    assert (
        repr(seqs)
        == "2x (ENSMUSG00000039616[GCCCTTCAAA...], ENSMUSG00000056468[GCCAGGGGGA...]) dna seqcollection"
    )

    data = {
        "ENSMUSG00000039616": "GCCCTTCAAATTT",
        "ENSMUSG00000056468": "GCCAGGGGGAAAAGGGAGAA",
    }
    seqs = new_aln.make_unaligned_seqs(data, moltype="dna")
    assert (
        repr(seqs)
        == "2x (ENSMUSG00000039616[GCCCTTCAAA...], ENSMUSG00000056468[GCCAGGGGGA...]) dna seqcollection"
    )

    data = {
        "a": "TCGAT",
    }
    seqs = new_aln.make_unaligned_seqs(data, moltype="dna")
    assert repr(seqs) == "1x (a[TCGAT]) dna seqcollection"

    data = {
        "a": "TCGAT" * 2,
    }
    seqs = new_aln.make_unaligned_seqs(data, moltype="dna")
    assert repr(seqs) == "1x (a[TCGATTCGAT]) dna seqcollection"


def test_sequence_collection_set_wrap_affects_repr_html():
    """the wrap argument affects the number of columns"""
    # indirectly tested via counting number of occurrences of 'class="label"'
    seqs = new_aln.make_unaligned_seqs({"a": "AAAAA", "b": "AAA--"}, moltype="dna")
    orig = seqs._repr_html_()
    seqs.set_repr_policy(wrap=3)  # break alignment into 2
    got = seqs._repr_html_()
    token = 'class="label"'
    assert got.count(token) == 2 * orig.count(token)

    # using environment variable
    env_name = "COGENT3_ALIGNMENT_REPR_POLICY"
    os.environ[env_name] = "wrap=2"
    seqs = new_aln.make_unaligned_seqs({"a": "AAAAA", "b": "AAA--"}, moltype="dna")
    got = seqs._repr_html_()
    os.environ.pop(env_name, None)
    assert got.count(token) == 3 * orig.count(token)


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


def test_sequence_collection_take_seqs(ragged_padded):
    """SequenceCollection take_seqs should return new SequenceCollection with selected seqs."""
    a = ragged_padded.take_seqs(list("bc"))
    assert isinstance(a, new_aln.SequenceCollection)
    assert a.names == ["b", "c"]
    assert str(a.seqs["b"]) == "AAA---"
    assert str(a.seqs["c"]) == "AAAA--"
    assert a.num_seqs == 2
    # should be able to negate
    a = ragged_padded.take_seqs(list("bc"), negate=True)
    assert a.names == ["a"]
    assert str(a.seqs["a"]) == "AAAAAA"


def test_sequence_collection_take_seqs_str(ragged_padded):
    """string arg to SequenceCollection take_seqs should work."""
    a = ragged_padded.take_seqs("a")
    assert a.names == ["a"]
    assert a.num_seqs == 1

    # should be able to negate
    a = ragged_padded.take_seqs("a", negate=True)
    assert isinstance(a, new_aln.SequenceCollection)
    assert a.names == ["b", "c"]
    assert a.num_seqs == 2


def test_sequence_collection_take_seqs_info():
    """take_seqs should preserve info attribute"""
    moltype = new_moltype.get_moltype("dna")
    orig = new_aln.make_unaligned_seqs(
        {"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"},
        moltype=moltype,
        info={"key": "value"},
    )
    subset = orig.take_seqs(list("ab"))
    assert set(subset.info) == set(orig.info)


def test_sequence_collection_take_seqs_moltype():
    """take_seqs should preserve the MolType"""
    moltype = new_moltype.get_moltype("dna")
    orig = new_aln.make_unaligned_seqs(
        {"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"}, moltype=moltype
    )
    subset = orig.take_seqs(list("ab"))
    assert set(subset.moltype) == set(orig.moltype)


def test_sequence_collection_take_seqs_empty_names():
    moltype = new_moltype.get_moltype("dna")
    orig = new_aln.make_unaligned_seqs(
        {"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"}, moltype=moltype
    )
    with pytest.raises(ValueError):
        _ = orig.take_seqs([])


def test_sequence_collection_num_seqs():
    """SequenceCollection.num_seqs should count seqs."""

    seqs = new_aln.make_unaligned_seqs(
        {"seq1": "ACGU", "seq2": "CGUA", "seq3": "CCGU"}, moltype="rna"
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


def test_get_seq_names_if(ragged_padded):
    """SequenceCollection.get_seq_names_if should return names of seqs where f(row) is True"""

    assert ragged_padded.get_seq_names_if(is_long) == []
    assert ragged_padded.get_seq_names_if(is_med) == ["a", "c"]
    # return order should reflect names when updated
    ragged_padded.names = ["b", "c", "a"]
    assert ragged_padded.get_seq_names_if(is_med) == ["c", "a"]
    assert ragged_padded.get_seq_names_if(is_med, negate=True) == ["b"]
    assert ragged_padded.get_seq_names_if(is_any) == ["b", "c", "a"]
    assert ragged_padded.get_seq_names_if(is_any, negate=True) == []


def test_take_seqs_if(ragged_padded):
    """SequenceCollection take_seqs_if should return seqs where f(seq) is True"""

    with pytest.raises(ValueError):
        ragged_padded.take_seqs_if(is_long)
    with pytest.raises(ValueError):
        ragged_padded.take_seqs_if(is_any, negate=True)
    got = ragged_padded.take_seqs_if(is_any)
    assert got.names == ragged_padded.names
    assert got.num_seqs == ragged_padded.num_seqs
    assert isinstance(ragged_padded.take_seqs_if(is_med), new_aln.SequenceCollection)

    got = ragged_padded.take_seqs_if(is_med, negate=True)
    assert got.names == ["b"]
    assert got.num_seqs == 1


def test_sequence_collection_to_dict():
    """SequenceCollection.to_dict should return dict of strings (not obj)"""
    data = {"seq1": "GATTTT", "seq2": "GATC??"}
    seqs = new_aln.make_unaligned_seqs(data, moltype="dna")
    assert seqs.to_dict() == data
    for i in list(seqs.to_dict().values()):
        assert isinstance(i, str)


def test_sequence_collection_get_seq():
    """SequenceCollection.get_seq should return specified seq"""
    seqs = new_aln.make_unaligned_seqs(
        {"seq1": "GATTTT", "seq2": "GATC??"}, moltype="dna"
    )
    assert seqs.get_seq("seq1") == "GATTTT"
    with pytest.raises(KeyError):
        seqs.get_seq("seqx")


@pytest.mark.xfail(reason="AttributeError: 'MolType' object has no attribute 'degap'")
def test_sequence_collection_degap():
    """SequenceCollection.degap should strip gaps from each seq"""
    seqs = new_aln.make_unaligned_seqs({"s1": "ATGRY?", "s2": "T-AG??"}, moltype="dna")
    assert seqs.degap().to_dict() == {"s1": "ATGRY", "s2": "TAG"}


def test_sequence_collection_to_fasta():
    """SequenceCollection should return correct FASTA string"""
    seqs = new_aln.make_unaligned_seqs(["AAA", "CCC"], moltype="dna")
    assert seqs.to_fasta() == ">seq_0\nAAA\n>seq_1\nCCC\n"
    assert seqs.to_fasta(block_size=2) == ">seq_0\nAA\nA\n>seq_1\nCC\nC\n"

    seqs = new_aln.make_unaligned_seqs(["GCATGCAT", "TCAGACGT"], moltype="dna")
    assert seqs.to_fasta() == ">seq_0\nGCATGCAT\n>seq_1\nTCAGACGT\n"
    assert seqs.to_fasta(block_size=4) == ">seq_0\nGCAT\nGCAT\n>seq_1\nTCAG\nACGT\n"
    assert seqs.to_fasta(block_size=3) == ">seq_0\nGCA\nTGC\nAT\n>seq_1\nTCA\nGAC\nGT\n"


def test_sequence_collection_to_phylip():
    """SequenceCollection should return PHYLIP string format correctly"""
    align_norm = new_aln.make_unaligned_seqs(
        [
            "ACDEFGHIKLMNPQRSTUVWY-",
            "ACDEFGHIKLMNPQRSUUVWF-",
            "ACDEFGHIKLMNPERSKUVWC-",
            "ACNEFGHIKLMNPQRS-UVWP-",
        ],
        moltype="protein",
    )

    assert (
        align_norm.to_phylip()
        == """4  22\nseq_0     ACDEFGHIKLMNPQRSTUVWY-\nseq_1     ACDEFGHIKLMNPQRSUUVWF-\nseq_2     ACDEFGHIKLMNPERSKUVWC-\nseq_3     ACNEFGHIKLMNPQRS-UVWP-\n"""
    )


@pytest.mark.parametrize(
    "gc,seqs",
    (
        (1, ("TCCTGA", "GATTT?")),
        (1, ("ACGTAA---", "ACGAC----", "ACGCAATGA")),
        (2, ("GATTTT", "TCCAGG")),
    ),
)
def test_has_terminal_stop_true(gc, seqs):
    data = {f"s{i}": s for i, s in enumerate(seqs)}
    seqs = new_aln.make_unaligned_seqs(data, moltype="dna")
    assert seqs.has_terminal_stop(gc=gc)


@pytest.mark.parametrize(
    "gc,seqs",
    ((1, ("TCCTCA", "GATTTT")), (2, ("GATTTT", "TCCCGG")), (1, ("CCTCA", "ATTTT"))),
)
def test_has_terminal_stop_false(gc, seqs):
    data = {f"s{i}": s for i, s in enumerate(seqs)}
    seqs = new_aln.make_unaligned_seqs(data, moltype="dna")
    assert not seqs.has_terminal_stop(gc=gc)


def test_get_similar():
    data = {
        "a": "AAAAAAAAAA",  # target
        "b": "AAAAAAAAAA",  # identical
        "c": "AAAAAAAAAT",  # 90% identical
        "d": "AAAAAAAATT",  # 80% identical
        "e": "AATTTTTTTT",  # 20% identical
        "f": "TTTTTTTTTT",  # 0% identical
    }
    seqs = new_aln.make_unaligned_seqs(data, moltype="dna")
    target = seqs.get_seq("a")
    got = seqs.get_similar(target, min_similarity=0.8)
    assert got.names == ["a", "b", "c", "d"]
    got = seqs.get_similar(target, min_similarity=0.81)
    assert got.names == ["a", "b", "c"]
    got = seqs.get_similar(target, min_similarity=0.75, max_similarity=0.9)
    assert got.names == ["c", "d"]
    got = seqs.get_similar(target, min_similarity=0.75, max_similarity=0.89)
    assert got.names == ["d"]


def test_sequence_collection_init_annotated_seqs():
    """correctly construct from list with annotated seq"""
    seq = new_moltype.DNA.make_seq(seq="GCCAGGGGGGAAAG-GGAGAA", name="seq1")
    _ = seq.add_feature(biotype="exon", name="name", spans=[(4, 10)])
    coll = new_aln.make_unaligned_seqs([seq], moltype="dna")
    features = list(coll.get_features(biotype="exon"))
    assert len(features) == 1


def test_sequence_collection_make_unaligned_seqs_annotated():
    """annotate_from_gff should work on data from gff3 files"""
    fasta_path = os.path.join("data/c_elegans_WS199_dna_shortened.fasta")
    gff3_path = os.path.join("data/c_elegans_WS199_shortened_gff.gff3")
    seq = load_seq(fasta_path, moltype="dna")
    seq = new_moltype.DNA.make_seq(seq=str(seq), name="seq1")
    seq.annotate_from_gff(gff3_path)


def test_sequence_collection_annotation_db_assign_none():
    """assigning None to annotation_db breaks conection"""
    seq_coll = new_aln.make_unaligned_seqs(
        {"seq1": "ACGU", "seq2": "CGUA", "test_seq": "CCGU"}, moltype="rna"
    )
    seq_coll.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])
    seq_coll.add_feature(seqid="seq2", biotype="xyzzz", name="abc", spans=[(1, 2)])
    seq_coll.annotation_db = None
    assert seq_coll.annotation_db is None


def test_sequence_collection_get_annotations_from_any_seq():
    """get_annotations_from_any_seq returns correct annotations"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = new_aln.make_unaligned_seqs(data, moltype="dna")
    db = GffAnnotationDb()
    db.add_feature(seqid="seq1", biotype="exon", name="annotation1", spans=[(3, 8)])
    db.add_feature(seqid="seq2", biotype="exon", name="annotation2", spans=[(1, 2)])
    db.add_feature(seqid="seq3", biotype="exon", name="annotation3", spans=[(3, 6)])
    seqs.annotation_db = db
    got = list(seqs.get_features())
    assert len(got) == 3
    assert "biotype='exon', name='annotation1', map=[3:8]/9" in str(got[0])
    assert "biotype='exon', name='annotation2', map=[1:2]/9" in str(got[1])
    assert "biotype='exon', name='annotation3', map=[3:6]/9" in str(got[2])

    got = list(seqs.get_features(name="annotation1"))
    assert len(got) == 1
    assert "biotype='exon', name='annotation1', map=[3:8]/9" in str(got[0])

    got = list(seqs.get_features(biotype="exon", name="annotation2"))
    assert len(got) == 1
    assert "biotype='exon', name='annotation2', map=[1:2]/9" in str(got[0])

    got = list(seqs.get_features(name="annotation3"))
    assert len(got) == 1
    assert "biotype='exon', name='annotation3', map=[3:6]/9" in str(got[0])


def test_sequence_collection_copy_annotations(gff_db):
    """copy_annotations copies records from annotation db"""
    seq_coll = new_aln.make_unaligned_seqs(
        {"seq1": "ACGU", "seq2": "CGUA", "test_seq": "CCGU"}, moltype="rna"
    )
    seq_coll.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])
    seq_coll.add_feature(seqid="seq2", biotype="xyzzz", name="abc", spans=[(1, 2)])
    expect = seq_coll.annotation_db.num_matches() + gff_db.num_matches()
    seq_coll.copy_annotations(gff_db)
    assert seq_coll.annotation_db.num_matches() == expect


def test_sequence_collection_get_seq_annotated():
    """SequenceCollection.get_seq should return specified seq"""
    seqs = new_aln.make_unaligned_seqs(
        {"seq1": "GATTTT", "seq2": "GATC??"}, moltype="dna"
    )
    seqs.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])

    with_annos = seqs.get_seq("seq1", copy_annotations=True)
    assert len(with_annos.annotation_db) == 1

    without_annos = seqs.get_seq("seq1", copy_annotations=False)
    assert len(without_annos.annotation_db) == 0


def _make_seq(name):
    raw_seq = "AACCCAAAATTTTTTGGGGGGGGGGCCCC"
    cds = (15, 25)
    utr = (12, 15)
    seq = new_moltype.DNA.make_seq(seq=raw_seq, name=name)
    seq.add_feature(biotype="CDS", name="CDS", spans=[cds])
    seq.add_feature(biotype="5'UTR", name="5' UTR", spans=[utr])
    return seq


@pytest.mark.xfail(
    reason="construction of Sequences from SequenceCollection does not propagate \
    annotations to get a sequence with annotations, use get_seq with copy_annotations=True"
)
def test_sequence_collection_init_seqs_have_annotations():
    """annotations on input seqs correctly merged and propagated"""

    seq_coll = new_aln.make_unaligned_seqs(
        {"seq1": _make_seq("seq1"), "seq2": _make_seq("seq2")}, moltype="dna"
    )
    coll_db = seq_coll.annotation_db
    assert len(coll_db) == 4
    for seq in seq_coll.seqs:
        db = seq.annotation_db
        assert db is coll_db


@pytest.mark.xfail(
    reason="annotating sequences does not propagate annotations to the SequenceCollection"
)
def test_add_to_seq_updates_coll():
    """annotating a seq updates the db of the propagated"""
    seq_coll = new_aln.make_unaligned_seqs(
        {
            "x": "AACCCAAAATTTTTTGGGGGGGGGGCCCC",
            "y": "AACCCAAAATTTTTTGGGGGGGGGGCCCC",
        },
        moltype="dna",
    )
    x = seq_coll.get_seq("x")
    assert len(seq_coll.annotation_db) == len(x.annotation_db) == 0
    x.add_feature(biotype="exon", name="E1", spans=[(3, 8)])
    assert len(seq_coll.annotation_db) == len(x.annotation_db) == 1


@pytest.mark.xfail(reason="load_seq creating old-style Sequence")
def test_copy_annotations_incompat_fails(seqcoll_db, gb_db):
    """copy_annotations copies records from annotation db"""
    db = seqcoll_db.annotation_db
    seqcoll_db.annotation_db = db
    with pytest.raises(TypeError):
        seqcoll_db.copy_annotations(gb_db)


@pytest.mark.xfail(reason="load_seq creating old-style Sequence")
def test_copy_annotations_incompat_type_fails(seqcoll_db, seqs):
    """copy_annotations copies records from annotation db"""

    with pytest.raises(TypeError):
        seqcoll_db.copy_annotations(seqs)
