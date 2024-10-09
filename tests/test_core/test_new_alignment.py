import json
import os
import pathlib
import re
from warnings import catch_warnings, filterwarnings

import numpy
import pytest

from cogent3 import get_app, load_unaligned_seqs, open_
from cogent3._version import __version__
from cogent3.core import new_alignment, new_alphabet, new_moltype, new_sequence
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import GffAnnotationDb, load_annotations
from cogent3.util.deserialise import deserialise_object
from cogent3.util.misc import get_object_provenance


@pytest.mark.parametrize(
    "func", (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs)
)
@pytest.mark.parametrize("index", (True, False))
def test_indexing_seqs_prop(func, index):
    names = ["seq1", "seq2", "seq3"]
    seqs = "GGGTAC", "GTTTGC", "ACGTAC"
    raw = dict(zip(names, seqs))
    obj = func(raw, moltype="dna")
    got = obj.seqs[0 if index else "seq1"]
    assert str(got) == raw["seq1"]


@pytest.mark.parametrize(
    "func", (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs)
)
def test_indexing_seqs_repr(func):
    names = ["seq1", "seq2", "seq3"]
    seqs = "GGGTAC", "GTTTGC", "ACGTAC"
    raw = dict(zip(names, seqs))
    obj = func(raw, moltype="dna")
    got = repr(obj.seqs)
    print(got, type(obj.seqs))


@pytest.fixture(scope="session")
def tmp_path(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp_path")


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
def dna_alphabet():
    moltype = new_moltype.get_moltype("dna")
    return moltype.degen_gapped_alphabet


@pytest.fixture
def dna_moltype():
    return new_moltype.get_moltype("dna")


@pytest.fixture
def dna_make_seq():
    return new_moltype.get_moltype("dna").make_seq


@pytest.fixture
def rna_moltype():
    return new_moltype.get_moltype("rna")


@pytest.fixture
def dna_sd(str_seqs_dict: dict[str, str], dna_alphabet):
    return new_alignment.SeqsData(data=str_seqs_dict, alphabet=dna_alphabet)


@pytest.fixture
def int_arr():
    return numpy.arange(17, dtype=numpy.uint8)


@pytest.fixture
def sdv_s2(dna_sd: new_alignment.SeqsData) -> new_alignment.SeqDataView:
    return dna_sd.get_view("seq2")


@pytest.fixture(scope="function")
def seqs() -> new_alignment.SequenceCollection:
    data = {"seq1": "AAAAAA", "seq2": "TTTT", "seq3": "ATTCCCC"}
    return new_alignment.make_unaligned_seqs(data, moltype="dna")


@pytest.fixture
def ragged_padded_dict():
    return {"a": "AAAAAA", "b": "AAA---", "c": "AAAA--"}


@pytest.fixture
def ragged_padded(ragged_padded_dict):
    return new_alignment.make_unaligned_seqs(ragged_padded_dict, moltype="dna")


@pytest.fixture
def ragged():
    return new_alignment.make_unaligned_seqs(
        {"a": "AAAAAA", "b": "AAA", "c": "AAAA"}, moltype="dna"
    )


@pytest.fixture
def unordered():
    return new_alignment.make_unaligned_seqs(
        {"a": "AAAAA", "c": "CCCCC"}, moltype="dna"
    )


@pytest.fixture
def ordered1():
    seqs = new_alignment.make_unaligned_seqs(
        {"a": "AAAAA", "c": "CCCCC"}, moltype="dna"
    )
    seqs.names = ["a", "c"]
    return seqs


@pytest.fixture
def ordered2():
    seqs = new_alignment.make_unaligned_seqs(
        {"a": "AAAAA", "c": "CCCCC"}, moltype="dna"
    )
    seqs.names = ["c", "a"]
    return seqs


@pytest.fixture(scope="function")
def gb_db(DATA_DIR):
    return load_annotations(path=DATA_DIR / "annotated_seq.gb")


@pytest.fixture(scope="function")
def gff_db(DATA_DIR):
    return load_annotations(path=DATA_DIR / "simple.gff")


@pytest.fixture(scope="function")
def seqcoll_db(DATA_DIR):
    from cogent3.parse.fasta import iter_fasta_records

    # NOTE: using iter_fasta_records directly because the c elegans
    # file has a lower case sequence and the iter_fasta_records function
    # coerces by default

    fasta_path = DATA_DIR / "c_elegans_WS199_dna_shortened.fasta"
    data = dict(iter_fasta_records(fasta_path))
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    seq_coll.annotation_db = load_annotations(
        path=DATA_DIR / "c_elegans_WS199_shortened_gff.gff3"
    )
    return seq_coll


def make_typed(seq, data_type, moltype):
    if data_type is numpy.ndarray:
        seq = moltype.most_degen_alphabet().to_indices(seq)
    elif data_type is bytes:
        seq = seq.encode("utf-8")
    return seq


def test_seqs_data_default_attributes(dna_sd: new_alignment.SeqsData):
    assert dna_sd.names == ["seq1", "seq2", "seq3"]
    assert isinstance(dna_sd.alphabet, new_alphabet.CharAlphabet)


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
def test_seqs_data_view_repr_default(dna_sd: new_alignment.SeqsData, seqid: str):
    seq = dna_sd.get_seq_str(seqid=seqid)
    got = dna_sd.get_view(seqid)
    expect = (
        f"SeqDataView(seqid='{seqid}', parent={seq}, slice_record={got.slice_record!r})"
    )
    assert repr(got) == expect


def test_seqs_data_view_repr_default_long(dna_alphabet):
    longseq = "CGATCGTAGTACGTGTCAAGTCTGAC"
    trunc = f"{longseq[:10]}...{longseq[-5:]}"

    d = {"long": longseq}
    sd = new_alignment.SeqsData(data=d, alphabet=dna_alphabet)
    got = sd.get_view("long")
    expect = (
        f"SeqDataView(seqid='long', parent={trunc}, slice_record={got.slice_record!r})"
    )
    assert repr(got) == expect


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
@pytest.mark.parametrize("sliced", (False, True))
@pytest.mark.parametrize("step", (1, -1))
def test_seqs_data_view_copy(dna_alphabet, seqid, sliced, step):
    seq1 = "ATGTTCTC"
    seq2 = "ATGAACTCATT"
    start, stop = 2, 6
    data = {"seq1": seq1, "seq2": seq2}

    sd = new_alignment.SeqsData(data=data, alphabet=dna_alphabet)
    sdv = sd.get_view(seqid)
    sliced_sdv = sdv[start:stop:step]
    copied_sdv = sliced_sdv.copy(sliced=sliced)

    assert copied_sdv.str_value == data[seqid][start:stop:step]


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
def test_seqs_data_get_seq_view(str_seqs_dict, dna_alphabet, seqid):
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=dna_alphabet)
    seq = str_seqs_dict[seqid]
    parent_len = len(seq)
    got = sd.get_view(seqid)
    assert isinstance(got, new_alignment.SeqDataView)
    assert got.parent == sd
    assert got.seqid == seqid
    assert got.parent_len == parent_len


@pytest.mark.parametrize("seq", ("seq1", "seq2"))
@pytest.mark.parametrize("start", (None, -1, 0, 1, 4))
@pytest.mark.parametrize("stop", (None, -1, 0, 1, 4))
def test_seqs_data_get_seq_str(str_seqs_dict, dna_alphabet, seq, start, stop):
    # slicing should be tested in test_get_seq_array
    expect = str_seqs_dict[seq][start:stop]
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=dna_alphabet)
    got = sd.get_seq_str(seqid=seq, start=start, stop=stop)
    assert expect == got


def test_seqs_data_get_seq_str_empty(dna_sd: new_alignment.SeqsData):
    with pytest.raises(TypeError):
        dna_sd.get_seq_str()


def test_seqs_data_names(str_seqs_dict, dna_alphabet):
    expect = str_seqs_dict.keys()
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=dna_alphabet)
    got = sd.names
    # returns iterator
    assert list(got) == list(expect)


def test_seqs_data_seq_lengths(str_seqs_dict, arr_seqs_dict, dna_alphabet):
    expect = {k: len(v) for k, v in str_seqs_dict.items()}
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=dna_alphabet)
    got = sd.seq_lengths()
    assert got == expect

    expect = {k: len(v) for k, v in arr_seqs_dict.items()}
    sd = new_alignment.SeqsData(data=arr_seqs_dict, alphabet=dna_alphabet)
    got = sd.seq_lengths()
    assert got == expect


def test_seqs_data_get_seq_array(str_seqs_dict, dna_alphabet):
    # seq1
    expect = numpy.array([2, 1, 3, 0], dtype="uint8")
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=dna_alphabet)
    got = sd.get_seq_array(seqid="seq1")
    assert numpy.array_equal(got, expect)

    # seq2
    expect = numpy.array([3, 0, 0, 0, 3, 1, 2], dtype="uint8")
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=dna_alphabet)
    got = sd.get_seq_array(seqid="seq2")
    assert numpy.array_equal(got, expect)


def test_seqs_data_get_seq_bytes(dna_sd: new_alignment.SeqsData):
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
    assert got.parent == dna_sd
    assert got.seqid == seq


@pytest.mark.parametrize("idx", (0, 1))
def test_seqs_data_getitem_int(str_seqs_dict, dna_alphabet, idx):
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=dna_alphabet)
    got = sd[idx]
    assert got.parent == sd
    assert got.seqid == list(str_seqs_dict)[idx]


@pytest.mark.parametrize("seqid", ("seq1", "seq2"))
def test_seqs_data_getitem_str(
    dna_sd,
    seqid,
):
    got = dna_sd[seqid]
    assert got.parent == dna_sd
    assert got.seqid == seqid


@pytest.mark.parametrize("make_seq", (True, False))
def test_seqs_data_getitem_raises(dna_sd, make_seq):
    ms = new_moltype.get_moltype("dna").make_seq if make_seq else None
    dna_sd.make_seq = ms
    invalid_index = ["this", "shouldn't", "work"]
    with pytest.raises(NotImplementedError):
        _ = dna_sd[invalid_index]


def test_seqs_data_subset(dna_sd):
    got = dna_sd.subset("seq1")
    assert isinstance(got, new_alignment.SeqsData)
    assert got.names == ["seq1"]
    assert got.get_seq_str(seqid="seq1") == dna_sd.get_seq_str(seqid="seq1")

    got = dna_sd.subset(["seq1", "seq2"])
    assert isinstance(got, new_alignment.SeqsData)
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


def test_seqs_data_to_alphabet():
    ASCII = new_moltype.ASCII.alphabet
    DNA = new_moltype.DNA.degen_gapped_alphabet
    RNA = new_moltype.RNA.degen_gapped_alphabet
    seqs = new_alignment.SeqsData(
        data={"a": "AAA", "b": "TTT", "c": "CCC"}, alphabet=ASCII
    )
    dna_seqs = seqs.to_alphabet(DNA)

    assert dna_seqs.get_seq_str(seqid="b") == "TTT"
    assert dna_seqs.alphabet == DNA

    rna_seqs = dna_seqs.to_alphabet(RNA)
    assert rna_seqs.get_seq_str(seqid="b") == "UUU"
    assert rna_seqs.alphabet == RNA


def test_seqs_data_to_alphabet_invalid():
    ASCII = new_moltype.ASCII.alphabet
    DNA = new_moltype.DNA.degen_gapped_alphabet
    seqs = new_alignment.SeqsData(
        data={"a": "AAA", "b": "TTT", "c": "LLL"}, alphabet=ASCII
    )
    with pytest.raises(ValueError):
        _ = seqs.to_alphabet(DNA)


@pytest.mark.parametrize("reverse", (False, True))
def test_seqs_data_round_trip(reverse, dna_alphabet):
    seqs_data = new_alignment.SeqsData(
        data={"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}, alphabet=dna_alphabet
    )
    seqs_data = seqs_data.reverse() if reverse else seqs_data

    rd = seqs_data.to_rich_dict()
    got = deserialise_object(rd)
    assert isinstance(got, new_alignment.SeqsData)
    assert got.to_rich_dict() == seqs_data.to_rich_dict()
    expect = "ACGG"[::-1] if reverse else "ACGG"
    assert str(got.get_view("seq1")) == expect


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
def test_seq_data_view_slice_returns_self(seq1: str, index: slice, dna_alphabet):
    sdv = new_alignment.SeqDataView(
        parent=seq1,
        seqid="seq1",
        alphabet=dna_alphabet,
        parent_len=len(seq1),
        slice_record=new_sequence.SliceRecord(parent_len=len(seq1)),
    )
    got = sdv[index]
    assert isinstance(got, new_alignment.SeqDataView)


# SeqDataView tests for value properties
@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seq_data_view_value(str_seqs_dict: dict, dna_alphabet, start, stop, step):
    seq = "seq2"
    expect = str_seqs_dict[seq][start:stop:step]
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=dna_alphabet)
    # Get SeqDataView on seq
    sdv = sd.get_view(seq)
    sdv2 = sdv[start:stop:step]
    got = sdv2.str_value
    assert got == expect


@pytest.mark.parametrize("rev", (False, True))
def test_seq_data_view_to_rich_dict(rev):
    data = {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}
    alpha = new_moltype.DNA.degen_gapped_alphabet
    sd = new_alignment.SeqsData(data=data, alphabet=alpha)
    sdv = sd.get_view("seq1")
    sdv = sdv[::-1] if rev else sdv
    got = sdv.to_rich_dict()
    expect = {
        "init_args": {
            "seqid": "seq1",
            "parent": sdv.str_value,
            "alphabet": alpha.to_rich_dict(),
            "slice_record": sdv.slice_record.to_rich_dict(),
        },
        "type": get_object_provenance(sdv),
        "version": __version__,
    }

    assert got == expect


@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seqs_data_array_value(arr_seqs_dict: dict, dna_alphabet, start, stop, step):
    seq = "seq2"
    expect = arr_seqs_dict[seq][start:stop:step]
    sd = new_alignment.SeqsData(data=arr_seqs_dict, alphabet=dna_alphabet)
    # Get SeqDataView on seq
    sdv = sd.get_view(seq)
    got = sdv.array_value[start:stop:step]
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seqs_data_bytes_value(str_seqs_dict: dict, dna_alphabet, start, stop, step):
    seq = "seq2"
    expect = str_seqs_dict[seq][start:stop:step]
    expect = expect.encode("utf8")
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=dna_alphabet)
    # Get SeqDataView on seq
    sdv = sd.get_view(seq)
    got = sdv.bytes_value[start:stop:step]
    assert expect == got


# SeqDataView tests for special methods that access "value" properties
def test_array(sdv_s2: new_alignment.SeqDataView):
    expect = sdv_s2.array_value
    got = numpy.array(sdv_s2)
    assert numpy.array_equal(expect, got)


def test_bytes(sdv_s2: new_alignment.SeqDataView):
    expect = sdv_s2.bytes_value
    got = bytes(sdv_s2)
    assert expect == got


@pytest.mark.parametrize("moltype", ("dna", "rna", "protein", "protein_with_stop"))
def test_make_unaligned_seqs_dict(moltype):
    """test SequenceCollection constructor utility function"""
    data = {"a": "AGGCCC", "b": "AGAAAA"}
    got = new_alignment.make_unaligned_seqs(data, moltype=moltype)
    assert isinstance(got, new_alignment.SequenceCollection)
    assert got._seqs_data.make_seq == new_moltype.get_moltype(moltype).make_seq

    # should also work if seqs are arrays
    data = {"a": numpy.array([0, 2, 1, 3, 2, 3]), "b": numpy.array([0, 2, 1, 3, 2, 1])}
    got = new_alignment.make_unaligned_seqs(data, moltype=moltype)
    assert isinstance(got, new_alignment.SequenceCollection)
    assert got._seqs_data.make_seq == new_moltype.get_moltype(moltype).make_seq


@pytest.mark.parametrize("moltype", ("dna", "rna", "protein", "protein_with_stop"))
@pytest.mark.parametrize("cls", (list, set))
def test_make_unaligned_seqs_collection(moltype, cls):
    """test SequenceCollection constructor utility function"""
    data = cls(["AGGCCC", "AGAAAA"])
    got = new_alignment.make_unaligned_seqs(data, moltype=moltype)
    assert isinstance(got, new_alignment.SequenceCollection)
    assert got._seqs_data.make_seq == new_moltype.get_moltype(moltype).make_seq
    assert got.num_seqs == 2

    # should also work if seqs are arrays
    data = [numpy.array([0, 2, 1, 3, 2, 3]), numpy.array([0, 2, 1, 3, 2, 1])]
    got = new_alignment.make_unaligned_seqs(data, moltype=moltype)
    assert isinstance(got, new_alignment.SequenceCollection)
    assert got._seqs_data.make_seq == new_moltype.get_moltype(moltype).make_seq

    # also when seqs are list(name: seq) pairs
    data = [["seq1", "AGGCCC"], ["seq2", "AGAAAA"]]
    got = new_alignment.make_unaligned_seqs(data, moltype=moltype)
    assert isinstance(got, new_alignment.SequenceCollection)

    # or tuples
    data = [("seq1", "AGGCCC"), ("seq2", "AGAAAA")]
    got = new_alignment.make_unaligned_seqs(data, moltype=moltype)
    assert isinstance(got, new_alignment.SequenceCollection)


def test_make_unaligned_seqs_label_to_name():
    """test SequenceCollection constructor utility function"""
    data = {"a": "AGGCCC", "b": "AGAAAA"}

    def f(x):
        return x.upper()

    got = new_alignment.make_unaligned_seqs(data, moltype="dna", label_to_name=f)
    assert list(got.names) == ["A", "B"]
    assert got._seqs_data.make_seq == new_moltype.get_moltype("dna").make_seq


def test_make_unaligned_seqs_raises():
    data = "AGTCCTGA"
    with pytest.raises(NotImplementedError):
        new_alignment.make_unaligned_seqs(data, moltype="dna")


def test_make_unaligned_seqs_incompatible_moltype(dna_sd):
    with pytest.raises(ValueError):
        _ = new_alignment.make_unaligned_seqs(dna_sd, moltype="rna")


def test_make_unaligned_seqs_no_seqs():
    data = {}
    with pytest.raises(ValueError):
        new_alignment.make_unaligned_seqs(data, moltype="dna")


@pytest.mark.parametrize(
    "collection_maker",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_sequence_collection_names_is_list(collection_maker):
    """expected to be a list"""
    seqs = collection_maker({"a": b"AAAAA", "b": b"TTTTT"}, moltype="dna")
    assert isinstance(seqs.names, list)


@pytest.mark.parametrize(
    "collection_maker",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_sequence_collection_names(collection_maker):
    seqs = {"seq1": "AAAAAA", "seq2": "TTTT--", "seq3": "ATTCCC"}
    seq_coll = collection_maker(seqs, moltype="dna")
    assert seq_coll.names == ["seq1", "seq2", "seq3"]
    seq_coll.names = ["seq2", "seq3", "seq1"]
    assert seq_coll.names == ["seq2", "seq3", "seq1"]
    seq_coll.names = ["seq1", "seq2"]
    assert seq_coll.names == ["seq1", "seq2"]
    with pytest.raises(ValueError):
        seq_coll.names = ["seq1", "seq2", "seq3", "seq4"]


def test_make_unaligned_seqs_list_str():
    """SequenceCollection init from list of sequences should use indices as keys"""
    seqs = ["TTTTT", "CCCCC", "GGGGG"]
    a = new_alignment.make_unaligned_seqs(seqs, moltype="dna")
    assert a.seqs["seq_0"] == "TTTTT"
    assert a.seqs["seq_1"] == "CCCCC"
    assert a.seqs["seq_2"] == "GGGGG"
    assert a.names == ["seq_0", "seq_1", "seq_2"]


def test_make_unaligned_seqs_pairs():
    """SequenceCollection init from list of (key,val) pairs should work correctly"""
    seqs = [["a", "AAA"], ["t", "TTT"], ["c", "CCC"]]
    a = new_alignment.make_unaligned_seqs(seqs, moltype="dna")
    assert a.seqs["a"] == "AAA"
    assert a.seqs["t"] == "TTT"
    assert a.seqs["c"] == "CCC"
    assert a.names == ["a", "t", "c"]


def test_make_unaligned_seqs_list_sequences():
    """correctly construct from list of sequences of length 2"""
    seq1 = new_moltype.DNA.make_seq(seq="AC", name="seq1")
    seq2 = new_moltype.DNA.make_seq(seq="AC", name="seq2")
    coll = new_alignment.make_unaligned_seqs([seq1, seq2], moltype="dna")
    assert isinstance(coll, new_alignment.SequenceCollection)


def test_sequence_collection_init_ordered(ordered1, ordered2):
    """SequenceCollection should iterate over seqs correctly even if ordered"""
    first = ordered1
    sec = ordered2

    assert first.names == ["a", "c"]
    assert sec.names == ["c", "a"]


def test_sequence_collection_info_source():
    """info.source exists if load seqs given a filename"""
    path = pathlib.Path("data/brca1.fasta")
    seqs = load_unaligned_seqs(path, moltype="dna", new_type=True)
    assert seqs.info.source == str(path)


@pytest.mark.parametrize(
    "collection_maker",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_sequence_collection_init_ambig(collection_maker):
    """SequenceCollection and Alignment should tolerate ambiguous chars"""
    _ = collection_maker({"s0": "AAA", "s1": "CCC"}, moltype="dna")
    _ = collection_maker({"s0": "ANS", "s1": "CWC"}, moltype="dna")
    _ = collection_maker({"s0": "A-A", "s1": "CC-"}, moltype="dna")
    _ = collection_maker({"s0": "A?A", "s1": "CC-"}, moltype="dna")


@pytest.mark.parametrize(
    "collection_maker",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_sequence_collection_iter_seqs_ragged_padded(
    ragged_padded_dict, collection_maker
):
    """SequenceCollection.iter_seqs() method should support reordering of seqs"""
    coll = collection_maker(ragged_padded_dict, moltype="dna")
    seqs = list(coll.iter_seqs())
    assert list(map(str, seqs)) == ["AAAAAA", "AAA---", "AAAA--"]
    seqs = list(coll.iter_seqs(seq_order=["b", "a", "a"]))
    assert list(map(str, seqs)) == ["AAA---", "AAAAAA", "AAAAAA"]
    assert str(seqs[1]) == str(seqs[2])
    assert str(seqs[0]) == str(coll.seqs["b"])


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
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    assert (
        repr(seqs)
        == "2x (ENSMUSG00000039616[GCCCTTCAAA...], ENSMUSG00000056468[GCCAGGGGGA...]) dna seqcollection"
    )

    data = {
        "ENSMUSG00000039616": "GCCCTTCAAATTT",
        "ENSMUSG00000056468": "GCCAGGGGGAAAAGGGAGAA",
    }
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    assert (
        repr(seqs)
        == "2x (ENSMUSG00000039616[GCCCTTCAAA...], ENSMUSG00000056468[GCCAGGGGGA...]) dna seqcollection"
    )

    data = {
        "a": "TCGAT",
    }
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    assert repr(seqs) == "1x (a[TCGAT]) dna seqcollection"

    data = {
        "a": "TCGAT" * 2,
    }
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    assert repr(seqs) == "1x (a[TCGATTCGAT]) dna seqcollection"


def test_sequence_collection_set_wrap_affects_repr_html():
    """the wrap argument affects the number of columns"""
    # indirectly tested via counting number of occurrences of 'class="label"'
    seqs = new_alignment.make_unaligned_seqs(
        {"a": "AAAAA", "b": "AAA--"}, moltype="dna"
    )
    orig = seqs._repr_html_()
    seqs.set_repr_policy(wrap=3)  # break alignment into 2
    got = seqs._repr_html_()
    token = 'class="label"'
    assert got.count(token) == 2 * orig.count(token)

    # using environment variable
    env_name = "COGENT3_ALIGNMENT_REPR_POLICY"
    os.environ[env_name] = "wrap=2"
    seqs = new_alignment.make_unaligned_seqs(
        {"a": "AAAAA", "b": "AAA--"}, moltype="dna"
    )
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
    seqs.set_repr_policy(num_seqs=seqs.num_seqs)
    got = seqs._repr_html_()
    assert got.count("</tr>") == seqs.num_seqs + 1


def test_sequence_collection_repr_html_correct_num_seqs(seqs):
    # tests _repr_html_ displays correct number of sequences
    got = seqs._repr_html_()
    seq_lens = numpy.array([len(seqs.seqs[name]) for name in seqs.names])
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


@pytest.mark.parametrize(
    "collection_maker",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
@pytest.mark.parametrize(
    "sample", [["a"], ["b"], ["c"], ["a", "b"], ["a", "c"], ["b", "c"]]
)
def test_sequence_collection_take_seqs(ragged_padded_dict, collection_maker, sample):
    """take_seqs should return new SequenceCollection/Alignment with selected seqs."""
    orig = collection_maker(ragged_padded_dict, moltype="dna")
    subset = orig.take_seqs(sample)
    assert subset.names == sample
    assert subset.num_seqs == len(sample)
    # should be able to negate
    neg_subset = orig.take_seqs(sample, negate=True)
    assert neg_subset.names == [name for name in orig.names if name not in sample]


@pytest.mark.parametrize(
    "collection_maker",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_take_seqs_str(ragged_padded_dict, collection_maker):
    """string arg to SequenceCollection take_seqs should work."""
    orig = collection_maker(ragged_padded_dict, moltype="dna")
    subset = orig.take_seqs("a")
    assert subset.names == ["a"]
    assert subset.num_seqs == 1

    # should be able to negate
    neg_subset = orig.take_seqs("a", negate=True)
    assert neg_subset.names == ["b", "c"]
    assert neg_subset.num_seqs == 2


@pytest.mark.parametrize(
    "collection_maker",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_take_seqs_info(ragged_padded_dict, collection_maker):
    """take_seqs should preserve info attribute"""
    orig = collection_maker(
        ragged_padded_dict,
        moltype="dna",
        info={"key": "value"},
    )
    subset = orig.take_seqs(list("ab"))
    assert set(subset.info) == set(orig.info)


@pytest.mark.parametrize(
    "collection_maker",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_take_seqs_moltype(ragged_padded_dict, collection_maker):
    """take_seqs should preserve the MolType"""
    orig = collection_maker(ragged_padded_dict, moltype="dna")
    subset = orig.take_seqs(list("ab"))
    assert set(subset.moltype) == set(orig.moltype)


@pytest.mark.parametrize(
    "collection_maker",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_take_seqs_empty_names(
    ragged_padded_dict, collection_maker
):
    """take_seqs should raise ValueError if no seqs are selected."""
    orig = new_alignment.make_unaligned_seqs(ragged_padded_dict, moltype="dna")
    with pytest.raises(ValueError):
        _ = orig.take_seqs([])


@pytest.mark.parametrize(
    "collection_maker",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_take_seqs_copy_annotations(gff_db, collection_maker):
    data = {"test_seq": "ACGT--", "test_seq2": "CGTTTA"}
    seq_coll = collection_maker(data, moltype="dna", annotation_db=gff_db)
    seq_coll.annotation_db = gff_db
    # if copy annotations is true, only the annotations for the selected seqs should be copied
    just_seq2 = seq_coll.take_seqs(names="test_seq2", copy_annotations=True)
    assert len(just_seq2.annotation_db) == 1

    # if copy annotations is false, the annotation db should be the same as original
    just_seq1 = seq_coll.take_seqs(names="test_seq", copy_annotations=False)
    assert just_seq1.annotation_db == gff_db


def test_sequence_collection_num_seqs():
    """SequenceCollection.num_seqs should count seqs."""
    data = {"seq1": "ACGU", "seq2": "CGUA", "seq3": "CCGU"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="rna")
    assert seqs.num_seqs == 3


@pytest.mark.parametrize("index", [(0, "seq1"), (1, "seq2"), (2, "seq3")])
def test_sequence_collection_getitem(seqs, index):
    got1 = seqs.seqs[index[0]]
    got2 = seqs.seqs[index[1]]

    assert isinstance(got1, new_sequence.Sequence)
    assert isinstance(got2, new_sequence.Sequence)
    assert got1 == got2


# Functions for testing SequenceCollection methods that accept a function as an argument
def is_long(x):
    return len(x) > 10


def is_med(x):
    return len(str(x).replace("-", "")) > 3


def is_any(x):
    return len(x) > 0


def test_sequence_collection_get_seq_names_if(ragged_padded):
    """SequenceCollection.get_seq_names_if should return names of seqs where f(row) is True"""

    assert ragged_padded.get_seq_names_if(is_long) == []
    assert ragged_padded.get_seq_names_if(is_med) == ["a", "c"]
    # return order should reflect names when updated
    ragged_padded.names = ["b", "c", "a"]
    assert ragged_padded.get_seq_names_if(is_med) == ["c", "a"]
    assert ragged_padded.get_seq_names_if(is_med, negate=True) == ["b"]
    assert ragged_padded.get_seq_names_if(is_any) == ["b", "c", "a"]
    assert ragged_padded.get_seq_names_if(is_any, negate=True) == []


def test_sequence_collection_take_seqs_if(ragged_padded):
    """SequenceCollection take_seqs_if should return seqs where f(seq) is True"""

    with pytest.raises(ValueError):
        ragged_padded.take_seqs_if(is_long)
    with pytest.raises(ValueError):
        ragged_padded.take_seqs_if(is_any, negate=True)
    got = ragged_padded.take_seqs_if(is_any)
    assert got.names == ragged_padded.names
    assert got.num_seqs == ragged_padded.num_seqs
    assert isinstance(
        ragged_padded.take_seqs_if(is_med), new_alignment.SequenceCollection
    )

    got = ragged_padded.take_seqs_if(is_med, negate=True)
    assert got.names == ["b"]
    assert got.num_seqs == 1


def test_sequence_collection_to_dict():
    """SequenceCollection.to_dict should return dict of strings (not obj)"""
    data = {"seq1": "GATTTT", "seq2": "GATC??"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    assert seqs.to_dict() == data
    for i in list(seqs.to_dict().values()):
        assert isinstance(i, str)


def test_sequence_collection_get_seq():
    """SequenceCollection.get_seq should return specified seq"""
    data = {"seq1": "GATTTT", "seq2": "GATC??"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    assert seqs.get_seq("seq1") == "GATTTT"
    with pytest.raises(KeyError):
        seqs.get_seq("seqx")


def test_sequence_collection_degap():
    """SequenceCollection.degap should strip gaps from each seq"""
    data = {"s1": "ATGRY?", "s2": "T-AG??"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = seqs.degap().to_dict()
    expect = {"s1": "ATGRY", "s2": "TAG"}
    assert got == expect


def test_sequence_collection_degap_info():
    """.degap should preserve info attributes"""
    data = {"s1": "ATGRY?", "s2": "T-AG??"}
    aln = new_alignment.make_unaligned_seqs(data, moltype="dna")
    aln.info.path = "blah"
    got = aln.degap()
    assert got.info.path == "blah"


def test_sequence_collection_to_fasta():
    """SequenceCollection should return correct FASTA string"""
    data = {"seq_0": "AAA", "seq_1": "CCC"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    assert seqs.to_fasta() == ">seq_0\nAAA\n>seq_1\nCCC\n"
    assert seqs.to_fasta(block_size=2) == ">seq_0\nAA\nA\n>seq_1\nCC\nC\n"

    data = {"seq_0": "GCATGCAT", "seq_1": "TCAGACGT"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    assert seqs.to_fasta() == ">seq_0\nGCATGCAT\n>seq_1\nTCAGACGT\n"
    assert seqs.to_fasta(block_size=4) == ">seq_0\nGCAT\nGCAT\n>seq_1\nTCAG\nACGT\n"
    assert seqs.to_fasta(block_size=3) == ">seq_0\nGCA\nTGC\nAT\n>seq_1\nTCA\nGAC\nGT\n"


def test_sequence_collection_is_ragged(ragged, ragged_padded):
    """SequenceCollection is_ragged should return true if ragged alignment"""
    assert ragged.is_ragged()
    assert not ragged_padded.is_ragged()


def test_sequence_collection_ragged(ragged):
    """SequenceCollection seqs should work on ragged alignment"""
    ragged.names = "bac"
    assert [ragged.seqs[name] for name in ragged.names] == [
        "AAA",
        "AAAAAA",
        "AAAA",
    ]


def test_sequence_collection_to_phylip():
    """SequenceCollection should return PHYLIP string format correctly"""
    data = {
        "seq_0": "ACDEFGHIKLMNPQRSTUVWY-",
        "seq_1": "ACDEFGHIKLMNPQRSUUVWF-",
        "seq_2": "ACDEFGHIKLMNPERSKUVWC-",
        "seq_3": "ACNEFGHIKLMNPQRS-UVWP-",
    }
    align_norm = new_alignment.make_unaligned_seqs(data, moltype="protein")

    assert (
        align_norm.to_phylip()
        == """4  22\nseq_0     ACDEFGHIKLMNPQRSTUVWY-\nseq_1     ACDEFGHIKLMNPQRSUUVWF-\nseq_2     ACDEFGHIKLMNPERSKUVWC-\nseq_3     ACNEFGHIKLMNPQRS-UVWP-\n"""
    )


def test_sequence_collection_to_phylip_ragged():
    """SequenceCollection should refuse to convert ragged seqs to phylip"""
    data = {"seq1": "KGA-", "seq2": "KGA"}
    align_rag = new_alignment.make_unaligned_seqs(data, moltype="protein")

    with pytest.raises(ValueError):
        align_rag.to_phylip()


@pytest.mark.parametrize(
    "gc,seqs",
    (
        (1, ("TCCTGA", "GATTT?")),
        (1, ("ACGTAA---", "ACGAC----", "ACGCAATGA")),
        (2, ("GATTTT", "TCCAGG")),
    ),
)
def test_sequence_collection_has_terminal_stop_true(gc, seqs):
    data = {f"s{i}": s for i, s in enumerate(seqs)}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    assert seq_coll.has_terminal_stop(gc=gc)


@pytest.mark.parametrize(
    "gc,seqs",
    ((1, ("TCCTCA", "GATTTT")), (2, ("GATTTT", "TCCCGG")), (1, ("CCTCA", "ATTTT"))),
)
def test_sequence_collection_has_terminal_stop_false(gc, seqs):
    data = {f"s{i}": s for i, s in enumerate(seqs)}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    assert not seq_coll.has_terminal_stop(gc=gc)


def test_sequence_collection_has_terminal_stop_strict():
    data = {f"s{i}": s for i, s in enumerate(("CCTCA", "ATTTT"))}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    with pytest.raises(new_alphabet.AlphabetError):
        seq_coll.has_terminal_stop(gc=1, strict=True)


def test_sequence_collection_get_translation_trim_stop():
    data = {"seq1": "GATTCCTAG", "seq2": "GATTCCTCC"}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = seq_coll.get_translation(trim_stop=True)
    expect = {"seq1": "DS", "seq2": "DSS"}
    assert got.to_dict() == expect


def test_sequence_collection_get_translation_raises():
    """should raise error if self.moltype is not a nucleic acid"""
    data = {"seq1": "PAR", "seq2": "PQR"}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="protein")
    with pytest.raises(new_alphabet.AlphabetError):
        _ = seq_coll.get_translation(trim_stop=True)


@pytest.mark.parametrize(
    "seqs",
    (
        {"seq1": "GATTTT", "seq2": "GATC??"},
        {"seq1": "GAT---", "seq2": "GATCTT"},
        {"seq1": "GAT-T-", "seq2": "GATCTT"},
        {"seq1": "GATTTT", "seq2": "?GATCT"},
    ),
)
def test_sequence_collection_get_translation(seqs):
    """SequenceCollection.get_translation translates each seq"""
    seq_coll = new_alignment.make_unaligned_seqs(seqs, moltype="dna")
    got = seq_coll.get_translation(incomplete_ok=True)
    assert got.num_seqs == 2
    assert got.moltype == new_moltype.PROTEIN


def test_sequence_collection_get_translation_with_stop():
    data = {"seq1": "?GATCT", "seq2": "GATTAG"}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = seq_coll.get_translation(
        incomplete_ok=True, include_stop=True, trim_stop=False
    )
    assert got.to_dict() == {"seq1": "XS", "seq2": "D*"}
    assert got.moltype == new_moltype.PROTEIN_WITH_STOP


def test_sequence_collection_get_translation_non_div_3():
    data = {"seq1": "?GATCTA", "seq2": "GATTAGGG"}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = seq_coll.get_translation(incomplete_ok=True, include_stop=True)
    assert got.to_dict() == {"seq1": "XS", "seq2": "D*"}


@pytest.mark.parametrize(
    "data", ({"seq1": "GATTTT", "seq2": "GATC??"}, {"seq1": "GAT---", "seq2": "?GATCT"})
)
def test_sequence_collection_get_translation_error(data):
    """SequenceCollection.get_translation translates each seq"""
    # check for a failure when no moltype specified
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    with pytest.raises(TypeError):
        seq_coll.get_translation()


@pytest.mark.parametrize(
    "data",
    (
        {"seq1": "GATTTT", "seq2": "GATCTT"},
        {"seq1": "GAT---", "seq2": "GGATCT"},
        {"seq1": "GATTTT", "seq2": "?GATCT"},
        {"seq1": "GATTTT", "seq2": "GATC??"},
        {"seq1": "GAT-T-", "seq2": "GATCTT"},
    ),
)
def test_sequence_collection_get_translation_info(data):
    """SequenceCollection.get_translation preserves info attribute"""
    seq_coll = new_alignment.make_unaligned_seqs(
        data, moltype="dna", info={"key": "value"}
    )
    got = seq_coll.get_translation(incomplete_ok=True)
    assert got.info["key"] == "value"


def test_sequence_collection_get_translation_incomplete():
    """get translation works on incomplete codons"""
    data = {"seq1": "GATN--", "seq2": "?GATCT"}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = seq_coll.get_translation(incomplete_ok=True)
    assert got.to_dict() == {"seq1": "DX", "seq2": "XS"}
    with pytest.raises(new_alphabet.AlphabetError):
        _ = seq_coll.get_translation(incomplete_ok=False)


@pytest.mark.parametrize(
    "gc,seqs",
    (
        (1, ("--AT-CTGA", "GATAAATT?")),
        (1, ("ACGTGA---", "ACGAC----", "ACGCAATGA")),
        (1, ("CCTCA-", "ATTTTA")),
        (2, ("GATTTT", "TCCAGG")),
    ),
)
def test_sequence_collection_trim_stop_codons(gc, seqs):
    data = {f"s{i}": s for i, s in enumerate(seqs)}

    expect = {}
    for k, v in data.items():
        if "-" in v:
            v = re.sub("(TGA|AGG)(?=[-]*$)", "---", v)
        else:
            v = re.sub("(TGA|AGG)", "", v)
        expect[k] = v
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = seqs.trim_stop_codons(gc=gc).to_dict()

    assert got == expect


@pytest.mark.parametrize(
    "gc,seqs",
    ((1, ("T-CTGC", "GATAA?")), (2, ("GATTTT", "TCCCGG")), (1, ("CCTGC", "GATAA"))),
)
def test_sequence_collection_trim_stop_codons_no_stop(gc, seqs):
    data = {f"s{i}": s for i, s in enumerate(seqs)}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = seqs.trim_stop_codons(gc=gc)
    assert got is seqs


@pytest.mark.parametrize(
    "data", ({"s1": "CCTCA", "s2": "ATTTT"}, {"s1": "CCTCA-", "s2": "ATTTTA"})
)
def test_sequence_collection_trim_stop_codons_strict(data):
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    with pytest.raises(new_alphabet.AlphabetError):
        seqs.trim_stop_codons(gc=1, strict=True)


def test_sequence_collection_trim_stop_codons_info():
    """trim_stop_codons should preserve info attribute"""
    data = {"seq1": "ACGTAA", "seq2": "ACGACG", "seq3": "ACGCGT"}
    seq_coll = new_alignment.make_unaligned_seqs(
        data,
        moltype="dna",
        info={"key": "value"},
    )
    seq_coll = seq_coll.trim_stop_codons()
    assert seq_coll.info["key"] == "value"


def test_sequence_collection_trim_stop_codons_annotation_db(gff_db):
    """trim_stop_codons should preserve info attribute"""
    data = {"seq_1": "ACGTAA", "seq_2": "ACGACG", "seq_3": "ACGCGT"}
    seq_coll = new_alignment.make_unaligned_seqs(
        data,
        moltype="dna",
        info={"key": "value"},
    )
    seq_coll.annotation_db = gff_db
    trimmed = seq_coll.trim_stop_codons()
    assert trimmed.annotation_db == seq_coll.annotation_db


@pytest.fixture(scope="function")
def dna_seqs_with_dupes():
    data = {
        "a": "ACGG",
        "b": "ACGG",  # strict identical
        "c": "ACGN",  # non-strict matches above
        "d": "ACGT",
        "e": "ACGT",
        "k": "ACGT",  # strict identical
        "f": "RAAA",
        "g": "YAAA",  # non-strict identical
    }

    return new_alignment.make_unaligned_seqs(data, moltype="dna")


def test_sequence_collection_get_identical_sets_dna(dna_seqs_with_dupes):
    """correctly identify sets of identical sequences"""
    got = dna_seqs_with_dupes.get_identical_sets(
        mask_degen=False
    )  # a strict comparison
    # convert to frozenset, so we can do a comparison robust to set order
    got = frozenset(frozenset(s) for s in got)
    expect = [{"a", "b"}, {"d", "e", "k"}]
    expect = frozenset(frozenset(s) for s in expect)
    assert got == expect


def test_sequence_collection_get_identical_sets_dna_mask_degen(dna_seqs_with_dupes):
    got = dna_seqs_with_dupes.get_identical_sets(mask_degen=True)
    got = frozenset(frozenset(s) for s in got)
    expect = [{"a", "b", "c"}, {"d", "e", "k"}, {"f", "g"}]
    expect = frozenset(frozenset(s) for s in expect)
    assert got == expect


@pytest.fixture(scope="function")
def protein_seqs_with_dupes():
    data = {
        "a": "ACGT",
        "b": "ACGT",  # strict identical
        "c": "ACGX",  # non-strict matches above
        "d": "TTTT",
        "e": "TTTT",
        "k": "TTTT",  # strict identical
        "f": "BAAA",
        "g": "ZAAA",  # non-strict identical
    }

    return new_alignment.make_unaligned_seqs(data, moltype="protein")


def test_sequence_collection_get_identical_sets_protein(protein_seqs_with_dupes):
    got = protein_seqs_with_dupes.get_identical_sets(
        mask_degen=False
    )  # a strict comparison
    # convert to frozenset, so we can do a comparison robust to set order
    got = frozenset(frozenset(s) for s in got)
    expect = [{"a", "b"}, {"d", "e", "k"}]
    expect = frozenset(frozenset(s) for s in expect)
    assert got == expect


def test_sequence_collection_get_identical_sets_protein_mask_degen(
    protein_seqs_with_dupes,
):
    got = protein_seqs_with_dupes.get_identical_sets(mask_degen=True)
    got = frozenset(frozenset(s) for s in got)
    expect = [{"a", "b", "c"}, {"d", "e", "k"}, {"f", "g"}]
    expect = frozenset(frozenset(s) for s in expect)
    assert got == expect


def test_sequence_collection_get_identical_sets_text():
    data = {
        "a": "ACGT",
        "b": "ACGT",  # strict identical
        "c": "ACGX",  # non-strict matches above
        "d": "TTTT",
        "e": "TTTT",
        "k": "TTTT",  # strict identical
        "f": "BAAA",
        "g": "ZAAA",  # non-strict identical
    }
    seqs = new_alignment.make_unaligned_seqs(data, moltype="text")
    got = seqs.get_identical_sets(mask_degen=False)
    # convert to frozenset, so we can do a comparison robust to set order
    got = frozenset(frozenset(s) for s in got)
    expect = [{"a", "b"}, {"d", "e", "k"}]
    expect = frozenset(frozenset(s) for s in expect)
    assert got == expect

    with catch_warnings():
        filterwarnings("ignore", category=UserWarning)
        got = seqs.get_identical_sets(mask_degen=True)

    got = frozenset(frozenset(s) for s in got)
    assert got == expect


def test_sequence_collection_get_identical_sets_ragged_raises(ragged):
    with pytest.raises(ValueError):
        _ = ragged.get_identical_sets()


def rev(seq):
    return seq[::-1]


@pytest.mark.parametrize("transform", [rev, lambda x: x, None])
def test_sequence_collection_get_similar(transform):
    data = {
        "a": "AAAAAAAAAA",  # target
        "b": "AAAAAAAAAA",  # identical
        "c": "AAAAAAAAAT",  # 90% identical
        "d": "AAAAAAAATT",  # 80% identical
        "e": "AATTTTTTTT",  # 20% identical
        "f": "TTTTTTTTTT",  # 0% identical
    }
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    target = seqs.get_seq("a")
    got = seqs.get_similar(target, min_similarity=0.8, transform=transform)
    assert got.names == ["a", "b", "c", "d"]
    got = seqs.get_similar(target, min_similarity=0.81, transform=transform)
    assert got.names == ["a", "b", "c"]
    got = seqs.get_similar(
        target, min_similarity=0.75, max_similarity=0.9, transform=transform
    )
    assert got.names == ["c", "d"]
    got = seqs.get_similar(
        target, min_similarity=0.75, max_similarity=0.89, transform=transform
    )
    assert got.names == ["d"]


# todo: kath: add test for Alignment when it can be constructed with Sequence objects
@pytest.mark.parametrize("coll_maker", [new_alignment.make_unaligned_seqs])
def test_sequence_collection_init_annotated_seqs(coll_maker):
    """correctly construct from list with annotated seq"""
    seq = new_moltype.DNA.make_seq(seq="GCCAGGGGGGAAAG-GGAGAA", name="seq1")
    seq.add_feature(biotype="exon", name="name", spans=[(4, 10)])
    coll = coll_maker({"seq1": seq}, moltype="dna")
    features = list(coll.get_features(biotype="exon"))
    assert len(features) == 1


def test_sequence_collection_get_motif_probs():
    data = {"a": "AC-", "b": "AC", "c": "AC"}
    aln = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = aln.get_motif_probs()
    expect = {"A": 0.5, "C": 0.5, "G": 0.0, "T": 0.0}
    assert got == expect

    # exclude unobserved
    data = {"a": "AC-", "b": "AC", "c": "AC"}
    aln = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = aln.get_motif_probs(exclude_unobserved=True)
    expect = {"A": 0.5, "C": 0.5}
    assert got == expect

    # allow gap
    data = {"a": "----", "b": "ACGT", "c": "ACGT"}
    aln = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = aln.get_motif_probs(allow_gap=True)
    expect = {"A": 2 / 12, "C": 2 / 12, "G": 2 / 12, "T": 2 / 12, "-": 4 / 12}
    assert got == expect

    # add pseudocounts
    data = {"a": "AC-", "b": "AC", "c": "AC"}
    aln = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = aln.get_motif_probs(pseudocount=1)
    expect = {"A": 4 / 10, "C": 4 / 10, "G": 1 / 10, "T": 1 / 10}
    assert got == expect

    # pseudocount and allow gap
    data = {"a": "AC-", "b": "AC", "c": "AC"}
    aln = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = aln.get_motif_probs(pseudocount=1, allow_gap=True)
    expect = {"A": 4 / 12, "C": 4 / 12, "G": 1 / 12, "T": 1 / 12, "-": 2 / 12}
    assert got == expect


def test_sequence_collection_get_motif_probs_protein():
    data = {"a": "MVSB", "b": "MVS", "c": "MVP"}
    aln = new_alignment.make_unaligned_seqs(data, moltype="protein")
    got = aln.get_motif_probs()
    alphabet = aln.moltype.alphabet
    expect = {k: 0 for k in alphabet}
    expect["M"] = 1 / 3
    expect["V"] = 1 / 3
    expect["S"] = 2 / 9
    expect["P"] = 1 / 9
    assert got == expect

    got = aln.get_motif_probs(exclude_unobserved=True)
    expect = {"M": 1 / 3, "V": 1 / 3, "S": 2 / 9, "P": 1 / 9}
    assert got == expect

    got = aln.get_motif_probs(include_ambiguity=True, exclude_unobserved=True)
    expect = {"M": 3 / 10, "V": 3 / 10, "S": 2 / 10, "P": 1 / 10, "B": 1 / 10}


def test_sequence_collection_counts_per_seq():
    """SequenceCollection.counts_per_seq handles motif length, allow_gaps etc.."""
    data = {"a": "AAAA??????", "b": "CCCGGG--NN", "c": "CCGGTTCCAA"}
    coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    mtype = coll.moltype
    got = coll.counts_per_seq()
    assert got["a", "A"] == 4
    assert len(got.motifs) == len(mtype.alphabet)
    got = coll.counts_per_seq(include_ambiguity=True, allow_gap=True)
    # N, -, ? are the additional states
    assert len(got.motifs) == 7
    expect = {"-": 2, "?": 0, "A": 0, "C": 3, "G": 3, "N": 2, "T": 0}
    b = got["b"].to_dict()
    for k in expect:
        assert b[k] == expect[k]

    got = coll.counts_per_seq(motif_length=2)
    assert len(got.motifs), len(mtype.alphabet) ** 2
    assert got["a", "AA"] == 2
    assert got["b", "GG"] == 1
    got = coll.counts_per_seq(exclude_unobserved=True)
    expect = {"C": 4, "G": 2, "T": 2, "A": 2}
    c = got["c"].to_dict()
    for k in expect:
        assert c[k] == expect[k]


@pytest.mark.parametrize(
    "coll_maker", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_sequence_collection_probs_per_seq(coll_maker):
    data = {"seq1": "AA??", "seq2": "CG-N", "seq3": "CGAA"}
    coll = coll_maker(data, moltype="dna")
    got = coll.probs_per_seq()
    assert got["seq1", "A"] == 1.0
    assert got["seq2", "C"] == 0.5
    assert got["seq2", "G"] == 0.5
    assert got["seq3", "C"] == 0.25
    assert got["seq3", "G"] == 0.25
    assert got["seq3", "A"] == 0.5

    got = coll.probs_per_seq(allow_gap=True)
    assert got["seq2", "-"] == 1 / 3

    got = coll.probs_per_seq(include_ambiguity=True, allow_gap=True)
    assert got["seq2", "N"] == 0.25


def test_sequence_collection_counts():
    """SequenceCollection.counts handles motif length, allow_gaps etc.."""
    data = {"a": "AAAA??????", "b": "CCCGGG--NN"}
    coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = coll.counts()
    expect = dict(A=4, C=3, G=3)
    assert all(got[k] == v for k, v in expect.items())

    got = coll.counts(motif_length=2)
    expect = dict(AA=2, CC=1, CG=1, GG=1)
    assert all(got[k] == v for k, v in expect.items())

    got = coll.counts(motif_length=2, allow_gap=True)
    expect["--"] = 1
    assert all(got[k] == v for k, v in expect.items())

    got = coll.counts(motif_length=2, include_ambiguity=True, allow_gap=True)
    expect = dict(AA=2, CC=1, CG=1, GG=1, NN=1)
    expect.update({"??": 3, "--": 1})
    assert all(got[k] == v for k, v in expect.items())


def test_sequence_collection_annotation_db_assign_none():
    """assigning None to annotation_db breaks conection"""
    data = {"seq1": "ACGU", "seq2": "CGUA", "test_seq": "CCGU"}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="rna")
    seq_coll.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])
    seq_coll.add_feature(seqid="seq2", biotype="xyzzz", name="abc", spans=[(1, 2)])
    assert seq_coll.annotation_db is not None
    seq_coll.annotation_db = None
    assert seq_coll.annotation_db is None


def test_sequence_collection_annotation_db_assign_same(gff_db):
    """assigning the same annotation_db"""
    data = {"test_seq": "ACGT", "test_seq2": "CGTTTA"}
    seq_coll = new_alignment.make_unaligned_seqs(
        data, moltype="dna", annotation_db=gff_db
    )
    seq_coll.annotation_db = gff_db
    assert seq_coll.annotation_db is gff_db


@pytest.mark.parametrize("coll_maker", [new_alignment.make_unaligned_seqs])
def test_sequence_collection_get_annotations_from_any_seq(coll_maker):
    """get_annotations_from_any_seq returns correct annotations"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = coll_maker(data, moltype="dna")
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


def test_sequence_collection_add_feature(seqs):
    _ = seqs.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])
    assert seqs.annotation_db.num_matches() == 1

    # if seqid is not in seqs, raise error
    with pytest.raises(ValueError):
        _ = seqs.add_feature(seqid="bad_seq", biotype="xyz", name="abc", spans=[(1, 2)])


def test_sequence_collection_make_feature(seqs):
    data = {"seqid": "seq1", "biotype": "xyz", "name": "abc", "spans": [(1, 2)]}
    feat = seqs.make_feature(feature=data)
    assert isinstance(feat, Feature)
    assert seqs.annotation_db.num_matches() == 0


def test_sequence_collection_get_features(seqs):
    # no annotation db should return empty list
    got = seqs.get_features(seqid="seq1")
    assert not list(got)

    seqs.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 3)])
    seqs.add_feature(seqid="seq2", biotype="xyz", name="abc", spans=[(1, 3)])
    seqs.add_feature(seqid="seq3", biotype="xyz", name="abcde", spans=[(3, 5)])
    got = seqs.get_features(seqid="seq1")
    assert len(list(got)) == 1

    got = seqs.get_features(biotype="xyz")
    assert len(list(got)) == 3

    got = seqs.get_features(name="abc")
    assert len(list(got)) == 2

    got = seqs.get_features(biotype="xyz", name="abc")
    assert len(list(got)) == 2

    got = seqs.get_features(biotype="xyz", start=0, stop=5, allow_partial=True)
    assert len(list(got)) == 3

    got = seqs.get_features(start=2, stop=6, allow_partial=False)
    assert len(list(got)) == 1

    got = seqs.get_features(biotype="xyz", start=2, stop=4, allow_partial=True)
    assert len(list(got)) == 3

    with pytest.raises(ValueError):
        _ = list(
            seqs.get_features(
                seqid="bad_seq", biotype="xyz", start=2, stop=4, allow_partial=True
            )
        )


def test_sequence_collection_copy_annotations(gff_db):
    """copy_annotations copies records from annotation db"""
    data = {"seq1": "ACGU", "seq2": "CGUA", "test_seq": "CCGU"}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="rna")
    seq_coll.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])
    seq_coll.add_feature(seqid="seq2", biotype="xyzzz", name="abc", spans=[(1, 2)])
    expect = seq_coll.annotation_db.num_matches() + gff_db.num_matches()
    seq_coll.copy_annotations(gff_db)
    assert seq_coll.annotation_db.num_matches() == expect

    # copy annotations with no current annotations
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="rna")
    seq_coll.copy_annotations(gff_db)
    assert seq_coll.annotation_db.num_matches() == gff_db.num_matches()


def test_sequence_collection_copy_annotations_same_annotations(gff_db):
    data = {"seq1": "ACGU", "seq2": "CGUA", "test_seq": "CCGU"}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="rna")
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="rna")

    # copy annotations with the same annotation_db
    seq_coll.annotation_db = gff_db
    seq_coll.copy_annotations(gff_db)

    assert seq_coll.annotation_db.num_matches() == gff_db.num_matches()


def test_sequence_collection_copy_annotations_none_matching(gff_db):
    """copy annotations should old copy annotations for matching seqids"""
    data = {"name1": "ACGU", "name2": "CGUA", "name_3": "CCGU"}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="rna")
    seq_coll.add_feature(seqid="name1", biotype="xyz", name="abc", spans=[(1, 2)])
    seq_coll.add_feature(seqid="name2", biotype="xyzzz", name="abc", spans=[(1, 2)])
    expect = seq_coll.annotation_db.num_matches()
    assert gff_db.num_matches() > 0
    seq_coll.copy_annotations(gff_db)
    assert seq_coll.annotation_db.num_matches() == expect


def test_sequence_collection_copy_annotations_no_db(gff_db):
    data = {"seq1": "ACGU", "seq2": "CGUA", "test_seq": "CCGU"}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="rna")

    seq_coll.copy_annotations(gff_db)
    assert seq_coll.annotation_db.num_matches() == gff_db.num_matches()


def test_sequence_collection_get_seq_annotated():
    """SequenceCollection.get_seq should return specified seq"""

    data = {"seq1": "GATTTT", "seq2": "GATC??"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    seqs.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])
    seqs.add_feature(seqid="seq2", biotype="xy", name="ab", spans=[(1, 2)])

    with_annos = seqs.get_seq("seq1", copy_annotations=True)
    assert len(with_annos.annotation_db) == 1

    without_annos = seqs.get_seq("seq1", copy_annotations=False)
    assert len(without_annos.annotation_db) == 2


def _make_seq(name):
    raw_seq = "AACCCAAAATTTTTTGGGGGGGGGGCCCC"
    cds = (15, 25)
    utr = (12, 15)
    seq = new_moltype.DNA.make_seq(seq=raw_seq, name=name)
    seq.add_feature(biotype="CDS", name="CDS", spans=[cds])
    seq.add_feature(biotype="5'UTR", name="5' UTR", spans=[utr])
    return seq


def test_sequence_collection_init_seqs_have_annotations():
    """annotations on input seqs correctly merged and propagated"""
    data = {"seq1": _make_seq("seq1"), "seq2": _make_seq("seq2")}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    coll_db = seq_coll.annotation_db
    assert len(coll_db) == 4
    for name in seq_coll.names:
        seq = seq_coll.get_seq(name)
        db = seq.annotation_db
        assert db is coll_db


def test_sequence_collection_add_to_seq_updates_coll():
    """annotating a seq updates the db of the propagated"""
    data = {
        "x": "AACCCAAAATTTTTTGGGGGGGGGGCCCC",
        "y": "AACCCAAAATTTTTTGGGGGGGGGGCCCC",
    }
    seq_coll = new_alignment.make_unaligned_seqs(
        data,
        moltype="dna",
    )
    x = seq_coll.get_seq("x")
    assert len(seq_coll.annotation_db) == len(x.annotation_db) == 0
    x.add_feature(biotype="exon", name="E1", spans=[(3, 8)])
    assert len(seq_coll.annotation_db) == len(x.annotation_db) == 1


def test_sequence_collection_copy_annotations_incompat_type_fails(seqcoll_db, seqs):
    with pytest.raises(TypeError):
        seqcoll_db.copy_annotations(seqs)


@pytest.mark.parametrize("moltype", ("dna", "rna", "protein"))
def test_sequence_collection_to_moltype(moltype):
    data = dict(s1="ACGAA-", s2="ACCCAA")
    seqs = new_alignment.make_unaligned_seqs(data, moltype="text")
    mt_seqs = seqs.to_moltype(moltype=moltype)
    got = mt_seqs.seqs["s1"]
    assert got.moltype.label == moltype
    assert isinstance(got, new_sequence.Sequence)


@pytest.mark.parametrize("moltype", ("dna", "protein"))
def test_sequence_collection_to_moltype_same_moltype(moltype):
    data = dict(s1="ACGTT", s2="ACCTT")
    seqs = new_alignment.make_unaligned_seqs(data, moltype=moltype)
    got = seqs.to_moltype(moltype=moltype)
    assert got is seqs


def test_sequence_collection_to_moltype_with_gaps():
    """correctly convert to specified moltype"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="text")
    dna_seqs = seqs.to_moltype("dna")
    assert dna_seqs.moltype.label == "dna"
    assert dna_seqs.to_dict() == data

    # we can convert from text to dna to rna
    rna_seqs = dna_seqs.to_moltype("rna")
    assert rna_seqs.moltype.label == "rna"

    # but not from text to rna directly
    with pytest.raises(ValueError):
        seqs.to_moltype("rna")


def test_sequence_collection_to_moltype_info():
    """correctly convert to specified moltype"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = new_alignment.make_unaligned_seqs(
        data, moltype="text", info={"key": "value"}
    )
    dna = seqs.to_moltype("dna")
    assert dna.info["key"] == "value"


def test_sequence_collection_to_moltype_annotation_db():
    """correctly convert to specified moltype"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="text")
    db = GffAnnotationDb()
    db.add_feature(seqid="seq1", biotype="exon", name="annotation1", spans=[(3, 8)])
    db.add_feature(seqid="seq2", biotype="exon", name="annotation2", spans=[(1, 2)])
    db.add_feature(seqid="seq3", biotype="exon", name="annotation3", spans=[(3, 6)])
    seqs.annotation_db = db
    dna = seqs.to_moltype("dna")
    assert len(dna.annotation_db) == 3


def test_sequence_collection_to_dna():
    """correctly convert to dna"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="text")

    # convert from text to dna
    dna_seqs = seqs.to_dna()
    assert dna_seqs.moltype.label == "dna"
    assert dna_seqs.to_dict() == data

    # convert from rna to dna
    data = {"seq1": "ACGUACGUA", "seq2": "ACCGAA---", "seq3": "ACGUACGUU"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="RNA")
    dna_seqs = seqs.to_dna()
    assert dna_seqs.moltype.label == "dna"
    assert dna_seqs.to_dict() == {
        name: seq.replace("U", "T") for name, seq in data.items()
    }


def test_sequence_collection_to_rna():
    data = {"seq1": "ACGUACGUA", "seq2": "ACCGAA---", "seq3": "ACGUACGUU"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="text")

    # convert from text to rna
    rna_seqs = seqs.to_rna()
    assert rna_seqs.moltype.label == "rna"
    assert rna_seqs.to_dict() == data

    # convert from dna to rna
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    rna_seqs = seqs.to_rna()
    assert rna_seqs.moltype.label == "rna"
    assert rna_seqs.to_dict() == {
        name: seq.replace("T", "U") for name, seq in data.items()
    }


@pytest.fixture(scope="function")
def dotplot_seqs():
    data = {
        "Human": "CAGATTTGGCAGTT-",
        "Mouse": "CAGATTCAGCAGGTG",
        "Rat": "CAGATTCAGCAGG--G",
    }
    return new_alignment.make_unaligned_seqs(data, moltype="dna")


def test_sequence_collection_dotplot(dotplot_seqs):
    """exercising dotplot method"""
    # with provided names
    plot = dotplot_seqs.dotplot(name1="Human", name2="Mouse")
    assert str(plot.seq1) != str(plot.seq2)

    # without providing names
    plot = dotplot_seqs.dotplot()
    assert str(plot.seq1) != str(plot.seq2)

    # providing only one name
    plot = dotplot_seqs.dotplot(name1="Human")
    assert str(plot.seq1) != str(plot.seq2)

    # a collection of one sequence should make dotplot with itself
    less_seqs = dotplot_seqs.take_seqs("Human")
    plot = less_seqs.dotplot()
    assert str(plot.seq1) == str(plot.seq2)

    # k larger than window should raise an error
    with pytest.raises(AssertionError):
        dotplot_seqs.dotplot(window=5, k=11)

    # names not in the collection should raise an error
    with pytest.raises(ValueError):
        dotplot_seqs.dotplot(name1="Human", name2="Dog")


def test_sequence_collection_dotplot_annotated():
    """exercising dotplot method with annotated sequences"""
    db = GffAnnotationDb()
    db.add_feature(seqid="Human", biotype="exon", name="fred", spans=[(10, 15)])

    data = {"Human": "CAGATTTGGCAGTT-", "Mouse": "CAGATTCAGCAGGTG"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    seqs.annotation_db = db
    seqs = seqs.take_seqs(["Human", "Mouse"], copy_annotations=True)
    _ = seqs.dotplot(show_progress=False)


@pytest.mark.parametrize(
    "collection_maker",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_add_seqs(collection_maker):
    data = dict(
        [("name1", "AAA"), ("name2", "A--"), ("name3", "AAA"), ("name4", "AAA")]
    )
    data2 = dict([("name5", "TTT"), ("name6", "---")])
    aln = collection_maker(data, moltype="dna", info={"key": "foo"})
    out_aln = aln.add_seqs(data2)
    assert len(out_aln.names) == 6


@pytest.mark.parametrize(
    "collection_maker",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_add_seqs_duplicate_raises(collection_maker):
    """add_seqs should raise an error if duplicate names"""
    data = dict(
        [("name1", "AAA"), ("name2", "AAA"), ("name3", "AAA"), ("name4", "AAA")]
    )
    data2 = dict([("name5", "BBB"), ("name1", "CCC")])
    aln = collection_maker(data, moltype="text", info={"key": "foo"})
    with pytest.raises(ValueError):
        _ = aln.add_seqs(data2)

    # but it should work if we allow duplicates
    out_aln = aln.add_seqs(data2, force_unique_keys=False)
    assert len(out_aln.names) == 5


@pytest.mark.parametrize(
    "collection_maker",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_add_seqs_info(collection_maker):
    """add_seqs should preserve info attribute"""
    data = dict(
        [("name1", "AAA"), ("name2", "AAA"), ("name3", "AAA"), ("name4", "AAA")]
    )
    data2 = dict([("name5", "BBB"), ("name6", "CCC")])
    aln = collection_maker(data, moltype="text", info={"key": "foo"})
    aln2 = collection_maker(data2, moltype="text", info={"key": "bar"})
    out_aln = aln.add_seqs(aln2)
    assert out_aln.info["key"] == "foo"


@pytest.mark.xfail(reason="design of strand attribute is not finalized")
@pytest.mark.parametrize(
    "collection_maker",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_add_seqs_strand(collection_maker):
    """add_seqs should preserve strand attribute, and accept new strand information
    for the new sequences"""
    data = dict(
        [("name1", "AAA"), ("name2", "AAA"), ("name3", "AAA"), ("name4", "AAA")]
    )
    strand = {"name2": -1}
    data2 = dict([("name5", "TTT"), ("name6", "CCC")])
    strand2 = {
        "name5": -1,
    }
    aln = collection_maker(data, moltype="dna", info={"key": "foo"}, strand=strand)
    out_aln = aln.add_seqs(data2, strand=strand2)
    assert out_aln.seqs.strand["name2"] == -1
    assert out_aln.seqs.strand["name5"] == -1


@pytest.mark.parametrize(
    "collection_maker",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_write(collection_maker, tmp_path):
    """SequenceCollection.write should write in correct format"""
    data = {"a": "AAAA", "b": "TTTT", "c": "CCCC"}
    seqs = collection_maker(data, moltype="dna")
    fn = tmp_path / "seqs.fasta"
    seqs.write(fn)
    with open(fn, newline=None) as infile:
        result = infile.read()
    assert result == ">a\nAAAA\n>b\nTTTT\n>c\nCCCC\n"


@pytest.mark.parametrize(
    "collection_maker",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_write_gapped(collection_maker, tmp_path):
    data = {"a": "AAA--", "b": "TTTTT", "c": "CCCCC"}
    seqs = collection_maker(data, moltype="dna")
    fn = tmp_path / "seqs.fasta"
    seqs.write(fn)
    with open(fn, newline=None) as infile:
        result = infile.read()
    assert result == ">a\nAAA--\n>b\nTTTTT\n>c\nCCCCC\n"


def test_get_ambiguous_positions():
    """SequenceCollection.get_ambiguous_positions should return pos"""
    aln = new_alignment.make_unaligned_seqs(
        {"s1": "ATGRY?", "s2": "T-AG??"}, moltype="dna"
    )
    assert aln.get_ambiguous_positions() == {
        "s2": {4: "?", 5: "?"},
        "s1": {3: "R", 4: "Y", 5: "?"},
    }


def test_sequence_collection_consistent_gap_degen_handling():
    """gap degen character should be treated consistently"""
    # the degen character '?' can be a gap, so when we strip gaps it should
    # be gone too
    raw_seq = "---??-??TC-GGCG-GCA-G-GC-?-C-TAN-GCGC-CCTC-AGGA?-???-??--"
    raw_ungapped = re.sub("[-?]", "", raw_seq)
    re.sub("[N?]+", "", raw_seq)
    dna = new_moltype.DNA.make_seq(seq=raw_seq)

    aln = new_alignment.make_unaligned_seqs({"a": dna, "b": dna}, moltype="dna")
    expect = new_alignment.make_unaligned_seqs(
        {"a": raw_ungapped, "b": raw_ungapped}, moltype="dna"
    ).to_fasta()
    assert aln.degap().to_fasta() == expect


def test_sequence_collection_pad_seqs(ragged):
    # pad to max length
    padded1 = ragged.pad_seqs()
    seqs1 = list(padded1.iter_seqs(seq_order=["a", "b", "c"]))
    assert list(map(str, seqs1)) == ["AAAAAA", "AAA---", "AAAA--"]

    # pad to alternate length
    padded2 = ragged.pad_seqs(pad_length=10)
    seqs2 = list(padded2.iter_seqs(seq_order=["a", "b", "c"]))
    assert list(map(str, seqs2)) == ["AAAAAA----", "AAA-------", "AAAA------"]

    # assertRaises error when pad_length is less than max seq length
    with pytest.raises(ValueError):
        _ = ragged.pad_seqs(pad_length=5)


@pytest.fixture
def ambigs_coll():
    data = {"a": "AAAA??????", "b": "CCCGGG--NN"}
    return new_alignment.make_unaligned_seqs(data, moltype="dna")


def test_sequence_collection_get_lengths(ambigs_coll):
    got = ambigs_coll.get_lengths()
    expect = {"a": 4, "b": 6}
    assert got == expect


def test_sequence_collection_get_lengths_include_ambiguity(ambigs_coll):
    # NOTE '?' is excluded as it could be a gap
    got = ambigs_coll.get_lengths(include_ambiguity=True)
    expect = {"a": 4, "b": 8}
    assert got == expect


def test_sequence_collection_get_lengths_allow_gap(ambigs_coll):
    got = ambigs_coll.get_lengths(include_ambiguity=True, allow_gap=True)
    expect = {"a": 10, "b": 10}
    assert got == expect


def test_sequence_collection_strand_symmetry():
    """exercising strand symmetry test"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    result = seqs.strand_symmetry()
    assert numpy.allclose(result["seq1"].observed.array, [[3, 2], [2, 2]])
    assert numpy.allclose(result["seq2"].observed.array, [[3, 0], [2, 1]])


def test_sequence_collection_rename_seqs():
    """successfully rename sequences"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    new = seqs.rename_seqs(lambda x: x.upper())
    expect = {n.upper() for n in data}
    assert set(new.names) == expect


def test_sequence_collection_apply_pssm():
    """should successfully produce pssm scores"""
    from cogent3.parse import jaspar

    _, pwm = jaspar.read("data/sample.jaspar")
    data = {
        "ENSMUSG00000056468": "GCCAGGGGGGAAAGGGAGAA",
        "ENSMUSG00000039616": "GCCCTTCAAATTTGGTTTCT",
        "ENSMUSG00000024091": "TTTCCAGGGCAGACAAGACG",
        "ENSMUSG00000024056": "ACAATAATGCCGAGAGCCAG",
        "ENSMUSG00000054321": "TATGAAAATTTTTGCCAGGC",
        "ENSMUSG00000052469": "CCTGTTTGCCTTTAAATATT",
        "ENSMUSG00000024261": "CAGACAAGAAACCAGCAACA",
        "ENSMUSG00000052031": "AGCGAGTATNCACGCACAGA",
        "ENSMUSG00000067872": "ACACAGCTCTGACAACTCAT",
        "ENSMUSG00000023892": "GTAACATCAGTACAGCACAG",
    }
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    max_seq_len = max(seqs.get_lengths())

    scores = seqs.apply_pssm(path="data/sample.jaspar", show_progress=False)
    assert scores.shape == (len(data), max_seq_len - pwm.shape[0] + 1)
    scores = seqs.apply_pssm(pssm=pwm.to_pssm(), show_progress=False)
    assert scores.shape == (len(data), max_seq_len - pwm.shape[0] + 1)

    # using the names argument works to return scores in the correct order
    seqs = new_alignment.make_unaligned_seqs(
        {"ENSMUSG00000056468": "GCCAGGGGGGAAAGGGAGAA"}, moltype="dna"
    )
    expect = []
    expect.extend(seqs.apply_pssm(pssm=pwm.to_pssm(), show_progress=False))
    seqs = new_alignment.make_unaligned_seqs(
        {"ENSMUSG00000039616": "GCCCTTCAAATTTGGTTTCT"}, moltype="dna"
    )
    expect.extend(seqs.apply_pssm(pssm=pwm.to_pssm(), show_progress=False))
    expect = numpy.array(expect)
    seqs = new_alignment.make_unaligned_seqs(
        {
            "ENSMUSG00000056468": "GCCAGGGGGGAAAGGGAGAA",
            "ENSMUSG00000039616": "GCCCTTCAAATTTGGTTTCT",
        },
        moltype="dna",
    )
    got = seqs.apply_pssm(
        pssm=pwm.to_pssm(), show_progress=False, names="ENSMUSG00000056468"
    )
    assert numpy.allclose(got, expect[:1])
    got = seqs.apply_pssm(
        pssm=pwm.to_pssm(), show_progress=False, names=["ENSMUSG00000039616"]
    )
    assert numpy.allclose(got, expect[1:])
    got = seqs.apply_pssm(
        pssm=pwm.to_pssm(),
        show_progress=False,
        names=["ENSMUSG00000056468", "ENSMUSG00000039616"],
    )
    assert numpy.allclose(got, expect)
    got = seqs.apply_pssm(
        pssm=pwm.to_pssm(),
        show_progress=False,
        names=["ENSMUSG00000039616", "ENSMUSG00000056468"],
    )
    assert numpy.allclose(got, expect[::-1])


def test_sequence_collection_apply_pssm2():
    """apply_pssm fail if ragged sequences"""
    data = {
        "ENSMUSG00000056468": "GCCAGGGGGGAAAGGGAGAA",
        "ENSMUSG00000039616": "GCCCTTCAAATTT",
    }
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    with pytest.raises(AssertionError):
        _ = seqs.apply_pssm(path="data/sample.jaspar", show_progress=False)


def test_sequence_collection_get_seq_entropy():
    """get_seq_entropy should get entropy of each seq"""
    a = new_alignment.make_unaligned_seqs(dict(a="ACCC", b="AGTA"), moltype="dna")
    entropy = a.entropy_per_seq()
    e = 0.81127812445913283  # sum(p log_2 p) for p = 0.25, 0.75
    assert numpy.allclose(entropy, numpy.array([e, 1.5]))


def test_sequence_collection_write_to_json(tmp_path):
    # test writing to json file
    data = {"a": "AAAA", "b": "TTTT", "c": "CCCC"}
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    path = str(tmp_path / "sample.json")
    seq_coll.write(path)
    with open_(path) as fn:
        got = json.loads(fn.read())
        assert got == seq_coll.to_rich_dict()


def test_sequence_collection_to_rich_dict():
    """to_rich_dict produces correct dict"""
    data = {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")

    got = seqs.to_rich_dict()
    seqs_data = {
        "init_args": {
            "data": {
                name: seqs._seqs_data.get_seq_str(seqid=name) for name in seqs.names
            },
            "alphabet": seqs.moltype.most_degen_alphabet().to_rich_dict(),
            "strand": seqs._seqs_data._strand,
            "offset": seqs._seqs_data._offset,
            "reversed": seqs._seqs_data.is_reversed,
        },
        "type": get_object_provenance(seqs._seqs_data),
        "version": __version__,
    }
    expect = {
        "seqs_data": seqs_data,
        "type": get_object_provenance(seqs),
        "version": __version__,
        "init_args": {
            "moltype": seqs.moltype.label,
            "names": seqs.names,
            "info": seqs.info,
        },
    }
    assert got == expect


def test_sequence_collection_to_rich_dict_annotation_db():
    data = {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    seqs.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])

    got = seqs.to_rich_dict()
    db = seqs.annotation_db.to_rich_dict()
    seqs_data = {
        "init_args": {
            "data": {
                name: seqs._seqs_data.get_seq_str(seqid=name) for name in seqs.names
            },
            "alphabet": seqs.moltype.most_degen_alphabet().to_rich_dict(),
            "strand": seqs._seqs_data._strand,
            "offset": seqs._seqs_data._offset,
            "reversed": seqs._seqs_data.is_reversed,
        },
        "type": get_object_provenance(seqs._seqs_data),
        "version": __version__,
    }
    expect = {
        "seqs_data": seqs_data,
        "init_args": {
            "moltype": seqs.moltype.label,
            "names": seqs.names,
            "info": seqs.info,
            "annotation_db": db,
        },
        "type": get_object_provenance(seqs),
        "version": __version__,
    }
    assert got == expect


def test_sequence_collection_to_rich_dict_reversed_seqs():
    data = {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna", strand={"seq1": -1})
    reversed_seqs = seqs.reverse_complement()

    got = reversed_seqs.to_rich_dict()
    seqs_data = {
        "init_args": {
            "data": {
                name: seqs._seqs_data.get_seq_str(seqid=name) for name in seqs.names
            },
            "alphabet": seqs.moltype.most_degen_alphabet().to_rich_dict(),
            "strand": {"seq1": -1},
            "offset": seqs._seqs_data._offset,
            "reversed": True,
        },
        "type": get_object_provenance(seqs._seqs_data),
        "version": __version__,
    }
    expect = {
        "seqs_data": seqs_data,
        "init_args": {
            "moltype": seqs.moltype.label,
            "names": seqs.names,
            "info": seqs.info,
        },
        "type": get_object_provenance(seqs),
        "version": __version__,
    }
    assert got == expect


def test_sequence_collection_to_json():
    """roundtrip of to_json produces correct dict"""
    seq_coll = new_alignment.make_unaligned_seqs(
        {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}, moltype="dna"
    )
    txt = seq_coll.to_json()
    got = json.loads(txt)
    expect = seq_coll.to_rich_dict()
    assert got == expect


@pytest.mark.parametrize("reverse", (False, True))
def test_sequence_collection_round_trip(reverse):
    seq_coll = new_alignment.make_unaligned_seqs(
        {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}, moltype="dna"
    )
    seq_coll = seq_coll.reverse_complement() if reverse else seq_coll

    rd = seq_coll.to_rich_dict()
    got = deserialise_object(rd)
    assert isinstance(got, new_alignment.SequenceCollection)
    assert got.to_rich_dict() == seq_coll.to_rich_dict()


@pytest.mark.parametrize("moltype", ("dna", "rna"))
def test_sequence_collection_distance_matrix_singleton_collection(moltype):
    """SequenceCollection.distance_matrix() should raise error if collection
    only contains a single sequence"""
    collection = new_alignment.make_unaligned_seqs(
        {"s1": "ACGTACGTAGTCGCG"}, moltype=moltype
    )
    with pytest.raises(ValueError):
        _ = collection.distance_matrix()


@pytest.mark.parametrize("moltype", ("dna", "rna"))
def test_sequence_collection_distance_matrix_same_seq(moltype):
    """Identical seqs should return distance measure of 0.0"""
    data = dict([("s1", "ACGACGAGCGCG"), ("s2", "GGACGACGCG"), ("s3", "GGACGACGCG")])
    collection = new_alignment.make_unaligned_seqs(data, moltype=moltype)
    with catch_warnings():
        filterwarnings("ignore", category=DeprecationWarning)
        dists = collection.distance_matrix(calc="pdist")

    # all comparison of a sequence to itself should be zero
    for seq in collection.names:
        assert dists[(seq, seq)] == 0.0

    # s2 and s3 are identical, so should be zero
    assert dists[("s2", "s3")] == 0.0
    assert dists[("s3", "s2")] == 0.0

    # s1 is different from s2 and s3
    assert dists[("s1", "s2")] > 0.0
    assert dists[("s1", "s3")] > 0.0


@pytest.mark.parametrize("moltype", ("protein", "text", "bytes"))
def test_sequence_collection_distance_matrix_raises_wrong_moltype(moltype):
    data = {"s1": "ACGTA", "s2": "ACGTA"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype=moltype)
    with pytest.raises(NotImplementedError):
        seqs.distance_matrix()


@pytest.mark.parametrize("moltype", ("dna", "rna"))
def test_sequence_collection_distance_matrix_passes_correct_moltype(moltype):
    data = {"s1": "ACGTA", "s2": "ACGTA"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype=moltype)

    seqs.distance_matrix()


@pytest.mark.parametrize(
    "collection_maker",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_sequence_collection_reverse_complement(collection_maker):
    data = {"s1": "AACC", "s2": "GGTT"}
    seqs = collection_maker(data, moltype="dna")
    got = seqs.reverse_complement()
    expect = {"s1": "GGTT", "s2": "AACC"}
    assert got.to_dict() == expect
    assert str(got.get_seq(seqname="s1")) == "GGTT"

    got = got.reverse_complement()
    assert got.to_dict() == data


# Alignment tests


@pytest.fixture
def aligned_dict():
    return dict(seq1="ACG--T", seq2="-CGAAT", seq3="------", seq4="--GA--")


@pytest.fixture
def gapped_seqs_dict():
    return dict(
        seq1="AA--CC",
        seq2="--AA--",
        seq3="-A-A-A",
        seq4="AAA---",
        seq5="---AAA",
        seq6="------",
        seq7="AAAAAA",
    )


@pytest.fixture
def aligned_array_dict():
    return dict(
        seq1=numpy.array([1, 2, 3, 4, 4], dtype=numpy.uint8),
        seq2=numpy.array([4, 2, 4, 3, 4], dtype=numpy.uint8),
        seq3=numpy.array([4, 4, 4, 4, 4], dtype=numpy.uint8),
        seq4=numpy.array([4, 4, 1, 4, 4], dtype=numpy.uint8),
    )


@pytest.fixture
def aligned_seqs_data(aligned_dict, dna_alphabet, dna_make_seq):
    return new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )


@pytest.fixture
def gap_seqs():
    return [
        ("A---CTG-C", [[1, 3], [4, 4]]),
        ("-GTAC----", [[0, 1], [4, 5]]),
        ("---A--T--", [[0, 3], [1, 5], [2, 7]]),
    ]


@pytest.mark.parametrize("i", range(3))
def test_seq_to_gap_coords_str(gap_seqs, i, dna_alphabet):
    seq, gap_coords = gap_seqs[i]
    got_ungapped, got_map = new_alignment.seq_to_gap_coords(seq, alphabet=dna_alphabet)
    expect = dna_alphabet.to_indices(seq.replace("-", ""))
    assert numpy.array_equal(got_ungapped, expect)
    assert numpy.array_equal(got_map, gap_coords)


def test_seq_to_gap_coords_str_all_gaps(dna_alphabet):
    parent_seq = "-----"
    expect_gaplen = numpy.array([len(parent_seq)])
    got_ungap, got_map = new_alignment.seq_to_gap_coords(
        parent_seq, alphabet=dna_alphabet
    )
    expect = numpy.array([])
    assert numpy.array_equal(got_ungap, expect)
    assert got_map[:, 1] == expect_gaplen


def test_seq_to_gap_coords_str_no_gaps(dna_alphabet):
    parent_seq = "ACTGC"
    got_ungap, got_map = new_alignment.seq_to_gap_coords(
        parent_seq, alphabet=dna_alphabet
    )
    expect = dna_alphabet.to_indices(parent_seq)
    assert numpy.array_equal(got_ungap, expect)
    assert got_map.size == 0


def test_seq_to_gap_coords_arr_all_gaps(dna_alphabet):
    parent_seq = dna_alphabet.to_indices("-----")
    got_ungap, got_map = new_alignment.seq_to_gap_coords(
        parent_seq, alphabet=dna_alphabet
    )
    assert got_ungap.size == 0
    assert numpy.array_equal(got_map, numpy.array([[0, 5]]))


def test_seq_to_gap_coords_arr_no_gaps(dna_alphabet):
    parent_seq = dna_alphabet.to_indices("ACTGC")
    got_ungap, got_empty_arr = new_alignment.seq_to_gap_coords(
        parent_seq, alphabet=dna_alphabet
    )
    assert numpy.array_equal(got_ungap, parent_seq)
    assert got_empty_arr.size == 0


@pytest.mark.parametrize("i", range(3))
def test_seq_to_gap_coords_arr(gap_seqs, i, dna_alphabet):
    seq, gap_coords = gap_seqs[i]
    seq = dna_alphabet.to_indices(seq)
    got_ungapped, got_map = new_alignment.seq_to_gap_coords(seq, alphabet=dna_alphabet)
    expect = seq[seq != 4]  # gap_char = 4
    assert numpy.array_equal(got_ungapped, expect)
    assert numpy.array_equal(got_map, gap_coords)


@pytest.mark.parametrize("i", range(3))
def test_seq_to_gap_coords_arr_dispatch_equal(gap_seqs, i, dna_alphabet):
    """seq_to_gap_coords should return the same gap coords when input is a string/array/bytes"""
    seq_str, _ = gap_seqs[i]
    seq_array = dna_alphabet.to_indices(seq_str)
    seq_bytes = dna_alphabet.array_to_bytes(seq_array)
    seq_from_str, gaps_from_str = new_alignment.seq_to_gap_coords(
        seq_str, alphabet=dna_alphabet
    )
    seq_from_array, gaps_from_arr = new_alignment.seq_to_gap_coords(
        seq_array, alphabet=dna_alphabet
    )
    seq_from_bytes, gaps_from_bytes = new_alignment.seq_to_gap_coords(
        seq_bytes, alphabet=dna_alphabet
    )
    assert numpy.array_equal(gaps_from_str, gaps_from_arr)
    assert numpy.array_equal(gaps_from_arr, gaps_from_bytes)
    assert numpy.array_equal(seq_from_str, seq_from_array)
    assert numpy.array_equal(seq_from_array, seq_from_bytes)


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_str(aligned_dict, seqid, dna_moltype):
    """str() of an Aligned instance should return gapped sequence cast to a string"""
    aln = new_alignment.make_aligned_seqs(aligned_dict, moltype=dna_moltype)
    aligned = aln.seqs[seqid]
    got = str(aligned)
    expect = aligned_dict[seqid]
    assert got == expect


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_array(aligned_dict, seqid, dna_moltype):
    """array() of an Aligned instance should return gapped sequence as an array"""
    aln = new_alignment.make_aligned_seqs(aligned_dict, moltype=dna_moltype)
    aligned = aln.seqs[seqid]
    got = numpy.array(aligned)
    expect = dna_moltype.degen_gapped_alphabet.to_indices(aligned_dict[seqid])
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_bytes(aligned_dict, seqid, dna_moltype):
    aln = new_alignment.make_aligned_seqs(aligned_dict, moltype=dna_moltype)
    aligned = aln.seqs[seqid]
    got = bytes(aligned)
    expect = dna_moltype.degen_gapped_alphabet.array_to_bytes(
        dna_moltype.degen_gapped_alphabet.to_indices(aligned_dict[seqid])
    )
    assert got == expect


@pytest.mark.parametrize("start", range(6))
@pytest.mark.parametrize("stop", range(6))
@pytest.mark.parametrize("step", range(1, 4))
@pytest.mark.parametrize(
    "seqid", ("seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7")
)
def test_aligned_getitem_slice(gapped_seqs_dict, seqid, start, stop, step):
    aln = new_alignment.make_aligned_seqs(gapped_seqs_dict, moltype="dna")
    al = aln.seqs[seqid]
    got = al[start:stop:step]
    expect = gapped_seqs_dict[seqid][start:stop:step]
    assert str(got) == expect

    got1 = got.gapped_seq
    assert got1 == expect

    # Aligned.seq should return the sequence without gaps
    got2 = got.seq
    expect = expect.replace("-", "")
    assert got2 == expect


@pytest.mark.parametrize(
    "seqid", ("seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7")
)
@pytest.mark.parametrize("i", range(6))
def test_aligned_getitem_int(gapped_seqs_dict, seqid, i):
    aln = new_alignment.make_aligned_seqs(gapped_seqs_dict, moltype="dna")
    al = aln.seqs[seqid]
    assert isinstance(al, new_alignment.Aligned)
    got = al[i]
    expect = gapped_seqs_dict[seqid][i]
    assert str(got) == expect


def test_aligned_getitem_raises(gapped_seqs_dict):
    aln = new_alignment.make_aligned_seqs(gapped_seqs_dict, moltype="dna")
    al = aln.seqs["seq1"]
    with pytest.raises(NotImplementedError):
        _ = al[(9.0, 20.0)]


@pytest.mark.parametrize(
    "seqid", ("seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7")
)
def test_aligned_iter(seqid, gapped_seqs_dict):
    aln = new_alignment.make_aligned_seqs(gapped_seqs_dict, moltype="dna")
    aligned = aln.seqs[seqid]
    for i, got in enumerate(aligned):
        expect = gapped_seqs_dict[seqid][i]  # directly index the sequence
        assert got == expect


@pytest.mark.parametrize("seqid", ("seq0", "seq1", "seq2"))
def test_aligned_seqs_data_init(seqid, gap_seqs, dna_alphabet, dna_make_seq):
    """test initialisation of AlignedSeqsData object with dictionary of ungapped
    sequences and dictionary of gap coordinates and cumulated gap lengths"""
    seqs = {f"seq{i}": seq.replace("-", "") for i, (seq, _) in enumerate(gap_seqs)}
    gaps = {f"seq{i}": numpy.array(gaps) for i, (_, gaps) in enumerate(gap_seqs)}
    ad = new_alignment.AlignedSeqsData(
        seqs=seqs, gaps=gaps, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    assert ad.get_seq_str(seqid=seqid) == seqs[seqid]
    assert numpy.array_equal(ad.get_gaps(seqid=seqid), gaps[seqid])


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4", "seq5"))
@pytest.mark.parametrize("data_type", (str, numpy.array, bytes))
def test_aligned_seqs_data_init_gapped(
    gapped_seqs_dict, seqid, data_type, dna_alphabet, dna_moltype, dna_make_seq
):
    """AlignedSeqsData should handle data from gapped sequences correctly,
    including edge cases (all gaps, no gaps, etc)"""

    typed_data = {
        name: make_typed(seq, data_type=data_type, moltype=dna_moltype)
        for name, seq in gapped_seqs_dict.items()
    }

    seq_data = {
        name: new_alignment.seq_to_gap_coords(seq, alphabet=dna_alphabet)[0]
        for name, seq in typed_data.items()
    }
    gap_data = {
        name: new_alignment.seq_to_gap_coords(seq, alphabet=dna_alphabet)[1]
        for name, seq in typed_data.items()
    }
    asd = new_alignment.AlignedSeqsData(
        seqs=seq_data, gaps=gap_data, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    assert asd.align_len == 6
    assert asd.get_gapped_seq_str(seqid=seqid) == gapped_seqs_dict[seqid]


@pytest.mark.parametrize("data_type", (str, numpy.array, bytes))
def test_aligned_seqs_data_unequal_seqlens_raises(
    data_type, dna_alphabet, dna_moltype, dna_make_seq
):
    """from_aligned_seqs should raise an error if sequences are of unequal length"""
    data = dict(
        seq1=make_typed("A-A", data_type=data_type, moltype=dna_moltype),
        seq2=make_typed("AAAAAAA--", data_type=data_type, moltype=dna_moltype),
    )
    with pytest.raises(ValueError):
        _ = new_alignment.AlignedSeqsData.from_aligned_seqs(
            data=data, alphabet=dna_alphabet, make_seq=dna_make_seq
        )
    # directly creating an AlignedSeqsData object should also raise an error
    seq_data = {
        name: new_alignment.seq_to_gap_coords(seq, alphabet=dna_alphabet)[0]
        for name, seq in data.items()
    }
    gap_data = {
        name: new_alignment.seq_to_gap_coords(seq, alphabet=dna_alphabet)[1]
        for name, seq in data.items()
    }
    with pytest.raises(ValueError):
        _ = new_alignment.AlignedSeqsData(
            seqs=seq_data, gaps=gap_data, alphabet=dna_alphabet, make_seq=dna_make_seq
        )


def test_aligned_seqs_data_diff_keys_raises(dna_alphabet, dna_make_seq):
    """AlignedSeqsData expect identical keys in seqs and gaps"""
    seqs = dict(seq1=numpy.array([2, 1]), seq2=numpy.array([2, 0, 3, 1]))
    gaps = dict(seq1=numpy.array([[1, 3]]), seq3=numpy.array([[0, 1]]))

    with pytest.raises(ValueError):
        _ = new_alignment.AlignedSeqsData(
            seqs=seqs, gaps=gaps, alphabet=dna_alphabet, make_seq=dna_make_seq
        )
    # assert that it would work if we indeed had the same keys
    gaps["seq2"] = gaps.pop("seq3")
    asd = new_alignment.AlignedSeqsData(
        seqs=seqs, gaps=gaps, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    assert asd.get_seq_str(seqid="seq1") == "AC"


def test_aligned_seqs_data_omit_seqs_gaps_raises(dna_alphabet, dna_make_seq):
    with pytest.raises(ValueError):
        _ = new_alignment.AlignedSeqsData(
            seqs={}, gaps={}, alphabet=dna_alphabet, make_seq=dna_make_seq
        )


def test_aligned_seqs_data_names(aligned_dict, dna_alphabet, dna_make_seq):
    got = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    assert isinstance(got, new_alignment.AlignedSeqsData)
    assert got.names == ["seq1", "seq2", "seq3", "seq4"]


def test_aligned_seqs_data_len(aligned_dict, dna_alphabet, dna_make_seq):
    got = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    assert len(got) == len(aligned_dict["seq1"])
    assert got.align_len == len(aligned_dict["seq1"])


@pytest.mark.parametrize(
    "seqid, index", [("seq1", 0), ("seq2", 1), ("seq3", 2), ("seq4", 3)]
)
@pytest.mark.parametrize("moltype", ("dna_moltype", "rna_moltype"))
def test_aligned_seqs_data_getitem(seqid, index, aligned_array_dict, moltype, request):
    moltype = request.getfixturevalue(moltype)
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_array_dict,
        alphabet=moltype.degen_gapped_alphabet,
        make_seq=moltype.make_seq,
    )
    got_with_seqid = numpy.array(ad[seqid])
    got_with_index = numpy.array(ad[index])

    expect = aligned_array_dict[seqid][aligned_array_dict[seqid] != 4]

    assert numpy.array_equal(got_with_seqid, expect)
    assert numpy.array_equal(got_with_index, expect)


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
@pytest.mark.parametrize("moltype", ("dna_moltype", "rna_moltype"))
def test_aligned_seqs_data_get_seq_array(aligned_array_dict, seqid, moltype, request):
    moltype = request.getfixturevalue(moltype)
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_array_dict,
        alphabet=moltype.degen_gapped_alphabet,
        make_seq=moltype.make_seq,
    )
    got = ad.get_seq_array(seqid=seqid)
    expect = aligned_array_dict[seqid]  # get original data
    expect = expect[expect != 4]  # remove gaps
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
@pytest.mark.parametrize("moltype", ("dna_moltype", "rna_moltype"))
def test_aligned_seqs_data_get_gapped_seq_array(
    aligned_array_dict, seqid, moltype, request
):
    moltype = request.getfixturevalue(moltype)
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_array_dict,
        alphabet=moltype.degen_gapped_alphabet,
        make_seq=moltype.make_seq,
    )
    got = ad.get_gapped_seq_array(seqid=seqid)  # get original data
    expect = aligned_array_dict[seqid]
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
@pytest.mark.parametrize("moltype", ("dna_moltype", "rna_moltype"))
def test_aligned_seqs_data_get_seq_str(aligned_array_dict, seqid, moltype, request):
    moltype = request.getfixturevalue(moltype)
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_array_dict,
        alphabet=moltype.degen_gapped_alphabet,
        make_seq=moltype.make_seq,
    )
    got = ad.get_seq_str(seqid=seqid)
    expect = moltype.degen_gapped_alphabet.from_indices(aligned_array_dict[seqid])
    expect = expect.replace("-", "")  # remove gaps
    assert got == expect


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
@pytest.mark.parametrize("moltype", ("dna_moltype", "rna_moltype"))
def test_aligned_seqs_data_get_gapped_seq_str(
    aligned_array_dict, seqid, moltype, request
):
    moltype = request.getfixturevalue(moltype)
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_array_dict,
        alphabet=moltype.degen_gapped_alphabet,
        make_seq=moltype.make_seq,
    )
    got = ad.get_gapped_seq_str(seqid=seqid)
    expect = moltype.degen_gapped_alphabet.from_indices(aligned_array_dict[seqid])
    assert got == expect


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
@pytest.mark.parametrize("moltype", ("dna_moltype", "rna_moltype"))
def test_aligned_seqs_data_get_seq_bytes(aligned_array_dict, seqid, moltype, request):
    moltype = request.getfixturevalue(moltype)
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_array_dict,
        alphabet=moltype.degen_gapped_alphabet,
        make_seq=moltype.make_seq,
    )
    got = ad.get_seq_bytes(seqid=seqid)
    raw = aligned_array_dict[seqid]
    raw = raw[raw != 4]
    expect = moltype.degen_gapped_alphabet.array_to_bytes(raw)
    assert got == expect


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
@pytest.mark.parametrize("moltype", ("dna_moltype", "rna_moltype"))
def test_aligned_seqs_data_get_gapped_seq_bytes(
    aligned_array_dict, seqid, moltype, request
):
    moltype = request.getfixturevalue(moltype)
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_array_dict,
        alphabet=moltype.degen_gapped_alphabet,
        make_seq=moltype.make_seq,
    )
    got = ad.get_gapped_seq_bytes(seqid=seqid)
    raw = aligned_array_dict[seqid]
    expect = moltype.degen_gapped_alphabet.array_to_bytes(raw)
    assert got == expect


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_seqs_data_seq_lengths(seqid, aligned_dict, dna_alphabet, dna_make_seq):
    """AlignedSeqsData.seq_lengths should return the length of the ungapped sequences"""
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    got = ad.seq_lengths()[seqid]
    expect = len(aligned_dict[seqid].replace("-", ""))
    assert got == expect


def test_aligned_seqs_data_add_seqs(dna_alphabet, dna_make_seq):
    data = {"seq1": "ACGT-", "seq2": "ACG-T", "seq3": "A---T"}
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=data, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    new_data = {"seq4": "ACGTT", "seq5": "ACG--", "seq6": "-----"}
    new_ad = ad.add_seqs(new_data)
    assert new_ad.names == ["seq1", "seq2", "seq3", "seq4", "seq5", "seq6"]


def test_aligned_seqs_data_add_seqs_diff_lengths_raises(dna_alphabet, dna_make_seq):
    """adding sequences of different lengths should raise an error"""

    data = {"seq1": "ACGT-", "seq2": "ACG-T"}
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=data, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    new_data = {"seq3": "AC", "seq4": "A-"}
    with pytest.raises(ValueError):
        _ = ad.add_seqs(new_data)

    new_data = {"seq3": "AAA", "seq4": "A-"}
    with pytest.raises(ValueError):
        _ = ad.add_seqs(new_data)


def test_aligned_seqs_data_add_seqs_diff_moltype_raises(dna_alphabet, dna_make_seq):
    """adding sequences of different moltype should raise an error"""
    data = {"seq1": "ACGT-", "seq2": "ACG-T"}  # DNA
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=data, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    new_data = {"seq3": "ACGU-", "seq4": "ACG-U"}  # RNA
    with pytest.raises(ValueError):
        _ = ad.add_seqs(new_data)


def test_aligned_seqs_data_add_seqs_duplicate_keys_raises(dna_alphabet, dna_make_seq):
    """AlignedSeqsData.add_seqs should raise an error if their are duplicated
    sequence names"""
    data = {"seq1": "ACGT-", "seq2": "ACG-T"}
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=data, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    new_data = {"seq2": "ACGT-", "seq3": "ACT-T"}  # seq2 is duplicated
    with pytest.raises(ValueError):
        _ = ad.add_seqs(new_data)

    # if we set force_unique_keys=False, it should not raise an error
    new_ad = ad.add_seqs(new_data, force_unique_keys=False)
    assert new_ad.names == ["seq1", "seq2", "seq3"]


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_seqs_data_get_aligned_view(
    aligned_dict, seqid, dna_alphabet, dna_make_seq
):
    # str on an ADV should return the ungapped sequence
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    got = ad.get_view(seqid)
    assert got.parent == ad
    assert got.parent_len == ad.align_len
    assert str(got) == aligned_dict[seqid].replace("-", "")


def test_aligned_seqs_data_subset_raises(aligned_dict, dna_alphabet, dna_make_seq):
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    with pytest.raises(ValueError):
        _ = ad.subset(["seq99"])


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_data_view_array(aligned_array_dict, dna_alphabet, dna_make_seq, seqid):
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_array_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    view = ad.get_view(seqid)
    got = numpy.array(view)
    expect = aligned_array_dict[seqid][aligned_array_dict[seqid] != 4]  # remove gaps
    assert numpy.array_equal(got, expect)

    # directly accessing .array_value property should return the same result
    got = view.array_value
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_data_view_gapped_array(
    aligned_array_dict, dna_alphabet, dna_make_seq, seqid
):
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_array_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    view = ad.get_view(seqid)
    got = view.gapped_array_value
    expect = aligned_array_dict[seqid]
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_data_view_str(aligned_dict, dna_alphabet, dna_make_seq, seqid):
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    view = ad.get_view(seqid)
    got = str(view)
    expect = aligned_dict[seqid].replace("-", "")
    assert got == expect

    # directly accessing .str_value property should return the same result
    got = view.str_value
    assert got == expect


def test_aligned_data_view_gapped_str_value(aligned_dict, dna_alphabet, dna_make_seq):
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    view = ad.get_view("seq1")
    got = view.gapped_str_value
    expect = aligned_dict["seq1"]
    assert got == expect


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_data_view_bytes(aligned_array_dict, dna_alphabet, dna_make_seq, seqid):
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_array_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    view = ad.get_view(seqid)
    got = bytes(view)
    expect = aligned_array_dict[seqid][aligned_array_dict[seqid] != 4]  # remove gaps
    expect = dna_alphabet.array_to_bytes(expect)  # convert to bytes
    assert numpy.array_equal(got, expect)

    # directly accessing .bytes_value property should return the same result
    got = view.bytes_value
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_data_view_gapped_bytes_value(
    aligned_array_dict, dna_alphabet, dna_make_seq, seqid
):
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_array_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    view = ad.get_view(seqid)
    got = view.gapped_bytes_value
    expect = dna_alphabet.array_to_bytes(aligned_array_dict[seqid])
    assert numpy.array_equal(got, expect)


def test_alignment_init(aligned_dict, dna_moltype, dna_alphabet, dna_make_seq):
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    aln = new_alignment.Alignment(seqs_data=ad, moltype=dna_moltype)
    assert aln.moltype == dna_moltype


def test_make_aligned_seqs_dict(aligned_dict):
    aln = new_alignment.make_aligned_seqs(aligned_dict, moltype="dna")
    assert isinstance(aln, new_alignment.Alignment)
    # if we index a seq, it should be an Aligned instance
    assert isinstance(aln.seqs["seq1"], new_alignment.Aligned)
    # if we use .get_seq, it should be a Sequence instance
    assert isinstance(aln.get_seq("seq1"), new_sequence.Sequence)
    assert aln.names == ["seq1", "seq2", "seq3", "seq4"]
    assert aln.to_dict() == aligned_dict


def test_make_aligned_seqs_aligned_seqs_data(aligned_dict, dna_alphabet, dna_make_seq):
    aligned_seqs_data = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    aln = new_alignment.make_aligned_seqs(aligned_seqs_data, moltype="dna")
    assert isinstance(aln, new_alignment.Alignment)
    assert aln.moltype.label == "dna"
    assert aln.names == ["seq1", "seq2", "seq3", "seq4"]
    assert aligned_dict == aln.to_dict()


def test_make_aligned_seqs_incompatible_moltype(
    aligned_dict, dna_alphabet, dna_make_seq
):
    ad = new_alignment.AlignedSeqsData.from_aligned_seqs(
        data=aligned_dict, alphabet=dna_alphabet, make_seq=dna_make_seq
    )
    with pytest.raises(ValueError):
        _ = new_alignment.make_aligned_seqs(ad, moltype="rna")


def test_make_aligned_seqs_incompatible_type():
    with pytest.raises(NotImplementedError):
        _ = new_alignment.make_aligned_seqs("TGCA", moltype="dna")


def test_alignment_get_gapped_seq():
    """Alignment.get_gapped_seq should return seq, with gaps"""
    aln = new_alignment.make_aligned_seqs(
        {"seq1": "--TTT?", "seq2": "GATC??"}, moltype="dna"
    )
    got = aln.get_gapped_seq(seqname="seq1")
    expect = "--TTT?"
    assert got == expect


def test_alignment_iter_positions():
    data = {"a": "AAAAAA", "b": "AAA---", "c": "AAAA--"}
    r = new_alignment.make_aligned_seqs({k: data[k] for k in "cb"}, moltype="dna")
    assert list(r.iter_positions(pos_order=[5, 1, 3])) == list(
        map(list, ["--", "AA", "A-"])
    )
    # reorder names
    r = new_alignment.make_aligned_seqs(data, moltype="dna")
    cols = list(r.iter_positions())
    assert cols == list(map(list, ["AAA", "AAA", "AAA", "A-A", "A--", "A--"]))


def test_alignment_to_dict(gapped_seqs_dict):
    aln = new_alignment.make_aligned_seqs(gapped_seqs_dict, moltype="protein")
    assert aln.to_dict() == gapped_seqs_dict


def test_alignment_get_lengths(gapped_seqs_dict):
    aln = new_alignment.make_aligned_seqs(gapped_seqs_dict, moltype="dna")
    got = aln.get_lengths()
    expect = {name: len(seq.replace("-", "")) for name, seq in gapped_seqs_dict.items()}
    assert got == expect


@pytest.mark.parametrize(
    "moltype",
    ["rna", "dna", "protein"],
)
def test_upac_consensus_allow_gaps(moltype):
    aln = new_alignment.make_aligned_seqs(
        {"s1": "ACGG", "s2": "ACGG", "s3": "-CGG"},
        moltype=moltype,
    )
    # default behaviour
    iupac = aln.iupac_consensus()
    assert iupac == "?CGG"

    # allow_gaps
    iupac = aln.iupac_consensus(allow_gap=False)
    assert iupac == "ACGG"


def test_alignment_to_pretty():
    """produce correct pretty print formatted text"""
    seqs = {"seq1": "ACGAANGA", "seq2": "-CGAACGA", "seq3": "ATGAACGA"}
    expect = ["seq1    ACGAANGA", "seq2    -....C..", "seq3    .T...C.."]

    aln = new_alignment.make_aligned_seqs(seqs, moltype="dna")
    got = aln.to_pretty(name_order=["seq1", "seq2", "seq3"])
    assert got == "\n".join(expect)

    got = aln.to_pretty(name_order=["seq1", "seq2", "seq3"], wrap=4)
    expect = [
        "seq1    ACGA",
        "seq2    -...",
        "seq3    .T..",
        "",
        "seq1    ANGA",
        "seq2    .C..",
        "seq3    .C..",
    ]
    assert got == "\n".join(expect)


def test_alignment_to_html():
    """produce correct html formatted text"""
    seqs = {"seq1": "ACG", "seq2": "-CT"}

    aln = new_alignment.make_aligned_seqs(seqs, moltype="dna")
    got = aln.to_html(ref_name="longest")  # name_order=['seq1', 'seq2'])
    # ensure balanced tags are in the txt
    for tag in ["<style>", "</style>", "<div", "</div>", "<table>", "</table>"]:
        assert tag in got

    ref_row = (
        '<tr><td class="label">seq1</td>'
        '<td><span class="A_dna">A</span>'
        '<span class="C_dna">C</span>'
        '<span class="G_dna">G</span></td></tr>'
    )
    other_row = (
        '<tr><td class="label">seq2</td>'
        '<td><span class="ambig_dna">-</span>'
        '<span class="C_dna">.</span>'
        '<span class="T_dna">T</span></td></tr>'
    )

    assert ref_row in got
    assert other_row in got
    assert got.find(ref_row) < got.find(other_row)

    # using different ref sequence
    ref_row = (
        '<tr><td class="label">seq2</td>'
        '<td><span class="terminal_ambig_dna">-</span>'
        '<span class="C_dna">C</span>'
        '<span class="T_dna">T</span></td></tr>'
    )
    other_row = (
        '<tr><td class="label">seq1</td>'
        '<td><span class="A_dna">A</span>'
        '<span class="C_dna">.</span>'
        '<span class="G_dna">G</span></td></tr>'
    )
    got = aln.to_html(ref_name="seq2")
    # order now changes
    assert got.find(ref_row) < got.find(other_row)


def test_alignment_repr():
    data = {
        "ENSMUSG00000056468": "GCCAGGGGGAAAA",
        "ENSMUSG00000039616": "GCCCTTCAAATTT",
    }
    seqs = new_alignment.make_aligned_seqs(data, moltype="dna")
    assert (
        repr(seqs)
        == "2 x 13 dna alignment: ENSMUSG00000056468[GCCAGGGGGA...], ENSMUSG00000039616[GCCCTTCAAA...]"
    )

    data = {
        "ENSMUSG00000039616": "GCCCTTCAAATTT",
        "ENSMUSG00000056468": "GCCAGGGGGAAAA",
    }
    seqs = new_alignment.make_aligned_seqs(data, moltype="dna")

    assert (
        repr(seqs)
        == "2 x 13 dna alignment: ENSMUSG00000039616[GCCCTTCAAA...], ENSMUSG00000056468[GCCAGGGGGA...]"
    )

    data = {
        "a": "TCGAT",
    }
    seqs = new_alignment.make_aligned_seqs(data, moltype="dna")
    assert repr(seqs) == "1 x 5 dna alignment: a[TCGAT]"

    data = {
        "a": "TCGAT" * 2,
    }
    seqs = new_alignment.make_aligned_seqs(data, moltype="dna")
    assert repr(seqs) == "1 x 10 dna alignment: a[TCGATTCGAT]"

    data = {
        "a": "A" * 11,
        "b": "B" * 11,
        "c": "C" * 11,
        "d": "D" * 11,
        "e": "E" * 11,
    }
    seqs = new_alignment.make_aligned_seqs(data, moltype="text")
    assert (
        repr(seqs)
        == "5 x 11 text alignment: a[AAAAAAAAAA...], b[BBBBBBBBBB...], c[CCCCCCCCCC...], ..."
    )


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_alignment_getitem_slice(aligned_dict, seqid):
    """slicing an alignment should propogate the slice to aligned instances"""
    start = 1
    stop = 5
    step = 3
    aln = new_alignment.make_aligned_seqs(aligned_dict, moltype="dna")
    sliced_aln = aln[start:stop:step]
    got = sliced_aln.get_seq(seqid)
    expect = aligned_dict[seqid][start:stop:step].replace("-", "")
    assert got == expect


@pytest.mark.parametrize("index", (0, 1, 2, 3))
@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_alignment_getitem_int(aligned_dict, index, seqid):
    aln = new_alignment.make_aligned_seqs(aligned_dict, moltype="dna")
    got = aln[index].get_gapped_seq(seqid)
    expect = aligned_dict[seqid][index]
    assert got == expect


def test_alignment_getitem_raises(aligned_dict):
    aln = new_alignment.make_aligned_seqs(aligned_dict, moltype="dna")
    with pytest.raises(NotImplementedError):
        _ = aln[1.0]


@pytest.mark.parametrize("start", range(6))
@pytest.mark.parametrize("stop", range(6))
@pytest.mark.parametrize("step", range(1, 3))
@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_alignment_slice_pos_step_gapped(aligned_dict, start, stop, step, seqid):
    """slicing an alignment should propogate the slice to aligned instances"""
    aln = new_alignment.make_aligned_seqs(aligned_dict, moltype="dna")
    sliced_aln = aln[start:stop:step]
    got = sliced_aln.seqs[seqid].gapped_seq
    expect = aligned_dict[seqid][start:stop:step]
    assert got == expect


@pytest.mark.parametrize("start", range(6))
@pytest.mark.parametrize("stop", range(6))
@pytest.mark.parametrize("step", range(1, 3))
@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_alignment_slice_pos_step_ungapped(aligned_dict, start, stop, step, seqid):
    """slicing an alignment should propogate the slice to aligned instances"""
    aln = new_alignment.make_aligned_seqs(aligned_dict, moltype="dna")
    sliced_aln = aln[start:stop:step]
    got = sliced_aln.seqs[seqid].seq
    expect = aligned_dict[seqid][start:stop:step].replace("-", "")
    assert got == expect


@pytest.mark.parametrize("start", range(6))
@pytest.mark.parametrize("stop", range(6))
@pytest.mark.parametrize("step", range(-3, 0))
@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_alignment_slice_neg_step_gapped(aligned_dict, start, stop, step, seqid):
    """slicing an alignment should propogate the slice to aligned instances"""
    aln = new_alignment.make_aligned_seqs(aligned_dict, moltype="dna")
    sliced_aln = aln[start:stop:step]
    got = sliced_aln.seqs[seqid].data.gapped_str_value
    expect = aligned_dict[seqid][start:stop:step]
    assert got == expect


@pytest.mark.parametrize("start", range(6))
@pytest.mark.parametrize("stop", range(6))
@pytest.mark.parametrize("step", range(-3, 0))
@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_alignment_slice_neg_step_ungapped(aligned_dict, start, stop, step, seqid):
    """slicing an alignment should propogate the slice to aligned instances"""
    aln = new_alignment.make_aligned_seqs(aligned_dict, moltype="dna")
    dna = aln.moltype
    sliced_aln = aln[start:stop:step]
    got = sliced_aln.seqs[seqid].seq
    expect = dna.complement(aligned_dict[seqid][start:stop:step].replace("-", ""))
    assert got == expect


@pytest.mark.xfail(reason="not implemented yet")
@pytest.mark.parametrize(
    "coll_maker", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_alignment_strand_invalid(aligned_dict, coll_maker):
    with pytest.raises(ValueError):
        _ = coll_maker(aligned_dict, moltype="dna", strand={"seq1": "invalid"})

    with pytest.raises(ValueError):
        _ = coll_maker(aligned_dict, moltype="dna", strand={"seq1": 3})


@pytest.mark.xfail(reason="design of strand attribute is not finalized")
@pytest.mark.parametrize(
    "coll_maker", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_alignment_strand(aligned_dict, coll_maker):
    aln = coll_maker(aligned_dict, moltype="dna", strand={"seq1": -1})
    assert aln.seqs.strand["seq1"] == -1
    # if we don't set the strand information, it should be set to 1
    assert aln.seqs.strand["seq2"] == 1
    # reverse complement should flip the strand orientation
    rc_aln = aln.rc()
    assert rc_aln.seqs.strand["seq1"] == 1
    assert rc_aln.seqs.strand["seq2"] == -1


@pytest.mark.xfail(reason="not implemented for new style alignments")
@pytest.mark.parametrize("method", ("ic_score", "cogent3_score", "sp_score"))
def test_alignment_quality_methods(method):
    data = {
        "DogFaced": "TG----AATATGT------GAAAGAG",
        "FreeTaile": "TTGAAGAATATGT------GAAAGAG",
        "LittleBro": "CTGAAGAACCTGTGAAAGTGAAAGAG",
    }
    expected_score = dict(
        cogent3_score=-123.0, ic_score=32.93032499, sp_score=-25.25687109
    )[method]
    aln = new_alignment.make_aligned_seqs(
        data,
        moltype="dna",
        info=dict(align_params=dict(lnL=-123.0)),
    )
    app = get_app(method)
    score = app(aln)
    assert numpy.allclose(score, expected_score)


def test_alignment_counts_per_pos():
    """correctly count motifs"""
    exp = numpy.array(
        [
            [1, 1, 1, 0],
            [0, 2, 0, 1],
            [0, 0, 3, 0],
            [1, 1, 0, 1],
            [0, 0, 3, 0],
            [1, 1, 0, 1],
        ]
    )

    exp_gap = numpy.array(
        [
            [1, 1, 0, 1, 0],
            [0, 2, 0, 0, 1],
            [0, 0, 3, 0, 0],
            [0, 2, 0, 1, 0],
            [0, 1, 2, 0, 0],
            [0, 2, 0, 1, 0],
        ]
    )

    aln = new_alignment.make_aligned_seqs(
        {"s1": "TCAGAG", "s2": "CCACAC", "s3": "AGATAT"}, moltype="dna"
    )
    obs = aln.counts_per_pos()
    assert numpy.array_equal(obs.array, exp)
    assert numpy.array_equal(obs.motifs, tuple(aln.moltype.alphabet))
    obs = aln.counts_per_pos(motif_length=2)
    assert numpy.array_equal(obs[0, "TC"], 1)
    assert numpy.array_equal(obs[1, "AC"], 1)
    assert numpy.array_equal(obs[2, "AC"], 1)
    aln = new_alignment.make_aligned_seqs(
        {"s1": "TCAGAG", "s2": "CCACAC", "s4": "G-ACCC"}, moltype="dna"
    )
    obs = aln.counts_per_pos(allow_gap=True)
    assert numpy.array_equal(obs.array, exp_gap)
    aln = new_alignment.make_aligned_seqs(
        {"s1": "-RAT", "s2": "ACCT", "s3": "GTGT"}, moltype="dna"
    )
    c = aln.counts_per_pos(include_ambiguity=False, allow_gap=True)
    assert numpy.array_equal(set(c.motifs), set("ACGT-"))


def test_alignment_get_gap_array():
    data = {
        "DogFaced": "TG-",
        "FreeTaile": "T-G",
        "LittleBro": "---",
        "BigBro": "GTT",
        "MiddleBro": "-T-",
    }
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    got = aln.get_gap_array()
    expect = numpy.array(
        [
            [False, False, True],
            [False, True, False],
            [True, True, True],
            [False, False, False],
            [True, False, True],
        ],
        dtype=bool,
    )
    assert numpy.allclose(got, expect)


def test_get_position_indices():
    """get_position_indices should return names of cols where f(col)"""

    def gap_1st(x):
        return x[0] == "-"

    def gap_2nd(x):
        return x[1] == "-"

    def gap_3rd(x):
        return x[2] == "-"

    def is_list(x):
        return isinstance(x, list)

    gaps = new_alignment.make_aligned_seqs(
        {"a": "AAAAAAA", "b": "A--A-AA", "c": "AA-----"}, moltype="dna"
    )

    assert gaps.get_position_indices(gap_1st) == []
    assert gaps.get_position_indices(gap_2nd) == [1, 2, 4]
    assert gaps.get_position_indices(gap_3rd) == [2, 3, 4, 5, 6]
    assert gaps.get_position_indices(is_list) == [0, 1, 2, 3, 4, 5, 6]
    # should be able to negate
    assert gaps.get_position_indices(gap_2nd, negate=True) == [0, 3, 5, 6]
    assert gaps.get_position_indices(gap_1st, negate=True) == [0, 1, 2, 3, 4, 5, 6]
    assert gaps.get_position_indices(is_list, negate=True) == []


def test_iupac_consensus_rna():
    """Alignment iupac_consensus should use RNA IUPAC symbols correctly"""
    aln = new_alignment.make_aligned_seqs(
        {
            "seq1": "UCAGN-UCAGN-UCAGN-UCAGAGCAUN-",
            "seq2": "UUCCAAGGNN--UUCCAAGGNNAGCAG--",
            "seq3": "UUCCAAGGNN--UUCCAAGGNNAGCUA--",
            "seq4": "UUUUCCCCAAAAGGGGNNNN--AGCUA--",
            "seq5": "UUUUCCCCAAAAGGGGNNNN--AGCUA--",
        },
        moltype="rna",
    )

    # following IUPAC consensus calculated by hand
    assert aln.iupac_consensus() == "UYHBN?BSNN??KBVSN?NN??AGCWD?-"


def test_iupac_consensus_dna():
    """Alignment iupac_consensus should use DNA IUPAC symbols correctly"""
    aln = new_alignment.make_aligned_seqs(
        {
            "seq1": "TCAGN-TCAGN-TCAGN-TCAGAGCATN-",
            "seq2": "TTCCAAGGNN--TTCCAAGGNNAGCAG--",
            "seq3": "TTCCAAGGNN--TTCCAAGGNNAGCTA--",
            "seq4": "TTTTCCCCAAAAGGGGNNNN--AGCTA--",
            "seq5": "TTTTCCCCAAAAGGGGNNNN--AGCTA--",
        },
        moltype="dna",
    )
    # following IUPAC consensus calculated by hand
    assert aln.iupac_consensus() == "TYHBN?BSNN??KBVSN?NN??AGCWD?-"


def test_iupac_consensus_protein():
    """Alignment iupac_consensus should use protein IUPAC symbols correctly"""
    aln = new_alignment.make_aligned_seqs(
        {
            "seq1": "ACDEFGHIKLMNPQRSTUVWY-",
            "seq2": "ACDEFGHIKLMNPQRSUUVWF-",
            "seq3": "ACDEFGHIKLMNPERSKUVWC-",
            "seq4": "ACNEFGHIKLMNPQRS-UVWP-",
        },
        moltype="protein",
    )
    # following IUPAC consensus calculated by hand
    # Test all uppper
    assert aln.iupac_consensus() == "ACBEFGHIKLMNPZRS?UVWX-"


def test_majority_consensus():
    """Alignment.majority_consensus should return commonest symbol per column"""
    # Check the exact strings expected from string transform
    aln = new_alignment.make_aligned_seqs(
        {
            "seq1": "ACG",
            "seq2": "ACG",
            "seq3": "TTT",
        },
        moltype="dna",
    )
    assert aln.majority_consensus() == "ACG"


def test_probs_per_pos():
    """Alignment.probs_per_pos should find Pr(symbol) in each
    column"""
    # 4 seqs (easy to calculate probabilities)
    align = new_alignment.make_aligned_seqs(
        {"seq1": "AAA", "seq2": "ACA", "seq3": "GGG", "seq4": "GUC"}, moltype="rna"
    )
    got = align.probs_per_pos()
    # check that the column probs match the counts we expect
    expect = [
        {"A": 0.5, "G": 0.5},
        {"A": 0.25, "C": 0.25, "G": 0.25, "U": 0.25},
        {"A": 0.5, "G": 0.25, "C": 0.25},
    ]
    for pos, probs in enumerate(expect):
        for char, prob in probs.items():
            assert numpy.allclose(got[pos, char], prob)


def test_entropy_per_pos():
    """SequenceCollection.entropy_per_pos should match hand-calculated values"""
    aln = new_alignment.make_aligned_seqs({"seq1": "ATA", "seq2": "AAA"}, moltype="dna")
    obs = aln.entropy_per_pos()
    assert numpy.allclose(obs, [0, 1, 0])
    # check what happens with only one input sequence
    aln = new_alignment.make_aligned_seqs({"seq1": "TGC"}, moltype="dna")
    obs = aln.entropy_per_pos()
    assert numpy.allclose(obs, [0, 0, 0])


@pytest.mark.parametrize(
    "coll_maker", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_entropy_excluding_unobserved(coll_maker):
    """omitting unobserved motifs should not affect entropy calculation"""
    a = coll_maker(dict(a="ACAGGG", b="AGACCC", c="GGCCTA"), moltype="dna")
    entropy_excluded = a.entropy_per_seq(exclude_unobserved=True)
    entropy_unexcluded = a.entropy_per_seq(exclude_unobserved=False)
    assert numpy.allclose(entropy_excluded, entropy_unexcluded)


def test_seq_entropy_just_gaps():
    """get_seq_entropy should get entropy of each seq"""
    aln = new_alignment.make_aligned_seqs(dict(a="A---", b="----"), moltype="dna")
    got = aln.entropy_per_seq()
    expect = numpy.array([0, numpy.nan], dtype=numpy.float64)

    assert numpy.allclose(got, expect, equal_nan=True)

    aln = new_alignment.make_aligned_seqs(dict(a="----", b="----"), moltype="dna")
    entropy = aln.entropy_per_seq()
    assert entropy is None


def test_get_gap_array():
    aln = new_alignment.make_aligned_seqs(
        {"seq1": "A-GN", "seq2": "TG--", "seq3": "----"}, moltype="dna"
    )
    got = aln.get_gap_array()
    expect = numpy.array(
        [
            [False, True, False, True],
            [False, False, True, True],
            [True, True, True, True],
        ]
    )
    assert numpy.allclose(got, expect)

    got = aln.get_gap_array(include_ambiguity=False)
    expect = numpy.array(
        [
            [False, True, False, False],
            [False, False, True, True],
            [True, True, True, True],
        ]
    )


def test_count_gaps_per_pos():
    """correctly compute the number of gaps"""
    data = {"a": "AAAA---GGT", "b": "CCC--GG?GT"}
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    # per position
    got = aln.count_gaps_per_pos(include_ambiguity=False)
    assert numpy.array_equal(got.array, [0, 0, 0, 1, 2, 1, 1, 0, 0, 0])
    got = aln.count_gaps_per_pos(include_ambiguity=True)
    assert numpy.array_equal(got.array, [0, 0, 0, 1, 2, 1, 1, 1, 0, 0])


def test_count_gaps_per_seq():
    """correctly compute the number of gaps"""
    data = {"a": "AAAA---GGT", "b": "CCC--GG?GT"}
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    got = aln.count_gaps_per_seq(include_ambiguity=False)
    assert numpy.array_equal(got.array, [3, 2])
    assert numpy.array_equal(got["b"], 2)
    got = aln.count_gaps_per_seq(include_ambiguity=True)
    assert numpy.array_equal(got.array, [3, 3])
    assert numpy.array_equal(got["b"], 3)
    # per seq, unique
    got = aln.count_gaps_per_seq(include_ambiguity=False, unique=True)
    assert numpy.array_equal(got.array, [1, 2])
    got = aln.count_gaps_per_seq(include_ambiguity=True, unique=True)
    assert numpy.array_equal(got.array, [2, 2])

    data = {"a": "AAAGGG", "b": "------", "c": "------"}
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    got = aln.count_gaps_per_seq(include_ambiguity=False, unique=True)
    assert numpy.array_equal(got.array, [6, 0, 0])
    assert numpy.array_equal(got["a"], 6)
    assert numpy.array_equal(got["b"], 0)

    # per_seq, induced_by
    data = {"a": "--ACGT---GTAC", "b": "--ACGTA--GT--", "c": "--ACGTA-AGT--"}
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    got = aln.count_gaps_per_seq(unique=False, induced_by=True)
    assert numpy.array_equal(got.array, [2, 1, 2])
    assert numpy.array_equal(got["b"], 1)


def test_omit_bad_seqs():
    """omit_bad_seqs should return alignment w/o seqs causing most gaps"""
    data = {
        "s1": "---ACC---TT-",
        "s2": "---ACC---TT-",
        "s3": "---ACC---TT-",
        "s4": "--AACCG-GTT-",
        "s5": "--AACCGGGTTT",
        "s6": "AGAACCGGGTT-",
    }

    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    # with defaults, excludes s6
    expect = data.copy()
    del expect["s6"]
    result = aln.omit_bad_seqs()
    assert result.to_dict() == expect
    # with quantile 0.5, just s1, s2, s3
    expect = data.copy()
    for key in ("s6", "s5"):
        del expect[key]
    result = aln.omit_bad_seqs(0.5)
    assert result.to_dict() == expect


def test_matching_ref():
    """Alignment.matching_ref returns new aln with well-aln to temp"""
    data = {
        "s1": "UC-----CU---C",
        "s2": "UC------U---C",
        "s3": "UUCCUUCUU-UUC",
        "s4": "UU-UUUU-UUUUC",
        "s5": "-------------",
    }

    aln = new_alignment.make_aligned_seqs(data, moltype="rna")
    result = aln.matching_ref("s3", 0.9, 5)
    assert result.to_dict() == {"s3": "UUCCUUCUU-UUC", "s4": "UU-UUUU-UUUUC"}
    result2 = aln.matching_ref("s4", 0.9, 4)
    assert result2.to_dict() == {"s3": "UUCCUUCUU-UUC", "s4": "UU-UUUU-UUUUC"}
    result3 = aln.matching_ref("s1", 0.9, 4)
    assert result3.to_dict() == {
        "s2": "UC------U---C",
        "s1": "UC-----CU---C",
        "s5": "-------------",
    }

    result4 = aln.matching_ref("s3", 0.5, 13)
    assert result4.to_dict() == {"s3": "UUCCUUCUU-UUC", "s4": "UU-UUUU-UUUUC"}


def test_sliding_windows():
    """sliding_windows should return slices of alignments."""
    alignment = new_alignment.make_aligned_seqs(
        {"seq1": "ACGTACGT", "seq2": "ACGTACGT", "seq3": "ACGTACGT"}, moltype="dna"
    )
    result = []

    for bit in alignment.sliding_windows(5, 2):
        result += [bit]
    assert result[0].to_dict() == {"seq3": "ACGTA", "seq2": "ACGTA", "seq1": "ACGTA"}
    assert result[1].to_dict() == {"seq3": "GTACG", "seq2": "GTACG", "seq1": "GTACG"}

    result = []
    for bit in alignment.sliding_windows(5, 1):
        result += [bit]
    assert result[0].to_dict() == {"seq3": "ACGTA", "seq2": "ACGTA", "seq1": "ACGTA"}
    assert result[1].to_dict() == {"seq3": "CGTAC", "seq2": "CGTAC", "seq1": "CGTAC"}
    assert result[2].to_dict() == {"seq3": "GTACG", "seq2": "GTACG", "seq1": "GTACG"}
    assert result[3].to_dict() == {"seq3": "TACGT", "seq2": "TACGT", "seq1": "TACGT"}


@pytest.mark.xfail(reason="todo: get_features not properly implemented so test hangs")
def test_get_feature():
    aln = new_alignment.make_aligned_seqs(
        {"x": "-AAAAAAAAA", "y": "TTTT--CCCT"}, moltype="dna"
    )

    db = aln.annotation_db
    db.add_feature(seqid="y", biotype="exon", name="A", spans=[(5, 8)])
    feat = list(aln.get_features(seqid="y", biotype="exon", on_alignment=False))[0]
    assert feat.get_slice().to_dict() == dict(x="AAA", y="CCT")


@pytest.fixture
def alignment():
    data = {
        "seq1": "A-GT",
        #         012
        "seq2": "-GGT",
    }
    return new_alignment.make_aligned_seqs(data, moltype="dna")


def test_aligned_view_parent_coords(alignment):
    seqid = "seq2"
    a1 = alignment.seqs[seqid]
    got = a1.parent_coordinates()
    assert got == (seqid, 0, 4, 1)

    a2 = alignment[2:]
    a1_2 = a2.seqs[seqid]
    assert a1_2.parent_coordinates() == (seqid, 2, 4, 1)

    # now getting the sequence coordinates
    s2 = a2.get_seq(seqid)
    expect = s2.parent_coordinates()
    got = a1_2.parent_coordinates(seq_coords=True)
    assert got == expect


def test_aligned_view_parent_coords_reversed(alignment):
    seqid = "seq2"
    a2 = alignment[2:].rc()
    a1_2 = a2.seqs[seqid]
    assert a1_2.parent_coordinates() == (seqid, 2, 4, -1)

    s2 = a2.get_seq(seqid)
    expect = s2.parent_coordinates()
    got = a1_2.parent_coordinates(seq_coords=True)
    assert got == expect


@pytest.mark.parametrize("rced", [True, False])
def test_get_seq_from_slice(alignment, rced):
    seqid = "seq2"
    raw = str(alignment.seqs[seqid])
    dna = alignment.moltype
    expect = dna.rc(raw[2:]) if rced else raw[2:]

    a2 = alignment[2:].rc() if rced else alignment[2:]
    seq = a2.get_seq(seqid)
    got = str(seq)
    assert got == expect


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
def test_alignment_indexing_string(alignment, seqid):
    # when indexing with a string, should return an Aligned instance
    got = alignment.seqs[seqid]
    assert isinstance(got, new_alignment.Aligned)
    assert str(got) == alignment.to_dict()[seqid]
