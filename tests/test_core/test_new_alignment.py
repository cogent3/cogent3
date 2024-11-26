import json
import os
import pathlib
import re
from warnings import catch_warnings, filterwarnings

import numpy
import pytest

from cogent3 import get_app, load_aligned_seqs, load_unaligned_seqs, open_
from cogent3._version import __version__
from cogent3.core import new_alignment, new_alphabet, new_moltype, new_sequence
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import GffAnnotationDb, load_annotations
from cogent3.core.location import FeatureMap, LostSpan, Span
from cogent3.util.deserialise import deserialise_object
from cogent3.util.misc import get_object_provenance


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
def rna_alphabet():
    moltype = new_moltype.get_moltype("rna")
    return moltype.degen_gapped_alphabet


@pytest.fixture
def rna_moltype():
    return new_moltype.get_moltype("rna")


@pytest.fixture
def dna_sd(str_seqs_dict: dict[str, str], dna_alphabet):
    return new_alignment.SeqsData.from_seqs(data=str_seqs_dict, alphabet=dna_alphabet)


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
    return new_alignment.make_unaligned_seqs(
        {"a": "AAAAA", "c": "CCCCC"}, moltype="dna"
    )


@pytest.fixture
def ordered2():
    return new_alignment.make_unaligned_seqs(
        {"c": "CCCCC", "a": "AAAAA"}, moltype="dna"
    )


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


def test_seqs_data_construction(str_seqs_dict, dna_alphabet):
    """SeqsData can be constructed from a dict and alphabet, either directly or via from_seqs"""
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=dna_alphabet)
    assert isinstance(sd, new_alignment.SeqsData)

    sd = new_alignment.SeqsData.from_seqs(data=str_seqs_dict, alphabet=dna_alphabet)
    assert isinstance(sd, new_alignment.SeqsData)


def test_seqs_data_construction_wrong_alphabet(str_seqs_dict, rna_alphabet):
    """SeqsData should raise ValueError if alphabet is incompatible with data"""
    with pytest.raises(new_alphabet.AlphabetError):
        _ = new_alignment.SeqsData(data=str_seqs_dict, alphabet=rna_alphabet)


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
    sd = new_alignment.SeqsData.from_seqs(data=d, alphabet=dna_alphabet)
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
@pytest.mark.parametrize("start", (None, 0, 1))
@pytest.mark.parametrize("stop", (None, 3, 4))
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


@pytest.fixture(params=["str_seqs_dict", "arr_seqs_dict"])
def seqs_dicts(request):
    return request.getfixturevalue(request.param)


@pytest.mark.parametrize("seqid", ["seq1", "seq2", "seq3"])
def test_seqs_data_seq_lengths(seqs_dicts, dna_alphabet, seqid):
    sd = new_alignment.SeqsData(data=seqs_dicts, alphabet=dna_alphabet)
    expect = len(seqs_dicts[seqid])
    got = sd.get_seq_length(seqid)
    assert got == expect


@pytest.mark.parametrize("seqid", ["seq1", "seq2", "seq3"])
def test_seqs_data_get_seq_array(arr_seqs_dict, dna_alphabet, seqid):
    expect = arr_seqs_dict[seqid]
    sd = new_alignment.SeqsData(data=arr_seqs_dict, alphabet=dna_alphabet)
    got = sd.get_seq_array(seqid=seqid)
    assert numpy.array_equal(got, expect)


def test_seqs_data_get_seq_bytes(dna_sd: new_alignment.SeqsData):
    # getting seqid and slicing tested in test_get_seq_str
    got = dna_sd.get_seq_bytes(seqid="seq1")
    assert isinstance(got, bytes)


@pytest.mark.parametrize("seq", ("seq1", "seq2"))
def test_seqs_data_getitem_str(dna_sd, seq):
    got = dna_sd[seq]
    assert isinstance(got, new_alignment.SeqDataView)
    assert got.parent == dna_sd
    assert got.seqid == seq


@pytest.mark.parametrize("idx", (0, 1))
def test_seqs_data_getitem_int(str_seqs_dict, dna_sd, idx):
    got = dna_sd[idx]
    assert isinstance(got, new_alignment.SeqDataView)
    assert got.parent == dna_sd
    assert got.seqid == list(str_seqs_dict)[idx]


def test_seqs_data_getitem_raises(dna_sd):
    invalid_index = ["this", "shouldn't", "work"]
    with pytest.raises(NotImplementedError):
        _ = dna_sd[invalid_index]


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
    with pytest.raises(new_moltype.MolTypeError):
        _ = seqs.to_alphabet(DNA)


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
            "parent_len": len(sdv.str_value),
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


# tests of SequenceCollection and Alignment constructor utility function


@pytest.mark.parametrize(
    "mk_cls, cls",
    [
        (new_alignment.make_unaligned_seqs, new_alignment.SequenceCollection),
        (new_alignment.make_aligned_seqs, new_alignment.Alignment),
    ],
)
@pytest.mark.parametrize("moltype", ("dna", "rna", "protein", "protein_with_stop"))
def test_make_seqs(moltype, mk_cls, cls):
    """SequenceCollection and Alignment constructor functions should handle
    dict/list/set as input"""
    # dict of strings
    data = {"a": "AGGCCC", "b": "AGAAAA"}
    got = mk_cls(data, moltype=moltype)
    assert isinstance(got, cls)

    # dict of arrays
    data = {"a": numpy.array([0, 2, 1, 3, 2, 3]), "b": numpy.array([0, 2, 1, 3, 2, 1])}
    got = mk_cls(data, moltype=moltype)
    assert isinstance(got, cls)

    # list of str
    data = ["AGGCCC", "AGAAAA"]
    got = mk_cls(data, moltype=moltype)
    assert isinstance(got, cls)
    assert got.num_seqs == 2

    # set of str
    data = {"AGGCCC", "AGAAAA"}
    got = mk_cls(data, moltype=moltype)
    assert isinstance(got, cls)
    assert got.num_seqs == 2

    # list of array
    data = [numpy.array([0, 2, 1, 3, 2, 3]), numpy.array([0, 2, 1, 3, 2, 1])]
    got = mk_cls(data, moltype=moltype)
    assert isinstance(got, cls)

    # list(name: seq) pairs
    data = [["seq1", "AGGCCC"], ["seq2", "AGAAAA"]]
    got = mk_cls(data, moltype=moltype)
    assert isinstance(got, cls)

    # tuples
    data = [("seq1", "AGGCCC"), ("seq2", "AGAAAA")]
    got = mk_cls(data, moltype=moltype)
    assert isinstance(got, cls)


@pytest.mark.parametrize(
    "mk_cls, cls",
    [
        (new_alignment.make_unaligned_seqs, new_alignment.SeqsData),
        (new_alignment.make_aligned_seqs, new_alignment.AlignedSeqsData),
    ],
)
def test_make_seqs_label_to_name(mk_cls, cls, dna_alphabet):
    """SequenceCollection and Alignment constructor functions should convert names
    using label_to_name function if provided"""

    def f(x):
        return x.upper()

    # for data as a dict
    data = {"a": "AGGCCC", "b": "AGAAAA"}
    got = mk_cls(data, moltype="dna", label_to_name=f)
    assert list(got.names) == ["A", "B"]

    # for data as a SeqsData object
    data = cls.from_seqs(data=data, alphabet=dna_alphabet)
    got = mk_cls(data, moltype="dna", label_to_name=f)
    assert list(got.names) == ["A", "B"]


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_make_seqs_raises(mk_cls):
    """cannot construct SequenceCollection or Alignment from a string"""
    data = "AGTCCTGA"
    with pytest.raises(NotImplementedError):
        mk_cls(data, moltype="dna")


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_make_seqs_no_seqs(mk_cls):
    """cannot construct SequenceCollection or Alignment from an empty dict"""
    data = {}
    with pytest.raises(ValueError):
        mk_cls(data, moltype="dna")


@pytest.mark.parametrize(
    "mk_cls, seqs_data_cls",
    [
        (new_alignment.make_unaligned_seqs, new_alignment.SeqsData),
        (new_alignment.make_aligned_seqs, new_alignment.AlignedSeqsData),
    ],
)
def test_make_seqs_incompatible_moltype(mk_cls, seqs_data_cls, dna_alphabet):
    """SequenceCollection and Alignment constructor functions should raise an error
    if the provided moltype is incompatible with the alphabet of the seqs data"""

    data = {"a": "AGGCCC", "b": "AGAAAA"}
    seqs_data = seqs_data_cls.from_seqs(data=data, alphabet=dna_alphabet)

    with pytest.raises(ValueError):
        _ = mk_cls(seqs_data, moltype="rna")


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_make_seqs_from_list_generates_correct_names(mk_cls):
    """SequenceCollection/Alignment init from list of sequences should use indices as keys"""
    seqs = ["TTTTT", "CCCCC", "GGGGG"]
    a = mk_cls(seqs, moltype="dna")
    assert str(a.seqs["seq_0"]) == "TTTTT"
    assert str(a.seqs["seq_1"]) == "CCCCC"
    assert str(a.seqs["seq_2"]) == "GGGGG"
    assert a.names == ["seq_0", "seq_1", "seq_2"]


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_make_seqs_from_pairs(mk_cls):
    """SequenceCollection/Alignment init from list of (key,val) pairs should work correctly"""
    seqs = [["a", "AAA"], ["t", "TTT"], ["c", "CCC"]]
    a = mk_cls(seqs, moltype="dna")
    assert str(a.seqs["a"]) == "AAA"
    assert str(a.seqs["t"]) == "TTT"
    assert str(a.seqs["c"]) == "CCC"
    assert a.names == ["a", "t", "c"]


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_make_seqs_from_sequences(mk_cls):
    """SequenceCollection and Alignment constructor functions can be provided with
    a list of Sequence objects"""
    # if no names, they should be generated via the standard naming convention,
    # seq_0, seq_1, etc.
    seq1 = new_moltype.DNA.make_seq(seq="AC")
    seq2 = new_moltype.DNA.make_seq(seq="AC")
    coll = mk_cls([seq1, seq2], moltype="dna")
    assert isinstance(coll, new_alignment.SequenceCollection)
    assert coll.names == ["seq_0", "seq_1"]

    # if the sequences have names, they should be used
    seq1 = new_moltype.DNA.make_seq(seq="AC", name="seq1")
    seq2 = new_moltype.DNA.make_seq(seq="AC", name="seq2")
    coll = mk_cls([seq1, seq2], moltype="dna")
    assert isinstance(coll, new_alignment.SequenceCollection)
    assert coll.names == ["seq1", "seq2"]

    # if the data dict has different names to the seq names,
    # the names from data should be used
    coll = mk_cls({"s1": seq1, "s2": seq2}, moltype="dna")
    assert isinstance(coll, new_alignment.SequenceCollection)
    assert coll.names == ["s1", "s2"]
    assert coll.get_seq("s1") == seq1


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
@pytest.mark.parametrize(
    "seq_name, parent_name", [("seq_1", "parent_1"), ("seq_2", "parent_2")]
)
def test_make_seqs_renamed_seqs(mk_cls, seq_name, parent_name, dna_alphabet):
    # the parent_name should persist from sequence object to SequenceCollection
    seq_view_1 = new_sequence.SeqView(
        parent="AAAA", parent_len=4, seqid="parent_1", alphabet=dna_alphabet
    )
    seq_view_2 = new_sequence.SeqView(
        parent="TTTT", parent_len=4, seqid="parent_2", alphabet=dna_alphabet
    )

    seq_1 = new_moltype.DNA.make_seq(seq=seq_view_1, name="seq_1")
    seq_2 = new_moltype.DNA.make_seq(seq=seq_view_2, name="seq_2")

    seqs = mk_cls([seq_1, seq_2], moltype="dna")
    assert set(seqs.names) == {"seq_1", "seq_2"}
    assert seqs.get_seq(seq_name)._seq.seqid == parent_name


@pytest.mark.parametrize(
    "mk_cls, data_cls",
    [
        (new_alignment.make_unaligned_seqs, new_alignment.SeqsData),
        (new_alignment.make_aligned_seqs, new_alignment.AlignedSeqsData),
    ],
)
@pytest.mark.parametrize("seq", ("a", "b"))
def test_make_seqs_offset(mk_cls, data_cls, seq, dna_alphabet):
    """SequenceCollection and Alignment constructor functions should handle
    offset argument"""
    data = {"a": "AGGCCC", "b": "AGAAAA"}
    offset = {"a": 1, "b": 2}
    seqs = mk_cls(data, moltype="dna", offset=offset)
    got = seqs.get_seq(seq)
    assert got._seq.offset == offset[seq]

    # if data is a SeqsData object, this should fail
    data = data_cls.from_seqs(data=data, alphabet=new_moltype.DNA.degen_gapped_alphabet)
    with pytest.raises(ValueError):
        _ = mk_cls(data, moltype="dna", offset=offset)

    # if provided with sequence objects with offsets, they should be propogated
    seq_view_1 = new_sequence.SeqView(
        parent="AAAA", parent_len=4, alphabet=dna_alphabet, offset=1
    )
    seq_view_2 = new_sequence.SeqView(
        parent="TTTT", parent_len=4, alphabet=dna_alphabet, offset=2
    )
    seq_1 = new_moltype.DNA.make_seq(seq=seq_view_1, name="seq_1")
    seq_2 = new_moltype.DNA.make_seq(seq=seq_view_2, name="seq_2")

    seqs = mk_cls([seq_1, seq_2], moltype="dna")

    got = seqs.get_seq("seq_1")
    assert got._seq.offset == 1

    got = seqs.get_seq("seq_2")
    assert got._seq.offset == 2


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_make_seqs_invalid_chars(mk_cls):
    data = {"seq1": "AGT1CCT", "seq2": "AGT$CCC"}
    with pytest.raises(new_alphabet.AlphabetError):
        mk_cls(data, moltype="dna")


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_sequence_collection_names_is_list(mk_cls):
    """expected to be a list"""
    seqs = mk_cls({"a": b"AAAAA", "b": b"TTTTT"}, moltype="dna")
    assert isinstance(seqs.names, list)
    assert seqs.names == ["a", "b"]


def test_sequence_collection_init_ordered(ordered1, ordered2):
    """SequenceCollection should iterate over seqs correctly even if ordered"""
    first = ordered1
    sec = ordered2

    assert first.names == ["a", "c"]
    assert sec.names == ["c", "a"]


@pytest.mark.parametrize("load_cls", [load_aligned_seqs, load_unaligned_seqs])
def test_sequence_collection_info_source(load_cls):
    """info.source exists if load seqs given a filename"""
    path = pathlib.Path("data/brca1.fasta")
    seqs = load_cls(path, moltype="dna", new_type=True)
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


def test_sequence_collection_iter_seqs_renamed(ragged_padded_dict):
    """SequenceCollection iter_seqs() method should support renaming of seqs"""
    coll = new_alignment.make_unaligned_seqs(ragged_padded_dict, moltype="dna")
    renamed = coll.rename_seqs(renamer=lambda x: f"seq_{x}")
    seqs = list(renamed.iter_seqs())
    assert seqs == ["AAAAAA", "AAA---", "AAAA--"]


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
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
@pytest.mark.parametrize(
    "sample", [["a"], ["b"], ["c"], ["a", "b"], ["a", "c"], ["b", "c"]]
)
def test_sequence_collection_take_seqs(ragged_padded_dict, mk_cls, sample):
    """take_seqs should return new SequenceCollection/Alignment with selected seqs."""
    orig = mk_cls(ragged_padded_dict, moltype="dna")
    subset = orig.take_seqs(sample)
    assert subset.names == sample
    assert subset.num_seqs == len(sample)
    # should be able to negate
    neg_subset = orig.take_seqs(sample, negate=True)
    assert neg_subset.names == [name for name in orig.names if name not in sample]


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_take_seqs_rc(mk_cls):
    data = {"a": "ACGT", "b": "CGTA", "c": "TTTT"}
    orig = mk_cls(data, moltype="dna")
    rc_seqs = orig.rc().take_seqs(["a", "b"])
    got = rc_seqs.to_dict()
    expect = {"a": "ACGT", "b": "TACG"}
    assert got == expect


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_take_seqs_str(ragged_padded_dict, mk_cls):
    """string arg to SequenceCollection take_seqs should work."""
    orig = mk_cls(ragged_padded_dict, moltype="dna")
    subset = orig.take_seqs("a")
    assert subset.names == ["a"]
    assert subset.num_seqs == 1

    # should be able to negate
    neg_subset = orig.take_seqs("a", negate=True)
    assert neg_subset.names == ["b", "c"]
    assert neg_subset.num_seqs == 2


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_take_seqs_info(ragged_padded_dict, mk_cls):
    """take_seqs should preserve info attribute"""
    orig = mk_cls(
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
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_take_seqs_copy_annotations(gff_db, mk_cls):
    data = {"test_seq": "ACGT--", "test_seq2": "CGTTTA"}
    seq_coll = mk_cls(data, moltype="dna", annotation_db=gff_db)
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
    assert ragged_padded.get_seq_names_if(is_med, negate=True) == ["b"]
    assert ragged_padded.get_seq_names_if(is_any) == ["a", "b", "c"]
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


@pytest.fixture
def gap_ambig_seqs():
    return {"s1": "ATGRY?", "s2": "T-AG??"}


@pytest.mark.parametrize("rced", [False, True])
@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_degap(mk_cls, gap_ambig_seqs, rced):
    """SequenceCollection.degap should strip gaps from each seq"""

    seqs = mk_cls(gap_ambig_seqs, moltype="dna")
    # Test normal case
    degapped = seqs.rc().degap() if rced else seqs.degap()
    got = degapped.to_dict()

    expect = {"s1": "RYCAT", "s2": "CTA"} if rced else {"s1": "ATGRY", "s2": "TAG"}
    assert got == expect

    # Test empty sequences case
    empty_seqs = mk_cls({"empty1": "", "empty2": ""}, moltype="dna")
    empty_seqs = empty_seqs.degap()
    got_empty = empty_seqs.to_dict()
    expect_empty = {"empty1": "", "empty2": ""}
    assert got_empty == expect_empty


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_degap_info(mk_cls, gap_ambig_seqs):
    """.degap should preserve info attributes"""
    aln = mk_cls(gap_ambig_seqs, moltype="dna")
    aln.info.path = "blah"
    got = aln.degap()
    assert got.info.path == "blah"


def test_alignment_degap_sliced(gap_ambig_seqs):
    """degap should apply slice_record to alignment"""
    aln = new_alignment.make_aligned_seqs(gap_ambig_seqs, moltype="dna")
    sliced = aln[:3]
    got = sliced.degap()
    expect = {"s1": "ATG", "s2": "TA"}
    assert got.to_dict() == expect

    # ATGRY?
    # T-AG??
    #  * * *
    sliced = aln[1::2]
    got = sliced.degap()
    expect = {"s1": "TR", "s2": "G"}
    assert got.to_dict() == expect

    sliced = aln[::-1]
    got = sliced.degap()
    expect = {"s1": "YRGTA", "s2": "GAT"}


def test_get_degapped_relative_to():
    """should remove all columns with a gap in sequence with given name"""
    aln = new_alignment.make_aligned_seqs(
        [
            ["name1", "-AC-DEFGHI---"],
            ["name2", "XXXXXXXXXXXXX"],
            ["name3", "YYYY-YYYYYYYY"],
            ["name4", "-KL---MNPR---"],
        ],
        moltype="protein",
    )
    expect = dict(
        [
            ["name1", "ACDEFGHI"],
            ["name2", "XXXXXXXX"],
            ["name3", "YY-YYYYY"],
            ["name4", "KL--MNPR"],
        ]
    )
    result = aln.get_degapped_relative_to("name1")
    assert result.to_dict() == expect

    with pytest.raises(ValueError):
        aln.get_degapped_relative_to("nameX")


def test_get_degapped_relative_to_no_or_all_gaps():
    """should handle case where no gaps are present in reference
    or when the reference is all gaps"""
    aln = new_alignment.make_aligned_seqs(
        [
            ["name1", "-AC-DEFGHI---"],
            ["name2", "XXXXXXXXXXXXX"],
            ["name3", "YYYY-YYYYYYYY"],
            ["name4", "-------------"],
        ],
        moltype="protein",
    )
    expect = dict(
        [
            ["name1", "-AC-DEFGHI---"],
            ["name2", "XXXXXXXXXXXXX"],
            ["name3", "YYYY-YYYYYYYY"],
            ["name4", "-------------"],
        ]
    )
    got = aln.get_degapped_relative_to("name2").to_dict()
    assert got == expect

    got = aln.get_degapped_relative_to("name4").to_dict()
    assert got == {"name1": "", "name2": "", "name3": "", "name4": ""}


@pytest.mark.parametrize("rc", [False, True])
def test_get_degapped_relative_to_sliced(rc):
    aln = new_alignment.make_aligned_seqs(
        [
            ["name1", "-AC-TTT"],
            ["name2", "AAAAAA-"],
            ["name3", "YYYY-YY"],
            ["name4", "-AC---G"],
        ],
        moltype="dna",
    )
    sliced = aln[:5].rc() if rc else aln[:5]
    expect = dict(
        [
            ["name1", "ACT"],
            ["name2", "AAA"],
            ["name3", "YY-"],
            ["name4", "AC-"],
        ]
    )
    expect = (
        {name: aln.moltype.complement(seq[::-1]) for name, seq in expect.items()}
        if rc
        else expect
    )
    result = sliced.get_degapped_relative_to("name1")
    assert result.to_dict() == expect


def test_get_degapped_relative_to_info():
    """should remove all columns with a gap in sequence with given name
    while preserving info attribute"""
    aln = new_alignment.make_aligned_seqs(
        [
            ["name1", "-AC-DEFGHI---"],
            ["name2", "XXXXXX--XXXXX"],
            ["name3", "YYYY-YYYYYYYY"],
            ["name4", "-KL---MNPR---"],
        ],
        moltype="protein",
        info={"key": "foo"},
    )
    out_aln = new_alignment.make_aligned_seqs(
        [
            ["name1", "ACDEFGHI"],
            ["name2", "XXXX--XX"],
            ["name3", "YY-YYYYY"],
            ["name4", "KL--MNPR"],
        ],
        moltype="protein",
        info={"key": "bar"},
    )
    gdrt = aln.get_degapped_relative_to("name1")
    assert gdrt.info["key"] == "foo"


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_to_fasta(mk_cls):
    """SequenceCollection and Alignment should return correct FASTA string"""
    data = {"seq_0": "AAA", "seq_1": "CCC"}
    seqs = mk_cls(data, moltype="dna")
    assert seqs.to_fasta() == ">seq_0\nAAA\n>seq_1\nCCC\n"
    assert seqs.to_fasta(block_size=2) == ">seq_0\nAA\nA\n>seq_1\nCC\nC\n"

    data = {"seq_0": "GCATGCAT", "seq_1": "TCAGACGT"}
    seqs = mk_cls(data, moltype="dna")
    assert seqs.to_fasta() == ">seq_0\nGCATGCAT\n>seq_1\nTCAGACGT\n"
    assert seqs.to_fasta(block_size=4) == ">seq_0\nGCAT\nGCAT\n>seq_1\nTCAG\nACGT\n"
    assert seqs.to_fasta(block_size=3) == ">seq_0\nGCA\nTGC\nAT\n>seq_1\nTCA\nGAC\nGT\n"


def test_sequence_collection_is_ragged(ragged, ragged_padded):
    """SequenceCollection is_ragged should return true if ragged alignment"""
    assert ragged.is_ragged()
    assert not ragged_padded.is_ragged()


def test_sequence_collection_ragged(ragged):
    assert [ragged.seqs[name] for name in ragged.names] == [
        "AAAAAA",
        "AAA",
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
    "mk_cls", [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs]
)
@pytest.mark.parametrize(
    "gc,seqs",
    (
        (1, ("TCCTGA", "GATTT?")),
        (1, ("ACGTAA---", "ACGAC----", "ACGCAATGA")),
        (2, ("GATTTT", "TCCAGG")),
    ),
)
def test_sequence_collection_has_terminal_stop_true(gc, seqs, mk_cls):
    data = {f"s{i}": s for i, s in enumerate(seqs)}
    seq_coll = mk_cls(data, moltype="dna")
    assert seq_coll.has_terminal_stop(gc=gc)


@pytest.mark.parametrize(
    "mk_cls", [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs]
)
@pytest.mark.parametrize(
    "gc,seqs",
    ((1, ("TCCTCA", "GATTTT")), (2, ("GATTTT", "TCCCGG")), (1, ("CCTCA", "ATTTT"))),
)
def test_sequence_collection_has_terminal_stop_false(gc, seqs, mk_cls):
    data = {f"s{i}": s for i, s in enumerate(seqs)}
    seq_coll = mk_cls(data, moltype="dna")
    assert not seq_coll.has_terminal_stop(gc=gc)


@pytest.mark.parametrize(
    "mk_cls", [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs]
)
def test_sequence_collection_has_terminal_stop_strict(mk_cls):
    data = {f"s{i}": s for i, s in enumerate(("CCTCA", "ATTTT"))}
    seq_coll = mk_cls(data, moltype="dna")
    with pytest.raises(new_alphabet.AlphabetError):
        seq_coll.has_terminal_stop(gc=1, strict=True)


@pytest.mark.parametrize(
    "mk_cls,expect",
    (
        (new_alignment.make_unaligned_seqs, {"seq1": "DS", "seq2": "DSS"}),
        (new_alignment.make_aligned_seqs, {"seq1": "DS-", "seq2": "DSS"}),
    ),
)
def test_get_translation_trim_stop(mk_cls, expect):
    data = {"seq1": "GATTCCTAG", "seq2": "GATTCCTCC"}
    seqs = mk_cls(data, moltype="dna")
    got = seqs.get_translation(trim_stop=True)
    assert got.to_dict() == expect


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_get_translation_raises(mk_cls):
    """should raise error if self.moltype is not a nucleic acid"""
    data = {"seq1": "PAR", "seq2": "PQR"}
    seqs = mk_cls(data, moltype="protein")
    with pytest.raises(new_moltype.MolTypeError):
        _ = seqs.get_translation(trim_stop=True)


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
@pytest.mark.parametrize(
    "seqs",
    (
        {"seq1": "GATTTT", "seq2": "GATC??"},
        {"seq1": "GAT---", "seq2": "GATCTT"},
        {"seq1": "GAT-T-", "seq2": "GATCTT"},
        {"seq1": "GATTTT", "seq2": "?GATCT"},
    ),
)
def test_get_translation(seqs, mk_cls):
    """SequenceCollection.get_translation translates each seq"""
    seqs = mk_cls(seqs, moltype="dna")
    got = seqs.get_translation(incomplete_ok=True)
    assert got.num_seqs == 2
    assert got.moltype == new_moltype.PROTEIN


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_get_translation_renamed(seqs, mk_cls):
    seqs = mk_cls({"s1": "GATTTT", "s2": "GATCTT"}, moltype="dna")
    renamed = seqs.rename_seqs(lambda x: f"{x}_renamed")
    got = renamed.get_translation(incomplete_ok=True)
    assert got.names == ["s1_renamed", "s2_renamed"]


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_get_translation_with_stop(mk_cls):
    data = {"seq1": "?GATAG", "seq2": "GATTAG"}
    seqs = mk_cls(data, moltype="dna")
    got = seqs.get_translation(incomplete_ok=True, include_stop=True, trim_stop=False)
    assert got.to_dict() == {"seq1": "X*", "seq2": "D*"}
    assert got.moltype == new_moltype.PROTEIN_WITH_STOP


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_get_translation_non_div_3(mk_cls):
    data = {"seq1": "?GATCTA", "seq2": "GATTAGG"}
    seqs = mk_cls(data, moltype="dna")
    got = seqs.get_translation(incomplete_ok=True, include_stop=True)
    assert got.to_dict() == {"seq1": "XS", "seq2": "D*"}


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
@pytest.mark.parametrize(
    "data", ({"seq1": "GATTTT", "seq2": "GATC??"}, {"seq1": "GAT---", "seq2": "?GATCT"})
)
def test_get_translation_error(data, mk_cls):
    seqs = mk_cls(data, moltype="dna")
    with pytest.raises(TypeError):
        seqs.get_translation()


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
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
def test_get_translation_info(data, mk_cls):
    """SequenceCollection.get_translation preserves info attribute"""
    seqs = mk_cls(data, moltype="dna", info={"key": "value"})
    got = seqs.get_translation(incomplete_ok=True)
    assert got.info["key"] == "value"


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_get_translation_incomplete(mk_cls):
    """get translation works on incomplete codons"""
    data = {"seq1": "GATN--", "seq2": "?GATCT"}
    seqs = mk_cls(data, moltype="dna")
    got = seqs.get_translation(incomplete_ok=True)
    assert got.to_dict() == {"seq1": "DX", "seq2": "XS"}
    with pytest.raises(new_alphabet.AlphabetError):
        _ = seqs.get_translation(incomplete_ok=False)


@pytest.mark.parametrize(
    "gc,seqs",
    (
        (1, ("--AT-CTGA", "GATAAATT?")),
        (1, ("ACGTGA---", "ACGAC----", "ACGCAATGA")),
        (1, ("CCTCA-", "ATTTTA")),
        (2, ("GATTTT", "TCCAGG")),
    ),
)
@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_trim_stop_codons(gc, seqs, mk_cls):
    data = {f"s{i}": s for i, s in enumerate(seqs)}

    expect = {}
    for k, v in data.items():
        if "-" in v or mk_cls == new_alignment.make_aligned_seqs:
            v = re.sub("(TGA|AGG)(?=[-]*$)", "---", v)
        else:
            v = re.sub("(TGA|AGG)", "", v)
        expect[k] = v
    seqs = mk_cls(data, moltype="dna")
    got = seqs.trim_stop_codons(gc=gc).to_dict()

    assert got == expect


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
@pytest.mark.parametrize(
    "gc,seqs",
    ((1, ("T-CTGC", "GATAA?")), (2, ("GATTTT", "TCCCGG")), (1, ("CCTGC", "GATAA"))),
)
def test_trim_stop_codons_no_stop(gc, seqs, mk_cls):
    data = {f"s{i}": s for i, s in enumerate(seqs)}
    seqs = mk_cls(data, moltype="dna")
    got = seqs.trim_stop_codons(gc=gc)
    assert got is seqs


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
@pytest.mark.parametrize(
    "data", ({"s1": "CCTCA", "s2": "ATTTT"}, {"s1": "CCTCA-", "s2": "ATTTTA"})
)
def test_trim_stop_codons_strict(data, mk_cls):
    seqs = mk_cls(data, moltype="dna")
    with pytest.raises(new_alphabet.AlphabetError):
        seqs.trim_stop_codons(gc=1, strict=True)


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_trim_stop_codons_info(mk_cls):
    """trim_stop_codons should preserve info attribute"""
    data = {"seq1": "ACGTAA", "seq2": "ACGACG", "seq3": "ACGCGT"}
    seqs = mk_cls(
        data,
        moltype="dna",
        info={"key": "value"},
    )
    seqs = seqs.trim_stop_codons()
    assert seqs.info["key"] == "value"


@pytest.mark.parametrize(
    "mk_cls",
    (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
)
def test_trim_stop_codons_annotation_db(gff_db, mk_cls):
    """trim_stop_codons should preserve info attribute"""
    data = {"seq_1": "ACGTAA", "seq_2": "ACGACG", "seq_3": "ACGCGT"}
    seqs = mk_cls(
        data,
        moltype="dna",
        info={"key": "value"},
    )
    seqs.annotation_db = gff_db
    trimmed = seqs.trim_stop_codons()
    assert trimmed.annotation_db == seqs.annotation_db


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


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_sequence_collection_counts_per_seq(mk_cls):
    """SequenceCollection.counts_per_seq handles motif length, allow_gaps etc.."""
    data = {"a": "AAAA??????", "b": "CCCGGG--NN", "c": "CCGGTTCCAA"}
    coll = mk_cls(data, moltype="dna")
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


def test_counts_per_seq_text_moltype():
    """produce correct counts per seq with text moltypes"""
    data = {"a": "AAAA??????", "b": "CCCGGG--NN", "c": "CCGGTTCCAA"}
    coll = new_alignment.make_aligned_seqs(data, moltype="text")
    got = coll.counts_per_seq(include_ambiguity=True, allow_gap=True)
    assert got.col_sum()["-"] == 2
    assert got.col_sum()["?"] == 6
    assert got.col_sum()["T"] == 2


def test_counts_per_pos_text_moltype():
    """produce correct counts per pos with default moltypes"""
    data = {"a": "AAAA??????", "b": "CCCGGG--NN", "c": "CCGGTTCCAA"}
    coll = new_alignment.make_aligned_seqs(data, moltype="text")
    got = coll.counts_per_pos()
    # should not include gap character
    assert "-" not in got.motifs
    # allowing gaps
    got = coll.counts_per_pos(allow_gap=True)
    # should include gap character
    assert got[5, "-"] == 0
    assert got[6, "-"] == 1

    # now with motif-length 2
    got = coll.counts_per_pos(motif_length=2)
    found_motifs = set()
    lengths = set()
    for m in got.motifs:
        lengths.add(len(m))
        found_motifs.update(m)
    assert "-" not in found_motifs
    assert lengths == {2}


@pytest.mark.parametrize(
    "mk_cls", (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs)
)
def test_sequence_collection_probs_per_seq(mk_cls):
    data = {"seq1": "AA??", "seq2": "CG-N", "seq3": "CGAA"}
    coll = mk_cls(data, moltype="dna")
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


@pytest.mark.parametrize(
    "mk_cls", [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs]
)
def test_sequence_collection_count_ambiguous_per_seq(mk_cls):
    data = {
        "a": "AGGGT",
        "b": "AGGN?",
        "c": "?????",
        "d": "-----",
    }
    coll = mk_cls(data, moltype="dna")
    got = coll.count_ambiguous_per_seq()
    expect = {"a": 0, "b": 2, "c": 5, "d": 0}
    assert got == expect


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


@pytest.mark.parametrize("rc", (True, False))
def test_sequence_collection_init_seqs_rc(rc):
    """the rc status of the input seqs is propagated and annotations maintained
    when all seqs have the same rc status"""
    data = {"seq1": _make_seq("seq1"), "seq2": _make_seq("seq2")}
    data = {k: v.rc() for k, v in data.items()} if rc else data
    seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    assert seq_coll._is_reversed == rc
    assert len(seq_coll.annotation_db) > 0
    got = seq_coll.to_dict()
    expect = {k: str(v) for k, v in data.items()}
    assert got == expect

    # when is_reversed is specified, it should override the rc status of the input seqs
    seq_coll = new_alignment.make_unaligned_seqs(
        data, moltype="dna", is_reversed=not rc
    )
    assert seq_coll._is_reversed == (not rc)
    assert len(seq_coll.annotation_db) > 0
    got = seq_coll.to_dict()
    expect = {k: str(v.rc()) for k, v in data.items()}
    assert got == expect


def test_sequence_collection_init_seqs_mixed_rc():
    """annotations on input seqs with mixed rc means dropping annotations, the
    rc of the collection is set to False"""
    data = {"seq1": _make_seq("seq1"), "seq2": _make_seq("seq2").rc()}
    with catch_warnings():
        filterwarnings("ignore", category=UserWarning)
        seq_coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
    coll_db = seq_coll.annotation_db
    assert not len(coll_db)
    assert not seq_coll._is_reversed


def test_sequence_collection_copy_annotations_incompat_type_fails(seqcoll_db, seqs):
    with pytest.raises(TypeError):
        seqcoll_db.copy_annotations(seqs)


@pytest.mark.parametrize(
    "mk_cls, seq_cls",
    [
        (new_alignment.make_unaligned_seqs, new_sequence.Sequence),
        (new_alignment.make_aligned_seqs, new_alignment.Aligned),
    ],
)
@pytest.mark.parametrize("moltype", ("dna", "rna", "protein"))
def test_sequence_collection_to_moltype(moltype, mk_cls, seq_cls):
    data = dict(s1="ACGAA-", s2="ACCCAA")
    seqs = mk_cls(data, moltype="text")
    mt_seqs = seqs.to_moltype(moltype=moltype)
    got = mt_seqs.seqs["s1"]
    assert got.moltype.label == moltype
    assert isinstance(got, seq_cls)

    # should also work with moltype objects
    mt = new_moltype.get_moltype(moltype)
    mt_seqs = seqs.to_moltype(moltype=mt)
    got = mt_seqs.seqs["s1"]
    assert got.moltype.label == moltype
    assert isinstance(got, seq_cls)


@pytest.mark.parametrize("rc", (True, False))
@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_to_moltype_rc(rc, mk_cls):
    data = dict(s1="TCAG", s2="TCA-")
    seqs = mk_cls(data, moltype="dna")
    seqs = seqs.rc() if rc else seqs
    got = seqs.to_moltype("rna").to_dict()
    expect = {"s1": "CUGA", "s2": "-UGA"} if rc else dict(s1="UCAG", s2="UCA-")
    assert got == expect


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
@pytest.mark.parametrize("moltype", ("dna", "protein"))
def test_sequence_collection_to_moltype_same_moltype(moltype, mk_cls):
    data = dict(s1="ACGTT", s2="ACCTT")
    seqs = mk_cls(data, moltype=moltype)
    got = seqs.to_moltype(moltype=moltype)
    assert got is seqs


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_to_moltype_with_gaps(mk_cls):
    """correctly convert to specified moltype"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = mk_cls(data, moltype="text")
    dna_seqs = seqs.to_moltype("dna")
    assert dna_seqs.moltype.label == "dna"
    assert dna_seqs.to_dict() == data

    # we can convert from text to dna to rna
    rna_seqs = dna_seqs.to_moltype("rna")
    assert rna_seqs.moltype.label == "rna"

    # but not from text to rna directly
    with pytest.raises(new_moltype.MolTypeError):
        seqs.to_moltype("rna")


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_to_moltype_info(mk_cls):
    """correctly convert to specified moltype"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = mk_cls(data, moltype="text", info={"key": "value"})
    dna = seqs.to_moltype("dna")
    assert dna.info["key"] == "value"


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_to_moltype_annotation_db(mk_cls):
    """correctly convert to specified moltype"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = mk_cls(data, moltype="text")
    db = GffAnnotationDb()
    db.add_feature(seqid="seq1", biotype="exon", name="annotation1", spans=[(3, 8)])
    db.add_feature(seqid="seq2", biotype="exon", name="annotation2", spans=[(1, 2)])
    db.add_feature(seqid="seq3", biotype="exon", name="annotation3", spans=[(3, 6)])
    seqs.annotation_db = db
    dna = seqs.to_moltype("dna")
    assert len(dna.annotation_db) == 3


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_to_dna(mk_cls):
    """correctly convert to dna"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = mk_cls(data, moltype="text")

    # convert from text to dna
    dna_seqs = seqs.to_dna()
    assert dna_seqs.moltype.label == "dna"
    assert dna_seqs.to_dict() == data

    # convert from rna to dna
    data = {"seq1": "ACGUACGUA", "seq2": "ACCGAA---", "seq3": "ACGUACGUU"}
    seqs = mk_cls(data, moltype="RNA")
    dna_seqs = seqs.to_dna()
    assert dna_seqs.moltype.label == "dna"
    assert dna_seqs.to_dict() == {
        name: seq.replace("U", "T") for name, seq in data.items()
    }


def test_to_dna():
    """alignment cast to DNA works"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    aln = new_alignment.make_aligned_seqs(data, moltype="text")
    dna = aln.to_dna()
    assert set(dna.names) == set(aln.names)
    assert dna.moltype.label == "dna"
    # should fail if invalid character set
    paln = dna.get_translation()
    with pytest.raises(new_moltype.MolTypeError):
        _ = paln.to_dna()


def test_to_dna_info():
    """alignment cast to DNA preserves info attribute"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    aln = new_alignment.make_aligned_seqs(data, info={"key": "value"}, moltype="text")
    dna = aln.to_dna()
    assert dna.info["key"] == "value"


def test_to_rna():
    """alignment cast to RNA works"""
    data = {"seq1": "ACGUACGUA", "seq2": "ACCGAA---", "seq3": "ACGUACGUU"}
    aln = new_alignment.make_aligned_seqs(data, moltype="text")
    rna = aln.to_rna()
    assert set(rna.names) == set(aln.names)
    assert rna.moltype.label == "rna"


def test_to_rna_info():
    """alignment cast to RNA preserves info attribute"""
    data = {"seq1": "ACGUACGUA", "seq2": "ACCGAA---", "seq3": "ACGUACGUU"}
    aln = new_alignment.make_aligned_seqs(data, info={"key": "value"}, moltype="text")
    rna = aln.to_rna()
    assert rna.info["key"] == "value"


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_to_rna(mk_cls):
    data = {"seq1": "ACGUACGUA", "seq2": "ACCGAA---", "seq3": "ACGUACGUU"}
    seqs = mk_cls(data, moltype="text")

    # convert from text to rna
    rna_seqs = seqs.to_rna()
    assert rna_seqs.moltype.label == "rna"
    assert rna_seqs.to_dict() == data

    # convert from dna to rna
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = mk_cls(data, moltype="dna")
    rna_seqs = seqs.to_rna()
    assert rna_seqs.moltype.label == "rna"
    assert rna_seqs.to_dict() == {
        name: seq.replace("T", "U") for name, seq in data.items()
    }


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_add_seqs(mk_cls):
    data = dict(
        [("name1", "AAA"), ("name2", "A--"), ("name3", "AAA"), ("name4", "AAA")]
    )
    data2 = dict([("name5", "TTT"), ("name6", "---")])
    aln = mk_cls(data, moltype="dna", info={"key": "foo"})
    out_aln = aln.add_seqs(data2)
    assert len(out_aln.names) == 6


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_add_seqs_reversed(mk_cls):
    data = dict(
        [("name1", "AAA"), ("name2", "A--"), ("name3", "AAA"), ("name4", "AAA")]
    )
    data2 = dict([("name5", "TTT"), ("name6", "---")])
    aln = mk_cls(data, moltype="dna", info={"key": "foo"})
    out_aln = aln.add_seqs(data2)
    assert len(out_aln.names) == 6


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_add_seqs_duplicate_raises(mk_cls):
    """add_seqs should raise an error if duplicate names"""
    data = dict(
        [("name1", "AAA"), ("name2", "AAA"), ("name3", "AAA"), ("name4", "AAA")]
    )
    data2 = dict([("name5", "BBB"), ("name1", "CCC")])
    aln = mk_cls(data, moltype="text", info={"key": "foo"})
    with pytest.raises(ValueError):
        _ = aln.add_seqs(data2)

    # but it should work if we allow duplicates
    out_aln = aln.add_seqs(data2, force_unique_keys=False)
    assert len(out_aln.names) == 5


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_add_seqs_info(mk_cls):
    """add_seqs should preserve info attribute"""
    data = dict(
        [("name1", "AAA"), ("name2", "AAA"), ("name3", "AAA"), ("name4", "AAA")]
    )
    data2 = dict([("name5", "BBB"), ("name6", "CCC")])
    aln = mk_cls(data, moltype="text", info={"key": "foo"})
    aln2 = mk_cls(data2, moltype="text", info={"key": "bar"})
    out_aln = aln.add_seqs(aln2)
    assert out_aln.info["key"] == "foo"


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
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_write_gapped(mk_cls, tmp_path):
    data = {"a": "AAA--", "b": "TTTTT", "c": "CCCCC"}
    seqs = mk_cls(data, moltype="dna")
    fn = tmp_path / "seqs.fasta"
    seqs.write(fn)
    with open(fn, newline=None) as infile:
        result = infile.read()
    assert result == ">a\nAAA--\n>b\nTTTTT\n>c\nCCCCC\n"


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_get_ambiguous_positions(mk_cls):
    aln = mk_cls({"s1": "ATGRY?", "s2": "T-AG??"}, moltype="dna")
    assert aln.get_ambiguous_positions() == {
        "s2": {4: "?", 5: "?"},
        "s1": {3: "R", 4: "Y", 5: "?"},
    }


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_consistent_gap_degen_handling(mk_cls):
    """gap degen character should be treated consistently"""
    # the degen character '?' can be a gap, so when we strip gaps it should
    # be gone too
    raw_seq = "---??-??TC-GGCG-GCA-G-GC-?-C-TAN-GCGC-CCTC-AGGA?-???-??--"
    raw_ungapped = re.sub("[-?]", "", raw_seq)
    re.sub("[N?]+", "", raw_seq)
    dna = new_moltype.DNA.make_seq(seq=raw_seq)

    aln = mk_cls({"a": dna, "b": dna}, moltype="dna")
    expect = mk_cls({"a": raw_ungapped, "b": raw_ungapped}, moltype="dna").to_fasta()
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

    # will pad to max sequence length if pad_length less than max length
    padded3 = ragged.pad_seqs(pad_length=5)
    seqs3 = list(padded3.iter_seqs(seq_order=["a", "b", "c"]))
    assert list(map(str, seqs3)) == ["AAAAAA", "AAA---", "AAAA--"]


def test_sequence_collection_pad_seqs_reversed():
    mk_seq = new_moltype.DNA.make_seq
    data = {
        "a": mk_seq(seq="T", name="A"),
        "b": mk_seq(seq="TCG", name="B"),
        "c": mk_seq(seq="TTGG", name="C"),
    }
    ragged = new_alignment.make_unaligned_seqs(data, moltype="dna")
    rc = ragged.rc()
    padded = rc.pad_seqs()
    got = padded.to_dict()
    expect = {"a": "A---", "b": "CGA-", "c": "CCAA"}
    assert got == expect


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


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_rename_seqs(mk_cls):
    """successfully rename sequences"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = mk_cls(data, moltype="dna")
    new = seqs.rename_seqs(lambda x: x.upper())
    expect = {n.upper() for n in data}
    assert set(new.names) == expect
    # the names should not change in the seqsdata
    assert set(new._seqs_data.names) == set(data)


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_subsequent_rename(mk_cls):
    """sequences can be renamed multiple times"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = mk_cls(data, moltype="dna")
    new = seqs.rename_seqs(lambda x: x.upper())
    new_again = new.rename_seqs(lambda x: f"{x[0]}{x[-1]}")
    expect = {"S1", "S2", "S3"}
    assert set(new_again.names) == expect
    # the names should not change in the seqsdata
    assert set(new._seqs_data.names) == set(data)


@pytest.mark.parametrize(
    "mk_cls",
    [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs],
)
def test_sequence_collection_rename_non_unique_fails(mk_cls):
    """renaming to non-unique names should raise an error"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = mk_cls(data, moltype="dna")
    with pytest.raises(ValueError):
        _ = seqs.rename_seqs(lambda x: x[:1])


def test_alignment_rename_sliced():
    """a sliced alignment that is renmed will retain slice information"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    sliced = aln[:3]
    renamed = sliced.rename_seqs(lambda x: x.upper())
    expect = {n.upper(): seq[:3] for n, seq in data.items()}
    assert renamed.to_dict() == expect


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
    expect = {
        "seqs": data,
        "type": get_object_provenance(seqs),
        "version": __version__,
        "init_args": {
            "moltype": seqs.moltype.label,
            "name_map": seqs._name_map,
            "info": seqs.info,
        },
    }
    assert got == expect


def test_sequence_collection_to_rich_dict_reversed_seqs():
    data = {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    reversed_seqs = seqs.rc()

    got = reversed_seqs.to_rich_dict()
    expect = {
        "seqs": reversed_seqs.to_dict(),
        "init_args": {
            "moltype": seqs.moltype.label,
            "name_map": seqs._name_map,
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


@pytest.mark.parametrize("rc", (False, True))
def test_sequence_collection_round_trip(rc):
    seq_coll = new_alignment.make_unaligned_seqs(
        {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}, moltype="dna"
    )
    seq_coll = seq_coll.rc() if rc else seq_coll

    rd = seq_coll.to_rich_dict()
    got = deserialise_object(rd)
    assert isinstance(got, new_alignment.SequenceCollection)
    assert got.to_rich_dict() == seq_coll.to_rich_dict()


def test_sequence_collection_distance_matrix_singleton_collection(dna_moltype):
    """SequenceCollection.distance_matrix() should raise error if collection
    only contains a single sequence"""
    collection = new_alignment.make_unaligned_seqs(
        {"s1": "ACGTACGTAGTCGCG"}, moltype=dna_moltype
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
    data = {"s1": "ACGA", "s2": "ACGA"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype=moltype)
    got = seqs.distance_matrix()
    assert got[("s1", "s2")] == 0.0


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
def aligned_seqs_data(aligned_dict, dna_alphabet):
    return new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_dict, alphabet=dna_alphabet
    )


@pytest.fixture
def gap_seqs():
    return [
        ("A---CTG-C", [[1, 3], [4, 4]]),
        ("-GTAC----", [[0, 1], [4, 5]]),
        ("---A--T--", [[0, 3], [1, 5], [2, 7]]),
    ]


@pytest.mark.parametrize("i", range(3))
def test_decompose_gapped_seq_sequences(gap_seqs, i, dna_alphabet):
    seq, gap_coords = gap_seqs[i]
    dna = new_moltype.get_moltype("dna")
    got_ungapped, got_map = new_alignment.decompose_gapped_seq(
        dna.make_seq(seq=seq), alphabet=dna_alphabet
    )
    expect = dna_alphabet.to_indices(seq.replace("-", ""))
    assert numpy.array_equal(got_ungapped, expect)
    assert numpy.array_equal(got_map, gap_coords)


@pytest.mark.parametrize("i", range(3))
def test_decompose_gapped_seq_str(gap_seqs, i, dna_alphabet):
    seq, gap_coords = gap_seqs[i]
    got_ungapped, got_map = new_alignment.decompose_gapped_seq(
        seq, alphabet=dna_alphabet
    )
    expect = dna_alphabet.to_indices(seq.replace("-", ""))
    assert numpy.array_equal(got_ungapped, expect)
    assert numpy.array_equal(got_map, gap_coords)


def test_decompose_gapped_seq_str_all_gaps(dna_alphabet):
    parent_seq = "-----"
    expect_gaplen = numpy.array([len(parent_seq)])
    got_ungap, got_map = new_alignment.decompose_gapped_seq(
        parent_seq, alphabet=dna_alphabet
    )
    expect = numpy.array([])
    assert numpy.array_equal(got_ungap, expect)
    assert got_map[:, 1] == expect_gaplen


def test_decompose_gapped_seq_str_no_gaps(dna_alphabet):
    parent_seq = "ACTGC"
    got_ungap, got_map = new_alignment.decompose_gapped_seq(
        parent_seq, alphabet=dna_alphabet
    )
    expect = dna_alphabet.to_indices(parent_seq)
    assert numpy.array_equal(got_ungap, expect)
    assert got_map.size == 0


def test_decompose_gapped_seq_arr_all_gaps(dna_alphabet):
    parent_seq = dna_alphabet.to_indices("-----")
    got_ungap, got_map = new_alignment.decompose_gapped_seq(
        parent_seq, alphabet=dna_alphabet
    )
    assert got_ungap.size == 0
    assert numpy.array_equal(got_map, numpy.array([[0, 5]]))


def test_decompose_gapped_seq_arr_no_gaps(dna_alphabet):
    parent_seq = dna_alphabet.to_indices("ACTGC")
    got_ungap, got_empty_arr = new_alignment.decompose_gapped_seq(
        parent_seq, alphabet=dna_alphabet
    )
    assert numpy.array_equal(got_ungap, parent_seq)
    assert got_empty_arr.size == 0


@pytest.mark.parametrize("i", range(3))
def test_decompose_gapped_seq_arr(gap_seqs, i, dna_alphabet):
    seq, gap_coords = gap_seqs[i]
    seq = dna_alphabet.to_indices(seq)
    got_ungapped, got_map = new_alignment.decompose_gapped_seq(
        seq, alphabet=dna_alphabet
    )
    expect = seq[seq != 4]  # gap_char = 4
    assert numpy.array_equal(got_ungapped, expect)
    assert numpy.array_equal(got_map, gap_coords)


@pytest.mark.parametrize("i", range(3))
def test_decompose_gapped_seq_arr_dispatch_equal(gap_seqs, i, dna_alphabet):
    """decompose_gapped_seq should return the same gap coords when input is a string/array/bytes"""
    seq_str, _ = gap_seqs[i]
    seq_array = dna_alphabet.to_indices(seq_str)
    seq_bytes = dna_alphabet.array_to_bytes(seq_array)
    seq_from_str, gaps_from_str = new_alignment.decompose_gapped_seq(
        seq_str, alphabet=dna_alphabet
    )
    seq_from_array, gaps_from_arr = new_alignment.decompose_gapped_seq(
        seq_array, alphabet=dna_alphabet
    )
    seq_from_bytes, gaps_from_bytes = new_alignment.decompose_gapped_seq(
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


@pytest.mark.parametrize(
    "raw_seq,coords",
    (("ACGGTAAAG", ((2, 4), (5, 8))), ("CCC---CCC", ((0, 3), (6, 9)))),
)
def test_aligned_getitem_featuremap(raw_seq, coords):
    dna = new_moltype.get_moltype("dna")
    im, seq = dna.make_seq(seq=raw_seq).parse_out_gaps()
    gaps = numpy.array([im.gap_pos, im.cum_gap_lengths]).T
    asd = new_alignment.AlignedSeqsData.from_seqs_and_gaps(
        seqs={"seq1": numpy.array(seq)},
        gaps={"seq1": gaps},
        alphabet=dna.most_degen_alphabet(),
    )
    aln = new_alignment.make_aligned_seqs(asd, moltype=dna)
    ia = aln.seqs["seq1"]
    length = len(raw_seq)
    fmap = FeatureMap(spans=[Span(s, e) for s, e in coords], parent_length=length)
    expect = "".join(raw_seq[s:e] for s, e in fmap.get_coordinates())
    got = ia[fmap]
    assert str(got) == expect


@pytest.fixture
def aligned():
    data = {
        "seq1": "AAAGG--GGG-AACCCT",
    }
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    return aln.seqs["seq1"]


def test_aligned_getitem_featuremap_allgap(aligned):
    fmap = FeatureMap(spans=[LostSpan(4)], parent_length=0)
    sliced = aligned[fmap]
    assert not sliced


def test_aligned_getitem_featuremap_multi_spans(aligned):
    #                      1111111
    #            01234567890123456
    #             ***   **    ***
    # raw_seq = "AAAGG--GGG-AACCCT"
    #            01234  567 890123
    #                         1111

    fmap = FeatureMap.from_locations(
        locations=[(1, 4), (7, 9), (13, 16)], parent_length=17
    )
    sliced = aligned[fmap]
    assert str(sliced) == "AAGGGCCC"


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
def test_aligned_seqs_data_init(seqid, gap_seqs, dna_alphabet):
    """test initialisation of AlignedSeqsData object with dictionary of ungapped
    sequences and dictionary of gap coordinates and cumulated gap lengths"""
    seqs = {f"seq{i}": seq.replace("-", "") for i, (seq, _) in enumerate(gap_seqs)}
    gaps = {
        f"seq{i}": numpy.array(gaps, dtype=numpy.int32)
        for i, (_, gaps) in enumerate(gap_seqs)
    }
    ad = new_alignment.AlignedSeqsData.from_seqs_and_gaps(
        seqs={**seqs}, gaps=gaps, alphabet=dna_alphabet
    )
    got = ad.get_seq_str(seqid=seqid)
    assert got == seqs[seqid]
    assert numpy.array_equal(ad.get_gaps(seqid=seqid), gaps[seqid])


@pytest.mark.parametrize(
    "kwargs",
    (
        dict(offset=dict(s4=1)),
        dict(gaps=dict(s4=numpy.array([[0, 1]]))),
        dict(ungapped_seqs=dict(s4=numpy.array([1, 2, 3]))),
    ),
)
def test_aligned_seqs_data_init_check_raises(dna_alphabet, kwargs):
    raw = (
        "A---CTG-C",
        "-GTAC----",
        "---A--T--",
    )
    gapped = numpy.array([dna_alphabet.to_indices(seq) for seq in raw])
    names = "s1", "s2", "s3"
    with pytest.raises(ValueError):
        new_alignment.AlignedSeqsData(
            gapped_seqs=gapped, names=names, alphabet=dna_alphabet, **kwargs
        )


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4", "seq5"))
@pytest.mark.parametrize("data_type", (str, numpy.array, bytes))
def test_aligned_seqs_data_init_gapped(
    gapped_seqs_dict, seqid, data_type, dna_alphabet, dna_moltype
):
    """AlignedSeqsData should handle data from gapped sequences correctly,
    including edge cases (all gaps, no gaps, etc)"""

    typed_data = {
        name: make_typed(seq, data_type=data_type, moltype=dna_moltype)
        for name, seq in gapped_seqs_dict.items()
    }

    seq_data = {
        name: new_alignment.decompose_gapped_seq(seq, alphabet=dna_alphabet)[0]
        for name, seq in typed_data.items()
    }
    gap_data = {
        name: new_alignment.decompose_gapped_seq(seq, alphabet=dna_alphabet)[1]
        for name, seq in typed_data.items()
    }
    asd = new_alignment.AlignedSeqsData.from_seqs_and_gaps(
        seqs=seq_data, gaps=gap_data, alphabet=dna_alphabet
    )
    assert asd.align_len == 6
    assert asd.get_gapped_seq_str(seqid=seqid) == gapped_seqs_dict[seqid]


@pytest.mark.parametrize("data_type", (str, numpy.array, bytes))
def test_aligned_seqs_data_unequal_seqlens_raises(data_type, dna_alphabet, dna_moltype):
    """from_seqs for AlignedSeqsData should raise an error if sequences are of unequal length"""
    data = dict(
        seq1=make_typed("A-A", data_type=data_type, moltype=dna_moltype),
        seq2=make_typed("AAAAAAA--", data_type=data_type, moltype=dna_moltype),
    )
    with pytest.raises(ValueError):
        _ = new_alignment.AlignedSeqsData.from_seqs(data=data, alphabet=dna_alphabet)
    # directly creating an AlignedSeqsData object should also raise an error
    seq_data = {
        name: new_alignment.decompose_gapped_seq(seq, alphabet=dna_alphabet)[0]
        for name, seq in data.items()
    }
    gap_data = {
        name: new_alignment.decompose_gapped_seq(seq, alphabet=dna_alphabet)[1]
        for name, seq in data.items()
    }
    with pytest.raises(ValueError):
        _ = new_alignment.AlignedSeqsData.from_seqs_and_gaps(
            seqs=seq_data, gaps=gap_data, alphabet=dna_alphabet
        )


def test_from_seqs_and_gaps(dna_alphabet):
    # AlignedSeqsData should be able to be constructed from sequences and gap maps
    seqs = {"seq1": "ACCTA", "seq2": ""}
    gaps = {"seq1": numpy.array([[0, 1]]), "seq2": numpy.array([[0, 6]])}
    asd = new_alignment.AlignedSeqsData.from_seqs_and_gaps(
        seqs=seqs, gaps=gaps, alphabet=dna_alphabet
    )
    assert asd.get_seq_str(seqid="seq1") == "ACCTA"
    assert asd.get_gapped_seq_str(seqid="seq1") == "-ACCTA"
    assert asd.get_seq_str(seqid="seq2") == ""
    assert asd.get_gapped_seq_str(seqid="seq2") == "------"


def test_from_seqs_and_gaps_diff_seq_lens_raises(dna_alphabet):
    seqs = {"seq1": "ACCTA", "seq2": "A"}
    gaps = {"seq1": numpy.array([[0, 1]]), "seq2": numpy.array([[0, 1]])}
    with pytest.raises(ValueError):
        _ = new_alignment.AlignedSeqsData.from_seqs_and_gaps(
            seqs=seqs, gaps=gaps, alphabet=dna_alphabet
        )


def test_from_seqs_and_gaps_diff_keys_raises(dna_alphabet):
    seqs = {"seq1": "ACCTA", "seq2": "A"}
    gaps = {"seq1": numpy.array([[0, 1]]), "seq3": numpy.array([[0, 1]])}
    with pytest.raises(ValueError):
        _ = new_alignment.AlignedSeqsData.from_seqs_and_gaps(
            seqs=seqs, gaps=gaps, alphabet=dna_alphabet
        )


@pytest.mark.parametrize("seqid, i", (("seq1", 0), ("seq2", 1), ("seq3", 2)))
def test_from_names_and_array(dna_alphabet, seqid, i):
    names = ["seq1", "seq2", "seq3"]
    data = numpy.array([[0, 1, 2, 3], [3, 2, 1, 0], [4, 4, 4, 4]])
    asd = new_alignment.AlignedSeqsData.from_names_and_array(
        names=names, data=data, alphabet=dna_alphabet
    )
    assert asd.names == ("seq1", "seq2", "seq3")
    got = asd.get_gapped_seq_array(seqid=seqid)
    expect = data[i]
    assert numpy.array_equal(got, expect)


def test_from_names_and_array_empty_raises(dna_alphabet):
    names = []
    data = numpy.array([]).reshape(0, 0)

    with pytest.raises(ValueError):
        _ = new_alignment.AlignedSeqsData.from_names_and_array(
            names=names, data=data, alphabet=dna_alphabet
        )


def test_from_names_and_array_mismatched_length(dna_alphabet):
    names = ["seq1", "seq2"]
    data = numpy.array([[1, 0, 1], [0, 1, 0], [1, 1, 1]])

    with pytest.raises(ValueError):
        _ = new_alignment.AlignedSeqsData.from_names_and_array(
            names=names, data=data, alphabet=dna_alphabet
        )


def test_aligned_seqs_data_diff_keys_raises(dna_alphabet):
    """AlignedSeqsData expect identical keys in seqs and gaps"""
    seqs = dict(seq1=numpy.array([2, 1]), seq2=numpy.array([2, 0, 3, 1]))
    gaps = dict(
        seq1=numpy.array([[1, 3]], dtype=numpy.int32),
        seq3=numpy.array([[0, 1]], dtype=numpy.int32),
    )

    with pytest.raises(ValueError):
        _ = new_alignment.AlignedSeqsData.from_seqs_and_gaps(
            seqs=seqs, gaps=gaps, alphabet=dna_alphabet
        )
    # assert that it would work if we indeed had the same keys
    gaps["seq2"] = gaps.pop("seq3")
    asd = new_alignment.AlignedSeqsData.from_seqs_and_gaps(
        seqs=seqs, gaps=gaps, alphabet=dna_alphabet
    )
    got = asd.get_seq_str(seqid="seq1")
    assert got == "AC"


def test_aligned_seqs_data_omit_seqs_gaps_raises(dna_alphabet):
    with pytest.raises(ValueError):
        _ = new_alignment.AlignedSeqsData.from_seqs_and_gaps(
            seqs={}, gaps={}, alphabet=dna_alphabet
        )


def test_aligned_seqs_data_names(aligned_dict, dna_alphabet):
    got = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_dict, alphabet=dna_alphabet
    )
    assert isinstance(got, new_alignment.AlignedSeqsData)
    assert got.names == ("seq1", "seq2", "seq3", "seq4")


def test_aligned_seqs_data_len(aligned_dict, dna_alphabet):
    got = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_dict, alphabet=dna_alphabet
    )
    assert len(got) == len(aligned_dict["seq1"])
    assert got.align_len == len(aligned_dict["seq1"])


@pytest.mark.parametrize(
    "seqid, index", [("seq1", 0), ("seq2", 1), ("seq3", 2), ("seq4", 3)]
)
@pytest.mark.parametrize("moltype", ("dna_moltype", "rna_moltype"))
def test_aligned_seqs_data_getitem(seqid, index, aligned_array_dict, moltype, request):
    moltype = request.getfixturevalue(moltype)
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_array_dict, alphabet=moltype.degen_gapped_alphabet
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
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_array_dict, alphabet=moltype.degen_gapped_alphabet
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
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_array_dict, alphabet=moltype.degen_gapped_alphabet
    )
    got = ad.get_gapped_seq_array(seqid=seqid)  # get original data
    expect = aligned_array_dict[seqid]
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
@pytest.mark.parametrize("moltype", ("dna_moltype", "rna_moltype"))
def test_aligned_seqs_data_get_seq_str(aligned_array_dict, seqid, moltype, request):
    moltype = request.getfixturevalue(moltype)
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_array_dict, alphabet=moltype.degen_gapped_alphabet
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
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_array_dict, alphabet=moltype.degen_gapped_alphabet
    )
    got = ad.get_gapped_seq_str(seqid=seqid)
    expect = moltype.degen_gapped_alphabet.from_indices(aligned_array_dict[seqid])
    assert got == expect


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
@pytest.mark.parametrize("moltype", ("dna_moltype", "rna_moltype"))
def test_aligned_seqs_data_get_seq_bytes(aligned_array_dict, seqid, moltype, request):
    moltype = request.getfixturevalue(moltype)
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_array_dict, alphabet=moltype.degen_gapped_alphabet
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
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_array_dict, alphabet=moltype.degen_gapped_alphabet
    )
    got = ad.get_gapped_seq_bytes(seqid=seqid)
    raw = aligned_array_dict[seqid]
    expect = moltype.degen_gapped_alphabet.array_to_bytes(raw)
    assert got == expect


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_seqs_data_get_seq_length(seqid, aligned_dict, dna_alphabet):
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_dict, alphabet=dna_alphabet
    )
    raw = aligned_dict[seqid].replace("-", "")
    assert ad.get_seq_length(seqid) == len(raw)


def test_aligned_seqs_data_add_seqs(dna_alphabet):
    data = {"seq1": "ACGT-", "seq2": "ACG-T", "seq3": "A---T"}
    ad = new_alignment.AlignedSeqsData.from_seqs(data=data, alphabet=dna_alphabet)
    new_data = {"seq4": "ACGTT", "seq5": "ACG--", "seq6": "-----"}
    new_ad = ad.add_seqs(new_data)
    assert new_ad.names == ("seq1", "seq2", "seq3", "seq4", "seq5", "seq6")


def test_aligned_seqs_data_add_seqs_diff_lengths_raises(dna_alphabet):
    """adding sequences of different lengths should raise an error"""

    data = {"seq1": "ACGT-", "seq2": "ACG-T"}
    ad = new_alignment.AlignedSeqsData.from_seqs(data=data, alphabet=dna_alphabet)
    new_data = {"seq3": "AC", "seq4": "A-"}
    with pytest.raises(ValueError):
        _ = ad.add_seqs(new_data)

    new_data = {"seq3": "AAA", "seq4": "A-"}
    with pytest.raises(ValueError):
        _ = ad.add_seqs(new_data)


def test_aligned_seqs_data_add_seqs_diff_moltype_raises(dna_alphabet):
    """adding sequences of different moltype should raise an error"""
    data = {"seq1": "ACGT-", "seq2": "ACG-T"}  # DNA
    ad = new_alignment.AlignedSeqsData.from_seqs(data=data, alphabet=dna_alphabet)
    new_data = {"seq3": "ACGU-", "seq4": "ACG-U"}  # RNA
    with pytest.raises(new_alphabet.AlphabetError):
        _ = ad.add_seqs(new_data)


def test_aligned_seqs_data_add_seqs_duplicate_keys_raises(dna_alphabet):
    """AlignedSeqsData.add_seqs should raise an error if their are duplicated
    sequence names"""
    data = {"seq1": "ACGT-", "seq2": "ACG-T"}
    ad = new_alignment.AlignedSeqsData.from_seqs(data=data, alphabet=dna_alphabet)
    new_data = {"seq2": "ACGT-", "seq3": "ACT-T"}  # seq2 is duplicated
    with pytest.raises(ValueError):
        _ = ad.add_seqs(new_data)

    # if we set force_unique_keys=False, it should not raise an error
    new_ad = ad.add_seqs(new_data, force_unique_keys=False)
    assert new_ad.names == ("seq1", "seq2", "seq3")


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_seqs_data_get_aligned_view(aligned_dict, seqid, dna_alphabet):
    # str on an ADV should return the ungapped sequence
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_dict, alphabet=dna_alphabet
    )
    got = ad.get_view(seqid)
    assert got.parent == ad
    assert got.parent_len == ad.align_len
    assert str(got) == aligned_dict[seqid].replace("-", "")


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_data_view_array(aligned_array_dict, dna_alphabet, seqid):
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_array_dict, alphabet=dna_alphabet
    )
    view = ad.get_view(seqid)
    got = numpy.array(view)
    expect = aligned_array_dict[seqid][aligned_array_dict[seqid] != 4]  # remove gaps
    assert numpy.array_equal(got, expect)

    # directly accessing .array_value property should return the same result
    got = view.array_value
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_data_view_gapped_array(aligned_array_dict, dna_alphabet, seqid):
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_array_dict, alphabet=dna_alphabet
    )
    view = ad.get_view(seqid)
    got = view.gapped_array_value
    expect = aligned_array_dict[seqid]
    assert numpy.array_equal(got, expect)


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_data_view_str(aligned_dict, dna_alphabet, seqid):
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_dict, alphabet=dna_alphabet
    )
    view = ad.get_view(seqid)
    got = str(view)
    expect = aligned_dict[seqid].replace("-", "")
    assert got == expect

    # directly accessing .str_value property should return the same result
    got = view.str_value
    assert got == expect


def test_aligned_data_view_gapped_str_value(aligned_dict, dna_alphabet):
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_dict, alphabet=dna_alphabet
    )
    view = ad.get_view("seq1")
    got = view.gapped_str_value
    expect = aligned_dict["seq1"]
    assert got == expect


@pytest.mark.parametrize("seqid", ("seq1", "seq2", "seq3", "seq4"))
def test_aligned_data_view_bytes(aligned_array_dict, dna_alphabet, seqid):
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_array_dict, alphabet=dna_alphabet
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
def test_aligned_data_view_gapped_bytes_value(aligned_array_dict, dna_alphabet, seqid):
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_array_dict, alphabet=dna_alphabet
    )
    view = ad.get_view(seqid)
    got = view.gapped_bytes_value
    expect = dna_alphabet.array_to_bytes(aligned_array_dict[seqid])
    assert numpy.array_equal(got, expect)


@pytest.fixture
def simple_aln():
    return new_alignment.make_aligned_seqs(
        {"a": "T-C", "b": "---", "c": "AAA"}, moltype="dna"
    )


def test_alignment_array_seqs(simple_aln):
    got = simple_aln.array_seqs
    expect = numpy.array([[0, 4, 1], [4, 4, 4], [2, 2, 2]])
    assert numpy.array_equal(got, expect)


def test_alignment_array_seqs_renamed(simple_aln):
    # if the attribute ._array_seqs has been created (which happens when the
    # property .array_seqs is accessed), it should not be recreated when the
    # alignment is renamed
    orig_arr_seqs = simple_aln.array_seqs
    renamed = simple_aln.rename_seqs(renamer=lambda x: x.upper())
    assert renamed.array_seqs is orig_arr_seqs


def test_alignment_array_seqs_take_seqs(simple_aln):
    """an alignment which has been subset should return the correct array_seqs"""
    subset = simple_aln.take_seqs(["a", "c"])
    got = subset.array_seqs
    expect = numpy.array([[0, 4, 1], [2, 2, 2]])
    assert numpy.array_equal(got, expect)


def test_alignment_array_seqs_sliced(simple_aln):
    """an alignment which has been sliced should be reflected in the array_seqs"""
    # T-C
    # ---
    # AAA
    # **

    sliced = simple_aln[:2]
    got = sliced.array_seqs
    expect = numpy.array([[0, 4], [4, 4], [2, 2]])
    assert numpy.array_equal(got, expect)

    # T-C
    # ---
    # AAA
    # * *
    sliced = simple_aln[::2]
    got = sliced.array_seqs
    expect = numpy.array([[0, 1], [4, 4], [2, 2]])
    assert numpy.array_equal(got, expect)


def test_alignment_array_seqs_reverse_complement(simple_aln):
    """reverse complementing an alignment should be reflected in the array_seqs"""
    rc = simple_aln.rc()
    got = rc.array_seqs
    # grab the raw array data and reverse complement it
    raw = simple_aln.to_dict(as_array=True)
    expect = numpy.array([rc.moltype.complement(seq[::-1]) for seq in raw.values()])
    assert numpy.array_equal(got, expect)


def test_alignment_array_positions(simple_aln):
    got = simple_aln.array_positions
    expect = numpy.array([[0, 4, 2], [4, 4, 2], [1, 4, 2]])
    assert numpy.array_equal(got, expect)


def test_alignment_array_positions_take_positions(simple_aln):
    """an alignment which has been subset should return the correct array_positions"""
    subset = simple_aln.take_positions([0])
    got = subset.array_positions
    expect = numpy.array([[0, 4, 2]])
    assert numpy.array_equal(got, expect)


def test_alignment_array_positions_sliced(simple_aln):
    """an alignment which has been sliced should be reflected in the array_positions"""
    sliced = simple_aln[:2]
    got = sliced.array_positions
    expect = numpy.array([[0, 4, 2], [4, 4, 2]])
    assert numpy.array_equal(got, expect)

    sliced = simple_aln[::2]
    expect = numpy.array([[0, 4, 2], [1, 4, 2]])
    got = sliced.array_positions
    assert numpy.array_equal(got, expect)


def test_array_positions_reverse_complement(simple_aln):
    # orig:
    # T-C   041
    # ---   444
    # AAA   222

    # rc:
    # G-A   342
    # ---   444
    # TTT   000

    # iter columns should iter over columns in rc

    sliced = simple_aln[::-1]
    got = sliced.array_positions
    expect = numpy.array([[3, 4, 0], [4, 4, 0], [2, 4, 0]])
    assert numpy.array_equal(got, expect)


def test_alignment_init(aligned_dict, dna_moltype, dna_alphabet):
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_dict, alphabet=dna_alphabet
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


def test_make_aligned_seqs_aligned_seqs_data(aligned_dict, dna_alphabet):
    aligned_seqs_data = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_dict, alphabet=dna_alphabet
    )
    aln = new_alignment.make_aligned_seqs(aligned_seqs_data, moltype="dna")
    assert isinstance(aln, new_alignment.Alignment)
    assert aln.moltype.label == "dna"
    assert aln.names == ["seq1", "seq2", "seq3", "seq4"]
    assert aligned_dict == aln.to_dict()


def test_make_aligned_seqs_incompatible_moltype(aligned_dict, dna_alphabet):
    ad = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_dict, alphabet=dna_alphabet
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


def test_alignment_to_html_text_moltype():
    """exercising producing html for text moltype"""
    seqs = {"seq1": "ACG", "seq2": "-CT"}

    aln = new_alignment.make_aligned_seqs(seqs, moltype="text")
    got = aln.to_html(ref_name="longest")
    ref_row = (
        '<tr><td class="label">seq1</td>'
        '<td><span class="A_text">A</span>'
        '<span class="C_text">C</span>'
        '<span class="G_text">G</span></td></tr>'
    )
    other_row = (
        '<tr><td class="label">seq2</td>'
        '<td><span class="ambig_text">-</span>'
        '<span class="C_text">.</span>'
        '<span class="T_text">T</span></td></tr>'
    )

    assert ref_row in got
    assert other_row in got


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


def test_slice_align():
    """slicing alignment should work correctly"""
    data = {"seq1": "ACGACGACG", "seq2": "ACGACGACG", "seq3": "ACGACGACG"}
    alignment = new_alignment.make_aligned_seqs(data, moltype="dna")
    sub_align = alignment[2:5]
    assert isinstance(sub_align, new_alignment.Alignment)
    expect = {"seq1": "GAC", "seq2": "GAC", "seq3": "GAC"}
    assert sub_align.to_dict() == expect
    # slice third positions
    sub_align = alignment[2::3]
    expect = {"seq1": "GGG", "seq2": "GGG", "seq3": "GGG"}
    assert sub_align.to_dict() == expect


def test_slice_align_info():
    """slicing alignment preserves info attribute"""
    data = {"seq1": "ACGACGACG", "seq2": "ACGACGACG", "seq3": "ACGACGACG"}
    alignment = new_alignment.make_aligned_seqs(
        data, info={"key": "value"}, moltype="dna"
    )
    sub_align = alignment[2:5]
    assert len(sub_align) == 3
    assert sub_align.info["key"] == "value"


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


def test_alignment_quality():
    """check alignment method correctly invokes underlying app"""
    aln = new_alignment.make_aligned_seqs(
        ["AAAC", "ACGC", "AGCC", "A-TC"], moltype="dna"
    )
    got = aln.alignment_quality(equifreq_mprobs=False)
    expect = (
        2 * numpy.log2(1 / 0.4)
        + numpy.log2(1 / (4 * 0.4))
        + (1 / 2) * numpy.log2(1 / (8 / 15))
        + (1 / 4) * numpy.log2(1 / (4 / 15))
    )
    assert numpy.all(got == expect)


def test_variable_positions():
    """correctly identify variable positions"""
    new_seqs = {"A": "-CG-C", "B": "ACAA?", "C": "GCGAC"}
    aln = new_alignment.make_aligned_seqs(new_seqs, moltype="dna")
    assert aln.variable_positions(include_gap_motif=True) == (0, 2, 3, 4)
    assert aln.variable_positions(include_gap_motif=False) == (0, 2)
    new_seqs = {"A": "GCGAC", "B": "GCGAC", "C": "GCGAC"}
    aln = new_alignment.make_aligned_seqs(new_seqs, moltype="dna")
    assert aln.variable_positions(include_gap_motif=True) == ()
    assert aln.variable_positions(include_gap_motif=False) == ()
    new_seqs = {"A": "-CG?C-", "B": "ACAAYA", "C": "GCGACA"}
    aln = new_alignment.make_aligned_seqs(new_seqs, moltype="dna")
    assert aln.variable_positions(include_gap_motif=False, include_ambiguity=True) == (
        0,
        2,
        3,
        4,
    )
    assert aln.variable_positions(include_gap_motif=True, include_ambiguity=True) == (
        0,
        2,
        3,
        4,
        5,
    )


def test_variable_positions_motif_length():
    """correctly identify variable positions"""
    #                 *  - ?
    new_seqs = {"A": "-CG-CA", "B": "ACGACA", "C": "GCGACY"}
    aln = new_alignment.make_aligned_seqs(new_seqs, moltype="dna")
    # only the first dinucleotide is variable if gaps and ambigs disallowed
    assert aln.variable_positions(
        motif_length=2, include_gap_motif=False, include_ambiguity=False
    ) == (0, 1)
    # if gaps are allowed, the second dinucleotide is also variable
    assert aln.variable_positions(
        motif_length=2, include_gap_motif=True, include_ambiguity=False
    ) == (0, 1, 2, 3)
    # if no gaps but ambig allowed, the third dinucleotide is also variable
    assert aln.variable_positions(
        motif_length=2, include_gap_motif=False, include_ambiguity=True
    ) == (0, 1, 4, 5)
    # if gaps and ambig allowed all are variable
    assert aln.variable_positions(
        motif_length=2, include_gap_motif=True, include_ambiguity=True
    ) == (0, 1, 2, 3, 4, 5)
    # variable remains modulo motif_length
    new_seqs = {"A": "-CG-CAA", "B": "ACGACAC", "C": "GCGACYG"}
    aln = new_alignment.make_aligned_seqs(new_seqs, moltype="dna")
    assert aln.variable_positions(
        motif_length=2, include_gap_motif=True, include_ambiguity=True
    ) == (0, 1, 2, 3, 4, 5)


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

    # if we slice the alignment, the gap array should reflect the slice
    sliced = aln[:2]
    got = sliced.get_gap_array()
    expect = expect[:, :2]
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


def test_alignment_get_feature():
    aln = new_alignment.make_aligned_seqs(
        {"x": "-AAAAAAAAA", "y": "TTTT--CCCT"}, moltype="dna"
    )

    db = aln.annotation_db
    db.add_feature(seqid="y", biotype="exon", name="A", spans=[(5, 8)])
    feat = list(aln.get_features(seqid="y", biotype="exon", on_alignment=False))[0]
    assert feat.get_slice().to_dict() == dict(x="AAA", y="CCT")


def test_alignment_distance_matrix():
    """Alignment distance_matrix should produce correct scores"""
    data = dict([("s1", "ACGTACGTA"), ("s2", "GTGTACGTA")])
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    dists = aln.distance_matrix(calc="hamming", show_progress=False)
    assert dists == {("s1", "s2"): 2.0, ("s2", "s1"): 2.0}
    # and for protein
    aa = aln.get_translation()
    dists = aa.distance_matrix(calc="hamming")
    assert dists == {("s1", "s2"): 1.0, ("s2", "s1"): 1.0}

    # when there are invalid data
    data = dict(
        seq1="AGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
        seq2="TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
        seq3="TACAAAAAAAAGGGGCGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
    )

    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    with pytest.raises(ArithmeticError):
        # default settings cause an exception
        _ = aln.distance_matrix(calc="paralinear")
    # but setting drop_invalid=False allows calc
    dists = aln.distance_matrix(calc="paralinear", drop_invalid=True)
    assert dists is None


def test_alignment_sample_with_replacement():
    # test with replacement -- just verify that it rnus
    alignment = new_alignment.make_aligned_seqs(
        {"seq1": "GATC", "seq2": "GATC"}, moltype="dna"
    )
    sample = alignment.sample(n=100, with_replacement=True)
    assert len(sample) == 100
    # ensure that sampling with replacement works on single col alignment
    alignment1 = new_alignment.make_aligned_seqs(
        {"seq1": "A", "seq2": "A"}, moltype="dna"
    )
    result = alignment1.sample(with_replacement=True)
    assert len(result) == 1


def test_alignment_sample_tuples():
    ##### test with motif size != 1 #####
    alignment = new_alignment.make_aligned_seqs(
        {
            "seq1": "AACCDDEEFFGGHHIIKKLLMMNNPP",
            "seq2": "AACCDDEEFFGGHHIIKKLLMMNNPP",
        },
        moltype="protein",
    )
    # ensure length correct
    sample = alignment.sample(n=10, motif_length=2)
    assert len(sample), 20
    # test columns alignment preserved
    seqs = list(sample.to_dict().values())
    assert seqs[0] == seqs[1]
    # ensure each char occurs twice as sampling dinucs without replacement
    for char in seqs[0]:
        assert seqs[0].count(char) == 2

    assert all(len(set(seqs[0][i : i + 2])) == 1 for i in range(0, len(seqs[0]), 2))


def test_alignment_sample_sliced():
    """sample should only return characters from current view"""
    # not a deterministic test, but has a small false positive rate of 1/1001
    data = {"seq1": "A" * 1000 + "C"}
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    sliced = aln[-1:]
    got = sliced.sample(n=1)
    assert got.to_dict() == {"seq1": "C"}

    rc = aln.rc()
    sliced = rc[:1]
    got = sliced.sample(n=1)
    assert got.to_dict() == {"seq1": "G"}


def test_alignment_take_positions():
    """SequenceCollection take_positions should return new alignment w/ specified pos"""
    gaps = new_alignment.make_aligned_seqs(
        {"a": "AAAAAAA", "b": "A--A-AA", "c": "AA-----"}, moltype="dna"
    )
    assert gaps.take_positions([5, 4, 0]).to_dict() == {
        "a": "AAA",
        "b": "A-A",
        "c": "--A",
    }

    # should be able to negate
    assert gaps.take_positions([5, 4, 0], negate=True).to_dict() == {
        "a": "AAAA",
        "b": "--AA",
        "c": "A---",
    }


def test_alignment_take_positions_sliced():
    #       012345
    #       -TAA--
    #       TTAGGG
    # slice: ***
    # pos:    **

    data = {
        "seq1": "-TAA--",
        "seq2": "TTAGGG",
    }
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    sliced = aln[1:4]
    got = sliced.take_positions([1, 2]).to_dict()
    expect = {"seq1": "AA", "seq2": "AG"}
    assert got == expect

    rev = aln.rc()
    got = rev.take_positions([1, 2]).to_dict()
    expect = {"seq1": "-T", "seq2": "CC"}
    assert got == expect


def test_alignment_take_positions_info():
    aln = new_alignment.make_aligned_seqs(
        {"a": "AAAAAAA", "b": "A--A-AA", "c": "AA-----"},
        moltype="dna",
        info={"key": "value"},
    )
    tps = aln.take_positions([5, 4, 0])
    assert tps.info["key"] == "value"


def test_take_positions_if():
    """take_positions_if should return cols where f(col) is True"""

    gaps = new_alignment.make_aligned_seqs(
        {"a": "AAAAAAA", "b": "A--A-AA", "c": "AA-----"}, moltype="dna"
    )

    def gap_1st(x):
        return x[0] == "-"

    def gap_2nd(x):
        return x[1] == "-"

    def gap_3rd(x):
        return x[2] == "-"

    def is_list(x):
        return isinstance(x, list)

    got = gaps.take_positions_if(gap_1st).to_dict()
    assert got == {"a": "", "b": "", "c": ""}

    got = gaps.take_positions_if(gap_2nd).to_dict()
    assert got == {"a": "AAA", "b": "---", "c": "A--"}

    got = gaps.take_positions_if(gap_3rd).to_dict()
    assert got == {"a": "AAAAA", "b": "-A-AA", "c": "-----"}

    got = gaps.take_positions_if(is_list).to_dict()
    assert got == gaps.to_dict()

    # should be able to negate
    got = gaps.take_positions_if(gap_2nd, negate=True).to_dict()
    assert got == {"a": "AAAA", "b": "AAAA", "c": "A---"}

    got = gaps.take_positions_if(gap_3rd, negate=True).to_dict()
    assert got == {"a": "AA", "b": "A-", "c": "AA"}


def _make_and_filter(raw, expected, motif_length, drop_remainder):
    # a simple filter func
    aln = new_alignment.make_aligned_seqs(raw, moltype="dna", info={"key": "value"})

    def func(x):
        return "-" not in "".join([str(s) for s in x])

    result = aln.filtered(
        func,
        motif_length=motif_length,
        warn=False,
        drop_remainder=drop_remainder,
    )
    assert result.to_dict() == expected
    assert result.info["key"] == "value"


def test_filtered():
    """filtered should return new alignment with positions consistent with
    provided callback function"""
    # a simple filter option
    raw = {"a": "ACGACGACG", "b": "CCC---CCC", "c": "AAAA--AAA"}
    _make_and_filter(raw, {"a": "ACGACG", "b": "CCCCCC", "c": "AAAAAA"}, 1, True)
    # check with motif_length = 2
    _make_and_filter(raw, {"a": "ACAC", "b": "CCCC", "c": "AAAA"}, 2, True)
    # check with motif_length = 3
    _make_and_filter(raw, {"a": "ACGACG", "b": "CCCCCC", "c": "AAAAAA"}, 3, True)


def test_filtered_drop_remainder():
    """filter allows dropping"""
    raw = {"a": "ACGACGACG", "b": "CCC---CCC", "c": "AAAA--AAA"}
    aln = new_alignment.make_aligned_seqs(raw, moltype="dna")

    def func(x):
        return "-" not in "".join([str(s) for s in x])

    got = aln.filtered(func, motif_length=1, warn=False)
    assert len(got) == 6
    # raises an assertion if the length is not modulo

    with pytest.raises(ValueError):
        # because alignment not modulo 2
        got = aln.filtered(func, motif_length=2, drop_remainder=False)
    got = aln.filtered(func, motif_length=2, drop_remainder=True, warn=False)
    assert len(got) == 4


def test_no_degenerates():
    """no_degenerates correctly excludes columns containing IUPAC ambiguity codes"""
    data = {
        "s1": "AAA CCC GGG TTT".replace(" ", ""),
        "s2": "CCC GGG T-T AAA".replace(" ", ""),
        "s3": "GGR YTT AAA CCC".replace(" ", ""),
    }
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")

    # motif length of 1, defaults - no gaps allowed
    result = aln.no_degenerates().to_dict()
    expect = {
        "s1": "AA CC GG TTT".replace(" ", ""),
        "s2": "CC GG TT AAA".replace(" ", ""),
        "s3": "GG TT AA CCC".replace(" ", ""),
    }
    assert result == expect

    # allow gaps
    result = aln.no_degenerates(allow_gap=True).to_dict()
    expect = {
        "s1": "AA CC GGG TTT".replace(" ", ""),
        "s2": "CC GG T-T AAA".replace(" ", ""),
        "s3": "GG TT AAA CCC".replace(" ", ""),
    }
    assert result == expect

    # motif length of 3, defaults - no gaps allowed
    result = aln.no_degenerates(motif_length=3, allow_gap=False).to_dict()
    expect = {
        "s1": "TTT".replace(" ", ""),
        "s2": "AAA".replace(" ", ""),
        "s3": "CCC".replace(" ", ""),
    }
    assert result == expect

    # allow gaps
    result = aln.no_degenerates(motif_length=3, allow_gap=True).to_dict()
    expect = {
        "s1": "GGG TTT".replace(" ", ""),
        "s2": "T-T AAA".replace(" ", ""),
        "s3": "AAA CCC".replace(" ", ""),
    }
    assert result == expect

    # for length non-divisible by motif_length
    data = {
        "s1": "AAA CCC GGG TTT T".replace(" ", ""),
        "s2": "CCC GGG T-T AAA A".replace(" ", ""),
        "s3": "GGR YTT AAA CCC C".replace(" ", ""),
    }
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    result = aln.no_degenerates(motif_length=3, allow_gap=False).to_dict()
    expect = {
        "s1": "TTT".replace(" ", ""),
        "s2": "AAA".replace(" ", ""),
        "s3": "CCC".replace(" ", ""),
    }


def test_no_degenerates_bad_moltype_raises():
    """no_degenerates should raise an error if moltype has no degenerate symbols"""
    aln = new_alignment.make_aligned_seqs({"s1": "ACGT", "s2": "ACGT"}, moltype="text")
    with pytest.raises(new_moltype.MolTypeError):
        _ = aln.no_degenerates()


def test_omit_gap_pos_motif_length():
    """consistency with different motif_length values"""
    data = {
        "seq1": "CAGGTCGACCTCGGC---------CACGAC",
        "seq2": "CAGATCGACCTCGGC---------CACGAC",
        "seq3": "CAGATCGACCTCGGT---------CACGAT",
        "seq4": "CAGATCGACCTCGGCGAACACGGCCATGAT",
        "seq5": "CCGATCGACATGGGC---------CACGAT",
        "seq6": "GCC---------------------------",
    }
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    got1 = aln.omit_gap_pos(motif_length=1)
    got3 = aln.omit_gap_pos(motif_length=3)
    assert len(got3) == len(got1)
    assert got3.to_dict() == got1.to_dict()


def test_omit_gap_pos():
    """Alignment omit_gap_pos should return alignment w/o positions of gaps"""
    aln = new_alignment.make_aligned_seqs(
        {"a": "--A-BC-", "b": "-CB-A--", "c": "--D-EF-"}, moltype="protein"
    )
    # first, check behavior when we're just acting on the cols (and not
    # trying to delete the naughty seqs).

    # default should strip out cols that are 100% gaps
    result = aln.omit_gap_pos()
    assert result.to_dict() == {"a": "-ABC", "b": "CBA-", "c": "-DEF"}
    # if allowed_gap_frac is 1, shouldn't delete anything
    assert aln.omit_gap_pos(1).to_dict() == {
        "a": "--A-BC-",
        "b": "-CB-A--",
        "c": "--D-EF-",
    }

    # if allowed_gap_frac is 0, should strip out any cols containing gaps
    assert aln.omit_gap_pos(0).to_dict() == {"a": "AB", "b": "BA", "c": "DE"}
    # intermediate numbers should work as expected
    assert aln.omit_gap_pos(0.4).to_dict() == {"a": "ABC", "b": "BA-", "c": "DEF"}
    assert aln.omit_gap_pos(0.7).to_dict() == {"a": "-ABC", "b": "CBA-", "c": "-DEF"}

    # when we increase the number of sequences to 6, more differences
    # start to appear.
    new_aln = aln.add_seqs({"d": "-------", "e": "XYZXYZX", "f": "AB-CDEF"})

    # if no gaps are allowed, we get None
    got = new_aln.omit_gap_pos(0)
    assert not len(got)


def test_omit_gap_pos_no_gap_moltype():
    """if moltype does not support gap, just returns self"""
    aln = new_alignment.make_aligned_seqs(
        {"a": "--A-BC-", "b": "-CB-A--", "c": "--D-EF-"}, moltype="bytes"
    )
    got = aln.omit_gap_pos()
    assert got is aln


@pytest.fixture(scope="session")
def brca1_data():
    return load_aligned_seqs("data/brca1.fasta", moltype="dna", new_type=True).to_dict()


@pytest.mark.parametrize("calc", ("hamming", None))
def test_alignment_quick_tree(calc, brca1_data):
    """quick tree method returns tree"""
    aln = new_alignment.make_aligned_seqs(brca1_data, moltype="dna")[:100]
    aln = aln.take_seqs(["Human", "Rhesus", "HowlerMon", "Galago", "Mouse"])
    kwargs = dict(show_progress=False)
    if calc:
        kwargs["calc"] = calc

    # bootstrap
    tree = aln.quick_tree(bootstrap=2, **kwargs)
    assert set(tree.get_tip_names()) == set(aln.names)
    types = {
        type(float(edge.params["support"]))
        for edge in tree.preorder()
        if not edge.is_root()
    }
    assert types == {float}


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


@pytest.mark.parametrize("rc", [True, False])
@pytest.mark.parametrize(
    "func", (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs)
)
def test_alignment_offset_propagation(aligned_dict, func, rc):
    # providing an offset should set the offset on precisely the specified seq
    aln = func(aligned_dict, moltype="dna", offset={"seq1": 10})
    seq = aln.get_seq("seq1").rc() if rc else aln.get_seq("seq1")
    assert seq._seq.offset == 10
    assert seq.annotation_offset == 10

    seq = aln.get_seq("seq2").rc() if rc else aln.get_seq("seq2")
    assert seq._seq.offset == 0
    assert seq.annotation_offset == 0


def test_alignment_offset_sliced(aligned_dict):
    aln = new_alignment.make_aligned_seqs(
        aligned_dict, moltype="dna", offset={"seq1": 10}
    )
    sliced = aln[2:]
    seq = sliced.get_seq("seq1")
    assert seq._seq.offset == 10
    assert seq.annotation_offset == 12


@pytest.mark.parametrize("name", ("s1", "s2", "s3"))
def test_get_seq_with_sliced_aln(name):
    seqs = {
        "s1": "GTTGAAGTAGTAGAAGTTCCAAATAATGAA",
        "s2": "GTG------GTAGAAGTTCCAAATAATGAA",
        "s3": "GCTGAAGTAGTGGAAGTTGCAAAT---GAA",
    }
    aln = new_alignment.make_aligned_seqs(seqs, moltype="dna")
    start, stop = 1, 5
    a1 = aln[start:stop]

    seq = a1.get_seq(name)
    assert isinstance(seq, new_sequence.Sequence), seq

    got = str(seq)
    expect = seqs[name][start:stop].replace("-", "")
    assert got == expect, (got, expect)


@pytest.mark.parametrize("name", ("s1", "s2", "s3"))
def test_get_seq_with_sliced_rced_aln(name):
    seqs = {
        "s1": "GTTGAAGTAGTAGAAGTTCCAAATAATGAA",
        "s2": "GTG------GTAGAAGTTCCAAATAATGAA",
        "s3": "GCTGAAGTAGTGGAAGTTGCAAAT---GAA",
    }
    aln = new_alignment.make_aligned_seqs(seqs, moltype="dna")
    start, stop = 1, 5
    a1 = aln[start:stop]
    a1 = a1.rc()
    got = str(a1.get_seq(name))

    dna = new_moltype.get_moltype("dna")
    expect = dna.complement(seqs[name][start:stop].replace("-", ""))[::-1]
    assert got == expect, (got, expect)


@pytest.mark.parametrize("name", ("s1", "s2", "s3"))
def test_get_seq_with_sliced_aln_multiple_spans(name):
    seqs = {  # the sliced seq has:
        "s1": "GTTGA--TAGTAGAAGTTCCAAATAATGAA",  # span gap span
        "s2": "G----TT------AAGTTCCAAATAATGAA",  # gap span gap
        "s3": "G--GA--TA--GGAAGTTGCAAAT---GAA",  # gap span gap span gap
    }
    aln = new_alignment.make_aligned_seqs(seqs, moltype="dna")
    start, stop = 1, 10
    a1 = aln[start:stop]
    seq = a1.get_seq(name)
    assert isinstance(seq, new_sequence.Sequence), seq

    expect = seqs[name][start:stop].replace("-", "")
    got = str(seq)
    assert got == expect, (got, expect)


@pytest.mark.parametrize("name", ("s1", "s2", "s3"))
def test_get_seq_with_sliced_rced_aln_multiple_spans(name):
    seqs = {  # the sliced seq has:
        "s1": "GTTGA--TAGTAGAAGTTCCAAATAATGAA",  # span gap span
        "s2": "G----TT------AAGTTCCAAATAATGAA",  # gap span gap
        "s3": "G--GA--TA--GGAAGTTGCAAAT---GAA",  # gap span gap span gap
    }
    aln = new_alignment.make_aligned_seqs(seqs, moltype="dna")
    start, stop = 1, 10
    a1 = aln[start:stop]
    a1 = a1.rc()
    got = str(a1.get_seq(name))
    dna = new_moltype.get_moltype("dna")
    expect = dna.complement(seqs[name][start:stop].replace("-", ""))[::-1]
    assert got == expect, (got, expect)


@pytest.mark.parametrize("name", ("s1", "s2", "s3"))
def test_get_gapped_seq_with_sliced_aln(name):
    seqs = {
        "s1": "G-TG---TAGTAGAAGTTCCAAATAATGAA",
        "s2": "GTG------GTAGAAGTTCCAAATAATGAA",
        "s3": "GC--AAGTAGTGGAAGTTGCAAAT---GAA",
    }
    aln = new_alignment.make_aligned_seqs(seqs, moltype="dna")
    start, stop = 1, 10
    a1 = aln[start:stop]

    seq = a1.get_gapped_seq(name)
    assert isinstance(seq, new_sequence.Sequence), seq

    got = str(seq)
    expect = seqs[name][start:stop]
    assert got == expect, (got, expect)


@pytest.mark.parametrize("name", ("s1", "s2", "s3"))
def test_aln_rev_slice(name):
    seqs = {
        "s1": "AAGGTTCC",
        "s2": "AAGGTTCC",
        "s3": "AAGGTTCC",
    }

    aln = new_alignment.make_aligned_seqs(seqs, moltype="dna")
    got = aln[5:1]
    assert not got

    seq = got.get_gapped_seq(name)
    assert not str(seq)

    seq = got.get_seq(name)
    assert not str(seq)


## Tests of _IndexableSeqs


@pytest.fixture
def names_seqs():
    names = ["seq1", "seq2", "seq3"]
    seqs = ["GGGTAC", "GTTTGC", "ACGTAC"]
    return names, seqs


@pytest.mark.parametrize(
    "func", (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs)
)
@pytest.mark.parametrize("index", (True, False))
def test_indexing_seqs_prop(names_seqs, func, index):
    names, seqs = names_seqs
    raw = dict(zip(names, seqs))
    obj = func(raw, moltype="dna")
    got = obj.seqs[0 if index else "seq1"]
    assert str(got) == raw["seq1"]


def test_sequence_collection_indexing_seqs_repr(names_seqs):
    names, seqs = names_seqs
    raw = dict(zip(names, seqs))
    obj = new_alignment.make_unaligned_seqs(raw, moltype="dna")
    got = repr(obj.seqs)
    class_name = obj.seqs[0].__class__.__name__
    expect = f"({class_name}({seqs[0]}), + {len(names)-1} seqs)"
    assert got == expect


def test_alignment_indexing_seqs_repr(names_seqs):
    names, seqs = names_seqs
    raw = dict(zip(names, seqs))
    obj = new_alignment.make_aligned_seqs(raw, moltype="dna")
    got = repr(obj.seqs)
    class_name = "Aligned"
    expect = f"({class_name}(name={names[0]!r}, seq={seqs[0]!r}, moltype={obj.moltype.name!r}), + {len(names)-1} seqs)"
    assert got == expect


def test_aligned_repr():
    seqs = {
        "s1": "G-TG---TAGTAGAAGTTCCAAATAATGAA",
        "s2": "GTG------GTAGAAGTTCCAAATAATGAA",
        "s3": "GC--AAGTAGTGGAAGTTGCAAAT---GAA",
    }
    aln = new_alignment.make_aligned_seqs(seqs, moltype="dna")
    got = repr(aln.seqs["s1"])
    expect = f"Aligned(name='s1', seq='G-TG---... {len(seqs['s1'])}', moltype='dna')"
    assert got == expect


@pytest.mark.parametrize(
    "func", (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs)
)
def test_indexing_seqs_iter(names_seqs, func):
    names, seqs = names_seqs
    raw = dict(zip(names, seqs))
    obj = func(raw, moltype="dna")
    # exercise __iter__
    got = list(map(str, obj.seqs))
    expect = list(seqs)
    assert got == expect


@pytest.mark.parametrize(
    "gapped_seq",
    (
        "ACGGCTA",  # No gaps
        "AC--GT-AGC",  # With gaps
        "------",  # All gaps
        # Edge cases
        "-A-C-G-T",
        "A-C-G-T-",
    ),
)
def test_gapped_seq_round_trip(gapped_seq):
    gapped_seq = new_moltype.DNA.gapped_alphabet.to_indices(gapped_seq)
    # split into components
    ungapped_seq, gaps = new_alignment.decompose_gapped_seq_array(gapped_seq, 4)

    # Recreate the gapped sequence from the ungapped sequence and gaps
    recreated_gapped_seq = new_alignment.compose_gapped_seq(ungapped_seq, gaps, 4)

    # Test the output of gapped_seq_from_components against the original sequence
    numpy.testing.assert_array_equal(recreated_gapped_seq, gapped_seq)


def test_asd_get_gapped_seq(aligned_dict, dna_alphabet):
    a_slice = slice(None, 4, 2)
    seqid = "seq4"
    orig_seq = aligned_dict[seqid][a_slice]
    orig_array = dna_alphabet.to_indices(orig_seq)
    orig_revd = orig_array[::-1]

    asd = new_alignment.AlignedSeqsData.from_seqs(
        data=aligned_dict, alphabet=new_moltype.DNA.most_degen_alphabet()
    )
    view = asd.get_view(seqid)
    fwd = view[a_slice]
    rev = fwd[::-1]
    fwd_gapped = fwd.gapped_array_value
    numpy.testing.assert_array_equal(fwd_gapped, orig_array)
    rev_gapped = rev.gapped_array_value
    numpy.testing.assert_array_equal(rev_gapped, orig_revd)


def test_make_gap_filter():
    """make_gap_filter returns f(seq) -> True if aligned ok w/ query"""
    RNA = new_moltype.get_moltype("rna")
    s1 = RNA.make_seq(seq="UC-----CU---C")
    s3 = RNA.make_seq(seq="UUCCUUCUU-UUC")
    s4 = RNA.make_seq(seq="UU-UUUU-UUUUC")
    # check that the behavior is ok for gap runs
    f1 = new_alignment.make_gap_filter(s1, 0.9, 5)
    f3 = new_alignment.make_gap_filter(s3, 0.9, 5)
    # Should return False since s1 has gap run >= 5 with respect to s3
    assert f3(s1) == False
    # Should return False since s3 has an insertion run >= 5 to s1
    assert f1(s3) == False
    # Should retun True since s4 does not have a long enough gap or ins run
    assert f3(s4) == True
    f3 = new_alignment.make_gap_filter(s3, 0.9, 6)
    assert f3(s1) == True

    # Check that behavior is ok for gap_fractions
    f1 = new_alignment.make_gap_filter(s1, 0.5, 6)
    f3 = new_alignment.make_gap_filter(s3, 0.5, 6)
    # Should return False since 0.53% of positions are diff for gaps
    assert f3(s1) == False
    assert f1(s3) == False
    assert f3(s4) == True


@pytest.fixture(scope="session")
def codon_and_aa_alns():
    import cogent3

    data = dict(s1="ATG --- --- GAT --- AAA", s2="ATG CAA TCG AAT GAA ATA")
    dna = cogent3.make_aligned_seqs(
        {n: s.replace(" ", "") for n, s in data.items()}, moltype="dna", new_type=True
    )
    aa = dna.get_translation()
    return dna, aa


def test_alignment_apply_scaled_gaps_aa_to_codon(codon_and_aa_alns):
    codon, aa = codon_and_aa_alns
    ungapped = codon.degap()
    scaled = aa.apply_scaled_gaps(ungapped, aa_to_codon=True)
    assert scaled.to_dict() == codon.to_dict()


def test_alignment_apply_scaled_gaps_codon_to_aa(codon_and_aa_alns):
    codon, aa = codon_and_aa_alns
    ungapped = aa.degap()
    scaled = codon.apply_scaled_gaps(ungapped, aa_to_codon=False)
    assert scaled.to_dict() == aa.to_dict()


def test_alignment_apply_scaled_gaps_invalid_seqlen(codon_and_aa_alns):
    codon, aa = codon_and_aa_alns
    codon = codon[:-2]
    ungapped = codon.degap()
    with pytest.raises(ValueError):
        aa.apply_scaled_gaps(ungapped, aa_to_codon=True)


def test_alignment_apply_scaled_gaps_aa2codon_invalid_moltype(codon_and_aa_alns):
    codon, aa = codon_and_aa_alns
    codon = codon.to_moltype("text")
    ungapped = codon.degap()
    with pytest.raises(ValueError):
        aa.apply_scaled_gaps(ungapped, aa_to_codon=True)


def test_alignment_apply_scaled_gaps_codon2aa_invalid_moltype(codon_and_aa_alns):
    codon, aa = codon_and_aa_alns
    ungapped = aa.degap()
    codon = codon.to_moltype("text")
    with pytest.raises(ValueError):
        codon.apply_scaled_gaps(ungapped, aa_to_codon=False)


def test_alignment_copy(simple_aln):
    got = simple_aln.copy()
    # mutable data structures should be different IDs
    assert got._name_map is not simple_aln._name_map
    assert got.info is not simple_aln.info
    assert got.annotation_db is not simple_aln.annotation_db
    # immutable data structures should be the same object
    assert got._slice_record is simple_aln._slice_record
    assert got._seqs_data is simple_aln._seqs_data
    # sliced data should have same underlying length
    sl = simple_aln[:2].copy()
    assert sl._seqs_data.align_len == simple_aln._seqs_data.align_len

    # renamed seqs should be propagated
    renamed = simple_aln.rename_seqs(renamer=lambda x: x.upper())
    copied = renamed.copy()
    assert set(renamed.names) != set(simple_aln.names)
    assert set(copied.names) == set(renamed.names)
    # and can get a sequence with a new name
    assert str(copied.seqs["A"]) == str(renamed.seqs["A"])


def test_alignment_deepcopy(simple_aln):
    got = simple_aln.deepcopy()
    # all data structures should be different IDs
    assert got._name_map is not simple_aln._name_map
    assert got.info is not simple_aln.info
    assert got.annotation_db is not simple_aln.annotation_db
    assert got._slice_record is not simple_aln._slice_record
    assert got._seqs_data is not simple_aln._seqs_data

    # sliced data should have different underlying data length
    sl = simple_aln[:2].deepcopy()
    assert sl._seqs_data.align_len == 2 != simple_aln._seqs_data.align_len

    # renamed seqs should be propagated
    renamed = simple_aln.rename_seqs(renamer=lambda x: x.upper())
    copied = renamed.deepcopy()
    assert set(renamed.names) != set(simple_aln.names)
    assert set(copied.names) == set(renamed.names)
    # and can get a sequence with a new name
    assert str(copied.seqs["A"]) == str(renamed.seqs["A"])


@pytest.mark.parametrize(
    "mk_cls",
    [
        new_alignment.make_unaligned_seqs,
        new_alignment.make_aligned_seqs,
    ],
)
def test_collections_equal(aligned_dict, mk_cls):
    coll1 = mk_cls(aligned_dict, moltype="dna")
    coll2 = mk_cls(aligned_dict, moltype="dna")
    assert coll1 is not coll2
    assert coll1 == coll2


@pytest.mark.parametrize(
    "mk_cls",
    [
        new_alignment.make_unaligned_seqs,
        new_alignment.make_aligned_seqs,
    ],
)
@pytest.mark.parametrize(
    "kwargs",
    [dict(moltype="protein"), dict(is_reversed=True), dict(name_map=dict(A="a")), {}],
)
def test_collections_not_equal(aligned_dict, mk_cls, kwargs):
    coll1 = mk_cls(aligned_dict, moltype="dna")
    aligned_dict = aligned_dict if kwargs else {**aligned_dict, **{"seq1": "TTTTTT"}}
    kwargs = {**dict(moltype="dna"), **kwargs}
    coll2 = mk_cls(aligned_dict, **kwargs)
    assert coll1 != coll2


def test_alignment_not_equal_sliced(aligned_dict):
    coll1 = new_alignment.make_aligned_seqs(aligned_dict, moltype="dna")
    coll2 = coll1[:2]
    assert coll1 != coll2


@pytest.mark.parametrize(
    "type1,type2",
    [
        (new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs),
        (new_alignment.make_aligned_seqs, new_alignment.make_unaligned_seqs),
    ],
)
def test_alignment_not_equal_types(aligned_dict, type1, type2):
    coll1 = type1(aligned_dict, moltype="dna")
    coll2 = type2(aligned_dict, moltype="dna")
    assert coll1 != coll2


def test_alignment_to_rich_dict_sliced():
    # slices are realised
    data = {
        "seq1": "ATG",
        "seq2": "TGA",
    }
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    sliced = aln[1:]
    assert sliced._slice_record.parent_len == 3
    rd = sliced.to_rich_dict()
    got = deserialise_object(rd)
    assert got.to_dict() == {"seq1": "TG", "seq2": "GA"}
    assert got._slice_record.parent_len == 2


def test_alignment_to_rich_dict_round_trip():
    data = {
        "seq1": "ATG",
        "seq2": "TGA",
    }
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    rd = aln.to_rich_dict()
    got = deserialise_object(rd)
    assert isinstance(got, new_alignment.Alignment)
    assert got.to_dict() == aln.to_dict()
    assert got is not aln


@pytest.mark.parametrize(
    "mk_cls", [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs]
)
def test_alignment_to_rich_dict_round_trip_rc(mk_cls):
    data = {
        "seq1": "ATG",
        "seq2": "TGA",
    }
    aln = mk_cls(data, moltype="dna").rc()
    rd = aln.to_rich_dict()
    got = deserialise_object(rd)
    assert got.to_dict() == aln.to_dict()
    assert got._is_reversed == False


@pytest.mark.parametrize(
    "mk_cls", [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs]
)
def test_alignment_to_rich_dict_round_trip_renamed(mk_cls):
    # name_map is preserved
    data = {
        "seq1": "ATG",
        "seq2": "TGA",
    }
    aln = mk_cls(data, moltype="dna")
    renamed_aln = aln.rename_seqs(renamer=lambda x: x.upper())
    rd = renamed_aln.to_rich_dict()
    got = deserialise_object(rd)
    assert got._name_map == {"SEQ1": "seq1", "SEQ2": "seq2"}
    assert got.to_dict() == renamed_aln.to_dict()
    assert got is not renamed_aln


@pytest.mark.parametrize(
    "mk_cls", [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs]
)
def test_alignment_to_rich_dict_round_trip_offset(mk_cls):
    # offset is not preserved
    data = {
        "seq1": "ATG",
        "seq2": "TGA",
    }
    offsets = {"seq1": 10, "seq2": 20}
    aln = mk_cls(data, moltype="dna", offset=offsets)
    rd = aln.to_rich_dict()
    got = deserialise_object(rd)._seqs_data.offset
    expect = {"seq1": 0, "seq2": 0}
    assert got == expect


@pytest.mark.parametrize(
    "mk_cls", [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs]
)
def test_alignment_to_rich_dict_round_trip_info(mk_cls):
    # info is preserved
    data = {
        "seq1": "ATG",
        "seq2": "TGA",
    }
    aln = mk_cls(data, moltype="dna", info={"key": "value"})
    rd = aln.to_rich_dict()
    got = deserialise_object(rd)
    assert got.info["key"] == "value"
    assert got.to_dict() == aln.to_dict()
    assert got is not aln


@pytest.mark.parametrize(
    "mk_cls", [new_alignment.make_unaligned_seqs, new_alignment.make_aligned_seqs]
)
def test_alignment_to_rich_dict_round_trip_annotation_db(gff_db, mk_cls):
    # serialisation will drop the annotation_db
    data = {
        "seq1": "ATG",
        "seq2": "TGA",
    }
    aln = mk_cls(data, moltype="dna")
    aln.annotation_db = gff_db
    assert len(aln.annotation_db) > 0
    rd = aln.to_rich_dict()
    got = deserialise_object(rd)
    assert len(got.annotation_db) == 0


def test_deserialise_alignment():
    data = {
        "init_args": {
            "moltype": "dna",
            "name_map": {"new_seq1": "seq1", "new_seq2": "seq2"},
            "info": {},
        },
        "type": "cogent3.core.new_alignment.Alignment",
        "version": "2023.10",
        "seqs": {"seq1": "ATCG", "seq2": "TAGC"},
    }

    aln = deserialise_object(data)

    assert isinstance(aln, new_alignment.Alignment)
    assert aln.names == ["new_seq1", "new_seq2"]
    assert str(aln.get_seq("new_seq1")) == "ATCG"
    assert str(aln.get_seq("new_seq2")) == "TAGC"
