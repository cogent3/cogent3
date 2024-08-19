import json
import os
import pathlib
import re
from warnings import catch_warnings, filterwarnings

import numpy
import pytest

from cogent3 import load_unaligned_seqs, open_
from cogent3._version import __version__
from cogent3.core import new_alignment, new_alphabet, new_moltype, new_sequence
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import GffAnnotationDb, load_annotations
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
def alpha(moltype="dna"):
    moltype = new_moltype.get_moltype(moltype)
    return moltype.degen_gapped_alphabet


@pytest.fixture
def dna_sd(str_seqs_dict: dict[str, str], alpha):
    return new_alignment.SeqsData(data=str_seqs_dict, alphabet=alpha)


@pytest.fixture
def int_arr():
    return numpy.arange(17, dtype=numpy.uint8)


@pytest.fixture
def sdv_s2(dna_sd: new_alignment.SeqsData) -> new_alignment.SeqDataView:
    return dna_sd.get_seq_view(seqid="seq2")


@pytest.fixture(scope="function")
def seqs() -> new_alignment.SequenceCollection:
    data = {"seq1": "AAAAAA", "seq2": "TTTT", "seq3": "ATTCCCC"}
    return new_alignment.make_unaligned_seqs(data, moltype="dna")


@pytest.fixture
def ragged_padded():
    return new_alignment.make_unaligned_seqs(
        {"a": "AAAAAA", "b": "AAA---", "c": "AAAA--"}, moltype="dna"
    )


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


def test_seqs_data_default_attributes(dna_sd: new_alignment.SeqsData):
    assert dna_sd.names == ["seq1", "seq2", "seq3"]
    assert isinstance(dna_sd.alphabet, new_alphabet.CharAlphabet)


def test_seqs_data_view_zero_step_raises(dna_sd):
    with pytest.raises(ValueError):
        new_alignment.SeqDataView(seq=dna_sd, seqid=seq1, step=0, seq_len=4)


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
def test_seqs_data_view_repr_default(dna_sd: new_alignment.SeqsData, seqid: str):
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
    sd = new_alignment.SeqsData(data=d, alphabet=alpha)
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

    sd = new_alignment.SeqsData(data=data, alphabet=alpha)
    sdv = sd.get_seq_view(seqid=seqid)
    sliced_sdv = sdv[start:stop:step]

    copied_sdv = sliced_sdv.copy(sliced=sliced)

    assert copied_sdv.str_value == data[seqid][start:stop:step]


@pytest.mark.parametrize("seqid", ["seq1", "seq2"])
def test_seqs_data_get_seq_view(str_seqs_dict, alpha, seqid):
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=alpha)
    seq = str_seqs_dict[seqid]
    seq_len = len(seq)
    got = sd.get_seq_view(seqid)
    assert isinstance(got, new_alignment.SeqDataView)
    assert got.seq == sd
    assert got.stop == seq_len
    assert got.seqid == seqid
    assert got.seq_len == seq_len


@pytest.mark.parametrize("seq", ("seq1", "seq2"))
@pytest.mark.parametrize("start", (None, -1, 0, 1, 4))
@pytest.mark.parametrize("stop", (None, -1, 0, 1, 4))
def test_seqs_data_get_seq_str(str_seqs_dict, alpha, seq, start, stop):
    # slicing should be tested in test_get_seq_array
    expect = str_seqs_dict[seq][start:stop]
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=alpha)
    got = sd.get_seq_str(seqid=seq, start=start, stop=stop)
    assert expect == got


def test_seqs_data_get_seq_str_empty(dna_sd: new_alignment.SeqsData):
    with pytest.raises(TypeError):
        dna_sd.get_seq_str()


def test_seqs_data_names(str_seqs_dict, alpha):
    expect = str_seqs_dict.keys()
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=alpha)
    got = sd.names
    # returns iterator
    assert list(got) == list(expect)


def test_seqs_data_seq_lengths(str_seqs_dict, arr_seqs_dict, alpha):
    expect = {k: len(v) for k, v in str_seqs_dict.items()}
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=alpha)
    got = sd.seq_lengths()
    assert got == expect

    expect = {k: len(v) for k, v in arr_seqs_dict.items()}
    sd = new_alignment.SeqsData(data=arr_seqs_dict, alphabet=alpha)
    got = sd.seq_lengths()
    assert got == expect


def test_seqs_data_get_seq_array(str_seqs_dict, alpha):
    # seq1
    expect = numpy.array([2, 1, 3, 0], dtype="uint8")
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=alpha)
    got = sd.get_seq_array(seqid="seq1")
    assert numpy.array_equal(got, expect)

    # seq2
    expect = numpy.array([3, 0, 0, 0, 3, 1, 2], dtype="uint8")
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=alpha)
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
    assert got.seq == dna_sd
    assert got.seqid == seq


@pytest.mark.parametrize("idx", (0, 1))
def test_seqs_data_getitem_int(str_seqs_dict, alpha, idx):
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=alpha)
    got = sd[idx]
    assert got.seq == sd
    assert got.seqid == list(str_seqs_dict)[idx]


@pytest.mark.parametrize("seqid", ("seq1", "seq2"))
def test_seqs_data_getitem_str(
    dna_sd,
    seqid,
):
    got = dna_sd[seqid]
    assert got.seq == dna_sd
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
def test_seqs_data_round_trip(reverse, alpha):
    seqs_data = new_alignment.SeqsData(
        data={"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}, alphabet=alpha
    )
    seqs_data = seqs_data.reverse_seqs() if reverse else seqs_data

    rd = seqs_data.to_rich_dict()
    got = deserialise_object(rd)
    assert isinstance(got, new_alignment.SeqsData)
    assert got.to_rich_dict() == seqs_data.to_rich_dict()
    expect = "ACGG"[::-1] if reverse else "ACGG"
    assert str(got.get_seq_view("seq1")) == expect


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
    sdv = new_alignment.SeqDataView(seq=seq1, seqid="seq1", seq_len=len(seq1))
    got = sdv[index]
    assert isinstance(got, new_alignment.SeqDataView)


# SeqDataView tests for value properties
@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seq_data_view_value(str_seqs_dict: dict, alpha, start, stop, step):
    seq = "seq2"
    expect = str_seqs_dict[seq][start:stop:step]
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=alpha)
    # Get SeqDataView on seq
    sdv = sd.get_seq_view(seqid=seq)
    sdv2 = sdv[start:stop:step]
    got = sdv2.str_value
    assert got == expect


@pytest.mark.parametrize("rev", (False, True))
def test_seq_data_view_to_rich_dict(rev):
    data = {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}
    alpha = new_moltype.DNA.degen_gapped_alphabet
    sd = new_alignment.SeqsData(data=data, alphabet=alpha)
    sdv = sd.get_seq_view(seqid="seq1")
    sdv = sdv[::-1] if rev else sdv
    got = sdv.to_rich_dict()
    expect = {
        "init_args": {
            "seqid": "seq1",
            "step": sdv.step,
            "offset": sdv.offset,
            "seq_len": sdv.seq_len,
            "seq": sdv.str_value,
        },
        "type": get_object_provenance(sdv),
        "version": __version__,
    }

    assert got == expect


@pytest.mark.parametrize("start", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("stop", (None, 0, 1, 4, -1, -4))
@pytest.mark.parametrize("step", (None, 1, 2, 3, -1, -2, -3))
def test_seqs_data_array_value(arr_seqs_dict: dict, alpha, start, stop, step):
    seq = "seq2"
    expect = arr_seqs_dict[seq][start:stop:step]
    sd = new_alignment.SeqsData(data=arr_seqs_dict, alphabet=alpha)
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
    sd = new_alignment.SeqsData(data=str_seqs_dict, alphabet=alpha)
    # Get SeqDataView on seq
    sdv = sd.get_seq_view(seqid=seq)
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
    data = {"a": numpy.array(["AGGCCC"]), "b": numpy.array(["AGAAAA"])}
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
    data = [numpy.array(["AGGCCC"]), numpy.array(["AGAAAA"])]
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


def test_sequence_collection_init(seqs):
    assert isinstance(seqs.seqs, new_alignment.SeqsData)
    assert seqs.seqs.make_seq is not None


def test_sequence_collection_names_is_list():
    """expected to be a list"""
    seqs = new_alignment.make_unaligned_seqs(
        {"a": b"AAAAA", "b": b"TTTTT"}, moltype="dna"
    )
    assert isinstance(seqs.names, list)


def test_sequence_collection_names(seqs):
    assert seqs.names == ["seq1", "seq2", "seq3"]
    seqs.names = ["seq2", "seq3", "seq1"]
    assert seqs.names == ["seq2", "seq3", "seq1"]
    seqs.names = ["seq1", "seq2"]
    assert seqs.names == ["seq1", "seq2"]
    with pytest.raises(ValueError):
        seqs.names = ["seq1", "seq2", "seq3", "seq4"]


def test_make_unaligned_seqs_list_str():
    """SequenceCollection init from list of sequences should use indices as keys"""
    seqs = ["TTTTT", "CCCCC", "GGGGG"]
    a = new_alignment.make_unaligned_seqs(seqs, moltype="dna")
    assert len(a.seqs) == 3
    assert a.seqs["seq_0"] == "TTTTT"
    assert a.seqs["seq_1"] == "CCCCC"
    assert a.seqs["seq_2"] == "GGGGG"
    assert a.names == ["seq_0", "seq_1", "seq_2"]


def test_make_unaligned_seqs_pairs():
    """SequenceCollection init from list of (key,val) pairs should work correctly"""
    seqs = [["a", "AAA"], ["t", "TTT"], ["c", "CCC"]]
    a = new_alignment.make_unaligned_seqs(seqs, moltype="dna")
    assert len(a.seqs) == 3
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


def test_sequence_collection_init_ambig():
    """SequenceCollection should tolerate ambiguous chars"""
    _ = new_alignment.make_unaligned_seqs(["AAA", "CCC"], moltype="dna")
    _ = new_alignment.make_unaligned_seqs(["ANS", "CWC"], moltype="dna")
    _ = new_alignment.make_unaligned_seqs(["A-A", "CC-"], moltype="dna")
    _ = new_alignment.make_unaligned_seqs(["A?A", "CC-"], moltype="dna")


def test_sequence_collection_iter_seqs_ragged_padded(ragged_padded):
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
    a = ragged_padded.take_seqs(["b", "c"])
    assert isinstance(a, new_alignment.SequenceCollection)
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
    assert isinstance(a, new_alignment.SequenceCollection)
    assert a.names == ["b", "c"]
    assert a.num_seqs == 2


def test_sequence_collection_take_seqs_info():
    """take_seqs should preserve info attribute"""
    data = {"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"}
    orig = new_alignment.make_unaligned_seqs(
        data,
        moltype="dna",
        info={"key": "value"},
    )
    subset = orig.take_seqs(list("ab"))
    assert set(subset.info) == set(orig.info)


def test_sequence_collection_take_seqs_moltype():
    """take_seqs should preserve the MolType"""
    data = {"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"}
    orig = new_alignment.make_unaligned_seqs(data, moltype="dna")
    subset = orig.take_seqs(list("ab"))
    assert set(subset.moltype) == set(orig.moltype)


def test_sequence_collection_take_seqs_empty_names():
    """take_seqs should raise ValueError if no seqs are selected."""
    data = {"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"}
    orig = new_alignment.make_unaligned_seqs(data, moltype="dna")
    with pytest.raises(ValueError):
        _ = orig.take_seqs([])


def test_sequence_collection_take_seqs_copy_annotations(gff_db):
    data = {"test_seq": "ACGT", "test_seq2": "CGTTTA"}
    seq_coll = new_alignment.make_unaligned_seqs(
        data, moltype="dna", annotation_db=gff_db
    )
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
    assert [ragged.seqs.get_seq_str(seqid=name) for name in ragged.names] == [
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


def test_sequence_collection_init_annotated_seqs():
    """correctly construct from list with annotated seq"""
    seq = new_moltype.DNA.make_seq(seq="GCCAGGGGGGAAAG-GGAGAA", name="seq1")
    seq.add_feature(biotype="exon", name="name", spans=[(4, 10)])
    coll = new_alignment.make_unaligned_seqs([seq], moltype="dna")
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

    # pseudocount and allow gap
    data = {"a": "AC-", "b": "AC", "c": "AC"}
    aln = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = aln.get_motif_probs(pseudocount=1)
    expect = {"A": 4 / 11, "C": 4 / 11, "G": 1 / 11, "T": 1 / 11, "-": 1 / 11}


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


def test_sequence_collection_probs_per_seq():
    data = {"seq1": "AA??", "seq2": "CG-N", "seq3": "CGAA"}
    coll = new_alignment.make_unaligned_seqs(data, moltype="dna")
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


def test_sequence_collection_get_annotations_from_any_seq():
    """get_annotations_from_any_seq returns correct annotations"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
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


def test_sequence_collection_add_seqs():
    data = dict(
        [("name1", "AAA"), ("name2", "AAA"), ("name3", "AAA"), ("name4", "AAA")]
    )
    data2 = dict([("name5", "BBB"), ("name6", "CCC")])
    aln = new_alignment.make_unaligned_seqs(data, moltype="text", info={"key": "foo"})
    out_aln = aln.add_seqs(data2)
    assert len(out_aln.names) == 6


def test_sequence_collection_add_seqs_duplicate_raises():
    """add_seqs should raise an error if duplicate names"""
    data = dict(
        [("name1", "AAA"), ("name2", "AAA"), ("name3", "AAA"), ("name4", "AAA")]
    )
    data2 = dict([("name5", "BBB"), ("name1", "CCC")])
    aln = new_alignment.make_unaligned_seqs(data, moltype="text", info={"key": "foo"})
    with pytest.raises(ValueError):
        _ = aln.add_seqs(data2)

    # but it should work if we allow duplicates
    out_aln = aln.add_seqs(data2, force_unique_keys=False)
    assert len(out_aln.names) == 5


def test_sequence_collection_add_seqs_info():
    """add_seqs should preserve info attribute"""
    data = dict(
        [("name1", "AAA"), ("name2", "AAA"), ("name3", "AAA"), ("name4", "AAA")]
    )
    data2 = dict([("name5", "BBB"), ("name6", "CCC")])
    aln = new_alignment.make_unaligned_seqs(data, moltype="text", info={"key": "foo"})
    aln2 = new_alignment.make_unaligned_seqs(data2, moltype="text", info={"key": "bar"})
    out_aln = aln.add_seqs(aln2)
    assert out_aln.info["key"] == "foo"


def test_sequence_collection_write(tmp_path):
    """SequenceCollection.write should write in correct format"""
    data = {"a": "AAAA", "b": "TTTT", "c": "CCCC"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    fn = tmp_path / "seqs.fasta"
    seqs.write(fn)
    with open(fn, newline=None) as infile:
        result = infile.read()
    assert result == ">a\nAAAA\n>b\nTTTT\n>c\nCCCC\n"


def test_sequence_collection_write_gapped(tmp_path):
    data = {"a": "AAA--", "b": "TTTT", "c": "CCCC"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    fn = tmp_path / "seqs.fasta"
    seqs.write(fn)
    with open(fn, newline=None) as infile:
        result = infile.read()
    assert result == ">a\nAAA--\n>b\nTTTT\n>c\nCCCC\n"


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
    max_seq_len = max(seqs.seqs.seq_lengths().values())

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
            "data": {name: seqs.seqs.get_seq_str(seqid=name) for name in seqs.names},
            "alphabet": seqs.moltype.most_degen_alphabet().to_rich_dict(),
            "reversed_seqs": seqs.seqs.reversed,
        },
        "type": get_object_provenance(seqs.seqs),
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
            "data": {name: seqs.seqs.get_seq_str(seqid=name) for name in seqs.names},
            "alphabet": seqs.moltype.most_degen_alphabet().to_rich_dict(),
            "reversed_seqs": seqs.seqs.reversed,
        },
        "type": get_object_provenance(seqs.seqs),
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
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    reversed_seqs = seqs.reverse_complement()

    got = reversed_seqs.to_rich_dict()
    seqs_data = {
        "init_args": {
            "data": {name: seqs.seqs.get_seq_str(seqid=name) for name in seqs.names},
            "alphabet": seqs.moltype.most_degen_alphabet().to_rich_dict(),
            "reversed_seqs": {"seq1": True, "seq2": True, "seq3": True},
        },
        "type": get_object_provenance(seqs.seqs),
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


def test_sequence_collection_reverse_complement():
    data = {"s1": "AACC", "s2": "GGTT"}
    seqs = new_alignment.make_unaligned_seqs(data, moltype="dna")
    got = seqs.reverse_complement()
    expect = {"s1": "GGTT", "s2": "AACC"}
    assert got.to_dict() == expect
    assert str(got.get_seq(seqname="s1")) == "GGTT"

    got = got.reverse_complement()
    assert got.to_dict() == data
