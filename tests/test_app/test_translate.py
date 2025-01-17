import os

import pytest

import cogent3
from cogent3.app.composable import NotCompleted
from cogent3.app.translate import (
    best_frame,
    get_fourfold_degenerate_sets,
    select_translatable,
    translate_frames,
    translate_seqs,
)

DNA = cogent3.get_moltype("dna")

_NEW_TYPE = "COGENT3_NEW_TYPE" in os.environ


def test_best_frame():
    """correctly identify best frame with/without allowing rc"""
    seq = DNA.make_seq(seq="ATGCTAACATAAA", name="fake1")
    assert best_frame(seq) == 1
    assert best_frame(seq, require_stop=True) == 1

    # a challenging seq, translatable in 1 and 3 frames, ending on stop in
    # frame 1. Should return frame 1 irrespective of require_stop
    seq = DNA.make_seq(seq="ATGTTACGGACGATGCTGAAGTCGAAGATCCACCGCGCCACGGTGACCTGCTGA")
    assert best_frame(seq) == 1

    # a rc seq
    seq = DNA.make_seq(
        seq="AATATAAATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTTCATAAAGTCATA",
        name="fake2",
    )
    assert best_frame(seq, allow_rc=True) == 1
    with pytest.raises(ValueError):
        best_frame(seq, allow_rc=True, require_stop=True)

    rc = seq.rc()
    assert best_frame(rc, allow_rc=True) == -1


@pytest.mark.skipif(
    _NEW_TYPE,
    reason="new_type does not yet support mixed strand collections",
)
def test_select_translatable():
    """correctly get translatable seqs"""
    data = {
        "a": "AATATAAATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTTCATAAAGTCATA",
        "rc": "TATGACTTTATGAAGTAATAAACTGCTGTTCTCATGCTGTAATGAGCTGGCATTTATATT",
    }
    seqs = cogent3.make_unaligned_seqs(data=data, moltype=DNA)
    trans = select_translatable(allow_rc=False)
    tr = trans(seqs)  # pylint: disable=not-callable
    ex = data.copy()
    ex.pop("rc")
    assert tr.to_dict() == ex

    trans = select_translatable(allow_rc=True)
    tr = trans(seqs)  # pylint: disable=not-callable
    ex = data.copy()
    ex["rc"] = data["a"]
    assert tr.to_dict() == ex

    # if seqs not translatable returns NotCompletedResult
    data = {"a": "TAATTGATTAA", "b": "GCATAATTA"}
    seqs = cogent3.make_unaligned_seqs(data=data, moltype=DNA)
    got = select_translatable(allow_rc=False, frame=1)(seqs)  # pylint: disable=not-callable
    assert isinstance(got, NotCompleted)


def test_translate_frames():
    """returns translated sequences"""
    seq = DNA.make_seq(seq="ATGCTGACATAAA", name="fake1")
    tr = translate_frames(seq)
    assert tr == ["MLT*", "C*HK", "ADI"]
    # with the bacterial nuclear and plant plastid code
    tr = translate_frames(seq, gc="Euplotid Nuclear")
    assert tr == ["MLT*", "CCHK", "ADI"]


def test_translate_seqcoll():
    """correctly translate a sequence collection"""
    seqs = {"a": "ATGAGG", "b": "ATGTAA"}
    seqs = cogent3.make_unaligned_seqs(seqs, moltype="dna")
    # trim terminal stops
    translater = translate_seqs()
    aa = translater(seqs)  # pylint: disable=not-callable
    assert aa.to_dict() == {"a": "MR", "b": "M"}
    assert aa.moltype.label == "protein"
    # don't trim terminal stops, returns NotCompleted
    translater = translate_seqs(trim_terminal_stop=False)
    aa = translater(seqs)  # pylint: disable=not-callable
    assert isinstance(aa, NotCompleted)


def test_translate_aln():
    """correctly translates alignments"""
    data = {"a": "ATGAGGCCC", "b": "ATGTTT---"}
    # an array alignment
    aln = cogent3.make_aligned_seqs(data, moltype="dna")
    translater = translate_seqs()
    aa = translater(aln)  # pylint: disable=not-callable
    assert aa.to_dict() == {"a": "MRP", "b": "MF-"}
    assert aa.moltype.label == "protein"
    assert isinstance(aa, type(aln))
    # Alignment
    aln = aln.to_type(array_align=True)
    aa = translater(aln)  # pylint: disable=not-callable
    assert aa.to_dict() == {"a": "MRP", "b": "MF-"}
    assert aa.moltype.label == "protein"
    assert isinstance(aa, type(aln))


@pytest.mark.parametrize("code_id", range(1, 3))
def test_get_fourfold_degenerate_sets(code_id):
    """correctly identify 4-fold degenerate codons"""
    # using straight characters
    expect = set()
    for di in "GC", "GG", "CT", "CC", "TC", "CG", "AC", "GT":
        expect.update([frozenset(di + n for n in "ACGT")])

    got = get_fourfold_degenerate_sets(cogent3.get_code(code_id), as_indices=False)
    assert got == expect

    with pytest.raises(AssertionError):
        # as_indices requires an alphabet
        get_fourfold_degenerate_sets(cogent3.get_code(1), as_indices=True)

    expect = set()
    for di in "GC", "GG", "CT", "CC", "TC", "CG", "AC", "GT":
        codons = [tuple(DNA.alphabet.to_indices(x)) for x in [di + n for n in "ACGT"]]
        expect.update([frozenset(codons)])

    got = get_fourfold_degenerate_sets(
        cogent3.get_code(code_id),
        alphabet=DNA.alphabet,
        as_indices=True,
    )
    assert got == expect


@pytest.fixture(params=(None, 0, 1, 2))
def framed_seqs(DATA_DIR, request):
    # sample sequences with terminating stop codon
    # using valid values for frame
    data = {
        "NineBande": "GCAAGGCGCCAACAGAGCAGATGGGCTGAAAGTAAGGAAACATGTAATGATAGGCAGACTTAA",
        "Mouse": "GCAGTGAGCCAGCAGAGCAGATGGGCTGCAAGTAAAGGAACATGTAACGACAGGCAGGTTTAA",
        "Human": "GCAAGGAGCCAACATAACAGATGGGCTGGAAGTAAGGAAACATGTAATGATAGGCGGACTTAA",
        "HowlerMon": "GCAAGGAGCCAACATAACAGATGGGCTGAAAGTGAGGAAACATGTAATGATAGGCAGACTTAA",
        "DogFaced": "GCAAGGAGCCAGCAGAACAGATGGGTTGAAACTAAGGAAACATGTAATGATAGGCAGACTTAA",
    }
    prefix = "A" * (request.param or 0)
    frame = None if request.param is None else request.param + 1
    for k, s in data.items():
        data[k] = prefix + s
    return cogent3.make_unaligned_seqs(data=data, moltype="dna", info={"frame": frame})


def test_select_translatable_with_frame_terminal_stop(framed_seqs):
    frame = framed_seqs.info.frame
    sl = slice(None, None) if frame is None else slice(frame - 1, None)
    expect = {s.name: str(s[sl]) for s in framed_seqs.seqs}
    app = select_translatable(frame=frame, trim_terminal_stop=False)
    got = app(framed_seqs)  # pylint: disable=not-callable
    assert got.to_dict() == expect


def test_select_translatable_with_frame_no_stop(framed_seqs):
    frame = framed_seqs.info.frame
    sl = slice(None, -3) if frame is None else slice(frame - 1, -3)
    expect = {s.name: str(s[sl]) for s in framed_seqs.seqs}
    app = select_translatable(frame=frame, trim_terminal_stop=True)
    got = app(framed_seqs)  # pylint: disable=not-callable
    assert got.to_dict() == expect


def test_select_translatable_exclude_internal_stop():
    aln = cogent3.make_unaligned_seqs(
        {
            "internal_stop": "AATTAAATGTGA",
            "s2": "TATGACTAA",
        },
        moltype="dna",
    )
    app = select_translatable(frame=1)
    result = app(aln)  # pylint: disable=not-callable
    expect = {"s2": "TATGAC"}
    assert result.to_dict() == expect


@pytest.mark.parametrize("frame", [-1, 0, 4])
def test_select_translatable_invalid_frame(frame):
    with pytest.raises(AssertionError):
        _ = select_translatable(frame=frame)
