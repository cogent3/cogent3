import os

import numpy
import pytest
from numpy import log2
from numpy.testing import assert_allclose

from cogent3 import (
    get_app,
    get_moltype,
    load_aligned_seqs,
    make_aligned_seqs,
    make_tree,
    make_unaligned_seqs,
)
from cogent3.align.align import (
    local_pairwise,
    make_dna_scoring_dict,
    make_generic_scoring_dict,
)
from cogent3.app import align as align_app
from cogent3.app.align import (
    _combined_refseq_gaps,
    _gap_difference,
    _gap_union,
    _GapOffset,
    _gaps_for_injection,
    _merged_gaps,
    pairwise_to_multiple,
    smith_waterman,
)
from cogent3.app.composable import NotCompleted
from cogent3.core.location import gap_coords_to_map

DNA = get_moltype("dna")
_NEW_TYPE = "COGENT3_NEW_TYPE" in os.environ

if _NEW_TYPE:
    from cogent3.core.new_alignment import Aligned
else:
    from cogent3.core.alignment import Aligned

_seqs = {
    "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
    "Bandicoot": "NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
    "Rhesus": "GCCAGCTCATTACAGCATGAGAACAGTTTGTTACTCACT",
    "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAGTTTATTATACACT",
}

_nucleotide_models = [
    "JC69",
    "K80",
    "F81",
    "HKY85",
    "TN93",
    "GTR",
    "ssGN",
    "GN",
]

_codon_models = [
    "CNFGTR",
    "CNFHKY",
    "MG94HKY",
    "MG94GTR",
    "GY94",
    "H04G",
    "H04GK",
    "H04GGK",
    "GNC",
]


def make_pairwise(data, refseq_name, moltype="dna", array_align=False):
    """returns series of refseq, [(n, pwise aln),..]. All alignments are to ref_seq"""
    aln = make_aligned_seqs(
        data,
        array_align=array_align,
        moltype=moltype,
    )
    refseq = aln.get_seq(refseq_name)
    pwise = [
        (n, aln.take_seqs([refseq_name, n]).omit_gap_pos())
        for n in aln.names
        if n != refseq_name
    ]
    return refseq, pwise


def make_aligned(gaps_lengths, seq, name="seq1"):
    seq = seq.moltype.make_seq(seq=seq, name=name)
    seq.name = name
    return Aligned.from_map_and_seq(gap_coords_to_map(gaps_lengths, len(seq)), seq)


@pytest.fixture
def refalignment_seqs():
    """Fixture providing seqs for refalignment tests"""
    return make_unaligned_seqs(_seqs, moltype=DNA)


def test_align_to_ref(refalignment_seqs):
    """correctly aligns to a reference"""
    aligner = align_app.align_to_ref(ref_seq="Human")
    aln = aligner(refalignment_seqs)
    expect = {
        "Bandicoot": "---NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
        "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAG---TTTATTATACACT",
        "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
        "Rhesus": "GCCAGCTCATTACAGCATGAGAAC---AGTTTGTTACTCACT",
    }
    assert aln.to_dict() == expect


@pytest.mark.parametrize(
    "test_moltype",
    ["text", "rna", "protein", "protein_with_stop"],
)
def test_align_to_ref_generic_moltype(test_moltype):
    """tests when the moltype is generic"""
    aligner = align_app.align_to_ref(moltype=test_moltype)
    assert aligner._moltype.label == test_moltype
    assert aligner._kwargs["S"] == make_generic_scoring_dict(
        10,
        get_moltype(test_moltype),
    )


def test_align_to_ref_result_has_moltype(refalignment_seqs):
    """aligned object has correct moltype"""
    aligner = align_app.align_to_ref(moltype="dna")
    got = aligner(refalignment_seqs)
    assert got.moltype.label == "dna"


def test_merged_gaps(refalignment_seqs):
    """correctly merges gaps"""
    a = {2: 3, 4: 9}
    b = {2: 6, 8: 5}
    # omitting one just returns the other
    assert _merged_gaps(a, {}) is a
    assert _merged_gaps({}, b) is b
    got = _merged_gaps(a, b)
    assert got == [(2, 6), (4, 9), (8, 5)]


def test_aln_to_ref_known(refalignment_seqs):
    """correctly recapitulates known case"""
    orig = make_aligned_seqs(
        {
            "Ref": "CAG---GAGAACAGAAACCCAT--TACTCACT",
            "Qu1": "CAG---GAGAACAG---CCCGTGTTACTCACT",
            "Qu2": "CAGCATGAGAACAGAAACCCGT--TA---ACT",
            "Qu3": "CAGCATGAGAACAGAAACCCGT----CTCACT",
            "Qu4": "CAGCATGAGAACAGAAACCCGTGTTACTCACT",
            "Qu5": "CAG---GAGAACAG---CCCAT--TACTCACT",
            "Qu6": "CAG---GA-AACAG---CCCAT--TACTCACT",
            "Qu7": "CAG---GA--ACAGA--CCCGT--TA---ACT",
        },
        moltype="dna",
    )
    expect = orig.to_dict()
    aligner = align_app.align_to_ref(ref_seq="Ref")
    aln = aligner.main(orig.degap())
    assert aln.to_dict() == expect


def test_gap_union(refalignment_seqs):
    """correctly identifies the union of all gaps"""
    # fails if not all sequences same
    seq = DNA.make_seq(seq="AACCCGTT")
    all_gaps = {0: 3, 2: 1, 5: 3, 6: 3}
    make_aligned(all_gaps, seq)
    gap_sets = [
        {5: 1, 6: 3},
        {2: 1, 5: 3},
        {2: 1, 5: 1, 6: 2},
        {0: 3},
    ]
    seqs = [make_aligned(gaps, seq) for gaps in gap_sets]
    got = _gap_union(seqs)
    assert got == dict(all_gaps)

    # must all be Aligned instances
    with pytest.raises(TypeError):
        _gap_union([*seqs, "GGGGGGGG"])

    # must all have the same name
    with pytest.raises(ValueError):
        _gap_union([*seqs, make_aligned({}, seq, name="blah")])


def test_gap_difference(refalignment_seqs):
    """correctly identifies the difference in gaps"""
    seq = DNA.make_seq(seq="AACCCGTT")
    {0: 3, 2: 1, 5: 3, 6: 3}
    gap_sets = [
        {5: 1, 6: 3},
        {2: 1, 5: 3},
        {2: 1, 5: 1, 6: 2},
        {0: 3},
    ]
    seqs = [make_aligned(gaps, seq) for gaps in gap_sets]
    union = _gap_union(seqs)
    expects = [
        [{0: 3, 2: 1}, {5: 2}],
        [{0: 3, 6: 3}, {}],
        [{0: 3}, {5: 2, 6: 1}],
        [{2: 1, 5: 3, 6: 3}, {}],
    ]
    for seq, (plain, overlap) in zip(seqs, expects, strict=False):
        seq_gaps = dict(seq.map.get_gap_coordinates())
        got_plain, got_overlap = _gap_difference(seq_gaps, union)
        assert got_plain == dict(plain)
        assert got_overlap == dict(overlap)


def test_merged_gaps(refalignment_seqs):
    """correctly handles gap values"""
    a_gaps = {0: 2}
    b_gaps = {2: 2}
    assert _merged_gaps(a_gaps, {}) == a_gaps
    assert _merged_gaps({}, b_gaps) == b_gaps


def test_combined_refseq_gaps(refalignment_seqs):
    union = {0: 3, 2: 1, 5: 3, 6: 3}
    gap_sets = [
        [(5, 1), (6, 3)],
        [(2, 1), (5, 3)],
        [(2, 1), (5, 1), (6, 2)],
        [(0, 3)],
    ]
    # for subset gaps, their alignment position is the
    # offset + their position + their gap length
    expects = [
        {6: 2, 0: 3, 2: 1},
        {0: 3, 10: 3},
        {0: 3, 5 + 1 + 1: 2, 6 + 2 + 2: 1},
        {2 + 3: 1, 5 + 3: 3, 6 + 3: 3},
    ]
    for i, gap_set in enumerate(gap_sets):
        got = _combined_refseq_gaps(dict(gap_set), union)
        assert got == expects[i]

    # if union gaps equals ref gaps
    got = _combined_refseq_gaps({2: 2}, {2: 2})
    assert got == {}


def test_gaps_for_injection(refalignment_seqs):
    # for gaps before any otherseq gaps, alignment coord is otherseq coord
    oseq_gaps = {2: 1, 6: 2}
    rseq_gaps = {0: 3}
    expect = {0: 3, 2: 1, 6: 2}
    seqlen = 50
    got = _gaps_for_injection(oseq_gaps, rseq_gaps, seqlen)
    assert got == expect
    # for gaps after otherseq gaps seq coord is align coord minus gap
    # length totals
    got = _gaps_for_injection(oseq_gaps, {4: 3}, seqlen)
    expect = {2: 1, 3: 3, 6: 2}
    assert got == expect
    got = _gaps_for_injection(oseq_gaps, {11: 3}, seqlen)
    expect = {2: 1, 6: 2, 8: 3}
    assert got == expect
    # gaps beyond sequence length added to end of sequence
    got = _gaps_for_injection({2: 1, 6: 2}, {11: 3, 8: 3}, 7)
    expect = {2: 1, 6: 2, 7: 6}
    assert got == expect


def test_pairwise_to_multiple(refalignment_seqs):
    """the standalone function constructs a multiple alignment"""
    expect = {
        "Ref": "CAG---GAGAACAGAAACCCAT--TACTCACT",
        "Qu1": "CAG---GAGAACAG---CCCGTGTTACTCACT",
        "Qu2": "CAGCATGAGAACAGAAACCCGT--TA---ACT",
        "Qu3": "CAGCATGAGAACAGAAACCCGT----CTCACT",
        "Qu7": "CAG---GA--ACAGA--CCCGT--TA---ACT",
        "Qu4": "CAGCATGAGAACAGAAACCCGTGTTACTCACT",
        "Qu5": "CAG---GAGAACAG---CCCAT--TACTCACT",
        "Qu6": "CAG---GA-AACAG---CCCAT--TACTCACT",
    }
    aln = make_aligned_seqs(expect, moltype="dna").omit_gap_pos()
    expect = aln.to_dict()
    for refseq_name in ["Qu3"]:
        refseq, pwise = make_pairwise(expect, refseq_name)
        got = pairwise_to_multiple(pwise, ref_seq=refseq, moltype=refseq.moltype)
        assert len(got) == len(aln)
        orig = dict(pwise)
        _, pwise = make_pairwise(got.to_dict(), refseq_name)
        got = dict(pwise)
        # should be able to recover the original pairwise alignments
        for key, value in got.items():
            assert value.to_dict() == orig[key].to_dict(), refseq_name

        with pytest.raises(TypeError):
            pairwise_to_multiple(pwise, "ACGG", DNA)


def make_pwise_from_dict(data, ref_name):
    """Helper function to create pairwise alignments from a dictionary of sequences

    Parameters
    ----------
    data : dict
        Nested dictionary of sequences with structure {name: {ref_name: seq, name: seq}}
    ref_name : str
        Name of the reference sequence

    Returns
    -------
    list, cogent3.core.sequence.Sequence
        List of [name, alignment] pairs and reference sequence
    """
    result = []
    for n, seqs in data.items():
        result.append(
            [n, make_aligned_seqs(data=seqs, moltype="dna", array_align=False)],
        )
    ref_seq = result[0][1].get_seq(ref_name)
    return result, ref_seq


def test_pairwise_to_multiple_long_seqs():
    """correctly handle alignments with long sequences"""
    pwise = {
        "Platypus": {
            "Opossum": "-----------------GTGC------GAT-------------------------------CCAAAAACCTGTGTC--ACCGT--------GCC----CAGAGCCTCC----CTCAGGCCGCTCGGGGAG---TG-------GCCCCCCG--GC-GGAGGGCAGGGATGGGGAGT-AGGGGTGGCAGTC----GGAACTGGAAGAGCTT-TACAAACC---------GA--------------------GGCT-AGAGGGTC-TGCTTAC-------TTTTTACCTTGG------------GTTTG-CCAGGAGGTAG----------AGGATGA-----------------CTAC--ATCAAG----AGC------------TGGG-------------",
            "Platypus": "CAGGATGACTACATCAAGAGCTGGGAAGATAACCAGCAAGGAGATGAAGCTCTGGACACTACCAAAGACCCCTGCCAGAACGTGAAGTGCAGCCGACACAAGGTCTGCATCGCTCAGGGCTACCAGAGAGCCATGTGTATCAGCCGCAAGAAGCTGGAGCACAGGATCAAGCAGCCAGCCCTGAAACTCCATGGAAACAGAGAGAGCTTCTGCAAGCCTTGTCACATGACCCAGCTGGCCTCTGTCTGCGGCTCGGACGGACACACTTACAGCTCCGTGTGCAAACTGGAGCAGCAGGCCTGTCTGACCAGCAAGCAGCTGACAGTCAAGTGTGAAGGCCAGTGCCCGTGCCCCACCGATCATGTTCCAGCCTCCACCGCTGATGGAAAACAAGAGACCT",
        },
        "Wombat": {
            "Opossum": "GTGCGATCCAAAAACCTGTGTCACCGTGCCCAGAGCCTCCCTCAGGCCGCTCGG-GGAGTGGCCCCCCGGCGGAGGGCAGGGATGGGGAGTAGGGGTGGCAGTCGGAACTGGAAGAGCTTTACAAACCGAGGCTAGAGGGTCTGCTTACTTTTTACCTTGG------GTTT--GC-CAGGA---GGT----AGAGGATGACTACATCAAGAGCTGGG---------------------------",
            "Wombat": "--------CA----------TCACCGC-CCCTGCACC---------CGGCTCGGCGGAGGGGGATTCTAA-GGGGGTCAAGGATGGCGAG-ACCCCTGGCAATTTCA--TGGAGGA------CGAGCAATGGCT-----GTC-GTCCATCTCCCAGTATAGCGGCAAGATCAAGCACTGGAACCGCTTCCGAGACGATGACTACATCAAGAGCTGGGAGGACAGTCAGCAAGGAGATGAAGCGC",
        },
    }
    pwise, ref_seq = make_pwise_from_dict(pwise, "Opossum")
    aln = pairwise_to_multiple(pwise, ref_seq, ref_seq.moltype)
    assert not isinstance(aln, NotCompleted)


def test_pairwise_to_multiple_short_seqs():
    """correctly handle alignments with short sequences"""
    pwise = {
        "Platypus": {
            "Opossum": "-----------------GTGC------GAT-------------------------------CCAAAAACCTGTGTC",
            "Platypus": "CAGGATGACTACATCAAGAGCTGGGAAGATAACCAGCAAGGAGATGAAGCTCTGGACACTACCAAAGACCCCTGCC",
        },
        "Wombat": {
            "Opossum": "GTGCGATCCAAAAACCTGTGTC",
            "Wombat": "--------CA----------TC",
        },
    }
    pwise, ref_seq = make_pwise_from_dict(pwise, "Opossum")
    aln = pairwise_to_multiple(pwise, ref_seq, ref_seq.moltype)
    assert not isinstance(aln, NotCompleted)


@pytest.fixture
def progressive_seqs():
    """fixture providing seqs for progressive alignment tests"""
    return make_unaligned_seqs(_seqs, moltype=DNA)


@pytest.fixture
def progressive_treestring():
    """fixture providing tree string for progressive alignment tests"""
    return "(Bandicoot:0.4,FlyingFox:0.05,(Rhesus:0.06,Human:0.0):0.04);"


def test_progressive_align_protein_moltype():
    """tests guide_tree is None and moltype is protein"""
    from cogent3 import load_aligned_seqs

    seqs = load_aligned_seqs("data/nexus_aa.nxs", moltype="protein")
    seqs = seqs.degap()
    seqs = seqs.take_seqs(["Rat", "Cow", "Human", "Mouse", "Whale"])
    aligner = align_app.progressive_align(model="WG01")
    got = aligner(seqs)
    assert not isinstance(got, NotCompleted)
    aligner = align_app.progressive_align(model="protein")
    got = aligner(seqs)
    assert not isinstance(got, NotCompleted)


def test_progressive_align_nuc(progressive_seqs):
    """progressive alignment with nuc models"""
    aligner = align_app.progressive_align(model="TN93", distance="TN93")
    aln = aligner(progressive_seqs)
    # TODO: revert to isinstance when new_type is merged
    assert aln.__class__.__name__.endswith("Alignment")
    assert len(aln) == 42
    assert aln.moltype == aligner._moltype
    # TODO the following is not robust across operating systems
    # so commenting out for now, but needs to be checked
    # expect = {'Human': 'GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT',
    #           'Rhesus': 'GCCAGCTCATTACAGCATGAGAA---CAGTTTGTTACTCACT',
    #           'Bandicoot': 'NACTCATTAATGCTTGAAACCAG---CAGTTTATTGTCCAAC',
    #           'FlyingFox': 'GCCAGCTCTTTACAGCATGAGAA---CAGTTTATTATACACT'}
    # got = aln.to_dict()
    # assert got == expect


def test_progressive_fails():
    """should return NotCompletedResult along with message"""
    # Bandicoot has an inf-frame stop codon
    seqs = make_unaligned_seqs(
        data={"Human": "GCCTCA", "Rhesus": "GCCAGCTCA", "Bandicoot": "TGATCATTA"},
        moltype="dna",
    )
    aligner = align_app.progressive_align(model="codon")
    got = aligner(seqs)
    assert isinstance(got, NotCompleted)


@pytest.mark.parametrize("model", _nucleotide_models)
def test_progressive_align_nuc_model(progressive_seqs, model):
    """progressive alignment with all nuc models"""
    aligner = align_app.progressive_align(model=model)
    aln = aligner(progressive_seqs)
    assert not isinstance(aln, NotCompleted)


@pytest.mark.parametrize("model", _codon_models)
def test_progressive_align_codon_model(progressive_seqs, model):
    """progressive alignment with all codon models"""
    aligner = align_app.progressive_align(model=model)
    aln = aligner(progressive_seqs)
    assert not isinstance(aln, NotCompleted)


def test_progressive_align_guide_tree(progressive_seqs, progressive_treestring):
    """progressive alignment works with guide tree"""
    tree = make_tree(treestring=progressive_treestring)
    aligner = align_app.progressive_align(model="TN93", guide_tree=tree)
    aln = aligner(progressive_seqs)
    assert not isinstance(aln, NotCompleted)
    assert len(aln) == 42
    assert aln.moltype == aligner._moltype


def test_progressive_align_model_guide_tree(progressive_seqs, progressive_treestring):
    """progressive alignment works with model guide tree"""
    tree = make_tree(treestring=progressive_treestring)
    aligner = align_app.progressive_align(model="TN93", guide_tree=tree)
    aln = aligner(progressive_seqs)
    assert not isinstance(aln, NotCompleted)
    assert len(aln) == 42
    assert aln.moltype == aligner._moltype


def test_gap_offset_empty():
    """create an empty offset"""
    goff = _GapOffset({})
    for i in range(4):
        assert goff[i] == 0

    goff = _GapOffset({}, invert=True)
    for i in range(4):
        assert goff[i] == 0


def test_gap_offset_repr_str():
    """repr and str work"""
    goff = _GapOffset({}, invert=True)
    for func in (str, repr):
        assert func(goff) == "{}"


def test_gap_offset():
    goff = _GapOffset({1: 2, 3: 4})
    assert goff.min_pos == 1
    assert goff.max_pos == 3
    assert goff.total == 6
    assert goff[0] == 0
    assert goff[1] == 0
    assert goff[2] == 2
    assert goff[3] == 2
    assert goff[4] == 6


def test_gap_offset_invert():
    aln2seq = _GapOffset({2: 1, 5: 2, 7: 2}, invert=True)
    assert aln2seq._store == {3: 1, 2: 0, 8: 3, 6: 1, 12: 5, 10: 3}
    assert aln2seq.max_pos == 12
    assert aln2seq.min_pos == 2
    assert aln2seq[11] == 3
    seq2aln = _GapOffset({2: 1, 5: 2, 7: 2})
    for seq_pos in range(20):
        aln_pos = seq_pos + seq2aln[seq_pos]
        assert aln_pos - aln2seq[aln_pos] == seq_pos


@pytest.mark.parametrize("array_align", [True, False])
def test_information_content_score(array_align):
    """Tests that the alignment_quality generates the right alignment quality
    value based on the Hertz-Stormo metric. expected values are hand calculated
    using the formula in the paper."""
    app_equifreq = get_app("ic_score", equifreq_mprobs=True)
    app_not_equifreq = get_app("ic_score", equifreq_mprobs=False)

    aln = make_aligned_seqs(
        data=["AATTGA", "AGGTCC", "AGGATG", "AGGCGT"],
        moltype="dna",
        array_align=array_align,
    )
    got = app_equifreq(aln)
    expect = log2(4) + (3 / 2) * log2(3) + (1 / 2) * log2(2) + (1 / 2) * log2(2)
    assert_allclose(got, expect)
    # should be the same with the default moltype too
    aln = make_aligned_seqs(
        ["AATTGA", "AGGTCC", "AGGATG", "AGGCGT"],
        moltype="text",
        array_align=array_align,
    )
    got = app_equifreq(aln)
    assert_allclose(got, expect)

    aln = make_aligned_seqs(
        data=["AAAC", "ACGC", "AGCC", "A-TC"],
        moltype="dna",
        array_align=array_align,
    )
    got = app_not_equifreq(aln)
    expect = (
        2 * log2(1 / 0.4)
        + log2(1 / (4 * 0.4))
        + (1 / 2) * log2(1 / (8 / 15))
        + (1 / 4) * log2(1 / (4 / 15))
    )
    assert_allclose(got, expect)

    # 1. Alignment just gaps - alignment_quality returns 0.0
    aln = make_aligned_seqs(
        data=["----", "----"],
        moltype="dna",
        array_align=array_align,
    )
    got = app_equifreq(aln)
    assert_allclose(got, 0.0)

    # 2 Just one sequence - alignment_quality returns 0.0
    aln = make_aligned_seqs(data=["AAAC"], moltype="dna", array_align=array_align)
    got = app_equifreq(aln)
    assert_allclose(got, 0.0)

    # 3.1 Two seqs, one all gaps. (equifreq_mprobs=True)
    aln = make_aligned_seqs(
        data=["----", "ACAT"],
        moltype="dna",
        array_align=array_align,
    )
    got = app_equifreq(aln)
    assert_allclose(got, 1.1699250014423124)

    # 3.2 Two seqs, one all gaps. (equifreq_mprobs=False)
    aln = make_aligned_seqs(
        data=["----", "AAAA"],
        moltype="dna",
        array_align=array_align,
    )
    got = app_not_equifreq(aln)
    assert_allclose(got, -2)


@pytest.fixture
def aln():
    aligner = align_app.progressive_align(model="TN93", distance="TN93")
    seqs = make_unaligned_seqs(_seqs, moltype=DNA)
    return aligner(seqs)


@pytest.fixture
def seqs():
    return make_unaligned_seqs(_seqs, moltype=DNA)


def test_cogent3_score(aln):
    get_score = get_app("cogent3_score")
    score = get_score(aln)
    assert score < -100


@pytest.mark.parametrize("del_all_params", [True, False])
def test_cogent3_score_missing(aln, del_all_params):
    get_score = get_app("cogent3_score")
    if del_all_params:
        aln.info.pop("align_params")
    else:
        aln.info["align_params"].pop("lnL")
    score = get_score(aln)
    assert isinstance(score, NotCompleted)


def test_sp_score_exclude_gap():
    # no gap penalty
    app = get_app("sp_score", calc="pdist", gap_extend=0, gap_insert=0)
    data = {"s1": "AAGAA-A", "s2": "-ATAATG", "s3": "C-TGG-G"}
    # prop unchanged s1-s2, s1-s3
    expect = sum([6 * 3 / 6, 0, 5 * 2 / 5])
    aln = make_aligned_seqs(data, moltype="dna", array_align=False)
    got = app.main(aln)
    assert_allclose(got, expect)


def test_sp_fail():
    aln = make_aligned_seqs(
        data={"a": "ATG---------AATCGAAGA", "b": "GTG---------GAAAAGCAG"},
        moltype="dna",
    )
    app = get_app("sp_score")
    got = app.main(aln)
    assert isinstance(got, NotCompleted)
    assert "NaN" in got.message


def test_sp_score_additive_gap():
    # additive gap score
    app = get_app("sp_score", calc="pdist", gap_extend=1, gap_insert=0)
    data = {"s1": "AAGAA-A", "s2": "-ATAATG", "s3": "C-TGG-G"}
    # match score
    mscore = numpy.array([6 * 3 / 6, 0, 5 * 2 / 5])
    # gap score
    gscore = numpy.array([2, 1, 3])
    aln = make_aligned_seqs(data, moltype="dna")
    got = app.main(aln)
    assert_allclose(got, (mscore - gscore).sum())


def test_sp_score_affine_gap():
    # affine gap score
    app = get_app("sp_score", calc="pdist", gap_extend=1, gap_insert=2)
    data = {"a": "AAGAA-A", "b": "-ATAATG", "c": "C-TGG-G"}
    # match score
    mscore = numpy.array([6 * 3 / 6, 0, 5 * 2 / 5])
    # gap score
    gscore = numpy.array([2 + 4, 2 + 1, 3 + 6])
    aln = make_aligned_seqs(data, moltype="dna")
    got = app.main(aln)
    assert_allclose(got, (mscore - gscore).sum())


def test_progressive_align_one_seq(seqs):
    """progressive alignment with no provided tree and approx_dists=False
    will use a quick alignment to build the tree"""
    aligner = align_app.progressive_align(model="TN93", approx_dists=True)
    seqs = seqs.take_seqs(seqs.names[0])
    got = aligner(seqs)
    assert isinstance(got, NotCompleted)


def test_progressive_align_tree_from_reference(seqs):
    """progressive alignment with no provided tree and approx_dists=False
    will use a quick alignment to build the tree"""
    aligner = align_app.progressive_align(model="TN93", approx_dists=False)
    aln = aligner(seqs)
    # TODO: revert to isinstance when new_type is merged
    assert aln.__class__.__name__.endswith("Alignment")
    assert len(aln) == 42
    assert aln.moltype == aligner._moltype


def test_progressive_align_tree_from_approx_dist(seqs):
    """progressive alignment with no provided tree and approx_dists=True
    will use an approximated distance measure to build the tree"""
    aligner = align_app.progressive_align(model="TN93", approx_dists=True)
    aln = aligner(seqs)
    # TODO: revert to isinstance when new_type is merged
    assert aln.__class__.__name__.endswith("Alignment")
    assert len(aln) == 42
    assert aln.moltype == aligner._moltype


def test_progressive_align_iters(seqs):
    """progressive alignment works with iters>1"""
    aligner = align_app.progressive_align(model="TN93")
    aln = aligner(seqs)
    # TODO: revert to isinstance when new_type is merged
    assert aln.__class__.__name__.endswith("Alignment")
    assert len(aln) == 42
    assert aln.moltype == aligner._moltype


def test_smith_waterman_matches_local_pairwise(seqs):
    aligner = smith_waterman()
    coll = make_unaligned_seqs(
        data=[seqs.get_seq("Human"), seqs.get_seq("Bandicoot")],
        moltype="dna",
    )
    got = aligner(coll)
    s = make_dna_scoring_dict(10, -1, -8)
    insertion = 20
    extension = 2
    expect = local_pairwise(
        seqs.get_seq("Human"),
        seqs.get_seq("Bandicoot"),
        s,
        insertion,
        extension,
        return_score=False,
    )
    assert got.to_dict() == expect.to_dict()


def test_smith_waterman_score(seqs):
    aligner = smith_waterman()
    coll = make_unaligned_seqs(
        data=[seqs.get_seq("Human"), seqs.get_seq("Bandicoot")],
        moltype="dna",
    )
    aln = aligner(coll)
    got = aln.info["align_params"]["sw_score"]
    s = make_dna_scoring_dict(10, -1, -8)
    insertion = 20
    extension = 2
    _, expect = local_pairwise(
        seqs.get_seq("Human"),
        seqs.get_seq("Bandicoot"),
        s,
        insertion,
        extension,
        return_score=True,
    )
    assert got == expect


@pytest.mark.parametrize(
    "moltype",
    ["text", "rna", "protein", "protein_with_stop"],
)
def test_smith_waterman_generic_moltype(moltype):
    """tests when the moltype is generic"""
    aligner = smith_waterman(moltype=moltype)
    assert aligner._score_matrix == make_generic_scoring_dict(10, get_moltype(moltype))


def test_smith_waterman_no_moltype(seqs):
    """If no moltype is provided and the SequenceCollection has no specified moltype, the
    default moltype ('dna') should be used.
    """
    aligner = smith_waterman()
    coll = make_unaligned_seqs(
        data=[seqs.get_seq("Human"), seqs.get_seq("Bandicoot")],
        moltype="dna",
    )
    aln = aligner(coll)
    assert aln.moltype.label == "dna"


@pytest.mark.parametrize("moltype_1", ["text", "dna", "rna", "protein"])
@pytest.mark.parametrize("moltype_2", ["text", "dna", "rna", "protein"])
def test_smith_waterman_wrong_moltype(moltype_1, moltype_2):
    """If the moltypes differ between SW app and SequenceCollection,
    the SW moltype should be used
    """
    aligner = smith_waterman(moltype=moltype_1)
    coll = make_unaligned_seqs(
        data={"Human": "AUUCGAUGG", "Bandicoot": "AUUGCCCGAUGG"},
        moltype=moltype_2,
    )
    aln = aligner(coll)
    assert aln.moltype.label == moltype_1


def test_smith_waterman_raises(seqs):
    """SW should fail when given a SequenceCollection that deos not contain 2 seqs"""
    aligner = smith_waterman()
    coll = make_unaligned_seqs(
        data=[seqs.get_seq("Human"), seqs.get_seq("Bandicoot"), seqs.get_seq("Rhesus")],
        moltype="dna",
    )
    aln = aligner(coll)
    assert isinstance(aln, NotCompleted)

    coll = make_unaligned_seqs(data=[seqs.get_seq("Human")], moltype="dna")
    aln = aligner(coll)
    assert isinstance(aln, NotCompleted)


def test_aln_two():
    """correctly recapitulates known case"""
    orig = make_aligned_seqs(
        {
            "Ref": "CAGGAGAACAGAAACCCATTACTCACT",
            "Qu7": "CAGGA--ACAGA--CCCGTTA---ACT",
        },
        moltype="dna",
    )
    expect = orig.to_dict()
    aligner = align_app.align_to_ref(ref_seq="Ref")
    seqs = orig.degap()
    aln = aligner.main(seqs)
    assert aln.to_dict() == expect


def test_codon_incomplete(DATA_DIR):
    names = ["FlyingFox", "DogFaced", "FreeTaile"]
    aln = load_aligned_seqs(DATA_DIR / "brca1.fasta", moltype="dna")
    seqs = aln.take_seqs(names)[2700:3000].degap()
    aligner = align_app.progressive_align("codon")
    aln = aligner(seqs)
    assert aln  # will fail if aln is a NotCompleted instance
    # now make sure the resulting ungapped sequences are modulo 3
    seqs = aln.degap().to_dict().values()
    assert {len(s) % 3 for s in seqs} == {0}
