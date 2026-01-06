import pytest

import cogent3
from cogent3.core.location import FeatureMap, Span

DNA = cogent3.get_moltype("dna")


def close_dbs(*objs):
    for obj in objs:
        if not hasattr(obj, "has_annotation_db"):
            db = obj
        elif obj.has_annotation_db():
            db = obj.annotation_db
        else:
            continue

        db.close()


@pytest.mark.parametrize(
    ("coord", "expected"),
    [
        ((0, 4), "[0:4]"),
        ((0, 5), "[0:5]"),
        ((0, 6), "[0:5, 5:6]"),
        ((5, 10), "[5:10]"),
        ((-1, 10), "[-1-, 0:5, 5:10]"),
        ((5, 11), "[5:10, -1-]"),
        ((0, 10), "[0:5, 5:10]"),
        ((10, 0), "[10:5, 5:0]"),
    ],
)
def test_spans(coord, expected):
    # a simple two part map of length 10
    fmap = FeatureMap.from_locations(locations=[(0, 5), (5, 10)], parent_length=10)
    start, end = coord
    # try different spans on the above map
    r = repr(Span(start, end, reverse=start > end).remap_with(fmap))
    assert r == expected, repr((r, expected))


def test_get_by_annotation():
    seq = DNA.make_seq(seq="ATCGATCGAT" * 5, name="base")
    seq.add_feature(biotype="test_type", name="test_label", spans=[(5, 10)])
    seq.add_feature(biotype="test_type", name="test_label2", spans=[(15, 18)])

    answer = list(seq.get_features(biotype="test_type"))
    assert len(answer) == 2
    assert str(seq[answer[0]]) == "TCGAT"
    assert str(seq[answer[1]]) == "TCG"

    answer = list(seq.get_features(biotype="test_type", name="test_label"))
    assert len(answer) == 1
    assert str(seq[answer[0]]) == "TCGAT"

    # test ignoring of a partial annotation
    sliced_seq = seq[:17]
    answer = list(sliced_seq.get_features(biotype="test_type", allow_partial=False))
    assert len(answer) == 1
    assert str(sliced_seq[answer[0]]) == "TCGAT"
    close_dbs(seq, sliced_seq)


def test_get_by_seq_annotation_allow_gaps():
    aln = cogent3.make_aligned_seqs(
        {"a": "ATCGAAATCGAT", "b": "ATCGA--TCGAT"},
        moltype="dna",
    )
    # original version was putting annotation directly on seq
    f = aln.add_feature(
        seqid="b",
        biotype="test_type",
        name="test_label",
        spans=[(4, 6)],
        on_alignment=False,
    )
    # in the original get_by_seq_annotations(), inside the method it
    # used self[feature.map.start : feature.map.end],
    # so what was returned was an alignment slice compared to below,
    # which is only non-gap positions of "b"
    # when we set allow_gaps, we get gaps in an alignment
    f = next(iter(aln.get_features(seqid="b", biotype="test_type")))
    got = f.get_slice(allow_gaps=True).to_dict()
    assert got == {"b": "A--T", "a": "AAAT"}
    # otherwise we only get the exact positions represented
    got = f.get_slice(allow_gaps=False).to_dict()
    assert got == {"b": "AT", "a": "AT"}
    close_dbs(aln)
