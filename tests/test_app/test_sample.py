from unittest import TestCase

import pytest

import cogent3
from cogent3.app import composable, sample
from cogent3.app.composable import NotCompleted
from cogent3.core.moltype import MolTypeError

DNA = cogent3.get_moltype("dna")


class TranslateTests(TestCase):
    def _codon_positions(self):
        """correctly return codon positions"""
        aln = cogent3.make_aligned_seqs(
            [("a", "ACGACGACG"), ("b", "GATGATGAT")],
            moltype=DNA,
        )
        one = sample.take_codon_positions(1)
        got = one(aln)
        assert got.to_dict() == {"a": "AAA", "b": "GGG"}

        two = sample.take_codon_positions(2)
        got = two(aln)
        assert got.to_dict() == {"a": "CCC", "b": "AAA"}
        three = sample.take_codon_positions(3)
        got = three(aln)
        assert got.to_dict() == {"a": "GGG", "b": "TTT"}

        one_two = sample.take_codon_positions(1, 2)
        got = one_two(aln)
        assert got.to_dict() == {"a": "ACACAC", "b": "GAGAGA"}
        one_three = sample.take_codon_positions(1, 3)
        got = one_three(aln)
        assert got.to_dict() == {"a": "AGAGAG", "b": "GTGTGT"}
        two_three = sample.take_codon_positions(2, 3)
        got = two_three(aln)
        assert got.to_dict() == {"a": "CGCGCG", "b": "ATATAT"}

    def test_take_codon_positions_alignment(self):
        """correctly return codon positions from Alignment"""
        self._codon_positions()

    def test_take_named_3(self):
        """3 named seqs"""
        select = sample.take_named_seqs("a", "b", "c")
        assert select._init_vals == {"names": tuple("abc"), "negate": False}

    def test_take_named(self):
        """returns collections containing named seqs"""
        select = sample.take_named_seqs("a", "b")
        alns = [
            cogent3.make_aligned_seqs(
                [
                    ("a", "GCAAGCGTTTAT"),
                    ("b", "GCTTTTGTCAAT"),
                    ("c", "GC--GCGTTTAT"),
                    ("d", "GCAAGCNNTTAT"),
                ],
                moltype=DNA,
            ),
            cogent3.make_aligned_seqs(
                [
                    ("a", "GGAAGCGTTTAT"),
                    ("b", "GCTTTTGTCAAT"),
                    ("c", "GC--GCGTTTAT"),
                    ("d", "GCAAGCNNTTAT"),
                ],
                moltype=DNA,
            ),
        ]
        got = [aln.to_dict() for aln in map(select, alns) if aln]
        expected = [
            {"a": "GCAAGCGTTTAT", "b": "GCTTTTGTCAAT"},
            {"a": "GGAAGCGTTTAT", "b": "GCTTTTGTCAAT"},
        ]
        assert got == expected
        # return False if a named seq absent
        aln = cogent3.make_aligned_seqs(
            [("c", "GC--GCGTTTAT"), ("d", "GCAAGCNNTTAT")],
            moltype=DNA,
        )
        got = select(aln)
        assert not got
        assert type(got) == composable.NotCompleted

        # using negate
        select = sample.take_named_seqs("c", negate=True)
        alns = [
            cogent3.make_aligned_seqs(
                [("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
                moltype=DNA,
            ),
            cogent3.make_aligned_seqs(
                [
                    ("a", "GGAAGCGTTTAT"),
                    ("b", "GCTTTTGTCAAT"),
                    ("c", "GC--GCGTTTAT"),
                ],
                moltype=DNA,
            ),
        ]
        got = [aln.to_dict() for aln in map(select, alns) if aln]
        expect = [
            dict(d)
            for d in [
                [("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
                [("a", "GGAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
            ]
        ]
        assert got == expect

    def test_fixedlength(self):
        """correctly returns data with specified length"""
        aln = cogent3.make_aligned_seqs(
            [("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
            moltype=DNA,
        )

        fl = sample.fixed_length(4)
        got = fl(aln)
        assert len(got) == 4
        fl = sample.fixed_length(9, moltype="dna")
        got = fl(aln)
        assert len(got) == 9
        assert list(got.moltype) == list(DNA)

        alns = [
            cogent3.make_aligned_seqs(
                [("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
                moltype=DNA,
            ),
            cogent3.make_aligned_seqs(
                [("a", "GGAAGCGT"), ("b", "GCTTT-GT")],
                moltype=DNA,
            ),
        ]
        fl = sample.fixed_length(9)
        got = [a for a in map(fl, alns) if a]
        assert len(got[0]) == 9
        expected = {"a": "GCAAGCGTT", "b": "GCTTTTGTC"}
        assert got[0].to_dict() == expected

        fl = sample.fixed_length(600)
        got = [a for a in map(fl, alns) if a]
        expected = []
        assert got == expected
        # returns NotCompletedResult if nothing satisifies
        got = fl(alns[0])
        assert type(got) == sample.NotCompleted

        fl = sample.fixed_length(9, random=True)
        got = fl(aln)
        assert len(got) == 9
        assert set(aln.names) == set("ab")

        # these will be just a subset as sampling one triplet
        fl = sample.fixed_length(3, random=True, motif_length=3)
        d = cogent3.make_aligned_seqs(
            [("a", "GCAAGCGTGTAT"), ("b", "GCTACTGTCAAT")],
            moltype=DNA,
        )
        expect = d.to_dict()
        got = fl(d)
        assert len(got) == 3
        for name, seq in got.to_dict().items():
            assert seq in expect[name]
        # as above, but with moltype defined
        fl = sample.fixed_length(3, random=True, motif_length=3, moltype="dna")
        got = fl(d)
        assert len(got) == 3
        for name, seq in got.to_dict().items():
            assert seq in expect[name]

        fl = sample.fixed_length(9, start=2)
        got = fl(aln)
        assert len(got) == 9
        assert got.to_dict() == aln[2:11].to_dict()

        fl = sample.fixed_length(4, start="random")
        expect = aln.to_dict()
        got = fl(aln)
        assert len(got) == 4
        for name, seq in got.to_dict().items():
            assert seq in expect[name]

    def test_omit_duplicated(self):
        """correctly drop duplicated sequences"""
        # strict omit_duplicated
        data = {
            "a": "ACGT",
            "b": "ACG-",  # identical excepting -
            "c": "ACGN",  # non-strict matches above
            "d": "ACGG",
            "e": "ACGG",
            "k": "ACGG",  # strict identical
            "f": "RAAA",
            "g": "YAAA",  # non-strict identical
            "h": "GGGG",
        }  # unique!
        seqs = cogent3.make_unaligned_seqs(data, moltype=DNA)

        # mask_degen = True : [{'a', 'c', 'b'}, {'k', 'd', 'e'},
        # {'g', 'f'}] are dupe sets. Only 'h' unique
        drop = sample.omit_duplicated(mask_degen=True, choose=None, moltype="dna")
        got = drop(seqs)
        assert got.to_dict() == {"h": "GGGG"}
        # mask_degen = False : [{'a', 'b'}, {'k', 'd', 'e'}]
        # c, f, g, h
        drop = sample.omit_duplicated(mask_degen=False, choose=None, moltype="dna")
        got = drop(seqs)
        expect = {
            "a": "ACGT",
            "b": "ACG-",
            "c": "ACGN",
            "f": "RAAA",
            "g": "YAAA",
            "h": "GGGG",
        }
        assert got.to_dict() == expect

    def test_omit_duplicated_aligned(self):
        """omit_duplicated works on aligned sequences"""
        data = {
            "a": "ACGT",
            "b": "ACG-",  # identical excepting -
            "c": "ACGN",  # non-strict matches above
            "d": "ACGG",
            "e": "ACGG",
            "k": "ACGG",  # strict identical
            "f": "RAAA",
            "g": "YAAA",  # non-strict identical
            "h": "GGGG",
        }  # unique!
        # choose longest
        seqs = cogent3.make_aligned_seqs(data, moltype=DNA)
        drop = sample.omit_duplicated(mask_degen=True, choose="longest", moltype="dna")
        got = drop(seqs)
        expect = {"a": "ACGT", "k": "ACGG", "g": "YAAA", "h": "GGGG"}
        assert got.to_dict() == expect

        # choose random
        drop = sample.omit_duplicated(mask_degen=True, choose="random", moltype="dna")
        got1 = drop(seqs)
        seqnames = set(got1.names)
        duplicates = [{"a", "c", "b"}, {"k", "d", "e"}, {"g", "f"}]
        # should only be one of each group
        for dupes in duplicates:
            assert len(dupes & seqnames) == 1

    def test_concat(self):
        """returns concatenated alignment"""
        alns = [
            cogent3.make_aligned_seqs(d, moltype=DNA)
            for d in [
                {"seq1": "AAA", "seq2": "AAA", "seq3": "AAA"},
                {"seq1": "TTT", "seq2": "TTT", "seq3": "TTT", "seq4": "TTT"},
                {"seq1": "CC", "seq2": "CC", "seq3": "CC"},
            ]
        ]
        ccat = sample.concat(intersect=True)
        got = ccat(alns)
        assert got.to_dict() == {
            "seq1": "AAATTTCC",
            "seq2": "AAATTTCC",
            "seq3": "AAATTTCC",
        }

        ccat = sample.concat(intersect=False)
        got = ccat(alns)
        assert got.to_dict() == {
            "seq1": "AAATTTCC",
            "seq2": "AAATTTCC",
            "seq3": "AAATTTCC",
            "seq4": "???TTT??",
        }

    def test_concat_handles_moltype(self):
        """coerces to type"""
        alns = [
            cogent3.make_aligned_seqs(d, moltype=DNA)
            for d in [
                {"seq1": "AAA", "seq2": "AAA", "seq3": "AAA"},
                {"seq1": "TTT", "seq2": "TTT", "seq3": "TTT", "seq4": "TTT"},
                {"seq1": "CC", "seq2": "CC", "seq3": "CC"},
            ]
        ]
        ccat = sample.concat()
        got = ccat(alns)
        assert isinstance(got.moltype, type(DNA))

    def test_concat_validates_type(self):
        """raises TypeError if not known alignment type"""
        data = [
            {"seq1": "AAA", "seq2": "AAA", "seq3": "AAA"},
            cogent3.make_aligned_seqs(
                {"seq1": "TTT", "seq2": "TTT", "seq3": "TTT", "seq4": "TTT"},
                moltype=DNA,
            ),
        ]
        ccat = sample.concat()
        # triggered by first record
        got = ccat(data)
        assert isinstance(got, composable.NotCompleted)

        # triggered by second record
        got = ccat(data[::-1])
        assert isinstance(got, composable.NotCompleted)

    def test_trim_stop_codons(self):
        """trims stop codons using the specified genetic code"""
        trimmer = sample.trim_stop_codons()  # defaults to standard code
        seqs = cogent3.make_unaligned_seqs(
            {"seq1": "AAATTTCCC", "seq2": "AAATTTTAA"},
            moltype="dna",
        )
        got = trimmer(seqs)
        expect = {"seq1": "AAATTTCCC", "seq2": "AAATTT"}
        assert got.to_dict() == expect

        trimmer = sample.trim_stop_codons(gc=1)  # standard code
        seqs = cogent3.make_unaligned_seqs(
            {"seq1": "AAATTTCCC", "seq2": "AAATTTTAA"},
            moltype="dna",
        )
        got = trimmer(seqs)
        expect = {"seq1": "AAATTTCCC", "seq2": "AAATTT"}
        assert got.to_dict() == expect
        trimmer = sample.trim_stop_codons(gc=1)  # standard code
        aln = cogent3.make_aligned_seqs(
            {"seq1": "AAATTTCCC", "seq2": "AAATTTTAA"},
            moltype="dna",
        )
        got = trimmer(aln)
        expect = {"seq1": "AAATTTCCC", "seq2": "AAATTT---"}
        assert got.to_dict() == expect

        # different genetic code
        trimmer = sample.trim_stop_codons(gc=2)  # mt code
        seqs = cogent3.make_unaligned_seqs(
            {"seq1": "AAATTTCCC", "seq2": "AAATTTAGA"},
            moltype="dna",
        )
        got = trimmer(seqs)
        expect = {"seq1": "AAATTTCCC", "seq2": "AAATTT"}
        assert got.to_dict() == expect

    def test_take_n_seqs(self):
        """select specified number of sequences from a collection"""
        seqs1 = cogent3.make_unaligned_seqs(
            {
                "a": "ACGT",
                "b": "ACG-",
                "c": "ACGN",
                "d": "ACGG",
                "e": "ACGG",
                "k": "ACGG",
                "f": "RAAA",
                "g": "YAAA",
                "h": "GGGG",
            },
            moltype="dna",
        )
        seqs2 = seqs1.take_seqs(["a", "c", "e", "g", "h"])

        # by order, fixed
        take = sample.take_n_seqs(3, fixed_choice=True)
        got = take(seqs1)
        assert len(got.names) == 3
        # this should return NotCompleted because it applies the names present in 1 to the next one
        got = take(seqs2)
        assert isinstance(got, NotCompleted)

        take = sample.take_n_seqs(30)
        # this should fail because too few seqs
        got = take(seqs1)
        assert isinstance(got, NotCompleted)

        # by order, not fixed
        take = sample.take_n_seqs(3, fixed_choice=False)
        got1 = take(seqs1)
        got2 = take(seqs2)
        assert set(got1.names) != set(got2.names)

        # random choice, fixed
        take = sample.take_n_seqs(3, random=True, fixed_choice=True)
        assert take._fixed_choice is True

        got1 = take(seqs2)
        got2 = take(seqs1)
        assert got1.names == got2.names

        # random choice, not fixed
        take = sample.take_n_seqs(2, random=True, fixed_choice=False)
        assert take._fixed_choice is False
        # testing this is hard, we simply expect the labels to differ on subsequent call
        # the probability of drawing a specific pair of names on one call is 1/(9 choose 2) = 1/36
        # at n = 11, the probability all the pairs will be identical is ~=0
        first_call = take(seqs1)
        for _ in range(11):
            got = take(seqs1)
            different = first_call.names != got.names
            if different:
                break

        assert different, "failed to generate different random sample"

        # try setting the seed
        take = sample.take_n_seqs(2, random=True, seed=123)
        got = take(seqs1)
        assert not isinstance(got, NotCompleted)


def test_concat_coerced_moltype():
    # moltype of final result is the first one seen
    concat = sample.concat()
    aln1 = cogent3.make_aligned_seqs(
        {"s1": "AAA", "s2": "CAA", "s3": "AAA"},
        moltype="dna",
    )
    aln2 = cogent3.make_aligned_seqs(
        {"s1": "GCG", "s2": "GGG", "s3": "GGT"},
        moltype=DNA,
    )
    result = concat([aln1, aln2])
    assert result.moltype.label == "dna"


@pytest.mark.parametrize(
    "data",
    [[], [cogent3.make_aligned_seqs({"s1": "", "s2": ""}, moltype=DNA)]],
)
def test_concat_empty(data):
    # triggered by empty alignment
    ccat = sample.concat()
    got = ccat(data)
    assert isinstance(got, NotCompleted)


@pytest.fixture
def bad_gap_data():
    return {
        "s1": "---ACC---TT-",
        "s2": "---ACC---TT-",
        "s3": "---ACC---TT-",
        "s4": "--AACCG-GTT-",
        "s5": "--AACCGGGTTT",
        "s6": "AGAACCGGGTT-",
        "s7": "------------",
    }


@pytest.mark.parametrize(
    ("gap_fraction", "quantile", "expected_keys"),
    [
        (
            1,
            None,
            {"s1", "s2", "s3", "s4", "s5", "s6"},
        ),  # default just eliminates strict gap sequences
        (0.5, None, {"s4", "s5", "s6"}),  # providing a more stringent gap_frac
        (
            1,
            6 / 7,
            {"s1", "s2", "s3", "s4", "s5"},
        ),  # setting quantile drops additional sequences
    ],
)
def test_omit_bad_seqs(bad_gap_data, gap_fraction, quantile, expected_keys):
    """correctly omit bad sequences from an alignment"""

    dropbad = sample.omit_bad_seqs(gap_fraction=gap_fraction, quantile=quantile)
    aln = cogent3.make_aligned_seqs(bad_gap_data, moltype=DNA)
    got = dropbad(aln)
    expected = {k: bad_gap_data[k] for k in expected_keys}
    assert got.to_dict() == expected


def test_omit_bad_seqs_error():
    with pytest.raises(MolTypeError):
        sample.omit_bad_seqs(moltype="text")


@pytest.fixture
def bad_ambig_gap_data():
    return {
        "s1": "---ACC---TT-",
        "s2": "---ACC---TT-",
        "s3": "SSSACCSSSTTS",
        "s4": "SSAACCGSGTTS",
        "s5": "--AACCGGGTTT",
        "s6": "AGAACCGGGTT-",
        "s7": "SSSSSSSSSSSS",
    }


def test_omit_bad_seqs_ambigs(bad_ambig_gap_data):
    """correctly omit sequences from a new style alignment via fraction of ambiguous data"""
    from cogent3.core import alignment as c3_alignment

    aln = c3_alignment.make_aligned_seqs(bad_ambig_gap_data, moltype="dna")

    # drop sequences with all ambiguous data
    dropbad = sample.omit_bad_seqs(ambig_fraction=1)
    got = dropbad(aln)
    assert set(got.to_dict().keys()) == {"s1", "s2", "s3", "s4", "s5", "s6"}

    # drop sequences with more than 50% ambiguous data
    dropbad = sample.omit_bad_seqs(ambig_fraction=0.5)
    got = dropbad(aln)
    assert set(got.to_dict().keys()) == {"s1", "s2", "s4", "s5", "s6"}

    # drop sequences with any ambiguous data
    dropbad = sample.omit_bad_seqs(ambig_fraction=0.01)
    got = dropbad(aln)
    assert set(got.to_dict().keys()) == {"s1", "s2", "s5", "s6"}

    # drop sequences with more than 50% ambiguous data and 50% gaps
    dropbad = sample.omit_bad_seqs(gap_fraction=0.5, ambig_fraction=0.5)
    got = dropbad(aln)
    assert set(got.to_dict().keys()) == {"s4", "s5", "s6"}


def test_filter_degen():
    """just_nucs correctly identifies data with only nucleotides"""
    aln = cogent3.make_aligned_seqs(
        [("a", "ACGA-GACG"), ("b", "GATGATGYT")],
        moltype=DNA,
    )
    degen = sample.omit_degenerates(moltype="dna")
    got = degen(aln)  # pylint: disable=not-callable
    assert got.to_dict() == {"a": "ACGAGAG", "b": "GATGTGT"}
    # not doing isinstance during transition to new_type
    assert got.__class__.__name__.endswith("Alignment")

    # no ungapped columns
    aln = cogent3.make_aligned_seqs(
        [("a", "-C-A-G-C-"), ("b", "G-T-A-G-T")],
        moltype=DNA,
    )
    got = degen(aln)
    assert isinstance(got, composable.NotCompleted)

    # we get back the alignment type we passed in
    aln = cogent3.make_aligned_seqs(
        [("a", "ACGA-GACG"), ("b", "GATGATGYT")],
        moltype=DNA,
    )
    got = degen(aln)
    # not doing isinstance during transition to new_type
    assert got.__class__.__name__.endswith("Alignment")

    # motif length exludes tuples with a degenerate site
    aln = cogent3.make_aligned_seqs(
        {"a": "ACGA-GACG", "b": "GATGATGYT"},
        moltype=DNA,
    )
    degen = sample.omit_degenerates(moltype="dna", motif_length=2)
    got = degen(aln)
    expect = cogent3.make_aligned_seqs({"a": "ACGA", "b": "GATG"}, moltype="dna")
    assert got.to_dict() == expect.to_dict()


def test_omit_gapped():
    """omit_gap_pos correctly drops aligned columns"""
    # array alignment
    data = [("a", "ACGA-GA-CG"), ("b", "GATGATG-AT")]
    aln = cogent3.make_aligned_seqs(data, moltype=DNA)
    nogaps = sample.omit_gap_pos(moltype="dna", allowed_frac=0)  # default
    got = nogaps(aln)  # pylint: disable=not-callable
    # not doing isinstance during transition to new_type
    assert got.__class__.__name__.endswith("Alignment")
    expect = {"a": "ACGAGACG", "b": "GATGTGAT"}
    assert got.to_dict() == expect
    # standard alignment
    aln = cogent3.make_aligned_seqs(data, moltype=DNA)
    got = nogaps(aln)  # pylint: disable=not-callable
    assert got.__class__.__name__.endswith("Alignment")
    assert got.to_dict() == expect
    # non-exclusive gaps
    not_all_gaps = sample.omit_gap_pos(moltype="dna")  # default
    expect = {"a": "ACGA-GACG", "b": "GATGATGAT"}
    aln = cogent3.make_aligned_seqs(data, moltype=DNA)
    got = not_all_gaps(aln)  # pylint: disable=not-callable
    assert got.to_dict() == expect
    aln = cogent3.make_aligned_seqs(data, moltype=DNA)
    got = not_all_gaps(aln)  # pylint: disable=not-callable
    assert got.to_dict() == expect
    # with motif length
    not_all_gaps = sample.omit_gap_pos(
        moltype="dna",
        allowed_frac=0,
        motif_length=2,
    )
    aln = cogent3.make_aligned_seqs(data, moltype=DNA)
    expect = {"a": "ACGACG", "b": "GATGAT"}
    got = not_all_gaps(aln)  # pylint: disable=not-callable
    assert got.to_dict() == expect

    # no ungapped columns returns NotCompleted
    aln = cogent3.make_aligned_seqs(
        [("a", "-C-A-G-C-"), ("b", "G-T-A-G-T")],
        moltype=DNA,
    )
    got = nogaps(aln)  # pylint: disable=not-callable
    assert isinstance(got, composable.NotCompleted)


def test_minlength():
    """correctly identifies data with minimal length"""
    aln = cogent3.make_aligned_seqs(
        [("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
        moltype="text",
    )

    # if using subtract_degen, fails if incorect moltype
    ml = sample.min_length(9, subtract_degen=True)
    got = ml(aln)  # pylint: disable=not-callable
    assert isinstance(got, composable.NotCompleted)
    assert got.type == "ERROR"

    # but works if subtract_degen is False
    ml = sample.min_length(9, subtract_degen=False)
    aln = ml(aln)  # pylint: disable=not-callable
    assert len(aln) == 12
    # or if moltype provided
    ml = sample.min_length(9, subtract_degen=True, moltype="dna")
    aln = ml(aln)  # pylint: disable=not-callable
    assert len(aln) == 12

    alns = [
        cogent3.make_aligned_seqs(
            [("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
            moltype=DNA,
        ),
        cogent3.make_aligned_seqs(
            [("a", "GGAAGCGT"), ("b", "GCTTT-GT")],
            moltype=DNA,
        ),
    ]
    ml = sample.min_length(9)
    got = [aln.to_dict() for aln in map(ml, alns) if aln]  # pylint: disable=not-callable
    expected = [{"a": "GCAAGCGTTTAT", "b": "GCTTTTGTCAAT"}]
    assert got == expected

    # returns NotCompletedResult if nothing satisifies
    got = ml(alns[1])  # pylint: disable=not-callable
    assert isinstance(got, sample.NotCompleted)

    alns = [
        cogent3.make_unaligned_seqs(
            [("a", "GGAAGCGT"), ("b", "GCTTNGT")],
            moltype=DNA,
        ),
    ]
    ml = sample.min_length(6)
    got = [aln.to_dict() for aln in map(ml, alns) if aln]  # pylint: disable=not-callable
    expected = [{"a": "GGAAGCGT", "b": "GCTTNGT"}]
    assert got == expected

    ml = sample.min_length(7)
    got = [aln.to_dict() for aln in map(ml, alns) if aln]  # pylint: disable=not-callable
    expected = []
    assert got == expected


def test_codon_positions_4fold_degen():
    """codon_positions correctly return fourfold degenerate bases"""
    #                           **4---**4---
    aln = cogent3.make_aligned_seqs(
        [("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")],
        moltype=DNA,
    )
    expect = {"a": "AT", "b": "TC"}
    ffold = sample.take_codon_positions(fourfold_degenerate=True)
    got = ffold(aln)
    assert got.to_dict() == expect
    # error if no moltype
    with pytest.raises(AssertionError):
        _ = sample.take_codon_positions(moltype=None)


def test_fourfold_empty_alignment():
    aln = cogent3.make_aligned_seqs(
        {"s1": "ACGACGACG", "s2": "GATGATGAT"}, moltype="dna"
    )
    take_fourfold = cogent3.get_app(
        "take_codon_positions",
        fourfold_degenerate=True,
        moltype="dna",
    )
    result = take_fourfold.main(aln)
    assert isinstance(result, NotCompleted)
    assert result.message == "result is empty"
