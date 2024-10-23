import itertools
import json
import os
import pathlib
import re
from os import remove
from tempfile import TemporaryDirectory, mktemp
from unittest import TestCase

import numpy
import pytest
from numpy import array, log2, nan
from numpy.testing import assert_allclose, assert_equal

from cogent3 import (
    get_app,
    get_code,
    get_moltype,
    load_aligned_seqs,
    load_seq,
    load_unaligned_seqs,
    make_aligned_seqs,
    make_seq,
    make_unaligned_seqs,
    open_,
)
from cogent3._version import __version__
from cogent3.core.alignment import (
    Aligned,
    Alignment,
    ArrayAlignment,
    DataError,
    SequenceCollection,
    _coerce_to_unaligned_seqs,
    _construct_unaligned_seq,
    _SequenceCollectionBase,
    make_gap_filter,
)
from cogent3.core.alphabet import AlphabetError
from cogent3.core.annotation_db import GffAnnotationDb
from cogent3.core.moltype import AB
from cogent3.core.sequence import ArraySequence, Sequence, SeqView
from cogent3.maths.util import safe_p_log_p
from cogent3.parse.fasta import MinimalFastaParser
from cogent3.util.misc import get_object_provenance

DNA = get_moltype("dna")
RNA = get_moltype("rna")
PROTEIN = get_moltype("protein")
ASCII = get_moltype("text")
BYTES = get_moltype("bytes")


class SequenceCollectionBaseTests(object):
    """base class for testing the SequenceCollection object.

    Unlike Alignments, SequenceCollections can have sequences that are not equal
    length. This module contains all the code that _doesn't_ depend on being
    able to look at "ragged" SequenceCollections. It is intended that all
    classes that inherit from SequenceCollection should have test classes that
    inherit from this class, but that the SequenceCollection tests themselves
    will additionally contain code to deal with SequenceCollections of unequal
    length.

    set self.Class in subclasses to generate the right constructor.
    """

    Class = SequenceCollection
    brca1_data = load_aligned_seqs("data/brca1.fasta").to_dict()

    def setUp(self):
        """Define some standard SequenceCollection objects."""
        self.one_seq = self.Class({"a": "AAAAA"})
        self.ragged_padded = self.Class({"a": "AAAAAA", "b": "AAA---", "c": "AAAA--"})
        self.identical = self.Class({"a": "AAAA", "b": "AAAA"})
        self.gaps = self.Class({"a": "AAAAAAA", "b": "A--A-AA", "c": "AA-----"})
        self.gaps_rna = self.Class(
            {
                "a": RNA.make_seq(seq="AAAAAAA"),
                "b": RNA.make_seq(seq="A--A-AA"),
                "c": RNA.make_seq(seq="AA-----"),
            }
        )
        self.unordered = self.Class({"a": "AAAAA", "b": "BBBBB"})
        self.ordered1 = self.Class({"a": "AAAAA", "b": "BBBBB"}, names=["a", "b"])
        self.ordered2 = self.Class({"a": "AAAAA", "b": "BBBBB"}, names=["b", "a"])
        self.mixed = self.Class({"a": "ABCDE", "b": "LMNOP"})
        self.end_gaps = self.Class(
            {"a": "--A-BC-", "b": "-CB-A--", "c": "--D-EF-"}, names=["a", "b", "c"]
        )
        self.many = self.Class(
            {
                "a": RNA.make_seq(seq="UCAGUCAGUU"),
                "b": RNA.make_seq(seq="UCCGUCAAUU"),
                "c": RNA.make_seq(seq="ACCAUCAGUC"),
                "d": RNA.make_seq(seq="UCAAUCGGUU"),
                "e": RNA.make_seq(seq="UUGGUUGGGU"),
                "f": RNA.make_seq(seq="CCGGGCGGCC"),
                "g": RNA.make_seq(seq="UCAACCGGAA"),
            }
        )
        # Additional SequenceCollections for tests added 6/4/04 by Jeremy
        # Widmann
        self.sequences = self.Class(
            [
                RNA.make_seq(seq="UCAG"),
                RNA.make_seq(seq="UCAG"),
                RNA.make_seq(seq="UCAG"),
            ]
        )
        # Additional SequenceCollection for tests added 1/30/06 by Cathy
        # Lozupone
        self.omitSeqsTemplate_aln = self.Class(
            {
                "s1": RNA.make_seq(seq="UC-----CU---C"),
                "s2": RNA.make_seq(seq="UC------U---C"),
                "s3": RNA.make_seq(seq="UUCCUUCUU-UUC"),
                "s4": RNA.make_seq(seq="UU-UUUU-UUUUC"),
                "s5": RNA.make_seq(seq="-------------"),
            }
        )

        self.a = ArrayAlignment(["AAA", "AAA"])
        self.b = Alignment(["AAA", "AAA"])
        self.c = SequenceCollection(["AAA", "AAA"])

    def test_init_aln(self):
        """SequenceCollection should init from existing alignments"""
        exp = self.Class(["AAA", "AAA"])
        x = self.Class(self.a)
        y = self.Class(self.b)
        z = self.Class(self.c)
        assert x == exp
        assert z == exp
        assert y == exp

    test_init_aln.__doc__ = Class.__name__ + test_init_aln.__doc__

    def test_names_attribute(self):  # ported
        """expected to be a list"""
        seqs = self.Class({"a": b"AAAAA", "b": b"BBBBB"}, names=("a", "b"))
        self.assertIsInstance(seqs.names, list)

    def test_init_name_mapped(self):  # ported
        """SequenceCollection init should allow name mapping function"""
        d = {"a": "AAAAA", "b": "BBBBB"}

        def f(x):
            return x.upper()

        a = self.Class(d, label_to_name=f)
        self.assertNotEqual(a, d)
        self.assertNotEqual(sorted(a.named_seqs.items()), sorted(d.items()))
        d_upper = {"A": "AAAAA", "B": "BBBBB"}
        self.assertEqual(a, d_upper)
        self.assertEqual(sorted(a.named_seqs.items()), sorted(d_upper.items()))

    def test_init_seq(self):  # ported
        """SequenceCollection init from list of sequences should use indices as keys"""
        seqs = ["AAAAA", "BBBBB", "CCCCC"]
        a = self.Class(seqs)
        self.assertEqual(len(a.named_seqs), 3)
        self.assertEqual(a.named_seqs["seq_0"], "AAAAA")
        self.assertEqual(a.named_seqs["seq_1"], "BBBBB")
        self.assertEqual(a.named_seqs["seq_2"], "CCCCC")
        self.assertEqual(a.names, ["seq_0", "seq_1", "seq_2"])
        self.assertEqual(list(a.seqs), ["AAAAA", "BBBBB", "CCCCC"])

    def test_init_pairs(self):  # ported
        """SequenceCollection init from list of (key,val) pairs should work correctly"""
        seqs = [["x", "XXX"], ["b", "BBB"], ["c", "CCC"]]
        a = self.Class(seqs)
        self.assertEqual(len(a.named_seqs), 3)
        self.assertEqual(a.named_seqs["x"], "XXX")
        self.assertEqual(a.named_seqs["b"], "BBB")
        self.assertEqual(a.named_seqs["c"], "CCC")
        self.assertEqual(a.names, ["x", "b", "c"])
        self.assertEqual(list(a.seqs), ["XXX", "BBB", "CCC"])

    def test_init_ordered(self):  # ported
        """SequenceCollection should iterate over seqs correctly even if ordered"""
        first = self.ordered1
        sec = self.ordered2
        un = self.unordered

        self.assertEqual(first.names, ["a", "b"])
        self.assertEqual(sec.names, ["b", "a"])
        self.assertEqual(set(un.names), set(un.named_seqs.keys()))

        first_list = list(first.seqs)
        sec_list = list(sec.seqs)
        un_list = list(un.seqs)

        self.assertEqual(first_list, ["AAAAA", "BBBBB"])
        self.assertEqual(sec_list, ["BBBBB", "AAAAA"])

        # check that the unordered seq matches one of the lists
        self.assertTrue((un_list == first_list) or (un_list == sec_list))
        self.assertNotEqual(first_list, sec_list)

    def test_init_ambig(self):  # ported
        """SequenceCollection should tolerate ambiguous chars"""
        aln = self.Class(["AAA", "CCC"], moltype=DNA)
        aln = self.Class(["ANS", "CWC"], moltype=DNA)
        aln = self.Class(["A-A", "CC-"], moltype=DNA)
        aln = self.Class(["A?A", "CC-"], moltype=DNA)

    def test_seq_len_get(self):  # ported
        """SequenceCollection seq_len should return length of longest seq"""
        self.assertEqual(self.one_seq.seq_len, 5)
        self.assertEqual(self.identical.seq_len, 4)
        self.assertEqual(self.gaps.seq_len, 7)

    def test_Seqs(self):  # ported
        """SequenceCollection seqs property should return seqs in correct order."""
        first = self.ordered1
        sec = self.ordered2
        un = self.unordered

        first_list = list(first.seqs)
        sec_list = list(sec.seqs)
        un_list = list(un.seqs)

        self.assertEqual(first_list, ["AAAAA", "BBBBB"])
        self.assertEqual(sec_list, ["BBBBB", "AAAAA"])

        # check that the unordered seq matches one of the lists
        self.assertTrue((un_list == first_list) or (un_list == sec_list))
        self.assertNotEqual(first_list, sec_list)

    def test_iter_seqs(self):  # ported
        """SequenceCollection iter_seqs() method should support reordering of seqs"""
        self.ragged_padded = self.Class(
            self.ragged_padded.named_seqs, names=["a", "b", "c"]
        )
        seqs = list(self.ragged_padded.iter_seqs())
        self.assertEqual(seqs, ["AAAAAA", "AAA---", "AAAA--"])
        seqs = list(self.ragged_padded.iter_seqs(seq_order=["b", "a", "a"]))
        self.assertEqual(seqs, ["AAA---", "AAAAAA", "AAAAAA"])
        self.assertIs(seqs[1], seqs[2])
        self.assertIs(seqs[0], self.ragged_padded.named_seqs["b"])

    def test_Items(self):  # will not port
        """SequenceCollection iter_selected should iterate over items in specified order."""
        # should work if one row
        self.assertEqual(list(self.one_seq.iter_selected()), ["A"] * 5)
        # should take order into account
        self.assertEqual(list(self.ordered1.iter_selected()), ["A"] * 5 + ["B"] * 5)
        self.assertEqual(list(self.ordered2.iter_selected()), ["B"] * 5 + ["A"] * 5)

    def test_take_seqs(self):  # ported
        """SequenceCollection take_seqs should return new SequenceCollection with selected seqs."""
        a = self.ragged_padded.take_seqs(list("bc"))
        self.assertTrue(isinstance(a, _SequenceCollectionBase))
        self.assertEqual(a, {"b": "AAA---", "c": "AAAA--"})
        # should be able to negate
        a = self.ragged_padded.take_seqs(list("bc"), negate=True)
        self.assertEqual(a, {"a": "AAAAAA"})

    def test_take_seqs_str(self):  # ported
        """string arg to SequenceCollection take_seqs should work."""
        a = self.ragged_padded.take_seqs("a", negate=True)
        self.assertTrue(isinstance(a, _SequenceCollectionBase))
        self.assertEqual(a, {"b": "AAA---", "c": "AAAA--"})
        # should be able to negate
        a = self.ragged_padded.take_seqs("a")
        self.assertEqual(a, {"a": "AAAAAA"})

    def test_take_seqs_info(self):  # ported
        """take_seqs should preserve info attribute"""
        orig = self.Class(
            data={"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"}, info={"key": "value"}
        )
        subset = orig.take_seqs(list("ab"))
        self.assertEqual(set(subset.info), set(orig.info))

    def test_take_seqs_moltype(self):  # ported
        """take_seqs should preserve the MolType"""
        orig = self.Class(
            data={"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"}, moltype=DNA
        )
        subset = orig.take_seqs(list("ab"))
        self.assertEqual(set(subset.moltype), set(orig.moltype))

    def test_get_seq_indices(self):  # ported
        """SequenceCollection get_seq_indices should return names of seqs where f(row) is True"""
        srp = self.ragged_padded

        def is_long(x):
            return len(x) > 10

        def is_med(x):
            return len(str(x).replace("-", "")) > 3  # strips gaps

        def is_any(x):
            return len(x) > 0

        self.assertEqual(srp.get_seq_indices(is_long), [])
        srp.names = "cba"
        self.assertEqual(srp.get_seq_indices(is_med), ["c", "a"])
        srp.names = "bac"
        self.assertEqual(srp.get_seq_indices(is_med), ["a", "c"])
        self.assertEqual(srp.get_seq_indices(is_any), ["b", "a", "c"])
        # should be able to negate
        self.assertEqual(srp.get_seq_indices(is_med, negate=True), ["b"])
        self.assertEqual(srp.get_seq_indices(is_any, negate=True), [])

    def test_take_seqs_if(self):  # ported
        """SequenceCollection take_seqs_if should return seqs where f(row) is True"""

        def is_long(x):
            return len(x) > 10

        def is_med(x):
            return len(str(x).replace("-", "")) > 3

        def is_any(x):
            return len(x) > 0

        srp = self.ragged_padded
        self.assertEqual(srp.take_seqs_if(is_long), {})
        srp.names = "cba"
        self.assertEqual(srp.take_seqs_if(is_med), {"c": "AAAA--", "a": "AAAAAA"})
        srp.names = list(srp.named_seqs.keys())
        self.assertEqual(srp.take_seqs_if(is_med), {"c": "AAAA--", "a": "AAAAAA"})
        self.assertEqual(srp.take_seqs_if(is_any), srp)
        self.assertTrue(isinstance(srp.take_seqs_if(is_med), _SequenceCollectionBase))
        # should be able to negate
        self.assertEqual(srp.take_seqs_if(is_med, negate=True), {"b": "AAA---"})

    def test_get_identical_sets(self):  # ported
        """correctly identify sets of identical sequences"""
        from warnings import catch_warnings, filterwarnings

        # for DNA
        data = {
            "a": "ACGG",
            "b": "ACGG",  # strict identical
            "c": "ACGN",  # non-strict matches above
            "d": "ACGT",
            "e": "ACGT",
            "k": "ACGT",  # strict identical
            "f": "RAAA",
            "g": "YAAA",
        }  # non-strict identical

        seqs = self.Class(data=data, moltype=DNA)
        got = seqs.get_identical_sets(mask_degen=False)  # a strict comparison
        # convert to frozenset, so we can do a comparison robust to set order
        got = frozenset(frozenset(s) for s in got)
        expect = [{"a", "b"}, {"d", "e", "k"}]
        expect = frozenset(frozenset(s) for s in expect)
        self.assertEqual(got, expect)

        got = seqs.get_identical_sets(mask_degen=True)
        got = frozenset(frozenset(s) for s in got)
        expect = [{"a", "b", "c"}, {"d", "e", "k"}, {"f", "g"}]
        expect = frozenset(frozenset(s) for s in expect)
        self.assertEqual(got, expect)

        # for PROTEIN
        data = {
            "a": "ACGT",
            "b": "ACGT",  # strict identical
            "c": "ACGX",  # non-strict matches above
            "d": "TTTT",
            "e": "TTTT",
            "k": "TTTT",  # strict identical
            "f": "BAAA",
            "g": "ZAAA",
        }  # non-strict identical

        seqs = self.Class(data=data, moltype=PROTEIN)
        got = seqs.get_identical_sets(mask_degen=False)  # a strict comparison
        # convert to frozenset, so we can do a comparison robust to set order
        got = frozenset(frozenset(s) for s in got)
        expect = [{"a", "b"}, {"d", "e", "k"}]
        expect = frozenset(frozenset(s) for s in expect)
        self.assertEqual(got, expect)

        got = seqs.get_identical_sets(mask_degen=True)
        got = frozenset(frozenset(s) for s in got)
        expect = [{"a", "b", "c"}, {"d", "e", "k"}, {"f", "g"}]
        expect = frozenset(frozenset(s) for s in expect)
        self.assertEqual(got, expect)

        # if the moltype has no degen characters, just return for mask_degen
        seqs = self.Class(data=data, moltype=ASCII)
        got = seqs.get_identical_sets(mask_degen=False)
        # convert to frozenset, so we can do a comparison robust to set order
        got = frozenset(frozenset(s) for s in got)
        expect = [{"a", "b"}, {"d", "e", "k"}]
        expect = frozenset(frozenset(s) for s in expect)
        self.assertEqual(got, expect)

        with catch_warnings():
            filterwarnings("ignore", category=UserWarning)
            got = seqs.get_identical_sets(mask_degen=True)

        got = frozenset(frozenset(s) for s in got)
        self.assertEqual(got, expect)

    def test_is_ragged(self):  # ported
        """SequenceCollection is_ragged should return true if ragged alignment"""
        assert not self.identical.is_ragged()
        assert not self.gaps.is_ragged()

    def test_to_phylip(self):  # ported
        """SequenceCollection should return PHYLIP string format correctly"""
        align_norm = self.Class(
            [
                "ACDEFGHIKLMNPQRSTUVWY-",
                "ACDEFGHIKLMNPQRSUUVWF-",
                "ACDEFGHIKLMNPERSKUVWC-",
                "ACNEFGHIKLMNPQRS-UVWP-",
            ]
        )

        self.assertEqual(
            align_norm.to_phylip(),
            """4  22\nseq_0     ACDEFGHIKLMNPQRSTUVWY-\nseq_1     ACDEFGHIKLMNPQRSUUVWF-\nseq_2     ACDEFGHIKLMNPERSKUVWC-\nseq_3     ACNEFGHIKLMNPQRS-UVWP-\n""",
        )

    def test_to_fasta(self):  # ported
        """SequenceCollection should return correct FASTA string"""
        aln1 = self.Class(["AAA", "CCC"])
        self.assertEqual(aln1.to_fasta(), ">seq_0\nAAA\n>seq_1\nCCC\n")
        self.assertEqual(aln1.to_fasta(block_size=2), ">seq_0\nAA\nA\n>seq_1\nCC\nC\n")

        aln2 = self.Class(["GCATGCAT", "TCAGACGT"])
        self.assertEqual(aln2.to_fasta(), ">seq_0\nGCATGCAT\n>seq_1\nTCAGACGT\n")
        self.assertEqual(
            aln2.to_fasta(block_size=4), ">seq_0\nGCAT\nGCAT\n>seq_1\nTCAG\nACGT\n"
        )
        self.assertEqual(
            aln2.to_fasta(block_size=3), ">seq_0\nGCA\nTGC\nAT\n>seq_1\nTCA\nGAC\nGT\n"
        )

    def test_to_nexus(self):  # will not port
        """SequenceCollection should return correct Nexus string format"""
        align_norm = self.Class(
            [
                "ACDEFGHIKLMNPQRSTUVWY-",
                "ACDEFGHIKLMNPQRSUUVWF-",
                "ACDEFGHIKLMNPERSKUVWC-",
                "ACNEFGHIKLMNPQRS-UVWP-",
            ]
        )
        expect = (
            "#NEXUS\n\nbegin data;\n"
            "    dimensions ntax=4 nchar=22;\n"
            "    format datatype=protein interleave=yes missing=? gap=-;\n"
            "    matrix\n"
            "    seq_0    ACDEFGHIKLMNPQRSTUVWY-\n"
            "    seq_1    ACDEFGHIKLMNPQRSUUVWF-\n"
            "    seq_2    ACDEFGHIKLMNPERSKUVWC-\n"
            "    seq_3    ACNEFGHIKLMNPQRS-UVWP-\n\n    ;\nend;"
        )
        got = align_norm.to_nexus("protein")
        self.assertEqual(got, expect)

    def test_num_seqs(self):  # ported
        """SequenceCollection.num_seqs should count seqs."""
        aln = self.Class({"seq1": "ACGU", "seq2": "CGUA", "seq3": "CCGU"})
        self.assertEqual(aln.num_seqs, 3)

    def test_add_info(self):  # will not port
        """__add__ should preserve info attribute"""
        align1 = self.Class(
            {"a": "AAAA", "b": "TTTT", "c": "CCCC"}, info={"key": "foo"}
        )
        align2 = self.Class(
            {"a": "GGGG", "b": "----", "c": "NNNN"}, info={"key": "bar"}
        )
        align = align1 + align2
        self.assertEqual(align.info["key"], "foo")

    def test_add_seqs_info(self):  # ported
        """add_seqs should preserve info attribute"""
        data = [("name1", "AAA"), ("name2", "AAA"), ("name3", "AAA"), ("name4", "AAA")]
        data2 = [("name5", "BBB"), ("name6", "CCC")]
        aln = self.Class(data, info={"key": "foo"})
        aln2 = self.Class(data2, info={"key": "bar"})
        out_aln = aln.add_seqs(aln2)
        self.assertEqual(out_aln.info["key"], "foo")

    def test_write(self):  # ported
        """SequenceCollection.write should write in correct format"""
        aln = self.Class([("a", "AAAA"), ("b", "TTTT"), ("c", "CCCC")])
        fn = mktemp(suffix=".fasta")
        aln.write(fn)
        with open(fn, newline=None) as infile:
            result = infile.read()
        self.assertEqual(result, ">a\nAAAA\n>b\nTTTT\n>c\nCCCC\n")
        remove(fn)

    def test_len(self):  # ported
        """len(SequenceCollection) returns length of longest sequence"""
        aln = self.Class([("a", "AAAA"), ("b", "TTTT"), ("c", "CCCC")])
        self.assertEqual(len(aln), 4)

    def test_get_seq(self):  # ported
        """SequenceCollection.get_seq should return specified seq"""
        aln = self.Class({"seq1": "GATTTT", "seq2": "GATC??"})
        self.assertEqual(aln.get_seq("seq1"), "GATTTT")
        self.assertRaises(KeyError, aln.get_seq, "seqx")

    def test_to_dict(self):  # ported
        """SequenceCollection.to_dict should return dict of strings (not obj)"""
        aln = self.Class({"seq1": "GATTTT", "seq2": "GATC??"})
        self.assertEqual(aln.to_dict(), {"seq1": "GATTTT", "seq2": "GATC??"})
        for i in list(aln.to_dict().values()):
            assert isinstance(i, str)

    def test_get_ambiguous_positions(self):  # ported
        """SequenceCollection.get_ambiguous_positions should return pos"""
        aln = self.Class({"s1": "ATGRY?", "s2": "T-AG??"}, moltype=DNA)
        self.assertEqual(
            aln.get_ambiguous_positions(),
            {"s2": {4: "?", 5: "?"}, "s1": {3: "R", 4: "Y", 5: "?"}},
        )

    def test_degap(self):  # ported
        """SequenceCollection.degap should strip gaps from each seq"""
        aln = self.Class({"s1": "ATGRY?", "s2": "T-AG??"}, moltype=DNA)
        self.assertEqual(aln.degap(), {"s1": "ATGRY", "s2": "TAG"})

    def test_degap_info(self):  # ported
        """.degap should preserve info attributes"""
        aln = self.Class({"s1": "ATGRY?", "s2": "T-AG??"}, moltype=DNA)
        aln.info.path = "blah"
        got = aln.degap()
        self.assertEqual(got.info.path, "blah")

    def test_with_modified_termini(self):  # will not port
        """SequenceCollection.with_modified_termini should code trailing gaps as ?"""
        aln = self.Class({"s1": "AATGR--", "s2": "-T-AG?-"}, moltype=DNA)
        self.assertEqual(
            aln.with_modified_termini(), {"s1": "AATGR??", "s2": "?T-AG??"}
        )

    def test_make_gap_filter(self):
        """make_gap_filter returns f(seq) -> True if aligned ok w/ query"""
        s1 = RNA.make_seq(seq="UC-----CU---C")
        s3 = RNA.make_seq(seq="UUCCUUCUU-UUC")
        s4 = RNA.make_seq(seq="UU-UUUU-UUUUC")
        # check that the behavior is ok for gap runs
        f1 = make_gap_filter(s1, 0.9, 5)
        f3 = make_gap_filter(s3, 0.9, 5)
        # Should return False since s1 has gap run >= 5 with respect to s3
        self.assertEqual(f3(s1), False)
        # Should return False since s3 has an insertion run >= 5 to s1
        self.assertEqual(f1(s3), False)
        # Should retun True since s4 does not have a long enough gap or ins run
        self.assertEqual(f3(s4), True)
        f3 = make_gap_filter(s3, 0.9, 6)
        self.assertEqual(f3(s1), True)

        # Check that behavior is ok for gap_fractions
        f1 = make_gap_filter(s1, 0.5, 6)
        f3 = make_gap_filter(s3, 0.5, 6)
        # Should return False since 0.53% of positions are diff for gaps
        self.assertEqual(f3(s1), False)
        self.assertEqual(f1(s3), False)
        self.assertEqual(f3(s4), True)

    def test_omit_gap_seqs(self):  # will not port
        """SequenceCollection omit_gap_seqs should return alignment w/o seqs with gaps"""
        # check default params
        self.assertEqual(self.gaps.omit_gap_seqs(), self.gaps.omit_gap_seqs(0))
        # check for boundary effects
        self.assertEqual(self.gaps.omit_gap_seqs(-1), {})
        self.assertEqual(self.gaps.omit_gap_seqs(0), {"a": "AAAAAAA"})
        self.assertEqual(self.gaps.omit_gap_seqs(0.1), {"a": "AAAAAAA"})
        self.assertEqual(self.gaps.omit_gap_seqs(3.0 / 7 - 0.01), {"a": "AAAAAAA"})
        self.assertEqual(
            self.gaps.omit_gap_seqs(3.0 / 7), {"a": "AAAAAAA", "b": "A--A-AA"}
        )
        self.assertEqual(
            self.gaps.omit_gap_seqs(3.0 / 7 + 0.01), {"a": "AAAAAAA", "b": "A--A-AA"}
        )
        self.assertEqual(
            self.gaps.omit_gap_seqs(5.0 / 7 - 0.01), {"a": "AAAAAAA", "b": "A--A-AA"}
        )
        self.assertEqual(self.gaps.omit_gap_seqs(5.0 / 7 + 0.01), self.gaps)
        self.assertEqual(self.gaps.omit_gap_seqs(0.99), self.gaps)
        # check new object creation
        self.assertIsNot(self.gaps.omit_gap_seqs(0.99), self.gaps)
        self.assertTrue(
            isinstance(self.gaps.omit_gap_seqs(3.0 / 7), _SequenceCollectionBase)
        )
        # repeat tests for object that supplies its own gaps
        self.assertEqual(self.gaps_rna.omit_gap_seqs(-1), {})
        self.assertEqual(self.gaps_rna.omit_gap_seqs(0), {"a": "AAAAAAA"})
        self.assertEqual(self.gaps_rna.omit_gap_seqs(0.1), {"a": "AAAAAAA"})
        self.assertEqual(self.gaps_rna.omit_gap_seqs(3.0 / 7 - 0.01), {"a": "AAAAAAA"})
        self.assertEqual(
            self.gaps_rna.omit_gap_seqs(3.0 / 7), {"a": "AAAAAAA", "b": "A--A-AA"}
        )
        self.assertEqual(
            self.gaps_rna.omit_gap_seqs(3.0 / 7 + 0.01),
            {"a": "AAAAAAA", "b": "A--A-AA"},
        )
        self.assertEqual(
            self.gaps_rna.omit_gap_seqs(5.0 / 7 - 0.01),
            {"a": "AAAAAAA", "b": "A--A-AA"},
        )
        self.assertEqual(self.gaps_rna.omit_gap_seqs(5.0 / 7 + 0.01), self.gaps_rna)
        self.assertEqual(self.gaps_rna.omit_gap_seqs(0.99), self.gaps_rna)
        self.assertIsNot(self.gaps_rna.omit_gap_seqs(0.99), self.gaps_rna)
        self.assertTrue(
            isinstance(self.gaps_rna.omit_gap_seqs(3.0 / 7), _SequenceCollectionBase)
        )

    def test_omit_gap_runs(self):  # will not port
        """SequenceCollection omit_gap_runs should return alignment w/o runs of gaps"""
        # negative value will still let through ungapped sequences
        self.assertEqual(self.gaps.omit_gap_runs(-5), {"a": "AAAAAAA"})
        # test edge effects
        self.assertEqual(self.gaps.omit_gap_runs(0), {"a": "AAAAAAA"})
        self.assertEqual(self.gaps.omit_gap_runs(1), {"a": "AAAAAAA"})
        self.assertEqual(self.gaps.omit_gap_runs(2), {"a": "AAAAAAA", "b": "A--A-AA"})
        self.assertEqual(self.gaps.omit_gap_runs(3), {"a": "AAAAAAA", "b": "A--A-AA"})
        self.assertEqual(self.gaps.omit_gap_runs(4), {"a": "AAAAAAA", "b": "A--A-AA"})
        self.assertEqual(self.gaps.omit_gap_runs(5), self.gaps)
        self.assertEqual(self.gaps.omit_gap_runs(6), self.gaps)
        self.assertEqual(self.gaps.omit_gap_runs(1000), self.gaps)
        # test new object creation
        self.assertIsNot(self.gaps.omit_gap_runs(6), self.gaps)
        self.assertTrue(isinstance(self.gaps.omit_gap_runs(6), _SequenceCollectionBase))

    def test_consistent_gap_degen_handling(self):  # ported
        """gap degen character should be treated consistently"""
        # the degen character '?' can be a gap, so when we strip gaps it should
        # be gone too
        raw_seq = "---??-??TC-GGCG-GCA-G-GC-?-C-TAN-GCGC-CCTC-AGGA?-???-??--"
        raw_ungapped = re.sub("[-?]", "", raw_seq)
        re.sub("[N?]+", "", raw_seq)
        dna = DNA.make_seq(seq=raw_seq)

        aln = self.Class(data=[("a", dna), ("b", dna)])
        expect = self.Class(data=[("a", raw_ungapped), ("b", raw_ungapped)]).to_fasta()
        self.assertEqual(aln.degap().to_fasta(), expect)
        seqs = self.Class(data=[("a", dna), ("b", dna)])
        self.assertEqual(seqs.degap().to_fasta(), expect)

    def test_pad_seqs(self):  # ported
        """SequenceCollection pad_seqs should work on alignment."""
        # pad to max length
        padded1 = self.ragged_padded.pad_seqs()
        seqs1 = list(padded1.iter_seqs(seq_order=["a", "b", "c"]))
        self.assertEqual(list(map(str, seqs1)), ["AAAAAA", "AAA---", "AAAA--"])

        # pad to alternate length
        padded1 = self.ragged_padded.pad_seqs(pad_length=10)
        seqs1 = list(padded1.iter_seqs(seq_order=["a", "b", "c"]))
        self.assertEqual(
            list(map(str, seqs1)), ["AAAAAA----", "AAA-------", "AAAA------"]
        )

        # assertRaises error when pad_length is less than max seq length
        self.assertRaises(ValueError, self.ragged_padded.pad_seqs, 5)

    def test_to_moltype_info(self):  # ported
        """correctly convert to specified moltype"""
        data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        seqs = self.Class(data=data, info={"key": "value"})
        dna = seqs.to_moltype(DNA)
        self.assertEqual(dna.info["key"], "value")

    def test_get_lengths(self):  # ported
        """get_lengths handles motif length, allow_gaps etc.."""
        data = {"a": "AAAA??????", "b": "CCCGGG--NN"}
        coll = self.Class(data=data, moltype=DNA)
        got = coll.get_lengths()
        expect = {"a": 4, "b": 6}
        self.assertEqual(got, expect)
        got = coll.get_lengths(include_ambiguity=True)
        expect = {"a": 4, "b": 8}  # note ? is excluded as it could be a gap
        self.assertEqual(got, expect)

        got = coll.get_lengths(include_ambiguity=True, allow_gap=True)
        # note ? is excluded as it could be a gap
        expect = {"a": 10, "b": 10}
        self.assertEqual(got, expect)

    def test_strand_symmetry(self):  # ported
        """exercising strand symmetry test"""
        data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        seqs = self.Class(data, moltype=DNA)
        result = seqs.strand_symmetry()
        assert_allclose(result["seq1"].observed.array, [[3, 2], [2, 2]])
        assert_allclose(result["seq2"].observed.array, [[3, 0], [2, 1]])

    def test_rename_seqs(self):  # ported
        """successfully rename sequences"""
        data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        seqs = self.Class(data, moltype=DNA)
        new = seqs.rename_seqs(lambda x: x.upper())
        expect = {n.upper() for n in data}
        self.assertEqual(set(new.names), expect)
        # renamed contains the name map as info attribute
        self.assertEqual(new.info.name_map, {k.upper(): k for k in data})
        # if Alignment class, make sure data attribute does not have gaps
        if self.Class == Alignment:
            for seq in new.seqs:
                self.assertFalse("-" in str(seq.data))

    def test_apply_pssm(self):  # ported
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
        seqs = self.Class(data=data, moltype=DNA)
        scores = seqs.apply_pssm(path="data/sample.jaspar", show_progress=False)
        self.assertEqual(scores.shape, (len(data), len(seqs) - pwm.shape[0] + 1))
        scores = seqs.apply_pssm(pssm=pwm.to_pssm(), show_progress=False)
        self.assertEqual(scores.shape, (len(data), len(seqs) - pwm.shape[0] + 1))

        # using the names argument works to return scores in the correct order
        seqs = self.Class(
            data={"ENSMUSG00000056468": "GCCAGGGGGGAAAGGGAGAA"}, moltype=DNA
        )
        expect = []
        expect.extend(seqs.apply_pssm(pssm=pwm.to_pssm(), show_progress=False))
        seqs = self.Class(
            data={"ENSMUSG00000039616": "GCCCTTCAAATTTGGTTTCT"}, moltype=DNA
        )
        expect.extend(seqs.apply_pssm(pssm=pwm.to_pssm(), show_progress=False))
        expect = numpy.array(expect)
        seqs = self.Class(
            data={
                "ENSMUSG00000056468": "GCCAGGGGGGAAAGGGAGAA",
                "ENSMUSG00000039616": "GCCCTTCAAATTTGGTTTCT",
            },
            moltype=DNA,
        )
        got = seqs.apply_pssm(
            pssm=pwm.to_pssm(), show_progress=False, names="ENSMUSG00000056468"
        )
        assert_allclose(got, expect[:1])
        got = seqs.apply_pssm(
            pssm=pwm.to_pssm(), show_progress=False, names=["ENSMUSG00000039616"]
        )
        assert_allclose(got, expect[1:])
        got = seqs.apply_pssm(
            pssm=pwm.to_pssm(),
            show_progress=False,
            names=["ENSMUSG00000056468", "ENSMUSG00000039616"],
        )
        assert_allclose(got, expect)
        got = seqs.apply_pssm(
            pssm=pwm.to_pssm(),
            show_progress=False,
            names=["ENSMUSG00000039616", "ENSMUSG00000056468"],
        )
        assert_allclose(got, expect[::-1])

    def test_set_repr_policy_no_input(self):  # ported
        """repr_policy should remain unchanged"""
        seqs = self.Class({"a": "AAAAA"})
        seqs.set_repr_policy(num_seqs=None, num_pos=None)
        self.assertEqual(
            seqs._repr_policy,
            dict(num_seqs=10, num_pos=60, ref_name="longest", wrap=60),
        )

    def test_set_repr_policy_invalid_input(self):  # ported
        """repr_policy should remain unchanged"""
        seqs = self.Class({"a": "AAAAA"})
        invalid_args = (
            dict(num_seqs="foo", err=TypeError),
            dict(num_pos=4.2, err=TypeError),
            dict(ref_name="blah", err=ValueError),
            dict(wrap=3.1, err=TypeError),
        )
        for arg in invalid_args:
            err = arg.pop("err")
            with self.assertRaises(err):
                seqs.set_repr_policy(**arg)
            self.assertEqual(
                seqs._repr_policy,
                dict(num_seqs=10, num_pos=60, ref_name="longest", wrap=60),
            )

    def test_set_wrap_affects_repr_html(self):  # ported for seqcollection
        """the wrap argument affects the number of columns"""
        # indirectly tested via counting number of occurrences of 'class="label"'
        seqs = self.Class({"a": "AAAAA", "b": "AAA--"})
        orig = seqs._repr_html_()
        seqs.set_repr_policy(wrap=3)  # break alignment into 2
        got = seqs._repr_html_()
        token = 'class="label"'
        self.assertEqual(got.count(token), 2 * orig.count(token))

        # using environment variable
        env_name = "COGENT3_ALIGNMENT_REPR_POLICY"
        os.environ[env_name] = "wrap=2"
        seqs = self.Class({"a": "AAAAA", "b": "AAA--"})
        got = seqs._repr_html_()
        os.environ.pop(env_name, None)
        self.assertEqual(got.count(token), 3 * orig.count(token))

    def test_get_seq_entropy(self):  # ported
        """get_seq_entropy should get entropy of each seq"""
        a = self.Class(dict(a="ACCC", b="AGTA"), moltype=DNA)
        entropy = a.entropy_per_seq()
        e = 0.81127812445913283  # sum(p log_2 p) for p = 0.25, 0.75
        assert_allclose(entropy, array([e, 1.5]))

    def test_write_to_json(self):  # ported
        # test writing to json file
        aln = self.Class([("a", "AAAA"), ("b", "TTTT"), ("c", "CCCC")])
        with TemporaryDirectory(".") as dirname:
            path = str(pathlib.Path(dirname) / "sample.json")
            aln.write(path)
            with open_(path) as fn:
                got = json.loads(fn.read())
                self.assertEqual(got, aln.to_rich_dict())


class SequenceCollectionTests(SequenceCollectionBaseTests, TestCase):
    """Tests of the SequenceCollection object. Includes ragged collection tests.

    Should not test alignment-specific features.
    """

    def setUp(self):
        """Adds self.ragged for ragged collection tests."""
        self.ragged = SequenceCollection({"a": "AAAAAA", "b": "AAA", "c": "AAAA"})
        super(SequenceCollectionTests, self).setUp()

    def test_seq_len_get_ragged(self):  # ported
        """SequenceCollection seq_len get should work for ragged seqs"""
        self.assertEqual(self.ragged.seq_len, 6)

    def test_is_ragged_ragged(self):  # ported
        """SequenceCollection is_ragged should return True if ragged"""
        self.assertTrue(self.ragged.is_ragged())

    def test_Seqs_ragged(self):  # ported
        """SequenceCollection seqs should work on ragged alignment"""
        self.ragged.names = "bac"
        self.assertEqual(list(self.ragged.seqs), ["AAA", "AAAAAA", "AAAA"])

    def test_iter_seqs_ragged(self):  # ported
        """SequenceCollection iter_seqs() method should support reordering of seqs"""
        self.ragged.names = ["a", "b", "c"]
        seqs = list(self.ragged.iter_seqs())
        self.assertEqual(seqs, ["AAAAAA", "AAA", "AAAA"])
        seqs = list(self.ragged.iter_seqs(seq_order=["b", "a", "a"]))
        self.assertEqual(seqs, ["AAA", "AAAAAA", "AAAAAA"])
        self.assertIs(seqs[1], seqs[2])
        self.assertIs(seqs[0], self.ragged.named_seqs["b"])

    def test_toPHYLIP_ragged(self):  # ported
        """SequenceCollection should refuse to convert ragged seqs to phylip"""
        align_rag = self.Class(
            [
                "ACDEFGHIKLMNPQRSTUVWY-",
                "ACDEFGHIKLMNPQRSUUVWF-",
                "ACDEFGHIKLMNPERSKUVWC-",
                "ACNEFGHIKLMNUVWP-",
            ]
        )

        # no longer applicable in new implementation
        with self.assertRaises(ValueError):
            align_rag.to_phylip()

    def test_pad_seqs_ragged(self):  # ported
        """SequenceCollection pad_seqs should work on ragged alignment."""
        # pad to max length
        padded1 = self.ragged.pad_seqs()
        seqs1 = list(padded1.iter_seqs(seq_order=["a", "b", "c"]))
        self.assertEqual(list(map(str, seqs1)), ["AAAAAA", "AAA---", "AAAA--"])

        # pad to alternate length
        padded1 = self.ragged.pad_seqs(pad_length=10)
        seqs1 = list(padded1.iter_seqs(seq_order=["a", "b", "c"]))
        self.assertEqual(
            list(map(str, seqs1)), ["AAAAAA----", "AAA-------", "AAAA------"]
        )

        # assertRaises error when pad_length is less than max seq length
        self.assertRaises(ValueError, self.ragged.pad_seqs, 5)

    def test_info_source(self):  # ported
        """info.source exists if load seqs given a filename"""
        path = pathlib.Path("data/brca1.fasta")
        seqs = load_unaligned_seqs(path, moltype="dna")
        self.assertEqual(seqs.info.source, str(path))

    def test_apply_pssm2(self):  # ported
        """apply_pssm fail if ragged sequences"""
        data = {
            "ENSMUSG00000056468": "GCCAGGGGGGAAAGGGAGAA",
            "ENSMUSG00000039616": "GCCCTTCAAATTT",
        }
        seqs = self.Class(data=data, moltype=DNA)
        with self.assertRaises(AssertionError):
            _ = seqs.apply_pssm(path="data/sample.jaspar", show_progress=False)

    def test_construction(self):  # ported
        """correctly construct from list of sequences of length 2"""
        seq1 = make_seq(seq="AC", name="seq1")
        seq2 = make_seq(seq="AC", name="seq2")
        coll = SequenceCollection(data=[seq1, seq2])


def _make_filter_func(aln):
    array_align = type(aln) == ArrayAlignment
    if array_align:
        gap = aln.alphabet.with_gap_motif().to_indices("-")[0]
    else:
        gap = "-"

    def func_str(x):
        return gap not in "".join([str(s) for s in x])

    def func_arr(x):
        return (x != gap).all()

    return func_arr if array_align else func_str


class AlignmentBaseTests(SequenceCollectionBaseTests):
    """Tests of basic Alignment functionality. All Alignments should pass these.

    Note that this is not a TestCase: need to subclass to test each specific
    type of Alignment. Override self.Constructor with your alignment class
    as a constructor.
    """

    def test_alignment_quality(self):
        """check alignment method correctly invokes underlying app"""
        aln = self.Class(["AAAC", "ACGC", "AGCC", "A-TC"], moltype="dna")
        got = aln.alignment_quality(equifreq_mprobs=False)
        expect = (
            2 * log2(1 / 0.4)
            + log2(1 / (4 * 0.4))
            + (1 / 2) * log2(1 / (8 / 15))
            + (1 / 4) * log2(1 / (4 / 15))
        )
        assert_allclose(got, expect)

    def test_filter_drop_remainder(self):
        """filter allows dropping"""
        raw = {"a": "ACGACGACG", "b": "CCC---CCC", "c": "AAAA--AAA"}
        aln = self.Class(raw)
        func = _make_filter_func(aln)
        got = aln.filtered(func, motif_length=1, warn=False)
        self.assertEqual(len(got), 6)
        # raises an assertion if the length is not modulo
        with self.assertRaises(ValueError):
            # because alignment not modulo 2
            got = aln.filtered(func, motif_length=2, drop_remainder=False)
        got = aln.filtered(func, motif_length=2, drop_remainder=True, warn=False)
        self.assertEqual(len(got), 4)

    def test_take_positions_info(self):
        aln = self.Class(
            {"a": "AAAAAAA", "b": "A--A-AA", "c": "AA-----"}, info={"key": "value"}
        )
        tps = aln.take_positions([5, 4, 0])
        self.assertEqual(tps.info["key"], "value")

    def test_take_positions_if(self):
        """SequenceCollection take_positions_if should return cols where f(col) is True"""

        def gap_1st(x):
            return x[0] == "-"

        def gap_2nd(x):
            return x[1] == "-"

        def gap_3rd(x):
            return x[2] == "-"

        def is_list(x):
            return isinstance(x, list)

        self.gaps.names = "abc"

        self.assertEqual(
            self.gaps.take_positions_if(gap_1st), {"a": "", "b": "", "c": ""}
        )
        self.assertEqual(
            self.gaps.take_positions_if(gap_2nd), {"a": "AAA", "b": "---", "c": "A--"}
        )
        self.assertEqual(
            self.gaps.take_positions_if(gap_3rd),
            {"a": "AAAAA", "b": "-A-AA", "c": "-----"},
        )
        self.assertEqual(self.gaps.take_positions_if(is_list), self.gaps)
        self.assertEqual(self.gaps.take_positions_if(is_list), self.gaps)

        self.assertTrue(
            isinstance(self.gaps.take_positions_if(gap_1st), _SequenceCollectionBase)
        )
        # should be able to negate
        self.assertEqual(self.gaps.take_positions_if(gap_1st, negate=True), self.gaps)
        self.assertEqual(
            self.gaps.take_positions_if(gap_2nd, negate=True),
            {"a": "AAAA", "b": "AAAA", "c": "A---"},
        )
        self.assertEqual(
            self.gaps.take_positions_if(gap_3rd, negate=True),
            {"a": "AA", "b": "A-", "c": "AA"},
        )

    def test_omit_gap_pos2(self):
        """consistency with different motif_length values"""
        data = {
            "seq1": "CAGGTCGACCTCGGC---------CACGAC",
            "seq2": "CAGATCGACCTCGGC---------CACGAC",
            "seq3": "CAGATCGACCTCGGT---------CACGAT",
            "seq4": "CAGATCGACCTCGGCGAACACGGCCATGAT",
            "seq5": "CCGATCGACATGGGC---------CACGAT",
            "seq6": "GCC---------------------------",
        }
        aln = self.Class(data, moltype=DNA)
        got1 = aln.omit_gap_pos(motif_length=1)
        got3 = aln.omit_gap_pos(motif_length=3)
        self.assertEqual(len(got3), len(got1))
        self.assertEqual(got3.to_dict(), got1.to_dict())

    def test_omit_bad_seqs(self):
        """omit_bad_seqs should return alignment w/o seqs causing most gaps"""
        data = {
            "s1": "---ACC---TT-",
            "s2": "---ACC---TT-",
            "s3": "---ACC---TT-",
            "s4": "--AACCG-GTT-",
            "s5": "--AACCGGGTTT",
            "s6": "AGAACCGGGTT-",
        }

        aln = self.Class(data, moltype=DNA)
        # with defaults, excludes s6
        expect = data.copy()
        del expect["s6"]
        result = aln.omit_bad_seqs()
        self.assertEqual(result.to_dict(), expect)
        # with quantile 0.5, just s1, s2, s3
        expect = data.copy()
        for key in ("s6", "s5"):
            del expect[key]
        result = aln.omit_bad_seqs(0.5)
        self.assertEqual(result.to_dict(), expect)

    def test_matching_ref(self):
        """Alignment.matching_ref returns new aln with well-aln to temp"""
        aln = self.omitSeqsTemplate_aln
        result = aln.matching_ref("s3", 0.9, 5)
        self.assertEqual(result, {"s3": "UUCCUUCUU-UUC", "s4": "UU-UUUU-UUUUC"})
        result2 = aln.matching_ref("s4", 0.9, 4)
        self.assertEqual(result2, {"s3": "UUCCUUCUU-UUC", "s4": "UU-UUUU-UUUUC"})
        result3 = aln.matching_ref("s1", 0.9, 4)
        self.assertEqual(
            result3,
            {"s2": "UC------U---C", "s1": "UC-----CU---C", "s5": "-------------"},
        )
        result4 = aln.matching_ref("s3", 0.5, 13)
        self.assertEqual(result4, {"s3": "UUCCUUCUU-UUC", "s4": "UU-UUUU-UUUUC"})

    def test_iupac_consensus_RNA(self):
        """SequenceCollection iupac_consensus should use RNA IUPAC symbols correctly"""
        alignmentUpper = self.Class(
            [
                "UCAGN-UCAGN-UCAGN-UCAGAGCAUN-",
                "UUCCAAGGNN--UUCCAAGGNNAGCAG--",
                "UUCCAAGGNN--UUCCAAGGNNAGCUA--",
                "UUUUCCCCAAAAGGGGNNNN--AGCUA--",
                "UUUUCCCCAAAAGGGGNNNN--AGCUA--",
            ],
            moltype=RNA,
        )

        # following IUPAC consensus calculated by hand
        # Test all uppper
        self.assertEqual(
            alignmentUpper.iupac_consensus(), "UYHBN?BSNN??KBVSN?NN??AGCWD?-"
        )

    def test_iupac_consensus_DNA(self):
        """SequenceCollection iupac_consensus should use DNA IUPAC symbols correctly"""
        alignmentUpper = self.Class(
            [
                "TCAGN-TCAGN-TCAGN-TCAGAGCATN-",
                "TTCCAAGGNN--TTCCAAGGNNAGCAG--",
                "TTCCAAGGNN--TTCCAAGGNNAGCTA--",
                "TTTTCCCCAAAAGGGGNNNN--AGCTA--",
                "TTTTCCCCAAAAGGGGNNNN--AGCTA--",
            ]
        )
        # following IUPAC consensus calculated by hand
        # Test all uppper
        self.assertEqual(
            alignmentUpper.iupac_consensus(DNA), "TYHBN?BSNN??KBVSN?NN??AGCWD?-"
        )

    def test_iupac_consensus_Protein(self):
        """SequenceCollection iupac_consensus should use protein IUPAC symbols correctly"""
        alignmentUpper = self.Class(
            [
                "ACDEFGHIKLMNPQRSTUVWY-",
                "ACDEFGHIKLMNPQRSUUVWF-",
                "ACDEFGHIKLMNPERSKUVWC-",
                "ACNEFGHIKLMNPQRS-UVWP-",
            ]
        )
        # following IUPAC consensus calculated by hand
        # Test all uppper
        self.assertEqual(
            alignmentUpper.iupac_consensus(PROTEIN), "ACBEFGHIKLMNPZRS?UVWX-"
        )

    def test_is_ragged(self):
        """SequenceCollection is_ragged should return true if ragged alignment"""
        assert not self.identical.is_ragged()
        assert not self.gaps.is_ragged()

    def test_probs_per_pos(self):
        """SequenceCollection.probs_per_pos should find Pr(symbol) in each
        column"""
        # make an alignment with 4 seqs (easy to calculate probabilities)
        align = self.Class(["AAA", "ACA", "GGG", "GUC"])
        got = align.probs_per_pos()
        # check that the column probs match the counts we expect
        expect = [
            {"A": 0.5, "G": 0.5},
            {"A": 0.25, "C": 0.25, "G": 0.25, "U": 0.25},
            {"A": 0.5, "G": 0.25, "C": 0.25},
        ]
        for pos, probs in enumerate(expect):
            for char, prob in probs.items():
                assert_allclose(got[pos, char], prob)

    def test_majority_consensus(self):
        """SequenceCollection.majority_consensus should return commonest symbol per column"""
        # Check the exact strings expected from string transform
        self.assertEqual(self.sequences.majority_consensus(), "UCAG")

    def test_uncertainties(self):
        """SequenceCollection.uncertainties should match hand-calculated values"""
        aln = self.Class(["ABC", "AXC"])
        obs = aln.entropy_per_pos()
        assert_allclose(obs, [0, 1, 0])
        # check what happens with only one input sequence
        aln = self.Class(["ABC"])
        obs = aln.entropy_per_pos()
        assert_allclose(obs, [0, 0, 0])

    def test_sample_with_replacement(self):
        # test with replacement -- just verify that it rnus
        alignment = self.Class({"seq1": "gatc", "seq2": "gatc"})
        sample = alignment.sample(1000, with_replacement=True)
        self.assertEqual(len(sample), 1000)
        # ensure that sampling with replacement works on single col alignment
        alignment1 = self.Class({"seq1": "A", "seq2": "A"})
        result = alignment1.sample(with_replacement=True)
        self.assertEqual(len(result), 1)

    def test_sample_tuples(self):
        ##### test with motif size != 1 #####
        alignment = self.Class(
            {
                "seq1": "AABBCCDDEEFFGGHHIIJJKKLLMMNNOOPP",
                "seq2": "AABBCCDDEEFFGGHHIIJJKKLLMMNNOOPP",
            }
        )
        shuffled = alignment.sample(motif_length=2)
        # ensure length correct
        sample = alignment.sample(10, motif_length=2)
        self.assertEqual(len(sample), 20)
        # test columns alignment preserved
        seqs = list(sample.to_dict().values())
        self.assertEqual(seqs[0], seqs[1])
        # ensure each char occurs twice as sampling dinucs without replacement
        for char in seqs[0]:
            self.assertEqual(seqs[0].count(char), 2)

    def test_copy(self):
        """correctly copy an alignment"""
        aln = self.Class(data=[("a", "AC-GT"), ("b", "ACCGT")])
        copied = aln.copy()
        self.assertTrue(type(aln), type(copied))
        self.assertEqual(aln.to_dict(), copied.to_dict())
        self.assertEqual(id(aln.moltype), id(copied.moltype))
        aln = self.Class(data=[("a", "AC-GT"), ("b", "ACCGT")], info={"check": True})
        copied = aln.copy()
        self.assertEqual(aln.info, copied.info)

    def test_to_pretty(self):
        """produce correct pretty print formatted text"""
        seqs = {"seq1": "ACGAANGA", "seq2": "-CGAACGA", "seq3": "ATGAACGA"}
        expect = ["seq1    ACGAANGA", "seq2    -....C..", "seq3    .T...C.."]

        aln = self.Class(data=seqs, moltype=DNA)
        got = aln.to_pretty(name_order=["seq1", "seq2", "seq3"])
        self.assertEqual(got, "\n".join(expect))

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
        self.assertEqual(got, "\n".join(expect))

    def test_to_html(self):
        """produce correct html formatted text"""
        seqs = {"seq1": "ACG", "seq2": "-CT"}

        aln = self.Class(data=seqs, moltype=DNA)
        got = aln.to_html(ref_name="longest")  # name_order=['seq1', 'seq2'])
        # ensure balanced tags are in the txt
        for tag in ["<style>", "</style>", "<div", "</div>", "<table>", "</table>"]:
            self.assertTrue(tag in got)

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

        self.assertTrue(ref_row in got)
        self.assertTrue(other_row in got)
        self.assertTrue(got.find(ref_row) < got.find(other_row))

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
        self.assertTrue(got.find(ref_row) < got.find(other_row))

    def test_variable_positions(self):
        """correctly identify variable positions"""
        new_seqs = {"A": "-CG-C", "B": "ACAA?", "C": "GCGAC"}
        aln = self.Class(data=new_seqs, moltype=DNA)
        self.assertEqual(aln.variable_positions(include_gap_motif=True), [0, 2, 3, 4])
        self.assertEqual(aln.variable_positions(include_gap_motif=False), [0, 2])
        new_seqs = {"A": "GCGAC", "B": "GCGAC", "C": "GCGAC"}
        aln = self.Class(data=new_seqs, moltype=DNA)
        self.assertEqual(aln.variable_positions(include_gap_motif=True), [])
        self.assertEqual(aln.variable_positions(include_gap_motif=False), [])

    def test_to_type_info(self):
        """interconverting between alignment types preserves info attribute"""
        new_seqs = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        array_align = self.Class == ArrayAlignment
        # when array_align arg matches instance class, no conversion
        # and get back self
        aln = self.Class(data=new_seqs, info={"key": "value"})
        new = aln.to_type(array_align=array_align)
        self.assertEqual(id(aln), id(new))
        self.assertEqual(new.info["key"], "value")

    def test_to_dna(self):
        """alignment cast to DNA works"""
        data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        aln = self.Class(data=data)
        dna = aln.to_dna()
        self.assertEqual(set(dna.names), set(aln.names))
        self.assertTrue(dna.moltype == DNA)
        # should fail if invalid character set
        paln = dna.get_translation()
        self.assertRaises(AlphabetError, paln.to_dna)

    def test_to_dna_info(self):
        """alignment cast to DNA preserves info attribute"""
        data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        aln = self.Class(data=data, info={"key": "value"})
        dna = aln.to_dna()
        self.assertEqual(dna.info["key"], "value")

    def test_to_rna(self):
        """alignment cast to RNA works"""
        data = {"seq1": "ACGUACGUA", "seq2": "ACCGAA---", "seq3": "ACGUACGUU"}
        aln = self.Class(data=data)
        rna = aln.to_rna()
        self.assertEqual(set(rna.names), set(aln.names))
        self.assertTrue(rna.moltype == RNA)

    def test_to_rna_info(self):
        """alignment cast to RNA preserves info attribute"""
        data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        aln = self.Class(data=data, info={"key": "value"})
        rna = aln.to_rna()
        self.assertEqual(rna.info["key"], "value")

    def test_to_protein(self):
        """alignment cast to protein works"""
        data = {"seq1": "TYV", "seq3": "TYV", "seq2": "TE-"}
        aln = self.Class(data=data)
        paln = aln.to_protein()
        self.assertEqual(set(paln.names), set(aln.names))
        self.assertTrue(paln.moltype == PROTEIN)
        # should fail if invalid character set
        self.assertRaises(AlphabetError, paln.to_dna)

    def test_to_protein_info(self):
        """alignment cast to protein preserves info attribute"""
        data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        aln = self.Class(data=data, info={"key": "value"})
        dna = aln.to_dna()
        self.assertEqual(dna.info["key"], "value")

    def test_replace_seqs_info(self):
        """replace_seqs should preserve info attribute"""
        a = self.Class(
            {"seq1": "ACGU", "seq2": "C-UA", "seq3": "C---"}, info={"key": "value"}
        )
        seqs = {"seq1": "AAACCCGGGUUU", "seq2": "CCCUUUAAA", "seq3": "CCC"}
        result = a.replace_seqs(seqs)  # default behaviour
        self.assertEqual(result.info["key"], "value")

    def test_counts(self):  # ported
        """SequenceCollection.counts handles motif length, allow_gaps etc.."""
        data = {"a": "AAAA??????", "b": "CCCGGG--NN"}
        coll = self.Class(data=data, moltype=DNA)
        got = coll.counts()
        expect = dict(A=4, C=3, G=3)
        for k, v in expect.items():
            self.assertEqual(got[k], v)

        got = coll.counts(motif_length=2)
        expect = dict(AA=2, CC=1, CG=1, GG=1)
        for k, v in expect.items():
            self.assertEqual(got[k], v)

        got = coll.counts(motif_length=2, allow_gap=True)
        expect.update({"--": 1})
        for k, v in expect.items():
            self.assertEqual(got[k], v)

        got = coll.counts(motif_length=2, include_ambiguity=True, allow_gap=True)
        expect = dict(AA=2, CC=1, CG=1, GG=1, NN=1)
        expect.update({"??": 3, "--": 1})
        for k, v in expect.items():
            self.assertEqual(got[k], v)

    def test_counts_per_pos(self):
        """correctly count motifs"""
        exp = array(
            [
                [1, 1, 1, 0],
                [0, 2, 0, 1],
                [0, 0, 3, 0],
                [1, 1, 0, 1],
                [0, 0, 3, 0],
                [1, 1, 0, 1],
            ]
        )

        exp_gap = array(
            [
                [1, 1, 0, 1, 0],
                [0, 2, 0, 0, 1],
                [0, 0, 3, 0, 0],
                [0, 2, 0, 1, 0],
                [0, 1, 2, 0, 0],
                [0, 2, 0, 1, 0],
            ]
        )

        s1 = DNA.make_seq(seq="TCAGAG", name="s1")
        s2 = DNA.make_seq(seq="CCACAC", name="s2")
        s3 = DNA.make_seq(seq="AGATAT", name="s3")
        s4 = DNA.make_seq(seq="G-ACCC", name="s4")
        aln = self.Class([s1, s2, s3], moltype=DNA)
        obs = aln.counts_per_pos()
        assert_equal(obs.array, exp)
        assert_equal(obs.motifs, tuple(DNA.alphabet))
        obs = aln.counts_per_pos(motif_length=2)
        assert_equal(obs[0, "TC"], 1)
        assert_equal(obs[1, "AC"], 1)
        assert_equal(obs[2, "AC"], 1)
        aln = self.Class([s1, s2, s4], moltype=DNA)
        obs = aln.counts_per_pos(allow_gap=True)
        assert_equal(obs.array, exp_gap)
        aln = self.Class(["-RAT", "ACCT", "GTGT"], moltype="dna")
        c = aln.counts_per_pos(include_ambiguity=False, allow_gap=True)
        assert_equal(set(c.motifs), set("ACGT-"))

    def test_counts_per_seq_default_moltype(self):
        """produce correct counts per seq with default moltypes"""
        data = {"a": "AAAA??????", "b": "CCCGGG--NN", "c": "CCGGTTCCAA"}
        coll = self.Class(data=data)
        got = coll.counts_per_seq()
        try:
            self.assertEqual(got.col_sum()["-"], 0)
        except KeyError:
            pass  # text moltype in Alignment excludes '-'
        got = coll.counts_per_seq(include_ambiguity=True, allow_gap=True)
        self.assertEqual(got.col_sum()["-"], 2)

    def test_counts_per_pos_default_moltype(self):  # will not port
        """produce correct counts per pos with default moltypes"""
        data = {"a": "AAAA??????", "b": "CCCGGG--NN", "c": "CCGGTTCCAA"}
        coll = self.Class(data=data)
        got = coll.counts_per_pos()
        # should not include gap character
        self.assertNotIn("-", got.motifs)
        # allowing gaps
        got = coll.counts_per_pos(allow_gap=True)
        # should include gap character
        self.assertEqual(got[5, "-"], 0)
        self.assertEqual(got[6, "-"], 1)

        # now with motif-length 2
        got = coll.counts_per_pos(motif_length=2)
        found_motifs = set()
        lengths = set()
        for m in got.motifs:
            lengths.add(len(m))
            found_motifs.update(m)
        self.assertTrue("-" not in found_motifs)
        self.assertEqual(lengths, {2})

    def test_entropy_per_pos_just_gaps(self):
        """pos with just gaps have nan"""
        a = self.Class(dict(a="A---", b="C---", c="C---"), moltype=DNA)
        entropy = a.entropy_per_pos()
        assert_allclose(entropy, [0.91829583, numpy.nan, numpy.nan, numpy.nan])
        a = self.Class(dict(a="---", b="---", c="---"), moltype=DNA)
        entropy = a.entropy_per_pos()
        assert_allclose(entropy, [numpy.nan, numpy.nan, numpy.nan])

    def test_entropy_excluding_unobserved(self):
        """omitting unobserved motifs should not affect entropy calculation"""
        a = self.Class(dict(a="ACAGGG", b="AGACCC", c="GGCCTA"), moltype=DNA)
        entropy_excluded = a.entropy_per_seq(exclude_unobserved=True)
        entropy_unexcluded = a.entropy_per_seq(exclude_unobserved=False)
        assert_allclose(entropy_excluded, entropy_unexcluded)

    def test_distance_matrix(self):
        """Alignment distance_matrix should produce correct scores"""
        data = dict([("s1", "ACGTACGTA"), ("s2", "GTGTACGTA")])
        aln = self.Class(data=data, moltype="dna")
        dists = aln.distance_matrix(calc="hamming", show_progress=False)
        self.assertEqual(dists, {("s1", "s2"): 2.0, ("s2", "s1"): 2.0})
        # and for protein
        aa = aln.get_translation()
        dists = aa.distance_matrix(calc="hamming")
        self.assertEqual(dists, {("s1", "s2"): 1.0, ("s2", "s1"): 1.0})

        # when there are invalid data
        data = dict(
            seq1="AGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
            seq2="TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
            seq3="TACAAAAAAAAGGGGCGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
        )

        aln = self.Class(data=data, moltype="dna")
        with self.assertRaises(ArithmeticError):
            # default settings cause an exception
            dists = aln.distance_matrix(calc="paralinear")
        # but setting drop_invalid=False allows calc
        dists = aln.distance_matrix(calc="paralinear", drop_invalid=True)

    def test_get_gapped_seq(self):
        """Alignment.get_gapped_seq should return seq, with gaps"""
        aln = self.Class({"seq1": "--TTT?", "seq2": "GATC??"})
        self.assertEqual(str(aln.get_gapped_seq("seq1")), "--TTT?")

    def test_count_gaps_per_seq(self):
        """correctly compute the number of gaps"""
        data = {"a": "AAAA---GGT", "b": "CCC--GG?GT"}
        aln = self.Class(data=data, moltype=DNA)
        got = aln.count_gaps_per_seq(include_ambiguity=False)
        assert_equal(got.array, [3, 2])
        assert_equal(got["b"], 2)
        got = aln.count_gaps_per_seq(include_ambiguity=True)
        assert_equal(got.array, [3, 3])
        assert_equal(got["b"], 3)
        # per seq, unique
        got = aln.count_gaps_per_seq(include_ambiguity=False, unique=True)
        assert_equal(got.array, [1, 2])
        got = aln.count_gaps_per_seq(include_ambiguity=True, unique=True)
        assert_equal(got.array, [2, 2])

        data = {"a": "AAAGGG", "b": "------", "c": "------"}
        aln = self.Class(data=data, moltype=DNA)
        got = aln.count_gaps_per_seq(include_ambiguity=False, unique=True)
        assert_equal(got.array, [6, 0, 0])
        assert_equal(got["a"], 6)
        assert_equal(got["b"], 0)

        # per_seq, induced_by
        data = {"a": "--ACGT---GTAC", "b": "--ACGTA--GT--", "c": "--ACGTA-AGT--"}
        aln = self.Class(data=data, moltype=DNA)
        got = aln.count_gaps_per_seq(unique=False, induced_by=True)
        assert_equal(got.array, [2, 1, 2])
        assert_equal(got["b"], 1)

    def test_info_source(self):
        """info.source exists if load_aligned_seqs given a filename"""
        array_align = self.Class == ArrayAlignment
        seqs = load_aligned_seqs("data/brca1.fasta", array_align=array_align)
        self.assertEqual(seqs.info.source, "data/brca1.fasta")

    def test_seq_entropy_just_gaps(self):
        """get_seq_entropy should get entropy of each seq"""
        a = self.Class(dict(a="A---", b="----"), moltype=DNA)
        entropy = a.entropy_per_seq()
        assert_allclose(entropy, [0, numpy.nan])
        a = self.Class(dict(a="----", b="----"), moltype=DNA)
        entropy = a.entropy_per_seq()
        self.assertIs(entropy, None)

    def test_repr_html(self):  # ported
        """exercises method normally invoked in notebooks"""
        aln = self.Class({"a": "AAAAA", "b": "AAA--"})
        aln.set_repr_policy(num_seqs=5, num_pos=40)
        self.assertEqual(aln[:3]._repr_policy, aln._repr_policy)
        row_a = '<tr><td class="label">a</td>'
        row_b = '<tr><td class="label">b</td>'
        # default order is longest sequence at top
        got = aln._repr_html_()
        self.assertTrue(got.find(row_a) < got.find(row_b))
        # change order, a should now be last
        aln.set_repr_policy(num_seqs=5, num_pos=40, ref_name="b")
        got = aln._repr_html_()
        self.assertTrue(got.find(row_a) > got.find(row_b))
        # tests repr policy has been successfully applied
        aln = load_aligned_seqs("data/brca1.fasta", moltype="dna")
        aln.set_repr_policy(num_seqs=2)
        got = aln._repr_html_()
        self.assertEqual(got.count("</tr>"), 3)
        aln.set_repr_policy(num_seqs=3)
        got = aln._repr_html_()
        self.assertEqual(got.count("</tr>"), 4)
        aln.set_repr_policy(num_seqs=len(aln.seqs))
        got = aln._repr_html_()
        self.assertEqual(got.count("</tr>"), len(aln.seqs) + 1)
        # tests _repr_html_ displays correct number of sequences
        aln = load_aligned_seqs("data/brca1.fasta", moltype="dna")
        got = aln._repr_html_()
        self.assertIn("%d x %d" % (aln.num_seqs, aln.seq_len), got.splitlines()[-2])

    def test_seqlogo(self):
        """exercise producing a seq logo"""
        data = {
            "seq1": "CAGGTCGACCTCGGC---------CACGAC",
            "seq2": "CAGATCGACCTCGGC---------CACGAC",
            "seq3": "CAGATCGACCTCGGT---------CACGAT",
            "seq4": "CAGATCGACCTCGGCGAACACGGCCATGAT",
            "seq5": "CCGATCGACATGGGC---------CACGAT",
            "seq6": "GCC---------------------------",
        }
        # with a defined moltype
        aln = self.Class(data, moltype=DNA)
        logo = aln.seqlogo()
        # using wrap argument
        logo = aln.seqlogo(wrap=20)
        # should work for protein too
        aa = aln.get_translation()
        aa.seqlogo()

        # without a defined moltype
        aln = self.Class(data)
        aln.seqlogo()


class ArrayAlignmentTests(AlignmentBaseTests, TestCase):
    Class = ArrayAlignment

    def test_slice_align(self):
        """slicing alignment should work correctly"""
        data = {"seq1": "ACGACGACG", "seq2": "ACGACGACG", "seq3": "ACGACGACG"}
        alignment = self.Class(data=data)
        sub_align = alignment[2:5]
        self.assertTrue(isinstance(sub_align, self.Class))
        expect = {"seq1": "GAC", "seq2": "GAC", "seq3": "GAC"}
        self.assertEqual(sub_align.to_dict(), expect)
        # slice third positions
        sub_align = alignment[2::3]
        expect = {"seq1": "GGG", "seq2": "GGG", "seq3": "GGG"}
        self.assertEqual(sub_align.to_dict(), expect)

    def test_slice_align_info(self):
        """slicing alignment preserves info attribute"""
        data = {"seq1": "ACGACGACG", "seq2": "ACGACGACG", "seq3": "ACGACGACG"}
        alignment = self.Class(data=data, info={"key": "value"})
        sub_align = alignment[2:5]
        self.assertTrue(len(sub_align) == 3)
        self.assertEqual(sub_align.info["key"], "value")


class AlignmentTests(AlignmentBaseTests, TestCase):
    Class = Alignment

    def test_sliding_windows(self):
        """sliding_windows should return slices of alignments."""
        alignment = self.Class(
            {"seq1": "ACGTACGT", "seq2": "ACGTACGT", "seq3": "ACGTACGT"}
        )
        result = []
        for bit in alignment.sliding_windows(5, 2):
            result += [bit]
        self.assertEqual(
            result[0].to_dict(), {"seq3": "ACGTA", "seq2": "ACGTA", "seq1": "ACGTA"}
        )
        self.assertEqual(
            result[1].to_dict(), {"seq3": "GTACG", "seq2": "GTACG", "seq1": "GTACG"}
        )

        result = []
        for bit in alignment.sliding_windows(5, 1):
            result += [bit]
        self.assertEqual(
            result[0].to_dict(), {"seq3": "ACGTA", "seq2": "ACGTA", "seq1": "ACGTA"}
        )
        self.assertEqual(
            result[1].to_dict(), {"seq3": "CGTAC", "seq2": "CGTAC", "seq1": "CGTAC"}
        )
        self.assertEqual(
            result[2].to_dict(), {"seq3": "GTACG", "seq2": "GTACG", "seq1": "GTACG"}
        )
        self.assertEqual(
            result[3].to_dict(), {"seq3": "TACGT", "seq2": "TACGT", "seq1": "TACGT"}
        )

    def test_with_gaps_from(self):
        """with_gaps_from should overwrite with gaps."""
        gapless = self.Class({"seq1": "TCG", "seq2": "TCG"})
        pregapped = self.Class({"seq1": "-CG", "seq2": "TCG"})
        template = self.Class({"seq1": "A-?", "seq2": "ACG"})
        r1 = gapless.with_gaps_from(template).to_dict()
        r2 = pregapped.with_gaps_from(template).to_dict()
        self.assertEqual(r1, {"seq1": "T-G", "seq2": "TCG"})
        self.assertEqual(r2, {"seq1": "--G", "seq2": "TCG"})

    def test_get_degapped_relative_to_info(self):
        """should remove all columns with a gap in sequence with given name
        while preserving info attribute"""
        aln = self.Class(
            [
                ["name1", "-AC-DEFGHI---"],
                ["name2", "XXXXXX--XXXXX"],
                ["name3", "YYYY-YYYYYYYY"],
                ["name4", "-KL---MNPR---"],
            ],
            info={"key": "foo"},
        )
        out_aln = self.Class(
            [
                ["name1", "ACDEFGHI"],
                ["name2", "XXXX--XX"],
                ["name3", "YY-YYYYY"],
                ["name4", "KL--MNPR"],
            ],
            info={"key": "bar"},
        )
        gdrt = aln.get_degapped_relative_to("name1")
        self.assertEqual(gdrt.info["key"], "foo")

    def test_deepcopy2(self):
        """ "Aligned.deepcopy correctly handles gapped sequences"""
        seqs = self.Class(
            data={
                "a": "CAGATTTGGCAGTT-",
                "b": "-AGATTCAGCAGGTG",
                "c": "CAGAT-CAGCAGGTG",
                "d": "CAGATTCAGCAGGTG",
            },
            moltype="dna",
        )
        lengths = {len(s.deepcopy(sliced=True)) for s in seqs.seqs}
        self.assertEqual(lengths, {len(seqs)})
        rc = seqs.rc()
        lengths = {len(s.deepcopy(sliced=True)) for s in rc.seqs}
        self.assertEqual(lengths, {len(seqs)})

    def test_to_moltype_annotations(self):
        """correctly convert to specified moltype with proper sequence annotations"""

        s1 = Sequence(seq="TTTTTTAAAA", name="test_seq1")
        s2 = Sequence(seq="AAAATTTTTT", name="test_seq2")
        s3 = Sequence(seq="AATTTTTAAA", name="test_seq3")
        s1.add_feature(biotype="exon", name="fred", spans=[(0, 6)])
        s2.add_feature(biotype="exon", name="fred", spans=[(4, 10)])
        s3.add_feature(biotype="exon", name="fred", spans=[(2, 7)])
        data = {"seq1": s1, "seq2": s2, "seq3": s3}
        aln = self.Class(data=data)
        aln.add_feature(biotype="demo", name="one", spans=[(0, 1), (2, 4)])
        rna = aln.to_moltype("rna")
        for name in rna.names:
            orig_seq = aln.get_seq(name)
            new_seq = rna.get_seq(name)
        # check the sequence moltypes
        self.assertEqual({s.data.moltype.label for s in rna.seqs}, {"rna"})
        self.assertEqual(rna.moltype.label, "rna")

    def test_construction(self):
        """correctly construct from list of sequences of length 2"""
        seq1 = make_seq(seq="AC-", name="seq1")
        seq1 = Aligned(*seq1.parse_out_gaps())
        seq2 = make_seq(seq="ACG", name="seq2")
        seq2 = Aligned(*seq2.parse_out_gaps())
        coll = self.Class(data=[seq1, seq2])


class ArrayAlignmentSpecificTests(TestCase):
    """Tests of the ArrayAlignment object and its methods"""

    def setUp(self):
        """Define some standard alignments."""
        self.a2 = ArrayAlignment(["ABC", "DEF"], names=["x", "y"])
        seqs = []
        for s in ["abaa", "abbb"]:
            seqs.append(AB.make_seq(seq=s, preserve_case=True))
        self.a = ArrayAlignment(seqs, moltype=AB)
        self.b = Alignment(["ABC", "DEF"])
        self.c = SequenceCollection(["ABC", "DEF"])

    def test_init(self):
        """ArrayAlignment init should work from a sequence"""
        a = ArrayAlignment(array([[0, 1, 2], [3, 4, 5]]).T)
        assert_equal(a.seq_data, array([[0, 3], [1, 4], [2, 5]], "B"))
        assert_equal(a.array_positions, array([[0, 1, 2], [3, 4, 5]], "B"))
        assert_equal(a.names, ["seq_0", "seq_1", "seq_2"])

    def test_init_seqs(self):
        """ArrayAlignment init should work from ArraySequence objects."""
        s = list(map(ArraySequence, ["abc", "def"]))
        a = ArrayAlignment(s)
        assert_equal(a.seq_data, array(["abc", "def"], "c").view("B"))

    def test_init_generic(self):
        """ArrayAlignment init should work from generic objects."""
        s = ["abc", "def"]
        a = ArrayAlignment(s)
        assert_equal(a.seq_data, array(["abc", "def"], "c").view("B"))

    def test_init_aln(self):
        """ArrayAlignment init should work from another alignment."""
        s = ["abc", "def"]
        a = ArrayAlignment(s)
        b = ArrayAlignment(a)
        self.assertIsNot(a.seq_data, b.seq_data)
        assert_equal(b.seq_data, array(["abc", "def"], "c").view("B"))

    def test_init_dict(self):
        """ArrayAlignment init should work from dict."""
        s = {"abc": "AAACCC", "xyz": "GCGCGC"}
        a = ArrayAlignment(s, names=["abc", "xyz"])
        assert_equal(a.seq_data, array(["AAACCC", "GCGCGC"], "c").view("B"))
        self.assertEqual(tuple(a.names), ("abc", "xyz"))

    def test_iter(self):
        """ArrayAlignment iter should iterate over positions"""
        result = list(iter(self.a2))
        for i, j in zip(result, ["AD", "BE", "CF"]):
            self.assertEqual(i, list(j))

    def test_getitem(self):
        """ArrayAlignment getitem act like standard alignment slice"""
        a2 = self.a2
        expect = {"x": "B", "y": "E"}
        got = a2[1]
        self.assertEqual(got.to_dict(), expect)
        expect = {"x": "BC", "y": "EF"}
        got = a2[1:]
        self.assertEqual(got.to_dict(), expect)

    def test_get_sub_alignment(self):
        """ArrayAlignment get_sub_alignment should get requested part of alignment"""
        a = ArrayAlignment({"x": "ABCE", "y": "FGHI", "z": "JKLM"})
        # passing in positions should keep all seqs, but just selected
        # positions
        b = ArrayAlignment({"x": "BC", "y": "GH", "z": "KL"})
        a_1 = a.get_sub_alignment(pos=[1, 2])
        self.assertEqual(a_1.names, b.names)

        self.assertEqual(a_1.seqs, b.seqs)
        # ...and with negate_pos, should keep all except the positions passed in
        a_2 = a.get_sub_alignment(pos=[0, 3], negate_pos=True)
        self.assertEqual(a_2.seqs, b.seqs)
        self.assertEqual(a_2.names, b.names)
        # passing in seqs should keep all positions, but just selected seqs
        c = ArrayAlignment({"x": "ABCE", "z": "JKLM"})
        a_3 = a.get_sub_alignment(seqs=[0, 2])
        self.assertEqual(a_3.seqs, c.seqs)
        # check that labels were updates as well...
        self.assertEqual(a_3.names, c.names)
        # ...and should work with negate_seqs to exclude just selected seqs
        a_4 = a.get_sub_alignment(seqs=[1], negate_seqs=True)
        self.assertEqual(a_4.seqs, c.seqs)
        self.assertEqual(a_4.names, c.names)
        # should be able to do both seqs and positions simultaneously
        d = ArrayAlignment({"x": "BC", "z": "KL"})
        a_5 = a.get_sub_alignment(seqs=[0, 2], pos=[1, 2])
        self.assertEqual(a_5.seqs, d.seqs)
        self.assertEqual(a_5.names, d.names)

    def test_get_sub_alignment_info(self):
        """ArrayAlignment get_sub_alignment should preserve info attribute"""
        a = ArrayAlignment({"x": "ABCE", "y": "FGHI", "z": "JKLM"}, info={"key": "foo"})
        # passing in positions should keep all seqs, but just selected
        # positions
        b = ArrayAlignment({"x": "BC", "y": "GH", "z": "KL"}, info={"key": "bar"})
        a_1 = a.get_sub_alignment(pos=[1, 2])
        self.assertEqual(a_1.names, b.names)
        self.assertEqual(a.info["key"], "foo")

    def test_str(self):
        """ArrayAlignment str should return FASTA representation of aln"""
        self.assertEqual(str(self.a2), ">x\nABC\n>y\nDEF\n")
        # should work if labels diff length
        self.a2.names[-1] = "yyy"
        self.assertEqual(str(self.a2), ">x\nABC\n>yyy\nDEF\n")

    def test_entropy_per_pos(self):
        """entropy_per_pos should get entropy of each pos"""
        a = self.a
        f = a.entropy_per_pos()
        e = array([0, 0, 1, 1])
        assert_allclose(f, e)
        f = a.entropy_per_pos(motif_length=2)
        e = array([0, 1])
        assert_allclose(f, e)
        seqs = []
        for s in ["-GAT", "ACCT", "GAGT"]:
            seqs.append(make_seq(seq=s, moltype="dna"))
        a = ArrayAlignment(seqs)
        f = a.entropy_per_pos(allow_gap=True)
        e = array([1.584962500721156, 1.584962500721156, 1.584962500721156, 0])
        assert_allclose(f, e)

        seqs = []
        for s in ["-RAT", "ACCT", "GTGT"]:
            seqs.append(make_seq(seq=s, moltype="dna"))
        a = ArrayAlignment(seqs)

        # "-RAT"
        # "ACCT"
        # "GTGT"
        f = a.entropy_per_pos(allow_gap=False, include_ambiguity=False)
        e = [
            2 * safe_p_log_p(array([1 / 2])).sum(),
            2 * safe_p_log_p(array([1 / 2])).sum(),
            3 * safe_p_log_p(array([1 / 3])).sum(),
            0,
        ]
        assert_allclose(f, e)

        f = a.entropy_per_pos(include_ambiguity=True)
        e = [
            2 * safe_p_log_p(array([1 / 2])).sum(),
            3 * safe_p_log_p(array([1 / 3])).sum(),
            3 * safe_p_log_p(array([1 / 3])).sum(),
            0,
        ]
        assert_allclose(f, e)

        f = a.entropy_per_pos(allow_gap=True)
        e = [
            3 * safe_p_log_p(array([1 / 3])).sum(),
            2 * safe_p_log_p(array([1 / 2])).sum(),
            3 * safe_p_log_p(array([1 / 3])).sum(),
            0,
        ]
        assert_allclose(f, e)


class IntegrationTests(TestCase):
    """Test for integration between regular and model seqs and alns"""

    def setUp(self):
        """Intialize some standard sequences"""
        self.r1 = RNA.make_seq(seq="AAA", name="x")
        self.r2 = RNA.make_seq(seq="CCC", name="y")
        self.m1 = RNA.make_array_seq(seq="AAA", name="xx")
        self.m2 = RNA.make_array_seq(seq="CCC", name="yy")

    def test_model_to_model(self):
        """Model seq should work with dense alignment"""
        a = ArrayAlignment([self.m1, self.m2])
        self.assertEqual(a.to_dict(), {"xx": "AAA", "yy": "CCC"})
        a = ArrayAlignment([self.m1, self.m2], moltype=DNA)
        self.assertEqual(a.to_dict(), {"xx": "AAA", "yy": "CCC"})
        self.assertEqual(self.m1.name, "xx")

    def test_regular_to_model(self):
        """Regular seq should work with dense alignment"""
        a = ArrayAlignment([self.r1, self.r2])
        self.assertEqual(str(a), ">x\nAAA\n>y\nCCC\n")
        a = ArrayAlignment([self.r1, self.r2], moltype=DNA)
        self.assertEqual(str(a), ">x\nAAA\n>y\nCCC\n")
        self.assertEqual(self.r1.name, "x")

    def test_model_to_regular(self):
        """Model seq should work with regular alignment"""
        a = Alignment([self.m1, self.m2])
        self.assertEqual(str(a), ">xx\nAAA\n>yy\nCCC\n")
        a = Alignment([self.m1, self.m2], moltype=DNA)
        self.assertEqual(str(a), ">xx\nAAA\n>yy\nCCC\n")
        self.assertEqual(self.m1.name, "xx")

    def test_regular_to_regular(self):
        """Regular seq should work with regular alignment"""
        a = Alignment([self.r1, self.r2])
        self.assertEqual(str(a), ">x\nAAA\n>y\nCCC\n")
        a = Alignment([self.r1, self.r2], moltype=DNA)
        self.assertEqual(str(a), ">x\nAAA\n>y\nCCC\n")
        self.assertEqual(self.r1.name, "x")

    def test_model_aln_to_regular_aln(self):
        """Dense aln should convert to regular aln"""
        a = ArrayAlignment([self.r1, self.r2])
        d = Alignment(a)
        self.assertEqual(str(d), ">x\nAAA\n>y\nCCC\n")
        d = Alignment(a, moltype=DNA)
        self.assertEqual(str(d), ">x\nAAA\n>y\nCCC\n")
        self.assertEqual(self.r1.name, "x")

    def test_regular_aln_to_model_aln(self):
        """Regular aln should convert to model aln"""
        a = Alignment([self.r1, self.r2])
        d = ArrayAlignment(a)
        self.assertEqual(str(d), ">x\nAAA\n>y\nCCC\n")
        d = ArrayAlignment(a, moltype=DNA)
        self.assertEqual(str(d), ">x\nAAA\n>y\nCCC\n")
        self.assertEqual(self.r1.name, "x")

    def test_regular_aln_to_regular_aln(self):
        """Regular aln should convert to regular aln"""
        a = Alignment([self.r1, self.r2])
        d = Alignment(a)
        self.assertEqual(str(d), ">x\nAAA\n>y\nCCC\n")
        d = Alignment(a, moltype=DNA)
        self.assertEqual(str(d), ">x\nAAA\n>y\nCCC\n")
        self.assertEqual(self.r1.name, "x")

    def test_model_aln_to_model_aln(self):
        """Model aln should convert to model aln"""
        a = Alignment([self.r1, self.r2])
        d = Alignment(a)
        self.assertEqual(str(d), ">x\nAAA\n>y\nCCC\n")
        d = Alignment(a, moltype=DNA)
        self.assertEqual(str(d), ">x\nAAA\n>y\nCCC\n")
        self.assertEqual(self.r1.name, "x")


@pytest.mark.parametrize(
    "raw_seq,coords",
    (("ACGGTAAAG", ((2, 4), (5, 8))), ("CCC---CCC", ((0, 3), (6, 9)))),
)
def test_featuremap_slice_aligned(raw_seq, coords):
    from cogent3.core.location import FeatureMap, Span

    im, seq = DNA.make_seq(seq=raw_seq).parse_out_gaps()
    ia = Aligned(im, seq)
    length = len(raw_seq)
    fmap = FeatureMap(spans=[Span(s, e) for s, e in coords], parent_length=length)
    expect = "".join(raw_seq[s:e] for s, e in fmap.get_coordinates())
    got = ia[fmap]
    assert str(got) == expect


@pytest.mark.parametrize("cls", (ArrayAlignment, Alignment))
def test_to_type(cls):
    """correctly interconvert between alignment types"""
    new_seqs = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    array_align = cls == ArrayAlignment
    # when array_align arg matches instance class, no conversion
    # and get back self
    aln = cls(data=new_seqs)
    new = aln.to_type(array_align=array_align)
    assert id(aln) == id(new)

    # when array_align arg does not match, should get back the opposite type
    new = aln.to_type(array_align=not array_align)
    assert not isinstance(new, cls)

    # we should be able to specify moltype and alignment
    new = aln.to_type(array_align=not array_align, moltype=DNA)
    assert new.to_dict() == new_seqs
    # and translate
    assert new.get_translation().to_dict() == {
        "seq1": "TYV",
        "seq3": "TYV",
        "seq2": "TE-",
    }

    # should work on ArrayAlign when just moltype changes
    new_seqs = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---"}
    aln = cls(data=new_seqs)
    new = aln.to_type(array_align=array_align, moltype=DNA)
    new = new.no_degenerates()  # this should not fail!
    assert len(new) == len(aln) - 3

    # should correctly apply to existing moltype
    aln = cls(data=new_seqs, moltype=DNA)
    new = aln.to_type(array_align=not array_align)
    assert aln.moltype == new.moltype


@pytest.mark.parametrize(
    "moltype",
    ["rna", "dna", "protein"],
)
def test_upac_consensus_allow_gaps(moltype):
    aln = make_aligned_seqs(
        data={"s1": "ACGG", "s2": "ACGG", "s3": "-CGG"},
        moltype=moltype,
        array_align=False,
    )
    # default behaviour
    iupac = aln.iupac_consensus()
    assert iupac == "?CGG"

    # allow_gaps
    iupac = aln.iupac_consensus(allow_gap=False)
    assert iupac == "ACGG"


@pytest.mark.parametrize(
    "moltype",
    ["rna", "dna", "protein"],
)
def test_upac_consensus_allow_gaps_array_alignment(moltype):
    aln = make_aligned_seqs(
        data={"s1": "ACGG", "s2": "ACGG", "s3": "-CGG"},
        moltype=moltype,
        array_align=False,
    )
    # default behaviour
    iupac = aln.iupac_consensus()
    assert iupac == "?CGG"

    # allow_gaps
    iupac = aln.iupac_consensus(allow_gap=False)
    assert iupac == "ACGG"


def test_rc_iter():
    dna = DNA.make_seq(seq="ACG", name="seq1")
    rc = dna.rc()

    got = [x for x in rc]
    expect = ["C", "G", "T"]
    assert got == expect


def test_to_dna_raises():
    seq = make_seq(seq="ETV", moltype="protein")
    with pytest.raises(AlphabetError):
        seq.to_moltype("dna")


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment))
def test_annotate_from_gff3(cls):  # will not port
    """annotate_from_gff should work on data from gff3 files"""
    fasta_path = os.path.join("data/c_elegans_WS199_dna_shortened.fasta")
    gff3_path = os.path.join("data/c_elegans_WS199_shortened_gff.gff3")
    seq = load_seq(fasta_path, moltype="dna")

    # using annotate_from_gff will nest annotations
    seq_coll = cls({seq.name: seq})
    seq_coll.annotate_from_gff(gff3_path)
    member_seq = seq_coll.get_seq(seq.name)

    matches = list(member_seq.get_features())
    # 11 features
    assert len(matches) == 11
    matches = list(member_seq.get_features(biotype="gene"))
    assert len(matches) == 1
    matches = list(matches[0].get_children(biotype="mRNA"))
    assert len(matches) == 1
    matches = list(matches[0].get_children(biotype="exon"))
    assert len(matches) == 3


@pytest.fixture(scope="session")
def seqcoll_db():
    fasta_path = os.path.join("data/c_elegans_WS199_dna_shortened.fasta")
    gff3_path = os.path.join("data/c_elegans_WS199_shortened_gff.gff3")
    seq = load_seq(fasta_path, moltype="dna")
    seq_coll = SequenceCollection({seq.name: seq})
    seq_coll.annotate_from_gff(gff3_path)
    return seq_coll


def test_seqcoll_query(seqcoll_db):
    """querying from a SequenceCollection produces features bound to their seqs"""
    matches = list(seqcoll_db.get_features(seqid="I"))
    # 11 features
    assert len(matches) == 11
    matches = list(seqcoll_db.get_features(seqid="I", biotype="gene"))
    assert len(matches) == 1
    matches = list(matches[0].get_children(biotype="mRNA"))
    assert len(matches) == 1
    matches = list(matches[0].get_children(biotype="exon"))
    assert len(matches) == 3


def test_align_get_features():
    #                    0123456789   the positions
    seq1 = make_seq(seq="ACG--ACCGT", moltype="dna", name="seq1")
    seq2 = make_seq(seq="ACGGGCCCGT", moltype="dna", name="seq2")
    #                      *****      the CDS feature
    seq2.add_feature(biotype="CDS", name="fake01", spans=[(2, 7)], strand="+")
    aln = make_aligned_seqs(data=[seq1, seq2], array_align=False)
    feat = list(aln.get_features(biotype="CDS"))[0]
    sl = aln[feat]
    # slice is correct length
    assert len(sl) == (7 - 2)
    # returns correct value
    assert sl.to_dict() == dict(seq1="G--AC", seq2="GGGCC")


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment))
def test_init_annotated_seqs(cls):  # ported for SequenceCollection
    """correctly construct from list with annotated seq"""
    seq = make_seq(seq="GCCAGGGGGGAAAG-GGAGAA", name="seq1")
    _ = seq.add_feature(biotype="exon", name="name", spans=[(4, 10)])
    coll = cls(data=[seq])
    features = list(coll.get_features(biotype="exon"))
    assert len(features) == 1


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment))
def test_get_annotations_from_any_seq(cls):  # ported for SequenceCollection
    """get_annotations_from_any_seq returns correct annotations"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = cls(data, moltype=DNA)
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


def test_annotate_matches_to():
    """Aligned.annotate_matches_to correctly delegates to sequence"""

    aln = Alignment(dict(x="TTCCACTTCCGCTT"), moltype="dna")
    aln.annotation_db = GffAnnotationDb()
    seq = aln.named_seqs["x"]
    pattern = "CCRC"
    annot = seq.annotate_matches_to(
        pattern=pattern, biotype="domain", name="fred", allow_multiple=True
    )
    got = [a.get_slice() for a in annot]
    matches = ["CCAC", "CCGC"]
    assert got == matches
    annot = seq.annotate_matches_to(
        pattern=pattern, biotype="domain", name="fred", allow_multiple=False
    )
    got = [a.get_slice() for a in annot]
    assert got == matches[:1]

    # handles regex from aa
    aln = Alignment(dict(x="TTCCACTTCCGCTT"), moltype="dna")
    aln.annotation_db = GffAnnotationDb()
    gc = get_code(1)
    aa_regex = gc.to_regex("FHF")
    s = aln.named_seqs["x"].annotate_matches_to(
        aa_regex, "domain", "test", allow_multiple=False
    )
    a = list(aln.get_features(seqid="x"))[0]
    assert str(aln[a].named_seqs["x"]) == "TTCCACTTC"


@pytest.fixture(scope="function")
def gb_db(DATA_DIR):
    from cogent3.core.annotation_db import load_annotations

    return load_annotations(path=DATA_DIR / "annotated_seq.gb")


@pytest.fixture(scope="function")
def gff_db(DATA_DIR):
    from cogent3.core.annotation_db import load_annotations

    return load_annotations(path=DATA_DIR / "simple.gff")


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment))
def test_copy_annotations(cls, gff_db):  # ported for SequenceCollection
    """copy_annotations copies records from annotation db"""

    seq_coll = cls({"seq1": "ACGU", "seq2": "CGUA", "test_seq": "CCGU"})
    seq_coll.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])
    seq_coll.add_feature(seqid="seq2", biotype="xyzzz", name="abc", spans=[(1, 2)])
    expect = seq_coll.annotation_db.num_matches() + gff_db.num_matches()
    seq_coll.copy_annotations(gff_db)
    assert seq_coll.annotation_db.num_matches() == expect


def _make_seq(name):
    raw_seq = "AACCCAAAATTTTTTGGGGGGGGGGCCCC"
    cds = (15, 25)
    utr = (12, 15)
    # name is required for creating annotations
    seq = DNA.make_seq(seq=raw_seq, name=name)
    seq.add_feature(biotype="CDS", name="CDS", spans=[cds])
    seq.add_feature(biotype="5'UTR", name="5' UTR", spans=[utr])
    return seq


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment))
def test_init_seqs_have_annotations(cls, gff_db):  # ported for SequenceCollection
    """annotations on input seqs correctly merged and propagated"""

    seq_coll = cls({"seq1": _make_seq("seq1"), "seq2": _make_seq("seq2")})
    coll_db = seq_coll.annotation_db
    assert len(coll_db) == 4
    for seq in seq_coll.seqs:
        if cls == Alignment:
            db = seq.data.annotation_db
        else:
            db = seq.annotation_db
        assert db is coll_db


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment))
def test_add_to_seq_updates_coll(cls, gff_db):  # ported for SequenceCollection
    """annotating a seq updates the db of the propagated"""
    seq_coll = cls(
        {"x": "AACCCAAAATTTTTTGGGGGGGGGGCCCC", "y": "AACCCAAAATTTTTTGGGGGGGGGGCCCC"}
    )
    x = seq_coll.get_seq("x")
    assert len(seq_coll.annotation_db) == len(x.annotation_db) == 0
    x.add_feature(biotype="exon", name="E1", spans=[(3, 8)])
    assert len(seq_coll.annotation_db) == len(x.annotation_db) == 1


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment))
def test_assign_none(cls, gff_db):  # ported for SequenceCollection
    """assigning None to annotation_db breaks conection"""

    seq_coll = cls({"seq1": "ACGU", "seq2": "CGUA", "test_seq": "CCGU"})
    seq_coll.add_feature(seqid="seq1", biotype="xyz", name="abc", spans=[(1, 2)])
    seq_coll.add_feature(seqid="seq2", biotype="xyzzz", name="abc", spans=[(1, 2)])
    seq_coll.annotation_db = None
    assert seq_coll.annotation_db is None


def test_copy_annotations_incompat_fails(seqcoll_db, gb_db):  # ported
    """copy_annotations copies records from annotation db"""
    db = seqcoll_db.annotation_db
    seqcoll = seqcoll_db.rename_seqs(
        lambda x: "AE017341" if x == seqcoll_db.names[0] else seqcoll_db.names[0]
    )
    seqcoll.annotation_db = db
    with pytest.raises(TypeError):
        seqcoll.copy_annotations(gb_db)


def test_copy_annotations_incompat_type_fails(seqcoll_db, gb_db):  # ported
    """copy_annotations copies records from annotation db"""
    with pytest.raises(TypeError):
        seqcoll_db.copy_annotations({"a": "ACGGT"})


def test_aligned_deepcopy_sliced():
    a = Aligned(*DNA.make_seq(seq="GCAAGGCGCCAA").parse_out_gaps())
    sliced = a[2:3]
    sl_cp = sliced.deepcopy(sliced=True)
    assert str(sl_cp) == str(sliced)


def test_deepcopy_with_features(DATA_DIR):
    """correctly deepcopy Aligned objects in an alignment"""
    path = DATA_DIR / "brca1_5.paml"
    # generates an annotatable Alignment object
    aln = load_aligned_seqs(path, array_align=False, moltype="dna")
    db = GffAnnotationDb()
    # when the annotation is outside(before) boundary of the slice
    db.add_feature(seqid="NineBande", biotype="exon", name="annot1", spans=[(0, 10)])
    # when the annotation is across boundary of the slice
    db.add_feature(seqid="Mouse", biotype="exon", name="annot2", spans=[(10, 21)])
    # when the annotation is within boundary of the slice
    db.add_feature(seqid="Human", biotype="exon", name="annot3", spans=[(20, 25)])
    # when the annotation is across boundary of the slice
    db.add_feature(seqid="HowlerMon", biotype="exon", name="annot4", spans=[(25, 32)])
    # when the annotation is outside(after) boundary of the slice
    db.add_feature(seqid="DogFaced", biotype="exon", name="annot5", spans=[(40, 45)])
    aln.annotation_db = db
    aln = aln[20:30]
    # no slice
    copied = aln.deepcopy(sliced=False)
    feats = list(copied.get_features(biotype="exon", allow_partial=True))
    assert len(feats) == 3  # overlap is Mouse, Human, HowlerMon
    # with slice
    copied = aln.deepcopy(sliced=True)
    feats = list(copied.get_features(biotype="exon", allow_partial=True))
    assert len(feats) == 3
    # rc drops annotations only when sliced is True
    rced = aln.rc()
    copied = rced.deepcopy(sliced=False)
    feats = list(copied.get_features(biotype="exon", allow_partial=True))
    assert len(feats) == 3  # overlap is Mouse, Human, HowlerMon


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment))
def test_deepcopy_aligned(cls):
    """correctly deep copy aligned objects in an alignment"""
    data = {"seq1": "ACGACGACG", "seq2": "ACGACGACG"}
    seqs = cls(data)
    copied = seqs.deepcopy(sliced=True)
    orig_rd = seqs.to_rich_dict()
    cpy_rd = copied.to_rich_dict()
    assert orig_rd == cpy_rd
    assert id(copied) != id(seqs)
    for name in seqs.names:
        assert id(copied.named_seqs[name]) != copied.named_seqs[name]


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment))
def test_seq_rename_preserves_annotations(cls):
    """rename seqs discards all annotations"""
    data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    seqs = cls(data, moltype=DNA)
    seqs.add_feature(seqid="seq1", biotype="exon", name="fred", spans=[(3, 8)])
    assert seqs.annotation_db is not None
    new = seqs.rename_seqs(lambda x: x.upper())
    assert len(new.annotation_db) == 1
    assert len(list(new.get_features(biotype="exon")))
    # using original seq name should also work
    assert len(list(new.get_features(seqid="seq1")))


def test_to_rich_dict_not_alignment():  # ported
    """to_rich_dict produces correct dict"""
    data = {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}
    aln = SequenceCollection(data, moltype="dna")
    try:
        seq_type = get_object_provenance(aln.seqs[0].data)
    except AttributeError:
        seq_type = get_object_provenance(aln.seqs[0])

    got = aln.to_rich_dict()

    data = {k: SeqView(seq=s, seqid=k).to_rich_dict() for k, s in data.items()}

    expect = {
        "seqs": {
            "seq1": {
                "name": "seq1",
                "seq": data["seq1"],
                "info": None,
                "type": seq_type,
                "moltype": aln.moltype.label,
                "version": __version__,
                "annotation_offset": 0,
            },
            "seq2": {
                "name": "seq2",
                "seq": data["seq2"],
                "info": None,
                "type": seq_type,
                "moltype": aln.moltype.label,
                "version": __version__,
                "annotation_offset": 0,
            },
            "seq3": {
                "name": "seq3",
                "seq": data["seq3"],
                "info": None,
                "type": seq_type,
                "moltype": aln.moltype.label,
                "version": __version__,
                "annotation_offset": 0,
            },
        },
        "moltype": aln.moltype.label,
        "info": None,
        "type": get_object_provenance(aln),
        "version": __version__,
    }
    assert got == expect


def test_array_alignment_to_rich_dict_not_alignment():  # will not port
    data = {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}
    aln = ArrayAlignment(data, moltype="dna")
    try:
        seq_type = get_object_provenance(aln.seqs[0].data)
    except AttributeError:
        seq_type = get_object_provenance(aln.seqs[0])

    got = aln.to_rich_dict()

    data = {k: SeqView(seq=s, seqid=k).to_rich_dict() for k, s in data.items()}

    expect = {
        "seqs": {
            "seq1": {
                "name": "seq1",
                "seq": data["seq1"],
                "info": None,
                "type": seq_type,
                "moltype": aln.moltype.label,
                "version": __version__,
                "annotation_offset": 0,
            },
            "seq2": {
                "name": "seq2",
                "seq": data["seq2"],
                "info": None,
                "type": seq_type,
                "moltype": aln.moltype.label,
                "version": __version__,
                "annotation_offset": 0,
            },
            "seq3": {
                "name": "seq3",
                "seq": data["seq3"],
                "info": None,
                "type": seq_type,
                "moltype": aln.moltype.label,
                "version": __version__,
                "annotation_offset": 0,
            },
        },
        "moltype": aln.moltype.label,
        "info": None,
        "type": get_object_provenance(aln),
        "version": __version__,
    }
    assert got == expect


def test_to_rich_dict_alignment():
    """to_rich_dict produces correct dict"""
    data = {"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"}
    aln = Alignment(data, moltype="dna")
    got = aln.to_rich_dict()

    data = {
        k: Aligned(
            *make_seq(seq=s, name=k, moltype="dna").parse_out_gaps()
        ).to_rich_dict()
        for k, s in data.items()
    }

    expect = {
        "seqs": {
            "seq1": data["seq1"],
            "seq2": data["seq2"],
            "seq3": data["seq3"],
        },
        "moltype": aln.moltype.label,
        "info": None,
        "type": get_object_provenance(aln),
        "version": __version__,
    }
    assert got == expect


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment))
def test_dotplot_annotated(cls):  # ported for SequenceCollection
    """exercising dotplot method with annotated sequences"""
    db = GffAnnotationDb()
    db.add_feature(seqid="Human", biotype="exon", name="fred", spans=[(10, 15)])

    seqs = cls(data={"Human": "CAGATTTGGCAGTT-", "Mouse": "CAGATTCAGCAGGTG"})
    seqs.annotation_db = db
    seqs = seqs.take_seqs(["Human", "Mouse"])
    _ = seqs.dotplot(show_progress=False)


def test_array_alignment_get_seq_entropy():
    """ArrayAlignment get_seq_entropy should get entropy of each seq"""
    seqs = [AB.make_seq(seq=s, preserve_case=True) for s in ["abab", "bbbb", "abbb"]]
    a = ArrayAlignment(seqs, moltype=AB)
    entropy = a.entropy_per_seq()
    e = 0.81127812445913283  # sum(p log_2 p) for p = 0.25, 0.75
    assert_allclose(entropy, array([1, 0, e]))


def test_get_seq_entropy():
    """Alignment get_seq_entropy should get entropy of each seq"""
    seqs = [AB.make_seq(seq=s, preserve_case=True) for s in ["abab", "bbbb", "abbb"]]
    a = Alignment(seqs, moltype=AB)
    entropy = a.entropy_per_seq()
    e = 0.81127812445913283  # sum(p log_2 p) for p = 0.25, 0.75
    assert_allclose(entropy, array([1, 0, e]))


@pytest.mark.parametrize("moltype", ("dna", "rna"))
def test_distance_matrix_singleton_collection(moltype):  # ported
    """SequenceCollection.distance_matrix() should raise error if collection
    only contains a single sequence"""
    collection = make_unaligned_seqs(data={"s1": "ACGTACGTAGTCGCG"}, moltype=moltype)
    with pytest.raises(ValueError):
        _ = collection.distance_matrix()


@pytest.mark.parametrize("moltype", ("dna", "rna"))
def test_collection_distance_matrix_same_seq(moltype):  # ported
    """Identical seqs should return distance measure of 0.0"""
    data = dict(
        [("s1", "ACGTACGTAGTCGCG"), ("s2", "GTGTACGTATCGCG"), ("s3", "GTGTACGTATCGCG")]
    )
    collection = make_unaligned_seqs(data=data, moltype=moltype)
    dists = collection.distance_matrix(calc="pdist")

    # all comparison of a sequence to itself should be zero
    for seq in collection.names:
        assert dists[(seq, seq)] == 0.0

    # s2 and s3 are identical, so should be zero
    assert dists[("s2", "s3")] == 0.0
    assert dists[("s3", "s2")] == 0.0


@pytest.mark.parametrize("moltype", ("protein", "text", "bytes"))
def test_distance_matrix_fails_wrong_moltype(moltype):  # ported
    data = [("s1", "ACGTA"), ("s2", "ACGTA")]
    seqs = make_unaligned_seqs(data=data, moltype=moltype)
    with pytest.raises(NotImplementedError):
        seqs.distance_matrix()


@pytest.mark.parametrize("moltype", ("dna", "rna"))
def test_distance_matrix_passes_correct_moltype(moltype):  # ported
    data = [("s1", "ACGTA"), ("s2", "ACGTA")]
    seqs = make_unaligned_seqs(data=data, moltype=moltype)
    seqs.distance_matrix()


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
@pytest.mark.parametrize(
    "seqs", ({"seq1": "GATTTT", "seq2": "GATC??"}, {"seq1": "GAT---", "seq2": "?GATCT"})
)
def test_get_translation2(cls, seqs):  # ported for SequenceCollection
    """SequenceCollection.get_translation translates each seq"""
    alignment = cls(data=seqs, moltype=DNA)
    got = alignment.get_translation()
    assert len(got) == 2
    assert got.moltype == PROTEIN


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
def test_get_translation_with_stop(cls):  # ported for SequenceCollection
    seqs = {"seq1": "GATTAG", "seq2": "?GATCT"}
    alignment = cls(data=seqs, moltype=DNA)
    got = alignment.get_translation(include_stop=True)
    assert got.to_dict() == {"seq1": "D*", "seq2": "XS"}


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment, SequenceCollection))
def test_get_translation_trim_stop(cls):  # ported for SequenceCollection
    seqs = {"seq1": "GATTCCTAG", "seq2": "GATTCCTCC"}
    alignment = cls(data=seqs, moltype=DNA)
    expect = {"seq1": "DS", "seq2": "DSS"}
    if cls != SequenceCollection:
        expect = {"seq1": "DS-", "seq2": "DSS"}

    got = alignment.get_translation(trim_stop=True)
    assert got.to_dict() == expect


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
@pytest.mark.parametrize(
    "seqs", ({"seq1": "GATTTT", "seq2": "GATC??"}, {"seq1": "GAT---", "seq2": "?GATCT"})
)
def test_get_translation_error(cls, seqs):  # ported for SequenceCollection
    """SequenceCollection.get_translation translates each seq"""
    # check for a failure when no moltype specified
    alignment = cls(data=seqs)
    with pytest.raises(TypeError):
        alignment.get_translation()


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
@pytest.mark.parametrize("method", ("ic_score", "cogent3_score", "sp_score"))
def test_alignment_quality_methods(cls, method):
    data = {
        "DogFaced": "TG----AATATGT------GAAAGAG",
        "FreeTaile": "TTGAAGAATATGT------GAAAGAG",
        "LittleBro": "CTGAAGAACCTGTGAAAGTGAAAGAG",
    }
    expected_score = dict(
        cogent3_score=-123.0, ic_score=32.93032499, sp_score=-25.25687109
    )[method]
    aln = cls(
        data,
        moltype="dna",
        info=dict(align_params=dict(lnL=-123.0)),
    )
    app = get_app(method)
    score = app(aln)
    assert_allclose(score, expected_score)


@pytest.mark.parametrize("method", ("ic_score", "cogent3_score", "sp_score"))
def test_alignment_quality_methods_oneseq(method):
    data = {
        "DogFaced": "TG----AATATGT------GAAAGAG",
    }
    aln = make_aligned_seqs(
        data,
        moltype="dna",
        info=dict(align_params=dict(lnL=-123.0)),
    )
    app = get_app(method)
    score = app(aln)
    assert_allclose(score, 0.0)


@pytest.mark.parametrize("method", ("ic_score", "cogent3_score", "sp_score"))
def test_alignment_quality_methods_zero_length(method):
    data = {
        "a": "",
        "b": "",
        "c": "",
    }
    aln = make_aligned_seqs(
        data,
        moltype="dna",
        info=dict(align_params=dict(lnL=-123.0)),
    )
    app = get_app(method)
    score = app(aln)
    assert_allclose(score, 0.0)


def test_get_gap_array_equivalence():
    # make sure produced gap arrays are identical between the
    # two Alignment classes
    data = {
        "DogFaced": "TG----AATATGT------GAAAGAG",
        "FreeTaile": "TTGAAGAATATGT------GAAAGAG",
        "LittleBro": "CTGAAGAACCTGTGAAAGTGAAAGAG",
    }
    array_aln = make_aligned_seqs(data, moltype="dna", array_align=True)
    aln = make_aligned_seqs(data, moltype="dna", array_align=False)
    assert_allclose(array_aln.get_gap_array(), aln.get_gap_array())


@pytest.mark.parametrize("reverse", (False, True))
def test_aligned_rich_dict(reverse):
    map_, s = make_seq(
        seq="TTGAAGAATATGT------GAAAGAG", name="s1", moltype="dna"
    ).parse_out_gaps()
    seq = Aligned(map_, s)
    if reverse:
        seq = seq.rc()

    rd = seq.to_rich_dict()
    got = Aligned.from_rich_dict(rd)
    assert str(got) == str(seq)


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
@pytest.mark.parametrize(
    "seqs",
    (
        {"seq1": "GATTTT", "seq2": "GATC??"},
        {"seq1": "GAT---", "seq2": "?GATCT"},
    ),
)
def test_get_translation_info(cls, seqs):  # ported for SequenceCollection
    """SequenceCollection.get_translation preserves info attribute"""
    alignment = cls(data=seqs, moltype=DNA, info={"key": "value"})
    got = alignment.get_translation()
    assert got.info["key"] == "value"


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
@pytest.mark.parametrize(
    "gc,seqs",
    (
        (1, ("TCCTGA", "GATTT?")),
        (1, ("ACGTAA---", "ACGAC----", "ACGCAATGA")),
        (2, ("GATTTT", "TCCAGG")),
    ),
)
def test_has_terminal_stop_true(cls, gc, seqs):  # ported for SequenceCollection
    gc = get_code(gc)
    data = {f"s{i}": s for i, s in enumerate(seqs)}
    seqs = cls(data=data, moltype="dna")
    assert seqs.has_terminal_stop(gc=gc)


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
@pytest.mark.parametrize(
    "gc,seqs",
    ((1, ("TCCTCA", "GATTTT")), (2, ("GATTTT", "TCCCGG")), (1, ("CCTCA", "ATTTT"))),
)
def test_has_terminal_stop_false(cls, gc, seqs):  # ported for SequenceCollection
    gc = get_code(gc)
    data = {f"s{i}": s for i, s in enumerate(seqs)}
    seqs = cls(data=data, moltype="dna")
    assert not seqs.has_terminal_stop(gc=gc)


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
def test_has_terminal_stop_strict(cls):  # ported for SequenceCollection
    gc = get_code(1)
    data = {f"s{i}": s for i, s in enumerate(("CCTCA", "ATTTT"))}
    seqs = cls(data=data, moltype="dna")
    with pytest.raises(AlphabetError):
        seqs.has_terminal_stop(gc=gc, strict=True)


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
@pytest.mark.parametrize(
    "gc,seqs",
    (
        (1, ("--AT-CTGA", "GATAAATT?")),
        (1, ("ACGTGA---", "ACGAC----", "ACGCAATGA")),
        (1, ("CCTCA-", "ATTTTA")),
        (2, ("GATTTT", "TCCAGG")),
    ),
)
def test_trim_stops_true(cls, gc, seqs):  # ported for SequenceCollection
    gc = get_code(gc)
    data = {f"s{i}": s for i, s in enumerate(seqs)}

    expect = {}
    for k, v in data.items():
        if cls != SequenceCollection or "-" in v:
            v = re.sub("(TGA|AGG)(?=[-]*$)", "---", v)
        else:
            v = re.sub("(TGA|AGG)", "", v)
        expect[k] = v

    seqs = cls(data=data, moltype="dna")
    got = seqs.trim_stop_codons(gc=gc).to_dict()

    assert got == expect


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
@pytest.mark.parametrize(
    "gc,seqs",
    ((1, ("T-CTGC", "GATAA?")), (2, ("GATTTT", "TCCCGG")), (1, ("CCTGC", "GATAA"))),
)
def test_trim_terminal_stops_nostop(cls, gc, seqs):  # ported for SequenceCollection
    gc = get_code(gc)
    data = {f"s{i}": s for i, s in enumerate(seqs)}
    seqs = cls(data=data, moltype="dna")
    got = seqs.trim_stop_codons(gc=gc)
    assert got is seqs


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
@pytest.mark.parametrize("seqs", (("CCTCA", "ATTTT"), ("CCTCA-", "ATTTTA")))
def test_trim_terminal_stops_strict(cls, seqs):  # ported for SequenceCollection
    gc = get_code(1)
    data = {f"s{i}": s for i, s in enumerate(seqs)}
    seqs = cls(data=data, moltype="dna")
    with pytest.raises(AlphabetError):
        seqs.trim_stop_codons(gc=gc, strict=True)


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
def test_trim_stop_codons_info(cls):  # ported for SequenceCollection
    """trim_stop_codons should preserve info attribute"""
    coll = cls(
        data={"seq1": "ACGTAA", "seq2": "ACGACG", "seq3": "ACGCGT"},
        moltype=DNA,
        info={"key": "value"},
    )
    coll = coll.trim_stop_codons()
    assert coll.info["key"] == "value"


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
def test_get_translation_incomplete(cls):  # ported for SequenceCollection
    """get translation works on incomplete codons"""
    alignment = cls(data={"seq1": "GATN--", "seq2": "?GATCT"}, moltype=DNA)
    got = alignment.get_translation(incomplete_ok=True)
    assert got.to_dict() == {"seq1": "D?", "seq2": "XS"}
    with pytest.raises(AlphabetError):
        _ = alignment.get_translation(incomplete_ok=False)


@pytest.mark.parametrize("name", ("s1", "s2", "s3"))
def test_get_seq_with_sliced_aln(name):
    seqs = {
        "s1": "GTTGAAGTAGTAGAAGTTCCAAATAATGAA",
        "s2": "GTG------GTAGAAGTTCCAAATAATGAA",
        "s3": "GCTGAAGTAGTGGAAGTTGCAAAT---GAA",
    }
    aln = make_aligned_seqs(data=seqs, moltype="dna", array_align=False)
    start, stop = 1, 5
    a1 = aln[start:stop]

    seq = a1.get_seq(name)
    assert isinstance(seq, Sequence), seq

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
    aln = make_aligned_seqs(data=seqs, moltype="dna", array_align=False)
    start, stop = 1, 5
    a1 = aln[start:stop]
    a1 = a1.rc()
    got = str(a1.get_seq(name))

    dna = get_moltype("dna")
    expect = dna.complement(seqs[name][start:stop].replace("-", ""))[::-1]
    assert got == expect, (got, expect)


@pytest.mark.parametrize("name", ("s1", "s2", "s3"))
def test_get_seq_with_sliced_aln_multiple_spans(name):
    seqs = {  # the sliced seq has:
        "s1": "GTTGA--TAGTAGAAGTTCCAAATAATGAA",  # span gap span
        "s2": "G----TT------AAGTTCCAAATAATGAA",  # gap span gap
        "s3": "G--GA--TA--GGAAGTTGCAAAT---GAA",  # gap span gap span gap
    }
    aln = make_aligned_seqs(data=seqs, moltype="dna", array_align=False)
    start, stop = 1, 10
    a1 = aln[start:stop]
    seq = a1.get_seq(name)
    assert isinstance(seq, Sequence), seq

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
    aln = make_aligned_seqs(data=seqs, moltype="dna", array_align=False)
    start, stop = 1, 10
    a1 = aln[start:stop]
    a1 = a1.rc()
    got = str(a1.get_seq(name))
    dna = get_moltype("dna")
    expect = dna.complement(seqs[name][start:stop].replace("-", ""))[::-1]
    assert got == expect, (got, expect)


@pytest.mark.parametrize("name", ("s1", "s2", "s3"))
def test_get_gapped_seq_with_sliced_aln(name):
    seqs = {
        "s1": "G-TG---TAGTAGAAGTTCCAAATAATGAA",
        "s2": "GTG------GTAGAAGTTCCAAATAATGAA",
        "s3": "GC--AAGTAGTGGAAGTTGCAAAT---GAA",
    }
    aln = make_aligned_seqs(data=seqs, moltype="dna", array_align=False)
    start, stop = 1, 10
    a1 = aln[start:stop]

    seq = a1.get_gapped_seq(name)
    assert isinstance(seq, Sequence), seq

    got = str(seq)
    expect = seqs[name][start:stop]
    assert got == expect, (got, expect)


@pytest.mark.parametrize("name", ("s1", "s2", "s3"))
@pytest.mark.parametrize("array_align", (True, False))
def test_aln_rev_slice(name, array_align):
    seqs = {
        "s1": "AAGGTTCC",
        "s2": "AAGGTTCC",
        "s3": "AAGGTTCC",
    }

    aln = make_aligned_seqs(data=seqs, moltype="dna", array_align=array_align)
    got = aln[5:1]
    assert not got

    seq = got.get_gapped_seq(name)
    assert not str(seq)

    seq = got.get_seq(name)
    assert not str(seq)


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
def test_init_duplicate_keys_raises(cls):
    """SequenceCollection init from (key, val) pairs should fail on dup. keys"""
    data = [["x", "XXX"], ["b", "BBB"], ["x", "CCC"], ["d", "DDD"], ["a", "AAA"]]
    with pytest.raises(ValueError):
        cls(data)


@pytest.mark.parametrize(
    "seq",
    (
        "ACGG",
        DNA.make_seq(seq="ACGG"),
        numpy.array([2, 1, 3, 3], dtype=int),
        Aligned(*DNA.make_seq(seq="AC-GG").parse_out_gaps()),
    ),
)
def test_construct_unaligned_seq(seq):
    moltype = get_moltype("dna")
    got = _construct_unaligned_seq(seq, name="seq1", moltype=moltype)
    assert isinstance(got, (Sequence, ArraySequence))


@pytest.mark.parametrize(
    "data",
    (
        [["s1", "ACGG"]],
        {"s1": "ACGG"},
        {"s1": DNA.make_seq(seq="ACGG")},
        (
            Aligned(
                *DNA.make_seq(seq="AC-GG", name="s1").parse_out_gaps(),
            ),
        ),
        iter([["s1", "ACGG"]]),
    ),
)
def test_to_unaligned_seqs(data):
    moltype = get_moltype("dna")
    got = _coerce_to_unaligned_seqs(data, names=None, moltype=moltype)
    assert isinstance(got[0], list)
    assert isinstance(got[0][0], (Sequence, ArraySequence))


@pytest.mark.parametrize("moltype", ("dna", "protein"))
@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
def test_same_moltype(cls, moltype):  # ported for SequenceCollection
    moltype = get_moltype(moltype)
    data = dict(s1="ACGTT", s2="ACCTT")
    seqs = cls(data, moltype=moltype)
    got = seqs.to_moltype(moltype=moltype)
    assert got is seqs


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
def test_add(cls):  # will not port for SequenceCollection
    """__add__ should concatenate sequence data, by name"""
    align1 = cls({"a": "AAAA", "b": "TTTT", "c": "CCCC"})
    align2 = cls({"a": "GGGG", "b": "----", "c": "NNNN"})
    align = align1 + align2
    concatdict = align.to_dict()
    assert concatdict == {"a": "AAAAGGGG", "b": "TTTT----", "c": "CCCCNNNN"}


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_count_gaps_per_pos(cls):
    """correctly compute the number of gaps"""
    data = {"a": "AAAA---GGT", "b": "CCC--GG?GT"}
    aln = cls(data=data, moltype=DNA)
    # per position
    got = aln.count_gaps_per_pos(include_ambiguity=False)
    assert_equal(got.array, [0, 0, 0, 1, 2, 1, 1, 0, 0, 0])
    got = aln.count_gaps_per_pos(include_ambiguity=True)
    assert_equal(got.array, [0, 0, 0, 1, 2, 1, 1, 1, 0, 0])


def _make_and_filter(cls, raw, expected, motif_length, drop_remainder):
    # a simple filter func
    aln = cls(raw, info={"key": "value"})
    func = _make_filter_func(aln)
    result = aln.filtered(
        func,
        motif_length=motif_length,
        warn=False,
        drop_remainder=drop_remainder,
    )
    assert result.to_dict() == expected
    assert result.info["key"] == "value"


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_filtered(cls):
    """filtered should return new alignment with positions consistent with
    provided callback function"""
    # a simple filter option
    raw = {"a": "ACGACGACG", "b": "CCC---CCC", "c": "AAAA--AAA"}
    _make_and_filter(cls, raw, {"a": "ACGACG", "b": "CCCCCC", "c": "AAAAAA"}, 1, True)
    # check with motif_length = 2
    _make_and_filter(cls, raw, {"a": "ACAC", "b": "CCCC", "c": "AAAA"}, 2, True)
    # check with motif_length = 3
    _make_and_filter(cls, raw, {"a": "ACGACG", "b": "CCCCCC", "c": "AAAAAA"}, 3, True)


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
def test_counts_per_seq(cls):  # ported
    """SequenceCollection.counts_per_seq handles motif length, allow_gaps etc.."""
    data = {"a": "AAAA??????", "b": "CCCGGG--NN", "c": "CCGGTTCCAA"}
    coll = cls(data=data, moltype="dna")
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


def test_counts_per_seq_arr_alignment():
    """ArrayAlignment counts_per_seq should return motif counts each seq"""
    seqs = []
    for s in ["abaa", "abbb"]:
        seqs.append(AB.make_seq(seq=s, preserve_case=True))
    a = ArrayAlignment(seqs, moltype=AB)
    f = a.counts_per_seq()
    assert_equal(f.array, array([[3, 1], [1, 3]]))
    f = a.counts_per_seq(motif_length=2, exclude_unobserved=True)
    assert_equal(f.array, array([[1, 1, 0], [0, 1, 1]]))


def test_coevolution_segments():
    """specifying coordinate segments produces matrix with just those"""
    aln = load_aligned_seqs("data/brca1.fasta", moltype="dna")
    aln = aln.take_seqs(aln.names[:20])
    aln = aln.no_degenerates()[:20]
    coevo = aln.coevolution(segments=[(4, 6), (11, 13)], show_progress=False)
    assert coevo.template.names[0] == [4, 5, 11, 12]


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_iter_positions(cls):
    data = {"a": "AAAAAA", "b": "AAA---", "c": "AAAA--"}
    r = cls({k: data[k] for k in "cb"})
    assert list(r.iter_positions(pos_order=[5, 1, 3])) == list(
        map(list, ["--", "AA", "A-"])
    )
    # reorder names
    r = cls(data)
    cols = list(r.iter_positions())
    assert cols == list(map(list, ["AAA", "AAA", "AAA", "A-A", "A--", "A--"]))


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_add_from_ref_aln(cls):
    """should add or insert seqs based on align to reference"""
    aln1 = cls(
        [
            ["name1", "-AC-DEFGHI---"],
            ["name2", "XXXXXX--XXXXX"],
            ["name3", "YYYY-YYYYYYYY"],
        ]
    )

    aln2 = cls(
        [
            ["name1", "ACDEFGHI"],
            ["name4", "KL--MNPR"],
            ["name5", "KLACMNPR"],
            ["name6", "KL--MNPR"],
        ]
    )

    aligned_to_ref_out_aln_inserted = cls(
        [
            ["name1", "-AC-DEFGHI---"],
            ["name4", "-KL---MNPR---"],
            ["name5", "-KL-ACMNPR---"],
            ["name6", "-KL---MNPR---"],
            ["name2", "XXXXXX--XXXXX"],
            ["name3", "YYYY-YYYYYYYY"],
        ]
    )

    aln2_wrong_refseq = cls((("name1", "ACDXFGHI"), ("name4", "KL--MNPR")))

    aln2_wrong_refseq_name = cls([["nameY", "ACDEFGHI"], ["name4", "KL--MNPR"]])

    aln2_different_aln_class = ArrayAlignment(
        [["name1", "ACDEFGHI"], ["name4", "KL--MNPR"]]
    )

    aligned_to_ref_out_aln = cls(
        [
            ["name1", "-AC-DEFGHI---"],
            ["name2", "XXXXXX--XXXXX"],
            ["name3", "YYYY-YYYYYYYY"],
            ["name4", "-KL---MNPR---"],
        ]
    )

    out_aln = aln1.add_from_ref_aln(aln2, after_name="name1")
    assert str(aligned_to_ref_out_aln_inserted) == str(out_aln)  # test insert_after

    out_aln = aln1.add_from_ref_aln(aln2, before_name="name2")
    assert aligned_to_ref_out_aln_inserted == out_aln  # test insert_before

    with pytest.raises(ValueError):
        aln1.add_from_ref_aln(aln2_wrong_refseq_name)

    aln = aln1.add_from_ref_aln(aln2_different_aln_class)
    # test_align_to_refseq_different_aln_class
    assert aligned_to_ref_out_aln == aln

    with pytest.raises(ValueError):
        # test wrong_refseq
        aln1.add_from_ref_aln(aln2_wrong_refseq)


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_replace_seqs(cls):
    """replace_seqs should replace 1-letter w/ 3-letter seqs"""
    a = cls({"seq1": "ACGU", "seq2": "C-UA", "seq3": "C---"})
    seqs = {"seq1": "AAACCCGGGUUU", "seq2": "CCCUUUAAA", "seq3": "CCC"}
    result = a.replace_seqs(seqs)  # default behaviour
    assert result.to_dict() == {
        "seq1": "AAACCCGGGUUU",
        "seq2": "CCC---UUUAAA",
        "seq3": "CCC---------",
    }

    result = a.replace_seqs(seqs, aa_to_codon=True)  # default behaviour
    assert result.to_dict() == {
        "seq1": "AAACCCGGGUUU",
        "seq2": "CCC---UUUAAA",
        "seq3": "CCC---------",
    }

    # should correctly gap the same sequences with same length
    result = a.replace_seqs(a.degap(), aa_to_codon=False)  # default behaviour
    assert result.to_dict() == {"seq1": "ACGU", "seq2": "C-UA", "seq3": "C---"}

    # should fail when not same length if aa_to_codon is False
    new = SequenceCollection(
        [(n, s.replace("-", "")) for n, s in list(a[:3].to_dict().items())],
        moltype="bytes",
    )
    with pytest.raises(ValueError):
        a.replace_seqs(new, aa_to_codon=False)

    # check the gaps are changed
    aln1 = cls(data={"a": "AC-CT", "b": "ACGCT"})
    aln2 = cls(data={"a": "ACC-T", "b": "ACGCT"})

    result = aln1.replace_seqs(aln2, aa_to_codon=False)
    assert id(aln1) != id(aln2)
    assert aln1.to_dict() == result.to_dict()


def test_get_alphabet_and_moltype():
    """ArrayAlignment should figure out correct alphabet and moltype"""
    s1 = "A"
    s2 = RNA.make_seq(seq="AA")

    d = ArrayAlignment(s1)
    assert d.moltype is BYTES
    assert d.alphabet is BYTES.alphabet

    d = ArrayAlignment(s1, moltype=RNA)
    assert d.moltype is RNA
    assert d.alphabet is RNA.alphabets.degen_gapped

    d = ArrayAlignment([s2])
    assert d.moltype is RNA
    assert d.alphabet is RNA.alphabets.degen_gapped

    d = ArrayAlignment(s2, moltype=DNA)
    assert d.moltype is DNA
    assert d.alphabet is DNA.alphabets.degen_gapped
    # checks for containers
    d = ArrayAlignment([s2])
    assert d.moltype is RNA
    d = ArrayAlignment({"x": s2})
    assert d.moltype is RNA
    d = ArrayAlignment(set([s2]))
    assert d.moltype is RNA


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment, SequenceCollection))
def test_aln_from_fasta_parser(cls):
    """aln_from_fasta_parser should init from iterator"""
    s = ">aa\nAC\n>bb\nAA\n>c\nGG\n".splitlines()
    p = MinimalFastaParser(s)
    aln = cls(p, moltype=DNA)
    assert aln.named_seqs["aa"] == "AC"
    assert aln.to_dict() == {"aa": "AC", "bb": "AA", "c": "GG"}
    s2 = ">aa\nAC\n>bb\nAA\n>c\nGG\n"
    d = cls(MinimalFastaParser(s2.splitlines()))
    assert d.to_dict() == aln.to_dict()


@pytest.mark.parametrize(
    "d", ({"a": "AAAAA", "b": "BBBBB"}, {"a": b"AAAAA", "b": b"BBBBB"})
)
@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment, SequenceCollection))
def test_init_dict(cls, d):  # will not port for SequenceCollection
    """SequenceCollection init from dict should work as expected"""
    # from bytes strings
    a = cls(d, names=["a", "b"])
    expect = {k: v.decode("utf8") if isinstance(v, bytes) else v for k, v in d.items()}
    assert a.to_dict() == expect


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment, SequenceCollection))
def test_init_seq_info(cls):
    """SequenceCollection init from seqs w/ info should preserve data"""
    a = Sequence(seq="AAA", name="a", info={"x": 3})
    b = Sequence(seq="CCC", name="b", info={"x": 4})
    c = Sequence(seq="GGG", name="c", info={"x": 5})
    seqs = [c, b, a]
    a = cls(seqs)
    assert list(a.names) == ["c", "b", "a"]
    assert list(map(str, a.seqs)) == ["GGG", "CCC", "AAA"]
    if cls is not ArrayAlignment:
        # ArrayAlignment is allowed to strip info objects
        assert [i.info.x for i in a.seqs] == [5, 4, 3]
    # check it still works if constructed from same class
    b = cls(a)
    assert list(b.names) == ["c", "b", "a"]
    assert list(map(str, b.seqs)) == ["GGG", "CCC", "AAA"]
    if cls is not ArrayAlignment:
        # ArrayAlignment is allowed to strip Info objects
        assert [i.info.x for i in b.seqs] == [5, 4, 3]


@pytest.mark.parametrize("val", (None, 3))
def test_init_empty(val):
    """ArrayAlignment init should fail if empty."""
    with pytest.raises((ValueError, TypeError)):
        ArrayAlignment(val)
        ArrayAlignment()


def test_make_case():
    data = {"a": "tata", "b": "tgtc", "c": "gcga", "d": "gaac", "e": "gagc"}
    got = make_aligned_seqs(data=data, moltype="dna")
    assert got == {k: v.upper() for k, v in data.items()}


@pytest.fixture(scope="session")
def brca1_data():
    return load_aligned_seqs("data/brca1.fasta").to_dict()


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment, SequenceCollection))
def test_dotplot(cls, brca1_data):  # ported for SequenceCollection
    """exercising dotplot method"""
    seqs = cls(data=brca1_data, moltype=DNA)
    _ = seqs.dotplot()
    with pytest.raises(AssertionError):
        seqs.dotplot(window=5, k=11)


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_dotplot_alignment(cls):
    """exercising dotplot method"""
    aln = cls([["name1", "TTTTTTAAAA"], ["name2", "AAAATTTTTT"]])
    aln = aln[2:8]
    draw = aln.dotplot(show_progress=False)
    expected = {("name1", "TTTTAA"), ("name2", "AATTTT")}
    got = {(s.name, str(s)) for s in (draw.seq1, draw.seq2)}
    assert got == expected


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_get_degapped_relative_to(cls):
    """should remove all columns with a gap in sequence with given name"""
    aln = cls(
        [
            ["name1", "-AC-DEFGHI---"],
            ["name2", "XXXXXX--XXXXX"],
            ["name3", "YYYY-YYYYYYYY"],
            ["name4", "-KL---MNPR---"],
        ]
    )
    expect = dict(
        [
            ["name1", "ACDEFGHI"],
            ["name2", "XXXX--XX"],
            ["name3", "YY-YYYYY"],
            ["name4", "KL--MNPR"],
        ]
    )
    result = aln.get_degapped_relative_to("name1")
    assert result.to_dict() == expect

    with pytest.raises(ValueError):
        aln.get_degapped_relative_to("nameX")


@pytest.mark.parametrize(
    "data",
    (
        {"seq1": "ACAACGACG", "seq2": "ACGACGACG"},
        {"seq1": "-CAACGACG", "seq2": "ACGACGACG"},
        {"seq1": "ACAACGAC-", "seq2": "ACGACGACG"},
        {"seq1": "AC-ACGACG", "seq2": "ACGACGACG"},
        {"seq1": "ACA-CGAC-", "seq2": "ACGACGACG"},
    ),
)
@pytest.mark.parametrize("name", ("seq1", "seq2"))
@pytest.mark.parametrize("rev", (False, True))
def test_sliced_deepcopy(data, name, rev):
    """correctly deep copy aligned objects in an alignment"""

    orig = Alignment(data, moltype="dna")
    slice_start, slice_end = 3, 5
    aln = orig[slice_start:slice_end]
    if rev:
        aln = aln.rc()

    notsliced = aln.deepcopy(sliced=False)
    # the annotation offsets should match original object
    assert {s.data.annotation_offset for s in notsliced.seqs} == {
        s.data.annotation_offset for s in aln.seqs
    }
    sliced = aln.deepcopy(sliced=True)
    assert sliced.to_dict() == notsliced.to_dict()
    sliced_seq = sliced.named_seqs[name]
    notsliced_seq = notsliced.named_seqs[name]
    assert str(sliced_seq) == str(notsliced_seq)

    # the annotation offsets should match original object
    assert {s.data.annotation_offset for s in sliced.seqs} == {
        s.data.annotation_offset for s in aln.seqs
    }
    # if not sliced in copy, underlying seq data is identical to original
    assert (
        notsliced.named_seqs[name].data._seq.seq is orig.named_seqs[name].data._seq.seq
    )
    # but not the same for sliced
    assert (
        sliced.named_seqs[name].data._seq.seq is not orig.named_seqs[name].data._seq.seq
    )

    assert sliced.named_seqs[name].map.parent_length == len(
        str(sliced_seq).replace("-", "")
    )
    assert notsliced.named_seqs[name].map.parent_length == len(
        str(notsliced_seq).replace("-", "")
    )
    # and map.parent_length and len(data) should match
    assert sliced.named_seqs[name].map.parent_length == len(
        sliced.named_seqs[name].data
    )

    # if sliced, seq data should be < orig
    assert len(sliced.named_seqs[name].data) < len(orig.named_seqs[name].data)


def test_aligned_deepcopy_sliced_map_matches_data():
    m, seq = DNA.make_seq(seq="ACAACGACG", name="seq1").parse_out_gaps()
    aligned = Aligned(m, seq)
    sliced = aligned[3:5]
    sliced_copy = sliced.deepcopy(sliced=True)
    assert sliced_copy.map.parent_length == len(sliced_copy.data)


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
def test_iter_selected(cls):  # ported for SequenceCollection
    """SequenceCollection iter_selected() should iterate over items in correct order"""
    # should work if one row
    one_seq = cls({"a": "AAAAA"})
    ragged_padded = cls({"a": "AAAAAA", "b": "AAA---", "c": "AAAA--"})
    ordered1 = cls({"a": "AAAAA", "b": "BBBBB"}, names=["a", "b"])
    ordered2 = cls({"a": "AAAAA", "b": "BBBBB"}, names=["b", "a"])

    assert list(one_seq.iter_selected()) == (["A"] * 5)
    # should take order into account
    assert list(ordered1.iter_selected()) == (["A"] * 5 + ["B"] * 5)
    assert list(ordered2.iter_selected()) == (["B"] * 5 + ["A"] * 5)
    # should allow row and/or col specification
    r = ragged_padded
    assert list(r.iter_selected(seq_order=["c", "b"], pos_order=[5, 1, 3])) == list(
        "-AA-A-"
    )
    # should not interfere with superclass iteritems()
    i = list(r.named_seqs.items())
    i.sort()
    assert i == ([("a", "AAAAAA"), ("b", "AAA---"), ("c", "AAAA--")])


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_take_positions(cls):
    """SequenceCollection take_positions should return new alignment w/ specified pos"""
    gaps = cls({"a": "AAAAAAA", "b": "A--A-AA", "c": "AA-----"})
    assert gaps.take_positions([5, 4, 0]) == {"a": "AAA", "b": "A-A", "c": "--A"}
    assert isinstance(gaps.take_positions([0]), _SequenceCollectionBase)

    # should be able to negate
    assert gaps.take_positions([5, 4, 0], negate=True) == {
        "a": "AAAA",
        "b": "--AA",
        "c": "A---",
    }


@pytest.mark.parametrize("array_align", (True, False))
@pytest.mark.parametrize(
    "seq_moltype", (("ACCT--GT", "dna"), ("ACCT--GU", "rna"), ("MTST--T", "protein"))
)
def test_alignment_propogates_seqid_to_seqview(array_align, seq_moltype):
    data = {"seq1": seq_moltype[0], "seq2": seq_moltype[0], "seq3": seq_moltype[0]}
    aln = make_aligned_seqs(data, moltype=seq_moltype[1], array_align=array_align)

    if array_align:
        assert aln.seqs[0]._seq.seqid == "seq1"
    else:
        assert aln.seqs[0].data._seq.seqid == "seq1"


@pytest.mark.parametrize("cls", (str, bytes))
def test_construct_unaligned_seq_propogates_seqid(cls):
    data = "ACGT"
    if cls is bytes:
        seq = cls(data, "utf8")
    else:
        seq = cls(data)
    got = _construct_unaligned_seq(seq, name="seq1", moltype=DNA)
    assert got._seq.seqid == "seq1"


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_sample_info(cls):
    """Alignment.sample should preserver info attribute"""
    alignment = cls(
        {"seq1": "ABCDEFGHIJKLMNOP", "seq2": "ABCDEFGHIJKLMNOP"},
        info={"key": "value"},
    )
    # effectively permute columns, preserving length
    shuffled = alignment.sample()
    assert shuffled.info["key"] == "value"
    # ensure length correct
    sample = alignment.sample(10)
    assert sample.info["key"] == "value"


@pytest.mark.parametrize("cls", (SequenceCollection, Alignment, ArrayAlignment))
def test_to_json(cls):  # ported for SequenceCollection
    """roundtrip of to_json produces correct dict"""
    aln = cls({"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"})
    txt = aln.to_json()
    got = json.loads(txt)
    expect = aln.to_rich_dict()
    assert got == expect


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_get_position_indices(cls):
    """get_position_indices should return names of cols where f(col)"""

    def gap_1st(x):
        return x[0] == "-"

    def gap_2nd(x):
        return x[1] == "-"

    def gap_3rd(x):
        return x[2] == "-"

    def is_list(x):
        return isinstance(x, list)

    gaps = cls({"a": "AAAAAAA", "b": "A--A-AA", "c": "AA-----"}, names=["a", "b", "c"])

    assert gaps.get_position_indices(gap_1st) == []
    assert gaps.get_position_indices(gap_2nd) == [1, 2, 4]
    assert gaps.get_position_indices(gap_3rd) == [2, 3, 4, 5, 6]
    assert gaps.get_position_indices(is_list) == [0, 1, 2, 3, 4, 5, 6]
    # should be able to negate
    assert gaps.get_position_indices(gap_2nd, negate=True) == [0, 3, 5, 6]
    assert gaps.get_position_indices(gap_1st, negate=True) == [0, 1, 2, 3, 4, 5, 6]
    assert gaps.get_position_indices(is_list, negate=True) == []


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_positions(cls):
    """positions property should iterate over positions, using self.names"""
    r = cls({"a": "AAAAAA", "b": "AAA---", "c": "AAAA--"})
    r.names = ["a", "b", "c"]
    expect = [list(v) for v in ("AAA", "AAA", "AAA", "A-A", "A--", "A--")]
    assert list(r.positions) == expect


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment, SequenceCollection))
def test_set_repr_policy_valid_input(cls):  # ported for SequenceCollection
    """repr_policy should be set to new values"""
    seqs = cls({"a": "AAAAA", "b": "AAA--"})
    seqs.set_repr_policy(num_seqs=5, num_pos=40, ref_name="a", wrap=10)
    assert seqs._repr_policy == dict(num_seqs=5, num_pos=40, ref_name="a", wrap=10)


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_set_repr_policy_valid_input_slices(cls):
    """repr_policy should be set to new values"""
    seqs = cls({"a": "AAAAA", "b": "AAA--"})
    seqs.set_repr_policy(num_seqs=5, num_pos=40, ref_name="a", wrap=10)
    # should persist in slicing
    sliced = seqs[:2]
    assert sliced._repr_policy == dict(num_seqs=5, num_pos=40, ref_name="a", wrap=10)


def test_array_align_error_with_mixed_length():
    data = dict(s1="ACGG", s2="A-G")
    with pytest.raises(ValueError, match=".* not all the same length.*"):
        make_aligned_seqs(data=data)


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment, SequenceCollection))
def test_add_seqs(cls):
    """add_seqs should return an alignment with the new sequences appended or inserted"""
    data = [("name1", "AAA"), ("name2", "AAA"), ("name3", "AAA"), ("name4", "AAA")]
    data1 = [("name1", "AAA"), ("name2", "AAA")]
    data2 = [("name3", "AAA"), ("name4", "AAA")]
    data3 = [("name5", "BBB"), ("name6", "CCC")]
    aln = cls(data)
    aln3 = cls(data3)

    out_aln = aln.add_seqs(aln3)
    # test append at the end
    assert str(out_aln) == str(cls(data + data3))

    out_aln = aln.add_seqs(aln3, before_name="name3")
    assert str(out_aln) == str(cls(data1 + data3 + data2))

    # test insert before

    out_aln = aln.add_seqs(aln3, after_name="name2")
    assert str(out_aln) == str(cls(data1 + data3 + data2))  # test insert after

    out_aln = aln.add_seqs(aln3, before_name="name1")
    # test if insert before first seq works
    assert str(out_aln) == str(cls(data3 + data))

    out_aln = aln.add_seqs(aln3, after_name="name4")
    # test if insert after last seq works
    assert str(out_aln) == str(cls(data + data3))

    with pytest.raises(ValueError):
        # wrong after/before name
        aln.add_seqs(aln3, before_name="name5")

    with pytest.raises(ValueError):
        # wrong after/before name
        aln.add_seqs(aln3, after_name="name5")

    if isinstance(aln, Alignment) or isinstance(aln, ArrayAlignment):
        with pytest.raises((DataError, ValueError)):
            aln.add_seqs(aln3 + aln3)
    else:
        exp = set([seq for name, seq in data])
        exp.update([seq + seq for name, seq in data3])
        got = set()
        for seq in aln.add_seqs(aln3 + aln3).seqs:
            got.update([str(seq).strip()])
        assert got == exp


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_omit_gap_pos(cls):
    """Alignment omit_gap_pos should return alignment w/o positions of gaps"""
    aln = cls({"a": "--A-BC-", "b": "-CB-A--", "c": "--D-EF-"}, names=["a", "b", "c"])
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
    new_aln_data = aln.named_seqs.copy()
    new_aln_data["d"] = "-------"
    new_aln_data["e"] = "XYZXYZX"
    new_aln_data["f"] = "AB-CDEF"
    aln = cls(new_aln_data)

    # if no gaps are allowed, we get None
    assert aln.omit_gap_pos(0) is None


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
def test_no_degenerates(cls):
    """no_degenerates correctly excludes columns containing IUPAC ambiguity codes"""
    data = {
        "s1": "AAA CCC GGG TTT".replace(" ", ""),
        "s2": "CCC GGG T-T AAA".replace(" ", ""),
        "s3": "GGR YTT AAA CCC".replace(" ", ""),
    }
    aln = cls(data=data, moltype=DNA)

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


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
@pytest.mark.parametrize("moltype", ("bytes", "text"))
def test_no_degenerates_invalid_moltype(cls, moltype):
    # raises ValueError if a default moltype -- with no
    # degen characters -- is used
    data = {
        "s1": "AAA CCC GGG TTT".replace(" ", ""),
        "s2": "CCC GGG T-T AAA".replace(" ", ""),
        "s3": "GGR YTT AAA CCC".replace(" ", ""),
    }

    aln = cls(data=data, moltype=moltype)
    with pytest.raises(ValueError):
        aln.no_degenerates()


@pytest.mark.parametrize("cls", (Alignment, ArrayAlignment))
@pytest.mark.parametrize("calc", ("hamming", None))
def test_quick_tree(cls, calc, brca1_data):
    """quick tree method returns tree"""
    aln = cls(brca1_data, moltype=DNA)[:100]
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


@pytest.mark.parametrize("raw", ("-AAAGGGGGAACCCT", "AAAGGGGGAACCCT"))
def test_slice_aligned(raw):
    imap, seq = DNA.make_seq(seq=raw, name="x").parse_out_gaps()
    al = Aligned(imap, seq)
    sliced = al[:-3]
    assert str(sliced) == raw[:-3]


def test_slice_aligned_featuremap_allgap():
    from cogent3.core.location import FeatureMap, LostSpan

    imap, seq = DNA.make_seq(seq="AAAGGGGGAACCCT", name="x").parse_out_gaps()
    al = Aligned(imap, seq)
    fmap = FeatureMap(spans=[LostSpan(4)], parent_length=0)
    sliced = al[fmap]
    assert not sliced


def test_slice_aligned_featuremap_multi_spans():
    from cogent3.core.location import FeatureMap

    #                    1111111
    #          01234567890123456
    #           ***   **    ***
    raw_seq = "AAAGG--GGG-AACCCT"
    #          01234  567 890123
    #                       1111
    imap, seq = DNA.make_seq(seq=raw_seq, name="x").parse_out_gaps()
    al = Aligned(imap, seq)
    fmap = FeatureMap.from_locations(
        locations=[(1, 4), (7, 9), (13, 16)], parent_length=len(raw_seq)
    )
    sliced = al[fmap]
    assert str(sliced) == "AAGGGCCC"


def test_sequence_collection_repr():  # ported
    data = {
        "ENSMUSG00000056468": "GCCAGGGGGAAAAGGGAGAA",
        "ENSMUSG00000039616": "GCCCTTCAAATTT",
    }
    seqs = SequenceCollection(data=data, moltype=DNA)
    assert (
        repr(seqs)
        == "2x (ENSMUSG00000039616[GCCCTTCAAA...], ENSMUSG00000056468[GCCAGGGGGA...]) dna seqcollection"
    )

    data = {
        "ENSMUSG00000039616": "GCCCTTCAAATTT",
        "ENSMUSG00000056468": "GCCAGGGGGAAAAGGGAGAA",
    }
    seqs = SequenceCollection(data=data, moltype=DNA)
    assert (
        repr(seqs)
        == "2x (ENSMUSG00000039616[GCCCTTCAAA...], ENSMUSG00000056468[GCCAGGGGGA...]) dna seqcollection"
    )

    data = {
        "a": "TCGAT",
    }
    seqs = SequenceCollection(data=data, moltype=DNA)
    assert repr(seqs) == "1x (a[TCGAT]) dna seqcollection"

    data = {
        "a": "TCGAT" * 2,
    }
    seqs = SequenceCollection(data=data, moltype=DNA)
    assert repr(seqs) == "1x (a[TCGATTCGAT]) dna seqcollection"

    data = {
        "a": "A" * 11,
        "b": "B" * 3,
        "c": "C" * 3,
        "d": "D" * 11,
        "e": "E" * 8,
    }
    seqs = SequenceCollection(data=data, moltype=ASCII)
    assert repr(seqs) == "5x (b[BBB], ..., d[DDDDDDDDDD...]) text seqcollection"


@pytest.mark.parametrize("cls", (ArrayAlignment, Alignment))
def test_alignment_repr(cls):
    data = {
        "ENSMUSG00000056468": "GCCAGGGGGAAAA",
        "ENSMUSG00000039616": "GCCCTTCAAATTT",
    }
    seqs = cls(data=data, moltype=DNA)
    assert (
        repr(seqs)
        == "2 x 13 dna alignment: ENSMUSG00000056468[GCCAGGGGGA...], ENSMUSG00000039616[GCCCTTCAAA...]"
    )

    data = {
        "ENSMUSG00000039616": "GCCCTTCAAATTT",
        "ENSMUSG00000056468": "GCCAGGGGGAAAA",
    }
    seqs = cls(data=data, moltype=DNA)
    assert (
        repr(seqs)
        == "2 x 13 dna alignment: ENSMUSG00000039616[GCCCTTCAAA...], ENSMUSG00000056468[GCCAGGGGGA...]"
    )

    data = {
        "a": "TCGAT",
    }
    seqs = cls(data=data, moltype=DNA)
    assert repr(seqs) == "1 x 5 dna alignment: a[TCGAT]"

    data = {
        "a": "TCGAT" * 2,
    }
    seqs = cls(data=data, moltype=DNA)
    assert repr(seqs) == "1 x 10 dna alignment: a[TCGATTCGAT]"

    data = {
        "a": "A" * 11,
        "b": "B" * 11,
        "c": "C" * 11,
        "d": "D" * 11,
        "e": "E" * 11,
    }
    seqs = cls(data=data, moltype=ASCII)
    assert (
        repr(seqs)
        == "5 x 11 text alignment: a[AAAAAAAAAA...], b[BBBBBBBBBB...], c[CCCCCCCCCC...], ..."
    )


@pytest.mark.parametrize("cls", (ArrayAlignment, Alignment, SequenceCollection))
def test_empty_data(cls):
    with pytest.raises(ValueError) as e:
        _ = cls(())
    assert str(e.value) == f"{cls.__name__} must take at least one sequence."


@pytest.mark.parametrize("cls", (ArrayAlignment, Alignment))
def test_coevolution(cls):
    """correctly produces matrix of coevo measures"""
    data = {"s0": "AA", "s1": "AA", "s2": "GG", "s3": "GG", "s4": "GC"}
    aln = cls(data=data, moltype=DNA)
    coevo = aln.coevolution(stat="rmi", show_progress=False)
    expect = array([[nan, nan], [0.78333333, nan]])
    assert_allclose(coevo.array, expect)
    coevo = aln.coevolution(stat="nmi", show_progress=False)
    assert coevo[1, 1] != expect[1, 1]
    # now check invoking drawable produces a result object with a drawable
    # attribute
    coevo = aln.coevolution(stat="nmi", drawable="box", show_progress=False)
    assert hasattr(coevo, "drawable")
    aln = load_aligned_seqs("data/brca1.fasta", moltype="dna")
    aln = aln.take_seqs(aln.names[:20])
    aln = aln.no_degenerates()[:20]
