#!/usr/bin/env python
import json
import os
import pathlib
import re

from os import remove
from tempfile import TemporaryDirectory, mktemp
from unittest import TestCase, main

import numpy

from numpy import array, log2, nan, transpose
from numpy.testing import assert_allclose, assert_equal

from cogent3 import load_aligned_seqs, load_unaligned_seqs, make_seq, open_
from cogent3.core.alignment import (
    Aligned,
    Alignment,
    ArrayAlignment,
    DataError,
    SequenceCollection,
    _SequenceCollectionBase,
    aln_from_array,
    aln_from_array_aln,
    aln_from_array_seqs,
    aln_from_collection,
    aln_from_empty,
    aln_from_fasta,
    aln_from_generic,
    make_gap_filter,
    seqs_from_aln,
    seqs_from_array,
    seqs_from_array_seqs,
    seqs_from_empty,
    seqs_from_fasta,
    seqs_from_generic,
    seqs_from_kv_pairs,
)
from cogent3.core.alphabet import AlphabetError
from cogent3.core.annotation import Feature
from cogent3.core.moltype import AB, ASCII, BYTES, DNA, PROTEIN, RNA
from cogent3.core.sequence import ArraySequence, RnaSequence, Sequence
from cogent3.maths.util import safe_p_log_p
from cogent3.parse.fasta import MinimalFastaParser
from cogent3.util.misc import get_object_provenance


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Jeremy Widmann",
    "Catherine Lozuopone",
    "Gavin Huttley",
    "Rob Knight",
    "Daniel McDonald",
    "Jan Kosinski",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class alignment_tests(TestCase):
    """Tests of top-level functions."""

    def test_seqs_from_array(self):
        """seqs_from_array should return chars, and successive indices."""
        a = array([[0, 1, 2], [2, 1, 0]])  # three 2-char seqs
        obs_a, obs_labels = seqs_from_array(a)
        # note transposition
        assert_equal(obs_a, [array([0, 2]), array([1, 1]), array([2, 0])])
        self.assertEqual(obs_labels, None)

    def test_seqs_from_array_seqs(self):
        """seqs_from_array_seqs should return model seqs + names."""
        s1 = ArraySequence("ABC", name="a")
        s2 = ArraySequence("DEF", name="b")
        obs_a, obs_labels = seqs_from_array_seqs([s1, s2])
        self.assertEqual(obs_a, [s1, s2])  # seq -> numbers
        self.assertEqual(obs_labels, ["a", "b"])

    def test_seqs_from_generic(self):
        """seqs_from_generic should initialize seqs from list of lists, etc."""
        s1 = "ABC"
        s2 = "DEF"
        obs_a, obs_labels = seqs_from_generic([s1, s2])
        self.assertEqual(obs_a, ["ABC", "DEF"])
        self.assertEqual(obs_labels, [None, None])

    def test_seqs_from_fasta(self):
        """seqs_from_fasta should initialize seqs from fasta-format string"""
        s = ">aa\nAB\nC\n>bb\nDE\nF\n"
        obs_a, obs_labels = seqs_from_fasta(s)
        self.assertEqual(obs_a, ["ABC", "DEF"])
        self.assertEqual(obs_labels, ["aa", "bb"])

    def test_seqs_from_aln(self):
        """seqs_from_aln should initialize from existing alignment"""
        c = SequenceCollection(["abc", "def"])
        obs_a, obs_labels = seqs_from_aln(c)
        self.assertEqual(obs_a, ["abc", "def"])
        self.assertEqual(obs_labels, ["seq_0", "seq_1"])

    def test_seqs_from_kv_pairs(self):
        """seqs_from_kv_pairs should initialize from key-value pairs"""
        c = [["a", "abc"], ["b", "def"]]
        obs_a, obs_labels = seqs_from_kv_pairs(c)
        self.assertEqual(obs_a, ["abc", "def"])
        self.assertEqual(obs_labels, ["a", "b"])

    def test_seqs_from_empty(self):
        """seqs_from_empty should always raise ValueError"""
        self.assertRaises(ValueError, seqs_from_empty, "xyz")

    def test_aln_from_array(self):
        """aln_from_array should return same array, and successive indices."""
        a = array([[0, 1, 2], [3, 4, 5]])  # three 2-char seqs
        obs_a, obs_labels = aln_from_array(a)
        assert_equal(obs_a, transpose(a))
        assert_equal(obs_labels, None)

    def test_aln_from_array_seqs(self):
        """aln_from_array_seqs should initialize aln from sequence objects."""
        s1 = ArraySequence("ACC", name="a", alphabet=RNA.alphabet)
        s2 = ArraySequence("GGU", name="b", alphabet=RNA.alphabet)
        obs_a, obs_labels = aln_from_array_seqs([s1, s2], alphabet=BYTES.alphabet)
        assert_equal(obs_a, array([[2, 1, 1], [3, 3, 0]], "b"))
        # seq -> numbers
        assert_equal(obs_labels, ["a", "b"])

    def test_aln_from_generic(self):
        """aln_from_generic should initialize aln from list of lists, etc."""
        s1 = "AAA"
        s2 = "GGG"
        obs_a, obs_labels = aln_from_generic(
            [s1, s2], "b", alphabet=RNA.alphabet
        )  # specify array type
        assert_equal(obs_a, array([[2, 2, 2], [3, 3, 3]], "b"))  # str -> chars
        assert_equal(obs_labels, [None, None])

    def test_aln_from_fasta(self):
        """aln_from_fasta should initialize aln from fasta-format string"""
        s = ">aa\nAB\nC\n>bb\nDE\nF\n"
        obs_a, obs_labels = aln_from_fasta(s.splitlines())
        assert_equal(obs_a, array(["ABC", "DEF"], "c").view("B"))  # seq -> numbers
        assert_equal(obs_labels, ["aa", "bb"])

    def test_aln_from_array_aln(self):
        """aln_from_array_aln should initialize from existing alignment"""
        a = ArrayAlignment(array([[0, 1, 2], [3, 4, 5]]), conversion_f=aln_from_array)
        obs_a, obs_labels = aln_from_array_aln(a)
        assert_equal(obs_a, a.seq_data)
        assert_equal(obs_labels, a.names)

    def test_aln_from_collection(self):
        """aln_from_collection should initialize from existing alignment"""
        a = SequenceCollection(["AAA", "GGG"])
        obs_a, obs_labels = aln_from_collection(a, alphabet=RNA.alphabet)
        assert_equal(a.to_fasta(), ">seq_0\nAAA\n>seq_1\nGGG\n")
        assert_equal(obs_a, array([[2, 2, 2], [3, 3, 3]]))

    def test_aln_from_empty(self):
        """aln_from_empty should always raise ValueError"""
        self.assertRaises(ValueError, aln_from_empty, "xyz")


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
                "a": RnaSequence("AAAAAAA"),
                "b": RnaSequence("A--A-AA"),
                "c": RnaSequence("AA-----"),
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
                "a": RnaSequence("UCAGUCAGUU"),
                "b": RnaSequence("UCCGUCAAUU"),
                "c": RnaSequence("ACCAUCAGUC"),
                "d": RnaSequence("UCAAUCGGUU"),
                "e": RnaSequence("UUGGUUGGGU"),
                "f": RnaSequence("CCGGGCGGCC"),
                "g": RnaSequence("UCAACCGGAA"),
            }
        )
        # Additional SequenceCollections for tests added 6/4/04 by Jeremy
        # Widmann
        self.sequences = self.Class(list(map(RnaSequence, ["UCAG", "UCAG", "UCAG"])))
        # Additional SequenceCollection for tests added 1/30/06 by Cathy
        # Lozupone
        self.omitSeqsTemplate_aln = self.Class(
            {
                "s1": RnaSequence("UC-----CU---C"),
                "s2": RnaSequence("UC------U---C"),
                "s3": RnaSequence("UUCCUUCUU-UUC"),
                "s4": RnaSequence("UU-UUUU-UUUUC"),
                "s5": RnaSequence("-------------"),
            }
        )

        self.a = ArrayAlignment(["AAA", "AAA"])
        self.b = Alignment(["AAA", "AAA"])
        self.c = SequenceCollection(["AAA", "AAA"])

    def test_deepcopy(self):
        """correctly deep copy aligned objects in an alignment"""
        data = {"seq1": "ACGACGACG", "seq2": "ACGACGACG"}
        seqs = self.Class(data)
        copied = seqs.deepcopy(sliced=True)
        assert_equal(seqs.to_rich_dict(), copied.to_rich_dict())
        self.assertNotEqual(id(copied), id(seqs))
        for name in seqs.names:
            self.assertNotEqual(id(copied.named_seqs[name]), copied.named_seqs[name])

    def test_guess_input_type(self):
        """SequenceCollection  _guess_input_type should figure out data type correctly"""
        git = self.a._guess_input_type
        self.assertEqual(git(self.a), "array_aln")
        self.assertEqual(git(self.b), "aln")
        self.assertEqual(git(self.c), "collection")
        self.assertEqual(git(">ab\nabc"), "fasta")
        self.assertEqual(git([">ab", "abc"]), "fasta")
        self.assertEqual(git(["abc", "def"]), "generic")
        # precedence over generic
        self.assertEqual(git([[1, 2], [4, 5]]), "kv_pairs")
        self.assertEqual(git([[1, 2, 3], [4, 5, 6]]), "generic")
        self.assertEqual(git([ArraySequence("abc")]), "array_seqs")
        self.assertEqual(git(array([[1, 2, 3], [4, 5, 6]])), "array")
        self.assertEqual(git({"a": "aca"}), "dict")
        self.assertEqual(git([]), "empty")

    def test_init_aln(self):
        """SequenceCollection should init from existing alignments"""
        exp = self.Class(["AAA", "AAA"])
        x = self.Class(self.a)
        y = self.Class(self.b)
        z = self.Class(self.c)
        assert_equal(x, exp)
        assert_equal(z, exp)
        assert_equal(y, exp)

    test_init_aln.__doc__ = Class.__name__ + test_init_aln.__doc__

    def test_init_dict(self):
        """SequenceCollection init from dict should work as expected"""
        d = {"a": "AAAAA", "b": "BBBBB"}
        a = self.Class(d, names=["a", "b"])
        self.assertEqual(a, d)
        self.assertEqual(list(a.named_seqs.items()), list(d.items()))

        # from bytes strings
        a = self.Class({"a": b"AAAAA", "b": b"BBBBB"}, names=["a", "b"])
        self.assertEqual(a, d)
        self.assertEqual(list(a.named_seqs.items()), list(d.items()))

    def test_names_attribute(self):
        """expected to be a list"""
        seqs = self.Class({"a": b"AAAAA", "b": b"BBBBB"}, names=("a", "b"))
        self.assertIsInstance(seqs.names, list)

    def test_init_name_mapped(self):
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

    def test_init_seq(self):
        """SequenceCollection init from list of sequences should use indices as keys"""
        seqs = ["AAAAA", "BBBBB", "CCCCC"]
        a = self.Class(seqs)
        self.assertEqual(len(a.named_seqs), 3)
        self.assertEqual(a.named_seqs["seq_0"], "AAAAA")
        self.assertEqual(a.named_seqs["seq_1"], "BBBBB")
        self.assertEqual(a.named_seqs["seq_2"], "CCCCC")
        self.assertEqual(a.names, ["seq_0", "seq_1", "seq_2"])
        self.assertEqual(list(a.seqs), ["AAAAA", "BBBBB", "CCCCC"])

    def test_init_seq_info(self):
        """SequenceCollection init from seqs w/ info should preserve data"""
        a = Sequence("AAA", name="a", info={"x": 3})
        b = Sequence("CCC", name="b", info={"x": 4})
        c = Sequence("GGG", name="c", info={"x": 5})
        seqs = [c, b, a]
        a = self.Class(seqs)
        self.assertEqual(list(a.names), ["c", "b", "a"])
        self.assertEqual(list(map(str, a.seqs)), ["GGG", "CCC", "AAA"])
        if self.Class is not ArrayAlignment:
            # ArrayAlignment is allowed to strip info objects
            self.assertEqual([i.info.x for i in a.seqs], [5, 4, 3])
        # check it still works if constructed from same class
        b = self.Class(a)
        self.assertEqual(list(b.names), ["c", "b", "a"])
        self.assertEqual(list(map(str, b.seqs)), ["GGG", "CCC", "AAA"])
        if self.Class is not ArrayAlignment:
            # ArrayAlignment is allowed to strip Info objects
            self.assertEqual([i.info.x for i in b.seqs], [5, 4, 3])

    def test_init_annotated_seqs(self):
        """correctly construct from list with annotated seq"""
        if self.Class == ArrayAlignment:
            # this class cannot be annotated
            return
        seq = make_seq("GCCAGGGGGGAAAG-GGAGAA", name="seq1")
        _ = seq.add_feature("exon", "name", [(4, 10)])
        coll = self.Class(data=[seq])
        got_seq = coll.get_seq("seq1")
        ann = got_seq.annotations[0]
        self.assertEqual(str(got_seq[ann]), "GGGGGG")

    def test_init_pairs(self):
        """SequenceCollection init from list of (key,val) pairs should work correctly"""
        seqs = [["x", "XXX"], ["b", "BBB"], ["c", "CCC"]]
        a = self.Class(seqs)
        self.assertEqual(len(a.named_seqs), 3)
        self.assertEqual(a.named_seqs["x"], "XXX")
        self.assertEqual(a.named_seqs["b"], "BBB")
        self.assertEqual(a.named_seqs["c"], "CCC")
        self.assertEqual(a.names, ["x", "b", "c"])
        self.assertEqual(list(a.seqs), ["XXX", "BBB", "CCC"])

    def test_init_duplicate_keys(self):
        """SequenceCollection init from (key, val) pairs should fail on dup. keys"""
        seqs = [["x", "XXX"], ["b", "BBB"], ["x", "CCC"], ["d", "DDD"], ["a", "AAA"]]
        self.assertRaises(ValueError, self.Class, seqs)
        aln = self.Class(seqs, remove_duplicate_names=True)
        self.assertEqual(
            str(self.Class(seqs, remove_duplicate_names=True)),
            ">x\nXXX\n>b\nBBB\n>d\nDDD\n>a\nAAA\n",
        )

    def test_init_ordered(self):
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

    def test_init_ambig(self):
        """SequenceCollection should tolerate ambiguous chars"""
        aln = self.Class(["AAA", "CCC"], moltype=DNA)
        aln = self.Class(["ANS", "CWC"], moltype=DNA)
        aln = self.Class(["A-A", "CC-"], moltype=DNA)
        aln = self.Class(["A?A", "CC-"], moltype=DNA)

    def test_aln_from_fasta_parser(self):
        """aln_from_fasta_parser should init from iterator"""
        s = ">aa\nAC\n>bb\nAA\n>c\nGG\n".splitlines()
        p = MinimalFastaParser(s)
        aln = self.Class(p, moltype=DNA)
        self.assertEqual(aln.named_seqs["aa"], "AC")
        self.assertEqual(aln.to_fasta(), ">aa\nAC\n>bb\nAA\n>c\nGG\n")
        s2_ORIG = ">x\nCA\n>b\nAA\n>>xx\nGG"
        s2 = ">aa\nAC\n>bb\nAA\n>c\nGG\n"
        d = ArrayAlignment(MinimalFastaParser(s2.splitlines()))
        d.to_fasta()
        self.assertEqual(d.to_fasta(), aln.to_fasta())

    def test_aln_from_fasta(self):
        """SequenceCollection should init from fasta-format string"""
        s = ">aa\nAC\n>bb\nAA\n>c\nGG\n"
        aln = self.Class(s)
        self.assertEqual(aln.to_fasta(), s)

    def test_seq_len_get(self):
        """SequenceCollection seq_len should return length of longest seq"""
        self.assertEqual(self.one_seq.seq_len, 5)
        self.assertEqual(self.identical.seq_len, 4)
        self.assertEqual(self.gaps.seq_len, 7)

    def test_Seqs(self):
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

    def test_iter_seqs(self):
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

    def test_Items(self):
        """SequenceCollection iter_selected should iterate over items in specified order."""
        # should work if one row
        self.assertEqual(list(self.one_seq.iter_selected()), ["A"] * 5)
        # should take order into account
        self.assertEqual(list(self.ordered1.iter_selected()), ["A"] * 5 + ["B"] * 5)
        self.assertEqual(list(self.ordered2.iter_selected()), ["B"] * 5 + ["A"] * 5)

    def test_iter_selected(self):
        """SequenceCollection iter_selected() should iterate over items in correct order"""
        # should work if one row
        self.assertEqual(list(self.one_seq.iter_selected()), ["A"] * 5)
        # should take order into account
        self.assertEqual(list(self.ordered1.iter_selected()), ["A"] * 5 + ["B"] * 5)
        self.assertEqual(list(self.ordered2.iter_selected()), ["B"] * 5 + ["A"] * 5)
        # should allow row and/or col specification
        r = self.ragged_padded
        self.assertEqual(
            list(r.iter_selected(seq_order=["c", "b"], pos_order=[5, 1, 3])),
            list("-AA-A-"),
        )
        # should not interfere with superclass iteritems()
        i = list(r.named_seqs.items())
        i.sort()
        self.assertEqual(i, [("a", "AAAAAA"), ("b", "AAA---"), ("c", "AAAA--")])

    def test_take_seqs(self):
        """SequenceCollection take_seqs should return new SequenceCollection with selected seqs."""
        a = self.ragged_padded.take_seqs(list("bc"))
        self.assertTrue(isinstance(a, _SequenceCollectionBase))
        self.assertEqual(a, {"b": "AAA---", "c": "AAAA--"})
        # should be able to negate
        a = self.ragged_padded.take_seqs(list("bc"), negate=True)
        self.assertEqual(a, {"a": "AAAAAA"})

    def test_take_seqs_str(self):
        """string arg to SequenceCollection take_seqs should work."""
        a = self.ragged_padded.take_seqs("a", negate=True)
        self.assertTrue(isinstance(a, _SequenceCollectionBase))
        self.assertEqual(a, {"b": "AAA---", "c": "AAAA--"})
        # should be able to negate
        a = self.ragged_padded.take_seqs("a")
        self.assertEqual(a, {"a": "AAAAAA"})

    def test_take_seqs_info(self):
        """take_seqs should preserve info attribute"""
        orig = self.Class(
            data={"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"}, info={"key": "value"}
        )
        subset = orig.take_seqs(list("ab"))
        self.assertEqual(set(subset.info), set(orig.info))

    def test_take_seqs_moltype(self):
        """take_seqs should preserve the MolType"""
        orig = self.Class(
            data={"a": "CCCCCC", "b": "AAA---", "c": "AAAA--"}, moltype=DNA
        )
        subset = orig.take_seqs(list("ab"))
        self.assertEqual(set(subset.moltype), set(orig.moltype))

    def test_get_seq_indices(self):
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

    def test_take_seqs_if(self):
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

    def test_get_identical_sets(self):
        """correctly identify sets of identical sequences"""
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

        got = seqs.get_identical_sets(mask_degen=True)
        got = frozenset(frozenset(s) for s in got)
        self.assertEqual(got, expect)

    def test_get_similar(self):
        """SequenceCollection get_similar should get all sequences close to target seq"""
        aln = self.many
        RnaSequence("GGGGGGGGGG")
        RnaSequence("----------")
        # test min and max similarity ranges
        result = aln.get_similar(
            aln.named_seqs["a"], min_similarity=0.4, max_similarity=0.7
        )
        for seq in "cefg":
            self.assertIn(seq, result.named_seqs)
            self.assertEqual(result.named_seqs[seq], aln.named_seqs[seq])
        self.assertEqual(len(result.named_seqs), 4)

        result = aln.get_similar(
            aln.named_seqs["a"], min_similarity=0.95, max_similarity=1
        )
        for seq in "a":
            self.assertIn(seq, result.named_seqs)
            self.assertEqual(result.named_seqs[seq], aln.named_seqs[seq])
        self.assertEqual(len(result.named_seqs), 1)

        result = aln.get_similar(
            aln.named_seqs["a"], min_similarity=0.75, max_similarity=0.85
        )
        for seq in "bd":
            self.assertIn(seq, result.named_seqs)
            self.assertEqual(result.named_seqs[seq], aln.named_seqs[seq])
        self.assertEqual(len(result.named_seqs), 2)

        result = aln.get_similar(
            aln.named_seqs["a"], min_similarity=0, max_similarity=0.2
        )
        self.assertEqual(result, {})

        # test some sequence transformations
        def transform(s):
            return s[1:4]

        result = aln.get_similar(
            aln.named_seqs["a"], min_similarity=0.5, transform=transform
        )
        for seq in "abdfg":
            self.assertIn(seq, result.named_seqs)
            self.assertEqual(result.named_seqs[seq], aln.named_seqs[seq])
        self.assertEqual(len(result.named_seqs), 5)

        def transform(s):
            return s[-3:]

        result = aln.get_similar(
            aln.named_seqs["a"], min_similarity=0.5, transform=transform
        )
        for seq in "abcde":
            self.assertIn(seq, result.named_seqs)
            self.assertEqual(result.named_seqs[seq], aln.named_seqs[seq])
        self.assertEqual(len(result.named_seqs), 5)

        # test a different distance metric
        def metric(x, y):
            return str(x).count("G") + str(y).count("G")

        result = aln.get_similar(
            aln.named_seqs["a"], min_similarity=5, max_similarity=10, metric=metric
        )
        for seq in "ef":
            self.assertIn(seq, result.named_seqs)
            self.assertEqual(result.named_seqs[seq], aln.named_seqs[seq])
        self.assertEqual(len(result.named_seqs), 2)

        # test the combination of a transform and a distance metric
        aln = self.Class(
            dict(enumerate(map(RnaSequence, ["GA-GU", "A-GAC", "GG-GG"]))), moltype=RNA
        )

        def transform(s):
            return RnaSequence(str(s).replace("G", "A").replace("U", "C"))

        metric = RnaSequence.frac_same_non_gaps

        def null_transform(s):
            return RnaSequence(str(s))

        # first, do it without the transformation
        try:
            result = aln.get_similar(
                aln.named_seqs[0], min_similarity=0.5, metric=metric
            )
        except TypeError:  # need to coerce to RNA seq w/ null_transform
            result = aln.get_similar(
                aln.named_seqs[0],
                min_similarity=0.5,
                metric=metric,
                transform=null_transform,
            )
        for seq in [0, 2]:
            self.assertIn(seq, result.named_seqs)
            self.assertEqual(result.named_seqs[seq], aln.named_seqs[seq])
        self.assertEqual(len(result.named_seqs), 2)
        # repeat with higher similarity
        try:
            result = aln.get_similar(
                aln.named_seqs[0], min_similarity=0.8, metric=metric
            )
        except TypeError:  # need to coerce to RNA
            result = aln.get_similar(
                aln.named_seqs[0],
                min_similarity=0.8,
                metric=metric,
                transform=null_transform,
            )
        for seq in [0]:
            self.assertIn(seq, result.named_seqs)
            self.assertEqual(result.named_seqs[seq], aln.named_seqs[seq])
        self.assertEqual(len(result.named_seqs), 1)
        # then, verify that the transform changes the results
        result = aln.get_similar(
            aln.named_seqs[0], min_similarity=0.5, metric=metric, transform=transform
        )
        for seq in [0, 1, 2]:
            self.assertIn(seq, result.named_seqs)
            self.assertEqual(result.named_seqs[seq], aln.named_seqs[seq])
        self.assertEqual(len(result.named_seqs), 3)

        result = aln.get_similar(
            aln.named_seqs[0], min_similarity=0.8, metric=metric, transform=transform
        )
        for seq in [0, 1]:
            self.assertIn(seq, result.named_seqs)
            self.assertEqual(result.named_seqs[seq], aln.named_seqs[seq])
        self.assertEqual(len(result.named_seqs), 2)

    def test_is_ragged(self):
        """SequenceCollection is_ragged should return true if ragged alignment"""
        assert not self.identical.is_ragged()
        assert not self.gaps.is_ragged()

    def test_to_phylip(self):
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

    def test_to_fasta(self):
        """SequenceCollection should return correct FASTA string"""
        aln = self.Class(["AAA", "CCC"])
        self.assertEqual(aln.to_fasta(), ">seq_0\nAAA\n>seq_1\nCCC\n")

        # NOTE THE FOLLOWING SURPRISING BEHAVIOR BECAUSE OF THE TWO-ITEM
        # SEQUENCE RULE:
        aln = self.Class(["AA", "CC"])
        self.assertEqual(aln.to_fasta(), ">A\nA\n>C\nC\n")

    def test_to_nexus(self):
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

    def test_to_rich_dict(self):
        """to_rich_dict produces correct dict"""
        aln = self.Class({"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"})
        try:
            seq_type = get_object_provenance(aln.seqs[0].data)
        except AttributeError:
            seq_type = get_object_provenance(aln.seqs[0])

        got = aln.to_rich_dict()
        expect = {
            "seqs": {
                "seq1": {
                    "name": "seq1",
                    "seq": "ACGG",
                    "info": None,
                    "type": seq_type,
                    "moltype": aln.moltype.label,
                    "version": __version__,
                },
                "seq2": {
                    "name": "seq2",
                    "seq": "CGCA",
                    "info": None,
                    "type": seq_type,
                    "moltype": aln.moltype.label,
                    "version": __version__,
                },
                "seq3": {
                    "name": "seq3",
                    "seq": "CCG-",
                    "info": None,
                    "type": seq_type,
                    "moltype": aln.moltype.label,
                    "version": __version__,
                },
            },
            "moltype": aln.moltype.label,
            "info": None,
            "type": get_object_provenance(aln),
            "version": __version__,
        }
        self.assertEqual(got, expect)

    def test_to_json(self):
        """roundtrip of to_json produces correct dict"""
        aln = self.Class({"seq1": "ACGG", "seq2": "CGCA", "seq3": "CCG-"})
        got = json.loads(aln.to_json())
        expect = aln.to_rich_dict()
        self.assertEqual(got, expect)

    def test_num_seqs(self):
        """SequenceCollection.num_seqs should count seqs."""
        aln = self.Class({"seq1": "ACGU", "seq2": "CGUA", "seq3": "CCGU"})
        self.assertEqual(aln.num_seqs, 3)

    def test_copy_annotations(self):
        """SequenceCollection copy_annotations should copy from seq objects"""
        if self.Class == ArrayAlignment:
            return
        aln = self.Class({"seq1": "ACGU", "seq2": "CGUA", "seq3": "CCGU"})
        seq_1 = Sequence("ACGU", name="seq1")
        seq_1.add_feature("xyz", "abc", [(1, 2)])
        seq_5 = Sequence("ACGUAAAAAA", name="seq5")
        seq_5.add_feature("xyzzz", "abc", [(1, 2)])
        annot = {"seq1": seq_1, "seq5": seq_5}
        aln.copy_annotations(annot)
        aln_seq_1 = aln.named_seqs["seq1"]
        if not hasattr(aln_seq_1, "annotations"):
            aln_seq_1 = aln_seq_1.data
        aln_seq_2 = aln.named_seqs["seq2"]
        if not hasattr(aln_seq_2, "annotations"):
            aln_seq_2 = aln_seq_2.data
        self.assertEqual(len(aln_seq_1.annotations), 1)
        self.assertEqual(aln_seq_1.annotations[0].name, "abc")
        self.assertEqual(len(aln_seq_2.annotations), 0)

    def test_annotate_from_gff(self):
        """SequenceCollection.annotate_from_gff should read gff features"""
        aln = self.Class({"seq1": "ACGU", "seq2": "CGUA", "seq3": "C-GU"})
        gff = [
            ["seq1", "prog1", "snp", "1", "2", "1.0", "+", "1", '"abc"'],
            ["seq3", "prog2", "del", "1", "3", "1.0", "+", "1", '"xyz"'],
            ["seq5", "prog2", "snp", "2", "3", "1.0", "+", "1", '"yyy"'],
        ]
        gff = list(map("\t".join, gff))
        if self.Class == ArrayAlignment:
            return

        aln.annotate_from_gff(gff)
        aln_seq_1 = aln.named_seqs["seq1"]
        if not hasattr(aln_seq_1, "annotations"):
            aln_seq_1 = aln_seq_1.data
        aln_seq_2 = aln.named_seqs["seq2"]
        if not hasattr(aln_seq_2, "annotations"):
            aln_seq_2 = aln_seq_2.data
        self.assertEqual(len(aln_seq_1.annotations), 1)
        self.assertEqual(aln_seq_1.annotations[0].name, "abc")
        self.assertEqual(len(aln_seq_2.annotations), 0)

        if self.Class == Alignment:
            aln_seq_3 = aln.get_seq("seq3")
            matches = [m for m in aln_seq_3.get_annotations_matching("*")]
            self.assertFalse("-" in matches[0].get_slice())

    def test_annotate_from_gff3(self):
        """annotate_from_gff should work on data from gff3 files"""
        from cogent3.parse.fasta import FastaParser

        if self.Class == ArrayAlignment:
            return

        fasta_path = os.path.join("data/c_elegans_WS199_dna_shortened.fasta")
        gff3_path = os.path.join("data/c_elegans_WS199_shortened_gff.gff3")
        name, seq = next(FastaParser(fasta_path))

        # using annotate_from_gff will nest annotations
        aln = self.Class({name: seq})
        aln.annotate_from_gff(gff3_path)
        aln_seq = aln.named_seqs[name]
        if not hasattr(aln_seq, "annotations"):
            aln_seq = aln_seq.data
        matches = [m for m in aln_seq.get_annotations_matching("*", extend_query=True)]
        # 13 features with one having 2 parents, so 14 instances should be found
        self.assertEqual(len(matches), 14)
        matches = [m for m in aln_seq.get_annotations_matching("gene")]
        self.assertEqual(len(matches), 1)
        matches = matches[0].get_annotations_matching("mRNA")
        self.assertEqual(len(matches), 1)
        matches = matches[0].get_annotations_matching("exon")
        self.assertEqual(len(matches), 3)

    def test_add(self):
        """__add__ should concatenate sequence data, by name"""
        align1 = self.Class({"a": "AAAA", "b": "TTTT", "c": "CCCC"})
        align2 = self.Class({"a": "GGGG", "b": "----", "c": "NNNN"})
        align = align1 + align2
        concatdict = align.to_dict()
        self.assertEqual(
            concatdict, {"a": "AAAAGGGG", "b": "TTTT----", "c": "CCCCNNNN"}
        )

    def test_add_info(self):
        """__add__ should preserve info attribute"""
        align1 = self.Class(
            {"a": "AAAA", "b": "TTTT", "c": "CCCC"}, info={"key": "foo"}
        )
        align2 = self.Class(
            {"a": "GGGG", "b": "----", "c": "NNNN"}, info={"key": "bar"}
        )
        align = align1 + align2
        self.assertEqual(align.info["key"], "foo")

    def test_add_seqs(self):
        """add_seqs should return an alignment with the new sequences appended or inserted"""
        data = [("name1", "AAA"), ("name2", "AAA"), ("name3", "AAA"), ("name4", "AAA")]
        data1 = [("name1", "AAA"), ("name2", "AAA")]
        data2 = [("name3", "AAA"), ("name4", "AAA")]
        data3 = [("name5", "BBB"), ("name6", "CCC")]
        aln = self.Class(data)
        aln3 = self.Class(data3)

        out_aln = aln.add_seqs(aln3)
        # test append at the end
        self.assertEqual(str(out_aln), str(self.Class(data + data3)))

        out_aln = aln.add_seqs(aln3, before_name="name3")
        self.assertEqual(
            str(out_aln), str(self.Class(data1 + data3 + data2))
        )  # test insert before

        out_aln = aln.add_seqs(aln3, after_name="name2")
        self.assertEqual(
            str(out_aln), str(self.Class(data1 + data3 + data2))
        )  # test insert after

        out_aln = aln.add_seqs(aln3, before_name="name1")
        # test if insert before first seq works
        self.assertEqual(str(out_aln), str(self.Class(data3 + data)))

        out_aln = aln.add_seqs(aln3, after_name="name4")
        # test if insert after last seq works
        self.assertEqual(str(out_aln), str(self.Class(data + data3)))

        self.assertRaises(
            ValueError, aln.add_seqs, aln3, before_name="name5"
        )  # wrong after/before name
        self.assertRaises(
            ValueError, aln.add_seqs, aln3, after_name="name5"
        )  # wrong after/before name

        if isinstance(aln, Alignment) or isinstance(aln, ArrayAlignment):
            self.assertRaises((DataError, ValueError), aln.add_seqs, aln3 + aln3)
        else:
            exp = set([seq for name, seq in data])
            exp.update([seq + seq for name, seq in data3])
            got = set()
            for seq in aln.add_seqs(aln3 + aln3).seqs:
                got.update([str(seq).strip()])
            self.assertEqual(got, exp)

    def test_add_seqs_info(self):
        """add_seqs should preserve info attribute"""
        data = [("name1", "AAA"), ("name2", "AAA"), ("name3", "AAA"), ("name4", "AAA")]
        data2 = [("name5", "BBB"), ("name6", "CCC")]
        aln = self.Class(data, info={"key": "foo"})
        aln2 = self.Class(data2, info={"key": "bar"})
        out_aln = aln.add_seqs(aln2)
        self.assertEqual(out_aln.info["key"], "foo")

    def test_write(self):
        """SequenceCollection.write should write in correct format"""
        aln = self.Class([("a", "AAAA"), ("b", "TTTT"), ("c", "CCCC")])
        fn = mktemp(suffix=".fasta")
        aln.write(fn)
        with open(fn, newline=None) as infile:
            result = infile.read()
        self.assertEqual(result, ">a\nAAAA\n>b\nTTTT\n>c\nCCCC\n")
        remove(fn)

    def test_len(self):
        """len(SequenceCollection) returns length of longest sequence"""
        aln = self.Class([("a", "AAAA"), ("b", "TTTT"), ("c", "CCCC")])
        self.assertEqual(len(aln), 4)

    def test_get_translation(self):
        """SequenceCollection.get_translation translates each seq"""
        for seqs in [
            {"seq1": "GATTTT", "seq2": "GATC??"},
            {"seq1": "GAT---", "seq2": "?GATCT"},
        ]:
            alignment = self.Class(data=seqs, moltype=DNA)
            got = alignment.get_translation()
            self.assertEqual(len(got), 2)
            self.assertEqual(got.moltype, PROTEIN)
            # check for a failure when no moltype specified
            alignment = self.Class(data=seqs)
            try:
                alignment.get_translation()
            except AttributeError:
                pass

    def test_get_translation_info(self):
        """SequenceCollection.get_translation preserves info attribute"""
        for seqs in [
            {"seq1": "GATTTT", "seq2": "GATC??"},
            {"seq1": "GAT---", "seq2": "?GATCT"},
        ]:
            alignment = self.Class(data=seqs, moltype=DNA, info={"key": "value"})
            got = alignment.get_translation()
            self.assertEqual(got.info["key"], "value")

    def test_get_translation_incomplete(self):
        """get translation works on incomplete codons"""
        alignment = self.Class(data={"seq1": "GATN--", "seq2": "?GATCT"}, moltype=DNA)
        got = alignment.get_translation(incomplete_ok=True)
        self.assertEqual(got.to_dict(), {"seq1": "D?", "seq2": "XS"})
        with self.assertRaises(AlphabetError):
            got = alignment.get_translation(incomplete_ok=False)

    def test_get_seq(self):
        """SequenceCollection.get_seq should return specified seq"""
        aln = self.Class({"seq1": "GATTTT", "seq2": "GATC??"})
        self.assertEqual(aln.get_seq("seq1"), "GATTTT")
        self.assertRaises(KeyError, aln.get_seq, "seqx")

    def test_to_dict(self):
        """SequenceCollection.to_dict should return dict of strings (not obj)"""
        aln = self.Class({"seq1": "GATTTT", "seq2": "GATC??"})
        self.assertEqual(aln.to_dict(), {"seq1": "GATTTT", "seq2": "GATC??"})
        for i in list(aln.to_dict().values()):
            assert isinstance(i, str)

    def test_get_ambiguous_positions(self):
        """SequenceCollection.get_ambiguous_positions should return pos"""
        aln = self.Class({"s1": "ATGRY?", "s2": "T-AG??"}, moltype=DNA)
        self.assertEqual(
            aln.get_ambiguous_positions(),
            {"s2": {4: "?", 5: "?"}, "s1": {3: "R", 4: "Y", 5: "?"}},
        )

    def test_degap(self):
        """SequenceCollection.degap should strip gaps from each seq"""
        aln = self.Class({"s1": "ATGRY?", "s2": "T-AG??"}, moltype=DNA)
        self.assertEqual(aln.degap(), {"s1": "ATGRY", "s2": "TAG"})

    def test_degap_info(self):
        """.degap should preserve info attributes"""
        aln = self.Class({"s1": "ATGRY?", "s2": "T-AG??"}, moltype=DNA)
        aln.info.path = "blah"
        got = aln.degap()
        self.assertEqual(got.info.path, "blah")

    def test_with_modified_termini(self):
        """SequenceCollection.with_modified_termini should code trailing gaps as ?"""
        aln = self.Class({"s1": "AATGR--", "s2": "-T-AG?-"}, moltype=DNA)
        self.assertEqual(
            aln.with_modified_termini(), {"s1": "AATGR??", "s2": "?T-AG??"}
        )

    def test_make_gap_filter(self):
        """make_gap_filter returns f(seq) -> True if aligned ok w/ query"""
        s1 = RnaSequence("UC-----CU---C")
        s3 = RnaSequence("UUCCUUCUU-UUC")
        s4 = RnaSequence("UU-UUUU-UUUUC")
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

    def test_omit_gap_seqs(self):
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

    def test_omit_gap_runs(self):
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

    def test_consistent_gap_degen_handling(self):
        """gap degen character should be treated consistently"""
        # the degen character '?' can be a gap, so when we strip gaps it should
        # be gone too
        raw_seq = "---??-??TC-GGCG-GCA-G-GC-?-C-TAN-GCGC-CCTC-AGGA?-???-??--"
        raw_ungapped = re.sub("[-?]", "", raw_seq)
        re.sub("[N?]+", "", raw_seq)
        dna = DNA.make_seq(raw_seq)

        aln = self.Class(data=[("a", dna), ("b", dna)])
        expect = self.Class(data=[("a", raw_ungapped), ("b", raw_ungapped)]).to_fasta()
        self.assertEqual(aln.degap().to_fasta(), expect)
        seqs = self.Class(data=[("a", dna), ("b", dna)])
        self.assertEqual(seqs.degap().to_fasta(), expect)

    def test_pad_seqs(self):
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

    def test_to_moltype(self):
        """correctly convert to specified moltype"""
        data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        seqs = self.Class(data=data)
        dna = seqs.to_moltype("dna")
        self.assertEqual(dna.moltype.label, "dna")
        rc = dna.rc().to_dict()
        expect = {"seq1": "TACGTACGT", "seq2": "---TTCGGT", "seq3": "AACGTACGT"}
        self.assertEqual(rc, expect)

        data = {"seq1": "TTTTTTAAAA", "seq2": "AAAATTTTTT", "seq3": "AATTTTTAAA"}
        seqs = self.Class(data=data)
        rna = seqs.to_moltype("rna")
        self.assertEqual(rna.moltype.label, "rna")
        rc = rna.rc().to_dict()
        expect = {"seq1": "UUUUAAAAAA", "seq2": "AAAAAAUUUU", "seq3": "UUUAAAAAUU"}
        self.assertEqual(rc, expect)
        # calling with a null object should raise an exception
        with self.assertRaises(ValueError):
            seqs.to_moltype(None)

        with self.assertRaises(ValueError):
            seqs.to_moltype("")

    def test_to_moltype_info(self):
        """correctly convert to specified moltype"""
        data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        seqs = self.Class(data=data, info={"key": "value"})
        dna = seqs.to_moltype(DNA)
        self.assertEqual(dna.info["key"], "value")

    def test_get_lengths(self):
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

    def test_strand_symmetry(self):
        """exercising strand symmetry test"""
        data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        seqs = self.Class(data, moltype=DNA)
        result = seqs.strand_symmetry()
        assert_allclose(result["seq1"].observed.array, [[3, 2], [2, 2]])
        assert_allclose(result["seq2"].observed.array, [[3, 0], [2, 1]])

    def test_dotplot(self):
        """exercising dotplot method"""
        seqs = self.Class(data=self.brca1_data, moltype=DNA)
        _ = seqs.dotplot(show_progress=False)

    def test_dotplot_annotated(self):
        """exercising dotplot method with annotated sequences"""
        seqs = self.Class(data={"Human": "CAGATTTGGCAGTT-", "Mouse": "CAGATTCAGCAGGTG"})

        seqs = seqs.take_seqs(["Human", "Mouse"])

        if type(self.Class) != ArrayAlignment:
            # we annotated Human
            seq = seqs.get_seq("Human")
            _ = seq.add_feature("exon", "fred", [(10, 15)])

        _ = seqs.dotplot(show_progress=False)

    def test_rename_seqs(self):
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

    def test_apply_pssm(self):
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

    def test_set_repr_policy_no_input(self):
        """repr_policy should remain unchanged"""
        seqs = self.Class({"a": "AAAAA"})
        seqs.set_repr_policy(num_seqs=None, num_pos=None)
        self.assertEqual(
            seqs._repr_policy,
            dict(num_seqs=10, num_pos=60, ref_name="longest", wrap=60),
        )

    def test_set_repr_policy_invalid_input(self):
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

    def test_set_repr_policy_valid_input(self):
        """repr_policy should be set to new values"""
        seqs = self.Class({"a": "AAAAA", "b": "AAA--"})
        seqs.set_repr_policy(num_seqs=5, num_pos=40, ref_name="a", wrap=10)
        self.assertEqual(
            seqs._repr_policy, dict(num_seqs=5, num_pos=40, ref_name="a", wrap=10)
        )

        if self.Class == SequenceCollection:
            # this class cannot slice
            return True

        # should persist in slicing
        self.assertEqual(
            seqs[:2]._repr_policy, dict(num_seqs=5, num_pos=40, ref_name="a", wrap=10)
        )

    def test_set_wrap_affects_repr_html(self):
        """the wrap argument affects the number of columns"""
        if self.Class == SequenceCollection:
            # this class does not have this method
            return True

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
        self.assertEqual(got.count(token), 3 * orig.count(token))
        os.environ.pop(env_name, None)

    def test_get_seq_entropy(self):
        """get_seq_entropy should get entropy of each seq"""
        a = self.Class(dict(a="ACCC", b="AGTA"), moltype=DNA)
        entropy = a.entropy_per_seq()
        e = 0.81127812445913283  # sum(p log_2 p) for p = 0.25, 0.75
        assert_allclose(entropy, array([e, 1.5]))

    def test_write_to_json(self):
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

    def test_seq_len_get_ragged(self):
        """SequenceCollection seq_len get should work for ragged seqs"""
        self.assertEqual(self.ragged.seq_len, 6)

    def test_is_ragged_ragged(self):
        """SequenceCollection is_ragged should return True if ragged"""
        self.assertTrue(self.ragged.is_ragged())

    def test_Seqs_ragged(self):
        """SequenceCollection seqs should work on ragged alignment"""
        self.ragged.names = "bac"
        self.assertEqual(list(self.ragged.seqs), ["AAA", "AAAAAA", "AAAA"])

    def test_iter_seqs_ragged(self):
        """SequenceCollection iter_seqs() method should support reordering of seqs"""
        self.ragged.names = ["a", "b", "c"]
        seqs = list(self.ragged.iter_seqs())
        self.assertEqual(seqs, ["AAAAAA", "AAA", "AAAA"])
        seqs = list(self.ragged.iter_seqs(seq_order=["b", "a", "a"]))
        self.assertEqual(seqs, ["AAA", "AAAAAA", "AAAAAA"])
        self.assertIs(seqs[1], seqs[2])
        self.assertIs(seqs[0], self.ragged.named_seqs["b"])

    def test_toPHYLIP_ragged(self):
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

    def test_pad_seqs_ragged(self):
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

    def test_info_source(self):
        """info.source exists if load seqs given a filename"""
        seqs = load_unaligned_seqs("data/brca1.fasta")
        self.assertEqual(seqs.info.source, "data/brca1.fasta")

    def test_apply_pssm2(self):
        """apply_pssm fail if ragged sequences"""
        data = {
            "ENSMUSG00000056468": "GCCAGGGGGGAAAGGGAGAA",
            "ENSMUSG00000039616": "GCCCTTCAAATTT",
        }
        seqs = self.Class(data=data, moltype=DNA)
        with self.assertRaises(AssertionError):
            _ = seqs.apply_pssm(path="data/sample.jaspar", show_progress=False)

    def test_construction(self):
        """correctly construct from list of sequences of length 2"""
        seq1 = make_seq("AC", name="seq1")
        seq2 = make_seq("AC", name="seq2")
        coll = SequenceCollection(data=[seq1, seq2])


def _make_filter_func(aln):
    array_align = type(aln) == ArrayAlignment
    if array_align:
        gap = aln.alphabet.with_gap_motif().to_indices("-")[0]
    else:
        gap = "-"

    def func_str(x):
        return gap not in "".join(x)

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
        """Tests that the alignment_quality generates the right alignment quality
        value based on the Hertz-Stormo metric. expected values are hand calculated
        using the formula in the paper."""
        aln = self.Class(["AATTGA", "AGGTCC", "AGGATG", "AGGCGT"], moltype="dna")
        got = aln.alignment_quality(equifreq_mprobs=True)
        expect = log2(4) + (3 / 2) * log2(3) + (1 / 2) * log2(2) + (1 / 2) * log2(2)
        assert_allclose(got, expect)
        # should be the same with the default moltype too
        aln = self.Class(["AATTGA", "AGGTCC", "AGGATG", "AGGCGT"])
        got = aln.alignment_quality(equifreq_mprobs=True)
        assert_allclose(got, expect)

        aln = self.Class(["AAAC", "ACGC", "AGCC", "A-TC"], moltype="dna")
        got = aln.alignment_quality(equifreq_mprobs=False)
        expect = (
            2 * log2(1 / 0.4)
            + log2(1 / (4 * 0.4))
            + (1 / 2) * log2(1 / (8 / 15))
            + (1 / 4) * log2(1 / (4 / 15))
        )
        assert_allclose(got, expect)

        # 1. Alignment just gaps - alignment_quality returns None
        aln = self.Class(["----", "----"])
        got = aln.alignment_quality(equifreq_mprobs=True)
        self.assertIsNone(got)

        # 2 Just one sequence - alignment_quality returns None
        aln = self.Class(["AAAC"])
        got = aln.alignment_quality(equifreq_mprobs=True)
        self.assertIsNone(got)

        # 3.1 Two seqs, one all gaps. (equifreq_mprobs=True)
        aln = self.Class(["----", "ACAT"])
        got = aln.alignment_quality(equifreq_mprobs=True)
        assert_allclose(got, 1.1699250014423124)

        # 3.2 Two seqs, one all gaps. (equifreq_mprobs=False)
        aln = self.Class(["----", "AAAA"])
        got = aln.alignment_quality(equifreq_mprobs=False)
        assert_allclose(got, -2)

    def make_and_filter(self, raw, expected, motif_length, drop_remainder):
        # a simple filter func
        aln = self.Class(raw, info={"key": "value"})
        func = _make_filter_func(aln)
        result = aln.filtered(
            func,
            motif_length=motif_length,
            log_warnings=False,
            drop_remainder=drop_remainder,
        )
        self.assertEqual(result.to_dict(), expected)
        self.assertEqual(result.info["key"], "value")

    def test_filtered(self):
        """filtered should return new alignment with positions consistent with
        provided callback function"""
        # a simple filter option
        raw = {"a": "ACGACGACG", "b": "CCC---CCC", "c": "AAAA--AAA"}
        self.make_and_filter(
            raw, {"a": "ACGACG", "b": "CCCCCC", "c": "AAAAAA"}, 1, True
        )
        # check with motif_length = 2
        self.make_and_filter(raw, {"a": "ACAC", "b": "CCCC", "c": "AAAA"}, 2, True)
        # check with motif_length = 3
        self.make_and_filter(
            raw, {"a": "ACGACG", "b": "CCCCCC", "c": "AAAAAA"}, 3, True
        )

    def test_filter_drop_remainder(self):
        """filter allows dropping"""
        raw = {"a": "ACGACGACG", "b": "CCC---CCC", "c": "AAAA--AAA"}
        aln = self.Class(raw)
        func = _make_filter_func(aln)
        got = aln.filtered(func, motif_length=1, log_warnings=False)
        self.assertEqual(len(got), 6)
        # raises an assertion if the length is not modulo
        with self.assertRaises(ValueError):
            # because alignment not modulo 2
            got = aln.filtered(func, motif_length=2, drop_remainder=False)
        got = aln.filtered(
            func, motif_length=2, drop_remainder=True, log_warnings=False
        )
        self.assertEqual(len(got), 4)

    def test_positions(self):
        """SequenceCollection positions property should iterate over positions, using self.names"""
        r = self.Class({"a": "AAAAAA", "b": "AAA---", "c": "AAAA--"})
        r.names = ["a", "b", "c"]
        self.assertEqual(
            list(r.positions),
            list(map(list, ["AAA", "AAA", "AAA", "A-A", "A--", "A--"])),
        )

    def test_iter_positions(self):
        # """SequenceCollection iter_positions() method should support reordering of #cols"""
        r = self.Class(self.ragged_padded.named_seqs, names=["c", "b"])
        self.assertEqual(
            list(r.iter_positions(pos_order=[5, 1, 3])),
            list(map(list, ["--", "AA", "A-"])),
        )
        # reorder names
        r = self.Class(self.ragged_padded.named_seqs, names=["a", "b", "c"])
        cols = list(r.iter_positions())
        self.assertEqual(
            cols, list(map(list, ["AAA", "AAA", "AAA", "A-A", "A--", "A--"]))
        )

    def test_take_positions(self):
        """SequenceCollection take_positions should return new alignment w/ specified pos"""
        self.assertEqual(
            self.gaps.take_positions([5, 4, 0]), {"a": "AAA", "b": "A-A", "c": "--A"}
        )
        self.assertTrue(
            isinstance(self.gaps.take_positions([0]), _SequenceCollectionBase)
        )
        # should be able to negate
        self.assertEqual(
            self.gaps.take_positions([5, 4, 0], negate=True),
            {"a": "AAAA", "b": "--AA", "c": "A---"},
        )

    def test_take_positions_info(self):
        aln = self.Class(
            {"a": "AAAAAAA", "b": "A--A-AA", "c": "AA-----"}, info={"key": "value"}
        )
        tps = aln.take_positions([5, 4, 0])
        self.assertEqual(tps.info["key"], "value")

    def test_get_position_indices(self):
        """SequenceCollection get_position_indices should return names of cols where f(col)"""

        def gap_1st(x):
            return x[0] == "-"

        def gap_2nd(x):
            return x[1] == "-"

        def gap_3rd(x):
            return x[2] == "-"

        def is_list(x):
            return isinstance(x, list)

        self.gaps = self.Class(self.gaps.named_seqs, names=["a", "b", "c"])

        self.assertEqual(self.gaps.get_position_indices(gap_1st), [])
        self.assertEqual(self.gaps.get_position_indices(gap_2nd), [1, 2, 4])
        self.assertEqual(self.gaps.get_position_indices(gap_3rd), [2, 3, 4, 5, 6])
        self.assertEqual(self.gaps.get_position_indices(is_list), [0, 1, 2, 3, 4, 5, 6])
        # should be able to negate
        self.assertEqual(
            self.gaps.get_position_indices(gap_2nd, negate=True), [0, 3, 5, 6]
        )
        self.assertEqual(
            self.gaps.get_position_indices(gap_1st, negate=True), [0, 1, 2, 3, 4, 5, 6]
        )
        self.assertEqual(self.gaps.get_position_indices(is_list, negate=True), [])

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

    def test_no_degenerates(self):
        """no_degenerates correctly excludes columns containing IUPAC ambiguity codes"""
        data = {
            "s1": "AAA CCC GGG TTT".replace(" ", ""),
            "s2": "CCC GGG T-T AAA".replace(" ", ""),
            "s3": "GGR YTT AAA CCC".replace(" ", ""),
        }
        aln = self.Class(data=data, moltype=DNA)

        # motif length of 1, defaults - no gaps allowed
        result = aln.no_degenerates().to_dict()
        expect = {
            "s1": "AA CC GG TTT".replace(" ", ""),
            "s2": "CC GG TT AAA".replace(" ", ""),
            "s3": "GG TT AA CCC".replace(" ", ""),
        }
        self.assertEqual(result, expect)

        # allow gaps
        result = aln.no_degenerates(allow_gap=True).to_dict()
        expect = {
            "s1": "AA CC GGG TTT".replace(" ", ""),
            "s2": "CC GG T-T AAA".replace(" ", ""),
            "s3": "GG TT AAA CCC".replace(" ", ""),
        }
        self.assertEqual(result, expect)

        # motif length of 3, defaults - no gaps allowed
        result = aln.no_degenerates(motif_length=3).to_dict()
        expect = {
            "s1": "TTT".replace(" ", ""),
            "s2": "AAA".replace(" ", ""),
            "s3": "CCC".replace(" ", ""),
        }
        self.assertEqual(result, expect)

        # allow gaps
        result = aln.no_degenerates(motif_length=3, allow_gap=True).to_dict()
        expect = {
            "s1": "GGG TTT".replace(" ", ""),
            "s2": "T-T AAA".replace(" ", ""),
            "s3": "AAA CCC".replace(" ", ""),
        }
        self.assertEqual(result, expect)

        # raises ValueError if a default moltype -- with no
        # degen characters -- is used
        aln = self.Class(data=data)
        self.assertRaises(ValueError, aln.no_degenerates)

    def test_omit_gap_pos(self):
        """Alignment omit_gap_pos should return alignment w/o positions of gaps"""
        aln = self.end_gaps
        # first, check behavior when we're just acting on the cols (and not
        # trying to delete the naughty seqs).

        # default should strip out cols that are 100% gaps
        result = aln.omit_gap_pos()
        self.assertEqual(result.to_dict(), {"a": "-ABC", "b": "CBA-", "c": "-DEF"})
        # if allowed_gap_frac is 1, shouldn't delete anything
        self.assertEqual(
            aln.omit_gap_pos(1).to_dict(),
            {"a": "--A-BC-", "b": "-CB-A--", "c": "--D-EF-"},
        )
        # if allowed_gap_frac is 0, should strip out any cols containing gaps
        self.assertEqual(
            aln.omit_gap_pos(0).to_dict(), {"a": "AB", "b": "BA", "c": "DE"}
        )
        # intermediate numbers should work as expected
        self.assertEqual(
            aln.omit_gap_pos(0.4).to_dict(), {"a": "ABC", "b": "BA-", "c": "DEF"}
        )
        self.assertEqual(
            aln.omit_gap_pos(0.7).to_dict(), {"a": "-ABC", "b": "CBA-", "c": "-DEF"}
        )

        # when we increase the number of sequences to 6, more differences
        # start to appear.
        new_aln_data = aln.named_seqs.copy()
        new_aln_data["d"] = "-------"
        new_aln_data["e"] = "XYZXYZX"
        new_aln_data["f"] = "AB-CDEF"
        aln = self.Class(new_aln_data)

        # if no gaps are allowed, we get None
        self.assertEqual(aln.omit_gap_pos(0), None)

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

    def test_sample(self):
        """Alignment.sample should permute alignment by default"""
        alignment = self.Class({"seq1": "ABCDEFGHIJKLMNOP", "seq2": "ABCDEFGHIJKLMNOP"})
        # effectively permute columns, preserving length
        shuffled = alignment.sample()
        self.assertEqual(len(shuffled), len(alignment))
        # ensure length correct
        sample = alignment.sample(10)
        self.assertEqual(len(sample), 10)
        # test columns alignment preserved
        seqs = list(sample.to_dict().values())
        self.assertEqual(seqs[0], seqs[1])
        # ensure each char occurs once as sampling without replacement
        for char in seqs[0]:
            self.assertEqual(seqs[0].count(char), 1)

    def test_sample_info(self):
        """Alignment.sample should preserver info attribute"""
        alignment = self.Class(
            {"seq1": "ABCDEFGHIJKLMNOP", "seq2": "ABCDEFGHIJKLMNOP"},
            info={"key": "value"},
        )
        # effectively permute columns, preserving length
        shuffled = alignment.sample()
        self.assertEqual(shuffled.info["key"], "value")
        # ensure length correct
        sample = alignment.sample(10)
        self.assertEqual(sample.info["key"], "value")

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

    def test_to_type(self):
        """correctly interconvert between alignment types"""
        new_seqs = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        array_align = self.Class == ArrayAlignment
        # when array_align arg matches instance class, no conversion
        # and get back self
        aln = self.Class(data=new_seqs)
        new = aln.to_type(array_align=array_align)
        self.assertEqual(id(aln), id(new))

        # when array_align arg does not match, should get back the opposite type
        new = aln.to_type(array_align=not array_align)
        self.assertFalse(isinstance(new, self.Class))

        # we should be able to specify moltype and alignment
        new = aln.to_type(array_align=not array_align, moltype=DNA)
        self.assertEqual(new.to_dict(), new_seqs)
        # and translate
        self.assertEqual(
            new.get_translation().to_dict(),
            {"seq1": "TYV", "seq3": "TYV", "seq2": "TE-"},
        )

        # should work on ArrayAlign when just moltype changes
        new_seqs = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---"}
        aln = self.Class(data=new_seqs)
        new = aln.to_type(array_align=array_align, moltype=DNA)
        new = new.no_degenerates()  # this should not fail!
        self.assertEqual(len(new), len(aln) - 3)

        # should correctly apply to existing moltype
        aln = self.Class(data=new_seqs, moltype=DNA)
        new = aln.to_type(array_align=not array_align)
        self.assertEqual(aln.moltype, new.moltype)

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

    def test_replace_seqs(self):
        """replace_seqs should replace 1-letter w/ 3-letter seqs"""
        a = self.Class({"seq1": "ACGU", "seq2": "C-UA", "seq3": "C---"})
        seqs = {"seq1": "AAACCCGGGUUU", "seq2": "CCCUUUAAA", "seq3": "CCC"}
        result = a.replace_seqs(seqs)  # default behaviour
        self.assertEqual(
            result.to_fasta(),
            ">seq1\nAAACCCGGGUUU\n>seq2\nCCC---UUUAAA\n>seq3\nCCC---------\n",
        )

        result = a.replace_seqs(seqs, aa_to_codon=True)  # default behaviour
        self.assertEqual(
            result.to_fasta(),
            ">seq1\nAAACCCGGGUUU\n>seq2\nCCC---UUUAAA\n>seq3\nCCC---------\n",
        )

        # should correctly gap the same sequences with same length
        result = a.replace_seqs(a.degap(), aa_to_codon=False)  # default behaviour
        self.assertEqual(
            result.to_dict(), {"seq1": "ACGU", "seq2": "C-UA", "seq3": "C---"}
        )

        # should fail when not same length if aa_to_codon is False
        new = SequenceCollection(
            [(n, s.replace("-", "")) for n, s in list(a[:3].to_dict().items())]
        )
        self.assertRaises(ValueError, a.replace_seqs, new, aa_to_codon=False)

        # check the gaps are changed
        aln1 = self.Class(data={"a": "AC-CT", "b": "ACGCT"})
        aln2 = self.Class(data={"a": "ACC-T", "b": "ACGCT"})

        result = aln1.replace_seqs(aln2, aa_to_codon=False)
        self.assertTrue(id(aln1) != id(aln2))
        self.assertEqual(aln1.to_dict(), result.to_dict())

    def test_replace_seqs_info(self):
        """replace_seqs should preserve info attribute"""
        a = self.Class(
            {"seq1": "ACGU", "seq2": "C-UA", "seq3": "C---"}, info={"key": "value"}
        )
        seqs = {"seq1": "AAACCCGGGUUU", "seq2": "CCCUUUAAA", "seq3": "CCC"}
        result = a.replace_seqs(seqs)  # default behaviour
        self.assertEqual(result.info["key"], "value")

    def test_counts(self):
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

    def test_counts_per_seq(self):
        """SequenceCollection.counts_per_seq handles motif length, allow_gaps etc.."""
        data = {"a": "AAAA??????", "b": "CCCGGG--NN", "c": "CCGGTTCCAA"}
        coll = self.Class(data=data, moltype="dna")
        mtype = coll.moltype
        got = coll.counts_per_seq()
        self.assertEqual(got["a", "A"], 4)
        self.assertEqual(len(got.motifs), len(mtype.alphabet))
        got = coll.counts_per_seq(include_ambiguity=True, allow_gap=True)
        # N, -, ? are the additional states
        self.assertEqual(len(got.motifs), 7)
        expect = {"-": 2, "?": 0, "A": 0, "C": 3, "G": 3, "N": 2, "T": 0}
        b = got["b"].to_dict()
        for k in expect:
            self.assertEqual(b[k], expect[k])

        got = coll.counts_per_seq(motif_length=2)
        self.assertEqual(len(got.motifs), len(mtype.alphabet) ** 2)
        self.assertEqual(got["a", "AA"], 2)
        self.assertEqual(got["b", "GG"], 1)
        got = coll.counts_per_seq(exclude_unobserved=True)
        expect = {"C": 4, "G": 2, "T": 2, "A": 2}
        c = got["c"].to_dict()
        for k in expect:
            self.assertEqual(c[k], expect[k])

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

        s1 = DNA.make_seq("TCAGAG", name="s1")
        s2 = DNA.make_seq("CCACAC", name="s2")
        s3 = DNA.make_seq("AGATAT", name="s3")
        s4 = DNA.make_seq("G-ACCC", name="s4")
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

    def test_counts_per_pos_default_moltype(self):
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

    def test_get_seq_entropy(self):
        """ArrayAlignment get_seq_entropy should get entropy of each seq"""
        seqs = [AB.make_seq(s, preserve_case=True) for s in ["abab", "bbbb", "abbb"]]
        a = self.Class(seqs, alphabet=AB.alphabet)
        entropy = a.entropy_per_seq()
        e = 0.81127812445913283  # sum(p log_2 p) for p = 0.25, 0.75
        assert_allclose(entropy, array([1, 0, e]))

    def test_seq_entropy_just_gaps(self):
        """ArrayAlignment get_seq_entropy should get entropy of each seq"""
        a = self.Class(dict(a="A---", b="----"), moltype=DNA)
        entropy = a.entropy_per_seq()
        assert_allclose(entropy, [0, numpy.nan])
        a = self.Class(dict(a="----", b="----"), moltype=DNA)
        entropy = a.entropy_per_seq()
        self.assertIs(entropy, None)

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

    def test_quick_tree(self):
        """quick tree method returns tree"""
        aln = self.Class(self.brca1_data, moltype=DNA)
        aln = aln.take_seqs(aln.names[:5])[:100]
        tree = aln.quick_tree(calc="hamming", show_progress=False)
        # bootstrap
        tree = aln.quick_tree(calc="hamming", bootstrap=2, show_progress=False)
        self.assertEqual(set(tree.get_tip_names()), set(aln.names))
        for edge in tree.preorder():
            if edge.is_root():
                continue
            self.assertIsInstance(edge.params["support"], float)

    def test_get_gapped_seq(self):
        """Alignment.get_gapped_seq should return seq, with gaps"""
        aln = self.Class({"seq1": "--TTT?", "seq2": "GATC??"})
        self.assertEqual(str(aln.get_gapped_seq("seq1")), "--TTT?")

    def test_count_gaps_per_pos(self):
        """correctly compute the number of gaps"""
        data = {"a": "AAAA---GGT", "b": "CCC--GG?GT"}
        aln = self.Class(data=data, moltype=DNA)
        # per position
        got = aln.count_gaps_per_pos(include_ambiguity=False)
        assert_equal(got.array, [0, 0, 0, 1, 2, 1, 1, 0, 0, 0])
        got = aln.count_gaps_per_pos(include_ambiguity=True)
        assert_equal(got.array, [0, 0, 0, 1, 2, 1, 1, 1, 0, 0])

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

    def test_coevolution(self):
        """correctly produces matrix of coevo measures"""
        data = {"s0": "AA", "s1": "AA", "s2": "BB", "s3": "BB", "s4": "BC"}
        aln = self.Class(data=data)
        coevo = aln.coevolution(method="rmi", show_progress=False)
        expect = array([[nan, nan], [0.78333333, nan]])
        assert_allclose(coevo.array, expect)
        coevo = aln.coevolution(method="nmi", show_progress=False)
        self.assertNotEqual(coevo[1, 1], expect[1, 1])
        # now check invoking drawable produces a result object with a drawable
        # attribute
        coevo = aln.coevolution(method="nmi", drawable="box", show_progress=False)
        self.assertTrue(hasattr(coevo, "drawable"))
        aln = load_aligned_seqs("data/brca1.fasta", moltype="dna")
        aln = aln.take_seqs(aln.names[:20])
        aln = aln.no_degenerates()[:20]

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

    def test_repr_html(self):
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

    def test_sliced_deepcopy(self):
        """correctly deep copy aligned objects in an alignment"""

        def eval_data(data):
            orig = self.Class(data)
            aln = orig[2:5]

            notsliced = aln.deepcopy(sliced=False)
            sliced = aln.deepcopy(sliced=True)
            for name in orig.names:
                # if not sliced, underlying seq data should be same length
                # as original
                self.assertEqual(
                    len(notsliced.named_seqs[name].data),
                    len(orig.named_seqs[name].data),
                )
                # and the map.parent_length attributes should match
                self.assertEqual(
                    notsliced.named_seqs[name].map.parent_length,
                    orig.named_seqs[name].map.parent_length,
                )
                # and map.parent_length and len(data) should match
                self.assertEqual(
                    sliced.named_seqs[name].map.parent_length,
                    len(sliced.named_seqs[name].data),
                )

                # if sliced, seq data should be < orig
                self.assertLess(
                    len(sliced.named_seqs[name].data), len(orig.named_seqs[name].data)
                )
                # and map.parent_length and len(data) should match
                self.assertEqual(
                    sliced.named_seqs[name].map.parent_length,
                    len(sliced.named_seqs[name].data),
                )

        eval_data({"seq1": "ACAACGACG", "seq2": "ACGACGACG"})
        eval_data({"seq1": "-CAACGACG", "seq2": "ACGACGACG"})
        eval_data({"seq1": "ACAACGAC-", "seq2": "ACGACGACG"})
        eval_data({"seq1": "AC-ACGACG", "seq2": "ACGACGACG"})
        eval_data({"seq1": "ACA-CGAC-", "seq2": "ACGACGACG"})

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

    def test_get_degapped_relative_to(self):
        """should remove all columns with a gap in sequence with given name"""
        aln = self.Class(
            [
                ["name1", "-AC-DEFGHI---"],
                ["name2", "XXXXXX--XXXXX"],
                ["name3", "YYYY-YYYYYYYY"],
                ["name4", "-KL---MNPR---"],
            ]
        )
        out_aln = self.Class(
            [
                ["name1", "ACDEFGHI"],
                ["name2", "XXXX--XX"],
                ["name3", "YY-YYYYY"],
                ["name4", "KL--MNPR"],
            ]
        )
        self.assertEqual(aln.get_degapped_relative_to("name1"), out_aln)

        self.assertRaises(ValueError, aln.get_degapped_relative_to, "nameX")

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

    def test_add_from_ref_aln(self):
        """should add or insert seqs based on align to reference"""
        aln1 = self.Class(
            [
                ["name1", "-AC-DEFGHI---"],
                ["name2", "XXXXXX--XXXXX"],
                ["name3", "YYYY-YYYYYYYY"],
            ]
        )

        aln2 = self.Class(
            [
                ["name1", "ACDEFGHI"],
                ["name4", "KL--MNPR"],
                ["name5", "KLACMNPR"],
                ["name6", "KL--MNPR"],
            ]
        )

        aligned_to_ref_out_aln_inserted = self.Class(
            [
                ["name1", "-AC-DEFGHI---"],
                ["name4", "-KL---MNPR---"],
                ["name5", "-KL-ACMNPR---"],
                ["name6", "-KL---MNPR---"],
                ["name2", "XXXXXX--XXXXX"],
                ["name3", "YYYY-YYYYYYYY"],
            ]
        )

        aln2_wrong_refseq = self.Class((("name1", "ACDXFGHI"), ("name4", "KL--MNPR")))

        aln2_wrong_refseq_name = self.Class(
            [["nameY", "ACDEFGHI"], ["name4", "KL--MNPR"]]
        )

        aln2_different_aln_class = ArrayAlignment(
            [["name1", "ACDEFGHI"], ["name4", "KL--MNPR"]]
        )

        aln2_list = [["name1", "ACDEFGHI"], ["name4", "KL--MNPR"]]

        aligned_to_ref_out_aln = self.Class(
            [
                ["name1", "-AC-DEFGHI---"],
                ["name2", "XXXXXX--XXXXX"],
                ["name3", "YYYY-YYYYYYYY"],
                ["name4", "-KL---MNPR---"],
            ]
        )

        out_aln = aln1.add_from_ref_aln(aln2, after_name="name1")
        self.assertEqual(
            str(aligned_to_ref_out_aln_inserted), str(out_aln)
        )  # test insert_after

        out_aln = aln1.add_from_ref_aln(aln2, before_name="name2")
        self.assertEqual(aligned_to_ref_out_aln_inserted, out_aln)  # test insert_before

        self.assertRaises(
            ValueError, aln1.add_from_ref_aln, aln2_wrong_refseq_name
        )  # test wrong_refseq_name

        aln = aln1.add_from_ref_aln(aln2_different_aln_class)
        self.assertEqual(
            aligned_to_ref_out_aln, aln
        )  # test_align_to_refseq_different_aln_class

        aln = aln1.add_from_ref_aln(aln2_list)
        self.assertEqual(aligned_to_ref_out_aln, aln)  # test from_list

        self.assertRaises(
            ValueError, aln1.add_from_ref_aln, aln2_wrong_refseq
        )  # test wrong_refseq

    def test_annotate_matches_to(self):
        """Aligned.annotate_matches_to correctly delegates to sequence"""
        from cogent3 import get_code

        aln = Alignment(dict(x="TTCCACTTCCGCTT"), moltype="dna")
        seq = aln.named_seqs["x"]
        pattern = "CCRC"
        annot = seq.annotate_matches_to(
            pattern=pattern, annot_type="domain", name="fred", allow_multiple=True
        )
        got = [a.get_slice() for a in annot]
        matches = ["CCAC", "CCGC"]
        self.assertEqual(got, matches)
        annot = seq.annotate_matches_to(
            pattern=pattern, annot_type="domain", name="fred", allow_multiple=False
        )
        got = [a.get_slice() for a in annot]
        self.assertEqual(got, matches[:1])

        # handles regex from aa
        aln = Alignment(dict(x="TTCCACTTCCGCTT"), moltype="dna")
        gc = get_code(1)
        aa_regex = gc.to_regex("FHF")
        s = aln.named_seqs["x"].annotate_matches_to(
            aa_regex, "domain", "test", allow_multiple=False
        )
        a = list(aln.get_annotations_from_seq("x"))[0]
        self.assertEqual(str(aln[a].named_seqs["x"]), "TTCCACTTC")

    def test_deepcopy(self):
        """correctly deepcopy Aligned objects in an alignment"""
        path = "data/brca1_5.paml"
        # generates an annotatable Alignment object
        aln = load_aligned_seqs(path, array_align=False, moltype="dna")
        # when the annotation is outside(before) boundary of the slice
        aln.named_seqs["NineBande"].data.add_annotation(
            Feature, "exon", "annot1", [(0, 10)]
        )
        # when the annotation is across boundary of the slice
        aln.named_seqs["Mouse"].data.add_annotation(
            Feature, "exon", "annot2", [(10, 21)]
        )
        # when the annotation is within boundary of the slice
        aln.named_seqs["Human"].data.add_annotation(
            Feature, "exon", "annot3", [(20, 25)]
        )
        # when the annotation is across boundary of the slice
        aln.named_seqs["HowlerMon"].data.add_annotation(
            Feature, "exon", "annot4", [(25, 32)]
        )
        # when the annotation is outside(after) boundary of the slice
        aln.named_seqs["DogFaced"].data.add_annotation(
            Feature, "exon", "annot5", [(40, 45)]
        )
        aln = aln[20:30]

        # for these species, each has an annotation spanning slice boundary or within it
        for name in ["Mouse", "Human", "HowlerMon"]:
            new_seq = aln.named_seqs[name].deepcopy(sliced=True)
            seq = aln.named_seqs[name]
            self.assertNotEqual(new_seq.map.parent_length, seq.map.parent_length)

            self.assertEqual(len(new_seq.data), 10)
            self.assertTrue(new_seq.data.is_annotated())
            self.assertEqual(len(new_seq.data.annotations), 1)
            # tests the case when sliced argument if False
            new_seq = aln.named_seqs[name].deepcopy(sliced=False)
            self.assertEqual(new_seq.map.parent_length, seq.map.parent_length)
            self.assertEqual(len(new_seq.data), len(aln.named_seqs[name].data))
            self.assertTrue(new_seq.data.is_annotated())
        # for these species, each has an annotation outside slice
        for name in ["NineBande", "DogFaced"]:
            new_seq = aln.named_seqs[name].deepcopy(sliced=True)
            self.assertEqual(len(new_seq.data), 10)
            self.assertFalse(new_seq.data.is_annotated())
            # tests the case when sliced argument if False
            new_seq = aln.named_seqs[name].deepcopy(sliced=False)
            self.assertEqual(len(new_seq.data), len(aln.named_seqs[name].data))
            self.assertTrue(new_seq.data.is_annotated())
            self.assertEqual(len(new_seq.data.annotations), 1)

        # add another human annotation that is outside slice
        aln.named_seqs["Human"].data.add_annotation(
            Feature, "exon", "annot6", [(40, 45)]
        )
        # tests the case when sliced argument if False regarding the Human sequence
        new_seq = aln.named_seqs["Human"].deepcopy(sliced=False)
        self.assertEqual(len(new_seq.data), len(aln.named_seqs["Human"].data))
        self.assertTrue(new_seq.data.is_annotated())
        self.assertEqual(len(new_seq.data.annotations), 2)

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

    def test_dotplot(self):
        """exercising dotplot method"""
        aln = self.Class([["name1", "TTTTTTAAAA"], ["name2", "AAAATTTTTT"]])
        aln = aln[2:8]
        draw = aln.dotplot(show_progress=False)
        expected = set([("name1", "TTTTAA"), ("name2", "AATTTT")])
        got = {(s.name, s._seq) for s in (draw.seq1, draw.seq2)}
        self.assertEqual(got, expected)

    def test_to_moltype_annotations(self):
        """correctly convert to specified moltype with proper sequence annotations"""
        s1 = Sequence("TTTTTTAAAA", name="test_seq1")
        s2 = Sequence("AAAATTTTTT", name="test_seq2")
        s3 = Sequence("AATTTTTAAA", name="test_seq3")
        s1.add_annotation(Feature, "exon", "fred", [(0, 6)])
        s2.add_annotation(Feature, "exon", "fred", [(4, 10)])
        s3.add_annotation(Feature, "exon", "fred", [(2, 7)])
        data = {"seq1": s1, "seq2": s2, "seq3": s3}
        aln = self.Class(data=data)
        aln.add_feature("demo", "one", [(0, 1), (2, 4)])
        rna = aln.to_moltype("rna")
        for name in rna.names:
            orig_seq = aln.get_seq(name)
            new_seq = rna.get_seq(name)
            self.assertEqual(len(orig_seq.annotations), len(new_seq.annotations))
            for src, dest in zip(orig_seq.annotations, new_seq.annotations):
                self.assertEqual(src.get_coordinates(), dest.get_coordinates())
                self.assertIsInstance(src, dest.__class__)
                self.assertIs(dest.parent, new_seq)
        # check the sequence moltypes
        self.assertEqual({s.data.moltype.label for s in rna.seqs}, {"rna"})
        self.assertEqual(rna.moltype.label, "rna")

    def test_get_annotations_from_any_seq(self):
        """get_annotations_from_any_seq returns correct annotations"""
        data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        seqs = self.Class(data, moltype=DNA)
        seqs.get_seq("seq1").add_annotation(Feature, "exon", "annotation1", [(3, 8)])
        seqs.get_seq("seq2").add_annotation(Feature, "exon", "annotation2", [(1, 2)])
        seqs.get_seq("seq3").add_annotation(Feature, "exon", "annotation3", [(3, 6)])
        got = seqs.get_annotations_from_any_seq()
        self.assertEqual(len(got), 3)
        self.assertEqual(str(got[0]), 'exon "annotation1" at [3:8]/9')
        self.assertEqual(str(got[1]), 'exon "annotation2" at [1:2]/9')
        self.assertEqual(str(got[2]), 'exon "annotation3" at [3:6]/9')

        got = seqs.get_annotations_from_any_seq(annotation_type="*", name="annotation1")
        self.assertEqual(len(got), 1)
        self.assertEqual(str(got[0]), 'exon "annotation1" at [3:8]/9')

        got = seqs.get_annotations_from_any_seq(
            annotation_type="exon", name="annotation2"
        )
        self.assertEqual(len(got), 1)
        self.assertEqual(str(got[0]), 'exon "annotation2" at [1:2]/9')

        got = seqs.get_annotations_from_any_seq(annotation_type="*", name="annotation3")
        self.assertEqual(len(got), 1)
        self.assertEqual(str(got[0]), 'exon "annotation3" at [3:6]/9')

    def test_rename_handles_annotations(self):
        """rename seqs on Alignment preserves annotations"""
        from cogent3.core.annotation import Feature

        data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
        seqs = self.Class(data, moltype=DNA)
        x = seqs.get_seq("seq1").add_annotation(Feature, "exon", "fred", [(3, 8)])
        expect = str(x.get_slice())
        new = seqs.rename_seqs(lambda x: x.upper())
        got = list(new.get_annotations_from_any_seq("exon"))[0]
        self.assertEqual(got.get_slice().to_dict()["SEQ1"], expect)

    def test_construction(self):
        """correctly construct from list of sequences of length 2"""
        seq1 = make_seq("AC-", name="seq1")
        seq1 = Aligned(*seq1.parse_out_gaps())
        seq2 = make_seq("ACG", name="seq2")
        seq2 = Aligned(*seq2.parse_out_gaps())
        coll = self.Class(data=[seq1, seq2])


class ArrayAlignmentSpecificTests(TestCase):
    """Tests of the ArrayAlignment object and its methods"""

    def setUp(self):
        """Define some standard alignments."""
        self.a = ArrayAlignment(
            array([[0, 1, 2], [3, 4, 5]]), conversion_f=aln_from_array
        )
        self.a2 = ArrayAlignment(["ABC", "DEF"], names=["x", "y"])
        seqs = []
        for s in ["abaa", "abbb"]:
            seqs.append(AB.make_seq(s, preserve_case=True))
        self.a = ArrayAlignment(seqs, alphabet=AB.alphabet)
        self.b = Alignment(["ABC", "DEF"])
        self.c = SequenceCollection(["ABC", "DEF"])

    def test_init(self):
        """ArrayAlignment init should work from a sequence"""
        a = ArrayAlignment(array([[0, 1, 2], [3, 4, 5]]), conversion_f=aln_from_array)
        assert_equal(a.seq_data, array([[0, 3], [1, 4], [2, 5]], "B"))
        assert_equal(a.array_positions, array([[0, 1, 2], [3, 4, 5]], "B"))
        assert_equal(a.names, ["seq_0", "seq_1", "seq_2"])

    def test_guess_input_type(self):
        """ArrayAlignment _guess_input_type should figure out data type correctly"""
        git = self.a._guess_input_type
        self.assertEqual(git(self.a), "array_aln")
        self.assertEqual(git(self.b), "aln")
        self.assertEqual(git(self.c), "collection")
        self.assertEqual(git(">ab\nabc"), "fasta")
        self.assertEqual(git([">ab", "abc"]), "fasta")
        self.assertEqual(git(["abc", "def"]), "generic")
        # precedence over generic
        self.assertEqual(git([[1, 2], [4, 5]]), "kv_pairs")
        self.assertEqual(git([[1, 2, 3], [4, 5, 6]]), "generic")
        self.assertEqual(git([ArraySequence("abc")]), "array_seqs")
        self.assertEqual(git(array([[1, 2, 3], [4, 5, 6]])), "array")
        self.assertEqual(git({"a": "aca"}), "dict")
        self.assertEqual(git([]), "empty")

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

    def test_init_empty(self):
        """ArrayAlignment init should fail if empty."""
        self.assertRaises(TypeError, ArrayAlignment)
        self.assertRaises(ValueError, ArrayAlignment, 3)

    def test_get_alphabet_and_moltype(self):
        """ArrayAlignment should figure out correct alphabet and moltype"""
        s1 = "A"
        s2 = RNA.make_seq("AA")

        d = ArrayAlignment(s1)
        self.assertIs(d.moltype, BYTES)
        self.assertIs(d.alphabet, BYTES.alphabet)

        d = ArrayAlignment(s1, moltype=RNA)
        self.assertIs(d.moltype, RNA)
        self.assertIs(d.alphabet, RNA.alphabets.degen_gapped)

        d = ArrayAlignment(s1, alphabet=RNA.alphabet)
        self.assertIs(d.moltype, RNA)
        self.assertIs(d.alphabet, RNA.alphabet)

        d = ArrayAlignment(s2)
        self.assertIs(d.moltype, RNA)
        self.assertIs(d.alphabet, RNA.alphabets.degen_gapped)

        d = ArrayAlignment(s2, moltype=DNA)
        self.assertIs(d.moltype, DNA)
        self.assertIs(d.alphabet, DNA.alphabets.degen_gapped)
        # checks for containers
        d = ArrayAlignment([s2])
        self.assertIs(d.moltype, RNA)
        d = ArrayAlignment({"x": s2})
        self.assertIs(d.moltype, RNA)
        d = ArrayAlignment(set([s2]))
        self.assertIs(d.moltype, RNA)

    def test_iter(self):
        """ArrayAlignment iter should iterate over positions"""
        result = list(iter(self.a2))
        for i, j in zip(result, [list(i) for i in ["AD", "BE", "CF"]]):
            self.assertEqual(i, j)

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
        a = ArrayAlignment(">x ABCE >y FGHI >z JKLM".split())
        # passing in positions should keep all seqs, but just selected
        # positions
        b = ArrayAlignment(">x BC >y GH >z KL".split())
        a_1 = a.get_sub_alignment(pos=[1, 2])
        self.assertEqual(a_1.names, b.names)

        self.assertEqual(a_1.seqs, b.seqs)
        # ...and with invert_pos, should keep all except the positions passed in
        a_2 = a.get_sub_alignment(pos=[0, 3], invert_pos=True)
        self.assertEqual(a_2.seqs, b.seqs)
        self.assertEqual(a_2.names, b.names)
        # passing in seqs should keep all positions, but just selected seqs
        c = ArrayAlignment(">x ABCE >z JKLM".split())
        a_3 = a.get_sub_alignment(seqs=[0, 2])
        self.assertEqual(a_3.seqs, c.seqs)
        # check that labels were updates as well...
        self.assertEqual(a_3.names, c.names)
        # ...and should work with invert_seqs to exclude just selected seqs
        a_4 = a.get_sub_alignment(seqs=[1], invert_seqs=True)
        self.assertEqual(a_4.seqs, c.seqs)
        self.assertEqual(a_4.names, c.names)
        # should be able to do both seqs and positions simultaneously
        d = ArrayAlignment(">x BC >z KL".split())
        a_5 = a.get_sub_alignment(seqs=[0, 2], pos=[1, 2])
        self.assertEqual(a_5.seqs, d.seqs)
        self.assertEqual(a_5.names, d.names)

    def test_get_sub_alignment_info(self):
        """ArrayAlignment get_sub_alignment should preserve info attribute"""
        a = ArrayAlignment(">x ABCE >y FGHI >z JKLM".split(), info={"key": "foo"})
        # passing in positions should keep all seqs, but just selected
        # positions
        b = ArrayAlignment(">x BC >y GH >z KL".split(), info={"key": "bar"})
        a_1 = a.get_sub_alignment(pos=[1, 2])
        self.assertEqual(a_1.names, b.names)
        self.assertEqual(a.info["key"], "foo")

    def test_str(self):
        """ArrayAlignment str should return FASTA representation of aln"""
        self.assertEqual(str(self.a2), ">x\nABC\n>y\nDEF\n")
        # should work if labels diff length
        self.a2.names[-1] = "yyy"
        self.assertEqual(str(self.a2), ">x\nABC\n>yyy\nDEF\n")

    def test_counts_per_seq(self):
        """ArrayAlignment counts_per_seq should return motif counts each seq"""
        a = self.a
        f = a.counts_per_seq()
        assert_equal(f.array, array([[3, 1], [1, 3]]))
        f = a.counts_per_seq(motif_length=2, exclude_unobserved=True)
        assert_equal(f.array, array([[1, 1, 0], [0, 1, 1]]))

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
            seqs.append(make_seq(s, moltype="dna"))
        a = ArrayAlignment(seqs)
        f = a.entropy_per_pos(allow_gap=True)
        e = array([1.584962500721156, 1.584962500721156, 1.584962500721156, 0])
        assert_allclose(f, e)

        seqs = []
        for s in ["-RAT", "ACCT", "GTGT"]:
            seqs.append(make_seq(s, moltype="dna"))
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

    def test_coevolution_segments(self):
        """specifying coordinate segments produces matrix with just those"""
        aln = load_aligned_seqs("data/brca1.fasta", moltype="dna")
        aln = aln.take_seqs(aln.names[:20])
        aln = aln.no_degenerates()[:20]
        coevo = aln.coevolution(segments=[(4, 6), (11, 13)], show_progress=False)
        self.assertEqual(coevo.template.names[0], [4, 5, 11, 12])


class IntegrationTests(TestCase):
    """Test for integration between regular and model seqs and alns"""

    def setUp(self):
        """Intialize some standard sequences"""
        self.r1 = RNA.make_seq("AAA", name="x")
        self.r2 = RNA.make_seq("CCC", name="y")
        self.m1 = RNA.make_array_seq("AAA", name="xx")
        self.m2 = RNA.make_array_seq("CCC", name="yy")

    def test_model_to_model(self):
        """Model seq should work with dense alignment"""
        a = ArrayAlignment([self.m1, self.m2])
        self.assertEqual(str(a), ">xx\nAAA\n>yy\nCCC\n")
        a = ArrayAlignment([self.m1, self.m2], moltype=DNA)
        self.assertEqual(str(a), ">xx\nAAA\n>yy\nCCC\n")
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


# run tests if invoked from command line
if __name__ == "__main__":
    main()
