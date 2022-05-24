#!/usr/bin/env python


from unittest import TestCase, main

from numpy import alltrue, array, transpose

from cogent3.core.alignment import Alignment, ArrayAlignment
from cogent3.core.moltype import RNA
from cogent3.core.sequence import ArraySequence, RnaSequence


__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Sandra Smit", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

from numpy.testing import assert_equal


class AllTests(TestCase):
    def setUp(self):
        """setUp method for all tests"""
        # named sequences
        self.rna1 = RnaSequence("UCAGGG", name="rna1")
        self.rna2 = RnaSequence("YCU-RG", name="rna2")
        self.rna3 = RnaSequence("CAA-NR", name="rna3")
        self.model1 = ArraySequence(
            "UCAGGG", name="rna1", alphabet=RNA.alphabets.degen_gapped
        )
        self.model2 = ArraySequence(
            "YCU-RG", name="rna2", alphabet=RNA.alphabets.degen_gapped
        )
        self.model3 = ArraySequence(
            "CAA-NR", name="rna3", alphabet=RNA.alphabets.degen_gapped
        )

        self.aln = Alignment([self.rna1, self.rna2, self.rna3], moltype=RNA)
        self.da = ArrayAlignment(
            [self.model1, self.model2, self.model3],
            moltype=RNA,
            alphabet=RNA.alphabets.degen_gapped,
        )

        # seqs no name
        self.nn_rna1 = RnaSequence("UCAGGG")
        self.nn_rna2 = RnaSequence("YCU-RG")
        self.nn_rna3 = RnaSequence("CAA-NR")

        self.nn_model1 = ArraySequence("UCAGGG", alphabet=RNA.alphabets.degen_gapped)
        self.nn_model2 = ArraySequence("YCU-RG", alphabet=RNA.alphabets.degen_gapped)
        self.nn_model3 = ArraySequence("CAA-NR", alphabet=RNA.alphabets.degen_gapped)

        self.nn_aln = Alignment([self.nn_rna1, self.nn_rna2, self.nn_rna3], moltype=RNA)
        self.nn_da = ArrayAlignment(
            [self.nn_model1, self.nn_model2, self.nn_model3],
            moltype=RNA,
            alphabet=RNA.alphabets.degen_gapped,
        )

    def test_printing_named_seqs(self):
        """Printing named seqs should work the same on Aln and DenseAln"""
        # Note: the newline trailing each sequence is intentional, because
        # we want each FASTA-format record to be separated.
        exp_lines_general = [">rna1", "UCAGGG", ">rna2", "YCU-RG", ">rna3", "CAA-NR"]
        self.assertEqual(str(self.aln), "\n".join(exp_lines_general) + "\n")
        self.assertEqual(str(self.da), "\n".join(exp_lines_general) + "\n")

    def test_printing_unnamed_seqs(self):
        """Printing unnamed sequences should work the same on Aln and DenseAln"""
        exp_lines_gen = [">seq_0", "UCAGGG", ">seq_1", "YCU-RG", ">seq_2", "CAA-NR\n"]
        self.assertEqual(str(self.nn_aln), "\n".join(exp_lines_gen))
        self.assertEqual(str(self.nn_da), "\n".join(exp_lines_gen))

    def test_ArrayAlignment_without_moltype(self):
        """Expect MolType to be picked up from the sequences."""

        m1 = ArraySequence("UCAG", alphabet=RNA.alphabets.degen_gapped, name="rna1")
        m2 = ArraySequence("CCCR", alphabet=RNA.alphabets.degen_gapped, name="rna2")
        da = ArrayAlignment([m1, m2])
        exp_lines = [">rna1", "UCAG", ">rna2", "CCCR"]
        self.assertEqual(str(da), "\n".join(exp_lines) + "\n")

    def test_names(self):
        # Should both alignments handle names the same way?
        self.assertEqual(self.aln.names, ["rna1", "rna2", "rna3"])
        self.assertEqual(self.da.names, ["rna1", "rna2", "rna3"])
        # On unnamed sequences the behavior is now the same.
        self.assertEqual(self.nn_aln.names, ["seq_0", "seq_1", "seq_2"])
        self.assertEqual(self.nn_da.names, ["seq_0", "seq_1", "seq_2"])

    def test_seqFreqs(self):
        """seqFreqs should work the same on Alignment and ArrayAlignment"""
        RNA.alphabets.degen_gapped.index
        # 'UCAGGG'
        # 'YCU-RG'
        # 'CAA-NR'

        expected_counts = {
            0: {"U": 1, "C": 1, "A": 1, "G": 3},
            1: {"Y": 1, "C": 1, "U": 1, "-": 1, "R": 1, "G": 1},
            2: {"C": 1, "A": 2, "-": 1, "N": 1, "R": 1},
        }
        got1 = self.da.counts_per_seq(allow_gap=True, include_ambiguity=True)
        got2 = self.aln.counts_per_seq(allow_gap=True, include_ambiguity=True)
        for pos, counts in expected_counts.items():
            for char in counts:
                self.assertEqual(got1[pos, char], expected_counts[pos][char])
                self.assertEqual(got2[pos, char], expected_counts[pos][char])

    def test_subset_positions_ArrayAlignment(self):
        # because dict order volatile, need to grab the
        # the index for ambig characters from the object
        # The full data comes from these seqs
        # 'UCAGGG'
        # 'YCU-RG'
        # 'CAA-NR'
        get_index = RNA.alphabets.degen_gapped.index
        G = get_index("-")
        N = get_index("N")
        R = get_index("R")
        Y = get_index("Y")
        full_data = array([[0, 1, 2, 3, 3, 3], [Y, 1, 0, G, R, 3], [1, 2, 2, G, N, R]])

        model1 = ArraySequence("UCG", name="rna1", alphabet=RNA.alphabets.degen_gapped)
        model2 = ArraySequence("YCG", name="rna2", alphabet=RNA.alphabets.degen_gapped)
        model3 = ArraySequence("CAR", name="rna3", alphabet=RNA.alphabets.degen_gapped)
        sub_da = ArrayAlignment(
            [model1, model2, model3], moltype=RNA, alphabet=RNA.alphabets.degen_gapped
        )

        sub_data = array([[0, 1, 3], [Y, 1, 3], [1, 2, R]])

        # First check some data
        assert_equal(self.da.array_seqs, full_data)
        assert_equal(self.da.array_positions, transpose(full_data))
        assert_equal(sub_da.array_seqs, sub_data)
        assert_equal(sub_da.array_positions, transpose(sub_data))

        obs_sub_da_TP = self.da.take_positions([0, 1, 5])
        obs_sub_da_SA = self.da.get_sub_alignment(pos=[0, 1, 5])

        # When using the get_sub_alignment method the data is right
        self.assertEqual(obs_sub_da_SA, sub_da)
        self.assertNotEqual(obs_sub_da_SA, self.da)
        assert_equal(obs_sub_da_SA.array_seqs, sub_data)
        assert_equal(obs_sub_da_SA.array_positions, transpose(sub_data))

        # For the take_positions method: Why does this work
        self.assertEqual(obs_sub_da_TP, sub_da)
        self.assertNotEqual(obs_sub_da_TP, self.da)
        # If the data doesn't match?
        assert_equal(obs_sub_da_TP.array_seqs, sub_data)
        assert_equal(obs_sub_da_TP.array_positions, transpose(sub_data))
        # Shouldn't the __eq__ method check the data at least?

    def test_subset_positions_Alignment(self):
        rna1 = RnaSequence("UCG", name="rna1")
        rna2 = RnaSequence("YCG", name="rna2")
        rna3 = RnaSequence("CAR", name="rna3")

        sub_aln = Alignment([rna1, rna2, rna3], moltype=RNA)

        obs_sub_aln = self.aln.take_positions([0, 1, 5])
        self.assertEqual(obs_sub_aln, sub_aln)
        self.assertNotEqual(obs_sub_aln, self.aln)
        # string representations should be the same. This fails right
        # now, because sequence order is not maintained. See separate test.
        self.assertEqual(str(obs_sub_aln), str(sub_aln))

    def test_take_positions_sequence_order(self):
        """Alignment take_positions should maintain seq order"""
        # This works
        self.assertEqual(self.da.names, ["rna1", "rna2", "rna3"])
        sub_da = self.da.get_sub_alignment(pos=[0, 1, 5])
        self.assertEqual(sub_da.names, ["rna1", "rna2", "rna3"])
        # seq order not maintained in Alignment
        self.assertEqual(self.aln.names, ["rna1", "rna2", "rna3"])
        sub_aln = self.aln.take_positions([0, 1, 5])
        self.assertEqual(sub_aln.names, ["rna1", "rna2", "rna3"])

    def test_subset_seqs_Alignment(self):
        rna1 = RnaSequence("UCG", name="rna1")
        rna2 = RnaSequence("YCG", name="rna2")
        rna3 = RnaSequence("CAR", name="rna3")

        sub_aln = Alignment([rna2, rna3], moltype=RNA)
        aln = Alignment([rna1, rna2, rna3], moltype=RNA)
        obs_sub_aln = aln.take_seqs(["rna2", "rna3"])

        self.assertEqual(obs_sub_aln, sub_aln)
        self.assertEqual(str(obs_sub_aln), str(sub_aln))

        # Selected sequences should be in specified order?
        obs_sub_aln_1 = self.aln.take_seqs(["rna3", "rna2"])
        obs_sub_aln_2 = self.aln.take_seqs(["rna2", "rna3"])
        self.assertNotEqual(str(obs_sub_aln_1), str(obs_sub_aln_2))

    def test_subset_seqs_ArrayAlignment(self):
        model1 = ArraySequence("UCG", name="rna1", alphabet=RNA.alphabets.degen_gapped)
        model2 = ArraySequence("YCG", name="rna2", alphabet=RNA.alphabets.degen_gapped)
        model3 = ArraySequence("CAR", name="rna3", alphabet=RNA.alphabets.degen_gapped)
        sub_da = ArrayAlignment(
            [model1, model2, model3], moltype=RNA, alphabet=RNA.alphabets.degen_gapped
        )

        # take_seqs by name should have the same effect as
        # get_sub_alignment by seq idx?
        obs_sub_da_TS = self.da.take_seqs(["rna1"])
        obs_sub_da_SA = self.da.get_sub_alignment(seqs=[0])

        # These two are now the same. Fixed mapping of key to char array.
        self.assertEqual(obs_sub_da_TS, obs_sub_da_SA)
        self.assertEqual(str(obs_sub_da_TS), str(obs_sub_da_SA))

    def test_aln_equality(self):
        # When does something compare equal?
        self.assertEqual(self.da == self.da, True)
        # one sequence less
        other_da1 = ArrayAlignment(
            [self.model1, self.model2], moltype=RNA, alphabet=RNA.alphabets.degen_gapped
        )
        self.assertEqual(self.da == other_da1, False)
        # seqs in different order -- doesn't matter
        other_da2 = ArrayAlignment(
            [self.model1, self.model3, self.model2],
            moltype=RNA,
            alphabet=RNA.alphabets.degen_gapped,
        )
        self.assertEqual(self.da == other_da2, True)
        # seqs in different encoding -- doesn't matter, only looks at data
        other_da3 = ArrayAlignment([self.model1, self.model2, self.model3])
        # Should this compare False even though the data is exactly the same?
        # The moltype is different...
        self.assertEqual(self.da == other_da3, True)
        assert alltrue(list(map(alltrue, self.da.array_seqs == other_da3.array_seqs)))

    def test_seq_equality(self):
        model1 = ArraySequence("UCG", name="rna1", alphabet=RNA.alphabets.degen_gapped)
        model2 = ArraySequence("UCG", name="rna1", alphabet=RNA.alphabets.degen_gapped)
        # Shouldn't the above two sequences be equal?
        self.assertEqual(model1, model2)
        # string comparison is True
        self.assertEqual(str(model1), str(model2))

    def test_seq_ungapping(self):
        rna1 = RnaSequence("U-C-A-G-", name="rna1")
        model1 = ArraySequence(
            "U-C-A-G-", name="rna1", alphabet=RNA.alphabets.degen_gapped
        )

        self.assertEqual(rna1, "U-C-A-G-")
        self.assertEqual(rna1.degap(), "UCAG")

        # check is produces the right string from the beginning
        self.assertEqual(str(model1), "U-C-A-G-")
        assert_equal(model1._data, [0, 4, 1, 4, 2, 4, 3, 4])
        # ArraySequence should maybe have the same degap method as normal seq
        self.assertEqual(str(model1.degap()), "UCAG")

    def test_the_rest_of_ModelSequence(self):
        """The class ArraySequence has 14 methods, but only 2 unittests.
        You might want to add some tests there..."""
        # note: mostly these are tested in derived classes, for convenience.
        pass


if __name__ == "__main__":
    main()
