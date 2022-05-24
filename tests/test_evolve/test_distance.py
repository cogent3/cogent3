#!/usr/bin/env python
import os
import warnings

from unittest import TestCase, main

import numpy

from numpy.testing import assert_allclose, assert_equal

from cogent3 import (
    DNA,
    PROTEIN,
    RNA,
    load_aligned_seqs,
    make_aligned_seqs,
    make_unaligned_seqs,
)
from cogent3.evolve.distance import EstimateDistances
from cogent3.evolve.fast_distance import (
    DistanceMatrix,
    HammingPair,
    JC69Pair,
    LogDetPair,
    ParalinearPair,
    PercentIdentityPair,
    TN93Pair,
    _calculators,
    _fill_diversity_matrix,
    _hamming,
    _jc69_from_matrix,
    available_distances,
    get_distance_calculator,
    get_moltype_index_array,
    seq_to_indices,
)
from cogent3.evolve.models import F81, HKY85, JC69
from cogent3.evolve.pairwise_distance_numba import (
    fill_diversity_matrix as numba_fill_diversity_matrix,
)


warnings.filterwarnings("ignore", "Not using MPI as mpi4py not found")


# hides the warning from taking log of -ve determinant
numpy.seterr(invalid="ignore")


__author__ = "Gavin Huttley, Yicheng Zhu and Ben Kaehler"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Yicheng Zhu", "Ben Kaehler"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class TestPair(TestCase):
    dna_char_indices = get_moltype_index_array(DNA)
    rna_char_indices = get_moltype_index_array(RNA)
    alignment = make_aligned_seqs(
        data=[("s1", "ACGTACGTAC"), ("s2", "GTGTACGTAC")], moltype=DNA
    )

    ambig_alignment = make_aligned_seqs(
        data=[("s1", "RACGTACGTACN"), ("s2", "AGTGTACGTACA")], moltype=DNA
    )

    diff_alignment = make_aligned_seqs(
        data=[("s1", "ACGTACGTTT"), ("s2", "GTGTACGTAC")], moltype=DNA
    )

    def test_char_to_index(self):
        """should correctly recode a DNA & RNA seqs into indices"""
        seq = "TCAGRNY?-"
        expected = [0, 1, 2, 3, -9, -9, -9, -9, -9]
        indices = seq_to_indices(seq, self.dna_char_indices)
        assert_equal(indices, expected)
        seq = "UCAGRNY?-"
        indices = seq_to_indices(seq, self.rna_char_indices)
        assert_equal(indices, expected)

    def test_fill_diversity_matrix_all(self):
        """make correct diversity matrix when all chars valid"""
        s1 = seq_to_indices("ACGTACGTAC", self.dna_char_indices)
        s2 = seq_to_indices("GTGTACGTAC", self.dna_char_indices)
        matrix = numpy.zeros((4, 4), float)
        # self-self should just be an identity matrix
        _fill_diversity_matrix(matrix, s1, s1)
        assert_equal(matrix.sum(), len(s1))
        assert_equal(
            matrix,
            numpy.array(
                [[2, 0, 0, 0], [0, 3, 0, 0], [0, 0, 3, 0], [0, 0, 0, 2]], float
            ),
        )

        # small diffs
        matrix.fill(0)
        _fill_diversity_matrix(matrix, s1, s2)
        assert_equal(
            matrix,
            numpy.array(
                [[2, 0, 0, 0], [1, 2, 0, 0], [0, 0, 2, 1], [0, 0, 0, 2]], float
            ),
        )

    def test_fill_diversity_matrix_some(self):
        """make correct diversity matrix when not all chars valid"""
        s1 = seq_to_indices("RACGTACGTACN", self.dna_char_indices)
        s2 = seq_to_indices("AGTGTACGTACA", self.dna_char_indices)
        matrix = numpy.zeros((4, 4), float)
        # small diffs
        matrix.fill(0)
        _fill_diversity_matrix(matrix, s1, s2)
        assert_equal(
            matrix,
            numpy.array(
                [[2, 0, 0, 0], [1, 2, 0, 0], [0, 0, 2, 1], [0, 0, 0, 2]], float
            ),
        )

    def test_python_vs_numba_fill_matrix(self):
        """python & cython fill_diversity_matrix give same answer"""
        s1 = seq_to_indices("RACGTACGTACN", self.dna_char_indices)
        s2 = seq_to_indices("AGTGTACGTACA", self.dna_char_indices)
        matrix1 = numpy.zeros((4, 4), float)
        _fill_diversity_matrix(matrix1, s1, s2)
        matrix2 = numpy.zeros((4, 4), float)
        numba_fill_diversity_matrix(matrix2, s1, s2)
        assert_allclose(matrix1, matrix2)

    def test_hamming_from_matrix(self):
        """compute hamming from diversity matrix"""
        s1 = seq_to_indices("ACGTACGTAC", self.dna_char_indices)
        s2 = seq_to_indices("GTGTACGTAC", self.dna_char_indices)
        matrix = numpy.zeros((4, 4), float)
        _fill_diversity_matrix(matrix, s1, s2)
        total, p, dist, var = _hamming(matrix)
        self.assertEqual(total, 10.0)
        self.assertEqual(dist, 2)
        self.assertEqual(p, 0.2)

    def test_hamming_pair(self):
        """get distances dict"""
        calc = HammingPair(DNA, alignment=self.alignment)
        calc.run(show_progress=False)
        dists = calc.get_pairwise_distances()
        dists = dists.to_dict()
        dist = 2.0
        expect = {("s1", "s2"): dist, ("s2", "s1"): dist}
        self.assertEqual(list(dists.keys()), list(expect.keys()))
        assert_allclose(list(dists.values()), list(expect.values()))

    def test_percent_pair(self):
        """get distances dict"""
        calc = PercentIdentityPair(DNA, alignment=self.alignment)
        calc.run(show_progress=False)
        dists = calc.get_pairwise_distances()
        dists = dists.to_dict()
        dist = 0.2
        expect = {("s1", "s2"): dist, ("s2", "s1"): dist}
        self.assertEqual(list(dists.keys()), list(expect.keys()))
        assert_allclose(list(dists.values()), list(expect.values()))

    def test_jc69_from_matrix(self):
        """compute JC69 from diversity matrix"""
        s1 = seq_to_indices("ACGTACGTAC", self.dna_char_indices)
        s2 = seq_to_indices("GTGTACGTAC", self.dna_char_indices)
        matrix = numpy.zeros((4, 4), float)
        _fill_diversity_matrix(matrix, s1, s2)
        total, p, dist, var = _jc69_from_matrix(matrix)
        self.assertEqual(total, 10.0)
        self.assertEqual(p, 0.2)

    def test_wrong_moltype(self):
        """specifying wrong moltype raises ValueError"""
        with self.assertRaises(ValueError):
            _ = JC69Pair(PROTEIN, alignment=self.alignment)

    def test_jc69_from_alignment(self):
        """compute JC69 dists from an alignment"""
        calc = JC69Pair(DNA, alignment=self.alignment)
        calc.run(show_progress=False)
        self.assertEqual(calc.lengths["s1", "s2"], 10)
        self.assertEqual(calc.proportions["s1", "s2"], 0.2)
        # value from OSX MEGA 5
        assert_allclose(calc.dists["s1", "s2"], 0.2326161962)
        # value**2 from OSX MEGA 5
        assert_allclose(calc.variances["s1", "s2"], 0.029752066125078681)
        # value from OSX MEGA 5
        assert_allclose(calc.stderr["s1", "s2"], 0.1724878724)

        # same answer when using ambiguous alignment
        calc.run(self.ambig_alignment, show_progress=False)
        assert_allclose(calc.dists["s1", "s2"], 0.2326161962)

        # but different answer if subsequent alignment is different
        calc.run(self.diff_alignment, show_progress=False)
        self.assertTrue(calc.dists["s1", "s2"] != 0.2326161962)

    def test_tn93_from_matrix(self):
        """compute TN93 distances"""
        calc = TN93Pair(DNA, alignment=self.alignment)
        calc.run(show_progress=False)
        self.assertEqual(calc.lengths["s1", "s2"], 10)
        self.assertEqual(calc.proportions["s1", "s2"], 0.2)
        # value from OSX MEGA 5
        assert_allclose(calc.dists["s1", "s2"], 0.2554128119)
        # value**2 from OSX MEGA 5
        assert_allclose(calc.variances["s1", "s2"], 0.04444444445376601)
        # value from OSX MEGA 5
        assert_allclose(calc.stderr["s1", "s2"], 0.2108185107)

        # same answer when using ambiguous alignment
        calc.run(self.ambig_alignment, show_progress=False)
        assert_allclose(calc.dists["s1", "s2"], 0.2554128119)

        # but different answer if subsequent alignment is different
        calc.run(self.diff_alignment, show_progress=False)
        self.assertTrue(calc.dists["s1", "s2"] != 0.2554128119)

    def test_distance_pair(self):
        """get distances dict"""
        calc = TN93Pair(DNA, alignment=self.alignment)
        calc.run(show_progress=False)
        dists = calc.get_pairwise_distances()
        dists = dists.to_dict()
        dist = 0.2554128119
        expect = {("s1", "s2"): dist, ("s2", "s1"): dist}
        self.assertEqual(list(dists.keys()), list(expect.keys()))
        assert_allclose(list(dists.values()), list(expect.values()))

    def test_logdet_pair_dna(self):
        """logdet should produce distances that match MEGA"""
        aln = load_aligned_seqs("data/brca1_5.paml", moltype=DNA)
        logdet_calc = LogDetPair(moltype=DNA, alignment=aln)
        logdet_calc.run(use_tk_adjustment=True, show_progress=False)
        dists = logdet_calc.get_pairwise_distances().to_dict()
        all_expected = {
            ("Human", "NineBande"): 0.075336929999999996,
            ("NineBande", "DogFaced"): 0.0898575452,
            ("DogFaced", "Human"): 0.1061747919,
            ("HowlerMon", "DogFaced"): 0.0934480008,
            ("Mouse", "HowlerMon"): 0.26422862920000001,
            ("NineBande", "Human"): 0.075336929999999996,
            ("HowlerMon", "NineBande"): 0.062202897899999998,
            ("DogFaced", "NineBande"): 0.0898575452,
            ("DogFaced", "HowlerMon"): 0.0934480008,
            ("Human", "DogFaced"): 0.1061747919,
            ("Mouse", "Human"): 0.26539976700000001,
            ("NineBande", "HowlerMon"): 0.062202897899999998,
            ("HowlerMon", "Human"): 0.036571181899999999,
            ("DogFaced", "Mouse"): 0.2652555144,
            ("HowlerMon", "Mouse"): 0.26422862920000001,
            ("Mouse", "DogFaced"): 0.2652555144,
            ("NineBande", "Mouse"): 0.22754789210000001,
            ("Mouse", "NineBande"): 0.22754789210000001,
            ("Human", "Mouse"): 0.26539976700000001,
            ("Human", "HowlerMon"): 0.036571181899999999,
        }
        for pair in dists:
            got = dists[pair]
            expected = all_expected[pair]
            assert_allclose(got, expected)

    def test_slice_dmatrix(self):
        data = {
            ("ABAYE2984", "Atu3667"): 0.25,
            ("ABAYE2984", "Avin_42730"): 0.638,
            ("ABAYE2984", "BAA10469"): None,
            ("Atu3667", "ABAYE2984"): 0.25,
            ("Atu3667", "Avin_42730"): 2.368,
            ("Atu3667", "BAA10469"): 0.25,
            ("Avin_42730", "ABAYE2984"): 0.638,
            ("Avin_42730", "Atu3667"): 2.368,
            ("Avin_42730", "BAA10469"): 1.85,
            ("BAA10469", "ABAYE2984"): 0.25,
            ("BAA10469", "Atu3667"): 0.25,
            ("BAA10469", "Avin_42730"): 1.85,
        }
        darr = DistanceMatrix(data)
        names = darr.template.names[0][:3]
        got = darr[:3, :3]
        self.assertEqual(list(got.template.names[0]), names)

    def test_logdet_tk_adjustment(self):
        """logdet using tamura kumar differs from classic"""
        aln = load_aligned_seqs("data/brca1_5.paml", moltype=DNA)
        logdet_calc = LogDetPair(moltype=DNA, alignment=aln)
        logdet_calc.run(use_tk_adjustment=True, show_progress=False)
        tk = logdet_calc.get_pairwise_distances()
        logdet_calc.run(use_tk_adjustment=False, show_progress=False)
        not_tk = logdet_calc.get_pairwise_distances()
        self.assertNotEqual(tk, not_tk)

    def test_logdet_pair_aa(self):
        """logdet shouldn't fail to produce distances for aa seqs"""
        aln = load_aligned_seqs("data/brca1_5.paml", moltype=DNA)
        aln = aln.get_translation()
        logdet_calc = LogDetPair(moltype=PROTEIN, alignment=aln)
        logdet_calc.run(use_tk_adjustment=True, show_progress=False)
        logdet_calc.get_pairwise_distances()

    def test_logdet_missing_states(self):
        """should calculate logdet measurement with missing states"""
        data = [
            (
                "seq1",
                "GGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
            ),
            (
                "seq2",
                "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTNTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
            ),
        ]
        aln = make_aligned_seqs(data=data, moltype=DNA)
        logdet_calc = LogDetPair(moltype=DNA, alignment=aln)
        logdet_calc.run(use_tk_adjustment=True, show_progress=False)

        dists = logdet_calc.get_pairwise_distances().to_dict()
        self.assertTrue(list(dists.values())[0] is not None)

        logdet_calc.run(use_tk_adjustment=False, show_progress=False)
        dists = logdet_calc.get_pairwise_distances().to_dict()
        self.assertTrue(list(dists.values())[0] is not None)

    def test_logdet_variance(self):
        """calculate logdet variance consistent with hand calculation"""
        data = [
            (
                "seq1",
                "GGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
            ),
            (
                "seq2",
                "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
            ),
        ]
        aln = make_aligned_seqs(data=data, moltype=DNA)
        logdet_calc = LogDetPair(moltype=DNA, alignment=aln)
        logdet_calc.run(use_tk_adjustment=True, show_progress=False)
        self.assertEqual(logdet_calc.variances[1, 1], None)

        index = dict(list(zip("ACGT", list(range(4)))))
        J = numpy.zeros((4, 4))
        for p in zip(data[0][1], data[1][1]):
            J[index[p[0]], index[p[1]]] += 1
        for i in range(4):
            if J[i, i] == 0:
                J[i, i] += 0.5
        J /= J.sum()
        M = numpy.linalg.inv(J)
        var = 0.0
        for i in range(4):
            for j in range(4):
                var += M[j, i] ** 2 * J[i, j] - 1
        var /= 16 * len(data[0][1])

        logdet_calc.run(use_tk_adjustment=False, show_progress=False)
        logdet_calc.get_pairwise_distances()
        assert_allclose(logdet_calc.variances[1, 1], var, atol=1e-3)

    def test_logdet_for_determinant_lte_zero(self):
        """returns distance of None if the determinant is <= 0"""
        data = dict(
            seq1="AGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
            seq2="TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
        )
        aln = make_aligned_seqs(data=data, moltype=DNA)

        logdet_calc = LogDetPair(moltype=DNA, alignment=aln)
        logdet_calc.run(use_tk_adjustment=True, show_progress=False)
        dists = logdet_calc.get_pairwise_distances().to_dict()
        self.assertTrue(numpy.isnan(list(dists.values())[0]))
        logdet_calc.run(use_tk_adjustment=False, show_progress=False)
        dists = logdet_calc.get_pairwise_distances().to_dict()
        self.assertTrue(numpy.isnan(list(dists.values())[0]))

        # but raises ArithmeticError if told to
        logdet_calc = LogDetPair(moltype=DNA, alignment=aln, invalid_raises=True)
        with self.assertRaises(ArithmeticError):
            logdet_calc.run(use_tk_adjustment=True, show_progress=False)

    def test_paralinear_pair_aa(self):
        """paralinear shouldn't fail to produce distances for aa seqs"""
        aln = load_aligned_seqs("data/brca1_5.paml", moltype=DNA)
        aln = aln.get_translation()
        paralinear_calc = ParalinearPair(moltype=PROTEIN, alignment=aln)
        paralinear_calc.run(show_progress=False)
        paralinear_calc.get_pairwise_distances()

    def test_paralinear_distance(self):
        """calculate paralinear variance consistent with hand calculation"""
        data = [
            (
                "seq1",
                "GGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
            ),
            (
                "seq2",
                "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
            ),
        ]
        aln = make_aligned_seqs(data=data, moltype=DNA)
        paralinear_calc = ParalinearPair(moltype=DNA, alignment=aln)
        paralinear_calc.run(show_progress=False)

        index = dict(list(zip("ACGT", list(range(4)))))
        J = numpy.zeros((4, 4))
        for p in zip(data[0][1], data[1][1]):
            J[index[p[0]], index[p[1]]] += 1
        for i in range(4):
            if J[i, i] == 0:
                J[i, i] += 0.5
        J /= J.sum()
        numpy.linalg.inv(J)
        f = J.sum(1), J.sum(0)
        dist = -0.25 * numpy.log(
            numpy.linalg.det(J) / numpy.sqrt(f[0].prod() * f[1].prod())
        )

        assert_allclose(paralinear_calc.dists["seq1", "seq2"], dist)

    def test_paralinear_variance(self):
        """calculate paralinear variance consistent with hand calculation"""
        data = [
            (
                "seq1",
                "GGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
            ),
            (
                "seq2",
                "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
            ),
        ]
        aln = make_aligned_seqs(data=data, moltype=DNA)
        paralinear_calc = ParalinearPair(moltype=DNA, alignment=aln)
        paralinear_calc.run(show_progress=False)

        index = dict(list(zip("ACGT", list(range(4)))))
        J = numpy.zeros((4, 4))
        for p in zip(data[0][1], data[1][1]):
            J[index[p[0]], index[p[1]]] += 1
        for i in range(4):
            if J[i, i] == 0:
                J[i, i] += 0.5
        J /= J.sum()
        M = numpy.linalg.inv(J)
        f = J.sum(1), J.sum(0)
        var = 0.0
        for i in range(4):
            for j in range(4):
                var += M[j, i] ** 2 * J[i, j]
            var -= 1 / numpy.sqrt(f[0][i] * f[1][i])
        var /= 16 * len(data[0][1])

        assert_allclose(paralinear_calc.variances[1, 1], var, atol=1e-3)

    def test_paralinear_for_determinant_lte_zero(self):
        """returns distance of None if the determinant is <= 0"""
        data = dict(
            seq1="AGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
            seq2="TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
        )
        aln = make_aligned_seqs(data=data, moltype=DNA)

        paralinear_calc = ParalinearPair(moltype=DNA, alignment=aln)
        paralinear_calc.run(show_progress=False)
        dists = paralinear_calc.get_pairwise_distances().to_dict()
        self.assertTrue(numpy.isnan(list(dists.values())[0]))
        paralinear_calc.run(show_progress=False)
        dists = paralinear_calc.get_pairwise_distances().to_dict()
        self.assertTrue(numpy.isnan(list(dists.values())[0]))

    def test_paralinear_pair_dna(self):
        """calculate paralinear distance consistent with logdet distance"""
        data = [
            (
                "seq1",
                "TAATTCATTGGGACGTCGAATCCGGCAGTCCTGCCGCAAAAGCTTCCGGAATCGAATTTTGGCA",
            ),
            (
                "seq2",
                "AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGG",
            ),
        ]
        aln = make_aligned_seqs(data=data, moltype=DNA)
        paralinear_calc = ParalinearPair(moltype=DNA, alignment=aln)
        paralinear_calc.run(show_progress=False)
        logdet_calc = LogDetPair(moltype=DNA, alignment=aln)
        logdet_calc.run(show_progress=False)
        self.assertEqual(logdet_calc.dists[1, 1], paralinear_calc.dists[1, 1])
        self.assertEqual(paralinear_calc.variances[1, 1], logdet_calc.variances[1, 1])

    def test_duplicated(self):
        """correctly identifies duplicates"""

        def get_calc(data):
            aln = make_aligned_seqs(data=data, moltype=DNA)
            calc = ParalinearPair(moltype=DNA, alignment=aln)
            calc(show_progress=False)
            return calc

        # no duplicates
        data = [
            (
                "seq1",
                "GGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
            ),
            (
                "seq2",
                "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
            ),
        ]
        calc = get_calc(data)
        self.assertEqual(calc.duplicated, None)
        data = [
            (
                "seq1",
                "GGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
            ),
            (
                "seq2",
                "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
            ),
            (
                "seq3",
                "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
            ),
        ]
        calc = get_calc(data)
        self.assertTrue(
            {"seq2": ["seq3"]} == calc.duplicated
            or {"seq3": ["seq2"]} == calc.duplicated
        )
        # default to get all pairwise distances
        pwds = calc.get_pairwise_distances().to_dict()
        self.assertEqual(pwds[("seq2", "seq3")], 0.0)
        self.assertEqual(pwds[("seq2", "seq1")], pwds[("seq3", "seq1")])

        # only unique seqs when using include_duplicates=False

        pwds = calc.get_pairwise_distances(include_duplicates=False).to_dict()
        present = list(calc.duplicated.keys())[0]
        missing = calc.duplicated[present][0]
        self.assertEqual(set([(present, missing)]), set([("seq2", "seq3")]))
        self.assertTrue((present, "seq1") in pwds)
        self.assertFalse((missing, "seq1") in pwds)


class TestGetDisplayCalculators(TestCase):
    def test_get_calculator(self):
        """exercising getting specified calculator"""
        for key in _calculators:
            get_distance_calculator(key)
            get_distance_calculator(key.upper())

        with self.assertRaises(ValueError):
            get_distance_calculator("blahblah")

    def test_available_distances(self):
        """available_distances has correct content"""
        content = available_distances()
        self.assertEqual(content.shape, (6, 2))
        self.assertEqual(content["tn93", 1], "dna, rna")


class TestDistanceMatrix(TestCase):
    def test_to_dict(self):
        """distance matrix correctly produces a 1D dict"""
        data = {("s1", "s2"): 0.25, ("s2", "s1"): 0.25}
        dmat = DistanceMatrix(data)
        got = dmat.to_dict()
        self.assertEqual(got, data)

    def test_matrix_dtype(self):
        """tests DistanceMatrix correctly accepts the data with proper dtype"""
        data = {
            ("ABAYE2984", "Atu3667"): None,
            ("ABAYE2984", "Avin_42730"): 0.638,
            ("ABAYE2984", "BAA10469"): None,
            ("Atu3667", "ABAYE2984"): None,
            ("Atu3667", "Avin_42730"): 2.368,
            ("Atu3667", "BAA10469"): None,
            ("Avin_42730", "ABAYE2984"): 0.638,
            ("Avin_42730", "Atu3667"): 2.368,
            ("Avin_42730", "BAA10469"): 1.85,
            ("BAA10469", "ABAYE2984"): None,
            ("BAA10469", "Atu3667"): None,
            ("BAA10469", "Avin_42730"): 1.85,
        }
        names = set()
        for p in data:
            names.update(p)

        # tests when data has None values and DistanceMatrix using dtype('float')
        darr = DistanceMatrix(data)
        self.assertEqual(darr.shape, (4, 4))
        self.assertEqual(set(darr.names), names)
        for (a, b), dist in data.items():
            if dist is None:
                assert numpy.isnan(darr[a, b])
            else:
                assert_allclose(dist, darr[a, b])

        data = {
            ("ABAYE2984", "Atu3667"): "None",
            ("ABAYE2984", "Avin_42730"): 0.638,
            ("ABAYE2984", "BAA10469"): None,
            ("Atu3667", "ABAYE2984"): None,
            ("Atu3667", "Avin_42730"): 2.368,
            ("Atu3667", "BAA10469"): "None",
            ("Avin_42730", "ABAYE2984"): 0.638,
            ("Avin_42730", "Atu3667"): 2.368,
            ("Avin_42730", "BAA10469"): 1.85,
            ("BAA10469", "ABAYE2984"): None,
            ("BAA10469", "Atu3667"): None,
            ("BAA10469", "Avin_42730"): 1.85,
        }

        # tests when data has str values and DistanceMatrix using dtype('float')
        with self.assertRaises(ValueError):
            darr = DistanceMatrix(data)

    def test_dropping_from_matrix(self):
        """pairwise distances should have method for dropping invalid data"""
        data = {
            ("ABAYE2984", "Atu3667"): None,
            ("ABAYE2984", "Avin_42730"): 0.638,
            ("ABAYE2984", "BAA10469"): None,
            ("Atu3667", "ABAYE2984"): None,
            ("Atu3667", "Avin_42730"): 2.368,
            ("Atu3667", "BAA10469"): None,
            ("Avin_42730", "ABAYE2984"): 0.638,
            ("Avin_42730", "Atu3667"): 2.368,
            ("Avin_42730", "BAA10469"): 1.85,
            ("BAA10469", "ABAYE2984"): None,
            ("BAA10469", "Atu3667"): None,
            ("BAA10469", "Avin_42730"): 1.85,
        }

        darr = DistanceMatrix(data)
        new = darr.drop_invalid()
        self.assertEqual(new, None)

        data = {
            ("ABAYE2984", "Atu3667"): 0.25,
            ("ABAYE2984", "Avin_42730"): 0.638,
            ("ABAYE2984", "BAA10469"): None,
            ("Atu3667", "ABAYE2984"): 0.25,
            ("Atu3667", "Avin_42730"): 2.368,
            ("Atu3667", "BAA10469"): 0.25,
            ("Avin_42730", "ABAYE2984"): 0.638,
            ("Avin_42730", "Atu3667"): 2.368,
            ("Avin_42730", "BAA10469"): 1.85,
            ("BAA10469", "ABAYE2984"): 0.25,
            ("BAA10469", "Atu3667"): 0.25,
            ("BAA10469", "Avin_42730"): 1.85,
        }
        darr = DistanceMatrix(data)
        new = darr.drop_invalid()
        self.assertEqual(new.shape, (2, 2))

    def test_take_dists(self):
        """subsets the distance matrix"""
        data = {
            ("ABAYE2984", "Atu3667"): 0.25,
            ("ABAYE2984", "Avin_42730"): 0.638,
            ("ABAYE2984", "BAA10469"): None,
            ("Atu3667", "ABAYE2984"): 0.25,
            ("Atu3667", "Avin_42730"): 2.368,
            ("Atu3667", "BAA10469"): 0.25,
            ("Avin_42730", "ABAYE2984"): 0.638,
            ("Avin_42730", "Atu3667"): 2.368,
            ("Avin_42730", "BAA10469"): 1.85,
            ("BAA10469", "ABAYE2984"): 0.25,
            ("BAA10469", "Atu3667"): 0.25,
            ("BAA10469", "Avin_42730"): 1.85,
        }
        darr = DistanceMatrix(data)
        got1 = darr.take_dists(["ABAYE2984", "Atu3667", "Avin_42730"])
        got2 = darr.take_dists("BAA10469", negate=True)
        assert_allclose(got1.array.astype(float), got2.array.astype(float))

    def test_build_phylogeny(self):
        """build a NJ tree"""
        from cogent3 import make_tree

        dists = {
            ("DogFaced", "FlyingFox"): 0.05,
            ("DogFaced", "FreeTaile"): 0.14,
            ("DogFaced", "LittleBro"): 0.16,
            ("DogFaced", "TombBat"): 0.15,
            ("FlyingFox", "DogFaced"): 0.05,
            ("FlyingFox", "FreeTaile"): 0.12,
            ("FlyingFox", "LittleBro"): 0.13,
            ("FlyingFox", "TombBat"): 0.14,
            ("FreeTaile", "DogFaced"): 0.14,
            ("FreeTaile", "FlyingFox"): 0.12,
            ("FreeTaile", "LittleBro"): 0.09,
            ("FreeTaile", "TombBat"): 0.1,
            ("LittleBro", "DogFaced"): 0.16,
            ("LittleBro", "FlyingFox"): 0.13,
            ("LittleBro", "FreeTaile"): 0.09,
            ("LittleBro", "TombBat"): 0.12,
            ("TombBat", "DogFaced"): 0.15,
            ("TombBat", "FlyingFox"): 0.14,
            ("TombBat", "FreeTaile"): 0.1,
            ("TombBat", "LittleBro"): 0.12,
        }
        dists = DistanceMatrix(dists)
        got = dists.quick_tree(show_progress=False)
        expect = make_tree(
            treestring="((TombBat,(DogFaced,FlyingFox)),LittleBro,FreeTaile)"
        )
        self.assertTrue(expect.same_topology(got))

    def test_names(self):
        """names property works"""
        data = {
            ("ABAYE2984", "Atu3667"): 0.25,
            ("ABAYE2984", "Avin_42730"): 0.638,
            ("ABAYE2984", "BAA10469"): None,
            ("Atu3667", "ABAYE2984"): 0.25,
            ("Atu3667", "Avin_42730"): 2.368,
            ("Atu3667", "BAA10469"): 0.25,
            ("Avin_42730", "ABAYE2984"): 0.638,
            ("Avin_42730", "Atu3667"): 2.368,
            ("Avin_42730", "BAA10469"): 1.85,
            ("BAA10469", "ABAYE2984"): 0.25,
            ("BAA10469", "Atu3667"): 0.25,
            ("BAA10469", "Avin_42730"): 1.85,
        }
        names = set()
        for p in data:
            names.update(p)
        darr = DistanceMatrix(data)
        self.assertEqual(set(darr.names), names)
        darr = darr.drop_invalid()
        for n in ("ABAYE2984", "BAA10469"):
            names.remove(n)
        self.assertEqual(set(darr.names), names)


class DistancesTests(TestCase):
    def setUp(self):
        self.al = make_aligned_seqs(
            data={
                "a": "GTACGTACGATC",
                "b": "GTACGTACGTAC",
                "c": "GTACGTACGTTC",
                "e": "GTACGTACTGGT",
            }
        )
        self.collection = make_unaligned_seqs(
            data={
                "a": "GTACGTACGATC",
                "b": "GTACGTACGTAC",
                "c": "GTACGTACGTTC",
                "e": "GTACGTACTGGT",
            }
        )

    def assertDistsAlmostEqual(self, expected, observed, precision=4):
        observed = dict([(frozenset(k), v) for (k, v) in list(observed.items())])
        expected = dict([(frozenset(k), v) for (k, v) in list(expected.items())])
        for key in expected:
            self.assertAlmostEqual(expected[key], observed[key], precision)

    def test_EstimateDistances(self):
        """testing (well, exercising at least), EstimateDistances"""
        d = EstimateDistances(self.al, JC69())
        d.run(show_progress=False)
        canned_result = {
            ("b", "e"): 0.440840,
            ("c", "e"): 0.440840,
            ("a", "c"): 0.088337,
            ("a", "b"): 0.188486,
            ("a", "e"): 0.440840,
            ("b", "c"): 0.0883373,
        }
        result = d.get_pairwise_distances().to_dict()
        self.assertDistsAlmostEqual(canned_result, result)

        # excercise writing to file
        d.write("junk.txt")
        try:
            os.remove("junk.txt")
        except OSError:
            pass  # probably parallel

    def test_EstimateDistancesWithMotifProbs(self):
        """EstimateDistances with supplied motif probs"""
        motif_probs = {"A": 0.1, "C": 0.2, "G": 0.2, "T": 0.5}
        d = EstimateDistances(self.al, HKY85(), motif_probs=motif_probs)
        d.run(show_progress=False)
        canned_result = {
            ("a", "c"): 0.07537,
            ("b", "c"): 0.07537,
            ("a", "e"): 0.39921,
            ("a", "b"): 0.15096,
            ("b", "e"): 0.39921,
            ("c", "e"): 0.37243,
        }
        result = d.get_pairwise_distances().to_dict()
        self.assertDistsAlmostEqual(canned_result, result)

    def test_EstimateDistances_fromThreeway(self):
        """testing (well, exercising at least), EsimateDistances fromThreeway"""
        d = EstimateDistances(self.al, JC69(), threeway=True)
        d.run(show_progress=False)
        canned_result = {
            ("b", "e"): 0.495312,
            ("c", "e"): 0.479380,
            ("a", "c"): 0.089934,
            ("a", "b"): 0.190021,
            ("a", "e"): 0.495305,
            ("b", "c"): 0.0899339,
        }
        result = d.get_pairwise_distances(summary_function="mean").to_dict()
        self.assertDistsAlmostEqual(canned_result, result)

    def test_EstimateDistances_fromUnaligned(self):
        """Excercising estimate distances from unaligned sequences"""
        d = EstimateDistances(
            self.collection, JC69(), do_pair_align=True, rigorous_align=True
        )
        d.run(show_progress=False)
        canned_result = {
            ("b", "e"): 0.440840,
            ("c", "e"): 0.440840,
            ("a", "c"): 0.088337,
            ("a", "b"): 0.188486,
            ("a", "e"): 0.440840,
            ("b", "c"): 0.0883373,
        }
        result = d.get_pairwise_distances().to_dict()
        self.assertDistsAlmostEqual(canned_result, result)

        d = EstimateDistances(
            self.collection, JC69(), do_pair_align=True, rigorous_align=False
        )
        d.run(show_progress=False)
        canned_result = {
            ("b", "e"): 0.440840,
            ("c", "e"): 0.440840,
            ("a", "c"): 0.088337,
            ("a", "b"): 0.188486,
            ("a", "e"): 0.440840,
            ("b", "c"): 0.0883373,
        }
        result = d.get_pairwise_distances().to_dict()
        self.assertDistsAlmostEqual(canned_result, result)

    def test_EstimateDistances_other_model_params(self):
        """test getting other model params from EstimateDistances"""
        d = EstimateDistances(self.al, HKY85(), est_params=["kappa"])
        d.run(show_progress=False)
        # this will be a Number object with Mean, Median etc ..
        kappa = d.get_param_values("kappa")
        self.assertAlmostEqual(kappa.mean, 0.8939, 4)
        # this will be a dict with pairwise instances, it's called by the above
        # method, so the correctness of it's values is already checked
        kappa = d.get_pairwise_param("kappa")

    def test_EstimateDistances_modify_lf(self):
        """tests modifying the lf"""

        def constrain_fit(lf):
            lf.set_param_rule("kappa", is_constant=True)
            lf.optimise(local=True)
            return lf

        d = EstimateDistances(self.al, HKY85(), modify_lf=constrain_fit)
        d.run(show_progress=False)
        result = d.get_pairwise_distances().to_dict()
        d = EstimateDistances(self.al, F81())
        d.run(show_progress=False)
        expect = d.get_pairwise_distances().to_dict()
        self.assertDistsAlmostEqual(expect, result)

    def test_get_raw_estimates(self):
        """correctly return raw result object"""
        d = EstimateDistances(self.al, HKY85(), est_params=["kappa"])
        d.run(show_progress=False)
        expect = {
            ("a", "b"): {
                "kappa": 1.0000226766004808e-06,
                "length": 0.18232155856115662,
            },
            ("a", "c"): {
                "kappa": 1.0010380037049357e-06,
                "length": 0.087070406623635604,
            },
            ("a", "e"): {"kappa": 2.3965871843412687, "length": 0.4389176272584539},
            ("b", "e"): {"kappa": 2.3965871854366592, "length": 0.43891762729173389},
            ("b", "c"): {
                "kappa": 1.0010380037049357e-06,
                "length": 0.087070406623635604,
            },
            ("c", "e"): {"kappa": 0.57046787478038707, "length": 0.43260232210282784},
        }
        got = d.get_all_param_values()
        for pair in expect:
            for param in expect[pair]:
                self.assertAlmostEqual(got[pair][param], expect[pair][param], places=6)

    def test_no_calc(self):
        """returns None if no calculation done"""
        al = load_aligned_seqs("data/brca1_5.paml")
        d = EstimateDistances(al, submodel=HKY85())
        self.assertEqual(d.get_pairwise_distances(), None)

    def test_to_table(self):
        """converts a distance matrix to a Table"""
        data = {
            ("A", "B"): 2,
            ("A", "C"): 3,
            ("B", "C"): 1,
            ("B", "A"): 2,
            ("C", "A"): 3,
            ("C", "B"): 1,
        }
        darr = DistanceMatrix(data)
        table = darr.to_table()
        self.assertEqual(table.shape, (3, 4))
        self.assertEqual(table.columns["names"].tolist(), list(darr.names))
        self.assertEqual(table["A", "B"], 2)
        self.assertEqual(table["A", "A"], 0)


if __name__ == "__main__":
    main()
