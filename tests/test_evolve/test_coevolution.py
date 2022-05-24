from os import environ, remove
from tempfile import NamedTemporaryFile, mktemp
from unittest import TestCase, main

from numpy import (
    arange,
    array,
    e,
    greater_equal,
    less_equal,
    log,
    nan,
    sqrt,
    zeros,
)

from cogent3 import (
    DNA,
    PROTEIN,
    RNA,
    load_aligned_seqs,
    make_aligned_seqs,
    make_tree,
)
from cogent3.core.alignment import ArrayAlignment
from cogent3.core.alphabet import CharAlphabet
from cogent3.evolve.coevolution import (
    DEFAULT_NULL_VALUE,
    AAGapless,
    aln_position_pairs_cmp_threshold,
    aln_position_pairs_ge_threshold,
    aln_position_pairs_le_threshold,
    ancestral_state_alignment,
    ancestral_state_pair,
    ancestral_state_position,
    ancestral_states_input_validation,
    build_coevolution_matrix_filepath,
    calc_pair_scale,
    coevolution_matrix_to_csv,
    coevolve_alignment,
    coevolve_alignments,
    coevolve_alignments_validation,
    coevolve_pair,
    coevolve_position,
    count_ge_threshold,
    count_le_threshold,
    csv_to_coevolution_matrix,
    filter_exclude_positions,
    filter_non_parsimony_informative,
    filter_threshold_based_multiple_interdependency,
    freqs_from_aln,
    freqs_to_array,
    get_allowed_perturbations,
    get_ancestral_seqs,
    get_dg,
    get_dgg,
    get_positional_frequencies,
    get_positional_probabilities,
    get_subalignments,
    identify_aln_positions_above_threshold,
    ignore_excludes,
    is_parsimony_informative,
    join_positions,
    ltm_to_symmetric,
    make_weights,
    merge_alignments,
    mi,
    mi_alignment,
    mi_pair,
    mi_position,
    n_random_seqs,
    nmi,
    nmi_alignment,
    nmi_pair,
    nmi_position,
    normalized_mi,
    parse_coevolution_matrix_filepath,
    pickle_coevolution_result,
    probs_from_dict,
    protein_dict,
    resampled_mi_alignment,
    sca_alignment,
    sca_input_validation,
    sca_pair,
    sca_position,
    unpickle_coevolution_result,
    validate_alignment,
    validate_alphabet,
    validate_ancestral_seqs,
    validate_position,
    validate_tree,
)
from cogent3.maths.stats.number import CategoryCounter


__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Greg Caporaso"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Beta"

from numpy.testing import assert_allclose, assert_equal


class CoevolutionTests(TestCase):
    """Tests of coevolution.py"""

    def setUp(self):
        """Set up variables for us in tests"""
        self.run_slow_tests = int(environ.get("TEST_SLOW_APPC", 0))
        # Data used in SCA tests
        self.dna_aln = ArrayAlignment(
            data=list(zip(list(range(4)), ["ACGT", "AGCT", "ACCC", "TAGG"])),
            moltype=DNA,
        )
        self.rna_aln = ArrayAlignment(
            data=list(zip(list(range(4)), ["ACGU", "AGCU", "ACCC", "UAGG"])),
            moltype=RNA,
        )
        self.protein_aln = ArrayAlignment(
            data=list(zip(list(range(4)), ["ACGP", "AGCT", "ACCC", "TAGG"])),
            moltype=PROTEIN,
        )
        self.dna_aln_gapped = ArrayAlignment(
            data=list(zip(list(range(4)), ["A-CGT", "AGC-T", "-ACCC", "TAGG-"])),
            moltype=DNA,
        )
        self.freq = ArrayAlignment(
            data=list(
                zip(
                    list(range(20)),
                    [
                        "TCT",
                        "CCT",
                        "CCC",
                        "CCC",
                        "CCG",
                        "CC-",
                        "AC-",
                        "AC-",
                        "AA-",
                        "AA-",
                        "GA-",
                        "GA-",
                        "GA-",
                        "GA-",
                        "GA-",
                        "G--",
                        "G--",
                        "G--",
                        "G--",
                        "G--",
                    ],
                )
            ),
            moltype=PROTEIN,
        )
        self.two_pos = ArrayAlignment(
            data=list(
                zip(
                    list(map(str, list(range(20)))),
                    [
                        "TC",
                        "CC",
                        "CC",
                        "CC",
                        "CC",
                        "CC",
                        "AC",
                        "AC",
                        "AA",
                        "AA",
                        "GA",
                        "GA",
                        "GA",
                        "GA",
                        "GA",
                        "GT",
                        "GT",
                        "GT",
                        "GT",
                        "GT",
                    ],
                )
            ),
            moltype=PROTEIN,
        )
        self.tree20 = make_tree(treestring=tree20_string)
        self.gpcr_aln = gpcr_aln
        self.myos_aln = myos_aln
        # a made-up dict of base frequencies to use as the natural freqs
        # for SCA calcs on DNA seqs
        self.dna_base_freqs = dict(list(zip("ACGT", [0.25] * 4)))
        self.rna_base_freqs = dict(list(zip("ACGU", [0.25] * 4)))
        self.protein_aln4 = ArrayAlignment(
            [("A1", "AACF"), ("A12", "AADF"), ("A123", "ADCF"), ("A111", "AAD-")],
            moltype=PROTEIN,
        )
        self.rna_aln4 = ArrayAlignment(
            [("A1", "AAUU"), ("A12", "ACGU"), ("A123", "UUAA"), ("A111", "AAA-")],
            moltype=RNA,
        )
        self.dna_aln4 = ArrayAlignment(
            [("A1", "AATT"), ("A12", "ACGT"), ("A123", "TTAA"), ("A111", "AAA?")],
            moltype=DNA,
        )
        self.tree4 = make_tree(
            treestring="((A1:0.5,A111:0.5):0.5,(A12:0.5,A123:0.5):0.5);"
        )

    def test_alignment_analyses_moltype_protein(self):
        """alignment methods work with moltype = PROTEIN"""

        r = mi_alignment(self.protein_aln4)
        self.assertEqual(r.shape, (4, 4))
        r = nmi_alignment(self.protein_aln4)
        self.assertEqual(r.shape, (4, 4))
        r = sca_alignment(self.protein_aln4, cutoff=0.75)
        self.assertEqual(r.shape, (4, 4))

        r = ancestral_state_alignment(self.protein_aln4, self.tree4)
        self.assertEqual(r.shape, (4, 4))

    def test_alignment_analyses_moltype_rna(self):
        """alignment methods work with moltype = RNA"""

        r = mi_alignment(self.rna_aln4)
        self.assertEqual(r.shape, (4, 4))
        r = nmi_alignment(self.rna_aln4)
        self.assertEqual(r.shape, (4, 4))
        r = sca_alignment(
            self.rna_aln4,
            cutoff=0.75,
            alphabet="ACGU",
            background_freqs=self.rna_base_freqs,
        )
        self.assertEqual(r.shape, (4, 4))

        r = ancestral_state_alignment(self.rna_aln4, self.tree4)
        self.assertEqual(r.shape, (4, 4))

    def test_alignment_analyses_moltype_dna(self):
        """alignment methods work with moltype = DNA"""

        r = mi_alignment(self.dna_aln4)
        self.assertEqual(r.shape, (4, 4))
        r = nmi_alignment(self.dna_aln4)
        self.assertEqual(r.shape, (4, 4))
        r = sca_alignment(
            self.dna_aln4,
            cutoff=0.75,
            alphabet="ACGT",
            background_freqs=self.dna_base_freqs,
        )
        self.assertEqual(r.shape, (4, 4))

        r = ancestral_state_alignment(self.dna_aln4, self.tree4)
        self.assertEqual(r.shape, (4, 4))

    def test_join_positions(self):
        """join_positions functions as expected"""
        self.assertEqual(
            join_positions(list("ABCD"), list("WXYZ")), ["AW", "BX", "CY", "DZ"]
        )
        self.assertEqual(join_positions(list("AAA"), list("BBB")), ["AB", "AB", "AB"])
        self.assertEqual(join_positions([], []), [])

    def test_mi(self):
        """mi calculations function as expected with valid data"""
        assert_allclose(mi(1.0, 1.0, 1.0), 1.0)
        assert_allclose(mi(1.0, 1.0, 2.0), 0.0)
        assert_allclose(mi(1.0, 1.0, 1.5), 0.5)

    def test_normalized_mi(self):
        """normalized mi calculations function as expected with valid data"""
        assert_allclose(normalized_mi(1.0, 1.0, 1.0), 1.0)
        assert_allclose(normalized_mi(1.0, 1.0, 2.0), 0.0)
        assert_allclose(normalized_mi(1.0, 1.0, 1.5), 0.3333, 3)

    def test_mi_pair(self):
        """mi_pair calculates mi from a pair of columns"""
        aln = ArrayAlignment(data={"1": "AB", "2": "AB"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1), 0.0)
        aln = ArrayAlignment(data={"1": "AB", "2": "BA"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1), 1.0)
        # order of positions doesn't matter (when it shouldn't)
        aln = ArrayAlignment(data={"1": "AB", "2": "AB"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1), mi_pair(aln, pos1=1, pos2=0))
        aln = ArrayAlignment(data={"1": "AB", "2": "BA"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1), mi_pair(aln, pos1=1, pos2=0))

    def test_wrapper_functions_handle_invalid_parameters(self):
        """coevolve_*: functions error on missing parameters"""
        # missing cutoff
        aln = ArrayAlignment(data={"1": "AC", "2": "AC"}, moltype=PROTEIN)
        self.assertRaises(ValueError, coevolve_pair, sca_pair, aln, 0, 1)
        self.assertRaises(ValueError, coevolve_position, sca_position, aln, 0)
        self.assertRaises(ValueError, coevolve_alignment, sca_alignment, aln)
        self.assertRaises(ValueError, coevolve_alignments, sca_alignment, aln, aln)

    def test_coevolve_pair(self):
        """coevolve_pair: returns same as pair methods called directly"""
        aln = ArrayAlignment(data={"1": "AC", "2": "AC"}, moltype=PROTEIN)
        t = make_tree(treestring="(1:0.5,2:0.5);")
        cutoff = 0.50
        # mi_pair == coevolve_pair(mi_pair,...)
        assert_allclose(
            coevolve_pair(mi_pair, aln, pos1=0, pos2=1), mi_pair(aln, pos1=0, pos2=1)
        )
        assert_allclose(
            coevolve_pair(nmi_pair, aln, pos1=0, pos2=1), nmi_pair(aln, pos1=0, pos2=1)
        )
        assert_allclose(
            coevolve_pair(ancestral_state_pair, aln, pos1=0, pos2=1, tree=t),
            ancestral_state_pair(aln, pos1=0, pos2=1, tree=t),
        )
        assert_allclose(
            coevolve_pair(sca_pair, aln, pos1=0, pos2=1, cutoff=cutoff),
            sca_pair(aln, pos1=0, pos2=1, cutoff=cutoff),
        )

    def test_coevolve_position(self):
        """coevolve_position: returns same as position methods called directly"""
        aln = ArrayAlignment(data={"1": "AC", "2": "AC"}, moltype=PROTEIN)
        t = make_tree(treestring="(1:0.5,2:0.5);")
        cutoff = 0.50
        # mi_position == coevolve_position(mi_position,...)
        assert_allclose(
            coevolve_position(mi_position, aln, position=0),
            mi_position(aln, position=0),
        )
        assert_allclose(
            coevolve_position(nmi_position, aln, position=0),
            nmi_position(aln, position=0),
        )
        assert_allclose(
            coevolve_position(ancestral_state_position, aln, position=0, tree=t),
            ancestral_state_position(aln, position=0, tree=t),
        )
        assert_allclose(
            coevolve_position(sca_position, aln, position=0, cutoff=cutoff),
            sca_position(aln, position=0, cutoff=cutoff),
        )

    def test_coevolve_alignment(self):
        """coevolve_alignment: returns same as alignment methods"""
        aln = ArrayAlignment(data={"1": "AC", "2": "AC"}, moltype=PROTEIN)
        t = make_tree(treestring="(1:0.5,2:0.5);")
        cutoff = 0.50
        # mi_alignment == coevolve_alignment(mi_alignment,...)
        assert_allclose(coevolve_alignment(mi_alignment, aln), mi_alignment(aln))
        assert_allclose(coevolve_alignment(mip_alignment, aln), mip_alignment(aln))
        assert_allclose(coevolve_alignment(mia_alignment, aln), mia_alignment(aln))
        assert_allclose(coevolve_alignment(nmi_alignment, aln), nmi_alignment(aln))
        assert_allclose(
            coevolve_alignment(ancestral_state_alignment, aln, tree=t),
            ancestral_state_alignment(aln, tree=t),
        )
        assert_allclose(
            coevolve_alignment(sca_alignment, aln, cutoff=cutoff),
            sca_alignment(aln, cutoff=cutoff),
        )

    def test_coevolve_alignments_validation_idenifiers(self):
        """coevolve_alignments_validation: seq/tree validation functions"""
        method = sca_alignment
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AD"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1": "EFW", "2": "EGY"}, moltype=PROTEIN)
        t = make_tree(treestring="(1:0.5,2:0.5);")
        # OK w/ no tree
        coevolve_alignments_validation(method, aln1, aln2, 2, None)
        # OK w/ tree
        coevolve_alignments_validation(method, aln1, aln2, 2, None, tree=t)
        # If there is a plus present in identifiers, we only care about the
        # text before the colon
        aln1 = ArrayAlignment(data={"1+a": "AC", "2+b": "AD"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1 + c": "EFW", "2 + d": "EGY"}, moltype=PROTEIN)
        t = make_tree(treestring="(1+e:0.5,2 + f:0.5);")
        # OK w/ no tree
        coevolve_alignments_validation(method, aln1, aln2, 2, None)
        # OK w/ tree
        coevolve_alignments_validation(method, aln1, aln2, 2, None, tree=t)

        # mismatch b/w alignments seq names
        aln1 = ArrayAlignment(data={"3": "AC", "2": "AD"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1": "EFW", "2": "EGY"}, moltype=PROTEIN)
        t = make_tree(treestring="(1:0.5,2:0.5);")
        self.assertRaises(
            AssertionError,
            coevolve_alignments_validation,
            method,
            aln1,
            aln2,
            2,
            None,
            tree=t,
        )

        # mismatch b/w alignments and tree seq names
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AD"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1": "EFW", "2": "EGY"}, moltype=PROTEIN)
        t = make_tree(treestring="(3:0.5,2:0.5);")
        self.assertRaises(
            AssertionError,
            coevolve_alignments_validation,
            method,
            aln1,
            aln2,
            2,
            None,
            tree=t,
        )

        # mismatch b/w alignments in number of seqs
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AD", "3": "AA"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1": "EFW", "2": "EGY"}, moltype=PROTEIN)
        t = make_tree(treestring="(1:0.5,2:0.5);")
        self.assertRaises(
            AssertionError, coevolve_alignments_validation, method, aln1, aln2, 2, None
        )
        self.assertRaises(
            AssertionError,
            coevolve_alignments_validation,
            method,
            aln1,
            aln2,
            2,
            None,
            tree=t,
        )

        # mismatch b/w alignments & tree in number of seqs
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AD"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1": "EFW", "2": "EGY"}, moltype=PROTEIN)
        t = make_tree(treestring="(1:0.5,(2:0.5,3:0.25));")
        self.assertRaises(
            AssertionError,
            coevolve_alignments_validation,
            method,
            aln1,
            aln2,
            2,
            None,
            tree=t,
        )

    def test_coevolve_alignments_validation_min_num_seqs(self):
        """coevolve_alignments_validation: ValueError on fewer than min_num_seqs"""
        method = mi_alignment
        # too few sequences -> ValueError
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AD"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1": "EFW", "2": "EGY"}, moltype=PROTEIN)
        coevolve_alignments_validation(method, aln1, aln2, 1, None)
        coevolve_alignments_validation(method, aln1, aln2, 2, None)
        self.assertRaises(
            ValueError, coevolve_alignments_validation, method, aln1, aln2, 3, None
        )

    def test_coevolve_alignments_validation_max_num_seqs(self):
        """coevolve_alignments_validation: min_num_seqs <= max_num_seqs"""
        method = mi_alignment
        # min_num_seqs > max_num_seqs-> ValueError
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AD"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1": "EFW", "2": "EGY"}, moltype=PROTEIN)
        coevolve_alignments_validation(method, aln1, aln2, 1, None)
        coevolve_alignments_validation(method, aln1, aln2, 1, 3)
        coevolve_alignments_validation(method, aln1, aln2, 2, 3)
        self.assertRaises(
            ValueError, coevolve_alignments_validation, method, aln1, aln2, 3, 2
        )

    def test_coevolve_alignments_validation_moltypes(self):
        """coevolve_alignments_validation: valid for acceptable moltypes"""
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AU"}, moltype=RNA)
        aln2 = ArrayAlignment(data={"1": "EFW", "2": "EGY"}, moltype=PROTEIN)
        # different moltype
        coevolve_alignments_validation(mi_alignment, aln1, aln2, 2, None)
        coevolve_alignments_validation(nmi_alignment, aln1, aln2, 2, None)
        coevolve_alignments_validation(resampled_mi_alignment, aln1, aln2, 2, None)
        self.assertRaises(
            AssertionError,
            coevolve_alignments_validation,
            sca_alignment,
            aln1,
            aln2,
            2,
            None,
        )
        self.assertRaises(
            AssertionError,
            coevolve_alignments_validation,
            ancestral_state_alignment,
            aln1,
            aln2,
            2,
            None,
        )

    def test_coevolve_alignments(self):
        """coevolve_alignments: returns correct len(aln1) x len(aln2) matrix"""
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AD"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1": "EFW", "2": "EGY"}, moltype=PROTEIN)
        combined_aln = ArrayAlignment(
            data={"1": "ACEFW", "2": "ADEGY"}, moltype=PROTEIN
        )
        t = make_tree(treestring="(1:0.5,2:0.5);")
        cutoff = 0.50
        # MI
        m = mi_alignment(combined_aln)
        expected = array([[m[2, 0], m[2, 1]], [m[3, 0], m[3, 1]], [m[4, 0], m[4, 1]]])
        assert_allclose(coevolve_alignments(mi_alignment, aln1, aln2), expected)
        # MI (return_full=True)
        assert_allclose(
            coevolve_alignments(mi_alignment, aln1, aln2, return_full=True), m
        )
        # NMI
        m = nmi_alignment(combined_aln)
        expected = array([[m[2, 0], m[2, 1]], [m[3, 0], m[3, 1]], [m[4, 0], m[4, 1]]])
        assert_allclose(coevolve_alignments(nmi_alignment, aln1, aln2), expected)
        # AS
        m = ancestral_state_alignment(combined_aln, tree=t)
        expected = array([[m[2, 0], m[2, 1]], [m[3, 0], m[3, 1]], [m[4, 0], m[4, 1]]])
        assert_allclose(
            coevolve_alignments(ancestral_state_alignment, aln1, aln2, tree=t), expected
        )
        # SCA
        m = sca_alignment(combined_aln, cutoff=cutoff)
        expected = array([[m[2, 0], m[2, 1]], [m[3, 0], m[3, 1]], [m[4, 0], m[4, 1]]])
        assert_allclose(
            coevolve_alignments(sca_alignment, aln1, aln2, cutoff=cutoff), expected
        )

    def test_coevolve_alignments_watches_min_num_seqs(self):
        """coevolve_alignments: error on too few sequences"""
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AD"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1": "EFW", "2": "EGY"}, moltype=PROTEIN)

        coevolve_alignments(mi_alignment, aln1, aln2)
        coevolve_alignments(mi_alignment, aln1, aln2, min_num_seqs=0)
        coevolve_alignments(mi_alignment, aln1, aln2, min_num_seqs=1)
        coevolve_alignments(mi_alignment, aln1, aln2, min_num_seqs=2)
        self.assertRaises(
            ValueError, coevolve_alignments, mi_alignment, aln1, aln2, min_num_seqs=3
        )
        self.assertRaises(
            ValueError, coevolve_alignments, mi_alignment, aln1, aln2, min_num_seqs=50
        )

    def test_coevolve_alignments_watches_max_num_seqs(self):
        """coevolve_alignments: filtering or error on too many sequences"""
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AD", "3": "YP"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(
            data={"1": "ACP", "2": "EAD", "3": "PYP"}, moltype=PROTEIN
        )

        # keep all seqs
        tmp_filepath = NamedTemporaryFile(
            prefix="tmp_test_coevolution", suffix=".fasta"
        ).name
        coevolve_alignments(
            mi_alignment, aln1, aln2, max_num_seqs=3, merged_aln_filepath=tmp_filepath
        )
        self.assertEqual(load_aligned_seqs(tmp_filepath).num_seqs, 3)

        # keep 2 seqs
        coevolve_alignments(
            mi_alignment, aln1, aln2, max_num_seqs=2, merged_aln_filepath=tmp_filepath
        )
        seqs = load_aligned_seqs(tmp_filepath)
        self.assertEqual(seqs.num_seqs, 2)

        # error if no sequence filter
        self.assertRaises(
            ValueError,
            coevolve_alignments,
            mi_alignment,
            aln1,
            aln2,
            max_num_seqs=2,
            merged_aln_filepath=tmp_filepath,
            sequence_filter=None,
        )

        # clean up the temporary file
        remove(tmp_filepath)

    def test_coevolve_alignments_different_MolType(self):
        """coevolve_alignments: different MolTypes supported"""
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AU"}, moltype=RNA)
        aln2 = ArrayAlignment(data={"1": "EFW", "2": "EGY"}, moltype=PROTEIN)
        combined_aln = ArrayAlignment(data={"1": "ACEFW", "2": "AUEGY"})
        t = make_tree(treestring="(1:0.5,2:0.5);")
        # MI
        m = mi_alignment(combined_aln)
        expected = array([[m[2, 0], m[2, 1]], [m[3, 0], m[3, 1]], [m[4, 0], m[4, 1]]])
        assert_allclose(coevolve_alignments(mi_alignment, aln1, aln2), expected)
        # MI (return_full=True)
        assert_allclose(
            coevolve_alignments(mi_alignment, aln1, aln2, return_full=True), m
        )
        # NMI
        m = nmi_alignment(combined_aln)
        expected = array([[m[2, 0], m[2, 1]], [m[3, 0], m[3, 1]], [m[4, 0], m[4, 1]]])
        assert_allclose(coevolve_alignments(nmi_alignment, aln1, aln2), expected)

    def test_mi_pair_cols_default_exclude_handling(self):
        """mi_pair returns null_value on excluded by default"""
        aln = ArrayAlignment(data={"1": "AB", "2": "-B"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1), DEFAULT_NULL_VALUE)
        aln = ArrayAlignment(data={"1": "-B", "2": "-B"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1), DEFAULT_NULL_VALUE)
        aln = ArrayAlignment(data={"1": "AA", "2": "-B"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1), DEFAULT_NULL_VALUE)
        aln = ArrayAlignment(data={"1": "AA", "2": "PB"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1, excludes="P"), DEFAULT_NULL_VALUE)

    def test_mi_pair_cols_non_default_exclude_handling(self):
        """mi_pair uses non-default exclude_handler when provided"""
        aln = ArrayAlignment(data={"1": "A-", "2": "A-"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1), DEFAULT_NULL_VALUE)
        assert_allclose(
            mi_pair(aln, pos1=0, pos2=1, exclude_handler=ignore_excludes), 0.0
        )

    def test_mi_pair_cols_and_entropies(self):
        """mi_pair calculates mi from a pair of columns and precalc entropies"""
        aln = ArrayAlignment(data={"1": "AB", "2": "AB"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1, h1=0.0, h2=0.0), 0.0)
        aln = ArrayAlignment(data={"1": "AB", "2": "BA"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1, h1=1.0, h2=1.0), 1.0)
        # incorrect positional entropies provided to ensure that the
        # precalculated values are used, and that entorpies are not
        # caluclated on-the-fly.
        aln = ArrayAlignment(data={"1": "AB", "2": "AB"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1, h1=1.0, h2=1.0), 2.0)

    def test_mi_pair_alt_calculator(self):
        """mi_pair uses alternate mi_calculator when provided"""
        aln = ArrayAlignment(data={"1": "AB", "2": "AB"}, moltype=PROTEIN)
        assert_allclose(mi_pair(aln, pos1=0, pos2=1), 0.0)
        assert_allclose(
            mi_pair(aln, pos1=0, pos2=1, mi_calculator=normalized_mi),
            DEFAULT_NULL_VALUE,
        )

    def test_mi_position_valid_input(self):
        """mi_position functions with varied valid input"""
        aln = ArrayAlignment(data={"1": "ACG", "2": "GAC"}, moltype=PROTEIN)
        assert_allclose(mi_position(aln, 0), array([1.0, 1.0, 1.0]))
        aln = ArrayAlignment(data={"1": "ACG", "2": "ACG"}, moltype=PROTEIN)
        assert_allclose(mi_position(aln, 0), array([0.0, 0.0, 0.0]))
        aln = ArrayAlignment(data={"1": "ACG", "2": "ACG"}, moltype=PROTEIN)
        assert_allclose(mi_position(aln, 2), array([0.0, 0.0, 0.0]))

    def test_mi_position_from_alignment_nmi(self):
        """mi_position functions w/ alternate mi_calculator"""
        aln = ArrayAlignment(data={"1": "ACG", "2": "ACG"}, moltype=PROTEIN)
        assert_allclose(mi_position(aln, 0), array([0.0, 0.0, 0.0]))
        aln = ArrayAlignment(data={"1": "ACG", "2": "ACG"}, moltype=PROTEIN)
        assert_allclose(
            mi_position(aln, 0, mi_calculator=normalized_mi),
            array([DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE]),
        )

    def test_mi_position_from_alignment_default_exclude_handling(self):
        """mi_position handles excludes by setting to null_value"""
        aln = ArrayAlignment(data={"1": "ACG", "2": "G-C"}, moltype=PROTEIN)
        assert_allclose(mi_position(aln, 0), array([1.0, DEFAULT_NULL_VALUE, 1.0]))
        aln = ArrayAlignment(data={"1": "ACG", "2": "GPC"}, moltype=PROTEIN)
        assert_allclose(
            mi_position(aln, 0, excludes="P"), array([1.0, DEFAULT_NULL_VALUE, 1.0])
        )

    def test_mi_position_from_alignment_non_default_exclude_handling(self):
        """mi_position handles excludes w/ non-default method"""
        aln = ArrayAlignment(data={"1": "ACG", "2": "G-C"}, moltype=PROTEIN)
        assert_allclose(
            mi_position(aln, 0, exclude_handler=ignore_excludes), array([1.0, 1.0, 1.0])
        )

    def test_mi_alignment_excludes(self):
        """mi_alignment handles excludes properly"""
        expected = array(
            [
                [0.0, DEFAULT_NULL_VALUE, 0.0],
                [DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE],
                [0.0, DEFAULT_NULL_VALUE, 0.0],
            ]
        )
        # gap in second column
        aln = ArrayAlignment(data={"1": "ACG", "2": "A-G"}, moltype=PROTEIN)
        assert_allclose(mi_alignment(aln), expected)

        # excludes = 'P'
        aln = ArrayAlignment(data={"1": "ACG", "2": "APG"}, moltype=PROTEIN)
        assert_allclose(mi_alignment(aln, excludes="P"), expected)

        # gap in first column
        expected = array(
            [
                [DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE],
                [DEFAULT_NULL_VALUE, 0.0, 0.0],
                [DEFAULT_NULL_VALUE, 0.0, 0.0],
            ]
        )
        aln = ArrayAlignment(data={"1": "-CG", "2": "ACG"}, moltype=PROTEIN)
        assert_allclose(mi_alignment(aln), expected)

    def test_mi_alignment_high(self):
        """mi_alignment detected perfectly correlated columns"""
        expected = [[1.0, 1.0], [1.0, 1.0]]
        aln = ArrayAlignment(data={"1": "AG", "2": "GA"}, moltype=PROTEIN)
        assert_allclose(mi_alignment(aln), expected)

    def test_mi_alignment_low(self):
        """mi_alignment detected in perfectly uncorrelated columns"""
        expected = [[0.0, 0.0], [0.0, 1.0]]
        aln = ArrayAlignment(data={"1": "AG", "2": "AC"}, moltype=PROTEIN)
        assert_allclose(mi_alignment(aln), expected)

    def test_resampled_mi_alignment(self):
        """resampled_mi_alignment returns without error"""
        aln = ArrayAlignment(
            data={"1": "ACDEF", "2": "ACFEF", "3": "ACGEF"}, moltype=PROTEIN
        )
        resampled_mi_alignment(aln)
        aln = ArrayAlignment(
            data={"1": "ACDEF", "2": "ACF-F", "3": "ACGEF"}, moltype=PROTEIN
        )
        resampled_mi_alignment(aln)

    def test_coevolve_alignment(self):
        """coevolve_alignment functions as expected with varied input"""
        aln1 = ArrayAlignment(
            data={"1": "ACDEF", "2": "ACFEF", "3": "ACGEF"}, moltype=PROTEIN
        )
        # no kwargs passed
        assert_allclose(coevolve_alignment(mi_alignment, aln1), mi_alignment(aln1))
        # different method passed
        assert_allclose(coevolve_alignment(nmi_alignment, aln1), nmi_alignment(aln1))
        # kwargs passed
        assert_allclose(
            coevolve_alignment(mi_alignment, aln1, mi_calculator=nmi),
            nmi_alignment(aln1),
        )

    def test_build_coevolution_matrix_filepath(self):
        """build_coevolution_matrix_filepath functions w/ varied input"""
        self.assertEqual(build_coevolution_matrix_filepath("./blah.fasta"), "./blah")
        self.assertEqual(build_coevolution_matrix_filepath("blah.fasta"), "./blah")
        self.assertEqual(build_coevolution_matrix_filepath("blah"), "./blah")
        self.assertEqual(build_coevolution_matrix_filepath("./blah"), "./blah")

        self.assertEqual(
            build_coevolution_matrix_filepath(
                "./blah.fasta", output_dir="./duh/", method="xx", alphabet="yyy"
            ),
            "./duh/blah.yyy.xx",
        )
        self.assertEqual(
            build_coevolution_matrix_filepath(
                "./blah.fasta",
                output_dir="./duh/",
                method="xx",
                alphabet="yyy",
                parameter=0.25,
            ),
            "./duh/blah.yyy.xx",
        )
        self.assertEqual(
            build_coevolution_matrix_filepath(
                "./blah.fasta", output_dir="./duh/", method="xx"
            ),
            "./duh/blah.xx",
        )
        self.assertEqual(
            build_coevolution_matrix_filepath(
                "./blah.fasta", output_dir="./duh/", method="sca", parameter=0.25
            ),
            "./duh/blah.sca_25",
        )
        self.assertEqual(
            build_coevolution_matrix_filepath(
                "./blah.fasta",
                output_dir="./duh/",
                method="sca",
                parameter=0.25,
                alphabet="xx",
            ),
            "./duh/blah.xx.sca_25",
        )
        # no trailing / to output_dir
        self.assertEqual(
            build_coevolution_matrix_filepath(
                "./blah.fasta",
                output_dir="./duh",
                method="sca",
                parameter=0.25,
                alphabet="xx",
            ),
            "./duh/blah.xx.sca_25",
        )

        self.assertRaises(
            ValueError,
            build_coevolution_matrix_filepath,
            "./blah.fasta",
            "./duh/",
            "sca",
        )
        self.assertRaises(
            ValueError,
            build_coevolution_matrix_filepath,
            "./blah.fasta",
            "./duh/",
            "sca",
            "xx",
        )

    def test_pickle_coevolution_result_error(self):
        """pickle matrix: IOError handled correctly"""
        m = array([[1, 2], [3, 4]])
        self.assertRaises(IOError, pickle_coevolution_result, m, "")

    def test_unpickle_coevolution_result_error(self):
        """unpickle matrix: IOError handled correctly"""
        self.assertRaises(IOError, unpickle_coevolution_result, "invalid/file/path.pkl")

    def test_pickle_and_unpickle(self):
        """unpickle(pickle(matrix)) == matrix"""
        for expected in [4.5, array([1.2, 4.3, 5.5]), array([[1.4, 2.2], [3.0, 0.4]])]:
            filepath = mktemp()
            pickle_coevolution_result(expected, filepath)
            actual = unpickle_coevolution_result(filepath)
            assert_allclose(actual, expected)
            remove(filepath)

    def test_csv_coevolution_result_error(self):
        """matrix -> csv: IOError handled correctly"""
        m = array([[1, 2], [3, 4]])
        self.assertRaises(IOError, coevolution_matrix_to_csv, m, "")

    def test_uncsv_coevolution_result_error(self):
        """csv -> matrix: IOError handled correctly"""
        self.assertRaises(IOError, csv_to_coevolution_matrix, "invalid/file/path.pkl")

    def test_csv_and_uncsv(self):
        """converting to/from csv matrix results in correct coevolution matrix"""
        expected = array([[1.4, 2.2], [DEFAULT_NULL_VALUE, 0.4]])
        filepath = mktemp()
        coevolution_matrix_to_csv(expected, filepath)
        actual = csv_to_coevolution_matrix(filepath)
        assert_allclose(actual, expected)
        remove(filepath)

    def test_parse_coevolution_matrix_filepath(self):
        """Parsing matrix filepaths works as expected."""
        expected = ("myosin_995", "a1_4", "nmi")
        self.assertEqual(
            parse_coevolution_matrix_filepath("pkls/myosin_995.a1_4.nmi.pkl"), expected
        )
        self.assertEqual(
            parse_coevolution_matrix_filepath("pkls/myosin_995.a1_4.nmi.csv"), expected
        )
        expected = ("p53", "orig", "mi")
        self.assertEqual(parse_coevolution_matrix_filepath("p53.orig.mi.pkl"), expected)
        self.assertEqual(parse_coevolution_matrix_filepath("p53.orig.mi.csv"), expected)

    def test_parse_coevolution_matrix_filepath_error(self):
        """Parsing matrix file paths handles invalid filepaths"""
        self.assertRaises(
            ValueError, parse_coevolution_matrix_filepath, "pkls/myosin_995.nmi.pkl"
        )
        self.assertRaises(
            ValueError, parse_coevolution_matrix_filepath, "pkls/myosin_995.pkl"
        )
        self.assertRaises(
            ValueError, parse_coevolution_matrix_filepath, "pkls/myosin_995"
        )
        self.assertRaises(ValueError, parse_coevolution_matrix_filepath, "")

    def test_identify_aln_positions_above_threshold(self):
        """Extracting scores above threshold works as expected"""
        m = array(
            [
                [
                    DEFAULT_NULL_VALUE,
                    DEFAULT_NULL_VALUE,
                    DEFAULT_NULL_VALUE,
                    DEFAULT_NULL_VALUE,
                ],
                [0.3, 1.0, DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE],
                [0.25, 0.75, 1.0, DEFAULT_NULL_VALUE],
                [0.9, 0.751, 0.8, 1.0],
            ]
        )
        self.assertEqual(identify_aln_positions_above_threshold(m, 0.75, 0), [])
        self.assertEqual(identify_aln_positions_above_threshold(m, 0.75, 1), [1])
        self.assertEqual(identify_aln_positions_above_threshold(m, 0.75, 2), [1, 2])
        self.assertEqual(
            identify_aln_positions_above_threshold(m, 0.75, 3), [0, 1, 2, 3]
        )

        m = ltm_to_symmetric(m)
        self.assertEqual(identify_aln_positions_above_threshold(m, 0.75, 0), [3])
        self.assertEqual(identify_aln_positions_above_threshold(m, 0.75, 1), [1, 2, 3])
        self.assertEqual(identify_aln_positions_above_threshold(m, 0.75, 2), [1, 2, 3])
        self.assertEqual(
            identify_aln_positions_above_threshold(m, 0.75, 3), [0, 1, 2, 3]
        )

        self.assertEqual(identify_aln_positions_above_threshold(m, 1.1, 0), [])
        self.assertEqual(identify_aln_positions_above_threshold(m, -5.0, 0), [1, 2, 3])
        self.assertEqual(
            identify_aln_positions_above_threshold(m, -5.0, 1), [0, 1, 2, 3]
        )

    def test_count_ge_threshold(self):
        """count_ge_threshold works as expected"""
        m = array([[DEFAULT_NULL_VALUE] * 3] * 3)
        self.assertEqual(count_ge_threshold(m, 1.0), (0, 0))
        self.assertEqual(
            count_ge_threshold(m, DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE), (0, 0)
        )
        self.assertEqual(count_ge_threshold(m, 1.0, 42), (0, 9))

        m = array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        self.assertEqual(count_ge_threshold(m, 4), (5, 9))
        self.assertEqual(count_ge_threshold(m, 8), (1, 9))
        self.assertEqual(count_ge_threshold(m, 9), (0, 9))

        m = array(
            [
                [0, DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE],
                [DEFAULT_NULL_VALUE, 4, 5],
                [6, 7, 8],
            ]
        )
        self.assertEqual(count_ge_threshold(m, 4), (5, 6))
        self.assertEqual(count_ge_threshold(m, 8), (1, 6))
        self.assertEqual(count_ge_threshold(m, 9), (0, 6))

    def test_count_le_threshold(self):
        """count_le_threshold works as expected"""
        m = array([[DEFAULT_NULL_VALUE] * 3] * 3)
        self.assertEqual(count_le_threshold(m, 1.0), (0, 0))
        self.assertEqual(
            count_le_threshold(m, DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE), (0, 0)
        )
        self.assertEqual(count_le_threshold(m, 1.0, 42), (0, 9))

        m = array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        self.assertEqual(count_le_threshold(m, 4), (5, 9))
        self.assertEqual(count_le_threshold(m, 8), (9, 9))
        self.assertEqual(count_le_threshold(m, 9), (9, 9))

        m = array(
            [
                [0, DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE],
                [DEFAULT_NULL_VALUE, 4, 5],
                [6, 7, 8],
            ]
        )
        self.assertEqual(count_le_threshold(m, 4), (2, 6))
        self.assertEqual(count_le_threshold(m, 8), (6, 6))
        self.assertEqual(count_le_threshold(m, 9), (6, 6))

    def test_count_ge_threshold_symmetric_ignore_diagonal(self):
        """count_ge_threshold works with symmetric and/or ignoring diag = True"""
        # no good scores, varied null value
        m = array([[DEFAULT_NULL_VALUE] * 3] * 3)
        self.assertEqual(count_ge_threshold(m, 1.0, symmetric=True), (0, 0))
        self.assertEqual(count_ge_threshold(m, 1.0, symmetric=True), (0, 0))
        self.assertEqual(count_ge_threshold(m, 1.0, 42, symmetric=True), (0, 6))
        self.assertEqual(count_ge_threshold(m, 1.0, ignore_diagonal=True), (0, 0))
        self.assertEqual(count_ge_threshold(m, 1.0, ignore_diagonal=True), (0, 0))
        self.assertEqual(count_ge_threshold(m, 1.0, 42, ignore_diagonal=True), (0, 6))
        self.assertEqual(
            count_ge_threshold(m, 1.0, ignore_diagonal=True, symmetric=True), (0, 0)
        )
        self.assertEqual(
            count_ge_threshold(m, 1.0, ignore_diagonal=True, symmetric=True), (0, 0)
        )
        self.assertEqual(
            count_ge_threshold(m, 1.0, 42, ignore_diagonal=True, symmetric=True), (0, 3)
        )

        #  no null values, varied other values
        m = array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        self.assertEqual(count_ge_threshold(m, 4), (5, 9))
        self.assertEqual(count_ge_threshold(m, 4, symmetric=True), (4, 6))
        self.assertEqual(count_ge_threshold(m, 4, ignore_diagonal=True), (3, 6))
        self.assertEqual(
            count_ge_threshold(m, 4, symmetric=True, ignore_diagonal=True), (2, 3)
        )

        # null and mixed values
        m = array(
            [
                [0, DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE],
                [3, 4, DEFAULT_NULL_VALUE],
                [DEFAULT_NULL_VALUE, 7, 8],
            ]
        )
        self.assertEqual(count_ge_threshold(m, 4), (3, 5))
        self.assertEqual(count_ge_threshold(m, 4, symmetric=True), (3, 5))
        self.assertEqual(count_ge_threshold(m, 4, ignore_diagonal=True), (1, 2))
        self.assertEqual(
            count_ge_threshold(m, 4, symmetric=True, ignore_diagonal=True), (1, 2)
        )

    def test_count_le_threshold_symmetric_ignore_diagonal(self):
        """count_le_threshold works with symmetric and/or ignoring diag = True"""
        # varied null value
        m = array([[DEFAULT_NULL_VALUE] * 3] * 3)
        self.assertEqual(count_le_threshold(m, 1.0, symmetric=True), (0, 0))
        self.assertEqual(count_le_threshold(m, 1.0, symmetric=True), (0, 0))
        self.assertEqual(count_le_threshold(m, 1.0, 42, symmetric=True), (0, 6))
        self.assertEqual(count_le_threshold(m, 1.0, ignore_diagonal=True), (0, 0))
        self.assertEqual(count_le_threshold(m, 1.0, ignore_diagonal=True), (0, 0))
        self.assertEqual(count_le_threshold(m, 1.0, 42, ignore_diagonal=True), (0, 6))
        self.assertEqual(
            count_le_threshold(m, 1.0, ignore_diagonal=True, symmetric=True), (0, 0)
        )
        self.assertEqual(
            count_le_threshold(m, 1.0, ignore_diagonal=True, symmetric=True), (0, 0)
        )
        self.assertEqual(
            count_le_threshold(m, 1.0, 42, ignore_diagonal=True, symmetric=True), (0, 3)
        )

        #  no null values, varied other values
        m = array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        self.assertEqual(count_le_threshold(m, 4), (5, 9))
        self.assertEqual(count_le_threshold(m, 4, symmetric=True), (3, 6))
        self.assertEqual(count_le_threshold(m, 4, ignore_diagonal=True), (3, 6))
        self.assertEqual(
            count_le_threshold(m, 4, symmetric=True, ignore_diagonal=True), (1, 3)
        )

        # null and mixed values
        m = array(
            [
                [0, DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE],
                [3, 4, DEFAULT_NULL_VALUE],
                [DEFAULT_NULL_VALUE, 7, 8],
            ]
        )
        self.assertEqual(count_le_threshold(m, 4), (3, 5))
        self.assertEqual(count_le_threshold(m, 4, symmetric=True), (3, 5))
        self.assertEqual(count_le_threshold(m, 4, ignore_diagonal=True), (1, 2))
        self.assertEqual(
            count_le_threshold(m, 4, symmetric=True, ignore_diagonal=True), (1, 2)
        )

    def test_aln_position_pairs_cmp_threshold_intramolecular(self):
        """aln_position_pairs_ge_threshold: intramolecular matrix"""
        m = array(
            [
                [0, DEFAULT_NULL_VALUE, DEFAULT_NULL_VALUE],
                [3, 4, DEFAULT_NULL_VALUE],
                [DEFAULT_NULL_VALUE, 7, 8],
            ]
        )
        # cmp_function = ge
        self.assertEqual(
            aln_position_pairs_cmp_threshold(m, 3.5, greater_equal),
            [(1, 1), (2, 1), (2, 2)],
        )
        # cmp_function = greater_equal, alt null_value
        self.assertEqual(
            aln_position_pairs_cmp_threshold(m, 3.5, greater_equal, null_value=4),
            [(2, 1), (2, 2)],
        )
        # cmp_function = le
        self.assertEqual(
            aln_position_pairs_cmp_threshold(m, 3.5, less_equal), [(0, 0), (1, 0)]
        )

        # results equal results with wrapper functions
        self.assertEqual(
            aln_position_pairs_cmp_threshold(m, 3.5, greater_equal),
            aln_position_pairs_ge_threshold(m, 3.5),
        )
        self.assertEqual(
            aln_position_pairs_cmp_threshold(m, 3.5, less_equal),
            aln_position_pairs_le_threshold(m, 3.5),
        )
        self.assertEqual(
            aln_position_pairs_cmp_threshold(m, 3.5, greater_equal, null_value=4),
            aln_position_pairs_ge_threshold(m, 3.5, null_value=4),
        )
        self.assertEqual(
            aln_position_pairs_cmp_threshold(m, 3.5, less_equal, null_value=0),
            aln_position_pairs_le_threshold(m, 3.5, null_value=0),
        )

    def test_aln_position_pairs_ge_threshold_intermolecular(self):
        """aln_position_pairs_ge_threshold: intermolecular matrix"""
        m = array([[1.0, 10.0, 4.0, 3.0], [9.0, 18.0, 5.0, 6.0]])
        # error if failed to specify intermolecular_data_only=True
        self.assertRaises(
            AssertionError, aln_position_pairs_cmp_threshold, m, 3.5, greater_equal
        )
        # cmp_function = ge
        self.assertEqual(
            aln_position_pairs_cmp_threshold(
                m, 3.5, greater_equal, intermolecular_data_only=True
            ),
            [(1, 4), (2, 4), (0, 5), (1, 5), (2, 5), (3, 5)],
        )
        # cmp_function = greater_equal, alt null_value
        self.assertEqual(
            aln_position_pairs_cmp_threshold(
                m, 3.5, greater_equal, null_value=18.0, intermolecular_data_only=True
            ),
            [(1, 4), (2, 4), (0, 5), (2, 5), (3, 5)],
        )
        # cmp_function = le
        self.assertEqual(
            aln_position_pairs_cmp_threshold(
                m, 3.5, less_equal, intermolecular_data_only=True
            ),
            [(0, 4), (3, 4)],
        )

        # results equal results with wrapper functions
        self.assertEqual(
            aln_position_pairs_cmp_threshold(
                m, 3.5, greater_equal, intermolecular_data_only=True
            ),
            aln_position_pairs_ge_threshold(m, 3.5, intermolecular_data_only=True),
        )
        self.assertEqual(
            aln_position_pairs_cmp_threshold(
                m, 3.5, less_equal, intermolecular_data_only=True
            ),
            aln_position_pairs_le_threshold(m, 3.5, intermolecular_data_only=True),
        )

        self.assertEqual(
            aln_position_pairs_cmp_threshold(
                m, 3.5, greater_equal, null_value=4.0, intermolecular_data_only=True
            ),
            aln_position_pairs_ge_threshold(
                m, 3.5, null_value=4.0, intermolecular_data_only=True
            ),
        )
        self.assertEqual(
            aln_position_pairs_cmp_threshold(
                m, 3.5, less_equal, null_value=18.0, intermolecular_data_only=True
            ),
            aln_position_pairs_le_threshold(
                m, 3.5, null_value=18.0, intermolecular_data_only=True
            ),
        )

    def test_is_parsimony_informative_strict(self):
        """is_parsimony_informative functions as expected with strict=True"""
        freqs = {"A": 25}
        self.assertFalse(is_parsimony_informative(freqs, strict=True))
        freqs = {"A": 25, "-": 25}
        self.assertFalse(is_parsimony_informative(freqs, strict=True))
        freqs = {"A": 25, "?": 25}
        self.assertFalse(is_parsimony_informative(freqs, strict=True))
        freqs = {"A": 25, "B": 1}
        self.assertFalse(is_parsimony_informative(freqs, strict=True))
        freqs = {"A": 1, "B": 1, "C": 1, "D": 1, "E": 1}
        self.assertFalse(is_parsimony_informative(freqs, strict=True))
        freqs = {"A": 2, "B": 1, "C": 1, "D": 1, "E": 1}
        self.assertFalse(is_parsimony_informative(freqs, strict=True))
        freqs = {"A": 2, "B": 2, "C": 1, "D": 1, "E": 1}
        self.assertFalse(is_parsimony_informative(freqs, strict=True))

        freqs = {"A": 25, "B": 2}
        self.assertTrue(is_parsimony_informative(freqs, strict=True))
        freqs = {"A": 2, "B": 2, "C": 2, "D": 2, "E": 2}
        self.assertTrue(is_parsimony_informative(freqs, strict=True))

    def test_is_parsimony_informative_non_strict(self):
        """is_parsimony_informative functions as expected with strict=False"""
        freqs = {"A": 25}
        self.assertFalse(is_parsimony_informative(freqs, strict=False))
        freqs = {"A": 25, "-": 25}
        self.assertFalse(is_parsimony_informative(freqs, strict=False))
        freqs = {"A": 25, "?": 25}
        self.assertFalse(is_parsimony_informative(freqs, strict=False))
        freqs = {"A": 25, "B": 1}
        self.assertFalse(is_parsimony_informative(freqs, strict=False))
        freqs = {"A": 1, "B": 1, "C": 1, "D": 1, "E": 1}
        self.assertFalse(is_parsimony_informative(freqs, strict=False))
        freqs = {"A": 2, "B": 1, "C": 1, "D": 1, "E": 1}
        self.assertFalse(is_parsimony_informative(freqs, strict=False))

        freqs = {"A": 2, "B": 2, "C": 1, "D": 1, "E": 1}
        self.assertTrue(is_parsimony_informative(freqs, strict=False))
        freqs = {"A": 25, "B": 2}
        self.assertTrue(is_parsimony_informative(freqs, strict=False))
        freqs = {"A": 2, "B": 2, "C": 2, "D": 2, "E": 2}
        self.assertTrue(is_parsimony_informative(freqs, strict=False))

    def test_is_parsimony_informative_non_default(self):
        """is_parsimony_informative functions w non default paramters"""
        # NEED TO UPDATE THESE TESTS BASED ON MY ERROR IN THE
        # DEFINITION OF PARSIMONY INFORMATIVE.
        # changed minimum_count
        freqs = {"A": 25, "B": 2}
        self.assertFalse(is_parsimony_informative(freqs, minimum_count=3, strict=False))
        freqs = {"A": 25, "B": 1}
        self.assertTrue(is_parsimony_informative(freqs, minimum_count=1, strict=False))

        # different value of strict yields different results
        freqs = {"A": 25, "B": 2, "C": 3}
        self.assertTrue(is_parsimony_informative(freqs, minimum_count=3, strict=False))
        self.assertFalse(is_parsimony_informative(freqs, minimum_count=3, strict=True))

        # changed minimum_differences
        freqs = {"A": 25, "B": 25}
        self.assertFalse(
            is_parsimony_informative(freqs, minimum_differences=3, strict=False)
        )
        freqs = {"A": 25}
        self.assertTrue(
            is_parsimony_informative(freqs, minimum_differences=1, strict=False)
        )

        # changed ignored
        freqs = {"A": 25, "-": 25, "?": 25}
        self.assertTrue(is_parsimony_informative(freqs, ignored=None, strict=False))
        freqs = {"A": 25, "?": 25}
        self.assertTrue(is_parsimony_informative(freqs, ignored="", strict=False))
        freqs = {"A": 25, "-": 25}
        self.assertTrue(is_parsimony_informative(freqs, ignored=None, strict=False))
        freqs = {"A": 25, "C": 25}
        self.assertFalse(is_parsimony_informative(freqs, ignored="A", strict=False))

    def test_filter_non_parsimony_informative_intramolecular(self):
        """non-parsimony informative sites in intramolecular matrix -> null"""
        aln = make_aligned_seqs(
            data={"1": "ACDE", "2": "ACDE", "3": "ACDE", "4": "ACDE"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        expected = array([[DEFAULT_NULL_VALUE] * 4] * 4)
        filter_non_parsimony_informative(aln, m)
        assert_allclose(m, expected)

        aln = make_aligned_seqs(
            data={"1": "ACDE", "2": "FCDE", "3": "ACDE", "4": "FCDE"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array(
            [
                [42.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        expected = array([[DEFAULT_NULL_VALUE] * 4] * 4)
        expected[0, 0] = 42.0
        filter_non_parsimony_informative(aln, m)
        assert_allclose(m, expected)

    def test_filter_non_parsimony_informative_intermolecular(self):
        """non-parsimony informative sites in intermolecular matrix -> null"""
        # all non-parsimony informative
        aln = make_aligned_seqs(
            data={"1": "ACDEWQ", "2": "ACDEWQ", "3": "ACDEWQ", "4": "ACDEWQ"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array([[1.0, 10.0, 4.0, 3.0], [9.0, 18.0, 5.0, 6.0]])
        expected = array([[DEFAULT_NULL_VALUE] * 4] * 2)
        filter_non_parsimony_informative(aln, m, intermolecular_data_only=True)
        assert_allclose(m, expected)
        # one non-parsimony informative pair of positions
        aln = make_aligned_seqs(
            data={"1": "FCDEWD", "2": "ACDEWQ", "3": "ACDEWD", "4": "FCDEWQ"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array([[1.0, 10.0, 4.0, 3.0], [9.0, 18.0, 5.0, 6.0]])
        expected = array([[DEFAULT_NULL_VALUE] * 4] * 2)
        expected[1, 0] = 9.0
        filter_non_parsimony_informative(aln, m, intermolecular_data_only=True)
        assert_allclose(m, expected)
        # all parsimony informative
        aln = make_aligned_seqs(
            data={"1": "FFFFFF", "2": "FFFFFF", "3": "GGGGGG", "4": "GGGGGG"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array([[1.0, 10.0, 4.0, 3.0], [9.0, 18.0, 5.0, 6.0]])
        expected = array([[1.0, 10.0, 4.0, 3.0], [9.0, 18.0, 5.0, 6.0]])
        filter_non_parsimony_informative(aln, m, intermolecular_data_only=True)
        assert_allclose(m, expected)

    def test_filter_exclude_positions_intramolecular(self):
        """filter_exclude_positions: functions for intramolecular data"""
        # filter zero positions (no excludes)
        aln = make_aligned_seqs(
            data={"1": "WCDE", "2": "ACDE", "3": "ACDE", "4": "ACDE"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        expected = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        filter_exclude_positions(aln, m)
        assert_allclose(m, expected)
        # filter zero positions (max_exclude_percentage = percent exclude)
        aln = make_aligned_seqs(
            data={"1": "-CDE", "2": "A-DE", "3": "AC-E", "4": "ACD-"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        expected = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        filter_exclude_positions(aln, m, max_exclude_percent=0.25)
        assert_allclose(m, expected)
        # filter zero positions (max_exclude_percentage too high)
        aln = make_aligned_seqs(
            data={"1": "-CDE", "2": "A-DE", "3": "AC-E", "4": "ACD-"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        expected = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        filter_exclude_positions(aln, m, max_exclude_percent=0.5)
        assert_allclose(m, expected)
        # filter one position (defualt max_exclude_percentage)
        aln = make_aligned_seqs(
            data={"1": "-CDE", "2": "ACDE", "3": "ACDE", "4": "ACDE"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        expected = array(
            [
                [DEFAULT_NULL_VALUE] * 4,
                [DEFAULT_NULL_VALUE, 18.0, 5.0, 6.0],
                [DEFAULT_NULL_VALUE, 1.0, 3.0, 2.0],
                [DEFAULT_NULL_VALUE, 0.0, 1.0, 33.0],
            ]
        )
        filter_exclude_positions(aln, m)
        assert_allclose(m, expected)
        # filter one position (non-defualt max_exclude_percentage)
        aln = make_aligned_seqs(
            data={"1": "-CDE", "2": "ACDE", "3": "ACDE", "4": "-CDE"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        expected = array(
            [
                [DEFAULT_NULL_VALUE] * 4,
                [DEFAULT_NULL_VALUE, 18.0, 5.0, 6.0],
                [DEFAULT_NULL_VALUE, 1.0, 3.0, 2.0],
                [DEFAULT_NULL_VALUE, 0.0, 1.0, 33.0],
            ]
        )
        filter_exclude_positions(aln, m, max_exclude_percent=0.49)
        assert_allclose(m, expected)
        # filter all positions (defualt max_exclude_percentage)
        aln = make_aligned_seqs(
            data={"1": "----", "2": "ACDE", "3": "ACDE", "4": "ACDE"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        expected = array([[DEFAULT_NULL_VALUE] * 4] * 4)
        filter_exclude_positions(aln, m)
        assert_allclose(m, expected)
        # filter all positions (non-defualt max_exclude_percentage)
        aln = make_aligned_seqs(
            data={"1": "----", "2": "A-DE", "3": "AC--", "4": "-CDE"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 3.0],
            ]
        )
        expected = array([[DEFAULT_NULL_VALUE] * 4] * 4)
        filter_exclude_positions(aln, m, max_exclude_percent=0.49)
        assert_allclose(m, expected)

        # filter one position (defualt max_exclude_percentage,
        # non-defualt excludes)
        aln = make_aligned_seqs(
            data={"1": "WCDE", "2": "ACDE", "3": "ACDE", "4": "ACDE"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        expected = array(
            [
                [DEFAULT_NULL_VALUE] * 4,
                [DEFAULT_NULL_VALUE, 18.0, 5.0, 6.0],
                [DEFAULT_NULL_VALUE, 1.0, 3.0, 2.0],
                [DEFAULT_NULL_VALUE, 0.0, 1.0, 33.0],
            ]
        )
        filter_exclude_positions(aln, m, excludes="W")
        assert_allclose(m, expected)

        # filter one position (defualt max_exclude_percentage,
        # non-defualt null_value)
        aln = make_aligned_seqs(
            data={"1": "-CDE", "2": "ACDE", "3": "ACDE", "4": "ACDE"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array(
            [
                [1.0, 10.0, 4.0, 3.0],
                [9.0, 18.0, 5.0, 6.0],
                [4.0, 1.0, 3.0, 2.0],
                [21.0, 0.0, 1.0, 33.0],
            ]
        )
        expected = array(
            [
                [999.0] * 4,
                [999.0, 18.0, 5.0, 6.0],
                [999.0, 1.0, 3.0, 2.0],
                [999.0, 0.0, 1.0, 33.0],
            ]
        )
        filter_exclude_positions(aln, m, null_value=999.0)
        assert_allclose(m, expected)

    def test_filter_exclude_positions_intermolecular(self):
        """filter_exclude_positions: functions for intermolecular data"""
        # these tests correspond to alignments of length 4 and 2 positions
        # respectively, hence a coevolution_matrix with shape = (2,4)

        # filter zero positions (no excludes)
        merged_aln = make_aligned_seqs(
            data={"1": "WCDEDE", "2": "ACDEDE", "3": "ACDEDE", "4": "ACDEDE"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array([[1.0, 10.0, 4.0, 3.0], [9.0, 18.0, 5.0, 6.0]])
        expected = array([[1.0, 10.0, 4.0, 3.0], [9.0, 18.0, 5.0, 6.0]])
        filter_exclude_positions(merged_aln, m, intermolecular_data_only=True)
        assert_allclose(m, expected)

        # filter one position (aln1)
        merged_aln = make_aligned_seqs(
            data={"1": "WC-EDE", "2": "ACDEDE", "3": "ACDEDE", "4": "ACDEDE"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array([[1.0, 10.0, 4.0, 3.0], [9.0, 18.0, 5.0, 6.0]])
        expected = array(
            [[1.0, 10.0, DEFAULT_NULL_VALUE, 3.0], [9.0, 18.0, DEFAULT_NULL_VALUE, 6.0]]
        )
        filter_exclude_positions(merged_aln, m, intermolecular_data_only=True)
        assert_allclose(m, expected)
        # filter one position (aln2)
        merged_aln = make_aligned_seqs(
            data={"1": "WCEEDE", "2": "ACDEDE", "3": "ACDEDE", "4": "ACDED-"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array([[1.0, 10.0, 4.0, 3.0], [9.0, 18.0, 5.0, 6.0]])
        expected = array([[1.0, 10.0, 4.0, 3.0], [DEFAULT_NULL_VALUE] * 4])
        filter_exclude_positions(merged_aln, m, intermolecular_data_only=True)
        assert_allclose(m, expected)

        # filter two positions (aln1 & aln2)
        merged_aln = make_aligned_seqs(
            data={"1": "-CEEDE", "2": "ACDEDE", "3": "ACDEDE", "4": "ACDED-"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array([[1.0, 10.0, 4.0, 3.0], [9.0, 18.0, 5.0, 6.0]])
        expected = array(
            [[DEFAULT_NULL_VALUE, 10.0, 4.0, 3.0], [DEFAULT_NULL_VALUE] * 4]
        )
        filter_exclude_positions(merged_aln, m, intermolecular_data_only=True)
        assert_allclose(m, expected)

        # filter two positions (aln1 & aln2, alt excludes)
        merged_aln = make_aligned_seqs(
            data={"1": "WCEEDE", "2": "ACDEDE", "3": "ACDEDE", "4": "ACDEDW"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array([[1.0, 10.0, 4.0, 3.0], [9.0, 18.0, 5.0, 6.0]])
        expected = array(
            [[DEFAULT_NULL_VALUE, 10.0, 4.0, 3.0], [DEFAULT_NULL_VALUE] * 4]
        )
        filter_exclude_positions(
            merged_aln, m, intermolecular_data_only=True, excludes="W"
        )
        assert_allclose(m, expected)

        # filter two positions (aln1 & aln2, alt null_value)
        merged_aln = make_aligned_seqs(
            data={"1": "-CEEDE", "2": "ACDEDE", "3": "ACDEDE", "4": "ACDED-"},
            moltype=PROTEIN,
            array_align=True,
        )
        m = array([[1.0, 10.0, 4.0, 3.0], [9.0, 18.0, 5.0, 6.0]])
        expected = array([[999.0, 10.0, 4.0, 3.0], [999.0] * 4])
        filter_exclude_positions(
            merged_aln, m, intermolecular_data_only=True, null_value=999.0
        )
        assert_allclose(m, expected)

    def test_filter_threshold_based_multiple_interdependency_intermolecular(self):
        "multiple interdependency filter functions with intermolecular data"
        ## cmp_function = ge
        # lower boundary
        null = DEFAULT_NULL_VALUE
        m = array(
            [
                [0.63, 0.00, null],
                [0.75, 0.10, 0.45],
                [0.95, 0.32, 0.33],
                [1.00, 0.95, 0.11],
            ]
        )
        expected = array(
            [
                [null, null, null],
                [null, null, 0.45],
                [null, null, null],
                [null, null, null],
            ]
        )
        actual = filter_threshold_based_multiple_interdependency(
            None, m, 0.95, 0, greater_equal, True
        )
        assert_allclose(actual, expected)
        # realisitic test case
        m = array(
            [
                [0.63, 0.00, null],
                [0.75, 0.10, 0.45],
                [0.95, 0.32, 0.33],
                [1.00, 0.95, 0.11],
            ]
        )
        expected = array(
            [
                [null, 0.00, null],
                [null, 0.10, 0.45],
                [null, 0.32, 0.33],
                [null, null, null],
            ]
        )
        actual = filter_threshold_based_multiple_interdependency(
            None, m, 0.95, 1, greater_equal, True
        )
        assert_allclose(actual, expected)
        # upper boundary, nothing filtered
        null = DEFAULT_NULL_VALUE
        m = array(
            [
                [0.63, 0.00, null],
                [0.75, 0.10, 0.45],
                [0.95, 0.32, 0.33],
                [1.00, 0.95, 0.11],
            ]
        )
        expected = m
        actual = filter_threshold_based_multiple_interdependency(
            None, m, 0.95, 5, greater_equal, True
        )
        assert_allclose(actual, expected)

        # cmp_function = less_equal, realistic test case
        m = array(
            [
                [0.63, 0.00, null],
                [0.75, 0.10, 0.45],
                [0.95, 0.32, 0.33],
                [1.00, 0.95, 0.11],
            ]
        )
        expected = array(
            [
                [0.63, null, null],
                [0.75, null, null],
                [null, null, null],
                [1.00, null, null],
            ]
        )
        actual = filter_threshold_based_multiple_interdependency(
            None, m, 0.35, 1, less_equal, True
        )
        assert_allclose(actual, expected)

    def test_filter_threshold_based_multiple_interdependency_intramolecular(self):
        "multiple interdependency filter functions with intramolecular data"
        null = DEFAULT_NULL_VALUE
        ## cmp_function = ge
        # lower bound, everything filtered
        m = array(
            [
                [0.63, 0.75, 0.95, 1.00],
                [0.75, 0.10, null, 0.95],
                [0.95, null, 0.33, 0.11],
                [1.00, 0.95, 0.11, 1.00],
            ]
        )
        expected = array(
            [
                [null, null, null, null],
                [null, null, null, null],
                [null, null, null, null],
                [null, null, null, null],
            ]
        )
        actual = filter_threshold_based_multiple_interdependency(
            None, m, 0.95, 0, greater_equal
        )
        assert_allclose(actual, expected)

        # realistic test case
        m = array(
            [
                [0.63, 0.75, 0.95, 1.00],
                [0.75, 0.10, null, 0.95],
                [0.95, null, 0.33, 0.11],
                [1.00, 0.95, 0.11, 1.00],
            ]
        )
        expected = array(
            [
                [null, null, null, null],
                [null, 0.10, null, null],
                [null, null, 0.33, null],
                [null, null, null, null],
            ]
        )
        actual = filter_threshold_based_multiple_interdependency(
            None, m, 0.95, 1, greater_equal
        )
        assert_allclose(actual, expected)

        # upper boundary, nothing filtered
        m = array(
            [
                [0.63, 0.75, 0.95, 1.00],
                [0.75, 0.10, null, 0.95],
                [0.95, null, 0.33, 0.11],
                [1.00, 0.95, 0.11, 1.00],
            ]
        )
        expected = m
        actual = filter_threshold_based_multiple_interdependency(
            None, m, 0.95, 5, greater_equal
        )
        assert_allclose(actual, expected)

        ## cmp_function = le
        # realistic test case
        m = array(
            [
                [0.63, 0.75, 0.95, 1.00],
                [0.75, 0.10, null, 0.95],
                [0.95, null, 0.33, 0.11],
                [1.00, 0.95, 0.11, 1.00],
            ]
        )
        expected = array(
            [
                [0.63, 0.75, null, 1.00],
                [0.75, 0.10, null, 0.95],
                [null, null, null, null],
                [1.00, 0.95, null, 1.00],
            ]
        )
        actual = filter_threshold_based_multiple_interdependency(
            None, m, 0.33, 1, less_equal
        )
        assert_allclose(actual, expected)

    def test_probs_from_dict(self):
        """probs_from_dict: dict of probs -> list of probs in alphabet's order"""
        d = {"A": 0.25, "D": 0.52, "C": 0.23}
        a = list("ACD")
        assert_allclose(probs_from_dict(d, a), [0.25, 0.23, 0.52])
        a = list("ADC")
        assert_allclose(probs_from_dict(d, a), [0.25, 0.52, 0.23])
        a = list("DCA")
        assert_allclose(probs_from_dict(d, a), [0.52, 0.23, 0.25])
        a = CharAlphabet("DCA")
        assert_allclose(probs_from_dict(d, a), [0.52, 0.23, 0.25])

        # protein natural probs
        l = probs_from_dict(protein_dict, AAGapless)
        for i in range(20):
            assert_allclose(l[i], protein_dict[AAGapless[i]], 0.001)

    def test_freqs_from_aln(self):
        """freqs_from_aln: freqs of alphabet chars in aln is calc'ed correctly"""
        # non-default scaled_aln_size
        aln = ArrayAlignment(
            data=list(zip(list(range(4)), ["ACGT", "AGCT", "ACCC", "TAGG"])),
            moltype=PROTEIN,
        )
        alphabet = "ACGT"
        expected = [4, 5, 4, 3]
        assert_equal(freqs_from_aln(aln, alphabet, 16), expected)
        # change the order of the alphabet
        alphabet = "TGCA"
        expected = [3, 4, 5, 4]
        assert_equal(freqs_from_aln(aln, alphabet, 16), expected)
        # default scaled_aln_size, sums of freqs == 100
        alphabet = "ACGT"
        expected = [25.0, 31.25, 25, 18.75]
        assert_allclose(freqs_from_aln(aln, alphabet), expected)
        # alphabet char which doesn't show up gets zero freq
        alphabet = "ACGTW"
        expected = [25.0, 31.25, 25, 18.75, 0]
        assert_allclose(freqs_from_aln(aln, alphabet), expected)
        # alignment char which doesn't show up is silently ignored
        aln = ArrayAlignment(
            data=list(zip(list(range(4)), ["ACGT", "AGCT", "ACCC", "TWGG"])),
            moltype=PROTEIN,
        )
        alphabet = "ACGT"
        expected = [18.75, 31.25, 25, 18.75]
        assert_allclose(freqs_from_aln(aln, alphabet), expected)

    def test_freqs_to_array(self):
        """freqs_to_array: should convert CategoryCounter object to array"""
        # should work with empty object
        f = CategoryCounter()
        f2a = freqs_to_array
        assert_allclose(f2a(f, AAGapless), zeros(20))
        # should work with full object, omitting unwanted keys
        f = CategoryCounter({"A": 20, "Q": 30, "X": 20})
        expected = zeros(20)
        expected[AAGapless.index("A")] = 20
        expected[AAGapless.index("Q")] = 30
        assert_allclose(f2a(f, AAGapless), expected)

        # should work for normal dict and any alphabet
        d = {"A": 3, "D": 1, "C": 5, "E": 2}
        alpha = "ABCD"
        exp = array([3, 0, 5, 1])
        assert_allclose(f2a(d, alpha), exp)

    def test_get_allowed_perturbations(self):
        """get_allowed_perturbations: should work for different cutoff values"""
        counts = [50, 40, 10, 0]
        a = list("ACGT")
        self.assertEqual(get_allowed_perturbations(counts, 1.0, a), [])
        self.assertEqual(get_allowed_perturbations(counts, 0.51, a), [])
        self.assertEqual(get_allowed_perturbations(counts, 0.5, a), ["A"])
        self.assertEqual(get_allowed_perturbations(counts, 0.49, a), ["A"])
        self.assertEqual(get_allowed_perturbations(counts, 0.401, a), ["A"])
        self.assertEqual(get_allowed_perturbations(counts, 0.40, a), ["A", "C"])
        self.assertEqual(get_allowed_perturbations(counts, 0.399, a), ["A", "C"])
        self.assertEqual(get_allowed_perturbations(counts, 0.10, a), ["A", "C", "G"])
        self.assertEqual(get_allowed_perturbations(counts, 0.0, a), a)

    def test_get_subalignments(self):
        """get_subalignments: works with different alignment sizes and cutoffs"""
        aln = ArrayAlignment(
            data={1: "AAAA", 2: "AAAC", 3: "AACG", 4: "ACCT", 5: "ACG-"},
            moltype=PROTEIN,
        )
        sub_aln_0A = ArrayAlignment(
            data={1: "AAAA", 2: "AAAC", 3: "AACG", 4: "ACCT", 5: "ACG-"},
            moltype=PROTEIN,
        )
        sub_aln_0C = {}
        sub_aln_1A = ArrayAlignment(
            data={1: "AAAA", 2: "AAAC", 3: "AACG"}, moltype=PROTEIN
        )
        sub_aln_1C = ArrayAlignment(data={4: "ACCT", 5: "ACG-"}, moltype=PROTEIN)
        sub_aln_2G = ArrayAlignment(data={5: "ACG-"}, moltype=PROTEIN)

        self.assertEqual(get_subalignments(aln, 0, ["A"]), [sub_aln_0A])
        self.assertEqual(get_subalignments(aln, 0, ["C"]), [sub_aln_0C])
        self.assertEqual(get_subalignments(aln, 1, ["A"]), [sub_aln_1A])
        self.assertEqual(get_subalignments(aln, 1, ["C"]), [sub_aln_1C])
        self.assertEqual(
            get_subalignments(aln, 1, ["A", "C"]), [sub_aln_1A, sub_aln_1C]
        )
        self.assertEqual(get_subalignments(aln, 2, ["G"]), [sub_aln_2G])
        self.assertEqual(get_subalignments(aln, 3, ["-"]), [sub_aln_2G])

    def test_get_positional_frequencies_w_scale(self):
        """get_positional_frequencies: works with default scaled_aln_size"""
        aln = ArrayAlignment(
            data={1: "ACDE", 2: "ADDE", 3: "AEED", 4: "AFEF"}, moltype=PROTEIN
        )
        expected_0 = array([100.0, 0.0, 0.0, 0.0, 0.0])
        expected_1 = array([0.0, 25.0, 25.0, 25.0, 25.0])
        expected_2 = array([0.0, 0.0, 50.0, 50.0, 0.0])
        expected_3 = array([0.0, 0.0, 25.0, 50.0, 25.0])
        assert_allclose(get_positional_frequencies(aln, 0, "ACDEF"), expected_0)
        assert_allclose(get_positional_frequencies(aln, 1, "ACDEF"), expected_1)
        assert_allclose(get_positional_frequencies(aln, 2, "ACDEF"), expected_2)
        assert_allclose(get_positional_frequencies(aln, 3, "ACDEF"), expected_3)
        # extra characters (W) are silently ignored -- is this the desired
        # behavior?
        aln = ArrayAlignment(
            data={1: "WCDE", 2: "ADDE", 3: "AEED", 4: "AFEF"}, moltype=PROTEIN
        )
        expected_0 = array([75.0, 0.0, 0.0, 0.0, 0.0])
        assert_allclose(get_positional_frequencies(aln, 0, "ACDEF"), expected_0)
        # 20 residue amino acid alphabet
        aln = ArrayAlignment(
            data={1: "ACDE", 2: "ADDE", 3: "AEED", 4: "AFEF"}, moltype=PROTEIN
        )
        expected = array([100.0] + [0.0] * 19)
        assert_allclose(get_positional_frequencies(aln, 0, AAGapless), expected)

    def test_get_positional_frequencies(self):
        """get_positional_frequencies: works with non-default scaled_aln_size"""
        aln = ArrayAlignment(
            data={1: "ACDE", 2: "ADDE", 3: "AEED", 4: "AFEF"}, moltype=PROTEIN
        )
        expected_0 = array([4.0, 0.0, 0.0, 0.0, 0.0])
        expected_1 = array([0.0, 1.0, 1.0, 1.0, 1.0])
        expected_2 = array([0.0, 0.0, 2.0, 2.0, 0.0])
        expected_3 = array([0.0, 0.0, 1.0, 2.0, 1.0])
        assert_allclose(get_positional_frequencies(aln, 0, "ACDEF", 4), expected_0)
        assert_allclose(get_positional_frequencies(aln, 1, "ACDEF", 4), expected_1)
        assert_allclose(get_positional_frequencies(aln, 2, "ACDEF", 4), expected_2)
        assert_allclose(get_positional_frequencies(aln, 3, "ACDEF", 4), expected_3)
        # extra characters (W) are silently ignored -- is this the desired
        # behavior?
        aln = ArrayAlignment(
            data={1: "WCDE", 2: "ADDE", 3: "AEED", 4: "AFEF"}, moltype=PROTEIN
        )
        expected_0 = array([3.0, 0.0, 0.0, 0.0, 0.0])
        assert_allclose(get_positional_frequencies(aln, 0, "ACDEF", 4), expected_0)
        # 20 residue amino acid alphabet
        aln = ArrayAlignment(
            data={1: "ACDE", 2: "ADDE", 3: "AEED", 4: "AFEF"}, moltype=PROTEIN
        )
        expected = array([4.0] + [0.0] * 19)
        assert_allclose(get_positional_frequencies(aln, 0, AAGapless, 4), expected)

    def test_validate_alphabet_invalid(self):
        """validate_alphabet: raises error on incompatible alpabet and freqs"""
        # len(alpha) > len(freqs)
        self.assertRaises(ValueError, validate_alphabet, "ABC", {"A": 0.5, "B": 0.5})
        self.assertRaises(ValueError, validate_alphabet, "ABCD", {"A": 0.5, "B": 0.5})
        # len(alpha) == len(freqs)
        self.assertRaises(ValueError, validate_alphabet, "AC", {"A": 0.5, "B": 0.5})
        # len(alpha) < len(freqs)
        self.assertRaises(ValueError, validate_alphabet, "A", {"A": 0.5, "B": 0.5})
        self.assertRaises(ValueError, validate_alphabet, "", {"A": 0.5, "B": 0.5})
        # different values, len(alpha) > len(freqs)
        self.assertRaises(ValueError, validate_alphabet, [1, 42, 3], {42: 0.5, 1: 0.5})
        self.assertRaises(
            ValueError, validate_alphabet, CharAlphabet("ABC"), {"A": 0.5, "C": 0.5}
        )

    def test_validate_alphabet_valid(self):
        """validate_alphabet: does nothing on compatible alpabet and freqs"""
        validate_alphabet("AB", {"A": 0.5, "B": 0.5})
        validate_alphabet(CharAlphabet("AB"), {"A": 0.5, "B": 0.5})
        validate_alphabet([1, 42, 8], {1: 0.5, 42: 0.25, 8: 0.25})

    def test_validate_position_invalid(self):
        """validate_position: raises error on invalid position"""
        self.assertRaises(ValueError, validate_position, self.dna_aln, 4)
        self.assertRaises(ValueError, validate_position, self.dna_aln, 42)
        self.assertRaises(ValueError, validate_position, self.dna_aln, -1)
        self.assertRaises(ValueError, validate_position, self.dna_aln, -199)

    def test_validate_position_valid(self):
        """validate_position: does nothing on valid position"""
        validate_position(self.dna_aln, 0)
        validate_position(self.dna_aln, 1)
        validate_position(self.dna_aln, 2)
        validate_position(self.dna_aln, 3)

    def test_validate_alignment(self):
        """validate_alignment: ValueError on bad alignment characters"""
        # ambiguous characters
        aln = ArrayAlignment(
            data={0: "BA", 1: "AC", 2: "CG", 3: "CT", 4: "TA"}, moltype=PROTEIN
        )
        self.assertRaises(ValueError, validate_alignment, aln)
        aln = ArrayAlignment(
            data={0: "NA", 1: "AC", 2: "CG", 3: "CT", 4: "TA"}, moltype=DNA
        )
        self.assertRaises(ValueError, validate_alignment, aln)
        aln = ArrayAlignment(
            data={0: "YA", 1: "AC", 2: "CG", 3: "CU", 4: "UA"}, moltype=RNA
        )
        self.assertRaises(ValueError, validate_alignment, aln)

        aln = ArrayAlignment(
            data={0: "AA", 1: "AC", 2: "CG", 3: "CT", 4: "TA"}, moltype=PROTEIN
        )
        validate_alignment(aln)
        aln = ArrayAlignment(
            data={0: "AA", 1: "AC", 2: "CG", 3: "CT", 4: "TA"}, moltype=DNA
        )
        validate_alignment(aln)
        aln = ArrayAlignment(
            data={0: "AA", 1: "AC", 2: "CG", 3: "CU", 4: "UA"}, moltype=RNA
        )
        validate_alignment(aln)

    def test_coevolve_functions_validate_alignment(self):
        """coevolve_*: functions run validate alignment"""
        aln = ArrayAlignment(
            data={"0": "BA", "1": "AC", "2": "CG", "3": "CT", "4": "TA"},
            moltype=PROTEIN,
        )
        self.assertRaises(ValueError, coevolve_pair, mi_pair, aln, 0, 1)
        self.assertRaises(ValueError, coevolve_position, mi_position, aln, 0)
        self.assertRaises(ValueError, coevolve_alignment, mi_alignment, aln)
        self.assertRaises(ValueError, coevolve_alignments, mi_alignment, aln, aln)

    def test_get_positional_probabilities_w_non_def_num_seqs(self):
        """get_positional_probabilities: works w/ non-def num_seqs"""
        freqs = [1.0, 2.0, 0.0]
        probs = [0.33, 0.33, 0.33]
        expected = array([0.444411, 0.218889, 0.300763])
        assert_allclose(get_positional_probabilities(freqs, probs, 3), expected)

    def test_get_dg(self):
        """get_dg: returns delta_g vector"""
        p = [0.1, 0.2, 0.3]
        a = [0.5, 0.6, 0.7]
        expected = [log(0.1 / 0.5), log(0.2 / 0.6), log(0.3 / 0.7)]
        assert_allclose(get_dg(p, a), expected)

    def test_get_dgg(self):
        """get_dgg: returns delta_delta_g value given two delta_g vectors"""
        v1 = array([0.05, 0.5, 0.1])
        v2 = array([0.03, 0.05, 0.1])
        expected = sqrt(sum((v1 - v2) * (v1 - v2))) / 100 * e
        assert_allclose(get_dgg(v1, v2), expected)

    def test_get_positional_probabilities_w_def_num_seqs(self):
        """get_positional_probabilities: works w/ num_seqs scaled to 100 (def)"""
        freqs = [15.0, 33.0, 52.0]
        probs = [0.33, 0.33, 0.33]
        expected = array([2.4990e-5, 0.0846, 3.8350e-5])
        assert_allclose(get_positional_probabilities(freqs, probs), expected, 0.001)

    def test_get_positional_probs_handles_rounding_error_in_freqs(self):
        """get_positional_probabilities: works w/ rounding error in freqs"""
        # Since freqs are scaled to scaled_aln_size, rounding error can cause
        # errors for positions that are perfectly controled. Testing here that
        # that value error is handled.
        # default scaled_aln_size
        freqs = [100.0000000001, 0.0, 0.0]
        probs = [0.33, 0.33, 0.33]
        expected = array([7.102218e-49, 4.05024e-18, 4.05024e-18])
        assert_allclose(get_positional_probabilities(freqs, probs), expected, rtol=1e-5)
        # value that is truely over raises an error
        freqs = [101.0000000001, 0.0, 0.0]
        probs = [0.33, 0.33, 0.33]
        self.assertRaises(ValueError, get_positional_probabilities, freqs, probs)
        # non-default scaled_aln_size
        freqs = [50.0000000001, 0.0, 0.0]
        probs = [0.33, 0.33, 0.33]
        expected = array([8.42747e-25, 2.01252e-9, 2.01252e-9])
        assert_allclose(
            get_positional_probabilities(freqs, probs, 50), expected, rtol=1e-5
        )
        # value that is truely over raises an error
        freqs = [51.0000000001, 0.0, 0.0]
        probs = [0.33, 0.33, 0.33]
        self.assertRaises(ValueError, get_positional_probabilities, freqs, probs, 50)

    def test_sca_input_validation(self):
        """sca_input_validation: handles sca-specific validation steps"""
        # moltype != PROTEIN makes background freqs required
        self.assertRaises(ValueError, sca_input_validation, self.dna_aln, cutoff=0.4)
        self.assertRaises(ValueError, sca_input_validation, self.rna_aln, cutoff=0.4)
        # no cutoff -> ValueError
        self.assertRaises(ValueError, sca_input_validation, self.protein_aln)
        # low cutoff -> ValueError
        self.assertRaises(
            ValueError, sca_input_validation, self.protein_aln, cutoff=-0.001
        )
        # high cutoff -> ValueError
        self.assertRaises(
            ValueError, sca_input_validation, self.protein_aln, cutoff=1.001
        )
        # good cut-off -> no error
        sca_input_validation(self.protein_aln, cutoff=0.50)
        sca_input_validation(self.protein_aln, cutoff=0.0)
        sca_input_validation(self.protein_aln, cutoff=1.0)

        # only bad alphabet -> ValueError
        self.assertRaises(
            ValueError, sca_input_validation, self.dna_aln, cutoff=0.5, alphabet="ABC"
        )
        # only bad background_freqs -> ValueError
        self.assertRaises(
            ValueError,
            sca_input_validation,
            self.dna_aln,
            cutoff=0.5,
            background_freqs={"A": 0.25, "C": 0.75},
        )
        # incompatible background_freqs & alphabet provided -> ValueError
        self.assertRaises(
            ValueError,
            sca_input_validation,
            self.dna_aln,
            cutoff=0.5,
            alphabet="ABC",
            background_freqs={"A": 0.25, "C": 0.75},
        )

        # default alphabet, background_freqs -> no error
        sca_input_validation(self.protein_aln, cutoff=0.50)
        # compatible non-default alphabet, backgorund_freqs -> no error
        sca_input_validation(
            self.dna_aln, cutoff=0.50, alphabet="A", background_freqs={"A": 1.0}
        )

        # Note: don't need a full set of tests of validate_alphabet here --
        # it's tested on it's own.

    def test_sca_pair_no_error(self):
        """sca_pair: returns w/o error"""
        r = sca_pair(
            self.dna_aln,
            1,
            0,
            cutoff=0.50,
            alphabet="ACGT",
            background_freqs=self.dna_base_freqs,
        )
        r = coevolve_pair(
            sca_pair,
            self.dna_aln,
            1,
            0,
            cutoff=0.50,
            alphabet="ACGT",
            background_freqs=self.dna_base_freqs,
        )

    def test_sca_pair_return_all(self):
        """sca_pair: handles return_all by returning lists of proper length"""
        # two allowed_perturbations
        a = "ACGT"
        aln = ArrayAlignment(
            data={0: "AA", 1: "AC", 2: "CG", 3: "CT", 4: "TA"}, moltype=DNA
        )
        actual = sca_pair(
            aln,
            0,
            1,
            cutoff=0.33,
            return_all=True,
            alphabet=a,
            background_freqs=self.dna_base_freqs,
        )
        self.assertEqual(len(actual), 2)
        self.assertEqual(actual[0][0], "A")
        self.assertEqual(actual[1][0], "C")
        # one allowed_perturbations
        a = "ACGT"
        aln = ArrayAlignment(
            data={0: "AA", 1: "AC", 2: "AG", 3: "CT", 4: "TA"}, moltype=DNA
        )
        actual = sca_pair(
            aln,
            0,
            1,
            0.33,
            return_all=True,
            alphabet=a,
            background_freqs=self.dna_base_freqs,
        )
        self.assertEqual(len(actual), 1)
        self.assertEqual(actual[0][0], "A")
        # zero allowed_perturbations
        actual = sca_pair(
            aln,
            0,
            1,
            1.0,
            return_all=True,
            alphabet=a,
            background_freqs=self.dna_base_freqs,
        )
        # expected = [('A',-1),('C',-1)]
        expected = DEFAULT_NULL_VALUE
        assert_allclose(actual, expected)

        # pos1 == pos2
        actual = sca_pair(
            aln,
            0,
            0,
            0.33,
            return_all=True,
            alphabet=a,
            background_freqs=self.dna_base_freqs,
        )

        actual = list(zip(*actual))
        expected = list(zip(*[("A", 2.40381185618)]))
        assert_equal(actual[0], expected[0])
        assert_allclose(actual[1], expected[1])

    def test_sca_pair_error(self):
        """sca_pair:returns w/ error when appropriate"""
        a = "ACGT"
        # pos1 out of range
        self.assertRaises(
            ValueError,
            coevolve_pair,
            sca_pair,
            self.dna_aln,
            100,
            1,
            cutoff=0.50,
            alphabet=a,
            background_freqs=self.dna_base_freqs,
        )
        # pos2 out of range
        self.assertRaises(
            ValueError,
            coevolve_pair,
            sca_pair,
            self.dna_aln,
            0,
            100,
            cutoff=0.50,
            alphabet=a,
            background_freqs=self.dna_base_freqs,
        )
        # pos1 & pos2 out of range
        self.assertRaises(
            ValueError,
            coevolve_pair,
            sca_pair,
            self.dna_aln,
            100,
            100,
            cutoff=0.50,
            alphabet=a,
            background_freqs=self.dna_base_freqs,
        )

        # bad cut-off
        self.assertRaises(
            ValueError,
            coevolve_pair,
            sca_pair,
            self.dna_aln,
            0,
            1,
            cutoff=1.2,
            alphabet=a,
            background_freqs=self.dna_base_freqs,
        )

        # incompatible alphabet and background freqs
        self.assertRaises(
            ValueError,
            coevolve_pair,
            sca_pair,
            self.dna_aln,
            0,
            1,
            cutoff=0.2,
            alphabet=a,
        )
        self.assertRaises(
            ValueError,
            coevolve_pair,
            sca_pair,
            self.dna_aln,
            0,
            1,
            cutoff=0.2,
            alphabet="ACGTBC",
            background_freqs=self.dna_base_freqs,
        )

    def test_sca_position_no_error(self):
        """sca_position: returns w/o error"""
        r = sca_position(
            self.dna_aln, 1, 0.50, alphabet="ACGT", background_freqs=self.dna_base_freqs
        )
        # sanity check -- coupling w/ self
        assert_allclose(r[1], 3.087, 0.01)
        r = sca_position(
            self.dna_aln_gapped,
            1,
            0.50,
            alphabet="ACGT",
            background_freqs=self.dna_base_freqs,
        )
        assert_allclose(r[1], 3.387, 0.01)

        # same tests, but called via coevolve_position
        r = coevolve_position(
            sca_position,
            self.dna_aln,
            1,
            cutoff=0.50,
            alphabet="ACGT",
            background_freqs=self.dna_base_freqs,
        )
        # sanity check -- coupling w/ self
        assert_allclose(r[1], 3.087, 0.01)
        r = coevolve_position(
            sca_position,
            self.dna_aln_gapped,
            1,
            cutoff=0.50,
            alphabet="ACGT",
            background_freqs=self.dna_base_freqs,
        )
        # sanity check -- coupling w/ self
        assert_allclose(r[1], 3.387, 0.01)

    def test_sca_position_error(self):
        """sca_position: returns w/ error when appropriate"""
        a = "ACGT"
        # position out of range
        self.assertRaises(
            ValueError,
            coevolve_position,
            sca_position,
            self.dna_aln,
            100,
            cutoff=0.50,
            alphabet=a,
            background_freqs=self.dna_base_freqs,
        )
        # bad cutoff
        self.assertRaises(
            ValueError,
            coevolve_position,
            sca_position,
            self.dna_aln,
            1,
            cutoff=-8.2,
            alphabet=a,
            background_freqs=self.dna_base_freqs,
        )

        # incompatible alphabet and background freqs
        self.assertRaises(
            ValueError,
            coevolve_position,
            sca_position,
            self.dna_aln,
            0,
            cutoff=0.2,
            alphabet=a,
        )
        self.assertRaises(
            ValueError,
            coevolve_position,
            sca_position,
            self.dna_aln,
            0,
            cutoff=0.2,
            alphabet="ACGTBC",
            background_freqs=self.dna_base_freqs,
        )

    def test_sca_position_returns_same_as_sca_pair(self):
        """sca_position: returns same as sca_pair called on each pos"""
        expected = []
        for i in range(len(self.dna_aln)):
            expected.append(
                sca_pair(
                    self.dna_aln,
                    1,
                    i,
                    0.50,
                    alphabet="ACGT",
                    background_freqs=self.dna_base_freqs,
                )
            )
        actual = sca_position(
            self.dna_aln, 1, 0.50, alphabet="ACGT", background_freqs=self.dna_base_freqs
        )
        assert_allclose(actual, expected)
        # change some of the defaults to make sure they make it through
        bg_freqs = {"A": 0.50, "C": 0.50}
        expected = []
        for i in range(len(self.dna_aln)):
            expected.append(
                sca_pair(
                    self.dna_aln,
                    1,
                    i,
                    0.50,
                    alphabet="AC",
                    null_value=52.0,
                    scaled_aln_size=20,
                    background_freqs=bg_freqs,
                )
            )
        actual = sca_position(
            self.dna_aln,
            1,
            0.50,
            alphabet="AC",
            null_value=52.0,
            scaled_aln_size=20,
            background_freqs=bg_freqs,
        )
        assert_allclose(actual, expected)

    def test_sca_alignment_no_error(self):
        """sca_alignment: returns w/o error"""
        r = sca_alignment(
            self.dna_aln, 0.50, alphabet="ACGT", background_freqs=self.dna_base_freqs
        )
        # sanity check -- coupling w/ self
        assert_allclose(r[0][0], 2.32222608171)

        # same test, but called via coevolve_alignment
        r = coevolve_alignment(
            sca_alignment,
            self.dna_aln,
            cutoff=0.50,
            alphabet="ACGT",
            background_freqs=self.dna_base_freqs,
        )
        # sanity check -- coupling w/ self
        assert_allclose(r[0][0], 2.32222608171)

    def test_sca_alignment_error(self):
        """sca_alignment: returns w/ error when appropriate"""
        a = "ACGT"
        # incompatible alphabet and background freqs
        self.assertRaises(
            ValueError,
            coevolve_position,
            sca_position,
            self.dna_aln,
            0,
            cutoff=0.2,
            alphabet=a,
        )
        self.assertRaises(
            ValueError,
            coevolve_position,
            sca_position,
            self.dna_aln,
            0,
            cutoff=0.2,
            alphabet="ACGTBC",
            background_freqs=self.dna_base_freqs,
        )

    def test_sca_alignment_returns_same_as_sca_position(self):
        """sca_alignment: returns same as sca_position on every position"""
        expected = []
        for i in range(len(self.dna_aln)):
            expected.append(
                sca_position(
                    self.dna_aln,
                    i,
                    0.50,
                    alphabet="ACGT",
                    background_freqs=self.dna_base_freqs,
                )
            )
        actual = sca_alignment(
            self.dna_aln, 0.50, alphabet="ACGT", background_freqs=self.dna_base_freqs
        )
        assert_allclose(actual, expected)
        # change some of the defaults to make sure they make it through
        bg_freqs = {"A": 0.50, "C": 0.50}
        expected = []
        for i in range(len(self.dna_aln)):
            expected.append(
                sca_position(
                    self.dna_aln,
                    i,
                    0.50,
                    alphabet="AC",
                    null_value=52.0,
                    scaled_aln_size=20,
                    background_freqs=bg_freqs,
                )
            )
        actual = sca_alignment(
            self.dna_aln,
            0.50,
            alphabet="AC",
            null_value=52.0,
            scaled_aln_size=20,
            background_freqs=bg_freqs,
        )
        assert_allclose(actual, expected)

    def test_sca_pair_gpcr(self):
        """sca_pair: reproduces several GPCR data from Suel et al., 2003"""
        assert_allclose(sca_pair(self.gpcr_aln, 295, 18, 0.32), 0.12, 0.1)
        assert_allclose(sca_pair(self.gpcr_aln, 295, 124, 0.32), 1.86, 0.1)
        assert_allclose(sca_pair(self.gpcr_aln, 295, 304, 0.32), 0.3, 0.1)
        # covariation w/ self
        assert_allclose(sca_pair(self.gpcr_aln, 295, 295, 0.32), 7.70358628)

    def test_sca_position_gpcr(self):
        """sca_position: reproduces several GPCR data from Suel et al., 2003"""
        if not self.run_slow_tests:
            return
        vector = sca_position(self.gpcr_aln, 295, 0.32)
        assert_allclose(vector[18], 0.12, 0.1)
        assert_allclose(vector[124], 1.86, 0.1)
        assert_allclose(vector[304], 0.3, 0.1)
        # covariation w/ self == null_value
        assert_allclose(vector[295], nan)

    def test_ltm_to_symmetric(self):
        """ltm_to_symmetric: making ltm matrices symmetric functions"""
        m = arange(9).reshape((3, 3))
        expected = [[0, 3, 6], [3, 4, 7], [6, 7, 8]]
        assert_equal(ltm_to_symmetric(m), expected)
        # non-square matrices not supported
        self.assertRaises(AssertionError, ltm_to_symmetric, arange(10).reshape(5, 2))
        self.assertRaises(AssertionError, ltm_to_symmetric, arange(10).reshape(2, 5))

    def test_merge_alignments(self):
        """merging alignments of same moltype functions as expected"""
        # PROTEIN
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AD"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1": "EF", "2": "EG"}, moltype=PROTEIN)
        combined_aln = ArrayAlignment(data={"1": "ACEF", "2": "ADEG"}, moltype=PROTEIN)
        actual = merge_alignments(aln1, aln2)
        self.assertEqual(actual, combined_aln)
        self.assertEqual(actual.moltype, PROTEIN)
        # RNA
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AU"}, moltype=RNA)
        aln2 = ArrayAlignment(data={"1": "GG", "2": "UG"}, moltype=RNA)
        combined_aln = ArrayAlignment(data={"1": "ACGG", "2": "AUUG"}, moltype=RNA)
        actual = merge_alignments(aln1, aln2)
        self.assertEqual(actual, combined_aln)
        self.assertEqual(actual.moltype, RNA)
        # DNA
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AT"}, moltype=DNA)
        aln2 = ArrayAlignment(data={"1": "GG", "2": "TG"}, moltype=DNA)
        combined_aln = ArrayAlignment(data={"1": "ACGG", "2": "ATTG"}, moltype=DNA)
        actual = merge_alignments(aln1, aln2)
        self.assertEqual(actual, combined_aln)
        self.assertEqual(actual.moltype, DNA)

    def test_merge_alignments_ignores_id_following_plus(self):
        """merge_alignments ignores all seq id characters after '+'"""
        aln1 = ArrayAlignment(data={"1+a": "AC", "2+b": "AD"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1 + c": "EFW", "2 + d": "EGY"}, moltype=PROTEIN)
        combined_aln = ArrayAlignment(
            data={"1": "ACEFW", "2": "ADEGY"}, moltype=PROTEIN
        )
        self.assertEqual(merge_alignments(aln1, aln2), combined_aln)
        # not all ids have a +
        aln1 = ArrayAlignment(data={"1": "AC", "2+b": "AD"}, moltype=PROTEIN)
        aln2 = ArrayAlignment(data={"1+c": "EFW", "2": "EGY"}, moltype=PROTEIN)
        combined_aln = ArrayAlignment(
            data={"1": "ACEFW", "2": "ADEGY"}, moltype=PROTEIN
        )
        self.assertEqual(merge_alignments(aln1, aln2), combined_aln)

    def test_merge_alignments_different_moltype(self):
        """merging alignments of different moltype functions as expected"""
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AU"}, moltype=RNA)
        aln2 = ArrayAlignment(data={"1": "EF", "2": "EG"}, moltype=PROTEIN)
        combined_aln = ArrayAlignment(data={"1": "ACEF", "2": "AUEG"})
        self.assertEqual(merge_alignments(aln1, aln2), combined_aln)
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AT"}, moltype=DNA)
        aln2 = ArrayAlignment(data={"1": "EF", "2": "EG"}, moltype=PROTEIN)
        combined_aln = ArrayAlignment(data={"1": "ACEF", "2": "ATEG"})
        self.assertEqual(merge_alignments(aln1, aln2), combined_aln)
        aln1 = ArrayAlignment(data={"1": "AC", "2": "AT"}, moltype=DNA)
        aln2 = ArrayAlignment(data={"1": "UC", "2": "UG"}, moltype=RNA)
        combined_aln = ArrayAlignment(data={"1": "ACUC", "2": "ATUG"})
        self.assertEqual(merge_alignments(aln1, aln2), combined_aln)

    def test_n_random_seqs(self):
        """n_random_seqs: functions as expected"""
        aln1 = make_aligned_seqs(
            data=list(zip(list("abcd"), ["AA", "AC", "DD", "GG"])),
            moltype=PROTEIN,
            array_align=True,
        )
        # Number of returned sequences correct
        self.assertEqual(n_random_seqs(aln1, 1).num_seqs, 1)
        self.assertEqual(n_random_seqs(aln1, 2).num_seqs, 2)
        self.assertEqual(n_random_seqs(aln1, 3).num_seqs, 3)
        self.assertEqual(n_random_seqs(aln1, 4).num_seqs, 4)

        # Sequences are correct
        new_aln = n_random_seqs(aln1, 3)
        self.assertEqual(new_aln.num_seqs, 3)
        for n in new_aln.names:
            self.assertEqual(new_aln.get_seq(n), aln1.get_seq(n))

        # Objects are equal when all are requested
        self.assertEqual(n_random_seqs(aln1, 4), aln1)

        # Objects are not equal when subset are requested
        self.assertNotEqual(n_random_seqs(aln1, 3), aln1)

        # In 1000 iterations, we get at least one different alignment --
        # this tests the random selection
        different = False
        new_aln = n_random_seqs(aln1, 2)
        for i in range(1000):
            new_aln2 = n_random_seqs(aln1, 2)
            if new_aln != new_aln2:
                different = True
                break
        self.assertTrue(different)


class AncestorCoevolve(TestCase):
    """Tests of the ancestral state method for detecting coevolution"""

    def setUp(self):
        """ """
        # t1, ancestral_states1, and aln1_* are used to test that when
        # alternate seqs are used with the same tree and ancestral_states,
        # the results vary when appropriate
        self.t1 = make_tree(
            treestring="((A:0.5,B:0.5):0.5,(C:0.5,(D:0.5,E:0.5):0.5):0.5);"
        )
        self.ancestral_states1 = ArrayAlignment(
            data={"root": "AAA", "edge.0": "AAA", "edge.1": "AAA", "edge.2": "AAA"},
            moltype=PROTEIN,
        )
        self.ancestral_states1_w_gaps = ArrayAlignment(
            data={"root": "AAA", "edge.0": "AAA", "edge.1": "A-A", "edge.2": "AA-"},
            moltype=PROTEIN,
        )

        # no correlated changes count
        self.aln1_1 = ArrayAlignment(
            data={"A": "AAC", "B": "AAD", "C": "AAA", "D": "AAE", "E": "AFA"},
            moltype=PROTEIN,
        )
        # 1 correlated change count
        self.aln1_2 = ArrayAlignment(
            data={"A": "AAC", "B": "AAD", "C": "AAA", "D": "AEE", "E": "AFF"},
            moltype=PROTEIN,
        )
        # 1 different correlated change count
        self.aln1_3 = ArrayAlignment(
            data={"A": "AAC", "B": "AAD", "C": "AAA", "D": "AGE", "E": "AFH"},
            moltype=PROTEIN,
        )
        # 3 correlated change counts
        self.aln1_4 = ArrayAlignment(
            data={"A": "AAC", "B": "AGD", "C": "AAA", "D": "AGE", "E": "AFH"},
            moltype=PROTEIN,
        )
        # 8 correlated change counts
        self.aln1_5 = ArrayAlignment(
            data={"A": "YYC", "B": "HGD", "C": "AAA", "D": "AGE", "E": "AFH"},
            moltype=PROTEIN,
        )
        self.aln1_w_gaps = ArrayAlignment(
            data={"A": "AAC", "B": "AAD", "C": "AAA", "D": "AG-", "E": "A-H"},
            moltype=PROTEIN,
        )

        # t2, ancestral_states2_*, and aln2 are used to test that when
        # alternate ancestral states are used with the same aln and tree,
        # the results vary when appropriate
        self.t2 = make_tree(treestring="(A:0.5,B:0.5,C:0.5);")
        self.ancestral_states2_1 = ArrayAlignment(data={"root": "AA"}, moltype=PROTEIN)
        self.ancestral_states2_2 = ArrayAlignment(data={"root": "CC"}, moltype=PROTEIN)
        self.ancestral_states2_3 = ArrayAlignment(data={"root": "EF"}, moltype=PROTEIN)
        self.aln2 = ArrayAlignment(
            data={"A": "AA", "B": "CC", "C": "CA"}, moltype=PROTEIN
        )

        # t3_*, ancestral_states3, and aln3 are used to test that when
        # alternate trees are used with the same aln and ancestral_states,
        # the results vary when appropriate
        self.t3_1 = make_tree(treestring="(A:0.5,(B:0.5,C:0.5):0.5);")
        self.t3_2 = make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);")
        self.ancestral_states3 = ArrayAlignment(
            data={"root": "CC", "edge.0": "AD"}, moltype=PROTEIN
        )
        self.aln3 = ArrayAlignment(
            data={"A": "AC", "B": "CA", "C": "CC"}, moltype=PROTEIN
        )

    def test_validate_ancestral_seqs_invalid(self):
        """validate_ancestral_seqs: ValueError on incompatible anc. seqs & tree"""
        # edge missing
        aln = ArrayAlignment(data={"A": "AC", "B": "CA", "C": "CC"}, moltype=PROTEIN)
        self.assertRaises(
            ValueError,
            validate_ancestral_seqs,
            aln,
            tree=make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);"),
            ancestral_seqs=ArrayAlignment(data={"root": "AA"}, moltype=PROTEIN),
        )
        # root missing
        self.assertRaises(
            ValueError,
            validate_ancestral_seqs,
            aln,
            tree=make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);"),
            ancestral_seqs=ArrayAlignment(data={"edge.0": "AA"}, moltype=PROTEIN),
        )
        # correct numSeqs but wrong names
        self.assertRaises(
            ValueError,
            validate_ancestral_seqs,
            aln,
            tree=make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);"),
            ancestral_seqs=ArrayAlignment(
                data={"root": "AA", "edge.1": "AA"}, moltype=PROTEIN
            ),
        )
        self.assertRaises(
            ValueError,
            validate_ancestral_seqs,
            aln,
            tree=make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);"),
            ancestral_seqs=ArrayAlignment(
                data={"r": "AA", "edge.0": "AA"}, moltype=PROTEIN
            ),
        )
        self.assertRaises(
            ValueError,
            validate_ancestral_seqs,
            aln,
            tree=make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);"),
            ancestral_seqs=ArrayAlignment(data={"r": "AA", "e": "AA"}, moltype=PROTEIN),
        )
        # different tree: invalid
        aln = ArrayAlignment(
            data={"A": "AC", "B": "CA", "C": "CC", "D": "DD"}, moltype=PROTEIN
        )
        self.assertRaises(
            ValueError,
            validate_ancestral_seqs,
            aln,
            tree=make_tree(treestring="((A:0.5,B:0.5):0.5,(C:0.5,D:0.5):0.5);"),
            ancestral_seqs=ArrayAlignment(
                data={"root": "AA", "e": "AA", "edge.1": "AA"}, moltype=PROTEIN
            ),
        )

    def test_validate_ancestral_seqs_valid(self):
        """validate_ancestral_seqs: does nothing on compatible anc. seqs & tree"""
        aln = ArrayAlignment(data={"A": "AC", "B": "CA", "C": "CC"}, moltype=PROTEIN)
        # valid data -> no error
        validate_ancestral_seqs(
            aln,
            tree=make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);"),
            ancestral_seqs=ArrayAlignment(
                data={"root": "AA", "edge.0": "AA"}, moltype=PROTEIN
            ),
        )
        # different tree: valid
        aln = ArrayAlignment(
            data={"A": "AC", "B": "CA", "C": "CC", "D": "DD"}, moltype=PROTEIN
        )
        validate_ancestral_seqs(
            aln,
            tree=make_tree(treestring="((A:0.5,B:0.5):0.5,(C:0.5,D:0.5):0.5);"),
            ancestral_seqs=ArrayAlignment(
                data={"root": "AA", "edge.0": "AA", "edge.1": "AA"}, moltype=PROTEIN
            ),
        )

    def test_ancestral_states_input_validation(self):
        """ancestral_states_input_validation: all validation steps performed"""
        aln = ArrayAlignment(
            data={"A": "AC", "B": "CA", "C": "CC", "D": "DD"}, moltype=PROTEIN
        )
        # incompatible tree and ancestral states (more thorough testing in
        # test_validate_ancestral_seqs)
        self.assertRaises(
            ValueError,
            ancestral_states_input_validation,
            aln,
            tree=make_tree(treestring="((A:0.5,B:0.5):0.5,(C:0.5,D:0.5):0.5);"),
            ancestral_seqs=ArrayAlignment(
                data={"root": "AA", "e": "AA", "edge.1": "AA"}, moltype=PROTEIN
            ),
        )
        # no tree provided
        self.assertRaises(
            ValueError,
            ancestral_states_input_validation,
            aln,
            ancestral_seqs=ArrayAlignment(
                data={"root": "AA", "e": "AA", "edge.1": "AA"}, moltype=PROTEIN
            ),
        )
        # incompatible tree and alignment (more tests in test_validate_tree)
        aln = ArrayAlignment(data={"A": "AC", "B": "CA", "C": "CC"}, moltype=PROTEIN)
        self.assertRaises(
            ValueError,
            ancestral_states_input_validation,
            aln,
            tree=make_tree(treestring="((A:0.5,B:0.5):0.5,(C:0.5,D:0.5):0.5);"),
        )

    def test_validate_tree_valid(self):
        """validate_tree: does nothing on compatible tree and aln"""
        t = make_tree(treestring="((A:0.5,B:0.5):0.5,(C:0.5,D:0.5):0.5);")
        aln = ArrayAlignment(
            data={"A": "AC", "B": "CA", "C": "CC", "D": "DD"}, moltype=PROTEIN
        )
        validate_tree(aln, t)
        t = make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);")
        aln = ArrayAlignment(data={"A": "AC", "B": "CA", "C": "CC"}, moltype=PROTEIN)
        validate_tree(aln, t)

    def test_validate_tree_invalid(self):
        """validate_tree: raises ValueError on incompatible tree and aln"""
        # different scale tree and aln
        t = make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);")
        aln = ArrayAlignment(
            data={"A": "AC", "B": "CA", "C": "CC", "D": "DD"}, moltype=PROTEIN
        )
        self.assertRaises(ValueError, validate_tree, aln, t)
        t = make_tree(treestring="((A:0.5,B:0.5):0.5,(C:0.5,D:0.5):0.5);")
        aln = ArrayAlignment(data={"A": "AC", "B": "CA", "C": "CC"}, moltype=PROTEIN)
        self.assertRaises(ValueError, validate_tree, aln, t)
        # same scale tree and aln, but different names
        t = make_tree(treestring="((A:0.5,B:0.5):0.5,(C:0.5,Dee:0.5):0.5);")
        aln = ArrayAlignment(
            data={"A": "AC", "B": "CA", "C": "CC", "D": "DD"}, moltype=PROTEIN
        )
        self.assertRaises(ValueError, validate_tree, aln, t)

    def test_get_ancestral_seqs(self):
        """get_ancestral_seqs: returns valid collection of ancestral seqs"""
        t = make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);")
        aln = ArrayAlignment(data={"A": "AA", "B": "AA", "C": "AC"}, moltype=PROTEIN)
        expected = ArrayAlignment(data={"root": "AA", "edge.0": "AA"}, moltype=PROTEIN)
        self.assertEqual(get_ancestral_seqs(aln, t, optimise=False), expected)
        t = make_tree(treestring="(A:0.5,B:0.5,C:0.5);")
        aln = ArrayAlignment(data={"A": "AA", "B": "AA", "C": "AC"}, moltype=PROTEIN)
        expected = ArrayAlignment(data={"root": "AA"}, moltype=PROTEIN)
        self.assertEqual(get_ancestral_seqs(aln, t, optimise=False), expected)

        t = make_tree(
            treestring="(((A1:0.5,A2:0.5):0.5,B:0.5):0.5,\
            (C:0.5,D:0.5):0.5);"
        )
        aln = ArrayAlignment(
            data={"A1": "AD", "A2": "AD", "B": "AC", "C": "AC", "D": "AC"},
            moltype=PROTEIN,
        )
        expected = ArrayAlignment(
            data={"root": "AC", "edge.0": "AD", "edge.1": "AC", "edge.2": "AC"},
            moltype=PROTEIN,
        )
        self.assertEqual(get_ancestral_seqs(aln, t, optimise=False), expected)

    def test_get_ancestral_seqs_handles_gaps(self):
        """get_ancestral_seqs: handles gaps"""
        # gaps handled OK
        t = make_tree(treestring="(A:0.5,B:0.5,C:0.5);")
        aln = ArrayAlignment(data={"A": "A-", "B": "AA", "C": "AA"}, moltype=PROTEIN)
        expected = ArrayAlignment(data={"root": "AA"}, moltype=PROTEIN)
        self.assertEqual(get_ancestral_seqs(aln, t, optimise=False), expected)

    def test_get_ancestral_seqs_handles_ambiguous_residues(self):
        """get_ancestral_seqs: handles ambiguous residues"""
        # Non-canonical residues handled OK
        t = make_tree(treestring="(A:0.5,B:0.5,C:0.5);")
        aln = ArrayAlignment(data={"A": "AX", "B": "Z-", "C": "BC"}, moltype=PROTEIN)
        actual = get_ancestral_seqs(aln, t, optimise=False)
        self.assertEqual(len(actual), 2)
        self.assertEqual(actual.num_seqs, 1)

    def test_ancestral_state_alignment_handles_ancestral_state_calc(self):
        """ancestral_state_alignment: functions when calc'ing ancestral states"""
        t = make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);")
        aln = ArrayAlignment(data={"A": "AA", "B": "AA", "C": "AC"}, moltype=PROTEIN)
        assert_equal(ancestral_state_alignment(aln, t), [[0, 0], [0, 2]])
        # non-bifurcating tree
        t = make_tree(treestring="(A:0.5,B:0.5,C:0.5);")
        aln = ArrayAlignment(data={"A": "AA", "B": "AA", "C": "AC"}, moltype=PROTEIN)
        assert_equal(ancestral_state_alignment(aln, t), [[0, 0], [0, 2]])

    def test_ancestral_state_position_handles_ancestral_state_calc(self):
        """ancestral_state_position: functions when calc'ing ancestral states"""
        t = make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);")
        aln = ArrayAlignment(data={"A": "AA", "B": "AA", "C": "AC"}, moltype=PROTEIN)
        assert_equal(ancestral_state_position(aln, t, 0), [0, 0])
        assert_equal(ancestral_state_position(aln, t, 1), [0, 2])

    def test_ancestral_state_pair_handles_ancestral_state_calc(self):
        """ancestral_state_position: functions when calc'ing ancestral states"""
        t = make_tree(treestring="((A:0.5,B:0.5):0.5,C:0.5);")
        aln = ArrayAlignment(data={"A": "AA", "B": "AA", "C": "AC"}, moltype=PROTEIN)
        self.assertEqual(ancestral_state_pair(aln, t, 0, 0), 0)
        self.assertEqual(ancestral_state_pair(aln, t, 0, 1), 0)
        self.assertEqual(ancestral_state_pair(aln, t, 1, 1), 2)
        self.assertEqual(ancestral_state_pair(aln, t, 1, 0), 0)

    def test_ancestral_state_alignment_no_error_on_gap(self):
        """ancestral_state_alignment: return w/o error with gapped seqs"""
        ancestral_state_alignment(self.aln1_w_gaps, self.t1, self.ancestral_states1)
        ancestral_state_alignment(self.aln1_1, self.t1, self.ancestral_states1_w_gaps)

    def test_ancestral_state_methods_handle_bad_ancestor_aln(self):
        """ancestral state methods raise error on bad ancestor alignment"""
        # bad length and seq names
        self.assertRaises(
            ValueError,
            coevolve_alignment,
            ancestral_state_alignment,
            self.aln1_2,
            tree=self.t1,
            ancestral_seqs=self.ancestral_states2_1,
        )
        self.assertRaises(
            ValueError,
            coevolve_position,
            ancestral_state_position,
            self.aln1_2,
            0,
            tree=self.t1,
            ancestral_seqs=self.ancestral_states2_1,
        )
        self.assertRaises(
            ValueError,
            coevolve_pair,
            ancestral_state_pair,
            self.aln1_2,
            0,
            1,
            tree=self.t1,
            ancestral_seqs=self.ancestral_states2_1,
        )
        # bad seq names
        self.assertRaises(
            ValueError,
            coevolve_alignment,
            ancestral_state_alignment,
            self.aln1_2,
            tree=self.t1,
            ancestral_seqs=self.aln1_2,
        )
        self.assertRaises(
            ValueError,
            coevolve_position,
            ancestral_state_position,
            self.aln1_2,
            0,
            tree=self.t1,
            ancestral_seqs=self.aln1_2,
        )
        self.assertRaises(
            ValueError,
            coevolve_pair,
            ancestral_state_pair,
            self.aln1_2,
            0,
            1,
            tree=self.t1,
            ancestral_seqs=self.aln1_2,
        )
        # bad length
        a = ArrayAlignment(
            data={"root": "AC", "edge.0": "AD", "edge.1": "AA", "edge.2": "EE"}
        )
        self.assertRaises(
            ValueError,
            coevolve_alignment,
            ancestral_state_alignment,
            self.aln1_2,
            tree=self.t1,
            ancestral_seqs=a,
        )
        self.assertRaises(
            ValueError,
            coevolve_position,
            ancestral_state_position,
            self.aln1_2,
            0,
            tree=self.t1,
            ancestral_seqs=a,
        )
        self.assertRaises(
            ValueError,
            coevolve_pair,
            ancestral_state_pair,
            self.aln1_2,
            0,
            1,
            tree=self.t1,
            ancestral_seqs=a,
        )

    def test_ancestral_states_methods_handle_bad_position_numbers(self):
        """coevolve_* w/ ancestral_states raise ValueError on bad position"""

        self.assertRaises(
            ValueError,
            coevolve_position,
            ancestral_state_position,
            self.aln1_2,
            42,
            tree=self.t1,
            ancestral_states=self.ancestral_states2_1,
        )
        self.assertRaises(
            ValueError,
            coevolve_pair,
            ancestral_state_pair,
            self.aln1_2,
            0,
            42,
            tree=self.t1,
            ancestral_states=self.ancestral_states2_1,
        )
        self.assertRaises(
            ValueError,
            coevolve_pair,
            ancestral_state_pair,
            self.aln1_2,
            42,
            0,
            tree=self.t1,
            ancestral_states=self.ancestral_states2_1,
        )

    def test_ancestral_state_alignment_non_bifurcating_tree(self):
        """ancestral_state_alignment: handles non-bifurcating tree correctly"""
        assert_equal(
            ancestral_state_alignment(self.aln2, self.t2, self.ancestral_states2_3),
            [[9, 9], [9, 9]],
        )

    def test_ancestral_state_alignment_bifurcating_tree(self):
        """ancestral_state_alignment: handles bifurcating tree correctly"""
        assert_allclose(
            ancestral_state_alignment(self.aln1_5, self.t1, self.ancestral_states1),
            [[5, 5, 5], [5, 11.6, 11.6], [5, 11.6, 11.6]],
        )

    def test_ancestral_state_alignment_ancestor_difference(self):
        """ancestral_state_alignment: different ancestor -> different result"""
        # ancestral_states2_1
        assert_equal(
            ancestral_state_alignment(self.aln2, self.t2, self.ancestral_states2_1),
            [[5, 2], [2, 2]],
        )
        # ancestral_states2_2
        assert_equal(
            ancestral_state_alignment(self.aln2, self.t2, self.ancestral_states2_2),
            [[2, 2], [2, 5]],
        )
        # ancestral_states2_3
        assert_equal(
            ancestral_state_alignment(self.aln2, self.t2, self.ancestral_states2_3),
            [[9, 9], [9, 9]],
        )

    def test_ancestral_state_position_ancestor_difference(self):
        """ancestral_state_position: difference_ancestor -> different result"""
        # ancestral_states2_1
        assert_equal(
            ancestral_state_position(self.aln2, self.t2, 0, self.ancestral_states2_1),
            [5, 2],
        )
        assert_equal(
            ancestral_state_position(self.aln2, self.t2, 1, self.ancestral_states2_1),
            [2, 2],
        )
        # ancestral_states2_2
        assert_equal(
            ancestral_state_position(self.aln2, self.t2, 0, self.ancestral_states2_2),
            [2, 2],
        )
        assert_equal(
            ancestral_state_position(self.aln2, self.t2, 1, self.ancestral_states2_2),
            [2, 5],
        )
        # ancestral_states2_3
        assert_equal(
            ancestral_state_position(self.aln2, self.t2, 0, self.ancestral_states2_3),
            [9, 9],
        )
        assert_equal(
            ancestral_state_position(self.aln2, self.t2, 1, self.ancestral_states2_3),
            [9, 9],
        )

    def test_ancestral_state_pair_ancestor_difference(self):
        """ancestral_state_pair: difference_ancestor -> different result"""
        # ancestral_states2_1
        self.assertEqual(
            ancestral_state_pair(self.aln2, self.t2, 0, 0, self.ancestral_states2_1), 5
        )
        self.assertEqual(
            ancestral_state_pair(self.aln2, self.t2, 0, 1, self.ancestral_states2_1), 2
        )
        self.assertEqual(
            ancestral_state_pair(self.aln2, self.t2, 1, 1, self.ancestral_states2_1), 2
        )
        self.assertEqual(
            ancestral_state_pair(self.aln2, self.t2, 1, 0, self.ancestral_states2_1), 2
        )
        # ancestral_states2_2
        self.assertEqual(
            ancestral_state_pair(self.aln2, self.t2, 0, 0, self.ancestral_states2_2), 2
        )
        self.assertEqual(
            ancestral_state_pair(self.aln2, self.t2, 0, 1, self.ancestral_states2_2), 2
        )
        self.assertEqual(
            ancestral_state_pair(self.aln2, self.t2, 1, 1, self.ancestral_states2_2), 5
        )
        self.assertEqual(
            ancestral_state_pair(self.aln2, self.t2, 1, 0, self.ancestral_states2_2), 2
        )
        # ancestral_states2_3
        self.assertEqual(
            ancestral_state_pair(self.aln2, self.t2, 0, 0, self.ancestral_states2_3), 9
        )
        self.assertEqual(
            ancestral_state_pair(self.aln2, self.t2, 0, 1, self.ancestral_states2_3), 9
        )
        self.assertEqual(
            ancestral_state_pair(self.aln2, self.t2, 1, 1, self.ancestral_states2_3), 9
        )
        self.assertEqual(
            ancestral_state_pair(self.aln2, self.t2, 1, 0, self.ancestral_states2_3), 9
        )

    def test_ancestral_state_alignment_tree_difference(self):
        """ancestral_state_alignment: different result on different tree"""
        # tree: t3_1
        assert_equal(
            ancestral_state_alignment(self.aln3, self.t3_1, self.ancestral_states3),
            [[7, 5], [5, 5]],
        )
        # tree: t3_2
        assert_equal(
            ancestral_state_alignment(self.aln3, self.t3_2, self.ancestral_states3),
            [[2, 2], [2, 5]],
        )

    def test_ancestral_state_position_tree_difference(self):
        """ancestral_state_position: different result on different tree"""
        # tree: t3_1
        assert_equal(
            ancestral_state_position(self.aln3, self.t3_1, 0, self.ancestral_states3),
            [7, 5],
        )
        assert_equal(
            ancestral_state_position(self.aln3, self.t3_1, 1, self.ancestral_states3),
            [5, 5],
        )
        # tree: t3_2
        assert_equal(
            ancestral_state_position(self.aln3, self.t3_2, 0, self.ancestral_states3),
            [2, 2],
        )
        assert_equal(
            ancestral_state_position(self.aln3, self.t3_2, 1, self.ancestral_states3),
            [2, 5],
        )

    def test_ancestral_state_pair_tree_difference(self):
        """ancestral_state_pair: different result on different tree"""
        # tree: t3_1
        assert_allclose(
            ancestral_state_pair(self.aln3, self.t3_1, 0, 1, self.ancestral_states3), 5
        )
        assert_allclose(
            ancestral_state_pair(self.aln3, self.t3_1, 1, 0, self.ancestral_states3), 5
        )
        assert_allclose(
            ancestral_state_pair(self.aln3, self.t3_1, 0, 0, self.ancestral_states3), 7
        )
        assert_allclose(
            ancestral_state_pair(self.aln3, self.t3_1, 1, 1, self.ancestral_states3), 5
        )
        # tree: t3_2
        assert_allclose(
            ancestral_state_pair(self.aln3, self.t3_2, 0, 1, self.ancestral_states3), 2
        )
        assert_allclose(
            ancestral_state_pair(self.aln3, self.t3_2, 1, 0, self.ancestral_states3), 2
        )
        assert_allclose(
            ancestral_state_pair(self.aln3, self.t3_2, 0, 0, self.ancestral_states3), 2
        )
        assert_allclose(
            ancestral_state_pair(self.aln3, self.t3_2, 1, 1, self.ancestral_states3), 5
        )

    def test_ancestral_state_alignment_aln_difference(self):
        """ancestral_state_alignment: difference aln -> different result"""
        expected = [[0, 0, 0], [0, 2, 0], [0, 0, 7.8]]
        actual = ancestral_state_alignment(self.aln1_1, self.t1, self.ancestral_states1)
        assert_allclose(actual, expected)

        expected = [[5, 5, 5], [5, 11.6, 11.6], [5, 11.6, 11.6]]
        actual = ancestral_state_alignment(self.aln1_5, self.t1, self.ancestral_states1)
        assert_allclose(actual, expected)

    def test_ancestral_state_position_aln_difference(self):
        """ancestral_state_position: difference aln -> different result"""

        expected = [0, 0, 0]
        actual = ancestral_state_position(
            self.aln1_1, self.t1, 0, self.ancestral_states1
        )
        assert_allclose(actual, expected)
        expected = [0, 2, 0]
        actual = ancestral_state_position(
            self.aln1_1, self.t1, 1, self.ancestral_states1
        )
        assert_allclose(actual, expected)
        expected = [0, 0, 7.8]
        actual = ancestral_state_position(
            self.aln1_1, self.t1, 2, self.ancestral_states1
        )
        assert_allclose(actual, expected)

        expected = [5, 5, 5]
        actual = ancestral_state_position(
            self.aln1_5, self.t1, 0, self.ancestral_states1
        )
        assert_allclose(actual, expected)
        expected = [5, 11.6, 11.6]
        actual = ancestral_state_position(
            self.aln1_5, self.t1, 1, self.ancestral_states1
        )
        assert_allclose(actual, expected)
        expected = [5, 11.6, 11.6]
        actual = ancestral_state_position(
            self.aln1_5, self.t1, 2, self.ancestral_states1
        )
        assert_allclose(actual, expected)

    def test_ancestral_state_pair_aln_difference(self):
        """acestral_state_pair: different aln -> different result"""
        assert_allclose(
            ancestral_state_pair(self.aln1_1, self.t1, 0, 0, self.ancestral_states1), 0
        )
        assert_allclose(
            ancestral_state_pair(self.aln1_1, self.t1, 1, 1, self.ancestral_states1), 2
        )
        assert_allclose(
            ancestral_state_pair(self.aln1_1, self.t1, 2, 2, self.ancestral_states1),
            7.8,
        )

        assert_allclose(
            ancestral_state_pair(self.aln1_5, self.t1, 0, 1, self.ancestral_states1), 5
        )
        assert_allclose(
            ancestral_state_pair(self.aln1_5, self.t1, 0, 2, self.ancestral_states1), 5
        )
        assert_allclose(
            ancestral_state_pair(self.aln1_5, self.t1, 1, 2, self.ancestral_states1),
            11.6,
        )

    def test_ancestral_state_pair_symmetry(self):
        """ancestral_state_pair: value[i,j] == value[j,i]"""
        assert_allclose(
            ancestral_state_pair(self.aln1_5, self.t1, 0, 1, self.ancestral_states1),
            ancestral_state_pair(self.aln1_5, self.t1, 1, 0, self.ancestral_states1),
        )
        assert_allclose(
            ancestral_state_pair(self.aln1_5, self.t1, 0, 2, self.ancestral_states1),
            ancestral_state_pair(self.aln1_5, self.t1, 2, 0, self.ancestral_states1),
        )
        assert_allclose(
            ancestral_state_pair(self.aln1_5, self.t1, 1, 2, self.ancestral_states1),
            ancestral_state_pair(self.aln1_5, self.t1, 2, 1, self.ancestral_states1),
        )

    def est_ancestral_state_methods_handle_alt_null_value(self):
        """ancetral state methods handle non-default null value"""
        # need to rewrite a test of this -- right now there's no way to get
        # null values into the ancestral states result, but that will change
        # when I fix the exclude handling
        pass


# following are support funcs for ResampledMiTests


def make_freqs(c12):
    c1, c2 = CategoryCounter(), CategoryCounter()
    for a, b in c12.expand():
        c1[a] += 1
        c2[b] += 1
    return c1, c2


def make_sample(freqs):
    d = []
    for i, s in enumerate(freqs.expand()):
        d += [("s%d" % i, s)]
    return make_aligned_seqs(data=d)


def _calc_mi():
    """one mutual info hand calc"""
    from math import log

    i = 37 / 42 * -log(37 / 42, 2) - (5 / 42 * log(5 / 42, 2))
    j = 39 / 42 * -log(39 / 42, 2) - (3 / 42 * log(3 / 42, 2))
    k = (
        34 / 42 * -log(34 / 42, 2)
        - (3 / 42 * log(3 / 42, 2))
        - (5 / 42 * log(5 / 42, 2))
    )
    return i + j - k


class ResampledMiTests(TestCase):
    def setUp(self):
        self.c12 = CategoryCounter(["AA", "AA", "BB", "BB", "BC"])
        self.c1, self.c2 = make_freqs(self.c12)
        self.aln = make_sample(self.c12)

    def test_calc_weights(self):
        """resampled mi weights should be correctly computed"""
        w1 = make_weights(self.c1, 5)
        w2 = make_weights(self.c2, 5)
        e = [
            ("A", {"C": 0.033333333333333333, "B": 0.066666666666666666}),
            ("B", {"A": 0.066666666666666666, "C": 0.033333333333333333}),
            ("C", {"A": 0.050000000000000003, "B": 0.050000000000000003}),
        ]

        weights = []
        for w in w1, w2:
            for k, d in w:
                weights += list(d.values())
        assert_allclose(sum(weights), 0.5)
        w2.sort()
        self.assertEqual(w2, e)

    def test_scaled_mi(self):
        """resampled mi should match hand calc"""

        def calc_scaled(data, expected_smi):
            col_i, col_j = CategoryCounter(), CategoryCounter()
            for i, j in data:
                col_i[i] += 1
                col_j[j] += 1
            pair_freqs = CategoryCounter(data)
            weights_i = make_weights(col_i, col_i.sum)
            weights_j = make_weights(col_j, col_j.sum)
            entropy = mi(col_i.entropy, col_j.entropy, pair_freqs.entropy)
            assert_allclose(entropy, _calc_mi())
            scales = calc_pair_scale(data, col_i, col_j, weights_i, weights_j)
            scaled_mi = 1 - sum(
                [w * pair_freqs[pr] for pr, e, w in scales if entropy <= e]
            )
            assert_allclose(scaled_mi, expected_smi)

        data = [
            "BN",
            "BN",
            "BP",
            "BN",
            "PN",
            "BN",
            "BN",
            "BN",
            "BN",
            "BN",
            "BN",
            "BN",
            "BN",
            "PN",
            "BN",
            "PN",
            "BN",
            "BN",
            "BN",
            "BN",
            "BN",
            "BP",
            "BN",
            "BN",
            "BN",
            "BN",
            "BP",
            "BN",
            "BN",
            "BN",
            "BN",
            "PN",
            "PN",
            "BN",
            "BN",
            "BN",
            "BN",
            "BN",
            "BN",
            "BN",
            "BN",
            "BN",
        ]
        calc_scaled(data, 8 / 42)

    def test_resampled_mi_interface(self):
        """resampled_mi_alignment should correctly compute statistic from
        alignment"""
        arr = resampled_mi_alignment(self.aln)
        # expected value from hand calculation
        assert_allclose(arr.tolist(), [[1.0, 0.78333333], [0.78333333, 1.0]])


ALN_FILE = """Seq_1   ACDEFG
Seq_2   STVWY-
Seq_3 WY.ZBX"""

# J in here
ALN_FILE_WRONG_KEY = """Seq_1   ACDKLM
Seq_2   JINCK-
Seq_3 VX.MAB"""

# last seq too long
ALN_FILE_INC_SHAPE = """Seq_1   ACDKLM
Seq_2   LINCK-
Seq_3 VX.MABN"""


gpcr_ungapped = """>OPSD_SPAAU
MNGTEGPFFYVPMVNTSGIVRSPYEYPQYYLVNPAAYAALGAYMFLLILVGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVLGRLGCNIEGFFATLGGEIALWSLVVLAIERWVVVCKPISNFRFGENHAIMGLAFTWIMAMACAAPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFICHFSIPLTIVFFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVIMMVIAFLVCWLPYAGVAWWIFTHQGSEFGPVFMTIPAFFAKSSSIYNPMIYICLNKQFRHCMITTLCCGKNPFEEEEGASTSVSSSSVSPAA-
>B1AR_CANFA
LPDGAATAARLLVPASPSASPLAPTSEGPAPLSQQWTAGIGLLMALIVLLIVAGNVLVIAAIAKTPRLQTLTNLFIMSLASADLVMGLLVVPFGATIVMRGRWEYGSFLCELWTSVDVLCVTASIETLCVIALDRYLAITAPFYQSLLTRARARALVCTVWAISALVSFLPILMHWWRAR-R--KCCDFV--------TNRAYAIASSVVSFYVPLCIMAFVYLRVFREAQKQNGRRRPSRLVALREQKALKTLGIIMGVFTLCWLPFFLANVVKAFHRDL-VPDRLFVFFNWLGYANSAFNPIIYCRSP-DFRRAFQRLLCCARRAARGSHGAAG------PPPSPG
>OPSB_HUMAN
---MRKMSEEEFYLFKNISSVGPWDGPQYHIAPVWAFYLQAAFMGTVFLIGFPLNAMVLVATLRYKKLRQPLNYILVNVSFGGFLLCIFSVFPVFVASCNGYFVFGRHVCALEGFLGTVAGLVTGWSLAFLAFERYIVICKPFGNFRFSSKHALTVVLATWTIGIGVSIPPFFGWSRFIPEGLQCSCGPDWYTVGTKYRSESYTWFLFIFCFIVPLSLICFSYTQLLRALKAVAAQQQESATTQKAEREVSRMVVVMVGSFCVCYVPYAAFAMYMVNNRNHGLDLRLVTIPSFFSKSACIYNPIIYCFMNKQFQACIMKMVCGKAMTD---ESDTCSTVSSTQVGPN-
>OPSD_PROCL
NP-Y-GNFTVVDMAPKDILHMIHPHWYQYPPMNPMMYPLLLIFMLFTGILCLAGNFVTIWVFMNTKSLRTPANLLVVNLAMSDFLMMFTMFPPMMVTCYYHTWTLGPTFCQVYAFLGNLCGCASIWTMVFITFDRYNVIVKGVAGEPLSTKKASLWILTIWVLSITWCIAPFFGWNRYVPEGNLTGCGTD--YLSEDILSRSYLYDYSTWVYYLPLLP-IYCYVSIIKAVAAHMGIRNEEAQKTSAECRLAKIAMTTVALWFIAWTPYLLINWVGMFARSY-LSPVYTIWGYVFAKANAVYNPIVYAISHPKYRAAMEKKLPCLSCKTESDDVSESAEEKAESA----
>BRB1_RABIT
ASQGPLELQPSNQSQLAPPNATSC--SGAPDAWDLLHRLLPTFIIAIFTLGLLGNSFVLSVFLLARRRLSVAEIYLANLAASDLVFVLGLPFWAENVRNQFDWPFGAALCRIVNGVIKANLFISIFLVVAISQDRYSVLVHPMSRRGRRRRQAQATCALIWLAGGLLSTPTFVLRSVRA-LN--SACILL---LPHEAWHWLRMVELNLLGFLLPLAAILFFNCHILASLRRR---RVPSRCGGPRDSKSTALILTLVASFLVCWAPYHFFAFLECLWQVHEFTDLGLQLSNFSAFVNSCLNPVIYVFVGRLFRTKVWELCQQCSPR---------LAPV-------S
>5H7_.ENLA
NLLPSEFMTERPLNTTEQDLTKPDCGKELLLYGDTEKIVIGVVLSIITLFTIAGNALVIISVCIVKKLRQPSNYLVVSLAAADLSVAVAVMPFVIITDLVGGWLFGKVFCNVFIAMDVMCCTASIMTLCVISVDRYLGITRPLYPARQNGKLMAKMVFIVWLLSASITLPPLFGWAK--N-V--RVCLIS--------QDFGYTVYSTAVAFYIPMTVMLVMYQRIFVAAKISSKLDRKNISIFKREQKAARTLGIIVGAFTFCWLPFFLLSTARPFICGICMPLRLERTLLWLGYTNSLINPLIYAFFNRDLRTTFWNLLRCKYTNINRRLSAASTERHEGIL----
>B3AR_FELCA
MAPWPHGNGSLASWPDAPTLTPNTANTSGLPGVPWAVALAGALLALAVLATVGGNLLVIVAIARTPRLQTMTNVFVTSLATADLVVGLLVVPPGATLALTGHWPLGATGCELWTSVDVLCVTASIETLCALAVDRYLAVTNPLYGALVTKRRARAAVVLVWVVSAAVSFAPIMSKWWRVQ-R--HCCAFA--------SNIPYALLSSSVSFYLPLLVMLFVYARVFVVATRQGVPRRPARLLPLREHRALRTLGLIMGTFSLCWLPFFVANVVRALGGPS-VPSPAFLALNWLGYANSAFNPLIYCRSP-DFRSAFRRLLCRCRLEERHAAASGAAALTRPAESGLP
>TLR1_DROME
IIDNRDNLESINEAKDFLTECLFPSPTRPYELPWEQKTIWAIIFGLMMFVAIAGNGIVLWIVTGHRSMRTVTNYFLLNLSIADLLMSSLNCVFNFIFMLNSDWPFGSIYCTINNFVANVTVSTSVFTLVAISFDRYIAIVDPL-KRRTSRRKVRIILVLIWALSCVLSAPCLLYSS----SR--TVCFMMDGRYPTSMADYAYNLIILVLTTGIPMIVMLICYSLMGRVPGGSSIGTDRQMESMKSKRKVVRMFIAIVSIFAICWLPYHLFFIYAYHNNQVKYVQHMYLGFYWLAMSNAMVNPLIYYWMNKRFRMYFQRIICCCCVGLTRHRFDSPNSSNRHTRAETK
>ITR_CATCO
EQDFWSFNESSRNSTVGNETF-GNQTVNPLKRNEEVAKVEVTVLALVLFLALAGNLCVLIAIYTAKHTQSRMYYLMKHLSIADLVVAVFQVLPQLIWDITFRFYGPDFLCRLVKYLQTVGMFASTYMLVLMSIDRCIAICQPL--RSLHKRKDRCYVIVSWALSLVFSVPQVYIFSLRE----VYDCWGD---FVQPWGAKAYITWISLTIYIIPVAILGGCYGLISFKIWQNANAVSSVKLVSKAKITTVKMTFVIVLAYIVCWTPFFFVQMWSAWDPEA-REAMPFIISMLLASLNSCCNPWIYMFFAGHLFHDLKQSLLCCSTLYLKSSQCRCKSNCSTYVIKST
>NK1R_RANCA
----MNSNISAQNDSALNSTIQNGTKINQFIQPPWQIALWSVAYSIIVIVSLVGNIIVMWIIIAHKRMRTVTNYFLVNLAFAEASMSAFNTVINFTYAIHNHWYYGLIYCKFHNFFPISAVFTSIYSMTAIALDRYMAIIHPL-KPRLSATATKIVICVIWSFSFCMAFPLGYYAD---GG---DICYLNPDSEENRKYEQVYQVLVFCLIYILPLLVIGCAYTFIGMTLWAS--PSDRYHEQVVAKRKVVKMMIVVVCTFAICWLPFHIFFLLQTLHEMTKFYQQFYLAIMWLAMSSTMYNPIIYCCLNDRFRIGFKHVFRWCPFIRAGEY----STRYLQTQSSMY
>A2AA_MOUSE
----MGSLQPDAGNSSWNGTEAPGGGTRATPYSLQVTLTLVCLAGLLMLFTVFGNVLVIIAVFTSRALKAPQNLFLVSLASADILVATLVIPFSLANEVMGYWYFGKVWCEIYLALDVLFCTSSIVHLCAISLDRYWSITQAIYNLKRTPRRIKAIIVTVWVISAVISFPPLISIEKKGQ-P--PSCKIN--------DQKWYVISSSIGSFFAPCLIMILVYVRIYQIAKRRRGGASRWRGRQNREKRFTFVLAVVIGVFVVCWFPFFFTYTLIAVGCP--VPSQLFNFFFWFGYCNSSLNPVIYTIFNHDFRRAFKKILCRGDRKRIV------------------
>MSHR_VULVU
SGQGPQRRLLGSPNATSPTTPHFKLAANQTGPRCLEVSIPNGLFLSLGLVSVVENVLVVAAIAKNRNLHSPMYYFIGCLAVSDLLVSVTNVLETAVMLLVEAAAVVQQLDDIIDVLICGSMVSSLCFLGAIAVDRYLSIFYALYHSIVTLPRAWRAISAIWVASVLSSTLFIAYYY----------------------NNHTAVLLCLVSFFVAMLVLMAVLYVHMLARARQHIARKRQHSVHQGFGLKGAATLTILLGIFFLCWGPFFLHLSLMVLCPQHGCVFQNFNLFLTLIICNSIIDPFIYAFRSQELRKTLQEVVLCSW-----------------------
>OPSD_BOVIN
MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA
>OPSG_MOUSE
DHYEDSTHASIFTYTNSNSTKGPFEGPNYHIAPRWVYHLTSTWMILVVVASVFTNGLVLAATMRFKKLRHPLNWILVNLAVADLAETIIASTISVVNQIYGYFVLGHPLCVIEGYIVSLCGITGLWSLAIISWERWLVVCKPFGNVRFDAKLATVGIVFSWVWAAIWTAPPIFGWSRYWPYGLKTSCGPDVFSGTSYPGVQSYMMVLMVTCCIFPLSIIVLCYLQVWLAIRAVAKQQKESESTQKAEKEVTRMVVVMVFAYCLCWGPYTFFACFATAHPGYAFHPLVASLPSYFAKSATIYNPIIYVFMNRQFRNCILHLFGKKVDDS-----SELVSSV--SSVSPA
>TA2R_RAT
--MWLNS-----------TSLGACFRPVNITLQERRAIASPWFAASFCALGLGSNLLALSVLAGARPPRSSFLALLCGLVLTDFLGLLVTGAVVASQHAALLTDPGCRLCHFMGAAMVFFGLCPLLLGAAMAAERFVGITRPFSRPAATSRRAWATVGLVWVGAGTLGLLPLLGLGRYSVQYPGSWCFLT----LGAERGDVAFGLMFALLGSVSVGLSLLLNTVSVATLCRVYHAREATQRPRDCEVEMMVQLVGIMVVATVCWMPLLVFILQTLLQTLPRTTERQLLIYLRVATWNQILDPWVYILFRRSVLRRLHPRFTSQLQAVSLHSPPTQ------------
>CCR3_HUMAN
LENFSSSYDYGENESDSCCTSPPC---PQDFSLNFDRAFLPALYSLLFLLGLLGNGAVAAVLLSRRTALSSTDTFLLHLAVADTLLVLTLPLWAVDAAVQ--WVFGSGLCKVAGALFNINFYAGALLLACISFDRYLNIVHATLYRRGPPARVTLTCLAVWGLCLLFALPDFIFLSAHH-RL-ATHCQYN----FPQVGRTALRVLQLVAGFLLPLLVMAYCYAHILAVL---------LVSRGQRRLRAMRLVVVVVVAFALCWTPYHLVVLVDILMDLGSRVDVAKSVTSGLGYMHCCLNPLLYAFVGVKFRERMWMLLLRLGCPN----QRGL--------PSSS
>THRR_CRILO
LPEGRAIYLNKSHSPPAPLAPFISEDASGYLTSPWLRLFIPSVYTFVFVVSLPLNILAIAVFVLKMKVKKPAVVYMLHLAMADVLFVSVLPLKISYYFSGSDWQFGSGMCRFATAAFYCNMYASIMLMTVISIDRFLAVVYPISLSWRTLGRANFTCLVIWVMAIMGVVPLLLKEQTTR--N--TTCHDVLNETLLQGFYSYYFSAFSAVFFLVPLIISTICYMSIIRCL------SSSSVANRSKKSRALFLSAAVFCVFIVCFGPTNVLLIMHYLLLSDEKAYFAYLLCVCVSSVSCCIDPLIYYYASSECQRHLYGILCCKESSDPNSYNSTGDTCS--------
>AA2A_CAVPO
----------------------------------MSSSVYITVELVIAVLAILGNVLVCWAVWINSNLQNVTNYFVVSLAAADIAVGVLAIPFAITISTG--FCAACHGCLFFACFVLVLTQSSIFSLLTITIDRYIAIRIPLYNGLVTCTRAKGIIAICWVLSFAIGLTPMLGWNNCSS-E-QVTCLFE-----DVVPMNYMVYYNFFAFVLVPLLLMLGIYLRIFLAARRQESQGERTRSTLQKEVHPAKSLAIIVGLFALCCLPLNIINCFTFFCPECHAPPWLMYLTIILSHGNSVVNPLIYAYRIREFRQTFRKIIRSHILRRRELFKAGGAHSPEGEQVSLR
>C3AR_MOUSE
--------------MESFDADTNSTDLHSRPLFQPQDIASMVILGLTCLLGLLGNGLVLWVAGVKMK-TTVNTVWFLHLTLADFLCCLSLPFSLAHLILQGHWPYGLFLCKLIPSIIILNMFASVFLLTAISLDRCLIVHKPICQNHRNVRTAFAICGCVWVVAFVMCVPVFVYRDLFI-ED-DYVDQFT-YDNHVPTPLMAITITRLVVGFLVPFFIMVICYSLIVFRM--------RKTNFTKSRNKTFRVAVAVVTVFFICWTPYHLVGVLLLITDPEEAVMSWDHMSIALASANSCFNPFLYALLGKDFRKKARQSIKGILEAAFSEELTHSASS---------
>TA2R_BOVIN
--MWPNA-----------SSLGPCFRPMNITLEERRLIASPWFAASFCLVGLASNLLALSVLMGARQSRSSFLTFLCGLVLTDFMGLLVTGAIVVTQHFVLFVDPGCSLCHFMGVIMVFFGLCPLLLGAAMASERFLGITRPFRPATASQRRAWTTVGLVWASALALGLLPLLGVGHYTVQYPGSWCFLT----LGTDPGDVAFGLLFALLGSISVGMSFLLNTISVATLCHVYHGATAQQRPRDCEVEMMVQLMGIMVVASICWMPLLVFIAQTVLQSPPRLTERQLLIYLRVATWNQILDPWVYILFRRAVIQRFYPRLSTRSRSLSLQPQLTR------------
>OPSD_RABIT
MNGTEGPDFYIPMSNQTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTTTLYTSLHGYFVFGPTGCNVEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWIMALACAAPPLVGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPLIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSSSIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATASKTETSQVAPA
>OPSD_PHOVI
MNGTEGPNFYVPFSNKTGVVRSPFEFPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVGFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTLPAFFAKAASIYNPVIYIMMNKQFRTCMITTLCCGKNPLGDDEVSASASKTETSQVAPA
>D3DR_RAT
--------MAPLSQISTHLNSTCGAENSTGVNRARPHAYYALSYCALILAIIFGNGLVCAAVLRERALQTTTNYLVVSLAVADLLVATLVMPWVVYLEVTGGWNFSRICCDVFVTLDVMMCTASILNLCAISIDRYTAVVMPVGTGQSSCRRVALMITAVWVLAFAVSCPLLFGFN---D-P--SICSI---------SNPDFVIYSSVVSFYVPFGVTVLVYARIYIVLRQRTSLPLQPRGVPLREKKATQMVVIVLGAFIVCWLPFFLTHVLNTHCQACHVSPELYRATTWLGYVNSALNPVIYTTFNVEFRKAFLKILSC-------------------------
>DADR_.ENLA
---------------MTFNITSMDEDVLLTERESSFRVLTGCFLSVLILSTLLGNTLVCAAVIRFRHRSKVTNFFVISLAVSDLLVAVLVMPWKAVAEIAGFWPFG-TFCNIWVAFDIMCSTASILNLCVISVDRYWAISSPFYERKMTPKVAFIMIGVAWTLSVLISFIPVQLNWHKALN-TMDNCDSS--------LNRTYAISSSLISFYIPVAIMIVTYTRIYRIAAKQLDCESSLKTSFKRETKVLKTLSVIMGVFVCCWLPFFILNCIVPFCDPSCISSTTFDVFVWFGWANSSLNPIIYAFNA-DFRKAFSNLLGCYRLCPTSNNIIETAVVYSCQ-----
>PAR2_HUMAN
---VDGTSHVTGKGVTVETVFSVDEFSASVLTGKLTTVFLPIVYTIVFVVGLPSNGMALWVFLFRTKKKHPAVIYMANLALADLLSVIWFPLKIAYHIHGNNWIYGEALCNVLIGFFYGNMYCSILFMTCLSVQRYWVIVNPMGHSRKKANIAIGISLAIWLLILLVTIPLYVVKQTIF--N--TTCHDVLPEQLLVGDMFNYFLSLAIGVFLFPAFLTASAYVLMIRML------SAMDENSEKKRKRAIKLIVTVLAMYLICFTPSNLLLVVHYFLIKSSHVYALYIVALCLSTLNSCIDPFVYYFVSHDFRDHAKNALLCRSVRTVKQMQVSLKSSS--------
>GPRA_RAT
TPANQSAEASESNVSATVPRAAAVTPFQSLQLVHQLKGLIVMLYSIVVVVGLVGNCLLVLVIARVRRLHNVTNFLIGNLALSDVLMCAACVPLTLAYAFEPRWVFGGGLCHLVFFLQPVTVYVSVFTLTTIAVDRYVVLVHPL-RRRISLKLSAYAVLGIWALSAVLALPAAVHTYHVEHD--VRLCEEF--WGSQERQRQIYAWGLLLGTYLLPLLAILLSYVRVSVKLRNRPGS-SQADWDRARRRRTFCLLVVVVVVFALCWLPLHIFNLLRDLDPRAYAFGLVQLLCHWLAMSSACYNPFIYAWLHDSFREELRKMLLSWPRKIVPHGQNMT------------
>ET1R_PIG
EFSLVVTTHRPTNLALPSNGSMHNYCPQQTKITSAFKYINTVISCTIFIVGMVGNATLLRIIYQNKCMRNGPNALIASLALGDLIYVVIDLPINVFKLLAGRHDFGVFLCKLFPFLQKSSVGITVLNLCALSVDRYRAVASWSVQGIGIPLVTAIEIVSIWILSFILAIPEAIGFVMVPKT--HKTCMLNATSKFMEFYQDVKDWWLFGFYFCMPLVCTAIFYTLMTCEMLNRGSL-IALSEHLKQRREVAKTVFCLVVIFALCWFPLHLSRILKKTVYDESFLLLMDYIGINLATMNSCINPIALYFVSKKFKNCFQSCLCCCCYQSKSLMTSVPWKNHEQNNHNTE
>NTR2_MOUSE
-----METSSLWPPRPSPSAGLSLEARLGVDTRLWAKVLFTALYSLIFALGTAGNALSVHVVLKARTRPGRLRYHVLSLALSALLLLLISVPMELYNFVWSHWVFGDLGCRGYYFVRELCAYATVLSVASLSAERCLAVCQPLARRLLTPRRTCRLLSLVWVASLGLALPMAVIMGQKH--AASRVCTVL---VSRASSRSTFQVKRAGLLRSPLWELTAILNGITVNHLVALVQARHKDASQIRSLQHSAQVLRAIVAVYVICWLPYHARRLMYCYIPDDDFYHYFYMVTNTLFYVSSAVTPVLYNAVSSSFRKLFLESLSSLCGEQRSVVPLPQSTYSFRLWGSPR
>CKR7_MOUSE
QDEVTDDYIGENTTVDYTLYESVC---FKKDVRNFKAWFLPLMYSVICFVGLLGNGLVILTYIYFKRLKTMTDTYLLNLAVADILFLLILPFWAYSEAKS--WIFGVYLCKGIFGIYKLSFFSGMLLLLCISIDRYVAIVQAVRHRARVLLISKLSCVGIWMLALFLSIPELLYSGLQK-GE-TLRCSLV---SAQVEALITIQVAQMVFGFLVPMLAMSFCYLIIIRTL---------LQARNFERNKAIKVIIAVVVVFIVFQLPYNGVVLAQTVANFNKQLNIAYDVTYSLASVRCCVNPFLYAFIGVKFRSDLFKLFKDLGCLSQERLRHWS--------HVRN
>NY5R_HUMAN
-------------NTAATRNSDFPVWDDYKSSVDDLQYFLIGLYTFVSLLGFMGNLLILMALMKKRNQKTTVNFLIGNLAFSDILVVLFCSPFTLTSVLLDQWMFGKVMCHIMPFLQCVSVLVSTLILISIAIVRYHMIKHPI-SNNLTANHGYFLIATVWTLGFAICSPLPVFHSLVESS--RYLCVES---WPSDSYRIAFTISLLLVQYILPLVCLTVSHTSVCRSISCGVHEKRSVTRIKKRSRSVFYRLTILILVFAVSWMPLHLFHVVTDFNDNLRHFKLVYCICHLLGMMSCCLNPILYGFLNNGIKADLVSLIHCLHM----------------------
>VG74_HSVSA
VKLDFSSEDFSNYSYNYSGDIYYGDVAPCVVNFLISESALAFIYVLMFLCNAIGNSLVLRTFLKYRA-QAQSFDYLMMGFCLNSLFLAGYLLMRLLRMFE--IFMNTELCKLEAFFLNLSIYWSPFILVFISVLRCLLIFCATRLWVKKTLIGQVFLCCSFVLACFGALPHVMVTSYYE----PSSCIEE---VLTEQLRTKLNTFHTWYSFAGPLFITVICYSMSCYKL---------FKTKLSKRAEVVTIITMTTLLFIVFCIPYYIMESIDTLLRVGSAIVYGIQCTYMLLVLYYCMLPLMFAMFGSLFRQRMAAWCKTICHC---------------------
>CCKR_HUMAN
DSLLVNGSNITPPCELGLENETLFCLDQPRPSKEWQPAVQILLYSLIFLLSVLGNTLVITVLIRNKRMRTVTNIFLLSLAVSDLMLCLFCMPFNLIPNLLKDFIFGSAVCKTTTYFMGTSVSVSTFNLVAISLERYGAICKPLSRVWQTKSHALKVIAATWCLSFTIMTPYPIYN----NNQTANMCRFL---LPNDVMQQSWHTFLLLILFLIPGIVMMVAYGLISLELYQGRANSNSSAANLMAKKRVIRMLIVIVVLFFLCWMPIFSANAWRAYDTASRLSGTPISFILLLSYTSSCVNPIIYCFMNKRFRLGFMATFPCCPNPGPPGARGEVTTGASLSRFSYS
>OPSR_CARAU
GD--ETTRESMFVYTNSNNTRDPFEGPNYHIAPRWVYNLATVWMFFVVVASTFTNGLVLVATAKFKKLRHPLNWILVNLAVADLAETLLASTISVTNQFFGYFILGHPMCIFEGFTVSVCGIAGLWSLTVISWERWVVVCKPFGNVKFDAKWASAGIIFSWVWSAIWCAPPIFGWSRFWPHGLKTSCGPDVFSGSEDPGVQSYMIVLMITCCIIPLAIIILCYIAVWLAIRTVAQQQKDSESTQKAEKEVSRMVVVMIFAYCFCWGPYTFCACFAAANPGYAFHPLAAAMPAYFAKSATIYNPIIYVFMNRQFRVCIMQLFGKKVDDG-----SEVSS------VAPA
>OPS2_LIMPO
PN-----ASVVDTMPKEMLYMIHEHWYAFPPMNPLWYSILGVAMIILGIICVLGNGMVIYLMMTTKSLRTPTNLLVVNLAFSDFCMMAFMMPTMASNCFAETWILGPFMCEVYGMAGSLFGCASIWSMVMITLDRYNVIVRGMAAAPLTHKKATLLLLFVWIWSGGWTILPFFGWSRYVPEGNLTSCTVD--YLTKDWSSASYVIIYGLAVYFLPLITMIYCYFFIVHAVAEHNVAANADQQKQSAECRLAKVAMMTVGLWFMAWTPYLIIAWAGVFSSGTRLTPLATIWGSVFAKANSCYNPIVYGISHPRYKAALYQRFPSLACGSGESGSDVKTMEEKPKSPEA-
>O5I1_HUMAN
-----------MEFTDRNYTLVTEFILLGFPTRPELQIVLFLMFLTLYAIILIGNIGLMLLIRIDPHLQTPMYFFLSNLSFVDLCYFSDIVPKMLVNFLSENKSISYYGCALQFYFFCTFADTESFILAAMAYDRYVAICNPLYTVVMSRGICMRLIVLSYLGGNMSSLVHTSFAFIL---KNHFFCDLPKLSCTDTTINEWLLSTYGSSVEIICFIIIIISYFFILLSV--------LKIRSFSGRKKTFSTCASHLTSVTIYQGTLLFIYSRPSYLY---SPNTDKIISVFYTIFIPVLNPLIYSLRNKDVKDAAEKVLRSKVDSS--------------------
>OAR_HELVI
--TEEVIEDDRDACAVADDPKYPSSFGITLAVPEWEAICTAIVLTLIIISTIVGNILVILSVFTYKPLRIVQNFFIVSLAVADLTVAILVLPLNVAYSILGQWVFGIYVCKMWLTCDIMCCTSSILNLCAIALDRYWAITDPIYAQKRTLERVLLMIGVVWVLSLIISSPPLLGWNDW-E-P--TPCRLT--------SQPGFVIFSSSGSFYIPLVIMTVVYFEIYLATKKRAVYEEKQRISLTRERRAARTLGIIMGVFVVCWLPFFVIYLVIPFCASCCLSNKFINFITWLGYCNSALNPLIYTIFNMDFRRAFKKLLCMKP-----------------------
>O1D4_HUMAN
-------------MDGDNQSENSQFLLLGISESPEQQQILFWMFLSMYLVTVLGNVLIILAISSDSHLHTPMYFFLANLSFTDLFFVTNTIPKMLVNFQSQNKAISYAGCLTQLYFLVSLVTLDNLILAVMAYDRYVAICCPLYVTAMSPGLCVLLLSLCWGLSVLYGLLLTFLLTRV---THYLFCDMYWLACSNTHIIHTALIATGWFIFLTLLGFMTTSYVRIVRTI--------LQMPSASKKYKTFSTCASHLGVVSLFYGTLAMVYLQPLHTY----SMKDSVATVMYAVLTPMMNPFIYSLRNKDMHGAPGRVLWRPFQRP--------------------
>OPSD_APIME
AR-F-NNQTVVDKVPPDMLHLIDANWYQYPPLNPMWHGILGFVIGMLGFVSAMGNGMVVYIFLSTKSLRTPSNLFVINLAISNFLMMFCMSPPMVINCYYETWVLGPLFCQIYAMLGSLFGCGSIWTMTMIAFDRYNVIVKGLSGKPLSINGALIRIIAIWLFSLGWTIAPMFGWNRYVPEGNMTACGTD--YFNRGLLSASYLVCYGIWVYFVPLFLIIYSYWFIIQAVAAHMNVRSSENQNTSAECKLAKVALMTISLWFMAWTPYLVINFSGIFNLVK-ISPLFTIWGSLFAKANAVYNPIVYGISHPKYRAALFAKFPSLACAAEPSSDAVSVTDNEKSNA---
>TRFR_CHICK
----------MENGTGDEQNHTGLLLSSQEFVTAEYQVVTILLVLLICGLGIVGNIMVVLVVLRTKHMRTPTNCYLVSLAVADLMVLVAAGLPNITESLYKSWVYGYVGCLCITYLQYLGINASSFSITAFTIERYIAICHPIAQFLCTFSRAKKIIIFVWSFASVYCMLWFFLLDLN--DT-VVSCGYK------RSYYSPIYMMDFGIFYVLPMVLATVLYGLIARILFLNVNSNKSFNSTIASRRQVTKMLAVVVVLFAFLWMPYRTLVVVNSFLSSPFQENWFLLFCRICIYLNSAINPVIYNLMSQKFRAAFRKLCNCHLKRDKKPANYSVKESDHFSSEIED
>OPSR_HUMAN
DSYEDSTQSSIFTYTNSNSTRGPFEGPNYHIAPRWVYHLTSVWMIFVVTASVFTNGLVLAATMKFKKLRHPLNWILVNLAVADLAETVIASTISIVNQVSGYFVLGHPMCVLEGYTVSLCGITGLWSLAIISWERWLVVCKPFGNVRFDAKLAIVGIAFSWIWSAVWTAPPIFGWSRYWPHGLKTSCGPDVFSGSSYPGVQSYMIVLMVTCCIIPLAIIMLCYLQVWLAIRAVAKQQKESESTQKAEKEVTRMVVVMIFAYCVCWGPYTFFACFAAANPGYAFHPLMAALPAYFAKSATIYNPVIYVFMNRQFRNCILQLFGKKVDDG-----SELVSSV--SSVSPA
>CB2R_MOUSE
----MEGCRETEVTNGSNGGLEFNPMKEYMILSSGQQIAVAVLCTLMGLLSALENMAVLYIILSSRRRRKPSYLFISSLAGADFLASVIFACNFVIFHVFHG-VDSNAIFLLKIGSVTMTFTASVGSLLLTAVDRYLCLCYPPYKALVTRGRALVALCVMWVLSALISYLPLMGWTC-----CPSPCSEL------FPLIPNDYLLGWLLFIAILFSGIIYTYGYVLWKAHRHAEHQVPGIARMRLDVRLAKTLGLVLAVLLICWFPALALMGHSLVTTLSDQVKEAFAFCSMLCLVNSMVNPIIYALRSGEIRSAAQHCLIGWKKYLQGLGPEGKVTETEADVKTT-
>PE21_HUMAN
--MSPCGPLNLSLAGEATTCAA---PWVPNTSAVPPSGASPALPIFSMTLGAVSNLLALALLAQAAGSATTFLLFVASLLATDLAGHVIPGALVLRLYTAG-RAPAGGACHFLGGCMVFFGLCPLLLGCGMAVERCVGVTRPLHAARVSVARARLALAAVAAVALAVALLPLARVGRYELQYPGTWCFIG--LGPPGGWRQALLAGLFASLGLVALLAALVCNTLSGLALHRSRRRAHGPRRARAHDVEMVGQLVGIMVVSCICWSPMLVLVALAVGGWSSTSLQRPLFLAVRLASWNQILDPWVYILLRQAVLRQLLRLLPPRAGAKGGPAGLGLSSLRSSRHSGLS
>OPR._CAVPO
GSHLQGNLSLLSPNHSGLPPHLLLNASHSAFLPLGLKVTIVGLYLAVCIGGLLGNCLVMYVILRHTKMKTATNIYIFNLALADTLVLLTLPFQATDILLGF-WPFGNTLCKTVIAIDYYNMFTSTFTLTAMSVDRYVAICHPIALDVRTSSKAQAVNVAIWALALVVGVPVAIMGSAQVEE---IECLVE-IPDPQDYWGPVFAVSIFLFSFIIPVLIISVCYSLMIRRLHGVRLL-SGSREKDRNLRRITRLVLVVVAVFVGCWTPVQVFVLVQGLGVQPETTVAILRFCTALGYVNSCLNPILYAFLDENFKACFRKFCCASALHREMQVSDRVALGCKTTETVPR
>OPSV_CHICK
-----MSSDDDFYLFTNGSVPGPWDGPQYHIAPPWAFYLQTAFMGIVFAVGTPLNAVVLWVTVRYKRLRQPLNYILVNISASGFVSCVLSVFVVFVASARGYFVFGKRVCELEAFVGTHGGLVTGWSLAFLAFERYIVICKPFGNFRFSSRHALLVVVATWLIGVGVGLPPFFGWSRYMPEGLQCSCGPDWYTVGTKYRSEYYTWFLFIFCFIVPLSLIIFSYSQLLSALRAVAAQQQESATTQKAEREVSRMVVVMVGSFCLCYVPYAALAMYMVNNRDHGLDLRLVTIPAFFSKSACVYNPIIYCFMNKQFRACIMETVCGKPLTD--DSDASTSSVSSSQVGPT-
>YKR5_CAEEL
MNSENGLDSVTQIMYDMKKYNIVNDVLPPPNHEDLHVVIMAVSYLLLFLLGTCGNVAVLTTIYHVIRTLDNTLIYVIVLSCVDFGVCLSLPITVIDQILGF-WMFGKIPCKLHAVFENFGKILSALILTAMSFDRYAGVC---------HPQRKRLRSRNFAITILLAPGMLTRM-------KIEKCTVD----IDSQMFTAFTIYQFILCYCTPLVLIAFFYTKLLSKLRE----TRTFKSSQIPFLHISLYTLAVACFYFLCWTPFWMATLFAVYLENSPVFVYIMYFIHALPFTNSAINWILYGRVFLETVS---------------------------------
>5H1B_MOUSE
CAPPPPAASQTGVPLTNLSHNSADGYIYQDSIALPWKVLLVALLALITLATTLSNAFVIATVYRTRKLHTPANYLIASLAVTDLLVSILVMPISTMYTVTGRWTLGQVVCDFWLSSDITCCTASIMHLCVIALDRYWAITDAVYSAKRTPKRAAIMIVLVWVFSISISLPPFFWR----E-E--LDCFVN-------TDHVLYTVYSTVGAFYLPTLLLIALYGRIYVEARSRRVSLEKKKLMAARERKATKTLGIILGAFIVCWLPFFIISLVMPICKDAWFHMAIFDFFNWLGYLNSLINPIIYTMSNEDFKQAFHKLIRFKCAG---------------------
>P2YR_HUMAN
AAFLAGPGSSWGNSTVASTAAVSSSFKCALTKTGFQFYYLPAVYILVFIIGFLGNSVAIWMFVFHMKPWSGISVYMFNLALADFLYVLTLPALIFYYFNKTDWIFGDAMCKLQRFIFHVNLYGSILFLTCISAHRYSGVVYPLSLGRLKKKNAICISVLVWLIVVVAISPILFYSGTG----KTITCYDT-TSDEYLRSYFIYSMCTTVAMFCVPLVLILGCYGLIVRALIY------KDLDNSPLRRKSIYLVIIVLTVFAVSYIPFHVMKTMNLRARLDDRVYATYQVTRGLASLNSCVDPILYFLAGDTFRRRLSRATRKASRRSEANLQSKSLPEFKQNGDTSL
>5H4_CAVPO
------------------MDKLDANVSSKEGFGSVEKVVLLTFLSAVILMAILGNLLVMVAVCRDRQRKIKTNYFIVSLAFADLLVSVLVMPFGAIELVQDIWVYGEMFCLVRTSLDVLLTTASIFHLCCISLDRYYAICCQPYRNKMTPLRIALMLGGCWVIPMFISFLPIMQGWNNIRK-NSTYCVFM--------VNKPYAITCSVVAFYIPFLLMVLAYYRIYVTAKEHGRPDQHSTHRMRTETKAAKTLCIIMGCFCLCWAPFFVTNIVDPFIDYT-VPGQLWTAFLWLGYINSGLNPFLYAFLNKSFRRAFLIILCCDDERYRRPSILGQTINGSTHVLR--
>MC4R_PIG
GMHTSLHFWNRSTYGLHSNASEPLGKGYSEGGCYEQLFVSPEVFVTLGVISLLENILVIVAIAKNKNLHSPMYFFICSLAVADMLVSVSNGSETIVITLLNSQSFTVNIDNVIDSVICSSLLASICSLLSIAVDRYFTIFYALYHNIMTVKRVGIIISCIWAVCTVSGVLFIIYYS----------------------DDSSAVIICLITVFFTMLALMASLYVHMFLMARLH-RIPGTGTIRQGANMKGAITLTILIGVFVVCWAPFFLHLIFYISCPQNVCFMSHFNLYLILIMCNSIIDPLIYALRSQELRKTFKEIICCYPLGGLCDLSSRY------------
>SSR1_MOUSE
GEGACSRGPGSGAADGMEEPGRNASQNGTLSEGQGSAILISFIYSVVCLVGLCGNSMVIYVILRYAKMKTATNIYILNLAIADELLMLSVPFLVTSTLLRH-WPFGALLCRLVLSVDAVNMFTSIYCLTVLSVDRYVAVVHPIAARYRRPTVAKVVNLGVWVLSLLVILPIVVFSRTAADG--TVACNML-MPEPAQRWLVGFVLYTFLMGFLLPVGAICLCYVLIIAKMRMVALK-AGWQQRKRSERKITLMVMMVVMVFVICWMPFYVVQLVNVFAEQ--DDATVSQLSVILGYANSCANPILYGFLSDNFKRSFQRILCLS-----WMDNAAETALKSRAYSVED
>CML1_HUMAN
EDEDYNTSISYGDEYPDYLDSIVVLEDLSPLEARVTRIFLVVVYSIVCFLGILGNGLVIIIATFKMK-KTVNMVWFLNLAVADFLFNVFLPIHITYAAMDYHWVFGTAMCKISNFLLIHNMFTSVFLLTIISSDRCISVLLPVSQNHRSVRLAYMACMVIWVLAFFLSSPSLVFRDTAN-SS--WPTHSQ-MDPVGYSRHMVVTVTRFLCGFLVPVLIITACYLTIVCKL---------HRNRLAKTKKPFKIIVTIIITFFLCWCPYHTLNLLELHHTAMSVFSLGLPLATALAIANSCMNPILYVFMGQDFKK-FKVALFSRLVNALSEDTGHSFTKMSSMNERTS
>ACM1_RAT
------------MNTSVPPAVSPNITVLAPGKGPWQVAFIGITTGLLSLATVTGNLLVLISFKVNTELKTVNNYFLLSLACADLIIGTFSMNLYTTYLLMGHWALGTLACDLWLALDYVASNASVMNLLLISFDRYFSVTRPLYRAKRTPRRAALMIGLAWLVSFVLWAPAILFWQYLV-VL-AGQCYIQ------FLSQPIITFGTAMAAFYLPVTVMCTLYWRIYRETENRRGKAKRKTFSLVKEKKAARTLSAILLAFILTWTPYNIMVLVSTFCKDC-VPETLWELGYWLCYVNSTVNPMCYALCNKAFRDTFRLLLLCRWDKRRWRKIPKRPSRQC-------
>MRG_HUMAN
QNPNLVSQLCGVFLQNETNETIHMQMSMAVGQQALPLNIIAPKAVLVSLCGVLLNGTVFWLLCCGAT--NPYMVYILHLVAADVIYLCCSAVGFLQVTLLTYHGVVFFIPDFLAILSPFSFEVCLCLLVAISTERCVCVLFPIYRCHRPKYTSNVVCTLIWGLPFCINIVKSLFLTYWK-------------------KACVIFLKLSGLFHAILSLVMCVSSLTLLIRFL--------CCSQQQKATRVYAVVQISAPMFLLWALPLSVAPLITDF----KMFVTTSYLISLFLIINSSANPIIYFFVGSLRKKRLKESLRVILQRALADKPEVGIDPMEQPHSTQH
>A2AC_RAT
AEGPNGSDAGEWGSGGGANASGTDWGPPPGQYSAGAVAGLAAVVGFLIVFTVVGNVLVVIAVLTSRALRAPQNLFLVSLASADILVATLVMPFSLANELMAYWYFGQVWCGVYLALDVLFCTSSIVHLCAISLDRYWSVTQAVYNLKRTPRRVKATIVAVWLISAVISFPPLVSFYR-------PQCGLN--------DETWYILSSCIGSFFAPCLIMGLVYARIYRVAKLRRRAVCRRKVAQAREKRFTFVLAVVMGVFVLCWFPFFFSYSLYGICREAQLPEPLFKFFFWIGYCNSSLNPVIYTVFNQDFRRSFKHILFRRRRRGFRQ-----------------
>OPSD_SEPOF
---MGRDIPDNETWWYNPTMEVHPHWKQFNQVPDAVYYSLGIFIGICGIIGCTGNGIVIYLFTKTKSLQTPANMFIINLAFSDFTFSLVNGFPLMTISCFIKWVFGMAACKVYGFIGGIFGLMSIMTMSMISIDRYNVIGRPMASKKMSHRRAFLMIIFVWMWSTLWSIGPIFGWGAYVLEGVLCNCSFD--YITRDSATRSNIVCMYIFAFCFPILIIFFCYFNIVMAVSNHRLNLRKAQAGASAEMKLAKISIVIVTQFLLSWSPYAVVALLAQFGPIEWVTPYAAQLPVMFAKASAIHNPLIYSVSHPKFREAIAENFPWIITCCQFDEKEVEEIPATEQS-GGE
>SSR3_HUMAN
SVSTTSEPENASSAWPPDATLGNVSAGPSPAGLAVSGVLIPLVYLVVCVVGLLGNSLVIYVVLRHTASPSVTNVYILNLALADELFMLGLPFLAAQNALSY-WPFGSLMCRLVMAVDGINQFTSIFCLTVMSVDRYLAVVHPTSARWRTAPVARTVSAAVWVASAVVVLPVVVFSGVPR-----STCHMQ-WPEPAAAWRAGFIIYTAALGFFGPLLVICLCYLLIVVKVRSARVWAPSCQRRRRSERRVTRMVVAVVALFVLCWMPFYVLNIVNVVCPLPPAFFGLYFLVVALPYANSCANPILYGFLSYRFKQGFRRVLLRPSRRVRSQEPTVGEDEEEEDG---E
>O.1R_HUMAN
MGVPPGSREPSPVPPDYED-EFLRYLWRDYLYPKQYEWVLIAAYVAVFVVALVGNTLVCLAVWRNHHMRTVTNYFIVNLSLADVLVTAICLPASLLVDITESWLFGHALCKVIPYLQAVSVSVAVLTLSFIALDRWYAICHPL-LFKSTARRARGSILGIWAVSLAIMVPQAAVMECSSRTRLFSVCDER---WADDLYPKIYHSCFFIVTYLAPLGLMAMAYFQIFRKLWGRQPRFLAEVKQMRARRKTAKMLMVVLLVFALCYLPISVLNVLKRVFGMFEAVYACFTFSHWLVYANSAANPIIYNFLSGKFREQFKAAFSCCLPGLGPCGSLKASHKS---LSLQS
>O2C1_HUMAN
-------------MDGVNDSSLQGFVLMSISDHPQLEMIFFIAILFSYLLTLLGNSTIILLSRLEARLHTPMYFFLSNLSSLDLAFATSSVPQMLINLWGPGKTISYGGCITQLYVFLWLGATECILLVVMAFDRYVAVCRPLYTAIMNPQLCWLLAVIAWLGGLGNSVIQSTFTLQL---PEGFLCEVPKLACGDTSLNQAVLNGVCTFFTAVPLSIIVISYCLIAQAV--------LKIHSAEGRRKAFNTCLSHLLVVFLFYGSASYGYLLPAKNS---KQDQGKFISLFYSLVTPMVNPLIYTLRNMEVKGALRRLLGKGREVG--------------------
>OPSD_MOUSE
MNGTEGPNFYVPFSNVTGVGRSPFEQPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVVFTWIMALACAAPPLVGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIFFLICWLPYASVAFYIFTHQGSNFGPIFMTLPAFFAKSSSIYNPVIYIMLNKQFRNCMLTTLCCGKNPLGDDDASATASKTETSQVAPA
>CB1R_TARGR
EFFNRSVSTFKENDDNLKCGENFMDMECFMILTASQQLIIAVLSLTLGTFTVLENFLVLCVILQSRTRCRPSYHFIGSLAVADLLGSVIFVYSFLDFHVFHR-KDSSNVFLFKLGGVTASFTASVGSLFLTAIDRYISIHRPLYKRIVTRTKAVIAFCVMWTIAIIIAVLPLLGWNCK-K--LKSVCSDI------FPLIDENYLMFWIGVTSILLLFIVYAYVYILWKAHSHSEDQITRPEQTRMDIRLAKTLVLILVVLIICWGPLLAIMVYDVFGKMNNPIKTVFAFCSMLCLMDSTVNPIIYALRSQDLRHAFLEQCPPCEGTSQPLDNSMEGNN-AGNVHRAA
>B2AR_HUMAN
MGQPGNGSAFLLAPNRSHAPDH----DVTQQRDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGAAHILMKMWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRAI-N--TCCDFF--------TNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQGRTLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNL-IRKEVYILLNWIGYVNSGFNPLIYCRSP-DFRIAFQELLCLRRSSLKAYGNGYSGNTGEQSGYHVE
>OPSD_ZOSOP
MNGTEGPFFYIPMVNTTGIVRSPYEYPQYYLVNPAAYACLGAYMFFLILVGFPVNFLTLYVTLEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVLGRLGCNIEGFFATLGGEIALWSLVVLAIERWVVVCKPISNFRFGENHAIMGVAFTWFMASACAVPPLVGWSRYIPEGMQCSCGVDYYTRAEGFNNESFVIYMFIVHFCIPLAVVGFCYGRLLCAVKEAAAAQQESETTQRAEREVSRMVVIMVIGFLVCWLPYASVAWYIFTHQGSEFGPPFMTVPAFFAKSSSIYNPMIYICMNKQFRHCMITTLCCGKNPFEEEEGASTVSSSSVSPAA--
>GP37_HUMAN
LAQNGSLGEGIHEPGGPRRGNRLKNPFYPLTQESYGAYAVMCLSVVIFGTGIIGNLAVMCIVCHNYYMRSISNSLLANLAFWDFLIIFFCLPLVIFHELTKKWLLEDFSCKIVPYIEVASLGVTTFTLCALCIDRFRAATNVQYEMIENCSSTTAKLAVIWVGALLLALPEVVLRQLSKI-K--KISPDLTIYVLALTYDSARLWWYFGCYFCLPTLFTITCSLVTARKIRKAEKATRGNKRQIQLESQMNCTVVALTILYGFCIIPENICNIVTAYMATGQTMDLLNIISQFLLFFKSCVTPVLLFCLCKPFSRAFMECCCCCCEECIQKSSTVTYTTELELSPFST
>OPSR_CHICK
HEEEDTTRDSVFTYTNSNNTRGPFEGPNYHIAPRWVYNLTSVWMIFVVAASVFTNGLVLVATWKFKKLRHPLNWILVNLAVADLGETVIASTISVINQISGYFILGHPMCVVEGYTVSACGITALWSLAIISWERWFVVCKPFGNIKFDGKLAVAGILFSWLWSCAWTAPPIFGWSRYWPHGLKTSCGPDVFSGSSDPGVQSYMVVLMVTCCFFPLAIIILCYLQVWLAIRAVAAQQKESESTQKAEKEVSRMVVVMIVAYCFCWGPYTFFACFAAANPGYAFHPLAAALPAYFAKSATIYNPIIYVFMNRQFRNCILQLFGKKVDDG-----SEVVSSVSNSSVSPA
>IL8B_BOVIN
EGFEDEFGNYSGTPPTEDYDYSPC----EISTETLNKYAVVVIDALVFLLSLLGNSLVMLVILYSRIGRSVTDVYLLNLAMADLLFAMTLPIWTASKAKG--WVFGTPLCKVVSLLKEVNFYSGILLLACISMDRYLAIVHATRTLTQKWHWVKFICLGIWALSVILALPIFIFREAYQ-YS-DLVCYED-LGANTTKWRMIMRVLPQTFGFLLPLLVMLFCYGFTLRTL---------FSAQMGHKHRAMRVIFAVVLVFLLCWLPYNLVLIADTLMRAHNDIGRALDATEILGFLHSCLNPLIYVFIGQKFRHGLLKIMAIHGLISKEFLAKDG------------
>OPSD_MACFA
MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNAEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLFGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEARAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSASIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA
>5H2B_RAT
ILQKTCDHLILTDRSGLKAESAAEEMKQTAENQGNTVHWAALLIFAVIIPTIGGNILVILAVSLEKRLQYATNYFLMSLAVADLLVGLFVMPIALLTIMFEAWPLPLALCPAWLFLDVLFSTASIMHLCAISLDRYIAIKKPIANQCNSRTTAFVKITVVWLISIGIAIPVPIKGIEA-NA---ITCELT------KDRFGSFMLFGSLAAFFAPLTIMIVTYFLTIHALRKKRRMGKKPAQTISNEQRASKVLGIVFLFFLLMWCPFFITNVTLALCDSCTTLKTLLQIFVWVGYVSSGVNPLIYTLFNKTFREAFGRYITCNYQATKSVKVLRKGNSMVENSKFFT
>OPSR_ANOCA
NDEDDTTRDSLFTYTNSNNTRGPFEGPNYHIAPRWVYNITSVWMIFVVIASIFTNGLVLVATAKFKKLRHPLNWILVNLAIADLGETVIASTISVINQISGYFILGHPMCVLEGYTVSTCGISALWSLAVISWERWVVVCKPFGNVKFDAKLAVAGIVFSWVWSAVWTAPPVFGWSRYWPHGLKTSCGPDVFSGSDDPGVLSYMIVLMITCCFIPLAVILLCYLQVWLAIRAVAAQQKESESTQKAEKEVSRMVVVMIIAYCFCWGPYTVFACFAAANPGYAFHPLAAALPAYFAKSATIYNPIIYVFMNRQFRNCIMQLFGKKVDDG-----SELVSSVSNSSVSPA
>5H1A_MOUSE
-MDMFSLGQGNNTTTSLEPFGTGGNDTGLSNVTFSYQVITSLLLGTLIFCAVLGNACVVAAIALERSLQNVANYLIGSLAVTDLMVSVLVLPMAALYQVLNKWTLGQVTCDLFIALDVLCCTSSILHLCAIALDRYWAITDPIYVNKRTPRRAAALISLTWLIGFLISIPPMLGWRA---NP--NECTIS--------KDHGYTIYSTFGAFYIPLLLMLVLYGRIFRAARFRKNEEAKRKMALARERKTVKTLGIIMGTFILCWLPFFIVALVLPFCESSHMPELLGAIINWLGYSNSLLNPVIYAYFNKDFQNAFKKIIKCKFCR---------------------
>CCR4_HUMAN
IYTSDNYTEEMG-SGDYDSMKEPC---FREENANFNKIFLPTIYSIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVADLLFVITLPFWAVDAVAN--WYFGNFLCKAVHVIYTVNLYSSVLILAFISLDRYLAIVHATSQRPRKLLAEKVVYVGVWIPALLLTIPDFIFANVSE-DD-RYICDRF---YPNDLWVVVFQFQHIMVGLILPGIVILSCYCIIISKL---------SHSKGHQKRKALKTTVILILAFFACWLPYYIGISIDSFILLENTVHKWISITEALAFFHCCLNPILYAFLGAKFKTSAQHALTSVSRGSSLKILSKG----------GH
>CKR3_MACMU
MTTSLDTVETFGPTSYDDDMGLLC---EKADVGALIAQFVPPLYSLVFMVGLLGNVVVVMILIKYRRLRIMTNIYLLNLAISDLLFLFTLPFWIHYVRERN-WVFSHGMCKVLSGFYHTGLYSEIFFIILLTIDRYLAIVHAVALRARTVTFGVITSIVTWGLAVLAALPEFIFYGTEK-LF-KTLCSAIYPQDTVYSWRHFHTLKMTILCLALPLLVMAICYTGIIKTL---------LRCPSKKKYKAIRLIFVIMAVFFIFWTPYNVAILISTYQSVLKHLDLFVLATEVIAYSHCCVNPVIYAFVGERFRKYLRHFFHRHVLMHLGKYIPFLT-------SSVS
>OPSD_GLOME
MNGTEGLNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSVLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYIPLNLAVANLFMVFGGFTTTLYTSLHAYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLALTWVMAMACAAPPLVGWSRYIPEGMQCSCGIDYYTSRQEVNNESFVIYMFVVHFTIPLVIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVVAFLICWVPYASVAFYIFTHQGSDFGPIFMTIPSFFAKSSSIYNPVIYIMMNKQLRNCMLTTLCCGRNPLGDDEASTTASKTETSQVAPA
>FMLR_MOUSE
---MDTNMSLLMNKSAVNLMNVSGSTQSVSAGYIVLDVFSYLIFAVTFVLGVLGNGLVIWVAGFRMK-HTVTTISYLNLAIADFCFTSTLPFYIASMVMGGHWPFGWFMCKFIYTVIDINLFGSVFLIALIALDRCICVLHPVAQNHRTVSLAKKVIIVPWICAFLLTLPVIIRLTTVPSPWPVEKRKVA------VTMLTVRGIIRFIIGFSTPMSIVAICYGLITTKI---------HRQGLIKSSRPLRVLSFVVAAFFLCWCPFQVVALISTIQVREPGIVTALKITSPLAFFNSCLNPMLYVFMGQDFRERLIHSLPASLERALTEDSAQTGT----------
>OPSD_LAMJA
MNGTEGDNFYVPFSNKTGLARSPYEYPQYYLAEPWKYSALAAYMFFLILVGFPVNFLTLFVTVQHKKLRTPLNYILLNLAMANLFMVLFGFTVTMYTSMNGYFVFGPTMCSIEGFFATLGGEVALWSLVVLAIERYIVICKPMGNFRFGNTHAIMGVAFTWIMALACAAPPLVGWSRYIPEGMQCSCGPDYYTLNPNFNNESYVVYMFVVHFLVPFVIIFFCYGRLLCTVKEAAAAQQESASTQKAEKEVTRMVVLMVIGFLVCWVPYASVAFYIFTHQGSDFGATFMTLPAFFAKSSALYNPVIYILMNKQFRNCMITTLCCGKNPLGDDESGASVSSVSTSPVSPA
>C5AR_CAVPO
---MMVTVSYDYDYNSTFLPDGFVD--NYVERLSFGDLVAVVIMVVVFLVGVPGNALVVWVTACEAR-RHINAIWFLNLAAADLLSCLALPILLVSTVHLNHWYFGDTACKVLPSLILLNMYTSILLLATISADRLLLVLSPICQRFRGGCLAWTACGLAWVLALLLSSPSFLYRRTHN-SF--VYCVTD-YG-RDISKERAVALVRLLVGFIVPLITLTACYTFLLLRT---------WSRKATRSAKTVKVVVAVVSSFFVFWLPYQVTGILLAWHSPNRNTKALDAVCVAFAYINCCINPIIYVVAGHGFQGRLLKSLPSVLRNVLTEESLDKSTVD--------
>P2Y8_.ENLA
ATSYPTFLTTPYLPMKLLMNLTNDTEDICVFDEGFKFLLLPVSYSAVFMVGLPLNIAAMWIFIAKMRPWNPTTVYMFNLALSDTLYVLSLPTLVYYYADKNNWPFGEVLCKLVRFLFYANLYSSILFLTCISVHRYRGVCHPISLRRMNAKHAYVICALVWLSVTLCLVPNLIFVTVS-----NTICHDT-TRPEDFARYVEYSTAIMCLLFGIPCLIIAGCYGLMTRELMKP---SGNQQTLPSYKKRSIKTIIFVMIAFAICFMPFHITRTLYYYARLLNVINVTYKVTRPLASANSCIDPILYFLANDRYRRRLIRTVRRRSSVPNRRCMHTNMTAGPLPVISAE
>5H1D_FUGRU
ELDNNSLDYFSSNFTDIPSN--TTVAHWTEATLLGLQISVSVVLAIVTLATMLSNAFVIATIFLTRKLHTPANFLIGSLAVTDMLVSILVMPISIVYTVSKTWSLGQIVCDIWLSSDITFCTASILHLCVIALDRYWAITDALYSKRRTMRRAAVMVAVVWVISISISMPPLFWR----H-E--KECMVN-------TDQISYTLYSTFGAFYVPTVLLIILYGRIYVAARSRKLALERKRLCAARERKATKTLGIILGAFIICWLPFFVVTLVWAICKEC-FDPLLFDVFTWLGYLNSLINPVIYTVFNDEFKQAFQKLIKFRR-----------------------
>5HT_LYMST
TGQFINGSHSSRSRDNASANDTDDRYWSLTVYSHEHLVLTSVILGLFVLCCIIGNCFVIAAVMLERSLHNVANYLILSLAVADLMVAVLVMPLSVVSEISKVWFLHSEVCDMWISVDVLCCTASILHLVAIAMDRYWAVTSI-YIRRRSARRILLMIMVVWIVALFISIPPLFGWRD-PD-K--GTCIIS--------QDKGYTIFSTVGAFYLPMLVMMIIYIRIWLVARSRNDTRTREKLELKRERKAARTLAIITGAFLICWLPFFIIALIGPFVDPEGIPPFARSFVLWLGYFNSLLNPIIYTIFSPEFRSAFQKILFGKYRRGHR------------------
>5H6_HUMAN
-----------MVPEPGPTANSTPAWGAGPPSAPGGSGWVAAALCVVIALTAAANSLLIALICTQPARNTS-NFFLVSLFTSDLMVGLVVMPPAMLNALYGRWVLARGLCLLWTAFDVMCCSASILNLCLISLDRYLLILSPLYKLRMTPLRALALVLGAWSLAALASFLPLLLGWH---E-VPGQCRLL--------ASLPFVLVASGLTFFLPSGAICFTYCRILLAARKQVESRRLATKHSRKALKASLTLGILLGMFFVTWLPFFVANIVQAVCDC--ISPGLFDVLTWLGYCNSTMNPIIYPLFMRDFKRALGRFLPCPRCPRERQASLASSGPRPG-LSLQQ
>AG2R_MERUN
------MALNSSADDGIKRIQDDC---PKAGRHSYIFVMIPTLYSIIFVVGIFGNSLVVIVIYFYMKLKTVASVFLLNLALADLCFLLTLPVWAVYTAMEYRWPFGNHLCKIASAGISFNLYASVFLLTCLSIDRYLAIVHPMSRLRRTMLVAKVTCVVIWLLAGLASLPAVIHRNVYF-TN--TVCAFH-YESQNSTLPVGLGLTKNILGFMFPFLIILTSYTLIWKALKKA----YEIQKNKPRNDDIFRIIMAIVLFFFFSWIPHQIFTFLDVLIQLGDVVDTAMPITICIAYFNNCLNPLFYGFLGKKFKKYFLQLLKYIPPKAKSHSSLSTRPSD-------N
>GRHR_HORSE
--MANSDSLEQDPNHCSAINNSIPLIQGKLPTLTVSGKIRVTVTFFLFLLSTAFNASFLLKLQKWTQKLSRMKVLLKHLTLANLLETLIVMPLDGMWNITVQWYAGEFLCKVLSYLKLFSMYAPAFMMVVISLDRSLAITQPL-AVQSNSKLEQSMISLAWILSIVFAGPQLYIFRMIYTV--FSQCVTH--SFPQWWHQAFYNFFTFGCLFIIPLLIMLICNAKIIFALTRVPRKNQSKNNIPRARLRTLKMTVAFATSFVVCWTPYYVLGIWYWFDPEMRVSDPVNHFFFLFAFLNPCFDPLIYGYFSL-------------------------------------
>REIS_TODPA
---------------------MFGNPAMTGLHQFTMWEHYFTGSIYLVLGCVVFSLCGMCIIFLARQKPRRKYAILIHVLITAMAVNGGDPAHASSSIVGR-WLYGSVGCQLMGFWGFFGGMSHIWMLFAFAMERYMAVCHREFYQQMPSVYYSIIVGLMYTFGTFWATMPLLGWASYGLEVHGTSCTIN--YSVSDESYQSYVFFLAIFSFIFPMVSGWYAISKAWSGLSAIPD-EKEKDKDILSEEQLTALAGAFILISLISWSGFGYVAIYSALTHGGQLSHLRGHVPPIMSKTGCALFPLLIFLLTARSLPKSDTKKP--------------------------
>GRPR_RAT
NCSHLNLEVDPFLSCNNTFNQTLSPPKMDNWFHPGIIYVIPAVYGLIIVIGLIGNITLIKIFCTVKSMRNVPNLFISSLALGDLLLLVTCAPVDASKYLADRWLFGRIGCKLIPFIQLTSVGVSVFTLTALSADRYKAIVRPMIQASHALMKICLKAALIWIVSMLLAIPEAVFSDLHPQT--FISCAPY---HSNELHPKIHSMASFLVFYIIPLSIISVYYYFIARNLIQSLPVNIHVKKQIESRKRLAKTVLVFVGLFAFCWLPNHVIYLYRSYHYSEMLHFITSICARLLAFTNSCVNPFALYLLSKSFRKQFNTQLLCCQPSLLNR--SHSMTSFKSTNP-SA
>CKR5_MACMU
----MDYQVSSPTYDIDYYTSEPC---QKINVKQIAARLLPPLYSLVFIFGFVGNILVVLILINCKRLKSMTDIYLLNLAISDLLFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWVVAVFASLPGIIFTRSQR-EG-HYTCSSHFPYSQYQFWKNFQTLKMVILGLVLPLLVMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCSIFA-------SSVY
>B4AR_MELGA
--------------MTPLPAGNGSVPNCSWAAVLSRQWAVGAALSITILVIVAGNLLVIVAIAKTPRLQTMTNVFVTSLACADLVMGLLVVPPGATILLSGHWPYGTVVCELWTSLDVLCVTASIETLCAIAVDRYLAITAPLYEALVTKGRAWAVVCMVWAISAFISFLPIMNHWWRDV-R--RCCDFV--------TNMTYAIVSSTVSFYVPLLVMIFVYVRVFAVATRHSRGRRPSRLLAIKEHKALKTLGIIMGTFTLCWLPFFVANIIKVFCRPL-VPDQLFLFLNWLGYVNSAFNPIIYCRSP-DFRSAFRKLLCCPRRADRRLHAAPQAFSPRGDPMEDS
>OPSV_.ENLA
-----MLEEEDFYLFKNVSNVSPFDGPQYHIAPKWAFTLQAIFMGMVFLIGTPLNFIVLLVTIKYKKLRQPLNYILVNITVGGFLMCIFSIFPVFVSSSQGYFFFGRIACSIDAFVGTLTGLVTGWSLAFLAFERYIVICKPMGNFNFSSSHALAVVICTWIIGIVVSVPPFLGWSRYMPEGLQCSCGPDWYTVGTKYRSEYYTWFIFIFCFVIPLSLICFSYGRLLGALRAVAAQQQESASTQKAEREVSRMVIFMVGSFCLCYVPYAAMAMYMVTNRNHGLDLRLVTIPAFFSKSSCVYNPIIYSFMNKQFRGCIMETVCGRPMSD--DSSVSSSTVSSSQVSPA-
>OPSG_ORYLA
ENGTEGKNFYIPMNNRTGLVRSPYEYPQYYLADPWQFKLLGIYMFFLILTGFPINALTLVVTAQNKKLRQPLNFILVNLAVAGLIMVCFGFTVCIYSCMVGYFSLGPLGCTIEGFMATLGGQVSLWSLVVLAIERYIVVCKPMGSFKFTATHSAAGCAFTWIMASSCAVPPLVGWSRYIPEGIQVSCGPDYYTLAPGFNNESFVMYMFSCHFCVPVFTIFFTYGSLVMTVKAAAAQQQDSASTQKAEKEVTRMCFLMVLGFLLAWVPYASYAAWIFFNRGAAFSAMSMAIPSFFSKSSALFNPIIYILLNKQFRNCMLATIGMGG-----------VSTSKTEVSTAA
>NY6R_RABIT
--MEVSLNDPASNKTSAKSNSSAFFYFESCQSPSLALLLLLIAYTVVLIMGICGNLSLITIIFKKQRAQNVTNILIANLSLSDILVCVMCIPFTAIYTLMDRWIFGNTMCKLTSYVQSVSISVSIFSLVLIAIERYQLIVNPR-GWKPSASHAYWGIMLIWLFSLLLSIPLLLSYHLTDSH--HVVCVEH---WPSKTNQLLYSTSLIMLQYFVPLGFMFICYLKIVICLHKRNSKRRENESRLTENKRINTMLISIVVTFAACWLPLNTFNVIFDWYHEVCHHDLVFAICHLVAMVSTCINPLFYGFLNRNFQKDLVVLIHHCLCFALRERY---TLHTDESKGSLR
>CCR4_BOVIN
IFTSDNYTEDDLGSGDYDSMKEPC---FREENAHFNRIFLPTVYSIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVADLLFVLTLPFWAVDAVAN--WYFGKFLCKAVHVIYTVNLYSSVLILAFISLDRYLAIVHATSQKPRKLLAEKVVYVGVWLPAVLLTIPDLIFADIKE-DE-RYICDRF---YPSDLWLVVFQFQHIVVGLLLPGIVILSCYCIIISKL---------SHSKGYQKRKALKTTVILILTFFACWLPYYIGISIDSFILLESTVHKWISITEALAFFHCCLNPILYAFLGAKFKTSAQHALTSVSRGSSLKILSKG----------GH
>OPS2_SCHGR
YESSVGLPLLGWNVPTEHLDLVHPHWRSFQVPNKYWHFGLAFVYFMLMCMSSLGNGIVLWIYATTKSIRTPSNMFIVNLALFDVLMLLEMPMLVVSSLFYQR-PVWELGCDIYAALGSVAGIGSAINNAAIAFDRYRTISCPI-DGRLTQGQVLALIAGTWVWTLPFTLMPLLRIWSRFAEGFLTTCSFD--YLTDDEDTKVFVGCIFAWSYAFPLCLICCFYYRLIGAVREHNVKSNADTEAQSAEIRIAKVALTIFFLFLCSWTPYAVVAMIGAFGNRAALTPLSTMIPAVTAKIVSCIDPWVYAINHPRFRAEVQKRMKWLHLGEDARSSKSDRTVGNVSASA--
>DADR_HUMAN
--------------MRTLNTSAMDGTGLVVERDFSVRILTACFLSLLILSTLLGNTLVCAAVIRFRHRSKVTNFFVISLAVSDLLVAVLVMPWKAVAEIAGFWPFG-SFCNIWVAFDIMCSTASILNLCVISVDRYWAISSPFYERKMTPKAAFILISVAWTLSVLISFIPVQLSWHKAGN-TIDNCDSS--------LSRTYAISSSVISFYIPVAIMIVTYTRIYRIAQKQVECESSFKMSFKRETKVLKTLSVIMGVFVCCWLPFFILNCILPFCGSGCIDSNTFDVFVWFGWANSSLNPIIYAFNA-DFRKAFSTLLGCYRLCPATNNAIETAAMFSSH-----
>THRR_HUMAN
LTEYRLVSINKSSPLQKQLPAFISEDASGYLTSSWLTLFVPSVYTGVFVVSLPLNIMAIVVFILKMKVKKPAVVYMLHLATADVLFVSVLPFKISYYFSGSDWQFGSELCRFVTAAFYCNMYASILLMTVISIDRFLAVVYPMSLSWRTLGRASFTCLAIWALAIAGVVPLVLKEQTIQ--N--TTCHDVLNETLLEGYYAYYFSAFSAVFFFVPLIISTVCYVSIIRCL------SSSAVANRSKKSRALFLSAAVFCIFIICFGPTNVLLIAHYSFLSHEAAYFAYLLCVCVSSISSCIDPLIYYYASSECQRYVYSILCCKESSDPSSYNSSGDTCS--------
>V2R_HUMAN
MLMASTTSAVPGHPSLPSLPSNSSQERPLDTRDPLLARAELALLSIVFVAVALSNGLVLAALARRGRHWAPIHVFIGHLCLADLAVALFQVLPQLAWKATDRFRGPDALCRAVKYLQMVGMYASSYMILAMTLDRHRAICRPMAYRHGSGAHWNRPVLVAWAFSLLLSLPQLFIFAQRNSG--VTDCWAC---FAEPWGRRTYVTWIALMVFVAPTLGIAACQVLIFREIHASGRRPGEGAHVSAAVAKTVRMTLVIVVVYVLCWAPFFLVQLWAAWDPEA-LEGAPFVLLMLLASLNSCTNPWIYASFSSSVSSELRSLLCCARGRTPPSLGPQDSSLAKDTSS---
>GALR_HUMAN
----MELAVGNLSEGNASCPEPPAPEPGPLFGIGVENFVTLVVFGLIFALGVLGNSLVITVLARSKPPRSTTNLFILNLSIADLAYLLFCIPFQATVYALPTWVLGAFICKFIHYFFTVSMLVSIFTLAAMSVDRYVAIVHSRSSSLRVSRNALLGVGCIWALSIAMASPVAYHQGLF-SN--QTFCWEQ---WPDPRHKKAYVVCTFVFGYLLPLLLICFCYAKVLNHLHKKLK--NMSKKSEASKKKTAQTVLVVVVVFGISWLPHHIIHLWAEFGVFPPASFLFRITAHCLAYSNSSVNPIIYAFLSENFRKAYKQVFKCHIRKDSHLSDTKEPPSTNCTHV---
>5H7_RAT
SSWMPHLLSGFLEVTASPAPTNVSGCGEQINYGRVEKVVIGSILTLITLLTIAGNCLVVISVCFVKKLRQPSNYLIVSLALADLSVAVAVMPFVSVTDLIGGWIFGHFFCNVFIAMDVMCCTASIMTLCVISIDRYLGITRPLYPVRQNGKCMAKMILSVWLLSASITLPPLFGWAQ--N-D--KVCLIS--------QDFGYTIYSTAVAFYIPMSVMLFMYYQIYKAARKSSRLERKNISIFKREQKAATTLGIIVGAFTVCWLPFFLLSTARPFICGTCIPLWVERTCLWLGYANSLINPFIYAFFNRDLRTTYRSLLQCQYRNINRKLSAAGAERPERSEFVLQ
>UL33_RCMVM
----MDVLLGTEELEDELHQLHFNYTCVPSLGLSVARDAETAVNFLIVLVGGPMNFLVLATQMLSNRSVSTPTLYMTNLYLANLLTVATLPFLMLSNRGL--VGSSPEGCKIAALAYYATCTAGFATLMLIAINRYR-VIHQRRSGAGSKRQTYAVLAVTWLASLMCASPAPLYATVMAA-DAFETCIIYSYDQVK-TVLATFKILITMIWGITPVVMMSWFYVFFYRRL---------KLTSYRRRSQTLTFVTTLMLSFLVVQTPFVAIMSYDSYGVLNNKRDAVSMLARVVPNFHCLLNPVLYAFLGRDFNKRFILCISGKLFSRRRALRERAGPVCALP---SK
>O.YR_MOUSE
GTPAANWSIELDLGSGVPPGAEGNLTAGPPRRNEALARVEVAVLCLILFLALSGNACVLLALRTTRHKHSRLFFFMKHLSIADLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPL--RSLRRRTDRLAVLATWLGCLVASVPQVHIFSLRE----VFDCWAV---FIQPWGPKAYVTWITLAVYIVPVIVLAACYGLISFKIWQNRAAVSSVKLISKAKIRTVKMTFIIVLAFIVCWTPFFFVQMWSVWDVNA-KEASAFIIAMLLASLNSCCNPWIYMLFTGHLFHELVQRFLCCSARYLKGSRPGENSSTFVLSRCSS
>PF2R_RAT
--MSINS---------SKQPASSAAGLIANTTCQTENRLSVFFSIIFMTVGIVSNSLAIAILMKAYQSKASFLLLASGLVITDFFGHLINGGIAVFVYASDKFDQSNILCSVFGISMVFSGLCPLFLGSTMAIERCIGVTNPLHSTKITSKHVKMILSGVCMFAVFVALLPILGHRDYQIQASRTWCFYN--TEHIEDWEDRFYLLFFSSLGLLALGISFSCNAVTGVTLLRVKFRSQQHRQGRSHHLEMVIQLLAIMCVSCVCWSPFLVTMANIAINGNNPVTCETTLFALRMATWNQILDPWVYILLRKAVLRNLYKLASRCCGVNIISLHIWELKVAAISESPAA
>D1DR_FUGRU
--------------MAQNFSTVGDGKQMLLERDSSKRVLTGCFLSLLIFTTLLGNTLVCVAVTKFRHRSKVTNFFVISLAISDLLVAILVMPWKAATEIMGFWPFG-EFCNIWVAFDIMCSTASILNLCVISVDRYWAISSPFYERKMTPKVACLMISVAWTLSVLISFIPVQLNWHKALN-PPDNCDSS--------LNRTYAISSSLISFYIPVAIMIVTYTRIYRIAQKQSLSECSFKMSFKRETKVLKTLSVIMGVFVCCWLPFFILNCMVPFCEADCISSTTFDVFVWFGWANSSLNPIIYAFNA-DFRKAFSILLGCHRLCPGNS-AIEIAPLSNPSCQYQP
>CKR2_RAT
HSLFPRSIQELDEGATTPYDYDDGEPCHKTSVKQIGAWILPPLYSLVFIFGFVGNMLVIIILISCKKLKSMTDIYLFNLAISDLLFLLTLPFWAHYAANE--WVFGNIMCKLFTGLYHIGYFGGIFFIILLTIDRYLAIVHAVALKARTVTFGVITSVVTWVVAVFASLPGIIFTKSEQ-ED-QHTCGPY----FPTIWKNFQTIMRNILSLILPLLVMVICYSGILHTL--------FRCRNEKKRHRAVRLIFAIMIVYFLFWTPYNIVLFLTTFQEFLMHLDQAMQVTETLGMTHCCVNPIIYAFVGEKFRRYLSIFFRKHIAKNLCKQCPVFVSSTFTPSTGEQ
>OPRK_HUMAN
PPNSSAWFPGWAEPDSNGSAGSEDAQLEPAHISPAIPVIITAVYSVVFVVGLVGNSLVMFVIIRYTKMKTATNIYIFNLALADALVTTTMPFQSTVYLMNS-WPFGDVLCKIVISIDYYNMFTSIFTLTMMSVDRYIAVCHPVALDFRTPLKAKIINICIWLLSSSVGISAIVLGGTKVDVD-VIECSLQFPDDDYSWWDLFMKICVFIFAFVIPVLIIIVCYTLMILRLKSVRLL-SGSREKDRNLRRITRLVLVVVAVFVVCWTPIHIFILVEALGSTSTAALSSYYFCIALGYTNSSLNPILYAFLDENFKRCFRDFCFPLKMRMERQSTSRVAYLRDIDGMNKP
>ACM3_BOVIN
N---------ISRAAGNLSSPNGTTSDPLGGHTIWQVVFIAFLTGVLALVTIIGNILVIVAFKVNKQLKTVNNYFLLSLACADLIIGVISMNLFTTYIIMNRWALGNLACDLWLSIDYVASNASVMNLLVISFDRYFSITRPLYRAKRTTKRAGVMIGLAWVISFILWAPAILFWQYFV-VP-PGECFIQ------FLSEPTITFGTAIAAFYMPVTIMTILYWRIYKETEKRKTRTKRKRMSLIKEKKAAQTLSAILLAFIITWTPYNIMVLVNTFCDSC-IPKTYWNLGYWLCYINSTVNPVCYALCNKTFRNTFKMLLLCQCDKRKRRKQQYQHKRVPEQAL---
>C3.1_HUMAN
---MDQFPESVTENFEYDDLAEAC---YIGDIVVFGTVFLSIFYSVIFAIGLVGNLLVVFALTNSKKPKSVTDIYLLNLALSDLLFVATLPFWTHYLINEKG--LHNAMCKFTTAFFFIGFFGSIFFITVISIDRYLAIVLAASMNNRTVQHGVTISLGVWAAAILVAAPQFMFTKQKE-----NECLGDYPEVLQEIWPVLRNVETNFLGFLLPLLIMSYCYFRIIQTL---------FSCKNHKKAKAIKLILLVVIVFFLFWTPYNVMIFLETLKLYDKDLRLALSVTETVAFSHCCLNPLIYAFAGEKFRRYLYHLYGKCLAVLCGRSVHVDRSRHGSVLSSNF
>5H4_HUMAN
------------------MDKLDANVSSEEGFGSVEKVVLLTFLSTVILMAILGNLLVMVAVCWDRQRKIKTNYFIVSLAFADLLVSVLVMPFGAIELVQDIWIYGEVFCLVRTSLDVLLTTASIFHLCCISLDRYYAICCQPYRNKMTPLRIALMLGGCWVIPTFISFLPIMQGWNNIRK-NSTYCVFM--------VNKPYAITCSVVAFYIPFLLMVLAYYRIYVTAKEHSRPDQHSTHRMRTETKAAKTLCIIMGCFCLCWAPFFVTNIVDPFIDYT-VPGQVWTAFLWLGYINSGLNPFLYAFLNKSFRRAFLIILCCDDERYRRPSILGQTINGSTHVLRDA
>5H2A_PIG
NTSDAFNWTVDSENRTNLSCEGCLSPPCFSLLHLQEKNWSALLTAVVIILTIAGNILVIMAVSLEKKLQNATNYFLMSLAIADMLLGFLVMPVSMLTILYGYWPLPSKLCAVWIYLDVLFSTASIMHLCAISLDRYVAIQNPIHRRFNSRTKAFLKIIAVWTISVGISMPIPVFGLQD-VF---GSCLL---------ADDNFVLIGSFVSFFIPLTIMVITYFLTIKSLQKEREPGRRTMQSISNEQKACKVLGIVFFLFVVMWCPFFITNIMAVICKESDVIGALLNVFVWIGYLSSAVNPLVYTLFNKTYRSAFSRYIQCQYKENKKPLQLILAYKSSQLQTGQK
>OPSD_DIPVU
MNGTEGPYFYVPMVNTSGIVRSPYEYPQYYLVNPAAYAALGAYMFLLILVGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVLGRLGCNIEGFFATLGGEIALWSLVVLAIERWVVVCKPISNFRFGENHAIMGLAFTWLMALACAAPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFICHFSIPLLVVFFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVIMMVIAFLVCWLPYASVAWWIFTHQGSDFGPVFMTIPAFFAKSSSIYNPMIYICLNKQFRHCMITTLCCGKNPFEEEEGASTSVSSSSVSPAA-
>HH1R_HUMAN
----------MSLPNSSCLLEDKMCEGNKTTMASPQLMPLVVVLSTICLVTVGLNLLVLYAVRSERKLHTVGNLYIVSLSVADLIVGAVVMPMNILYLLMSKWSLGRPLCLFWLSMDYVASTASIFSVFILCIDRYRSVQQPLYLKYRTKTRASATILGAWFLSFLWVIPILGWNHFM--VR-EDKCETD------FYDVTWFKVMTAIINFYLPTLLMLWFYAKIYKAVRQHLRSQYVSGLHMNRERKAAKQLGFIMAAFILCWIPYFIFFMVIAFCKNC-CNEHLHMFTIWLGYINSTLNPLIYPLCNENFKKTFKRILHIRS-----------------------
>OPSD_DELDE
MNGTEGLNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSVLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVANLFMVFGGFTTTLYTSLHAYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLALTWIMAMACAAPPLVGWSRYIPEGMQCSCGIDYYTLSPEVNNESFVIYMFVVHFTIPLVIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVVAFLICWVPYASVAFYIFTHQGSDFGPIFMTIPSFFAKSSSIYNPVIYIMMNKQFRNCMLTTLCCGRNPLGDDEASTTASKTETSQVAPA
>OPS1_DROME
PLS---NGSVVDKVTPDMAHLISPYWNQFPAMDPIWAKILTAYMIMIGMISWCGNGVVIYIFATTKSLRTPANLLVINLAISDFGIMITNTPMMGINLYFETWVLGPMMCDIYAGLGSAFGCSSIWSMCMISLDRYQVIVKGMAGRPMTIPLALGKIAYIWFMSSIWCLAPAFGWSRYVPEGNLTSCGID--YLERDWNPRSYLIFYSIFVYYIPLFLICYSYWFIIAAVSAHMNVRSSEDAEKSAEGKLAKVALVTITLWFMAWTPYLVINCMGLFKFEG-LTPLNTIWGACFAKSAACYNPIVYGISHPKYRLALKEKCPCCVFGKVDDGKSSDSEAESKA-----
>RGR_BOVIN
---------------------MAESGTLPTGFGELEVLAVGTVLLVEALSGLSLNILTILSFCKTPELRTPSHLLVLSLALADSGIS-LNALVAATSSLLRRWPYGSEGCQAHGFQGFVTALASICSSAAVAWGRYHHFCTRS---RLDWNTAVSLVFFVWLSSAFWAALPLLGWGHYDYEPLGTCCTLD--YSRGDRNFTSFLFTMAFFNFLLPLFITVVSYRLME--------------QKLGKTSRPPVNTVLPARTLLLGWGPYALLYLYATIADATSISPKLQMVPALIAKAVPTVNAMNYALGSEMVHRGIWQCLSPQRREHSREQ----------------
>VU51_HSV6U
----------------MEKETKSLAWPATAEFYGWVFIFSSIQLCTMVLLTVRFNSFKVGR-E--------YAVFTFAGMSFNCFLLPIKMGLLSGH-----WSLPRDFCAILLYIDDFSIYFSSWSLVFMAIERINHFCYSTLLNENSKALAKVCFPIVWIISGVQALQMLNNYKATA---ETPQCFLA--------FRSGYDMWLMLVYSVMIPVMLVFIYIYSKNFM-----------LLKDELSTVTTYLCIYLLLGTIAHLPKAGLSEIESD----KIFYGLRDIFMALPVLKVYYIPVMAYCMACDDHTVPVRLCSIWLVNLCKKCFSCTLEVGIKMLK---
>A1AA_ORYLA
-----------MTPSSVTLNCSNCSHVLAPELNTVKAVVLGMVLGIFILFGVIGNILVILSVVCHRHLQTVTYYFIVNLAVADLLLSSTVLPFSAIFEILDRWVFGRVFCNIWAAVDVLCCTASIMSLCVISVDRYIGVSYPLYPAIMTKRRALLAVMLLWVLSVIISIGPLFGWKEP--AP--TVCKIT--------EEPGYAIFSAVGSFYLPLAIILAMYCRVYVVAQKELRSFALRLLKFSREKKAAKTLGIVVGCFVLCWLPFFLVLPIGSIFPAYRPSDTVFKITFWLGYFNSCINPIIYLCSNQEFKKAFQSLLGVHCLRMTPRAHHHHTQGHSLT-----
>B1AR_BOVIN
VPDGAATAARLLVPASPPASLLTSASEGPPLPSQQWTAGMGLLMAFIVLLIVVGNVLVLVAIAKTPRLQTLTNLFIMSLASADLVMGLLVVPFGATIVVWGRWEYGSFFCELWTSVDVLCVTASIETLCVIALDRYLAITSPFYQSLLTRARARALVCTVWAISALVSFLPIFMQWWGDS-R--ECCDFI--------INEGYAITSSVVSFYVPLCIMAFVYLRVFREAQKQPGRRRPPRLVALREQKALKTLGIIMGVFTLCWLPFFLANVVKAFHRDL-VPDRLFVFFNWLGYANSAFNPIIYCRSP-DFRKAFQRLLCCARRAACGSHAAAGCLAVARPSPSPG
>OPSH_ASTFA
NE--DTTRESAFVYTNANNTRDPFEGPNYHIAPRWVYNVSSLWMIFVVIASVFTNGLVIVATAKFKKLRHPLNWILVNLAIADLGETVLASTISVINQIFGYFILGHPMCVFEGWTVSVCGITALWSLTIISWERWVVVCKPFGNVKFDGKWAAGGIIFSWVWAIIWCTPPIFGWSRYWPHGLKTSCGPDVFSGSEDPGVASYMITLMLTCCILPLSIIIICYIFVWSAIHQVAQQQKDSESTQKAEKEVSRMVVVMILAFIVCWGPYASFATFSAVNPGYAWHPLAAAMPAYFAKSATIYNPIIYVFMNRQFRSCIMQLFGKKVEDA-----SEVSTAS--------
>OL15_MOUSE
-------------MEVDSNSSSGTFILMGVSDHPHLEIIFFAVILASYLLTLVGNLTIILLSRLDARLHTPMYFFLSNLSSLDLAFTTSSVPQMLKNLWGPDKTISYGGCVTQLYVFLWLGATECILLVVMAFDRYVAVCRPLYMTVMNPRLCWGLAAISWLGGLGNSVIQSTFTLQL---PDNFLCEVPKLACGDTSLNEAVLNGVCTFFTVVPVSVILVSYCFIAQAV--------MKIRSVEGRRKAFNTCVSHLVVVFLFYGSAIYGYLLPAKSS---NQSQGKFISLFYSVVTPMVNPLIYTLRNKEVKGALGRLLGKGRGAS--------------------
>ETBR_HUMAN
SLAPAEVPKGDRTAGSPPRTISPPPCQGPIEIKETFKYINTVVSCLVFVLGIIGNSTLLRIIYKNKCMRNGPNILIASLALGDLLHIVIDIPINVYKLLAEDWPFGAEMCKLVPFIQKASVGITVLSLCALSIDRYRAVASWSIKGIGVPKWTAVEIVLIWVVSVVLAVPEAIGFDIITRI--LRICLLHQKTAFMQFYKTAKDWWLFSFYFCLPLAITAFFYTLMTCEMLRKSGM-IALNDHLKQRREVAKTVFCLVLVFALCWLPLHLSRILKLTLYNQSFLLVLDYIGINMASLNSCINPIALYLVSKRFKNCFKSCLCCWCQSFE-EKQSLEFKANDHGYDNFR
>PE23_MOUSE
---------MASMWAPEHSAEAHSNLS---STTDDCGSVSVAFPITMMVTGFVGNALAMLLVSRSYRRKKSFLLCIGWLALTDLVGQLLTSPVVILVYLSQRLDPSGRLCTFFGLTMTVFGLSSLLVASAMAVERALAIRAPHYASHMKTRAT-PVLLGVWLSVLAFALLPVLGVG----RYSGTWCFISNETDPAREPGSVAFASAFACLGLLALVVTFACNLATIKALVSRAKASQSSAQWGRITTETAIQLMGIMCVLSVCWSPLLIMMLKMIFNQMSEKECNSFLIAVRLASLNQILDPWVYLLLRKILLRKFCQIRDHTN-YASSSTSLPCWSDQLER-----
>TA2R_MOUSE
--MWPNG-----------TSLGACFRPVNITLQERRAIASPWFAASFCALGLGSNLLALSVLAGARPPRSSFLALLCGLVLTDFLGLLVTGAIVASQHAALLTDPSCRLCYFMGVAMVFFGLCPLLLGAAMASERFVGITRPFSRPTATSRRAWATVGLVWVAAGALGLLPLLGLGRYSVQYPGSWCFLT----LGTQRGDVVFGLIFALLGSASVGLSLLLNTVSVATLCRVYHTREATQRPRDCEVEMMVQLVGIMVVATVCWMPLLVFIMQTLLQTPPRATEHQLLIYLRVATWNQILDPWVYILFRRSVLRRLHPRFSSQLQAVSLRRPPAQ------------
>ACM3_HUMAN
N---------VSRAAGNFSSPDGTTDDPLGGHTVWQVVFIAFLTGILALVTIIGNILVIVSFKVNKQLKTVNNYFLLSLACADLIIGVISMNLFTTYIIMNRWALGNLACDLWLAIDYVASNASVMNLLVISFDRYFSITRPLYRAKRTTKRAGVMIGLAWVISFVLWAPAILFWQYFV-VP-PGECFIQ------FLSEPTITFGTAIAAFYMPVTIMTILYWRIYKETEKRKTRTKRKRMSLVKEKKAAQTLSAILLAFIITWTPYNIMVLVNTFCDSC-IPKTFWNLGYWLCYINSTVNPVCYALCNKTFRTTFKMLLLCQCDKKKRRKQQYQHKRAPEQAL---
>CCR4_FELCA
IYPSDNYTEDDLGSGDYDSMKEPC---FREENAHFNRIFLPTVYSIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVADLLFVLTLPFWAVDAVAN--WYFGKFLCKAVHVIYTVNLYSSVLILAFISLDRYLAIVHATSQRPRKLLAEKVVYVGVWIPALLLTIPDFIFANVRE-DG-RYICDRF---YPSDSWLVVFQFQHIMVGLILPGIVILSCYCIIISKL---------SHSKGYQKRKALKTTVILILAFFACWLPYYIGISIDSFILLESTVHKWISITEALAFFHCCLNPILYAFLGAKFKTSAQHALTSVSRGSSLKILSKG----------GH
>GPRJ_MOUSE
AEAAEALLPHGLMGLHEEHSWMSNRTELQYELNPGEVATASIFFGALWLFSIFGNSLVCLVIHRSRRTQSTTNYFVVSMACADLLISVASTPFVVLQFTTGRWTLGSAMCKVVRYFQYLTPGVQIYVLLSICIDRFYTIVYPL-SFKVSREKAKKMIAASWILDAAFVTPVFFFYG---SNW-HCNYFLP-----PSWEGTAYTVIHFLVGFVIPSILIILFYQKVIKYIWRIDGR-RTMNIVPRTKVKTVKMFLLLNLVFLFSWLPFHVAQLWHPHEQDYKKSSLVFTAVTWVSFSSSASKPTLYSIYNANFRRGMKETFCMSSMKCYRSNAYTIKRNYVGISEIPP
>AA2B_MOUSE
------------------------------MQLETQDALYVALELVIAALAVAGNVLVCAAVGASSALQTPTNYFLVSLATADVAVGLFAIPFAITISLG--FCTDFHGCLFLACFVLVLTQSSIFSLLAVAVDRYLAIRVPLYKGLVTGTRARGIIAVLWVLAFGIGLTPFLGWNSKDI-A-PLTCLFE-----NVVPMSYMVYFNFFGCVLPPLLIMLVIYIKIFMVACKQ---MDHSRTTLQREIHAAKSLAMIVGIFALCWLPVHAINCITLFHPALDKPKWVMNVAILLSHANSVVNPIVYAYRNRDFRYSFHKIISRYVLCQAETKGGSGLSLGL-------
>MAS_HUMAN
-----MDGSNVTSFVVEEPTNISTGRNASVGNAHRQIPIVHWVIMSISPVGFVENGILLWFLCFRMR-RNPFTVYITHLSIADISLLFCIFILSIDYALDYESSGHYYTIVTLSVTFLFGYNTGLYLLTAISVERCLSVLYPIYRCHRPKYQSALVCALLWALSCLVTTMEYVMCI-DRHS--RNDC-----------RAVIIFIAILSFLVFTPLMLVSSTILVVKIRK----------NTWASHSSKLYIVIMVTIIIFLIFAMPMRLLYLLYYEYW--STFGNLHHISLLFSTINSSANPFIYFFVGSSKKKRFKESLKVVLTRAFKDEMQPRTVTVETVV----
>TSHR_SHEEP
LQAFDNHYDYTVCGGSEEMVCTPKSDEFNPCEDIMGYKFLRIVVWFVSLLALLGNVFVLVILLTSHYKLTVPRFLMCNLAFADFCMGLYLLLIASVDLYTQSWQTG-PGCNTAGFFTVFASELSVYTLTVITLERWYAITFAMLDRKIRLWHAYVIMLGGWVCCFLLALLPLVGISSY----KVSICLPM-----TETPLALAYIILVLLLNIIAFIIVCACYVKIYITVRNP------HYNPGDKDTRIAKRMAVLIFTDFMCMAPISFYALSALMNKPLITVTNSKILLVLFYPLNSCANPFLYAIFTKAFQRDVFMLLSKFGICKRQAQAYRGSTGIRVQKVPPD
>NK2R_BOVIN
----MGACVVMTDINISSGLDSNATGITAFSMPGWQLALWTAAYLALVLVAVMGNATVIWIILAHQRMRTVTNYFIVNLALADLCMAAFNAAFNFVYASHNIWYFGRAFCYFQNLFPITAMFVSIYSMTAIAADRYMAIVHPF-QPRLSAPGTRAVIAGIWLVALALAFPQCFYST---GA---TKCVVAWPEDSGGKMLLLYHLIVIALIYFLPLVVMFVAYSVIGLTLWRRGHQHGANLRHLQAKKKFVKTMVLVVVTFAICWLPYHLYFILGTFQEDIKFIQQVYLALFWLAMSSTMYNPIIYCCLNHRFRSGFRLAFRCCPWVTPTEE----YTPSLSTRVNRC
>5H2C_MOUSE
LLVWQFDISISPVAAIVTDTFNSSDGGRLFQFPDGVQNWPALSIVVIIIMTIGGNILVIMAVSMEKKLHNATNYFLMSLAIADMLVGLLVMPLSLLAILYDYWPLPRYLCPVWISLDVLFSTASIMHLCAISLDRYVAVRSPVHSRFNSRTKAIMKIAIVWAISIGVSVPIPVIGLRD-VF--NTTCVL---------NDPNFVLIGSFVAFFIPLTIMVITYFLTIYVLRRQKKKPRGTMQAINNEKKASKVLGIVFFVFLIMWCPFFITNILSVLCGKAKLMEKLLNVFVWIGYVCSGINPLVYTLFNKIYRRAFSKYLRCDYKPDKKPPVRQIALSGRELNVNIY
>OLF1_CHICK
-------------MASGNCTTPTTFILSGLTDNPGLQMPLFMVFLAIYTITLLTNLGLIALISVDLHLQTPMYIFLQNLSFTDAAYSTVITPKMLATFLEERKTISYVGCILQYFSFVLLTVTESLLLAVMAYDRYVAICKPLYPSIMTKAVCWRLVESLYFLAFLNSLVHTSGLLKL---SNHFFCDISQISSSSIAISELLVIISGSLFVMSSIIIILISYVFIILTV--------VMIRSKDGKYKAFSTCTSHLMAVSLFHGTVIFMYLRPVKLF---SLDTDKIASLFYTVVIPMLNPLIYSWRNKEVKDALRRLTATTFGFIDSKAVQ--------------
>O.2R_RAT
ASELNETQEPFLNPTDYDDEEFLRYLWREYLHPKEYEWVLIAGYIIVFVVALIGNVLVCVAVWKNHHMRTVTNYFIVNLSLADVLVTITCLPATLVVDITETWFFGQSLCKVIPYLQTVSVSVSVLTLSCIALDRWYAICHPL-MFKSTAKRARNSIVVIWIVSCIIMIPQAIVMERSSKTTLFTVCDER---WGGEVYPKMYHICFFLVTYMAPLCLMVLAYLQIFRKLWCRKARVAAEIKQIRARRKTARMLMVVLLVFAICYLPISILNVLKRVFGMFETVYAWFTFSHWLVYANSAANPIIYNFLSGKFREEFKAAFSCCLG-VHRRQGDRLESRKSLTTQISN
>CKR8_MOUSE
DYTMEPNVTMT--DYYPDFFTAPC---DAEFLLRGSMLYLAILYCVLFVLGLLGNSLVILVLVGCKKLRSITDIYLLNLAASDLLFVLSIPFQTHNLLDQ--WVFGTAMCKVVSGLYYIGFFSSMFFITLMSVDRYLAIVHAVAIKVRTASVGTALSLTVWLAAVTATIPLMVFYQVAS-ED-MLQCFQF-YEEQSLRWKLFTHFEINALGLLLPFAILLFCYVRILQQL---------RGCLNHNRTRAIKLVLTVVIVSLLFWVPFNVALFLTSLHDLHQRLALAIHVTEVISFTHCCVNPVIYAFIGEKFKKHLMDVFQKS-CSHIFLYLGRQRQ------LSSN
>GPRF_MACNE
----MDPEETSVYLDYYYATSPNPDIRETHSHVPYTSVFLPVFYTAVFLTGVLGNLVLMGALHFKPGSRRLIDIFIINLAASDFIFLVTLPLWVDKEASLGLWRTGSFLCKGSSYMISVNMHCSVFLLTCMSVDRYLAIVCPVSRKFRRTDCAYVVCASIWFISCLLGLPTLLSRELT-IDD-KPYCAEK----KATPLKLIWSLVALIFTFFVPLLSIVTCYCCIARKLCAH---YQQSGKHNKKLKKSIKIIFIVVAAFLVSWLPFNTSKLLAIVSGLQAILQLGMEVSGPLAFANSCVNPFIYYIFDSYIRRAIVHCLCPCLKNYDFGSSTETALSTFIHAEDFT
>MC3R_MOUSE
NSSCCLSSVSPMLPNLSEHPAAPPASNRSGSGFCEQVFIKPEVFLALGIVSLMENILVILAVVRNGNLHSPMYFFLCSLAAADMLVSLSNSLETIMIAVINSDQFIQHMDNIFDSMICISLVASICNLLAIAIDRYVTIFYALYHSIMTVRKALTLIGVIWVCCGICGVMFIIYYS----------------------EESKMVIVCLITMFFAMVLLMGTLYIHMFLFARLHIAVAGVVAPQQHSCMKGAVTITILLGVFIFCWAPFFLHLVLIITCPTNICYTAHFNTYLVLIMCNSVIDPLIYAFRSLELRNTFKEILCGCNSMNLG------------------
>CCR5_RAT
DDLYKELAIYSNSTEIPLQDSIFCSTEEGPLLTSFKTIFMPVAYSLIFLLGMMGNILVLVILERHRHTRSSTETFLFHLAVADLLLVFILPFAVAEGSVG--WVLGTFLCKTVIALHKINFYCSSLLLACIAVDRYLAIVHAVAYRRRRLLSIHITCSTIWLAGFLFALPELLFAKVVQ-NE-LPQCIFSQENEAETRAWFASRFLYHTGGFLLPMLVMAWCYVGVVHRL--------LQAQRRPQRQKAVRVAILVTSIFLLCWSPYHIVIFLDTLERLKGYLSVAITLCEFLGLAHCCLNPMLYTFAGVKFRSDLSRLLTKLGCAG---PASLC--------PGWR
>OPSD_LOLFO
---MGRDIPDNETWWYNPYMDIHPHWKQFDQVPAAVYYSLGIFIAICGIIGCVGNGVVIYLFTKTKSLQTPANMFIINLAFSDFTFSLVNGFPLMTISCFMKWVFGNAACKVYGLIGGIFGLMSIMTMTMISIDRYNVIGRPMASKKMSHRKAFIMIIFVWIWSTIWAIGPIFGWGAYTLEGVLCNCSFD--YITRDTTTRSNILCMYIFAFMCPIVVIFFCYFNIVMSVSNHRLNLRKAQAGANAEMKLAKISIVIVTQFLLSWSPYAVVALLAQFGPIEWVTPYAAQLPVMFAKASAIHNPMIYSVSHPKFRERIASNFPWILTCCQYDEKEIEEIPAGEQS-GGE
>OPSD_TODPA
----GRDLRDNETWWYNPSIVVHPHWREFDQVPDAVYYSLGIFIGICGIIGCGGNGIVIYLFTKTKSLQTPANMFIINLAFSDFTFSLVNGFPLMTISCFLKWIFGFAACKVYGFIGGIFGFMSIMTMAMISIDRYNVIGRPMASKKMSHRRAFIMIIFVWLWSVLWAIGPIFGWGAYTLEGVLCNCSFD--YISRDSTTRSNILCMFILGFFGPILIIFFCYFNIVMSVSNHRLNLRKAQAGANAEMRLAKISIVIVSQFLLSWSPYAVVALLAQFGPLEWVTPYAAQLPVMFAKASAIHNPMIYSVSHPKFREAISQTFPWVLTCCQFDDKETEEIPAGESSDAAP
>OPRK_RAT
LPNSSSWFPNWAESDSNGSVGSEDQQLEPAHISPAIPVIITAVYSVVFVVGLVGNSLVMFVIIRYTKMKTATNIYIFNLALADALVTTTMPFQSAVYLMNS-WPFGDVLCKIVISIDYYNMFTSIFTLTMMSVDRYIAVCHPVALDFRTPLKAKIINICIWLLASSVGISAIVLGGTKVDVD-VIECSLQFPDDEYSWWDLFMKICVFVFAFVIPVLIIIVCYTLMILRLKSVRLL-SGSREKDRNLRRITKLVLVVVAVFIICWTPIHIFILVEALGSTSTAVLSSYYFCIALGYTNSSLNPVLYAFLDENFKRCFRDFCFPIKMRMERQSTNRVASMRDVGGMNKP
>ACM4_RAT
------M-NFTPVNGSSANQSVRLVTAAHNHLETVEMVFIATVTGSLSLVTVVGNILVMLSIKVNRQLQTVNNYFLFSLGCADLIIGAFSMNLYTLYIIKGYWPLGAVVCDLWLALDYVVSNASVMNLLIISFDRYFCVTKPLYPARRTTKMAGLMIAAAWVLSFVLWAPAILFWQFVV-VP-DNQCFIQ------FLSNPAVTFGTAIAAFYLPVVIMTVLYIHISLASRSRSIAVRKKRQMAARERKVTRTIFAILLAFILTWTPYNVMVLVNTFCQSC-IPERVWSIGYWLCYVNSTINPACYALCNATFKKTFRHLLLCQYRNIGTAR----------------
>SSR3_RAT
SVPTTLDPGNASSAWPLDTSLGNASAGTSLAGLAVSGILISLVYLVVCVVGLLGNSLVIYVVLRHTSSPSVTSVYILNLALADELFMLGLPFLAAQNALSY-WPFGSLMCRLVMAVDGINQFTSIFCLTVMSVDRYLAVVHPTSARWRTAPVARMVSAAVWVASAVVVLPVVVFSGVPR-----STCHMQ-WPEPAAAWRTAFIIYTAALGFFGPLLVICLCYLLIVVKVRSTSCQAPACQRRRRSERRVTRMVVAVVALFVLCWMPFYLLNIVNVVCPLPPAFFGLYFLVVALPYANSCANPILYGFLSYRFKQGFRRILLRPSRRVRSQEPGSGEEDEEEEERREE
>OPSD_HUMAN
MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA
>OPSG_GECGE
RDDDDTTRGSVFTYTNTNNTRGPFEGPNYHIAPRWVYNLVSFFMIIVVIASCFTNGLVLVATAKFKKLRHPLNWILVNLAFVDLVETLVASTISVFNQIFGYFILGHPLCVIEGYVVSSCGITGLWSLAIISWERWFVVCKPFGNIKFDSKLAIIGIVFSWVWAWGWSAPPIFGWSRYWPHGLKTSCGPDVFSGSVELGCQSFMLTLMITCCFLPLFIIIVCYLQVWMAIRAVAAQQKESESTQKAEREVSRMVVVMIVAFCICWGPYASFVSFAAANPGYAFHPLAAALPAYFAKSATIYNPVIYVFMNRQFRNCIMQLFGKKVDDG-----SEAVSSVSNSSVAPA
>OAR_BOMMO
TEIYDVIEDEKDVCAVADEPNIPCSFGISLAVPEWEAICTAIILTMIIISTVVGNILVILSVFTYKPLRIVQNFFIVSLAVADLTVAILVLPLNVAYSILGQWVFGIYVCKMWLTCDIMCCTSSILNLCAIALDRYWAITDPIYAQKRTLERVLFMIGIVWILSLVISSPPLLGWNDW-E-P--TPCRLT--------SQPGFVIFSSSGSFYIPLVIMTVVYFEIYLATKKRAVYEEKQRISLTRERRAARTLGIIMGVFVVCWLPFFVIYLVIPFCVSCCLSNKFINFITWLGYVNSALNPLIYTIFNMDFRRAFKKLLFIKC-----------------------
>YQH2_CAEEL
EQSTPARENLPNREIYQIFQFTLVYALPLSNHDNSSLMLIAGFYALLFMFGTCGNAAILAVVHHVKGRHNTTLTYICILSIVDFLSMLPIPMTIIDQILGF-WMFDTFACKLFRLLEHIGKIFSTFILVAFSIDRYCAVCHPLQVRVRNQRTVFVFLGIMFFVTCVMLSPILLYAHSKVTRMHLYKCVDD----LGRELFVVFTLYSFVLAYLMPLLFMIYFYYEMLIRLFKQLVGGGEEKKLTIPVGHIAIYTLAICSFHFICWTPYWISILYSLYEELYYAFIYFMYGVHALPYINSASNFILYGLLNRQLHNAPERKYTRNGVGGRQMSHALTSELIAIPSSSCR
>UL33_HSV7J
----MICYSFAKNVTFAFLIILQNFFSQHDEEYKYNYTCITPTVRKAQRLESVINGIMLTLILPVSTKQTITSPYLITLFISDSLHSLTVLLLTLNREAL--TNLNQALCQCVLFVYSASCTYSLCMLAVISTIRYR-TLQRRTLNDKNNNHIKRNVGILFLSSAMCAIPAVLYVQVEKK-KNYGKCNIHSTQKAY-DLFIGIKIVYCFLWGIFPTVIFSYFYVIFGKTL---------RALTQSKHNKTLSFISLLILSFLCIQIPNLLVMSVEIFFLYIIQREIVQIISRLMPEIHCLSNPLVYAFTRTDFRLRFYDFIKCNLCNSSLKRKRNP------------
>C5AR_RAT
DPISNDSSEITYDYSDGTPNPDMPADGVYIPKMEPGDIAALIIYLAVFLVGVTGNALVVWVTAFEAK-RTVNAIWFLNLAVADLLSCLALPILFTSIVKHNHWPFGDQACIVLPSLILLNMYSSILLLATISADRFLLVFKPICQKFRRPGLAWMACGVTWVLALLLTIPSFVFRRIHK-SD--ILCNID-YSKGPFFIEKAIAILRLMVGFVLPLLTLNICYTFLLIRT---------WSRKATRSTKTLKVVMAVVTCFFVFWLPYQVTGVILAWLPRSQSVERLNSLCVSLAYINCCVNPIIYVMAGQGFHGRLRRSLPSIIRNVLSEDSLGRSTMD--------
>CKR6_HUMAN
DSSEDYFVSVNTSYYSVDSEMLLC---SLQEVRQFSRLFVPIAYSLICVFGLLGNILVVITFAFYKKARSMTDVYLLNMAIADILFVLTLPFWAVSHATGA-WVFSNATCKLLKGIYAINFNCGMLLLTCISMDRYIAIVQATRLRSRTLPRSKIICLVVWGLSVIISSSTFVFNQKYN-TQ-SDVCEPKQTVSEPIRWKLLMLGLELLFGFFIPLMFMIFCYTFIVKTL---------VQAQNSKRHKAIRVIIAVVLVFLACQIPHNMVLLVTAANLGKKLIGYTKTVTEVLAFLHCCLNPVLYAFIGQKFRNYFLKILKDLWCVRRKYKSSGFEN------ISRQ
>B2AR_MOUSE
MGPHGNDSDFLLAPNGSRAPDH----DVTQERDEAWVVGMAILMSVIVLAIVFGNVLVITAIAKFERLQTVTNYFIISLACADLVMGLAVVPFGASHILMKMWNFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYVAITSPFYQSLLTKNKARVVILMVWIVSGLTSFLPIQMHWYRAI-D--TCCDFF--------TNQAYAIASSIVSFYVPLCVMVFVYSRVFQVAKRQGRSLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIRDNL-IPKEVYILLNWLGYVNSAFNPLIYCRSP-DFRIAFQELLCLRRSSSKTYGNGYSDYTGEPNTCQLG
>DADR_DIDMA
---------------MPLNDTTMDRRGLVVERDFSFRILTACFLSLLILSTLLGNTLVCAAVIRFRHRSKVTNFFVISLAVSDLLVAVLVMPWKAVAEIAGFWPFG-SFCNIWVAFDIMCSTASILNLCVISVDRYWAISSPFYERKMTPKAAFILISVAWTLSVLISFIPVQLNWHKAGN-TMDNCDSS--------LSRTYAISSSLISFYIPVAIMIVTYTRIYRIAQKQVECESSFKMSFKRETKVLKTLSVIMGVFVCCWLPFFILNCMVPFCESDCIDSITFDVFVWFGWANSSLNPIIYAFNA-DFRKAFSTLLGCYRLCPTANNAIETGAVFSSH-----
>O.YR_MACMU
GELAANWSTEAVNSSAAPPGAEGNCTAGPPRRNEALARVEVAVLCLILFLALSGNACVLLALRTTRHKHSRLFFFMKHLSIADLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPL--RSLRRRTDRLAVLATWLGCLVASAPQVHIFSLRE----VFDCWAV---FIQPWGPKAYITWITLAVYIVPVIVLAACYGLISFKIWQNRMAVSSVKLISKAKIRTVKMTFIIVLAFIVCWTPFFFVQMWSVWDANA-KEASAFIIVMLLASLNSCCNPWIYMLFTGHLFHELVQRFLCCSASYLKGNRLGENSSSFVLSHRSS
>AG2R_MELGA
------MVPNYSTEETVKRIHVDC---PVSGRHSYIYIMVPTVYSIIFIIGIFGNSLVVIVIYCYMKLKTVASIFLLNLALADLCFLITLPLWAAYTAMEYQWPFGNCLCKLASAGISFNLYASVFLLTCLSIDRYLAIVHPVSRIRRTMFVARVTCIVIWLLAGVASLPVIIHRNIFF-LN--TVCGFR-YDNNNTTLRVGLGLSKNLLGFLIPFLIILTSYTLIWKTLKKA----YQIQRNKTRNDDIFKMIVAIVFFFFFSWIPHQVFTFLDVLIQLHDIVDTAMPFTICIAYFNNCLNPFFYVFFGKNFKKYFLQLIKYIPPNVSTHPSLTTRPPE-------N
>TSHR_MOUSE
LQAFESHYDYTVCGDNEDMVCTPKSDEFNPCEDIMGYRFLRIVVWFVSLLALLGNIFVLLILLTSHYKLTVPRFLMCNLAFADFCMGVYLLLIASVDLYTHSWQTG-PGCNTAGFFTVFASELSVYTLTVITLERWYAITFAMLDRKIRLRHAYTIMAGGWVSCFLLALLPMVGISSY----KVSICLPM-----TDTPLALAYIVLVLLLNVVAFVVVCSCYVKIYITVRNP------QYNPRDKDTKIAKRMAVLIFTDFMCMAPISFYALSALMNKPLITVTNSKILLVLFYPLNSCANPFLYAIFTKAFQRDVFILLSKFGICKRQAQAYQGSTGIQIQKIPQD
>CCR3_MOUSE
FAFLLENSTSPYDYGENESDFSDSPPCPQDFSLNFDRTFLPALYSLLFLLGLLGNGAVAAVLLSQRTALSSTDTFLLHLAVADVLLVLTLPLWAVDAAVQ--WVFGPGLCKVAGALFNINFYAGAFLLACISFDRYLSIVHATIYRRDPRVRVALTCIVVWGLCLLFALPDFIYLSANY-RL-ATHCQYN----FPQVGRTALRVLQLVAGFLLPLLVMAYCYAHILAVL---------LVSRGQRRFRAMRLVVVVVAAFAVCWTPYHLVVLVDILMDVGSHVDVAKSVTSGMGYMHCCLNPLLYAFVGVKFREKMWMLFTRLGRSDQRGPQRQPSWSETTEASYLG
>AA2B_RAT
------------------------------MQLETQDALYVALELVIAALAVAGNVLVCAAVGASSALQTPTNYFLVSLATADVAVGLFAIPFAITISLG--FCTDFHSCLFLACFVLVLTQSSIFSLLAVAVDRYLAIRVPLYKGLVTGTRARGIIAVLWVLAFGIGLTPFLGWNSKDI-T-PVKCLFE-----NVVPMSYMVYFNFFGCVLPPLLIMMVIYIKIFMVACKQ---MEHSRTTLQREIHAAKSLAMIVGIFALCWLPVHAINCITLFHPALDKPKWVMNVAILLSHANSVVNPIVYAYRNRDFRYSFHRIISRYVLCQTDTKGGSGFSLSL-------
>MSHR_RANTA
PVLGSQRRLLGSLNCTPPATFPLMLAPNRTGPQCLEVSIPNGLFLSLGLVSLVENVLVVAAIAKNSNLHSPMYYFICCLAVSDLLVSVSNVLETAVMLLLEAAAVVQQLDNVIDVLICGSMVSSLCFLGAIAVDRYISIFYALYHSVVTLPRAWRIIAAIWVASILTSLLFITYYY----------------------NNHTVVLLCLVGFFIAMLALMAVLYVHMLARACQHIARKRQRPIHRGFGLKGAATLTILLGVFFLCWGPFFLHLSLIVLCPQHGCIFKNFNLFLALIICNAIVDPLIYAFRSQELRKTLQEVLQCSW-----------------------
>CKR3_RAT
EEELKTVVETFETTPYEYEWAPPC---EKVSIRELGSWLLPPLYSLVFIVGLLGNMMVVLILIKYRKLQIMTNIYLLNLAISDLLFLFTVPFWIHYVLWNE-WGFGHCMCKMLSGLYYLALYSEIFFIILLTIDRYLAIVHAVALRARTVTFATITSIITWGFAVLAALPEFIFHESQD-NF-DLSCSPRYPEGEEDSWKRFHALRMNIFGLALPLLIMVICYSGIIKTL---------LRCPNKKKHKAIQLIFVVMIVFFIFWTPYNLVLLLSAFHSTFIHLDLAMQVTEVITHTHCCINPIIYAFVGERFRKHLRLFFHRNVAIYLRKYISFLT-------SSVS
>OPS2_HEMSA
DFGYPEGVSIVDFVRPEIKPYVHQHWYNYPPVNPMWHYLLGVIYLFLGTVSIFGNGLVIYLFNKSAALRTPANILVVNLALSDLIMLTTNVPFFTYNCFSGGWMFSPQYCEIYACLGAITGVCSIWLLCMISFDRYNIICNGFNGPKLTTGKAVVFALISWVIAIGCALPPFFGWGNYILEGILDSCSYD--YLTQDFNTFSYNIFIFVFDYFLPAAIIVFSYVFIVKAIFAHMNVRSNEADAQRAEIRIAKTALVNVSLWFICWTPYALISLKGVMGDTSGITPLVSTLPALLAKSCSCYNPFVYAISHPKYRLAITQHLPWFCVHETETKSNDDAQDKA-------
>AG2R_RABIT
------MMLNSSTEDGIKRIQDDC---PKAGRHNYIFVMIPTLYSIIFVVGIFGNSLAVIVIYFYMKLKTVASVFLLNLALADLCFLLTLPLWAVYTAMEYRWPFGNYLCKIASASVSFNLYASVFLLTCLSIDRYLAIVHPMSRLRRTMLVAKVTCIIIWLLAGLASLPAIIHRNVFF-TN--TVCAFH-YESQNSTLPIGLGLTKNILGFLFPFLIILTSYTLIWKALKKA----YEIQKNKPRNDDIFKIIMAIVLFFFFSWVPHQIFTFLDVLIQLGDIVDTAMPITICIAYFNNCLNPLFYGFLGKKFKKYFLQLLKYIPPKAKSHSNLSTRPSD-------N
>5H2A_HUMAN
NTSDAFNWTVDSENRTNLSCEGCLSPSCLSLLHLQEKNWSALLTAVVIILTIAGNILVIMAVSLEKKLQNATNYFLMSLAIADMLLGFLVMPVSMLTILYGYWPLPSKLCAVWIYLDVLFSTASIMHLCAISLDRYVAIQNPIHSRFNSRTKAFLKIIAVWTISVGISMPIPVFGLQD-VF---GSCLL---------ADDNFVLIGSFVSFFIPLTIMVITYFLTIKSLQKEEPGGRRTMQSISNEQKACKVLGIVFFLFVVMWCPFFITNIMAVICKESDVIGALLNVFVWIGYLSSAVNPLVYTLFNKTYRSAFSRYIQCQYKENKKPLQLILAYKSSQLQMGQK
>AA2A_HUMAN
-------------------------------MPIMGSSVYITVELAIAVLAILGNVLVCWAVWLNSNLQNVTNYFVVSLAAADIAVGVLAIPFAITISTG--FCAACHGCLFIACFVLVLTQSSIFSLLAIAIDRYIAIRIPLYNGLVTGTRAKGIIAICWVLSFAIGLTPMLGWNNCGS-Q-QVACLFE-----DVVPMNYMVYFNFFACVLVPLLLMLGVYLRIFLAARRQESQGERARSTLQKEVHAAKSLAIIVGLFALCWLPLHIINCFTFFCPDCHAPLWLMYLAIVLSHTNSVVNPFIYAYRIREFRQTFRKIIRSHVLRQQEPFKAAGAHGSDGEQVSLR
>HH2R_CAVPO
-------------------MAFNGTVPSFCMDFTVYKVTISVILIILILVTVAGNVVVCLAVGLNRRLRSLTNCFIVSLAVTDLLLGLLVLPFSAIYQLSCKWSFSKVFCNIYTSLDVMLCTASILNLFMISLDRYCAVTDPLYPVLITPARVAISLVFIWVISITLSFLSIHLGWN--RN-TIVKCKVQ--------VNEVYGLVDGLVTFYLPLLIMCITYFRIFKIAREQR--IGSWKAATIREHKATVTLAAVMGAFIICWFPYFTVFVYRGLKGDD-VNEVFEDVVLWLGYANSALNPILYAALNRDFRTAYHQLFCCRLASHNSHETSLRRSQCQEPRW---
>IL8B_HUMAN
KGEDLSNYSYSSTLPPFLLDAAPC----EPESLEINKYFVVIIYALVFLLSLLGNSLVMLVILYSRVGRSVTDVYLLNLALADLLFALTLPIWAASKVNG--WIFGTFLCKVVSLLKEVNFYSGILLLACISVDRYLAIVHATRTLTQKRYLVKFICLSIWGLSLLLALPVLLFRRTVY-NV-SPACYED-MGNNTANWRMLLRILPQSFGFIVPLLIMLFCYGFTLRTL---------FKAHMGQKHRAMRVIFAVVLIFLLCWLPYNLVLLADTLMRTQNHIDRALDATEILGILHSCLNPLIYAFIGQKFRHGLLKILAIHGLISKDSLPKDS------------
>OLF5_RAT
-------------MSSTNQSSVTEFLLLGLSRQPQQQQLLFLLFLIMYLATVLGNLLIILAIGTDSRLHTPMYFFLSNLSFVDVCFSSTTVPKVLANHILGSQAISFSGCLTQLYFLAVFGNMDNFLLAVMSYDRFVAICHPLYTTKMTRQLCVLLVVGSWVVANMNCLLHILLMARL---SPHFFCDGTKLSCSDTHLNELMILTEGAVVMVTPFVCILISYIHITCAV--------LRVSSPRGGWKSFSTCGSHLAVVCLFYGTVIAVYFNPSSSH---LAGRDMAAAVMYAVVTPMLNPFIYSLRNSDMKAALRKVLAMRFPSKQ-------------------
>OPSB_CHICK
TDLPEDFYIPMALDAPNITALSPFLVPQTHLGSPGLFRAMAAFMFLLIALGVPINTLTIFCTARFRKLRSHLNYILVNLALANLLVILVGSTTACYSFSQMYFALGPTACKIEGFAATLGGMVSLWSLAVVAFERFLVICKPLGNFTFRGSHAVLGCVATWVLGFVASAPPLFGWSRYIPEGLQCSCGPDWYTTDNKWHNESYVLFLFTFCFGVPLAIIVFSYGRLLITLRAVARQQEQSATTQKADREVTKMVVVMVLGFLVCWAPYTAFALWVVTHRGRSFEVGLASIPSVFSKSSTVYNPVIYVLMNKQFRSCMLKLLFCGRSPFGDDEDVSGSSVSSSHVAPA-
>AA2A_RAT
----------------------------------MGSSVYITVELAIAVLAILGNVLVCWAVWINSNLQNVTNFFVVSLAAADIAVGVLAIPFAITISTG--FCAACHGCLFFACFVLVLTQSSIFSLLAIAIDRYIAIRIPLYNGLVTGVRAKGIIAICWVLSFAIGLTPMLGWNNCST-K-RVTCLFE-----DVVPMNYMVYYNFFAFVLLPLLLMLAIYLRIFLAARRQESQGERTRSTLQKEVHAAKSLAIIVGLFALCWLPLHIINCFTFFCSTCHAPPWLMYLAIILSHSNSVVNPFIYAYRIREFRQTFRKIIRTHVLRRQEPFQAGGAHSTEGEQVSLR
>OLF6_CHICK
-------------MASGNCTTPTTFILSGLTDNPGLQMPLFMVFLAIYTITLLTNLGLIALISIDLQLQTPMYIFLQNLSFTDAVYSTVITPKMLATFLEETKTISYVGCILQYFSFVLLTVRECLLLAVMAYDRYAAICKPLYPAIMTKAVCWRLVKGLYSLAFLNFLVHTSGLLKL---SNHFFCDNSQISSSSTALNELLVFIFGSLFVMSSIITILISYVFIILTV--------VRIRSKERKYKAFSTCTSHLMAVSLFHGTIVFMYFQPANNF---SLDKDKIMSLFYTVVIPMLNPLIYSWRNKEVKDALHRAIATAVLFH--------------------
>P2Y6_HUMAN
-----------MEWDNGTGQALGLPPTTCVYRENFKQLLLPPVYSAVLAAGLPLNICVITQICTSRRALTRTAVYTLNLALADLLYACSLPLLIYNYAQGDHWPFGDFACRLVRFLFYANLHGSILFLTCISFQRYLGICHPLWHKRGGRRAAWLVCVAVWLAVTTQCLPTAIFAATG-----RTVCYDL-SPPALATHYMPYGMALTVIGFLLPFAALLACYCLLACRLCRQ---GPAEPVAQERRGKAARMAVVVAAAFAISFLPFHITKTAYLAVRSTEAFAAAYKGTRPFASANSVLDPILFYFTQKKFRRRPHELLQKLTAKWQRQGR---------------
>OPS5_DROME
SLGDGSVFPMGHGYPAEYQHMVHAHWRGFREAPIYYHAGFYIAFIVLMLSSIFGNGLVIWIFSTSKSLRTPSNLLILNLAIFDLFMCTN-MPHYLINATVGYIVGGDLGCDIYALNGGISGMGASITNAFIAFDRYKTISNPI-DGRLSYGQIVLLILFTWLWATPFSVLPLFQIWGRYPEGFLTTCSFD--YLTNTDENRLFVRTIFVWSYVIPMTMILVSYYKLFTHVRVHNVKANANADNMSVELRIAKAALIIYMLFILAWTPYSVVALIGCFGEQQLITPFVSMLPCLACKSVSCLDPWVYATSHPKYRLELERRLPWLGIREKHATSGTSSVSGDTLALSVQ
>YR13_CAEEL
--------------------MSNIFSVPLDPMSVAVGIPYVCFFIILSVVGIIGNVIVIYAIAGDRNRKSVMNILLLNLAVADLANLIFTIPEWIPPVFFGSWLFPSFLCPVCRYLECVFLFASISTQMIVCIERYIAIVLPMARQLCSRRNVLITVLVDWIFVACFASPYAVWHSVK-LFQLSATCSNT---VGKSTWWQGYKLTEFLAFYFVPCFIITVVYTKVAKCLWCKCLDSSRSSDALRTRRNVVKMLIACVAVYFVCYSPIQVIFLSKAVLNVTHPPYDFILLMNALAMTCSASNPLLYTLFSQKFRRRLRDVLYCPSDVENETKTYYSGPRASF------
>5H1B_CRIGR
CAPPPPAASQTGVPLVNLSHNSAESHIYQDSIALPWKVLLVALLALITLATTLSNAFVIATVYRTRKLHTPANYLIASLAVTDLLVSILVMPVSTMYTVTGRWTLGQVVCDFWLSSDITCCTASIMHLCVIALDRYWAITDAVYAAKRTPKRAAIMIALVWVFSISISLPPFFWR----E-E--LTCLVN-------TDHVLYTVYSTGGAFYLPTLLLIALYGRIYVEARSRRVSLEKKKLMAARERKATKTLGIILGAFIVCWLPFFIISLVMPICKDAWFHMATLDFFNWLGYLNSLINPIIYTMSNEDFKQAFHKLIRFKCAG---------------------
>D2DR_BOVIN
---MDPLNLSWYDDDPESRNWSRPFNGSEGKADRPPYNYYAMLLTLLIFVIVFGNVLVCMAVSREKALQTTTNYLIVSLAVADLLVATLVMPWVVYLEVVGEWKFSRIHCDIFVTLDVMMCTASILNLCAISIDRYTAVAMPMNTRYSSKRRVTVMIAIVWVLSFTISCPMLFGLN---Q----NECII---------ANPAFVVYSSIVSFYVPFIVTLLVYIKIYIVLRRRRTSMSRRKLSQQKEKKATQMLAIVLGVFIICWLPFFITHILNIHCDCN-IPPVLYSAFTWLGYVNSAVNPIIYTTFNIEFRKAFLKILHC-------------------------
>VU51_HSV6Z
----------------MEKETKSLAWPATAEFYGWVFIFSSIQLCTVVFLTVRFNGFKVGR-E--------YAVFTFAGMSFNCFLLPIKMGLLSGH-----WTLPRDFCAILLYIDDFSAYFSSWSLVFMAIERINYFCYSTLLNENSKALAKVCFPIVWVVSGVQALQMLNNYKATA---ETGQCFLA--------FRSGHDMWLMLVYSVVIPVMLVFFYLYSKNFM-----------LLKDELSSVTTYLCIYLLLGTIAHLPKAALSEIESD----KIFYGLRDIFMALPVLKVYYISAMAYCMACDDHTVPVRLCSIWLVNLCKKCFSCTLEVGIKMLK---
>D2DR_MOUSE
---MDPLNLSWYDDDLERQNWSRPFNGSEGKADRPHYNYYAMLLTLLIFIIVFGNVLVCMAVSREKALQTTTNYLIVSLAVADLLVATLVMPWVVYLEVVGEWKFSRIHCDIFVTLDVMMCTASILNLCAISIDRYTAVAMPMNTRYSSKRRVTVMIAIVWVLSFTISCPLLFGLN---Q----NECII---------ANPAFVVYSSIVSFYVPFIVTLLVYIKIYIVLRKRRTSMSRRKLSQQKEKKATQMLAIVLGVFIICWLPFFITHILNIHCDCN-IPPVLYSAFTWLGYVNSAVNPIIYTTFNIEFRKAFMKILHC-------------------------
>CCR4_MACFA
IYTSDNYTEEMG-SGDYDSIKEPC---FREENAHFNRIFLPTIYSIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVADLLYVITLPFWAVDAVAN--WYFGNFLCKAVHVIYTVNLYSSVLILAFISLDRYLAIVHATSQRPRKLLAEKVVYVGVWIPALLLTIPDFIFASVSE-DD-RYICDRF---YPNDLWVVVFQFQHIMVGLILPGIVILSCYCIIISKL---------SHSKGHQKRKALKTTVILILAFFACWLPYYIGISIDSFILLENTVHKWISITEALGFFHCCLNPILYAFLGAKFKTSAQHALTSVSRGSSLKILSKG----------GH
>HH1R_BOVIN
---------MTCPNSSCVFEDKMCQGNKTAPANDAQLTPLVVVLSTISLVTVGLNLLVLYAVRSERKLHTVGNLYIVSLSVADLIVGVVVMPMNILYLLMSRWSLGRPLCLFWLSMDYVASTASIFSVFILCIDRYRSVQQPLYLRYRTKTRASITILAAWFLSFLWIIPILGWRHFQ--EP-EDKCETD------FYNVTWFKVMTAIINFYLPTLLMLWFYAKIYKAVRQHLRSQYVSGLHMNRERKAAKQLGFIMAAFIICWIPYFIFFMVIAFCESC-CNQHVHMFTIWLGYINSTLNPLIYPLCNENFKKTFKKILHIRS-----------------------
>OPS2_DROPS
AQTG-GNRSVLDNVLPDMAPLVNPHWSRFAPMDPTMSKILGLFTLVILIISCCGNGVVVYIFGGTKSLRTPANLLVLNLAFSDFCMMASQSPVMIINFYYETWVLGPLWCDIYAACGSLFGCVSIWSMCMIAFDRYNVIVKGINGTPMTIKTSIMKIAFIWMMAVFWTIMPLIGWSSYVPEGNLTACSID--YMTRQWNPRSYLITYSLFVYYTPLFMICYSYWFIIATVAAHMNVRSSEDCDKSAENKLAKVALTTISLWFMAWTPYLIICYFGLFKIDG-LTPLTTIWGATFAKTSAVYNPIVYGISHPNDRLVLKEKCPMCVCGTTDEPKPDATSEAESKD----
>OPSD_SCYCA
MNGTEGENFYIPMSNKTGVVRSPFDYPQYYLAEPWKFSVLAAYMFFLIIAGFPVNFLTLYVTIQHKKLRQPLNYILLNLAVADLFMIFGGFPSTMITSMNGYFVFGPSGCNFEGFFATLGGEIGLWSLVVLAIERYVVVCKPMSNFRFGSQHAFMGVGLTWIMAMACAFPPLVGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFSIPLTIIFFCYGRLVCTVKEAAAQQQESETTQRAEREVTRMVIIMVIAFLICWLPYASVAFFIFCNQGSEFGPIFMTIPAFFAKAASLYNPLIYILMNKQFRNCMITTICCGKNPFEEEESTSASSVSSSQVAPAA
>OPSD_BUFBU
MNGTEGPNFYIPMSNKTGVVRSPFEYPQYYLAEPWQYSILCAYMFLLILLGFPINFMTLYVTIQHKKLRTPLNYILLNLAFANHFMVLCGFTVTMYSSMNGYFILGATGCYVEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFSENHAVMGVAFTWIMALSCAVPPLLGWSRYIPEGMQCSCGVDYYTLKPEVNNESFVIYMFVVHFTIPLIIIFFCYGRLVCTVKEAAAQQQESATTQKAEKEVTRMVIIMVVFFLICWVPYASVAFFIFSNQGSEFGPIFMTVPAFFAKSSSIYNPVIYIMLNKQFRNCMITTLCCGKNPFGEDDASSAASSVSSSQVSPA
>ACM4_CHICK
AQPWQAKMANLTYDNVTLSNRSEVAIQPPTNYKTVELVFIATVTGSLSLVTVVGNILVMLSIKVNRQLQTVNNYFLFSLACADLIIGVFSMNLYTVYIIKGYWPLGAVVCDLWLALDYVVSNASVMNLLIISFDRYFCVTKPLYPARRTTKMAGLMIAAAWILSFILWAPAILFWQFIV-VH-ERECYIQ------FLSNPAVTFGTAIAAFYLPVVIMTVLYIHISLASRSRSIAVRKKRQMAAREKKVTRTIFAILLAFILTWTPYNVMVLINTFCETC-VPETVWSIGYWLCYVNSTINPACYALCNATFKKTFKHLLMCQYRNIGTAR----------------
>GRHR_RAT
--MANNASLEQDQNHCSAINNSIPLTQGKLPTLTLSGKIRVTVTFFLFLLSTAFNASFLVKLQRWTQKLSRMKVLLKHLTLANLLETLIVMPLDGMWNITVQWYAGEFLCKVLSYLKLFSMYAPAFMMVVISLDRSLAVTQPL-AVQSKSKLERSMTSLAWILSIVFAGPQLYIFRMIYAV--FSQCVTH--SFPQWWHEAFYNFFTFSCLFIIPLLIMLICNAKIIFALTRVPRKNQSKNNIPRARLRTLKMTVAFGTSFVICWTPYYVLGIWYWFDPEMRVSEPVNHFFFLFAFLNPCFDPLIYGYFSL-------------------------------------
>C5AR_MOUSE
---MNSSFEINYDHY-GTMDPNIPADGIHLPKRQPGDVAALIIYSVVFLVGVPGNALVVWVTAFEPD-GPSNAIWFLNLAVADLLSCLAMPVLFTTVLNHNYWYFDATACIVLPSLILLNMYASILLLATISADRFLLVFKPICQKVRGTGLAWMACGVAWVLALLLTIPSFVYREAYK-SE--TVCGIN-YGGGSFPKEKAVAILRLMVGFVLPLLTLNICYTFLLLRT---------WSRKATRSTKTLKVVMAVVICFFIFWLPYQVTGVMIAWLPPSKRVEKLNSLCVSLAYINCCVNPIIYVMAGQGFHGRLLRSLPSIIRNALSEDSVGRSTDD--------
>5H2A_MACMU
NTSDAFNWTVESENRTNLSCEGCLSPSCLSLLHLQEKNWSALLTAVVIILTIAGNILVIMAVSLEKKLQNATNYFLMSLAIADMLLGFLVMPVSMLTILYGYWPLPSKLCAVWIYLDVLFSTASIMHLCAISLDRYVAIQNPIHSRFNSRTKAFLKIIAVWTISVGISMPIPVFGLQD-VF---GSCLL---------ADDNFVLIGSFVSFFIPLTIMVITYFLTIKSLQKEDPGGRRTMQSISNEQKACKVLGIVFFLFVVMWCPFFITNIMAVICKESDVIGALLNVFVWIGYLSSAVNPLVYTLFNKTYRSAFSRYIQCQYKENKKPLQLILAYKSSQLQMGQK
>LSHR_BOVIN
ESELSGWDYDYGFCLPKTLQCAPEPDAFNPCEDIMGYNFLRVLIWLINILAITGNVTVLFVLLTSRYKLTVPRFLMCNLSFADFCMGLYLLLIASVDAQTKGWQTG-SGCSAAGFFTVFASELSVYTLTVITLERWHTITYAILDQKLRLKHAIPVMLGGWLFSTLIAVLPLVGVSNY----KVSICLPM-----VESTLSQVYILTILILNVMAFIIICACYIKIYFAVQNP------ELMATNKDTKIAKKMAVLIFTDFTCMAPISFFAISAAFKVPLITVTNSKVLLVLFYPVNSCANPFLYAIFTKAFQRDFFLLLSKFGCCKYRAELYRRSNCKNGFTGSNK
>SSR5_HUMAN
LFPASTPSWNASSPGAASGGGDNRTLVGPAPSAGARAVLVPVLYLLVCAAGLGGNTLVIYVVLRFAKMKTVTNIYILNLAVADVLYMLGLPFLATQNAASF-WPFGPVLCRLVMTLDGVNQFTSVFCLTVMSVDRYLAVVHPLSARWRRPRVAKLASAAAWVLSLCMSLPLLVFADVQE-----GTCNAS-WPEPVGLWGAVFIIYTAVLGFFAPLLVICLCYLLIVVKVRAAGV--RVGCVRRRSERKVTRMVLVVVLVFAGCWLPFFTVNIVNLAVALPPASAGLYFFVVILSYANSCANPVLYGFLSDNFRQSFQKVLCLRKGSGAKDADATEQQQEATRPRTAA
>B2AR_PIG
MGQPGNRSVFLLAPNGSHAPDQ----DVPQERDEAWVVGMAIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGASHILMKMWTFGSFWCEFWISIDVLCVTASIETLCVIAVDRYLAITSPFYQCLLTKNKARVVILMVWVVSGLISFLPIKMHWYQAL-N--ACCDFF--------TNQPYAIASSIVSFYLPLVVMVFVYSRVFQVARRQGRSHRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHGIHDNL-IPKEVYILLNWVGYVNSAFNPLIYCRSP-DFRMAFQELLCLHRSSLKAYGNGCSDYTGEQSGCYLG
>FSHR_RAT
SDMMYNEFDYDLCNEVVDVTCSPKPDAFNPCEDIMGYNILRVLIWFISILAITGNTTVLVVLTTSQYKLTVPRFLMCNLAFADLCIGIYLLLIASVDIHTKSWQTG-AGCDAAGFFTVFASELSVYTLTAITLERWHTITHAMLECKVQLRHAASVMVLGWTFAFAAALFPIFGISSY----KVSICLPM-----IDSPLSQLYVMALLVLNVLAFVVICGCYTHIYLTVRNP------TIVSSSSDTKIAKRMATLIFTDFLCMAPISFFAISASLKVPLITVSKAKILLVLFYPINSCANPFLYAIFTKNFRRDFFILLSKFGCYEMQAQIYRTNFHARKSHCSSA
>D4DR_MOUSE
---MGNSSATEDGGLLAGRGP---ESLGTGAGLGGAGAAALVGGVLLIGLVLAGNSLVCVSVASERTLQTPTNYFIVSLAAADLLLAVLVLPLFVYSEVQGGWLLSPRLCDTLMAMDVMLCTASIFNLCAISVDRFVAVTVPL-RYNQQGQCQLLLIAATWLLSAAVASPVVCGLN---G-R--AVCCL---------ENRDYVVYSSVCSFFLPCPLMLLLYWATFRGLRRWPEPRRRGAKITGRERKAMRVLPVVVGAFLVCWTPFFVVHITRALCPACFVSPRLVSAVTWLGYVNSALNPIIYTIFNAEFRSVFRKTLRLRC-----------------------
>GALT_MOUSE
--------------------MADIQNISLDSPGSVGAVAVPVVFALIFLLGMVGNGLVLAVLLQPGPPGSTTDLFILNLAVADLCFILCCVPFQAAIYTLDAWLFGAFVCKTVHLLIYLTMYASSFTLAAVSVDRYLAVRHPLSRALRTPRNARAAVGLVWLLAALFSAPYLSYYG---GA--LELCVPA----WEDARRRALDVATFAAGYLLPVTVVSLAYGRTLCFLWAAGP-AAAAEARRRATGRAGRAMLAVAALYALCWGPHHALILCFWYGRFAPATYACRLASHCLAYANSCLNPLVYSLASRHFRARFRRLWPCGHRRHRHHHHRLHPASSGPAGYPGD
>AG2R_RAT
------MALNSSAEDGIKRIQDDC---PKAGRHSYIFVMIPTLYSIIFVVGIFGNSLVVIVIYFYMKLKTVASVFLLNLALADLCFLLTLPLWAVYTAMEYRWPFGNHLCKIASASVSFNLYASVFLLTCLSIDRYLAIVHPMSRLRRTMLVAKVTCIIIWLMAGLASLPAVIHRNVYF-TN--TVCAFH-YESRNSTLPIGLGLTKNILGFLFPFLIILTSYTLIWKALKKA----YEIQKNKPRNDDIFRIIMAIVLFFFFSWVPHQIFTFLDVLIQLGDIVDTAMPITICIAYFNNCLNPLFYGFLGKKFKKYFLQLLKYIPPKAKSHSSLSTRPSD-------N
>NMBR_RAT
PNLSLPTEASESELEPEVWENDFLPDSDGTTAELVIRCVIPSLYLIIISVGLLGNIMLVKIFLTNSTMRSVPNIFISNLAAGDLLLLLTCVPVDASRYFFDEWVFGKLGCKLIPAIQLTSVGVSVFTLTALSADRYRAIVNPMMQTSGVVLWTSLKAVGIWVVSVLLAVPEAVFSEVARSS--FTACIPY---QTDELHPKIHSVLIFLVYFLIPLVIISIYYYHIAKTLIRSLPGNEHTKKQMETRKRLAKIVLVFVGCFVFCWFPNHILYLYRSFNYKELGHMIVTLVARVLSFSNSCVNPFALYLLSESFRKHFNSQLCCGQKSYPERSTSYLMTSLKSNAKNVV
>NY1R_CANFA
TSFSQVENHSIFCNFSE-NSQFLAFESDDCHLPLAMIFTLALAYGAVIILGVTGNLALIMIILKQKEMRNVTNILIVNLSFSDLLVAIMCLPFTFVYTLMDHWVFGEAMCKLNPFVQCVSITVSIFSLVLIAVERHQLIINPR-GWRPNNRHAYVGIAVIWVLAVVSSLPFLIYQVLTDKD--KYVCFDK---FPSDSHRLSYTTLLLMLQYFGPLCFIFICYFKIYIRLKRRNMMMRDNKYRSSETKRINIMLLSIVVAFAVCWLPLTIFNTVFDWNHQICNHNLLFLLCHLTAMISTCVNPIFYGFLNKNFQRDLQFFFNFCDFRSRDDDY---TMHTDVSKTSLK
>A2AD_HUMAN
AAGPNASGAGERGSGGVANASGASWGPPRGQYSAGAVAGLAAVVGFLIVFTVVGNVLVVIAVLTSRALRAPQNLFLVSLASADILVATLVMPFSLANELMAYWYFGQVWCGVYLALDVLFCTSSIVHLCAISLDRYWSVTQAVYNLKRTPRRVKATIVAVWLISAVISFPPLVSLYR-------PQCGLN--------DETWYILSSCIGSFFAPCLIMGLVYARIYRVAKLRRRAVCRRKVAQAREKRFTFVLAVVMGVFVLCWFPFFFSYSLYGICREAQVPGPLFKFFFWIGYCNSSLNPVIYTVFNQDFRRSFKHILFRRRRRGFRQ-----------------
>V1AR_SHEEP
SSRWWPLDAGDANTSGDLAGLGEDGGPQADTRNEELAKLEIAVLAVIFVVAVLGNSSVLLALHRTPRKTSRMHLFIRHLSLADLAVAFFQVLPQLGWDITYRFRGPDGLCRVVKHMQVFAMFASAYMLVVMTADRYIAVCHPLKTLQQPARRSRLMIAAAWVLSFVLSTPQYFVFSMVETK--TYDCWAN---FIHPWGLPAYVTWMTGSVFVAPVVILGTCYGFICYHIWRKVLHVSSVKTISRAKIRTVKMTFVIVTAYIVCWAPFFIIQMWSAWDKNFESENPATAIPALLASLNSCCNPWIYMFFSGHLLQDCAQSFPCCQNVKRTFTREGSTSFTNNRSPTNS
>FML1_MOUSE
-----------MESNYSIHLNGSEVVVYDSTISRVLWILSMVVVSITFFLGVLGNGLVIWVAGFRMP-HTVTTIWYLNLALADFSFTATLPFLLVEMAMKEKWPFGWFLCKLVHIAVDVNLFGSVFLIAVIALDRCICVLHPVAQNHRTVSLARNVVVGSWIFALILTLPLFLFLTTVRVSWVEERLNTA------ITFVTTRGIIRFIVSFSLPMSFVAICYGLITTKI---------HKKAFVNSSRPFRVLTGVVASFFICWFPFQLVALLGTVWLKEKIIGRLVNPTSSLAFFNSCLNPILYVFMGQDFQERLIHSLSSRLQRALSEDSGHIAS----------
>OPSD_PETMA
MNGTEGENFYIPFSNKTGLARSPFEYPQYYLAEPWKYSVLAAYMFFLILVGFPVNFLTLFVTVQHKKLRTPLNYILLNLAVANLFMVLFGFTLTMYSSMNGYFVFGPTMCNFEGFFATLGGEMSLWSLVVLAIERYIVICKPMGNFRFGSTHAYMGVAFTWFMALSCAAPPLVGWSRYLPEGMQCSCGPDYYTLNPNFNNESFVIYMFLVHFIIPFIVIFFCYGRLLCTVKEAAAAQQESASTQKAEKEVTRMVVLMVIGFLVCWVPYASVAFYIFTHQGSDFGATFMTVPAFFAKTSALYNPIIYILMNKQFRNCMITTLCCGKNPLGDEDSGASVSSVSTSQVSPA
>B2AR_MACMU
MGQPGNGSAFLLAPNGSHAPDH----DVTQERDEAWVVGMGIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGAAHILMKMWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRAI-N--TCCDFF--------TNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQGRTLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNL-IPKEVYILLNWVGYVNSGFNPLIYCRSP-DFRIAFQELLCLRRSSLKACGNGYSGNTGEQSGYHLE
>AG2R_SHEEP
------MILNSSTEDGIKRIQDDC---PKAGRHNYIFIMIPTLYSIIFVVGLFGNSLVVIVIYFYMKLKTVASVFLLNLALADLCFLLTLPLWAVYTAMEYRWPFGNYLCKIASGSVSFNLYASVFLLTCLSIDRYLAIVHPMSRLRRTMLVAKVTCIIIWLLAGLASLPTIIHRNVFF-TN--TVCAFHVYESQNSTLPVGLGLTKNILGFLFPFLIILTSYTLIWKTLKKA----YEIQKNKPRKDDIFKIILAIVLFFFFSWVPHQIFTFMDVLIQLGDIVDTAMPITICLAYFNNCLNPPFYGFLGKKFKKYFLQLLKYIPPKAKSHSNLSTRPSE-------N
>AG2R_BOVIN
------MILNSSTEDGIKRIQDDC---PKAGRHNYIFIMIPTLYSIIFVVGIFGNSLVVIVIYFYMKLKTVASVFLLNLALADLCFLLTLPLWAVYTAMEYRWPFGNYLCKIASASVSFNLYASVFLLTCLSIDRYLAIVHPMSRLRRTMLVAKVTCIIIWLLAGLASLPTIIHRNVFF-TN--TVCAFH-YESQNSTLPVGLGLTKNILGFLFPFLIILTSYTLIWKTLKKA----YEIQKNKPRKDDIFKIILAIVLFFFFSWVPHQIFTFMDVLIQLGDIVDTAMPITICLAYFNNCLNPLFYGFLGKKFKKYFLQLLKYIPPKAKSHSNLSTRPSE-------N
>ACM2_HUMAN
--------------MNNSTNSSNNSLALTSPYKTFEVVFIVLVAGSLSLVTIIGNILVMVSIKVNRHLQTVNNYFLFSLACADLIIGVFSMNLYTLYTVIGYWPLGPVVCDLWLALDYVVSNASVMNLLIISFDRYFCVTKPLYPVKRTTKMAGMMIAAAWVLSFILWAPAILFWQFIV-VE-DGECYIQ------FFSNAAVTFGTAIAAFYLPVIIMTVLYWHISRASKSRVKMPAKKKPPPSREKKVTRTILAILLAFIITWAPYNVMVLINTFCAPC-IPNTVWTIGYWLCYINSTINPACYALCNATFKKTFKHLLMCHYKNIGATR----------------
>PF2R_HUMAN
--MSMNN---------SKQLVSPAAALLSNTTCQTENRLSVFFSVIFMTVGILSNSLAIAILMKAYQSKASFLLLASGLVITDFFGHLINGAIAVFVYASDKFDQSNVLCSIFGICMVFSGLCPLLLGSVMAIERCIGVTKPIHSTKITSKHVKMMLSGVCLFAVFIALLPILGHRDYKIQASRTWCFYN--TEDIKDWEDRFYLLLFSFLGLLALGVSLLCNAITGITLLRVKFKSQQHRQGRSHHLEMVIQLLAIMCVSCICWSPFLVTMANIGINGNHLETCETTLFALRMATWNQILDPWVYILLRKAVLKNLYKLASQCCGVHVISLHIWELKVAAISESPVA
>US28_HCMVA
----MTPTTTTAELTTEFDYDEDATPCVFTDVLNQSKPVTLFLYGVVFLFGSIGNFLVIFTITWRRRIQCSGDVYFINLAAADLLFVCTLPLWMQYLLDH--NSLASVPCTLLTACFYVAMFASLCFITEIALDRYYAIVYMR---YRPVKQACLFSIFWWIFAVIIAIPHFMVVTKK-----DNQCMTD-YDYLEVSYPIILNVELMLGAFVIPLSVISYCYYRISRIV---------AVSQSRHKGRIVRVLIAVVLVFIIFWLPYHLTLFVDTLKLLKRSLKRALILTESLAFCHCCLNPLLYVFVGTKFRQELHCLLAEFRQRLFSRDVSWYRSSPSRRETSSD
>NK4R_HUMAN
NLTSSPAPTASPSPAPSWTPSPRPGPAHPFLQPPWAVALWSLAYGAVVAVAVLGNLVVIWIVLAHKRMRTVTNSFLVNLAFADAAMAALNALVNFIYALHE-WYFGANYCRFQNFFPITAVFASIYSMTAIAVDRYMAIIDPL-KPRLSATATRIVIGSIWILAFLLAFPQCLYSK---GR---TLCYVQ--WPEGSRQHFTYHMIVIVLVYCFPLLIMGITYTIVGITLWGG--PCDKYQEQLKAKRKVVKMMIIVVVTFAICWLPYHIYFILTAIYQQLKYIQQVYLASFWLAMSSTMYNPIIYCCLNKRFRAGFKRAFRWCPFIHVSSY----ATRLHPMRQSSL
>MC3R_RAT
NSSCCPSSSYPTLPNLSQHPAAPSASNRSGSGFCEQVFIKPEVFLALGIVSLMENILVILAVVRNGNLHSPMYFFLLSLLQADMLVSLSNSLETIMIVVINSDQFIQHMDNIFDSMICISLVASICNLLAIAVDRYVTIFYALYHSIMTVRKALSLIVAIWVCCGICGVMFIVYYS----------------------EESKMVIVCLITMFFAMVLLMGTLYIHMFLFARLHIAAADGVAPQQHSCMKGAVTITILLGVFIFCWAPFFLHLVLIITCPTNICYTAHFNTYLVLIMCNSVIDPLIYAFRSLELRNTFKEILCGCNGMNVG------------------
>OPS1_HEMSA
TFGYPEGMTVADFVPDRVKHMVLDHWYNYPPVNPMWHYLLGVVYLFLGVISIAGNGLVIYLYMKSQALKTPANMLIVNLALSDLIMLTTNFPPFCYNCFSGGWMFSGTYCEIYAALGAITGVCSIWTLCMISFDRYNIICNGFNGPKLTQGKATFMCGLAWVISVGWSLPPFFGWGSYTLEGILDSCSYD--YFTRDMNTITYNICIFIFDFFLPASVIVFSYVFIVKAIFAHMNVRSNEAETQRAEIRIAKTALVNVSLWFICWTPYAAITIQGLLGNAEGITPLLTTLPALLAKSCSCYNPFVYAISHPKFRLAITQHLPWFCVHEKDPNDVEETQEKS-------
>AA3R_RAT
-----------------------MKANNTTTSALWLQITYITMEAAIGLCAVVGNMLVIWVVKLNRTLRTTTFYFIVSLALADIAVGVLVIPLAIAVSLE--VQMHFYACLFMSCVLLVFTHASIMSLLAIAVDRYLRVKLTVYRTVTTQRRIWLFLGLCWLVSFLVGLTPMFGWN-RKE-L-TLSCHFR-----SVVGLDYMVFFSFITWILIPLVVMCIIYLDIFYIIRNK--NFRETRAFYGREFKTAKSLFLVLFLFALCWLPLSIINFVSYFNVK--IPEIAMCLGILLSHANSMMNPIVYACKIKKFKETYFVILRACRLCQTSDSLDSN------------
>SSR4_RAT
GED----TTWTPGINASWAPDEEEDAVRSDGTGTAGMVTIQCIYALVCLVGLVGNALVIFVILRYAKMKTATNIYLLNLAVADELFMLSVPFVASAAALRH-WPFGAVLCRAVLSVDGLNMFTSVFCLTVLSVDRYVAVVHPLAATYRRPSVAKLINLGVWLASLLVTLPIAVFADTRPGGE-AVACNLH---WPHPAWSAVFVIYTFLLGFLLPVLAIGLCYLLIVGKMRAVALR-AGWQQRRRSEKKITRLVLMVVTVFVLCWMPFYVVQLLNLFVTS--LDATVNHVSLILSYANSCANPILYGFLSDNFRRSFQRVLCLRCCLLETTGGAEETALKSRGGPGCI
>NK1R_RAT
-----MDNVLPMDSDLFPNISTNTSESNQFVQPTWQIVLWAAAYTVIVVTSVVGNVVVIWIILAHKRMRTVTNYFLVNLAFAEACMAAFNTVVNFTYAVHNVWYYGLFYCKFHNFFPIAALFASIYSMTAVAFDRYMAIIHPL-QPRLSATATKVVIFVIWVLALLLAFPQGYYST---SR---VVCMIEWPEHPNRTYEKAYHICVTVLIYFLPLLVIGYAYTVVGITLWAS-IPSDRYHEQVSAKRKVVKMMIVVVCTFAICWLPFHVFFLLPYINPDLKFIQQVYLASMWLAMSSTMYNPIIYCCLNDRFRLGFKHAFRCCPFISAGDY----STRYLQTQSSVY
>OPSD_SARDI
MNGTEGPFFYIPMVNTTGIVRSPYEYPQYYLVNPAAYAILGAYMFFLIIVGFPVNFMTLYVTLEHKKLRTPLNYILLNLAVADLFMVIGGFTTTMYTSMHGYFVLGRLGCNLEGFFATLGGMISLWSLAVLAIERWVVVCKPISNFRFGENHAIMGVSLTWGMALACTVPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVLYMFFCHFTIPLTIIFFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVIIMVIGFLVCWLPYASVAWFIFTHQGSEFGPLFMTIPAFFAKSSSIYNPMIYICMNKQFRHCMITTLFCGKNPF---EGEEETEASSASSVSPA
>FSHR_HORSE
FDMMYSEFEYDLCNEVVDVTCSPKPDAFNPCEDIMGYDILRVLIWFISILAITGNIIVLVILITSQYKLTVPRFLMCNLAFADLCIGIYLLLIASVDIHTKSWQTG-AGCDAAGFFTVFASELSVYTLTAITLERWHTITHAMLECKVQLRHAASVMLVGWIFAFAVALLPIFGISTY----KVSICLPM-----IDSPLSQLYVMSLLVLNVLAFVVICGCYIHIYLTVRNP------NIVSSSSDTKIAKRMAILIFTDFLCMAPISFFAISASLKVPLITVSKSKILLVLFYPINSCANPFLYAIFTKNFRRDFFILLSKFGCYEMQAQLYRTISHPRNGHCPPT
>OPSG_CARAU
MNGTEGKNFYVPMSNRTGLVRSPFEYPQYYLAEPWQFKILALYLFFLMSMGLPINGLTLVVTAQHKKLRQPLNFILVNLAVAGTIMVCFGFTVTFYTAINGYFVLGPTGCAVEGFMATLGGEVALWSLVVLAIERYIVVCKPMGSFKFSSSHAFAGIAFTWVMALACAAPPLFGWSRYIPEGMQCSCGPDYYTLNPDYNNESYVIYMFVCHFILPVAVIFFTYGRLVCTVKAAAAQQQDSASTQKAEREVTKMVILMVFGFLIAWTPYATVAAWIFFNKGADFSAKFMAIPAFFSKSSALYNPVIYVLLNKQFRNCMLTTIFCGKNPLGDDESS--TSKTEVSSVSPA
>AA3R_RABIT
------------------------MPDNSTTLFLAIRASYIVFEIVIGVCAVVGNVLVIWVIKLNPSLKTTTFYFIFSLALADIAVGFLVMPLAIVISLG--ITIGFYSCLVMSCLLLVFTHASIMSLLAIAVDRYLRVKLTVYRRVTTQRRIWLALGLCWVVSLLVGFTPMFGWN-MKE-S-DFQCKFD-----SVIPMEYMVFFSFFTWILIPLLLMCALYVYIFYIIRNK--SFKETGAFYRREFKTAKSLFLVLALFAGCWLPLSIINCVTYFKCK--VPDVVLLVGILLSHANSMMNPIVYACKIQKFKETYLLIFKARVTCQPSDSLDPS------------
>5H1D_HUMAN
SPLNQSAEGLPQEASNRSLNATETSEAWDPRTLQALKISLAVVLSVITLATVLSNAFVLTTILLTRKLHTPANYLIGSLATTDLLVSILVMPISIAYTITHTWNFGQILCDIWLSSDITCCTASILHLCVIALDRYWAITDALYSKRRTAGHAATMIAIVWAISICISIPPLFWR----Q-E--SDCLVN-------TSQISYTIYSTCGAFYIPSVLLIILYGRIYRAARNRKLALERKRISAARERKATKILGIILGAFIICWLPFFVVSLVLPICRDSWIHPALFDFFTWLGYLNSLINPIIYTVFNEEFRQAFQKIVPFRKAS---------------------
>CCR4_MACMU
IYTSDNYTEEMG-SGDYDSIKEPC---FREENAHFNRIFLPTIYSIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVADLLFVITLPFWAVDAVAN--WYFGNFLCKAVHVIYTVNLYSSVLILAFISLDRYLAIVHATSQKPRKLLAEKVVYVGVWIPALLLTIPDFIFASVSE-DD-RYICDRF---YPNDLWVVVFQFQHIMVGLILPGIDILSCYCIIISKL---------SHSKGHQKRKALKTTVILILAFFACWLPYYIGISIDSFILLENTVHKWISITEALAFFHCCLNPILYAFLGAKFKTSAQHALTSVSRGSSLKILSKG----------GH
>OPSP_CHICK
------MSSNSSQAPPNGTPGPFDGPQWPYQAPQSTYVGVAVLMGTVVACASVVNGLVIVVSICYKKLRSPLNYILVNLAVADLLVTLCGSSVSLSNNINGFFVFGRRMCELEGFMVSLTGIVGLWSLAILALERYVVVCKPLGDFQFQRRHAVSGCAFTWGWALLWSAPPLLGWSSYVPEGLRTSCGPN--WYTGGSNNNSYILSLFVTCFVLPLSLILFSYTNLLLTLRAAAAQQKEADTTQRAEREVTRMVIVMVMAFLLCWLPYSTFALVVATHKGIIIQPVLASLPSYFSKTATVYNPIIYVFMNKQFQSCLLEMLCCGYQPQRTGKASPGVTAAGLRNKVMP
>O7C1_HUMAN
-------------METGNQTHAQEFLLLGFSATSEIQFILFGLFLSMYLVTFTGNLLIILAICSDSHLHTPMYFFLSNLSFADLCFTSTTVPKMLLNILTQNKFITYAGCLSQIFFFTSFGCLDNLLLTVMAYDRFVAVCHPLYTVIMNPQLCGLLVLGSWCISVMGSLLETLTVLRL---SPHFFCDLLKLACSDTFINNVVIYFATGVLGVISFTGIFFSYYKIVFSI--------LRISSAGRKHKAFSTCGSHLSVVTLFYGTGFGVYLSSAATP---SSRTSLVASVMYTMVTPMLNPFIYSLRNTDMKRALGRLLSRATFFNGDITAGLS------------
>OPR._PIG
GSPLQGNLSLLSPNHSLLPPHLLLNASHGAFLPLGLKVTIVGLYLAVCVGGLLGNCLVMYVILRHTKMKTATNIYIFNLALADTAVLLTLPFQGTDVLLGF-WPFGNALCKAVIAIDYYNMFTSAFTLTAMSVDRYVAICHPIALDVRTSSKAQAVNVAIWALASIVGVPVAIMGSAQVEE---IECLVE-IPAPQDYWGPVFAVCIFLFSFVIPVLIISVCYSLMVRRLRGVRLL-SGSREKDRNLRRITRLVLVVVAVFVGCWTPVQVFVLVQGLGVQPETAVAVLRFCTALGYVNSCLNPILYAFLDENFKACFRKFCCAPTRRREMQVSDRVALACKTSETVPR
>OPSD_RAT
MNGTEGPNFYVPFSNITGVVRSPFEQPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIGLWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIFFLICWLPYASVAMYIFTHQGSNFGPIFMTLPAFFAKTASIYNPIIYIMMNKQFRNCMLTSLCCGKNPLGDDEASATASKTETSQVAPA
>P2UR_HUMAN
----MAADLGPWNDTINGTWDGDELGYRCRFNEDFKYVLLPVSYGVVCVLGLCLNAVALYIFLCRLKTWNASTTYMFHLAVSDALYAASLPLLVYYYARGDHWPFSTVLCKLVRFLFYTNLYCSILFLTCISVHRCLGVLRPLSLRWGRARYARRVAGAVWVLVLACQAPVLYFVTTS-----RVTCHDT-SAPELFSRFVAYSSVMLGLLFAVPFAVILVCYVLMARRLLKP---YGTSGGLPRAKRKSVRTIAVVLAVFALCFLPFHVTRTLYYSFRSLNAINMAYKVTRPLASANSCLDPVLYFLAGQRLVRFARDAKPPTGPSPATPARRRLTDMQRIGDVLGS
>D1DR_OREMO
-MEIFTTTRGTSAGPEPAPGGHGGTDSPR-TSDLSLRALTGCVLCILIVSTLLGNALVCAAVIKFRHRSKVTNAFVISLAVSDLFVAVLVMPWRAVSEVAGVWLFG-AFCDTWVAFDIMCSTASILHLCIISMDRYWAISSPFYERRMTPRFGCVMIGVAWTLSVLISFIPVQLNWHARR--DPGDCNAS--------LNRTYAISSSLISFYIPVLIMVGTYTRIFRIGRTQPALESSLKTSFRRETKVLKTLSVIMGVFVFCWLPFFVLNCMVPFCRLECVSDTTFSVFVWFGWANSSLNPVIYAFNA-DFRKAFSTILGCSRYCRTSAVEAVDYHHDTTLQK---
>OPSP_PETMA
HHSLPSALPSATGGNGTVATMHNPFERPLEGIAPWNFTMLAALMGTITALSLGENFAVIVVTARFRQLRQPLNYVLVNLAAADLLVSAIGGSVSFFTNIKGYFFLGVHACVLEGFAVTYFGVVALWSLALLAFERYFVICRPLGNFRLQSKHAVLGLAVVWVFSLACTLPPVLGWSSYRPSMIGTTCEPN--WYSGELHDHTFILMFFSTCFIFPLAVIFFSYGKLIQKLKKASETQRGLESTRRAEQQVTRMVVVMILAFLVCWMPYATFSIVVTACPTI---PLLAAVPAFFSKTATVYNPVIYIFMNKQFRDCFVQVLPCKGLKKVSATQTAGASVNTQSPGNRH
>DADR_RAT
---------------MAPNTSTMDEAGLPAERDFSFRILTACFLSLLILSTLLGNTLVCAAVIRFRHRSKVTNFFVISLAVSDLLVAVLVMPWKAVAEIAGFWPLG-PFCNIWVAFDIMCSTASILNLCVISVDRYWAISSPFYERKMTPKAAFILISVAWTLSVLISFIPVQLSWHKAGN-EDDNCDTR--------LSRTYAISSSLISFYIPVAIMIVTYTSIYRIAQKQVECESSFKMSFKRETKVLKTLSVIMGVFVCCWLPFFISNCMVPFCGSECIDSITFDVFVWFGWANSSLNPIIYAFNA-DFQKAFSTLLGCYRLCPTTNNAIETAVVFSSH-----
>GALS_RAT
------------MNGSGSQGAENTSQEGGSGGWQPEAVLVPLFFALIFLVGTVGNALVLAVLLRGGQAVSTTNLFILNLGVADLCFILCCVPFQATIYTLDDWVFGSLLCKAVHFLIFLTMHASSFTLAAVSLDRYLAIRYPLSRELRTPRNALAAIGLIWGLALLFSGPYLSYYR---AN--LTVCHPA----WSAPRRRAMDLCTFVFSYLLPVLVLSLTYARTLRYLWRTDP-VTAGSGSQRAKRKVTRMIIIVAVLFCLCWMPHHALILCVWFGRFPRATYALRILSHLVSYANSCVNPIVYALVSKHFRKGFRKICAGLLRPAPRRASGRVHSGSMLEQESTD
>5H7_MOUSE
SSWMPHLLSGFPEVTASPAPTNVSGCGEQINYGRVEKVVIGSILTLITLLTIAGNCLVVISVCFVKNVRQPSNYLIVSLALADLSVAVAVMPFVSVTDLIGGWIFGHFFCNVFIAMDVMCCTASIMTLCVISIDRYLGITRPLYPVRQNGKCMAKMILSVWPLSASITLPPLFGWAQ--N-D--KVCLIS--------QDFGYTIYSTAVAFYIPMSVMLFMYYQIYKAARKSSRLERKNISSFKREQKAATTLGIIVGAFTVCWLPFFLLSTARPFICGTCIPLWVERTCLWLGYANSLINPFIYSFFNRDLRTTYRSLLQCQYRNINRKLSAAGAERPERSEFVLQ
>EBI2_HUMAN
------MDIQMANNFTPPSATPQGNDCDLYAHHSTARIVMPLHYSLVFIIGLVGNLLALVVIVQNRKKINSTTLYSTNLVISDILFTTALPTRIAYYAMGFDWRIGDALCRITALVFYINTYAGVNFMTCLSIDRFIAVVHPLYNKIKRIEHAKGVCIFVWILVFAQTLPLLINPMS--KQEERITCMEY--NFEETKSLPWILLGACFIGYVLPLIIILICYSQICCKLFRTAK-QNPLTEKSGVNKKALNTIILIIVVFVLCFTPYHVAIIQHMIKKLRHSFQISLHFTVCLMNFNCCMDPFIYFFACKGYKRKVMRMLKRQVSVSISSAVKSAMTETQMMIHSKS
>TSHR_CANFA
LQAFDSHYDYTVCGGNEDMVCTPKSDEFNPCEDIMGYKFLRIVVWFVSLLALLGNVFVLIVLLTSHYKLTVPRFLMCNLAFADFCMGMYLLLIASVDLYTHSWQTG-PGCNTAGFFTVFASELSVYTLTVITLERWYAITFAMLDRKIRLRHAYAIMVGGWVCCFLLALLPLVGISSY----KVSICLPM-----TETPLALAYIILVLLLNIVAFIIVCSCYVKIYITVRNP------QYNPGDKDTKIAKRMAVLIFTDFMCMAPISFYALSALMNKPLITVTNSKILLVLFYPLNSCANPFLYAIFTKAFQRDVFILLSKFGICKRQAQAYRGSAGIQIQKVTRD
>NY4R_RAT
HLMASLSPAFLQGKNGTNPLDSLYNLSDGCQDSADLLAFIITTYSVETVLGVLGNLCLIFVTTRQKEKSNVTNLLIANLAFSDFLMCLICQPLTVTYTIMDYWIFGEVLCKMLTFIQCMSVTVSILSLVLVALERHQLIINPT-GWKPSISQAYLGIVVIWFISCFLSLPFLANSILNDED--KVVCFVS---WSSDHHRLIYTTFLLLFQYCVPLAFILVCYMRIYQRLQRQRA-THTCSSRVGQMKRINGMLMAMVTAFAVLWLPLHVFNTLEDWYQEACHGNLIFLMCHLFAMASTCVNPFIYGFLNINFKKDIKALVLTCRCRPPQGEP---TVHTDLSKGSMR
>CKR5_PAPHA
----MDYQVSSPTYDIDYYTSEPC---QKINVKQIAARLLPPLYSLVFIFGFVGNILVVLILINCKRLKSMTDIYLLNLAISDLLFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWVVAVFASLPGIIFTRSQR-EG-HYTCSSHFPYSQYQFWKNFQTLKIVILGLVLPLLVMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCSIFA-------SSVY
>OLF2_RAT
------------MESGNSTRRFSSFFLLGFTENPQLHFLIFALFLSMYLVTVLGNLLIIMAIITQSHLHTPMYFFLANLSFVDICFTSTTIPKMLVNIYTQSKSITYEDCISQMCVFLVFAELGNFLLAVMAYDRYVA-CHPLYTVIVNHRLCILLLLLSWVISIFHAFIQSLIVLQL---TPHFFCELNQLTCSDNFPSHLIMNLVPVMLAAISFSGILYSYFKIVSSI--------HSISTVQGKYKAFSTCASHLSIVSLFYSTGLGVYVSSAVVQ---SSHSAASASVMYTVVTPMLNPFIYSLRNKDVKRALERLLEGNCKVHHWTG----------------
>HH2R_HUMAN
-------------------MAPNGTASSFCLDSTACKITITVVLAVLILITVAGNVVVCLAVGLNRRLRNLTNCFIVSLAITDLLLGLLVLPFSAIYQLSCKWSFGKVFCNIYTSLDVMLCTASILNLFMISLDRYCAVMDPLYPVLVTPVRVAISLVLIWVISITLSFLSIHLGWN--RN-TTSKCKVQ--------VNEVYGLVDGLVTFYLPLLIMCITYYRIFKVARDQR--ISSWKAATIREHKATVTLAAVMGAFIICWFPYFTAFVYRGLRGDD-INEVLEAIVLWLGYANSALNPILYAALNRDFRTGYQQLFCCRLANRNSHKTSLRRTQSREPRQ---
>CKR5_HYLLE
----MDYQVSSPTYDIDYDTSEPC---QKINVKQIAARLLPPLYSLVFIFGFVGNMLVILVLINCKRLKSMTDIYLLNLAISDLFFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWVVAVFASLPGIIFTRSQK-EG-HYTCSSHFPYSQYQFWKNFQTLKIVILGLVLPLLVMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKHFCKCCSIFA-------SSVY
>B3AR_HUMAN
MAPWPHENSSLAPWPDLPTLAPNTANTSGLPGVPWEAALAGALLALAVLATVGGNLLVIVAIAWTPRLQTMTNVFVTSLAAADLVMGLLVVPPAATLALTGHWPLGATGCELWTSVDVLCVTASIETLCALAVDRYLAVTNPLYGALVTKRCARTAVVLVWVVSAAVSFAPIMSQWWRVQ-R--RCCAFA--------SNMPYVLLSSSVSFYLPLLVMLFVYARVFVVATRQGVPRRPARLLPLREHRALCTLGLIMGTFTLCWLPFFLANVLRALGGPS-VPGPAFLALNWLGYANSAFNPLIYCRSP-DFRSAFRRLLCRCGRRLPPEPCAAAGVPAARSSPAQP
>NY6R_MOUSE
--MEVLTNQPTPNKTSGKSNNSAFFYFESCQPPFLAILLLLIAYTVILIMGIFGNLSLIIIIFKKQRAQNVTNILIANLSLSDILVCVMCIPFTVIYTLMDHWVFGNTMCKLTSYVQSVSVSVSIFSLVLIAIERYQLIVNPR-GWKPRVAHAYWGIILIWLISLTLSIPLFLSYHLTNTH--QVACVEI---WPSKLNQLLFSTSLFMLQYFVPLGFILICYLKIVLCLRKRTRQRKENKSRLNENKRVNVMLISIVVTFGACWLPLNIFNVIFDWYHEMCHHDLVFVVCHLIAMVSTCINPLFYGFLNKNFQKDLMMLIHHCWCGEPQESY---TMHTDESKGSLK
>OLF4_RAT
-------------MTGNNQTLILEFLLLGLPIPSEYHLLFYALFLAMYLTIILGNLLIIVLVRLDSHLHMPMYLFLSNLSFSDLCFSSVTMPKLLQNMQSQVPSISYTGCLTQLYFFMVFGDMESFLLVVMAYDRYVAICFPLYTTIMSTKFCASLVLLLWMLTMTHALLHTLLIARL---SLHFFCDISKLSCSDIYVNELMIYILGGLIIIIPFLLIVMSYVRIFFSI--------LKFPSIQDIYKVFSTCGSHLSVVTLFYGTIFGIYLCPSGNN---STVKEIAMAMMYTVVTPMLNPFIYSLRNRDMKRALIRVICTKKISL--------------------
>OPSD_ZEUFA
MNGTEGPDFYVPMVNTTGIVRSPYDYPQYYLVNPAAFSMLAAYMFFLILVGFPVNFLTLYVTMEHKKLRTPLNYILLNLAVANLFMVIGGFTTTMYTSMHGYFVLGRTGCNLEGFFATLGGEIALWSLVVLAVERWVVVCKPISNFRFGENHAVMGVSFTWLMACACSVPPLFGWSRYIPEGMQCSCGIDYYTRAPGYNNESFVIYMFVCHFSIPLTIIFFCYGRLLCAVKDAAAAQQESETTQRAEREVSRMVVIMVIGFLICWLPYASVAWFIFTHQGSEFGPVFMTIPAFFAKSSAIYNPMIYICMNKQFRHCMITTLCCGKNPFEEEEGASTASSVSSSHVSPA
>OPSD_CYPCA
MNGTEGPMFYVPMSNATGVVKSPYDYPQYYLVAPWAYGCLAAYMFFLIITGFPINFLTLYVTIEHKKLRTPLNYILLNLAISDLFMVFGGFTTTMYTSLHGYFVFGRIGCNLEGFFATLGGEMGLWSLVVLAFERWMVVCKPVSNFRFGENHAIMGVVFTWFMACTCAVPPLVGWSRYIPEGMQCSCGVDYYTRAPGYNNESFVIYMFLVHFIIPLIVIFFCYGRLVCTVKDAAAQQQESETTQRAEREVTRMVVIMVIGFLICWIPYASVAWYIFTHQGSEFGPVFMTVPAFFAKSAAVYNPCIYICMNKQFRHCMITTLCCGKNPFEEEEGASTASSVSSSSVSPA
>A2AR_LABOS
-----MDPLNATGMDAFTAIHLNASWSADSGYSLAAIASIAALVSFLILFTVVGNILVVIAVLTSRALKAPQNLFLVSLATADILVATLVMPFSLANELMGYWYFGKVWCGIYLALDVLFCTSSIVHLCAISLDRYWSVTQAVYNLKRTPKRVKCIIVIVWLISAFISSPPLLSIDS--I-S--PQCMLN--------DDTWYILSSSMASFFAPCLIMILVYIRIYQVAKTRKRRIAEKKVSQAREKRFTFVLAVVMGVFVVCWFPFFFSYSLHAVCRDYKIPDTLFK-FFWIGYCNSSLNPAIYTIFNRDFRRAFQKILCKSWKKSF-------------------
>D3DR_HUMAN
--------MASLSQLSSHLNYTCGAENSTGASQARPHAYYALSYCALILAIVFGNGLVCMAVLKERALQTTTNYLVVSLAVADLLVATLVMPWVVYLEVTGGWNFSRICCDVFVTLDVMMCTASILNLCAISIDRYTAVVMPVGTGQSSCRRVALMITAVWVLAFAVSCPLLFGFN---D-P--TVCSI---------SNPDFVIYSSVVSFYLPFGVTVLVYARIYVVLKQRTSLPLQPRGVPLREKKATQMVAIVLGAFIVCWLPFFLTHVLNTHCQTCHVSPELYSATTWLGYVNSALNPVIYTTFNIEFRKAFLKILSC-------------------------
>D2DR_FUGRU
-----MDVFTQYAYNDSIFDNGTWSANETTKDETHPYNYYAMLLTLLIFVIVFGNVLVCMAVSREKALQTTTNYLIVSLAVADLLVATLVMPWVVYLEVVGEWRFSKIHCDIFVTLDVMMCTASILNLCAISIDRYTAVAMPMNTRYSSRRRVTVMISVVWVLSFAISCPLLFGLN---T-R--SLCFI---------ANPAFVVYSSIVSFYVPFIVTLLVYVQIYVVLRKRQTSLSKRKISQQKEKKATQMLAIVLGVFIICWLPFFITHILNTHCTRCKVPAEMYNAFTWLGYVNSAVNPIIYTTFNVEFRKAFIKILHC-------------------------
>P2Y4_HUMAN
ASTESSLLRSLGLSPGPGS---SEVELDCWFDEDFKFILLPVSYAVVFVLGLGLNAPTLWLFIFRLRPWDATATYMFHLALSDTLYVLSLPTLIYYYAAHNHWPFGTEICKFVRFLFYWNLYCSVLFLTCISVHRYLGICHPLALRWGRPRLAGLLCLAVWLVVAGCLVPNLFFVTTS-----TVLCHDT-TRPEEFDHYVHFSSAVMGLLFGVPCLVTLVCYGLMARRLYQP----LPGSAQSSSRLRSLRTIAVVLTVFAVCFVPFHITRTIYYLARLLNIVNVVYKVTRPLASANSCLDPVLYLLTGDKYRRQLRQLCGGGKPQPRTAASSLASSCRWAATPQDS
>GALR_MOUSE
----MELAMVNLSEGNGSDPEPPAPESRPLFGIGVENFITLVVFGLIFAMGVLGNSLVITVLARSKPPRSTTNLFILNLSIADLAYLLFCIPFQATVYALPTWVLGAFICKFIHYFFTVSMLVSIFTLAAMSVDRYVAIVHSRSSSLRVSRNALLGVGFIWALSIAMASPVAYHQRLF-SN--QTFCWEQ---WPNKLHKKAYVVCTFVFGYLLPLLLICFCYAKVLNHLHKKLK--NMSKKSEASKKKTAQTVLVVVVVFGISWLPHHVVHLWAEFGAFPPASFFFRITAHCLAYSNSSVNPIIYAFLSENFRKAYKQVFKCHVCDESPRSETKEPPSTNCTHV---
>OAR1_LOCMI
SSAAEEPQDALVGGDACGGRRPPSVLGVRLAVPEWEVAVTAVSLSLIILITIVGNVLVVLSVFTYKPLRIVQNFFIVSLAVADLTVAVLVMPFNVAYSLIQRWVFGIVVCKMWLTCDVLCCTASILNLCAIALDRYWAITDPIYAQKRTLRRVLAMIAGVWLLSGVISSPPLIGWNDW-N-D--TPCQLT--------EEQGYVIYSSLGSFFIPLFIMTIVYVEIFIATKRRPVYEEKQRISLSKERRAARTLGIIMGVFVVCWLPFFLMYVIVPFCNPSKPSPKLVNFITWLGYINSALNPIIYTIFNLDFRRAFKKLLHFKT-----------------------
>ML1A_SHEEP
GSPGGTPKGNGSSALLNVSQAAPGAGDGVRPRPSWLAATLASILIFTIVVDIVGNLLVVLSVYRNKKLRNAGNVFVVSLAVADLLVAVYPYPLALASIVNNGWSLSSLHCQLSGFLMGLSVIGSVFSITGIAINRYCCICHSLYGKLYSGTNSLCYVFLIWTLTLVAIVPNLCVGT-LQYDP-IYSCTFT------QSVSSAYTIAVVVFHFIVPMLVVVFCYLRIWALVLQVWK-PDNKPKLKPQDFRNFVTMFVVFVLFAICWAPLNFIGLVVASDPASRIPEWLFVASYYMAYFNSCLNAIIYGLLNQNFRQEYRKIIVSLCTTKMFFVDSSNRKPSPLIANHNL
>OPSD_SOLSO
MNGTEGPYFYIPMLNTTGIVRSPYEYPQYYLVNPAAYAALCAYMFLLILLGFPINFLTLYVTIEHKKLRTPLNYILLNLAVANLFMVFGGFTTTMYTSMHGYFVLGRLGCNLEGFFATLGGEIGLWSLVVLAVERWMVVCKPISNFRFTENHAIMGLGFTWFAASACAVPPLVGWSRYIPEGMQCSCGVDYYTRAEGFNNESFVVYMFVCHFLIPLIVVFFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVVIMVIAFLICWCPYAGVAWYIFSNQGSEFGPLFMTIPAFFAKSSSIYNPLIYIFMNKQFRHCMITTLCCGKNPFEEEEGSTTSASSSSVSPAA-
>ML1A_MOUSE
-------MKGNVSELLNATQQAPGGGEGGRPRPSWLASTLAFILIFTIVVDILGNLLVILSVYRNKKLRNSGNIFVVSLAVADLVVAVYPYPLVLTSILNNGWNLGYLHCQVSAFLMGLSVIGSIFNITGIAMNRYCYICHSLYDKIYSNKNSLCYVFLIWMLTLIAIMPNLQTGT-LQYDP-IYSCTFT------QSVSSAYTIAVVVFHFIVPMIIVIFCYLRIWVLVLQVRR-PDNKPKLKPQDFRNFVTMFVVFVLFAICWAPLNLIGLIVASDPATRIPEWLFVASYYLAYFNSCLNAIIYGLLNQNFRKEYKKIIVSLCTAKMFFVESSNCKPSPLIPNNNL
>OPSD_POERE
MNGTEGPYFYVPMVNTTGIVRSPYEYPQYYLVSPAAYACLGAYMFFLILVGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTIYTSMHGYFVLGRLGCNLEGYFATLGGEIGLWSLVVLAVERWLVVCKPISNFRFSENHAIMGLVFTWIMANSCAAPPLLGWSRYIPEGMQCSCGVDYYTRAEGFNNESFVIYMFICHFCIPLIVVFFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVVIMVIGFLVCWIPYASVAWYIFTHQGSEFGPLFMTVPAFFAKSASIYNPLIYICMNKQFRHCMITTLCCGKNPFEEEEGASTASSVSSSSVSPA
>PE22_HUMAN
----------------MGNASNDSQSEDCETRQWLPPGESPAISSVMFSAGVLGNLIALALLARRWRSLSLFHVLVTELVFTDLLGTCLISPVVLASYARNQLAPESRACTYFAFAMTFFSLATMLMLFAMALERYLSIGHPYYQRRVSASGGLAVLPVIYAVSLLFCSLPLLDYGQYVQYCPGTWCFIR--------HGRTAYLQLYATLLLLLIVSVLACNFSVILNLIRMGGPRRGERVSMAEETDHLILLAIMTITFAVCSLPFTIFAYMNETSS---RKEKWDLQALRFLSINSIIDPWVFAILRPPVLRLMRSVLCCRISLRTQDATQTSSKQADL------
>CKR3_HUMAN
MTTSLDTVETFGTTSYYDDVGLLC---EKADTRALMAQFVPPLYSLVFTVGLLGNVVVVMILIKYRRLRIMTNIYLLNLAISDLLFLVTLPFWIHYVRGHN-WVFGHGMCKLLSGFYHTGLYSEIFFIILLTIDRYLAIVHAVALRARTVTFGVITSIVTWGLAVLAALPEFIFYETEE-LF-ETLCSALYPEDTVYSWRHFHTLRMTIFCLVLPLLVMAICYTGIIKTL---------LRCPSKKKYKAIRLIFVIMAVFFIFWTPYNVAILLSSYQSILKHLDLVMLVTEVIAYSHCCMNPVIYAFVGERFRKYLRHFFHRHLLMHLGRYIPFLT-------SSVS
>GPRA_HUMAN
TPANQSAEASAGNGSVAGADAPAVTPFQSLQLVHQLKGLIVLLYSVVVVVGLVGNCLLVLVIARVPRLHNVTNFLIGNLALSDVLMCTACVPLTLAYAFEPRWVFGGGLCHLVFFLQPVTVYVSVFTLTTIAVDRYVVLVHPL--RRASRCASAYAVLAIWALSAVLALPPAVHTYHVEHD--VRLCEEF--WGSQERQRQLYAWGLLLVTYLLPLLVILLSYVRVSVKLRNRPGC-SQADWDRARRRRTFCLLVVVVVVFAVCWLPLHVFNLLRDLDPHAYAFGLVQLLCHWLAMSSACYNPFIYAWLHDSFREELRKLLVAWPRKIAPHGQNMT------------
>OPSD_LIZAU
MNGTEGPYFYIPMVNTTGIVRSPYEYPQYYLVNPAAYAALGAYMFLLILIGFPVNFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVLGRLGCNLEGFFATLGGEIALWSLVVLAVERWMVVCKPISNFRFGEDHAIMGLAFTWVMAAACAVPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFVCHFLIPLVVVFFCYGRLLCAVKEAAAAQQESETTQRAEREVSRMVVIMVVAFLVCWCPYAGVAWYIFTHQGSEFGPLFMTFPAFFAKSSSIYNPMIYICMNKQFRQCMITTLCCGKNPFEEEEGASTSVSSSSVSPAA-
>PF2R_BOVIN
--MSTNS---------SIQPVSPESELLSNTTCQLEEDLSISFSIIFMTVGILSNSLAIAILMKAYQYKSSFLLLASALVITDFFGHLINGTIAVFVYASDKFDKSNILCSIFGICMVFSGLCPLFLGSLMAIERCIGVTKPIHSTKITTKHVKMMLSGVCFFAVFVALLPILGHRDYKIQASRTWCFYK--TDEIKDWEDRFYLLLFAFLGLLALGISFVCNAITGISLLKVKFRSQQHRQGRSHHFEMVIQLLGIMCVSCICWSPFLVTMASIGMNIQDKDSCERTLFTLRMATWNQILDPWVYILLRKAVLRNLYVCTRRCCGVHVISLHVWELKVAAISDLPVT
>GALR_RAT
----MELAPVNLSEGNGSDPEPP-AEPRPLFGIGVENFITLVVFGLIFAMGVLGNSLVITVLARSKPPRSTTNLFILNLSIADLAYLLFCIPFQATVYALPTWVLGAFICKFIHYFFTVSMLVSIFTLAAMSVDRYVAIVHSRSSSLRVSRNALLGVGFIWALSIAMASPVAYYQRLF-SN--QTFCWEH---WPNQLHKKAYVVCTFVFGYLLPLLLICFCYAKVLNHLHKKLK--NMSKKSEASKKKTAQTVLVVVVVFGISWLPHHVIHLWAEFGAFPPASFFFRITAHCLAYSNSSVNPIIYAFLSENFRKAYKQVFKCRVCNESPHGDAKEPPSTNCTHV---
>OPSU_BRARE
MNGTEGPAFYVPMSNATGVVRSPYEYPQYYLVAPWAYGFVAAYMFFLIITGFPVNFLTLYVTIEHKKLRTPLNYILLNLAIADLFMVFGGFTTTMYTSLHGYFVFGRLGCNLEGFFATLGGEMGLKSLVVLAIERWMVVCKPVSNFRFGENHAIMGVAFTWVMACSCAVPPLVGWSRYIPEGMQCSCGVDYYTRTPGVNNESFVIYMFIVHFFIPLIVIFFCYGRLVCTVKEAARQQQESETTQRAEREVTRMVIIMVIAFLICWLPYAGVAWYIFTHQGSEFGPVFMTLPAFFAKTSAVYNPCIYICMNKQFRHCMITTLCCGKNPFEEEEGASTASSVSSSSVSPA
>GPR6_RAT
GPPAASAALGGGGGPNGSLELSSQLPAGPSGLLLSAVNPWDVLLCVSGTVIAGENALVVALIASTPALRTPMFVLVGSLATADLLAG-CGLILHFVFQ-Y--VVPSETVSLLMVGFLVASFAASVSSLLAITVDRYLSLYNALYYSRRTLLGVHLLLAATWTVSLGLGLLPVLGWNCL-A--DRASCSVV-------RPLTRSHVALLSTSFFVVFGIMLHLYVRICQVVWRHIALHCLAPPHLAATRKGVGTLAVVLGTFGASWLPFAIYCVVGSQ----EDPAIYTYATLLPATYNSMINPIIYAFRNQEIQRALWLLFCGCFQSKVPFRSRSP------------
>ACM3_CHICK
DSPETTESFPFSTVETTNSSLNATIKDPLGGHAVWQVVLIAFLTGIIALVTIIGNILVIVSFKVNKQLKTVNNYFLLSLACADLIIGVISMNLFTTYIIMGHWALGNLACDLWLSIDYVASNASVMNLLVISFDRYFSITRPLYRAKRTTKRAGVMIGLAWIISFVLWAPAILFWQYFV-VP-LDECFIQ------FLSEPIITFGTAIAAFYLPVTIMSILYWRIYKETEKRKTRTKRKRMSLIKEKKAAQTLSAILFAFIITWTPYNIMVLVNTFCDC--VPKTVWNLGYWLCYINSTVNPVCYALCNKMFRNTFKMLLLCQCDKRKRRKQQYQHKRIPREAS---
>NY5R_CANFA
-------------NTAATRNSDFPVWDDYKSSVDDLQYFLIGLYTFVSLLGFMGNLLILMALMRKRNQKTMVNFLIGNLAFSDILVVLFCSPFTLTSVLLDQWMFGKVMCHIMPFLQCVSVLVSTLILISIAIVRYHMIKHPI-SNNLTANHGYFLIATVWTLGFAICSPLPVFHSLVESS--RYLCVES---WPSDSYRIAFTISLLLVQYILPLVCLTVSHTSVCRSISCGVHDNRSIMRIKKRSRSVFYRLTILILVFAVSWMPLHLFHVVTDFNDNLRHFKLVYCICHLLGMMSCCLNPILYGFLNNGIKADLISLIQCLHMS---------------------
>MC5R_SHEEP
-MNSSFHLHFLDLGLNATEGNLSGLSVRNASSPCEDMGIAVEVFLALGLISLLENILVIGAIVRNRNLHIPMYFFVGSLAVADMLVSLSNFWETITIYLLTNDASVRHLDNVFDSMICISVVASMCSLLAIAVDRYVTIFCRLYQRIMTGRRSGAIIAGIWAFCTSCGTVFIVYYY----------------------EESTYVVVCLIAMFLTMLLLMASLYTHMFLLARTH-RIPGHSSVRQRTGVKGAITLAMLLGVFIICWAPFFLHLILMISCPQNSCFMSHFNMYLILIMCNSVIDPLIYAFRSQEMRKTFKEIVCFQGFRTPCRFPSTY------------
>V2R_RAT
MLLVSTVSAVPGLFSPPSSPSNSSQEELLDDRDPLLVRAELALLSTIFVAVALSNGLVLGALIRRGRRWAPMHVFISHLCLADLAVALFQVLPQLAWDATDRFHGPDALCRAVKYLQMVGMYASSYMILAMTLDRHRAICRPMAYRHGGGARWNRPVLVAWAFSLLLSLPQLFIFAQRDSG--VFDCWAR---FAEPWGLRAYVTWIALMVFVAPALGIAACQVLIFREIHASGRRPSEGAHVSAAMAKTVRMTLVIVIVYVLCWAPFFLVQLWAAWDPEA-LERPPFVLLMLLASLNSCTNPWIYASFSSSVSSELRSLLCCAQRHTTHSLGPQDSSLMKDTPS---
>EDG2_MOUSE
QPQFTAMNEQQCFYNESIAFFYNRSGKYLATEWNTVSKLVMGLGITVCVFIMLANLLVMVAIYVNRRFHFPIYYLMANLAAADFFAG-LAYFYLMFNTGPNTRRLTVSTWLLRQGLIDTSLTASVANLLAIAIERHITVFRMQLHTRMSNRRVVVVIVVIWTMAIVMGAIPSVGWNCI-C--DIDHCSNM------APLYSDSYLVFWAIFNLVTFVVMVVLYAHIFGYVRQRRMSSSGPRRNRDTMMSLLKTVVIVLGAFIVCWTPGLVLLLLDVCCPQC-DVLAYEKFFLLLAEFNSAMNPIIYSYRDKEMSATFRQILCCQRNENPNGPTEGSNHTILAGVHSND
>EDG1_RAT
VKALRSQVSDYGNYDIIVRHYNYTGKLNIGVEKDHGIKLTSVVFILICCLIILENIFVLLTIWKTKKFHRPMYYFIGNLALSDLLAG-VAYTANLLLSGATTYKLTPAQWFLREGSMFVALSASVFSLLAIAIERYITMLKMKLHNGSNSSRSFLLISACWVISLILGGLPIMGWNCI-S--SLSSCSTV------LPLYHKHYILFCTTVFTLLLLSIVILYCRIYSLVRTRLTFISKASRSSEKSLALLKTVIIVLSVFIACWAPLFILLLLDVGCKAKCDILYKAEYFLVLAVLNSGTNPIIYTLTNKEMRRAFIRIISCCKCPNGDSAGKFKEFSRSKSDNSSH
>CKR8_HUMAN
DYTLDLSVTTVTDYYYPDIFSSPC---DAELIQTNGKLLLAVFYCLLFVFSLLGNSLVILVLVVCKKLRSITDVYLLNLALSDLLFVFSFPFQTYYLLDQ--WVFGTVMCKVVSGFYYIGFYSSMFFITLMSVDRYLAVVHAVALKVRTIRMGTTLCLAVWLTAIMATIPLLVFYQVAS-ED-VLQCYSF-YNQQTLKWKIFTNFKMNILGLLIPFTIFMFCYIKILHQL---------KRCQNHNKTKAIRLVLIVVIASLLFWVPFNVVLFLTSLHSMHQQLTYATHVTEIISFTHCCVNPVIYAFVGEKFKKHLSEIFQKS-CSQIFNYLGRQKS------SSCQ
>B1AR_MACMU
LPDGVATAARLLVPASPPASLLPPASEGPEPLSQQWTAGMGLLMALIVLLIVAGNVLVIVAIAKTPRLQTLTNLFIMSLASADLVMGLLVVPFGATIVVWGRWEYGSFFCELWTSVDVLCVTASIETLCVIALDRYLAITSPFYQSLLTRARARGLVCTVWAISALVSFLPILMHWWRAR-R--KCCDFV--------TNRAYAIASSVVSFYVPLCIMAFVYLRVFREAQKQNGRRRPSRLVALREQKALKTLGIIMGVFTLCWLPFFLANVVKAFHREL-VPDRLFVFFNWLGYANSAFNPIIYCRSP-DFRNAFQRLLCCARRAARRRHAAHGCLARPGPPPSPG
>DBDR_.ENLA
FQHLDSDQVASWQSPEMLMNKSVSRESQRRKELVAGQIVTGSLLLLLIFWTLFGNILVCTAVMRFRHRSRVTNIFIVSLAVSDLLVALLVMPWKAVAEVAGHWPFG-AFCDIWVAFDIMCSTASILNLCVISVDRYWAISSPFYERKMTQRVALLMISTAWALSVLISFIPVQLSWHKSDH-STGNCDSS--------LNRTYAISSSLISFYIPVAIMIVTYTRIYRIAQIQSCSQTSLRTSIKKETKVLKTLSIIMGVFVCCWLPFFILNCMVPFCDRSCVSETTFDIFVWFGWANSSLNPIIYAFNA-DFRKVFSSLLGCGHWCSTTPVETVNYNQDTLFHK---
>O00155
PTEPWSPSPGSAPWDYSGLDGLEELELCPAGDLPYGYVYIPALYLAAFAVGLLGNAFVVWLLAGRRGPRRLVDTFVLHLAAADLGFVLTLPLWAAAAARRP-WPFGDGLCKLSTFALAGTRSAGALLLAGMSVDRYLAVVKLLARPLRTPRCAVASCCGVWAVALLAGLPSLVYRGLQPLPG-DSQCGE-----EPSHAFQGLSLLLLLLTFVLPLVVTLFCYCRISRRL-------RRPPHVGRARRNSLRIIFAIESTFVGSWLPFSALRAVFHLARLGLALRWGLTIATCLAFVNSCANPLIYLLLDRSFRARALDGACGRTGRLARRISSASSVFRCRAQAANT
>NY4R_HUMAN
HLLALLLPKSPQGENRSKPLGTPYNFSEHCQDSVDVMVFIVTSYSIETVVGVLGNLCLMCVTVRQKEKANVTNLLIANLAFSDFLMCLLCQPLTAVYTIMDYWIFGETLCKMSAFIQCMSVTVSILSLVLVALERHQLIINPT-GWKPSISQAYLGIVLIWVIACVLSLPFLANSILENAD--KVVCTES---WPLAHHRTIYTTFLLLFQYCLPLGFILVCYARIYRRLQRQGRVKGTYSLRAGHMKQVNVVLVVMVVAFAVLWLPLHVFNSLEDWHHEACHGNLIFLVCHLLAMASTCVNPFIYGFLNTNFKKEIKALVLTCQQSAPLEES---TVHTEVSKGSLR
>SSR2_MOUSE
QLNGSQVWVSSPFDLNGSLGPSNGSNQTEPYYDMTSNAVLTFIYFVVCVVGLCGNTLVIYVILRYAKMKTITNIYILNLAIADELFMLGLPFLAMQVALVH-WPFGKAICRVVMTVDGINQFTSIFCLTVMSIDRYLAVVHPISAKWRRPRTAKMINVAVWCVSLLVILPIMIYAGLRSWG--RSSCTIN-WPGESGAWYTGFIIYAFILGFLVPLTIICLCYLFIIIKVKSSGIR-VGSSKRKKSEKKVTRMVSIVVAVFIFCWLPFYIFNVSSVSVAISPALKGMFDFVVILTYANSCANPILYAFLSDNFKKSFQNVLCLVKVSGTEDGERSDLNETTETQRTLL
>EDG2_SHEEP
QPQFTAMNEPQCFYNESIAFFYNRSGKYLATEWNTVSKLVMGLGITVCIFIMLANLLVMVAIYVNRRFHFPIYYLMANLAAADFFAG-LAYFYLMFNTGPNTRRLTVSTWLLRQGLIDTTVTASVANLLAIAIERHITVFRMQLHTRMSNRRVVVVIVVIWTMAIVMGAIPSVGWNCI-C--DIENCSNM------APLYSDSYLVFWAIFNLVTFVVMVVLYAHIFGYVRQRRMSSSGPRRNRDTMMSLLKTVVIVLGAFIICWTPGLVLLLLDVCCPQC-DVLAYEKFFLLLAEFNSAMNPIIYSYRDKEMSATFRQILCCQRSENTSGPTEGSNHTILAGVHSND
>5H6_MOUSE
-----------MVPEPGPVNSSTPAWGPGPPPAPGGSGWVAAALCVVIVLTAAANSLLIALICTQPALRNTSNFFLVSLFTSDLMVGLVVMPPAMLNALYGRWVLARGLCLLWTAFDVMCCSASILNLCLISLDRYLLILSPLYKLRMTAPRALALILGAWSLAALASFLPLLLGWH---E-APGQCRLL--------ASLPYVLVASGVTFFLPSGAICFTYCRILLAARKQMESRRLTTKHSRKALKASLTLGILLSMFFVTWLPFFVASIAQAVCDC--ISPGLFDVLTWLGYCNSTMNPIIYPLFMRDFKRALGRFVPCVHCPPEHRASPASSGARPGLSLQQV
>ACTR_CAVPO
-------------MKHIIHASGNVNGTARNNSDCPHVALPEEIFFIISITGVLENLIIILAVIKNKNLQFPMYFFICSLAISDMLGSLYKILESILIMFRNMGSFETTTDDIIDTMFILSLLGSIFSLLAIAVDRYITIFHALYHSIVTMHRTIAVLSIIWTFCIGSGITMVLFFS----------------------HHHVPTVLTFTSLFPLMLVFILCLYVHMFLMARSH---ARNISTLPRGNMRGAITLTILLGVFIFCWAPFILHILLVTFCPNNTCYISLFHVNGMLIMCNAVIDPFIYAFRSPELRSAFRRMISYSKCL---------------------
>OPSD_ALLMI
MNGTEGPDFYIPFSNKTGVVRSPFEYPQYYLAEPWKYSALAAYMFMLIILGFPINFLTLYVTVQHKKLRSPLNYILLNLAVADLFMVLGGFTTTLYTSMNGYFVFGVTGCYFEGFFATLGGEVALWCLVVLAIERYIVVCKPMSNFRFGENHAIMGVVFTWIMALTCAAPPLVGWSRYIPEGMQCSCGVDYYTLKPEVNNESFVIYMFVVHFAIPLAVIFFCYGRLVCTVKEAAAQQQESATTQKAEKEVTRMVIIMVVSFLICWVPYASVAFYIFSNQGSDFGPVFMTIPAFFAKSSAIYNPVIYIVMNKQFRNCMITTLCCGKNPLGDDETATGTSSVSTSQVSPA
>OPSD_RANTE
MNGTEGPNFYIPMSNKTGVVRSPFEYPQYYLAEPWKYSILAAYMFLLILLGFPINFMTLYVTIQHKKLRTPLNYILLNLAFANHFMVLCGFTITLYTSLHGYFVFGQSGCYFEGFFATLGGEIALWSLVALAIERYIVVCKPMSNFRFGENHAMMGVAFTWIMALACAVPPLFGWSRYIPEGMQCSCGVDYYTLKPEINNESFVIYMFVVHFLIPLIIITFCYGRLVCTVKEAAAQQQESATTQKAEKEVTRMVIIMVIFFLICWVPYAYVAFYIFCNQGSEFGPIFMTVPAFFAKSSAIYNPVIYIMLNKQFRNCMITTLCCGKNPFGDDDASSAATSVSTSQVSPA
>A2AC_CAVPO
GPNASGAGEG----GGGVNASGAVWGPPPSQYSAGAVAGLAAVVGFLIVFTVVGNVLVVIAVLTSRALRAPQNLFLVSLASADILVATLVMPFSLANELMAYWYFGQVWCGVYLALDVLFCTSSIVHLCAISLDRYWSVTQAVYNLKRTPRRVKATIVAVWLISAIISFPPLVSFYR-------PRCGLN--------DETWYILSSCIGSFFAPCLIMGLVYARIYRVAKLRRRAVCRRKVAQAREKRFTFVLAVVMGVFVLCWFPFFFSYSLYGICREAQLPTPLFKFFFWIGYCNSSLNPVIYTIFNQDFRRSFKHILFRRRRRGFRQ-----------------
>PE23_RAT
---------MAGVWAPEHSVEAHSNQS---SAADGCGSVSVAFPITMMVTGFVGNALAMLLVVRSYRRKKSFLLCIGWLALTDLVGQLLTSPVVILVYLSQRLDPSGRLCTFFGLTMTVFGLSSLLVASAMAVERALAIRAPHYASHMKTRAT-PVLLGVWLSVLAFALLPVLGVG----RYSGTWCFISNETDSAREPGSVAFASAFACLGLLALVVTFACNLATIKALVSRAKASQSSAQWGRITTETAIQLMGIMCVLSVCWSPLLIMMLKMIFNQMSEKECNSFLIAVRLASLNQILDPWVYLLLRKILLRKFCQIRDHTN-YASSSTSLPCWSDQLER-----
>OPSF_ANGAN
MNGTEGPNFYVPMSNVTGVVRSPFEYPQYYLAEPWAYSALAAYMFFLIIAGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVFGPTGCNIEGFFATLGGEIALWCLVVLAVERWMVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLFGWSRYIPEGMQCSCGMDHYAPNPETYNESFVIYMFICHFTIPLTVISFCYGRLVCTVKEATAQQQESETTQRAEREVTRMVIIMVISFLVCWVPYASVAWYIFTHQGSSFGPIFMTIPAFFAKSSSLYNPLIYICMNKQSRNCMITTLCCGKNPFEEEEGASTASSVSS--VSPA
>OPSR_ORYLA
NE--DTTRGSAFTYTNSNHTRDPFEGPNYHIAPRWVYNLATLWMFFVVVLSVFTNGLVLVATAKFKKLRHPLNWILSNLAIADLGETVFASTISVCNQFFGYFILGHPMCVFEGYVVSTCGIAALWSLTIISWERWVVVCKPFGNVKFDAKWAIGGIVFSWVWSAVWCAPPVFGWSRYWPHGLKTSCGPDVFSGSDDPGVQSYMIVLMITCCIIPLAIIILCYLAVWLAIRAVAMQQKESESTQKAEREVSRMVVVMIVAYCVCWGPYTFFACFAAANPGYAFHPLAAAMPAYFAKSATIYNPVIYVFMNRQFRTCIMQLFGKQVDDG-----SEVSS------VAPA
>CKR5_CERAE
----MDYQVSSPTYDIDNYTSEPC---QKINVKQIAARLLPPLYSLVFIFGFVGNILVVLILINCKRLKSMTDIYLLNLAISDLLFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWVVAVFASLPRIIFTRSQR-EG-HYTCSSHFPYSQYQFWKNFQTLKIVILGLVLPLLVMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCSIFA-------SSVY
>D2DR_CERAE
---MDPLNLSWYDDDLERQNWSRPFNGSDGKADRPHYNYYATLLTLLIAVIVFGNVLVCMAVSREKALQTTTNYLIVSLAVADLLVATLVMPWVVYLEVVGEWKFSKIHCDIFVTLDVMMCTASILNLCAISIDRYTAVAMPMNTRYSSKRRVTVMIAIVWVLSFTISCPLLFGLN---Q----NECII---------ANPAFVVYSSIVSFYVPFIVTLLVYIKIYIVLRRRRTSMSRRKLSQQKEKKATQMLAIVLGVFIICWLPFFITHILNIHCDCN-IPPVLYSAFTWLGYVNSAVNPIIYTTFNIEFRKAFLKILHC-------------------------
>AA3R_SHEEP
-------------------------MPVNSTAVSWTSVTYITVEILIGLCAIVGNVLVIWVVKLNPSLQTTTFYFIVSLALADIAVGVLVMPLAIVISLG--VTIHFYSCLFMTCLMLIFTHASIMSLLAIAVDRYLRVKLTVYRRVTTQRRIWLALGLCWLVSFLVGLTPMFGWN-MKA-D-FLPCRFR-----SVMRMDYMVYFSFFLWILVPLVVMCAIYFDIFYIIRNR--SSRETGAFYGREFKTAKSLLLVLFLFALCWLPLSIINCILYFDGQ--VPQTVLYLGILLSHANSMMNPIVYAYKIKKFKETYLLILKACVMCQPSKSMDPS------------
>IL8B_RAT
SGDIDSYN-YSSDPPFTLSDAAPC----PSANLDINRYAVVVIYVLVTLLSLVGNSLVMLVILYNRSTCSVTDVYLLNLAIADLFFALTLPVWAASKVNG--WIFGSFLCKVFSFLQEITFYSSVLLLACISMDRYLAIVHATSTLIQKRHLVKFVCITMWFLSLVLSLPIFILRTTVK-PS-TVVCYEN-IGNNTSKWRVVLRILPQTYGFLLPLLIMLFCYGFTLRTL---------FKAHMGQKHRAMRVIFAVVLVFLLCWLPYNIVLFTDTLMRTKNEINKALEATEILGFLHSCLNPIIYAFIGQKFRHGLLKIMANYGLVSKEFLAKEG------------
>NY5R_MOUSE
ASPAWEDYRGTENNTSAARNTAFPVWEDYRGSVDDLQYFLIGLYTFVSLLGFMGNLLILMAVMKKRNQKTTVNFLIGNLAFSDILVVLFCSPFTLTSVLLDQWMFGKAMCHIMPFLQCVSVLVSTLILISIAIVRYHMIKHPI-SNNLTANHGYFLIATVWTLGFAICSPFPVFHSLVESS--KYLCVES---WPSDSYRIAFTISLLLVQYILPLVCLTVSHTSVCRSISCGAQEKRSLTRIKKRSRSVFYRLTILILVFAVSWMPLHVFHVVTDFNDNLRHFKLVYCICHLLGMMSCCLNPILYGFLNNGIKADLRALIHCLHMS---------------------
>OPS3_DROPS
ARLSAESRLLGWNVPPDELRHIPEHWLIYPEPPESMNYLLGTLYIFFTVISMIGNGLVMWVFSAAKSLRTPSNILVINLAFCDFMMMIKTPIFIYNSFHQG-YALGHLGCQIFGVIGSYTGIAAGATNAFIAYDRYNVITRPM-EGKMTHGKAIAMIIFIYLYATPWVVACYTESWGRFPEGYLTSCTFD--YLTDNFDTRLFVACIFFFSFVCPTTMITYYYSQIVGHVFSHNVDSNVDKSKEAAEIRIAKAAITICFLFFASWTPYGVMSLIGAFGDKTLLTPGATMIPACTCKMVACIDPFVYAISHPRYRMELQKRCPWLAISEKAPESRAAEQQQTTAA----
>B1AR_MOUSE
LPDGAATAARLLVLASPPASLLPPASEGSAPLSQQWTAGMGLLVALIVLLIVVGNVLVIVAIAKTPRLQTLTNLFIMSLASADLVMGLLVVPFGATIVVWGRWEYGSFFCELWTSVDVLCVTASIETLCVIALDRYLAITSPFYQSLLTRARARALVCTVWAISALVSFLPILMHWWRAR-R--KCCDFV--------TNRAYAIASSVVSFYVPLCIMAFVYLRVFREAQKQNGRRRPSRLVALREQKALKTLGIIMGVFTLCWLPFFLANVVKAFHRDL-VPDRLFVFFNWLGYANSAFNPIIYCRSP-DFRKAFQRLLCCARRAACRRRAAHGCLARAGPPPSPG
>O6A1_HUMAN
------------MEWRNHSGRVSEFVLLGFPAPVPIQVILFALLLLAYVLVLTENTLIIMAIRNHSTLHKPMYFFLANMSFLEIWYVTVTIPKMLAGFVGSKQLISFEGCMTQLYFFLGLGCTECVLLAVMANDRYMAICYLLNPVIVSGRLCVQMAAGSWAGGFGISMVKVFLISGL---SNHFFCDVSNLSCTDMSTAELTDFILAIFILLGPLSVTGASYVAITGAV--------MHIPSAAGRYKAFSTCASHFNVVIIFYAASIFIYARPKALS---AFDTNKLVSVLYAVIVPLLNPIIYCLRNQEVKRALCCILHLYQHQDPDPKKGSR------------
>GPR8_HUMAN
PLDSRGSFSLPTMGANVSQDNGTGHNATFSEPLPFLYVLLPAVYSGICAVGLTGNTAVILVILRAPKMKTVTNVFILNLAVADGLFTLVLPVNIAEHLLQY-WPFGELLCKLVLAVDHYNIFSSIYFLAVMSVDRYLVVLATVHMPWRTYRGAKVASLCVWLGVTVLVLPFFSFAGVYS-LQ-VPSCGLS-FPWPERVWFKASRVYTLVLGFVLPVCTICVLYTDLLRRLRAVR--RSGAKALGKARRKVTVLVLVVLAVCLLCWTPFHLASVVALTTDLPPLVISMSYVITSLTYANSCLNPFLYAFLDDNFRKNFRSILRC-------------------------
>CCKR_.ENLA
SSTNGTHNLTTANWPPWNLNCTPILDRKKPSPSDLNLWVRIVMYSVIFLLSVFGNTLIIIVLVMNKRLRTITNSFLLSLALSDLMVAVLCMPFTLIPNLMENFIFGEVICRAAAYFMGLSVSVSTFNLVAISIERYSAICNPLSRVWQTRSHAYRVIAATWVLSSIIMIPYLVYK----DRRVGHQCRLV---WPSKQVQQAWYVLLLTILFFIPGVVMIVAYGLISRELYRGKMDINNSEAKLMAKKRVIRMLIVIVAMFFICWMPIFVANTWKAFDELSTLTGAPISFIHLLSYTSACVNPLIYCFMNKRFRKAFLGTFSSCIKP----CRNFRATGASLSKFSYT
>ETBR_RAT
SSAPAEVTKGGRVAGVPPRS-FPPPCQRKIEINKTFKYINTIVSCLVFVLGIIGNSTLLRIIYKNKCMRNGPNILIASLALGDLLHIIIDIPINAYKLLAGDWPFGAEMCKLVPFIQKASVGITVLSLCALSIDRYRAVASWSIKGIGVPKWTAVEIVLIWVVSVVLAVPEAIGFDVITRV--LRVCMLNQKTAFMQFYKTAKDWWLFSFYFCLPLAITAIFYTLMTCEMLRKSGM-IALNDHLKQRREVAKTVFCLVLVFALCWLPLHLSRILKLTLYDQSFLLVLDYIGINMASLNSCINPIALYLVSKRFKNCFKSCLCCWCQTFE-EKQSLEFKANDHGYDNFR
>EDG2_HUMAN
QPQFTAMNEPQCFYNESIAFFYNRSGKHLATEWNTVSKLVMGLGITVCIFIMLANLLVMVAIYVNRRFHFPIYYLMANLAAADFFAG-LAYFYLMFNTGPNTRRLTVSTWLLRQGLIDTSLTASVANLLAIAIERHITVFRMQLHTRMSNRRVVVVIVVIWTMAIVMGAIPSVGWNCI-C--DIENCSNM------APLYSDSYLVFWAIFNLVTFVVMVVLYAHIFGYVRQRRMSSSGPRRNRDTMMSLLKTVVIVLGAFIICWTPGLVLLLLDVCCPQC-DVLAYEKFFLLLAEFNSAMNPIIYSYRDKEMSATFRQILCCQRSENPTGPTESSNHTILAGVHSND
>5H1D_RAT
SLPNQSLEGLPQEASNRSLNAT---GAWDPEVLQALRISLVVVLSIITLATVLSNAFVLTTILLTKKLHTPANYLIGSLATTDLLVSILVMPISIAYTTTRTWNFGQILCDIWVSSDITCCTASILHLCVIALDRYWAITDALYSKRRTAGHAAAMIAAVWAISICISIPPLFWR----H-E--SDCLVN-------TSQISYTIYSTCGAFYIPSILLIILYGRIYVAARSRKLALERKRISAARERKATKTLGIILGAFIICWLPFFVVSLVLPICRDSWIHPALFDFFTWLGYLNSLINPVIYTVFNEDFRQAFQRVVHFRKAS---------------------
>ET1R_BOVIN
ELSFVVTTHQPTNLALPSNGSMHNYCPQQTKITSAFKYINTVISCTIFIVGMVGNATLLRIIYQNKCMRNGPNALIASLALGDLIYVVIDLPINVFKLLAGRNDFGVFLCKLFPFLQKSSVGITVLNLCALSVDRYRAVASWSVQGIGIPLVTAIEIVSIWILSFILAIPEAIGFVMVPRT--HRTCMLNATSKFMEFYQDVKDWWLFGFYFCMPLVCTAIFYTLMTCEMLNRGSL-IALSEHLKQRREVAKTVFCLVVIFALCWFPLHLSRILKKTVYDESFLLLMDYIGINLATMNSCINPIALYFVSKKFKNCFQSCLCCCCYQSKSLMTSVPWKNHEQNNHNTE
>GRE2_BALAM
----MSGGEASITGRTAPELN-ASAAPLDDERELGETVAATALLLAIILVTIVGNSLVIISVFTYRPLRSVQNFFVVSLAVADLTVALFVLPLNVAYRLLNQWLLGSYLCQMWLTCDILCCTSSILNLCVIALDRYWAITDPIYAQKRTIRRVNTMIAAVWALSLVISVPPLLGWNDW-T-E--TPCTLT--------QR-LFVVYSSSGSFFIPLIIMSVVYAKIFFATKRRSVHEEKQRISLSKERKAARVLGVIMGVFVVCWLPFFLMYAIVPFCTNCPPSQRVVDFVTWLGYVNSSLNPIIYTIYNKDFRTAFSRLLRCDRRMSA-------------------
>LSHR_MOUSE
ENELSGWDYDYDFCSPKTLQCTPEPDAFNPCEDIMGYAFLRVLIWLINILAIFGNLTVLFVLLTSRYKLTVPRFLMCNLSFADFCMGLYLLLIASVDSQTKGWQTG-SGCSAAGFFTVFASELSVYTLTVITLERWHTITYAVLDQKLRLRHAIPIMLGGWIFSTLMATLPLVGVSSY----KVSICLPM-----VESTLSQVYILSILLLNAVAFVVICACYVRIYFAVQNP------ELTAPNKDTKIAKKMAILIFTDFTCMAPISFFAISAAFKVPLITVTNSKVLLVLFYPVNSCANPFLYAVFTKAFQRDFFLLLSRFGCCKHRAELYRRFNSKNGFPRSSK
>MSHR_MOUSE
STQEPQKSLLGSLNSN--ATSHLGLATNQSEPWCLYVSIPDGLFLSLGLVSLVENVLVVIAITKNRNLHSPMYYFICCLALSDLMVSVSIVLETTIILLLEVVALVQQLDNLIDVLICGSMVSSLCFLGIIAIDRYISIFYALYHSIVTLPRARRAVVGIWMVSIVSSTLFITYYY----------------------KKHTAVLLCLVTFFLAMLALMAILYAHMFTRACQHIAQKRRRSIRQGFCLKGAATLTILLGIFFLCWGPFFLHLLLIVLCPQHSCIFKNFNLFLLLIVLSSTVDPLIYAFRSQELRMTLKEVLLCSW-----------------------
>5H5A_RAT
LPINLTSFSLSTPSTLEPNRSDTEALRTSQSFLSAFRVLVLTLLGFLAAATFTWNLLVLATILRVRTFHRVPHNLVASMAISDVLVAVLVMPLSLVHELSGRWQLGRRLCQLWIACDVLCCTASIWNVTAIALDRYWSITRHLYTLRARKRVSNVMILLTWALSAVISLAPLLFGWGE-S-E-SEECQVS--------REPSYTVFSTVGAFYLPLCVVLFVYWKIYKAAKFRATVTEGDTWREQKEQRAALMVGILIGVFVLCWFPFFVTELISPLCSW-DIPALWKSIFLWLGYSNSFFNPLIYTAFNRSYSSAFKVFFSKQQ-----------------------
>BRS3_HUMAN
QTLISITNDTESSSSVVSNDNTNKGWSGDNSPGIEALCAIYITYAVIISVGILGNAILIKVFFKTKSMQTVPNIFITSLAFGDLLLLLTCVPVDATHYLAEGWLFGRIGCKVLSFIRLTSVGVSVFTLTILSADRYKAVVKPLRQPSNAILKTCVKAGCVWIVSMIFALPEAIFSNVYTMT--FESCTSY---VSKKLLQEIHSLLCFLVFYIIPLSIISVYYSLIARTLYKSIPTQSHARKQIESRKRIARTVLVLVALFALCWLPNHLLYLYHSFTSQTAMHFIFTIFSRVLAFSNSCVNPFALYWLSKSFQKHFKAQLFCCKAERPEPPVADTMGTVPGTGSIQM
>OPSD_CATBO
GGGF-GNQTVVDKVPPEMLHLVDAHWYQFPPMNPLWHAILGFVIGILGMISVIGNGMVIYIFTTTKSLRTPSNLLVINLAISDFLMMLSMSPAMVINCYYETWVLGPLVCELYGLTGSLFGCGSIWTMTMIAFDRYNVIVKGLSAKPMTINGALLRILGIWFFSLGWTIAPMFGWNRYVPEGNMTACGTD--YLTKDLLSRSYILVYSFFCYFLPLFLIIYSYFFIIQAVAAHMNVRSAENQSTSAECKLAKVALMTISLWFMAWTPYLVINYAGIFETVK-INPLFTIWGSLFAKANAVYNPIVYGISHPKYRAALFQRFPSLACSSGPAG-ADTTEGTEKPAA---
>GASR_PRANA
GSSLCHPGVSLLNSSSAGNLSCEPPRIRGTGTRELELAIRITLYAVIFLMSIGGNMLIIVVLGLSRRLRTVTNAFLLSLAVSDLLLAVACMPFTLLPNLMGTFIFGTVICKAVSYLMGVSVSVSTLNLVAIALERYSAICRPLARVWQTRSHAARVILATWLLSGLLMVPYPVYTVVQP---V-LQCMHR---WPSARVRQTWSVLLLMLLFFIPGVVMAVAYGLISRELYLGTPGASANQAKLLAKKRVVRMLLVIVLLFFLCWLPIYSANTWCAFDGPGALSGAPISFIHLLSYASACVNPLVYCFMHRRFRQACLDTCARCCPRPPRARPRPLPSIASLSRLSYT
>P2Y3_CHICK
----------------MSMANFTGGRNSCTFHEEFKQVLLPLVYSVVFLLGLPLNAVVIGQIWLARKALTRTTIYMLNLAMADLLYVCSLPLLIYNYTQKDYWPFGDFTCKFVRFQFYTNLHGSILFLTCISVQRYMGICHPLWHKKKGKKLTWLVCAAVWFIVIAQCLPTFVFASTG-----RTVCYDL-SPPDRSTSYFPYGITLTITGFLLPFAAILACYCSMARILCQK---ELIGLAVHKKKDKAVRMIIIVVIVFSISFFPFHLTKTIYLIVRSSQAFAIAYKCTRPFASMNSVLDPILFYFTQRKFRESTRYLLDKMSSKWRQDHCISY------------
>AG2R_.ENLA
-----MSNASTVETSDVERIAVNC---SKSGMHNYIFIAIPIIYSTIFVVGVFGNSMVVIVIYSYMKMKTVASIFLMNLALSDLCFVITLPLWAAYTAMHYHWPFGNFLCKVASTAITLNLYTTVFLLTCLSIDRYSAIVHPMSRIWRTAMVARLTCVGIWLVAFLASMPSIIYRQIYL-TN--TVCAIV-YDSGHIYFMVGMSLAKNIVGFLIPFLIILTSYTLIGKTLKEV------YRAQRARNDDIFKMIVAVVLLFFFCWIPYQVFTFLDVLIQMDDIVDTGMPITICIAYFNSCLNPFLYGFFGKNFRKHFLQLIKYIPPKMRTHASVNTSLSD-------T
>FSHR_PIG
FDTMYSEFDYDLCNEVVDVICSPEPDTFNPCEDIMGHDILRVLIWFISILAITGNIIVLVILITSQYKLTVPRFLMCNLAFADLCIGIYLLLIASVDIHTKTWQTG-AGCDAAGFFTVFASELSVYTLTAITLERWHTITHAMLQCKVQLRHAASIMLVGWIFAFTVALFPIFGISSY----KVSICLPM-----IDSPLSQLYVVSLLVLNVLAFVVICGCYTHIYLTVRNP------NIMSSSSDTKIAKRMAMLIFTDFLCMAPISFFAISASLKVPLITVSKSKILLVLFYPINSCANPFLYAIFTKNFRRDVFILLSKFGCYEMQAQTYRTNIHPRNGHCPPA
>MSHR_DAMDA
PVLGSQRRLLGSLNCTPPATFPLTLAPNRTGPQCLEVSIPDGLFLSLGLVSLVENVLVVAAIAKNRNLHSPMYYFICCLAMSDLLVSVSNVLETAVMLLLEAAAVVQQLDNVIDVLICGSMVSSLCFLGAIAVDRYISIFYALYHSVVTLPRAWRIIAAIWVASILTSLLFITYYY----------------------NNHTVVLLCLVGFFIAMLALMAVLYVHMLARACQHIARKRQRPIHQGFGLKGAATLTILLGVFFLCWGPFFLHLSLIVLCPQHGCIFKNFNLFLALIICNAIVDPLIYAFRSQELRKTLQEVLQCSW-----------------------
>ETBR_HORSE
LSAPPQMPKAGRTAGAQRRTLPPPPCERTIEIKETFKYINTVVSCLVFVLGIIGNSTLLRIIYKNKCMRNGPNILIASLALGDLLHIIIDIPINVYKLLAEDWPFGVEMCKLVPFIQKASVGITVLSLCALSIDRYRAVASWSIKGIGVPKWTAVEIVLIWVVSVVLAVPEAVGFDMITRI--LRICLLHQKTAFMQFYKNAKDWWLFSFYFCLPLAITAFFYTLETCEMLRKSGM-IALNDHLKQRREVAKTVFCLVLVFALCWLPLHLSRILKHTLYDQSFLLVLEYIGINMASLNSCINPIALYLVSKRFKNCFKWCLCCWCQSFE-EKQSLEFKANDHGYDNFR
>MSHR_CAPHI
PALGSPRRLLGSLNCTPPATLPLTLAPNRTGPQCLEVSIPDGLFLSLGLVSLVENVLVVAAIAKNRNLHSPMYYFICCLAMSDLLVSVSNVLETAVMLLLEAAAVVQQLDNVIDVLICSSMVSSLCFLGAIAVDRYISIFYALYHSVVTLPRAWRIIAAIWVASILTSVLSITYYY----------------------NNHTVVLLCLVGFFIAMLALMAVLYVHMLARACQHIARKRQRPIHQGFGLKGAATLTILLGVFFLCWGPFFLHLSLIVLCPQHGCIFKNFNLFLALIICNAIVDPLIYAFRSQELRKTLQEVLQCSW-----------------------
>AG2S_RAT
------MTLNSSTEDGIKRIQDDC---PKAGRHNYIFVMIPTLYSIIFVVGIFGNSLVVIVIYFYMKLKTVASVFLLNLALADLCFLLTLPLWAVYTAMEYRWPFGNHLCKIASASVSFNLYASVFLLTCLSIDRYLAIVHPMSRLRRTMLVAKVTCIIIWLMAGLASLPAVIYRNVYF-TN--TVCAFH-YESQNSTLPIGLGLTKNILGFVFPFLIILTSYTLIWKALKKA----YKIQKNTPRNDDIFRIIMAIVLFFFFSWVPHQIFTFLDVLIQLGDIVDTAMPITICIAYFNNCLNPLFYGFLGKKFKKYFLQLLKYIPPTAKSHAGLSTRPSD-------N
>D5DR_FUGRU
NFYNETEPTEPRGGVDPLRVVTAAEDVPAPVGGVSVRALTGCVLCALIVSTLLGNTLVCAAVIKFRHRSKVTNAFVVSLAVSDLFVAVLVMPWRAVSEVAGVWLFG-RFCDTWVAFDIMCSTASILNLCVISMDRYWAISNPFYERRMTRRFAFLMIAVAWTLSVLISFIPVQLNWHRASS-EQGDCNAS--------LNRTYAISSSLISFYIPVLIMVGTYTRIFRIAQTQRASESALKTSFKRETKVLKTLSVIMGVFVFCWLPFFVLNCVVPFCDVDCVSDTTFNIFVWFGWANSSLNPVIYAFNA-DFRKAFTTILGCSKFCSSSAVQAVDYHHDTTLQK---
>PI2R_BOVIN
------------------------MADSCRNLTYVRDSVGPATSTLMFVAGVVGNGLALGILGARRHRPSAFAVLVTGLGVTDLLGTCFLSPAVFAAYARNSARGRPALCDAFAFAMTFFGLASTLILFAMAVERCLALSHPYYAQLDGPRRARLALPAIYAFCTIFCSLPFLGLGQHQQYCPGSWCFIR---MRSAEPGGCAFLLAYASLVALLVAAIVLCNGSVTLSLCRMQRRRCPRPRAGEDEVDHLILLALMTGIMAVCSLPLTPQIRGFTQAIAPDSSEMGDLLAFRFNAFNPILDPWVFILFRKSVFQRLKLWFCCLYSRPAQGDSRTSRKDSSAPPALEG
>A1AA_RABIT
----------MVFLSGNASDSSNCT-HPPAPVNISKAILLGVILGGLILFGVLGNILVILSVACHRHLHSVTHYYIVNLAVADLLLTSTVLPFSAIFEILGYWAFGRVFCNIWAAVDVLCCTASIISLCVISIDRYIGVSYPLYPTIVTQRRGLRALLCVWAFSLVISVGPLFGWRQP--AP--TICQIN--------EEPGYVLFSALGSFYVPLTIILAMYCRVYVVAKREAKNFSVRLLKFSREKKAAKTLGIVVGCFVLCWLPFFLVMPIGSFFPDFKPPETVFKIVFWLGYLNSCINPIIYPCSSQEFKKAFQNVLKIQCLRRKQSSKHALSQ----------
>OPS4_DROPS
SSGSDELQFLGWNVPPDQIQYIPEHWLTQLEPPASMHYMLGVFYIFLFFASTLGNGMVIWIFSTSKSLRTPSNMFVLNLAVFDLIMCLKAPIFIYNSFHRG-FALGNTWCQIFASIGSYSGIGAGMTNAAIGYDRYNVITKPM-NRNMTFTKAVIMNIIIWLYCTPWVVLPLTQFWDRFPEGYLTSCSFD--YLSDNFDTRLFVGTIFLFSFVVPTLMILYYYSQIVGHVFNHNVESNVDKSKETAEIRIAKAAITICFLFFVSWTPYGVMSLIGAFGDKSLLTPGATMIPACTCKLVACIEPFVYAISHPRYRMELQKRCPWLGVNEKSGEASSAQTQQTSAA----
>5H1F_RAT
-------------MDFLNSSD-QNLTSEELLNRMPSKILVSLTLSGLALMTTTINCLVITAIIVTRKLHHPANYLICSLAVTDFLVAVLVMPFSIVYIVRESWIMGQGLCDLWLSVDIICCTCSILHLSAIALDRYRAITDAVYARKRTPRHAGITITTVWVISVFISVPPLFWR----S-R--DQCIIK-------HDHIVSTIYSTFGAFYIPLVLILILYYKIYRAARTLLKHWRRQKISGTRERKAATTLGLILGAFVICWLPFFVKELVVNICEKCKISEEMSNFLAWLGYLNSLINPLIYTIFNEDFKKAFQKLVRCRN-----------------------
>NK2R_CAVPO
----MGACVIVTNTNISSGLESNTTGITAFSMPTWQLALWATAYLALVLVAVTGNATVTWIILAHQRMRTVTNYFIVNLALADLCMAAFNAAFNFVYASHNIWYFGRAFCYFQNLFPITAMFVSIYSMTAIAIDRYMAIVHPF-QPRLSAPSTKAVIGGIWLVALALAFPQCFYST---GA---TKCVVAWPEDSRDKSLLLYHLVVIVLIYLLPLTVMFVAYSIIGLTLWRRRHQHGANLRHLQAKKKFVKTMVLVVVTFAICWLPYHLYFILGSFQEDIKFIQQVYLALFWLAMSSTMYNPIIYCCLNRRFRSGFRLAFRCCPWVTPTEE----HTPSFSLRVNRC
>5H1A_RAT
-MDVFSFGQGNNTTASQEPFGTGGNVTSISDVTFSYQVITSLLLGTLIFCAVLGNACVVAAIALERSLQNVANYLIGSLAVTDLMVSVLVLPMAALYQVLNKWTLGQVTCDLFIALDVLCCTSSILHLCAIALDRYWAITDPIYVNKRTPRRAAALISLTWLIGFLISIPPMLGWRT---DP--DACTIS--------KDHGYTIYSTFGAFYIPLLLMLVLYGRIFRAARFRKNEEAKRKMALARERKTVKTLGIIMGTFILCWLPFFIVALVLPFCESSHMPALLGAIINWLGYSNSLLNPVIYAYFNKDFQNAFKKIIKCKFCRR--------------------
>ACM2_PIG
--------------MNNSTNSSNSGLALTSPYKTFEVVFIVLVAGSLSLVTIIGNILVMVSIKVNRHLQTVNNYFLFSLACADLIIGVFSMNLYTLYTVIGYWPLGPVVCDLWLALDYVVSNASVMNLLIISFDRYFCVTKPLYPVKRTTKMAGMMIAAAWVLSFILWAPAILFWQFIV-VE-DGECYIQ------FFSNAAVTFGTAIAAFYLPVIIMTVLYWHISRASKSRVKMPAKKKPPPSREKKVTRTILAILLAFIITWAPYNVMVLINTFCAPC-IPNTVWTIGYWLCYINSTINPACYALCNATFKKTFKHLLMCHYKNIGATR----------------
>GPRV_HUMAN
-----------------------MPFPNCSAPSTVVATAVGVLLGLECGLGLLGNAVALWTFLFRVRVWKPYAVYLLNLALADLLLAACLPFLAAFYLSLQAWHLGRVGCWALRFLLDLSRSVGMAFLAAVALDRYLRVVHPRKVNLLSPQAALGVSGLVWLLMVALTCPGLLISEAAQ--S--TRCHSF-YSRADGSFSIIWQEALSCLQFVLPFGLIVFCNAGIIRALQKR----LREPEKQPKLQRAQALVTLVVVLFALCFLPCFLARVLMHIFQNLCAVAHTSDVTGSLTYLHSVVNPVVYCFSSPTFRSSYRRVFHTLRGKGQAAEPPDF------------
>YT66_CAEEL
----------------------------------------MGVTFHPGIVGNITNLMVLASRR-------LRAMYLRALAVADLLCMLFVLVFVSTEYLAKNKLYQIYQCHLMLTLINWALGAGVYVVVALSLERYISIVFPMFRTWNSPQRATRAIVIAFLIPAIFYVPYAITRYKGK---VTIYSMDD---IYTTFYWQIYKWTREAILRFLPIIILTVLNIQIMIAFRKRMFQNKRKEQGTQKDDTLMYMLGGTVLMSLVCNIPAAINLLLIDETLKKLDYQIFRAVANILEITNHASQFYVFCACSTDYRTTFLQKFPCFKTDYANRDRLRSVIQKQGSVEHTT
>RGR_MOUSE
---------------------MAATRALPAGLGELEVLAVGTVLLMEALSGISLNGLTIFSFCKTPDLRTPSNLLVLSLALADTGISLNALVAAVSSLLRR-WPHGSEGCQVHGFQGFATALASICGSAAVAWGRYHHYCTRR---QLAWDTAIPLVLFVWMSSAFWASLPLMGWGHYDYEPVGTCCTLD--YSRGDRNFISFLFTMAFFNFLVPLFITHTSYRFME--------------QKFSRSGHLPVNTTLPGRMLLLGWGPYALLYLYAAIADVSFISPKLQMVPALIAKTMPTINAINYALHREMVCRGTWQCLSPQKSKKDRTQA---------------
>A2AA_PIG
----MGSLQPEAGNASWNGTEAPGGGARATPYSLQVTLTLVCLAGLLMLFTVFGNVLVIIAVFTSRALKAPQNLFLVSLASADILVATLVIPFSLANEVMGYWYFGKAWCEIYLALDVLFCTSSIVHLCAISLDRYWSITQAIYNLKRTPRRIKAIIVTVWVISAVISFPPLISIEKKAQ-P--PRCEIN--------DQKWYVISSCIGSFFAPCLIMILVYVRIYQIAKRRRGGASRWRGRQNREKRFTFVLAVVIGVFVVCWFPFFFTYTLTAVGCS--VPPTLFKFFFWFGYCNSSLNPVIYTIFNHDFRRAFKKILCRGDRKRIV------------------
>PI2R_HUMAN
------------------------MADSCRNLTYVRGSVGPATSTLMFVAGVVGNGLALGILSARRPRPSAFAVLVTGLAATDLLGTSFLSPAVFVAYARNSARGGPALCDAFAFAMTFFGLASMLILFAMAVERCLALSHPYYAQLDGPRCARLALPAIYAFCVLFCALPLLGLGQHQQYCPGSWCFLR---MRWAQPGGAAFSLAYAGLVALLVAAIFLCNGSVTLSLCRMKRHLGPRPRTGEDEVDHLILLALMTVVMAVCSLPLTIRCFTQAVAPDS-SSEMGDLLAFRFYAFNPILDPWVFILFRKAVFQRLKLWVCCLCLGPAHGDSQTPRRDPRAPSAPVG
>AA3R_CANFA
-------------------------MAVNGTALLLANVTYITVEILIGLCAIVGNVLVIWVVKLNPSLQTTTFYFIVSLALADIAVGVLVMPLAIVISLG--ITIQFYNCLFMTCLLLIFTHASIMSLLAIAVDRYLRVKLTVYRRVTTQRRIWLALGLCWLVSFLVGLTPMFGWN-MKE-H-FLSCQFS-----SVMRMDYMVYFSFFTWILIPLVVMCAIYLDIFYVIRNK--NSKETGAFYGREFKTAKSLFLVLFLFAFSWLPLSIINCITYFHGE--VPQIILYLGILLSHANSMMNPIVYAYKIKKFKETYLLIFKTYMICQSSDSLDSS------------
>CB1R_MOUSE
EFYNKSLSSFKENEDNIQCGENFMDMECFMILNPSQQLAIAVLSLTLGTFTVLENLLVLCVILHSRSRCRPSYHFIGSLAVADLLGSVIFVYSFVDFHVFHR-KDSPNVFLFKLGGVTASFTASVGSLFLTAIDRYISIHRPLYKRIVTRPKAVVAFCLMWTIAIVIAVLPLLGWNCK-K--LQSVCSDI------FPLIDETYLMFWIGVTSVLLLFIVYAYMYILWKAHSHSEDQVTRPDQARMDIRLAKTLVLILVVLIICWGPLLAIMVYDVFGKMNKLIKTVFAFCSMLCLLNSTVNPIIYALRSKDLRHAFRSMFPSCEGTAQPLDNSMGHANNTASMHRAA
>OPSD_CANFA
MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNVEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGMQCSCGIDYYTLKPEINNESFVIYMFVVHFAIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSDFGPIFMTLPAFFAKSSSIYNPVIYIMMNKQFRNCMITTLCCGKNPLGDDEASASASKTETSQVAPA
>BRS3_CAVPO
QTLISITNDTESSSSVVSNDTTNKGWTGDNSPGIEALCAIYITYAVIISVGILGNAILIKVFFKTKSMQTVPNIFITSLALGDLLLLLTCVPVDATHYLAEGWLFGRIGCKVLSFIRLTSVGVSVFTLTILSADRYKAVVKPLRQPSNAILKTCAKAGCIWIMSMIFALPEAIFSNVHTMT--SEWCAFY---VSEKLLQEIHALLSFLVFYIIPLSIISVYYSLIARTLYKSIPTQSHARKQVESRKRIAKTVLVLVALFALCWLPNHLLNLYHSFTHKAAIHFIVTIFSRVLAFSNSCVNPFALYWLSKTFQKQFKAQLFCCKGELPEPPLAATMGRVSGTENTHI
>OPSD_DIPAN
MNGTEGPFFYVPMVNTTGIVRSPYEYPQYYLVNPAAYAALGAYMFLLILVGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVLGGFTTTMYTSMHGYFVLGRLGCNIEGFFATLGGEIALWSLVVLAIERWVVVCKPISNFRFGENHAIMGLAFTWTMAMACAAPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFICHFTIPLTVVFFCYGRLLCAVKEAAAAQQESETTQRAEKEVTRMVIMMVIAFLVCWLPYASVAWYIFTHQGSEFGPVFMTIPAFFAKSSSIYNPMIYICLNKQFRHCMITTLCCGKNPFEEEEGASTSVSSSSVSPAA-
>AG2R_MOUSE
------MALNSSTEDGIKRIQDDC---PRAGRHSYIFVMIPTLYSIIFVVGIFGNSLVVIVIYFYMKLKTVASVFLLNLALADLCFLLTLPLWAVYTAMEYRWPFGNHLCKIASASVSFNLYASVFLLTCLSIDRYLAIVHPMSRLRRTMLVAKVTCIIIWLMAGLASLPAVIHRNVYF-TN--TVCAFH-YESRNSTLPIGLGLTKNILGFLFPFLIILTSYTLIWKALKKA----YEIQKNKPRNDDIFRIIMAIVLFFFFSWVPHQIFTFLDVLIQLGDIVDTAMPITICIAYFNNCLNPLFYGFLGKKFKKYFLQLLKYIPPKAKSHSSLSTRPSD-------N
>YMJC_CAEEL
NLRLNESPYKYVMSNNTTIPSCLTDRQMSLSVSSTEGVLIGTIIPILVLFGISGNILNLTVLLAPNL-RTRSNQLLACLAVADIVSLVVILPHSMAHYETFERKFYGKYKFQIIAMTNWSIATATWLVFVICLERLIIIKYPLPRNVVTIIVVTTFILTSYNHVSHACAEKLFCNGTQY--S-RWFRNEPPNSEFMKSVVRVAPQVNAIFVVLIPVVLVIIFNVMLILTLRQRKTISQFTQLQSKTEHKVTITVTAIVTCFTITQSPSAFVTFLSSYVH--RDWVTLSAICTILVVLGKALNFVLFCLSSASFRQRLLMQTKQGILRKSTRYTSVA------------
>OPSG_ASTFA
NE--ETTRESAFVYTNANNTRDPFEGPNYHIAPRWVYNLASLWMIIVVIASIFTNSLVIVATAKFKKLRHPLNWILVNLAIADLGETVLASTISVFNQVFGYFVLGHPMCIFEGWTVSVCGITALWSLTIISWERWVVVCKPFGNVKFDGKWAAGGIIFAWTWAIIWCTPPIFGWSRYWPHGLKTSCGPDVFSGSEDPGVASYMVTLLLTCCILPLSVIIICYIFVWNAIHQVAQQQKDSESTQKAEKEVSRMVVVMILAFILCWGPYASFATFSALNPGYAWHPLAAALPAYFAKSATIYNPIIYVFMNRQFRSCIMQLFGKKVEDA-----SEVSTAS--------
>DUFF_HUMAN
FEDVWNSSYGVNDSFPDGDYDANLEAAAPCHSCNLLDDSALPFFILTSVLGILASSTVLFMLFRPLFQLCPGWPVLAQLAVGSALFSIVVPVLAPGLGST--RSSALCSLGYCVWYGSAFAQALLLGCHASLGHRLGAGQ--------VPGLTLGLTVGIWGVAALLTLPVTLASG------SGGLCTLI-YSTELKALQATHTVACLAIFVLLPLGLFGAKGLKKA-------------------LGMGPGPWMNILWAWFIFWWPHGVVLGLDFLVRSKQALDLLLNLAEALAILHCVATPLLLALFCHQATRTLLPSLPLPEGWSSHLDTLGS------------
>OPR._RAT
GSHFQGNLSLLN---ETVPHHLLLNASHSAFLPLGLKVTIVGLYLAVCIGGLLGNCLVMYVILRHTKMKTATNIYIFNLALADTLVLLTLPFQGTDILLGF-WPFGNALCKTVIAIDYYNMFTSTFTLTAMSVDRYVAICHPIALDVRTSSKAQAVNVAIWALASVVGVPVAIMGSAQVEE---IECLVE-IPAPQDYWGPVFAICIFLFSFIIPVLIISVCYSLMIRRLRGVRLL-SGSREKDRNLRRITRLVLVVVAVFVGCWTPVQVFVLVQGLGVQPETAVAILRFCTALGYVNSCLNPILYAFLDENFKACFRKFCCASSLHREMQVSDRVGLGCKTSETVPR
>TRFR_SHEEP
------------MENETGSELNQTQLQPRAVVALEYQVVTILLVLIICGLGIVGNIMVVLVVMRTKHMRTPTNCYLVSLAVADLMVLVAAGLPNITDSIYGSWVYGYVGCLCITYLQYLGINASSCSITAFTIERYIAICHPIAQFLCTFSRAKKIIIFVWAFTSIYCMLWFFLLDLN--DA-SCGYKIS------RNYYSPIYLMDFGVFYVVPMILATVLYGFIARILFLSLNSNRYFNSTVSSRKQVTKMLAVVVILFALLWMPYRTLVVVNSFLSSPFQENWFLLFCRICIYLNSAINPVIYNLMSQKFRAAFRKLCNCKQKPVEKPANYSVKESDHFSTELDD
>ACTR_MESAU
-------------MKHIITPYEHTNDTARNNSDCPDVVLPEEIFFTISIIGVLENLIVLLAVVKNKNLQCPMYFFICSLAISDMLGSLYKILENILIMFRNRGNFESTADDIIDCMFILSLLGSIFSLSVIAADRYITIFHALYHSIVTMRRTIITLTVIWIFCTGSGIAMVIFFS----------------------HHHVPTVLTFTSLFPLMLVFILCLYIHMFLLARSH---ARKISTLPRANMKGAITLTILLGVFIFCWAPFILHVLLMTFCPNNVCYMSLFQINGMLIMCNAVIDPFIYAFRSPELRDAFKKMFSCHRYQ---------------------
>CKR7_HUMAN
QDEVTDDYIGDNTTVDYTLFESLC---SKKDVRNFKAWFLPIMYSIICFVGLLGNGLVVLTYIYFKRLKTMTDTYLLNLAVADILFLLTLPFWAYSAAKS--WVFGVHFCKLIFAIYKMSFFSGMLLLLCISIDRYVAIVQAVRHRARVLLISKLSCVGIWILATVLSIPELLYSDLQR-SE-AMRCSLI---TEHVEAFITIQVAQMVIGFLVPLLAMSFCYLVIIRTL---------LQARNFERNKAIKVIIAVVVVFIVFQLPYNGVVLAQTVANFNKQLNIAYDVTYSLACVRCCVNPFLYAFIGVKFRNDLFKLFKDLGCLSQEQLRQWS--------HIRR
>O.YR_SHEEP
GAFAANWSAEAVNGSAAPPGTEGNRTAGPPQRNEALARVEVAVLSLILFLALSGNACVLLALRTTRHKHSRLFFFMKHLSIADLAVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPL--RSLRRRTDRLAVLATWLGCLVASAPQVHIFSLRE----VFDCWAV---FIQPWGPKAYITWITLAVYIVPVIVLAACYGLISFKIWQNRAAVSNVKLISKAKIRTVKMTFIVVLAFIVCWTPFFFKQMWSVWDADA-KEASAFIIAMLLASLNSCCNPWIYMLFTGHLFQDLVQRFLCCSFRRLKGSQLGEHSYTFVLSRHSS
>GPR7_HUMAN
MDNASFSEPWPANASGPDPALSCSNASTLAPLPAPLAVAVPVVYAVICAVGLAGNSAVLYVLLRAPRMKTVTNLFILNLAIADELFTLVLPINIADFLLRQ-WPFGELMCKLIVAIDQYNTFSSLYFLTVMSADRYLVVLATARVAGRTYSAARAVSLAVWGIVTLVVLPFAVFARLD--QG-RRQCVLV-FPQPEAFWWRASRLYTLVLGFAIPVSTICVLYTTLLCRLHAMR--DSHAKALERAKKRVTFLVVAILAVCLLCWTPYHLSTVVALTTDLPPLVIAISYFITSLTYANSCLNPFLYAFLDASFRRNLRQLITCRAAA---------------------
>GPRJ_HUMAN
TETATPLPSQYLMELSEEHSWMSNQTDLHYVLKPGEVATASIFFGILWLFSIFGNSLVCLVIHRSRRTQSTTNYFVVSMACADLLISVASTPFVLLQFTTGRWTLGSATCKVVRYFQYLTPGVQIYVLLSICIDRFYTIVYPL-SFKVSREKAKKMIAASWIFDAGFVTPVLFFYG---SNW-HCNYFLP-----SSWEGTAYTVIHFLVGFVIPSVLIILFYQKVIKYIWRIDGR-RTMNIVPRTKVKTIKMFLILNLLFLLSWLPFHVAQLWHPHEQDYKKSSLVFTAITWISFSSSASKPTLYSIYNANFRRGMKETFCMSSMKCYRSNAYTIKKNYVGISEIPS
>OAJ1_HUMAN
--MLLCFRFGNQSMKRENFTLITDFVFQGFSSFHEQQITLFGVFLALYILTLAGNIIIVTIIRIDLHLHTPMYFFLSMLSTSETVYTLVILPRMLSSLVGMSQPMSLAGCATQMFFFVTFGITNCFLLTAMGYDRYVAICNPLYMVIMNKRLRIQLVLGACSIGLIVAITQVTSVFRL---PPHFFCDIRKLSCIDTTVNEILTLIISVLVLVVPMGLVFISYVLIISTI--------LKIASVEGRKKAFATCASHLTVVIVHYSCASIAYLKPKSEN---TREHDQLISVTYTVITPLLNPVVYTLRNKEVKDALCRAVGGKFS----------------------
>RDC1_CANFA
YAEPGNFSDISWPCNSSDCIVVDTVLCPNMPNKSVLLYTLSFIYIFIFVIGMIANSVVVWVNIQAKTTGYDTHCYILNLAIADLWVVVTIPVWVVSLVQHNQWPMGELTCKITHLIFSINLFGSIFFLTCMSVDRYLSITYFATSSRRKKVVRRAVCVLVWLLAFCVSLPDTYYLKTVTNNE--TYCRSFYPEHSVKEWLISMELVSVVLGFAIPFCVIAVFYCLLARAI---------SASSDQEKQSSRKIIFSYVVVFLVCWLPYHVVVLLDIFSILHNFLFTALHVTQCLSLVHCCVNPVLYSFINRNYRYELMKAFIFKYSAKTGLTKLIDEYSALEQNAK--
>CML1_RAT
EYEGYNDSSIYGEEYSDGSDYIVDLEEAGPLEAKVAEVFLVVIYSLVCFLGILGNGLVIVIATFKMK-KTVNTVWFVNLAVADFLFNIFLPIHITYAAMDYHWVFGKAMCKISSFLLSHNMYTSVFLLTVISFDRCISVLLPVSQNHRSVRLAYMTCVVVWVWLSSESPPSLVFGHVST-SF--HSTHPR-TDPVGYSRHVAVTVTRFLCGFLIPVFIITACYLTIVFKL---------QRNRQAKTKKPFKIIITIIITFFLCWCPYHTLYLLELHHTAVSVFSLGLPLATAVAIANSCMNPILYVFMGHDFKK-FKVALFSRLVNALSEDTGPSFTKMSSLIEKAS
>FSHR_EQUAS
--MMYSEFDYDLCNEVVDVTCSPKPDAFNPCEDIMGYDILRVLIWFISILAITGNIIVLVILITSQYKLTVPRFLMCNLAFADLCIGIYLLLIASVDIHTKSWQTG-AGCDAAGFFTVFGSELSVYTLTAITLERWHTITHAMLECKVQLRHAASVMLVGWIFGFGVGLLPIFGISTY----KVSICLPM-----IDSPLSQLYVMSLLVLNVLAFVVICGCYTHIYLTVRNP------NIVSSSSDTKIAKRMGILIFTDFLCMAPISFFGISASLKVALITVSKSKILLVLFYPINSCANPFLYAIFTKNFRRDFFILLSKFGCYEMQAQTYRTISHPKNGPCPPT
>OPSD_MULSU
MNGTEGPYFYIPMVNTTGIVRSPYDYPQYYLVNPAAYAALGAYMFFLILVGFPINFLTLYVTIEHKKLRTPLNYILLNLAVANLFMVFGGFTTTMYTSMHGYFVLGRLGCNLEGFFATLGGEIALWSLVVLAVERWMVVCKPISNFRFGENHAIMGLAMTWLMASACAVPPLVGWSRYIPEGMQCSCGVDYYTRAEGFNNESFVVYMFCCHFMIPLIIVFFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVVIMVIAFLVCWLPYASVAWWIFTHQGSEFGPVFMTIPAFFAKSSSIYNPMIYICMNKQFRNCMITTLCCGKNPFEEEEGASSSSVSSSSVSPAA
>NK3R_HUMAN
GNLSSSPSALGLPVASPAPSQPWANLTNQFVQPSWRIALWSLAYGVVVAVAVLGNLIVIWIILAHKRMRTVTNYFLVNLAFSDASMAAFNTLVNFIYALHSEWYFGANYCRFQNFFPITAVFASIYSMTAIAVDRYMAIIDPL-KPRLSATATKIVIGSIWILAFLLAFPQCLYSK---GR---TLCFVQ--WPEGPKQHFTYHIIVIILVYCFPLLIMGITYTIVGITLWGG--PCDKYHEQLKAKRKVVKMMIIVVMTFAICWLPYHIYFILTAIYQQLKYIQQVYLASFWLAMSSTMYNPIIYCCLNKRFRAGFKRAFRWCPFIKVSSY----TTRFHPNRQSSM
>CKR5_CERTO
----MDYQVSSPTYDIDYYTSEPC---QKINVKQIAARLLPPLYSLVFIFGFVGNILVVLILINCKRLKSMTDIYLLNLAISDLLFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWVVAVFASLPGIIFTRSQR-EG-HYTCSPHFPYSQYQFWKNFQTLKIVILGLVLPLLVMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCSIFA-------SSVY
>PI2R_MOUSE
DGHPGPPSVTPGSPLSAGGREWQGMAGSCWNITYVQDSVGPATSTLMFVAGVVGNGLALGILGARRRHPSAFAVLVTGLAVTDLLGTCFLSPAVFVAYARNSAHGGTMLCDTFAFAMTFFGLASTLILFAMAVERCLALSHPYYAQLDGPRCARFALPSIYAFCCLFCSLPLLGLGEHQQYCPGSWCFIR---MRSAQPGGCAFSLAYASLMALLVTSIFFCNGSVTLSLYHMRRHFVPTSRAREDEVYHLILLALMTVIMAVCSLPLMIRGFTQAIAPDS--REMGDLLAFRFNAFNPILDPWVFILFRKAVFQRLKFWLCCLCARSVHGDLQAPRRDPPAPTSLQA
>FSHR_CHICK
FGPVENEFDYGLCNEVVDFVCSPKPDAFNPCEDIMGYNVLRVLIWFINILAITGNTTVLIILISSQYKLTVPRFLMCNLAFADLCIGIYLLFIASVDIQTKSWQTG-AGCNAAGFFTVFASELSVYTLTVITLERWHTITYAMLNRKVRLRHAVIIMVFGWMFAFTVALLPIFGISSY----KVSICLPM-----IETPFSQAYVIFLLVLNVLAFVIICICYICIYFTVRNP------NVISSNSDTKIAKRMAILIFTDFLCMAPISFFAISASLRVPLITVSKSKILLVLFYPINSCANPFLYAIFTKTFRRDFFILLSKFGCCEMQAQIYRTNFHTRNGHYPTA
>OPRM_BOVIN
FSHLEGNLSDPCGPNRTELGGSDRLCPSAGSPSMITAIIIMALYSIVCVVGLFGNFLVMYVIVRYTKMKTATNIYIFNLALADALATSTLPFQSVNYLMGT-WPFGTILCKIVISIDYYNMFTSIFTLCTMSVDRYIAVCHPVALDLRTPRNAKIINICNWILSSAIGLPVMFMATTKYGS---IDCTLT-FSHPTWYWENLLKICVFIFAFIMPILIITVCYGLMILRLKSVRML-SGSKEKDRNLRRITRMVLVVVAVFIVCWTPIHIYVIIKALITIPTFQTVSWHFCIALGYTNSCLNPVLYAFLDENFKRCFREFCIPTSSTIEQQNSTRIPSTANTVDRTNH
>P2Y5_HUMAN
--------------------MVSVNSSHCFYNDSFKYTLYGCMFSMVFVLGLVSNCVAIYIFICVLKVRNETTTYMINLAMSDLLFVFTLPFRIFYFTTRN-WPFGDLLCKISVMLFYTNMYGSILFLTCISVDRFLAIVYPFSKTLRTKRNAKIVCTGVWLTVIGGSAPAVFVQSTHS---ASEACFENFPEATWKTYLSRIVIFIEIVGFFIPLILNVTCSSMVLKTLTKP----VTLSRSKINKTKVLKMIFVHLIIFCFCFVPYNINLILYSLVRTQAAVRTMYPITLCIAVSNCCFDPIVYYFTSDTIQNSIKMKNWSVRRSDFRFSEVHGNLQTLKSKIFDN
>HM74_HUMAN
----------MNRHHLQDHFLEIDKKNCCVFRDDFIAKVLPPVLGLEFIFGLLGNGLALWIFCFHLKSWKSSRIFLFNLAVADFLLIICLPFVMDYYVRRSDWNFGDIPCRLVLFMFAMNRQGSIIFLTVVAVDRYFRVVHPHALNKISNWTAAIISCLLWGITVGLTVHLLKKKLLIQ-PA--NVCISF-----SICHTFRWHEAMFLLEFLLPLGIILFCSARIIWSLRQR------QMDRHAKIKRAITFIMVVAIVFVICFLPSVVVRIRIFWLLHTRSVDLAFFITLSFTYMNSMLDPVVYYFSSPSFPNFFSTLINRCLQRKMTGEPDNNTGDPNKTRGAPE
>OPSD_DICLA
MNGTEGPFFYVPMVNTTGIVRSPYDYPQYYLVSPAAYAALGAYMFLLILLGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVLGRLGCNMEGFFATLGGEIGLWSLVVLAVERWLVVCKPISNFRFGENHAIMGLAFTWVMACSCAVPPLVGWSRYIPEGMQCSCGVDYYTRAEGFNNESFVIYMFACHFIIPMCVVFFCYGRLLCAVKEAAAAQQESETTQRAEKEVTRMVVIMGIAFLICWCPYASVAWYIFTHQGSEFGPVFMTLPAFFAKTSSVYNPLIYILMNKQFRHCMITTLCCGKNPFEEEEGASTSVSSSSVSPAA-
>GPRW_HUMAN
GTRGCSDRQPGVLTRDRSCSRKMNSSGCLSEEVGSLRPLTVVILSASIVVGVLGNGLVLWMTVFRMA-RTVSTVCFFHLALADFMLSLSLPIAMYYIVSRQ-WLLGEWACKLYITFVFLSYFASNCLLVFISVDRCISVLYPVALNHRTVQRASWLAFGVWLLAAALCSAHLKFRTTR-FNSNETAQIWI-----VVEGHIIGTIGHFLLGFLGPLAIIGTCAHLIRAKL---------LREGWVHANRPKRLLLVLVSAFFIFWSPFNVVLLVHLWRRV-PRMLLILQASFALGCVNSSLNPFLYVFVGRDFQEKFFQSLTSALARAFGEEEFLSPRE---------
>UL33_HSV6U
--------------------MDTVIELSKLQFKGNASCTSTPTLKTARIMESAVTGITLTTSIPMIIKHNATSFYVITLFASDFVLMWCVFFMTVNRKQL--FSFNRFFCQLVYFIYHAVCSYSISMLAIIATIRYK-TLHRRKKTESKTSSTGRNIGILLLASSMCAIPTALFVKTNGM-KKTGKCVVYSSKKAY-ELFLAVKIVFSFIWGVLPTMVFSFFYVIFCKAL---------HDVTEKKYKKTLFFIRILLLSFLLIQIPYIAILICEIAFLYMARVEILQLIIRLMPQVHCFSNPLVYAFTGGELRNRFTACFQSFFPKTLCSTQKRKDQNSKSKASVEK
>A1AB_HUMAN
HNTSAPAHWGELKNANFTGPNQTSSNSTLPQLDITRAISVGLVLGAFILFAIVGNILVILSVACNRHLRTPTNYFIVNLAMADLLLSFTVLPFSAALEVLGYWVLGRIFCDIWAAVDVLCCTASILSLCAISIDRYIGVRYSLYPTLVTRRKAILALLSVWVLSTVISIGPLLGWKEP--AP--KECGVT--------EEPFYALFSSLGSFYIPLAVILVMYCRVYIVAKRTHNPIAVKLFKFSREKKAAKTLGIVVGMFILCWLPFFIALPLGSLFSTLKPPDAVFKVVFWLGYFNSCLNPIIYPCSSKEFKRAFVRILGCQCRGRRRRRRRRRTYRPWTRGGSLE
>OPSD_SPHSP
GG-Y-GNQTVVDKVLPEMLHLIDPHWYQFPPMNPLWHGLLGFVIGCLGFVSVVGNGMVIYIFSTTKGLRTPSNLLVVNLAFSDFLMMLSMSPPMVINCYYETWVLGPFMCELYALLGSLFGCGSIWTMVMIALDRYNVIVKGLAAKPMTNKTAMLRILGIWAMSIAWTVFPLFGWNRYVPEGNMTACGTD--YLNKEWVSRSYILVYSVFVYFLPLATIIYSYWFIVQAVSAHMNVRSAENANTSAECKLAKVALMTISLWFFAWTPYLVTDFSGIFEWGK-ISPLATIWCSLFAKANAVYNPIVYGISHPKYRAALNKKFPSLACASEPDDTASQSDEKSASA----
>ET1R_RAT
EFNFLGTTLQPPNLALPSNGSMHGYCPQQTKITTAFKYINTVISCTIFIVGMVGNATLLRIIYQNKCMRNGPNALIASLALGDLIYVVIDLPINVFKLLAGRNDFGVFLCKLFPFLQKSSVGITVLNLCALSVDRYRAVASWSVQGIGIPLITAIEIVSIWILSFILAIPEAIGFVMVPRT--HRTCMLNATTKFMEFYQDVKDWWLFGFYFCMPLVCTAIFYTLMTCEMLNRGSL-IALSEHLKQRREVAKTVFCLVVIFALCWFPLHLSRILKKTVYDESFLLLMDYIGINLATMNSCINPIALYFVSKKFKNCFQSCLCCCCHQSKSLMTSVPWKNQEQN-HNTE
>B1AR_PIG
LPDGAATAARLLVPASPPASLLTPASEGSVQLSQQWTAGMGLLMALIVLLIVAGNVLVIVAIAKTPRLQTLTNLFIMSLASADLVMGLLVVPFGATIVVWGRWEYGSFFCELWTSVDVLCVTASIETLCVIALDRYLAITSPFYQSLLTRAAR-ALVCTVWAISALVSFLPILMHWWRDR-R--KCCDFV--------TNRAYAIASSVVSFYVPLCIMAFVYLRVFREAQKQNGRRRPSRLVALREQKALKTLGIIMGVFTLCWLPFFLANVVKAFHRDL-VPDRLFVFFNWLGYANSAFNPIIYCRSP-DFRKAFQRLLCCARRVARGSCAAAGCLAVARPPPSPG
>OLF3_RAT
-------------MDSSNRTRVSEFLLLGFVENKDLQPLIYGLFLSMYLVTVIGNISIIVAIISDPCLHTPMYFFLSNLSFVDICFISTTVPKMLVNIQTQNNVITYAGCITQIYFFLLFVELDNFLLTIMAYDRYVAICHPMYTVIMNYKLCGFLVLVSWIVSVLHALFQSLMMLAL---PPHYFCEPNQLTCSDAFLNDLVIYFTLVLLATVPLAGIFYSYFKIVSSI--------CAISSVHGKYKAFSTCASHLSVVSLFYCTGLGVYLSSAANN---SSQASATASVMYTVVTPMVNPFIYSLRNKDVKSVLKKTLCEEVIRSPPSLLHFFCFIFCY------
>5H5A_HUMAN
LPVNLTSFSLSTPSPLETNHSGKDDLRPSSPLLSVFGVLILTLLGFLVAATFAWNLLVLATILRVRTFHRVPHNLVASMAVSDVLVAALVMPLSLVHELSGRWQLGRRLCQLWIACDVLCCTASIWNVTAIALDRYWSITRHMYTLRTRKCVSNVMIALTWALSAVISLAPLLFGWGE-S-E-SEECQVS--------REPSYAVFSTVGAFYLPLCVVLFVYWKIYKAAKFRATVPEGDTWREQKEQRAALMVGILIGVFVLCWIPFFLTELISPLCSC-DIPAIWKSIFLWLGYSNSFFNPLIYTAFNKNYNSAFKNFFSRQH-----------------------
>OPRD_RAT
FSLLANVSDTFPSAFPSASANASGSPGARSASSLALAIAITALYSAVCAVGLLGNVLVMFGIVRYTKLKTATNIYIFNLALADALATSTLPFQSAKYLMET-WPFGELLCKAVLSIDYYNMFTSIFTLTMMSVDRYIAVCHPVALDFRTPAKAKLINICIWVLASGVGVPIMVMAVTQPGA---VVCTLQ-FPSPSWYWDTVTKICVFLFAFVVPILIITVCYGLMLLRLRSVRLL-SGSKEKDRSLRRITRMVLVVVGAFVVCWAPIHIFVIVWTLVDINPLVVAALHLCIALGYANSSLNPVLYAFLDENFKRCFRQLCRAPCGGQEPGSLRRPRVTACTPSDGPG
>CKRV_MOUSE
EIPAVTEPSYNTVAKNDFMSGFLC---FSINVRAFGITVPTPLYSLVFIIGVIGHVLVVLVLIQHKRLRNMTSIYLFNLAISDLVFLSTLPFWVDYIMKGD-WIFGNAMCKFVSGFYYLGLYSDMFFITLLTIDRYLAVVHVVALRARTVTFGIISSIITWVLAALVSIPCLYVFKSQM-EF-YHTCRAILPRKSLIRFLRFQALTMNILGLILPLLAMIICYTRIINVL---------HRRPNKKKAKVMRLIFVITLLFFLLLAPYYLAAFVSAFEDVLQQVDLSLMITEALAYTHCCVNPVIYVFVGKRFRKYLWQLFRRHTAITLPQWLPFLA-------SARL
>THRR_RAT
PLEGRAVYLNKSRFPPMPPPPFISEDASGYLTSPWLTLFIPSVYTFVFIVSLPLNILAIAVFVFRMKVKKPAVVYMLHLAMADVLFVSVLPFKISYYFSGTDWQFGSGMCRFATAACYCNMYASIMLMTVISIDRFLAVVYPISLSWRTLGRANFTCVVIWVMAIMGVVPLLLKEQTTQ--N--TTCHDVLNETLLHGFYSYYFSAFSAIFFLVPLIISTVCYTSIIRCL------SSSAVANRSKKSRALFLSAAVFCIFIVCFGPTNVLLIVHYLLLSDETAYFAYLLCVCVTSVASCIDPLIYYYASSECQKHLYSILCCRESSDSNSCNSTGDTCS--------
>DADR_PIG
--------------MRTLNTSTMDGTGLVVERDFSFRILTACFLSLLILSTLLGNTLVCAAVIRFRHRSKVTNFFVISLAVSDLLVAVLVMPWKAVAEIAGFWPFG-SFCNIWVAFDIMCSTASILNLCVISVDRYWAISSPFYERKMTPKAAFILISVAWTLSVLISFIPVQLSWHKAGN-TTHNCDSS--------LSRTYAISSSLISFYIPVAIMIVTYTRIYRIAQKQAECESSFKMSFKRETKVLKTLSVIMGVFVCCWLPFFILNCMVPFCGSGCIDSITFDVFVWFGWANSSLNPIIYAFNA-DFRKAFSTLLGCYRLCPTSTNAIETAVVFSSH-----
>GRHR_BOVIN
--MANSDSPEQNENHCSAINSSIPLTPGSLPTLTLSGKIRVTVTFFLFLLSTIFNTSFLLKLQNWTQKLSRMKLLLKHLTLANLLETLIVMPLDGMWNITVQWYAGELLCKVLSYLKLFSMYAPAFMMVVISLDRSLAITKPL-AVKSNSKLGQFMIGLAWLLSSIFAGPQLYIFGMIHEG--FSQCVTH--SFPQWWHQAFYNFFTFSCLFIIPLLIMVICNAKIIFTLTRVPHKNQSKNNIPRARLRTLKMTVAFATSFTVCWTPYYVLGIWYWFDPDMRVSDPVNHFFFLFAFLNPCFDPLIYGYFSL-------------------------------------
>OLF3_CHICK
-------------MASGNCTTPTTFILSGLTDNPGLQMPLFMVFLAIYTITLLTNLGLIRLISVDLHLQTPMYIFLQNLSFTDAAYSTVITPKMLATFLEERKTISYVGCILQYFSFVLLTTSECLLLAVMAYDRYVAICKPLYPAIMTKAVCWRLVESLYFLAFLNSLVHTCGLLKL---SNHFFCDISQISSSSIAISELLVIISGSLFVMSSIIIILISYVFIILTV--------VMIRSKDGKYKAFSTCTSHLMAVSLFHGTVIFMYLRPVKLF---SLDTDKIASLFYTVVIPMLNPLIYSWRNKEVKDALRRLTATTFGFIDSKAVQ--------------
>ML1A_PHOSU
-------MKGNGSTLLNASQQAPGVGEGGGPRPSWLASTLAFILIFTIVVDILGNLLVILSVYRNKKLRNAGNIFVVSLAIADLVVAIYPYPLVLTSIFNNGWNLGYLHCQISAFLMGLSVIGSIFNITGIAINRYCYICHSLYDRLYSNKNSLCYVFLIWVLTLVAIMPNLQTGT-LQYDP-IYSCTFT------QSVSSAYTIAVVVFHFIVPMIIVIFCYLRIWILVLQVRR-PDSKPRLKPQDFRNFVTMFVVFVLFAICWAPLNFIGLIVASDPATRIPEWLFVASYYMAYFNSCLNAIIYGLLNQNFRQEYKRILVSLFTAKMCFVDSSNCKPAPLIANNNL
>NY2R_HUMAN
EEMKVEQYGP-QTTPRGELVPDPEPELIDSTKLIEVQVVLILAYCSIILLGVIGNSLVIHVVIKFKSMRTVTNFFIANLAVADLLVNTLCLPFTLTYTLMGEWKMGPVLCHLVPYAQGLAVQVSTITLTVIALDRHRCIVYHL-ESKISKRISFLIIGLAWGISALLASPLAIFREYSLFE--IVACTEKWPGEEKSIYGTVYSLSSLLILYVLPLGIISFSYTRIWSKLKNHSPG-AANDHYHQRRQKTTKMLVCVVVVFAVSWLPLHAFQLAVDIDSQVKEYKLIFTVFHIIAMCSTFANPLLYGWMNSNYRKAFLSAFRCEQRLDAIHSE---AKKNLEVRKNSG
>FMLR_HUMAN
-----------METNSSLPTNISGGTPAVSAGYLFLDIITYLVFAVTFVLGVLGNGLVIWVAGFRMT-HTVTTISYLNLAVADFCFTSTLPFFMVRKAMGGHWPFGWFLCKFLFTIVDINLFGSVFLIALIALDRCVCVLHPVTQNHRTVSLAKKVIIGPWVMALLLTLPVIIRVTTVPSPWPKERINVA------VAMLTVRGIIRFIIGFSAPMSIVAVSYGLIATKI---------HKQGLIKSSRPLRVLSFVAAAFFLCWSPYQVVALIATVRIREKEIGIAVDVTSALAFFNSCLNPMLYVFMGQDFRERLIHALPASLERALTEDSTQTTL----------
>NK1R_CAVPO
-----MDNVLPVDSDLFPNISTNTSEPNQFVQPAWQIVLWAAAYTVIVVTSVVGNVVVMWIILAHKRMRTVTNYFLVNLAFAEASMAAFNTVVNFTYAVHNEWYYGLFYCKFHNFFPIAAVFASIYSMTAVAFDRYMAIIHPL-QPRLSATATKVVICVIWVLALLLAFPQGYYST---GR---VVCMIEWPSHPDKIYEKVYHICVTVLIYFLPLLVIGYAYTVVGITLE---IPSDRYHEQVSAKRKVVKMMIVVVCTFAICWLPFHIFFLLPYINPDLKFIQQVYLAIMWLAMSSTMYNPIIYCCLNDRFRLGFKHAFRCCPFISAADY----STRYFQTQGSVY
>ML1B_HUMAN
NGSFANCCEAGGWAVRPGWSGAGSARPSRTPRPPWVAPALSAVLIVTTAVDVVGNLLVILSVLRNRKLRNAGNLFLVSLALADLVVAFYPYPLILVAIFYDGWALGEEHCKASAFVMGLSVIGSVFNITAIAINRYCYICHSMYHRIYRRWHTPLHICLIWLLTVVALLPNFFVGS-LEYDP-IYSCTFI------QTASTQYTAAVVVIHFLLPIAVVSFCYLRIWVLVLQARK-PESRLCLKPSDLRSFLTMFVVFVIFAICWAPLNCIGLAVAINPQEQIPEGLFVTSYLLAYFNSCLNAIVYGLLNQNFRREYKRILLALWNPRHCIQDASKQSPAPPIIGVQH
>C3AR_RAT
--------------MESFTADTNSTDLHSRPLFKPQDIASMVILSLTCLLGLPGNGLVLWVAGVKMK-RTVNTVWFLHLTLADFLCCLSLPFSVAHLILRGHWPYGLFLCKLIPSVIILNMFASVFLLTAISLDRCLMVHKPICQNHRSVRTAFAVCGCVWVVTFVMCIPVFVYRDLLV-ED-DYFDQLM-YGNHAWTPQVAITISRLVVGFLVPFFIMITCYSLIVFRM--------RKTNLTKSRNKTLRVAVAVVTVFFVCWIPYHIVGILLVITDQEEVVLPWDHMSIALASANSCFNPFLYALLGKDFRKKARQSVKGILEAAFSEELTHSAPS---------
>OPSB_GECGE
MNGTEGINFYVPLSNKTGLVRSPFEYPQYYLADPWKFKVLSFYMFFLIAAGMPLNGLTLFVTFQHKKLRQPLNYILVNLAAANLVTVCCGFTVTFYASWYAYFVFGPIGCAIEGFFATIGGQVALWSLVVLAIERYIVICKPMGNFRFSATHAIMGIAFTWFMALACAGPPLFGWSRFIPEGMQCSCGPDYYTLNPDFHNESYVIYMFIVHFTVPMVVIFFSYGRLVCKVREAAAQQQESATTQKAEKEVTRMVILMVLGFLLAWTPYAATAIWIFTNRGAAFSVTFMTIPAFFSKSSSIYNPIIYVLLNKQFRNCMVTTICCGKNPFGDEDVSSSVSSVSSSQVAPA
>SSR5_RAT
LSLASTPSWNAS---AASSGNHNWSLVGSASPMGARAVLVPVLYLLVCTVGLSGNTLVIYVVLRHAKMKTVTNVYILNLAVADVLFMLGLPFLATQNAVVSYWPFGSFLCRLVMTLDGINQFTSIFCLMVMSVDRYLAVVHPLSARWRRPRVAKMASAAVWVFSLLMSLPLLVFADVQE-----GTCNLS-WPEPVGLWGAAFITYTSVLGFFGPLLVICLCYLLIVVKVKAAMR--VGSSRRRRSEPKVTRMVVVVVLVFVGCWLPFFIVNIVNLAFTLPPTSAGLYFFVVVLSYANSCANPLLYGFLSDNFRQSFRKVLCLRRGYGMEDADAIERPQATLPTRSCE
>B2AR_CANFA
MGQPANRSVFLLAPNGSHAPDQ----GDSQERSEAWVVGMGIVMSLIVLAIVFGNVLVITAIARFERLQTVTNYFITSLACADLVMGLAVVPFGASHILMKMWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFYQSLLTKNKARVVILMVWIVSGLTSFLPIQMHWYRAI-N--TCCDFF--------TNQAYAIASSIVSFYLPLVVMVFVYSRVFQVAQRQGRSHRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNL-IPKEVYILLNWVGYVNSAFNPLIYCRSP-DFRIAFQELLCLRRSSLKAYGNGYSDYAGEHSGCHLG
>CKR6_MOUSE
FGTDDYDN---TEYYSIPPDHGPC---SLEEVRNFTKVFVPIAYSLICVFGLLGNIMVVMTFAFYKKARSMTDVYLLNMAITDILFVLTLPFWAVTHATNT-WVFSDALCKLMKGTYAVNFNCGMLLLACISMDRYIAIVQATRVRSRTLTHSKVICVAVWFISIIISSPTFIFNKKYE-LQ-RDVCEPRRSVSEPITWKLLGMGLELFFGFFTPLLFMVFCYLFIIKTL---------VQAQNSKRHRAIRVVIAVVLVFLACQIPHNMVLLVTAVNTGKKVLAYTRNVAEVLAFLHCCLNPVLYAFIGQKFRNYFMKIMKDVWCMRRKNKMPGFESY-----ISRQ
>NYR_DROME
TLSGLQFETYNITVMMNFSCDDYDLLSEDMWSSAYFKIIVYMLYIPIFIFALIGNGTVCYIVYSTPRMRTVTNYFIASLAIGDILMSFFCEPSSFISLFILNWPFGLALCHFVNYSQAVSVLVSAYTLVAISIDRYIAIMWPL-KPRITKRYATFIIAGVWFIALATALPIPIVSGLDI---WHTKCEKYREMWPSRSQEYYYTLSLFALQFVVPLGVLIFTYARITIRVWAKPGETNRDQRMARSKRKMVKMMLTVVIVFTCCWLPFNILQLLLNDEEFADPLPYVWFAFHWLAMSHCCYNPIIYCYMNARFRSGFVQLMHRMPGLRRWCCLRSVSGTGPALPLNRM
>HH2R_MOUSE
-------------------MEPNGTVHSCCLDSIALKVTISVVLTTLIFITVAGNVVVCLAVSLNRRLRSLTNCFIVSLAATDLLLGLLVMPFSAIYQLSFKWRFGQVFCNIYTSLDVMLCTASILNLFMISLDRYCAVTDPLYPVLVTPVRVAISLVFIWVISITLSFLSIHLGWN--RN-TF-KCKVQ--------VNEVYGLVDGMVTFYLPLLIMCVTYYRIFKIAREQR--ISSWKAATIREHKATVTLAAVMGAFIVCWFPYFTAFVYRGLRGDD-VNEVVEGIVLWLGYANSALNPILYATLNRDFRMAYQQLFHCKLASHNSHKTSLRRSQSREGRW---
>TRFR_BOVIN
-----------MENETGSELN-QTQLQPRAVVALEYQVVTILLVLIICGLGIVGNIMVVLVVMRTKHMRTPTNCYLVSLAVADLMVLVAAGLPNITDSIYGSWVYGYVGCLCITYLQYLGINASSCSITAFTIERYIAICHPIAQFLCTFSRAKKIIIFVWAFTSIYCMLWFFLLDLN--DA-SCGYKIS------RNYYSPIYLMDFGVFYVVPMILATVLYGFIARILFLNLNSNRYFNSTVSSRKQVTKMLAVVVILFALLWMPYRTLVVVNSFLSSPFQENWFLLFCRICIYLNSAINPVIYNLMSQKFRAAFRKLCNCKQKPVEKPANYSVKESDRFSTELDD
>OPSP_ICTPU
-------MASIILINFSETDTLHLGSVNDHIMPRIGYTILSIIMALSSTFGIILNMVVIIVTVRYKQLRQPLNYALVNLAVADLGCPVFGGLLTAVTNAMGYFSLGRVGCVLEGFAVAFFGIAGLCSVAVIAVDRYMVVCRPLGAVMFQTKHALAGVVFSWVWSFIWNTPPLFGWGSYQLEGVMTSCAPN--WYRRDPVNVSYILCYFMLCFALPFATIIFSYMHLLHTLWQVAKLVADSGSTAKVEVQVARMVVIMVMAFLLTWLPYAAFALTVIIDSNIYINPVIGTIPAYLAKSSTVFNPIIYIFMNRQFRDYALPCLLCGKNPWAAKEGRDSTVSKNTSVSPL-
>OPSR_.ENLA
NDDDDTTRSSVFTYTNSNNTRGPFEGPNYHIAPRWVYNLTSIWMIFVVFASVFTNGLVIVATLKFKKLRHPLNWILVNMAIADLGETVIASTISVFNQIFGYFILGHPMCVLEGFTVSTCGITALWSLTVIAWERWFVVCKPFGNIKFDEKLAATGIIFSWVWSAGWCAPPMFGWSRFWPHGLKTSCGPDVFSGSSDPGVQSYMLVLMITCCIIPLAIIILCYLHVWWTIRQVAQQQKESESTQKAEREVSRMVVVMIVAYIFCWGPYTFFACFAAFSPGYSFHPLAAALPAYFAKSATIYNPIIYVFMNRQFRNCIYQMFGKKVDDG-----SEVVSSVSNSSVSPA
>ML1._MOUSE
---------MGPTKAVPTPFGCIGCKLPKPDYPPALIIFMFCAMVITVVVDLIGNSMVILAVTKNKKLRNSGNIFVASLSVADMLVAIYPYPLMLYAMSVGGWDLSQLQCQMVGLVTGLSVVGSIFNITAIAINRYCYICHSLYKRIFSLRNTCIYLVVTWVMTVLAVLPNMYIGT-IEYDP-TYTCIFN------YVNNPAFTVTIVCIHFVLPLIIVGYCYTKIWIKVLAARD-AGQNPDNQFAEVRNFLTMFVIFLLFAVCWCPVNVLTVLVAVIPKEKIPNWLYLAAYCIAYFNSCLNAIIYGILNESFRREYWTIFHAMRHPILFISHLISTRALTRARVRAR
>OPSG_CHICK
MNGTEGINFYVPMSNKTGVVRSPFEYPQYYLAEPWKYRLVCCYIFFLISTGLPINLLTLLVTFKHKKLRQPLNYILVNLAVADLFMACFGFTVTFYTAWNGYFVFGPVGCAVEGFFATLGGQVALWSLVVLAIERYIVVCKPMGNFRFSATHAMMGIAFTWVMAFSCAAPPLFGWSRYMPEGMQCSCGPDYYTHNPDYHNESYVLYMFVIHFIIPVVVIFFSYGRLICKVREAAAQQQESATTQKAEKEVTRMVILMVLGFMLAWTPYAVVAFWIFTNKGADFTATLMAVPAFFSKSSSLYNPIIYVLMNKQFRNCMITTICCGKNPFGDEDVSSTVSSVSSSQVSPA
>GRHR_HUMAN
--MANSASPEQNQNHCSAINNSIPLMQGNLPTLTLSGKIRVTVTFFLFLLSATFNASFLLKLQKWTQKLSRMKLLLKHLTLANLLETLIVMPLDGMWNITVQWYAGELLCKVLSYLKLFSMYAPAFMMVVISLDRSLAITRPL-ALKSNSKVGQSMVGLAWILSSVFAGPQLYIFRMIHKV--FSQCVTH--SFSQWWHQAFYNFFTFSCLFIIPLFIMLICNAKIIFTLTRVPHENQSKNNIPRARLKTLKMTVAFATSFTVCWTPYYVLGIWYWFDPEMRLSDPVNHFFFLFAFLNPCFDPLIYGYFSL-------------------------------------
>5HT2_APLCA
----------------MLCGRLRHTMNSTTCFFSHRTVLIGIVGSLIIAVSVVGNVLVCLAIFTEPISHSKSKFFIVSLAVADLLLALLVMTFALVNSLYGYWLFGETFCFIWMSADVMCETASIFSICVISYNRLKQVQKPLYEEFMTTTRALLIIASLWICSFVVSFVPFFLEWHELGD-PKPECLFD--------VHFIYSVIYSLFCFYIPCTLMLRNYLRLFLIAKKH-RIHRLHRNQGTQGSKAARTLTIITGTFLACWLPFFIINPIEAVDEHL-IPLECFMVTIWLGYFNSCVNPIIYGTSNSKFRAAFQRLLRCRSVKSTVSSISPVSWIRPSLLDGP-
>OPSD_ASTFA
MNGTEGPYFYVPMSNATGVVRSPYEYPQYYLAPPWAYACLAAYMFFLILVGFPVNFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSLNGYFVFGRLGCNLEGFFATFGGINSLWCLVVLSIERWVVVCKPMSNFRFGENHAIMGVAFTWFMALACTVPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFVVHFLTPLFVITFCYGRLVCTVKEAAAQQQESETTQRAEREVTRMVILMFIAYLVCWLPYASVSWWIFTNQGSEFGPIFMTVPAFFAKSSSIYNPVIYICLNKQFRHCMITTLCCGKNPF--EEEEGATEASSVSSVSPA
>5H2C_HUMAN
GLLVWQCDISVSPVAAIVTDIFNTSDGGRFKFPDGVQNWPALSIVIIIIMTIGGNILVIMAVSMEKKLHNATNYFLMSLAIADMLVGLLVMPLSLLAILYDYWPLPRYLCPVWISLDVLFSTASIMHLCAISLDRYVAIRNPIHSRFNSRTKAIMKIAIVWAISIGVSVPIPVIGLRD-VF--NTTCVL---------NDPNFVLIGSFVAFFIPLTIMVITYCLTIYVLRRQKKKPRGTMQAINNERKASKVLGIVFFVFLIMWCPFFITNILSVLCEKSKLMEKLLNVFVWIGYVCSGINPLVYTLFNKIYRRAFSNYLRCNYKVEKKPPVRQIALSGRELNVNIY
>DBDR_RAT
PGRNRTAQPARLGLQRQLAQVDAPAG--SATPLGPAQVVTAGLLTLLIVWTLLGNVLVCAAIVRSRHRAKMTNIFIVSLAVSDLFVALLVMPWKAVAEVAGYWPFG-TFCDIWVAFDIMCSTASILNLCIISVDRYWAISRPFYERKMTQRVALVMVGLAWTLSILISFIPVQLNWHRDEG-RTENCDSS--------LNRTYAISSSLISFYIPVAIMIVTYTRIYRIAQVQRGADPSLRASIKKETKVFKTLSMIMGVFVCCWLPFFILNCMVPFCSSGCVSETTFDIFVWFGWANSSLNPIIYAFNA-DFRKVFAQLLGCSHFCFRTPVQTVNYNQDTVFHK---
>EDG1_MOUSE
VKALRSSVSDYGNYDIIVRHYNYTGKLNIGAEKDHGIKLTSVVFILICCFIILENIFVLLTIWKTKKFHRPMYYFIGNLALSDLLAG-VAYTANLLLSGATTYKLTPAQWFLREGSMFVALSASVFSLLAIAIERYITMLKMKLHNGSNSSRSFLLISACWVISLILGGLPSMGWNCI-S--SLSSCSTV------LPLYHKHYILFCTTVFTLLLLSIAILYCRIYSLVRTRLTFISKGSRSSEKSLALLKTVIIVLSVFIACWAPLFILLLLDVGCKAKCDILYKAEYFLVLAVLNSGTNPIIYTLTNKEMRRAFIRIVSCCKCPNGDSAGKFKEFSRSKSDNSSH
>B3AR_BOVIN
MAPWPPGNSSLTPWPDIPTLAPNTANASGLPGVPWAVALAGALLALAVLATVGGNLLVIVAIARTPRLQTMTNVFVTSLATADLVVGLLVVPPGATLALTGHWPLGVTGCELWTSVDVLCVTASIETLCALAVDRYLAVTNPLYGALVTKRRALAAVVLVWVVSAAVSFAPIMSKWWRIQ-R--RCCTFA--------SNMPYALLSSSVSFYLPLLVMLFVYARVFVVATRQGVPRRPARLLPLREHRALRTLGLIMGTFTLCWLPFFVVNVVRALGGPS-VSGPTFLALNWLGYANSAFNPLIYCRSP-DFRSAFRRLLCRCR---PEEHLAAAGAPTALTSPAGP
>NY2R_CAVPO
EEIKVEPYGPGHTTPRGELAPDPEPELIDSTKLTEVRVVLILAYCSIILLGVVGNSLVIHVVIKFKSMRTVTNFFIANLAVADLLVNTLCLPFTLTYTLMGEWKMGPVLCHLVPYAQGLAVQVSTVTLTVIALDRHRCIVYHL-DSKISKQNSFLIIGLAWGISALLASPLAIFREYSLFE--IVACTEKWPGEEKSIYGTVYSLSSLLILYVLPLGIISVSYVRIWSKLKNHSPG-AANDHYHQRRQKTTKMLVFVVVVFAVSWLPLHAFQLAVDIDSQVKEYKLIFTVFHIIAMCSTFANPLLYGWMNSNYRKAFLSAFRCQQRLDAIQSE---AKTNVEVEKNHG
>PE22_RAT
---------------MDNSFNDSRRVENCESRQYLLSDESPAISSVMFTAGVLGNLIALALLARRWRSISLFHVLVTELVLTDLLGTCLISPVVLASYSRNQLAPESRACTYFAFTMTFFSLATMLMLFAMALERYLAIGHPYYRRRVSRRGGLAVLPAIYGVSLLFCSLPLLNYGEYVQYCPGTWCFIQ--------HGRTAYLQLYATVLLLLIVAVLGCNISVILNLIRMRGPRRGERTSMAEETDHLILLAIMTITFAVCSLPFTIFAYMDETSS---RKEKWDLRALRFLSVNSIIDPWVFVILRPPVLRLMRSVLCCRTSLRAPEAPGAS--QTDLCGQL--
>A2AB_MOUSE
--------------------MSGPAMVHQEPYSVQATAAIASAITFLILFTIFGNALVILAVLTSRSLRAPQNLFLVSLAAADILVATLIIPFSLANELLGYWYFWRAWCEVYLALDVLFCTSSIVHLCAISLDRYWAVSRALYNSKRTPRRIKCIILTVWLIAAVISLPPLIYKGD-Q-----PQCELN--------QEAWYILASSIGSFFAPCLIMILVYLRIYVIAKRSGVAWWRRRTQLSREKRFTFVLAVVIGVFVVCWFPFFFSYSLGAICPQHKVPHGLFQFFFWIGYCNSSLNPVIYTIFNQDFRRAFRRILCRQWTQTGW------------------
>AG2S_HUMAN
------MILNSSTEDGIKRIQDDC---PKAGRHNYIFVMIPTLYSIIFVVGIFGNSLVVIVIYFYMKLKTVASVFLLNLALADLCFLLTLPLWAVYTAMEYRWPFGNYLCKIASASVSFNLYASVFLLTCLSIDRYLAIVHPMSRLRRTMLVAKVTCIIIWLLAGLASLPAIIHRNVFF-TN--TVCAFH-YESRNSTLPIGLGLTKNILGSCFPFLIILTSYTLIWKALKKA----YEIQKNNPRNDDIFRIIMAIVLFFFFSWIPHQIFTFLDVLIQQGDIVDTAMPITIWIAYFNNCLNPLFYGFLGKKFKKDILQLLKYIPPKAKSHSNLSTRPSD-------N
>B2AR_BOVIN
MGQPGNRSVFLLAPNASHAPDQ----NVTLERDEAWVVGMGILMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGACHILMKMWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYLAITSPFYQCLLTKNKARVVILMVWIVSGLTSFLPIQMHWYRAI-N--TCCDFF--------TNQPYAIASSIVSFYLPLVVMVFVYSRVFQVAKRQGRSQRRTSKFYLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIKDNL-IRKEIYILLNWLGYINSAFNPLIYCRSP-DFRIAFQELLCLRRSSLKAYGNGCSDYTGEQSGYHLG
>IL8B_RABIT
G--DFSNYSYSTDLPPTLLDSAPC----RSESLETNSYVVLITYILVFLLSLLGNSLVMLVILYSRSTCSVTDVYLLNLAIADLLFATTLPIWAASKVHG--WTFGTPLCKVVSLVKEVNFYSGILLLACISVDRYLAIVHATRTMIQKRHLVKFICLSMWGVSLILSLPILLFRNAIF-NS-SPVCYED-MGNSTAKWRMVLRILPQTFGFILPLLVMLFCYVFTLRTL---------FQAHMGQKHRAMRVIFAVVLIFLLCWLPYNLVLLTDTLMRTHNDIDRALDATEILGFLHSCLNPIIYAFIGQKFRYGLLKILAAHGLISKEFLAKES------------
>GASR_RAT
GSSLCRPGVSLLNSSSAGNLSCDPPRIRGTGTRELEMAIRITLYAVIFLMSVGGNVLIIVVLGLSRRLRTVTNAFLLSLAVSDLLLAVACMPFTLLPNLMGTFIFGTVICKAISYLMGVSVSVSTLNLVAIALERYSAICRPLARVWQTRSHAARVILATWLLSGLLMVPYPVYTMVQP---V-LQCMHR---WPSARVQQTWSVLLLLLLFFIPGVVIAVAYGLISRELYLGPGPPRPNQAKLLAKKRVVRMLLVIVLLFFLCWLPVYSVNTWRAFDGPGALSGAPISFIHLLSYVSACVNPLVYCFMHRRFRQACLDTCARCCPRPPRARPQPLPSIASLSRLSYT
>A1AB_MOUSE
HNTSAPAHWGELKDANFTGPNQTSSNSTLPQLDVTRAISVGCLG-AFILFAIVGNILVILSVACNRHLRTPTNYFIVNLAIADLLLSFTDLPFSATLEVLGYWVLGRIFCDIWAAVDVLCCTASILSLCAISIDRYIGVRYSLYPTLVTRRKAILALLSVWVLSTVISIGPLLGWKEP--AP--KECGVT--------EEPFYALFSSLGSFYIPLAVILVMYCRVYIVAKRTHNPIAVKLFKFSREKKAAKTLGIVVGMFILCWLPFFIALPLGSLFSTLKPPDAVFKVVFWLGYFNSCLNPIIYPCSSKEFKRAFMRILGCQCRGGRRRRRRRRTYRPWTRGGSLE
>OPSR_FELCA
AGLEDSTRASIFTYTNSNATRGPFEGPNYHIAPRWVYHVTSAWMIFVVIASVFTNGLVLAATMKFKKLRHPLNWILVNLAVADLAETIIASTISVVNQIYGYFVLGHPMCVLEGYTVSLCGITGLWSLAIISWERWLVVCKPFGNVRFDAKLAIAGIAFSWIWAAVWTAPPIFGWSRYWPHGLKTSCGPDVFSGSSYPGVQSYMIVLMITCCIIPLSVIVLCYLQVWLAIRAVAKQQKESESTQKAEKEVTRMVMVMIFAYCVCWGPYTFFACFAAAHPGYAFHPLVAALPAYFAKSATIYNPIIYVFMNRQFRNCIMQLFGKKVDDG-----SELASSV--SSVSPA
>P2Y7_HUMAN
-------------------MNTTSSAAPPSLGVEFISLLAIILLSVALAVGLPGNSFVVWSILKRMQKRSVTALMVLNLALADLAVLLTAPFFLHFLAQGT-WSFGLAGCRLCHYVCGVSMYASVLLITAMSLDRSLAVARPFSQKLRTKAMARRVLAGIWVLSFLLATPVLAYRTVVPK--NMSLCFPR---YPSEGHRAFHLIFEAVTGFLLPFLAVVASYSDIGRRL---------QARRFRRSRRTGRLVVLIILTFAAFWLPYHVVNLAEAGRALAKRLSLARNVLIALAFLSSSVNPVLYACAGGGLLRSAGVGFVAKLLEGTGSEASSTQTARSGPAALEP
>GPRL_HUMAN
---------MNSTLDGNQSSHPFCLLAFGYLETVNFCLLEVLIIVFLTVLIISGNIIVIFVFHCAPLNHHTTSYFIQTMAYADLFVGVSCVVPSLSLLHHPLPVEESLTCQIFGFVVSVLKSVSMASLACISIDRYIAITKPLYNTLVTPWRLRLCIFLIWLYSTLVFLPSFFHWG---K-P-VFQWCAE-----SWHTDSYFTLFIVMMLYAPAALIVCFTYFNIFRICQQHRFSGETGEVQACPDKRYAMVLFRITSVFYILWLPYIIYFLLESSTGHS--NRFASFLTTWLAISNSFCNCVIYSLSNSVFQRGLKRLSGAMCTSCASQTTANDGPLNGCHI----
>CKR3_MOUSE
TDEIKTVVESFETTPYEYEWAPPC---EKVRIKELGSWLLPPLYSLVFIIGLLGNMMVVLILIKYRKLQIMTNIYLFNLAISDLLFLFTVPFWIHYVLWNE-WGFGHYMCKMLSGFYYLALYSEIFFIILLTIDRYLAIVHAVALRARTVTFATITSIITWGLAGLAALPEFIFHESQD-SF-EFSCSPRYPEGEEDSWKRFHALRMNIFGLALPLLVMVICYSGIIKTL---------LRCPNKKKHKAIRLIFVVMIVFFIFWTPYNLVLLFSAFHRTFKHLDLAMQVTEVIAYTHCCVNPVIYAFVGERFRKHLRLFFHRNVAVYLGKYIPFLT-------SSVS
>HH1R_MOUSE
----------MRLPNTSSASEDKMCEGNRTAMASPQLLPLVVVLSSISLVTVGLNLGVLYAVRSERKLHTVGNLYIVSLSVADLIVGAIVMPMNILYLIMTKWSLGRPLCLFWLSMDYVASTASIFSVFILCIDRYRSVQQPLYLRYRTKTRASATILGAWFLSFLWVIPILGWHHFT--EL-EDKCETD------FYNVTWFKIMTAIINFYLPTLLMLWFYVKIYNGVRRHLRSQYVSGLHLNRERKAAKQLGCIMAAFILCWIPYFIFFMVIAFCNSC-CSEPVHMFTIWLGYINSTLNPLIYPLCNENFKKTFKKILHIRS-----------------------
>SSR1_HUMAN
GEGGGSRGPGAGAADGMEEPGRNASQNGTLSEGQGSAILISFIYSVVCLVGLCGNSMVIYVILRYAKMKTATNIYILNLAIADELLMLSVPFLVTSTLLRH-WPFGALLCRLVLSVDAVNMFTSIYCLTVLSVDRYVAVVHPIAARYRRPTVAKVVNLGVWVLSLLVILPIVVFSRTAADG--TVACNML-MPEPAQRWLVGFVLYTFLMGFLLPVGAICLCYVLIIAKMRMVALK-AGWQQRKRSERKITLMVMMVVMVFVICWMPFYVVQLVNVFAEQ--DDATVSQLSVILGYANSCANPILYGFLSDNFKRSFQRILCLS-----WMDNAAETALKSRAYSVED
>MC5R_MOUSE
-MNSSSTLTVLNLTLNASEDGILGSNVKNKSLACEEMGIAVEVFLTLGLVSLLENILVIGAIVKNKNLHSPMYFFVGSLAVADMLVSMSNAWETVTIYLLNNDTFVRHIDNVFDSMICISVVASMCSLLAIAVDRYITIFYALYHHIMTARRSGVIIACIWTFCISCGIVFIIYYY----------------------EESKYVIICLISMFFTMLFFMVSLYIHMFLLARNH-RIPRYNSVRQRTSMKGAITLTMLLGIFIVCWSPFFLHLILMISCPQNSCFMSYFNMYLILIMCNSVIDPLIYALRSQEMRRTFKEIVCCHGFRRPCRLLGGY------------
>OPSD_CHELA
MNGTEGPYFYIPMVNTTGIVRSPYEYPQYYLVNPAAYAALGAYMFLLILVGFPVNFLTLYVTLEHKKLRTPLNYILLNLAVADLFMVLGGFTTTMYTSMHGYFVLGRLGCNVEGFFATLGGEIALWSLVVLAIERWVVVCKPISNFRFSEDHAIMGLAFTWVMASACAVPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFVCHFLIPLVVVFFCYGRLLCAVKEAAAAQQESETTQRAEREVSRMVVIMVVAFLVCWCPYAGVAWYIFTHQGSEFGPLFMTFPAFFAKSSSIYNPMIYICMNKQFRHCMITTLCCGKNPFEEEEGASTSVSSSSVSPAA-
>OLF6_RAT
----------MAWSTGQNLSTPGPFILLGFPGPRSMRIGLFLLFLVMYLLTVVGNLAIISLVGAHRCLQTPMYFFLCNLSFLEIWFTTACVPKTLATFAPRGGVISLAGCATQMYFVFSLGCTEYFLLAVMAYDRYLAICLPLYGGIMTPGLAMRLALGSWLCGFSAITVPATLIARL---SNHFFCDISVLSCTDTQVVELVSFGIAFCVILGSCGITLVSYAYIITTI--------IKIPSARGRHRAFSTCSSHLTVVLIWYGSTIFLHVRTSVES---SLDLTKAITVLNTIVTPVLNPFIYTLRNKDVKEALRRTVKGK------------------------
>D4DR_HUMAN
---MGNRSTADADGLLAGRGPAAGASAGASAGLAGQGAAALVGGVLLIGAVLAGNSLVCVSVATERALQTPTNSFIVSLAAADLLLALLVLPLFVYSEVQGGWLLSPRLCDALMAMDVMLCTASIFNLCAISVDRFVAVAVPLYNRQGGSRRQLLLIGATWLLSAAVAAPVLCGLN---G-R--AVCRL---------EDRDYVVYSSVCSFFLPCPLMLLLYWATFRGLQRWPPPRRRRAKITGRERKAMRVLPVVVGAFLLCWTPFFVVHITQALCPACSVPPRLVSAVTWLGYVNSALNPVIYTVFNAEFRNVFRKALRACC-----------------------
>NY2R_MOUSE
VEVKVEPYGPGHTTPRGELPPDPEPELIDSTKLVEVQVILILAYCSIILLGVVGNSLVIHVVIKFKSMRTVTNFFIANLAVADLLVNTLCLPFTLTYTLMGEWKMGPVLCHLVPYAQGLAVQVSTITLTVIALDRHRCIVYHL-ESKISKRISFLIIGLAWGISALLASPLAIFREYSLFE--IVACTEKWPGEEKSVYGTVYSLSTLLILYVLPLGIISFSYTRIWSKLRNHSPG-AASDHYHQRRHKMTKMLVCVVVVFAVSWLPLHAFQLAVDIDSHVKEYKLIFTVFHIIAMCSTFANPLLYGWMNSNYRKAFLSAFRCEQRLDAIHSE---AKKNLEVKKNNG
>ETBR_MOUSE
SSAPAEVTKGGRGAGVPPRS-FPPPCQRNIEISKTFKYINTIVSCLVFVLGIIGNSTLLRIIYKNKCMRNGPNILIASLALGDLLHIIIDIPINTYKLLAEDWPFGAEMCKLVPFIQKASVGITVLSLCALSIDRYRAVASWSIKGIGVPKWTAVEIVLIWVVSVVLAVPEAIGFDMITRV--LRVCMLNQKTAFMQFYKTAKDWWLFSFYFCLPLAITAVFYTLMTCEMLRKSGM-IALNDHLKQRREVAKTVFCLVLVFALCWLPLHLSRILKLTLYDQSFLLVLDYIGINMASLNSCINPIALYLVSKRFKNCFKSCLCCWCQTFE-EKQSLEFKANDHGYDNFR
>NK2R_RAT
----MGTRAIVSDANILSGLESNATGVTAFSMPGWQLALWATAYLALVLVAVTGNATVIWIILAHERMRTVTNYFIINLALADLCMAAFNATFNFIYASHNIWYFGRAFCYFQNLFPITAMFVSIYSMTAIAADRYMAIVHPF-QPRLSAPSTKAIIAGIWLVALALASPQCFYST---GA---TKCVVAWPNDNGGKMLLLYHLVVFVLIYFLPLLVMFGAYSVIGLTLWKRPRHHGANLRHLQAKKKFVKAMVLVVLTFAICWLPYHLYFILGTFQEDIKFIQQVYLALFWLAMSSTMYNPIIYCCLNHRFRSGFRLAFRCCPWVTPTEE----HTPSLSRRVNRC
>OPR._MOUSE
GSHFQGNLSLLN---ETVPHHLLLNASHSAFLPLGLKVTIVGLYLAVCIGGLLGNCLVMYVILRHTKMKTATNIYIFNLALADTLVLLTLPFQGTDILLGF-WPFGNALCKTVIAIDYYNMFTSTFTLTAMSVDRYVAICHPIALDVRTSSKAQAVNVAIWALASVVGVPVAIMGSAQVEE---IECLVE-IPAPQDYWGPVFAICIFLFSFIIPVLIISVCYSLMIRRLRGVRLL-SGSREKDRNLRRITRLVLVVVAVFVGCWTPVQVFVLVQGLGVQPETAVAILRFCTALGYVNSCLNPILYAFLDENFKACFRKFCCASALHREMQVSDRVGLGCKTSETVPR
>YLD1_CAEEL
--------------------------MIIFYLYVATQVFVAIAFVLLMATAIIGNSVVMWIIYQHKVMHYGFNYFLFNMAFADLLIALFNVGTSWTYNLYYDWWYG-DLCTLTSFFGIAPTTVSVCSMMALSWDRCQAVVNPLQKRPLSRKRSVIAILIIWVVSTVTALPFAIAASVNSVTSKAHVCSAP--------VNTFFEKVLFGIQYALPIIILGSTFTRIAVAFRATEATSSLKNNHTRAKSKAVKMLFLMVVAFVVCWLPYHIYHAFALEEFFDARGKYAYLLIYWIAMSSCAYNPIIYCFANERFRIGFRYVFRWIPVIDCKKEQYEYMRSMAISLQKGR
>5HT_BOMMO
LPLQNCSWNSTGWEPNWNVTVWQASAPFDTPAALVRAAAKAVVLGLLILATVVGNVFVIAAILLERHLRSAANNLILSLAVADLLVACLVMPLGAVYEVVQRWTLGPELCDMWTSGDVLCCTASILHLVAIALDRYWAVTNI-YIHASTAKRVGMMIACVWTVSFFVCIAQLLGWKDPDS-E--LRCVVS--------QDVGYQIFATASSFYVPVLIILILYWRIYQTARKRPSLKPKEAADSKRERKAAKTLAIITGAFVACWLPFFVLAILVPTCDC-EVSPVLTSLSLWLGYFNSTLNPVIYTVFSPEFRHAFQRLLCGRRVRRRRAPQ---------------
>OPSD_SALPV
MNGTEGPYFYIPMVNTTGIVRSPYEYPQYYLVNPAAYAALGAYMFFLILLGFPINFLTLYVTLEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVLGRLGCNLEGFFATLGGEIGLWSLVVLAIERWVVVCKPISNFRFGENHAIMGLAFTWIMACACAVPPLVGWSRYIPEGMQCSCGVDYYTRAEGFNNESFVVYMFTCHFCIPLTIIGFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVILMVVGFLVCWLPYASVAWYIFSNQGSQFGPLFMTIPAFFAKSSSVYNPMIYICMNKQFRHCMITTLCCGKNPFEEEEGASTSSVSSSSVSPAA
>OPSD_LITMO
MNGTEGPYFYVPMVNTSGIVRSPYEYPQYYLVNPAAYAALGAYMFLLILVGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVLGRLGCNIEGFFATLGGEIALWSLVVLAIERWVVVCKPISNFRFGENHAIMGLAFTWLMAMACAAPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFVCHFLIPLMVVFFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVVIMVIAFLICWCPYAGVAWWIFTHQGSDFGPVFMTIPAFFAKSSSIYNPMIYICLNKQFRHCMITTLCCGKNPFEEEEGASTSVSSSSVSPAA-
>P2UR_RAT
----MAAGLDSWNSTINGTWEGDELGYKCRFNEDFKYVLLPVSYGVVCVLGLCLNVVALYIFLCRLKTWNASTTYMFHLAVSDSLYAASLPLLVYYYAQGDHWPFSTVLCKLVRFLFYTNLYCSILFLTCISVHRCLGVLRPLSLSWGHARYARRVAAVVWVLVLACQAPVLYFVTTS-----RITCHDT-SARELFSHFVAYSSVMLGLLFAVPFSIILVCYVLMARRLLKP---AYGTTGLPRAKRKSVRTIALVLAVFALCFLPFHVTRTLYYSFRSLNAINMAYKITRPLASANSCLDPVLYFLAGQRLVRFARDAKPATEPTPSPQARRKLTDTVRKDLSISS
>OLF4_CHICK
-------------MASGNCTTPTTFILSGLTDNPGLQMPLFMVFLAIYTITLLTNLGLIALISVDLHLQTPMYIFLQNLSFTDAAYSTVITPKMLATFLEERKTISYIGCILQYFSFVLLTVTESLLLAVMAYDRYVAICKPLYPSIMTKAVCWRLVKGLYSLAFLNSLVHTSGLLKL---SNHFFCDNSQISSSSTTLNELLVFIFGSLFAMSSIITILISYVFIILTV--------VRIRSKDGKYKAFSTCTSHLMAVSLFHGTVIFMYLRPVKLF---SLDTDKIASLFYTVVIPMLNPLIYSWRNKEVKDALRRVIATNVWIH--------------------
>LSHR_RAT
ENELSGWDYDYGFCSPKTLQCAPEPDAFNPCEDIMGYAFLRVLIWLINILAIFGNLTVLFVLLTSRYKLTVPRFLMCNLSFADFCMGLYLLLIASVDSQTKGWQTG-SGCGAAGFFTVFASELSVYTLTVITLERWHTITYAVLDQKLRLRHAIPIMLGGWLFSTLIATMPLVGISNY----KVSICLPM-----VESTLSQVYILSILILNVVAFVVICACYIRIYFAVQNP------ELTAPNKDTKIAKKMAILIFTDFTCMAPISFFAISAAFKVPLITVTNSKILLVLFYPVNSCANPFLYAIFTKAFQRDFLLLLSRFGCCKRRAELYRRSNCKNGFPGASK
>MSHR_BOVIN
PALGSQRRLLGSLNCTPPATLPFTLAPNRTGPQCLEVSIPDGLFLSLGLVSLVENVLVVAAIAKNRNLHSPMYYFICCLAVSDLLVSVSNVLETAVMLLLEAAAVVQQLDNVIDVLICGSMVSSLCFLGAIAVDRYISIFYALYHSVVTLPRAWRIIAAIWVASILTSLLFITYYY----------------------NNHKVILLCLVGLFIAMLALMAVLYVHMLARACQHIARKRQRPIHQGFGLKGAATLTILLGVFFLCWGPFFLHLSLIVLCPQHGCIFKNFNLFLALIICNAIVDPLIYAFRSQELRKTLQEVLQCSW-----------------------
>PE24_RAT
--------------------MSIPGVNASFSSTPERLNSPVTIPAVMFIFGVVGNLVAIVVLCKSRKKETTFYTLVCGLAVTDLLGTLLVSPVTIATYMKGQWPGDQALCDYSTFILLFFGLSGLSIICAMSIERYLAINHAYYSHYVDKRLAGLTLFAVYASNVLFCALPNMGLGRSERQYPGTWCFID---WTTNVTAYAAFSYMYAGFSSFLILATVLCNVLVCGALLRMAAAVASFRRIAGAEIQMVILLIATSLVVLICSIPLVVRVFINQLYQPSDISRNPDLQAIRIASVNPILDPWIYILLRKTVLSKAIEKIKCLFCRIGGSGRDGSRRTSSAMSGHSR
>5H1B_FUGRU
----MEGTNNTTGWTHFDSTSNRTSKSFDEEVKLSYQVVTSFLLGALILCSIFGNACVVAAIALERSLQNVANYLIGSLAVTDLMVSVLVLPMAALYQVLNRWTLGQIPCDIFISLDMLCCTSSILHLCVIALDRYWAITEPIYMKKRTPRRAAVLISVTWLVGFSISIPPMLIMRSQPA-N--KQCKIT--------QDPWYTIYSTFGAFYIPLTLMLVLYGRIFKAARFRRHEETKRKIALARERKTVKTLGIIMGTFILCWLPFFIVALVMPFCQESFMPHWLKDVINWLGYSNSLLNPIIYAYFNKDFQSAFKKIIKCHFCRA--------------------
>GP39_HUMAN
-------MASPSLPGSDCSQIIDHSHVPEFEVATWIKITLILVYLIIFVMGLLGNSATIRVTQVLQKLQKEVTDHMVSLACSDILVFLIGMPMEFYSIIWNPTSSYTLSCKLHTFLFEACSYATLLHVLTLSFERYIAICHPFYKAVSGPCQVKLLIGFVWVTSALVALPLLFAMGTEYETSNMSICTNL-----SRWTVFQSSIFGAFVVYLVVLLSVAFMCWNMMQVLMKSRPPKSESEESRTARRQTIIFLRLIVVTLAVCWMPNQIRRIMAAAKPKHRAYMILLPFSETFFYLSSVINPLLYTVSSQQFRRVFVQVLCCRLSLQHANHEKRLTDSARFVQRPLL
>BRB2_HUMAN
FSADMLNVTLQGPTLNGTFAQSKC---PQVEWLGWLNTIQPPFLWVLFVLATLENIFVLSVFCLHKSSCTVAEIYLGNLAAADLILACGLPFWAITISNNFDWLFGETLCRVVNAIISMNLYSSICFLMLVSIDRYLALVKTMMGRMRGVRWAKLYSLVIWGCTLLLSSPMLVFRTMKE-HN--TACVIS---YPSLIWEVFTNMLLNVVGFLLPLSVITFCTMQIMQVLRNN---EMQKFKEIQTERRATVLVLVVLLLFIICWLPFQISTFLDTLHRLGRIIDVITQIASFMAYSNSCLNPLVYVIVGKRFRKKSWEVYQGVCQKGGCRSEPIQLRTS-------I
>RDC1_HUMAN
YAEPGNFSDISWPCNSSDCIVVDTVMCPNMPNKSVLLYTLSFIYIFIFVIGMIANSVVVWVNIQAKTTGYDTHCYILNLAIADLWVVLTIPVWVVSLVQHNQWPMGELTCKVTHLIFSINLFSGIFFLTCMSVDRYLSITYFTTPSSRKKMVRRVVCILVWLLAFCVSLPDTYYLKTVTNNE--TYCRSFYPEHSIKEWLIGMELVSVVLGFAVPFSIIAVFYFLLARAI---------SASSDQEKHSSRKIIFSYVVVFLVCWLPYHVAVLLDIFSILHHALFTALHVTQCLSLVHCCVNPVLYSFINRNYRYELMKAFIFKYSAKTGLTKLIDEYSALEQNAK--
>OPSB_ASTFA
QEFQEDFYIPIPLDTNNITALSPFLVPQDHLGGSGIFMIMTVFMLFLFIGGTSINVLTIVCTVQYKKLRSHLNYILVNLAISNLLVSTVGSFTAFVSFLNRYFIFGPTACKIEGFVATLGGMVSLWSLSVVAFERWLVICKPVGNFSFKGTHAIIGCALTWFFALLASTPPLFGWSRYIPEGLQCSCGPDWYTTENKYNNESYVMFLFCFCFGFPFTVILFCYGQLLFTLKSAAKAQADSASTQKAEREVTKMVVVMVMGFLVCWLPYASFALWVVFNRGQSFDLRLGTIPSCFSKASTVYNPVIYVFMNKQFRSCMMKLIFCGKSPFGDDEEASSSSVGPEK-----
>CB1B_FUGRU
TNASDFPLSNGSGEATQCGEDIVDNMECFMILTPAQQLVIVILAITLGTFTVLENFVVLCVILHSHTRSRPSYHFIGSLAVADLIGSIIFVYSFLDFHVLHR-KDSPSIFLFKLAGVIASFTASVGSLFLTAIDRYVSIHRPMYKRIITKTKAVIAFSVMWAISIEFSLLPLLGWNCK-R--LHSVCSDI------FPLIDEKYLMFWIGMTTVLLLFIIYAYMFILWKSHHHSEGQTVRPEQARMDLRLAKTLVLILVALIICWGPLLAIMVYDLFGRVNDFIKTVFAFCSMLCLLNSTINPVIYAMRSKDLRRAFVNICHMCRGTTQSLDSSAEVRSTGGRAGKDR
>5H2B_HUMAN
ILQSTFVHVISSNWSGLQTESIPEEMKQIVEEQGNKLHWAALLILMVIIPTIGGNTLVILAVSLEKKLQYATNYFLMSLAVADLLVGLFVMPIALLTIMFEAWPLPLVLCPAWLFLDVLFSTASIMHLCAISVDRYIAIKKPIANQYNSRATAFIKITVVWLISIGIAIPVPIKGIET-NP---ITCVLT------KERFGDFMLFGSLAAFFTPLAIMIVTYFLTIHALQKKRRTGKKSVQTISNEQRASKVLGIVFFLFLLMWCPFFITNITLVLCDSCTTLQMLLEIFVWIGYVSSGVNPLVYTLFNKTFRDAFGRYITCNYRATKSVKTLRKRNPMAENSKFFK
>OAR2_LYMST
RNFSVSADVWLCGANFSQEWQLMQPVCSTKYDSITIFITVAVVLTLITLWTILGNFFVLMALYRYGTLRTMSNCLIGNLAISDLLLAVTVLPISTVHDLLGYWVFGEFTCTLWLCMDVLYCTASIWGLCTVAFDRYLATVYPVYHDQRSVRKAVGCIVFVWIFSIVISFAPFIGWQHM-S-F--YQCILF--------TSSSYVLYSSMGSFVIPAILMAFMYVRIFVVLHNQRNKLSMKRRFELREQRATKRMLLIMACFCVCWMPFLFMYILRSVCDTCHMNQHFVAAIIWLGYVNSSLNPVLYTLFNDDFKVAFKRLIGARSPSAYRSPGPRR------------
>GRPR_HUMAN
DCFLLNLEVDHFMHCN-ISSHSADLPVNDDWSHPGILYVIPAVYGVIILIGLIGNITLIKIFCTVKSMRNVPNLFISSLALGDLLLLITCAPVDASRYLADRWLFGRIGCKLIPFIQLTSVGVSVFTLTALSADRYKAIVRPMIQASHALMKICLKAAFIWIISMLLAIPEAVFSDLHPQT--FISCAPY---HSNELHPKIHSMASFLVFYVIPLSIISVYYYFIAKNLIQSLPVNIHVKKQIESRKRLAKTVLVFVGLFAFCWLPNHVIYLYRSYHYSEMLHFVTSICARLLAFTNSCVNPFALYLLSKSFRKQFNTQLLCCQPGLIIR--SHSMTSLKSTNPSVA
>DBDR_HUMAN
PGSNGTAYPGQFALYQQLAQGNAVGGSAGAPPLGPSQVVTACLLTLLIIWTLLGNVLVCAAIVRSRHRANMTNVFIVSLAVSDLFVALLVMPWKAVAEVAGYWPFG-AFCDVWVAFDIMCSTASILNLCVISVDRYWAISRPFYKRKMTQRMALVMVGLAWTLSILISFIPVQLNWHRDED-NAENCDSS--------LNRTYAISSSLISFYIPVAIMIVTYTRIYRIAQVQSAADTSLRASIKKETKVLKTLSVIMGVFVCCWLPFFILNCMVPFCSGHCVSETTFDVFVWFGWANSSLNPVIYAFNA-DFQKVFAQLLGCSHFCSRTPVETVNYNQDIVFHK---
>B2AR_MESAU
MGPPGNDSDFLLTTNGSHVPDH----DVTEERDEAWVVGMAILMSVIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGASHILMKMWNFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYIAITSPFYQSLLTKNKARMVILMVWIVSGLTSFLPIQMHWYRAI-D--TCCDFF--------TNQAYAIASSIVSFYVPLVVMVFVYSRVFQVAKRQGRSLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNL-IPKEVYILLNWLGYVNSAFNPLIYCRSP-DFRIAFQELLCLRRSSSKAYGNGYSDYMGEASGCQLG
>5H1B_HUMAN
PAGSETWVPQANLSSAPSQNCSAKDYIYQDSISLPWKVLLVMLLALITLATTLSNAFVIATVYRTRKLHTPANYLIASLAVTDLLVSILVMPISTMYTVTGRWTLGQVVCDFWLSSDITCCTASILHLCVIALDRYWAITDAVYSAKRTPKRAAVMIALVWVFSISISLPPFFWR----E-E--SECVVN-------TDHILYTVYSTVGAFYFPTLLLIALYGRIYVEARSRRVSLEKKKLMAARERKATKTLGIILGAFIVCWLPFFIISLVMPICKDAWFHLAIFDFFTWLGYLNSLINPIIYTMSNEDFKQAFHKLIRFKCTS---------------------
>GPRH_HUMAN
------MNGLEVAPPGLITNFSLATAEQCGQETPLENMLFASFYLLDFILALVGNTLALWLFIRDHKSGTPANVFLMHLAVADLSCVLVLPTRLVYHFSGNHWPFGEIACRLTGFLFYLNMYASIYFLTCISADRFLAIVHPVSLKLRRPLYAHLACAFLWVVVAVAMAPLLVSPQTVQVCL-QLYRE----------KASHHALVSLAVAFTFPFITTVTCYLLIIRSLR------QGLRVEKRLKTKAVRMIAIVLAIFLVCFVPYHVNRSVYVLHYRSRILALANRITSCLTSLNGALDPIMYFFVAEKFRHALCNLLCGKRLKGPPPSFEGKAKSEL-------
>ACM1_MOUSE
------------MNTSVPPAVSPNITVLAPGKGPWQVAFIGSTTGLLSLATVTGNLLVLISIKVNTELKTVNNYFLLSLACADLIIGTFSMNLYTTYLLMGHWALGTLACDLWLALDYVASNASVMNLLLISFDRYFSVTRPLYRAKRTPRRAALMIGLAWLVSFVLWAPAILFWQYLV-VL-AGQCYIQ------FLSQPIITFGTAMAAFYLPVTVMCTLYWRIYRETENRRGKAKRKTFSLVKEKKAARTLSAILLAFILTWTPYNIMVLVSTFCKDC-VPETLWELGYWLCYVNSTVNPMCYASCNKAFRDHFRLLLLCRWDKRRWRKIPKRPSRQC-------
>B3AR_SHEEP
MAPWPPGNSFLTPWPDIPTLAPNTANASGLPGVPWAVALAGALLALAVLATVGGNLLVIVAIARTPRLQTMTNVFVTSLATADLVVGLLVVPPGATLALTGHWPLGVTGCELWTSVDVLCVTASIETLCALAVDRYLAVTNPLYGALVTKRRARAAVVLVWVVSAAVSFAPIMSKWWRVQ-R--RCCTFA--------SNMPYALLSSSVSFYLPLLVMLFVYARVFVVDTRQGVPRRPARLLPLREHRALRTLGLIMGTFTLCWLPFFVVNVVRALGGPS-VSGPTFLALNWLGYANSAFNPLIYCRSP-DFRSAFRRLLCRCPPEEHLAAASPPTVLTSPAGPRQP
>A1AD_MOUSE
TGSGEDNQSSTAEAGAAASGEVNGSAAVGGLVVSAQGVGVGVFLAAFILTAVAGNLLVILSVACNRHLQTVTNYFIVNLAVADLLLSAAVLPFSATMEVLGFWPFGRTFCDVWAAVDVLCCTASILSLCTISVDRYVGVRHSLYPAIMTERKAAAILALLWAVALVVSVGPLLGWKEP--VP--RFCGIT--------EEVGYAIFSSVCSFYLPMAVIVVMYCRVYVVARSTHTLLSVRLLKFSREKKAAKTLAIVVGVFVLCWFPFFFVLPLGSLFPQLKPSEGVFKVIFWLGYFNSCVNPLIYPCSSREFKRAFLRLLRCQCRRRRRR--LWPSLDRR-PALRLC
>C5AR_CANFA
SMNFSPPEYPDYG-TATLDPNIFVDESLNTPKLSVPDMIALVIFVMVFLVGVPGNFLVVWVTGFEVR-RTINAIWFLNLAVADLLSCLALPILFSSIVQQGYWPFGNAACRILPSLILLNMYASILLLTTISADRFVLVFNPICQNYRGPQLAWAACSVAWAVALLLTVPSFIFRGVHT-PF--MTCGVD-YSGVGVLVERGVAILRLLMGFLGPLVILSICYTFLLIRT---------WSRKATRSTKTLKVVVAVVVSFFVLWLPYQVTGMMMALFYKHRRVSRLDSLCVAVAYINCCINPIIYVLAAQGFHSRFLKSLPARLRQVLAEESVGRSTVD--------
>CKRA_HUMAN
QVSWGHYSGDEEDAYSAEPLPELC---YKADVQAFSRAFQPSVSLTVAALGLAGNGLVLATHLAARRARSPTSAHLLQLALADLLLALTLPFAAAGALQG--WSLGSATCRTISGLYSASFHAGFLFLACISADRYVAIARALGPRPSTPGRAHLVSVIVWLLSLLLALPALLFS--QD-RE-QRRCRLIFPEGLTQTVKGASAVAQVALGFALPLGVMVACYALLGRTL---------LAARGPERRRALRVVVALVAAFVVLQLPYSLALLLDTADLLAKRKDVALLVTSGLALARCGLNPVLYAFLGLRFRQDLRRLLRGGSSPSGPQPRRGC--------PRLS
>MTR_BUFMA
LNLDCSELPNSSWVNSSMENQSSNSTRDPLKRNEEVAKVEVTVLALILFLALAGNICVLLGIYINRHKHSRMYFFMKHLSIADLVVAIFQVLPQLIWDITFRFYAPDLVCRLVTYLQVVGMFASTYMLLLMSLDRCLAICQPL--RSLHRRSDCVYVLFTWILSFLLSTPQTVIFSLTE----VYDCRAD---FIQPWGPKAYITWITLAVYIIPVMILSVCYGLISYKIWQNRATVSSVRLISKAKIRTVKMTFIIVLAYIVCWTPFFFVQMWSVWDPNP-KEASLFIIAMLLGSLNSCCNPWIYMLFTGHLFHDLLQSFLCCSARYLKTQQQGSSNSSTFVLSRKS
>B1AR_HUMAN
LPDGAATAARLLVPASPPASLLPPASESPEPLSQQWTAGMGLLMALIVLLIVAGNVLVIVAIAKTPRLQTLTNLFIMSLASADLVMGLLVVPFGATIVVWGRWEYGSFFCELWTSVDVLCVTASIETLCVIALDRYLAITSPFYQSLLTRARARGLVCTVWAISALVSFLPILMHWWRAR-R--KCCDFV--------TNRAYAIASSVVSFYVPLCIMAFVYLRVFREAQKQNGRRRPSRLVALREQKALKTLGIIMGVFTLCWLPFFLANVVKAFHREL-VPDRLFVFFNWLGYANSAFNPIIYCRSP-DFRKAFQGLLCCARRAARRRHATHGCLARPGPPPSPG
>ACM2_RAT
--------------MNNSTNSSNNGLAITSPYKTFEVVFIVLVAGSLSLVTIIGNILVMVSIKVSRHLQTVNNYFLFSLACADLIIGVFSMNLYTLYTVIGYWPLGPVVCDLWLALDYVVSNASVMNLLIISFDRYFCVTKPLYPVKRTTKMAGMMIAAAWVLSFILWAPAILFWQFIV-VE-DGECYIQ------FFSNAAVTFGTAIAAFYLPVIIMTVLYWHISRASKSRVKMPAKKKPPPSREKKVTRTILAILLAFIITWAPYNVMVLINTFCAPC-IPNTVWTIGYWLCYINSTINPACYALCNATFKKTFKHLLMCHYKNIGATR----------------
>AG22_MERUN
AATSRNITSSLPFVNLNMSGTNDLIFNCSHKPSDKHLEAIPVLYYLIFVIGFAVNIIVVSLFCCQKGPKKVSSIYIFNLAVADLLLLATLPLWATYYSYRYDWLFGPVMCKVFGSFLTLNMFASIFFITCMSVDRYQSVIYPFLSQRRNPWQASYVVPLVWCMACLSSLPTFYFRDVRT-LG--NACVMAFPPEKYAQWSAGIALMKNVLGFIIPLIFIATCYFGIRKHLLKT----NSYGKNRITRDQVLKMAAAVVLAFIICWLPFHVLTFLDALSWMGAVIDLALPFAILLGFTNSCVNPFLYCFVGNRFQQKLRSMFRVPITWLQGKRETMSREMDTFVS----
>NY5R_PIG
-------------NTVATRNSGFPVWEDYKGSVDDLQYFLIGLYTFVSLLGFMGNLLILMAVMRKRNQKTTVNFLIGNLAFSDILVVLFCSPFTLTSVLLDQWMFGKVMCHIMPFLQCVTVLVSTLILISIAIVRYHMIKHPV-SNNLTANHGYFLIATVWTLGLAICSPLPVFHSLVESS--RYLCVES---WPSDSYRIAFTISLLLVQYILPLVCLTVSHTSVCRTISCGVPESRSIMRLRKRSRSVFYRLTVLILVFAVSWMPLHLFHVVTDFNDNLRHFKLVYCICHLLGMMSCCLNPILYGFLNNGIKADLMSLIHCLHVS---------------------
>CB1R_HUMAN
EFYNKSLSSFKENEENIQCGENFMDIECFMVLNPSQQLAIAVLSLTLGTFTVLENLLVLCVILHSRSRCRPSYHFIGSLAVADLLGSVIFVYSFIDFHVFHR-KDSRNVFLFKLGGVTASFTASVGSLFLTAIDRYISIHRPLYKRIVTRPKAVVAFCLMWTIAIVIAVLPLLGWNCE-K--LQSVCSDI------FPHIDETYLMFWIGVTSVLLLFIVYAYMYILWKAHSHSEDQVTRPDQARMDIRLAKTLVLILVVLIICWGPLLAIMVYDVFGKMNKLIKTVFAFCSMLCLLNSTVNPIIYALRSKDLRHAFRSMFPSCEGTAQPLDNSMGKHANNAASVHRA
>D3DR_MOUSE
--------MAPLSQISSHINSTCGAENSTGVNRARPHAYYALSYCALILAIIFGNGLVCAAVLRERALQTTTNYLVVSLAVADLLVATLVMPWVVYLEVTGGWNFSRICCDVFVTLDVMMCTASILNLCAISIDRYTAVVMPVGTGQSSCRRVALMITAVWVLAFAVSCPLLFGFN---D-P--SICSI---------SNPDFVIYSSVVSFYVPFGVTVLVYARIYMVLRQRTSLPLQPRGVPLREKKATQMVVIVLGAFIVCWLPFFLTHVLNTHCQACHVSPELYRATTWLGYVNSALNPVIYTTFNIEFRKAFLKILSC-------------------------
>GPRF_MACMU
----MDPEETSVYLDYYYATSPNPDIRETHSHVPYTSVFLPVFYTAVFLTGVLGNLVLMGALHFKPGSRRLIDIFIINLAASDFIFLVTLPLWVDKEASLGLWRTGSFLCKGSSYMISVNMHCSVFLLTCMSVDRYLAIVCPVSRKFRRTDCAYVVCASIWFISCLLGLPTLLSRELT-IDD-KPYCAEK----KATPLKLIWSLVALIFTFFVPLLSIVTCYCCIARKLCAH---YQQSGKHNKKLKKSIKIIFIVVAAFLVSWLPFNTFKLLAIVSGLQAMLQLGMEVSGPLAFANSCVNPFIYYIFDSYIRRAIVHCLCPCLKNYDFGSSTETALSTFIHAEDFT
>CB1R_RAT
EFYNKSLSSFKENEENIQCGENFMDMECFMILNPSQQLAIAVLSLTLGTFTVLENLLVLCVILHSRSRCRPSYHFIGSLAVADLLGSVIFVYSFVDFHVFHR-KDSPNVFLFKLGGVTASFTASVGSLFLTAIDRYISIHRPLYKRIVTRPKAVVAFCLMWTIAIVIAVLPLLGWNCK-K--LQSVCSDI------FPLIDETYLMFWIGVTSVLLLFIVYAYMYILWKAHSHSEDQVTRPDQARMDIRLAKTLVLILVVLIICWGPLLAIMVYDVFGKMNKLIKTVFAFCSMLCLLNSTVNPIIYALRSKDLRHAFRSMFPSCEGTAQPLDNSMGHANNTASMHRAA
>CKR4_HUMAN
IADTTLDESIYSNYYLYESIPKPC---TKEGIKAFGELFLPPLYSLVFVFGLLGNSVVVLVLFKYKRLRSMTDVYLLNLAISDLLFVFSLPFWGYYAADQ--WVFGLGLCKMISWMYLVGFYSGIFFVMLMSIDRYLAIVHAVSLRARTLTYGVITSLATWSVAVFASLPGFLFSTCYT-ER-HTYCKTK-YSLNSTTWKVLSSLEINILGLVIPLGIMLFCYSMIIRTL---------QHCKNEKKNKAVKMIFAVVVLFLGFWTPYNIVLFLETLVELERYLDYAIQATETLAFVHCCLNPIIYFFLGEKFRKYILQLFKTCRGLFVLCQYCGLTP------SSSY
>PD2R_MOUSE
----------------------MNESYRCQTSTWVERGSSATMGAVLFGAGLLGNLLALVLLARSGLPPSVFYVLVCGLTVTDLLGKCLISPMVLAAYAQNQPASGNQLCETFAFLMSFFGLASTLQLLAMAVECWLSLGHPFYQRHVTLRRGVLVAPVVAAFCLAFCALPFAGFGKFVQYCPGTWCFIQ-MIHKERSFSVIGFSVLYSSLMALLVLATVVCNLGAMYNLYDMAQSYRHGSLHPLEELDHFVLLALMTVLFTMCSLPLIYRAYYGAFKL--AEGDSEDLQALRFLSVISIVDPWIFIIFRTSVFRMLFHKVFTRPLIYRNWSSHSQL-----------
>OPSD_BUFMA
MNGTEGPNFYIPMSNKTGVVRSPFEYPQYYLAEPWQYSVLCAYMFLLILLGFPINFMTLYVTIQHKKLRTPLNYILLNLAFANHFMVLCGFTVTMYSSMNGYFVFGQTGCYVEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFSENHAIMGVAFTWIMALACAAPPLFGWSRYIPEGMQCSCGVDYYTLKPEVNNESFVIYMFVVHFLIPLIIIFFCYGRLVCTVKEAAAQQQESATTQKAEKEVTRMVIIMVVFFLICWVPYASVAFFIFTHQGSEFGPVFMTIPAFFAKSSSIYNPVIYIMLNKQFRNCMITTLCCGKNPFGDEDASSAASSVSSSQVSPA
>5H1F_MOUSE
-------------MDFLNASD-QNLTSEELLNRMPSKILVSLTLSGLALMTTTINSLVIAAIIVTRKLHHPANYLICSLAVTDFLVAVLVMPFSIVYIVRESWIMGQVLCDIWLSVDIICCTCSILHLSAIALDRYRAITDAVYARKRTPRHAGIMITIVWVISVFISMPPLFWR----S-R--DECVIK-------HDHIVSTIYSTFGAFYIPLVLILILYYKIYRAARTLLKHWRRQKISGTRERKAATTLGLILGAFVICWLPFFVKELVVNVCEKCKISEEMSNFLAWLGYLNSLINPLIYTIFNEDFKKAFQKLVRCRY-----------------------
>B2AR_RAT
MEPHGNDSDFLLAPNGSRAPGH----DITQERDEAWVVGMAILMSVIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGASHILMKMWNFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYVAITSPFYQSLLTKNKARVVILMVWIVSGLTSFLPIQMHWYRAI-D--TCCDFF--------TNQAYAIASSIVSFYVPLVVMVFVYSRVFQVAKRQGRSLRSSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIRANL-IPKEVYILLNWLGYVNSAFNPLIYCRSP-DFRIAFQELLCLRRSSSKTYGNGYSDYTGEQSAYQLG
>O.YR_PIG
GVLAANWSAEAVNSSAAPPEAEGNRTAGPPQRNEALARVEVAVLCLILFLALSGNACVLLALRTTRHKHSRLFFFMKHLSIADLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPL--RALRRPADRLAVLATWLGCLVASAPQVHIFSLRE----VFDCWAV---FIQPWGPKAYITWITLAVYIVPVIVLAACYGLISFKIWQNRAAVSSVKLISKAKIRTVKMTFIIVLAFIVCWTPFFFVQMWSVWDADA-KEASAFIIAMLLASLNSCCNPWIYMLFTGHLFHELVQRFLCCSSSHLKTSRPGENSSTFVLSQHSS
>OPSB_CARAU
PEFHEDFYIPIPLDINNLSAYSPFLVPQDHLGNQGIFMAMSVFMFFIFIGGASINILTILCTIQFKKLRSHLNYILVNLSIANLFVAIFGSPLSFYSFFNRYFIFGATACKIEGFLATLGGMVGLWSLAVVAFERWLVICKPLGNFTFKTPHAIAGCILPWISALAASLPPLFGWSRYIPEGLQCSCGPDWYTTNNKYNNESYVMFLFCFCFAVPFGTIVFCYGQLLITLKLAAKAQADSASTQKAEREVTKMVVVMVLGFLVCWAPYASFSLWIVSHRGEEFDLRMATIPSCLSKASTVYNPVIYVLMNKQFRSCMMKMVCGKNIEE---DEASTSSVAPEK-----
>5H1B_RAT
CAPPPPATSQTGVPLANLSHNSADDYIYQDSIALPWKVLLVALLALITLATTLSNAFVIATVYRTRKLHTPANYLIASLAVTDLLVSILVMPISTMYTVTGRWTLGQVVCDFWLSSDITCCTASIMHLCVIALDRYWAITDAVYSAKRTPKRAAIMIVLVWVFSISISLPPFFWR----E-E--LDCFVN-------TDHVLYTVYSTVGAFYLPTLLLIALYGRIYVEARSRRVSLEKKKLMAARERKATKTLGIILGAFIVCWLPFFIISLVMPICKDAWFHMAIFDFFNWLGYLNSLINPIIYTMSNEDFKQAFHKLIRFKCTG---------------------
>H963_HUMAN
-----------------------MTNSSFFCPVYKDLEPFTYFFYLVFLVGIIGSCFATWAFIQKNTNHRCVSIYLINLLTADFLLTLALPVKIVVDLGVAPWKLKIFHCQVTACLIYINMYLSIIFLAFVSIDRCLQLTHSCIYRIQEPGFAKMISTVVWLMVLLIMVPNMMIPIKD-SN---VGCMEF--KKEFGRNWHLLTNFICVAIFLNFSAIILISNCLVIRQLYR-----NKDNENYPNVKKALINILLVTTGYIICFVPYHIVRIPYTLSQTEISLFKAKEATLLLAVSNLCFDPILYYHLSKAFRSKVTETFASPKETKAQKEKLRC------------
>OPS3_DROME
ARLSAETRLLGWNVPPEELRHIPEHWLTYPEPPESMNYLLGTLYIFFTLMSMLGNGLVIWVFSAAKSLRTPSNILVINLAFCDFMMMVKTPIFIYNSFHQG-YALGHLGCQIFGIIGSYTGIAAGATNAFIAYDRFNVITRPM-EGKMTHGKAIAMIIFIYMYATPWVVACYTETWGRFPEGYLTSCTFD--YLTDNFDTRLFVACIFFFSFVCPTTMITYYYSQIVGHVFSHNVESNVDKNKETAEIRIAKAAITICFLFFCSWTPYGVMSLIGAFGDKTLLTPGATMIPACACKMVACIDPFVYAISHPRYRMELQKRCPWLALNEKAPESSAVEPQQTTAA----
>CCR5_MOUSE
DDLYKELAFYSNSTEIPLQDSNFCSTVEGPLLTSFKAVFMPVAYSLIFLLGMMGNILVLVILERHRHTRSSTETFLFHLAVADLLLVFILPFAVAEGSVG--WVLGTFLCKTVIALHKINFYCSSLLVACIAVDRYLAIVHAVAYRRRRLLSIHITCTAIWLAGFLFALPELLFAKVGQ-ND-LPQCTFSQENEAETRAWFTSRFLYHIGGFLLPMLVMGWCYVGVVHRL--------LQAQRRPQRQKAVRVAILVTSIFFLCWSPYHIVIFLDTLERLKGYLSVAITLCEFLGLAHCCLNPMLYTFAGVKFRSDLSRLLTKLGCAG---PASLC--------PNWR
>5H1F_CAVPO
-------------MDFLNSSD-QNLTSEELLHRMPSKILVSLTLSGLALMTTTINSLVIAAIIVTRKLHHPANYLICSLAVTDFLVAVLVMPFSIVYIVRESWIMGQVLCDIWLSVDIICCTCSILHLSAIALDRYRAITDAVYARKRTPKQAGIMITIVWIISVFISMPPLFWR----S-R--DECIIK-------HDHIVSTIYSTFGAFYIPLVLILILYYKIYKAAKTLLRHWRRQKISGTRERKAATTLGLILGAFVICWLPFFVKELVVNVCEKCKISEEMANFLAWLGYLNSLINPLIYTIFNEDFKKAFQKLVRCQY-----------------------
>OPSD_MESBI
MNGTEGLNFYVPFSNHTGVVRSPFEYPQYYLAEPWQFSVLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVANLFMVLGGFTTTLYTSMHAYFIFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLALTWIMALACAAPPLVGWSRYIPEGMQCSCGVDYYTPSPEVNNESFVVYMFVVHFSIPMVIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVVIMVVAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPSFFAKSSAIYNPVIYIMMNKQFRNCMLTTLCCGRNPLGDDEVSTTASKTETSQVAPA
>PI2R_RAT
GRPDGPPSITPESPLIVGGREWQGMAGSCWNITYVQDSVGPATSTLMFVAGVVGNGLALGILGARRRHPSAFAVLVTGLAVTDLLGTCFLSPAVFVAYARNSAHGGTMLCDTFAFAMTFFGLASTLILFAMAVERCLALSHPYYAQLDGPRCARLALPAIYAFCCLFCSLPLLGLGEHQQYCPGSWCFIR---MRSPQPGGCAFSLAYASLMALLVTSIFFCNGSVTLSLCHMRRHFVPTSRAREDEVYHLILLALMTGIMAVCSLPLTIRGFTQAIAPDS--REMGDLHAFRFNAFNPILDPWVFILFRKAVFQRLKFWLCCLCARSVHGDLQTPRRDTLAPDSLQA
>IL8A_PANTR
PQMWDFDDLNFTGMPPTDEGYSPC----RLETETLNKYVVIITYALVFLLSLLGNSLVMLVILYSRVGRSVTDVYLLNLALADLLFALTLPIWAASKVNG--WIFGTFLCKVVSLLKEVNFYSGILLLACISVDRYLAIVHATRTLTQKRHLVKFVCLGCWGLSMNLSLPFFLFRQAYH-NS-SPVCYEV-LGNDTAKWRMVLRILPHTFGFIVPLFVMLFCYGFTLRTL---------FKAHMGQKHRAMRVIFAVVLIFLLCWLPYNLVLLADTLMRTQNNIGRALDATEILGFLHSCLNPIIYAFIGQNFRHGFLKILAMHGLVSKEFLARHR------------
>CCR5_HUMAN
EDLFWELDRLDNYNDTSLVENHLCPATEGPLMASFKAVFVPVAYSLIFLLGVIGNVLVLVILERHRQTRSSTETFLFHLAVADLLLVFILPFAVAEGSVG--WVLGTFLCKTVIALHKVNFYCSSLLLACIAVDRYLAIVHAVAYRHRRLLSIHITCGTIWLVGFLLALPEILFAKVSQ-NN-LPRCTFSQENQAETHAWFTSRFLYHVAGFLLPMLVMGWCYVGVVHRL--------RQAQRRPQRQKAVRVAILVTSIFFLCWSPYHIVIFLDTLARLKGSLPVAITMCEFLGLAHCCLNPMLYTFAGVKFRSDLSRLLTKLGCTG---PASLC--------PSWR
>CB2R_HUMAN
----MEECWVTEIANGSKDGLDSNPMKDYMILSGPQKTAVAVLCTLLGLLSALENVAVLYLILSSHQRRKPSYLFIGSLAGADFLASVVFACSFVNFHVFHG-VDSKAVFLLKIGSVTMTFTASVGSLLLTAIDRYLCLRYPPYKALLTRGRALVTLGIMWVLSALVSYLPLMGWTC-----CPRPCSEL------FPLIPNDYLLSWLLFIAFLFSGIIYTYGHVLWKAHQHSGHQVPGMARMRLDVRLAKTLGLVLAVLLICWFPVLALMAHSLATTLSDQVKKAFAFCSMLCLINSMVNPVIYALRSGEIRSSAHHCLAHWKKCVRGLGSEAKVTETEADGKITP
>AA2B_HUMAN
------------------------------MLLETQDALYVALELVIAALSVAGNVLVCAAVGTANTLQTPTNYFLVSLAAADVAVGLFAIPFAITISLG--FCTDFYGCLFLACFVLVLTQSSIFSLLAVAVDRYLAICVPLYKSLVTGTRARGVIAVLWVLAFGIGLTPFLGWNSKDT-T-LVKCLFE-----NVVPMSYMVYFNFFGCVLPPLLIMLVIYIKIFLVACRQ---MDHSRTTLQREIHAAKSLAMIVGIFALCWLPVHAVNCVTLFQPAQNKPKWAMNMAILLSHANSVVNPIVYAYRNRDFRYTFHKIISRYLLCQADVKSGNGLGVGL-------
>MSHR_HUMAN
AVQGSQRRLLGSLNSTPTAIPQLGLAANQTGARCLEVSISDGLFLSLGLVSLVENALVVATIAKNRNLHSPMYCFICCLALSDLLVSGTNVLETAVILLLEAAAVLQQLDNVIDVITCSSMLSSLCFLGAIAVDRYISIFYALYHSIVTLPRAPRAVAAIWVASVVFSTLFIAYYY----------------------DDHVAVLLCLVVFFLAMLVLMAVLYVHMLARACQHIARKRQRPVHQGFGLKGAVTLTILLGIFFLCWGPFFLHLTLIVLCPEHGCIFKNFNLFLALIICNAIIDPLIYAFHSQELRRTLKEVLTCSW-----------------------
>ETBR_PIG
SSSPPQMPKGGRMAGPPARTLTPPPCEGPIEIKDTFKYINTVVSCLVFVLGIIGNSTLLRIIYKNKCMRNGPNILIASLALGDLLHIIIDIPINVYKLLAEDWPFGVEMCKLVPFIQKASVGITVLSLCALSIDRYRAVASWSIKGIGVPKWTAVEIVLIWVVSVVLAVPEALGFDMITRI--LRICLLHQKTAFMQFYKTAKDWWLFSFYFCLPLAITAFFYTLMTCEMLRKSGM-IALNDHLKQRREVAKTVFCLVLVFALCWLPLHLSRILKLTLYDQSFLLVLDYIGINMASLNSCINPIALYLVSKRFKNCFKSCLCCWCQSFE-EKQSLEFKANDHGYDNFR
>OAR2_LOCMI
SSAAEEPQDALVGGDACGGRRPPSVLGVRLAVPEWEVAVTAVSLSLIILITIVGNVLVVLSVFTYKPLRIVQNFFIVSLAVADLTVAVLVMPFNVAYSLIQRWVFGIVVCKMWLTCDVLCCTASILNLCAIALDRYWAITDPIYAQKRTLRRVLAMIAGVWLLSGVISSPPLIGWNDW-N-D--TPCQLT--------EEQGYVIYSSLGSFFIPLFIMTIVYVEIFIATKRRPVYEEKQRISLSKERRAARTLGIIMGVFVVCWLPFFLMYVIVPFCNPSKPSPKLVNFITWLGYINSALNPIIYTIFNLDFRRAFKKLLHFKT-----------------------
>MC5R_HUMAN
-MNSSFHLHFLDLNLNATEGNLSGPNVKNKSSPCEDMGIAVEVFLTLGVISLLENILVIGAIVKNKNLHSPMYFFVCSLAVADMLVSMSSAWETITIYLLNNDAFVRHIDNVFDSMICISVVASMCSLLAIAVDRYVTIFYALYHHIMTARRSGAIIAGIWAFCTGCGIVFILYYS----------------------EESTYVILCLISMFFAMLFLLVSLYIHMFLLARTH-RIPGASSARQRTSMQGAVTVTMLLGVFTVCWAPFFLHLTLMLSCPQNSRFMSHFNMYLILIMCNSVMDPLIYAFRSQEMRKTFKEIICCRGFRIACSFPRRD------------
>5H6_RAT
-----------MVPEPGPVNSSTPAWGPGPPPAPGGSGWVAAALCVVIVLTAAANSLLIVLICTQPARNTS-NFFLVSLFTSDLMVGLVVMPPAMLNALYGRWVLARGLCLLWTAFDVMCCSASILNLCLISLDRYLLILSPLYKLRMTAPRALALILGAWSLAALASFLPLLLGWH---E-APGQCRLL--------ASLPFVLVASGVTFFLPSGAICFTYCRILLAARKQMESRRLATKHSRKALKASLTLGILLGMFFVTWLPFFVANIAQAVCDC--ISPGLFDVLTWLGYCNSTMNPIIYPLFMRDFKRALGRFLPCVHCPPEHRPALPPAVPDQASACSRC
>5H1D_CAVPO
SPPNQSEEGLPQEASNRSLNATETPGDWDPGLLQALKVSLVVVLSIITLATVLSNAFVLTTILLTRKLHTPANYLIGSLATTDLLVSILVMPISIAYTTTRTWNFGQILCDIWVSSDITCCTASILHLCVIALDRYWAITDALYSKRRTAGHAGAMIAAVWVISICISIPPLFWR----Q-E--SDCLVN-------TSQISYTIYSTCGAFYIPSVLLIILYSRIYRAARSRKLALERKRISAARERKATKTLGIILGAFIVCWLPFFVVSLVLPICRDSWIHPALFDFFTWLGYLNSLINPIIYTVFNEDFRQAFQKVVHFRKAS---------------------
>NMBR_MOUSE
SNLSFPTEANESELVPEVWEKDFLPDSDGTTAELVIRCVIPSLYLIIISVGLLGNIMLVKIFLTNSAMRNVPNIFISNLAAGDLLLLLTCVPVDASRYFFDEWVFGKLGCKLIPAIQLTSVGVSVFTLTALSADRYRAIVNPMMQTSGVLLWTSLKAVGIWVVSVLLAVPEAVFSEVARSS--FTACIPY---QTDELHPKIHSVLIFLVYFLIPLVIISIYYYHIAKTLIKSLPGNEHTKKQMETRKRLAKIVLVFVGCFVFCWFPNHVLYLYRSFNYKELGHMIVTLVARVLSFSNSCVNPFALYLLSESFRKHFNSQLCCGRKSYPERSTSYLMTSLKSNTKNVV
>OPSD_SHEEP
MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPQGMQCSCGALYFTLKPEINNESFVIYMFVVHFSIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKSSSVYNPVIYIMMNKQFRNCMLTTLCCGKNPLGDDEASTTVSKTETSQVAPA
>OPSD_GALML
MNGTEGENFYVPMSNKTGVVRNPFEYPQYYLADHWMFAVLAAYMFFLIITGFPVNFLTLFVTIQNKKLRQPLNYILLNLAVANLFMVFGGFTTTLITSMNGYFVFGSTGCNLEGFFATLGGEISLWSLVVLAIERYVVVCKPMSNFRFGSQHAIAGVSLTWVMAMACAAPPLVGWSRYIPEGLQCSCGIDYYTPKPEINNVSFVIYMFVVHFSIPLTIIFFCYGRLVCTVKAAAAQQQESETTQRAEREVTRMVVIMVIGFLICWLPYASVALYIFNNQGSEFGPVFMTIPSFFAKSSALYNPLIYILMNKQFRNCMITTLCCGKNPFEEEESTSASSVSSSQVSPAA
>A1AA_RAT
----------MVLLSENASEGSNCT-HPPAPVNISKAILLGVILGGLIIFGVLGNILVILSVACHRHLHSVTHYYIVNLAVADLLLTSTVLPFSAIFEILGYWAFGRVFCNIWAAVDVLCCTASIMGLCIISIDRYIGVSYPLYPTIVTQRRGVRALLCVWVLSLVISIGPLFGWRQP--AP--TICQIN--------EEPGYVLFSALGSFYVPLAIILVMYCRVYVVAKREAKNFSVRLLKFSREKKAAKTLGIVVGCFVLCWLPFFLVMPIGSFFPDFKPSETVFKIVFWLGYLNSCINPIIYPCSSQEFKKAFQNVLRIQCLRRRQSSKHALSQ----------
>GHSR_RAT
SEEPEP-NVTLDLDWDASPGNDSLPDELLPLFPAPLLAGVTATCVALFVVGISGNLLTMLVVSRFRELRTTTNLYLSSMAFSDLLIFLCMPLDLVRLWQYRPWNFGDLLCKLFQFVSESCTYATVLTITALSVERYFAICFPLAKVVVTKGRVKLVILVIWAVAFCSAGPIFVLVG---T-D-TNECRAT---FAVRSGLLTVMVWVSSVFFFLPVFCLTVLYSLIGRKLWRRD--AVGASLRDQNHKQTVKMLAVVVFAFILCWLPFHVGRYLFSKSFEPQISQYCNLVSFVLFYLSAAINPILYNIMSKKYRVAVFKLLGFESFSQRKLSTLKDKSSINT------
>O1E2_HUMAN
-------------MMGQNQTSISDFLLLGLPIQPEQQNLCYALFLAMYLTTLLGNLLIIVLIRLDSHLHTPVYLFLSNLSFSDLCFSSVTMPKLLQNMQNQDPSIPYADCLTQMYFFLYFSDLESFLLVAMAYDRYVAICFPMYTAIMSPMLCLSVVALSWVLTTFHAMLHTLLMARL---CPHFFCDMSKLACSDTRVNEWVIFIMGGLILVIPFLLILGSYARIVSSI--------LKVPSSKGICKAFSTCGSHLSVVSLFYGTVIGLYLCPSANS---STLKDTVMAMMYTVVTPMLTPFIYSLRNRDMKGALERVICKRKNPFLL------------------
>OPSD_SARSL
MNGTEGPYFYVPMVNTSGIVRSPYEYPQYYLVNPAAYARLGAYMFLLILVGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVLGRLGCNIEGFFATLGGEIALWSLVVLAIERWVVVCKPISNFRFGENHAIMGLAFTWLMALACAAPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFVCHFTVPLMVVFFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVIMMVVAFLVCWLPYASVAWWIFTHQGSEFGPVFMTIPAFFAKSSSIYNPMIYICLNKQFRHCMITTLCCGKNPFEEEEGASTSVSSSSVSPAA-
>DOP2_DROME
SATSATLSPAMVATGGGGTTTPEPDLSEFLEALPNDRVGLLAFLFLFSFATVFGNSLVILAVIRERYLHTATNYFITSLAVADCLVGLVVMPFSALYEVLENWFFGTDWCDIWRSLDVLFSTASILNLCVISLDRYWAITDPFYPMRMTVKRAAGLIAAVWICSSAISFPAIVWWRAARP----YKCTFT--------EHLGYLVFSSTISFYLPLLVMVFTYCRIYRAAVIQMGKLSRKLAKFAKEKKAAKTLGIVMGVFIICWLPFFVVNLLSGFCIECEHEEIVSAIVTWLGWINSCMNPVIYACWSRDFRRAFVRLLCMCCPRKIRRKYQPTFATRRCYSTCSL
>GASR_RABIT
VASLCRPGGPLLNNSGTGNLSCEPPRIRGAGTRELELAIRVTLYAVIFLMSVGGNILIIVVLGLSRRLRTVTNAFLLSLAVSDLLLAVACMPFTLLPNLMGTFIFGTVICKAVSYLMGVSVSVSTLSLVAIALERYSAICRPLARVWQTRSHAARVILATWLLSGLLMVPYPVYTAVQP---V-LQCVHR---WPSARVRQTWSVLLLLLLFFVPGVVMAVAYGLISRELYLGSGPPRPAQAKLLAKKRVVRMLLVIVVLFFMCWLPVYSANTWRAFDGPGALSGAPISFIHLLSYASACVNPLVYCFMHRRFRQACLDTCARCCPRPPRARPRPLPSIASLSRLSYT
>IL8B_CANFA
G--DIDNYTYNTEMPIIPADSAPC----RPESLDINKYAVVVIYVLVFVLNLLGNSLVIMVVLYSRVSHSVTDVYLLNLAIADLLFALTLPIWAVSKVKG--WIFGTPLCKIVSLLKEVNFYSGILLLASISMDRYLAIVHATRRLTQKKHWVKFICLGIWALSLILSLPIFVFRRAIN-YS-SPVCYED-MGTNTTKLRIVMRALPQTFGFIVPLMIMLFCYGLTLRTL---------FEAHMGQKHRAMRVIFAVVLVFLLCWLPYNLVADTLMRLQAINDIGRALDATEILGFFHSCLNPLIYAFIGQKFRHGLLKIMAFHGLISKEYLPKDS------------
>KI01_RAT
---------------MDNTTTTEPPKQPCTRNTLITQQIIPMLYCVVFITGVLLNGISGWIFFYVPS-SKSFIIYLKNIVVADFLMGLTFPFKVLSDSGLGPWQLNVFVFRVSAVIFYVNMYVSIAFFGLISFDRYYKIVKPLVSIVQSVNYSKVLSVLVWVLMLLLAVPNIILTNQS-TN---IQCMEL--KNELGRKWHKASNYVFVSIFWIVFLLLTVFYMAITRKIFK----LKSRKNSISVKRKSSRNIFSIVLAFVACFAPYHVARIPYTKSQTEETLLYTKEFTLLLSAANVCLDPISISSYASRLEKS--------------------------------
>O18793
RSRFIRNTNGSGEEVTTFFDYDYGAPCHKFDVKQIGAQLLPPLYSLVFIFGFVGNMLVVLILINCKKLKSLTDIYLLNLAISDLLFLITLPLWAHSAANE--WVFGNAMCKLFTGLYHIGYLGGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWLVAVFASVPGIIFTKCQE-ED-VYICGPY----FPRGWNNFHTIMRNILGLVLPLLIMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWTPYNIVILLNTFQEFFRQLDQATQVTETLGMTHCCINPIIYAFVGEKFRRYLSMFFRKYITKRFCKQCPVFVTSTNTPSTAEQ
>NY1R_HUMAN
TLFSQVENHSVHSNFSEKNAQLLAFENDDCHLPLAMIFTLALAYGAVIILGVSGNLALIIIILKQKEMRNVTNILIVNLSFSDLLVAIMCLPFTFVYTLMDHWVFGEAMCKLNPFVQCVSITVSIFSLVLIAVERHQLIINPR-GWRPNNRHAYVGIAVIWVLAVASSLPFLIYQVMTDKD--KYVCFDQ---FPSDSHRLSYTTLLLVLQYFGPLCFIFICYFKIYIRLKRRNMMMRDNKYRSSETKRINIMLLSIVVAFAVCWLPLTIFNTVFDWNHQICNHNLLFLLCHLTAMISTCVNPIFYGFLNKNFQRDLQFFFNFCDFRSRDDDY---TMHTDVSKTSLK
>B2AR_FELCA
----MGQPGNRSVFLLAPNGSHAPDQDGTQERNDAWVVGMGIVMSLIVLAIVFGNVLVITAIARFERLQTVTNYFITSLACADLVMGLAVVPFGASHILMKMWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFYQSLLTKNKARVVILMVWIVSGLTSFLPIQMHWYRAI-N--TCCDFF--------TNQAYAIASSIVSFYLPLVVMVFVYSRVFQVAQRQDGRHRRASKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNL-IPKEVYILLNWVGYVNSAFNPLIYCRSP-DFRIAFQELLCLRRSSLKAYGNGYSDYAGEHSGGPLG
>OLF1_RAT
-------------MTEENQTVISQFLLLFLPIPSEHQHVFYALFLSMYLTTVLGNLIIIILIHLDSHLHTPMYLFLSNLSFSDLCFSSVTMPKLLQNMQSQVPSIPFAGCLTQLYFYLYFADLESFLLVAMAYDRYVAICFPLYMSIMSPKLCVSLVVLSWVLTTFHAMLHTLLMARL---SPHFFCDISKLSCSDTHVNELVIFVMGGLVIVIPFVLIIVSYARVVASI--------LKVPSVRGIHKIFSTCGSHLSVVSLFYGTIIGLYLCPSANN---STVKETVMAMMYTVVTPMLNPFIYSLRNRDMKEALIRVLCKKKITFCL------------------
>GASR_HUMAN
GASLCRPGAPLLNSSSVGNLSCEPPRIRGAGTRELELAIRITLYAVIFLMSVGGNMLIIVVLGLSRRLRTVTNAFLLSLAVSDLLLAVACMPFTLLPNLMGTFIFGTVICKAVSYLMGVSVSVSTLSLVAIALERYSAICRPLARVWQTRSHAARVIVATWLLSGLLMVPYPVYTVVQP---V-LQCVHR---WPSARVRQTWSVLLLLLLFFIPGVVMAVAYGLISRELYLGPGPSRPTQAKLLAKKRVVRMLLVIVVLFFLCWLPVYSANTWRAFDGPGALSGAPISFIHLLSYASACVNPLVYCFMHRRFRQACLETCARCCPRPPRARPRALPSIASLSRLSYT
>CKR2_MOUSE
FTRSIQELDEGATTPYDYDDGEPC---HKTSVKQIGAWILPPLYSLVFIFGFVGNMLVIIILIGCKKLKSMTDIYLLNLAISDLLFLLTLPFWAHYAANE--WVFGNIMCKVFTGLYHIGYFGGIFFIILLTIDRYLAIVHAVALKARTVTFGVITSVVTWVVAVFASLPGIIFTKSKQ-DD-HYTCGPY----FTQLWKNFQTIMRNILSLILPLLVMVICYSGILHTL--------FRCRNEKKRHRAVRLIFAIMIVYFLFWTPYNIVLFLTTFQESLKHLDQAMQVTETLGMTHCCINPVIYAFVGEKFRRYLSIFFRKHIAKRLCKQCPVFV-------SSTF
>CKRA_MOUSE
PTEQVSWGLYSGYDEEAYSVGPLPELCYKADVQAFSRAFQPSVSLMVAVLGLAGNGLVLATHLAARRTRSPTSVHLLQLALADLLLALTLPFAAAGALQG--WNLGSTTCRAISGLYSASFHAGFLFLACINADRYVAIARALGQRPSTPSRAHLVSVFVWLLSLFLALPALLFS--RD-RE-QRRCRLIFPESLTQTVKGASAVAQVVLGFALPLGVMAACYALLGRTL---------LAARGPERRRALRVVVALVVAFVVLQLPYSLALLLDTADLLAKRKDLALLVTGGLTLVRCSLNPVLYAFLGLRFRRDLRRLLQGGGCSPKPNPRGRCSCSAPTETHSLS
>OPSD_GOBNI
MNGTEGPFFYIPMVNTTGVVRSPYEYPQYYLVNPAAYACLGAYMFFLILVGFPVNFLTLYVTLEHKKLRTPLNYILLNLAVADLFMVFGGFTTTIYTSMHGYFVLGRLGCNIEGFFATLGGEIALWSLVVLAIERWVVVCKPISNFRFGENHAIMGVAFTWFMASACAVPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFTVHFCIPLAVVGFCYGRLLCAVKEAAAAQQESETTQRAEREVSRMVVIMVIGFLVCWLPYASVAWYIFTHQGSEFGPLFMTIPAFFAKSSSIYNPMIYICMNKQFRHCMITTLCCGKNPFEEEEGASTVSSSSVSPAA--
>AG2R_PIG
------MILNSSTEDSIKRIQDDC---PKAGRHNYIFVMIPTLYSIIFVVGIFGNSLVVIVIYFYMKLKTVASVFLLNLALADLCFLLTLPLWAVYTAMEYRWPFGNYLCKIASASVSFNLYASVFLLTCLSIDRYLAIVHPMSRLRRTMLVAKVTCIIIWLLAGLASLPTIIHRNVFF-TN--TVCAFH-YESQNSTLPVGLGLTKNILGFLFPFLIILTSYTLIWKALKKA----YEIQKNKPRNDDIFKIIMAIVLFFFFSWVPHQIFTFLDVLIQLGDIVDTAMPITICLAYFNNCLNPLFYGFLGKKFKKYFLQLLKYIPPKAKSHSSLSTRPSE-------N
>BONZ_MACMU
-----MAEYDHYEDDGFLNSFNDSSQEEHQDFLQFRKVFLPCMYLVVFVCGLVGNSLVLVISIFYHKLQSLTDVFLVNLPLADLVFVCTLPFWAYAGIHE--WIFGQVMCKTLLGVYTINFYTSMLILTCITVDRFIVVVKATNQQAKRMTWGKVICLLIWVISLLVSLPQIIYGNVFNLD--KLICGYH-----DEEISTVVLATQMTLGFFLPLLAMIVCYSVIIKTL---------LHAGGFQKHRSLKIIFLVMAVFLLTQTPFNLVKLIRSTHWEYTSFHYTIIVTEAIAYLRACLNPVLYAFVSLKFRKNFWKLVKDIGCLPYLGVSHQWKTFSASHNVEAT
>PE24_RABIT
--------------------MSTPVANASASSMPELLNNPVTIPAVMFIFGVVGNLVAIVVLCKSRKKETTFYTLVCGLAVTDLLGTLLVSPVTIATYMKGQWPGGQALCDYSTFILLFFGLSGLSIICAMSIERYLAINHAYYSHYVDKRLAGLTLFAVYASNVLFCALPNMGLGRSRLQFPDTWCFID---WRTNVTAHAAFSYMYAGFSSFLILATVLCNVLVCGALLRMAAAAASFRRIAGAEIQMVILLIATSLVVLICSIPLVVRVFINQLYQPDEISQNPDLQAIRIASVNPILDPWIYILLRKTVLSKAIEKIKCLFCRIGGSRRDRSRRTSSAMSTHSR
>FSHR_MACFA
FDMTYAEFDYDLCNEVVDVTCSPKPDAFNPCEDILGYNILRVLIWFISILAITGNIIVLVTLTTSQYKLTVPRFLMCNLAFADLCIGIYLLLIASVDIHTKSWQTG-AGCDAAGFFTVFASELSVYTLTAITLERWHTITHAMLDCKVHVRHAASVMVMGWIFAFAAALFPIFGISSY----KVSICLPM-----IDSPLSQLYVMSLLVLNVLAFVVICGCYTHIYLTVRNP------NIVSSSSDTRIAKRMAMLIFTDFLCMAPISFFAISASLKVPLITVSKAKILLVLFYPINSCANPFLYAIFTKNFRRDFFILLSKFGCYEMQAQIYRTNSHPRNGHCSSA
>ACM1_DROME
-------------MYGNQTNGTIGFETKGPRYSLASMVVMGFVAAILSTVTVAGNVMVMISFKIDKQLQTISNYFLFSLAIADFAIGTISMPLFAVTTILGYWPLGPIVCDTWLALDYLASNASVLNLLIINFDRYFSVTRPLYRAKGTTNRAAVMIG-AWGISLLLWPPWIYSWPYIE-VP-KDECYIQ-----FIETNQYITFGTALAAFYFPVTIMCFLYWRIWRETKKRNPNKKKKSQEKRQESKAAKTLSAILLSFIITWTPYNILVLIKPLTTCSCIPTELWDFFYALCYINSTINPMSYALCNATFRRTYVRILTCKWHTRNREGMVRG------------
>BRS3_MOUSE
QTLISITNDTETSSSVVSNDTTHKGWTGDNSPGIEALCAIYITYAGIISVGILGNAILIKVFFKTKSMQTVPNIFITSLAFGDLLLLLTCVPVDATHYLAEGWLFGKVGCKVLSFIRLTSVGVSVFTLTILSADRYKAVVKPLRQPPNAILKTCAKAGGIWIVSMIFALPEAIFSNVYTVT--FESCNSY---ISERLLQEIHSLLCFLVFYIIPLSIISVYYSLIARTLYKSIPTQSHARKQIESRKRIAKTVLVLVALFALCWLPNHLLYLYHSFTYESDVPFVIIIFSRVLAFSNSCVNPFALYWLSKTFQQHFKAQLCCLKAEQPEPPLGDIMGRVPATGSAHV
>OPSD_MUGCE
MNGTEGPYFYIPMVNTTGIVRSPYEYPQYYLVNPAAYAALGAYMFLLILLGFPINFLTLYVTIEHKKLRTPLNYILLNLAVANLFMVFGGFTTTMYTSMHGYFVLGRLGCNLEGFFATLGGEIALWSLVVLAVERWMVVCKPISNFRFGENHAIMGLAFTWVMASACAVPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFVCHFLIPLVVVFFCYGRLLCAVKEAAAAQQESETTQRAEREVSRMVVIMVVAFLICWCPYAGVAWYIFTHQGSEFGPLFMTFPAFFAKSSSIYNPMIYICMNKQFRHCMITTLCCGKNPFEEEEGASTSVSSSSVSPAA-
>C5AR_HUMAN
SFNYTTPDYGHYDDKDTLDLNTPVD--KTSNTLRVPDILALVIFAVVFLVGVLGNALVVWVTAFEAK-RTINAIWFLNLAVADFLSCLALPILFTSIVQHHHWPFGGAACSILPSLILLNMYASILLLATISADRFLLVFKPICQNFRGAGLAWIACAVAWGLALLLTIPSFLYRVVRE-PP--VLCGVD-YS-HDKRRERAVAIVRLVLGFLWPLLTLTICYTFILLRT---------WSRRATRSTKTLKVVVAVVASFFIFWLPYQVTGIMMSFLEPSLLLNKLDSLCVSFAYINCCINPIIYVVAGQGFQGRLRKSLPSLLRNVLTEESVVRSTVD--------
>MAS_RAT
------MDQSNMTSFAEEKAMNTSSRNASLGTSHPPIPIVHWVIMSISPLGFVENGILLWFLCFRMR-RNPFTVYITHLSIADISLLFCIFILSIDYALDYESSGHYYTIVTLSVTFLFGYNTGLYLLTAISVERCLSVLYPIYRCHRPKHQSAFVCALLWALSCLVTTMEYVMCI-DSHS--QSDC-----------RAVIIFIAILSFLVFTPLMLVSSTILVVKIRK----------NTWASHSSKLYIVIMVTIIIFLIFAMPMRVLYLLYYEYW--STFGNLHNISLLFSTINSSANPFIYFFVGSSKKKRFRESLKVVLTRAFKDEMQPRTVSIETVV----
>GPRC_RAT
NLSGLPRDCIEAGTPENISAAVPSQGSVVESEPELVVNPWDIVLCSSGTLICCENAVVVLIIFHSPSLRAPMFLLIGSLALADLLAG-LGLIINFVFA-Y--LLQSEATKLVTIGLIVASFSASVCSLLAITVDRYLSLYYALYHSERTVTFTYVMLVMLWGTSTCLGLLPVMGWNCL-R--DESTCSVV-------RPLTKNNAAILSISFLFMFALMLQLYIQICKIVMRHIALHFLATSHYVTTRKGISTLALILGTFAACWMPFTLYSLIADY----TYPSIYTYATLLPATYNSIINPVIYAFRNQEIQKALCLICCGCIPNTLSQRARSP------------
>V1AR_HUMAN
SSPWWPLATGAGNTSREAEALGEGNGPPRDVRNEELAKLEIAVLAVTFAVAVLGNSSVLLALHRTPRKTSRMHLFIRHLSLADLAVAFFQVLPQMCWDITYRFRGPDWLCRVVKHLQVFGMFASAYMLVVMTADRYIAVCHPLKTLQQPARRSRLMIAAAWVLSFVLSTPQYFVFSMIETK--ARDCWAT---FIQPWGSRAYVTWMTGGIFVAPVVILGTCYGFICYNIWCNFLLVSSVKSISRAKIRTVKMTFVIVTAYIVCWAPFFIIQMWSVWDPMSESENPTITITALLGSLNSCCNPWIYMFFSGHLLQDCVQSFPCCQNMKEKFNKEDTTFYSNNRSPTNS
>A1AA_HUMAN
----------MVFLSGNASDSSNCT-QPPAPVNISKAILLGVILGGLILFGVLGNILVILSVACHRHLHSVTHYYIVNLAVADLLLTSTVLPFSAIFEVLGYWAFGRVFCNIWAAVDVLCCTASIMGLCIISIDRYIGVSYPLYPTIVTQRRGLMALLCVWALSLVISIGPLFGWRQP--AP--TICQIN--------EEPGYVLFSALGSFYLPLAIILVMYCRVYVVAKREAKTFSVRLLKFSREKKAAKTLGIVVGCFVLCWLPFFLVMPIGSFFPDFKPSETVFKIVFWLGYLNSCINPIIYPCSSQEFKKAFQNVLRIQCLCRKQSSKHALSQ----------
>CML2_RAT
PSNSTPLALNLSLALREDAPGNLTGDLSEHQQYVIALFLSCLYTIFLFPIGFVGNILILVVNISFREKMTIPDLYFINLAAADLILVADSLIEVFNLDEQ--YYDIAVLCTFMSLFLQINMYSSVFFLTWMSFDRYLALAKAMCGLFRTKHHARLSCGLIWMASVSATLVPFTAVHLR-A----CFCFA---------DVREVQWLEVTLGFIVPFAIIGLCYSLIVRALIR----AHRHRGLRPRRQKALRMIFAVVLVFFICWLPENVFISVHLLQWAQHAYPLTGHIVNLAAFSNSCLSPLIYSFLGETFRDKLRLYVAQKTSLPALNRFCHADSTEQSDVKFSS
>GPRF_HUMAN
----MDPEETSVYLDYYYATSPNSDIRETHSHVPYTSVFLPVFYTAVFLTGVLGNLVLMGALHFKPGSRRLIDIFIINLAASDFIFLVTLPLWVDKEASLGLWRTGSFLCKGSSYMISVNMHCSVLLLTCMSVDRYLAIVWPVSRKFRRTDCAYVVCASIWFISCLLGLPTLLSRELT-IDD-KPYCAEK----KATPIKLIWSLVALIFTFFVPLLSIVTCYCCIARKLCAH---YQQSGKHNKKLKKSIKIIFIVVAAFLVSWLPFNTFKFLAIVSGLRAILQLGMEVSGPLAFANSCVNPFIYYIFDSYIRRAIVHCLCPCLKNYDFGSSTETALSTFIHAEDFA
>CKR5_TRAFR
----MDYQVSSPTYDIDYYTSEPC---QKVNVKQIAARLLPPLYSLVFIFGFVGNILVVLILINCKRLKSMTDIYLLNLAISDLFFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWVVAVFASLPGIIFTRSQR-EG-HYTCSSHFPYSQYQFWKNFQTLKIVILGLVLPLLVMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCSIFA-------SSVY
>GP68_HUMAN
----------------MGNITADNSSMSCTIDHTIHQTLAPVVYVTVLVVGFPANCLSLYFGYLQIKARNELGVYLCNLTVADLFYICSLPFWLQYVLQHDNWSHGDLSCQVCGILLYENIYISVGFLCCISVDRYLAVAHPFFHQFRTLKAAVGVSVVIWAKELLTSIYFLMHEEV--QH---RVCFEHPIQAWQRAIQRAINYYRFLVGFLFPICLLLASYQGILRAVRR------SHGTQKSRKDQIQRLVLSTVVIFLACFLPYHVLLLVRSVWEASKGVFNAYHFSLLLTSFNCVADPVLYCFVSETTHRDLARLRGACLAFLTCSRTGRAAPEASGKSGAQG
>YYI3_CAEEL
SDPNAEDLYITMTPSVSTENDTTVWATEEPAAIVWRHPLLAIALFSICLLTVAGNCLVVIAVCTKKYIWVTRLYLIISLAIADLIVGVIVMPMNSLFEIANHWLFGLMMCDVFHAMDILASTASIWNLCVISLDRYMAGQDPIYRDKVSKRRILMAILSVWVLSAILSFPGIIWWWRTSP-H--SQCLFT--------DSKMYVSFSSLVSFYIPLFLILFAYGKVYIIATRHKLNKSRQMMRYVHEQRAARTLSIVVGAFILCWTPFFVFTPLTAFCESCSNKETIFTFVTWAGHLNSMLNPLIYSRFSRDFRRAFKQILTCQRQQKVKTAFKTPLISVTQMAPRFS
>CKR1_MACMU
METPNTTEDYDMITEFDYGDATPC---HKVNERAILAQLLPPLYSLVFVIGVVGNLLVVLVLVQYKRLKNMTNIYLLNLAISDLLFLFTLPFLIYYKSTDD-WIFGDAMCKILSGFYYTGLYSEIFFIILLTIDRYLAIVHAVALRARTVTFGVITSIIIWALAILASSPLMYFSKTQW-NI-RHSCNIHFPYESFQQWKLFQALKLNLFGLVLPLLVMIVCYTGIIKIL---------LRRPNEKKSKAVRLIFVIMIIFFLFWTPYNLTELISVFQEFLRQLDLAMEVTEVIANMHCCVNPVIYAFAGERFRKYLRQLFHRRVAVHLVKWLPFLV-------SSTS
>BRB2_RAT
IEMFNITTQALGSAHNGTFSEVNC---PDTEWWSWLNAIQAPFLWVLFLLAALENIFVLSVFCLHKTNCTVAEIYLGNLAAADLILACGLPFWAITIANNFDWLFGEVLCRVVNTMIYMNLYSSICFLMLVSIDRYLALVKTMMGRMRGVRWAKLYSLVIWSCTLLLSSPMLVFRTMKD-HN--TACVIV---YPSRSWEVFTNMLLNLVGFLLPLSIITFCTVRIMQVLRNN---EMKKFKEVQTEKKATVLVLAVLGLFVLCWFPFQISTFLDTLLRLGRAVDIVTQISSYVAYSNSCLNPLVYVIVGKRFRKKSREVYQAICRKGGCMGESVQLRTS-------I
>A1AD_RAT
TGSGEDNQSSTGEPGAAASGEVNGSAAVGGLVVSAQGVGVGVFLAAFILTAVAGNLLVILSVACNRHLQTVTNYFIVNLAVADLLLSAAVLPFSATMEVLGFWAFGRTFCDVWAAVDVLCCTASILSLCTISVDRYVGVRHSLYPAIMTERKAAAILALLWAVALVVSVGPLLGWKEP--VP--RFCGIT--------EEVGYAIFSSVCSFYLPMAVIVVMYCRVYVVARSTHTLLSVRLLKFSREKKAAKTLAIVVGVFVLCWFPFFFVLPLGSLFPQLKPSEGVFKVIFWLGYFNSCVNPLIYPCSSREFKRAFLRLLRCQCRRRRRR--LWAASTGD-ARSDCA
>TSHR_BOVIN
LQAFDSHYDYTVCGGSEDMVCTPKSDEFNPCEDIMGYKFLRIVVWFVSLLALLGNVFVLVILLTSHYKLTVPRFLMCNLAFADFCMGLYLLLIASVDLYTQSWQTG-PGCNTAGFFTVFASELSVYTLTVITLERWHAITFAMLDRKIRLWHAYVIMLGGWVCCFLLALLPLVGISSY----KVSICLPM-----TETPLALAYIILVLLLNIIAFIIVCACYVKIYITVRNP------HYNPGDKDTRIAKRMAVLIFTDFMCMAPISFYALSALMNKPLITVTNSKILLVLFYPLNSCANPFLYAIFTKAFQRDVFMLLSKFGICKRQAQAYRGSTGIRVQKVPPD
>MC3R_HUMAN
NASCCLPSVQPTLPNGSEHLQAPFFSNQSSSAFCEQVFIKPEIFLSLGIVSLLENILVILAVVRNGNLHSPMYFFLCSLAVADMLVSVSNALETIMIAIVHSDQFIQHMDNIFDSMICISLVASICNLLAIAVDRYVTIFYALYHSIMTVRKALTLIVAIWVCCGVCGVVFIVYYS----------------------EESKMVIVCLITMFFAMMLLMGTLYVHMFLFARLHIAAADGVAPQQHSCMKGAVTITILLGVFIFCWAPFFLHLVLIITCPTNICYTAHFNTYLVLIMCNSVIDPLIYAFRSLELRNTFREILCGCNGMNLG------------------
>MSHR_SHEEP
PVLGSQRRLLGSLNCTPPATLPLTLAPNRTGPQCLEVSIPDGLFLSLGLVSLVENVLVVAAIAKNRNLHSPMYYFICCLAMSDLLVSVSNVLETAVMLLLEAAAVVQQLDNVIDVLICSSMVSSLCFLGAIAVDRYISIFYALYHSVVTLPRAWRIIAAIWVASILTSVLSITYYY----------------------NNHTVVLLCLVGFFIAMLALMAVLYVHMLARACQHIARKRQRPIHQGFGLKGAATLTILLGVFFLCWGPFFLHLSLIVLCPQHGCIFKNFNLFLALIICNAIVDPLIYAFRSQELRKTLQEVLQCSW-----------------------
>CKR1_MOUSE
MEISDFTEAYPTTTEFDYGDSTPC---QKTAVRAFGAGLLPPLYSLVFIIGVVGNVLMILVLMQHRRLQSMTSIYLFNLAVSDLVFLFTLPFWIDYKLKDD-WIFGDAMCKLLSGFYYLGLYSEIFFIILLTIDRYLAIVHAVALRARTVTLGIITSIITWALAILASMPALYFFKAQW-EF-HRTCSPHFPYKSLKQWKRFQALKLNLLGLILPLLVMIICYAGIIRIL---------LRRPSEKKVKAVRLIFAITLLFFLLWTPYNLSVFVSAFQDVLKHLDLAMQVTEVIAYTHCCVNPIIYVFVGERFWKYLRQLFQRHVAIPLAKWLPFLT-------SSIS
>A2AA_BOVIN
----MGSLQPDAGNASWNGTEAPGGGARATPYSLQVTLTLVCLAGLLMLFTVFGNVLVIIAVFTSRALKAPQNLFLVSLASADILVATLVIPFSLANEVMGYWYFGKAWCEIYLALDVLFCTSSIVHLCAISLDRYWSITQAIYNLKRTPRRIKAIIVTVWVISAVISFPPLISFEKKRP-S--PRCEIN--------DQKWYVISSSIGSFFAPCLIMILVYVRIYQIAKRRSGGASRWRGRQNREKRFTFVLAVVIGVFVVCWFPFFFTYTLTAIGCP--VPPTLFKFFFWFGYCNSSLNPVIYTIFNHDFRRAFKKILCRGDRKRIV------------------
>SSR3_MOUSE
SEPMTLDPGNTSSTWPLDTTLGNTSAGASLTGLAVSGILISLVYLVVCVVGLLGNSLVIYVVLRHTSSPSVTSVYILNLALADELFMLGLPFLAAQNALSY-WPFGSLMCRLVMAVDGINQFTSIFCLTVMSVDRYLAVVHPTSARWRTAPVARTVSRAVWVASAVVVLPVVVFSGVPR-----STCHMQ-WPEPAAAWRTAFIIYMAALGFFGPLLVICLCYLLIVVKVRSTSCQAPACQRRRRSERRVTRMVVAVVALFVLCWMPFYLLNIVNVVCPLPPAFFGLYFLVVALPYANSCANPILYGFLSYRFKQGFRRILLRPSRRIRSQEPGSGEEDEEEEERREE
>P2YR_BOVIN
TAFLADPGSPWGNSTVTSTAAVASPFKCALTKTGFQFYYLPAVYILVFIIGFLGNSVAIWMFVFHMKPWSGISVYMFNLALADFLYVLTLPALIFYYFNKTDWIFGDAMCKLQRFIFHVNLYGSILFLTCISAHRYSGVVYPLSLGRLKKKNAVYISVLVWLIVVVGISPILFYSGTG----KTITCYDT-TSDEYLRSYFIYSMCTTVAMFCVPLVLILGCYGLIVRALIY------KDLDNSPLRRKSIYLVIIVLTVFAVSYIPFHVMKTMNLRARLDDRVYATYQVTRGLASLNSCVDPILYFLAGDTFRRRLSRATRKASRRSEANLQSKSLSEFKQNGDTSL
>5H5B_RAT
TPGIAFPPGPESCSDSPSSGRGGLILSGREPPFSAFTVLVVTLLVLLIAATFLWNLLVLVTILRVRAFHRVPHNLVASTAVSDVLVAALVMPLSLVSELSAGWQLGRSLCHVWISFDVLCCTASIWNVAAIALDRYWTITRHLYTLRTRRRASALMIAITWALSALIALAPLLFGWGE-D-A-LQRCQVS--------QEPSYAVFSTCGAFYVPLAVVLFVYWKIYKAAKFRATVTSGDSWREQKEKRAAMMVGILIGVFVLCWIPFFLTELVSPLCAC-SLPPIWKSIFLWLGYSNSFFNPLIYTAFNKNYNNAFKSLFTKQR-----------------------
>OPSD_CHICK
MNGTEGQDFYVPMSNKTGVVRSPFEYPQYYLAEPWKFSALAAYMFMLILLGFPVNFLTLYVTIQHKKLRTPLNYILLNLVVADLFMVFGGFTTTMYTSMNGYFVFGVTGCYIEGFFATLGGEIALWSLVVLAVERYVVVCKPMSNFRFGENHAIMGVAFSWIMAMACAAPPLFGWSRYIPEGMQCSCGIDYYTLKPEINNESFVIYMFVVHFMIPLAVIFFCYGNLVCTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTNQGSDFGPIFMTIPAFFAKSSAIYNPVIYIVMNKQFRNCMITTLCCGKNPLGDEDTSAGTSSVSTSQVSPA
>B3AR_MOUSE
MAPWPHRNGSLALWSDAPTLDPSAAN---TSGLPGVPWAAALAGALLALATVGGNLLVIIAIARTPRLQTITNVFVTSLAAADLVVGLLVMPPGATLALTGHWPLGETGCELWTSVDVLCVTASIETLCALAVDRYLAVTNPLYGTLVTKRRARAAVVLVWIVSAAVSFAPIMSQWWRVQ-E--RCCSFA--------SNMPYALLSSSVSFYLPLLVMLFVYARVFVVAKRQGVPRRPARLLPLREHRALRTLGLIMGIFSLCWLPFFLANVLRALAGPS-VPSGVFIALNWLGYANSAFNPVIYCRSP-DFRDAFRRLLCSYGGRGPEEPRAVTARQSPPLNRFDG
>SSR1_RAT
GEGVCSRGPGSGAADGMEEPGRNSSQNGTLSEGQGSAILISFIYSVVCLVGLCGNSMVIYVILRYAKMKTATNIYILNLAIADELLMLSVPFLVTSTLLRH-WPFGALLCRLVLSVDAVNMFTSIYCLTVLSVDRYVAVVHPIAARYRRPTVAKVVNLGVWVLSLLVILPIVVFSRTAADG--TVACNML-MPEPAQRWLVGFVLYTFLMGFLLPVGAICLCYVLIIAKMRMVALK-AGWQQRKRSERKITLMVMMVVMVFVICWMPFYVVQLVNVFAEQ--DDATVSQLSVILGYANSCANPILYGFLSDNFKRSFQRILCLS-----WMDNAAETALKSRAYSVED
>AA1R_HUMAN
----------------------------MPPSISAFQAAYIGIEVLIALVSVPGNVLVIWAVKVNQALRDATFCFIVSLAVADVAVGALVIPLAILINIG--PQTYFHTCLMVACPVLILTQSSILALLAIAVDRYLRVKIPLYKMVVTPRRAAVAIAGCWILSFVVGLTPMFGWNNLSN-G-VIKCEFE-----KVISMEYMVYFNFFVWVLPPLLLMVLIYLEVFYLIRKQK-VSGDPQKYYGKELKIAKSLALILFLFALSWLPLHILNCITLFCPSCHKPSILTYIAIFLTHGNSAMNPIVYAFRIQKFRVTFLKIWNDHFRCQPAPPID--PDD---------
>PE23_HUMAN
PFCTRLNHSYTGMWAPERSAEARGNLTRPPGSGEDCGSVSVAFPITMLLTGFVGNALAMLLVSRSYRRKKSFLLCIGWLALTDLVGQLLTTPVVIVVYLSKQIDPSGRLCTFFGLTMTVFGLSSLFIASAMAVERALAIRAPHYASHMKTRATRAVLLGVWLAVLAFALLPVLGVG----QYTGTWCFISNGTSSSHNWGNLFFASAFAFLGLLALTVTFSCNLATIKALVSRAKASQSSAQWGRITTETAIQLMGIMCVLSVCWSPLLIMMLKMIFNQTSQKECNFFLIAVRLASLNQILDPWVYLLLRKILLRKFCQIRYHTNNYASSSTSLPCWSDHLER-----
>OLF0_RAT
---------------MNNQTFITQFLLLGLPIPEEHQHLFYALFLVMYLTTILGNLLIIVLVQLDSQLHTPMYLFLSNLSFSDLCFSSVTMPKLLQNMRSQDTSIPYGGCLAQTYFFMVFGDMESFLLVAMAYDRYVAICFPLYTSIMSPKLCTCLVLLLWMLTTSHAMMHTLLAARL---SLNFFCDLFKLACSDTYINELMIFIMSTLLIIIPFFLIVMSYARIISSI--------LKVPSTQGICKVFSTCGSHLSVVSLFYGTIIGLYLCPAGNN---STVKEMVMAMMYTVVTPMLNPFIYSLRNRDMKRALIRVICSMKITL--------------------
>BRB2_RABIT
-MLNITSQVLAPALNGSVSQSSGC---PNTEWSGWLNVIQAPFLWVLFVLATLENLFVLSVFCLHKSSCTVAEVYLGNLAAADLILACGLPFWAVTIANHFDWLFGEALCRVVNTMIYMNLYSSICFLMLVSIDRYLALVKTMIGRMRRVRWAKLYSLVIWGCTLLLSSPMLVFRTMKD-YN--TACIID---YPSRSWEVFTNVLLNLVGFLLPLSVITFCTVQILQVLRNN---EMQKFKEIQTERRATVLVLAVLLLFVVCWLPFQVSTFLDTLLKLGHVIDVITQVGSFMGYSNSCLNPLVYVIVGKRFRKKSREVYRAACPKAGCVLEPVQLRTS-------I
>CCKR_MOUSE
DSLLMNGSNITPPCELGLENETLFCLDQPQPSKEWQSAVQILLYSFIFLLSVLGNTLVITVLIRNKRMRTVTNIFLLSLAVSDLMLCLFCMPFNLIPNLLKDFIFGSAVCKTTTYFMGTSVSVSTFNLVAISLERYGAICRPLSRVWQTKSHALKVIAATWCLSFTIMTPYPIYN----NNQTANMCRFL---LPSDAMQQSWQTFLLLILFLIPGVVMVVAYGLISLELYQGRINSSGSAANLIAKKRVIRMLIVIVVLFFLCWMPIFSANAWRAYDTVSHLSGTPISFILLLSYTSSCVNPIIYCFMNKRFRLGFMATFPCCPNPGPTGVRGEVTIRASLSRYSYS
>ACM3_RAT
N---------ISQETGNFSS-NDTSSDPLGGHTIWQVVFIAFLTGFLALVTIIGNILVIVAFKVNKQLKTVNNYFLLSLACADLIIGVISMNLFTTYIIMNRWALGNLACDLWLSIDYVASNASVMNLLVISFDRYFSITRPLYRAKRTTKRAGVMIGLAWVISFVLWAPAILFWQYFV-VP-PGECFIQ------FLSEPTITFGTAIAAFYMPVTIMTILYWRIYKETEKRKTRTKRKRMSLIKEKKAAQTLSAILLAFIITWTPYNIMVLVNTFCDSC-IPKTYWNLGYWLCYINSTVNPVCYALCNKTFRTTFKTLLLCQCDKRKRRKQQYQHKRVPEQAL---
>5H1E_HUMAN
-------------MNITNCT--TEASMAIRPKTITEKMLICMTLVVITTLTTLLNLAVIMAIGTTKKLHQPANYLICSLAVTDLLVAVLVMPLSIIYIVMDRWKLGYFLCEVWLSVDMTCCTCSILHLCVIALDRYWAITNAIYARKRTAKRAALMILTVWTISIFISMPPLFWRS---S-P--SQCTIQ-------HDHVIYTIYSTLGAFYIPLTLILILYYRIYHAAKSLNDLGERQQISSTRERKAARILGLILGAFILSWLPFFIKELIVGLSIYT-VSSEVADFLTWLGYVNSLINPLLYTSFNEDFKLAFKKLIRCREHT---------------------
>OPS1_CALVI
ALT---NGSVTDKVTPDMAHLVHPYWNQFPAMEPKWAKFLAAYMVLIATISWCGNGVVIYIFSTTKSLRTPANLLVINLAISDFGIMITNTPMMGINLFYETWVLGPLMCDIYGGLGSAFGCSSILSMCMISLDRYNVIVKGMAGQPMTIKLAIMKIALIWFMASIWTLAPVFGWSRYVPEGNLTSCGID--YLERDWNPRSYLIFYSIFVYYLPLFLICYSYWFIIAAVSAHMNVRSSEDADKSAEGKLAKVALVTISLWFMAWTPYTIINTLGLFKYEG-LTPLNTIWGACFAKSAACYNPIVYGISHPKYGIALKEKCPCCVFGKVDDGK-ASNNESETKA----
>V1BR_RAT
---MNSEPSWTATPSPGGTLPVPNATTPWLGRDEELAKVEIGILATVLVLATGGNLAVLLTLGRHGHKRSRMHLFVLHLALTDLGVALFQVLPQLLWDITYRFQGSDLLCRAVKYLQVLSMFASTYMLLAMTLDRYLAVCHPLRSLRQPSQSTYPLIAAPWLLAAILSLPQVFIFSLRESG--VLDCWAD---FYFSWGPRAYITWTTMAIFVLPVAVLSACYGLICHEIYKNRGLVSSISTISRAKIRTVKMTFVIVLAYIACWAPFFSVQMWSVWDENADSTNVAFTISMLLGNLSSCCNPWIYMGFNSRLLPRSLSHHACCTGSKPQVHRQLSRTTLLTHACGSP
>OPSD_SARPI
MNGTEGPFFYIPMSNATGLVRSPYDYPQYYLVPPWGYACLAAYMFLLILTGFPVNFLTLYVTIEHKKLRSPLNYILLNLAVADLFMVIGGFTTTMWTSLNGYFVFGRMGCNIEGFFATLGGEIALWSLVVLSMERWIVVCKPISNFRFGENHAVMGVAFSWFMAAACAVPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFVVHFTCPLTIITFCYGRLVCTVKEAAAQQQESETTQRAEREVTRMVIIMFVAFLACWVPYASVAWYIFTHQGSEFGPVFMTIPAFFAKSSAVYNPVIYICLNKQFRHCMITTLCCGKNPFEEEEGSTTSVCSVSPAA---
>V1BR_MOUSE
---MDSEPSWTATPSPGGTLFVPNTTTPWLGRDEELAKVEIGILATVLVLATGGNLAVLLILGLQGHKRSRMHLFVLHLALTDLGVALFQVLPQLLWDITYRFQGSDLLCRAVKYLQVLSMFASTYMLLAMTLDRYLAVCHPLRSLQQPSQSTYPLIAAPWLLAAILSLPQVFIFSLRESG--VLDCWAD---FYFSWGPRAYITWTTMAIFVLPVVVLTACYGLICHEIYKNATLVSSISTISRAKIRTVKMTFVIVLAYIACWAPFFSVQMWSVWDENADSTNVAFTISMLLGNLSSCCNPWIYMGFNSHLLPRSLSHRACCRGSKPRVHRQLSRTTLLTHTCGPS
>CCKR_CAVPO
DSLFVNGSNITSACELGFENETLFCLDRPRPSKEWQPAVQILLYSLIFLLSVLGNTLVITVLIRNKRMRTVTNIFLLSLAVSDLMLCLFCMPFNLIPSLLKDFIFGSAVCKTTTYFMGTSVSVSTFNLVAISLERYGAICKPLSRVWQTKSHALKVIAATWCLSFTIMTPYPIYN----NNQTGNMCRFL---LPNDVMQQTWHTFLLLILFLIPGIVMMVAYGLISLELYQGRINSSSSTANLMAKKRVIRMLIVIVVLFFLCWMPIFSANAWRAYDTVSHLSGTPISFILLLSYTSSCVNPIIYCFMNKRFRLGFMATFPCCPNPGTPGVRGEMTTGASLSRYSYS
>OPRM_RAT
LSHVDGNQSDPCGLNRTGLGGNDSLCPQTGSPSMVTAITIMALYSIVCVVGLFGNFLVMYVIVRYTKMKTATNIYIFNLALADALATSTLPFQSVNYLMGT-WPFGTILCKIVISIDYYNMFTSIFTLCTMSVDRYIAVCHPVALDFRTPRNAKIVNVCNWILSSAIGLPVMFMATTKYGS---IDCTLT-FSHPTWYWENLLKICVFIFAFIMPVLIITVCYGLMILRLKSVRML-SGSKEKDRNLRRITRMVLVVVAVFIVCWTPIHIYVIIKALITIPTFQTVSWHFCIALGYTNSCLNPVLYAFLDENFKRCFREFCIPTSSTIEQQNSTRVPSTANTVDRTNH
>FMLR_RABIT
-----------MDSNASLPLNVSGGTQATPAGLVVLDVFSYLILVVTFVLGVLGNGLVIWVTGFRMT-HTVTTISYLNLALADFSFTSTLPFFIVTKALGGHWPFGWFLCKFVFTIVDINLFGSVFLIALIALDRCICVLHPVAQNHRNVSLAKKVIVGPWICALLLTLPVIIRVTTLSSPWPAEKLKVA------ISMFMVRGIIRFIIGFSTPMSIVAVCYGLIATKI---------HRQGLIKSSRPLRVLSFVVASFLLCWSPYQIAALIATVRIREKDLRIVLDVTSFVAFFNSCLNPMLYVFMGQDFRERLIHSLPASLERALSEDSAQTTS----------
>OPSP_COLLI
----MDPTNSPQEPPHTSTPGPFDGPQWPHQAPRGMYLSVAVLMGIVVISASVVNGLVIVVSIRYKKLRSPLNYILVNLAMADLLVTLCGSSVSFSNNINGFFVFGKRLCELEGFMVSLTGIVGLWSLAILALERYVVVCRPLGDFRFQHRHAVTGCAFTWVWSLLWTTPPLLGWSSYVPEGLRTSCGPN--WYTGGSNNNSYILTLFVTCFVMPLSLILFSYANLLMTLRAAAAQQQESDTTQQAERQVTRMVVAMVMAFLICWLPYTTFALVVATNKDIAIQPALASLPSYFSKTATVYNPIIYVFMNKQFQSCLLKMLCCGHHPRGTGRTAPAGLRNKVTPSHPV
>O1D2_HUMAN
-------------MDGGNQSEGSEFLLLGMSESPEQQQILFWMFLSMYLVTVVGNVLIILAISSDSRLHTPVYFFLANLSFTDLFFVTNTIPKMLVNLQSHNKAISYAGCLTQLYFLVSLVALDNLILAVMAYDRYVAICCPLYTTAMSPKLCILLLSLCWVLSVLYGLIHTLLMTRV---THYIFCEMYRMACSNIQINHTVLIATGCFIFLIPFGFVIISYVLIIRAI--------LRIPSVSKKYKAFSTCASHLGAVSLFYGTLCMVYLKPLHTY----SVKDSVATVMYAVVTPMMNPFIYSLRNKDMHGALGRLLDKHFKRLT-------------------
>IL8A_HUMAN
PQMWDFDDLNFTGMPPADEDYSPC----MLETETLNKYVVIIAYALVFLLSLLGNSLVMLVILYSRVGRSVTDVYLLNLALADLLFALTLPIWAASKVNG--WIFGTFLCKVVSLLKEVNFYSGILLLACISVDRYLAIVHATRTLTQKRHLVKFVCLGCWGLSMNLSLPFFLFRQAYH-NS-SPVCYEV-LGNDTAKWRMVLRILPHTFGFIVPLFVMLFCYGFTLRTL---------FKAHMGQKHRAMRVIFAVVLIFLLCWLPYNLVLLADTLMRTQNNIGRALDATEILGFLHSCLNPIIYAFIGQNFRHGFLKILAMHGLVSKEFLARHR------------
>EBP2_HUMAN
NPDKDGGTPDSGQELRGNLTGQIQNPLYPVTESSYSAYAIMLLALVVFAVGIVGNLSVMCIVWHSYYLKSAWNSILASLALWDFLVLFFCLPIVIFNEITKQRLLGDVSCRAVPFMEVSSLGVTTFSLCALGIDRFHVATSTLVRPIERCQSILAKLAVIWVGSMTLAVPELLLWQLAQM----KPSASLSLYSLVMTYQNARMWWYFGCYFCLPILFTVTCQLVTWRVRGPPRKS-CRASKHEQCESQLNSTVVGLTVVYAFCTLPENVCNIVVAYLSTEQTLDLLGLINQFSTFFKGAITPVLLLCICRPLGQAFLDCCCCCCCEECGGASEASKLKTEVSSSIYF
>5H5B_MOUSE
TPGLAFPPGPESCSDSPSSGRGGLILPGREPPFSAFTVLVVTLLVLLIAATFLWNLLVLVTILRVRAFHRVPHNLVASTAVSDVLVAVLVMPLSLVSELSAGWQLGRSLCHVWISFDVLCCTASIWNVAAIALDRYWTITRHLYTLRTRSRASALMIAITWALSALIALAPLLFGWGE-D-A-LQRCQVS--------QEPSYAVFSTCGAFYLPLAVVLFVYWKIYKAAKFRATVTSGDSWREQKEKRAAMMVGILIGVFVLCWIPFFLTELISPLCAC-SLPPIWKSIFLWLGYSNSFFNPLIYTAFNKNYNNAFKSLFTKQR-----------------------
>V2R_BOVIN
MFMASTTSAVPWHLSQPTPAGNGSEGELLTARDPLLAQAELALLSTVFVAVALSNGLVLGALVRRGRRWAPMHVFIGHLCLADLAVALFQVLPQLAWDATDRFRGPDALCRAVKYLQMVGMYASSYMILAMTLDRHRAICRPMAHRHGGGTHWNRPVLLAWAFSLLFSLPQLFIFAQRDSG--VLDCWAR---FAEPWGLRAYVTWIALMVFVAPALGIAACQVLIFREIHASGCRPAEGARVSAAVAKTVKMTLVIVIVYVLCWAPFFLVQLWAAWDPEA-REGPPFVLLMLLASLNSCTNPWIYASFSSSISSELRSLLCCTWRRAPPSPGPQESFLAKDTPS---
>GRHR_MOUSE
--MANNASLEQDPNHCSAINNSIPLIQGKLPTLTVSGKIRVTVTFFLFLLSTAFNASFLLKLQKWTQKLSRMKVLLKHLTLANLLETLIVMPLDGMWNITVQWYAGEFLCKVLSYLKLFSMYAPAFMMVVISLDRSLAITQPL-AVQSNSKLEQSMISLAWILSIVFAGPQLYIFRMIYTV--FSQCVTH--SFPQWWHQAFYNFFTFGCLFIIPLLIMLICNAKIIFALTRVPRKNQSKNNIPRARLRTLKMTVAFATSFVVCWTPYYVLGIWYWFDPEMRVSEPVNHFFFLFAFLNPCFDPLIYGYFSL-------------------------------------
>OLF1_CANFA
--------------MDGNYTLVTEFILLGFPTRPELQIVLFLVFLTLYGIILTGNIGLMMLIRTDPHLQTPMYFFLSNLSFADLCFSSAIVPKMLVNFLSENKSISLYGCALQFYFSCAFADTESFILAAMAYDRYVAICNPLYTVVMSRGICVWLIVLSYIGGNMSSLVHTSFAFIL---KNHFFCDLPKLSCTDTSVNEWLLSTYGSSVEIFCFIVIVISYYFILRSV--------LRIRSSSGRKKTFSTCASHLTSVAIYQGTLLFIYSRPTYLY---TPNTDKIISVFYTIIIPVLNPLIYSLRNKDVKDAAKRAVRLKVDSS--------------------
>GPRY_HUMAN
SVSSWPYSSHRMRFITNHSDQATPNVTTCPMDEKLLSTVLTTSYSVIFIVGLVGNIIALYVFLGIHRKRNSIQIYLLNVAIADLLLIFCLPFRIMYHINQNKWTLGVILCKVVGTLFYMNMYISIILLGFISLDRYIKINRSIQRKAITTKQSIYVCCIVWMLALGGFLTMIILTLKK-----STMCFHY--RDKHNAKGEAIFNFILVVMFWLIFLLIILSYIKIGKNLLRISKR--SKFPNSGKYATTARNSFIVLIIFTICFVPYHAFRFIYISSQLNEIVHKTNEIMLVLSSFNSCLDPVMYFLMSSNIRKIMCQLLFRRFQGEPSRSESTSLHDTSVAVKIQS
>OPSD_RANPI
MNGTEGPNFYIPMSNKTGVVRSPFDYPQYYLAEPWKYSVLAAYMFLLILLGLPINFMTLYVTIQHKKLRTPLNYILLNLGVCNHFMVLCGFTITMYTSLHGYFVFGQTGCYFEGFFATLGGEIALWSLVVLAIERYIVVCKPMSNFRFGENHAMMGVAFTWIMALACAVPPLFGWSRYIPEGMQCSCGVDYYTLKPEVNNESFVIYMFVVHFLIPLIIISFCYGRLVCTVKEAAAQQQESATTQKAEKEVTRMVIIMVIFFLICWVPYAYVAFYIFTHQGSEFGPIFMTVPAFFAKSSAIYNPVIYIMLNKQFRNCMITTLCCGKNPFGDDDASSAATSVSTSQVSPA
>GP38_HUMAN
GSPWNGSDGPEGAREPPWPALPPCDERRCSPFPLGALVPVTAVCLCLFVVGVSGNVVTVMLIGRYRDMRTTTNLYLGSMAVSDLLILLGLPFDLYRLWRSRPWVFGPLLCRLSLYVGEGCTYATLLHMTALSVERYLAICRPLARVLVTRRRVRALIAVLWAVALLSAGPFLFLVG---AEAFSRECRPS----PAQLGALRVMLWVTTAYFFLPFLCLSILYGLIGRELWSSPLR--AASGRERGHRQTVRVLLVVVLAFIICWLPFHVGRIIYINTEDSYFSQYFNIVALQLFYLSASINPILYNLISKKYRAAAFKLLLARKSRPRGFHRSRDDTGGDTVGYTET
>GPRM_HUMAN
PILEINMQSESNITVRDDIDDINTNMYQPLSYPLSFQVSLTGFLMLEIVLGLGSNLTVLVLYCMKSNINSVSNIITMNLHVLDVIICVGCIPLTIVILLLSLESNTALICCFHEACVSFASVSTAINVFAITLDRYDISVKPA-NRILTMGRAVMLMISIWIFSFFSFLIPFIEVNFFS--TKTLLCVST---EYYTELGMYYHLLVQIPIFFFTVVVMLITYTKILQALNIRTISQHEARERRERQKRVFRMSLLIISTFLLCWTPISVLNTTILCLGPSDLLVKLRLCFLVMAYGTTIFHPLLYAFTRQKFQKVLKSKMKKRVVSIVEADPLPNWIDPKRNKKITF
>GPR3_HUMAN
GAGSPLAWLSAGSGNVNVSSVGPAEGPTGPAAPLPSPKAWDVVLCISGTLVSCENALVVAIIVGTPAFRAPMFLLVGSLAVADLLAG-LGLVLHFAAV-F--CIGSAEMSLVLVGVLAMAFTASIGSLLAITVDRYLSLYNALYYSETTVTRTYVMLALVWGGALGLGLLPVLAWNCL-D--GLTTCGVV-------YPLSKNHLVVLAIAFFMVFGIMLQLYAQICRIVCRHIALHLLPASHYVATRKGIATLAVVLGAFAACWLPFTVYCLLGDA----HSPPLYTYLTLLPATYNSMINPIIYAFRNQDVQKVLWAVCCCCSSSKIPFRSRSP------------
>CCKR_RAT
DSLLMNGSNITPPCELGLENETLFCLDQPQPSKEWQSALQILLYSIIFLLSVLGNTLVITVLIRNKRMRTVTNIFLLSLAVSDLMLCLFCMPFNLIPNLLKDFIFGSAVCKTTTYFMGTSVSVSTFNLVAISLERYGAICRPLSRVWQTKSHALKVIAATWCLSFTIMTPYPIYN----NNQTANMCRFL---LPSDAMQQSWQTFLLLILFLLPGIVMVVAYGLISLELYQGRLNSSSSAANLIAKKRVIRMLIVIVVLFFLCWMPIFSANAWRAYDTVSHLSGTPISFILLLSYTSSCVNPIIYCFMNKRFRLGFMATFPCCPNPGPPGVRGEVTIRALLSRYSYS
>RDC1_MOUSE
YAEPGNYSDINWPCNSSDCIVVDTVQCPTMPNKNVLLYTLSFIYIFIFVIGMIANSVVVWVNIQAKTTGYDTHCYILNLAIADLWVVITTPVWVVSLVQHNQWPMGELTCKITHLIFSINLFGSIFFLACMSVDRYLSITYFTTSSYKKKMVRRVVCILVWLLAFFVSLPDTYYLKAVTNNE--TYCRSFYPEHSIKEWLIGMELVSVILGFAVPFTIIAIFYFLLARAM---------SASGDQEKHSSRKIIFSYVVVFLVCWLPYHFVVLLDIFSILHNVLFTALHVTQCLSLVHCCVNPVLYSFINRNYRYELMKAFIFKYSAKTGLTKLIDEYSALEQNTK--
>OLF9_RAT
-------------MTRRNQTAISQFFLLGLPFPPEYQHLFYALFLAMYLTTLLGNLIIIILILLDSHLHTPMYLFLSNLSFADLCFSSVTMPKLLQNMQSQVPSIPYAGCLAQIYFFLFFGDLGNFLLVAMAYDRYVAICFPLYMSIMSPKLCVSLVVLSWVLTTFHAMLHTLLMARL---SPHYFCDMSKVACSDTHDNELAIFILGGPIVVLPFLLIIVSYARIVSSI--------FKVPSSQSIHKAFSTCGSHLSVVSLFYGTVIGLYLCPSANN---STVKETVMSLMYTMVTPMLNPFIYSLRNRDIKDALEKIMCKKQIPSFL------------------
>5H1D_RABIT
SPSNQSAEGLPQEAANRSLNATGTPEAWDPGTLQALKISLAVVLSIITVATVLSNTFVLTTILLTRKLHTPANYLIGSLATTDLLVSILVMPISIAYTITHTWNFGQVLCDIWVSSDITCCTASILHLCVIALDRYWAITDALYSKRRTAGHAAAMIAVVWAISICISIPPLFWR----H-E--SDCLVN-------TSQISYTIYSTCGAFYIPSVLLIVLYGRIYMAARNRKLALERKRISAARERKATKTLGIILGAFIGCWLPFFVASLVLPICRDSWMPPGLFDFFTWLGYLNSLINPIIYTVFNEDFRQAFQRVIHFRKAF---------------------
>ML1C_.ENLA
-----MMEVNSTCLDCRTPGTIRTEQDAQDSASQGLTSALAVVLIFTIVVDVLGNILVILSVLRNKKLQNAGNLFVVSLSIADLVVAVYPYPVILIAIFQNGWTLGNIHCQISGFLMGLSVIGSVFNITAIAINRYCYICHSLYDKLYNQRSTWCYLGLTWILTIIAIVPNFFVGS-LQYDP-IFSCTFA------QTVSSSYTITVVVVHFIVPLSVVTFCYLRIWVLVIQVHR-QDFKQKLTQTDLRNFLTMFVVFVLFAVCWAPLNFIGLAVAINPFHKIPEWLFVLSYFMAYFNSCLNAVIYGVLNQNFRKEYKRILMSLLTPRLLFLDTSRSKPSPAVTNNNQ
>SSR2_PIG
LLNGSQPWLSSPFDLNGSVATANSSNQTEPYYDLTSNAVLTFIYFVVCIIGLCGNTLVIYVILRYAKMKTITNIYILNLAIADELFMLGLPFLAMQVALVH-WPFGKAICRVVMTVDGINQFTSIFCLTVMSIDRYLAVVHPISAKWRRPRTAKMINVAVWGVSLLVILPIMIYAGLRSWG--RSSCTIN-WPGESGAWYTGFIIYAFILGFLVPLTIICLCYLFIIIKVKSSGIR-VGSSKRKKSEKKVTRMVSIVVAVFIFCWLPFYIFNVSSVSVAISPALKGMFDFVVVLTYANSCANPILYAFLSDNFKKSFQNVLCLVKVSGTDDGERSDLNETTETQRTLL
>NTR1_HUMAN
EEALLAPGFGNASGNASERVLAAPSSELDVNTDIYSKVLVTAVYLALFVVGTVGNTVTAFTLARKKSLQSTVHYHLGSLALSDLLTLLLAMPVELYNFIWVHWAFGDAGCRGYYFLRDACTYATALNVASLSVERYLAICHPFAKTLMSRSRTKKFISAIWLASALLTVPMLFTMGEQN--SGGLVCTPT---IHTATVKVVIQVN-TFMSFIFPMVVISVLNTIIANKLTVMEHSMAIEPGRVQALRHGVRVLRAVVIAFVVCWLPYHVRRLMFCYISDEDFYHYFYMVTNALFYVSSTINPILYNLVSANFRHIFLATLACLCPVWRRRRK-RPSVSSNHTLSSNA
>AG2R_CAVPO
------MILNSSTEDGIKRIQDDC---PKAGRHSYIFVMIPTLYSIIFVVGIFGNSLVVIVIYFYMKLKTVASVFLLNLALADICFLLTLPLWAVYTAMEYRWPFGNYLCKIASASVSFNLYASVFLLTCLSIDRYLAIVHPMSRLRRTMLVAKVTCVIIWLMAGLASLPAVIHRNVFF-TN--TVCAFH-YESQNSTLPIGLGLTKNILGFMFPFLIILTSYTLIWKALKKA----YEIQKNKPRNDDIFKIIMAIVLFFFFSWVPHQIFTFLDVLIQLGDIVDTAMPITICIAYFNNCLNPLFYGFLGKKFKKYFLQLLKYIPPKAKSHSTLSTRPSD-------N
>PE23_RABIT
PFCTRLNHSYPGMWAP----EARGNLTRPPGPGEDCGSVSVAFPITMLITGFVGNALAMLLVSRSYRRKKSFLLCIGWLALTDLVGQLLTSPVVILVYLSKQLDPSGRLCTFFGLTMTVFGLSSLFIASAMAVERALAIRAPHYASHMKTRATRAVLLGVWLAVLAFALLPVLGVG----QYTGTWCFISNGTSSSHNWGNLFFASTFAFLGLLALAITFTCNLATIKALVSRAKASQSSAQWGRITTETAIQLMGIMCVLSVCWSPLLIMMLKMIFNQTSQKECNFFLIAVRLASLNQILDPWVYLLLRKILLRKFCQVIHENNEQKDEIQRENRHEEARDSEKSKT
>CB1A_FUGRU
TDLFGNRNTTRD-ENSIQCGENFMDMECFMILTPSQQLAVAVLSLTLGTFTVLENLVVLCVIFQSRTRCRPSYHFIGSLAVADLLGSVIFVYSFLDFHVFHK-KDSPNVFLFKLGGVTASFTASVGSLFLTAIDRYISIHRPLYRRIVTRTKAVIAFCMMWTISIIIAVLPLLGWNCK-R--LNSVCSDI------FPLIDENYLMFWIGVTSVLVLFIIYAYIYILWKAHHHAEGQTTRPEQTRMDIRLAKTLVLILAVLVICWGPLLAIMVYDLFWKMDDNIKTVFAFCSMLCLLNSTVNPIIYALRSRDLRHAFLSSCHACRGSAQQLDNSLENVN--ISANRAA
>CKR9_HUMAN
DDYGSESTSSMEDYVNFNFTDFYC---EKNNVRQFASHFLPPLYWLVFIVGALGNSLVILVYWYCTRVKTMTDMFLLNLAIADLLFLVTLPFWAIAAADQ--WKFQTFMCKVVNSMYKMNFYSCVLLIMCISVDRYIAIAQAMTWREKRLLYSKMVCFTIWVLAAALCIPEILYSQIKE-G--IAICTMVYPSDESTKLKSAVLTLKVILGFFLPFVVMACCYTIIIHTL---------IQAKKSSKHKALKVTITVLTVFVLSQFPYNCILLVQTIDAYATNIDICFQVTQTIAFFHSCLNPVLYVFVGERFRRDLVKTLKNLGCISQAQWVSFT--------GSLK
>V1BR_HUMAN
---MDSGPLWDANPTPRGTLSAPNATTPWLGRDEELAKVEIGVLATVLVLATGGNLAVLLTLGQLGRKRSRMHLFVLHLALTDLAVALFQVLPQLLWDITYRFQGPDLLCRAVKYLQVLSMFASTYMLLAMTLDRYLAVCHPLRSLQQPGQSTYLLIAAPWLLAAIFSLPQVFIFSLRESG--VLDCWAD---FGFPWGPRAYLTWTTLAIFVLPVTMLTACYSLICHEICKNRGLVSSINTISRAKIRTVKMTFVIVLAYIACWAPFFSVQMWSVWDKNADSTNVAFTISMLLGNLNSCCNPWIYMGFNSHLLPRPLRHLACCGGPQPRMRRRLSHTTLLTRSSCPA
>5H4_RAT
------------------MDRLDANVSSNEGFGSVEKVVLLTFFAMVILMAILGNLLVMVAVCRDRQRKIKTNYFIVSLAFADLLVSVLVNAFGAIELVQDIWFYGEMFCLVRTSLDVLLTTASIFHLCCISLDRYYAICCQPYRNKMTPLRIALMLGGCWVIPMFISFLPIMQGWNNIRK-NSTFCVFM--------VNKPYAITCSVVAFYIPFLLMVLAYYRIYVTAKEHSRPDQHSTHRMRTETKAAKTLCVIMGCFCFCWAPFFVTNIVDPFIDYT-VPEKVWTAFLWLGYINSGLNPFLYAFLNKSFRRAFLIILCCDDERYKRPPILGQTINGSTHVLR--
>5HTA_DROME
DSSLFGEMLANRSGHLDLINGTTSKVAEDDFTQLLRMAVTSVLLGLMILVTIIGNVFVIAAIILERNLQNVANYLVASLAVADLFVACLVMPLGAVYEISQGWILGPELCDIWTSCDVLCCTASILHLVAIAVDRYWAVTNI-YIHSRTSNRVFMMIFCVWTAAVIVSLAPQFGWKDPDE-Q--QKCMVS--------QDVSYQVFATCCTFYVPLMVILALYWKIYQTARKRPMQKRKETLEAKRERKAAKTLAIITGAFVVCWLPFFVMALTMPLCAACQISDSVASLFLWLGYFNSTLNPVIYTIFSPEFRQAFKRILFGGHRPVHYRSGKL-------------
>HH1R_CAVPO
-MSFLPGMTPVTLSNFSWALEDRMLEGNSTTTPTRQLMPLVVVLSSVSLVTVALNLLVLYAVRSERKLHTVGNLYIVSLSVADLIVGAVVMPMSILYLHRSAWILGRPLCLFWLSMDYVASTASIFSVFILCIDRYRSVQQPLYLRYRTKTRASATILGAWLLSFLWVIPILGWHHFM--EP-EKKCETD------FYDVTWFKVMTAIINFYLPTLLMLWFYIRIYKAVRRHLRSQYTSGLHLNRERKAAKQLGCIMAAFILCWIPYFVFFMVIAFCKSC-SNEPVHMFTIWLGYLNSTLNPLIYPLCNENFRKTFKRILRIPP-----------------------
>OPSD_PHOGR
MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVGLTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGFNFGPIFMTLPAFFAKAAAIYNPVIYIMMNKQFRTCMITTLCCGKNPLGDDEVSASASKTETSQVAPA
>MSHR_CANFA
VWQGPQRRLLGSLNGTSPATPHFELAANQTGPRCLEVSIPNGLFLSLGLVSVVENVLVVAAIAKNRNLHSPMYYFIGCLAVSDLLVSVTNVLETAVMLLVEAAAVVQQLDDIIDVLICGSMVSSLCFLGAIAVDRYLSIFYALYHSIVTLPRAWRAISAIWVASVLSSTLFIAYYY----------------------NNHTAVLLCLVSFFVAMLVLMAVLYVHMLARARQH-IAKRQHSVHQGFGLKGAATLTILLGIFFLCWGPFFLHLSLMVLCPQHGCVFQNFNLFLTLIICNSIIDPFIYAFRSQELRKTLQEVVLCSWA----------------------
>A2AC_HUMAN
AAGPNASGAGERGSGGVANASGASWGPPRGQYSAGAVAGLAAVVGFLIVFTVVGNVLVVIAVLTSRALRAPQNLFLVSLASADILVATLVMPFSLANELMAYWYFGQVWCGVYLALDVLFCTSSIVHLCAISLDRYWSVTQAVYNLKRTPRRVKATIVAVWLISAVISFPPLVSLYR-------PQCGLN--------DETWYILSSCIGSFFAPCLIMGLVYARIYRVAKRRRRAVCRRKVAQAREKRFTFVLAVVMGVFVLCWFPFFFIYSLYGICREAQVPGPLFKFFFWIGYCNSSLNPVIYTVFNQDFRPSFKHILFRRRRRGFRQ-----------------
>B1AR_SHEEP
VPDGAATAARLLVP-SPLRLAADLGQRGTPLLSQQWTVGMGLLMAFIVLLIVAGNVLVIVAIAKTPRLQTLTNLFIMSLASADLVMGLLVVPFGATIVVWGRWEYGSFFCELWTSVDVLCVTASIETLCVIALDRYLAITSPFYQSLLTRARARALVCTVWAISALVSFLPIFMQWWGDS-R--ECCDFI--------INEGYAITSSVVSFYVPLCIMAFVYLRVFREAQKQNGRRRPSRLVALREQKALKTLGIIMGVFTLCWLPFFLANVVKAFHRDL-VPDRLFVFFNWLGYANSAFNPIIYCRSP-DFRKAFQRLLCCARRAACGSHGAAGCLAVARPSPSPG
>OPS4_DROME
S-GNGDLQFLGWNVPPDQIQYIPEHWLTQLEPPASMHYMLGVFYIFLFCASTVGNGMVIWIFSTSKSLRTPSNMFVLNLAVFDLIMCLKAPIFNSFHRGFA-IYLGNTWCQIFASIGSYSGIGAGMTNAAIGYDRYNVITKPM-NRNMTFTKAVIMNIIIWLYCTPWVVLPLTQFWDRFPEGYLTSCSFD--YLSDNFDTRLFVGTIFFFSFVCPTLMILYYYSQIVGHVFSHNVESNVDKSKETAEIRIAKAAITICFLFFVSWTPYGVMSLIGAFGDKSLLTQGATMIPACTCKLVACIDPFVYAISHPRYRLELQKRCPWLGVNEKSGEISSAQQQTTAA-----
>O1E1_HUMAN
-------------MMGQNQTSISDFLLLGLPIQPEQQNLCYALFLAMYLTTLLGNLLIIVLIRLDSHLHTPMYLFLSNLSFSDLCFSSVTIPKLLQNMQNQDPSIPYADCLTQMYFFLLFGDLESFLLVAMAYDRYVAICFPLYTAIMSPMLCLALVALSWVLTTFHAMLHTLLMARL---CPHFFCDMSKLAFSDTRVNEWVIFIMGGLILVIPFLLILGSYARIVSSI--------LKVPSSKGICKAFSTCGSHLSVVSLFYGTVIGLYLCSSANS---STLKDTVMAMMYTVVTPMLNPFIYSLRNRDMKGALSRVIHQKKTFFSL------------------
>MSHR_CHICK
-MSMLAPLRLVREPWNASEGNQSNATAGAGGAWCQGLDIPNELFLTLGLVSLVENLLVVAAILKNRNLHSPTYYFICCLAVSDMLVSVSNLAKTLFMLLMEHASIVRHMDNVIDMLICSSVVSSLSFLGVIAVDRYITIFYALYHSIMTLQRAVVTMASVWLASTVSSTVLITYYY----------------------RRNNAILLCLIGFFLFMLVLMLVLYIHMFALACHHSISQKQPTIYRTSSLKGAVTLTILLGVFFICWGPFFFHLILIVTCPTNTCFFSYFNLFLILIICNSVVDPLIYAFRSQELRRTLREVVLCSW-----------------------
>A1AA_BOVIN
----------MVFLSGNASDSSNCT-HPPPPVNISKAILLGVILGGLILFGVLGNILVILSVACHRHLHSVTHYYIVNLAVADLLLTSTVLPFSAIFEILGYWAFGRVFCNVWAAVDVLCCTASIMGLCIISIDRYIGVSYPLYPTIVTQKRGLMALLCVWALSLVISIGPLFGWRQP--AP--TICQIN--------EEPGYVLFSALGSFYVPLTIILVMYCRVYVVAKREAKNFSVRLLKFSREKKAAKTLGIVVGCFVLCWLPFFLVMPIGSFFPDFRPSETVFKIAFWLGYLNSCINPIIYPCSSQEFKKAFQNVLRIQCLRRKQSSKHTLSH----------
>CKR5_PYGBI
----MDYQVSSPTYDIDYYTSEPC---QKVNVKQIAARLLPPLYSLVFIFGFVGNILVVLILINCKRLKSMTDIYLLNLAISDLFFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWVVAVFASLPGIIFTRSQR-EG-HYTCSSHFPYSQYQFWKNFQTLKIVILGLVLPLLVMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCYIFA-------SSVY
>OPSO_SALSA
TLRIAVNGVSYNEASEIYKPHADPFTGPITNLAPWNFAVLATLMFVITSLSLFENFTVMLATYKFKQLRQPLNYIIVNLSLADFLVSLTGGTISFLTNARGYFFLGNWACVLEGFAVTYFGIVAMWSLAVLSFERYFVICRPLGNVRLRGKHAALGLLFVWTFSFIWTIPPVFGWCSYTVSKIGTTCEPN--WYSNNIWNHTYIITFFVTCFIMPLGMIIYCYGKLLQKLRKVSHD--RLGNAKKPERQVSRMVVVMIVAYLVGWTPYAAFSIIVTACPTIYLDPRLAAAPAFFSKTAAVYNPVIYVFMNKQVSTQLNWGFWSRA-----------------------
>NY1R_RAT
TLFSRVENYSVHYNVSE-NSPFLAFENDDCHLPLAVIFTLALAYGAVIILGVSGNLALIIIILKQKEMRNVTNILIVNLSFSDLLVAVMCLPFTFVYTLMDHWVFGETMCKLNPFVQCVSITVSIFSLVLIAVERHQLIINPR-GWRPNNRHAYIGITVIWVLAVASSLPFVIYQILTDKD--KYVCFDK---FPSDSHRLSYTTLLLVLQYFGPLCFIFICYFKIYIRLKRRNMMIRDSKYRSSETKRINVMLLSIVVAFAVCWLPLTIFNTVFDWNHQICNHNLLFLLCHLTAMISTCVNPIFYGFLNKNFQRDLQFFFNFCDFRSRDDDY---TMHTDVSKTSLK
>OLF3_CANFA
-------------MGTGNQTWVREFVLLGLSSDWDTEVSLFVLFLITYMVTVLGNFLIILLIRLDSRLHTPMYFFLTNLSLVDVSYATSIIPQMLAHLLAAHKAIPFVSCAAQLFFSLGLGGIEFVLLAVMAYDRYVAVCDPLYSVIMHGGLCTRLAITSWVSGSMNSLMQTVITFQL---PDHISCELLRLACVDTSSNEIAIMVSSIVLLMTPFCLVLLSYIQIISTI--------LKIQSTEGRKKAFHTCASHLTVVVLCYGMAIFTYIQPRSSP---SVLQEKLISLFYSVLTPMLNPMIYSVRNKEVKGAWQKLLGQLTGITSKLAT---------------
>P2Y9_HUMAN
DRRFIDFQFQDSNSSLRPRLGNATANNTCIVDDSFKYNLNGAVYSVVFILGLITNSVSLFVFCFRMKMRSETAIFITNLAVSDLLFVCTLPFKIFYNFNRH-WPFGDTLCKISGTAFLTNIYGSMLFLTCISVDRFLAIVYPFSRTIRTRRNSAIVCAGVWILVLSGGISASLFSTTN----ATTTCFEGFSKRVWKTYLSKITIFIEVVGFIIPLILNVSCSSVVLRTLRKP----ATLSQIGTNKKKVLKMITVHMAVFVVCFVPYNSVLFLYALVRSQRFAKIMYPITLCLATLNCCFDPFIYYFTLESFQKSFYINAHIRMESLFKTETPLTIQEEVSDQTTNN
>MC5R_BOVIN
-MNSSFHLHFLDLGLNTTDGNLSGLSVQNASSLCEDMGIAVEVFLALGLISLLENILVIGAIVRNRNLHTPMYFFVGSLAVADMLVSLSNSWETITIYLLTNDASVRHLDNVFDSMICISVVASMCSLLAIAVDRYVTIFCALYQRIMTGRRSGAIIGGIWAFCASCGTVFIVYYY----------------------EESTYVVICLIAMFLTMLLLMASLYTHMFLLARTH-RIPGHSSVRQRTGVKGAITLAMLLGVFIVCWAPFFLHLILMISCPHNSCFMSHFNMYLILIMCNSVIDPLIYAFRSQEMRKTFKEIVCFQSFRTPCRFPSRY------------
>CKR5_TRAPH
----MDYQVSSPTYDIDYYTSEPC---QKVNVKQIAARLLPPLYSLVFIFGFVGNILVVLILINCKRLKSMTDIYLLNLAISDLFFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWVVAVFASLPGIIFTRSQR-EG-HYTCSSHFPYSQYQFWKNFQTLKIVILGLVLPLLVMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCSIFA-------SSVY
>OPSB_SAIBB
--MSKMPEEEEFYLFKNISSVGPWDGPQYHIAPVWAFQLQAAFMGIVFLAGLPLNSMVLVATVRYKKLRHPLNYVLVNVSVGGFLLCIFSVLPVFVNSCNGYFVFGRHVCALEGFLGTVAGLVTGWSLAFLAFERYIVICKPFGNFRFSSKHALMVVLTTWTIGIGVSIPPFFGWSRYIAEGLQCSCGPDWYTVGTKYRSEYYTWFLFIFCFIVPLSLICFSYAQLLRALKAVAAQQQESATTQKAEREVSRMVVVMVGSFCVCYVPYAALAMYMVNNRNHGLDLRLVSIPAFFSKSSCIYNPIIYCFMNKQFRACIMEMVCGKAMTD---ESDISSTVSSSQVGPN-
>OAR_DROME
LSTAQADKDSAGECEGAVEELHASILGLQLAVPEWEALLTALVLSVIIVLTIIGNILVILSVFTYKPLRIVQNFFIVSLAVADLTVALLVLPFNVAYSILGRWEFGIHLCKLWLTCDVLCCTSSILNLCAIALDRYWAITDPIYAQKRTVGRVLLLISGVWLLSLLISSPPLIGWNDW-T-S--TPCELT--------SQRGYVIYSSLGSFFIPLAIMTIVYIEIFVATRRRGVNEEKQKISLSKERRAARTLGIIMGVFVICWLPFFLMYVILPFCQTCCPTNKFKNFITWLGYINSGLNPVIYTIFNLDYRRAFKRLLGLN------------------------
>PAR3_HUMAN
-WTGATITVKIKCPEESASHLHVKNATMGYLTSSLSTKLIPAIYLLVFVVGVPANAVTLWMLFFRTR-SICTTVFYTNLAIADFLFCVTLPFKIAYHLNGNNWVFGEVLCRATTVIFYGNMYCSILLLACISINRYLAIVHPFYRGLPKHTYALVTCGLVWATVFLYMLPFFILKQEYY--D--TTCHDVNTCESSSPFQLYYFISLAFFGFLIPFVLIIYCYAAIIRTL----------NAYDHRWLWYVKASLLILVIFTICFAPSNIILIIHHANYYYDGLYFIYLIALCLGSLNSCLDPFLYFLMSKTRNHSTAYLTK--------------------------
>NMBR_HUMAN
SNLSVTTGANESGSVPEGWERDFLPASDGTTTELVIRCVIPSLYLLIITVGLLGNIMLVKIFITNSAMRSVPNIFISNLAAGDLLLLLTCVPVDASRYFFDEWMFGKVGCKLIPVIQLTSVGVSVFTLTALSADRYRAIVNPMMQTSGALLRTCVKAMGIWVVSVLLAVPEAVFSEVARSS--FTACIPY---QTDELHPKIHSVLIFLVYFLIPLAIISIYYYHIAKTLIKSLPGNEHTKKQMETRKRLAKIVLVFVGCFIFCWFPNHILYMYRSFNYNELGHMIVTLVARVLSFGNSCVNPFALYLLSESFRRHFNSQLCCGRKSYQERGTSYLMTSLKSNAKNMV
>AA1R_BOVIN
----------------------------MPPSISAFQAAYIGIEVLIALVSVPGNVLVIWAVKVNQALRDATFCFIVSLAVADVAVGALVIPLAILINIG--PRTYFHTCLKVACPVLILTQSSILALLAMAVDRYLRVKIPLYKTVVTPRRAVVAITGCWILSFVVGLTPMFGWNNLSN-G-VIECQFE-----KVISMEYMVYFNFFVWVLPPLLLMVLIYMEVFYLIRKQK-VSGDPQKYYGKELKIAKSLALILFLFALSWLPLHILNCITLFCPSCHMPRILIYIAIFLSHGNSAMNPIVYAFRIQKFRVTFLKIWNDHFRCQPAPPID--PDD---------
>C.C1_MOUSE
----------MESSTAFYDYHDKLSLLCENNVIFFSTISTIVLYSLVFLLSLVGNSLVLWVLVKYENLESLTNIFILNLCLSDLMFSCLLPVLISAQWS---WFLGDFFCKFFNMIFGISLYSSIFFLTIMTIHRYLSVVSPITLGIHTLRCRVLVTSCVWAASILFSIPDAVFHKVIS-----LNCKYS------EHHGFLASVYQHNIFFLLSMGIILFCYVQILRTL---------FRTRSRQRHRTVRLIFTVVVAYFLSWAPYNLTLFLKTGIIQQQQLDIAMIICRHLAFSHCCFNPVLYVFVGIKFRRHLKHLFQQVWLCRKTSSTVPCEGPSFY------
>OPS1_PATYE
NGTLNRSMTPNTGWEGPYDMSVHLHWTQFPPVTEEWHYIIGVYITIVGLLGIMGNTTVVYIFSNTKSLRSPSNLFVVNLAVSDLIFSAVNGFPLLTVSSFHQWIFGSLFCQLYGFVGGVFGLMSINTLTAISIDRYVVITKPLASQTMTRRKVHLMIVIVWVLSILLSIPPFFGWGAYIPEGFQTSCTFD--YLTKTARTRTYIVVLYLFGFLIPLIIIGVCYVLIIRGVRRHSMKARANNKRARSELRISKIAMTVTCLFIISWSPYAIIALIAQFGPAHWITPLVSELPMMLAKSSSMHNPVVYALSHPKFRKALYQRVPWLFCCCKPKEKADFRSVTRTESVNSD
>5H2B_MOUSE
ILQKTCDHLILTNRSGLETDSVAEEMKQTVEGQGHTVHWAALLILAVIIPTIGGNILVILAVALEKRLQYATNYFLMSLAIADLLVGLFVMPIALLTIMFEAWPLPLALCPAWLFLDVLFSTASIMHLCAISLDRYIAIKKPIANQCNTRATAFIKITVVWLISIGIAIPVPIKGIET-NP---VTCELT------KDRFGSFMVFGSLAAFFVPLTIMVVTYFLTIHTLQKKRRMGKRSAQTISNEQRASKALGVVFFLFLLMWCPFFITNLTLALCDSCTTLKTLLEIFVWIGYVSSGVNPLIYTLFNKTFREAFGRYITCNYRATKSVKALRKGNSMVENSKFFT
>OPSG_HUMAN
DSYEDSTQSSIFTYTNSNSTRGPFEGPNYHIAPRWVYHLTSVWMIFVVIASVFTNGLVLAATMKFKKLRHPLNWILVNLAVADLAETVIASTISVVNQVYGYFVLGHPMCVLEGYTVSLCGITGLWSLAIISWERWMVVCKPFGNVRFDAKLAIVGIAFSWIWAAVWTAPPIFGWSRYWPHGLKTSCGPDVFSGSSYPGVQSYMIVLMVTCCITPLSIIVLCYLQVWLAIRAVAKQQKESESTQKAEKEVTRMVVVMVLAFCFCWGPYAFFACFAAANPGYPFHPLMAALPAFFAKSATIYNPVIYVFMNRQFRNCILQLFGKKVDDG-----SELVSSV--SSVSPA
>ACM4_MOUSE
------MANFTPVNGSSANQSVRLVTTAHNHLETVEMVFIATVTGSLSLVTVVGNILVMLSIKVNRQLQTVNNYFLFSLACADLIIGAFSMNLYTLYIIKGYWPLGAVVCDLWLALDYVVSNASVMNLLIISFDRYFCVTKPLYPARRTTKMAGLMIAAAWVLSFVLWAPAILFWQFVV-VP-DNQCFIQ------FLSNPAVTFGTAIAAFYLPVVIMTVLYIHISLASRSRSIAVRKKRQMAARERKVTRTIFAILLAFILTWTPYNVMVLVNTFCQSC-IPERVWSIGYWLCYVNSTINPACYALCNATFKKTFRHLLLCQYRNIGTAR----------------
>AA1R_RAT
----------------------------MPPYISAFQAAYIGIEVLIALVSVPGNVLVIWAVKVNQALRDATFCFIVSLAVADVAVGALVIPLAILINIG--PQTYFHTCLMVACPVLILTQSSILALLAIAVDRYLRVKIPLYKTVVTQRRAAVAIAGCWILSLVVGLTPMFGWNNLSN-G-VIKCEFE-----KVISMEYMVYFNFFVWVLPPLLLMVLIYLEVFYLIRKQK-VSGDPQKYYGKELKIAKSLALILFLFALSWLPLHILNCITLFCPTCQKPSILIYIAIFLTHGNSAMNPIVYAFRIHKFRVTFLKIWNDHFRCQPKPPID--AED---------
>AA1R_CHICK
----------------------------MAQSVTAFQAAYISIEVLIALVSVPGNILVIWAVKMNQALRDATFCFIVSLAVADVAVGALVIPLAIIINIG--PQTEFYSCLMMACPVLILTESSILALLAIAVDRYLRVKIPVYKSVVTPRRAAVAIACCWIVSFLVGLTPMFGWNNLNN-V-VIKCQFE-----TVISMEYMVYFNFFVWVLPPLLLMLLIYLEVFNLIRTQK-VSNDPQKYYGKELKIAKSLALVLFLFALSWLPLHILNCITLFCPSCKTPHILTYIAIFLTHGNSAMNPIVYAFRIKKFRTAFLQIWNQYFCCKTNKSSS--VN----------
>P2YR_MELGA
PELLAG-----------GWAAGNASTKCSLTKTGFQFYYLPTVYILVFITGFLGNSVAIWMFVFHMRPWSGISVYMFNLALADFLYVLTLPALIFYYFNKTDWIFGDVMCKLQRFIFHVNLYGSILFLTCISVHRYTGVVHPLSLGRLKKKNAVYVSSLVWALVVAVIAPILFYSGTG----KTITCYDT-TADEYLRSYFVYSMCTTVFMFCIPFIVILGCYGLIVKALIY------KDLDNSPLRRKSIYLVIIVLTVFAVSYLPFHVMKTLNLRARLDDKVYATYQVTRGLASLNSCVDPILYFLAGDTFRRRLSRATRKSSRRSEPNVQSKSLTEYKQNGDTSL
>GRE1_BALAM
----MEGPPLSPAPADNVTLNVSCGRPATLFDWADHRLISLLALAFLNLMVVAGNLLVVMAVFVHSKLRTVTNLFIVSLACADLLVGMLVLPFSATLEVLDVWLYGDVWCSVWLAVDVWMCTSSILNLCAISLDRYLAVSQPIYPSLMSTRRAKQLIAAVWVLSFVICFPPLVGWNDR-G-S--LTCELT--------NERGYVIYSALGSFFLPSTVMLFFYGRIYRTAVSTRVSVRHQARRFRMETKAAKTVGIIVGLFILCWLPFFVCYLVRGFCADC-VPPLLFSVFFWLGYCNSAVNPCVYALCSRDFRFAFSSILCKCVCRRGAMERRFRRSQTEEDCEVAD
>TDA8_MOUSE
---------------------MAMNSMCIEEQRHLEHYLFPVVYIIVFIVSVPANIGSLCVSFLQAKKENELGIYLFSLSLSDLLYALTLPLWINYTWNKDNWTFSPTLCKGSVFFTYMNFYSSTAFLTCIALDRYLAVVYPLFSFLRTRRFAFITSLSIWILESFFNSMLLWKDETS-DKSNFTLCYDK---YPLEKWQINLNLFRTCMGYAIPLITIMICNHKVYRAV------RHNQATENSEKRRIIKLLASITLTFVLCFTPFHVMVLIRCVLERDWQTFTVYRVTVALTSLNCVADPILYCFVTETGRADMWNILKLCTRKHNRHQGKKRRDAVELEIID--
>MSHR_CEREL
PVLGSQRRLLGSLNCTPPATFPLTLAPNRTGPQCLEVAIPDGLFLSLGLVSLVENVLVVAAIAKNRNLQSPMYYFICCLAMSDLLVSVSNVLETAVMLLLEAAAVVQQLDNVIDVLICGSMVSSLCFLGAIAVDRYISIFYALYHSVVTLPRAWRIIAAIWVASILTSLLFITYYY----------------------NNHTVVLLCLVGFFIAMLALMAVLYVHMLARACQHIARKRQRPIHQGFGLKGAATLTILLGVFFLCWGPFFLHLSLIVLCPQHGCIFKNFNLFLALIICNAIVDPLIYAFRSQELRKTLQEVLQCSW-----------------------
>NY1R_MOUSE
TLFSKVENHSIHYNASE-NSPLLAFENDDCHLPLAVIFTLALAYGAVIILGVSGNLALIIIILKQKEMRNVTNILIVNLSFSDLLVAVMCLPFTFVYTLMDHWVFGETMCKLNPFVQCVSITVSIFSLVLIAVERHQLIINPR-GWRPNNRHAYIGITVIWVLAVASSLPFVIYQILTDKD--KYVCFDK---FPSDSHRLSYTTLLLVLQYFGPLCFIFICYFKIYIRLKRRNMMIRDSKYRSSETKRINIMLLSIVVAFAVCWLPLTIFNTVFDWNHQICNHNLLFLLCHLTAMISTCVNPIFYGFLNKNFQRDLQFFFNFCDFRSRDDDY---TMHTDVSKTSLK
>GPR1_MACMU
EDLEETLFEEFENYSYALDYYSLESDLEEKVQLGVVHWVSLVLYCLSFVLGIPGNAIVIWFTGFKWKKTVS-TLWFLNLAIADFIFLLFLPLYISYVVMNFHWPFGIWLCKANSFTAQLNMFASVFFLTVISLDHYIHLIHPVSHRHRTLKNSLIVIIFIWLLASLIGGPALYFR--D--NN-HTLCYNNHDPDLTVIRHHVLTWVKFIVGYLFPLLTMSICYLCLIFKV---------KKRSILISSRHFWTILAVVVAFVVCWTPYHLFSIWELTIHHNHVMQAGIPLSTGLAFLNSCLNPILYVLISKKFQARFRSSVAEILKYTLWEVSCSGNSETKNLCLLET
>ACM1_PIG
------------MNTSAPPAVSPNITVLAPGKGPWQVAFIGITTGLLSLATVTGNLLVLISFKVNTELKTVNNYFLLSLACADLIIGTFSMNLYTTYLLMGHWALGTLACDLWLALDYVASNASVMNLLLISFDRYFSVTRPLYRAKRTPRRAALMIGLAWLVSFVLWAPAILFWQYLV-VL-AGQCYIQ------FLSQPIITFGTAMAAFYLPVTVMCTLYWRIYRETENRRGKAKRKTFSLVKEKKAARTLSAILLAFIVTWTPYNIMVLVSTFCKDC-VPETLWELGYWLCYVNSTINPMCYALCNKAFRDTFRLLLLCRWDKRRWRKIPKRPSRQC-------
>RTA_RAT
EAHSTNQNKMCPGMSEALELYSRGFLTIEQIATLPPPAVTNYIFLLLCLCGLVGNGLVLWFFGFSIK-RTPFSIYFLHLASADGIYLFSKAVIALLNMGTFLGSFPDYVRRVSRIVGLCTFFAGVSLLPAISIERCVSVIFPMYWRRRPKRLSAGVCALLWLLSFLVTSIHNYFCM-FL---SGTACL-----------NMDISLGILLFFLFCPLMVLPCLALILHVE---------CRARRRQRSAKLNHVVLAIVSVFLVSSIYLGIDWFLFWVFQ--IPAPFPEYVTDLCICINSSAKPIVYFLAGRDKSQRLWEPLRVVFQRALRDGAEPGNTVTMEMQCPSG
>GP43_HUMAN
------------------------------MLPDWKSSLILMAYIIIFLTGLPANLLALRAFVGRIRQPAPVHILLLSLTLADLLLLLLLPFKIIEAASNFRWYLPKVVCALTSFGFYSSIYCSTWLLAGISIERYLGVAFPVYKLSRRPLYGVIAALVAWVMSFGHCTIVIIVQYLNT--N--ITCYEN-FTDNQLDVVLPVRLELCLVLFFIPMAVTIFCYWRFVWIMLSQP------LVGAQRRRRAVGLAVVTLLNFLVCFGPYNVSHLVGYHQR---KSPWWRSIAVVFSSLNASLDPLLFYFSSSVVRRAFGRGLQVLRNQGSSLLGRRG---------TAE
>CKR1_HUMAN
METPNTTEDYDTTTEFDYGDATPC---QKVNERAFGAQLLPPLYSLVFVIGLVGNILVVLVLVQYKRLKNMTSIYLLNLAISDLLFLFTLPFWIDYKLKDD-WVFGDAMCKILSGFYYTGLYSEIFFIILLTIDRYLAIVHAVALRARTVTFGVITSIIIWALAILASMPGLYFSKTQW-EF-HHTCSLHFPHESLREWKLFQALKLNLFGLVLPLLVMIICYTGIIKIL---------LRRPNEKKSKAVRLIFVIMIIFFLFWTPYNLTILISVFQDFLRHLDLAVQVTEVIAYTHCCVNPVIYAFVGERFRKYLRQLFHRRVAVHLVKWLPFLV-------SSTS
>VK02_SPVKA
YEYSTITDYYNTINNDITSSSVIKAFDNNCTFLEDTKYHIIVIHIILFLLGSIGNIFVVSLIAFKRN-KSITDIYILNLSMSDCIFVFQIPFIVYSKLDQ--WIFGNILCKIMSVLYYVGFFSNMFIITLMSIDRYFAIVHPIRQPYRTKRIGILMCCSAWLLSLILSSPVSKLYENIP--MDIYQCTLTENDSIIAFIKRLMQIEITILGFLIPIIIFVYCYYRIFTTV---------VRLRNRRKYKSIKIVLMIVVCSLICWIPLYIVLMIATIVSLYLNLAYAITFSETISLARCCINPIIYTLIGEHVRSRISSICSCIYRDNRIRKKLFSNII---------
>AA1R_CAVPO
----------------------------MPHSVSAFQAAYIGIEVLIALVSVPGNVLVIWAVKVNQALRDATFCFIASLAVADVAVGALVIPLAILINIG--PQTYFHTCLMVACPVLILTQSSILALLAIAVDRYLRVKIPLYKTVVTPRRAAVAIAGCWILSLVVGLTPMFGWNNLSN-G-VIKCEFE-----KVISMEYMVYFNFFVWVLPPLLLMVLIYLEVFYLIRKQK-VSGDPQKYYGKELKIAKSLALILFLFALSWLPLHILNCITLFCPTCHKPTILTYIAIFLTHGNSAMNPIVYAFRIQKFRVTFLKIWNDHFRCQPEPPID--VDD---------
>P2Y6_RAT
-----------MERDNGTIQAPGLPPTTCVYREDFKRLLLPPVYSVVLVVGLPLNVCVIAQICASRRTLTRSAVYTLNLALADLLYACSLPLLIYNYARGDHWPFGDLACRLVRFLFYANLHGSILFLTCISFQRYLGICHPLWHKRGGRRAAWVVCGVVWLVVTAQCLPTAVFAATG-----RTVCYDL-SPPILSTRYLPYGMALTVIGFLLPFTALLACYCRMARRLCRQ---GPAGPVAQERRSKAARMAVVVAAVFVISFLPFHITKTAYLAVRSTETFAAAYKGTRPFASANSVLDPILFYFTQQKFRRQPHDLLQKLTAKWQRQRV---------------
>O2H3_HUMAN
---------------MDNQSSTPGFLLLGFSEHPGLGRTLFVDVITSYLLTLVGNTLIILLSALDTKLHSPMYFFLSNLSFLDLCFTTSCVPQMLANLWGPKKTISFLDCSVQIFIFLSLGTTECILMKVMAFDRYVAVCQPLYATIIHPRLCWQLASVAWVIGLVGSVVQTPSTLHL---PDDFVCEVPRLSCEDTSYNEIQVAVASVFILVVPLSLILVSYGAITWAV--------LRINSATAWRKAFGTCSSHLTVVTLFYSSVIAVYLQPKNPY---AQGRGKFFGLFYAVGTPSLNPLVYTLRNKEIKRALRRLLGKERDSRESWRAA--------------
>OPRK_MOUSE
LPNSSSWFPNWAESDSNGSVGSEDQQLESAHISPAIPVIITAVYSVVFVVGLVGNSLVMFVIIRYTKMKTATNIYIFNLALADALVTTTMPFQSAVYLMNS-WPFGDVLCKIVISIDYYNMFTSIFTLTMMSVDRYIAVCHPVALDFRTPLKAKIINICIWLLASSVGISAIVLGGTKVDVD-VIECSLQFPDDEYSWWDLFMKICVFVFAFVIPVLIIIVCYTLMILRLKSVRLL-SGSREKDRNLRRITKLVLVVVAVFIICWTPIHIFILVEALGSTSTAALSSYYFCIALGYTNSSLNPVLYAFLDENFKRCFRDFCFPIKMRMERQSTNRVASMRDVGGMNKP
>5HTB_DROME
TTSNLSQIVWNRSVNGNGNSNDEQERAAVEFWLLVKMIAMAVVLGLMILVTIIGNVFVIAAIILERNLQNVANYLVASLAVADLFVACLVMPLGAVYEISNGWILGPELCDIWTSCDVLCCTASILHLVAIAADRYWTVTNI-YNNLRTPRRVFLMIFCVWFAALIVSLAPQFGWKDPDE-E--QHCMVS--------QDVGYQIFATCCTFYVPLLVILFLYWKIYIIARKRPHQKRRQLLEAKRERKAAQTLAIITGAFVICWLPFFVMALTMSLCKECEIHTAVASLFLWLGYFNSTLNPVIYTIFNPEFRRAFKRILFGRKAAARARSAKI-------------
>OPR._HUMAN
GSHLQGNLSLLSPNHSLLPPHLLLNASHGAFLPLGLKVTIVGLYLAVCVGGLLGNCLVMYVILRHTKMKTATNIYIFNLALADTLVLLTLPFQGTDILLGF-WPFGNALCKTVIAIDYYNMFTSTFTLTAMSVDRYVAICHPIALDVRTSSKAQAVNVAIWALASVVGVPVAIMGSAQVEE---IECLVE-IPTPQDYWGPVFAICIFLFSFIVPVLVISVCYSLMIRRLRGVRLL-SGSREKDRNLRRITRLVLVVVAVFVGCWTPVQVFVLAQGLGVQPETAVAILRFCTALGYVNSCLNPILYAFLDENFKACFRKFCCASALRRDVQVSDRVALACKTSETVPR
>VG74_KSHV
LDDDESWNETLNMSGYDYSGNFSLEVSVCEMTTVVPYTWNVGILSLIFLINVLGNGLVTYIFCKHRS-RAGAIDILLLGICLNSLCLSISLLAEVLM--FLFNIISTGLCRLEIFFYYLYVYLDIFSVVCVSLVRYLLVAYSTSWPKKQSLGWVLTSAALLIALVLSGDACRHRSRVVD-PVKQAMCYEN---NMTADWRLHVRTVSVTAGFLLPLALLILFYALTWCVV---------RRTKLQARRKVRGVIVAVVLLFFVFCFPYHVLNLLDTLLRRRGLINVGLAVTSLLQALYSAVVPLIYSCLGSLFRQRMYGLFQSLRQSFMSGATT--------------
>PE22_MOUSE
---------------MDNFLNDSKLMEDCKSRQWLLSGESPAISSVMFSAGVLGNLIALALLARRWRSISLFHVLVTELVLTDLLGTCLISPVVLASYSRNQLAPESHACTYFAFTMTFFSLATMLMLFAMALERYLSIGYPYYRRHLSRRGGLAVLPVIYGASLLFCSLPLLNYGEYVQYCPGTWCFIR--------HGRTAYLQLYATMLLLLIVAVLACNISVILNLIRMRGPRRGERTSMAEETDHLILLAIMTITFAICSLPFTIFAYMDETSS---LKEKWDLRALRFLSVNSIIDPWVFAILRPPVLRLMRSVLCCRTSLRTQEAQQTSSKQTDLCGQL--
>A1AD_RABIT
GSGEDNRSSAGEPGGAGGGGEVNGTAAVGGLVVSAQSVGVGVFLAAFILTAVAGNLLVILSVACNRHLQTVTNYFIVNLAVADLLLSATVLPFSATMEVLGFWAFGRAFCDVWAAVDVLCCTASILSLCTISVDRYVGVRHSLYPAIMTERKAAAILALLWAVALVVSMGPLLGWKEP--VP--RFCGIT--------EEVGYAVFSSLCSFYLPMAVIVVMYCRVYVVARSTHTFLSVRLLKFSREKKAAKTLAIVVGVFVLCWFPFFFVLPLGSLFPQLKPSEGVFKVIFWLGYFNSCVNPLIYPCSSREFKRAFLRLLRCQCRRRRRRRPLWRASAGGGPHPDCA
>OPSD_CRIGR
MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVICKPMSNFRFGENHAIMGVVFTWIMALACAAPPLVGWSRYIPEGMQCSCGVDYYTLKPEVNNESFVIYMFVVHFTIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVILMVVFFLICWFPYAGVAFYIFTHQGSNFGPIFMTLPAFFAKSSSIYNPVIYIMMNKQFRNCMLTTLCCGKNILGDDEASATASKTETSQVAPA
>AA1R_CANFA
----------------------------MPPAISAFQAAYIGIEVLIALVSVPGNVLVIWAVKVNQALRDATFCFIVSLAVADVAVGALVIPLAILINIG--PRTYFHTCLMVACPVLILTQSSILALLAIAVDRYLRVKIPLYKTVVTPRRAAVAIAGCWILSFVVGLTPLFGWNRLGN-G-VIKCEFE-----KVISMEYMVYFNFFVWVLPPLLLMVLIYLEVFYLIRRQK-VSGDPQKYYGKELKIAKSLALILFLFALSWLPLHILNCITLFCPSCRKPSILMYIAIFLTHGNSAMNPIVYAFRIQKFRVTFLKIWNDHFRCQPTPPVD--PHD---------
>DCDR_.ENLA
---------MENFSIFNVTVNVWHADLDVGNSDLSLRALTGLLLSLLILSTLLGNTLVCLAVIKFRHRSKVTNFFVISLAVSDLFVALLVMPWKAVTEVAGFWVFG-DFCDTWVAFDIMCSTASILNLCIISLDRYWAIASPFYERKMTQRVAFIMIGVAWTLSILISFIPVQLSWHKSEE-HTENCDSS--------LNRTYAISSSLISFYIPVVIMIGTYTRIYRIAQTQSSRENSLKTSFRKETKVLKTLSIIMGVFVFCWLPFFVLNCMIPFCHMNCVSETTFNIFVWFGWANSSLNPVIYAFNA-DFRKAFTTILGCNRFCSSNNVEAVNYHHDTTFQK---
>O.YR_BOVIN
GAFAANWSAEAVNGSAAPPGTEGNRTAGPPQRNEALARVEVAVLCLILFLALSGNACVLLALRTTRHKHSRLFFFMKHLSIADLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPL--RSLSRRTDRLAVLVTWLGCLVASAPQVHIFSLRE----VFDCWAV---FIQPWGPKAYITWITLAVYIVPVIVLATCYGLISFKIWQNRAIVSNVKLISKAKIRTVKMTFIVVLAFIVCWTPFFFVQMWSVWDADA-KEASPFIIAMLLASLNSCCNPWIYMLFTGHLFQELVQRFLCCSFRRLKGSRPGENSSTFVLSQYSS
>YQNJ_CAEEL
RDSVINASSAVSTTTLPPLDIPMTSMKPPSIIPTVELVLGTITYLVIIAMTVVGNTLVVVAVFSYRPLKKVQNYFLVSLAASDLAVAIFVMPLHVVTFLAGQWLLGVTVCQFFTTADILLCTSSILNLCAIALDRYWAIHNPIYAQKRTTKFVCIVIVIVWILSMLISVPPIIGWNNW-M---EDSCGLS--------TEKAFVVFSAAGSFFLPLLVMVVVYVKIFISARQRNPTRKREKISVAKEKRAAKTIAVIIFVFSFCWLPFFVAYVIRPFCETCTLVQQVEQAFTWLGYINSSLNPFLYGILNLEFRRAFKKILCPKAVLEQRRRRMSA------------
>OPS2_PATYE
------------------MPFPLNRTDTALVISPSEFRIIGIFISICCIIGVLGNLLIIIVFAKRRSVRRPINFFVLNLAVSDLIVALLGYPMTAASAFSNRWIFDNIGCKIYAFLCFNSGVISIMTHAALSFCRYIIICQYGYRKKITQTTVLRTLFSIWSFAMFWTLSPLFGWSSYVIEVVPVSCSVN--WYGHGLGDVSYTISVIVAVYVFPLSIIVFSYGMILQEKVCGIRARYTPRFIQDIEQRVTFISFLMMAAFMVAWTPYAIMSALAIGSFN--VENSFAALPTLFAKASCAYNPFIYAFTNANFRDTVVEIMAPWTTRRVGVSTLPWRRRTSAVNTTDI
>RGR_HUMAN
---------------------MAETSALPTGFGELEVLAVGMVLLVEALSGLSLNTLTIFSFCKTPELRTPCHLLVLSLALADSGIS-LNALVAATSSLLRRWPYGSDGCQAHGFQGFVTALASICSSAAIAWGRYHHYCTRS---QLAWNSAVSLVLFVWLSSAFWAALPLLGWGHYDYEPLGTCCTLD--YSKGDRNFTSFLFTMSFFNFAMPLFITITSYSLME--------------QKLGKSGHLQVNTTLPARTLLLGWGPYAILYLYAVIADVTSISPKLQMVPALIAKMVPTINAINYALGNEMVCRGIWQCLSPQKREKDRTK----------------
>PE23_PIG
------------MWAPERSAEEQGNLTRSLGSSEDCGSVSVVFPMTMLITGFVGNALAMLLVSQSYRRKKSFLLCIGWLALTDMVGQLLTSPVVIVLYLSHQLDPSGRLCTFFGLTMTAFGLSSLFIASAMAVERALAIRAPHYSSHMKTSATRAVLLGVWLAVLAFALLPVLGVG----QYTGTWCFISNETSSENNWGNIFFASAFSFLGLSALVVTFACNLATIKALVSRAKASQSSAQWGRITTETAIQLMGIMCVLSVCWSPLLIMMLKTIFNQTSQNECNFFLIAVRLASLNQILDPWVYLLLRKILLQKFCQAVSQKQREEAATLIFTHPGEARVLFSKSK
>ACM1_MACMU
------------MNTSAPPAVSPNITVLAPGKGPWQVAFIGITTGLLSLATVTGNLLVLISFKVNTELKTVNNYFLLSLACADLIIGTFSMNLYTTYLLMGHWALGTLACDLWLALDYVASNASVMNLLLISFDRYFSVTRPLYRAKRTPRRAALMIGLAWLVSFVLWAPAILFWQYLV-VL--GQCYIQ------FLSQPIITFGTAMAAFYLPVTVMCTLYWRIYRETENRRGKAKRKTFSLVKEKKAARTLSAILLAFILTWTPYNIMVLVSTFCKDC-VPETLWELGYWLCYVNSTINPMCYALCNKAFRDTFRLLLLCRWDKRRWRKIPKRPSRQC-------
>GPRF_CERAE
----MDPEETSVYLDYYYATSPNPDIRETHSHVPYTSVFLPVFYIAVFLTGVLGNLVLMGALHFKPGSRRLIDIFIINLAASDFIFLVTLPLWVDKEASLGLWRTGSFLCKGSSYMISVNMHCSVFLLTCMSVDRYLAIVCPVSRKFRRTDCAYVVCASIWFISCLLGLPTLLSRELT-IDD-KPYCAEK----KATPLKLIWSLVALIFTFFVPLLSIVTCYCRIARKLCAH---YQQSGKHNKKLKKSIKIIFIVVAAFLVSWLPFNTSKLLAIVSGLQAILQLGMEVSGPLAFANSCVNPFIYYIFDSYIRRAIVHCLCPCLKNYDFGSSTETALSTFIHAEDFT
>NK3R_RAT
-G--NFSSALGLPATTQAPSQVRANLTNQFVQPSWRIALWSLAYGLVVAVAVFGNLIVIWIILAHKRMRTVTNYFLVNLAFSDASVAAFNTLINFIYGLHSEWYFGANYCRFQNFFPITAVFASIYSMTAIAVDRYMAIIDPL-KPRLSATATKIVIGSIWILAFLLAFPQCLYSK---GR---TLCYVQ--WPEGPKQHFTYHIIVIILVYCFPLLIMGVTYTIVGITLWGG--PCDKYHEQLKAKRKVVKMMIIVVVTFAICWLPYHVYFILTAIYQQLKYIQQVYLASFWLAMSSTMYNPIIYCCLNKRFRAGFKRAFRWCPFIQVSSY----TTRFHPTRQSSL
>CKR3_CERAE
MTTSLYTVETFGPTSYDDDMGLLC---EKADVGALIAQFVPPLYSLVFTVGLLGNVVVVMILIKYRRLRIMTNIYLLNLAISDLLFLFTLPFWIHYVREHN-WVFSHGMCKVLSGFYHTGLYSEIFFIILLTIDRYLAIVHAVALRARTVTFGVITSIVTWGLAVLVALPEFIFYGTEE-LF-ETLCSAIYPQDTVYSWRHFHTLKMTILCLALPLLVMAICYTGIIKTL---------LKCPSKKKYKAIRLIFVIMAVFFIFWTPYNVAILISTYQSILKHVDLVVLVTEVIAYSHCCVNPVIYAFVGERFRKYLRHFFHRHVLMHLGRYIPFLT-------SSVS
>PE24_HUMAN
--------------------MSTPGVNSSASLSPDRLNSPVTIPAVMFIFGVVGNLVAIVVLCKSRKKETTFYTLVCGLAVTDLLGTLLVSPVTIATYMKGQWPGGQPLCEYSTFILLFFSLSGLSIICAMSVERYLAINHAYYSHYVDKRLAGLTLFAVYASNVLFCALPNMGLGSSRLQYPDTWCFID---WTTNVTAHAAYSYMYAGFSSFLILATVLCNVLVCGALLRMAAAASSFRRIAGAEIQMVILLIATSLVVLICSIPLVVRVFVNQLYQPSEVSKNPDLQAIRIASVNPILDPWIYILLRKTVLSKAIEKIKCLFCRIGGSRRERSQRTSSAMSGHSR
>GPR4_PIG
--------------------MGNGTWEGCHVDSRVDHLFPPSLYIFVIGVGLPTNCLRLWAAYRQVRQRNELGVYLMNLSIADLLYICTLPLWVDYFLHHDNWIHGPGSCKLFGFIFYTNIYISIAFLCCISVDRYLAVAHPLFARLRRVKTAVAVSSVVWATELGANSVPLFHDEL--NH---TFCFEKPMEGWVAWMVAWMNLYRVFVGFLFPWALMLLSYRGILRAVRG------SVSTERQEKAKIKRLALSLIAIVLVCFAPYHVLLLSRSAVYLGERVFSAYHSSLAFTSLNCVADPILYCLVNEGARSDVAKALHNLLRFLTSDKPQEMDTPLTSKRNSMA
>PE23_BOVIN
PFCTRFNHSDPGIWAAERAVEAPNNLTLPPEPSEDCGSVSVAFSMTMMITGFVGNALAITLVSKSYRRKKSFLLCIGWLALTDMVGQLLTSPVVIVLYLSHQLDPSGRLCTFFGLTMTVFGLSSLFIASAMAVERALATRAPHYSSHMKTSVTRAVLLGVWLAVLAFALLPVLGVG----QYTGTWCFISNGTNSRQNWGNVFFASAFAILGLSALVVTFACNLATIKALVSRAKASQSSAQWGRITTETAIQLMGIMCVLSVCWSPLLIMMLKMIFNHTSQDECNFFLIAVRLASLNQILDPWVYLLLRKILLQKFCQLLKGHSYGLDTEGGTENNLYISNLSRFFI
>OPSD_RANCA
MNGTEGPNFYVPMSNKTGIVRSPFEYPQYYLAEPWKYSVLAAYMFLLILLGLPINFMTLYVTIQHKKLRTPLNYILLNLAFANHFMVLCGFTITMYTSLHGYFVFGQTGCYFEGFFATLGGEIALWSLVVLAIERYIVVCKPMSNFRFGENHAMMGVAFTWIMALACAVPPLFGWSRYIPEGMQCSCGVDYYTLKPEVNNESFVIYMFVVHFLIPLIIISFCYGRLVCTVKEAAAQQQESATTQKAEKEVTRMVVIMVIFFLICWVPYAYVAFYIFTHQGSEFGPIFMTVPAFFAKSSAIYNPVIYIMLNKQFRNCMITTLCCGKNPFGDEDASSAATSVSTSQVSPA
>OPS6_DROME
GR----NLSLAESVPAEIMHMVDPYWYQWPPLEPMWFGIIGFVIAILGTMSLAGNFIVMYIFTSSKGLRTPSNMFVVNLAFSDFMMMFTMFPPVVLNGFYGTWIMGPFLCELYGMFGSLFGCVSIWSMTLIAYDRYCVIVKGMARKPLTATAAVLRLMVVWTICGAWALMPLFGWNRYVPEGNMTACGTD--YFAKDWWNRSYIIVYSLWVYLTPLLTIIFSYWHIMKAVAAHNVANSEADKSKAIEIKLAKVALTTISLWFFAWTPYTIINYAGIFESMH-LSPLSTICGSVFAKANAVCNPIVYGLSHPKYKQVLREKMPCLACGKDDLTSDSRSESQA-------
>OPRM_HUMAN
LSHLDGNLSDPCGPNRTDLGGRDSLCPPTGSPSMITAITIMALYSIVCVVGLFGNFLVMYVIVRYTKMKTATNIYIFNLALADALATSTLPFQSVNYLMGT-WPFGTILCKIVISIDYYNMFTSIFTLCTMSVDRYIAVCHPVALDFRTPRNAKIINVCNWILSSAIGLPVMFMATTKYGS---IDCTLT-FSHPTWYWENLLKICVFIFAFIMPVLIITVCYGLMILRLKSVRML-SGSKEKDRNLRRITRMVLVVVAVFIVCWTPIHIYVIIKALVTIPTFQTVSWHFCIALGYTNSCLNPVLYAFLDENFKRCFREFCIPTSSNIEQQNSTRIPSTANTVDRTNH
>AA2A_MOUSE
----------------------------------MGSSVYIMVELAIAVLAILGNVLVCWAVWINSNLQNVTNFFVVSLAAADIAVGVLAIPFAITISTG--FCAACHGCLFIACFVLVLTQSSIFSLLAIAIDRYIAIRIPLYNGLVTGMKAKGIIAICWVLSFAIGLTPMLGWNNCST-K-RVTCLFE-----DVVPMNYMVYYNFFAFVLLPLLLMLAIYLRIFLAARRQESQGERTRSTLQKEVHAAKSLAIIVGLFALCWLPLHIINCFTFFCSTCHAPPWLMYLAIILSHSNSVVNPFIYAYRIREFRQTFRKIIRTHVLRRQEPFRAGGAHSTEGEQVSLR
>OPSG_CAVPO
DAYEDSTQASLFTYTNSNNTRGPFEGPNYHIAPRWVYHLTSAWMTIVVIASIFTNGLVLVATMRFKKLRHPLNWILVNLAVADLAETVIASTISVVNQVYGYFVLGHPLCVVEGYTVSLCGITGLWSLAIISWERWLVVCKPFGNVRFDAKLAIVGIVFSWVWSAVWTAPPIFGWSRYWPYGLKTSCGPDVFSGTSYPGVQSYMMVLMVTCCITPLSIIVLCYLHVWLAIRAVAKQQKESESTQKAEKEVTRMVVVMVLAYCLCWGPYAFFACFATANPGYSFHPLVAALPAYFAKSATIYNPIIYVFMNRQFRNCILQLFGKKVEDSSELSSTSRSVSPAA------
>CCR4_CERTO
IYTSDNYTEEMG-SGDYDSIKEPC---FREKNAHFNRIFLPTIYSIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVADLLFVITLPFWAVDAVAN--WYFGNFLCKAVHVIYTVNLYSSVLILAFISLDRYLAIVHATSQKPRKLLAEKVVYVGVWIPALLLTIPGFIFASVSE-DD-RFICDRF---YPNDLWVVVFQFQHIMVGLILPGIVILSCYCIIISKL---------SHSKGHQKRKALKTTVILILAFFACWLPYYIGISIDSFILLENTVHKWISITEALAFFHCCLNPILYAFLGAKFKTSAQHALTSVSRGSSLKILSKG----------GH
>OPSB_ANOCA
MNGTEGINFYVPLSNKTGLVRSPFEYPQYYLAEPWKYKVVCCYIFFLIFTGLPINILTLLVTFKHKKLRQPLNYILVNLAVADLFMACFGFTVTFYTAWNGYFIFGPIGCAIEGFFATLGGQVALWSLVVLAIERYIVVCKPMGNFRFSATHALMGISFTWFMSFSCAAPPLLGWSRYIPEGMQCSCGPDYYTLNPDYHNESYVLYMFGVHFVIPVVVIFFSYGRLICKVREAAAQQQESASTQKAEREVTRMVILMVLGFLLAWTPYAMVAFWIFTNKGVDFSATLMSVPAFFSKSSSLYNPIIYVLMNKQFRNCMITTICCGKNPFGDEDVSSSVSSVSSSQVSPA
>O1A2_HUMAN
-------------MKKENQSFNLDFILLGVTSQQEQNNVFFVIFLCIYPITLTGNLLIILAICADIRLHNPMYFLLANLSLVDIIFSSVTIPKVLANHLLGSKFISFGGCLMQMYFMIALAKADSYTLAAMAYDRAVAISCPLYTTIMSPRSCILLIAGSWVIGNTSALPHTLLTASL---SANFYCDIMKLSCSDVVFFNVKMMYLGVGVFSLPLLCIIVSYVQVFSTV--------FQVPSTKSLFKAFCTCGSHLTVVFLYYGTTMGMYFRPLTSY----SPKDAVITVMYVAVTPALNPFIYSLRNWDMKAALQKLFSKRISS---------------------
>A2AC_MOUSE
AEGPNGSDAGEWGSGGGANASGTDWVPPPGQYSAGAVAGLAAVVGFLIVFTVVGNVLVVIAVLTSRALRAPQNLFLVSLASADILVATLVMPFSLANELMAYWYFGQVWCGVYLALDVLFCTSSIVHLCAISLDRYWSVTQAVYNLKRTPRRVKATIVAVWLISAVISFPPLVSFYR-------PQCGLN--------DETWYILSSCIGSFFAPCLIMGLVYARIYRVAKLRRRAVCRRKVAQAREKRFTFVLAVVMGVFVLCWFPFFFSYSLYGICREAQLPEPLFKFFFWIGYCNSSLNPVIYTVFNQDFRRSFKHILFRRRRRGFRQ-----------------
>5H2C_RAT
LLVWQFDISISPVAAIVTDTFNSSDGGRLFQFPDGVQNWPALSIVVIIIMTIGGNILVIMAVSMEKKLHNATNYFLMSLAIADMLVGLLVMPLSLLAILYDYWPLPRYLCPVWISLDVLFSTASIMHLCAISLDRYVAIRNPIHSRFNSRTKAIMKIAIVWAISIGVSVPIPVIGLRD-VF--NTTCVL---------NDPNFVLIGSFVAFFIPLTIMVITYFLTIYVLRRQKKKPRGTMQAINNEKKASKVLGIVFFVFLIMWCPFFITNILSVLCGKAKLMEKLLNVFVWIGYVCSGINPLVYTLFNKIYRRAFSKYLRCDYKPDKKPPVRQIALSGRELNVNIY
>OPRM_PIG
FSHLEGNLSDPCIRNRTELGGSDSLCPPTGSPSMVTAITIMALYSIVCVVGLFGNFLVMYVIVRYTKMKTATNIYIFNLALADALATSTLPFQSVNYLMGT-WPFGTILCKIVISIDYYNMFTSIFTLCTMSVDRYIAVCHPVALDFRTPRNAKIINVCNWILSSAIGLPVMFMATTKYGS---IDCALT-FSHPTWYWENLLKICVFIFAFIMPVLIITVCYGLMILRLKSVRML-SGSKEKDRNLRRITRMVLVVVAVFIVCWTPIHIYVIIKALITIPTFQTVSWHFCIALGYTNSCLNPVLYAFLDENFKRCFREFCIPTSSTIEQQNSARIPSTANTVDRTNH
>B1AR_FELCA
LPDGAATAARLLVPASPSASPLTPTSEGPAPLSQQWTAGIGLLMALIVLLIVAGNVLVIVAIAKTPRLQTLTNLFIMSLASADLVMGLLVVPFGATIVMRGRWEYGSFFCELWTSVDVLCVTASIETLCVIALDRYLAITSPFYQSLLTRARARALVCTVWAISALVSFLPILMHWWRAR-R--KCCDFV--------TNRAYAIASSVVSFYVPLCIMAFVYLRVFREAQKQNGRRRPSRLVALREQKALKTLGIIMGVFTLCWLPFFLANVVKAFHRDL-VPDRLFVFFNWLGYANSAFNPIIYCRSP-DFRKAFQRLLCFARRAARGGHAAAGCLPGTRPPPSPG
>BONZ_HUMAN
------MAEHDYHEDYGFSSFNDSSQEEHQDFLQFSKVFLPCMYLVVFVCGLVGNSLVLVISIFYHKLQSLTDVFLVNLPLADLVFVCTLPFWAYAGIHE--WVFGQVMCKSLLGIYTINFYTSMLILTCITVDRFIVVVKATNQQAKRMTWGKVTSLLIWVISLLVSLPQIIYGNVFNLD--KLICGYH-----DEAISTVVLATQMTLGFFLPLLTMIVCYSVIIKTL---------LHAGGFQKHRSLKIIFLVMAVFLLTQMPFNLMKFIRSTHWEYTSFHYTIMVTEAIAYLRACLNPVLYAFVSLKFRKNFWKLVKDIGCLPYLGVSHQWKTFSASHNVEAT
>CH23_HUMAN
EDEDYNTSISYGDEYPDYLDSIVVLEDLSPLEARVTRIFLVVVYSIVCFLGILGNGLVIIIATFKMK-KTVNMVWFLNLAVADFLFNVFLPIHITYAAMDYHWVFGTAMCKISNFLLIHNMFTSVFLLTIISSDRCISVLLPVSQNHRSVRLAYMACMVIWVLAFFLSSPSLVFRDTAN-SS--WPTHSQ-MDPVGYSRHMVVTVTRFLCGFLVPVLIITACYLTIVCKL---------QRNRLAKTKKPFKIIVTIIITFFLCWCPYHTLNLLELHHTAMSVFSLGLPLATALAIANSCMNPILYVFMGQDFKK-FKVALFSRLVNALSEDTGHSFTKMSSMNERTS
>TSHR_HUMAN
LQAFDSHYDYTICGDSEDMVCTPKSDEFNPCEDIMGYKFLRIVVWFVSLLALLGNVFVLLILLTSHYKLNVPRFLMCNLAFADFCMGMYLLLIASVDLYTHSWQTG-PGCNTAGFFTVFASELSVYTLTVITLERWYAITFAMLDRKIRLRHACAIMVGGWVCCFLLALLPLVGISSY----KVSICLPM-----TETPLALAYIVFVLTLNIVAFVIVCCCHVKIYITVRNP------QYNPGDKDTKIAKRMAVLIFTDFICMAPISFYALSAILNKPLITVSNSKILLVLFYPLNSCANPFLYAIFTKAFQRDVFILLSKFGICKRQAQAYRGSTDIQVQKVTHD
>GHSR_PIG
SEEPGPNLTLPDLGWDAPPENDSLVEELLPLFPTPLLAGVTATCVALFVVGIAGNLLTMLVVSRFREMRTTTNLYLSSMAFSDLLIFLCMPLDLFRLWQYRPWNLGNLLCKLFQFVSESCTYATVLTITALSVERYFAICFPLAKVVVTKGRVKLVILVIWAVAFCSAGPIFVLVG---T-D-TNECRAT---FAVRSGLLTVMVWVSSVFFFLPVFCLTVLYSLIGRKLWRRGE-AVGSSLRDQNHKQTVKMLAVVVFAFILCWLPFHVGRYLFSKSLEPQISQYCNLVSFVLFYLSAAINPILYNIMSKKYRVAVFKLLGFEPFSQRKLSTLKDESSINT------
>B1AR_RAT
LPDGAATAARLLVLASPPASLLPPASEGSAPLSQQWTAGMGLLLALIVLLIVVGNVLVIVAIAKTPRLQTLTNLFIMSLASADLVMGLLVVPFGATIVVWGRWEYGSFFCELWTSVDVLCVTASIETLCVIALDRYLAITLPFYQSLLTRARARALVCTVWAISALVSFLPILMHWWRAR-R--KCCDFV--------TNRAYAIASSVVSFYVPLCIMAFVYLRVFREAQKQNGRRRPSRLVALREQKALKTLGIIMGVFTLCWLPFFLANVVKAFHRDL-VPDRLFVFFNWLGYANSAFNPIIYCRSP-DFRKAFQRLLCCARRAACRRRAAHGCLARAGPPPSPG
>ACM2_CHICK
-----------MNNSTYINSSSENVIALESPYKTIEVVFIVLVAGSLSLVTIIGNILVMVSIKVNRHLQTVNNYFLFSLACADLIIGIFSMNLYTLYTVIGYWPLGPVVCDLWLALDYVVSNASVMNLLIISFDRYFCVTKPLYPVKRTTKMAGMMIAAAWVLSFILWAPAILFWQFIV-VP-DKDCYIQ------FFSNPAVTFGTAIAAFYLPVIIMTVLYWQISRASKSRVKMPAKKKPPPSREKKVTRTILAILLAFIITWTPYNVMVLINSFCASC-IPGTVWTIGYWLCYINSTINPACYALCNATFKKTFKHLLMCHYKNIGATR----------------
>GPRK_HUMAN
ATAVTTVRTNASGLEVPLFHLFARLDEELHGTFPGLCVALMAVHGAIFLAGLVLNGLALYVFCCRTRAKTPSVIYTINLVVTDLLVGLSLPTRFAVYYGA---RGCLRCAFPHVLGYFLNMHCSILFLTCICVDRYLAIVRPEPAACRQPACARAVCAFVWLAAGAVTLSVLGVTG----------S-----------RPCCRVFALTVLEFLLPLLVISVFTGRIMCALSRP----GLLHQGRQRRVRAMQLLLTVLIIFLVCFTPFHARQVAVALWPDMHTSLVVYHVAVTLSSLNSCMDPIVYCFVTSGFQATVRGLFGQHGEREPSSGDVVSSGRHHILSAGPH
>ACM1_HUMAN
------------MNTSAPPAVSPNITVLAPGKGPWQVAFIGITTGLLSLATVTGNLLVLISFKVNTELKTVNNYFLLSLACADLIIGTFSMNLYTTYLLMGHWALGTLACDLWLALDYVASNASVMNLLLISFDRYFSVTRPLYRAKRTPRRAALMIGLAWLVSFVLWAPAILFWQYLV-VL-AGQCYIQ------FLSQPIITFGTAMAAFYLPVTVMCTLYWRIYRETENRRGKAKRKTFSLVKEKKAARTLSAILLAFILTWTPYNIMVLVSTFCKDC-VPETLWELGYWLCYVNSTINPMCYALCNKAFRDTFRLLLLCRWDKRRWRKIPKRPSRQC-------
>ML1A_HUMAN
-------MQGNGSALPNASQ---PVLRGDGARPSWLASALACVLIFTIVVDILGNLLVILSVYRNKKLRNAGNIFVVSLAVADLVVAIYPYPLVLMSIFNNGWNLGYLHCQVSGFLMGLSVIGSIFNITGIAINRYCYICHSLYDKLYSSKNSLCYVLLIWLLTLAAVLPNLRAGT-LQYDP-IYSCTFA------QSVSSAYTIAVVVFHFLVPMIIVIFCYLRIWILVLQVQR-PDRKPKLKPQDFRNFVTMFVVFVLFAICWAPLNFIGLAVASDPASRIPEWLFVASYYMAYFNSCLNAIIYGLLNQNFRKEYRRIIVSLCTARVFFVDSSNWKPSPLMTNNNV
>TSHR_RAT
LQAFDSHYDYTVCGDNEDMVCTPKSDEFNPCEDIMGYKFLRIVVWFVSPMALLGNVFVLFVLLTSHYKLTVPRFLMCNLAFADFCMGVYLLLIASVDLYTHTWQTG-PGCNTAGFFTVFASELSVYTLTVITLERWYAITFAMLDRKIRLRHAYTIMAGGWVSCFLLALLPMVGISSY----KVSICLPM-----TDTPLALAYIALVLLLNVVAFVIVCSCYVKIYITVRNP------QYNPRDKDTKIAKRMAVLIFTDFMCMAPISFYALSALMNKPLITVTNSGVLLVLFYPLNSCANPFLYAIFTKAFQRDVFILLSKFGLCKHQAQAYQANTGIQIQKIPQD
>OPS2_DROME
AQSS-GNGSVLDNVLPDMAHLVNPYWSRFAPMDPMMSKILGLFTLAIMIISCCGNGVVVYIFGGTKSLRTPANLLVLNLAFSDFCMMASQSPVMIINFYYETWVLGPLWCDIYAGCGSLFGCVSIWSMCMIAFDRYNVIVKGINGTPMTIKTSIMKILFIWMMAVFWTVMPLIGWSAYVPEGNLTACSID--YMTRMWNPRSYLITYSLFVYYTPLFLICYSYWFIIAAVAAHMNVRSSEDCDKSAEGKLAKVALTTISLWFMAWTPYLVICYFGLFKIDG-LTPLTTIWGATFAKTSAVYNPIVYGISHPKYRIVLKEKCPMCVFGNTDEPKPDATSEADSKA----
>OPRD_HUMAN
PPLFANASDAYPSACPSAGANASGPPGARSASSLALAIAITALYSAVCAVGLLGNVLVMFGIVRYTKMKTATNIYIFNLALADALATSTLPFQSAKYLMET-WPFGELLCKAVLSIDYYNMFTSIFTLTMMSVDRYIAVCHPVALDFRTPAKAKLINICIWVLASGVGVPIMVMAVTRPGA---VVCMLQ-FPSPSWYWDTVTKICVFLFAFVVPILIITVCYGLMLLRLRSVRLL-SGSKEKDRSLRRITRMVLVVVGAFVVCWAPIHIFVIVWTLVDIDPLVVAALHLCIALGYANSSLNPVLYAFLDENFKRCFRQLCRKPCGRPDPSSFSRARVTACTPSDGPG
>OPSD_.ENLA
MNGTEGPNFYVPMSNKTGVVRSPFDYPQYYLAEPWQYSALAAYMFLLILLGLPINFMTLFVTIQHKKLRTPLNYILLNLVFANHFMVLCGFTVTMYTSMHGYFIFGPTGCYIEGFFATLGGEVALWSLVVLAVERYIVVCKPMANFRFGENHAIMGVAFTWIMALSCAAPPLFGWSRYIPEGMQCSCGVDYYTLKPEVNNESFVIYMFIVHFTIPLIVIFFCYGRLLCTVKEAAAQQQESLTTQKAEKEVTRMVVIMVVFFLICWVPYAYVAFYIFTHQGSNFGPVFMTVPAFFAKSSAIYNPVIYIVLNKQFRNCLITTLCCGKNPFGDEDGSSAASSVSSSQVSPA
>PAFR_HUMAN
----------------------MEPHDSSHMDSEFRYTLFPIVYSIIFVLGVIANGYVLWVFARLYPKFNEIKIFMVNLTMADMLFLITLPLWIVYYQNQGNWILPKFLCNVAGCLFFINTYCSVAFLGVITYNRFQAVTRPITAQANTRKRGISLSLVIWVAIVGAASYFLILDSTN-G--NVTRCFEH---YEKGSVPVLIIHIFIVFSFFLVFLIILFCNLVIIRTLLMQP---VQQQRNAEVKRRALWMVCTVLAVFIICFVPHHVVQLPWTLAELGQAINDAHQVTLCLLSTNCVLDPVIYCFLTKKFRKHLTEKFYSMRSSRKCSRATTDPFNQIPGNSLKN
>OPSG_SCICA
DSHEDSTQSSIFTYTNSNATRGPFEGPNYHIAPRWVYHITSTWMIIVVIASVFTNGLVLVATMKFKKLRHPLNWILVNLAIADLAETVIASTISVVNQLYGYFVLGHPLCVVEGYTVSVCGITGLWSLAIISWERWLVVCKPFGNMRFDAKLAIVGIAFSWIWSAVWTAPPIFGWSRYWPYGLKTSCGPDVFSGTSYPGVQSYMMVLMVTCCIIPLSIIILCYLQVWLAIRAVAKQQKESESTQKAEKEVTRMVVVMVFAYCLCWGPYTFFACFATANPGYAFHPLVAALPAYFAKSATIYNPIIYVFMNRQFRNCILQLFGKKVDDTSELSSASKSVSPAA------
>HH2R_RAT
-------------------MEPNGTVHSCCLDSMALKVTISVVLTTLILITIAGNVVVCLAVSLNRRLRSLTNCFIVSLAATDLLLGLLVLPFSAIYQLSFTWSFGHVFCNIYTSLDVMLCTASILNLFMISLDRYCAVTDPLYPVLVTPVRVAISLVFIWVISITLSFLSIHLGWN--RN-TF-KCKVQ--------VNEVYGLVDGLVTFYLPLLIMCVTYYRIFKIAREQR--ISSWKAATIREHKATVTLAAVMGAFIICWFPYFTAFVYRGLRGDD-INEAVEGIVLWLGYANSALNPILYAALNRDFRTAYQQLFHCKFASHNSHKTSLRRSQSREGRW---
>V1AR_RAT
SSPWWPLTTEGSNGSQEAARLGEGDSPLGDVRNEELAKLEIAVLAVIFVVAVLGNSSVLLALHRTPRKTSRMHLFIRHLSLADLAVAFFQVLPQLCWDITSSFRGPDWLCRVVKHLQVFAMFASAYMLVVMTADRYIAVCHPLKTLQQPARRSRLMIATSWVLSFILSTPQYFIFSVIETK--TQDCWAT---FIQPWGTRAYVTWMTSGVFVAPVVVLGTCYGFICYHIWRNLLVVSSVKSISRAKIRTVKMTFVIVSAYILCWAPFFIVQMWSVWDENFDSENPSITITALLASLNSCCNPWIYMFFSGHLLQDCVQSFPCCHSMAQKFAKDDSTSYSNNRSPTNS
>GPR1_HUMAN
EDLEETLFEEFENYSYDLDYYSLESDLEEKVQLGVVHWVSLVLYCLAFVLGIPGNAIVIWFTGLKWK-KTVTTLWFLNLAIADFIFLLFLPLYISYVAMNFHWPFGIWLCKANSFTAQLNMFASVFFLTVISLDHYIHLIHPVSHRHRTLKNSLIVIIFIWLLASLIGGPALYFR--D--NN-HTLCYNNHDPDLTLIRHHVLTWVKFIIGYLFPLLTMSICYLCLIFKV---------KKRTVLISSRHFWTILVVVVAFVVCWTPYHLFSIWELTIHHNHVMQAGIPLSTGLAFLNSCLNPILYVLISKKFQARFRSSVAEILKYTLWEVSCSGNSETKNLCLLET
>DADR_MACMU
--------------MRTLNTSAMDGTGLVVERDFSVRILTACFLSLLILSTLLGNTLVCAAVIRFRHRSKVTNFFVISLAVSDLLVAVLVMPWKAVAEIAGFWPFG-SFCNIWVAFDIMCSTASILNLCVISVDRYWAISSPFYERKMTPKAAFILISVAWTLSVLISFIPVQLSWHKAGN-TIDNCDSS--------LSRTYAISSSVISFYIPVAIMIVTYTRIYRIAQKQVECESSFKMSFKRETKVLKTLSVIMGVFVCCWLPFFILNCILPFCGSGCIDSITFDVFVWFGWANSSLNPIIYAFNA-DFRKAFSTLLGCYRLCPATNNAIETAAMFSSH-----
>AG2R_CANFA
------MILNSSTEDGIKRIQDDC---PKAGRHNYIFVMIPTLYSIIFVVGIFGNSLVVIVIYFYMKLKTVASVFLLNLALADLCFLLTLPLWAVYTAMEYRWPFGNYLCKIASASVSFNLYASVFLLTCLSIDRYVAIVHPMSPVRRTMLMAKVTCIIIWLLAGLASLPTIIHRNVFF-TN--TVCAFH-YESQNSTLPIGLGLTKNILGFLFPFLIILTSYTLIWKTLKRA----YEIQKNKPRNDDIFKIIMAIVLFFFFSWVPHQIFTFLDVLIQLGDIVDTAMPITICIAYFNNCLNPLFYGFLGKKFKKYFLQLLKYIPPKAKSHSSLSTRPSD-------H
>OPSD_ATHBO
MNGTEGPYFYIPMLNTTGVVRSPYEYPQYYLVNPAAYAVLGAYMFFLILVGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTIYTSMHGYFVLGRLGCNVEGFSATLGGEIALWSLVVLAIERWVVVCKPISNFRFGENHAIMGVAFTWFMAAACAVPPLFGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFTCHFCIPLMVVFFCYGRLVCAVKEAAAAQQESETTQRAEREVTRMVIIMVVSFLVSWVPYASVAWYIFTHQGSEFGPLFMTIPAFFAKSSSIYNPMIYICMNKQFRHCMITTLCCGKNPFEEEEGASSSSVSSSSVSPAA
>RDC1_RAT
YVEPGNYSDSNWPCNSSDCIVVDTVQCPAMPNKNVLLYTLSFIYIFIFVIGMIANSVVVWVNIQAKTTGYDTHCYILNLAIADLWVVITIPVWVVSLVQHNQWPMGELTCKITHLIFSINLFGSIFFLACMSVDRYLSITYFTTSSYKKKMVLRVVCVLVWLLAFFVSLPDTYYLKTVTNNE--TYCRSFYPEHSIKEWLIGMELVSVILGFAVPFTIIAIFYFLLARAM---------SASGDQEKHSSRKIIFSYVVVFLVCWLPYHFVVLLDIFSILHNVLFTALHVTQCLSLVHCCVNPVLYSFINRNYRYELMKAFIFKYSAKTGLTKLIDEYSALEQNAKA-
>IL8A_RAT
EGDFEEEFGNITRMLPTGEYFSPC----KR-VPMTNRQAVVVFYALVFLLSLLGNSLVMLVILYRRRTRSVTDVYVLNLAIADLLFSLTLPFLAVSKWKG--WIFGTPLCKMVSLLKEVNFFSGILLLACISVDRYLAIVHATRTLTRKRYLVKFVCMGTWGLSLVLSLPFAIFRQAYK-RS-GTVCYEV-LGEATADLRITLRGLSHIFGFLLPLFIMLVCYGLTLRTL---------FKAHMRQKRRAMWVIFAVVLVFLLCCLPYNLVLLSDTLLGAHNNIDQALYITEILGFSHSCLNPVIYAFVGQSFRHEFLKILAN--LVHKEVLTHHS------------
>CB2R_RAT
----MEGCRELELTNGSNGGLEFNPMKEYMILSDAQQIAVAVLCTLMGLLSALENVAVLYLILSSQRRRKPSYLFIGSLAGADFLASVIFACNFVIFHVFHG-VDSRNIFLLKIGSVTMTFTASVGSLLLTAVDRYLCLCYPPYKALVTRGRALVALGVMWVLSALISYLPLMGWTC-----CPSPCSEL------FPLIPNDYLLGWLLFIAILFSGIIYTYGYVLWKAHQHTEHQVPGIARMRLDVRLAKTLGLVMAVLLICWFPALALMGHSLVTTLSDKVKEAFAFCSMLCLVNSMVNPIIYALRSGEIRSAAQHCLTGWKKYLQGLGSEGKVTETEAEVKTTT
>LSHR_HUMAN
ESELSGWDYEYGFCLPKTPRCAPEPDAFNPCEDIMGYDFLRVLIWLINILAIMGNMTVLFVLLTSRYKLTVPRFLMCNLSFADFCMGLYLLLIASVDSQTKGWQTG-SGCSTAGFFTVFASELSVYTLTVITLERWHTITYAILDQKLRLRHAILIMLGGWLFSSLIAMLPLVGVSNY----KVSICFPM-----VETTLSQVYILTILILNVVAFFIICACYIKIYFAVRNP------ELMATNKDTKIAKKMAILIFTDFTCMAPISFFAISAAFKVPLITVTNSKVLLVLFYPINSCANPFLYAIFTKTFQRDFFLLLSKFGCCKRRAELYRRSNCKNGFTGSNK
>OPSD_LIZSA
MNGTEGPYFYIPMVNTTGIVRSPYEYPQYYLVNPAAYAALGAYMFLLILVGFPINFLTLYVTIEHKKLRTPLNYILLNLAVANLFMVFGGFTTTMYTSMHGYFVLGRLGCNLEGFFATLGGEIALWSLVVLAIERWMVVCKPISNFRFGEDHAIMGLAFTWVMAAACAVPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVIYMFVCHFLIPLVVVFFCYGRLLCAVKEAAAAQQESETTQRAEREVSRMVVIMVVAFLICWCPYAGVAWYIFTHQGSEFGPLFMTFPAFFAKSSSIYNPMIYICMNKQFRHCMITTLCCGKNPFEEEEGASTSVSSSSVSPAA-
>CB1R_FELCA
EFYNKSLSSYKENEENIQCGENFMDMECFMILNPSQQLAIAVLSLTLGTFTVLENLLVLCVILHSRSRCRPSYHFIGSLAVADLLGSVIFVYSFVDFHVFHR-KDSPNVFLFKLGGVTASFTASVGSLFLTAIDRYISIHRPLYKKIVTRPKAVVAFCLMWTIAIVIAVLPLLGWNCK-K--LQSVCSDI------FPLIDETYLMFWIGVTSVLLLFIVYAYMYILWKAHIHSEDQVTRPDQARMDIRLAKTLVLILVVLIICWGPLLAIMVYDVFGKMNKLIKTVFAFCSMLCLLNSTVNPIIYALRSKDLRHAFRSMFPSCEGTAQPLDNSMGHANNTANVHRAA
>EDG3_HUMAN
TALPPRLQPVRGNETLREHYQYVGKLAGRLKEASEGSTLTTVLFLVICSFIVLENLMVLIAIWKNNKFHNRMYFFIGNLALCDLLAG-IAYKVNILMSGKKTFSLSPTVWFLREGSMFVALGASTCSLLAIAIERHLTMIKMRPYDANKRHRVFLLIGMCWLIAFTLGALPILGWNCL-H--NLPDCSTI------LPLYSKKYIAFCISIFTAILVTIVILYARIYFLVKSS---KVANHNNSERSMALLRTVVIVVSVFIACWSPLFILFLIDVACRVQCPILFKAQWFIVLAVLNSAMNPVIYTLASKEMRRAFFRLVCNCLVRGRGARASPIRSKSSSSNNSSH
>UR2R_HUMAN
AATGSSVPEPPGGPNATLNSSWASPTEPSSLEDLVATGTIGTLLSAMGVVGVVGNAYTLVVTCRSLRAVASMYVYVVNLALADLLYLLSIPFIVATYVTKE-WHFGDVGCRVLFGLDFLTMHASIFTLTVMSSERYAAVLRPLDTVQRPKGYRKLLALGTWLLALLLTLPVMLAMR---GP--KSLCLPA----WGPRAHRAYLTLLFATSIAGPGLLIGLLYARLARAYRRSQR--ASFKRARRPGARALRLVLGIVLLFWACFLPFWLWQLLAQYHQAPRTARIVNYLTTCLTYGNSCANPFLYTLLTRNYRDHLRGRVRGPGSGGGRGPVPSLRCSGRSLSSCSP
>AA2B_CHICK
--------------------------------MNTMKTTYIVLELIIAVLSIAGNVLVCWAVAINSTLKNATNYFLVSLAVADIAVGLLAIPFAITISIG--FQVDFHSCLFFACFVLVLTQSSIFSLLAVAIDRYLAIKIPLYNSLVTGKRARGLIAVLWLLSFVIGLTPLMGWNKAMG-A-FISCLFE-----NVVTMSYMVYFNFFGCVLLPLIIMLGIYIKIFMVACKQ---MGNSRTTLQKEVHAAKSLAIIVGLFAFCWLPLHILNCITHFHEEFSKPEWVMYVAIILSHANSVINPIIYAYRIRDFRYTFHKIISKILCKTDDFPKCTTVTNVNAPAASVT
>5H2A_MOUSE
NTSEASNWTIDAENRTNLSCEGYLPPTCLSILHLQEKNWSALLTTVVIILTIAGNILVIMAVSLEKKLQNATNYFLMSLAIADMLLGFLVMPVSMLTILYGYWPLPSKLCAVWIYLDVLFSTASIMHLCAISLDRYVAIQNPIHSRFNSRTKAFLKIIAVWTISVGISMPIPVFGLQD-VF---GSCLL---------ADDNFVLIGSFVAFFIPLTIMVITYFLTIKSLQKEEPGGRRTMQSISNEQKACKVLGIVFFLFVVMWCPFFITNIMAVICKESNVIGALLNVFVWIGYLSSAVNPLVYTLFNKTYRSAFSRYIQCQYKENRKPLQLILAYKSSQLQVGQK
>US27_HCMVA
-------MTTSTNNQTLTQVSNMTNHTLNSTEIYQLFEYTRLGVWLMCIVGTFLNVLVITTILYYRRKKSPSDTYICNLAVADLLIVVGLPFFLEYAKHHP-KLSREVVCSGLNACFYICLFAGVCFLINLSMDRYCVIVWGVLNRVRNNKRATCWVVIFWILAVLMGMPHYLMYSHT-----NNECVGE-ANETSGWFPVFLNTKVNICGYLAPIALMAYTYNRMVRFI---------INYVGKWHMQTLHVLLVVVVSFASFWFPFNLALFLESIRLLANVIIFCLYVGQFLAYVRACLNPGIYILVGTQMRKDMWTTLRVFACCCVKQEIPYQKDIQRRAKHTKR
>NK1R_HUMAN
-----MDNVLPVDSDLSPNISTNTSEPNQFVQPAWQIVLWAAAYTVIVVTSVVGNVVVMWIILAHKRMRTVTNYFLVNLAFAEASMAAFNTVVNFTYAVHNEWYYGLFYCKFHNFFPIAAVFASIYSMTAVAFDRYMAIIHPL-QPRLSATATKVVICVIWVLALLLAFPQGYYST---SR---VVCMIEWPEHPNKIYEKVYHICVTVLIYFLPLLVIGYAYTVVGITLE---IPSDRYHEQVSAKRKVVKMMIVVVCTFAICWLPFHIFFLLPYINPDLKFIQQVYLAIMWLAMSSTMYNPIIYCCLNDRFRLGFKHAFRCCPFISAGDY----STRYLQTQGSVY
>5H1A_HUMAN
-MDVLSPGQGNNTTSPPAPFETGGNTTGISDVTVSYQVITSLLLGTLIFCAVLGNACVVAAIALERSLQNVANYLIGSLAVTDLMVSVLVLPMAALYQVLNKWTLGQVTCDLFIALDVLCCTSSILHLCAIALDRYWAITDPIYVNKRTPRRAAALISLTWLIGFLISIPPMLGW-----DP--DACTIS--------KDHGYTIYSTFGAFYIPLLLMLVLYGRIFRAARFRKNEEAKRKMALARERKTVKTLGIIMGTFILCWLPFFIVALVLPFCESSHMPTLLGAIINWLGYSNSLLNPVIYAYFNKDFQNAFKKIIKCKFCRQ--------------------
>IL8B_MOUSE
SGDLDIFN-YSSGMPSILPDAVPC----HSENLEINSYAVVVIYVLVTLLSLVGNSLVMLVILYNRSTCSVTDVYLLNLAIADLFFALTLPVWAASKVNG--WTFGSTLCKIFSYVKEVTFYSSVLLLACISMDRYLAIVHATSTLIQKRHLVKFVCIAMWLLSVILALPILILRNPVK-LS-TLVCYED-VGNNTSRLRVVLRILPQTFGFLVPLLIMLFCYGFTLRTL---------FKAHMGQKHRAMRVIFAVVLVFLLCWLPYNLVLFTDTLMRTKDDIDKALNATEILGFLHSCLNPIIYAFIGQKFRHGLLKIMATYGLVSKEFLAKEG------------
>GPRY_MOUSE
LCSSHGMHFITNYSDQASQNFGVPNVTSCPMDEKLLSTVLTTFYSVIFLVGLVGNIIALYVFLGIHRKRNSIQIYLLNVAVADLLLIFCLPFRIMYHINQNKWTLGVILCKVVGTLFYMNMYISIILLGFISLDRYIKINRSIQRRAITTKQSIYVCCIVWTVALAGFLTMIILTLKK-----STMCFHY--RDRHNAKGEAIFNFVLVVMFWLIFLLIILSYIKIGKNLLRISKR--SKFPNSGKYATTARNSFIVLIIFTICFVPYHAFRFIYISSQLNEIIHKTNEIMLVFSSFNSCLDPVMYFLMSSNIRKIMCQLLFRRFQSEASRSESTSLHDLSVTVKMPQ
>PAR2_RAT
---LDTPPPITGKGAPVEPGFSVDEFSASVLTGKLTTVFLPVIYIIVFVIGLPSNGMALWVFFFRTKKKHPAVIYMANLALADLLSVIWFPLKISYHLHGNDWTYGDALCKVLIGFFYGNMYCSILFMTCLSVQRYWVIVNPMGHSRKRANIAVGVSLAIWLLIFLVTIPLYVMRQTIY--N--TTCHDVLPEEVLVGDMFSYFLSLAIGVFLFPALLTASAYVLMIKTL------SAMDEHSEKKRRRAIRLIITVLSMYFICFAPSNVLLVVHYFLIKSSHVYALYLVALCLSTLNSCIDPFVYYFVSKDFRDQARNALLCRSVRTVKRMQISLKSSS--------
>5H5A_MOUSE
LPVNLTSFSLSTPSSLEPNRSDTEVLRPSRPFLSAFRVLVLTLLGFLAAATFTWNLLVLATILKVRTFHRVPHNLVASMAISDVLVAVLVMPLSLVHELSGRWQLGRRLCQLWIACDVLCCTASIWNVTAIALDRYWSITRHLYTLRTRKRVSNVMILLTWALSTVISLAPLLFGWGE-S-E-SEECQVS--------REPSYTVFSTVGAFYLPLWLVLFVYWKIYRAAKFRATVTEGDTWREQKEQRAALMVGILIGVFVLCWFPFFVTELISPLCSW-DVPAIWKSIFLWLGYSNSFFNPLIYTAFNRSYSSAFKVFFSKQQ-----------------------
>IL8A_GORGO
PQMWDFDDLNFTGMPPIDEDYSPC----RLETETLNKYVVIITYALAFLLSLLGNSLVMLVILYSRGGRSVTDVYLLNLALADLLFALTLPIWAASKVNG--WIFGTFLCKVVSLLKEVNFYSGILLLACISVDRYLAIVHATRTLTQKRHLVKFVCLGCWGLSMILSLPFFLFRQAYH-NS-SPVCYEV-LGNDTAKWRMVLRILPHTFGFIVPLFVMLFCYGFTLRTL---------FKAHMGQKHRAMRVIFAVVLIFLLCWLPYNLVLLADTLMRTQNNVSLALDATEILGFLHSCLNPIIYAFIGQNFRHGFLKILAMHGLVSKEFLARHR------------
>5H2A_RAT
NTSEASNWTIDAENRTNLSCEGYLPPTCLSILHLQEKNWSALLTTVVIILTIAGNILVIMAVSLEKKLQNATNYFLMSLAIADMLLGFLVMPVSMLTILYGYWPLPSKLCAIWIYLDVLFSTASIMHLCAISLDRYVAIQNPIHSRFNSRTKAFLKIIAVWTISVGISMPIPVFGLQD-VF---GSCLL---------ADDNFVLIGSFVAFFIPLTIMVITYFLTIKSLQKEEPGGRRTMQSISNEQKACKVLGIVFFLFVVMWCPFFITNIMAVICKESNVIGALLNVFVWIGYLSSAVNPLVYTLFNKTYRSAFSRYIQCQYKENRKPLQLILAYKSSQLQVGQK
>SSR2_BOVIN
ELNETQPWLTTPFDLNGSVGAANISNQTEPYYDLASNVVLTFIYFVVCIIGLCGNTLVIYVILRYAKMKTITNIYILNLAIADELFMLGLPFLAMQVALVH-WPFGKAICRVVMTVDGINQFTSIFCLTVMSIDRYLAVVHPISAKWRRPRTAKMINVAVWGVSLLVILPIMIYAGLRSWG--RSSCTIN-WPGESGAWYTGFIIYAFILGFLVPLTIICLCYLFIIIKVKSSGIR-VGSSKRKKSEKKVTRMVSIVVAVFIFCWLPFYIFNVSSVSVAISPALKGMFDFVVVLTYANSCANPILYAFLSDNFKKSFQNVLCLVKVSGTDDGERSDLNETTETQRTLL
>GP42_HUMAN
-----------------------MDTGPDQSYFSGNHWFVFSVYLLTFLVGLPLNLLALVVFVGKLRRPVAVDVLLLNLTASDLLLLLFLPFRMVEAANGMHWPLPFILCPLSGFIFFTTIYLTALFLAAVSIERFLSVAHPLYKTRPRLGQAGLVSVACWLLASAHCSVVYVIEFSGD--T--GTCYLE-FWKDQLAILLPVRLEMAVVLFVVPLIITSYCYSRLVWILGR--------GGSHRRQRRVAGLVAATLLNFLVCFGPYNVSHVVGYICG---ESPVWRIYVTLLSTLNSCVDPFVYYFSSSGFQADFHELLRRLCGLWGQWQQESSGGEEQRADRPAE
>APJ_MACMU
---------MEEGGDFDNYYGADNQSECEYTDWKSSGALIPAIYMLVFLLGTTGNGLVLWTVFRSSRKRRSADIFIASLAVADLTFVVTLPLWATYTYRDYDWPFGTFSCKLSSYLIFVNMYASVFCLTGLSFDRYLAIVRPVNARLRLRVSGAVATAVLWVLAALLAMPVMVFRTTGDQCY-MDYSMVA-TVSSDWAWEVGLGVSSTTVGFVVPFTIMLTCYFFIAQTIAGHFR--KERIEGLRKRRRLLSIIVVLVVTFALCWMPYHLVKTLYMLGSLLLFLMNVFPYCTCISYVNSCLNPFLYAFFDPRFRQACTSMLCCGQSRCAGTSHSSSSSGHSQGPGPNM
>AG22_MOUSE
RNITSSRPFDNLNATGTNESAFNC----SHKPSDKHLEAIPVLYYMIFVIGFAVNIVVVSLFCCQKGPKKVSSIYIFNLALADLLLLATLPLWATYYSYRYDWLFGPVMCKVFGSFLTLNMFASIFFITCMSVDRYQSVIYPFLSQRRNPWQASYVVPLVWCMACLSSLPTFYFRDVRT-LG--NACIMAFPPEKYAQWSAGIALMKNILGFIIPLIFIATCYFGIRKHLLKT----NSYGKNRITRDQVLKMAAAVVLAFIICWLPFHVLTFLDALTWMGAVIDLALPFAILLGFTNSCVNPFLYCFVGNRFQQKLRSVFRVPITWLQGKRETMSREMD-------T
>FSHR_SHEEP
FDMMYSEFDYDLCSEVVDVTCSPEPDAFNPCEDIMGYDILRVLIWFISILAITGNILVLVILITSQYKLTVPRFLMCNLAFADLCIGIYLLLIASVDVHTKSWQTG-AGCDAAGFFTVFASELSVYTLTAITLERWHTITHAMLECKVHVRHAASIMLVGWVFAFAVALFPIFGISSY----KVSICLPM-----IDSPLSQLYVMSLLVLNVLAFVVICGCYTHIYLTVRNP------NITSSSSDTKIAKRMAMLIFTDFLCMAPISFFAISASLKVPLITVSKSKILLVLFYPINSCANPFLYAIFTRNFRRDFFILLSKFGCYEVQAQTYRSNFHPRNGHCPPA
>YN84_CAEEL
EVFHHISTTNKIFQKMFDKRNFSTDYTFNPKTFPGYRTYVASTYISFNVVGFVINAWVLYVVAPLLFVPKSILFYIFALCVGDLMTMIAMLLLVIELVFG--TWQFS-SMVCTSYLIFDSMNKFMAPMIVFLISRTCYSTVCLGEKAATLKYAIIQFCIAFAFVMILLWPVFAYSQVFTQEVVMRKCGFF----PPPQIEFWFNLIACITSYAVPLFGIIYWYVSVPFFLKRR---LVASSSMDAALRKVITTVLLLTVIYVLCWTPYWVSMFANRIWIMEKSIIIISYFIHLLPYISCVAYPLIFTLLNRGIRSAHAKIVADQRRRFRSLTDEASRTIPGTKMKKNE
>5H1D_CANFA
SPPNQSLEGLLQEASNRSLNATETPEAWGPETLQALKISLALLLSIITMATALSNAFVLTTIFLTRKLHTPANYLIGSLAMTDLLVSILVMPISIAYTTTRTWSFGQILCDIWLSSDITCCTASILHLCVIALDRYWAITDALYSKRRTAGRAAVMIATVWVISICISIPPLFWR----Q-E--SDCQVN-------TSQISYTIYSTCGAFYIPSVLLIILYGRIYVAARNRKLALERKRISAARERKATKTLGIILGAFIVCWLPFFVASLVLPICRASWLHPALFDFFTWLGYLNSLINPIIYTVFNEEFRQAFQRVVHVRKAS---------------------
>GASR_BOVIN
GASLCRSGGPLLNGSGTGNLSCEPPRIRGAGTRELELAIRVTLYAVIFLMSVGGNVLIIVVLGLSRRLRTVTNAFLLSLAVSDLLLAVACMPFTLLPNLMGTFIFGTVVCKAVSYFMGVSVSVSTLSLVAIALERYSAICRPLARVWQTRSHAARVIVATWMLSGLLMVPYPVYTAVQP---V-LQCMHR---WPSARVRQTWSVLLLLLLFFVPGVVMAVAYGLISRELYLGPGPTRPAQAKLLAKKRVVRMLLVIVVLFFLCWLPVYSANTWRAFDGPGALSGAPISFIHLLTYASACVNPLVYCFMHRRFRQACLDTCTRCCPRPPRARPRPLPSIASLSRLSYT
>BRS4_BOMOR
QTLPSAISSIAHLESLNDSFILGAKQSEDVSPGLEILALISVTYAVIISVGILGNTILIKVFFKIKSMQTVPNIFITSLAFGDLLLLLTCVPVDASRYIVDTWMFGRAGCKIISFIQLTSVGVSVFTLTVLSADRYRAIVKPLLQTSDAVLKTCGKAVCVWIISMLLAAPEAVFSDLYETT--FEACAPY---VSEKILQETHSLICFLVFYIVPLSIISAYYFLIAKTLYKSMPAHTHARKQIESRKRVAKTVLVLVALFAVCWLPNHMLYLYRSFTYHSAFHLSATIFARVLAFSNSCVNPFALYWLSRSFRQHFKKQVYCCKTEPPAS--QQSTGITAVKGNIQM
>OPSB_RAT
-----MSGE-EFYLFQNISSVGPWDGPQYHIAPVWAFHLQAAFMGFVFFAGTPLNATVLVATLHYKKLRQPLNYILVNVSLGGFLFCIFSVFTVFIASCHGYFLFGRHVCALEAFLGSVAGLVTGWSLAFLAFERYLVICKPFGNIRFNSKHALTVVLITWTIGIGVSIPPFFGWSRFIPEGLQCSCGPDWYTVGTKYRSEHYTWFLFIFCFIIPLSLICFSYFQLLRTLRAVAAQQQESATTQKAEREVSHMVVVMVGSFCLCYVPYAALAMYMVNNRNHGLYLRLVTIPAFFSKSSCVYNPIIYCFMNKQFRACILEMVCRKPMTD---ESDMSSTVSSSKVGPH-
>NK2R_RABIT
----MGACDIVTEANISSDIDSNATGVTAFSMPGWQLALWATAYLALVLVAVVGNATVIWIILAHRRMRTVTNYFIVNLALADLCMATFNAAFNFVYASHNIWYFGRAFCYFQNLFPITAMFVSIYSMTAIAADRYMAIVHPF-QPRLSGPGTKAVIAGIWLVALALAFPQCFYST---GA---TKCVVAWPEDSGGKMLLLYHLTVIALIYFLPLVVMFVAYSVIGFKLWRRPGHHGANLRHLRAKKKFVKTMVLVVVTFAVCWLPYHLYFLLGHFQDDIKFIQQVYLVLFWLAMSSTMYNPIIYCCLNHRFRSGFRLAFRCCPWVTPTEE----HTPSLSVRVNRC
>PAR2_MOUSE
---LETQPPITGKGVPVEPGFSIDEFSASILTGKLTTVFLPVVYIIVFVIGLPSNGMALWIFLFRTKKKHPAVIYMANLALADLLSVIWFPLKISYHLHGNNWVYGEALCKVLIGFFYGNMYCSILFMTCLSVQRYWVIVNPMGHPRKKANIAVGVSLAIWLLIFLVTIPLYVMKQTIY--N--TTCHDVLPEEVLVGDMFNYFLSLAIGVFLFPALLTASAYVLMIKTL------SAMDEHSEKKRQRAIRLIITVLAMYFICFAPSNLLLVVHYFLIKTSHVYALYLVALCLSTLNSCIDPFVYYFVSKDFRDHARNALLCRSVRTVNRMQISLKSGS--------
>OPSD_CARAU
MNGTEGDMFYVPMSNATGIVRSPYDYPQYYLVAPWAYACLAAYMFFLIITGFPVNFLTLYVTIEHKKLRTPLNYILLNLAISDLFMVFGGFTTTMYTSLHGYFVFGRVGCNPEGFFATLGGEMGLWSLVVLAFERWMVVCKPVSNFRFGENHAIMGVVFTWFMACTCAVPPLVGWSRYIPEGMQCSCGVDYYTRPQAYNNESFVIYMFIVHFIIPLIVIFFCYGRLVCTVKEAAAQHEESETTQRAEREVTRMVVIMVIGFLICWIPYASVAWYIFTHQGSEFGPVFMTLPAFFAKTAAVYNPCIYICMNKQFRHCMITTLCCGKNPFEEEEGASTASSVSSSSVSPA
>TRFR_HUMAN
------------MENETVSELNQTQLQPRAVVALEYQVVTILLVLIICGLGIVGNIMVVLVVMRTKHMRTPTNCYLVSLAVADLMVLVAAGLPNITDSIYGSWVYGYVGCLCITYLQYLGINASSCSITAFTIERYIAICHPIAQFLCTFSRAKKIIIFVWAFTSLYCMLWFFLLDLN--DA-SCGYKIS------RNYYSPIYLMDFGVFYVVPMILATVLYGFIARILFLNLNVNRCFNSTVSSRKQVTKMLAVVVILFALLWMPYRTLVVVNSFLSSPFQENWFLLFCRICIYLNSAINPVIYNLMSQKFRAAFRKLCNCKQKPTEKPANYSVKESDHFSTELDD
>PE24_MOUSE
GTIPRSNRELQRCVLLTTTIMSIPGVNASFSSTPERLNSPVTIPAVMFIFGVVGNLVAIVVLCKSRKKETTFYTLVCGLAVTDLLGTLLVSPVTIATYMKGQWPGDQALCDYSTFILLFFGLSGLSIICAMSIERYLAINHAYYSHYVDKRLAGLTLFAIYASNVLFCALPNMGLGRSERQYPGTWCFID---WTTNVTAYAAFSYMYAGFSSFLILATVLCNVLVCGALLRMAAAVASFRRIAGAEIQMVILLIATSLVVLICSIPLVVRVFINQLYQPNDISRNPDLQAIRIASVNPILDPWIYILLRKTVLSKAIEKIKCLFCRIGGSGRDSSRRTSSAMSGHSR
>C3AR_CAVPO
--------------MESSSAETNSTGLHLEPQYQPETILAMAILGLTFVLGLPGNGLVLWVAGLKMR-RTVNTVWFLHLTVADFVCCLSLPFSMAHLALRGYWPYGEILCKFIPTVIIFNMFASVFLLTAISLDRCLMVLKPICQNHRNVRTACIICGCIWLVAFVLCIPVFVYRETFT-EE-DDLSPFT-HEYRTPRLLKVITFTRLVVGFLLPMIIMVACYTLIIFRM--------RRVRVVKSWNKALHLAMVVVTIFLICWAPYHVFGVLILFINPEAALLSWDHVSIALASANSCFNPFLYALLGRDLRKRVRQSMKGILEAAFSEDISKSAFS---------
>GPR3_MOUSE
GAGSSMAWFSAGSGSVNVSSVDPVEEPTGPATLLPSPRAWDVVLCISGTLVSCENALVVAIIVGTPAFRAPMFLLVGSLAVADLLAG-LGLVLHFAAD-F--CIGSPEMSLMLVGVLAMAFTASIGSLLAITVDRYLSLYNALYYSETTVTRTYVMLALVWVGALGLGLVPVLAWNCR-D--GLTTCGVV-------YPLSKNHLVVLAIAFFMVFGIMLQLYAQICRIVCRHIALHLLPASHYVATRKGIATLAVVLGAFAACWLPFTVYCLLGDA----DSPRLYTYLTLLPATYNSMINPVIYAFRNQDVQKVLWAICCCCSTSKIPFRSRSP------------
>OPSD_OCTDO
--MVESTTLVNQTWWYNPTVDIHPHWAKFDPIPDAVYYSVGIFIGVVGIIGILGNGVVIYLFSKTKSLQTPANMFIINLAMSDLSFSAINGFPLKTISAFMKWIFGKVACQLYGLLGGIFGFMSINTMAMISIDRYNVIGRPMASKKMSHRRAFLMIIFVWMWSIVWSVGPVFNWGAYVPEGILTSCSFD--YLSTDPSTRSFILCMYFCGFMLPIIIIAFCYFNIVMSVSNHRLNLRKAQAGASAEMKLAKISMVIITQFMLSWSPYAIIALLAQFGPAEWVTPYAAELPVLFAKASAIHNPIVYSVSHPKFREAIQTTFPWLLTCCQFDEKECEEVVASERG-GES
>B1AR_MELGA
WLPPDCGPHNRSGGGGATAAPTGSRQVSAELLSQQWEAGMSLLMALVVLLIVAGNVLVIAAIGRTQRLQTLTNLFITSLACADLVMGLLVVPFGATLVVRGTWLWGSFLCECWTSLDVLCVTASIETLCVIAIDRYLAITSPFYQSLMTRARAKVIICTVWAISALVSFLPIMMHWWRDL-K--GCCDFV--------TNRAYAIASSIISFYIPLLIMIFVYLRVYREAKEQNGRRKTSRVMAMREHKALKTLGIIMGVFTLCWLPFFLVNIVNVFNRDL-VPDWLFVFFNWLGYANSAFNPIIYCRSP-DFRKAFKRLLCFPRKADRRLHAGGQFISTLGSPEHSP
>5H1B_CAVPO
PAVLGSQTGLPHANVSAPPNNAP-SHIYQDSIALPWKVLLVVLLALITLATTLSNAFVIATVYRTRKLHTPANYLIASLAFTDLLVSILVMPISTMYTVTGRWTLGQALCDFWLSSDITCCTASIMHLCVIALDRYWAITDAVYSAKRTPRRAAGMIALVWVFSICISLPPFFWR----E-E--LDCLVN-------TDHVLYTVYSTGGAFYLPTLLLIALYGRIYVEARSRRVSLEKKKLMAARERKATKTLGVILGAFIVCWLPFFIISLVMPICKDAWFHMAIFDFFTWLGYLNSLINPIIYTMSNEDFKQAFHKLIRFKCTT---------------------
>BRB2_MOUSE
IEMFNVTTQVLGSALNGTLSKDNC---PDTEWWSWLNAIQAPFLWVLFLLAALENLFVLSVFFLHKNSCTVAEIYLGNLAAADLILACGLPFWAITIANNFDWVFGEVLCRVVNTMIYMNLYSSICFLMLVSIDRYLALVKTMMGRMRGVRWAKLYSLVIWGCTLLLSSPMLVFRTMRE-HN--TACVIV---YPSRSWEVFTNVLLNLVGFLLPLSVITFCTVRILQVLRNN---EMKKFKEVQTERKATVLVLAVLGLFVLCWVPFQISTFLDTLLRLGHAVDVITQISSYVAYSNSGLNPLVYVIVGKRFRKKSREVYRVLCQKGGCMGEPVQLRTS-------I
>OPSB_APIME
YVPSMREKFLGWNVPPEYSDLVRPHWRAFPAPGKHFHIGLAIIYSMLLIMSLVGNCCVIWIFSTSKSLRTPSNMFIVSLAIFDIIMAFEMPMLVISSFMERM--GWEIGCDVYSVFGSISGMGQAMTNAAIAFDRYRTISCPI-DGRLNSKQAAVIIAFTWFWVTPFTVLPLLKVWGRYTEGFLTTCSFD--FLTDDEDTKVFVTCIFIWAYVIPLIFIILFYSRLLSSIRNHNVKSNQDKER-SAEVRIAKVAFTIFFLFLLAWTPYATVALIGVYGNRELLTPVSTMLPAVFAKTVSCIDPWIYAINHPRYRQELQKRCKWMGIHEPETTSDATKTDE--------
>O.YR_HUMAN
GALAANWSAEAANASAAPPGAEGNRTAGPPRRNEALARVEVAVLCLILLLALSGNACVLLALRTTRQKHSRLFFFMKHLSIADLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPL--RSLRRRTDRLAVLATWLGCLVASAPQVHIFSLRE----VFDCWAV---FIQPWGPKAYITWITLAVYIVPVIVLATCYGLISFKIWQNRVAVSSVKLISKAKIRTVKMTFIIVLAFIVCWTPFFFVQMWSVWDANA-KEASAFIIVMLLASLNSCCNPWIYMLFTGHLFHELVQRFLCCSASYLKGRRLGENSSSFVLSHRSS
>OPSR_CAPHI
ANFEESTQGSIFTYTNSNSTRDPFEGPNYHIAPRWVYHLTSAWMVFVVIASVFTNGLVLAATMRFKKLRHPLNWILVNLAIADLAETIIASTISVVNQMYGYFVLGHPLCVVEGYTVSLCGITGLWSLAIISWERWMVVCKPFGNVRFDAKLATAGIAFSWIWAAVWTAPPIFGWSRYWPHGLKTSCGPDVFSGSSYPGVQSYMIVLMITCCFIPLSVIILCYLQVWLAIRAVAKQQKESESTQKAEKEVTRMVMVMIFAYCLCWGPYTFFACFAAAHPGYAFHPLVAALPAYFAKSATIYNPIIYVFMNRQFRNCILQLFGKKVDDS-----SELASSV--SSVSPA
>APJ_HUMAN
---------MEEGGDFDNYYGADNQSECEYTDWKSSGALIPAIYMLVFLLGTTGNGLVLWTVFRSSRKRRSADIFIASLAVADLTFVVTLPLWATYTYRDYDWPFGTFFCKLSSYLIFVNMYASVFCLTGLSFDRYLAIVRPVNARLRLRVSGAVATAVLWVLAALLAMPVMVLRTTGDQCY-MDYSMVA-TVSSEWAWEVGLGVSSTTVGFVVPFTIMLTCYFFIAQTIAGHFR--KERIEGLRKRRRLLSIIVVLVVTFALCWMPYHLVKTLYMLGSLLLFLMNIFPYCTCISYVNSCLNPFLYAFFDPRFRQACTSMLCCGQSRCAGTSHSSSSSGHSQGPGPNM
>NY1R_PIG
TLSSQVENHSIYYNFSEKNSQFLAFENDDCHLPLAMIFTLALAYGAVIILGVSGNLALIIIILKQKEMRNVTNILIVNLSFSDLLVAIMCLPFTFVYTLMDHWVFGEVMCKLNPFVQCVSITVSIFSLVLIAVERHQLIINPR-GWRPSNRHAYVGIAVIWVLAVASSLPFLIYQVLTDKD--KYVCFDK---FLSDSHRLSYTTLLLVLQYFGPLCFIFICYFKIYIRLKRRNMMMRDNKYRSSETKRINVMLLSIVVAFAVCWLPLTIFNTVFDWNHQICNHNLLFLLCHLTAMISTCINPIFYGFLNKNFQRDLQFFFNFCDFRSRDDDY---TMHTDVSKTSLK
>OAR1_LYMST
---------MSRDIFMKRLRLHLLFDEVAMVTHIVGDVLSSVLLCAVVLLVLVGNTLVVAAVATSRKLRTVTNVFIVNLACADLLLGVLVLPFSAVNEIKDVWIFGHVWCQVWLAVDVWLCTASILNLCCISLDRYLAITRPIYPGLMSAKRAKTLVAGVWLFSFVICCPPLIGWNDGGT-Y--TTCELT--------NSRGYRIYAALGSFFIPMLVMVFFYLQIYRAAVKTHKPMRLHMQKFNREKKAAKTLAIIVGAFIMCWMPFFTIYLVGAFCENC-ISPIVFSVAFWLGYCNSAMNPCVYALFSRDFRFAFRKLLTCSCKAWSKNRSFRPIQLHCATQDDAK
>OLF5_CHICK
-------------MALGNCTTPTTFILSGLTDNPRLQMPLFMVFLAIYTITLLANLGLIALISVDFHLQTPMYIFLQNLSFTDAAYSTVITPKMLATFLEERRTISYVGCILQYFSFVLLTSSECLLLAVMAYDRYVAICKPLYPAIMTKAVCWRLVEGLYSLAFLNSLVHTSGLLKL---SNHFFCDNSQISSSSTTLNELLVFIFGSWFAMSSIITTPISYVFIILTV--------VRIRSKDGKYKAFSTCTSHLMAVSLFHGTVIFMYLRPVKLF---SLDTDKIASLFYTVVIPMLNPLIYSWRNKEVKDALRRVIATNVWIH--------------------
>A2AR_CARAU
----------MDVTQSNATKDDANITVTPWPYTETAAAFIILVVSVIILVSIVGNVLVIVAVLTSRALRAPQNLFLVSLACADILVATLVIPFSLANEIMGYWFFGSTWCAFYLALDVLFCTSSIVHLCAISLDRYWSVTKAVYNLKRTPKRIKSMIAVVWVISAVISFPPLIMTKH-------KECLIN--------DETWYILSSSLVSFFAPGFIMITVYCKIYRVAKQRSKQASKTKVAQMREKRFTFVLTVVMGVFVLCWFPFFFTYSLHAICGDSEPPEALFKLFFWIGYCNSSVNPIIYTIFNRDFRKAFKKICLLDCAAHLRDSCLGTCIFECHQKSNQE
>VQ3L_CAPVK
SNYTTAYNTTYYSDDYDDYEVSIVDIPHCDDGVDTTSFGLITLYSTIFFLGLFGNIIVLTVLRKYKI-KTIQDMFLLNLTLSDLIFVLVFPFNLYDSIAKQ-WSLGDCLCKFKAMFYFVGFYNSMSFITLMSIDRYLAVVHPVSMPIRTKRYGIVLSMVVWIVSTIESFPIMLFYETKK-VY-ITYCHVF-YNDNAKIWKLFINFEINIFGMIIPLTILLYCYYKILNTL----------KTSQTKNKKAIKMVFLIVICSVLFLLPFSVTVFVSSLYLLNRFVNLAVHVAEIVSLCHCFINPLIYAFCSREFTKKLLRLRTTSSAGSISIG----------------
>GP40_HUMAN
--------------------------------MDLPPQLSFGLYVAAFALGFPLNVLAIRGATAHARRLTPSLVYALNLGCSDLLLTVSLPLKAVEALASGAWPLPASLCPVFAVAHFFPLYAGGGFLAALSAGRYLGAAFPLYQAFRRPCYSWGVCAAIWALVLCHLGLVFGLEAPGGTPVGSPVCLEA----WDPASAGPARFSLSLLLFFLPLAITAFCYVGCLRAL-------ARSGLTHRRKLRAAWVAGGALLTLLLCVGPYNASNVASFLYP--NLGGSWRKLGLITGAWSVVLNPLVTGYLGRGPGLKTVCAARTQGGKSQK------------------
>CKR8_MACMU
--MDYTLDPSMTTMTDYYYPDSLSSPCDGELIQRNDKLLLAVFYCLLFVFSLLGNSLVILVLVVCKKLRNITDIYLLNLALSDLLFVFSFPFQTYYQLDQ--WVFGTVMCKVVSGFYYIGFYSSMFFITLMSVDRYLAVVHAVIKVRTIRMGTTTLSLLVWLTAIMATIPLLVFYQVAS-ED-VLQCYSF-YNQQTLKWKIFTNFEMNILGLLIPFTIFMFCYIKILHQL---------KRCQNHNKTKAIRLVLIVVIASLLFWVPFNVVLFLTSLHSMHQQLNYATHVTEIISFTHCCVNPVIYAFVGEKFKKHLSEIFQKSCSHIFIYLGRQMSSSCQQHSFRSS
>NTR1_RAT
EATFLALSLSNGSGNTSESDTAGPNSDLDVNTDIYSKVLVTAIYLALFVVGTVGNSVTAFTLARKKSLQSTVHYHLGSLALSDLLILLLAMPVELYNFIWVHWAFGDAGCRGYYFLRDACTYATALNVASLSVERYLAICHPFAKTLMSRSRTKKFISAIWLASALLAIPMLFTMGLQN--SGGLVCTPI---VDTATVKVVIQVN-TFMSFLFPMLVISILNTVIANKLTVMEHSMTIEPGRVQALRHGVLVLRAVVIAFVVCWLPYHVRRLMFCYISDEDFYHYFYMLTNALFYVSSAINPILYNLVSANFRQVFLSTLACLCPGWRHRRKKRPSMSSNHAFSTSA
>B3AR_MACMU
MAPWPHGNSSLVPWPDVPTLAPNTANTSGLPGVPWAAALAGALLALAVLATVGGNLLVIVAITRTPRLQTMTNVFVTSLAAADLVMGLLVVPPAATLVLTGHWPLGATGCELWTSVDVLCVTASIETLCALAVDRYLAVTNPLYGALVTKRRARAAVVLVWVVSAAVSFAPIMSQWWRVQ-R--RCCAFA--------SNMPYVLLSSSVSFYLPLLVMLFVYARVFVVATRQGVPRRPARLLPLREHRALCTLGLIMGTFTLCWLPFFLANVLRALGGPS-VPDPAFLALNWLGYANSAFNPLIYCRSP-DFRSAFRRLLCHCGGRLPREPCAADAPLRPGPAPRSP
>GALS_MOUSE
------------MNGSDSQGAEDSSQE-GGGGWQPEAVLVPLFFALIFLVGAVGNALVLAVLLRGGQAVSTTNLFILNLGVADLCFILCCVPFQATIYTLDDWVFGSLLCKAVHFLIFLTMHASSFTLAAVSLDRYLAIRYPMSRELRTPRNALAAIGLIWGLALLFSGPYLSYYS---AN--LTVCHPA----WSAPRRRAMDLCTFVFSYLLPVLVLSLTYARTLHYLWRTDP-VAAGSGSQRAKRKVTRMIVIVAVLFCLCWMPHHALILCVWFGRFPRATYALRILSHLVSYANSCVNPIVYALVSKHFRKGFRKICAGLLRRAPRRASGRVHSGGMLEPESTD
>P2Y3_MELGA
----------------MSMANFTAGRNSCTFQEEFKQVLLPLVYSVVFLLGLPLNAVVIGQIWLARKALTRTTIYMLNLATADLLYVCSLPLLIYNYTQKDYWPFGDFTCKFVRFQFYTNLHGSILFLTCISVQRYMGICHPLWHKKKGKKLTWLVCAAVWFIVIAQCLPTFVFASTG-----RTVCYDL-SPPDRSASYFPYGITLTITGFLLPFAAILACYCSMARILCQK---ELIGLAVHKKKDKAVRMIIIVVIVFSISFFPFHLTKTIYLIVRSSQAFAIAYKCTRPFASMNSVLDPILFYFTQRKFRESTRYLLDKMSSKWRHDHCITY------------
>OPSD_TURTR
MNGTEGLNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSVLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVANLFMVFGGFTTTLYTSLHAYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLALTWIMAMACAAAPLVGWSRYIPEGMQCSCGIDYYTSRQEVNNESFVIYMFVVHFTIPLVIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVVAFLICWVPYASVAFYIFTHQGSDFGPIFMTIPSFFAKSSSIYNPVIYIMMNKQFRNCMLTTLCCGRNPLGDDEASTTASKTETSQVAPA
>A2AB_CAVPO
-------------------------MDHQEPYSVQATAAIAAVITFLILFTIFGNALVILAVLTSRSLPAPQNLFLVSLAAADILVATLIIPFSLANELLGYWYFWRTWCEVYLALDVLFCTSSIVHLCAISLDRYWAVSRALYNSKRTPRRIKCIILTVWLIAAVISLPPLIYKGD-Q-----PQCKIN--------QEAWYILASSIGSFFAPCLIMILVYLRIYLIAKRSGAVWWRRRTQMTREKRFTFVLAVVIGVFVLCWFPFFFTYSLGAICPQHKVPHGLFQFFFWIGYCNSSLNPVIYTIFNQDFRRAFRRILCRQWTQTAW------------------
>CML1_MOUSE
EYDAYNDSGIYDDEYSDGFGYFVDLEEASPWEAKVAPVFLVVIYSLVCFLGLLGNGLVIVIATFKMK-KTVNTVWFVNLAVADFLFNIFLPMHITYAAMDYHWVFGKAMCKISNFLLSHNMYTSVFLLTVISFDRCISVLLPVSQNHRSIRLAYMTCSAVWVLAFFLSSPSLVFRDTAN-SS--HPAHSQ-VVSTGYSRHVAVTVTRFLCGFLIPVFIITACYLTIVFKL---------QRNRLAKNKKPFKIIITIIITFFLCWCPYHTLYLLELHHTAVSVFSLGLPLATAVAIANSCMNPILYVFMGHDFRK-FKVALFSRLANALSEDTGPSFTKMSSLNEKAS
>5HT1_APLCA
-MKSLKSSTHDVPHPEHVVWAPPAYDEQHHLFFSHGTVLIGIVGSLIITVAVVGNVLVCLAIFTEPISHSKSNFFIVSLAVADLLLALLVMTFALVNDMYGYWLFGETFCFIWMSADVMCETASIFSICVISYDRLKQVQKPLYEEFMTTTRALLIIACLWICSFVLSFVPIFLEWHELGD-AKHVCLFD--------VHFTYSVIYSFICFYVPCTLMLTNYLRLFLIAQTHQLRASSYRNQGTQGSKAARTLTIITGTFLACWLPFFIINPIAAADEHL-IPLECFMVTIWLGYFNSSVNPIIYGTSNSKFRAAFKRLLRCRSVKSVVGSISPVSWIRPSRLDLSS
>5H1A_FUGRU
NDSNATSGYSDTAAVDWDEGENATGSGSLPDPELSYQIITSLFLGALILCSIFGNSCVVAAIALERSLQNVANYLIGSLAVTDLMVSVLVLPMAALYQVLNKWTLGQDICDLFIALDVLCCTSSILHLCAIALDRYWAITDPIYVNKRTPRRAAVLISVTWLIGFSISIPPMLGWRS---NP--DACIIS--------QDPGYTIYSTFGAFYIPLILMLVLYGRIFKAARFRINEGTRRKIALARERKTVKTLGIIMGTFIFCWLPFFIVALVLPFCAENYMPEWLGAVINWLGYSNSLLNPIIYAYFNKDFQSAFKKILRCKFHRH--------------------
>ACTR_BOVIN
-------------MKHILNLYENINSTARNNSDCPAVILPEEIFFTVSIVGVLENLMVLLAVAKNKSLQSPMYFFICSLAISDMLGSLYKILENVLIMFKNMGSFESTADDVVDSLFILSLLGSICSLSVIAADRYITIFHALYHRIMTPHRALVILTVLWAGCTGSGITIVTFFS----------------------HHHVPTVIAFTALFPLMLAFILCLYVHMFLLARSH---TRRTPSLPKANMRGAVTLTVLLGVFIFCWAPFVLHVLLMTFCPADACYMSLFQVNGVLIMCNAIIDPFIYAFRSPELRVAFKKMVICNCYQ---------------------
>CKR9_MOUSE
DDFSYDSTASTDDYMNLNFSSFFC---KKNNVRQFASHFLPPLYWLVFIVGTLGNSLVILVYWYCTRVKTMTDMFLLNLAIADLLFLATLPFWAIAAAGQ--WMFQTFMCKVVNSMYKMNFYSCVLLIMCISVDRYIAIVQAMVWRQKRLLYSKMVCITIWVMAAVLCTPEILYSQVSG-G--IATCTMVYPKDKNAKLKSAVLILKVTLGFFLPFMVMAFCYTIIIHTL---------VQAKKSSKHKALKVTITVLTVFIMSQFPYNSILVVQAVDAYATNIDICFQVTQTIAFFHSCLNPVLYVFVGERFRRDLVKTLKNLGCISQAQWVSFT--------GSLK
>5H1B_SPAEH
CAPPPPAGSQTQTPSSNLSHNSADSYIYQDSIALPWKVLLVALLALITLATTLSNAFVIATVYRTRKLHTPANYLIASLAVTDLLVSILVMPISTMYTVTGRWTLGQVVCDFWLSSDITCCTASIMHLCVIALDRYWAITDAVYSAKRTPRRAAVMIALVWVFSISISLPRFFWR----E-E--LDCLVN-------TDHVLYTVYSTVGAFYLPTLLLIALYGRIYVEARSRRVSLEKKKLMAARERKATKTLGIILGAFIVCWLPFFIISLVMPICKDAWFHMAIFDFFNWLGYLNSLINPIIYTMPNEDFKQAFHKLIRFKCTG---------------------
>OPS4_DROVI
VSGNGDLQFLGWNVPPDQIQHIPEHWLTQLEPPASMHYMLGVFYIFLFCASTVGNGMVIWIFSTSKALRTPSNMFVLNLAVFDFIMCLKAPIFIYNSFHRG-FALGNTGCQIFAAIGSYSGIGAGMTNAAIGYDRLNVITKPM-NRNMTFTKAIIMNVIIWLYCTPWVVLPLTQFWDRFPEGYLTSCTFD--YLTDNFDTRLFVGTIFFFSFVCPTLMIIYYYSQIVGHVFSHNVESNVDKSKDTAEIRIAKAAITICFLFFVSWTPYGVMSLIGAFGDKSLLTPGATMIPACTCKLVACIDPFVYAISHPRYRMELQKRCPWLAIDEKAPESSSAEQQQTTAA----
>OPSG_RABIT
ESHEDSTQASIFTYTNSNSTRGPFEGPNFHIAPRWVYHLTSAWMILVVIASVFTNGLVLVATMRFKKLRHPLNWILVNLAVADLAETVIASTISVVNQFYGYFVLGHPLCVVEGYTVSLCGITGLWSLAIISWERWLVVCKPFGNVRFDAKLAIAGIAFSWIWAAVWTAPPIFGWSRYWPYGLKTSCGPDVFSGTSYPGVQSYMMVLMVTCCIIPLSVIVLCYLQVWMAIRTVAKQQKESESTQKAEKEVTRMVVVMVFAYCLCWGPYTFFACFATAHPGYSFHPLVAAIPSYFAKSATIYNPIIYVFMNRQFRNCILQLFGKKVEDS-----SELASSV--SSVSPA
>NK3R_RABIT
LVQAGNLSSSLPSSVPGLPTTPRANLTNQFVQPSWRIALWSLAYGVVVAVAVFGNLIVIWIILAHKRMRTVTNYFLVNLAFSDASMAAFNTLVNFIYALHSEWYFGANYCRFQNFFPITAVFASIYSMTAIAVDRYMAIIDPL-KPRLSATATKIVIGSIWILAFLLALPQCLYSK---GR---TLCYVQ--WPEGPKQHFIYHIIVIILVYCFPLLIMGITYTIVGITLWGG--PCDKYHEQLKAKRKVVKMMIIVVVTFAICWLPYHIYFILTAIYQQLKYIQQVYLASFWLAMSSTMYNPIIYCCLNKRFRAGFKRAFRWCPFIQVSSYDELEPTRQSSLYTVTR
>ACM5_HUMAN
-------MEGDSYHNATTVNGTPVNHQPLERHRLWEVITIAAVTAVVSLITIVGNVLVMISFKVNSQLKTVNNYYLLSLACADLIIGIFSMNLYTTYILMGRWALGSLACDLWLALDYVASNASVMNLLVISFDRYFSITRPLYRAKRTPKRAGIMIGLAWLISFILWAPAILCWQYLV-VP-LDECQIQ------FLSEPTITFGTAIAAFYIPVSVMTILYCRIYRETEKRNPSTKRKRVVLVKERKAAQTLSAILLAFIITWTPYNIMVLVSTFCDKC-VPVTLWHLGYWLCYVNSTVNPICYALCNRTFRKTFKMLLLCRWKKKKVEEKLYW------------
>O3A1_HUMAN
----------MQPESGANGTVIAEFILLGLLEAPGLQPVVFVLFLFAYLVTVRGNLSILAAVLVEPKLHTPMYFFLGNLSVLDVGCISVTVPSMLSRLLSRKRAVPCGACLTQLFFFHLFVGVDCFLLTAMAYDQFLAICRPLYSTRMSQTVQRMLVAASWACAFTNALTHTVAMSTL---NNHFYCDLPQLSCSSTQLNELLLFAVGFIMAGTPMALIVISYIHVAAAV--------LRIRSVEGRKKAFSTCGSHLTVVAIFYGSGIFNYMRLGSTK---LSDKDKAVGIFNTVINPMLNPIIYSFRNPDVQSAIWRMLTGRRSLA--------------------
>CKRB_BOVIN
YNQSTDYYYEENEMNDTHDYSQYEVICIKEEVRKFAKVFLPAFFTIAFIIGLAGNSTVVAIYAYYKKRRTKTDVYILNLAVADLFLLFTLPFWAVNAVHG--WVLGKIMCKVTSALYTVNFVSGMQFLACISTDRYWAVTKAP-SQSGVGKPCWVICFCVWVAAILLSIPQLVFYTVN----HKARCVPI--YHLGTSMKASIQILEICIGFIIPFLIMAVCYFITAKTL---------IKMPNIKKSQPLKVLFTVVIVFIVTQLPYNIVKFCQAIDIIYKRMDVAIQITESIALFHSCLNPVLYVFMGTSFKNYIMKVAKKYGSW---------NVEEIPFESEDA
>TA2R_HUMAN
--MWPNG-----------SSLGPCFRPTNITLEERRLIASPWFAASFCVVGLASNLLALSVLAGARQTRSSFLTFLCGLVLTDFLGLLVTGTIVVSQHAALFVDPGCRLCRFMGVVMIFFGLSPLLLGAAMASERYLGITRPFRPAVASQRRAWATVGLVWAAALALGLLPLLGVGRYTVQYPGSWCFLT----LGAESGDVAFGLLFSMLGGLSVGLSFLLNTVSVATLCHVYHGEAAQQRPRDSEVEMMAQLLGIMVVASVCWLPLLVFIAQTVLRNPPRTTEKELLIYLRVATWNQILDPWVYILFRRAVLRRLQPRLSTRPRRVSLCGPAWSTATSASRVQAIL
>AG2R_CHICK
------MVPNYSTEETVKRIHVDC---PVSGRHSYIYIMVPTVYSIIFIIGIFGNSLVVIVIYCYMKLKTVASIFLLNLALADLCFLITLPLWAAYTAMEYQWPFGNCLCKLASAGISFNLYASVFLLTCLSIDRYLAIVHPVSRIRRTMFVARVTCIVIWLLAGVASLPVIIHRNIFF-LN--TVCGFR-YDNNNTTLRVGLGLSKNLLGFLIPFLIILTSYTLIWKTLKKA----YQIQRNKTRNDDIFKMIVAIVFFFFFSWIPHQVFTFLDVLIQLHDIVDTAMPFTICIAYFNNCLNPFFYVFFGKNFKKYFLQLIKYIPPNVSTHPSLTTRPPE-------N
>P2YR_MOUSE
AAFLAGLGSLWGNSTVASTAAVSSSFQCALTKTGFQFYYLPAVYILVFIIGFLGNSVAIWMFVFHMKPWSGISVYMFNLALADFLYVLTLPALIFYYFNKTDWIFGDAMCKLQRFIFHVNLYGSILFLTCISAHRYSGVVYPLSLGRLKKKNAIYVSVLVWLIVVVAISPILFYSGTG----KTVTCYDT-TSNDYLRSYFIYSMCTTVAMFCIPLVLILGCYGLIVKALIY------NDLDNSPLRRKSIYLVIIVLTVFAVSYIPFHVMKTMNLRARLDDRVYATYQVTRGLASLNSCVDPILYFLAGDTFRRRLSRATRKASRRSEANLQSKSLSEFKQNGDTSL
>NY5R_RAT
-------------NTAAARNAAFPAWEDYRGSVDDLQYFLIGLYTFVSLLGFMGNLLILMAVMKKRNQKTTVNFLIGNLAFSDILVVLFCSPFTLTSVLLDQWMFGKAMCHIMPFLQCVSVLVSTLILISIAIVRYHMIKHPI-SNNLTANHGYFLIATVWTLGFAICSPLPVFHSLVESS--KYLCVES---WPSDSYRIAFTISLLLVQYILPLVCLTVSHTSVCRSISCGAHEKRSITRIKKRSRSVFYRLTILILVFAVSWMPLHVFHVVTDFNDNLRHFKLVYCICHLLGMMSCCLNPILYGFLNNGIKADLRALIHCLHMS---------------------
>THRR_PAPHA
LTEYRLVSINKSSPLQKPLPAFISEDASGYLTSSWLTLFVPSVYTGVFVVSLPVNIMAIVVFILKMKVKKPAVVYMLHLATADVLFVSVLPFKISYYLSGSDWQFGSELCRFVTAAFYCNMYASILLMTVISIDRFLAVVYPMSLSWRTLGRASFTCLAIWALAIAGVVPLLLKEQTIQ--N--TTCHDVLNETLLEGYYAYYFSAFSAVFFFVPLIISTVCYVSIIRCL------SSSTVANRSKKSRALFLSAAVFCIFIICFGPTNILLIAHYSFLSHEAAYFAYLLCVCVSSISCCIDPLIYYYASSECQRYVYSILCCKESSDPSSSNSSGDTCS--------
>FSHR_HUMAN
FDMTYTEFDYDLCNEVVDVTCSPKPDAFNPCEDIMGYNILRVLIWFISILAITGNIIVLVILTTSQYKLTVPRFLMCNLAFADLCIGIYLLLIASVDIHTKSWQTG-AGCDAAGFFTVFASELSVYTLTAITLERWHTITHAMLDCKVQLRHAASVMVMGWIFAFAAALFPIFGISSY----KVSICLPM-----IDSPLSQLYVMSLLVLNVLAFVVICGCYIHIYLTVRNP------NIVSSSSDTRIAKRMAMLIFTDFLCMAPISFFAISASLKVPLITVSKAKILLVLFHPINSCANPFLYAIFTKNFRRDFFILLSKCGCYEMQAQIYRTNTHPRNGHCSSA
>OPSB_CONCO
MNGTEGPNFYVPMSNATGVVRSPFEYPQYYLAEPWAFSILAAYMFFLIITGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVFGETGCNLEGYFATLGGEISLWSLVVLAIERWVVVCKPISNFRFGENHAIMGLTLTWVMANACAMPPLFGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFLVHFTIPLTIISFCYGRLVCAVKEAAAQQQESETTQRAEREVTRMVVIMVISFLVCWIPYASVAWYIFTHQGSTFGPIFMTVPSFFAKSSSIYNPMIYICMNKQFRNCMITTLFCGKNPFEGEE--EGASAVSS--VSPA
>A2AA_HUMAN
----MGSLQPDAGNASWNGTEAPGGGARATPYSLQVTLTLVCLAGLLMLLTVFGNVLVIIAVFTSRALKAPQNLFLVSLASADILVATLVIPFSLANEVMGYWYFGKAWCEIYLALDVLFCTSSIVHLCAISLDRYWSITQAIYNLKRTPRRIKAIIITVWVISAVISFPPLISIEKKGQ-P--PRCEIN--------DQKWYVISSCIGSFFAPCLIMILVYVRIYQIAKRRRVGASRWRGRQNREKRFTFVLAVVIGVFVVCWFPFFFTYTLTAVGCS--VPRTLFKFFFWFGYCNSSLNPVIYTIFNHDFRRAFKKILCRGDRKRIV------------------
>OPSD_TRIMA
MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNVEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTLPAFFAKSASIYNPVIYIMMNKQFRNCMLTTICCGKNPFAEEEGATTVSKTETSQVAPA
>AA3R_HUMAN
-------------------------MPNNSTALSLANVTYITMEIFIGLCAIVGNVLVICVVKLNPSLQTTTFYFIVSLALADIAVGVLVMPLAIVVSLG--ITIHFYSCLFMTCLLLIFTHASIMSLLAIAVDRYLRVKLTVYKRVTTHRRIWLALGLCWLVSFLVGLTPMFGWN-MKE-Y-FLSCQFV-----SVMRMDYMVYFSFLTWIFIPLVVMCAIYLDIFYIIRNK--NSKETGAFYGREFKTAKSLFLVLFLFALSWLPLSIINCIIYFNGE--VPQLVLYMGILLSHANSMMNPIVYAYKIKKFKETYLLILKACVVCHPSDSLDTS------------
>ET1R_HUMAN
ELSFLVTTHQPTNLVLPSNGSMHNYCPQQTKITSAFKYINTVISCTIFIVGMVGNATLLRIIYQNKCMRNGPNALIASLALGDLIYVVIDLPINVFKLLAGRNDFGVFLCKLFPFLQKSSVGITVLNLCALSVDRYRAVASWSVQGIGIPLVTAIEIVSIWILSFILAIPEAIGFVMVPKT--HKTCMLNATSKFMEFYQDVKDWWLFGFYFCMPLVCTAIFYTLMTCEMLNRGSL-IALSEHLKQRREVAKTVFCLVVIFALCWFPLHLSRILKKTVYNESFLLLMDYIGINLATMNSCINPIALYFVSKKFKNCFQSCLCCCCYQSKSLMTSVPWKNHDQNNHNTD
>GALS_HUMAN
------------MNVSGCPGAGNASQAGGGGGWHPEAVIVPLLFALIFLVGTVGNTLVLAVLLRGGQAVSTTNLFILNLGVADLCFILCCVPFQATIYTLDGWVFGSLLCKAVHFLIFLTMHASSFTLAAVSLDRYLAIRYPLSRELRTPRNALAAIGLIWGLSLLFSGPYLSYYR---AN--LTVCHPA----WSAPRRRAMDICTFVFSYLLPVLVLGLTYARTLRYLWRADP-VAAGSGARRAKRKVTRMILIVAALFCLCWMPHHALILCVWFGQFPRATYALRILSHLVSYANSCVNPIVYALVSKHFRKGFRTICAGLLGRAPGRASGRVHSGSVLERESSD
>OPSV_ORYLA
-------MGKYFYLYENISKVGPYDGPQYYLAPTWAFYLQAAFMGFVFFVGTPLNFVVLLATAKYKKLRVPLNYILVNITFAGFIFVTFSVSQVFLASVRGYYFFGQTLCALEAAVGAVAGLVTSWSLAVLSFERYLVICKPFGAFKFGSNHALAAVIFTWFMGV-VRCPPFFGWSRYIPEGLGCSCGPDWYTNCEEFSCASYSKFLLVTCFICPITIIIFSYSQLLGALRAVAAQQAESASTQKAEKEVSRMIIVMVASFVTCYGPYALTAQYYAYSQDENKDYRLVTIPAFFSKSSCVYNPLIYAFMNKQFNGCIMEMVFGKKMEE-----ASESTDS--------
>YWO1_CAEEL
----------------------------MNHVDYVAHVIVMPIVLSIGMINQCLNVCTLLHIRT------SIFLYLKASAIADILSIVAFIPFLFRHAKLIDELGMFYHAHLELPLINALISASALNIVAMTVDRYVSVCHPINETKPSRRRTMLIIVMIYFIALMIYFPSVFQKKLG-VTTIYTIVRNE---VEALQVFKFYLIVRECICRWGPVLLLVILNMCVVRGLRKIRTEQRQLRSPRDDRSRISVLLFVTSATFIICNIPASVISFFVRRVSGSLFWQIFRAIANLLQVTSYLYNFYLYALCSSEYRHAFLRLFGCRSSLSPTSTGDSPGKRCHQAVVLLG
>CKR5_RAT
--MDFQGSIPTYIYDIDYSMSAPC---QKVNVKQIAAQLLPPLYSLVFIFGFVGNMMVFLILISCKKLKSMTDIYLFNLAISDLLFLLTLPFWAHYAANE--WVFGNIMCKLFTGIYHIGYFGGIFFIILLTIDRYLAIVHAVAIKARTVNFGVITSVVTWVVAVFVSLPEIIFMRSQK-EG-HYTCSPHFLHIQYRFWKHFQTLKMVILSLILPLLVMVICYSGILNTL--------FRCRNEKKRHRAVRLIFAIMIVYFLFWTPYNIVLLLTTFQEYFNRLDQAMQVTETLGMTHCCLNPVIYAFVGEKFRNYLSVFFRKHIVKRFCKHCSIFV-------SSVY
>5HT_HELVI
TSSEWFDGSNCSWVDAVSWGCTSTNATSTDVTSFVLMAVTSVVLALIILATIVGNVFVIAAIIIERNLQNVANYLVASLAVADLMVACLVMPLGAVYEVSQGWILGPELCDMWTSSDVLCSSASILHLVAIATDRYWAVTDV-YIHIRNEKRIFTMIVLVWGAALVVSLAPQLGWKDPDT-Q--QKCLVS--------QDLAYQIFATMSTFYVPLAVILILYWKIFQTARRRPAPEKKESLEAKRERKAAKTLAIITGAFVFCWLPFFIMALVMPICQTCVISDYLASFFLWLGYFNSTLNPVIYTIFSPDFRQAFARILFGTHRRRRYKKF---------------
>NTR2_HUMAN
-----METSSPRPPRPSSNPGLSLDARLGVDTRLWAKVLFTALYALIWALGAAGNALSVHVVLKARARAGRLRHHVLSLALAGLLLLLVGVPVELYSFVWFHWVFGDLGCRGYYFVHELCAYATVLSVAGLSAERCLAVCQPLARSLLTPRRTRWLVALSWAASLGLALPMAVIMGQKH--AASRVCTVL---VVSRTALQVFIQVNVLVSFVLPLALTAFLNGVTVSHLLALGQVRHKDVRRIRSLQRSVQVLRAIVVMYVICWLPYHARRLMYCYVPDDNFYHYFYMVTNTLFYVSSAVTPLLYNAVSSSFRKLFLEAVSSLCGEHHPMKRLPPMDTASGFGDPPE
>NY4R_MOUSE
HFLAPLFPGSLQGKNGTNPLDSPYNFSDGCQDSAELLAFIITTYSIETILGVLGNLCLIFVTTRQKEKSNVTNLLIANLAFSDFLMCLICQPLTVTYTIMDYWIFGEVLCKMLTFIQCMSVTVSILSLVLVALERHQLIINPT-GWKPSIFQAYLGIVVIWFISCFLSLPFLANSTLNDED--KVVCFVS---WSSDHHRLIYTTFLLLFQYCIPLAFILVCYIRIYQRLQRQKHVAHACSSRAGQMKRINSMLMTMVTAFAVLWLPLHVFNTLEDWYQEACHGNLIFLMCHLLAMASTCVNPFIYGFLNINFKKDIKALVLTCHCRSPQGES---TVHTDLSKGSMR
>PE21_RAT
--MSPYG-LNLSLVDEATTCVTPRVPNTSVVLPTGGNGTSPALPIFSMTLGAVSNVLALALLAQVAGSTATFLLFVASLLAIDLAGHVIPGALVLRLYTAG-RAPAGGACHFLGGCMVFFGLCPLLLGCGMAVERCVGVTQPLHAARVSVARARLALALLAAMALAVALLPLVHVGHYELQYPGTWCFIS--LGPPGGWRQALLAGLFAGLGLAALLAALVCNTLSGLALLRDRRRSRGLRRVHAHDVEMVGQLVGIMVVSCICWSPLLVLVVLAIGGWNSNSLQRPLFLAVRLASWNQILDPWVYILLRQAMLRQLLRLLPLRVSAKGGPTELSLSSLRSSRHSGFS
>OPSG_RAT
DHYEDSTQASIFTYTNSNSTRGPFEGPNYHIAPRWVYHLTSTWMILVVIASVFTNGLVLAATMRFKKLRHPLNWILVNLAVADLAETIIASTISVVNQIYGYFVLGHPLCVIEGYIVSLCGITGLWSLAIISWERWLVVCKPFGNVRFDAKLATVGIVFSWVWAAVWTAPPIFGWSRYWPYGLKTSCGPDVFSGTSYPGVQSYMMVLMVTCCIFPLSIIVLCYLQVWLAIRAVAKQQKESESTQKAEKEVTRMVVVMVFAYCLCWGPYTFFACFATAHPGYAFHPLVASLPSYFAKSATIYNPIIYVFMNRQFRNCILQLFGKKVDDS-----SELVSSV--SSVSPA
>P2UR_MOUSE
----MAADLEPWNSTINGTWEGDELGYKCRFNEDFKYVLLPVSYGVVCVLGLCLNVVALYIFLCRLKTWNASTTYMFHLAVSDSLYAASLPLLVYYYARGDHWPFSTVLCKLVRFLFYTNLYCSILFLTCISVHRCLGVLRPLSLRWGRARYARRVAAVVWVLVLACQAPVLYFVTTS-----RITCHDT-SARELFSHFVAYSSVMLGLLFAVPFSVILVCYVLMARRLLKP---YGTTGGLPRAKRKSVRTIALVLAVFALCFLPFHVTRTLYYSFRSLNAINMAYKITRPLASANSCLDPVLYFLAGQRLVRFARDAKPPTEPTPSPQARRKL--TVRKDLSVSS
>A2AC_DIDMA
-MDLQLTTNSTDSGDRGGSSNESLQRQPPSQYSPAEVAGLAAVVSFLIVFTIVGNVLVVIPVLTSRALKAPQNLFLVSLASADILVATLVMPFSLANELMNYWYFGKVWCDIYLALDVLFCTSSIVHLCAISLDRYWSVTQAVYNLKRTPRRIKGIIVTVWLISAVISFPPLISLYR-------PQCELN--------DETWYILSSCIGSFFAPCIIMVLVYVRIYRVAKLRRRKLCRRKVTQAREKRFTFVLAVVMGVFVVCWFPFFFTYSLYGICREAQVPETLFKFFFWFGYCNSSLNPVIYTIFNQDFRRSFKHILFKKKKKTSLQ-----------------
>GPR1_RAT
EVSREMLFEELDNYSYALEYYSQEPDAEENVYPGIVHWISLLLYALAFVLGIPGNAIVIWFMGFKWK-KTVTTLWFLNLAIADFIFVLFLPLYISYVALSFHWPFGRWLCKLNSFIAQLNMFSSVFFLTVISLDRYIHLIHPGSHPHRTLKNSLLVVLFVWLLASLLGGPTLYFR--D--NN-RIICYNNQEYELTLMRHHVLTWVKFLFGYLLPLLTMSSCYLCLIFKT---------KKQNILISSKHLWMILSVVIAFMVCWTPFHLFSIWELSIHHNNVLQGGIPLSTGLAFLNSCLNPILYVIISKKFQARFRASVAEVLKRSLWEASCSGSAETKSLSLLET
>CCR4_PAPAN
IYTSDNYTEEMG-SGDYDSIKEPC---FREENAHFNRIFLPTIYSIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVADLLFVITLPFWAVDAVAN--WYFGNFLCKAVHVIYTVNLYSSVLILAFISLDRYLAIVHATSQRPRKLLAEKVVYVGVWIPALLLTIPDFIFASVSE-DD-RYICDRF---YPNDLWVVVFQFQHIMVGLILPGIVILSCYCIIISKL---------SHSKGHQKRKALKTTVILILAFFACWLPYYIGISIDSFILLENTVHKWISITEALAFFHCCLNPILYAFLGAKFKTSAQHALTSVSRGSSLKILSKG----------GH
>D2D1_.ENLA
---MDPQNLSMYNDD------INNGTNGTAVDQKPHYNYYAMLLTLLVFVIVFGNVLVCIAVSREKALQTTTNYLIVSLAVADLLVATLVMPWAVYMEVVGEWRFSRIHCDIFVTLDVMMCTASILNLCAISIDRYTAVAMPMNTRYSSKRRVTVMISVVWVLSFAISCPLLFGLN---S----KVCII---------DNPAFVIYSSIVSFYVPFIVTLLVYVQIYIVLRKRRTSMSKKKLSQHKEKKATQMLAIVLGVFIICWLPFFIIHILNMHCNCN-IPQALYSAFTWLGYVNSAVNPIIYTTFNVEFRKAFIKILHC-------------------------
>OLF8_RAT
---------------MNNKTVITHFLLLGLPIPPEHQQLFFALFLIMYLTTFLGNLLIVVLVQLDSHLHTPMYLFLSNLSFSDLCFSSVTMLKLLQNIQSQVPSISYAGCLTQIFFFLLFGYLGNFLLVAMAYDRYVAICFPLYTNIMSHKLCTCLLLVFWIMTSSHAMMHTLLAARL---SLNFFCDLFKLACSDTYVNELMIHIMGVIIIVIPFVLIVISYAKIISSI--------LKVPSTQSIHKVFSTCGSHLSVVSLFYGTIIGLYLCPSGDN---FSLKGSAMAMMYTVVTPMLNPFIYSLRNRDMKQALIRVTCSKKISLPW------------------
>GASR_MOUSE
GSSLCRPGVSLLNSSSAGNLSCETPRIRGTGTRELELTIRITLYAVIFLMSVGGNVLIIVVLGLSRRLRTVTNAFLLSLAVSDLLLAVACMPFTLLPNLMGTFIFGTVICKAVSYLMGVSVSVSTLNLAAIALERYSAICRPLARVWQTRSHAARVILATWLLSGLLMVPYPVYTVVQP---I-LQCMHL---WPSERVQQMWSVLLLILLFFIPGVVMAVAYGLISRELYLGTGPPRPNQAKLLAKKRVVRMLLVIVLLFFVCWLPVYSANTWRAFDGPGALAGAPISFIHLLSYTSACANPLVYCFMHRRFRQACLDTCARCCPRPPRARPRPLPSIASLSRLSYT
>OPRK_CAVPO
LPNGSAWLPGWAEPDGNGSAGPQDEQLEPAHISPAIPVIITAVYSVVFVVGLVGNSLVMFVIIRYTKMKTATNIYIFNLALADALVTTTMPFQSTVYLMNS-WPFGDVLCKIVISIDYYNMFTSIFTLTMMSVDRYIAVCHPVALDFRTPLKAKIINICIWLLSSSVGISAIILGGTKVDVD-IIECSLQFPDDDYSWWDLFMKICVFVFAFVIPVLIIIVCYTLMILRLKSVRLL-SGSREKDRNLRRITRLVLVVVAVFIICWTPIHIFILVEALGSTSTAALSSYYFCIALGYTNSSLNPILYAFLDENFKRCFRDFCFPIKMRMERQSTSRVAYMRNVDGVNKP
>OLFD_CANFA
-------------MTEKNQTVVSEFVLLGLPIDPDQRDLFYALFLAMYVTTILGNLLIIVLIQLDSHLHTPMYLFLSNLSFSDLCFSSVTMPKLLQNMQSQVPSIPYAGCLTQMYFFLFFGDLESFLLVAMAYDRYVAICFPLYTTIMSPKLCFSLLVLSWVLTMFHAVLHTLLMARL---CPHFFCDMSKLACSDTQVNELVIFIMGGLILVIPFLLIITSYARIVSSI--------LKVPSAIGICKVFSTCGSHLSVVSLFYGTVIGLYLCPSANN---STVKETIMAMMYTVVTPMLNPFIYSLRNKDMKGALRRVICRKKITFSV------------------
>O2F1_HUMAN
-------------MGTDNQTWVSEFILLGLSSDWDTRVSLFVLFLVMYVVTVLGNCLIVLLIRLDSRLHTPMYFFLTNLSLVDVSYATSVVPQLLAHFLAEHKAIPFQSCAAQLFFSLALGGIEFVLLAVMAYDRYVAVCDALYSAIMHGGLCARLAITSWVSGFISSPVQTAITFQL---PDHISCELLRLACVDTSSNEVTIMVSSIVLLMTPFCLVLLSYIQIISTI--------LKIQSREGRKKAFHTCASHLTVVALCYGVAIFTYIQPHSSP---SVLQEKLFSVFYAILTPMLNPMIYSLRNKEVKGAWQKLLWKFSGLTSKLAT---------------
>OPS._MOUSE
------------MLSEASDFNSSGSRSEGSVFSRTEHSVIAAYLIVAGITSILSNVVVLGIFIKYKELRTPTNAVIINLAFTDIGVSSIGYPMSAASDLHGSWKFGHAGCQIYAGLNIFFGMVSIGLLTVVAMDRYLTISCPDVGRRMTTNTYLSMILGAWINGLFWALMPIIGWASYA--PTGATCTIN--WRNNDTSFVSYTMMVIVVNFIVPLTVMFYCYYHVSRSLRLYAS-TAHLHRDWADQADVTKMSVIMILMFLLAWSPYSIVCLWACFGNPKKIPPSMAIIAPLFAKSSTFYNPCIYVAAHKKFRKAMLAMFKCQPHLAVPEPSTLPLAPVRI------
>A1AB_RAT
HNTSAPAHWGELKDDNFTGPNQTSSNSTLPQLDVTRAISVGLVLGAFILFAIVGNILVILSVACNRHLRTPTNYFIVNLAIADLLLSFTVLPFSATLEVLGYWVLGRIFCDIWAAVDVLCCTASILSLCAISIDRYIGVRYSLYPTLVTRRKAILALLSVWVLSTVISIGPLLGWKEP--AP--KECGVT--------EEPFYALFSSLGSFYIPLAVILVMYCRVYIVAKRTHNPIAVKLFKFSREKKAAKTLGIVVGMFILCWLPFFIALPLGSLFSTLKPPDAVFKVVFWLGYFNSCLNPIIYPCSSKEFKRAFMRILGCQCRGGRRRRRRRRTYRPWTRGGSLE
>D4DR_RAT
---MGNSSATGDGGLLAGRGP---ESLGTGTGLGGAGAAALVGGVLLIGMVLAGNSLVCVSVASERILQTPTNYFIVSLAAADLLLAVLVLPLFVYSEVQGGWLLSPRLCDTLMAMDVMLCTASIFNLCAISVDRFVAVTVPL-RYNQQGQCQLLLIAATWLLSAAVAAPVVCGLN---G-R--TVCCL---------EDRDYVVYSSICSFFLPCPLMLLLYWATFRGLRRWPAPRKRGAKITGRERKAMRVLPVVVGAFLMCWTPFFVVHITRALCPACFVSPRLVSAVTWLGYVNSALNPIIYTIFNAEFRSVFRKTLRLRC-----------------------
>GP52_HUMAN
SRWTEWRILNMSSGIVNASERHSCPLGFGHYSVVDVCIFETVVIVLLTFLIIAGNLTVIFAFHCAPLHHYTTSYFIQTMAYADLFVGVSCLVPTLSLLHYSTGVHESLTCRVFGYIISVLKSVSMACLACISVDRYLAITKPLYNQLVTPCRLRICIILIWIYSCLIFLPSFFGWG---PGY-IFEWCAT-----SWLTSAYFTGFIVCLLYAPAAFVVCFTYFHIFKICRQHFPSDSSRETGHSPDRRYAMVLFRITSVFYMLWLPYIIYFLLESSRV--LDNPTLSFLTTWLAVSNSFCNCVIYSLSNGVFRLGLRRLFETMCTSCMCVKDQEARANSCSI-----
>ACTR_MOUSE
-------------MKHIINSYEHTNDTARNNSDCPDVVLPEEIFFTISVIGILENLIVLLAVIKNKNLQSPMYFFICSLAISDMLGSLYKILENILIMFRNMGSFESTADDIIDCMFILSLLGSIFSLSVIAADRYITIFHALYHSIVTMRRTIITLTIIWMFCTGSGITMVIFFS----------------------HHHIPTVLTFTSLFPLMLVFILCLYIHMFLLARSH---ARKISTLPRTNMKGAMTLTILLGVFIFCWAPFVLHVLLMTFCPNNVCYMSLFQVNGMLIMCNAVIDPFIYAFRSPELRDAFKRMLFCNRY----------------------
>SSR4_HUMAN
GEEGLGTAWPSAANASSAPAEAEEAVAGPGDARAAGMVAIQCIYALVCLVGLVGNALVIFVILRYAKMKTATNIYLLNLAVADELFMLSVPFVASSAALRH-WPFGSVLCRAVLSVDGLNMFTSVFCLTVLSVDRYVAVVHPLAATYRRPSVAKLINLGVWLASLLVTLPIAIFADTRPGGQ-AVACNLQ---WPHPAWSAVFVVYTFLLGFLLPVLAIGLCYLLIVGKMRAVALR-AGWQQRRRSEKKITRLVLMVVVVFVLCWMPFYVVQLLNLVVTS--LDATVNHVSLILSYANSCANPILYGFLSDNFRRSFQRVLCLRCCLLEGAGGAEETALKSKGGAGCM
>OPSD_PIG
MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFMLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLALTWVMALACAAPPLVGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFSIPLVIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVVAFLICWLPYASVAFYIFTHQGSDFGPIFMTIPAFFAKSASIYNPVIYIMMNKQFRNCMLTTLCCGKNPLGDDEASTTTSKTETSQVAPA
>GRHR_SHEEP
--MANGDSPDQNENHCSAINSSILLTPGSLPTLTLSGKIRVTVTFFLFLLSTIFNTSFLLKLQNWTQKLSKMKVLLKHLTLANLLETLIVMPLDGMWNITVQWYAGELLCKVLSYLKLFSMYAPAFMMVVISLDRSLAITRPL-AVKSNSKLGQFMIGLAWLLSSIFAGPQLYIFGMIHEG--FSQCVTH--SFPQWWHQAFYNFFTFSCLFIIPLLIMLICNAKIIFTLTRVPHKNQSKNNIPQARLRTLKMTVAFATSFTVCWTPYYVLGIWYWFDPDMRVSDPVNHFFFLFAFLNPCFDPLIYGYFSL-------------------------------------
>DOP1_DROME
YDGTTLTSFYNESSWTNASEMDTIVGEEPEPLSLVSIVVVGIFLSVLIFLSVAGNILVCLAIYTDGSPRIG-NLFLASLAIADLFVASLVMTFAGVNDLLGYWIFGAQFCDTWVAFDVMCSTASILNLCAISMDRYIHIKDPLYGRWVTRRVAVITIAAIWLLAAFVSFVPISLGIHRPLI-KYPTCALD--------LTPTYAVVSSCISFYFPCVVMIGIYCRLYCYAQKHFKVHTHSSPYHVSDHKAAVTVGVIMGVFLICWVPFFCVNITAAFCKTC-IGGQTFKILTWLGYSNSAFNPIIYSIFNKEFRDAFKRILTMRNP----------------------
>AA1R_RABIT
----------------------------MPPSISAFQAAYIGIEVLIALVSVPGNVLVIWAVKVNQALRDATFCFIVSLAVADVAVGALVIPLAILINIG--PETYFHTCLMVACPVLILTQSSILALLAIAVDRYLRVKIPLYKAVVTPRRAAVAIAGCWILSLVVGLTPMFGWNNLRN-G-VIKCEFE-----KVISMEYMVYFNFFVWVLPPLLLMVLIYLEVFYLIRRQK-ASGDPHKYYGKELKIAKSLALILFLFALSWLPLHILNCVTLFCPSCQKPSILVYTAIFLTHGNSAMNPIVYAFRIHKFRVTFLKIWNDHFRCRPAPAGDGDPND---------
>OPSI_ASTFA
LNGFEGDNFYIPMSNRTGLVRDPFVYEQYYLAEPWQFKLLACYMFFLICLGLPINGFTLFVTAQHKKLQQPLNFILVNLAVAGMIMVCFGFTITISSAVNGYFYFGPTACAIEGFMATLGGEVALWSLVVLAIERYIVVCKPMGSFKFSASHALGGIGFTWFMAMTCAAPPLVGWSRYIPEGLQCSCGPDYYTLNPKYNNESYVIYMFVVHFIVPVTVIFFTYGRLVCTVKSAAAAQQDSASTQKAEKEVTRMVILMVVGFLVAWTPYATVAAWIFFNKGAAFTAQFMAVPAFFSKSSALFNPIIYVLLNKQFRNCMLTTLFCGKNPLGDEESSTVTVSS----VSPA
>ACTR_HUMAN
-------------MKHIINSYENINNTARNNSDCPRVVLPEEIFFTISIVGVLENLIVLLAVFKNKNLQAPMYFFICSLAISDMLGSLYKILENILIILRNMGSFETTADDIIDSLFVLSLLGSIFSLSVIAADRYITIFHALYHSIVTMRRTVVVLTVIWTFCTGTGITMVIFFS----------------------HHHVPTVITFTSLFPLMLVFILCLYVHMFLLARSH---TRKISTLPRANMKGAITLTILLGVFIFCWAPFVLHVLLMTFCPSNACYMSLFQVNGMLIMCNAVIDPFIYAFRSPELRDAFKKMIFCSRYW---------------------
>CCR4_RAT
IYTSDNYSEEVG-SGDYDSNKEPC---FRDENENFNRIFLPTIYFIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVADLLFVITLPFWAVDAMAD--WYFGKFLCKAVHIIYTVNLYSSVLILAFISLDRYLAIVHATSQSARKLLAEKAVYVGVWIPALLLTIPDIIFADVSQ-DG-RYICDRL---YPDSLWMVVFQFQHIMVGLILPGIVILSCYCIIISKL---------SHSKGHQKRKALKTTVILILAFFACWLPYYVGISIDSFILLESVVHKWISITEALAFFHCCLNPILYAFLGAKFKSSAQHALNSMSRGSSLKILSKG----------GH
>OPSD_ANOCA
MNGTEGQNFYVPMSNKTGVVRNPFEYPQYYLADPWQFSALAAYMFLLILLGFPINFLTLFVTIQHKKLRTPLNYILLNLAVANLFMVLMGFTTTMYTSMNGYFIFGTVGCNIEGFFATLGGEMGLWSLVVLAVERYVVICKPMSNFRFGETHALIGVSCTWIMALACAGPPLLGWSRYIPEGMQCSCGVDYYTPTPEVHNESFVIYMFLVHFVTPLTIIFFCYGRLVCTVKAAAAQQQESATTQKAEREVTRMVVIMVISFLVCWVPYASVAFYIFTHQGSDFGPVFMTIPAFFAKSSAIYNPVIYILMNKQFRNCMIMTLCCGKNPLGDEETSAGTSTVSTSQVSPA
>ML1A_CHICK
-------MRANGSELNGTVLPRDPPAEGSPRRPPWVTSTLATILIFTIVVDLLGNLLVILSVYRNKKLRNAGNIFVVSLAIADLVVAIYPYPLVLTSVFHNGWNLGYLHCQISGFLMGLSVIGSIFNITGIAINRYCYICHSLYDKLYSDKNSLCYVGLIWVLTVVAIVPNLFVGS-LQYDP-IYSCTFA------QSVSSAYTIAVVFFHFILPIAIVTYCYLRIWILVIQVRR-PDNNPRLKPHDFRNFVTMFVVFVLFAVCWAPLNFIGLAVAVDPETRIPEWLFVSSYYMAYFNSCLNAIIYGLLNQNFRREYKKIVVSFCTAKAFFQDSSNSKPSPLITNNNQ
>EDG1_HUMAN
VKAHRSSVSDYVNYDIIVRHYNYTGKLNISADKENSIKLTSVVFILICCFIILENIFVLLTIWKTKKFHRPMYYFIGNLALSDLLAG-VAYTANLLLSGATTYKLTPAQWFLREGSMFVALSASVFSLLAIAIERYITMLKMKLHNGSNNFRLFLLISACWVISLILGGLPIMGWNCI-S--ALSSCSTV------LPLYHKHYILFCTTVFTLLLLSIVILYCRIYSLVRTRLTFISKASRSSE-NVALLKTVIIVLSVFIACWAPLFILLLLDVGCKVKCDILFRAEYFLVLAVLNSGTNPIIYTLTNKEMRRAFIRIMSCCKCPSGDSAGKFKEFSRSKSDNSSH
>OPS1_DROPS
APS---NGSVVDKVTPDMAHLISPYWDQFPAMDPIWAKILTAYMIIIGMISWCGNGVVIYIFATTKSLRTPANLLVINLAISDFGIMITNTPMMGINLYFETWVLGPMMCDIYAGLGSAFGCSSIWSMCMISLDRYQVIVKGMAGRPMTIPLALGKIAYIWFMSTIWCCLAPVFWSRYVPEGNLTSCGID--YLERDWNPRSYLIFYSIFVYYIPLFLICYSYWFIIAAVSAHMNVRSSEDADKSAEGKLAKVALVTISLWFMAWTPYLVINCMGLFKFEG-LTPLNTIWGACFAKSAACYNPIVYGISHPKYRLALKEKCPCCVFGKVDDGK-SSTSEAESKA----
>PAFR_CAVPO
----------------------MELNSSSRVDSEFRYTLFPIVYSIIFVLGIIANGYVLWVFARLYPKLNEIKIFMVNLTVADLLFLITLPLWIVYYSNQGNWFLPKFLCNLAGCLFFINTYCSVAFLGVITYNRFQAVKYPITAQATTRKRGIALSLVIWVAIVAAASYFLVMDSTN-G--NITRCFEH---YEKGSKPVLIIHICIVLGFFIVFLLILFCNLVIIHTLLRQP---VKQQRNAEVRRRALWMVCTVLAVFVICFVPHHMVQLPWTLAELGQAINDAHQVTLCLLSTNCVLDPVIYCFLTKKFRKHLSEKLNIMRSSQKCSRVTTDPINHTPVNPIKN
>OPSB_MOUSE
-----MSGEDDFYLFQNISSVGPWDGPQYHLAPVWAFRLQAAFMGFVFFVGTPLNAIVLVATLHYKKLRQPLNYILVNVSLGGFLFCIFSVFTVFIASCHGYFLFGRHVCALEAFLGSVAGLVTGWSLAFLAFERYVVICKPFGSIRFNSKHALMVVLATWIIGIGVSIPPFFGWSRFIPEGLQCSCGPDWYTVGTKYRSEYYTWFLFIFCFIIPLSLICFSYSQLLRTLRAVAAQQQESATTQKAEREVSHMVVVMVGSFCLCYVPYAALAMYMVNNRNHGLDLRLVTIPAFFSKSSCVYNPIIYCFMNKQFRACILEMVCRKPMAD---ESDVSSTVSSSKVGPH-
>HH2R_CANFA
-------------------MISNGTGSSFCLDSPPCRITVSVVLTVLILITIAGNVVVCLAVGLNRRLRSLTNCFIVSLAITDLLLGLLVLPFSAFYQLSCRWSFGKVFCNIYTSLDVMLCTASILNLFMISLDRYCAVTDPLYPVLITPVRVAVSLVLIWVISITLSFLSIHLGWN--RN-TIPKCKVQ--------VNLVYGLVDGLVTFYLPLLVMCITYYRIFKIARDQR--MGSWKAATIGEHKATVTLAAVMGAFIICWFPYFTVFVYRGLKGDD-INEAFEAVVLWLGYANSALNPILYATLNRDFRTAYQQLFRCRPASHNAQETSLRRNQSREPMR---
>C3.1_MOUSE
--MSTSFPELDLENFEYDDSAEAC---YLGDIVAFGTIFLSVFYALVFTFGLVGNLLVVLALTNSRKPKSITDIYLLNLALSDLLFVATLPFWTHYLISHEG--LHNAMCKLTTAFFFIGFFGGIFFITVISIDRYLAIVLAASMNNRTVQHGVTISLGVWAAAILVASPQFMFTKRKD-----NECLGDYPEVLQEMWPVLRNSEVNILGFALPLLIMSFCYFRIIQTL---------FSCKNRKKARAVRLILLVVFAFFLFWTPYNIMIFLETLKFYNRDLRLALSVTETVAFSHCCLNPFIYAFAGEKFRRYLGHLYRKCLAVLCGHPVHTGRSRQDSILSS-F
>ML1._SHEEP
---------MGRTLAVPTPYGCIGCKLPQPDYPPALIVFMFCAMVITIVVDLIGNSMVILAVSKNKKLRNSGNVFVVSLSVADMLVAIYPYPLMLHAMAIGGWDLSKLQCQMVGFITGLSVVGSIFNIMAIAINRYCYICHSLYERIFSVRNTCIYLAVTWIMTVLAVLPNMYIGT-IEYDP-TYTCIFN------YVNNPAFAVTIVCIHFVLPLLIVGFCYVKIWTKVLAARD-AGQNPDNQLAEVRNFLTMFVIFLLFAVCWCPINALTVLVAVNPKEKIPNWVYLAAYFIAYFNSCLNAVIYGVLNENFRREYWTIFHAMRHPVLFLSGLLTAQAHTHARARAR
>GRPR_MOUSE
NCSHLNLDVDPFLSCNDTFNQSLSPPKMDNWFHPGFIYVIPAVYGLIIVIGLIGNITLIKIFCTVKSMRNVPNLFISSLALGDLLLLVTCAPVDASKYLADRWLFGRIGCKLIPFIQLTSVGVSVFTLTALSADRYKAIVRPMIQASHALMKICLKAALIWIVSMLLAIPEAVFSDLHPQT--FISCAPY---HSNELHPKIHSMASFLVFYVIPLAIISVYYYFIARNLIQSLPVNIHVKKQIESRKRLAKTVLVFVGLFAFCWLPNHVIYLYRSYHYSEMLHFVTSICARLLAFTNSCVNPFALYLLSKSFRKQFNTQLLCCQPGLMNR--SHSMTSFKSTNP-SA
>OPS._HUMAN
------------MLRNNLGNSSDSKNEDGSVFSQTEHNIVATYLIMAGMISIISNIIVLGIFIKYKELRTPTNAIIINLAVTDIGVSSIGYPMSAASDLYGSWKFGYAGCQVYAGLNIFFGMASIGLLTVVAVDRYLTICLPDVGRRMTTNTYIGLILGAWINGLFWALMPIIGWASYA--PTGATCTIN--WRKNDRSFVSYTMTVIAINFIVPLTVMFYCYYHVTLSIKHHSD-TESLNRDWSDQIDVTKMSVIMICMFLVAWSPYSIVCLWASFGDPKKIPPPMAIIAPLFAKSSTFYNPCIYVVANKKFRRAMLAMFKCQTHQTMPVTSILPLASGRI------
>5HT1_DROME
VSYQGITSSNLGDSNTTLVPLLEEFAAGEFVLPPLTSIFVSIVLLIVILGTVVGNVLVCIAVCMVRKLRRPCNYLLVSLALSDLCVALLVMPMALLYEVLEKWNFGPLLCDIWVSFDVLCCTASILNLCAISVDRYLAITKPLYGVKRTPRRMMLCVGIVWLAAACISLPPLLILGN--E-G--PICTVC--------QNFAYQIYATLGSFYIPLSVMLFVYYQIFRAARRILLGHKKLRFQLAKEKKASTTLGIIMSAFTVCWLPFFILALIRPFETM-HVPASLSSLFLWLGYANSLLNPIIYATLNRDFRKPFQEILYFRCSSLNTMMRENYPPSQRVMLGDER
>CRH2_MOUSE
------MANVTLKPLCPLLEEMVQLPNHSNSSLRYIDHVSVLLHGLASLLGLVENGLILFVVGCRMR-QTVVTTWVLHLALSDLLAAASLPFFTYFLAVGHSWELGTTFCKLHSSVFFLNMFASGFLLSAISLDRCLQVVRPVAQNHRTVAVAHRVCLMLWALAVLNTIPYFVFRDTIPWNPRDTTCDY---------RQKALAVSKFLLAFMVPLAIIASSHVAVSLRL---------HHRGRQRTGRFVRLVAAIVVAFVLCWGPYHIFSLLEARAHS-QLASRGLPFVTSLAFFNSVVNPLLYVFTCPDMLYKLRRSLRAVLESVLVEDSDQSRRASSTATPAST
>ML1C_CHICK
-------------MERPGSNGSCSGCRLEGGPAARAASGLAAVLIVTIVVDVLGNALVILSVLRNKKLRNAGNIFVVSLSVADLVVAVYPYPLILSAIFHNGWTMGNIHCQISGFLMGLSVIGSIFNITAIAINRYCYICHSLYDKLFNLKNTCCYICLTWTLTVVAIVPNFFVGS-LQYDP-IYSCTFA------QTVSTSYTITVVVVHFIVPLSIVTFCYLRIWILVIQVHR-QDCKQKIRAADIRNFLTMFVVFVLFAVCWGPLNFIGLAVSINPSKHIPEWLFVLSYFMAYFNSCLNAVIYGLLNQNFRKEYKRILLMLRTPRLLFIDVSKSKPSPAVTNNNQ
>5H1F_HUMAN
-------------MDFLNSSD-QNLTSEELLNRMPSKILVSLTLSGLALMTTTINSLVIAAIIVTRKLHHPANYLICSLAVTDFLVAVLVMPFSIVYIVRESWIMGQVVCDIWLSVDITCCTCSILHLSAIALDRYRAITDAVYARKRTPKHAGIMITIVWIISVFISMPPLFWR----S-R--DECIIK-------HDHIVSTIYSTFGAFYIPLALILILYYKIYRAAKTLFKHWRRQKISGTRERKAATTLGLILGAFVICWLPFFVKELVVNVCDKCKISEEMSNFLAWLGYLNSLINPLIYTIFNEDFKKAFQKLVRCRC-----------------------
>GRHR_PIG
--MANSASPEQNQNHCSAINSSILLTQGNLPTLTLSPNIRVTVTFFLFLLSTAFNASFLLKLQKWTQKLSRMKVLLKHLTLANLLETLIVMPLDGMWNITVQWYAGEFLCKVLSYLKLFSMYAPAFMMVVISLDRSLAITRPL-AVKSNSRLGRFMIGLAWLLSSIFAGPQLYIFRMIHEG--FSQCVTH--SFPQWWHQAFYDFFTFSCLFIIPLLIMLICNAKIMFTLTRVPHNNQSKNNIPRARLRTLKMTVAFAASFIVCWTPYLVLGIWYWFDPEMRVSDPVNHFFFLFAFLNPCFDPLIYGYFSL-------------------------------------
>MAS_MOUSE
------MDQSNMTSLAEEKAMNTSSRNASLGSSHPPIPIVHWVIMSISPLGFVENGILLWFLCFRMR-RNPFTVYITHLSMADISLLFCIFILSIDYALDYESSGHHYTIVTLSVTFLFGYNTGLYLLTAISVERCLSVLYPIYTSHRPKHQSAFVCALLCALSCLVTTMEYVMCI-DSHS--RSDC-----------RAVIIFIAILSFLVFTPLMLVSSSILVVKIRK----------NTWASHSSKLYIVIMVTIIIFLIFAMPMRVLYLLYYEYW--SAFGNLHNISLLFSTINSSANPFIYFFVGSSKKKRFRESLKVVLTRAFKDEMQPRTVSIETVV----
>UL33_HCMVA
----------------------------MTGPLFAIRTTEAVLNTFIIFVGGPLNAIVLITQLLTNRGYSTPTIYMTNLYSTNFLTLTVLPFIVLSNQ-WL-LPAGVASCKFLSVIYYSSCTVGFATVALIAADRYR-VLHKRTYARQSYRSTYMILLLTWLAGLIFSVPAAVYTTVVMN-NGHATCVLYVAEEVH-TVLLSWKVLLTMVWGAAPVIMMTWFYAFFYSTV---------QRTSQKQRSRTLTFVSVLLISFVALQTPYVSLMIFNSYATTATLRRTIGTLARVVPHLHCLINPILYALLGHDFLQRMRQCFRGQLLDRRAFLRSQQTNLAAGNNSQSV
>CKR5_GORGO
----MDYQVSSPTYDIDYYTSEPC---QKTNVKQIAARLLPPLYSLVFIFGFVGNMLVILILINCKRLKSMTDIYLLNLAISDLFFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWVVAVFASLPGIIFTRSQK-EG-HYTCSSHFPYSQYQFWKNFQTLKIVILGLVLPLLVMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCSIFA-------SSVY
>PE21_MOUSE
--MSPCG-LNLSLADEAATCATPRLPNTSVVLPTGDNGTSPALPIFSMTLGAVSNVLALALLAQVAGSAATFLLFVASLLAIDLAGHVIPGALVLRLYTAG-RAPAGGACHFLGGCMVFFGLCPLLLGCGMAVERCVGVTQPLHAARVSVARARLALAVLAAMALAVALLPLVHVGRYELQYPGTWCFIS--LGPRGGWRQALLAGLFAGLGLAALLAALVCNTLSGLALLRDRRRSRGPRRVHAHDVEMVGQLVGIMVVSCICWSPLLVLVVLAIGGWNSNSLQRPLFLAVRLASWNQILDPWVYILLRQAMLRQLLRLLPLRVSAKGGPTELGLSSLRSSRHSGFS
>A1AB_MESAU
HNTSAPAQWGELKDANFTGPNQTSSNSTLPQLDVTRAISVGLVLGAFILFAIVGNILVILSVACNRHLRTPTNYFIVNLAIADLLLSFTVLPFSATLEVLGYWVLGRIFCDIWAAVDVLCCTASILSLCAISIDRYIGVRYSLYPTLVTRRKAILALLSVWVLSTVISIGPLLGWKEP--AP--KECGVT--------EEPFYALFSSLGSFYIPLAVILVMYCRVYIVAKRTHNPIAVKLFKFSREKKAAKTLGIVVGMFILCWLPFFIALPLGSLFSTLKPPDAVFKVVFWLGYFNSCLNPIIYPCSSKEFKRAFMRILGCQCRSGRRRRRRRRTYRPWTRGGSLE
>O3A2_HUMAN
----------MEPEAGTNRTAVAEFILLGLVQTEEMQPVVFVLLLFAYLVTTGGNLSILAAVLVEPKLHAPMYFFLGNLSVLDVGCITVTVPAMLGRLLSHKSTISYDACLSQLFFFHLLAGMDCFLLTAMAYDRLLAICQPLYSTRMSQTVQRMLVAASLACAFTNALTHTVAMSTL---NNHFYCDLPQLSCSSTQLNELLLFAVGFIMAGTPLVLIITAYSHVAAAV--------LRIRSVEGRKKAFSTCGSHLTVVCLFFGRGIFNYMRLGSEE---ASDKDKGVGVFNTVINPMLNPLIYSLRNPDVQGALWQIFLGRRSLTA-------------------
>OPSU_CARAU
-------MDAWTYQFGNLSKISPFEGPQYHLAPKWAFYLQAAFMGFVFFVGTPLNAIVLFVTMKYKKLRQPLNYILVNISLGGFIFDTFSVSQVFFSALRGYYFFGYTLCAMEAAMGSIAGLVTGWSLAVLAFERYVVICKPFGSFKFGQSQALGAVALTWIIGIGCATPPFWGWSRYIPEGIGTACGPDWYTKNEEYNTESYTYFLLVSCFMMPIMIITFSYSQLLGALRAVAAQQAESASTQKAEKEVSRMVVVMVGSFVVCYGPYAITALYFSYAEDSNKDYRLVAIPSLFSKSSCVYNPLIYAFMNKQFNACIMETVFGKKIDE-----SSESSVSA-------
>PAFR_RAT
----------------------MEQNGSFRVDSEFRYTLFPIVYSVIFVLGVVANGYVLWVFATLYPKLNEIKIFMVNLTVADLLFLMTLPLWIVYYSNEGDWIVHKFLCNLAGCLFFINTYCSVAFLGVITYNRYQAVAYPITAQATTRKRGITLSLVIWISIAATASYFLATDSTN-G--NITRCFEH---YEPYSVPILVVHIFITSCFFLVFFLIFYCNMVIIHTLLTRP---VRQQRKPEVKRRALWMVCTVLAVFVICFVPHHVVQLPWTLAELGQAINDAHQITLCLLSTNCVLDPVIYCFLTKKFRKHLSEKFYSMRSSRKCSRATSDPANQTPVLPLKN
>YYO1_CAEEL
ISPNASNYLTYPFDGLCLQKFFYQLQTSLRRFTPYEEIIYTTVYIIISVAAVIGNGLVIMAVVRKKTMRTNRNVLILNLALSNLILAITNIPFLWLPSIDFEFPYSRFFCKFANVLPGSNIYCSTLTISVMAIDRYYSVKKLKASNRKQCFHAVLVSLAIWIVSFILSLPLLLYYETS--MQEVRQCRLVRLPDITQSIQLLMSILQVAFLYIVPLFVLSIFNVKLTRFLKTNDSHLKNNNKTNQRTNRTTSLLIAMAGSYAALWFPFTLITFLIDFELIINLVERIDQTCKMVSMLSICVNPFLYGFLNTNFRHEFSDIYYRYIRCETKSQPAGRIAHHRQDSVYND
>OPS1_LIMPO
PN-----ASVVDTMPKEMLYMIHEHWYAFPPMNPLWYSILGVAMIILGIICVLGNGMVIYLMMTTKSLRTPTNLLVVNLAFSDFCMMAFMMPTMTSNCFAETWILGPFMCEVYGMAGSLFGCASIWSMVMITLDRYNVIVRGMAAAPLTHKKATLLLLFVWIWSGGWTILPFFGWSRYVPEGNLTSCTVD--YLTKDWSSASYVVIYGLAVYFLPLITMIYCYFFIVHAVAEHNVAANADQQKQSAECRLAKVAMMTVGLWFMAWTPYLIISWAGVFSSGTRLTPLATIWGSVFAKANSCYNPIVYGISHPRYKAALYQRFPSLACGSGESGSDVKTMEEKPKIPEA-
>OPSV_APIME
YLPAGPPRLLGWNVPAEELIHIPEHWLVYPEPNPSLHYLLALLYILFTFLALLGNGLVIWIFCAAKSLRTPSNMFVVNLAICDFFMMIKTPIFIYNSFNTG-FALGNLGCQIFAVIGSLTGIGAAITNAAIAYDRYSTIARPL-DGKLSRGQVILFIVLIWTYTIPWALMPVMGVWGRFPEGFLTSCSFD--YLTDTNEIRIFVATIFTFSYCIPMILIIYYYSQIVSHVVNHNVDSNANTSSQSAEIRIAKAAITICFLYVLSWTPYGVMSMIGAFGNKALLTPGVTMIPACTCKAVACLDPYVYAISHPKYRLELQKRLPWLELQEKPISDSTSTPPASS------
>5H1B_DIDMA
PPASGSLTSSQTNHSTFPNPNAPDLEPYQDSIALPWKVLLATFLGLITLGTTLSNAFVIATVSRTRKLHTPANYLIASLAVTDLLVSILVMPISTMYTVTGRWTLGQVVCDFWLSSDITCCTASILHLCVIALDRYWAITDAVYSAKRTPKRAAGMIIMVWVFSVSISMPPLFWR------E--ADCSVN-------TDHILYTVYSTVGAFYFPTLLLIALYGRIYVEARSRKVSLEKKKLMAARERKATRTLGIILGAFIVCWLPFFIISLALPICDDAWFHLAIFDFFNWLGYLNSLINPIIYTKSNDDFKQAFQKLMRFRRTS---------------------
>NY2R_PIG
EEMKMEPSGPGHTTPRGELAPDSEPELKDSTKLIEVQIILILAYCSIILLGVVGNSLVIHVVIKFKSMRTVTNFFIANLAVADLLVNTLCLPFTLTYTLMGEWKMGPVLCHLVPYAQGLAVQVSTITLTVIALDRHRCIVYHL-ESKISKRISFLIIGLAWGISALLASPLAIFREYS-FE--IVACTEKWPGEEKSIYGTVYSLSSLLILYVLPLGIISFSYARIWSKLKNHSPG-GVNDHYHQRRQKTTKMLVCVVVVFAVSWLPLHAFQLAVDIDSQVKEYKLIFTVFHIIAMCSTFANPLLYGWMNSNYRKAFLSAFRCEQRLDAIHSE---AKKNLEATKNGG
>B1AR_.ENLA
WGPMEC--RNRS----GTPTTVPSPMHPLPELTHQWTMGMTMFMAAIILLIVMGNIMVIVAIGRNQRLQTLTNVFITSLACADLIMGLFVVPLGATLVVSGRWLYGSIFCEFWTSVDVLCVTASIETLCVISIDRYIAITSPFYQSLLTKGRAKGIVCSVWGISALVSFLPIMMHWWRDM-K--GCCDFV--------TNRAYAIASSIISFYFPLIIMIFVYIRVFKEAQKQHGRRILSKILVAKEQKALKTLGIIMGTFTLCWLPFFLANVVNVFYRNL-IPDKLFLFLNWLGYANSAFNPIIYCRSP-DFRKAFKRLLCCPKKADRHLHTTGEFVNSLDTN---A
>D2DR_HUMAN
---MDPLNLSWYDDDLERQNWSRPFNGSDGKADRPHYNYYATLLTLLIAVIVFGNVLVCMAVSREKALQTTTNYLIVSLAVADLLVATLVMPWVVYLEVVGEWKFSRIHCDIFVTLDVMMCTASILNLCAISIDRYTAVAMPMNTRYSSKRRVTVMISIVWVLSFTISCPLLFGLN---Q----NECII---------ANPAFVVYSSIVSFYVPFIVTLLVYIKIYIVLRRRRTSMSRRKLSQQKEKKATQMLAIVLGVFIICWLPFFITHILNIHCDCN-IPPVLYSAFTWLGYVNSAVNPIIYTTFNIEFRKAFLKILHC-------------------------
>FML2_HUMAN
-----------METNFSIPLNETEEVLPEPAGHTVLWIFSLLVHGVTFVFGVLGNGLVIWVAGFRMT-RTVNTICYLNLALADFSFSAILPFRMVSVAMREKWPFASFLCKLVHVMIDINLFVSVYLITIIALDRCICVLHPAAQNHRTMSLAKRVMTGLWIFTIVLTLPNFIFWTTISAFWAVERLNVF------ITMAKVFLILHFIIGFTVPMSIITVCYGIIAAKI---------HRNHMIKSSRPLRVFAAVVASFFICWFPYELIGILMAVWLKEKIILVLINPTSSLAFFNSCLNPILYVFMGRNFQERLIRSLPTSLERALTEVPDSATSAS--------
>ACM5_MACMU
-------MEGDSYHNATTVNGTPVYHQPLERHRLWEVISIAAVTAVVSLITIVGNVLVMISFKVNSQLKTVNNYYLLSLACADLIIGIFSMNLYTTYILMGRWALGSLACDLWLALDYVASNASVMNLLVISFDRYFSITRPLYRAKRTPKRAGVMIGLAWLISFILWAPAILCWQYLV-VP-LDECQIQ------FLSEPTITFGTAIAAFYIPVSVMTILYCRIYRETEKRNPSTKRKRMVLVKERKAAQTLSAILLAFIITWTPYNIMVLVSTFCDKC-VPVTLWHLGYWLCYVNSTVNPICYALCNRTFRKTFKMLLLCRWKKKKVEEKLYW------------
>PE22_CANFA
----------------MGSISNNSGSEDCESREWLPSGESPAISSAMFSAGVLGNLIALALLARRWRSISLFHVLVTELVFTDLLGTCLISPVVLASYARNQLEPERRACTYFAFAMTFFSLATMLMLFAMALERYLSIGRPYYQRHVTRRGGLAVLPTIYTVSLLFCSLPLLGYGQYVQYCPGTWCFIR--------HGRTAYLQLYATLLLLLIVAVLACNFSVILNLIRMDGSRRGERVSVAEETDHLILLAIMTITFAICSLPFTIFAYMNETSS---RREKWDLQALRFLSINSIIDPWVFAIFRPPVLRLMRSVLCCRVSLRAQDATQTSSRLTFVDTS---
>5H1B_RABIT
AAGSQIAVPQANLSAAHSHNCSAEGYIYQDSIALPWKVLLVLLLALFTLATTLSNAFVVATVYRTRKLHTPANYLIASLAVTDLLVSILVMPISTMYTVTGRWTLGQVVCDLWLSSDITCCTASIMHLCVIALDRYWAITDAVYSAKRTPKRAAIMIRLVWVFSICISLPPFFWR----E-E--SECLVN-------TDHVLYTVYSTVGAFYLPTLLLIALYGRIYVEARSRRVSLEKKKLMAARERKATKTLGIILGVFIVCWLPFFIISLVMPICKDAWFHQAIFDFFTWLGYVNSLINPIIYTMSNEDFKQAFHKLIRFKCTS---------------------
>GLHR_ANTEL
HSNHTPNGTQFHQCSKIPVQCVPKSDAFHPCEDIMGYVWLTVVSFMVGAVALVANLVVALVLLTSQRRLNVTRFLMCNLAFADFILGLYIFILTSVSAVTRGWQNG-AGCKILGFLAVFSSELSLFTLVMMTIERFYAIVHAMMNARLSFRKTVRFMIGGWIFALVMAVVPLTGVSGY----KVAICLPF-----VSDATSTAYVAFLLLVNGASFISVMYLYSRMLYVVVSG---GDMEGAPKRNDSKVAKRMAILVFTDMLCWAPIAFFGLLAAFGQTLLTVTQSKILLVFFFPINSICNPFLYAFFTKAFKRELFTALSRIGFCKFRALKYNGSRSRRHHSTVNA
>MSHR_OVIMO
PALGSQRRLLGSLNCTPPATLPLTLAPNRTGPQCLEVSIPNGLFLSLGLVSLVENVLVVAAIAKNRNLHSPMYYFVCCLAMSDLLVSVSNVLETAVMLLLEAAAVVQQLDNVIDVLICSSMVSSLCFLGAIAVDRYISIFYALYHSVVTLPRAWRIIAAIWVASILTSVLSITYYY----------------------NNHTVVLLCLVGFFIAMLALMAVLYVHMLARACRHIARKRQRPIHQGFGLKGAATLTILLGVFFLCWGPFFLHLSLIVLCPQHGCIFKNFNLFLALIICNAIVDPLIYAFRSQELRKTLQEVLQCSW-----------------------
>P2YR_CHICK
PELLAG-----------GWAAGNATTKCSLTKTGFQFYYLPTVYILVFITGFLGNSVAIWMFVFHMRPWSGISVYMFNLALADFLYVLTLPALIFYYFNKTDWIFGDVMCKLQRFIFHVNLYGSILFLTCISVHRYTGVVHPLSLGRLKKKNAVYVSSLVWALVVAVIAPILFYSGTG----KTITCYDT-TADEYLRSYFVYSMCTTVFMFCIPFIVILGCYGLIVKALIY------KDLDNSPLRRKSIYLVIIVLTVFAVSYLPFHVMKTLNLRARLDDKVYATYQVTRGLASLNSCVDPILYFLAGDTFRRRLSRATRKSSRRSEPNVQSKSLTEYKQNGDTSL
>ADMR_MOUSE
LEPDNDFRDIHNWTELLHLFNQTFTDCHIEFNENTKHVVLFVFYLAIFVVGLVENVLVICVNCRRSGRVGMLNLYILNMAIADLGIILSLPVWMLEVMLYETWLWGSFSCRFIHYFYLVNMYSSIFFLTCLSIDRYVTLTNTSSWQRHQHRIRRAVCAGVWVLSAIIPLPEVVHIQLLD--E--PMCLFLAPFETYSAWALAVALSATILGFLLPFLLIAVFNILTACRL---------RRQRQTESRRHCLLMWAYIVVFAICWLPYQVTMLLLTLHGTHNLLYFFYEIIDCFSMLHCVANPILYNFLSPSFRGRLLSLVVRYLPKEQARAAGGRQHSIIITKEGSL
>FSHR_BOVIN
FDVMYSEFDYDLCNEVVDVTCSPEPDAFNPCEDIMGDDILRVLIWFISILAITGNILVLVILITSQYKLTVPRFLMCNLAFADLCIGIYLLLIASVDVHTKTWQTG-AGCDAAGFFTVFASELSVYTLTAITLERWHTITHAMLECKVQLRHAASIMLVGWIFAFAVALFPIFGISSY----KVSICLPM-----IDSPLSQLYVMSLLVLNVLAFVVICGCYTHIYLTVRNP------NITSSSSDTKIAKRMAMLIFTDFLCMAPISFFAISASLKVPLITVSKSKILLVLFYPINSCANPFLYAIFTKNFRRDFFILLSKFGCYEVQAQTYRSNFHPRNGHCPPA
>B3AR_RAT
MAPWPHKNGSLAFWSDAPTLDPSAAN---TSGLPGVPWAAALAGALLALATVGGNLLVITAIARTPRLQTITNVFVTSLATADLVVGLLVMPPGATLALTGHWPLGATGCELWTSVDVLCVTASIETLCALAVDRYLAVTNPLYGTLVTKRRARAAVVLVWIVSATVSFAPIMSQWWRVQ-E--RCCSFA--------SNMPYALLSSSVSFYLPLLVMLFVYARVFVVAKRQGVPRRPARLLPLGEHRALRTLGLIMGIFSLCWLPFFLANVLRALVGPS-VPSGVFIALNWLGYANSAFNPLIYCRSP-DFRDAFRRLLCSYGGRGPEEPRVVTSRQNSPLNRFDG
>D3DR_CERAE
--------MAPLSQLSGHLNYTCGVENSTGASQARPHAYYALSYCALILAIVFGNGLVCMAVLKERALQTTTNYLVVSLAVADLLVATLVMPWVVYLEVTGGWNFSRVCCDVFVTLDVMMCTASILNLCAISIDRYTAVVMPVGTGQSSCRRVTLMITAVWVLAFAVSCPLLFGFN---D-P--TVCSI---------SNPDFVIYSSVVSFYLPFGVTVLVYARIYVVLKQRTSLPLQPRGVPLREKKATQMVAIVLGAFIVCWLPFFLTHVLNTHCQTCHVSPELYSATTWLGYVNSALNPVIYTTFNIEFRKAFLKILSC-------------------------
>NK2R_MOUSE
----MGAHASVTDTNILSGLESNATGVTAFSMPGWQLALWATAYLALVLVAVTGNATVIWIILAHERMRTVTNYFIINLALADLCMAAFNATFNFIYASHNIWYFGSTFCYFQNLFPVTAMFVSIYSMTAIAADRYMAIVHPF-QPRLSAPSTKAVIAVIWLVALALASPQCFYST---GA---TKCVVAWPNDNGGKMLLLYHLVVFVLIYFLPLVVMFAAYSVIGLTLWKRPRHHGANLRHLQAKKKFVKAMVLVVVTFAICWLPYHLYFILGTFQEDIKFIQQVYLALFWLAMSSTMYNPIIYCCLNHRFRSGFRLAFRCCPWGTPTEE----HTPSISRRVNRC
>BONZ_MACNE
------MAEHDYHEDYGLNSFNDSSQEEHQDFLQFRKVFLPCMYLVVFVCGLVGNSLVLVISIFYHKLQSLTDVFLVNLPLADLVFVCTLPFWAYAGIHE--WIFGQVMCKTLLGVYTINFYTSMLILTCITVDRFIVVVKATNQQAKRMTWGKVICLLIWVISLLVSLPQIIYGNVFNLD--KLICGYH-----DKEISTVVLATQMTLGFFLPLLAMIVCYSVIIKTL---------LHAGGFQKHRSLKIIFLVMAVFLLTQTPFNLVKLIRSTHWEYTSFHYTIIVTEAIAYLRACLNPVLYAFVSLKFRKNFWKLVKDIGCLPYLGVSHQWKTFSASHNVEAT
>NK1R_MOUSE
-----MDNVLPVDSDLFPNTSTNTSESNQFVQPTWQIVLWAAAYTVIVVTSVVGNVVVIWIILAHKRMRTVTNYFLVNLAFAEACMAAFNTVVNFTYAVHNVWYYGLFYCKFHNFFPIAALFASIYSMTAVAFDRYMAIIHPL-QPRLSATATKVVIFVIWVLALLLAFPQGYYST---SR---VVCMIEWPEHPNRTYEKAYHICVTVLIYFLPLLVIGYAYTVVGITLE---IPSDRYHEQVSAKRKVVKMMIVVVCTFAICWLPFHIFFLLPYINPDLKFIQQVYLASMWLAMSSTMYNPIIYCCLNDRFRLGFKHAFRCCPFISAGDY----STRYLQTQSSVY
>O.2R_HUMAN
ASELNETQEPFLNPTDYDDEEFLRYLWREYLHPKEYEWVLIAGYIIVFVVALIGNVLVCVAVWKNHHMRTVTNYFIVNLSLADVLVTITCLPATLVVDITETWFFGQSLCKVIPYLQTVSVSVSVLTLSCIALDRWYAICHPL-MFKSTAKRARNSIVIIWIVSCIIMIPQAIVMECSTKTTLFTVCDER---WGGEIYPKMYHICFFLVTYMAPLCLMVLAYLQIFRKLWCRKSRVAAEIKQIRARRKTARMLMVVLLVFAICYLPISILNVLKRVFGMFETVYAWFTFSHWLVYANSAANPIIYNFLSGKFREEFKAAFSCCCLGVHHRQEDRLESRKSLTTQISN
>IL8A_RABIT
WTWFEDEFANATGMPPVEKDYSPC----LVVTQTLNKYVVVVIYALVFLLSLLGNSLVMLVILYSRSNRSVTDVYLLNLAMADLLFALTMPIWAVSKEKG--WIFGTPLCKVVSLVKEVNFYSGILLLACISVDRYLAIVHATRTLTQKRHLVKFICLGIWALSLILSLPFFLFRQVFS-NS-SPVCYED-LGHNTAKWRMVLRILPHTFGFILPLLVMLFCYGFTLRTL---------FQAHMGQKHRAMRVIFAVVLIFLLCWLPYNLVLLADTLMRTHNDIDRALDATEILGFLHSCLNPIIYAFIGQNFRNGFLKMLAARGLISKEFLTRHR------------
>AG2S_MOUSE
------MILNSSIEDGIKRIQDDC---PKAGRHNYIFVMIPTLYSIIFVVGIFGNSLVVIVIYFYMKLKTVASVFLLNLALADLCFLLTLPLWAVYTAMEYQWPFGNHLCKIASASVSFNLYASVFLLTCLSIDRYLAIVHPMSRLRRTMLVAKVTCIIIWLMAGLASLPAVIHRNVYF-TN--TVCAFH-YESQNSTLPIGLGLTKNILGFVFPFVIILTSYTLIWKALKKA----YKIQKNTPRNDDIFRIIMAIVLFFFFSWVPHQIFSFLDVLIQLGDVVDTAMPITICIAYFNNCLNPLFYGFLGKKFKRYFLQLLKYIPPKARSHAGLSTRPSD-------N
>OPSD_ANGAN
MNGTEGPNFYIPMSNITGVVRSPFEYPQYYLAEPWAYTILAAYMFTLILLGFPVNFLTLYVTIEHKKLRTPLNYILLNLAVANLFMVFGGFTTTVYTSMHGYFVFGETGCNLEGYFATLGGEISLWSLVVLAIERWVVVCKPMSNFRFGENHAIMGLAFTWIMANSCAMPPLFGWSRYIPEGMQCSCGVDYYTLKPEVNNESFVIYMFIVHFSVPLTIISFCYGRLVCTVKEAAAQQQESETTQRAEREVTRMVVIMVIAFLVCWVPYASVAWYIFTHQGSTFGPVFMTVPSFFAKSSAIYNPLIYICLNSQFRNCMITTLFCGKNPFQEEEGASTASSVSS--VSPA
>O3A3_HUMAN
----MSLQKLMEPEAGTNRTAVAEFILLGLVQTEEMQPVVFVLLLFAYLVTIGGNLSILAAVLVEPKLHAPMYFFLGNLSVLDVGCITVTVPAMLGRLLSHKSTISYDACLSQLFFFHLLAGMDCFLLTAMAYDRLLAICQPLYSTRMSQTVQRMLVAASWACAFTNALTHTVAMSTL---NNHFYCDLPQLSCSSTQLNELLLFVAAAFMAVAPLVFISVSYAHVVAAV--------LQIRSAEGRKKAFSTCGSHLTVVGIFYGTGVFSYMRLGSVE---SSDKDKGVGVFMTVINPMLNPLIYSLRNTDVQGALCQLLVGERSLT--------------------
>GPRC_MOUSE
NLSGLPRDCIDAGAPENISAAVPSQGSVAESEPELVVNPWDIVLCSSGTLICCENAVVVLIIFHSPSLRAPMFLLIGSLALADLLAG-LGLIINFVFA-Y--LLQSEATKLVTIGLIVASFSASVCSLLAITVDRYLSLYYALYHSERTVTFTYVMLVMLWGTSICLGLLPVMGWNCL-R--DESTCSVV-------RPLTKNNAAILSISFLFMFALMLQLYIQICKIVMRHIALHFLATSHYVTTRKGVSTLALILGTFAACWMPFTLYSLIADY----TYPSIYTYATLLPATYNSIINPVIYAFRNQEIQKALCLICCGCIPSSLSQRARSP------------
>D1DR_CARAU
-------------MAVLDLNLTTVIDSGFMESDRSVRVLTGCFLSVLILSTLLGNTLVCAAVTKFRHRSKVTNFFVISLAVSDLLVAVLVMPWKAVTEVAGFWPFG-AFCDIWVAFDIMCSTASILNLCVISVDRYWAISSPFYERKMTPRVAFVMISGAWTLSVLISFIPVQLKWHKAVN-PTDNCDSS--------LNRTYAISSSLISFYIPVAIMIVTYTQIYRIAQKQGSNESSFKLSFKRETKVLKTLSVIMGVFVCCWLPFFILNCMVPFCKRTCISPTTFDVFVWFGWANSSLNPIIYAFNA-DFRRAFAILLGCQRLCPGSI-SMET------------
>OPS1_SCHGR
GG-F-ANQTVVDKVPPEMLYLVDPHWYQFPPMNPLWHGLLGFVIGVLGVISVIGNGMVIYIFSTTKSLRTPSNLLVVNLAFSDFLMMFTMSAPMGINCYYETWVLGPFMCELYALFGSLFGCGSIWTMTMIALDRYNVIVKGLSAKPMTNKTAMLRILFIWAFSVAWTIMPLFGWNRYVPEGNMTACGTD--YLTKDWVSRSYILVYSFFVYLLPLGTIIYSYFFILQAVSAHMNVRSAEASQTSAECKLAKVALMTISLWFFGWTPYLIINFTGIFETMK-ISPLLTIWGSLFAKANAVFNPIVYGISHPKYRAALEKKFPSLACASSSDDNTSVSDEKSEKSASA-
>MSHR_CAPCA
PVLGSQRRLLGSLNCTPPATFPLTLAPNRTGPQCLEVSIPDGLFLSLGLVSLVENVLVVAAIAKNRNLHSPMYYFICCLAVSDLLVSVSNVLETAVMLLLEAAAVVQQLDNVIDMLICGSMVSSLCFLGAIAVDRYISIFYALYHSVVTLPRAWRIIAAIWVASILTSLLFITYYY----------------------NNHTVVLLCLVGFFIAMLALMAVLYVHMLARACQHIARKRQRPIHQGFGLKGAATLTILLGVFFLCWGPFFLHLSLIVLCPQHGCIFKNFNLFLALIICNAIVDPLIYAFRSQELRKTLQEVLQCSW-----------------------
>O3A4_HUMAN
-----------MDLGNSGNDSVTKFVLLGLTETAALQPILFVIFLLAYVTTIGGTLSILAAILMETKLHSPMYFFLGNLSLPDVGCVSVTVPAMLSHFISNDRSIPYKACLSELFFFHLLAGADCFLLTIMAYDRYLAICQSLYSSRMSWGIQQALVGMSWVFSFTNALTQTVALSPL---NNHFYCDLPQLSCASVHLNGQLLFVAAAFMGVAPLVLITVSYAHVAAAV--------LRIRSAEGKKKAFSTCSSHLTVVGIFYGTGVFSYTRLGSVE---SSDKDKGIGILNTVISPMLNPLIYWTSLLDVGCISHCSSDAGVSPGPPVQSSLCLSPPPGWGGLSP
>P70658
IYTSDNYSEEVG-SGDYDSNKEPC---FRDENVHFNRIFLPTIYFIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVADLLFVITLPFWAVDAMAD--WYFGKFLCKAVHIIYTVNLYSSVLILAFISLDRYLAIVHATSQRPRKLLAEKAVYVGVWIPALLLTIPDFIFADVSQGDD-RYICDRL---YPDSLWMVVFQFQHIMVGLILPGIVILSCYCIIISKL---------SHSKGHQKRKALKTTVILILAFFACWLPYYVGISIDSFILLGSIVHKWISITEALAFFHCCLNPILYAFLGAKFKSSAQHALNSMSRGSSLKILSKG----------GH
>OPSD_AMBTI
MNGTEGPNFYVPFSNKSGVVRSPFEYPQYYLAEPWQYSVLAAYMFLLILLGFPVNFLTLYVTIQHKKLRTPLNYILLNLAFANHFMVFGGFPVTMYSSMHGYFVFGQTGCYIEGFFATMGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVMMTWIMALACAAPPLFGWSRYIPEGMQCSCGVDYYTLKPEVNNESFVIYMFLVHFTIPLMIIFFCYGRLVCTVKEAAAQQQESATTQKAEKEVTRMVIIMVVAFLICWVPYASVAFYIFSNQGTDFGPIFMTVPAFFAKSSAIYNPVIYIVLNKQFRNCMITTICCGKNPFGDDETTSAASSVSSSQVSPA
>V2R_PIG
-MLRATTSAVPRALSWPAAPGNGSEREPLDDRDPLLARVELALLSTVFVAVALSNGLVLGALVRRGRRWAPMHVFIGHLCLADLAVALFQVLPQLAWDATYRFRGPDALCRAVKYLQMVGMYASSYMILAMTLDRHRAICRPMAYRHGGGARWNRPVLVAWAFSLLLSLPQLFIFAQRDSG--VLDCWAS---FAEPWGLRAYVTWIALMVFVAPALGIAACQVLIFREIHTSGRRPREGARVSAAMAKTARMTLVIVAVYVLCWAPFFLVQLWSVWDPKA-REGPPFVLLMLLASLNSCTNPWIYASFSSSISSELRSLLCCPRRRTPPSLRPQESFSARDTSS---
>5H7_HUMAN
GSWAPHLLS---EVTASPAPTNASGCGEQINYGRVEKVVIGSILTLITLLTIAGNCLVVISVCFVKKLRQPSNYLIVSLALADLSVAVAVMPFVSVTDLIGGWIFGHFFCNVFIAMDVMCCTASIMTLCVISIDRYLGITRPLYPVRQNGKCMAKMILSVWLLSASITLPPLFGWAQ--N-D--KVCLIS--------QDFGYTIYSTAVAFYIPMSVMLFMYYQIYKAARKSSRLERKNISIFKREQKAATTLGIIVGAFTVCWLPFFLLSTARPFICGTCIPLWVERTFLWLGYANSLINPFIYAFFNRDLRTTYRSLLQCQYRNINRKLSAAGAERPERPEFVLR
>NTR2_RAT
-----METSSPWPPRPSPSAGLSLEARLGVDTRLWAKVLFTALYSLIFAFGTAGNALSVHVVLKARARPGRLRYHVLSLALSALLLLLVSMPMELYNFVWSHWVFGDLGCRGYYFVRELCAYATVLSVASLSAERCLAVCQPLARRLLTPRRTRRLLSLVWVASLGLALPMAVIMGQKH--AASRVCTVL---VSRATLQVFIQVN-VLVSFALPLALTAFLNGITVNHLMALVQARHKDASQIRSLQHSAQVLRAIVAVYVICWLPYHARRLMYCYIPDDDFYHYFYMVTNTLFYVSSAVTPILYNAVSSSFRKLFLESLGSLCGEQHSLVPLPQSTYSFRLWGSPR
>OPSD_RAJER
MNGTEGENFYVPMSNKTGVVRSPFDYPQYYLGEPWMFSALAAYMFFLILTGLPVNFLTLFVTIQHKKLRQPLNYILLNLAVSDLFMVFGGFTTTIITSMNGYFIFGPAGCNFEGFFATLGGEVGLWCLVVLAIERYMVVCKPMANFRFGSQHAIIGVVFTWIMALSCAGPPLVGWSRYIPEGLQCSCGVDYYTMKPEVNNESFVIYMFVVHFTIPLIVIFFCYGRLVCTVKEAAAQQQESESTQRAEREVTRMVIIMVVAFLICWVPYASVAFYIFINQGCDFTPFFMTVPAFFAKSSAVYNPLIYILMNKQFRNCMITTICLGKNPFEEEESTSAASSVSSSQVAPA
>OPSD_CAMAB
GGGF-GNQTVVDKVPPEMLHMVDAHWYQFPPMNPLWHALLGFVIGVLGVISVIGNGMVIYIFTTTKSLRTPSNLLVVNLAISDFLMMLCMSPAMVINCYYETWVLGPLFCELYGLAGSLFGCASIWTMTMIAFDRYNVIVKGLSAKPMTINGALIRILTIWFFTLAWTIAPMFGWNRYVPEGNMTACGTD--YLTKDLFSRSYILIYSIFVYFTPLFLIIYSYFFIIQAVAAHMNVRSAENQSTSAECKLAKVALMTISLWFMAWTPYLVINYSGIFETTK-ISPLFTIWGSLFAKANAVYNPIVYGISHPKYRAALFQKFPSLACTTEPTG-ADTTEGNEKPAA---
>SSR4_MOUSE
VED----TTWTPGINASWAPEQEEDAMGSDGTGTAGMVTIQCIYALVCLVGLVGNALVIFVILRYAKMKTATNIYLLNLAVADELFMLSVPFVRSAAALRH-WPFGAVLCRAVLSVDGLNMFTSVFCLTVLSVDRYVAVVHPLTATYRRPSVAKLINLGVWLASLLVTLPIAVFADTRPGGE-AVACNLH---WPHPAWSAVFVIYTFLLGFLPPVLAIGLCYLLIVGKMRAVALR-GGWQQRRRSEKKITRLVLMVVTVFVLCWMPFYVVQLLNLFVTS--LDATVNHVSLILSYANSCANPILYGFLSDNFRRSFQRVLCLRCCLLETTGGAEETALKSRGGAGCI
>NY1R_.ENLA
NFSTYFENLSVPNNISG---NITFPISEDCALPLPMIFTLALAYGAVIILGLSGNLALIIIILKQKEMRNVTNILIVNLSFSDLLATIMCLPFTLIYTLMDHWIFGEVMCKLNEYIQCVSVTVSIFSLVLIAIERHQLIINPR-GWRPNNRHACFGITVIWGFAMACSTPLMMYSVLTDIG--KYVCLED---FPEDKFRLSYTTLLFILQYLGPLCFIFVCYTKIFLRLKRRNMMIRDNKYRSSETKRINIMLLSIVVGFALCWLPFFIFNLVFDWNHEACNHNLLFLICHLTAMISTCVNPIFYGFLNKNFQRDLQFFFNFCDFRSREDDY---TMHTDVSKTSLK
>LSHR_PIG
ESELSDWDYDYGFCSPKTLQCAPEPDAFNPCEDIMGYDFLRVLIWLINILAIMGNVTVLFVLLTSHYKLTVPRFLMCNLSFADFCMGLYLLLIASVDAQTKGWQTG-NGCSVAGFFTVFASELSVYTLTVITLERWHTITYAILDQKLRLRHAIPIMLGGWLFSTLIAMLPLVGVSSY----KVSICLPM-----VETTLSQVYILTILILNVVAFIIICACYIKIYFAVQNP------ELMATNKDTKIAKKMAVLIFTDFTCMAPISFFAISAALKVPLITVTNSKVLLVLFYPVNSCANPFLYAIFTKAFRRDFFLLLSKSGCCKHQAELYRRK---NGFTGSNK
>GPR4_HUMAN
--------------------MGNHTWEGCHVDSRVDHLFPPSLYIFVIGVGLPTNCLALWAAYRQVQQRNELGVYLMNLSIADLLYICTLPLWVDYFLHHDNWIHGPGSCKLFGFIFYTNIYISIAFLCCISVDRYLAVAHPLFARLRRVKTAVAVSSVVWATELGANSAPLFHDEL--NH---TFCFEKPMEGWVAWMVAWMNLYRVFVGFLFPWALMLLSYRGILRAVRG------SVSTERQEKAKIKRLALSLIAIVLVCFAPYHVLLLSRSAIYLGERVFSAYHSSLAFTSLNCVADPILYCLVNEGARSDVAKALHNLLRFLASDKPQEMETPLTSKRNSTA
>OPSD_POMMI
MNGTEGPFFYIPMVNTTGIVRSPYEYPQYYLVNPAAYAALGAYMFFLILTGFPINFLTLYVTLEHKKLRTALNLILLNLAVADLFMVFGGFTTTMYTSMHGYFVLGRLGCNVEGFFATLGGEIALWSLVVLAVERWVVVCKPISNFRFTENHAIMGVAFSWIMAATCAVPPLVGWSRYIPEGMQCSCGVDYYTRAEGFNNESFVIYMFIVHFLAPLIVIFFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVIIMVIGFLTSWLPYASVAWYIFTHQGTEFGPLFMTIPAFFAKSSALYNPMIYICMNKQFRHCMITTLCCGKNPFEEEEGASTASSVSSSSVSPA
>CB1R_POEGU
EFYNKSLSTFKDNEENIQCGENFMDMECFMILNPSQQLAIAVLSLTLGTFTVLENLLVLCVILHSRSRCRPSYHFIGSLAVADLLGSVIFVYSFVDFHVFHR-KDSPNVFLFKLGGVTASFTASVGSLFLTAIDRYISIHRPLYKRIVTRPKAVVAFCVMWTIAIVIAVLPLLGWNCK-K--LNSVCSDI------FPLIDETYLMFWIGVTSILLLFIVYAYMYILWKAHSHTEDQITRPDQTRMDIRLAKTLVLILVVLIICWGPLLAIMVYDVFGKMNKLIKTIFAFCSMLCLLNSTVNPIIYALRSKDLRHAFRSMFPTCEGTAQPLDNSMEANN-AGNVHRAA
>ADMR_HUMAN
AVPTSDLGEIHNWTELLDLFNHTLSECHVELSQSTKRVVLFALYLAMFVVGLVENLLVICVNWRGSGRAGLMNLYILNMAIADLGIVLSLPVWMLEVTLDYTWLWGSFSCRFTHYFYFVNMYSSIFFLVCLSVDRYVTLTSASSWQRYQHRVRRAMCAGIWVLSAIIPLPEVVHIQLVE--E--PMCLFMAPFETYSTWALAVALSTTILGFLLPFPLITVFNVLTACRL---------RQPGQPKSRRHCLLLCAYVAVFVMCWLPYHVTLLLLTLHGTHHLLYFFYDVIDCFSMLHCVINPILYNFLSPHFRGRLLNAVVHYLPKDQTKAGTCAQHSIIITKGDSQ
>GPR6_HUMAN
G-PPAAAALGAGGGANGSLELSSQLSAGPPGLLLPAVNPWDVLLCVSGTVIAGENALVVALIASTPALRTPMFVLVGSLATADLLAG-CGLILHFVFQ-Y--LVPSETVSLLTVGFLVASFAASVSSLLAITVDRYLSLYNALYYSRRTLLGVHLLLAATWTVSLGLGLLPVLGWNCL-A--ERAACSVV-------RPLARSHVALLSAAFFMVFGIMLHLYVRICQVVWRHIALHCLAPPHLAATRKGVGTLAVVLGTFGASWLPFAIYCVVGSH----EDPAVYTYATLLPATYNSMINPIIYAFRNQEIQRALWLLLCGCFQSKVPFRSRSP------------
>5H1D_MOUSE
SLPNQSLEGLPQEASNRSLNAT---GAWDPEVLQALRISLVVVLSVITLATVLSNAFVLTTILLTKKLHTPANYLIGSLATTDLLVSILVMPISIAYTTTRTWNFGQILCDIWVSSDITCCTASILHLCVIALDRYWAITDALYSKRRTAGHAAAMIAAVWIISICISIPPLFWR----H-E--SDCLVN-------TSQISYTIYSTCGAFYIPSILLIILYGRIYVAARSRKLALERKRISAARERKATKTLGIILGAFIICWLPFFVVSLVLPICRDSWIHPALFDFFTWLGYLNSLINPVIYTVFNEDFRQAFQKVVHFRKAS---------------------
>OPSD_ORYLA
MNGTEGPYFNVPMVNTTGIVRSPYEYPQYYLVSPAAYAALGAYMFFLILVGFPINFLTLYVTLEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVLGRLGCNLEGFFATLGGEIGLWSLVVLAIERWVVVCKPISNFRFGENHAIMGLVFTWIMAASCAVPPLVGWSRYIPEGMQCSCGVDYYTRAEGFNNESFVVYMFVCHFLIPLIVVFFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVVIMVIGFLVCWLPYASVAWYIFTNQGSEFGPLFMTIPAFFAKSSSIYNPAIYICMNKQFRNCMITTLCCGKNPFEEEEGASTASSVSSSSVSPA
>YDBM_CAEEL
EGAGEDVDHHSLFCPKKLVGNLKGFIRNQYHQHETIQILKGSALFLLVLWTIFANSLVFIVLYKNPRLQTVPNLLVGNLAFSDLALGLIVLPLSSVYAIAGEWVFPDALCEVFVSADILCSTASIWNLSIVGLDRYWAITSPVYMSKRNKRTAGIMILSVWISSALISLAPLLGWKQ--Q-TTVRQCTFL--------DLPSYTVYSATGSFFIPTLLMFFVYFKIYQAFAKHYRKRKPKAISAAKERRGVKVLGIILGCFTVCWAPFFTMYVLVQFCKDCSPNAHIEMFITWLGYSNSAMNPIIYTVFNRDYQIALKRLFTSEKKPSSTSRV---------------
>THRR_MOUSE
LLEGRAVYLNISLPPHTPPPPFISEDASGYLTSPWLTLFMPSVYTIVFIVSLPLNVLAIAVFVLRMKVKKPAVVYMLHLAMADVLFVSVLPFKISYYFSGTDWQFGSGMCRFATAAFYGNMYASIMLMTVISIDRFLAVVYPISLSWRTLGRANFTCVVIWVMAIMGVVPLLLKEQTTR--N--TTCHDVLSENLMQGFYSYYFSAFSAIFFLVPLIVSTVCYTSIIRCL------SSSAVANRSKKSRALFLSAAVFCIFIVCFGPTNVLLIVHYLFLSDEAAYFAYLLCVCVSSVSCCIDPLIYYYASSECQRHLYSILCCKESSDPNSCNSTGDTCS--------
>H218_RAT
----MGGLYSEYLNPEKVQEHYNYTKETLDMQETPSRKVASAFIIILCCAIVVENLLVLIAVARNSKFHSAMYLFLGNLAASDLLAG-VAFVANTLLSGPVTLSLTPLQWFAREGSAFITLSASVFSLLAIAIERQVAIAKVKLYGSDKSCRMLMLIGASWLISLILGGLPILGWNCL-D--HLEACSTV------LPLYAKHYVLCVVTIFSVILLAIVALYVRIYFVVRSS-----HADVAGPQTLALLKTVTIVLGVFIICWLPAFSILLLDSTCPVRCPVLYKAHYFFAFATLNSLLNPVIYTWRSRDLRREVLRPLLCWRQGKGATGRRGGPLRSSSSLERGL
>CKR5_MOUSE
--MDFQGSVPTYIYDIDYGMSAPC---QKINVKQIAAQLLPPLYSLVFIFGFVGNMMVFLILISCKKLKSVTDIYLLNLAISDLLFLLTLPFWAHYAANE--WIFGNIMCKVFTGVYHIGYFGGIFFIILLTIDRYLAIVHAVALKVRTVNFGVITSVVTWVVAVFASLPEIIFTRSQK-EG-HYTCSPHFPHTQYHFWKSFQTLKMVILSLILPLLVMIICYSGILHTL--------FRCRNEKKRHRAVRLIFAIMIVYFLFWTPYNIVLLLTTFQEFFNRLDQAMQATETLGMTHCCLNPVIYAFVGEKFRSYLSVFFRKHIVKRFCKRCSIFV-------SSVY
>B3AR_CANFA
MAPWPHGNGSVASWPAAPTPTPDAANTSGLPGAPWAVALAGALLALEVLATVGGNLLVIVAIARTPRLQTMTNVFVTSLATADLVVGLLVVPPGATLALTGRWPLGATGCELWTSVDVLCVTASIETLCALAVDRYLAVTNPLYGALVTKRRARAAVVLVWVVSAAVSFAPIMSKWWRVQ-R--HCCAFA--------SNIPYALLSSSVSFYLPLLVMLFVYARVFLVATRQAVPLRPARLLPLREHRALRTLGLIVGTFTLCWLPFFVANVMRALGGPS-VPSPALLALNWLGYANSAFNPLIYCRSP-DFRSAFRRLLCRCR---REEHRAAAAAPAALTSPAES
>P2Y5_CHICK
-----------------------MVSSNCSTEDSFKYTLYGCVFSMVFVLGLIANCVAIYIFTFTLKVRNETTTYMLNLAISDLLFVFTLPFRIYYFVVRN-WPFGDVLCKISVTLFYTNMYGSILFLTCISVDRFLAIVHPFSKTLRTKRNARIVCVAVWITVLAGSTPASFFQSTNR---EQRTCFENFPESTWKTYLSRIVIFIEIVGFFIPLILNVTCSTMVLRTLNKP----LTLSRNKLSKKKVLKMIFVHLVIFCFCFVPYNITLILYSLMRTQTAVRTMYPVTLCIAVSNCCFDPIVYYFTSDTNSELDKKQQVHQNT----------------------
>GPRO_HUMAN
HASRMSVLRAKPMSNSQRLLLLSPGSPPRTGSISYINIIMPSVFGTICLLGIIGNSTVIFAVVKKSKCNNVPDIFIINLSVVDLLFLLGMPFMIHQLMGNGVWHFGETMCTLITAMDANSQFTSTYILTAMAIDRYLATVHPISTKFRKPSVATLVICLLWALSFISITPVWLYARLIP-FP-AVGCGIR--LPNPDTDLYWFTLYQFFLAFALPFVVITAAYVRILQRMTSSV--PASQRSIRLRTKRVTRTAIAICLVFFVCWAPYYVLQLTQLSISRPTFVYLYNAAIS-LGYANSCLNPFVYIVLCETFRKRLVLSVKPAAQGQLRAVSNAQESKGT-------
>THRR_.ENLA
SGEGSGDQAPVSRSARKPIRRNITKEAEQYLSSQWLTKFVPSLYTVVFIVGLPLNLLAIIIFLFKMKVRKPAVVYMLNLAIADVFFVSVLPFKIAYHLSGNDWLFGPGMCRIVTAIFYCNMYCSVLLIASISVDRFLAVVYPMSLSWRTMSRAYMACSFIWLISIASTIPLLVTEQTQK--D--TTCHDVLDLKDLKDFYIYYFSSFCLLFFFVPFIITTICYIGIIRSL------SSSSIENSCKKTRALFLAVVVLCVFIICFGPTNVLFLTHYLQEANEFLYFAYILSACVGSVSCCLDPLIYYYASSQCQRYLYSLLCCRKVSEPGSSTGQLDNCS--------
>GHSR_HUMAN
SEEPGFNLTLADLDWDASPGNDSLGDELLQLFPAPLLAGVTATCVALFVVGIAGNLLTMLVVSRFRELRTTTNLYLSSMAFSDLLIFLCMPLDLVRLWQYRPWNFGDLLCKLFQFVSESCTYATVLTITALSVERYFAICFPLAKVVVTKGRVKLVIFVIWAVAFCSAGPIFVLVG---T-D-TNECRPT---FAVRSGLLTVMVWVSSIFFFLPVFCLTVLYSLIGRKLWRRGD-VVGASLRDQNHKQTVKMLAVVVFAFILCWLPFHVGRYLFSKSFEPQISQYCNLVSFVLFYLSAAINPILYNIMSKKYRVAVFRLLGFEPFSQRKLSTLKDESSINT------
>CKR3_CAVPO
PEEAELETEFPGTTFYDYEFAQPC---FKVSITDLGAQFLPSLFSLVFIVGLLGNITVIVVLTKYQKLKIMTNIYLLNLAISDLLFLFTLPFWTYYVHWNK-WVFGHFMCKIISGLYYVGLFSEIFFIILLTIDRYLAIVHAVALRTRTVTFGIITSVITWVLAVLAALPEFMFYGTQG-HF-VLFCGPSYPEKKEHHWKRFQALRMNIFGLALPLLIMIICYTGIIKTL---------LRCPSKKKYKAIRLIFVIMVVFFVFWTPYNLLLLFSAFDLSFKQLDMAKHVTEVIAHTHCCINPIIYAFVGERFQKYLRHFLHRNVTMHLSKYIPFFS-------SSIS
>ADMR_RAT
LAPDNDFREIHNWTELLHLFNQTFSDCRMELNENTKQVVLFVFYLAIFVVGLVENVLVICVNCRRSGRVGMLNLYILNMAVADLGIILSLPVWMLEVMLYETWLWGSFSCRFIHYFYLANMYSSIFFLTCLSIDRYVTLTNTSSWQRHQHRIRRAVCAGVWVLSAIIPLPEVVHIQLLD--E--PMCLFLAPFETYSAWALAVALSATILGFLLPFPLIAVFNILSACRL---------RRQGQTESRRHCLLMWAYIVVFAICWLPYHVTMLLLTLHTTHNFLYFFYEITDCFSMLHCVANPILYNFLSPSFRGRLLSLVVRYLPKEQARAAGGRQHSIIITKEGSL
>GPR._MOUSE
--------MDLINSSTHVINVSTSLTNSTGVPTPAPKTIIAASLFMAFIIGVISNGLYLWMLQFKMQ-RTVNTLLFFHLILSYFISTLILPFMATSFLQDNHWVFGSVLCKAFNSTLSVSMFASVFFLSAISVARYYLILHPVSQQHRTPHWASRIALQIWISATILSIPYLVFRTTHD-IS--DWESKE-HQTLGQWIHAACFVGRFLLGFLLPFLVIIFCYKRVATKM---------KEKGLFKSSKPFKVMVTAVISFFVCWMPYHVHSGLVLTKSQP-PLHLTLGLAVVTISFNTVVSPVLYLFTGENFKV-FKKSILALFNSTFSDISSTEETEI--------
>CML2_HUMAN
APNTTSPELNLSHPLLGTALANGTGELSEHQQYVIGLFLSCLYTIFLFPIGFVGNILILVVNISFREKMTIPDLYFINLAVADLILVADSLIEVFNLHER--YYDIAVLCTFMSLFLQVNMYSSVFFLTWMSFDRYIALARAMCSLFRTKHHARLSCGLIWMASVSATLVPFTAVHLQ-A----CFCFA---------DVREVQWLEVTLGFIVPFAIIGLCYSLIVRVLVR----AHRHRGLRPRRQKALRMILAVVLVFFVCWLPENVFISVHLLQRTQHAHPLTGHIVNLAAFSNSCLNPLIYSFLGETFRDKLRLYIEQKTNLPALNRFCHADSTEQSDVRFSS
>AA2A_CANFA
-------------------------------MSTMGSWVYITVELAIAVLAILGNVLVCWAVWLNSNLQNVTNYFVVSLAAADIAVGVLAIPFAITISTG--FCAACHNCLFFACFVLVLTQSSIFSLLAIAIDRYIAIRIPLYNGLVTGTRAKGIIAVCWVLSFAIGLTPMLGWNNCSS-Q-QVACLFE-----DVVPMNYMVYYNFFAFVLVPLLLMLGVYLRIFLAARRQESQGERARSTLQKEVHAAKSLAIIVGLFALCWLPLHIINCFTFFCPECHAPLWLMYLTIVLSHTNSVVNPFIYAYRIREFRQTFRKIIRSHVLRRREPFKAGGAHGSDGEQISLR
>NK2R_MESAU
----MGGRAIVTDTNIFSGLESNTTGVTAFSMPAWQLALWATAYLGLVLVAVTGNATVIWIILAHERMRTVTNYFIINLALADLCMAAFNATFNFVYASHNIWYFGRAFCYFQNLFPITAMFVSIYSMTAIAADRYMAIVHPF-QPRLSAPITKATIAGIWLVALALASPQCFYST---GA---TKCVVAWPNDNGGKMLLLYHLVVFVLVYFLPLVVMFVAYSVIGLTLWKRPRHHGANLRHLHAKKKFVKAMVLVVLTFAICWLPYHLYFILGSFQKDIKFIQQVYLALFWLAMSSTMYNPIIYCCLNHRFRSGFRLAFRCCPWVTPTEE----RTPSLSRRVNRC
>PAFR_MOUSE
----------------------MEHNGSFRVDSEFRYTLFPIVYSVIFILGVVANGYVLWVFANLYPKLNEIKIFMVNLTMADLLFLITLPLWIVYYYNEGDWILPNFLCNVAGCLFFINTYCSVAFLGVITYNRYQAVAYPITAQATTRKRGISLSLIIWVSIVATASYFLATDSTN-G--NITRCFEH---YEPYSVPILVVHVFIAFCFFLVFFLIFYCNLVIIHTLLTQP---MRQQRKAGVKRRALWMVCTVLAVFIICFVPHHVVQLPWTLAELGQAINDAHQITLCLLSTNCVLDPVIYCFLTKKFRKHLSEKFYSMRSSRKCSRATSDPANQTPIVSLKN
>GU27_RAT
--------------------------------MILNCNPFSGLFLSMYLVTVLGNLLIILAVSSNSHLHNLMYFFLSNLSFVDICFISTTIPKMLVNIHSQTKDISYIECLSQVYFLTTFGGMDNFLLTLMACDRYVAICHPLYTVIMNLQLCALLILMFWLIMFCVSLIHVLLMNEL---NPHFFCELAKVANSDTHINNVFMYVVTSLLGLIPMTGILMSYSQIASSL--------LKMSSSVSKYKAFSTCGSHLCVVSLFYGSATIVYFCSSVLHS---THKKMIASLMYTVISPMLNPFIYSLRNKDVKGALGKLFIRVASCPLWSKDFRPRQSL--------
>ET3R_.ENLA
GNVLNMSPP------------PPSPCLSRAKIRHAFKYVTTILSCVIFLVGIVGNSTLLRIIYKNKCMRNGPNVLIASLALGDLFYILIAIPIISISFWLS----TGHSEYIYQLVHLYRARVYSLSLCALSIDRYRAVASWNIRSIGIPVRKAIELTLIWAVAIIVAVPEAIAFNLVELV--ILVCMLPQTSDFMRFYQEVKVWWLFGFYFCLPLACTGVFYTLMSCEMLSINGM-IALNDHMKQRREVAKTVFCLVVIFALCWLPLHVSSIFVRLSATVQLLMVMNYTGINMASLNSCIGPVALYFVSRKFKNCFQSCLCCWCHRPTLTITPMDWKANGHDLDLDR
>B3AR_CAPHI
MAPWPPRNSSLTPWPDIPTLAPNTANASGLPGVPWAVALAGALLALAVLAIVGGNLLVIVAIARTPRLQTMTNVFVTSLATADLVVGLLVVPPGATLALTGHWPLGVTGCELWTSVDVLCVTASIETLCALAVDRYLAVTNPLYGALVTKRRARAAVVLVWVVSAAVSFAPIMSKWWRVQ-R--RCCTFA--------SNMPYALLSSSVSFYLPLLVMLFVYARVFVVATRQGVPRRPARLLPLREHRALRTLGLIMGTFTLCWLPFFVVNVVRALGGPS-VSGPTFLALNWLGYANSAFNPLIYCRSP-DFQSAFRRLLCRCRPEEHLAAASPPRVLTSPAGPRQP
>GP41_HUMAN
-----------------------MDTGPDQSYFSGNHWFVFSVYLLTFLVGLPLNLLALVVFVGKLQRPVAVDVLLLNLTASDLLLLLFLPFRMVEAANGMHWPLPFILCPLSGFIFFTTIYLTALFLAAVSIERFLSVAHPLYKTRPRLGQAGLVSVACWLLASAHCSVVYVIEFSGD--T--GTCYLE-FRKDQLAILLPVRLEMAVVLFVVPLIITSYCYSRLVWILGR--------GGSHRRQRRVAGLLAATLLNFLVCFGPYNVSHVVGYICG---ESPAWRIYVTLLSTLNSCVDPFVYYFSSSGFQADFHELLRRLCGLWGQWQQESSGGEEQRADRPAE
>O.YR_RAT
GTPAANWSVELDLGSGVPPGEEGNRTAGPPQRNEALARVEVAVLCLILFLALSGNACVLLALRTTRHKHSRLFFFMKHLSIRDLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPL--RSLRRRTDRLAVLGTWLGCLVASAPQVHIFSLRE----VFDCWAV---FIQPWGPKAYVTWITLAVYIVPVIVLAACYGLISFKIWQNRAAVSSVKLISKAKIRTVKMTFIIVLAFIVCWTPFFFVQMWSVWDVNA-KEASAFIIAMLLASLNSCCNPWIYMLFTGHLFHELVQRFFCCSARYLKGSRPGENSSTFVLSRRSS
>OLF2_CHICK
-------------MASGNCTTPTTFILSGLTDNPRLQMPLFMVFLVIYTTTLLTNLGLIALIGMDLHLQTPMYIFLQNLSFTDAAYSTVITPKMLATFLEERRTISYVGCILQYFSFVLLTTSEWLLLAVMAYDRYVAICKPLYPSIMTKAVCWRLVKGLYSLAFLNSLVHTSGLLKL---SNHFFCDNRQISSSSTTLNELLVIISGSLFVMSSIITILISYVFIILTV--------VMIRSKDGKYKAFSTCTSHLMAVSLFHGTVIFMYLRSVKLF---SLDTDKIASLFYTVVIPMLNPLIYSWRNKEVKDALRRLTATSVWLH--------------------
>NTR1_MOUSE
PHPQFGLETMLLALSLSNGSGLEPNSNLDVNTDIYSKVLVTAVYLALFVVGTVGNSVTAFTLARKKSLQSTVHYHLGSLALSDLLILLLAMPVELYNFIWVHWAFGDAGCRGYYFLRDACTYATALNVASLSVERYLAICHPFAKTLMSRSRTKKFISAIWLASALLAVPMLFTMGLQN--SGGLVCTPT---VVDTATVKVVIQVNTFMSFLFPMLIISILNTVIANKLTVMEHSMSIEPGRVQALRHGVLVLRAVVIAFVVCWLPYHVRRLMFCYISDEDFYHYFYMLTNALFYVSSAINPILYNLVSANFRQVFLSTLACLCPGWRRRRKKRPSMSSNHAFSTSA
>OLF7_RAT
------------MERRNHSGRVSEFVLLGFPAPAPLRVLLFFLSLL-YVLVLTENMLIIIAIRNHPTLHKPMYFFLANMSFLEIWYVTVTIPKMLAGFIGSKQLISFEACMTQLYFFLGLGCTECVLLAVMAYDRYVAICHPLYPVIVSSRLCVQMAAGSWAGGFGISMVKVFLISRL---SNHFFCDVSNLSCTDMSTAELTDFVLAIFILLGPLSVTGASYMAITGAV--------MRIPSAAGRHKAFSTCASHLTVVIIFYAASIFIYARPKALS---AFDTNKLVSVLYAVIVPLFNPIIYCLRNQDVKRALRRTLHLAQDQEANTNKGSK------------
>MC5R_RAT
-MNSSSHLTLLDLTLNASEDNILGQNVNNKSSACEDMGIAVEVFLTLGLVSLLENILVIGAIVKNKNLHSPMYFFVGSLAVADMLVSMSNAWETITIYLINNDTFVRHIDNVFDSMICISVVASMCSLLAIAVDRYITIFYALYHHIMTARRSGVIIACIWTFCISCGIVFIIYYY----------------------EESKYVIVCLISMFFTMLFFMVSLYIHMFLLARNH-RIPRYNSVRQRASMKGAITLTMLLGIFIVCWSPFFLHLILMISCPQNACFMSYFNMYLILIMCNSVIDPLIYALRSQEMRRTFKEIICCHGFRRTCTLLGRY------------
>BONZ_CERAE
------MAEYDHYEDNGFNSFNDSSQEEHQDFLQFSKVFLPCMYLVVFVCGLVGNSLVLVISIFYHKLQSLTDVFLVNLPLADLVFVCTLPFWAYAGIHE--WIFGQVMCKTLLGIYTINFYTSMLILTCITVDRFIVVVKATNQQAKKMTWGKVICLLIWVISLLVSLPQIIYGNVFNLD--KLICGYH-----DEEISTVVLATQMTLGFFLPLLAMIVCYSVIIKTL---------LHAGGFQKHRSLKIIFLVMAVFLLTQTPFNLVKLIRSTHWEYTSFHYTIIVTEAIAYLRACLNPVLYAFVSLKFRKNFWKLVKDIGCLPYLGVSHQWKTFSASHNVEAT
>SSR5_MOUSE
LSLASTPSWNAS---AASSGSHNWSLVDPVSPMGARAVLVPVLYLLVCTVGLGGNTLVIYVVLRYAKMKTVTNVYILNLAVADVLFMLGLPFLATQNAVSY-WPFGSFLCRLVMTLDGINQFTSIFCLMVMSVDRYLAVVHPLSARWRRPRVAKLASAAVWVFSLLMSLPLLVFADVQE-----GTCNLS-WPEPVGLWGAAFITYTSVLGFFGPLLVICLCYLLIVVKVKAAMR--VGSSRRRRSERKVTRMVVVVVLVFVGCWLPFFIVNIVNLAFTLPPTSAGLYFFVVVLSYANSCANPLLYGFLSDNFRQSFRKALCLRRGYGVEDADAIERPQTTLPTRSCE
>OPSR_ASTFA
GD--DTTREAAFTYTNSNNTKDPFEGPNYHIAPRWVYNLATCWMFFVVVASTVTNGLVLVASAKFKKLRHPLNWILVNLAIADLLETLLASTISVCNQFFGYFILGHPMCVFEGFTVATCGIAGLWSLTVISWERWVVVCKPFGNVKFDGKMATAGIVFTWVWSAVWCAPPIFGWSRYWPHGLKTSCGPDVFSGSEDPGVQSYMIVLMITCCFIPLGIIILCYIAVWWAIRTVAQQQKDSESTQKAEKEVSRMVVVMIMAYCFCWGPYTFFACFAAANPGYAFHPLAAAMPAYFAKSATIYNPVIYVFMNRQFRVCIMQLFGKKVDDG-----SEVSS------VAPA
>C.C1_HUMAN
------MESSGNPESTTFFYYDLQSQPCENQAWVFATLATTVLYCLVFLLSLVGNSLVLWVLVKYESLESLTNIFILNLCLSDLVFACLLPVWISPYHWG--WVLGDFLCKLLNMIFSISLYSSIFFLTIMTIHRYLSVVSPLTLRVPTLRCRVLVTMAVWVASILSSILDTIFHKVLS-----SGCDYS------ELTWYLTSVYQHNLFFLLSLGIILFCYVEILRTL---------FRSRSKRRHRTVKLIFAIVVAYFLSWGPYNFTLFLQTLFRTQQQLEYALLICRNLAFSHCCFNPVLYVFVGVKFRTHLKHVLRQFWFCRLQAPSPASFAYEGASFY---
>GALT_RAT
--------------------MADIQNISLDSPGSVGAVAVPVIFALIFLLGMVGNGLVLAVLLQPGPPRSTTDLFILNLAVADLCFILCCVPFQAAIYTLDAWLFGAFVCKTVHLLIYLTMYASSFTLAAVSLDRYLAVRHPLSRALRTPRNARAAVGLVWLLAALFSAPYLSYYG---GA--LELCVPA----WEDARRRALDVATFAAGYLLPVAVVSLAYGRTLCFLWAAP--AAAAEARRRATGRAGRAMLAVAALYALCWGPHHALILCFWYGRFAPATYACRLASHCLAYANSCLNPLVYSLASRHFRARFRRLWPCGRRRHRHHHR-AHPASSGPAGYPGD
>ACM4_HUMAN
-----MANFTPVNGSSGNQSVRLVTSSSHNRYETVEMVFIATVTGSLSLVTVVGNILVMLSIKVNRQLQTVNNYFLFSLACADLIIGAFSMNLYTVYIIKGYWPLGAVVCDLWLALDYVVSNASVMNLLIISFDRYFCVTKPLYPARRTTKMAGLMIAAAWVLSFVLWAPAILFWQFVV-VP-DNHCFIQ------FLSNPAVTFGTAIAAFYLPVVIMTVLYIHISLASRSRSIAVRKKRQMAARERKVTRTIFAILLAFILTWTPYNVMVLVNTFCQSC-IPDTVWSIGYWLCYVNSTINPACYALCNATFKKTFRHLLLCQYRNIGTAR----------------
>AG22_RAT
RNITSSLPFDNLNATGTNESAFNC----SHKPADKHLEAIPVLYYMIFVIGFAVNIVVVSLFCCQKGPKKVSSIYIFNLAVADLLLLATLPLWATYYSYRYDWLFGPVMCKVFGSFLTLNMFASIFFITCMSVDRYQSVIYPFLSQRRNPWQASYVVPLVWCMACLSSLPTFYFRDVRT-LG--NACIMAFPPEKYAQWSAGIALMKNILGFIIPLIFIATCYFGIRKHLLKT----NSYGKNRITRDQVLKMAAAVVLAFIICWLPFHVLTFLDALTWMGAVIDLALPFAILLGFTNSCVNPFLYCFVGNRFQQKLRSVFRVPITWLQGKRETMSREMD-------T
>OPRM_MOUSE
LSHVDGNQSDPCGPNRTGLGGSHSLCPQTGSPSMVTAITIMALYSIVCVVGLFGNFLVMYVIVRYTKMKTATNIYIFNLALADALATSTLPFQSVNYLMGT-WPFGNILCKIVISIDYYNMFTSIFTLCTMSVDRYIAVCHPVALDFRTPRNAKIVNVCNWILSSAIGLPVMFMATTKYGS---IDCTLT-FSHPTWYWENLLKICVFIFAFIMPVLIITVCYGLMILRLKSVRML-SGSKEKDRNLRRITRMVLVVVAVFIVCWTPIHIYVIIKALITIPTFQTVSWHFCIALGYTNSCLNPVLYAFLDENFKRCFREFCIPTSSTIEQQNSARIPSTANTVDRTNH
>MC4R_HUMAN
GMHTSLHLWNRSSYRLHSNASESLGKGYSDGGCYEQLFVSPEVFVTLGVISLLENILVIVAIAKNKNLHSPMYFFICSLAVADMLVSVSNGSETIIITLLNSQSFTVNIDNVIDSVICSSLLASICSLLSIAVDRYFTIFYALYHNIMTVKRVGIIISCIWAACTVSGILFIIYYS----------------------DDSSAVIICLITMFFTMLALMASLYVHMFLMARLH-RIPGTGAIRQGANMKGAITLTILIGVFVVCWAPFFLHLIFYISCPQNVCFMSHFNLYLILIMCNSIIDPLIYALRSQELRKTFKEIICCYPLGGLCDLSSRY------------
>5H7_CAVPO
STWTPRLLSGVPEVAASPSPSNVSGCGEQINYGRAEKVVIGSILTLITLLTIAGNCLVVISVCFVKKLRQPSNYLIVSLALADLSVAVAVIPFVSVTDLIGGWIFGHFFCNVFIAMDVMCCTASIMTLCVISIDRYLGITRPLYPVRQNGKCMPKMILSVWLLSASITLPPLFGWAQ--N-D--KVCLIS--------QDFGYTIYSTAVAFYIPMSVMLFMYYRIYKAARKSSRLERKNISIFKREQKAATTLGIIVGAFTVCWLPFFLLSTARPFICGTCIPLWVERTCLWLGYANSLINPFIYAFFNRDLRTTYRSLLQCQYRNINRKLSAAGAERPERPECVLQ
>NK2R_HUMAN
----MGTCDIVTEANISSGPESNTTGITAFSMPSWQLALWAPAYLALVLVAVTGNAIVIWIILAHRRMRTVTNYFIVNLALADLCMAAFNAAFNFVYASHNIWYFGRAFCYFQNLFPITAMFVSIYSMTAIAADRYMAIVHPF-QPRLSAPSTKAVIAGIWLVALALASPQCFYST---GA---TKCVVAWPEDSGGKTLLLYHLVVIALIYFLPLAVMFVAYSVIGLTLWRRPGHHGANLRHLQAKKKFVKTMVLVVLTFAICWLPYHLYFILGSFQEDIKFIQQVYLALFWLAMSSTMYNPIIYCCLNHRFRSGFRLAFRCCPWVTPTKE----PTTSLSTRVNRC
>NY1R_CAVPO
TSFSQLENHSVHYNLSEEKPSFFAFENDDCHLPLAVIFTLALAYGAVIILGVSGNLALILIILKQKEMRNVTNILIVNLSFSDLLVAIMCLPFTFVYTLMDHWIFGEIMCKLNPFVQCVSITVSIFSLVLIAVERHQLIINPR-GWRPNNRHAYIGIAVIWVLAVASSLPFMIYQVLTDKD--KLVCFDQ---FPSDSHRLSYTTLLLVLQYFGPLCFIFICYFKIYIRLKRRNMMMRDSKYRSSESKRINIMLLSIVVAFAVCWLPLTIFNTVFDWNHQICNHNLLFLLCHLTAMISTCVNPIFYGFLNKNFQRDLQFFFNFCDFRSRDDDY---TMHTDVSKTSLK
>O.1R_RAT
PGVPTSSGEPFHLPPDYED-EFLRYLWRDYLYPKQYEWVLIAAYVAVFLIALVGNTLVCLAVWRNHHMRTVTNYFIVNLSLADVLVTAICLPASLLVDITESWLFGHALCKVIPYLQAVSVSVAVLTLSFIALDRWYAICHPL-LFKSTARRARGSILGIWAVSLAVMVPQAAVMECSSRTRLFSVCDER---WADELYPKIYHSCFFFVTYLAPLGLMGMAYFQIFRKLWGPQPRFLAEVKQMRARRKTAKMLMVVLLVFALCYLPISVLNVLKRVFGMFEAVYACFTFSHWLVYANSAANPIIYNFLSGKFREQFKAAFSCCLPGLGPS-----RHKS---LSLQS
>YTJ5_CAEEL
-------------MPNYTVPPDPADTSWDSPYSIPVQIVVWIIIIVLSLETIIGNAMVVMAYRIERNSKQVSNRYIVSLAISDLIIGIEGFPFFTVYVLNGDWPLGWVACQTWLFLDYTLCLVSILTVLLITADRYLSVCHTAYLKWQSPTKTQLLIVMSWLLPAIIFGIMIYGWQAMTQSTSGAECSAP------FLSNPYVNMGMYVAYYWTTLVAMLILYKVFSSGYQKKSQPDRLAPPNKTDTFLSASGTITFIVGFFAILWSPYYIMATVYGFCKG--IPSFLYTLSYYMCYLNSSGNPFAYALANRQFRSAFMRMFRGNFNKVA------------------
>KI01_HUMAN
---------------MINSTSTQPPDESCSQNLLITQQIIPVLYCMVFIAGILLNGVSGWIFFYVPS-SKSFIIYLKNIVIADFVMSLTFPFKILGDSGLGPWQLNVFVCRVSAVLFYVNMYVSIVFFGLISFDRYYKIVKPLTSFIQSVSYSKLLSVIVWMLMLLLAVPNIILTNQS-TQ---IKCIEL--KSELGRKWHKASNYIFVAIFWIVFLLLIVFYTAITKKIFK----LKSSRNSTSVKKKSSRNIFSIVFVFFVCFVPYHIARIPYTKSQTEEILRYMKEFTLLLSAANVCLDPIIYFFLCQPFREILCKKLHIPLKAQNDLDISRIESTDTL------
>PD2R_HUMAN
---------------------MKSPFYRCQNTTSVEKGNSAVMGGVLFSTGLLGNLLALGLLARSGLLPSVFYMLVCGLTVTDLLGKCLLSPVVLAAYAQNRPALDNSLCQAFAFFMSFFGLSSTLQLLAMALECWLSLGHPFYRRHITLRLGALVAPVVSAFSLAFCALPFMGFGKFVQYCPGTWCFIQ-MVHEEGSLSVLGYSVLYSSLMALLVLATVLCNLGAMRNLYAMAEPGREASPQPLEELDHLLLLALMTVLFTMCSLPVIYRAYYGAFKDV-TSEEAEDLRALRFLSVISIVDPWIFIIFRSPVFRIFFHKIFIRPLRYRSRCSNST------------
>C3.1_RAT
--MPTSFPELDLENFEYDDSAEAC---YLGDIVAFGTIFLSIFYSLVFTFGLVGNLLVVLALTNSRKSKSITDIYLLNLALSDLLFVATLPFWTHYLISHEG--LHNAMCKLTTAFFFIGFFGGIFFITVISIDRYLAIVLAASMNNRTVQHGVTISLGVWAAAILVASPQFMFTKRKD-----NECLGDYPEVLQEIWPVLRNSEVNILGFVLPLLIMSFCYFRIVRTL---------FSCKNRKKARAIRLILLVVVVFFLFWTPYNIVIFLETLKFYNRDLRWALSVTETVAFSHCCLNPFIYAFAGEKFRRYLRHLYNKCLAVLCGRPVHAGRSRQDSILSS-L
>UR2R_RAT
TVSGSTVTELPGDSNVSLNSSWSGPTDPSSLKDLVATGVIGAVLSAMGVVGMVGNVYTLVVMCRFLRASASMYVYVVNLALADLLYLLSIPFIIATYVTKD-WHFGDVGCRVLFSLDFLTMHASIFTLTIMSSERYAAVLRPLDTVQRSKGYRKLLVLGTWLLALLLTLPMMLAIQ---GS--KSLCLPA----WGPRAHRTYLTLLFGTSIVGPGLVIGLLYVRLARAYWLS---ASFKQTRRLPNPRVLYLILGIVLLFWACFLPFWLWQLLAQYHEAMETARIVNYLTTCLTYGNSCINPFLYTLLTKNYREYLRGRQRSLGSSCHSPGSPGSLQQDSGRSLSSS
>5H4_MOUSE
------------------MDKLDANVSSNEGFRSVEKVVLLTFLAVVILMAILGNLLVMVAVCRDRQRKIKTNYFIVSLAFADLLVSVLVMPFGAIELVQDIWAYGEMFCLVRTSLDVLLTTASIFHLCCISLDRYYAICCQPYRNKMTPLRIALMLGGCWVLPMFISFLPIMQGWNNIRK-NSTWCVFM--------VNKPYAITCSVVAFYIPFLLMVLAYYRIYVTAKEHSRPDQHSTHRMRTETKAAKTLCVIMGCFCFCWAPFFVTNIVDPFIDYT-VPEQVWTAFLWLGYINSGLNPFLYAFLNKSFRRAFLIILCCDDERYKRPPILGQTINGSTHVLR--
>MC4R_RAT
GMYTSLHLWNRSSHGLHGNASESLGKGHSDGGCYEQLFVSPEVFVTLGVISLLENILVIVAIAKNKNLHSPMYFFICSLAVADMLVSVSNGSETIVITLLNSQSFTVNIDNVIDSVICSSLLASICSLLSIAVDRYFTIFYALYHNIMTVRRVGIIISCIWAACTVSGVLFIIYYS----------------------DDSSAVIICLITMFFTMLVLMASLYVHMFLMARLH-RIPGTGTIRQGANMKGAITLTILIGVFVVCWAPFFLHLLFYISCPQNVCFMSHFNLYLILIMCNAVIDPLIYALRSQELRKTFKEIICFYPLGGICELPGRY------------
>OLF2_CANFA
-------------MDGKNCSSVNEFLLVGISNKPGVKVTLFITFLIVYLIILVANLGMIILIRMDSQLHTPMYFFLSHLSFSDARYSTAVGPRMLVGFIAKNKSIPFYSCAMQWLVFCTFVDSECLLLAVMAFDRYKAISHPLYTVSMSSRVCSLLMAGVYLVGIMDASVNTILTFRL---CNHFFCDVPLLSCSDTQVNELVIFTIFGFIELITLSGLFVSYCYIILAV--------RKINSAEGRFKAFSTCTSHLTAVAIFQGTMLFMYFRPSSSY---SLDQDKIISLFYSLVIPMLNPLIYSLRNKDVKEALKKLKNKKWFH---------------------
>AG2S_.ENLA
----MLSNISAGENSEVEKIVVKC---SKSGMHNYIFITIPIIYSTIFVVGVFGNSLVVIVIYSYMKMKTMASVFLMNLALSDLCFVITLPLWAVYTAMHYHWPFGDLLCKIASTAITLNLYTTVFLLTCLSIDRYSAIVHPMSRIRRTVMVARLTCVGIWLVAFLASLPSVIYRQIFI-TN--TVCALV-YHSGHIYFMVGMSLVKNIVGFFIPFVIILTSYTLIGKTLKEV------YRAQRARNDDIFKMIVAVVLLFFFCWIPHQVFTFLDVLIQMDDIVDTGMPITICIAYFNSCLNPFLYGFFGKKFRKHFLQLIKYIPPKMRTHASVNTRLSD-------T
>OPRD_MOUSE
SSPLVNLSDAFPSAFPSAGANASGSPGARSASSLALAIAITALYSAVCAVGLLGNVLVMFGIVRYTKLKTATNIYIFNLALADALATSTLPFQSAKYLMET-WPFGELLCKAVLSIDYYNMFTSIFTLTMMSVDRYIAVCHPVALDFRTPAKAKLINICIWVLASGVGVPIMVMAVTQPGA---VVCMLQ-FPSPSWYWDTVTKICVFLFAFVVPILIITVCYGLMLLRLRSVRLL-SGSKEKDRSLRRITRMVLVVVGAFVVCWAPIHIFVIVWTLVDINPLVVAALHLCIALGYANSSLNPVLYAFLDENFKRCFRQLCRTPCGRQEPGSLRRPRVTACTPSDGPG
>CKR5_PYGNE
----MDYQVSSPTYDIDYYTSEPC---QKVNVKQIAARLLPPLYSLVFIFGFVGNILVVLILINCKRLKSMTDIYLLNLAISDLFFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWVVAVFASLPGIIFTRSQR-EG-HYTCSSHFPYSQYQFWKNFQTLKIVILGLVLPLLIMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCSIFA-------SSVY
>OPSD_GAMAF
MNGTEGPYFYVPMVNTTGIVRSPYEYPQYYLVSPAAYACLGAYMFFLILVGFPVNFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTIYTSMHGYFVLGRLGCNLEGYFATLGGEIGLWSLVVLAVERWLVVCKPISNFRFTENHAIMGLVFTWIMANACAAPPLLGWSRYIPEGMQCSCGVDYYTRAEGFNNESFVIYMFICHFCIPLVVVFFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVVILVIGFLVCWTPYASVAWYIFSNQGSEFGPLFMTIPAFFAKSSSIYNPMIYICMNKQFRHCMITTLCCGKNPFEEEEGASTASSVSSSSVSPA
>LSHR_CALJA
ESGQSGWDYDYGFHLPKTPRCAPEPDAFNPCEDIMGYDFLRVLIWLINILAIMGNMTVLFVLLTSRYKLTVPRFLMCNLSFADFCMGLYLLLIASVDSQTKGWQTG-SGCNTAGFFTVFASELSVYTLTVITLERWHTITYAILDQKLRLRHAILIMLGGWLFSSLIAMLPLVGVSNY----KVSICFPM-----VETTLSQIYILTILILNVVAFIIICACYIKIYFAVRNP------ELMATNKDTKIAKKMAILIFTDFTCMAPISFFAISAAFKMPLITVTNSKVLLVLFYPINSCANPFLYAIFTKTFRRDFFLLLGKFGCCKHRAELYRRSNYKNGFTGSSK
>O1F1_HUMAN
-------------MSGTNQSSVSEFLLLGLSRQPQQQHLLFVFFLSMYLATVLGNLLIILSVSIDSCLHTPMYFFLSNLSFVDICFSFTTVPKMLANHILETQTISFCGCLTQMYFVFMFVDMDNFLLAVMAYDHFVAVCHPLYTAKMTHQLCALLVAGLWVVANLNVLLHTLLMAPL---STHFFCDVTKLSCSDTHLNEVIILSEGALVMITPFLCILASYMHITCTV--------LKVPSTKGRWKAFSTCGSHLAVVLLFYSTIIAVYFNPLSSH---SAEKDTMATVLYTVVTPMLNPFIYSLRNRYLKGALKKVVGRVVFSV--------------------
>GALT_HUMAN
--------------------MADAQNISLDSPGSVGAVAVPVVFALIFLLGTVGNGLVLAVLLQPGPPGSTTDLFILNLAVADLCFILCCVPFQATIYTLDAWLFGALVCKAVHLLIYLTMYASSFTLAAVSVDRYLAVRHPLSRALRTPRNARAAVGLVWLLAALFSAPYLSYYG---GA--LELCVPA----WEDARRRALDVATFAAGYLLPVAVVSLAYGRTLRFLWAAP--AAAAEARRRATGRAGRAMLAVAALYALCWGPHHALILCFWYGRFAPATYACRLASHCLAYANSCLNPLVYALASRHFRARFRRLWPCGRRRRHRARRALRGPPGCPGDARPS
>A2AA_CAVPO
----MGSLQPDSGNASWNGTEGPGGGTRATPYSLQVTVTLVCLVGLLILLTVFGNVLVIIAVFTSRALKAPQNLFLVSLASADILVATLVIPFSLANEVMGYWYFGKAWCEIYLALDVLFCTSSIVHLCAISLDRYWSITQAIYNLKRTPRRIKAIIVTVWVISAVISFPPLISFEK-AQ-P--PRCEIN--------DQKWYVISSSIGSFFAPCLIMILVYVRIYQIAKRRRGGASRWRGRQNREKRFTFVLAVVIGVFVVCWFPFFFTYTLTAVGCS--VPRTLFKFFFWFGYCNSSLNPVIYTIFNHDFRRAFKKILCRGDRKRIV------------------
>A2AA_RAT
----MGSLQPDAGNSSWNGTEAPGGGTRATPYSLQVTLTLVCLAGLLMLFTVFGNVLVIIAVFTSRALKAPQNLFLVSLASADILVATLVIPFSLANEVMGYWYFGKVWCEIYLALDVLFCTSSIVHLCAISLDRYWSITQAIYNLKRTPRRIKAIIVTVWVISAVISFPPLISIEKKGQ-P--PSCKIN--------DQKWYVISSSIGSFFAPCLIMILVYVRIYQIAKRRRAGASRWRGRQNREKRFTFVLAVVIGVFVVCWFPFFFTYTLIAVGCP--VPYQLFNFFFWFGYCNSSLNPVIYTIFNHDFRRAFKKILCRGDRKRIV------------------
>OLF4_CANFA
-------------MELENDTRIPEFLLLGFSEEPKLQPFLFGLFLSMYLVTILGNLLLILAVSSDSHLHTPMYFFLANLSFVDICFTCTTIPKMLVNIQTQRKVITYESCIIQMYFFELFAGIDNFLLTVMAYDRYMAICYPLYMVIMNPQLCSLLLLVSWIMSALHSLLQTLMVLRL---SPHFFCELNQLACSDTFLNNMMLYFAAILLGVAPLVGVLYSYFKIVSSI--------RGISSAHSKYKAFSTCASHLSVVSLFYCTSLGVYLSSAAPQ---STHTSSVASVMYTVVTPMLNPFIYSLRNKDIKGALNVFFRGKP-----------------------
>SSR2_RAT
QFNGSQVWIPSPFDLNGSLGPSNGSNQTEPYYDMTSNAVLTFIYFVVCVVGLCGNTLVIYVILRYAKMKTITNIYILNLAIADELFMLGLPFLAMQVALVH-WPFGKAICRVVMTVDGINQFTSIFCLTVMSIDRYLAVVHPISAKWRRPRTAKMINVAVWGVSLLVILPIMIYAGLRSWG--RSSCTIN-WPGESGAWYTGFIIYAFILGFLVPLTIICLCYLFIIIKVKSSGIR-VGSSKRKKSEKKVTRMVSIVVAVFIFCWLPFYIFNVSSVSVAISPALKGMFDFVVILTYANSCANPILYAFLSDNFKKSFQNVLCLVKVSGAEDGERSDLNETTETQRTLL
>BRS3_SHEEP
QTLISTTNDTESSSSVVPNDSTNKRRTGDNSPGIEALCAIYITYAVIISVGILGNAILIKVFFKTKSMQTVPNIFITSLAFGDLLLLLTCVPVDVTHYLAEGWLFGRIGCKVLSFIRLTSVGVSVFTLTILSADRYKAVVKPLRQPPNAILKTCAKAGCIWIMSMIIALPEAIFSNVYTVT--FKACASY---VSERLLQEIHSLLCFLVFYIIPLSIISVYYSLIARTLYKSIPTQRHARKQIESRKRIAKTVLVLVALFALCWLPNHLLYLYRSFTSQTTVHLFVTIISRILAFSNSCVNPFALYWLSNTFQQHFKAQLFCCKAGRPDPTAANTMGRVPGAASTQM
>ETBR_BOVIN
SSATPQIPRGGRMAGIPPR--TPPPCDGPIEIKETFKYINTVVSCLVFVLGIIGNSTLLRIIYKNKCMRNGPNILIASLALGDLLHIIIDIPINTYKLLAKDWPFGVEMCKLVPFIQKASVGITVLSLCALSIDRYRAVASWSIKGIGVPKWTAVEIVLIWVVSVVLAVPEAVGFDIITRI--LRICLLHQKTAFMQFYKTAKDWWLFSFYFCLPLAITALFYTLMTCEMLRKSGM-IALNDHLKQRREVAKTVFCLVLVFALCWLPLHLSRILKLTLYDQSFLLVLDYIGINMASLNSCINPIALYLVSKRFKNCFKSCLCCWCQSFE-EKQSLEFKANDHGYDNFR
>FML1_HUMAN
-----------METNFSTPLNEYEEVSYESAGYTVLRILPLVVLGVTFVLGVLGNGLVIWVAGFRMT-RTVTTICYLNLALADFSFTATLPFLIVSMAMGEKWPFGWFLCKLIHIVVDINLFGSVFLIGFIALDRCICVLHPVAQNHRTVSLAMKVIVGPWILALVLTLPVFLFLTTVTASWPEERLKVA------ITMLTARGIIRFVIGFSLPMSIVAICYGLIAAKI---------HKKGMIKSSRPLRVLTAVVASFFICWFPFQLVALLGTVWLKEKIIDILVNPTSSLAFFNSCLNPMLYVFVGQDFRERLIHSLPTSLERALSEDSAPTAS----------
>V1AR_MOUSE
SSPWWPLTTEGANSSREAAGLGEGGSPPGDVRNEELAKLEVTVLAVIFVVAVLGNSSVLLALHRTPRKTSRMHLFIRHLSLADLAVAFFQVLPQLCWDITYRFRGPDWLCRVVKHLQVFAMFASSYMLVVMTADRYIAVCHPLKTLQQPARRSRLMIAASWGLSFVLSIPQYFIFSVIETK--AQDCWAT---FIPPWGTRAYVTWMTSGVFVVPVIILGTCYGFICYHIWRNLLVVSSVKSISRAKIRTVKMTFVIVSAYILCWTPFFIVQMWSVWDTNFDSENPSTTITALLASLNSCCNPWIYMFFSGHLLQDCVQSFPCCQSIAQKFAKDDSTSYSNNRSPTNS
>OPSB_BOVIN
--MSKMSEEEEFLLFKNISLVGPWDGPQYHLAPVWAFHLQAVFMGFVFFVGTPLNATVLVATLRYRKLRQPLNYILVNVSLGGFIYCIFSVFIVFITSCYGYFVFGRHVCALEAFLGCTAGLVTGWSLAFLAFERYIIICKPFGNFRFSSKHALMVVVATWTIGIGVSIPPFFGWSRFVPEGLQCSCGPDWYTVGTKYYSEYYTWFLFIFCYIVPLSLICFSYSQLLGALRAVAAQQQESASTQKAEREVSHMVVVMVGSFCLCYTPYAALAMYIVNNRNHGVDLRLVTIPAFFSKSACVYNPIIYCFMNKQFRACIMEMVCGKPMTD---ESELSSTVSSSQVGPN-
>APJ_MOUSE
-----------MEDDGYNYYGADNQSECDYADWKPSGALIPAIYMLVFLLGTTGNGLVLWTVFRTSRKRRSADIFIASLAVADLTFVVTLPLWATYTYREFDWPFGTFSCKLSSYLIFVNMYASVFCLTGLSFDRYLAIVRPVNARLRLRVSGAVATAVLWVLAALLAVPVMVFRSTDAQCY-MDYSMVA-TSNSEWAWEVGLGVSSTAVGFVVPFTIMLTCYFFIAQTIAGHFR--KERIEGLRKRRRLLSIIVVLVVTFALCWMPYHLVKTLYMLGSLLIFLMNVFPYCTCISYVNSCLNPFLYAFFDPRFRQACTSMLCCDQSGCKGTPHSSSSSGHSQGPGPNM
>O1G1_HUMAN
-------------MEGKNLTSISECFLLGFSEQLEEQKPLFGSFLFMYLVTVAGNLLIILVIITDTQLHTPMYFFLANLSLADACFVSTTVPKMLANIQIQSQAISYSGCLLQLYFFMLFVMLEAFLLAVMAYDCYVAICHPLYILIMSPGLCIFLVSASWIMNALHSLLHTLLMNSL---SPHFFCDINSLSCTDPFTNELVIFITGGLTGLICVLCLIISYTNVFSTI--------LKIPSAQGKRKAFSTCSSHLSVVSLFFGTSFCVDFSSPSTH---SAQKDTVASVMYTVVTPMLNPFIYSLRNQEIKSSLRKLIWVRKIHSP-------------------
>CKR2_HUMAN
FIRNTNESGEEVTTFFDYDYGAPC---HKFDVKQIGAQLLPPLYSLVFIFGFVGNMLVVLILINCKKLKCLTDIYLLNLAISDLLFLITLPLWAHSAANE--WVFGNAMCKLFTGLYHIGYFGGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWLVAVFASVPGIIFTKCQK-ED-VYVCGPY----FPRGWNNFHTIMRNILGLVLPLLIMVICYSGILKTL--------LRCRNEKKRHRAVRVIFTIMIVYFLFWTPYNIVILLNTFQEFFSQLDQATQVTETLGMTHCCINPIIYAFVGEKFRSLFHIALGCRIAPLQKPVCGGPV-------KVTT
>5H2A_CRIGR
NSSDASNWTIDGENRTNLSFEGYLPPTCLSILHLQEKNWSALLTAVVIILTIAGNILVIMAVSLEKKLQNATNYFLMSLAIADMLLGFLVMPVSMLTILYGYWPLPSKLCAVWIYLDVLFSTASIMHLCAISLDRYVAIQNPIHSRFNSRTKAFLKIIAVWTISVGVSMPIPVFGLQD-VF---GSCLL---------ADDNFVLIGSFVAFFIPLTIMVITYFLTIKSLQKEEPGGRRTMQSISNEQKACKVLGIVFFLFVVMWCPFFITNIMAVICKESHVIGALLNVFVWIGYLSSAVNPLVYTLFNKTYRSAFSRYIQCQYKENRKPLQLILAYKSSQLQAGQN
>TRFR_MOUSE
------------MENDTVSEMNQTELQPQAAVALEYQVVTILLVVIICGLGIVGNIMVVLVVMRTKHMRTPTNCYLVSLAVADLMVLVAAGLPNITDSIYGSWVYGYVGCLCITYLQYLGINASSCSITAFTIERYIAICHPIAQFLCTFSRAKKIIIFVWAFTSIYCMLWFFLLDLN--NA-SCGYKIS------RNYYSPIYLMDFGVFYVVPMILATVLYGFIARILFLNLNLNRCFNSTVSSRKQVTKMLAVVVILFALLWMPYRTLVVVNSFLSSPFQENWFLLFCRICIYLNSAINPVIYNLMSQKFRAAFRKLCNCKQKPTEKAANYSVKESDRFSTELED
>OPSH_CARAU
MNGTEGNNFYVPLSNRTGLVRSPFEYPQYYLAEPWQFKLLAVYMFFLICLGLPINGLTLICTAQHKKLRQPLNFILVNLAVAGAIMVCFGFTVTFYTAINGYFALGPTGCAVEGFMATLGGEVALWSLVVLAIERYIVVCKPMGSFKFSSTHASAGIAFTWVMAMACAAPPLVGWSRYIPEGIQCSCGPDYYTLNPEYNNESYVLYMFICHFILPVTIIFFTYGRLVCTVKAAAAQQQDSASTQKAEREVTKMVILMVLGFLVAWTPYATVAAWIFFNKGAAFSAQFMAIPAFFSKTSALYNPVIYVLLNKQFRSCMLTTLFCGKNPLGDEESSTVSS------VSPA
>A2AB_RAT
--------------------MSGPTMDHQEPYSVQATAAIASAITFLILFTIFGNALVILAVLTSRSLRAPQNLFLVSLAAADILVATLIIPFSLANELLGYWYFWRAWCEVYLALDVLFCTSSIVHLCAISLDRYWAVSRALYNSKRTPCRIKCIILTVWLIAAVISLPPLIYKGD-Q-----PQCELN--------QEAWYILASSIGSFFAPCLIMILVYLRIYVIAKRSGVAWWRRRTQLSREKRFTFVLAVVIGVFVVCWFPFFFSYSLGAICPQHKVPHGLFQFFFWIGYCNSSLNPVIYTVFNQDFRRAFRRILCRPWTQTGW------------------
>TA2R_CERAE
--MWPNG-----------SSLGPCFRPTNITLEERRLIASPWFAASFCVVGLASNLLALSVLAGARQTRSSFLTFLCGLVLTDFLGLLVTGAIVVSQHAALFVDPGCRLCRFMGVVMIFFGLSPLLLGATMASERFLGITRPFRPVVTSQRRAWATVGLVWAAALALGLLPLLGLGRYTVQYPGSWCFLT----LGAESGDVAFGLLFSMLGGLSVGLSFLLNTVSVATLCHVYHGEAAQQRPRDSEVEMMAQLLGIMLVASVCWLPLLVFIAQTVLRNPPRATEQELLIYLRVATWNQILDPWVYILFRRAVLRRLQPRLSTRPRSLSLQPQLTQ------------
>SSR2_HUMAN
PLNGSHTWLSIPFDLNGSVVSTNTSNQTEPYYDLTSNAVLTFIYFVVCIIGLCGNTLVIYVILRYAKMKTITNIYILNLAIADELFMLGLPFLAMQVALVH-WPFGKAICRVVMTVDGINQFTSIFCLTVMSIDRYLAVVHPISAKWRRPRTAKMITMAVWGVSLLVILPIMIYAGLRSWG--RSSCTIN-WPGESGAWYTGFIIYTFILGFLVPLTIICLCYLFIIIKVKSSGIR-VGSSKRKKSEKKVTRMVSIVVAVFIFCWLPFYIFNVSSVSMAISPALKGMFDFVVVLTYANSCANPILYAFLSDNFKKSFQNVLCLVKVSGTDDGERSDLNETTETQRTLL
>ACM4_.ENLA
----MENDTWENESSASNHSIDETIVEIPGKYQTMEMIFIATVTGSLSLVTVVGNILVMLSIKVNRQLQTVNNYFLFSLACADLIIGVFSMNLYSLYIIKGYWPLGPIVCDLWLALDYVVSNASVMNLLIISLER-FCVTKPLYPARRTTKMAGLMIAAAWLLSFELWAPAILFWQFIV-VP-SGECYIQ------FLSNPAVTFGTAIAAFYLPVVIMTILYIHISLASRSRSIAVRKKRQMAAREKKVTRTIFAILLAFIITWTPYNVMVLINTFCQTC-IPETIWYIGYWLCYVNSTINPACYALCNATFKKTFKHLLMCQYKSIGTAR----------------
>Y..5_CAEEL
SVNESCDNYVEIFNKINYFFRDDQVINGTEYSPKEFGYFITFAYMLIILFGAIGNFLTIIVVILNPAMRTTRNFFILNLALSDFFVCIVTAPTTLYTVLYMFWPFSRTLCKIAGSLQGFNIFLSTFSIASIAVDRYVLIIFPT-KRERQQNLSFCFFIMIWVISLILAVPLLQASDLTPCDLALYICHEQEIWEKMIISKGTYTLAVLITQYAFPLFSLVFAYSRIAHRMKLRTTNSQRRRSVVERQRRTHLLLVCVVAVFAVAWLPLNVFHIFNTFELVN-FSVTTFSICHCLAMCSACLNPLIYAFFNHNFRIEFMHLFDRVGLRSLRVVIFGEMRTEFRSRGGCK
>GRHR_CLAGA
TLLLSNPTNVLDNSSVLNVSVSPPVLKWETPTFTTAARFRVAATLVLFVFRAASNLSVLLSVTRGRGLASHLRPLIASLASADLVMTFVVMPLDAVWNVTVQWYAGDAMCKLMCFLKLFAMHSAAFILVVVSLDRHHAILHPL-DTLDAGRRNRRMLLTAWILSLLLASPQLFIFRAIKVD--FVQCATH--SFQQHWQETAYNMFHFVTLYVFPLLVMSLCYTRILVEINRQGEPRSGTDMIPKARMKTLKMTIIIVASFVICWTPYYLLGIWYWFQPQMVIPDYVHHVFFVFGNLNTCCDPVIYGFFTPSFRADLSRCFCWRNQNASAKSLPHFSGEAESDLGSGD
>CKR5_HUMAN
----MDYQVSSPIYDINYYTSEPC---QKINVKQIAARLLPPLYSLVFIFGFVGNMLVILILINCKRLKSMTDIYLLNLAISDLFFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAVVHAVALKARTVTFGVVTSVITWVVAVFASLPGIIFTRSQK-EG-HYTCSSHFPYSQYQFWKNFQTLKIVILGLVLPLLVMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCSIFA-------SSVY
>OPSD_NEOSA
MNGTEGPYFYVPMVNTTGVVRSPYEYPQYYLVNPAAFAVLGAYMFFLIIFGFPINFLTLYVTLEHKKLRTPLNYILLNLAVADLFMVIGGFTTTMYSSMHGYFVLGRLGCNLEGFSATLGGMISLWSLAVLAIERWVVVCKPTSNFRFGENHAIMGVSLTWTMALACTVPPLVGWSRYIPEGMQCSCGIDYYTRAEGFNNESFVLYMFFCHFMVPLIIIFFCYGRLLCAVKEAAAAQQESETTQRAEREVTRMVILMVIGYLVCWLPYASVAWFIFTHQGSEFGPLFMTIPAFFAKSSSIYNPVIYICMNKQFRNCMITTLFCGKNPF---EGEEETEASSASSVSPA
>A2AB_HUMAN
-------------------------MDHQDPYSVQATAAIAAAITFLILFTIFGNALVILAVLTSRSLRAPQNLFLVSLAAADILVATLIIPFSLANELLGYWYFRRTWCEVYLALDVLFCTSSIVHLCAISLDRYWAVSRALYNSKRTPRRIKCIILTVWLIAAVISLPPLIYKGD-Q-----PQCKLN--------QEAWYILASSIGSFFAPCLIMILVYLRIYLIAKRSGAIWWRRRAHVTREKRFTFVLAVVIGVFVLCWFPFFFSYSLGAICPKHKVPHGLFQFFFWIGYCNSSLNPVIYTIFNQDFRRAFRRILCRPWTQTAW------------------
>CKR4_MOUSE
VTDTTQDETVYNSYYFYESMPKPC---TKEGIKAFGEVFLPPLYSLVFLLGLFGNSVVVLVLFKYKRLKSMTDVYLLNLAISDLLFVLSLPFWGYYAADQ--WVFGLGLCKIVSWMYLVGFYSGIFFIMLMSIDRYLAIVHAVSLKARTLTYGVITSLITWSVAVFASLPGLLFSTCYT-EH-HTYCKTQ-YSVNSTTWKVLSSLEINVLGLLIPLGIMLFWYSMIIRTL---------QHCKNEKKNRAVRMIFGVVVLFLGFWTPYNVVLFLETLVELERYLDYAIQATETLGFIHCCLNPVIYFFLGEKFRKYITQLFRTCRGPLVLCKHCDFMS------SSSY
>A1AD_HUMAN
GSGEDNRSSAGEPGSAGAGGDVNGTAAVGGLVVSAQGVGVGVFLAAFILMAVAGNLLVILSVACNRHLQTVTNYFIVNLAVADLLLSATVLPFSATMEVLGFWAFGRAFCDVWAAVDVLCCTASILSLCTISVDRYVGVRHSLYPAIMTERKAAAILALLWVVALVVSVGPLLGWKEP--VP--RFCGIT--------EEAGYAVFSSVCSFYLPMAVIVVMYCRVYVVARSTHTFLSVRLLKFSREKKAAKTLAIVVGVFVLCWFPFFFVLPLGSLFPQLKPSEGVFKVIFWLGYFNSCVNPLIYPCSSREFKRAFLRLLRCQCRRRRRRRPLWRASTSG-LRQDCA
>C3AR_HUMAN
--------------MASFSAETNSTDLLSQPWNEPPVILSMVILSLTFLLGLPGNGLVLWVAGLKMQ-RTVNTIWFLHLTLADLLCCLSLPFSLAHLALQGQWPYGRFLCKLIPSIIVLNMFASVFLLTAISLDRCLVVFKPICQNHRNVGMACSICGCIWVVAFVMCIPVFVYREIFT-ED-YNLGQFT-DDDQVPTPLVAITITRLVVGFLLPSVIMIACYSFIVFRM--------QRGRFAKSQSKTFRVAVVVVAVFLVCWTPYHIFGVLSLLTDPEKTLMSWDHVCIALASANSCFNPFLYALLGKDFRKKARQSIQGILEAAFSEELTRSVIS---------
>AG2R_HUMAN
------MILNSSTEDGIKRIQDDC---PKAGRHNYIFVMIPTLYSIIFVVGIFGNSLVVIVIYFYMKLKTVASVFLLNLALADLCFLLTLPLWAVYTAMEYRWPFGNYLCKIASASVSFNLYASVFLLTCLSIDRYLAIVHPMSRLRRTMLVAKVTCIIIWLLAGLASLPAIIHRNVFF-TN--TVCAFH-YESQNSTLPIGLGLTKNILGFLFPFLIILTSYTLIWKALKKA----YEIQKNKPRNDDIFKIIMAIVLFFFFSWIPHQIFTFLDVLIQLGDIVDTAMPITICIAYFNNCLNPLFYGFLGKKFKRYFLQLLKYIPPKAKSHSNLSTRPSD-------N
>CKR5_PANTR
----MDYQVSSPIYDIDYYTSEPC---QKINVKQIAARLLPPLYSLVFIFGFVGNMLVILILINCKRLKSMTDIYLLNLAISDLFFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWVVAVFASLPGIIFTRSQK-EG-HYTCSSHFPYSQYQFWKNFQTLKIVILGLVLPLLVMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCSIFA-------SSVY
>PF2R_MOUSE
--MSMNS---------SKQPVSPAAGLIANTTCQTENRLSVFFSIIFMTVGILSNSLAIAILMKAYQSKASFLLLASGLVITDFFGHLINGGIAVFVYASDKFDQSNILCSIFGISMVFSGLCPLFLGSAMAIERCIGVTNPIHSTKITSKHVKMILSGVCMFAVFVAVLPILGHRDYQIQASRTWCFYN--TEHIEDWEDRFYLLFFSFLGLLALGVSFSCNAVTGVTLLRVKFRSQQHRQGRSHHLEMIIQLLAIMCVSCVCWSPFLVTMANIAINGNNPVTCETTLFALRMATWNQILDPWVYILLRKAVLRNLYKLASRCCGVNIISLHIWELKVAAISESPAA
>OPSB_ORYLA
VEFPDDFWIPIPLDTNNVTALSPFLVPQDHLGSPTIFYSMSALMFVLFVAGTAINLLTIACTLQYKKLRSHLNYILVNMAVANLIVASTGSSTCFVCFAFKYMVLGPLGCKIEGFTAALGGMVSLWSLAVIAFERWLVICKPLGNFVFKSEHALLCCALTWVCGLCASVPPLVGWSRYIPEGMQCSCGPDWYTTGNKFNNESFVMFLFCFCFAVPFSIIVFCYSQLLFTLKMAAKAQADSASTQKAEKEVTRMVVVMVVAFLVCYVPYASFALWVINNRGQTFDLRLATIPSCVSKASTVYNPVIYVLLNKQFRLCMKKMLGMSADED---EESSTSKVGPS------
>CCKR_RABIT
ASLLGNASGIPPPCELGLDNETLFCLDQPPPSKEWQPAVQILLYSLIFLLSVLGNTLVITVLIRNKRMRTVTNIFLLSLAISDLMLCLFCMPFNLIPNLLKDFIFGSALCKTTTYLMGTSVSVSTLNLVAISLERYGAICKPLSRVWQTKSHALKVIAATWCLSFAIMTPYPIYN----NNQTANMCRFL---LPSDVMQQAWHTFLLLILFLIPGIVMMVAYGMISLELYQGRVSSSSSAATLMAKKRVIRMLMVIVVLFFLCWMPIFSANAWRAYDTVSRLSGTPISFILLLSYTSSCVNPIIYCFMNRRFRLGFMATFPCCPNPGPPGPRAEATTRASLSRYSYS
>BRB1_HUMAN
SSWPPLELQSSNQSQLFPQNATAC--DNAPEAWDLLHRVLPTFIISICFFGLLGNLFVLLVFLLPRRQLNVAEIYLANLAASDLVFVLGLPFWAENIWNQFNWPFGALLCRVINGVIKANLFISIFLVVAISQDRYRVLVHPMSGRQQRRRQARVTCVLIWVVGGLLSIPTFLLRSIQA-LN--TACILL---LPHEAWHFARIVELNILGFLLPLAAIVFFNYHILASLRTR---VSRTRVRGPKDSKTTALILTLVVAFLVCWAPYHFFAFLEFLFQVQDFIDLGLQLANFFAFTNSSLNPVIYVFVGRLFRTKVWELYKQCTPK---------LAPI-------S
>TRFR_RAT
------------MENETVSELNQTELPPQVAVALEYQVVTILLVVVICGLGIVGNIMVVLVVMRTKHMRTATNCYLVSLAVADLMVLVAAGLPNITDSIYGSWVYGYVGCLCITYLQYLGINASSCSITAFTIERYIAICHPIAQFLCTFSRAKKIIIFVWAFTSIYCMLWFFLLDLN--DA-SCGYKIS------RNYYSPIYLMDFGVFYVMPMILATVLYGFIARILFLNMNLNRCFNSTVSSRKQVTKMLAVVVILFALLWMPYRTLVVVNSFLSSPFQENWFLLFCRICIYLNSAINPVIYNLMSQKFRAAFRKLCNCKQKPTEKAANYSVKESDRFSTELDD
>AG22_HUMAN
KNITSGLHFGLVNISGNNESTLNC----SQKPSDKHLDAIPILYYIIFVIGFLVNIVVVTLFCCQKGPKKVSSIYIFNLAVADLLLLATLPLWATYYSYRYDWLFGPVMCKVFGSFLTLNMFASIFFITCMSVDRYQSVIYPFLSQRRNPWQASYIVPLVWCMACLSSLPTFYFRDVRT-LG--NACIMAFPPEKYAQWSAGIALMKNILGFIIPLIFIATCYFGIRKHLLKT----NSYGKNRITRDQVLKMAAAVVLAFIICWLPFHVLTFLDALAWMGAVIDLALPFAILLGFTNSCVNPFLYCFVGNRFQQKLRSVFRVPITWLQGKRESMSREME-------T
>VU51_HSV7J
-----------------MKNIDLTNWKLLAEIYEYLFFFSFFFLCLLVIIVVKFNNSTVGR-E--------YTFSTFSGMLVYILLLPVKMGMLTKM-----WDVSTDYCIILMFLSDFSFIFSSWALTLLALERINNFSFSEIKVNETKILKQMSFPIIWVTSIFQAVQISMKYKKSQ---EDDYCLLA--------IRSAEEAWILLMYTVVIPTFIVFFYVLNKRFL-----------FLERDLNSIVTHLSLFLFFGALCFFPASVLNEFNCN----RLFYGLHELLIVCLELKIFYVPTMTYIISCENYRLAAKAFFCKCFKPCFLMPSLRSTQF--------
>PAR3_MOUSE
-WTGATTTIKAECPEDSISTLHVNNATIGYLRSSLSTQVIPAIYILLFVVGVPSNIVTLWKLSLRTK-SISLVIFHTNLAIADLLFCVTLPFKIAYHLNGNNWVFGEVMCRITTVVFYGNMYCAILILTCMGINRYLATAHPFYQKLPKRSFSLLMCGIVWVMVFLYMLPFVILKQEYH--E--TTCHDVDACESPSSFRFYYFVSLAFFGFLIPFVIIIFCYTTLIHKL----------KSKDRIWLGYIKAVLLILVIFTICFAPTNIILVIHHANYYYDSLYFMYLIALCLGSLNSCLDPFLYFVMSKVVDQLNP------------------------------
>GPRC_HUMAN
NLSGLPRDYLDAAAAENISAAVSSRVPAVEPEPELVVNPWDIVLCTSGTLISCENAIVVLIIFHNPSLRAPMFLLIGSLALADLLAG-IGLITNFVFA-Y--LLQSEATKLVTIGLIVASFSASVCSLLAITVDRYLSLYYALYHSERTVTFTYVMLVMLWGTSICLGLLPVMGWNCL-R--DESTCSVV-------RPLTKNNAAILSVSFLFMFALMLQLYIQICKIVMRHIALHFLATSHYVTTRKGVSTLAIILGTFAACWMPFTLYSLIADY----TYPSIYTYATLLPATYNSIINPVIYAFRNQEIQKALCLICCGCIPSSLAQRARSP------------
>BRB2_CAVPO
---MFNITSQVSALNATLAQGNSC---LDAEWWSWLNTIQAPFLWVLFVLAVLENIFVLSVFFLHKSSCTVAEIYLGNLAVADLILAFGLPFWAITIANNFDWLFGEVLCRMVNTMIQMNMYSSICFLMLVSIDRYLALVKTMMGRMRGVRWAKLYSLVIWGCALLLSSPMLVFRTMKD-HN--TACLII---YPSLTWQVFTNVLLNLVGFLLPLSIITFCTVQIMQVLRNN---EMQKFKEIQTERRATVLVLAVLLLFVVCWLPFQIGTFLDTLRLLGHVIDLITQISSYLAYSNSCLNPLVYVIVGKRFRKKSREVYHGLCRSGGCVSEPAQLRTS-------I
>P2YR_RAT
AAFLAGLGSLWGNSTIASTAAVSSSFRCALIKTGFQFYYLPAVYILVFIIGFLGNSVAIWMFVFHMKPWSGISVYMFNLALADFLYVLTLPALIFYYFNKTDWIFGDVMCKLQRFIFHVNLYGSILFLTCISAHRYSGVVYPLSLGRLKKKNAIYVSVLVWLIVVVAISPILFYSGTG----KTVTCYDS-TSDEYLRSYFIYSMCTTVAMFCIPLVLILGCYGLIVRALIY------KDLDNSPLRRKSIYLVIIVLTVFAVSYIPFHVMKTMNLRARLDDRVYATYQVTRGLASLNSCVDPILYFLAGDTFRRRLSRATRKASRRSEANLQSKSLSEFKQNGDTSL
>PF2R_SHEEP
--MSTNN---------SVQPVSPASELLSNTTCQLEEDLSISFSIIFMTVGILSNSLAIAILMKAYQYKSSFLLLASALVITDFFGHLINGTIAVFVYASDKFDKSNILCSIFGICMVFSGLCPLFLGSLMAIERCIGVTKPIHSTKITTKHVKMMLSGVCFFAVFVALLPILGHRDYKIQASRTWCFYK--TDQIKDWEDRFYLLLFAFLGLLALGISFVCNAITGISLLKVKFRSQQHRQGRSHHFEMVIQLLGIMCVSCICWSPFLVTMASIGMNIQDKDSCERTLFTLRMATWNQILDPWVYILLRKAVLRNLYVCTRRCCGVHVISLHVWELKVAAISDLPVT
>GPRO_RAT
QTSLLSTGPNASNISDGQDNLTLPGSPPRTGSVSYINIIMPSVFGTICLLGIVGNSTVIFAVVKKSKCSNVPDIFIINLSVVDLLFLLGMPFMIHQLMGNGVWHFGETMCTLITAMDANSQFTSTYILTAMTIDRYLATVHPISTKFRKPSMATLVICLLWALSFISITPVWLYARLIP-FP-AVGCGIR--LPNPDTDLYWFTLYQFFLAFALPFVVITAAYVKILQRMTSSV--PASQRSIRLRTKRVTRTAIAICLVFFVCWAPYYVLQLTQLSISRPTFVYLYNAAIS-LGYANSCLNPFVYIVLCETFRKRLVLSVKPAAQGQLRTVSNAQESKGT-------
>ML1._HUMAN
---------MGPTLAVPTPYGCIGCKLPQPEYPPALIIFMFCAMVITIVVDLIGNSMVILAVTKNKKLRNSGNIFVVSLSVADMLVAIYPYPLMLHAMSIGGWDLSQLQCQMVGFITGLSVVGSIFNIVAIAINRYCYICHSLYERIFSVRNTCIYLVITWIMTVLAVLPNMYIGT-IEYDP-TYTCIFN------YLNNPVFTVTIVCIHFVLPLLIVGFCYVRIWTKVLAARD-AGQNPDNQLAEVRNFLTMFVIFLLFAVCWCPINVLTVLVAVSPKEKIPNWLYLAAYFIAYFNSCLNAVIYGLLNENFRREYWTIFHAMRHPIIFFPGLISARTLARARAHAR
>AVT_CATCO
-------------MGRIANQTTASNDTDPFGRNEEVAKMEITVLSVTFFVAVIGNLSVLLAMHNTKKKSSRMHLFIKHLSLADMVVAFFQVLPQLCWEITFRFYGPDFLCRIVKHLQVLGMFASTYMMVMMTLDRYIAICHPLKTLQQPTQRAYIMIGSTWLCSLLLSTPQYFIFSLSESY--VYDCWGH---FIEPWGIRAYITWITVGIFLIPVIILMICYGFICHSIWKNMIGVSSVTIISRAKLRTVKMTLVIVLAYIVCWAPFFIVQMWSVWDENFDSENAAVTLSALLASLNSCCNPWIYMLFSGHLLYDFLRCFPCCKKPRNMLQKEDSTLLTKLAAGRMT
>ACM3_PIG
N---------ISQAAGNFSSPNGTTSDPLGGHTIWQVVFIAFLTGILALVTIIGNILVIVAFKVNKQLKTVNNYFLLSLACADLIIGVISMNLFTTYIIMNRWALGNLACDLWLSIDYVASNASVMNLLVISFDRYFSITRPLYRAKRTTKRAGVMIGLAWVISFILWAPAILFWQYFV-VP-PGECFIQ------FLSEPTITFGTAIAAFYMPVTIMTILYWRIYKETEKRKTRTKRKRMSLIKEKKAAQTLSAILLAFIITWTPYNIMVLVNTFCDSC-IPKTYWNLGYWLCYINSTVNPVCYALCNKTFRTTFKMLLLCQCDKRKRRKQQYQHKRVPEQAL---
>MSHR_ALCAA
PVLGSQRRLLGSLNCTPPATFSLTLAPNRTGPQCLEVSIPDGLFLSLGLVSLVENVLVVAAIAKNRNLHSPMYYFICCLAVSDLLVSVSNVLETAVMLLLEAAAVVQQLDNVIDVLICGSMVSSLCFLGAIAMDRYISIFYALYHSVVTLPRAWRIIAAIWVASILTSLLFITYYY----------------------NNHTVVLLCLVGFFIAMLALMAILYVHMLARACQHIARKRQHPIHQGFGLKGAATLTILLGVFFLCWGPFFLHLSLIVLCPQHGCIFKNFNLFLALIICNAIVDPLIYAFRSQELRKTLQEVLQCSW-----------------------
>GP72_MOUSE
TGPNASSHFWANYTFSDWQNFVGRRRYGAESQNPTVKALLIVAYSFTIVFSLFGNVLVCHVIFKNQRMHSATSLFIVNLAVADIMITLLNTPFTLVRFVNSTWVFGKGMCHVSRFAQYCSLHVSALTLTAIAVDRHQVIMHPL-KPRISITKGVIYIAVIWVMATFFSLPHAICQKLFT--EVRSLCLPD-FPEPADLFWKYLDLATFILLYLLPLFIISVAYARVAKKLWLCGDVTEQYLALRRKKKTTVKMLVLVVVLFALCWFPLNCYVLLLSSKAI-HTNNALYFAFHWFAMSSTCYNPFIYCWLNENFRVELKALLSMCQRPPKPQEDRLPVAWTEKSHGRRA
>ACM5_RAT
--------MEGESYNESTVNGTPVNHQALERHGLWEVITIAVVTAVVSLMTIVGNVLVMISFKVNSQLKTVNNYYLLSLACADLIIGIFSMNLYTTYILMGRWVLGSLACDLWLALDYVASNASVMNLLVISFDRYFSITRPLYRAKRTPKRAGIMIGLAWLVSFILWAPAILCWQYLV-VP-PDECQIQ------FLSEPTITFGTAIAAFYIPVSVMTILYCRIYRETEKRNLSTKRKRMVLVKERKAAQTLSAILLAFIITWTPYNIMVLVSTFCDKC-VPVTLWHLGYWLCYVNSTINPICYALCNRTFRKTFKLLLLCRWKKKKVEEKLYW------------
>CKR5_PONPY
----MDYQVSSPTYDIDYYTSEPC---QKINVKQIAARLLPPLYSLVFIFGFVGNMLVILILINCKRLKSMTDIYLLNLAISDLFFLLTVPFWAHYAAAQ--WDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDRYLAIVHAVALKARTVTFGVVTSVITWVVAVFASLPGIIFTRSQK-EG-HYTCSSHFPYSQYQFWKNFQTLKIVILGLVLPLLVMVICYSGILKTL--------LRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNIVLLLNTFQEFFNRLDQAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCSIFA-------SSVY
>GASR_CANFA
GASLCRAGGALLNSSGAGNLSCEPPRLRGAGTRELELAIRVTLYAVIFLMSVGGNVLIIVVLGLSRRLRTVTNAFLLSLAVSDLLLAVACMPFTLLPNLMGTFIFGTVVCKAVSYLMGVSVSVSTLSLVAIALERYSAICRPLARVWQTRSHAARVIIATWMLSGLLMVPYPVYTAVQP---A-LQCVHR---WPSARVRQTWSVLLLLLLFFVPGVVMAVAYGLISRELYLGPGPPRPYQAKLLAKKRVVRMLLVIVVLFFLCWLPLYSANTWRAFDSSGALSGAPISFIHLLSYASACVNPLVYCFMHRRFRQACLETCARCCPRPPRARPRPLPSIASLSRLSYT
>TLR2_DROME
TLSTDQPAVGDVEDAAEDAAASMETGSFAFVVPWWRQVLWSILFGGMVIVATGGNLIVVWIVMTTKRMRTVTNYFIVNLSIADAMVSSLNVTFNYYYMLDSDWPFGEFYCKLSQFIAMLSICASVFTLMAISIDRYVAIIRPL-QPRMSKRCNLAIAAVIWLASTLISCPMMIIYR----NR--TVCYPEDGPTNHSTMESLYNILIIILTYFLPIVSMTVTYSRVGIELWGSTIGTPRQVENVRSKRRVVKMMIVVVLIFAICWLPFHSYFIITSCYPAIPFIQELYLAIYWLAMSNSMYNPIIYCWMNSRFRYGFKMVFRWCLFVRVGTEPFSRYSCSGSPDHNRI
>D2DR_MELGA
------MDPLNLSWYNTGDRNWSEPVNESSADQKPQYNYYAVLLTLLIFVIVFGNVLVCMAVSREKALQTTTNYLIVSLAVADLLVATLVMPWVVYLEVVGEWRFSRIHCDIFVTLDVMMCTASILNLCAISIDRYTAAAMPMNTRYSSKRRVTVMIACVWVLSFAISSPILFGLN---E----RECII---------ANPAFVVYSSVVSFYVPFIVTLLVYVQIYMVLRRRSTLMNRRKLSQQKEKKATQMLAIVLGVFIICWLPFFITHILNMHCDCN-IPPAMYSAFTWLGYVNSAVNPIIYTTFNIEFRKAFMKILHC-------------------------
>EDG2_BOVIN
QPQFTAMNEQQCFSNESIAFFYNRSGKYLATEWNTVTKLVMGLGITVCIFIMLANLLVMVAIYVNRRFHFPIYYLMANLAAADFFAG-LAYFYLMFNTGPNTRRLTVSTWLLRQGLIDTSLTVSVANLLAIAIERHITVFRMQLHARMSNRRVVVVIVVIWTMAIVMGAIPSVGWNCI-C--DIENCSNM------APLYSDSYLVFWAIFNLVTFVVMVVLYAHIFGYVRQRRMSSSGPRRNRDTMMSLLKTVVIVLGAFIICWTPGLVLLLLDVCCPQC-DVLAYEKFFLLLAEFNSAMNPIIYSYRDKEMSATFRQILCCQRSENTSGPTEGSNHTILAGVHSND
>NY2R_BOVIN
EEMKVDQFGPGHTTLPGELAPDSEPELIDSTKLIEVQVVLILAYCSIILLGVIGNSLVIHVVIKFKSMRTVTNFFIANLAVADLLVNTLCLPFTLTYTLMGEWKMGPVLCHLVPYAQGLAVQVSTITLTVIALDRHRCIVYHL-ESKISKQISFLIIGLAWGVSALLASPLAIFREYSLFE--IVACTEKWPGEEKGIYGTIYSLSSLLILYVLPLGIISFSYTRIWSKLKNHSPG-AAHDHYHQRRQKTTKMLVCVVVVFAVSWLPLHAFQLAVDIDSHVKEYKLIFTVFHIIAMCSTFANPLLYGWMNSNYRKAFLSAFRCEQRLDAIHSE---AKKHLQVTKNNG
>HH1R_RAT
----------MSFANTSSTFEDKMCEGNRTAMASPQLLPLVVVLSSISLVTVGLNLLVLYAVHSERKLHTVGNLYIVSLSVADLIVGAVVMPMNILYLIMTKWSLGRPLCLFWLSMDYVASTASIFSVFILCIDRYRSVQQPLYLRYRTKTRASATILGAWFFSFLWVIPILGWHHFM--EL-EDKCETD------FYNVTWFKIMTAIINFYLPTLLMLWFYVKIYKAVRRHLRSQYVSGLHLNRERKAAKQLGFIMAAFILCWIPYFIFFMVIAFCKSC-CSEPMHMFTIWLGYINSTLNPLIYPLCNENFKKTFKKILHIRS-----------------------"""

gpcr_aln = ArrayAlignment(data=gpcr_ungapped.split("\n"), moltype=PROTEIN)

myos_data = """>gi|107137|pir||A37102
LSRIITRIQA
>gi|11024712|ref|NP-060003.1|
LAQLITRTQA
>gi|11276950|pir||A59286
LSRIITRIQA
>gi|11276952|pir||A59293
LAQLITRTQA
>gi|11276954|pir||A59234
LSLIITRIQA
>gi|11276955|pir||A59236
LSSIFKLIQA
>gi|11321579|ref|NP-003793.1|
LVTLMTSTQA
>gi|11342672|ref|NP-002461.1|
LAKLITRTQA
>gi|1197168|dbj|BAA08111.1|
LSSIFKLIQA
>gi|12003423|gb|AAG43570.1|AF21
LVTLMTRTQA
>gi|12003425|gb|AAG43571.1|AF21
LVTLMTRTQA
>gi|12003427|gb|AAG43572.1|AF21
LVTLMTRTQA
>gi|12053672|emb|CAC20413.1|
LSRIITRIQA
>gi|12060489|dbj|BAB20630.1|
LSRIITRIQA
>gi|12657350|emb|CAC27776.1|
LASLVTLTQA
>gi|12657354|emb|CAC27778.1|
LAKLVTMTQA
>gi|127741|sp|P02563|MYH6-RAT
LSRIITRIQA
>gi|127748|sp|P02564|MYH7-RAT
LSRIITRIQA
>gi|127755|sp|P12847|MYH3-RAT
LAKLITRTQA
>gi|1289512|gb|AAC59911.1|
LSLIITRIQA
>gi|1289514|gb|AAC59912.1|
LSLIITRIQA
>gi|13431707|sp|Q28641|MYH4-RAB
LAQLITRTQA
>gi|13431711|sp|Q90339|MYSS-CYP
LALLVTMTQA
>gi|13431716|sp|Q9UKX2|MYH2-HUM
LAQLITRTQA
>gi|13431717|sp|Q9UKX3|MYHD-HUM
LVTLMTSTQA
>gi|13431724|sp|Q9Y623|MYH4-HUM
LAQLITRTQA
>gi|1346637|sp|P02565|MYH3-CHIC
LAQLITRTQA
>gi|13560269|dbj|BAB40920.1|
LAQLMTRTQA
>gi|13560273|dbj|BAB40922.1|
LSRIITRIQA
>gi|13638390|sp|P12882|MYH1-HUM
LAQLITRTQA
>gi|14017756|dbj|BAB47399.1|
LSLIITRIQA
>gi|15384839|emb|CAC59753.1|
LATLVTMTQA
>gi|1581130|prf||2116354A
LSRIITRIQA
>gi|1619328|emb|CAA27817.1|
LAKLITRTQA
>gi|16508127|gb|AAL17913.1|
LSRIITRIQA
>gi|1698895|gb|AAB37320.1|
LSRIITRIQA
>gi|17907763|dbj|BAB79445.1|
LAKILTMLQA
>gi|179508|gb|AAA51837.1|
LSRIITRIQA
>gi|179510|gb|AAA62830.1|
LSRIITRIQA
>gi|18859641|ref|NP-542766.1|
LSRIITRIQA
>gi|191618|gb|AAA37159.1|
LSRIITRIQA
>gi|191620|gb|AAA37160.1|
LSRIITRIQA
>gi|191622|gb|AAA37161.1|
LSRIITRIQA
>gi|191624|gb|AAA37162.1|
LSRIITRIQA
>gi|2119306|pir||I49464
LSRIITRIQA
>gi|2119307|pir||I48175
LSRIITRIQA
>gi|2119308|pir||I48153
LSRIITRIQA
>gi|212376|gb|AAA48972.1|
LAQLITRTQA
>gi|21623523|dbj|BAC00871.1|
LAALVGMVQA
>gi|21743235|dbj|BAB40921.2|
LAQLITRTQA
>gi|21907898|dbj|BAC05679.1|
LAQIITRTQA
>gi|21907900|dbj|BAC05680.1|
LAQIITRTQA
>gi|21907902|dbj|BAC05681.1|
LSRIITRIQA
>gi|219524|dbj|BAA00791.1|
LSRIITRIQA
>gi|22121649|gb|AAM88909.1|
LAQIITRTQA
>gi|23379831|gb|AAM88910.1|
LAQLITRTQA
>gi|2351219|dbj|BAA22067.1|
LSHLVTMTQA
>gi|2351221|dbj|BAA22068.1|
LVNLVTMTQA
>gi|2351223|dbj|BAA22069.1|
LALLVTMTQA
>gi|27764861|ref|NP-002462.1|
LSRIITRMQA
>gi|297024|emb|CAA79675.1|
LSRIITRMQA
>gi|3024204|sp|Q02566|MYH6-MOUS
LSRIITRIQA
>gi|3041706|sp|P13533|MYH6-HUMA
LSRIITRMQA
>gi|3041708|sp|P13540|MYH7-MESA
LSRIITRIQA
>gi|3043372|sp|P11055|MYH3-HUMA
LAKLITRTQA
>gi|34870884|ref|XP-213345.2|
LAQLITRTQA
>gi|34870892|ref|XP-340820.1|
LAQIITRTQA
>gi|37720046|gb|AAN71741.1|
LARILTGIQA
>gi|38091410|ref|XP-354614.1|
LAKLITRTQA
>gi|38091413|ref|XP-354615.1|
LAQLITRTQA
>gi|38177589|gb|AAF00096.2|AF11
LSLIISGIQA
>gi|38347761|dbj|BAD01606.1|
LSLLLTRTQA
>gi|38347763|dbj|BAD01607.1|
LSLLLTRTQA
>gi|38488753|ref|NP-942118.1|
LARILTGIQA
>gi|3915779|sp|P13539|MYH6-MESA
LSRIITRIQA
>gi|402372|gb|AAA62313.1|
LSRIITRIQA
>gi|402374|gb|AAB59701.1|
LSRIITRIQA
>gi|41350446|gb|AAS00505.1|
LAALVTMTQA
>gi|41386691|ref|NP-776542.1|
LAQLITRTQA
>gi|41386711|ref|NP-777152.1|
LSRIITRIQA
>gi|42476190|ref|NP-060004.2|
LAQLITRTQA
>gi|42662294|ref|XP-371398.2|
LAKVLTLLQA
>gi|45382109|ref|NP-990097.1|
LAKILTMIQA
>gi|45383005|ref|NP-989918.1|
LAKILTMLQA
>gi|45383668|ref|NP-989559.1|
LAQLITRTQA
>gi|4557773|ref|NP-000248.1|
LSRIITRIQA
>gi|45595719|gb|AAH67305.1|
LAQLITRTQA
>gi|476355|pir||A46762
LSRIITRIQA
>gi|4808809|gb|AAD29948.1|
LVTLMTSTQA
>gi|4808811|gb|AAD29949.1|
LAQLITRTQA
>gi|4808813|gb|AAD29950.1|
LAQLITRTQA
>gi|4808815|gb|AAD29951.1|
LAQLITRTQA
>gi|5360746|dbj|BAA82144.1|
LAQLITRTQA
>gi|5360748|dbj|BAA82145.1|
LAQLITRTQA
>gi|5360750|dbj|BAA82146.1|
LAQLITRTQA
>gi|547966|sp|P12883|MYH7-HUMAN
LSRIITRIQA
>gi|56655|emb|CAA34064.1|
LSRIITRIQA
>gi|6093461|sp|P79293|MYH7-PIG
LSRIITRIQA
>gi|6683485|dbj|BAA89233.1|
LAQLITRTQA
>gi|6708502|gb|AAD09454.2|
LAKIMTMLQC
>gi|7209643|dbj|BAA92289.1|
LATLVTMTQA
>gi|7248371|dbj|BAA92710.1|
LAKILTMIQA
>gi|7669506|ref|NP-005954.2|
LAQLITRTQA
>gi|8393804|ref|NP-058935.1|
LSRIITRIQA
>gi|8393807|ref|NP-058936.1|
LSRIITRIQA
>gi|86358|pir||A29320
LAQLITRTQA
>gi|88201|pir||S04090
LAKLITRTQA
>gi|92498|pir||S06005
LSRIITRIQA
>gi|92499|pir||S06006
LSRIITRIQA
>gi|92509|pir||A24922
LAKLITRTQA
>gi|940233|gb|AAA74199.1|
LAQLITRTQA
>gi|9800486|gb|AAF99314.1|AF272
LAQLITRTQA
>gi|9800488|gb|AAF99315.1|AF272
LAQLITRTQA
>gi|9971579|dbj|BAB12571.1|
LAALVTMTQA"""

myos_aln = ArrayAlignment(data=myos_data.split("\n"), moltype=PROTEIN)

# a randomly generated tree to use in tests
tree20_string = "(((0:0.5,1:0.5):0.5,(((2:0.5,3:0.5):0.5,(4:0.5,(5:0.5,6:0.5):0.5):0.5):0.5,((7:0.5,8:0.5):0.5,((9:0.5,((10:0.5,11:0.5):0.5,12:0.5):0.5):0.5,13:0.5):0.5):0.5):0.5):0.5,(((14:0.5,(15:0.5,16:0.5):0.5):0.5,17:0.5):0.5,(18:0.5,19:0.5):0.5):0.5);"

if __name__ == "__main__":
    main()
