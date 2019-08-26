"""testing the default import"""
import os

from unittest import TestCase, main

from cogent3 import available_apps
from cogent3.app import align, evo, io, sample, translate, tree


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.8.23a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestAvalableApps(TestCase):
    def test_available_apps(self):
        """available_apps returns a table"""
        from cogent3.util.table import Table

        apps = available_apps()
        self.assertIsInstance(apps, Table)
        self.assertTrue(apps.shape[0] > 10)

    def test_composable_pairwise_applications(self):
        """Properly compose two composable applications"""
        from cogent3.app.composable import Composable

        # Initialise some data for the application configurations.
        test_model1 = evo.model("HKY85")
        test_model2 = evo.model("HKY85", name="hky85")
        test_hyp = evo.hypothesis(test_model1, test_model2)
        test_num_reps = 1000

        applications = [
            align.align_to_ref(),
            align.progressive_align(model="nucleotide"),
            evo.ancestral_states(),
            evo.bootstrap(hyp=test_hyp, num_reps=test_num_reps),
            evo.hypothesis(test_model1, test_model2),
            evo.model("HKY85"),
            evo.tabulate_stats(),
            io.write_db(os.getcwd()),
            io.write_json(os.getcwd()),
            sample.fixed_length(100),
            sample.min_length(100),
            io.write_seqs(os.getcwd()),
            sample.omit_bad_seqs(),
            sample.omit_degenerates(),
            sample.omit_duplicated(),
            sample.take_codon_positions(1),
            sample.take_named_seqs(),
            sample.trim_stop_codons(gc=1),
            translate.select_translatable(),
            tree.quick_tree(),
            tree.scale_branches(),
            tree.uniformize_tree(),
        ]

        for app in applications:
            self.assertIsInstance(app, Composable)

        composable_application_tuples = [
            (fir_element, sec_element)
            for fir_element in applications
            for sec_element in applications
            if fir_element != sec_element
            and fir_element._output_type & sec_element._input_type != set()
        ]

        for composable_application_tuple in composable_application_tuples:
            composable_application_tuple[0].disconnect()
            composable_application_tuple[1].disconnect()
            # Compose two composable applications, there should not be exceptions.
            res = composable_application_tuple[0] + composable_application_tuple[1]

    def test_uncomposable_pairwise_applications(self):
        """Properly identify two uncomposable applications"""
        from cogent3.app.composable import Composable

        # Initialise some data for the application configurations.
        test_model1 = evo.model("HKY85")
        test_model2 = evo.model("GN")
        test_hyp = evo.hypothesis(test_model1, test_model2)
        test_num_reps = 100

        applications = [
            align.align_to_ref(),
            align.progressive_align(model="GY94"),
            evo.ancestral_states(),
            evo.bootstrap(hyp=test_hyp, num_reps=test_num_reps),
            evo.hypothesis(test_model1, test_model2),
            evo.model("GN"),
            evo.tabulate_stats(),
            io.write_db(os.getcwd()),
            io.write_json(os.getcwd()),
            sample.fixed_length(100),
            sample.min_length(100),
            io.write_seqs(os.getcwd()),
            sample.omit_bad_seqs(),
            sample.omit_degenerates(),
            sample.omit_duplicated(),
            sample.take_codon_positions(1),
            sample.take_named_seqs(),
            sample.trim_stop_codons(gc=1),
            translate.select_translatable(),
            tree.quick_tree(),
            tree.scale_branches(),
            tree.uniformize_tree(),
        ]

        for app in applications:
            self.assertIsInstance(app, Composable)

        uncomposable_application_tuples = [
            (fir_element, sec_element)
            for fir_element in applications
            for sec_element in applications
            if fir_element != sec_element
            and fir_element._output_type & sec_element._input_type == set()
        ]

        for uncomposable_application_tuple in uncomposable_application_tuples:
            uncomposable_application_tuple[0].disconnect()
            uncomposable_application_tuple[1].disconnect()

            # Compose two uncomposable applications, there should be exceptions.
            with self.assertRaises(TypeError):
                res = (
                    uncomposable_application_tuple[0]
                    + uncomposable_application_tuple[1]
                )


if __name__ == "__main__":
    main()
