"""testing the default import"""
import os

from tempfile import TemporaryDirectory
from unittest import TestCase, main

from cogent3 import available_apps
from cogent3.app import align, dist, evo, io, sample, translate, tree


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def _get_all_composables(tmp_dir_name):
    test_model1 = evo.model("HKY85")
    test_model2 = evo.model("GN")
    test_hyp = evo.hypothesis(test_model1, test_model2)
    test_num_reps = 100

    applications = [
        align.align_to_ref(),
        align.progressive_align(model="GY94"),
        dist.fast_slow_dist(moltype="dna", fast_calc="hamming"),
        evo.ancestral_states(),
        evo.bootstrap(hyp=test_hyp, num_reps=test_num_reps),
        evo.hypothesis(test_model1, test_model2),
        evo.model("GN"),
        evo.tabulate_stats(),
        sample.fixed_length(100),
        sample.min_length(100),
        io.write_db(tmp_dir_name, create=True),
        io.write_json(tmp_dir_name, create=True),
        io.write_seqs(tmp_dir_name, create=True),
        sample.omit_bad_seqs(),
        sample.omit_degenerates(),
        sample.omit_duplicated(),
        sample.take_codon_positions(1),
        sample.take_named_seqs(),
        sample.take_n_seqs(2),
        sample.trim_stop_codons(gc=1),
        translate.select_translatable(),
        tree.quick_tree(),
        tree.scale_branches(),
        tree.uniformize_tree(),
    ]
    return applications


class TestAvailableApps(TestCase):
    def test_available_apps(self):
        """available_apps returns a table"""
        from cogent3.util.table import Table

        apps = available_apps()
        self.assertIsInstance(apps, Table)
        self.assertTrue(apps.shape[0] > 10)

    def test_composable_pairwise_applications(self):
        """Properly compose two composable applications"""
        from cogent3.app.composable import Composable

        with TemporaryDirectory(dir=".") as dirname:
            applications = _get_all_composables(os.path.join(dirname, "delme"))

            for app in applications:
                self.assertIsInstance(app, Composable)

            composable_application_tuples = [
                (app1, app2)
                for app1 in applications
                for app2 in applications
                if app1 != app2 and app1._output_types & app2._input_types != set()
            ]

            for composable_application_tuple in composable_application_tuples:
                composable_application_tuple[0].disconnect()
                composable_application_tuple[1].disconnect()
                # Compose two composable applications, there should not be exceptions.
                composable_application_tuple[0] + composable_application_tuple[1]

            for app in applications:
                if hasattr(app, "data_store"):
                    app.data_store.close()

    def test_incompatible_pairwise_applications(self):
        """Properly identify two incompatible applications"""
        from cogent3.app.composable import Composable

        with TemporaryDirectory(dir=".") as dirname:
            applications = _get_all_composables(os.path.join(dirname, "delme"))

            for app in applications:
                self.assertIsInstance(app, Composable)

            incompatible_application_tuples = [
                (app1, app2)
                for app1 in applications
                for app2 in applications
                if app1 != app2 and app1._output_types & app2._input_types == set()
            ]

            for incompatible_application_tuple in incompatible_application_tuples:
                incompatible_application_tuple[0].disconnect()
                incompatible_application_tuple[1].disconnect()

                # Compose two incompatible applications, there should be exceptions.
                with self.assertRaises(TypeError):
                    res = (
                        incompatible_application_tuple[0]
                        + incompatible_application_tuple[1]
                    )

            for app in applications:
                if hasattr(app, "data_store"):
                    app.data_store.close()


if __name__ == "__main__":
    main()
