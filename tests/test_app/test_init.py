"""testing the default import"""
import os

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase, main

import pytest

from cogent3 import app_help, available_apps, get_app, open_data_store
from cogent3.app.composable import (
    LOADER,
    WRITER,
    __app_registry,
    define_app,
    is_composable,
)
from cogent3.util.misc import get_object_provenance
from cogent3.util.table import Table


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.10.31a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def _get_all_composables(tmp_dir_name):
    tmp_dir_name = Path(tmp_dir_name)
    test_model1 = get_app("model", "HKY85")
    test_model2 = get_app("model", "GN")
    test_hyp = get_app("hypothesis", test_model1, test_model2)
    test_num_reps = 100
    return [
        get_app(
            "align_to_ref",
        ),
        get_app("progressive_align", model="GY94"),
        get_app("fast_slow_dist", moltype="dna", fast_calc="hamming"),
        get_app("ancestral_states"),
        get_app("bootstrap", hyp=test_hyp, num_reps=test_num_reps),
        get_app("hypothesis", test_model1, test_model2),
        get_app("model", "GN"),
        get_app("tabulate_stats"),
        get_app("fixed_length", 100),
        get_app("sample.min_length", 100),
        get_app("write_db", open_data_store(tmp_dir_name / "delme.sqlitedb", mode="w")),
        get_app(
            "write_json",
            open_data_store(tmp_dir_name / "json", suffix="json", mode="w"),
        ),
        get_app(
            "write_seqs",
            open_data_store(tmp_dir_name / "fasta", suffix="fasta", mode="w"),
        ),
        get_app("omit_bad_seqs"),
        get_app("omit_degenerates"),
        get_app("omit_duplicated"),
        get_app("take_codon_positions", 1),
        get_app("take_named_seqs"),
        get_app("take_n_seqs", 2),
        get_app("trim_stop_codons", gc=1),
        get_app("select_translatable"),
        get_app("quick_tree"),
        get_app("scale_branches"),
        get_app("uniformize_tree"),
    ]


class TestAvailableApps(TestCase):
    def test_available_apps(self):
        """available_apps returns a table"""
        from cogent3.util.table import Table

        apps = available_apps()
        self.assertIsInstance(apps, Table)
        self.assertTrue(apps.shape[0] > 10)

    def test_composable_pairwise_applications(self):
        """Properly compose two composable applications"""

        with TemporaryDirectory(dir=".") as dirname:
            applications = _get_all_composables(os.path.join(dirname, "delme"))
            for app in applications:
                self.assertTrue(is_composable(app), msg=app)

            composable_application_tuples = [
                (app1, app2)
                for app1 in applications
                for app2 in applications
                if app1 != app2
                and (
                    app1._return_types & app2._data_types
                    or app1._return_types & {"SerialisableType", "IdentifierType"}
                )
                and app1.app_type is not WRITER
                and app2.app_type is not LOADER
            ]

            for app_a, app_b in composable_application_tuples:
                app_a.disconnect()
                app_b.disconnect()
                # Compose two composable applications, there should not be exceptions.
                app_a + app_b

    def test_incompatible_pairwise_applications(self):
        """Properly identify two incompatible applications"""

        with TemporaryDirectory(dir=".") as dirname:
            applications = _get_all_composables(os.path.join(dirname, "delme"))
            for app in applications:
                self.assertTrue(is_composable(app))

            incompatible_application_tuples = [
                (app1, app2)
                for app1 in applications
                for app2 in applications
                if app1.app_type is WRITER
                or app2.app_type is LOADER
                and app1 != app2
                and not app1._return_types & app2._data_types
                and not app1._return_types & {"SerialisableType", "IdentifierType"}
            ]

            for app_a, app_b in incompatible_application_tuples:
                err_type = ValueError if app_a is app_b else TypeError
                app_a.disconnect()
                app_b.disconnect()

                # Compose two incompatible applications, there should be exceptions.
                with self.assertRaises(err_type):
                    app_a + app_b


def test_available_apps_local():
    """available_apps robust to local scope apps"""

    @define_app
    def dummy(val: int) -> int:
        return val

    apps = available_apps()
    assert isinstance(apps, Table)
    __app_registry.pop(get_object_provenance(dummy), None)


@pytest.mark.parametrize("name", ("sample.min_length", "min_length"))
def test_get_app(name):
    __app_registry.pop(get_object_provenance(min_length), None)
    app = get_app(name, 500)
    assert app.__class__.__name__.endswith(name.split(".")[-1])


@define_app
def min_length(val: int) -> int:
    return val


def test_get_app_fail():

    __app_registry[get_object_provenance(min_length)] = True

    with pytest.raises(NameError):
        _ = get_app("min_length", 500)

    __app_registry.pop(get_object_provenance(min_length), None)


def test_app_help(capsys):
    app_help("compress")
    got = capsys.readouterr().out
    assert "Options" in got
    assert got.count("bytes") >= 2  # both input and output types are bytes


@define_app
class blah:
    def __init__(self):
        ...

    def main(self, val: int) -> int:
        return val


@pytest.mark.parametrize(
    "app_doc,init_doc", ((None, None), ("text", None), (None, "text"), ("text", "text"))
)
def test_app_help_no_docs(capsys, app_doc, init_doc):
    blah.__doc__ = app_doc
    blah.__init__.__doc__ = init_doc
    app_help("blah")
    got = capsys.readouterr().out
    if app_doc:
        assert "Overview" in got

    if init_doc:
        assert "Options" in got


if __name__ == "__main__":
    main()
