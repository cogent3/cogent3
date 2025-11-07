from pathlib import Path

import pytest

from cogent3 import app_help, available_apps, get_app, open_data_store
from cogent3.app.composable import LOADER, WRITER, is_app
from cogent3.core.table import Table

try:
    import piqtree  # noqa: F401

    has_piqtree = True
except ImportError:
    has_piqtree = False


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


def test_available_apps():
    """available_apps returns a table"""
    from cogent3.core.table import Table

    apps = available_apps()
    assert isinstance(apps, Table)
    assert apps.shape[0] > 10


def _get_incompat_app_pairs(tmp_path):
    """Generate all incompatible application pairs"""
    applications = _get_all_composables(tmp_path / "delme")
    return [
        (app1, app2)
        for app1 in applications
        for app2 in applications
        if app1.app_type is WRITER
        or (
            app2.app_type is LOADER
            and app1 != app2
            and not app1._return_types & app2._data_types
            and not app1._return_types & {"SerialisableType", "IdentifierType"}
        )
    ]


def _get_compat_app_pairs(tmp_path):
    """Generate all incompatible application pairs"""
    applications = _get_all_composables(tmp_path / "delme")
    return [
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


def pytest_generate_tests(metafunc):
    """Dynamically generate test parameters"""
    if "incompat_app_pair" in metafunc.fixturenames:
        tmp_path = metafunc.config._tmp_path_factory.mktemp("test_apps")
        pairs = _get_incompat_app_pairs(tmp_path)
        metafunc.parametrize(
            "incompat_app_pair",
            pairs,
            ids=[
                f"{app1.__class__.__name__}-{app2.__class__.__name__}"
                for app1, app2 in pairs
            ],
        )
    elif "compat_pair" in metafunc.fixturenames:
        tmp_path = metafunc.config._tmp_path_factory.mktemp("test_apps")
        pairs = _get_compat_app_pairs(tmp_path)
        metafunc.parametrize(
            "compat_pair",
            pairs,
            ids=[
                f"{app1.__class__.__name__}-{app2.__class__.__name__}"
                for app1, app2 in pairs
            ],
        )


def test_incompatible_pairwise_applications(incompat_app_pair):
    """Properly identify two incompatible applications"""
    app_a, app_b = incompat_app_pair
    assert is_app(app_a)
    assert is_app(app_b)
    err_type = ValueError if app_a is app_b else TypeError
    # make sure they're not connected
    app_a.disconnect()
    app_b.disconnect()

    # Compose two incompatible applications, there should be exceptions.
    with pytest.raises(err_type):
        app_a + app_b


def test_composable_pairwise_applications(compat_pair):
    """Properly compose two composable applications"""
    app_a, app_b = compat_pair
    assert is_app(app_a)
    assert is_app(app_b)

    app_a.disconnect()
    app_b.disconnect()
    # Compose two composable applications, there should not be exceptions.
    _ = app_a + app_b


@pytest.mark.parametrize("name", ["sample.min_length", "min_length"])
def test_get_app(name):
    app = get_app(name, 500)
    assert app.__class__.__name__.endswith(name.split(".")[-1])


def test_get_app_kwargs():
    # when an app has a name kwarg
    # we should still be able to use get_app!
    _ = get_app("model", "F81", name="F81-model")


def test_app_help(capsys):
    app_help("concat")
    got = capsys.readouterr().out
    assert "Options" in got
    assert got.count("SerialisableType") == 1  # output type


@pytest.mark.parametrize(
    "app_name",
    ("bootstrap", "from_primitive", "load_db", "take_named_seqs")[:1],
)
def test_app_help_signature(capsys, app_name):
    from cogent3.app import _get_app_matching_name, _make_signature

    with pytest.raises(ValueError, match="app cannot be None"):
        _make_signature(None)

    got = _make_signature(_get_app_matching_name(app_name))
    # app name is in quotes
    assert f"{app_name!r}" in got
    # check split across multiple lines if long signature
    if len(got) > 70:
        assert got.count("\n") > 1


def test_available_apps_filter():
    """available apps can be filtered by name"""
    app_name_filter: str = "load"
    filtered_apps = available_apps(app_name_filter)
    assert isinstance(filtered_apps, Table)
    assert len(filtered_apps) > 0
    # check every returned table row 'name' has filter in it
    assert sum(app_name_filter in n for n in filtered_apps.columns["name"]) == len(
        filtered_apps,
    )


def test_available_apps_package():
    """available apps lists just the package"""
    filtered_apps = available_apps("model")
    assert "cogent3" in filtered_apps.columns["package"]
    assert all("." not in pkg for pkg in filtered_apps.columns["package"])


_package_licenses = [("typing_extensions", "PSF-2"), ("numpy", "BSD")]
if has_piqtree:
    _package_licenses.append(("piqtree", "GPL"))


@pytest.mark.parametrize(("package", "license"), _package_licenses)
def test_available_apps_license(package, license):
    # table contains the correct license
    from cogent3.app import _get_licenses

    # "typing_extensions" has no trove classifier
    # rich has license but not in trove
    got = _get_licenses(package)
    assert license in got


def test_available_apps_license_col():
    available = available_apps()
    assert "licenses" in available.columns
    assert "BSD" in available.columns["licenses"]
