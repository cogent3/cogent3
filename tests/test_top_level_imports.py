import importlib
import subprocess
import sys

import pytest

import cogent3


@pytest.mark.parametrize(
    "name",
    sorted(cogent3._import_mapping),
)
def test_lazy_import(name):
    """all names in _import_mapping are accessible via cogent3.<name>"""
    attr = getattr(cogent3, name)
    module_path = cogent3._import_mapping[name]
    module = importlib.import_module(f"cogent3.{module_path}")
    assert attr is getattr(module, name)


@pytest.mark.parametrize(
    "name,expected_type",
    [
        ("Alignment", type),
        ("SequenceCollection", type),
        ("Sequence", type),
        ("Table", type),
        ("MolType", type),
        ("GeneticCode", type),
        ("PhyloNode", type),
    ],
)
def test_toplevel_types(name, expected_type):
    """core return types are accessible as cogent3.<Type>"""
    attr = getattr(cogent3, name)
    assert isinstance(attr, expected_type), f"cogent3.{name} is not a {expected_type}"


def test_getattr_returns_cached_module():
    # os is bound at the top of cogent3/__init__.py, so __getattr__
    # short-circuits via the globals() cache without consulting _import_mapping
    import os as _os

    assert cogent3.__getattr__("os") is _os


def test_getattr_fallback_imports_stdlib_module():
    # names not in _import_mapping should fall through to __import__
    import json as _json

    assert "json" not in cogent3._import_mapping
    assert cogent3.__getattr__("json") is _json


def test_getattr_unknown_name_raises():
    with pytest.raises(AttributeError):
        cogent3.__getattr__("definitely_not_a_real_module_xyz")


def test_warnings_env_var_applied():
    """COGENT3_WARNINGS env var registers a simplefilter on import."""
    import os

    env = {**os.environ, "COGENT3_WARNINGS": "error"}
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import cogent3\n"
            "import warnings\n"
            "warnings.warn('triggered', UserWarning)\n",
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
        env=env,
    )
    assert result.returncode != 0
    assert "UserWarning" in result.stderr


def test_profile_first_import_no_circular_error():
    """Importing cogent3.core.profile first must not trigger a circular import."""
    result = subprocess.run(
        [sys.executable, "-c", "import cogent3.core.profile"],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, (
        f"Importing cogent3.core.profile failed:\n{result.stderr}"
    )


def test_profile_import_then_make_seq():
    """After importing profile first, MolType.make_seq must still work."""
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import cogent3.core.profile\n"
            "from cogent3.core.moltype import DNA\n"
            "seq = DNA.make_seq(seq='ACGT', name='t')\n"
            "assert str(seq) == 'ACGT'\n",
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, (
        f"make_seq after profile-first import failed:\n{result.stderr}"
    )


@pytest.mark.slow
@pytest.mark.parametrize(
    "module",
    [
        "cogent3.core.alignment",
        "cogent3.core.sequence",
        "cogent3.core.tree",
        "cogent3.core.table",
        "cogent3.core.moltype",
        "cogent3.core.profile",
        "cogent3.evolve.models",
        "cogent3.evolve.distance",
        "cogent3.evolve.fast_distance",
        "cogent3.evolve.parameter_controller",
        "cogent3.align.pairwise",
        "cogent3.align.progressive",
        "cogent3.app",
        "cogent3.app.io",
        "cogent3.app.evo",
        "cogent3.parse.fasta",
        "cogent3.parse.genbank",
        "cogent3.parse.cogent3_json",
        "cogent3.draw.dotplot",
        "cogent3.draw.dendrogram",
        "cogent3.phylo.nj",
        "cogent3.maths.optimisers",
    ],
)
def test_no_circular_imports(module):
    """Each subpackage must be importable in a fresh process without circular import errors."""
    result = subprocess.run(
        [sys.executable, "-c", f"import {module}"],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, f"Importing {module} failed:\n{result.stderr}"
