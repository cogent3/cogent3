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


def test_profile_first_import_no_circular_error():
    """Importing cogent3.core.profile first must not trigger a circular import."""
    result = subprocess.run(
        [sys.executable, "-c", "import cogent3.core.profile"],
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
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, (
        f"make_seq after profile-first import failed:\n{result.stderr}"
    )
