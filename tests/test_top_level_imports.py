import importlib

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
