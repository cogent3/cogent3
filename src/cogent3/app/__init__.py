import importlib
import inspect

from warnings import filterwarnings


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

__all__ = ["align", "composable", "dist", "evo", "io", "sample", "translate", "tree"]


def _get_app_attr(name, obj, mod, is_composable):
    """returns app details for display"""
    _input_types = getattr(obj, "_input_types", [])
    _output_types = getattr(obj, "_output_types", [])
    _data_types = getattr(obj, "_data_types", [])

    if _input_types:
        if type(_input_types) == str:
            _input_types = [_input_types]
        in_type = [{None: ""}.get(e, e) for e in _input_types]
    else:
        in_type = []

    if _output_types:
        if type(_output_types) == str:
            _output_types = [_output_types]
        out_type = [{None: ""}.get(e, e) for e in _output_types]
    else:
        out_type = []

    if _data_types:
        if type(_data_types) == str:
            _data_types = [_data_types]
        data_type = [{None: ""}.get(e, e) for e in _data_types]
    else:
        data_type = []

    row = [
        mod.__name__,
        name,
        is_composable,
        obj.__doc__,
        ", ".join(in_type),
        ", ".join(out_type),
        ", ".join(data_type),
    ]
    return row


def available_apps():
    """returns table of all available apps"""
    from cogent3.util.table import Table
    from .composable import Composable, user_function

    # excluding composable, find all class
    rows = []
    for m in __all__:
        if m == "composable":
            continue
        mod = importlib.import_module(f"{__name__}.{m}")
        for name, obj in inspect.getmembers(mod, inspect.isclass):
            if name.startswith("_"):
                continue
            if obj.__module__ == mod.__name__:
                is_composable = issubclass(obj, Composable)
                rows.append(_get_app_attr(name, obj, mod, is_composable))

    mod = importlib.import_module(f"{__name__}.composable")
    rows.append(
        _get_app_attr(
            "user_function", user_function, mod, issubclass(user_function, Composable)
        )
    )
    header = ["module", "name", "composable", "doc", "inputs", "outputs", "data type"]
    table = Table(header, rows)
    # todo the inputs/outputs/data type columns no longer being populated as
    # these attributes are instance attributes and not resolved at present
    table = table.get_columns(header)
    return table
