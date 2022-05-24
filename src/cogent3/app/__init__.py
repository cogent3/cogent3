import importlib
import inspect


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

__all__ = ["align", "composable", "dist", "evo", "io", "sample", "translate", "tree"]


def _get_app_attr(name, obj, mod, is_composable):
    """returns app details for display"""

    _types = {"_input_types": [], "_output_types": [], "_data_types": []}

    for tys in _types.keys():
        types = getattr(obj, tys, None) or []
        types = [types] if type(types) == str else types
        _types[tys] = [{None: ""}.get(e, e) for e in types]

    row = [
        mod.__name__,
        name,
        is_composable,
        obj.__doc__,
        ", ".join(_types["_input_types"]),
        ", ".join(_types["_output_types"]),
        ", ".join(_types["_data_types"]),
    ]
    return row


def available_apps():
    """returns Table listing the available apps"""
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
    return Table(header, rows)
