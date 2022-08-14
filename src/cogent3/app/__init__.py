import importlib
import inspect


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

__all__ = ["align", "composable", "dist", "evo", "io", "sample", "translate", "tree"]


def _get_app_attr(name):
    """returns app details for display"""
    from .composable import Composable, is_composable, user_function

    modname, name = name.rsplit(".", maxsplit=1)
    mod = importlib.import_module(modname)
    obj = getattr(mod, name)
    iscomposable = (
        issubclass(user_function, Composable)
        if name == "user_function"
        else is_composable(f"{obj.__module__}.{name}")
    )

    _types = {"_data_types": [], "_return_types": []}

    for tys in _types:
        types = getattr(obj, tys, None) or []
        types = [types] if type(types) == str else types
        _types[tys] = [{None: ""}.get(e, e) for e in types]

    return [
        mod.__name__,
        name,
        iscomposable,
        obj.__doc__,
        _types["_data_types"],
        _types["_return_types"],
    ]


def available_apps():
    """returns Table listing the available apps"""
    from cogent3.util.table import Table

    # excluding composable, find all class
    rows = []
    for m in __all__:
        if m == "composable":
            continue
        mod = importlib.import_module(f"{__name__}.{m}")
        rows.extend(
            _get_app_attr(f"{obj.__module__}.{name}")
            for name, obj in inspect.getmembers(mod, inspect.isclass)
            if not name.startswith("_") and obj.__module__ == mod.__name__
        )

    mod = importlib.import_module(f"{__name__}.composable")
    rows.append(_get_app_attr(f"{mod.__name__}.user_function"))
    header = ["module", "name", "composable", "doc", "outputs", "data type"]
    return Table(header, rows)
