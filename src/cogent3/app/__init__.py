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


def _get_app_attr(name, is_composable):
    """returns app details for display"""

    modname, name = name.rsplit(".", maxsplit=1)
    mod = importlib.import_module(modname)
    obj = getattr(mod, name)

    _types = {"_data_types": [], "_return_types": []}

    for tys in _types:
        types = getattr(obj, tys, None) or []
        types = [types] if type(types) == str else types
        _types[tys] = [{None: ""}.get(e, e) for e in types]

    return [
        mod.__name__,
        name,
        is_composable,
        obj.__doc__,
        _types["_data_types"],
        _types["_return_types"],
    ]


def available_apps():
    """returns Table listing the available apps"""
    from cogent3.util.table import Table

    from .composable import Composable, __app_registry, user_function

    rows = [_get_app_attr(app, True) for app in __app_registry]

    mod = importlib.import_module(f"{__name__}.composable")
    rows.append(
        _get_app_attr(
            f"{mod.__name__}.user_function", issubclass(user_function, Composable)
        )
    )
    header = ["module", "name", "composable", "doc", "outputs", "data type"]
    return Table(header, rows)
