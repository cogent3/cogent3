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
        ", ".join(sorted(_types["_data_types"])),
        ", ".join(sorted(_types["_return_types"])),
    ]


def available_apps():
    """returns Table listing the available apps"""
    import cogent3.app as apps

    from cogent3.util.table import Table

    from .composable import Composable, __app_registry, user_function

    # registration of apps does not happen until their modules are imported
    for name in __all__:
        importlib.import_module(f"cogent3.app.{name}")

    rows = [_get_app_attr(app, is_comp) for app, is_comp in __app_registry.items()]
    header = ["module", "name", "composable", "doc", "input type", "output type"]
    return Table(header=header, data=rows)
