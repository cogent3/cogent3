import importlib
import inspect

from warnings import filterwarnings


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.07.10a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

__all__ = ["align", "composable", "dist", "io", "sample", "translate", "tree"]


def available_apps():
    """returns table of all available apps"""
    from cogent3.util.table import Table
    from .composable import Composable

    # exclude composable, find all class
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
                rows.append([mod.__name__, name, is_composable, obj.__doc__])
    header = ["module", "name", "composable", "doc"]
    table = Table(header, rows)
    return table
