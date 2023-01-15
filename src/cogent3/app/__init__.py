from __future__ import annotations

import contextlib
import importlib
import inspect
import re
import textwrap

from .io_new import open_data_store


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.10.31a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

__all__ = [
    "align",
    "composable",
    "dist",
    "evo",
    "io_new",
    "sample",
    "translate",
    "tree",
]


def _get_app_attr(name, is_composable):
    """returns app details for display"""

    modname, name = name.rsplit(".", maxsplit=1)
    mod = importlib.import_module(modname)
    obj = getattr(mod, name)

    _types = _make_types(obj)

    return [
        mod.__name__,
        name,
        is_composable,
        _doc_summary(obj.__doc__ or ""),
        ", ".join(sorted(_types["_data_types"])),
        ", ".join(sorted(_types["_return_types"])),
    ]


def _make_types(app) -> dict:
    """returns type hints for the input and output"""
    _types = {"_data_types": [], "_return_types": []}
    for tys in _types:
        types = getattr(app, tys, None) or []
        types = [types] if type(types) == str else types
        _types[tys] = [{None: ""}.get(e, e) for e in types]
    return _types


def available_apps():
    """returns Table listing the available apps"""
    from cogent3.util.table import Table

    from .composable import __app_registry

    # registration of apps does not happen until their modules are imported
    for name in __all__:
        importlib.import_module(f"cogent3.app.{name}")

    # exclude apps from deprecated modules
    deprecated = ["cogent3.app.io."]

    rows = []
    for app, is_comp in __app_registry.items():
        if any(app.startswith(d) for d in deprecated):
            continue

        with contextlib.suppress(AttributeError):
            # probably a local scope issue in testing!
            rows.append(_get_app_attr(app, is_comp))

    header = ["module", "name", "composable", "doc", "input type", "output type"]
    return Table(header=header, data=rows)


_get_param = re.compile('(?<=").+(?=")')


def _make_signature(app: type) -> str:
    init_sig = inspect.signature(app.__init__)
    params = ", ".join(
        [f'"{app.__name__}"']
        + [
            _get_param.findall(repr(v))[0]
            for k, v in init_sig.parameters.items()
            if k != "self"
        ]
    )
    sig = f"{app.__name__} = get_app({params})"
    subsequent_indent = " " * (sig.find("(") + 1)
    return "\n".join(textwrap.wrap(sig, subsequent_indent=subsequent_indent))


def _doc_summary(doc):
    """return first para of docstring"""
    result = []
    for line in doc.splitlines():
        if line := line.strip():
            result.append(line)
        else:
            break
    return " ".join(result)


def _get_app_matching_name(name: str):
    """name can include module name"""
    table = available_apps()
    if "." in name:
        modname, name = name.rsplit(".", maxsplit=1)
        table = table.filtered(
            lambda x: x[0].endswith(modname) and x[1] == name,
            columns=["module", "name"],
        )
    else:
        table = table.filtered(lambda x: name == x, columns="name")

    if table.shape[0] == 0:
        raise NameError(f"no app matching name {name!r}")
    elif table.shape[0] > 1:
        raise NameError(f"too many apps matching name {name!r},\n{table}")
    modname, _ = table.tolist(columns=["module", "name"])[0]
    mod = importlib.import_module(modname)
    return getattr(mod, name)


def get_app(name: str, *args, **kwargs):
    """returns app instance, use app_help() to display arguments

    Parameters
    ----------
    name
        app name, e.g. 'minlength', or can include module information,
        e.g. 'cogent3.app.sample.minlength' or 'sample.minlength'. Use the
        latter (qualified class name) style when there are multiple matches
        to name.
    *args, **kwargs
        positional and keyword arguments passed through to the app

    Returns
    -------
    cogent3 app instance

    Raises
    ------
    NameError when multiple apps have the same name. In that case use a
    qualified class name, as shown above.
    """
    return _get_app_matching_name(name)(*args, **kwargs)


def _make_head(text: str) -> list[str]:
    """makes a restructured text formatted header"""
    return [text, "-" * len(text)]


def _clean_params_docs(text: str) -> str:
    """remove unnecessary indentation"""
    text = text.splitlines(keepends=False)
    prefix = re.compile(r"^\s{8}")  # expected indentation of constructor doc
    doc = []
    for line in text:
        line = prefix.sub("", line)
        if line.strip():
            doc.append(line)

    return "\n".join(doc)


def _clean_overview(text: str) -> str:
    text = text.split()
    return "\n".join(textwrap.wrap(" ".join(text), break_long_words=False))


def app_help(name: str):
    """displays help for the named app

    Parameters
    ----------
    name
        app name, e.g. 'minlength', or can include module information,
        e.g. 'cogent3.app.sample.minlength' or 'sample.minlength'. Use the
        latter (qualified class name) style when there are multiple matches
        to name.
    """
    app = _get_app_matching_name(name)
    docs = []
    app_doc = app.__doc__ or ""
    if app_doc.strip():
        docs.extend(_make_head("Overview") + [_clean_overview(app_doc), ""])

    docs.extend(_make_head("Options for making the app") + [_make_signature(app)])
    init_doc = app.__init__.__doc__ or ""
    if init_doc.strip():
        docs.extend(["", _clean_params_docs(init_doc)])

    types = _make_types(app)
    docs.extend([""] + _make_head("Input type") + [", ".join(types["_data_types"])])
    docs.extend([""] + _make_head("Output type") + [", ".join(types["_return_types"])])

    print("\n".join(docs))
