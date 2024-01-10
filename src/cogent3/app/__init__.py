from __future__ import annotations

import contextlib
import importlib
import inspect
import re
import textwrap

from cogent3.util.table import Table

from .io import open_data_store


__all__ = [
    "align",
    "composable",
    "dist",
    "evo",
    "io",
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


def available_apps(name_filter: str | None = None) -> Table:
    """returns Table listing the available apps

    Parameters
    ----------
    name_filter
        include apps whose name includes name_filter
    """
    from cogent3.util.table import Table

    from .composable import __app_registry

    # registration of apps does not happen until their modules are imported
    for name in __all__:
        importlib.import_module(f"cogent3.app.{name}")

    # exclude apps from deprecated modules
    deprecated = []

    rows = []
    for app, is_comp in __app_registry.items():
        if any(app.startswith(d) for d in deprecated):
            continue

        if name_filter and name_filter not in app:
            continue

        with contextlib.suppress(AttributeError):
            # probably a local scope issue in testing!
            rows.append(_get_app_attr(app, is_comp))

    header = ["module", "name", "composable", "doc", "input type", "output type"]
    return Table(header=header, data=rows)


_get_param = re.compile('(?<=").+(?=")')


def _make_signature(app: type) -> str:
    from cogent3.util.misc import get_object_provenance

    init_sig = inspect.signature(app.__init__)
    app_name = app.__name__
    params = [f'"{app_name}"']
    empty_default = inspect._empty
    for k, v in init_sig.parameters.items():
        if k == "self":
            continue

        txt = repr(v).replace("\n", " ")
        # clean up text when callable() used  as a type hint
        txt = txt.replace("<built-in function callable>", "callable")

        val = _get_param.findall(txt)[0]
        if v.default is not empty_default and callable(v.default):
            val = val.split("=", maxsplit=1)
            if hasattr(v.default, "app_type"):
                val[-1] = f" {str(v.default)}"
            else:
                val[-1] = f" {get_object_provenance(v.default)}"
            val = "=".join(val)
        params.append(val.replace("\n", " "))

    sig_prefix = f"{app_name}_app = get_app"
    sig_prefix_length = len(sig_prefix)

    if len(", ".join(params)) + sig_prefix_length <= 68:  # plus 2 parentheses makes 70
        params = ", ".join(params)
        return f"{sig_prefix}({params})"

    indent = " " * 4
    params = ",\n".join([f"{indent}{p}" for p in params])

    return f"{sig_prefix}(\n{params},\n)"


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
    modname, _ = table.to_list(columns=["module", "name"])[0]
    mod = importlib.import_module(modname)
    return getattr(mod, name)


def get_app(_app_name: str, *args, **kwargs):
    """returns app instance, use app_help() to display arguments

    Parameters
    ----------
    _app_name
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
    return _get_app_matching_name(_app_name)(*args, **kwargs)


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
        doc.append(line)

    # Remove empty lines at the beginning and end of docstrings
    if doc[0] == "":
        doc.pop(0)
    if doc[-1] == "":
        doc.pop()

    return "\n".join(doc)


def _clean_overview(text: str) -> str:
    text = text.split()
    return "\n".join(textwrap.wrap(" ".join(text), break_long_words=False))


def app_help(name: str):
    """displays help for the named app

    Parameters
    ----------
    name
        app name, e.g. 'min_length', or can include module information,
        e.g. 'cogent3.app.sample.min_length' or 'sample.min_length'. Use the
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
