import json
import re
from importlib import import_module

from cogent3.util.io import open_, path_exists

_deserialise_func_map = {}


class register_deserialiser:
    """
    registration decorator for functions to inflate objects that were
    serialised using json.

    Functions are added to a dict which is used by the deserialise_object()
    function. The type string(s) must uniquely identify the appropriate
    value for the dict 'type' entry, e.g. 'cogent3.core.table.Table'.

    Parameters
    ----------
    args: str or sequence of str
        must be unique
    """

    def __init__(self, *args) -> None:
        for type_str in args:
            if not isinstance(type_str, str):
                msg = f"{type_str!r} is not a string"
                raise TypeError(msg)
            assert type_str not in _deserialise_func_map, (
                f"{type_str!r} already in {list(_deserialise_func_map)}"
            )
        self._type_str = args

    def __call__(self, func):
        for type_str in self._type_str:
            _deserialise_func_map[type_str] = func
        return func


def get_class(provenance: str) -> type:
    index = provenance.rfind(".")
    assert index > 0
    klass = provenance[index + 1 :]
    nc = "NotCompleted"
    klass = nc if nc in klass else klass
    mod = import_module(provenance[:index])
    return getattr(mod, klass)


_pat = re.compile("[a-z]")


def str_to_version(v):
    letter = _pat.search(v)
    return tuple(f"{v[: letter.start()]}.{letter.group()}.{letter.end():}".split("."))


@register_deserialiser("cogent3.evolve.parameter_controller")
def deserialise_likelihood_function(data):
    """returns a cogent3 likelihood function instance"""
    data.pop("version", None)
    model = deserialise_object(data.pop("model"))
    tree = deserialise_object(data.pop("tree"))
    constructor_args = data.pop("likelihood_construction")
    motif_probs = data.pop("motif_probs")
    param_rules = data.pop("param_rules")
    name = data.pop("name", None)
    lf = model.make_likelihood_function(tree, **constructor_args)
    lf.set_name(name)
    lf = model.make_likelihood_function(tree, **constructor_args)

    if isinstance(constructor_args["loci"], list):
        locus_names = constructor_args["loci"]
        align = data["alignment"]
        aln = [deserialise_object(align[k]) for k in locus_names]
        if locus_names[0] in motif_probs:
            mprobs = [motif_probs[k] for k in motif_probs]
        else:
            mprobs = [motif_probs]
    else:
        aln = deserialise_object(data.pop("alignment"))
        mprobs = [motif_probs]

    lf.set_alignment(aln)
    with lf.updates_postponed():
        for motif_probs in mprobs:
            lf.set_motif_probs(motif_probs)
        for rule in param_rules:
            lf.set_param_rule(**rule)
    return lf


def deserialise_object(data):
    """
    deserialises from json

    Parameters
    ----------
    data
        path to json file, json string or a dict

    Returns
    -------
    If the dict from json.loads does not contain a "type" key, the object will
    be returned as is. Otherwise, it will be deserialised to a cogent3 object.

    Notes
    -----
    The value of the "type" key is used to identify the specific function for recreating
    the original instance.
    """
    if path_exists(data):
        with open_(data) as infile:
            data = json.load(infile)

    if isinstance(data, str):
        data = json.loads(str(data))

    type_ = data.get("type", None) if hasattr(data, "get") else None
    if type_ is None:
        return data

    for type_str, func in _deserialise_func_map.items():
        if type_str in type_:
            break
    else:
        msg = f"deserialising '{type_}' from json"
        raise NotImplementedError(msg)

    return func(data)
