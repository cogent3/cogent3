#!/usr/bin/env python
import json

from importlib import import_module

import cogent3

from cogent3.core.alignment import Aligned
from cogent3.core.genetic_code import get_code
from cogent3.core.moltype import _CodonAlphabet, get_moltype
from cogent3.util.io import open_, path_exists


__author__ = ["Gavin Huttley"]
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

_deserialise_func_map = {}


class register_deserialiser:
    """
    registration decorator for functions to inflate objects that were
    serialised using json.

    Functions are added to a dict which is used by the deserialise_object()
    function. The type string(s) must uniquely identify the appropriate
    value for the dict 'type' entry, e.g. 'cogent3.util.table.Table'.

    Parameters
    ----------
    args: str or sequence of str
        must be unique
    """

    def __init__(self, *args):
        for type_str in args:
            if not isinstance(type_str, str):
                raise TypeError(f"{type_str!r} is not a string")
            assert (
                type_str not in _deserialise_func_map
            ), f"{type_str!r} already in {list(_deserialise_func_map)}"
        self._type_str = args

    def __call__(self, func):
        for type_str in self._type_str:
            _deserialise_func_map[type_str] = func
        return func


def _get_class(provenance):
    index = provenance.rfind(".")
    assert index > 0
    klass = provenance[index + 1 :]
    nc = "NotCompleted"
    klass = nc if nc in klass else klass
    mod = import_module(provenance[:index])
    klass = getattr(mod, klass)
    return klass


@register_deserialiser(
    "cogent3.util.table.Table",
    "cogent3.util.dict_array.DictArray",
    "cogent3.evolve.fast_distance.DistanceMatrix",
)
def deserialise_tabular(data):
    """deserialising DictArray, Table instances"""
    data.pop("version", None)
    type_ = data.pop("type")
    klass = _get_class(type_)
    if type_.endswith("Table"):
        if "init_table" in data:
            result = klass()
            result.__setstate__(data)
        else:
            result = klass(**data)
    elif "dictarray" in type_.lower():
        named_dims = data.pop("names")
        array = data.pop("array")
        template = klass(*named_dims)
        result = template.wrap(array)
    else:  # DistanceMatrix
        # dists is a list of simple dists from which we reconstruct a 1D dict
        dists = {}
        for element in data["dists"]:
            key = tuple(element[:2])
            value = element[2]
            dists[key] = value
        data["dists"] = dists
        result = klass(**data)

    return result


@register_deserialiser("cogent3.app.composable.NotCompleted")
def deserialise_not_completed(data):
    """deserialising NotCompletedResult"""
    data.pop("version", None)
    klass = _get_class(data.pop("type"))
    init = data.pop("not_completed_construction")
    args = init.pop("args")
    kwargs = init.pop("kwargs")
    return klass(*args, **kwargs)


def deserialise_map_spans(map_element):
    map_element.pop("version", None)
    map_klass = _get_class(map_element.pop("type"))
    spans = []
    for element in map_element["spans"]:
        element.pop("version", None)
        klass = _get_class(element.pop("type"))
        instance = klass(**element)
        spans.append(instance)

    map_element["spans"] = spans
    return map_klass(**map_element)


def deserialise_annotation(data, parent):
    annots = []
    for element in data:
        element.pop("version", None)
        annotations = element.pop("annotations", None)
        klass = _get_class(element.pop("type"))
        kwargs = element.pop("annotation_construction")
        kwargs["map"] = deserialise_map_spans(kwargs["map"])
        instance = klass(parent, **kwargs)
        if annotations:
            deserialise_annotation(annotations, instance)
        annots.append(instance)
    parent.annotations += tuple(annots)


@register_deserialiser("cogent3.app.result")
def deserialise_result(data):
    """returns a result object"""
    data.pop("version", None)
    klass = _get_class(data.pop("type"))
    kwargs = data.pop("result_construction")
    result = klass(**kwargs)
    items = data.pop("items") if "items" in data else data.items()
    for key, value in items:
        # only deserialise the result object, other attributes loaded as
        # required
        if type(value) == dict and "app.result" in str(value.get("type")):
            value = deserialise_object(value)
        try:
            result[key] = value
        except TypeError:
            result[tuple(key)] = value
    return result


@register_deserialiser("cogent3.core.moltype")
def deserialise_moltype(data):
    """returns a cogent3 MolType instance, or a CodonAlphabet"""
    data.pop("version", None)
    label = data["moltype"]
    data["moltype"] = get_moltype(label)
    klass = _get_class(data.pop("type"))
    if klass == _CodonAlphabet:
        gc = get_code(data.pop("genetic_code"))
        result = _CodonAlphabet(**data)
        result._gc = gc
    else:
        result = data["moltype"]

    return result


@register_deserialiser("cogent3.core.alphabet")
def deserialise_alphabet(data):
    """returns a cogent3 Alphabet instance"""
    data.pop("version", None)
    if _get_class(data.get("type")) == _CodonAlphabet:
        result = deserialise_moltype(data)
        return result

    label = data["moltype"]
    data["moltype"] = get_moltype(label)
    key = "data" if "data" in data else "motifset"
    motifs = data.pop(key)
    klass = _get_class(data.pop("type"))
    result = klass(motifs, **data)
    return result


@register_deserialiser("cogent3.core.sequence")
def deserialise_seq(data, aligned=False):
    """deserialises sequence and any annotations

    Parameters
    ----------
    data : dict
        a result of json.loads of a to_rich_dict()
    aligned
        whether sequence type is for an Alignment, in which case an Aligned
        instance will be returned
    Returns
    -------

    """
    from cogent3.core.moltype import get_moltype

    data.pop("version", None)
    data["moltype"] = get_moltype(data.pop("moltype"))
    annotations = data.pop("annotations", None)
    make_seq = data["moltype"].make_seq
    _ = data.pop("type")
    if "-" in data["seq"]:
        aligned = True

    data.pop("moltype")
    result = make_seq(**data)
    if aligned:
        map_, result = result.parse_out_gaps()

    if annotations:
        deserialise_annotation(annotations, result)

    if aligned:
        result = Aligned(map_, result)

    return result


@register_deserialiser("cogent3.core.alignment")
def deserialise_seq_collections(data):
    """returns a cogent3 sequence/collection/alignment instance"""
    # We first try to load moltype/alphabet using get_moltype
    from cogent3.core.moltype import get_moltype

    data.pop("version", None)
    data["moltype"] = get_moltype(data.pop("moltype"))
    annotations = data.pop("annotations", None)
    type_ = data.pop("type")
    klass = _get_class(type_)
    assert "alignment" in type_.lower(), "not alignment type"
    aligned = not type_.endswith("SequenceCollection")
    seqs = []
    for v in data.pop("seqs").values():
        v["moltype"] = data["moltype"]
        seq = deserialise_seq(v, aligned=aligned)
        seqs.append(seq)

    result = klass(seqs, **data)

    if annotations:
        deserialise_annotation(annotations, result)

    return result


@register_deserialiser("cogent3.core.tree")
def deserialise_tree(data):
    """returns a cogent3 PhyloNode instance"""
    data.pop("version", None)
    # we load tree using make_tree, then populate edge attributes
    newick = data.pop("newick")
    edge_attr = data.pop("edge_attributes")
    tree = cogent3.make_tree(treestring=newick)
    for edge in tree.get_edge_vector():
        params = edge_attr.get(edge.name, {})
        edge.params.update(params)
    return tree


@register_deserialiser(
    "cogent3.evolve.substitution_model", "cogent3.evolve.ns_substitution_model"
)
def deserialise_substitution_model(data):
    """returns a cogent3 substitution model instance"""
    from cogent3.evolve.models import get_model

    data.pop("version", None)
    kw = {} if "kw" not in data else data.pop("kw")
    sm = None
    if kw and "name" in kw:
        name = kw.pop("name")
        try:
            sm = get_model(name, **kw)
        except ValueError:  # user defined sm?
            pass

    if sm is None:
        alphabet = deserialise_alphabet(data.pop("alphabet"))
        klass = _get_class(data.pop("type"))
        sm = klass(alphabet, **data)

    return sm


@register_deserialiser("cogent3.evolve.parameter_controller")
def deserialise_likelihood_function(data):
    """returns a cogent3 likelihood function instance"""
    data.pop("version", None)
    model = deserialise_substitution_model(data.pop("model"))
    tree = deserialise_tree(data.pop("tree"))
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
        aln = [deserialise_seq_collections(align[k]) for k in locus_names]
        if locus_names[0] in motif_probs:
            mprobs = [motif_probs[k] for k in motif_probs]
        else:
            mprobs = [motif_probs]
    else:
        aln = deserialise_seq_collections(data.pop("alignment"))
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

    if type(data) is str:
        data = json.loads(data)

    type_ = data.get("type", None)
    if type_ is None:
        return data

    for type_str, func in _deserialise_func_map.items():
        if type_str in type_:
            break
    else:
        msg = f"deserialising '{type_}' from json"
        raise NotImplementedError(msg)

    return func(data)
