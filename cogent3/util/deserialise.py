#!/usr/bin/env python
import os
import json
from importlib import import_module

import cogent3
from cogent3.core.moltype import get_moltype, _CodonAlphabet
from cogent3.core.genetic_code import get_code
from cogent3.util.misc import open_

__author__ = ["Gavin Huttley"]
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


def _get_class(provenance):
    index = provenance.rfind('.')
    assert index > 0
    klass = provenance[index + 1:]
    mod = import_module(provenance[:index])
    klass = getattr(mod, klass)
    return klass


def deserialise_moltype(data):
    """returns a cogent3 MolType instance, or a CodonAlphabet"""
    label = data['moltype']
    data['moltype'] = get_moltype(label)
    klass = _get_class(data.pop('type'))
    if klass == _CodonAlphabet:
        gc = get_code(data.pop('genetic_code'))
        result = _CodonAlphabet(**data)
        result._gc = gc
    else:
        result = data['moltype']

    return result


def deserialise_alphabet(data):
    """returns a cogent3 Alphabet instance"""
    if _get_class(data.get('type')) == _CodonAlphabet:
        result = deserialise_moltype(data)
        return result

    label = data['moltype']
    data['moltype'] = get_moltype(label)
    key = 'data' if 'data' in data else 'motifset'
    motifs = data.pop(key)
    klass = _get_class(data.pop('type'))
    result = klass(motifs, **data)
    return result


def deserialise_seq(data):
    """returns a cogent3 sequence/collection/alignment instance"""
    # We first try to load moltype/alphabet using get_moltype
    from cogent3.core.moltype import get_moltype

    data['moltype'] = get_moltype(data.pop('moltype'))
    make_seq = data['moltype'].make_seq
    type_ = data.pop('type')
    klass = _get_class(type_)

    if 'alignment' not in type_.lower():
        data.pop('moltype')
        result = make_seq(**data)
    else:
        seqs = [make_seq(**v) for v in data.pop('seqs').values()]
        result = klass(seqs, **data)

    return result


def deserialise_tree(data):
    """returns a cogent3 PhyloNode instance"""
    # we load tree using LoadTree, then populate edge attributes
    newick = data.pop('newick')
    edge_attr = data.pop('edge_attributes')
    tree = cogent3.LoadTree(treestring=newick)
    for edge in tree.get_edge_vector():
        params = edge_attr.get(edge.name, {})
        edge.params.update(params)
    return tree


def deserialise_substitution_model(data):
    """returns a cogent3 substitution model instance"""
    from cogent3.evolve.models import get_model
    kw = {} if 'kw' not in data else data.pop('kw')
    sm = None
    if kw and 'name' in kw:
        name = kw.pop('name')
        try:
            sm = get_model(name, **kw)
        except ValueError:  # user defined sm?
            pass

    if sm is None:
        alphabet = deserialise_alphabet(data.pop('alphabet'))
        klass = _get_class(data.pop('type'))
        sm = klass(**data)

    return sm


def deserialise_likelihood_function(data):
    """returns a cogent3 likelihood function instance"""
    model = deserialise_substitution_model(data.pop('model'))
    aln = deserialise_seq(data.pop('alignment'))
    tree = deserialise_tree(data.pop('tree'))
    constructor_args = data.pop('likelihood_construction')
    motif_probs = data.pop('motif_probs')
    param_rules = data.pop('param_rules')
    lf = model.make_likelihood_function(tree, **constructor_args)
    lf.set_alignment(aln)
    with lf.updates_postponed():
        lf.set_motif_probs(motif_probs)
        for rule in param_rules:
            lf.set_param_rule(**rule)
    return lf


def deserialise_object(data):
    if type(data) is str and os.path.exists(data):
        with open_(data) as infile:
            data = json.load(infile)

    if type(data) is str:
        data = json.loads(data)

    type_ = data['type']
    if 'core.sequence' in type_ or 'core.alignment' in type_:
        func = deserialise_seq
    elif 'core.tree' in type_:
        func = deserialise_tree
    elif ('evolve.substitution_model' in type_ or
          'evolve.ns_substitution_model' in type_):
        func = deserialise_substitution_model
    elif 'evolve.parameter_controller' in type_:
        func = deserialise_likelihood_function
    elif 'core.moltype' in type_:
        func = deserialise_moltype
    elif 'core.alphabet' in type_:
        func = deserialise_alphabet
    else:
        msg = "deserialising '%s' from json" % type_
        raise NotImplementedError(msg)
    return func(data)
