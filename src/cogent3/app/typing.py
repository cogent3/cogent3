"""defined type hints for app composability"""
# todo write more extensive docstring explaining limited use of these types
from __future__ import annotations

from typing import ForwardRef, Iterable, TypeVar, Union

from typing_extensions import get_args, get_origin


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

AlignedSeqsType = TypeVar("AlignedSeqsType", "Alignment", "ArrayAlignment")
UnalignedSeqsType = TypeVar("UnalignedSeqsType", bound="SequenceCollection")
SeqsCollectionType = Union[AlignedSeqsType, UnalignedSeqsType]
SeqType = TypeVar("SeqType", bound="Sequence")
PairwiseDistanceType = TypeVar("PairwiseDistanceType", bound="DistanceMatrix")
TabularType = TypeVar("TabularType", "Table", "DictArray", "DistanceMatrix")
TreeType = TypeVar("TreeType", "TreeNode", "PhyloNode")
SerialisableType = TypeVar("SerialisableType")
BootstrapResultType = TypeVar("BootstrapResultType", bound="bootstrap_result")
HypothesisResultType = TypeVar("HypothesisResultType", bound="hypothesis_result")
ModelCollectionResultType = TypeVar(
    "ModelCollectionResultType", bound="model_collection_result"
)
ModelResultType = TypeVar("ModelResultType", bound="model_result")
TabularResultType = TypeVar("TabularResultType", bound="tabular_result")
GenericResultType = TypeVar("GenericResultType", bound="generic_result")
ResultType = Union[
    GenericResultType,
    BootstrapResultType,
    HypothesisResultType,
    ModelResultType,
    TabularResultType,
]

# todo when move to python 3.8 define protocols for IdentifierType and SerialisableType
IdentifierType = TypeVar("IdentifierType")

# the following constants are deprecated
ALIGNED_TYPE = "aligned"
IDENTIFIER_TYPE = "identifier"
PAIRWISE_DISTANCE_TYPE = "pairwise_distances"
SEQUENCE_TYPE = "sequences"
SERIALISABLE_TYPE = "serialisable"
TABULAR_TYPE = "tabular"
TREE_TYPE = "tree"
BOOTSTRAP_RESULT_TYPE = "bootstrap_result"
HYPOTHESIS_RESULT_TYPE = "hypothesis_result"
MODEL_RESULT_TYPE = "model_result"
RESULT_TYPE = "result"
TABULAR_RESULT_TYPE = "tabular_result"

_mappings = {
    TABULAR_TYPE: TabularType,
    ALIGNED_TYPE: AlignedSeqsType,
    SEQUENCE_TYPE: SeqsCollectionType,
    IDENTIFIER_TYPE: IdentifierType,
    PAIRWISE_DISTANCE_TYPE: PairwiseDistanceType,
    SERIALISABLE_TYPE: SerialisableType,
    TREE_TYPE: TreeType,
    BOOTSTRAP_RESULT_TYPE: BootstrapResultType,
    HYPOTHESIS_RESULT_TYPE: HypothesisResultType,
}


def get_constraint_names(*hints) -> set[str, ...]:
    """returns the set of named constraints of a type hint"""
    all_hints = set()
    for hint in hints:
        if hint in (SerialisableType, IdentifierType):
            all_hints.add(hint.__name__)
            continue

        if getattr(hint, "__bound__", None):
            all_hints.add(hint.__bound__)
            continue

        if getattr(hint, "__constraints__", None):
            all_hints.update(hint.__constraints__)
            continue

        if get_origin(hint) == Union:
            all_hints.update(get_constraint_names(*get_args(hint)))

        if type(hint) == type:
            all_hints.add(hint.__name__)

    all_hints = {h.__forward_arg__ if type(h) == ForwardRef else h for h in all_hints}
    return all_hints


def hints_from_strings(*strings: Iterable[str]) -> list:
    """returns list of type hints corresponding to string values"""
    types = []
    for string in strings:
        string = string.lower()
        if string not in _mappings:
            raise ValueError(f"{string!r} not a known type constant")
        types.append(_mappings[string])
    return types
