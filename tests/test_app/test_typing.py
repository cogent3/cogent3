import pytest

from cogent3.app.typing import (
    AlignedSeqsType,
    SeqsCollectionType,
    SerialisableType,
    TabularType,
    UnalignedSeqsType,
    get_constraint_names,
    hints_from_strings,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


def test_get_constraint_names():
    """returns the correct names"""
    from cogent3.core.alignment import (
        Alignment,
        ArrayAlignment,
        SequenceCollection,
    )
    from cogent3.evolve.fast_distance import DistanceMatrix
    from cogent3.util.dict_array import DictArray
    from cogent3.util.table import Table

    got = get_constraint_names(AlignedSeqsType)
    assert got == {obj.__name__ for obj in (Alignment, ArrayAlignment)}
    got = get_constraint_names(UnalignedSeqsType)
    assert got == {SequenceCollection.__name__}

    got = get_constraint_names(SeqsCollectionType)
    assert got == {
        obj.__name__ for obj in (Alignment, ArrayAlignment, SequenceCollection)
    }

    got = get_constraint_names(TabularType)

    assert got == {obj.__name__ for obj in (Table, DictArray, DistanceMatrix)}


def test_get_constraint_names_builtins():
    """handles built-in types"""
    from typing import Union

    got = get_constraint_names(Union[str, bytes])
    assert got == {"str", "bytes"}


def test_get_constraint_names_no_constraint():
    """SerialisableType does not define any compatible types"""

    got = get_constraint_names(SerialisableType)
    assert got == set()


def test_hints_from_strings_invalid():
    """raise an exception if unknown string"""
    with pytest.raises(ValueError):
        hints_from_strings("abcde")


def test_hints_from_strings():
    got = hints_from_strings("serialisable", "aligned")
    assert got == [SerialisableType, AlignedSeqsType]
