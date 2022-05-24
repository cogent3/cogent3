"""COmparative GENomics Toolkit 3: providing a first-class genomic sequence
analysis experience within Jupyter notebooks plus supporting parallel
execution on compute systems with 1000s of CPUs."""

import os
import pathlib
import pickle
import sys
import warnings

from typing import Callable, Optional, Union

from cogent3.app import available_apps
from cogent3.core.alignment import (
    Alignment,
    ArrayAlignment,
    Sequence,
    SequenceCollection,
)
from cogent3.core.genetic_code import available_codes, get_code
# note that moltype has to be imported last, because it sets the moltype in
# the objects created by the other modules.
from cogent3.core.moltype import (
    ASCII,
    DNA,
    PROTEIN,
    RNA,
    STANDARD_CODON,
    CodonAlphabet,
    available_moltypes,
    get_moltype,
)
from cogent3.core.tree import PhyloNode, TreeBuilder, TreeError, TreeNode
from cogent3.evolve.fast_distance import (
    available_distances,
    get_distance_calculator,
)
from cogent3.evolve.models import available_models, get_model
from cogent3.parse.cogent3_json import load_from_json
from cogent3.parse.newick import parse_string as newick_parse_string
from cogent3.parse.sequence import FromFilenameParser
from cogent3.parse.table import load_delimited
from cogent3.parse.tree_xml import parse_string as tree_xml_parse_string
from cogent3.util.io import get_format_suffixes, open_
from cogent3.util.table import Table as _Table
from cogent3.util.table import cast_str_to_array


__author__ = ""
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Gavin Huttley",
    "Rob Knight",
    "Peter Maxwell",
    "Jeremy Widmann",
    "Catherine Lozupone",
    "Matthew Wakefield",
    "Edward Lang",
    "Greg Caporaso",
    "Mike Robeson",
    "Micah Hamady",
    "Sandra Smit",
    "Zongzhi Liu",
    "Andrew Butterfield",
    "Amanda Birmingham",
    "Brett Easton",
    "Hua Ying",
    "Jason Carnes",
    "Raymond Sammut",
    "Helen Lindsay",
    "Daniel McDonald",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


_min_version = (3, 7)
if (sys.version_info.major, sys.version_info.minor) < _min_version:
    PY_VERSION = ".".join([str(n) for n in sys.version_info])
    _min_version = ".".join(str(e) for e in _min_version)
    raise RuntimeError(
        f"Python-{_min_version} or greater is required, Python-{PY_VERSION} used."
    )
elif (sys.version_info.major, sys.version_info.minor) == (3, 7):
    warnings.warn(
        "The minimum supported python version will change to 3.8 at 2022.10",
        category=DeprecationWarning,
    )

version = __version__
version_info = tuple([int(v) for v in version.split(".") if v.isdigit()])


warn_env = "COGENT3_WARNINGS"

if warn_env in os.environ:
    warnings.simplefilter(os.environ[warn_env])


def make_seq(seq, name=None, moltype=None):
    """
    Parameters
    ----------
    seq : str
        raw string to be converted to sequence object
    name : str
        sequence name
    moltype
        name of a moltype or moltype instance

    Returns
    -------
    returns a sequence object
    """
    moltype = moltype or "text"
    moltype = get_moltype(moltype)
    seq = moltype.make_seq(seq)
    if name is not None:
        seq.name = name
    return seq


def _make_seq_container(
    klass, data, moltype=None, label_to_name=None, info=None, source=None, **kw
):
    """utility function for creating the different sequence collection/alignment instances"""
    if moltype is not None:
        moltype = get_moltype(moltype)

    info = info or {}
    for other_kw in ("constructor_kw", "kw"):
        other_kw = kw.pop(other_kw, None) or {}
        kw.update(other_kw)
    assert isinstance(info, dict), "info must be a dict"
    source = source or info.get("source", "unknown")
    info["source"] = str(source)

    return klass(
        data=data, moltype=moltype, label_to_name=label_to_name, info=info, **kw
    )


def make_unaligned_seqs(
    data, moltype=None, label_to_name=None, info=None, source=None, **kw
):
    """Initialize an unaligned collection of sequences.

    Parameters
    ----------
    data
        sequences
    moltype
        the moltype, eg DNA, PROTEIN, 'dna', 'protein'
    label_to_name
        function for converting original name into another name.
    info
        a dict from which to make an info object
    source
        origins of this data, defaults to 'unknown'. Converted to a string
        and added to info["source"].
    **kw
        other keyword arguments passed to SequenceCollection
    """

    return _make_seq_container(
        SequenceCollection,
        data,
        moltype=moltype,
        label_to_name=label_to_name,
        info=info,
        source=source,
        **kw,
    )


def make_aligned_seqs(
    data,
    moltype=None,
    array_align=True,
    label_to_name=None,
    info=None,
    source=None,
    **kw,
):
    """Initialize an aligned collection of sequences.

    Parameters
    ----------
    data
        sequences
    moltype
        the moltype, eg DNA, PROTEIN, 'dna', 'protein'
    array_align : bool
        if True, returns ArrayAlignment, otherwise an annotatable Alignment
    label_to_name
        function for converting original name into another name.
    info
        a dict from which to make an info object
    source
        origins of this data, defaults to 'unknown'. Converted to a string
        and added to info["source"].
    **kw
        other keyword arguments passed to alignment class
    """
    klass = ArrayAlignment if array_align else Alignment
    return _make_seq_container(
        klass,
        data,
        moltype=moltype,
        label_to_name=label_to_name,
        info=info,
        source=source,
        **kw,
    )


def _load_seqs(file_format, filename, fmt, kw, parser_kw):
    """utility function for loading sequences"""
    fmt = fmt or file_format
    if not fmt:
        msg = "could not determined file format, set using the format argument"
        raise ValueError(msg)
    parser_kw = parser_kw or {}
    for other_kw in ("constructor_kw", "kw"):
        other_kw = kw.pop(other_kw, None) or {}
        kw.update(other_kw)
    return list(FromFilenameParser(filename, fmt, **parser_kw))


def load_seq(
    filename: Union[str, pathlib.Path],
    format: Optional[str] = None,
    moltype: Optional[str] = None,
    label_to_name: Optional[Callable] = None,
    parser_kw: dict = None,
    info: dict = None,
    **kw,
) -> Sequence:
    """
    loads unaligned sequences from file

    Parameters
    ----------
    filename : str
        path to sequence file
    format : str
        sequence file format, if not specified tries to guess from the path suffix
    moltype : str
        the moltype, eg DNA, PROTEIN, 'dna', 'protein'
    label_to_name : callable
        function for converting original name into another name.
    parser_kw : dict
        optional arguments for the parser
    info : dict
        a dict from which to make an info object
    **kw
        other keyword arguments passed to SequenceCollection

    Notes
    -----
    Returns **one** sequence from a file. Use load_aligned_seqs or
    load_unaligned_seqs to get a collection.

    Returns
    -------
    ``Sequence``
    """
    info = info or {}
    info["source"] = str(filename)
    file_format, _ = get_format_suffixes(filename)
    if file_format == "json":
        seq = load_from_json(filename, (Sequence,))
        seq.name = label_to_name(seq.name) if label_to_name else seq.name
        return seq

    data = _load_seqs(file_format, filename, format, kw, parser_kw)
    name, seq = data[0]
    name = label_to_name(name) if label_to_name else name
    result = make_seq(seq, name, moltype=moltype)
    result.info.update(info)
    return result


def load_unaligned_seqs(
    filename: Union[str, pathlib.Path],
    format=None,
    moltype=None,
    label_to_name=None,
    parser_kw=None,
    info=None,
    **kw,
):
    """
    loads unaligned sequences from file

    Parameters
    ----------
    filename : str
        path to sequence file
    format : str
        sequence file format, if not specified tries to guess from the path suffix
    moltype
        the moltype, eg DNA, PROTEIN, 'dna', 'protein'
    label_to_name
        function for converting original name into another name.
    parser_kw : dict
        optional arguments for the parser
    info
        a dict from which to make an info object
    **kw
        other keyword arguments passed to SequenceCollection

    Returns
    -------
    ``SequenceCollection``
    """
    file_format, _ = get_format_suffixes(filename)
    if file_format == "json":
        return load_from_json(filename, (SequenceCollection,))

    data = _load_seqs(file_format, filename, format, kw, parser_kw)
    return make_unaligned_seqs(
        data,
        label_to_name=label_to_name,
        moltype=moltype,
        source=filename,
        info=info,
        **kw,
    )


def load_aligned_seqs(
    filename: Union[str, pathlib.Path],
    format=None,
    array_align=True,
    moltype=None,
    label_to_name=None,
    parser_kw=None,
    info=None,
    **kw,
):
    """
    loads aligned sequences from file

    Parameters
    ----------
    filename : str
        path to sequence file
    format : str
        sequence file format, if not specified tries to guess from the path suffix
    moltype
        the moltype, eg DNA, PROTEIN, 'dna', 'protein'
    array_align : bool
        if True, returns ArrayAlignment, otherwise an annotatable Alignment
    label_to_name
        function for converting original name into another name.
    parser_kw : dict
        optional arguments for the parser

    Returns
    -------
    ``ArrayAlignment`` or ``Alignment`` instance
    """
    file_format, _ = get_format_suffixes(filename)
    if file_format == "json":
        return load_from_json(filename, (Alignment, ArrayAlignment))

    data = _load_seqs(file_format, filename, format, kw, parser_kw)
    return make_aligned_seqs(
        data,
        array_align=array_align,
        label_to_name=label_to_name,
        moltype=moltype,
        source=filename,
        info=info,
        **kw,
    )


def make_table(
    header=None,
    data=None,
    row_order=None,
    digits=4,
    space=4,
    title="",
    max_width=1e100,
    index_name=None,
    legend="",
    missing_data="",
    column_templates=None,
    data_frame=None,
    format="simple",
    **kwargs,
):
    """

    Parameters
    ----------
    header
        column headings
    data
        a 2D dict, list or tuple. If a dict, it must have column
        headings as top level keys, and common row labels as keys in each
        column.
    row_order
        the order in which rows will be pulled from the twoDdict
    digits
        floating point resolution
    space
        number of spaces between columns or a string
    title
        as implied
    max_width
        maximum column width for printing
    index_name
        column name with values to be used as row identifiers and keys
        for slicing. All column values must be unique.
    legend
        table legend
    missing_data
        replace missing data with this
    column_templates
        dict of column headings
        or a function that will handle the formatting.
    limit
        exits after this many lines. Only applied for non pickled data
        file types.
    data_frame
        a pandas DataFrame, supersedes header/rows
    format
        output format when using str(Table)

    """
    if any(isinstance(a, str) for a in (header, data)):
        raise TypeError("str type invalid, if it's a path use load_table()")

    data = kwargs.get("rows", data)
    if data_frame is not None:
        from pandas import DataFrame

        if not isinstance(data_frame, DataFrame):
            raise TypeError(f"expecting a DataFrame, got{type(data_frame)}")

        data = {c: data_frame[c].to_numpy() for c in data_frame}

    return _Table(
        header=header,
        data=data,
        digits=digits,
        row_order=row_order,
        title=title,
        column_templates=column_templates,
        space=space,
        missing_data=missing_data,
        max_width=max_width,
        index_name=index_name,
        legend=legend,
        data_frame=data_frame,
        format=format,
    )


def load_table(
    filename: Union[str, pathlib.Path],
    sep=None,
    reader=None,
    digits=4,
    space=4,
    title="",
    missing_data="",
    max_width=1e100,
    index_name=None,
    legend="",
    column_templates=None,
    static_column_types=False,
    limit=None,
    format="simple",
    skip_inconsistent=False,
    **kwargs,
):
    """

    Parameters
    ----------
    filename
        path to file containing a tabular data
    sep
        the delimiting character between columns
    reader
        a parser for reading filename. This approach assumes the first
        row returned by the reader will be the header row.
    static_column_types
        if True, and reader is None, identifies columns
        with a numeric/bool data types from the first non-header row.
        This assumes all subsequent entries in that column are of the same type.
        Default is False.
    digits
        floating point resolution
    space
        number of spaces between columns or a string
    title
        as implied
    missing_data
        character assigned if a row has no entry for a column
    max_width
        maximum column width for printing
    index_name
        column name with values to be used as row identifiers and keys
        for slicing. All column values must be unique.
    legend
        table legend
    column_templates
        dict of column headings
        or a function that will handle the formatting.
    limit
        exits after this many lines. Only applied for non pickled data
        file types.
    format
        output format when using str(Table)
    skip_inconsistent
        skips rows that have different length to header row
    """
    if not any(isinstance(filename, t) for t in (str, pathlib.PurePath)):
        raise TypeError(
            "filename must be string or Path, perhaps you want make_table()"
        )

    sep = sep or kwargs.pop("delimiter", None)
    file_format, compress_format = get_format_suffixes(filename)

    if file_format == "json":
        return load_from_json(filename, (_Table,))
    elif file_format in ("pickle", "pkl"):
        with open_(filename, mode="rb") as f:
            loaded_table = pickle.load(f)

        r = _Table()
        r.__setstate__(loaded_table)
        return r

    if reader:
        with open_(filename, newline=None) as f:
            data = list(reader(f))
            header = data[0]
            data = {column[0]: column[1:] for column in zip(*data)}
    else:
        if file_format == "csv":
            sep = sep or ","
        elif file_format == "tsv":
            sep = sep or "\t"

        header, rows, loaded_title, legend = load_delimited(
            filename, sep=sep, limit=limit, **kwargs
        )
        if skip_inconsistent:
            num_fields = len(header)
            rows = [r for r in rows if len(r) == num_fields]
        else:
            lengths = set(map(len, [header] + rows))
            if len(lengths) != 1:
                msg = f"inconsistent number of fields {lengths}"
                raise ValueError(msg)

        title = title or loaded_title
        data = {column[0]: column[1:] for column in zip(header, *rows)}

    for key, value in data.items():
        data[key] = cast_str_to_array(value, static_type=static_column_types)

    return make_table(
        header=header,
        data=data,
        digits=digits,
        title=title,
        column_templates=column_templates,
        space=space,
        missing_data=missing_data,
        max_width=max_width,
        index_name=index_name,
        legend=legend,
        format=format,
    )


def make_tree(treestring=None, tip_names=None, format=None, underscore_unmunge=False):
    """Initialises a tree.

    Parameters
    ----------
    treestring
        a newick or xml formatted tree string
    tip_names
        a list of tip names, returns a "star" topology tree
    format : str
        indicates treestring is either newick or xml formatted, default
        is newick
    underscore_unmunge : bool
        replace underscores with spaces in all names read, i.e. "sp_name"
        becomes "sp name"

    Notes
    -----
    Underscore unmunging is turned off by default, although it is part
    of the Newick format.

    Returns
    -------
    PhyloNode
    """
    assert treestring or tip_names, "must provide either treestring or tip_names"
    if tip_names:
        tree_builder = TreeBuilder().create_edge
        tips = [tree_builder([], tip_name, {}) for tip_name in tip_names]
        tree = tree_builder(tips, "root", {})
        return tree

    if format is None and treestring.startswith("<"):
        format = "xml"
    parser = tree_xml_parse_string if format == "xml" else newick_parse_string
    tree_builder = TreeBuilder().create_edge
    # FIXME: More general strategy for underscore_unmunge
    if parser is newick_parse_string:
        tree = parser(treestring, tree_builder, underscore_unmunge=underscore_unmunge)
    else:
        tree = parser(treestring, tree_builder)
    if not tree.name_loaded:
        tree.name = "root"

    return tree


def load_tree(
    filename: Union[str, pathlib.Path], format=None, underscore_unmunge=False
):
    """Constructor for tree.

    Parameters
    ----------
    filename : str
        a file path containing a newick or xml formatted tree.
    format : str
        either xml or json, all other values default to newick. Overrides
        file name suffix.
    underscore_unmunge : bool
        replace underscores with spaces in all names read, i.e. "sp_name"
        becomes "sp name".

    Notes
    -----
    Underscore unmunging is turned off by default, although it is part
    of the Newick format. Only the cogent3 json and xml tree formats are
    supported.

    Returns
    -------
    PhyloNode
    """
    file_format, _ = get_format_suffixes(filename)
    format = format or file_format
    if format == "json":
        return load_from_json(filename, (TreeNode, PhyloNode))

    with open_(filename) as tfile:
        treestring = tfile.read()

    return make_tree(treestring, format=format, underscore_unmunge=underscore_unmunge)
