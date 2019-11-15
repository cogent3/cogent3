"""The most commonly used constructors are available from this toplevel module.
The rest are in the subpackages: core, draw, evolve, format, maths, parse
and phylo.
"""

import os
import pickle
import re
import sys
import warnings

import numpy

from cogent3.app import available_apps
from cogent3.core.alignment import (
    Alignment,
    ArrayAlignment,
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
from cogent3.core.tree import TreeBuilder, TreeError
from cogent3.evolve.fast_distance import (
    available_distances,
    get_distance_calculator,
)
from cogent3.evolve.models import available_models, get_model
from cogent3.parse.newick import parse_string as newick_parse_string
from cogent3.parse.sequence import FromFilenameParser
from cogent3.parse.table import autogen_reader, load_delimited
from cogent3.parse.tree_xml import parse_string as tree_xml_parse_string
from cogent3.util.misc import get_format_suffixes, open_
from cogent3.util.table import Table as _Table


__author__ = ""
__copyright__ = "Copyright 2007-2019, The Cogent Project"
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
__version__ = "2019.11.15.a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


if sys.version_info < (3, 6):
    PY_VERSION = ".".join([str(n) for n in sys.version_info])
    raise RuntimeError(
        "Python-3.6 or greater is required, Python-%s used." % PY_VERSION
    )

NUMPY_VERSION = re.split(r"[^\d]", numpy.__version__)
numpy_version_info = tuple([int(i) for i in NUMPY_VERSION if i.isdigit()])
if numpy_version_info < (1, 3):
    raise RuntimeError("Numpy-1.3 is required, %s found." % NUMPY_VERSION)

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
        origins of this data, defaults to 'unknown'
    **kw
        other keyword arguments passed to SequenceCollection
    """

    if moltype is not None:
        moltype = get_moltype(moltype)

    info = info or {}
    for other_kw in ("constructor_kw", "kw"):
        other_kw = kw.pop(other_kw, None) or {}
        kw.update(other_kw)
    assert isinstance(info, dict), "info must be a dict"
    info["source"] = source or "unknown"

    return SequenceCollection(
        data=data, moltype=moltype, label_to_name=label_to_name, info=info, **kw
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
        origins of this data, defaults to 'unknown'
    **kw
        other keyword arguments passed to SequenceCollection
    """
    if moltype is not None:
        moltype = get_moltype(moltype)

    info = info or {}
    for other_kw in ("constructor_kw", "kw"):
        other_kw = kw.pop(other_kw, None) or {}
        kw.update(other_kw)
    assert isinstance(info, dict), "info must be a dict"
    info["source"] = source or "unknown"
    klass = ArrayAlignment if array_align else Alignment
    return klass(
        data=data, moltype=moltype, label_to_name=label_to_name, info=info, **kw
    )


def load_unaligned_seqs(
    filename,
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

    Returns
    -------
    ``SequenceCollection``
    """
    parser_kw = parser_kw or {}
    for other_kw in ("constructor_kw", "kw"):
        other_kw = kw.pop(other_kw, None) or {}
        kw.update(other_kw)
    data = list(FromFilenameParser(filename, format, **parser_kw))
    return make_unaligned_seqs(
        data,
        label_to_name=label_to_name,
        moltype=moltype,
        source=filename,
        info=info,
        **kw,
    )


def load_aligned_seqs(
    filename,
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
    parser_kw = parser_kw or {}
    for other_kw in ("constructor_kw", "kw"):
        other_kw = kw.pop(other_kw, None) or {}
        kw.update(other_kw)
    data = list(FromFilenameParser(filename, format, **parser_kw))
    return make_aligned_seqs(
        data,
        array_align=array_align,
        label_to_name=label_to_name,
        moltype=moltype,
        source=filename,
        info=info,
        **kw,
    )


def LoadSeqs(
    filename=None,
    format=None,
    data=None,
    moltype=None,
    name=None,
    aligned=True,
    label_to_name=None,
    parser_kw=None,
    constructor_kw=None,
    array_align=True,
    **kw,
):
    """
    .. deprecated:: 2019.8.30a

        ``LoadSeqs`` will be removed in ``cogent3`` 2020.1.1. It's replaced by
        ``load_unaligned_seqs``, ``load_aligned_seqs``, ``make_unaligned_seqs``
        and ``make_aligned_seqs``.
    """
    kwargs = locals()
    from cogent3.util.warning import deprecated

    if filename and aligned:
        deprecated("function", "LoadSeqs", "load_aligned_seqs", "2020.1.1", 1)
        for key in ("aligned", "data"):
            del kwargs[key]
        return load_aligned_seqs(**kwargs)
    elif filename:
        deprecated("function", "LoadSeqs", "load_unaligned_seqs", "2020.1.1", 1)
        for key in ("aligned", "data", "array_align"):
            del kwargs[key]
        return load_unaligned_seqs(**kwargs)
    elif aligned:
        deprecated("function", "LoadSeqs", "make_aligned_seqs", "2020.1.1", 1)
        for key in ("filename", "format", "aligned", "parser_kw"):
            del kwargs[key]
        return make_aligned_seqs(**kwargs)
    else:
        deprecated("function", "LoadSeqs", "make_unaligned_seqs", "2020.1.1", 1)
        for key in ("filename", "format", "aligned", "array_align", "parser_kw"):
            del kwargs[key]
        return make_unaligned_seqs(**kwargs)


def make_table(
    header=None,
    rows=None,
    row_order=None,
    digits=4,
    space=4,
    title="",
    max_width=1e100,
    row_ids=None,
    legend="",
    missing_data="",
    column_templates=None,
    dtype=None,
    data_frame=None,
    format="simple",
):
    """

    Parameters
    ----------
    header
        column headings
    rows
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
    row_ids
        if True, the 0'th column is used as row identifiers and keys
        for slicing.
    legend
        table legend
    column_templates
        dict of column headings
        or a function that will handle the formatting.
    dtype
        optional numpy array typecode.
    limit
        exits after this many lines. Only applied for non pickled data
        file types.
    data_frame
        a pandas DataFrame, supersedes header/rows
    format
        output format when using str(Table)

    """
    table = _Table(
        header=header,
        rows=rows,
        digits=digits,
        row_order=row_order,
        title=title,
        dtype=dtype,
        column_templates=column_templates,
        space=space,
        missing_data=missing_data,
        max_width=max_width,
        row_ids=row_ids,
        legend=legend,
        data_frame=data_frame,
        format=format,
    )

    return table


def load_table(
    filename,
    sep=None,
    reader=None,
    digits=4,
    space=4,
    title="",
    missing_data="",
    max_width=1e100,
    row_ids=None,
    legend="",
    column_templates=None,
    dtype=None,
    static_column_types=False,
    limit=None,
    format="simple",
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
    header
        column headings
    rows
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
    missing_data
        character assigned if a row has no entry for a column
    max_width
        maximum column width for printing
    row_ids
        if True, the 0'th column is used as row identifiers and keys
        for slicing.
    legend
        table legend
    column_templates
        dict of column headings
        or a function that will handle the formatting.
    dtype
        optional numpy array typecode.
    limit
        exits after this many lines. Only applied for non pickled data
        file types.
    data_frame
        a pandas DataFrame, supersedes header/rows
    format
        output format when using str(Table)

    """
    sep = sep or kwargs.pop("delimiter", None)
    file_format, compress_format = get_format_suffixes(filename)

    if not (reader or static_column_types):
        if file_format == "pickle":
            f = open_(filename, mode="rb")
            loaded_table = pickle.load(f)
            f.close()
            return _Table(**loaded_table)
        elif file_format == "csv":
            sep = sep or ","
        elif file_format == "tsv":
            sep = sep or "\t"

        header, rows, loaded_title, legend = load_delimited(
            filename, delimiter=sep, limit=limit, **kwargs
        )
        title = title or loaded_title
    else:
        f = open_(filename, newline=None)
        if not reader:
            if file_format == "csv":
                sep = sep or ","
            elif file_format == "tsv":
                sep = sep or "\t"
            elif not sep:
                raise ValueError(
                    "static_column_types option requires a value " "for sep"
                )

            reader = autogen_reader(
                f, sep, limit=limit, with_title=kwargs.get("with_title", False)
            )

        rows = [row for row in reader(f)]
        f.close()
        header = rows.pop(0)
    return make_table(
        header=header,
        rows=rows,
        digits=digits,
        title=title,
        dtype=dtype,
        column_templates=column_templates,
        space=space,
        missing_data=missing_data,
        max_width=max_width,
        row_ids=row_ids,
        legend=legend,
        format=format,
    )

    return table


def LoadTable(
    filename=None,
    sep=None,
    reader=None,
    header=None,
    rows=None,
    row_order=None,
    digits=4,
    space=4,
    title="",
    missing_data="",
    max_width=1e100,
    row_ids=None,
    legend="",
    column_templates=None,
    dtype=None,
    static_column_types=False,
    limit=None,
    data_frame=None,
    format="simple",
    **kwargs,
):
    """
    .. deprecated:: 2019.8.30a

        ``LoadTable`` will be removed in ``cogent3`` 2020.1.1. It's replaced by
        ``load_table`` and ``make_table``.
    """
    from cogent3.util.warning import deprecated

    args = {k: v for k, v in locals().items() if k != "deprecated"}
    kwargs = args.pop("kwargs", {})
    args.update(kwargs)
    if filename:
        deprecated("function", "LoadTable", "load_table", "2020.1.1", 1)

        return load_table(**args)

    deprecated("function", "LoadTable", "make_table", "2020.1.1", 1)
    for skip in ("filename", "sep", "reader", "static_column_types", "limit"):
        del args[skip]
    return make_table(**args)


def make_tree(treestring=None, tip_names=None, format=None, underscore_unmunge=False):
    """Initialises a tree.

    Parameters
    ----------
    treestring
        a newick or xml formatted tree string.
    tip_names
        a list of tip names.

    Notes
    -----
    Underscore unmunging is turned off by default, although it is part
    of the Newick format. Set ``underscore_unmunge=True`` to replace underscores
    with spaces in all names read.
    """
    assert treestring or tip_names, "must provide either treestring or tip_names"
    if tip_names:
        tree_builder = TreeBuilder().create_edge
        tips = [tree_builder([], tip_name, {}) for tip_name in tip_names]
        tree = tree_builder(tips, "root", {})
        return tree

    if format is None and treestring.startswith("<"):
        format = "xml"
    if format == "xml":
        parser = tree_xml_parse_string
    else:
        parser = newick_parse_string
    tree_builder = TreeBuilder().create_edge
    # FIXME: More general strategy for underscore_unmunge
    if parser is newick_parse_string:
        tree = parser(treestring, tree_builder, underscore_unmunge=underscore_unmunge)
    else:
        tree = parser(treestring, tree_builder)
    if not tree.name_loaded:
        tree.name = "root"

    return tree


def load_tree(filename, format=None, underscore_unmunge=False):
    """Constructor for tree.

    Parameters
    ----------
    filename
        a file containing a newick or xml formatted tree.

    Notes
    -----
    Underscore unmunging is turned off by default, although it is part
    of the Newick format. Set ``underscore_unmunge=True`` to replace underscores
    with spaces in all names read.
    """

    with open_(filename) as tfile:
        treestring = tfile.read()
        if format is None and filename.endswith(".xml"):
            format = "xml"
    tree = make_tree(treestring, format=format, underscore_unmunge=underscore_unmunge)
    return tree


def LoadTree(
    filename=None,
    treestring=None,
    tip_names=None,
    format=None,
    underscore_unmunge=False,
):
    """
    .. deprecated:: 2019.8.30a

        ``LoadTree`` will be removed in ``cogent3`` 2020.1.1. It's replaced by
        ``load_tree`` and ``make_tree``.
    """
    from cogent3.util.warning import deprecated

    if filename:
        deprecated("function", "LoadTree", "load_tree", "2020.1.1", 1)
        return load_tree(filename, format=format, underscore_unmunge=underscore_unmunge)

    deprecated("function", "LoadTree", "make_tree", "2020.1.1", 1)
    tree = make_tree(
        treestring=treestring,
        tip_names=tip_names,
        format=format,
        underscore_unmunge=underscore_unmunge,
    )
    return tree
