"""COmparative GENomics Toolkit 3: providing a first-class genomic sequence
analysis experience within Jupyter notebooks plus supporting parallel
execution on compute systems with 1000s of CPUs."""

import os
import pathlib
import pickle
import warnings
from typing import Callable, Optional, Union

from cogent3._version import __version__
from cogent3.app import app_help, available_apps, get_app, open_data_store  # noqa
from cogent3.core import annotation_db as _anno_db
from cogent3.core.alignment import (
    Alignment,
    ArrayAlignment,
    Sequence,
    SequenceCollection,
)
from cogent3.core.genetic_code import available_codes, get_code  # noqa

# note that moltype has to be imported last, because it sets the moltype in
# the objects created by the other modules.
from cogent3.core.moltype import (
    ASCII,  # noqa
    DNA,  # noqa
    PROTEIN,  # noqa
    RNA,  # noqa
    available_moltypes,  # noqa
    get_moltype,  # noqa
)
from cogent3.core.tree import PhyloNode, TreeBuilder, TreeError, TreeNode
from cogent3.evolve.fast_distance import (
    available_distances,  # noqa
    get_distance_calculator,  # noqa
)
from cogent3.evolve.models import available_models, get_model  # noqa
from cogent3.parse.cogent3_json import load_from_json
from cogent3.parse.newick import parse_string as newick_parse_string
from cogent3.parse.sequence import get_parser, is_genbank
from cogent3.parse.table import load_delimited
from cogent3.parse.tree_xml import parse_string as tree_xml_parse_string
from cogent3.util import warning as c3warn
from cogent3.util.io import get_format_suffixes, is_url, open_
from cogent3.util.progress_display import display_wrap
from cogent3.util.table import Table as _Table
from cogent3.util.table import cast_str_to_array

__copyright__ = "Copyright 2007-2023, The Cogent Project"
__credits__ = "https://github.com/cogent3/cogent3/graphs/contributors"
__license__ = "BSD-3"


version = __version__
version_info = tuple([int(v) for v in version.split(".") if v.isdigit()])


warn_env = "COGENT3_WARNINGS"

if warn_env in os.environ:
    warnings.simplefilter(os.environ[warn_env])


import logging

# suppress numba warnings
__numba_logger = logging.getLogger("numba")
__numba_logger.setLevel(logging.WARNING)

load_annotations = _anno_db.load_annotations


def make_seq(
    seq,
    name: str = None,
    moltype=None,
    new_type: bool = False,
    annotation_offset: int = 0,
    annotation_db: Optional[_anno_db.SupportsFeatures] = None,
    **kw: dict,
):  # refactor: type hinting, need to capture optional args and the return type
    """
    Parameters
    ----------
    seq
        raw string to be converted to sequence object
    name
        sequence name
    moltype
        name of a moltype or moltype instance
    new_type
        if True, returns a new type Sequence (cogent3.core.new_sequence.Sequence).
        The default will be changed to True in 2024.12. Support for the old
        style will be removed as of 2025.6.
    annotation_offset
        integer indicating start position relative to annotations
    **kw
        other keyword arguments passed to Sequence

    Returns
    -------
    returns a sequence object
    """
    moltype = moltype or "text"
    if new_type or "COGENT3_NEW_TYPE" in os.environ:
        from cogent3.core import new_moltype

        moltype = new_moltype.get_moltype(moltype)
    else:
        moltype = get_moltype(moltype)
    seq = moltype.make_seq(
        seq=seq,
        name=name,
        annotation_offset=annotation_offset,
        **kw,
    )
    if annotation_db:
        seq.annotation_db = annotation_db
    return seq


def _make_seq_container(
    klass,
    data,
    moltype=None,
    label_to_name=None,
    info=None,
    source=None,
    **kw,
):
    """utility function for creating the different sequence collection/alignment instances"""
    if moltype is not None:
        moltype = get_moltype(moltype)

    info = info or {}
    for other_kw in ("constructor_kw", "kw"):
        other_kw = kw.pop(other_kw, None) or {}
        kw |= other_kw
    assert isinstance(info, dict), "info must be a dict"
    source = source or info.get("source", "unknown")
    info["source"] = str(source)

    return klass(
        data=data,
        moltype=moltype,
        label_to_name=label_to_name,
        info=info,
        **kw,
    )


def make_unaligned_seqs(
    data,
    moltype=None,
    label_to_name=None,
    info=None,
    source=None,
    new_type=False,
    **kw,
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
    new_type
        if True, the returned SequenceCollection will be of the new type,
        (cogent3.core.new_sequence.SequenceCollection). The default will be
        changed to True in 2024.12. Support for the old style will be removed
        as of 2025.6.
    **kw
        other keyword arguments passed to SequenceCollection
    """

    if new_type or "COGENT3_NEW_TYPE" in os.environ:
        if moltype is None:
            raise ValueError("Argument 'moltype' is required when 'new_type=True'")

        from cogent3.core import new_alignment

        return new_alignment.make_unaligned_seqs(
            data,
            moltype=moltype,
            label_to_name=label_to_name,
            info=info,
            source=source,
            **kw,
        )
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
    new_type=False,
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
    new_type
        if True, the returned Alignment will be of the new type,
        (cogent3.core.new_sequence.Alignment). The default will be
        changed to True in 2025.6. Support for the old style will be removed
        as of 2025.12.
    **kw
        other keyword arguments passed to alignment class
    """
    if new_type or "COGENT3_NEW_TYPE" in os.environ:
        if moltype is None:
            raise ValueError("Argument 'moltype' is required when 'new_type=True'")

        from cogent3.core import new_alignment

        return new_alignment.make_aligned_seqs(
            data,
            moltype=moltype,
            label_to_name=label_to_name,
            info=info,
            source=source,
            **kw,
        )

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


def _load_files_to_unaligned_seqs(
    *,
    path: os.PathLike,
    format: Optional[str] = None,
    moltype: Optional[str] = None,
    label_to_name: Optional[Callable] = None,
    parser_kw: Optional[dict] = None,
    info: Optional[dict] = None,
    new_type: bool = False,
    ui=None,
) -> SequenceCollection:
    """loads multiple files and returns as a sequence collection"""

    file_names = list(path.parent.glob(path.name))
    seqs = [
        load_seq(
            fn,
            format=format,
            moltype=moltype,
            label_to_name=label_to_name,
            parser_kw=parser_kw,
        )
        for fn in ui.series(file_names)
    ]

    return make_unaligned_seqs(
        seqs,
        label_to_name=label_to_name,
        moltype=moltype,
        source=path,
        info=info,
        new_type=new_type,
    )


def _load_seqs(file_format, filename, fmt, parser_kw):
    """utility function for loading sequences"""
    if not is_url(filename):
        filename = pathlib.Path(filename).expanduser()
    fmt = fmt or file_format
    if not fmt:
        msg = "could not determined file format, set using the format argument"
        raise ValueError(msg)
    parser_kw = parser_kw or {}
    parser = get_parser(fmt)
    return list(parser(filename, **parser_kw))


T = Optional[_anno_db.SupportsFeatures]


def _load_genbank_seq(
    filename: os.PathLike,
    parser_kw: dict,
    just_seq: bool = False,
) -> tuple[str, str, T]:
    """utility function for loading sequences"""
    from cogent3.parse.genbank import iter_genbank_records

    for name, seq, features in iter_genbank_records(filename, **parser_kw):
        break
    else:
        raise ValueError(f"No sequences found in {filename}")

    db = (
        None
        if just_seq
        else _anno_db.GenbankAnnotationDb(
            data=features.pop("features", None),
            seqid=name,
        )
    )
    return name, seq, db


def load_seq(
    filename: os.PathLike,
    annotation_path: Optional[os.PathLike] = None,
    format: Optional[str] = None,
    moltype: Optional[str] = None,
    label_to_name: Optional[Callable] = None,
    parser_kw: Optional[dict] = None,
    info: Optional[dict] = None,
    new_type: bool = False,
    annotation_offset: int = 0,
    **kw: dict,
) -> Sequence:
    """
    loads unaligned sequences from file

    Parameters
    ----------
    filename
        path to sequence file
    annotation_path
        path to annotation file, ignored if format is genbank
    format
        sequence file format, if not specified tries to guess from the path suffix
    moltype
        the moltype, eg DNA, PROTEIN, 'dna', 'protein'
    label_to_name
        function for converting original name into another name.
    parser_kw
        optional arguments for the parser
    info
        a dict from which to make an info object
    new_type
        if True, returns a new type Sequence (cogent3.core.new_sequence.Sequence)
        The default will be changed to True in 2024.12. Support for the old
        style will be removed as of 2025.6.
    annotation_offset
        integer indicating start position relative to annotations
    **kw
        other keyword arguments passed to sequence loader

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
        seq = load_from_json(filename, (Sequence,))  # need to support new seq here
        seq.name = label_to_name(seq.name) if label_to_name else seq.name
        return seq

    if is_genbank(format or file_format):
        name, seq, db = _load_genbank_seq(
            filename,
            parser_kw or {},
            just_seq=annotation_path is not None,
        )
    else:
        db = None
        data = _load_seqs(file_format, filename, format, parser_kw)
        name, seq = data[0]

    if annotation_path is not None:
        db = load_annotations(path=annotation_path, seqids=[name])

    name = label_to_name(name) if label_to_name else name

    result = make_seq(
        seq,
        name,
        moltype=moltype,
        new_type=new_type,
        annotation_offset=annotation_offset,
        annotation_db=db,
        **kw,
    )
    result.info.update(info)

    return result


@display_wrap
def load_unaligned_seqs(
    filename: Union[str, pathlib.Path],
    format=None,
    moltype=None,
    label_to_name=None,
    parser_kw: Optional[dict] = None,
    info: Optional[dict] = None,
    new_type: bool = False,
    **kw,
) -> SequenceCollection:
    """
    loads unaligned sequences from file

    Parameters
    ----------
    filename
        path to sequence file or glob pattern. If a glob we assume a single
        sequence per file. All seqs returned in one SequenceCollection.
    format
        sequence file format, if not specified tries to guess from the path suffix
    moltype
        the moltype, eg DNA, PROTEIN, 'dna', 'protein'
    label_to_name
        function for converting original name into another name.
    parser_kw
        optional arguments for the parser
    info
        a dict from which to make an info object
    new_type
        if True, the returned SequenceCollection will be of the new type,
        (cogent3.core.new_sequence.SequenceCollection). The default will be
        changed to True in 2024.12. Support for the old style will be removed
        as of 2025.6.
    **kw
        other keyword arguments passed to SequenceCollection, or show_progress.
        The latter induces a progress bar for number of files processed when
        filename is a glob pattern.

    Returns
    -------
    ``SequenceCollection``
    """
    ui = kw.pop("ui")
    filename = pathlib.Path(filename)
    file_format, _ = get_format_suffixes(filename)

    if "*" in filename.name:
        return _load_files_to_unaligned_seqs(
            path=filename,
            format=file_format,
            moltype=moltype,
            label_to_name=label_to_name,
            parser_kw=parser_kw,
            info=info,
            new_type=new_type,
            ui=ui,
        )

    if file_format == "json":
        return load_from_json(filename, (SequenceCollection,))

    data = _load_seqs(file_format, filename, format, parser_kw)

    return make_unaligned_seqs(
        data,
        label_to_name=label_to_name,
        moltype=moltype,
        source=filename,
        info=info,
        new_type=new_type,
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
    new_type: bool = False,
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
    new_type
        if True, the returned Alignment will be of the new type,
        (cogent3.core.new_alignment.Alignment). The default will be
        changed to True in 2024.12. Support for the old style will be removed
        as of 2025.6.
    kw
        passed to make_aligned_seqs

    Returns
    -------
    ``ArrayAlignment`` or ``Alignment`` instance
    """
    file_format, _ = get_format_suffixes(filename)
    if file_format == "json":
        return load_from_json(filename, (Alignment, ArrayAlignment))

    data = _load_seqs(file_format, filename, format, parser_kw)
    return make_aligned_seqs(
        data,
        array_align=array_align,
        label_to_name=label_to_name,
        moltype=moltype,
        source=filename,
        info=info,
        new_type=new_type,
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
            "filename must be string or Path, perhaps you want make_table()",
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
            filename,
            sep=sep,
            limit=limit,
            **kwargs,
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


@c3warn.deprecated_args(
    version="2024.12",
    reason="argument has no effect",
    discontinued="name_nodes",
)
def make_tree(
    treestring=None,
    tip_names=None,
    format=None,
    underscore_unmunge=False,
):
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
        tips = [tree_builder([], str(tip_name), {}) for tip_name in tip_names]
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
    filename: Union[str, pathlib.Path],
    format=None,
    underscore_unmunge=False,
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
