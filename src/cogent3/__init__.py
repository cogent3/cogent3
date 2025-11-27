"""COmparative GENomics Toolkit 3: providing a first-class genomic sequence
analysis experience within Jupyter notebooks plus supporting parallel
execution on compute systems with 1000s of CPUs."""

import logging
import os
import pathlib
import typing
import warnings
from collections.abc import Callable
from importlib import import_module

from cogent3._version import __version__

if typing.TYPE_CHECKING:  # pragma: no cover
    import numpy
    import numpy.typing as npt

    from cogent3.core.alignment import Alignment, SequenceCollection
    from cogent3.core.annotation_db import SupportsFeatures
    from cogent3.core.moltype import MolTypeLiteral
    from cogent3.core.sequence import Sequence

__copyright__ = "Copyright 2007-date, The Cogent Project"
__credits__ = "https://github.com/cogent3/cogent3/graphs/contributors"
__license__ = "BSD-3"


def __getattr__(name: str) -> typing.Any:  # noqa: ANN401
    if (attr := globals().get(name)) is not None:
        return attr

    if name not in _import_mapping:
        try:
            attr = __import__(name)
        except ImportError as err:
            raise AttributeError(name) from err

        globals()[name] = attr
        return attr

    module_name = _import_mapping[name]
    module = import_module(f".{module_name}", package=__name__)
    attr = getattr(module, name)
    globals()[name] = attr
    return attr


_import_mapping = {
    "available_distances": "evolve.fast_distance",
    "get_distance_calculator": "evolve.fast_distance",
    "available_datasets": "_dataset",
    "get_dataset": "_dataset",
    "set_storage_defaults": "_plugin",
    "make_aligned_seqs": "core.alignment",
    "make_unaligned_seqs": "core.alignment",
    "load_annotations": "core.annotation_db",
    "load_tree": "core.tree",
    "make_tree": "core.tree",
    "PhyloNode": "core.tree",
    "TreeError": "core.tree",
    "load_table": "core.table",
    "make_table": "core.table",
    "load_delimited": "parse.table",
    "ASCII": "core.moltype",
    "DNA": "core.moltype",
    "RNA": "core.moltype",
    "PROTEIN": "core.moltype",
    "available_moltypes": "core.moltype",
    "get_moltype": "core.moltype",
    "MolTypeLiteral": "core.moltype",
    "open_": "util.io",
    "available_models": "evolve.models",
    "get_model": "evolve.models",
    "available_codes": "core.genetic_code",
    "get_code": "core.genetic_code",
    "app_help": "app",
    "available_apps": "app",
    "get_app": "app",
    "open_data_store": "app.io",
}


def __dir__() -> list[str]:
    return list(_import_mapping.keys()) + list(globals().keys())


__all__ = list(_import_mapping.keys())

version = __version__
version_info = tuple(int(v) for v in version.split(".") if v.isdigit())


warn_env = "COGENT3_WARNINGS"

if warn := os.environ.get(warn_env):
    warnings.simplefilter(warn)


# suppress numba warnings
__numba_logger = logging.getLogger("numba")
__numba_logger.setLevel(logging.WARNING)


def make_seq(
    seq,
    name: str | None = None,
    moltype: "MolTypeLiteral | None" = None,
    annotation_offset: int = 0,
    annotation_db: "SupportsFeatures | None" = None,
    **kwargs: typing.Any,  # noqa: ANN401
):  # refactor: type hinting, need to capture optional args and the return type
    """
    Parameters
    ----------
    seq
        raw string to be converted to sequence object
    name
        sequence name
    moltype
        name of a moltype or moltype instance. If None, defaults to 'text'.
    annotation_offset
        integer indicating start position relative to annotations
    **kwargs
        other keyword arguments passed to Sequence

    Returns
    -------
    returns a sequence object
    """
    from cogent3.core.moltype import get_moltype

    mtype = get_moltype(moltype)

    seq = mtype.make_seq(
        seq=seq,
        name=name,
        annotation_offset=annotation_offset,
        **kwargs,
    )
    if annotation_db:
        seq.annotation_db = annotation_db
    return seq


def _load_files_to_unaligned_seqs(
    *,
    path: os.PathLike,
    format_name: str | None = None,
    moltype: "MolTypeLiteral | None" = None,
    label_to_name: Callable | None = None,
    parser_kw: dict | None = None,
    info: dict | None = None,
    **kwargs: typing.Any,  # noqa: ANN401
) -> "SequenceCollection":
    """loads multiple files and returns as a sequence collection"""
    from cogent3.core.alignment import make_unaligned_seqs

    ui = kwargs.pop("ui")
    file_names = list(path.parent.glob(path.name))
    seqs = [
        load_seq(
            fn,
            format_name=format_name,
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
    )


def _load_genbank_seq(
    filename: os.PathLike,
    parser_kw: dict,
    just_seq: bool = False,
) -> tuple[str, "str | bytes | npt.NDArray[numpy.integer]", "SupportsFeatures | None"]:
    """utility function for loading sequences"""
    from cogent3.core.annotation_db import GenbankAnnotationDb
    from cogent3.parse.genbank import iter_genbank_records

    for name, seq, features in iter_genbank_records(filename, **parser_kw):
        break
    else:
        msg = f"No sequences found in {filename}"
        raise ValueError(msg)

    db = (
        None
        if just_seq
        else GenbankAnnotationDb(
            data=features.pop("features", None),
            seqid=name,
        )
    )
    return name, seq, db


def load_seq(
    filename: os.PathLike | str,
    annotation_path: os.PathLike | str | None = None,
    format_name: str | None = None,
    moltype: "MolTypeLiteral | None" = None,
    label_to_name: Callable | None = None,
    parser_kw: dict | None = None,
    info: dict | None = None,
    annotation_offset: int = 0,
    **kwargs: typing.Any,  # noqa: ANN401
) -> "Sequence":
    """
    loads unaligned sequences from file

    Parameters
    ----------
    filename
        path to sequence file
    annotation_path
        path to annotation file, ignored if format is genbank
    format_name
        sequence file format, if not specified tries to guess from the path suffix
    moltype
        the moltype, eg DNA, PROTEIN, 'dna', 'protein'
    label_to_name
        function for converting original name into another name.
    parser_kw
        optional arguments for the parser
    info
        a dict from which to make an info object
    annotation_offset
        integer indicating start position relative to annotations
    **kwargs
        other keyword arguments passed to sequence loader

    Notes
    -----
    Returns **one** sequence from a file. Use load_aligned_seqs or
    load_unaligned_seqs to get a collection.

    Returns
    -------
    ``Sequence``
    """
    from cogent3._plugin import get_seq_format_parser_plugin
    from cogent3.core.annotation_db import load_annotations
    from cogent3.parse.cogent3_json import load_from_json
    from cogent3.parse.sequence import is_genbank
    from cogent3.util.io import get_format_suffixes, is_url

    if not is_url(filename):
        filename = pathlib.Path(filename).expanduser()

    info = info or {}
    info["source"] = str(filename)
    file_suffix, _ = get_format_suffixes(filename)
    parser_kw = parser_kw or {}
    if file_suffix == "json":
        from cogent3.core.sequence import Sequence
        from cogent3.core.sequence import Sequence as OldSeq

        seq = load_from_json(filename, (Sequence, OldSeq))
        seq.name = label_to_name(seq.name) if label_to_name else seq.name
        return seq

    if is_genbank(format_name or file_suffix):
        name, seq, db = _load_genbank_seq(
            pathlib.Path(filename),
            parser_kw,
            just_seq=annotation_path is not None,
        )
    else:
        db = None
        parser = get_seq_format_parser_plugin(
            format_name=format_name,
            file_suffix=file_suffix,
            unaligned_seqs=True,
        )
        if parser.result_is_storage:
            msg = (
                "Cannot reliably derive single sequence from multi-sequence storage. "
                "Use either load_unaligned_seqs() or load_aligned_seqs()."
            )
            raise ValueError(msg)

        data = list(parser.loader(filename, **parser_kw))
        name, seq = data[0]

    name = label_to_name(name) if label_to_name else name

    if annotation_path is not None:
        db = load_annotations(path=annotation_path, seqids=[name])

    result = make_seq(
        seq,
        name,
        moltype=moltype,
        annotation_offset=annotation_offset,
        annotation_db=db,
        **kwargs,
    )
    result.info.update(info)

    return result


def load_unaligned_seqs(
    filename: str | pathlib.Path,
    format_name: str | None = None,
    moltype: "MolTypeLiteral | None" = None,
    label_to_name: typing.Callable[[str], str] | None = None,
    parser_kw: dict | None = None,
    info: dict | None = None,
    **kwargs: typing.Any,  # noqa: ANN401
) -> "SequenceCollection":
    """
    loads unaligned sequences from file

    Parameters
    ----------
    filename
        path to sequence file or glob pattern. If a glob we assume a single
        sequence per file. All seqs returned in one SequenceCollection.
    format_name
        sequence file format, if not specified tries to guess from the path suffix
    moltype
        the moltype, eg DNA, PROTEIN, 'dna', 'protein'
    label_to_name
        function for converting original name into another name.
    parser_kw
        optional arguments for the parser
    info
        a dict from which to make an info object
    **kwargs
        other keyword arguments passed to SequenceCollection, or show_progress.
        The latter induces a progress bar for number of files processed when
        filename is a glob pattern.

    Returns
    -------
    ``SequenceCollection``
    """
    from cogent3._plugin import get_seq_format_parser_plugin
    from cogent3.core.alignment import SequenceCollection, make_unaligned_seqs
    from cogent3.parse.cogent3_json import load_from_json
    from cogent3.util.io import get_format_suffixes, is_url

    if not is_url(filename):
        filename = pathlib.Path(filename).expanduser()

    file_suffix, _ = get_format_suffixes(filename)
    if "*" in filename.name:
        from cogent3.util.progress_display import display_wrap

        func = display_wrap(_load_files_to_unaligned_seqs)
        return func(
            path=filename,
            format_name=format_name or file_suffix,
            moltype=moltype,
            label_to_name=label_to_name,
            parser_kw=parser_kw,
            info=info,
            **kwargs,
        )

    if file_suffix == "json":
        return load_from_json(filename, (SequenceCollection,))

    if not (file_suffix or format_name):
        msg = "could not determined file format, set using the format argument"
        raise ValueError(msg)

    parser = get_seq_format_parser_plugin(
        format_name=format_name,
        file_suffix=file_suffix,
        unaligned_seqs=True,
    )
    parser_kw = parser_kw or {}
    if parser.result_is_storage:
        data = parser.loader(path=filename, **parser_kw)
    else:
        data = list(parser.loader(filename, **parser_kw))

    return make_unaligned_seqs(
        data,
        label_to_name=label_to_name,
        moltype=moltype,
        source=filename,
        info=info,
        **kwargs,
    )


def load_aligned_seqs(
    filename: str | pathlib.Path,
    format_name: str | None = None,
    moltype: "MolTypeLiteral | None" = None,
    label_to_name: typing.Callable[[str], str] | None = None,
    parser_kw: dict | None = None,
    info: dict | None = None,
    **kwargs: typing.Any,  # noqa: ANN401
) -> "Alignment":
    """
    loads aligned sequences from file

    Parameters
    ----------
    filename
        path to sequence file
    format_name
        sequence file format, if not specified tries to guess from the path suffix
    moltype
        the moltype, eg DNA, PROTEIN, 'dna', 'protein'
    label_to_name
        function for converting original name into another name.
    parser_kw
        optional arguments for the parser
    kwargs
        passed to make_aligned_seqs

    Returns
    -------
    ``Alignment`` instance
    """
    from cogent3._plugin import get_seq_format_parser_plugin
    from cogent3.core.alignment import Alignment, make_aligned_seqs
    from cogent3.parse.cogent3_json import load_from_json
    from cogent3.util.io import get_format_suffixes, is_url

    if not is_url(filename):
        filename = pathlib.Path(filename).expanduser()

    file_suffix, _ = get_format_suffixes(filename)
    if file_suffix == "json":
        return load_from_json(filename, (Alignment,))

    parser = get_seq_format_parser_plugin(
        format_name=format_name,
        file_suffix=file_suffix,
        unaligned_seqs=False,
    )
    parser_kw = parser_kw or {}
    if parser.result_is_storage:
        data = parser.loader(path=filename, **parser_kw)
    else:
        data = list(parser.loader(filename, **parser_kw))

    return make_aligned_seqs(
        data,
        label_to_name=label_to_name,
        moltype=moltype,
        source=filename,
        info=info,
        **kwargs,
    )
