import abc
import contextlib
import os
import pathlib
import typing
import warnings
from collections.abc import Callable, Mapping

from cogent3.format import clustal, fasta, gde, paml, phylip
from cogent3.parse.record import FileFormatError

SeqsTypes = "CollectionBase[typing.Any] | Sequence"

FORMATTERS: dict[str, Callable[..., str]] = {
    "phylip": phylip.alignment_to_phylip,
    "paml": paml.alignment_to_paml,
    "fasta": fasta.seqs_to_fasta,
    "gde": gde.alignment_to_gde,
    "clustal": clustal.clustal_from_alignment,
}


class SequenceWriterBase(abc.ABC):
    """Base class for sequence format parsers."""

    @property
    @abc.abstractmethod
    def name(self) -> str:
        """name of the format"""
        ...

    @property
    def supports_unaligned(self) -> bool:
        """True if the writer supports unaligned sequences"""
        return True

    @property
    def supports_aligned(self) -> bool:
        """True if the writer supports aligned sequences"""
        return True

    @property
    def saves_metadata(self) -> bool:
        """True if the writer saves collection metadata"""
        return False

    @property
    @abc.abstractmethod
    def supported_suffixes(self) -> set[str]:
        """Return list of file suffixes this parser supports"""
        ...

    def formatted(self, seqcoll: SeqsTypes, **kwargs: typing.Any) -> str:
        """returns a string representation of the sequence collection

        Parameters
        ----------
        seqcoll
            sequence, or a sequence collection, to format, must have a to_dict() method
        """
        formatter = FORMATTERS[self.name]
        return formatter(seqcoll.to_dict(), **kwargs)

    def write(
        self,
        *,
        path: pathlib.Path | str,
        seqcoll: SeqsTypes,
        **kwargs: typing.Any,
    ) -> pathlib.Path | str:
        """returns a path after writing the sequence collection to a file
        Parameters
        ----------
        path
            path to the file to write
        seqcoll
            sequence collection to write, must have a to_dict() method
        kwargs
            additional arguments to pass to the formatter
        """
        from scinexus.io_util import atomic_write

        output = self.formatted(seqcoll, **kwargs)
        with atomic_write(path, mode="wt") as f:
            f.write(output)
        return path


def warn_if_metadata_unsaved(
    metadata: Mapping[str, typing.Any], writer: typing.Any, *, stacklevel: int
) -> None:
    """Emit a UserWarning if writer does not save the non-empty metadata.

    A writer without a ``saves_metadata`` attribute is treated as not saving
    metadata. ``stacklevel`` is passed to ``warnings.warn``.
    """
    if not metadata or getattr(writer, "saves_metadata", False):
        return
    msg = (
        f"the {writer.name!r} format does not save metadata, "
        f"so these keys are not written: {sorted(metadata)}"
    )
    warnings.warn(msg, UserWarning, stacklevel=stacklevel)


class FastaWriter(SequenceWriterBase):
    @property
    def name(self) -> str:
        return "fasta"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"fasta", "fa", "fna", "faa", "mfa"}


class GdeWriter(SequenceWriterBase):
    @property
    def name(self) -> str:
        return "gde"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"gde"}


class PhylipWriter(SequenceWriterBase):
    """Parser for PHYLIP format sequence files."""

    @property
    def name(self) -> str:
        return "phylip"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"phylip", "phy"}


class PamlWriter(SequenceWriterBase):
    """Parser for PAML format sequence files."""

    @property
    def name(self) -> str:
        return "paml"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"paml"}


def save_to_filename(alignment, filename, format, **kw) -> None:
    """Arguments:
    - alignment: to be written
    - filename: name of the sequence alignment file
    - format: the multiple sequence file format
    """
    from scinexus.io_util import atomic_write

    if format is None:
        msg = "format not known"
        raise FileFormatError(msg)

    with atomic_write(filename, mode="wt") as f:
        try:
            write_alignment_to_file(f, alignment, format, **kw)
        except Exception:
            with contextlib.suppress(Exception):
                os.unlink(filename)
            raise


def write_alignment_to_file(f, alignment, format, **kw) -> None:
    format = format.lower()
    if format not in FORMATTERS:
        msg = f"Unsupported file format {format}"
        raise FileFormatError(msg)
    contents = FORMATTERS[format](alignment, **kw)
    f.write(contents)
    f.close()
