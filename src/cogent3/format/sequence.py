import abc
import contextlib
import os
import pathlib
import typing

from cogent3.format import clustal, fasta, gde, paml, phylip
from cogent3.parse.record import FileFormatError
from cogent3.util.io import atomic_write

if typing.TYPE_CHECKING:
    from cogent3.core.new_alignment import Alignment, SequenceCollection


SeqsTypes = typing.Union["SequenceCollection", "Alignment"]

FORMATTERS = {
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
    @abc.abstractmethod
    def supported_suffixes(self) -> set[str]:
        """Return list of file suffixes this parser supports"""
        ...

    def formatted(self, seqcoll: SeqsTypes, **kwargs) -> str:
        """returns a string representation of the sequence collection

        Parameters
        ----------
        seqcoll
            sequence collection to format, must have a to_dict() method
        """
        formatter = FORMATTERS[self.name]
        return formatter(seqcoll.to_dict(), **kwargs)

    def write(
        self,
        *,
        path: pathlib.Path,
        seqcoll: SeqsTypes,
        **kwargs,
    ) -> pathlib.Path:
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
        output = self.formatted(seqcoll, **kwargs)
        with atomic_write(path, mode="wt") as f:
            f.write(output)
        return path


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
