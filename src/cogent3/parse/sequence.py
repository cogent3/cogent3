"""Delegator for sequence data format parsers."""

import abc
import functools
import pathlib
import typing

from cogent3.parse import (
    clustal,
    fasta,
    gbseq,
    gcg,
    genbank,
    nexus,
    paml,
    phylip,
    tinyseq,
)
from cogent3.util.io import iter_splitlines

ParserOutputType = typing.Iterable[tuple[str, str] | dict]

SeqParserInputTypes = str | pathlib.Path | tuple | list


class SequenceParserBase(abc.ABC):
    """Base class for sequence format parsers."""

    @property
    @abc.abstractmethod
    def name(self) -> str:
        """name of the format"""
        ...

    @property
    @abc.abstractmethod
    def supported_suffixes(self) -> set[str]:
        """Return list of file suffixes this parser supports"""
        ...

    @property
    def result_is_storage(self) -> bool:
        """True if the loader directly returns SeqsDataABC or AlignedSeqdDataABC"""
        return False

    @property
    def supports_unaligned(self) -> bool:
        """True if the loader supports unaligned sequences"""
        return True

    @property
    def supports_aligned(self) -> bool:
        """True if the loader supports aligned sequences"""
        return True

    @property
    @abc.abstractmethod
    def loader(self) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
        """a callable for loading from a file"""
        ...


class LineBasedParser:
    """wrapper class to standardise input for line-based sequence format parsers"""

    def __init__(self, parser: typing.Callable[[typing.Any], ParserOutputType]) -> None:
        self._parse = parser

    @functools.singledispatchmethod
    def __call__(self, data, **kwargs) -> ParserOutputType:
        msg = f"Unsupported data type {type(data)}"
        raise TypeError(msg)

    @__call__.register
    def _(self, data: str, **kwargs) -> ParserOutputType:
        yield from self._parse(iter_splitlines(data), **kwargs)

    @__call__.register
    def _(self, data: pathlib.Path, **kwargs) -> ParserOutputType:
        if not data.exists():
            msg = f"File '{data}' does not exist"
            raise FileNotFoundError(msg)
        yield from self._parse(iter_splitlines(data), **kwargs)

    @__call__.register
    def _(self, data: tuple, **kwargs) -> ParserOutputType:
        # we're assuming this is already split by lines
        yield from self._parse(data, **kwargs)

    @__call__.register
    def _(self, data: list, **kwargs) -> ParserOutputType:
        # we're assuming this is already split by lines
        yield from self._parse(data, **kwargs)


PARSERS = {
    "phylip": LineBasedParser(phylip.MinimalPhylipParser),
    "phy": LineBasedParser(phylip.MinimalPhylipParser),
    "paml": LineBasedParser(paml.PamlParser),
    "fasta": fasta.iter_fasta_records,
    "mfa": fasta.iter_fasta_records,
    "fa": fasta.iter_fasta_records,
    "faa": fasta.iter_fasta_records,
    "fna": fasta.iter_fasta_records,
    "xmfa": LineBasedParser(fasta.MinimalXmfaParser),
    "gde": LineBasedParser(fasta.MinimalGdeParser),
    "aln": LineBasedParser(clustal.ClustalParser),
    "clustal": LineBasedParser(clustal.ClustalParser),
    "gb": genbank.rich_parser,
    "gbk": genbank.rich_parser,
    "gbff": genbank.rich_parser,
    "genbank": genbank.rich_parser,
    "msf": LineBasedParser(gcg.MsfParser),
    "nex": LineBasedParser(nexus.MinimalNexusAlignParser),
    "nxs": LineBasedParser(nexus.MinimalNexusAlignParser),
    "nexus": LineBasedParser(nexus.MinimalNexusAlignParser),
}


class FastaParser(SequenceParserBase):
    """Parser for FASTA format sequence files."""

    @property
    def name(self) -> str:
        return "fasta"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"fasta", "fa", "fna", "faa", "mfa"}

    @property
    def loader(self) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
        return fasta.iter_fasta_records


class GdeParser(SequenceParserBase):
    """Parser for GDE format sequence files."""

    @property
    def name(self) -> str:
        return "gde"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"gde"}

    @property
    def loader(self) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
        return LineBasedParser(fasta.MinimalGdeParser)


class PhylipParser(SequenceParserBase):
    """Parser for PHYLIP format sequence files."""

    @property
    def name(self) -> str:
        return "phylip"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"phylip", "phy"}

    @property
    def loader(self) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
        return LineBasedParser(phylip.MinimalPhylipParser)


class ClustalParser(SequenceParserBase):
    """Parser for Clustal format sequence files."""

    @property
    def name(self) -> str:
        return "clustal"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"clustal", "aln"}

    @property
    def loader(self) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
        return LineBasedParser(clustal.ClustalParser)


class PamlParser(SequenceParserBase):
    """Parser for PAML format sequence files."""

    @property
    def name(self) -> str:
        return "paml"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"paml"}

    @property
    def loader(self) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
        return LineBasedParser(paml.PamlParser)


class NexusParser(SequenceParserBase):
    """Parser for Nexus format sequence files."""

    @property
    def name(self) -> str:
        return "nexus"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"nexus", "nex", "nxs"}

    @property
    def loader(self) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
        return LineBasedParser(nexus.MinimalNexusAlignParser)


class GenbankParser(SequenceParserBase):
    """Parser for GenBank format sequence files."""

    @property
    def name(self) -> str:
        return "genbank"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"gb", "gbk", "gbff", "genbank"}

    @property
    def loader(self) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
        return genbank.rich_parser


class MsfParser(SequenceParserBase):
    """Parser for MSF format sequence files."""

    @property
    def name(self) -> str:
        return "msf"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"msf"}

    @property
    def loader(self) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
        return LineBasedParser(gcg.MsfParser)


class XmfaParser(SequenceParserBase):
    """Parser for XMFA format sequence files."""

    @property
    def name(self) -> str:
        return "xmfa"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"xmfa"}

    @property
    def loader(self) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
        return LineBasedParser(fasta.MinimalXmfaParser)


class TinyseqParser(SequenceParserBase):
    """Parser for Tinyseq format sequence files."""

    @property
    def name(self) -> str:
        return "tinyseq"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"tinyseq"}

    @property
    def loader(self) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
        return LineBasedParser(tinyseq.TinyseqParser)


class GbSeqParser(SequenceParserBase):
    """Parser for GbSeq format sequence files."""

    @property
    def name(self) -> str:
        return "gbseq"

    @property
    def supported_suffixes(self) -> set[str]:
        return {"gbseq"}

    @property
    def loader(self) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
        return gbseq.GbSeqXmlParser


XML_PARSERS = {"gbseq": gbseq.GbSeqXmlParser, "tseq": tinyseq.TinyseqParser}


def get_parser(fmt: str) -> typing.Callable[[SeqParserInputTypes], ParserOutputType]:
    """returns a sequence format parser"""
    try:
        return PARSERS[fmt]
    except KeyError:
        msg = f"Unsupported format {fmt!r}"
        raise ValueError(msg)


def is_genbank(fmt: str) -> bool:
    """whether the provided format is a genbank format"""
    try:
        return get_parser(fmt).__module__.endswith("genbank")
    except ValueError:
        return False
