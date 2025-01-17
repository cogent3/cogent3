"""Delegator for sequence data format parsers."""

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

_lc_to_wc = "".join([[chr(x), "?"]["A" <= chr(x) <= "Z"] for x in range(256)])


ParserOutputType = typing.Iterable[tuple[str, str]]


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

XML_PARSERS = {"gbseq": gbseq.GbSeqXmlParser, "tseq": tinyseq.TinyseqParser}

SeqParserInputTypes = typing.Union[str, pathlib.Path, tuple, list]


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
