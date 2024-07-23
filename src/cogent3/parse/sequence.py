"""Delegator for sequence data format parsers."""

import functools
import pathlib
import typing
import xml.dom.minidom

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
from cogent3.parse.record import FileFormatError
from cogent3.util import warning as c3warn
from cogent3.util.io import get_format_suffixes, iter_splitlines, open_

_lc_to_wc = "".join([[chr(x), "?"]["A" <= chr(x) <= "Z"] for x in range(256)])


@c3warn.deprecated_callable(
    "2024.9",
    reason="allow more customised parser implementations",
    is_discontinued=True,
)
def FromFilenameParser(filename, format=None, **kw):  # pragma: no cover
    """Arguments:
    - filename: name of the sequence alignment file
    - format: the multiple sequence file format
    """
    if format is None:
        format, _ = get_format_suffixes(filename)

    with open_(filename, newline=None, mode="rt") as f:
        data = f.read()

    return FromFileParser(data.splitlines(), format, **kw)


@c3warn.deprecated_callable(
    "2024.9",
    reason="allow more customised parser implementations",
    is_discontinued=True,
)
def FromFileParser(f, format, dialign_recode=False, **kw):  # pragma: no cover
    format = format.lower()
    if format in XML_PARSERS:
        doctype = format
        format = "xml"
    else:
        doctype = None
    if format == "xml":
        source = dom = xml.dom.minidom.parse(f)
        if doctype is None:
            doctype = str(dom.doctype.name).lower()
        if doctype not in XML_PARSERS:
            raise FileFormatError(f"Unsupported XML doctype {doctype}")
        parser = XML_PARSERS[doctype]
    else:
        if format not in PARSERS:
            raise FileFormatError(f"Unsupported file format {format}")
        parser = PARSERS[format]
        source = f
    for name, seq in parser(source, **kw):
        if isinstance(seq, str):
            if dialign_recode:
                seq = seq.translate(_lc_to_wc)
            if not seq.isupper():
                seq = seq.upper()
        yield name, seq


ParserOutputType = typing.Iterable[typing.Tuple[str, str]]


class LineBasedParser:
    """wrapper class to standardise input for line-based sequence format parsers"""

    def __init__(self, parser: typing.Callable[[typing.Any], ParserOutputType]) -> None:
        self._parse = parser

    @functools.singledispatchmethod
    def __call__(self, data, **kwargs) -> ParserOutputType:
        raise TypeError(f"Unsupported data type {type(data)}")

    @__call__.register
    def _(self, data: str, **kwargs) -> ParserOutputType:
        yield from self._parse(iter_splitlines(data), **kwargs)

    @__call__.register
    def _(self, data: pathlib.Path, **kwargs) -> ParserOutputType:
        if not data.exists():
            raise FileNotFoundError(f"File '{data}' does not exist")
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
    "paml": LineBasedParser(paml.PamlParser),
    "fasta": LineBasedParser(fasta.MinimalFastaParser),
    "mfa": LineBasedParser(fasta.MinimalFastaParser),
    "fa": LineBasedParser(fasta.MinimalFastaParser),
    "faa": LineBasedParser(fasta.MinimalFastaParser),
    "fna": LineBasedParser(fasta.MinimalFastaParser),
    "xmfa": LineBasedParser(fasta.MinimalXmfaParser),
    "gde": LineBasedParser(fasta.MinimalGdeParser),
    "aln": LineBasedParser(clustal.ClustalParser),
    "clustal": LineBasedParser(clustal.ClustalParser),
    "gb": LineBasedParser(genbank.RichGenbankParser),
    "gbk": LineBasedParser(genbank.RichGenbankParser),
    "gbff": LineBasedParser(genbank.RichGenbankParser),
    "genbank": LineBasedParser(genbank.RichGenbankParser),
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
        raise ValueError(f"Unsupported format {fmt!r}")
