#!/usr/bin/env python
"""Classes for reading multiple sequence alignment files in different formats."""

import json
import re
import xml.dom.minidom

from cogent3.app.data_store import load_record_from_json
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
from cogent3.util.deserialise import deserialise_object
from cogent3.util.misc import get_object_provenance, open_


__author__ = "Cath Lawrence"
__copyright__ = "Copyright 2007-2020, The Cogent Project"
__credits__ = [
    "Cath Lawrence",
    "Gavin Huttley",
    "Peter Maxwell",
    "Matthew Wakefield",
    "Rob Knight",
]
__license__ = "BSD-3"
__version__ = "2020.7.2a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

_lc_to_wc = "".join([[chr(x), "?"]["A" <= chr(x) <= "Z"] for x in range(256)])
_compression = re.compile(r"\.(gz|bz2)$")


def FromFilenameParser(filename, format=None, **kw):
    """Arguments:
            - filename: name of the sequence alignment file
            - format: the multiple sequence file format
    """
    format = format_from_filename(filename, format)
    f = open_(filename, newline=None, mode="rt")
    return FromFileParser(f, format, **kw)


def FromFileParser(f, format, dialign_recode=False, **kw):
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
            raise FileFormatError("Unsupported XML doctype %s" % doctype)
        parser = XML_PARSERS[doctype]
    else:
        if format not in PARSERS:
            raise FileFormatError("Unsupported file format %s" % format)
        parser = PARSERS[format]
        source = f
    for (name, seq) in parser(source, **kw):
        if isinstance(seq, str):
            if dialign_recode:
                seq = seq.translate(_lc_to_wc)
            if not seq.isupper():
                seq = seq.upper()
        yield (name, seq)

    f.close()


def format_from_filename(filename, format=None):
    """Detects format based on filename."""
    if format:
        return format
    else:
        r = _compression.search(filename)
        if r:
            filename = filename[: r.start()]

        return filename[filename.rfind(".") + 1 :]


def load_from_json(filename, classes):
    """Loads objects from json files.

    Parameters
    ----------
    filename: name of the json file
    classes: A series of the Cogent3 types, for example: (Alignment, ArrayAlignment)

    """
    assert all(
        (isinstance(klass, type) for klass in classes)
    ), "classes should be a series of Cogent3 types, for example: (Alignment, ArrayAlignment)"

    with open_(filename) as f:
        content = json.loads(f.read())
    try:
        _, data, completed = load_record_from_json(content)
        if not completed:
            raise TypeError("json file is a record for type NotCompleted.")
    except KeyError:
        data = content

    type_ = data.get("type", None)
    if type_ is None:
        raise TypeError("json does not contain 'type' key")

    valid_types = {get_object_provenance(klass) for klass in classes}
    if type_ not in valid_types:
        raise TypeError(f"Invalid data type: {type_} is not one of {valid_types}")

    return deserialise_object(data)


PARSERS = {
    "phylip": phylip.MinimalPhylipParser,
    "paml": paml.PamlParser,
    "fasta": fasta.MinimalFastaParser,
    "mfa": fasta.MinimalFastaParser,
    "fa": fasta.MinimalFastaParser,
    "faa": fasta.MinimalFastaParser,
    "fna": fasta.MinimalFastaParser,
    "xmfa": fasta.MinimalXmfaParser,
    "gde": fasta.MinimalGdeParser,
    "aln": clustal.ClustalParser,
    "clustal": clustal.ClustalParser,
    "gb": genbank.RichGenbankParser,
    "gbk": genbank.RichGenbankParser,
    "genbank": genbank.RichGenbankParser,
    "msf": gcg.MsfParser,
    "nex": nexus.MinimalNexusAlignParser,
    "nxs": nexus.MinimalNexusAlignParser,
}

XML_PARSERS = {"gbseq": gbseq.GbSeqXmlParser, "tseq": tinyseq.TinyseqParser}
