#!/usr/bin/env python
"""Classes for reading multiple sequence alignment files in different formats."""

import xml.dom.minidom

from cogent.parse import fasta, phylip, paml, clustal, genbank
from cogent.parse import gbseq, tinyseq, macsim, gcg
from cogent.parse.record import FileFormatError

__author__ = "Cath Lawrence"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Cath Lawrence", "Gavin Huttley", "Peter Maxwell",
                    "Matthew Wakefield", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

_lc_to_wc = ''.join([[chr(x),'?']['A' <= chr(x) <= 'Z'] for x in range(256)])

def FromFilenameParser(filename, format=None, **kw):
    """Arguments:
            - filename: name of the sequence alignment file
            - format: the multiple sequence file format
    """
    format = format_from_filename(filename, format)
    f = open(filename, 'U')
    return FromFileParser(f, format, **kw)

def FromFileParser(f, format, dialign_recode=False, **kw):
    format = format.lower()
    if format in XML_PARSERS:
        doctype = format
        format = 'xml'
    else:
        doctype = None
    if format == 'xml':
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
        if isinstance(seq, basestring):
            if dialign_recode:
                seq = seq.translate(_lc_to_wc)
            if not seq.isupper():
                seq = seq.upper()
        yield (name, seq)

def format_from_filename(filename, format=None):
    """Detects format based on filename."""
    if format:
        return format
    else:
        return filename[filename.rfind('.')+1:]

PARSERS  = {
        'phylip': phylip.MinimalPhylipParser,
        'paml':  paml.PamlParser,
        'fasta': fasta.MinimalFastaParser,
        'mfa': fasta.MinimalFastaParser,
        'fa': fasta.MinimalFastaParser,
        'faa': fasta.MinimalFastaParser,
        'fna': fasta.MinimalFastaParser,
        'xmfa': fasta.MinimalXmfaParser,
        'gde': fasta.MinimalGdeParser,
        'aln': clustal.ClustalParser,
        'clustal': clustal.ClustalParser,
        'gb': genbank.RichGenbankParser,
        'gbk': genbank.RichGenbankParser,
        'genbank': genbank.RichGenbankParser,
        'msf': gcg.MsfParser,
        }

XML_PARSERS = {
        'gbseq': gbseq.GbSeqXmlParser,
        'tseq': tinyseq.TinyseqParser,
        'macsim': macsim.MacsimParser,
        }

