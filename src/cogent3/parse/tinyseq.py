"""Parser for NCBI Tiny Seq XML format.
DOCTYPE TSeqSet PUBLIC "-//NCBI//NCBI TSeq/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd"
"""

import io
import xml.dom.minidom

from cogent3.core import moltype

"""
CAUTION:
This XML PARSER uses minidom. This means a bad performance for
big files (>5MB), and huge XML files will for sure crash the program!
"""


def TinyseqParser(doc):
    """Parser for NCBI Tiny Seq XML format.
    DOCTYPE TSeqSet PUBLIC "-//NCBI//NCBI TSeq/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd"

    Parameters
    ----------
    doc
        An xml.dom.minidom.Document, file object of string

    Returns
    -------
    name, cogent sequence

    CAUTION:
    This XML PARSER uses minidom. This means a bad performance for
    big files (>5MB), and huge XML files will for sure crash the program!
    """
    if isinstance(doc, xml.dom.minidom.Document):
        dom_obj = doc
    elif isinstance(doc, io.IOBase):
        dom_obj = xml.dom.minidom.parse(doc)
    elif isinstance(doc, str):
        dom_obj = xml.dom.minidom.parseString(doc)
    else:
        raise TypeError
    for record in dom_obj.getElementsByTagName("TSeq"):
        raw_seq = (
            record.getElementsByTagName("TSeq_sequence")[0].childNodes[0].nodeValue
        )
        name = record.getElementsByTagName("TSeq_accver")[0].childNodes[0].nodeValue

        # cast as string to de-unicode
        raw_string = str(raw_seq).upper()
        name = str(name)

        if (
            record.getElementsByTagName("TSeq_seqtype")[0].getAttribute("value")
            == "protein"
        ):
            alphabet = moltype.PROTEIN
        else:
            alphabet = moltype.DNA

        seq = alphabet.make_seq(raw_string, name=name)

        seq.add_feature(biotype="genbank_id", name=name, spans=[(0, len(seq))])

        organism = str(
            record.getElementsByTagName("TSeq_orgname")[0].childNodes[0].nodeValue
        )

        seq.add_feature(biotype="organism", name=organism, spans=[(0, len(seq))])

        yield (name, seq)
