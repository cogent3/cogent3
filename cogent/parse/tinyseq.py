#!/usr/bin/env python
"""Parser for NCBI Tiny Seq XML format.
DOCTYPE TSeqSet PUBLIC "-//NCBI//NCBI TSeq/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd"
"""
import xml.dom.minidom
from cogent.core import annotation, moltype

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Matthew Wakefield", "Peter Maxwell", "Gavin Huttley",
                    "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production"

"""
CAUTION:
This XML PARSER uses minidom. This means a bad performance for 
big files (>5MB), and huge XML files will for sure crash the program!
"""

def TinyseqParser(doc):
    """Parser for NCBI Tiny Seq XML format.
    DOCTYPE TSeqSet PUBLIC "-//NCBI//NCBI TSeq/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd"
    Arguments:
        - doc: An xml.dom.minidom.Document, file object of string
    Yields:
        - name, cogent sequence
    
    CAUTION:
    This XML PARSER uses minidom. This means a bad performance for 
    big files (>5MB), and huge XML files will for sure crash the program!
    """
    if isinstance(doc,xml.dom.minidom.Document):
        dom_obj = doc
    elif isinstance(doc,file):
        dom_obj = xml.dom.minidom.parse(doc)
    elif isinstance(doc,str):
        dom_obj = xml.dom.minidom.parseString(doc)
    else:
        raise TypeError
    for record in dom_obj.getElementsByTagName('TSeq'):
        raw_seq = record.getElementsByTagName(
                        'TSeq_sequence')[0].childNodes[0].nodeValue
        name = record.getElementsByTagName(
                        'TSeq_accver')[0].childNodes[0].nodeValue
        
        #cast as string to de-unicode
        raw_string = str(raw_seq).upper()
        name=str(name)
        
        if record.getElementsByTagName(
                        'TSeq_seqtype')[0].getAttribute('value') == u'protein':
            alphabet = moltype.PROTEIN
        else:
            alphabet = moltype.DNA
        
        seq = alphabet.makeSequence(raw_string, Name=name)
        
        seq.addAnnotation(annotation.Feature, "genbank_id", name, [(0,len(seq))])
        
        organism = str(record.getElementsByTagName(
                                    'TSeq_orgname')[0].childNodes[0].nodeValue)
        
        seq.addAnnotation(annotation.Feature, "organism", organism, [(0,len(seq))])
        
        yield (name, seq)


def parse(*args):
    return TinyseqParser(*args).next()[1]
