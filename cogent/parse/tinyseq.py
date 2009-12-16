#!/usr/bin/env python
#from xml.dom.minidom import parse as xmlparse
from cogent.core import annotation, moltype

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Matthew Wakefield", "Peter Maxwell", "Gavin Huttley",
                    "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@anu.edu.au"
__status__ = "Production"

def TinyseqParser(doc):
    for record in doc.getElementsByTagName('TSeq'):
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
