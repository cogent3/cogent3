#!/usr/bin/env python
"""Parser for NCBI Sequence Set XML format.
DOCTYPE Bioseq-set PUBLIC "-//NCBI//NCBI Seqset/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_Seqset.dtd"
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

def GbSeqXmlParser(doc):
    """Parser for NCBI Sequence Set XML format.
    DOCTYPE Bioseq-set PUBLIC "-//NCBI//NCBI Seqset/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_Seqset.dtd"
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
    for record in dom_obj.getElementsByTagName('GBSeq'):
        raw_seq = record.getElementsByTagName(
                        'GBSeq_sequence')[0].childNodes[0].nodeValue
        name = record.getElementsByTagName(
                        'GBSeq_accession-version')[0].childNodes[0].nodeValue
        
        #cast as string to de-unicode
        raw_string = str(raw_seq).upper()
        name=str(name)
        
        if record.getElementsByTagName(
                        'GBSeq_moltype')[0].childNodes[0].nodeValue == u'9':
            alphabet = moltype.PROTEIN
        else:
            alphabet = moltype.DNA

        seq = alphabet.makeSequence(raw_string, Name=name)
        
        all = annotation.Map([(0,len(seq))], parent_length=len(seq))
        seq.addAnnotation(annotation.Source, all, name, all)
        
        organism = str(record.getElementsByTagName(
                                'GBSeq_organism')[0].childNodes[0].nodeValue)
        
        seq.addAnnotation(annotation.Feature, "organism", organism, [(0,len(seq))])
        
        features = record.getElementsByTagName('GBFeature')
        for feature in features:
            key = str(feature.getElementsByTagName(
                                    'GBFeature_key')[0].childNodes[0].nodeValue)
            
            if key == 'source':
                continue
            
            spans = []
            feature_name = ""
            
            for interval in feature.getElementsByTagName("GBInterval"):
                try:
                    start = int(interval.getElementsByTagName(
                                    "GBInterval_from")[0].childNodes[0].nodeValue)
                    end= int(interval.getElementsByTagName(
                                    "GBInterval_to")[0].childNodes[0].nodeValue)
                    spans.append((start-1, end))
                except IndexError:
                    point = int(interval.getElementsByTagName(
                                "GBInterval_point")[0].childNodes[0].nodeValue)
                    spans.append((point-1, point))
            if spans == []:
                spans = [(0,len(seq))]
            for qualifier in feature.getElementsByTagName("GBQualifier"):
                qname = qualifier.getElementsByTagName(
                                    "GBQualifier_name")[0].childNodes[0].nodeValue
                if qname == u'gene':
                    feature_name = qualifier.getElementsByTagName(
                                "GBQualifier_value")[0].childNodes[0].nodeValue
            seq.addAnnotation(annotation.Feature, key, feature_name, spans)
        yield (name, seq)

def parse(*args):
    return GbSeqXmlParser(*args).next()[1]
