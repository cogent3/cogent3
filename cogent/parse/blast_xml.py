#!/usr/bin/env python
"""Parsers for XML output of blast, psi-blast and blat.
"""

__author__ = "Kristian Rother"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Micah Hamady"]
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Prototype"

import xml.dom.minidom

"""
CAUTION:
This XML BLAST PARSER uses minidom. This means a bad performance for 
big files (>5MB), and huge XML files will for sure crash the program!
(06/2009 Kristian)

Possible improvements:
- convert some values into floats automatically (feature request)
- MH recommends sax.* for faster processing.
- test against nt result
- test really big file.
- consider high speed parser for standard output
"""

from cogent.parse.blast import BlastResult

# field names used to parse tags and create dict.
HIT_XML_FIELDNAMES = ['QUERY ID','SUBJECT_ID','HIT_DEF','HIT_ACCESSION',\
                      'HIT_LENGTH']

HSP_XML_FIELDS = (
    ('PERCENT_IDENTITY','Hsp_identity'),
    ('ALIGNMENT_LENGTH','Hsp_align-len'),
    ('MISMATCHES',''),
    ('GAP_OPENINGS','Hsp_gaps'),
    ('QUERY_START','Hsp_query-from'),
    ('QUERY_END','Hsp_query-to'),
    ('SUBJECT_START','Hsp_hit-from'),
    ('SUBJECT_END','Hsp_hit-to'),
    ('E_VALUE','Hsp_evalue'),
    ('BIT_SCORE','Hsp_bit-score'),
    ('SCORE','Hsp_score'),
    ('POSITIVE','Hsp_positive'),
    ('QUERY_ALIGN','Hsp_qseq'),
    ('SUBJECT_ALIGN','Hsp_hseq'),
    ('MIDLINE_ALIGN','Hsp_midline'),
    )

HSP_XML_FIELDNAMES = [x[0] for x in HSP_XML_FIELDS]
HSP_XML_TAGNAMES = [x[1] for x in HSP_XML_FIELDS]
     

def get_tag(record, name, default=None):
    """
    Loks in the XML tag 'record' for 
    other tags named 'name', and returns the value of the first one.
    If none is found, it returns 'default'.
    """
    tag = record.getElementsByTagName(name)
    if len(tag) and len(tag[0].childNodes):
        return tag[0].childNodes[0].nodeValue
    else:
        return default

def parse_hit(hit_tag,query_id=1):
    """
    Parses a 'Hit' dom object.
    Returns a list of lists with HSP data.
    """
    result = []
    # parse elements from hit tag
    hit_id = get_tag(hit_tag,'Hit_id')
    hit_def = get_tag(hit_tag,'Hit_def')
    accession = get_tag(hit_tag,'Hit_accession')
    length = int(get_tag(hit_tag,'Hit_len'), 0)
    hit_data = [query_id,hit_id, hit_def, accession, length]
    # process HSPS in this hit.
    for hsp_tag in hit_tag.getElementsByTagName('Hsp'):
        result.append(hit_data + parse_hsp(hsp_tag))
    return result

def parse_hsp(hsp_tag):
    """
    Parses a 'Hsp' XML dom object. Returns a list of values,
    according to the items in HSP_XML_FIELDS.
    """
    result = []
    for tag_name in HSP_XML_TAGNAMES:
        result.append(get_tag(hsp_tag,tag_name,0))
    # what about these?
    #    self.identity = int(self.get_tag(record,'Hsp_identity', 0))
    #    self.positive = int(self.get_tag(record, 'Hsp_positive', 0))
    return result


def parse_header(tag):
    """
    Parses a 'BlastOutput' dom object.
    Returns a dict with information from the blast header
    """
    result = {}
    result['application'] = get_tag(tag,'BlastOutput_program')
    result['version'] = get_tag(tag,'BlastOutput_version')
    result['reference'] = get_tag(tag,'BlastOutput_reference')
    result['query'] = get_tag(tag,'BlastOutput_query-def')
    result['query_letters'] = int(get_tag(tag,'BlastOutput_query-len'))
    result['database'] = get_tag(tag,'BlastOutput_db')
    # add data fro Parameters tag
    for param_tag in tag.getElementsByTagName('BlastOutput_param'):
        #for param_tag in tag.getElementsByTagName('Parameters'):
        data = parse_parameters(param_tag)
        for k in data:
            result[k] = data[k]
    return result

def parse_parameters(tag):
    """Parses a 'BlastOutput_param' dom object."""
    result = {}
    result['matrix'] = get_tag(tag,'Parameters_matrix')
    result['expect'] = get_tag(tag,'Parameters_expect')
    result['gap_open_penalty'] = float(get_tag(tag,'Parameters_gap-open'))
    result['gap_extend_penalty'] = float(get_tag(tag,'Parameters_gap-extend'))
    result['filter'] = get_tag(tag,'Parameters_filter')
    return result
    
def MinimalBlastParser7(lines, include_column_names=False, format='xml'):
    """Yields succesive records from lines (props, data list).

    lines must be XML BLAST output format.

    output:
    props is a dict of {UPPERCASE_KEY:value}.
    data_list is a list of list of strings, optionally with header first.

    LIST CONTAINS [HIT][HSP][strings], FIRST ENTRY IS LIST OF LABELS!
    """
    doc = ''.join(lines)
    dom_obj = xml.dom.minidom.parseString(doc)
    query_id = 1
    for record in dom_obj.getElementsByTagName('BlastOutput'):
        props = parse_header(record)
        hits = [HIT_XML_FIELDNAMES + HSP_XML_FIELDNAMES]
        for hit in record.getElementsByTagName('Hit'):
            hits += parse_hit(hit,query_id)
        yield props,hits


class BlastXMLResult(BlastResult):
    """the BlastResult objects have the query sequence as keys,
    and the values are lists of lists of dictionaries.
    The FIELD NAMES given are the keys of the dict.
    """
    # FIELD NAMES
    QUERY_ALIGN = 'HSP QSEQ'
    SUBJECT_ALIGN = 'HSP HSEQ'
    MIDLINE_ALIGN = 'HSP MIDLINE'
    HIT_DEF = 'HIT_DEF'
    HIT_ACCESSION = 'HIT_ACCESSION'
    HIT_LENGTH = 'HIT_LENGTH'
    SCORE = 'SCORE'
    POSITIVE = 'POSITIVE'

    #FieldComparisonOperators = (
    #    BlastResult.FieldComparisonOperators = {
    #        HIT_DEF:(_gt, float)
    #        }
    # .. to be done

    # .. extend HitKeys
    HitKeys = BlastResult.HitKeys.union(
                                  set([ HIT_DEF,
                                        HIT_ACCESSION,
                                        HIT_LENGTH,
                                        SCORE,
                                        POSITIVE,
                                        QUERY_ALIGN,
                                        SUBJECT_ALIGN,
                                        MIDLINE_ALIGN ]))

    def __init__(self, data, psiblast=False, parser=None, xml=False):
        # iterate blast results, generate data structure
        """
        Init using blast 7 or blast 9 results

        data: blast output from the m = 9 output option
        psiblast: if True, will expect psiblast output, else expects 
            blast output

        """
        # further improvement:
        # add XML option to BlastResult __init__ instead of
        # using a separate class.

        if not parser:
            if xml:
                parser = MinimalBlastParser7
            elif psiblast:
                parser = MinimalPsiBlastParser9
            else:
                parser = MinimalBlastParser9

        # code below copied from BlastResult, unchanged.
        mp = parser(data, True)
                
        for props, rec_data in mp:
         
            iteration = 1
            if self.ITERATION in props:
                iteration = int(props[self.ITERATION])

            hits = []
            # check if found any hits
            if len(rec_data) > 1:
                for h in rec_data[1:]:
                    hits.append(dict(zip(rec_data[0], h)))
            else:
                hits.append(dict(zip(rec_data[0], ['' for x in rec_data[0]])))
            
            # get blast version of query id
            query_id = hits[0][self.QUERY_ID]

            if query_id not in self: 
                self[query_id] = [] 
            self[query_id].append(hits)
        

