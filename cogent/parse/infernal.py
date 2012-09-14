#!/usr/bin/env python
#file cogent.parse.infernal.py
"""Parses various Infernal output formats for the commandline version of:
Infernal 1.0 and 1.0.2 only."""

from cogent.parse.table import ConvertFields, SeparatorFormatParser,is_empty

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

def CmsearchParser(lines):
    """Parser for tabfile format cmsearch result.
    
        - IMPORTANT: Will not parse standard output from cmsearch.  You must
            use --tabfile with cmsearch to get correct format to use this
            parser.
        
        - NOTE: Will only work with search result files with a single CM
            as a query.  Will not work with multiple search result files
            that have been concatenated.
        
        - Result will be list of hits with following order:
        [target name, target start, target stop, query start, query stop,
            bit score, E-value, GC%]
        
    """
    # Converting indices and %GC to integers and bit score to float.
    # Since E-value is only present if CM is calibrated, leaving as string.
    conversion_fields = [(2,int),(3,int),(4,int),(5,int),(6,float),(8,int)]
    cmsearch_converter = ConvertFields(conversion_fields)
    #Ignore hash characters
    good_lines = []
    for l in lines:
        if not l.startswith('#'):
            good_lines.append(l)
    #make parser
    cmsearch_parser = SeparatorFormatParser(with_header=False,\
                                            converter=cmsearch_converter,\
                                            ignore=None,\
                                            sep=None)
    
    return cmsearch_parser(good_lines)

def CmalignScoreParser(lines):
    """Parser for tabfile format cmalign score result.
    
        - IMPORTANT: Will only parse standard output from cmalign.
                
        - NOTE: Will only work with search result files with a single CM
            as a query.  Will not work with multiple alignment result files
            that have been concatenated.
        
        - Result will be list of hits with following order:
        [seq idx, seq name, seq len, total bit score, struct bit score,
            avg prob, elapsed time]
        
    """
    # Converting indices and %GC to integers and bit score to float.
    # Since E-value is only present if CM is calibrated, leaving as string.
    conversion_fields = [(0,int),(2,int),(3,float),(4,float),(5,float)]
    cmalign_score_converter = ConvertFields(conversion_fields)
    #Ignore hash characters
    good_lines = []
    for l in lines:
        line = l.strip()
        if line.startswith('# STOCKHOLM 1.0'):
            break
        if line and (not line.startswith('#')):
            good_lines.append(l)
    #make parser
    cmalign_score_parser = SeparatorFormatParser(with_header=False,\
                                            converter=cmalign_score_converter,\
                                            ignore=None,\
                                            sep=None)
    
    return cmalign_score_parser(good_lines)

    

