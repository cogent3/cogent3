#!/usr/bin/env python
#file cogent.parse.infernal.py
"""Parses various Infernal output formats for the commandline version of:
Infernal 1.0 and 1.0.2 only."""

from cogent.parse.table import ConvertFields, SeparatorFormatParser,is_empty

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
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

    
