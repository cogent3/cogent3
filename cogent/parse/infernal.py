#!/usr/bin/env python
#file cogent.parse.infernal.py
"""Parses various Infernal output formats."""

from cogent.parse.table import ConvertFields, SeparatorFormatParser,is_empty

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.4"
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
    conversion_fields = [(1,int),(2,int),(3,int),(4,int),(5,float),(7,int)]
    cmsearch_converter = ConvertFields(conversion_fields)
    #Ignore hash characters and empty lines
    ignore_funct = lambda x: x.startswith('#') or is_empty(x)
    #make parser
    cmsearch_parser = SeparatorFormatParser(with_header=False,\
                                            converter=cmsearch_converter,\
                                            ignore=ignore_funct,\
                                            sep=None)
    
    return cmsearch_parser(lines)

    