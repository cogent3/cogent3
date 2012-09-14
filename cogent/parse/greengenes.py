#!/usr/bin/env python

"""Parse the Greengenes formatted sequence data records

The script is intended to be used with the following input:
http://greengenes.lbl.gov/Download/Sequence_Data/Greengenes_format/greengenes16SrRNAgenes.txt.gz
"""

from cogent.parse.record_finder import DelimitedRecordFinder
from cogent.parse.record import DelimitedSplitter, GenericRecord

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Daniel McDonald"] 
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Prototype"

def make_ignore_f(start_line):
    """Make an ignore function that ignores bad gg lines"""
    def ignore(line):
        """Return false if line is bad"""
        return not line or ['',''] == line or [start_line,''] == line
    return ignore

def DefaultDelimitedSplitter(delimiter):
    """Wraps delimited splitter to handle empty records"""
    parser = DelimitedSplitter(delimiter=delimiter)
    def f(line):
        parsed = parser(line)
        if len(parsed) == 1:
            parsed.append('')
        return parsed
    return f

def MinimalGreengenesParser(lines,LineDelim="=",RecStart="BEGIN",RecEnd="END"):
    """Parses raw Greengeens 16S rRNA Gene records
   
    lines  :  open records file
    LineDelim  :  individual line delimiter, eg foo=bar
    RecStart  :  start identifier for a record
    RecEnd  :  end identifier for a record
    """
    line_parser = DefaultDelimitedSplitter(delimiter=LineDelim)

    # parse what the ending record looks like so it can match after being split
    RecordDelim = line_parser(RecEnd)

    # make sure to ignore the starting record
    ignore = make_ignore_f(RecStart)

    parser = DelimitedRecordFinder(RecordDelim, constructor=line_parser, 
                                   keep_delimiter=False, ignore=ignore)

    for record in parser(lines):
        yield GenericRecord(record)

all_ids = lambda x,y: True
specific_ids = lambda x,y: x in y
def SpecificGreengenesParser(lines, fields, ids=None, **kwargs):
    """Yield specific fields from successive Greengenes records
    
    If ids are specified, only the records for the set of ids passed in will
    be returned. Parser will silently ignore ids that are not present in the
    set of ids as well as silently ignore ids in the set that are not present
    in the records file.

    ids : must either test True or be an iterable with prokMSA_ids

    Returns tuples in 'fields' order
    """
    parser = MinimalGreengenesParser(lines, **kwargs)

    if ids:
        ids = set(ids)
        id_lookup = specific_ids
    else:
        id_lookup = all_ids

    for record in parser:
        if id_lookup(record['prokMSA_id'], ids):
            yield tuple([record[field] for field in fields])

def main():
    from optparse import make_option
    from cogent.util.misc import parse_command_line_parameters
    from sys import exit, stdout
    
    script_info = {}
    script_info['brief_description'] = "Parse raw Greengenes 16S records"
    script_info['script_description'] = """Parse out specific fields from raw Greengenes 16S records. These records are rich but often only a subset of each record is required for downstream processing."""
    script_info['script_usage'] = []
    script_info['script_usage'].append(("""Example:""","""Greengenes taxonomy and raw sequences are needed:""","""python greengenes.py -i greengenes16SrRNAgenes.txt -o gg_seq_and_tax.txt -f prokMSA_id,greengenes_tax_string,aligned_seq"""))
    script_info['script_usage'].append(("""Example:""","""Spitting out the available fields from Greengenes:""","""python greengenes.py -i greengenes16SrRNAgenes.txt --print-fields"""))
    script_info['output_description'] = """The resulting output file will contain a header that is prefixed with a # and delimited by the specified delimiter (default is tab). All records will follow in the same order with the same delimiter. It is possible for some key/value pairs within a record to lack a value. In this case, the value placed will be ''"""
    script_info['required_options']=[make_option('--input','-i',dest='input',\
                  help='Greengenes Records')]
    script_info['optional_options']=[\
               make_option('--output','-o',dest='output',help='Output file'),
               make_option('--fields','-f',dest='fields',\
                  help='Greengenes fields to keep'),
               make_option('--delim','-d',dest='delim',help='Output delimiter',\
                       default="\t"),
               make_option('--list-of-ids','-l',dest='ids',default=None,\
                   help='File with a single column list of ids to retrieve'),
               make_option('--print-fields','-p',dest='print_fields',\
                  help='Prints available fields from first Greengenes Record',\
                  action='store_true',default=False)]
    script_info['version'] = __version__
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.print_fields:
        gg_parser = MinimalGreengenesParser(open(opts.input))
        rec = gg_parser.next()
        print '\n'.join(sorted(rec.keys()))
        exit(0)

    if not opts.fields:
        print option_parser.usage()
        print
        print "Greengenes fields must be specified!"
        exit(1)

    if not opts.output:
        output = stdout
    else:
        output = open(opts.output,'w')

    fields = opts.fields.split(',')
    output.write("#%s\n" % opts.delim.join(fields))

    if opts.ids:
        ids = set([l.strip() for l in open(opts.ids, 'U')])
    else:
        ids = None

    gg_parser = SpecificGreengenesParser(open(opts.input), fields, ids)

    for record in gg_parser:
        output.write(opts.delim.join(record))
        output.write('\n')

    if opts.output:
        output.close()
    
if __name__ == '__main__':
    main()
