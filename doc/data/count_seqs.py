#!/usr/bin/env python
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The PyCogent project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from glob import glob
from cogent.util.option_parsing import (
 parse_command_line_parameters, 
 make_option)
from cogent.parse.fasta import MinimalFastaParser

script_info = {}
script_info['brief_description'] = "Count sequences in one or more fasta files."
script_info['script_description'] = "This script counts the number of sequences in one or more fasta files and prints the results to stdout."
script_info['script_usage'] = [\
 ("Count sequences in one file",
  "Count the sequences in a fasta file and write results to stdout.",
  "%prog -i in.fasta"),
 ("Count sequences in two file",
  "Count the sequences in two fasta files and write results to stdout.",
  "%prog -i in1.fasta,in2.fasta"),
  ("Count the sequences in many fasta files",
   "Count the sequences all .fasta files in current directory and write results to stdout. Note that -i option must be quoted.",
   "%prog -i \"*.fasta\"")]
script_info['output_description']= "Tabular data is written to stdout."
script_info['required_options'] = [
 make_option('-i','--input_fps',
        help='the input filepaths (comma-separated)'),
]
script_info['optional_options'] = [
 make_option('--suppress_errors',action='store_true',\
        help='Suppress warnings about missing files [default: %default]',
        default=False)
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    suppress_errors = opts.suppress_errors

    input_fps = []
    for input_fp in opts.input_fps.split(','):
        input_fps.extend(glob(input_fp))

    for input_fp in input_fps:
        i = 0
        try:
            input_f = open(input_fp,'U')
        except IOError,e:
            if suppress_errors:
                continue
            else:
                print input_fp, e
        for s in MinimalFastaParser(input_f):
            i += 1
        print input_fp, i

if __name__ == "__main__":
    main()