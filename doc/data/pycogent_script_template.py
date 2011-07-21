#!/usr/bin/env python
# File created on 15 Jul 2011
from __future__ import division

__author__ = "AUTHOR_NAME"
__copyright__ = "COPYRIGHT_INFORMATION"
__credits__ = ["AUTHOR_NAME"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "AUTHOR_NAME"
__email__ = "AUTHOR_EMAIL"
__status__ = "Development"



from cogent.util.option_parsing import parse_command_line_parameters, make_option

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # Example required option
 #make_option('-i','--input_dir',type="existing_filepath",help='the input directory'),\
]
script_info['optional_options'] = [\
 # Example optional option
 #make_option('-o','--output_dir',type="new_dirpath",help='the output directory [default: %default]'),\
]
script_info['version'] = __version__



def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)


if __name__ == "__main__":
    main()