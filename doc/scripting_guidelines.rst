Scripting guidelines
====================
Developing command line interface script is a convenient way to create easily reproducible analyses with PyCogent. This document covers the support for developing standardized interfaces in PyCogent. In addition to make your code easier to distribute (or to return to months after it was originally written), several GUI generators make use of the PyCogent option_parsing module (such as the PyCogent GUI) -- this means that defining your interfaces with PyCogent's option_parsing module will allow your code to be easily wrapped in graphical interfaces.


The script_info dictionary
--------------------------
The script_info dictionary is the central piece of information required to define a cogent script. ``script_info`` is passed to ``parse_command_line_parameters`` to define the command line interface for your script. Additional several tools have been developed to import and use this object to define other types of interfaces (e.g., script form in t\he PyCogent GUI) or to auto-generate script documentation (e.g., for the QIIME project). This section covers the values that can be defined in your ``script_info`` dictionaries, what they do, and their default values.


Core option types used in PyCogent command line interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are the core options used by the PyCogent ``option_parsing`` module.

+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
|        key                    |  Description                                                                                       |    Default   |
+===============================+====================================================================================================+==============+
| script_description            | a paragraph description of the script's functionality                                              |    REQUIRED  |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| script_usage                  | a list of tuples illustrating example usages of the script                                         |       []     |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| output_description            | a paragraph description of the script's output                                                     |       ""     |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| version                       | a version number for the script                                                                    |   REQUIRED   |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| required_options              | a list of optparse Option objects that are required for the script to run                          |        []    |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| optional_options              | a list of optparse Option objects that are optional for the script to run                          |        []    |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| disallow_positional_arguments | do not allow positional arguments to be passed to the script                                       |  True        |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| help_on_no_arguments          | print help text if the script is called with no options or arguments                               |   True       |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| suppress_verbose              | do not auto-generate a verbose option for the script                                               |    False     |  
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+


Options known to be used by the tools outside of the PyCogent codebase
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These options are known to be used by tools outside of the primary PyCogent code base. It's best to not name new options with these names to avoid conflicts. 

+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
|        key                    |  Description                                                                                       |    Used by   |
+===============================+====================================================================================================+==============+
|  brief_description            | a one-sentence description of the script, used by some document generators                         |    Q         |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| script_type                   | a definition of the type of script, used by some graphical interfaces                              |      Q,P     |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| script_name                   | a brief "human readable" name for the script, used in some graphical interfaces                    |       Q,P    |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| output_type                   | a list of tuples noting the type (in a controlled vocabulary) of each possible output              |       Q      |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| option_label                  | a dictionary matching option names to "human readable" names, used in some graphical interfaces    |   Q          |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+

* "Used by" key : Q: `QIIME <http://www.qiime.org>`_; P: PyCogent GUI.




Example of a simple cogent script
---------------------------------

::
	
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
	 ("",
	  "Count the sequences in a fasta file and write results to stdout.",
	  "%prog -i in.fasta"),
	 ("",
	  "Count the sequences in two fasta files and write results to stdout.",
	  "%prog -i in1.fasta,in2.fasta"),
	  ("",
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
	