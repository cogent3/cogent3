Scripting guidelines
====================
Developing command line interfaces for your scripts is a convenient way to create easily reproducible analyses with PyCogent. This document covers the support for developing standardized interfaces in PyCogent. In addition to making your code easier to distribute (or to return to months after it was originally written), several GUI generators that are currently in development make use of the PyCogent ``option_parsing`` module -- this means that defining your interfaces with PyCogent's ``option_parsing`` module will allow your code to be easily wrapped in graphical interfaces.

PyCogent command line interfaces are based on the ``optparse`` module in the Python Standard Library, using ``cogent.util.option_parsing.parse_command_line_parameters`` as a convenience wrapper for defining interfaces based on certain information provided by the script developer. A fully functional `example script <./scripting_guidelines.html#countseqs>`_ and a `script template <./scripting_guidelines.html#scripttemplate>`_ are provided below.

This page will help you develop your command line interfaces after you've familiarized yourself with ``optparse``, but will not be a replacement to a general understanding of how to interact with ``optparse``. You should therefore refer both to the `optparse documentation <http://docs.python.org/library/optparse.html>`_ as well as the `PyCogent Scripting Guidelines <./scripting_guidelines.html>`_. As support for optparse will not continue into Python 3.0, we will be switching to ``argparse`` when we transition PyCogent to Python 3. We'll do our best to minimize the work in transitioning from ``optparse`` to ``argparse`` for PyCogent scripts by changing the interface to ``cogent.util.option_parsing.parse_command_line_parameters`` as little as possible.

This document starts with a basic example of creating a script. It then covers guidelines for defining scripts which all PyCogent developers must adhere to, and which we suggest all PyCogent script developers adhere to. Finally it covers some details of the ``script_info`` object and custom command line option types defined in PyCogent.

You should also review the PyCogent `coding guidelines <./coding_guidelines.html>`_. The scripting guidelines presented here are an extension of those.

Quickstart : how to create your first PyCogent script in 5 minutes
------------------------------------------------------------------

These steps show you how to quickly create a working PyCogent script from the PyCogent script template. Some very basic familiarity with python will help here (e.g., identifying the ``main()`` function and adding to it with the correct indentation).

1. Copy and paste the `script template <./scripting_guidelines.html#scripttemplate>`_ to a new file on your system. The next steps will assume that you're saving that file as ``/home/pycogent_user/my_cogent_script.py``.
2. Open a command terminal and enter the following command::

	chmod 755 /home/pycogent_user/my_cogent_script.py

3. Open ``my_cogent_script.py`` with a text editor (such as emacs, pico, or TextWrangler) and add the following to the bottom of the main() function (be sure you have proper indentation)::

	print "Hello World!"

Next, after the line::

	script_info['version'] = __version__

add the following line::

	script_info['help_on_no_arguments'] = False

4. Run and print help text::

	python /home/pycogent_user/my_cogent_script.py -h

5. Run script -- ``Hello World!`` will be printed to the terminal::

	python /home/pycogent_user/my_cogent_script.py

You've now got a working PyCogent script. You can continue working with this script to add the functionality you're interested in. See the `optparse documentation <http://docs.python.org/library/optparse.html>`_ documentation for discussion of how to create options, as well as the `example script <./scripting_guidelines.html#countseqs>`_ provided below. 

General notes on designing command line interfaces
--------------------------------------------------

Design convenient command line interfaces. The goal of your interface is to make things easy for the user (who is often you). This section covers some guidelines for how to do that.

If your script is difficult to work with, or has requirements that are not intuitive for users who frequently work with command line applications, people won't use your code. Have people who are better programmers than you interact with your command line interface and give you feedback on it. 

If there are tasks that are automatable, automate them. For example, if you can make a good guess at what an output file should be named from an input file and a parameter choice, do that and use it as the default output path (but allow the user to overwrite it with a command line option).

Define sensible default values for your command line options. If most of the time that a script is used it will require a parameter to be set to a certain value, make that value the default to simplify the interface.

Have the user specify named options rather than positional arguments. The latter are more difficult to work with as users need to remember the order that they need to be passed. PyCogent scripts do not allow positional arguments by default, but if you must use them you can override this behavior by setting ``script_info['disallow_positional_arguments'] = False``. Note that this contradicts what the ``optparse`` docs say - we disagree with their comment that all required options should be passed as positional arguments. 

Avoid making assumptions about how a script will be run. Perhaps most importantly, don't assume that the script will be run from the same directory that the script lives in. Users often want to copy executables into a centralized directory on their system (e.g., ``/usr/local/bin``). Facilitate that by not requiring that the script is run from a specific location. If you rely on data files, you have other options such as having users set an environment variable that defines where data files live on the system. Test your script from multiple locations on the file system!

Designing PyCogent command line interfaces
------------------------------------------

This section covers guidelines for how to build PyCogent command line interfaces using the ``script_info`` dictionary and the ``cogent.util.option_parsing.parse_command_line_parameters`` function. Some of this is general to ``optparse`` and some is specific to PyCogent.

Flag options
^^^^^^^^^^^^

Flags are boolean options to your script. ``optparse`` supports these directly, so you should never have to define an option that explicitly takes ``True`` or ``False`` on the command line.

Flags to your script should always be either ``action='store_true'`` or ``action='store_false'``, and do not need to define a type. The names of these options should suggest whether the option enables something (e.g., ``--print_to_stdout``) which would be defined with ``action='store_true'`` (i.e., default is False), or whether the option disables something (e.g., ``--suppress_stdout``) which would be defined with ``action='store_false'`` (i.e., the default is True). A bad name for a flag is ``--stdout`` as it's not clear what this option does.

Always define ``default`` for boolean options to set the default option for your script. If ``action='store_true'`` you should *always* pass ``default=False``. If ``action='store_false'`` you should *always* pass ``default=True``.

Choice options
^^^^^^^^^^^^^^
Use ``type=choice`` when an option is passed as a string and can be one of several acceptable values. This saves you from having to check that the user passed an acceptable value. This is done by ``optparse``, so saves you lines of code that you'd need to test, and standardizes how errors are handled. The acceptable choices are defined with ``choices=``. An example choice option definition is::

	alignment_method_choices = ['pynast','mafft','muscle']
	o = make_option('-m','--alignment_method',type='choice',
	                help='Method for aligning sequences. Valid choices are: '+\
	                ', '.join(alignment_method_choices) + ' [default: %default]',
	                choices=alignment_method_choices, default='pynast')

Note that the help text here includes the list of acceptable options. This is generally a good idea as it's convenient for the user. It's not a good idea however if this is a big list (say, more than 5 or so options). If the user passes something invalid (such as ``raxml`` in this example) the list of acceptable options will be included in the error text.

Defining where output will be stored
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a single file is created, allow the user to define that file name. If multiple files are created, allow the user to define a directory name and store all of the output files in that directory. Use the ``new_filepath`` and ``new_dirpath``, respectively, to define these output types. These will raise errors if the file or directory already exists, which is generally good as it avoids overwriting results that may have taken a long time to generate.

Defining options
----------------

Use ``make_option`` (`described here <http://docs.python.org/library/optparse.html#populating-the-parser>`_) to create options. As in that example, you'll define these in lists that get set as ``script_info['required_options']`` and ``script_info['optional_options']``.

Use the PyCogent custom option types when specifying input file or directory paths. These standardize error handling in the case of input files which don't exist or aren't readable and output files which already exist.

Don't define ``dest=``. By default this gets set to the long-form parameter option (e.g. ``dest='input_fp'`` is implied if your option is ``--input_fp``). Defining this as something else will confuse other people who may end up doing maintenance work on your scripts in the future.

Always define ``default=`` for optional options, and never define ``default=`` for required options. The default value for all options is ``None``, but it's convenient to explicitly define that for readability.

Always define ``help=``, and provide useful information in this string. Include ``[default: %default]`` for optional options, but not for required options (as there can be no default for a required option, or it'd be optional). The ``%default`` gets replaced with the value provided for ``default=``. It sometimes makes sense to include additional information in the ``[default:%default]`` text if the option on it's own is not informative. For example::

	make_option("--output_fp",default=None,help="output filepath [default:%default; print to stdout]")

``action=store`` and ``type=string`` are defaults, and therefore do not need to be included. Leave these values out to keep your code cleaner.

If you need to pass multiple paths or strings to a single option, do this by passing a comma-separated string. The ``existing_filepaths`` option type expects strings in this format and takes care of splitting them on commas and returning a list, so if you're passing multiple input filepaths set ``type='existing_filepaths'``.

Naming options
--------------

``optparse`` allows for users to define short-form (e.g., ``-i``) and long-form (``--input_fp``) option names. For options that are commonly used, define both a long-form and a short-form parameter name::

	make_option('-i','--input_dir',type="existing_filepath",help='the input directory')

For options that are infrequently used define only a long-form parameter name::

	make_option('--output_file_type',help='the file type for graphical output',default='pdf')

This helps with reducing clutter and saving convenient short-form parameter names for future options that may be added.

Make paths to files end with ``_fp`` and paths to directories end with ``_dir``. This helps users understand exactly what must be passed to a script.

Some standard names for common options are listed below. You should use these whenever possible.

+-------------------------------+----------------------------------------------------------------------------------------------------+
|        Description            | Option name                                                                                        |
+===============================+====================================================================================================+
|  path to an input file        | ``-i``, ``--input_fp``                                                                             |
+-------------------------------+----------------------------------------------------------------------------------------------------+
|  path to an output file       | ``-o``, ``--output_fp``                                                                            |
+-------------------------------+----------------------------------------------------------------------------------------------------+
|  path to an input directory   | ``-i``, ``--input_dir``                                                                            |
+-------------------------------+----------------------------------------------------------------------------------------------------+
|  path to an output dir        | ``-o``, ``--output_dir``                                                                           |
+-------------------------------+----------------------------------------------------------------------------------------------------+
|  path to a log file           | ``-l``, ``--log_fp``                                                                               |
+-------------------------------+----------------------------------------------------------------------------------------------------+

What documentation should be included in my scripts?
----------------------------------------------------

The ``script_documentation`` entry in ``script_info`` should describe the basic functionality of your script. This entry is typically one to several sentences. Be sure not to add line breaks yourself - ``optparse`` will take care of this for you, and the formatting will look better than if you try to do it yourself.

The ``usage_examples`` entry in ``script_info`` should list one or more examples of commands that need to be run to execute your script. These should be actual calls to commands. A user should be able to copy this and paste it on the command line and have the script run (provided they put the right input files in place). See the `example script <./scripting_guidelines.html#countseqs>`_ for instances of what good usage examples look like. ``script_info['usage_examples']`` must be a list of tuples with three string entries each where the first entry is a concise title for the example, the second entry is a description of the example and why certain parameter settings are being made, and the third entry should be the exact command that needs to be run. Start these examples with ``%prog`` - this gets replaced with the name of your script and is convenient so you don't have to remember to update the usage examples if the name of your script changes.

The ``output_description`` entry in ``script_info`` should describe the output generated by the script. This entry is typically one to several sentences. Again, don't add line breaks yourself.

The script_info dictionary
--------------------------
The ``script_info`` dictionary is the central piece of information required to define a cogent script. ``script_info`` is passed to ``parse_command_line_parameters`` to define the command line interface for your script. Additionally several tools have been developed to import and use this object to define other types of interfaces (e.g., script form in the PyCogent beta GUI) or to auto-generate script documentation (e.g., for the QIIME project). This section covers the values that can be defined in your ``script_info`` dictionaries, what they do, and their default values.


Core values defined in PyCogent command line interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are the core values defined in the ``script_info`` dictionary used by the PyCogent ``option_parsing`` module.

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


Values known to be used by the tools outside of the PyCogent codebase
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These values are known to be used by tools outside of the PyCogent code base in ``script_info`` objects. It's best to not name new values with these names to avoid conflicts. 

+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
|        key                    |  Description                                                                                       |    Used by   |
+===============================+====================================================================================================+==============+
|  brief_description            | a one-sentence description of the script, used by some document generators                         |    Q         |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| script_type                   | a definition of the type of script, used by some graphical interfaces                              |      Q,PG    |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| script_name                   | a brief "human readable" name for the script, used in some graphical interfaces                    |       Q,PG   |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| output_type                   | a list of tuples noting the type (in a controlled vocabulary) of each possible output              |       Q      |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+
| option_label                  | a dictionary matching option names to "human readable" names, used in some graphical interfaces    |   Q          |
+-------------------------------+----------------------------------------------------------------------------------------------------+--------------+

* "Used by" key : Q: `QIIME <http://www.qiime.org>`_; PG: PyCogent beta GUI.

Setting values in script_info
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``script_info`` object is simply a dict, so the standard method for setting and working with dict entries applies. Some examples are::

	script_info['brief_description'] = "Count sequences in one or more fasta files."
	script_info['required_options'] = [
	 make_option('-i','--input_fps',
	        help='the input filepaths (comma-separated)'),
	]

Custom command line option types
--------------------------------
Several custom option types are defined in PyCogent. These are:

* ``existing_path`` : Specify a path to a directory or file. Path must exist or an error is raised.

* ``new_path`` : Specify a path to a directory or file. Path must not exist or an error is raised.

* ``existing_filepath`` : Specify a path to a file.  Path must exist or an error is raised.

* ``existing_filepaths`` : Specify a comma-separated list of file paths. All paths must exist or an error is raised. These are returned as a list split on commas.

* ``new_filepath`` :  Specify a path to a file.  Path must not exist or an error is raised.

* ``existing_dirpath`` :  Specify a path to a directory.  Path must exist or an error is raised.

* ``new_dirpath`` :  Specify a path to a directory.  Path must not exist or an error is raised.

.. _scripttemplate:

Template for a new PyCogent script
----------------------------------
The following is a template for a PyCogent script. You can cut/paste this into a new file to form the basis of your new script. This template forms a fully functional PyCogent script, so on copying this you should be able to run the script to confirm that it is working. If you name your script ``my_cogent_script.py``, you should be able to execute this as follows::

	python my_cogent_script.py

This will print help text and exit.

PyCogent script template::
	
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

.. _countseqs:

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
	
