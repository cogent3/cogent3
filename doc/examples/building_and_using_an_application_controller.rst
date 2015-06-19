.. _new-app-controller:

***********************************************
Building and using a new application controller
***********************************************

.. sectionauthor:: Greg Caporaso

.. note: The code in this file is specifically not doctested because it is describing how to define a new application controller, which won't be in cogent. This file won't be available to import in most/any cases, so shouldn't be tested.

Overview
========
This document provides an example for defining and using a new application controller. We'll look at wrapping the ``formatdb`` application from the BLAST 2.2.20 package `available from NCBI <http://www.ncbi.nlm.nih.gov/BLAST/download.shtml>`_. (Note this is what NCBI now refers to as *legacy BLAST*, not BLAST+.) 

This document was developed in the process of writing the full ``formatdb`` application controller in PyCogent. You can find that file in your PyCogent directory at: ``cogent/app/formatdb.py``. After you work through this example, you should refer to that file to see what the full application controller and convenience wrappers look like.

A more complete reference on `PyCogent Application Controllers can be found in :ref:`app-controllers`.

Building a formatdb application controller
==========================================

Step 0. Decide which application you want to support, and which version.
------------------------------------------------------------------------

Decide what version of the program you want to support and install the application. Check out the features, and decide what functionality you want to support::
	
	formatdb --help

For the sake of brevity in this example, we'll support only the basic functionality: creating nucleic acid or protein blast databases in the default format, and with default names (i.e., named based on the names of the input file). 

So, we'll support the ``-i``, ``-o``, and ``-p`` parameters.

Step 1. Define the class and the class variables.
-------------------------------------------------
First, create a new file called ``minimal_formatdb.py``. Open this in a text editor, and add the following header lines::

	#!/usr/bin/python
	
	from cogent.app.util import CommandLineApplication, ResultPath
	from cogent.app.parameters import ValuedParameter

This imports some classes that we'll be using from PyCogent. For this to function correctly, ``cogent`` must be in your ``$PYTHONPATH``. Next define your class, ``MinimalFormatDb``::

	class MinimalFormatDb(CommandLineApplication):
	    """ Example APPC for the formatdb program in the blast package """

Note that you're defining this class to inherit from ``CommandLineApplication``. This is where most of the heavy lifting is done.

We'll next want to define the following class variables:

	* ``_command``: the command called to start the application
	* ``_parameters``: the parameters to be passed (we'll select a few that we're going to support)
	* ``_input_handler``: function describing how to convert a single parameter passed to the app controller object into input for the function
	* ``_suppress_stdout``: (if other than False)
	* ``_suppress_stderr``: (if other than False)


This is done by adding the following lines::

    _command = "formatdb"
    _parameters = {\
        '-i':ValuedParameter(Prefix='-',Name='i',Delimiter=' ',IsPath=True),
        '-o':ValuedParameter(Prefix='-',Name='o',Delimiter=' ',Value='T'),
        '-p':ValuedParameter(Prefix='-',Name='p',Delimiter=' ',Value='F')
    }
    _input_handler = "_input_as_parameter"

You'll see here that we're only defining the ``-i``, ``-o``, and ``-p`` parameters, hence the name of this call being ``MinimalFormatDb``. An important variable to note here is ``_input_handler``. We'll come back to that next. 

An addition thing to note here is that I'm setting Value for ``-p`` to ``F`` instead of the default of ``T``. This is because I usually build nucleotide databases (specified by ``-p F``), not protein databases (which is the default, and specified by ``-p T``).


Step 2. Overwrite methods as necessary.
---------------------------------------
We'll next create a non-default input handler. The input handler takes the input that you'll eventually provide when calling an instance of ``MinimalFormatDb``, and prepares it be passed to the actual command line application that you're calling. ``formatdb`` requires that you pass a fasta files via the ``-i`` parameters, so we'll define a new function ``_input_as_parameter``, here::

    def _input_as_parameter(self,data):
        """ Set -i to data
        """
        self.Parameters['-i'].on(data)
        return ''

Input handlers return a string that gets appended to the command, but turning parameters on also caused them to be added to the command. For that reason, this input handler returns an empty string -- otherwise ``-i`` would end up being passed twice to ``formatdb``.

Finally, we'll define the ``_get_result_paths``, which is the function that tells the application controller what files it should expect to be written, and under what circumstances. Our ``MinimalFormatDb`` application controller writes a log file under all circumstances, and nucleotide databases when ``-p`` is set to ``F`` or protein databases when ``-p`` is set to ``T``. We return the resulting files as a dict of ResultPath objects::

    def _get_result_paths(self,data):
        """
        """
        result = {}
        
        result['log'] = ResultPath(\
         Path=self.WorkingDir+'formatdb.log',IsWritten=True)

        if self.Parameters['-p'].Value == 'F':
            extensions = ['nhr','nin','nsi','nsq','nsd']
        else:
            extensions = ['phr','pin','psi','psq','psd']
            
        for extension in extensions:
            result[extension] = ResultPath(\
             Path=data + '.' + extension,\
             IsWritten=True)
        return result

At this stage, you've created an application controller which supports interacting with a few features of the ``formatdb`` command line application controller. In the next step, we'll look at how to use your new application controller.

Using the new formatdb application controller
=============================================

Next we'll import the new ``minimal_formatdb`` application controller, and test it out. For the following examples, you need to access some files that are in your ``cogent/doc/data`` directory. For simplicity, we'll assume that on your system this directory is ``/home/pycogent_user/PyCogent/cogent/doc/data``. You should always replace this directory with the path as it exists on your machine.

Open a python interpreter in the directory where you created your ``minimal_formatdb.py`` and enter the following commands::

	>>> import minimal_formatdb
	>>> fdb = minimal_formatdb.MinimalFormatDb()
	>>> res = fdb('/home/pycogent_user/PyCogent/doc/data/refseqs.fasta')
	>>> res
	
You'll see that you've created a new protein BLAST database -- you can tell because you have the nucleotide database files in the result object (i.e., they begin with ``n``).

Next clean up your the files that were created::

	>>> res.cleanUp()

Next we'll change some parameters settings, and confirm the changes::

	>>> fdb = minimal_formatdb.MinimalFormatDb()
	>>> fdb.Parameters['-p'].on('T')
	>>> fdb.Parameters['-p'].isOn()
	True
	>>> fdb.Parameters['-p'].Value
	'T'
	>>> str(fdb.Parameters['-p'])
	'-p T'
	
We've just set the -p parameter to F, indicating that a protein database should be built instead of a nucleotide database. Note that the database files now begin with ``p``. Run the appc and investigate the results::

	>>> res = fdb('/home/pycogent_user/PyCogent/doc/data/refseqs.fasta')
	>>> res

Next clean up your the files that were created::

	>>> res.cleanUp()


Tips and tricks when writing applications controllers
=======================================================
One of the most useful features of application controller object when building and debugging is the HALT_EXEC parameter that can be passed to the constructor. This will cause the program to halt just before executing, and print the command that would have been run. For example:

	>>> fdb = minimal_formatdb.MinimalFormatDb(HALT_EXEC=True)
	>>> res = fdb('/home/pycogent_user/PyCogent/doc/data/refseqs.fasta')
	[Traceback omitted]
	Halted exec with command:
	cd "/home/pycogent_user/"; formatdb -o T -i "/home/pycogent_user/PyCogent/doc/data/refseqs.fasta" -p F > "/tmp/tmpBpMUXE0ksEhzIZA1SSbS.txt" 2> "/tmp/tmpSKc0PRhTl47SZfkxY0g1.txt"

You can then leave the interpreter and paste this command onto the command line, and see directly what happens if this command is called. It's usually useful to remove the stdout and stderr redirects (i.e., everything after and including the first ``>``). For example::

	cd "/home/pycogent_user/"; formatdb -o T -i "/home/pycogent_user/PyCogent/doc/data/refseqs.fasta"

