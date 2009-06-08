********************************************************
  Application Controller Documentation  
********************************************************

.. sectionauthor:: Greg Caporaso, Sandra Smit

About this document
===================

This document was written by Greg Caporaso and Sandra Smit, the developers of
the application controller framework in PyCogent. Since its beginnings, the 
application controller framework has been
further developed and improved by multiple people involved in the PyCogent
project, and of course many application controllers have been implemented.
Thanks!

The purpose of this document (together with the code in the PyCogent library) is
to serve as a guideline for new users who want to implement application
controllers. It should also be helpful as a reference document for more
experienced coders. To meet the different needs of these people, this document
is split up into multiple sections that discuss the framework in different
levels of detail. We first provide a summary, then documentation for application
controller users, then documentation for application controller developers, and
finally a developer's reference detailing the underlying objects.

Feedback from users of this document will help us to improve it, so please don't
hesitate to contact us: gregcaporaso@gmail.com or sandra_smit at users.sourceforge.net.

.. % ============================================================================


The application controller framework: brief overview
====================================================


Summary
-------

An application performs a certain task. It uses parameters to give the user
control over the possible settings. In general an application requires some
input data and produces some output. To run RNAfold on (i.e. predict an RNA
secondary structure for) the sequence in seqs.txt, using some settings for
temperature and GU pairs), and to save the output in the file result.txt, the
following command would have to be executed::

   RNAfold -T 25 -noGU <seqs.txt >result.txt

Suppose you want to fold the sequences at two different temperatures (25 and 37
degrees Celcius) and investigate the difference between the structures. This
would require a lot of manual work (run RNA fold twice, store all the files, use
python to read in the files, parse the output, and calculate the difference).
This task becomes a lot easier if we can control the application RNAfold from
within the python environment. The following example folds two sequences at two
different temperatures and reports the symmetric difference between the
resulting structures. It works by creating a single application controller, and
running  it with different settings on the same input. The output is parsed and
compared. The code snippet doesn't leave any files on the system.

.. % This command would run RNAfold on (i.e. predict an RNA secondary structure for) the sequences in seqs.txt. It would write the output to the file result.txt. And it would use the specified settings (for temperature and GU pairs).

::

   seqs = ['GCCCGGAUAGCUCAGUUAAAAUCCCCGUGUCGGUGGUUCAAUUCCGCCUCCGGGCCCA',\
   'CCCCUGGUGGCCAUAGCGGGGAUUAACCCGGCCUGUGCCGGCAGCGGCACGGGAGGUCGCUGCCAGGGG']
   r = RNAfold(WorkingDir="/tmp")

   r.Parameters['-T'].on(25)
   print "COMMAND AT 25C:", r.BaseCommand
   res = r(seqs)
   fold_25 = just_seqs_parser(res['StdOut'].readlines())
   res.cleanUp()

   r.Parameters['-T'].on(37)
   print "COMMAND AT 37C:", r.BaseCommand
   res = r(seqs)
   fold_37 = just_seqs_parser(res['StdOut'].readlines())
   res.cleanUp()

   for seq_25, struct_25 in fold_25:
       for seq_37, struct_37 in fold_37:
           if seq_25 == seq_37:
               print 'SEQ:', seq_25
               print '25C:', struct_25
               print '37C:', struct_37
               print 'BASE PAIR SIMILARITY: %.3f'%\
                   (compare_pairs(struct_25.toPairs(), struct_37.toPairs()))

Running the script results in the following output::

   COMMAND AT 25C: cd "/tmp/"; RNAfold -d1 -T 25 -S 1.07
   COMMAND AT 37C: cd "/tmp/"; RNAfold -d1 -T 37 -S 1.07
   SEQ: GCCCGGAUAGCUCAGUUAAAAUCCCCGUGUCGGUGGUUCAAUUCCGCCUCCGGGCCCA
   25C: (((((((..((..((((.....(((((...))).))...))))..)).)))))))... (-26.25)
   37C: ((((((.((((...)))).............(((((.......))))).))))))... (-19.6)
   BASE PAIR SIMILARITY: 0.222
   SEQ: CCCCUGGUGGCCAUAGCGGGGAUUAACCCGGCCUGUGCCGGCAGCGGCACGGGAGGUCGCUGCCAGGGG
   25C: (((((((((((.((..((((......)))).(((((((((....)))))))))..)).))))))))))) (-55.52)
   37C: (((((((..((.....((((......)))).(((((((((....))))))))).....))..))))))) (-46.3)
   BASE PAIR SIMILARITY: 0.846

It has tremendous benefits to have the capability to run an application in an
integrative fashion from the Python environment, using Python object as input,
catching and directly using the output in a downstream analysis. It is
timesaving and allows for the full automation of an analysis (prevents manual
intervention). This is extremely convenient, for example, if you wish to perform
a comparison of a large variety of parameters passed to a single application.
Rather than manually running the application many times from the command line,
you can script the runs and automatically process the output. A piece of python
code that handles this is called an application controller (APPC). The purpose
of the application controller framework is to support various application
controllers in general. The framework is generic and supports actions that all
applications have in common. Each application requires a specific APPC that
defines all the variables and methods for that unique application.

At the moment the framework supports only one type of application: the command
line application. In the (near) future we want to support web applications as
well.

This section addresses the different components of an application and how these
components are represented in the framework. The framework is set up in two
files: cogent/app/util.py and cogent/app/parameters.py


Command line applications
-------------------------

Every application requires a base command (such as RNAfold, muscle, or clustalw)
and some parameters. These two properties are stored in the high level class
Application.

The class CommandlineApplication inherits from Application and has additional
features. It deals with the working directory, input, and output. The central
method is the __call__ method in which the full command-to-execute is built up
and executed. The result of running the application controller is returned to
the user.

The util.py file contains one other class: ApplicationError. This class is used
to raise exceptions for application controllers.


Parameters
----------

Most applications allow you to specify a certain set of parameters to control
how the program runs. Parameters can control many different features of an
application, such as the temperature at which RNA is folded, the number of gaps
allowed in an alignment, or the name of an output file. They come in many forms
as well, some are simply flags, some always require a value, some can have
optional values.

.. % Parameter
.. % -- FlagParameter
.. % -- ValuedParameter
.. % -- MixedParameter
.. % Parameters
.. % ParameterError
.. % FilePath

The application controller framework supports three types of parameters, which
will be discussed below. Subclassing to specify new types of parameters or to
make certain attributes fixed, is very easy.

The abstract Parameter class defines the basic functionality of a parameter: it
initializes all the Parameter attributes and it defines a Parameter ID which is
a unique identifier for each parameter. In general a parameter has a prefix
(usually a dash) and a name. Some parameters have values. The Parameter object
is discussed in more detail in section :ref:`sec:build`.

There are three subclasses from the class Parameter. FlagParameter is used for
parameters that don't have values (e.g. allow GU pairs or not). ValuedParameters
are used for paramaters that specify some value (e.g. the temperature or some
input file). MixedParameters are parameters that might or might not have a value
(e.g. the -d parameter in RNAfold). All parameters of an application are grouped
in a Paramters object. The class Parameters is a special type of dictionary that
allows lookups by parameter ID or synonyms.

The parameters.py file contains two more classes. ParameterError is used to
raise exceptions in the parameter framework. The class FilePath defines paths on
a system, it can print itself in a special way and add other parts of a path.


Input
-----

Input can be very diverse between applications. Most often it requires a file or
some data directly from the command line. Application input is handled by "input
handlers". There are a few generic input handlers in CommandLineApplication
object. Specific APPCs can use these methods directly or overwrite them. The
methods process the input data for the application. They might for example write
a certain Python object to a temporary file, and change some application
parameters to use this file.


Output
------

All applications produce some form of output. It can be limited to information
on "standard output" (stdout) and "standard error" (stderr). Many applications
produce additional output files. Most (unfortunately not all) applications
report a meaningful exit status that inform the user on whether the execution of
the program was succesful. The class CommandLineAppResult handles all aspects of
application output: stdout, stderr, exit status, and the additional output
files. Access to all the available files is handles by the class ResultPath.
More technical aspect of these classes is discussed in section :ref:`sec:build`.

.. % ResultPath
.. % CommandLineAppResult

.. % ============================================================================
.. % \newpage


Using an application controller
===============================


Summary
-------

#. Create an instance of some app controller

#. Turn parameters on and off

#. Optionally change the working directory

#. Optionally check the base command which is built-up from the above
   information

#. Set the input handler

#. Possibly redirect StdOut and StdErr

#. Apply the instance to the input data, store the results

#. Use the results as you like

#. Possibly clean up files created by the program and the APPC


Creating an instance with basic settings (parameters, working directory)
------------------------------------------------------------------------

The first step toward running an application is creating an instance of the
APPC. Two basic settings are the parameters and the working directory. Below are
some examples on how to do this. Note that the working directory must be an
absolute path.

All parameters have the methods isOn and isOff to check whether the parameter is
active or not. Parameters can be turned on and off (with or without a value)
with the  on() and off() methods. The values can also be set during
initialization of the APPC. When specifying the parameters upon initialization
the __init__ parameter params should be a dictionary of parameters that should
be turned on, keyed by either the Parameter ID or a synonym. The values in
params should be the values to turn the parameters on with for Valued or Mixed
Parameters, or None for Flag or Mixed Parameters.

It is useful to check the BaseCommand to see if all the parameters have the
correct settings and if the working directory is correct. During debugging it is
useful to check whether the command runs on the normal command line.

.. % \subsection{Setting/changing parameters}
.. % \subsection{Changing the working directory}

::

   Initialization without params, only defaults are on.
   >>> from cogent.app.vienna_package import RNAfold
   >>> r = RNAfold()

   Initialization with params, set new values for this instance
   >>> r = RNAfold(params={'-T':25,'-d':None,'-4':None,'-S':1.2})

   Initialization changing the Working directory (must be absolute path!)
   >>> r = RNAfold(WorkingDir='/tmp')
   >>> print r.BaseCommand
   cd "/tmp/"; RNAfold -d1 -T 37 -S 1.07

   Changing the working directory after initialization (must be absolute path!)
   >>> r = RNAfold()
   >>> r.WorkingDir = '/tmp'
   >>> print r.BaseCommand
   cd "/tmp/"; RNAfold -d1 -T 37 -S 1.07

   Checking the parameters
   >>> r = RNAfold()
   >>> print r.Parameters['-P'].isOn()
   False
   >>> print r.Parameters['-P']
   <BLANKLINE>
   >>> print r.Parameters['-T'].isOn()
   True
   >>> print r.Parameters['-T']
   -T 37


Other settings on initialization
--------------------------------

The input handler could be set (if not, the default is used)  ::

   On initialization
   >>> r = RNAfold(InputHandler="_input_as_string")
   >>> print r.InputHandler
   _input_as_string

   After initialization
   >>> r = RNAfold()
   >>> print r.InputHandler # default
   _input_as_lines
   >>> r.InputHandler = "_input_as_path"
   >>> print r.InputHandler
   _input_as_path

Standard out and standard error can be suppressed. If SuppressStderr or
SuppressStdout are set to True, stdout and stderr will be routed to /dev/null.
The default is to store these results in a temporary file. Redirecting StdErr
might be useful for programs that write a lot of useless information to this
filestream.

.. % input handler

Some parameters concerning the creation of temporary files can be changed.
TmpDir: default is /tmp. TmpNameLen is the length of the filenames, default is
20.

HALT_EXEC is a parameter that can be set to True for debugging purposes. It
stops the process right before execution of the system call, it leaves all the
input files (incl. temporary) in place. This allows the user to check whether
the input is generated correctly. See Section :ref:`sec:haltexec` for more
details.


Running the application, using the output, and cleaning up
----------------------------------------------------------

When calling the instance of the APPC on some data the __call__ method is
invoked. The call method has to optional parameters: data (the input data) and
remove_tmp (if True the temporary files are removed). The call method returns a
CommandLineAppResult object, containing all the application output information.

The output dictionary can be used to access the resulting files. All the
information can be incorporated in a downstream analysis. In the example below
the aligned sequences in clustalw format are parsed and printed.

Additionally CommandLineAppResult contains one public method: cleanUp() which
takes no parameters.  The method cleanUp() should be used when you want to
delete the files that were created by the CommandLineApplication from disk. Note
that after cleanUp() you may still have access to your files, but these are not
reliable. You will only have access to what has already been loaded into memory
(ie. only a fraction of your file typically), so you should only run cleanUp()
after you are done accessing you files. Also note that running cleanUp() is not
required. If you want the result files to remain on disc you should not run
cleanUp() and they will be left in place. This is useful for running an
application for later analysis of results. ::

   >>> from cogent import PROTEIN
   >>> from cogent.app.clustalw import Clustalw
   >>> from cogent.parse.clustal import ClustalParser
   >>> s1 = PROTEIN.Sequence('MHSSIVLATVLFVAIASASKTRELCMKSL')
   >>> s2 = PROTEIN.Sequence('MALAEADDGAVVFGEEQEALVLKSWAVMKKDA')
   >>> s3 = PROTEIN.Sequence('MSTVEGREFSEDQEALVVKSWTVMKLNAGELALKF')

   >>> c = Clustalw(InputHandler="_input_as_seqs")
   >>> result = c([s1,s2,s3])
   >>> print result['ExitStatus']
   0
   >>> aln_txt = result['Align'].readlines()
   >>> for label, seq in ClustalParser(aln_txt): print "%s: %s"%(label, seq)
   2: MALAEADDGAVVFGEEQEALVLKSWAVMKKDA-------
   3: MSTVEGRE----FSEDQEALVVKSWTVMKLNAGELALKF
   1: MHSSIVLAT-VLFVAIASASKTRELCMKSL---------
   >>> result.cleanUp()

.. % ============================================================================
.. % \newpage


.. _sec:build:

Designing and implementing a new type of application controller
===============================================================

Each specific application that you wish to control through PyCogent requires an
application controller, i.e., a subclass of CommandLineApplication. Building the
new application controller consists of three steps:

#. Creating the application controller class: Overwrite CommandLineApplication
   to define your new application controller, and define the class data. (Section
   :ref:`sec:step1`.)

#. Input handing: Determine whether the built-in input handlers (in
   CommandLineApplication) are sufficient. If not, write one or more input handling
   methods. (Section :ref:`sec:step2`.)

#. Output handling: Determine whether the program writes any output files to
   disk. If so, implement the _get_result_paths method. (Section :ref:`sec:step3`.)


.. _sec:step1:

Step 1: Creating the application controller class and defining class data
-------------------------------------------------------------------------

All of these class variables are discussed in detail in Sections
:ref:`sec:application` and :ref:`sec:commandlineapplication`. ---

**The following class data must be overwritten:**

_command:
   The command used to run the command (a string).

**The following class data can be overwritten:**

_parameters:
   A dictionary of Parameter objects. Keys should be the identifiers of the
   parameters, and values should be the Parameter objects.

_command_delimiter:
   String that specifies the delimiter between the components of a full command,
   e.g. the command, parameters, and arguments.

_synonyms:
   A dictionary of parameter synonyms. Keys should be the alternative keys to
   lookup a parameter, and values should be the identifiers used in the _parameters
   dictionary.

_input_handler:
   The name of the input handler method that should be used by default. The value
   should be a string (see CommandLineApplication.__call__ for how it's used).

_working_dir:
   Specifies where the command should be run (string). Default is current working
   directory.

_suppress_stdout:
   Boolean value that specifies what happens with standard output (stdout) by
   default.

_suppress_stderr:
   Boolean value that specifies what happens with standard error (stderr) by
   default.

Defining parameters
^^^^^^^^^^^^^^^^^^^

All parameters should be one of the three built-in types: FlagParamater,
ValuedParameter, or MixedParameter. (We don't know of any types that wouldn't
fit into this framework, but if you come across any, please let us know.)
Examples illustrating how to define the three different parameter types can be
found in Section :ref:`sec:parameters`. The _parameters dict is a mapping of
parameter identifiers, or Prefix and Name joined by the empty string, to
parameter objects. All parameters which can be passed to an application should
be defined in the parameters dict. Usually you can get this list by reviewing
the application's documentation. See Section :ref:`sec:rnafoldexample` for an
example including the definition of the _parameters dict. Note: if for a given
ValuedParameter or MixedParameter, the value is intended to be a path to a
directory or file, ``IsPath=True`` must be passed when initializing those
parameters. ---  **Defining a new Parameter type** ---  If the application
you're working with uses a type of parameter that is not supported by the
framework yet, you might want to write your own subclass. To subclass Parameter,
the following methods will need to be implemented: __str__, isOn(), isOff(),
on(), off(). These methods cover the two important characteristics of each
parameter: knowing how to print itself, based on its status, and knowing how to
be turned on or off. It is unlikely that you will need to subclass parameter if
working with CommandLineApplication subclasses. If you think you do, please let
us know. ---  **Writing constructor functions/wrappers** ---  There might be
several reasons, such as to make some attribute of the parameter fixed, to write
a wrapper around or constructor function for a parameter. For example to fixate
the prefix of the FlagParameter, one might write this::

   >>> from cogent.app.parameters import FlagParameter
   >>> def DashedFlag(name):
   ...   return FlagParameter('-',name)
   ...
   >>> tree = DashedFlag('tree')
   >>> tree
   <cogent.app.parameters.FlagParameter object at ...
   >>> tree.on()
   >>> print tree
   -tree


.. _sec:rnafoldexample:

A complete Command-LineApplication subclass example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subclass of CommandLineApplication might look something like this::

   class RNAfold(CommandLineApplication):
       """Application controller for RNAfold (in the Vienna RNA package)
       """
       _command = 'RNAfold'
       _parameters = {
       '-p':MixedParameter(Prefix='-',Name='p',Delimiter='',Value=False),
       '-C':FlagParameter(Prefix='-',Name='C'),
       '-T':ValuedParameter(Prefix='-',Name='T',Value=37,Delimiter=' '),
       '-4':FlagParameter(Prefix='-',Name=4),
       '-d':MixedParameter(Prefix='-',Name='d',Delimiter='',Value=1),
       '-noLP':FlagParameter(Prefix='-',Name='noLP'),
       '-noGU':FlagParameter(Prefix='-',Name='noGU'),
       '-noCloseGU':FlagParameter(Prefix='-',Name='noCloseGU'),
       '-e':ValuedParameter(Prefix='-',Name='e',Delimiter=' '),
       '-P':ValuedParameter(Prefix='-',Name='P',Delimiter=' '),
       '-nsp':ValuedParameter(Prefix='-',Name='nsp',Delimiter=' '),
       '-S':ValuedParameter(Prefix='-',Name='S',Value=1.07,Delimiter=' ')}
       _synonyms = {'Temperature':'-T','Temp':'-T','Scale':'-S'}
       _input_handler = '_input_as_lines'
       _suppress_stderr = True 

If the built-in input handlers are sufficient, and no output to disk is written
by the program, this would complete the application controller.


.. _sec:step2:

Step 2: Input handling
----------------------

Not all applications handle their input in the same way. The input might be
specified as a filename on the command line, as a list of values on the command
line,   or an input file might be specified through parameters. Some input data
might also require processing before it is used by the application.

To give the user control over how input is handled without having to overwrite
__call__(), small input handling methods can be specified in the application
controller. In most cases, the CommandLineApplication input handlers can
probably be used (e.g., passing data via stdin or a temp file), but for more
complicated input formats, custom input handlers may need to be written for a
CommandLineApplication subclass. Every input handling method should take one
parameter, data, and return a string that will be appended to the command, e.g.
``/path/to/input/file.txt``, if a path is passed to the application. (In this
example, you would want to use CommandLineApplication._input_as_path as the
input handler.)

By writing multiple input handling methods, multiple types of input can be
handled by one application. The user can specify which one they want to use in a
certain instance by setting the _input_handler class variable, or the
InputHandler initialization variable.

For example, RNAfold takes a list of sequences from stdin. In this case, none of
the built-in input handlers provides this functionaloty. The following input
handler (from cogent.app.rnafold.Rnafold) writes the sequences (data) to a
temporary file and redirects them to stdin. ::

   def _input_as_lines(self,data):
       """Returns '<temp_filename to redirect input to stdin"""
       return ''.join(['<',super(RNAfold,self)._input_as_lines(data)])

Clustalw requires the input filename be passes via the -infile paramter. This
custom input handler from cogent.app.clustalw.Clustalw performs that function.
Note that the empty string is returned, as input handlers are required to return
a string that should be appended to the command line.  ::

   def _input_as_string(self,data):
       """Makes data the value of a specific parameter
       This method returns the empty string. The parameter will be printed
       automatically once set.
       """
       if data:
           self.Parameters['-infile'].on(data)
       return ''

The default input handler should be set (as a string) via the class variable.
See the example in Section :ref:`sec:rnafoldexample`.


.. _sec:step3:

Step 3: Output handling
-----------------------

Stdout and the exit status of any program are caught automatically. Stderr is
accessible as well, unless suppressed via the _suppress_stderr class variable or
the SupressStderr instance varaible. Any other files that are written should be
made accessible by specifying their paths in the method _get_result_paths(). If
you don't overwrite this method, it is assumed that the program doesn't create
additional output files, so if it does, they will be written, but won't be
accessible through the CommandLineAppResult object, and won't be cleaned up upon
program termination!

Names and locations of output files may be fixed, but they can also be created
on the fly based on things such as input file name, data the application is
called on, a combination of values of parameters, or specified filename plus a
fixed suffix. Since the generation of output files is so application specific
and may be very complex, each application controller should handle its own
output.

The _get_result_paths method should take data (as passed to __call__) as an
argument. This is necessary to allow access to any possible variable used by the
program. The user has access to data, self._input_filename (for an on the fly
generated input file), all parameter values, and all public attributes of an
Application.

_get_result_paths() should return a dictionary of ResultPath objects. The file
streams resulting from a run of the application (in the CommandLineAppResult)
will be accessed by the keys in the dictionary. The ResultPath specifies the
*absolute* path of a file and whether the file has been written. This dictionary
is used as input for the CommandLineAppResult which will handle opening the
files etc.

As an example we show the output handling method of RNAfold. For a more complex
example, see RnaView. ::

   def _get_result_paths(self,data):
           """Specifies the paths of output files generated by the application

           data: the data the instance of the application is called on

           You always get back: StdOut,StdErr, and ExitStatus
           RNAfold can produce two additional output files:
               a secondary structure structure graph. Default name: rna.ps
               a dot plot of the base pairing matrix. Default name: dp.ps
           The default names are used for unnamed sequences. Files are created
               in the current working directory.
           You can make a sequence named by inserting a line '>name' above it in
               your input file (or list of sequences). The ss and dp files for 
               named sequences will be written to name_ss.ps and name_dp.ps
           """
           result = {}
           name_counter = 0
           seq_counter = 0
           if not isinstance(data,list):
               #means data is file
               data = open(data).readlines()
           for item in data:
               if item.startswith('>'):
                   name_counter += 1
                   name = item.strip('>\n')
                   result[(name+'_ss')] =\
                       ResultPath(Path=(self.WorkingDir+name+'_ss.ps'))
                   result[(name+'_dp')] =\
                       ResultPath(Path=(self.WorkingDir+name+'_dp.ps'),\
                       IsWritten=self.Parameters['-p'].isOn())
               else:
                   seq_counter += 1

           result['SS'] = ResultPath(Path=self.WorkingDir+'rna.ps',\
               IsWritten=seq_counter - name_counter > 0) #Secondary Structure
           result['DP'] = ResultPath(Path=self.WorkingDir+'dot.ps',
               IsWritten=(self.Parameters['-p'].isOn() and\
               seq_counter - name_counter > 0)) #DotPlot
           return result


Tips and tricks for creating application controllers
----------------------------------------------------


.. _sec:haltexec:

HALT_EXEC is your friend
^^^^^^^^^^^^^^^^^^^^^^^^

The __init__ method takes a boolean parameter, HALT_EXEC, which is False by
default. Setting HALT_EXEC=True will cause __call__ to exit before the system
call, print out the complete command that was about to be run, and leave all
temporary files in place. This is extremely useful for debugging, because it
allows you to run the application directly with the input that was generated by
the application controller. You can therefore run the command and look directly
at stdout and stderr, debug any temporary files that were created, etc. If the
application you are controlling is slow, this can also allow you to debug
earlier steps without having to wait for the application to run. HALT_EXEC is
your friend.

.. % ============================================================================
.. % \newpage


Application controller base classes: Developer's reference
==========================================================


Command line applications
-------------------------


.. _sec:application:

Application: cogent.app. util.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Application is an abstract class that contains some data used by for all
application handlers that could be written. Private class data of Application
consists of:

_command:
   The command used to run the command (a string). If the command is in your path
   (in a directory listed in the environment variable ``$PATH``, or found by typing
   ``which`` followed by the command) you can provide only the command, e.g.
   ``RNAfold``. If the command is not in your path, you must specify the absolute
   path it, e.g. ``/some/other/bin/rnaview``. The use of absolute paths here is
   *not* recommended, because the location of the installation might be different
   on every machine. Instead, consider setting your ``$PATH`` environment variable
   to include the directory where the application is installed.

   For example if you are writing an application controller for ``ls`` where you
   might run: ---   ls -al \*.jpg ---  _command should be set to ``ls``.

_parameters:
   A dictionary of Parameter objects. Keys should be the identifiers of the
   parameters, and values should be the Parameter objects. This dictionary defines
   which parameters are available to the application. No values are specified,
   except for occasional default values. The default value for _parameters is the
   empty dictionary. If the application takes any command line parameters, this
   must be overwritten. This is almost always the case. See
   cogent.app.clustalw.Clustalw._parameters for an example of when parameters is
   overwritten.  ::

      _parameters = {'-T':ValuedParameter('-','T',Delimiter='=')}

_synonyms:
   A dictionary of parameter synonyms. Keys should be the alternative keys to
   lookup a parameter, and values should be the identifiers used in the _parameters
   dictionary. It probably a good idea to comment on the available synonyms in the
   docstring of the application controller, so users that haven't read the manual
   know what they can use to control the parameters. The default value for
   _synonyms is the empty dictionary. See
   cogent.app.vienna_package.ViennaPackage._synonyms for an example of this being
   overwritten.  ::

      _synonyms = {'Temperature':'-T', 'Temp':'-T'}

_command_delimiter:
   String that specifies the delimiter between the components of a full command,
   e.g. the command, parameters, and arguments. The default value is ' ' (a single
   space). This delimiter will work for any Unix application, so it is usually not
   overwritten. (We are interested in hearing about any circumstances where this
   might be overwritten. Please let us know if you come across any. One example
   might be if the command being constructed is a URL.)

   In the above 'ls' example, a single space (' ') spearates the command
   componenets: the base command ``ls``, the parameters ``-al``, and the argument
   ``*.jpg``.

The only method that is defined by Application is __init__, which takes one
optional argument, params. The value of params should be a dictionary of
parameters that should be turned on. Keys should be either the Parameter ID or a
synonym. The values in params should be the values to turn the parameters on
with for Valued or Mixed Parameters, or None for Flag or Mixed Parameters.

Application is never directly instantiated, but is instead inherited (either
directly or indirectly) by all application controllers. It is necessary that
Application.__init__() be called somewhere during the initialization of your
class, but if you are inheriting from a higher level class (such as
CommandLineApplication) this should already be handled.


.. _sec:commandlineapplication:

Command-LineApplication: cogent.app.util.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CommandLineApplication is an abstract class for command line application
controllers. Several variables are class data to facilitate subclassing
CommandLineApplication and to allow definition of defaults for Application
Controller subclasses. This class was designed to be easily and minimally
subclassed.

CommandLineApplication inherits from Application. CommandLineApplication
contains the following additional class data:

_input_handler:
   The name of the input handler method that should be used by default. The value
   should be a string (see CommandLineApplication.__call__ for how it's used). The
   input handling methods are private, so they should start with an underscore. The
   default value for _input_handler is '_input_as_string'. The input handler can be
   changed on instance level via the InputHandler initialization parameter.

_working_dir:
   Specifies the default working directory (string). The working_dir is where many
   applications write out their output. Setting this value gives you control over
   where output is written. The value of _working_dir should be an *absolute* path.
   If the value of _working_dir is None (the default) the current working directory
   will be used. The working directory can be changed at instance level via the
   WorkingDir initialization parameter.

_suppress_stdout:
   Boolean value that specifies what happens with standard output (stdout) by
   default. If the value is False (default), stdout is caught and accessible in the
   result object. If the value is True, stdout is routed to /dev/null and won't be
   accessible. Suppression of stdout can also be controlled at instance level via
   the SuppressStdout initialization parameter.

_suppress_stderr:
   Boolean value that specifies what happens with standard error (stderr) by
   default. Some programs write a lot to stderr which you might want to ignore. If
   the value is False (default), stderr is caught and accessible in the result
   object. If the value is True, stderr is routed to /dev/null and won't be
   accessible. Suppression of stderr can also be controlled at instance level via
   the SuppressStderr initialization parameter.


Class data can be overruled on the instance level by passing alternate data in
as parameters to __init__(). These parameters are InputHandler, SupressStderr,
and WorkingDir.  Note that _working_dir and WorkingDir must always be an
absolute path, although no explicit checking is done for this. You *will* get
weird results in many cases if you use relative paths. WorkingDir, InputHandler,
and SupressStderr are all public attributes of CommandLineApplication, and can
be modified at anytime. You should (obviously) not modify the private versions
of these attributes. Note that if WorkingDir does not exist on the system it
will be created, and it will not be removed after the program runs.

There is an additional private variable _input_filename. This is set to the
string containing the absolute path to an input file when the input file is a
python generated temporary file. This should not be accessed from outside of the
program, but may be useful at times when subclassing.

CommandLineApplication defines several methods. These include::

   __init__(), __call__(), _input_as_string(), _input_as_multiline_string(),
   _input_as_lines(), _input_as_path(), _input_as_paths(), _absolute(),
    _get_base_command(), _get_WorkingDir(), _set_WorkingDir(), _accept_exit_status(),
    _get_result_paths(), getTmpFilename()

We will go over these in differing depths, because for most cases, these are
background methods that should never be called directly, or overwritten.

__init__():
   Initializes the object, taking as parameters params (see Application),
   InputHandler, WorkingDir, SupressStderr (discussed above). This method *must* be
   called by subclasses in their __init__() if they have one. For most purposes,
   you will never need to overwrite this method.

__call__():
   This is the method that does most of the work in the CommandLineApplication.
   Most of a users interaction with CommandLineApplications will be through this
   method, which takes data as a parameter. data is the data that should be passed
   as input to the application when it is called, default is None. Note that before
   data is appended to the command the InputHandler function is called on it. If
   data=None, no data is passed into the function, and the input handler will not
   be called. You should at all costs avoid overwriting __call__() as a lot is
   going on here.

_input_as_string():
   The default input handler. This acts on one parameter, data, that is passed in.
   It type casts data to a string, and returns the string.

_input_as_lines():
   An alternate input handler. In this case, data is a a sequence of lines to be
   written in a temporary file. This allows you interact with programs which only
   takes files as input, when you have created a data file on the fly. The return
   value of this function is a string representing the absolute path to the
   filename, which will be created with in self.WorkingDir.

_input_as_multiline_string():
   Input handler, similar to _input_as_lines, except data is a single string which
   should be written to a temporary file. The temporary file's path is passed as
   input to the application as input.

_input_as_path():
   Another alternate input handler. This is similar to _input_as_string, but casts
   the input to a FilePath object rather than a string. If the input is a path,
   this input handler should be used.

_input_as_paths():
   Yet another alternate input handler. This is similar to _input_as_path, but
   operates on a list of paths.

_absolute():
   Converts a filename to an absolute path if it is not already. The path that is
   appended is self.WorkingDir. The result is a FilePath object.

_get_base_command():
   Appends the necessary parameters to self._command and returns the full command
   as a string (without input and output).

_get_WorkingDir() and _set_WorkingDir():
   accessor methods for the WorkingDir attribute.

_accept_exit_status():
   This function takes a string containing the return value of the application that
   was run. It is meant to be overwritten when necessary. It's purpose is to
   analyze the exit_status of the application being run to determine if an
   ApplicationError should be raised. By default, no ApplicationError is raised
   regardless of the exit_status. In a subclass this is handy because you can
   customize what exit statuses are acceptable to you, and which are not, or you
   can not define the function in your subclass and accept all exit statuses.

_get_result_paths():
   This method is used to initialize the CommandLineAppResult class (see Section
   :ref:`sec:commandlineappresult`). This method should be overwritten if the
   application creates output other than stdout and stderr.  A dict should be
   returned with ResultPath objects keyed by the names that you'd like to access
   their data by in the CommandLineAppResult object. When building the ResultPath
   objects, you will need to construct the names of all of the files that are being
   created. For this reason, you will need access to all of the data that the
   application has access to in the case of dynamic filenames. In order to
   construct these file name you have access to the Parameters object, data (which
   is passed in to the function) in the case where, for example, the output
   filename is specified as input to the program. The name of the input filename,
   when generated as a temporary file is available as self._input_filename, for
   cases where the output file name is based on the name of the input filename.
   This, in addition to system calls if necessary, should provide all of the
   information needed to build the names and paths of output files.

getTmpFilename():
   Generates a random filename using ``TmpLenName`` random alphanumeric (upper and
   lowercase) characters. The result will be an absolute path (presuming that
   ``TmpDir`` is absolute, which it should be), and the filename will begin with
   ``prefix``, end with ``suffix``, and be in
   ``tmp_dir`` or ``TmpDir``. The ``tmp_dir`` parameter
   overrides the class/object-level default. Note that this function does not
   actually created the file, just the filename. The result is a ``FilePath``
   object.

Two module level functions are also implemented::

   get_tmp_filename, guess_input_handler

get_tmp_filename:
   A module level implementation of ``CommandLineApplication.getTmpFilename().``

guess_input_handler:
   This is a module-level function intended to pick the right input_handler in case
   the input is a set of sequences. It will return one of four input handlers: ---
   _input_as_multiline_string, _input_as_path, _input_as_seqs, _input_as_lines.


Parameters
----------


Parameter: cogent.app. parameters.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The class Parameter is an abstract class. Every Parameter object has six
attributes: Prefix, Name, Value, Delimiter, Quote, and IsPath. All attributes
may have any value, as long as it can be type casted into a string.

The Prefix of a parameter specifies the character that precedes the name of the
parameter. It is mandatory to specify a prefix for a parameter, although it may
be the empty string. For example: '-' is the prefix in '-T=37', and '\*' is the
prefix in '\*d'. Note that some characters may have to be escaped (e.g.
`\backslash`).

Name is the second mandatory attribute of Parameter. The combination of the
prefix and name of a parameter should form a unique combination that identifies
the parameter. This ID is a public property of Parameter and will function later
on as the key in the dictionary of parameters.

The attribute Value specifies the value of a parameter. It will be clear that
not all parameters, such as flags, require a value. Therefore this field is
optional in the __init__ method. For example, the value in '-T=37' is 37, the
value in '-d1' is 1.

The Delimiter specifies what separates the name from the value when a parameter
is printed. For example: '=' in 'T=37' or ' ' (single space) in '#r 14'.

The Quote is an optional attribute that determines which characters will
surround the value when the parameter is printed. Be alert on escaping quotes,
since most quote-values will have a special meaning in python. At the moment
only symmetrical quotes are supported, such as " ' " (single quote) in " -p='a'
". Asymmetrical quotes are not possible, e.g. 'd=[4]'. *Is this something that
should get supported?*

IsPath should be set to true if the Value of the Parameter object is intended to
be a path to a directory or file. Paths require special handling when printing,
and Value is therefore cast to a cogent.app.parameters.FilePath object. IsPath
is only used by ValuedParameter and MixedParameter objects, and has no effect on
FlagParameters.

Every type of parameter prints itself differently. A flag will only print a
combination of its prefix and name; another parameter may include everything.
Therefore, the __str__ has to be specified in each specific subclass of
Parameter. Whether a parameter is printed is determined by its value. This is
also subclass specific and will be explained in the following sections.


FlagParameter: cogent.app. parameters.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FlagParameter inherits from Parameter. A flag can't have a value, it is just on
or off. For example: if '-tree' is set, a dendrogram is calculated; if '-tree'
is not set, the tree is not calculated. Since a flag can never have a value, we
can easily use the value to specify whether the flag will be printed or not. If
Value=True, the parameter will print itself; if Value=False, it won't.

A FlagParameter can be initialized with three things. Prefix, Name (mandatory),
and Value (optional). The default for Value is False to indicate that the
parameter is off (i.e. not printed) by default. The only thing that counts for a
flag is whether its value evaluates to True or to False.

If a FlagParameter has to print itself, it checks first whether it is on or off
(Value=True or Value=False). If it is off, it will return the empty string. If
it is on it will return the combination of its prefix and name.

The methods isOn() and isOff() will return True or False depending on the Value
of the FlagParameter. These methods can be used to see whether the parameter
will be printed on the command line or not. With the methods on() and off() the
parameter can be turned on or off. These methods don't take a value, because a
flag can't have a value. Internally, they'll set parameter.Value to True or
False.

.. % Example can be removed (b/c is in section 2?)?

::

   >>> from cogent.app.parameters import FlagParameter
   >>> tree = FlagParameter(Prefix='-',Name='tree')
   >>> tree.isOn()
   False
   >>> print tree
   <BLANKLINE>
   >>> tree.on()
   >>> print tree
   -tree


ValuedParameter: cogent.app. parameters.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ValuedParameter also inherits from Parameter. In addition to setting all
attributes of the parameter during initialization, a default value is set. This
is a private property of a ValuedParameter and will be set to the value with
which the parameter.Value is initialized. The Default value is available for
inspection through parameter.Default. The default value should not be changed by
the user. With the method reset() the Value of the parameter will be reset to
the default value.

Like in FlagParameter, the value is used to control whether the parameter will
print itself or not. If the Value is None, the parameter is off and __str__ will
return the empty string. If the Value is anything else, the parameter will be
printed in full glory: prefix, name, value and optionally delimiter and quotes.
If IsPath is True, the value will be wrapped in double quotes when printed
allowing for spaces in paths.

The methods isOn() and isOff() can be used to check whether the parameter will
be printed or not. If parameter.Value is not None, the parameter is on and will
be printed. If parameter.Value is None the parameter is off and won't be
printed. By using the method on(value) the Value of the parameter is set to the
specified value. If you accidentally try to turn the parameter on with the value
None, an error will be raised. Calling off() will set the Value of the parameter
to None.

.. % Example can be removed (b/c is in section 2?)?

::

   >>> from cogent.app.parameters import ValuedParameter
   >>> temp = ValuedParameter(Prefix='-',Name='T',Delimiter="=")
   >>> temp.isOn()
   False
   >>> print temp
   <BLANKLINE>
   >>> temp.on(37)
   >>> print temp
   -T=37
   >>> temp_def = ValuedParameter(Prefix='-',Name='T',Value=100,Delimiter="=")
   >>> temp_def.Default
   100
   >>> print temp_def
   -T=100
   >>> temp_def.on(15)
   >>> print temp_def
   -T=15
   >>> temp_def.reset()
   >>> print temp_def
   -T=100


MixedParameter: cogent.app. parameters.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MixedParameter is a subclass of ValuedParameter, because they share many
features. A MixedParameter is a parameter that has an optional value; sometimes
it behaves like a FlagParameter, sometimes like a ValuedParameter. An example
is: '-d[0\ `\mid`\ 1\ `\mid`\ 2]'.

During initialization the Default value is set like in ValuedParameter. The
method reset() is available to reset the parameter value to the default.

This type of Parameter has the most complicated control over 'on' or 'off'. If
the Value is False, the parameter is off. If the Value is None, the parameter is
on, but behaves like a flag (only prefix and name will be printed), if the Value
is anything else, the parameter is on and behaves like a ValuedParameter.

The methods isOn() and isOff() have the same functionality as in the other
parameter types. When using on(val=None) it is optional to specify the value. If
a MixedParameter is turned on without a value it will behave like flag. When
turned on with a value, it will behave like a ValuedParameter. The method off()
sets the Value to False, which indicates that the parameter should not be
printed.

.. % Example can be removed (b/c is in section 2?)?

::

   >>> from cogent.app.parameters import MixedParameter
   >>> d = MixedParameter(Prefix='-',Name='d',Delimiter='')
   >>> d.isOff()
   True
   >>> d.on()
   >>> print d
   -d
   >>> d.on(2)
   >>> print d
   -d2


FilePath: cogent.app. parameters.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The FilePath object inherits from string, and should be used to wrap all strings
that represent paths. Examples inlcude::

   my_file.txt
   /path/to/my/file.txt
   /path/to/my/dir/

Wrapping paths in a FilePath object wraps the path in quotes when it should be,
for example when passed to a system call, and doesn't wrap it in quotes when it
shouldn't be, for example when performing operations on strings. The following
example illustrates how this fails with a simple string, but performs as it
should with a FilePath. In this example, p1 and p2 are simple strings, and p3,
p4, p5, and p6 are FilePath objects. Since the example path contains spaces, a
system call would not generate the desired result if the path is not wrapped in
quotes. The FilePath object will wrap it in quotes when it is cast to a string,
but will not wrap it in quotes when performing other string operations. The
string object, on the other hand, does not differentiate, and joining p1 and p2
results in quotes placed in the middle of the string. ::

   >>> p1 = '"/path to/"'
   >>> p2 = '"my_file.txt"'
   >>> print str(p1 + p2)
   "/path to/""my_file.txt"
   >>> from cogent.app.parameters import FilePath
   >>> p3 = FilePath("/path to/")
   >>> p4 = FilePath("my_file.txt")
   >>> print str(p3+p4)
   "/path to/my_file.txt"
   >>> p5 = FilePath("/path to/")
   >>> p6 = FilePath("my_file.txt")
   >>> print str(p5+p6)
   "/path to/my_file.txt"

The FilePath object is used by MixedParameter and ValuedParameter when their
IsPath attribute is set to True. This causes the Value attribute to be cast to a
FilePath object, and it is wrapped in quotes when used in a system call. The
_input_as_path input handler also casts the input to a FilePath object. In
general, if you are working with a string in an application controller that
represents a path to a file or directory, for example in a custom input handler,
that string should be cast to a FilePath object. Failure to do this will result
in errors if users pass a path that contains spaces.


.. _sec:parameters:

Parameters: cogent.app. parameters.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For most applications multiple parameters can be set. An application controller
should have a set of known parameters with optional default values. All
parameters of an application are organized in a dictionary where keys are
parameter identifiers (combination of prefix and name) and values are Parameter
objects.

Sometimes it might be hard to remember what the identifiers of the parameters in
a specific application are. Lets look at an example. Suppose the user knows that
the temperature can be set in three applications. In application *A* with '-T',
in application *B* with '\*temp', and in application *C* with '--t'. It is very
likely that he/she doesn't remember what parameter is used in which application.
If every application controller has an synonyms dictionary which maps 'temp' to
the identifier of the application specific parameter, the user can always lookup
the temperature with parameters['temp'].

To support the lookup of parameters by synonyms, the class Parameters is not a
simple dictionary, but a MappedDict (in cogent.util.misc). A MappedDict is a
dictionary that can apply some function to a lookup value, before it looks it up
in the dictionary. This function is called a mask. In the Parameters object the
mask allows users to look up parameters in the dictionary by synonyms.

The Parameters object uses the private function _find_synonym() internally to
determine by what key the parameter will be looked up. If the key, given by the
user, appears in the synonyms dictionary the key to use for the parameters
dictionary is looked up. Otherwise, it is assumed that the user used an existing
key in the parameters dictionary. ::

   >>> from cogent.app.parameters import FlagParameter
   >>> a = FlagParameter('-','a')
   >>> from cogent.app.parameters import  ValuedParameter
   >>> b = ValuedParameter('-','T',Value=37,Delimiter='=')
   >>> from cogent.app.parameters import MixedParameter
   >>> c = MixedParameter('-','d',Value=0)
   >>> params = {'-a':a,'-T':b,'-d':c}
   >>> synonyms = {'temp':'-T','distance':'-d'}
   >>> from cogent.app.parameters import Parameters
   >>> p = Parameters(params,synonyms)
   >>> p['-a'].isOn()
   False
   >>> print p['temp']
   -T=37
   >>> p['distance'].on(2)
   >>> print p['-d']
   -d2


Input
-----

CommandLineApplication subclasses are called using the __call__ method with a
single variable, data. This is the value passed to the application on the
command line. The value of data will differ based on the application you are
interfacing. Controlling for this without having to overwrite __call__ for every
CommandLineApplication is the purpose of the _input_handlers discussed in
Section :ref:`sec:commandlineapplication`.

Some examples of data that might be passed to CommandLineApplications are
strings, via the _input_as_string input handler, a list of lines that should be
written to file and then passed to the application, via the _input_as_lines
input handler, or a path to a file or directory, via the _input_as_path input
handler. To define the input handler that should be used, the class data
_input_handler should be set. If one of the default input handlers is not
applicable for a new CommandLineApplication, you will need to write a custom
input handler. See ``cogent.app.raxml.Raxml._input_as_seqs`` for an example of a
custom input handler.

For a discussion of the predefined input handlers, see Section
:ref:`sec:commandlineapplication`. For a discussion on defining custom input
handlers, see Section :ref:`sec:step2`.

.. % I think we're better off just pointing to the relevant discussion, rather than including there here.
.. % Is there a good example that we could put here? I hesitate to put a real example, b/c if we're going to
.. % rewrite this as executable documentation, it will fail on any system that didn't have the application
.. % we use as an example.


Output
------


ResultPath:  cogent.app.util.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ResultPath object is intended to hold the important information pertaining
to output files created by an application, namely the path to where the file can
be found, and whether the file was written or not. ResultPath is a very simple
container class. It has no methods aside from __init__().

To initialize a ResultPath object you must specify the path to the output file
by setting the Path parameter. This must be a string, and it is *highly*
recommended that this be an absolute path, though relative paths will also work
in many cases. You can also optionally specify a boolean value specifying
whether the file has been written or not by setting the IsWritten parameter.
This is True by default. ::

   >>> from cogent.app.util import ResultPath
   >>> rp = ResultPath(Path='/tmp/my_output.txt',IsWritten=True)
   >>> rp.Path
   '/tmp/my_output.txt'
   >>> rp.IsWritten
   True
   >>> rp
   <cogent.app.util.ResultPath object at ...


.. _sec:commandlineappresult:

CommandLine-AppResult: cogent.app.util.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CommandLineApp result object is a dictionary, and a container class for the
results from a run of a CommandLineApplication.  CommandLineAppResult has three
default items in it:

StdOut:
   The stdout from the application, a file object

StdErr:
   The stderr from the application, a file object (or None when suppressing stderr,
   more on this in CommandLineApplication)

ExitStatus:
   The exit status of the application, and int that is returned by the
   CommandLineApplication, usually indicating the success or failure of the run of
   the application.

The __init__() method takes several required parameters:

out:
   a file handle which will be assigned to 'StdOut' in the CommandLineAppResult

err:
   a file handle or None which will be assigned to 'StdErr' in the
   CommandLineAppResult

exit_status:
   an int or None which will be assigned to 'ExitStatus' in the
   CommandLineAppResult

result_paths:
   a dictionary of ResultPath objects keyed by the name by which you want to access
   the output. For each item in result_paths, an entry will be created in the
   CommandLineAppResult and if ResultPath.IsWritten is True, the file specified by
   ResultPath.Path will be opened. If the file can not be opened (due to not being
   found at the specified path, or inadequate read access) an ApplicationError will
   be raised.

Note that as a user, you will never instantiate a CommandLineAppResult, it is
taken care of upon calling the CommandLineApplication. You may at times however
be responsible for creating the result_paths dict, but more on this in
CommandLineApplication.

