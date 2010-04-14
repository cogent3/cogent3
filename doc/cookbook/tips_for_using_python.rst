.. authors: Doug Wendel

*********************
Tips for using python
*********************

If you are new to Python, this is the right place for you. This is not a comprehensive guide, rather this section is intended help you install an appropriate version of Python and get you started using the language with PyCogent.

Checking you version
====================

Before starting, it's recommended that you update your version of python to the latest 2.x stable build. At a minimum, you need to be using Python 2.5.x. If you have Python installed, open a terminal window and type "python". The first few lines of text should tell you what version you're running. For example: 

| ``$ python``
| ``Python 2.6.4 (r264:75706, Mar 11 2010, 12:48:01)``
| ``[GCC 4.2.1 (Apple Inc. build 5646) (dot 1)] on darwin``
| ``Type "help", "copyright", "credits" or "license" for more information.``

In this case, Python 2.6.4 is installed. If you don't have Python installed or your version is older than 2.5.x, download the most recent 2.x release of Python here:

http://www.python.org/download/

**DO NOT** install Python 3.x or above. This version of the language was significantly restructured and will not work with PyCogent.

Getting help
============

Python comes with a built-in help utility. To access it, open a terminal window and type ``python`` to enter Python's interactive mode:

| ``>>> help()``

| ``Welcome to Python 2.6!  This is the online help utility.``

| ``If this is your first time using Python, you should definitely check out``
| ``the tutorial on the Internet at http://docs.python.org/tutorial/.``

| ``Enter the name of any module, keyword, or topic to get help on writing``
| ``Python programs and using Python modules.  To quit this help utility and``
| ``return to the interpreter, just type "quit".``

| ``To get a list of available modules, keywords, or topics, type "modules",``
| ``"keywords", or "topics".  Each module also comes with a one-line summary``
| ``of what it does; to list the modules whose summaries contain a given word``
| ``such as "spam", type "modules spam".``

| ``help>``

Note that you are now in the interactive help mode as noted by the prompt ``help>``. In this mode you can type in the names of various objects and get additional information on them. For example, if I want to know something about the ``map`` function:

| ``help> map``

| ``Help on built-in function map in module __builtin__:``

| ``map(...)``
| ``map(function, sequence[, sequence, ...]) -> list``

| ``Return a list of the results of applying the function to the items of``
| ``the argument sequence(s).  If more than one sequence is given, the``
| ``function is called with an argument list consisting of the corresponding``
| ``item of each sequence, substituting None for missing values when not all``
| ``sequences have the same length.  If the function is None, return a list of``
| ``the items of the sequence (or a list of tuples if more than one sequence).``
| ``(END)``
	
To exit the interactive help mode, simply enter a blank line at the ``help>`` prompt:

| ``help>``

| ``You are now leaving help and returning to the Python interpreter.``
| ``If you want to ask for help on a particular object directly from the``
| ``interpreter, you can type "help(object)".  Executing "help('string')"``
| ``has the same effect as typing a particular string at the help> prompt.``
	
As the parting message suggests, you can also invoke help on a specific object directly:

| ``>>> help(abs)``

| ``Help on built-in function abs in module __builtin__:``

| ``abs(...)``
| ``abs(number) -> number``

| ``Return the absolute value of the argument.``
| ``(END)``

Using the dir() function
========================

Another useful built-in function is ``dir()``. As the name implies, it's use is to list defined names in the current scope. To list the currently defined names:

| ``>>> dir()``
| ``['__builtins__', '__doc__', '__name__', '__package__']``
	
The list shows which names are currently defined. This list includes all imported modules and variable names. For example, if I define a new variable, it will also show up in this list:

| ``>>> my_variable = 'Just testing'``
| ``>>> dir()``
| ``['__builtins__', '__doc__', '__name__', '__package__', 'my_variable']``
	
Imported modules will also be reflected in this list:

| ``>>> import os``
| ``>>> import sys``
| ``>>> dir()``
| ``['__builtins__', '__doc__', '__name__', '__package__', 'my_variable', 'os', 'sys']``

``dir()`` can also be used to list the names defined within a module:

| ``>>> import sys``
| ``>>> dir(sys)``
| ``['__displayhook__', '__doc__', '__excepthook__', '__name__', '__package__', '__stderr__', '__stdin__', '__stdout__',...``

It also works on variable types. For example, let's see what attributes the string class has as defined:

| ``>>> dir(str)``
| ``['__add__', '__class__', '__contains__', '__delattr__', '__doc__', '__eq__', '__format__', '__ge__',...``

You can also use ``dir()`` on a defined variable. It will inspect the variable's type and report the attributes for that type. In this case, we defined a variable ``my_variable`` of type ``str``. Calling ``dir(my_variable)`` will product the same result as calling ``dir(str)``:

| ``>>> my_variable = 'Just testing'``
| ``>>> dir(my_variable)``
| ``['__add__', '__class__', '__contains__', '__delattr__', '__doc__', '__eq__', '__format__', '__ge__',...``
	
Hello PyCogent!
===============

Now that we've gotten our feet wet, let's write a simple function that returns a friendly message. This is a simple function which takes in one parameter, ``your_name``, and outputs the user's name prefixed with a standard message. Calling your new function is as simple as typing the name of the function and supplying the appropriate variables::

	>>> def hello_pycogent(your_name):
	...     message = 'PyCogent bids you welcome ' + your_name
	...     print message
	... 
	>>> hello_pycogent('John Smith')
	PyCogent bids you welcome John Smith
	
Enter each line as you see it and note that white space is important! There are no brackets or keywords to signal blocks of code. Instead, indentation is used to designate related lines of code.

Further Python documentation
============================

Now that you've got Python up and running and know a few commands, it might be useful to browse the official documentation. There is a comprehensive list of information and some excellent tutorials to work though:

http://docs.python.org/

There are also many code examples to be found in the Python cookbook:

http://code.activestate.com/recipes/langs/python/