The Readme
==========

:Download: `From sourceforge <http://sourceforge.net/projects/pycogent>`_ or follow the :ref:`quick-install` instructions.
:Registration: To be informed of bugs, new releases please subscribe to the `mailing lists at sourceforge <http://sourceforge.net/projects/pycogent>`_.

Dependencies
------------

The toolkit requires Python 2.5.1 or greater, and Numpy 1.3 or greater. Aside from these the dependencies below are optional and the code will work as is. A C compiler, however, will allow external C module's responsible for the likelihood and matrix exponentiation calculations to be compiled, resulting in significantly improved performance.

.. _required:

Required
^^^^^^^^

- Python_: the language the toolkit is primarily written in, and in which the user writes control scripts.
- Numpy_: This is a python module used for speeding up matrix computations. It is available as source code for \*nix.
- zlib_: This is a compression library which is available for all platforms and comes pre-installed on most too. If, by chance, your platform doesn't have this installed then download the source from the zlib_ site and follow the install instructions, or refer to the instructions for `compiling matplotlib`_.

.. note:: On some linux platforms (like Ubuntu), you must specifically install a ``python-dev`` package so that the Python_ header files required for building some external dependencies are available.

Optional
^^^^^^^^

- C compiler: This is standard on most \*nix platforms. On Macos X this is available for free in the Developer tools which, if you don't already have them, can be obtained from Apple_.
- Matplotlib_: used to plot several kinds of graphs related to codon usage. For installation, see these instructions for `compiling matplotlib`_.
- Cython_: This module is only necessary if you are a developer who wants to modify the \*.pyx files.
- mpi4py_: Message Passing Interface interface, required for parallel computation.
- SQLAlchemy_ and `MySQL-python`_: These are required for the Ensembl querying code.

If you use the :ref:`quick-install` approach, these are all specified in the requirements file.

Installation
------------

If you don't wish to use the :ref:`quick-install` approach then the conventional \*nix platform (including MacOS X) python package installation procedure applies. Download the software from `here <http://sourceforge.net/projects/pycogent>`_. Uncompress the archive and change into the ``PyCogent`` directory and type::

$ python setup.py build

This automatically compiles the modules. If you have administrative privileges type::

$ sudo python setup.py install

This then places the entire package into your python/site-packages folder.

If you do not have administrator privileges on your machine you can change the build approach to::

$ python setup.py build_ext -if

which compiles the extensions in place (the ``i`` option) forcibly (the ``f`` option, ie even if they've already been compiled). Then move the cogent directory to where you want it (or leave it in place) and add this location to your python path using ``sys.path.insert(0, "/your/path/to/PyCogent")`` in each script, or by setting shell environment variables (e.g. ``$ export PYTHONPATH=/your/path/to/PyCogent:$PYTHONPATH``)

Testing
-------

``PyCogent/tests`` contains all the tests (currently >3100). You can most readily run the tests using the ``PyCogent/run_tests`` shell script. This is done by typing:

.. code-block:: guess
    
    $ sh run_tests

which will automatically build extensions in place, set up the PYTHONPATH and run ``PyCogent/tests/alltests.py``. Note that if certain optional applications are not installed this will be indicated in the output as "can't find" or "not installed". A "`.`" will be printed to screen for each test and if they all pass, you'll see output like:

.. code-block:: guess
    
    Ran 3299 tests in 58.128s
    
    OK

Tips for usage
--------------

A good IDE can greatly simplify writing control scripts. Features such as code completion and definition look-up are extremely useful. For a complete list of `editors go here`_.

To get help on attributes of an object in python, use

.. code-block:: python
    
    >>> dir(myalign)

to list the attributes of ``myalign`` or 

.. code-block:: python
    
    >>> help(myalign.writeToFile)

to figure out how to use the ``myalign.writeToFile`` method. Also note that the directory structure of the package is similar to the import statements required to use a module -- to see the contents of ``alignment.py`` or ``sequence.py`` you need to look in the directory ``cogent/core`` path, to use the classes in those files you specify ``cogent.core`` for importing.

.. _Python: http://www.python.org
.. _Cython: http://www.cython.org/
.. _Numpy: http://numpy.scipy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net
.. _Apple: http://www.apple.com
.. _Pyrex: http://www.cosc.canterbury.ac.nz/~greg/python/Pyrex/
.. _`editors go here`: http://www.python.org/cgi-bin/moinmoin/PythonEditors
.. _mpi4py: http://code.google.com/p/mpi4py
.. _`restructured text`: http://docutils.sourceforge.net/rst.html
.. _gcc: http://gcc.gnu.org/
.. _SQLAlchemy: http://www.sqlalchemy.org
.. _`MySQL-python`: http://sourceforge.net/projects/mysql-python
.. _zlib: http://www.zlib.net/
.. _`compiling matplotlib`: http://bioinformatics.anu.edu.au/groups/huttleylab/wiki/da9fe/Building_matplotlib_for_Snow_Leopard.html
