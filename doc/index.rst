.. _contents:

####################################
Welcome to PyCogent's documentation!
####################################

**Contents:**

.. toctree::
   :maxdepth: 1

   README
   coding_guidelines
   examples/index
   licenses
   cookbook/index

.. todolist::

.. todo::
    
    **FOR RELEASE:** unlink cookbook until it's complete?

.. _quick-install:

******************
Quick installation
******************

We have just implemented support for installing PyCogent using a python package installer called pip_. This is the most most robust approach we've found for building all the dependencies and, better yet, requires almost no modifications of our code! I'm going to assume you have ``easy_install`` on your system (this comes standard with new Mac's for instance), that you have administrator privileges and that you're connected to the internet.

The key steps are then to:

#. Download the :download:`requirements file <../cogent-requirements.txt>`. Delete lines for the dependencies you're not interested in (or if one causes trouble during installation).
    .. note:: Numpy is an absolute requirement.
#. Install pip
    
    .. code-block:: guess
    
        $ sudo easy_install -U pip

#. Use pip to download PyCogent plus all dependencies, build and install them.

    .. code-block:: guess
    
        $ sudo pip install -r path/to/cogent-requirements.txt


If you want to grab the development version of PyCogent, just replace the first line of the requirements file with ``https://pycogent.svn.sourceforge.net/svnroot/pycogent/trunk``.

.. _pip: http://pypi.python.org/pypi/pip

.. todo::
    
    **FOR RELEASE:** update the cogent-requirements.txt file to point to stable release tarball

Citation
--------

If you use this software for published work please cite either -- `Knight et al., 2007, Genome Biol, 8, R171 <http://genomebiology.com/2007/8/8/R171>`_; or, `Butterfield et al., 2004, BMC Bioinformatics, 5, 1 <http://www.biomedcentral.com/1471-2105/5/1>`_.

Contacts
--------

If you find a bug or have feature requests, please post comments on the corresponding `tracker`_ page at sourceforge. If you have any questions please post on the appropriate projects `forums`_ page at sourceforge.

**************
For Developers
**************

Grabbing from the subversion repository
---------------------------------------

To grab PyCogent from the sourceforge subversion repository, do the following::

    $ svn co https://pycogent.svn.sourceforge.net/svnroot/pycogent/trunk PyCogent

Building/testing the documentation
----------------------------------

To build the documentation or ``doctest`` the contents, you'll need to install Sphinx. Assuming you have ``easy_install`` configured on your machine (and if not, download ``setuptools`` from http://pypi.python.org/pypi/setuptools and follow the instructions). Then grab Sphinx::

    $ sudo easy_install -U sphinx

Generating the html form of the documentation requires changing into the ``doc`` directory and executing a ``make`` command::

    $ cd path/to/PyCogent/doc
    $ make html
    ... # bunch of output
    Build finished. The HTML pages are in _build/html.

This prints a bunch of output to screen and creates a directory ``PyCogent/doc/_build`` and within that ``html``. The index file ``PyCogent/doc/_build/html/index.html`` is the root of the documentation.

One can also generate a pdf file, using the Sphinx latex generation capacity. This is slightly more involved. (It also requires that you have an installation of TeTex_.)

.. _TeTex: http://www.tug.org/texlive/

First generate the latex ::

    $ make latex
    ... # bunch of output
    Build finished; the LaTeX files are in _build/latex.
    Run `make all-pdf' or `make all-ps' in that directory to run these through (pdf)latex.

then change into the ``latex`` dir and build the pdf ::

    $ cd _build/latex
    $ make all-pdf

You can now open ``PyCogent.pdf``.

To actually test the documentation, you need to be in the ``doc`` directory and then execute another ``make`` command::

    $ cd path/to/PyCogent/doc
    $ make doctest

The results are in ``_build/doctest/output.txt``.

.. note:: The documentation does not test for presence of 3rd party dependencies (such as applications or python modules) like the PyCogent ``unittest`` test suite. If you don't have all the 3rd party applications installed you will see failures. At this point **no effort** is being expended to hide such failures.

Adding to the documentation
---------------------------

For new use-case examples. Look at any of the existing examples. The restructured text format is pretty easy to write (for overview see the Sphinx `rest overview`_). The conventions adopted by PyCogent are using heading levels to be consistent with the Python.org standard (taken from `Sphinx headings`_). They are

- # with overline, for parts
- \* with overline, for chapters
- =, for sections
- -, for subsections
- ^, for subsubsections
- ", for paragraphs
- +, added for sub-paragraphs (non-standard)

Create your file in the ``examples`` directory, giving it a ``.rst`` suffix. Link it into the documentation tree, adding a line into the ``examples/index.rst`` file.

Then test that it works (rather than testing the entire suite, you can use the convenience script within doc). For instance, the following is a single test of one file::

    $ cd path/to/PyCogent/doc
    $ python doctest_rsts.py examples/reverse_complement.rst

Add todo's into the rst files using the ``todo`` directive as in

::
    
    .. todo::
    
        some task

To see the list of todo's in the project, uncomment the line that sets ``todo_include_todos=True`` in ``doc/conf.py``, then cd into the ``doc/`` and make the html docs again. The todo's are listed on the main page.

.. warning:: Be sure to revert the conf.py file back to it's original state so you don't change everyone else's documentation too!

If you have any questions, contact gavin_.

.. _`rest overview`: http://sphinx.pocoo.org/rest.html
.. _`Sphinx headings`: http://sphinx.pocoo.org/rest.html#sections
.. _gavin: Gavin.Huttley@anu.edu.au
.. _tracker: http://sourceforge.net/tracker2/?group_id=186234
.. _forums: http://sourceforge.net/forum/?group_id=186234

******
Search
******

* :ref:`search`

