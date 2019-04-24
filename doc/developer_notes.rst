How to make a contribution
==========================

Look at the `PyCogent3 issues <https://bitbucket.org/pycogent3/cogent3/issues>`_. Pick something you think you can tackle which is not already assigned and have a go! (Following the process outlined below.)

For Developers
==============

You first need to `install and configure Mercurial <https://confluence.atlassian.com/get-started-with-bitbucket/install-and-set-up-mercurial-860009660.html>`_ on your machine.

After that, the process is:

#. `fork the PyCogent3 repository <https://confluence.atlassian.com/bitbucketserver/using-forks-in-bitbucket-server-776639958.html>`_
#. clone this fork to your local machine
#. follow the :ref:`dev-install` instructions to install
#. :ref:`run-tests` to make sure the install was correct
#. `create a new branch <https://confluence.atlassian.com/bitbucket/branching-a-repository-223217999.html#BranchingaRepository-CreateaMercurialbranch>`_
#. make your changes, and add tests to the test-suite
#. Keep your repository in sync with the upstream PyCogent3 repository
#. run the test suite
#. commit your changes and push to your Bitbucket repo
#. make a pull request

.. _dev-install:

Installing in developer mode
----------------------------

Assuming you are in either a `virtualenv` or a conda environment, do the following::

    $ cd <to local repo>
    $ pip install -e .

That will also build the library.

.. _run-tests:

Run the test suite
------------------

We use the `unittest` framework for testing. The `tests/` directory largely mirrors the structure of `cogent3`, so finding the place to put tests should be pretty straighforward.

To run the tests, either::

    $ cd <to local repo>
    $ ./run_tests

or, you can also run these tests directly within the `tests/` directory.::

    $ cd <to local repo>/tests
    $ python alltests.py

In both cases, allowed options are ``--output-ok`` and ``--verbose``.

Building/testing the documentation
----------------------------------

To build the documentation or ``doctest`` the contents, you'll need to install Sphinx::

    $ pip install sphinx

Generating the html form of the documentation requires changing into the ``doc`` directory and executing a ``make`` command::

    $ cd path/to/PyCogent3/doc
    $ make html
    ... # bunch of output
    Build finished. The HTML pages are in _build/html.

The index file czan be found in ``PyCogent3/doc/_build/html/index.html`` is the root of the documentation.

You can also generate a pdf file, using the Sphinx latex generation capacity. This is slightly more involved. (It also requires that you have an installation of TeTex_.)

.. _TeTex: http://www.tug.org/texlive/

First generate the latex ::

    $ make latex
    ... # bunch of output
    Build finished; the LaTeX files are in _build/latex.
    Run `make all-pdf' or `make all-ps' in that directory to run these through (pdf)latex.

then change into the ``latex`` dir and build the pdf ::

    $ cd _build/latex
    $ make all-pdf

You can now open ``PyCogent3.pdf``.

To actually test the documentation, you need to be in the ``doc`` directory and then execute another ``make`` command::

    $ cd path/to/PyCogent3/doc
    $ make doctest

The results are in ``_build/doctest/output.txt``.

Adding to the documentation
---------------------------

You can maximise the cogent3 user experience for yourself and others by contributing to the documentation. If you solve a problem that you think might prove useful to others then fork, add it into the documentation and do a pull request. If you can think of ways to improve the existing documents let us know via a `ticket <https://bitbucket.org/pycogent3/cogent3/issues>`_.

For guidance on adding documentation, look at any of the existing examples. The restructured text format is pretty easy to write (for overview see the Sphinx `rest overview`_). The conventions adopted by PyCogent3 are to use heading levels to be consistent with the Python.org standard (taken from `Sphinx headings`_). They are

- # with overline, for parts
- \* with overline, for chapters
- =, for sections
- -, for subsections
- ^, for subsubsections
- ", for paragraphs
- +, added for sub-paragraphs (non-standard)

If it's a use-case, create your file in the ``examples`` directory, giving it a ``.rst`` suffix. Link it into the documentation tree, adding a line into the ``examples/index.rst`` file. If it's something you think should be added into the cookbook, add it into the appropriate cookbook document.

The new documentation checklist
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Things you should check before committing your new document:

- Add a line at the beginning with yourself as author (``.. sectionauthor:: My Name``) so people can contact you with feedback.
- Add any data files used in your documentation under ``PyCogent3/doc/data/``
- Add a download link to those files to ``PyCogent3/doc/data_file_links.rst`` following the style employed in that file.
- Spellcheck!!
- Check what you wrote is valid restructured text by building the documents for both html and latex. If your document isn't connected into the table of contents, Sphinx will print a warning to screen.
- Check you have correctly marked up the content and that it looks OK. Make sure that python code and shell commands are correctly highlighted and that literals are marked up as literals. In particular, check the latex build since it is common for text to span beyond the page margins. If the latter happens, revise your document!
- Check that it works (rather than testing the entire suite, you can use the convenience script within doc). For instance, the following is a single test of one file::

   $ cd path/to/PyCogent3/doc
   $ python doctest_rsts.py examples/reverse_complement.rst

Adding TODOs
^^^^^^^^^^^^

Add todo's into the rst files using the ``todo`` directive as in

::

    .. todo::

        some task

To see the list of todo's in the project, uncomment the line that sets ``todo_include_todos=True`` in ``doc/conf.py``, then cd into the ``doc/`` and make the html docs again. The todo's are listed on the main page.

.. warning:: Be sure to revert the conf.py file back to it's original state so you don't accidentally commit the change as this affects everyone else's documentation too!

Developing C-extensions
-----------------------

Extensions for PyCogent3 should be written in `Cython <http://www.cython.org/>`_.

If you have any questions, contact Gavin_.

.. _`rest overview`: http://sphinx.pocoo.org/rest.html
.. _`Sphinx headings`: http://sphinx.pocoo.org/rest.html#sections
.. _Gavin: Gavin.Huttley@anu.edu.au
.. _PyCogent3: https://bitbucket.org/pycogent3/cogent3
