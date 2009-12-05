For Developers
==============

Anyone can contribute to the development of PyCogent, not just registered developers. If you figure out a solution to something using PyCogent that you'd like to share please post comments on the forums_ page at sourceforge and we'll look at including your tips! If you have ideas for improving the current documentation, then we want to here from you so please post on the tracker_ page!

.. _forums: http://sourceforge.net/forum/?group_id=186234
.. _tracker: http://sourceforge.net/tracker2/?group_id=186234

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

You can maximise the cogent user experience for yourself and others by contributing to the documentation. If you solve a problem that you think might prove useful to others then add it into the documentation. If you can think of ways to improve the existing documents let us know.

For guidance on adding documentation, look at any of the existing examples. The restructured text format is pretty easy to write (for overview see the Sphinx `rest overview`_). The conventions adopted by PyCogent are to use heading levels to be consistent with the Python.org standard (taken from `Sphinx headings`_). They are

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
- Add any data files used in your documentation under ``PyCogent/doc/data/``
- Add a download link to those files to ``PyCogent/doc/data_file_links.rst`` following the style employed in that file.
- Spellcheck!!
- Check what you wrote is valid restructured text by building the documents for both html and latex. If your document isn't connected into the table of contents, Sphinx will print a warning to screen.
- Check you have correctly marked up the content and that it looks OK. Make sure that python code and shell commands are correctly highlighted and that literals are marked up as literals. In particular, check the latex build since it is common for text to span beyond the page margins. If the latter happens, revise your document!
- Check that it works (rather than testing the entire suite, you can use the convenience script within doc). For instance, the following is a single test of one file::

   $ cd path/to/PyCogent/doc
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

Extensions for PyCogent should be written in `Cython <http://www.cython.org/>`_.

If you have any questions, contact Gavin_.

.. _`rest overview`: http://sphinx.pocoo.org/rest.html
.. _`Sphinx headings`: http://sphinx.pocoo.org/rest.html#sections
.. _Gavin: Gavin.Huttley@anu.edu.au
