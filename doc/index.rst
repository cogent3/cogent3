.. _contents:

#####################################
Welcome to Cogent3's documentation!
#####################################

**Contents**

.. toctree::
    :maxdepth: 1
    
    install
    coding_guidelines
    data_file_links
    examples/index
    cookbook/index
    developer_notes
    licenses
    ChangeLog

.. note:: All the code presented in this documentation is working, but does not reflect the current best-practice. The documentation will be significantly revised in the next couple of months.


.. todolist::

********************
WARNING & DISCLAIMER
********************

Be warned that Cogent3_ is a significantly changed library from the original `PyCogent <http://www.pycogent.org>`_. The renaming has been done to emphasise these differences and to make the project name and import statement consistent (``cogent`` was always the import name, originating in the ``pyevolve`` project from 2004!).

Most of the changes from PyCogent involved elimination of modules, using `black <https://github.com/psf/black>`_ and `isort <https://github.com/timothycrosley/isort>`_ for coding style, rationalisation of interfaces and the addition of new features. For instance, we have an experimental ``cogent3.app`` module (documentation still being written) that is intended to present a functional programming style interface to ``cogent3`` capabilities.

The rewrite has been a massive amount of work and unfortunately the changes to the API are only partially documented at present -- please see the `wiki pages <https://bitbucket.org/Cogent3/cogent3/wiki/Home>`_. We will endeavour to make that up-to-date when we can. But if you want to see what happened to something, search in the repository history.

*************
YOU CAN HELP!
*************

Posting Bugs
============

If you discover a bug please post a ticket at the issues_ page!

********
Overview
********

Cogent3 is a software library for genomic biology. It is a fully integrated and thoroughly tested framework for: conducting novel probabilistic analyses of biological sequence evolution; and generating publication quality graphics. It is distinguished by many unique built-in capabilities (such as true codon alignment) and the frequent addition of entirely new methods for the analysis of genomic data.

Our primary goal is to provide a collection of rigorously validated tools for the manipulation and analysis of genome biology data sets.

.. todo::

    update the literature cited

*******
Support
*******

We acknowledge the provision by `Wingware <http://wingware.com>`_ of free licenses for their professional IDE. Their commitment, over more than a decade, to supporting our development of Open Source software for science is greatly appreciated. (The port to Python 3 was made much easier by using Wing's refactor tools!)

********
Citation
********

If you use this software for published work please cite -- `Knight et al., 2007, Genome Biol, 8, R171 <http://genomebiology.com/2007/8/8/R171>`_.

******
Search
******

* :ref:`search`

.. _Cogent3: https://bitbucket.org/Cogent3/cogent3

*******
Support
*******

We acknowledge the provision by `Wingware <http://wingware.com>`_ of free licenses for their professional IDE. Their commitment, over more than a decade, to supporting our development of Open Source software for science is greatly appreciated. (The port of PyCogent to Python 3 was made much easier by using Wing's refactor tools!)

.. _issues: https://bitbucket.org/pycogent3/cogent3/issues/