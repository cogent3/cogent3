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
   cookbook/index
   developer_notes
   licenses
   ChangeLog

.. todolist::

.. todo::

    **FOR RELEASE:** unlink cookbook until it's complete

********
Overview
********

PyCogent is a software library for genomic biology. It is a fully integrated and thoroughly tested framework for: controlling third-party applications; devising workflows; querying databases; conducting novel probabilistic analyses of biological sequences; and generating publication quality graphics. It is distinguished by many unique built-in capabilities (such as true codon alignment) and the frequent addition of entirely new methods for the analysis of genomic data.

Our primary goal is to provide a collection of rigourously validated tools for the manipulation and analysis of genome biology data sets. The project is routinely employed in numerous labs across the world and has provided essential capabilities for many high profile publications.

.. _quick-install:

******************
Quick installation
******************

For systems with ``easy_install``
=================================

The following assumes you have ``easy_install`` on your system (this comes standard with new Mac's for instance), that you have administrator privileges and that you're connected to the internet.

The key steps are then to:

1. Download the :download:`requirements file <../cogent-requirements.txt>`. Delete lines for the dependencies you're not interested in (or if one causes trouble during installation).

.. note:: Numpy is an absolute requirement.

2. Install pip ::

    $ sudo easy_install -U pip

3. Use pip to download, build and install PyCogent plus all dependencies. ::

    $ sudo pip install -r path/to/cogent-requirements.txt

If you want to grab the development version of PyCogent, just replace the first line of the requirements file with ``https://pycogent.svn.sourceforge.net/svnroot/pycogent/trunk``.

.. note:: Although there is a dependency on matplotlib_ for some of the drawing code, significant changes to the matplotlib API mean that those PyCogent capabilities are broken with recent matplotlib releases. We are waiting until matplotlib has been updated to use numpy 1.3, at which point we will be updating our support.

.. todo::

    **FOR RELEASE:** specify the supported matplotlib version

.. _pip: http://pypi.python.org/pypi/pip
.. _matplotlib: http://matplotlib.sourceforge.net/

.. todo::

    **FOR RELEASE:** update the cogent-requirements.txt file to point to stable release tarball

Installing ``easy_install``
===========================

If your system doesn't have ``easy_install``, then execute the following::

    $ sudo curl http://peak.telecommunity.com/dist/ez_setup.py | python

and following the instructions for the pip based installation.

********
Citation
********

If you use this software for published work please cite -- `Knight et al., 2007, Genome Biol, 8, R171 <http://genomebiology.com/2007/8/8/R171>`_.

*************************
Contacts and contributing
*************************

If you find a bug, have feature or documentation requests then please post comments on the corresponding `tracker`_ page at sourceforge. If you have any questions please post on the appropriate projects `forums`_ page at sourceforge. We appreciate your input!

.. _tracker: http://sourceforge.net/tracker2/?group_id=186234
.. _forums: http://sourceforge.net/forum/?group_id=186234

******
Search
******

* :ref:`search`

