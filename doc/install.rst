.. _quick-install:

Quick installation
==================

For systems with ``easy_install``
---------------------------------

For the list of dependencies see the :ref:`required` software list.

The following assumes you have ``easy_install`` on your system (this comes standard with new Mac's for instance), that you have administrator privileges and that you're connected to the internet. See below if you don't have ``easy_install``.

The key steps for the minimal install are:

1. Download the :download:`requirements file <../cogent-requirements.txt>`.

2. Install pip ::

    $ sudo easy_install -U pip

3. Use pip to download, build and install PyCogent plus the numpy dependency. ::

    $ DONT_USE_PYREX=1 sudo pip install -r path/to/cogent-requirements.txt

.. note:: The ``DONT_USE_PYREX=1`` statement is required if you have Pyrex installed due to a conflict between setuptools and later versions of Pyrex. If you don't have Pyrex, this will still work.

If the above fails to download PyCogent you can `download the tarball <http://sourceforge.net/projects/pycogent>`_ to your hard drive and replace the first line of the :download:`requirements file <../cogent-requirements.txt>` with the full path to the tarball, e.g. ``/Users/my_user_name/Downloads/cogent-1.4.tgz``.

Optional installs
^^^^^^^^^^^^^^^^^

To use the Ensembl querying code
""""""""""""""""""""""""""""""""

Add the following lines to the requirements file ::

    MySQL-python>=1.2.2
    SQLAlchemy>=0.5

.. note:: The MySQL-python module requires that you have MySQL installed.

To use the parallel capabilities
""""""""""""""""""""""""""""""""

Add the following to the requirements file ::

    mpi4py>=1.0

To build the documentation
""""""""""""""""""""""""""

Add the following to the requirements file ::

    Sphinx>=0.6

To use the development version of PyCogent
""""""""""""""""""""""""""""""""""""""""""

Just replace the first line of the requirements file with ``https://pycogent.svn.sourceforge.net/svnroot/pycogent/trunk``.

To use the graphics capabilities
""""""""""""""""""""""""""""""""

You need to install matplotlib_ (version 0.99+) to use the drawing code. However, compiling matplotlib can be a challenge. We therefore suggest you obtain a prebuilt binary for your platform from the matplotlib_ project page rather than modify the requirements file. For OSX, we suggest reading the following instructions on how to `compiling matplotlib`_.

.. _pip: http://pypi.python.org/pypi/pip
.. _matplotlib: http://matplotlib.sourceforge.net/
.. _`compiling matplotlib`: http://bioinformatics.anu.edu.au/groups/huttleylab/wiki/da9fe/Building_matplotlib_for_Snow_Leopard.html

.. todo::

    **FOR RELEASE:** update the tarball name for the version.

Installing ``easy_install``
---------------------------

If your system doesn't have ``easy_install``, you can execute the following::

    $ sudo curl http://peak.telecommunity.com/dist/ez_setup.py | python

or, if you are on a linux system that has a package manager, you may only need to do something like::

    $ sudo apt-get install python-setuptools

Use the approach to getting ``easy_install`` that best suites your system, then follow the (above) instructions for the ``pip`` based installation.
