.. _quick-install:

Quick installation
==================

PyCogent.app for OSX 10.6
-------------------------

`Download PyCogent.app <http://sourceforge.net/projects/pycogent/files/PyCogent.app/>`_. This native OSX app comes bundled with all required dependencies and is a download, decompress and go experience! It also implements the novel script form system that controls command line scripts via a form based input mechanism.

By virtual machine
------------------

One way to install PyCogent is to install the QIIME virtual machine using VirtualBox. The installation instructions can be found `here <http://qiime.sourceforge.net/install/virtual_box.html>`_.

Please, note that this is the only installation method supported for Windows and that natively Windows does not support gz files properly so to uncompress a gz file in Windows use `7-zip <http://www.7-zip.org/>`_.

For systems with ``easy_install``
---------------------------------

For the list of dependencies see the :ref:`required` software list.

The following assumes you have ``easy_install`` on your system (this comes standard with new Macs for instance), that you have administrator privileges and that you're connected to the internet. See below if you don't have ``easy_install``.

The key steps for the minimal install are:

1. Download the :download:`requirements file <../cogent-requirements.txt>`.

2. Install pip ::

    $ sudo easy_install -U pip

3. Use pip to download, build and install PyCogent plus the numpy dependency. ::

    $ DONT_USE_PYREX=1 sudo pip install -r path/to/cogent-requirements.txt

.. note:: The ``DONT_USE_PYREX=1`` statement is required if you have Pyrex installed due to a conflict between setuptools and later versions of Pyrex. If you don't have Pyrex, this will still work.

If the above fails to download PyCogent you can `download the tarball <http://sourceforge.net/projects/pycogent>`_ to your hard drive and replace the first line of the :download:`requirements file <../cogent-requirements.txt>` with the full path to the tarball, e.g. ``/Users/my_user_name/Downloads/cogent-1.5.2.tgz``.

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

You need to install matplotlib_ (version 1.1.0+) to use the drawing code. However, compiling matplotlib can be a challenge. We therefore suggest you obtain a prebuilt binary for your platform from the matplotlib_ project page rather than modify the requirements file. For OSX, we suggest reading the following instructions on `compiling matplotlib`_.

.. _pip: http://pypi.python.org/pypi/pip
.. _matplotlib: http://matplotlib.sourceforge.net/
.. _`compiling matplotlib`: http://sourceforge.net/projects/pycogent/forums/forum/651121/topic/5635916

.. todo::

    **FOR RELEASE:** update the tarball name for the version.

Installing ``easy_install``
---------------------------

If your system doesn't have ``easy_install``, you can execute the following::

    $ sudo curl http://peak.telecommunity.com/dist/ez_setup.py | python

or, if you are on a linux system that has a package manager, you may only need to do something like::

    $ sudo apt-get install python-setuptools

Use the approach to getting ``easy_install`` that best suites your system, then follow the (above) instructions for the ``pip`` based installation.
