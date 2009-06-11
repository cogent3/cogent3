.. _quick-install:

Quick installation
==================

For systems with ``easy_install``
---------------------------------

The following assumes you have ``easy_install`` on your system (this comes standard with new Mac's for instance), that you have administrator privileges and that you're connected to the internet.

The key steps are then to:

1. Download the :download:`requirements file <../cogent-requirements.txt>`. Delete lines for the dependencies you're not interested in (or if one causes trouble during installation).

.. note:: Numpy is an absolute requirement.

2. Install pip ::

    $ sudo easy_install -U pip

3. Use pip to download, build and install PyCogent plus all dependencies. ::

    $ sudo pip install -r path/to/cogent-requirements.txt

If you want to grab the development version of PyCogent, just replace the first line of the requirements file with ``https://pycogent.svn.sourceforge.net/svnroot/pycogent/trunk``.

.. note:: Although there is a dependency on matplotlib_ for some of the drawing code, significant changes to the matplotlib API mean that those PyCogent capabilities are broken with recent matplotlib releases. We are waiting until matplotlib has been updated to use numpy 1.3, at which point we will be updating our support. Until then, the code requires matplotlib version 0.87.6.

.. _pip: http://pypi.python.org/pypi/pip
.. _matplotlib: http://matplotlib.sourceforge.net/

.. todo::

    **FOR RELEASE:** update the cogent-requirements.txt file to point to stable release tarball

Installing ``easy_install``
---------------------------

If your system doesn't have ``easy_install``, then execute the following::

    $ sudo curl http://peak.telecommunity.com/dist/ez_setup.py | python

and following the instructions for the pip based installation.
