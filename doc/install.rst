.. _quick-install:

Quick installation
==================

Using ``pip``
-------------

The key steps for the minimal install are:

1. Install pip ::

    $ sudo easy_install -U pip

2. Use pip to download, build and install PyCogent plus the numpy dependency. ::

    $ DONT_USE_PYREX=1 sudo pip install -r path/to/cogent-requirements.txt

.. note:: The ``DONT_USE_PYREX=1`` statement is required if you have Pyrex installed due to a conflict between setuptools and later versions of Pyrex. If you don't have Pyrex, this will still work.

If the above fails to download PyCogent you can `download the tarball <https://bitbucket.org/pycogent3/pycogent3>`_ to your hard drive and pip install as.

::

    $ pip installs /path/to/downloaded/archive.zip

Optional installs
^^^^^^^^^^^^^^^^^

To use the drawing and parallel computing capabilities do::

    $ pip install cogent3[all]

Using conda
-----------

Support not yet in place, but it will be...

.. TODO Write conda instructions
