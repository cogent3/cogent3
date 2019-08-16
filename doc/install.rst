.. _quick-install:

Quick installation
==================

Until the library is officially released, these installation instructions explicitly reference the `bitbucket repository <https://bitbucket.org/pycogent3/cogent3>`_.

Using ``pip``
-------------

Assuming you have `pip <https://pypi.python.org/pypi/pip/>`_ installed on your system, the key steps for the minimal install are:

1. Install numpy ::

    $ pip install numpy

2. Install PyCogent:
    
    a) If you have `mercurial <https://pypi.python.org/pypi/Mercurial/3.9.1>`_ installed::

        $ DONT_USE_CYTHON=1 pip install hg+https://bitbucket.org/pycogent3/cogent3

    b) If you don't have mercurial installed, `download a zip file <https://bitbucket.org/pycogent3/cogent3/downloads>`_ to your hard drive and pip install as::
    
        $ DONT_USE_CYTHON=1 pip install /path/to/downloaded/archive.zip

.. note:: Use the ``DONT_USE_CYTHON=1`` if you want to be sure and use the ``*.c`` files we generated. If you don't have Cython installed, it has no effect.


Optional installs
^^^^^^^^^^^^^^^^^

To use the drawing and parallel computing capabilities, you will need to download a zip archive as indicated above. Then do::

    $ pip install /path/to/downloaded/archive.zip[all]

Using conda
-----------

Support not yet in place, but it will be...

.. TODO Write conda instructions

.. todo::

    complete this after writing conda package
