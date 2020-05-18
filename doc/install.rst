************
Installation
************

Minimal installation
====================

Suitable for when you don't need plotting.

.. code-block:: bash

   $ pip install cogent3

Install with graphing tools
===========================

.. code-block:: bash

   $ pip install cogent3[extra]

.. note:: Writing image files requires you install ``plotly-orca``. That's best done using ``conda`` (see below).

Install with developer tools
============================

.. code-block:: bash

   $ pip install cogent3[dev]

.. note:: Installs all dependencies that can be installed using ``pip``.

Installing the development version
==================================

.. code-block:: bash

   $ pip install git+https://github.com/cogent3/cogent3.git@develop#egg=cogent3

.. warning:: This is the bleeding edge!

Installing using ``conda``
==========================

Strictly speaking, this is just a ``pip install`` into a conda_ environment. Do this is if you want all the capabilities that come with Jupyterlab_ and Plotly_, particularly the ability to save image files.

Manual creation of the ``conda`` environment
--------------------------------------------

Once you've downloaded and installed conda_, create a custom environment, e.g.

.. code-block:: bash

    $ conda create -n c3-env python=3
    $ conda activate c3-env

Once the environment is active

.. code-block:: bash

    (c3-env) $ conda install jupyterlab # or pip install jupyterlab if you prefer

then follow the instructions posted on the Plotly_ pipy page to install and configure Plotly for jupyter.

.. _conda: https://docs.conda.io/en/latest/miniconda.html
.. _Plotly: https://pypi.org/project/plotly/
.. _Jupyterlab: https://jupyter.org

Then install ``cogent3`` using the variant of ``pip install`` from above that you want.
