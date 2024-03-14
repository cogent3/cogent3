.. _install:

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

The graphing dependency Plotly_ is installed too.

.. code-block:: bash

   $ pip install "cogent3[extra]"

Install with developer tools
============================

Everything we use for ``cogent3`` development.

.. code-block:: bash

   $ pip install "cogent3[dev]"

.. note:: Installs all dependencies that can be installed using ``pip``.

Installing the development version
==================================

.. code-block:: bash

   $ pip install git+https://github.com/cogent3/cogent3.git@develop#egg=cogent3

.. warning:: The interface can change without warning.

Installing using ``conda`` / ``mamba``
======================================

Activate your conda_ environment, then

.. code-block:: bash

    (myenv) $ conda install cogent3

.. _conda: https://docs.conda.io/en/latest/miniconda.html
.. _Plotly: https://pypi.org/project/plotly/
