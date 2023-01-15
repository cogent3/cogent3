.. _get_app:

Getting an app
==============

The ``get_app()`` function is a top-level import. It allows you to create an app instance by passing the name of the app you want as a string along with the values for any positional or keyword arguments. (See :ref:`app_help`.) Below, I create an ``omit_degenerates`` app for a DNA moltype.

.. jupyter-execute::

    from cogent3 import get_app

    just_nucs = get_app("omit_degenerates", moltype="dna")
    just_nucs
