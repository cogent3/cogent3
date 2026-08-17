.. _union_dict:

``UnionDict`` -- a dict with set like operations and keys as attributes
=======================================================================

This object combines the key-element storage of a ``dict`` with the union operation of a ``set()`` object. It is used in the ``cogent3.draw`` module, primarily for the ``figure`` and ``layout`` attributes.

Accessing elements of a ``UnionDict``
-------------------------------------

Keys in a ``UnionDict`` can be accessed like attributes

.. jupyter-execute::

    from cogent3.util.union_dict import UnionDict

    data = UnionDict(a=2, b={"c": 24, "d": [25]})
    data.a

.. jupyter-execute::

    data["a"]

.. jupyter-execute::

    data.b.d

Updating a ``UnionDict``
------------------------

If you use the ``|`` bitwise operator to compare two dicts and the left one is a ``UnionDict``, a union operation is performed.

.. jupyter-execute::

    from cogent3.util.union_dict import UnionDict

    data = UnionDict(a=2, b={"c": 24, "d": [25]})
    data.b |= {"d": 25}
    data.b

This can also be done using the ``union`` method.

.. jupyter-execute::

    data.b.union({"d": [25]})

.. jupyter-execute::

    data.b
    {"c": 24, "d": [25]}

Accessing a non-existent ``UnionDict`` key
------------------------------------------

.. jupyter-execute::
    :raises: KeyError

    from cogent3.util.union_dict import UnionDict

    data = UnionDict(a=2, b={"c": 24, "d": [25]})
    data["k"]

But if accessing as an attribute, you get an attribute error.

.. jupyter-execute::
    :raises: AttributeError

    data.k
