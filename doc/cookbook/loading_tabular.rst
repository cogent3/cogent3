.. jupyter-execute::
    :hide-code:

    import set_working_directory

Loading a csv file
==================

We load a tab separated data file using the ``load_table()`` function. The format is inferred from the filename suffix and you will note, in this case, it's not actually a `csv` file.

.. jupyter-execute::

    from cogent3 import load_table

    table = load_table("data/stats.tsv")
    table

.. note:: The known filename suffixes for reading are ``.csv``, ``.tsv`` and ``.pkl`` or ``.pickle`` (Python's pickle format).

.. note:: If you invoke the static column types argument, i.e.``load_table(..., static_column_types=True)`` and the column data are not static, those columns will be left as a string type.

Loading delimited specifying the format
=======================================

Although unnecessary in this case, it's possible to override the suffix by specifying the delimiter using the ``sep`` argument.

.. jupyter-execute::

    from cogent3 import load_table

    table = load_table("data/stats.tsv", sep="\t")
    table

Loading delimited data without a header line
============================================

To create a table from the follow examples, you specify your header and use ``make_table()``.

Using ``load_delimited()``
--------------------------

This is just a standard parsing function which does not do any filtering or converting elements to non-string types.

.. jupyter-execute::

    from cogent3.parse.table import load_delimited

    header, rows, title, legend = load_delimited("data/CerebellumDukeDNaseSeq.pk", header=False, sep="\t")
    rows[:4]

Using ``FilteringParser``
-------------------------

.. jupyter-execute::

    from cogent3.parse.table import FilteringParser
    
    reader = FilteringParser(with_header=False, sep="\t")
    rows = list(reader("data/CerebellumDukeDNaseSeq.pk"))
    rows[:4]

Selectively loading parts of a big file
=======================================

Loading a set number of lines from a file
-----------------------------------------

The ``limit`` argument specifies the number of lines to read.

.. jupyter-execute::

    from cogent3 import load_table

    table = load_table("data/stats.tsv", limit=2)
    table

Loading only some rows
----------------------

If you only want a subset of the contents of a file, use the ``FilteringParser``. This allows skipping certain lines by using a callback function. We illustrate this with ``stats.tsv``, skipping any rows with ``"Ratio"`` > 10.

.. jupyter-execute::

    from cogent3.parse.table import FilteringParser

    reader = FilteringParser(
        lambda line: float(line[2]) <= 10, with_header=True, sep="\t"
    )
    table = load_table("data/stats.tsv", reader=reader, digits=1)
    table

You can also ``negate`` a condition, which is useful if the condition is complex. In this example, it means keep the rows for which ``Ratio > 10``.

.. jupyter-execute::

    reader = FilteringParser(
        lambda line: float(line[2]) <= 10, with_header=True, sep="\t", negate=True
    )
    table = load_table("data/stats.tsv", reader=reader, digits=1)
    table

Loading only some columns
-------------------------

Specify the columns by their names.

.. jupyter-execute::

    from cogent3.parse.table import FilteringParser

    reader = FilteringParser(columns=["Locus", "Ratio"], with_header=True, sep="\t")
    table = load_table("data/stats.tsv", reader=reader)
    table

Or, by their index.

.. jupyter-execute::

    from cogent3.parse.table import FilteringParser

    reader = FilteringParser(columns=[0, -1], with_header=True, sep="\t")
    table = load_table("data/stats.tsv", reader=reader)
    table

.. note:: The ``negate`` argument does not affect the columns evaluated.

Load raw data as a list of lists of strings
-------------------------------------------

We just use ``FilteringParser``.

.. jupyter-execute::

    from cogent3.parse.table import FilteringParser

    reader = FilteringParser(with_header=True, sep="\t")
    data = list(reader("data/stats.tsv"))

We just display the first two lines.

.. jupyter-execute::

    data[:2]

.. note:: The individual elements are all ``str``.

Make a table from header and rows
=================================

.. jupyter-execute::

    from cogent3 import make_table

    header = ["A", "B", "C"]
    rows = [range(3), range(3, 6), range(6, 9), range(9, 12)]
    table = make_table(header=["A", "B", "C"], data=rows)
    table

Make a table from a ``dict``
============================

For a ``dict`` with key's as column headers.

.. jupyter-execute::

    from cogent3 import make_table

    data = dict(A=[0, 3, 6], B=[1, 4, 7], C=[2, 5, 8])
    table = make_table(data=data)
    table

Specify the column order when creating from a ``dict``.
=======================================================

.. jupyter-execute::

    table = make_table(header=["C", "A", "B"], data=data)
    table

Create the table with an index
==============================

A ``Table`` can be indexed like a dict if you designate a column as the index (and that column has a unique value for every row).

.. jupyter-execute::

    table = load_table("data/stats.tsv", index_name="Locus")
    table["NP_055852"]

.. jupyter-execute::

    table["NP_055852", "Region"]

.. note:: The ``index`` argument also applies when using ``make_table()``.

Create a table from a ``pandas.DataFrame``
==========================================

.. jupyter-execute::

    from pandas import DataFrame

    from cogent3 import make_table

    data = dict(a=[0, 3], b=["a", "c"])
    df = DataFrame(data=data)
    table = make_table(data_frame=df)
    table

Create a table from header and rows
===================================

.. jupyter-execute::

    from cogent3 import make_table

    table = make_table(header=["a", "b"], data=[[0, "a"], [3, "c"]])
    table

Create a table from dict
========================

``make_table()`` is the utility function for creating ``Table`` objects from standard python objects.

.. jupyter-execute::

    from cogent3 import make_table

    data = dict(a=[0, 3], b=["a", "c"])
    table = make_table(data=data)
    table

Create a table from a 2D dict
=============================

.. jupyter-execute::

    from cogent3 import make_table

    d2D = {
        "edge.parent": {
            "NineBande": "root",
            "edge.1": "root",
            "DogFaced": "root",
            "Human": "edge.0",
        },
        "x": {
            "NineBande": 1.0,
            "edge.1": 1.0,
            "DogFaced": 1.0,
            "Human": 1.0,
        },
        "length": {
            "NineBande": 4.0,
            "edge.1": 4.0,
            "DogFaced": 4.0,
            "Human": 4.0,
        },
    }
    table = make_table(
        data=d2D,
    )
    table

Create a table that has complex python objects as elements
==========================================================

.. jupyter-execute::

    from cogent3 import make_table

    table = make_table(
        header=["abcd", "data"],
        data=[[range(1, 6), "0"], ["x", 5.0], ["y", None]],
        missing_data="*",
        digits=1,
    )
    table

Create an empty table
=====================

.. jupyter-execute::

    from cogent3 import make_table

    table = make_table()
    table
