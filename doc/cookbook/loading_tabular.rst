Loading a csv file
^^^^^^^^^^^^^^^^^^

We load a tab separated data file using the ``load_table()`` function. The format is inferred from the filename suffix and you will note, in this case, it's not actually a `csv` file.

.. jupyter-execute::
    :linenos:

    from cogent3 import load_table

    table = load_table("data/stats.tsv")
    table

.. note:: The known filename suffixes are ``.csv``, ``.tsv`` and ``.pkl`` or ``.pickle`` (Python's pickle format).

Loading delimited specifying the format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Although unnecessary in this case, it's possible to override the suffix by specifying the delimiter using the ``sep`` argument.

.. jupyter-execute::
    :linenos:

    from cogent3 import load_table

    table = load_table("data/stats.tsv", sep="\t")
    table

Selectively loading parts of a big file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you only want a subset of the contents of a file, use the ``FilteringParser``. This allows skipping certain lines by using a callback function. We illustrate this with ``stats.tsv``, skipping any rows with ``"Ratio"`` > 10.

.. jupyter-execute::
    :linenos:

    from cogent3.parse.table import FilteringParser

    reader = FilteringParser(
        lambda line: float(line[2]) <= 10, with_header=True, sep="\t"
    )
    table = load_table("data/stats.tsv", reader=reader, digits=1)
    table

.. note:: You can also ``negate`` a condition, which is useful if the condition is complex.

Make a table from header and rows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::
    :linenos:

    from cogent3 import make_table

    header = ["A", "B", "C"]
    rows = [range(3), range(3, 6), range(6, 9), range(9, 12)]
    table = make_table(header=["A", "B", "C"], data=rows)
    table

Make a table from a ``dict``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a ``dict`` with key's as column headers.

.. jupyter-execute::
    :linenos:

    from cogent3 import make_table

    data = dict(A=[0, 3, 6], B=[1, 4, 7], C=[2, 5, 8])
    table = make_table(data=data)
    table

Specify the column order when creating from a ``dict``.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::
    :linenos:

    table = make_table(header=["C", "A", "B"], data=data)
    table

Create the table with an index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``Table`` can be indexed like a dict if you designate a column as the index (and that column has a unique value for every row).

.. jupyter-execute::

    table = load_table("data/stats.tsv", index="Locus")
    table["NP_055852"]

.. jupyter-execute::

    table["NP_055852", "Region"]

.. note:: The ``index`` argument also applies when using ``make_table()``.

Create a table from a ``pandas`` data frame
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::
    :linenos:

    from pandas import DataFrame

    df = DataFrame(data=[[0, 1], [3, 7]], columns=["a", "b"])
    table = make_table(data_frame=df)
    table
