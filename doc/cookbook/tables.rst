.. jupyter-execute::
    :hide-code:

    import set_working_directory

************
Tabular data
************

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul

``Table`` handles tabular data, storing as columns in a, you guessed it, ``columns`` attribute. The latter acts like a dictionary, with the column names as the keys and the column values being  ``numpy.ndarray`` instances. The table itself is iterable over rows.

.. note:: ``Table`` is immutable at the level of the individual ``ndarray`` not being writable.

.. include:: ./loading_tabular.rst

Adding a new column
===================

.. jupyter-execute::

    from cogent3 import make_table

    table = make_table()
    table.columns["a"] = [1, 3, 5]
    table.columns["b"] = [2, 4, 6]
    table

Add a title and a legend to a table
===================================

This can be done when you create the table.

.. jupyter-execute::

    from cogent3 import make_table

    data = dict(a=[0, 3], b=["a", "c"])
    table = make_table(data=data, title="Sample title", legend="a legend")
    table

It can be done by directly assigning to the corresponding attributes.

.. jupyter-execute::

    data = dict(a=[0, 3], b=["a", "c"])
    table = make_table(data=data)
    table.title = "My title"
    table

Iterating over table rows
=========================

``Table`` is a row oriented object. Iterating on the table returns each row as a new ``Table`` instance.

.. jupyter-execute::

    from cogent3 import load_table

    table = load_table("data/stats.tsv")
    for row in table:
        print(row)
        break

The resulting rows can be indexed using their column names.

.. jupyter-execute::

    for row in table:
        print(row["Locus"])

How many rows are there?
========================

The ``Table.shape`` attribute is like that of a ``numpy`` ``array``. The first element (``Table.shape[0]``) is the number of rows.

.. jupyter-execute::

    from cogent3 import make_table

    data = dict(a=[0, 3, 5], b=["a", "c", "d"])
    table = make_table(data=data)
    table.shape[0] == 3

How many columns are there?
===========================

``Table.shape[1]`` is the number of columns. Using the table from above.

.. jupyter-execute::

    table.shape[1] == 2

Iterating over table columns
============================

The ``Table.columns`` attribute is a ``Columns`` instance, an object with ``dict`` attributes.

.. jupyter-execute::

    from cogent3 import load_table

    table = load_table("data/stats.tsv")
    table.columns

.. jupyter-execute::

    table.columns["Region"]

So iteration is the same as for dicts.

.. jupyter-execute::

    for name in table.columns:
        print(name)

Table slicing using column names
================================

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    table

Slice using the column name.

.. jupyter-execute::

    table[:2, "Region":]

Table slicing using indices
===========================

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    table[:2, :1]

Changing displayed numerical precision
======================================

We change the ``Ratio`` column to using scientific notation.

.. jupyter-execute::

    from cogent3 import load_table

    table = load_table("data/stats.tsv")
    table.format_column("Ratio", "%.1e")
    table

Change digits or column spacing
===============================

This can be done on table loading,

.. jupyter-execute::

    table = load_table("data/stats.tsv", digits=1, space=2)
    table

or, for spacing at least, by modifying the attributes

.. jupyter-execute::

    table.space = "    "
    table

Wrapping tables for display
===========================

Wrapping generates neat looking tables whether or not you index the table rows. We demonstrate here

.. jupyter-execute::

    from cogent3 import make_table

    h = ["name", "A/C", "A/G", "A/T", "C/A"]
    rows = [["tardigrade", 0.0425, 0.1424, 0.0226, 0.0391]]
    wrap_table = make_table(header=h, data=rows, max_width=30)
    wrap_table

.. jupyter-execute::

    wrap_table = make_table(header=h, data=rows, max_width=30, index_name="name")
    wrap_table

Display the top of a table using ``head()``
===========================================

.. jupyter-execute::

    table = make_table(data=dict(a=list(range(10)), b=list(range(10))))
    table.head()

You change how many rows are displayed.

.. jupyter-execute::

    table.head(2)

The table shape is that of the original table.

Display the bottom of a table using ``tail()``
==============================================

.. jupyter-execute::

    table.tail()

You change how many rows are displayed.

.. jupyter-execute::

    table.tail(1)

Display random rows from a table
================================

.. jupyter-execute::

    table.set_repr_policy(random=3)
    table

Change the number of rows displayed by ``repr()``
=================================================

.. jupyter-execute::

    table.set_repr_policy(head=2, tail=3)
    table

.. note:: The ``...`` indicates the break between the top and bottom rows.

Changing column headings
========================

The table ``header`` is immutable. Changing column headings is done as follows.

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    print(table.header)
    table = table.with_new_header("Ratio", "Stat")
    print(table.header)

Adding a new column
===================

.. jupyter-execute::

    from cogent3 import make_table

    table = make_table()
    table

.. jupyter-execute::

    table.columns["a"] = [1, 3, 5]
    table.columns["b"] = [2, 4, 6]
    table

Create a new column from existing ones
======================================

This can be used to take a single, or multiple columns and generate a new column of values. Here we'll take 2 columns and return True/False based on a condition.

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    table = table.with_new_column(
        "LargeCon",
        lambda r_v: r_v[0] == "Con" and r_v[1] > 10.0,
        columns=["Region", "Ratio"],
    )
    table

Get table data as a numpy array
===============================

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    table.array

Get a table column as a list
============================

Via the ``Table.tolist()`` method.

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    locus = table.tolist("Locus")
    locus

Or directly from the column array object.

.. jupyter-execute::

    table.columns["Locus"].tolist()

Get multiple table columns as a list
====================================

This returns a row oriented list.

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    rows = table.tolist(["Region", "Locus"])
    rows

.. note:: column name order dictates the element order per row

Get the table as a row oriented ``dict``
========================================

Keys in the resulting dict are the row indices, the value is a dict of column name, value pairs.

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    table.to_dict()

Get the table as a column oriented ``dict``
===========================================

Keys in the resulting dict are the column names, the value is a list.

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    table.columns.to_dict()

Get the table as a ``pandas.DataFrame``
=======================================

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    df = table.to_dataframe()
    df

You can also specify column(s) are categories

.. jupyter-execute::

    df = table.to_dataframe(categories="Region")

Get a table of counts as a contingency table
============================================

If our table consists of counts data, the ``Table`` can convert it into a ``CategoryCount`` instance that can be used for performing basic contingency table statistical tests, e.g. chisquare, G-test of independence, etc.. To do this, we must specify which column contains the row names using the ``index_name`` argument.

.. jupyter-execute::

    table = make_table(data={"Ts": [31, 58], "Tv": [36, 138], "": ["syn", "nsyn"]}, index_name="")
    table

.. jupyter-execute::

    contingency = table.to_categorical(["Ts", "Tv"])
    contingency

.. jupyter-execute::

    g_test = contingency.G_independence()
    g_test

Appending tables
================

.. warning:: Only for tables with the same columns.

Can be done without specifying a new column (set the first argument to ``appended`` to be ``None``). Here we simply use the same table data.

.. jupyter-execute::

    table1 = load_table("data/stats.tsv")
    table2 = load_table("data/stats.tsv")
    table = table1.appended(None, table2)
    table

Specifying with a new column. In this case, the value of the ``table.title`` becomes the value for the new column.

.. jupyter-execute::

    table1.title = "Data1"
    table2.title = "Data2"
    table = table1.appended("Data#", table2, title="")
    table

.. note:: We assigned an empty string to ``title``, otherwise the resulting table has the same ``title`` attribute as that of ``table1``.

Summing a single column
=======================

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    table.summed("Ratio")

Because each column is just a ``numpy.ndarray``, this also can be done directly via the array methods.

.. jupyter-execute::

    table.columns["Ratio"].sum()

Summing multiple columns or rows - strictly numerical data
==========================================================

We define a strictly numerical table,

.. jupyter-execute::

    from cogent3 import make_table

    all_numeric = make_table(
        header=["A", "B", "C"], data=[range(3), range(3, 6), range(6, 9), range(9, 12)]
    )
    all_numeric

and sum all columns (default condition)

.. jupyter-execute::

    all_numeric.summed()

and all rows

.. jupyter-execute::

    all_numeric.summed(col_sum=False)

Summing multiple columns or rows with mixed non-numeric/numeric data
====================================================================

We define a table with mixed data, like a distance matrix.

.. jupyter-execute::

    mixed = make_table(
        header=["A", "B", "C"], data=[["*", 1, 2], [3, "*", 5], [6, 7, "*"]]
    )
    mixed

and sum all columns (default condition), ignoring non-numerical data

.. jupyter-execute::

    mixed.summed(strict=False)

and all rows

.. jupyter-execute::

    mixed.summed(col_sum=False, strict=False)

Filtering table rows
====================

We can do this by providing a reference to an external function

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    sub_table = table.filtered(lambda x: x < 10.0, columns="Ratio")
    sub_table

or using valid python syntax within a string, which is executed

.. jupyter-execute::

    sub_table = table.filtered("Ratio < 10.0")
    sub_table

You can also filter for values in multiple columns

.. jupyter-execute::

    sub_table = table.filtered("Ratio < 10.0 and Region == 'NonCon'")
    sub_table

Filtering table columns
=======================

We select only columns that have a sum > 20 from the ``all_numeric`` table constructed above.

.. jupyter-execute::

    big_numeric = all_numeric.filtered_by_column(lambda x: sum(x) > 20)
    big_numeric

Standard sorting
================

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    table.sorted(columns="Ratio")

Reverse sorting
===============

.. jupyter-execute::

    table.sorted(columns="Ratio", reverse="Ratio")

Sorting involving multiple columns, one reversed
================================================

.. jupyter-execute::

    table.sorted(columns=["Region", "Ratio"], reverse="Ratio")

Getting raw data for a single column
====================================

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    raw = table.tolist("Region")
    raw

Getting raw data for multiple columns
=====================================

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    raw = table.tolist(["Locus", "Region"])
    raw

Getting distinct values
=======================

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    assert table.distinct_values("Region") == set(["NonCon", "Con"])

Counting occurrences of values
==============================

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    assert table.count("Region == 'NonCon' and Ratio > 1") == 1

Counting unique values
======================

This returns a ``CategoryCounter``, a dict like class.

.. jupyter-execute::

    from cogent3 import make_table

    table = make_table(
        data=dict(A=["a", "b", "b", "b", "a"], B=["c", "c", "c", "c", "d"])
    )
    unique = table.count_unique("A")
    type(unique)

.. jupyter-execute::

    unique

For multiple columns.

.. jupyter-execute::

    unique = table.count_unique(["A", "B"])
    unique

.. jupyter-execute::

    r = unique.to_table()
    r

Joining or merging tables
=========================

We do a standard inner join here for a restricted subset. We must specify the columns that will be used for the join. Here we just use ``Locus``.

.. jupyter-execute::

    rows = [
        ["NP_004893", True],
        ["NP_005079", True],
        ["NP_005500", False],
        ["NP_055852", False],
    ]
    region_type = make_table(header=["Locus", "LargeCon"], data=rows)
    stats_table = load_table("data/stats.tsv")
    new = stats_table.joined(region_type, columns_self="Locus")
    new

.. note:: If the tables have titles, column names are prefixed with those instead of ``right_``.

.. note:: The ``joined()`` method is just a wrapper for the ``inner_join()`` and ``cross_join()`` (row cartesian product) methods, which you can use directly.

Transpose a table
=================

.. jupyter-execute::

    from cogent3 import make_table

    header = ["#OTU ID", "14SK041", "14SK802"]
    rows = [
        [-2920, "332", 294],
        [-1606, "302", 229],
        [-393, 141, 125],
        [-2109, 138, 120],
    ]
    table = make_table(header=header, rows=rows)
    table

We require a new column heading for the current header data. We also need to specify which existing column will become the header.

.. jupyter-execute::

    tp = table.transposed(new_column_name="sample", select_as_header="#OTU ID")
    tp

Specify markdown as the ``str()`` format
========================================

Using the method provides finer control over formatting.

.. jupyter-execute::

    from cogent3 import load_table

    table = load_table("data/stats.tsv", format="md")
    print(table)

Specify latex as the ``str()`` format
=====================================

Using the method provides finer control over formatting.

.. jupyter-execute::

    from cogent3 import load_table

    table = load_table("data/stats.tsv", format="tex")
    print(table)

Get a table as a markdown formatted string
==========================================

We use the ``justify`` argument to indicate the column justification.

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    print(table.to_markdown(justify="ccr"))

Get a table as a latex formatted string
=======================================

.. jupyter-execute::

    table = load_table(
        "data/stats.tsv", title="Some stats.", legend="Derived from something."
    )
    print(table.to_latex(justify="ccr", label="tab:table1"))

Get a table as a restructured text csv-table
============================================

.. jupyter-execute::

    table = load_table(
        "data/stats.tsv", title="Some stats.", legend="Derived from something."
    )
    print(table.to_rst(csv_table=True))

Get a table as a restructured text grid table
=============================================

.. jupyter-execute::

    table = load_table(
        "data/stats.tsv", title="Some stats.", legend="Derived from something."
    )
    print(table.to_rst())

Getting a latex format table with ``to_string()``
=================================================

It is also possible to specify column alignment, table caption and other arguments.

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    print(table.to_string(format="latex"))

Getting a bedGraph format with ``to_string()``
==============================================

This format allows display of annotation tracks on genome browsers. A small sample of a bigger table.

.. jupyter-execute::
    :hide-code:

    rows = [
        ["1", 100, 101, 1.123],
        ["1", 101, 102, 1.123],
        ["1", 102, 103, 1.123],
        ["1", 103, 104, 1.123],
        ["1", 104, 105, 1.123],
        ["1", 105, 106, 1.123],
        ["1", 106, 107, 1.123],
        ["1", 107, 108, 1.123],
        ["1", 108, 109, 1],
        ["1", 109, 110, 1],
        ["1", 110, 111, 1],
        ["1", 111, 112, 1],
        ["1", 112, 113, 1],
        ["1", 113, 114, 1],
        ["1", 114, 115, 1],
        ["1", 115, 116, 1],
        ["1", 116, 117, 1],
        ["1", 117, 118, 1],
        ["1", 118, 119, 2],
        ["1", 119, 120, 2],
        ["1", 120, 121, 2],
        ["1", 150, 151, 2],
        ["1", 151, 152, 2],
        ["1", 152, 153, 2],
        ["1", 153, 154, 2],
        ["1", 154, 155, 2],
        ["1", 155, 156, 2],
        ["1", 156, 157, 2],
        ["1", 157, 158, 2],
        ["1", 158, 159, 2],
        ["1", 159, 160, 2],
        ["1", 160, 161, 2],
    ]
    bgraph = make_table(header=["chrom", "start", "end", "value"], rows=rows)

.. jupyter-execute::

    bgraph.head()

Then converted.

.. jupyter-execute::

    print(
        bgraph.to_string(
            format="bedgraph",
            name="test track",
            description="test of bedgraph",
            color=(255, 0, 0),
            digits=0,
        )
    )

Getting a table as html
=======================

.. jupyter-execute::

    from cogent3 import load_table

    table = load_table("data/stats.tsv")
    straight_html = table.to_rich_html(compact=True)

We can provide customised formatting via a callback function.

.. jupyter-execute::

    def format_cell(value, row_num, col_num):
        style = 'style="background: rgba(176, 245, 102, 0.25);"' if value else ""
        return f"<td {style}>{value}</td>"

    rich_html = table.to_rich_html(row_cell_func=format_cell, compact=False)

Which produces the following...

.. jupyter-execute::
    :hide-code:

    from IPython.core.display import HTML
    HTML(rich_html)

We could also use control html element format.

.. jupyter-execute::

    element_format = dict(thead=f'<thead style="background: rgba(0, 250, 0, 0.1);">')
    rich_html = table.to_rich_html(element_formatters=element_format)

Which produces the following...

.. jupyter-execute::
    :hide-code:

    HTML(rich_html)

What formats can be written?
============================

Appending any of the following to a filename will cause that format to be used for writing.

.. jupyter-execute::

    from cogent3.format.table import known_formats

    known_formats

Writing a latex formmated file
==============================

.. jupyter-execute::

    table.write("stats_tab.tex", justify="ccr", label="tab:table1")

Writing delimited formats
=========================

The delimiter can be specified explicitly using the ``sep`` argument or implicitly via the file name suffix.

.. jupyter-execute::

    table = load_table("data/stats.tsv")
    table.write("stats_tab.txt", sep="\t")

..  cleanup

.. jupyter-execute::
    :hide-code:

    import pathlib
    
    for name in ("stats_tab.txt", "stats_tab.tex"):
        p = pathlib.Path(name)
        if p.exists():
            p.unlink()