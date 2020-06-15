Data Manipulation using ``Table``
=================================

.. sectionauthor:: Gavin Huttley

..
    Copyright 2007-2009, The Cogent Project
    Credits Gavin Huttley, Felix Schill
    License, GPL
    version, 1.3.0.dev
    Maintainer, Gavin Huttley
    Email, gavin.huttley@anu.edu.au
    Status, Production

The toolkit has a ``Table`` object that can be used for manipulating tabular data. It's properties can be considered like an ordered 2 dimensional dictionary or tuple with flexible output format capabilities of use for exporting data for import into external applications. Importantly, via the restructured text format one can generate html or latex formatted tables. The ``table`` module is located within ``cogent3.util``. The ``load_table`` and ``make_table`` convenience functions are provided as top-level ``cogent3`` imports.

Table creation
--------------

Tables can be created directly using the Table object itself, or a convenience function that handles loading from files. We import both here:

.. jupyter-execute::
    :linenos:

    from cogent3 import load_table, make_table
    from cogent3.util.table import Table

You can create a ``Table`` with no data.

.. jupyter-execute::
    :linenos:

    t = Table(header=["col 1", "col 2"], data=[])
    t.shape == (0, 2)

Let's create a very simple, rather nonsensical, table first. To create a table requires a header series, and a 2D series (either of type ``tuple``, ``list``, ``dict``) or a `pandas DataFrame <https://pandas.pydata.org/>`_..

.. jupyter-execute::
    :linenos:

    column_headings = ["chrom", "stableid", "length"]

The string "chrom" will become the first column heading, "stableid" the second column heading, etc. The data are,

.. jupyter-execute::
    :linenos:

    rows = [
        ["X", "ENSG00000005893", 1353],
        ["A", "ENSG00000019485", 1827],
        ["A", "ENSG00000019102", 999],
        ["X", "ENSG00000012174", 1599],
        ["X", "ENSG00000010671", 1977],
        ["A", "ENSG00000019186", 1554],
        ["A", "ENSG00000019144", 4185],
        ["X", "ENSG00000008056", 2307],
        ["A", "ENSG00000018408", 1383],
        ["A", "ENSG00000019169", 1698],
    ]

We create the simplest of tables.

.. jupyter-execute::
    :linenos:

    t = Table(header=column_headings, data=rows)
    print(t)

The format above is referred to as 'simple' format in the documentation. Notice that the numbers in this table have 4 decimal places, despite the fact the original data were largely strings and had ``max`` of 3 decimal places precision. ``Table`` converts string representations of numbers to their appropriate form when you do ``str(table)`` or print the table.

We have several things we might want to specify when creating a table: the precision and or format of floating point numbers (integer argument - ``digits``), the spacing between columns (integer argument or actual string of whitespace - ``space``), title (argument - ``title``), and legend (argument - ``legend``). Lets modify some of these and provide a title and legend.

.. jupyter-execute::
    :linenos:

    t = Table(
        header=column_headings,
        data=rows,
        title="Alignment lengths",
        legend="Some analysis",
        digits=2,
        space="        ",
    )
    print(t)

.. note:: The ``repr()`` of a table gives a quick summary.

.. jupyter-execute::
    :linenos:

    t

The Table class cannot handle arbitrary python objects, unless they are passed in as strings. Note in this case we now directly pass in the column headings list and the handling of missing data can be explicitly specified..

.. jupyter-execute::
    :linenos:

    t2 = Table(
        header=["abcd", "data"],
        data=[[str(list(range(1, 6))), "0"], ["x", 5.0], ["y", None]],
        missing_data="*",
        digits=1,
    )
    print(t2)

Table column headings can be assessed from the ``table.header`` property

.. jupyter-execute::
    :linenos:

    assert t2.header == ("abcd", "data")

this cannot be changed.

.. jupyter-execute::
    :linenos:
    :raises: TypeError

    t2.header[1] = "Data"

If you want to change the header, use the ``with_new_header`` method. This can be done one column at a time, or as a batch. The returned Table is identical aside from the modified column labels.

.. jupyter-execute::
    :linenos:

    mod_header = t2.with_new_header("abcd", "ABCD")
    assert mod_header.header == ("ABCD", "data")
    mod_header = t2.with_new_header(["abcd", "data"], ["ABCD", "DATA"])
    print(mod_header)

Tables may also be created from 2-dimensional dictionaries. In this case, special capabilities are provided to enforce printing rows in a particular order.

.. jupyter-execute::
    :linenos:

    d2D = {
        "edge.parent": {
            "NineBande": "root",
            "edge.1": "root",
            "DogFaced": "root",
            "Human": "edge.0",
            "edge.0": "edge.1",
            "Mouse": "edge.1",
            "HowlerMon": "edge.0",
        },
        "x": {
            "NineBande": 1.0,
            "edge.1": 1.0,
            "DogFaced": 1.0,
            "Human": 1.0,
            "edge.0": 1.0,
            "Mouse": 1.0,
            "HowlerMon": 1.0,
        },
        "length": {
            "NineBande": 4.0,
            "edge.1": 4.0,
            "DogFaced": 4.0,
            "Human": 4.0,
            "edge.0": 4.0,
            "Mouse": 4.0,
            "HowlerMon": 4.0,
        },
        "y": {
            "NineBande": 3.0,
            "edge.1": 3.0,
            "DogFaced": 3.0,
            "Human": 3.0,
            "edge.0": 3.0,
            "Mouse": 3.0,
            "HowlerMon": 3.0,
        },
        "z": {
            "NineBande": 6.0,
            "edge.1": 6.0,
            "DogFaced": 6.0,
            "Human": 6.0,
            "edge.0": 6.0,
            "Mouse": 6.0,
            "HowlerMon": 6.0,
        },
        "edge.name": [
            "Human",
            "HowlerMon",
            "Mouse",
            "NineBande",
            "DogFaced",
            "edge.0",
            "edge.1",
        ],
    }
    row_order = d2D["edge.name"]
    d2D["edge.name"] = dict(zip(row_order, row_order))
    t3 = Table(
        ["edge.name", "edge.parent", "length", "x", "y", "z"],
        d2D,
        row_order=row_order,
        missing_data="*",
        space=8,
        max_width=50,
        index="edge.name",
        title="My title",
        legend="legend: this is a nonsense example.",
    )
    print(t3)

In the above we specify a maximum width of the table, and also specify row identifiers (using ``index``, the name to use as row identifiers). This has the effect of forcing the table to wrap when the simple text format is used, but wrapping does not occur for any other format. The ``index`` is a column containing data for slicing the table by row, and as identifiers are presented in each wrapped sub-table.

Wrapping generates neat looking tables whether or not you index the table rows. We demonstrate here

.. jupyter-execute::
    :linenos:

    from cogent3 import make_table

    h = ["A/C", "A/G", "A/T", "C/A"]
    rows = [[0.0425, 0.1424, 0.0226, 0.0391]]
    wrap_table = make_table(header=h, data=rows, max_width=30)
    print(wrap_table)
    wrap_table = make_table(header=h, data=rows, max_width=30, index="A/C")
    print(wrap_table)

We can also customise the formatting of individual columns.

.. jupyter-execute::
    :linenos:

    rows = (
        ("NP_003077_hs_mm_rn_dna", "Con", 2.5386013224378985),
        ("NP_004893_hs_mm_rn_dna", "Con", 0.12135142635634111e06),
        ("NP_005079_hs_mm_rn_dna", "Con", 0.95165949788861326e07),
        ("NP_005500_hs_mm_rn_dna", "Con", 0.73827030202664901e-07),
        ("NP_055852_hs_mm_rn_dna", "Con", 1.0933217708952725e07),
    )

We first create a table and show the default formatting behaviour for ``Table``.

.. jupyter-execute::
    :linenos:

    t46 = Table(["Gene", "Type", "LR"], rows)
    print(t46)

We then format the ``LR`` column to use a scientific number format.

.. jupyter-execute::
    :linenos:

    t46 = Table(["Gene", "Type", "LR"], rows)
    t46.format_column("LR", "%.4e")
    print(t46)

It is safe to directly modify certain attributes, such as the title, legend and white space separating columns, which we do for the ``t46``.

.. jupyter-execute::
    :linenos:

    t46.title = "A new title"
    t46.legend = "A new legend"
    t46.space = "  "
    print(t46)

We can provide settings for multiple columns.

.. jupyter-execute::
    :linenos:

    t3 = Table(
        ["edge.name", "edge.parent", "length", "x", "y", "z"], d2D, row_order=row_order
    )
    t3.format_column("x", "%.1e")
    t3.format_column("y", "%.2f")
    print(t3)

In some cases, the contents of a column can be of different types. In this instance, rather than passing a column template we pass a reference to a function that will handle this complexity. To illustrate this we will define a function that formats floating point numbers, but returns everything else as is.

.. jupyter-execute::
    :linenos:

    def formatcol(value):
        if isinstance(value, float):
            val = "%.2f" % value
        else:
            val = str(value)
        return val

We apply this to a table with mixed string, integer and floating point data.

.. jupyter-execute::
    :linenos:

    t6 = Table(
        ["ColHead"],
        [["a"], [1], [0.3], ["cc"]],
        column_templates=dict(ColHead=formatcol),
    )
    print(t6)

Creating a Table from a pandas DataFrame
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assign the ``DataFrame`` instance to the ``data_frame`` argument.

.. jupyter-execute::
    :linenos:

    from pandas import DataFrame

    df = DataFrame(data=[[0, 1], [3, 7]], columns=["a", "b"])
    print(df)
    df_as_table = make_table(data_frame=df)
    print(df_as_table)

Representation of tables
^^^^^^^^^^^^^^^^^^^^^^^^

The representation formatting provides a quick overview of a table's dimensions and it's contents. We show this for a table with 3 columns and multiple rows

.. jupyter-execute::
    :linenos:

    t46

and larger

.. jupyter-execute::
    :linenos:

    t3

.. note:: within a script use ``print(repr(t3))`` to get the same representation.

Table output
------------

Table can output in multiple formats, including restructured text or 'rest' and delimited. These can be obtained using the ``to_string`` method and ``format`` argument as follows. Using table ``t`` from above,

.. jupyter-execute::
    :linenos:

    print(t.to_string(format="rest"))

or Markdown format

.. jupyter-execute::
    :linenos:

    print(t.to_string(format="md"))

which can also take an optional `justify` argument. The latter must be a series with a value for each column. (It only affects the html display of a Markdown table.)

.. jupyter-execute::
    :linenos:

    print(t.to_string(format="md", justify="lcr"))

where the values `lcr` correspond to left, centre and right justification.

In the case of Markdown, the pipe character (``|``) is special and so cells containing it must be escaped.

.. jupyter-execute::
    :linenos:

    md_table = make_table(
        header=["a", "b"], data=[["val1", "val2"], ["has | symbol", "val4"]]
    )
    print(md_table.to_string(format="md"))

Arguments such as ``space`` have no effect in this case. The table may also be written to file in any of the available formats (latex, simple text, html, pickle) or using a custom separator (such as a comma or tab). This makes it convenient to get data into other applications (such as R or a spreadsheet program).

The display format can be specified for a ``Table`` using any valid argument to ``to_string()``. For instance, we can make a ``Table`` instance that defaults to Markdown display.

.. jupyter-execute::
    :linenos:

    md_table = make_table(
        header=["a", "b"],
        data=[["val1", "val2"], ["has | symbol", "val4"]],
        format="md",
    )
    print(md_table)

This can be changed by modifying the `format` attribute, for example

.. jupyter-execute::
    :linenos:

    md_table.format = "rst"
    print(md_table)

Here is the latex format, note how the title and legend are joined into the latex table caption. We also provide optional arguments for the column alignment (fist column left aligned, second column right aligned and remaining columns centred) and a label for table referencing.

.. jupyter-execute::
    :linenos:

    print(t3.to_string(format="tex", justify="lrcccc", label="table:example"))

More complex latex table justifying is also possible. Specifying the width of individual columns requires passing in a series (list or tuple) of justification commands. In the following we introduce the command for specific columns widths.

.. jupyter-execute::
    :linenos:

    print(t3.to_string(format="tex", justify=["l", "p{3cm}", "c", "c", "c", "c"]))
    print(t3.to_string(sep=","))

You can specify any standard text character that will work with your desired target. Useful separators are tabs (``\t``), or pipes (``|``). If ``Table`` encounters the specified separator character within a cell, it wraps the cell in quotes -- a standard approach to facilitate import by other applications. We will illustrate this with ``t2``.

.. jupyter-execute::
    :linenos:

    print(t2.to_string(sep=","))

Test the writing of phylip distance matrix format.

.. jupyter-execute::
    :linenos:

    rows = [
        ["a", "", 0.088337278874079342, 0.18848582712597683, 0.44084000179091454],
        ["c", 0.088337278874079342, "", 0.088337278874079342, 0.44083999937417828],
        ["b", 0.18848582712597683, 0.088337278874079342, "", 0.44084000179090932],
        ["e", 0.44084000179091454, 0.44083999937417828, 0.44084000179090932, ""],
    ]
    header = ["seq1/2", "a", "c", "b", "e"]
    dist = Table(header=header, data=rows, index="seq1/2")
    print(dist.to_string(format="phylip"))

The ``to_string`` method also provides generic html generation via the restructured text format. The ``to_rich_html`` method can be used to generate the html table element by itself, with greater control over formatting. Specifically, users can provide custom callback functions to the ``row_cell_func`` and ``header_cell_func`` arguments to control in detail the formatting of table elements, or use the simpler dictionary based ``element_formatters`` approach. We use the above ``dist`` table to provide a specific callback that will set the background color for diagonal cells. We first write a function that takes the cell value and coordinates, returning the html formmatted text.

.. jupyter-execute::
    :linenos:

    def format_cell(value, row_num, col_num):
        bgcolor = ["", ' bgcolor="#0055ff"'][value == ""]
        return "<td%s>%s</td>" % (bgcolor, value)

We then call the method, without this argument, then with it.

.. jupyter-execute::
    :linenos:

    straight_html = dist.to_rich_html(compact=True)
    print(straight_html)
    rich_html = dist.to_rich_html(row_cell_func=format_cell, compact=False)
    print(rich_html)

Convert Table to pandas DataFrame
---------------------------------

If you have ``pandas`` installed, you can convert a ``Table`` instance to a ``DataFrame``.

.. jupyter-execute::
    :linenos:

    tbl = Table(header=["a", "b"], data=[[0, 1], [3, 7]])
    df = tbl.to_dataframe()
    type(df)
    print(df)

Exporting bedGraph format
-------------------------

One export format available is bedGraph_. This format can be used for viewing data as annotation track in a genome browser. This format allows for unequal spans and merges adjacent spans with the same value. The format has many possible arguments that modify the appearance in the genome browser. For this example we just create a simple data set.

.. jupyter-execute::
    :linenos:

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

    bgraph = make_table(header=["chrom", "start", "end", "value"], data=rows)

    print(
        bgraph.to_string(
            format="bedgraph",
            name="test track",
            graphType="bar",
            description="test of bedgraph",
            color=(255, 0, 0),
        )
    )

The bedgraph formatter defaults to rounding values to 2 decimal places. You can adjust that precision using the ``digits`` argument.

.. jupyter-execute::
    :linenos:

    print(
        bgraph.to_string(
            format="bedgraph",
            name="test track",
            graphType="bar",
            description="test of bedgraph",
            color=(255, 0, 0),
            digits=0,
        )
    )

.. note:: Writing files in bedgraph format is done using the ``write(format='bedgraph', name='test track', description='test of bedgraph', color=(255,0,0))``.

.. _bedGraph: https://genome.ucsc.edu/goldenPath/help/bedgraph.html

Saving a table for reloading
----------------------------

Saving a table object to file for later reloading can be done using the standard ``write`` method and ``filename`` argument to the ``Table`` constructor, specifying any of the formats supported by ``to_string``. The table loading will recreate a table from raw data located at ``filename``. To illustrate this, we first write out the table ``t3`` in ``pickle`` format, then the table ``t2`` in a csv (comma separated values format). We then remove it's header and write/reload as a tsv (tab separated values format).

.. jupyter-execute::
    :linenos:

    t3 = Table(
        ["edge.name", "edge.parent", "length", "x", "y", "z"],
        d2D,
        row_order=row_order,
        missing_data="*",
        space=8,
        max_width=50,
        index="edge.name",
        title="My title",
        legend="legend: this is a nonsense example.",
    )
    t3.write("t3.pickle")
    t3_loaded = load_table("t3.pickle")
    print(t3_loaded)
    t2 = Table(
        ["abcd", "data"],
        [[str([1, 2, 3, 4, 5]), "0"], ["x", 5.0], ["y", None]],
        missing_data="*",
        title="A \ntitle",
    )
    t2.write("t2.csv")
    t2_loaded = load_table("t2.csv", header=True, with_title=True)
    print(t2_loaded)
    t2.title = ""
    t2.write("t2.tsv")
    t2_loaded = load_table("t2.tsv")
    print(t2_loaded)

Note the ``missing_data`` attribute is not saved in the delimited format, but is in the ``pickle`` format. In the next case, I'm going to override the digits format on reloading of the table.

.. jupyter-execute::
    :linenos:

    t2 = Table(
        ["abcd", "data"],
        [[str([1, 2, 3, 4, 5]), "0"], ["x", 5.0], ["y", None]],
        missing_data="*",
        title="A \ntitle",
        legend="And\na legend too",
    )
    t2.write("t2.csv", sep=",")
    t2_loaded = load_table(
        "t2.csv", header=True, with_title=True, with_legend=True, sep=",", digits=2
    )
    print(t2_loaded)

A few things to note about the delimited file saving: formatting arguments are lost in saving to a delimited format; the ``header`` argument specifies whether the first line of the file should be treated as the header; the ``with_title`` and ``with_legend`` arguments are necessary if the file contains them, otherwise they become the header or part of the table. Importantly, if you wish to preserve numerical precision use the ``pickle`` format.

``pickle`` can load a useful object from the pickled ``Table`` by itself, without needing to know anything about the ``Table`` class.

.. jupyter-execute::
    :linenos:

    import pickle

    f = open("t3.pickle", "rb")
    pickled = pickle.load(f)
    f.close()
    sorted(pickled.keys())
    pickled["data"]["columns"]["length"]

We can read in a delimited format using a custom reader. There are two approaches. The first one allows specifying different type conversions for different columns. The second allows specifying a whole line-based parser.

You can also read and write tables in gzip compressed format. This can be done simply by ending a filename with '.gz' or specifying ``compress=True``. We write a compressed file the two different ways and read it back in.

.. jupyter-execute::
    :linenos:

    t2.write("t2.csv.gz", sep=",")
    t2_gz = load_table("t2.csv.gz", sep=",", with_title=True, with_legend=True)
    t2_gz.shape == t2.shape
    t2.write("t2.csv", sep=",", compress=True)
    t2_gz = load_table("t2.csv.gz", sep=",", with_title=True, with_legend=True)
    t2_gz.shape == t2.shape

Filtering lines on reading
--------------------------

If you only want a subset of the contents of a file, use the ``FilteringParser``. This allows skipping certain lines by using a callback function. We illustrate this using the above data, skipping any rows with ``edge.name`` starting with ``edge``.

.. jupyter-execute::
    :hide-code:

    _t = make_table(header=t3.header, data=t3.columns.to_dict())
    _t.write("t3.tab", sep="\t")


.. jupyter-execute::
    :linenos:

    from cogent3.parse.table import FilteringParser

    reader = FilteringParser(
        lambda line: not line[0].startswith("edge"), with_header=True, sep="\t"
    )
    tips = load_table("t3.tab", reader=reader, digits=1, space=2)
    print(tips)

You can also ``negate`` the condition, useful if the condition is complex (which is not really the case here).

.. jupyter-execute::
    :linenos:

    reader = FilteringParser(
        lambda line: line[0].startswith("edge"), negate=True, with_header=True, sep="\t"
    )

We can also limit the amount of data to be read in, very handy for checking large files.

.. jupyter-execute::
    :linenos:

    t3a = load_table("t3.tab", sep="\t", limit=3)
    print(t3a)

Limiting also works when ``static_column_types`` is invoked

.. jupyter-execute::
    :linenos:

    t3a = load_table("t3.tab", sep="\t", limit=3, static_column_types=True)
    t3a.shape[0] == 3

In the above example, the data type in a column is static, e.g. all values in ``x`` are floats. Rather than providing a custom reader, you can get the ``Table`` to construct such a reader based on the first data row using the ``static_column_types`` argument.

.. jupyter-execute::
    :linenos:

    t3a = load_table("t3.tab", static_column_types=True, digits=1, sep="\t")
    print(t3a)

If you invoke the ``static_column_types`` argument and the column data are not static, you'll get back a string type.

.. jupyter-execute::
    :linenos:

    t3b = make_table(header=["A", "B"], data=[[1, 1], ["a", 2]])
    print(t3b)
    t3b.write("test3b.txt", sep="\t")
    t3b = load_table("test3b.txt", sep="\t", static_column_types=True)
    t3b.columns["A"]

Table slicing and iteration
---------------------------

The Table class is capable of slicing by row, range of rows, column or range of columns headings or used to identify a single cell. Slicing using the method ``get_columns`` can also be used to reorder columns. In the case of columns, either the string headings or their position integers can be used. For rows, if ``index`` was specified, the cell values in that column can also be used.

.. jupyter-execute::
    :linenos:

    t4 = Table(
        ["edge.name", "edge.parent", "length", "x", "y", "z"],
        d2D,
        row_order=row_order,
        index="edge.name",
        title="My title",
    )

We subset ``t4`` by column and reorder them.

.. jupyter-execute::
    :linenos:

    new = t4.get_columns(["z", "y"])
    print(new)

We use the column position indexes to do get the same table.

.. jupyter-execute::
    :linenos:

    new = t4.get_columns([5, 4])
    print(new)

We can also using more general slicing, by both rows and columns. The following returns all rows from 4 on, and columns up to (but excluding) 'y':

.. jupyter-execute::
    :linenos:

    k = t4[4:, :"y"]
    print(k)

We can explicitly reference individual cells, in this case using both row and column keys.

.. jupyter-execute::
    :linenos:

    val = t4["HowlerMon", "y"]
    print(val)

We slice a single row,

.. jupyter-execute::
    :linenos:

    new = t4[3]
    print(new)

and range of rows.

.. jupyter-execute::
    :linenos:

    new = t4[3:6]
    print(new)

You can iterate over the table one row at a time and slice the rows. We illustrate this for slicing a single column,

.. jupyter-execute::
    :linenos:

    for row in t:
        print(row["stableid"])

and for multiple columns.

.. jupyter-execute::
    :linenos:

    for row in t:
        print(row["stableid"], row["length"])

The numerical slice equivalent to the first case above would be ``row[0]``, to the second case either ``row[:]``, ``row[:2]``.

Filtering tables - selecting subsets of rows/columns
----------------------------------------------------

We want to be able to slice a table, based on some condition(s), to produce a new subset table. For instance, we construct a table with type and probability values.

.. jupyter-execute::
    :linenos:

    header = ["Gene", "type", "LR", "df", "Prob"]
    rows = (
        ("NP_003077_hs_mm_rn_dna", "Con", 2.5386, 1, 0.1111),
        ("NP_004893_hs_mm_rn_dna", "Con", 0.1214, 1, 0.7276),
        ("NP_005079_hs_mm_rn_dna", "Con", 0.9517, 1, 0.3293),
        ("NP_005500_hs_mm_rn_dna", "Con", 0.7383, 1, 0.3902),
        ("NP_055852_hs_mm_rn_dna", "Con", 0.0000, 1, 0.9997),
        ("NP_057012_hs_mm_rn_dna", "Unco", 34.3081, 1, 0.0000),
        ("NP_061130_hs_mm_rn_dna", "Unco", 3.7986, 1, 0.0513),
        ("NP_065168_hs_mm_rn_dna", "Con", 89.9766, 1, 0.0000),
        ("NP_065396_hs_mm_rn_dna", "Unco", 11.8912, 1, 0.0006),
        ("NP_109590_hs_mm_rn_dna", "Con", 0.2121, 1, 0.6451),
        ("NP_116116_hs_mm_rn_dna", "Unco", 9.7474, 1, 0.0018),
    )
    t5 = Table(header, rows)
    print(t5)

We then seek to obtain only those rows that contain probabilities < 0.05. We use valid python code within a string. **Note:** Make sure your column headings could be valid python variable names or the string based approach will fail (you could use an external function instead, see below).

.. jupyter-execute::
    :linenos:

    sub_table1 = t5.filtered(callback="Prob < 0.05")
    print(sub_table1)

Using the above table we test the function to extract the raw data for a single column,

.. jupyter-execute::
    :linenos:

    raw = sub_table1.tolist("LR")
    raw

and from multiple columns.

.. jupyter-execute::
    :linenos:

    raw = sub_table1.tolist(columns=["df", "Prob"])
    raw

We can also do filtering using an external function, in this case we use a ``lambda`` to obtain only those rows of type 'Unco' that contain probabilities < 0.05, modifying our callback function.

.. jupyter-execute::
    :linenos:

    sub_table2 = t5.filtered(
        lambda ty_pr: ty_pr[0] == "Unco" and ty_pr[1] < 0.05, columns=("type", "Prob")
    )
    print(sub_table2)

This can also be done using the string approach.

.. jupyter-execute::
    :linenos:

    sub_table2 = t5.filtered("type == 'Unco' and Prob < 0.05")
    print(sub_table2)

We can also filter table columns using ``filtered_by_column``. Say we only want the numerical columns, we can write a callback that returns ``False`` if some numerical operation fails, ``True`` otherwise.

.. jupyter-execute::
    :linenos:

    def is_numeric(values):
        try:
            sum(values)
        except TypeError:
            return False
        return True

    print(t5.filtered_by_column(callback=is_numeric))

Appending tables
----------------

Tables may also be appended to each other, to make larger tables. We'll construct two simple tables to illustrate this.

.. jupyter-execute::
    :linenos:

    geneA = Table(
        ["edge.name", "edge.parent", "z"],
        [["Human", "root", 6.0], ["Mouse", "root", 6.0], ["Rat", "root", 6.0]],
        title="Gene A",
    )
    geneB = Table(
        ["edge.name", "edge.parent", "z"],
        [["Human", "root", 7.0], ["Mouse", "root", 7.0], ["Rat", "root", 7.0]],
        title="Gene B",
    )
    print(geneB)

we now use the ``appended`` Table method to create a new table, specifying that we want a new column created (by passing the ``new_column`` argument a heading) in which the table titles will be placed.

.. jupyter-execute::
    :linenos:

    new = geneA.appended("Gene", geneB, title="Appended tables")
    print(new)

We repeat this without adding a new column.

.. jupyter-execute::
    :linenos:

    new = geneA.appended(None, geneB, title="Appended, no new column")
    print(new)

Miscellaneous
-------------

Tables have a ``shape`` attribute, which specifies *x* (number of columns) and *y* (number of rows). The attribute is a tuple and we illustrate it for the above ``sub_table`` tables. Combined with the ``filtered`` method, this attribute can tell you how many rows satisfy a specific condition.

.. jupyter-execute::
    :linenos:

    t5.shape
    sub_table1.shape
    sub_table2.shape

For instance, 3 of the 11 rows in ``t`` were significant and belonged to the ``Unco`` type.

For completeness, we generate a table with no rows and assess its shape.

.. jupyter-execute::
    :linenos:

    sub_table3 = t5.filtered(
        lambda ty_pr: ty_pr[0] == "Unco" and ty_pr[1] > 0.1, columns=("type", "Prob")
    )
    sub_table3.shape

The distinct values can be obtained for a single column,

.. jupyter-execute::
    :linenos:

    distinct = new.distinct_values("edge.name")
    assert distinct == set(["Rat", "Mouse", "Human"]), distinct

or multiple columns

.. jupyter-execute::
    :linenos:

    distinct = new.distinct_values(["edge.parent", "z"])
    assert distinct == set([("root", 6.0), ("root", 7.0)]), distinct

We can compute column sums. Assuming only numerical values in a column.

.. jupyter-execute::
    :linenos:

    assert new.summed("z") == 39.0, new.summed("z")

We construct an example with mixed numerical and non-numerical data. We now compute the column sum with mixed non-numerical/numerical data.

.. jupyter-execute::
    :linenos:

    mix = make_table(header=["A", "B"], data=[[0, ""], [1, 2], [3, 4]])
    print(mix)
    mix.summed("B", strict=False)

We also compute row sums for the pure numerical and mixed non-numerical/numerical rows. For summing across rows we must specify the actual row index as an ``int``.

.. jupyter-execute::
    :linenos:

    mix.summed(0, col_sum=False, strict=False)
    mix.summed(1, col_sum=False)

We can compute the totals for all columns or rows too.

.. jupyter-execute::
    :linenos:

    mix.summed(strict=False)
    mix.summed(col_sum=False, strict=False)

We test these for a strictly numerical table.

.. jupyter-execute::
    :linenos:

    non_mix = make_table(header=["A", "B"], data=[[0, 1], [1, 2], [3, 4]])
    non_mix.summed()
    non_mix.summed(col_sum=False)

We can normalise a numerical table by row,

.. jupyter-execute::
    :linenos:

    print(non_mix.normalized(by_row=True))

or by column, such that the row/column sums are 1.

.. jupyter-execute::
    :linenos:

    print(non_mix.normalized(by_row=False))

We normalize by an arbitrary function (maximum value) by row,

.. jupyter-execute::
    :linenos:

    print(non_mix.normalized(by_row=True, denominator_func=max))

by column.

.. jupyter-execute::
    :linenos:

    print(non_mix.normalized(by_row=False, denominator_func=max))

Extending tables
----------------

In some cases it is desirable to compute an additional column from existing column values. This is done using the ``with_new_column`` method. We'll use t4 from above, adding two of the columns to create an additional column.

.. jupyter-execute::
    :linenos:

    t7 = t4.with_new_column("Sum", callback="z+x", digits=2)
    print(t7)

We test this with an externally defined function.

.. jupyter-execute::
    :linenos:

    func = lambda x_y: x_y[0] * x_y[1]
    t7 = t4.with_new_column("Sum", callback=func, columns=("y", "z"), digits=2)
    print(t7)
    func = lambda x: x ** 3
    t7 = t4.with_new_column("Sum", callback=func, columns="y", digits=2)
    print(t7)

Sorting tables
--------------

We want a table sorted according to values in a column.

.. jupyter-execute::
    :linenos:

    sorted = t5.sorted(columns="LR")
    print(sorted)

We want a table sorted according to values in a subset of columns, note the order of columns determines the sort order.

.. jupyter-execute::
    :linenos:

    sorted = t5.sorted(columns=("LR", "type"))
    print(sorted)

We now do a sort based on 2 columns.

.. jupyter-execute::
    :linenos:

    sorted = t5.sorted(columns=("type", "LR"))
    print(sorted)

Reverse sort a single column

.. jupyter-execute::
    :linenos:

    sorted = t5.sorted("LR", reverse="LR")
    print(sorted)

Sort by just specifying the ``reverse`` column

.. jupyter-execute::
    :linenos:

    sorted = t5.sorted(reverse="LR")
    print(sorted)

Reverse sort one column but not another

.. jupyter-execute::
    :linenos:

    sorted = t5.sorted(columns=("type", "LR"), reverse="LR")
    print(sorted)

Reverse sort both columns.

.. jupyter-execute::
    :linenos:

    sorted = t5.sorted(columns=("type", "LR"), reverse=("type", "LR"))
    print(sorted)

Joining Tables
--------------

The Table object is capable of joins or merging of records in two tables. There are two fundamental types of joins -- inner and outer -- with there being different sub-types. We demonstrate these first constructing some simple tables.

.. jupyter-execute::
    :linenos:

    a = Table(
        header=["index", "col2", "col3"],
        data=[[1, 2, 3], [2, 3, 1], [2, 6, 5]],
        title="A",
    )
    print(a)
    b = Table(
        header=["index", "col2", "col3"],
        data=[[1, 2, 3], [2, 2, 1], [3, 6, 3]],
        title="B",
    )
    print(b)
    c = Table(header=["index", "col_c2"], rows=[[1, 2], [3, 2], [3, 5]], title="C")
    print(c)

For a natural inner join, only 1 copy of columns with the same name are retained. So we expect the headings to be identical between the table ``a``/``b`` and the result of ``a.joined(b)`` or ``b.joined(a)``.

.. jupyter-execute::
    :linenos:

    assert a.joined(b).header == b.header
    assert b.joined(a).header == a.header

For a standard inner join, the joined table should contain all columns from ``a`` and ``b`` excepting the index column(s). Simply providing a column name (or index) selects this behaviour. Note that in this case, column names from the second table are made unique by prefixing them with that tables title. If the right table does not have a title, a default value `right` is used.

.. jupyter-execute::
    :linenos:

    b.title = None
    c.joined(b)
    b.title = "B"
    assert a.joined(b, "index").header == ("index", "col2", "col3", "B_col2", "B_col3")

Note that the table title's were used to prefix the column headings from the second table. We further test this using table ``c`` which has different dimensions.

.. jupyter-execute::
    :linenos:

    assert a.joined(c, "index").header == ("index", "col2", "col3", "C_col_c2")

It's also possible to specify index columns using numerical values, the results of which should be the same.

.. jupyter-execute::
    :linenos:

    r1 = a.joined(b, [0, 2])
    r2 = a.joined(b, ["index", "col3"])
    assert r1.tolist() == r2.tolist()

Additionally, it's possible to provide two series of indices for the two tables. Here, they have identical values.

.. jupyter-execute::
    :linenos:

    assert (
        a.joined(b, ["index", "col3"], ["index", "col3"]).tolist()
        == a.joined(b, ["index", "col3"]).tolist()
    )

The results of a standard join between tables ``a`` and ``b`` are

.. jupyter-execute::
    :linenos:

    print(a.joined(b, ["index"], title="A&B"))

We demo the table specific indices.

.. jupyter-execute::
    :linenos:

    print(a.joined(c, ["col2"], ["index"], title='A&C by "col2/index"'))

Tables ``a`` and ``c`` share a single row with the same value in the ``index`` column, hence a join by that index should return a table with just that row.

.. jupyter-execute::
    :linenos:

    print(a.joined(c, "index", title='A&C by "index"'))

A natural join of tables ``a`` and ``b`` results in a table with only rows that were identical between the two parents.

.. jupyter-execute::
    :linenos:

    print(a.joined(b, title="A&B Natural Join"))

We test the outer join by defining an additional table with different dimensions, and conducting a join specifying ``inner_join=False``.

.. jupyter-execute::
    :linenos:

    d = Table(header=["index", "col_c2"], data=[[5, 42], [6, 23]], title="D")
    print(d)
    print(c.joined(d, inner_join=False, title="C&D Outer join"))

We establish the ``joined`` method works for mixtures of character and numerical data, setting some indices and some cell values to be strings.

.. jupyter-execute::
    :linenos:

    a = Table(
        header=["index", "col2", "col3"],
        data=[[1, 2, "3"], ["2", 3, 1], [2, 6, 5]],
        title="A",
    )
    b = Table(
        header=["index", "col2", "col3"],
        data=[[1, 2, "3"], ["2", 2, 1], [3, 6, 3]],
        title="B",
    )
    assert (
        a.joined(b, ["index", "col3"], ["index", "col3"]).tolist()
        == a.joined(b, ["index", "col3"]).tolist()
    )

We test that the ``joined`` method works when the column index orders differ.

.. jupyter-execute::
    :linenos:

    t1_header = ["a", "b"]
    t1_rows = [(1, 2), (3, 4)]
    t2_header = ["b", "c"]
    t2_rows = [(3, 6), (4, 8)]
    t1 = Table(t1_header, data=t1_rows, title="t1")
    t2 = Table(t2_header, data=t2_rows, title="t2")
    t3 = t1.joined(t2, columns_self=["b"], columns_other=["b"])
    print(t3)

We then establish that a join with no values does not cause a failure, just returns an empty ``Table``.

.. jupyter-execute::
    :linenos:

    t4_header = ["b", "c"]
    t4_rows = [(5, 6), (7, 8)]
    t4 = make_table(header=t4_header, data=t4_rows)
    t4.title = "t4"
    t5 = t1.joined(t4, columns_self=["b"], columns_other=["b"])
    print(t5)

Whose representation looks like

.. jupyter-execute::
    :linenos:

    t5

Transposing a table
-------------------

Tables can be transposed.

.. jupyter-execute::
    :linenos:

    from cogent3 import make_table

    title = "#Full OTU Counts"
    header = ["#OTU ID", "14SK041", "14SK802"]
    rows = [
        [-2920, "332", 294],
        [-1606, "302", 229],
        [-393, 141, 125],
        [-2109, 138, 120],
        [-5439, 104, 117],
        [-1834, 70, 75],
        [-18588, 65, 47],
        [-1350, 60, 113],
        [-2160, 57, 52],
        [-11632, 47, 36],
    ]
    table = make_table(header=header, rows=rows, title=title)
    print(table)

We now transpose this. We require a new column heading for header data and an identifier for which existing column will become the header (default is index 0).

.. jupyter-execute::
    :linenos:

    tp = table.transposed(new_column_name="sample", select_as_header="#OTU ID", space=2)
    print(tp)

We test transposition with default value is the same.

.. jupyter-execute::
    :linenos:

    tp = table.transposed(new_column_name="sample", space=2)
    print(tp)

We test transposition selecting a different column to become the header.

.. jupyter-execute::
    :linenos:

    tp = table.transposed(new_column_name="sample", select_as_header="14SK802", space=2)
    print(tp)

Counting rows
-------------

We can count the number of rows for which a condition holds. This method uses the same arguments as ``filtered`` but returns an integer result only.

.. jupyter-execute::
    :linenos:

    print(c.count("col_c2 == 2"))
    print(c.joined(d, inner_join=False).count("index==3 and D_index==5"))

Testing a sub-component
-----------------------

Before using ``Table``, we exercise some formatting code:

.. jupyter-execute::
    :linenos:

    from cogent3.format.table import formatted_cells, phylip_matrix, latex

We check we can format an arbitrary 2D list, without a header, using the ``formatted_cells`` function directly.

.. jupyter-execute::
    :linenos:

    data = [[230, "acdef", 1.3], [6, "cc", 1.9876]]
    head = ["one", "two", "three"]
    header, formatted = formatted_cells(data, header=head)
    print(formatted)
    print(header)

We directly test the latex formatting.

.. jupyter-execute::
    :linenos:

    print(
        latex(formatted, header, justify="lrl", caption="A legend", label="table:test")
    )

..
    Import the ``os`` module so some file cleanup can be done at the end. To check the contents of those files, just delete the following prior to running the test. The try/except clause below is aimed at case where ``junk.pdf`` wasn't created due to ``reportlab`` not being present.

.. jupyter-execute::
    :hide-code:

    import os

    to_delete = ["t3.pickle", "t2.csv", "t2.tsv", "t2.csv.gz", "t3.tab", "test3b.txt"]
    for f in to_delete:
        try:
            os.remove(f)
        except OSError:
            pass
