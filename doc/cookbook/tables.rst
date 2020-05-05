Tabular data
------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul

Table handles tabular data, storing as columns in a, you guessed it, ``columns`` attribute. The latter acts like a dictionary, with the column names as the keys and the column values being  ``numpy.ndarray`` instances. The table itself is iterable over rows.

.. note:: ``Table`` is immutable at the level of the individual ``ndarray`` not being writable.

.. note:: ``Table`` can be interconverted from a ``pandas.DataFrame``.

.. include:: ./loading_tabular.rst

Iterating over table rows
^^^^^^^^^^^^^^^^^^^^^^^^^

``Table`` is a row oriented object. Iterating on the table returns each row as a new ``Table`` instance.

.. doctest::

    >>> from cogent3 import load_table
    >>> table = load_table("data/stats.tsv")
    >>> for row in table:
    ...     print(row)
    ...     break
    ...
    =============================
        Locus    Region     Ratio
    -----------------------------
    NP_003077       Con    2.5386
    -----------------------------

The resulting rows can be indexed using their column names.

.. doctest::

    >>> for row in table:
    ...     print(row["Locus"])
    NP_003077
    NP_004893
    NP_005079
    NP_005500
    NP_055852

Iterating over table columns
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``Table.columns`` attribute is a ``Columns`` instance, an object with ``dict`` attributes.

.. doctest::

    >>> from cogent3 import load_table
    >>> table = load_table("data/stats.tsv")
    >>> table.columns
    Columns('Locus': <U9, 'Region': <U6, 'Ratio': float64)
    >>> print(table.columns["Region"])
    ['Con' 'Con' 'Con' 'NonCon' 'NonCon']

So iteration is the same as for dicts.

.. doctest::

    >>> for name in table.columns:
    ...     print(name)
    ...
    Locus
    Region
    Ratio

Table slicing using column names
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> print(table[:2, :'Region'])
    =========
        Locus
    ---------
    NP_003077
    NP_004893
    ---------

Table slicing using column indices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> print(table[:2,: 1])
    =========
        Locus
    ---------
    NP_003077
    NP_004893
    ---------

Changing displayed numerical precision
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We change the ``Ratio`` column to using scientific notation.

.. doctest::

    >>> from cogent3 import load_table
    >>> table = load_table("data/stats.tsv")
    >>> table.format_column('Ratio', '%.1e')
    >>> print(table)
    ==============================
        Locus    Region      Ratio
    ------------------------------
    NP_003077       Con    2.5e+00
    NP_004893       Con    1.2e+05
    NP_005079       Con    9.5e+06
    NP_005500    NonCon    7.4e-08
    NP_055852    NonCon    1.1e+07
    ------------------------------

Change digits or column spacing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This can be done on table loading,

.. doctest::

    >>> table = load_table('data/stats.tsv', digits=1, space=2)
    >>> print(table)
    =============================
        Locus  Region       Ratio
    -----------------------------
    NP_003077     Con         2.5
    NP_004893     Con    121351.4
    NP_005079     Con   9516595.0
    NP_005500  NonCon         0.0
    NP_055852  NonCon  10933217.7
    -----------------------------

or, for spacing at least, by modifying the attributes

.. doctest::

    >>> table.space = '    '
    >>> print(table)
    =================================
        Locus    Region         Ratio
    ---------------------------------
    NP_003077       Con           2.5
    NP_004893       Con      121351.4
    NP_005079       Con     9516595.0
    NP_005500    NonCon           0.0
    NP_055852    NonCon    10933217.7
    ---------------------------------

Changing column headings
^^^^^^^^^^^^^^^^^^^^^^^^

The table ``header`` is immutable. Changing column headings is done as follows.

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> print(table.header)
    ('Locus', 'Region', 'Ratio')
    >>> table = table.with_new_header('Ratio', 'Stat')
    >>> print(table.header)
    ('Locus', 'Region', 'Stat')

Adding a new column
^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import make_table
    >>> table = make_table()
    >>> table
    0 rows x 0 columns
    >>> table.columns["a"] = [1, 3, 5]
    >>> table.columns["b"] = [2, 4, 6]
    >>> table
    ======
    a    b
    ------
    1    2
    3    4
    5    6
    ------
    <BLANKLINE>
    3 rows x 2 columns

Create a new column from existing ones
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This can be used to take a single, or multiple columns and generate a new column of values. Here we'll take 2 columns and return True/False based on a condition.

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> table = table.with_new_column('LargeCon',
    ...                     lambda r_v: r_v[0] == 'Con' and r_v[1]>10.0,
    ...                     columns=['Region', 'Ratio'])
    >>> print(table)
    ================================================
        Locus    Region            Ratio    LargeCon
    ------------------------------------------------
    NP_003077       Con           2.5386       False
    NP_004893       Con      121351.4264        True
    NP_005079       Con     9516594.9789        True
    NP_005500    NonCon           0.0000       False
    NP_055852    NonCon    10933217.7090       False
    ------------------------------------------------

Get table data as a numpy array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> table.array
    array([['NP_003077', 'Con', 2.5386013224378985],
           ['NP_004893', 'Con', 121351.42635634111],
           ['NP_005079', 'Con', 9516594.978886133],
           ['NP_005500', 'NonCon', 7.382703020266491e-08],
           ['NP_055852', 'NonCon', 10933217.708952725]], dtype=object)

Get a table column as a list
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Via the ``Table.tolist()`` method.

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> locus = table.tolist("Locus")
    >>> locus
    ['NP_003077', 'NP_004893', 'NP_005079', 'NP_005500', 'NP_055852']

Or directly from the column array object.

.. doctest::

    >>> table.columns["Locus"].tolist()
    ['NP_003077', 'NP_004893', 'NP_005079', 'NP_005500', 'NP_055852']

Get multiple table columns as a list
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This returns a row oriented list.

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> rows = table.tolist(["Region", "Locus"])
    >>> rows
    [['Con', 'NP_003077'], ['Con', 'NP_004893'], ...

.. note:: column name order dictates the element order per row

Get the table as a row oriented ``dict``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Keys in the resulting dict are the row indices, the value is a dict of column name, value pairs.

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> table.to_dict()
    {0: {'Locus': 'NP_003077', 'Region': 'Con', 'Ratio': ...


Get the table as a column oriented ``dict``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Keys in the resulting dict are the column names, the value is a list.

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> table.columns.to_dict()
    {'Locus': ['NP_003077', 'NP_004893',...

Get the table as a ``pandas.DataFrame``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> df = table.to_dataframe()
    >>> df
           Locus  Region         Ratio
    0  NP_003077     Con  2.538601e+00
    1  NP_004893     Con  1.213514e+05
    2  NP_005079     Con  9.516595e+06
    3  NP_005500  NonCon  7.382703e-08
    4  NP_055852  NonCon  1.093322e+07

You can also specify column(s) are categories

.. doctest::
    
    >>> df = table.to_dataframe(categories="Region")

Appending tables
^^^^^^^^^^^^^^^^

Can be done without specifying a new column. Here we simply use the same table data.

.. doctest::

    >>> table1 = load_table("data/stats.tsv")
    >>> table2 = load_table("data/stats.tsv")
    >>> table = table1.appended(None, table2)
    >>> print(table)
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_003077       Con           2.5386
    NP_004893       Con      121351.4264
    NP_005079       Con     9516594.9789
    NP_005500    NonCon           0.0000
    NP_055852    NonCon    10933217.7090
    NP_003077       Con           2.5386
    NP_004893       Con      121351.4264
    NP_005079       Con     9516594.9789
    NP_005500    NonCon           0.0000
    NP_055852    NonCon    10933217.7090
    ------------------------------------

or with a new column

.. doctest::

    >>> table1.title = 'Data1'
    >>> table2.title = 'Data2'
    >>> table = table1.appended('Data#', table2, title='')
    >>> print(table)
    =============================================
    Data#        Locus    Region            Ratio
    ---------------------------------------------
    Data1    NP_003077       Con           2.5386
    Data1    NP_004893       Con      121351.4264
    Data1    NP_005079       Con     9516594.9789
    Data1    NP_005500    NonCon           0.0000
    Data1    NP_055852    NonCon    10933217.7090
    Data2    NP_003077       Con           2.5386
    Data2    NP_004893       Con      121351.4264
    Data2    NP_005079       Con     9516594.9789
    Data2    NP_005500    NonCon           0.0000
    Data2    NP_055852    NonCon    10933217.7090
    ---------------------------------------------

.. note:: We assigned an empty string to ``title``, otherwise the resulting table has the same ``title`` attribute as that of ``table1``.

Summing a single column
^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> table.summed('Ratio')
    20571166.652...

Because each column is just a ``numpy.ndarray``, this also can be done directly via the array methods.

.. doctest::

    >>> table.columns["Ratio"].sum()
    20571166.652...

Summing multiple columns or rows - strictly numerical data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We define a strictly numerical table,

.. doctest::

    >>> from cogent3 import make_table
    >>> all_numeric = make_table(header=['A', 'B', 'C'], data=[range(3),
    ...                                 range(3,6), range(6,9), range(9,12)])
    >>> print(all_numeric)
    =============
    A     B     C
    -------------
    0     1     2
    3     4     5
    6     7     8
    9    10    11
    -------------

and sum all columns (default condition)

.. doctest::

    >>> all_numeric.summed()
    [18, 22, 26]

and all rows

.. doctest::

    >>> all_numeric.summed(col_sum=False)
    [3, 12, 21, 30]

Summing multiple columns or rows with mixed non-numeric/numeric data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We define a table with mixed data, like a distance matrix.

.. doctest::

    >>> mixed = make_table(header=['A', 'B', 'C'], data=[['*',1,2], [3,'*', 5],
    ...                                                 [6,7,'*']])
    >>> print(mixed)
    ===========
    A    B    C
    -----------
    *    1    2
    3    *    5
    6    7    *
    -----------

and sum all columns (default condition), ignoring non-numerical data

.. doctest::

    >>> mixed.summed(strict=False)
    [9, 8, 7]

and all rows

.. doctest::

    >>> mixed.summed(col_sum=False, strict=False)
    [3, 8, 13]


Filtering table rows
^^^^^^^^^^^^^^^^^^^^

We can do this by providing a reference to an external function

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> sub_table = table.filtered(lambda x: x < 10.0, columns='Ratio')
    >>> print(sub_table)
    =============================
        Locus    Region     Ratio
    -----------------------------
    NP_003077       Con    2.5386
    NP_005500    NonCon    0.0000
    -----------------------------

or using valid python syntax within a string, which is executed

.. doctest::

    >>> sub_table = table.filtered("Ratio < 10.0")
    >>> print(sub_table)
    =============================
        Locus    Region     Ratio
    -----------------------------
    NP_003077       Con    2.5386
    NP_005500    NonCon    0.0000
    -----------------------------

You can also filter for values in multiple columns

.. doctest::

    >>> sub_table = table.filtered("Ratio < 10.0 and Region == 'NonCon'")
    >>> print(sub_table)
    =============================
        Locus    Region     Ratio
    -----------------------------
    NP_005500    NonCon    0.0000
    -----------------------------

Filtering table columns
^^^^^^^^^^^^^^^^^^^^^^^

We select only columns that have a sum > 20 from the ``all_numeric`` table constructed above.

.. doctest::

    >>> big_numeric = all_numeric.filtered_by_column(lambda x: sum(x)>20)
    >>> print(big_numeric)
    ========
     B     C
    --------
     1     2
     4     5
     7     8
    10    11
    --------

Standard sorting
^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> print(table.sorted(columns='Ratio'))
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_005500    NonCon           0.0000
    NP_003077       Con           2.5386
    NP_004893       Con      121351.4264
    NP_005079       Con     9516594.9789
    NP_055852    NonCon    10933217.7090
    ------------------------------------

Reverse sorting
^^^^^^^^^^^^^^^

.. doctest::

    >>> print(table.sorted(columns='Ratio', reverse='Ratio'))
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_055852    NonCon    10933217.7090
    NP_005079       Con     9516594.9789
    NP_004893       Con      121351.4264
    NP_003077       Con           2.5386
    NP_005500    NonCon           0.0000
    ------------------------------------

Sorting involving multiple columns, one reversed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> print(table.sorted(columns=['Region', 'Ratio'], reverse='Ratio'))
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_005079       Con     9516594.9789
    NP_004893       Con      121351.4264
    NP_003077       Con           2.5386
    NP_055852    NonCon    10933217.7090
    NP_005500    NonCon           0.0000
    ------------------------------------

Getting raw data for a single column
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> raw = table.tolist('Region')
    >>> print(raw)
    ['Con', 'Con', 'Con', 'NonCon', 'NonCon']

Getting raw data for multiple columns
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> raw = table.tolist(['Locus', 'Region'])
    >>> print(raw)
    [['NP_003077', 'Con'], ['NP_004893', 'Con'], ...

Getting distinct values
^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> assert table.distinct_values('Region') == set(['NonCon', 'Con'])

Counting occurrences of values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> assert table.count("Region == 'NonCon' and Ratio > 1") == 1

Joining or merging tables
^^^^^^^^^^^^^^^^^^^^^^^^^

We do a standard inner join here for a restricted subset. We must specify the columns that will be used for the join. Here we just use ``Locus``.

.. doctest::

    >>> rows = [['NP_004893', True], ['NP_005079', True],
    ...         ['NP_005500', False], ['NP_055852', False]]
    >>> region_type = make_table(header=['Locus', 'LargeCon'], data=rows)
    >>> stats_table = load_table("data/stats.tsv")
    >>> new = stats_table.joined(region_type, columns_self='Locus')
    >>> print(new)
    ======================================================
        Locus    Region            Ratio    right_LargeCon
    ------------------------------------------------------
    NP_004893       Con      121351.4264              True
    NP_005079       Con     9516594.9789              True
    NP_005500    NonCon           0.0000             False
    NP_055852    NonCon    10933217.7090             False
    ------------------------------------------------------

.. note:: If the tables have titles, column names are prefixed with those instead of ``right_``.

.. note:: The ``joined()`` method is just a wrapper for the ``inner_join()`` and ``cross_join()`` (row cartesian product) methods, which you can use directly.

Specify markdown as the ``str()`` format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using the method provides finer control over formatting.

.. doctest::

    >>> from cogent3 import load_table
    >>> table = load_table("data/stats.tsv", format="md")
    >>> print(table)
    |     Locus | Region |         Ratio |
    |-----------|--------|---------------|
    | NP_003077 |    Con |        2.5386 |
    | NP_004893 |    Con |   121351.4264 |
    | NP_005079 |    Con |  9516594.9789 |
    | NP_005500 | NonCon |        0.0000 |
    | NP_055852 | NonCon | 10933217.7090 |

Specify latex as the ``str()`` format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using the method provides finer control over formatting.

.. doctest::

    >>> from cogent3 import load_table
    >>> table = load_table("data/stats.tsv", format="tex")
    >>> print(table)
    \begin{table}[htp!]
    \centering
    \begin{tabular}{ r r r }
    \hline
    \bf{Locus} & \bf{Region} & \bf{Ratio} \\
    \hline
    \hline
    NP_003077 &    Con &        2.5386 \\...

Get a table as a markdown formatted string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use the ``justify`` argument to indicate the column justification.

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> print(table.to_markdown(justify="ccr"))
    |     Locus | Region |         Ratio |
    |:---------:|:------:|--------------:|
    | NP_003077 |    Con |        2.5386 |
    | NP_004893 |    Con |   121351.4264 |
    | NP_005079 |    Con |  9516594.9789 |
    | NP_005500 | NonCon |        0.0000 |
    | NP_055852 | NonCon | 10933217.7090 |


Get a table as a latex formatted string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv", title="Some stats.",
    ...                    legend="Derived from something.")
    >>> print(table.to_latex(justify="ccr", label="tab:table1"))
    \begin{table}[htp!]
    \centering
    \begin{tabular}{ c c r }...

Get a table as a restructured text csv-table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv", title="Some stats.",
    ...                    legend="Derived from something.")
    >>> print(table.to_rst(csv_table=True))
    .. csv-table:: Some stats.
        :header: "Locus", "Region", "Ratio"
    <BLANKLINE>
        NP_003077, Con, 2.5386...

Get a table as a restructured text grid table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv", title="Some stats.",
    ...                    legend="Derived from something.")
    >>> print(table.to_rst())
    +------------------------------------+
    |            Some stats.             |
    +-----------+--------+---------------+
    |     Locus | Region |         Ratio |
    +===========+========+===============+
    | NP_003077 |    Con |        2.5386 |
    +-----------+--------+---------------...


What formats can be written?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3.format.table import known_formats
    >>> known_formats
    ('bedgraph', 'phylip', 'rest', 'rst', 'markdown', 'md', 'latex', 'tex', 'html', 'simple', 'csv', 'tsv')

Writing delimited formats
^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> table.write('stats_tab.txt', sep='\t')

Writing latex format
^^^^^^^^^^^^^^^^^^^^

It is also possible to specify column alignment, table caption and other arguments.

.. doctest::

    >>> table = load_table("data/stats.tsv")
    >>> print(table.to_string(format='latex'))
    \begin{table}[htp!]
    \centering
    \begin{tabular}{ r r r }
    \hline
    \bf{Locus} & \bf{Region} & \bf{Ratio} \\
    \hline
    \hline
    NP_003077 &    Con &        2.5386 \\
    NP_004893 &    Con &   121351.4264 \\
    NP_005079 &    Con &  9516594.9789 \\
    NP_005500 & NonCon &        0.0000 \\
    NP_055852 & NonCon & 10933217.7090 \\
    \hline
    \end{tabular}
    \end{table}

Writing bedGraph format
^^^^^^^^^^^^^^^^^^^^^^^

This format allows display of annotation tracks on genome browsers.

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> rows = [['1', 100, 101, 1.123], ['1', 101, 102, 1.123],
    ...         ['1', 102, 103, 1.123], ['1', 103, 104, 1.123],
    ...         ['1', 104, 105, 1.123], ['1', 105, 106, 1.123],
    ...         ['1', 106, 107, 1.123], ['1', 107, 108, 1.123],
    ...         ['1', 108, 109, 1], ['1', 109, 110, 1],
    ...         ['1', 110, 111, 1], ['1', 111, 112, 1],
    ...         ['1', 112, 113, 1], ['1', 113, 114, 1],
    ...         ['1', 114, 115, 1], ['1', 115, 116, 1],
    ...         ['1', 116, 117, 1], ['1', 117, 118, 1],
    ...         ['1', 118, 119, 2], ['1', 119, 120, 2],
    ...         ['1', 120, 121, 2], ['1', 150, 151, 2],
    ...         ['1', 151, 152, 2], ['1', 152, 153, 2],
    ...         ['1', 153, 154, 2], ['1', 154, 155, 2],
    ...         ['1', 155, 156, 2], ['1', 156, 157, 2],
    ...         ['1', 157, 158, 2], ['1', 158, 159, 2],
    ...         ['1', 159, 160, 2], ['1', 160, 161, 2]]
    ...
    >>> bgraph = make_table(header=['chrom', 'start', 'end', 'value'],
    ...                   rows=rows)
    ...
    >>> print(bgraph.to_string(format='bedgraph', name='test track',
    ...     description='test of bedgraph', color=(255,0,0))) # doctest: +NORMALIZE_WHITESPACE
    track type=bedGraph name="test track" description="test of bedgraph" color=255,0,0
    1	100	108	1.12
    1	108	118	1.00
    1	118	161	2.00

    The bedgraph formatter defaults to rounding values to 2 decimal places. You can adjust that precision using the ``digits`` argument.

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print(bgraph.to_string(format='bedgraph', name='test track',
    ...   description='test of bedgraph', color=(255,0,0), digits=0)) # doctest: +NORMALIZE_WHITESPACE
    track type=bedGraph name="test track" description="test of bedgraph" color=255,0,0
    1	100	118	1.00
    1	118	161	2.00
