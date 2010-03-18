Tabular data
------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul

.. doctest::
    :hide:

    >>> # just saving some tabular data for subsequent data
    >>> from cogent import LoadTable
    >>> rows = (('NP_003077', 'Con', 2.5386013224378985),
    ... ('NP_004893', 'Con', 0.12135142635634111e+06),
    ... ('NP_005079', 'Con', 0.95165949788861326e+07),
    ... ('NP_005500', 'NonCon', 0.73827030202664901e-07),
    ... ('NP_055852', 'NonCon', 1.0933217708952725e+07))
    >>> table = LoadTable(header=['Locus', 'Region', 'Ratio'], rows=rows)
    >>> table.writeToFile('stats.txt', sep=',')

Loading delimited formats
^^^^^^^^^^^^^^^^^^^^^^^^^

We load a comma separated data file using the generic ``LoadTable`` function.

.. doctest::

    >>> from cogent import LoadTable
    >>> table = LoadTable('stats.txt', sep=',')
    >>> print table
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_003077       Con           2.5386
    NP_004893       Con      121351.4264
    NP_005079       Con     9516594.9789
    NP_005500    NonCon           0.0000
    NP_055852    NonCon    10933217.7090
    ------------------------------------

Reading large files
^^^^^^^^^^^^^^^^^^^

For really large files the automated conversion used by the standard read mechanism can be quite slow. If the data within a column is consistently of one type, set the ``LoadTable`` argument ``static_column_types=True``. This causes the ``Table`` object to create a custom reader.

.. doctest::

    >>> table = LoadTable('stats.txt', static_column_types=True)
    >>> print table
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_003077       Con           2.5386
    NP_004893       Con      121351.4264
    NP_005079       Con     9516594.9789
    NP_005500    NonCon           0.0000
    NP_055852    NonCon    10933217.7090
    ------------------------------------

Formatting
^^^^^^^^^^

Changing displayed numerical precision
""""""""""""""""""""""""""""""""""""""

We change the ``Ratio`` column to using scientific notation.

.. doctest::

    >>> table.setColumnFormat('Ratio', '%.1e')
    >>> print table
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
"""""""""""""""""""""""""""""""

This can be done on table loading,

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',', digits=1, space=2)
    >>> print table
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

    >>> table.Space = '    '
    >>> print table
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

The table ``Header`` is immutable. Changing column headings is done as follows.

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> print table.Header
    ['Locus', 'Region', 'Ratio']
    >>> table = table.withNewHeader('Ratio', 'Stat')
    >>> print table.Header
    ['Locus', 'Region', 'Stat']

Creating new columns from existing ones
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This can be used to take a single, or multiple columns and generate a new column of values. Here we'll take 2 columns and return True/False based on a condition.

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> table = table.withNewColumn('LargeCon',
    ...                     lambda (r,v): r == 'Con' and v>10.0,
    ...                     columns=['Region', 'Ratio'])
    >>> print table
    ================================================
        Locus    Region            Ratio    LargeCon
    ------------------------------------------------
    NP_003077       Con           2.5386       False
    NP_004893       Con      121351.4264        True
    NP_005079       Con     9516594.9789        True
    NP_005500    NonCon           0.0000       False
    NP_055852    NonCon    10933217.7090       False
    ------------------------------------------------

Appending tables
^^^^^^^^^^^^^^^^

Can be done without specifying a new column. Here we simply use the same table data.

.. doctest::

    >>> table1 = LoadTable('stats.txt', sep=',')
    >>> table2 = LoadTable('stats.txt', sep=',')
    >>> table = table1.appended(None, table2)
    >>> print table
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

    >>> table1.Title = 'Data1'
    >>> table2.Title = 'Data2'
    >>> table = table1.appended('Data#', table2, title='')
    >>> print table
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

.. note:: We assigned an empty string to ``title``, otherwise the resulting table has the same ``Title`` attribute as that of ``table1``.

Summing a single column
^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> table.summed('Ratio')
    20571166.652847398

Summing multiple columns or rows - strictly numerical data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We define a strictly numerical table,

.. doctest::

    >>> all_numeric = LoadTable(header=['A', 'B', 'C'], rows=[range(3),
    ...                                 range(3,6), range(6,9), range(9,12)])
    >>> print all_numeric
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

    >>> mixed = LoadTable(header=['A', 'B', 'C'], rows=[['*',1,2], [3,'*', 5],
    ...                                                 [6,7,'*']])
    >>> print mixed
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

    >>> table = LoadTable('stats.txt', sep=',')
    >>> sub_table = table.filtered(lambda x: x < 10.0, columns='Ratio')
    >>> print sub_table
    =============================
        Locus    Region     Ratio
    -----------------------------
    NP_003077       Con    2.5386
    NP_005500    NonCon    0.0000
    -----------------------------

or using valid python syntax within a string, which is executed

.. doctest::

    >>> sub_table = table.filtered("Ratio < 10.0")
    >>> print sub_table
    =============================
        Locus    Region     Ratio
    -----------------------------
    NP_003077       Con    2.5386
    NP_005500    NonCon    0.0000
    -----------------------------

You can also filter for values in multiple columns

.. doctest::

    >>> sub_table = table.filtered("Ratio < 10.0 and Region == 'NonCon'")
    >>> print sub_table
    =============================
        Locus    Region     Ratio
    -----------------------------
    NP_005500    NonCon    0.0000
    -----------------------------

Filtering table columns
^^^^^^^^^^^^^^^^^^^^^^^

We select only columns that have a sum > 20 from the ``all_numeric`` table constructed above.

.. doctest::

    >>> big_numeric = all_numeric.filteredByColumn(lambda x: sum(x)>20)
    >>> print big_numeric
    ========
     B     C
    --------
     1     2
     4     5
     7     8
    10    11
    --------

Sorting
^^^^^^^

Standard sorting
""""""""""""""""

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> print table.sorted(columns='Ratio')
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
"""""""""""""""

.. doctest::

    >>> print table.sorted(columns='Ratio', reverse='Ratio')
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
""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> print table.sorted(columns=['Region', 'Ratio'], reverse='Ratio')
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_005079       Con     9516594.9789
    NP_004893       Con      121351.4264
    NP_003077       Con           2.5386
    NP_055852    NonCon    10933217.7090
    NP_005500    NonCon           0.0000
    ------------------------------------

Getting raw data
^^^^^^^^^^^^^^^^

For a single column
"""""""""""""""""""

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> raw = table.getRawData('Region')
    >>> print raw
    ['Con', 'Con', 'Con', 'NonCon', 'NonCon']

For multiple columns
""""""""""""""""""""

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> raw = table.getRawData(['Locus', 'Region'])
    >>> print raw
    [['NP_003077', 'Con'], ['NP_004893', 'Con'], ...

Iterating over table rows
^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> for row in table:
    ...     print row['Locus']
    ...
    NP_003077
    NP_004893
    NP_005079
    NP_005500
    NP_055852

Table slicing
^^^^^^^^^^^^^

Using column names
""""""""""""""""""

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> print table[:2, :'Region']
    =========
        Locus
    ---------
    NP_003077
    NP_004893
    ---------

Using column indices
""""""""""""""""""""

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> print table[:2,: 1]
    =========
        Locus
    ---------
    NP_003077
    NP_004893
    ---------

SQL-like capabilities
^^^^^^^^^^^^^^^^^^^^^

Distinct values
"""""""""""""""

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> assert table.getDistinctValues('Region') == set(['NonCon', 'Con'])

Counting
""""""""

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> assert table.count("Region == 'NonCon' and Ratio > 1") == 1

Joining tables
""""""""""""""

SQL-like join operations requires tables have different ``Title`` attributes which are not ``None``. We do a standard inner join here for a restricted subset. We must specify the columns that will be used for the join. Here we just use ``Locus`` but multiple columns can be used, and their names can be different between the tables. Note that the second table's title becomes a part of the column names.

.. doctest::

    >>> rows = [['NP_004893', True], ['NP_005079', True],
    ...         ['NP_005500', False], ['NP_055852', False]]
    >>> region_type = LoadTable(header=['Locus', 'LargeCon'], rows=rows,
    ...                 title='RegionClass')
    >>> stats_table = LoadTable('stats.txt', sep=',', title='Stats')
    >>> new = stats_table.joined(region_type, columns_self='Locus')
    >>> print new
    ============================================================
        Locus    Region            Ratio    RegionClass_LargeCon
    ------------------------------------------------------------
    NP_004893       Con      121351.4264                    True
    NP_005079       Con     9516594.9789                    True
    NP_005500    NonCon           0.0000                   False
    NP_055852    NonCon    10933217.7090                   False
    ------------------------------------------------------------

Exporting
^^^^^^^^^

Writing delimited formats
"""""""""""""""""""""""""

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> table.writeToFile('stats_tab.txt', sep='\t')

Writing latex format
""""""""""""""""""""

It is also possible to specify column alignment, table caption and other arguments.

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> print table.tostring(format='latex')
    \begin{longtable}[htp!]{ r r r }
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
    \end{longtable}

.. we remove the table data

.. doctest::
    :hide:

    >>> import os
    >>> os.remove('stats.txt')
    >>> os.remove('stats_tab.txt')
