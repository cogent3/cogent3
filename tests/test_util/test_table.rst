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

The toolkit has a ``Table`` object that can be used for manipulating tabular data. It's properties can be considered like an ordered 2 dimensional dictionary or tuple with flexible output format capabilities of use for exporting data for import into external applications. Importantly, via the restructured text format one can generate html or latex formatted tables. The ``table`` module is located within ``cogent.util``. The ``LoadTable`` convenience function is provided as a top-level ``cogent`` import.

Table creation
--------------

Tables can be created directly using the Table object itself, or a convenience function that handles loading from files. We import both here:

.. doctest::
    
    >>> from cogent import LoadTable
    >>> from cogent.util.table import Table

First, if you try and create a ``Table`` without any data, it raises a ``RuntimeError``.

.. doctest::
    
    >>> t = Table()
    Traceback (most recent call last):
    RuntimeError: header and rows must be provided to Table
    >>> t = Table(header=[], rows=[])
    Traceback (most recent call last):
    RuntimeError: header and rows must be provided to Table

Let's create a very simple, rather nonsensical, table first. To create a table requires a header series, and a 2D series (either of type ``tuple``, ``list``, ``dict``).

.. doctest::
    
    >>> column_headings = ['Journal', 'Impact']

The string "Journal" will become the first column heading, "Impact" the second column heading. The data are,

.. doctest::
    
    >>> rows = [['INT J PARASITOL', 2.9],
    ... ['J MED ENTOMOL', 1.4],
    ... ['Med Vet Entomol', 1.0],
    ... ['INSECT MOL BIOL', 2.85],
    ... ['J AM MOSQUITO CONTR', 0.811],
    ... ['MOL PHYLOGENET EVOL', 2.8],
    ... ['HEREDITY', 1.99e+0],
    ... ['AM J TROP MED HYG', 2.105],
    ... ['MIL MED', 0.605],
    ... ['MED J AUSTRALIA', 1.736]]

We create the simplest of tables.

.. doctest::
    
    >>> t = Table(header = column_headings, rows = rows)
    >>> print t
    =============================
                Journal    Impact
    -----------------------------
        INT J PARASITOL    2.9000
          J MED ENTOMOL    1.4000
        Med Vet Entomol    1.0000
        INSECT MOL BIOL    2.8500
    J AM MOSQUITO CONTR    0.8110
    MOL PHYLOGENET EVOL    2.8000
               HEREDITY    1.9900
      AM J TROP MED HYG    2.1050
                MIL MED    0.6050
        MED J AUSTRALIA    1.7360
    -----------------------------

The format above is referred to as 'simple' format in the documentation. Notice that the numbers in this table have 4 decimal places, despite the fact the original data were largely strings and had ``max`` of 3 decimal places precision. ``Table`` converts string representations of numbers to their appropriate form when you do ``str(table)`` or print the table.

We have several things we might want to specify when creating a table: the precision and or format of floating point numbers (integer argument - ``digits``), the spacing between columns (integer argument or actual string of whitespace - ``space``), title (argument - ``title``), and legend (argument - ``legend``). Lets modify some of these and provide a title and legend.

.. doctest::
    
    >>> t = Table(column_headings, rows, title='Journal impact factors', legend='From ISI',
    ...     digits=2, space='        ')
    >>> print t
    Journal impact factors
    =================================
                Journal        Impact
    ---------------------------------
        INT J PARASITOL          2.90
          J MED ENTOMOL          1.40
        Med Vet Entomol          1.00
        INSECT MOL BIOL          2.85
    J AM MOSQUITO CONTR          0.81
    MOL PHYLOGENET EVOL          2.80
               HEREDITY          1.99
      AM J TROP MED HYG          2.10
                MIL MED          0.60
        MED J AUSTRALIA          1.74
    ---------------------------------
    From ISI

The Table class cannot handle arbitrary python objects, unless they are passed in as strings. Note in this case we now directly pass in the column headings list and the handling of missing data can be explicitly specified..

.. doctest::
    
    >>> t2 = Table(['abcd', 'data'], [[str(range(1,6)), '0'],
    ...                               ['x', 5.0], ['y', None]],
    ...           missing_data='*')
    >>> print t2
    =========================
               abcd      data
    -------------------------
    [1, 2, 3, 4, 5]         0
                  x    5.0000
                  y         *
    -------------------------

Table column headings can be assessed from the ``table.Header`` property

.. doctest::
    
    >>> assert t2.Header == ['abcd', 'data']

and this is immutable (cannot be changed).

.. doctest::
    
    >>> t2.Header[1] = 'Data'
    Traceback (most recent call last):
    RuntimeError: Table Header is immutable, use withNewHeader

If you want to change the Header, use the ``withNewHeader`` method. This can be done one column at a time, or as a batch. The returned Table is identical aside from the modified column labels.

.. doctest::
    
    >>> mod_header = t2.withNewHeader('abcd', 'ABCD')
    >>> assert mod_header.Header == ['ABCD', 'data']
    >>> mod_header = t2.withNewHeader(['abcd', 'data'], ['ABCD', 'DATA'])
    >>> print mod_header
    =========================
               ABCD      DATA
    -------------------------
    [1, 2, 3, 4, 5]         0
                  x    5.0000
                  y         *
    -------------------------

Tables may also be created from 2-dimensional dictionaries. In this case, special capabilities are provided to enforce printing rows in a particular order.

.. doctest::
    
    >>> d2D={'edge.parent': {'NineBande': 'root', 'edge.1': 'root',
    ... 'DogFaced': 'root', 'Human': 'edge.0', 'edge.0': 'edge.1',
    ... 'Mouse': 'edge.1', 'HowlerMon': 'edge.0'}, 'x': {'NineBande': 1.0,
    ... 'edge.1': 1.0, 'DogFaced': 1.0, 'Human': 1.0, 'edge.0': 1.0,
    ... 'Mouse': 1.0, 'HowlerMon': 1.0}, 'length': {'NineBande': 4.0,
    ... 'edge.1': 4.0, 'DogFaced': 4.0, 'Human': 4.0, 'edge.0': 4.0,
    ... 'Mouse': 4.0, 'HowlerMon': 4.0}, 'y': {'NineBande': 3.0, 'edge.1': 3.0,
    ... 'DogFaced': 3.0, 'Human': 3.0, 'edge.0': 3.0, 'Mouse': 3.0,
    ... 'HowlerMon': 3.0}, 'z': {'NineBande': 6.0, 'edge.1': 6.0,
    ... 'DogFaced': 6.0, 'Human': 6.0, 'edge.0': 6.0, 'Mouse': 6.0,
    ... 'HowlerMon': 6.0},
    ... 'edge.name': ['Human', 'HowlerMon', 'Mouse', 'NineBande', 'DogFaced',
    ... 'edge.0', 'edge.1']}
    >>> row_order = d2D['edge.name']
    >>> d2D['edge.name'] = dict(zip(row_order, row_order))
    >>> t3 = Table(['edge.name', 'edge.parent', 'length', 'x', 'y', 'z'], d2D,
    ... row_order = row_order, missing_data='*', space=8, max_width = 50,
    ... row_ids = True, title = 'My Title',
    ... legend = 'Legend: this is a nonsense example.')
    >>> print t3
    My Title
    ==========================================
    edge.name        edge.parent        length
    ------------------------------------------
        Human             edge.0        4.0000
    HowlerMon             edge.0        4.0000
        Mouse             edge.1        4.0000
    NineBande               root        4.0000
     DogFaced               root        4.0000
       edge.0             edge.1        4.0000
       edge.1               root        4.0000
    ------------------------------------------
    <BLANKLINE>
    continued: My Title
    =====================================
    edge.name             x             y
    -------------------------------------
        Human        1.0000        3.0000
    HowlerMon        1.0000        3.0000
        Mouse        1.0000        3.0000
    NineBande        1.0000        3.0000
     DogFaced        1.0000        3.0000
       edge.0        1.0000        3.0000
       edge.1        1.0000        3.0000
    -------------------------------------
    <BLANKLINE>
    continued: My Title
    =======================
    edge.name             z
    -----------------------
        Human        6.0000
    HowlerMon        6.0000
        Mouse        6.0000
    NineBande        6.0000
     DogFaced        6.0000
       edge.0        6.0000
       edge.1        6.0000
    -----------------------
    <BLANKLINE>
    Legend: this is a nonsense example.

In the above we specify a maximum width of the table, and also specify row identifiers (using ``row_ids``, the integer corresponding to the column at which data begin, preceding columns are taken as the identifiers). This has the effect of forcing the table to wrap when the simple text format is used, but wrapping does not occur for any other format. The ``row_ids`` are keys for slicing the table by row, and as identifiers are presented in each wrapped sub-table.

We can also customise the formatting of individual columns.

.. doctest::
    
    >>> rows = (('NP_003077_hs_mm_rn_dna', 'Con', 2.5386013224378985),
    ... ('NP_004893_hs_mm_rn_dna', 'Con', 0.12135142635634111e+06),
    ... ('NP_005079_hs_mm_rn_dna', 'Con', 0.95165949788861326e+07),
    ... ('NP_005500_hs_mm_rn_dna', 'Con', 0.73827030202664901e-07),
    ... ('NP_055852_hs_mm_rn_dna', 'Con', 1.0933217708952725e+07))

We first create a table and show the default formatting behaviour for ``Table``.

.. doctest::
    
    >>> t46 = Table(['Gene', 'Type', 'LR'], rows)
    >>> print t46
    ===============================================
                      Gene    Type               LR
    -----------------------------------------------
    NP_003077_hs_mm_rn_dna     Con           2.5386
    NP_004893_hs_mm_rn_dna     Con      121351.4264
    NP_005079_hs_mm_rn_dna     Con     9516594.9789
    NP_005500_hs_mm_rn_dna     Con           0.0000
    NP_055852_hs_mm_rn_dna     Con    10933217.7090
    -----------------------------------------------

We then format the ``LR`` column to use a scientific number format.

.. doctest::
    
    >>> t46 = Table(['Gene', 'Type', 'LR'], rows)
    >>> t46.setColumnFormat('LR', "%.4e")
    >>> print t46
    ============================================
                      Gene    Type            LR
    --------------------------------------------
    NP_003077_hs_mm_rn_dna     Con    2.5386e+00
    NP_004893_hs_mm_rn_dna     Con    1.2135e+05
    NP_005079_hs_mm_rn_dna     Con    9.5166e+06
    NP_005500_hs_mm_rn_dna     Con    7.3827e-08
    NP_055852_hs_mm_rn_dna     Con    1.0933e+07
    --------------------------------------------

It is safe to directly modify certain attributes, such as the title, legend and white space separating columns, which we do for the ``t46``.

.. doctest::
    
    >>> t46.Title = "A new title"
    >>> t46.Legend = "A new legend"
    >>> t46.Space = '  '
    >>> print t46
    A new title
    ========================================
                      Gene  Type          LR
    ----------------------------------------
    NP_003077_hs_mm_rn_dna   Con  2.5386e+00
    NP_004893_hs_mm_rn_dna   Con  1.2135e+05
    NP_005079_hs_mm_rn_dna   Con  9.5166e+06
    NP_005500_hs_mm_rn_dna   Con  7.3827e-08
    NP_055852_hs_mm_rn_dna   Con  1.0933e+07
    ----------------------------------------
    A new legend

We can provide settings for multiple columns.

.. doctest::
    
    >>> t3 = Table(['edge.name', 'edge.parent', 'length', 'x', 'y', 'z'], d2D,
    ... row_order = row_order)
    >>> t3.setColumnFormat('x', "%.1e")
    >>> t3.setColumnFormat('y', "%.2f")
    >>> print t3
    ===============================================================
    edge.name    edge.parent    length          x       y         z
    ---------------------------------------------------------------
        Human         edge.0    4.0000    1.0e+00    3.00    6.0000
    HowlerMon         edge.0    4.0000    1.0e+00    3.00    6.0000
        Mouse         edge.1    4.0000    1.0e+00    3.00    6.0000
    NineBande           root    4.0000    1.0e+00    3.00    6.0000
     DogFaced           root    4.0000    1.0e+00    3.00    6.0000
       edge.0         edge.1    4.0000    1.0e+00    3.00    6.0000
       edge.1           root    4.0000    1.0e+00    3.00    6.0000
    ---------------------------------------------------------------

In some cases, the contents of a column can be of different types. In this instance, rather than passing a column template we pass a reference to a function that will handle this complexity. To illustrate this we will define a function that formats floating point numbers, but returns everything else as is.

.. doctest::
    
    >>> def formatcol(value):
    ...     if isinstance(value, float):
    ...         val = "%.2f" % value
    ...     else:
    ...         val = str(value)
    ...     return val

We apply this to a table with mixed string, integer and floating point data.

.. doctest::
    
    >>> t6 = Table(['ColHead'], [['a'], [1], [0.3], ['cc']],
    ... column_templates = dict(ColHead=formatcol))
    >>> print t6
    =======
    ColHead
    -------
          a
          1
       0.30
         cc
    -------

Table output
------------

Table can output in multiple formats, including restructured text or 'rest' and delimited. These can be obtained using the ``tostring`` method and ``format`` argument as follows. Using table ``t`` from above,

.. doctest::
    
    >>> print t.tostring(format='rest')
    +------------------------------+
    |    Journal impact factors    |
    +---------------------+--------+
    |             Journal | Impact |
    +=====================+========+
    |     INT J PARASITOL |   2.90 |
    +---------------------+--------+
    |       J MED ENTOMOL |   1.40 |
    +---------------------+--------+
    |     Med Vet Entomol |   1.00 |
    +---------------------+--------+
    |     INSECT MOL BIOL |   2.85 |
    +---------------------+--------+
    | J AM MOSQUITO CONTR |   0.81 |
    +---------------------+--------+
    | MOL PHYLOGENET EVOL |   2.80 |
    +---------------------+--------+
    |            HEREDITY |   1.99 |
    +---------------------+--------+
    |   AM J TROP MED HYG |   2.10 |
    +---------------------+--------+
    |             MIL MED |   0.60 |
    +---------------------+--------+
    |     MED J AUSTRALIA |   1.74 |
    +---------------------+--------+
    | From ISI                     |
    +------------------------------+

Arguments such as ``space`` have no effect in this case. The table may also be written to file in any of the available formats (latex, simple text, html, pickle) or using a custom separator (such as a comma or tab). This makes it convenient to get data into other applications (such as R or a spreadsheet program).

Here is the latex format, note how the title and legend are joined into the latex table caption. We also provide optional arguments for the column alignment (fist column left aligned, second column right aligned and remaining columns centred) and a label for table referencing.

.. doctest::
    
    >>> print t3.tostring(format='tex', justify="lrcccc", label="table:example")
    \begin{longtable}[htp!]{ l r c c c c }
    \hline
    \bf{edge.name} & \bf{edge.parent} & \bf{length} & \bf{x} & \bf{y} & \bf{z} \\
    \hline
    \hline
        Human &      edge.0 & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
    HowlerMon &      edge.0 & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
        Mouse &      edge.1 & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
    NineBande &        root & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
     DogFaced &        root & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
       edge.0 &      edge.1 & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
       edge.1 &        root & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
    \hline
    \label{table:example}
    \end{longtable}

More complex latex table justifying is also possible. Specifying the width of individual columns requires passing in a series (list or tuple) of justification commands. In the following we introduce the command for specific columns widths.

.. doctest::
    
    >>> print t3.tostring(format='tex', justify=["l","p{3cm}","c","c","c","c"])
    \begin{longtable}[htp!]{ l p{3cm} c c c c }
    \hline
    \bf{edge.name} & \bf{edge.parent} & \bf{length} & \bf{x} & \bf{y} & \bf{z} \\
    \hline
    \hline
        Human &      edge.0 & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
    HowlerMon &      edge.0 & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
        Mouse &      edge.1 & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
    NineBande &        root & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
     DogFaced &        root & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
       edge.0 &      edge.1 & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
       edge.1 &        root & 4.0000 & 1.0e+00 & 3.00 & 6.0000 \\
    \hline
    \end{longtable}
    >>> print t3.tostring(sep=',')
    edge.name,edge.parent,length,      x,   y,     z
        Human,     edge.0,4.0000,1.0e+00,3.00,6.0000
    HowlerMon,     edge.0,4.0000,1.0e+00,3.00,6.0000
        Mouse,     edge.1,4.0000,1.0e+00,3.00,6.0000
    NineBande,       root,4.0000,1.0e+00,3.00,6.0000
     DogFaced,       root,4.0000,1.0e+00,3.00,6.0000
       edge.0,     edge.1,4.0000,1.0e+00,3.00,6.0000
       edge.1,       root,4.0000,1.0e+00,3.00,6.0000

You can specify any standard text character that will work with your desired target. Useful separators are tabs ('\\t'), or pipes ('\|'). If ``Table`` encounters any of these characters within a cell, it wraps the cell in quotes -- a standard approach to facilitate import by other applications. We will illustrate this with ``t2``.

.. doctest::
    
    >>> print t2.tostring(sep=', ')
               abcd,   data
    "[1, 2, 3, 4, 5]",      0
                  x, 5.0000
                  y,      *

Note that I introduced an extra space after the column just to make the result more readable in this example.

Test the writing of phylip distance matrix format.

.. doctest::
    
    >>> rows = [['a', '', 0.088337278874079342, 0.18848582712597683,
    ...  0.44084000179091454], ['c', 0.088337278874079342, '',
    ...  0.088337278874079342, 0.44083999937417828], ['b', 0.18848582712597683,
    ...  0.088337278874079342, '', 0.44084000179090932], ['e',
    ...  0.44084000179091454, 0.44083999937417828, 0.44084000179090932, '']]
    >>> header = ['seq1/2', 'a', 'c', 'b', 'e']
    >>> dist = Table(rows = rows, header = header,
    ...  row_ids = True)
    >>> print dist.tostring(format = 'phylip')
       4
    a           0.0000  0.0883  0.1885  0.4408
    c           0.0883  0.0000  0.0883  0.4408
    b           0.1885  0.0883  0.0000  0.4408
    e           0.4408  0.4408  0.4408  0.0000

The ``tostring`` method also provides generic html generation via the restructured text format. The ``toRichHtmlTable`` method can be used to generate the html table element by itself, with greater control over formatting. Specifically, users can provide custom callback functions to the ``row_cell_func`` and ``header_cell_func`` arguments to control in detail the formatting of table elements, or use the simpler dictionary based ``element_formatters`` approach. We use the above ``dist`` table to provide a specific callback that will set the background color for diagonal cells. We first write a function that takes the cell value and coordinates, returning the html formmatted text.

.. doctest::
    
    >>> def format_cell(value, row_num, col_num):
    ...     bgcolor=['', ' bgcolor="#0055ff"'][value=='']
    ...     return '<td%s>%s</td>' % (bgcolor, value)

We then call the method, without this argument, then with it.

.. doctest::
    
    >>> straight_html = dist.toRichHtmlTable()
    >>> print straight_html
    <table><tr><th>seq1/2</th><th>a...
    >>> rich_html = dist.toRichHtmlTable(row_cell_func=format_cell,
    ...                                  compact=False)
    >>> print rich_html
    <table>
    <tr>
    <th>seq1/2</th>
    <th>a</th>
    <th>c</th>
    <th>b</th>
    <th>e</th>
    </tr>
    <tr>
    <td>a</td>
    <td bgcolor="#0055ff"></td>
    <td>0.0883</td>...

Saving a table for reloading
----------------------------

Saving a table object to file for later reloading can be done using the standard ``writeToFile`` method and ``filename`` argument to the ``Table`` constructor, specifying any of the formats supported by ``tostring``. The table loading will recreate a table from raw data located at ``filename``. To illustrate this, we first write out the table ``t3`` in ``pickle`` format and then the table ``t2`` in a csv (comma separated values format).

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> t3 = Table(['edge.name', 'edge.parent', 'length', 'x', 'y', 'z'], d2D,
    ... row_order = row_order, missing_data='*', space=8, max_width = 50,
    ... row_ids = True, title = 'My Title',
    ... legend = 'Legend: this is a nonsense example.')
    >>> t3.writeToFile("t3.pickle")
    >>> t3_loaded = LoadTable(filename = "t3.pickle")
    >>> print t3_loaded
    My Title
    ==========================================
    edge.name        edge.parent        length
    ------------------------------------------
        Human             edge.0        4.0000
    HowlerMon             edge.0        4.0000
        Mouse             edge.1        4.0000
    NineBande               root        4.0000
     DogFaced               root        4.0000
       edge.0             edge.1        4.0000
       edge.1               root        4.0000
    ------------------------------------------
    <BLANKLINE>
    continued: My Title
    =====================================
    edge.name             x             y
    -------------------------------------
        Human        1.0000        3.0000
    HowlerMon        1.0000        3.0000
        Mouse        1.0000        3.0000
    NineBande        1.0000        3.0000
     DogFaced        1.0000        3.0000
       edge.0        1.0000        3.0000
       edge.1        1.0000        3.0000
    -------------------------------------
    <BLANKLINE>
    continued: My Title
    =======================
    edge.name             z
    -----------------------
        Human        6.0000
    HowlerMon        6.0000
        Mouse        6.0000
    NineBande        6.0000
     DogFaced        6.0000
       edge.0        6.0000
       edge.1        6.0000
    -----------------------
    <BLANKLINE>
    Legend: this is a nonsense example.
    >>> t2 = Table(['abcd', 'data'], [[str(range(1,6)), '0'], ['x', 5.0],
    ... ['y', None]], missing_data='*', title = 'A \ntitle')
    >>> t2.writeToFile('t2.csv', sep=',')
    >>> t2_loaded = LoadTable(filename = 't2.csv', header = True, with_title = True,
    ...  sep = ',')
    >>> print t2_loaded
    A 
    title
    =========================
               abcd      data
    -------------------------
    [1, 2, 3, 4, 5]         0
                  x    5.0000
                  y          
    -------------------------

Note the ``missing_data`` attribute is not saved in the delimited format, but is in the ``pickle`` format. In the next case, I'm going to override the digits format on reloading of the table.

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> t2 = Table(['abcd', 'data'], [[str(range(1,6)), '0'], ['x', 5.0],
    ... ['y', None]], missing_data='*', title = 'A \ntitle',
    ... legend = "And\na legend too")
    >>> t2.writeToFile('t2.csv', sep=',')
    >>> t2_loaded = LoadTable(filename = 't2.csv', header = True,
    ... with_title = True, with_legend = True, sep = ',', digits = 2)
    >>> print t2_loaded
    A 
    title
    =======================
               abcd    data
    -----------------------
    [1, 2, 3, 4, 5]       0
                  x    5.00
                  y        
    -----------------------
    And
    a legend too

A few things to note about the delimited file saving: formatting arguments are lost in saving to a delimited format; the ``header`` argument specifies whether the first line of the file should be treated as the header; the ``with_title`` and ``with_legend`` arguments are necessary if the file contains them, otherwise they become the header or part of the table. Importantly, if you wish to preserve numerical precision use the ``pickle`` format.

``cPickle`` can load a useful object from the pickled ``Table`` by itself, without needing to know anything about the ``Table`` class.

.. doctest::
    
    >>> import cPickle
    >>> f = file("t3.pickle")
    >>> pickled = cPickle.load(f)
    >>> f.close()
    >>> pickled.keys()
    ['digits', 'row_ids', 'rows', 'title', 'space', 'max_width', 'header',...
    >>> pickled['rows'][0]
    ['Human', 'edge.0', 4.0, 1.0, 3.0, 6.0]

We can read in a delimited format using a custom reader, which we'll now import. We convert columns 2-5 to floats by specifying a field convertor. We then create a reader, specifying the data (below a list but can be a file) properties. Note that if no convertor is provided all data are returned as strings. We can also provide this reader to the ``Table`` constructor for a more direct way of opening such files. In this case, ``Table`` assumes there is a header row and nothing else.

.. doctest::
    
    >>> from cogent.parse.table import ConvertFields, SeparatorFormatParser
    >>> t3.Title = t3.Legend = None
    >>> comma_sep = t3.tostring(sep=",").splitlines()
    >>> print comma_sep
    ['edge.name,edge.parent,length,     x,     y,     z', '    Human,    ...
    >>> converter = ConvertFields([(2,float), (3,float), (4,float), (5, float)])
    >>> reader = SeparatorFormatParser(with_header=True,converter=converter,
    ...      sep=",")
    >>> comma_sep = [line for line in reader(comma_sep)]
    >>> print comma_sep
    [['edge.name', 'edge.parent', 'length', 'x', 'y', 'z'], ['Human',...
    >>> t3.writeToFile("t3.tab", sep="\t")
    >>> reader = SeparatorFormatParser(with_header=True,converter=converter,
    ...      sep="\t")
    >>> t3a = LoadTable(filename="t3.tab", reader=reader, title="new title",
    ...       space=2)
    >>> print t3a
    new title
    ======================================================
    edge.name  edge.parent  length       x       y       z
    ------------------------------------------------------
        Human       edge.0  4.0000  1.0000  3.0000  6.0000
    HowlerMon       edge.0  4.0000  1.0000  3.0000  6.0000
        Mouse       edge.1  4.0000  1.0000  3.0000  6.0000
    NineBande         root  4.0000  1.0000  3.0000  6.0000
     DogFaced         root  4.0000  1.0000  3.0000  6.0000
       edge.0       edge.1  4.0000  1.0000  3.0000  6.0000
       edge.1         root  4.0000  1.0000  3.0000  6.0000
    ------------------------------------------------------

In the above example, the data type in a column is static, e.g. all values in ``x`` are floats. Rather than providing a custom reader, you can get the ``Table`` to construct such a reader based on the first data row using the ``static_column_types`` argument.

.. doctest::
    
    >>> t3a = LoadTable(filename="t3.tab", static_column_types=True, digits=1,
    ...                 sep='\t')
    >>> print t3a
    =======================================================
    edge.name    edge.parent    length      x      y      z
    -------------------------------------------------------
        Human         edge.0       4.0    1.0    3.0    6.0
    HowlerMon         edge.0       4.0    1.0    3.0    6.0
        Mouse         edge.1       4.0    1.0    3.0    6.0
    NineBande           root       4.0    1.0    3.0    6.0
     DogFaced           root       4.0    1.0    3.0    6.0
       edge.0         edge.1       4.0    1.0    3.0    6.0
       edge.1           root       4.0    1.0    3.0    6.0
    -------------------------------------------------------

If you invoke the ``static_column_types`` argument and the column data are not static, you'll get a ``ValueError``. We show this by first creating a simple table with mixed data types in a column, write to file and then try to load with  ``static_column_types=True``.

.. doctest::
    
    >>> t3b = LoadTable(header=['A', 'B'], rows=[[1,1], ['a', 2]], sep=2)
    >>> print t3b
    ======
    A    B
    ------
    1    1
    a    2
    ------
    >>> t3b.writeToFile('test3b.txt', sep='\t')
    >>> t3b = LoadTable('test3b.txt', sep = '\t', static_column_types=True)
    Traceback (most recent call last):
    ValueError: invalid literal for int() with base 10: 'a'

We also test the reader function for a tab delimited format with missing data at the end.

.. doctest::
    
    >>> data = ['ab\tcd\t', 'ab\tcd\tef']
    >>> tab_reader = SeparatorFormatParser(sep='\t')
    >>> for line in tab_reader(data):
    ...     assert len(line) == 3, line

We can likewise specify a writer, using a custom field formatter and provide this to the ``Table`` directly for writing. We first illustrate how the writer works to generate output. We then use it to escape some text fields in quotes. In order to read that back in, we define a custom reader that strips these quotes off.

.. doctest::
    
    >>> from cogent.format.table import FormatFields, SeparatorFormatWriter
    >>> formatter = FormatFields([(0,'"%s"'), (1,'"%s"')])
    >>> writer = SeparatorFormatWriter(formatter=formatter, sep=" | ")
    >>> for formatted in writer(comma_sep, has_header=True):
    ...      print formatted
    edge.name | edge.parent | length | x | y | z
    "Human" | "edge.0" | 4.0 | 1.0 | 3.0 | 6.0
    "HowlerMon" | "edge.0" | 4.0 | 1.0 | 3.0 | 6.0
    "Mouse" | "edge.1" | 4.0 | 1.0 | 3.0 | 6.0
    "NineBande" | "root" | 4.0 | 1.0 | 3.0 | 6.0
    "DogFaced" | "root" | 4.0 | 1.0 | 3.0 | 6.0
    "edge.0" | "edge.1" | 4.0 | 1.0 | 3.0 | 6.0
    "edge.1" | "root" | 4.0 | 1.0 | 3.0 | 6.0
    >>> t3.writeToFile(filename="t3.tab", writer=writer)
    >>> strip = lambda x: x.replace('"', '')
    >>> converter = ConvertFields([(0,strip), (1, strip)])
    >>> reader = SeparatorFormatParser(with_header=True, converter=converter,
    ...       sep="|", strip_wspace=True)
    >>> t3a = LoadTable(filename="t3.tab", reader=reader, title="new title",
    ...       space=2)
    >>> print t3a
    new title
    =============================================
    edge.name  edge.parent  length    x    y    z
    ---------------------------------------------
        Human       edge.0     4.0  1.0  3.0  6.0
    HowlerMon       edge.0     4.0  1.0  3.0  6.0
        Mouse       edge.1     4.0  1.0  3.0  6.0
    NineBande         root     4.0  1.0  3.0  6.0
     DogFaced         root     4.0  1.0  3.0  6.0
       edge.0       edge.1     4.0  1.0  3.0  6.0
       edge.1         root     4.0  1.0  3.0  6.0
    ---------------------------------------------

.. note:: There are performance issues for large files. Pickling has proven very slow for saving very large files and introduces significant file size bloat. A simple delimited format is much more efficient both storage wise and, if you use a custom reader (or specify ``static_column_types=True``), to generate and read. A custom reader was approximately 6 fold faster than the standard delimited file reader.

Table slicing and iteration
---------------------------

The Table class is capable of slicing by row, range of rows, column or range of columns headings or used to identify a single cell. Slicing using the method ``getColumns`` can also be used to reorder columns. In the case of columns, either the string headings or their position integers can be used. For rows, if ``row_ids`` was specified as ``True`` the 0'th cell in each row can also be used.

.. doctest::
    
    >>> t4 = Table(['edge.name', 'edge.parent', 'length', 'x', 'y', 'z'], d2D,
    ... row_order = row_order, row_ids = True, title = 'My Title')

We subset ``t4`` by column and reorder them.

.. doctest::
    
    >>> new = t4.getColumns(['z', 'y'])
    >>> print new
    My Title
    =============================
    edge.name         z         y
    -----------------------------
        Human    6.0000    3.0000
    HowlerMon    6.0000    3.0000
        Mouse    6.0000    3.0000
    NineBande    6.0000    3.0000
     DogFaced    6.0000    3.0000
       edge.0    6.0000    3.0000
       edge.1    6.0000    3.0000
    -----------------------------

We use the column position indexes to do get the same table.

.. doctest::
    
    >>> new = t4.getColumns([5, 4])
    >>> print new
    My Title
    =============================
    edge.name         z         y
    -----------------------------
        Human    6.0000    3.0000
    HowlerMon    6.0000    3.0000
        Mouse    6.0000    3.0000
    NineBande    6.0000    3.0000
     DogFaced    6.0000    3.0000
       edge.0    6.0000    3.0000
       edge.1    6.0000    3.0000
    -----------------------------

We can also using more general slicing, by both rows and columns. The following returns all rows from 4 on, and columns up to (but excluding) 'y':

.. doctest::
    
    >>> k = t4[4:, :'y']
    >>> print k
    My Title
    ============================================
    edge.name    edge.parent    length         x
    --------------------------------------------
     DogFaced           root    4.0000    1.0000
       edge.0         edge.1    4.0000    1.0000
       edge.1           root    4.0000    1.0000
    --------------------------------------------

We can explicitly reference individual cells, in this case using both row and column keys.

.. doctest::
    
    >>> val = t4['HowlerMon', 'y']
    >>> print val
    3.0

We slice a single row,

.. doctest::
    
    >>> new = t4[3]
    >>> print new
    My Title
    ================================================================
    edge.name    edge.parent    length         x         y         z
    ----------------------------------------------------------------
    NineBande           root    4.0000    1.0000    3.0000    6.0000
    ----------------------------------------------------------------

and range of rows.

.. doctest::
    
    >>> new = t4[3:6]
    >>> print new
    My Title
    ================================================================
    edge.name    edge.parent    length         x         y         z
    ----------------------------------------------------------------
    NineBande           root    4.0000    1.0000    3.0000    6.0000
     DogFaced           root    4.0000    1.0000    3.0000    6.0000
       edge.0         edge.1    4.0000    1.0000    3.0000    6.0000
    ----------------------------------------------------------------

You can get disjoint rows.

.. doctest::
    
    >>> print t4.getDisjointRows(['Human', 'Mouse', 'DogFaced'])
    My Title
    ================================================================
    edge.name    edge.parent    length         x         y         z
    ----------------------------------------------------------------
        Human         edge.0    4.0000    1.0000    3.0000    6.0000
        Mouse         edge.1    4.0000    1.0000    3.0000    6.0000
     DogFaced           root    4.0000    1.0000    3.0000    6.0000
    ----------------------------------------------------------------

You can iterate over the table one row at a time and slice the rows. We illustrate this for slicing a single column,

.. doctest::
    
    >>> for row in t:
    ...     print row['Journal']
    INT J PARASITOL
    J MED ENTOMOL
    Med Vet Entomol
    INSECT MOL BIOL
    J AM MOSQUITO CONTR
    MOL PHYLOGENET EVOL
    HEREDITY
    AM J TROP MED HYG
    MIL MED
    MED J AUSTRALIA

and for multiple columns.

.. doctest::
    
    >>> for row in t:
    ...     print row['Journal'], row['Impact']
    INT J PARASITOL 2.9
    J MED ENTOMOL 1.4
    Med Vet Entomol 1.0
    INSECT MOL BIOL 2.85
    J AM MOSQUITO CONTR 0.811
    MOL PHYLOGENET EVOL 2.8
    HEREDITY 1.99
    AM J TROP MED HYG 2.105
    MIL MED 0.605
    MED J AUSTRALIA 1.736

The numerical slice equivalent to the first case above would be ``row[0]``, to the second case either ``row[:]``, ``row[:2]``.

Filtering tables - selecting subsets of rows/columns
----------------------------------------------------

We want to be able to slice a table, based on some condition(s), to produce a new subset table. For instance, we construct a table with type and probability values.

.. doctest::
    
    >>> header = ['Gene', 'type', 'LR', 'df', 'Prob']
    >>> rows = (('NP_003077_hs_mm_rn_dna', 'Con', 2.5386, 1, 0.1111),
    ...         ('NP_004893_hs_mm_rn_dna', 'Con', 0.1214, 1, 0.7276),
    ...         ('NP_005079_hs_mm_rn_dna', 'Con', 0.9517, 1, 0.3293),
    ...         ('NP_005500_hs_mm_rn_dna', 'Con', 0.7383, 1, 0.3902),
    ...         ('NP_055852_hs_mm_rn_dna', 'Con', 0.0000, 1, 0.9997),
    ...         ('NP_057012_hs_mm_rn_dna', 'Unco', 34.3081, 1, 0.0000),
    ...         ('NP_061130_hs_mm_rn_dna', 'Unco', 3.7986, 1, 0.0513),
    ...         ('NP_065168_hs_mm_rn_dna', 'Con', 89.9766, 1, 0.0000),
    ...         ('NP_065396_hs_mm_rn_dna', 'Unco', 11.8912, 1, 0.0006),
    ...         ('NP_109590_hs_mm_rn_dna', 'Con', 0.2121, 1, 0.6451),
    ...         ('NP_116116_hs_mm_rn_dna', 'Unco', 9.7474, 1, 0.0018))
    >>> t5 = Table(header, rows)
    >>> print t5
    =========================================================
                      Gene    type         LR    df      Prob
    ---------------------------------------------------------
    NP_003077_hs_mm_rn_dna     Con     2.5386     1    0.1111
    NP_004893_hs_mm_rn_dna     Con     0.1214     1    0.7276
    NP_005079_hs_mm_rn_dna     Con     0.9517     1    0.3293
    NP_005500_hs_mm_rn_dna     Con     0.7383     1    0.3902
    NP_055852_hs_mm_rn_dna     Con     0.0000     1    0.9997
    NP_057012_hs_mm_rn_dna    Unco    34.3081     1    0.0000
    NP_061130_hs_mm_rn_dna    Unco     3.7986     1    0.0513
    NP_065168_hs_mm_rn_dna     Con    89.9766     1    0.0000
    NP_065396_hs_mm_rn_dna    Unco    11.8912     1    0.0006
    NP_109590_hs_mm_rn_dna     Con     0.2121     1    0.6451
    NP_116116_hs_mm_rn_dna    Unco     9.7474     1    0.0018
    ---------------------------------------------------------

We then seek to obtain only those rows that contain probabilities < 0.05. We use valid python code within a string. **Note:** Make sure your column headings could be valid python variable names or the string based approach will fail (you could use an external function instead, see below).

.. doctest::
    
    >>> sub_table1 = t5.filtered(callback = "Prob < 0.05")
    >>> print sub_table1
    =========================================================
                      Gene    type         LR    df      Prob
    ---------------------------------------------------------
    NP_057012_hs_mm_rn_dna    Unco    34.3081     1    0.0000
    NP_065168_hs_mm_rn_dna     Con    89.9766     1    0.0000
    NP_065396_hs_mm_rn_dna    Unco    11.8912     1    0.0006
    NP_116116_hs_mm_rn_dna    Unco     9.7474     1    0.0018
    ---------------------------------------------------------

Using the above table we test the function to extract the raw data for a single column,

.. doctest::
    
    >>> raw = sub_table1.getRawData('LR')
    >>> raw
    [34.308100000000003, 89.976600000000005, 11.8912, 9.7474000000000007]

and from multiple columns.

.. doctest::
    
    >>> raw = sub_table1.getRawData(columns = ['LR', 'df', 'Prob'])
    >>> raw
    [[34.308100000000003, 1, 0.0], [89.976600000000005, 1, 0.0],...

We can also do filtering using an external function, in this case we use a ``lambda`` to obtain only those rows of type 'Unco' that contain probabilities < 0.05, modifying our callback function.

.. doctest::
    
    >>> func = lambda (ty, pr): ty == 'Unco' and pr < 0.05
    >>> sub_table2 = t5.filtered(columns = ('type', 'Prob'), callback = func)
    >>> print sub_table2
    =========================================================
                      Gene    type         LR    df      Prob
    ---------------------------------------------------------
    NP_057012_hs_mm_rn_dna    Unco    34.3081     1    0.0000
    NP_065396_hs_mm_rn_dna    Unco    11.8912     1    0.0006
    NP_116116_hs_mm_rn_dna    Unco     9.7474     1    0.0018
    ---------------------------------------------------------

This can also be done using the string approach.

.. doctest::
    
    >>> sub_table2 = t5.filtered(callback = "type == 'Unco' and Prob < 0.05")
    >>> print sub_table2
    =========================================================
                      Gene    type         LR    df      Prob
    ---------------------------------------------------------
    NP_057012_hs_mm_rn_dna    Unco    34.3081     1    0.0000
    NP_065396_hs_mm_rn_dna    Unco    11.8912     1    0.0006
    NP_116116_hs_mm_rn_dna    Unco     9.7474     1    0.0018
    ---------------------------------------------------------

We can also filter table columns using ``filteredByColumn``. Say we only want the numerical columns, we can write a callback that returns ``False`` if some numerical operation fails, ``True`` otherwise.

.. doctest::
    
    >>> def is_numeric(values):
    ...     try:
    ...         sum(values)
    ...     except TypeError:
    ...         return False
    ...     return True
    >>> print t5.filteredByColumn(callback=is_numeric)
    =======================
         LR    df      Prob
    -----------------------
     2.5386     1    0.1111
     0.1214     1    0.7276
     0.9517     1    0.3293
     0.7383     1    0.3902
     0.0000     1    0.9997
    34.3081     1    0.0000
     3.7986     1    0.0513
    89.9766     1    0.0000
    11.8912     1    0.0006
     0.2121     1    0.6451
     9.7474     1    0.0018
    -----------------------

Appending tables
----------------

Tables may also be appended to each other, to make larger tables. We'll construct two simple tables to illustrate this.

.. doctest::
    
    >>> geneA = Table(['edge.name', 'edge.parent', 'z'], [['Human','root',
    ... 6.0],['Mouse','root', 6.0], ['Rat','root', 6.0]],
    ... title='Gene A')
    >>> geneB = Table(['edge.name', 'edge.parent', 'z'], [['Human','root',
    ... 7.0],['Mouse','root', 7.0], ['Rat','root', 7.0]],
    ... title='Gene B')
    >>> print geneB
    Gene B
    ==================================
    edge.name    edge.parent         z
    ----------------------------------
        Human           root    7.0000
        Mouse           root    7.0000
          Rat           root    7.0000
    ----------------------------------

we now use the ``appended`` Table method to create a new table, specifying that we want a new column created (by passing the ``new_column`` argument a heading) in which the table titles will be placed.

.. doctest::
    
    >>> new = geneA.appended('Gene', geneB, title='Appended tables')
    >>> print new
    Appended tables
    ============================================
      Gene    edge.name    edge.parent         z
    --------------------------------------------
    Gene A        Human           root    6.0000
    Gene A        Mouse           root    6.0000
    Gene A          Rat           root    6.0000
    Gene B        Human           root    7.0000
    Gene B        Mouse           root    7.0000
    Gene B          Rat           root    7.0000
    --------------------------------------------

We repeat this without adding a new column.

.. doctest::
    
    >>> new = geneA.appended(None, geneB, title="Appended, no new column")
    >>> print new
    Appended, no new column
    ==================================
    edge.name    edge.parent         z
    ----------------------------------
        Human           root    6.0000
        Mouse           root    6.0000
          Rat           root    6.0000
        Human           root    7.0000
        Mouse           root    7.0000
          Rat           root    7.0000
    ----------------------------------

Miscellaneous
-------------

Tables have a ``Shape`` attribute, which specifies *x* (number of columns) and *y* (number of rows). The attribute is a tuple and we illustrate it for the above ``sub_table`` tables. Combined with the ``filtered`` method, this attribute can tell you how many rows satisfy a specific condition.

.. doctest::
    
    >>> t5.Shape
    (11, 5)
    >>> sub_table1.Shape
    (4, 5)
    >>> sub_table2.Shape
    (3, 5)

For instance, 3 of the 11 rows in ``t`` were significant and belonged to the ``Unco`` type.

For completeness, we generate a table with no rows and assess its shape.

.. doctest::
    
    >>> func = lambda (ty, pr): ty == 'Unco' and pr > 0.1
    >>> sub_table3 = t5.filtered(columns = ('type', 'Prob'), callback = func)
    >>> sub_table3.Shape
    (0, 5)

The distinct values can be obtained for a single column,

.. doctest::
    
    >>> distinct = new.getDistinctValues("edge.name")
    >>> assert distinct == set(['Rat', 'Mouse', 'Human'])

or multiple columns

.. doctest::
    
    >>> distinct = new.getDistinctValues(["edge.parent", "z"])
    >>> assert distinct == set([('root', 6.0), ('root', 7.0)]), distinct

We can compute column sums. Assuming only numerical values in a column.

.. doctest::
    
    >>> assert new.summed('z') == 39., new.summed('z')

We construct an example with mixed numerical and non-numerical data. We now compute the column sum with mixed non-numerical/numerical data.

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> mix = LoadTable(header=['A', 'B'], rows=[[0,''],[1,2],[3,4]])
    >>> print mix
    ======
    A    B
    ------
    0     
    1    2
    3    4
    ------
    >>> mix.summed('B', strict=False)
    6

We also compute row sums for the pure numerical and mixed non-numerical/numerical rows. For summing across rows we must specify the actual row index as an ``int``.

.. doctest::
    
    >>> mix.summed(0, col_sum=False, strict=False)
    0
    >>> mix.summed(1, col_sum=False)
    3

We can compute the totals for all columns or rows too.

.. doctest::
    
    >>> mix.summed(strict=False)
    [4, 6]
    >>> mix.summed(col_sum=False, strict=False)
    [0, 3, 7]

It is not currently possible to do a subset of columns/rows. We show this for rows here.

.. doctest::
    
    >>> mix.summed([0, 2], col_sum=False, strict=False)
    Traceback (most recent call last):
    RuntimeError: unknown indices type: [0, 2]

We test these for a strictly numerical table.

.. doctest::
    
    >>> non_mix = LoadTable(header=['A', 'B'], rows=[[0,1],[1,2],[3,4]])
    >>> non_mix.summed()
    [4, 7]
    >>> non_mix.summed(col_sum=False)
    [1, 3, 7]

We can normalise a numerical table by row,

.. doctest::
    
    >>> print non_mix.normalized(by_row=True)
    ================
         A         B
    ----------------
    0.0000    1.0000
    0.3333    0.6667
    0.4286    0.5714
    ----------------

or by column, such that the row/column sums are 1.

.. doctest::
    
    >>> print non_mix.normalized(by_row=False)
    ================
         A         B
    ----------------
    0.0000    0.1429
    0.2500    0.2857
    0.7500    0.5714
    ----------------

We normalize by an arbitrary function (maximum value) by row,

.. doctest::
    
    >>> print non_mix.normalized(by_row=True, denominator_func=max)
    ================
         A         B
    ----------------
    0.0000    1.0000
    0.5000    1.0000
    0.7500    1.0000
    ----------------

by column.

.. doctest::
    
    >>> print non_mix.normalized(by_row=False, denominator_func=max)
    ================
         A         B
    ----------------
    0.0000    0.2500
    0.3333    0.5000
    1.0000    1.0000
    ----------------

Extending tables
----------------

In some cases it is desirable to compute an additional column from existing column values. This is done using the ``withNewColumn`` method. We'll use t4 from above, adding two of the columns to create an additional column.

.. doctest::
    
    >>> t7 = t4.withNewColumn('Sum', callback="z+x", digits=2)
    >>> print t7
    My Title
    ==================================================================
    edge.name    edge.parent    length       x       y       z     Sum
    ------------------------------------------------------------------
        Human         edge.0      4.00    1.00    3.00    6.00    7.00
    HowlerMon         edge.0      4.00    1.00    3.00    6.00    7.00
        Mouse         edge.1      4.00    1.00    3.00    6.00    7.00
    NineBande           root      4.00    1.00    3.00    6.00    7.00
     DogFaced           root      4.00    1.00    3.00    6.00    7.00
       edge.0         edge.1      4.00    1.00    3.00    6.00    7.00
       edge.1           root      4.00    1.00    3.00    6.00    7.00
    ------------------------------------------------------------------

We test this with an externally defined function.

.. doctest::
    
    >>> func = lambda (x, y): x * y
    >>> t7 = t4.withNewColumn('Sum', callback=func, columns=("y","z"),
    ... digits=2)
    >>> print t7
    My Title
    ===================================================================
    edge.name    edge.parent    length       x       y       z      Sum
    -------------------------------------------------------------------
        Human         edge.0      4.00    1.00    3.00    6.00    18.00
    HowlerMon         edge.0      4.00    1.00    3.00    6.00    18.00
        Mouse         edge.1      4.00    1.00    3.00    6.00    18.00
    NineBande           root      4.00    1.00    3.00    6.00    18.00
     DogFaced           root      4.00    1.00    3.00    6.00    18.00
       edge.0         edge.1      4.00    1.00    3.00    6.00    18.00
       edge.1           root      4.00    1.00    3.00    6.00    18.00
    -------------------------------------------------------------------
    >>> func = lambda x: x**3
    >>> t7 = t4.withNewColumn('Sum', callback=func, columns="y", digits=2)
    >>> print t7
    My Title
    ===================================================================
    edge.name    edge.parent    length       x       y       z      Sum
    -------------------------------------------------------------------
        Human         edge.0      4.00    1.00    3.00    6.00    27.00
    HowlerMon         edge.0      4.00    1.00    3.00    6.00    27.00
        Mouse         edge.1      4.00    1.00    3.00    6.00    27.00
    NineBande           root      4.00    1.00    3.00    6.00    27.00
     DogFaced           root      4.00    1.00    3.00    6.00    27.00
       edge.0         edge.1      4.00    1.00    3.00    6.00    27.00
       edge.1           root      4.00    1.00    3.00    6.00    27.00
    -------------------------------------------------------------------

Sorting tables
--------------

We want a table sorted according to values in a column.

.. doctest::
    
    >>> sorted = t5.sorted(columns = 'LR')
    >>> print sorted
    =========================================================
                      Gene    type         LR    df      Prob
    ---------------------------------------------------------
    NP_055852_hs_mm_rn_dna     Con     0.0000     1    0.9997
    NP_004893_hs_mm_rn_dna     Con     0.1214     1    0.7276
    NP_109590_hs_mm_rn_dna     Con     0.2121     1    0.6451
    NP_005500_hs_mm_rn_dna     Con     0.7383     1    0.3902
    NP_005079_hs_mm_rn_dna     Con     0.9517     1    0.3293
    NP_003077_hs_mm_rn_dna     Con     2.5386     1    0.1111
    NP_061130_hs_mm_rn_dna    Unco     3.7986     1    0.0513
    NP_116116_hs_mm_rn_dna    Unco     9.7474     1    0.0018
    NP_065396_hs_mm_rn_dna    Unco    11.8912     1    0.0006
    NP_057012_hs_mm_rn_dna    Unco    34.3081     1    0.0000
    NP_065168_hs_mm_rn_dna     Con    89.9766     1    0.0000
    ---------------------------------------------------------

We want a table sorted according to values in a subset of columns, note the order of columns determines the sort order.

.. doctest::
    
    >>> sorted = t5.sorted(columns=('LR', 'type'))
    >>> print sorted
    =========================================================
                      Gene    type         LR    df      Prob
    ---------------------------------------------------------
    NP_055852_hs_mm_rn_dna     Con     0.0000     1    0.9997
    NP_004893_hs_mm_rn_dna     Con     0.1214     1    0.7276
    NP_109590_hs_mm_rn_dna     Con     0.2121     1    0.6451
    NP_005500_hs_mm_rn_dna     Con     0.7383     1    0.3902
    NP_005079_hs_mm_rn_dna     Con     0.9517     1    0.3293
    NP_003077_hs_mm_rn_dna     Con     2.5386     1    0.1111
    NP_061130_hs_mm_rn_dna    Unco     3.7986     1    0.0513
    NP_116116_hs_mm_rn_dna    Unco     9.7474     1    0.0018
    NP_065396_hs_mm_rn_dna    Unco    11.8912     1    0.0006
    NP_057012_hs_mm_rn_dna    Unco    34.3081     1    0.0000
    NP_065168_hs_mm_rn_dna     Con    89.9766     1    0.0000
    ---------------------------------------------------------

We now do a sort based on 2 columns.

.. doctest::
    
    >>> sorted = t5.sorted(columns=('type', 'LR'))
    >>> print sorted
    =========================================================
                      Gene    type         LR    df      Prob
    ---------------------------------------------------------
    NP_055852_hs_mm_rn_dna     Con     0.0000     1    0.9997
    NP_004893_hs_mm_rn_dna     Con     0.1214     1    0.7276
    NP_109590_hs_mm_rn_dna     Con     0.2121     1    0.6451
    NP_005500_hs_mm_rn_dna     Con     0.7383     1    0.3902
    NP_005079_hs_mm_rn_dna     Con     0.9517     1    0.3293
    NP_003077_hs_mm_rn_dna     Con     2.5386     1    0.1111
    NP_065168_hs_mm_rn_dna     Con    89.9766     1    0.0000
    NP_061130_hs_mm_rn_dna    Unco     3.7986     1    0.0513
    NP_116116_hs_mm_rn_dna    Unco     9.7474     1    0.0018
    NP_065396_hs_mm_rn_dna    Unco    11.8912     1    0.0006
    NP_057012_hs_mm_rn_dna    Unco    34.3081     1    0.0000
    ---------------------------------------------------------

Reverse sort a single column

.. doctest::
    
    >>> sorted = t5.sorted('LR', reverse = 'LR')
    >>> print sorted
    =========================================================
                      Gene    type         LR    df      Prob
    ---------------------------------------------------------
    NP_065168_hs_mm_rn_dna     Con    89.9766     1    0.0000
    NP_057012_hs_mm_rn_dna    Unco    34.3081     1    0.0000
    NP_065396_hs_mm_rn_dna    Unco    11.8912     1    0.0006
    NP_116116_hs_mm_rn_dna    Unco     9.7474     1    0.0018
    NP_061130_hs_mm_rn_dna    Unco     3.7986     1    0.0513
    NP_003077_hs_mm_rn_dna     Con     2.5386     1    0.1111
    NP_005079_hs_mm_rn_dna     Con     0.9517     1    0.3293
    NP_005500_hs_mm_rn_dna     Con     0.7383     1    0.3902
    NP_109590_hs_mm_rn_dna     Con     0.2121     1    0.6451
    NP_004893_hs_mm_rn_dna     Con     0.1214     1    0.7276
    NP_055852_hs_mm_rn_dna     Con     0.0000     1    0.9997
    ---------------------------------------------------------

Reverse sort one column but not another

.. doctest::
    
    >>> sorted = t5.sorted(columns=('type', 'LR'), reverse = 'LR')
    >>> print sorted
    =========================================================
                      Gene    type         LR    df      Prob
    ---------------------------------------------------------
    NP_065168_hs_mm_rn_dna     Con    89.9766     1    0.0000
    NP_003077_hs_mm_rn_dna     Con     2.5386     1    0.1111
    NP_005079_hs_mm_rn_dna     Con     0.9517     1    0.3293
    NP_005500_hs_mm_rn_dna     Con     0.7383     1    0.3902
    NP_109590_hs_mm_rn_dna     Con     0.2121     1    0.6451
    NP_004893_hs_mm_rn_dna     Con     0.1214     1    0.7276
    NP_055852_hs_mm_rn_dna     Con     0.0000     1    0.9997
    NP_057012_hs_mm_rn_dna    Unco    34.3081     1    0.0000
    NP_065396_hs_mm_rn_dna    Unco    11.8912     1    0.0006
    NP_116116_hs_mm_rn_dna    Unco     9.7474     1    0.0018
    NP_061130_hs_mm_rn_dna    Unco     3.7986     1    0.0513
    ---------------------------------------------------------

Reverse sort both columns.

.. doctest::
    
    >>> sorted = t5.sorted(columns=('type', 'LR'), reverse = ('type', 'LR'))
    >>> print sorted
    =========================================================
                      Gene    type         LR    df      Prob
    ---------------------------------------------------------
    NP_057012_hs_mm_rn_dna    Unco    34.3081     1    0.0000
    NP_065396_hs_mm_rn_dna    Unco    11.8912     1    0.0006
    NP_116116_hs_mm_rn_dna    Unco     9.7474     1    0.0018
    NP_061130_hs_mm_rn_dna    Unco     3.7986     1    0.0513
    NP_065168_hs_mm_rn_dna     Con    89.9766     1    0.0000
    NP_003077_hs_mm_rn_dna     Con     2.5386     1    0.1111
    NP_005079_hs_mm_rn_dna     Con     0.9517     1    0.3293
    NP_005500_hs_mm_rn_dna     Con     0.7383     1    0.3902
    NP_109590_hs_mm_rn_dna     Con     0.2121     1    0.6451
    NP_004893_hs_mm_rn_dna     Con     0.1214     1    0.7276
    NP_055852_hs_mm_rn_dna     Con     0.0000     1    0.9997
    ---------------------------------------------------------

Joining Tables
--------------

The Table object is capable of joins or merging of records in two tables. There are two fundamental types of joins -- inner and outer -- with there being different sub-types. We demonstrate these first constructing some simple tables.

.. doctest::
    
    >>> a=Table(header=["index", "col2","col3"],
    ...         rows=[[1,2,3],[2,3,1],[2,6,5]], title="A")
    >>> print a
    A
    =====================
    index    col2    col3
    ---------------------
        1       2       3
        2       3       1
        2       6       5
    ---------------------
    >>> b=Table(header=["index", "col2","col3"],
    ...         rows=[[1,2,3],[2,2,1],[3,6,3]], title="B")
    >>> print b
    B
    =====================
    index    col2    col3
    ---------------------
        1       2       3
        2       2       1
        3       6       3
    ---------------------
    >>> c=Table(header=["index","col_c2"],rows=[[1,2],[3,2],[3,5]],title="C")
    >>> print c
    C
    ===============
    index    col_c2
    ---------------
        1         2
        3         2
        3         5
    ---------------

For a natural inner join, only 1 copy of columns with the same name are retained. So we expect the headings to be identical between the table ``a``/``b`` and the result of ``a.joined(b)`` or ``b.joined(a)``.

.. doctest::
    
    >>> assert a.joined(b).Header == b.Header
    >>> assert b.joined(a).Header == a.Header

For a standard inner join, the joined table should contain all columns from ``a`` and ``b`` excepting the index column(s). Simply providing a column name (or index) selects this behaviour. Note that in this case, column names from the second table are made unique by prefixing them with that tables title. If the provided tables do not have a title, a ``RuntimeError`` is raised.

.. doctest::
    
    >>> b.Title = None
    >>> try:
    ...     a.joined(b)
    ... except RuntimeError:
    ...     pass
    >>> b.Title = 'B'
    >>> assert a.joined(b, "index").Header == ["index", "col2", "col3",
    ...                                        "B_col2", "B_col3"]
    ...                                         

Note that the table title's were used to prefix the column headings from the second table. We further test this using table ``c`` which has different dimensions.

.. doctest::
    
    >>> assert a.joined(c,"index").Header == ["index","col2","col3",
    ...                                       "C_col_c2"]

It's also possible to specify index columns using numerical values, the results of which should be the same.

.. doctest::
    
    >>> assert a.joined(b,[0, 2]).getRawData() ==\
    ...                          a.joined(b,["index","col3"]).getRawData()

Additionally, it's possible to provide two series of indices for the two tables. Here, they have identical values.

.. doctest::
    
    >>> assert a.joined(b, ["index", "col3"],["index", "col3"]).getRawData()\
    ...         == a.joined(b,["index","col3"]).getRawData()

The results of a standard join between tables ``a`` and ``b`` are

.. doctest::
    
    >>> print a.joined(b, ["index"], title='A&B')
    A&B
    =========================================
    index    col2    col3    B_col2    B_col3
    -----------------------------------------
        1       2       3         2         3
        2       3       1         2         1
        2       6       5         2         1
    -----------------------------------------

We demo the table specific indices.

.. doctest::
    
    >>> print a.joined(c, ["col2"], ["index"], title='A&C by "col2/index"')
    A&C by "col2/index"
    =================================
    index    col2    col3    C_col_c2
    ---------------------------------
        2       3       1           2
        2       3       1           5
    ---------------------------------

Tables ``a`` and ``c`` share a single row with the same value in the ``index`` column, hence a join by that index should return a table with just that row.

.. doctest::
    
    >>> print a.joined(c, "index", title='A&C by "index"')
    A&C by "index"
    =================================
    index    col2    col3    C_col_c2
    ---------------------------------
        1       2       3           2
    ---------------------------------

A natural join of tables ``a`` and ``b`` results in a table with only rows that were identical between the two parents.

.. doctest::
    
    >>> print a.joined(b, title='A&B Natural Join')
    A&B Natural Join
    =====================
    index    col2    col3
    ---------------------
        1       2       3
    ---------------------

We test the outer join by defining an additional table with different dimensions, and conducting a join specifying ``inner_join=False``.

.. doctest::
    
    >>> d=Table(header=["index", "col_c2"], rows=[[5,42],[6,23]], title="D")
    >>> print d
    D
    ===============
    index    col_c2
    ---------------
        5        42
        6        23
    ---------------
    >>> print c.joined(d,inner_join=False, title='C&D Outer join')
    C&D Outer join
    ======================================
    index    col_c2    D_index    D_col_c2
    --------------------------------------
        1         2          5          42
        1         2          6          23
        3         2          5          42
        3         2          6          23
        3         5          5          42
        3         5          6          23
    --------------------------------------

We establish the ``joined`` method works for mixtures of character and numerical data, setting some indices and some cell values to be strings.

.. doctest::
    
    >>> a=Table(header=["index", "col2","col3"],
    ...         rows=[[1,2,"3"],["2",3,1],[2,6,5]], title="A")
    >>> b=Table(header=["index", "col2","col3"],
    ...         rows=[[1,2,"3"],["2",2,1],[3,6,3]], title="B")
    >>> assert a.joined(b, ["index", "col3"],["index", "col3"]).getRawData()\
    ...         == a.joined(b,["index","col3"]).getRawData()

We test that the ``joined`` method works when the column index orders differ.

.. doctest::
    
    >>> t1_header = ['a', 'b']
    >>> t1_rows = [(1,2),(3,4)]
    >>> t2_header = ['b', 'c']
    >>> t2_rows = [(3,6),(4,8)]
    >>> t1 = Table(header = t1_header, rows = t1_rows, title='t1')
    >>> t2 = Table(header = t2_header, rows = t2_rows, title='t2')
    >>> t3 = t1.joined(t2, columns_self = ["b"], columns_other = ["b"])
    >>> print t3
    ==============
    a    b    t2_c
    --------------
    3    4       8
    --------------

We then establish that a join with no values does not cause a failure, just returns an empty ``Table``.

.. doctest::
    
    >>> t4_header = ['b', 'c']
    >>> t4_rows = [(5,6),(7,8)]
    >>> t4 = LoadTable(header = t4_header, rows = t4_rows)
    >>> t4.Title = 't4'
    >>> t5 = t1.joined(t4, columns_self = ["b"], columns_other = ["b"])

Transposing a table
-------------------

Tables can be transposed.

.. doctest::
    
    >>> from cogent import LoadTable
    >>> title='#Full OTU Counts'
    >>> header = ['#OTU ID', '14SK041', '14SK802']
    >>> rows = [[-2920, '332', 294], 
    ...         [-1606, '302', 229], 
    ...         [-393, 141, 125], 
    ...         [-2109, 138, 120], 
    ...         [-5439, 104, 117], 
    ...         [-1834, 70, 75], 
    ...         [-18588, 65, 47], 
    ...         [-1350, 60, 113], 
    ...         [-2160, 57, 52], 
    ...         [-11632, 47, 36]]
    >>> table = LoadTable(header=header,rows=rows,title=title)
    >>> print table
    #Full OTU Counts
    =============================
    #OTU ID    14SK041    14SK802
    -----------------------------
      -2920        332        294
      -1606        302        229
       -393        141        125
      -2109        138        120
      -5439        104        117
      -1834         70         75
     -18588         65         47
      -1350         60        113
      -2160         57         52
     -11632         47         36
    -----------------------------

We now transpose this. We require a new column heading for header data and an identifier for which existing column will become the header (default is index 0).

.. doctest::
    
    >>> tp = table.transposed(new_column_name='sample',
    ...             select_as_header='#OTU ID', space=2)
    ...             
    >>> print tp
    ==============================================================================
     sample  -2920  -1606  -393  -2109  -5439  -1834  -18588  -1350  -2160  -11632
    ------------------------------------------------------------------------------
    14SK041    332    302   141    138    104     70      65     60     57      47
    14SK802    294    229   125    120    117     75      47    113     52      36
    ------------------------------------------------------------------------------

We test transposition with default value is the same.

.. doctest::
    
    >>> tp = table.transposed(new_column_name='sample', space=2)
    ...             
    >>> print tp
    ==============================================================================
     sample  -2920  -1606  -393  -2109  -5439  -1834  -18588  -1350  -2160  -11632
    ------------------------------------------------------------------------------
    14SK041    332    302   141    138    104     70      65     60     57      47
    14SK802    294    229   125    120    117     75      47    113     52      36
    ------------------------------------------------------------------------------

We test transposition selecting a different column to become the header.

.. doctest::
    
    >>> tp = table.transposed(new_column_name='sample',
    ...             select_as_header='14SK802', space=2)
    ...             
    >>> print tp
    ==============================================================================
     sample    294    229   125    120    117     75      47    113     52      36
    ------------------------------------------------------------------------------
    #OTU ID  -2920  -1606  -393  -2109  -5439  -1834  -18588  -1350  -2160  -11632
    14SK041    332    302   141    138    104     70      65     60     57      47
    ------------------------------------------------------------------------------

Counting rows
-------------

We can count the number of rows for which a condition holds. This method uses the same arguments as ``filtered`` but returns an integer result only.

.. doctest::
    
    >>> print c.count("col_c2 == 2")
    2
    >>> print c.joined(d,inner_join=False).count("index==3 and D_index==5")
    2

Testing a sub-component
-----------------------

Before using ``Table``, we exercise some formatting code:

.. doctest::
    
    >>> from cogent.format.table import formattedCells, phylipMatrix, latex

We check we can format an arbitrary 2D list, without a header, using the ``formattedCells`` function directly.

.. doctest::
    
    >>> data = [[230, 'acdef', 1.3], [6, 'cc', 1.9876]]
    >>> head = ['one', 'two', 'three']
    >>> header, formatted = formattedCells(data, header = head)
    >>> print formatted
    [['230', 'acdef', '1.3000'], ['  6', '   cc', '1.9876']]
    >>> print header
    ['one', '  two', ' three']

We directly test the latex formatting.

.. doctest::
    
    >>> print latex(formatted, header, justify='lrl', caption='A legend',
    ...             label="table:test")
    \begin{longtable}[htp!]{ l r l }
    \hline
    \bf{one} & \bf{two} & \bf{three} \\
    \hline
    \hline
    230 & acdef & 1.3000 \\
      6 &    cc & 1.9876 \\
    \hline
    \caption{A legend}
    \label{table:test}
    \end{longtable}

..
    Import the ``os`` module so some file cleanup can be done at the end. To check the contents of those files, just delete the following prior to running the test. The try/except clause below is aimed at case where ``junk.pdf`` wasn't created due to ``reportlab`` not being present.

.. doctest::
    :hide:
    
    >>> import os
    >>> to_delete = ['t3.pickle', 't2.csv', 't3.tab', 'test3b.txt']
    >>> for f in to_delete:
    ...     try:
    ...         os.remove(f)
    ...     except OSError:
    ...         pass

