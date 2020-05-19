Loading a csv file
^^^^^^^^^^^^^^^^^^

.. authors, Gavin Huttley

We load a tab separated data file using the ``load_table()`` function. The format is inferred from the filename suffix.

.. doctest::

    >>> from cogent3 import load_table
    >>> table = load_table("data/stats.tsv")
    >>> print(table)
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_003077       Con           2.5386
    NP_004893       Con      121351.4264
    NP_005079       Con     9516594.9789
    NP_005500    NonCon           0.0000
    NP_055852    NonCon    10933217.7090
    ------------------------------------

.. note:: The known filename suffixes are ``.csv``, ``.tsv`` and ``.pkl`` or ``.pickle`` (Python's pickle format).

Loading delimited specifying the format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Although unecessary in this case, it's possible to override the suffix by specifying the delimiter using the ``seq`` argument.

.. doctest::

    >>> from cogent3 import load_table
    >>> table = load_table("data/stats.tsv", sep="\t")
    >>> print(table)
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_003077       Con           2.5386
    NP_004893       Con      121351.4264
    NP_005079       Con     9516594.9789
    NP_005500    NonCon           0.0000
    NP_055852    NonCon    10933217.7090
    ------------------------------------

Make a table from header and rows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import make_table
    >>> header = ['A', 'B', 'C']
    >>> rows = [range(3), range(3,6), range(6,9), range(9,12)]
    >>> table = make_table(header=['A', 'B', 'C'], data=rows)
    >>> print(table)
    =============
    A     B     C
    -------------
    0     1     2
    3     4     5
    6     7     8
    9    10    11
    -------------

Make a table from a ``dict``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a ``dict`` with key's as column headers.

.. doctest::

    >>> from cogent3 import make_table
    >>> data = dict(A=[0, 3, 6], B=[1, 4, 7], C=[2, 5, 8])
    >>> table = make_table(data=data)
    >>> print(table)
    ===========
    A    B    C
    -----------
    0    1    2
    3    4    5
    6    7    8
    -----------

Specify the column order when creating from a ``dict``.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = make_table(header=["C", "A", "B"], data=data)
    >>> print(table)
    ===========
    C    A    B
    -----------
    2    0    1
    5    3    4
    8    6    7
    -----------

Create the table with an index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``Table`` can be indexed like a dict if you designate a column as the index (and that column has a unique value for every row).

.. doctest::

    >>> table = load_table("data/stats.tsv", index="Locus")
    >>> print(table["NP_055852"])
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_055852    NonCon    10933217.7090
    ------------------------------------
    >>> table["NP_055852", "Region"]
    'NonCon'

.. note:: The ``index`` argument also applies when using ``make_table()``.

Create a table from a ``pandas`` data frame
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> from pandas import DataFrame
    >>> df = DataFrame(data=[[0, 1], [3, 7]], columns=["a", "b"])
    >>> table = make_table(data_frame=df)
    >>> print(table)
    ======
    a    b
    ------
    0    1
    3    7
    ------
