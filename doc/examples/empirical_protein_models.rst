Use an empirical protein substitution model
===========================================

.. sectionauthor:: Gavin Huttley

This file contains an example of importing an empirically determined protein substitution matrix such as Dayhoff et al 1978 and using it to create a substitution model. The globin alignment is from the PAML distribution.

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree, PROTEIN
    >>> from cogent.evolve.substitution_model import EmpiricalProteinMatrix
    >>> from cogent.parse.paml_matrix import PamlMatrixParser

Make a tree object.  In this case from a string.

.. doctest::

    >>> treestring="(((rabbit,rat),human),goat-cow,marsupial);"
    >>> t = LoadTree(treestring=treestring)

Import the alignment, explicitly setting the ``moltype`` to be protein

.. doctest::

    >>> al = LoadSeqs('data/abglobin_aa.phylip',
    ...                interleaved=True,
    ...                moltype=PROTEIN,
    ...                )

Open the file that contains the empirical matrix and parse the matrix and frequencies.

.. doctest::

    >>> matrix_file = open('data/dayhoff.dat')

The ``PamlMatrixParser`` will import the matrix and frequency from files designed for Yang's PAML package.  This format is the lower half of the matrix in three letter amino acid name order white space delineated followed by motif frequencies in the same order.

.. doctest::

    >>> empirical_matrix, empirical_frequencies = PamlMatrixParser(matrix_file)

Create an Empirical Protein Matrix Substitution model object.  This will take the unscaled empirical matrix and use it and the motif frequencies to create a scaled Q matrix.

.. doctest::

    >>> sm = EmpiricalProteinMatrix(empirical_matrix, empirical_frequencies)

Make a parameter controller, likelihood function object and optimise.

.. doctest::

    >>> lf = sm.makeLikelihoodFunction(t)
    >>> lf.setAlignment(al)
    >>> lf.optimise(show_progress = False)
    >>> print lf.getLogLikelihood()
    -1706...
    >>> print lf
    Likelihood Function Table
    =============================
         edge    parent    length
    -----------------------------
       rabbit    edge.0    0.0785
          rat    edge.0    0.1750
       edge.0    edge.1    0.0324
        human    edge.1    0.0545
       edge.1      root    0.0269
     goat-cow      root    0.0972
    marsupial      root    0.2424
    -----------------------------
    ===============
    motif    mprobs
    ---------------
        A    0.0871
        C    0.0335
        D    0.0469
        E    0.0495
        F    0.0398
        G    0.0886
        H    0.0336
        I    0.0369
        K    0.0805
        L    0.0854
        M    0.0148
        N    0.0404
        P    0.0507
        Q    0.0383
        R    0.0409
        S    0.0696
        T    0.0585
        V    0.0647
        W    0.0105
        Y    0.0299
    ---------------
