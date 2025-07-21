.. jupyter-execute::
    :hide-code:

    import set_working_directory

Use an empirical protein substitution model
===========================================

.. sectionauthor:: Gavin Huttley

This file contains an example of importing an empirically determined protein substitution matrix such as Dayhoff et al 1978 and using it to create a substitution model. The globin alignment is from the PAML distribution.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, make_tree
    from cogent3.evolve.substitution_model import EmpiricalProteinMatrix
    from cogent3.parse.paml_matrix import PamlMatrixParser

Make a tree object.  In this case from a string.

.. jupyter-execute::

    treestring = "(((rabbit,rat),human),goat-cow,marsupial);"
    t = make_tree(treestring)

Import the alignment, explicitly setting the ``moltype`` to be protein

.. jupyter-execute::

    al = load_aligned_seqs("data/abglobin_aa.phylip", moltype="protein")

Open the file that contains the empirical matrix and parse the matrix and frequencies.

.. jupyter-execute::

    matrix_file = open("data/dayhoff.dat")

The ``PamlMatrixParser`` will import the matrix and frequency from files designed for Yang's PAML package.  This format is the lower half of the matrix in three letter amino acid name order white space delineated followed by motif frequencies in the same order.

.. jupyter-execute::

    empirical_matrix, empirical_frequencies = PamlMatrixParser(matrix_file)

Create an Empirical Protein Matrix Substitution model object.  This will take the unscaled empirical matrix and use it and the motif frequencies to create a scaled Q matrix.

.. jupyter-execute::

    sm = EmpiricalProteinMatrix(empirical_matrix, empirical_frequencies)

Make a parameter controller, likelihood function object and optimise.

.. jupyter-execute::

    lf = sm.make_likelihood_function(t)
    lf.set_alignment(al)
    lf.optimise(show_progress=False)
    lf
