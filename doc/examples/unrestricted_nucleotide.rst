Specifying and using an unrestricted nucleotide substitution model
==================================================================

.. sectionauthor:: Gavin Huttley

Do standard ``cogent`` imports.

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree, DNA
    >>> from cogent.evolve.predicate import MotifChange
    >>> from cogent.evolve.substitution_model import Nucleotide

.. don't pollute screen during execution with uninteresting warning

.. doctest::
    :hide:
    
    >>> import warnings
    >>> warnings.filterwarnings("ignore", "Model not reversible")

To specify substitution models we use the ``MotifChange`` class from predicates. In the case of an unrestricted nucleotide model, we specify 11 such ``MotifChanges``, the last possible change being ignored (with the result it is constrained to equal 1, thus calibrating the matrix). Also note that this is a non-reversible model and thus we can't assume the nucleotide frequencies estimated from the alignments are reasonable estimates for the root frequencies. We therefore specify they are to be optimised using ``optimise_motif_probs`` argument.

.. doctest::

    >>> ACTG = list('ACTG')
    >>> preds = [MotifChange(i, j, forward_only=True) for i in ACTG for j in ACTG if i != j]
    >>> del(preds[-1])
    >>> preds
    [A>C, A>T, A>G, C>A, C>T, C>G, T>A, T>C, T>G, G>A, G>C]
    >>> sm = Nucleotide(predicates=preds, recode_gaps=True,
    ...                 optimise_motif_probs=True)
    >>> print sm
    <BLANKLINE>
    Nucleotide ( name = ''; type = 'None'; params = ['A>T', 'C>G', 'T>G',...

We'll illustrate this with a sample alignment and tree in ``data/primate_cdx2_promoter.fasta``.

.. doctest::

    >>> al = LoadSeqs("data/primate_cdx2_promoter.fasta", moltype=DNA)
    >>> al
    3 x 1525 dna alignment: human[AGCGCCCGCGG...], macaque[AGC...
    >>> tr = LoadTree(tip_names=al.Names)
    >>> print tr
    (human,macaque,chimp)root;

We now construct the parameter controller with each predicate constant across the tree, and get the likelihood function calculator.

.. doctest::

    >>> lf = sm.makeLikelihoodFunction(tr)
    >>> lf.setAlignment(al)
    >>> lf.setTablesFormat(digits=2, space=3)
    >>> lf.setName('Unrestricted model')

We want to make the most general continuous time Markov model, which requires the predicates be independent for every edge.

.. doctest::
    
    >>> lf.optimise(local=True, show_progress=False)

In the output from the ``optimise`` call you'll see progress from the simulated annealing optimiser which is used first, and the Powell optimiser which finishes things off.

.. doctest::

    >>> print lf
    Unrestricted model
    ==========================================================================
     A>C    A>G    A>T    C>A    C>G    C>T    G>A    G>C    T>A    T>C    T>G
    --------------------------------------------------------------------------
    0.49   4.88   1.04   2.04   0.99   7.89   9.00   1.55   0.48   5.53   1.57
    --------------------------------------------------------------------------
    =========================
       edge   parent   length
    -------------------------
      human     root     0.00
    macaque     root     0.04
      chimp     root     0.01
    -------------------------
    ==============
    motif   mprobs
    --------------
        T     0.26
        C     0.26
        A     0.24
        G     0.24
    --------------

This data set is very small, so the parameter estimates are poor and hence doing something like allowing the parameters to differ between edges is silly. If you have lots of data it makes sense to allow parameters to differ between edges, which can be specified by modifying the ``lf`` as follows.

.. doctest::

    >>> for pred in preds:
    ...     lf.setParamRule(str(pred), is_independent=True)

You would then re-optimise the model as above.
