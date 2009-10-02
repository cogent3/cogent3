Likelihood analysis of multiple loci
====================================

.. sectionauthor:: Gavin Huttley

We want to know whether an exchangeability parameter is different between alignments. We will specify a null model, under which each alignment get's it's own motif probabilities and all alignments share branch lengths and the exchangeability parameter kappa (the transition / transversion ratio). We'll split the example alignment into two-pieces.

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree, LoadTable
    >>> from cogent.evolve.models import HKY85
    >>> from cogent.recalculation.scope import EACH, ALL
    >>> from cogent.maths.stats import chisqprob
    >>> aln = LoadSeqs("data/long_testseqs.fasta")
    >>> half = len(aln)/2
    >>> aln1 = aln[:half]
    >>> aln2 = aln[half:]

We provide names for those alignments, then construct the tree, model instances.

.. doctest::

    >>> loci_names = ["1st-half", "2nd-half"]
    >>> loci = [aln1, aln2]
    >>> tree = LoadTree(tip_names=aln.getSeqNames())
    >>> mod = HKY85()

To make a likelihood function with multiple alignments we provide the list of loci names. We can then specify a parameter (other than length) to be the same across the loci (using the imported ``ALL``) or different for each locus (using ``EACH``). We conduct a LR test as before.

.. doctest::

    >>> lf = mod.makeLikelihoodFunction(tree,loci=loci_names,digits=2,space=3)
    >>> lf.setParamRule("length", is_independent=False)
    >>> lf.setParamRule('kappa', loci = ALL)
    >>> lf.setAlignment(loci)
    >>> lf.optimise(local=True, show_progress=False)
    >>> print lf
    Likelihood Function Table
    =========================
       locus   motif   mprobs
    -------------------------
    1st-half       T     0.22
    1st-half       C     0.18
    1st-half       A     0.38
    1st-half       G     0.21
    2nd-half       T     0.24
    2nd-half       C     0.19
    2nd-half       A     0.35
    2nd-half       G     0.22
    -------------------------
    ==============
    kappa   length
    --------------
     3.98     0.13
    --------------
    >>> all_lnL = lf.getLogLikelihood()
    >>> all_nfp = lf.getNumFreeParams()
    >>> lf.setParamRule('kappa', loci = EACH)
    >>> lf.optimise(local=True, show_progress=False)
    >>> print lf
    Likelihood Function Table
    ================
       locus   kappa
    ----------------
    1st-half    4.33
    2nd-half    3.74
    ----------------
    =========================
       locus   motif   mprobs
    -------------------------
    1st-half       T     0.22
    1st-half       C     0.18
    1st-half       A     0.38
    1st-half       G     0.21
    2nd-half       T     0.24
    2nd-half       C     0.19
    2nd-half       A     0.35
    2nd-half       G     0.22
    -------------------------
    ======
    length
    ------
      0.13
    ------
    >>> each_lnL = lf.getLogLikelihood()
    >>> each_nfp = lf.getNumFreeParams()
    >>> LR = 2 * (each_lnL - all_lnL)
    >>> df = each_nfp - all_nfp

Just to pretty up the result display, I'll print a table consisting of the test statistics created on the fly.

    >>> print LoadTable(header=['LR', 'df', 'p'],
    ...             rows=[[LR, df, chisqprob(LR, df)]], digits=2, space=3)
    ================
      LR   df      p
    ----------------
    1.59    1   0.21
    ----------------
