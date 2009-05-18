Analysis of rate heterogeneity
==============================

A simple example for analyses involving rate heterogeneity among sites. In this case we will simulate an alignment with two rate categories and then try to recover the rates from the alignment.

.. doctest::

    >>> from cogent.evolve.substitution_model import Nucleotide
    >>> from cogent import LoadTree

Make an alignment with equal split between rates 0.6 and 0.2, and then concatenate them to create a new alignment.

.. doctest::

    >>> model = Nucleotide(equal_motif_probs=True)
    >>> tree = LoadTree("data/test.tree")
    >>> lf = model.makeLikelihoodFunction(tree)
    >>> lf.setParamRule('length', value=0.6, is_const=True)
    >>> aln1 = lf.simulateAlignment(sequence_length=1000)
    >>> lf.setParamRule('length', value=0.2, is_const=True)
    >>> aln2 = lf.simulateAlignment(sequence_length=1000)
    >>> aln3 = aln1 + aln2

Start from scratch, optimising only rates and the rate probability ratio.

.. doctest::

    >>> model = Nucleotide(equal_motif_probs=True, ordered_param="rate",
    ...                    distribution="free")
    >>> lf = model.makeLikelihoodFunction(tree, bins=2)
    >>> lf.setAlignment(aln3)
    >>> lf.optimise(local=True, max_restarts=2, show_progress = False)

We want to know the bin probabilities and the posterior probabilities.

.. doctest::

    >>> bprobs = lf.getParamValue('bprobs')
    >>> pp = lf.getBinProbs()
    >>> for bin in [0,1]:
    ...     p = pp[bin]
    ...     rate = lf.getParamValue('rate', bin='bin%s'%bin)
    ...     print '%.2f of sites have rate %.2f' % (bprobs[bin], rate)
    ...     print 'Avg probs over the fast (%.2f) and slow (%.2f) halves' % \
    ...        (sum(p[:1000])/1000, sum(p[1000:])/1000)
    0.12 of sites have rate 0.22
    Avg probs over the fast (0.05) and slow (0.18) halves
    0.88 of sites have rate 1.10
    Avg probs over the fast (0.95) and slow (0.82) halves

We'll now use a gamma distribution on the sample alignment, specifying the number of bins as 4. We specify that the bins have equal density using the ``lf.setParamRule('bprobs', is_const=True)`` command.

.. doctest::

    >>> model = Nucleotide(equal_motif_probs=True, ordered_param="rate",
    ...                    distribution="gamma")
    >>> lf = model.makeLikelihoodFunction(tree, bins=4)
    >>> lf.setParamRule('bprobs', is_const=True)
    >>> lf.setAlignment(aln3)
    >>> lf.optimise(local=True, max_restarts=2, show_progress = False)
