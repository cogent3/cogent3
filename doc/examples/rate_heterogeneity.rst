Analysis of rate heterogeneity
==============================

.. sectionauthor:: Gavin Huttley

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
    >>> aln1 = lf.simulateAlignment(sequence_length=10000)
    >>> lf.setParamRule('length', value=0.2, is_const=True)
    >>> aln2 = lf.simulateAlignment(sequence_length=10000)
    >>> aln3 = aln1 + aln2

Start from scratch, optimising only rates and the rate probability ratio.

.. doctest::

    >>> model = Nucleotide(equal_motif_probs=True, ordered_param="rate",
    ...                    distribution="free")
    >>> lf = model.makeLikelihoodFunction(tree, bins=2, digits=2, space=3)
    >>> lf.setAlignment(aln3)
    >>> lf.optimise(local=True, max_restarts=2, show_progress = False)

We want to know the bin probabilities and the posterior probabilities.

.. doctest::
    
    >>> bprobs = [t for t in lf.getStatistics() if 'bin' in t.Title][0]

Print the ``bprobs`` sorted by ``'rate'`` will generate a table like

.. code-block:: python
    
    bin params
    ====================
     bin   bprobs   rate
    --------------------
    bin0     0.49   0.49
    bin1     0.51   1.48
    --------------------

We'll now use a gamma distribution on the sample alignment, specifying the number of bins as 4. We specify that the bins have equal density using the ``lf.setParamRule('bprobs', is_const=True)`` command.

.. doctest::

    >>> model = Nucleotide(equal_motif_probs=True, ordered_param="rate",
    ...                    distribution="gamma")
    >>> lf = model.makeLikelihoodFunction(tree, bins=4)
    >>> lf.setParamRule('bprobs', is_const=True)
    >>> lf.setAlignment(aln3)
    >>> lf.optimise(local=True, max_restarts=2, show_progress = False)
