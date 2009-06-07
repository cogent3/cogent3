Estimate parameter values using a sampling from a dataset
=========================================================

.. sectionauthor:: Gavin Huttley

This script uses the ``sample`` method of the alignment class to provide an estimate for a two stage optimisation. This allows rapid optimisation of long alignments and complex models with a good chance of arriving at the global maximum for the model and data. Local optimisation of the full alignment may end up in local maximum and for this reason results from this strategy my be inaccurate.

From cogent import all the components we need.

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree
    >>> from cogent.evolve import  substitution_model

Load your alignment, note that if your file ends with a suffix that is the same as it's format (assuming it's a supported format) then you can just give the filename. Otherwise you can specify the format using the format argument.

.. doctest::

    >>> aln = LoadSeqs(filename = "data/long_testseqs.fasta")

Get your tree

.. doctest::

    >>> t = LoadTree(filename = "data/test.tree")

Get the substitution model (defaults to Felsensteins 1981 model)

.. doctest::

    >>> sm = substitution_model.Nucleotide()

Make a likelihood function from a sample of the alignment the ``sample`` method selects the chosen number of bases at random.

.. doctest::

    >>> lf = sm.makeLikelihoodFunction(t)
    >>> lf.setMotifProbsFromData(aln)
    >>> lf.setAlignment(aln.sample(20))

Optimise with the local optimiser

.. doctest::

    >>> lf.optimise(local=True, show_progress=False)

Next use the whole alignment

.. doctest::

    >>> lf.setAlignment(aln)

and the faster Powell optimiser that will only find the best result near the provided starting point

.. doctest::
    
    >>> lf.optimise(local=True, show_progress=False)

.. doctest::
    
    >>> print lf
    Likelihood Function Table
    =============================
         edge    parent    length
    -----------------------------
        Human    edge.0    0.0309
    HowlerMon    edge.0    0.0412
       edge.0    edge.1    0.0359
        Mouse    edge.1    0.2666
       edge.1      root    0.0226
    NineBande      root    0.0895
     DogFaced      root    0.1095
    -----------------------------
    ===============
    motif    mprobs
    ---------------
        T    0.2317
        C    0.1878
        A    0.3681
        G    0.2125
    ---------------
