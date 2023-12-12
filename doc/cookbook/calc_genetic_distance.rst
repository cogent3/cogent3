.. jupyter-execute::
    :hide-code:

    import set_working_directory

****************************
Genetic distance calculation
****************************

Fast pairwise distance estimation
=================================

For a limited number of evolutionary models a fast implementation is
available.

.. jupyter-execute::

    from cogent3 import available_distances

    available_distances()

Computing genetic distances using the ``Alignment`` object
==========================================================

Abbreviations listed from ``available_distances()`` can be used as values for the ``distance_matrix(calc=<abbreviation>)``.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
    dists = aln.distance_matrix(calc="tn93", show_progress=False)
    dists

Using the distance calculator directly
======================================

.. jupyter-execute::

    from cogent3 import get_distance_calculator, load_aligned_seqs

    aln = load_aligned_seqs("data/primate_brca1.fasta")
    dist_calc = get_distance_calculator("tn93", alignment=aln)
    dist_calc

.. jupyter-execute::

    dist_calc.run(show_progress=False)
    dists = dist_calc.get_pairwise_distances()
    dists

The distance calculation object can provide more information. For instance, the standard errors.

.. jupyter-execute::

    dist_calc.stderr

Likelihood based pairwise distance estimation
=============================================

The standard ``cogent3`` likelihood function can also be used to estimate distances. Because these require numerical optimisation they can be significantly slower than the fast estimation approach above.

The following will use the F81 nucleotide substitution model and perform numerical optimisation.

.. jupyter-execute::

    from cogent3 import get_model, load_aligned_seqs
    from cogent3.evolve import distance

    aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
    d = distance.EstimateDistances(aln, submodel=get_model("F81"))
    d.run(show_progress=False)
    dists = d.get_pairwise_distances()
    dists

*****************************************************
Get the names of sequences with max pairwise distance
*****************************************************

Given a ``DistanceMatrix`` object, finding the sequences that have the maximum pairwise distance is achieved through the ``max_pair`` method. 

.. jupyter-execute::
    :raises:
    
    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
    dists = aln.distance_matrix(calc="tn93", show_progress=False)
    dists.max_pair()

    
To find the maximum distance, index the ``DistanceMatrix`` with the result of ``max_pair``.

.. jupyter-execute::
    :raises:
    
    dists[dists.max_pair()]


*****************************************************
Get the names of sequences with min pairwise distance
*****************************************************

Given a ``DistanceMatrix`` object, finding the sequences that have the minimum pairwise distance is achieved through the ``min_pair`` method. 


.. note:: As the distance between a sequence and itself is zero, and this is not informative, ``min_pair`` will return the smallest distance not on the diagonal.
    
 
.. jupyter-execute::
    :raises:
    
    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
    dists = aln.distance_matrix(calc="tn93", show_progress=False)
    dists.min_pair()

To find the minimum distance, index the ``DistanceMatrix`` with the result of ``min_pair``.

.. jupyter-execute::
    :raises:

    dists[dists.min_pair()]
