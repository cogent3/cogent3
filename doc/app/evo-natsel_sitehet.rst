.. jupyter-execute::
    :hide-code:

    import set_working_directory

``natsel_sitehet`` – a test of site heterogeneity
=================================================

This app evaluates evidence for whether sites differ in their mode of
natural selection (`Nielsen and Yang
1998 <https://www.ncbi.nlm.nih.gov/pubmed/9539414>`__).

.. jupyter-execute::

    from cogent3.app import evo, io

    loader = io.load_aligned(format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")

    sites_differ = evo.natsel_sitehet(
        "GNC", tree="data/primate_brca1.tree", optimise_motif_probs=False
    )

    result = sites_differ(aln)
    result

The models have been constructed such that site-class bins have names indicating the mode of natural selection: -ve is purifying (omega<1); neutral (omega=1); and +ve is positive natural selection (omega>1). The two parameters of interest relating to these are the ``bprobs`` (the maximum likelihood estimate of the frequency of the site-class) and the corresponding value of omega.

.. jupyter-execute::

    result.alt.lf

Getting the individual site posterior probabilities
---------------------------------------------------

I’m just displaying the posterior-probabilities from the first 20 positions only.

.. jupyter-execute::

    bprobs = result.alt.lf.get_bin_probs()
    bprobs[:, :20]
