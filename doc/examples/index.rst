#####################
Cogent Usage Examples
#####################

**************************************
A Note on the Computable Documentation
**************************************

The following examples are all available as standalone text files which can be computed using the Python doctest module. One caveat with these tests is a subset will fail sporadically (or even consistently), although there is nothing 'wrong' with the software. These failures arise because of the typically very small data-sets we use in order for the documentation to compute in a short time period. As a result of their small size, the results from numerical optimisations are volatile and can change from one run to another -- leading to 'failures'. Specific examples that are prone to these problems involve the HMM models, the test of Neutrality, rate heterogeneity, unrestricted nucleotide substitution model and even the simplest example.

Examples not explicitly attributed were authored by Gavin Huttley.

*****************
Data manipulation
*****************

.. toctree::
    :maxdepth: 1

    translate_dna
    seq_features
    complete_seq_features
    reverse_complement
    align_codons_to_protein
    aln_profile
    genetic_code_aa_index
    handling_tabular_data
    query_ncbi
    manipulating_tree_nodes

**********************************
Controlling 3rd party applications
**********************************

.. toctree::
    :maxdepth: 1

    application_controller_framework
    alignment_app_controllers
    phylogeny_app_controllers
    generating_app_commandlines

*********************
General data analysis
*********************

.. toctree::
    :maxdepth: 1

    perform_PCoA_analysis
    perform_nmds
    motif_results

******************
Data Visualisation
******************

.. toctree::
    :maxdepth: 1

    draw_dendrogram
    draw_dotplot

*******************
Modelling Evolution
*******************

.. toctree::
    :maxdepth: 1

    simple
    relative_rate
    neutral_test
    empirical_protein_models
    rate_heterogeneity
    testing_multi_loci
    reuse_results
    unrestricted_nucleotide
    simulate_alignment
    parametric_bootstrap
    estimate_startingpoint
    coevolution
    hmm_par_heterogeneity

***************************
Phylogenetic Reconstruction
***************************

.. toctree::
    :maxdepth: 1

    calculate_pairwise_distances
    calculate_neigbourjoining_tree
    calculate_UPGMA_cluster
    phylo_by_ls
    maketree_from_proteinseqs

