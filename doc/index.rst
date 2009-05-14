#####################
Cogent Usage Examples
#####################

.. Contents::

**********
The Readme
**********

.. include:: ../README

**************************************
A Note on the Computable Documentation
**************************************

The following examples are all available as standalone text files which can be computed using the Python doctest module. One caveat with these tests is a subset will fail sporadically (or even consistently), although there is nothing 'wrong' with the software. These failures arise because of the typically very small data-sets we use in order for the documentation to compute in a short time period. As a result of their small size, the results from numerical optimisations are volatile and can change from one run to another -- leading to 'failures'. Specific examples that are prone to these problems involve the HMM models, the test of Neutrality, rate heterogeneity, unrestricted nucleotide substitution model and even the simplest example.

Examples not explicitly attributed were authored by Gavin Huttley.

*****************
Data manipulation
*****************

.. include:: translate_dna.rst
.. include:: seq_features.rst
.. include:: complete_seq_features.rst
.. include:: reverse_complement.rst
.. include:: align_codons_to_protein.rst
.. include:: aln_profile.rst
.. include:: genetic_code_aa_index.rst
.. include:: handling_tabular_data.rst
.. include:: query_ncbi.rst

**********************************
Controlling 3rd party applications
**********************************

.. include:: alignment_app_controllers.rst
.. include:: phylogeny_app_controllers.rst
.. include:: generating_app_commandlines.rst

*********************
General data analysis
*********************

.. include:: perform_PCoA_analysis.rst
.. include:: motif_results.rst

******************
Data Visualisation
******************

.. include:: draw_dendrogram.rst
.. include:: draw_dotplot.rst

*******************
Modelling Evolution
*******************

.. include:: simple.rst
.. include:: relative_rate.rst
.. include:: neutral_test.rst
.. include:: empirical_protein_models.rst
.. include:: rate_heterogeneity.rst
.. include:: testing_multi_loci.rst
.. include:: reuse_results.rst
.. include:: unrestricted_nucleotide.rst
.. include:: simulate_alignment.rst
.. include:: parametric_bootstrap.rst
.. include:: estimate_startingpoint.rst
.. include:: coevolution.rst

***************************
Phylogenetic Reconstruction
***************************

.. include:: calculate_pairwise_distances.rst
.. include:: calculate_neigbourjoining_tree.rst
.. include:: calculate_UPGMA_cluster.rst
.. include:: phylo_by_ls.rst
.. include:: maketree_from_proteinseqs.rst

************************
Python Coding Guidelines
************************

.. include:: coding_guidelines.rst
