.. jupyter-execute::
    :hide-code:

    import set_working_directory


.. _nonstationary-model-aa-inference:


Infering Amino Acid frequencies with a nonstationary model
==========================================================

.. sectionauthor:: Peter Goodman, Andrew Wheeler

This example demonstrates how to apply a non-stationary, branch-heterogeneous maximum likelihood model to infer state frequency vectors across different branches or branch groups on a phylogenetic tree, highlighting a practical use case. The model used is a maximum likelihood, time-heterogeneous, and time-reversible model.

Non-stationarity means that amino acid frequencies are not fixed over time. This characteristic allows a non-stationary model to infer new frequencies along different branches or clades of a tree.

In this example, we infer an amino acid frequency vector for the primate and rodent clades in a mammal phylogeny, as well as a frequency vector for the root of these two clades. The purpose of this inference was to produce three amino acid frequency vectors for comparison. The rodent vector represented the amino acid frequencies associated with clade reported to exhibit a higher effective population size (:math:`N_{e}`), while the primate vector showed frequencies associated with clade reported to exhibit a lower :math:`N_{e}`. These clade specific measures allow us to measure associations with other known properties of the clades, such as amino acid fitness metrics.

Load the necessary functions from Cogent3.

.. jupyter-execute::

    from cogent3 import load_tree, load_aligned_seqs
    from cogent3.evolve.substitution_model import EmpiricalProteinMatrix
    from cogent3.parse.paml_matrix import PamlMatrixParser

Load in the fasta formatted alignment file and designate molecule type, in this example we use protein.

.. jupyter-execute::

    aln = load_aligned_seqs("data/primate_rodent.fasta", moltype="protein")

Load in a substitution model (or select from models available in Cogent3). We use a single substitution model to start, and later apply branch heterogeneous frequencies.

.. jupyter-execute::

    with open("data/Q.mammal") as matrix_file:
        empirical_matrix , empirical_frequencies = PamlMatrixParser(matrix_file)
        sm = EmpiricalProteinMatrix(empirical_matrix, empirical_frequencies, with_rate=True, distribution="free")

Load in a phylogenetic tree file in Newick format.

.. jupyter-execute::

    tree = load_tree("data/primate_rodent.tre")

Designate clades based in your tree for which you want to infer different frequency vectors. Use two ingroup species and the closest possible outgroup species to define each clade. All internal branches and nodes are included in the designated clade. ``stem=True`` sets the root of the clade as part of the clade.

.. jupyter-execute::

    primate_edges = tree.get_edge_names("CHLOR_SAB","GALEO_VAR", outgroup_name="OCHOT_PRI", clade=True, stem=True)
    rodent_edges = tree.get_edge_names("MICRO_OCH","OCHOT_PRI", outgroup_name="GALEO_VAR", clade=True, stem=True)


We display the phylogenetic tree, using edge coloring to visualise the tree scopes represented by the model.

.. jupyter-execute::

    dnd = tree.get_figure()
    dnd.style_edges(primate_edges, line=dict(color="red"), legendgroup="Primates")
    dnd.style_edges(rodent_edges, line=dict(color="blue"), legendgroup="Rodents")
    dnd.scale_bar = None
    dnd.show(width=600, height=700)

Create a likelihood function.

.. jupyter-execute::

    lf = sm.make_likelihood_function(tree)
    lf.set_alignment(aln)

Tell the model to infer mprobs (frequency vectors) for each set of edges you defined above. Setting ``clade = True`` will give a single amino acid frequency vector for all the edges within the clade, rather than a new vector for every branch. Single edges, such as the root, can be designated to be optimized separately. Any edge not included in an optimization function will retain the amino acid frequencies of the substitution model you used.

.. jupyter-execute::

    lf.set_param_rule("mprobs", edges=primate_edges, clade=True, stem=True)
    lf.set_param_rule("mprobs", edges=rodent_edges, clade=True, stem=True)
    lf.set_param_rule("mprobs", edge="root")

Optimize the likelihood function.

.. jupyter-execute::

    lf.optimise(max_restarts=5, tolerance=1e-9, show_progress=False)
    lf
