.. jupyter-execute::
    :hide-code:

    import set_working_directory

**************************************
Evolutionary Analysis Using Likelihood
**************************************

Specifying substitution models
==============================

The available pre-defined substitution models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In cases where code takes a substitution model as an argument, you can use the value under “Abbreviation” as a string.

.. jupyter-execute::

    from cogent3 import available_models

    available_models()

Getting a substitution model with ``get_model()``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    from cogent3.evolve.models import get_model

    hky = get_model("HKY85")
    print(hky)

Rate heterogeneity models
~~~~~~~~~~~~~~~~~~~~~~~~~

We illustrate this for the gamma distributed case using examples of the canned models displayed above. Creating rate heterogeneity variants of the canned models can be done by using optional arguments that get passed to the substitution model class.

For nucleotide
--------------

We specify a general time reversible nucleotide model with gamma distributed rate heterogeneity.

.. jupyter-execute::

    from cogent3.evolve.models import get_model

    sub_mod = get_model("GTR", with_rate=True, distribution="gamma")
    print(sub_mod)

For codon
---------

We specify a conditional nucleotide frequency codon model with nucleotide general time reversible parameters and a parameter for the ratio of nonsynonymous to synonymous substitutions (omega) with gamma distributed rate heterogeneity.

.. jupyter-execute::

    from cogent3.evolve.models import get_model

    sub_mod = get_model("CNFGTR", with_rate=True, distribution="gamma")
    print(sub_mod)

For protein
-----------

We specify a Jones, Taylor and Thornton 1992 empirical protein substitution model with gamma distributed rate heterogeneity.

.. jupyter-execute::

    from cogent3.evolve.models import get_model

    sub_mod = get_model("JTT92", with_rate=True, distribution="gamma")
    print(sub_mod)

Making a likelihood function
============================

You start by specifying a substitution model and use that to construct a likelihood function for a specific tree.

.. jupyter-execute::

    from cogent3 import make_tree
    from cogent3.evolve.models import get_model

    sub_mod = get_model("F81")
    tree = make_tree("(a,b,(c,d))")
    lf = sub_mod.make_likelihood_function(tree)

Providing an alignment to a likelihood function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You need to load an alignment and then provide it a likelihood function. I construct very simple trees and alignments for this example.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs, make_tree
    from cogent3.evolve.models import get_model

    sub_mod = get_model("F81")
    tree = make_tree("(a,b,(c,d))")
    lf = sub_mod.make_likelihood_function(tree)
    aln = make_aligned_seqs(
        [("a", "ACGT"), ("b", "AC-T"), ("c", "ACGT"), ("d", "AC-T")]
    )
    lf.set_alignment(aln)

Scoping parameters on trees – time heterogeneous models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For many evolutionary analyses, it’s desirable to allow different branches on a tree to have different values of a parameter. We show this for a simple codon model case here where we want the great apes (the clade that includes human and orangutan) to have a different value of the ratio of nonsynonymous to synonymous substitutions. This parameter is identified in the precanned ``CNFGTR`` model as ``omega``.

.. jupyter-execute::

    from cogent3 import load_tree
    from cogent3.evolve.models import get_model

    tree = load_tree("data/primate_brca1.tree")
    print(tree.ascii_art())

.. jupyter-execute::

    sm = get_model("CNFGTR")
    lf = sm.make_likelihood_function(tree, digits=2)
    lf.set_param_rule(
        "omega",
        tip_names=["Human", "Orangutan"],
        outgroup_name="Galago",
        clade=True,
        init=0.5,
    )

We’ve set an *initial* value for this clade so that the edges affected by this rule are evident below.

.. jupyter-execute::

    lf

A more extensive description of capabilities is in :ref:`scope-params-on-trees`.

Specifying a parameter as constant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This means the parameter will not be modified during likelihood maximisation. We show this here by making the ``omega`` parameter constant at the value 1 – essentially the condition of selective neutrality.

.. jupyter-execute::

    from cogent3 import load_tree
    from cogent3.evolve.models import get_model

    tree = load_tree("data/primate_brca1.tree")
    sm = get_model("CNFGTR")
    lf = sm.make_likelihood_function(tree, digits=2)
    lf.set_param_rule("omega", is_constant=True)

Providing a starting value for a parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This can be useful to improve performance, the closer you are to the maximum likelihood estimator the quicker optimisation will be.

.. jupyter-execute::

    from cogent3 import load_tree
    from cogent3.evolve.models import get_model

    tree = load_tree("data/primate_brca1.tree")
    sm = get_model("CNFGTR")
    lf = sm.make_likelihood_function(tree, digits=2)
    lf.set_param_rule("omega", init=0.1)

Setting parameter bounds for optimisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This can be useful for stopping optimisers from getting stuck in a bad part of parameter space. The following is for ``omega`` in a codon model. I’m also providing an initial guess for the parameter (``init=0.1``) as well as a lower bound. An initial guess that is close to the maximum likelihood estimate will speed up optimisation.

.. jupyter-execute::

    from cogent3 import load_tree
    from cogent3.evolve.models import get_model

    tree = load_tree("data/primate_brca1.tree")
    sm = get_model("CNFGTR")
    lf = sm.make_likelihood_function(tree, digits=2)
    lf.set_param_rule("omega", init=0.1, lower=1e-9, upper=20.0)

Setting an upper bound for branch length
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the branch length estimates seem too large, setting just an upper bound can be sensible. This will apply to all edges on the tree.

.. jupyter-execute::

    from cogent3 import load_tree
    from cogent3.evolve.models import get_model

    tree = load_tree("data/primate_brca1.tree")
    sm = get_model("F81")
    lf = sm.make_likelihood_function(tree)
    lf.set_param_rule("length", upper=1.0)

.. note:: If, after optimising, the branch lengths equal to the upper value you set then the function has not been fully maximised and you should consider adjusting the boundary again.

Specifying rate heterogeneity functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We extend the simple gamma distributed rate heterogeneity case for nucleotides from above to construction of the actual likelihood function. We do this for 4 bins and constraint the bin probabilities to be equal.

.. jupyter-execute::

    from cogent3 import load_tree
    from cogent3.evolve.models import get_model

    sm = get_model("GTR", with_rate=True, distribution="gamma")
    tree = load_tree("data/primate_brca1.tree")
    lf = sm.make_likelihood_function(tree, bins=4, digits=2)
    lf.set_param_rule("bprobs", is_constant=True)

For more detailed discussion of defining and using these models see :ref:`rate-heterogeneity`.

Specifying Phylo-HMMs
~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    from cogent3 import load_tree
    from cogent3.evolve.models import get_model

    sm = get_model("GTR", with_rate=True, distribution="gamma")
    tree = load_tree("data/primate_brca1.tree")
    lf = sm.make_likelihood_function(tree, bins=4, sites_independent=False, digits=2)
    lf.set_param_rule("bprobs", is_constant=True)

For more detailed discussion of defining and using these models see :ref:`rate-heterogeneity-hmm`.

Fitting likelihood functions - Choice of optimisers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are 2 types of optimiser: simulated annealing, a *global* optimiser; and Powell, a *local* optimiser. The simulated annealing method is slow compared to Powell and in general Powell is an adequate choice. I setup a simple nucleotide model to illustrate these.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree
    from cogent3.evolve.models import get_model

    tree = load_tree("data/primate_brca1.tree")
    aln = load_aligned_seqs("data/primate_brca1.fasta")
    sm = get_model("F81")
    lf = sm.make_likelihood_function(tree, digits=3, space=2)
    lf.set_alignment(aln)
    lf.optimise(show_progress=False)

The default is to use Powell. For Powell, it’s recommended to set the ``max_restarts`` argument since this provides a mechanism for Powell to attempt restarting the optimisation from a slightly different spot which can help in overcoming local maxima.

.. jupyter-execute::

    lf.optimise(local=True, max_restarts=5, show_progress=False)

We might want to do crude simulated annealing following by more rigorous Powell. To do this we first need to use the global optimiser, setting ``local=False`` setting a large value for ``global_tolerance``.

.. jupyter-execute::

    lf.optimise(local=False, global_tolerance=1.0, show_progress=False)

Followed by a standard call to ``optimise()``.

.. jupyter-execute::

    lf.optimise(show_progress=False, max_restarts=5, tolerance=1e-8)

How to check your optimisation was successful
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is no guarantee that an optimised function has achieved a global maximum. We can, however, be sure that a maximum was achieved by validating that the optimiser stopped because the specified tolerance condition was met, rather than exceeding the maximum number of evaluations. The latter number is set to ensure optimisation doesn’t proceed endlessly. If the optimiser exited because this limit was exceeded you can be sure that the function **has not** been successfully optimised.

We can monitor this situation using the ``limit_action`` argument to ``optimise``. Providing the value ``raise`` causes an exception to be raised if this condition occurs, as shown below. Providing ``warn`` (default) instead will cause a warning message to be printed to screen but execution will continue. The value ``ignore`` hides any such message.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree
    from cogent3.evolve.models import get_model

    tree = load_tree("data/primate_brca1.tree")
    aln = load_aligned_seqs("data/primate_brca1.fasta")
    sm = get_model("F81")
    lf = sm.make_likelihood_function(tree, digits=3, space=2)
    lf.set_alignment(aln)
    try:
        lf.optimise(
            show_progress=False,
            limit_action="raise",
            max_evaluations=10,
            return_calculator=True,
        )
    except ArithmeticError as err:
        print(err)

.. note:: We recommend using ``limit_action='raise'`` and catching the ``ArithmeticError`` error explicitly (as demonstrated above). You really shouldn't be using results from such an optimisation run.

Overview of the fitted likelihood function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Jupyter, the likelihood function object presents a representation of the main object features.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree
    from cogent3.evolve.models import get_model

    sm = get_model("GTR")
    tree = load_tree("data/primate_brca1.tree")
    lf = sm.make_likelihood_function(tree)
    aln = load_aligned_seqs("data/primate_brca1.fasta")
    lf.set_alignment(aln)
    lf.optimise(local=True, show_progress=False)
    lf

Log likelihood and number of free parameters
--------------------------------------------

Reusing the optimised ``lf`` object from above, we can get the log-likelihood and the number of free parameters.

.. jupyter-execute::

    lnL = lf.lnL
    lnL

.. jupyter-execute::

    nfp = lf.nfp
    nfp

.. warning:: The number of free parameters (nfp) refers only to the number of parameters that were modifiable by the optimiser. Typically, the degrees-of-freedom of a likelihood ratio test statistic is computed as the difference in nfp between models. This will not be correct for models in which a boundary conditions exist (rate heterogeneity models where a parameter value boundary is set between bins).

Aikake Information Criterion
----------------------------

Reusing the optimised ``lf`` object from above.

.. jupyter-execute::

    lf.get_aic()

We can also get the second-order AIC.

.. jupyter-execute::

    lf.get_aic(second_order=True)

Bayesian Information Criterion
------------------------------

Reusing the optimised ``lf`` object from above.

.. jupyter-execute::

    lf.get_bic()

Getting maximum likelihood estimates
------------------------------------

Reusing the optimised ``lf`` object from above.

One at a time
'''''''''''''

We get the statistics out individually. We get the ``length`` for the Human edge and the exchangeability parameter ``A/G``.

.. jupyter-execute::

    a_g = lf.get_param_value("A/G")
    a_g

.. jupyter-execute::

    human = lf.get_param_value("length", "Human")
    human

Just the motif probabilities
''''''''''''''''''''''''''''

.. jupyter-execute::

    mprobs = lf.get_motif_probs()
    mprobs

As tables
'''''''''

.. jupyter-execute::

    tables = lf.get_statistics(with_motif_probs=True, with_titles=True)
    tables[0]  # just displaying the first

Testing Hypotheses - Using Likelihood Ratio Tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We test the molecular clock hypothesis for human and chimpanzee lineages. The null has these two branches constrained to be equal.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree
    from cogent3.evolve.models import get_model

    tree = load_tree("data/primate_brca1.tree")
    aln = load_aligned_seqs("data/primate_brca1.fasta")
    sm = get_model("F81")
    lf = sm.make_likelihood_function(tree, digits=3, space=2)
    lf.set_alignment(aln)
    lf.set_param_rule(
        "length",
        tip_names=["Human", "Chimpanzee"],
        outgroup_name="Galago",
        clade=True,
        is_independent=False,
    )
    lf.set_name("Null Hypothesis")
    lf.optimise(local=True, show_progress=False)
    null_lnL = lf.lnL
    null_nfp = lf.nfp
    lf

The alternate allows the human and chimpanzee branches to differ by just setting all lengths to be independent.

.. jupyter-execute::

    lf.set_param_rule("length", is_independent=True)
    lf.set_name("Alt Hypothesis")
    lf.optimise(local=True, show_progress=False)
    alt_lnL = lf.lnL
    alt_nfp = lf.nfp
    lf

We import the function for computing the probability of a chi-square test statistic, compute the likelihood ratio test statistic, degrees of freedom and the corresponding probability.

.. jupyter-execute::

    from scipy.stats.distributions import chi2

    LR = 2 * (alt_lnL - null_lnL)  # the likelihood ratio statistic
    df = alt_nfp - null_nfp  # the test degrees of freedom
    p = chi2.sf(LR, df)
    print(f"LR={LR:.4f} ; df={df}; p={df:.4f}")

Testing Hypotheses - By parametric bootstrapping
------------------------------------------------

If we can't rely on the asymptotic behaviour of the LRT, e.g. due to small alignment length, we can use a parametric bootstrap. Convenience functions for that are described in more detail here :ref:`parametric-bootstrap`.

In general, however, this capability derives from the ability of any defined ``evolve`` likelihood function to simulate an alignment. This property is provided as ``simulate_alignment`` method on likelihood function objects.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree
    from cogent3.evolve.models import get_model

    tree = load_tree("data/primate_brca1.tree")
    aln = load_aligned_seqs("data/primate_brca1.fasta")

    sm = get_model("F81")
    lf = sm.make_likelihood_function(tree, digits=3, space=2)
    lf.set_alignment(aln)
    lf.set_param_rule(
        "length",
        tip_names=["Human", "Chimpanzee"],
        outgroup_name="Galago",
        clade=True,
        is_independent=False,
    )
    lf.set_name("Null Hypothesis")
    lf.optimise(local=True, show_progress=False)
    sim_aln = lf.simulate_alignment()
    sim_aln[:60]

Determining confidence intervals on MLEs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The profile method is used to calculate a confidence interval for a named parameter. We show it here for a global substitution model exchangeability parameter (*kappa*, the ratio of transition to transversion rates) and for an edge specific parameter (just the human branch length).

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree
    from cogent3.evolve.models import get_model

    tree = load_tree("data/primate_brca1.tree")
    aln = load_aligned_seqs("data/primate_brca1.fasta")
    sm = get_model("HKY85")
    lf = sm.make_likelihood_function(tree)
    lf.set_alignment(aln)
    lf.optimise(local=True, show_progress=False)
    kappa_lo, kappa_mle, kappa_hi = lf.get_param_interval("kappa")
    print(f"lo={kappa_lo:.2f} ; mle={kappa_mle:.2f} ; hi={kappa_hi:.2f}")
    human_lo, human_mle, human_hi = lf.get_param_interval("length", "Human")
    print(f"lo={human_lo:.2f} ; mle={human_mle:.2f} ; hi={human_hi:.2f}")

Saving results
~~~~~~~~~~~~~~

The best approach is to use the json string from the ``to_json()`` method. The saved data can be later reloaded using ``cogent3.util.deserialise.deserialise_object()``. The ``json`` data contains the alignment, tree topology, substitution model, parameter values, etc..

To illustrate this, I create a very simple likelihood function. The ``json`` variable below is just a string that can be saved to disk.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree
    from cogent3.evolve.models import get_model

    aln = make_aligned_seqs(data=dict(a="ACGG", b="ATAG", c="ATGG"))
    tree = make_tree(tip_names=aln.names)
    sm = get_model("F81")
    lf = sm.make_likelihood_function(tree)
    lf.set_alignment(aln)
    json = lf.to_json()
    json[:60]  # just truncating the displayed string

We deserialise the object from the string.

.. jupyter-execute::

    from cogent3.util.deserialise import deserialise_object

    newlf = deserialise_object(json)
    newlf

Reconstructing ancestral sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We first fit a likelihood function.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree
    from cogent3.evolve.models import get_model

    tree = load_tree("data/primate_brca1.tree")
    aln = load_aligned_seqs("data/primate_brca1.fasta")
    sm = get_model("F81")
    lf = sm.make_likelihood_function(tree, digits=3, space=2)
    lf.set_alignment(aln)
    lf.optimise(show_progress=False)

We then get the most likely ancestral sequences.

.. jupyter-execute::

    ancestors = lf.likely_ancestral_seqs()
    ancestors[:60]

Or we can get the posterior probabilities (returned as a ``DictArray``) of sequence states at each node.

.. jupyter-execute::

    ancestral_probs = lf.reconstruct_ancestral_seqs()
    ancestral_probs["root"][:5]

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There’s nothing that improves performance quite like being close to the maximum likelihood values. So using the ``set_param_rule`` method to provide good starting values can be very useful. As this can be difficult to do one easy way is to build simpler models that are nested within the one you’re interested in. Fitting those models and then relaxing constraints until you’re at the parameterisation of interest can markedly improve optimisation speed.
