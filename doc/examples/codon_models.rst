.. jupyter-execute::
    :hide-code:

    import set_working_directory

Using codon models
==================

.. sectionauthor:: Gavin Huttley

The basic paradigm for evolutionary modelling is:

#. construct the codon substitution model
#. constructing likelihood function
#. modify likelihood function (setting rules)
#. providing the alignment(s)
#. optimisation
#. get results out

.. note:: In the following, a result followed by '...' just means the output has been truncated for the sake of a succinct presentation.

Constructing the codon substitution model
-----------------------------------------

For the time-reversible category, Cogent3 implements 4 basic rate matrices: i) NF models, these are nucleotide frequency weighted rate matrices and were initially described by Muse and Gaut (1994); ii) a variant of (i) where position specific nucleotide frequencies are used; iii) TF models, these are tuple (codon in this case) frequency weighted rate matrices and were initially described by Goldman and Yang (1994); iv) CNF, these use the conditional nucleotide frequency and have developed by Yap, Lindsay, Easteal and Huttley. These different models can be created using provided convenience functions which will be the case here, or specified by directly calling the ``TimeReversibleCodon`` substitution model class and setting the argument ``mprob_model`` equal to:

- NF, ``mprob_model='monomer'``
- NF with position specific nucleotide frequencies, ``mprob_model='monomers'``
- TF, ``mprob_model=None``
- CNF, ``mprob_model='conditional'``

In the following I will construct GTR variants of i and iv and a HKY variant of iii.

We import these explicitly from the ``cogent3.evolve.models`` module.

.. jupyter-execute::

    from cogent3.evolve.models import get_model

These are functions and calling them returns the indicated substitution model with default behaviour of recoding gap characters into N's.

.. jupyter-execute::

    tf = get_model("GY94")
    nf = get_model("MG94GTR")
    cnf = get_model("CNFGTR")

In the following demonstration I will use only the CNF form (``cnf``).

For our example we load a sample alignment and tree as per usual. To reduce the computational overhead for this example we will limit the number of sampled taxa.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree

    aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
    tree = load_tree("data/primate_brca1.tree")

Standard test of neutrality
---------------------------

We construct a likelihood function and constrain omega parameter (the ratio of nonsynonymous to synonymous substitutions) to equal 1. We also set some display formatting parameters.

.. jupyter-execute::

    lf = cnf.make_likelihood_function(tree, digits=2, space=3)
    lf.set_param_rule("omega", is_constant=True, value=1.0)

We then provide an alignment and optimise the model. In the current case we just use the local optimiser (hiding progress to keep this document succinct). We then print(the results.)

.. note:: I'm going to specify a set of conditions that will be used for all optimiser steps. For those new to python, one can construct a dictionary with the following form ``{'argument_name': argument_value}``, or alternatively ``dict(argument_name=argument_value)``. I'm doing the latter. This dictionary is then passed to functions/methods by prefacing it with ``**``.

.. jupyter-execute::

    optimiser_args = dict(
        local=True, max_restarts=5, tolerance=1e-8, show_progress=False
    )
    lf.set_alignment(aln)
    lf.optimise(**optimiser_args)
    lf

In the above output, the first table shows the maximum likelihood estimates (MLEs) for the substitution model parameters that are 'global' in scope. For instance, the ``C/T=4.58`` MLE indicates that the relative rate of substitutions between C and T is nearly 5 times the background substitution rate.

The above function has been fit using the default counting procedure for estimating the motif frequencies, i.e. codon frequencies are estimated as the average of the observed codon frequencies. If you wanted to numerically optimise the motif probabilities, then modify the likelihood function creation line to

.. code-block:: python

    lf = cnf.make_likelihood_function(tree, optimise_motif_probs=True)

We can then free up the omega parameter, but before we do that we'll store the log-likelihood and number of free parameters for the current model form for reuse later.

.. jupyter-execute::

    neutral_lnL = lf.get_log_likelihood()
    neutral_nfp = lf.get_num_free_params()
    lf.set_param_rule("omega", is_constant=False)
    lf.optimise(**optimiser_args)
    non_neutral_lnL = lf.get_log_likelihood()
    non_neutral_nfp = lf.get_num_free_params()
    lf

We then conduct a likelihood ratio test whether the MLE of omega significantly improves the fit over the constraint it equals 1. We import the convenience function from the ``cogent3`` stats module.

.. jupyter-execute::

    from scipy.stats.distributions import chi2

    LR = 2 * (non_neutral_lnL - neutral_lnL)
    df = non_neutral_nfp - neutral_nfp
    print(chi2.sf(LR, df))

Not surprisingly, this is significant. We then ask whether the Human and Chimpanzee edges have a value of omega that is significantly different from the rest of the tree.

.. jupyter-execute::

    lf.set_param_rule(
        "omega", tip_names=["Chimpanzee", "Human"], outgroup_name="Galago", clade=True
    )
    lf.optimise(**optimiser_args)
    lf
    chimp_human_clade_lnL = lf.get_log_likelihood()
    chimp_human_clade_nfp = lf.get_num_free_params()

.. jupyter-execute::

    LR = 2 * (chimp_human_clade_lnL - non_neutral_lnL)
    df = chimp_human_clade_nfp - non_neutral_nfp
    print(chi2.sf(LR, df))

This is basically a replication of the original Huttley et al (2000) result for *BRCA1*.

Rate-heterogeneity model variants
---------------------------------

It is also possible to specify rate-heterogeneity variants of these models. In the first instance we'll create a likelihood function where these rate-classes are global across the entire tree. Because fitting these models can be time consuming I'm going to recreate the non-neutral likelihood function from above first, fit it, and then construct the rate-heterogeneity likelihood function. By doing this I can ensure that the richer model starts with parameter values that produce a log-likelihood the same as the null model, ensuring the subsequent optimisation step improves the likelihood over the null.

.. jupyter-execute::

    lf = cnf.make_likelihood_function(tree, digits=2, space=3)
    lf.set_alignment(aln)
    lf.optimise(**optimiser_args)
    non_neutral_lnL = lf.get_log_likelihood()
    non_neutral_nfp = lf.get_num_free_params()

Now, we have a null model which we know (from having fit it above) has a MLE < 1. We will construct a rate-heterogeneity model with just 2 rate-classes (neutral and adaptive) that are separated by the boundary of omega=1. These rate-classes are specified as discrete bins in Cogent3 and the model configuration steps for a bin or bins are done using the ``set_param_rule`` method. To ensure the alternate model starts with a likelihood at least as good as the previous we need to make the probability of the neutral site-class bin ~= 1 (these are referenced by the ``bprobs`` parameter type) and assign the null model omega MLE to this class.

To get all the parameter MLEs (branch lengths, GTR terms, etc ..) into the alternate model we get an annotated tree from the null model which will have these values associated with it.

.. jupyter-execute::

    annot_tree = lf.get_annotated_tree()
    omega_mle = lf.get_param_value("omega")

We can then construct a new likelihood function, specifying the rate-class properties.

.. jupyter-execute::

    rate_lf = cnf.make_likelihood_function(
        annot_tree, bins=["neutral", "adaptive"], digits=2, space=3
    )

We define a very small value (``epsilon``) that is used to specify the starting values.

.. jupyter-execute::

    epsilon = 1e-6

We now provide starting parameter values for ``omega`` for the two bins, setting the boundary

.. jupyter-execute::

    rate_lf.set_param_rule("omega", bin="neutral", upper=1, init=omega_mle)
    rate_lf.set_param_rule(
        "omega", bin="adaptive", lower=1 + epsilon, upper=100, init=1 + 2 * epsilon
    )

and provide the starting values for the bin probabilities (``bprobs``).

.. jupyter-execute::

    rate_lf.set_param_rule("bprobs", init=[1 - epsilon, epsilon])

The above statement essentially assigns a probability of nearly 1 to the 'neutral' bin. We now set the alignment and fit the model.

.. jupyter-execute::

    rate_lf.set_alignment(aln)
    rate_lf.optimise(**optimiser_args)
    rate_lnL = rate_lf.get_log_likelihood()
    rate_nfp = rate_lf.get_num_free_params()
    LR = 2 * (rate_lnL - non_neutral_lnL)
    df = rate_nfp - non_neutral_nfp
    rate_lf

.. jupyter-execute::

    print(chi2.sf(LR, df))

We can get the posterior probabilities of site-classifications out of this model as

.. jupyter-execute::

    pp = rate_lf.get_bin_probs()

This is a ``DictArray`` class which stores the probabilities as a ``numpy.array``.

Mixing branch and site-heterogeneity
------------------------------------

The following implements a modification of the approach of Zhang, Nielsen and Yang (Mol Biol Evol, 22:2472–9, 2005). For this model class, there are groups of branches for which all positions are evolving neutrally but some proportion of those neutrally evolving sites change to adaptively evolving on so-called foreground edges. For the current example, we'll define the Chimpanzee and Human branches as foreground and everything else as background. The following table defines the parameter scopes.

.. jupyter-execute::
    :hide-code:

    from IPython.core.display import HTML
    from numpy import array

    from cogent3 import make_table

    header = ['Site Class', 'Proportion', 'Background Edges', 'Foreground Edges']
    data = {'Site Class': array(['0', '1', '2a', '2b'], dtype='<U2'), 'Proportion': array(['p0', 'p1', 'p2', 'p3'], dtype='<U2'), 'Background Edges': array(['0 < omega0 < 1', 'omega1 = 1', '0 < omega0 < 1', 'omega1 = 1'],
      dtype='<U14'), 'Foreground Edges': array(['0 < omega0 < 1', 'omega1 = 1', '0 < omega2 > 1', '0 < omega0 < 1'],
      dtype='<U14')}
    data = {k: array(data[k], dtype='U') for k in data}
    table = make_table(header, data=data)
    HTML(table.set_repr_policy(show_shape=False))

.. note:: Our implementation is not as parametrically succinct as that of Zhang et al, we have 1 additional bin probability.

After Zhang et al, we first define a null model that has 2 rate classes '0' and '1'. We also get all the MLEs out using ``get_statistics``, just printing out the bin parameters table in the current case.

.. jupyter-execute::

    rate_lf = cnf.make_likelihood_function(tree, bins=["0", "1"], digits=2, space=3)
    rate_lf.set_param_rule("omega", bin="0", upper=1.0 - epsilon, init=1 - epsilon)
    rate_lf.set_param_rule("omega", bins="1", is_constant=True, value=1.0)
    rate_lf.set_alignment(aln)
    rate_lf.optimise(**optimiser_args)
    tables = rate_lf.get_statistics(with_titles=True)
    for table in tables:
        if "bin" in table.title:
            print(table)

We're also going to use the MLEs from the ``rate_lf`` model, since that nests within the more complex branch by rate-class model. This is unfortunately quite ugly compared with just using the annotated tree approach described above. It is currently necessary, however, due to a bug in constructing annotated trees for models with binned parameters.

.. jupyter-execute::

    globals = [t for t in tables if "global" in t.title][0]
    globals = dict(zip(globals.header, globals.to_list()[0]))
    bin_params = [t for t in tables if "bin" in t.title][0]
    rate_class_omegas = dict(bin_params.to_list(["bin", "omega"]))
    rate_class_probs = dict(bin_params.to_list(["bin", "bprobs"]))
    lengths = [t for t in tables if "edge" in t.title][0]
    lengths = dict(lengths.to_list(["edge", "length"]))

We now create the more complex model,

.. jupyter-execute::

    rate_branch_lf = cnf.make_likelihood_function(
        tree, bins=["0", "1", "2a", "2b"], digits=2, space=3
    )

and set from the nested null model the branch lengths,

.. jupyter-execute::

    for branch, length in lengths.items():
        rate_branch_lf.set_param_rule("length", edge=branch, init=length)

GTR term MLES,

.. jupyter-execute::

    for param, mle in globals.items():
        rate_branch_lf.set_param_rule(param, init=mle)

binned parameter values,

.. jupyter-execute::

    rate_branch_lf.set_param_rule(
        "omega", bins=["0", "2a"], upper=1.0, init=rate_class_omegas["0"]
    )
    rate_branch_lf.set_param_rule(
        "omega", bins=["1", "2b"], is_constant=True, value=1.0
    )
    rate_branch_lf.set_param_rule(
        "omega",
        bins=["2a", "2b"],
        edges=["Chimpanzee", "Human"],
        init=99,
        lower=1.0,
        upper=100.0,
        is_constant=False,
    )

and the bin probabilities.

.. jupyter-execute::

    rate_branch_lf.set_param_rule(
        "bprobs",
        init=[
            rate_class_probs["0"] - epsilon,
            rate_class_probs["1"] - epsilon,
            epsilon,
            epsilon,
        ],
    )

The result of these steps is to create a rate/branch model with initial parameter values that result in likelihood the same as the null.

.. jupyter-execute::

    rate_branch_lf.set_alignment(aln)

.. code-block:: python

    rate_branch_lf.optimise(**optimiser_args)
    print(rate_branch_lf)
    Likelihood function statistics
    log-likelihood = -6753.4561
    number of free parameters = 21
    =========================
          edge   bin    omega
    -------------------------
        Galago     0     0.00
        Galago     1     1.00
        Galago    2a     0.00
        Galago    2b     1.00
     HowlerMon     0     0.00
     HowlerMon     1     1.00
     HowlerMon    2a     0.00
     HowlerMon    2b     1.00
        Rhesus     0     0.00
        Rhesus     1     1.00
        ...
