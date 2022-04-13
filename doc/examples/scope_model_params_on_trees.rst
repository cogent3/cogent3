.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _scope-params-on-trees:

Allowing substitution model parameters to differ between branches
=================================================================

.. sectionauthor:: Gavin Huttley

A common task concerns assessing how substitution model exchangeability parameters differ between evolutionary lineages. This is most commonly of interest for the case of testing for natural selection. Here I'll demonstrate the different ways of scoping parameters across trees for the codon model case and how these can be used for evolutionary modelling.

We start with the standard imports, plus using a canned codon substitution model and then load the sample data set.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_tree
    from cogent3.evolve.models import MG94HKY

    aln = load_aligned_seqs("data/long_testseqs.fasta")
    tree = load_tree("data/test.tree")

We construct the substitution model and likelihood function and set the alignment.

.. jupyter-execute::

    sm = MG94HKY()
    lf = sm.make_likelihood_function(tree, digits=2, space=3)
    lf.set_alignment(aln)

At this point we have a likelihood function with two exchangeability parameters from the substitution model (``kappa`` the transition/transversion ratio; ``omega`` the nonsynonymous/synonymous ratio) plus branch lengths for all tree edges. To facilitate subsequent discussion I now display the tree

.. jupyter-execute::

    print(tree.ascii_art())

In order to scope a parameter on a tree (meaning specifying a subset of edges for which the parameter is to be treated differently to the remainder of the tree) requires uniquely identifying the edges. We do this using the following arguments to the likelihood function ``set_param_rule`` method:

- ``tip_names``: the name of two tips
- ``outgroup_name``: the name of a tip that is not part of the clade of interest
- ``clade``: if ``True``, all lineages descended from the tree node identified by the ``tip_names`` and ``outgroup_name`` argument are affected by the other arguments. If ``False``, then the ``stem`` argument must apply.
- ``stem``: Whether the edge leading to the node is included.

The next concepts include exactly what can be scoped and how. In the case of testing for distinctive periods of natural selection it is common to specify distinct values for ``omega`` for an edge. I'll first illustrate some possible uses for the arguments above by setting ``omega`` to be distinctive for specific edges. I will set a value for omega so that printing the likelihood function illustrates what edges have been effected, but I won't actually do any model fitting.

Specifying a clade
------------------

I'm going to cause ``omega`` to attain a different value for all branches aside from the primate clade and stem (``HowlerMon``, ``Human``, ``edge.0``).

.. jupyter-execute::

    lf.set_param_rule(
        "omega",
        tip_names=["DogFaced", "Mouse"],
        outgroup_name="Human",
        init=2.0,
        clade=True,
    )
    lf

As you can see ``omega`` for the primate edges I listed above have the default parameter value (1.0), while the others have what I've assigned. In fact, you could omit the ``clade`` argument as this is the default, but I think for readability of scripts it's best to be explicit.

Specifying a stem
-----------------

This time I'll specify the stem leading to the primates as the edge of interest.

.. note:: I need to reset the ``lf`` so all edges have the default value again. I'll show this only for this example, but rest assured I'm doing it for all others too.

.. jupyter-execute::

    lf.set_param_rule("omega", init=1.0)
    lf.set_param_rule(
        "omega",
        tip_names=["Human", "HowlerMon"],
        outgroup_name="Mouse",
        init=2.0,
        stem=True,
        clade=False,
    )
    lf

Specifying clade and stem
-------------------------

I'll specify that both the primates and their stem are to be considered.

.. jupyter-execute::
    :hide-code:

    lf.set_param_rule("omega", init=1.0)

.. jupyter-execute::

    lf.set_param_rule(
        "omega",
        tip_names=["Human", "HowlerMon"],
        outgroup_name="Mouse",
        init=2.0,
        stem=True,
        clade=True,
    )
    lf

Alternate arguments for specifying edges
----------------------------------------

The likelihood function ``set_param_rule`` method also has the arguments of ``edge`` and ``edges``. These allow specific naming of the tree edge(s) to be affected by a rule. In general, however, the ``tip_names`` + ``outgroup_name`` combo is more robust.

Applications of scoped parameters
---------------------------------

The general use-cases for which a tree scope can be applied are:

1. constraining all edges identified by a rule to have a specific value which is constant and not modifiable

.. code-block:: python

    lf.set_param_rule(
        "omega",
        tip_names=["Human", "HowlerMon"],
        outgroup_name="Mouse",
        clade=True,
        is_constant=True,
    )

2. all edges identified by a rule have the same but different value to the rest of the tree

.. code-block:: python

    lf.set_param_rule(
        "omega", tip_names=["Human", "HowlerMon"], outgroup_name="Mouse", clade=True
    )

3. allowing all edges identified by a rule to have different values of the parameter with the remaining tree edges having the same value

.. code-block:: python

    lf.set_param_rule(
        "omega",
        tip_names=["Human", "HowlerMon"],
        outgroup_name="Mouse",
        clade=True,
        is_independent=True,
    )

4. allowing all edges to have a different value

.. code-block:: python

    lf.set_param_rule("omega", is_independent=True)

I'll demonstrate these cases sequentially as they involve gradually increasing the degrees of freedom in the model. First we'll constrain ``omega`` to equal 1 on the primate edges. I'll then optimise the model.

.. note:: here I'm specifying a constant value for the parameter and so I **must** use the argument ``value`` to set it. This not to be confused with the argument ``init`` that is used for providing initial (starting) values for fitting.

.. jupyter-execute::
    :hide-code:

    lf.set_param_rule("omega", init=1.0)

.. jupyter-execute::

    lf.set_param_rule(
        "omega",
        tip_names=["Human", "HowlerMon"],
        outgroup_name="Mouse",
        clade=True,
        value=1.0,
        is_constant=True,
    )
    lf.optimise(local=True, show_progress=False)
    lf

I'll now free up ``omega`` on the primate clade, but making it a single value shared by all primate lineages.

.. jupyter-execute::

    lf.set_param_rule(
        "omega",
        tip_names=["Human", "HowlerMon"],
        outgroup_name="Mouse",
        clade=True,
        is_constant=False,
    )
    lf.optimise(local=True, show_progress=False)
    lf

Finally I'll allow all primate edges to have different values of ``omega``.

.. jupyter-execute::

    lf.set_param_rule(
        "omega",
        tip_names=["Human", "HowlerMon"],
        outgroup_name="Mouse",
        clade=True,
        is_independent=True,
    )
    lf.optimise(local=True, show_progress=False)
    lf

We now allow ``omega`` to be different on all edges.

.. jupyter-execute::

    lf.set_param_rule("omega", is_independent=True)
    lf.optimise(local=True, show_progress=False)
    lf
