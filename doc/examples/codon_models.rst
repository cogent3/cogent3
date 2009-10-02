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

PyCogent implements 4 basic rate matrices, described in a recently accepted manuscript: i) NF models, these are nucleotide frequency weighted rate matrices and were initially described by Muse and Gaut (1994); ii) a variant of (i) where position specific nucleotide frequencies are used; iii) TF models, these are tuple (codon in this case) frequency weighted rate matrices and were initially described by Goldman and Yang (1994); iv) CNF, these use the conditional nucleotide frequency and have developed by Yap, Lindsay, Easteal and Huttley. These different models can be created using provided convenience functions which will be the case here, or specified by directly calling the ``Codon`` substitution model class and setting the argument ``mprob_model`` equal to:

- NF, ``mprob_model='monomer'``
- NF with position specific nucleotide frequencies, ``mprob_model='monomers'``
- TF, ``mprob_model=None``
- CNF, ``mprob_model='conditional'``

.. warning:: The TF form is the currently the default, but this will be changed in the near future.

In the following I will construct GTR variants of i and iv and a HKY variant of iii.

We import these explicitly from the ``cogent.evolve.models`` module.

.. doctest::

    >>> from cogent.evolve.models import CNFGTR, MG94GTR, GY94

These are functions and calling them returns the indicated substitution model with default behaviour of recoding gap characters into N's.

.. doctest::

    >>> tf = GY94()
    >>> nf = MG94GTR()
    >>> cnf = CNFGTR()

In the following demonstration I will use only the CNF form (``cnf``).

For our example we load a sample alignment and tree as per usual. To reduce the computational overhead for this example we will limit the number of sampled taxa.

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree, DNA
    >>> aln = LoadSeqs('data/primate_brca1.fasta', moltype=DNA)
    >>> tree = LoadTree('data/primate_brca1.tree')

Standard test of neutrality
---------------------------

We construct a likelihood function and constrain omega parameter (the ratio of nonsynonymous to synonymous substitutions) to equal 1. We also set some display formatting parameters.

.. doctest::

    >>> lf = cnf.makeLikelihoodFunction(tree, digits=2, space=3)
    >>> lf.setParamRule('omega', is_const=True, value=1.0)

We then provide an alignment and optimise the model. In the current case we just use the local optimiser (hiding progress to keep this document succinct). We then print the results.

.. note:: I'm going to specify a set of conditions that will be used for all optimiser steps. For those new to python, one can construct a dictionary with the following form ``{'argument_name': argument_value}``, or alternatively ``dict(argument_name=argument_value)``. I'm doing the latter. This dictionary is then passed to functions/methods by prefacing it with ``**``.

.. doctest::
    
    >>> optimiser_args = dict(local=True, max_restarts=5,
    ...                     show_progress=False, tolerance=1e-8)
    >>> lf.setAlignment(aln)
    >>> lf.optimise(**optimiser_args)
    >>> print lf
    Likelihood Function Table
    ========================================
     A/C    A/G    A/T    C/G    C/T   omega
    ----------------------------------------
    1.10   4.07   0.84   1.95   4.58    1.00
    ----------------------------------------
    ============================
          edge   parent   length
    ----------------------------
        Galago     root     0.53
     HowlerMon     root     0.14
        Rhesus   edge.3     0.07
     Orangutan   edge.2     0.02
       Gorilla   edge.1     0.01
         Human   edge.0     0.02
    Chimpanzee   edge.0     0.01
        edge.0   edge.1     0.00
        edge.1   edge.2     0.01
        edge.2   edge.3     0.04
        edge.3     root     0.02
    ----------------------------
    ==============
    motif   mprobs
    --------------
      CTT     0.01
      ACC     0.00...

In the above output, the first table shows the maximum likelihood estimates (MLEs) for the substitution model parameters that are 'global' in scope. For instance, the ``C/T=4.58`` MLE indicates that the relative rate of substitutions between C and T is nearly 5 times the background substitution rate.

The above function has been fit using the default counting procedure for estimating the motif frequencies, i.e. codon frequencies are estimated as the average of the observed codon frequencies. If you wanted to numerically optimise the motif probabilities, then modify the likelihood function creation line to

.. code-block:: python

    lf = cnf.makeLikelihoodFunction(tree,optimise_motif_probs=True)

We can then free up the omega parameter, but before we do that we'll store the log-likelihood and number of free parameters for the current model form for reuse later.

.. doctest::

    >>> neutral_lnL = lf.getLogLikelihood()
    >>> neutral_nfp = lf.getNumFreeParams()
    >>> lf.setParamRule('omega', is_const=False)
    >>> lf.optimise(**optimiser_args)
    >>> print lf
    Likelihood Function Table
    ========================================
     A/C    A/G    A/T    C/G    C/T   omega
    ----------------------------------------
    1.08   3.86   0.78   1.96   4.08    0.75
    ----------------------------------------
    ============================
          edge   parent   length
    ----------------------------
        Galago     root     0.53
     HowlerMon     root     0.14...
    >>> non_neutral_lnL = lf.getLogLikelihood()
    >>> non_neutral_nfp = lf.getNumFreeParams()

We then conduct a likelihood ratio test whether the MLE of omega significantly improves the fit over the constraint it equals 1. We import the convenience function from the cogent stats module.

    >>> from cogent.maths.stats import chisqprob
    >>> LR = 2*(non_neutral_lnL-neutral_lnL)
    >>> df = non_neutral_nfp - neutral_nfp
    >>> print chisqprob(LR, df)
    0.0026...

Not surprisingly, this is significant. We then ask whether the Human and Chimpanzee edges have a value of omega that is significantly different from the rest of the tree.

.. doctest::

    >>> lf.setParamRule('omega', tip_names=['Chimpanzee', 'Human'],
    ...                          outgroup_name='Galago', is_clade=True)
    >>> lf.optimise(**optimiser_args)
    >>> print lf
    Likelihood Function Table
    ================================
     A/C    A/G    A/T    C/G    C/T
    --------------------------------
    1.08   3.86   0.78   1.96   4.07
    --------------------------------
    ====================================
          edge   parent   length   omega
    ------------------------------------
        Galago     root     0.53    0.73
     HowlerMon     root     0.14    0.73
        Rhesus   edge.3     0.07    0.73
     Orangutan   edge.2     0.02    0.73
       Gorilla   edge.1     0.01    0.73
         Human   edge.0     0.02    2.39
    Chimpanzee   edge.0     0.01    2.39
        edge.0   edge.1     0.00    0.73...
    >>> chimp_human_clade_lnL = lf.getLogLikelihood()
    >>> chimp_human_clade_nfp = lf.getNumFreeParams()
    >>> LR = 2*(chimp_human_clade_lnL-non_neutral_lnL)
    >>> df = chimp_human_clade_nfp-non_neutral_nfp
    >>> print chisqprob(LR, df)
    0.028...

This is basically a replication of the original Huttley et al (2000) result for *BRCA1*.

Rate-heterogeneity model variants
---------------------------------

It is also possible to specify rate-heterogeneity variants of these models. In the first instance we'll create a likelihood function where these rate-classes are global across the entire tree. Because fitting these models can be time consuming I'm going to recreate the non-neutral likelihood function from above first, fit it, and then construct the rate-heterogeneity likelihood function. By doing this I can ensure that the richer model starts with parameter values that produce a log-likelihood the same as the null model, ensuring the subsequent optimisation step improves the likelihood over the null.

.. doctest::

    >>> lf = cnf.makeLikelihoodFunction(tree, digits=2, space=3)
    >>> lf.setAlignment(aln)
    >>> lf.optimise(**optimiser_args)
    >>> non_neutral_lnL = lf.getLogLikelihood()
    >>> non_neutral_nfp = lf.getNumFreeParams()

Now, we have a null model which we know (from having fit it above) has a MLE < 1. We will construct a rate-heterogeneity model with just 2 rate-classes (neutral and adaptive) that are separated by the boundary of omega=1. These rate-classes are specified as discrete bins in PyCogent and the model configuration steps for a bin or bins are done using the ``setParamRule`` method. To ensure the alternate model starts with a likelihood at least as good as the previous we need to make the probability of the neutral site-class bin ~= 1 (these are referenced by the ``bprobs`` parameter type) and assign the null model omega MLE to this class.

To get all the parameter MLEs (branch lengths, GTR terms, etc ..) into the alternate model we get an annotated tree from the null model which will have these values associated with it.

.. doctest::

    >>> annot_tree = lf.getAnnotatedTree()
    >>> omega_mle = lf.getParamValue('omega')

We can then construct a new likelihood function, specifying the rate-class properties.

.. doctest::

    >>> rate_lf = cnf.makeLikelihoodFunction(annot_tree,
    ...                     bins = ['neutral', 'adaptive'], digits=2, space=3)

We define a very small value (``epsilon``) that is used to specify the starting values.

.. doctest::

    >>> epsilon=1e-6

We now provide starting parameter values for ``omega`` for the two bins, setting the boundary

.. doctest::

    >>> rate_lf.setParamRule('omega', bin='neutral', upper=1, init=omega_mle)
    >>> rate_lf.setParamRule('omega', bin='adaptive', lower=1+epsilon,
    ...         upper=100, init=1+2*epsilon)

and provide the starting values for the bin probabilities (``bprobs``).

.. doctest::

    >>> rate_lf.setParamRule('bprobs', init=[1-epsilon, epsilon])

The above statement essentially assigns a probability of nearly 1 to the 'neutral' bin. We now set the alignment and fit the model.

.. doctest::

    >>> rate_lf.setAlignment(aln)
    >>> rate_lf.optimise(**optimiser_args)
    >>> rate_lnL = rate_lf.getLogLikelihood()
    >>> rate_nfp = rate_lf.getNumFreeParams()
    >>> LR = 2*(rate_lnL-non_neutral_lnL)
    >>> df = rate_nfp-non_neutral_nfp
    >>> print rate_lf
    Likelihood Function Table
    ============================
          edge   parent   length
    ----------------------------
        Galago     root     0.56
     HowlerMon     root     0.14
        Rhesus   edge.3     0.07
     Orangutan   edge.2     0.02
       Gorilla   edge.1     0.01
         Human   edge.0     0.02
    Chimpanzee   edge.0     0.01
        edge.0   edge.1     0.00
        edge.1   edge.2     0.01
        edge.2   edge.3     0.03
        edge.3     root     0.02
    ----------------------------
    =========================
         bin   bprobs   omega
    -------------------------
     neutral     0.14    0.01
    adaptive     0.86    1.17
    -------------------------
    ================================
     A/C    A/G    A/T    C/G    C/T
    --------------------------------
    1.07   3.96   0.78   1.96   4.20
    --------------------------------
    ==============
    motif   mprobs
    --------------
      CTT     0.01...
    >>> print chisqprob(LR, df)
    0.000...

We can get the posterior probabilities of site-classifications out of this model as

.. doctest::

    >>> pp = rate_lf.getBinProbs()

This is a ``DictArray`` class which stores the probabilities as a ``numpy.array``.

Mixing branch and site-heterogeneity
------------------------------------

The following implements a modification of the approach of Zhang, Nielsen and Yang (Mol Biol Evol, 22:2472–9, 2005). For this model class, there are groups of branches for which all positions are evolving neutrally but some proportion of those neutrally evolving sites change to adaptively evolving on so-called foreground edges. For the current example, we'll define the Chimpanzee and Human branches as foreground and everything else as background. The following table defines the parameter scopes.

+--------------+----------------+----------------------+---------------------+
|  Site class  |   Proportion   |   Background edges   |  Foreground edges   |
+==============+================+======================+=====================+
|           0  |          p_0   |      0 < omega_0 < 1 |   0 < omega_0 < 1   |
+--------------+----------------+----------------------+---------------------+
|           1  |          p_1   |          omega_1=1   |         omega_1=1   |
+--------------+----------------+----------------------+---------------------+
|          2a  |          p_2   |    0 < omega_0 < 1   |       omega_2 > 1   |
+--------------+----------------+----------------------+---------------------+
|          2b  |          p_3   |          omega_1=1   |       omega_2 > 1   |
+--------------+----------------+----------------------+---------------------+

.. note:: Our implementation is not as parametrically succinct as that of Zhang et al, we have 1 additional bin probability.

After Zhang et al, we first define a null model that has 2 rate classes '0' and '1'. We also get all the MLEs out using ``getStatistics``, just printing out the bin parameters table in the current case.

.. doctest::
    
    >>> rate_lf = cnf.makeLikelihoodFunction(tree, bins = ['0', '1'],
    ...                              digits=2, space=3)
    >>> rate_lf.setParamRule('omega', bin='0', upper=1.0-epsilon,
    ...                      init=1-epsilon)
    >>> rate_lf.setParamRule('omega', bins='1', is_const=True, value=1.0)
    >>> rate_lf.setAlignment(aln)
    >>> rate_lf.optimise(**optimiser_args)
    >>> tables = rate_lf.getStatistics(with_titles=True)
    >>> for table in tables:
    ...     if 'bin' in table.Title:
    ...         print table
    bin params
    ====================
    bin   bprobs   omega
    --------------------
      0     0.11    0.00
      1     0.89    1.00
    --------------------

We're also going to use the MLEs from the ``rate_lf`` model, since that nests within the more complex branch by rate-class model. This is unfortunately quite ugly compared with just using the annotated tree approach described above. It is currently necessary, however, due to a bug in constructing annotated trees for models with binned parameters.

.. doctest::
    
    >>> globals = [t for t in tables if 'global' in t.Title][0]
    >>> globals = dict(zip(globals.Header, globals.getRawData()[0]))
    >>> bin_params = [t for t in tables if 'bin' in t.Title][0]
    >>> rate_class_omegas = dict(bin_params.getRawData(['bin', 'omega']))
    >>> rate_class_probs = dict(bin_params.getRawData(['bin', 'bprobs']))
    >>> lengths = [t for t in tables if 'edge' in t.Title][0]
    >>> lengths = dict(lengths.getRawData(['edge', 'length']))

We now create the more complex model,

.. doctest::
    
     >>> rate_branch_lf = cnf.makeLikelihoodFunction(tree,
     ...             bins = ['0', '1', '2a', '2b'], digits=2, space=3)

and set from the nested null model the branch lengths,

.. doctest::
    
    >>> for branch, length in lengths.items():
    ...     rate_branch_lf.setParamRule('length', edge=branch, init=length)

GTR term MLES,

.. doctest::
    
    >>> for param, mle in globals.items():
    ...     rate_branch_lf.setParamRule(param, init=mle)

binned parameter values,

.. doctest::
    
    >>> rate_branch_lf.setParamRule('omega', bins=['0', '2a'], upper=1.0,
    ...                 init=rate_class_omegas['0'])
    >>> rate_branch_lf.setParamRule('omega', bins=['1', '2b'], is_const=True,
    ...                 value=1.0)
    >>> rate_branch_lf.setParamRule('omega', bins=['2a', '2b'],
    ...                    edges=['Chimpanzee', 'Human'], init=99,
    ...                    lower=1.0, upper=100.0, is_const=False)

and the bin probabilities.

.. doctest::
    
    >>> rate_branch_lf.setParamRule('bprobs',
    ...         init=[rate_class_probs['0']-epsilon,
    ...               rate_class_probs['1']-epsilon, epsilon, epsilon])

The result of these steps is to create a rate/branch model with initial parameter values that result in likelihood the same as the null.

.. doctest::
    
    >>> rate_branch_lf.setAlignment(aln)

This function can then be optimised as before. The results of one such optimisation are shown below. As you can see, the ``omega`` value for the '2a' and '2b' bins is at the upper bounds, indicating the model is not maximised in this case.

.. code-block:: python
    
    rate_branch_lf.optimise(**optimiser_args)
    print rate_branch_lf
    Likelihood Function Table
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
        Rhesus    2a     0.00
        Rhesus    2b     1.00
     Orangutan     0     0.00
     Orangutan     1     1.00
     Orangutan    2a     0.00
     Orangutan    2b     1.00
       Gorilla     0     0.00
       Gorilla     1     1.00
       Gorilla    2a     0.00
       Gorilla    2b     1.00
         Human     0     0.00
         Human     1     1.00
         Human    2a   100.00
         Human    2b   100.00
    Chimpanzee     0     0.00
    Chimpanzee     1     1.00
    Chimpanzee    2a   100.00
    Chimpanzee    2b   100.00...
