**************************************
Evolutionary analysis using likelihood
**************************************

Specifying substitution models
==============================

Canned models
-------------

Many standard evolutionary models come pre-defined in the ``cogent.evolve.models`` module.

The available nucleotide, codon and protein models are

.. doctest::
    
    >>> from cogent.evolve import models
    >>> print models.nucleotide_models
    ['JC69', 'K80', 'F81', 'HKY85', 'TN93', 'GTR']
    >>> print models.codon_models
    ['CNFGTR', 'CNFHKY', 'MG94HKY', 'MG94GTR', 'GY94', 'H04G', 'H04GK', 'H04GGK']
    >>> print models.protein_models
    ['DSO78', 'AH96', 'AH96_mtmammals', 'JTT92', 'WG01']

While those values are strings, a function of the same name exists within the module so creating the substitution models requires only calling that function. I demonstrate that for a nucleotide model here.

.. doctest::
    
    >>> from cogent.evolve.models import F81
    >>> sub_mod = F81()

We'll be using these for the examples below.

Rate heterogeneity models
-------------------------

We illustrate this for the gamma distributed case using examples of the canned models displayed above. Creating rate heterogeneity variants of the canned models can be done by using optional arguments that get passed to the substitution model class.

For nucleotide
^^^^^^^^^^^^^^

We specify a general time reversible nucleotide model with gamma distributed rate heterogeneity.

.. doctest::
    
    >>> from cogent.evolve.models import GTR
    >>> sub_mod = GTR(with_rate=True, distribution='gamma')
    >>> print sub_mod
    <BLANKLINE>
    Nucleotide ( name = 'GTR'; type = 'None'; params = ['A/G', 'A/T', 'A/C', 'C/T', 'C/G']; number of motifs = 4; motifs = ['T', 'C', 'A', 'G'])
    <BLANKLINE>

For codon
^^^^^^^^^

We specify a conditional nucleotide frequency codon model with nucleotide general time reversible parameters and a parameter for the ratio of nonsynonymous to synonymous substitutions (omega) with gamma distributed rate heterogeneity.

.. doctest::
    
    >>> from cogent.evolve.models import CNFGTR
    >>> sub_mod = CNFGTR(with_rate=True, distribution='gamma')
    >>> print sub_mod
    <BLANKLINE>
    Codon ( name = 'CNFGTR'; type = 'None'; params = ['A/G', 'A/C', 'C/T', 'A/T', 'C/G', 'omega']; ...

For protein
^^^^^^^^^^^

We specify a Jones, Taylor and Thornton 1992 empirical protein substitution model with gamma distributed rate heterogeneity.

.. doctest::
    
    >>> from cogent.evolve.models import JTT92
    >>> sub_mod = JTT92(with_rate=True, distribution='gamma')
    >>> print sub_mod
    <BLANKLINE>
    Empirical ( name = 'JTT92'; type = 'None'; number of motifs = 20; motifs = ['A', 'C'...

Specifying likelihood functions
===============================

Making a likelihood function
----------------------------

You start by specifying a substitution model and use that to construct a likelihood function for a specific tree.

.. doctest::
    
    >>> from cogent import LoadTree
    >>> from cogent.evolve.models import F81
    >>> sub_mod = F81()
    >>> tree = LoadTree(treestring='(a,b,(c,d))')
    >>> lf = sub_mod.makeLikelihoodFunction(tree)

Providing an alignment to a likelihood function
-----------------------------------------------

You need to load an alignment and then provide it a likelihood function. I construct very simple trees and alignments for this example.

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import F81
    >>> sub_mod = F81()
    >>> tree = LoadTree(treestring='(a,b,(c,d))')
    >>> lf = sub_mod.makeLikelihoodFunction(tree)
    >>> aln = LoadSeqs(data=[('a', 'ACGT'), ('b', 'AC-T'), ('c', 'ACGT'),
    ...                      ('d', 'AC-T')])
    ...                     
    >>> lf.setAlignment(aln)

Scoping parameters on trees
---------------------------

For many evolutionary analyses, it's desirable to allow different branches on a tree to have different values of a parameter. We show this for a simple codon model case here where we want the great apes (the clade that includes human and orangutan) to have a different value of the ratio of nonsynonymous to synonymous substitutions. This parameter is identified in the precanned ``CNFGTR`` model as ``omega``.

.. doctest::
    
    >>> from cogent import LoadTree
    >>> from cogent.evolve.models import CNFGTR
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> print tree.asciiArt()
              /-Galago
             |
    -root----|--HowlerMon
             |
             |          /-Rhesus
              \edge.3--|
                       |          /-Orangutan
                        \edge.2--|
                                 |          /-Gorilla
                                  \edge.1--|
                                           |          /-Human
                                            \edge.0--|
                                                      \-Chimpanzee
    >>> sm = CNFGTR()
    >>> lf = sm.makeLikelihoodFunction(tree, digits=2)
    >>> lf.setParamRule('omega', tip_names=['Human', 'Orangutan'], outgroup_name='Galago', is_clade=True, init=0.5)

We've set an *initial* value for this clade so that the edges affected by this rule are evident below.

.. doctest::
    
    >>> print lf
    Likelihood Function Table
    ====================================
     A/C     A/G     A/T     C/G     C/T
    ------------------------------------
    1.00    1.00    1.00    1.00    1.00
    ------------------------------------
    =======================================
          edge    parent    length    omega
    ---------------------------------------
        Galago      root      1.00     1.00
     HowlerMon      root      1.00     1.00
        Rhesus    edge.3      1.00     1.00
     Orangutan    edge.2      1.00     0.50
       Gorilla    edge.1      1.00     0.50
         Human    edge.0      1.00     0.50
    Chimpanzee    edge.0      1.00     0.50
        edge.0    edge.1      1.00     0.50
        edge.1    edge.2      1.00     0.50
        edge.2    edge.3      1.00     1.00
        edge.3      root      1.00     1.00
    ---------------------------------------...

A more extensive description of capabilities is in :ref:`scope-params-on-trees`.

Specifying parameter values
---------------------------

Specifying a parameter as constant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This means the parameter will not be modified during likelihood maximisation. We show this here by making the ``omega`` parameter constant at the value 1 -- essentially the condition of selective neutrality.

.. doctest::
    
    >>> from cogent import LoadTree
    >>> from cogent.evolve.models import CNFGTR
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> sm = CNFGTR()
    >>> lf = sm.makeLikelihoodFunction(tree, digits=2)
    >>> lf.setParamRule('omega', is_constant=True)

Providing a starting value for a parameter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This can be useful to improve performance, the closer you are to the maximum likelihood estimator the quicker optimisation will be.

.. doctest::
    
    >>> from cogent import LoadTree
    >>> from cogent.evolve.models import CNFGTR
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> sm = CNFGTR()
    >>> lf = sm.makeLikelihoodFunction(tree, digits=2)
    >>> lf.setParamRule('omega', init=0.1)

Setting bounds for optimising a function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This can be useful for stopping optimisers from getting stuck in a bad part of parameter space.

.. doctest::
    
    >>> from cogent import LoadTree
    >>> from cogent.evolve.models import CNFGTR
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> sm = CNFGTR()
    >>> lf = sm.makeLikelihoodFunction(tree, digits=2)
    >>> lf.setParamRule('omega', init=0.1, lower=1e-9, upper=20.0)

Specifying rate heterogeneity functions
---------------------------------------

We extend the simple gamma distributed rate heterogeneity case for nucleotides from above to construction of the actual likelihood function. We do this for 4 bins and constraint the bin probabilities to be equal.

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import GTR
    >>> sm = GTR(with_rate=True, distribution='gamma')
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> lf = sm.makeLikelihoodFunction(tree, bins=4, digits=2)
    >>> lf.setParamRule('bprobs', is_constant=True)

For more detailed discussion of defining and using these models see :ref:`rate-heterogeneity`.

Specifying Phylo-HMMs
---------------------

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import GTR
    >>> sm = GTR(with_rate=True, distribution='gamma')
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> lf = sm.makeLikelihoodFunction(tree, bins=4, sites_independent=False,
    ...                                 digits=2)
    >>> lf.setParamRule('bprobs', is_constant=True)

For more detailed discussion of defining and using these models see :ref:`rate-heterogeneity-hmm`.

Fitting likelihood functions
============================

Choice of optimisers
--------------------

There are 2 types of optimiser: simulated annealing, a *global* optimiser; and Powell, a *local* optimiser. The simulated annealing method is slow compared to Powell and in general Powell is an adequate choice. I setup  a simple nucleotide model to illustrate these.

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import F81
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> sm = F81()
    >>> lf = sm.makeLikelihoodFunction(tree, digits=3, space=2)
    >>> lf.setAlignment(aln)

The default is to use the simulated annealing optimiser followed by Powell.

.. doctest::
    
    >>> lf.optimise(show_progress=False)

We can specify just using the local optimiser. To do so, it's recommended to set the ``max_restarts`` argument since this provides a mechanism for Powell to attempt restarting the optimisation from slightly different sport which can help in overcoming local maxima.

.. doctest::
    
    >>> lf.optimise(local=True, max_restarts=5, show_progress=False)

We might want to do crude simulated annealing following by more rigorous Powell.

.. doctest::
    
    >>> lf.optimise(show_progress=False, global_tolerance=1.0, tolerance=1e-8,
    ...              max_restarts=5)

Checkpointing runs
------------------

See :ref:`checkpointing-optimisation`.

How to check your optimisation was successful.
----------------------------------------------

There is no guarantee that an optimised function has achieved a global maximum. We can, however, be sure that a maximum was achieved by validating that the optimiser stopped because the specified tolerance condition was met, rather than exceeding the maximum number of evaluations. The latter number is set to ensure optimisation doesn't proceed endlessly. If the optimiser exited because this limit was exceeded you can be sure that the function **has not** been successfully optimised.

We can monitor this situation using the ``limit_action`` argument to ``optimise``. Providing the value ``raise`` causes an exception to be raised if this condition occurs, as shown below. Providing ``warn`` (default) instead will cause a warning message to be printed to screen but execution will continue. The value ``ignore`` hides any such message.

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import F81
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> sm = F81()
    >>> lf = sm.makeLikelihoodFunction(tree, digits=3, space=2)
    >>> lf.setAlignment(aln)
    >>> max_evals = 10
    >>> lf.optimise(show_progress=False, limit_action='raise',
    ...              max_evaluations=max_evals, return_calculator=True)
    ... 
    Traceback (most recent call last):
    ArithmeticError: FORCED EXIT from optimiser after 10 evaluations

.. note:: We recommend using ``limit_action='raise'`` and catching the ``ArithmeticError`` error explicitly. You really shouldn't be using results from such an optimisation run.

Getting statistics out of likelihood functions
==============================================

Model fit statistics
--------------------

Log likelihood and number of free parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import GTR
    >>> sm = GTR()
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> lf = sm.makeLikelihoodFunction(tree)
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> lf.setAlignment(aln)

We get the log-likelihood and the number of free parameters.

.. doctest::
    
    >>> lnL = lf.getLogLikelihood()
    >>> print lnL
    -24601.9...
    >>> nfp = lf.getNumFreeParams()
    >>> print nfp
    16

.. warning:: The number of free parameters (nfp) refers only to the number of parameters that were modifiable by the optimiser. Typically, the degrees-of-freedom of a likelihood ratio test statistic is computed as the difference in nfp between models. This will not be correct for models in which boundary conditions exist (rate heterogeneity models where a parameter value boundary is set between bins).

Information theoretic measures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Aikake Information Criterion
""""""""""""""""""""""""""""

.. note:: this measure only makes sense when the model has been optimised, a step I'm skipping here in the interests of speed.

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import GTR
    >>> sm = GTR()
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> lf = sm.makeLikelihoodFunction(tree)
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> lf.setAlignment(aln)
    >>> AIC = lf.getAic()
    >>> AIC
    49235.869...

We can also get the second-order AIC.

.. doctest::
    
    >>> AICc = lf.getAic(second_order=True)
    >>> AICc
    49236.064...

Bayesian Information Criterion
""""""""""""""""""""""""""""""

.. note:: this measure only makes sense when the model has been optimised, a step I'm skipping here in the interests of speed.

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import GTR
    >>> sm = GTR()
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> lf = sm.makeLikelihoodFunction(tree)
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> lf.setAlignment(aln)
    >>> BIC = lf.getBic()
    >>> BIC
    49330.9475...

Getting maximum likelihood estimates
------------------------------------

We fit the model defined in the previous section and use that in the following.

One at a time
^^^^^^^^^^^^^

We get the statistics out individually. We get the ``length`` for the Human edge and the exchangeability parameter ``A/G``.

.. doctest::
    
    >>> lf.optimise(local=True, show_progress=False)
    >>> a_g = lf.getParamValue('A/G')
    >>> print a_g
    5.25...
    >>> human = lf.getParamValue('length', 'Human')
    >>> print human
    0.006...

Just the motif probabilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> mprobs = lf.getMotifProbs()
    >>> print mprobs
    ====================================
         T         C         A         G
    ------------------------------------
    0.2406    0.1742    0.3757    0.2095
    ------------------------------------

On the tree object
^^^^^^^^^^^^^^^^^^

If written to file in xml format, then model parameters will be saved. This can be useful for later plotting or recreating likelihood functions.

.. doctest::
    
    >>> annot_tree = lf.getAnnotatedTree()
    >>> print annot_tree.getXML() #doctest: +SKIP
    <?xml version="1.0"?>
    <clade>
      <clade>
         <name>Galago</name>
         <param><name>A/G</name><value>5.25342689214</value></param>
         <param><name>A/C</name><value>1.23159157151</value></param>
         <param><name>C/T</name><value>5.97001104267</value></param>
         <param><name>length</name><value>0.173114172705</value></param>...

.. warning:: This method fails for some rate-heterogeneity models.

As tables
^^^^^^^^^

.. doctest::
    
    >>> tables = lf.getStatistics(with_motif_probs=True, with_titles=True)
    >>> for table in tables:
    ...     if 'global' in table.Title:
    ...         print table
    global params
    ==============================================
       A/C       A/G       A/T       C/G       C/T
    ----------------------------------------------
    1.2316    5.2534    0.9585    2.3159    5.9700
    ----------------------------------------------

Testing hypotheses
==================

Using likelihood ratio tests
----------------------------

We test the molecular clock hypothesis for human and chimpanzee lineages. The null has these two branches constrained to be equal.

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import F81
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> sm = F81()
    >>> lf = sm.makeLikelihoodFunction(tree, digits=3, space=2)
    >>> lf.setAlignment(aln)
    >>> lf.setParamRule('length', tip_names=['Human', 'Chimpanzee'],
    ...         outgroup_name='Galago', is_clade=True, is_independent=False)
    ...                 
    >>> lf.setName('Null Hypothesis')
    >>> lf.optimise(local=True, show_progress=False)
    >>> null_lnL = lf.getLogLikelihood()
    >>> null_nfp = lf.getNumFreeParams()
    >>> print lf
    Null Hypothesis
    ==========================
          edge  parent  length
    --------------------------
        Galago    root   0.167
     HowlerMon    root   0.044
        Rhesus  edge.3   0.021
     Orangutan  edge.2   0.008
       Gorilla  edge.1   0.002
         Human  edge.0   0.004
    Chimpanzee  edge.0   0.004
        edge.0  edge.1   0.000...

The alternate allows the human and chimpanzee branches to differ by just setting all lengths to be independent.

.. doctest::
    
    >>> lf.setParamRule('length', is_independent=True)
    >>> lf.setName('Alt Hypothesis')
    >>> lf.optimise(local=True, show_progress=False)
    >>> alt_lnL = lf.getLogLikelihood()
    >>> alt_nfp = lf.getNumFreeParams()
    >>> print lf
    Alt Hypothesis
    ==========================
          edge  parent  length
    --------------------------
        Galago    root   0.167
     HowlerMon    root   0.044
        Rhesus  edge.3   0.021
     Orangutan  edge.2   0.008
       Gorilla  edge.1   0.002
         Human  edge.0   0.006
    Chimpanzee  edge.0   0.003
        edge.0  edge.1   0.000...

We import the function for computing the probability of a chi-square test statistic, compute the likelihood ratio test statistic, degrees of freedom and the corresponding probability.

.. doctest::
    
    >>> from cogent.maths.stats import chisqprob
    >>> LR = 2 * (alt_lnL - null_lnL) # the likelihood ratio statistic
    >>> df = (alt_nfp - null_nfp) # the test degrees of freedom
    >>> p = chisqprob(LR, df)
    >>> print 'LR=%.4f ; df = %d ; p=%.4f' % (LR, df, p)
    LR=3.3294 ; df = 1 ; p=0.0681

By parametric bootstrapping
---------------------------

If we can't rely on the asymptotic behaviour of the LRT, e.g. due to small alignment length, we can use a parametric bootstrap. Convenience functions for that are described in more detail here :ref:`parametric-bootstrap`.

In general, however, this capability derives from the ability of any defined ``evolve`` likelihood function to simulate an alignment. This property is provided as ``simulateAlignment`` method on likelihood function objects.

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import F81
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> sm = F81()
    >>> lf = sm.makeLikelihoodFunction(tree, digits=3, space=2)
    >>> lf.setAlignment(aln)
    >>> lf.setParamRule('length', tip_names=['Human', 'Chimpanzee'],
    ...         outgroup_name='Galago', is_clade=True, is_independent=False)
    ...                 
    >>> lf.setName('Null Hypothesis')
    >>> lf.optimise(local=True, show_progress=False)
    >>> sim_aln = lf.simulateAlignment()
    >>> print repr(sim_aln)
    7 x 2814 dna alignment: Gorilla...

Determining confidence intervals on MLEs
========================================

The profile method is used to calculate a confidence interval for a named parameter. We show it here for a global substitution model exchangeability parameter (*kappa*, the ratio of transition to transversion rates) and for an edge specific parameter (just the human branch length).

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import HKY85
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> sm = HKY85()
    >>> lf = sm.makeLikelihoodFunction(tree)
    >>> lf.setAlignment(aln)
    >>> lf.optimise(local=True, show_progress=False)
    >>> kappa_lo, kappa_mle, kappa_hi = lf.getParamInterval('kappa')
    >>> print "lo=%.2f ; mle=%.2f ; hi = %.2f" % (kappa_lo, kappa_mle, kappa_hi)
    lo=3.78 ; mle=4.44 ; hi = 5.22
    >>> human_lo, human_mle, human_hi = lf.getParamInterval('length', 'Human')
    >>> print "lo=%.2f ; mle=%.2f ; hi = %.2f" % (human_lo, human_mle, human_hi)
    lo=0.00 ; mle=0.01 ; hi = 0.01

Saving results
==============

Use either the annotated tree or statistics tables to obtain objects that can easily be written to file.

Visualising statistics on trees
===============================

We look at the distribution of ``omega`` from the CNF codon model family across different primate lineages. We allow each edge to have an independent value for ``omega``.

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import CNFGTR
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> sm = CNFGTR()
    >>> lf = sm.makeLikelihoodFunction(tree, digits=2, space=2)
    >>> lf.setParamRule('omega', is_independent=True, upper=10.0)
    >>> lf.setAlignment(aln)
    >>> lf.optimise(show_progress=False, local=True)
    >>> print lf
    Likelihood Function Table
    ============================
     A/C   A/G   A/T   C/G   C/T
    ----------------------------
    1.07  3.88  0.79  1.96  4.09
    ----------------------------
    =================================
          edge  parent  length  omega
    ---------------------------------
        Galago    root    0.53   0.85
     HowlerMon    root    0.14   0.71
        Rhesus  edge.3    0.07   0.58
     Orangutan  edge.2    0.02   0.49
       Gorilla  edge.1    0.01   0.43
         Human  edge.0    0.02   2.44
    Chimpanzee  edge.0    0.01   2.28
        edge.0  edge.1    0.00   0.01
        edge.1  edge.2    0.01   0.55
        edge.2  edge.3    0.04   0.33
        edge.3    root    0.02   1.10...

We need an annotated tree object to do the drawing, we write this out to an XML formatted file so it can be reloaded for later reuse.

.. doctest::
    
    >>> annot_tree = lf.getAnnotatedTree()
    >>> annot_tree.writeToFile('result_tree.xml')

We first import an unrooted dendrogram and then generate a heat mapped image to file where edges are colored red by the magnitude of ``omega`` with maximal saturation when ``omega=1``.

.. doctest::
    
    >>> from cogent.draw.dendrogram import ContemporaneousDendrogram
    >>> dend = ContemporaneousDendrogram(annot_tree)
    >>> fig = dend.makeFigure(height=6, width=6, shade_param='omega',
    ...                      max_value=1.0, stroke_width=2)
    >>> fig.savefig('omega_heat_map.png')

Reconstructing ancestral sequences
==================================

We first fit a likelihood function.

.. doctest::
    
    >>> from cogent import LoadTree, LoadSeqs
    >>> from cogent.evolve.models import F81
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> sm = F81()
    >>> lf = sm.makeLikelihoodFunction(tree, digits=3, space=2)
    >>> lf.setAlignment(aln)
    >>> lf.optimise(show_progress=False, local=True)

We then get the most likely ancestral sequences.

.. doctest::
    
    >>> ancestors = lf.likelyAncestralSeqs()
    >>> print ancestors
    >root
    TGTGGCACAAATACTCATGCCAGCTCATTACAGCA...

Or we can get the posterior probabilities (returned as a ``DictArray``) of sequence states at each node.

.. doctest::
    
    >>> ancestral_probs = lf.reconstructAncestralSeqs()
    >>> print ancestral_probs['root']
    ============================================
                 T         C         A         G
    --------------------------------------------
       0    0.1816    0.0000    0.0000    0.0000
       1    0.0000    0.0000    0.0000    0.1561
       2    0.1816    0.0000    0.0000    0.0000
       3    0.0000    0.0000    0.0000    0.1561...

Tips for improved performance
=============================

Sequentially build the fitting
------------------------------

There's nothing that improves performance quite like being close to the maximum likelihood values. So using the ``setParamRule`` method to provide good starting values can be very useful. As this can be difficult to do one easy way is to build simpler models that are nested within the one you're interested in. Fitting those models and then relaxing constraints until you’re at the parameterisation of interest can markedly improve optimisation speed.

Being able to save results to file allows you to do this between sessions.

Sampling
--------

If you're dealing with a very large alignment, another approach is to use a subset of the alignment to fit the model then try fitting the entire alignment. The alignment method does have an method to facilitate this approach. The following samples 99 codons without replacement.

.. doctest::
    
    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> smpl = aln.sample(n=99, with_replacement=False, motif_length=3)
    >>> len(smpl)
    297

While this samples 99 nucleotides without replacement.

.. doctest::
    
    >>> smpl = aln.sample(n=99, with_replacement=False)
    >>> len(smpl)
    99

.. following cleans up files

.. doctest::
    :hide:

    >>> from cogent.util.misc import remove_files
    >>> remove_files(['result_tree.xml', 'omega_heat_map.png'],
    ...               error_on_missing=False)
