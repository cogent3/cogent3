*************************************************
Using likelihood to perform evolutionary analyses
*************************************************

Specifying substitution models
==============================

Canned models
-------------

MotifChange and predicates
--------------------------

Rate heterogeneity models
-------------------------

Specifying likelihood functions
===============================

Scoping parameters on trees
---------------------------

Specifying parameter values
---------------------------

.. constant, bounds, initial

Specifying rate heterogeneity functions
---------------------------------------

Specifying Phylo-HMMs
---------------------

Fitting likelihood functions
============================

Choice of optimisers
--------------------

Checkpointing runs
------------------

How to check your optimisation was successful.
----------------------------------------------

.. Try again, use global optimisation, check maximum numbers of calculations not exceeded.

Getting statistics out of likelihood functions
==============================================

.. the annotated tree, the tables, getParamValue

Testing hypotheses
==================

.. LRTs, assuming chisq, bootstrapping, randomisation

Determining confidence intervals on MLEs
========================================

Saving results
==============

Visualising statistics on trees
===============================

Reconstructing ancestral sequences
==================================

.. most likely ancestors, the complete posterior probabilities

Tips for improved performance
=============================

Sequentially build the fitting
------------------------------

.. start with null, then modify lf to alternate. Don't forget to record the values you need.

.. how to specify the alt so it is the null for rate heterogeneity models

Sampling
--------

.. using a subset of data

