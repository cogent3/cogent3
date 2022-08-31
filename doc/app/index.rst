############
The ``apps``
############

``cogent3`` comes with pre-defined "apps". These are capabilities that can be used by themselves, or added together to define a pipeline. Apps that can be combined into a pipeline are called *composable*.

.. The goal is to allow using apps without requiring detailed understanding of programming. In other words, being able to use them as functions without writing loops or conditionals.

Overview
========

.. toctree::
    :maxdepth: 1

    warning
    app-overview
    dstore
    not-completed
    available-apps

The ``progressive_align`` App
=============================

``cogent3`` has a built-in progressive aligner. This is a tree-based HMM aligner that can use any ``cogent3`` substitution model.

.. toctree::
    :maxdepth: 1

    align-nucleotide
    align-codon
    align-protein


The ``model`` App
=================

``cogent3`` possesses unique capabilities for molecular evolutionary analyses. Of particular significance are our non-stationary model classes and the associated measurement of genetic distance from these. That said, we also implement most of the conventional substitution models. See ``cogent3.available_models()`` for a display of the built-in models.

You can use these models either through the general ``model`` app, or some task specific apps. Both are described below.

.. toctree::
    :maxdepth: 1

    evo-model
    evo-model-with-tree
    evo-model-timehet
    evo-extract-model-stats
    evo-ancestral-states
    evo-dt-nuc-model
    evo-tr-nuc-model
    evo-nstat-codon-model
    evo-tr-codon-model

The ``hypothesis`` App
======================

Takes two (or more) ``model`` instances and performs hypothesis test.

.. toctree::
    :maxdepth: 1

    evo-hypothesis

The ``natsel`` Apps
===================

These apps all employ codon models to evaluate different types of hypotheses regarding the mode of natural selection.

.. toctree::
    :maxdepth: 1

    evo-natsel_neutral
    evo-natsel_timehet
    evo-natsel_sitehet
    evo-natsel_zhang


Write your own Apps
===================

You can write your own apps.

.. toctree::
    :maxdepth: 1

    user_function
