########
The apps
########

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

Sequence Alignment App
======================

``cogent3`` has bult-in alignment capabilities. The progressive alignment is a tree-based HMM aligner that can use any ``cogent3`` substitution model. There are default settings for the nucleotide, codon and protein cases, or you can customise the aligner settings.

.. toctree::
    :maxdepth: 2

    align-nucleotide
    align-codon
    align-protein
