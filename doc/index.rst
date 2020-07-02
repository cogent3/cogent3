#############
Documentation
#############


cogent3_ is a software library for genomic biology. It is a comprehensive framework for manipulation of genomic sequence data and for conducting molecular evolutionary analyses. It is distinguished by many unique built-in capabilities. Of particular note are our non-stationary Markov substitution models. These include `non-stationary nucleotide <https://www.ncbi.nlm.nih.gov/pubmed/25503772>`_, and `non-stationary codon <https://www.ncbi.nlm.nih.gov/pubmed/28175284>`_ models. 

.. dropdown:: Click to see an animation showing testing a hypothesis involving a non-stationary nucleotide process.

    .. raw:: html

        <object type="image/gif" data="../_static/gif/demo-fit-ns.gif" style="height: 400px; width: auto;"></object>

Plus, the ability to manipulate biological sequences by their annotations.

.. dropdown:: Click to see an animation showing defining and using sequence annotations.

    .. raw:: html

        <object type="image/gif" data="../_static/gif/demo-annotate.gif" style="height: 400px; width: auto;"></object>

.. panels::
    :header: bg-primary
    :footer: text-right
    
    ---

    ``cogent3`` apps
    ^^^^^^^^^^^^^^^^
    
    ``cogent3`` comes with pre-defined "apps" that simplify otherwise complex tasks. They provide capabilities that can be used by themselves, or added together to define a pipeline. They also simplify parallel execution of pipelines.

    +++++++++++

    .. link-button:: app/index
        :type: ref
        :text: …
        :classes: stretched-link

    ---

    Cookbook
    ^^^^^^^^

    The cookbook presents short code recipes targeted at specific problems.

    +++++++++++

    .. link-button:: cookbook/index
        :type: ref
        :text: …
        :classes: stretched-link

    ---

    Tutorials
    ^^^^^^^^^

    The tutorials present code for solving more extensive problems.
    
    +++++++++++

    .. link-button:: examples/index
        :type: ref
        :text: …
        :classes: stretched-link

    ---

    API
    ^^^

    The API for major ``cogent3`` objects
    
    +++++++++++

    .. link-button:: api/index
        :type: ref
        :text: …
        :classes: stretched-link

.. toctree::
    :hidden:
    :maxdepth: 3

    app/index
    cookbook/index
    examples/index
    api/index

.. _cogent3: https://cogent3.org
