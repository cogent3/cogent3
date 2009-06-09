Drawing dendrograms and saving to PDF
=====================================

.. sectionauthor:: Gavin Huttley

From cogent import all the components we need.

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree
    >>> from cogent.evolve.models import MG94HKY
    >>> from cogent.draw import dendrogram

Do a model, see the neutral test example for more details of this

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> aln = LoadSeqs("data/long_testseqs.fasta")
    >>> t = LoadTree("data/test.tree")
    >>> sm = MG94HKY()
    >>> nonneutral_lf = sm.makeLikelihoodFunction(t)
    >>> nonneutral_lf.setParamRule("omega", is_independent = True)
    >>> nonneutral_lf.setAlignment(aln)
    >>> nonneutral_lf.optimise(show_progress=False)

We will draw two different dendrograms -- one with branch lengths contemporaneous, the other where length is scaled.

Specify the dimensions of the canvas in pixels

.. doctest::

    >>> height, width = 500, 700

Dendrogram with branch lengths not proportional
-----------------------------------------------

.. doctest::

    >>> np = dendrogram.ContemporaneousDendrogram(nonneutral_lf.tree)
    >>> np.drawToPDF('tree-unscaled.pdf' , width, height, stroke_width=2.0,
    ... show_params = ['r'], label_template = "%(r).2g", shade_param = 'r',
    ... max_value = 1.0, show_internal_labels=False, font_size = 10,
    ... scale_bar = None, use_lengths=False)

Dendrogram with branch lengths proportional
-------------------------------------------

.. doctest::

    >>> p = dendrogram.SquareDendrogram(nonneutral_lf.tree)
    >>> p.drawToPDF('tree-scaled.pdf', width, height, stroke_width=2.0,
    ... shade_param = 'r', max_value = 1.0, show_internal_labels=False,
    ... font_size = 10)

Separating the analysis and visualisation steps
-----------------------------------------------

It's typically better to not have the analysis and drawing code in the same script, since drawing involves frequent iterations. This requires saving a tree for later reuse. This can be done using an annotated tree, which looks just like a tree, but has the maximum-likelihood parameter estimates attached to each tree edge. The tree must be saved in xml format to preserve the parameter estimates. The annotated tree is obtained from the likelihood function and saved to file specifying the format with the .xml suffix. This file can then be loaded using the standard ``LoadTree`` method in a separate script and used for drawing.

.. doctest::

    >>> at = nonneutral_lf.getAnnotatedTree()
    >>> at.writeToFile('annotated_tree.xml')

.. we clean up after ourselves, deleting the file

.. doctest::
    :hide:

    >>> import os
    >>> for file_name in 'tree-scaled.pdf', 'tree-unscaled.pdf', 'annotated_tree.xml':
    ...     os.remove(file_name)
