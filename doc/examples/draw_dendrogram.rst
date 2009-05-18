Drawing dendrograms and saving to PDF
=====================================

From cogent import all the components we need.

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree
    >>> from cogent.evolve.models import Y98
    >>> from cogent.draw import dendrogram

Do a model, see the neutral test example for more details of this

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> al = LoadSeqs("data/test.paml")
    >>> t = LoadTree("data/test.tree")
    >>> sm = Y98()
    >>> nonneutral_lf = sm.makeLikelihoodFunction(t)
    >>> nonneutral_lf.setParamRule("omega", is_independent = 1)
    >>> nonneutral_lf.setAlignment(al)
    >>> nonneutral_lf.optimise(tolerance = 1.0)
    Outer loop = 0...
    >>> nonneutral_lf.optimise(local = True)
        Number of function evaluations = 1; current F = 139...

We will draw two different dendrograms -- one with branch lengths contemporaneous, the other where length is scaled.

Specify the dimensions of the canvas in pixels

.. doctest::

    >>> height, width = 500, 700

Dendrogram with branch lengths not proportional

.. doctest::

    >>> np = dendrogram.ContemporaneousDendrogram(nonneutral_lf.tree)
    >>> np.drawToPDF('tree-unscaled.pdf' , width, height, stroke_width=2.0,
    ... show_params = ['r'], label_template = "%(r).2g", shade_param = 'r',
    ... max_value = 1.0, show_internal_labels=False, font_size = 10,
    ... scale_bar = None, use_lengths=False)

Dendrogram with branch lengths proportional

.. doctest::

    >>> p = dendrogram.SquareDendrogram(nonneutral_lf.tree)
    >>> p.drawToPDF('tree-scaled.pdf', width, height, stroke_width=2.0,
    ... shade_param = 'r', max_value = 1.0, show_internal_labels=False,
    ... font_size = 10)

To save a tree for later reuse, either for analysis of drawing can be done using an annotated tree, which looks just like a tree, but has the maximum-likelihood parameter estimates attached to each tree edge. This tree can be saved in xml format, which preserve these parameter estimates. The annotated tree is obtained from the likelihood function with following command.

.. doctest::

    >>> at = nonneutral_lf.getAnnotatedTree()

Saving this to file is done using the normal ``writeToFile`` method, specifying a filename with the .xml suffix.
