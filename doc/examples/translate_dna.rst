Translating DNA into protein
============================

.. sectionauthor:: Gavin Huttley

To translate a DNA alignment, read it in assigning the DNA alphabet. Note setting ``aligned = False`` is critical for loading sequences of unequal length. Different genetic codes are available in ``cogent.core.genetic_code``

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> al = LoadSeqs('data/test2.fasta', moltype=DNA, aligned = False)
    >>> pal = al.getTranslation()
    >>> print pal.toFasta()
    >DogFaced
    ARSQQNRWVETKETCNDRQT
    >HowlerMon
    ARSQHNRWAESEETCNDRQT
    >Human
    ARSQHNRWAGSKETCNDRRT
    >Mouse
    AVSQQSRWAASKGTCNDRQV
    >NineBande
    RQQSRWAESKETCNDRQT

To save this result to a file, use the ``writeToFile`` method.
