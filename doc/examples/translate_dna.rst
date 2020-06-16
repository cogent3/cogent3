Translating DNA into protein
============================

.. sectionauthor:: Gavin Huttley

To translate a DNA alignment, read it in assigning the DNA alphabet. Note setting ``aligned=False`` is critical for loading sequences of unequal length. Different genetic codes are available in ``cogent3.core.genetic_code``

.. jupyter-execute::
    :linenos:

    from cogent3 import load_unaligned_seqs

    al = load_unaligned_seqs("data/test2.fasta", moltype="dna")
    pal = al.get_translation()
    pal.to_fasta()
