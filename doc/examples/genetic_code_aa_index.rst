Compute the effect of a nucleotide substitution on residue polarity in two different genetic codes using GeneticCode and AAIndex
================================================================================================================================

.. sectionauthor:: Greg Caporaso

This document illustrates how to work with a genetic code object, and compare two different genetic codes. Here we compare the change in residue polarity, as judged by the Woese Polarity Requirement index (Woese 1973), resulting from a nucleotide substitution if the sequence is translated with the standard nuclear genetic code, or the vertebrate mitochondrial genetic code. 

First, we load the genetic code objects and look at how they differ from one another. 

.. doctest::

    >>> from cogent.core.genetic_code import GeneticCode
    >>> code = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    >>> standard_nuclear_genetic_code = GeneticCode(code)
    >>> code = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG'
    >>> vertebrate_mitochondrial_genetic_code = GeneticCode(code)
    >>> standard_nuclear_genetic_code == vertebrate_mitochondrial_genetic_code
    False

We'll make some synonyms for the objects for simplicity, and then look at the differences between the two codes:

.. doctest::

    >>> ngc = standard_nuclear_genetic_code
    >>> mgc = vertebrate_mitochondrial_genetic_code

    >>> differences = ngc.changes(mgc).items()
    >>> differences.sort()
    >>> differences
    [('AGA', 'R*'), ('AGG', 'R*'), ('ATA', 'IM'), ('TGA', '*W')]



Next, let's load the Woese Polar Requirement ``AAIndex`` data, and find the effect of an ATA to ATG substitution with each of the two ``GeneticCode`` objects.

.. doctest::
    
    >>> from cogent.parse.aaindex import getWoeseDistanceMatrix
    >>> woese_distance_matrix = getWoeseDistanceMatrix()
    >>> woese_distance_matrix[ngc['ATA']][ngc['ATG']]
    0.39999999999999947
    >>> woese_distance_matrix[mgc['ATA']][mgc['ATG']]
    0.0

This illustrates that there is a difference in residue polarity associated with substitution only in the standard nuclear code (where ATA to ATG translates to an isoleucine to methionine substitution). In the vertebrate mitochondrial code, ATA to ATG is a synonymous substitution. Calculations of this type were central to [1]_ which presents the study that these modules were initially developed for.

``GeneticCode`` objects can also be used to translate DNA sequences (where asterisks in the results refer to stop-translation characters):

.. doctest::
    
    >>> dna = "AAACGCTGTGTGTGAGATGAAAAA"
    >>> ngc.translate(dna)
    'KRCV*DEK'
    >>> mgc.translate(dna)
    'KRCVWDEK'

The standard nuclear genetic code can also be loaded as ``DEFAULT``:

.. doctest::
    
    >>> from cogent.core.genetic_code import DEFAULT
    >>> DEFAULT == standard_nuclear_genetic_code
    True

**Citations**

.. [1] Caporaso, Yarus, and Knight. *Error minimization and coding triplet/binding site associations are independent features of the canonical genetic code.* J Mol Evol, 61(5):597-607, 2005.

