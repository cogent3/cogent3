*************************************
Loading nucleotide, protein sequences
*************************************

.. author, Tony Walters

Using `LoadSeqs`
================

Simple case of loading a list of aligned nucleotide sequences in FASTA format.  **Note**: `LoadSeqs` by default expects aligned data to be passed to it (aligned=True) and generates an `Alignment` object.

.. doctest::
    
    >>> from cogent import LoadSeqs
    >>> seqs = ['>seq1','AATCG-A','>seq2','AATCGGA']
    >>> seqs_loaded = LoadSeqs(data=seqs)
    >>> print seqs_loaded
    >seq1
    AATCG-A
    >seq2
    AATCGGA
    <BLANKLINE>

Simple case of loading a list of aligned amino acid sequences in FASTA format, with and without molecule type specification.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> from cogent import DNA, PROTEIN
    >>> protein_seqs = ['>seq1','DEKQL-RG','>seq2','DDK--SRG']
    >>> proteins_loaded = LoadSeqs(data=protein_seqs)
    >>> print proteins_loaded
    >seq1
    DEKQL-RG
    >seq2
    DDK--SRG
    <BLANKLINE>
    >>> proteins_loaded = LoadSeqs(data=protein_seqs, moltype=PROTEIN)
    >>> print proteins_loaded
    >seq1
    DEKQL-RG
    >seq2
    DDK--SRG
    <BLANKLINE>

Simple case of loading aligned sequences from a file.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> fasta_file_seqs = LoadSeqs(filename="data/primate_cdx2_promoter.fasta")

Loading unaligned sequences from a file, specifying phylip format and unaligned sequences.
**Note**: With aligned=False, `LoadSeqs` creates a `SequenceCollection` object, not an `Alignment` object

.. doctest::

    >>> from cogent import LoadSeqs
    >>> clustal_file_seqs = LoadSeqs(filename="data/abglobin_aa.phylip", format="phylip", aligned=False)

Load a list of aligned nucleotide sequences, while specifying the DNA molecule type and stripping the comments from the label.  Stripping is accomplished by passing a function that removes everything after the first whitespace to the label_to_name parameter.

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> DNA_seqs = ['>sample1 Mus musculus','AACCTGC--C','>sample2 Gallus gallus','AAC-TGCAAC']
    >>> loaded_seqs = LoadSeqs(data=DNA_seqs, moltype=DNA, label_to_name=lambda x: x.split()[0])
    >>> print loaded_seqs
    >sample1
    AACCTGC--C
    >sample2
    AAC-TGCAAC
    <BLANKLINE>

Using alternative constructors for the `Alignment` object
=========================================================

An example of using an alternative constructor is given below.  A constructor is passed to the aligned parameter in lieu of True or False.

.. doctest::
    
    >>> from cogent import LoadSeqs
    >>> from cogent.core.alignment import DenseAlignment
    >>> seqs = ['>seq1','AATCG-A','>seq2','AATCGGA']
    >>> seqs_loaded = LoadSeqs(data=seqs,aligned=DenseAlignment)
    >>> print seqs_loaded
    >seq1
    AATCG-A
    >seq2
    AATCGGA
    <BLANKLINE>
