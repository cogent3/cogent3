*******************
Accessing databases
*******************

.. Gavin Huttley, Kristian Rother, Patrick Yannul

NCBI
====

EUtils is a web service offered by the NCBI to access the sequence, literature and other databases by a special format of URLs. PyCogent offers an interface to construct the URLs and retrieve the results in text format.

For sequence
------------

Fetching FASTA or Genpept sequences from NCBI using GI's
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.db.ncbi import EFetch
    >>> fasta = EFetch(id='1234567,459567',rettype='fasta').read()
    >>> genpept = EFetch(id='1234567,459567',rettype='gp').read()

Retrieving GenPept files from NCBI via Eutils
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We query for just one accession to illustrate the process. A more general query can be executed by replacing ``'BAB52044`` with ``'"lysyl tRNA-synthetase"[ti] AND bacteria[orgn]'`` in the snippet below.

.. doctest::

    >>> from cogent.db.ncbi import EUtils
    >>> e = EUtils(numseqs=100, db='protein', rettype='gp')
    >>> result = e['BAB52044']
    >>> print result.read()
    LOCUS       BAB52044                 548 aa            linear   BCT 16-MAY-2009
    DEFINITION  lysyl tRNA synthetase [Mesorhizobium loti MAFF303099].
    ACCESSION   BAB52044
    VERSION     BAB52044.1  GI:14025444
    DBSOURCE    accession BA000012.4
    KEYWORDS    .
    SOURCE      Mesorhizobium loti MAFF303099...

Retrieving and parsing GenBank entries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.parse.genbank import RichGenbankParser
    >>> from cogent.db.ncbi import EUtils
    >>> e = EUtils(numseqs=100, db='protein', rettype='gp')
    >>> result = e['"lysyl tRNA-synthetase"[ti] AND bacteria[orgn]']
    >>> parser = RichGenbankParser(result.readlines())
    >>> gb = [(accession, seq) for accession, seq in parser]

Printing the resulting list (``gb``) will generate output like:

.. code-block:: python
    
    [('AAA83071', Sequence(MSEQHAQ... 505)), ('ACS40931', ...

For other
---------

.. OMIM, ??

Retrieving PubMed abstracts from NCBI via EUtils
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> from cogent.db.ncbi import EUtils
    >>> e = EUtils(db='pubmed',rettype='brief')
    >>> result = e['Simon Easteal'].read()
    >>> print result
    <BLANKLINE>
    1: Yap VB et al. Estimates of the effect of na...[PMID: 19815689] 
    <BLANKLINE>
    2: Cherbuin N et al. Risk factors of transition fr...[PMID: 19628940] ...
    >>> e = EUtils(db='pubmed',rettype='abstract')
    >>> result = e['19815689'].read()
    >>> print result
    <BLANKLINE>
    1: Mol Biol Evol. 2009 Oct 8; [Epub ahead of print] 
    <BLANKLINE>
    Estimates of the effect of natural selection on protein coding content.
    <BLANKLINE>
    Yap VB, Lindsay H, Easteal S, Huttley G.
    <BLANKLINE>
    Department of Statistics and Applied Probability, National University of
    Singapore, Singapore.
    <BLANKLINE>
    Analysis of natural selection is key to understanding many core biological
    processes, including the emergence of competition, co-operation, and complexity,...

Retrieving PubMed abstracts via PMID
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.db.ncbi import EUtils
    >>> e = EUtils(db='pubmed',rettype='abstract')
    >>> result = e['14983078'].read()

KEGG
====

Complete genomes
----------------

*To be written.*

Orthologs
---------

*To be written.*

Functional assignments
----------------------

*To be written.*

Pathway assignments
-------------------

*To be written.*

Ensembl
=======

Connecting
----------

*To be written.*

.. Hosts and species

Get genomic features
--------------------

*To be written.*

Get alignments
--------------

*To be written.*

Get SNPs
--------

*To be written.*

PDB
===

For structures
--------------

*To be written.*

Rfam
====

For RNA secondary structures, alignments, functions
---------------------------------------------------

*To be written.*

GoldenPath (not yet implemented)
================================

*To be written.*

Whole-genome alignments, orthologs, annotation tracks
-----------------------------------------------------

*To be written.*
