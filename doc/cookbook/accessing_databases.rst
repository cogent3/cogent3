*******************
Accessing databases
*******************

.. sectionauthor:: Kristian Rother, Patrick Yannul, Gavin Huttley


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

.. doctest::

    >>> from cogent.db.ncbi import EUtils
    >>> e = EUtils(numseqs=100, db='protein', rettype='gp')
    >>> result = e['"lysyl tRNA-synthetase"[ti] AND bacteria[orgn]']
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
    >>> print gb
    [('BAB52044', Sequence(MAGSNTI... 548)),...

For other
---------

.. OMIM, PUBMED, ??

Retrieving PubMed abstracts from NCBI via EUtils
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.db.ncbi import EUtils
    >>> e = EUtils(db='pubmed',rettype='brief')
    >>> result = e['Michalsky Preissner'].read()
    >>> print result
    <BLANKLINE>
    1: Hildebrand PW et al. SuperLooper--a prediction ser...[PMID: 19429894] 
    <BLANKLINE>
    2: Fest S et al. Survivin minigene DNA vaccina...[PMID: 19291796]...
    >>> e = EUtils(db='pubmed',rettype='abstract')
    >>> result = e['Michalsky Preissner'].read()
    >>> print result
    <BLANKLINE>
    1: Nucleic Acids Res. 2009 May 8; [Epub ahead of print] 
    <BLANKLINE>
    SuperLooper--a prediction server for the...

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

Orthologs
---------

Functional assignments
----------------------

Pathway assignments
-------------------

Ensembl
=======

Connecting
----------

.. Hosts and species

get genomic features
--------------------

Get alignments
--------------

Get SNPs
--------

PDB
===

for structures
--------------

Rfam
====

for rna secondary structures, alignments, functions
---------------------------------------------------

GoldenPath (not yet implemented)
================================

whole-genome alignments, orthologs, annotation tracks
-----------------------------------------------------

