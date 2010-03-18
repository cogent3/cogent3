*******************
Accessing databases
*******************

.. Gavin Huttley, Kristian Rother, Patrick Yannul

NCBI
====

EUtils is a web service offered by the NCBI to access the sequence, literature and other databases by a special format of URLs. PyCogent offers an interface to construct the URLs and retrieve the results in text format.

For sequence
------------

Fetching FASTA or Genpept sequences from NCBI using EFetch with GI's
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.db.ncbi import EFetch
    >>> ef = EFetch(id='459567', rettype='fasta')
    >>> lines = ef.read().splitlines()
    >>> for line in lines:
    ...     print line[:40]
    ... 
    >gi|459567|dbj|D28543.1|HPCNS5PC Hepatit
    GAGCACGACATCTACCAATGTTGCCAACTGAACCCAGAGG
    GGCTTTACCTTGGTGGTCCCATGTTTAACTCGCGAGGTCA
    CGGGGTTCTTCCAACCAGCATGGGCAATACCCTCACATGT
    GCAGGCCTCACCAATTCTGACATGTTGGTTTGCGGAGATG
    TC
    <BLANKLINE>

.. TODO need commentary on the below

.. doctest::

    >>> genpept = EFetch(id='1234567,459567', rettype='gp').read()

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


Parsing in more detail:  a single GenBank entry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. TODO you could select these from each sequence using the getFeaturesMatching

.. doctest::

    >>> from cogent.db.ncbi import EUtils
    >>> from cogent.parse.genbank import RichGenbankParser
    >>> e = EUtils(db="nucleotide", rettype="gb")
    >>> record = e['154102'].readlines()
    >>> parser = RichGenbankParser(record)
    >>> accession, seq = [record for record in parser][0]
    >>> accession
    'STYHEMAPRF'
    >>> type(seq)
    <class 'cogent.core.sequence.DnaSequence'>
    >>> def gene_and_cds(f):
    ...     return f['type'] == 'CDS' and 'gene' in f
    ... 
    >>> cds_features = [f for f in seq.Info.features if gene_and_cds(f)]
    >>> for cds_feature in cds_features:
    ...     print cds_feature['gene'], cds_feature['location']
    ... 
    ['hemA'] 732..1988
    ['prfA'] 2029..3111

Retrieving a bacterial genome file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To obtain a full bacterial genome, run the following to get the complete *Salmonella typhimurium* genome sequence (Genbank) file. (For this documentation, we include a partial file for illustration purposes.)

.. code-block:: python
    
    from cogent.db.ncbi import EUtils
    e = EUtils(db="nucleotide", rettype="gb")
    outfile = open('data/ST.genome.gb','w')
    outfile.write(e['AE006468'].read())
    outfile.close()

Now do the analysis

.. doctest::
    
    >>> from cogent.parse.genbank import RichGenbankParser
    >>> infile = open('data/ST_genome_part.gb', 'r')
    >>> parser = RichGenbankParser(infile)
    >>> accession, seq = [record for record in parser][0]
    >>> gene_and_cds = lambda f: f['type'] == 'CDS' and 'gene' in f
    >>> gene_name = lambda f: f['gene']
    >>> all_cds = [f for f in seq.Info.features if gene_and_cds(f)]
    >>> for cds in sorted(all_cds, key=gene_name):
    ...     print cds['gene'][0].ljust(6),
    ...     print cds['protein_id'], cds['location']
    ... 
    mog    ['AAL18972.1'] 8729..9319
    talB   ['AAL18971.1'] 7665..8618
    thrA   ['AAL18966.1'] 337..2799
    thrB   ['AAL18967.1'] 2801..3730
    thrC   ['AAL18968.1'] 3734..5020
    thrL   ['AAL18965.1'] 190..255
    yaaA   ['AAL18969.1'] complement(5114..5887)
    yaaH   ['AAL18973.1'] complement(9376..9942)
    yaaJ   ['AAL18970.1'] complement(5966..7396)

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
    1: Chipman P et al. No association between the se...[PMID: 19997044] 
    <BLANKLINE>
    2: Yap VB et al. Estimates of the effect of na...[PMID: 19815689] 
    <BLANKLINE>
    3: Cherbuin N et al. Risk factors of transition fr...[PMID: 19628940]...
    >>> e = EUtils(db='pubmed',rettype='abstract')
    >>> result = e['19815689'].read()
    >>> print result
    <BLANKLINE>
    1: Mol Biol Evol. 2010 Mar;27(3):726-34. Epub 2009 Oct 8. 
    <BLANKLINE>
    Estimates of the effect of natural selection on protein-coding content.
    <BLANKLINE>
    Yap VB, Lindsay H, Easteal S, Huttley G.
    <BLANKLINE>
    Department of Statistics and Applied Probability, National University of
    Singapore, Singapore, Singapore. stayapvb@nus.edu.sg
    <BLANKLINE>
    Analysis of natural selection is key to understanding many core biological
    processes, including the emergence of competition, cooperation, and complexity...

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

.. following cleans up files

.. doctest::
    :hide:
    
    >>> from cogent.util.misc import remove_files
    >>> remove_files('ST.genome.gb', error_on_missing=False)
