.. _blast-usage:

*****************
Controlling BLAST
*****************

.. authors, Gavin Huttley, Tom Elliott

Preliminaries
-------------

In order to run BLAST locally (from a program running on your computer) you will need to do three things:

- download the BLAST "executables" from NCBI
- make sure these programs are available on your ``PATH``
- construct and format a database to search against

NCBI has recently changed the BLAST programs, and as yet PyCogent does not support the new versions. The "legacy" versions are available from `here <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_ (login as guest).

Detailed instructions are beyond the scope of this example, but after downloading the programs and setting up your ``PATH``, you should test BLAST by doing this from the command line:

::
    
    $ blastall --help

Which should give this:

::
    
    blastall 2.2.22   arguments:
    
      -p  Program Name [String]
      -d  Database [String]
        default = nr...

The file ``refseqs.fasta`` contains some short sequences for use in the following examples. It is available from :download:`here <../data/refseqs.fasta>`.

.. TODO add to data_file_links.rst

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> seqs = LoadSeqs('data/refseqs.fasta',
    ...     moltype=DNA, aligned=False)
    >>> for seq in seqs.Seqs:
    ...     print seq
    ... 
    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC
    TGCAGCTTGAGCCACAGGAGAGAGCCTTC
    TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC
    ACCGATGAGATATTAGCACAGGGGAATTAGAACCA
    TGTCGAGAGTGAGATGAGATGAGAACA
    ACGTATTTTAATTTGGCATGGT
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    CCAGAGCGAGTGAGATAGACACCCAC

These sequences can be formatted as a database by running the ``formatdb`` program contained in the NCBI download (this assumes that ``formatdb`` is on your ``PATH``)

.. doctest::
    
    >>> from cogent.app.formatdb import build_blast_db_from_fasta_path
    >>> result = build_blast_db_from_fasta_path('data/refseqs.fasta')
    >>> result[0]
    'data/refseqs.fasta'...

The function ``build_blast_db_from_fasta_path()`` returns a tuple containing the path to the database, and a list of paths for all the new files written by ``formatdb``.

The final requirement is to know the path to the ``data`` directory that comes with your BLAST download.  This directory contains substitution matrices and other files.
 
BLAST with text output
----------------------

In this example, we load a DNA sequence from a file in the data directory and BLAST against our formatted database. The parameters dictionary contains two flagged arguments: -p for the program to use, and -m for the type of output we want. '-m':'9' requests "tabular with comment lines." See ``blastall --help`` for more details.

Also, the application controller is set up to require a path to the data directory even though we don't use a substitution matrix for DNA. Here we can just pass any string.

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> from cogent.app.blast import blast_seqs, Blastall
    >>> seqs = LoadSeqs('data/inseqs.fasta', moltype=DNA, aligned=False)
    >>> seq = seqs.getSeq('s2_like_seq')
    >>> seq
    DnaSequence(TGCAGCT... 28)
    >>> params={'-p':'blastn','-m':'9'}
    >>> result = blast_seqs([seq], 
    ...    Blastall, 
    ...    blast_db = 'data/refseqs.fasta',
    ...    blast_mat_root = 'xxx',
    ...    params = params)
    >>> data = result['StdOut'].read()
    >>> print data.split('\n')[:1]
    ['# BLASTN 2.2...

We save the results for further processing 

.. doctest::
    
    >>> outfile = open('data/blast_test.txt','w')
    >>> outfile.write(data)
    >>> outfile.close()

The simplest way to organize the results is to use a parser. The BLAST parsers operate on a file object.

.. doctest::

    >>> from cogent.parse.blast import MinimalBlastParser9
    >>> blastfile = open('data/blast_test.txt', 'r')
    >>> blast_results = MinimalBlastParser9(blastfile)
    >>> type(blast_results)
    <type 'generator'>
    >>> for result in blast_results:
    ...     print result
    ... 
    ({'QUERY': '1', 'FIELDS': 'Query id...

The results include one item for each query sequence. Each result consists of a tuple whose first item is a dictionary of metadata. The second item is a list of hits. For example, you might do this

.. doctest::

    >>> from cogent.parse.blast import MinimalBlastParser9
    >>> blastfile = open('data/blast_test.txt', 'r')
    >>> blast_results = MinimalBlastParser9(blastfile)
    >>> for result in blast_results:
    ...     meta_data, hit_list = result
    ...     fields = meta_data['FIELDS'].split(',')
    ...     for key, value in zip(fields, hit_list[0]):
    ...         print key.strip().ljust(20), value
    Query id             1
    Subject id           s2
    % identity           89.66
    alignment length     29
    mismatches           2
    gap openings         1
    q. start             1
    q. end               28
    s. start             1
    s. end               29
    e-value              6e-05
    bit score            26.3

NCBI recommends that you use XML as the output for BLAST. (They reserve the right to change the format for other output types). XML is the output when we pass '-m':'7'.

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> from cogent.app.blast import blast_seqs, Blastall
    >>> seqs = LoadSeqs('data/inseqs.fasta', moltype=DNA, aligned=False)
    >>> seq = seqs.getSeq('s2_like_seq')
    >>> params={'-p':'blastn','-m':'7'}
    >>> result = blast_seqs([seq], 
    ...    Blastall, 
    ...    blast_db = 'data/refseqs.fasta',
    ...    blast_mat_root = 'xxx',
    ...    params = params)
    >>> data = result['StdOut'].read()
    >>> outfile = open('data/blast_test.xml','w')
    >>> outfile.write(data)
    >>> outfile.close()

One nice thing about this format is that it includes the alignment. The organization of the results from this parser is slightly different. Each result is still a tuple, but the first item of the tuple is metadata about the BLAST settings (``meta_meta_data``). The keys for the fields in the output are contained as the first element in the list that is the second item in the result tuple.

.. doctest::

    >>> from cogent.parse.blast_xml import MinimalBlastParser7
    >>> blastfile = open('data/blast_test.xml', 'r')
    >>> blast_results = MinimalBlastParser7(blastfile)
    >>> for result in blast_results:
    ...     meta_meta_data, hit_list = result
    ...     key_list = hit_list[0]
    ...     for key, value in zip(key_list, hit_list[1]):
    ...         if 'ALIGN' in key:  
    ...             continue
    ...         print key.ljust(20), value
    QUERY ID             1
    SUBJECT_ID           lcl|s2
    HIT_DEF              No definition line found
    HIT_ACCESSION        s2
    HIT_LENGTH           29
    PERCENT_IDENTITY     26
    MISMATCHES           0
    GAP_OPENINGS         1
    QUERY_START          1
    QUERY_END            28
    SUBJECT_START        1
    SUBJECT_END          29
    E_VALUE              6.00825e-05
    BIT_SCORE            26.2635
    SCORE                13
    POSITIVE             26
    >>> from cogent.parse.blast_xml import MinimalBlastParser7
    >>> blastfile = open('data/blast_test.xml', 'r')
    >>> blast_results = MinimalBlastParser7(blastfile)
    >>> for result in blast_results:
    ...     meta_meta_data, hit_list = result
    ...     key_list = hit_list[0]
    ...     for key in ('QUERY_ALIGN','MIDLINE_ALIGN','SUBJECT_ALIGN'):
    ...         i = key_list.index(key)
    ...         print hit_list[1][i][:40]
    TGCAGCTTGAG-CACAGGTTAGAGCCTTC
    ||||||||||| ||||||  |||||||||
    TGCAGCTTGAGCCACAGGAGAGAGCCTTC

.. doctest::
    :hide:
    
    >>> from cogent.util.misc import remove_files
    >>> remove_files(['data/blast_test.txt', 'data/blast_test.xml'],
    ...              error_on_missing=False)

