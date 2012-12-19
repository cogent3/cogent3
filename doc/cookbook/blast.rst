.. _blast-usage:

*****************
Controlling BLAST
*****************

.. authors, Gavin Huttley, Tom Elliott, Jeremy Widmann

Preliminaries
-------------

In order to run BLAST locally (from a program running on your computer) you will need to do three things:

- Download the BLAST "executables" from NCBI
- Make sure these programs are available on your ``PATH``
- Construct and format a database to search against
- Tested on version 2.2.13

NCBI has recently changed the BLAST programs, and as yet PyCogent does not support the new versions. The "legacy" versions are available from `here <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_ (login as guest).


Detailed installation instructions are beyond the scope of this example, but are available at `NCBI's website <http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/unix_setup.html>`_ .
After downloading the programs and setting up your ``PATH``, you should test BLAST by doing this from the command line:

::
    
    $ blastall --help

Which should give this:

::
    
    blastall 2.2.22   arguments:
    
      -p  Program Name [String]
      -d  Database [String]
        default = nr...

The file ``refseqs.fasta`` contains some short sequences for use in the following examples.

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

BLAST with XML output
---------------------

In this example, we load a DNA sequence from a file in the data directory and BLAST against our formatted database as above.

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

BLAST with protein sequences
----------------------------

In this example, we load a protein sequence from a file in the data directory and BLAST against a new protein database we will create.  Since we want to BLAST protein sequences instead of DNA, we will have to construct a new BLAST database.

The file ``refseqs_protein.fasta`` contains some short sequences for use in the following examples.

.. doctest::
    
    >>> from cogent.app.formatdb import build_blast_db_from_fasta_path
    >>> result = build_blast_db_from_fasta_path('data/refseqs_protein.fasta', is_protein=True)
    >>> result[0]
    'data/refseqs_protein.fasta'...

Notice that we set the parameter ``is_protein`` to ``True`` since our database consists of protein sequences this time.  This was not necessary in the previous example, because ``is_protein`` is set to ``False`` by default.

Now that we have built our protein BLAST database, we can load our sequence and BLAST against this database.

.. doctest::

    >>> from cogent import LoadSeqs, PROTEIN
    >>> from cogent.app.blast import blast_seqs, Blastall
    >>> seqs = LoadSeqs('data/inseqs_protein.fasta', moltype=PROTEIN, aligned=False)
    >>> seq = seqs.getSeq('1091044_fragment')
	>>> seq
	ProteinSequence(IPLDFDK... 26)
	
Notice we need to use '-p':'blastp' in the parameters dictionary, since ``blastp`` is used for protein.

.. doctest::
	
	>>> params={'-p':'blastp','-m':'9'}
	>>> result = blast_seqs([seq], 
	...    Blastall, 
	...    blast_db = 'data/refseqs_protein.fasta',
	...    params = params)
	>>> data = result['StdOut'].read()
	>>> print data.split('\n')[:1]
	['# BLASTP 2.2...

We save the results for further processing 

.. doctest::
    
    >>> outfile = open('data/blast_protein_test.txt','w')
    >>> outfile.write(data)
    >>> outfile.close()

Now we will explore some of the convenience methods of the ``BlastResult`` object.

.. doctest::

	>>> from cogent.parse.blast import BlastResult
	>>> blast_results = BlastResult(open('data/blast_protein_test.txt','r'))

Suppose we want to filter our results based on various criteria.  In many cases you may want to only keep the top '3' matches with the longest 'ALIGNMENT LENGTH' for the query sequence to the target.

.. doctest::

	>>> best_hits = dict(blast_results.bestHitsByQuery(field='ALIGNMENT LENGTH', n=3))
	>>> query_1_best_hits = best_hits['1']
	>>> for hit in query_1_best_hits:
	...     for key, value in hit.items():
	...             print key.ljust(20), value
	...     print
	... 
	MISMATCHES           0
	ALIGNMENT LENGTH     26
	Q. END               26
	BIT SCORE            56.2
	% IDENTITY           100.00
	Q. START             1
	S. START             30
	S. END               55
	GAP OPENINGS         0
	QUERY ID             1
	E-VALUE              5e-12
	SUBJECT ID           1091044
	<BLANKLINE>
	MISMATCHES           10
	ALIGNMENT LENGTH     27
	Q. END               25
	BIT SCORE            33.5
	% IDENTITY           55.56
	Q. START             1
	S. START             32
	S. END               58
	GAP OPENINGS         1
	QUERY ID             1
	E-VALUE              3e-05
	SUBJECT ID           5326864
	<BLANKLINE>
	MISMATCHES           16
	ALIGNMENT LENGTH     24
	Q. END               25
	BIT SCORE            22.3
	% IDENTITY           33.33
	Q. START             2
	S. START             19
	S. END               42
	GAP OPENINGS         0
	QUERY ID             1
	E-VALUE              0.077
	SUBJECT ID           14286173
	<BLANKLINE>

The fist of the top 3 hits for alignment length has 0 MISMATCHES and a % IDENTITY of 100.00.  The next 2 hits have many MISMATCHES and a much lower % IDENTITY.  Lets filter the results again, but by E-VALUE this time:

.. doctest::

	>>> best_hits = dict(blast_results.bestHitsByQuery(field='E-VALUE', n=3))
	>>> query_1_best_hits = best_hits['1']
	>>> for hit in query_1_best_hits:
	...     for key, value in hit.items():
	...             print key.ljust(20), value
	...     print
	... 
	MISMATCHES           0
	ALIGNMENT LENGTH     26
	Q. END               26
	BIT SCORE            56.2
	% IDENTITY           100.00
	Q. START             1
	S. START             30
	S. END               55
	GAP OPENINGS         0
	QUERY ID             1
	E-VALUE              5e-12
	SUBJECT ID           1091044
	<BLANKLINE>
	MISMATCHES           10
	ALIGNMENT LENGTH     27
	Q. END               25
	BIT SCORE            33.5
	% IDENTITY           55.56
	Q. START             1
	S. START             32
	S. END               58
	GAP OPENINGS         1
	QUERY ID             1
	E-VALUE              3e-05
	SUBJECT ID           5326864
	<BLANKLINE>
	MISMATCHES           6
	ALIGNMENT LENGTH     18
	Q. END               26
	BIT SCORE            30.4
	% IDENTITY           66.67
	Q. START             9
	S. START             31
	S. END               48
	GAP OPENINGS         0
	QUERY ID             1
	E-VALUE              3e-04
	SUBJECT ID           15964668
	<BLANKLINE>

You can filter the BLAST results by any of the fields you like.  You can also use the ``BlastResult`` object to do a quick assessment of your BLAST results looking only at the fields you like:

.. doctest::

	>>> fields = ['SUBJECT ID', 'BIT SCORE', 'E-VALUE']
	>>> for query, results in blast_results.items():
	...     print ''.join([f.ljust(20) for f in fields])
	...     for result in results[-1]:
	...             print ''.join(map(str,[result[field].ljust(20) for field in fields]))
	SUBJECT ID          BIT SCORE           E-VALUE             
	1091044             56.2                5e-12               
	5326864             33.5                3e-05               
	15964668            30.4                3e-04               
	17229033            29.6                5e-04               
	21112072            28.1                0.001               
	4704732             25.8                0.007               
	13541117            24.6                0.016               
	15826629            24.3                0.020               
	14286173            22.3                0.077               
	6323138             21.9                0.10                
	18313548            20.8                0.22                
	21674812            20.0                0.38                
	14600438            20.0                0.38                
	4996210             18.5                1.1                 
	15605963            17.3                2.5                 
	15615431            16.5                4.2                 

.. doctest::
	    :hide:

	    >>> from cogent.util.misc import remove_files
	    >>> remove_files(['data/blast_protein_test.txt'],
	    ...              error_on_missing=False)
