*********************
Dealing with HTS data
*********************

FASTQ formatted files
=====================

Parsing
-------

FASTQ format can be exported by Illumina's pipeline software.

.. doctest::
    
    >>> from cogent.parse.fastq import MinimalFastqParser
    >>> for label, seq, qual in MinimalFastqParser('data/fastq.txt'):
    ...     print label
    ...     print seq
    ...     print qual
    GAPC_0015:6:1:1259:10413#0/1
    AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
    ````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
    GAPC_0015:6:1:1283:11957#0/1
    TATGTATATATAACATATACATATATACATACATA
    ]KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb...


Converting quality scores to numeric data
-----------------------------------------

In FASTQ format, ASCII characters are used to represent base-call quality. Unfortunately, vendors differ in the range of characters used. According to their documentation, Illumina uses the character range from 64-104. We parse the sequence file and convert the characters into integers on the fly.

.. doctest::
    
    >>> from cogent.parse.fastq import MinimalFastqParser
    >>> for label, seq, qual in MinimalFastqParser('data/fastq.txt'):
    ...     qual = map(lambda x: ord(x)-64, qual)
    ...     print label
    ...     print seq
    ...     print qual
    GAPC_0015:6:1:1259:10413#0/1
    AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
    [32, 32, 32, 32, 25, 30, 20, 29, 32, 29, 35, 30, 35, 33, 34, 35, 33, ...
