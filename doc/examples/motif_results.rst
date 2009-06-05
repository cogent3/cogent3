Motif results example
=====================

.. sectionauthor:: Jeremy Widmann

In this example we will be parsing a motif results file and doing some basic operations highlighting the features of the various core motif handling objects. We first want to import the necessary modules.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> from cogent.parse.meme import MemeParser

Now we want to parse the MEME (http://meme.sdsc.edu) motif results file and the fasta file that we passed to MEME originally. This will construct a ``MotifResults`` object and a ``SequenceCollection`` object respectively, then add the ``SequenceCollection`` to the ``MotifResults``.

.. doctest::

    >>> results = MemeParser(open('data/motif_example_meme_results.txt','U'))
    >>> seqs = LoadSeqs('data/motif_example.fasta',aligned=False)
    >>> results.Alignment = seqs

Lets quickly look at an overview of the ``MotifResults``. First, we can check the ``MolType`` of the sequences, how many sequences were searched, and how many distinct motifs were found.

.. doctest::

    >>> results.MolType.label
    'protein'
    >>> len(results.Alignment.NamedSeqs)
    96
    >>> len(results.Motifs)
    10

Here 10 unique motifs were found searching 96 protein sequences. Now lets look in more detail at the motifs that were found by MEME. The ``MotifResults`` object has a list of ``Motif`` objects. Each ``Motif`` object contains a list of ``Module`` objects that make up the motif. In this example, each ``Motif`` has only one ``Module``. Show the module ID, Evalue, and number of instances of the module in the set of sequences.

.. doctest::

    >>> for motif in results.Motifs:
    ...     module = motif.Modules[0]
    ...     print module.ID, module.Evalue, len(module.NamedSeqs)
    ... 
    1 0.0 50
    2 7.3e-239 41
    3 2.2e-254 45
    4 7.8e-153 37
    5 2.9e-120 13
    6 2e-99 29
    7 4.2e-54 20
    8 6.1e-41 29
    9 3.5e-15 9
    10 5.3e-14 7

Module 1 has the smallest Evalue and the most instances, so we'll look at this one in more detail. The ``Module`` object is a subclass of the core ``Alignment`` object, so it shares much of this functionality. We can look at the consensus sequence for Module 1 calculated in different ways. Lets compare the consensus that MEME provides, the calculated majority consensus, and the calculated IUPAC consensus.

.. doctest::

    >>> module_1 = results.Motifs[0].Modules[0]
    >>> module_1.ID
    '1'
    >>> module_1.ConsensusSequence
    'GKPVVVDFWATWCGPCRxEAPILEELAKE'
    >>> module_1.majorityConsensus(transform=str)
    'GKPVVVDFWATWCGPCRAEAPILEELAKE'
    >>> module_1.IUPACConsensus()
    'XXXXXXXXXXXXCXXCXXXXXXXXXXXXX'

Here we can see that the consensus sequences provided by MEME and the calculated majority consensus are about the same. The IUPAC consensus is an easy way to see if any positions in the module are absolutely conserved. To get a better idea of the conservation of the module, we can calculate the uncertainties for every position in the module.

.. doctest::

    >>> iupac = module_1.IUPACConsensus()
    >>> majority = module_1.majorityConsensus()
    >>> uncertainty = module_1.uncertainties()
    >>> for i,m,u in zip(iupac,majority,uncertainty):
    ...     print i,m,u
    ... 
    X G 2.69585768303
    X K 2.29582593843
    X P 2.96578451217
    X V 1.61117952123
    X V 1.91067699662
    X V 2.01512726036
    X D 1.57769736083
    X F 0.777268500731
    X W 2.0045407601
    X A 0.522179190202
    X T 2.70369216641
    X W 0.282292189082
    C C 0.0
    X G 1.96072818839
    X P 0.937268500731
    C C 0.0
    X R 2.03875770182
    X A 3.68637013016
    X E 2.60359082041
    X A 2.9672863748
    X P 0.282292189082
    X I 3.49915032218
    X L 2.19664948376
    X E 2.71346937346
    X E 2.49058231553
    X L 1.94895812367
    X A 2.71230564207
    X K 2.85533775047
    X E 2.36191706121

The first column is the IUPAC consensus symbol, the second column is the majority consensus symbol, and the third column is the uncertainty at a given position in the module. The smaller the number, the less uncertainty, and therefore the more conserved the majority residue is at that position. Now that we have examined Module 1 in some detail, lets do some more simple tasks. How many different sequences is Module 1 located in?

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> module_1.LocationDict
    {'18309723': [284], '15614085': [58], '15966937': [59], 
    '18406743': [42, 362], '30021713': [75], '15988313': [59], 
    '16123427': [51], '15899007': [47], '15805225': [70], '16761507': [51],
    '7290567': [19], '16804867': [17], '4200327': [77], '27375582': [79],
    '3323237': [17], '17531233': [26], '267116': [18], '2822332': [32],
    '23098307': [51], '16759994': [58], '1651717': [32], '7109697': [16],
    '4155972': [15], '1174686': [20], '11499727': [15], '19746502': [90],
    '15599256': [24], '6687568': [18], '15597673': [155], '15615431': [70],
    '19705357': [37], '17537401': [25], '12044976': [17], '1633495': [19],
    '20092028': [5], '16078864': [58], '21222859': [101], '17547503': [61],
    '15805374': [76], '15614140': [61], '17229859': [20], '15218394': [21],
    '13358154': [19], '15605725': [59], '15791337': [0], '135765': [15],
    '140543': [42], '1388082': [27], '1729944': [16]}
    >>> len(module_1.LocationDict)
    49

The ``LocationDict`` property of the ``Module`` object is a dictionary of sequences IDs and indices in the sequence where the module was found. Here we see that Module 1 was found in 49 different sequences, which means that it was found twice in one sequence. We can find what other modules were found to have more than one instance in a given sequence.

.. doctest::

    >>> for motif in results.Motifs:
    ...     module = motif.Modules[0]
    ...     for seq_id, indices in module.LocationDict.items():
    ...             if len(indices) > 1:
    ...                     print module.ID, seq_id, indices
    ... 
    1 18406743 [42, 362]
    3 18406743 [104, 264, 424]

We see that Module 1 and Module 3 have more than one instance in sequence 18406743. Since this sequence is the only one to contain multiple instances of the same module, lets quickly examine some statistics of the alignment.

.. doctest::

    >>> len(results.Alignment.NamedSeqs['18406743'])
    578
    >>> lengths = [len(seq) for seq in results.Alignment.Seqs]
    >>> min(lengths)
    89
    >>> max(lengths)
    578
    >>> sum(lengths)/float(len(lengths))
    169.86458333333334

This sequence is the longest of all the sequences searched and more than 3 times longer than the average sequence.
