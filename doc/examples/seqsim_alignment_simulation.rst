Seqsim Simple Alignment Simulation Example
==========================================

.. sectionauthor:: Julia Goodrich

This is a very simple example of how to use the ``seqsim`` module to simulate
an alignment for a tree starting with a random sequence and substitution rate 
matrix (q). The rate matrix gives the rate constant of going from one character 
in the sequence to another character in the sequence, the Q matrix determines 
the rate of change of the sequence.

First we will perform the necessary imports:

* ``Rates`` is an object that stores the rate matrix data, it can also be used 
    to generate a random rate matrix given an ``Alphabet``.

* ``DnaUsage`` is a ``Usage`` object that stores the usage of each nucleotide.

* ``DnaPairs`` is an Alphabet it stores the DNA pairs (AA,AT,AC,...), it can
    be passed into the ``Rates`` object, defining the rate matrix pairs for DNA.

* ``DNA`` is a ``MolType`` object for DNA.

* ``RangeNode`` is the main ``seqsim`` Node object, it allows for the easy 
    evolution of sequences.

* ``DndParser`` is a parser for a newick format tree.

.. doctest::

    >>> from cogent.seqsim.usage import Rates, DnaUsage
    >>> from cogent.core.usage import DnaPairs
    >>> from cogent.core.moltype import DNA
    >>> from cogent.core.alignment import Alignment
    >>> from cogent.seqsim.tree import RangeNode
    >>> from cogent.parse.tree import DndParser
    
Now, lets specify a 4 taxon tree:

.. doctest::

    >>> t = DndParser('(a:0.4,b:0.3,(c:0.15,d:0.2)edge.0:0.1);', 
    ... constructor = RangeNode)
    
To generate a random DNA sequence, we first specify nucleotide frequencies 
with the ``DnaUsage`` object. Then we create a random DNA sequence that is 
five bases long.

.. doctest::

    >>> u = DnaUsage({'A':0.5,'T':0.2,'C':0.15,'G':0.25})
    >>> s = DNA.ModelSeq(u.randomIndices(5))
    >>> q = Rates.random(DnaPairs)

Set q at the base of the tree and propagate it to all nodes in the tree,

.. doctest::

    >>> t.Q = q
    >>> t.propagateAttr('Q')

Set a P matrix (probability matrix) from every Q matrix on each node: 
    P(t) = e^(Qt),

.. doctest::

    >>> t.assignP()

Use ``evolve`` to evolve sequences for each tip, Note: must evolve sequence
data, not sequence object itself (for speed)

.. doctest::
    >>> t.evolve(s._data)

Build alignment,

.. doctest::

    >>> seqs = {}
    >>> for n in t.tips():
    ...     seqs[n.Name] = DNA.ModelSeq(n.Sequence)
    >>> aln = Alignment(seqs)

The result is a Cogent ``Alignment`` object, which can be used the same way as
any other alignment object.

``evolveSeqs`` can be used instead of evolve to evolve multiple sequences
according to the same tree (can model either different genes, or different rate
categories within a gene that you then combine, etc...),

.. doctest::

    >>> from numpy import concatenate

First you need to use ``assignPs`` to assign the proper P matricies given rates:

.. doctest::

    >>> t.assignPs([.5, .75, 1])

There needs to be the same number of random sequences as there are rate 
catigories so we create a list of 3 random sequences,

.. doctest::

    >>> s = [DNA.ModelSeq(u.randomIndices(5))._data for i in range(0,3)]

Then use ``evolveSeqs`` to evolve a sequence for every tip with every rate.

.. doctest::

    >>> t.evolveSeqs(s)

Now to concatenate the sequences,

.. doctest::

    >>> seqs = {}
    >>> for n in t.tips():
    ...     for s in n.Sequences:
    ...         seqs[n.Name] = DNA.ModelSeq(concatenate(tuple(n.Sequences)))
    >>> aln = Alignment(seqs)

