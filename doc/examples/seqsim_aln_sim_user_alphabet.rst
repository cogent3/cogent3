Seqsim Alignment Simulation Example with Non-standard alphabet
==============================================================

.. sectionauthor:: Julia Goodrich

This is an example of how to use PyCogent's ``seqsim`` module to simulate an 
alignment where the alphabet is defined by the user for a simple tree starting 
with a random sequence and a random substitution rate matrix. 

First we will perform the necessary imports.


.. doctest::

    >>> from cogent.seqsim.usage import Rates
    >>> from cogent.core.alignment import DenseAlignment
    >>> from cogent.seqsim.tree import RangeNode
    >>> from cogent.parse.tree import DndParser
    >>> from cogent.core.alphabet import CharAlphabet
    >>> from cogent.seqsim.usage import Usage
    >>> from cogent.core.sequence import ModelSequence

Now, lets specify a 4 taxon tree:

.. doctest::

    >>> t = DndParser('(a:0.4,b:0.3,(c:0.15,d:0.2)edge.0:0.1);', 
    ... constructor = RangeNode)

Create the alphabet by passing in the characters to ``CharAlphabet`` then create
tuples of all the possible pairs using ** operator

.. doctest::

    >>> Bases = CharAlphabet('ABCD')
    >>> Pairs = Bases**2

Generate a random sequence with the new alphabet and a random rate matrix,
``Usage`` is being used to define character frequencies for the random
sequence. Then we create a random sequence of length five.

.. doctest::

    >>> u = Usage({'A':0.5,'B':0.2,'C':0.15,'D':0.25}, Alphabet = Bases)    
    >>> s = ModelSequence(u.randomIndices(5))
    >>> q = Rates.random(Pairs)

Set q at the base of the tree and propagate it to all nodes in the tree,

.. doctest::

    >>> t.Q = q
    >>> t.propagateAttr('Q')

Set a P matrix from every Q matrix on each node,

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
    ...     seqs[n.Name] = ModelSequence(n.Sequence,Bases)
    >>> aln = DenseAlignment(seqs,Alphabet=Bases)

The result is a Cogent ``Alignment`` object, which can be used the same way as
any other alignment object.

