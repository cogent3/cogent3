************************
Structural Data Advanced
************************

This section covers more advanced structural entity handling tasks.

What happens if I try to add an entity (e.g atom) to a another entity (e.g. residue) if it is already there?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A PyCogent ``Entity`` is a subclass of a dictionary. Adding children is
essentially the same as updating a dictionary, with a minimal amount of
book-keeping. It is equivalent to the following:

::

    child.setParent(parent)
    child_id = child.getId()
    parent[child_id] = child
    parent.setModified(True, False)
 
This points the child entity to it's new parent (line 1) and adds the child to the
parent dictionary (line 3). The call to ``setModified`` notifies all parents
of the parent of the modification. A dictionary has unique keys and so a parent
has children with unique ids. If you try to add a child which has an id clash
it will update the parent and override the previous child, just like you would
update a dictionary.
 

Why are the short ids inside a tuple? 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Short ids are parts of a long id. The long id is a tuple. Short ids can be
concatenated to form a long id. This would not be possible if short ids were not
within a tuple initially. For example:
 
>>> (0,) + ('A',) + (('GLY', 209, ' '),) + (('C', ' '),)
(0, 'A', ('GLY', 209, ' '), ('C', ' '))
 
The output here is a valid long id of an atom for use in ``AtomHolder`` instances.


How do I select children of a ``MultiEntity`` instance by some feature.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Selection is a common task and ``PyCogent`` has a unified syntax for this via the 
``selectChildren`` method. The idea behind it is as follows:
 
    #. gather "requested data" from all children.
    #. compare each returned child value to the template "value" using the
       "operator"
    #. return children for which the comparison is ``True``
 
The signature of this method is selectChildren("value", "operator", "requested data").
In the first step all children return the "requested data", the request might be
an attribute, a value corresponding to a key in the ``parent.xtra`` dictionary or
any other query supported by the ``getData`` method.

.. doctest::

    >>> from cogent.parse.pdb import PDBParser
    >>> pdb_fh = open('data/1HQF.pdb')
    >>> pdb_structure = PDBParser(pdb_fh)
    >>> model = pdb_structure[(0,)]
    >>> chainA = model[('A',)]
 
Example 1: select all alanines from a chain.

.. doctest::

    >>> alanines = chainA.selectChildren('ALA', 'eq', 'name')
 
This requests the "name" attribute from all children in chain A and uses the
"eq" (equals) operator to compare this to "ALA". It returns a list of residues
which have this name.
 
Example 2: select all residues, which are not amino acids or nucleic acids.


.. doctest::

    >>> selection = chainA.selectChildren('H', 'eq', 'h_flag')
 
This requests the "h_flag" i.e. hetero-atom flag from all residues. For amino
acids and nucleic acids this should be "" for all other molecular entities "H",
so the function returns only ligands, waters etc.
 
Example 3: What if some children have data to return?
 
First we pick out a residue and modify it's xtra dictionary to contain some
custom data. We mark lys39 as a catalytic residue.
 
.. doctest::

    >>> lys39 = chainA[(('LYS', 39, ' '),)]
    >>> lys39.xtra['CATALYTIC'] = True
 
All other residues do not have a value corresponding to the "CATALYTIC" key. But
we still can select all "CATALYTIC" residues in chain A.
 

.. doctest::

    >>> catalytic = chainA.selectChildren(True, 'eq', 'CATALYTIC', xtra=True)
    >>> catalytic
    {(('LYS', 39, ' '),): <Residue LYS resseq=39 icode= >}
 
The difference is that we have requested a value from the "xtra" dictionary
instead of a hypothetical "CATALYTIC" attribute.


What comparison "operators" are supported for the ``selectChildren`` method?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The "operator" can be either a) a string corresponding to a function from the
``operator`` module from the python standard library. The list of currently
supported operators is: "gt", "ge", "lt", "le", "eq", "ne", "or_", "and_", 
"contains", "is_", "is_not" or alternatively it can be a a custom function,
which has the following signature operator(value, got), where "got" is the value
returned by the child and "value" is what it is compared to.


How can I copy, deep copy or serialize pickle an entity.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PyCogent ``MutltiEntity`` and ``Entity`` are Python objects and they
support the copy and deepcopy protocol.


.. doctest::

    >>> import cPickle
    >>> pickledA = cPickle.dumps(chainA)
    >>> unpickledA = cPickle.loads(pickledA)
    >>> unpickledA is chainA
    False
    >>> unpickledA == chainA
    True

In the above we have pickled and unpickled a ``MultiEntity`` instance.
This results in a new instance "unpickledA" which is the same as "chainA",
but has a different id (different objects, identity fails).

If you are only interested in obtaining a copy of an ``Entity`` instance
and not being able to share entities between python sessions. You can use
the functions from the ``copy`` module. Please note that copies and 
deep copies are the same

.. doctest::

    >>> from copy import copy, deepcopy
    >>> otherA = copy(chainA)
    >>> otherA is chainA
    False
    >>> otherA == chainA
    True
    >>> cys119 = chainA[(('CYS', 119, ' '),)]
    >>> cys119_other = otherA[(('CYS', 119, ' '),)]
    >>> cys119 is cys119_other
    False
    >>> cys119 == cys119_other
    True

