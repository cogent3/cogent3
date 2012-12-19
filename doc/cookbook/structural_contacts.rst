Intramolecular contacts
-----------------------

This section of the documentation explains how to use several PyCogent modules to find contacts between macromolecular entities e.g. residues or chains both within and outside the context of a crystal. It assumes that the following boilerplate code is entered into the interpreter. This includes the necessary imports and sets up a working set of ``Entity`` instances.

.. doctest::

    >>> from cogent.struct import contact
    >>> from cogent.parse.pdb import PDBParser
    >>> from cogent.struct import selection, manipulation, annotation
    >>> pdb_fh = open('data/1HQF.pdb')
    >>> pdb_structure = PDBParser(pdb_fh)
    >>> model = pdb_structure[(0,)]
    >>> chainA = model[('A',)]
    >>> chainB = model[('B',)]
    >>> chainC = model[('C',)]

Find contacts within a asymmetric unit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To find contacts and save the result in the xtra dictionary of an ``Entity`` instance we have to use the ``contact.contacts_xtra`` function. This function calls is a wrapper around a low level function to find the contacts and annotates the input with the result of this call, but also returns a dictionary of the annotations.

In this simple case we do not look for contacts in the crystal lattice, but for contacts between chains in the asymmetric unit, which corresponds to the default parameters of the ``contacts_xtra`` function. We seek for all contacts between all atoms in the asymmetric unit, which are at most 5.0A apart and belong to different chains.

.. doctest::

    >>> conts = contact.contacts_xtra(model, xtra_key='CNT', search_limit=5.0)

Let's explore what the output is all about. The keys in the dictionary are all the atoms, which have are involved in some contact(s). The length of the dictionary is not the total number of contacts.

.. doctest::

    >>> atom_id = ('1HQF', 0, 'B', ('VAL', 182, ' '), ('C', ' '))
    >>> atom_id_conts = conts[atom_id]

The value for this is a dictionary, where the contacts are store in the given ``xtra_key``.

.. doctest::

    >>> atom_id_conts
    {'CNT': {('1HQF', 0, 'A', ('GLY', 310, ' '), ('CA', ' ')): (4.5734119648245102, 0, 0)}}

The value is a dictionary of contacts with keys being ids of the involved atoms and values tuples defining the distance in A, symmetry operation id, and unit cell id. For contacts within the asymmetric unit there is no symmetry operation and no unit cell translation so both are ``0``.
