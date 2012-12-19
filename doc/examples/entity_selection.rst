Selecting and grouping entities
===============================

.. sectionauthor:: Marcin Cieslik

The feature, which distinguishes PyCogent's approach to the handling of macromolecular structures is the flexible and concise way of selecting, grouping and retrieving data from entities. The concepts of entity and hierarchy are similar.

Overview of methods and functions covered.
------------------------------------------

The methods covered in this section of the manual deal with selecting entities for purposes like: "select all hydrogen atoms from chain B", "mask all hetero atoms", "remove all water molecules" etc. We start with the high-level functions first, which are concise and standard to low-level methods for fine grained manipulations

Selection based on hierarchy.
-----------------------------

Let's start by accessing a PDB file and creating a structure entity. We establish a connection to the PDB file server download a file and parse it.

.. doctest::

    >>> from cogent.parse.pdb import PDBParser
    >>> from cogent.db.pdb import Pdb
    >>> pdb = Pdb()
    >>> socket_handle = pdb['2E1F']
    >>> structure = PDBParser(socket_handle)

Let's see what we got

.. doctest::

    >>> print structure.header['name']
    HYDROLASE
    >>> print structure.header['experiment_type']
    X-RAY DIFFRACTION

WOW, thats descriptive. At least we know it is an X-Ray structure. Now how many chains does it have?

.. doctest::
    
    >>> structure[(0,)].getChildren()
    [<Chain id=A>]

We found the 'A' chain of the first (0-based indexing) model. We can dig deeper

.. doctest::

    >>> structure[(0,)][('A',)].sortedkeys()[0:2]
    [(('H_HOH', 1, ' '),), (('H_HOH', 2, ' '),)]

Only waters? Probably not. You can see what is inside a chain by looking inside the dictionary to get the list of short ids and child entities:

.. doctest::

    >>> chain_A = structure[(0,)][('A',)]
    >>> # chain_A.keys() # get the short_ids
    >>> # chain_A.values() # get the children
    >>> len(chain_A)
    147

This number is too high because we counted water molecules not only amino acids. But typing ``structure[(0,)][('A',)]`` is pretty boring and it requires to inspect the number of models and chain ids first. The function which allows to select entities from the hierarchy based on their identity is called ``einput``

.. doctest::

    >>> from cogent.struct.selection import einput
    >>> all_residues = einput(structure, 'R', 'my_residues')
    >>> all_atoms = einput(structure, 'A')
    >>> len(all_residues)
    147

Still waters are included.

Selection based on properties.
------------------------------

We already have a collection of entities ``all_residues`` which contains all residues in the structure regardless of the number of chains and models. Our task is to determine the number of non-water residues. The property which allows us to distinguish a water molecule from an amino acid is the name, which is stored as the ``name`` attribute.

.. doctest::

    >>> chain_A.name
    'A'
    >>> first_child = chain_A.sortedvalues()[0]
    >>> first_child.name
    'H_HOH'

We could write a loop to select those residues we can either loop over the residues in ``chain_A`` or ``all_residues`` as they are the same:

.. doctest::

    >>> non_water = []
    >>> for residue in chain_A:
    ...     if residue.name != 'H_HOH':
    ...          non_water.append(residue)
    ...
    >>> len(non_water)
    95

To make this more convenient each entity e.g. a ``Chain`` instance has a method to select children based on a property ``selectChildren``. The equivalent of the above expression is:

.. doctest::
    
    >>> non_water = chain_A.selectChildren('H_HOH', 'ne', 'name').values()

or

.. doctest::

    >>> non_water = all_residues.selectChildren('H_HOH', 'ne', 'name').values()
    >>> len(non_water)
    95

The first argument is a value, the second an operator name from the ``operator`` module, here 'ne' is for 'Not Equal'. The last argument 'name' is resolved by the ``data_children`` method which allows the user to retrieve data from a child entities attributes, xtra dictionary or methods. Here we get the data from the 'name' attribute. The ``selectChildren`` method returns a dictionary, where keys are the short ids and values are the child entities. The result can be put into a new entity holder.

.. doctest::

    >>> non_water_holder = einput(non_water, 'R')

But having to first group the entities via ``einput`` then select them only to put them into a new container seems awkward. It can be done in one step using the ``select`` function.

.. doctest::

    >>> from cogent.struct.selection import select
    >>> non_water_holder = select(structure, 'R', 'H_HOH', 'ne', 'name')
    >>> len(non_water_holder)
    95

Is there a serine(s) in the sequence?

.. doctest::

    >>> serines = select(structure, 'R', 'SER', 'eq', 'name')
    >>> serines.sortedkeys()[0]
    ('2E1F', 0, 'A', ('SER', 1146, ' '))

The function raises a ``ValueError`` if no entities can be selected.
