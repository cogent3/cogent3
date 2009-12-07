Working with macromolecular structures
======================================

.. sectionauthor:: Marcin Cieslik

This part of the documentation presents examples on how to work with macromolecular structures i.e. coordinate files. This functionality has originates from ZenPDB.

At the current stage of development the input and output is limited to PDB files and some derivatives, but the "internal" representation and construction is file-format agnostic and other parsers can be easily added.

Hierarchy and identity
----------------------

A common way to describe macromolecular structures is the hierarchical representation. The hierarchy is made from entities (atoms, residues, chains, models, structures and crystals). In this hierarchical representation models are within structures (e.g. several NMR models of one protein structure or residues which are a collection of atoms). Hierarchical representations are unique i.e. each entity has to be uniquely identified for a given structure, which means that it has a unique identifier. We will refer to this identifier as the *full_id* in contrast with the *short_id* or just *id* which defines an entity uniquely only within it's parent. Each entity has only a single parent, but can have multiple children e.g. a residue is part of only one peptide chain, but it will contain multiple atoms.

An example of a *full_id*:

::

    # for an atom
    ('4TSV', 0, 'A', ('ARG', 131, ' '), ('O', ' '))
    # for a residue
    ('4TSV', 0, 'A', ('ARG', 131, ' '))
    # for a chain
    ('4TSV', 0, 'A')

The first is the identifier for the oxygen atom from the peptide bond of the ARG131 residue in the 'A' chain of the first model (0) in the structure available from the PDB as '4TSV'. A short version of an *id* which is specific only within a parent (i.e. for an atom within a residue) looks similar. In the example below we see the short id of the same atom. Of course this information is not enough to pin-point an atom in the structure but it is enough to identify different atoms within the same residue.

::

    (('O', ' '),)

As you can see the *full_id* is a linear tuple of short id's which can be either a tuple (e.g. ('O', ' ') for an oxygen atom) or a string (e.g. 'A' for chain A). All strings within a short id have some special meaning for example the id of a residue has the following structure ('three letter AA name', residue_id, 'insertion code'). It should be noted that according to the RCSB the ``residue_id`` is the integer number which should be 1 for the first natural residue in a protein. Residues can have negative ``residue_ids`` e.g. residues of a N-terminal affinity tag.

What is an Entity?
------------------

``Entity`` is the most basic class to provide methods specific to macromolecular structures. The ``Atom``, ``Residue``, ``Chain``, ``Model`` and ``Structure`` classes all inherit from the ``Entity`` class, yet there is some distinction between them. Only the ``Atom`` entity cannot contain other entities e.g. an instance of the ``Residue`` class can (and should) contain some ``Atom`` instances. The methods common to container entities are within the ``MultiEntity`` class which they all inherit from. The ``MultiEntity`` is also a subclass of the ``Entity`` class. It is important not to use the ``Entity`` and ``MultiEntity`` classes directly as some attributes (e.g. their position within the SMCRA hierarchy) have to be provided. In fact each entity is just a Python dictionary with almost all dictionary methods left untouched.

Parsing a macromolecular structure e.g. a PDB file means to create a ``Structure`` entity or in other words to recursively fill it with atoms, residues and chains.

Working with entities
---------------------

Our first task will be to parse and write a structure in a PDB file into an ``Entity`` hierarchy. This is quite easy and if you are familiar with the internals of PyCogent I hope it will also be obvious. You can use any PDB file the examples use the ``4TSV.pdb`` file in the doc/data directory.

The easy way, but implicit way:

.. doctest::

    >>> import cogent
    >>> structure = cogent.LoadStructure('data/4TSV.pdb')

This code involves quite a bit of magic, so let's do everything manually.

The ``cogent.LoadStructure`` method is a convenience method to get a structure object from a file in any of the supported formats. Right now PyCogent supports only the PDB file-format. Now let's read and write the same PDB file by using the ``PDBParser`` and ``PDBWriter`` function directly. The new_structure argument can be any ``Entity`` (e.g. a ``Structure`` entity) or a container of entities (e.g. a list of ``Atom`` and ``Residue`` entities):

.. doctest::

    >>> from cogent.parse.pdb import PDBParser
    >>> from cogent.format.pdb import PDBWriter
    >>> import tempfile, os
    >>> pdb_file = open('data/4TSV.pdb')
    >>> new_structure = PDBParser(pdb_file)
    >>> open_handle, file_name = tempfile.mkstemp()
    >>> os.close(open_handle)
    >>> new_pdb_file = open(file_name,'wb')
    >>> PDBWriter(new_pdb_file, new_structure)
    >>> new_structure
    <Structure id=4TSV>

In this code-listing we first import the PDB parser and PDB writer, open a PDB file and parse the structure. You can verify that the ``PDBParser`` does not close the open ``pdb_file``:

.. doctest::

    >>> assert not pdb_file.closed
    >>> assert not new_pdb_file.closed

Currently the ``PDBParser`` parses quite a lot information from the header of the PDB file and the atomic coordinates. It omits the anisotropic b-factors. Additional information is stored in the ``header`` attribute which is a dictionary.

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> structure.id # the static id tuple.
    ('4TSV',)
    >>> structure.getId() # the dynamic id tuple, use calls to get_id whenever possible.
    ('4TSV',)
    >>> structure.getFull_id() # only for the structure entity is the full_id identical to the id.
    ('4TSV',)
    >>> structure.header.keys() # the pdb header is parsed in to a dictionary as the header attribute
    ['bio_cmx', 'uc_mxs', 'name', 'solvent_content', 'expdta', 'bio_mxs',...
    >>> structure.header['id'] # this is the 4-char PDB ID parsed from the header and used to construct the structure.id
    '4TSV'
    >>> structure.header['expdta'] # if this is 'X-RAY' we probably deal with a x-ray structre and thus a lot crystallografic data is store in the header.
    'X-RAY'

Not all information from the PDB header is currently parsed, If you are interested in some special data you can access the unparsed header through the ``raw_header`` attribute, the same is true for the trailer. If you manage to extract the data from the ``raw_header`` you are ready to modify the modular code of the ``PDBParser`` class, please submit a patch!

::

    structure.raw_header
    structure.raw_trailer

The structure entity is a container for model entities, as you already know the structure is just a dictionary of models.

.. doctest::

    >>> structure.items()
    [((0,), <Model id=0>)]
    >>> structure.values()
    [<Model id=0>]
    >>> structure.keys()
    [(0,)]
    >>> first_model = structure.values()[0] # we name the first(and only) model in the structure
    >>> first_model_id = first_model.getId()

But PyCogent provides more specific methods to work with entities. The one which is useful to access the contents of an entity is ``getChildren``. The optional argument to the ``getChildren`` methods is a list of ids (e.g. to access only a subset of children) more concise and sophisticated methods to work with children will be introduced later

.. doctest::

    >>> structure.getChildren() # the output should be the same as structure.values()
    [<Model id=0>]
    >>> children_list = structure.getChildren([first_model_id])

A typical way to change a property of all children in a MultiEntity would be to write a loop. In this example we change the name of every residue to 'UNK'.

.. doctest::

   >>> some_model = structure.values()[0]
   >>> some_chain = some_model.values()[0]
   >>> for residue in some_chain.values():
   ...     residue.setName('UNK')
   ...

Pycogent allows to make it much shorter. Whenever a structure is created the top-level entity(i.e. the structure) gets pointer list to all the entities it contains stored as the ``table`` attribute. For example the structure entity will have a table with a list of all models, chains, residues and atoms that it contains. The keys of this table are *full_ids* the values the actual entities. The table is divided into sections based on the hierarchy i.e. there is a separate dictionary for residues, atoms, chains and models.

.. doctest::

   >>> sorted(structure.table.keys()) # all the different entity levels in the table (which is a normal dictionary)
   ['A', 'C', 'M', 'R']
   >>> structure.table['C'] # this is a full_id to entity mapping for all chains inside the structures
   {('4TSV', 0, ' '): <Chain id= >, ('4TSV', 0, 'A'): <Chain id=A>}

The creation of such a table is quite expensive so it is created for the structure entity, but there is no reason why you should not create a table for e.g. a chain if you need it.

.. doctest::

   >>> some_model = structure.values()[0]
   >>> some_chain = some_model.values()[0]
   >>> some_chain.setTable()
   >>> # some_chain.table['R'] # all the residues

There is however a catch. Tables are not dynamic, this means that they are not magically updated whenever a child changes it's id. This can be easily seen in following example where a new chain is created a residue moved into it. A table is created for the chain, but it does not update the key after the child changes it's name.

.. doctest::

    >>> from cogent.core.entity import Chain # the chain entity
    >>> new_chain = Chain('J') # an ampty chain named 'J'
    >>> new_chain.getId()
    ('J',)
    >>> some_residue = structure.table['R'].values()[0] # a semi-random residue from structure
    >>> # a possible output: <Residue UNK resseq=39 icode= >
    >>> some_residue.setName('001') # change the name to '001'
    >>> # some_residue.getId() # should return e.g. (('001', 39, ' '),)
    >>> # some_residue.getFull_id() # should return ('4TSV', 0, 'A', ('001', 39, ' '))
    >>> new_chain.addChild(some_residue) # move from chain 'A' in 4TSV into chain 'J'
    >>> # new_chain.keys() # should return: [(('001', 39, ' '),)]
    >>> new_chain.setTable()
    >>> # new_chain.table['R'].keys() # should return: [('J', ('001', 39, ' '))]
    >>> some_residue.setName('002') # change the name to '002'
    >>> # new_chain.keys() # should return: [(('002', 39, ' '),)] # updated!
    >>> # new_chain.table['R'].keys() # should return [('J', ('001', 39, ' '))] not updated
    >>> new_chain.setTable(force =True) # update table
    >>> # new_chain.table['R'].keys() # should return [('J', ('002', 39, ' '))] updated

It is important to realize that Python dictionaries are not sorted so the order of two equal dictionaries is not the same. Each time a child is changed in a way that affects the parent e.g. a part of it's id changes the parent dictionary will be updated and the order might also. You should **never** assume that an entity has a particular order.

.. doctest::

   >>> some_residue = some_chain.values()[0]
   >>> old_id = some_residue.getId() # e.g. (('ILE', 154, ' '),)
   >>> some_residue.setName('VAL')
   >>> new_id = some_residue.getId() # e.g. (('VAL', 154, ' '),)
   >>> some_chain.getChildren([old_id]) # nothin... not valid anymore
   []
   >>> # some_chain.getChildren([new_id]) # e.g. [<Residue VAL resseq=154 icode= >]

But the the table of an entity is static and does not get updated.

.. doctest::

   >>> some_full_id = some_residue.getFull_id() # entities in tables are stored using their full ids!!
   >>> # some_chain.table['R'][some_full_id] # should raise a KeyError
   >>> some_chain.setTable() # we make a new table
   >>> some_chain.table['R'][some_full_id]
   <Residue VAL resseq=131 icode= >

It is important to note that the table is a simple dictionary and the entity specific methods like ``getChildren`` are not available. You can figure out whether the table is up-to-date (or at least I hope I managed to code it right) using the ``modified`` attribute.

.. doctest::

   >>> some_chain.modified
   False

If the result were ``True`` the residue has been modified and might require to use the ``setTable``, or in some cases ``updateIds()`` methods.

.. doctest::

   >>> some_chain.setTable()
   >>> some_chain.updateIds()

Do not run those methods if you do not need to as they take some time.

The loop to run a child method can be implicitly omitted by using the dispatch method. It calls the method for every child.

.. doctest::

   >>> some_model = structure.values()[0]
   >>> some_chain = some_model.values()[1]
   >>> some_chain.dispatch('setName', 'UNK')
   >>> some_chain.modified
   True

The above method has exactly the same effect as the loop. All residues within the chain will have the name set to 'UNK'. You can verify that the id's and dictionary keys got updated.:

.. code-block:: python

    some_chain.keys()[0] # output random e.g. (('UNK', 260, ' '),)
    some_chain.values()[0] # e.g. <Residue UNK resseq=260 icode= >
