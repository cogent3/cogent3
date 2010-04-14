***************
Structural data
***************

.. sectionauthor:: Kristian Rother, Patrick Yannul

Protein structures
==================

Reading Protein structures
--------------------------

How to retrieve a structure from the PDB database?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> from cogent.db.pdb import Pdb
    >>> p = Pdb()
    >>> # KR commented because this takes long.
    >>> #pdb_file = p['1cse']
    >>> #pdb = pdb_file.read()
    >>> #len(pdb)
    283986

This example will retrieve the structure as a PDB file string.

How to parse a PDB file?
^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> from cogent.parse.pdb import PDBParser
    >>> struc = PDBParser(open('data/4TSV.pdb'))
    >>> struc
    <Structure id=4TSV>


What happens if a PDB file does not exist
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the wrong filename is given to the PDBParser function, you will 
get an *IOError*.


How to parse a PDB entry directly from the web
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> from cogent.parse.pdb import PDBParser
    >>> # struc = PDBParser(p['1cse'])
    >>> # KR commented because this takes long.


How to access PDB header information?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> struc.header['id']
    '4TSV'

    >>> struc.header['resolution']
    '1.80'

    >>> struc.header['r_free']
    '0.262'

    >>> struc.header['space_group']
    'H 3'

Navigating structure objects
----------------------------

What does a structure object contain?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A **cogent.parse.pdb.Structure** object as returned by PDBParser
contains a tree-like hierarchy of **Entity** objects. They are organized 
such that **Structures** that contain **Models** that contain **Chains** 
that contain ''Residues** that in turn contain **Atoms**. 
You can read more about the entity model on
[URL of Marcins example page].


How to access a model from a structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To get the first model out of a structure:

.. doctest::
    
    >>> model = struc[(0,)]
    >>> model
    <Model id=0>

The key contains the model number as a tuple.


How to access a chain from a model?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To get a particular chain:

.. doctest::
    
    >>> chain = model[('A',)]
    >>> chain
    <Chain id=A>


How to access a residue from a chain?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To get a particular residue:

.. doctest::
    
    >>> resi = chain[('ILE', 154, ' '),]
    >>> resi
    <Residue ILE resseq=154 icode= >



What properties does a residue have?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> resi.res_id
    154

    >>> resi.name
    'ILE'

    >>> resi.h_flag
    ' '

    >>> resi.seg_id
    '    '

How to access an atom from a residue?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To get a particular atom:

.. doctest::
    
    >>> atom = resi[("N", ' '),]
    >>> atom
    <Atom ('N', ' ')>


What properties does an atom have?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> atom.name
    ' N  '

    >>> atom.element
    ' N'

    >>> atom.coords
    array([ 142.986,   36.523,    6.838])

    >>> atom.bfactor
    13.35

    >>> atom.occupancy
    1.0


What happens if a model/chain/residue/atom does not exist?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If any entity does not exist, you will get a **KeyError**.


Is there something special about heteroatoms to consider?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yes, they have the **h_flag** attribute set in residues.


How are Altlocs/insertion codes represented?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both are part of the residue/atom ID.


Useful methods to access Structure objects
------------------------------------------


How to access all atoms, residues etc via a dictionary 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The **table** property of a structure returns a two-dimensional
dictionary containing all atoms. The keys are 1) the entity level
(any of 'A','R','C','M') and
2) the combined IDs of Structure, Model, Chain, Residue, Atom
as a tuple.

.. doctest::
    
    >>> struc.table['A'][('4TSV', 0, 'A', ('HIS', 73, ' '), ('O', ' '))]
    <Atom ('O', ' ')>



How to calculate the center of mass of a model or chain
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> # NEEDS TO BE CHECKED WITH MARCIN
    >>> model.coords
    array([ 147.35930713,   35.30383834,   -3.48538525])

    >>> chain.coords
    array([ 145.42204284,   34.6970624 ,   -3.82628478])


How to get a list of all residues in a chain?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> chain.values()[0]
    <Residue ILE resseq=154 icode= >


How to get a list of all atoms in a chain?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> resi.values()[0]
    <Atom ('N', ' ')>


Constructing structures
-----------------------

How to create a new entity?
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Structure/Model/Chain/Residue/Atom objects can be created as follows:

.. doctest::
    
    >>> from cogent.core.entity import Structure,Model,Chain,Residue,Atom
    >>> from numpy import array
    >>> s = Structure('my_struc')
    >>> m = Model((0),)
    >>> c = Chain(('A'),)
    >>> r = Residue(('ALA', 1, ' ',),False,' ')
    >>> a = Atom(('C  ',' ',), 'C', 1, array([0.0,0.0,0.0]), 1.0, 0.0, 'C')


How to add entities to each other?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> s.addChild(m)
    >>> m.addChild(c)
    >>> c.addChild(r)
    >>> r.addChild(a)
    >>> s.table
    >>> # SHOULD NOT BE {'A': {}, 'C': {}, 'R': {}, 'M': {}}

How to remove a residue from a chain?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> c.delChild(r.id)
    >>> s.table
    {'A': {}, 'C': {}, 'R': {}, 'M': {}}


Geometrical analyses
--------------------

Calculating euclidean distances between atoms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> from cogent.maths.geometry import distance
    >>> atom1 = resi[('N', ' '),]
    >>> atom2 = resi[('CA', ' '),]
    >>> distance(atom1.coords, atom2.coords)
    1.4691967192993618


Calculating euclidean distances between coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> from numpy import array
    >>> from cogent.maths.geometry import distance
    >>> a1 = array([1.0, 2.0, 3.0])
    >>> a2 = array([1.0, 4.0, 9.0])
    >>> distance(a1,a2)
    6.324...


Calculating flat angles from atoms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.struct.dihedral import angle
    >>> atom3 = resi[('C', ' '),]
    >>> a12 = atom2.coords-atom1.coords
    >>> a23 = atom2.coords-atom3.coords
    >>> angle(a12,a23)
    1.8568181080150168


Calculates the angle in radians.


Calculating flat angles from coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.struct.dihedral import angle
    >>> a1 = array([0.0, 0.0, 1.0])
    >>> a2 = array([0.0, 0.0, 0.0])
    >>> a3 = array([0.0, 1.0, 0.0])    
    >>> a12 = a2-a1
    >>> a23 = a2-a3
    >>> angle(a12,a23)
    1.5707963267948966


Calculates the angle in radians.



Calculating dihedral angles from atoms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.struct.dihedral import dihedral
    >>> atom4 = resi[('CG1', ' '),]
    >>> dihedral(atom1.coords,atom2.coords,atom3.coords, atom4.coords)
    259.49277688244217

Calculates the torsion in degrees.


Calculating dihedral angles from coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> from cogent.struct.dihedral import dihedral
    >>> a1 = array([0.0, 0.0, 1.0])
    >>> a2 = array([0.0, 0.0, 0.0])
    >>> a3 = array([0.0, 1.0, 0.0])
    >>> a4 = array([1.0, 1.0, 0.0])
    >>> dihedral(a1,a2,a3,a4)
    90.0


Calculates the torsion in degrees.


Other stuff
-----------

How to count the atoms in a structure?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> len(struc.table['A'].values())
    1187

How to iterate over chains in canonical PDB order?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In PDB, the chain with space as ID comes last, the 
others in alphabetical order.

.. doctest::
    
    >>> for chain in model.sortedvalues(): print chain
    <Chain id=A>
    <Chain id= >


How to iterate over chains in alphabetical order?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want the chains in purely alphabetical order:

.. doctest::
    
    >>> # KR 2 ROB: Is this what you requested or is the above example enough?
    >>> keys = model.keys()
    >>> keys.sort()
    >>> for chain in [model[id] for id in keys]: print chain
    <Chain id= >
    <Chain id=A>


How to iterate over all residues in a chain?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> residues = [resi for resi in chain.values()]
    >>> len(residues)
    148


How to remove all water molecules from a structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> water = [r for r in struc.table['R'].values() if r.name=='H_HOH']
    >>> for resi in water: resi.parent.delChild(resi.id)
    >>> # KR: I THINK struc.table IS NOT UPDATED AUTOMATICALLY
    >>> len(struc.table['A'].values())
    1117




