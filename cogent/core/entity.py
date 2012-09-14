#!/usr/bin/env python
"""Provides the entities, the building blocks of the SMRCA hierachy 
representation of a macromolecular structure.
    
The MultiEntity class is a special Entity class to hold multiple instances of
other entities. All Entities apart from the Atom can hold others and inherit
from the MultiEntity. The Entity is the most basic class to deal with 
structural and molecular data. Do not use it directly since some functions
depend on methods provided by sub-classes. Classes inheriting from MultiEntity 
have to provide some attributes during init e.g: self.level = a valid string 
inside the SMCRA hierarchy). Holders of entities are like normal MultiEntities,
but are temporary and are outside the parent-children axes.
"""

import cogent
from cogent.core.annotation import SimpleVariable
from numpy import (sqrt, arctan2, power, array, mean, sum)
from cogent.data.protein_properties import AA_NAMES, AA_ATOM_BACKBONE_ORDER, \
                                   AA_ATOM_REMOTE_ORDER, AREAIMOL_VDW_RADII, \
                                   DEFAULT_AREAIMOL_VDW_RADIUS, AA_NAMES_3to1
from cogent.data.ligand_properties import HOH_NAMES, LIGAND_AREAIMOL_VDW_RADII
from operator import itemgetter, gt, ge, lt, le, eq, ne, or_, and_, contains, \
                                                                 is_, is_not
from collections import defaultdict
from itertools import izip
from copy import copy, deepcopy

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

ALPHABET = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ_ '
HIERARCHY = ['H', 'S', 'M', 'C', 'R', 'A']
AREAIMOL_VDW_RADII.update(LIGAND_AREAIMOL_VDW_RADII)

# error while creating a structure (non-recoverable error)
class ConstructionError(Exception):
    """Cannot unambiguously create a structure."""
    pass

# warning while creating a structure 
# (something wrong with the input, but recoverable)
class ConstructionWarning(Exception):
    """Input violates some construction rules (contiguity)."""
    pass

def sort_id_list(id_list, sort_tuple):
    """Sorts lists of id tuples. The order is defined by the PDB file 
    specification."""
    (hol_loc, str_loc, mod_loc, chn_loc, res_loc, at_loc) = sort_tuple
    # even a simple id is a tuple, this makes sorting general
    def space_last(ch_id1, ch_id2):             # this is for chain sorting
        if ch_id1 == ' '  and ch_id2 != ' ':
            return 1
        if ch_id2 == ' '  and ch_id1 != ' ':
            return - 1
        if ch_id1 == ' '  and ch_id2 == ' ':
            return 0
        return cmp(ch_id1, ch_id2)

    def atom(at_id1, at_id2):
        # hydrogen atoms come last
        is_hydrogen1 = (at_id1[0] == 'H')
        is_hydrogen2 = (at_id2[0] == 'H')
        diff = cmp(is_hydrogen1, is_hydrogen2)

        # back bone come first
        if not diff:
            order1 = AA_ATOM_BACKBONE_ORDER.get(at_id1)
            order2 = AA_ATOM_BACKBONE_ORDER.get(at_id2)
            diff = cmp(order2, order1)

        # (B)eta, (D)elta, (G)amma, .... o(X)t
        if not diff:
            remote1 = AA_ATOM_REMOTE_ORDER.get(at_id1[1:2])
            remote2 = AA_ATOM_REMOTE_ORDER.get(at_id2[1:2])
            diff = cmp(remote1, remote2)

        # branching comes last
        if not diff:
            diff = cmp(at_id1[2:4], at_id2[2:4])
        return diff

        # SE vs CE - selenium first
        if not diff:
            alpha1 = ALPHABET.index(at_id1[0:1])
            alpha2 = ALPHABET.index(at_id2[0:1])
            diff = cmp(alpha2, alpha1)

    def residue(res_id1, res_id2):
        r1, r2 = 1, 1
        if res_id1 in AA_NAMES: r1 = 2
        if res_id1 in HOH_NAMES: r1 = 0
        if res_id2 in AA_NAMES: r2 = 2
        if res_id2 in HOH_NAMES: r2 = 0
        if r1 is r2:
            return cmp(res_id1, res_id2)
        else:
            return cmp(r2, r1)
    # this assumes that the implementation of sorting is stable.
    # does it work for others then cPython.
    if res_loc or res_loc is 0:
        id_list.sort(key=itemgetter(res_loc), cmp=lambda x, y: residue(x[0], y[0]))  # by res_name
    if at_loc or at_loc is 0:
        id_list.sort(key=itemgetter(at_loc), cmp=lambda x, y: space_last(x[1], y[1]))  # by alt_loc
    if at_loc or at_loc is 0:
        id_list.sort(key=itemgetter(at_loc), cmp=lambda x, y: atom(x[0], y[0]))  # by at_id
    if res_loc or res_loc is 0:
        id_list.sort(key=itemgetter(res_loc), cmp=lambda x, y: cmp(x[2], y[2]))  # by res_ic
    if res_loc or res_loc is 0:
        id_list.sort(key=itemgetter(res_loc), cmp=lambda x, y: cmp(x[1], y[1]))  # by res_id
    if chn_loc or chn_loc is 0:
        id_list.sort(key=itemgetter(chn_loc), cmp=space_last) # by chain
    if mod_loc or mod_loc is 0:
        id_list.sort(key=itemgetter(mod_loc))             # by model
    if str_loc or str_loc is 0:
        id_list.sort(key=itemgetter(str_loc))             # by structure
    return id_list

def merge(dicts):
    """Merges multiple dictionaries into a new one."""
    master_dict = {}
    for dict_ in dicts:
        master_dict.update(dict_)
    return master_dict

def unique(lists):
    """Merges multiple iterables into a unique sorted tuple (sorted set)."""
    master_set = set()
    for set_ in lists:
        master_set.update(set_)
    return tuple(sorted(master_set))


class Entity(dict):
    """Container object all entities inherit from it. Inherits from dict."""

    def __init__(self, id, name=None, *args):
        # This class has to be sub-classed!
        # the masked attribute has to be set before the __init__ of an Entity
        # because during __setstate__, __getstate__ sub-entities are iterated
        # by .values(), which relies on the attribute masked. to decide which
        # children should be omitted.
        self.masked = False
        self.parent = None                      # mandatory parent attribute
        self.modified = True                    # modified on creation

        self.id = (id,)                         # ids are non-zero lenght tuples
        self.name = (name or id)                # prefer name over duplicate id
        self.xtra = {}                          # mandatory xtra dict attribute
        # Dictionary that keeps additional properties
        dict.__init__(self, *args)              # finish init as dictionary

    def __copy__(self):
        return deepcopy(self)

    def __deepcopy__(self, memo):
        new_state = self.__getstate__()
        new_instance = self.__new__(type(self))
        new_instance.__setstate__(new_state)
        return new_instance

    def __getstate__(self):
        new_state = copy(self.__dict__)          # shallow
        new_state['parent'] = None
        return new_state

    def __setstate__(self, new_state):
        self.__dict__.update(new_state)

    def __repr__(self):
        """Default representation."""
        # mandatory getLevel from sub-class
        return "<Entity id=%s, level=%s>" % (self.getId(), self.getLevel())

    def __sub__(self, entity):
        """Override "-" as Euclidean distance between coordinates."""
        return sqrt(sum(pow(self.coords - entity.coords, 2)))

    def _setId(self, id):
        self.name = id[0]

    def _getId(self):
        return (self.name,)

    def getId(self):
        """Return the id."""
        return self._getId()

    def getFull_id(self):
        """Return the full id."""
        parent = self.getParent()
        if parent:
            full_id = parent.getFull_id()
        else:
            full_id = ()                    # we create a tuple on the top
        full_id = full_id + self.getId()   # merge tuples from the left
        return full_id

    def setId(self, id_=None):
        """Set the id. Calls the ``_setId`` method."""
        if (id_ and id_ != self.id) or (not id_ and (self.getId() != self.id)):
            self.id = (id_ or self.getId())
            self.setModified(True, True)
            self._setId(self.id)
            if self.parent:
                self.parent.updateIds()

    def _setMasked(self, masked, force=False):
        if masked != self.masked or force:
            self.masked = masked             # mask or unmask
            self.setModified(True, False)   # set parents as modified

    def setMasked(self, *args, **kwargs):
        """Set masked flag (``masked``) ``True``."""
        self._setMasked(True, *args, **kwargs)

    def setUnmasked(self, *args, **kwargs):
        """Set masked flag (``masked``) ``False``."""
        self._setMasked(False, *args, **kwargs)

    def setModified(self, up=True, down=False):
        """Set modified flag (``modified``) ``True``."""
        self.modified = True
        if up and self.parent:
            self.parent.setModified(True, False)

    def setUnmodified(self, up=False, down=False):
        """Set modified flag (``modified``) ``False``."""
        self.modified = False
        if up and self.parent:
            self.parent.setUnmodified(True, False)

    def setParent(self, entity):
        """Set the parent ``Entity`` and adds oneself as the child."""
        if self.parent != entity:
            # delete old parent
            self.delParent()
            # add new parent
            self.parent = entity
            self.parent.addChild(self)
            self.setModified(False, True)

    def delParent(self):
        """Detach mutually from the parent. Sets both child and parent modified 
        flags (``modified``) as ``True``."""
        if self.parent:
            self.parent.pop(self.getId())
            self.parent.setModified(True, False)
        self.parent = None
        self.setModified(False, True)

    def getModified(self):
        """Return value of the modified flag (``modified``)."""
        return self.modified

    def getMasked(self):
        """Return value of the masked flag (``masked``)."""
        return self.masked

    def setLevel(self, level):
        """Set level (``level``)."""
        self.level = level

    def getLevel(self):
        """Return level (``level``)in the hierarchy."""
        return self.level

    def setName(self, name):
        """Set name."""
        self.name = name
        self.setId()

    def getName(self):
        """Return name."""
        return self.name

    def getParent(self, level=None):
        """Return the parent ``Entity`` instance."""
        if not level:
            return self.parent
        elif level == self.level:
            return self
        return self.parent.getParent(level)

    def move(self, origin):
        """Subtract the origin coordinates from the coordintats (``coords``)."""
        self.coords = self.coords - origin

    def setCoords(self, coords):
        """Set the entity coordinates. Coordinates should be a 
        ``numpy.array``."""
        self.coords = coords

    def getCoords(self):
        """Get the entity coordinates."""
        return self.coords

    def getScoords(self):
        """Return spherical (r, theta, phi) coordinates."""
        x, y, z = self.coords
        x2, y2, z2 = power(self.coords, 2)
        scoords = array([sqrt(x2 + y2 + z2), \
                   arctan2(sqrt(x2 + y2), z), \
                   arctan2(y, x)])
        return scoords

    def getCcoords(self):
        """Return redundant, polar, clustering-coordinates on the unit-sphere. 
        This is only useful for clustering."""
        x, y, z = self.coords
        x2, y2, z2 = power(self.coords, 2)
        ccoords = array([arctan2(sqrt(y2 + z2), x), \
                          arctan2(sqrt(x2 + z2), y), \
                          arctan2(sqrt(x2 + y2), z)
                          ])
        return ccoords

    def setScoords(self):
        """Set ``entity.scoords``, see: getScoords."""
        self.scoords = self.getScoords()

    def setCcoords(self):
        """Set ``entity.ccoords``, see: getCcoords."""
        self.ccoords = self.getCcoords()


class MultiEntity(Entity):
    """The ``MultiEntity`` contains other ``Entity`` or ``MultiEntity`` 
    instances."""
    def __init__(self, long_id, short_id=None, *args):
        self.index = HIERARCHY.index(self.level)   # index corresponding to the hierarchy level
        self.table = dict([(level, {}) for level in HIERARCHY[self.index + 1:]])  # empty table
        Entity.__init__(self, long_id, short_id, *args)

    def __repr__(self):
        id_ = self.getId()
        return "<MultiEntity id=%s, holding=%s>" % (id_, len(self))

    def _link(self):
        """Recursively adds a parent pointer to children."""
        for child in self.itervalues(unmask=True):
            child.parent = self
            try:
                child._link()
            except AttributeError:
                pass

    def _unlink(self):
        """Recursively deletes the parent pointer from children."""
        for child in self.itervalues(unmask=True):
            child.parent = None
            try:
                child._unlink()
            except AttributeError:
                pass

    def __getstate__(self):
        new_dict = copy(self.__dict__)          # shallow copy
        new_dict['parent'] = None               # remove recursion
        new_children = []
        for child in self.itervalues(unmask=True):
            new_child_instance = deepcopy(child)
            new_children.append(new_child_instance)
        return (new_children, new_dict)

    def __setstate__(self, new_state):
        new_children, new_dict = new_state
        self.__dict__.update(new_dict)
        for child in new_children:
            self.addChild(child)

    def __copy__(self):
        return deepcopy(self)

    def __deepcopy__(self, memo):
        new_state = self.__getstate__()
        new_instance = self.__new__(type(self))
        new_instance.__setstate__(new_state)
        return new_instance

    def __iter__(self):
        return self.itervalues()

    def setSort_tuple(self, sort_tuple=None):
        """Set the ``sort_tuple attribute``. The ``sort_tuple`` is a tuple 
        needed by the ``sort_id_list`` function to correctly sort items within 
        entities."""
        if sort_tuple:
            self.sort_tuple = sort_tuple
        else:   # making the sort tuple, ugly, uggly, uaughhlly ble
            sort_tuple = [None, None, None, None, None, None]
            key_lenght = len(self.keys()[0])
            stop_i = self.index + 2                    # next level, open right [)
            start_i = stop_i - key_lenght               # before all nones
            indexes = range(start_i, stop_i)            # Nones to change
            for value, index in enumerate(indexes):
                sort_tuple[index] = value
            self.sort_tuple = sort_tuple

    def getSort_tuple(self):
        """Return the ``sort_tuple`` attribute. If not set calls the 
        ``setSort_tuple`` method first. See: ``setSort_tuple``."""
        if not hasattr(self, 'sort_tuple'):
            self.setSort_tuple()
        return self.sort_tuple

    def itervalues(self, unmask=False):
        return (v for v in super(MultiEntity, self).itervalues() if not v.masked or unmask)

    def iteritems(self, unmask=False):
        return ((k, v) for k, v in super(MultiEntity, self).iteritems() if not v.masked or unmask)

    def iterkeys(self, unmask=False):
        return (k for k, v in super(MultiEntity, self).iteritems() if not v.masked or unmask)

    def values(self, *args, **kwargs):
        return list(self.itervalues(*args, **kwargs))

    def items(self, *args, **kwargs):
        return list(self.iteritems(*args, **kwargs))

    def keys(self, *args, **kwargs):
        return list(self.iterkeys(*args, **kwargs))

    def __contains__(self, key, *args, **kwargs):
        return key in self.keys(*args, **kwargs)

    def sortedkeys(self, *args, **kwargs):
        list_ = sort_id_list(self.keys(*args, **kwargs), self.getSort_tuple())
        return list_

    def sortedvalues(self, *args, **kwargs):
        values = [self[i] for i in self.sortedkeys(*args, **kwargs)]
        return values

    def sorteditems(self, *args, **kwargs):
        items = [(i, self[i]) for i in self.sortedkeys()]
        return items

    def _setMasked(self, masked, force=False):
        """Set the masked flag (``masked``) recursively. If forced proceed even
        if the flag is already set correctly."""
        if masked != self.masked or force:   # the second condition is when
            if masked:                       # an entity has all children masked
                # we have to mask children        # but is not masked itself
                for child in self.itervalues():   # only unmasked children
                    child.setMasked()
                    child.setModified(False, False)
            else:
                # we have to unmask children
                for child in self.itervalues(unmask=True):
                    if child.masked or force:     # only masked children
                        child.setUnmasked(force=force)
                        child.setModified(False, False)
            self.masked = masked
            self.setModified(True, False)   # set parents as modified

    def setModified(self, up=True, down=True):
        """Set the modified flag (``modified``) ``True``. If down proceeds 
        recursively for all children. If up proceeds recursively for all 
        parents."""
        self.modified = True
        if up and self.parent:
            self.parent.setModified(True, False)
        if down:
            for child in self.itervalues(unmask=True):
                child.setModified(False, True)

    def setUnmodified(self, up=False, down=False):
        """Set the modified (``modified``) flag ``False``. If down proceeds 
        recursively for all children. If up proceeds recursively for all 
        parents."""
        self.modified = False
        if up and self.parent:
            self.parent.setUnmodified(True, False)
        if down:
            for child in self.itervalues(unmask=True):
                child.setUnmodified(False, True)

    def _initChild(self, child):
        """Initialize a child (during construction)."""
        child.parent = self
        self[child.getId()] = child

    def addChild(self, child):
        """Add a child."""
        child.setParent(self)
        child_id = child.getId()
        self[child_id] = child
        self.setModified(True, False)

    def delChild(self, child_id):
        """Remove a child."""
        child = self.get(child_id)
        if child:
            child.delParent()
            self.setModified(True, False)

    def getChildren(self, ids=None, **kwargs):
        """Return a copy of the list of children."""
        if ids:
            children = []
            for (id_, child) in self.iteritems(**kwargs):
                if id_ in ids:
                    children.append(child)
        else:
            children = self.values(**kwargs)
        return children

    def _setTable(self, entity):
        """Recursive helper method for ``entity.setTable``."""
        for e in entity.itervalues():
            self.table[e.getLevel()].update({e.getFull_id():e})
            self._setTable(e)

    def setTable(self, force=True, unmodify=True):
        """Populate the children table (``table``) recursively with all children
        grouped into hierarchy levels. If forced is ``True`` the table will be
        updated even if the ``Entity`` instance is not modified. If unmodify is 
        ``True`` the ``Entity`` modified flag (``modified``) will be set 
        ``False`` afterwards."""
        if self.modified or force:
            # a table is accurate as long as the contents of a dictionary do not
            # change.
            self.delTable()
            self._setTable(self)
        if unmodify:
            self.setUnmodified()

    def delTable(self):
        """Delete all children from the children-table (``table``). This does
        not modify the hierarchy."""
        self.table = dict([(level, {}) for level in HIERARCHY[self.index + 1:]])
        self.modified = True

    def getTable(self, level):
        """Return children of given level from the children-table 
        (``table``)."""
        return self.table[level]

    def updateIds(self):
        """Update self with children ids."""
        ids = []
        for (id_, child) in self.iteritems():
            new_id = child.getId()
            if id_ != new_id:
                ids.append((id_, new_id))
        for (old_id, new_id)  in ids:
            child = self.pop(old_id)
            self.update(((new_id, child),))

    def getData(self, attr, xtra=False, method=False, forgiving=True, sorted=False):
        """Get data from children attributes, methods and xtra dicts as a list. 
        If is ``True`` forgiving remove ``None`` values from the output. 
        ``Nones`` are place-holders if a child does not have the requested data.
        If xtra is True the xtra dictionary (``xtra``) will be searched, if
        method is ``True`` the child attribute will be called."""
        values = self.sortedvalues() if sorted else self.values()
        if xtra:
            # looking inside the xtra of children
            data = [child.xtra.get(attr) for child in values] # could get None
        else:
            # looking at attributes
            data = []
            for child in values:
                try:
                    if not method:
                        data.append(getattr(child, attr))
                    else:
                        data.append(getattr(child, attr)())
                except AttributeError: #
                    data.append(None)
        if forgiving: # remove Nones
            data = [point for point in data if point is not None]
        return data

    def propagateData(self, function, level, attr, **kwargs):
        """Propagate data from child level to this ``Entity`` instance. The
        function defines how children data should be transformed to become
        the parents data e.g. summed."""
        if self.index <= HIERARCHY.index(level) - 2:
            for child in self.itervalues():
                child.propagateData(function, level, attr, **kwargs)
        datas = self.getData(attr, **kwargs)
        if isinstance(function, basestring):
             function = eval(function)        
        transformed_datas = function(datas)
        if kwargs.get('xtra'):
            self.xtra[attr] = transformed_datas
        else:
            setattr(self, attr, transformed_datas)
        return transformed_datas

    def countChildren(self, *args, **kwargs):
        """Count children based on ``getData``. Additional arguments and 
        keyworded arguments are passed to the ``getData`` method."""
        data = self.getData(*args, **kwargs)
        children = defaultdict(int) # by default returns 0
        for d in data:
            children[d] += 1
        return children

    def freqChildren(self, *args, **kwargs):
        """Frequency of children based on ``countChildren``. Additional 
        arguments and keyworded arguments are passed to the ``countChildren`` 
        method."""
        children_count = self.countChildren(*args, **kwargs)
        lenght = float(len(self))  # it could be len(children_count)?
        for (key_, value_) in children_count.iteritems():
            children_count[key_] = value_ / lenght
        return children_count

    def splitChildren(self, *args, **kwargs):
        """Splits children into groups children based on ``getData``.
        Additional arguments and keyworded arguments are passed to the 
        ``getData`` method."""
        kwargs['forgiving'] = False
        data = self.getData(*args, **kwargs)
        clusters = defaultdict(dict) # by default returns {}
        for (key, (id_, child)) in izip(data, self.iteritems()):
            clusters[key].update({id_:child})
        return clusters

    def selectChildren(self, value, operator, *args, **kwargs):
        """Generic method to select children, based on ``getData``. 
        Returns a dictionary of children indexed by ids. Compares the data item
        for each child using the operator name e.g. "eq" and a value e.g. 
        "H_HOH". Additional arguments and keyworded arguments are passed to the 
        ``getData`` method."""
        kwargs['forgiving'] = False
        data = self.getData(*args, **kwargs)
        children = {}
        for (got, (id_, child)) in izip(data, self.iteritems()):
            if isinstance(operator, basestring):
                operator = eval(operator)
            if operator(value, got):
                children.update({id_:child})
        return children

    def ornamentChildren(self, *args, **kwargs):
        """Return a list of (ornament, (id, child)) tuples, based on 
        ``getData``. Useful for sorting see: Schwartzian transform. 
        Forgiving is set False. Additional arguments and keyworded arguments are
        passed to the ``getData`` method."""
        kwargs['forgiving'] = False
        data = self.getData(*args, **kwargs)
        children = []
        for (got, (id_, child)) in izip(data, self.iteritems()):
            children.append((got, (id_, child)))
        return children

    def ornamentdictChildren(self, *args, **kwargs):
        """Return a dictionary of ornaments indexed by child ids, based on 
        ``getData``. Forgiving is set False. Additional arguments and 
        keyworded arguments are passed to the ``getData`` method."""
        kwargs['forgiving'] = False
        data = self.getData(*args, **kwargs)
        propertydict = {}
        for (got, id_) in izip(data, self.iterkeys()):
            propertydict.update(((id_, got),))
        return propertydict

    def stripChildren(self, *args, **kwargs):
        """Strips children based on selection criteria. See: 
        ``selectChildren``. Additional arguments and keyworded arguments are 
        passed to the ``selectChildren`` method."""
        children_ids = self.selectChildren(*args, **kwargs).keys()
        for id_ in children_ids:
            self.delChild(id_)

    def maskChildren(self, *args, **kwargs):
        """Mask children based on selection criteria. See: ``selectChildren``.
        Additional arguments and keyworded arguments are passed to the 
        ``selectChildren`` method."""
        children = self.selectChildren(*args, **kwargs).itervalues()
        for child in children:
            child.setMasked() # child.setModified child.parent.setModified

    def unmaskChildren(self, *args, **kwargs):
        """Unmask children based on selection criteria. See: 
        ``selectChildren``. Additional arguments and keyworded arguments are 
        passed to the ``selectChildren`` method."""
        children = self.selectChildren(*args, **kwargs).itervalues()
        for child in children:
            child.setUnmasked() # child.setModified child.parent.setModified

    def moveRecursively(self, origin):
        """Move ``Entity`` instance recursively to the origin."""
        for child in self.itervalues():
            try:
                child.moveRecursively(origin)
            except:
                # Atoms do not have this
                child.move(origin)
                pass
        self.setCoords()

    def setCoordsRecursively(self):
        """Set coordinates (``coords``) recursively. Useful if any child had its
        coordinates changed."""
        for child in self.itervalues():
            try:
                child.setCoordsRecursively()
            except:
                #Atoms do not have this
                pass
        self.setCoords()

    def setCoords(self, *args, **kwargs):
        """Set coordinates (``coords``) as a centroid of children coordinates. 
        A subset of children can be selected for the calculation. See: 
        ``Entity.selectChildren``. Additional arguments and keyworded arguments
        are passed to the ``getData`` method."""
        # select only some children
        if args or kwargs:
            children = self.selectChildren(*args, **kwargs).values()
        else:
            children = self
        coords = []
        for child in children:
            coords.append(child.getCoords())
        self.coords = mean(coords, axis=0)

    def getCoords(self):
        """Returns the current coordinates (``coords``). Raises an 
        ``AttributeError`` if not set."""
        try:
            return self.coords
        except AttributeError:
            raise AttributeError, "Entity has coordinates not set."

    def dispatch(self, method, *args, **kwargs):
        """Calls a method of all children with given arguments and keyworded 
        arguments."""
        for child in self.itervalues():
            getattr(child, method)(*args, **kwargs)


class Structure(MultiEntity):
    """The ``Structure`` instance contains ``Model`` instances."""
    def __init__(self, id, *args, **kwargs):
        self.level = 'S'
        MultiEntity.__init__(self, id, *args, **kwargs)

    def __repr__(self):
        return '<Structure id=%s>' % self.getId()

    def removeAltmodels(self):
        """Remove all models with an id != 0"""
        self.stripChildren((0,), 'ne', 'id', forgiving=False)

    def getDict(self):
        """See: ``Entity.getDict``."""
        return {'structure':self.getId()[0]}


class Model(MultiEntity):
    """The ``Model`` instance contains ``Chain`` instances."""
    def __init__(self, id, *args, **kwargs):
        self.level = 'M'
        MultiEntity.__init__(self, id, *args, **kwargs)

    def __repr__(self):
        return "<Model id=%s>" % self.getId()

    def getDict(self):
        """See: ``Entity.getDict``."""
        try:
            from_parent = self.parent.getDict()
        except AttributeError:
            # we are allowed to silence this becaus a structure id is not 
            # required to write a proper pdb line.            
            from_parent = {}
        from_parent.update({'model':self.getId()[0]})
        return from_parent


class Chain(MultiEntity):
    """The ``Chain`` instance contains ``Residue`` instances."""
    def __init__(self, id, *args, **kwargs):
        self.level = 'C'
        MultiEntity.__init__(self, id, *args, **kwargs)

    def __repr__(self):
        return "<Chain id=%s>" % self.getId()

    def removeHetero(self):
        """Remove residues with the hetero flag."""
        self.stripChildren('H', 'eq', 'h_flag', forgiving=False)

    def removeWater(self):
        """Remove water residues."""
        self.stripChildren('H_HOH', 'eq', 'name', forgiving=False)

    def residueCount(self):
        """Count residues based on ``name``."""
        return self.countChildren('name')

    def residueFreq(self):
        """Calculate residue frequency (based on ``name``)."""
        return self.freqChildren('name')

    def getSeq(self, moltype ='PROTEIN'):
        """Returns a Sequence object from the ordered residues.
        The "seq_type" determines allowed residue names."""
        if moltype == 'PROTEIN':
            valid_names = AA_NAMES
            moltype = cogent.PROTEIN
        elif moltype == 'DNA':
            raise NotImplementedError('The sequence type: %s is not implemented' % moltype)
        elif moltype == 'RNA':
            raise NotImplementedError('The sequence type: %s is not implemented' % moltype)
        else:
            raise ValueError('The \'seq_type\' is not supported.')
        aa = ResidueHolder('aa', self.selectChildren(valid_names, contains, 'name'))
        aa_noic = ResidueHolder('noic', aa.selectChildren(' ', eq, 'res_ic'))
        
        raw_seq = []
        full_ids = []
        for res in aa_noic.sortedvalues():
            raw_seq.append(AA_NAMES_3to1[res.name])
            full_ids.append(res.getFull_id()[1:]) 
        raw_seq = "".join(raw_seq)
        
        seq = cogent.Sequence(moltype, raw_seq, self.getName())
        seq.addAnnotation(SimpleVariable, 'entity_id', 'S_id', full_ids)
        return seq

    def getDict(self):
        """See: ``Entity.getDict``."""
        from_parent = self.parent.getDict()
        from_parent.update({'chain_id':self.getId()[0]})
        return from_parent


class Residue(MultiEntity):
    """The ``Residue`` instance contains ``Atom`` instances."""
    def __init__(self, res_long_id, h_flag, seg_id, *args, **kwargs):
        self.level = 'R'
        self.seg_id = seg_id
        self.h_flag = h_flag
        self.res_id = res_long_id[1] #ID number
        self.res_ic = res_long_id[2] #ID long   NAME
        MultiEntity.__init__(self, res_long_id, res_long_id[0], *args, **kwargs)

    def __repr__(self):
        res_name, res_id, res_ic = self.getId()[0]
        full_name = (res_name, res_id, res_ic)
        return "<Residue %s resseq=%s icode=%s>" % full_name

    def _getId(self):
        """Return the residue full id. ``(name, res_id, res_ic)``."""
        return ((self.name, self.res_id, self.res_ic),)

    def _setId(self, id):
        """Set the residue id ``res_id``, name ``name`` and insertion code 
        ``res_ic`` from a full id."""
        (self.name, self.res_id, self.res_ic) = id[0]

    def removeHydrogens(self):
        """Remove hydrogen atoms."""
        self.stripChildren(' H', 'eq', 'element', forgiving=False)

    def getSeg_id(self):
        """Return the segment id."""
        return self.seg_id

    def setSeg_id(self, seg_id):
        """Set the segment id. This does not change the id."""
        self.seg_id = seg_id

    def getIc(self):
        """Return the insertion code."""
        return self.res_ic

    def setIc(self, res_ic):
        """Set the insertion code."""
        self.res_ic = res_ic
        self.setId()

    def getRes_id(self):
        """Get the id."""
        return self.res_id

    def setRes_id(self, res_id):
        """Set the id."""
        self.res_id = res_id
        self.setId()

    def getH_flag(self):
        """Return the hetero flag."""
        return self.h_flag

    def setH_flag(self, h_flag):
        """Sets the hetero flag. A valid flag is ' ' or 'H'. If 'H' the flag
        becomes part of the residue name i.e. H_XXX."""
        if not h_flag in (' ', 'H'):
            raise AttributeError, "Only ' ' and 'H' hetero flags allowed."
        if len(self.name) == 3:
            self.name = "%s_%s" % (h_flag, self.name)
        elif len(self.name) == 5:
            self.name = "%s_%s" % (h_flag, self.name[2:])
        else:
            raise ValueError, 'Non-standard residue name'
        self.h_flag = h_flag
        self.setId()

    def getDict(self):
        """See: ``Entity.getDict``."""
        from_parent = self.parent.getDict()
        if self.h_flag != ' ':
            at_type = 'HETATM'
        else:
            at_type = 'ATOM  '
        from_parent.update({'at_type': at_type,
                            'h_flag': self.h_flag,
                            'res_name': self.name,
                            'res_long_id': self.getId()[0],
                            'res_id': self.res_id,
                            'res_ic': self.res_ic,
                            'seg_id': self.seg_id, })
        return from_parent


class Atom(Entity):
    """The ``Atom`` class contains no children."""
    def __init__(self, at_long_id, at_name, ser_num, coords, occupancy, bfactor, element):
        self.level = 'A'
        self.index = HIERARCHY.index(self.level)
        self.coords = coords
        self.bfactor = bfactor
        self.occupancy = occupancy
        self.ser_num = ser_num
        self.at_id = at_long_id[0]
        self.alt_loc = at_long_id[1]
        self.table = dict([(level, {}) for level in HIERARCHY[self.index + 1:]])
        self.element = element
        Entity.__init__(self, at_long_id, at_name)

    def __nonzero__(self):
        return bool(self.id)

    def __repr__(self):
        return "<Atom %s>" % self.getId()

    def _getId(self):
        """Return the full id. The id of an atom is not its ' XX ' name 
        but this string after left/right spaces striping. The full id is 
        ``(at_id, alt_loc)``."""
        return ((self.at_id, self.alt_loc),)

    def _setId(self, id):
        """Set the atom id ``at_id`` and alternate location ``alt_loc`` from a
        full id. See: ``_getId``."""
        (self.at_id, self.alt_loc) = id[0]

    def setElement(self, element):
        """Set the atom element ``element``."""
        self.element = element

    def setName(self, name):
        """Set name and update the id."""
        self.name = name
        self.setAt_id(name.strip())

    def setAt_id(self, at_id):
        """Set id. An atom id should be derived from the atom name. See:
        ``_getId``."""
        self.at_id = at_id
        self.setId()

    def setAlt_loc(self, alt_loc):
        """Set alternate location identifier."""
        self.alt_loc = alt_loc
        self.setId()

    def setSer_num(self, n):
        """Set serial number."""
        self.ser_num = n

    def setBfactor(self, bfactor):
        """Set B-factor."""
        self.bfactor = bfactor

    def setOccupancy(self, occupancy):
        """Set occupancy."""
        self.occupancy = occupancy

    def setRadius(self, radius=None, radius_type=AREAIMOL_VDW_RADII, \
                            default_radius=DEFAULT_AREAIMOL_VDW_RADIUS):
        """Set radius, defaults to the AreaIMol VdW radius."""
        if radius:
            self.radius = radius
        else:
            try:
                self.radius = radius_type[(self.parent.name, self.name)]
            except KeyError:
                self.radius = default_radius

    def getSer_num(self):
        """Return the serial number."""
        return self.ser_num

    def getBfactor(self):
        """Return the B-factor."""
        return self.bfactor

    def getOccupancy(self):
        """Return the occupancy."""
        return self.occupancy

    def getRadius(self):
        """Return the radius."""
        return self.radius

    def getDict(self):
        """See: ``Entity.getDict``."""
        from_parent = self.parent.getDict()
        from_parent.update({'at_name': self.name,
                            'ser_num': self.ser_num,
                            'coords': self.coords,
                            'occupancy': self.occupancy,
                            'bfactor': self.bfactor,
                            'alt_loc': self.alt_loc,
                            'at_long_id': self.getId()[0],
                            'at_id': self.at_id,
                            'element': self.element})
        return from_parent


class Holder(MultiEntity):
    """The ``Holder`` instance exists outside the SMCRA hierarchy. Elements in 
    a ``Holder`` instance are indexed by the full id."""
    def __init__(self, name, *args):
        if not hasattr(self, 'level'):
            self.level = name
        MultiEntity.__init__(self, name, name, *args)

    def __repr__(self):
        return '<Holder level=%s name=%s>' % (self.level, self.getName())

    def addChild(self, child):
        """Add a child."""
        child_id = child.getFull_id()
        self[child_id] = child

    def delChild(self, child_id):
        """Remove a child."""
        self.pop(child_id)

    def updateIds(self):
        """Update self with children long ids."""
        ids = []
        for (id_, child) in self.iteritems():
            new_id = child.getFull_id()
            if id_ != new_id:
                ids.append((id_, new_id))
        for (old_id, new_id)  in ids:
            child = self.pop(old_id)
            self.update(((new_id, child),))


class StructureHolder(Holder):
    """The ``StructureHolder`` contains ``Structure`` instances. See: 
    ``Holder``."""
    def __init__(self, *args):
        self.level = 'H'
        Holder.__init__(self, *args)

    def __repr__(self):
        return "<StructureHolder name=%s>" % self.getName()


class ModelHolder(Holder):
    """The ``ModelHolder`` contains ``Model`` instances. See: ``Holder``."""
    def __init__(self, *args):
        self.level = 'S'
        Holder.__init__(self, *args)

    def __repr__(self):
        return "<ModelHolder name=%s>" % self.getName()


class ChainHolder(Holder):
    """The ``ChainHolder`` contains ``Chain`` instances. See: ``Holder``."""
    def __init__(self, *args):
        self.level = 'M'
        Holder.__init__(self, *args)

    def __repr__(self):
        return "<ChainHolder name=%s>" % self.getName()


class ResidueHolder(Holder):
    """The ``ResidueHolder`` contains ``Residue`` instances. See: ``Holder``."""
    def __init__(self, *args):
        self.level = 'C'
        Holder.__init__(self, *args)

    def __repr__(self):
        return "<ResidueHolder name=%s>" % self.getName()


class AtomHolder(Holder):
    """The ``AtomHolder`` contains ``Atom`` instances. See: ``Holder``."""
    def __init__(self, *args):
        self.level = 'R'
        Holder.__init__(self, *args)

    def __repr__(self):
        return "<AtomHolder name=%s>" % self.getName()


class StructureBuilder(object):
    """Constructs a ``Structure`` object. The ``StructureBuilder`` class is used
    by a parser class to parse a file into a ``Structure`` object. An instance 
    of a ``StructureBuilder`` has methods to create ``Entity`` instances and add
    them into the SMCRA hierarchy``."""
    def __init__(self):
        self.structure = None

    def initStructure(self, structure_id):
        """Initialize a ``Structure`` instance."""
        self.structure = Structure(structure_id)

    def initModel(self, model_id):
        """Initialize a ``Model`` instance and add it as a child to the 
        ``Structure`` instance. If a model is defined twice a 
        ``ConstructionError`` is raised."""
        if not (model_id,) in self.structure:
            self.model = Model(model_id)
            self.model.junk = AtomHolder('junk')
            self.structure._initChild(self.model)
        else:
            raise ConstructionError

    def initChain(self, chain_id):
        """Initialize a ``Chain`` instance and add it as a child to the 
        ``Model`` instance. If a chain is defined twice a 
        ``ConstructionWarning`` is raised. This means that the model is not
        continuous."""
        if not (chain_id,) in self.model:
            self.chain = Chain(chain_id)
            self.model._initChild(self.chain)
        else:
            self.chain = self.model[(chain_id,)]
            raise ConstructionWarning, "Chain %s is not continous" % chain_id

    def initSeg(self, seg_id):
        """Does not create an ``Entity`` instance, but updates the segment id,
        ``seg_id`` which is used to initialize ``Residue`` instances."""
        self.seg_id = seg_id

    def initResidue(self, res_long_id, res_name):
        """Initialize a ``Residue`` instance and add it as a child to the 
        ``Chain`` instance. If a residue is defined twice a 
        ``ConstructionWarning`` is raised. This means that the chain is not
        continuous."""
        if not (res_long_id,) in self.chain:
            self.residue = Residue(res_long_id, res_name, self.seg_id)
            self.chain._initChild(self.residue)
        else:
            self.residue = self.chain[(res_long_id,)]
            raise ConstructionWarning, "Residue %s%s%s is not continuous" % \
                                                                    res_long_id

    def initAtom(self, at_long_id, at_name, ser_num, coord, occupancy, \
                  bfactor, element):
        """Initialize an ``Atom`` instance and add is as child to the 
        ``Residue`` instance. If an atom is defined twice a 
        ``ConstructionError`` is raised and the ``Atom`` instance is added to 
        the ``structure.model.junk`` ``Holder`` instance."""
        if not (at_long_id,) in self.residue:
            self.atom = Atom(at_long_id, at_name, ser_num, coord, occupancy, \
                             bfactor, element)
            self.residue._initChild(self.atom)
        else:
            full_id = (tuple(self.residue[(at_long_id,)].getFull_id()), \
                       ser_num)
            self.model.junk._initChild(Atom(full_id, at_name, ser_num, coord, \
                                             occupancy, bfactor, element))
            raise ConstructionError, 'Atom %s%s is defined twice.' % at_long_id

    def getStructure(self):
        """Update coordinates (``coords``), set the children-table (``table``) 
        and return the ``Structure`` instance."""
        self.structure.setTable()
        self.structure.setCoordsRecursively()
        return self.structure

