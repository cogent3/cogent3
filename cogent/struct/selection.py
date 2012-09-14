"""Contains functions to select and group structural entities."""

from cogent.core.entity import StructureHolder, ModelHolder, ChainHolder, \
                               ResidueHolder, AtomHolder, HIERARCHY


__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"


def select(entities, level, *args, **kwargs):
    """Shorthand for ``einput`` and subsequent ``selectChildren``. Returns
    
    Returns a ``Holder`` instance. The "name" can be specified.
    
    Additional arguments and keyworded arguments are passed to the 
    ``selectChildren`` method of the holder instance.
    """
    try:
        name = kwargs.pop('name')
    except KeyError:
        name = 'select'
    holder = einput(entities, level)
    selection = holder.selectChildren(*args, **kwargs)
    try:
        holder = einput(selection.values(), level, name)
    except ValueError:
        raise ValueError('No entities have been selected')
    return holder

def einput(entities, level, name=None):
    """Creates a ``XyzHolder`` instance of entities at the specified level. Where
    Xyz is 'Structure', 'Model', 'Chain', Residue' or 'Atom'.
    
    Arguments:
    
        - entities: ``Entity`` instance or sequence of entities.
        - level: one of 'H', 'S', 'M', 'C', 'R', 'A'
        - name: optional name of the ``XyzHolder`` instance.
    """
    # Keep it bug-free
    all = {}
    index = HIERARCHY.index(level)
    for entity in entities:                     # __iter__ override in Entity
        if index > HIERARCHY.index(entity.level):       # call for children
            all.update(get_children(entity, level))
        elif  index < HIERARCHY.index(entity.level):    # call for parents
            all.update(get_parent(entity, level))
        else:
            all.update({entity.getFull_id():entity})   # call for self
    higher_level = HIERARCHY[index - 1]                 # one up;)
    if all:
        name = name or higher_level
        if higher_level == 'C':
            holder = ResidueHolder(name, all)
        elif higher_level == 'R':
            holder = AtomHolder(name, all)
        elif higher_level == 'M':
            holder = ChainHolder(name, all)
        elif higher_level == 'S':
            holder = ModelHolder(name, all)
        elif higher_level == 'H':
            holder = StructureHolder(name, all)
    else:
        raise ValueError, "einput got no input entites."
    holder.setSort_tuple()
    return holder

def get_children(entity, level):
    """Return unique entities of lower or equal level
    
    Arguments:
    
        - entity: any ``Entity`` instance.
        - level: one of 'H', 'S', 'M', 'C', 'R', 'A'
    """
    entity.setTable()
    return entity.table[level]

def get_parent(entity, level):
    """Returns unique entities of higher level.
    
    Arguments:
    
        - entity: any ``Entity`` instance.
        - level: one of 'H', 'S', 'M', 'C', 'R', 'A'
    """
    parent = entity.getParent(level) # get the correct parent
    return {parent.getFull_id(): parent}
