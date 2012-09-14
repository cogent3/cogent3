"""Contains functions to annotate macromolecular entities."""

from cogent.core.entity import HIERARCHY

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"

def xtradata(data, entity):
    """Annotates an entity with data from a ``{full_id:data}`` dictionary. The
    ``data`` should also be a dictionary.
    
    Arguments:
        - data: a dictionary, which is a mapping of full_id's (keys) and data
                dictionaries.
        - entity: top-level entity, which contains the entities which will hold 
                  the data."""
    for full_id, data in data.iteritems():
        sub_entity = entity
        strip_full_id = [i for i in full_id if i is not None]
        for short_id in strip_full_id:
            sub_entity = sub_entity[(short_id,)]
        sub_entity.xtra.update(data)

