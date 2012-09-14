#!/usr/bin/env python
"""Provides Info, DbRef, DbRefs

Info is a dictionary and is the annotation object of a Sequence object. 
"""
from cogent.parse.record import MappedRecord
from cogent.util.misc import Delegator, FunctionWrapper, ConstrainedDict

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

class DbRef(object):
    """Holds a database accession, and optionally other data.

    Accession:      id in the database: str or int
    Db:             database name: str
    Name:           short name of the record: str
    Description:    description of the record, possibly lengthy: str
    Data:           any data associated with the record: arbitrary object

    str(DbRef) always returns the accession.
    """
    def __init__(self, Accession, Db='', Name='', Description='', \
        Data=None):
        """Returns new DbRef.

        str(DbRef) always returns the accession as a string.
        """
        self.Accession = Accession
        self.Db = Db
        self.Name = Name
        self.Description = Description
        self.Data = Data

    def __str__(self):
        """Returns accession."""
        return str(self.Accession)

    def __int__(self):
        """Tries to coerce accession to int."""
        return int(self.Accession)

    def __cmp__(self, other):
        """Compares by accession: tries numeric first, then alphabetic"""
        try:
            return cmp(int(self), int(other))
        except:
            return cmp(str(self), str(other))

def _make_list(obj):
    """Returns list corresponding to or containing obj, depending on type."""
    if isinstance(obj, list):
        return obj
    elif isinstance(obj, tuple):
        return list(obj)
    else:
        return [obj]

class DbRefs(MappedRecord, ConstrainedDict):
    """Holds Database -> [Accessions] mapping.
    
    The accessions for a particular database are always stored as a list.

    DbRefs will ultimately contain methods for actually getting the records
    from known databases.
    """
    ValueMask = FunctionWrapper(_make_list)
    DefaultValue = []

KnownDatabases = dict.fromkeys(['RefSeq', 'GenBank', 'GenNucl', 'GenPept', 
    'GI', 'SwissProt', 'PIR', 'EMBL', 'DDBJ',
    'NDB', 'PDB', 'Taxon', 'LocusLink', 'UniGene', 'OMIM', 'PubMed', 'COGS', 
    'CDD', 'Pfam', 'Rfam', 'GO', 'dbEST', 'IPI', 'rRNA', 'EC', 'HomoloGene', 
    'KEGG', 'BRENDA', 'EcoCyc', 'HumanCyc', 'BLOCKS'])

class Info(MappedRecord, Delegator):
    """Dictionary that stores attributes for Sequence objects.
    
    Delegates to DbRefs for database IDs.
    """
    Required = {'Refs':None}
    def __init__(self, *args, **kwargs):
        """Returns new Info object. Creates DbRefs if necessary."""
        temp = dict(*args, **kwargs)
        if 'Refs' in temp:
            refs = temp['Refs']
            if not isinstance(refs, DbRefs):
                refs = DbRefs(refs)
        else:
            refs = DbRefs()
        #move keys into refs if they belong there: allows init from flat dict 
        for key, val in temp.items():
            if key in KnownDatabases:
                refs[key] = val
                del temp[key]
                
        Delegator.__init__(self, refs)
        self['Refs'] = refs
        MappedRecord.__init__(self, temp)

    def __getattr__(self, attr):
        """Checks for attr in Refs first."""
        if attr in KnownDatabases:
            return getattr(self.Refs, attr)
        else:
            return super(Info, self).__getattr__(attr)

    def __setattr__(self, attr, val):
        """Try to set in Refs first."""
        if attr in KnownDatabases:
            return setattr(self.Refs, attr, val)
        else:
            return super(Info, self).__setattr__(attr, val)
            
    def __getitem__(self, item):
        """Checks for item in Refs first."""
        if item in KnownDatabases:
            return getattr(self.Refs, item)
        else:
            return super(Info, self).__getitem__(item)

    def __setitem__(self, item, val):
        """Try to set in Refs first."""
        if item in KnownDatabases:
            return setattr(self.Refs, item, val)
        else:
            return super(Info, self).__setitem__(item, val)
    
    def __contains__(self, item):
        """Checks for item in Refs first."""
        if item in KnownDatabases:
            return item in self.Refs
        else:
            return super(Info, self).__contains__(item)

