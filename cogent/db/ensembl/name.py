import re

from species import Species

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

_release = re.compile("\d+")
def get_version_from_name(name):
    """returns the release and build identifiers from an ensembl db_name"""
    name = _name_delim.split(name)
    # if last entry has a release and build
    release, build = None, None
    if _release.findall(name[-2]):
        release = _release.findall(name[-2])[0]
        build = name[-1]
    else:
        release = _release.findall(name[-1])[0]
    return release, build

_name_delim = re.compile("_")
def get_dbtype_from_name(name):
    """returns the data base type from the name"""
    try:
        name = _name_delim.split(name)
    except TypeError, msg:
        print "Error:"
        print name, type(name), msg
        raise
    dbtype = None
    if name[0] == "ensembl":
        dbtype = name[1]
    else:
        dbtype = "_".join(name[2:-2])
    return dbtype

def get_db_prefix(name):
    """returns the db prefix, typically an organism or `ensembl'"""
    name = _name_delim.split(name)
    if name[0] == "ensembl":
        prefix = "ensembl"
    elif len(name) > 4:
        prefix = "_".join(name[:2])
    else:
        raise RuntimeError("Unknown name structure: %s" % "_".join(name))
    return prefix

class EnsemblDbName(object):
    """container for a db name, inferring different attributes from the name,
    such as species, version, build"""
    def __init__(self, db_name):
        """db_name: and Emsembl database name"""
        if isinstance(db_name, EnsemblDbName):
            db_name = db_name.Name
        self.Name = db_name
        self.Type = get_dbtype_from_name(db_name)
        self.Prefix = get_db_prefix(db_name)
        release, build = get_version_from_name(db_name)
        self.Release = release
        self.Build = build
        self.Species = None
        self.Species = Species.getSpeciesName(self.Prefix)
    
    def __repr__(self):
        build = ['', "; Build='%s'" % self.Build][self.Build != None]
        s = "db(Prefix='%s'; Type='%s'; Release='%s'%s)" % (self.Prefix, self.Type,
                    self.Release, build)
        return s
    
    def __str__(self):
        return self.Name
    
    def __cmp__(self, other):
        if isinstance(other, type(self)):
            other = other.Name
        return cmp(self.Name, other)
    
    def __hash__(self):
        return hash(self.Name)
