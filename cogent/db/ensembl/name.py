import re

from species import Species

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

_release = re.compile("\d+")
def get_version_from_name(name):
    """returns the release and build identifiers from an ensembl db_name"""
    r = _release.search(name)
    if r is None:
        return None, None
    
    # first number run is release, followed by build
    # note, for the ensemblgenomes naming system, the second digit run is the
    # standard Ensembl release and the first is for the specified genome
    release = name[r.start(): r.end()]
    b = [s for s in _name_delim.split(name[r.end():]) if s]
    
    return release, b

_name_delim = re.compile("_")
def get_dbtype_from_name(name):
    """returns the data base type from the name"""
    try:
        name = _release.split(name)
        name = [s for s in _name_delim.split(name[0]) if s]
    except TypeError, msg:
        print "Error:"
        print name, type(name), msg
        raise
    dbtype = None
    if name[0] == "ensembl":
        dbtype = name[1]
    else:
        dbtype = "_".join(name[2:])
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
        self.GeneralRelease = self.Release
        
        if len(build) == 1:
            if self.Type != 'compara':
                self.Build = build[0]
            else:
                self.Build = None
                self.GeneralRelease = build[0]
        elif build:
            self.Build = build[1]
            self.GeneralRelease = build[0]
        else:
            self.Build  = None
        
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
