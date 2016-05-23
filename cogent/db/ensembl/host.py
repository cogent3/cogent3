#!/usr/bin/env python
import warnings
import sqlalchemy as sql

try:
    import mysql.connector as mysql_connect
    connect_template = 'mysql+mysqlconnector://%(account)s/%(db_name)s?raise_on_warnings=False'
    password_arg = 'password'
    sql_version = tuple([int(v) for v in sql.__version__.split(".") if v.isdigit()])
    if sql_version < (0, 9, 7):
        warnings.warn('mysql.connector requires sqlalchemy >= 0.9.7\n')
        raise ImportError
except ImportError:
    import pymysql as mysql_connect
    connect_template = 'mysql+pymysql://%(account)s/%(db_name)s'
    password_arg = 'passwd'
except ImportError:
    import MySQLdb as mysql_connect
    connect_template = 'mysql+mysqldb://%(account)s/%(db_name)s'
    password_arg = 'passwd'

from cogent.util.table import Table
from cogent.db.ensembl.species import Species
from cogent.db.ensembl.name import EnsemblDbName
from cogent.db.ensembl.util import asserted_one

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Jason Merkin"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

class HostAccount(object):
    """host account data"""
    def __init__(self, host, user, passwd, port=None):
        super(HostAccount, self).__init__()
        self.host = host
        self.user = user
        self.passwd = passwd
        self.port = port or 3306
        self._hash = hash((self.host, self.user, self.port))
    
    def __cmp__(self, other):
        return cmp(self._hash, other._hash)
    
    def __eq__(self, other):
        return self._hash == other._hash
    
    def __hash__(self):
        return self._hash
    
    def __str__(self):
        return '%s:%s@%s:%s' % (self.user,self.passwd,self.host,self.port)
    

def get_ensembl_account(release=None):
    """returns an HostAccount for ensembl.
    
    Arguments:
        - release: if not specified, returns for the ensembl MySQL server
          hosting releases from 48"""
    port = [None, 5306][release is None or int(release) > 47]
    return HostAccount("ensembldb.ensembl.org", "anonymous", "", port=port)

def _get_default_connection():
    return "ensembldb.ensembl.org", "anonymous", ""

class EngineCache(object):
    """storage of active connections, indexed by account, database name"""
    _db_account = {}
    def __call__(self, account, db_name=None, pool_recycle=None):
        """returns an active SQLAlchemy connection engine"""
        assert account and db_name,"Must provide an account and a db"
        pool_recycle = pool_recycle or 3600
        if account not in self._db_account.get(db_name, []):
            if db_name == "PARENT":
                args = {password_arg: account.passwd}
                engine = mysql_connect.connect(host=account.host, user=account.user,
                                    port=account.port, **args)
            else:
                engine = sql.create_engine(connect_template % dict(account=account, 
                                        db_name=db_name), pool_recycle=pool_recycle)
            if db_name not in self._db_account:
                self._db_account[db_name] = {}
            self._db_account[db_name][account] = engine
        return self._db_account[db_name][account]

DbConnection = EngineCache()

def make_db_name_pattern(species=None, db_type=None, release=None):
    """returns a pattern for matching the db name against"""
    sep = r"%"
    pattern = ""
    if species:
        species = Species.getEnsemblDbPrefix(species)
        pattern = "%s%s" % (sep, species)
    if db_type:
        pattern = "%s%s%s" % (pattern, sep, db_type)
    if release:
        pattern = "%s%s%s" % (pattern, sep, release)
    assert pattern
    
    return "'%s%s'" % (pattern, sep)

def get_db_name(account=None, species=None, db_type=None, release=None,
                division=None, DEBUG=False):
    """returns the listing of valid data-base names as EnsemblDbName objects"""
    if account is None:
        account = get_ensembl_account(release=release)
    
    if DEBUG:
        print "Connection To:", account
        print "Selecting For:", species, db_type, release
    
    server = DbConnection(account, db_name='PARENT')
    cursor = server.cursor()
    show = "SHOW DATABASES"
    if species or db_type or release:
        pattern = make_db_name_pattern(species, db_type, release)
        show = "%s LIKE %s" % (show, pattern)
    if DEBUG:
        print show
    cursor.execute(show)
    rows = cursor.fetchall()
    dbs = []
    for row in rows:
        try:
            if division is not None and division not in row[0]:
                continue
            name = EnsemblDbName(row[0])
            if (release is None or name.Release == str(release)) and\
                                (db_type is None or name.Type == db_type):
                dbs.append(name)
        except (IndexError, RuntimeError):
            if DEBUG:
                print "FAIL:", row[0]
            continue
    return dbs

def get_latest_release(account = None):
    """returns the number of the latest release based on the compara db"""
    names = get_db_name(account=account, db_type="compara")
    compara = []
    for name in names:
        compara += [int(name.Release)]
    return str(max(compara))

if __name__ == "__main__":
    eaccount = get_ensembl_account(release='48')
    print get_db_name(account=eaccount, release="48", db_type='compara')
    print get_latest_release()
