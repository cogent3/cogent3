import sqlalchemy as sql

from cogent.util import table as cogent_table
from cogent.db.ensembl.host import DbConnection, get_db_name
from cogent.util.misc import flatten

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Jason Merkin"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

class Database(object):
    """holds the data-base connection and table attributes"""
    def __init__(self, account, species=None, db_type=None, release=None,
            pool_recycle=None, division=None):
        self._tables = {}
        self.db_name = get_db_name(account=account, species=species,
                          release=release, db_type=db_type, division=division)
        if not self.db_name:
            raise RuntimeError, "%s db doesn't exist for '%s' on '%s'" % \
                                        (db_type, species, account.host)
        else:
            self.db_name = self.db_name[0]
        self._db = DbConnection(account=account, db_name=self.db_name,
                                pool_recycle=pool_recycle)
        self._meta = sql.MetaData(self._db)
        self.Type = db_type
        
    def __str__(self):
        return str(self.db_name)
    
    def __cmp__(self, other):
        return cmp(self._db, other._db)
    
    def getTable(self, name):
        """returns the SQLalchemy table instance"""
        table = self._tables.get(name, None)
        if table is None:
            c = self._db.execute("DESCRIBE %s" % name)
            custom_columns = []
            for r in c.fetchall():
                Field = r["Field"]
                Type = r["Type"]
                if "tinyint" in Type:
                    custom_columns.append(sql.Column(Field, sql.Integer))
            try:
                table = sql.Table(name, self._meta, autoload=True,
                                    extend_existing=True, *custom_columns)
            except TypeError:
                # new arg name not supported, try old
                table = sql.Table(name, self._meta, autoload=True,
                                    useexisting=True, *custom_columns)
            
            self._tables[name] = table
        return table
    
    def getDistinct(self, table_name, column):
        """returns the Ensembl data-bases distinct values for the named
        property_type.
        
        Arguments:
            - table_name: the data base table name
            - column: valid values are biotype, status"""
        table = self.getTable(table_name)
        query = sql.select([table.c[column]], distinct=True)
        records = set()
        string_types = str, unicode
        for record in query.execute():
            if type(record) not in string_types and \
                type(record[0]) not in string_types:
                # multi-dimensioned list/tuple
                record = flatten(record)
            elif type(record) not in string_types:
                # list/tuple of strings
                record = tuple(record)
            else:
                # a string
                record = [record]
            
            records.update(record)
        return records
    
    def tableHasColumn(self, table_name, column):
        """returns True if table has column"""
        table = self.getTable(table_name)
        return hasattr(table.c, column)
    
    def getTablesRowCount(self, table_name=None):
        """returns a cogent Table object with the row count for each table
        in the database
        
        Arguments:
            - table_name: database table name. If none, all database tables
              assessed."""
        if type(table_name) == str:
            table_name = (table_name,)
        elif table_name is None:
            self._meta.reflect()
            table_name = self._meta.tables.keys()
        rows = []
        for name in table_name:
            table = self.getTable(name)
            count = table.count().execute().fetchone()[0]
            rows.append(['%s.%s' % (self.db_name, name), count])
        
        return cogent_table.Table(header=['name', 'count'], rows=rows)
    


