import os

from cogent.util.unit_test import TestCase, main
from cogent.db.ensembl.name import EnsemblDbName
from cogent.db.ensembl.host import get_db_name, get_latest_release,\
                        DbConnection, HostAccount, get_ensembl_account
from cogent.db.ensembl.species import Species

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

Release = 56

if 'ENSEMBL_ACCOUNT' in os.environ:
    username, password = os.environ['ENSEMBL_ACCOUNT'].split()
    account = HostAccount('127.0.0.1', username, password)
else:
    account = get_ensembl_account(release=Release)

class TestEnsemblDbName(TestCase):
    def test_cmp_name(self):
        """should validly compare names by attributes"""
        n1 = EnsemblDbName('homo_sapiens_core_46_36h')
        n2 = EnsemblDbName('homo_sapiens_core_46_36h')
        self.assertEqual(n1, n2)
    
    def test_name_without_build(self):
        """should correctly handle a db name without a build"""
        n = EnsemblDbName("pongo_pygmaeus_core_49_1")
        self.assertEqual(n.Prefix, "pongo_pygmaeus")
        self.assertEqual(n.Type, "core")
        self.assertEqual(n.Build, '1')

class TestDBconnects(TestCase):
    
    def test_get_ensembl_account(self):
        """return an HostAccount with correct port"""
        for release in [48, '48', None]:
            act_new = get_ensembl_account(release=release)
            self.assertEqual(act_new.port, 5306)
        
        for release in [45, '45']:
            act_old = get_ensembl_account(release=45)
            self.assertEqual(act_old.port, 3306)
    
    def test_getdb(self):
        """should discover human entries correctly"""
        for name, db_name in [("human", "homo_sapiens_core_49_36k"),
                      ("mouse", "mus_musculus_core_49_37b"),
                      ("rat", "rattus_norvegicus_core_49_34s"),
                      ("platypus", "ornithorhynchus_anatinus_core_49_1f")]:
            result = get_db_name(species=name, db_type="core", release='49')
            self.assertEqual(len(result), 1)
            result = result[0]
            self.assertEqual(result.Name, db_name)
            self.assertEqual(result.Release, '49')
    
    def test_latest_release_number(self):
        """should correctly the latest release number"""
        self.assertGreaterThan(get_latest_release(), "53")
    
    def test_get_all_available(self):
        """should return a listing of all the available databases on the
        indicated server"""
        available = get_db_name()
        # make sure we have a compara db present -- a crude check on
        # correctness
        one_valid = False
        for db in available:
            if db.Type == "compara":
                one_valid = True
                break
        self.assertEqual(one_valid, True)
        # now check that when we request available under a specific version
        # that we only receive valid ones back
        available = get_db_name(release="46")
        for db in available:
            self.assertEqual(db.Release, '46')
    
    def test_active_connections(self):
        """connecting to a database on a specified server should be done once
        only, but same database on a different server should be done"""
        ensembl_acct = get_ensembl_account(release='46')
        engine1 = DbConnection(account=ensembl_acct,db_name="homo_sapiens_core_46_36h")
        engine2 = DbConnection(account=ensembl_acct,db_name="homo_sapiens_core_46_36h")
        self.assertEqual(engine1, engine2)

if __name__ == "__main__":
    main()