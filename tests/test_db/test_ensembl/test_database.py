import os

from cogent.util.unit_test import TestCase, main
from cogent.db.ensembl.host import HostAccount, get_ensembl_account
from cogent.db.ensembl.database import Database

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "GPL"
__version__ = "1.5.2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

Release = 68

if 'ENSEMBL_ACCOUNT' in os.environ:
    args = os.environ['ENSEMBL_ACCOUNT'].split()
    host, username, password = args[0:3]
    kwargs = {}
    if len(args) > 3:
        kwargs['port'] = int(args[3])
    account = HostAccount(host, username, password, **kwargs)
else:
    account = get_ensembl_account(release=Release)

class TestDatabase(TestCase):
    
    def test_connect(self):
        human = Database(account=account, release=Release,
                    species='human', db_type='core')
        gene = human.getTable('gene')
    
    def test_get_distinct(self):
        """should return list of strings"""
        db = Database(account=account, release=Release,
                    species='human', db_type='variation')
        tn, tc = 'variation_feature', 'consequence_types'
        expected = set(('3_prime_UTR_variant', 'splice_acceptor_variant',
                        '5_prime_UTR_variant'))
        got = db.getDistinct(tn, tc)
        self.assertNotEquals(set(got) & expected, set())
        
        db = Database(account=account, release=Release,
                    species='human', db_type='core')
        tn, tc = 'gene', 'biotype'
        expected = set(['protein_coding', 'pseudogene', 'processed_transcript', 
          'Mt_tRNA', 'Mt_rRNA', 'IG_V_gene', 'IG_J_gene',
          'IG_C_gene', 'IG_D_gene', 'miRNA', 'misc_RNA', 'snoRNA', 'snRNA', 'rRNA'])
        got = set(db.getDistinct(tn, tc))
        self.assertNotEquals(set(got) & expected, set())
        
        db = Database(account=account, release=Release, db_type='compara')
        got = set(db.getDistinct('homology', 'description'))
        expected = set(['apparent_ortholog_one2one', 'ortholog_many2many',
          'ortholog_one2many', 'ortholog_one2one', 'within_species_paralog'])
        self.assertEquals(len(got&expected), len(expected))
    
    def test_get_table_row_counts(self):
        """should return correct row counts for some tables"""
        expect = {'homo_sapiens_core_68_37.analysis': 61L,
                  'homo_sapiens_core_68_37.seq_region': 55616L,
                  'homo_sapiens_core_68_37.assembly': 102090L,
                  'homo_sapiens_core_68_37.qtl': 0L}
        human = Database(account=account, release=Release,
                    species='human', db_type='core')
        table_names = [n.split('.')[1] for n in expect]
        got = dict(human.getTablesRowCount(table_names).getRawData())
        for dbname in expect:
            self.assertTrue(got[dbname] >= expect[dbname])
    
    def test_table_has_column(self):
        """return correct values for whether a Table has a column"""
        account = get_ensembl_account(release=Release)
        var61 = Database(account=account, release=61, species='human',
            db_type='variation')
        
        var62 = Database(account=account, release=62, species='human',
            db_type='variation')
        
        self.assertTrue(var61.tableHasColumn('transcript_variation',
            'peptide_allele_string'))
        self.assertFalse(var61.tableHasColumn('transcript_variation',
            'pep_allele_string'))
        
        self.assertTrue(var62.tableHasColumn('transcript_variation',
            'pep_allele_string'))
        self.assertFalse(var62.tableHasColumn('transcript_variation',
            'peptide_allele_string'))

if __name__ == "__main__":
    main()
