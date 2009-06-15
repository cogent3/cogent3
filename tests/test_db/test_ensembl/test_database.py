import os

from cogent.util.unit_test import TestCase, main
from cogent.db.ensembl.host import HostAccount, get_ensembl_account
from cogent.db.ensembl.database import Database

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "GPL"
__version__ = "1.4.0.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

Release=52

if 'ENSEMBL_ACCOUNT' in os.environ:
    username, password = os.environ['ENSEMBL_ACCOUNT'].split()
    account = HostAccount('127.0.0.1', username, password)
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
        tn, tc = 'variation_feature', 'consequence_type'
        expected = set((('3PRIME_UTR', 'ESSENTIAL_SPLICE_SITE'),
                    ('3PRIME_UTR', 'SPLICE_SITE'),
                    ('5PRIME_UTR', 'ESSENTIAL_SPLICE_SITE')))
        self.assertNotEquals(set(db.getDistinct(tn, tc)) & expected, set())
        
        db = Database(account=account, release=Release,
                    species='human', db_type='core')
        tn, tc = 'gene', 'biotype'
        expected = set(['rRNA', 'IG_pseudogene', 'rRNA_pseudogene',
          'snoRNA_pseudogene', 'snRNA_pseudogene', 'IG_V_gene', 'misc_RNA',
          'Mt_tRNA', 'misc_RNA_pseudogene', 'IG_J_gene', 'scRNA',
          'IG_C_gene', 'protein_coding', 'Mt_rRNA', 'scRNA_pseudogene',
          'pseudogene', 'snRNA',
          'unprocessed_pseudogene', 'tRNA_pseudogene', 'miRNA', 'IG_D_gene',
          'retrotransposed', 'Mt_tRNA_pseudogene', 'snoRNA',
          'processed_pseudogene', 'miRNA_pseudogene'])
        got = set(db.getDistinct(tn, tc))
        self.assertEquals(len(got&expected), len(expected))
        
        db = Database(account=account, release=Release, db_type='compara')
        got = set(db.getDistinct('homology', 'description'))
        expected = set(['apparent_ortholog_one2one', 'between_species_paralog',
        'ortholog_many2many', 'ortholog_one2many', 'ortholog_one2one',
        'within_species_paralog'])
        self.assertEquals(len(got&expected), len(expected))
        

if __name__ == "__main__":
    main()
