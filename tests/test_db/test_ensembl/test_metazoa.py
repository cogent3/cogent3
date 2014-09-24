from cogent.db.ensembl.host import HostAccount, get_ensembl_account
from cogent.db.ensembl.compara import Compara, Genome
from cogent.util.unit_test import TestCase, main

__author__ = "Jason Merkin"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Hua Ying"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

Release = 23
account = HostAccount('mysql.ebi.ac.uk','anonymous', '', port=4157)

class MZ_ComparaTestBase(TestCase):
    comp = Compara(['D.grimshawi', 'D.melanogaster'], Release=Release,
                    account=account, division='metazoa')

class MZ_TestCompara(MZ_ComparaTestBase):
    def test_query_genome(self):
        """compara should attach valid genome attributes by common name"""
        brca2 = self.comp.Dmelanogaster.getGeneByStableId("FBgn0050169")
        self.assertEquals(brca2.Symbol.lower(), 'brca2')
    
    def test_get_related_genes(self):
        """should correctly return the related gene regions from each genome"""
        # using sc35, a splicing factor
        sc35 = self.comp.Dmelanogaster.getGeneByStableId("FBgn0040286")
        Orthologs = self.comp.getRelatedGenes(gene_region=sc35,
                Relationship="ortholog_one2one")
        self.assertEquals("ortholog_one2one", Orthologs.Relationships[0])
    
    def test_get_related_genes2(self):
        """should handle case where gene is absent from one of the genomes"""
        # here, it is brca2
        brca2 = self.comp.Dmelanogaster.getGeneByStableId(
                                        StableId='FBgn0050169')
        orthologs = self.comp.getRelatedGenes(gene_region=brca2,
                        Relationship='ortholog_one2one')
        self.assertEquals(len(orthologs.Members),2)
    
    def test_get_collection(self):
        sc35 = self.comp.Dmelanogaster.getGeneByStableId(StableId="FBgn0040286")
        Orthologs = self.comp.getRelatedGenes(gene_region=sc35,
                        Relationship="ortholog_one2one")
        collection = Orthologs.getSeqCollection()
        self.assertTrue(len(collection.Seqs[0])> 1000)
    
class MZ_Genome(TestCase):
    def test_get_general_release(self):
        """should correctly infer the general release"""
        rel_lt_65 = Genome('D.melanogaster', Release=22, account=account)
        self.assertEqual(rel_lt_65.GeneralRelease, 75)
        self.assertEqual(rel_lt_65.CoreDb.db_name, 'drosophila_melanogaster_core_22_75_546')
        
        rel_gt_65 = Genome('D.melanogaster', Release=23, account=account)
        self.assertEqual(rel_gt_65.GeneralRelease, 76)
        self.assertEqual(rel_gt_65.CoreDb.db_name, 'drosophila_melanogaster_core_23_76_546')



if __name__ == "__main__":
    main()
