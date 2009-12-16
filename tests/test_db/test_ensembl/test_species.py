from cogent.util.unit_test import TestCase, main
from cogent.db.ensembl.species import Species

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

class TestSpeciesNamemaps(TestCase):
    def test_get_name_type(self):
        """should return the (latin|common) name given a latin, common or ensembl
        db prefix names"""
        self.assertEqual(Species.getSpeciesName("human"), "Homo sapiens")
        self.assertEqual(Species.getSpeciesName("homo_sapiens"), "Homo sapiens")
        self.assertEqual(Species.getCommonName("Mus musculus"), "Mouse")
        self.assertEqual(Species.getCommonName("mus_musculus"), "Mouse")
    
    def test_get_ensembl_format(self):
        """should take common or latin names and return the corresponding
        ensembl db prefix"""
        self.assertEqual(Species.getEnsemblDbPrefix("human"), "homo_sapiens")
        self.assertEqual(Species.getEnsemblDbPrefix("mouse"), "mus_musculus")
        self.assertEqual(Species.getEnsemblDbPrefix("Mus musculus"),
                                                "mus_musculus")
    
    def test_add_new_species(self):
        """should correctly add a new species/common combination and infer the
        correct ensembl prefix"""
        species_name, common_name = "Otolemur garnettii", "Bushbaby"
        Species.amendSpecies(species_name, common_name)
        self.assertEqual(Species.getSpeciesName(species_name), species_name)
        self.assertEqual(Species.getSpeciesName("Bushbaby"), species_name)
        self.assertEqual(Species.getSpeciesName(common_name), species_name)
        self.assertEqual(Species.getCommonName(species_name), common_name)
        self.assertEqual(Species.getCommonName("Bushbaby"), common_name)
        self.assertEqual(Species.getEnsemblDbPrefix("Bushbaby"), "otolemur_garnettii")
        self.assertEqual(Species.getEnsemblDbPrefix(species_name), "otolemur_garnettii")
        self.assertEqual(Species.getEnsemblDbPrefix(common_name), "otolemur_garnettii")
    
    def test_amend_existing(self):
        """should correctly amend an existing species"""
        species_name = 'Ochotona princeps'
        common_name1 = 'american pika'
        common_name2 = 'pika'
        ensembl_pref = 'ochotona_princeps'
        Species.amendSpecies(species_name, common_name1)
        self.assertEqual(Species.getCommonName(species_name),common_name1)
        Species.amendSpecies(species_name, common_name2)
        self.assertEqual(Species.getSpeciesName(common_name2), species_name)
        self.assertEqual(Species.getSpeciesName(ensembl_pref), species_name)
        self.assertEqual(Species.getCommonName(species_name), common_name2)
        self.assertEqual(Species.getCommonName(ensembl_pref), common_name2)
        self.assertEqual(Species.getEnsemblDbPrefix(species_name),ensembl_pref)
        self.assertEqual(Species.getEnsemblDbPrefix(common_name2),ensembl_pref)

if __name__ == "__main__":
    main()