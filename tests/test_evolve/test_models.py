from cogent.util.unit_test import TestCase, main
from cogent.evolve.models import JC69, F81, HKY85, GTR, \
    MG94HKY, MG94GTR, GY94, H04G, H04GK, H04GGK, \
    DSO78, AH96, AH96_mtmammals, JTT92, WG01, CNFGTR, CNFHKY

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class CannedModelsTest(TestCase):
    """Check each canned model can actually be instantiated."""
    
    def _instantiate_models(self, models):
        for model in models:
            model()
    
    def test_nuc_models(self):
        """excercising nucleotide model construction"""
        self._instantiate_models([JC69, F81, HKY85, GTR])
    
    def test_codon_models(self):
        """excercising codon model construction"""
        self._instantiate_models([CNFGTR, CNFHKY, MG94HKY, MG94GTR, GY94,
                                  H04G, H04GK, H04GGK])
    
    def test_aa_models(self):
        """excercising aa model construction"""
        self._instantiate_models([DSO78, AH96, AH96_mtmammals, JTT92, WG01])
    
if __name__ == '__main__':
    unittest.main()
