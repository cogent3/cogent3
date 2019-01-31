from cogent3.util.unit_test import TestCase, main
from cogent3.evolve.models import (JC69, F81, HKY85, TN93, GTR, GN, ssGN,
                                   MG94HKY, MG94GTR, GY94, H04G, H04GK, H04GGK,
                                   DSO78, AH96, AH96_mtmammals, JTT92, WG01,
                                   CNFGTR, CNFHKY, GNC, WG01_matrix, WG01_freqs)

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class CannedModelsTest(TestCase):
    """Check each canned model can actually be instantiated."""

    def _instantiate_models(self, models, **kwargs):
        for model in models:
            model(**kwargs)

    def test_nuc_models(self):
        """excercising nucleotide model construction"""
        self._instantiate_models([JC69, F81, HKY85, GTR, GN, ssGN])

    def test_codon_models(self):
        """excercising codon model construction"""
        self._instantiate_models([CNFGTR, CNFHKY, MG94HKY, MG94GTR, GY94,
                                  H04G, H04GK, H04GGK, GNC])

    def test_aa_models(self):
        """excercising aa model construction"""
        self._instantiate_models([DSO78, AH96, AH96_mtmammals, JTT92, WG01])

    def test_bin_options(self):
        kwargs = dict(with_rate=True, distribution='gamma')
        model = WG01(**kwargs)
        model = GTR(**kwargs)

    def test_empirical_values_roundtrip(self):
        model = WG01()
        assert model.get_motif_probs() == WG01_freqs
        assert (model.calc_exchangeability_matrix('dummy_mprobs') ==
                WG01_matrix).all()

    def test_solved_models(self):
        for klass in [TN93, HKY85, F81]:
            for scaled in [True, False]:
                model = klass(rate_matrix_required=False, do_scaling=scaled)
                model.check_psub_calculations_match()


if __name__ == '__main__':
    main()
