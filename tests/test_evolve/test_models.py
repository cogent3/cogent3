from unittest import TestCase, main

from cogent3.evolve import models as models_module
from cogent3.evolve.models import (
    AH96,
    CNFGTR,
    CNFHKY,
    DSO78,
    F81,
    GN,
    GNC,
    GTR,
    GY94,
    H04G,
    H04GGK,
    H04GK,
    HKY85,
    JC69,
    JTT92,
    MG94GTR,
    MG94HKY,
    TN93,
    WG01,
    AH96_mtmammals,
    WG01_freqs,
    WG01_matrix,
    available_models,
    codon_models,
    get_model,
    models,
    nucleotide_models,
    protein_models,
    ssGN,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.8.28a"
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
        self._instantiate_models(
            [CNFGTR, CNFHKY, MG94HKY, MG94GTR, GY94, H04G, H04GK, H04GGK, GNC]
        )

    def test_aa_models(self):
        """excercising aa model construction"""
        self._instantiate_models([DSO78, AH96, AH96_mtmammals, JTT92, WG01])

    def test_bin_options(self):
        kwargs = dict(with_rate=True, distribution="gamma")
        model = WG01(**kwargs)
        model = GTR(**kwargs)

    def test_empirical_values_roundtrip(self):
        model = WG01()
        assert model.get_motif_probs() == WG01_freqs
        assert (model.calc_exchangeability_matrix("dummy_mprobs") == WG01_matrix).all()

    def test_solved_models(self):
        for klass in [TN93, HKY85, F81]:
            model = klass(rate_matrix_required=False)
            model.check_psub_calculations_match()

    def test_get_model(self):
        """get_models successfully creates model instances"""
        for name in models:
            model = get_model(name)

        # just returns query if it's already a substitution model
        for mod in (CNFGTR(), WG01(), GN()):
            got = get_model(mod)
            self.assertEqual(id(got), id(mod))

        with self.assertRaises(ValueError):
            # unknown model raises exception
            _ = get_model("blah")


def get_sample_model_types(mod_type=None):
    opts = dict(
        codon=codon_models, nucleotide=nucleotide_models, protein=protein_models
    )
    return opts.get(mod_type, models)


class AvailableModelsTest(TestCase):
    def test_model_abbreviation(self):
        """make sure getting model abbreviations that exist"""
        got = set(available_models().tolist("Abbreviation"))
        expect = set(["JC69", "CNFGTR", "DSO78"])
        self.assertTrue(expect < got)

    def test_model_by_type(self):
        """correctly obtain models by type"""
        for model_type in ["codon", "nucleotide", "protein"]:
            table = available_models(model_type)
            got = table.distinct_values("Model Type")
            self.assertEqual(got, {model_type})

    def test_model_description(self):
        """correctly grabs function descriptions"""
        all_available = available_models()
        for abbrev, desc in all_available.tolist(["Abbreviation", "Description"]):
            func = getattr(models_module, abbrev)
            doc = func.__doc__.split()
            self.assertEqual(desc.split(), doc)


if __name__ == "__main__":
    main()
