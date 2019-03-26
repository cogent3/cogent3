#!/usr/bin/env python
#
# suite of cogent package unit tests.
# run suite by executing this file
#
from warnings import filterwarnings, resetwarnings
import importlib
import doctest
import sys
import os

resetwarnings()
filterwarnings("ignore", category=ResourceWarning, append=True)
filterwarnings("ignore",
               message="using slow exponentiator.+",
               category=UserWarning, append=True)
filterwarnings("ignore",
               message="can't resolve package from.+",
               category=ImportWarning, append=True)


import cogent3.util.unit_test as unittest


__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
               "Hau Ying", "Helen Lindsay", "Jeremy Widmann",
               "Sandra Smit", "Greg Caporaso", "Matthew Wakefield",
               "Ben Kaehler"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def suite():
    modules_to_test = [
        'test_phylo',
        'test_align.test_align',
        'test_app.test_align',
        'test_app.test_evo',
        'test_app.test_composable',
        'test_app.test_io',
        'test_app.test_translate',
        'test_app.test_tree',
        'test_app.test_sample',
        'test_app.test_data_store',
        'test_cluster.test_UPGMA',
        'test_core.test_alphabet',
        'test_core.test_alignment',
        'test_core.test_annotation',
        'test_core.test_core_standalone',
        'test_core.test_genetic_code',
        'test_core.test_info',
        'test_core.test_location',
        'test_core.test_maps',
        'test_core.test_moltype',
        'test_core.test_profile',
        'test_core.test_seq_aln_integration',
        'test_core.test_sequence',
        'test_core.test_tree',
        'test_data.test_molecular_weight',
        'test_evolve.test_best_likelihood',
        'test_evolve.test_bootstrap',
        'test_evolve.test_coevolution',
        'test_evolve.test_ns_substitution_model',
        'test_evolve.test_motifchange',
        'test_evolve.test_substitution_model',
        'test_evolve.test_scale_rules',
        'test_evolve.test_likelihood_function',
        'test_evolve.test_newq',
        'test_evolve.test_pairwise_distance',
        'test_evolve.test_parameter_controller',
        'test_evolve.test_models',
        'test_format.test_bedgraph',
        'test_format.test_fasta',
        'test_format.test_clustal',
        'test_maths.test_fit_function',
        'test_maths.test_geometry',
        'test_maths.test_matrix_logarithm',
        'test_maths.test_matrix_exponential_integration',
        'test_maths.test_period',
        'test_maths.test_stats.test_distribution',
        'test_maths.test_stats.test_information_criteria',
        'test_maths.test_stats.test_period',
        'test_maths.test_stats.test_special',
        'test_maths.test_stats.test_test',
        'test_maths.test_stats.test_ks',
        'test_maths.test_stats.test_util',
        'test_maths.test_stats.test_jackknife',
        'test_maths.test_optimisers',
        'test_maths.test_distance_transform',
        'test_maths.test_svd',
        'test_parse.test_cigar',
        'test_parse.test_clustal',
        'test_parse.test_dialign',
        'test_parse.test_ebi',
        'test_parse.test_fasta',
        'test_parse.test_gbseq',
        'test_parse.test_genbank',
        'test_parse.test_gff',
        'test_parse.test_greengenes',
        'test_parse.test_locuslink',
        'test_parse.test_ncbi_taxonomy',
        'test_parse.test_nexus',
        'test_parse.test_psl',
        'test_parse.test_pamlmatrix',
        'test_parse.test_phylip',
        'test_parse.test_rdb',
        'test_parse.test_record',
        'test_parse.test_record_finder',
        'test_parse.test_tinyseq',
        'test_parse.test_tree',
        'test_parse.test_unigene',
        'test_parse.test_blast_xml',
        'test_util.test_unit_test',
        'test_util.test_array',
        'test_util.test_deserialise',
        'test_util.test_dict2d',
        'test_util.test_misc',
        'test_util.test_recode_alignment',
        'test_util.test_transform',
        'test_recalculation.rst',
        'test_dictarray.rst',
        'test_core.test_features.rst',
        'test_util.test_table.rst',
    ]

    try:
        import matplotlib
    except:
        print("No matplotlib so not running test_draw.py", file=sys.stderr)
    else:
        modules_to_test.append('test_draw')

    assert sys.version_info >= (2, 6)

    alltests = unittest.TestSuite()

    for module in modules_to_test:
        if module.endswith('.rst'):
            module = os.path.join(*module.split(".")[:-1]) + ".rst"
            test = doctest.DocFileSuite(module, optionflags=doctest.REPORT_ONLY_FIRST_FAILURE |
                                        doctest.ELLIPSIS)
        else:
            test = unittest.findTestCases(importlib.import_module(module))
        alltests.addTest(test)
    return alltests


class BoobyTrappedStream(object):

    def __init__(self, output):
        self.output = output

    def write(self, text):
        self.output.write(text)
        raise RuntimeError("Output not allowed in tests")

    def flush(self):
        pass

    def isatty(self):
        return False

if __name__ == '__main__':
    if '--debug' in sys.argv:
        s = suite()
        s.debug()
    else:
        orig = sys.stdout
        if '--output-ok' in sys.argv:
            sys.argv.remove('--output-ok')
        else:
            sys.stdout = BoobyTrappedStream(orig)
        try:
            unittest.main(defaultTest='suite', argv=sys.argv)
        finally:
            sys.stdout = orig
