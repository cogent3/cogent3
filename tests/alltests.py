#!/usr/bin/env python
#
# suite of cogent package unit tests.
# run suite by executing this file
#
import doctest, cogent.util.unit_test as unittest, sys, os
from cogent.util.misc import app_path

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
                    "Hau Ying", "Helen Lindsay", "Jeremy Widmann",
                    "Sandra Smit", "Greg Caporaso", "Matthew Wakefield",
                "Ben Kaehler"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def my_import(name):
    """Imports a module, possibly qualified with periods. Returns the module.
    
    __import__ only imports the top-level module.
    
    Recipe from python documentation at:
    http://www.python.org/doc/2.4/lib/built-in-funcs.html
    """
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod

def module_present(modules):
    """returns True if dependencies present"""
    if type(modules) == str:
        modules = [modules]
    try:
        for module in modules:
            mod = __import__(module)
    except ImportError:
        return False
    return True

def suite():
    modules_to_test = [
        'test_recalculation.rst',
        'test_phylo',
        'test_dictarray.rst',
        'test_align.test_align',
        'test_cluster.test_UPGMA',
        'test_core.test_alphabet',
        'test_core.test_alignment',
        'test_core.test_annotation',
        'test_core.test_bitvector',
        'test_core.test_core_standalone',
        'test_core.test_features.rst',
        'test_core.test_genetic_code',
        'test_core.test_info',
        'test_core.test_location',
        'test_core.test_maps',
        'test_core.test_moltype',
        'test_core.test_profile',
        'test_core.test_seq_aln_integration',
        'test_core.test_sequence',
        'test_core.test_tree',
        'test_core.test_usage',
        'test_data.test_molecular_weight',
        'test_evolve.test_best_likelihood',
        'test_evolve.test_bootstrap',
        'test_evolve.test_coevolution',
        'test_evolve.test_models',
        'test_evolve.test_motifchange',
        'test_evolve.test_substitution_model',
        'test_evolve.test_scale_rules',
        'test_evolve.test_likelihood_function',
        'test_evolve.test_newq',
        'test_evolve.test_pairwise_distance',
        'test_evolve.test_parameter_controller',
        'test_format.test_bedgraph',
        'test_format.test_fasta',
        'test_maths.test_fit_function',
        'test_maths.test_geometry',
        'test_maths.test_matrix_logarithm',
        'test_maths.test_matrix_exponential_integration',
        'test_maths.test_period',
        'test_maths.test_matrix.test_distance',
        'test_maths.test_stats.test_distribution',
        'test_maths.test_stats.test_histogram',
        'test_maths.test_stats.test_information_criteria',
        'test_maths.test_stats.test_period',
        'test_maths.test_stats.test_special',
        'test_maths.test_stats.test_test',
        'test_maths.test_stats.test_ks',
        'test_maths.test_stats.test_util',
        'test_maths.test_stats.test_cai.test_adaptor',
        'test_maths.test_stats.test_cai.test_get_by_cai',
        'test_maths.test_stats.test_cai.test_util',
        'test_maths.test_optimisers',
        'test_maths.test_distance_transform',
        'test_motif.test_util',
        'test_parse.test_aaindex',
        'test_parse.test_bowtie',
        'test_parse.test_cigar',
        'test_parse.test_clustal',
        'test_parse.test_cut',
        'test_parse.test_cutg',
        'test_parse.test_dialign',
        'test_parse.test_ebi',
        'test_parse.test_fasta',
        'test_parse.test_gbseq',
        'test_parse.test_gibbs',
        'test_parse.test_genbank',
        'test_parse.test_gff',
        'test_parse.test_greengenes',
        'test_parse.test_locuslink',
        'test_parse.test_meme',
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
        'test_util.test_unit_test',
        'test_util.test_array',
        'test_util.test_dict2d',
        'test_util.test_misc',
        'test_util.test_organizer',
        'test_util.test_recode_alignment',
        'test_util.test_table.rst',
        'test_util.test_transform',
        ]

    try:
        import matplotlib
    except:
        print >> sys.stderr, "No matplotlib so not running test_draw.py"
    else:
        modules_to_test.append('test_draw')

    # we now toggle the db tests, based on an environment flag
    if int(os.environ.get('TEST_DB', 0)):
        db_tests = ['test_db.test_ncbi', 'test_db.test_pdb',
                        'test_db.test_rfam', 'test_db.test_util']

        # we check for an environment flag for ENSEMBL
        # we expect this to have the username and account for a localhost
        # installation of the Ensembl MySQL databases
        if 'ENSEMBL_ACCOUNT' in os.environ:
            # check for cogent.db.ensembl dependencies
            test_ensembl = True
            if not (module_present('MySQLdb') or module_present('mysql')):
                test_ensembl = False
                print >> sys.stderr, \
                    "Module 'mysql-connector-python' and 'MySQL-python' not "\
                        "present: skipping test"
            
            if not module_present('sqlalchemy'):
                test_ensembl = False
                print >> sys.stderr, \
                    "Module 'sqlalchemy' not present: skipping test"

            if test_ensembl:
                db_tests += ['test_db.test_ensembl.test_assembly',
                     'test_db.test_ensembl.test_database',
                     'test_db.test_ensembl.test_compara',
                     'test_db.test_ensembl.test_genome',
                     'test_db.test_ensembl.test_host',
                     'test_db.test_ensembl.test_metazoa',
                     'test_db.test_ensembl.test_species',
                     'test_db.test_ensembl.test_feature_level']
        else:
            print >> sys.stderr, "Environment variable ENSEMBL_ACCOUNT not "\
            "set: skipping db.ensembl tests"

        for db_test in db_tests:
            modules_to_test.append(db_test)
    else:
        print >> sys.stderr, \
                "Environment variable TEST_DB=1 not set: skipping db tests"
    
    assert sys.version_info >= (2, 6)
    
    alltests = unittest.TestSuite()
    
    for module in modules_to_test:
        if module.endswith('.rst'):
            module = os.path.join(*module.split(".")[:-1]) + ".rst"
            test = doctest.DocFileSuite(module, optionflags=
                doctest.REPORT_ONLY_FIRST_FAILURE |
                doctest.ELLIPSIS)
        else:
            test = unittest.findTestCases(my_import(module))
        alltests.addTest(test)
    return alltests

class BoobyTrappedStream(object):
    def __init__(self, output):
        self.output = output

    def write(self, text):
        self.output.write(text)
        raise RuntimeError, "Output not allowed in tests"
        
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
