#!/usr/bin/env python
#
# suite of cogent package unit tests.
# run suite by executing this file
#

import doctest, cogent.util.unit_test as unittest, sys, os
from cogent.util.misc import app_path

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
                    "Hau Ying", "Helen Lindsay", "Jeremy Widmann",
                    "Sandra Smit", "Greg Caporaso", "Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.0.1"
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

def suite():
    modules_to_test = [
        'test_recalculation.rest',
        'test_phylo',
        'test_dictarray.rest',
        'test_align.test_align',
        'test_align.test_algorithm',
        'test_align.test_weights.test_methods',
        'test_align.test_weights.test_util',
        'test_app.test_parameters',
        'test_app.test_util',
        'test_cluster.test_metric_scaling',
        'test_cluster.test_UPGMA',
        'test_core.test_alphabet',
        'test_core.test_alignment',
        'test_core.test_annotation',
        'test_core.test_bitvector',
        'test_core.test_core_standalone',
        'test_core.test_features.rest',
        'test_core.test_genetic_code',
        'test_core.test_info',
        'test_core.test_location',
        'test_core.test_maps',
        'test_core.test_moltype',
        'test_core.test_profile',
        'test_core.test_seq_aln_integration',
        'test_core.test_sequence',
        'test_core.test_tree',
        'test_core.test_tree2',
        'test_core.test_usage',
        'test_data.test_molecular_weight',
        'test_evolve.test_best_likelihood',
        'test_evolve.test_bootstrap',
        'test_evolve.test_motifchange',
        'test_evolve.test_substitution_model',
        'test_evolve.test_scale_rules',
        'test_evolve.test_likelihood_function',
        'test_evolve.test_parameter_controller',
        'test_format.test_mage',
        'test_format.test_fasta',
        'test_format.test_pdb_color',
        'test_maths.test_geometry',
        'test_maths.test_matrix_logarithm',
        'test_maths.test_matrix.test_distance',
        'test_maths.test_stats.test_distribution',
        'test_maths.test_stats.test_histogram',
        'test_maths.test_stats.test_special',
        'test_maths.test_stats.test_test',
        'test_maths.test_stats.test_ks',
        'test_maths.test_stats.test_util',
        'test_maths.test_stats.test_cai.test_adaptor',
        'test_maths.test_stats.test_cai.test_get_by_cai',
        'test_maths.test_stats.test_cai.test_util',
        'test_maths.test_optimisers',
        'test_motif.test_util',
        'test_parse.test_aaindex',
        'test_parse.test_agilent_microarray',
        'test_parse.test_blast',
        'test_parse.test_bpseq',
        'test_parse.test_cigar',
        'test_parse.test_clustal',
        'test_parse.test_column',
        'test_parse.test_comrna',
        'test_parse.test_consan',
        'test_parse.test_cove',
        'test_parse.test_ct',
        'test_parse.test_cut',
        'test_parse.test_cutg',
        'test_parse.test_dialign',
        'test_parse.test_ebi',
        'test_parse.test_fasta',
        'test_parse.test_gibbs',
        'test_parse.test_genbank',
        'test_parse.test_gff',
        'test_parse.test_ilm',
        'test_parse.test_locuslink',
        'test_parse.test_mage',
        'test_parse.test_meme',
        'test_parse.test_ncbi_taxonomy',
        'test_parse.test_nexus',
        'test_parse.test_nupack',
        'test_parse.test_phylip',
        'test_parse.test_pknotsrg',
        'test_parse.test_rdb',
        'test_parse.test_record',
        'test_parse.test_record_finder',
        'test_parse.test_rfam',
        'test_parse.test_rnaalifold',
        'test_parse.test_rna_fold',
        'test_parse.test_rnaview',
        'test_parse.test_rnaforester',
        'test_parse.test_sprinzl',
        'test_parse.test_tree',
        'test_parse.test_unigene',
        'test_seqsim.test_analysis',
        'test_seqsim.test_birth_death',
        'test_seqsim.test_markov',
        'test_seqsim.test_microarray',
        'test_seqsim.test_microarray_normalize',
        'test_seqsim.test_randomization',
        'test_seqsim.test_searchpath',
        'test_seqsim.test_sequence_generators',
        'test_seqsim.test_tree',
        'test_seqsim.test_usage',
        'test_struct.test_knots',
        'test_struct.test_pairs_util',
        'test_struct.test_rna2d',
        'test_util.test_unit_test',
        'test_util.test_array',
        'test_util.test_dict2d',
        'test_util.test_misc',
        'test_util.test_organizer',
        'test_util.test_table.rest',
        'test_util.test_transform',
        ]
        
    try:
        import reportlab
    except:
        print >>sys.stderr, "No reportlab so not running test_draw.py"
    else:
        # modules_to_test.append('test_draw')
        pass
    #Try importing modules for app controllers
    apps = [('blastall','test_blast'),
            ('carnac', 'test_carnac'),
            ('clearcut','test_clearcut'),
            ('clustalw','test_clustalw'),
            ('cmfinder.pl', 'test_cmfinder'),
            ('comrna', 'test_comrna'),
            ('contrafold', 'test_contrafold'),
            ('covea', 'test_cove'),
            ('dialign2-2', 'test_dialign'),
            ('dynalign', 'test_dynalign'),
            ('foldalign', 'test_foldalign'),
            ('ilm', 'test_ilm'),
            ('knetfold.pl', 'test_knetfold'),
            ('mafft', 'test_mafft'),
            ('mfold', 'test_mfold'),
            ('muscle', 'test_muscle'),
            ('Fold.out', 'test_nupack'),
            ('findphyl', 'test_pfold'),
            ('pknotsRG-1.2-i386-linux-static', 'test_pknotsrg'),
            ('RNAalifold', 'test_rnaalifold'),
            ('rnaview', 'test_rnaview'),
            ('RNAfold', 'test_vienna_package'),
            ('raxmlHPC', 'test_raxml'),
            ('sfold.X86_64.LINUX', 'test_sfold'),
            ('hybrid-ss-min', 'test_unafold'),
            ('mlagan', 'test_lagan'),
            ]
    for app, test_name in apps:
        if app_path(app):
            modules_to_test.append('test_app.'+test_name)
        else:
            print >>sys.stderr,"Can't find %s executable: skipping test" % app

    if app_path('muscle'):
        modules_to_test.append('test_format.test_pdb_color')
    
    # we now toggle the db tests, based on an evironment flag
    if int(os.environ.get('TEST_DB', 0)):
        for db_test in ['test_db.test_ncbi', 'test_db.test_pdb',
                            'test_db.test_rfam', 'test_db.test_util']:
            modules_to_test.append(db_test)
    else:
        print >>sys.stderr, "Environment variable TEST_DB=1 not set: skipping db tests"
    
    assert sys.version_info >= (2, 4)
    
    alltests = unittest.TestSuite()
    
    for module in modules_to_test:
        if module.endswith('.rest'):
            module = os.path.join(*module.split(".")[:-1])+".rest"
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
