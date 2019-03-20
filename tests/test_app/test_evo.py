from unittest import TestCase, main
from cogent3 import LoadTree
from cogent3.app import evo as evo_app

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestModel(TestCase):
    basedir = 'data'

    def test_model_str(self):
        """correct str representation"""
        model = evo_app.model('HKY85', time_het='max')
        got = str(model)
        self.assertEqual(got, ("model(type='model', sm='HKY85', tree=None, "
                               "name=None, sm_args=None, lf_args=None, "
                               "time_het='max', param_rules=None, "
                               "opt_args=None, split_codons=False, "
                               "show_progress=False)"))

    def test_model_tree(self):
        """allows tree to be string, None or tree"""
        treestring = '(a,b,c)'
        for tree in (treestring, LoadTree(treestring=treestring), None):
            mod = evo_app.model('HKY85', tree=tree)
            expect = None if tree is None else LoadTree(treestring=treestring)
            self.assertIsInstance(mod._tree, expect.__class__)

    def test_hypothesis_str(self):
        """correct str representation"""
        model1 = evo_app.model('HKY85')
        model2 = evo_app.model('HKY85', time_het='max')
        hyp = evo_app.hypothesis(model1, model2)
        got = str(hyp)
        expect = ("hypothesis(type='hypothesis', null='HKY85', "
                  "alternates=(model(type='model', sm='HKY85', tree=None, "
                  "name=None, sm_args=None, lf_args=None, "
                  "time_het='max', param_rules=None, opt_args=None,"
                  " split_codons=False, show_progress=False),))")
        self.assertEqual(got, expect)


if __name__ == '__main__':
    main()
