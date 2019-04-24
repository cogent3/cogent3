import os
from cogent3 import LoadTree, LoadSeqs, DNA
from unittest import TestCase, main

from cogent3.app import tree as tree_app

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestTree(TestCase):
    def test_scale_tree_lengths(self):
        """correctly scales tree lengths"""
        with self.assertRaises(AssertionError):
            _ = tree_app.scale_branches(nuc_to_codon=True,
                                        codon_to_nuc=True)

        scale_to_codon = tree_app.scale_branches(nuc_to_codon=True)
        tree = LoadTree(treestring='(a:3,b:6,c:9)')
        scale_to_codon = tree_app.scale_branches(nuc_to_codon=True)
        d = scale_to_codon(tree)
        got = {e.name: e.length for e in d.get_edge_vector(include_root=False)}
        expect = {'a': 1.0, 'b': 2.0, 'c': 3.0}
        self.assertEqual(got, expect)

        scale_from_codon = tree_app.scale_branches(codon_to_nuc=True)
        d = scale_from_codon(d)
        got = {e.name: e.length for e in d.get_edge_vector(include_root=False)}
        expect = {'a': 3.0, 'b': 6.0, 'c': 9.0}
        self.assertEqual(got, expect)

        by_scalar = tree_app.scale_branches(scalar=0.5)
        d = by_scalar(tree)
        got = {e.name: e.length for e in d.get_edge_vector(include_root=False)}
        expect = {'a': 6.0, 'b': 12.0, 'c': 18.0}
        self.assertEqual(got, expect)

        # handle case where a length is not defined, setting to minimum
        min_length = tree_app.scale_branches(min_length=66)
        tree = LoadTree(treestring='(a:3,b:6,c)')
        new = min_length(tree)
        got = {e.name: e.length for e in new.get_edge_vector(
            include_root=False)}
        expect = {'a': 3.0, 'b': 6.0, 'c': 66.0}
        self.assertEqual(got, expect)

    def test_quick_tree(self):
        """correctly calc a nj tree"""
        path = os.path.join(os.path.abspath(__file__).split('test_app')[0],
                            'data/brca1_5.paml')
        aln = LoadSeqs(path, moltype=DNA)
        quick = tree_app.quick_tree()
        tree = quick(aln)
        self.assertEqual(set(tree.get_tip_names()), set(aln.names))

    def test_uniformize_tree(self):
        """equivalent topologies should be the same"""
        a = LoadTree(treestring="(a,(b,c),(d,e))")
        b = LoadTree(treestring="(e,d,(a,(b,c)))")
        make_uniform = tree_app.uniformize_tree(root_at='c',
                                                ordered_names=list('abcde'))
        u_a = make_uniform(a).get_newick()
        u_b = make_uniform(b).get_newick()
        self.assertTrue(u_a == u_b)
        # but different ones different
        c = LoadTree(treestring="(e,c,(a,(b,d)))")
        u_c = make_uniform(c).get_newick()
        self.assertFalse(u_a == u_c)


if __name__ == '__main__':
    main()
