#!/usr/bin/env python

import unittest
import os

from cogent import LoadTree
from cogent.parse.tree import DndParser

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Andrew Butterfield",
                    "Matthew Wakefield", "Daniel McDonald", "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

base_path = os.getcwd()
data_path = os.path.join(base_path, 'data')

# Missing tests: edge attributes (Name, Length, Children) only get tested
# in passing by some of these tests.  See also xxx's

class TreeInterfaceForLikelihoodFunction(unittest.TestCase):
    
    default_newick = "((A:1,B:2)ab:3,((C:4,D:5)cd,E:6)cde:7)"
    
    def _maketree(self, treestring=None):
        if treestring is None:
            treestring = self.default_newick
        return LoadTree(treestring=treestring, underscore_unmunge=True)
    
    def setUp(self):
        self.default_tree = self._maketree()
    
    def test_getEdgeNames(self):
        tree = self._maketree()
        for (a, b, outgroup, result) in [
                ('A', 'B', None, ['A', 'B']),
                ('E', 'C', None, ['C', 'D', 'cd', 'E']),
                ('C', 'D', 'E', ['C', 'D'])]:
            self.assertEqual(tree.getEdgeNames(
                a, b, True, False, outgroup), result)
    
    def test_parser(self):
        """nasty newick"""
        nasty = "( (A :1.0,'B (b)': 2) [com\nment]pair:3,'longer name''s':4)dash_ed;"
        nice = "((A:1.0,'B (b)':2.0)pair:3.0,'longer name''s':4.0)dash_ed;"
        tree = self._maketree(nasty)
        tidied = tree.getNewick(with_distances=1)
        self.assertEqual(tidied, nice)
    
    # Likelihood Function Interface
    
    def test_getEdgeNames(self):
        tree = self.default_tree
        clade = tree.getEdgeNames('C', 'E', getstem=0, getclade=1)
        clade.sort()
        self.assertEqual(clade, ['C', 'D', 'E', 'cd'])
        
        all = tree.getEdgeNames('C', 'E', getstem=1, getclade=1)
        all.sort()
        self.assertEqual(all, ['C', 'D', 'E', 'cd', 'cde'])
        
        stem = tree.getEdgeNames('C', 'E', getstem=1, getclade=0)
        self.assertEqual(stem, ['cde'])
    
    def test_getEdgeNamesUseOutgroup(self):
        t1 = LoadTree(treestring="((A,B)ab,(F,(C,D)cd)cdf,E)root;")
        # a, e, ogroup f
        t2 = LoadTree(treestring="((E,(A,B)ab)abe,F,(C,D)cd)root;")
        expected = ['A', 'B', 'E', 'ab']
        for t in [t1, t2]:
            edges = t.getEdgeNames('A', 'E', getstem = False, getclade = True,
                                        outgroup_name = "F")
            edges.sort()
            self.assertEqual(expected, edges)
    
    def test_getConnectingNode(self):
        tree = self.default_tree
        self.assertEqual(tree.getConnectingNode('A', 'B').Name, 'ab')
        self.assertEqual(tree.getConnectingNode('A', 'C').Name, 'root')
    
    def test_getNodeMatchingName(self):
        tree = self.default_tree
        for (name, expect_tip) in [('A', True), ('ab', False)]:
            edge = tree.getNodeMatchingName(name)
            self.assertEqual(edge.Name, name)
            self.assertEqual(edge.istip(), expect_tip)
    
    def test_getEdgeVector(self):
        tree = self.default_tree
        names = [e.Name for e in tree.getEdgeVector()]
        self.assertEqual(names,
            ['A', 'B', 'ab', 'C', 'D', 'cd', 'E', 'cde', 'root'])
    
    def test_getNewickRecursive(self):
        orig = "((A:1.0,B:2.0)ab:3.0,((C:4.0,D:5.0)cd:6.0,E:7.0)cde:8.0)all;"
        unlen = "((A,B)ab,((C,D)cd,E)cde)all;"
        tree = self._maketree(orig)
        self.assertEqual(tree.getNewickRecursive(with_distances=1), orig)
        self.assertEqual(tree.getNewickRecursive(), unlen)

        tree.Name = "a'l"
        ugly_name = "((A,B)ab,((C,D)cd,E)cde)a'l;"
        ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a''l';"
        self.assertEqual(tree.getNewickRecursive(escape_name=True), \
                         ugly_name_esc)
        self.assertEqual(tree.getNewickRecursive(escape_name=False), ugly_name)
   
        tree.Name = "a_l"
        ugly_name = "((A,B)ab,((C,D)cd,E)cde)a_l;"
        ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a_l';"
        self.assertEqual(tree.getNewickRecursive(escape_name=True), \
                         ugly_name_esc)
        self.assertEqual(tree.getNewickRecursive(escape_name=False), ugly_name)
   
        tree.Name = "a l"
        ugly_name = "((A,B)ab,((C,D)cd,E)cde)a l;"
        ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)a_l;"
        self.assertEqual(tree.getNewickRecursive(escape_name=True), \
                         ugly_name_esc)
        self.assertEqual(tree.getNewickRecursive(escape_name=False), ugly_name)
   
        tree.Name = "'a l'"
        quoted_name = "((A,B)ab,((C,D)cd,E)cde)'a l';"
        quoted_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a l';"
        self.assertEqual(tree.getNewickRecursive(escape_name=True), \
                         quoted_name_esc)
        self.assertEqual(tree.getNewickRecursive(escape_name=False),quoted_name)
   
    def test_getNewick(self):
        orig = "((A:1.0,B:2.0)ab:3.0,((C:4.0,D:5.0)cd:6.0,E:7.0)cde:8.0)all;"
        unlen = "((A,B)ab,((C,D)cd,E)cde)all;"
        tree = self._maketree(orig)
        self.assertEqual(tree.getNewick(with_distances=1), orig)
        self.assertEqual(tree.getNewick(), unlen)

        tree.Name = "a'l"
        ugly_name = "((A,B)ab,((C,D)cd,E)cde)a'l;"
        ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a''l';"
        self.assertEqual(tree.getNewick(escape_name=True), ugly_name_esc)
        self.assertEqual(tree.getNewick(escape_name=False), ugly_name)
   
        tree.Name = "a_l"
        ugly_name = "((A,B)ab,((C,D)cd,E)cde)a_l;"
        ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a_l';"
        self.assertEqual(tree.getNewickRecursive(escape_name=True), \
                         ugly_name_esc)
        self.assertEqual(tree.getNewickRecursive(escape_name=False), ugly_name)
   
        tree.Name = "a l"
        ugly_name = "((A,B)ab,((C,D)cd,E)cde)a l;"
        ugly_name_esc = "((A,B)ab,((C,D)cd,E)cde)a_l;"
        self.assertEqual(tree.getNewickRecursive(escape_name=True), \
                         ugly_name_esc)
        self.assertEqual(tree.getNewickRecursive(escape_name=False), ugly_name)
   
        tree.Name = "'a l'"
        quoted_name = "((A,B)ab,((C,D)cd,E)cde)'a l';"
        quoted_name_esc = "((A,B)ab,((C,D)cd,E)cde)'a l';"
        self.assertEqual(tree.getNewick(escape_name=True), \
                         quoted_name_esc)
        self.assertEqual(tree.getNewick(escape_name=False),quoted_name)
   
    def test_XML(self):
        # should add some non-length parameters
        orig = self.default_tree
        xml = orig.getXML()
        parsed = LoadTree(treestring=xml)
        self.assertEqual(str(orig), str(parsed))
    
    # Magic methods
    
    def test_str(self):
        """testing (well, exercising at least), __str__"""
        str(self.default_tree)
    
    def test_repr(self):
        """testing (well, exercising at least), __repr__"""
        repr(self.default_tree)
    
    def test_eq(self):
        """testing (well, exercising at least), __eq__"""
        # xxx not good enough!
        t1 = self._maketree()
        t2 = self._maketree()
        self.assertTrue(t1 == t1)
        self.assertFalse(t1 == t2)
    
    def test_balanced(self):
        """balancing an unrooted tree"""
        t = LoadTree(treestring='((a,b),((c1,(c2,(c3,(c4,(c5,(c6,c7)))))),(d,e)),f)')
        b = LoadTree(treestring='(c1,(c2,(c3,(c4,(c5,(c6,c7))))),((d,e),((a,b),f)))')
        self.assertEqual(str(t.balanced()), str(b))
    
    def test_params_merge(self):
        t = LoadTree(treestring='((((a,b)ab,c)abc),d)')
        for (label, length, beta) in [('a',1, 20),('b',3,2.0),('ab',4,5.0),]:
            t.getNodeMatchingName(label).params = {'length':length, 'beta':beta}
        t = t.getSubTree(['b', 'c', 'd'])
        self.assertEqual(t.getNodeMatchingName('b').params,
                                {'length':7, 'beta':float(2*3+4*5)/(3+4)})
        self.assertRaises(ValueError, t.getSubTree, ['b','c','xxx'])
        self.assertEqual(str(t.getSubTree(['b','c','xxx'],ignore_missing=True)),
            '(b:7,c)root;')
    
    def test_making_from_list(self):
        tipnames_with_spaces = ['a_b','a b',"T'lk"]
        tipnames_with_spaces.sort()
        t = LoadTree(tip_names=tipnames_with_spaces)
        result = t.getTipNames()
        result.sort()
        assert result == tipnames_with_spaces
    
    def test_getsetParamValue(self):
        """test getting, setting of param values"""
        t = LoadTree(treestring='((((a:.2,b:.3)ab:.1,c:.3)abc:.4),d:.6)')
        self.assertEqual(t.getParamValue('length', 'ab'), 0.1, 2)
        t.setParamValue('zz', 'ab', 4.321)
        node = t.getNodeMatchingName('ab')
        self.assertEqual(4.321, node.params['zz'], 4)

class SmallTreeReshapeTestClass(unittest.TestCase):
    def test_rootswaps(self):
        """testing (well, exercising at least), unrooted"""
        new_tree = LoadTree(treestring="((a,b),(c,d))")
        new_tree = new_tree.unrooted()
        self.assert_(len(new_tree.Children) > 2, 'not unrooted right')
    
    def test_reroot(self):
        tree = LoadTree(treestring="((a,b),(c,d),e)")
        tree2 = tree.rootedWithTip('b')
        self.assertEqual(tree2.getNewick(), "(a,b,((c,d),e));")
    
    def test_sameShape(self):
        """test topology assessment"""
        t1 = LoadTree(treestring="(((s1,s5),s3),s2,s4);")
        t2 = LoadTree(treestring="((s1,s5),(s2,s4),s3);")
        t3 = LoadTree(treestring="((s1,s4),(s2,s5),s3);")
        assert t1.sameTopology(t2), (t1, t2)
        assert not t1.sameTopology(t3), (t1, t3)
        assert not t2.sameTopology(t3), (t2, t3)


#=============================================================================
# these are tests involving tree manipulation methods
# hence, testing them for small and big trees
# the tests are written once for the small tree, the big tree
# tests are performed by inheriting from this class, but over-riding
# the setUp.

class TestTree(unittest.TestCase):
    """tests for a single tree-type"""
    
    def setUp(self):
        self.name = 'small tree - '
        self.otu_names = ['NineBande', 'Mouse', 'HowlerMon', 'DogFaced']
        self.otu_names.sort()
        self.newick = '(((Human,HowlerMon),Mouse),NineBande,DogFaced);'
        self.newick_sorted = '(DogFaced,((HowlerMon,Human),Mouse),NineBande);'
        self.newick_reduced = '((HowlerMon,Mouse),NineBande,DogFaced);'
        self.tree = LoadTree(treestring = self.newick)
    
    def test_sorttree(self):
        """testing (well, exercising at least) treesort"""
        new_tree = self.tree.sorted()
        if hasattr(self, 'newick_sorted'):
            self.assertEqual(
                self.newick_sorted,
                new_tree.getNewick(with_distances=0))
    
    def test_getsubtree(self):
        """testing getting a subtree"""
        subtree = self.tree.unrooted().getSubTree(self.otu_names)
        
        new_tree = LoadTree(treestring = self.newick_reduced).unrooted()
        
        # check we get the same names
        self.assertEqual(*[len(t.Children) for t in (subtree,new_tree)])
        self.assertEqual(str(subtree), str(new_tree))
    
    def test_ascii(self):
        self.tree.asciiArt()
        # unlabeled internal node
        tr = DndParser("(B:0.2,(C:0.3,D:0.4):0.6)F;")
        tr.asciiArt(show_internal=True, compact=False)
        tr.asciiArt(show_internal=True, compact=True)
        tr.asciiArt(show_internal=False, compact=False)

# the following class repeats the above tests but using a big tree and big data-set
class BigTreeSingleTests(TestTree):
    """using the big-tree for single-tree tests"""
    def setUp(self):
        self.name = 'big tree - '
        self.otu_names = ['Horse', 'TombBat', 'Rhino', 'Pig', 'AsianElep',
                     'SpermWhal', 'Cat', 'Gorilla', 'Orangutan',
                     'bandicoot', 'Hedgehog', 'Sloth', 'HairyArma',
                     'Manatee', 'GoldenMol', 'Pangolin']
        self.otu_names.sort()
        self.newick = '((((((((FlyingFox,DogFaced),((FreeTaile,LittleBro),(TombBat,RoundEare))),(FalseVamp,LeafNose)),(((Horse,Rhino),(Pangolin,(Cat,Dog))),(Llama,(Pig,(Cow,(Hippo,(SpermWhal,HumpbackW))))))),(Mole,Hedgehog)),(TreeShrew,(FlyingLem,((Jackrabbit,(FlyingSqu,(OldWorld,(Mouse,Rat)))),(Galago,(HowlerMon,(Rhesus,(Orangutan,(Gorilla,(Human,Chimpanzee)))))))))),(((NineBande,HairyArma),(Anteater,Sloth)),(((Dugong,Manatee),((AfricanEl,AsianElep),(RockHyrax,TreeHyrax))),(Aardvark,((GoldenMol,(Madagascar,Tenrec)),(LesserEle,GiantElep)))))),(caenolest,(phascogale,(wombat,bandicoot))));'
        self.newick_reduced = '(((((TombBat,(((Horse,Rhino),(Pangolin,Cat)),(Pig,SpermWhal))),Hedgehog),(Orangutan,Gorilla)),((HairyArma,Sloth),((Manatee,AsianElep),GoldenMol))),bandicoot);'
        self.tree = LoadTree(treestring = self.newick)
    
    def test_getEdgeNames(self):
        """testing (well, exercising at least), getedgenames"""
        # Fell over on small tree because "stem descended from root
        # joiner was a tip"
        a,b = self.otu_names[:2]
        clade = self.tree.getEdgeNames(a, b, True, False)
    
    def test_getTipNames(self):
        """testing (well, exercising at least), getTipNames"""
        a,b = self.otu_names[:2]
        tips = self.tree.getTipNames()
        self.assertEqual(len(tips), 55)

if __name__ == '__main__':
    unittest.main()
