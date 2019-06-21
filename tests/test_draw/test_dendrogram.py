from unittest import TestCase, main

from cogent3 import LoadTree
from cogent3.draw.dendrogram import Dendrogram, SquareTreeGeometry


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class TestDendro(TestCase):
    def test_geometry(self):
        """tree geometry class get_edge_names works"""
        tree = LoadTree(treestring="(a,b,(c,(d,e)e1)e2)")
        geom = SquareTreeGeometry(tree)
        series = [
            dict(tip1name="d", tip2name="c", clade=True, stem=False),
            dict(tip1name="d", tip2name="c", clade=True, stem=True),
            dict(tip1name="d", tip2name="c", clade=False, stem=True),
            dict(tip1name="d", tip2name="c", clade=True, stem=False, outgroup_name="e"),
            dict(tip1name="d", tip2name="c", clade=False, stem=True, outgroup_name="e"),
        ]
        for kwargs in series[-1:]:
            expect = tree.get_edge_names(**kwargs)
            got = geom.get_edge_names(**kwargs)
            self.assertEqual(got, expect)

    def test_get_edges(self):
        """returns edge names"""
        tree = LoadTree(treestring="(a,b,(c,(d,e)e1)e2)")
        dnd = Dendrogram(tree=tree)
        edges = dnd.get_edge_names("d", "c", clade=True, stem=False)
        self.assertEqual(set(edges), set(["c", "d", "e", "e1"]))


if __name__ == "__main__":
    main()
