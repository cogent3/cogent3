from unittest import TestCase, main

from numpy.testing import assert_allclose

from cogent3 import make_tree
from cogent3.draw.dendrogram import (
    CircularTreeGeometry,
    Dendrogram,
    SquareTreeGeometry,
)
from cogent3.util.deserialise import deserialise_object
from cogent3.util.union_dict import UnionDict


__author__ = "Gavin Huttley and Rahul Ghangas"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Rahul Ghangas"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class TestDendro(TestCase):
    def test_geometry(self):
        """tree geometry class get_edge_names works"""
        tree = make_tree(treestring="(a,b,(c,(d,e)e1)e2)")
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
        tree = make_tree(treestring="(a,b,(c,(d,e)e1)e2)")
        dnd = Dendrogram(tree=tree)
        edges = dnd.get_edge_names("d", "c", clade=True, stem=False)
        self.assertEqual(set(edges), {"c", "d", "e", "e1"})

    def test_min_max_x_y(self):
        """correctly compute the min and max of x and y"""
        tree = make_tree(treestring="(A:1,B:2,C:3)")
        geom = CircularTreeGeometry(tree)
        geom.propagate_properties()
        got = max(map(abs, [geom.min_x, geom.max_x]))
        expect = max(map(abs, [e.x for e in geom.postorder()]))
        assert_allclose(got, expect)
        got = max(map(abs, [geom.min_y, geom.max_y]))
        expect = max(map(abs, [e.y for e in geom.postorder()]))
        assert_allclose(got, expect)

    def test_length_attr_valid(self):
        """Tests whether setting a custom length attribute provides valid x values"""
        tree = make_tree(
            treestring="((a:0.1,b:0.25):0.1,(c:0.02,d:0.03, (e:0.035, f:0.04):0.15):0.3 , g:0.3)"
        )
        geom = SquareTreeGeometry(tree, length_attr="custom")
        geom.params["custom"] = 1

        for e in geom.preorder():
            if e.is_root():
                continue
            e.params["custom"] = e.parent.params.get("custom", 1) * 2
        geom.propagate_properties()

        # .x attribute is cumulative from the root, which we have set to 1
        # for 'custom', e.g. a.x == 2 + 4 == 6
        func = geom.get_node_matching_name
        actual_vals = [
            func("root").x,
            func("a").x,
            func("b").x,
            func("c").x,
            func("d").x,
            func("e").x,
            func("f").x,
            func("g").x,
        ]

        expected_vals = [0, 6, 6, 6, 6, 14, 14, 2]

        # Root x resets to 0 so any assigned value to root is always discarded

        assert_allclose(actual_vals, expected_vals)

    def test_square_dendrogram_regression(self):
        tree = make_tree(treestring="(a:0.1,b:0.1,(c:0.05,(d:0.01,e:0.02):0.01):0.1)")
        dendrogram = Dendrogram(tree, style="square", contemporaneous=False)
        func = dendrogram.tree.get_node_matching_name
        actual_vals = [
            (func("root").x, func("root").y),
            (func("a").x, func("a").y),
            (func("b").x, func("b").y),
            (func("c").x, func("c").y),
            (func("d").x, func("d").y),
            (func("e").x, func("e").y),
        ]

        expected_vals = [
            (0, 1.3),
            (0.1, 2.6),
            (0.1, 1.3),
            (0.15, -2.6),
            (0.12, 0),
            (0.13, -1.3),
        ]

        assert_allclose(actual_vals, expected_vals)

    def test_dendro_shape(self):
        """exercising using different values of shape parameter"""
        tree = make_tree(treestring="(a:0.1,b:0.1,(c:0.05,(d:0.01,e:0.02):0.01):0.1)")
        for style in ("square", "angular", "circular", "radial"):
            dnd = Dendrogram(tree, style=style)
            # the figure attribute should be a dict
            fig = dnd.figure
            self.assertIsInstance(fig, UnionDict)
            # should have a layout and a data key
            self.assertTrue("layout" in fig)
            self.assertTrue("data" in fig)
            # data traces should be of type "scatter"
            self.assertEqual({tr.type for tr in fig.data}, {"scatter"})

    def test_dendro_with_support(self):
        """exercising creating dendrograms with support measure"""
        data = {
            "newick": "(A,(B,C)edge.1,(D,E)edge.0)",
            "edge_attributes": {
                "A": {"support": 1.0, "length": 0.148},
                "B": {"support": 1.0, "length": 0.098},
                "C": {"support": 1.0, "length": 0.134},
                "edge.1": {"support": 0.8, "length": 0.016},
                "D": {"support": 1.0, "length": 0.087},
                "E": {"support": 1.0, "length": 0.048},
                "edge.0": {"support": 0.6, "length": 0.131},
                "root": {"length": None},
            },
            "type": "cogent3.core.tree.PhyloNode",
            "version": "2019.10.17a",
        }
        tree = deserialise_object(data)
        # currently not supported for "circular", "radial"
        for style in ("square", "angular"):
            y_support = Dendrogram(tree, style=style, show_support=True, threshold=0.8)
            n_support = Dendrogram(tree, style=style, show_support=False)
            self.assertEqual(
                len(y_support.figure.layout.annotations)
                - len(n_support.figure.layout.annotations),
                2,
                style,
            )

    def test_style_edges(self):
        """test style_edges only accepts edges present in tree"""
        tree = make_tree(treestring="(a,b,(c,(d,e)e1)e2)")
        dnd = Dendrogram(tree=tree)
        dnd.style_edges("a", line=dict(color="magenta"))
        with self.assertRaises(ValueError):
            dnd.style_edges("foo", line=dict(color="magenta"))

    def test_tip_font(self):
        """test tip_font settable"""
        tree = make_tree(treestring="(a,b,(c,(d,e)e1)e2)")
        dnd = Dendrogram(tree=tree)
        dnd.tip_font |= dict(size=18)
        self.assertEqual(dnd.tip_font.size, 18)
        dnd.tip_font.size = 10
        self.assertEqual(dnd.tip_font.size, 10)
        dnd.tip_font.color = "red"
        self.assertEqual(dnd.tip_font["color"], "red")


if __name__ == "__main__":
    main()
