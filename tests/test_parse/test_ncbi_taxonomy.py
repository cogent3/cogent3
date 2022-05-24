#!/usr/bin/env python
"""Tests of parsers for dealing with NCBI Taxonomy files.
"""

from unittest import TestCase, main

from cogent3.parse.ncbi_taxonomy import (
    MissingParentError,
    NcbiName,
    NcbiNameLookup,
    NcbiNameParser,
    NcbiTaxon,
    NcbiTaxonLookup,
    NcbiTaxonomyFromFiles,
    NcbiTaxonParser,
)


__author__ = "Jason Carnes"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Jason Carnes", "Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

good_nodes = """1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
2\t|\t1\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
6\t|\t2\t|\tgenus\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t|\t
7\t|\t6\t|\tspecies\t|\tAC\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
9\t|\t7\t|\tsubspecies\t|\tBA\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
10\t|\t6\t|\tspecies\t|\tAC\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
""".split(
    "\n"
)

bad_nodes = """1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
2\t|\t1\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
6\t|\t2\t|\tgenus\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t|\t
7\t|\t6\t|\tspecies\t|\tAC\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
9\t|\t777\t|\tsubspecies\t|\tBA\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
10\t|\t666\t|\tspecies\t|\tAC\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
""".split(
    "\n"
)

good_names = """1\t|\tall\t|\t\t|\tsynonym\t|
1\t|\troot\t|\t\t|\tscientific name\t|
2\t|\tBacteria\t|\tBacteria <bacteria>\t|\tscientific name\t|
2\t|\tMonera\t|\tMonera <Bacteria>\t|\tin-part\t|
2\t|\tProcaryotae\t|\tProcaryotae <#1>\t|\tin-part\t|
2\t|\tProkaryotae\t|\tProkaryotae <#1>\t|\tin-part\t|
2\t|\teubacteria\t|\t\t|\tgenbank common name\t|
6\t|\tAzorhizobium\t|\t\t|\tscientific name\t|
7\t|\tAzorhizobium caulinodans\t|\t\t|\tscientific name\t|
9\t|\tBuchnera aphidicola\t|\t\t|\tscientific name\t|
10\t|\tFakus namus\t|\t\t|\tscientific name\t|
""".split(
    "\n"
)


class NcbiTaxonTests(TestCase):
    """Tests proper parsing of NCBI node file, e.g. nodes.dmp"""

    def test_init(self):
        """NcbiTaxon init should return object containing taxonomy data"""
        good_1 = """1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"""
        good_2 = """2\t|\t1\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"""
        good_3 = """6\t|\t2\t|\tgenus\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t|\n"""
        good_4 = """7\t|\t6\t|\tspecies\t|\tAC\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\n"""
        node_1 = NcbiTaxon(good_1)  # make a NcbiTaxon object
        node_2 = NcbiTaxon(good_2)  # from the corresponding
        node_3 = NcbiTaxon(good_3)  # line.
        node_4 = NcbiTaxon(good_4)
        self.assertEqual(node_1.Rank, "no rank")  # confirm object holds
        self.assertEqual(node_1.RankId, 28)  # right data
        self.assertEqual(node_1.ParentId, 1)
        self.assertEqual(node_2.Rank, "superkingdom")
        self.assertEqual(node_2.RankId, 27)
        self.assertEqual(node_2.ParentId, 1)
        self.assertEqual(node_3.Rank, "genus")
        self.assertEqual(node_3.RankId, 8)
        self.assertEqual(node_3.ParentId, 2)
        self.assertEqual(node_4.Rank, "species")
        self.assertEqual(node_4.RankId, 4)
        self.assertEqual(node_4.ParentId, 6)
        # test some comparisons
        assert node_1 > node_2
        assert node_1 > node_3
        assert node_1 > node_4
        assert node_1 == node_1
        assert node_2 < node_1
        assert node_2 == node_2
        assert node_4 < node_1
        assert node_3 > node_4

    def test_str(self):
        """NcbiTaxon str should write data in input format from nodes.dmp"""
        good = """2\t|\t1\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"""
        node = NcbiTaxon(good)
        self.assertEqual(str(node), good)
        root = """1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|"""
        NcbiTaxon(root)
        self.assertEqual(str(root), root)

    def test_bad_input(self):
        """NcbiTaxon init should raise ValueError if nodes missing"""
        bad_node_taxid = """\t|\t6\t|\tspecies\t|\tAC\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\n"""  # contains no taxon_id; not valid
        bad_node_parentid = """7\t|\t\t|\tspecies\t|\tAC\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\n"""  # contains no parent_id; not valid
        self.assertRaises(ValueError, NcbiTaxon, bad_node_taxid)
        self.assertRaises(ValueError, NcbiTaxon, bad_node_parentid)


class NcbiNameTests(TestCase):
    """Tests proper parsing NCBI name file, e.g. names.dmp."""

    def test_init(self):
        """NcbiName should init OK with well-formed name line"""
        line_1 = """1\t|\tall\t|\t\t|\tsynonym\t|\n"""
        line_2 = """1\t|\troot\t|\t\t|\tscientific name\t|\n"""
        line_3 = """2\t|\tBacteria\t|\tBacteria <bacteria>\t|\tscientific name\t|\n"""
        line_4 = """7\t|\tAzorhizobium caulinodans\t|\t\t|\tscientific name\t|\n"""
        name_1 = NcbiName(line_1)  # make an NcbiName object
        name_2 = NcbiName(line_2)  # from the corresponding line
        name_3 = NcbiName(line_3)
        name_4 = NcbiName(line_4)
        self.assertEqual(name_1.TaxonId, 1)  # test that the data
        self.assertEqual(name_1.NameClass, "synonym")  # fields in the object
        self.assertEqual(name_2.TaxonId, 1)  # hold right data
        self.assertEqual(name_2.NameClass, "scientific name")
        self.assertEqual(name_3.TaxonId, 2)
        self.assertEqual(name_3.NameClass, "scientific name")
        self.assertEqual(name_4.TaxonId, 7)
        self.assertEqual(name_4.NameClass, "scientific name")

    def test_str(self):
        """NcbiName str should return line in original format"""
        line = """1\t|\troot\t|\t\t|\tscientific name|\n"""
        name = NcbiName(line)
        self.assertEqual(str(name), line)

    def test_bad_input(self):
        """NcbiName init should raise correct errors on bad data"""
        bad_name_taxid = """\t|\troot\t|\t\t|\tscientific name\t|\n"""  # no tax_id
        self.assertRaises(ValueError, NcbiName, bad_name_taxid)


class NcbiNameLookupTest(TestCase):
    """Tests of the NcbiNameLookup factory function."""

    def test_init(self):
        """NcbiNameLookup should map taxon ids to scientific names"""
        names = list(NcbiNameParser(good_names))  # list of objects
        sci_names = NcbiNameLookup(names)  # NcbiNameLookup object
        root = names[1]  # NcbiName object made from 2nd line of good_name_file
        bacteria = names[2]  # from 3rd line of good_name_file
        azorhizobium = names[7]
        caulinodans = names[8]
        assert sci_names[1] is root  # gets NcbiName object from the
        assert sci_names[2] is bacteria  # NcbiNameLookup object and
        assert sci_names[6] is azorhizobium  # asks if it is the original
        assert sci_names[7] is caulinodans  # NcbiName object
        self.assertEqual(sci_names[1].Name, "root")
        self.assertEqual(sci_names[2].Name, "Bacteria")
        self.assertEqual(sci_names[7].Name, "Azorhizobium caulinodans")
        self.assertEqual(sci_names[9].Name, "Buchnera aphidicola")


class NcbiTaxonLookupTest(TestCase):
    """Tests of the NcbiTaxonLookup factory function."""

    def setUp(self):
        """Sets up the class tests"""
        self.names = list(NcbiNameParser(good_names))
        self.nodes = list(NcbiTaxonParser(good_nodes))
        self.taxID_to_obj = NcbiTaxonLookup(self.nodes)
        self.names_to_obj = NcbiNameLookup(self.names)

    def test_init(self):
        """NcbiTaxonLookup should have correct fields for input NcbiTaxon"""
        line1_obj = self.nodes[0]  # NcbiTaxon objects made from lines of
        line2_obj = self.nodes[1]  # good_node_file
        line3_obj = self.nodes[2]
        line4_obj = self.nodes[3]
        line5_obj = self.nodes[4]
        # gets NcbiTaxon object from
        assert self.taxID_to_obj[1] is line1_obj
        assert self.taxID_to_obj[2] is line2_obj  # NcbiTaxonLookup object &
        # asks if it is the original
        assert self.taxID_to_obj[6] is line3_obj
        assert self.taxID_to_obj[7] is line4_obj  # NcbiTaxon object
        assert self.taxID_to_obj[9] is line5_obj
        self.assertEqual(self.taxID_to_obj[1].ParentId, 1)  # checking a few
        self.assertEqual(self.taxID_to_obj[2].ParentId, 1)  # individual
        self.assertEqual(self.taxID_to_obj[6].ParentId, 2)  # fields of the
        self.assertEqual(self.taxID_to_obj[7].ParentId, 6)  # NcbiTaxon objs
        self.assertEqual(self.taxID_to_obj[9].ParentId, 7)
        self.assertEqual(self.taxID_to_obj[1].Rank, "no rank")
        self.assertEqual(self.taxID_to_obj[2].Rank, "superkingdom")
        self.assertEqual(self.taxID_to_obj[6].Rank, "genus")
        self.assertEqual(self.taxID_to_obj[7].Rank, "species")
        self.assertEqual(self.taxID_to_obj[7].EmblCode, "AC")
        self.assertEqual(self.taxID_to_obj[7].DivisionId, "0")
        self.assertEqual(self.taxID_to_obj[7].DivisionInherited, 1)
        self.assertEqual(self.taxID_to_obj[7].TranslTable, 11)
        self.assertEqual(self.taxID_to_obj[7].TranslTableInherited, 1)
        self.assertEqual(self.taxID_to_obj[7].TranslTableMt, 0)
        self.assertEqual(self.taxID_to_obj[7].TranslTableMtInherited, 1)


class NcbiTaxonomyTests(TestCase):
    """Tests of the NcbiTaxonomy class."""

    def setUp(self):
        self.tx = NcbiTaxonomyFromFiles(good_nodes, good_names)

    def test_init_good(self):
        """NcbiTaxonomyFromFiles should pass spot-checks of resulting objects"""
        self.assertEqual(len(self.tx.ByName), 6)
        self.assertEqual(len(self.tx.ById), 6)
        self.assertEqual(self.tx[10].Name, "Fakus namus")
        self.assertEqual(self.tx["1"].Name, "root")
        self.assertEqual(self.tx["root"].parent, None)
        self.assertEqual(self.tx.Deadbeats, {})

    def test_init_bad(self):
        """NcbiTaxonomyFromFiles should produce deadbeats by default"""
        bad_tx = NcbiTaxonomyFromFiles(bad_nodes, good_names)
        self.assertEqual(len(bad_tx.Deadbeats), 2)
        assert 777 in bad_tx.Deadbeats
        assert 666 in bad_tx.Deadbeats
        assert bad_tx.Deadbeats[777] == bad_tx[9]

    def test_init_strict(self):
        """NcbiTaxonomyFromFiles should fail if strict and deadbeats exist"""
        tx = NcbiTaxonomyFromFiles(good_nodes, good_names, strict=True)
        self.assertRaises(
            MissingParentError,
            NcbiTaxonomyFromFiles,
            bad_nodes,
            good_names,
            strict=True,
        )

    def test_Ancestors(self):
        """NcbiTaxonomy should support Ancestors correctly, not incl. self"""
        result = self.tx["7"].ancestors()
        tax_ids = [taxon_obj.TaxonId for taxon_obj in result]
        self.assertEqual(tax_ids, [6, 2, 1])

    def test_Parent(self):
        """NcbiTaxonomy should support parent correctly"""
        assert self.tx[7].parent is self.tx[6]
        assert self.tx[6].parent is self.tx[2]
        assert self.tx[2].parent is self.tx[1]
        assert self.tx[1].parent is None

    def test_Siblings(self):
        """NcbiTaxonomy should support Siblings correctly"""
        sibs = self.tx[7].siblings()
        self.assertEqual(len(sibs), 1)
        assert sibs[0] is self.tx[10]

    def test_Children(self):
        """NcbiTaxonomy should support children correctly"""
        children = self.tx[6].children
        self.assertEqual(len(children), 2)
        assert children[0] is self.tx[7]
        assert children[1] is self.tx[10]
        root_kids = self.tx["root"]
        self.assertEqual(len(root_kids), 1)
        assert root_kids[0] is self.tx[2]
        self.assertEqual(len(self.tx[10].children), 0)

    def test_Names(self):
        """NcbiTaxonomy should fill in names correctly"""
        self.assertEqual(self.tx["6"].Name, "Azorhizobium")
        self.assertEqual(self.tx["1"].Name, "root")
        self.assertEqual(self.tx["2"].Name, "Bacteria")
        self.assertEqual(self.tx["7"].Name, "Azorhizobium caulinodans")

    def test_last_common_ancestor(self):
        """NcbiTaxonomy should support last_common_ancestor()"""
        assert self.tx[9].last_common_ancestor(self.tx[9]) is self.tx[9]
        assert self.tx[9].last_common_ancestor(self.tx[7]) is self.tx[7]
        assert self.tx[9].last_common_ancestor(self.tx[10]) is self.tx[6]
        assert self.tx[9].last_common_ancestor(self.tx[1]) is self.tx[1]


class NcbiTaxonNodeTests(TestCase):
    """Tests of the NcbiTaxonNode class.

    Note: only testing methods that differ from the TreeNode base class.

    Note: nested_species is explicitly designed to test the case where the nodes
    file does _not_ contain the root, and where the id of the de facto
    root is not 1, to make sure there's nothing special about a node
    called 'root' or with id 1.
    """

    def test_getRankedDescendants(self):
        """NcbiTaxonNode getRankedDescendants should return correct list"""
        nested_species = """3\t|\t3\t|\tsuperkingdom\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
        11\t|\t3\t|\tkingdom\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
        22\t|\t11\t|\tclass\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
        44\t|\t22\t|\torder\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
        66\t|\t22\t|\torder\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t|\t
        77\t|\t66\t|\tfamily\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
        99\t|\t66\t|\tfamily\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
        88\t|\t44\t|\tfamily\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
        101\t|\t77\t|\tgenus\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
        202\t|\t77\t|\tgenus\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
        606\t|\t99\t|\tgenus\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t|\t
        707\t|\t88\t|\tgenus\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
        909\t|\t88\t|\tgenus\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
        123\t|\t909\t|\tgroup\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
        1111\t|\t123\t|\tspecies\t|\tAT\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
        2222\t|\t707\t|\tspecies\t|\tTT\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
        6666\t|\t606\t|\tspecies\t|\tGG\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t|\t
        7777\t|\t606\t|\tspecies\t|\tAC\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
        9999\t|\t202\t|\tspecies\t|\tBA\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
        1010\t|\t101\t|\tspecies\t|\tAC\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
        5555\t|\t555\t|\tspecies\t|\tAC\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
        555\t|\t3\t|\tsuperclass\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|""".split(
            "\n"
        )
        nested_names = [
            "3|a||scientific name|",
            "11|b||scientific name|",
            "555|c||scientific name|",
            "22|d||scientific name|",
            "44|e||scientific name|",
            "66|f||scientific name|",
            "88|g||scientific name|",
            "77|h||scientific name|",
            "99|i||scientific name|",
            "707|j||scientific name|",
            "909|k||scientific name|",
            "101|l||scientific name|",
            "202|m||scientific name|",
            "606|n||scientific name|",
            "2222|o||scientific name|",
            "123|p||scientific name|",
            "1111|q||scientific name|",
            "1010|r||scientific name|",
            "9999|s||scientific name|",
            "7777|t||scientific name|",
            "6666|u||scientific name|",
            "5555|z||scientific name|",
        ]
        tx = NcbiTaxonomyFromFiles(nested_species, nested_names)
        dec = tx[3].getRankedDescendants("superclass")
        self.assertEqual(len(dec), 1)
        assert dec[0] is tx[555]
        sp = tx["f"].getRankedDescendants("species")
        self.assertCountEqual(sp, [tx[1010], tx[9999], tx[7777], tx[6666]])
        empty = tx[11].getRankedDescendants("superclass")
        self.assertEqual(empty, [])
        gr = tx[3].getRankedDescendants("group")
        self.assertEqual(gr, [tx[123]])
        assert tx[3] is tx["a"]


if __name__ == "__main__":
    main()
