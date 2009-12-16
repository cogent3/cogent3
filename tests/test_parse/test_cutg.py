#!/usr/bin/env python
"""Unit tests for the CUTG database parsers.
"""
from cogent.parse.cutg import CutgParser, CutgSpeciesParser, InfoFromLabel
from cogent.parse.record import RecordError
from cogent.util.unit_test import TestCase, main

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

sample_gene = r""">AB000406\AB000406\100..1506\1407\BAA19100.1\Xenopus laevis\Xenopus laevis mRNA for protein phosphatase 2A regulatory subunit,complete cds./gene="PP2A"/codon_start=1/product="protein phosphatase 2A regulatory subunit"/protein_id="BAA19100.1"/db_xref="GI:1783183"
1 7 4 2 15 6 4 4 11 4 1 7 6 7 4 7 14 7 10 7 4 6 5 2 3 8 3 4 3 6 5 4 6 0 8 5 8 8 20 10 21 13 2 9 5 7 22 16 20 10 10 8 7 3 6 15 12 5 14 11 6 0 1 0 
>AB000458#1\AB000458\105..623\519\BAA22881.1\Xenopus laevis\Xenopus laevis Xem1 mRNA for transmembrane protein, complete cds./gene="Xem1"/codon_start=1/product="transmembrane protein"/protein_id="BAA22881.1"/db_xref="GI:2554596"
0 1 1 0 2 1 6 7 10 4 0 1 4 6 1 5 2 1 3 1 1 5 2 5 0 5 2 4 0 2 2 1 1 1 2 5 4 4 2 2 1 3 1 6 4 2 2 5 1 2 3 2 1 2 2 9 3 5 2 6 4 1 0 0 
>AB000736\AB000736\27..557\531\BAA19174.1\Xenopus laevis\Xenopus laevis mRNA for myelin basic protein, complete cds./codon_start=1/product="myelin basic protein"/protein_id="BAA19174.1"/db_xref="GI:1816437"
1 1 1 2 5 11 1 0 3 2 0 1 8 4 2 11 3 1 2 1 0 0 5 2 0 5 3 3 0 1 14 3 4 2 1 0 1 4 3 7 3 0 2 3 9 7 1 4 3 3 5 4 0 0 5 2 1 1 0 4 1 0 0 1 """.split('\n')

sample_species = r"""Salmonella enterica: 332
432 1587 640 1808 586 275 758 825 3557 1668 2943 1663 1267 928 1002 1397 1359 1277 1096 1880 1305 1267 911 608 1752 1002 1684 1868 2960 1627 995 2399 1228 2103 1605 1289 2024 2017 4800 1681 1719 3181 1699 2857 700 1723 4347 2560 1577 4174 1031 2938 496 678 1354 3501 1713 1955 3743 2497 1366 214 22 96 
Salmonella enterica IIIb 50:k:z: 3
2 5 2 4 3 3 4 6 12 5 7 1 7 5 1 4 7 3 10 13 5 6 9 4 3 5 17 14 14 14 6 18 5 12 15 10 4 11 24 5 16 17 9 19 6 3 11 9 9 18 7 20 2 0 8 17 6 11 16 11 5 2 0 1 
Salmonella enterica subsp. VII: 5
4 1 6 7 6 0 8 15 10 9 18 9 15 6 6 19 21 3 21 14 11 22 16 10 11 6 38 13 20 10 8 8 9 9 7 5 7 22 58 21 25 40 24 22 11 16 41 21 6 14 9 11 2 12 10 31 16 23 22 22 3 0 1 4""".split('\n')

strange_db = r'''>AB001737\AB001737\1..696\696\BAA19944.1\Mus musculus\Mus musculus mRNA for anti-CEA scFv antibody, complete cds./codon_start=1/product="anti-CEA scFv antibody"/protein_id="BAA19944.1"/db_xref="GI:2094751"/db_xref="IMGT/LIGM:AB001737"'''

class InfoFromLabelTests(TestCase):
    """Tests of the InfoFromLabel constructor."""
    def test_init(self):
        """InfoFromLabel should handle a typical label line"""
        i = InfoFromLabel(sample_gene[0])
        sa = self.assertEqual
        sa(i.GenBank, ['AB000406'])
        sa(i.Locus, 'AB000406')
        sa(i.CdsNumber, '1'),
        sa(i.Location, '100..1506')
        sa(i.Length, '1407')
        sa(i.Species, 'Xenopus laevis')
        sa(i.Description, r'Xenopus laevis mRNA for protein phosphatase 2A regulatory subunit,complete cds./gene="PP2A"/codon_start=1/product="protein phosphatase 2A regulatory subunit"/protein_id="BAA19100.1"/db_xref="GI:1783183"')
        sa(i.Gene, 'PP2A')
        sa(i.CodonStart, '1')
        sa(i.Product, 'protein phosphatase 2A regulatory subunit')
        sa(i.GenPept, ['BAA19100.1'])
        sa(i.GI, ['1783183'])

        j = InfoFromLabel(sample_gene[2])
        assert j.Refs is not i.Refs
        assert j._handler is not i._handler
        assert j._handler is j.Refs
        assert j.Refs.GI is not i.Refs.GI
        assert j.GI is not i.GI
        sa(j.GenBank, ['AB000458']),
        sa(j.Locus, 'AB000458')
        sa(j.CdsNumber, '1')
        sa(j.Location, '105..623')
        sa(j.Length, '519'),
        sa(j.Species, 'Xenopus laevis')
        sa(j.Description, 'Xenopus laevis Xem1 mRNA for transmembrane protein, complete cds./gene="Xem1"/codon_start=1/product="transmembrane protein"/protein_id="BAA22881.1"/db_xref="GI:2554596"')
        sa(j.GenPept, ['BAA22881.1']),
        sa(j.GI, ['2554596'])
        sa(j.Product, 'transmembrane protein')

    def test_init_unknown_db(self):
        """InfoFromLabel should handle a line whose database is unknown"""
        i = InfoFromLabel(strange_db)
        self.assertEqual(i.Locus, 'AB001737')
   
class CutgSpeciesParserTests(TestCase):
    """Tests of the CutgSpeciesParser."""
    def test_init(self):
        """CutgSpeciesParser should read records one at a time from lines"""
        recs = list(CutgSpeciesParser(sample_species))
        self.assertEqual(len(recs), 3)
        a, b, c = recs
        self.assertEqual(a.Species, 'Salmonella enterica')
        self.assertEqual(a.NumGenes, 332)
        self.assertEqual(a['CGA'], 432)
        self.assertEqual(a['UGG'], 1366)
        self.assertEqual(b.Species, 'Salmonella enterica IIIb 50:k:z')
        self.assertEqual(b.NumGenes, 3)
        self.assertEqual(b['CGA'], 2)
        self.assertEqual(b['UGG'], 5)
        self.assertEqual(c.Species, 'Salmonella enterica subsp. VII')
        self.assertEqual(c.NumGenes, 5)
        self.assertEqual(c['CGA'], 4)
        self.assertEqual(c['UGG'], 3)
        #check that it won't work if we're missing any lines
        self.assertRaises(RecordError, list, 
            CutgSpeciesParser(sample_species[1:]))
        self.assertRaises(RecordError, list,
            CutgSpeciesParser(sample_species[:-1]))
        #...but that it does work if we only have some of them
        recs = list(CutgSpeciesParser(sample_species[2:]))
        self.assertEqual(recs[0], b)
        self.assertEqual(len(list(CutgSpeciesParser(sample_species[1:],
            strict=False))), 2)

class CutgParserTests(TestCase):
    """Tests of the CutgParser.
    
    Note: these are fairly incomplete at present since most of the work is in
    parsing the label line, which is tested by itself.
    """
    def test_init(self):
        """CutgParser should read records one at a time from lines"""
        recs = list(CutgParser(sample_gene))
        self.assertEqual(len(recs), 3)
        a, b, c = recs
        self.assertEqual(a.Species, 'Xenopus laevis')
        self.assertEqual(a['CGC'], 7)
        self.assertEqual(a.GI, ['1783183'])
        self.assertRaises(RecordError, list, CutgParser(sample_gene[1:]))
        self.assertEqual(len(list(CutgParser(sample_gene[1:],strict=False))), 2)
        
if __name__ == '__main__':
    main()
