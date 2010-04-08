#!/usr/bin/env python
"""Tests of data retrieval from NCBI."""
from cogent.util.unit_test import TestCase, main
from cogent.db.ncbi import EUtils, ESearch, EFetch, ELink, ESearchResultParser,\
    ELinkResultParser, get_primary_ids, ids_to_taxon_ids, \
    taxon_lineage_extractor, taxon_ids_to_lineages, taxon_ids_to_names, \
    taxon_ids_to_names_and_lineages, \
    get_unique_lineages, get_unique_taxa
from string import strip

__author__ = "Mike Robeson"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Mike Robeson", "Rob Knight", "Zongzhi Liu"]
__license__ = "GPL"
__version__ = "1.4.1"
__maintainer__ = "Mike Robeson"
__email__ = "mike.robeson@colorado.edu"
__status__ = "Production"

class EUtilsTests(TestCase):
    """Tests of the EUtils class."""
    def test_simple_get(self):
        """EUtils simple access of an item should work"""
        g = EUtils(db='protein',rettype='gp')
        result = g['NP_003320'].read()
        assert result.startswith('LOCUS')
        assert 'NP_003320' in result

    def test_get_slice(self):
        """EUtils access of a slice should work"""
        g = EUtils(db='protein',rettype='gp', retmax=1)
        result = g['NP_003320':'NP_003322'].read()
        lines = result.splitlines()
        is_locus = lambda x: x.startswith('LOCUS')
        loci = filter(is_locus, lines)
        self.assertEqual(len(loci), 3)

    def test_get_list(self):
        """EUtils access of a list should work"""
        g = EUtils(db='protein',rettype='gp')
        result = g['NP_003320','NP_003321','NP_003322'].read()
        lines = result.splitlines()
        is_locus = lambda x: x.startswith('LOCUS')
        loci = filter(is_locus, lines)
        self.assertEqual(len(loci), 3)

    def test_get_from_taxonomy_db(self):
        """EUtils access from taxonomy database should work"""
        #note: this is more fragile than the nucleotide databases
        g = EUtils(db='taxonomy', rettype='Brief', retmode='text')
        ids = '9606[taxid] OR 28901[taxid]'
        result = sorted(g[ids].read().splitlines())
        self.assertEqual(result, ['Homo sapiens', 'Salmonella enterica'])
        
    def test_query(self):
        """EUtils access via a query should work"""
        g = EUtils(db='protein', rettype='gi', retmax=100)
        result = g['homo[organism] AND erf1[ti]'].read().splitlines()
        assert '5499721' in result  #gi of human eRF1
        #note: successfully retrieved 841,821 ids on a query for 'rrna',
        #although it took about 40 min so not in the tests. RK 1/3/07.

    def test_query_retmax(self):
        """EUtils should join results taken retmax at a time"""
        g = EUtils(db='protein', rettype='gi', retmax=3, DEBUG=False)
        result = g['homo[organism] AND myh7'].read().splitlines()
        assert len(result) > 1
        assert '83304912' in result  #gi of human myh7

    def test_query_max_recs(self):
        """EUtils should stop query at max_recs when max_recs < retmax"""
        g = EUtils(db='protein', rettype='gi', max_recs=5, DEBUG=False,
            retmax=100)
        result = g['homo[organism] AND myh7'].read().splitlines()
        self.assertEqual(len(result), 5)

    def test_query_max_recs_gt_retmax(self):
        """EUtils should stop query at max_recs when max_recs > retmax"""
        g = EUtils(db='protein', rettype='gi', max_recs=5, DEBUG=False,
            retmax=3)
        result = g['homo[organism] AND myh7'].read().splitlines()
        self.assertEqual(len(result), 5)


class ESearchTests(TestCase):
    """Tests of the ESearch class: gets primary ids from search."""
    def test_simple_search(self):
        """ESearch Access via a query should return accessions"""
        s = ESearch(db='protein', rettype='gi', retmax=1000,
            term='homo[organism] AND myh7')
        result = s.read()
        parsed = ESearchResultParser(result)
        assert '83304912' in parsed.IdList  #gi of human cardiac beta myh7

class ELinkTests(TestCase):
    """Tests of the ELink class: converts ids between databases"""
    def test_simple_elink(self):
        """ELink should retrieve a link from a single id"""
        l = ELink(db='taxonomy', dbfrom='protein', id='83304912')
        result = l.read()
        parsed = ELinkResultParser(result)
        self.assertEqual(parsed, ['9606'])  #human sequence

    def test_multiple_elink(self):
        """ELink should find unique links in a set of ids"""
        l = ELink(db='taxonomy', dbfrom='protein', 
            id='83304912 115496169 119586556 111309484')
        result = l.read()
        parsed = ELinkResultParser(result)
        self.assertEqual(sorted(parsed), ['10090', '9606'])  
        #human and mouse sequences

class EFetchTests(TestCase):
    """Tests of the EFetch class: gets records using primary ids."""
    def test_simple_efetch(self):
        """EFetch should return records from list of ids"""
        f = EFetch(db='protein', rettype='fasta', retmode='text', 
            id='111309484')
        result = f.read().splitlines()
        assert result[0].startswith('>')
        assert result[1].startswith('madaemaafg'.upper())

class NcbiTests(TestCase):
    """Tests of top-level convenience wrappers."""
    def setUp(self):
        """Define some lengthy data."""
        self.mouse_taxonomy = map(strip, 'cellular organisms; Eukaryota; Fungi/Metazoa group; Metazoa; Eumetazoa; Bilateria; Coelomata; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Euarchontoglires; Glires; Rodentia; Sciurognathi; Muroidea; Muridae; Murinae; Mus'.split(';'))
        self.human_taxonomy = map(strip, 'cellular organisms; Eukaryota; Fungi/Metazoa group; Metazoa; Eumetazoa; Bilateria; Coelomata; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Euarchontoglires; Primates; Haplorrhini; Simiiformes; Catarrhini; Hominoidea; Hominidae; Homininae; Homo'.split(';'))

    def test_get_primary_ids(self):
        """get_primary_ids should get primary ids from query"""
        res = get_primary_ids('homo[orgn] AND myh7[ti]', retmax=5, max_recs=7)
        self.assertEqual(len(res), 7)
        res = get_primary_ids('homo[orgn] AND myh7[ti]', retmax=5, max_recs=2)
        self.assertEqual(len(res), 2)
        res = get_primary_ids('homo[orgn] AND myh7[ti]', retmax=100)
        assert '115496168' in res

    def test_ids_to_taxon_ids(self):
        """ids_to_taxonomy should get taxon ids from primary ids"""
        ids = ['83304912', '115496169', '119586556', '111309484']
        result = ids_to_taxon_ids(ids, db='protein')
        self.assertEqual(sorted(result), ['10090', '9606'])

    def test_taxon_lineage_extractor(self):
        """taxon_lineage_extractor should find lineage lines"""
        lines = """ignore
        <Lineage>xxx;yyy</Lineage>
        ignore

        <Lineage>aaa;bbb</Lineage>
        """
        self.assertEqual(list(taxon_lineage_extractor(lines.splitlines())),
            [['xxx','yyy'],['aaa','bbb']])

    def test_taxon_ids_to_lineages(self):
        """taxon_ids_to_lineages should return lineages from taxon ids"""
        taxon_ids = ['10090', '9606']
        result = [self.mouse_taxonomy, self.human_taxonomy]
        self.assertEqualItems(list(taxon_ids_to_lineages(taxon_ids)), result)
    
    def test_taxon_ids_to_names(self):
        """taxon_ids_to_names should return names from taxon ids"""
        taxon_ids = ['10090', '9606']
        result = set(['Mus musculus', 'Homo sapiens'])
        self.assertEqual(set(taxon_ids_to_names(taxon_ids)), result)

    def test_taxon_ids_to_names_and_lineages(self):
        """taxon_ids_to_names should return names/lineages from taxon ids"""
        taxon_ids = ['10090', '9606']
        exp = [('10090', 'Mus musculus', self.mouse_taxonomy),
                  ('9606', 'Homo sapiens', self.human_taxonomy)]
        obs = list(taxon_ids_to_names_and_lineages(taxon_ids))
        self.assertEqualItems(obs, exp)

    def test_get_unique_lineages(self):
        """get_unique_lineages should return all lineages from a query"""
        result = get_unique_lineages('angiotensin[ti] AND rodents[orgn]')
        
        assert tuple(self.mouse_taxonomy) in result
        assert len(result) > 2

    def test_get_unique_taxa(self):
        """get_unique_taxa should return all taxa from a query"""
        result = get_unique_taxa('angiotensin[ti] AND primate[orgn]')
        assert 'Homo sapiens' in result
        assert len(result) > 2

if __name__ == '__main__':
    main()
