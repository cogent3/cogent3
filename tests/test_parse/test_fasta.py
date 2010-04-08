#!/usr/bin/env python
"""Unit tests for FASTA and related parsers.
"""
from cogent.parse.fasta import FastaParser, MinimalFastaParser, \
    NcbiFastaLabelParser, NcbiFastaParser, RichLabel, LabelParser, GroupFastaParser
from cogent.core.sequence import DnaSequence, Sequence, ProteinSequence as Protein
from cogent.core.info import Info
from cogent.parse.record import RecordError
from cogent.util.unit_test import TestCase, main

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

def Dna(seq, *args, **kwargs):
    seq = seq.replace('u','t')
    seq = seq.replace('U','T')
    d = DnaSequence(seq, *args, **kwargs)
    return d

class GenericFastaTest(TestCase):
    """Setup data for all the various FASTA parsers."""
    def setUp(self):
        """standard files"""
        self.labels = '>abc\n>def\n>ghi\n'.split('\n')
        self.oneseq = '>abc\nUCAG\n'.split('\n')
        self.multiline = '>xyz\nUUUU\nCC\nAAAAA\nG'.split('\n')
        self.threeseq='>123\na\n> \t abc  \t \ncag\ngac\n>456\nc\ng'.split('\n')
        self.twogood='>123\n\n> \t abc  \t \ncag\ngac\n>456\nc\ng'.split('\n')
        self.oneX='>123\nX\n> \t abc  \t \ncag\ngac\n>456\nc\ng'.split('\n')
        self.nolabels = 'GJ>DSJGSJDF\nSFHKLDFS>jkfs\n'.split('\n')
        self.empty = []
 
class MinimalFastaParserTests(GenericFastaTest):
    """Tests of MinimalFastaParser: returns (label, seq) tuples."""
       
    def test_empty(self):
        """MinimalFastaParser should return empty list from 'file' w/o labels"""
        self.assertEqual(list(MinimalFastaParser(self.empty)), [])
        self.assertEqual(list(MinimalFastaParser(self.nolabels, strict=False)),
            [])
        self.assertRaises(RecordError, list, MinimalFastaParser(self.nolabels))

    def test_no_labels(self):
        """MinimalFastaParser should return empty list from file w/o seqs"""
        #should fail if strict (the default)
        self.assertRaises(RecordError, list, 
            MinimalFastaParser(self.labels,strict=True))
        #if not strict, should skip the records
        self.assertEqual(list(MinimalFastaParser(self.labels, strict=False)), 
            [])
        
    def test_single(self):
        """MinimalFastaParser should read single record as (label, seq) tuple"""
        f = list(MinimalFastaParser(self.oneseq))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('abc', 'UCAG'))

        f = list(MinimalFastaParser(self.multiline))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('xyz', 'UUUUCCAAAAAG'))

    def test_multiple(self):
        """MinimalFastaParser should read multiline records correctly"""
        f = list(MinimalFastaParser(self.threeseq))
        self.assertEqual(len(f), 3)
        a, b, c = f
        self.assertEqual(a, ('123', 'a'))
        self.assertEqual(b, ('abc', 'caggac'))
        self.assertEqual(c, ('456', 'cg'))

    def test_multiple_bad(self):
        """MinimalFastaParser should complain or skip bad records"""
        self.assertRaises(RecordError, list, MinimalFastaParser(self.twogood))
        f = list(MinimalFastaParser(self.twogood, strict=False))
        self.assertEqual(len(f), 2)
        a, b = f
        self.assertEqual(a, ('abc', 'caggac'))
        self.assertEqual(b, ('456', 'cg'))

class FastaParserTests(GenericFastaTest):
    """Tests of FastaParser: returns sequence objects."""
       
    def test_empty(self):
        """FastaParser should return empty list from 'file' w/o labels"""
        self.assertEqual(list(FastaParser(self.empty)), [])
        self.assertEqual(list(FastaParser(self.nolabels, strict=False)),
            [])
        self.assertRaises(RecordError, list, FastaParser(self.nolabels))

    def test_no_labels(self):
        """FastaParser should return empty list from file w/o seqs"""
        #should fail if strict (the default)
        self.assertRaises(RecordError, list, 
            FastaParser(self.labels,strict=True))
        #if not strict, should skip the records
        self.assertEqual(list(FastaParser(self.labels, strict=False)), [])
        
    def test_single(self):
        """FastaParser should read single record as seq object"""
        f = list(FastaParser(self.oneseq))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('abc', 'UCAG'))
        self.assertEqual(a[1].Name, 'abc')

        f = list(FastaParser(self.multiline))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('xyz', 'UUUUCCAAAAAG'))
        self.assertEqual(a[1].Name, 'xyz')

    def test_single_constructor(self):
        """FastaParser should use constructors if supplied"""
        f = list(FastaParser(self.oneseq, Dna))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('abc', 'TCAG'))
        self.assertEqual(a[1].Name, 'abc')

        def upper_abc(x):
            return None, {'ABC': x.upper()}

        f = list(FastaParser(self.multiline, Dna, upper_abc))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, (None, 'TTTTCCAAAAAG'))
        self.assertEqual(a[1].Name, None)
        self.assertEqual(a[1].Info.ABC, 'XYZ')

    def test_multiple(self):
        """FastaParser should read multiline records correctly"""
        f = list(FastaParser(self.threeseq))
        self.assertEqual(len(f), 3)
        for i in f:
            assert isinstance(i[1], Sequence)
        a, b, c = f
        self.assertEqual((a[1].Name, a[1]), ('123', 'a'))
        self.assertEqual((b[1].Name, b[1]), ('abc', 'caggac'))
        self.assertEqual((c[1].Name, c[1]), ('456', 'cg'))

    def test_multiple_bad(self):
        """Parser should complain or skip bad records"""
        self.assertRaises(RecordError, list, FastaParser(self.twogood))
        f = list(FastaParser(self.twogood, strict=False))
        self.assertEqual(len(f), 2)
        a, b = f
        a, b = a[1], b[1]   #field 0 is name
        self.assertEqual((a.Name, a), ('abc', 'caggac'))
        self.assertEqual((b.Name, b), ('456', 'cg'))

    def test_multiple_constructor_bad(self):
        """Parser should complain or skip bad records w/ constructor"""

        def dnastrict(x, **kwargs):
            try:
                return Dna(x, check=True, **kwargs)
            except Exception, e:
                raise RecordError, "Could not convert sequence"
        
        self.assertRaises(RecordError, list, FastaParser(self.oneX, dnastrict))
        f = list(FastaParser(self.oneX, dnastrict, strict=False))
        self.assertEqual(len(f), 2)
        a, b = f
        a, b = a[1], b[1]
        self.assertEqual((a.Name, a), ('abc', 'caggac'.upper()))
        self.assertEqual((b.Name, b), ('456', 'cg'.upper()))

class NcbiFastaLabelParserTests(TestCase):
    """Tests of the label line parser for NCBI's FASTA identifiers."""
    def test_init(self):
        """Labels from genpept.fsa should work as expected"""
        i = NcbiFastaLabelParser(
            '>gi|37549575|ref|XP_352503.1| similar to EST gb|ATTS1136')[1]
        self.assertEqual(i.GI, ['37549575'])
        self.assertEqual(i.RefSeq, ['XP_352503.1'])
        self.assertEqual(i.Description, 'similar to EST gb|ATTS1136')

        i = NcbiFastaLabelParser(
            '>gi|32398734|emb|CAD98694.1| (BX538350) dbj|baa86974.1, possible')[1]
        self.assertEqual(i.GI, ['32398734'])
        self.assertEqual(i.RefSeq, [])
        self.assertEqual(i.EMBL, ['CAD98694.1'])
        self.assertEqual(i.Description, '(BX538350) dbj|baa86974.1, possible')

        i = NcbiFastaLabelParser(
            '>gi|10177064|dbj|BAB10506.1| (AB005238)   ')[1]
        self.assertEqual(i.GI, ['10177064'])
        self.assertEqual(i.DDBJ, ['BAB10506.1'])
        self.assertEqual(i.Description, '(AB005238)')

class NcbiFastaParserTests(TestCase):
    """Tests of the NcbiFastaParser."""
    def setUp(self):
        """Define a few standard files"""
        self.peptide = [
'>gi|10047090|ref|NP_055147.1| small muscle protein, X-linked [Homo sapiens]',
'MNMSKQPVSNVRAIQANINIPMGAFRPGAGQPPRRKECTPEVEEGVPPTSDEEKKPIPGAKKLPGPAVNL',
'SEIQNIKSELKYVPKAEQ',
'>gi|10047092|ref|NP_037391.1| neuronal protein [Homo sapiens]',
'MANRGPSYGLSREVQEKIEQKYDADLENKLVDWIILQCAEDIEHPPPGRAHFQKWLMDGTVLCKLINSLY',
'PPGQEPIPKISESKMAFKQMEQISQFLKAAETYGVRTTDIFQTVDLWEGKDMAAVQRTLMALGSVAVTKD'
]
        self.nasty = [
'  ',                               #0  ignore leading blank line
'>gi|abc|ref|def|',                 #1  no description -- ok
'UCAG',                             #2  single line of sequence
'#comment',                         #3  comment -- skip
'  \t   ',                          #4  ignore blank line between records
'>gi|xyz|gb|qwe|  \tdescr   \t\t',  #5  desciption has whitespace
'UUUU',                             #6  two lines of sequence
'CCCC',                             #7  
'>gi|bad|ref|nonsense',             #8  missing last pipe -- error
'ACU',                              #9  
'>gi|bad|description',              #10 not enough fields -- error       
'AAA',                              #11
'>gi|bad|ref|stuff|label',          #12
'XYZ',                              #13 bad sequence -- error
'>gi|bad|gb|ignore| description',   #14 label without sequence -- error
'>  gi  |  123  | dbj  | 456 | desc|with|pipes| ',#15 label w/ whitespace -- OK
'ucag',                             #16
'  \t  ',                           #17 ignore blank line inside record
'UCAG',                             #18
'tgac',                             #19 lowercase should be OK
'# comment',                        #20 comment -- skip
'NNNN',                             #21 degenerates should be OK
'   ',                              #22 ignore trailing blank line
]
        self.empty = []
        self.no_label = ['ucag']

    def test_empty(self):
        """NcbiFastaParser should accept empty input"""
        self.assertEqual(list(NcbiFastaParser(self.empty)), [])
        self.assertEqual(list(NcbiFastaParser(self.empty, Protein)), [])

    def test_normal(self):
        """NcbiFastaParser should accept normal record if loose or strict"""
        f = list(NcbiFastaParser(self.peptide, Protein))
        self.assertEqual(len(f), 2)
        a, b = f
        a, b = a[1], b[1]   #field 0 is the name
        self.assertEqual(a, 'MNMSKQPVSNVRAIQANINIPMGAFRPGAGQPPRRKECTPEVEEGVPPTSDEEKKPIPGAKKLPGPAVNLSEIQNIKSELKYVPKAEQ')
        self.assertEqual(a.Info.GI, ['10047090'])
        self.assertEqual(a.Info.RefSeq, ['NP_055147.1'])
        self.assertEqual(a.Info.DDBJ, [])
        self.assertEqual(a.Info.Description, 
            'small muscle protein, X-linked [Homo sapiens]')

        self.assertEqual(b, 'MANRGPSYGLSREVQEKIEQKYDADLENKLVDWIILQCAEDIEHPPPGRAHFQKWLMDGTVLCKLINSLYPPGQEPIPKISESKMAFKQMEQISQFLKAAETYGVRTTDIFQTVDLWEGKDMAAVQRTLMALGSVAVTKD')
        self.assertEqual(b.Info.GI, ['10047092'])
        self.assertEqual(b.Info.RefSeq, ['NP_037391.1'])
        self.assertEqual(b.Info.Description, 'neuronal protein [Homo sapiens]')

    def test_bad(self):
        """NcbiFastaParser should raise error on bad records if strict"""
        #if strict, starting anywhere in the first 15 lines should cause errors
        for i in range(15):
            self.assertRaises(RecordError,list,NcbiFastaParser(self.nasty[i:]))
        #...but the 16th is OK.
        r = list(NcbiFastaParser(self.nasty[15:]))[0]
        self.assertEqual(r, ('123', 'ucagUCAGtgacNNNN'))
        #test that we get what we expect if not strict
        r = list(NcbiFastaParser(self.nasty, Sequence, strict=False))
        self.assertEqual(len(r), 4)
        a, b, c, d = r
        self.assertEqual((a[1], a[1].Info.GI, a[1].Info.RefSeq, \
            a[1].Info.Description), 
            ('UCAG', ['abc'], ['def'], ''))
        self.assertEqual((b[1], b[1].Info.GI, b[1].Info.GenBank, \
                b[1].Info.Description),
            ('UUUUCCCC', ['xyz'], ['qwe'], 'descr'))
        self.assertEqual((c[1], c[1].Info.GI, c[1].Info.RefSeq, \
            c[1].Info.Description),
            ('XYZ', ['bad'], ['stuff'], 'label'))
        self.assertEqual((d[1], d[1].Info.GI, d[1].Info.DDBJ, \
            d[1].Info.Description),
            ('ucagUCAGtgacNNNN'.upper(), ['123'], ['456'], 'desc|with|pipes|'))
        #...and when we explicitly supply a constructor
        r = list(NcbiFastaParser(self.nasty, Dna, strict=False))
        self.assertEqual(len(r), 3)
        a, b, c = r
        a, b, c = a[1], b[1], c[1]
        self.assertEqual((a, a.Info.GI, a.Info.RefSeq, a.Info.Description), 
            ('TCAG', ['abc'], ['def'], ''))
        self.assertEqual((b, b.Info.GI, b.Info.GenBank, b.Info.Description),
            ('TTTTCCCC', ['xyz'], ['qwe'], 'descr'))
        self.assertEqual((c, c.Info.GI, c.Info.DDBJ, c.Info.Description),
            ('tcagTCAGtgacNNNN'.upper(), ['123'], ['456'], 'desc|with|pipes|'))
    

class LabelParsingTest(TestCase):
    """Test generic fasta label parsing"""
    def test_rich_label(self):
        """rich label correctly constructs label strings"""
        # labels should be equal based on the result of applying their
        # attributes to their string template
        k = RichLabel(Info(species="rat"), "%(species)s")
        l = RichLabel(Info(species="rat", seq_id="xy5"), "%(species)s")
        self.assertEqual(k, l)
        
        # labels should construct from Info components correctly
        k = RichLabel(Info(species="rat", seq_id="xy5"),
                      "%(seq_id)s:%(species)s")
        self.assertEqual(k, "xy5:rat")
        k = RichLabel(Info(species="rat", seq_id="xy5"),
                      "%(species)s:%(seq_id)s")
        self.assertEqual(k, "rat:xy5")
        
        # extra components should be ignored
        k = RichLabel(Info(species="rat", seq_id="xy5"), "%(species)s")
        self.assertEqual(k, "rat")
        
        # the label should have Info object
        self.assertEqual(k.Info.species, "rat")
        self.assertEqual(k.Info.seq_id, "xy5")
        
        # label should be constructable just like a normal string
        self.assertEqual(RichLabel('a'), 'a')
    
    def test_label_parser(self):
        """label parser factory function cope with mixed structure labels"""
        # the label parser factory function should correctly handle label lines
        # with mixed separators
        make = LabelParser("%(species)s:%(accession)s",
                                [[0,"accession", str],
                                [2, "species", str]],
                                split_with=": ")
        for label, expect in [(">abcd:human:misc", "misc:abcd"),
                              ("abcd:human:misc", "misc:abcd"),
                              (">abcd:Human misc", "misc:abcd"),
                              (">abcd Human:misc", "misc:abcd"),
                              (">abcd:Human misc", "misc:abcd")]:
            self.assertEqual(make(label), expect)
        
        # should raise an assertion error if template doesn't match at least one field name
        self.assertRaises(AssertionError, LabelParser, "%s:%s",
                                    [[0,"accession", str],
                                    [2, "species", str]],
                                    split_with=": ")
    

class GroupFastaParsingTest(TestCase):
    """test parsing of grouped sequences in a collection"""
    def test_groups(self):
        """correctly yield grouped sequences from fasta formatted data"""
        data = [">group1:seq1_id:species1",
                "ACTG",
                ">group1:seq2_id:species2",
                "ACTG",
                ">group2:seq3_id:species1",
                "ACGT",
                ">group2:seq4_id:species2",
                "ACGT"]
        expected = [{"species1": "ACTG", "species2":"ACTG"},
                    {"species1":"ACGT", "species2":"ACGT"}]
        label_to_name = LabelParser("%(species)s", [(0,"Group",str),
                            (1,"seq_id",str),(2,"species",str)], split_with=":")
        parser = GroupFastaParser(data, label_to_name, aligned=True)
        count = 0
        for group in parser:
            got = group.todict()
            want = expected[count]
            self.assertEqual(got, want)
            self.assertEqual(group.Info.Group, "group%s" % (count+1))
            count += 1
        
        # check we don't return a done group
        done_groups = ["group1"]
        parser = GroupFastaParser(data, label_to_name, done_groups=done_groups,
                    aligned=True)
        for group in parser:
            got = group.todict()
            want = expected[1]
            self.assertEqual(got, want)
            self.assertEqual(group.Info.Group, "group2")
        
    

if __name__ == '__main__':
    main()
