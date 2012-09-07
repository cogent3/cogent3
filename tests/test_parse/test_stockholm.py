#!/usr/bin/env python
"""
Provides tests for StockholmParser and related classes and functions.
"""

from cogent.parse.stockholm import is_gf_line, is_gc_line, is_gs_line, \
    is_gr_line, is_seq_line, is_structure_line, GfToInfo, GcToInfo, GsToInfo, \
    GrToInfo, MinimalStockholmParser, StockholmFinder, \
    StockholmParser, Sequence, is_empty_or_html
from cogent.util.unit_test import TestCase, main
from cogent.parse.record import RecordError
from cogent.core.info import Info
from cogent.struct.rna2d import WussStructure
from cogent.core.alignment import Alignment
from cogent.core.moltype import BYTES

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.2-dev"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

Sequence = BYTES.Sequence

class StockholmParserTests(TestCase):
    """ Tests componenets of the stockholm parser, in the stockholm.py file """

    def setUp(self):
        """ Construct some fake data for testing purposes """

        self._fake_headers = []
        temp = list(fake_headers.split('\n'))
        for line in temp:
            self._fake_headers.append(line.strip())
        del temp
        
        self._fake_gc_annotation = []
        temp = list(fake_gc_annotation.split('\n'))
        for line in temp:
            self._fake_gc_annotation.append(line.strip())
        del temp

        self._fake_gs_annotation = []
        temp = list(fake_gs_annotation.split('\n'))
        for line in temp:
            self._fake_gs_annotation.append(line.strip())
        del temp

        self._fake_gr_annotation = []
        temp = list(fake_gr_annotation.split('\n'))
        for line in temp:
            self._fake_gr_annotation.append(line.strip())
        del temp
        
        self._fake_record_no_headers =\
            list(fake_record_no_headers.split('\n'))

        self._fake_record_no_sequences =\
            list(fake_record_no_sequences.split('\n'))

        self._fake_record_no_structure =\
            list(fake_record_no_structure.split('\n'))

        self._fake_two_records =\
            list(fake_two_records.split('\n'))
            
        self._fake_record =\
            list(fake_record.split('\n'))

        self._fake_record_bad_header_1 =\
            list(fake_record_bad_header_1.split('\n'))
            
        self._fake_record_bad_header_2 =\
            list(fake_record_bad_header_2.split('\n'))

        self._fake_record_bad_sequence_1 =\
            list(fake_record_bad_sequence_1.split('\n'))

        self._fake_record_bad_structure_1 =\
            list(fake_record_bad_structure_1.split('\n'))                                                    
        self._fake_record_bad_structure_2 =\
            list(fake_record_bad_structure_2.split('\n'))

        self.single_family = single_family.split('\n')
            
    def test_is_empty_or_html(self):
        """is_empty_or_html: should ignore empty and HTML line"""
        line = '        '
        self.assertEqual(is_empty_or_html(line), True)
        line = '\n\n'
        self.assertEqual(is_empty_or_html(line), True)
        line = '<pre>'
        self.assertEqual(is_empty_or_html(line), True)
        line = '</pre>\n\n'
        self.assertEqual(is_empty_or_html(line), True)
        line = '\t<//\n'
        self.assertEqual(is_empty_or_html(line), False)

    def test_is_gf_line(self):
        """is_gf_line: functions correctly w/ various lines """
        self.assertEqual(is_gf_line('#=GF'), True)
        self.assertEqual(is_gf_line('#=GF AC   RF00001'), True)
        self.assertEqual(is_gf_line('#=GF CC   until it is\
            required for transcription. '), True)

        self.assertEqual(is_gf_line(''), False)
        self.assertEqual(is_gf_line('X07545.1/505-619 '), False)
        self.assertEqual(is_gf_line('#=G'), False)
        self.assertEqual(is_gf_line('=GF'), False)
        self.assertEqual(is_gf_line('#=GC SS_cons'), False)
        
    def test_is_gc_line(self):
        """is_gc_line: functions correctly w/ various lines """
        self.assertEqual(is_gc_line('#=GC'), True)
        self.assertEqual(is_gc_line('#=GC SS_cons'), True)
        self.assertEqual(is_gc_line('#=GC RF'), True)

        self.assertEqual(is_gc_line(''), False)
        self.assertEqual(is_gc_line('X07545.1/505-619 '), False)
        self.assertEqual(is_gc_line('#=G'), False)
        self.assertEqual(is_gc_line('=GF'), False)
        self.assertEqual(is_gc_line('#=GR SS'), False)
        
    def test_is_gs_line(self):
        """is_gs_line: functions correctly w/ various lines """
        self.assertEqual(is_gs_line('#=GS'), True)
        self.assertEqual(is_gs_line('#=GS Seq1   AC'), True)
        self.assertEqual(is_gs_line('#=GS Seq1   DE'), True)

        self.assertEqual(is_gs_line(''), False)
        self.assertEqual(is_gs_line('X07545.1/505-619 '), False)
        self.assertEqual(is_gs_line('#=G'), False)
        self.assertEqual(is_gs_line('=GF'), False)
        self.assertEqual(is_gs_line('#=GC SS_cons'), False)
        
    def test_is_gr_line(self):
        """is_gr_line: functions correctly w/ various lines """
        self.assertEqual(is_gr_line('#=GR'), True)
        self.assertEqual(is_gr_line('#=GR SS   ..<<..>>..'), True)
        self.assertEqual(is_gr_line('#=GR RF   cGGacG'), True)

        self.assertEqual(is_gr_line(''), False)
        self.assertEqual(is_gr_line('X07545.1/505-619 '), False)
        self.assertEqual(is_gr_line('#=G'), False)
        self.assertEqual(is_gr_line('=GF'), False)
        self.assertEqual(is_gr_line('#=GC SS_cons'), False)

    def test_is_seq_line(self):
        """is_seq_line: functions correctly w/ various lines """
        s = 'X07545.1/505-619                     .\
            .ACCCGGC.CAUA...GUGGCCG.GGCAA.CAC.CCGG.U.C..UCGUU'
        assert is_seq_line('s')
        assert is_seq_line('X07545.1/505-619')
        assert is_seq_line('M21086.1/8-123')

        assert not is_seq_line('')
        assert not is_seq_line('#GF=')
        assert not is_seq_line('//blah')

    def test_is_structure_line(self):
        """is_structure_line: functions correctly w/ various lines """
        s = '#=GC SS_cons\
            <<<<<<<<<........<<.<<<<.<...<.<...<<<<.<.<.......'
        self.assertEqual(is_structure_line(s), True)
        self.assertEqual(is_structure_line('#=GC SS_cons'), False)
        self.assertEqual(is_structure_line('#=GC SS_cons2'), False)
        self.assertEqual(is_structure_line('#=GC SS_cons '), True)

        self.assertEqual(is_structure_line(''), False)
        self.assertEqual(is_structure_line(' '), False)
        self.assertEqual(is_structure_line('#=GF AC   RF00001'), False)
        self.assertEqual(is_structure_line('X07545.1/505-619'), False)
        self.assertEqual(is_structure_line('=GC SS_cons'), False)
        self.assertEqual(is_structure_line('#=GC'), False)
        self.assertEqual(is_structure_line('#=GC RF'), False)

    def test_GfToInfo(self):
        """GfToInfo: correctly builds info object from header information"""
        info = GfToInfo(self._fake_headers)
        self.assertEqual(info['AccessionNumber'], 'RF00001')
        self.assertEqual(info['Identification'], '5S_rRNA')
        self.assertEqual(info['Comment'], 'This is a short comment')
        self.assertEqual(info['Author'], 'Griffiths-Jones SR')
        self.assertEqual(info['Sequences'], '606')
        self.assertEqual(info['DatabaseReference'],\
            ['URL; http://oberon.fvms.ugent.be:8080/rRNA/ssu/index.html;',\
            'URL; http://rdp.cme.msu.edu/html/;'])
        self.assertEqual(info['PK'],'not real')

    def test_GfToInfo_invalid_data(self):
        """GfToInfo: correctly raises error when necessary """
        invalid_headers = [['#=GF ACRF00001'],['#=GFACRF00001']]
        for h in invalid_headers:
            self.assertRaises(RecordError, GfToInfo, h)
    
    def test_GcToInfo(self):
        """GcToInfo: correctly builds info object from header information"""
        info = GcToInfo(self._fake_gc_annotation)
        self.assertEqual(info['ConsensusSecondaryStructure'], \
            '..........<<<<<<<<<<.....>>>>>>>>>>..')
        self.assertEqual(info['ReferenceAnnotation'], \
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')

    def test_GcToInfo_invalid_data(self):
        """GcToInfo: correctly raises error when necessary """
        invalid_headers = [['#=GCSS_cons ..<<..>>..'],['#=GCSAxxxxxxx']]
        for h in invalid_headers:
            self.assertRaises(RecordError, GcToInfo, h)

    def test_GsToInfo(self):
        """GsToInfo: correctly builds info object from header information"""
        info = GsToInfo(self._fake_gs_annotation)
        self.assertEqual(info['BasePair'], \
            {'1N77_C':['0 70 cWW CC','1 69 cWW CC','2 68 cWW CC',\
                '3 67 cWW CC']})

    def test_GsToInfo_invalid_data(self):
        """GsToInfo: correctly raises error when necessary """
        invalid_headers = [['#=GSBPS 0 10 cwW CC'],['#=GSACRF00001']]
        for h in invalid_headers:
            self.assertRaises(RecordError, GsToInfo, h)

    def test_GrToInfo(self):
        """GrToInfo: correctly builds info object from header information"""
        info = GrToInfo(self._fake_gr_annotation)
        self.assertEqual(info['SecondaryStructure'], \
            {'1N77_C':'..........<<<<<<<<<<.....>>>>>>>>>>..'})
        self.assertEqual(info['ReferenceAnnotation'], \
            {'1N77_C':'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'})

    def test_GrToInfo_invalid_data(self):
        """GrToInfo: correctly raises error when necessary """
        invalid_headers = [['#=GRSS ..<<..>>..'],['#=GRSAxxxxxx']]
        for h in invalid_headers:
            self.assertRaises(RecordError, GrToInfo, h)

    def test_StockholmStockholmParser_strict_missing_fields(self):
        """MinimalStockholmParser: toggle strict functions w/ missing fields"""
        # strict = True
        
        self.assertRaises(RecordError,list,\
            MinimalStockholmParser(self._fake_record_no_sequences))

        # strict = False
        # no header shouldn't be a problem
        headers, aln, struct = \
            list(MinimalStockholmParser(self._fake_record_no_headers,\
                strict=False))[0]
        self.assertEqual((headers,aln.todict(),str(struct)), \
            ({'GS':[],'GF':[],'GR':[],\
                'GC':['#=GC SS_cons                         ............>>>']},\
                {'Z11765.1/1-89':'GGUC'},'............>>>'))
        # should get empty on missing sequence or missing structure
        self.assertEqual(list(MinimalStockholmParser(\
            self._fake_record_no_sequences,\
            strict=False)), [])

    def test_MinimalStockholmParser_strict_invalid_sequence(self):
        """MinimalStockholmParser: toggle strict functions w/ invalid seq
        """
        #strict = True
        self.assertRaises(RecordError,list,\
            MinimalStockholmParser(self._fake_record_bad_sequence_1))

        # strict = False
        # you expect to get back as much information as possible, also
        # half records or sequences
        result = MinimalStockholmParser(\
            self._fake_record_bad_sequence_1,strict=False)
        self.assertEqual(len(list(MinimalStockholmParser(\
            self._fake_record_bad_sequence_1,strict=False))[0][1].NamedSeqs), 3)            

    def test_StockholmParser_strict_invalid_structure(self):
        """StockholmParser: toggle strict functions w/ invalid structure
        """
        #strict = True
        self.assertRaises(RecordError,list,\
            StockholmParser(self._fake_record_bad_structure_1))

        # strict = False
        self.assertEqual(list(MinimalStockholmParser(\
            self._fake_record_bad_structure_1,strict=False))[0][2],None)                                

    def test_MinimalStockholmParser_w_valid_data(self):
        """MinimalStockholmParser: integrity of output """

        # Some ugly constructions here, but this is what the output of
        # parsing fake_two_records should be
        headers = ['#=GF AC   RF00014','#=GF AU   Mifsud W']
        sequences =\
        {'U17136.1/898-984':\
        ''.join(['AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA',\
            'AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU']),\
        'M15749.1/155-239':\
        ''.join(['AACGCAUCGGAUUUCCCGGUGUAACGAA-UUUUCAAGUGCUUCUUGCAUU',\
            'AGCAAGUUUGAUCCCGACUCCUG-CGAGUCGGGAUUU']),\
        'AF090431.1/222-139':\
        ''.join(['CUCACAUCAGAUUUCCUGGUGUAACGAA-UUUUCAAGUGCUUCUUGCAUA',\
            'AGCAAGUUUGAUCCCGACCCGU--AGGGCCGGGAUUU'])}

        structure = WussStructure(''.join(\
        ['...<<<<<<<.....>>>>>>>....................<<<<<...',\
        '.>>>>>....<<<<<<<<<<.....>>>>>>>>>>..']))
        
        data = []
        for r in MinimalStockholmParser(self._fake_two_records, strict=False):
            data.append(r)
        self.assertEqual(\
            (data[0][0]['GF'],data[0][1].todict(),\
                str(data[0][2])),(headers,sequences,structure))
        assert isinstance(data[0][1],Alignment)

        # This line tests that invalid entries are ignored when strict=False
        # Note, there are two records in self._fake_two_records, but 2nd is
        # invalid
        self.assertEqual(len(data),1)            
            
    def test_StockholmFinder(self):
        """StockholmFinder: integrity of output """
        fake_record = ['a','//','b','b','//']
        num_records = 0
        data = []
        for r in StockholmFinder(fake_record):
            data.append(r)
            num_records += 1
        self.assertEqual(num_records, 2)
        self.assertEqual(data[0], ['a','//'])
        self.assertEqual(data[1], ['b','b','//'])
                    
        
    def test_StockholmParser(self):
        """StockholmParser: integrity of output """

        expected_sequences =\
        [''.join(['AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA',\
            'AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU']),\
        ''.join(['AACGCAUCGGAUUUCCCGGUGUAACGAA-UUUUCAAGUGCUUCUUGCAUU',\
            'AGCAAGUUUGAUCCCGACUCCUG-CGAGUCGGGAUUU']),\
        ''.join(['CUCACAUCAGAUUUCCUGGUGUAACGAA-UUUUCAAGUGCUUCUUGCAUA',\
            'AGCAAGUUUGAUCCCGACCCGU--AGGGCCGGGAUUU'])]
        expected_structure = ''.join(\
        ['...<<<<<<<.....>>>>>>>....................<<<<<...',\
        '.>>>>>....<<<<<<<<<<.....>>>>>>>>>>..'])            
             
        for r in StockholmParser(self._fake_record):
            headers = r.Info
            sequences = r
            structure = r.Info['Struct']
            self.assertEqual(headers['GF']['AccessionNumber'], 'RF00014')
            self.assertEqual(headers['GF']['Author'], 'Mifsud W')
            self.assertEqualItems(sequences.values(), expected_sequences)
            assert isinstance(sequences, Alignment)
            self.assertEqual(structure, expected_structure)
            assert isinstance(structure,WussStructure)

    def test_StockholmParser_strict_missing_fields(self):
        """StockholmParser: toggle strict functions correctly """
        # strict = True
        self.assertRaises(RecordError,list,\
            StockholmParser(self._fake_record_no_headers))
                        
        # strict = False
        self.assertEqual(list(StockholmParser(self._fake_record_no_headers,\
            strict=False)), [])
        self.assertEqual(list(StockholmParser(self._fake_record_no_sequences,\
            strict=False)), [])

    def test_StockholmParser_strict_invalid_headers(self):
        """StockholmParser: functions when toggling strict record w/ bad header
        """
        self.assertRaises(RecordError,list,\
            StockholmParser(self._fake_record_bad_header_1))
            
        self.assertRaises(RecordError,list,\
            StockholmParser(self._fake_record_bad_header_2))
        
        # strict = False
        x =  list(StockholmParser(self._fake_record_bad_header_1, strict=False))
        obs = list(StockholmParser(self._fake_record_bad_header_1,\
            strict=False))[0].Info.GF.keys()
        self.assertEqual(len(obs),1)
        obs = list(StockholmParser(self._fake_record_bad_header_2,\
            strict=False))[0].Info.GF.keys()
        self.assertEqual(len(obs),1)

    def test_StockholmParser_strict_invalid_sequences(self):
        """StockholmParser: functions when toggling strict w/ record w/ bad seq
        """
        self.assertRaises(RecordError,list,
            MinimalStockholmParser(self._fake_record_bad_sequence_1))
            
        # strict = False
        # in 'False' mode you expect to get back as much as possible, also
        # parts of sequences
        self.assertEqual(len(list(StockholmParser(\
            self._fake_record_bad_sequence_1,\
            strict=False))[0].NamedSeqs), 3)           
                            
    def test_StockholmParser_strict_invalid_structure(self):
        """StockholmParser: functions when toggling strict record w/ bad struct
        """
        # strict 
        self.assertRaises(RecordError,list,\
            StockholmParser(self._fake_record_bad_structure_2))
        #not strict
        self.assertEqual(list(StockholmParser(\
            self._fake_record_bad_structure_2,\
        strict=False)),[])

    def test_StockholmParser_single_family(self):
        """StockholmParser: should work on a family in stockholm format"""
        exp_header = {}
        exp_aln = {'K02120.1/628-682':\
            'AUGGGAAAUUCCCCCUCCUAUAACCCCCCCGCUGGUAUCUCCCCCUCAGACUGGC',\
            'D00647.1/629-683':\
            'AUGGGAAACUCCCCCUCCUAUAACCCCCCCGCUGGCAUCUCCCCCUCAGACUGGC'}
        exp_struct = '<<<<<<.........>>>>>>.........<<<<<<.............>>>>>>'
        aln = list(StockholmParser(self.single_family))[0]
        h = aln.Info['GF']
        a = aln
        s = aln.Info['Struct']
        self.assertEqual(h,exp_header)
        self.assertEqual(a,exp_aln)
        self.assertEqual(s,exp_struct)
        
        

# This is an altered version of some header info from Rfam.seed modified to
# incorporate different cases for testing
fake_headers = """#=GF AC   RF00001
#=GF AU   Griffiths-Jones SR
#=GF ID   5S_rRNA
#=GF RT   5S Ribosomal RNA Database.
#=GF DR   URL; http://oberon.fvms.ugent.be:8080/rRNA/ssu/index.html;
#=GF DR   URL; http://rdp.cme.msu.edu/html/;
#=GF CC   This is a short
#=GF CC   comment
#=GF SQ   606
#=GF PK   not real"""

fake_gc_annotation = """#=GC SS_cons    ..........<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
"""

fake_gs_annotation = """#=GS 1N77_C BP 0 70 cWW CC#=GS 1N77_C BP 1 69 cWW CC#=GS 1N77_C BP 2 68 cWW CC#=GS 1N77_C BP 3 67 cWW CC
"""

fake_gr_annotation = """#=GR 1N77_C SS    ..........<<<<<<<<<<.....>>>>>>>>>>..
#=GR 1N77_C RF    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
"""

fake_record_no_headers ="""Z11765.1/1-89                        GGUC
#=GC SS_cons                         ............>>>
//"""

fake_record_no_sequences ="""#=GF AC   RF00006
#=GC SS_cons                         ............>
//"""

fake_record_no_structure ="""#=GF AC   RF00006

Z11765.1/1-89                        GGUCAGC
//"""

fake_two_records ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GF AU   Mifsud W

U17136.1/898-984               AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons                   ...<<<<<<<.....>>>>>>>....................<<<<<...
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//
#=GF AC   RF00015
//"""

fake_record ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GF AU   Mifsud W

U17136.1/898-984               AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons                   ...<<<<<<<.....>>>>>>>....................<<<<<...
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//"""

fake_record_bad_header_1 ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GF AUMifsudW

U17136.1/898-984               AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons                   ...<<<<<<<.....>>>>>>>....................<<<<<...
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//"""

fake_record_bad_header_2 ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GFAUMifsud W

U17136.1/898-984               AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons                   ...<<<<<<<.....>>>>>>>....................<<<<<...
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//"""

fake_record_bad_sequence_1 ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GF AU   Mifsud W

U17136.1/898-984AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons                   ...<<<<<<<.....>>>>>>>....................<<<<<...
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//"""

fake_record_bad_structure_1 ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GF AU   Mifsud W

U17136.1/898-984               AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons...<<<<<<<.....>>>>>>>....................<<<<<...
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//"""

fake_record_bad_structure_2 ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GF AU   Mifsud W

U17136.1/898-984               AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons                   ...<<<<<<<.....>>>>>>>....................<<<<<!!!
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//"""

single_family=\
"""K02120.1/628-682      AUGGGAAAUUCCCCCUCCUAUAACCCCCCCGCUGGUAUCUCCCCCUCAGA
D00647.1/629-683      AUGGGAAACUCCCCCUCCUAUAACCCCCCCGCUGGCAUCUCCCCCUCAGA
#=GC SS_cons          <<<<<<.........>>>>>>.........<<<<<<.............>

K02120.1/628-682      CUGGC
D00647.1/629-683      CUGGC
#=GC SS_cons          >>>>>
//"""

# Run tests if called from the command line
if __name__ == '__main__':
    main()   
