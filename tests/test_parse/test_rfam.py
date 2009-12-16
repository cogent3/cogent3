#!/usr/bin/env python
"""
Provides tests for RfamParser and related classes and functions.
"""

from cogent.parse.rfam import is_header_line, is_seq_line, is_structure_line,\
HeaderToInfo, MinimalRfamParser, RfamFinder, NameToInfo, RfamParser,\
ChangedSequence, is_empty_or_html
from cogent.util.unit_test import TestCase, main
from cogent.parse.record import RecordError
from cogent.core.info import Info
from cogent.struct.rna2d import WussStructure
from cogent.core.alignment import Alignment
from cogent.core.moltype import BYTES

__author__ = "Sandra Smit and Greg Caporaso"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Sandra Smit", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"

Sequence = BYTES.Sequence

class RfamParserTests(TestCase):
    """ Tests componenets of the rfam parser, in the rfam.py file """

    def setUp(self):
        """ Construct some fake data for testing purposes """

        self._fake_headers = []
        temp = list(fake_headers.split('\n'))
        for line in temp:
            self._fake_headers.append(line.strip())
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

    def test_is_header_line(self):
        """is_header_line: functions correctly w/ various lines """
        self.assertEqual(is_header_line('#=GF'), True)
        self.assertEqual(is_header_line('#=GF AC   RF00001'), True)
        self.assertEqual(is_header_line('#=GF CC   until it is\
            required for transcription. '), True)

        self.assertEqual(is_header_line(''), False)
        self.assertEqual(is_header_line('X07545.1/505-619 '), False)
        self.assertEqual(is_header_line('#=G'), False)
        self.assertEqual(is_header_line('=GF'), False)
        self.assertEqual(is_header_line('#=GC SS_cons'), False)

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
        self.assertEqual(is_structure_line('#=GC SS_cons'), True)
        self.assertEqual(is_structure_line('#=GC SS_cons '), True)

        self.assertEqual(is_structure_line(''), False)
        self.assertEqual(is_structure_line(' '), False)
        self.assertEqual(is_structure_line('#=GF AC   RF00001'), False)
        self.assertEqual(is_structure_line('X07545.1/505-619'), False)
        self.assertEqual(is_structure_line('=GC SS_cons'), False)
        self.assertEqual(is_structure_line('#=GC'), False)
        self.assertEqual(is_structure_line('#=GC RF'), False)

    def test_HeaderToInfo(self):
        """HeaderToInfo: correctly builds info object from header information"""
        info = HeaderToInfo(self._fake_headers)
        self.assertEqual(info['Identification'], '5S_rRNA')
        self.assertEqual(info['RT'], None)
        self.assertEqual(info['Comment'], 'This is a short comment')
        self.assertEqual(info['Author'], 'Griffiths-Jones SR')
        self.assertEqual(info['Sequences'], '606')
        self.assertEqual(info['DatabaseReference'],\
            ['URL; http://oberon.fvms.ugent.be:8080/rRNA/ssu/index.html;',\
            'URL; http://rdp.cme.msu.edu/html/;'])
        self.assertEqual(info['PK'],'not real')

        self.assertEqual(info['Rfam'], ['RF00001'])

    def test_HeaderToInfo_invalid_data(self):
        """HeaderToInfo: correctly raises error when necessary """
        invalid_headers = [['#=GF ACRF00001'],['#=GFACRF00001']]
        for h in invalid_headers:
            self.assertRaises(RecordError, HeaderToInfo, h)

    def test_MinimalRfamParser_strict_missing_fields(self):
        """MinimalRfamParser: toggle strict functions w/ missing fields"""
        # strict = True
        
        self.assertRaises(RecordError,list,\
            MinimalRfamParser(self._fake_record_no_sequences))
        
        self.assertRaises(RecordError,list,\
            MinimalRfamParser(self._fake_record_no_structure))

        # strict = False
        # no header shouldn't be a problem
        self.assertEqual(list(MinimalRfamParser(self._fake_record_no_headers,\
            strict=False)), [([],{'Z11765.1/1-89':'GGUC'},'............>>>')])
        # should get empty on missing sequence or missing structure
        self.assertEqual(list(MinimalRfamParser(self._fake_record_no_sequences,\
            strict=False)), [])
        self.assertEqual(list(MinimalRfamParser(self._fake_record_no_structure,\
            strict=False)), [])

    def test_MinimalRfamParser_strict_invalid_sequence(self):
        """MinimalRfamParser: toggle strict functions w/ invalid seq
        """
        #strict = True
        self.assertRaises(RecordError,list,\
            MinimalRfamParser(self._fake_record_bad_sequence_1))

        # strict = False
        # you expect to get back as much information as possible, also
        # half records or sequences
        result = MinimalRfamParser(self._fake_record_bad_sequence_1,strict=False)
        self.assertEqual(len(list(MinimalRfamParser(\
            self._fake_record_bad_sequence_1,strict=False))[0][1].NamedSeqs), 3)            

    def test_MinimalRfamParser_strict_invalid_structure(self):
        """MinimalRfamParser: toggle strict functions w/ invalid structure
        """
        #strict = True
        self.assertRaises(RecordError,list,\
            MinimalRfamParser(self._fake_record_bad_structure_1))

        # strict = False
        self.assertEqual(list(MinimalRfamParser(\
            self._fake_record_bad_structure_1,strict=False))[0][2],None)                                

    def test_MinimalRfamParser_w_valid_data(self):
        """MinimalRfamParser: integrity of output """

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
        for r in MinimalRfamParser(self._fake_two_records, strict=False):
            data.append(r)
        self.assertEqual(data[0],(headers,sequences,structure))
        assert isinstance(data[0][1],Alignment)

        # This line tests that invalid entries are ignored when strict=False
        # Note, there are two records in self._fake_two_records, but 2nd is
        # invalid
        self.assertEqual(len(data),1)            
            
    def test_RfamFinder(self):
        """RfamFinder: integrity of output """
        fake_record = ['a','//','b','b','//']
        num_records = 0
        data = []
        for r in RfamFinder(fake_record):
            data.append(r)
            num_records += 1
        self.assertEqual(num_records, 2)
        self.assertEqual(data[0], ['a','//'])
        self.assertEqual(data[1], ['b','b','//'])

    def test_ChangedSequence(self):
        """ChangedSequence: integrity of output"""
        # Made up input, based on a line that would look like:
        # U17136.1/898-984  AACA..CAU..CAGAUUUCCU..GGUGUAA.CGAA
        
        s_in = 'AACA..CAU..CAGAUUUCCU..GGUGUAA.CGAA'
        s_out = 'AACA--CAU--CAGAUUUCCU--GGUGUAA-CGAA'
        sequence = ChangedSequence(s_in)
        
        self.assertEqual(sequence, s_out)

        # test some extremes on the seq
        # sequence of all blanks
        s_in = '.' * 5
        s_out = '-' * 5
        sequence = ChangedSequence(s_in)

        self.assertEqual(sequence, s_out)

        # sequence of no blanks
        s_in = 'U' * 5
        s_out = 'U' * 5
        sequence = ChangedSequence(s_in)

        self.assertEqual(sequence, s_out)


    def test_NameToInfo(self):
        """NameToInfo: integrity of output """
        # Made up input, based on a line that would look like:
        # U17136.1/898-984  AACA..CAU..CAGAUUUCCU..GGUGUAA.CGAA
        
        s_in = 'AACA..CAU..CAGAUUUCCU..GGUGUAA.CGAA'
        #s_out = 'AACA--CAU--CAGAUUUCCU--GGUGUAA-CGAA'
        sequence = Sequence(s_in, Name='U17136.1/898-984')
        info = NameToInfo(sequence)
        
        #self.assertEqual(seq, s_out)
        self.assertEqual(info['Start'], 897)
        self.assertEqual(info['End'], 984)
        self.assertEqual(info['GenBank'], ['U17136.1'])


    def test_NameToInfo_invalid_label(self):
        """NameToInfo: raises error on invalid label """
        s = 'AA'
        invalid_labels = ['U17136.1898-984','U17136.1/898984']
        for l in invalid_labels:
            self.assertRaises(RecordError,NameToInfo,\
                Sequence(s, Name=l))
        a = 'U17136.1/' #missing start/end positions
        b = '/898-984' #missing genbank id
        obs_info = NameToInfo(Sequence(s,Name=a))
        exp = Info({'GenBank':'U17136.1','Start':None,'End':None})
        self.assertEqual(obs_info,exp)
        obs_info = NameToInfo(Sequence(s,Name=b))
        exp = Info({'GenBank':None,'Start':897,'End':984})
        self.assertEqual(obs_info,exp)

        #strict = False
        # in strict mode you want to get back as much info as possible
        lab1 = 'U17136.1898-984'
        lab2 = 'U17136.1/898984'
        obs_info = NameToInfo(Sequence(s,Name=lab1), strict=False)
        exp = Info({'GenBank':None,'Start':None,'End':None})
        self.assertEqual(obs_info,exp)
        obs_info = NameToInfo(Sequence(s,Name=lab2), strict=False)
        exp = Info({'GenBank':'U17136.1','Start':None,'End':None})
        self.assertEqual(obs_info,exp)
                    
        
    def test_RfamParser(self):
        """RfamParser: integrity of output """

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
             
        for r in RfamParser(self._fake_record):
            headers,sequences,structure = r
            
            self.assertEqual(headers['Refs']['Rfam'], ['RF00014'])
            self.assertEqual(headers['Author'], 'Mifsud W')
            self.assertEqualItems(sequences.values(), expected_sequences)
            assert isinstance(sequences, Alignment)
            self.assertEqualItems([s.Info.GenBank for s in sequences.Seqs],
                [['U17136.1'],['M15749.1'],['AF090431.1']])
            self.assertEqualItems([s.Info.Start for s in sequences.Seqs],
                [897,154,221])
            self.assertEqual(structure, expected_structure)
            assert isinstance(structure,WussStructure)

    def test_RfamParser_strict_missing_fields(self):
        """RfamParser: toggle strict functions correctly """
        # strict = True
        self.assertRaises(RecordError,list,\
            RfamParser(self._fake_record_no_headers))
        
        self.assertRaises(RecordError,list,\
            RfamParser(self._fake_record_no_sequences))
        
        self.assertRaises(RecordError,list,\
            RfamParser(self._fake_record_no_structure))
                        
        # strict = False
        self.assertEqual(list(RfamParser(self._fake_record_no_headers,\
            strict=False)), [])
        self.assertEqual(list(RfamParser(self._fake_record_no_sequences,\
            strict=False)), [])
        self.assertEqual(list(RfamParser(self._fake_record_no_structure,\
            strict=False)), [])            

    def test_RFamParser_strict_invalid_headers(self):
        """RfamParser: functions when toggling strict w/ record w/ bad header
        """
        self.assertRaises(RecordError,list,\
            RfamParser(self._fake_record_bad_header_1))
            
        self.assertRaises(RecordError,list,\
            RfamParser(self._fake_record_bad_header_2))
        
        # strict = False
        x =  list(RfamParser(self._fake_record_bad_header_1, strict=False))
        obs = list(RfamParser(self._fake_record_bad_header_1,\
            strict=False))[0][0].keys()
        self.assertEqual(len(obs),1)
        obs = list(RfamParser(self._fake_record_bad_header_2,\
            strict=False))[0][0].keys()
        self.assertEqual(len(obs),1)

    def test_RfamParser_strict_invalid_sequences(self):
        """RfamParser: functions when toggling strict w/ record w/ bad seq
        """
        self.assertRaises(RecordError,list,
            MinimalRfamParser(self._fake_record_bad_sequence_1))
            
        # strict = False
        # in 'False' mode you expect to get back as much as possible, also
        # parts of sequences
        self.assertEqual(len(list(RfamParser(self._fake_record_bad_sequence_1,\
            strict=False))[0][1].NamedSeqs), 3)           
                            
    def test_RfamParser_strict_invalid_structure(self):
        """RfamParser: functions when toggling strict w/ record w/ bad struct
        """
        # strict 
        self.assertRaises(RecordError,list,\
            RfamParser(self._fake_record_bad_structure_2))
        #not strict
        self.assertEqual(list(RfamParser(self._fake_record_bad_structure_2,\
        strict=False)),[])

    def test_RfamParser_single_family(self):
        """RfamParser: should work on a single family in stockholm format"""
        exp_header = Info()
        exp_aln = {'K02120.1/628-682':\
            'AUGGGAAAUUCCCCCUCCUAUAACCCCCCCGCUGGUAUCUCCCCCUCAGACUGGC',\
            'D00647.1/629-683':\
            'AUGGGAAACUCCCCCUCCUAUAACCCCCCCGCUGGCAUCUCCCCCUCAGACUGGC'}
        exp_struct = '<<<<<<.........>>>>>>.........<<<<<<.............>>>>>>'
        h, a, s = list(RfamParser(self.single_family))[0]
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
#=GF AUMifsud W

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
