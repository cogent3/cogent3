#!/usr/bin/env python

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Production"

"""
Test code for kegg_ko.py in cogent.parse.  
"""
from cogent.util.unit_test import TestCase, main

from cogent.parse.kegg_ko import kegg_label_fields,\
  parse_kegg_taxonomy, ko_record_iterator, ko_record_splitter,\
  ko_default_parser, ko_first_field_parser, delete_comments,\
  ko_colon_fields, ko_colon_delimited_parser, _is_new_kegg_rec_group,\
  group_by_end_char, class_lines_to_fields, ko_class_parser, parse_ko,\
  parse_ko_file, make_tab_delimited_line_parser



class ParseKOTests(TestCase):
    def make_tab_delimited_line_parser(self):
        """make_tab_delimited_line_parser should generate line parser"""
        line ="good\tbad:good\tgood\tgood\tbad:good\tgood"
        parse_fn = make_tab_delimited_line_parser([0,2,3,5])
        obs = parse_fn(line)
        exp = "good\tgood\tgood\tgood\tgood\tgood"
        self.assertEqual(obs,exp)

    def test_kegg_label_fields(self):
        """kegg_label_fields should return fields from line"""
        # Format is species:gene_id [optional gene_name]; description.
        # Note that the '>' should already be stripped by the Fasta Parser
        test1 = \
          """stm:STM0001  thrL; thr operon leader peptide ; K08278 thr operon leader peptide"""
        test2 = \
          """stm:STM0002  thrA; bifunctional aspartokinase I/homeserine dehydrogenase I (EC:2.7.2.4 1.1.1.13); K00003 homoserine dehydrogenase [EC:1.1.1.3]; K00928 aspartate kinase [EC:2.7.2.4]"""
        obs = kegg_label_fields(test1)
        exp = ('stm:STM0001','stm','STM0001',\
              'thrL','thr operon leader peptide ; K08278 thr operon leader peptide')
        self.assertEqual(obs,exp)

        obs = kegg_label_fields(test2)
        exp = ('stm:STM0002', 'stm', 'STM0002', 'thrA', \
            'bifunctional aspartokinase I/homeserine dehydrogenase I (EC:2.7.2.4 1.1.1.13); K00003 homoserine dehydrogenase [EC:1.1.1.3]; K00928 aspartate kinase [EC:2.7.2.4]')
        
        self.assertEqual(obs,exp)

    def test_ko_record_iterator(self):
        """ko_record_iterator should iterate over KO records"""
        recs = []
        for rec in ko_record_iterator(TEST_KO_LINES):
            recs.append(rec)
        
        self.assertEqual(len(recs),3)
        self.assertEqual(len(recs[0]),31)
        
        exp = 'ENTRY       K01559                      KO\n'
        self.assertEqual(recs[0][0],exp) 
          
        exp = '            RCI: RCIX1162 RCIX2396\n'
        self.assertEqual(recs[0][-1],exp) 
    
        exp = 'ENTRY       K01561                      KO\n'
        self.assertEqual(recs[-1][0],exp) 
    
        exp = '            MSE: Msed_1088\n'
        self.assertEqual(recs[-1][-1],exp) 
    
    def test_ko_record_splitter(self):
        """ko_record_splitter should split ko lines into a dict of groups"""
        
        recs=[rec for rec in ko_record_iterator(TEST_KO_LINES)]
        split_recs =  ko_record_splitter(recs[0])
        exp = ['GENES       AFM: AFUA_4G13070\n',\
               '            PHA: PSHAa2393\n',\
               '            ABO: ABO_0668\n',\
               '            BXE: Bxe_B0037 Bxe_C0683 Bxe_C1002 Bxe_C1023\n',\
               '            MPT: Mpe_A2274\n',\
               '            BBA: Bd0910(catD)\n',\
               '            GBE: GbCGDNIH1_0998 GbCGDNIH1_1171\n',\
               '            FNU: FN1345\n', \
               '            RBA: RB13257\n',\
               '            HMA: rrnAC1925(mhpC)\n',\
               '            RCI: RCIX1162 RCIX2396\n'] 
        self.assertEqual(exp,split_recs["GENES"])
        exp = ['CLASS       Metabolism; Biosynthesis of Secondary Metabolites; Limonene and\n', '            pinene degradation [PATH:ko00903]\n', '            Metabolism; Xenobiotics Biodegradation and Metabolism; Caprolactam\n', '            degradation [PATH:ko00930]\n', '            Metabolism; Xenobiotics Biodegradation and Metabolism;\n', '            1,1,1-Trichloro-2,2-bis(4-chlorophenyl)ethane (DDT) degradation\n', '            [PATH:ko00351]\n', '            Metabolism; Xenobiotics Biodegradation and Metabolism; Benzoate\n', '            degradation via CoA ligation [PATH:ko00632]\n', '            Metabolism; Xenobiotics Biodegradation and Metabolism; Benzoate\n', '            degradation via hydroxylation [PATH:ko00362]\n']
        
    def test_ko_default_parser(self):
        """ko_default parser should strip out newlines and join lines together"""
        
        # Applies to 'NAME' and 'DEFINITION' lines 
        
        default_line_1 = ['NAME        E3.8.1.2\n']
        obs = ko_default_parser(default_line_1)
        self.assertEqual(obs,'E3.8.1.2')
        
        default_line_2 = ['DEFINITION  2-haloacid dehalogenase [EC:3.8.1.2]\n']
        obs = ko_default_parser(default_line_2)
        self.assertEqual(obs,'2-haloacid dehalogenase [EC:3.8.1.2]')

    def test_ko_first_field_parser(self):
        """ko_first_field_parser should strip out newlines and join lines
        together (first field only)"""
        obs = ko_first_field_parser(\
         ['ENTRY       K01559                      KO\n'])
        exp = 'K01559'
        self.assertEqual(obs,exp)
    
    def test_delete_comments(self):
        """delete_comments should delete parenthetical comments from lines"""

        test_line = \
          "bifunctional aspartokinase I/homeserine dehydrogenase I (EC:2.7.2.4 1.1.1.13);"
        exp = "bifunctional aspartokinase I/homeserine dehydrogenase I ;"
        obs = delete_comments(test_line)
        self.assertEqual(obs,exp)

        nested_test_line = \
          "text(comment1(comment2));"
        exp = "text;"
        obs = delete_comments(nested_test_line)
        self.assertEqual(obs,exp)

    def test_ko_colon_fields(self):
        """ko_colon_fields should convert lines to (key, [list of values])"""
        test_lines =\
            ['            BXE: Bxe_B0037 Bxe_C0683 Bxe_C1002 Bxe_C1023\n']
        
        obs = ko_colon_fields(test_lines)
        exp = ('BXE', ['Bxe_B0037', 'Bxe_C0683', 'Bxe_C1002', 'Bxe_C1023'])
        self.assertEqual(obs,exp)

        test_lines = ['            HMA: rrnAC1925(mhpC)\n'] 
        obs = ko_colon_fields(test_lines, without_comments = True)
        exp = ('HMA', ['rrnAC1925'])
        self.assertEqual(obs,exp)
    
 
        test_lines = ['            HMA: rrnAC1925(mhpC)\n'] 
        obs = ko_colon_fields(test_lines, without_comments = False)
        exp = ('HMA', ['rrnAC1925(mhpC)'])
        self.assertEqual(obs,exp)
    
    def test_ko_colon_delimited_parser(self):
        """ko_colon_delimited_parser should return a dict of id: values for
        colon delimited lines"""
        test_lines =\
           ['GENES       AFM: AFUA_4G13070\n',\
            '            PHA: PSHAa2393\n',\
            '            ABO: ABO_0668\n',\
            '            BXE: Bxe_B0037 Bxe_C0683 Bxe_C1002 Bxe_C1023\n',\
            '            MPT: Mpe_A2274\n',\
            '            BBA: Bd0910(catD)\n',\
            '            GBE: GbCGDNIH1_0998 GbCGDNIH1_1171\n',\
            '            FNU: FN1345\n',\
            '            RBA: RB13257\n',\
            '            HMA: rrnAC1925(mhpC)\n',\
            '            RCI: RCIX1162 RCIX2396\n']
        
        obs = ko_colon_delimited_parser(test_lines, without_comments = True)
        self.assertEqual(obs['BXE'],['Bxe_B0037','Bxe_C0683', 'Bxe_C1002','Bxe_C1023'])
        self.assertEqual(obs['PHA'],['PSHAa2393'])
        # Check that comments are stripped
        self.assertEqual(obs['BBA'],['Bd0910'])
        
        obs = ko_colon_delimited_parser(test_lines, without_comments = False)
        # Lines without comments shouldn't be affected
        self.assertEqual(obs['BXE'],['Bxe_B0037','Bxe_C0683', 'Bxe_C1002','Bxe_C1023'])
        self.assertEqual(obs['PHA'],['PSHAa2393'])
        # Comments should be preserved 
        self.assertEqual(obs['BBA'],['Bd0910(catD)'])

    def test_is_new_kegg_rec_group(self):
        """_is_new_kegg_rec_group should check for irregular field terminators in KEGG"""
        pass 
        # Handle unusual KEGG fields.

    def test_group_by_end_char(self):
        """group_by_end_char should yield successive lines that end with a given
        char, plus the last group of lines"""
        class_lines=['CLASS       Metabolism; Xenobiotics Biodegradation and Metabolism;\n',\
                    '            gamma-Hexachlorocyclohexane degradation [PATH:ko00361]\n',\
                    '            Metabolism; Xenobiotics Biodegradation and Metabolism;\n',\
                    '             1,2-Dichloroethane degradation [PATH:ko00631]\n']
       
        exp =[['CLASS       Metabolism; Xenobiotics Biodegradation and Metabolism;\n',\
               '            gamma-Hexachlorocyclohexane degradation [PATH:ko00361]\n'],\
              ['            Metabolism; Xenobiotics Biodegradation and Metabolism;\n',\
               '             1,2-Dichloroethane degradation [PATH:ko00631]\n']]
        for i,group in enumerate(group_by_end_char(class_lines)):
            self.assertEqual(group, exp[i])
 
    def test_class_lines_to_fields(self):
        """class_lines_to_fields should split groups of lines for one KO class
        definition"""
        class_lines1=['CLASS      Metabolism; Xenobiotics Biodegradation and Metabolism;\n',\
                    '            gamma-Hexachlorocyclohexane degradation [PATH:ko00361]\n']
        
        class_lines2=['            Metabolism; Xenobiotics Biodegradation and Metabolism;\n',\
                    '            1,2-Dichloroethane degradation [PATH:ko00631]\n']
       
        obs = class_lines_to_fields(class_lines1)
        exp = ('PATH:ko00361',('Metabolism', 'Xenobiotics Biodegradation and Metabolism', 'gamma-Hexachlorocyclohexane degradation'))
        self.assertEqual(obs,exp)
         
        obs = class_lines_to_fields(class_lines2)
        exp = ('PATH:ko00631',('Metabolism', 'Xenobiotics Biodegradation and Metabolism','1,2-Dichloroethane degradation'))
        self.assertEqual(obs,exp)
            

    def test_ko_class_parser(self):
        """ko_class_parser should return fields"""
        class_lines='CLASS       Metabolism; Xenobiotics Biodegradation and Metabolism;\n',\
                    '            gamma-Hexachlorocyclohexane degradation [PATH:ko00361]\n',\
                    '            Metabolism; Xenobiotics Biodegradation and Metabolism;\n',\
                    '             1,2-Dichloroethane degradation [PATH:ko00631]\n'
        exp = [('PATH:ko00361',('Metabolism','Xenobiotics Biodegradation and Metabolism',\
           'gamma-Hexachlorocyclohexane degradation')),\
              ('PATH:ko00631',('Metabolism', 'Xenobiotics Biodegradation and Metabolism', '1,2-Dichloroethane degradation'))]
               
        for i,obs in enumerate(ko_class_parser(class_lines)):
            self.assertEqual(obs,exp[i])

    
    def test_parse_ko(self):
        """parse_ko should parse a ko record into fields """
        lines = TEST_KO_LINES
        r =  parse_ko(lines)
        results = []
        for result in r:
            results.append(result)
        # For each entry we expect a dict
        
        self.assertEqual(results[0]["ENTRY"], "K01559")
        self.assertEqual(results[1]["ENTRY"], "K01560")
        self.assertEqual(results[2]["ENTRY"], "K01561")
        
        self.assertEqual(results[0]["NAME"], "E3.7.1.-")
        self.assertEqual(results[1]["NAME"], "E3.8.1.2")
        self.assertEqual(results[2]["NAME"], "E3.8.1.3")

        self.assertEqual(results[0].get("DEFINITION"), None) #case 1 has no def
        self.assertEqual(results[1]["DEFINITION"],\
          "2-haloacid dehalogenase [EC:3.8.1.2]")
        self.assertEqual(results[2]["DEFINITION"],\
          "haloacetate dehalogenase [EC:3.8.1.3]")
        self.assertEqual(len(results[0]["CLASS"]), 5)
        self.assertEqual(results[0]["CLASS"][4], \
          ('PATH:ko00362', ('Metabolism', \
          'Xenobiotics Biodegradation and Metabolism',\
          'Benzoate degradation via hydroxylation'))) 
       
        self.assertEqual(results[0]["DBLINKS"], \
          {'RN': ['R04488', 'R05100', 'R05363', \
                  'R05365', 'R06371', 'R07515', \
                  'R07831']}) 
        
        self.assertEqual(results[1]["DBLINKS"], \
          {'GO': ['0018784'], 'RN': ['R05287'], 'COG': ['COG1011']}) 
        
        self.assertEqual(results[2]["DBLINKS"], \
          {'GO': ['0018785'], 'RN': ['R05287']})
        
        self.assertEqual(results[0]["GENES"], \
          {'AFM': ['AFUA_4G13070'], 'FNU': ['FN1345'],\
           'GBE': ['GbCGDNIH1_0998', 'GbCGDNIH1_1171'],\
           'PHA': ['PSHAa2393'], \
           'BBA': ['Bd0910'], \
           'ABO': ['ABO_0668'],\
           'MPT': ['Mpe_A2274'],\
           'RCI': ['RCIX1162', 'RCIX2396'], \
           'BXE': ['Bxe_B0037', 'Bxe_C0683', 'Bxe_C1002', 'Bxe_C1023'],\
           'HMA': ['rrnAC1925'], \
           'RBA': ['RB13257']})
        

               
TEST_KO_LINES = ['ENTRY       K01559                      KO\n', '\
NAME        E3.7.1.-\n', '\
PATHWAY     ko00351  1,1,1-Trichloro-2,2-bis(4-chlorophenyl)ethane (DDT)\n', '\
                            degradation\n', '\
            ko00362  Benzoate degradation via hydroxylation\n', '\
            ko00632  Benzoate degradation via CoA ligation\n', '\
            ko00903  Limonene and pinene degradation\n', '\
            ko00930  Caprolactam degradation\n', '\
CLASS       Metabolism; Biosynthesis of Secondary Metabolites; Limonene and\n', '\
            pinene degradation [PATH:ko00903]\n', '\
            Metabolism; Xenobiotics Biodegradation and Metabolism; Caprolactam\n', '\
            degradation [PATH:ko00930]\n', '\
            Metabolism; Xenobiotics Biodegradation and Metabolism;\n', '\
            1,1,1-Trichloro-2,2-bis(4-chlorophenyl)ethane (DDT) degradation\n', '\
            [PATH:ko00351]\n', '\
            Metabolism; Xenobiotics Biodegradation and Metabolism; Benzoate\n', '\
            degradation via CoA ligation [PATH:ko00632]\n', '\
            Metabolism; Xenobiotics Biodegradation and Metabolism; Benzoate\n', '\
            degradation via hydroxylation [PATH:ko00362]\n', '\
DBLINKS     RN: R04488 R05100 R05363 R05365 R06371 R07515 R07831\n', '\
GENES       AFM: AFUA_4G13070\n', '\
            PHA: PSHAa2393\n', '\
            ABO: ABO_0668\n', '\
            BXE: Bxe_B0037 Bxe_C0683 Bxe_C1002 Bxe_C1023\n', '\
            MPT: Mpe_A2274\n', '\
            BBA: Bd0910(catD)\n', '\
            GBE: GbCGDNIH1_0998 GbCGDNIH1_1171\n', '\
            FNU: FN1345\n', '\
            RBA: RB13257\n', '\
            HMA: rrnAC1925(mhpC)\n', '\
            RCI: RCIX1162 RCIX2396\n', '\
///\n', '\
ENTRY       K01560                      KO\n', '\
NAME        E3.8.1.2\n', '\
DEFINITION  2-haloacid dehalogenase [EC:3.8.1.2]\n', '\
PATHWAY     ko00361  gamma-Hexachlorocyclohexane degradation\n', '\
            ko00631  1,2-Dichloroethane degradation\n', '\
CLASS       Metabolism; Xenobiotics Biodegradation and Metabolism;\n', '\
            gamma-Hexachlorocyclohexane degradation [PATH:ko00361]\n', '\
            Metabolism; Xenobiotics Biodegradation and Metabolism;\n', '\
            1,2-Dichloroethane degradation [PATH:ko00631]\n', '\
DBLINKS     RN: R05287\n', '\
            COG: COG1011\n', '\
            GO: 0018784\n', '\
GENES       NCR: NCU03617\n', '\
            ANI: AN5830.2 AN7918.2\n', '\
            AFM: AFUA_2G07750 AFUA_5G14640 AFUA_8G05870\n', '\
            AOR: AO090001000019 AO090003001435 AO090011000921\n', '\
            PST: PSPTO_0247(dehII)\n', '\
            PSP: PSPPH_1747(dehII1) PSPPH_5028(dehII2)\n', '\
            ATU: Atu0797 Atu3405(hadL)\n', '\
            ATC: AGR_C_1458 AGR_L_2834\n', '\
            RET: RHE_CH00996(ypch00330) RHE_PF00342(ypf00173)\n', '\
            MSE: Msed_0732\n', '\
///\n', '\
ENTRY       K01561                      KO\n', '\
NAME        E3.8.1.3\n', '\
DEFINITION  haloacetate dehalogenase [EC:3.8.1.3]\n', '\
PATHWAY     ko00361  gamma-Hexachlorocyclohexane degradation\n', '\
            ko00631  1,2-Dichloroethane degradation\n', '\
CLASS       Metabolism; Xenobiotics Biodegradation and Metabolism;\n', '\
            gamma-Hexachlorocyclohexane degradation [PATH:ko00361]\n', '\
            Metabolism; Xenobiotics Biodegradation and Metabolism;\n', '\
            1,2-Dichloroethane degradation [PATH:ko00631]\n', '\
DBLINKS     RN: R05287\n', '\
            GO: 0018785\n', '\
GENES       RSO: RSc0256(dehH)\n', '\
            REH: H16_A0197\n', '\
            BPS: BPSL0329\n', '\
            BPM: BURPS1710b_0537(dehH)\n', '\
            BPD: BURPS668_0347\n', '\
            STO: ST2570\n', '\
            MSE: Msed_1088\n', '\
///\n']


if __name__=="__main__":
    main()
