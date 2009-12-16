#!/usr/bin/env python
"""Tests of the AAIndex parser.
"""

from cogent.util.unit_test import TestCase, main
from cogent.parse.aaindex import AAIndex1Parser, AAIndex2Parser,\
AAIndexRecord, AAIndex1Record, AAIndex2Record, AAIndex1FromFiles,\
AAIndex2FromFiles

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Greg Caporaso"
__email__ = "caporaso@colorado.edu"
__status__ = "Production"

class test_aaindex1_parser(TestCase):
    """ Tests aindex1_parser class """
    def setUp(self):
        """ Setup some variables """
        self._fake_file = list(fake_file_aaindex1.split('\n'))

        self.AAIndexObjects = AAIndex1FromFiles(self._fake_file)

    def test_init(self):
        """ AAI1: Test that init run w/o error """
        aa1p = AAIndex1Parser()
        
    def test_read_file_as_list(self):
        """AAI1: Test that a file is correctly opened as a list """
        aap = AAIndex1Parser()
        AAIndexObjects = aap(self._fake_file)

    def test_correct_num_of_records(self):
        """AAI1: Test that one object is created per record """
        self.assertEqual(6, len(self.AAIndexObjects))
       
    def test_ID_entries(self):
        """ AAI1: Test ID Entries """
        self.assertEqual(self.AAIndexObjects['ANDN920101'].ID, 'ANDN920101')
        self.assertEqual(self.AAIndexObjects['ARGP820103'].ID, 'ARGP820103')
        self.assertEqual(self.AAIndexObjects['JURD980101'].ID, 'JURD980101')

    def test_single_line_Description_entries(self):
        """ AAI1: Test Single Line Description Entries """
        self.assertEqual(self.AAIndexObjects['ANDN920101'].Description,\
        'alpha-CH chemical shifts (Andersen et al., 1992)')
        self.assertEqual(self.AAIndexObjects['ARGP820103'].Description,\
        'Membrane-buried preference parameters (Argos et al., 1982)')

    def test_multi_line_Description_entries(self):
        """ AAI1: Test Multi Line Description Entries """        
        self.assertEqual(self.AAIndexObjects['JURD980101'].Description,\
        'Modified Kyte-Doolittle hydrophobicity scale (Juretic et al., 1998)')
        
    def test_LITDB_entries(self):
        """ AAI1: Test LITDB Entries """
        self.assertEqual(self.AAIndexObjects['ANDN920101'].LITDBEntryNum,\
        'LIT:1810048b PMID:1575719')
        self.assertEqual(self.AAIndexObjects['ARGP820103'].LITDBEntryNum,\
        'LIT:0901079b PMID:7151796')
        self.assertEqual(self.AAIndexObjects['JURD980101'].LITDBEntryNum,\
        '')

    def test_Authors_entries(self):
        """ AAI1: Test Authors Entries """
        self.assertEqual(self.AAIndexObjects['ANDN920101'].Authors,\
        'Andersen, N.H., Cao, B. and Chen, C.')
        self.assertEqual(self.AAIndexObjects['ARGP820103'].Authors,\
        'Argos, P., Rao, J.K.M. and Hargrave, P.A.')
        self.assertEqual(self.AAIndexObjects['JURD980101'].Authors,\
        'Juretic, D., Lucic, B., Zucic, D. and Trinajstic, N.')                

    def test_mult_line_Title_entries(self):
        """ AAI1: Test Multi Line Title Entries """
        self.assertEqual(self.AAIndexObjects['ANDN920101'].Title,\
        'Peptide/protein structure analysis using the chemical shift index ' +\
        'method: upfield alpha-CH values reveal dynamic helices and aL sites')
        self.assertEqual(self.AAIndexObjects['JURD980101'].Title,\
        'Protein transmembrane structure: recognition and prediction by ' +\
        'using hydrophobicity scales through preference functions')

    def test_sing_line_Title_entries(self):
        """ AAI1: Test Single Line Title Entries """
        self.assertEqual(self.AAIndexObjects['ARGP820103'].Title,\
        'Structural prediction of membrane-bound proteins')     
        
    def test_Citation_entries(self):
        """ AAI1: Test Citation Entries """
        self.assertEqual(self.AAIndexObjects['ANDN920101'].Citation,\
        'Biochem. and Biophys. Res. Comm. 184, 1008-1014 (1992)')
        self.assertEqual(self.AAIndexObjects['ARGP820103'].Citation,\
        'Eur. J. Biochem. 128, 565-575 (1982)')
        self.assertEqual(self.AAIndexObjects['JURD980101'].Citation,\
        'Theoretical and Computational Chemistry, 5, 405-445 (1998)')
         
    def test_Comments_entries(self):
        """ AAI1: Test Comments Entries """
        self.assertEqual(self.AAIndexObjects['ANDN920101'].Comments,\
        '')
        self.assertEqual(self.AAIndexObjects['ARGP820103'].Comments,\
        '')
        self.assertEqual(self.AAIndexObjects['JURD980101'].Comments,\
        '')
        self.assertEqual(self.AAIndexObjects['TSAJ990102'].Comments,\
        '(Cyh 113.7)')

    def test_single_line_Correlating_entries(self):
         """ AAI1: Test single line Correlating Entries """
         self.assertEqual(self.AAIndexObjects['ANDN920101'].\
         Correlating['BUNA790102'], 0.949)

    def test_empty_Correlating_entries(self):
        """ AAI1: Test empty Correlating Entries """       
        self.assertEqual(self.AAIndexObjects['WILM950104'].Correlating, {})

    def test_multi_line_Correlating_entries(self):
         """ AAI1: Test multi line Correlating Entries """
         self.assertEqual(self.AAIndexObjects['ARGP820103'].\
         Correlating['ARGP820102'], 0.961)
         self.assertEqual(self.AAIndexObjects['ARGP820103'].\
         Correlating['MIYS850101'], 0.822)
         self.assertEqual(self.AAIndexObjects['ARGP820103'].\
         Correlating['JURD980101'], 0.800)

         self.assertEqual(self.AAIndexObjects['JURD980101'].\
         Correlating['KYTJ820101'], 0.996)
         self.assertEqual(self.AAIndexObjects['JURD980101'].\
         Correlating['NADH010101'], 0.925)
         self.assertEqual(self.AAIndexObjects['JURD980101'].\
         Correlating['OOBM770101'], -0.903)
                                    
    def test_Data_entries(self):
        """ AAI1: Test Data Entries """
        self.assertEqual(self.AAIndexObjects['ANDN920101'].Data['A'],\
        4.35)
        self.assertEqual(self.AAIndexObjects['ANDN920101'].Data['Q'],\
        4.37)
        self.assertEqual(self.AAIndexObjects['ANDN920101'].Data['V'],\
        3.95)
        self.assertEqual(self.AAIndexObjects['ARGP820103'].Data['A'],\
        1.56)
        self.assertEqual(self.AAIndexObjects['ARGP820103'].Data['Q'],\
        0.51)
        self.assertEqual(self.AAIndexObjects['ARGP820103'].Data['V'],\
        1.14)
        self.assertEqual(self.AAIndexObjects['JURD980101'].Data['A'],\
        1.10)
        self.assertEqual(self.AAIndexObjects['JURD980101'].Data['Q'],\
        -3.68)
        self.assertEqual(self.AAIndexObjects['JURD980101'].Data['V'],\
        4.2)                  

                   
class test_aaindex2_parser(TestCase):
    def setUp(self):
        """ Setup some variables """
        self._fake_file = list(fake_file_aaindex2.split('\n'))
        self.AAIndexObjects = AAIndex2FromFiles(self._fake_file)

    def test_init(self):
        """ AAI2: Test that init run w/o error """
        aa2p = AAIndex2Parser()

    def test_read_file_as_list(self):
        """AAI1: Test that a file is correctly opened as a list """
        aap = AAIndex2Parser()
        AAIndexObjects = aap(self._fake_file)

    def test_correct_num_of_records(self):
        """AAI2: Test that one object is created per record """
        self.assertEqual(6, len(self.AAIndexObjects))

    def test_ID_entries(self):
        """ AAI2: Test ID Entries """
        self.assertEqual(self.AAIndexObjects['ALTS910101'].ID, 'ALTS910101')
        self.assertEqual(self.AAIndexObjects['BENS940103'].ID, 'BENS940103')
        self.assertEqual(self.AAIndexObjects['QUIB020101'].ID, 'QUIB020101')

    def test_Description_entries(self):
        """ AAI2: Test Description Entries """
        self.assertEqual(self.AAIndexObjects['ALTS910101'].Description,\
        'The PAM-120 matrix (Altschul, 1991)')
        self.assertEqual(self.AAIndexObjects['BENS940103'].Description,\
        'Log-odds scoring matrix collected in 74-100 PAM (Benner et al., '+\
        '1994)')
        self.assertEqual(self.AAIndexObjects['QUIB020101'].Description,\
        'STROMA score matrix for the alignment of known distant homologs ' +\
        '(Qian-Goldstein, 2002)')

    def test_LITDB_entries(self):
        """ AAI2: Test LITDB Entries """
        self.assertEqual(self.AAIndexObjects['ALTS910101'].LITDBEntryNum,\
        'LIT:1713145 PMID:2051488')
        self.assertEqual(self.AAIndexObjects['BENS940103'].LITDBEntryNum,\
        'LIT:2023094 PMID:7700864')
        self.assertEqual(self.AAIndexObjects['QUIB020101'].LITDBEntryNum,\
        'PMID:12211027')

    def test_Authors_entries(self):
        """ AAI2: Test Atuthor Entries """
        self.assertEqual(self.AAIndexObjects['ALTS910101'].Authors,\
        'Altschul, S.F.')
        self.assertEqual(self.AAIndexObjects['BENS940103'].Authors,\
        'Benner, S.A., Cohen, M.A. and Gonnet, G.H.')
        self.assertEqual(self.AAIndexObjects['QUIB020101'].Authors,\
        'Qian, B. and Goldstein, R.A.')

    def test_Title_entries(self):
        """ AAI2: Test Title Entries """
        self.assertEqual(self.AAIndexObjects['ALTS910101'].Title,\
        'Amino acid substitution matrices from an information theoretic ' +\
        'perspective')         
        self.assertEqual(self.AAIndexObjects['BENS940103'].Title,\
        'Amino acid substitution during functionally constrained divergent ' +\
        'evolution of protein sequences')
        self.assertEqual(self.AAIndexObjects['QUIB020101'].Title,\
        'Optimization of a new score function for the generation of '+\
        'accurate alignments')
        
    def test_Citation_entries(self):
        """ AAI2: Test citation entries """
        self.assertEqual(self.AAIndexObjects['ALTS910101'].Citation,\
        'J. Mol. Biol. 219, 555-565 (1991)')
        self.assertEqual(self.AAIndexObjects['BENS940103'].Citation,\
        'Protein Engineering 7, 1323-1332 (1994)')
        self.assertEqual(self.AAIndexObjects['QUIB020101'].Citation,\
        'Proteins. 48, 605-610 (2002)')

    def test_Comments_entries(self):
        """ AAI2: Tests null, single line, multi line comments """
        self.assertEqual(self.AAIndexObjects['ALTS910101'].Comments,\
        '')
        self.assertEqual(self.AAIndexObjects['BENS940103'].Comments,\
        'extrapolated to 250 PAM')
        self.assertEqual(self.AAIndexObjects['QUIB020101'].Comments,\
        '')
        self.assertEqual(self.AAIndexObjects['HENS920104'].Comments,\
        '#  Matrix made by matblas from blosum50.iij ' +
        '* #  BLOSUM Clustered Scoring Matrix in 1/3 Bit Units ' +
        '* #  Blocks Database = /data/blocks_5.0/blocks.dat ' +
        '* #  Cluster Percentage: >= 50 ' +
        '* #  Entropy =   0.4808, Expected =  -0.3573')
      
    def test_Data_entries_20x20_LTM(self):
         """ AAI2: correct data entries when 20x20 LTM"""
         self.assertEqual(self.AAIndexObjects['ALTS910101'].Data['A']['A'],\
         3.)
         self.assertEqual(self.AAIndexObjects['ALTS910101'].Data['Y']['R'],\
         -6.)
         self.assertEqual(self.AAIndexObjects['ALTS910101'].Data['V']['V'],\
         5.)
         self.assertEqual(self.AAIndexObjects['BENS940103'].Data['A']['A'],\
         2.4)
         self.assertEqual(self.AAIndexObjects['BENS940103'].Data['Y']['R'],\
         -2.0)
         self.assertEqual(self.AAIndexObjects['BENS940103'].Data['V']['V'],\
         3.4)
         self.assertEqual(self.AAIndexObjects['QUIB020101'].Data['A']['A'],\
         2.5)
         self.assertEqual(self.AAIndexObjects['QUIB020101'].Data['Y']['R'],\
         -0.9)
         self.assertEqual(self.AAIndexObjects['QUIB020101'].Data['V']['V'],\
         4.2)
         
    def test_Data_entries_20x20_Square(self):
        """ AAI2: correct data entries when 20x20 squ matrix """
        self.assertEqual(self.AAIndexObjects['HENS920104'].Data['V']['Y'],\
        -1)
        self.assertEqual(self.AAIndexObjects['HENS920104'].Data['Q']['A'],\
        -1)
        self.assertEqual(self.AAIndexObjects['HENS920104'].Data['N']['N'],\
        7)

    def test_Data_entries_with_abnormal_fields(self):
        """ AAI2: test correct data entries when more than std fields present

            Some entires in AAIndex2 have more that 20 fields, this tests that
            that data is corrected parsed and identified.
        """
        # There are no entries that fit this category that are square
        # matrices, which is all we are concerned with at this point,
        # so this method should just serve as a reminder to test this
        # when we begin parsing data other than square matrices.
        pass
        
    def test_Data_entries_21x21_LTM(self):
        """ AAI2: correct data entries when 21x21 LTM"""
        self.assertEqual(self.AAIndexObjects['KOSJ950101'].Data['-']['-'],\
        55.7)
        self.assertEqual(self.AAIndexObjects['KOSJ950101'].Data['Y']['-'],\
        0.3)
        self.assertEqual(self.AAIndexObjects['KOSJ950101'].Data['N']['R'],\
        3.0)

    def test_Data_entries_22x21_square(self):
        """ AAI2: correct data entries when 22x21 square matrix """
        # It's not really a sqaure matrix, but it's fully populated ...
        self.assertEqual(self.AAIndexObjects['OVEJ920102'].Data['J']['D'],\
        0.001)
        self.assertEqual(self.AAIndexObjects['OVEJ920102'].Data['-']['I'],\
        0.022)
        self.assertEqual(self.AAIndexObjects['OVEJ920102'].Data['D']['E'],\
        0.109)

        
class AAIndexRecordTests(TestCase):
    """ AAIR: Tests AAIndexRecord class """

    def setUp(self):
        self.id = "5"
        self.description = "Some Info"
        self.LITDB_entry_num = "25"
        self.authors = "Greg"
        self.title = "A test"
        self.citation = "something"
        self.comments = "This is a test, this is only a test"
        self.data = {}


class AAIndex1RecordTests(AAIndexRecordTests):
    """ AAIR1: Tests AAIndex1Records class """

    def setUp(self):
        AAIndexRecordTests.setUp(self)

        self.correlating = [0.987, 0.783, 1., 0]

        values = []
        keys = 'ARNDCQEGHILKMFPSTWYV'
        for i in range(20):
            values += [float(i) + 0.15]
        self.data = dict(zip(keys,values))

        self.aar = AAIndex1Record(self.id, self.description,\
                self.LITDB_entry_num, self.authors, self.title,\
                self.citation, self.comments, self.correlating, self.data)

    def test_init(self):
        """ AAIR1: Tests init method returns with no errors"""

        test_aar = AAIndex1Record(self.id, self.description,\
                self.LITDB_entry_num, self.authors, self.title,\
                self.citation, self.comments, self.correlating, self.data)

    def test_general_init_data(self):
        """ AAIR1: Tests init correctly initializes data"""

        self.assertEqual(self.aar.ID, str(self.id))
        self.assertEqual(self.aar.Description, str(self.description))
        self.assertEqual(self.aar.LITDBEntryNum,\
                str(self.LITDB_entry_num))
        self.assertEqual(self.aar.Authors, str(self.authors))
        self.assertEqual(self.aar.Title, str(self.title))
        self.assertEqual(self.aar.Citation, str(self.citation))
        self.assertEqual(self.aar.Comments, str(self.comments))
        self.assertEqual(self.aar.Correlating, self.correlating)
        self.assertEqual(self.aar.Data,self.data)

    def test_toSquareDistanceMatrix(self):
        """ AAIR1: Tests that _toSquareDistanceMatrix runs without returning an error """
        square = self.aar._toSquareDistanceMatrix()

    def test_toSquareDistanceMatrix_data_integrity_diagonal(self):
        """ AAIR1: Tests that diag = 0 when square matrix is built """
        square = self.aar._toSquareDistanceMatrix()
        # Test diagonal
        keys = 'ARNDCQEGHILKMFPSTWYV'
        for k in keys:
            self.assertEqual(square[k][k], 0.)

    def test_toSquareDistanceMatrix_data_integrity(self):
        """ AAIR1: Tests that _toSquareDistanceMatrix works right w/o stops """
        square = self.aar._toSquareDistanceMatrix()
        self.assertFloatEqualAbs(square['R']['A'], square['A']['R'])
        self.assertFloatEqualAbs(square['A']['R'], 1.)
        self.assertFloatEqualAbs(square['D']['N'], square['N']['D'])
        self.assertFloatEqualAbs(square['D']['N'], 1.)
        self.assertFloatEqualAbs(square['A']['C'], square['C']['A'])
        self.assertFloatEqualAbs(square['A']['C'], 4.)
        self.assertFloatEqualAbs(square['V']['A'], square['A']['V'])
        self.assertFloatEqualAbs(square['V']['A'], 19.)
        self.assertFloatEqualAbs(square['V']['Y'], square['Y']['V'])
        self.assertFloatEqualAbs(square['V']['Y'], 1.)

    def test_toSquareDistanceMatrix_data_integrity_w_stops(self):
        """ AAIR1: Tests that _toSquareDistanceMatrix works right w/ stops """
        square = self.aar._toSquareDistanceMatrix(include_stops=1)
        self.assertFloatEqualAbs(square['R']['A'], square['A']['R'])
        self.assertFloatEqualAbs(square['A']['R'], 1.)
        self.assertFloatEqualAbs(square['D']['N'], square['N']['D'])
        self.assertFloatEqualAbs(square['D']['N'], 1.)
        self.assertFloatEqualAbs(square['A']['C'], square['C']['A'])
        self.assertFloatEqualAbs(square['A']['C'], 4.)
        self.assertFloatEqualAbs(square['V']['A'], square['A']['V'])
        self.assertFloatEqualAbs(square['V']['A'], 19.)
        self.assertFloatEqualAbs(square['V']['Y'], square['Y']['V'])
        self.assertFloatEqualAbs(square['V']['Y'], 1.)
        self.assertFloatEqualAbs(square['V']['*'], None)
        self.assertFloatEqualAbs(square['*']['Y'], None)
        self.assertFloatEqualAbs(square['*']['*'], None)
        self.assertFloatEqualAbs(square['*']['R'], None)

    def test_toDistanceMatrix(self):
        """ AAIR1: Tests that toDistanceMatrix functions as expected """
        dm = self.aar.toDistanceMatrix()
        self.assertFloatEqualAbs(dm['R']['A'], dm['A']['R'])
        self.assertFloatEqualAbs(dm['A']['R'], 1.)
        self.assertFloatEqualAbs(dm['D']['N'], dm['N']['D'])
        self.assertFloatEqualAbs(dm['D']['N'], 1.)
        self.assertFloatEqualAbs(dm['A']['C'], dm['C']['A'])
        self.assertFloatEqualAbs(dm['A']['C'], 4.)
        self.assertFloatEqualAbs(dm['V']['A'], dm['A']['V'])
        self.assertFloatEqualAbs(dm['V']['A'], 19.)
        self.assertFloatEqualAbs(dm['V']['Y'], dm['Y']['V'])
        self.assertFloatEqualAbs(dm['V']['Y'], 1.)

    def test_toDistanceMatrix_w_stops(self):
        """ AAIR1: Tests that toDistanceMatrix works right w/ stops """
        square = self.aar.toDistanceMatrix(include_stops=1)
        self.assertFloatEqualAbs(square['R']['A'], square['A']['R'])
        self.assertFloatEqualAbs(square['A']['R'], 1.)
        self.assertFloatEqualAbs(square['D']['N'], square['N']['D'])
        self.assertFloatEqualAbs(square['D']['N'], 1.)
        self.assertFloatEqualAbs(square['A']['C'], square['C']['A'])
        self.assertFloatEqualAbs(square['A']['C'], 4.)
        self.assertFloatEqualAbs(square['V']['A'], square['A']['V'])
        self.assertFloatEqualAbs(square['V']['A'], 19.)
        self.assertFloatEqualAbs(square['V']['Y'], square['Y']['V'])
        self.assertFloatEqualAbs(square['V']['Y'], 1.)
        self.assertFloatEqualAbs(square['V']['*'], None)
        self.assertFloatEqualAbs(square['*']['Y'], None)
        self.assertFloatEqualAbs(square['*']['*'], None)
        self.assertFloatEqualAbs(square['*']['R'], None)

class AAIndex2RecordTests(AAIndexRecordTests):
    """ AAIR2: Tests AAIndex2Records class """

    def setUp(self):
        AAIndexRecordTests.setUp(self)

        # Build LTM data
        values = range(210)
        keys = 'ARNDCQEGHILKMFPSTWYV'

        self.LTMdata = dict.fromkeys(keys)

        i = 0
        for r in keys:
            new_row = dict.fromkeys(keys)
            for c in keys:
                if keys.find(c) <= keys.find(r):
                    new_row[c] = values[i]
                    i +=1
            self.LTMdata[r] = new_row



        self.aarLTM = AAIndex2Record(self.id, self.description,\
                self.LITDB_entry_num, self.authors, self.title,\
                self.citation, self.comments, self.LTMdata)

        # Build Square matrix data
        values = range(400)

        self.SQUdata = dict.fromkeys(keys)

        i = 0
        for r in keys:
            new_row = dict.fromkeys(keys)
            for c in keys:
                new_row[c] = values[i]
                i +=1
            self.SQUdata[r] = new_row

        self.aarSquare = AAIndex2Record(self.id, self.description,\
                self.LITDB_entry_num, self.authors, self.title,\
                self.citation, self.comments, self.SQUdata)

    def test_init(self):
        """ AAIR2: Tests init method returns with no errors"""

        test_aar = AAIndex2Record(self.id, self.description,\
                self.LITDB_entry_num, self.authors, self.title,\
                self.citation, self.comments, self.SQUdata)

    def test_init_data(self):
        """ AAIR2: Tests init correctly initializes data"""

        self.assertEqual(self.aarLTM.ID, str(self.id))
        self.assertEqual(self.aarLTM.Description, str(self.description))
        self.assertEqual(self.aarLTM.LITDBEntryNum,\
                str(self.LITDB_entry_num))
        self.assertEqual(self.aarLTM.Authors, str(self.authors))
        self.assertEqual(self.aarLTM.Title, str(self.title))
        self.assertEqual(self.aarLTM.Citation, str(self.citation))
        self.assertEqual(self.aarLTM.Comments, str(self.comments))

#    def test_matrix_values_col_by_row(self):
 #       """ Tests that keys and values correctly correspond in data LTM
#

#
 #           Also tests that reverse keys are same as forward keys.
#
 #       """
#
 #       data_matrix = self.aarLTM.Data
  #      self.assertEqual(data_matrix['A']['A'], 0)
  #      self.assertEqual(data_matrix['A']['R'], 1)
   #     self.assertEqual(data_matrix['R']['R'], 2)
    #    self.assertEqual(data_matrix['C']['H'], 40)
#        self.assertEqual(data_matrix['I']['M'], 87)
 #       self.assertEqual(data_matrix['D']['P'], 108)
 #       self.assertEqual(data_matrix['W']['V'], 207)
  #      self.assertEqual(data_matrix['Y']['V'], 208)
  #      self.assertEqual(data_matrix['V']['V'], 209)

 #   def test_LTM_values_row_by_col(self):
  #      """ Tests that keys are correctly linked to values in a LTM
#
 #           This tests that some random places hold the correct values.
  #          These are some randomly selected keys with hand calculated
   #         values.  Also included are the extreme values. Technically if
#            the first and last three are correct all values should be
 #           correct.
#
 #       """
  #      data_matrix = self.aarLTM.Data
   #     self.assertEqual(data_matrix['R']['A'], 1)
   #     self.assertEqual(data_matrix['H']['C'], 40)
   #     self.assertEqual(data_matrix['M']['I'], 87)
   #     self.assertEqual(data_matrix['P']['D'], 108)
   #     self.assertEqual(data_matrix['V']['W'], 207)
   #     self.assertEqual(data_matrix['V']['Y'], 208)
   #     self.assertEqual(data_matrix['A']['A'], 0)
    #    self.assertEqual(data_matrix['R']['R'], 2)
     #   self.assertEqual(data_matrix['V']['V'], 209)

    def test_Square_Matrix_values_row_by_col(self):
        """ AAIR2: Tests that key -> value pair integrity in Square matrix
        """
        data_matrix = self.aarSquare.Data
        self.assertEqual(data_matrix['R']['A'], 20)
        #self.assertEqual(data_matrix['H']['C'], 40)
        #self.assertEqual(data_matrix['M']['I'], 87)
        #self.assertEqual(data_matrix['P']['D'], 108)

        #self.assertEqual(data_matrix['V']['W'], 207)
        #self.assertEqual(data_matrix['V']['Y'], 208)
        self.assertEqual(data_matrix['A']['A'], 0)
        self.assertEqual(data_matrix['R']['R'], 21)
        self.assertEqual(data_matrix['V']['V'], 399)

    def test_toSquareDistanceMatrix_data_integrity(self):
        """ AAIR2: Tests that _toSquareDistanceMatrix works right w/o stops

        """
        square = self.aarSquare._toSquareDistanceMatrix()
        self.assertEqual(square['R']['A'], 20)
        self.assertEqual(square['A']['A'], 0)
        self.assertEqual(square['R']['R'], 21)
        self.assertEqual(square['V']['V'], 399)

    def test_toSquareDistanceMatrix_data_integrity_w_stops(self):
        """ AAIR2: Tests that _toSquareDistanceMatrix works right with stops
        """
        square = self.aarSquare._toSquareDistanceMatrix(include_stops=1)
        self.assertEqual(square['R']['A'], 20)
        self.assertEqual(square['A']['A'], 0)
        self.assertEqual(square['R']['R'], 21)
        self.assertEqual(square['V']['V'], 399)
        self.assertEqual(square['V']['*'], None)
        self.assertEqual(square['*']['Y'], None)
        self.assertEqual(square['*']['*'], None)
        self.assertEqual(square['*']['R'], None)

# Data for parser tests
        
fake_file_aaindex1 =\
"""
H ANDN920101
D alpha-CH chemical shifts (Andersen et al., 1992)
R LIT:1810048b PMID:1575719
A Andersen, N.H., Cao, B. and Chen, C.
T Peptide/protein structure analysis using the chemical shift index method:
  upfield alpha-CH values reveal dynamic helices and aL sites
J Biochem. and Biophys. Res. Comm. 184, 1008-1014 (1992)
C BUNA790102    0.949
I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    4.35    4.38    4.75    4.76    4.65    4.37    4.29    3.97    4.63    3.95
    4.17    4.36    4.52    4.66    4.44    4.50    4.35    4.70    4.60    3.95
//
H ARGP820101
D Hydrophobicity index (Argos et al., 1982)
R LIT:0901079b PMID:7151796
A Argos, P., Rao, J.K.M. and Hargrave, P.A.
T Structural prediction of membrane-bound proteins
J Eur. J. Biochem. 128, 565-575 (1982)
C JOND750101    1.000  SIMZ760101    0.967  GOLD730101    0.936
  TAKK010101    0.906  MEEJ810101    0.891  CIDH920105    0.867
  LEVM760106    0.865  CIDH920102    0.862  MEEJ800102    0.855
  MEEJ810102    0.853  CIDH920103    0.827  PLIV810101    0.820
  CIDH920104    0.819  LEVM760107    0.806  NOZY710101    0.800
  PARJ860101   -0.835  WOLS870101   -0.838  BULH740101   -0.854
I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    0.61    0.60    0.06    0.46    1.07      0.    0.47    0.07    0.61    2.22
    1.53    1.15    1.18    2.02    1.95    0.05    0.05    2.65    1.88    1.32
//
H TSAJ990102
D Volumes not including the crystallographic waters using the ProtOr (Tsai et
  al., 1999)
R PMID:10388571
A Tsai, J., Taylor, R., Chothia, C. and Gerstein, M.
T The packing density in proteins: standard radii and volumes
J J Mol Biol. 290, 253-266 (1999)
* (Cyh 113.7)
C TSAJ990101    1.000  CHOC750101    0.996  BIGC670101    0.992
  GOLD730102    0.991  KRIW790103    0.987  FAUJ880103    0.985
  GRAR740103    0.978  CHAM820101    0.978  CHOC760101    0.972
  FASG760101    0.940  LEVM760105    0.928  LEVM760102    0.918
  ROSG850101    0.909  DAWD720101    0.905  CHAM830106    0.896
  FAUJ880106    0.882  RADA880106    0.864  LEVM760107    0.861
  LEVM760106    0.841  RADA880103   -0.879
I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    90.0   194.0   124.7   117.3   103.3   149.4   142.2    64.9   160.0   163.9
   164.0   167.3   167.0   191.9   122.9    95.4   121.5   228.2   197.0   139.0
//
H JURD980101
D Modified Kyte-Doolittle hydrophobicity scale (Juretic et al., 1998)
R
A Juretic, D., Lucic, B., Zucic, D. and Trinajstic, N.
T Protein transmembrane structure: recognition and prediction by using
  hydrophobicity scales through preference functions
J Theoretical and Computational Chemistry, 5, 405-445 (1998)
C KYTJ820101    0.996  CHOC760103    0.967  NADH010102    0.931
  JANJ780102    0.928  NADH010101    0.925  EISD860103    0.901
  DESM900102    0.900  NADH010103    0.900  EISD840101    0.895
  RADA880101    0.893  MANP780101    0.887  WOLR810101    0.881
  PONP800103    0.879  JANJ790102    0.879  NADH010104    0.873
  CHOC760104    0.870  PONP800102    0.869  JANJ790101    0.868
  MEIH800103    0.861  PONP800101    0.858  NAKH920108    0.858
  RADA880108    0.857  PONP800108    0.856  ROSG850102    0.854
  PONP930101    0.849  RADA880107    0.842  BIOV880101    0.840
  MIYS850101    0.837  FAUJ830101    0.833  CIDH920104    0.832
  DESM900101    0.829  WARP780101    0.827  KANM800104    0.826
  LIFS790102    0.824  RADA880104    0.824  NADH010105    0.821
  NISK800101    0.816  NISK860101    0.808  BIOV880102    0.805
  ARGP820102    0.802  ARGP820103    0.800  VHEG790101   -0.814
  KRIW790101   -0.824  CHOC760102   -0.851  ROSM880101   -0.851
  MONM990101   -0.853  JANJ780103   -0.853  RACS770102   -0.855
  PRAM900101   -0.862  JANJ780101   -0.862  GUYH850101   -0.864
  GRAR740102   -0.864  MEIH800102   -0.879  KUHL950101   -0.884
  ROSM880102   -0.894  OOBM770101   -0.903
I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    1.10   -5.10   -3.50   -3.60    2.50   -3.68   -3.20   -0.64   -3.20    4.50
    3.80   -4.11    1.90    2.80   -1.90   -0.50   -0.70   -0.46    -1.3     4.2
//
H WILM950104
D Hydrophobicity coefficient in RP-HPLC, C18 with 0.1%TFA/2-PrOH/MeCN/H2O
  (Wilce et al. 1995)
R
A Wilce, M.C., Aguilar, M.I. and Hearn, M.T.
T Physicochemical basis of amino acid hydrophobicity scales: evaluation of four
  new scales of amino acid hydrophobicity coefficients derived from RP-HPLC of
  peptides
J Anal Chem. 67, 1210-1219 (1995)
C
I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
   -2.34    1.60    2.81   -0.48    5.03    0.16    1.30   -1.06   -3.00    7.26
    1.09    1.56    0.62    2.57   -0.15    1.93    0.19    3.59   -2.58    2.06
//
H ARGP820103
D Membrane-buried preference parameters (Argos et al., 1982)
R LIT:0901079b PMID:7151796
A Argos, P., Rao, J.K.M. and Hargrave, P.A.
T Structural prediction of membrane-bound proteins
J Eur. J. Biochem. 128, 565-575 (1982)
C ARGP820102    0.961  MIYS850101    0.822  NAKH900106    0.810
  EISD860103    0.810  KYTJ820101    0.806  JURD980101    0.800
I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    1.56    0.45    0.27    0.14    1.23    0.51    0.23    0.62    0.29    1.67
    2.93    0.15    2.96    2.03    0.76    0.81    0.91    1.08    0.68    1.14
//
"""

fake_file_aaindex2 =\
"""
H ALTS910101
D The PAM-120 matrix (Altschul, 1991)
R LIT:1713145 PMID:2051488
A Altschul, S.F.
T Amino acid substitution matrices from an information theoretic perspective
J J. Mol. Biol. 219, 555-565 (1991)
M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV
      3.
     -3.      6.
      0.     -1.      4.
      0.     -3.      2.      5.
     -3.     -4.     -5.     -7.      9.
     -1.      1.      0.      1.     -7.      6.
      0.     -3.      1.      3.     -7.      2.      5.
      1.     -4.      0.      0.     -5.     -3.     -1.      5.
     -3.      1.      2.      0.     -4.      3.     -1.     -4.      7.
     -1.     -2.     -2.     -3.     -3.     -3.     -3.     -4.     -4.      6.
     -3.     -4.     -4.     -5.     -7.     -2.     -4.     -5.     -3.      1.      5.
     -2.      2.      1.     -1.     -7.      0.     -1.     -3.     -2.     -2.     -4.      5.
     -2.     -1.     -3.     -4.     -6.     -1.     -4.     -4.     -4.      1.      3.      0.      8.
     -4.     -4.     -4.     -7.     -6.     -6.     -6.     -5.     -2.      0.      0.     -6.     -1.      8.
      1.     -1.     -2.     -2.     -3.      0.     -1.     -2.     -1.     -3.     -3.     -2.     -3.     -5.      6.
      1.     -1.      1.      0.     -1.     -2.     -1.      1.     -2.     -2.     -4.     -1.     -2.     -3.      1.      3.
      1.     -2.      0.     -1.     -3.     -2.     -2.     -1.     -3.      0.     -3.     -1.     -1.     -4.     -1.      2.      4.
     -7.      1.     -5.     -8.     -8.     -6.     -8.     -8.     -5.     -7.     -5.     -5.     -7.     -1.     -7.     -2.     -6.     12.
     -4.     -6.     -2.     -5.     -1.     -5.     -4.     -6.     -1.     -2.     -3.     -6.     -4.      4.     -6.     -3.     -3.     -1.      8.
      0.     -3.     -3.     -3.     -2.     -3.     -3.     -2.     -3.      3.      1.     -4.      1.     -3.     -2.     -2.      0.     -8.     -3.      5.
//
H BENS940103
D Log-odds scoring matrix collected in 74-100 PAM (Benner et al., 1994)
R LIT:2023094 PMID:7700864
A Benner, S.A., Cohen, M.A. and Gonnet, G.H.
T Amino acid substitution during functionally constrained divergent
  evolution of protein sequences
J Protein Engineering 7, 1323-1332 (1994)
* extrapolated to 250 PAM
M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV
     2.4
    -0.8     4.8
    -0.2     0.3     3.6
    -0.3    -0.5     2.2     4.8
     0.3    -2.2    -1.8    -3.2    11.8
    -0.3     1.6     0.7     0.8    -2.6     3.0
    -0.1     0.3     1.0     2.9    -3.2     1.7     3.7
     0.6    -1.0     0.4     0.2    -2.0    -1.1    -0.5     6.6
    -1.0     1.0     1.2     0.4    -1.3     1.4     0.2    -1.6     6.1
    -0.8    -2.6    -2.8    -3.9    -1.2    -2.0    -2.9    -4.3    -2.3     4.0
    -1.4    -2.4    -3.1    -4.2    -1.6    -1.7    -3.1    -4.6    -1.9     2.8     4.2
    -0.4     2.9     0.9     0.4    -2.9     1.7     1.2    -1.1     0.6    -2.3    -2.4     3.4
    -0.8    -1.8    -2.2    -3.2    -1.2    -1.0    -2.2    -3.5    -1.5     2.6     2.9    -1.5     4.5
    -2.6    -3.5    -3.2    -4.7    -0.7    -2.8    -4.3    -5.4     0.0     0.9     2.1    -3.6     1.3     7.2
     0.4    -1.0    -1.0    -1.0    -3.1    -0.2    -0.7    -1.7    -1.0    -2.6    -2.2    -0.8    -2.4    -3.8     7.5
     1.1    -0.2     0.9     0.4     0.1     0.1     0.1     0.4    -0.3    -1.8    -2.2     0.0    -1.4    -2.6     0.5     2.1
     0.7    -0.3     0.4    -0.2    -0.6    -0.1    -0.2    -1.0    -0.5    -0.3    -1.1     0.1    -0.4    -2.2     0.1     1.4     2.5
    -4.1    -1.6    -4.0    -5.5    -0.9    -2.8    -4.7    -4.1    -1.0    -2.3    -0.9    -3.6    -1.3     3.0    -5.2    -3.4    -3.7    14.7
    -2.6    -2.0    -1.4    -2.8    -0.4    -1.8    -3.0    -4.3     2.5    -1.0    -0.1    -2.4    -0.5     5.3    -3.4    -1.9    -2.1     3.6     8.1
     0.1    -2.2    -2.2    -2.9    -0.2    -1.7    -2.1    -3.1    -2.1     3.2     1.9    -1.9     1.8     0.1    -1.9    -1.0     0.2    -2.9    -1.4     3.4
//
H QUIB020101
D STROMA score matrix for the alignment of known distant homologs
  (Qian-Goldstein, 2002)
R PMID:12211027
A Qian, B. and Goldstein, R.A.
T Optimization of a new score function for the generation of accurate
  alignments
J Proteins. 48, 605-610 (2002)
M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV
     2.5
     0.2     5.2
     1.1     0.7     2.5
       1     0.1     3.3     5.3
     1.2    -1.3    -1.9    -3.1    11.5
    -0.1       2     1.9     1.1    -2.5     3.6
     1.2     1.9     2.3     3.2    -2.4     1.7     3.7
     1.4    -0.2     0.7     0.9    -1.3    -0.3     0.5     7.5
    -1.4     1.5     1.4     0.5    -1.7     1.4     0.3    -1.7     6.8
     0.3    -1.9    -2.4    -2.9    -3.2    -0.9    -3.1    -3.7    -1.8     4.5
    -0.2    -1.5    -2.4    -3.4    -1.6    -1.2    -1.5    -3.8    -2.4     3.4     5.2
    -0.2     3.4     1.6     1.4      -3     2.2     1.2     0.4     1.1    -1.5      -2     3.9
    -0.2    -1.4    -2.1    -2.8    -1.3    -0.6      -2    -3.8    -0.8     2.2     3.1    -0.5     5.4
    -1.6    -3.2    -2.5    -3.7    -0.8    -1.7   -13.7    -4.7    -0.9     2.2     3.7    -2.8     1.7       7
     0.7    -0.6    -0.1    -0.2    -3.6       1       0    -0.8    -2.1    -2.4    -1.4     0.2    -1.9    -4.1     8.1
     1.7     0.2     1.4     1.7     0.7     0.9     1.1     1.6    -0.1    -1.1    -0.8     1.4    -1.1    -2.5       2     2.8
     1.7     0.2     1.4     0.1     0.3    -0.1     1.6    -0.6    -0.2       0     0.3       1    -0.3    -0.8     1.1     2.6     0.4
    -3.3    -1.5      -4    -5.7    -0.5    -2.9    -4.7    -4.2    -1.2    -1.8    -1.2      -3    -0.6     3.7      -5    -2.8    -2.9    14.9
    -1.8    -0.9    -0.8    -2.9    -0.3    -1.5    -2.2    -4.8     2.9     0.2     0.8    -1.5     0.5     5.2    -3.3    -0.9    -0.8     4.9     8.1
     1.9    -2.8    -0.9    -2.5     0.7    -1.5    -1.3    -1.4    -2.5     4.5     3.4      -1     1.7     0.9    -1.1      -3     1.5    -2.5     0.3     4.2
//
H HENS920104
D BLOSUM50 substitution matrix (Henikoff-Henikoff, 1992)
R LIT:1902106 PMID:1438297
A Henikoff, S. and Henikoff, J.G.
T Amino acid substitution matrices from protein blocks
J Proc. Natl. Acad. Sci. USA 89, 10915-10919 (1992)
* #  Matrix made by matblas from blosum50.iij
* #  BLOSUM Clustered Scoring Matrix in 1/3 Bit Units
* #  Blocks Database = /data/blocks_5.0/blocks.dat
* #  Cluster Percentage: >= 50
* #  Entropy =   0.4808, Expected =  -0.3573
M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV
       5      -2      -1      -2      -1      -1      -1       0      -2      -1      -2      -1      -1      -3      -1       1       0      -3      -2       0
      -2       7      -1      -2      -4       1       0      -3       0      -4      -3       3      -2      -3      -3      -1      -1      -3      -1      -3
      -1      -1       7       2      -2       0       0       0       1      -3      -4       0      -2      -4      -2       1       0      -4      -2      -3
      -2      -2       2       8      -4       0       2      -1      -1      -4      -4      -1      -4      -5      -1       0      -1      -5      -3      -4
      -1      -4      -2      -4      13      -3      -3      -3      -3      -2      -2      -3      -2      -2      -4      -1      -1      -5      -3      -1
      -1       1       0       0      -3       7       2      -2       1      -3      -2       2       0      -4      -1       0      -1      -1      -1      -3
      -1       0       0       2      -3       2       6      -3       0      -4      -3       1      -2      -3      -1      -1      -1      -3      -2      -3
       0      -3       0      -1      -3      -2      -3       8      -2      -4      -4      -2      -3      -4      -2       0      -2      -3      -3      -4
      -2       0       1      -1      -3       1       0      -2      10      -4      -3       0      -1      -1      -2      -1      -2      -3       2      -4
      -1      -4      -3      -4      -2      -3      -4      -4      -4       5       2      -3       2       0      -3      -3      -1      -3      -1       4
      -2      -3      -4      -4      -2      -2      -3      -4      -3       2       5      -3       3       1      -4      -3      -1      -2      -1       1
      -1       3       0      -1      -3       2       1      -2       0      -3      -3       6      -2      -4      -1       0      -1      -3      -2      -3
      -1      -2      -2      -4      -2       0      -2      -3      -1       2       3      -2       7       0      -3      -2      -1      -1       0       1
      -3      -3      -4      -5      -2      -4      -3      -4      -1       0       1      -4       0       8      -4      -3      -2       1       4      -1
      -1      -3      -2      -1      -4      -1      -1      -2      -2      -3      -4      -1      -3      -4      10      -1      -1      -4      -3      -3
       1      -1       1       0      -1       0      -1       0      -1      -3      -3       0      -2      -3      -1       5       2      -4      -2      -2
       0      -1       0      -1      -1      -1      -1      -2      -2      -1      -1      -1      -1      -2      -1       2       5      -3      -2       0
      -3      -3      -4      -5      -5      -1      -3      -3      -3      -3      -2      -3      -1       1      -4      -4      -3      15       2      -3
      -2      -1      -2      -3      -3      -1      -2      -3       2      -1      -1      -2       0       4      -3      -2      -2       2       8      -1
       0      -3      -3      -4      -1      -3      -3      -4      -4       4       1      -3       1      -1      -3      -2       0      -3      -1       5
//
H KOSJ950101
D Context-dependent optimal substitution matrices for exposed helix
  (Koshi-Goldstein, 1995)
R LIT:2124140 PMID:8577693
A Koshi, J.M. and Goldstein, R.A.
T Context-dependent optimal substitution matrices.
J Protein Engineering 8, 641-645 (1995)
M rows = -ARNDCQEGHILKMFPSTWYV, cols = -ARNDCQEGHILKMFPSTWYV
    55.7
     3.0     3.0
     3.0     3.0     0.4
     0.1     3.0     3.0     2.1
     3.0     3.0     3.0     0.1     1.9
     2.2     2.4     3.0     0.8     1.3     3.0
    25.6    47.2     1.5     1.0     0.7     0.3     1.9
     2.3     4.3     0.6     0.2     2.0     0.8     0.1     0.3
     3.1     2.8     3.7     0.4     0.1     2.0    14.8     0.9    62.7
     1.3     0.4     0.3     4.6     0.3     0.1     1.9     0.5     2.2     5.1
     0.6     0.2     0.4     1.9     1.5     0.4     0.2     0.3    15.2     0.2     0.5
    48.2     3.3     0.1     3.2     4.9     0.1     1.7     1.7     1.4     3.0     0.6     1.0
     0.1     9.7     2.7     0.7     1.1     1.5    15.9     3.9     1.4     7.3    52.1     0.3     0.9
    11.0     2.0     0.4     0.6     0.1     0.6     0.5     0.1     0.6     2.9     0.1     0.1     0.1     0.1
     9.4     1.5     0.1     1.5     1.6    73.6     0.1     2.6     0.1     0.1     2.1     4.0     0.1     0.1     0.8
     0.7     0.3     2.2     0.1     0.1     0.1     0.1     8.4     5.7     2.0     4.5     0.3    47.5     8.2     0.9     1.6
     0.1     3.4     7.8     0.5     0.1     0.7     5.3     2.2     0.2     0.7     0.5     5.2     5.3     1.0     1.5     8.6     0.1
     4.9    56.8     1.5     1.0     0.3     0.9     5.8     0.1     0.2     1.6     2.1     2.4     0.2     0.1     1.1    20.2     2.0     1.2
     2.3     3.3     0.1     0.4     0.1       6     4.8     0.8     0.1     0.1     1.4     0.3     0.6     0.1     1.2     0.6     0.1     0.5    13.3
     0.3     4.7     7.5     1.8     0.1     4.4     0.7     0.1    56.9     0.6     0.1     2.3     1.2     2.2     0.1     0.1     0.1     0.1     4.4     0.1
    18.4     0.1     0.1     0.1     0.1     0.1     0.4     0.1     0.1     0.1       5     2.6    10.8     1.2     3.5     1.3     0.1     0.1     3.4     0.1     0.1
//
H OVEJ920102
D Environment-specific amino acid substitution matrix for alpha residues
  (Overington et al., 1992)
R LIT:1811128 PMID:1304904
A Overington, J., Donnelly, D., Johnson, M.S., Sali, A. and Blundell, T.L.
T Environment-specific amino acid substitution tables: tertiary templates
  and prediction of protein folds
J Protein Science 1, 216-226 (1992)
M rows = ACDEFGHIKLMNPQRSTVWYJ-, cols = ACDEFGHIKLMNPQRSTVWYJ
   0.355   0.007   0.090   0.100   0.050   0.177   0.037   0.077   0.096   0.056   0.081   0.103   0.106   0.090   0.088   0.163   0.120   0.098   0.065   0.036   0.252
   0.001   0.901   0.000   0.000   0.000   0.000   0.000   0.004   0.001   0.000   0.000   0.003   0.000   0.006   0.006   0.004   0.002   0.000   0.007   0.000   0.000
   0.038   0.000   0.315   0.109   0.006   0.041   0.027   0.009   0.033   0.004   0.009   0.088   0.051   0.089   0.023   0.065   0.048   0.013   0.012   0.011   0.009
   0.044   0.011   0.111   0.305   0.011   0.048   0.026   0.011   0.059   0.013   0.009   0.068   0.069   0.086   0.053   0.033   0.045   0.017   0.012   0.018   0.000
   0.017   0.000   0.005   0.007   0.415   0.004   0.009   0.039   0.025   0.097   0.042   0.013   0.006   0.011   0.009   0.009   0.014   0.041   0.053   0.085   0.009
   0.065   0.000   0.070   0.042   0.006   0.370   0.017   0.022   0.029   0.013   0.015   0.036   0.043   0.031   0.013   0.068   0.049   0.014   0.009   0.021   0.045
   0.010   0.000   0.012   0.011   0.010   0.007   0.571   0.003   0.022   0.005   0.015   0.043   0.006   0.035   0.021   0.016   0.008   0.017   0.009   0.037   0.009
   0.029   0.014   0.009   0.008   0.048   0.021   0.004   0.325   0.017   0.076   0.107   0.018   0.007   0.007   0.015   0.014   0.033   0.112   0.016   0.030   0.018
   0.053   0.007   0.044   0.081   0.020   0.041   0.044   0.026   0.336   0.029   0.059   0.073   0.045   0.094   0.163   0.041   0.054   0.026   0.041   0.028   0.036
   0.038   0.000   0.006   0.018   0.210   0.019   0.004   0.139   0.033   0.415   0.225   0.033   0.016   0.041   0.028   0.029   0.026   0.133   0.037   0.057   0.036
   0.013   0.000   0.004   0.003   0.016   0.007   0.000   0.043   0.014   0.053   0.197   0.010   0.000   0.018   0.004   0.003   0.010   0.018   0.021   0.021   0.018
   0.031   0.007   0.057   0.035   0.010   0.026   0.054   0.012   0.034   0.012   0.013   0.195   0.015   0.066   0.026   0.037   0.046   0.012   0.002   0.048   0.000
   0.022   0.000   0.036   0.035   0.005   0.026   0.011   0.009   0.020   0.006   0.000   0.013   0.424   0.013   0.016   0.039   0.011   0.009   0.002   0.000   0.000
   0.025   0.011   0.045   0.039   0.011   0.021   0.031   0.004   0.045   0.015   0.035   0.059   0.015   0.183   0.029   0.030   0.030   0.008   0.007   0.025   0.009
   0.019   0.011   0.012   0.023   0.005   0.008   0.019   0.010   0.069   0.009   0.004   0.018   0.013   0.028   0.348   0.030   0.019   0.005   0.007   0.018   0.018
   0.086   0.021   0.075   0.047   0.012   0.079   0.033   0.020   0.041   0.020   0.009   0.089   0.082   0.069   0.063   0.264   0.096   0.028   0.005   0.020   0.054
   0.043   0.007   0.039   0.033   0.020   0.038   0.014   0.026   0.032   0.015   0.026   0.057   0.028   0.046   0.035   0.065   0.266   0.037   0.016   0.034   0.000
   0.055   0.000   0.018   0.021   0.069   0.022   0.044   0.178   0.025   0.111   0.016   0.018   0.025   0.017   0.015   0.129   0.060   0.350   0.012   0.043   0.162
   0.009   0.000   0.003   0.004   0.022   0.004   0.007   0.006   0.012   0.006   0.020   0.001   0.001   0.006   0.004   0.002   0.007   0.003   0.588   0.064   0.000
   0.009   0.000   0.006   0.006   0.046   0.006   0.029   0.014   0.007   0.013   0.031   0.033   0.003   0.020   0.010   0.007   0.017   0.016   0.078   0.377   0.027
   0.009   0.000   0.001   0.000   0.001   0.004   0.001   0.002   0.002   0.002   0.004   0.000   0.000   0.004   0.003   0.006   0.004   0.010   0.000   0.005   0.297
   0.028   0.004   0.041   0.074   0.010   0.029   0.017   0.022   0.050   0.031   0.033   0.031   0.045   0.039   0.028   0.047   0.034   0.032   0.002   0.021   0.000
//
"""

# Run tests if called from the command line
if __name__ == '__main__':
    main()


                                                                                                                          

