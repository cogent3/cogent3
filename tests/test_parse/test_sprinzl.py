#!/usr/bin/env python
#file evo/parsers/test_sprinzl.py
"""Unit tests for the Sprinzl tRNA database parser.
"""
from string import strip
from cogent.parse.sprinzl import OneLineSprinzlParser, GenomicSprinzlParser,\
    _fix_sequence, get_pieces, get_counts, sprinzl_to_vienna
from cogent.util.unit_test import TestCase, main

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Jeremy Widmann", "Sandra Smit"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

sample_file = """Accession@@@AA@@@Anticodon@@@Species@@@Strain@@@0@@@1@@@2@@@3@@@4@@@5@@@6@@@7@@@8@@@9@@@10@@@11@@@12@@@13@@@14@@@15@@@16@@@17@@@17A@@@18@@@19@@@20@@@20A@@@20B@@@21@@@22@@@23@@@24@@@25@@@26@@@27@@@28@@@29@@@30@@@31@@@32@@@33@@@34@@@35@@@36@@@37@@@38@@@39@@@40@@@41@@@42@@@43@@@44@@@45@@@e11@@@e12@@@e13@@@e14@@@e15@@@e16@@@e17@@@e1@@@e2@@@e3@@@e4@@@e5@@@e27@@@e26@@@e25@@@e24@@@e23@@@e22@@@e21@@@46@@@47@@@48@@@49@@@5.@@@51@@@52@@@53@@@54@@@55@@@56@@@57@@@58@@@59@@@60@@@61@@@62@@@63@@@64@@@65@@@66@@@67@@@68@@@69@@@7.@@@71@@@72@@@73@@@74@@@75@@@76
GA0000001@@@Ala@@@TGC@@@Haemophilus influenzae@@@Rd KW20@@@-@@@G@@@G@@@G@@@G@@@C@@@C@@@T@@@T@@@A@@@G@@@C@@@T@@@C@@@A@@@G@@@C@@@T@@@-@@@G@@@G@@@G@@@-@@@-@@@A@@@G@@@A@@@G@@@C@@@G@@@C@@@C@@@T@@@G@@@C@@@T@@@T@@@T@@@G@@@C@@@A@@@C@@@G@@@C@@@A@@@G@@@G@@@A@@@G@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@G@@@T@@@C@@@A@@@G@@@C@@@G@@@G@@@T@@@T@@@C@@@G@@@A@@@T@@@C@@@C@@@C@@@G@@@C@@@T@@@A@@@G@@@G@@@C@@@T@@@C@@@C@@@A@@@-@@@-@@@-
GA0000002@@@Ala@@@GGC@@@Chlamydia pneumoniae @@@AR39@@@-@@@G@@@G@@@G@@@G@@@T@@@A@@@T@@@T@@@A@@@G@@@C@@@T@@@C@@@A@@@G@@@T@@@T@@@-@@@G@@@G@@@T@@@-@@@-@@@A@@@G@@@A@@@G@@@C@@@G@@@C@@@A@@@A@@@C@@@A@@@A@@@T@@@G@@@G@@@C@@@A@@@T@@@T@@@G@@@T@@@T@@@G@@@A@@@G@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@G@@@T@@@C@@@A@@@G@@@C@@@G@@@G@@@T@@@T@@@C@@@G@@@A@@@C@@@C@@@C@@@C@@@G@@@C@@@T@@@A@@@T@@@G@@@C@@@T@@@C@@@C@@@-@@@-@@@-@@@-
GA0000003@@@Ala@@@TGC@@@Chlamydia pneumoniae @@@AR39@@@-@@@G@@@G@@@G@@@G@@@A@@@C@@@T@@@T@@@A@@@G@@@C@@@T@@@T@@@A@@@G@@@T@@@T@@@-@@@G@@@G@@@T@@@-@@@-@@@A@@@G@@@A@@@G@@@C@@@G@@@T@@@C@@@T@@@G@@@A@@@T@@@T@@@T@@@G@@@C@@@A@@@T@@@T@@@C@@@A@@@G@@@A@@@A@@@G@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@-@@@G@@@T@@@C@@@A@@@G@@@G@@@A@@@G@@@T@@@T@@@C@@@G@@@A@@@A@@@T@@@C@@@T@@@C@@@C@@@T@@@A@@@G@@@T@@@C@@@T@@@C@@@C@@@-@@@-@@@-@@@-"""

sample_lines = ['\t'.join(i.split('@@@')) for i in sample_file.split('\n')]

class OneLineSprinzlParserTests(TestCase):
    """Tests of OneLineSprinzlParser"""
    def setUp(self):
        """standard tRNA file"""
        self.tRNAs = sample_lines #open('data_sprinzl.txt').read().split('\n')

    def test_minimal(self):
        """OneLineSprinzlParser should work on a minimal 'file'"""
        small = ['acc\taa\tac\tsp\tst\ta\tb\tc','q\tw\te\tr\tt\tA\tC\tG']
        p = OneLineSprinzlParser(small)
        result = list(p)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], 'ACG')
        self.assertEqual(result[0].Info.Accession, 'q')
        self.assertEqual(result[0].Info.AA, 'w')
        self.assertEqual(result[0].Info.Anticodon, 'e')
        self.assertEqual(result[0].Info.Species, 'r')
        self.assertEqual(result[0].Info.Strain, 't')
        
    def test_init(self):
        """OneLineSprinzlParser should read small file correctly"""
        p = OneLineSprinzlParser(self.tRNAs)
        recs = list(p)
        self.assertEqual(len(recs), 3)
        first, second, third = recs
        assert first.Info.Labels is second.Info.Labels
        assert first.Info.Labels is third.Info.Labels

        expected_label_list = "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 17A 18 19 20 20A 20B 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 e11 e12 e13 e14 e15 e16 e17 e1 e2 e3 e4 e5 e27 e26 e25 e24 e23 e22 e21 46 47 48 49 5. 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 7. 71 72 73 74 75 76".split()
        exp_labels = {}
        for i, label in enumerate(expected_label_list):
            exp_labels[label] = i

        self.assertEqual(first.Info.Labels, exp_labels)

        self.assertEqual(first.Info.Accession, 'GA0000001')
        self.assertEqual(first.Info.AA, 'Ala')
        self.assertEqual(first.Info.Anticodon, 'TGC')
        self.assertEqual(first.Info.Species, 'Haemophilus influenzae')
        self.assertEqual(first.Info.Strain, 'Rd KW20')
        self.assertEqual(first, '-GGGGCCTTAGCTCAGCT-GGG--AGAGCGCCTGCTTTGCACGCAGGAG-------------------GTCAGCGGTTCGATCCCGCTAGGCTCCA---'.replace('T','U'))

        self.assertEqual(third.Info.Accession, 'GA0000003')
        self.assertEqual(third.Info.AA, 'Ala')
        self.assertEqual(third.Info.Anticodon, 'TGC')
        self.assertEqual(third.Info.Species, 'Chlamydia pneumoniae')
        self.assertEqual(third.Info.Strain, 'AR39')
        self.assertEqual(third, '-GGGGACTTAGCTTAGTT-GGT--AGAGCGTCTGATTTGCATTCAGAAG-------------------GTCAGGAGTTCGAATCTCCTAGTCTCC----'.replace('T','U'))
         
genomic_sample = """3\t5950\tsequences\t\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t17A\t18\t19\t20\t20A\t20B\t21\t22\t23\t24\t25\t26\t27\t28\t29\t30\t31\t32\t33\t34\t35\t36\t37\t38\t39\t40\t41\t42\t43\t44\t45\te11\te12\te13\te14\te15\te16\te17\te1\te2\te3\te4\te5\te27\te26\te25\te24\te23\te22\te21\t46\t47\t48\t49\t50\t51\t52\t53\t54\t55\t56\t57\t58\t59\t60\t61\t62\t63\t64\t65\t66\t67\t68\t69\t70\t71\t72\t73\t74\t75\t76\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t
2\tGA0000001\tAla\t\tTGC\t\tHaemophilus influenzae\t\t\t\t\t\t\t\t\t\tRd KW20\t\t\t\tBacteria; Proteobacteria; gamma subdivision; Pasteurellaceae; Haemophilus\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t
\t\t\t\t-\tG\tG\tG\tG\tC\tC\tT\tT\tA\tG\tC\tT\tC\tA\tG\tC\tT\t-\tG\tG\tG\t-\t-\tA\tG\tA\tG\tC\tG\tC\tC\tT\tG\tC\tT\tT\tT\tG\tC\tA\tC\tG\tC\tA\tG\tG\tA\tG\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tG\tT\tC\tA\tG\tC\tG\tG\tT\tT\tC\tG\tA\tT\tC\tC\tC\tG\tC\tT\tA\tG\tG\tC\tT\tC\tC\tA\t-\t-\t-\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t
\t\t\t\t\t=\t=\t*\t=\t=\t=\t=\t\t\t=\t=\t=\t=\t\t\t\t\t\t\t\t\t\t\t\t=\t=\t=\t=\t\t=\t=\t=\t=\t=\t\t\t\t\t\t\t\t=\t=\t=\t=\t=\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t=\t=\t=\t=\t=\t\t\t\t\t\t\t\t=\t=\t=\t=\t=\t=\t=\t=\t=\t*\t=\t=\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t
3\tGA0000002\tAla\t\tGGC\t\tChlamydia pneumoniae \t\t\t\t\t\t\t\t\t\tAR39\t\t\t\tBacteria; Chlamydiales; Chlamydiaceae; Chlamydophila; Chlamydophila pneumoniae\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t
\t\t\t\t-\tG\tG\tG\tG\tT\tA\tT\tT\tA\tG\tC\tT\tC\tA\tG\tT\tT\t-\tG\tG\tT\t-\t-\tA\tG\tA\tG\tC\tG\tC\tA\tA\tC\tA\tA\tT\tG\tG\tC\tA\tT\tT\tG\tT\tT\tG\tA\tG\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tG\tT\tC\tA\tG\tC\tG\tG\tT\tT\tC\tG\tA\tC\tC\tC\tC\tG\tC\tT\tA\tT\tG\tC\tT\tC\tC\t-\t-\t-\t-\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t

\t\t\t\t\t=\t=\t*\t=\t*\t=\t=\t\t\t=\t=\t=\t=\t\t\t\t\t\t\t\t\t\t\t\t=\t=\t=\t=\t\t=\t=\t=\t=\t=\t\t\t\t\t\t\t\t=\t=\t=\t=\t=\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t=\t=\t=\t=\t=\t\t\t\t\t\t\t\t=\t=\t=\t=\t=\t=\t=\t*\t=\t*\t=\t=\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t

4\tGA0000003\tAla\t\tTGC\t\tChlamydia pneumoniae \t\t\t\t\t\t\t\t\t\tAR39\t\t\t\tBacteria; Chlamydiales; Chlamydiaceae; Chlamydophila; Chlamydophila pneumoniae\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t
\t\t\t\t-\tG\tG\tG\tG\tA\tC\tT\tT\tA\tG\tC\tT\tT\tA\tG\tT\tT\t-\tG\tG\tT\t-\t-\tA\tG\tA\tG\tC\tG\tT\tC\tT\tG\tA\tT\tT\tT\tG\tC\tA\tT\tT\tC\tA\tG\tA\tA\tG\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tG\tT\tC\tA\tG\tG\tA\tG\tT\tT\tC\tG\tA\tA\tT\tC\tT\tC\tC\tT\tA\tG\tT\tC\tT\tC\tC\t-\t-\t-\t-\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t
\t\t\t\t\t=\t=\t*\t=\t=\t=\t=\t\t\t=\t=\t=\t*\t\t\t\t\t\t\t\t\t\t\t\t*\t=\t=\t=\t\t=\t=\t=\t=\t=\t\t\t\t\t\t\t\t=\t=\t=\t=\t=\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t=\t=\t=\t=\t=\t\t\t\t\t\t\t\t=\t=\t=\t=\t=\t=\t=\t=\t=\t*\t=\t=\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t""".split('\n')

class GenomicSprinzlParserTests(TestCase):
    """Tests of the GenomicSprinzlParser class."""
    
    def test_single(self):
        """GenomicSprinzlParser should work with single sequence"""
        seqs = list(GenomicSprinzlParser(genomic_sample[0:4]))
        self.assertEqual(len(seqs), 1)
        s = seqs[0]
        self.assertEqual(s, '-GGGGCCTTAGCTCAGCT-GGG--AGAGCGCCTGCTTTGCACGCAGGAG-------------------GTCAGCGGTTCGATCCCGCTAGGCTCCA---'.replace('T','U'))
        self.assertEqual(s.Info.Accession, 'GA0000001')
        self.assertEqual(s.Info.AA, 'Ala')
        self.assertEqual(s.Info.Anticodon, 'UGC')
        self.assertEqual(s.Info.Species, 'Haemophilus influenzae')
        self.assertEqual(s.Info.Strain, 'Rd KW20')
        self.assertEqual(s.Info.Taxonomy, ['Bacteria', 'Proteobacteria', \
            'gamma subdivision', 'Pasteurellaceae', 'Haemophilus'])
        self.assertEqual(s.Pairing, '.==*====..====...........====.=====.......=====........................=====.......=========*==....')

    def test_multi(self):
        """GenomicSprinzlParser should work with multiple sequences"""
        seqs = list(GenomicSprinzlParser(genomic_sample))
        self.assertEqual(len(seqs), 3)
        self.assertEqual([s.Info.Accession for s in seqs], \
            ['GA0000001', 'GA0000002', 'GA0000003'])
        self.assertEqual(seqs[2].Info.Anticodon, 'UGC')
        self.assertEqual(seqs[0].Info.Order, seqs[2].Info.Order)

class FixSequenceTests(TestCase):
    """Tests that _fix_structure functions properly."""
    
    def test_fix_sequence(self):
        """Fix sequence should properly replace terminal gaps with CCA"""
        seqs = ['','ACGUUCC-','ACGUUC--','ACGUU---','ACGU----']
        results = ['','ACGUUCCA','ACGUUCCA','ACGUUCCA','ACGU-CCA']
        for s,r in zip(seqs,results):
            self.assertEqual(_fix_sequence(s),r)

class SprinzlToViennaTests(TestCase):

    def setUp(self):
        """setUp function for SprinzlToViennaTests"""
        self.structures = map(strip,STRUCTURES.split('\n'))
        self.vienna_structs = map(strip,VIENNA.split('\n'))
        self.short_struct = '...===...===.'
        #structure too long
        self.incorrect1 = ''.join(['..=====*..*==.............==*...=.=...',
            '....=.=..........................=====.......=====*=====......'])
        #two halves don't match 
        self.incorrect2 = ''.join(['..=====*..*===............==*...=.=...',
            '....=.=..........................=====.......=====*=====.....'])

    def test_get_pieces(self):
        """get_pieces: should return the correct pieces"""
        splits = [0,3,7,-1,13]
        self.assertEqual(get_pieces(self.short_struct, splits),\
            ['...','===.','..===','.'])
        #will include empty strings for indices outside of the structure
        self.assertEqual(get_pieces(self.short_struct,[2,10,20,30]),\
            ['.===...=','==.',''])
        #will return empty list if no break-positions are given
        self.assertEqual(get_pieces(self.short_struct,[]), [])

    def test_get_counts(self):
        """get_counts: should return list of lengths of paired regions"""
        self.assertEqual(get_counts('.===.=..'),[3,1])
        self.assertEqual(get_counts('...====..'),[4])
        self.assertEqual(get_counts('...'),[])

    def test_sprinzl_to_vienna(self):
        """sprinzl_to_vienna: should give expected output"""
        #This should only work for correct 
        for sprinzl,vienna in zip(self.structures, self.vienna_structs):
            self.assertEqual(sprinzl_to_vienna(sprinzl),vienna)
        #Check two obvious errors
        self.assertRaises(AssertionError,sprinzl_to_vienna,self.incorrect1)
        self.assertRaises(AssertionError,sprinzl_to_vienna,self.incorrect2)
        


STRUCTURES="""\
.===*===..====...........====.=====.......=====........................=====.......========*===....
.=======..=*=.............=*=.=====.......=====........................=====.......============....
.=======..====...........====.*====.......====*........................=.===.......===.========....
.=====.=..====...........====..===.........===.........................**===.......===**=.=====....
....====..====...........====.*====.......====*........................=====.......=========.......
..=====*..*==.............==*...=.=.......=.=..........................=====.......=====*=====.....
.=.=.==*..*.=*...........*=.*.*=.==.......==.=*........................=====.......=====*==.=.=....
.====*.=..**.=...........=.**.*====.......====*........................=.*.=.......=.*.==.*====....
.====*.=..**.=............=**.*====.......====*........................=.*.=.......=.*.==.*====...."""

VIENNA="""\
.(((((((..((((...........)))).(((((.......)))))........................(((((.......))))))))))))....
.(((((((..(((.............))).(((((.......)))))........................(((((.......))))))))))))....
.(((((((..((((...........)))).(((((.......)))))........................(.(((.......))).))))))))....
.(((((.(..((((...........))))..(((.........))).........................(((((.......)))))).)))))....
....((((..((((...........)))).(((((.......)))))........................(((((.......))))))))).......
..((((((..(((.............)))...(.(.......).)..........................(((((.......))))))))))).....
.(.(.(((..(.((...........)).).((.((.......)).))........................(((((.......)))))))).).)....
.(((((.(..((.(...........).)).(((((.......)))))........................(.(.(.......).).)).)))))....
.(((((.(..((.(............))).(((((.......)))))........................(.(.(.......).).)).)))))...."""


if __name__ == '__main__':
    main()
