#!/usr/bin/env python

from string import split, strip
from os import popen, remove
from glob import glob
from cogent.util.unit_test import TestCase, main
from cogent.parse.blast import QMEBlast9
from cogent.app.blast import seqs_to_stream,\
    make_subject_match_scorer, make_shotgun_scorer, keep_everything_scorer, \
    ids_from_seq_lower_threshold, PsiBlast, psiblast_n_neighbors

__author__ = "Micah Hamady"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Micah Hamady", "Rob Knight", "Catherine Lozupone"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Micah Hamady"
__email__ = "hamady@colorado.edu"
__status__ = "Prototype"

class BlastTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """Define some standard data"""
        self.rec = """# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 1
# Query: ece:Z4181
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4181	ece:Z4181	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4181	ecs:ECs3717	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4181	cvi:CV2421	41.67	72	42	0	39	110	29	100	2e-06	52.8
# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 2
# Query: ece:Z4181
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4181	ece:Z4181	100.00	110	0	0	1	110	1	110	3e-54	 211
ece:Z4181	ecs:ECs3717	100.00	110	0	0	1	110	1	110	3e-54	 211
ece:Z4181	cvi:CV2421	41.67	72	42	0	39	110	29	100	2e-08	59.0
ece:Z4181	sfl:CP0138	33.98	103	57	2	8	110	6	97	6e-06	50.5
ece:Z4181	spt:SPA2730	37.50	72	45	0	39	110	30	101	1e-05	49.8
ece:Z4181	sec:SC2804	37.50	72	45	0	39	110	30	101	1e-05	49.8
ece:Z4181	stm:STM2872	37.50	72	45	0	39	110	30	101	1e-05	49.8
# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 1
# Query: ece:Z4182
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4182	ece:Z4182	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4182	ecs:ECs3718	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4182	cvi:CV2422	41.67	72	42	0	39	110	29	100	2e-06	52.8""".split('\n')

        self.rec2 = """# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 1
# Query: ece:Z4181
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4181	ece:Z4181	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4181	ecs:ECs3717	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4181	spt:SPA2730	37.50	72	45	0	39	110	30	101	1e-05	49.8
# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 2
# Query: ece:Z4181
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4181	ecs:ECs3717	100.00	110	0	0	1	110	1	110	3e-54	 211
ece:Z4181	cvi:CV2421	41.67	72	42	0	39	110	29	100	2e-08	59.0
# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 1
# Query: ece:Z4182
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4182	ece:Z4182	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4182	cvi:CV2421	41.67	72	42	0	39	110	29	100	2e-06	52.8""".split('\n')

        self.rec3 = """# BLASTP 2.2.10 [Oct-19-2004]
# BLASTP 2.2.10 [Oct-19-2004]
# Query: ece:Z4181
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4181	ece:Z4182	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4181	ecs:ECs3717	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4181	spt:SPA2730	37.50	72	45	0	39	110	30	101	1e-05	49.8
# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 1
# Query: ece:Z4182
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4182	ece:Z4182	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4182	cvi:CV2421	41.67	72	42	0	39	110	29	100	2e-06	52.8
# BLASTP 2.2.10 [Oct-19-2004]
# Query: ece:Z4183
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4183	ece:Z4182	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4183	ecs:ECs3717	100.00	110	0	0	1	110	1	110	3e-54	 211
ece:Z4183	cvi:CV2421	41.67	72	42	0	39	110	29	100	2e-08	59.0""".split('\n')

        self.query_1 = """>gi|100002553| Bd2556c Bd2556c two-component system sensor histidine kinase 3092017:3094158 reverse MW:81963
MRLKNRLNNWISIRMGMVIVIFLGVSCGSMRSSTPPPAKDRLTEIDSLERLLPDCPTIASTLPLLRRLAFLYQQQSEMKVYNERLYENAMAVDSISVAYLGLKNLAEYYYDQSVRDSLEYYCSLVDSIAKARHEYPNVLFDVKSLSSQDLLWLGNYELAMSEAMDLYRLASNLDHRYGLLRCSETLGLIYQRIRRDSDAVVSFQESLDLLKDIKDVPDIMDTKVRLTSYQLESSVRTKQYASTERILGQYMALLDEQYKIYQEKNDLLSIKREYWLLYSFYTSFYLSQGDLENAKRSLDQASSYADSNWVEGDYAINTYLTVKARYHKAAGDIPLALHCINEVLETERLPEDIQFKADILKEQGQLGEVMALYDELYSTLTKRRGTSFLRQVNQLRTLHELHEKELKETELKEAGQRIARKQDLLIFILSISVVLLILLYVLFLYYRHLRSLKNQLQREKELLLESQRQLIKEKTRAEEASLMKSAFLANMSHEVRTPLNAIVGFSGLLVEPSTDEEERKEYSSIIRNNTDLMLNLVNDVLDLSRMETGDLHFDIKDHLLLVCCQMALESVRHRIPDGVKLTFSPAGEPIVVHVDNLRLQQLLTNLLTNAAKFTEKGEINLSFQLEPDRKKVCIAVTDTGAGIPLEKQATIFNRFEKLDDYKPGVGLGLSICLLIAERLDGALFIDSSYTDGARFVLILSCEIDSSIYNPPIEV"""

        self.query_2 = """>gi|100002557| Bd2560c Bd2560c conserved hypothetical protein 3097971:3098210 reverse MW:8927
MGKNQLIHGNEFHLLKQAEIHKATGKLVESLNLAAGSTGGFDIYKVVEAYFTDLEKRKEINDLLGISEPCETRVTEECFS
"""

        self.fasta_recs = """>gi|100002550| Bd2553c Bd2553c conserved hypothetical protein 3090609:3091013 reverse MW:14682
MMDFISVPLVVGIVCAGIYGLFELFVRKRERLAIIEKIGDKLDTSAFDGKLGLPNYMRNFSFSSLKAGCLLAGIGLGLLVGFIINMCMATNSYYDDGWYRHEVAGTAYGASVLLFGGIGLIIAFVIELKLGKNNK
>gi|100002551| Bd2554 Bd2554 RNA polymerase ECF-type sigma factor 3091112:3091717 forward MW:23408
LLPQVVTYLPGLRPLSTMELYTDTYYIQRIQAGDVACFACLLDKYSRPIHSLILKVVRSQEEAEELAQDTFMKVFKNLASFKGDCSFSTWIYRIAYNTAISSVRKKRYEFLAIEETTLENVSEEEITNLFGQTESTEQVQRLEVALEQLLPDERALILLFYWKEKTIEELVSITGLTASNIKVKLHRIRKKLFVLLNGMDHE
>gi|100002552| Bd2555 Bd2555 conserved hypothetical protein 3091713:3092066 forward MW:13332
MSKINTNKEQPDLLGDLFKRIPEEELPASFRSNVMRQIMLESAKAKKRDERFSLLAAIVASLIMISLAIVSFVYMEIPKIAIPTISTSALAFYLYIGAITLILLLADYKLRNLFHKKG
>gi|100002553| Bd2556c Bd2556c two-component system sensor histidine kinase 3092017:3094158 reverse MW:81963
MRLKNRLNNWISIRMGMVIVIFLGVSCGSMRSSTPPPAKDRLTEIDSLERLLPDCPTIASTLPLLRRLAFLYQQQSEMKVYNERLYENAMAVDSISVAYLGLKNLAEYYYDQSVRDSLEYYCSLVDSIAKARHEYPNVLFDVKSLSSQDLLWLGNYELAMSEAMDLYRLASNLDHRYGLLRCSETLGLIYQRIRRDSDAVVSFQESLDLLKDIKDVPDIMDTKVRLTSYQLESSVRTKQYASTERILGQYMALLDEQYKIYQEKNDLLSIKREYWLLYSFYTSFYLSQGDLENAKRSLDQASSYADSNWVEGDYAINTYLTVKARYHKAAGDIPLALHCINEVLETERLPEDIQFKADILKEQGQLGEVMALYDELYSTLTKRRGTSFLRQVNQLRTLHELHEKELKETELKEAGQRIARKQDLLIFILSISVVLLILLYVLFLYYRHLRSLKNQLQREKELLLESQRQLIKEKTRAEEASLMKSAFLANMSHEVRTPLNAIVGFSGLLVEPSTDEEERKEYSSIIRNNTDLMLNLVNDVLDLSRMETGDLHFDIKDHLLLVCCQMALESVRHRIPDGVKLTFSPAGEPIVVHVDNLRLQQLLTNLLTNAAKFTEKGEINLSFQLEPDRKKVCIAVTDTGAGIPLEKQATIFNRFEKLDDYKPGVGLGLSICLLIAERLDGALFIDSSYTDGARFVLILSCEIDSSIYNPPIEV
>gi|100002554| Bd2557c Bd2557c two-component system sensor histidine kinase 3094158:3095507 reverse MW:51247
LERKYNGEGKIFPVKRHRCLMSCYYCELYTMKGNSGKAQAYLDQATAYLDSSFGDRVEAQYLRTKSFYYWKEKDYRHALSAVNLALKINRDLDKLEMKKAVLQSSGQLQEAVTIYEEIINKTETINTDAFDRQIEQLRVLNDLNDLEKQDRELKLKSEQEALKQKQIVVSIGLLLVLMGLLYMLWRIYMHTKRLRNELLQEKDSLTASEKQLRVVTKEAEAANKKKSAFIANISHEVRTPLNAIVGFSELLASSEYSEEEKIRFAGEVNHSSELLLNLVNDVLDLSRLESGKIKFSVKPNDLVACCQRALDSIRHRVKPGVRLTFTPSIESYTLNTDALRLQQLLTNLLSNAAKFTSEGEINLSFTVDEGKEEVCFSVTDTGCGIPEDKCEKIFERFEKLDDFIQGTGLGLSVCQIISEQLNGSLSVDISYKDGARFVFIHPTNLIETPI
>gi|100002555| Bd2558c Bd2558c hypothetical protein 3095527:3095985 reverse MW:17134
LRGKNIHLGRVGCNYGKLLIFIDIYFVSLRIVSDKSMSRGFLRKSSVNTFIGIVWILFAVGTSAQNAVSKFRADSIRQSLSRIQKPQDKIPLLKELIGLYWQLPEEVLALKEIIDIAMPLDSIGIVYDAMAGLSRYYPAIRTFVRVGGALETV
>gi|100002556| Bd2559 Bd2559 30S ribosomal protein S1 3096095:3097882 forward MW:67092
MENLKNIQPVEDFNWDAFEQGETYTEVSKDDLVKTYDETLNTVKDKEVVMGTVTSMNKREVVVNIGFKSDGVVPMSEFRYNPDLKIGDEVEVYIESQEDKKGQLILSHKKARATRSWDRVNEALEKDEIIKGYIKCRTKGGMIVDVFGIEAFLPGSQIDVKPIRDYDVFVGKTMEFKIVKINQEFKNVVVSHKALIEAELEQQKKDIISKLEKGQVLEGTVKNITSYGVFIDLGGVDGLIHITDLSWGRVSHPEEIVQLDQKINVVILDFDDEKKRIALGLKQLTPHPWDALDTNLKVGDKVKGKVVVMADYGAFIEIAPGVEGLIHVSEMSWTQHLRSAQDFMKVGDEIEAVILTLDRDERKMSLGIKQLKADPWENIEERFPVGSRHAAKVRNFTNFGVFVEIEEGVDGLIHISDLSWTKKIKHPSEFTQIGAEIEVQVLEIDKENRRLSLGHKQLEENPWDVFETIFTVGSIHEGTIIEVLDKGAVISLPYGVEGFATPKHLVKEDGSQAQVDEKLSFKVIEFNKEAKRIILSHSRIFEDEQKGAKATSEKKASSKRGGKKEEESGMVTGPVEKTTLGDIEELAALKEKLSGK
>gi|100002557| Bd2560c Bd2560c conserved hypothetical protein 3097971:3098210 reverse MW:8927
MGKNQLIHGNEFHLLKQAEIHKATGKLVESLNLAAGSTGGFDIYKVVEAYFTDLEKRKEINDLLGISEPCETRVTEECFS
>gi|100002558| Bd2561 Bd2561 phosphoglycolate phosphatase 3098389:3099033 forward MW:24182
MKKLVIFDLDGTLLNTIADLAHSTNHALRQNGFPTHDVKEYNFFVGNGINKLFERALPEGEKTAENILKVREEFLKHYDLHNTDRSVPYPGVPELLALLQERGIKLAVASNKYQAATRKLIAHFFPSIQFTEVLGQREGVKAKPDPSIVNEIVERASISKESTLYVGDSDVDMQTAINSEVTSCGVTWGFRPRTELEKYAPDHIAEKAEDILKFI
>gi|100002559| Bd2562 Bd2562 conserved hypothetical protein 3099382:3100299 forward MW:35872
MSGNIKKIVEPNSGIDYSLEKDFKIFTLSKELPITTYPSYIRLGIVIYCVKGNAKIDIYSNKHIITPKELIIILPGQLVALTDVSVDFQIRYFTITESFYSDILSGISRFSPHFFFYMRQHYYFKMEDVETLSFVDFFELLIRKAVDPENQYRRESVILLLRILFLDIYNHYKVNSLDSTATIDVHKKELTHKFFQLVMSNYKVNRSVTFYANSLCITPKYLTMVVKEVSGKSAKDWITEYMILELKGLLTNSTLNIQEIVEKTQFSNQSSLGRFFRRHTGLSPLQYRKKYLTTEQRTNFSKNNTI
"""

    def test_seqs_to_stream(self):
        """seqs_to_stream should iterate over seqs"""
        sts = seqs_to_stream
        self.assertEqual(list(sts('>a\nTG\n>b\nWW\n', \
            '_input_as_multiline_string')),\
            [['>a','TG'],['>b','WW']])
        #skipping test for file open
        self.assertEqual(list(sts(['TG','WW'], '_input_as_seqs')), \
            [['>0','TG'],['>1','WW']])
        self.assertEqual(list(sts(['>a','TG','>b','WW'], \
            '_input_as_lines')),\
            [['>a','TG'],['>b','WW']])
        self.assertRaises(TypeError, sts, 'abc', 'xyz')
   
    def test_make_subject_match_scorer(self):
        """make_subject_match_scorer should keep ids matching n queries"""
        qm1 = make_subject_match_scorer(1)
        qm3 = make_subject_match_scorer(3)
        qm5 = make_subject_match_scorer(5)
        qmes = wrap_qmes(QMEBlast9(self.rec3))
        self.assertEqualItems(qm1(qmes), ['ece:Z4181','ece:Z4182','ece:Z4183'])
        self.assertEqualItems(qm3(qmes), ['ece:Z4181','ece:Z4183'])
        self.assertEqualItems(qm5(qmes), [])
   
    def test_make_shotgun_scorer(self):
        """make_shotgun_scorer should keep ids matching n queries"""
        sg1 = make_shotgun_scorer(1)
        sg2 = make_shotgun_scorer(2)
        sg3 = make_shotgun_scorer(3)
        sg4 = make_shotgun_scorer(4)
        sg5 = make_shotgun_scorer(5)
        qmes = wrap_qmes(QMEBlast9(self.rec3))
        self.assertEqualItems(sg1(qmes), keep_everything_scorer(qmes))
        self.assertEqualItems(sg2(qmes), \
            ['ece:Z4181','ece:Z4182','ece:Z4183','cvi:CV2421','ecs:ECs3717'])
        self.assertEqualItems(sg3(qmes), \
            ['ece:Z4181','ece:Z4182','ece:Z4183'])
        self.assertEqualItems(sg4(qmes), \
            ['ece:Z4182'])
        self.assertEqualItems(sg5(qmes), [])
     
    def test_keep_everything_scorer(self):
        """keep_everything_scorer should keep all ids found."""
        k = keep_everything_scorer(wrap_qmes(QMEBlast9(self.rec2)))
        self.assertEqualItems(k, \
            ['ece:Z4181','ecs:ECs3717','spt:SPA2730','cvi:CV2421','ece:Z4182'])

    def test_ids_from_seq_lower_threshold(self):
        "ids_from_seq_lower_threshold returns psiblast hits, decreasing sens"
        bdb_seqs = self.fasta_recs
        f = open('test_bdb', 'w')
        f.write(bdb_seqs)
        f.close()
        temp = popen('formatdb -i test_bdb -o T -p T')
        params = {'-j':2,
                '-d':'test_bdb'}
        query = self.query_1.split('\n')
        app = PsiBlast(params=params, 
                InputHandler='_input_as_lines')
        #the command below should result in finding itself and 2554
        #it should run for max_iterations
        result = ids_from_seq_lower_threshold(query, n=12, \
                max_iterations=4, app=app, core_threshold=1e-50, \
                lower_threshold=1e-20, step=10000)
        self.assertEqual(result[0],\
                [('gi|100002553', '0.0'), ('gi|100002554', '0.0')])
        self.assertEqual(result[1], 4)
        #if n=2, it should find the same sequences but only run for 1 iteration
        #since it would hit n after the first blast search
        result = ids_from_seq_lower_threshold(query, n=2, \
                max_iterations=4, app=app, core_threshold=1e-50, \
                lower_threshold=1e-20, step=10000)
        self.assertEqual(result[0],\
                [('gi|100002553', '0.0'), ('gi|100002554', '0.0')])
        self.assertEqual(result[1], 1)
        query = self.query_2.split('\n')
        #query_2_s e-value for itself is 9e-47, it should not be found
        #with the lower_threshold set to 1e-48
        result = ids_from_seq_lower_threshold(query, n=12, \
                max_iterations=4, app=app, core_threshold=1e-50, \
                lower_threshold=1e-48, step=10000)
        self.assertEqual(result[0], [])
        #it also should not be found if the max_iterations is set to 1
        result = ids_from_seq_lower_threshold(query, n=12, \
                max_iterations=1, app=app, core_threshold=1e-50, \
                lower_threshold=1e-20, step=10000)
        self.assertEqual(result[0], [])
        for fname in ['formatdb.log'] + glob('test_bdb*'):
            remove(fname)

    def test_psiblast_n_neighbors(self):
        "psiblast_n_neighbors psiblasts and stops when n neighbors are reached"
        bdb_seqs = self.fasta_recs
        f = open('test_bdb', 'w')
        f.write(bdb_seqs)
        f.close()
        temp = popen('formatdb -i test_bdb -o T -p T')
        params = {'-j':11}
        lines = bdb_seqs.split('\n')
        results = psiblast_n_neighbors(lines, n=12, blast_db='test_bdb', \
                method='lower_threshold', params=params,\
                core_threshold=1e-50, step=10000)
        #there should be 10 result entries since there were 10 queries
        self.assertEqual(len(results), 10)
        for i in results:
            #each query should at least find itself
            self.failUnless(len(results[i][0]) >= 1)
            #each query should iterate 8 times since it can never reach max
            self.assertEqual(results[i][1], 11)
        for fname in ['formatdb.log'] + glob('test_bdb*'):
            remove(fname)
        

def wrap_qmes(qmes):
    """Converts qmes into a dict of {q:{m:e}}"""
    d = {}
    for q, m, e in qmes:
        if q not in d:
            d[q] = {}
        d[q][m] = e
    return d

if __name__ == "__main__":
    main()
