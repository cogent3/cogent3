#!/usr/bin/env python
from __future__ import division

import unittest, os, tempfile

from cogent import DNA, RNA, STANDARD_CODON as CODON, PROTEIN, Sequence, \
                LoadSeqs

__author__ = "Peter Maxwell, Gavin Huttley and Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

base_path = os.getcwd()
data_path = os.path.join(base_path, 'data')


class ReadingWritingFileFormats(unittest.TestCase):
    """Testing ability to read file formats."""
    def setUp(self):
        pass
        
    def _loadfromfile(self, filename, test_write=True, **kw):
        filename = os.path.join(data_path, filename)
        aln = LoadSeqs(filename=filename, **kw)
        if test_write:
            suffix = filename.split('.')[-1]
            fn = tempfile.mktemp(suffix='.'+suffix)
            aln.writeToFile(filename=fn)
            os.remove(fn)
        
    def test_fasta(self):
        self._loadfromfile("formattest.fasta")
        
    def test_phylipsequential(self):
        self._loadfromfile('formattest.phylip')
        
    def test_clustal(self):
        self._loadfromfile('formattest.aln', test_write=False)
        
    def test_phylip_interleaved(self):
        self._loadfromfile('interleaved.phylip', test_write=False, interleaved=True)
        
    def test_paml(self):
        self._loadfromfile('formattest.paml')
        
    def test_gde(self):
        self._loadfromfile('formattest.gde')

    def test_msf(self):
        self._loadfromfile('formattest.msf', test_write=False)
        
class AlignmentTestMethods(unittest.TestCase):
    """Testing Alignment methods"""
    def setUp(self):
        self.alignment = LoadSeqs(filename = os.path.join(data_path,'brca1_5.paml'))

    def test_picklability(self):
        """Pickle an alignment containing an annotated sequence"""
        # This depends on alignments, sequences, features, maps and spans
        # Doesn't test round trip result is correct, which should possibly
        # be done for maps/spans, but seqs/alignments are just simple
        # python classes without __getstate__ etc.
        import cPickle as pickle
        seq1 = DNA.makeSequence("aagaagaagaccccca")
        seq2 = DNA.makeSequence("aagaagaagaccccct")
        seq2.addFeature('exon', 'fred', [(10,15)])
        aln = LoadSeqs(data={'a':seq1, 'b':seq2})
        aln2 = pickle.loads(pickle.dumps(aln))

    def test_empty_seq(self):
        """test creation of an alignment from scratch, with one sequence pure gap"""
        new_seqs = {'seq1': 'ATGATG', 'seq2': '------'}
        align = LoadSeqs(moltype=DNA, data=new_seqs)
        assert len(align) == 6, align
        
    def test_getNumSeqs(self):
        self.assertEqual(self.alignment.getNumSeqs(), 5)
        
    def test_variablePositions(self):
        new_seqs = {'seq1': 'ACGTACGT', 'seq2': 'ACCGACGT', 'seq3': 'ACGTACGT'}
        align = LoadSeqs(data=new_seqs)
        self.assertEqual(align.variablePositions(), [2,3])
            
    def test_numberseqs(self):
        """testing the number of sequences"""
        assert len(self.alignment.getSeqNames()) == 5
        
    def test_alignlsength(self):
        """testing the alignment length"""
        assert len(self.alignment) == 60
            
    def test_init_from_strings(self):
        """testing constructing an alignment from a dictionary of strings"""
        new_seqs = {'seq1': 'ACGTACGT', 'seq2': 'ACGTACGT', 'seq3': 'ACGTACGT'}
        LoadSeqs(data=new_seqs)
        
    def test_getSubAlignment(self):
        """test slicing otus, and return of new alignment"""
        fullset = ['DogFaced', 'Human', 'HowlerMon', 'Mouse', 'NineBande']
        subset = ['DogFaced', 'Human', 'HowlerMon', 'Mouse']
        subset.sort()
        sub_align = self.alignment.takeSeqs(subset)
        new = sub_align.getSeqNames()
        new.sort()
        assert new == subset, "included subset didn't work %s, %s" % (new, subset)
        
        # testing exclusion of one
        to_exclude = ['NineBande']
        sub_align = self.alignment.takeSeqs(to_exclude, negate = True)
        new = sub_align.getSeqNames()
        new.sort()
        assert new == subset, "excluded subset didn't work %s, %s" % (new, subset)
        
        # testing exclusion of two
        subset = ['DogFaced', 'HowlerMon', 'NineBande']
        subset.sort()
        to_exclude = ['Human', 'Mouse']
        sub_align = self.alignment.takeSeqs(to_exclude, negate = True)
        new = sub_align.getSeqNames()
        new.sort()
        assert new == subset, "excluded subset didn't work %s, %s" % (new, subset)
        
    def test_slice_align(self):
        """test slicing of sequences"""
        alignment = LoadSeqs(data={'seq1': 'ACGTACGT', 'seq2': 'ACGTACGT', 'seq3': 'ACGTACGT'})
        sub_align = alignment[2: 5]          
        self.assertEqual(len(sub_align), 3)  
        self.assertEqual(len(sub_align.getSeqNames()), 3)
        self.assertEqual(sub_align.todict(), {'seq1': 'GTA', 'seq2': 'GTA', 'seq3': 'GTA'})
                                             
        sub_align = alignment[5: 20]         
        self.assertEqual(len(sub_align), 3)  
        self.assertEqual(len(sub_align.getSeqNames()), 3)
        self.assertEqual(sub_align.todict(), {'seq1': 'CGT', 'seq2': 'CGT', 'seq3': 'CGT'})
                                             
        sub_align = alignment[2]             
        self.assertEqual(len(sub_align), 1)  
        self.assertEqual(sub_align.todict(), {'seq1': 'G', 'seq2': 'G', 'seq3': 'G'})
                                             
        sub_align = alignment[0]             
        self.assertEqual(len(sub_align), 1)  
        self.assertEqual(sub_align.todict(), {'seq1': 'A', 'seq2': 'A', 'seq3': 'A'})
                                             
        sub_align = alignment[7]             
        self.assertEqual(len(sub_align), 1)  
        self.assertEqual(sub_align.todict(), {'seq1': 'T', 'seq2': 'T', 'seq3': 'T'})
                                             
    def test_slidingWindows(self):          
        """test slicing of sequences"""      
        alignment = LoadSeqs(data = {'seq1': 'ACGTACGT', 'seq2': 'ACGTACGT', 'seq3': 'ACGTACGT'})
        result = []                          
        for bit in alignment.slidingWindows(5,2):
            result+=[bit]                    
        self.assertEqual(result[0].todict(), {'seq3': 'ACGTA', 'seq2': 'ACGTA', 'seq1': 'ACGTA'})
        self.assertEqual(result[1].todict(), {'seq3': 'GTACG', 'seq2': 'GTACG', 'seq1': 'GTACG'})
                                             
        result = []                          
        for bit in alignment.slidingWindows(5,1):
            result+=[bit]                    
        self.assertEqual(result[0].todict(), {'seq3': 'ACGTA', 'seq2': 'ACGTA', 'seq1': 'ACGTA'})
        self.assertEqual(result[1].todict(), {'seq3': 'CGTAC', 'seq2': 'CGTAC', 'seq1': 'CGTAC'})
        self.assertEqual(result[2].todict(), {'seq3': 'GTACG', 'seq2': 'GTACG', 'seq1': 'GTACG'})
        self.assertEqual(result[3].todict(), {'seq3': 'TACGT', 'seq2': 'TACGT', 'seq1': 'TACGT'})
        
    def test_withoutRedundantGaps(self):
        """test removal of redundant gaps (all entries in alignment column are gaps)"""
        alignment = LoadSeqs(data={'seq1': '--ACGT--GT---', 'seq2': '--ACGTA-GT---', 'seq3': '--ACGTA-GT---'})
        align_dict = alignment.omitGapPositions().todict()
        self.assertEqual(align_dict, {'seq1':'ACGT-GT', 'seq2':'ACGTAGT', 'seq3':'ACGTAGT'})
            
    def test_withoutAnyGaps(self):
        """test removal of all gaps (any entries in alignment column are gaps)"""
        alignment = LoadSeqs(data={'seq1': '--ACGT--GT---', 'seq2': '--ACGTA-GT---', 'seq3': '--ACGTA-GT---'})
        align_dict = alignment.omitGapPositions(allowed_gap_frac=0).todict()
        self.assertEqual(align_dict, {'seq1':'ACGTGT', 'seq2':'ACGTGT', 'seq3':'ACGTGT'})
            
        alignment = LoadSeqs(data={'seq1': 'ACGT', 'seq2': '----', 'seq3': '----'})
        align_dict = alignment.omitGapPositions(allowed_gap_frac=0).todict()
        self.assertEqual(align_dict, {'seq1':'', 'seq2':'', 'seq3':''})
    
    def test_degap(self):
        """test stripping gaps from collections and alignments"""
        aln = LoadSeqs(data={'seq1': '--ACGT--GT---', 'seq2': '--ACGTA-GT---',
                    'seq3': '--ACGTA-GT---'})
        observed = aln.degap()
        expect = {'seq1': 'ACGTGT', 'seq2': 'ACGTAGT', 'seq3': 'ACGTAGT'}
        self.assertEqual(observed.todict(), expect)
        collection = LoadSeqs(data={'seq1': '--ACGT--GT---',
                    'seq2': '--ACGTA-GT---', 'seq3': '--ACGTA-GT---'},
                    aligned=False, moltype=DNA)
        observed = collection.degap()
        self.assertEqual(observed.todict(), expect)
        self.assertEqual(observed.MolType, DNA)
    
    def test_DnaRna_interconversion(self):
        """test interconversion between Rna and Dna by SequenceCollection and
        Alignment"""
        dna = {'seq1': '--ACGT--GT---', 'seq2': '--ACGTA-GT---',
                    'seq3': '--ACGTA-GT---'}
        rna = {'seq1': '--ACGU--GU---', 'seq2': '--ACGUA-GU---',
                    'seq3': '--ACGUA-GU---'}
        aln_Dna = LoadSeqs(data=dna, moltype=DNA)
        aln_Rna = LoadSeqs(data=dna, moltype=RNA)
        collect_Dna = LoadSeqs(data=dna, aligned=False, moltype=DNA)
        collect_Rna = LoadSeqs(data=rna, aligned=False, moltype=RNA)
        assert aln_Rna.toDna().todict() == dna, (aln_Rna.toDna().todict(), dna)
        assert aln_Dna.toRna().todict() == rna, (aln_Dna.toRna().todict(), rna)
        assert collect_Rna.toDna().todict() == dna, \
                                        (collect_Rna.toDna().todict(), dna)
        assert collect_Dna.toRna().todict() == rna, \
                                        (collect_Dna.toRna().todict(), rna)
    
    def test_reversecomplement(self):
        """test reverse complementing of Alignments and SequenceCollection."""
        dna = {'seq1': '--ACGT--GT---', 'seq2': 'TTACGTA-GT---',
                    'seq3': '--ACGTA-GCC--'}
        dna_rc = {'seq1': '---AC--ACGT--', 'seq2': '---AC-TACGTAA',
                    'seq3': '--GGC-TACGT--'}
        # alignment with gaps
        aln = LoadSeqs(data=dna, moltype=DNA)
        aln_rc = aln.rc()
        self.assertEqual(aln_rc.todict(), dna_rc)
        # check collection, with gaps
        coll = LoadSeqs(data=dna, moltype=DNA, aligned=False)
        coll_rc = coll.rc()
        self.assertEqual(coll_rc.todict(), dna_rc)
        self.assertEqual(coll_rc.todict(), coll.reversecomplement().todict())
        # collection with no gaps
        dna = {'seq1': 'ACGTGT', 'seq2': 'TTACGTAGT',
                    'seq3': 'ACGTAGCC'}
        dna_rc = {'seq1': 'ACACGT', 'seq2': 'ACTACGTAA',
                    'seq3': 'GGCTACGT'}
        coll = LoadSeqs(data=dna, moltype=DNA, aligned=False)
        coll_rc = coll.rc()
        self.assertEqual(coll_rc.todict(), dna_rc)
    
    def test_getasdict(self):
        """getting the alignment as a dictionary"""
        seqs={'seq1': 'ACGT--GT', 'seq2': 'ACGTACGT', 'seq3': 'ACGTACGT'}
        alignment = LoadSeqs(data=seqs)
        align_dict = alignment.todict()
        self.assertEqual(align_dict, seqs)
        
    def test_alignadd(self):
        """testing adding one alignment to another."""
        align1= LoadSeqs(data={'a': 'AAAA', 'b': 'TTTT', 'c': 'CCCC'})
        align2 = LoadSeqs(data={'a': 'GGGG', 'b': '----', 'c': 'NNNN'})
        align = align1 + align2
        concatdict = align.todict()
        self.assertEqual(concatdict, {'a': 'AAAAGGGG', 'b': 'TTTT----', 'c': 'CCCCNNNN'})

    def test_replaceSeqs(self):
        """synchronize gaps between protein seqs and codon seqs"""
        pd={'FlyingFox': 'C-TNAH',
            'DogFaced':  'CGTNT-',
            'FreeTaile': '-GTDTH',
            'LittleBro': 'C-TD-H',
            'TombBat':   'C--STH'}
        pal = LoadSeqs(moltype = PROTEIN, data = pd)

        cu={'TombBat':   'TGTAGTACTCAT',
            'FreeTaile': 'GGCACAGATACTCAT',
            'FlyingFox': 'TGTACAAATGCTCAT',
            'LittleBro': 'TGTACAGATCAT',
            'DogFaced':  'TGTGGCACAAATACT'}

        co = LoadSeqs(moltype = DNA, data = cu, aligned = False)
        cal = pal.replaceSeqs(co)
        result = cal.todict()
        for taxon, expected_sequence in [
                ('FlyingFox', 'TGT---ACAAATGCTCAT'),
                ('DogFaced',   'TGTGGCACAAATACT---'),
                ('FreeTaile', '---GGCACAGATACTCAT'),
                ('LittleBro', 'TGT---ACAGAT---CAT'),
                ('TombBat',  'TGT------AGTACTCAT')]:
            self.assertEqual(result[taxon], expected_sequence)
            
    def test_sample(self):
        """Test sample generation"""
        alignment = LoadSeqs(data={'seq1': 'ABCDEFGHIJKLMNOP',
                                    'seq2': 'ABCDEFGHIJKLMNOP'})
        # effectively permute columns, preserving length
        shuffled = alignment.sample()
        # ensure length correct
        sample = alignment.sample(10)
        self.assertEqual(len(sample), 10)
        # test columns alignment preserved
        seqs = sample.todict().values()
        self.assertEqual(seqs[0], seqs[1])
        # ensure each char occurs once as sampling without replacement
        for char in seqs[0]:
            self.assertEqual(seqs[0].count(char), 1)
                
    def test_sample_with_replacement(self):
        #test with replacement
        alignment = LoadSeqs(data={'seq1': 'gatc', 'seq2': 'gatc'})
        sample = alignment.sample(1000, with_replacement=True)

    def test_sample_tuples(self):
        ##### test with motif size != 1 #####
        alignment = LoadSeqs(data={'seq1': 'AABBCCDDEEFFGGHHIIJJKKLLMMNNOOPP',
                                    'seq2': 'AABBCCDDEEFFGGHHIIJJKKLLMMNNOOPP'})
        shuffled = alignment.sample(motif_length=2)
        # ensure length correct
        sample = alignment.sample(10,motif_length=2)
        self.assertEqual(len(sample), 20)
        # test columns alignment preserved
        seqs = sample.todict().values()
        self.assertEqual(seqs[0], seqs[1])
        # ensure each char occurs twice as sampling dinucs without replacement
        for char in seqs[0]:
            self.assertEqual(seqs[0].count(char), 2)

    def test_translate(self):
        for seqs in [
                {'seq1': 'GATTTT', 'seq2': 'GATC??'}, 
                {'seq1': 'GAT---', 'seq2': '?GATCT'}]:
            alignment = LoadSeqs(data=seqs, moltype=DNA)
            self.assertEqual(len(alignment.getTranslation()), 2)
            # check for a failure when no moltype specified
            alignment = LoadSeqs(data=seqs)
            try:
                peps = alignment.getTranslation()
            except AttributeError:
                pass
        
        
    def test_seqnames(self):
        s1 = self.alignment.getSeq('Mouse')
        self.assertEqual(s1.getName(), 'Mouse')
    
    def test_withoutTerminalStopCodons(self):
        """test without terminal stop handling"""
        seq_coll = LoadSeqs(data = {'seq1': 'ACGTAA', 'seq2': 'ACGACG',
                            'seq3': 'ACGCGT'}, moltype = DNA, aligned=False)
        seq_coll = seq_coll.withoutTerminalStopCodons()
        seqs = seq_coll.todict()
        self.assertEqual(seqs['seq1'], 'ACG')   # note: not 'acg---'
        self.assertEqual(seqs['seq2'], 'ACGACG')
        aln = LoadSeqs(data = {'seq1': 'ACGTAA', 'seq2': 'ACGTGA',
                        'seq3': 'ACGTAA'}, moltype = DNA)
        aln = aln.withoutTerminalStopCodons()
        seqs = aln.todict()
        self.assertEqual(seqs['seq1'], 'ACG')   # note: not 'acg---'
        self.assertEqual(seqs['seq2'], 'ACG')
        self.assertEqual(seqs['seq3'], 'ACG')
    
    def test_hasTerminalStops(self):
        """test truth values for terminal stops"""
        # seq collections
        seq_coll = LoadSeqs(data = {'seq1': 'ACGTAA', 'seq2': 'ACG',
                            'seq3': 'ACGCGT'}, moltype = DNA, aligned=False)
        assert seq_coll.hasTerminalStops() == True
        seq_coll = LoadSeqs(data = {'seq1': 'ACGTAC', 'seq2': 'ACGACG',
                            'seq3': 'ACGCGT'}, moltype = DNA, aligned=False)
        assert seq_coll.hasTerminalStops() == False
        # alignments
        aln = LoadSeqs(data = {'seq1': 'ACGTAA', 'seq2': 'ACGCAA',
                            'seq3': 'ACGCGT'}, moltype = DNA)
        assert aln.hasTerminalStops() == True
        aln = LoadSeqs(data = {'seq1': 'ACGTAA', 'seq2': 'ACGTAG',
                            'seq3': 'ACGTGA'}, moltype = DNA)
        assert aln.hasTerminalStops() == True
        aln = LoadSeqs(data = {'seq1': 'ACGCAA', 'seq2': 'ACGCAA',
                            'seq3': 'ACGCGT'}, moltype = DNA)
        assert aln.hasTerminalStops() == False
    
    def test_slice(self):
        seqs = {'seq1': 'ACGTANGT', 'seq2': 'ACGTACGT', 'seq3': 'ACGTACGT'}
        alignment = LoadSeqs(data = seqs)
        short = {'seq1':'A', 'seq2':'A', 'seq3':'A'}
        self.assertEqual(alignment[0:1].todict(), short)
    
    def test_get_motifprobs(self):
        """calculation of motif probs"""
        seqs = {'seq1': 'ACGTANGT', 'seq2': '-CGTACGT', 'seq3': 'ACGTACGT'}
        aln = LoadSeqs(data = seqs, moltype=DNA)
        mprobs = aln.getMotifProbs(allow_gap=False)
        expected = {'A':5/22, 'T':6/22, 'C':5/22, 'G':6/22}
        self.assertEqual(mprobs, expected)
        mprobs = aln.getMotifProbs(allow_gap=True)
        expected = {'A':5/23, 'T':6/23, 'C':5/23, 'G':6/23, '-':1/23}
        self.assertEqual(mprobs, expected)
        mprobs = aln.getMotifProbs(allow_gap=False, include_ambiguity=True)
        expected = {'A':5.25/23, 'T':6.25/23, 'C':5.25/23, 'G':6.25/23}
        self.assertEqual(mprobs, expected)
        mprobs = aln.getMotifProbs(allow_gap=True, include_ambiguity=True)
        expected = {'A':5.25/24, 'T':6.25/24, 'C':5.25/24, 'G':6.25/24, '-':
                    1/24}
        self.assertEqual(mprobs, expected)
        seqs = {'seq1': 'ACGAANGA', 'seq2': '-CGAACGA', 'seq3': 'ACGAACGA'}
        aln = LoadSeqs(data = seqs, moltype=DNA)
        mprobs = aln.getMotifProbs(exclude_unobserved=True)
        expected = {'A':11/22, 'C':5/22, 'G':6/22}
        self.assertEqual(mprobs, expected)
        
    

# fileformats doesn't catch an exception when the file has no data!
class SequenceTestMethods(unittest.TestCase):
    """Testing Sequence methods"""
    def setUp(self):
        self.seq = Sequence(DNA,  'ATGACGTTGCGTAGCATAGCTCGA')
        
    def test_getlength(self):
        """testing getting length"""
        assert len(self.seq) == 24
        
    def test_getInMotifSize(self):
        """test accuracy of chunking various sizes"""
        self.assertEqual(self.seq.getInMotifSize(2),
                ['AT','GA','CG','TT','GC','GT','AG','CA','TA','GC','TC','GA'])
        self.assertEqual(self.seq.getInMotifSize(3),
                ['ATG','ACG','TTG','CGT','AGC','ATA','GCT','CGA'])
            
    def test_translate(self):
        """test of translating seqs"""
        seq = Sequence(DNA, 'ATGACGTTGCGTAGCATAGCTCGA').getTranslation()
        self.assertEqual(str(seq), 'MTLRSIAR')

    def test_ambig_translate(self):
        """test of translating seqs"""
        seq = Sequence(DNA, 'CGNTGN???---').getTranslation()
        self.assertEqual(str(seq), 'RX?-')

    def test_slidingWindows(self):
        """test sliding window along sequences"""
        result = []
        for bit in self.seq.slidingWindows(5,2):
            result+=[bit]
        self.assertEqual([str(x) for x in result],
                          ['ATGAC', 'GACGT', 'CGTTG', 'TTGCG', 'GCGTA',
                          'GTAGC', 'AGCAT', 'CATAG', 'TAGCT', 'GCTCG'])
        
        result = []
        for bit in self.seq.slidingWindows(5,1):
            result+=[bit]
        self.assertEqual([str(x) for x in result],
                          ['ATGAC', 'TGACG', 'GACGT', 'ACGTT', 'CGTTG',
                          'GTTGC', 'TTGCG', 'TGCGT', 'GCGTA', 'CGTAG',
                          'GTAGC', 'TAGCA', 'AGCAT', 'GCATA', 'CATAG',
                          'ATAGC', 'TAGCT', 'AGCTC', 'GCTCG', 'CTCGA'])
    
    def test_reversecomplement(self):
        """testing reversal and complementing of a sequence"""
        seq = Sequence(DNA, seq='ACTGTAA')
        rev = seq.reversecomplement()
        self.assertEqual(str(rev), 'TTACAGT')
        seq = Sequence(DNA, seq='ACTG-TAA')
        rev = seq.reversecomplement()
        self.assertEqual(str(rev), 'TTA-CAGT')
        #try amigbuities
        seq = Sequence(DNA, seq='ACHNRTAA')
        rev = seq.reversecomplement()
        self.assertEqual(str(rev), 'TTAYNDGT')
        
    def test_withoutTerminalStopCodon(self):
        """testing deleting terminal stop"""
        # for standard code
        seq = Sequence(DNA, seq='ACTTAA')
        seq2 = seq.withoutTerminalStopCodon()
        self.assertEqual(str(seq2), "ACT")
    
    def test_hasTerminalStop(self):
        """test check for terminal stop codons"""
        seq = Sequence(DNA, seq='ACTTAA')
        assert seq.hasTerminalStop() == True
        seq = Sequence(DNA, seq='ACTTAT') == False
        try:
            # only sequences with length divisible by 3 should work
            seq = Sequence(DNA, seq='ACTTA')
            seq.hasTerminalStop()
        except AssertionError:
            pass

if __name__ == '__main__':
    unittest.main()

