import os

from cogent import DNA
from cogent.util.unit_test import TestCase, main
from cogent.db.ensembl.host import HostAccount, get_ensembl_account
from cogent.db.ensembl.util import convert_strand
from cogent.db.ensembl.genome import Genome
from cogent.db.ensembl.sequence import _assemble_seq

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

Release = 56

if 'ENSEMBL_ACCOUNT' in os.environ:
    username, password = os.environ['ENSEMBL_ACCOUNT'].split()
    account = HostAccount('127.0.0.1', username, password)
else:
    account = get_ensembl_account(release=Release)

class GenomeTestBase(TestCase):
    human = Genome(Species="human", Release=Release, account=account)
    mouse = Genome(Species="mouse", Release=Release, account=account)
    rat = Genome(Species="rat", Release=Release, account=account)
    macaq = Genome(Species="macaque", Release=Release, account=account)
    brca2 = list(human.getGenesMatching(StableId="ENSG00000139618"))[0]

class TestGenome(GenomeTestBase):
    
    def test_other_features(self):
        """should correctly return record for ENSESTG00000000020"""
        est = self.human.getEstMatching(StableId='ENSESTG00000000020')
        direct = list(est)[0]
        ests = self.human.getFeatures(feature_types='est', CoordName=8,
                                                Start=105711700, End=105715900)
        stable_ids = [est.StableId for est in ests]
        self.assertContains(stable_ids, direct.StableId)
    
    def test_genome_comparison(self):
        """different genome instances with same CoreDb connection are equal"""
        h2 = Genome(Species='human', Release=Release, account=account)
        self.assertEquals(self.human, h2)
    
    def test_make_location(self):
        """should correctly make a location for an entire chromosome"""
        loc = self.human.makeLocation(CoordName=1)
        self.assertEquals(len(loc), 249250621)
    
    def test_get_region(self):
        """should return a generic region that extracts correct sequence"""
        chrom = 1
        Start = 11137
        End = Start+20
        region = self.human.getRegion(CoordName=chrom, Start=Start, End=End,
                        ensembl_coord=True)
        self.assertEquals(region.Location.Start, Start-1)
        self.assertEquals(region.Location.End, End)
        self.assertEquals(region.Location.CoordName, str(chrom))
        self.assertEquals(region.Location.CoordType, 'chromosome')
        self.assertEquals(region.Seq, 'ACCTCAGTAATCCGAAAAGCC')
    
    def test_get_assembly_exception_region(self):
        """should return correct sequence for region with an assembly
        exception"""
        ##old:chrY:57767412-57767433; New: chrY:59358024-59358045
        region = self.human.getRegion(CoordName = "Y", Start = 59358024, 
                            End = 59358045, Strand = 1, ensembl_coord = True)
        self.assertEquals(str(region.Seq), 'CGAGGACGACTGGGAATCCTAG')
    
    def test_getting_annotated_seq(self):
        """a region should return a sequence with the correct annotation"""
        new_loc = self.brca2.Location.resized(-100, 100)
        region = self.human.getRegion(region=new_loc)
        annot_seq = region.getAnnotatedSeq(feature_types='gene')
        gene_annots = annot_seq.getAnnotationsMatching('gene')
        self.assertEquals(gene_annots[0].Name, self.brca2.Symbol)
    
    def test_correct_feature_type_id_cache(self):
        """should obtain the feature type identifiers without failure"""
        self.assertNotEquals(self.human._feature_type_ids.CpGisland, None)
    
    def test_strand_conversion(self):
        """should consistently convert strand info"""
        self.assertEquals(convert_strand(None), 1)
        self.assertEquals(convert_strand(-1), -1)
        self.assertEquals(convert_strand(1), 1)
        self.assertEquals(convert_strand('-'), -1)
        self.assertEquals(convert_strand('+'), 1)
        self.assertEquals(convert_strand(-1.0), -1)
        self.assertEquals(convert_strand(1.0), 1)
    

class TestGene(GenomeTestBase):
    def _eval_brca2(self, brca2):
        """tests all attributes correctly created"""
        self.assertEquals(brca2.Symbol.lower(), 'brca2')
        self.assertEquals(brca2.StableId, 'ENSG00000139618')
        self.assertEquals(brca2.BioType.lower(), 'protein_coding')
        self.assertContains(brca2.Description.lower(), 'breast cancer')
        self.assertEquals(brca2.Status, 'KNOWN')
        self.assertEquals(brca2.CanonicalTranscript.StableId,
                        'ENST00000380152')
        # note length can change between genome builds
        self.assertGreaterThan(len(brca2), 83700)
        transcript = brca2.getMember('ENST00000380152')
        self.assertEquals(transcript.getCdsLength(),len(transcript.Cds))
    
    def test_get_genes_by_stable_id(self):
        """if get gene by stable_id, attributes should be correctly
        constructed"""
        self._eval_brca2(self.brca2)
    
    def test_get_exons(self):
        """transcript should return correct exons for brca2"""
        transcript = self.brca2.getMember('ENST00000380152')
        self.assertEquals(len(transcript.TranslatedExons), 26)
        self.assertEquals(len(transcript.Cds), 3419*3)
        self.assertEquals(len(transcript.ProteinSeq), 3418)
    
    def test_translated_exons(self):
        """should correctly translate a gene with 2 exons but 1st exon
        transcribed"""
        gene = self.mouse.getGenesMatching(StableId='ENSMUSG00000036136')
        gene = list(gene)[0]
        transcript = gene.getMember('ENSMUST00000041133')
        self.assertTrue(len(transcript.ProteinSeq) > 0)
        # now one on the - strand
        gene = self.mouse.getGenesMatching(StableId='ENSMUSG00000045912')
        gene = list(gene)[0]
        transcript = gene.Transcripts[0]
        self.assertTrue(len(transcript.ProteinSeq) > 0)
    
    def test_failed_ensembl_annotation(self):
        """we demonstrate a failed annotation by ensembl"""
        # I'm including this to demonstrate that Ensembl coords aren't
        # always right. This case has a macaque gene which we correctly
        # infer the CDS boundaries for according to Ensembl, but the CDS
        # length is not divisible by 3, hence our code fails.
        gene = self.macaq.getGenesMatching(StableId='ENSMMUG00000001551')
        gene = list(gene)[0]
        transcript = gene.getMember('ENSMMUT00000002194')
        try:
            # this fails because of the above length issue
            prot_seq = transcript.ProteinSeq
        except AssertionError:
            pass
        # now if you slice the CDS to be divisible by 3, you can get the same
        # protein sequence as that at Ensembl
        l = 3 * (transcript.getCdsLength() / 3)
        trunc_cds = transcript.Cds[:l]
        prot_seq = trunc_cds.getTranslation()
        self.assertEquals(str(prot_seq),
            'MPSSPLRVAVVCSSNQNRSMEAHNILSKRGFSVRSFGTGTHVKLPGPAPDKPNVYDFKTT'\
               'YDQMYNDLLRKDKELYTQNGILHMLDRNKRIKPRPERFQNCKDLFDLILTCEERVY')
    
    def test_gene_transcripts(self):
        """should return multiple transcripts"""
        stable_id = 'ENSG00000012048'
        gene = list(self.human.getGenesMatching(StableId=stable_id))[0]
        self.assertTrue(len(gene.Transcripts) > 1)
        # .. and correctly construct the Cds and location
        for transcript in gene.Transcripts:
            self.assertTrue(transcript.getCdsLength()>0)
            self.assertEquals(transcript.Location.CoordName,'17')
    
    def test_get_longest_cds_transcript(self):
        """should correctly return transcript with longest cds"""
        stable_id = 'ENSG00000178591'
        gene = list(self.human.getGenesMatching(StableId=stable_id))[0]
        ts = gene.getLongestCdsTranscript()
        self.assertEquals(ts.getCdsLength(), max(gene.getCdsLengths()))
    
    def test_rna_transcript_cds(self):
        """should return a Cds for an RNA gene too"""
        rna_gene = self.human.getGenesMatching(StableId='ENSG00000210049')
        rna_gene = list(rna_gene)[0]
        self.assertTrue(rna_gene.Transcripts[0].getCdsLength() > 0)
    
    def test_gene_annotation(self):
        """should correctly annotated a sequence"""
        annot_seq = self.brca2.getAnnotatedSeq(feature_types='gene')
        gene_annots = annot_seq.getAnnotationsMatching('gene')
        self.assertEquals(gene_annots[0].Name, self.brca2.Symbol)
    
    def test_get_by_symbol(self):
        """selecting a gene by it's HGNC symbol should correctly populate all
        specified attributes"""
        results = self.human.getGenesMatching(Symbol="BRCA2")
        for snp in results:
            self._eval_brca2(snp)
    
    def test_get_by_description(self):
        """if get by description, all attributes should be correctly
        constructed"""
        description='breast cancer type 2'
        results = list(self.human.getGenesMatching(Description=description))
        self.assertEquals(len(results), 1)
        self._eval_brca2(results[0])
    
    def test_get_member(self):
        """should return correct exon and translated exon"""
        transcript = self.brca2.getMember('ENST00000380152')
        # just returns the first
        exon_id = 'ENSE00001484009'
        exon = transcript.getMember(exon_id)
        trans_exon = transcript.getMember(exon_id,'TranslatedExons')
        self.assertEquals(exon.StableId, exon_id)
        self.assertEquals(trans_exon.StableId, exon_id)
        # we check we got Exon in the first call and TranslatedExon in the
        # second using the fact that the Exons entry is longer than the
        # TranslatedExons one
        self.assertGreaterThan(len(exon), len(trans_exon))
    
    def test_get_by_biotype(self):
        results = list(self.human.getGenesMatching(BioType='Mt_tRNA', like=False))
        self.assertEquals(len(results), 22)
        results = list(self.human.getGenesMatching(BioType='Mt_tRNA', like=True))
        self.assertEquals(len(results), 602)
    
    def test_get_by_decsr_biotype(self):
        """combining the description and biotype should return a result"""
        results = list(self.human.getGenesMatching(BioType="protein_coding",
                    Description="cancer"))
        self.assertTrue(len(results) > 100)
    

class TestVariation(GenomeTestBase):
    snp_names =  ['rs34213141', 'rs12791610', 'rs10792769', 'rs11545807', 'rs11270496']
    snp_nt_alleles = ['G/A', 'C/T', 'A/G', 'C/A', 'CAGCTCCAGCTC/-']
    snp_aa_alleles = ['G/R', 'P/L', 'Y/C', "V/F", "GAGAV/V"]
    snp_effects = ['INTRONIC']*3+[['REGULATORY_REGION', 'NON_SYNONYMOUS_CODING']]+['NON_SYNONYMOUS_CODING']
    snp_nt_len = [1, 1, 1, 1, 12]
    map_weights = [1,1,1,1,1]
    snp_flanks = [
     ('CTGAGGTGAGCCAGCGTTGGAGCTGTTTTTCCTTTCAGTATGAATTCCACAAGGAAATCATCTCAGGAGGAAGGGCTCATACTTGGATCCAGAAAATATCAACATAGCCAAAGAAAAACAATCAAGACATACCTCCAGGAGCTGTGTAACAGCAACCGGAAAGAGAAACAATGGTGTGTTCCTATGTGGGATATAAAGAGCCGGGGCTCAGGGGGCTCCACACCTGCACCTCCTTCTCACCTGCTCCTCTACCTGCTCCACCCTCAATCCACCAGAACCATGGGCTGCTGTGGCTGCTCC',
      'GAGGCTGTGGCTCCAGCTGTGGAGGCTGTGACTCCAGCTGTGGGAGCTGTGGCTCTGGCTGCAGGGGCTGTGGCCCCAGCTGCTGTGCACCCGTCTACTGCTGCAAGCCCGTGTGCTGCTGTGTTCCAGCCTGTTCCTGCTCTAGCTGTGGCAAGCGGGGCTGTGGCTCCTGTGGGGGCTCCAAGGGAGGCTGTGGTTCTTGTGGCTGCTCCCAGTGCAGTTGCTGCAAGCCCTGCTGTTGCTCTTCAGGCTGTGGGTCATCCTGCTGCCAGTGCAGCTGCTGCAAGCCCTACTGCTCCC'),
     ('GAAAATATCAACATAGCCAAAGAAAAACAATCAAGACATACCTCCAGGAGCTGTGTAACAGCAACCGGAAAGAGAAACAATGGTGTGTTCCTATGTGGGATATAAAGAGCCGGGGCTCAGGGGGCTCCACACCTGCACCTCCTTCTCACCTGCTCCTCTACCTGCTCCACCCTCAATCCACCAGAACCATGGGCTGCTGTGGCTGCTCCGGAGGCTGTGGCTCCAGCTGTGGAGGCTGTGACTCCAGCTGTGGGAGCTGTGGCTCTGGCTGCAGGGGCTGTGGCCCCAGCTGCTGTGCAC',
      'CGTCTACTGCTGCAAGCCCGTGTGCTGCTGTGTTCCAGCCTGTTCCTGCTCTAGCTGTGGCAAGCGGGGCTGTGGCTCCTGTGGGGGCTCCAAGGGAGGCTGTGGTTCTTGTGGCTGCTCCCAGTGCAGTTGCTGCAAGCCCTGCTGTTGCTCTTCAGGCTGTGGGTCATCCTGCTGCCAGTGCAGCTGCTGCAAGCCCTACTGCTCCCAGTGCAGCTGCTGTAAGCCCTGTTGCTCCTCCTCGGGTCGTGGGTCATCCTGCTGCCAATCCAGCTGCTGCAAGCCCTGCTGCTCATCCTC'),
     ('ATCAACATAGCCAAAGAAAAACAATCAAGACATACCTCCAGGAGCTGTGTAACAGCAACCGGAAAGAGAAACAATGGTGTGTTCCTATGTGGGATATAAAGAGCCGGGGCTCAGGGGGCTCCACACCTGCACCTCCTTCTCACCTGCTCCTCTACCTGCTCCACCCTCAATCCACCAGAACCATGGGCTGCTGTGGCTGCTCCGGAGGCTGTGGCTCCAGCTGTGGAGGCTGTGACTCCAGCTGTGGGAGCTGTGGCTCTGGCTGCAGGGGCTGTGGCCCCAGCTGCTGTGCACCCGTCT',
      'CTGCTGCAAGCCCGTGTGCTGCTGTGTTCCAGCCTGTTCCTGCTCTAGCTGTGGCAAGCGGGGCTGTGGCTCCTGTGGGGGCTCCAAGGGAGGCTGTGGTTCTTGTGGCTGCTCCCAGTGCAGTTGCTGCAAGCCCTGCTGTTGCTCTTCAGGCTGTGGGTCATCCTGCTGCCAGTGCAGCTGCTGCAAGCCCTACTGCTCCCAGTGCAGCTGCTGTAAGCCCTGTTGCTCCTCCTCGGGTCGTGGGTCATCCTGCTGCCAATCCAGCTGCTGCAAGCCCTGCTGCTCATCCTCAGGCTG'), 
      ('GCTGAAGAAACCATTTCAAACAGGATTGGAATAGGGAAACCCGGCACTCAGCTCGGCGCAAGCCGGCGGTGCCTTCAGACTAGAGAGCCTCTCCTCCGGTGCGCTGCAAGTAGGGCCTCGGCTCGAGGTCAACATTCTAGTTGTCCAGCGCTCCCTCTCCGGCACCTCGGTGAGGCTAGTTGACCCGACAGGCGCGGATCATGAGCAGCTGCAGGAGAATGAAGAGCGGGGACGTAATGAGGCCGAACCAGAGCTCCCGAGTCTGCTCCGCCAGCTTCTGGCACAACAGCATCTCGAAGA',
'GAACTTGAGACTCAGGACCGTAAGTACCCAGAAAAGGCGGAGCACCGCCAGCCGCTTCTCTCCATCCTGGAAGAGGCGCACGGACACGATGGTGGTGAAGTAGGTGCTGAGCCCGTCAGCGGCGAAGAAAGGCACGAACACGTTCCACCAGGAGAGGCCCGGGACCAGGCCATCCACACGCAGTGCCAGCAGCACAGAGAACACCAACAGGGCCAGCAGGTGCACGAAGATCTCGAAGGTGGCGAAGCCTAGCCACTGCACCAGCTCCCGGAGCGAGAAGAGCATCGCGCCCGTTGAGCG')]
    def test_get_variation_by_symbol(self):
        """should return correct snp when query genome by symbol"""
        # supplement this test with some synonymous snp's, where they have no
        # peptide alleles
        for i in range(4):
            snp = list(self.human.getVariation(Symbol=self.snp_names[i]))[0]
            self.assertEquals(snp.Symbol, self.snp_names[i])
            self.assertEquals(snp.Effect, self.snp_effects[i])
            self.assertEquals(snp.Alleles, self.snp_nt_alleles[i])
            self.assertEquals(snp.MapWeight, self.map_weights[i])
    
    def test_num_alleles(self):
        """should correctly infer the number of alleles"""
        for i in range(4):
            snp = list(self.human.getVariation(Symbol=self.snp_names[i]))[0]
            self.assertEquals(len(snp), self.snp_nt_len[i])
    
    def test_get_peptide_alleles(self):
        """should correctly infer the peptide alles"""
        for i in range(4):
            snp = list(self.human.getVariation(Symbol=self.snp_names[i]))[0]
            if snp.Effect == 'INTRONIC':
                continue
            
            self.assertEquals(snp.PeptideAlleles, self.snp_aa_alleles[i])
    
    def test_validation_status(self):
        """should return correct validation status"""
        data = (('rs34213141', set(['freq']), lambda x: set([x])),
                ('rs12791610', set(['hapmap']), lambda x: set([x])),
                ('rs10792769', set(['cluster', 'freq', 'doublehit']), set))
        for name, status, conv in data:
            snp = list(self.human.getVariation(Symbol=name))[0]
            self.assertTrue(status <= conv(snp.Validation))
    
    def test_get_flanking_seq(self):
        """should correctly get the flanking sequence"""
        for i in range(4): # only have flanking sequence for 3
            snp = list(self.human.getVariation(Symbol=self.snp_names[i]))[0]
            self.assertEquals(snp.FlankingSeq, self.snp_flanks[i])
    
    def test_variation_seq(self):
        """should return the sequence for a Variation snp if asked"""
        snp = list(self.human.getVariation(Symbol=self.snp_names[0]))[0]
        self.assertContains(snp.Alleles, str(snp.Seq))
    
    def test_get_validation_condition(self):
        """simple test of SNP validation status"""
        snp_status = [('rs94', False), ('rs90', True)]
        for symbol, status in snp_status:
            snp = list(self.human.getVariation(Symbol=symbol, validated=True))
            self.assertEquals(snp!=[], status)
    

class TestFeatures(GenomeTestBase):
    def setUp(self):
        self.igf2 = list(
                  self.human.getGenesMatching(StableId='ENSG00000167244'))[0]
    
    def test_CpG_island(self):
        """should return correct CpG islands"""
        CpGislands = self.human.getFeatures(region=self.igf2,
                            feature_types='CpG')
        expected_stats = [(630, 757), (652, 537), (3254, 3533)]
        obs_stats = [(int(island.Score), len(island)) \
                                            for island in CpGislands]
        obs_stats.sort()
        self.assertTrue(set(expected_stats) & set(obs_stats) != set())
    
    def test_get_multiple_features(self):
        """should not fail to get multiple feature types"""
        regions =\
            self.human.getFeatures(feature_types=['repeat','gene','cpg'],
                                    CoordName=1, Start=869936,End=901867)
        for region in regions:
            pass
    
    def test_repeats(self):
        """should correctly return a repeat"""
        loc = self.igf2.Location.resized(-1000, 1000)
        repeats = list(self.human.getFeatures(
                                    region=loc, feature_types='repeat'))
        self.assertTrue(len(repeats) >= 4)
    
    def test_genes(self):
        """should correctly identify igf2 within a region"""
        loc = self.igf2.Location.resized(-1000, 1000)
        genes = self.human.getFeatures(region=loc, feature_types='gene')
        symbols = [g.Symbol.lower() for g in genes]
        self.assertContains(symbols, self.igf2.Symbol.lower())
    
    def test_other_genes(self):
        """docstring for test_other_genes"""
        mouse = self.mouse.getRegion(CoordName='5', Start=150791005,
                                        End=150838512, Strand='-')
        rat = self.rat.getRegion(CoordName='12', Start=4282534, End=4324019,
                                        Strand='+')
        for region in [mouse, rat]:
            features = region.getFeatures(feature_types=['gene'])
            ann_seq = region.getAnnotatedSeq(feature_types='gene')
            genes = ann_seq.getAnnotationsMatching('gene')
            self.assertTrue(genes != [])
    
    def test_get_variation_feature(self):
        """should correctly return variation features within a region"""
        snps = self.human.getFeatures(feature_types='variation',
                                        region=self.brca2)
        # snp coordname, start, end should satsify constraints of brca2 loc
        c = 0
        loc = self.brca2.Location
        for snp in snps:
            self.assertEquals(snp.Location.CoordName, loc.CoordName)
            self.assertTrue(loc.Start < snp.Location.Start < loc.End)
            c += 1
            if c == 2:
                break
    
    def test_gene_feature_data_correct(self):
        """should apply gene feature data in a manner consistent with strand
        and the Cogent sequence annotations slice should return the same
        result"""
        plus = list(self.human.getFeatures(feature_types='gene',
                                           CoordName=13,
                                           Start=31787610,
                                           End=31871820))[0]
        minus = plus.Location.copy()
        minus.Strand *= -1
        minus = self.human.getRegion(region = minus)
        # get Sequence
        plus_seq = plus.getAnnotatedSeq(feature_types='gene')
        minus_seq = minus.getAnnotatedSeq(feature_types='gene')
        # the seqs should be the rc of each other
        self.assertEquals(str(plus_seq), str(minus_seq.rc()))
        # the Cds, however, from the annotated sequences should be identical
        plus_cds = plus_seq.getAnnotationsMatching('CDS')[0]
        minus_cds = minus_seq.getAnnotationsMatching('CDS')[0]
        self.assertEquals(str(plus_cds.getSlice()),str(minus_cds.getSlice()))
    
    def test_other_feature_data_correct(self):
        """should apply CpG feature data in a manner consistent with strand"""
        human = self.human
        coord = dict(CoordName=11, Start=2165124,End=2165724)
        exp_coord = dict(CoordName=11, Start=2165136, End=2165672)
        exp_loc = human.getRegion(Strand=1, ensembl_coord=True, **exp_coord)
        exp = exp_loc.Seq
        
        ps_feat = human.getRegion(Strand=1, **coord)
        ms_feat = human.getRegion(Strand=-1, **coord)
        
        ps_seq = ps_feat.getAnnotatedSeq(feature_types='CpG')
        ps_cgi = ps_seq.getAnnotationsMatching('CpGisland')[0]
        
        self.assertEquals(ps_feat.Seq, ms_feat.Seq.rc())
        
        self.assertEquals(ps_cgi.getSlice().rc(), exp)
        ms_seq = ms_feat.getAnnotatedSeq(feature_types='CpG')
        ms_cgi = ms_seq.getAnnotationsMatching('CpGisland')[0]
        
        self.assertEquals(ms_cgi.getSlice(), ps_cgi.getSlice())
    
    def test_other_repeat(self):
        """should apply repeat feature data in a manner consistent with strand"""
        coord=dict(CoordName=13, Start=32890200, End=32890500)
        ps_repeat = self.human.getRegion(Strand=1, **coord)
        ms_repeat = self.human.getRegion(Strand=-1, **coord)
        exp = DNA.makeSequence('CTTACTGTGAGGATGGGAACATTTTACAGCTGTGCTG'\
          'TCCAAACCGGTGCCACTAGCCACATTAAGCACTCGAAACGTGGCTAGTGCGACTAGAGAAGAGGA'\
          'TTTTCATACGATTTAGTTTCAATCACGCTAACCAGTGACGCGTGGCTAGTGG')
        
        self.assertEquals(ms_repeat.Seq, ps_repeat.Seq.rc())
        
        ps_annot_seq = ps_repeat.getAnnotatedSeq(feature_types='repeat')
        ms_annot_seq = ms_repeat.getAnnotatedSeq(feature_types='repeat')
        ps_seq = ps_annot_seq.getAnnotationsMatching('repeat')[0]
        ms_seq = ms_annot_seq.getAnnotationsMatching('repeat')[0]
        self.assertEquals(ms_seq.getSlice(), ps_seq.getSlice())
        self.assertEquals(ps_seq.getSlice(), exp)
    
    def test_get_features_from_nt(self):
        """should correctly return the encompassing gene from 1nt"""
        snp = list(self.human.getVariation(Symbol='rs34213141'))[0]
        gene=list(self.human.getFeatures(feature_types='gene',region=snp))[0]
        self.assertEquals(gene.StableId, 'ENSG00000204572')
    

class TestAssembly(TestCase):
    
    def test_assemble_seq(self):
        """should correctly fill in a sequence with N's"""
        expect = DNA.makeSequence("NAAAAANNCCCCCNNGGGNNN")
        frags = ["AAAAA","CCCCC","GGG"]
        positions = [(11, 16), (18, 23), (25, 28)]
        self.assertEqual(_assemble_seq(frags, 10, 31, positions), expect)
        positions = [(1, 6), (8, 13), (15, 18)]
        self.assertEqual(_assemble_seq(frags, 0, 21, positions), expect)
        # should work with:
        # start matches first frag start
        expect = DNA.makeSequence("AAAAANNCCCCCNNGGGNNN")
        positions = [(0, 5), (7, 12), (14, 17)]
        self.assertEqual(_assemble_seq(frags, 0, 20, positions), expect)
        # end matches last frag_end
        expect = DNA.makeSequence("NAAAAANNCCCCCNNGGG")
        positions = [(11, 16), (18, 23), (25, 28)]
        self.assertEqual(_assemble_seq(frags, 10, 28, positions), expect)
        # both start and end matched
        expect = DNA.makeSequence("AAAAANNCCCCCNNGGG")
        positions = [(10, 15), (17, 22), (24, 27)]
        self.assertEqual(_assemble_seq(frags, 10, 27, positions), expect)
        # one frag
        expect = DNA.makeSequence(''.join(frags))
        positions = [(10, 23)]
        self.assertEqual(_assemble_seq([''.join(frags)],10,23,positions),
                                expect)
    
if __name__ == "__main__":
    main()
