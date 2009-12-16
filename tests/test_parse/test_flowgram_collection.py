__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2009, The Cogent Project"
__credits__ = ["Jens Reeder","Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Julia Goodrich"
__email__ = "julia.goodrich@colorado.edu"
__Status__ = "Development"


from cogent.util.unit_test import TestCase, main
from types import GeneratorType
from numpy import array, transpose
from cogent.core.sequence import Sequence

from cogent.parse.flowgram_collection import FlowgramCollection, flows_from_array,\
     flows_from_generic,flows_from_kv_pairs,flows_from_empty,flows_from_dict,\
     flows_from_sff,assign_sequential_names,flows_from_flowCollection,\
     pick_from_prob_density, seqs_to_flows

from cogent.parse.flowgram import Flowgram
from cogent.core.alignment import SequenceCollection
from tempfile import mktemp
from os import remove

class flowgram_tests(TestCase):
    """Tests of top-level functions."""
    
    def test_flows_from_array(self):
        """flows_from_array should return chars, and successive indices."""
        a = array([[0,1,2],[2,1,0]])    #three 2-char seqs
        obs_a, obs_labels, obs_info = flows_from_array(a)
        #note transposition
        self.assertEqual(obs_a, [array([0,2]), array([1,1]), array([2,0])])
        self.assertEqual(obs_labels, None)
        self.assertEqual(obs_info, None)
        
    def test_flows_from_generic(self):
        """flows_from_flow should initialize from list of flowgram objects"""
        c = Flowgram('0.0 1.1 3.0 1.0', Name='a')
        b = Flowgram('0.5 1.0 4.0 0.0', Name = 'b')
        obs_a, obs_labels, obs_info = flows_from_generic([c,b])
        self.assertEqual(map(str,obs_a), ['0.0\t1.1\t3.0\t1.0',
                                          '0.5\t1.0\t4.0\t0.0'])
        self.assertEqual(obs_labels, ['a','b'])
        self.assertEqual(obs_info, [None,None])

        f = ['0.0 1.1 3.0 1.0','0.5 1.0 4.0 0.0']
        obs_a, obs_labels, obs_info = flows_from_generic(f)
        self.assertEqual(map(str,obs_a), ['0.0 1.1 3.0 1.0',
                                          '0.5 1.0 4.0 0.0'])
        self.assertEqual(obs_labels, [None,None])
        self.assertEqual(obs_info, [None,None])

    def test_flows_from_flowCollection(self):
        """flows_from_flowCollection should init from existing collection"""
        c = FlowgramCollection({'a':'0.0 1.1 3.0 1.0','b':'0.5 1.0 4.0 0.0'})
        obs_a, obs_labels, obs_info = flows_from_flowCollection(c)
        self.assertEqual(map(str,obs_a), ['0.0\t1.1\t3.0\t1.0',
                                          '0.5\t1.0\t4.0\t0.0'])
        self.assertEqual(obs_labels, ['a','b'])
        self.assertEqual(obs_info, [None,None])

    def test_flows_from_kv_pairs(self):
        """seqs_from_kv_pairs should initialize from key-value pairs"""
        c = [['a','0.0 1.1 3.0 1.0'],['b','0.5 1.0 4.0 0.0']]
        obs_a, obs_labels, obs_info = flows_from_kv_pairs(c)
        self.assertEqual(map(str,obs_a), ['0.0 1.1 3.0 1.0','0.5 1.0 4.0 0.0'])
        self.assertEqual(obs_labels, ['a','b'])
        self.assertEqual(obs_info, [None,None])

        c =[['a',Flowgram('0.0 1.1 3.0 1.0')],['b',Flowgram('0.5 1.0 4.0 0.0')]]
        obs_a, obs_labels, obs_info = flows_from_kv_pairs(c)
        self.assertEqual(map(str,obs_a), ['0.0\t1.1\t3.0\t1.0','0.5\t1.0\t4.0\t0.0'])
        self.assertEqual(obs_labels, ['a','b'])
        self.assertEqual(obs_info, [None,None])

    def test_flows_from_empty(self):
        """flowss_from_empty should always raise ValueError"""
        self.assertRaises(ValueError, flows_from_empty, 'xyz')

    def test_flows_from_dict(self):
        """flows_from_dict should init from dictionary"""
        c = {'a':'0.0 1.1 3.0 1.0','b':'0.5 1.0 4.0 0.0'}
        obs_a, obs_labels, obs_info = flows_from_dict(c)
        self.assertEqual(map(str,obs_a), ['0.0 1.1 3.0 1.0','0.5 1.0 4.0 0.0'])
        self.assertEqual(obs_labels, ['a','b'])
        self.assertEqual(obs_info, [None,None])

        c ={'a':Flowgram('0.0 1.1 3.0 1.0'),'b':Flowgram('0.5 1.0 4.0 0.0')}
        obs_a, obs_labels, obs_info = flows_from_dict(c)
        self.assertEqual(map(str,obs_a), ['0.0\t1.1\t3.0\t1.0','0.5\t1.0\t4.0\t0.0'])
        self.assertEqual(obs_labels, ['a','b'])
        self.assertEqual(obs_info, [None,None])
        
    def test_pick_from_prob_density(self):
        """Should take bin probabilitys and bin size and return random"""
        i = pick_from_prob_density([0,1.0,0,0],1)
        self.assertEqual(i,1)

        i = pick_from_prob_density([1.0,0,0,0],.01)
        self.assertEqual(i,0.0)

    def test_seqs_to_flows(self):
        """seqs_to_flows should take a list of seqs and probs and return """
        seqs = [('a','ATCGT'), ('b','ACCCAG'), ('c','GTAATG')]
        a = SequenceCollection(seqs)

        flows = seqs_to_flows(a.items())
        assert isinstance(flows,FlowgramCollection)
        
        for f,i in zip(flows,['0.0 1.0 0.0 0.0 1.0 0.0 1.0 1.0 1.0 0.0 0.0 0.0',
                            '0.0 1.0 3.0 0.0 0.0 1.0 0.0 1.0',
                            '0.0 0.0 0.0 1.0 1.0 2.0 0.0 0.0 1.0 0.0 0.0 1.0']):
            self.assertEqual(f,i)

        probs ={0:[1.0,0,0,0,0],1:[0,1.0,0,0,0],2:[0,0,1.0,0,0],3:[0,0,0,1.0,0]}
        
        flows = seqs_to_flows(a.items(), probs = probs, bin_size = 1.0)
        assert isinstance(flows,FlowgramCollection)
        
        for f,i in zip(flows,['0.0 1.0 0.0 0.0 1.0 0.0 1.0 1.0 1.0 0.0 0.0 0.0',
                            '0.0 1.0 3.0 0.0 0.0 1.0 0.0 1.0',
                            '0.0 0.0 0.0 1.0 1.0 2.0 0.0 0.0 1.0 0.0 0.0 1.0']):
            self.assertEqual(f,i)

            
        
    
class FlowgramCollectionTests(TestCase):
    """Tests sff parser functions"""

    Class = FlowgramCollection
  

    def test_guess_input_type(self):
        """  _guess_input_type should figure out data type correctly"""
        git = self.unordered._guess_input_type
        self.assertEqual(git(self.unordered), 'flowcoll')
        self.assertEqual(git(['0.0 1.1 3.0 1.0','0.5 1.0 4.0 0.0']), 'generic')
        self.assertEqual(git([Flowgram('0.0 1.1 3.0 1.0'),
                              Flowgram('0.5 1.0 4.0 0.0')]), 'generic')
        self.assertEqual(git([[1,2],[4,5]]), 'kv_pairs') #precedence over generic
        self.assertEqual(git([('a',Flowgram('0.0 1.1 3.0 1.0')),
                              ('b',Flowgram('0.5 1.0 4.0 0.0'))]), 'kv_pairs')
        self.assertEqual(git([[1,2,3],[4,5,6]]), 'generic')
        self.assertEqual(git(array([[1,2,3],[4,5,6]])), 'array')
        self.assertEqual(git({'a':'0.0 1.1 3.0 1.0'}), 'dict')
        self.assertEqual(git({'a':Flowgram('0.0 1.1 3.0 1.0')}), 'dict')
        self.assertEqual(git([]), 'empty')
        self.assertEqual(git('Common Header'), 'sff')


    def test_init_pairs(self):
        """FlowgramCollection init from list of (key,val) should work"""
        Flows = [['a','0.0 1.1 3.0 1.0'],['b','0.5 1.0 4.0 0.0']]
        a = self.Class(Flows)
        self.assertEqual(len(a.NamedFlows), 2)
        self.assertEqual(a.NamedFlows['a'], '0.0 1.1 3.0 1.0')
        self.assertEqual(a.NamedFlows['b'], '0.5 1.0 4.0 0.0')
        self.assertEqual(a.Names, ['a','b'])
        self.assertEqual(list(a.flows), ['0.0 1.1 3.0 1.0','0.5 1.0 4.0 0.0'])


    def test_init_aln(self):
        """FlowgramCollection should init from existing Collections"""
        start = self.Class([['a','0.0 1.1 3.0 1.0'],['b','0.5 1.0 4.0 0.0']])
        exp = self.Class([['a','0.0 1.1 3.0 1.0'],['b','0.5 1.0 4.0 0.0']])
        f = self.Class(start)
        self.assertEqual(f, exp)
    test_init_aln.__doc__ = Class.__name__ + test_init_aln.__doc__
 
    def test_init_dict(self):
        """FlowgramCollection init from dict should work as expected"""
        d = {'a':'0.0 1.1 3.0 1.0','b':'0.5 1.0 4.0 0.0'}
        a = self.Class(d)
        self.assertEqual(a, d)
        self.assertEqual(a.NamedFlows.items(), d.items())

    def test_init_name_mapped(self):
        """FlowgramCollection init should allow name mapping function"""
        d = {'a':'0.0 1.1 3.0 1.0','b':'0.5 1.0 4.0 0.0'}
        f = lambda x: x.upper()
        a = self.Class(d, name_conversion_f=f)
        self.assertNotEqual(a, d)
        self.assertNotEqual(a.NamedFlows.items(), d.items())
        d_upper = {'A':'0.0 1.1 3.0 1.0','B':'0.5 1.0 4.0 0.0'}
        self.assertEqual(a, d_upper)
        self.assertEqual(a.NamedFlows.items(), d_upper.items())
        

    def test_init_flow(self):
        """FlowgramCollection init from list of flowgrams should use indices
        as keys"""
        f1 = Flowgram('0.0 1.1 3.0 1.0')
        f2 = Flowgram('0.5 1.0 4.0 0.0')
        flows = [f1,f2]
        a = self.Class(flows)
        self.assertEqual(len(a.NamedFlows), 2)
        self.assertEqual(a.NamedFlows['seq_0'], '0.0 1.1 3.0 1.0')
        self.assertEqual(a.NamedFlows['seq_1'], '0.5 1.0 4.0 0.0')
        self.assertEqual(a.Names, ['seq_0','seq_1'])
        self.assertEqual(list(a.Flows), ['0.0 1.1 3.0 1.0','0.5 1.0 4.0 0.0'])


    def test_flows_from_sff(self):
        """flow_from_sff should init from sff iterator"""
        s = self.rec
        f = self.Class(s)
        self.assertEqual(f.NamedFlows['FIQU8OX05GCVRO'], self.flow)


    def test_init_duplicate_keys(self):
        """FlowgramCollection init from kv pairs should fail on dup. keys"""
        f = [['a','0.0 1.1 3.0 1.0'],['b','0.5 1.0 4.0 0.0'],
             ['b','1.5 2.0 0.0 0.5']]
        self.assertRaises(ValueError, self.Class, f)
        self.assertEqual(self.Class(f, remove_duplicate_names=True).Names,
                         ['a','b'])

    def test_init_ordered(self):
        """FlowgramCollection should iter over flows correctly, ordered too"""
        first = self.ordered1
        sec = self.ordered2
        un = self.unordered

        self.assertEqual(first.Names, ['a','b'])
        self.assertEqual(sec.Names, ['b', 'a'])
        self.assertEqual(un.Names, un.NamedFlows.keys())

        first_list = list(first.flow_str)
        sec_list = list(sec.flow_str)
        un_list = list(un.flow_str)

        self.assertEqual(first_list, ['0.0 1.1 3.0 1.0','0.5 1.0 4.0 0.0'])
        self.assertEqual(sec_list, ['0.5 1.0 4.0 0.0', '0.0 1.1 3.0 1.0'])
    
        #check that the unordered seq matches one of the lists
        self.assertTrue((un_list == first_list) or (un_list == sec_list))
        self.assertNotEqual(first_list, sec_list)
        
    def test_flow_str(self):
        """FlowgramCollection flow_str prop returns flows in correct order."""
        first = self.ordered1
        sec = self.ordered2
        un = self.unordered

        first_list = list(first.flow_str)
        sec_list = list(sec.flow_str)
        un_list = list(un.flow_str)

        self.assertEqual(first_list, ['0.0 1.1 3.0 1.0','0.5 1.0 4.0 0.0'])
        self.assertEqual(sec_list, ['0.5 1.0 4.0 0.0', '0.0 1.1 3.0 1.0'])
    
        #check that the unordered seq matches one of the lists
        self.assertTrue((un_list == first_list) or (un_list == sec_list))
        self.assertNotEqual(first_list, sec_list)

    def test_iter(self):
        """FlowgramCollection __iter__ method should yield flows inorder"""
        f = self.Class(['0.0 1.1 3.0 1.0','0.5 1.0 4.0 0.0','1.5 0.0 2.0 1.0'], \
            Names=['a','b','c'])

        for i,b in zip(f,['0.0 1.1 3.0 1.0','0.5 1.0 4.0 0.0',
                          '1.5 0.0 2.0 1.0']):
            self.assertEqual(i,b)

    def test_str(self):
        """FlowgramCollection __str__ should return sff format"""
        a = [Flowgram('0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0', Name='a',
                      header_info = {'Bases':'TACCCCTTGG','Name Length':'14'}),
             Flowgram('1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0', Name = 'b',
              header_info = {'Bases':'TTATTTACCG','Name Length':'14'})]
        f = FlowgramCollection(a, header_info = {'Flow Chars':'TACG'})
        
        self.assertEqual(str(f), """Common Header:\n  Flow Chars:\tTACG\n\n>a\n  Name Length:\t14\nBases:\tTACCCCTTGG\nFlowgram:\t0.5\t1.0\t4.0\t0.0\t1.5\t0.0\t0.0\t2.0\n\n>b\n  Name Length:\t14\nBases:\tTTATTTACCG\nFlowgram:\t1.5\t1.0\t0.0\t0.0\t2.5\t1.0\t2.0\t1.0\n""")        

    def test_len(self):
        """len(FlowgramCollection) returns length of longest sequence"""
        a = [('a','0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0'),
             ('b','1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0'),
             ('c','2.5 0.0 4.0 0.0 0.5 1.0 0.0 1.0')]
        f = FlowgramCollection(a)
        self.assertEqual(len(f), 3)


    def test_writeToFile(self):
        """FlowgramCollection.writeToFile should write in correct format"""
        a = [Flowgram('0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0', Name='a',
                      header_info = {'Bases':'TACCCCTTGG','Name Length':'14'}),
             Flowgram('1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0', Name = 'b',
              header_info = {'Bases':'TTATTTACCG','Name Length':'14'})]
        f = FlowgramCollection(a, header_info = {'Flow Chars':'TACG'})
        fn = mktemp(suffix='.sff')
        f.writeToFile(fn)
        result = open(fn, 'U').read()
        self.assertEqual(result, """Common Header:\n  Flow Chars:\tTACG\n\n>a\n  Name Length:\t14\nBases:\tTACCCCTTGG\nFlowgram:\t0.5\t1.0\t4.0\t0.0\t1.5\t0.0\t0.0\t2.0\n\n>b\n  Name Length:\t14\nBases:\tTTATTTACCG\nFlowgram:\t1.5\t1.0\t0.0\t0.0\t2.5\t1.0\t2.0\t1.0\n""")
        remove(fn)

    def test_createCommonHeader(self):
        """create_commor_header should return lines for sff common header"""
        a = [Flowgram('0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0', Name='a',
                      header_info = {'Bases':'TACCCCTTGG','Name Length':'14'}),
             Flowgram('1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0', Name = 'b',
              header_info = {'Bases':'TTATTTACCG','Name Length':'14'})]
        f = FlowgramCollection(a, header_info = {'Flow Chars':'TACG'})

        self.assertEqual('\n'.join(f.createCommonHeader()),
                         """Common Header:\n  Flow Chars:\tTACG""")
    def test_toFasta(self):
        """FlowgramCollection should return correct FASTA string"""
        f = self.Class( [  '0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0',
                                     '1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0',
                                     '2.5 0.0 4.0 0.0 0.5 1.0 0.0 1.0',
                                     '0.0 1.0 0.0 3.0 1.5 1.0 1.0 2.0'                                  
                                     ], header_info = {'Flow Chars':'TACG'})
        self.assertEqual(f.toFasta(), '>seq_0\nTACCCCTTGG\n>seq_1\nTTATTTACCG\n>seq_2\nTTTCCCCTAG\n>seq_3\nAGGGTTACGG')

        #NOTE THE FOLLOWING SURPRISING BEHAVIOR BECAUSE OF THE TWO-ITEM
        #SEQUENCE RULE:
        aln = self.Class(['0.5 1.0 0.0 0.0','0.0 1.0 1.0 0.0'],
                         header_info = {'Flow Chars':'TACG'})
        self.assertEqual(aln.toFasta(), '>A\nC\n>T\nA')
        
    def test_toPhylip(self):
        """FlowgramCollection should return PHYLIP string format correctly"""
        f = self.Class( [  '0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0',
                                     '1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0',
                                     '2.5 0.0 4.0 0.0 0.5 1.0 0.0 1.0',
                                     '0.0 1.0 0.0 3.0 1.5 1.0 1.0 2.0'                                  
                                     ], header_info = {'Flow Chars':'TACG'})

        phylip_str, id_map =  f.toPhylip()

        self.assertEqual(phylip_str, """4 10\nseq0000001 TACCCCTTGG\nseq0000002 TTATTTACCG\nseq0000003 TTTCCCCTAG\nseq0000004 AGGGTTACGG""")
        self.assertEqual(id_map, {'seq0000004':'seq_3', 'seq0000001':'seq_0', \
            'seq0000003': 'seq_2', 'seq0000002': 'seq_1'})


    def test_toNexus(self):
        """FlowgramCollection should return correct Nexus string format"""
        f = self.Class( [  '0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0',
                                     '1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0',
                                     '2.5 0.0 4.0 0.0 0.5 1.0 0.0 1.0',
                                     '0.0 1.0 0.0 3.0 1.5 1.0 1.0 2.0'                                  
                                     ], header_info = {'Flow Chars':'TACG'})

        expect = '#NEXUS\n\nbegin data;\n    dimensions ntax=4 nchar=10;\n'+\
        '    format datatype=dna interleave=yes missing=? gap=-;\n'+\
        '    matrix\n    seq_1    TTATTTACCG\n    seq_0'+\
        '    TACCCCTTGG\n    seq_3    AGGGTTACGG\n   '+\
        ' seq_2    TTTCCCCTAG\n\n    ;\nend;'
        self.assertEqual(f.toNexus('dna'), expect)


    def test_toSequenceCollection(self):
        """toSequenceCollection should return sequence collection from flows"""
        f = self.Class( [  '0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0',
                                     '1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0',
                                     '2.5 0.0 4.0 0.0 0.5 1.0 0.0 1.0',
                                     '0.0 1.0 0.0 3.0 1.5 1.0 1.0 2.0'                                  
                                     ], header_info = {'Flow Chars':'TACG'})
        s = f.toSequenceCollection()
        assert isinstance(s,SequenceCollection)
        for i,j in zip(s.iterSeqs(),['TACCCCTTGG','TTATTTACCG','TTTCCCCTAG',
                                   'AGGGTTACGG']):
            self.assertEqual(i,j)
            
        a = [Flowgram('0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0', Name='a',
                      header_info = {'Bases':'TACTTGG','Name Length':'14'}),
             Flowgram('1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0', Name = 'b',
              header_info = {'Bases':'TTATTTG','Name Length':'14'})]

        f = self.Class(a)
        s = f.toSequenceCollection(Bases = True)
        assert isinstance(s,SequenceCollection)
        for i,j in zip(s.iterSeqs(),['TACTTGG','TTATTTG']):
            self.assertEqual(i,j)

    def test_addFlows(self):
        """addFlows should return an alignment with the new sequences appended"""
        a = [('s4', '0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0'),
             ('s3', '1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0')]
        b = [('s1','2.5 0.0 4.0 0.0 0.5 1.0 0.0 1.0'),
             ('s2', '0.0 1.0 0.0 3.0 1.5 1.0 1.0 2.0')]
        f1 = self.Class(a, header_info = {'Flow Chars':'TACG'})
        f2 = self.Class(b, header_info = {'Flow Chars':'TACG'})
        self.assertEqual(f1.addFlows(f2).toFasta(),
                self.Class(a+b, header_info = {'Flow Chars':'TACG'}).toFasta())

        
    def test_iterFlows(self):
        """FlowgramCollection iterFlows() method should support reordering"""
        f = self.Class(['0.0 1.1 3.0 1.0','0.5 1.0 4.0 0.0','1.5 0.0 2.0 1.0'], \
            Names=['a','b','c'])
        flows = map(str,list(f.iterFlows()))
        self.assertEqual(flows, ['0.0\t1.1\t3.0\t1.0',
                         '0.5\t1.0\t4.0\t0.0','1.5\t0.0\t2.0\t1.0'])
        flows = list(f.iterFlows(flow_order=['b','a','a']))
        self.assertEqual(map(str,flows), ['0.5\t1.0\t4.0\t0.0',
                                          '0.0\t1.1\t3.0\t1.0',
                                          '0.0\t1.1\t3.0\t1.0'])
        self.assertSameObj(flows[1], flows[2])
        self.assertSameObj(flows[0], f.NamedFlows['b'])

    def test_Items(self):
        """FlowgramCollection Items should iterate over items in specified order."""
        #should work if one row
        self.assertEqual(list(self.one_seq.Items), [0.0, 1.1, 3.0, 1.0])
        #should take order into account
        self.assertEqual(list(self.ordered1.Items),
                         [0.0, 1.1, 3.0, 1.0] + [0.5, 1.0, 4.0, 0.0])
        self.assertEqual(list(self.ordered2.Items),
                         [0.5, 1.0, 4.0, 0.0] + [0.0, 1.1, 3.0, 1.0])

    def test_takeFlows(self):
        """takeFlows should return new FlowgramCollection with selected seqs."""
        f = self.Class(['0.0 1.1 3.0 1.0','0.5 1.0 4.0 0.0','1.5 0.0 2.0 1.0'], \
            Names=['a','b','c'])
        a = f.takeFlows('bc')
        self.assertTrue(isinstance(a, FlowgramCollection))
        self.assertEqual(a, {'b':'0.5 1.0 4.0 0.0','c':'1.5 0.0 2.0 1.0'})
        #should be able to negate
        a = f.takeFlows('bc', negate=True)
        self.assertEqual(a, {'a':'0.0 1.1 3.0 1.0'})


    def test_getFlowIndices(self):
        """FlowgramCollection getSeqIndices should return names of seqs where f(row) is True"""
        f = self.ambiguous
        is_long = lambda x: len(x) > 10
        is_med = lambda x: len(str(x).replace('N','')) > 7 #strips gaps
        is_any = lambda x: len(x) > 0
        self.assertEqual(f.getFlowIndices(is_long,Bases = True), [])
        f.Names = 'cba'
        self.assertEqual(f.getFlowIndices(is_med,Bases = True), ['c','a'])
        f.Names = 'bac'
        self.assertEqual(f.getFlowIndices(is_med,Bases = True), ['a','c'])
        self.assertEqual(f.getFlowIndices(is_any,Bases = True), ['b','a','c'])
        #should be able to negate
        self.assertEqual(f.getFlowIndices(is_med,Bases = True, negate=True),
                         ['b'])
        self.assertEqual(f.getFlowIndices(is_any, Bases = True,negate=True), [])
        

    def test_takeFlowsIf(self):
        """FlowgramCollection takeFlowsIf should return flows where f(row) is True"""
        is_long = lambda x: len(x) > 10
        is_med = lambda x: len(str(x).replace('N','')) > 7
        is_any = lambda x: len(x) > 0
        
        f = self.ambiguous
        self.assertEqual(f.takeFlowsIf(is_long, Bases = True), {})
        self.assertEqual(f.takeFlowsIf(is_med, Bases = True), \
            {'a':'0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0',
            'c':'1.5 1.0 2.0 0.0 1.5 0.0 0.0 2.0'})
        self.assertEqual(f.takeFlowsIf(is_any, Bases = True), f)
        self.assertTrue(isinstance(f.takeFlowsIf(is_med, Bases = True),
                                   FlowgramCollection))
        #should be able to negate
        self.assertEqual(f.takeFlowsIf(is_med,  Bases = True,negate=True), \
            {'b':'0.0 0.0 0.0 0.0 2.0 1.0 2.0 2.0'})

    def test_getFlow(self):
        """FlowgramCollection.getFlow should return specified flow"""
        a = [('a','0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0'),
             ('b','1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0'),
             ('c','2.5 0.0 4.0 0.0 0.5 1.0 0.0 1.0')]
        f = FlowgramCollection(a)
        self.assertEqual(f.getFlow('a'), '0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0')
        self.assertRaises(KeyError, f.getFlow, 'd')


    def test_getIntMap(self):
        """FlowgramCollection.getIntMap should return correct mapping."""
        f = self.Class({'seq1':'0.5 1.0 2.0 0.0',
                          'seq2':'1.5 0.0 0.0 2.0','seq3':'0.0 3.0 0.1 1.0'})
        int_keys = {'seq_0':'seq1','seq_1':'seq2','seq_2':'seq3'}
        int_map = {'seq_0':'0.5 1.0 2.0 0.0','seq_1':'1.5 0.0 0.0 2.0',
                   'seq_2':'0.0 3.0 0.1 1.0'}
        im,ik = f.getIntMap()
        self.assertEqual(ik,int_keys)
        self.assertEqual(im,int_map)
        #test change prefix from default 'seq_'
        prefix='seqn_'
        int_keys = {'seqn_0':'seq1','seqn_1':'seq2','seqn_2':'seq3'}
        int_map = {'seqn_0':'0.5 1.0 2.0 0.0','seqn_1':'1.5 0.0 0.0 2.0',
                   'seqn_2':'0.0 3.0 0.1 1.0'}
        im,ik = f.getIntMap(prefix=prefix)
        self.assertEqual(ik,int_keys)
        self.assertEqual(im,int_map)


    def test_toDict(self):
        """FlowgramCollection.toDict should return dict of strings (not obj)"""
        f = self.Class({'a': '0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0'
                          , 'b': '1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0'})
        self.assertEqual(f.toDict(), {'a':'0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0'
                                      ,'b':'1.5 1.0 0.0 0.0 2.5 1.0 2.0 1.0'})
        for i in f.toDict().values():
            assert isinstance(i, str)
            
    def test_omitAmbiguousFlows(self):
        """FlowgramCollection omitAmbiguousFlows should return flows w/o N's"""
        
        self.assertEqual(self.ambiguous.omitAmbiguousFlows(Bases=True),
                         {'a':'0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0',
                          'c':'1.5 1.0 2.0 0.0 1.5 0.0 0.0 2.0'})

        self.assertEqual(self.ambiguous.omitAmbiguousFlows(Bases=False),
                         {'a':'0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0',
                          'c':'1.5 1.0 2.0 0.0 1.5 0.0 0.0 2.0'})

        #check new object creation
        self.assertNotSameObj(self.ambiguous.omitAmbiguousFlows(),
                              self.ambiguous)
        self.assertTrue(isinstance(self.ambiguous.omitAmbiguousFlows(
                        Bases = True), FlowgramCollection))

    def test_setBases(self):
        """FlowgramCollection setBases should set Bases property correctly"""
        f = self.Class([Flowgram('0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0', Name='a',
                                 header_info = {'Bases':'TACCCCTTGG'}),
                        Flowgram('0.0 1.0 0.0 0.0 2.0 1.0 2.0 2.0', Name='b',
                                 header_info = {'Bases':'ATTACCGG'}),
                        Flowgram('1.5 1.0 2.0 0.0 1.5 0.0 0.0 2.0', Name='c',
                                 header_info = {'Bases':'TTACCTTGG'})],
                       header_info = {'Flow Chars':'TACG'})

        f.setBases()

        for i,b in zip(f,['TACCCCTTGG','ATTACCGG','TTACCTTGG']):
            self.assertEqual(i.Bases,b)
    
    def setUp(self):
        """Define some standard data"""
        self.rec = """Common Header:
  Magic Number:  0x2E736666
  Version:       0001
  Index Offset:  96099976
  Index Length:  1158685
  # of Reads:    57902
  Header Length: 440
  Key Length:    4
  # of Flows:    400
  Flowgram Code: 1
  Flow Chars:    TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
  Key Sequence:  TCAG

>FIQU8OX05GCVRO
  Run Prefix:   R_2008_10_15_16_11_02_
  Region #:     5
  XY Location:  2489_3906

  Run Name:       R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford
  Analysis Name:  /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis
  Full Path:      /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis

  Read Header Len:  32
  Name Length:      14
  # of Bases:       104
  Clip Qual Left:   5
  Clip Qual Right:  85
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.06	0.08	1.04	0.08	0.05	0.94	0.10	2.01	0.10	0.07	0.96	0.09	1.04	1.96	1.07	0.10	1.01	0.13	0.08	1.01	1.06	1.83	2.89	0.18	0.96	0.13	0.99	0.11	1.94	0.12	0.13	1.92	0.21	0.07	0.94	0.17	0.03	0.97	2.76	0.15	0.05	1.02	1.14	0.10	0.98	2.54	1.13	0.96	0.15	0.21	1.90	0.16	0.07	1.78	0.22	0.07	0.93	0.22	0.97	0.08	2.02	0.15	0.19	1.02	0.19	0.09	1.02	0.17	0.99	0.09	0.18	1.84	0.16	0.91	0.10	1.10	1.00	0.20	0.09	1.11	3.01	1.07	1.98	0.14	0.22	1.09	0.17	1.99	0.15	0.20	0.92	0.17	0.07	1.01	2.96	0.15	0.07	1.06	0.20	1.00	0.10	0.12	1.00	0.15	0.08	1.90	0.19	0.10	0.99	0.18	0.09	0.99	1.08	0.15	0.07	1.06	0.14	1.84	0.13	0.11	0.95	1.05	0.13	1.04	1.10	0.18	0.94	0.14	0.10	0.97	1.08	0.12	1.08	0.18	0.08	1.00	0.13	0.98	0.15	0.87	0.13	0.19	1.01	3.06	0.17	0.11	1.04	0.09	1.03	0.10	0.11	2.02	0.16	0.11	1.04	0.04	0.09	1.87	0.13	2.09	0.13	0.10	0.97	0.17	0.08	0.08	0.04	0.12	0.05	0.08	0.07	0.08	0.05	0.07	0.06	0.07	0.03	0.05	0.04	0.09	0.04	0.07	0.04	0.07	0.06	0.03	0.06	0.06	0.06	0.06	0.07	0.09	0.04	0.05	0.08	0.05	0.04	0.09	0.06	0.03	0.02	0.08	0.04	0.06	0.05	0.08	0.03	0.08	0.05	0.05	0.05	0.10	0.05	0.05	0.07	0.06	0.04	0.06	0.05	0.03	0.04	0.05	0.06	0.04	0.04	0.07	0.04	0.04	0.05	0.05	0.04	0.07	0.06	0.05	0.03	0.08	0.05	0.06	0.04	0.06	0.05	0.04	0.04	0.04	0.05	0.06	0.04	0.05	0.04	0.05	0.05	0.06	0.05	0.06	0.04	0.06	0.07	0.06	0.05	0.05	0.05	0.06	0.06	0.04	0.05	0.06	0.03	0.06	0.04	0.06	0.05	0.03	0.06	0.06	0.05	0.06	0.04	0.03	0.06	0.06	0.06	0.03	0.04	0.05	0.05	0.07	0.04	0.05	0.06	0.07	0.07	0.05	0.07	0.06	0.05	0.06	0.05	0.07	0.06	0.05	0.06	0.07	0.05	0.06	0.04	0.06	0.05	0.05	0.06	0.04	0.06	0.04	0.03	0.06	0.05	0.05	0.04	0.05	0.05	0.04	0.04	0.05	0.06	0.06	0.04	0.04	0.05	0.06	0.04	0.04	0.04	0.05	0.05	0.04	0.05	0.05	0.03	0.06	0.06	0.06	0.04	0.07	0.05	0.05	0.04	0.06	0.06	0.05	0.05	0.07	0.04	0.06	0.06	0.06	0.04	0.06	0.03	0.06	0.04	0.06	0.04	0.09	0.05	0.05	0.05	0.07	0.06	0.05	0.05	0.06	0.05	0.05	0.05	0.04	0.04	0.06	0.05	0.05	0.05	0.05	0.04	0.05	0.05	0.06	0.04	0.05	0.05	0.05	0.05	0.05	0.04	0.06	0.04	0.05	0.05	0.04	0.05	0.05	0.05	0.04
Flow Indexes:	1	3	6	8	8	11	13	14	14	15	17	20	21	22	22	23	23	23	25	27	29	29	32	32	35	38	39	39	39	42	43	45	46	46	46	47	48	51	51	54	54	57	59	61	61	64	67	69	72	72	74	76	77	80	81	81	81	82	83	83	86	88	88	91	94	95	95	95	98	100	103	106	106	109	112	113	116	118	118	121	122	124	125	127	130	131	133	136	138	140	143	144	144	144	147	149	152	152	155	158	158	160	160	163
Bases:	tcagGCTAACTGTAACCCTCTTGGCACCCACTAAACGCCAATCTTGCTGGAGTGTTTACCAGGCACCCAGCAATGTGAATAGTCActgagcgggctggcaaggc
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	37	37	37	37	37	39	39	39	39	24	24	24	37	34	28	24	24	24	28	34	39	39	39	39	39	39	39	39	39	39	39	39	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37

>FIQU8OX05F8ILF
  Run Prefix:   R_2008_10_15_16_11_02_
  Region #:     5
  XY Location:  2440_0913

  Run Name:       R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford
  Analysis Name:  /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis
  Full Path:      /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis

  Read Header Len:  32
  Name Length:      14
  # of Bases:       206
  Clip Qual Left:   5
  Clip Qual Right:  187
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.04	0.00	1.01	0.00	0.00	1.00	0.00	1.00	0.00	1.05	0.00	0.91	0.10	1.07	0.95	1.01	0.00	0.06	0.93	0.02	0.03	1.06	1.18	0.09	1.00	0.05	0.90	0.11	0.07	1.99	0.11	0.02	1.96	1.04	0.13	0.01	2.83	0.10	1.97	0.06	0.11	1.04	0.13	0.03	0.98	1.15	0.07	1.00	0.07	0.08	0.98	0.11	1.92	0.05	0.04	2.96	1.02	1.02	0.04	0.93	1.00	0.13	0.04	1.00	1.03	0.08	0.97	0.13	0.11	1.88	0.09	0.05	1.02	1.89	0.07	0.11	0.98	0.05	0.07	1.01	0.08	0.05	1.01	0.13	1.00	0.07	0.10	1.04	0.10	0.04	0.98	0.12	1.03	0.96	0.11	0.07	1.00	0.09	0.03	1.03	0.11	1.95	1.06	0.13	0.05	1.00	0.13	0.11	1.00	0.09	0.03	2.89	0.08	0.95	0.09	1.03	1.02	1.05	1.07	0.08	0.12	2.81	0.08	0.08	1.00	1.07	0.07	0.05	1.86	0.12	0.98	0.06	2.00	0.11	1.02	0.11	0.08	1.88	0.13	1.03	0.13	0.98	0.15	0.11	1.03	1.03	1.04	0.18	0.98	0.13	0.15	1.04	0.11	1.01	0.13	0.06	1.01	0.06	1.02	0.08	0.99	0.14	0.99	0.09	0.05	1.09	0.04	0.07	2.96	0.09	2.03	0.13	2.96	1.13	0.08	1.03	0.07	0.99	0.11	0.05	1.05	1.04	0.09	0.07	1.00	1.03	0.09	0.06	1.06	1.04	2.94	0.18	0.06	0.93	0.10	1.10	0.11	2.02	0.17	1.00	1.03	0.06	0.11	0.96	0.04	3.00	0.11	0.07	1.99	0.10	2.03	0.12	0.97	0.16	0.01	2.09	0.14	1.04	0.16	0.06	1.03	0.14	1.12	0.12	0.05	0.96	1.01	0.10	0.14	0.94	0.03	0.12	1.10	0.92	0.09	1.10	1.04	1.02	0.12	0.97	2.00	0.15	1.08	0.04	1.03	1.04	0.03	0.09	5.16	1.02	0.09	0.13	2.66	0.09	0.05	1.06	0.07	0.89	0.05	0.12	1.10	0.16	0.06	1.01	0.13	1.00	0.14	0.98	0.09	2.92	1.28	0.03	2.95	0.98	0.16	0.08	0.95	0.96	1.09	0.08	1.07	1.01	0.16	0.06	4.52	0.12	1.03	0.07	0.09	1.03	0.14	0.03	1.01	1.99	1.05	0.14	1.03	0.13	0.03	1.10	0.10	0.96	0.11	0.99	0.12	0.05	0.94	2.83	0.14	0.12	0.96	0.00	1.00	0.11	0.14	1.98	0.08	0.11	1.04	0.01	0.11	2.03	0.15	2.05	0.10	0.03	0.93	0.01	0.08	0.12	0.00	0.16	0.05	0.07	0.08	0.11	0.07	0.05	0.04	0.10	0.05	0.05	0.03	0.07	0.03	0.04	0.04	0.06	0.03	0.05	0.04	0.09	0.03	0.08	0.03	0.07	0.02	0.05	0.02	0.06	0.01	0.05	0.04	0.06	0.02	0.04	0.04	0.04	0.03	0.03	0.06	0.06	0.03	0.02	0.02	0.08	0.03	0.01	0.01	0.06	0.03	0.01	0.03	0.04	0.02	0.00	0.02	0.05	0.00	0.02	0.02	0.03	0.00	0.02	0.02	0.04	0.01	0.00	0.01	0.05
Flow Indexes:	1	3	6	8	10	12	14	15	16	19	22	23	25	27	30	30	33	33	34	37	37	37	39	39	42	45	46	48	51	53	53	56	56	56	57	58	60	61	64	65	67	70	70	73	74	74	77	80	83	85	88	91	93	94	97	100	102	102	103	106	109	112	112	112	114	116	117	118	119	122	122	122	125	126	129	129	131	133	133	135	138	138	140	142	145	146	147	149	152	154	157	159	161	163	166	169	169	169	171	171	173	173	173	174	176	178	181	182	185	186	189	190	191	191	191	194	196	198	198	200	201	204	206	206	206	209	209	211	211	213	216	216	218	221	223	226	227	230	233	234	236	237	238	240	241	241	243	245	246	249	249	249	249	249	250	253	253	253	256	258	261	264	266	268	270	270	270	271	273	273	273	274	277	278	279	281	282	285	285	285	285	285	287	290	293	294	294	295	297	300	302	304	307	308	308	308	311	313	316	316	319	322	322	324	324	327
Bases:	tcagAGACGCACTCAATTATTTCCATAGCTTGGGTAGTGTCAATAATGCTGCTATGAACATGGGAGTACAAATATTCTTCAAGATACTGATCTCATTTCCTTTAGATATATACCCAGAAGTGAAATTCCTGGATCACATAGTAGTTCTATTTTTATTTGATGAGAAACTTTATACTATTTTTCATAActgagcgggctggcaaggc
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	34	34	34	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	36	36	36	38	25	25	25	38	37	37	37	37	37	37	33	33	34	37	37	37	37	37	37	37	38	34	20	20	26	26	20	34	38	37	37	37	37	37	37	37	37	37	38	38	38	37	37	37	37	37	37	37	37	37	37

""".split('\n')

        self.flow = """1.06	0.08	1.04	0.08	0.05	0.94	0.10	2.01	0.10	0.07	0.96	0.09	1.04	1.96	1.07	0.10	1.01	0.13	0.08	1.01	1.06	1.83	2.89	0.18	0.96	0.13	0.99	0.11	1.94	0.12	0.13	1.92	0.21	0.07	0.94	0.17	0.03	0.97	2.76	0.15	0.05	1.02	1.14	0.10	0.98	2.54	1.13	0.96	0.15	0.21	1.90	0.16	0.07	1.78	0.22	0.07	0.93	0.22	0.97	0.08	2.02	0.15	0.19	1.02	0.19	0.09	1.02	0.17	0.99	0.09	0.18	1.84	0.16	0.91	0.10	1.10	1.00	0.20	0.09	1.11	3.01	1.07	1.98	0.14	0.22	1.09	0.17	1.99	0.15	0.20	0.92	0.17	0.07	1.01	2.96	0.15	0.07	1.06	0.20	1.00	0.10	0.12	1.00	0.15	0.08	1.90	0.19	0.10	0.99	0.18	0.09	0.99	1.08	0.15	0.07	1.06	0.14	1.84	0.13	0.11	0.95	1.05	0.13	1.04	1.10	0.18	0.94	0.14	0.10	0.97	1.08	0.12	1.08	0.18	0.08	1.00	0.13	0.98	0.15	0.87	0.13	0.19	1.01	3.06	0.17	0.11	1.04	0.09	1.03	0.10	0.11	2.02	0.16	0.11	1.04	0.04	0.09	1.87	0.13	2.09	0.13	0.10	0.97	0.17	0.08	0.08	0.04	0.12	0.05	0.08	0.07	0.08	0.05	0.07	0.06	0.07	0.03	0.05	0.04	0.09	0.04	0.07	0.04	0.07	0.06	0.03	0.06	0.06	0.06	0.06	0.07	0.09	0.04	0.05	0.08	0.05	0.04	0.09	0.06	0.03	0.02	0.08	0.04	0.06	0.05	0.08	0.03	0.08	0.05	0.05	0.05	0.10	0.05	0.05	0.07	0.06	0.04	0.06	0.05	0.03	0.04	0.05	0.06	0.04	0.04	0.07	0.04	0.04	0.05	0.05	0.04	0.07	0.06	0.05	0.03	0.08	0.05	0.06	0.04	0.06	0.05	0.04	0.04	0.04	0.05	0.06	0.04	0.05	0.04	0.05	0.05	0.06	0.05	0.06	0.04	0.06	0.07	0.06	0.05	0.05	0.05	0.06	0.06	0.04	0.05	0.06	0.03	0.06	0.04	0.06	0.05	0.03	0.06	0.06	0.05	0.06	0.04	0.03	0.06	0.06	0.06	0.03	0.04	0.05	0.05	0.07	0.04	0.05	0.06	0.07	0.07	0.05	0.07	0.06	0.05	0.06	0.05	0.07	0.06	0.05	0.06	0.07	0.05	0.06	0.04	0.06	0.05	0.05	0.06	0.04	0.06	0.04	0.03	0.06	0.05	0.05	0.04	0.05	0.05	0.04	0.04	0.05	0.06	0.06	0.04	0.04	0.05	0.06	0.04	0.04	0.04	0.05	0.05	0.04	0.05	0.05	0.03	0.06	0.06	0.06	0.04	0.07	0.05	0.05	0.04	0.06	0.06	0.05	0.05	0.07	0.04	0.06	0.06	0.06	0.04	0.06	0.03	0.06	0.04	0.06	0.04	0.09	0.05	0.05	0.05	0.07	0.06	0.05	0.05	0.06	0.05	0.05	0.05	0.04	0.04	0.06	0.05	0.05	0.05	0.05	0.04	0.05	0.05	0.06	0.04	0.05	0.05	0.05	0.05	0.05	0.04	0.06	0.04	0.05	0.05	0.04	0.05	0.05	0.05	0.04"""

        self.unordered = self.Class({'a':'0.0 1.1 3.0 1.0',
                                     'b':'0.5 1.0 4.0 0.0'})
        self.ordered1 = self.Class({'a':'0.0 1.1 3.0 1.0',\
                                    'b':'0.5 1.0 4.0 0.0'}, Names=['a','b'])
        self.ordered2 = self.Class({'a':'0.0 1.1 3.0 1.0',\
                                    'b':'0.5 1.0 4.0 0.0'}, Names=['b','a'])
        self.one_seq = self.Class({'a':'0.0 1.1 3.0 1.0'})

        self.ambiguous = self.Class([Flowgram('0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0', Name='a',
                                              header_info = {'Bases':'TACCCCTTGG'}),
                                     Flowgram('0.0 0.0 0.0 0.0 2.0 1.0 2.0 2.0', Name = 'b',
                                              header_info = {'Bases':'NTTACCGG'}),
                                     Flowgram('1.5 1.0 2.0 0.0 1.5 0.0 0.0 2.0', Name='c',
                                              header_info = {'Bases':'TTACCTTGG'})],
                                    header_info = {'Flow Chars':'TACG'})
                            
if __name__ == "__main__":
    main()

