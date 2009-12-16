#!/usr/bin/env python

"""tests for Trie and compressed Trie class."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2009, the PyCogent Project" 
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Prototype"


from cogent.util.unit_test import TestCase, main
from cogent.util.trie import Trie, Compressed_Trie, build_prefix_map, build_trie, \
    _build_prefix_map

class TrieTests(TestCase):
    
    def setUp(self):
        self.data = dict({"0": "ab", "1":"abababa", "2":"abab",
                          "3":"baba", "4":"ababaa","5":"a", "6":"abababa",
                          "7":"bab", "8":"babba"})
    
    def test_init(self):
        """Trie init should create an empty trie."""
        
        t = Trie()
        self.assertEqual(t.root.labels, [])
        self.assertEqual(t.root.children, {})
        
    def test_insert_find(self):
        """An added key should be found by find."""

        data = self.data
        
        t = Trie()
        for (label, seq) in data.iteritems():
            t.insert(seq, label)
           
        for  (label, seq) in data.iteritems():
            self.assertEqual(label in t.find(seq), True)
        self.assertEqual(t.find("cacgchagc"), [])
        self.assertEqual(t.find("abababa"), ["1","6"])

    def test_insert_unique(self):
        """insert_unique should insert only unique words."""

        data = self.data
        
        t = Trie()
        for (label, seq) in data.iteritems():
            t._insert_unique(seq, label)
                       
        self.assertEqual(t.find("ab"), [])
        self.assertEqual(t.find("cacgchagc"), [])
        self.assertEqual(t.find("abababa"), ["1"])


    def test_build_prefix_map(self):
        """prefix_map should map prefix strings."""
        
        self.assertEqual(dict(_build_prefix_map(self.data.iteritems())),
                         {'1': ['0', '2', '5', '6'],
                          '8': [],
                          '3': ['7'],
                          '4': []})

class Compressed_Trie_Tests(TestCase):
    
    def setUp(self):
        self.data = dict({"0": "ab", "1":"abababa", "2":"abab",
                          "3":"baba", "4":"ababaa","5":"a", "6":"abababa",
                          "7":"bab", "8":"babba"})

        self.trie = build_trie(self.data.iteritems())
    
    def test_init(self):
        """Trie init should create an empty trie."""
        
        t = Compressed_Trie()
        self.assertEqual(t.root.labels, [])
        self.assertEqual(t.root.children, {})
        self.assertEqual(t.root.key, "")
        
    def test_non_zero(self):
        """__non_zero__ should cehck for any data in the trie."""
        t = Compressed_Trie()
        self.assertEqual(t.__nonzero__(), False)
        self.assertEqual(self.trie.__nonzero__(), True)

    def test_len(self):
        """__len__ should return the number of seqs in the trie."""
        
        self.assertEqual(len(self.trie), 9)
        t = Compressed_Trie()
        self.assertEqual(len(t), 0)

    def test_size(self):
        """size should return the number of nodes in the trie."""
        
        self.assertEqual(self.trie.size(), 10)

        #empty trie contins only root node
        t = Compressed_Trie()
        self.assertEqual(size(t), 1)

    def test_to_string(self):
        """_to_string should create a string representation."""

        string_rep = """
key 
{
\tkey a['5']
\t{
\t\tkey b['0']
\t\t{
\t\t\tkey ab['2']
\t\t\t{
\t\t\t\tkey a
\t\t\t\t{
\t\t\t\t\tkey a['4']
\t\t\t\t}
\t\t\t\t{
\t\t\t\t\tkey ba['1', '6']
\t\t\t\t}
\t\t\t}
\t\t}
\t}
}
{
\tkey bab['7']
\t{
\t\tkey a['3']
\t}
\t{
\t\tkey ba['8']
\t}
}
"""
        self.assertEqual(str(self.trie), string_rep)
        
    def test_insert_find(self):
        """An added key should be found by find."""

        data = self.data
        t = Compressed_Trie()
        for (label, seq) in data.iteritems():
            t.insert(seq, label)
    
        for  (label, seq) in data.iteritems():
            self.assertEqual(label in t.find(seq), True)
        self.assertEqual(t.find("abababa"), ["1","6"])
        self.assertEqual(t.find("cacgchagc"), [])


    def test_prefixMap(self):
        """prefix_map (Compressed_Trie) should map prefix strings."""
        
        self.assertEqual(self.trie.prefixMap(),
                         {'1': ['6', '2', '0', '5'],
                          '8': ['7'],
                          '3': [], 
                          '4': []})
    
    def test_init(self):
        """Trie init should create an empty trie."""
        
        t = Trie()
        self.assertEqual(t.root.labels, [])
        self.assertEqual(t.root.children, {})
        
    def test_insert_find(self):
        """An added key should be found by find."""

        data = self.data
        
        t = Trie()
        for (label, seq) in data.iteritems():
            t.insert(seq, label)
           
        for  (label, seq) in data.iteritems():
            self.assertEqual(label in t.find(seq), True)
        self.assertEqual(t.find("cacgchagc"), [])
        self.assertEqual(t.find("abababa"), ["1","6"])

    def test_insert_unique(self):
        """insert_unique should insert only unique words."""

        data = self.data
        
        t = Trie()
        for (label, seq) in data.iteritems():
            t._insert_unique(seq, label)
                       
        self.assertEqual(t.find("ab"), [])
        self.assertEqual(t.find("cacgchagc"), [])
        self.assertEqual(t.find("abababa"), ["1"])


    def test_build_prefix_map(self):
        """prefix_map should map prefix strings."""
        
        self.assertEqual(dict(_build_prefix_map(self.data.iteritems())),
                         {'1': ['0', '2', '5', '6'],
                          '8': [],
                          '3': ['7'],
                          '4': []})
   
    def test_build_trie(self):
        """Build_trie should build trie from seqs."""
        
        t = build_trie(self.data.iteritems(), Trie)
        self.assertTrue(isinstance(t, Trie))
        for (label, seq) in self.data.iteritems():
            self.assertContains(t.find(seq), label)

        self.assertEqual(t.find(""), [])
        self.assertEqual(t.find("ccc"), [])


class Compressed_Trie_Tests(TestCase):
    
    def setUp(self):
        self.data = dict({"0": "ab", "1":"abababa", "2":"abab",
                          "3":"baba", "4":"ababaa","5":"a", "6":"abababa",
                          "7":"bab", "8":"babba"})

        self.trie = build_trie(self.data.iteritems())
    
    def test_init(self):
        """Trie init should create an empty trie."""
        
        t = Compressed_Trie()
        self.assertEqual(t.root.labels, [])
        self.assertEqual(t.root.children, {})
        self.assertEqual(t.root.key, "")
        
    def test_non_zero(self):
        """__non_zero__ should cehck for any data in the trie."""
        t = Compressed_Trie()
        self.assertEqual(t.__nonzero__(), False)
        self.assertEqual(self.trie.__nonzero__(), True)

    def test_len(self):
        """__len__ should return the number of seqs in the trie."""
        
        self.assertEqual(len(self.trie), 9)

    def test_size(self):
        """size should return the number of nodes in the trie."""
        
        self.assertEqual(self.trie.size(), 10)

    def test_to_string(self):
        """_to_string should create a string representation."""

        string_rep = """
key 
{
\tkey a['5']
\t{
\t\tkey b['0']
\t\t{
\t\t\tkey ab['2']
\t\t\t{
\t\t\t\tkey a
\t\t\t\t{
\t\t\t\t\tkey a['4']
\t\t\t\t}
\t\t\t\t{
\t\t\t\t\tkey ba['1', '6']
\t\t\t\t}
\t\t\t}
\t\t}
\t}
}
{
\tkey bab['7']
\t{
\t\tkey a['3']
\t}
\t{
\t\tkey ba['8']
\t}
}
"""
        self.assertEqual(str(self.trie), string_rep)
        
    def test_insert_find(self):
        """An added key should be found by find."""

        data = self.data
        t = Compressed_Trie()
        for (label, seq) in data.iteritems():
            t.insert(seq, label)
    
        for  (label, seq) in data.iteritems():
            self.assertEqual(label in t.find(seq), True)
        self.assertEqual(t.find("abababa"), ["1","6"])
        self.assertEqual(t.find("cacgchagc"), [])


    def test_prefixMap(self):
        """prefix_map (Compressed_Trie) should map prefix strings."""
        
        self.assertEqual(self.trie.prefixMap(),
                         {'1': ['6', '2', '0', '5'],
                          '8': ['7'],
                          '3': [], 
                          '4': []})

    def test_build_trie(self):
        """Build_trie should build trie from seqs."""
        
        t = build_trie(self.data.iteritems())
        
        self.assertTrue(isinstance(t, Compressed_Trie))
        
        for (label, seq) in self.data.iteritems():
            self.assertContains(t.find(seq), label)

        self.assertEqual(t.find(""), [])
        self.assertEqual(t.find("ccc"), [])


if __name__ == "__main__":
    main()
