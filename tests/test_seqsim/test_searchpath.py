#!/usr/bin/env python
    
"""Tests private methods of SearchPath and SearchNode classes.
"""

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import NonnegIntError
from cogent.seqsim.searchpath import SearchPath, SearchNode

__author__ = "Amanda Birmingham"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Amanda Birmingham"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Amanda Birmingham"
__email__ = "amanda.birmingham@thermofisher.com"
__status__ = "Production"

class SearchPathHelper(object):
    """Contains data and function defs used by public AND private tests"""

    #all primers have certain forbidden sequences: no runs longer than 3 of
    #any purines or pyrimidines
    standard_forbid_seq = ['AAAA', 'GAAA', 'AGAA', 'GGAA', 'AAGA', 'GAGA', \
                            'AGGA', 'GGGA', 'AAAG', 'GAAG', 'AGAG', 'GGAG', \
                            'AAGG', 'GAGG', 'AGGG', 'GGGG', 'CCCC', 'TCCC', \
                            'CTCC', 'TTCC', 'CCTC', 'TCTC', 'CTTC', 'TTTC', \
                            'CCCT', 'TCCT', 'CTCT', 'TTCT', 'CCTT', 'TCTT', \
                            'CTTT', 'TTTT']
                            
    alphabets = {SearchPath.DEFAULT_KEY:"ACGT"}
#end SearchPathHelper

class SearchNodeHelper(object):
    """Contains data and function defs used by public AND private tests"""
    
    #list of possible bases at any position
    alphabet = ["A", "C", "G", "T"]
#end SearchNodeHelper

class SearchPathTests(TestCase):
    """Tests public SearchPath methods."""

    #-------------------------------------------------------------------
    #Tests of clearNodes
    
    def test_clearNodes(self):
        """Should empty path stack and variable forbidden"""
        
        #create a searchpath and add just one node
        test = SearchPath(SearchPathHelper.alphabets)
        test.generate(1)
        
        #now call clear and make sure path value is "" (empty)
        test.clearNodes()
        self.assertEquals(test.Value, "")
    #end test_clearNodes
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    #Tests of generate method
    
    def tripletGenerator(self, pathobj, path_length):
        """Helper function to generate primers and catalogue triplets"""
        
        found_triplets = {}
 
        #make a hundred random searchpaths 
        for i in xrange(100):
            curr_path_val = pathobj.generate(path_length)
            
            #now find all the triplets in this path
            for r in xrange(path_length-2):
                curr_triplet = curr_path_val[r:r+3]
                found_triplets[curr_triplet] = True
            #end if
            
            #clear out the path
            pathobj.clearNodes()
        #next rand sequence
        
        return found_triplets
    #end tripletGenerator
    
    def test_generate_fullCoverage(self):
        """With no constraints, should produce all possible triplets"""
        
        path_length = 20
        test = SearchPath(SearchPathHelper.alphabets)
        
        #make a hundred random searchpaths and see what triplets produced
        found_triplets = self.tripletGenerator(test, path_length)
    
        num_found = len(found_triplets.keys())
        self.assertEquals(num_found, 64)
    #end test_generate_fullCoverage
    
    def test_generate_withForbidden(self):
        """With 2 triplet constraints, should produce all others"""
        
        forbidden_triplet = ["ATG", "CCT"]
        path_length = 20
        test = SearchPath(SearchPathHelper.alphabets, forbidden_triplet)
        
        #make a hundred random searchpaths and see what triplets produced
        found_triplets = self.tripletGenerator(test, path_length)

        num_found = len(found_triplets.keys())
        self.assertEquals(num_found, 62) 
    #end test_generate_oneForbidden
    
    def test_generate_nonePossible(self):
        """Should return null if no path can match constraints"""
        
        alphabet = {SearchPath.DEFAULT_KEY:"AB"}
        #forbid all combinations of alphabet
        forbidden_seqs = ["AA", "AB", "BB", "BA"]
        
        test = SearchPath(alphabet, forbidden_seqs)
        output = test.generate(2)
        
        self.assertEquals(output, None)
    #end test_generate_nonePossible
    
    def test_generate_multiple(self):
        """Should be able to call generate multiple times to extend path"""
        
        test = SearchPath(SearchPathHelper.alphabets)
        output1 = test.generate(2)
        output2 = test.generate(3)
        
        #make sure that the length of the path is now three
        self.assertEquals(len(output2), 3)
        
        #make sure that the new path is a superset of the old one
        self.assertEquals(output1, output2[:2])
    #end test_generate_multiple
    
    def test_generate_correctAlph(self):
        """Should get correct alphabet even if node is popped then readded"""
        
        test_alphs = {0:"A",1:"BC",2:"D",3:"E",SearchPath.DEFAULT_KEY:"X"}
        forbidden_seqs = ["CD"]
        test = SearchPath(test_alphs, forbidden_seqs)

        #given these position alphabets and this forbidden seq,
        #the only legal 3-node searchpath should be ABD.  Make
        #a hundred searchpaths and make sure this is the only one
        #that actually shows up.
        found_paths = {}
        for i in xrange(100):
            curr_path = test.generate(3)
            found_paths[curr_path] = True
            test.clearNodes()
        #next 
        
        #make sure there is only one path found and that it is the right one
        found_path_str = str("".join(found_paths.keys()))
        self.assertEquals(len(found_paths), 1)
        self.assertEquals(found_path_str, "ABD")
    #end test_generate_correctAlph
    #-------------------------------------------------------------------
    
    #-------------------------------------------------------------------
    #Tests of findAllowedOption
    
    def test_findAllowedOption_currentAllowed(self):
        """Should return true when current option is allowed"""
        
        #searchpath with no forbidden seqs, so anything should work
        test = SearchPath(SearchPathHelper.alphabets)
        test._add_node(SearchNode(SearchNodeHelper.alphabet))
        allowed_found = test.findAllowedOption()
        
        self.assertEquals(allowed_found, True)
    #end test_findAllowedOption_currentAllowed
    
    def test_findAllowedOption_otherAllowed(self):
        """Should return true when curr option is bad but another is good"""
        
        node_vals = []
        #create a path and put in 2 nodes; since all the forbidden seqs I 
        #used to init the path have 4 entries, there should be no chance that
        #the path I just created has anything forbidden in it
        test = self._fill_path(2, node_vals)
        #add the existing path value to the forbidden list
        test._fixed_forbidden["".join(node_vals)] = True
        test._forbidden_lengths[2] = True
        
        #call findAllowedOption ... should find next available good option
        allowed_found = test.findAllowedOption()
        self.assertEquals(allowed_found, True)
    #end test_findAllowedOption_otherAllowed  
    
    def test_findAllowedOption_none(self):
        """Should return false if curr option is bad and no good exist"""
        
        test = self._fill_path(1)        
        self._empty_top(test)
        #get the value of the top node's only remaining option;
        #add to forbidden
        last_option = test._get_top().Options
        test._fixed_forbidden["".join(last_option)] = True
        test._forbidden_lengths[1] = True

        #now make sure we get back result that no options for path remain
        allowed_found = test.findAllowedOption()
        self.assertEquals(allowed_found, False)        
    #end test_findAllowedOption_none
    #-------------------------------------------------------------------  

    #-------------------------------------------------------------------  
    #Tests of removeOption
    
    #Helper function
    def _fill_path(self, num_nodes, node_vals = []):
        """create a searchpath and add searchnodes; return path"""
        
        test = SearchPath(SearchPathHelper.alphabets, \
                SearchPathHelper.standard_forbid_seq)
        for i in xrange(num_nodes):
            curr_node = SearchNode(SearchNodeHelper.alphabet)
            node_vals.append(curr_node.Value)
            test._add_node(curr_node)
        #next i  
        
        return test
    #end _fill_path
    
    #Helper function
    def _empty_top(self, spath):
        """remove all but one options from the top node"""
        
        top_node = spath._get_top()
        num_options = len(top_node.Options)
        for i in xrange(num_options-1): top_node.removeOption()    
    #end _empty_top

    def test_removeOption_simple(self):
        """Should correctly remove option from untapped node"""     
        
        #create a searchpath and a searchnode
        test = self._fill_path(1) 
        orig_len_minus1 = len(test._get_top().Options) - 1
        
        #also check that remove result is true: node still has options
        has_options = test.removeOption()
        self.assertEqual(has_options, True)
        
        #get the top node and make sure that it has fewer options
        top_node = test._get_top()
        option_len = len(top_node.Options)
        self.assertEqual(option_len, orig_len_minus1)
    #end test_removeOption_simple
    
    def test_removeOption_empty(self):
        """Should return False if removing option leads to empty stack"""
        
        #create a searchpath with just one searchnode, then (almost) empty it
        test = self._fill_path(1)
        self._empty_top(test)
        
        #now remove the last option, and make sure the stack is now empty
        some_left = test.removeOption()
        self.assertEquals(some_left, False)
    #end test_removeOption_empty
    
    def test_removeOption_recurse(self):
        """Should correctly remove empty node and curr option of next"""
        
        #put two nodes in the search path and almost empty the top one
        test = self._fill_path(2)        
        self._empty_top(test)
        
        test.removeOption()
        #make sure there's only one item left in the path
        self.assertEquals(len(test._path_stack), 1)
        #make sure that it has one fewer options
        top_node = test._get_top()
        self.assertEquals(len(top_node.Options), len(top_node.Alphabet)-1)
    #end test_removeOption_recurse 
    #-------------------------------------------------------------------
#end SearchPathTests

class SearchNodeTests(TestCase):
    """Tests public SearchNode methods."""
    
    #-------------------------------------------------------------------
    #Tests of removeOption method
    
    def test_removeOption_someLeft(self):
        """removeOption should cull options and return T when some left."""
        
        #create a search node and get its current value
        test = SearchNode(SearchNodeHelper.alphabet)
        last_val = test.Value
        
        some_left = test.removeOption()
        
        #new current value must be different from old
        #and return value must be true
        self.assertNotEqual(test.Value, last_val)
        self.assertEqual(some_left, True)
    #end test_removeOption_someLeft
    
    def test_removeOption_noneLeft(self):
        """removeOption should cull options and return F when none left."""
        
        test = SearchNode(SearchNodeHelper.alphabet)
        num_options = len(test.Options)
        
        #removeOption num_options times: that should get 'em all
        for i in xrange(num_options): some_left = test.removeOption()
        
        #return value should be false (no options remain)
        self.assertEqual(some_left, False)        
    #end test_removeOption_noneLeft
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    #Test of options property (and _get_options method)
    
    def test_options(self):
        """Should return a copy of real options"""
        
        test = SearchNode(SearchNodeHelper.alphabet)
        optionsA = test.Options
        del optionsA[0]
        optionsB = test.Options
        self.assertNotEqual(len(optionsA), len(optionsB))
    #end test_options
    #------------------------------------------------------------------- 
#end SearchNodeTests
class SearchPathTests_private(TestCase):
    """Tests for private SearchPath methods."""
    
    #No need to test __str__: just calls toString in general_tools
    #No need to test _accept_option or _remove_accepted_option: they
        #simply pass in the base class implementation
    #No need to test _in_extra_forbidden: just returns False in base class
        #implementation

    #-------------------------------------------------------------------        
    #Helper functions
    
    def _fill_path(self, num_nodes, node_vals = []):
        """create a searchpath and add searchnodes; return path"""
        
        test = SearchPath(SearchPathHelper.alphabets, \
                SearchPathHelper.standard_forbid_seq)
        for i in xrange(num_nodes):
            curr_node = SearchNode(SearchNodeHelper.alphabet)
            node_vals.append(curr_node.Value)
            test._add_node(curr_node)
        #next i  
        
        return test
    #end _fill_path
    
    def _empty_top(self, spath):
        """remove all but one options from the top node"""
        
        top_node = spath._get_top()
        num_options = len(top_node.Options)
        for i in xrange(num_options-1): top_node.removeOption()    
    #end _empty_top
    #-------------------------------------------------------------------
    
    #-------------------------------------------------------------------
    #Tests of __init__
    
    def test_init_noForbid(self):
        """Init should correctly set private properties w/o forbid list"""

        test = SearchPath(SearchPathHelper.alphabets)
        real_result = len(test._fixed_forbidden.keys())
        self.assertEquals(real_result, 0)
    #end test_init_noForbid
    
    def test_init_withForbid(self):
        """Init should correctly set private properties, w/forbid list"""
        user_input = SearchPathHelper.standard_forbid_seq[:]
        user_input.extend(["AUG", "aaaaccuag"])

        test = SearchPath(SearchPathHelper.alphabets, user_input)
        user_input = [i.upper() for i in user_input]
        user_input.sort()
        real_result = test._fixed_forbidden.keys()
        real_result.sort()
        self.assertEquals(str(real_result), str(user_input))        
    #end test_init_withForbid
    
    def test_init_badAlphabets(self):
        """Init should fail if alphabets param is not dictionary-like"""
        
        self.assertRaises(ValueError, SearchPath, "blue")
    #end test_init_badAlphabets
    
    def test_init_noDefault(self):
        """Init should fail if alphabets param has no 'default' key"""
        
        self.assertRaises(ValueError, SearchPath, {12:"A"})
    #end test_init_noDefault  
    #------------------------------------------------------------------- 
    
    #-------------------------------------------------------------------
    #Tests of value property (and _get_value method)
    
    def test_value_empty(self):
        """Should return empty string when path is empty"""
        
        test = SearchPath(SearchPathHelper.alphabets)
        self.assertEquals(test.Value, "")
    #end test_value_empty
    
    def test_value(self):
        """Should return string of node values when nodes exist"""
        
        node_vals = []
        test = self._fill_path(3, node_vals)
        self.assertEquals(test.Value, "".join(node_vals))
    #end test_value
    #-------------------------------------------------------------------  

    #-------------------------------------------------------------------
    #Tests of _top_index property (and _get_top_index method)
    
    def test_top_index(self):
        """Should return index of top node when one exists"""
        
        test = self._fill_path(3)
        top_index = test._top_index
        self.assertEquals(top_index,2)
    #end test_top_index
    
    def test_top_index_None(self):
        """Should return None when stack has no entries"""
        
        test = SearchPath(SearchPathHelper.alphabets)
        top_index = test._top_index
        self.assertEquals(top_index, None)
    #end test_top_index_None
    #-------------------------------------------------------------------    
    
    #-------------------------------------------------------------------     
    #Tests of _get_top method
    
    def test_get_top(self):
        """Should return a reference to top node on stack if there is one"""

        test = SearchPath(SearchPathHelper.alphabets)
        test._add_node(SearchNode(SearchNodeHelper.alphabet))       
        topnode = SearchNode(SearchNodeHelper.alphabet)
        test._add_node(topnode)
        
        resultnode = test._get_top()
        self.assertEquals(resultnode, topnode)
    #end test_get_top
    
    def test_get_top_None(self):
        """Should return None if stack is empty"""
        
        test = SearchPath(SearchPathHelper.alphabets)
        topnode = test._get_top()
        self.assertEquals(topnode, None)
    #end test_get_top_None
    #------------------------------------------------------------------- 
    
    #-------------------------------------------------------------------
    #Test of _get_forbidden_lengths
    
    def test_get_forbidden_lengths(self):
        """get_forbidden_lengths should return dict of forbidden seq lens"""
        
        correct_result = str([3, 4, 9])
        user_input = SearchPathHelper.standard_forbid_seq[:]
        user_input.extend(["AUG", "aaaaccuag"])
        test = SearchPath(SearchPathHelper.alphabets, user_input)

        real_dict = test._get_forbidden_lengths()
        real_list = real_dict.keys()
        real_list.sort()
        real_result = str(real_list)
        self.assertEquals(real_result, correct_result)
    #end test_get_forbidden_lengths
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    #Tests of _add_node
    
    def test_add_node_first(self):
        """add_node should correctly add first node and increase top index."""
        
        test = SearchPath(SearchPathHelper.alphabets, \
                SearchPathHelper.standard_forbid_seq)
        test_node = SearchNode(SearchNodeHelper.alphabet)
        test._add_node(test_node)
        self.assertEquals(len(test._path_stack), 1)
        self.assertEquals(test._top_index, 0)
    #end test_add_node_first
    
    def test_add_node_subsequent(self):
        """add_node should correctly add additional nodes and up top index."""
        
        test = SearchPath(SearchPathHelper.alphabets, \
                SearchPathHelper.standard_forbid_seq)
        test_node = SearchNode(SearchNodeHelper.alphabet)
        test_node2 = SearchNode(SearchNodeHelper.alphabet)
        test._add_node(test_node)
        test._add_node(test_node2)
        self.assertEquals(len(test._path_stack), 2)
        self.assertEquals(test._top_index, 1)        
    #end test_add_node_subsequent
    #-------------------------------------------------------------------
    
    #-------------------------------------------------------------------
    #Tests of _get_nmer
    
    def test_get_nmer(self):
        """get_nmer should return correct nmer for n <= length of stack"""
        
        node_values = []
        n = 4
        
        test_path = SearchPath(SearchPathHelper.alphabets, \
                SearchPathHelper.standard_forbid_seq)
        
        for i in xrange(n+1):
            curr_node = SearchNode(SearchNodeHelper.alphabet)
            test_path._add_node(curr_node)
            node_values.append(curr_node.Value)
        #next
        
        #get a nmer, and get the last n values that were put on stack; 
        #should be the same
        real_result = test_path._get_nmer(n)
        correct_result = "".join(node_values[-n:])
        self.assertEquals(real_result, correct_result)
    #end test_get_nmer
    
    def test_get_nmer_tooLong(self):
        """get_nmer should return None for n > length of stack"""
        
        test_path = SearchPath(SearchPathHelper.alphabets, \
                SearchPathHelper.standard_forbid_seq)
        test_node = SearchNode(SearchNodeHelper.alphabet)
        test_path._add_node(test_node)
        
        #stack is 1 long.  Ask for a 2 mer
        real_result = test_path._get_nmer(2)
        self.assertEquals(real_result, None)        
    #end test_get_nmer_tooLong
    
    def test_get_nmer_len1(self):
        """get_nmer should return correct result for nmer 1 on full stack"""
        
        test_path = SearchPath(SearchPathHelper.alphabets, \
                SearchPathHelper.standard_forbid_seq)
        test_node = SearchNode(SearchNodeHelper.alphabet)
        test_path._add_node(test_node)
        correct_result = test_node.Value
        real_result = test_path._get_nmer(1)
        self.assertEquals(real_result, correct_result)
    #end test_get_nmer_len1
    
    def test_get_nmer_len0(self):
        """get_nmer should return an empty string if n is 0"""
 
        #if n is zero, this should return "" even when stack is empty 
        test_path = SearchPath(SearchPathHelper.alphabets, \
                SearchPathHelper.standard_forbid_seq)   
        real_result = test_path._get_nmer(0)
        self.assertEquals(real_result, "")
    #end test_get_nmer_len0
    
    def test_get_nmer_badArg(self):
        """get_nmer should error if given a non integer-castable n"""
        
        test_path = SearchPath(SearchPathHelper.alphabets, \
                SearchPathHelper.standard_forbid_seq)
        self.assertRaises(NonnegIntError, test_path._get_nmer, "blue")
    #end test_get_nmer_badArg
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    #Tests of _check_forbidden_seqs
    
    def test_check_forbidden_seqs_fixed(self):
        """Should return True if path includes a fixed forbidden seq"""
        
        forbidden_seq = ["G", "U", "A"]     
        user_input = ["".join(forbidden_seq)]
        user_input.extend(SearchPathHelper.standard_forbid_seq)
        
        test = SearchPath(SearchPathHelper.alphabets, user_input)
        test._add_node(SearchNode(SearchNodeHelper.alphabet))
        
        #add more values, and cheat so as to make them something forbidden
        for bad_val in forbidden_seq:
            curr_node = SearchNode(SearchNodeHelper.alphabet)
            curr_node._options[0] = bad_val #torque the node's innards
            test._add_node(curr_node)
        #next bad_val
        
        real_result = test._check_forbidden_seqs()
        self.assertEquals(real_result, True)            
    #end test_check_forbidden_seqs_fixed
    
    def test_check_forbidden_seqs_none(self):
        """Should return False if path includes no forbidden seqs"""
        
        #a seq that isn't in the standard fixed forbidden lib
        allowed_seq = ["C", "U", "A", "T"]        
        
        test = SearchPath(SearchPathHelper.alphabets, \
                SearchPathHelper.standard_forbid_seq)
        test._add_node(SearchNode(SearchNodeHelper.alphabet))

        #add more values, and cheat so as to make them something known
        for known_val in allowed_seq:
            curr_node = SearchNode(SearchNodeHelper.alphabet)
            curr_node._options[0] = known_val #torque the node's innards
            test._add_node(curr_node)
        #next bad_val
        
        real_result = test._check_forbidden_seqs()
        self.assertEquals(real_result, False)            
    #end test_check_forbidden_seqs_fixed
    #-------------------------------------------------------------------   
    
    #------------------------------------------------------------------- 
    #Tests of _get_alphabet
    
    def test_get_alphabet_exists(self):
        """Should return alphabet for position when one exists"""
        
        alph1 = "G"
        alph2 = "ACGT"
        test_alphs = {0:alph1, 2:alph1, SearchPath.DEFAULT_KEY:alph2}
        test = SearchPath(test_alphs)        
        
        real_alph = test._get_alphabet(2)
        self.assertEquals(str(real_alph), alph1)
    #end test_get_alphabet_exists
    
    def test_get_alphabet_default(self):
        """Should return default alphabet if none defined for position"""
        
        #SearchPathHelper.alphabets has only a default entry
        test = SearchPath(SearchPathHelper.alphabets) 
        real_alph = test._get_alphabet(0)
        correct_alph = SearchPathHelper.alphabets[SearchPath.DEFAULT_KEY]
        self.assertEquals(str(real_alph), str(correct_alph))
    #end test_get_alphabet_default
    
    def test_get_alphabet_badPosition(self):
        """Should raise error if input isn't castable to nonneg int"""
        
        test = SearchPath(SearchPathHelper.alphabets) 
        self.assertRaises(NonnegIntError, test._get_alphabet, "blue")
    #end test_get_alphabet_badPosition
    #-------------------------------------------------------------------
#end SearchPathTests_private


class SearchNodeTests_private(TestCase):
    """Tests for private SearchNode methods."""
    
    #No need to test __str__: just calls toString in general_tools
    #No need to test _get_value and the value property: just references
        #an item in an array
    #No need to test _get_alphabet and alphabet property: ibid
        
    #-------------------------------------------------------------------
    #Tests of __init__
    
    def test_init_noArg(self):
        """Init should correctly set private properties w/no arg"""
        correct_result = str(SearchNodeHelper.alphabet)

        test = SearchNode(SearchNodeHelper.alphabet)
        options = test.Options
        options.sort()
        real_result = str(options)
        self.assertEquals(real_result, correct_result)
    #end test_init_noArg    
    #-------------------------------------------------------------------    
#end SearchNodeTests_private

if __name__ == '__main__':
    main()
