#!/usr/bin/env python
"""SearchPath and SearchNode classes for generating random w/rules strings.

SearchPath is a class that can generate strings that follow certain rules
re: what alphabet may appear at each position in the string and whether any
substrings are forbidden.  If it encounters a forbidden substring, it is able
to backtrack the minimum amount necessary to get around that.  The main 
method is "generate", which takes in the length of the string to generate
and can be called multiple times (without clearing state) to grow the path 
incrementally.  Generate returns the first good string of the correct length
that it finds--or, if none is possible given the constraints, it returns
None.

Internally, SearchPath treats each position as a SearchNode, which contains 
a number of randomized options and knows how to remove options that have
been determined to produce forbidden substrings.

Revision History:
09/22/03 Amanda Birmingham: created from generalized elements of what started
    as primer_builder.py
11/04/03 Amanda Birmingham: removed unnecessary increment of position variable
    in SearchPath.generate.  Renamed public properties with leading uppercase
    to match standards.
12/03/03 Amanda Birmingham: renamed _remove_option (in SearchPath) to 
    removeOption to indicate it is now public
12/05/03 Amanda Birmingham: made DEFAULT_KEY property of SearchPath public
12/15/03 Amanda Birmingham: renamed _find_allowed_option (in SearchPath) to
    findAllowedOption to indicate it is now public
12/16/03 Amanda Birmingham: added optional parameter to removeOption to
    allow removal of accepted options from a completed path; see comment 
    in removeOption for more details.  Necessary to fix bug when 
    backtracking from a completed primer.
01/05/04 Amanda Birmingham: Altered init of SearchPath so that the default
    value of forbidden_seqs is defined as None, not []--Rob says that 
    using a default mutable object is a big no-no.
01/20/04 Amanda Birmingham: updated to fit into cvsroot and use packages.
"""

from random import shuffle
from cogent.util.misc import toString, makeNonnegInt

__author__ = "Amanda Birmingham"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Amanda Birmingham"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Amanda Birmingham"
__email__ = "amanda.birmingham@thermofisher.com"
__status__ = "Production"

class SearchPath(object):
    """Represents one path through the a random w/rules search tree."""
    
    DEFAULT_KEY = "default"
    
    def __init__(self, alphabets_by_position, forbidden_seqs = None):
        """Prepopulate with input forbidden values.
        
        alphabets_by_position: a dictionary-like object keyed by position
            in searchpath, specifying the alphabet that should be used at
            that position.  Must include a "default" key and alphabet
        forbidden_seqs: a list containing sequences that may never occur.
            In the case of primers, this should include runs of 4 
            purines or pyrimidines, as well as any user-defined forbidden
            seqs.
        """
        
        #reminder: can't use mutable object as default!
        if forbidden_seqs is None: forbidden_seqs = []
        
        #store the alphabet dictionary and forbidden len for later use
        self._alphabets = alphabets_by_position
        try:
            if self.DEFAULT_KEY not in self._alphabets:
                raise ValueError, "alphabets_by_position param must " + \
                                    "contain a " + self.DEFAULT_KEY + "key"
            #end if
        except TypeError:
            raise ValueError, "alphabets_by_position param must function" + \
                                " as a dictionary"
        #end except
        
        #create the private dict of fixed, unacceptable values
        input_forbidden = [(i.upper(), True) for i in forbidden_seqs]
        self._fixed_forbidden = dict(input_forbidden)
        
        #create a dictionary holding the lengths of fixed forbidden seqs
        self._forbidden_lengths = self._get_forbidden_lengths()
        
        #create the variable_forbidden and path_stack properties
        self.clearNodes()
    #end __init__
    
    def __str__(self):
        """Create a human-readable representation of object."""
        return toString(self)
    #end __str__
    
    def _get_top_index(self):
        """Return the value of top index of path stack"""
        
        result = None
        path_length = len(self._path_stack)
        if path_length > 0: result = path_length - 1
        return result
    #end _get_top_index
    
    _top_index = property(_get_top_index)

    def _get_value(self):
        """Read the curr value of each node, concat, and return string"""
        
        node_vals = []
        for curr_node in self._path_stack:
            node_vals.append(curr_node.Value)
        #next node
        
        return "".join(node_vals)
    #end _get_value
    
    Value = property(_get_value)
    
    def clearNodes(self):
        """Clear the node stack"""
        
        #create an stack to hold the path (most recent selection will go
        #on top)
        self._path_stack = []
    #end clearNodes

    def _get_forbidden_lengths(self):
        """Return a dictionary of lengths of all fixed forbidden sequences"""
        
        lengths = {}
        
        #for each key in self._fixed_forbidden, say we need to check
        #nmers of that length.
        for seq in self._fixed_forbidden: lengths[len(seq)] = True
        
        return lengths
    #end _get_forbidden_lengths
    
    def _get_top(self):
        """Retrieve (but don't pop) the top item on the path stack."""
        result = None
        if self._top_index is not None: 
            result = self._path_stack[self._top_index]
        #end if there are actually any nodes on the stack
        return result
    #end _get_top
    
    def _add_node(self, new_searchnode):
        """Add new searchnode to top of path stack and increase top index.
        
        new_searchnode: a SearchNode object to add to the path stack.
            NOTE that it is the responsibility of the user of this function
            to pass a SearchNode; no checks are made, so GIGO
        """
        
        self._path_stack.append(new_searchnode)
    #end _add_node
    
    def _get_alphabet(self, position):
        """Return the alphabet for input position. If none, return default
        
        position: a nonnegative integer or something castable to it.  
            Positions are assumed to be ZERO based.
        """
        
        position = makeNonnegInt(position)
        
        if position in self._alphabets:
            result = self._alphabets[position]
        else:
            result = self._alphabets[self.DEFAULT_KEY]
        #end if
        
        return result
    #end _get_alphabet
    
    def generate(self, path_length):
        """Generate a valid path of required length and return its value.
        
        Returns None if no path is possible.
        
        path_length: a nonnegative integer or castable to it.  Indicates
            length of desired valid path.
        """
        
        path_length = makeNonnegInt(path_length)
        
        #while length of path stack < path_length
        while len(self._path_stack) < path_length:
            #always have to get the alphabet based on the current 
            #top index plus 1.  This is because, if we have to remove
            #a node because it contributes to a forbidden sequence,
            #we need to make sure that the position used to generate
            #the next node is the same
            position = self._top_index
            if position is None: position = -1
            position += 1
            #get the alphabet for the next node
            curr_alphabet = self._get_alphabet(position)
        
            #make new SearchNode and push it onto path
            new_node = SearchNode(curr_alphabet)
            self._add_node(new_node)
            
            #select the next available allowed option
            option_exists = self.findAllowedOption()
            
            #if no more options, no valid searchpath is possible
            #break out and return None
            if not option_exists: return None
        #end while
        
        #return path as string
        return self.Value
    #end generate
    
    def findAllowedOption(self):
        """Finds, sets next allowed option.  Returns false if none exists."""
        
        #initially, assume that the current option is forbidden
        #(means we always check at least once)
        is_forbidden = True
        
        while is_forbidden == True:
            #check whether the current option of the top node is invalid
            is_forbidden = self._check_forbidden_seqs()
            
            if is_forbidden:
                options_remain = self.removeOption()
                
                #if the path is now empty, break out and return false
                if not options_remain: return False
            #end if current option is forbidden
        #end while current option is forbidden
        
        #we found a good option, so do any necessary bookkeeping
        self._accept_option()
    
        return True
    #end findAllowedOption
    
    def _accept_option(self):
        """Bookkeeping to accept a good option. Default impl does nothing."""
        pass
    #end _accept_option
    
    def _remove_accepted_option(self):
        """Bookkeeping to remove previously accepted. Default does nothing."""
        pass
    #end _remove_accepted_option
    
    def removeOption(self, top_is_accepted = False):
        """Remove the current option and return false if no more exist"""
        
        result = True
        
        #if the top option is accepted, remove its nmer from
        #the accepted option list (before we remove the option)
        #This option is only used when calling removeOption 
        #explicitly from another module after a path has been
        #generated.
        if top_is_accepted: self._remove_accepted_option()
       
        #tell top search node to removeOption
        curr_top = self._get_top()
        still_viable = curr_top.removeOption()
        
        #if that node is now out of options
        if not still_viable:
            #pop it off the path
            self._path_stack.pop()
            
            #remove the previously accepted option
            #(in node under the one we just popped)
            self._remove_accepted_option()
            
            #if stack is now empty
            if len(self._path_stack) == 0:
                result = False
            else:
                #call recursively to remove dead-end option(s)
                result = self.removeOption()
            #end if
        #end if
        
        return result
    #end removeOption
    
    def _check_forbidden_seqs(self):
        """Return t/f for whether current path has any forbidden seqs in it"""
        
        result = False #assume current path has no forbidden nmers
        
        #for every length we need to check (bc there's a fixed forbidden
        #sequence with that length), get the current rightmost string
        #of that length from the path and see if it is in a forbidden list
        for length_to_test in self._forbidden_lengths:
            curr_nmer = self._get_nmer(length_to_test)    

            #check if current nmer is forbidden in fixed or variable
            if (self._fixed_forbidden.__contains__(curr_nmer) or \
                self._in_extra_forbidden(curr_nmer)): 
                result = True
                break
        #next length to check
        
        return result
    #end _check_forbidden_seqs    
    
    def _in_extra_forbidden(self, nmer):
        """Return True if nmer forbidden in anything besides fixed forbid"""
        
        #default implementation *has* no other forbidden dictionariesm so
        #the nmer can't be *in* one
        return False
    #end _in_extra_forbidden
    
    def _get_nmer(self, n):
        """Get the string of the last N bases on the path stack.
        
        n: aninteger or integer-castable value; nonnegative
        
        Returns None if n is greater than the number of searchnodes on the
            path stack.  Otherwise, returns a string of the values of those
            nodes, in the order they were pushed on the stack
        """
        
        nmer = []
        result = None
        n = makeNonnegInt(n)
        
        #only try to get the nmer if there are at least n items on stack
        if n <= len(self._path_stack):
            #for n to 0 -- start with the furthest-back entry
            for temp_length in xrange(n, 0, -1):
                #Note that we add one bc temp length is base 1, while
                #top index isn't.  Ex: a stack with length 5 has top
                #index 4.  If you want the last 4, you want to start
                #at index 1.  Thus, 4 - 4 + 1 = 1
                temp_index = self._top_index - temp_length +1

                #get value of node at that index and add to list
                nmer.append(self._path_stack[temp_index].Value)
            #next intermediate length

            #join list to get current nmer and return
            result = "".join(nmer)
        #end if

        return result
    #end _get_nmer                
#end SearchPath

class SearchNode(object):
    """Represents a single choice (base) in a search path (primer)."""
    
    def __init__(self, an_alphabet):
        """Create a node ready to be used."""
        
        self._alphabet = list(an_alphabet)
        self._options = self.Alphabet #note: this gets a COPY
        shuffle(self._options)
        
        #don't need a current option index variable: current option is
        #always zero; it is the # of available options that changes
    #end __init__
    
    def __str__(self):
        """Create a human-readable representation of object."""
        return toString(self)
    #end __str__ 
    
    def _get_alphabet(self):
        """Get a copy of the class's list of choices"""
        
        return self._alphabet[:]
    #end _get_alphabet
    
    def _get_value(self):
        """Get the value of the current option."""
        
        return self._options[0]
    #end _get_value
    
    def _get_options(self):
        """Get a copy of the list of available options"""
        
        return self._options[:]
    #end _get_options
    
    Value = property(_get_value)
    Options = property(_get_options)
    Alphabet = property(_get_alphabet)
    
    def removeOption(self):
        """Remove current option and return t/f whether any more are left."""
        
        #assume, by default, that we won't run out of options,
        #then remove the current option from the options array
        result = True
        del self._options[0]

        #if there are no more options return false; otherwise, true
        if len(self._options) == 0: result = False
        return result
    #end removeOption
#end SearchNode
