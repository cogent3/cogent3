#!/usr/bin/env python
"""Code for handling RNA secondary structure.

RNA secondary structures can be represented in many different ways. The
representation makes a large difference to the efficiency of different
algorithms, so several different structural representations (and the means
to interconvert them) are provided.

Provides following classes:  
    Stem: representation of a stem in a secondary structure
    Partners: list holding partner of each position.
    Pairs: list of base pairs in a structure.
    StructureString: string holding a secondary structure representation.
    ViennaStructure: representation of Vienna-format RNA structure.
    WussStructure: Wash U secondary structure format, handles pseudoknots.
    StructureNode: for tree representation of nested RNA structure.
"""
from numpy import zeros, put
from cogent.util.transform import make_trans, float_from_string
from cogent.core.tree import TreeNode, TreeError
from cogent.util.misc import not_none, flatten

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Sandra Smit", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class PairError(ValueError):
    """Base class for errors in pairing."""
    pass

class Stem(object):
    """Holds a Start, an End, and a Length.
    
    Note that the Start and the End pair with each other, thus spanning the 
    range. Successive pairs count up from the Start and down from the End.

    Stem is intended to be _very_ lightweight and does no error checking. In
    other words, it will let you specify a Start that's after the End, a Length
    that can't exist because there's insufficient separation between the Start
    and the End, and so on. The same principle applies to __getitem__ (and hence
    __iter__): you can iterate through pairs in a stem that can't exist, e.g.
    Stem(6, 7, 10) will give you the 10 items from Stem(6,7,1) through 
    Stem(15, -2, 1) -- often not what you want. 
    
    Similarly, when you initialize a Stem it detects whether Start and End are 
    both set, and makes a pair accordingly, but updating the Start or End will 
    not re-do the check to see if a pair has formed or broken. The alternative
    is to make Start, End and Length into properties that reset the state of
    the others when updated, but this makes Stem too slow for applications such
    as BayesFold.
    """
    __slots__ = ['Start', 'End', 'Length']
    def __init__(self, Start=None, End=None, Length=0):
        """Returns a new Stem object."""
        self.Start = Start
        self.End = End
        #set length if specified
        if Length:
            self.Length = Length
        #otherwise, set to 1 if paired and 0 if unpaired
        elif Start is None or End is None:
            self.Length = 0
        else:
            self.Length = 1

    def __len__(self):
        """Returns self.Length."""
        return self.Length

    def __getitem__(self, item):
        """Masquerades as a list of single-base stems."""
        length = self.Length
        #bounds check
        if item < 0:
            item = length - item
        if (item < 0) or (item >= length):
            raise IndexError, "Index %s out of range." % item
        #return appropriate base pair
        return Stem(self.Start + item, self.End - item)
        
    def __cmp__(self, other):
        """Sorts by start, then by end, then by length (if possible)."""
        return cmp(self.Start, other.Start) or cmp(self.End, other.End) \
            or cmp(self.Length, other.Length)
                
    def extract(self, seq):
        """Returns bases in pairs as list of tuples.
        
        Note: always returns list, even if only one base pair.
        """
        if self.Length > 1:
            return flatten([p.extract(seq) for p in self])
        else:
            if self.Start is not None:
                start = seq[self.Start]
            else:
                start = None
            if self.End is not None:
                end = seq[self.End]
            else:
                end = None
            return [(start, end)]
    
    def __hash__(self):
        """Hashes the same as a tuple of (start, end, length).
        
        WARNING: if you change any of the values once an object is in a
        dict, you'll get unpleasant results. Don't do it!
        """
        return hash((self.Start, self.End, self.Length))

    def __str__(self):
        """String representation contains Start,End,Length."""
        return '(%s,%s,%s)' % (self.Start, self.End, self.Length)
    
    def __nonzero__(self):
        """Nonzero if length > 0."""
        return self.Length > 0

class Partners(list):
    """Holds list p such that p[i] is the index of the partner of i, or None.
    
    Primarily useful for testing whether a specified base is paired and, if so,
    extracting its partner.

    Each base may have precisely 0 or 1 partners. If A pairs with B, B must
    pair with A.
    All inconsistencies will be removed by setting previous partners to None.
    Checking for conflicts and raising errors should be done in method that 
    constructs the Partners.

    If constructing by hand, should initialize with list of [None] * seq_length.
    Typically, Partners will be constructed by code from some other data. Use
    the EmptyPartners(n) factory function to get an empty Partners list of
    length n.
    """
    def __setitem__(self, index, item):
        """Sets self[index] to item, enforcing integrity constraints."""
        if index == item:
            raise ValueError, "Cannot set base %s to pair with itself." % item
        #if item already paired, raise Error or make partner unpaired
        if item and self[item]:
            self[self[item]] = None
        #if already paired, need to make partner unpaired
        curr_partner = self[index]
        if curr_partner is not None:
            list.__setitem__(self, curr_partner, None)
        #set self[index] to item    
        list.__setitem__(self, index, item)
        #if item is not None, set self[item] to index
        if item is not None:
            list.__setitem__(self, item, index)
                
    def toPairs(self):
        """Converts the partners to sorted list of pairs."""
        result = Pairs()
        for first, second in enumerate(self):
            if first < second:
                result.append((first, second))
        return result

    def _not_implemented(self, *args, **kwargs):
        """Raises NotImplementedError for 'naughty' methods.
        
        Not allowed any methods that insert/remove items or that change the
        order of the items, including things like sort or reverse.
        """
        raise NotImplementedError
    
    __delitem__ = __delslice__ = __iadd__ = __imul__ = __setslice__ = append \
    = extend = insert = pop = remove = reverse = sort = _not_implemented

def EmptyPartners(length):
    """Returns empty list of Partners with specified length."""
    return Partners([None] * length)


class Pairs(list):
    """Holds list of base pairs, each of which is a 2-element sequence.
    
    This is a very lightweight object for storing base pairs, and does not
    perform any validation. Useful as an intermediate in many different
    calculations.
    """
    def toPartners(self, length, offset=0, strict=True):
        """Returns a Partners object, if possible.
        
        length of resulting sequence must be specified.
        offset is optional, and is added to each index.
        strict specifies whether collisions cause fatal errors. if not strict
        conflicts will be removed by the Partners object.
        """
        result = EmptyPartners(length)
        for up, down in self:
            upstream = up + offset
            downstream = down + offset
            
            if result[upstream] or result[downstream]:
                if strict:
                    raise ValueError, "Pairs contain conflicting partners: %s"\
                        % self
            result[upstream] = downstream
        return result
            
    def toVienna(self, length, offset=0, strict=True):
        """Returns a Vienna structure string, if possible.

        length of resulting sequence must be specified.
            Instead of parsing the sequence length, you can also throw in
            an object that has the required length (such as the sequence that
            the structure corresponds to).
        offset is optional, and is added to each index.
        strict specifies whether collisions cause fatal errors.
        """
        if self.hasPseudoknots():
            raise PairError, "Pairs contains pseudoknots %s"%(self)
        try:
            length = int(length)
        except ValueError: #raised when length can't be converted to int
            length = len(length)
        
        p = self.directed()
        result = ['.'] * length
        for up, down in p:
            try:
                upstream = up + offset
                downstream = down + offset
            except TypeError:
                continue

            if strict:
                if (result[upstream] != '.') or (result[downstream] != '.'):
                    raise ValueError, "Pairs contain conflicting partners: %s"\
                        % self
            result[upstream] = '('
            result[downstream] = ')'
        return ViennaStructure(''.join(result))

    def tuples(self):
        """Converts all pairs in self to tuples, in place.

        Useful for constructing dicts and for sorting (otherwise, pairs of
        different types, e.g. lists and tuples, will sort according to type
        rather than to position).
        """
        self[:] = map(tuple, self)

    def unique(self):
        """Returns copy of self omitting duplicate pairs, preserving order.

        Keeps the first occurrence of each pair.
        """
        seen = {}
        result = []
        for p in map(tuple, self):
            if p not in seen:
                seen[p] = True
                result.append(p)
        return Pairs(result)

    def directed(self):
        """Returns copy of self where all pairs are (upstream, downstream).
        
        Omits any unpaired bases and any duplicates. Result is in arbitrary 
        order.
        """
        seen = {}
        for up, down in self:
            if (up is None) or (down is None):
                continue    #omit unpaired bases
            if up > down:
                up, down = down, up
            seen[(up, down)] = True
        result = seen.keys()
        return Pairs(result)

    def symmetric(self):
        """Retruns copy of self where  each up, down pair has a down, up pair.
        
        Result is in arbitrary order. Double pairs and pairs containing None 
        are left out.
        """
        result = self.directed()
        result.extend([(down, up) for up, down in result])
        return Pairs(result)

    def paired(self):
        """Returns copy of self omitting items where a 'partner' is None."""
        return Pairs(filter(not_none, self))

        
    def hasPseudoknots(self):
        """Returns True if the pair list contains pseudoknots.
        
        (good_up,good_down) <=> (checked_up,checked_down)
        pseudoknot if checked_up<good_down and checked_down>good_down
        """
        pairs = self.directed()
        seen = [] # list of pairs against which you compare each time
        pairs.sort()
        for pair in pairs:
            if not seen:
                seen.append(pair)
            else:
                lastseen_up, lastseen_down = seen[-1]
                while pair[0] > lastseen_down:
                    seen.pop()
                    if not seen:
                        break
                    else:
                        lastseen_up,lastseen_down = seen[-1]
                if not seen:
                    seen.append(pair)
                    continue
                if pair[1]>lastseen_down:
                    #pseudoknot found
                    return True
                else:
                    #good pair
                    seen.append(pair)
        return False


    def hasConflicts(self):
        """Returns True if the pair list contains conflicts.

        Conflict occurs when a base has two different partners, or is asserted
        to be both paired and unpaired.
        """
        partners = {}
        for first, second in self:
            if first is None:
                if second is None:
                    continue    #no pairing info
                else:
                    first, second = second, first   #swap order so None is 2nd
            if second is None: #check first isn't paired
                if partners.get(first, None) is not None:
                    return True
                else:
                    partners[first] = None
            else:   #first and second were both non-empty: check partners
                if first in partners:
                    if partners[first] != second:
                        return True
                if second in partners:
                    if partners[second] != first:
                        return True
                #add current pair to the list of constraints
                partners[first] = second
                partners[second] = first
        #can only get here if there weren't conflicts
        return False
            
    def mismatches(self, sequence, pairs=None):
        """Counts pairs that can't form in sequence.

        Sequence must have a Pairs property that acts like a dictionary
        containing a 2-element tuple for each valid pair. Can also pass in
        the pairs explicitly.
        """
        mismatches = 0
        if pairs is None:
            try:
                pairs = sequence.Alphabet.Pairs
            except AttributeError:
                pairs = sequence.Pairs
            
        for up, down in self.directed():
            curr = (sequence[up], sequence[down])
            if curr not in pairs:
                mismatches += 1
        return mismatches
    
        
wuss_to_vienna_table = make_trans('<([{>)]} ', '(((()))) ', '.')

def wuss_to_vienna(data):
    """Converts WUSS format string to Vienna format.

    Any pseudoknots or unrecognized chars will convert to unpaired bases.
    Spaces will be preserved.
    """
    return ViennaStructure(data.translate(wuss_to_vienna_table))


class StructureString(str):
    """Base class for ViennaStructure and WussStructure. Immutable.
    
    StructureString holds a structure and a energy. By default energy is 
    set to None. 
    If you compare two StructureStrings the structure is the only important
    thing, since the energy is relative to the associated sequence.
    """
    Alphabet=None
    StartSymbols = ''      #dict of symbols that start base pairs
    EndSymbols = ''        #dict of symbols that end base pairs
    
    def __new__(cls, Structure, Energy=None):
        """Returns new StructureString."""
        a = cls.Alphabet
        if a:
            for i in Structure:
                if i not in a:
                    raise ValueError,\
                    "Tried to include unknown symbol '%s'" % i
        
        return str.__new__(cls,Structure)

    def __init__(self, Structure='', Energy=None):
        """Initializes StructureString with Structure and optionally Energy."""
        self.Energy = Energy
        self.toPartners()

    def __str__(self):
        """Returns string representaion of structure and energy, if known.

        Energy = 0 is different from Energy = None. Latter case is not printed.
        """
        if not self.Energy == None:
            return self + ' (' + str(self.Energy) + ')'
        else:
            return str.__str__(self)

    def toPartners(self):
        """Makes list containing partner of each position.
        
        Constructs a list from 0 to the number of bases, where each position
        contains the index of its pair (or None if it is unpaired).

        Note that the numbering starts at 0 for the first position.

        The algorithm here relies on the fact that any completely nested
        base-paired structure (no pseudoknots!) can be formally represented
        as a tree. Consequently, when you hit a closed base pair, you know
        that it it must pair with the last base pair you opened.
        """
        num_bases = len(self) #number of bases
        result = [None] * len(self) #array of None, one for each base
        stack = []
        start = self.StartSymbols
        end = self.EndSymbols
        for i, symbol in enumerate(self):
           if symbol in start:       #open a pair
               stack.append(i)
           elif symbol in end:     #close a pair
               curr = stack.pop()  #return and delete last element
               result[i] = curr #make i pair with the last element...
               result[curr] = i #...and the last element pair with i
               
        #test whether there are any open pairs left unaccounted for        
        if stack:
           raise IndexError, \
           "Too many open pairs in structure:\n%s" % self
        return Partners(result)

    def toPairs(self):
        """Makes list of (upstream,downstream) partners.
        
        Note that the numbering starts at 0 for the first position.
        Key will always be smaller than value.

        Result is in arbitrary order.
        """
        result = {}
        stack = []
        start = self.StartSymbols
        end = self.EndSymbols
        for i, symbol in enumerate(self):
           if symbol in start:       #open a pair
               stack.append(i)
           elif symbol in end:     #close a pair
               result[stack.pop()] = i
        #test whether there are any open pairs left unaccounted for        
        if stack:
           raise IndexError, \
           "Too many open pairs in structure:\n%s" % self
        return Pairs([(key,result[key]) for key in result])

    def toTree(self):
        """Returns tree version of structure.
        
        Each node in the tree corresponds to a loop. Runs of nodes with single
        non-leaf children correspond to stems.
        """
        root = StructureNode()
        curr_node = root
        start = self.StartSymbols
        end = self.EndSymbols
        for index, symbol in enumerate(self):
            if symbol in end:
                curr_node.End = index
                curr_node.Length = 1
                curr_node = curr_node.Parent
            else:
                new_node = StructureNode()
                new_node.Start = index
                curr_node.Children.append(new_node)
                new_node.Parent = curr_node
                if symbol in start:
                    curr_node = new_node
        return root

class ViennaStructure(StructureString):
    """Contains a Vienna dot-bracket structure, possibly with energy."""
    Alphabet = dict.fromkeys('(.)')
    StartSymbols = {'(':None}      #dict of symbols that start base pairs
    EndSymbols =   {')':None}      #dict of symbols that end base pairs
 
class WussStructure(StructureString):
    """Contains a Wuss Structure."""
    Alphabet = dict.fromkeys(
        '(<{[.~-,:_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ)>}]')
    StartSymbols = dict.fromkeys('(<{[')#dict of symbols that start base pairs
    EndSymbols = dict.fromkeys(')>}]')   #dict of symbols that end base pairs
 
def Vienna(data,Energy=None):
    """Tries to extract structure and energy from string data.

    Returns (structure, energy) where energy might be None.

    structure is just anything before the first space: doesn't validate.
    """
    pieces = data.strip().split(None, 1)
    if not pieces:
        return ViennaStructure('', Energy)
    else:
        if not Energy:
            try:
                energy = float_from_string(pieces[1])
            except (TypeError, ValueError, IndexError):
                energy = Energy
        else: #energy given by user overrules the one in structure
            energy = Energy
    return ViennaStructure(pieces[0], energy)

class StructureNode(TreeNode):
    """Operations on StructureNodes, trees, and subtrees."""
    Types = ["Stem", "Loop", "Bulge", "Junction", "End", "Flexible"]

    def __init__(self, Data=None, Children=None, Parent=None):
        """Returns a new StructureNode object."""
        #coerce data into correct type
        if Data is None:
            Data = Stem()
        elif not isinstance(Data, Stem):
            Data = Stem(Data)
        #initalize as TreeNode: will forward attributes to Data
        super(StructureNode, self).__init__(None, Children, Parent)
        self.Data = Data

    def _get_stems(self):
        """Returns nodes corresponding to stems leaving self."""
        return [i for i in self if i.IsPaired]
    
    def _get_unpaired(self):
        """Returns nodes corresponding to unpaired bases leaving self."""
        return [i for i in self if not i.IsPaired]
    
    Stems = property(_get_stems)
    Unpaired = property(_get_unpaired)

    def _get_end(self):
        """Gets End from self's data."""
        return self.Data.End

    def _set_end(self, val):
        """Sets End in self's data."""
        self.Data.End = val

    def _get_start(self):
        """Gets Start from self's data."""
        return self.Data.Start
    
    def _set_start(self, val):
        """Sets Start in self's data."""
        self.Data.Start = val

    def _get_length(self):
        """Gets Length from self's data."""
        return self.Data.Length

    def _set_length(self, val):
        """Sets Length in self's data."""
        self.Data.Length = val

    def _get_index(self):
        """Accessor for index: returns position of self in parent's list."""
        if self.Parent is None:
            return None
        else:
            ids = map(id, self.Parent.Children)
            return ids.index(id(self))

    def _set_index(self, index):
        """Mutator for index: moves self to new location in parent's list.
        
        NOTE: index is relative to the new list after self is removed, not to
        the old list before self is removed.
        """
        if self.Parent is None:
            raise TreeError, "Can't set Index in node %s without parent." % self
        else:
            curr_parent = self.Parent
            curr_parent.removeNode(self)
            curr_parent.Children.insert(index, self)

    End = property(_get_end, _set_end)
    Start = property(_get_start, _set_start)
    Length = property(_get_length, _set_length)
    Index = property(_get_index, _set_index)
    

    def _is_paired(self):
        """Returns True if self is paired, false otherwise."""
        return self.End is not None

    IsPaired = property(_is_paired)

    def _get_type(self):
        """Returns type of the node, depending on stems in and out."""
        #easy cases first
        if self.Parent is None:
            return "Root"
        if self.IsPaired:
            return "Stem"
        #if grandparent is None (i.e. it is a child of the root node), either 
        #end or 'flexible' depending on whether the current node is outside
        #the first and last helices or between them.
        if self.Parent.Parent is None:
            stems = self.Parent.Stems
            if not stems:
                return "End"
            first, last = stems[0], stems[-1]
            index = self.Index
            if index < stems[0].Index or index > stems[-1].Index:
                return "End"
            else:
                return "Flexible"
        #otherwise, depends on number of stems coming out of parent: if none,
        #it's a loop, if one, it's a bulge, and otherwise it's a junction.
        else:
            out_stems = self.Parent.Stems
            if not out_stems:
                return "Loop"
            elif len(out_stems) == 1:
                return "Bulge"
            else:
                return "Junction"
        #should never get down here, since all cases are handled above
        raise ValueError, '_get_type failed on node with start %s, end %s' % \
            (self.Start, self.End)
    Type = property(_get_type)
        
    def renumber(self, start=0):
        """Renumbers self and all child nodes consecutively.
        
        Returns the next number to be used.
        """
        if self.Parent is None: #no number for root node
            curr = start
        else:
            self.Start = start
            curr = start + max(1, self.Length)  #still add 1 if it's unpaired
        for i in self:
            curr = i.renumber(curr)
        if self.IsPaired:
            curr += self.Length
            self.End = curr - 1
        return curr
            
    def classify(self, terminate=True):
        """Returns string containing site classification"""
        #check whether we're at the original node or in an internal node
        if not terminate:
            #if it's paired, need to add 'S' for self, plus handle children
            if self.IsPaired:
                result = ['S']
                for i in self:
                    result.extend(i.classify(False))
                result.append('S')
                return result
            #otherwise, it's unpaired and we need to figure out what it is
            else:
                return self.Type[0] #just want first letter
        else:   #if it's the root, just need to handle children
            result = []
            for i in self:
                result.extend(i.classify(False))
            return ''.join(result)

    def __str__(self):
        """Returns string representation of tree, in Vienna format."""
        if self.IsPaired:
            prefix = '(' * self.Length
            suffix = ')' * self.Length
            return prefix + ''.join(map(str, self)) + suffix
        elif self.Parent is None:   #root node
            return ''.join(map(str, self))
        else:                       #unpaired base
            return '.'

    def unpair(self):
        """Breaks the first pair represented by the current node, if any.

        Returns True if the node is changed, False if it wasn't.
        """
        if self.IsPaired:
            curr_idx = self.Index
            first = StructureNode(Data=Stem(self.Start))
            last = StructureNode(Data=Stem(self.End))
            
            if self.Length > 1: #not melting the whole helix
                self.Start += 1
                self.End -= 1
                self.Length -= 1
                result = [first, self, last]
            else:   #melting the whole helix
                result = [first] + self.Children + [last]
            #replace current record in parent with the result
            #note use of a slice assignment instead of an index! This is to
            #replace with the elements, not with a list of the elements.
            self.Parent[curr_idx:curr_idx+1] = result
            return True
        else:
            return False

    def _pair_before_indices(self):
        """Detects the pair of positions before the current node. 
        
        Returns a tuple of (upstream, downstream) if possible; None
        otherwise.
        """
        curr_idx = self.Index
        curr_parent = self.Parent
        #bail out if it's the first or last base
        if (curr_idx == 0) or (curr_idx == len(curr_parent) - 1):
            return None
        #bail out if the base before and after self aren't unpaired
        before = curr_parent[curr_idx - 1]
        after = curr_parent[curr_idx + 1]
        if before.IsPaired or after.IsPaired:
            return None
        return (before.Start, after.Start)

    def pairBefore(self):
        """Forms a pair before the current node, if possible.
        
        Returns True if the pair was successfully created, False otherwise.
        """
        #get the indices of the pair, if possible
        indices = self._pair_before_indices()
        if not indices:
            return False
        #cache the current index and parent
        curr_idx = self.Index
        curr_parent = self.Parent
        #create a new pair, and add self to it
        new_pair = StructureNode(Data=Stem(indices[0], indices[1], 1))
        new_pair.Children.append(self)
        self.Parent = new_pair
        #add the new pair to parent, and delete its unpaired siblings
        curr_parent.Children.insert(curr_idx, new_pair)
        new_pair.Parent = curr_parent
        del curr_parent.Children[curr_idx + 1]
        del curr_parent.Children[curr_idx - 1]
        return True

    def _pair_after_indices(self):
        """Returns indices of the pair after the current node, if possible.

        Returns tuple containing the indices if possible, None otherwise.
        """
        #can't work if one or fewer children
        if len(self) < 2:
            return None
        before = self[0]
        after = self[-1]
        if before.IsPaired or after.IsPaired:
            return None
        return (before.Start, after.Start)

    def pairAfter(self):
        """Forms a pair after the current node, if possible.

        Returns True if the pair was successfully created, False otherwise.
        Note that all children of self will become children of the newly-
        created pair.
        """
        indices = self._pair_after_indices()
        if not indices:
            return False
        #make the new pair and devolve children to it 
        new_pair = StructureNode(Data=Stem(indices[0], indices[1], 1))
        new_pair.Children[:] = self.Children[1:-1]    #all except start and end
        for c in new_pair.Children:
            c.Parent = new_pair
        self.Children[:] = [new_pair]
        new_pair.Parent = self
        return True

    def pairChildren(self, first, second):
        """Forms a pair using two of the children of self.

        first and second can be the node objects or the indices.

        Returns True if the pair was created, False otherwise.
        """
        #check that first and second are really instances of self, or
        #convert from indices
        if isinstance(first, int):
            if first >= 0:
                first_index = first
            else:
                first_index = len(self) + first
            first_node = self[first]
        else:
            first_index = first.Index
            first_node = first
        if isinstance(second, int):
            if second >= 0:
                second_index = second
            else:
                second_index = len(self) + second
            second_node = self[second]
        else:
            second_index = second.Index
            second_node = second
        #check that the nodes are really children of self, and that they're
        #not the same, and that they're not already paired
        if self is not first_node.Parent or self is not second_node.Parent:
            return False
        if first_node is second_node:
            return False
        if first_node.IsPaired or second_node.IsPaired:
            return False
        #create the new node and append the appropriate children
        new_node = StructureNode(Data=  \
            Stem(first_node.Start, second_node.Start, 1))
        del self.Children[second_index]
        new_node.Children[:] = self.Children[first_index + 1 : second_index]
        for c in new_node.Children:
            c.Parent = new_node
        self.Children[first_index] = new_node
        new_node.Parent = self
        return True

    def expand(self):
        """If self is a stem, expands self into n separate pair nodes.

        Returns True if the node was expanded, False otherwise.
        """
        if self.Length <= 1:
            return False
    
        children = self.Children[:]
        self.Children[:] = []
        curr_pair = self
        start = self.Start
        end = self.End
        for i in range(1, self.Length):
            new_pair = StructureNode(Data=Stem(start+i, end-i, 1))
            curr_pair.Children.append(new_pair)
            new_pair.Parent = curr_pair
            curr_pair = new_pair
        new_pair.Children[:] = children
        for c in new_pair.Children:
            c.Parent = new_pair
        self.Length = 1
        return True
            
    def expandAll(self):
        """Expands self and all children of self."""
        for i in self:
            i.expandAll()
        self.expand()

    def collapse(self):
        """If self is a stem, extends self with as many base pairs as possible.

        Returns True if the stem was extended, False otherwise.
        """
        #no effect if we're at the root, or if the position isn't paired
        if not self.IsPaired or self.Parent is None:
            return False
        extended = False
        while len(self) == 1 and self[0].IsPaired:
            self.Children[:] = self[0].Children[:]
            for c in self.Children:
                c.Parent = self
            self.Length += 1
            extended = True
        return extended

    def collapseAll(self):
        """Collapse self and all children of self."""
        self.collapse()
        for i in self:
            i.collapseAll()

    def breakBadPairs(self, seq):
        """Breaks all pairs in self that can't form in seq."""
        pairs = seq.Alphabet.Pairs
        self.expandAll()
        for i in self.traverse_recursive():
            if i.IsPaired:
                if not (seq[i.Start], seq[i.End]) in pairs:
                    i.unpair()

    def extendHelix(self, seq):
        """Extends helix represented by self as far as possible with seq.
        
        Assumes that helix has already been collapsed.
        """
        pairs = seq.Alphabet.Pairs
        #form as many contiguous pairs before as possible
        curr = self
        while 1:
            indices = curr._pair_before_indices()
            if indices and ((seq[indices[0]], seq[indices[1]]) in pairs):
                curr.pairBefore()
                #see if we can extend the newly-created pair
                curr = curr.Parent
            else:
                break
        #form as many contiguous pairs after as possible
        curr = self
        while 1:
            indices = curr._pair_after_indices()
            if indices and ((seq[indices[0]], seq[indices[1]]) in pairs) \
                and (curr.Stems or (len(curr.Unpaired) >= 5)):
                curr.pairAfter()
                #see if we can extend the newly-created pair
                curr = curr[0]
            else:
                break
    
    def extendHelices(self, seq):
        """Extends all helices in self and its children as far as possible."""
        for i in self.traverse_recursive():
            if i.IsPaired:
                i.extendHelix(seq)
    
    def fitSeq(self, seq):
        """Corrects the structure in self according to the specified sequence.

        seq must be a Sequence object with .Alphabet, etc.

        Note that fitSeq does not renumber the structure; it may be necessary
        to call renumber(), especially if the structure does not start at the
        beginning of the sequence.

        Assumes that it will be called starting with the root node; does not
        check backwards in the tree.
        """
        self.breakBadPairs(seq)
        self.extendHelices(seq)


def classify(struct, verbose=False):
    """Classifies a Vienna-format string into structural categories.
    
    struct may be a string or a real Vienna structure. It is tested on
    validity, because classifying invalid structures is very unreliable.
    If the number of closing brackets is larger than the number of opening
    brackets, an error would be raised, but if it's the other way around, 
    weird characters could be added to the string. For instance, classifying
    '.((.)' would give "\x00SSLS". This happens because the ends are not
    handled correctly if the stack is not back to its one-level state.
    """
    
    #Test whether structure is valid. Classifying invalid structures
    #is very unreliable.
    try:
        Vienna(struct)
    except IndexError:
        raise IndexError, "Trying to classify an invalid Vienna structure: %s"\
        %(struct)
    
    MAX_STEMS=1000

    #implement stack as three-item list
    PARENT = 0
    ITEMS = 1
    DEGREE = 2

    STEM, LOOP, BULGE, JUNCTION, END, FLEXIBLE = map(ord, 'SLBJEF')
    #WARNING: won't work if more than max_stems come off a junction
    LEVELS = [FLEXIBLE, LOOP, BULGE] + [JUNCTION]*MAX_STEMS

    length = len(struct)
    result = zeros((length,), 'B')
    stack = [None,[],0]
    curr_level = stack
    if verbose:
        print 'curr_level:', curr_level
        print 'result:', result

    for i, c in enumerate(struct):
        if verbose:
            print 'pos, char:',i,c
        #open parens add new level to stack
        if c == '(':
            curr_level = [curr_level,[],1]
            result[i] = STEM
        #unpaired base gets appended to current level
        elif c == '.':
            curr_level[ITEMS].append(i)
        #closed parens subtract level from stack and assign state
        elif c == ')':   #note: will handle end separately
            result[i] = STEM
            put(result, curr_level[ITEMS], LEVELS[curr_level[DEGREE]])
            curr_level = curr_level[PARENT]
            curr_level[DEGREE] += 1
        if verbose:
            print 'curr_level:', curr_level
            print 'result', result
            
    #handle ends and flexible bases
    end_items = curr_level[ITEMS]
    if end_items:
        first_start = struct.find('(')
        if first_start == -1:
            first_start = length+1
        last_end = struct.rfind(')')
        put(result, [i for i in end_items if first_start<i<last_end], FLEXIBLE)
        put(result, [i for i in end_items if not first_start<i<last_end], END)
    return result.tostring()

