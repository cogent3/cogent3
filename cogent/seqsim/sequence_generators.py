#!/usr/bin/env python
"""sequence_generators.py: various types of random and non-random generators.

Currently provides:

SequenceGenerator: fills in degenerate sequences, either by cycling through
    all possibilities or by jumping to a particular sequence. Supports indexing,
    iteration, and slicing.

Partition: generates all the ways of dividing n objects among b bins. Useful for
    stepping through a space of compositions, or dividing a sequence.

The SequencGenerators are fairly elaborate, and allow complex modeling of RNA.
The present implementation is based on Freqs, and is relatively slow. An
array-based implementation that uses seqsim.usage objects, cogent core
alphabets, etc. is in the works, and should have essentially the same interface.
However, this implementation is fairly well-tested and was used to generate
the data for the Knight et al. 2005 NAR paper on hammerhead/isoleucine
motif folding.
"""
from operator import mul
from types import SliceType
from sys import path
from random import choice, random, shuffle, randrange
from cogent.maths.stats.util import Freqs
from cogent.struct.rna2d import ViennaStructure
from cogent.app.vienna_package import RNAfold
from numpy import logical_and, fromstring, byte

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

IUPAC_DNA = {'T':'T','C':'C','A':'A','G':'G',
             'R':'AG','Y':'TC','W':'TA','S':'CG','M':'CA','K':'TG',
             'B':'TCG','D':'TAG','H':'TCA','V':'CAG','N':'TCAG'}
IUPAC_RNA = {'U':'U','C':'C','A':'A','G':'G',
             'R':'AG','Y':'UC','W':'UA','S':'CG','M':'CA','K':'UG',
             'B':'UCG','D':'UAG','H':'UCA','V':'CAG','N':'UCAG'}

def permutations(n, k):
    """Returns the number of ways of choosing k items from n, in order.

    Defined as n!/k!
    """
    #Validation: k must be be between 0 and n (inclusive), and n must be >=0.
    if k > n:
        raise IndexError, "can't choose %s items from %s" % (k, n)
    elif k < 0:
        raise IndexError, "can't choose negative number of items"
    elif n < 0:
        raise IndexError, "can't choose from negative number of items"
    product = 1
    for i in xrange(n-k+1, n+1):
        product *= i
    return product

def combinations(n, k):
    """Returns the number of ways of choosing k items from n.

    Defined as n!/(k!(n-k)!)
    """
    #Validation: k must be be between 0 and n (inclusive), and n must be >=0.
    if k > n:
        raise IndexError, "can't choose %s items from %s" % (k, n)
    elif k < 0:
        raise IndexError, "can't choose negative number of items"
    elif n < 0:
        raise IndexError, "can't choose from negative number of items"
    #permutations(n, k) = permutations(n, n-k), so reduce computation by 
    #figuring out which requires calculation of fewer terms.
    if k > (n - k):
        larger = k
        smaller = n - k
    else:
        larger = n - k
        smaller = k

    product = 1
    #compute n!/(n-larger)! by multiplying terms from n to (n-larger+1)
    for i in xrange(larger+1, n+1):
        product *= i

    #divide by (smaller)! by multiplying terms from 2 to smaller
    for i in xrange(2, smaller+1): #no need to divide by 1... 
        product /= i    #ok to use integer division: should always be factor

    return product


def _slice_support(the_slice, length):
    """Takes a slice and the length of an object; returns normalized version.

    Specifically, corrects start and end for negative indices.
    """
    start = the_slice.start
    stop = the_slice.stop
    step = the_slice.step or 1
    #fill in missing values for start and end
    if start is None:
        start = 0
    if stop is None:
        stop = length
    #convert end-relative values to start-relative values
    if start < 0:
        start = length - start
    if stop < 0:
        stop = length - stop
    return (start, stop, step)

    #ensure step is not zero, or we never move in the sequence!
    if not step:
        step = 1

class SequenceGenerator(object):
    """Generates all the possibilities for a degnerate template."""
    def __init__(self, template='', alphabet=None, start=None):
        """Returns a new SequenceGenerator based on template"""
        if alphabet:
            self.Alphabet = alphabet
        else:
            self.Alphabet = IUPAC_RNA
        self.Template = template
        if start:
            self.validate(start)
            self.Start = start
        else:
            self.Start = [0] * len(template)

    def validate(self, state):
        """Check that state is allowable given the template."""
        possibilities = map(len, map(self.Alphabet.__getitem__, self.Template))
        for max_allowed, curr in zip(possibilities, state):
            if curr >= max_allowed:
                raise ValueError, "Tried to set a state to too high an index."
        
    def __str__(self):
        """Returns data about current iterator's template"""
        return "<SequenceGenerator object with template: %s, and alphabet: %s>"\
            % (self.Template, self.Alphabet)

    def __len__(self):
        """Returns the number of elements in all possible expansions."""
        return self.numPossibilities()

    def numPossibilities(self):
        """Same as __len__, except Python doesn't coerce result to an int"""
        if self.Template:
            return reduce(mul, map(len, map(self.Alphabet.__getitem__, \
                    self.Template)))
        else:
            return 0

    def _index2state(self, index):
        """Takes an index and returns the corresponding state."""
        expansions = map(self.Alphabet.__getitem__, self.Template)
        num_items = len(expansions)
        lengths = map(len, expansions)
        indices = range(num_items)
        indices.reverse()     #want to traverse in reverse order
        states = [0] * num_items    #initialize with zero
        for i in indices:
            if not index:
                break   #if index is zero, so is everything to the left of i
            if lengths[i] == 1:
                continue #skip anything that can't vary
            else:
                choices = lengths[i]
                states[i] = index % choices
                index //= choices
        return states

    def __getitem__(self, index):
        """Supports indexing. Now uses constant-time algorithm (fast)."""
        if type(index) is SliceType:
            return self._handle_slice(index)
        else:
            if index < 0:
                index = self.__len__() + index
            iterator = self.items(self._index2state(index))
            return iterator.next()

    def _handle_slice(self, index):
        """Needs separate method, since __getitem__ can't yield _and_ return.
        """
        length = self.__len__() #might be too big to fit in an int
        start, stop, step = _slice_support(index, length)
        #quick check to se if we can't return any items
        if (stop - start < 1):
            raise StopIteration
        if step < 1:
            raise NotImplementedError, \
                "Can't support negative step in irreversible sequence."""
        else:
            index = start
            iterator = self.items(self._index2state(start))
            while index < stop:
                for i in range(step - 1):
                    if i >= stop - 1:
                        raise StopIteration
                    iterator.next()
                index += step
                yield iterator.next()

    def __iter__(self):
        """Iterator interface using self.Start as the start_state."""
        return self.items(self.Start)

    def items(self, start_state=None):
        """Acts like a sequence containing all the possibilities."""
        #shortcut if the template is empty
        if not self.Template:
            return
        #figure out how many possibilities there are at each position, and
        #what the choices are
        expansions = map(self.Alphabet.__getitem__, self.Template)
        num_positions = len(expansions)
        lengths = map(len, expansions)
        #set the starting state, i.e. the array of what the current choice is
        #at each position.
        if start_state is None:
            indices = [0] * num_positions
        else:
            self.validate(start_state)
            indices = start_state[:]
        seq = [expansions[i][indices[i]] for i in range(num_positions)]
        #always return the sequence for the first possibility: there might not
        #be any more...
        yield ''.join(seq)
        while 1: 
            #find rightmost element that can be incremented
            pos = num_positions - 1
            while indices[pos] == lengths[pos] - 1:
                pos -= 1
                if pos < 0: #ran off end
                    return
            indices[pos] += 1
            seq[pos] = expansions[pos][indices[pos]]
            #reset the rest of the elements, if there are any
            pos += 1
            while pos <= num_positions - 1:
                indices[pos] = 0
                seq[pos] = expansions[pos][0]
                pos += 1
            #seq should always contain a list of the current states for each pos
            yield ''.join(seq)

class Partition(object):
    """Generator behaving like a list of the partitions of a set of n items.

    Usage: p = Partition(num_items, num_pieces, min_occupancy=0)

    Requires each bin to have at least min_occupancy pieces.
    """

    def __init__(self, num_items, num_pieces, min_occupancy=0):
        """Returns new Partition object with first partitions initialized.
        
        Usage: p = Partition(num_items, num_pieces, min_occupancy=0)
        
        Default is min_occupanct items in each bin, with all the leftovers
        in the first bin.
        """
        self.NumItems = num_items
        if num_pieces:
            self.NumPieces = num_pieces
        else:
            raise ValueError, "Cannot divide items among zero bins."
        self.MinOccupancy = min_occupancy
        self._reset()

    def __str__(self):
        """Prints string representation with pieces, items, and occupancy."""
        return "Items: %s Pieces: %s Min Per Piece: %s" % \
        (self.NumItems, self.NumPieces, self.MinOccupancy)

    def _validate(self, states):
        """Verify that states has right sum and meets occupancy restrictions.
        
        Raises ValueError if there is any problem: does not return anything.
        """
        num_pieces = self.NumPieces
        min_occupancy = self.MinOccupancy #cache for efficiency
        #check the number of pieces
        if len(states) != num_pieces:
            raise ValueError, "Tried to set state %s, but need %s pieces." % \
            (states, num_pieces)
        #check that no piece has too few items
        sum = 0
        for state in states:
            if state < min_occupancy:
                raise ValueError, \
                "Tried to set state %s, but need at least %s items per bin." %\
                (states, min_occupancy)
            sum += state
        #check that we have the right number of pieces
        if sum != self.NumItems:
            raise ValueError, \
            "Tried to set state %s, but it has %s pieces instead of %s." %\
            (states, sum, self.NumItems)
        

    def _reset(self, states=None):
        """Resets to a particular state given by sequence of states in each bin.
        
        Default: go to first partition, with MinOccupancy items in each bin and
        any leftovers in the first bin.
        """
        min_occupancy = self.MinOccupancy
        num_items = self.NumItems
        num_pieces = self.NumPieces #cache for efficiency
        if states:
            #check that we're not trying to set a bad state
            self._validate(states)
            self._bins = states
        else:
            reserved = (num_pieces - 1) * min_occupancy
            #check that we can actually divide the pieces among the bins OK
            if reserved + min_occupancy > num_items:
                raise ValueError, \
                "Can't divide %s items into %s pieces with at least %s in each."\
                % (num_items, num_pieces, min_occupancy)
            #otherwise, fill the bins    
            bins = [min_occupancy] * num_pieces
            bins[0] = num_items - reserved
            self._bins = bins
            self._reserved = reserved

    def __iter__(self):
        """Defines iterator interface, starting with self._bins."""
        return self.items()

    def _transform(self, value):
        """Transformation to be applied to return values.
        
        Default behavior is to copy, but can be overridden in derived classes.
        """
        return value[:]

    def items(self, bin_states=None):
        """Defines iterator interface, supporting for i in self."""
        #always copy the array of states, since we will be mutating it
        if bin_states:
            self._validate(bin_states)
            bins = bin_states[:]
        else:
            self._reset()
            bins = self._bins[:]
        #cache local vars for efficiency
        delta = self.MinOccupancy
        num_items = self.NumItems
        num_pieces = self.NumPieces
        transform = self._transform
        end_state = num_items - (delta * (num_pieces - 1))
        #always return the first state
        yield transform(bins)
        while 1:
            #check if we're done: when the last bin has all the items
            if bins[-1] == end_state:
                return
            #need to adjust the bins to the correct state for next time
            #find rightmost non-delta except the last
            rightmost = sum = 0
            #figure out the sum of all the items to the right of the bin we're
            #going to decrement, and also which bin we're going to decrement
            for i in xrange(len(bins)-2, -1, -1):
                curr = bins[i]
                if curr != delta:
                    rightmost = i
                    break
                else:
                    sum += curr
            #bins[-1] excluded from count above: also need to add 1 for newly
            #incremented item from the rightmost decrementable bin
            sum += bins[-1] + 1
            bins[rightmost] -= 1
            #leftover_bins counts the number of bins more than one to the left
            #of the rightmost
            leftover_bins = num_pieces - rightmost - 2
            if leftover_bins:
                bins[rightmost+2:] = [delta] * leftover_bins
                sum -= delta * leftover_bins
            bins[rightmost + 1] = sum
            yield transform(bins)

    def __len__(self):
        """Calculates the number of possible parameters with current state.

        Specifically, only takes into account the number of objects, the
        number of bins, and the minimum per bin: does _not_ take into account
        a particular start point.
        """
        #NOTE: I don't know why the following works, but it seems to be
        #empirically true when compared to the lengths of the resulting lists.
        cuts = self.NumPieces 
        items = self.NumItems - (self.MinOccupancy - 1) * cuts
        product = 1
        for i in range(items - cuts + 1, items):
            product *= i
        for i in range(2, cuts):
            product /= i
        return product

class Composition(Partition):
    """Generates evenly spaced composition intervals over an alphabet.

    Usage:c=Composition(spacing,min_occupancy=0,alphabet='ACGU')

    spacing should be a float representing the percentage of the space 
    separating successive values (e.g. 5 for 5% steps). Note that the spacing 
    may be approximated: check self.Spacing to see what the recorded value is.

    alphabet should be a list, in order, of the possible characters.

    min_occupancy should typically be 0 (can miss some symbols) or 1 (always
    require at least one of each symbol).

    Always yields an un-normalized Freqs containing counts of
    each symbol at each step.

    For a given alphabet A, the possible compositions of that alphabet can
    be represented as a simplex in len(A)-1 dimensions. Composition returns
    a representation of evenly distributed compositions in that space, with
    distances along all dimensions represented by spacing (i.e. if spacing
    is 0.05, then the next point in any dimension will be 0.05 away if it
    exists.)

    Useful for generating sequences of specified composition that can then
    be randomized.
    """
    def __init__(self,spacing,min_occupancy=0,alphabet='ACGU'):
        """Initializes new generator with specified spacing, alphabet, etc.

        Usage: c = Composition(spacing, min_occupancy=0,
                   alphabet='ACGU')

        See class documentation for details.
        """
        self.Spacing = spacing #also sets self._num_items
        self.Alphabet = alphabet #also sets self._num_pieces
        self.MinOccupancy = min_occupancy

    def _get_spacing(self):
        """Accessor for self.Spacing."""
        return self._spacing

    def _set_spacing(self, spacing):
        """Mutator for self.Spacing. Sets self.NumItems to correct value."""
        num_items = int(round(100.0/spacing))
        self._num_items = num_items
        self._spacing = 100.0/num_items

    Spacing = property(_get_spacing, _set_spacing, \
        doc="Set spacing and calculate string length.")

    def _get_num_items(self):
        """Accessor for self.NumItems."""
        return self._num_items
    def _set_num_items(self, num_items):
        """Mutator for self.NumItems: recalculates self.Spacing."""
        self._num_items = num_items
        self._spacing = 100.0/num_items

    NumItems = property(_get_num_items, _set_num_items, \
        doc="Set NumItems and calculate Spacing.")
    
    def _get_alphabet(self):
        """Accessor for self.Alphabet."""
        return self._alphabet
    def _set_alphabet(self, alphabet):
        """Mutator for self.Alphabet."""
        self._alphabet = alphabet

    Alphabet = property(_get_alphabet, _set_alphabet, \
        doc="Set alphabet and calculate number of pieces.")

    def _get_num_pieces(self):
        """Accessor for self.NumPieces"""
        return len(self.Alphabet)

    NumPieces = property(_get_num_pieces, doc="Get number of pieces.")

    def _transform(self, value):
        """Override superclass transform to yield Freqs.
        """
        return Freqs(dict(zip(self.Alphabet, value)))
       
    def __iter__(self):
        """Defines iterator interface, starting with self._bins."""
        return self.items()

class MageFrequencies(object):
    """Takes a Freqs and optionally a label.

    Writes out a Mage-format string.

    This presentation class is standalone to avoid cluttering 
    Freqs.
    """
    def __init__(self, freqs, label=''):
        """Returns a new MageFrequencies object.

        This is basically a labeled Freqs that can write itself
        out as a Mage-format string.
        """
        self.Freqs = freqs#don't mutate original
        self.Label = label
    
    def __str__(self):
        """Returns the frequency string, suitable for MAGE."""
        pieces = []
        freqs = self.Freqs
        known_bases = Freqs({
            'A':freqs.get('A',0),
            'C':freqs.get('C',0),
            'U':freqs.get('U',0) + freqs.get('T',0),
            'G':freqs.get('G',0),
            })
        #frequencies should sum to 1 for MAGE display.
        known_bases.normalize()
        #Only append label field if there is one.
        label = self.Label or ''
        if label:
            pieces.append('{%s}' % self.Label)
        for item in 'ACG':
            pieces.append(str(known_bases[item]))
        return ' '.join(pieces)

class SequenceHandle(list):
    """Holds mutable sequence that can join itself together as string.
    
    Sequence cannot vary in length.
    """
    def __init__(self, data='', alphabet=None):
        """Initializes new list over an alphabet. Rejects invalid entries."""
        if alphabet:
            for d in data:
                if d not in alphabet:
                    raise ValueError, "Item %s not in alphabet %s." \
                    % (d, alphabet)
        super(SequenceHandle, self).__init__(data)
        self.Alphabet = alphabet

    def __setitem__(self, index, item):
        """Checks that the item is in the alphabet."""
        alphabet = self.Alphabet
        if alphabet:
            try:
                absent = item not in alphabet
            except TypeError:
                raise ValueError, "Item %s not in alphabet %s." \
                % (item, alphabet)
            else:
                if absent:
                    raise ValueError, "Item %s not in alphabet %s." \
                    % (item, alphabet)
        super(SequenceHandle, self).__setitem__(index, item)

    def __setslice__(self, start, stop, values):
        """Checks that items are in alphabet, and that slice is same length."""
        orig_length = len(self)
        alphabet = self.Alphabet
        if alphabet:
            for v in values:
                try:
                    absent = v not in alphabet
                except TypeError:
                    raise ValueError, "Item %s not in alphabet %s." \
                    % (v, alphabet)
                else:
                    if absent:
                        raise ValueError, "Item %s not in alphabet %s." \
                        % (v, alphabet)
        super(SequenceHandle, self).__setslice__(start, stop, values)
        if len(self) != orig_length:
            raise ValueError, "Cannot change length of SequenceHandle."

    def __str__(self):
        """Returns self as a string, symbols joined."""
        try:
            return ''.join(self)
        except: #use built-in conversion methods for lists
            return super(SequenceHandle, self).__str__()

    def _naughty_method(self, *args, **kwargs):
        """Prevent other methods that change the length or set items."""
        raise NotImplementedError, \
            "May not change length of SequenceHandle."
    #note how _many_ methods are naughty... 
    __delitem__ = __delslice__ = __iadd__ = __imul__ = append \
    = extend = insert = pop = remove = _naughty_method
 
class BaseFrequency(Freqs):
    RNA = ['U', 'C', 'A', 'G']
    DNA = ['T', 'C', 'A', 'G']

    """Holds information about base frequencies."""

    def __init__(self, freqs, RNA=True):
        """Returns new BaseFrequency object, ensuring a count for each base."""
        if RNA:
            alphabet = self.RNA
        else:
            alphabet = self.DNA
        super(BaseFrequency, self).__init__(freqs, alphabet)
        for k in alphabet:
            if k not in self:
                self[k] = 0.0

class PairFrequency(Freqs):
    """Makes a frequency distribution of pairs from freqs of single items."""

    def __init__(self, freqs, pairs=None):
        """Makes pair frequency distribution.

        Usage: p = PairFrequency(freqs, pairs)

        freqs is the single-item frequencies
        pairs is the list of valid pairs from which samples will be drawn.
        If pairs is None (the default), constructs all possible pairs.
        """
        symbol_freqs = BaseFrequency(freqs)
        if pairs is None:
            symbols = symbol_freqs.keys()
            pairs = [(i, j) for i in symbols for j in symbols]
        pair_freqs = {}
        for i, j in pairs:
            try:
                pair_freqs[(i,j)] = symbol_freqs[i]*symbol_freqs[j]
            except KeyError, e:
                print symbol_freqs
                print i, j
                raise e
        super(PairFrequency, self).__init__(pair_freqs, pairs)
        self.normalize()

class BasePairFrequency(PairFrequency):
    """Holds information about base pair frequencies."""
    WatsonCrick = [('A','U'), ('U','A'),('G','C'),('C','G')]
    Wobble = WatsonCrick + [('G','U'), ('U','G')]

    def __init__(self, freqs, GU=True):
        if GU:
            pairs = self.Wobble
        else:
            pairs = self.WatsonCrick
        super(BasePairFrequency, self).__init__(freqs, pairs)   

class RegionModel(object):
    """Holds probability model for constructing random or randomized sequences.

    Supports the following interface:

        Current:     Reference to current sequence, or tuple of references.
        Template:    Degenerate sequence specifying the class of sequences to
                     produce. Immutable.
        Length:      Length of the current sequence. Read-only.
        Composition: Composition used to generate sequences (e.g. pairs).
        refresh():   Generate the next, random sequence.
        monomers(f): Update internal frequencies using symbol frequencies in f.
        
    Base class RegionModel behavior is to model a constant region.
    """
    def __init__(self, template='', composition=None):
        """Return a new RegionModel object. See class for documentation."""
        self.Composition = composition
        self.Template = template    #will set self.Current

    def _get_template(self):
        """Accessor method for self.Template"""
        return self._template
    
    def _set_template(self, data):
        """Mutator method for self.Template"""
        self._template = data
        self.Current = SequenceHandle(data)
        self.refresh()

    Template = property(_get_template, _set_template)

    def _get_composition(self):
        """Accessor method for self.Composition"""
        return self._composition

    def _set_composition(self, composition):
        """Mutator method for self.Composition"""
        self._composition = composition
        self.refresh()
        
    Composition = property(_get_composition, _set_composition)

    def __len__(self):
        """Returns length of the current string."""
        return len(self.Current)

    def refresh(self):
        """Replaces the current sequence with a new string fitting the model.

        Does nothing unless overridden in derived classes.
        """
        pass

    def monomers(self, composition, **kwargs):
        """Replaces the current composition with new Freqs."""
        self.Composition = composition  #no effect unless overridden

class ConstantRegion(RegionModel):
    """Holds a constant string: this is default behavior."""
    pass

class UnpairedRegion(RegionModel):
    """Holds an unpaired region: this gets filled in from self.Composition."""
    def refresh(self):
        """Fills in a sequence drawn randomly from composition."""
        if hasattr(self, "Current") and self.Current and self.Composition:
            self.Current[:] = self.Composition.randomSequence(len(self))

class ShuffledRegion(RegionModel):
    """Holds a non-degenerate template that is randomized by shuffling."""
    def refresh(self):
        """Randomizes the template by permuting the elements."""
        if hasattr(self, "Current") and self.Current:
            shuffle(self.Current)

class PairedRegion(RegionModel):
    """Holds complementary upstream and downstream strands."""
    def refresh(self):
        """Fills in tuple of paired sequences drawn from self.Composition."""
        if hasattr(self, "Current") and self.Current and self.Composition:
            length = len(self)
            upstream = self.Current[0]
            downstream = self.Current[1]
            pairs = self.Composition.randomSequence(length)
            for i in xrange(length):
                upstream[i] = pairs[i][0]
                downstream[i] = pairs[i][1]
            #downstream has the complements, but need to reverse it as well
            downstream.reverse()

    def __len__(self):
        """Returns length of (half of) the current helix, not the tuple..."""
        return len(self.Current[0])

    def _set_template(self, data):
        """Mutator method for self.Template"""
        data = list(data)
        self._template = data
        self.Current = (SequenceHandle(data), SequenceHandle(data))
        self.refresh()

    #Override base class _set_template in the property
    Template = property(RegionModel._get_template, _set_template)
    
    def monomers(self, composition, **kwargs):
        """Calculates pair distribution from monomer frequencies."""
        GU = kwargs.get('GU', True)
        self.Composition = BasePairFrequency(composition,GU)

class DegenRegion(RegionModel):
    """Handles a string of degenexerate bases.
    
    WARNING: Not tested!
    """
    def refresh(self):
        """Fills in degen bases randomly according to possible symbols"""
        if hasattr(self, "Current") and self.Current and self.Composition:
            result = []
            non_degen = dict.fromkeys('UCAG')
            freqs = self.Composition
            for b in self.Template:
                if b in non_degen:
                    result.append(b)
                else:
                    allowed_bases = IUPAC_RNA[b]
                    composition = Freqs(dict([(i,freqs[i]) for i in allowed_bases]))
                    composition.normalize()
                    result.append(composition.choice(random()))
            self.Current[:] = ''.join(result)

###WARNING: MATCHINGREGION HAS NOT YET BEEN TESTED####
class MatchingRegion(RegionModel):
    """Fills in the complement to specified constant region."""
    WatsonCrick = {'A':'U', 'U':'A', 'C':'G', 'G':'C'}
    Wobble = {'A':'U', 'U':'AG', 'C':'G', 'G':'UC'}

    def __init__(self):
        raise NotImplementedError, "NOT YET TESTED"

    def _init_current(self):
        """Initializes Current and some private variables."""
        wc = self.WatsonCrick
        template = self.Template
        self.Complement = [wc[base] for base in template]
        self.Current = SequenceHandle(self.Complement)
        self._freqs = {}
        
    def refresh(self):
        """Returns new sequence that could pair with template."""
        if self.GU:
            freqs = self._freqs
            bases = [freqs[base].randomSequence(1)[0] for base in self.Template]
            self.Current[:] = bases
        else:
            self.Current[:] = self.Complement

    def monomers(self, freqs, **kwargs):
        """Calculates Freqs of possibilities for each base."""
        if kwargs.get('GU', False):
            pairs = self.Wobble.items()
            self.GU = True
            for base, complements in pairs:
                freqs = self.Composition.copy()
                freqs.subset(complements)
                freqs.normalize()
                self._freqs[base] = freqs
        else:
            self.GU = False

class SequenceModel(object):
    """Stores state associated with generating a randomized sequence."""
    def __init__(self, order, composition=None, GU=True,    \
        constants=[], unpaireds=[], helices=[], matches=[], degenerates=[]):
        """Returns a new SequenceModel.

        constants, unpaireds, and helices should all be lists of RegionModels.
        order should be a string of the following format:
            [label1][index1] [label2][index2] ...
            where label is C for constant, U for unpaired, or H for helix,
            or D for degen,
            and index is the index of the region within the appropriate list.

            Use hyphens to indicate cuts.

            For example:
                "C0 U3 H1 U5 - H2 C1 H2 H1 U0"
                ...means constants[0] followed by unpaireds[3], followed by the
                first part of helices[1], followed by unpaireds[5], followed
                by the first part of helices[2], followed by constants[1],
                followed by the second part of helices[2], followed by the
                second part of helices[1], followed by unpaireds[0].
                There is a cut in the sequence between U5 and H2.

        Although this is a very general mechanism (each piece can potentially
        have its own composition, etc.), typically the functionality will be
        accessed programatically through other classes.
        """
        self.Helices = helices
        self.Unpaireds = unpaireds
        self.Constants = constants
        self.Degenerates = degenerates
        self.Matches = matches
        self.GU = GU
        #Don't _require_ a composition to be passed in, but if it isn't passed
        #in, then all the pieces must be initialized with their own compositions
        #beforehand.
        self.Composition = composition  
        self.Order = order

    def __len__(self):
        """Figures out the total length of all the components."""
        length = 0
        for i in self.Unpaireds + self.Constants + self.Matches + self.Degenerates:
            length += len(i)
        for h in self.Helices:
            length += 2 * len(h)
        return length
    
    def refresh(self):
        """Delegates each region to refresh itself."""
        for i in self.Helices + self.Unpaireds + self.Matches + self.Degenerates:
            i.refresh()

    def _get_order(self):
        """Accessor for self.Order."""
        return self._order

    def _set_order(self, order):
        """Figure out the order to put pieces in, using string format."""
        result = []
        segments = order.split('-')
        helix_counts = [0] * len(self.Helices)
        for s in segments:
            pieces = []
            components = s.split()
            for c in components:
                label = c[0]
                index = int(c[1:])
                if label == 'C':    #constant
                    pieces.append(self.Constants[index].Current)
                elif label == 'U':  #unpaired random region
                    pieces.append(self.Unpaireds[index].Current)
                elif label == 'D':  #degenerate region
                    pieces.append(self.Degenerates[index].Current)
                elif label == 'H':  #helix
                    pieces.append(self.Helices[index].Current[\
                    helix_counts[index]])
                    helix_counts[index] += 1
                    #will give IndexError if the helix is added too many times
                else:
                    raise ValueError, \
                    "SequenceModel got unknown label: %s" % label
            result.append(pieces)
        self._order = result
        self.refresh()

    Order = property(_get_order, _set_order, \
        doc="Stores order for accessing the pieces of the template.")

    def _get_composition(self):
        """Accessor for Composition."""
        return self._composition

    def _set_composition(self, composition):
        """Sets the composition of each of the components to a global value."""
        if composition:
            for i in self.Helices + self.Unpaireds + self.Matches + self.Degenerates:
                i.monomers(composition, GU=self.GU)
        self._composition = composition

    Composition = property(_get_composition, _set_composition)

    def _get_GU(self):
        """Accessor for GU."""
        return self._GU

    def _set_GU(self, GU):
        """Mutator for GU. Recalculates composition."""
        self._GU = GU
        if hasattr(self, 'Composition') and self.Composition:
            self.Composition = self.Composition #recalculate with GU

    GU = property(_get_GU,_set_GU,doc="Controls whether GU pairs are allowed.")

    def __getitem__(self, index):
        """Returns the index'th segment of the sequence in its current state."""
        return ''.join([str(i) for i in self.Order[index]])

    def __str__(self):
        return '-'.join(self)

class Rule(object):
    """Holds information about pairing constraints on motifs."""
    def __init__(self, upstream_seq, upstream_pos, downstream_seq, \
        downstream_pos, length):
        """Initialize new Rule object."""
        self.UpstreamSequence = upstream_seq
        self.UpstreamPosition = upstream_pos
        self.DownstreamSequence = downstream_seq
        self.DownstreamPosition = downstream_pos
        self.Length = length
        self.validate()

    def validate(self):
        """Sanity checks on rule object."""
        if self.Length <= 0:
            raise ValueError, "Helix length must be at least 1."
        if self.Length > self.DownstreamPosition + 1:
            raise ValueError, \
            "Helix length cannot be more than 1 greater than downstream start."
        if min(self.UpstreamSequence, self.UpstreamPosition, \
            self.DownstreamSequence, self.DownstreamPosition) < 0:
            raise ValueError, \
            "All sequences and positions must be >= 0."
        if self.UpstreamSequence == self.DownstreamSequence:
            if self.UpstreamPosition >= self.DownstreamPosition:
                raise ValueError, \
                "Upstream position must have lower index than downstream."
            if self.DownstreamPosition-self.UpstreamPosition+1 < 2*self.Length:
                raise ValueError, "Helices can't overlap."
        if self.UpstreamSequence > self.DownstreamSequence:
            raise ValueError, "Upstream sequence must have the smaller index."

    def isCompatible(self, other):
        """Checks that the helices in self and other don't overlap.

        Has to try all possible combinations of upstream and downstream
        sequences, since any could conflict.
        """
        if self.UpstreamSequence == other.UpstreamSequence:
            diff = abs(self.UpstreamPosition - other.UpstreamPosition)
            if self.UpstreamPosition <= other.UpstreamPosition:
                if diff < self.Length:
                    return False
            elif diff < other.Length:
                return False
                
        if self.DownstreamSequence == other.DownstreamSequence:
            diff = abs(self.DownstreamPosition - other.DownstreamPosition)
            if self.DownstreamPosition >= other.DownstreamPosition:
                if diff < self.Length:
                    return False
            elif diff < other.Length:
                return False

        if self.UpstreamSequence == other.DownstreamSequence:
            diff = abs(self.UpstreamPosition - other.DownstreamPosition)
            #only need to check if position in self <= position in other
            if self.UpstreamPosition <= other.DownstreamPosition:
                if diff < (self.Length + other.Length - 1):
                    return False

        if self.DownstreamSequence == other.UpstreamSequence:
            diff = abs(self.DownstreamPosition - other.UpstreamPosition)
            #only need to check if position in self >= position in other
            if self.DownstreamPosition >= other.UpstreamPosition:
                if diff < (self.Length + other.Length - 1):
                    return False
        #if none of the checks failed, return True
        return True

    def fitsInSequence(self, upstream):
        """Checks whether upstream sequence is too short to hold helix.
        
        Note: downstream sequence length doesn't need to be checked because
        the index that can't be overlapped is always 0.
        """
        if self.UpstreamPosition + self.Length > len(upstream):
            return False
        else:
            return True
        
    def __str__(self):
        """Human-readable rule string."""
        return "Up Seq: %s Up Pos: %s Down Seq: %s Down Pos: %s Length: %s" % \
            (self.UpstreamSequence, self.UpstreamPosition, \
             self.DownstreamSequence, self.DownstreamPosition, self.Length)

class Module(object):
    """Holds information about a module's required sequence and structure."""
    def __init__(self, sequence, structure):
        """Returns a new Module object with specified sequence and structure."""
        self.Sequence = sequence
        self.Structure = structure
        len(self)   #will raise error if lengths out of sync

    def __len__(self):
        """Returns length of sequence and structure."""
        seq = self.Sequence
        struct = self.Structure
        seq_length = len(seq)
        if seq_length != len(struct):
            raise ValueError, \
            "Lengths of sequence '%s' and structure '%s' differ." % \
            (seq, struct)
        else:
            return seq_length

    def __str__(self):
        """Returns string containing sequence and structure."""
        return "Sequence:  %s\nStructure: %s" % (self.Sequence, self.Structure)

    def matches(self, other, index=None, alphabet=IUPAC_RNA):
        """Tests whether sequence/structure in self match other at index.

        other must be an object that has Sequence and Structure properties.
        If index is None, will search for matches anywhere in other.
        
        ###THIS METHOD NEEDS ATTENTION: move responsibility for finding matches
        to the sequence objects themselves?
        """
        length = len(self)
        if not length:  #zero-length pattern matches everywhere by definition
            return True 
        if index is not None:   #index might be 0...
            seq_match = True
            this_seq = self.Sequence
            other_seq = other.Sequence
            for i in range(length):
                curr = False
                try:
                    curr = curr or (other_seq[i+index] in alphabet[this_seq[i]])
                except:
                    pass
                if not curr:
                    try:
                        curr = curr or (this_seq[i] in \
                            alphabet[other_seq[i+index]])
                    except:
                        pass
                if not curr:
                    seq_match = False
                    break
            struct_match = self.structureMatches(other.Structure, index)
            return seq_match and struct_match[0]
        else:
            other_length = len(other)
            seq = self.Sequence
            struct = self.Structure
            other_struct = other.Structure
            other_seq = other.Sequence
            curr = 0    #current index
            while curr <= other_length - length:   #don't run off end
                try:
                    index = other_seq.index(seq, curr)
                except ValueError:
                    return False    #no more matches to try
                if struct == other_struct[index:index+length]:
                    return True     #found struct and seq matches at same place
                if curr == index:
                    curr += 1       #always make sure curr is incremented
                else:
                    curr = index
            return False    #must have been a seq match but no struct match
                            #at last window if we got here after the loop

    def structureMatches(self, structure, index):
        """Tests whether structure in self matches other at index.
        
        structure must have PairList property, e.g. ViennaStructure.
        """
        length = len(self)
        if not length:  #zero-length pattern matches everywhere by definition
            return (True, )
        else:
            ss = fromstring(self.Structure, byte)
            structure_mask = ss != ord('x')
            diffs = ss != fromstring(structure[index:index+length], byte)
            result = not logical_and(diffs, structure_mask).any()
            return result, ss, structure_mask, diffs
            
class Motif(object):
    """Holds sequences and structures for a motif."""
    def __init__(self, modules, rules):
        """Initializes motif with sequences, structures, and rules"""
        self.Modules = modules
        self.Rules = rules                      
        self.validate()
        
    def validate(self):
        """Checks that sequences and structures are equal length, and rules ok.

        Specifically, there must be the same number of sequences as structures;
        the length of each sequence must be the length of each structure; and
        the rules may not refer to any index outside the known sequences and
        structures.
        """
        self._check_helix_lengths()
        self._check_rule_overlaps()
    
    def _check_helix_lengths(self):
        """Check upstream sequence of each rule to make sure the helix can fit."""
        for r in self.Rules:
            if not r.fitsInSequence(self.Modules[r.UpstreamSequence].Sequence):
               raise ValueError, "Rule '%s' can't fit in sequence '%s'." \
               % (r, self.Modules[r.UpstreamSequence].Sequence)
                
    def _check_rule_overlaps(self):
        """Check every pair of rules for overlaps in covered regions."""
        rules = self.Rules  #cache reference for efficiency
        for first in range(len(rules)):
            first_rule = rules[first]
            for second in range(first):
                second_rule = rules[second]
                if not first_rule.isCompatible(second_rule):
                    raise ValueError, "Rules '%s' and '%s' incompatible." \
                    % (first, second)

    def _check_rule_match(self, rule, pairlist, locations):
        """Check whether rule matches pairlist given module locations.
        
        pairlist should be a list where, for each position in a longer sequence,
        parlist[i] should be the index of the partner of i, or None if i is
        not paired.

        locations should be a list of the locations of each module, in the order
        that the rule expects to find them.
        """
        start_up = locations[rule.UpstreamSequence]+rule.UpstreamPosition
        start_down = locations[rule.DownstreamSequence]+rule.DownstreamPosition
        for i in range(rule.Length):
            if pairlist[start_up + i] != start_down - i:
                return False
        return True     #if nothing failed, everything must be OK

    def _get_rule_match_pairs(self, rule, pairlist, locations):
        """Get the pairs that the rule will check."""
        start_up = locations[rule.UpstreamSequence]+rule.UpstreamPosition
        start_down = locations[rule.DownstreamSequence]+rule.DownstreamPosition
        return [(start_up+i,start_down-i) for i in range(rule.Length)]

            
    def matches(self, sequence, structure, positions):
        """Checks that sequence and structure matches motifs/rules.

        sequence needs to support the string interface (specifically, s.index)
        if it is necessary to search for matches anywhere; otherwise, arbitrary
        sequences should work.

        structure must have a PairList property (like ViennaStructure),
        which is a list the same length as the sequence where the value of
        each position is the index of its partner, or None if it is unpaired.
        As with sequence, must support string interface to find arbitrary
        matches; arbitrary sequences are ok otherwise.

        positions must be a list the same length as self.Modules, containing the
        index at which each successive module should be searched for.

        ###TO BE IMPLEMENTED: IF POSITIONS IS NONE, SEARCH FOR THE MODULES
        ANYWHERE IN THE SEQUENCE.###
        """
        full_length = Module(sequence, structure)   #more convenient as object
        modules = self.Modules
        if len(positions) != len(modules):
            raise ValueError, "len(positions) must match number of modules."
        for position, module in zip(positions, modules):
            if not module.matches(full_length, position):
                return False
        #can only get here if all the modules matched: need to check rules
        pairlist = structure.toPartners()
        for rule in self.Rules:
            if not self._check_rule_match(rule, pairlist, positions):
                return False
        #if we got here, all the modules matched and all the rules were OK
        return True
 
    def structureMatches(self, structure, positions, offsets=None,debug=False):
        """Checks that structure only matches motifs/rules.

        structure must have a PairList property (like ViennaStructure),
        which is a list the same length as the sequence where the value of
        each position is the index of its partner, or None if it is unpaired.
        As with sequence, must support string interface to find arbitrary
        matches; arbitrary sequences are ok otherwise.

        positions must be a list the same length as self.Modules, containing the
        index at which each successive module should be searched for.

        ###TO BE IMPLEMENTED: IF POSITIONS IS NONE, SEARCH FOR THE MODULES
        ANYWHERE IN THE SEQUENCE.###
        """
        modules = self.Modules
        if len(positions) != len(modules):
            raise ValueError, "len(positions) must match number of modules."
        if offsets:
            positions = [p+o for p, o in zip(positions, offsets)]
        result = True
        for position, module in zip(positions, modules):
            matched, ss, mask, diffs = \
                module.structureMatches(structure, position) 
            if debug:
                print 'STRUC:', structure[position:position+len(ss)]
                print 'SS   :', ss.tostring()
                print 'MASK :', ''.join(map(str, map(int, mask)))
                print 'DIFFS:', ''.join(map(str, map(int,diffs)))
                print 'WHERE:'
                all = ['.'] * len(structure)
                all[position:position+len(ss)] = ['x']*len(ss)
                print ''.join(all)
            if not matched: 
                if debug:
                    result = False
                else:
                    return False
        if not result:
            return False
        #can only get here if all the modules matched: need to check rules
        pairlist = structure.toPartners()
        for rule in self.Rules:
            if debug:
                pairs = self._get_rule_match_pairs(rule, pairlist, positions)
                for up, down in pairs:
                    all = ['.'] * len(structure)
                    all[up] = '('
                    all[down] = ')'
                    print ''.join(all)
                    if not pairlist[up] == down:
                        print structure
                        raise Exception, "Failed to find partner in pairlist"
            if not self._check_rule_match(rule, pairlist, positions):
                return False
        #if we got here, all the modules matched and all the rules were OK
        return True
        

class SequenceEmbedder(object):
    """Generates and analyzes set of modules embedded inside longer sequence."""
    def __init__(self, length, num_to_do, motif, model, composition, GU=True,\
                 with_replacement=False, positions=None, primer_5='',
                 primer_3='', match_offsets=None, debug=False, report_seqs=False
                 ):
        """Initializes with a specified length sequence model, composition.
        
        Note that sampling with replacement does NOT give all the outcomes
        equal frequencies, e.g. with two choices (0,1) will happen half the
        time because there are 2 ways to get it, but only one way to get
        (0,0) or (1,1).
        """
        self.Model = model
        self.Motif = motif
        self.NumToDo = long(num_to_do)
        self.Length = long(length)
        self.WithReplacement = with_replacement #allows adjacent modules
        self.GU = GU
        self.RandomRegion = UnpairedRegion('N'*(length - len(self.Model)), \
            composition)
        self.Composition = composition
        self._fixed_positions = positions
        self.Positions = positions
        self.Primer3 = primer_3
        self.Primer5 = primer_5
        self.MatchOffsets = match_offsets
        self.Debug = debug
        self.ReportSeqs = report_seqs

    def _get_composition(self):
        """Accessor for self.Composition."""
        return self._composition

    def _set_composition(self, composition):
        """Mutator for self.Composition."""
        self._composition = composition
        if composition:
            self.Model.GU = self.GU
            self.Model.Composition = composition
            self.RandomRegion.Composition = composition

    Composition = property(_get_composition, _set_composition)

    def _choose_locations(self):
        """Picks out places for the modules."""
        random_positions = self.Length - len(self.Model)
        num_modules = len(self.Motif.Modules)
        locations = []
        with_replacement = self.WithReplacement
        if (not with_replacement) and (random_positions < num_modules):
            raise ValueError, "Not enough positions to place modules."
        while len(locations) < num_modules:
            if with_replacement:
                curr = randrange(random_positions + 1)
                locations.append(curr)
            else:
                curr = randrange(random_positions)
                if curr not in locations:
                    locations.append(curr)
        locations.sort()
        return locations

    def __str__(self):
        """Makes a new sequence with inserts at correct positions.
        
        Note: no longer mutates self.Positions.
        """
        pieces = [str(self.Primer5)]
        random = str(self.RandomRegion.Current)
        modules = list(self.Model)
        added_positions = 0
        last_position = 0
        positions = self.Positions[:]
        for i in range(len(positions)):
            curr_module = modules[i]
            curr_position = positions[i]
            pieces.append(random[last_position:curr_position])
            pieces.append(curr_module)
            last_position = curr_position
            positions[i] += added_positions
            added_positions += len(curr_module)
        pieces.append(random[last_position:])   #add anything left over
        pieces.append(str(self.Primer3))
        return ''.join(pieces)
            
    def refresh(self):
        """Generates a new version of each module, incl. the random region."""
        self.RandomRegion.refresh()
        self.Model.GU = self.GU
        self.Model.refresh()

    def countMatches(self, verbose=False, temp=25):
        """Generates NumToDo sequences, folds them, and returns match count."""
        positions = []
        seqs = []
        structs = []
        orig_positions = self._fixed_positions
        self.Positions = orig_positions
        for i in xrange(self.NumToDo):
            self.refresh()
            if not orig_positions:
                self.Positions = self._choose_locations()
            curr_seq = str(self)
            #adjust positions to account for inserted modules
            curr_positions = self.Positions[:]
            insert_length = len(self.Primer5)
            module_lengths = map(len, list(self.Model))
            for i in range(len(curr_positions)):
                curr_positions[i] += insert_length
                insert_length += module_lengths[i]
            positions.append(curr_positions)
            seqs.append(curr_seq)
        folder = RNAfold(params={'-T':temp})
        struct_file = folder(seqs)['StdOut']
        odd = False
        for line in struct_file:
            if odd:
                structs.append(ViennaStructure(line.split()[0]))
            odd = not odd
        good_count = 0
        if self.Debug:
            print "DEBUGGING"   #debug code: prints seqs, structs, matches
        for seq, struct, position in zip(seqs, structs, positions):
            matched = self.Motif.structureMatches(struct, position, \
                self.MatchOffsets,debug=self.Debug)
            if matched:
                good_count += 1
            if self.Debug or (matched and self.ReportSeqs):
                module_lengths = map(len, list(self.Model))
                if self.Debug:
                    print "Module lengths:", module_lengths
                    print "Positions:", position
                print seq
                print struct
                temp = [' '] * len(seq)
                for l, p in zip(module_lengths, position):
                    temp[p:p+l] = ['*']*l
                print ''.join(temp)
                if self.Debug:
                    print "Offsets:", self.MatchOffsets
                print matched
        return good_count
