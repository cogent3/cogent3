#!/usr/bin/env python
"""Code for handling multiple sequence alignments. In particular:

    - SequenceCollection handles both aligned and unaligned sequences.
    - Alignment and its subclasses handle multiple sequence alignments, storing
      the raw sequences and a gap map. Useful for very long alignments, e.g.
      genomics data.
    - DenseAlignment and its subclasses handle multiple sequence alignments as
      arrays of characters. Especially useful for short alignments that contain
      many sequences.

    WARNING: The various alignment objects try to guess the input type from the
    input, but this behavior has a few quirks. In particular, if the input is
    a sequence of two-item sequences (e.g. a list of two-character strings), 
    each sequence will be unpacked and the first item will be used as the
    label, the second as the sequence. For example, Alignment(['AA','CC','AA'])
    produces an alignment of three 1-character strings labeled A, C and A
    respectively. The reason for this is that the common case is that you have
    passed in a stream of two-item label, sequence pairs. However, this can
    cause confusion when testing.
"""
from __future__ import division
from types import GeneratorType
from cogent.core.annotation import Map, _Annotatable
import cogent   #will use to get at cogent.parse.fasta.MinimalFastaParser,
                #which is a circular import otherwise.
from cogent.format.alignment import save_to_filename
from cogent.core.info import Info as InfoClass
from cogent.core.sequence import frac_same, ModelSequence
from cogent.maths.stats.util import Freqs
from cogent.format.fasta import fasta_from_alignment
from cogent.format.phylip import phylip_from_alignment
from cogent.format.nexus import nexus_from_alignment
from cogent.parse.gff import GffParser, parse_attributes
from numpy import nonzero, array, logical_or, logical_and, logical_not, \
    transpose, arange, zeros, ones, take, put, uint8, ndarray
from numpy.random import randint, permutation

from cogent.util.dict2d import Dict2D

from copy import copy
from cogent.core.profile import Profile

__author__ = "Peter Maxwell and Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Rob Knight", "Gavin Huttley",
               "Jeremy Widmann", "Catherine Lozupone", "Matthew Wakefield",
               "Micah Hamady", "Daniel McDonald", "Jan Kosinski"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class DataError(Exception):
    pass

eps = 1e-6  #small number: 1-eps is almost 1, and is used for things like the
            #default number of gaps to allow in a column.

def assign_sequential_names(ignored, num_seqs, base_name='seq', start_at=0):
    """Returns list of num_seqs sequential, unique names.
    
    First argument is ignored; expect this to be set as a class attribute.
    """
    return ['%s_%s' % (base_name,i) for i in range(start_at,start_at+num_seqs)]

class SeqLabeler(object):
    """Allows flexible seq labeling in toFasta()."""

    def __init__(self, aln, label_f=assign_sequential_names, **kwargs):
        """Initializes a new seq labeler."""
        self._aln = aln
        self._label_f = label_f
        self._map = dict(zip(aln.Names, label_f(len(aln.Names, **kwargs))))

    def __call__(self, s):
        """Returns seq name from seq id"""
        return self._map[s.Name]
    

def coerce_to_string(s):
    """Converts an arbitrary sequence into a string."""
    if isinstance(s, str):  #if it's a string, OK as is
        return s
    if isinstance(s, Aligned):  #if it's an Aligned object, convert to string
        return str(s)
    curr = str(s)           #if its string is the same length, return that
    if len(curr) == len(s):
        return curr
    try:
        return ''.join(s)   #assume it's a seq of chars
    except(TypeError, ValueError):
        return ''.join(map(str, s)) #general case (slow, might not be correct)

def seqs_from_array(a, Alphabet=None):
    """SequenceCollection from array of pos x seq: names are integers.
    
    This is an InputHandler for SequenceCollection. It converts an arbitrary
    array of numbers into Sequence objects using seq_constructor, and
    leaves the sequences unlabeled.
    """
    return list(transpose(a)), None

def seqs_from_model_seqs(seqs, Alphabet=None):
    """Alignment from ModelSequence objects: seqs -> array, names from seqs.
    
    This is an InputHandler for SequenceCollection. It converts a list of
    Sequence objects with _data and Name properties into a SequenceCollection
    that uses those sequences.
    """
    return seqs, [s.Name for s in seqs]

def seqs_from_generic(seqs, Alphabet=None):
    """SequenceCollection from generic seq x pos data: seq of seqs of chars.
    
    This is an InputHandler for SequenceCollection. It converts a generic list
    (each item in the list will be mapped onto an object using
    seq_constructor and assigns sequential integers (0-based) as names.
    """
    names = []
    for s in seqs:
        if hasattr(s, 'Name'):
            names.append(s.Name)
        else:
            names.append(None)
    return seqs, names

def seqs_from_fasta(seqs, Alphabet=None):
    """SequenceCollection from FASTA-format string or lines.
    
    This is an InputHandler for SequenceCollection. It converts a FASTA-format
    string or collection of lines into a SequenceCollection object, preserving
    order..
    """
    if isinstance(seqs, str):
        seqs = seqs.splitlines()
    names, seqs = zip(*list(cogent.parse.fasta.MinimalFastaParser(seqs)))
    return list(seqs), list(names)

def seqs_from_dict(seqs, Alphabet=None):
    """SequenceCollection from dict of {label:seq_as_str}.
    
    This is an InputHandler for SequenceCollection. It converts a dict in
    which the keys are the names and the values are the sequences
    (sequence only, no whitespace or other formatting) into a
    SequenceCollection. Because the dict doesn't preserve order, the result
    will not necessarily be in alphabetical order."""
    names, seqs = map(list, zip(*seqs.items()))
    return seqs, names

def seqs_from_kv_pairs(seqs, Alphabet=None):
    """SequenceCollection from list of (key, val) pairs.
    
    This is an InputHandler for SequenceCollection. It converts a dict in
    which the keys are the names and the values are the sequences
    (sequence only, no whitespace or other formatting) into a
    SequenceCollection. Because the dict doesn't preserve order, the result
    will be in arbitrary order."""
    names, seqs =  map(list, zip(*seqs))
    return seqs, names

def seqs_from_aln(seqs, Alphabet=None):
    """SequenceCollection from existing SequenceCollection object: copies data.
    
    This is relatively inefficient: you should really use the copy() method
    instead, which duplicates the internal data structures.
    """
    return seqs.Seqs, seqs.Names

def seqs_from_empty(obj, *args, **kwargs):
    """SequenceCollection from empty data: raise exception."""
    raise ValueError, "Cannot create empty SequenceCollection."

class SequenceCollection(object):
    """Base class for Alignment, but also just stores unaligned seqs.
    
    Handles shared functionality: detecting the input type, writing out the
    sequences as different formats, translating the sequences, chopping off
    stop codons, looking up sequences by name, etc.
    
    A SequenceCollection must support:
    
    - input handlers for different data types
    - SeqData: behaves like list of lists of chars, holds seq data
    - Seqs: behaves like list of Sequence objects, iterable in name order
    - Names: behaves like list of names for the sequence objects
    - NamedSeqs: behaves like dict of {name:seq}
    - MolType: specifies what kind of sequences are in the collection
    """
    InputHandlers = {   'array': seqs_from_array,
                        'model_seqs': seqs_from_model_seqs,
                        'generic': seqs_from_generic,
                        'fasta': seqs_from_fasta,
                        'collection': seqs_from_aln,
                        'aln': seqs_from_aln,
                        'dense_aln': seqs_from_aln,
                        'dict': seqs_from_dict,
                        'empty': seqs_from_empty,
                        'kv_pairs':seqs_from_kv_pairs,
                    }
    
    IsArray = set(['array', 'model_seqs'])
    
    DefaultNameFunction = assign_sequential_names
    
    def __init__(self, data, Names=None, Alphabet=None, MolType=None, \
            Name=None, Info=None, conversion_f=None, is_array=False, \
            force_same_data=False, \
            remove_duplicate_names=False, label_to_name=None,
            suppress_named_seqs=False):
        """Initialize self with data and optionally Info.
        
        We are always going to convert to characters, so Sequence objects
        in the collection will lose additional special attributes they have.
        This is somewhat inefficient, so it might be worth revisiting this
        decision later.
        
        The handling of sequence names requires special attention. Depending
        on the input data, we might get the names from the sequences themselves,
        or we might add them from Names that are passed in. However, the Names
        attribute controls the order that we examine the sequences in, so if
        it is passed in it should override the order that we got from the
        input data (e.g. you might pass in unlabeled sequences with the names
        ['b','a'], so that you want the first sequence to be called 'b' and
        the second to be called 'a', or you might pass in labeled sequences,
        e.g. as a dict, and the names ['b','a'], indicating that you want the
        sequence called b to be first and the sequence called a to be second
        despite the fact that they are in arbitrary order in the original
        input. In this second situation, it is imortant that the sequences not
        be relabeled.
        
        This is handled as followed. If the sequences are passed in using a
        method that does not carry the names with it, the Names that are passed
        in will be handed out to successive sequences. If the sequences are
        passed in using a method that does carry the names with it, the Names
        that are passed in will be used to order the sequences, but they will
        not be relabeled. Note that if you're passing in a data type that
        is already labeled (e.g. a list of Sequence objects) you _must_ have
        unique names beforehand.
        
        It's possible that this additional handling should be moved to a
        separate object; the motivation for having it on Alignment __init__
        is that it's easy for users to construct Alignment objects directly.
        
        Parameters:
        
        data:           Data to convert into a SequenceCollection
        
        Names:          Order of Names in the alignment. Should match the
                        names of the sequences (after processing by
                        label_to_name if present).
        
        Alphabet:       Alphabet to use for the alignment (primarily important
                        for DenseAlignment)
        
        MolType:        MolType to be applied to the Alignment and to each seq.
        
        Name:           Name of the SequenceCollection.
        
        Info:           Info object to be attached to the alignment itself.
        
        conversion_f:   Function to convert string into sequence.
        
        is_array:       True if input is an array, False otherwise.
        
        force_same_data: True if data will be used as the same object.
        
        remove_duplicate_names: True if duplicate names are to be silently
                        deleted instead of raising errors.
        
        label_to_name: if present, converts name into f(name).
        """
        
        #read all the data in if we were passed a generator
        if isinstance(data, GeneratorType):
            data = list(data)
        #set the Name
        self.Name = Name
        #figure out alphabet and moltype
        self.Alphabet, self.MolType = \
            self._get_alphabet_and_moltype(Alphabet, MolType, data)
        if not isinstance(Info, InfoClass):
            if Info:
                Info = InfoClass(Info)
            else:
                Info = InfoClass()
        self.Info = Info
        #if we're forcing the same data, skip the validation
        if force_same_data:
            self._force_same_data(data, Names)
            curr_seqs = data
        #otherwise, figure out what we got and coerce it into the right type
        else:
            per_seq_names, curr_seqs, name_order = \
                self._names_seqs_order(conversion_f, data, Names, is_array, \
                label_to_name, remove_duplicate_names, \
                Alphabet=self.Alphabet)
            self.Names = name_order
            
            #will take only the seqs and names that are in name_order
            if per_seq_names != name_order:
                good_indices = []
                for n in name_order:
                    good_indices.append(per_seq_names.index(n))
                if hasattr(curr_seqs, 'astype'): #it's an array
                    #much faster to check than to raise exception in this case
                    curr_seqs = take(curr_seqs, good_indices, axis=0)
                else:
                    curr_seqs = [curr_seqs[i] for i in good_indices]
                per_seq_names = name_order
            
            #create NamedSeqs dict for fast lookups
            if not suppress_named_seqs:
                self.NamedSeqs = self._make_named_seqs(self.Names, curr_seqs)
        #Sequence objects behave like sequences of chars, so no difference
        #between Seqs and SeqData. Note that this differs for Alignments,
        #so be careful which you use if writing methods that should work for
        #both SequenceCollections and Alignments.
        self._set_additional_attributes(curr_seqs)

    def __str__(self):
        """Returns self in FASTA-format, respecting name order."""
        return ''.join(['>%s\n%s\n' % (name, self.getGappedSeq(name))
                for name in self.Names])
    
    def _make_named_seqs(self, names, seqs):
        """Returns NamedSeqs: dict of name:seq."""
        name_seq_tuples = zip(names, seqs)
        for n, s in name_seq_tuples:
            s.Name = n
        return dict(name_seq_tuples)
    
    def _set_additional_attributes(self, curr_seqs):
        """Sets additional attributes based on current seqs: class-specific."""
        self.SeqData = curr_seqs
        self._seqs = curr_seqs
        try:
            self.SeqLen = max(map(len, curr_seqs))
        except ValueError: #got empty sequence, for some reason?
            self.SeqLen = 0
    
    def _force_same_data(self, data, Names):
        """Forces dict that was passed in to be used as self.NamedSeqs"""
        self.NamedSeqs = data
        self.Names = Names or data.keys()
    
    def copy(self):
        """Returns deep copy of self."""
        result = self.__class__(self, MolType=self.MolType, Info=self.Info)
        return result
    
    def _get_alphabet_and_moltype(self, Alphabet, MolType, data):
        """Returns Alphabet and MolType, giving MolType precedence."""
        if Alphabet is None and MolType is None:
            if hasattr(data, 'MolType'):
                MolType = data.MolType
            elif hasattr(data, 'Alphabet'):
                Alphabet = data.Alphabet
            #check for containers
            else:
                curr_item = self._get_container_item(data)
                if hasattr(curr_item, 'MolType'):
                    MolType = curr_item.MolType
                elif hasattr(curr_item, 'Alphabet'):
                    Alphabet = curr_item.Alphabet
                else:
                    MolType = self.MolType #will be BYTES by default
        if Alphabet is not None and MolType is None:
            MolType = Alphabet.MolType
        if MolType is not None and Alphabet is None:
            try:
                Alphabet = MolType.Alphabets.DegenGapped
            except AttributeError:
                Alphabet = MolType.Alphabet
        return Alphabet, MolType
    
    def _get_container_item(self, data):
        """Checks container for item with Alphabet or MolType"""
        curr_item = None
        if hasattr(data, 'itervalues'):
            curr_item = data.itervalues().next()
        else:
            try:
                curr_item = iter(data).next()
            except:
                pass
        return curr_item

    def _strip_duplicates(self, names, seqs):
        """Internal function to strip duplicates from list of names"""
        if len(set(names)) == len(names):
            return set(), names, seqs
        #if we got here, there are duplicates
        unique_names = {}
        duplicates = {}
        fixed_names = []
        fixed_seqs = []
        for n, s in zip(names, seqs):
            if n in unique_names:
                duplicates[n] = 1
            else:
                unique_names[n] = 1
                fixed_names.append(n)
                fixed_seqs.append(s)
        if type(seqs) is ndarray:
            fixed_seqs = array(fixed_seqs, seqs.dtype)
        return duplicates, fixed_names, fixed_seqs

    def _names_seqs_order(self, conversion_f, data, Names, is_array, \
            label_to_name, remove_duplicate_names, Alphabet=None):
        """Internal function to figure out names, seqs, and name_order."""
        #figure out conversion function and whether it's an array
        if not conversion_f:
            input_type = self._guess_input_type(data)
            is_array = input_type in self.IsArray
            conversion_f = self.InputHandlers[input_type]
        #set seqs and names as properties
        if Alphabet:
            seqs, names = conversion_f(data, Alphabet=Alphabet)
        else:
            seqs, names = conversion_f(data)
        if names and label_to_name:
            names = map(label_to_name, names)
        curr_seqs = self._coerce_seqs(seqs, is_array)

        #if no names were passed in as Names, if we obtained them from
        #the seqs we should use them, but otherwise we should use the
        #default names
        if Names is None:
            if (names is None) or (None in names):
                per_seq_names = name_order = \
                    self.DefaultNameFunction(len(curr_seqs))
            else:   #got names from seqs
                per_seq_names = name_order = names
        else:
            #otherwise, names were passed in as Names: use this as the order
            #if we got names from the sequences, but otherwise assign the
            #names to successive sequences in order
            if (names is None) or (None in names):
                per_seq_names = name_order = Names
            else: #got names from seqs, so assume name_order is in Names
                per_seq_names = names
                name_order = Names
        
        #check for duplicate names
        duplicates, fixed_names, fixed_seqs = \
            self._strip_duplicates(per_seq_names, curr_seqs)
        if duplicates:
            if remove_duplicate_names:
                per_seq_names, curr_seqs = fixed_names, fixed_seqs
                #if name_order doesn't have the same names as per_seq_names,
                #replace it with per_seq_names
                if (set(name_order) != set(per_seq_names)) or\
                    (len(name_order) != len(per_seq_names)):
                    name_order = per_seq_names
            else:
                raise ValueError, \
                "Some names were not unique. Duplicates are:\n" + \
                str(sorted(duplicates.keys()))

        return per_seq_names, curr_seqs, name_order
    
    def _coerce_seqs(self, seqs, is_array):
        """Controls how seqs are coerced in _names_seqs_order.
        
        Override in subclasses where this behavior should differ.
        """
        if is_array:
            seqs = map(str, map(self.MolType.ModelSeq, seqs))
        return map(self.MolType.Sequence, seqs)
    
    def _guess_input_type(self, data):
        """Guesses input type of data; returns result as key of InputHandlers.
        
        First checks whether data is an Alignment, then checks for some common
        string formats, then tries to do it based on string or array properties.
        
        Returns 'empty' if check fails, i.e. if it can't recognize the sequence
        as a specific type. Note that bad sequences are not guaranteed to
        return 'empty', and may be recognized as another type incorrectly.
        """
        if isinstance(data, DenseAlignment):
            return 'dense_aln'
        if isinstance(data, Alignment):
            return 'aln'
        if isinstance(data, SequenceCollection):
            return 'collection'
        if isinstance(data, dict):
            return 'dict'
        if isinstance(data, str):
            if data.startswith('>'):
                return 'fasta'
            else:
                return 'generic'
        first = None
        try:
            first = data[0]
        except (IndexError, TypeError):
            pass
        try:
            first = iter(data).next()
        except (IndexError, TypeError, StopIteration):
            pass
        if first is None:
            return 'empty'
        try:
            if isinstance(first, ModelSequence):     #model sequence base type
                return 'model_seqs'
            elif hasattr(first, 'dtype'):    #array object
                return 'array'
            elif isinstance(first, str) and first.startswith('>'):
                return 'fasta'
            else:
                try:
                    dict(data)
                    return 'kv_pairs'
                except (TypeError, ValueError):
                    pass
            return 'generic'
        except (IndexError, TypeError), e:
            return 'empty'
    
    def __cmp__(self, other):
        """cmp first tests as dict, then as str."""
        c = cmp(self.NamedSeqs, other)
        if not c:
            return 0
        else:
            return cmp(str(self), str(other))
    
    def keys(self):
        """keys uses self.Names, which defaults to known keys if None.
        
        Note: returns copy, not original.
        """
        return self.Names[:]
    
    def values(self):
        """values returns values corresponding to self.Names."""
        return [self.NamedSeqs[n] for n in self.Names]
    
    def items(self):
        """items returns (name, value) pairs."""
        return [(n, self.NamedSeqs[n]) for n in self.Names]
    
    def iterSeqs(self, seq_order=None):
        """Iterates over values (sequences) in the alignment, in order.
        
        seq_order: list of keys giving the order in which seqs will be returned.
        Defaults to self.Names. Note that only these sequences will be
        returned, and that KeyError will be raised if there are sequences
        in order that have been deleted from the Alignment. If self.Names
        is None, returns the sequences in the same order as
        self.NamedSeqs.values().
        
        Use map(f, self.seqs()) to apply the constructor f to each seq. f must
        accept a single list as an argument.
        
        Always returns references to the same objects that are values of the
        alignment.
        """
        ns = self.NamedSeqs
        get = ns.__getitem__
        for key  in seq_order or self.Names:
            yield get(key)
    
    def _take_seqs(self): return list(self.iterSeqs())
    
    Seqs = property(_take_seqs)   #access as attribute if using default order.
    
    def takeSeqs(self, seqs, negate=False, **kwargs):
        """Returns new Alignment containing only specified seqs.
        
        Note that the seqs in the new alignment will be references to the
        same objects as the seqs in the old alignment.
        """
        get = self.NamedSeqs.__getitem__
        result = {}
        if 'MolType' not in kwargs:
            kwargs['MolType'] = self.MolType
        
        if negate:
            #copy everything except the specified seqs
            negated_names = []
            row_lookup = dict.fromkeys(seqs)
            for r, row in self.NamedSeqs.items():
                if r not in row_lookup:
                    result[r] = row
                    negated_names.append(r)
            seqs = negated_names #remember to invert the list of names
        
        else:
            #copy only the specified seqs
            for r in seqs:
                result[r] = get(r)
        if result:
            return self.__class__(result, Names=seqs, **kwargs)
        else:
            return {}   #safe value; can't construct empty alignment
    
    def getSeqIndices(self, f, negate=False):
        """Returns list of keys of seqs where f(row) is True.
        
        List will be in the same order as self.Names, if present.
        """
        get = self.NamedSeqs.__getitem__
        #negate function if necessary
        if negate:
            new_f = lambda x: not f(x)
        else:
            new_f = f
        #get all the seqs where the function is True
        return [key for key in self.Names if new_f(get(key))]
    
    def takeSeqsIf(self, f, negate=False, **kwargs):
        """Returns new Alignment containing seqs where f(row) is True.
        
        Note that the seqs in the new Alignment are the same objects as the
        seqs in the old Alignment, not copies.
        """
        #pass negate to get SeqIndices
        return self.takeSeqs(self.getSeqIndices(f, negate), **kwargs)
    
    def iterItems(self, seq_order=None, pos_order=None):
        """Iterates over elements in the alignment.
        
        seq_order (names) can be used to select a subset of seqs.
        pos_order (positions) can be used to select a subset of positions.
        
        Always iterates along a seq first, then down a position (transposes
        normal order of a[i][j]; possibly, this should change)..
        
        WARNING: Alignment.iterItems() is not the same as alignment.iteritems()
        (which is the built-in dict iteritems that iterates over key-value
        pairs).
        """
        if pos_order:
            for row in self.iterSeqs(seq_order):
                for i in pos_order:
                    yield row[i]
        else:
            for row in self.iterSeqs(seq_order):
                for i in row:
                    yield i
    
    Items = property(iterItems)
    
    def getItems(self, items, negate=False):
        """Returns list containing only specified items.
        
        items should be a list of (row_key, col_key) tuples.
        """
        get = self.NamedSeqs.__getitem__
        if negate:
            #have to cycle through every item and check that it's not in
            #the list of items to return
            item_lookup = dict.fromkeys(map(tuple, items))
            result = []
            for r in self.Names:
                curr_row = get(r)
                for c in range(len(curr_row)):
                    if (r, c) not in items:
                        result.append(curr_row[c])
            return result
        #otherwise, just pick the selected items out of the list
        else:
            return [get(row)[col] for row, col in items]
    
    def getItemIndices(self, f, negate=False):
        """Returns list of (key,val) tuples where f(self.NamedSeqs[key][val])."""
        get = self.NamedSeqs.__getitem__
        if negate:
            new_f = lambda x: not f(x)
        else:
            new_f = f
        result = []
        for row_label in self.Names:
            curr_row = get(row_label)
            for col_idx, item in enumerate(curr_row):
                if new_f(item):
                    result.append((row_label, col_idx))
        return result
    
    def getItemsIf(self, f, negate=False):
        """Returns list of items where f(self.NamedSeqs[row][col]) is True."""
        return self.getItems(self.getItemIndices(f, negate))
    
    def getSimilar(self, target, min_similarity=0.0, max_similarity=1.0, \
        metric=frac_same, transform=None):
        """Returns new Alignment containing sequences similar to target.
        
        target: sequence object to compare to. Can be in the alignment.
        
        min_similarity: minimum similarity that will be kept. Default 0.0.
        
        max_similarity: maximum similarity that will be kept. Default 1.0.
        (Note that both min_similarity and max_similarity are inclusive.)
        
        metric: similarity function to use. Must be f(first_seq, second_seq).
        The default metric is fraction similarity, ranging from 0.0 (0%
        identical) to 1.0 (100% identical). The Sequence classes have lots
        of methods that can be passed in as unbound methods to act as the
        metric, e.g. fracSameGaps.
        
        transform: transformation function to use on the sequences before
        the metric is calculated. If None, uses the whole sequences in each
        case. A frequent transformation is a function that returns a specified
        range of a sequence, e.g. eliminating the ends. Note that the
        transform applies to both the real sequence and the target sequence.
        
        WARNING: if the transformation changes the type of the sequence (e.g.
        extracting a string from an RnaSequence object), distance metrics that
        depend on instance data of the original class may fail.
        """
        if transform:
            target = transform(target)
        m = lambda x: metric(target, x)
        
        if transform:
            def f(x):
                result = m(transform(x))
                return min_similarity <= result <= max_similarity
        else:
            def f(x):
                result = m(x)
                return min_similarity <= result <= max_similarity
        
        return self.takeSeqsIf(f)
    
    def distanceMatrix(self, f):
        """Returns Matrix containing pairwise distances between sequences.
        f is the distance function f(x,y) -> distance between x and y.
        
        It's often useful to pass an unbound method in as f.
        
        Does not assume that f(x,y) == f(y,x) or that f(x,x) == 0.
        """
        get = self.NamedSeqs.__getitem__
        seqs = self.NamedSeqs.keys()
        result = Dict2D()
        for i in seqs:
            for j in seqs:
                d = f(get(i), get(j))
                if i not in result:
                    result[i] = {}
                if j not in result:
                    result[j] = {}
                result[i][j] = d
                result[j][i] = d
        return result
    
    def isRagged(self):
        """Returns True if alignment has sequences of different lengths."""
        seqs = self.Seqs      #Get all sequences in alignment
        length = len(seqs[0])     #Get length of first sequence
        for seq in seqs:
            #If lengths differ
            if length != len(seq):
                return True
        #lengths were all equal
        return False
    
    def toPhylip(self, generic_label=True, make_seqlabel=None):
        """
        Return alignment in PHYLIP format and mapping to sequence ids
        
        raises exception if invalid alignment
        
        Arguments:
            - make_seqlabel: callback function that takes the seq object and
              returns a label str
        """
        return phylip_from_alignment(self, generic_label=generic_label,
                        make_seqlabel=make_seqlabel)
    
    def toFasta(self, make_seqlabel=None):
        """Return alignment in Fasta format
        
        Arguments:
            - make_seqlabel: callback function that takes the seq object and
              returns a label str
        """
        return fasta_from_alignment(self, make_seqlabel=make_seqlabel)
    
    def toNexus(self, seq_type, interleave_len=50):
        """
        Return alignment in NEXUS format and mapping to sequence ids
        
        **NOTE** Not that every sequence in the alignment MUST come from
            a different species!! (You can concatenate multiple sequences from
            same species together before building tree)
        
        seq_type: dna, rna, or protein
        
        Raises exception if invalid alignment
        """
        return nexus_from_alignment(self, seq_type,
                                interleave_len=interleave_len)
    
    def getIntMap(self,prefix='seq_'):
        """Returns a dict with names mapped to enumerates integer names.
            
            - prefix: prefix for sequence label. Default = 'seq_'
            - int_keys is a dict mapping int names to sorted original names.
        """
        get = self.NamedSeqs.__getitem__
        int_keys = dict([(prefix+str(i),k) for i,k in \
                enumerate(sorted(self.NamedSeqs.keys()))])
        int_map = dict([(k, copy(get(v))) for k,v in int_keys.items()])
        return int_map, int_keys
    
    def getNumSeqs(self):
        """Returns the number of sequences in the alignment."""
        return len(self.NamedSeqs)
    
    def copyAnnotations(self, unaligned):
        """Copies annotations from seqs in unaligned to self, matching by name. 
        
        Alignment programs like ClustalW don't preserve annotations,
        so this method is available to copy annotations off the unaligned
        sequences.  
        
        unaligned should be a dictionary of Sequence instances.
        
        Ignores sequences that are not in self, so safe to use on larger dict
        of seqs that are not in the current collection/alignment.
        """
        for name, seq in unaligned.items():
            if name in self.NamedSeqs:
                self.NamedSeqs[name].copyAnnotations(seq)
    
    def annotateFromGff(self, f):
        """Copies annotations from gff-format file to self.

        Matches by name of sequence. This method expects a file handle, not
        the name of a file.

        Skips sequences in the file that are not in self.
        """
        for (name, source, feature, start, end, score,
                strand, frame, attributes, comments) in GffParser(f):
            if name in self.NamedSeqs:
                self.NamedSeqs[name].addFeature(    feature, 
                                    parse_attributes(attributes), 
                                    [(start,end)])


            '''
            self.NamedSeqs[seqname].data.addFeature(
                                feature,
                                parse_attributes(attributes),
                                [(start, end)])
   ''' 
    def replaceSeqs(self, seqs):
        """Returns new alignment with same shape but with data taken from seqs.

        Primary use is for aligning codons from protein alignment, or, more
        generally, substituting in codons from a set of protein sequences (not
        necessarily aligned). For this reason, it takes characters from seqs
        three at a time rather than one at a time (i.e. 3 characters in seqs
        are put in place of 1 character in self).

        If seqs is an alignment, any gaps in it will be ignored.
        """
        if hasattr(seqs, 'NamedSeqs'):
            seqs = seqs.NamedSeqs
        else:
            seqs = SequenceCollection(seqs).NamedSeqs
        new_seqs = []
        for label in self.Names:
            aligned = self.NamedSeqs[label]
            seq = seqs[label]
            if isinstance(seq, Aligned):
                seq = seq.data
            new_seqs.append((label, Aligned(aligned.map * 3, seq)))
        return self.__class__(new_seqs)
    
    def getGappedSeq(self, seq_name, recode_gaps=False):
        """Return a gapped Sequence object for the specified seqname.

        Note: return type may depend on what data was loaded into the
        SequenceCollection or Alignment.
        """
        return self.NamedSeqs[seq_name]
  
    def __add__(self, other):
        """Concatenates sequence data for same names"""
        aligned = isinstance(self, Alignment)
        
        if len(self.NamedSeqs) != len(other.NamedSeqs):
            raise ValueError("Alignments don't have same number of sequences")
        
        concatenated = []
        for name in self.Names:
            if name not in other.Names:
                raise ValueError("Right alignment doesn't have a '%s'" % name)
            new_seq = self.NamedSeqs[name] + other.NamedSeqs[name]
            concatenated.append(new_seq)
        
        new = self.__class__(MolType=self.MolType,
                                data=zip(self.Names, concatenated))
        
        if aligned:
            left = [a for a in self._shiftedAnnotations(new, 0) \
                        if a.map.End <= len(self)]
            right = [a for a in other._shiftedAnnotations(new, len(self)) \
                        if a.map.Start >= len(self)]
            new.annotations = left + right
        return new
    
    def addSeqs(self, other, before_name=None, after_name=None):
        """Returns new object of class self with sequences from other added.
        
        By default the sequence is appended to the end of the alignment,
        this can be changed by using either before_name or after_name arguments.
        
        Arguments:
            - other: same class as self or coerceable to that class
            - before_name: str - [default:None] name of the sequence before
              which sequence is added
            - after_name: str - [default:None] name of the sequence after
              which sequence is added
              
        If both before_name and after_name are specified, the seqs will be
        inserted using before_name.
        """
        assert not isinstance(other, str), "Must provide a series of seqs "+\
                                            "or an alignment"
        self_seq_class = self.Seqs[0].__class__
        try:
            combined = self.Seqs + other.Seqs
        except AttributeError:
            combined = self.Seqs + list(other)
        
        for seq in combined:
            assert seq.__class__ == self_seq_class,\
                "Seq classes different: Expected %s, Got %s" % \
                    (seq.__class__, self_seq_class)
                    
        combined_aln = self.__class__(data=combined)
        
        if before_name is None and after_name is None:
            return combined_aln
        
        if (before_name and before_name not in self.Names) \
            or \
            (after_name and after_name not in self.Names):
            name = before_name or after_name
            raise ValueError("The alignment doesn't have a sequence named '{0}'" 
                .format(name))
        
        if before_name is not None: # someone might have seqname of int(0)
            index = self.Names.index(before_name)
        elif after_name is not None:
            index = self.Names.index(after_name) + 1
        
        names_before = self.Names[:index]
        names_after = self.Names[index:]
        new_names = combined_aln.Names[len(self.Names):]
        
        aln_new = combined_aln.takeSeqs(new_names)
        if len(names_before) > 0:
            aln_before = self.takeSeqs(names_before)
            combined_aln = aln_before
            combined_aln = combined_aln.addSeqs(aln_new)
        else:                
            combined_aln = aln_new
        
        if len(names_after) > 0:
            aln_after = self.takeSeqs(names_after)
            combined_aln = combined_aln.addSeqs(aln_after)
        
        return combined_aln
    
    def writeToFile(self, filename=None, format=None, **kwargs):
        """Write the alignment to a file, preserving order of sequences.
        
        Arguments:
        - filename: name of the sequence file
        - format: format of the sequence file
        
        If format is None, will attempt to infer format from the filename
        suffix.
        """
        
        if filename is None:
            raise DataError('no filename specified')
        
        # need to turn the alignment into a dictionary
        align_dict = {}
        for seq_name in self.Names:
            align_dict[seq_name] = str(self.NamedSeqs[seq_name])
        
        if format is None and '.' in filename:
            # allow extension to work if provided
            format = filename[filename.rfind(".")+1:]

        if 'order' not in kwargs:
            kwargs['order'] = self.Names
        save_to_filename(align_dict, filename, format, **kwargs)
    
    def __len__(self):
        """len of SequenceCollection returns length of longest sequence."""
        return self.SeqLen
    
    def getTranslation(self, gc=None, **kwargs):
        """Returns a new alignment object with the DNA sequences translated,
        using the current codon moltype, into an amino acid sequence.
        """
        translated = []
        aligned = isinstance(self, Alignment)
        # do the translation
        try:
            for seqname in self.Names:
                if aligned:
                    seq = self.getGappedSeq(seqname)
                else:
                    seq = self.NamedSeqs[seqname]
                pep = seq.getTranslation(gc)
                translated.append((seqname, pep))
            return self.__class__(translated, **kwargs)
        except AttributeError, msg:
            raise AttributeError, "%s -- %s" % (msg, "Did you set a DNA MolType?")
    
    def getSeq(self, seqname):
        """Return a sequence object for the specified seqname.
        """
        return self.NamedSeqs[seqname]
    
    def todict(self):
        """Returns the alignment as dict of names -> strings.

        Note: returns strings, NOT Sequence objects.
        """
        align_dict = {}
        
        for seq_name in self.Names:
            align_dict[seq_name] = str(self.NamedSeqs[seq_name])
        
        return align_dict
    
    def getPerSequenceAmbiguousPositions(self):
        """Returns dict of seq:{position:char} for ambiguous chars.
        
        Used in likelihood calculations.
        """
        result = {}
        for name in self.Names:
            result[name] = ambig = {}
            for (i, motif) in enumerate(self.getGappedSeq(name)):
                if self.MolType.isAmbiguity(motif):
                    ambig[i] = motif
        return result
    
    def degap(self, **kwargs):
        """Returns copy in which sequences have no gaps."""
        new_seqs = []
        aligned = isinstance(self, Alignment)
        for seq_name in self.Names:
            if aligned:
                seq = self.NamedSeqs[seq_name].data
            else:
                seq = self.NamedSeqs[seq_name]
            new_seqs.append((seq_name, seq.degap()))
        return SequenceCollection(MolType=self.MolType, data=new_seqs, **kwargs)
    
    def withModifiedTermini(self):
        """Changes the termini to include termini char instead of gapmotif.
        
        Useful to correct the standard gap char output by most
        alignment programs when aligned sequences have different ends.
        """
        seqs = []
        for name in self.Names:
            seq = self.NamedSeqs[name].withTerminiUnknown()
            seqs.append((name, seq))
        return self.__class__(MolType=self.MolType, data=seqs)
    
    def hasTerminalStops(self, gc=None):
        """Returns True if any sequence has a terminal stop codon."""
        stops = []
        aligned = isinstance(self, Alignment)
        for seq_name in self.Names:
            if aligned:
                seq = self.NamedSeqs[seq_name].data
            else:
                seq = self.NamedSeqs[seq_name]
            stops.append(seq.hasTerminalStop(gc=gc))
        return max(stops)
    
    def withoutTerminalStopCodons(self, gc=None, **kwargs):
        """Removes any terminal stop codons from the sequences"""
        new_seqs = []
        aligned = isinstance(self, Alignment)
        for seq_name in self.Names:
            old_seq = self.NamedSeqs[seq_name]
            if aligned:
                new_seq = old_seq.data.withoutTerminalStopCodon(gc=gc)
                new_seq = Aligned(old_seq.map, new_seq)
            else:
                new_seq = old_seq.withoutTerminalStopCodon(gc=gc)
            new_seqs.append((seq_name, new_seq))
        return self.__class__(MolType=self.MolType, data=new_seqs, **kwargs)
    
    def getSeqNames(self):
        """Return a list of sequence names."""
        return self.Names[:]
    
    def getMotifProbs(self, alphabet=None, include_ambiguity=False,
            exclude_unobserved=False, allow_gap=False, pseudocount=0):
        """Return a dictionary of motif probs.
        
        Arguments:
            - include_ambiguity: if True resolved ambiguous codes are
              included in estimation of frequencies, default is False.
            - exclude_unobserved: if True, motifs that are not present in
              the alignment are excluded from the returned dictionary,
              default is False.
            - allow_gap: allow gap motif
        """
        if alphabet is None:
            alphabet = self.MolType.Alphabet
            if allow_gap:
                alphabet = alphabet.Gapped
        
        counts = {}
        for seq_name in self.Names:
            sequence = self.NamedSeqs[seq_name]
            motif_len = alphabet.getMotifLen()
            if motif_len > 1:
                posns = range(0, len(sequence)+1-motif_len, motif_len)
                sequence = [sequence[i:i+motif_len] for i in posns]
            for motif in sequence:
                if not allow_gap:
                    if self.MolType.Gap in motif:
                        continue
                
                if motif in counts:
                    counts[motif] += 1
                else:
                    counts[motif] = 1
        
        probs = {}
        if not exclude_unobserved:
            for motif in alphabet:
                probs[motif] = pseudocount
        
        for (motif, count) in counts.items():
            motif_set = alphabet.resolveAmbiguity(motif)
            if len(motif_set) > 1:
                if include_ambiguity:
                    count = float(count) / len(motif_set)
                else:
                    continue
            for motif in motif_set:
                probs[motif] = probs.get(motif, pseudocount) + count
        
        total = float(sum(probs.values()))
        for motif in probs:
            probs[motif] /= total
        
        return probs
    
    def getGapCount(self, seq_name):
        return len(self.NamedSeqs[seq_name].map.gaps())
    
    def getSeqFreqs(self):
        """Returns Profile of counts: seq by character.
        
        See documentation for _get_freqs: this just wraps it and converts the
        result into a Profile object organized per-sequence (i.e. per row).
        """
        return Profile(self._get_freqs(0), self.Alphabet)
    
    def _make_gaps_ok(self, allowed_gap_frac):
        """Makes the gaps_ok function used by omitGapPositions and omitGapSeqs.
        
        Need to make the function because if it's a method of Alignment, it
        has unwanted 'self' and 'allowed_gap_frac' parameters that impede the
        use of map() in takeSeqsIf.
        
        WARNING: may not work correctly if component sequences have gaps that
        are not the Alignment gap character. This is because the gaps are
        checked at the position level (and the positions are lists), rather than
        at the sequence level. Working around this issue would probably cause a
        significant speed penalty.
        """
        def gaps_ok(seq):
            seq_len = len(seq)
            try:
                num_gaps = seq.countGaps()
            except AttributeError:
                num_gaps = len(filter(self.MolType.Gaps.__contains__, seq))
            return num_gaps / seq_len <= allowed_gap_frac
        
        return gaps_ok
    
    def omitGapPositions(self, allowed_gap_frac=1-eps, del_seqs=False, \
        allowed_frac_bad_cols=0, seq_constructor=None):
        """Returns new alignment where all cols have <= allowed_gap_frac gaps.
        
        allowed_gap_frac says what proportion of gaps is allowed in each
        column (default is 1-eps, i.e. all cols with at least one non-gap
        character are preserved).
        
        If del_seqs is True (default:False), deletes the sequences that don't
        have gaps where everything else does. Otherwise, just deletes the
        corresponding column from all sequences, in which case real data as
        well as gaps can be removed.
        
        Uses seq_constructor(seq) to make each new sequence object.
        
        Note: a sequence that is all gaps will not be deleted by del_seqs
        (even if all the positions have been deleted), since it has no non-gaps
        in positions that are being deleted for their gap content. Possibly,
        this decision should be revisited since it may be a surprising
        result (and there are more convenient ways to return the sequences
        that consist wholly of gaps).
        """
        if seq_constructor is None:
            seq_constructor = self.MolType.Sequence
        gaps_ok = self._make_gaps_ok(allowed_gap_frac)
        #if we're not deleting the 'naughty' seqs that contribute to the
        #gaps, it's easy...
        if not del_seqs:
            return self.takePositionsIf(f=gaps_ok, \
                seq_constructor=seq_constructor)
        #otherwise, we have to figure out which seqs to delete.
        #if we get here, we're doing del_seqs.
        cols_to_delete = dict.fromkeys(self.getPositionIndices(gaps_ok, \
            negate=True))
        default_gap_f = self.MolType.Gaps.__contains__
        
        bad_cols_per_row = {}
        for key, row in self.NamedSeqs.items():
            try:
                is_gap = row.Alphabet.Gaps.__contains__
            except AttributeError:
                is_gap = default_gap_f
            
            for col in cols_to_delete:
                if not is_gap(str(row)[col]):
                    if key not in bad_cols_per_row:
                        bad_cols_per_row[key] = 1
                    else:
                        bad_cols_per_row[key] += 1
        #figure out which of the seqs we're deleting
        get = self.NamedSeqs.__getitem__
        seqs_to_delete = {}
        for key, count in bad_cols_per_row.items():
            if float(count)/len(get(key)) >= allowed_frac_bad_cols:
                seqs_to_delete[key] = True
        #It's _much_ more efficient to delete the seqs before the cols.
        good_seqs = self.takeSeqs(seqs_to_delete, negate=True)
        cols_to_keep = dict.fromkeys(range(self.SeqLen))
        for c in cols_to_delete:
            del cols_to_keep[c]
        if good_seqs:
            return good_seqs.takePositions(cols=cols_to_keep.keys(), \
                seq_constructor=seq_constructor)
        else:
            return {}
    
    def omitGapSeqs(self, allowed_gap_frac=0):
        """Returns new alignment with seqs that have <= allowed_gap_frac.
        
        allowed_gap_frac should be a fraction between 0 and 1 inclusive.
        Default is 0.
        """
        gaps_ok = self._make_gaps_ok(allowed_gap_frac)
        
        return self.takeSeqsIf(gaps_ok)
    
    def omitGapRuns(self, allowed_run=1):
        """Returns new alignment where all seqs have runs of gaps <=allowed_run.
        
        Note that seqs with exactly allowed_run gaps are not deleted.
        Default is for allowed_run to be 1 (i.e. no consecutive gaps allowed).
        
        Because the test for whether the current gap run exceeds the maximum
        allowed gap run is only triggered when there is at least one gap, even
        negative values for allowed_run will still let sequences with no gaps
        through.
        """
        def ok_gap_run(x):
            try:
                is_gap = x.Alphabet.Gaps.__contains__
            except AttributeError:
                is_gap = self.MolType.Gaps.__contains__
            curr_run = max_run = 0
            for i in x:
                if is_gap(i):
                    curr_run += 1
                    if curr_run > allowed_run:
                        return False
                else:
                    curr_run = 0
            #can only get here if max_run was never exceeded (although this
            #does include the case where the sequence is empty)
            return True
        
        return self.takeSeqsIf(ok_gap_run)
    
    def omitSeqsTemplate(self, template_name, gap_fraction, gap_run):
        """Returns new alignment where all seqs are well aligned with template.
        
        gap_fraction = fraction of positions that either have a gap in the
            template but not in the seq or in the seq but not in the template
        gap_run = number of consecutive gaps tolerated in query relative to
            sequence or sequence relative to query
        """
        template = self.NamedSeqs[template_name]
        gap_filter = make_gap_filter(template, gap_fraction, gap_run)
        return self.takeSeqsIf(gap_filter)
    
    def toDna(self):
        """Returns the alignment as DNA."""
        seqs = [self.NamedSeqs[name].toDna() for name in self.Names]
        aln = self.__class__(data=seqs, Names=self.Names[:], Name=self.Name, Info=self.Info)
        if isinstance(self, _Annotatable) and self.annotations:
            aln.annotations = self.annotations[:]
        return aln
    
    def toRna(self):
        """Returns the alignment as RNA"""
        seqs = [self.NamedSeqs[name].toRna() for name in self.Names]
        aln = self.__class__(data=seqs, Names=self.Names[:], Name=self.Name, Info=self.Info)
        if isinstance(self, _Annotatable) and self.annotations:
            aln.annotations = self.annotations[:]
        return aln
    
    def rc(self):
        """Returns the reverse complement alignment"""
        seqs = [self.NamedSeqs[name].rc() for name in self.Names]
        rc = self.__class__(data=seqs, Names=self.Names[:], Name=self.Name, Info=self.Info)
        if isinstance(self, _Annotatable) and self.annotations:
            self._annotations_nucleic_reversed_on(rc)
        return rc
    
    def reversecomplement(self):
        """Returns the reverse complement alignment. A synonymn for rc."""
        return self.rc()
    
    def padSeqs(self, pad_length=None, **kwargs):
        """Returns copy in which sequences are padded to same length.
            
            pad_length: Length all sequences are to be padded to.  Will pad
                to max sequence length if pad_length is None or less than max
                length.
        """
        #get max length
        max_len = max([len(s) for s in self.Seqs])
        #If a pad_length was passed in, make sure it is valid
        if pad_length is not None:
            pad_length = int(pad_length)
            if pad_length < max_len:
                raise ValueError, \
        "pad_length must be at greater or equal to maximum sequence length: %s"\
                    %(str(max_len))
        #pad_length is max sequence length.
        else:
            pad_length = max_len
            
        #Get new sequence list
        new_seqs = []
        aligned = isinstance(self, Alignment)
        #for each sequence, pad gaps to end
        for seq_name in self.Names:
            if aligned:
                seq = self.NamedSeqs[seq_name].data
            else:
                seq = self.NamedSeqs[seq_name]
            padded_seq = seq + '-'*(pad_length-len(seq))
            new_seqs.append((seq_name, padded_seq))
            
        #return new SequenceCollection object
        return SequenceCollection(MolType=self.MolType, data=new_seqs, **kwargs)
        
        
class Aligned(object):
    """One sequence in an alignment, a map between alignment coordinates and
    sequence coordinates"""
    
    def __init__(self, map, data, length=None):
        #Unlike the normal map constructor, here we take a list of pairs of
        #alignment coordinates, NOT a list of pairs of sequence coordinates
        if isinstance(map, list):
            map = Map(map, parent_length=length).inverse()
        self.map = map
        self.data = data
        if hasattr(data, 'Info'):
            self.Info = data.Info
        if hasattr(data, 'Name'):
            self.Name = data.Name

    def _get_moltype(self):
        return self.data.MolType
    MolType = property(_get_moltype)
    
    def copy(self, memo=None, _nil=[], constructor='ignored'):
        """Returns a shallow copy of self

        WARNING: cogent.core.sequence.Sequence does NOT implement a copy method,
        as such, the data member variable of the copied object will maintain
        reference to the original object.

        WARNING: cogent.core.location.Map does NOT implement a copy method, as 
        such, the data member variable of the copied object will maintain
        reference to the original object.
        """
        return self.__class__(self.map, self.data)

    def __repr__(self):
        return '%s of %s' % (repr(self.map), repr(self.data))
    
    def withTerminiUnknown(self):
        return self.__class__(self.map.withTerminiUnknown(), self.data)
    
    def copyAnnotations(self, other):
        self.data.copyAnnotations(other)
    
    def annotateFromGff(self, f):
        self.data.annotate_from_gff(f)

    def addFeature(self, *args, **kwargs):
        self.data.addFeature(*args, **kwargs)
    
    def __str__(self):
        """Returns string representation of aligned sequence, incl. gaps."""
        return str(self.getGappedSeq())
    
    def __cmp__(self, other):
        """Compares based on string representations."""
        return cmp(str(self), str(other))
    
    def __iter__(self):
        """Iterates over sequence one motif (e.g. char) at a time, incl. gaps"""
        return self.data.gappedByMapMotifIter(self.map)
    
    def getGappedSeq(self, recode_gaps=False):
        """Returns sequence as an object, including gaps."""
        return self.data.gappedByMap(self.map, recode_gaps)
    
    def __len__(self):
        # these make it look like Aligned should be a subclass of Map,
        # but then you have to be careful with __getitem__, __init__ and inverse.
        return len(self.map)
    
    def __add__(self, other):
        if self.data is other.data:
            (map, seq) = (self.map + other.map, self.data)
        else:
            seq = self.getGappedSeq() + other.getGappedSeq()
            (map, seq) = seq.parseOutGaps()
        return Aligned(map, seq)
            
    def __getitem__(self, slice):
        return Aligned(self.map[slice], self.data)
    
    def rc(self):
        return Aligned(self.map.reversed(), self.data)
        
    def toRna(self):
        return Aligned(self.map, self.data.toRna())
        
    def toDna(self):
        return Aligned(self.map, self.data.toDna())
        
    def getTracks(self, policy):
        policy = policy.at(self.map.inverse())
        return self.data.getTracks(policy)
    
    def remappedTo(self, map):
        #assert map is self.parent_map or ... ?
        #print 'REMAP', self.map, self
        #print 'ONTO', map, map.inverse()
        result = Aligned(map[self.map.inverse()].inverse(), self.data)
        #print 'GIVES', result.map, result
        #print
        return result
    
    def getAnnotationsMatching(self, alignment, *args):
        for annot in self.data.getAnnotationsMatching(*args):
            yield annot.remappedTo(alignment, self.map.inverse())
    
    def gapVector(self):
        """Returns gapVector of GappedSeq, for omitGapPositions."""
        return self.getGappedSeq().gapVector()
    
    def _masked_annotations(self, annot_types, mask_char, shadow):
        """returns a new aligned sequence with regions defined by align_spans
        and shadow masked."""
        new_data = self.data.withMaskedAnnotations(annot_types, mask_char, shadow)
        # we remove the mask annotations from self and new_data
        return self.__class__(self.map, new_data)
    

class AlignmentI(object):
    """Alignment interface object. Contains methods shared by implementations.
    
    Note that subclasses should inherit both from AlignmentI and from
    SequenceCollection (typically).
    
    Alignments are expected to be immutable once created. No mechanism is
    provided for maintaining reference consistency if data in the alignment
    are modified.
    
    An Alignment is expected to be able to generate the following:
    - Seqs:         Sequence objects in the alignment, can turn themselves into
                    strings. These are usually thought of as "rows" in an
                    alignment.
    - Positions:    Vectors representing data in each position in the alignment
                    These are usually thought of as "columns" in an alignment.
    - SeqData:      Vectors representing data in each sequence in the alignment,
                    not necessarily guaranteed to turn themselves into a string
    - Items:        Iterator over the characters in the alignment
    - Names:        List of names of sequences in the alignment. Used for
                    display order. A cheap way to omit or reorder sequences is
                    to modify the list of names.
    - NamedSeqs:    Dict of name -> seq object, used for lookup.
    - MolType:      MolType of the alignment.
    """
    DefaultGap = '-'                #default gap character for padding
    GapChars = dict.fromkeys('-?')  #default gap chars for comparisons
    
    def iterPositions(self, pos_order=None):
        """Iterates over positions in the alignment, in order.
        
        pos_order refers to a list of indices (ints) specifying the column
        order. This lets you rearrange positions if you want to (e.g. to pull
        out individual codon positions).
        
        Note that self.iterPositions() always returns new objects, by default
        lists of elements. Use map(f, self.iterPositions) to apply the
        constructor or function f to the resulting lists (f must take a single
        list as a parameter). Note that some sequences (e.g. ViennaStructures)
        have rules that prevent arbitrary strings of their symbols from being
        valid objects.
        
        Will raise IndexError if one of the indices in order exceeds the
        sequence length. This will always happen on ragged alignments:
        assign to self.SeqLen to set all sequences to the same length.
        """
        get = self.NamedSeqs.__getitem__
        pos_order = pos_order or xrange(self.SeqLen)
        seq_order = self.Names
        for pos in pos_order:
            yield [get(seq)[pos] for seq in seq_order]
    
    Positions = property(iterPositions)
    
    def takePositions(self, cols, negate=False, seq_constructor=None):
        """Returns new Alignment containing only specified positions.
        
        By default, the seqs will be lists, but an alternative constructor
        can be specified.
        
        Note that takePositions will fail on ragged positions.
        """
        if seq_constructor is None:
            seq_constructor = self.MolType.Sequence
        result = {}
        #if we're negating, pick out all the positions except specified indices
        if negate:
            col_lookup = dict.fromkeys(cols)
            for key, row in self.NamedSeqs.items():
                result[key] = seq_constructor([row[i] for i in range(len(row)) \
                if i not in col_lookup])
        #otherwise, just get the requested indices
        else:
            for key, row in self.NamedSeqs.items():
                result[key] = seq_constructor([row[i] for i in cols])
        return self.__class__(result, Names=self.Names)
    
    def getPositionIndices(self, f, negate=False):
        """Returns list of column indices for which f(col) is True."""
        #negate f if necessary
        if negate:
            new_f = lambda x: not f(x)
        else:
            new_f = f
        return [i for i, col in enumerate(self.Positions) if new_f(col)]
    
    def takePositionsIf(self, f, negate=False, seq_constructor=None):
        """Returns new Alignment containing cols where f(col) is True.
        
        Note that the seqs in the new Alignment are always new objects. Default
        constructor is list(), but an alternative can be passed in.
        """
        if seq_constructor is None:
            seq_constructor = self.MolType.Sequence
        return self.takePositions(self.getPositionIndices(f, negate), \
            seq_constructor=seq_constructor)
    
    def IUPACConsensus(self, alphabet=None):
        """Returns string containing IUPAC consensus sequence of the alignment.
        """
        if alphabet is None:
            alphabet = self.MolType
        consensus = []
        degen = alphabet.degenerateFromSequence
        for col in self.Positions:
            consensus.append(degen(coerce_to_string(col)))
        return coerce_to_string(consensus)
    
    def columnFreqs(self, constructor=Freqs):
        """Returns list of Freqs with item counts for each column.
        """
        return map(constructor, self.Positions)
    
    def columnProbs(self, constructor=Freqs):
        """Returns FrequencyDistribuutions w/ prob. of each item per column.
        
        Implemented as a list of normalized Freqs objects.
        """
        freqs = self.columnFreqs(constructor)
        
        for fd in freqs:
            fd.normalize()
        return freqs
    
    def majorityConsensus(self, transform=None, constructor=Freqs):
        """Returns list containing most frequent item at each position.
        
        Optional parameter transform gives constructor for type to which result
        will be converted (useful when consensus should be same type as
        originals).
        """
        col_freqs = self.columnFreqs(constructor)
        
        consensus = [freq.Mode for freq in col_freqs]
        if transform == str:
            return coerce_to_string(consensus)
        elif transform:
            return transform(consensus)
        else:
            return consensus
    
    def uncertainties(self, good_items=None):
        """Returns Shannon uncertainty at each position.
        
        Usage: information_list = alignment.information(good_items=None)
        
        If good_items is supplied, deletes any symbols that are not in
        good_items.
        """
        uncertainties = []
        #calculate column probabilities if necessary
        if hasattr(self, 'PositionumnProbs'):
            probs = self.PositionumnProbs
        else:
            probs = self.columnProbs()
        #calculate uncertainty for each column
        for prob in probs:
            #if there's a list of valid symbols, need to delete everything else
            if good_items:
                prob = prob.copy()  #do not change original
                #get rid of any symbols not in good_items
                for symbol in prob.keys():
                    if symbol not in good_items:
                        del prob[symbol]
                #normalize the probabilities and add to the list
                prob.normalize()
            uncertainties.append(prob.Uncertainty)
        return uncertainties
    
    def scoreMatrix(self):
        """Returns a position specific score matrix for the alignment."""
        return Dict2D(dict([(i,Freqs(col)) for i, col in enumerate(self.Positions)]))
    
    def _get_freqs(self, index=None):
        """Gets array of freqs along index 0 (= positions) or 1 (= seqs).
        
        index: if 0, will calculate the frequency of each symbol in each
        position (=column) in the alignment. Will return 2D array where the
        first index is the position, and the second index is the index of the
        symbol in the alphabet. For example, for the TCAG DNA Alphabet,
        result[3][0] would store the count of T at position 3 (i.e. the 4th
        position in the alignment.
        
        if 1, does the same thing except that the calculation is performed for
        each sequence, so the 2D array has the sequence index as the first
        index, and the symbol index as the second index. For example, for the
        TCAG DNA Alphabet, result[3][0] would store the count of T in the
        sequence at index 3 (i.e. the 4th sequence).
        
        First an DenseAligment object is created, next the calculation is done
        on this object. It is important that the DenseAlignment is initialized
        with the same MolType and Alphabet as the original Alignment.
        """
        da = DenseAlignment(self, MolType=self.MolType, Alphabet=self.Alphabet)
        return da._get_freqs(index)

    def getPosFreqs(self):
        """Returns Profile of counts: position by character.
        
        See documentation for _get_freqs: this just wraps it and converts the
        result into a Profile object organized per-position (i.e. per column).
        """
        return Profile(self._get_freqs(1), self.Alphabet)

    def sample(self, n=None, with_replacement=False, motif_length=1, \
        randint=randint, permutation=permutation):
        """Returns random sample of positions from self, e.g. to bootstrap.

        Arguments:
            - n: the number of positions to sample from the alignment.
              Default is alignment length
            - with_replacement: boolean flag for determining if sampled
              positions
            - random_series: a random number generator with
              .randint(min,max) .random() methods
              
            
        Notes: 
            By default (resampling all positions without replacement), generates
            a permutation of the positions of the alignment.

            Setting with_replacement to True and otherwise leaving parameters
            as defaults generates a standard bootstrap resampling of the 
            alignment.
            """
        population_size = len(self) // motif_length
        if not n:
            n = population_size
        if with_replacement:
            locations = randint(0, population_size, n)
        else:
            assert n <= population_size, (n, population_size, motif_length)
            locations = permutation(population_size)[:n]
        positions = [(loc*motif_length, (loc+1)*motif_length)
                for loc in locations]
        sample = Map(positions, parent_length=len(self))
        return self.gappedByMap(sample, Info=self.Info)

    def slidingWindows(self, window, step, start=None, end=None):
        """Generator yielding new Alignments of given length and interval.

        Arguments:
            - window: The length of each returned alignment.
            - step: The interval between the start of the successive
              alignment objects returned.
            - start: first window start position
            - end: last window start position
        """
        start = [start, 0][start is None]
        end = [end, len(self)-window+1][end is None]
        end = min(len(self)-window+1, end)
        if start < end and len(self)-end >= window-1:
            for pos in xrange(start, end, step):
                yield self[pos:pos+window]
        
    
  
def aln_from_array(a, array_type=None, Alphabet=None):
    """Alignment from array of pos x seq: no change, names are integers.
    
    This is an InputHandler for Alignment. It converts an arbitrary array
    of numbers without change, but adds successive integer names (0-based) to
    each sequence (i.e. column) in the input a. Data type of input is
    unchanged.
    """
    if array_type is None:
        result = a.copy()
    else:
        result = a.astype(array_type)
    return transpose(result), None

def aln_from_model_seqs(seqs, array_type=None, Alphabet=None):
    """Alignment from ModelSequence objects: seqs -> array, names from seqs.
    
    This is an InputHandler for Alignment. It converts a list of Sequence
    objects with _data and Label properties into the character array Alignment
    needs. All sequences must be the same length.
    
    WARNING: Assumes that the ModelSeqs are already in the right alphabet. If
    this is not the case, e.g. if you are putting sequences on a degenerate
    alphabet into a non-degenerate alignment or you are putting protein
    sequences into a DNA alignment, there will be problems with the alphabet
    mapping (i.e. the resulting sequences may be meaningless).
    
    WARNING: Data type of return array is not guaranteed -- check in caller!
    """
    data, names = [], []
    for s in seqs:
        data.append(s._data)
        names.append(s.Name)
    result = array(data)
    if array_type:
        result = result.astype(array_type)
    return result, names

def aln_from_generic(data, array_type=None, Alphabet=None):
    """Alignment from generic seq x pos data: sequence of sequences of chars.
    
    This is an InputHandler for Alignment. It converts a generic list (each
    item in the list will be mapped onto an Array object, with character
    transformations, all items must be the same length) into a numpy array,
    and assigns sequential integers (0-based) as names.
    
    WARNING: Data type of return array is not guaranteed -- check in caller!
    """
    result = array(map(Alphabet.toIndices, data))
    names = []
    for d in data:
        if hasattr(d, 'Name'):
            names.append(d.Name)
        else:
            names.append(None)
    if array_type:
        result = result.astype(array_type)
    return result, names

def aln_from_collection(seqs, array_type=None, Alphabet=None):
    """Alignment from SequenceCollection object, or its subclasses."""
    names = seqs.Names
    data = [seqs.NamedSeqs[i] for i in names]
    result = array(map(Alphabet.toIndices, data))
    if array_type:
        result = result.astype(array_type)
    return result, names

def aln_from_fasta(seqs, array_type=None, Alphabet=None):
    """Alignment from FASTA-format string or lines.
    
    This is an InputHandler for Alignment. It converts a FASTA-format string
    or collection of lines into an Alignment object. All sequences must be the
    same length.
    
    WARNING: Data type of return array is not guaranteed -- check in caller!
    """
    if isinstance(seqs, str):
        seqs = seqs.splitlines()
    return aln_from_model_seqs([ModelSequence(s, Name=l, Alphabet=Alphabet)\
        for l, s in cogent.parse.fasta.MinimalFastaParser(seqs)], array_type)

def aln_from_dict(aln, array_type=None, Alphabet=None):
    """Alignment from dict of {label:seq_as_str}.
    
    This is an InputHandler for Alignment. It converts a dict in which the
    keys are the names and the values are the sequences (sequence only, no
    whitespace or other formatting) into an alignment. Because the dict
    doesn't preserve order, the result will be in alphabetical order."""
    names, seqs = zip(*sorted(aln.items()))
    result = array(map(Alphabet.toIndices, seqs), array_type)
    return result, list(names)

def aln_from_kv_pairs(aln, array_type=None, Alphabet=None):
    """Alignment from sequence of (key, value) pairs.
    
    This is an InputHandler for Alignment. It converts a list in which the
    first item of each pair is the label and the second item is the sequence
    (sequence only, no whitespace or other formatting) into an alignment.
    Because the dict doesn't preserve order, the result will be in arbitrary
    order."""
    names, seqs = zip(*aln)
    result = array(map(Alphabet.toIndices, seqs), array_type)
    return result, list(names)

def aln_from_dense_aln(aln, array_type=None, Alphabet=None):
    """Alignment from existing DenseAlignment object: copies data.
    
    Retrieves data from Positions field. Uses copy(), so array data type
    should be unchanged.
    """
    if array_type is None:
        result = aln.ArrayPositions.copy()
    else:
        result = aln.ArrayPositions.astype(array_type)
    return transpose(result), aln.Names[:]

def aln_from_empty(obj, *args, **kwargs):
    """Alignment from empty data: raise exception."""
    raise ValueError, "Cannot create empty alignment."

#Implementation of Alignment base class

class DenseAlignment(AlignmentI, SequenceCollection):
    """Holds a dense array representing a multiple sequence alignment.
    
    An Alignment is _often_, but not necessarily, an array of chars. You might
    want to use some other data type for the alignment if you have a large
    number of symbols. For example, codons on an ungapped DNA alphabet has
    4*4*4=64 entries so can fit in a standard char data type, but tripeptides
    on the 20-letter ungapped protein alphabet has 20*20*20=8000 entries so
    can _not_ fit in a char and values will wrap around (i.e. you will get an
    unpredictable, wrong value for any item whose index is greater than the
    max value, e.g. 255 for uint8), so in this case you would need to use
    UInt16, which can hold 65536 values. DO NOT USE SIGNED DATA TYPES FOR YOUR
    ALIGNMENT ARRAY UNLESS YOU LOVE MISERY AND HARD-TO-DEBUG PROBLEMS.
    
    Implementation: aln[i] returns position i in the alignment.
    
    aln.Positions[i] returns the same as aln[i] -- usually, users think of this
    as a 'column', because alignment editors such as Clustal typically display
    each sequence as a row so a position that cuts across sequences is a
    column.
    
    aln.Seqs[i] returns a sequence, or 'row' of the alignment in standard
    terminology.
    
    WARNING: aln.Seqs and aln.Positions are different views of the same array,
    so if you change one you will change the other. This will no longer be
    true if you assign to Seqs or Positions directly, so don't do it. If you
    want to change the data in the whole array, always assign to a slice so
    that both views update: aln.Seqs[:] = x instead of aln.Seqs = x. If you
    get the two views out of sync, you will get all sorts of exceptions. No
    validation is performed on aln.Seqs and aln.Positions for performance
    reasons, so this can really get you into trouble.
    
    Alignments are immutable, though this is not enforced. If you change the
    data after the alignment is created, all sorts of bad things might happen.
    
    Class properties:
    Alphabet: should be an Alphabet object. Must provide mapping between items
    (possibly, but not necessarily, characters) in the alignment and indices
    of those characters in the resulting Alignment object.
    
    SequenceType: Constructor to use when building sequences. Default: Sequence
    
    InputHandlers: dict of {input_type:input_handler} where input_handler is
    from the InputHandlers above and input_type is a result of the method
    self._guess_input_type (should always be a string).
    
    Creating a new array will always result in a new object unless you use
    the force_same_object=True parameter.
    
    WARNING: Rebinding the Names attribute in a DenseAlignment is not
    recommended because not all methods will use the updated name order. This
    is because the original sequence and name order are used to produce data
    structures that are cached for efficiency, and are not updated if you
    change the Names attribute.
    
    WARNING: DenseAlignment strips off Info objects from sequences that have
    them, primarily for efficiency.
    """
    MolType = None             #will be set to BYTES on moltype import
    Alphabet = None            #will be set to BYTES.Alphabet on moltype import
    
    InputHandlers = {   'array':aln_from_array,
                        'model_seqs':aln_from_model_seqs,
                        'generic':aln_from_generic,
                        'fasta':aln_from_fasta,
                        'dense_aln':aln_from_dense_aln,
                        'aln':aln_from_collection,
                        'collection':aln_from_collection,
                        'dict':aln_from_dict,
                        'kv_pairs':aln_from_kv_pairs,
                        'empty':aln_from_empty,
                    }
    
    def __init__(self, *args, **kwargs):
        """Returns new DenseAlignment object. Inherits from SequenceCollection.
        """
        kwargs['suppress_named_seqs'] = True
        super(DenseAlignment, self).__init__(*args, **kwargs)
        self.ArrayPositions = transpose(\
            self.SeqData.astype(self.Alphabet.ArrayType))
        self.ArraySeqs = transpose(self.ArrayPositions)
        self.SeqData = self.ArraySeqs
        self.SeqLen = len(self.ArrayPositions)


    def _force_same_data(self, data, Names):
        """Forces array that was passed in to be used as self.ArrayPositions"""
        if isinstance(data, DenseAlignment):
            data = data._positions
        self.ArrayPositions = data
        self.Names = Names or self.DefaultNameFunction(len(data[0]))
    
    def _get_positions(self):
        """Override superclass Positions to return positions as symbols."""
        return map(self.Alphabet.fromIndices, self.ArrayPositions)
    
    Positions = property(_get_positions)
    
    def _get_named_seqs(self):
        if not hasattr(self, '_named_seqs'):
            seqs = map(self.Alphabet.toString, self.ArraySeqs)
            if self.MolType:
                seqs = map(self.MolType.Sequence, seqs)
            self._named_seqs = self._make_named_seqs(self.Names, seqs)
        return self._named_seqs
    
    NamedSeqs = property(_get_named_seqs)
    
    def keys(self):
        """Supports dict-like interface: returns names as keys."""
        return self.Names
    
    def values(self):
        """Supports dict-like interface: returns seqs as Sequence objects."""
        return [self.Alphabet.MolType.ModelSeq(i, Alphabet=self.Alphabet) \
            for i in self.ArraySeqs]
    
    def items(self):
        """Supports dict-like interface; returns (name, seq) pairs."""
        return zip(self.keys(), self.values())
    
    def __iter__(self):
        """iter(aln) iterates over positions, returning array slices.
        
        Each item in the result is be a position ('column' in standard
        terminology) within the alignment, with the sequneces in the same
        order as in the names.
        
        The result shares data with the original array, so if you change
        the result you change the Alignment.
        """
        return iter(self.Positions)
    
    def __getitem__(self, item):
        """getitem delegates to self.Positions., returning array slices.
        
        The result is a column or slice of columns, supporting full slice
        functionality (including stride). Use this to get a selection of
        positions from the alignment.
        
        Result shares data with the original array, so if you change the
        result you change the Alignment.
        """
        return self.Positions[item]
    
    def _coerce_seqs(self, seqs, is_array):
        """Controls how seqs are coerced in _names_seqs_order.
        
        Override in subclasses where this behavior should differ.
        """
        return seqs
    
    def getSubAlignment(self, seqs=None, pos=None, invert_seqs=False, \
        invert_pos=False):
        """Returns subalignment of specified sequences and positions.
        
        seqs and pos can be passed in as lists of sequence indices to keep
        or positions to keep.
        
        invert_seqs: if True (default False), gets everything _except_ the
        specified sequences.
        
        invert_pos: if True (default False), gets everything _except_ the
        specified positions.
        
        Unlike most of the other code that gets things out of an alignment,
        this method returns a new alignment that does NOT share data with the
        original alignment.
        """
        #figure out which positions to keep, and keep them
        if pos is not None:
            if invert_pos:
                pos_mask = ones(len(self.ArrayPositions))
                put(pos_mask, pos, 0)
                pos = nonzero(pos_mask)[0]
            data = take(self.ArrayPositions, pos, axis=0)
        else:
            data = self.ArrayPositions
        #figure out which sequences to keep, and keep them
        if seqs is not None:
            if invert_seqs:
                seq_mask = ones(len(self.ArraySeqs))
                put(seq_mask, seqs, 0)
                seqs = nonzero(seq_mask)[0]
            data = take(data, seqs, 1)
            names = [self.Names[i] for i in seqs]
        else:
            names = self.Names
        return self.__class__(data, map(str,names), self.Alphabet, \
            conversion_f=aln_from_array)
    
    def __str__(self):
        """Returns FASTA-format string.
        
        Should be able to handle joint alphabets, e.g. codons.
        """
        result = []
        names = map(str, self.Names)
        max_label_length = max(map(len, names)) + 1
        seq2str = self.Alphabet.fromIndices
        for l, s in zip(self.Names, self.ArraySeqs):
            result.append('>'+str(l)+'\n'+''.join(seq2str(s)))
        return '\n'.join(result) + '\n'
    
    def _get_freqs(self, index=None):
        """Gets array of freqs along index 0 (= positions) or 1 (= seqs).
        
        index: if 0, will calculate the frequency of each symbol in each
        position (=column) in the alignment. Will return 2D array where the
        first index is the position, and the second index is the index of the
        symbol in the alphabet. For example, for the TCAG DNA Alphabet,
        result[3][0] would store the count of T at position 3 (i.e. the 4th
        position in the alignment.
        
        if 1, does the same thing except that the calculation is performed for
        each sequence, so the 2D array has the sequence index as the first
        index, and the symbol index as the second index. For example, for the
        TCAG DNA Alphabet, result[3][0] would store the count of T in the
        sequence at index 3 (i.e. the 4th sequence).
        """
        if index:
            a = self.ArrayPositions
        else:
            a = self.ArraySeqs
        count_f = self.Alphabet.counts
        return array(map(count_f, a))
    
    def getPosFreqs(self):
        """Returns Profile of counts: position by character.
        
        See documentation for _get_freqs: this just wraps it and converts the
        result into a Profile object organized per-position (i.e. per column).
        """
        return Profile(self._get_freqs(1), self.Alphabet)
    
    def getSeqEntropy(self):
        """Returns array containing Shannon entropy for each seq in self.
        
        Uses the profile object from getSeqFreqs (see docstring) to calculate
        the per-symbol entropy in each sequence in the alignment, i.e. the
        uncertainty about each symbol in each sequence (or row). This can be
        used to, for instance, filter low-complexity sequences.
        """
        p = self.getSeqFreqs()
        p.normalizePositions()
        return p.rowUncertainty()
    
    def getPosEntropy(self):
        """Returns array containing Shannon entropy for each pos in self.
        
        Uses the profile object from getPosFreqs (see docstring) to calculate
        the per-symbol entropy in each position in the alignment, i.e. the
        uncertainty about each symbol at each position (or column). This can
        be used to, for instance, detect the level of conservation at each
        position in an alignment.
        """
        p = self.getPosFreqs()
        p.normalizePositions()
        return p.rowUncertainty()
    
    def IUPACConsensus(self, alphabet=None):
        """Returns string containing IUPAC consensus sequence of the alignment.
        """
        if alphabet is None:
            alphabet = self.MolType
        consensus = []
        degen = alphabet.degenerateFromSequence
        for col in self.Positions:
            consensus.append(degen(str(alphabet.ModelSeq(col, \
                Alphabet=alphabet.Alphabets.DegenGapped))))
        return coerce_to_string(consensus)
    
    def _make_gaps_ok(self, allowed_gap_frac):
        """Makes the gaps_ok function used by omitGapPositions and omitGapSeqs.
        
        Need to make the function because if it's a method of Alignment, it
        has unwanted 'self' and 'allowed_gap_frac' parameters that impede the
        use of map() in takeSeqsIf.
        
        WARNING: may not work correctly if component sequences have gaps that
        are not the Alignment gap character. This is because the gaps are
        checked at the column level (and the positions are lists), rather than
        at the row level. Working around this issue would probably cause a
        significant speed penalty.
        """
        def gaps_ok(seq):
            seq_len = len(seq)
            if hasattr(seq, 'countGaps'):
                num_gaps = seq.countGaps()
            elif hasattr(seq, 'count'):
                num_gaps = seq.count(self.Alphabet.Gap)
            else:
                num_gaps = sum(seq==self.Alphabet.GapIndex)
            return num_gaps / seq_len <= allowed_gap_frac
        
        return gaps_ok
    
    def columnFreqs(self, constructor=Freqs):
        """Returns list of Freqs with item counts for each column.
        """
        return map(constructor, self.Positions)

    def sample(self, n=None, with_replacement=False, motif_length=1, \
        randint=randint, permutation=permutation):
        """Returns random sample of positions from self, e.g. to bootstrap.

        Arguments:
            - n: the number of positions to sample from the alignment.
              Default is alignment length
            - with_replacement: boolean flag for determining if sampled
              positions
            - randint and permutation: functions for random integer in a
              specified range, and permutation, respectively.
              
            
        Notes: 
            By default (resampling all positions without replacement), generates
            a permutation of the positions of the alignment.

            Setting with_replacement to True and otherwise leaving parameters
            as defaults generates a standard bootstrap resampling of the 
            alignment.
            """
        population_size = len(self) // motif_length
        if not n:
            n = population_size
        if with_replacement:
            locations = randint(0, population_size, n)
        else:
            assert n <= population_size, (n, population_size, motif_length)
            locations = permutation(population_size)[:n]
        #check if we need to convert coords for multi-width motifs
        if motif_length > 1:
            locations = (locations*motif_length).repeat(motif_length)
            wrapped_locations =locations.reshape((n,motif_length))
            wrapped_locations += arange(motif_length)
        positions = take(self.ArrayPositions, locations, 0)
        result = self.__class__(positions.T,force_same_data=True, \
            Info=self.Info, Names=self.Names)
        return result
 
def aln_from_fasta_codons(seqs, array_type=None, Alphabet=None):
    """Codon alignment from FASTA-format string or lines.
    
    This is an InputHandler for taking a FASTA-format string of individual
    bases and converting it into an array by way of a CodonSequence object
    that groups triples of bases together and converts them into symbols on
    the codon alphabet (i.e. each group of 3 bases together is coded by a
    single symbol). This needs to override the normal aln_from_fasta
    InputHandler, which asssumes that it can convert the string into the
    array directly without this grouping step.
    """
    if isinstance(seqs, str):
        seqs = seqs.split('\n')
    return aln_from_model_seqs([CodonSequenceGap(s, Label=l) for l, s \
        in cogent.parse.fasta.MinimalFastaParser(seqs)])

    def xsample(self, n=None, with_replacement=False, motif_length=1, \
        random_series=random):
        """Returns random sample of positions from self, e.g. to bootstrap.

        Arguments:
            - n: the number of positions to sample from the alignment.
              Default is alignment length
            - with_replacement: boolean flag for determining if sampled
              positions
            - random_series: a random number generator with
              .randint(min,max) .random() methods
              
            
        Notes: 
            By default (resampling all positions without replacement), generates
            a permutation of the positions of the alignment.

            Setting with_replacement to True and otherwise leaving parameters
            as defaults generates a standard bootstrap resampling of the 
            alignment.
            """
        population_size = len(self) // motif_length
        if not n:
            n = population_size
        if with_replacement:
            locations = [random_series.randint(0, population_size)
                    for samp in xrange(n)]
        else:
            assert n <= population_size, (n, population_size, motif_length)
            locations = random_series.sample(xrange(population_size), n)
        positions = [(loc*motif_length, (loc+1)*motif_length)
                for loc in locations]
        sample = Map(positions, parent_length=len(self))
        return self.gappedByMap(sample, Info=self.Info)
 
class CodonDenseAlignment(DenseAlignment):
    """Stores alignment of gapped codons, no degenerate symbols."""
    InputHandlers = {   'array':aln_from_array,
                        'seqs':aln_from_model_seqs,
                        'generic':aln_from_generic,
                        'fasta':aln_from_fasta_codons,
                        'dense_aln':aln_from_dense_aln,
                        'aln': aln_from_collection,
                        'collection':aln_from_collection,
                        'dict':aln_from_dict,
                        'empty':aln_from_empty,
                    }

def make_gap_filter(template, gap_fraction, gap_run):
    """Returns f(seq) -> True if no gap runs and acceptable gap fraction.
    
    Calculations relative to template.
    gap_run = number of consecutive gaps allowed in either the template or seq
    gap_fraction = fraction of positions that either have a gap in the template
        but not in the seq or in the seq but not in the template
    NOTE: template and seq must both be ModelSequence objects.
    """
    template_gaps = array(template.gapVector())
    def result(seq):
        """Returns True if seq adhers to the gap threshold and gap fraction."""
        seq_gaps = array(seq.gapVector())
        #check if gap amount bad
        if sum(seq_gaps!=template_gaps)/float(len(seq)) > gap_fraction:
            return False
        #check if gap runs bad
        if '\x01'*gap_run in logical_and(seq_gaps, \
                logical_not(template_gaps)).astype(uint8).tostring():
            return False
        #check if insertion runs bad
        elif '\x01'*gap_run in logical_and(template_gaps, \
                logical_not(seq_gaps)).astype(uint8).tostring():
            return False
        return True
    
    return result

class Alignment(_Annotatable, AlignmentI, SequenceCollection):
    MolType = None  #note: this is reset to ASCII in moltype module
    def __init__(self, *args, **kwargs):
        """Returns new Alignment object: see SequenceCollection."""
        
        SequenceCollection.__init__(self, *args, **kwargs)
        
        #need to convert seqs to Aligned objects
        seqs = self.SeqData
        names = self.Names
        
        self._motif_probs = {}
        self._type = self.MolType.gettype()
        lengths = map(len, self.SeqData)
        if lengths and (max(lengths) != min(lengths)):
            raise DataError, "Not all sequences are the same length:\n" + \
                "max is %s, min is %s" % (max(lengths), min(lengths))
        aligned_seqs = []
        for s, n in zip(seqs, names):
            if isinstance(s, Aligned):
                s.Name = n  #ensure consistency
                aligned_seqs.append(s)
            else:
                aligned_seqs.append(self._seq_to_aligned(s, n))
        self.NamedSeqs = self.AlignedSeqs = dict(zip(names, aligned_seqs))
        self.SeqData = self._seqs = aligned_seqs
    
    def _coerce_seqs(self, seqs, is_array):
        if not min([isinstance(seq, _Annotatable) or isinstance(seq, Aligned) for seq in seqs]):
            seqs = map(self.MolType.Sequence, seqs)
        return seqs
    
    def _seq_to_aligned(self, seq, key):
        """Converts seq to Aligned object -- override in subclasses"""
        (map, seq) = self.MolType.Sequence(seq, key).parseOutGaps()
        return Aligned(map, seq)
    
    def getTracks(self, policy):
        # drawing code related
        # same as sequence but annotations go below sequence tracks
        return policy.tracksForAlignment(self)
    
    def getChildTracks(self, policy):
        """The only Alignment method required for cogent.draw"""
        tracks =  []
        for label in self.Names:
            seq = self.NamedSeqs[label]
            tracks += seq.getTracks(policy.copy(seqname=label))
        return tracks
    
    
    def __repr__(self):
        seqs = []
        limit = 10
        delimiter = ''
        for (count, name) in enumerate(self.Names):
            if count == 3:
                seqs.append('...')
                break
            elts = list(self.getGappedSeq(name)[:limit+1])
            if len(elts) > limit:
                elts.append('...')
            seqs.append("%s[%s]" % (name, delimiter.join(elts)))
        seqs = ', '.join(seqs)
        
        return "%s x %s %s alignment: %s" % (len(self.Names),
                self.SeqLen, self._type, seqs)
    
    def _mapped(self, slicemap):
        align = []
        for name in self.Names:
            align.append((name, self.NamedSeqs[name][slicemap]))
        return self.__class__(MolType=self.MolType, data=align)
    
    def gappedByMap(self, keep, **kwargs):
        # keep is a Map
        seqs = []
        for seq_name in self.Names:
            aligned = self.NamedSeqs[seq_name]
            seqmap = aligned.map[keep]
            seq = aligned.data.gappedByMap(seqmap)
            seqs.append((seq_name, seq))
        return self.__class__(MolType=self.MolType, data=seqs, **kwargs)
    
    def projectAnnotation(self, seq_name, annot):
        target_aligned = self.NamedSeqs[seq_name]
        if annot.parent is not self:
            raise ValueError('Annotation does not belong to this alignment')
        return annot.remappedTo(target_aligned.data, target_aligned.map)
        
    def getProjectedAnnotations(self, seq_name, *args):
        aln_annots = self.getAnnotationsMatching(*args)
        return [self.projectAnnotation(seq_name, a) for a in aln_annots]
    
    def getAnnotationsFromSequence(self, seq_name, *args):
        aligned = self.NamedSeqs[seq_name]
        return aligned.getAnnotationsMatching(self, *args)
    
    def getAnnotationsFromAnySequence(self, *args):
        result = []
        for seq_name in self.Names:
            result.extend(self.getAnnotationsFromSequence(seq_name, *args))
        return result
    
    def getBySequenceAnnotation(self, seq_name, *args):
        result = []
        for feature in self.getAnnotationsFromSequence(seq_name, *args):
            segment = self[feature.map.Start:feature.map.End]
            segment.Name = '%s "%s" %s to %s of %s' % (
                feature.type, feature.Name,
                feature.map.Start, feature.map.End, self.Name or '')
            result.append(segment)
        return result
    
    def withMaskedAnnotations(self, annot_types, mask_char=None, shadow=False):
        """returns an alignment with annot_types regions replaced by mask_char
        if shadow is False, otherwise all other regions are masked.
        
        Arguments:
            - annot_types: annotation type(s)
            - mask_char: must be a character valid for the seq MolType. The
              default value is the most ambiguous character, eg. '?' for DNA
            - shadow: whether to mask the annotated regions, or everything but
              the annotated regions"""
        masked_seqs = []
        for seq in self.Seqs:
            # we mask each sequence using these spans
            masked_seqs += [seq._masked_annotations(annot_types,mask_char,shadow)]
        new = self.__class__(data=masked_seqs, Info=self.Info, Name=self.Name)
        return new
    
    def variablePositions(self, include_gap_motif = True):
        """Return a list of variable position indexes.
        
        Arguments:
            - include_gap_motif: if False, sequences with a gap motif in a
              column are ignored."""
        seqs = [self.getGappedSeq(n) for n in self.Names]
        seq1 = seqs[0]
        positions = zip(*seqs[1:])
        result = []
        for (position, (motif1, column)) in enumerate(zip(seq1,positions)):
            for motif in column:
                if motif != motif1:
                    if include_gap_motif:
                        result.append(position)
                        break
                    elif motif != '-' and motif1 != '-':
                        result.append(position)
                        break
        
        return result
    
    def filtered(self, predicate, motif_length=1, **kwargs):
        """The alignment positions where predicate(column) is true.
        
        Arguments:
            - predicate: a callback function that takes an tuple of motifs and
              returns True/False
            - motif_length: length of the motifs the sequences should be split
              into, eg. 3 for filtering aligned codons."""
        gv = []
        kept = False
        seqs = [self.getGappedSeq(n).getInMotifSize(motif_length,
                                    **kwargs) for n in self.Names]
        
        positions = zip(*seqs)
        for (position, column) in enumerate(positions):
            keep = predicate(column)
            if kept != keep:
                gv.append(position*motif_length)
                kept = keep
        
        if kept:
            gv.append(len(positions)*motif_length)
        
        locations = [(gv[i], gv[i+1]) for i in range(0, len(gv), 2)]
        
        keep = Map(locations, parent_length=len(self))
        return self.gappedByMap(keep, Info=self.Info)
    
    def getSeq(self, seqname):
        """Return a ungapped Sequence object for the specified seqname.

        Note: always returns Sequence object, not ModelSequence.
        """
        return self.NamedSeqs[seqname].data
    
    def getGappedSeq(self, seq_name, recode_gaps=False):
        """Return a gapped Sequence object for the specified seqname.

        Note: always returns Sequence object, not ModelSequence.
        """
        return self.NamedSeqs[seq_name].getGappedSeq(recode_gaps)
    
    def iterPositions(self, pos_order=None):
        """Iterates over positions in the alignment, in order.
        
        pos_order refers to a list of indices (ints) specifying the column
        order. This lets you rearrange positions if you want to (e.g. to pull
        out individual codon positions).
        
        Note that self.iterPositions() always returns new objects, by default
        lists of elements. Use map(f, self.iterPositions) to apply the
        constructor or function f to the resulting lists (f must take a single
        list as a parameter). Note that some sequences (e.g. ViennaStructures)
        have rules that prevent arbitrary strings of their symbols from being
        valid objects.
        
        Will raise IndexError if one of the indices in order exceeds the
        sequence length. This will always happen on ragged alignments:
        assign to self.SeqLen to set all sequences to the same length.
        """
        get = self.NamedSeqs.__getitem__
        pos_order = pos_order or xrange(self.SeqLen)
        seq_order = self.Names
        aligned_objs = [get(seq) for seq in seq_order]
        seqs = map(str, aligned_objs)
        for pos in pos_order:
            yield [seq[pos] for seq in seqs]
    
    Positions = property(iterPositions)
    
    def withGapsFrom(self, template):
        """Same alignment but overwritten with the gaps from 'template'"""
        if len(self) != len(template):
            raise ValueError("Template alignment must be same length")
        gap = self.Alphabet.Gap
        tgp = template.Alphabet.Gap
        result = {}
        for name in self.Names:
            seq = self.getGappedSeq(name)
            if name not in template.Names:
                raise ValueError("Template alignment doesn't have a '%s'" 
                        % name)
            gsq = template.getGappedSeq(name)
            assert len(gsq) == len(seq)
            combo = []
            for (s,g) in zip(seq, gsq):
                if g == tgp:
                    combo.append(gap)
                else:
                    combo.append(s)
            result[name] = combo
        return Alignment(result, Alphabet=self.Alphabet.withGapMotif())
         
    def getDegappedRelativeTo(self, name):
        """Remove all columns with gaps in sequence with given name.
        
        Returns Alignment object of the same class.
        Note that the seqs in the new Alignment are always new objects.
        
        Arguments:
            - name: sequence name
        """
        if name not in self.Names:
            raise ValueError("The alignment doesn't have a sequence named '{0}'"
                .format(name))
        
        gap = self.Alphabet.Gap
        non_gap_cols = [i for i, col in enumerate(self.getGappedSeq(name)) 
                                    if col != gap]
        
        return self.takePositions(non_gap_cols)
    
    def addFromReferenceAln(self, ref_aln, before_name=None, after_name=None):
        """
        Insert sequence(s) to self based on their alignment to a reference
        sequence. Assumes the first sequence in ref_aln.Names[0] is the
        reference.
        
        By default the sequence is appended to the end of the alignment,
        this can be changed by using either before_name or after_name
        arguments.
        
        Returns Alignment object of the same class.
        
        Arguments:
            - ref_aln: reference alignment (Alignment object/series) of
              reference sequence and sequences to add.
              New sequences in ref_aln (ref_aln.Names[1:] are sequences to add.
              If series is used as ref_aln, it must have the structure
              [['ref_name', SEQ], ['name', SEQ]]
            - before_name: name of the sequence before which
              sequence is added
            - after_name: name of the sequence after which sequence is added
              If both before_name and after_name are specified seqs will be
              inserted using before_name.
        
        Example:
        Aln1:
        -AC-DEFGHI (name: seq1) 
        XXXXXX--XX (name: seq2)
        YYYY-YYYYY (name: seq3) 
        
        Aln2:
        ACDEFGHI   (name: seq1) 
        KL--MNPR   (name: seqX)
        KLACMNPR   (name: seqY)
        KL--MNPR   (name: seqZ)
        
        Out:
        -AC-DEFGHI (name: seq1) 
        XXXXXX--XX (name: seq2)
        YYYY-YYYYY (name: seq3)
        -KL---MNPR (name: seqX)
        -KL-ACMNPR (name: seqY)
        -KL---MNPR (name: seqZ)
        """
        
        if type(ref_aln) != type(self): # let the seq class try and guess
            ref_aln = self.__class__(ref_aln)
        
        ref_seq_name = ref_aln.Names[0]
        
        if ref_seq_name not in self.Names:
            raise ValueError, "The name of reference sequence ({0})"\
            "not found in the alignment \n(names in the alignment:\n{1}\n)"\
            .format(ref_seq_name, "\n".join(self.Names))
        
        if str(ref_aln.getGappedSeq(ref_seq_name)) \
            != str(self.getSeq(ref_seq_name)):
            raise ValueError, "Reference sequences are unequal."\
            "The reference sequence must not contain gaps"
        
        temp_aln = None
        for seq_name in ref_aln.Names[1:]:
            if seq_name in self.Names:
                raise ValueError, "The name of a sequence being added ({0})"\
                        "is already present".format(seq_name)
            
            seq = ref_aln.getGappedSeq(seq_name)
            new_seq = Aligned(self.NamedSeqs[ref_seq_name].map, seq)
            if not temp_aln:
                temp_aln = self.__class__({new_seq.Name: str(new_seq)})
            else:
                temp_aln = temp_aln.addSeqs(self.__class__({new_seq.Name:
                                            str(new_seq)}))
        
        aln = self.addSeqs(temp_aln, before_name, after_name)
        
        return aln
