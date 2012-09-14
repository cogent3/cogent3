#!/usr/bin/env python
"""
moltype.py

MolType provides services for resolving ambiguities, or providing the
correct ambiguity for recoding. It also maintains the mappings between
different kinds of alphabets, sequences and alignments.

One issue with MolTypes is that they need to know about Sequence, Alphabet,
and other objects, but, at the same time, those objects need to know about
the MolType. It is thus essential that the connection between these other
types and the MolType can be made after the objects are created.
"""

__author__ = "Peter Maxwell, Gavin Huttley and Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight", \
        "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

from cogent.core.alphabet import CharAlphabet, Enumeration, Alphabet, \
    AlphabetError, _make_complement_array
from cogent.util.misc import FunctionWrapper, add_lowercase, iterable, if_
from cogent.util.transform import allchars, keep_chars
from cogent.data.molecular_weight import DnaMW, RnaMW, ProteinMW
from cogent.core.sequence import Sequence as DefaultSequence, RnaSequence, \
    DnaSequence, ProteinSequence, ABSequence, NucleicAcidSequence, \
    ByteSequence, ModelSequence, ModelNucleicAcidSequence, \
    ModelDnaSequence, ModelRnaSequence, ModelDnaCodonSequence, \
    ModelRnaCodonSequence, ModelProteinSequence, ProteinWithStopSequence,\
    ModelProteinWithStopSequence
from cogent.core.genetic_code import DEFAULT as DEFAULT_GENETIC_CODE, \
    GeneticCodes
from cogent.core.alignment import Alignment, DenseAlignment, \
    SequenceCollection
from random import choice

import re
import string

import numpy
from numpy import array, sum, transpose, remainder, zeros, arange, newaxis, \
    ravel, asarray, fromstring, take, uint8, uint16, uint32

Float = numpy.core.numerictypes.sctype2char(float)
Int = numpy.core.numerictypes.sctype2char(int)

from string import maketrans, translate

IUPAC_gap = '-'

IUPAC_missing = '?'

IUPAC_DNA_chars = ['T','C','A','G']
IUPAC_DNA_ambiguities = {
    'N': ('A','C','T','G'),
    'R': ('A','G'),
    'Y': ('C','T'),
    'W': ('A','T'),
    'S': ('C','G'),
    'K': ('T','G'),
    'M': ('C','A'),
    'B': ('C','T','G'),
    'D': ('A','T','G'),
    'H': ('A','C','T'),
    'V': ('A','C','G')
        }
IUPAC_DNA_ambiguities_complements = {
    'A':'T','C':'G','G':'C','T':'A', '-':'-',
    'M':'K', 'K':'M',
    'N':'N',
    'R':'Y', 'Y':'R',
    'W':'W',
    'S':'S',
    'X':'X', # not technically an IUPAC ambiguity, but used by repeatmasker
    'V':'B', 'B':'V',
    'H':'D', 'D':'H'
    }

IUPAC_DNA_complements = {
    'A':'T','C':'G','G':'C','T':'A', '-':'-',
    }

IUPAC_DNA_complements = {
    'A':'T','C':'G','G':'C','T':'A', '-':'-',
    }


IUPAC_RNA_chars = ['U','C','A','G']     #note change in standard order from DNA
IUPAC_RNA_ambiguities = {
    'N': ('A','C','U','G'),
    'R': ('A','G'),
    'Y': ('C','U'),
    'W': ('A','U'),
    'S': ('C','G'),
    'K': ('U','G'),
    'M': ('C','A'),
    'B': ('C','U','G'),
    'D': ('A','U','G'),
    'H': ('A','C','U'),
    'V': ('A','C','G')
        }

IUPAC_RNA_ambiguities_complements = {
    'A':'U','C':'G','G':'C','U':'A', '-':'-',
    'M':'K', 'K':'M',
    'N':'N',
    'R':'Y', 'Y':'R',
    'W':'W',
    'S':'S',
    'X':'X', # not technically an IUPAC ambiguity, but used by repeatmasker
    'V':'B', 'B':'V',
    'H':'D', 'D':'H'
    }

IUPAC_RNA_complements = {
    'A':'U','C':'G','G':'C','U':'A', '-':'-',
    }


#Standard RNA pairing: GU pairs count as 'weak' pairs
RnaStandardPairs = {
    ('A','U'): True,    #True vs False for 'always' vs 'sometimes' pairing
    ('C','G'): True,
    ('G','C'): True,
    ('U','A'): True,
    ('G','U'): False,
    ('U','G'): False,
}

#Watson-Crick RNA pairing only: GU pairs don't count as pairs
RnaWCPairs = {
    ('A','U'): True,
    ('C','G'): True,
    ('G','C'): True,
    ('U','A'): True,
}

#RNA pairing with GU counted as standard pairs
RnaGUPairs = {
    ('A','U'): True,
    ('C','G'): True,
    ('G','C'): True,
    ('U','A'): True,
    ('G','U'): True,
    ('U','G'): True,
}

#RNA pairing with GU, AA, GA, CA and UU mismatches allowed as weak pairs
RnaExtendedPairs = {
    ('A','U'): True,
    ('C','G'): True,
    ('G','C'): True,
    ('U','A'): True,
    ('G','U'): False,
    ('U','G'): False,
    ('A','A'): False,
    ('G','A'): False,
    ('A','G'): False,
    ('C','A'): False,
    ('A','C'): False,
    ('U','U'): False,
}
#Standard DNA pairing: only Watson-Crick pairs count as pairs
DnaStandardPairs = {
    ('A','T'): True,
    ('C','G'): True,
    ('G','C'): True,
    ('T','A'): True,
}


# protein letters & ambiguity codes
IUPAC_PROTEIN_chars = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
    'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
    'W', 'Y']

PROTEIN_WITH_STOP_chars = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
    'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
    'W', 'Y', '*']

IUPAC_PROTEIN_ambiguities = {
    'B': ['N', 'D'],
    'X': IUPAC_PROTEIN_chars,
    'Z': ['Q', 'E'],
    }

PROTEIN_WITH_STOP_ambiguities = {
    'B': ['N', 'D'],
    'X': PROTEIN_WITH_STOP_chars,
    'Z': ['Q', 'E'],
    }

class FoundMatch(Exception):
    """Raised when a match is found in a deep loop to skip many levels"""
    pass

def make_matches(monomers=None, gaps=None, degenerates=None):
    """Makes a dict of symbol pairs (i,j) -> strictness.
    
    Strictness is True if i and j always match and False if they sometimes
    match (e.g. A always matches A, but W sometimes matches R).
    """
    result = {}
    #allow defaults to be left blank without problems
    monomers = monomers or {}
    gaps = gaps or {}
    degenerates = degenerates or {}
    #all monomers always match themselves and no other monomers
    for i in monomers:
        result[(i,i)] = True
    #all gaps always match all other gaps
    for i in gaps:
        for j  in gaps:
            result[(i,j)] = True
    #monomers sometimes match degenerates that contain them
    for i in monomers:
        for j in degenerates:
            if i in degenerates[j]:
                result[(i,j)] = False
                result[(j,i)] = False
    #degenerates sometimes match degenerates that contain at least one of
    #the same monomers
    for i in degenerates:
        for j in degenerates:
            try:
                for i_symbol in degenerates[i]:
                    if i_symbol in degenerates[j]:
                        result[(i,j)] = False
                        raise FoundMatch
            except FoundMatch:
                pass    #flow control: break out of doubly nested loop
    return result

def make_pairs(pairs=None, monomers=None, gaps=None, degenerates=None):
    """Makes a dict of symbol pairs (i,j) -> strictness.
    
    Expands pairs into all possible pairs using degen symbols.
    Strictness is True if i and j always pair, and False if they 'weakly' pair
    (e.g. GU pairs or if it is possible that they pair).
    
    If you want to make GU pairs count as 'always matching', pass in pairs
    that have (G,U) and (U, G) mapped to True rather than False.
    """
    result = {}
    #allow defaults to be left blank without problems
    pairs = pairs or {}
    monomers = monomers or {}
    gaps = gaps or {}
    degenerates = degenerates or {}
    #add in the original pairs: should be complete monomer pairs
    result.update(pairs)
    #all gaps 'weakly' pair with each other
    for i in gaps:
        for j in gaps:
            result[(i,j)] = False
    #monomers sometimes pair with degenerates if the monomer's complement
    #is in the degenerate symbol
    for i in monomers:
        for j in degenerates:
            found = False
            try:
                for curr_j in degenerates[j]:
                    #check if (i,curr_j) and/or (curr_j,i) is a valid pair:
                    #not mutually required if pairs are not all commutative!
                    if (i, curr_j) in pairs:
                        result[(i,j)] = False
                        found = True
                    if (curr_j, i) in pairs:
                        result[(j,i)] = False
                        found = True
                    if found:
                        raise FoundMatch
            except FoundMatch:
                pass    #flow control: break out of nested loop
    #degenerates sometimes pair with each other if the first degenerate
    #contains the complement of one of the bases in the second degenerate
    for i in degenerates:
        for j in degenerates:
            try:
                for curr_i in degenerates[i]:
                    for curr_j in degenerates[j]:
                        if (curr_i, curr_j) in pairs:
                            result[(i,j)] = False
                            raise FoundMatch
            except FoundMatch:
                pass    #just using for flow control
    #don't forget the return value!
    return result

#RnaPairingRules is a dict of {name:(base_pairs,degen_pairs)} where base_pairs
#is a dict with the non-degenerate pairing rules and degen_pairs is a dict with
#both the degenerate and non-degenerate pairing rules.
#NOTE: uses make_pairs to augment the initial dict after construction.
RnaPairingRules = {
    'Standard': RnaStandardPairs,
    'WC':       RnaWCPairs,
    'GU':       RnaGUPairs,
    'Extended': RnaExtendedPairs,
}

for k, v in RnaPairingRules.items():
    RnaPairingRules[k] = (v, make_pairs(v))

class CoreObjectGroup(object):
    """Container relating gapped, ungapped, degen, and non-degen objects."""
    _types = ['Base', 'Degen', 'Gap', 'DegenGap']
    
    def __init__(self, Base, Degen=None, Gapped=None, DegenGapped=None):
        """Returns new CoreObjectGroup. Only Base is required"""
        self.Base = Base
        self.Degen = Degen
        self.Gapped = Gapped
        self.DegenGapped = DegenGapped
        self._items = [Base, Degen, Gapped, DegenGapped]
        self._set_relationships()
    
    def _set_relationships(self):
        """Sets relationships between the different "flavors"."""
        self.Base.Gapped = self.Gapped
        self.Base.Ungapped = self.Base
        self.Base.Degen = self.Degen
        self.Base.NonDegen = self.Base

        statements = [
        "self.Degen.Gapped = self.DegenGapped",
        "self.Degen.Ungapped = self.Degen",
        "self.Degen.Degen = self.Degen",
        "self.Degen.NonDegen = self.Base",
        
        "self.Gapped.Gapped = self.Gapped",
        "self.Gapped.Ungapped = self.Base",
        "self.Gapped.Degen = self.DegenGapped",
        "self.Gapped.NonDegen = self.Gapped",
        
        "self.DegenGapped.Gapped = self.DegenGapped",
        "self.DegenGapped.Ungapped = self.Degen",
        "self.DegenGapped.Degen = self.DegenGapped",
        "self.DegenGapped.NonDegen = self.Gapped",
        ]
        for s in statements:
            try:
                exec(s)
            except AttributeError:
                pass
    
    def __getitem__(self, i):
        """Allows container to be indexed into, by type of object (e.g. Gap)."""
        return self.__dict__[i]
    
    def whichType(self, a):
        """Returns the type of an alphabet in self, or None if not present."""
        return self._types[self._items.find(a)]
    

class AlphabetGroup(CoreObjectGroup):
    """Container relating gapped, ungapped, degen, and non-degen alphabets."""
    
    def __init__(self, chars, degens, gap=IUPAC_gap, missing=IUPAC_missing, \
        MolType=None, constructor=None):
        """Returns new AlphabetGroup."""
        if constructor is None:
            if max(map(len, chars)) == 1:
                constructor = CharAlphabet
                chars = ''.join(chars)
                degens = ''.join(degens)
            else:
                constructor = Alphabet  #assume multi-char
        self.Base = constructor(chars, MolType=MolType)
        self.Degen = constructor(chars+degens, MolType=MolType)
        self.Gapped = constructor(chars+gap, gap, MolType=MolType)
        self.DegenGapped = constructor(chars+gap+degens+missing, gap, \
            MolType=MolType)
        self._items = [self.Base, self.Degen, self.Gapped, self.DegenGapped]
        self._set_relationships()
        #set complements if MolType was specified
        if MolType is not None:
            comps = MolType.Complements
            for i in self._items:
                i._complement_array = _make_complement_array(i, comps)

class MolType(object):
    """MolType: Handles operations that depend on the sequence type (e.g. DNA).
    
    The MolType knows how to connect alphabets, sequences, alignments, and so
    forth, and how to disambiguate ambiguous symbols and perform base
    pairing (where appropriate).
    
    WARNING: Objects passed to a MolType become associated with that MolType,
    i.e. if you pass ProteinSequence to a new MolType you make up, all
    ProteinSequences will now be associated with the new MolType. This may
    not be what you expect. Use preserve_existing_moltypes=True if you
    don't want to reset the moltype.
    """
    def __init__(self, motifset, Gap=IUPAC_gap, Missing=IUPAC_missing,\
            Gaps=None,
            Sequence=None, Ambiguities=None,
            label=None, Complements=None, Pairs=None, MWCalculator=None, \
            add_lower=False, preserve_existing_moltypes=False, \
            make_alphabet_group=False, ModelSeq=None):
        """Returns a new MolType object. Note that the parameters are in flux.
        
        Currently:
            motifset: Alphabet or sequence of items in the default
                alphabet. Does not include degenerates.
            
            Gap: default gap symbol
            
            Missing: symbol for missing data
            
            Gaps: any other symbols that should be treated as gaps (doesn't have
                  to include Gap or Missing; they will be silently added)
            
            Sequence: Class for constructing sequences.
            
            Ambiguities: dict of char:tuple, doesn't include gaps (these are
                hard-coded as - and ?, and added later.
            
            label: text label, don't know what this is used for. Unnecessary?
            
            Complements: dict of symbol:symbol showing how the non-degenerate
                single characters complement each other. Used for constructing
                on the fly the complement table, incl. support for mustPair and
                canPair.
            
            Pairs: dict in which keys are pairs of symbols that can pair
                with each other, values are True (must pair) or False (might
                pair). Currently, the meaning of GU pairs as 'weak' is conflated
                with the meaning of degenerate symbol pairs (which might pair
                with each other but don't necessarily, depending on how the
                symbol is resolved). This should be refactored.
            
            MWCalculator: f(seq) -> molecular weight.
            
            add_lower: if True (default: False) adds the lowercase versions of
                everything into the alphabet. Slated for deletion.
            
            preserve_existing_moltypes: if True (default: False), does not
            set the MolType of the things added in **kwargs to self.
            
            make_alphabet_group: if True, makes an AlphabetGroup relating
            the various alphabets to one another.
            
            ModelSeq: sequence type for modeling
        
        Note on "Degenerates" versus "Ambiguities": self.Degenerates contains
        _only_ mappings for degenerate symbols, whereas self.Ambiguities
        contains mappings for both degenerate and non-degenerate symbols.
        Sometimes you want one, sometimes the other, so both are provided.
        """
        self.Gap = Gap
        self.Missing = Missing
        self.Gaps = frozenset([Gap, Missing])
        if Gaps:
            self.Gaps = self.Gaps.union(frozenset(Gaps))
        self.label = label
        #set the sequence constructor
        if Sequence is None:
            Sequence = ''.join     #safe default string constructor
        elif not preserve_existing_moltypes:
            Sequence.MolType = self
        self.Sequence = Sequence
        
        #set the ambiguities
        ambigs = {self.Missing:tuple(motifset)+(self.Gap,),self.Gap:(self.Gap,)}
        if Ambiguities:
            ambigs.update(Ambiguities)
        for c in motifset:
            ambigs[c] = (c,)
        self.Ambiguities = ambigs
        
        #set Complements -- must set before we make the alphabet group
        self.Complements = Complements or {}
        
        if make_alphabet_group: #note: must use _original_ ambiguities here
            self.Alphabets = AlphabetGroup(motifset, Ambiguities, \
                MolType=self)
            self.Alphabet = self.Alphabets.Base
        else:
            if isinstance(motifset, Enumeration):
                self.Alphabet = motifset
            elif max(len(motif) for motif in motifset) == 1:
                self.Alphabet = CharAlphabet(motifset, MolType=self)
            else:
                self.Alphabet = Alphabet(motifset, MolType=self)
        #set the other properties
        self.Degenerates = Ambiguities and Ambiguities.copy() or {}
        self.Degenerates[self.Missing] = ''.join(motifset)+self.Gap
        self.Matches = make_matches(motifset, self.Gaps, self.Degenerates)
        self.Pairs = Pairs and Pairs.copy() or {}
        self.Pairs.update(make_pairs(Pairs, motifset, self.Gaps, \
            self.Degenerates))
        self.MWCalculator = MWCalculator
        #add lowercase characters, if we're doing that
        if add_lower:
            self._add_lowercase()
        #cache various other data that make the calculations faster
        self._make_all()
        self._make_comp_table()
        # a gap can be a true gap char or a degenerate character, typically '?'
        # we therefore want to ensure consistent treatment across the definition
        # of characters as either gap or degenerate
        self.GapString = ''.join(self.Gaps)
        strict_gap = "".join(set(self.GapString) - set(self.Degenerates))
        self.stripDegenerate = FunctionWrapper(
            keep_chars(strict_gap+''.join(self.Alphabet)))
        self.stripBad = FunctionWrapper(keep_chars(''.join(self.All)))
        to_keep = set(self.Alphabet) ^ set(self.Degenerates) - set(self.Gaps)
        self.stripBadAndGaps = FunctionWrapper(keep_chars(''.join(to_keep)))
        
        #make inverse degenerates from degenerates
        #ensure that lowercase versions also exist if appropriate
        inv_degens = {}
        for key, val in self.Degenerates.items():
            inv_degens[frozenset(val)] = key.upper()
            if add_lower:
                inv_degens[frozenset(''.join(val).lower())] = key.lower()
        for m in self.Alphabet:
            inv_degens[frozenset(m)] = m
            if add_lower:
                inv_degens[frozenset(''.join(m).lower())] = m.lower()
        for m in self.Gaps:
            inv_degens[frozenset(m)] = m
        self.InverseDegenerates = inv_degens
        
        #set array type for modeling alphabets
        try:
            self.ArrayType = self.Alphabet.ArrayType
        except AttributeError:
            self.ArrayType = None
        
        #set modeling sequence
        self.ModelSeq = ModelSeq
    
    def __repr__(self):
        """String representation of MolType.
        
        WARNING: This doesn't allow you to reconstruct the object in its present
        incarnation.
        """
        return 'MolType(%s)' % (self.Alphabet,)
    
    def gettype(self):
        """Returns type, e.g. 'dna', 'rna', 'protein'. Delete?"""
        return self.label
    
    def makeSequence(self, Seq, Name=None, **kwargs):
        """Returns sequence of correct type. Replace with just self.Sequence?"""
        return self.Sequence(Seq, Name, **kwargs)
    
    def verifySequence(self, seq, gaps_allowed=True, wildcards_allowed=True):
        """Checks whether sequence is valid on the default alphabet.
        
        Has special-case handling for gaps and wild-cards. This mechanism is
        probably useful to have in parallel with the validation routines that
        check specifically whether the sequence has gaps, degenerate symbols,
        etc., or that explicitly take an alphabet as input.
        """
        alpha = frozenset(self.Ambiguities)
        if gaps_allowed:
            alpha = alpha.union(self.Gaps)
        if wildcards_allowed:
            alpha = alpha.union(self.Missing)
        try:
            nonalpha = re.compile('[^%s]' % re.escape(''.join(alpha)))
            badchar = nonalpha.search(seq)
            if badchar:
                motif = badchar.group()
                raise AlphabetError(motif)
        except TypeError:   #not alphabetic sequence: try slow method
            for motif in seq:
                if motif not in alpha:
                    raise AlphabetError(motif)
    
    def isAmbiguity(self, querymotif):
        """Return True if querymotif is an amibiguity character in alphabet.
        
        Arguments:
            - querymotif: the motif being queried."""
        
        return len(self.Ambiguities[querymotif]) > 1
    
    def _whatAmbiguity(self, motifs):
        """The code that represents all of 'motifs', and minimal others.
        
        Does this duplicate DegenerateFromSequence directly?
        """
        most_specific = len(self.Alphabet) + 1
        result = self.Missing
        for (code, motifs2) in self.Ambiguities.items():
            for c in motifs:
                if c not in motifs2:
                    break
            else:
                if len(motifs2) < most_specific:
                    most_specific = len(motifs2)
                    result = code
        return result
    
    def whatAmbiguity(self, motifs):
        """The code that represents all of 'motifs', and minimal others.
        
        Does this duplicate DegenerateFromSequence directly?
        """
        if not hasattr(self, '_reverse_ambiguities'):
            self._reverse_ambiguities = {}
        motifs = frozenset(motifs)
        if motifs not in self._reverse_ambiguities:
            self._reverse_ambiguities[motifs] = self._whatAmbiguity(motifs)
        return self._reverse_ambiguities[motifs]
    
    def _add_lowercase(self):
        """Adds lowercase versions of keys and vals to each internal dict."""
        for name in ['Alphabet', 'Degenerates', 'Gaps', 'Complements', 'Pairs',
            'Matches']:
            curr = getattr(self, name)
            #temp hack to get around re-ordering
            if isinstance(curr, Alphabet):
                curr = tuple(curr)
            new = add_lowercase(curr)
            setattr(self, name, new)
    
    def _make_all(self):
        """Sets self.All, which contains all the symbols self knows about.
        
        Note that the value of items in self.All will be the string containing
        the possibly degenerate set of symbols that the items expand to.
        """
        all = {}
        for i in self.Alphabet:
            curr = str(i)
            all[i] = i
        for key, val in self.Degenerates.items():
            all[key] = val
        for i in self.Gaps:
            all[i] = i
        self.All = all
    
    def _make_comp_table(self):
        """Sets self.ComplementTable, which maps items onto their complements.
        
        Note: self.ComplementTable is only set if self.Complements exists.
        """
        if self.Complements:
            self.ComplementTable = maketrans(''.join(self.Complements.keys()),
                                             ''.join(self.Complements.values()))
    
    def complement(self, item):
        """Returns complement of item, using data from self.Complements.
        
        Always tries to return same type as item: if item looks like a dict,
        will return list of keys.
        """
        if not self.Complements:
            raise TypeError, \
            "Tried to complement sequence using alphabet without complements."
        try:
            return item.translate(self.ComplementTable)
        except (AttributeError, TypeError):
            item = iterable(item)
            get = self.Complements.get
            return item.__class__([get(i, i) for i in item])
    
    def rc(self, item):
        """Returns reverse complement of item w/ data from self.Complements.
        
        Always returns same type as input.
        """
        comp = list(self.complement(item))
        comp.reverse()
        if isinstance(item, str):
            return item.__class__(''.join(comp))
        else:
            return item.__class__(comp)
    
    def __contains__(self, item):
        """A MolType contains every character it knows about."""
        return item in self.All
    
    def __iter__(self):
        """A MolType iterates only over the characters in its Alphabet.."""
        return iter(self.Alphabet)
    
    def isGap(self, char):
        """Returns True if char is a gap."""
        return char in self.Gaps
    
    def isGapped(self, sequence):
        """Returns True if sequence contains gaps."""
        return self.firstGap(sequence) is not None
    
    def isDegenerate(self, sequence):
        """Returns True if sequence contains degenerate characters."""
        return self.firstDegenerate(sequence) is not None
    
    def isValid(self, sequence):
        """Returns True if sequence contains no items that are not in self."""
        try:
            return self.firstInvalid(sequence) is None
        except:
            return False
    
    def isStrict(self, sequence):
        """Returns True if sequence contains only items in self.Alphabet."""
        try:
            return (len(sequence)==0) or (self.firstNonStrict(sequence) is None)
        except:
            return False
    
    def isValidOnAlphabet(self, sequence, alphabet=None):
        """Returns True if sequence contains only items in alphabet.
        
        Alphabet can actually be anything that implements __contains__.
        Defaults to self.Alphabet if not supplied.
        """
        if alphabet is None:
            alphabet = self.Alphabet
        return first_index_in_set(sequence, alphabet) is not None
    
    def firstNotInAlphabet(self, sequence, alphabet=None):
        """Returns index of first item not in alphabet, or None.
        
        Defaults to self.Alphabet if alphabet not supplied.
        """
        if alphabet is None:
            alphabet = self.Alphabet
        return first_index_in_set(sequence, alphabet)
    
    def firstGap(self, sequence):
        """Returns the index of the first gap in the sequence, or None."""
        gap = self.Gaps
        for i, s in enumerate(sequence):
            if s in gap:
                return i
        return None
    
    def firstDegenerate(self, sequence):
        """Returns the index of first degenerate symbol in sequence, or None."""
        degen = self.Degenerates
        for i, s in enumerate(sequence):
            if s in degen:
                return i
        return None
    
    def firstInvalid(self, sequence):
        """Returns the index of first invalid symbol in sequence, or None."""
        all = self.All
        for i, s in enumerate(sequence):
            if not s in all:
                return i
        return None
    
    def firstNonStrict(self, sequence):
        """Returns the index of first non-strict symbol in sequence, or None."""
        monomers = self.Alphabet
        for i, s in enumerate(sequence):
            if not s in monomers:
                return i
        return None
    
    def disambiguate(self, sequence, method='strip'):
        """Returns a non-degenerate sequence from a degenerate one.
        
        method can be 'strip' (deletes any characters not in monomers or gaps)
        or 'random'(assigns the possibilities at random, using equal
        frequencies).
        """
        if method == 'strip':
            try:
                return sequence.__class__(self.stripDegenerate(sequence))
            except:
                ambi = self.Degenerates
                def not_ambiguous(x):
                    return not x in ambi
                return sequence.__class__(filter(not_ambiguous, sequence))
        
        elif method == 'random':
            degen = self.Degenerates
            result = []
            for i in sequence:
                if i in degen:
                    result.append(choice(degen[i]))
                else:
                    result.append(i)
            if isinstance(sequence, str):
                return sequence.__class__(''.join(result))
            else:
                return sequence.__class__(result)
        else:
            raise NotImplementedError, "Got unknown method %s" % method
    
    def degap(self, sequence):
        """Deletes all gap characters from sequence."""
        try:
            return sequence.__class__(sequence.translate( \
            allchars, self.GapString))
        except AttributeError:
            gap = self.Gaps
            def not_gap(x):
                return not x in gap
            return sequence.__class__(filter(not_gap, sequence))
    
    def gapList(self, sequence):
        """Returns list of indices of all gaps in the sequence, or []."""
        gaps = self.Gaps
        return [i for i, s in enumerate(sequence) if s in gaps]
    
    def gapVector(self, sequence):
        """Returns list of bool indicating gap or non-gap in sequence."""
        return map(self.isGap, sequence)
    
    def gapMaps(self, sequence):
        """Returns tuple containing dicts mapping between gapped and ungapped.
        
        First element is a dict such that d[ungapped_coord] = gapped_coord.
        Second element is a dict such that d[gapped_coord] = ungapped_coord.
        
        Note that the dicts will be invalid if the sequence changes after the
        dicts are made.
        
        The gaps themselves are not in the dictionary, so use d.get() or test
        'if pos in d' to avoid KeyErrors if looking up all elements in a gapped
        sequence.
        """
        ungapped = {}
        gapped = {}
        num_gaps = 0
        for i, is_gap in enumerate(self.gapVector(sequence)):
            if is_gap:
                num_gaps += 1
            else:
                ungapped[i] = i - num_gaps
                gapped[i - num_gaps] = i
        return gapped, ungapped
    
    def countGaps(self, sequence):
        """Counts the gaps in the specified sequence."""
        gaps = self.Gaps
        gap_count = 0
        for s in sequence:
            if s in gaps:
                gap_count += 1
        return gap_count
    
    def countDegenerate(self, sequence):
        """Counts the degenerate bases in the specified sequence."""
        degen = self.Degenerates
        degen_count = 0
        for s in sequence:
            if s in degen:
                degen_count += 1
        return degen_count
    
    def possibilities(self, sequence):
        """Counts number of possible sequences matching the sequence.
        
        Uses self.Degenerates to decide how many possibilites there are at
        each position in the sequence.
        """
        degen = self.Degenerates
        count = 1
        for s in sequence:
            if s in degen:
                count *= len(degen[s])
        return count
    
    def MW(self, sequence, method='random', delta=None):
        """Returns the molecular weight of the sequence.
        
        If the sequence is ambiguous, uses method (random or strip) to
        disambiguate the sequence.
        
        if delta is present, uses it instead of the standard weight adjustment.
        """
        if not sequence:
            return 0
        try:
            return self.MWCalculator(sequence, delta)
        except KeyError:    #assume sequence was ambiguous
            return self.MWCalculator(self.disambiguate(sequence, method), delta)
    
    def canMatch(self, first, second):
        """Returns True if every pos in 1st could match same pos in 2nd.
        
        Truncates at length of shorter sequence.
        Gaps are only allowed to match other gaps.
        """
        m = self.Matches
        for pair in zip(first, second):
            if pair not in m:
                return False
        return True
    
    def canMismatch(self, first, second):
        """Returns True if any position in 1st could cause a mismatch with 2nd.
        
        Truncates at length of shorter sequence.
        Gaps are always counted as matches.
        """
        m = self.Matches
        if not first or not second:
            return False
        
        for pair in zip(first, second):
            if not m.get(pair, None):
                return True
        return False
    
    def mustMatch(self, first, second):
        """Returns True if all positions in 1st must match positions in second."""
        return not self.canMismatch(first, second)
    
    def canPair(self, first, second):
        """Returns True if first and second could pair.
        
        Pairing occurs in reverse order, i.e. last position of second with
        first position of first, etc.
        
        Truncates at length of shorter sequence.
        Gaps are only allowed to pair with other gaps, and are counted as 'weak'
        (same category as GU and degenerate pairs).
        
        NOTE: second must be able to be reverse
        """
        p = self.Pairs
        sec = list(second)
        sec.reverse()
        for pair in zip(first, sec):
            if pair not in p:
                return False
        return True
    
    def canMispair(self, first, second):
        """Returns True if any position in 1st could mispair with 2nd.
        
        Pairing occurs in reverse order, i.e. last position of second with
        first position of first, etc.
        
        Truncates at length of shorter sequence.
        Gaps are always counted as possible mispairs, as are weak pairs like GU.
        """
        p = self.Pairs
        if not first or not second:
            return False
        
        sec = list(second)
        sec.reverse()
        for pair in zip(first, sec):
            if not p.get(pair, None):
                return True
        return False
    
    def mustPair(self, first, second):
        """Returns True if all positions in 1st must pair with second.
        
        Pairing occurs in reverse order, i.e. last position of second with
        first position of first, etc.
        """
        return not self.canMispair(first, second)
    
    def degenerateFromSequence(self, sequence):
        """Returns least degenerate symbol corresponding to chars in sequence.
        
        First tries to look up in self.InverseDegenerates. Then disambiguates
        and tries to look up in self.InverseDegenerates. Then tries converting
        the case (tries uppercase before lowercase). Raises TypeError if
        conversion fails.
        """
        symbols = frozenset(sequence)
        #check if symbols are already known
        inv_degens = self.InverseDegenerates
        result = inv_degens.get(symbols, None)
        if result:
            return result
        #then, try converting the symbols
        degens = self.All
        converted = set()
        for sym in symbols:
            for char in degens[sym]:
                converted.add(char)
        symbols = frozenset(converted)
        result = inv_degens.get(symbols, None)
        if result:
            return result
        #then, try converting case
        symbols = frozenset([s.upper() for s in symbols])
        result = inv_degens.get(symbols, None)
        if result:
            return result
        symbols = frozenset([s.lower() for s in symbols])
        result = inv_degens.get(symbols, None)
        if result:
            return result
        #finally, try to find the minimal subset containing the symbols
        symbols = frozenset([s.upper() for s in symbols])
        lengths = {}
        for i in inv_degens:
            if symbols.issubset(i):
                lengths[len(i)] = i
        if lengths:  #found at least some matches
            sorted = lengths.keys()
            sorted.sort()
            return inv_degens[lengths[sorted[0]]]
        
        #if we got here, nothing worked
        raise TypeError, "Cannot find degenerate char for symbols: %s" \
                % symbols

ASCII = MolType(
    # A default type for text read from a file etc. when we don't
    # want to prematurely assume DNA or Protein.
    Sequence = DefaultSequence,
    motifset = string.letters,
    Ambiguities = {},
    label = 'text',
    ModelSeq = ModelSequence,
    )

DNA = MolType(
    Sequence = DnaSequence,
    motifset = IUPAC_DNA_chars,
    Ambiguities = IUPAC_DNA_ambiguities,
    label = "dna",
    MWCalculator = DnaMW,
    Complements = IUPAC_DNA_ambiguities_complements,
    Pairs = DnaStandardPairs,
    make_alphabet_group=True,
    ModelSeq = ModelDnaSequence,
    )

RNA = MolType(
    Sequence = RnaSequence,
    motifset = IUPAC_RNA_chars,
    Ambiguities = IUPAC_RNA_ambiguities,
    label = "rna",
    MWCalculator = RnaMW,
    Complements = IUPAC_RNA_ambiguities_complements,
    Pairs = RnaStandardPairs,
    make_alphabet_group=True,
    ModelSeq = ModelRnaSequence,
    )

PROTEIN = MolType(
    Sequence = ProteinSequence,
    motifset = IUPAC_PROTEIN_chars,
    Ambiguities = IUPAC_PROTEIN_ambiguities,
    MWCalculator = ProteinMW,
    make_alphabet_group=True,
    ModelSeq = ModelProteinSequence,
    label = "protein")

PROTEIN_WITH_STOP = MolType(
    Sequence = ProteinWithStopSequence,
    motifset = PROTEIN_WITH_STOP_chars,
    Ambiguities = PROTEIN_WITH_STOP_ambiguities,
    MWCalculator = ProteinMW,
    make_alphabet_group=True,
    ModelSeq = ModelProteinWithStopSequence,
    label = "protein_with_stop")

BYTES = MolType(
    # A default type for arbitrary chars read from a file etc. when we don't
    # want to prematurely assume _anything_ about the data.
    Sequence = ByteSequence,
    motifset = map(chr, range(256)),
    Ambiguities = {},
    ModelSeq = ModelSequence,
    label = 'bytes')

#following is a two-state MolType useful for testing
AB = MolType(
    Sequence = ABSequence,
    motifset = 'ab',
    Ambiguities={},
    ModelSeq = ModelSequence,
    label='ab')

class _CodonAlphabet(Alphabet):
    """Codon alphabets are DNA TupleAlphabets with a genetic code attribute and some codon-specific methods"""
    
    def _with(self, motifs):
        a = Alphabet._with(self, motifs)
        a.__class__ = type(self)
        a._gc = self._gc
        return a
    
    def isCodingCodon(self, codon):
        return not self._gc.isStop(codon)
    
    def isStopCodon(self, codon):
        return self._gc.isStop(codon)
    
    def getGeneticCode(self):
        return self._gc
    

def CodonAlphabet(gc=DEFAULT_GENETIC_CODE, include_stop_codons=False):
    if isinstance(gc, (int, basestring)):
        gc = GeneticCodes[gc]
    if include_stop_codons:
        motifset = list(gc.Codons)
    else:
        motifset = list(gc.SenseCodons)
    motifset = [codon.upper().replace('U', 'T') for codon in motifset]
    a = _CodonAlphabet(motifset, MolType=DNA)
    a._gc = gc
    return a

def _method_codon_alphabet(ignore, *args, **kwargs):
    """If CodonAlphabet is set as a property, it gets self as extra 1st arg."""
    return CodonAlphabet(*args, **kwargs)

STANDARD_CODON = CodonAlphabet()

#Modify NucleicAcidSequence to avoid circular import
NucleicAcidSequence.CodonAlphabet = _method_codon_alphabet
NucleicAcidSequence.PROTEIN = PROTEIN
ModelRnaSequence.MolType = RNA
ModelRnaSequence.Alphabet = RNA.Alphabets.DegenGapped

ModelDnaSequence.MolType = DNA
ModelDnaSequence.Alphabet = DNA.Alphabets.DegenGapped

ModelProteinSequence.MolType = PROTEIN
ModelProteinSequence.Alphabet = PROTEIN.Alphabets.DegenGapped

ModelProteinWithStopSequence.MolType = PROTEIN_WITH_STOP
ModelProteinWithStopSequence.Alphabet = PROTEIN_WITH_STOP.Alphabets.DegenGapped

ModelSequence.Alphabet = BYTES.Alphabet

DenseAlignment.Alphabet = BYTES.Alphabet
DenseAlignment.MolType = BYTES

ModelDnaCodonSequence.Alphabet = DNA.Alphabets.Base.Triples
ModelRnaCodonSequence.Alphabet = RNA.Alphabets.Base.Triples

#Modify Alignment to avoid circular import
Alignment.MolType = ASCII
SequenceCollection.MolType = BYTES
