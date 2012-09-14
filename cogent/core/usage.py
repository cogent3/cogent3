#!/usr/bin/env python
"""Classes for dealing with base, codon, and amino acid usage.
"""
from __future__ import division
from cogent.maths.stats.util import Freqs, Numbers, UnsafeFreqs
from cogent.maths.stats.special import fix_rounding_error
from cogent.util.array import euclidean_distance
from cogent.util.misc import Delegator, FunctionWrapper, InverseDict
from cogent.core.genetic_code import GeneticCodes, GeneticCode as GenCodeClass
from cogent.core.info import Info as InfoClass
from cogent.core.alphabet import CharAlphabet
from string import upper

from numpy import array, concatenate, sum, mean, isfinite, sqrt

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Sandra Smit", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

RnaBases = CharAlphabet('UCAG')
DnaBases = CharAlphabet('TCAG')
AminoAcids = CharAlphabet('ACDEFGHIKLMNPQRSTVWY*') # * denotes termination
AB = CharAlphabet('ab') #used for testing
Chars = CharAlphabet(''.join(map(chr, range(256))), '-')  #used for raw chars
RnaBasesGap = CharAlphabet('UCAG-', '-')
DnaBasesGap = CharAlphabet('TCAG-', '-')
AminoAcidsGap = CharAlphabet('ACDEFGHIKLMNPQRSTVWY*-', '-')
DnaIupac = CharAlphabet('TCAGNVBHDKSWMYR')
RnaIupac = CharAlphabet('UCAGNVBHDKSWMYR')
AminoAcidsIupac = CharAlphabet('ACDEFGHIKLMNPQRSTVWY*XBZ')
DnaIupacGap = CharAlphabet('TCAG-NVBHDKSWMYR', '-')
RnaIupacGap = CharAlphabet('UCAG-NVBHDKSWMYR', '-')
AminoAcidsIupacGap = CharAlphabet('ACDEFGHIKLMNPQRSTVWY*-XBZ', '-')

RnaPairs = RnaBases**2
DnaPairs = DnaBases**2
RnaGapPairs = RnaBasesGap**2
DnaGapPairs = DnaBasesGap**2
AminoAcidPairs = AminoAcids**2
ABPairs = AB**2

#RnaBases = 'UCAG'
#DnaBases = 'TCAG'
RnaCodons = [i+j+k for i in RnaBases for j in RnaBases for k in RnaBases]
DnaCodons = [i+j+k for i in DnaBases for j in DnaBases for k in DnaBases]
#AminoAcids = 'ACDEFGHIKLMNPQRSTVWY*'
SGC = GeneticCodes[1]

RnaDinucs = [i+j for i in RnaBases for j in RnaBases]

RnaToDna = dict(zip(RnaBases, DnaBases))
DnaToRna = dict(zip(DnaBases, RnaBases))

Bases = RnaBases        #by default
Codons = RnaCodons      #by default

_equal_bases = Freqs(Bases)
_equal_codons = Freqs(Codons)
_equal_amino_acids = Freqs(AminoAcids[:-1])     #exclude Stop
for i in (_equal_bases, _equal_codons, _equal_amino_acids):
    i.normalize()

empty_rna_codons = dict.fromkeys(RnaCodons, 0.0)
empty_dna_codons = dict.fromkeys(DnaCodons, 0.0)

def seq_to_codon_dict(seq, empty_codons=empty_dna_codons):
    """Converts sequence into codon dict."""
    leftover = len(seq) % 3
    if leftover:
        seq += 'A' * (3-leftover)
    result = empty_codons.copy()
    for i in range(0, len(seq), 3):
        curr = seq[i:i+3]
        if curr in result:  #ignore others
            result[curr] += 1
    return result

def UnsafeCodonsFromString(seq, rna=False, formatted=False, **kwargs):
    if rna:
        d = empty_rna_codons
        if not formatted:
            seq = seq.upper().replace('T','U')
    else:
        d = empty_dna_codons
        if not formatted:
            seq = seq.upper().replace('U','T')
    return UnsafeCodonUsage(seq_to_codon_dict(seq,d), **kwargs)

def key_to_rna(key):
    """Sets key to uppercase RNA."""
    return key.upper().replace('T', 'U')

def key_to_dna(key):
    """Sets key to uppercase DNA."""
    return key.upper().replace('U', 'T')

class InfoFreqs(Freqs, Delegator):
    """Like Freqs, but has an Info object storing additional data.
    
    Intended for holding base or codon frequencies that come from a 
    particular sequence, so that the Info of the sequence can be preserved
    even if the sequence is deleted to free up memory.
    """
    def __init__(self, data=None, Info=None, **kwargs):
        """Intializes BaseUsage with data, either sequence or dict of freqs.
        
        Ignores additional kwargs (e.g. to support copy).

        Makes the _handler for delegator accessible with the name Info.
        """
        if Info is None:
            if hasattr(data, 'Info'):
                Info = data.Info
            else:
                Info = InfoClass()
        Delegator.__init__(self, Info)
        Freqs.__init__(self, data or [], **kwargs)

    def _get_info(self):
        """Accessor for Info."""
        return self._handler
    def _set_info(self, obj):
        """Mutator for Info."""
        self._handler = obj
    Info = property(_get_info, _set_info)


class BaseUsageI(object):
    """Provides shared interface for BaseUsage classes.

    BaseUsage stores counts of the four DNA or RNA bases.
    """
    def bases(self):
        """Supports bases/codons/positionalBases/aminoAcids interface."""
        return self

    def codons(self):
        """Predicts codon frequencies from the base frequencies."""
        result = {}
        base_copy = self.__class__(self)
        base_copy.normalize()
        for c in Codons:
            curr = 1
            for i in c:
                curr *= base_copy[i]
            result[c] = curr
        return CodonUsage(result, self.Info)

    def positionalBases(self):
        """Returns PositionalBaseUsage with copy of self at each position."""
        return PositionalBaseUsage(self.__class__(self), self.__class__(self),
            self.__class__(self), self.Info)

    def aminoAcids(self, genetic_code=SGC):
        """Predicts amino acid frequencies from the base frequencies."""
        return self.codons().aminoAcids(genetic_code)
    
    def distance(self,other):
        """Calculates the distance between two BaseUsages.

        Distance is measured in three directions, CG-content, CU-content, and
        GU-content.
        """
        return euclidean_distance(array(self.toCartesian()),\
            array(other.toCartesian()))

    def content(self, string):
        """Gets the sum of bases specified in string.
        
        For example, self.content('GC') gets the GC content.
        """
        return sum([self.get(i, 0) for i in string], 0)

    def toCartesian(self):
        """Returns tuple of x, y, z coordinates from BaseUsage.

        x=u+c, y=u+g, z=u+a

        Doesn't alter original object.
        """
        return self['UC'], self['UG'], self['UA']

    def fromCartesian(cls, *coords):
        """Returns new BaseUsage with A,C,G,U coordinates from UC,UG,UA.

        From UC,UG,UA to A,C,G(,U).

        This will only work when the original coordinates come from a simplex,
        where U+C+A+G=1
        """
        result = cls()
        x,y,z = coords
        u = fix_rounding_error((1-x-y-z)/-2)
        a, c, g = z-u, x-u, y-u
        result['A'] = a
        result['C'] = c
        result['G'] = g
        result['U'] = u
        return result
    fromCartesian = classmethod(fromCartesian)

class BaseUsage(BaseUsageI, InfoFreqs):
    """Stores frequencies of the four bases, mapped to RNA.

    This class is convenient but inefficient, since it automatically maps any
    lookups to the uppercase RNA alphabet internally. Use UnsafeBaseUsage for 
    speed when necessary.
    """
    
    Mask = FunctionWrapper(key_to_rna)
    RequiredKeys = dict.fromkeys(Bases)
    
    def __getitem__(self, key):
        """Normalizes key and treats T=U."""
        key = self.Mask(key)
        if len(key) == 2:  #pair of bases, e.g. GC for GC content
            dup = BaseUsage(self)
            dup.normalize()
            return sum([dup.get(i,0) for i in key], 0)
        else:
            return super(BaseUsage, self).__getitem__(key)

class UnsafeBaseUsage(BaseUsageI, UnsafeFreqs):
    """Stores frequencies of the four bases. Does not do any validation.

    This class avoids overriding most built-ins, so is much faster than
    BaseFreqs (although it is often less convenient).
    """
    RequiredKeys = dict.fromkeys(Bases)
    Info = {}   # for interface compatibility with InfoFreqs-based class

class CodonUsageI(object):
    """Stores codon usage for a gene or species.
    
    Note that CodonUsage objects get their own reference to _default_code
    during creation, so changing CodonUsage._default_code will not change the
    GeneticCode of any CodonUsage object that has already been created.
    """
    _default_code = SGC

    BlockAbbreviations = \
        {'UU':'F/L', 'CU':'Leu', 'AU':'I/M', 'GU':'Val', \
        'UC':'Ser', 'CC':'Pro', 'AC':'Thr', 'GC':'Ala',\
        'UA':'Y/*','CA':'H/Q','AA':'N/K','GA':'D/E',\
        'UG':'C*W', 'CG':'Arg', 'AG':'S/R', 'GG':'Gly'}
    
    BlockNames = \
        {'UU':'Phe/Leu', 'CU':'Leucine', 'AU':'Ile/Met', 'GU':'Valine', \
        'UC':'Serine', 'CC':'Proline', 'AC':'Threonine', 'GC':'Alanine',\
        'UA':'Tyr/Ter','CA':'His/Gln','AA':'Asn/Lys','GA':'Asp/Glu',\
        'UG':'Cys/Ter/Trp', 'CG':'Arginine', 'AG':'Ser/Arg', 'GG':'Glycine'}

    Blocks = [i+j for i in 'UCAG' for j in 'UCAG'] #UCAG order
    SingleAABlocks = ['GC','CG','GG','CU','CC','UC','AC','GU'] #alpha by aa
    SplitBlocks = ['UU', 'CA','AU','AA','AG','GA'] #UCAG order]

    AbbreviationsToBlocks = InverseDict(BlockAbbreviations)
    NamesToBlocks = InverseDict(BlockNames)

    BaseUsageClass = None   #Must overrride in subclasses

    def bases(self, purge_unwanted=False):
        """Returns overall base counts."""
        result = {}
        if purge_unwanted:
            data = self._purged_data()
        else:
            data = self
        for codon, freq in data.items():
            for base in codon:
                if base in result:
                    result[base] += freq
                else:
                    result[base] = freq
        return self.BaseUsageClass(result, self.Info)

    def codons(self):
        """Supports codons/aminoAcids/bases/positionalBases interface."""
        return self

    def rscu(self):
        """Normalizes self in-place to RSCU, relative synonymous codon usage.

        RSCU divides the frequency of each codon to the sum of the freqs for
        that codon's amino acid.
        """
        gc = self.GeneticCode
        syn = gc.Synonyms
        aa_sums = {}
        for key, codons in syn.items():
            aa_sums[key] = sum([self[c] for c in codons], 0)
        for codon in self:
            try:
                curr = self[codon]
                res = curr/aa_sums[gc[codon]]
            except (KeyError, ZeroDivisionError, FloatingPointError):
                pass
            else:
                if isfinite(res):
                    self[codon] = res
        return self
        

    def _purged_data(self):
        """Copy of self's freqs after removing bad/stop codons and singlets."""
        good_codons = self.RequiredKeys
        code = self.GeneticCode
        data = dict(self)   #need copy, since we're deleting required keys
        #delete any bad codons
        for codon in self:
            if codon not in good_codons:
                del data[codon]
        #delete any stop codons in the current code
        for codon in code['*']:
            try:
                c = codon.upper().replace('T','U')
                del data[c]
            except KeyError:
                pass    #don't care if it's not there
        #delete any single-item blocks in the current code (i.e. leaving
        #only doublets and quartets).
        for group in code.Blocks:
            if len(group) == 1:
                try:
                    c = group[0].upper().replace('T','U')
                    del data[c]
                except KeyError:
                    pass    #don't care if already deleted
        return data

    def positionalBases(self, purge_unwanted=False):
        """Calculates positional base usage.
        
        purge_unwanted controls whether or not to purge 1-codon groups, stop
        codons, and any codons containing degnerate bases before calculating 
        the base usage (e.g. to get Noboru Sueoka's P3 measurement): default 
        is False. Deletion of unwanted codons happens on a copy, not the 
        original data.
        """
        first = {}
        second = {}
        third = {}

        if purge_unwanted:  #make a copy of the data and delete things from it
            data = self._purged_data()
        else:
            data = self

        for codon, freq in data.items():
            try:
                p1, p2, p3 = codon
            except ValueError:
                continue #skip any incomplete codons
            
            if p1 in first:
                first[p1] += freq
            else:
                first[p1] = freq
            if p2 in second:
                second[p2] += freq
            else:
                second[p2] = freq
            if p3 in third:
                third[p3] += freq
            else:
                third[p3] = freq
        return PositionalBaseUsage(self.BaseUsageClass(first), \
            self.BaseUsageClass(second), self.BaseUsageClass(third), self.Info)

    def positionalGC(self, purge_unwanted=True):
        """Returns GC, P1, P2 P3. Use purge_unwanted=False to get raw counts."""
        p = self.positionalBases(purge_unwanted)
        p.normalize()
        result = [i['G'] + i['C'] for i in p]
        average = sum(result, 0)/3
        return [average] + result
        
    def fingerprint(self, which_blocks='quartets', include_mean=True,\
        normalize=True):
        """Returns fingerprint data for fingerprint plots.
        
        which_blocks: whether to include only the usual 4-codon quartets (the
                      default), the split blocks only, or all blocks.
        include_mean: whether to include the mean (True).
        normalize:    whether to normalize so that the quartets sum to 1 (True)
        """ 
        if which_blocks == 'split':
            blocks = self.SplitBlocks
        elif which_blocks == 'quartets':
            blocks = self.SingleAABlocks
        elif which_blocks == 'all':
            blocks = self.Blocks
        else:
            raise "Got invalid option %s for which_blocks:\n" % which_blocks+\
                    "  (valid options: 'split', 'quartets', 'all')."
        result = []
        for b in blocks:    #iterates over doublet string
            U, C, A, G = [self[b+i] for i in 'UCAG']
            all = U+C+A+G
            if G+C:
                g_ratio = G/(G+C)
            else:
                g_ratio = 0.5

            if A+U:
                a_ratio = A/(A+U)
            else:
                a_ratio=0.5
                
            result.append([g_ratio, a_ratio, all])
        result = array(result)
        
        if normalize:   #make the shown bubbles sum to 1
            sum_ = sum(result[:,-1])
            if sum_:
                result[:,-1] /= sum_
        
        if include_mean: #calculate mean from all codons
            third = self.positionalBases().Third
            U, C, A, G = [third[i] for i in 'UCAG']
            if G+C:
                g_ratio = G/(G+C)
            else:
                g_ratio = 0.5
            if A+U:
                a_ratio = A/(A+U)
            else:
                a_ratio=0.5
            result = concatenate((result, array([[g_ratio,a_ratio,1]])))
        
        return result

    def pr2bias(self, block):
        """Calculates PR2 bias for a specified block, e.g. 'AC' or 'UU'.

        Returns tuple of:
        (G/G+C, A/A+T, G/G+A, C/C+T, G/G+T, C/C+A)

        If any pair sums to zero, will raise ZeroDivisionError.

        block: codon block, e.g. 'AC', 'UU', etc. Any of the 16 doublets.
        """
        U, C, A, G = [self[block+i] for i in 'UCAG']
        return G/(G+C), A/(A+U), G/(G+A), C/(C+U), G/(G+U), C/(C+A)
            
    def aminoAcids(self, genetic_code=None):
        """Calculates amino acid usage, optionally using a specified code."""
        if genetic_code is None:
            curr_code = self.GeneticCode
        elif isinstance(genetic_code, GenCodeClass):
            curr_code = genetic_code
        else:
            curr_code = GeneticCodes[genetic_code]
        aa = {}
        for codon, freq in self.items():
            curr_aa = curr_code[codon]
            if curr_aa in aa:
                aa[curr_aa] += freq
            else:
                aa[curr_aa] = freq
        return AminoAcidUsage(aa, self.Info)
    
class CodonUsage(CodonUsageI, InfoFreqs):
    """Stores frequencies of the 64 codons, mapped to RNA.

    This class is convenient but inefficient, since it automatically maps any
    lookups to the uppercase RNA alphabet internally. Use UnsafeBaseUsage for 
    speed when necessary.
    """
    
    Mask = FunctionWrapper(key_to_rna)
    RequiredKeys = RnaCodons
    BaseUsageClass = BaseUsage
   
    def __init__(self, data=None, Info=None, GeneticCode=None, \
        Mask=None, ValueMask=None, Constraint=None):
        """Initializes new CodonUsage with Info and frequency data.
        
        Note: Mask, ValueMask and Constraint are ignored, but must be present
        to support copy() because of the ConstrainedContainer interface.
        """
        #check if we have a sequence: if so, take it 3 bases at a time
        #this will (properly) fail on lists of tuples or anything else where 
        #the items don't behave like strings.
        try:
            codons = [''.join(data[i:i+3]) for i in xrange(0, len(data), 3)]
        except:
            codons = data
        super(CodonUsage, self).__init__(codons, Info)
        
        if GeneticCode:
            if isinstance(GeneticCode, GenCodeClass):
                curr_code = GeneticCode
            else:
                curr_code = GeneticCodes[GeneticCode]
        else:
            curr_code = self._default_code
        self.__dict__['GeneticCode'] = curr_code
 
    def __getitem__(self, key):
        """Normalizes key and treats T=U."""
        key = self.Mask(key)
        if len(key) == 2:  #pair of bases, e.g. GC for GC content
            dup = BaseUsage(self)
            dup.normalize()
            return sum([dup.get(i,0) for i in key], 0)
        else:
            return super(CodonUsage, self).__getitem__(key)

class UnsafeCodonUsage(CodonUsageI, UnsafeFreqs):
    """Stores frequencies of the four bases. Must access as RNA.

    This class avoids overriding most built-ins, so is much faster than
    CodonFreqs (although it is often less convenient).
    """
    RequiredKeys = RnaCodons
    Info = {}   # for interface compatibility with InfoFreqs-based class
    Gene=None       #for CUTG compatibility
    Species=None    # for CUTG compaitibility
    BaseUsageClass = UnsafeBaseUsage

    def __init__(self, data=None, Info=None, GeneticCode=None, \
        Mask=None, ValueMask=None, Constraint=None):
        """Initializes new CodonUsage with Info and frequency data.
        
        Note: Mask, ValueMask and Constraint are ignored, but must be present
        to support copy() because of the ConstrainedContainer interface.
        """
        #check if we have a sequence: if so, take it 3 bases at a time
        #this will (properly) fail on lists of tuples or anything else where 
        #the items don't behave like strings.
        try:
            codons = [''.join(data[i:i+3]) for i in xrange(0, len(data), 3)]
        except:
            codons = data or {}
        UnsafeFreqs.__init__(self, codons)
        #set required keys
        for k in self.RequiredKeys:
            if k not in self:
                self[k] = 0.0
        #flatten Info onto self directly for lookups
        if Info:
            self.__dict__.update(Info)
        self.Info = Info or {}
        
        if GeneticCode:
            if isinstance(GeneticCode, GenCodeClass):
                curr_code = GeneticCode
            else:
                curr_code = GeneticCodes[GeneticCode]
        else:
            curr_code = self._default_code
        self.GeneticCode = curr_code

    
class PositionalBaseUsage(Delegator):
    """Stores a BaseUsage for each of the three codon positions."""
    
    def __init__(self, First=None, Second=None, Third=None, Info=None):
        """Returns new PositionalBaseUsage with values for the 3 positions."""
        Delegator.__init__(self, Info)
        self.__dict__['First'] = First or BaseUsage()
        self.__dict__['Second'] = Second or BaseUsage()
        self.__dict__['Third'] = Third or BaseUsage()

    def _get_info(self):
        """Accessor for Info."""
        return self._handler
    def _set_info(self, obj):
        """Mutator for Info."""
        self._handler = obj
    Info = property(_get_info, _set_info)

    def __getitem__(self, index):
        """Supports lookups by index."""
        if index == 0 or index == -3:
            return self.First
        elif index == 1 or index == -2:
            return self.Second
        elif index == 2 or index == -1:
            return self.Third
        else:
            raise IndexError, "PositionalBaseUsage only has 3 positions."

    def normalize(self):
        """Normalizes each of the component base usages."""
        self.First.normalize()
        self.Second.normalize()
        self.Third.normalize()

    def bases(self):
        """Returns distribution of the four bases, summed over positions."""
        sum = BaseUsage(Info=self.Info)
        for i in self:
            sum +=i
        return sum

    def codons(self):
        """Returns codon distribution, calculated from positional freqs."""
        result = {}
        first_copy, second_copy, third_copy = map(Freqs, self)
        first_copy.normalize()
        second_copy.normalize()
        third_copy.normalize()
        for c in Codons:
            result[c] = first_copy[c[0]] * second_copy[c[1]] * third_copy[c[2]]
        return CodonUsage(result, self.Info)

    def positionalBases(self):
        """Supports bases/codons/positionalBases/aminoAcids interface."""
        return self

    def aminoAcids(self, genetic_code=None):
        """Returns amino acid distribution."""
        return self.codons().aminoAcids(genetic_code)
        
class AminoAcidUsage(InfoFreqs):
    """Stores counts ofthe 20 canonical amino acids."""
    Mask = FunctionWrapper(upper)
    RequiredKeys = dict.fromkeys(AminoAcids)

    def bases(self, genetic_code=SGC, codon_usage=_equal_codons):
        """Predicts most likely set of base frequencies.
        
        Optionally uses a genetic code (default: standard genetic code) and 
        codon usage (default: unbiased codon usage).
        """
        result = self.codons(genetic_code, codon_usage).bases()
        result.normalize()
        return result


    def codons(self, genetic_code=SGC, codon_usage=_equal_codons):
        """Predicts most likely set of codon frequencies.

        Optionally uses genetic_code (to figure out which codons belong
        with each amino acid), and codon_usage (to get most likely codons for 
        each amino acid). Defaults are the standard genetic code and unbiased 
        codon frequencies.
        """
        result = {}
        normalized = Freqs(self)
        normalized.normalize()
        for aa, aa_freq in normalized.items():
            curr_codons = [c.upper().replace('T','U') for c in genetic_code[aa]]
            if not curr_codons:
                continue    #code might be missing some amino acids?
            curr_codon_freqs = Numbers([codon_usage[c] for c in curr_codons])
            curr_codon_freqs.normalize()
            for codon, c_freq in zip(curr_codons, curr_codon_freqs):
                result[codon] = c_freq * aa_freq
        return CodonUsage(result, self.info, genetic_code)

    def positionalBases(self, genetic_code=SGC, codon_usage=_equal_codons):
        """Predicts most likely set of positional base frequencies.

        Optionally uses a genetic code (default: standard genetic code) and
        codon usage (default: unbiased codon usage).
        """
        return self.codons(genetic_code, codon_usage).positionalBases()

    def aminoAcids(self):
        """Supports bases/positionalBases/aminoAcids/codons interface."""
        return self

class DinucI(object):
    """Provides shared interface for DinucUsage classes.

    DinucUsage stores counts of the 16 doublets.
    """

    def distance(self, other):
        """Calculates distance between two DinucUsage objects."""
        result = 0
        for k in self.RequiredKeys:
            result += (self[k]-other[k])**2
        return sqrt(result)

class DinucUsage(DinucI, InfoFreqs):
    """Stores frequencies of the 16 dinucleotides, mapped to RNA.

    This class is convenient but inefficient, since it automatically maps any
    lookups to the uppercase RNA alphabet internally. Use UnsafeBaseUsage for 
    speed when necessary.
    """
    Mask = FunctionWrapper(key_to_rna)
    RequiredKeys = RnaDinucs

    def __init__(self, data=None, Info=None, Overlapping=True, \
        GeneticCode=None, Mask=None, ValueMask=None, Constraint=None):
        """Initializes new CodonUsage with Info and frequency data.
        
        Note: Mask, ValueMask and Constraint are ignored, but must be present
        to support copy() because of the ConstrainedContainer interface.
        """
        #check if we have a sequence: if so, take it 3 bases at a time
        #this will (properly) fail on lists of tuples or anything else where 
        #the items don't behave like strings.
        if Mask is not None:
            self.Mask = Mask
        if ValueMask is not None:
            self.ValueMask = ValueMask
        try:
            data = self.Mask(data)
            if Overlapping == '3-1':
                range_ = range(2, len(data)-1, 3)
            elif Overlapping:
                range_ = range(0, len(data)-1)
            else:
                range_ = range(0, len(data)-1, 2)
            dinucs = [''.join(data[i:i+2]) for i in range_]
        except:
            dinucs = data
        super(DinucUsage, self).__init__(dinucs, Info)
        
    def __getitem__(self, key):
        """Normalizes key and treats T=U."""
        key = self.Mask(key)
        return super(DinucUsage, self).__getitem__(key)


#some useful constants...

EqualBases = BaseUsage()
EqualBases = BaseUsage(_equal_bases)
EqualPositionalBases = PositionalBaseUsage(BaseUsage(_equal_bases),
    BaseUsage(_equal_bases), BaseUsage(_equal_bases))
EqualAminoAcids = AminoAcidUsage(_equal_amino_acids)   #excl. Stop
EqualCodons = CodonUsage(_equal_codons)
