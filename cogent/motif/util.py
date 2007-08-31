#!/usr/bin/env python
"""Utility classes for general motif and module API."""

from __future__ import division
import string
import re
from cogent.core.alignment import Alignment
from cogent.core.location import Span

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Prototype"

class Location(Span):
    """Object that stores location information for a module

    -Sequence refers to the original sequence the module came from
    -SeqId is the key of the sequence in the alignment
    -Start is the position in the sequence

    """

    def __init__(self, SeqId, Start, End=None):
        """Initializes location object"""
        self.SeqId = SeqId
        Span.__init__(self,Start, End)

    def __cmp__(self,other):
        """Overwriting __cmp__ for sorting purposes"""
        return cmp(self.SeqId, other.SeqId)

class ModuleInstanceI(object):
    """Object that stores individual module instance information.

    Contains sequence, location, Pvalue and Evalue of a module instance as well
    as some basic instance functions.
    """

    def __init__(self, Sequence, Location, Pvalue=None, Evalue=None):
        """Initializes ModuleInstance object"""
        self.Sequence = Sequence
        self.Location = Location    #Location Object
        self.Pvalue = Pvalue
        self.Evalue = Evalue

    def distance(self,other):
        """Calculates the distance between two ModuleInstances"""
        raise NotImplementedError

    def __cmp__(self,other):
        """Overwriting __cmp__ function to compare ModuleInstance objects"""
        if self is other:
            return 0
        return cmp(self.Pvalue,other.Pvalue) \
               or cmp(self.Evalue,other.Evalue) \
               or cmp(self.Location,other.Location) \
               or cmp(str(self),str(other))
    
    def __lt__(self, other):
        return cmp(self, other) == -1
    
    def __le__(self, other):
        return cmp(self, other) <= 0

    def __gt__(self, other):
        return cmp(self, other) == 1

    def __ge__(self, other):
        return cmp(self, other) >= 0

    def __eq__(self, other):
        return self.__cmp__(other) == 0
    
    def __ne__(self, other):
        return cmp(self, other) != 0
    
class ModuleInstanceStr(ModuleInstanceI, str):
    """Constructor for ModuleInstance inheriting from string."""
    def __new__(cls, data='', *args, **kwargs):
        return str.__new__(cls, data)
    def __init__(self, *args, **kwargs):
        return ModuleInstanceI.__init__(self, *args, **kwargs)

def ModuleInstance(data, Location, Pvalue=None, Evalue=None, constructor=None):
    """Creates ModuleInstance given a constructor."""
    if constructor is None:
        #maybe add code to try to figure out what to do from the data later
        constructor=ModuleInstanceStr
    return constructor(data, Location, Pvalue, Evalue)

def seqs_from_empty(obj, *args, **kwargs):
    """Allows empty initialization of Module, useful when data must be added."""
    return [], []
    
class Module(Alignment):
    """Object that stores module information.

    Module is an Alignment of ModuleInstances.  Constructed as a dict keyed by
    location with ModuleInstance sequence as the value:
        - {(SeqId, Start): ModuleInstance}
    """
    InputHandlers = Alignment.InputHandlers.copy()
    InputHandlers['empty'] = seqs_from_empty

    def __init__(self, data=None, Template=None, MolType=None,\
                 Locations=None, Pvalue=None, Evalue=None, Llr=None,\
                 ID=None,ConsensusSequence=None):
        """Initializes Module object"""
        self.Template = Template
        if MolType is not None:
            self.MolType = MolType
        self.Pvalue = Pvalue
        self.Evalue = Evalue
        self.Llr = Llr      #Log likelihood ratio
        self.ID = ID
        self.ConsensusSequence = ConsensusSequence
        if isinstance(data, dict):
            data = sorted(data.items())
        else:
            try:
                data = sorted(data)
            except TypeError:
                pass
        super(Module, self).__init__(data, MolType=MolType)

    def update(self, other):
        """Updates self with info in other, in-place. WARNING: No validation!"""
        self.Names += other.Names
        self.NamedSeqs.update(other.NamedSeqs)

    def __setitem__(self, item, val):
        """Replaces item in self.NamedSeqs. WARNING: No validation!"""
        if item not in self.NamedSeqs:
            self.Names.append(item)
        self.NamedSeqs[item] = val

    def __repr__(self):
        return str(self.NamedSeqs)


    def __str__(self):
        """Returns string representation of IUPAC consensus sequence"""
        if len(self.MolType.Alphabet) < 20:
            return str(self.IUPACConsensus(self.MolType))
        return str(''.join(self.majorityConsensus()))

    
    def distance(self,other):
        """Calculates the distance between two Modules"""
        raise NotImplementedError

    def __cmp__(self,other):
        """Overwriting __cmp__ function to compare Module objects"""
        return cmp(self.Pvalue,other.Pvalue) \
               or cmp(self.Evalue,other.Evalue)

    def __hash__(self):
        """overwriting __hash__ function to hash Module object"""
        return id(self)

    def _get_location_dict(self):
        """Returns a dict of module locations.

        Represented as a dict with SeqId as key and [indices] as values:
            {SeqId:[indices]}
        """
        location_dict = {}
        for key in self.Names:
            try:
                location_dict[key[0]].append(key[1])
            except:
                location_dict[key[0]]=[key[1]]
        return location_dict

    LocationDict = property(_get_location_dict)
    
    def _get_loose(self):
        """Returns a list of all ModuleInstances not in self.Strict.
        """
        loose_list = []
        strict = self.Strict[0].Sequence
        for instance in self.values():
            if instance.Sequence != strict:
                loose_list.append(instance)
        return loose_list
        
    Loose = property(_get_loose)
    
    def _get_strict(self):
        """Returns a list of ModuleInstances with the most common sequence.
        """
        strict_dict = {} #Dictionary to hold counts of instance strings.
        #For each ModuleInstance in self.
        for instance in self.values():
            #If instance already in strict_dict then increment and append.
            if instance.Sequence in strict_dict:
                strict_dict[instance.Sequence][0]+=1
                strict_dict[instance.Sequence][1].append(instance)
            #Else, add count and instance to dict.
            else:
                strict_dict[instance.Sequence]=[1,[instance]]
        #List with all counts and instances
        count_list = strict_dict.values()
        count_list.sort()
        count_list.reverse()
        #Set self.Template as the Strict ModuleInstance sequence.
        self.Template = count_list[0][1][0].Sequence
        #Return list of ModuleInstances with the most common sequence.
        return count_list[0][1]
    
    Strict = property(_get_strict)
    
    def basePossibilityCount(self,degenerate_dict=None):
        """Returns number of possible combinations to form a degenerate string.
        """
        if degenerate_dict is None:
            degenerate_dict = self.MolType.Degenerates
        #Get degenerate string representation of module
        degenerate_string = self.__str__()
        #Get length of first degenerate character
        combinations = len(degenerate_dict.get(degenerate_string[0],'-'))
        #Multiply number of possibilities for each degenerate character together
        for i in range(1, len(degenerate_string)):
            combinations *= len(degenerate_dict.get(degenerate_string[i],'-'))
        #Return total possible ways to make module
        return combinations

    def _coerce_seqs(self, seqs, is_array):
        """Override _coerce_seqs so we keep the orig objects."""
        return seqs

    def _seq_to_aligned(self, seq, key):
        """Override _seq_to_aligned so we keep the orig objects."""
        return seq

class ModuleFinder(object):
    """Object that constructs a dict of modules given an alignment"""

    def __call__(self, *args):
        """Call method for ModuleFinder"""
        raise NotImplementedError

class ModuleConsolidator(object):
    """Object that takes in a list of modules and returns a consolidated list.
    Modules that are very similar are considered the same module.
    """

    def __call__(self, *args):
        """Call method for ModuleConsolidator"""
        raise NotImplementedError

class Motif(object):
    """Object that stores modules that are considered the same motif
    """

    def __init__(self, Modules=None, Info=None):
        """Initializes Motif object"""
        self.Modules = []
        try:
            #only one module in motif
            self.Modules.append(Modules)
        except:
            #list of modules
            self.Modules.extend(Modules)
        self.Info = Info

class MotifFinder(object):
    """Object that takes modules and constructs motifs

        - Takes in a list of modules and constructs a list of Motifs"""

    def __call__(self, *args):
        """Call method for MotifFinder"""
        raise NotImplementedError

class MotifFormatter(object):
    """Object that takes a list of Motifs and formats them for output to browser

        - Takes in a list of motifs and generates specified output format.
    """
    COLORS = [ "#00FF00", "#FFFF00", "#00FFFF", "#FF00FF",
               "#C0C0C0", "#FAEBD7", "#8A2BE2", "#A52A2A", "#00CC00", "#FF6600",
               "#FF33CC", "#CC33CC", "#9933FF", "#FFCCCC", "#00CCCC", "#999999",
               "#CC6666", "#CCCC33", "#66CCFF" ]

    STYLES = ["", "font-weight: bold", "font-style: italic"]
    def getColorMapS0(self, module_ids):
        """ Standalone version - needed b/c of pickle problem """
        color_map = {}
        mod = len(MotifFormatter.COLORS)
        smod = len(MotifFormatter.STYLES)
        for module_id in module_ids: 
            ix = int(module_id)
            cur_color = ix % mod
            cur_style = int(round((ix / mod))) % smod
            style_str = """background-color: %s; %s; font-family: 'Courier New', Courier"""
            color_map[module_id] = style_str % (
                                MotifFormatter.COLORS[cur_color],
                                MotifFormatter.STYLES[cur_style])

        return color_map

    def getColorMap(self, motif_results):
        """ Return color mapping for motif_results """
        module_ids = []
        for motif in motif_results.Motifs:
            for module in motif.Modules:
                module_ids.append(module.ID)
        return self.getColorMapS0(sorted(module_ids))
    
    def getColorMapRgb(self, motif_results):
        """ Return color mapping for motif_results using RGB rather than hex.
        """
        module_ids = []
        for motif in motif_results.Motifs:
            for module in motif.Modules:
                module_ids.append(module.ID)
        color_map = {}
        mod = len(MotifFormatter.COLORS)
        for module_id in module_ids: 
            ix = int(module_id)
            cur_color = ix % mod
            color_map[module_id] = \
                html_color_to_rgb(MotifFormatter.COLORS[cur_color])

        return color_map


    def __call__(self, *args):
        """Call method for MotifFormatter"""
        raise NotImplementedError

def html_color_to_rgb(colorstring):
    """ convert #RRGGBB to an (R, G, B) tuple 
    
        - From Python Cookbook.
    """
    colorstring = colorstring.strip()
    if colorstring[0] == '#': colorstring = colorstring[1:]
    if len(colorstring) != 6:
        raise ValueError, "input #%s is not in #RRGGBB format" % colorstring
    r, g, b = colorstring[:2], colorstring[2:4], colorstring[4:]
    r, g, b = [int(n, 16) for n in (r, g, b)]
    return (r, g, b)
    
class MotifResults(object):
    """Object that holds a list of Modules, Motifs and a dict of Results.
    """

    def __init__(self,Modules=None, Motifs=None, Results=None, Parameters=None,
                 Alignment=None,MolType=None):
        """Initializes MotifResults object."""
        self.Modules = Modules or []
        self.Motifs = Motifs or []
        self.Results = Results or {} #Results not belonging to other categories.
        if Parameters:
            self.__dict__.update(Parameters)
        self.Alignment = Alignment
        self.MolType = MolType
