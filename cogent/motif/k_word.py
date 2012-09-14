#!/usr/bin/env python
"""MotifFinder that searches for over-represented k-words in an alignment."""

from __future__ import division
from cogent.motif.util import Motif, Module, ModuleFinder, ModuleConsolidator, \
    MotifFinder, Location, ModuleInstance, MotifResults
from cogent.core.bitvector import PackedBases
from cogent.maths.stats.test import combinations, multiple_comparisons
from cogent.maths.stats.distribution import poisson_high
from numpy import array, fromstring

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Gavin Huttley", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

class KWordModuleFinder(ModuleFinder):
    """ModuleFinder that finds all k-words in an alignment.
    """
    
    def __init__(self,Alignment,MolType):
        """Initializing KWordModuleFinder.
        """
        self.ModuleDict = {}
        self.ModuleOrder = []
        self.Alignment = Alignment
        self.MolType = MolType
        self._make_char_array_aln()
    
    def _make_char_array_aln(self):
        """Turns self.Alignment into a character array.
        """
        for k,v in self.Alignment.items():
            self.Alignment[k]=array(v,'c')
        
            
    def __call__(self,word_length):
        """Builds a dict of all Modules and a list of their order.
        
            - module_dict is {module pattern:Module object}
            - module_order is a list in descending order of their count.
        """
        #Dictionary keying k-word to Module
        self.ModuleDict = {}
        #For each sequence in the alignment
        for key,seq in self.Alignment.items():
            #For each position in seq till end - word_length
            for i in range(0,len(seq)-word_length+1):
                #Get the current k-word
                word = seq[i:i+word_length].tostring()
                #Create a location object
                location = Location(key,i,i+word_length)
                #Create a ModuleInstance
                curr_instance = ModuleInstance(word,location)
                #Check to see if pattern is already in dict
                if word in self.ModuleDict:
                    #Add instance to Module
                    self.ModuleDict[word][(key,i)]=curr_instance
                #Not in dict
                else:
                    #Create a new module and add to dict
                    self.ModuleDict[word]=Module({(key,i):curr_instance},\
                        MolType=self.MolType)
        #Get list of counts
        module_counts = \
            [(len(mod.Names),word) for word,mod in self.ModuleDict.items()]
        #Sort and put in descending order
        module_counts.sort()
        module_counts.reverse()
        #Get list of only the words in descending order
        self.ModuleOrder = [word for i,word in module_counts]

    
class KWordModuleConsolidatorNucleotide(ModuleConsolidator):
    """Consolidates module instances obtained from an alignment into modules.
    
        - Must be initialized with an instance of KWordModuleFinder
    """
    
    def __init__(self,KFinder): 
        """Initializing KWordModuleConsolidatorNucleotide.
        """
        self.Modules=[]
        self.KFinder=KFinder
    
    def __call__(self,mismatches):
        """Consolidates ModuleInstances in KFinder.ModuleInstances into Modules.
        
            - mismatches accounts for the difference between two ModuleInstances
                in terms of purine vs pyrimidine or strong vs weak bonding.
        """
        #New list to hold consolidated modules
        consolidated_list = []
        #Bind module_dict and module_order locally
        module_dict = self.KFinder.ModuleDict
        module_order = self.KFinder.ModuleOrder
        #Check that dict is not empty
        if module_dict:
            #Iterate through modules in order
            for pattern in module_order:
                #Create a Bitvector representation of pattern
                pat_vec = PackedBases(pattern)
                added=False
                #Iterate through consolidated_list
                for curr_module,curr_vec in consolidated_list:
                    #If pat_vec and curr_vec are in the allowed mismatch cutoff
                    if sum(pat_vec ^ curr_vec) <= mismatches:
                        #Add module information to curr_module
                        curr_module.update(module_dict[pattern])
                        added=True
                        break
                if not added:
                    #Add new module to consolidated list
                    consolidated_list.append([module_dict[pattern],pat_vec])
        #self.Modules should have only the list of modules
        self.Modules = [i[0] for i in consolidated_list]
        for mod in self.Modules:
            mod.Template = str(mod)

class KWordModuleConsolidatorProtein(ModuleConsolidator):
    """Consolidates module instances obtained from an alignment into modules.
    
        - Must be initialized with an instance of KWordModuleFinder
    """
    
    def __init__(self,KFinder): 
        """Initializing KWordModuleConsolidatorProtein.
        """
        self.Modules=[]
        self.KFinder=KFinder
    
    def __call__(self,mismatches):
        """Consolidates ModuleInstances in KFinder.ModuleInstances into Modules.
        
            - mismatches accounts for the difference between two ModuleInstances
                in terms of purine vs pyrimidine or strong vs weak bonding.
        """
        #New list to hold consolidated modules
        consolidated_list = []
        #Bind module_dict and module_order locally
        module_dict = self.KFinder.ModuleDict
        module_order = self.KFinder.ModuleOrder
        #Check that dict is not empty
        if module_dict:
            #Iterate through modules in order
            for pattern in module_order:
                #Create a Bitvector representation of pattern
                pat_array = fromstring(pattern,'c') 
                added=False
                #Iterate through consolidated_list
                for curr_module,curr_array in consolidated_list:
                    #If pat_vec and curr_vec are in the allowed mismatch cutoff
                    if sum(pat_array != curr_array) <= mismatches:
                        #Add module information to curr_module
                        curr_module.update(module_dict[pattern])
                        added=True
                        break
                if not added:
                    #Add new module to consolidated list
                    consolidated_list.append([module_dict[pattern],pat_array])
        #self.Modules should have only the list of modules
        self.Modules = [i[0] for i in consolidated_list]

class KWordModuleFilterProtein(ModuleConsolidator):
    """Filters list of Modules by number of Modules and seqs where Module is in.
    
        - This is a strict ModuleConsolidator.  Does not allow mismatches, but
        rather consolidates based on a minimum allowed modules and minimum
        number of sequences that a module must be present in.
    """
    
    def __init__(self,KFinder,Alignment): 
        """Initializing KWordModuleFilterProtein.
        """
        self.Modules=[]
        self.KFinder=KFinder
        self.Alignment=Alignment
    
    def __call__(self,min_modules,min_seqs):
        """Adds modules to self.Modules from self.KFinder.ModuleDict as follows:
        
            - modules that have at least min_modules and are present in at least
            min_seqs will be added to self.Modules.
        """
        #First filter by min_modules
        for k,v in self.KFinder.ModuleDict.items():
            #If number of modules is greater than the minimum
            if len(v) >= min_modules:
                #Check to see that module is present in at least min_seqs
                curr_seqs = {}
                for curr in v.keys():
                    curr_seqs[curr[0]]=True
                if len(curr_seqs) >= min_seqs:
                    self.Modules.append(self.fixModuleSequence(v))
    
    def fixModuleSequence(self,module):
        """Remaps original (non-reduced) sequence string for each ModuleInstance
        """
        module_len=len(str(module))
        module.Template=str(module)
        for k,v in module.items():
            seq_id, module_start = k
            module_end = module_start+module_len
            loc = Location(seq_id,module_start,module_end)
            curr_str = \
                self.Alignment[seq_id][module_start:module_end]
            curr_instance = ModuleInstance(curr_str,loc)
            module[k]=curr_instance
        return module


class KWordMotifFinder(MotifFinder):
    """Constructs a list of Motifs from a list of Modules.
    
        - Must be initialized with a list of Modules and an Alignment.
    """
    
    def __init__(self,Modules,Alignment,Mismatches,BaseFrequency):
        """Initializing KWordMotifFinder.
        """
        self.Modules=Modules        #List of modules from Alignment
        self.Alignment=Alignment    #Alignment used to find modules
        self.Mismatches=Mismatches  #Number of mismatches used for consolidation
        self.BaseFrequency=BaseFrequency
        

    def __call__(self,threshold,alphabet_len=None,use_loose=True):
        """Returns a MotifResults object containing a list of significant Motifs
        
            - When a Module is found that is less than the threshold, a new 
                Motif object is created and added to the list of significant
                motifs.
        """
        #Create new MotifResults object
        k_word_results = MotifResults()
        k_word_results.Modules=[]
        k_word_results.Motifs=[]
        #Add alignment to MotifResults
        k_word_results.Alignment = self.Alignment
        #Give each module a unique ID
        module_id = 0
        #For each module        
        for i,module in enumerate(self.Modules):
            
            if use_loose:
                p_curr = self.getPValueLoose(module,alphabet_len)
            else:
                p_curr = self.getPValueStrict(module,alphabet_len)
            #If the P value is less than or equal to the threshold
            if p_curr <= threshold:
                #Add the P value to the current module
                module.Pvalue = p_curr
                module.ID = module_id
                module_id+=1
                #Append the module to the k_word_results Modules list
                k_word_results.Modules.append(module)
                #Create a Motif from the module and append to Motifs list
                k_word_results.Motifs.append(Motif(module))
        return k_word_results
        
    
    def getPValueStrict(self,module, alphabet_len=None):
        """Returns the Pvalue of a module.
        
            - pass alphabet_len if alphabet other than module.MolType was used
            to find modules (i.e. using a protein reduced alphabet)
        """
        #Length of the module
        module_len = len(str(module))
        
        #if moltype length has not been passed, get it from module.MolType
        #
        if not alphabet_len:
            alphabet_len = len(module.MolType.Alphabet)
        
        #Total length of the alignment
        aln_len=sum(map(len,self.Alignment.values()))
        #Number of sequences in the alignment
        num_seqs = len(self.Alignment)
        #get module_p
        module_p = 1
        for char in module.Template:
            module_p *= self.BaseFrequency.get(char,1)
        
        #Mean for passing to poisson_high NOT using loose correction
        strict_mean = \
            (aln_len + num_seqs*(1-module_len))*(module_p)
            
        #Strict P value from poisson_high
        strict_p_value = poisson_high(len(module)-1,strict_mean)
        #Correct P value for multiple comparisons
        strict_p_corrected = \
            multiple_comparisons(strict_p_value,alphabet_len**module_len)

        return strict_p_corrected

    
    def getPValueLoose(self,module, alphabet_len=None):
        """Returns the Pvalue of a module.
        
            - pass alphabet_len if alphabet other than module.MolType was used
            to find modules (i.e. using a protein reduced alphabet)
        """
        #Length of the module
        module_len = len(module.Template)
        
        #if moltype length has not been passed, get it from module.MolType
        #
        if not alphabet_len:
            alphabet_len = len(module.MolType.Alphabet)
        
        #Total length of the alignment
        aln_len=sum(map(len,self.Alignment.values()))
        #Number of sequences in the alignment
        num_seqs = len(self.Alignment)
        

        #get module_p
        module_p = 1
        for char in module.Template:
            module_p *= self._degen_p(char,module.MolType)
        
        #Mean for passing to poisson_high using loose correction
        
        loose_mean = \
            (aln_len + num_seqs*(1-module_len))*(module_p)
        #Loose P value from poisson_high
        loose_p_value = poisson_high(len(module)-1,loose_mean)
        #Correct P value for multiple comparisons
        loose_p_corrected = \
            multiple_comparisons(loose_p_value,alphabet_len**module_len)
       
        return loose_p_corrected
    
    def _degen_p(self,char,alphabet):
        """Returns the sum of the probabilities of seeing each degenerate char.
        """
        all = alphabet.Degenerates.get(char,char)
        return sum([self.BaseFrequency.get(x,1) for x in all])
            
        
