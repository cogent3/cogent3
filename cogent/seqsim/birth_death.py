#!/usr/bin/env python
"""
Code for birth-death (Yule) processes for simulating phylogenetic trees.

Also contains a double birth-death model for simulating horizontal gene transfer
histories (not yet tested).

"Production" status only applies to the single birth-death model.
"""

from cogent.seqsim.tree import RangeNode
from random import random

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Mike Robeson"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class ExtinctionError(Exception): pass
class TooManyTaxaError(Exception): pass

class BirthDeathModel(object):
    """Creates trees using a birth-death model.

    Initialize with timestep, birth prob, death prob.

    This class only produces the trees (using RangeNode); these trees already
    know how to evolve sequences, set rate matrices, etc.

    The trees returned will include lengths for each branch.

    WARNING: Sometimes, the ancestral node will die off, or that all the nodes
    will die off later. If that happens, an ExtinctionError will be raised.
    """
    NodeClass = RangeNode

    def __init__(self, BirthProb, DeathProb, TimePerStep, ChangedBirthProb=None, \
            ChangedDeathProb=None, ChangedBirthStep=None, ChangedDeathStep=None, \
            MaxStep=1000, MaxTaxa=None):
        """Returns a new BirthDeathModel object.

        BirthProb: probability that a node will split in timestep.
        DeathProb: probability that a node will die in timestep.
        TimePerStep: branch length (sequence distance units) for each step.
        ChangedBirthProb: new BirthProb to be set at the time step specified
                            in ChangedBirthStep.
        ChangedDeathProb: new DeathProb to be set at the time step specified
                            in ChangedDeathStep.
        ChangedBirthStep: time step at which the ChangedBirthProb is set.
        ChangedDeathStep: time step at which the ChangedDeathProb is set.
        MaxStep: maximum time step before stopping. Default 1000.
        MaxTaxa: maximum taxa before stopping. Default None.
        
        Sets CurrStep to 0 at the beginning to measure elapsed time.

        Note: if both a birth and a death occur in the same timestep, they
        will be ignored.

        WARNING: If neither MaxStep nor MaxTaxa is set, the simulation will keep
        going until all nodes are extinct, or you run out of memory!
        """
        self.CurrStep = 0
        self.BirthProb = BirthProb
        self.DeathProb = DeathProb 
        self.TimePerStep = float(TimePerStep)
        self.ChangedBirthProb = self.prob_check(ChangedBirthProb)
        self.ChangedBirthStep = self.step_check(ChangedBirthStep)
        self.ChangedDeathProb = self.prob_check(ChangedDeathProb)
        self.ChangedDeathStep = self.step_check(ChangedDeathStep)
        self.CurrDeathProb = DeathProb
        self.CurrBirthProb = BirthProb
        self.MaxStep = MaxStep
        self.MaxTaxa = MaxTaxa
        if TimePerStep <= 0:
            raise ValueError, "TimePerStep must be greater than zero"
        if not (0 <= BirthProb <= 1) or not (0 <= DeathProb <= 1):
            raise ValueError, "Birth and death probs must be between 0 and 1"
        #self.CurrStep = 0
        self.Tree = self.NodeClass()
        self.CurrTaxa = [self.Tree]

    def prob_check(self,prob_value):
        """Checks if prabability value lies between 0 and 1"""
        if prob_value is not None:
            if not (0 <= prob_value <= 1):
                raise ValueError, "\'prob_value\'  must be between 0 and 1"
            else:
                return prob_value
        else:
            return None
        
    def step_check(self,step_value):
        """Checks to see if value is greater than zero"""
        if step_value is not None:
            if step_value <= 0:
                raise ValueError, "\'stop_value\' must be greater than zero"
            else:
                return step_value
        else:
            return None
    
    def timeOk(self):
        """Return True only if the maximum time has not yet been reached."""
        #If MaxStep is not set, never say that the maximum time was reached
        if self.MaxStep is None:
            return True
        else:
            return self.CurrStep < self.MaxStep

    def taxaOk(self):
        """Returns True if the number of taxa is > 0 and < self.MaxTaxa.
        
        Note: MaxTaxa is exclusive (i.e. if MaxTaxa is 32, taxaOk will return
        False when the number of taxa is exactly 32, allowing you to stop when
        this number is reached).
        """
        num_taxa = len(self.CurrTaxa)
        if num_taxa < 1:
            return False
        if self.MaxTaxa is not None:
            return num_taxa < self.MaxTaxa
        #otherwise, if self.MaxTaxa was not set, any number is OK since we 
        #know we have at least item in the list or we wouldn't have got here.
        else:
            return True
    
    def B_Prob(self):
        """Checks to see if Birth probability changes during a time step.

        If target time step is defined and met, the new birth prob will 
        take affect from that point onward.
        """
        if self.ChangedBirthStep is None:
            return self.BirthProb
        else:
            if self.CurrStep < self.ChangedBirthStep:
                return self.BirthProb
            else:
                try:
                    self.CurrBirthProb = self.ChangedBirthProb
                    return self.ChangedBirthProb
                except ValueError: 'ChangedBirthProb value is \'None\''

    def D_Prob(self):
        """Checks to see if Death probability changes during a time step.

        If target time step is defined and met, the new death prob will 
        take affect from that point onward.
"""
        if self.ChangedDeathStep is None:
            return self.DeathProb
        else:
            if self.CurrStep < self.ChangedDeathStep:
                return self.DeathProb
            else:
                try:
                    self.CurrDeathProb = self.ChangedDeathProb
                    return self.ChangedDeathProb
                except ValueError: 'ChangedDeathProb value is \'None\''

    def step(self, random_f=random):
        """Advances the state of the object by one timestep.

        Specifically:
        
        For each node in the current taxa, decides whether it's going to
        produce a birth or a death.
        
        If a node dies, delete it from the list of current taxa.

        If a node gives birth, add two child nodes to the list of current
        taxa each with branchlength equal to the timestep, and delete the
        original node from the list of taxa.

        Otherwise, add the timestep to the node's branchlength.
        """
        #create list of new current nodes
        b = self.B_Prob()
        d = self.D_Prob()
        nc = self.NodeClass
        ts = self.TimePerStep
        new_list = []
        for node in self.CurrTaxa:
            died = random_f() < d
            born = random_f() < b
            #need to duplicate only if it was born and one didn't die
            if (born and not died):
                first_child = nc()
                second_child = nc()
                children = [first_child, second_child]
                #remember, we need to take care of both parent and child 
                #refs manually unless the tree class does it for us
                node.Children = children
                first_child.Parent = node
                second_child.Parent = node
                new_list.extend(children)
            elif (died and not born):
                #don't add the dead node to the new list
                continue
            else:   #i.e. if born and died, or if nothing happened
                new_list.append(node)
        #update time steps
        for node in new_list:
            if hasattr(node, 'Length') and node.Length is not None:
                node.Length += ts
            else:
                node.Length = ts
        self.CurrStep += 1
        #set the list of current nodes to the new list
        self.CurrTaxa = new_list

    def __call__(self, filter=True, exact=False, random_f=random):
        """Returns a new tree using params in self.
        
        If filter is True (the default), gets rid of extinct lineages.

        If exact is True (default is False), raises exception if we didn't
        get the right number of taxa

        WARNING: Because multiple births can happen in the same timestep, you
        might get more than the number of taxa you specify. Check afterwards!
        """
        self.CurrStep = 0
        self.Tree = self.NodeClass()
        self.CurrTaxa = [self.Tree]
        while 1:
            self.step(random_f)
            if not(self.timeOk() and self.taxaOk()):
                break

        if not self.CurrTaxa:
            raise ExtinctionError, "All taxa are extinct."
        if filter:
            self.Tree.filter(self.CurrTaxa, keep=True)
        if exact and self.MaxTaxa and (len(self.CurrTaxa) != self.MaxTaxa):
            raise TooManyTaxaError, "Got %s taxa, not %s." % \
                (len(self.CurrTaxa), self.MaxTaxa)
        return self.Tree

class GeneNode(RangeNode):
    """Holds a phylogenetic node that corresponds to a gene.

    Specificially, needs Species property holding ref to its species.

    WARNING: the current implementation does not take Species in __init__,
    but assumes you will create it manually after you make the node.
    """
    pass

class SpeciesNode(RangeNode):
    """Holds a phylogenetic node that corresponds to a species.

    Specifically, needs a Genes property that holds refs to its genes.
    
    WARNING: the current implementation does not take Genes in __init__,
    but assumes you will create it manually after you make the node.
    """
    pass

class DoubleBirthDeathModel(object):
    """Creates species and gene trees using a double birth-death model.

    Initialize with timestep, birth prob, death prob.

    This class only produces the trees (using RangeNode); these trees already
    know how to evolve sequences, set rate matrices, etc.

    The trees returned will include lengths for each branch.

    WARNING: Sometimes, the ancestral node will die off, or that all the nodes
    will die off later. If that happens, an ExtinctionError will be raised.
    """
    GeneClass = GeneNode
    SpeciesClass = SpeciesNode

    def __init__(self, GeneBirth, GeneDeath, GeneTransfer, SpeciesBirth, \
            SpeciesDeath, SpeciesRateChange, TimePerStep, GenesAtStart, \
            MaxStep=1000, MaxGenes=None, MaxSpecies=None, MaxGenome=None,
            DEBUG=False):
        """Returns a new BirthDeathModel object.

        GeneBirth: f(val, gene) returning True if gene is born given seed val.
        
        GeneDeath: f(val, gene) returning True if gene dies given seed val.
        
        GeneTransfer: f(val, gene) returning species that gene transfers to (or 
        None if no trransfer.)

        SpeciesBirth: f(val, species) returning True if species splits.
        
        SpeciesDeath: f(val, species) returning True if species dies.
        
        SpeciesRateChange: f(val, species) resetting species rate matrix given
        val.
        
        NOTE: in current implementation, Q only changes when species duplicates.
        
        TimePerStep: branch length (sequence distance units) for each step.
        GenesAtStart: number of genes at the beginning of the simulation.
        
        MaxStep: maximum time step before stopping. Default 1000.
        MaxGenes: maximum genes before stopping. Default None.
        MaxSpecies: maximum species before stopping. Default None.
        MaxGenome: maximum number of genes in a genome. Default None.
        
        Sets CurrStep to 0 at the beginning to measure elapsed time.

        Note: if both a birth and a death occur in the same timestep, they
        will be ignored.

        WARNING: If neither MaxStep nor MaxTaxa is set, the simulation will keep
        going until all nodes are extinct, or you run out of memory!
        """
        self.GeneBirth = GeneBirth 
        self.GeneDeath = GeneDeath
        self.GeneTransfer = GeneTransfer
        self.SpeciesBirth = SpeciesBirth
        self.SpeciesDeath = SpeciesDeath
        self.SpeciesRateChange = SpeciesRateChange
        self.TimePerStep = TimePerStep
        self.GenesAtStart = GenesAtStart
        self.MaxStep = MaxStep
        self.MaxGenes = MaxGenes
        self.MaxSpecies = MaxSpecies
        self.MaxGenome = MaxGenome
        self.DEBUG = DEBUG
        if TimePerStep <= 0:
            raise ValueError, "TimePerStep must be greater than zero"
        self._init_vars()

    def _init_vars(self):
        """Initialize vars before running the simulation."""
        self.CurrStep = 0
        self.SpeciesTree = self.SpeciesClass()
        self.SpeciesTree.Length = 0
        self.SpeciesTree.BirthDeathModel = self
        self.CurrSpecies = [self.SpeciesTree]
        self.SpeciesTree.CurrSpecies = self.CurrSpecies #ref to same object
        self.GeneTrees = [self.GeneClass() for i in range(self.GenesAtStart)]
        for i in self.GeneTrees:
            i.Length = 0
            i.BirthDeathModel = self
        self.CurrGenes = self.GeneTrees[:]

        #set gene/species references
        for i in self.CurrGenes:
            i.Species = self.SpeciesTree
        self.SpeciesTree.Genes = self.CurrGenes[:]  
        #note: copy of CurrGenes list, not reference

    def timeOk(self):
        """Return True only if the maximum time has not yet been reached."""
        #If MaxStep is not set, never say that the maximum time was reached
        if self.MaxStep is None:
            return True
        else:
            return self.CurrStep < self.MaxStep

    def genesOk(self):
        """Returns True if the number of genes is > 0 and < self.MaxGenes.
        
        Note: MaxGenes is exclusive (i.e. if MaxGenes is 32, genesOk will return
        False when the number of genes is exactly 32, allowing you to stop when
        this number is reached).
        """
        num_taxa = len(self.CurrGenes)
        if num_taxa < 1:
            return False
        if self.MaxGenes is not None:
            return num_taxa < self.MaxGenes
        #otherwise, if self.MaxTaxa was not set, any number is OK since we 
        #know we have at least item in the list or we wouldn't have got here.
        else:
            return True
        
    def speciesOk(self):
        """Returns True if the number of species is > 0 and < self.MaxSpecies.
        
        Note: MaxSpecies is exclusive (i.e. if MaxSpecies is 32, speciesOk 
        will return False when the number of species is exactly 32, allowing 
        you to stop when this number is reached).
        """
        num_taxa = len(self.CurrSpecies)
        if num_taxa < 1:
            return False
        if self.MaxSpecies is not None:
            return num_taxa < self.MaxSpecies
        #otherwise, if self.MaxTaxa was not set, any number is OK since we 
        #know we have at least item in the list or we wouldn't have got here.
        else:
            return True

    def genomeOk(self):
        """Returns True if the max genome size is < self.MaxGenome.
        
        Note: MaxGenome is exclusive (i.e. if MaxGenome is 32, genomeOk will 
        return False when the max genome size is exactly 32, allowing you to 
        stop when this number is reached).
        """
        max_taxa = max([len(i.Genes) for i in self.CurrSpecies])
        if num_taxa < 1:
            return False
        if self.MaxGenome is not None:
            return num_taxa < self.MaxGenome
        #otherwise, if self.MaxTaxa was not set, any number is OK since we 
        #know we have at least item in the list or we wouldn't have got here.
        else:
            return True

    def geneStep(self, random_f=random):
        """Advances the state of the genes by one timestep (except speciation).

        Specifically:
        
        Decides whether each gene will die, duplicate, or transfer.
        
        If a gene dies, delete it from the list of current genes.

        If a gene gives birth, add two child nodes to the list of current
        genes each with zero branchlength (will increment later), and delete the
        original node from the list of genes.

        If a gene transfers, handle like birth but also change the species.

        WARNING: This method does not increment the branch length or the time 
        counter. Handle separately!
        """
        #create list of new current nodes
        #Too complex to do combinations of states. Use three-pass algorithm:
        #1. birth
        #2. transfer
        #3. death
        #i.e. each copy gets a separate chance at death after it is made.
        #note that this differs slightly from what we do in the single
        #birth-death model where each original gene gets a chance at death
        #and a death and a duplication just cancel. Is this a problem with
        #the original model?
        self._gene_birth_step(random_f)
        self._gene_transfer_step(random_f)
        self._gene_death_step(random_f)

    def _duplicate_gene(self, gene, orig_species, new_species=None, \
        new_species_2=None):
        """Duplicates a gene, optionally attaching to new species.

        When called with only orig_species, duplicates the gene in the same 
        species (killing the old gene and making two copies).
        
        When called with orig_species and new_species, kills the old gene and 
        puts one new child into each of the old and new species (i.e. for
        horizontal gene transfer).

        When called with orig_species, new_species, and new_species_2, kills
        the old gene and puts one new child into each of the two new species
        (i.e. for speciation where all genes duplicate into new species).
        
        WARNING: Does not update self.CurrGenes (so can use in loop, but
        must update self.CurrGenes manually)."""
        
        gc = self.GeneClass
        #make new children
        first_child, second_child = gc(), gc()
        children = [first_child, second_child]
        #update gene parent/child refs
        gene.Children = children
        first_child.Parent = gene
        second_child.Parent = gene
        #init branch lengths
        first_child.Length = 0
        second_child.Length = 0
        #update species refs
        #first, figure out which species to deal with
        if new_species is None:   #add both to orig species
            first_species = orig_species
            second_species = orig_species
        elif new_species_2 is None: #add first to orig, second to new_species
            first_species = orig_species
            second_species = new_species
        else:   #add to the two new species
            first_species = new_species
            second_species = new_species_2
        #then, update the refs
        first_child.Species = first_species
        second_child.Species = second_species
        orig_species.Genes.remove(gene)
        first_species.Genes.append(first_child)
        second_species.Genes.append(second_child)
        #return the new genes for appending or whatever
        return first_child, second_child

    def _gene_birth_step(self, random_f=random):
        """Implements gene birth sweep."""
        gb = self.GeneBirth
        new_genes = []
        for gene in self.CurrGenes:
            if gb(random_f(), gene):
                new_genes.extend(self._duplicate_gene(gene, gene.Species))
            else:
                new_genes.append(gene)
        self.CurrGenes[:] = new_genes[:]

    def _gene_transfer_step(self, random_f=random):
        """Implements gene transfer sweep."""
        gt = self.GeneTransfer
        new_genes = []
        #step 2: transfer
        for gene in self.CurrGenes:
            new_species = gt(random_f(), gene)
            if new_species:
                new_genes.extend(self._duplicate_gene(gene, gene.Species, \
                    new_species))
            else:
                new_genes.append(gene)
        self.CurrGenes[:] = new_genes[:]

    def _gene_death_step(self, random_f=random):
        """Implements gene death sweep."""
        gd = self.GeneDeath
        new_genes = []
        for gene in self.CurrGenes:
            if gd(random_f(), gene):
                gene.Species.Genes.remove(gene)
            else:
                new_genes.append(gene)
        self.CurrGenes[:] = new_genes[:]

    def speciesStep(self, random_f=random):
        """Advances the state of the species by one timestep.

        Specifically:
        
        For each species in the current species, decides whether it's going to
        produce a birth or a death.
        
        If a species dies, delete it from the list of current species.

        If a species gives birth, add two child nodes to the list of current
        species, duplicate all their genes, and delete the
        original node from the list of taxa.
        Otherwise, add the timestep to the node's branchlength.
        """
        #make the species that are going to duplicate
        self._species_birth_step(random_f)
        #kill the species that are going to die
        self._species_death_step(random_f)
        self._kill_orphan_genes_step()

    def _species_death_step(self, random_f):
        """Kills species in self."""
        sd = self.SpeciesDeath
        new_list = []
        for s in self.CurrSpecies:
            if not sd(random_f(), s):
                new_list.append(s)
        self.CurrSpecies[:] = new_list[:]

    def _kill_orphan_genes_step(self):
        """Kills genes whose species has been removed."""
        new_list = []
        species_dict = dict.fromkeys(map(id, self.CurrSpecies))
        for g in self.CurrGenes:
            if id(g.Species) in species_dict:
                new_list.append(g)
        self.CurrGenes[:] = new_list[:]

    def _species_birth_step(self, random_f):
        sb = self.SpeciesBirth
        new_list = []
        for s in self.CurrSpecies:
            if sb(random_f(), s):
                new_list.extend(self._duplicate_species(s))
            else:
                new_list.append(s)
        self.CurrSpecies[:] = new_list[:]

    def _duplicate_species(self, species):
        """Duplicates a species by duplicating all its genes.
        
        WARNING: Doesn't remove from self.CurrSpecies: must do
        outside function (so can use while iterating over self.CurrSpecies).
        """
        for i in self.CurrGenes:
            assert i.Species in self.CurrSpecies
        for i in self.CurrSpecies:
            assert not i.Children
        if self.DEBUG:
            print '*** DUPLICATING SPECIES'
            print "SPECIES GENES AT START: ", len(species.Genes)
        sc = self.SpeciesClass
        #make new species
        first_child, second_child = sc(), sc()
        children = [first_child, second_child]
        #update species parent/child refs
        species.Children = children
        first_child.Parent = species
        second_child.Parent = species
        #update other child properties
        first_child.CurrSpecies = species.CurrSpecies
        first_child.Genes = []
        first_child.Length = 0
        second_child.CurrSpecies = species.CurrSpecies
        second_child.Genes = []
        second_child.Length = 0
        #update gene references
        curr_genes = self.CurrGenes
        if self.DEBUG:
            print "GENES BEFORE SWEEP: ", len(curr_genes)
            print "NUM GENES IN SPECIES: ", len(species.Genes)
        
        gene_counter = 0
        for gene in species.Genes[:]:
            assert gene.Species is species
            if self.DEBUG:
                print "handling gene ", gene_counter
            gene_counter += 1
            curr_genes.remove(gene)
            assert gene not in curr_genes
            for i in curr_genes: assert (i.Species in self.CurrSpecies) or \
                i.Species in [first_child, second_child]
            curr_genes.extend(self._duplicate_gene(gene, \
                gene.Species, first_child, second_child))
            for i in curr_genes: assert (i.Species in self.CurrSpecies) or \
                i.Species in [first_child, second_child]
        if self.DEBUG:
            print "GENES IN FIRST CHILD: ", len(first_child.Genes)
            print "GENES IN SECOND CHILD: ", len(second_child.Genes)
            print "GENES AFTER SWEEP: ", len(curr_genes)
        self.SpeciesTree.assignIds()
        if self.DEBUG:
            print "SPECIES TREE: ", self.SpeciesTree
        if self.DEBUG:
            print "SPECIES ASSIGNMENTS FOR EACH GENE"
        for i in curr_genes: 
            if self.DEBUG:
                print i.Species.Id
            assert (i.Species in self.CurrSpecies) or i.Species in [first_child, second_child]
        return children

    def updateLengths(self):
        """Adds timestep to the branch lengths of surviving genes/species."""
        ts = self.TimePerStep
        for gene in self.CurrGenes:
            gene.Length += ts
        for species in self.CurrSpecies:
            species.Length += ts

    def __call__(self, filter=True, exact_species=False, exact_genes=False, \
            random_f=random):
        """Returns a new tree using params in self.
        
        If filter is True (the default), gets rid of extinct lineages.

        If exact is True (default is False), raises exception if we didn't
        get the right number of taxa

        WARNING: Because multiple births can happen in the same timestep, you
        might get more than the number of taxa you specify. Check afterwards!
        """
        self._init_vars()
        done = False
        while not done:
            if self.DEBUG:
                print "CURR STEP:", self.CurrStep
            for i in self.CurrGenes: assert i.Species in self.CurrSpecies
            self.geneStep(random_f)
            for i in self.CurrGenes: assert i.Species in self.CurrSpecies
            self.speciesStep(random_f)
            for i in self.CurrGenes: assert i.Species in self.CurrSpecies
            self.updateLengths()
            self.CurrStep += 1
            done = not(self.timeOk() and self.genesOk() and self.speciesOk \
                and self.genomeOk)
        #check if all the constraints were met
        if not (self.CurrSpecies or self.CurrGenes):
            raise ExtinctionError, "All taxa are extinct."
        if exact_species and self.MaxSpecies and \
            (len(self.CurrSpecies) != self.MaxSpecies):
            raise TooManyTaxaError, "Got %s species, not %s." % \
                (len(self.CurrSpecies), self.MaxSpecies)
        if exact_genes and self.MaxGenes and \
            (len(self.CurrGenes) != self.MaxGenes):
            raise TooManyTaxaError, "Got %s genes, not %s." % \
                (len(self.CurrGenes), self.MaxGenes)
        #filter if required
        if filter:
            if self.DEBUG:
                print "***FILTERING..."
            self.SpeciesTree.assignIds()
            if self.DEBUG:
                print "BEFORE PRUNE: ", self.SpeciesTree
            self.SpeciesTree.filter(self.CurrSpecies, keep=True)
            if self.DEBUG:
                print "AFTER PRUNE: ", self.SpeciesTree
            for i, t in enumerate(self.GeneTrees):
                t.filter(self.CurrGenes)
        return self.SpeciesTree, self.GeneTrees
