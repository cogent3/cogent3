#!/usr/bin/env python

import numpy
import warnings
import substitution_calculation
from cogent.evolve.likelihood_tree import makeLikelihoodTreeLeaf

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def makeModel(mprob_model, tuple_alphabet, mask):        
    if mprob_model == "monomers":
        return PosnSpecificMonomerProbModel(tuple_alphabet, mask)
    elif mprob_model == "monomer":
        return MonomerProbModel(tuple_alphabet, mask)
    elif mprob_model == "conditional":
        return ConditionalMotifProbModel(tuple_alphabet, mask)
    elif mprob_model in ["word", "tuple", None]:
        return SimpleMotifProbModel(tuple_alphabet)
    else:
        raise ValueError("Unknown mprob model '%s'" % str(mprob_model))
    
class MotifProbModel(object):    
    def __init__(self, *whatever, **kw):
        raise NotImplementedError
        
    def calcWordProbs(self, *monomer_probs):  
        assert len(monomer_probs) == 1
        return monomer_probs[0]
        
    def calcWordWeightMatrix(self, *monomer_probs):  
        assert len(monomer_probs) == 1
        return monomer_probs[0]

    def makeMotifProbsDefn(self):
        """Makes the first part of a parameter controller definition for this
        model, the calculation of motif probabilities"""
        return substitution_calculation.PartitionDefn(
                name="mprobs", default=None, dimensions = ('locus','edge'),
                dimension=('motif', tuple(self.getInputAlphabet())))
        
    def setParamControllerMotifProbs(self, pc, motif_probs, **kw):
        pc.setParamRule('mprobs', value=motif_probs, **kw)
    
    def countMotifs(self, alignment, include_ambiguity=False, recode_gaps=True):
        result = None
        for seq_name in alignment.getSeqNames():
            sequence = alignment.getGappedSeq(seq_name, recode_gaps)
            leaf = makeLikelihoodTreeLeaf(sequence, self.getCountedAlphabet(), 
                    seq_name)
            count = leaf.getMotifCounts(include_ambiguity=include_ambiguity)
            if result is None:
                result = count.copy()
            else:
                result += count
        return result
        
    def adaptMotifProbs(self, motif_probs, auto=False):
        motif_probs = self.getInputAlphabet().adaptMotifProbs(motif_probs)
        assert abs(sum(motif_probs)-1.0) < 0.0001, motif_probs
        return motif_probs

    def makeEqualMotifProbs(self):
        alphabet = self.getInputAlphabet()
        p = 1.0/len(alphabet)
        return dict([(m,p) for m in alphabet])
        
    def makeSampleMotifProbs(self):
        import random
        motif_probs = numpy.array(
            [random.uniform(0.2, 1.0) for m in self.getCountedAlphabet()])
        motif_probs /= sum(motif_probs)
        return motif_probs

        
class SimpleMotifProbModel(MotifProbModel):
    def __init__(self, alphabet):
        self.alphabet = alphabet
        
    def getInputAlphabet(self):
        return self.alphabet
        
    def getCountedAlphabet(self):
        return self.alphabet
        
    def makeMotifWordProbDefns(self):
        monomer_probs = self.makeMotifProbsDefn()
        return (monomer_probs, monomer_probs, monomer_probs)       
        

class ComplexMotifProbModel(MotifProbModel):
    def __init__(self, tuple_alphabet, mask):
        """Arguments:
            - tuple_alphabet: series of multi-letter motifs
            - monomers: the monomers from which the motifs are made
            - mask: instantaneous change matrix"""
        self.mask = mask
        self.tuple_alphabet = tuple_alphabet
        self.monomer_alphabet = monomers = tuple_alphabet.MolType.Alphabet
        self.word_length = length = tuple_alphabet.getMotifLen()
        size = len(tuple_alphabet)

        # m2w[AC, 1] = C
        # w2m[0, AC, A] = True
        # w2c[ATC, AT*] = 1
        self.m2w = m2w = numpy.zeros([size, length], int)
        self.w2m = w2m = numpy.zeros([length, size, len(monomers)], int)
        contexts = monomers.getWordAlphabet(length-1)
        self.w2c = w2c = numpy.zeros([size, length*len(contexts)], int)
        for (i, word) in enumerate(tuple_alphabet):
            for j in range(length):
                monomer = monomers.index(word[j])
                context = contexts.index(word[:j]+word[j+1:])
                m2w[i, j] = monomer
                w2m[j, i, monomer] = 1
                w2c[i, context*length+j] = 1
            
        self.mutated_posn = numpy.zeros(mask.shape, int)
        self.mutant_motif = numpy.zeros(mask.shape, int)
        self.context_indices = numpy.zeros(mask.shape, int)
        
        for (i, old_word, j, new_word, diff) in self._mutations():
            self.mutated_posn[i,j] = diff
            mutant_motif = new_word[diff]
            context = new_word[:diff]+new_word[diff+1:]
            self.mutant_motif[i,j] = monomers.index(mutant_motif)
            c = contexts.index(context)
            self.context_indices[i,j] = c * length + diff 
    
    def _mutations(self):
        diff_pos = lambda x,y: [i for i in range(len(x)) if x[i] != y[i]]
        num_states = len(self.tuple_alphabet)
        for i in range(num_states):
            old_word = self.tuple_alphabet[i]
            for j in range(num_states):
                new_word = self.tuple_alphabet[j]
                if self.mask[i,j]:
                    assert self.mask[i,j] == 1.0
                    diffs = diff_pos(old_word, new_word)
                    assert len(diffs) == 1, (old_word, new_word)
                    diff = diffs[0]
                    yield i, old_word, j, new_word, diff
        

class MonomerProbModel(ComplexMotifProbModel):
        
    def getInputAlphabet(self):
        return self.monomer_alphabet
    
    def getCountedAlphabet(self):
        return self.monomer_alphabet
        
    def calcMonomerProbs(self, word_probs):
        monomer_probs = numpy.dot(word_probs, self.w2m.sum(axis=0))
        monomer_probs /= monomer_probs.sum()
        return monomer_probs
        
    def calcWordProbs(self, monomer_probs):  
        result = numpy.product(monomer_probs.take(self.m2w), axis=-1)
        # maybe simpler but slower, works ok:
        #result = numpy.product(monomer_probs ** (w2m, axis=-1))
        result /= result.sum()
        return result
        
    def calcWordWeightMatrix(self, monomer_probs):  
        result = monomer_probs.take(self.mutant_motif) * self.mask
        return result

    def makeMotifWordProbDefns(self):
        monomer_probs = self.makeMotifProbsDefn()
        word_probs = substitution_calculation.CalcDefn(
                self.calcWordProbs, name="wprobs")(monomer_probs)
        mprobs_matrix = substitution_calculation.CalcDefn(
                self.calcWordWeightMatrix, name="mprobs_matrix")(monomer_probs)
        return (monomer_probs, word_probs, mprobs_matrix)

    def adaptMotifProbs(self, motif_probs, auto=False):
        try:
            motif_probs = self.monomer_alphabet.adaptMotifProbs(motif_probs)
        except ValueError:
            motif_probs = self.tuple_alphabet.adaptMotifProbs(motif_probs)
            if not auto:
                warnings.warn('Motif probs overspecified', stacklevel=5)
            motif_probs = self.calcMonomerProbs(motif_probs)
        return motif_probs
        

class PosnSpecificMonomerProbModel(MonomerProbModel):
    def getCountedAlphabet(self):
        return self.tuple_alphabet
        
    def calcPosnSpecificMonomerProbs(self, word_probs):
        monomer_probs = numpy.dot(word_probs, self.w2m)
        monomer_probs /= monomer_probs.sum(axis=1)[..., numpy.newaxis]
        return list(monomer_probs)

    def calcWordProbs(self, monomer_probs):
        positions = range(self.word_length)
        assert len(monomer_probs) == self.m2w.shape[1], (
            len(monomer_probs), type(monomer_probs), self.m2w.shape)
        result = numpy.product(
            [monomer_probs[i].take(self.m2w[:,i]) 
            for i in positions], axis=0)
        result /= result.sum()
        return result
    
    def calcWordWeightMatrix(self, monomer_probs):  
        positions = range(self.word_length)
        monomer_probs = numpy.array(monomer_probs) # so [posn, motif]
        size = monomer_probs.shape[-1]
        # should be constant
        extended_indices = self.mutated_posn * size + self.mutant_motif
        #print size, self.word_length
        #for a in [extended_indices, self.mutated_posn, self.mutant_motif, 
        #        monomer_probs]:
        #    print a.shape, a.max()
             
        result = monomer_probs.take(extended_indices) * self.mask
        return result

    def makeMotifWordProbDefns(self):
        monomer_probs = substitution_calculation.PartitionDefn(
                name="psmprobs", default=None, 
                dimensions = ('locus', 'position', 'edge'),
                dimension=('motif', tuple(self.getInputAlphabet())))
        monomer_probs3 = monomer_probs.acrossDimension('position', [
            str(i) for i in range(self.word_length)])
        monomer_probs3 = substitution_calculation.CalcDefn(
                lambda *x:numpy.array(x), name='mprobs')(*monomer_probs3)
        word_probs = substitution_calculation.CalcDefn(
                self.calcWordProbs, name="wprobs")(monomer_probs3)
        mprobs_matrix = substitution_calculation.CalcDefn(
                self.calcWordWeightMatrix, name="mprobs_matrix")(
                    monomer_probs3)
        return (monomer_probs, word_probs, mprobs_matrix)

    def setParamControllerMotifProbs(self, pc, motif_probs, **kw):
        assert len(motif_probs) == self.word_length
        for (i,m) in enumerate(motif_probs):
            pc.setParamRule('psmprobs', value=m, position=str(i), **kw)

    def adaptMotifProbs(self, motif_probs, auto=False):
        try:
            motif_probs = self.monomer_alphabet.adaptMotifProbs(motif_probs)
        except ValueError:
            motif_probs = self.tuple_alphabet.adaptMotifProbs(motif_probs)
            motif_probs = self.calcPosnSpecificMonomerProbs(motif_probs)
        else:
            motif_probs = [motif_probs] * self.word_length
        return motif_probs

class ConditionalMotifProbModel(ComplexMotifProbModel):
    def getInputAlphabet(self):
        return self.tuple_alphabet
        
    def getCountedAlphabet(self):
        return self.tuple_alphabet

    def calcWordWeightMatrix(self, motif_probs): 
        context_probs = numpy.dot(motif_probs, self.w2c)
        context_probs[context_probs==0.0] = numpy.inf
        result = motif_probs / context_probs.take(self.context_indices)
        return result
        
    def makeMotifWordProbDefns(self):
        mprobs = self.makeMotifProbsDefn()
        mprobs_matrix = substitution_calculation.CalcDefn(
                self.calcWordWeightMatrix, name="mprobs_matrix")(mprobs)
        return (mprobs, mprobs, mprobs_matrix)

