#!/usr/bin/env python
"""
substitution_model.py

Contains classes for defining Markov models of substitution.
These classes depend on an Alphabet class member for defining the set
of motifs that each represent a state in the Markov chain. Examples of
a 'dna' type alphabet motif is 'a', and of a 'codon' type motif is'atg'.

By default all models include the gap motif ('-' for a 'dna' alphabet or
'---' for a 'codon' alphabet). This differs from software such as PAML,
where gaps are treated as ambiguituous states (specifically, as 'n'). The gap
motif state can be excluded from the substitution model using the method
excludeGapMotif(). It is recommended that to ensure the alignment and the
substitution model are defined with the same alphabet that modifications
are done to the substitution model alphabet and this instance is then given
to the alignment.

The model's substitution rate parameters are represented as a dictionary
with the parameter names as keys, and predicate functions as the values.
These predicate functions compare a pair of motifs, returning True or False.
Many such functions are provided as methods of the class. For instance,
the istransition method is pertinent to dna based models. This method returns
True if an 'a'/'g' or 'c'/'t' pair is passed to it, False otherwise. In this
way the positioning of parameters in the instantaneous rate matrix (commonly
called Q) is determined.

>>> model = Nucleotide(equal_motif_probs=True)
>>> model.setparameterrules({'alpha': model.istransition})
>>> parameter_controller = model.makeParamController(tree)
"""

import numpy
from numpy.linalg import svd
import warnings
import inspect

from cogent.core import moltype
from cogent.evolve import parameter_controller, predicate, motif_prob_model
from cogent.evolve.substitution_calculation import (
    SubstitutionParameterDefn as ParamDefn, 
    RateDefn, LengthDefn, ProductDefn, CallDefn, CalcDefn,
    PartitionDefn, NonParamDefn, AlignmentAdaptDefn, ExpDefn, 
    ConstDefn, GammaDefn, MonotonicDefn, SelectForDimension, 
    WeightedPartitionDefn)
from cogent.evolve.discrete_markov import PsubMatrixDefn
from cogent.evolve.likelihood_tree import makeLikelihoodTreeLeaf
from cogent.maths.optimisers import ParameterOutOfBoundsError

__author__ = "Peter Maxwell, Gavin Huttley and Andrew Butterfield"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Gavin Huttley", "Andrew Butterfield", "Peter Maxwell",
                    "Matthew Wakefield", "Brett Easton", "Rob Knight",
                    "Von Bing Yap"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def predicate2matrix(alphabet, pred, mask=None):
    """From a test like istransition() produce an MxM boolean matrix"""
    M = len(alphabet)
    result = numpy.zeros([M,M], int)
    for i in range(M):
        for j in range(M):
            if mask is None or mask[i,j]:
                result[i,j] = pred(alphabet[i], alphabet[j])
    return result

def redundancyInPredicateMasks(preds):
    # Calculate the nullity of the predicates.  If non-zero
    # there is some redundancy and the model will be overparameterised.
    if len(preds) <= 1:
        return 0
    eqns = 1.0 * numpy.array([list(mask.flat) for mask in preds.values()])
    svs = svd(eqns)[1]
    # count non-duplicate non-zeros singular values
    matrix_rank = len([sv for sv in svs if abs(sv) > 1e-8])
    return len(preds) - matrix_rank

def _maxWidthIfTruncated(pars, delim, each):
    # 'pars' is an array of lists of strings, how long would the longest
    # list representation be if the strings were truncated at 'each'
    # characters and joined together with 'delim'.
    return max([
            sum([min(len(par), each) for par in par_list])
            + len(delim) * (len(par_list)-1)
        for par_list in pars.flat])

def _isSymmetrical(matrix):
    return numpy.alltrue(numpy.alltrue(matrix == numpy.transpose(matrix)))

def extend_docstring_from(cls, pre=False): 
    def docstring_inheriting_decorator(fn):
        parts = [getattr(cls,fn.__name__).__doc__, fn.__doc__ or '']
        if pre: parts.reverse()
        fn.__doc__ = ''.join(parts) 
        return fn 
    return docstring_inheriting_decorator

class _SubstitutionModel(object):
    # Subclasses must provide
    #  .makeParamControllerDefns()
    
    def __init__(self, alphabet, 
            motif_probs=None, optimise_motif_probs=False,
            equal_motif_probs=False, motif_probs_from_data=None,
            motif_probs_alignment=None, mprob_model=None,  
            model_gaps=False, recode_gaps=False, motif_length=None,
            name="", motifs=None):
        # subclasses can extend this incomplete docstring
        """
         
        Alphabet:
         - alphabet - An Alphabet object
         - motif_length: Use a tuple alphabet based on 'alphabet'.
         - motifs: Use a subalphabet that only contains those motifs.
         - model_gaps: Whether the gap motif should be included as a state.
         - recode_gaps: Whether gaps in an alignment should be treated as an 
           ambiguous state instead.
           
        Motif Probability:
         - motif_probs: Dictionary of probabilities.
         - equal_motif_probs: Flag to set alignment motif probs equal.
         - motif_probs_alignment: An alignment from which motif probs are set.
        
         If none of these options are set then motif probs will be derived 
         from the data: ie the particular alignment provided later.
          
         - optimise_motif_probs: Treat like other free parameters.  Any values
           set by the other motif_prob options will be used as initial values.
          
         - mprob_model: 'tuple', 'conditional' or 'monomer' to specify how
           tuple-alphabet (including codon) motif probs are used.
           
        """
            
        # MISC
        assert len(alphabet) < 65, "Alphabet too big. Try explicitly "\
            "setting alphabet to PROTEIN or DNA"
        
        self.name = name
        self._optimise_motif_probs = optimise_motif_probs
        
        # ALPHABET        
        if recode_gaps:
            if model_gaps:
                warnings.warn("Converting gaps to wildcards AND modeling gaps")
            else:
                model_gaps = False
        
        self.recode_gaps = recode_gaps
        
        self.MolType = alphabet.MolType
        if model_gaps:
            alphabet = alphabet.withGapMotif()
                
        if motif_length > 1:
            alphabet = alphabet.getWordAlphabet(motif_length)
        
        if motifs is not None:
            alphabet = alphabet.getSubset(motifs)
        self.alphabet = alphabet
        self.gapmotif = alphabet.getGapMotif()
        self._word_length = alphabet.getMotifLen()
        
        # MOTIF PROB ALPHABET MAPPING
        if mprob_model is None:
            mprob_model = 'tuple' if self._word_length==1 else 'conditional'
        elif mprob_model == 'word':
            mprob_model = 'tuple'
        
        if model_gaps and mprob_model != 'tuple':
            raise ValueError("mprob_model must be 'tuple' to model gaps")
        
        isinst = self._isInstantaneous
        self._instantaneous_mask = predicate2matrix(self.alphabet, isinst)
        self._instantaneous_mask_f = self._instantaneous_mask * 1.0
        self.mprob_model = motif_prob_model.makeModel(mprob_model, alphabet, 
                self._instantaneous_mask_f)
        
        # MOTIF PROBS        
        if equal_motif_probs:
            assert not (motif_probs or motif_probs_alignment), \
                    "Motif probs equal or provided but not both"
            motif_probs = self.mprob_model.makeEqualMotifProbs()
        elif motif_probs_alignment is not None:
            assert not motif_probs, \
                    "Motif probs from alignment or provided but not both"
            motif_probs = self.countMotifs(motif_probs_alignment)
            motif_probs = motif_probs.astype(float) / sum(motif_probs)
            assert len(alphabet) == len(motif_probs)
            motif_probs = dict(zip(alphabet, motif_probs))
        if motif_probs:
            self.adaptMotifProbs(motif_probs) # to check
            self.motif_probs = motif_probs
            if motif_probs_from_data is None:
                motif_probs_from_data = False
        else:
            self.motif_probs = None
            if motif_probs_from_data is None:
                motif_probs_from_data = True
        self.motif_probs_from_align = motif_probs_from_data
         
    def getParamList(self):
        return []
    
    def __str__(self):
        s = ["\n%s (" % self.__class__.__name__ ]
        s.append("name = '%s'; type = '%s';" %
                (getattr(self, "name", None), getattr(self, "type", None)))
        if hasattr(self, "predicate_masks"):
            parlist = self.predicate_masks.keys()
            s.append("params = %s;" % parlist)
        motifs = self.getMotifs()
        s.append("number of motifs = %s;" % len(motifs))
        s.append("motifs = %s)\n" % motifs)
        return " ".join(s)
    
    def getAlphabet(self):
        return self.alphabet
    
    def getMprobAlphabet(self):
        return self.mprob_model.getInputAlphabet()
    
    def getMotifs(self):
        return list(self.getAlphabet())
    
    def getWordLength(self):
        return self._word_length
    
    def getMotifProbs(self):
        """Return the dictionary of motif probabilities."""
        return self.motif_probs.copy()
    
    def setParamControllerMotifProbs(self, pc, mprobs, **kw):
        return self.mprob_model.setParamControllerMotifProbs(pc, mprobs, **kw)
    
    def makeLikelihoodFunction(self, tree, motif_probs_from_align=None,
            optimise_motif_probs=None, aligned=True, expm=None, digits=None,
            space=None, **kw):
        
        if motif_probs_from_align is None:
            motif_probs_from_align = self.motif_probs_from_align
        
        if optimise_motif_probs is None:
            optimise_motif_probs = self._optimise_motif_probs
        
        kw['optimise_motif_probs'] = optimise_motif_probs
        kw['motif_probs_from_align'] = motif_probs_from_align
        
        if aligned:
            klass = parameter_controller.AlignmentLikelihoodFunction
        else:
            alphabet = self.getAlphabet()
            assert alphabet.getGapMotif() not in alphabet
            klass = parameter_controller.SequenceLikelihoodFunction
        
        result = klass(self, tree, **kw)
        
        if self.motif_probs is not None:
            result.setMotifProbs(self.motif_probs, is_constant=
                not optimise_motif_probs, auto=True)
        
        if expm is None:
            expm = self._default_expm_setting
        if expm is not None:
            result.setExpm(expm)
        
        if digits or space:
            result.setTablesFormat(digits=digits, space=space)
        
        return result
    
    def makeParamController(self, tree, motif_probs_from_align=None,
            optimise_motif_probs=None, **kw):
        # deprecate
        return self.makeLikelihoodFunction(tree,
                motif_probs_from_align = motif_probs_from_align,
                optimise_motif_probs = optimise_motif_probs,
                **kw)
    
    def convertAlignment(self, alignment):
        # this is to support for everything but HMM
        result = {}
        for seq_name in alignment.getSeqNames():
            sequence = alignment.getGappedSeq(seq_name, self.recode_gaps)
            result[seq_name] = self.convertSequence(sequence, seq_name)
        return result
    
    def convertSequence(self, sequence, name):
        # makeLikelihoodTreeLeaf, sort of an indexed profile where duplicate
        # columns stored once, so likelihoods only calc'd once
        return makeLikelihoodTreeLeaf(sequence, self.getAlphabet(), name)
    
    def countMotifs(self, alignment, include_ambiguity=False):
        return self.mprob_model.countMotifs(alignment, 
                include_ambiguity, self.recode_gaps)
    
    def makeAlignmentDefn(self, model):
        align = NonParamDefn('alignment', ('locus',))
        # The name of this matters, it's used in likelihood_function.py
        # to retrieve the correct (adapted) alignment.
        return AlignmentAdaptDefn(model, align)
    
    def adaptMotifProbs(self, motif_probs, auto=False):
        return self.mprob_model.adaptMotifProbs(motif_probs, auto=auto)
    
    def calcMonomerProbs(self, word_probs):
        # Not presently used, always go monomer->word instead
        return self.mprob_model.calcMonomerProbs(word_probs)
    
    def calcWordProbs(self, monomer_probs):
        return self.mprob_model.calcWordProbs(monomer_probs)
    
    def calcWordWeightMatrix(self, monomer_probs):
        return self.mprob_model.calcWordWeightMatrix(monomer_probs)
    
    def makeParamControllerDefns(self, bin_names, endAtQd=False):
        (input_probs, word_probs, mprobs_matrix) = \
                self.mprob_model.makeMotifWordProbDefns()

        if len(bin_names) > 1:
            bprobs = PartitionDefn(
                [1.0/len(bin_names) for bin in bin_names], name = "bprobs",
                dimensions=['locus'], dimension=('bin', bin_names))
        else:
            bprobs = None
        
        defns = {
            'align': self.makeAlignmentDefn(ConstDefn(self, 'model')),
            'bprobs': bprobs,
            'word_probs': word_probs,
             }
        
        rate_params = self.makeRateParams(bprobs)
        if endAtQd:
            defns['Qd'] = self.makeQdDefn(word_probs, mprobs_matrix, rate_params)
        else:
            defns['psubs'] = self.makePsubsDefn(
                    bprobs, word_probs, mprobs_matrix, rate_params)
        return defns


class DiscreteSubstitutionModel(_SubstitutionModel):
    _default_expm_setting = None
    
    def _isInstantaneous(self, x, y):
        return True
    
    def getParamList(self):
        return []
        
    def makeRateParams(self, bprobs):
        return []

    def makePsubsDefn(self, bprobs, word_probs, mprobs_matrix, rate_params):
        assert len(rate_params) == 0
        assert word_probs is mprobs_matrix, "Must use simple mprob model"
        motifs = tuple(self.getAlphabet())
        return PsubMatrixDefn(
            name="psubs", dimension = ('motif', motifs), default=None, 
            dimensions=('locus', 'edge'))
        

class _ContinuousSubstitutionModel(_SubstitutionModel):
    # subclass must provide:
    #
    # - parameter_order: a list of parameter names corresponding to the 
    #   arguments of:
    # 
    # - calcExchangeabilityMatrix(*params)
    #   convert len(self.parameter_order) params to a matrix
    
    """A substitution model for which the rate matrix (P) is derived from an 
    instantaneous rate matrix (Q).  The nature of the parameters used to define 
    Q is up to the subclasses.  
    """

    # At some point this can be made variable, and probably
    # the default changed to False
    long_indels_are_instantaneous = True
    
    _scalableQ = True
    _exponentiator = None
    _default_expm_setting = 'either'
    
    @extend_docstring_from(_SubstitutionModel)
    def __init__(self, alphabet, with_rate=False, ordered_param=None, 
            distribution=None, partitioned_params=None, do_scaling=None, **kw):

        """
         - with_rate: Add a 'rate' parameter which varies by bin. 
         - ordered_param: name of a single parameter which distinguishes any bins.
         - distribution: choices of 'free' or 'gamma' or an instance of some
           distribution. Could probably just deprecate free
         - partitioned_params: names of params to be partitioned across bins
         - do_scaling: Scale branch lengths as the expected number of 
           substitutions.  Reduces the maximum substitution df by 1.
        """
        
        _SubstitutionModel.__init__(self, alphabet, **kw)
        alphabet = self.getAlphabet() # as may be altered by recode_gaps etc.
               
        if do_scaling is None:
            do_scaling = self._scalableQ
        if do_scaling and not self._scalableQ:
            raise ValueError("Can't autoscale a %s model" % type(self).__name__)
        self._do_scaling = do_scaling
        
        # BINS
        if not ordered_param:
            if ordered_param is not None:
                warnings.warn('ordered_param should be a string or None')
                ordered_param = None
            if distribution:
                if with_rate:
                    ordered_param = 'rate'
                else:
                    raise ValueError('distribution provided without ordered_param')      
        elif not isinstance(ordered_param, str):
            warnings.warn('ordered_param should be a string or None')
            assert len(ordered_param) == 1, 'More than one ordered_param'
            ordered_param = ordered_param[0]
            assert ordered_param, "False value hidden in list"
        self.ordered_param = ordered_param
        
        if distribution == "gamma":
            distribution = GammaDefn
        elif distribution in [None, "free"]:
            distribution = MonotonicDefn
        elif isinstance(distribution, basestring):
            raise ValueError('Unknown distribution "%s"' % distribution)
        self.distrib_class = distribution
        
        if not partitioned_params:
            partitioned_params = ()
        elif isinstance(partitioned_params, str):
            partitioned_params = (partitioned_params,)
        else:
            partitioned_params = tuple(partitioned_params)
        if self.ordered_param:
            if self.ordered_param not in partitioned_params:
                partitioned_params += (self.ordered_param,)
        self.partitioned_params = partitioned_params
        
        if 'rate' in partitioned_params:
            with_rate = True
        self.with_rate = with_rate  

        # CACHED SHORTCUTS
        self._exponentiator = None
        #self._ident = numpy.identity(len(self.alphabet), float)
        
    def checkParamsExist(self):
        """Raise an error if the parameters specified to be partitioned or
        ordered don't actually exist."""
        for param in self.partitioned_params:
            if param not in self.parameter_order and param != 'rate':
                desc = ['partitioned', 'ordered'][param==self.ordered_param]
                raise ValueError('%s param "%s" unknown' % (desc, param))

    def _isInstantaneous(self, x, y):
        diffs = sum([X!=Y for (X,Y) in zip(x,y)])
        return diffs == 1 or (diffs > 1 and
                self.long_indels_are_instantaneous and self._isAnyIndel(x, y))
    
    def _isAnyIndel(self, x, y):
        """An indel of any length"""
        # Things get complicated when a contigous indel of any length is OK:
        if x == y:
            return False
        gap_start = gap_end = gap_strand = None
        for (i, (X,Y)) in enumerate(zip(x,y)):
            G = self.gapmotif[i]
            if X != Y:
                if X != G and Y != G:
                    return False  # non-gap differences had their chance above
                elif gap_start is None:
                    gap_start = i
                    gap_strand = [X,Y].index(G)
                elif gap_end is not None or [X,Y].index(G) != gap_strand:
                    return False # can't start a second gap
                else:
                    pass # extend open gap
            elif gap_start is not None:
                gap_end = i
        return True
    
    def calcQ(self, word_probs, mprobs_matrix, *params):
        Q = self.calcExchangeabilityMatrix(word_probs, *params)
        Q *= mprobs_matrix
        row_totals = Q.sum(axis=1)
        Q -= numpy.diag(row_totals)
        if self._do_scaling:
            Q *= 1.0 / (word_probs * row_totals).sum()
        return Q
    
    def makeQdDefn(self, word_probs, mprobs_matrix, rate_params):
        """Diagonalized Q, ie: rate matrix prepared for exponentiation"""
        Q = CalcDefn(self.calcQ, name='Q')(word_probs, mprobs_matrix, *rate_params)
        expm = NonParamDefn('expm')
        exp = ExpDefn(expm)
        Qd = CallDefn(exp, Q, name='Qd')
        return Qd
    
    def _makeBinParamDefn(self, edge_par_name, bin_par_name, bprob_defn):
        # if no ordered param defined, behaves as old, everything indexed by and edge
        if edge_par_name not in self.partitioned_params:
            return ParamDefn(dimensions=['bin'], name=bin_par_name)
        
        if edge_par_name == self.ordered_param:
            whole = self.distrib_class(bprob_defn, bin_par_name)
        else:
            # this forces them to average to one, but no forced order
            # this means you can't force a param value to be shared across bins
            # so 1st above approach has to be used
            whole = WeightedPartitionDefn(bprob_defn, bin_par_name+'_partn')
        whole.bin_names = bprob_defn.bin_names
        return SelectForDimension(whole, 'bin', name=bin_par_name)
    
    def makeRateParams(self, bprobs):
        params = []
        for param_name in self.parameter_order:
            if bprobs is None or param_name not in self.partitioned_params:
                defn = ParamDefn(param_name)
            else:
                e_defn = ParamDefn(param_name, dimensions=['edge', 'locus'])
                # should be weighted by bprobs*rates not bprobs
                b_defn = self._makeBinParamDefn(
                        param_name, param_name+'_factor', bprobs)
                defn = ProductDefn(b_defn, e_defn, name=param_name+'_BE')
            params.append(defn)
        return params

    def makeFundamentalParamControllerDefns(self, bin_names):
        """Everything one step short of the psubs, because cogent.align code
        needs to handle Q*t itself."""
        defns = self.makeParamControllerDefns(bin_names, endAtQd=True)
        assert not 'length' in defns
        defns['length'] = LengthDefn()
        return defns
        
    def makePsubsDefn(self, bprobs, word_probs, mprobs_matrix, rate_params):
        distance = self.makeDistanceDefn(bprobs)
        P = self.makeContinuousPsubDefn(word_probs, mprobs_matrix, distance, rate_params)
        return P
        
    def makeDistanceDefn(self, bprobs):
        length = LengthDefn()
        if self.with_rate and bprobs is not None:
            b_rate = self._makeBinParamDefn('rate', 'rate', bprobs)
            distance = ProductDefn(length, b_rate, name="distance")
        else:
            distance = length
        return distance
        
    def makeContinuousPsubDefn(self, word_probs, mprobs_matrix, distance, rate_params):
        Qd = self.makeQdDefn(word_probs, mprobs_matrix, rate_params)
        P = CallDefn(Qd, distance, name='psubs')
        return P
    

class General(_ContinuousSubstitutionModel):
    """A continuous substitution model with one free parameter for each and 
    every possible instantaneous substitution."""
    
    # k = self.param_pick[i,j], 0<=k<=N+1
    # k==0:   not instantaneous, should be 0.0 in Q
    # k<=N:   apply Kth exchangeability parameter
    # k==N+1: no parameter, should be 1.0 in unscaled Q

    #@extend_docstring_from(_ContinuousSubstitutionModel)
    def __init__(self, alphabet, **kw):
        _ContinuousSubstitutionModel.__init__(self, alphabet, **kw)

        alphabet = self.getAlphabet() # as may be altered by recode_gaps etc.
        mask = self._instantaneous_mask
        N = len(alphabet)
        self.param_pick = numpy.zeros([N,N], int)
        self.parameter_order = []
        for (i,x) in enumerate(alphabet):
            for j in numpy.flatnonzero(mask[i]):
                y = alphabet[j]
                self.parameter_order.append('%s/%s'%(x,y))
                self.param_pick[i,j] = len(self.parameter_order)
        if self._do_scaling:
            const_param = self.parameter_order.pop()
        self.symmetric = False
        self.checkParamsExist()
        
    def calcExchangeabilityMatrix(self, mprobs, *params):
        return numpy.array((0.0,)+params+(1.0,)).take(self.param_pick)

class GeneralStationary(_ContinuousSubstitutionModel):
    """A continuous substitution model with one free parameter for each and 
    every possible instantaneous substitution, except the last in each column.
    As general as can be while still having stationary motif probabilities"""
    
    #@extend_docstring_from(_ContinuousSubstitutionModel)
    def __init__(self, alphabet, **kw):
        _ContinuousSubstitutionModel.__init__(self, alphabet, **kw)

        alphabet = self.getAlphabet() # as may be altered by recode_gaps etc.
        mask = self._instantaneous_mask
        N = len(alphabet)
        self.param_pick = numpy.zeros([N,N], int)
        self.parameter_order = []
        self.last_in_column = []
        for (d, (row, col)) in enumerate(zip(mask, mask.T)):
            row = list(numpy.flatnonzero(row[d:])+d)
            col = list(numpy.flatnonzero(col[d:])+d)
            if col:
                self.last_in_column.append((col.pop(), d))
            else:
                assert not row
            inst = [(d,j) for j in row] + [(i,d) for i in col]
            
            for (i, j) in inst:
                (x,y) = [alphabet[k] for k in [i,j]]
                self.parameter_order.append('%s/%s'%(x,y))
                self.param_pick[i,j] = len(self.parameter_order)
        if self._do_scaling:
            const_param = self.parameter_order.pop()
        self.symmetric = False
        self.checkParamsExist()

    def calcExchangeabilityMatrix(self, mprobs, *params):
        R = numpy.array((0.0,)+params+(1.0,)).take(self.param_pick)
        for (i,j) in self.last_in_column:
            assert i > j
            row_total = numpy.dot(mprobs, R[j])
            col_total = numpy.dot(mprobs, R[:,j])
            required = row_total - col_total
            if required < 0.0:
                raise ParameterOutOfBoundsError
            R[i,j] = required / mprobs[i]
        return R

class Empirical(_ContinuousSubstitutionModel):
    """A continuous substitution model with a predefined instantaneous rate 
    matrix."""
    
    @extend_docstring_from(_ContinuousSubstitutionModel)
    def __init__(self, alphabet, rate_matrix, **kw):
        """
         - rate_matrix: The instantaneous rate matrix
        """
        _ContinuousSubstitutionModel.__init__(self, alphabet, **kw)
    
        alphabet = self.getAlphabet()  # as may be altered by recode_gaps etc.
        N = len(alphabet)
        assert rate_matrix.shape == (N, N)
        assert numpy.alltrue(numpy.diagonal(rate_matrix) == 0)
        self._instantaneous_mask_f = rate_matrix * 1.0
        self._instantaneous_mask = (self._instantaneous_mask_f != 0.0)
        self.symmetric = _isSymmetrical(self._instantaneous_mask_f)
        self.parameter_order = []
        self.checkParamsExist()

    def calcExchangeabilityMatrix(self, mprobs):
        return self._instantaneous_mask_f.copy()


class SubstitutionModel(_ContinuousSubstitutionModel):
    """A continuous substitution model with only user-specified substitution
    parameters."""
    
    @extend_docstring_from(_ContinuousSubstitutionModel)
    def __init__(self, alphabet, predicates=None, scales=None, **kw):
        """
         - predicates: a dict of {name:predicate}. See cogent.evolve.predicate
         - scales: scale rules, dict with predicates
        """
        self._canned_predicates = None
        _ContinuousSubstitutionModel.__init__(self, alphabet, **kw)

        (predicate_masks, predicate_order) = self._adaptPredicates(predicates or [])

        # Check for redundancy in predicates, ie: 1 or more than combine
        # to be equivalent to 1 or more others, or the distance params.
        # Give a clearer error in simple cases like always false or true.
        for (name, matrix) in predicate_masks.items():
            if numpy.alltrue((matrix == 0).flat):
                raise ValueError("Predicate %s is always false." % name)
        predicates_plus_scale = predicate_masks.copy()
        predicates_plus_scale[None] = self._instantaneous_mask
        if self._do_scaling:
            for (name, matrix) in predicate_masks.items():
                if numpy.alltrue((matrix == self._instantaneous_mask).flat):
                    raise ValueError("Predicate %s is always true." % name)
            if redundancyInPredicateMasks(predicate_masks):
                raise ValueError("Redundancy in predicates.")
            if redundancyInPredicateMasks(predicates_plus_scale):
                raise ValueError("Some combination of predicates is"
                        " equivalent to the overall rate parameter.")
        else:
            if redundancyInPredicateMasks(predicate_masks):
                raise ValueError("Redundancy in predicates.")
            if redundancyInPredicateMasks(predicates_plus_scale):
                warnings.warn("do_scaling=True would be more efficient than"
                        " these overly general predicates")

        self.predicate_masks = predicate_masks
        self.parameter_order = []
        self.predicate_indices = []
        self.symmetric = _isSymmetrical(self._instantaneous_mask)
        for pred in predicate_order:
            mask = predicate_masks[pred]
            if not _isSymmetrical(mask):
                self.symmetric = False
            indices = numpy.nonzero(mask)
            assert numpy.alltrue(mask[indices] == 1)
            self.parameter_order.append(pred)
            self.predicate_indices.append(indices)
        if not self.symmetric:
            warnings.warn('Model not reversible')
        (self.scale_masks, scale_order) = self._adaptPredicates(scales or [])
        self.checkParamsExist()
            
    def calcExchangeabilityMatrix(self, mprobs, *params):
        assert len(params) == len(self.predicate_indices), self.parameter_order
        R = self._instantaneous_mask_f.copy()
        for (indices, par) in zip(self.predicate_indices, params):
            R[indices] *= par
        return R
    
    def asciiArt(self, delim='', delim2='|', max_width=70):
        """An ASCII-art table representing the model.  'delim' delimits
        parameter names, 'delim2' delimits motifs"""
        # Should be implemented with table module instead.
        
        pars = self.getMatrixParams()
        par_names = self.getParamList()
        longest = max([len(name) for name in (par_names+[' '])])
        if delim:
            all_names_len = _maxWidthIfTruncated(pars, delim, 100)
            min_names_len = _maxWidthIfTruncated(pars, delim, 1)
        else:
            all_names_len = sum([len(name) for name in par_names])
            min_names_len = len(par_names)
        
        # Find a width-per-motif that is as big as can be without being too big
        w = min_names_len
        while (w+1) * len(self.alphabet) < max_width and w < all_names_len:
            w += 1
        
        # If not enough width truncate parameter names
        if w < all_names_len:
            each = w / len(par_names)
            if delim:
                while _maxWidthIfTruncated(pars, delim, each+1) <= w:
                    each += 1
                w = _maxWidthIfTruncated(pars, delim, each)
            else:
                w = each * len(par_names)
        else:
            each = longest
        
        rows = []
        # Only show header if there is enough width for the motifs
        if self.alphabet.getMotifLen() <= w:
            header = [str(motif).center(w) for motif in self.alphabet]
            header = [' ' * self.alphabet.getMotifLen() + ' '] + header + ['']
            header = delim2.join(header)
            rows.append(header)
            rows.append(''.join([['-',delim2][c == delim2] for c in header]))
        
        # pars in sub-cols, should also offer pars in sub-rows?
        for (motif, row2) in zip(self.alphabet, pars):
            row = []
            for par_list in row2:
                elt = []
                for par in par_names:
                    if par not in par_list:
                        par = ''
                    par = par[:each]
                    if not delim:
                        par = par.ljust(each)
                    if par:
                        elt.append(par)
                elt = delim.join(elt).ljust(w)
                row.append(elt)
            rows.append(delim2.join(([motif+' '] + row + [''])))
        return '\n'.join(rows)
    
    def getMatrixParams(self):
        """Return the parameter assignment matrix."""
        dim = len(self.alphabet)
        Pars = numpy.zeros([dim, dim], object)
        for x, y in [(x, y) for x in range(dim) for y in range(dim)]:
            Pars[x][y] = []  # a limitation of numpy.  [x,y] = [] fails!
            if not self._instantaneous_mask[x, y]:
                continue
            for par in self.predicate_masks:
                if self.predicate_masks[par][x, y]:
                    Pars[x, y].append(par)
            # sort the matrix entry to facilitate scaling calculations
            Pars[x, y].sort()
        return Pars
    
    def getParamList(self):
        """Return a list of parameter names."""
        return self.predicate_masks.keys()
    
    def isInstantaneous(self, x, y):
        return self._isInstantaneous(x, y)
    
    def getSubstitutionRateValueFromQ(self, Q, motif_probs, pred):
        pred_mask = self._adaptPredicates([pred])[0].values()[0]
        pred_row_totals = numpy.sum(pred_mask * Q, axis=1)
        inst_row_totals = numpy.sum(self._instantaneous_mask * Q, axis=1)
        r = sum(pred_row_totals * motif_probs)
        t = sum(inst_row_totals * motif_probs)
        pred_size = numpy.sum(pred_mask.flat)
        inst_size = sum(self._instantaneous_mask.flat)
        return (r / pred_size) / ((t-r) / (inst_size-pred_size))
    
    def getScaledLengthsFromQ(self, Q, motif_probs, length):
        lengths = {}
        for rule in self.scale_masks:
            lengths[rule] = length * self.getScaleFromQs(
                    [Q], [1.0], motif_probs, rule)
        return lengths
    
    def getScaleFromQs(self, Qs, bin_probs, motif_probss, rule):
        rule = self.getPredicateMask(rule)
        weighted_scale = 0.0
        bin_probs = numpy.asarray(bin_probs)
        for (Q, bin_prob, motif_probs) in zip(Qs, bin_probs, motif_probss):
            row_totals = numpy.sum(rule * Q, axis=1)
            motif_probs = numpy.asarray(motif_probs)
            word_probs = self.calcWordProbs(motif_probs)
            scale = sum(row_totals * word_probs)
            weighted_scale += bin_prob * scale
        return weighted_scale
    
    def getPredefinedPredicates(self):
        # overridden in subclasses
        return {'indel': predicate.parse('-/?')}
    
    def getPredefinedPredicate(self, name):
        # Called by predicate parsing code
        if self._canned_predicates is None:
            self._canned_predicates = self.getPredefinedPredicates()
        return self._canned_predicates[name].interpret(self)
        
    def _adaptPredicates(self, rules):
        # dict or list of callables, predicate objects or predicate strings
        if isinstance(rules, dict):
            rules = rules.items()
        else:
            rules = [(None, rule) for rule in rules]
        predicate_masks = {}
        order = []
        for (key, pred) in rules:
            (label, mask) = self.adaptPredicate(pred, key)
            if label in predicate_masks:
                raise KeyError('Duplicate predicate name "%s"' % label)
            predicate_masks[label] = mask
            order.append(label)
        return predicate_masks, order
    
    def adaptPredicate(self, pred, label=None):
        if isinstance(pred, str):
            pred = predicate.parse(pred)
        elif callable(pred):
            pred = predicate.UserPredicate(pred)
        pred_func = pred.makeModelPredicate(self)
        label = label or repr(pred)
        mask = predicate2matrix(
            self.getAlphabet(), pred_func, mask=self._instantaneous_mask)
        return (label, mask)
    
    def getPredicateMask(self, pred):
        if pred in self.scale_masks:
            mask = self.scale_masks[pred]
        elif pred in self.predicate_masks:
            mask = self.predicate_masks[pred]
        else:
            (label, mask) = self.adaptPredicate(pred)
        return mask


class _Nucleotide(SubstitutionModel):
    def getPredefinedPredicates(self):
        return {
            'transition' : predicate.parse('R/R') | predicate.parse('Y/Y'),
            'transversion' : predicate.parse('R/Y'),
            'indel': predicate.parse('-/?'),
            }
    

class Nucleotide(_Nucleotide):
    """A nucleotide substitution model."""
    def __init__(self, **kw):
        SubstitutionModel.__init__(self, moltype.DNA.Alphabet, **kw)
    

class Dinucleotide(_Nucleotide):
    """A nucleotide substitution model."""
    def __init__(self, **kw):
        SubstitutionModel.__init__(self, moltype.DNA.Alphabet, motif_length=2, **kw)
    

class Protein(SubstitutionModel):
    """Base protein substitution model."""
    def __init__(self, with_selenocysteine=False, **kw):
        alph = moltype.PROTEIN.Alphabet
        if not with_selenocysteine:
            alph = alph.getSubset('U', excluded=True)
        SubstitutionModel.__init__(self, alph, **kw)
    

def EmpiricalProteinMatrix(matrix, motif_probs=None, optimise_motif_probs=False,
        recode_gaps=True, do_scaling=True, **kw):
    alph = moltype.PROTEIN.Alphabet.getSubset('U', excluded=True)
    return Empirical(alph, rate_matrix=matrix, motif_probs=motif_probs,
            model_gaps=False, recode_gaps=recode_gaps, do_scaling=do_scaling,
            optimise_motif_probs=optimise_motif_probs, **kw)
    

class Codon(_Nucleotide):
    """Core substitution model for codons"""
    long_indels_are_instantaneous = True
    
    def __init__(self, alphabet=None, gc=None, **kw):
        if gc is not None:
            alphabet = moltype.CodonAlphabet(gc = gc)
        alphabet = alphabet or moltype.STANDARD_CODON
        SubstitutionModel.__init__(self, alphabet, **kw)
    
    def _isInstantaneous(self, x, y):
        if x == self.gapmotif or y == self.gapmotif:
            return x != y
        else:
            ndiffs = sum([X!=Y for (X,Y) in zip(x,y)])
            return ndiffs == 1
    
    def getPredefinedPredicates(self):
        gc = self.getAlphabet().getGeneticCode()
        def silent(x, y):
            return x != '---' and y != '---' and gc[x] == gc[y]
        def replacement(x, y):
            return x != '---' and y != '---' and gc[x] != gc[y]
        
        preds = _Nucleotide.getPredefinedPredicates(self)
        preds.update({
            'indel' : predicate.parse('???/---'),
            'silent' : predicate.UserPredicate(silent),
            'replacement' : predicate.UserPredicate(replacement),
            })
        return preds

