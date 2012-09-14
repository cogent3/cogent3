#!/usr/bin/env python
"""markov.py: various types of random and non-random generators.

Currently provides:

MarkovGenerator: reads in k-word frequencies from a text, and generates random text based on those frequencies. Also calculates the average entropy of the (k+1)th symbol (for sufficiently large k, converges on the entropy of the text).

NOTE:  The text must be a list of strings (e.g. lines of text).  If 
a single string is passed into the constructor it should be 
put into a list (i.e. ['your_string']) or it will result in errors
when calculating kword frequencies. 
"""
from __future__ import division
from operator import mul
from random import choice, shuffle, randrange
from cogent.maths.stats.util import UnsafeFreqs as Freqs
from cogent.util.array import cartesian_product
from cogent.maths.stats.test import G_fit
from copy import copy,deepcopy
from numpy import ones, zeros, ravel, array, rank, put, argsort, searchsorted,\
                  take
from numpy.random import random

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Jesse Zaneveld", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

class MarkovGenerator(object):
    """Holds k-word probabilities read from file, and can generate text."""
    def __init__(self, text=None, order=1, linebreaks=False, \
        calc_entropy=False, freqs=None, overlapping=True, \
        array_pseudocounts=0, delete_bad_suffixes=True):
        """Sets text and generates k-word frequencies."""
        self.Text = text
        self.Linebreaks = linebreaks
        self.Order = order
        self.Frequencies = freqs or {}
        self.RawCounts = {}
        self.FrequencyArray=None
        self.CountArray=None
        self.ExcludedCounts=None
        self.ArrayPseudocounts=array_pseudocounts
        self._calc_entropy = calc_entropy
        self.Entropy = None
        self.Prior = None
        self.Overlapping=overlapping
        if self.Text:
            self.calcFrequencies(delete_bad_suffixes)

    def calcFrequencies(self, delete_bad_suffixes=True):
        """For order k, gets the (k-1)-word frequencies plus what follows."""
        #reset text if possible -- but it might just be a string, so don't
        #complain if the reset fails.
        overlapping=self.Overlapping
        try:
            self.Text.reset()
        except AttributeError:
            try:
                self.Text.seek(0)
            except AttributeError:
                pass
        k = self.Order
        if k < 1:   #must be 0 or '-1': just need to count single bases
            self._first_order_frequency_calculation()
        else:   #need to figure out what comes after the first k bases
            all_freqs = {}
            for line in self.Text:
                if not self.Linebreaks:
                    line = line.strip()
                #skip the line if it's blank
                if (not line):
                    continue
                #otherwise, make a frequency distribution of symbols
                end = len(line) - k
                if overlapping:
                    rang=xrange(end)
                else:
                    rang=xrange(0,end,(k+1))
                for i in rang:
                    word, next = line[i:i+k], line[i+k]
                    curr = all_freqs.get(word, None)
                    if curr is None:
                        curr = Freqs({next:1})
                        all_freqs[word] = curr
                    else:
                        curr += next
            if self._calc_entropy:
                self.Entropy = self._entropy(all_freqs)
            self.Frequencies = all_freqs
            if delete_bad_suffixes:
                self.deleteBadSuffixes()
            self.RawCounts=deepcopy(all_freqs)
            #preserve non-normalized freqs
            for dist in self.Frequencies.values():
                dist.normalize()

    def wordToUniqueKey\
            (self,word, conversion_dict={'a':0,'c':1,'t':2,'g':3}):
        #since conversion_dict values are used as array indices later, 
        #values of conversion dict should range from 0 to (n-1), 
        #where n=number of characters in your alphabet
        uniqueKey=0
        alpha_len = len(conversion_dict)
        for i in range(0,len(word)):
            uniqueKey += (conversion_dict[word[i]]*alpha_len**i)
        return uniqueKey

    def makeCountArray(self):
        """Generates a 1 column array with indices equal to the keys for each 
        kword + character and raw counts of the occurances of that key as values.
        
        This allows counts for many k+1 long strings to be 
        found simultaneously with evaluateArrayProbability 
        (which also normalizes the raw counts to frequencies)"""
        #print "makeCountArray: before replaceDegen self.Rawcounts=",\
        #        self.RawCounts
        self.replaceDegenerateBases()
        #print "makeCountArray:self.RawCounts=",self.RawCounts #debugging
        counts=self.RawCounts
        self.CountArray=zeros((4**(self.Order+1)),'f')
        #Order 0 --> 4 spots ('a','c','t','g') 1 --> 16 etc
        #TODO: may generate problems if Order = -1
            
        if self.Order==0:
            #print "makeCountArray:counts=",counts #debugging
            for key in counts['']:
                #print "attempting to put",float(counts[''][key]),"into index",\
                #        self.wordToUniqueKey(key),\
                #        "of array CountArray=",self.CountArray #debugging
                put(self.CountArray,self.wordToUniqueKey(key),\
                        float(counts[''][key]))
                self.CountArray[self.wordToUniqueKey(key)]=counts[''][key]
                #print "placement successful!" #debugging
        else:
            for kword in counts.keys():
                for key in counts[kword]:
                    index=self.wordToUniqueKey(kword+key)
                    #debugging
                    #print "attempting to put",counts[kword][key],"at index",\
                    #        index,"of self.CountArray, which =",self.CountArray
                    
                    put(self.CountArray,index,counts[kword][key])
                    #print "placement sucessful!" #debugging
        #print "makeCountArray:raveling self" #debugging 
        if self.ArrayPseudocounts:
            self.CountArray = self.CountArray + float(self.ArrayPseudocounts)
            # adds to each count, giving unobserved keys frequency 
            #pseudocounts/n 
            #n= number of observed counts (rather than 0 frequency)
            # When the number of pseudocounts added is one,
            # this is 'Laplace's rule'
            #(See 'Biological Sequence Analysis',Durbin et. al, p.115)
        self.CountArray=ravel(self.CountArray)
        #print "makeCountArray:final CountArray=",self.CountArray

    def updateFrequencyArray(self):
        """updates the frequency array by re-normalizing CountArray"""
        self.FrequencyArray=deepcopy(self.CountArray) #preserve raw counts
        total_counts=sum(self.FrequencyArray)
        self.FrequencyArray=self.FrequencyArray/total_counts
       
    def replaceDegenerateBases(self,normal_bases=['a','t','c','g']):
        """remove all characters from self.Text that aren't
        a,t,c or g and replace them with random characters
        (when degenerate characters are rare, this is useful
        because it avoids assigning all kwords with those 
        characters artificially low conditional probabilities)"""
        
        def normalize_character(base,bases=normal_bases):
            if base not in bases:
                base=choice(bases)
            return base
        
        text=self.Text
        for i in range(len(text)):
            text[i]=\
                    ''.join(map(normalize_character,text[i].lower()))
        self.Text=text

    def deleteBadSuffixes(self):
        """Deletes all suffixes that can't lead to prefixes.

        For example, with word size 3, if acg is present but cg* is not
        present, acg is not allowed.

        Need to repeat until no more suffixes are deleted.
        """
        f = self.Frequencies
        #loop until we make a pass where we don't delete anything
        deleted = True
        while deleted:
            deleted = False
            for k, v in f.items():
                suffix = k[1:]
                for last_char in v.keys():
                    #if we can't make suffix + last_char, can't select that char
                    if suffix + last_char not in f:
                        del v[last_char]
                        deleted=True
                if not v:   #if we deleted the last item, delete prefix
                     del f[k]
                     deleted = True

    def _entropy(self, frequencies):
        """Calcuates average entropy of the (k+1)th character for k-words."""
        sum_ = 0.
        sum_entropy = 0.
        count = 0.
        for i in frequencies.values():
            curr_entropy = i.Uncertainty
            curr_sum = sum(i.values())
            sum_ += curr_sum
            sum_entropy += curr_sum * curr_entropy
            count += 1
        return sum_entropy/sum_
        

    def _first_order_frequency_calculation(self):
        """Handles single-character calculations, which are independent.

        Specifically, don't need to take into account any other characters, and
        can just feed the whole thing into a single Freqs.
        """
        freqs = Freqs('')
        for line in self.Text:
            freqs += line
        #get rid of line breaks if necessary
        if not self.Linebreaks:
            for badkey in ['\r', '\n']:
                try:
                    del freqs[badkey]
                except KeyError:
                    pass    #don't care if there weren't any
        #if order is negative, equalize the frequencies
        if self.Order < 0:
            for key in freqs:
                freqs[key] = 1
        self.RawCounts= {'':deepcopy(freqs)}
        freqs.normalize()
        self.Frequencies = {'':freqs}

    def next(self, length=1, burn=0):
        """Generates random text of specified length with current freqs.
        burn specifies the number of iterations to throw away while the chain
        converges.
        """
        if self.Order < 1:
            return self._next_for_uncorrelated_model(length)
        freqs = self.Frequencies    #cache reference since it's frequently used
        #just pick one of the items at random, since calculating the weighted
        #frequencies is not possible without storing lots of extra info
        keys = freqs.keys()
        curr = choice(keys)
        result = []
        for i in range(burn +length):
            next = freqs[curr].choice(random())
            if i >= burn:
                result.append(next)
            curr = curr[1:] + next
        return ''.join(result)

    def _next_for_uncorrelated_model(self, length):
        """Special case for characters that don't depend on previous text."""
        return ''.join(self.Frequencies[''].randomSequence(length))

    def evaluateProbability(self,seq):
        """Evaluates the probability of generating a 
        user-specified sequence given the model."""
        conditional_prob=1
        order=self.Order
        for i in range(0,(len(seq)-(order)),1):
            k=seq[i:i+order+1]
            try:
                conditional_prob *= self.Frequencies[k[:-1]][k[-1]]
            except KeyError:
                #if key not in Frequencies 0 < Freq < 1/n
                #To be conservative in the exclusion of models, use 1/n 
                if conditional_prob:
                    conditional_prob *= 1.0/(float(len(self.Text)))
                else:
                    conditional_prob = 1.0/(float(len(self.Text)))
        return conditional_prob
        
    def evaluateWordProbability(self,word):
        k=word[:self.Order+1]
        try:
            conditional_prob= self.Frequencies[k[:-1]][k[-1]]
        except KeyError:
            conditional_prob = 1.0/(float(len(self.Text)))
        return conditional_prob
    
    def evaluateArrayProbability(self,id_array):
        #takes an array of unique integer keys
        #corresponding to (k+1) long strings
        #[can be generated by self.wordToUniqueKey()]
        #Outputs probability
        if self.FrequencyArray is None:
            if self.CountArray is None:
                self.makeCountArray()
            self.updateFrequencyArray()
        freqs=take(self.FrequencyArray,id_array)
        prob=reduce(mul,freqs)
        return float(prob)
 
    def evaluateInitiationFrequency(self,kword,\
            allowed_bases=['a','t','c','g']):
        # takes a unique key corresponding to a k long word
        # calculates the initiation frequency for that kword
        # which is equal to its relative frequency

        #TODO: add case where order is < 1
        if len(kword) != (self.Order):
            raise KwordError  #kword must be equal to markov model order
        if self.CountArray is None:
            self.makeCountArray()
        unique_keys=[]
        #add term for each possible letter
        for base in allowed_bases:
            unique_keys.append(self.wordToUniqueKey(kword+base))
        id_array=ravel(array(unique_keys))
        counts=take(self.CountArray,id_array)
        total_kword_counts=sum(counts)
        total_counts=sum(self.CountArray)
        prob=float(total_kword_counts)/float(total_counts)
        return prob

    def excludeContribution(self,excluded_texts):
        #"""Excludes the contribution of a set of texts 
        #from the markov model.  This can be useful, for example,
        #to prevent self-contribution of the data in a gene to 
        #the model under which that gene is evaluated.
        #
        #A Markov Model is made from the strings, converted to a CountArray,
        #and then that Count array (stored as ExcludedCounts) is subtracted
        #from the current CountArray, and FrequencyArray is updated.
        #
        #The data excluded with this function can be restored with
        #restoreContribution
        #
        #Only one list of texts can be excluded at any time.  If a list of 
        #texts is already excluded when excludeContribution is called, that 
        #data will be restored before the new data is excluded"""
        
        #print "excludeContribution:excluded_texts=",excluded_texts #debugging
        if self.CountArray is None:
            #print ".excludeContribution:missing countArray" #debugging
            self.makeCountArray()
        if self.ExcludedCounts:
            self.restoreContribution()
        #generate mm using same parameters as current model
        exclusion_model=MarkovGenerator(excluded_texts,order=self.Order,\
                overlapping=self.Overlapping)
        exclusion_model.makeCountArray()
        self.ExcludedCounts=\
                (exclusion_model.CountArray)
        self.CountArray = self.CountArray-(self.ExcludedCounts)
        self.updateFrequencyArray()
        
    def restoreContribution(self):
        """Restores data excluded using excludeContribution, and 
        renormalizes FrequencyArray"""
        if self.ExcludedCounts:
            self.CountArray += (self.ExcludedCounts)
            self.ExcludedCounts=None
            self.updateFrequencyArray()

        
    
def count_kwords(source, k, delimiter=''):
    """Makes dict of {word:count} for specified k."""
    result = {}
    #reset to beginning if possible
    if hasattr(source, 'seek'):
        source.seek(0)
    elif hasattr(source, 'reset'):
        source.reset()
    if isinstance(source, str):
        if delimiter:
            source = source.split(delimiter)
        else:
            source = [source]
    for s in source:
        for i in range(len(s) - k + 1):
            curr = s[i:i+k]
            if curr in result:
                result[curr] += 1
            else:
                result[curr] = 1
    return result

def extract_prefix(kwords):
    """Converts dict of {w:count} to {w[:-1]:{w[-1]:count}}"""
    result = {}
    for w, count in kwords.items():
        prefix = w[:-1]
        suffix = w[-1]
        if prefix not in result:
            result[prefix] = {}
        curr = result[prefix]
        if suffix not in curr:
            curr[suffix] = {}
        curr[suffix] = count
    return result

def _get_expected_counts(kwords, kminus1):
    """Gets expected counts from counts of 2 successive kword lengths.."""
    result = []
    total = sum(kminus1.values())   #shortest
    prefixes = extract_prefix(kminus1)
    for k in kwords:
        result.append(kminus1[k[:-1]] * prefixes[k[1:-1]]/kminus1[k[:-1]])

def _pair_product(p, i, j):
    """Return product of counts of i and j from data."""
    try:
        return sum(p[i].values())*sum(p[j].values())
    except KeyError:
        return 0

def markov_order(word_counts, k, alpha):
    """Estimates Markov order of a source, using G test for fit.

    Uses following procedure:

    A source depends on the previous k letters more than the previous (k-1)
    letters iff Pr(a|w_k) != Pr(a|w_{k-1}) for all words of length k. If we 
    know Pr(a|w) for all symbols a and words of length k and k-1, we would
    expect count(a|w_i) to equal count(a|w_i[1:]) * count(w)/count(w[1:]).
    We can compare these expected frequencies to observed frequencies using
    the G test.
    
    max_length: maximum correlation length to try
    """
    if k == 0:  #special case: test for unequal freqs
        obs = word_counts.values()
        total = sum(obs)    #will remain defined through loop
        exp = [total/len(word_counts)] * len(word_counts)
    elif k == 1: #special case: test for pair freqs
        prefix_counts = extract_prefix(word_counts)
        total = sum(word_counts.values())
        words = word_counts.keys()
        
        exp = [_pair_product(prefix_counts, w[0], w[1])/total for w in words]
        obs = word_counts.values()
    else:   # k >= 3: need to do general Markov chain
        #expect count(a_i.w.b_i) to be Pr(b_i|w)*count(a_i.w)
        cwb = {}    #count of word.b
        cw = {}     #count of word
        caw = {}    #count of a.word
        #build up counts of prefix, word, and suffix
        for word, count in word_counts.items():
            aw, w, wb = word[:-1], word[1:-1], word[1:]
            if not wb in cwb:
                cwb[wb] = 0
            cwb[wb] += count
            if not aw in caw:
                caw[aw] = 0
            caw[aw] += count
            if not w in cw:
                cw[w] = 0
            cw[w] += count

        obs = word_counts.values()
        exp = [cwb[w[1:]]/(cw[w[1:-1]])*caw[w[:-1]] for w in \
            word_counts.keys()]
    return G_fit(obs, exp)

def random_source(a, k, random_f=random):
    """Makes a random Markov source on alphabet a with memory k.
    
    Specifically, for all words k, pr(i|k) = rand().
    """
    result = dict.fromkeys(map(''.join, cartesian_product([a]*k)))
    for k in result:
        result[k] = Freqs(dict(zip(a, random_f(len(a)))))
    return result

def markov_order_tests(a, max_order=5, text_len=10000, verbose=False):
    """Tests of the Markov order inferrer using Markov chains of diff. orders.
    """
    result = []
    max_estimated_order = max_order + 2
    for real_order in range(max_order):
        print "Actual Markov order:", real_order
        s = random_source(a, real_order)
        m = MarkovGenerator(order=real_order, freqs=s)
        text = m.next(text_len)
        for word_length in range(1, max_estimated_order+1):
            words = count_kwords(text, word_length)
            g, prob = markov_order(words, word_length-1, a)
            if verbose:
                print "Inferred order: %s G=%s P=%s" % (word_length-1, g, prob)
            result.append([word_length-1, g, prob])
    return result
            
         
        

if __name__ == '__main__':
    """Makes text of specified # chars from training file.
    
    Note: these were tested from the command line and confirmed working by RK 
    on 8/6/07.
    """
    from sys import argv, exit
    if len(argv) == 2 and argv[1] == 'm':
        markov_order_tests('tcag', verbose=True)
    elif len(argv) == 3 and argv[1] == 'm':
        infilename = argv[2]
        max_estimated_order = 12
        text = open(infilename).read().split('\n')
        for word_length in range(1, max_estimated_order):
            words = count_kwords(text, word_length)
            g, prob = markov_order(words, word_length-1, 'ATGC')
            print "Inferred order: %s G=%s P=%s" % (word_length-1, g, prob)
            
    else:
        try:
            length = int(argv[1])
            max_order = int(argv[2])
            text = open(argv[3], 'U')
        except:
           print "Usage: python markov.py num_chars order training_file"
           print "...or python markov.py m training_file to check order"
           exit()
        for order in range(max_order + 1):
            m = MarkovGenerator(text, order, calc_entropy=True)
            print order,':', 'Entropy=', m.Entropy, m.next(length=length)
        
