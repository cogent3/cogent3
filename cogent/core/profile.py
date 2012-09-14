#!/usr/bin/env python
"""Provides Profile and ProfileError object and CharMeaningProfile

Owner: Sandra Smit (Sandra Smit)
"""
from __future__ import division
#SUPPORT2425
#from __future__ import with_statement
from string import maketrans, translate

from numpy import array, sum, transpose, reshape, ones, zeros,\
    take, float64, ravel, nonzero, log, put, concatenate, argmax, cumsum,\
    sort, argsort, searchsorted, logical_and, asarray, uint8, add, subtract,\
    multiply, divide, newaxis, alltrue, max, all, isfinite
#from numpy.oldnumeric import sum
from numpy.random import random
from cogent.util.array import euclidean_distance, row_degeneracy,\
    column_degeneracy, row_uncertainty, column_uncertainty, safe_log
from cogent.format.table import formattedCells
##SUPPORT2425
import numpy #from cogent.util.unit_test import numpy_err

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Sandra Smit", "Gavin Huttley", "Rob Knight"
                    "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

class ProfileError(Exception):
    """Error raised for exceptions occuring in the Profile object"""
    pass


class Profile(object):
    """Profile class
    """
    
    def __init__(self, Data, Alphabet, CharOrder=None):
        """Initializes a new Profile object.

        Data: numpy 2D array with the Profile data in it.
        Specifically, each row of the array corresponds to a position in the
        Alignment. Each column of the array corresponds to a character in the
        Alphabet.
        
        Alphabet: an Alphabet object or anything that can act as a list of
            characters
        CharOrder: optional list of characters to which the columns
            in the Data correspond.
        """
        self.Data = Data
        self.Alphabet = Alphabet
        if CharOrder is None:
            self.CharOrder = list(self.Alphabet)
        else:
            self.CharOrder = CharOrder
        #the translation table is needed for making consensus sequences,
        #but will fail if the alphabet isn't made of chars (in which case,
        #we'll just skip the translation table, and certain downstream
        #operations may fail).
        try:
            self._translation_table = self._make_translation_table()
        except:
            pass
    
    def __str__(self):
        """Returns string representation of self.Data"""
        return str(self.Data)
    
    def _make_translation_table(self):
        """Makes a translation tables between the CharOrder and indices
        """
        indices = ''.join(map(chr, range(len(self.CharOrder))))
        chars = ''.join(map(str,self.CharOrder))
        return maketrans(chars, indices)

    def hasValidData(self, err=1e-16):
        """Returns True if all rows in self.Data add up to one
        
        err -- float, maximum deviation from 1 allowed, default is 1e-16

        Rounding errors might occur, so a small deviation from 1 is allowed.
            The default tolerance is 1e-16.
        """
        obs_sums = sums = sum(self.Data,1)
        lower_bound = ones(len(self.Data)) - err
        upper_bound = ones(len(self.Data)) + err
        tfs = sum(self.Data,1) == ones(len(self.Data))
        if (lower_bound <= obs_sums).all() and (obs_sums <= upper_bound).all():
            return True
        return False

    def hasValidAttributes(self):
        """Checks Alphabet, CharOrder, and size of self.Data"""
        if not reduce(logical_and, [c in self.Alphabet\
            for c in self.CharOrder]):
                return False
        elif self.Data.shape[1] != len(self.CharOrder):
            return False
        return True

    def isValid(self):
        """Check whether everything in the Profile is valid"""
        vd = self.hasValidData()
        va = self.hasValidAttributes()
        return vd and va

    def dataAt(self, pos, character=None):
        """Return data for a certain position (row!) and character (column)

        pos -- int, position (row) in the profile
        character -- str, character from the CharacterOrder

        If character is None, all data for the position is returned.
        """
        if not 0 <= pos < len(self.Data):
            raise ProfileError(\
            "Position %s is not present in the profile"%(pos))
        if character is None:
            return self.Data[pos,:]
        else:
            if character not in self.CharOrder:
                raise ProfileError(\
                "Character %s is not present in the profile's CharacterOrder"\
                    %(character))
            return self.Data[pos, self.CharOrder.index(character)]
 
    def copy(self):
        """Returns a copy of the Profile object

        WARNING: the data, alphabet and the character order are the same
        object in the original and the copy. This means you can rebind the
        attributes, but modifying them will change them in both the original
        and the copy.
        """ 
        return self.__class__(self.Data, self.Alphabet, self.CharOrder)
    
    def normalizePositions(self):
        """Normalizes the data by position (the rows!) to one

        It does not make sense to normalize anything with negative
        numbers in there. However, the method does NOT check for that, 
        because it would slow down the calculations too much. It will work, 
        but you might get very unexpected results.

        The method will raise an error when one or more rows add up to
        one. It checks explicitly for that to avoid OverflowErrors, 
        ZeroDivisionErrors, and infinities in the results. 

        WARNING: this method works in place with respect to the Profile
        object, not with respect to the Data attribute. Normalization
        rebinds self.Data to a new array.
        """
        row_sums = sum(self.Data,1)
        if (row_sums == 0).any():
            zero_indices = nonzero(row_sums==0)[0].tolist()
            raise ProfileError,\
            "Can't normalize profile, rows at indices %s add up to zero"\
            %(zero_indices)
        else:
            self.Data = self.Data/row_sums[:,newaxis]

    def normalizeSequences(self):
        """Normalized the data by sequences (the columns) to one
        
        It does not make sense to normalize anything with negative
        numbers in there. However, the method does NOT check for that, 
        because it would slow down the calculations too much. It will work, 
        but you might get very unexpected results.

        The method will raise an error when one or more columns add up to
        one. It checks explicitly for that to avoid OverflowErrors, 
        ZeroDivisionErrors, and infinities in the results. 

        WARNING: this method works in place with respect to the Profile
        object, not with respect to the Data attribute. Normalization
        rebinds self.Data to a new array.
        """
        col_sums = sum(self.Data, axis=0)
        if (col_sums == 0).any():
            zero_indices = nonzero(col_sums==0)[0].tolist()
            raise ProfileError,\
            "Can't normalize profile, columns at indices %s add up to zero"\
            %(zero_indices)
        else:
            self.Data = self.Data/col_sums
            
    def prettyPrint(self, include_header=False, transpose_data=False,\
        column_limit=None, col_sep='\t'):
        """Returns a string method of the data and character order.

        include_header: include charcter order or not
        transpose_data: data as is (rows are positions) or transposed
            (rows are characters) to line it up with an alignment
        column_limit = int, maximum number of columns displayed
        col_sep = string, column separator
        """
        h = self.CharOrder
        d = self.Data
        if column_limit is None:
            max_col_idx = d.shape[1]
        else:
            max_col_idx = column_limit
        if include_header and not transpose_data:
            r = [h]+d.tolist()
        elif include_header and transpose_data:
            r =[[x] + y for x,y in zip(h,transpose(d).tolist())]
        elif transpose_data:
            r = transpose(d).tolist()
        else:
            r = d.tolist()
        # resize the result based on the column limit
        if column_limit is not None:
            r = [row[:column_limit] for row in r]
        # nicely format the table content, discard the header (already included) 
        if r:
            new_header, formatted_res = formattedCells(r)
        else:
            formatted_res = r
        return '\n'.join([col_sep.join(map(str,i)) for i in formatted_res])
        
    def reduce(self,other,op=add,normalize_input=True,normalize_output=True):
        """Reduces two profiles with some operator and returns a new Profile
        
        other: Profile object
        op: operator (e.g. add, subtract, multiply, divide)
        normalize_input: whether the input profiles will be normalized
            before collapsing. The default is True.
        normalize_output: whether the output profile will be normalized.
            The default is True

        This function is intented for use on normalized profiles. For
        safety it'll try to normalize the data before collapsing them.
        If you do not normalize your data and set normalize_input to 
        False, you might get unexpected results. 
        
        It does check whether self.Data and other.Data have the same shape
        It does not check whether self and other have the same 
        CharOrder. The resulting Profile gets the alphabet and
        char order from self.

        """
        if self.Data.shape != other.Data.shape:
            raise ProfileError,\
                "Cannot collapse profiles of different size: %s, %s"\
                %(self.Data.shape,other.Data.shape)
        if normalize_input:
            self.normalizePositions()
            other.normalizePositions()
        
        try:
            ##SUPPORT2425
            ori_err = numpy.geterr()
            numpy.seterr(divide='raise')
            try: new_data = op(self.Data, other.Data)
            finally: numpy.seterr(**ori_err)
            #with numpy_err(divide='raise'):
                #new_data = op(self.Data, other.Data)
        except (OverflowError, ZeroDivisionError, FloatingPointError):
            raise ProfileError, "Can't do operation on input profiles"
        result = Profile(new_data, self.Alphabet, self.CharOrder)
        
        if normalize_output:
            result.normalizePositions()
        return result

    def __add__(self,other):
        """Binary + operator: adds two profiles element-wise.
        
        Input and output are NOT normalized.
        """
        return self.reduce(other, op=add, normalize_input=False,\
            normalize_output=False)
    
    def __sub__(self,other):
        """Binary - operator: subtracts two profiles element-wise
        
        Input and output are NOT normalized.
        """
        return self.reduce(other, op=subtract, normalize_input=False,\
            normalize_output=False)

    def __mul__(self,other):
        """* operator: multiplies two profiles element-wise
        
        Input and output are NOT normalized.
        """
        return self.reduce(other, op=multiply, normalize_input=False,\
            normalize_output=False)

    def __div__(self,other):
        """/ operator for old-style division: divides 2 profiles element-wise.

        Used when __future__.divsion not imported

        Input and output are NOT normalized.
        """
        return self.reduce(other, op=divide, normalize_input=False,\
            normalize_output=False)

    def __truediv__(self,other):
        """/ operator for new-style division: divides 2 profiles element-wise.

        Used when __future__.division is in action.
        
        Input and output are NOT normalized.
        """
        return self.reduce(other, op=divide, normalize_input=False,\
            normalize_output=False)

    def distance(self, other, method=euclidean_distance):
        """Returns the distance between two profiles

        other: Profile object
        method: function used to calculated the distance between two
        arrays.

        WARNING: In principle works only on profiles of the same size. 
        However, when one of the two profiles is 1D (which shouldn't 
        happen) and can be aligned with the other profile the distance 
        is still calculated and may give unexpected results.
        """
        try:
            return method(self.Data, other.Data)
        except ValueError: #frames not aligned 
            raise ProfileError,\
            "Profiles have different size (and are not aligned): %s %s"\
            %(self.Data.shape,other.Data.shape)
    
    def toOddsMatrix(self, symbol_freqs=None):
        """Returns the OddsMatrix of a profile as a new Profile.

        symbol_freqs: per character array of background frequencies
        e.g. [.25,.25,.25,.25] for equal frequencies for each of the 
        four bases.

        If no symbol frequencies are provided, all symbols will get equal 
        freqs. The length of symbol freqs should match the number of 
        columns in the profile! If symbol freqs contains a zero entry,
        a ProfileError is raised. This is done to prevent either a
        ZeroDivisionError (raised when zero is an int) or 'inf' in the 
        resulting matrix (which happens when zero is a float).
        """
        pl = self.Data.shape[1] #profile length
        #if symbol_freqs is None, create an array with equal frequencies
        if symbol_freqs is None:
            symbol_freqs = ones(pl)/pl
        else:
            symbol_freqs = array(symbol_freqs)
            
        #raise error when symbol_freqs has wrong length
        if len(symbol_freqs) != pl:
            raise ProfileError,\
            "Length of symbol freqs should be %s, but is %s"\
            %(pl,len(symbol_freqs))
        
        #raise error when symbol freqs contains zero (to prevent 
        #ZeroDivisionError or 'inf' in the resulting matrix)
        if sum(symbol_freqs != 0, 0) != len(symbol_freqs):
            raise ProfileError,\
            "Symbol frequency is not allowed to be zero: %s"\
            %(symbol_freqs)

        #calculate the OddsMatrix
        log_odds = self.Data/symbol_freqs
        return Profile(log_odds, self.Alphabet, self.CharOrder)

    def toLogOddsMatrix(self, symbol_freqs=None):
        """Returns the LogOddsMatrix of a profile as a new Profile/

        symbol_freqs: per character array of background frequencies
        e.g. [.25,.25,.25,.25] for equal frequencies for each of the 
        four bases.

        See toOddsMatrix for more information.
        """
        odds = self.toOddsMatrix(symbol_freqs)
        log_odds = safe_log(odds.Data)
        return Profile(log_odds, self.Alphabet, self.CharOrder)

    def _score_indices(self, seq_indices, offset=0):
        """Returns score of the profile for each slice of the seq_indices
        
        seq_indices: translation of sequence into indices that match the 
        characters in the CharOrder of the profile
        offset: where to start the matching procedure

        This function doesn't do any input validation. That is done in 'score'
        See method 'score' for more information.
        """
        data = self.Data
        pl = len(data) #profile length (number of positions)
        sl = len(seq_indices)

        r = range(pl) #fixed range
        result = []
        for starting_pos in range(offset, len(seq_indices)-pl+1):
            slice = seq_indices[starting_pos:starting_pos+pl]
            result.append(sum(array([data[i] for i in zip(r,slice)]), axis=0))
        return array(result)
    
    def _score_profile(self, profile, offset=0):
        """Returns score of the profile against the input_profile.

        profile: Profile of a sequence or alignment that has to be scored
        offset: where to start the matching procedure

        This function doesn't do any input validation. That is done in 'score'
        See method 'score' for more information.
        """
        data = self.Data
        self_l = len(data) #profile length
        other_l = len(profile.Data) #other profile length
        result = []
        for start in range(offset,other_l-self_l+1):
            stop = start + self_l
            slice = profile.Data[start:stop,:]
            result.append(sum(self.Data*slice))
        return array(result)

    def score(self, input_data, offset=0):
        """Returns a score of the profile against input_data (Profile or Seq).

        seq: Profile or Sequence object (or string)
        offset: starting index for searching in seq/profile

        Returns the score of the profile against all possible subsequences/
        subprofiles of the input_data. 
    
        This method determines how well a profile fits at different places
        in the sequence. This is very useful when the profile is a motif and 
        you want to find the position in the sequence that best matches the
        profile/motif. 

        Sequence Example:
        =================
            T   C   A   G
        0   .2  .4  .4  0
        1   .1  0   .9  0
        2   .1  .2  .3  .4

        Sequence: TCAAGT

        pos 0: TCA -> 0.5
        pos 1: CAA -> 1.6
        pos 2: AAG -> 1.7
        pos 3: AGT -> 0.5

        So the subsequence starting at index 2 in the sequence has the 
        best match with the motif

        Profile Example:
        ================
        Profile: same as above
        Profile to score:
            T   C   A   G
        0   1   0   0   0
        1   0   1   0   0   
        2   0   0   .5  .5
        3   0   0   0   1   
        4   .25 .25 .25 .25
       
        pos 0: rows 0,1,2 -> 0.55
        pos 1: rows 1,2,3 -> 1.25
        pos 2: rows 2,3,4 -> 0.45
        """

        #set up some local variables
        data = self.Data
        pl = len(data) #profile length
        is_profile = False

        #raise error if profile is empty
        if not data.any():
            raise ProfileError,"Can't score an empty profile"
        
        #figure out what the input_data type is
        if isinstance(input_data,Profile):
            is_profile = True
            to_score_length = len(input_data.Data)
            #raise error if CharOrders don't match
            if self.CharOrder != input_data.CharOrder:
                raise ProfileError, "Profiles must have same character order"
        else: #assumes it get a sequence
            to_score_length = len(input_data)
       
        #Profile should fit at least once in the sequence/profile_to_score
        if to_score_length < pl:
            raise ProfileError,\
            "Sequence or Profile to score should be at least %s "%(pl)+\
            "characters long, but is %s."%(to_score_length)
        #offset should be valid
        if not offset <= (to_score_length - pl):
            raise ProfileError, "Offset must be <= %s, but is %s"\
            %((to_score_length-pl), offset)

        #call the apropriate scoring function
        if is_profile:
            return self._score_profile(input_data, offset)
        else:
            #translate seq to indices
            if hasattr(self, '_translation_table'):
                seq_indices = array(map(ord,translate(str(input_data),\
                    self._translation_table)))
            else:   #need to figure out where each item is in the charorder
                idx = self.CharOrder.index
                seq_indices = array(map(idx, input_data))
            #raise error if some sequence characters are not in the CharOrder
            if (seq_indices > len(self.CharOrder)).any():
                raise ProfileError,\
                "Sequence contains characters that are not in the "+\
                "CharOrder"
            #now the profile is scored against the list of indices   
            return self._score_indices(seq_indices,offset)
     
    def rowUncertainty(self):
        """Returns the uncertainty (Shannon's entropy) for each row in profile
        
        Entropy is returned in BITS (not in NATS).
        """
        if not self.Data.any():
            return array([])
        try:
            return row_uncertainty(self.Data)
        except ValueError:
            raise ProfileError,\
            "Profile has to be two dimensional to calculate rowUncertainty"
            
    def columnUncertainty(self):
        """Returns uncertainty (Shannon's entropy) for each column in profile

        Uncertainty is returned in BITS (not in NATS).
        """
        if not self.Data.any():
            return array([])
        try:
            return column_uncertainty(self.Data)
        except ValueError:
            raise ProfileError,\
            "Profile has to be two dimensional to calculate columnUncertainty"

    def rowDegeneracy(self, cutoff=0.5):
        """Returns how many chars are needed to cover the cutoff value.

        cutoff: value that should be covered in each row

        For example:
        pos 0: [.1,.2,.3,.4] char order=TCAG. 
        If cutoff=0.75 -> degeneracy = 3 (degenearate char for CAG)
        If cutoff=0.25 -> degeneracy = 1 (G alone covers this cutoff)
        If cutoff=0.5  -> degeneracy = 2 (degenerate char for AG)

        If the cutoff value is not reached in the row, the returned value
        will be clipped to the length of the character order (=the number
        of columns in the Profile).
        """
        try:
            return row_degeneracy(self.Data,cutoff)
        except ValueError:
            raise ProfileError,\
            "Profile has to be two dimensional to calculate rowDegeneracy"

    def columnDegeneracy(self, cutoff=0.5):
        """Returns how many chars are neede to cover the cutoff value

        See rowDegeneracy for more information.
        """
        try:
            return column_degeneracy(self.Data,cutoff)
        except ValueError:
            raise ProfileError,\
            "Profile has to be two dimensional to calculate columnDegeneracy"

    def rowMax(self):
        """Returns ara containing most frequent element in each row of the profile."""
        return max(self.Data, 1)

    def toConsensus(self, cutoff=None, fully_degenerate=False,\
        include_all=False):
        """Returns the consensus sequence from a profile.

        cutoff: cutoff value, determines how much should be covered in a
        position (row) of the profile. Example: pos 0 [.2,.1,.3,.4]
        (CharOrder: TCAG). To cover .65 (=cutoff) we need two characters:
        A and G, which results in the degenerate character R.
        
        fully_degenerate: determines whether the fully degenerate character
        is returned at a position. For the example above an 'N' would
        be returned.
       
        inlcude_all: all possibilities are included in the degenerate 
        character. Example: row = UCAG = [.1,.3,.3,.3] cutoff = .4, 
        consensus = 'V' (even though only 2 chars would be enough to 
        reach the cutoff value).

        The Alphabet of the Profile should implement degenerateFromSequence.
        
        Note that cutoff has priority over fully_degenerate. In other words,
        if you specify a cutoff value and set fully_degenerate to true, 
        the calculation will be done with the cutoff value. If nothing 
        gets passed in, the maximum argument is chosen. In the first example
        above G will be returned.
        """
        #set up some local variables
        co = array(self.CharOrder, 'c')
        alpha = self.Alphabet
        data = self.Data

        #determine the action. Cutoff takes priority over fully_degenerate
        if cutoff:
            result = []
            degen = self.rowDegeneracy(cutoff)
            sorted = argsort(data)
            if include_all:
                #if include_all include all possiblilities in the degen char 
                for row_idx, (num_to_keep, row) in enumerate(zip(degen,sorted)):
                    to_take = [item for item in row[-num_to_keep:]\
                    if item in nonzero(data[row_idx])[0]] +\
                    [item for item in nonzero(data[row_idx] ==\
                        data[row_idx,row[-num_to_keep]])[0] if item in\
                        nonzero(data[row_idx])[0]]
                    result.append(alpha.degenerateFromSequence(\
                    map(str,take(co, to_take, axis=0))))
            else:
                for row_idx, (num_to_keep, row) in enumerate(zip(degen,sorted)):
                    result.append(alpha.degenerateFromSequence(\
                        map(str,take(co, [item for item in row[-num_to_keep:]\
                        if item in nonzero(data[row_idx])[0]]))))
                                    
        elif not fully_degenerate: 
            result = take(co, argmax(self.Data, axis=-1), axis=0)
        else:
            result = []
            for row in self.Data:
                result.append(alpha.degenerateFromSequence(\
                map(str,take(co, nonzero(row)[0], axis=0))))
        return ''.join(map(str,result))


    def randomIndices(self, force_accumulate=False, random_f = random):
        """Returns random indices matching current probability matrix.

        Stores cumulative sum (sort of) of probability matrix in 
        self._accumulated; Use force_accumulate to reset if you change 
        the matrix in place (which you shouldn't do anyway).

        The returned indices correspond to the characters in the
        CharOrder of the Profile.
        """
        if force_accumulate or not hasattr(self, '_accumulated'):
            self._accumulated = cumsum(self.Data, 1)
        choices = random_f(len(self.Data))
        return array([searchsorted(v, c) for v, c in\
            zip(self._accumulated, choices)])

    def randomSequence(self, force_accumulate=False, random_f = random):
        """Returns random sequence matching current probability matrix.

        Stores cumulative sum (sort of) of probability matrix in 
        self._accumulated; Use force_accumulate to reset if you change 
        the matrix in place (which you shouldn't do anyway).
        """
        co = self.CharOrder
        random_indices = self.randomIndices(force_accumulate,random_f)
        return ''.join(map(str,take(co,random_indices)))



        
def CharMeaningProfile(alphabet, char_order=None, split_degenerates=False):
    """Returns a Profile with the meaning of each character in the alphabet
    
    alphabet: Alphabet object (should have 'Degenerates'if split_degenerates
    is set to True)
    char_order: string indicating the order of the characters in the profile
    split_degenerates: whether the meaning of degenerate symbols in the 
    alphabet should be split up among the characters in the char order,
    or ignored.

    The returned profile has 255 rows (one for each ascii character) and
    one column for each character in the character order. The profile
    specifies the meaning of each character in the alphabet. Chars in the
    character order will count as a full character by themselves, degenerate
    characters might split their 'value' over several other charcters in the
    character order. 

    Splitting up degenerates: only degenerate characters of which the full 
    set of symbols it maps onto are in the character order are split up, 
    others are ignored. E.g. in the DnaAlphabet, if the char order is 
    TACG, ? (which maps to TCAG-) wouldn't be split up, 'R' (which maps 
    to 'AG') would.
    
    Any degenerate characters IN the character order will NOT be split up.
    It doesn't make sense to split up a character that is in the char order
    because it would create an empty column in the profile, so it might
    as well be left out alltogether.
    
    Example 1:
    Alphabet = DnaAlphabet
    Character order = "TCAG"
    Split degenerates = False
    
    All the nonzero rows in the resulting profile are:
    65: [0,0,1,0] (A)
    67: [0,1,0,0] (C)
    71: [0,0,0,1] (G)
    84: [1,0,0,0] (T)
    All other rows will be [0,0,0,0].

    Example 2:
    Alphabet = DnaAlphabet
    Character order = "AGN"
    Split degenerates = True
    
    All the nonzero rows in the resulting profile are:
    65: [1,0,0] (A)
    71: [0,1,0] (G)
    78: [0,0,1] (N)
    82: [.5,.5,0] (R)
    All other rows will be [0,0,0].
    
    Errors are raised when the character order is empty or when there's a
    character in the character order that is not in the alphabet.
    """
    if not char_order: #both testing for None and for empty string
        char_order = list(alphabet)
    
    char_order = array(char_order, 'c')
    lc = len(char_order) #length char_order

    #initialize the profile. 255 rows (one for each ascii char), one column
    #for each character in the character order
    result = zeros([255,lc],float64)

    if split_degenerates:
        degen = alphabet.Degenerates
        for degen_char in degen:
            #if all characters that the degenerate character maps onto are
            #in the character order, split its value up according to the
            #alphabet
            curr_degens = degen[degen_char]
            if all(map(char_order.__contains__, curr_degens)):
                contains = map(curr_degens.__contains__, char_order)
                result[ord(degen_char)] = \
                        array(contains, float)/len(curr_degens)
    #for each character in the character order, make an entry of ones and 
    #zeros, matching the character order
    for c in char_order:
        c = str(c)
        if c not in alphabet:
            raise ValueError, "Found character in the character order "+\
            "that is not in the specified alphabet: %s"%(c) 
        result[ord(c)] = array(c*lc, 'c') == char_order           
    return Profile(Data=result,Alphabet=alphabet,CharOrder=char_order)

