#!/usr/bin/env python
# knots.py

"""Contains code related to RNA (secondary) structure and pseudoknots.

Specifically, this module contains several methods to remove pseudoknots
from RNA structures. Pseudoknot removal is discussed in the following paper:

S. Smit, K. Rother, J. Heringa, and R. Knight
Manuscript in preparation.

If you use this code in your work, please cite this publication (in addition
to the PyCogent publication). If you need to cite the paper before submission,
please contact the author of this module.

Six functions are provided (see documentation of each function for
a more detailed description):

opt_all -- optimization approach that calculates all nested structures
    that optimize some value (e.g. keep the maximum number of base pairs)
conflict_elimination -- Removes pseudoknots from a structure by eliminating
    conflicting base pairs one by one. Two functions to determine which 
    paired region should be removed next are provided: max_conlficts
    and min_gain.
inc_order -- creates a nested structure by adding non-conflicting paired
    regions one by one to the solution; paired regions are processed
    from 5' to 3' start point or from 3' to 5' end point.
inc_length -- creates a nested structure by adding non-conflicting paired
    regions one at the time, starting with the longest region working towards
    the shortest region.
inc_range -- generates a nested structure by adding non-conflicting paired
    retions one at the time, starting with short-range interactions working
    towards long-range interactions.

These six functions represent the core objective of this module.
Two convenience functions supporting the opt_all function are added:
opt_single_random, and opt_single_property
There is also a modified version of the original Nussinov-Jacobson
algorithm present which is restricted to the given list of base pairs:
nussinov_restricted

In addition, the following supporting objects and functions are present:

PairedRegion -- object that represents a paired region in an RNA structure.
    A paired region is an uninterrupted stretch of base pairs with positions
    [(i,j),(i+1, j-1),(i+2,j-2), ...].
PairedRegions -- object (basically a list) that stores a collection of 
    PairedRegion objects. This is an alternative way of representing an
    RNA structure, where basically stretches of base pairs are condensed
    into PairedRegion objects.
PairedRegionFromPairs -- Factory function to create a PairedRegion object
    from a Pairs object (cogent.struct.rna2d)
PairedRegionsFromPairs -- Factory function to create a PairedRegions object
    from a Pairs object (cogent.struct.rna2d)
ConflictMatrix -- object to store a matrix of conflicts between different 
    paired regions. Row and Column indices correspond to PairedRegion IDs, 
    values in the matrix are True (if the regions conflict), and False (if
    the regions don't conflict).

Other smaller helper functions are:
contains_true, empty_matrix, pick_multi_best, dp_matrix_multi,
matrix_solutions, and add_back_non_conflicting
See their docstrings for detailed documentation.

NOTE:
None of the methods provided can handle overlapping base pairs (i.e. base
A pairs with B, and A pairs with C), because in that case it is unclear
whether the first or second pair should be kept (also, you could get 
different results based on the order of the base pairs). In general,
one removes pseudoknots in order to represent the structure in dot-bracket
format. If the list contains conflict, this is impossible anyway, so removing
the pseudoknots would not help. It is the responsibility of the user to 
remove overlapping base pairs before trying to obtain a nested structure.

"""

from __future__ import division
from random import choice
from numpy import sum, average, zeros
from cogent.struct.rna2d import Pairs
from cogent.util.dict2d import Dict2D

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Sandra Smit, Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"


class PairedRegion(object):
    """Store an uninterrupted list of base pairs by start, end, and length
    
    A paired region (a.k.a. ladder or (helical) region) is a stretch of
        perfectly nested base pairs with positions:
        [(m,n), (m+1,n-1), (m+2,n-2),...]
    This object is very similar to the Stem object in cogent.struct.rna2d.
        In addition to the start, end, and length it stores the actual 
        base pairs, and it stores a region ID. It has many more methods
        than the Stem object.
    This object performs no error checking. You can specify and End 
        before a Start point, or a Start and End which are closer 
        to each other than 2 times the Length.
    """

    def __init__(self, Start, End, Length, Id=None):
        """Initialize a new PairedRegion object

        Start -- int, specifying the starting index of the paired region
            in the sequence. This is the 5' side of the 5' halfregion.
        End -- int, specifying the end index of the paired region
            in the sequence. This is the 3' side of the 3' halfregion.
        Length -- int, specifying the length of the paired region, i.e. the
            number of base pairs in the region. 
        Id -- string or int, unique identifier for this PairedRegion.

        During initialization a Pairs object is created. The first pair is
        always (Start,End), additional pairs add up from the Start point and
        down from the End point, so (Start+1,End-1), (Start+2, End-2) etc.
        """
        self.Start = Start
        self.End = End
        if Length < 1:
            raise ValueError(\
                "PairedRegion should contain at least one base pair")
        self.Length = Length
        self.Id = Id
        self.Pairs = Pairs()
        for i in range(self.Length):
            self.Pairs.append((self.Start+i, self.End-i)) 

    def __str__(self):
        """Return string representation of PairedRegion, list of Pairs
        """
        return str(self.Pairs)

    def __len__(self):
        """Return Length of the PairedRegion, i.e. the number of base pairs
        """
        return self.Length

    def __eq__(self, other):
        """Compares Pairs and IDs
       
        other -- PairedRegion object

        If IDs are not set (both None), Pairs is the only criterion. 
        """
        if self.Pairs == other.Pairs and self.Id == other.Id:
            return True
        return False

    def __ne__(self, other):
        """Return True if two PairedRegion objects differ

        other -- PairedRegion object
        """
        return not self == other

    def upstream(self):
        """Return list of upstream positions in self from 5' to 3'
        """
        return [i for (i,j) in self.Pairs]

    def downstream(self):
        """Return list of downstream positions in self from 5' to 3'
        """
        result = [j for (i,j) in self.Pairs]
        result.reverse()
        return result

    def paired(self):
        """Return sorted list of paired positions in this region
        """
        result = self.upstream() + self.downstream()
        result.sort()
        return result

    def range(self):
        """Return the range of this region
        
        The range of the region is the number of bases between the
        highest upstream and the lowest downstream position.
        (i.e. the number of unpaired bases in the hairpin if this
        were the only paired region in the structure)
        For example: the range of region with start=3, end=10, and len=2
            would be 4.
        Performs no error checking. If region overlaps with itself,
            a negative number will be returned...
        """
        return min(self.downstream()) - max(self.upstream()) - 1

    def overlapping(self, other):
        """Returns True if two regions overlap

        other -- PairedRegion object

        Two regions overlap if there is at least one base which
            is a member of both regions. (Definition from Studnicka 1978)
        Identical regions are overlapping.
        """
        ref_pos = dict.fromkeys(self.paired())
        for pos in other.paired():
            if pos in ref_pos:
                return True
        return False

    def conflicting(self, other):
        """Return True if the regions are conflicting, False if they are nested

        other -- PairedRegion object

        Two paired regions are conflicting if they are organized in a 
            knotted fashion. This means both other.Start and other.End 
            have to be between self.Start and self.End or both shouldn't be.
            See for example Studnicka 1978 for defintion, or any other
            paper with a general pseudoknot definiton in there.
        Overlapping regions cause an error, because you can't determine
            whether they are conflicting or not. However, two identical
            regions are defined as NOT conflicting, even though they are
            overlapping, and thus False is returned. For non-identical,
            but overlapping regions an error will be raised.
        """
        if self == other: # equal blocks
            return False
        # if not equal, but overlapping, raise error
        if self.overlapping(other):
            raise ValueError("Can only handle non-overlapping regions")
        
        if (other.Start > self.End and other.End > self.End) or\
            (self.Start > other.End and self.End > other.End):
            return False
        if (self.Start < other.Start < self.End and\
            self.Start < other.End < self.End) or\
            (other.Start < self.Start < other.End and\
            other.Start < self.End < other.End):
            return False
        return True

    def score(self, scoring_function):
        """Sets self.Score to value of scoring function applies to self

        scoring_function -- function that can be applied to a PairedRegion
            object and returns a numerical value

        Note: this method has no return value, it sets a property of the 
            instance.
        """
        self.Score = scoring_function(self)

def PairedRegionFromPairs(pairs, Id=None):
    """Return new PairedRegion object from Pairs object

    pairs -- Pairs object or list of tuples with up and downstream positions.
    Id -- string or int, unique identifier of this region.

    This is a factory function to create a PairedRegion object
        from a Pairs object. It assumes the pairs are fully nested
        with positions [(m,n), (m+1,n-1), (m+2,n-2),...].
    The function doesn't validate the input pairs, so if the assumtion
        does not hold, the Pairs in the resulting PairedRegion might
        differ from the input pairs.
    It makes the pairs directed and sorts them, then it extracts the 
        Start, End and Length and initialized a new PairedRegion with 
        those parameters.
    """
    if not pairs:
        raise ValueError("PairedRegion should contain at least one pair")
    # preprocess
    raw_pairs = Pairs(pairs)
    if raw_pairs.hasConflicts():
        raise ValueError("Cannot handle pairs with conflicts")
    p = raw_pairs.directed()
    p.sort()
    # initialize variables
    start = p[0][0]
    end = p[0][1]
    length = len(p)
    return PairedRegion(start, end, length, Id=Id) 


class PairedRegions(list):
    """Stores a list of PairedRegion objects
    
    A PairedRegions object is a condensed way of looking at an RNA structure,
        where continuous stretches of base pairs are collapsed into 
        PairedRegion objects. See the documentation on PairedRegion for more
        details.
    """
    
    def __init__(self, regions=None):
        """Initialize new PairedRegions object

        regions -- list of PairedRegion objects

        The object is meant to store a list of PairedRegion objects,
            however it is very light-weight and does not perform any
            validation on the input.
        """
        if regions is None:
            regions = []
        self[:] = list(regions)

    def __str__(self):
        """Return string representation of PairedRegions object

        Each PairedRegion is presented as Id:Start,End,Length;
            A PairedRegions object is presented as '(' + space-delimited
            list of PairedRegion objects + ')'
            For example: (A:2,10,2; B:12,20,3;)
        """
        result = []
        for i in self:
            result.append("%s:%s,%s,%s;"%(i.Id,i.Start,i.End,i.Length))
        return '('+' '.join(result)+')'
    
    def __eq__(self, other):
        """Return True if two PairedRegions objects are equal

        other -- PairedRegions object

        Two regions are equal if they have the same length and contain
            the same PairedRegion objects.
        """
        if len(self) != len(other):
            return False
        for i in self:
            if i not in other:
                return False
        return True

    def __ne__(self, other):
        """Return True if two PairedRegions objects are different

        other -- PairedRegions object
        """
        return not self == other

    def byId(self):
        """Return dict of {ID: PairedRegion}

        This function only works if the IDs for each region in self
            are unique. If multiple regions with the same ID are found, 
            an error is raised.
        """
        result = {}
        for pr in self: # for every paired region
            if pr.Id in result:
                raise ValueError("Duplicate key found")
            result[pr.Id] = pr
        return result

    def numberOfRegions(self):
        """Return the number of PairedRegion objects in the list
        """
        return len(self)

    def totalLength(self):
        """Return the cumulative length of all PairedRegion objects in self
        
        The totalLength is the total number of base pairs in this
            PairedRegions object. So, it adds the number of pairs in each
            PairedRegion in the list.
        """
        if self:
            return sum(map(len, self))
        else:
            return 0

    def totalScore(self):
        """Return sum of Score values of each PairedRegion in self

        This method simply adds all the Score attributes (!= None)
            for each PairedRegion in this PairedRegions object.
        """
        score = 0
        for pr in self:
            try:
                score += pr.Score
            except AttributeError:
                raise ValueError("Score not set for %s"%(str(self)))
            except TypeError:
                raise ValueError("Score should be numerical, but is %s"\
                    %(pr.Score))
        return score

    def toPairs(self):
        """Return Pairs object containing all the pairs in each PairedRegion

        This method does not validate the pairs. It simply adds all the 
            pairs in each PairedRegion to the result. Pairs might occur twice
            in the result. The resulting Pairs object is sorted.
        """
        result = Pairs()
        for pr in self:
            result.extend(pr.Pairs)
        result.sort()
        return result

    def byStartEnd(self):
        """Return dict of {(pr.Start, pr.End): pr}

        Keys in the dictionary are tuples of start and end positions,
            the values are the PairedRegion objects themselves. If a 
            Start/End combination is already in the dictionary, an error is 
            raised.
        """
        result = {}
        for pr in self:
            se = (pr.Start, pr.End)
            if se in result:
                raise ValueError("Duplicate key found: %s"%(str(se)))
            result[se] = pr
        return result
        #return dict([((pr.Start, pr.End), pr) for pr in self])

    def lowestStart(self):
        """Return lowest begin value of any PairedRegion in the list
        
        It lists all the Start values for all the PairedRegion objects
            in the list and returns the lowest value. If there are no
            regions in self, None is returned.
        """
        start_values = [pr.Start for pr in self]
        if not start_values:
            return None
        else:
            return min(start_values)

    def highestEnd(self):
        """Return highest end value of any PairedRegion in the list

        It lists all the End values for all the PairedRegion objects
            in the list and returns the highest value. If there are no
            regions in self, None is returned.
        """
        end_values = [pr.End for pr in self]
        if not end_values:
            return None
        else:
            return max(end_values)

    def sortedIds(self):
        """Return sorted list of region IDs
        """
        all_ids = [pr.Id for pr in self]
        all_ids.sort()
        return all_ids

    def upstream(self):
        """Return sorted list of upstream positions
        """
        result = []
        for pr in self:
            result.extend(pr.upstream())
        result.sort()
        return result

    def downstream(self):
        """Return sorted list of downstream positions
        """
        result = []
        for pr in self:
            result.extend(pr.downstream())
        result.sort()
        return result

    def pairedPos(self):
        """Return sorted list of all paired positions
        """
        result = self.upstream() + self.downstream()
        result.sort()
        return result

    def boundaries(self):
        """Return sorted list of all start and end points
        """
        result = []
        for pr in self:
            result.append(pr.Start)
            result.append(pr.End)
        result.sort()
        return result

    def enumeratedBoundaries(self):
        """Return dict of {boundary_index: boundary}

        Return value is dictionary created from tuples from the enumeration
        of all boundaries.
        """
        return dict(enumerate(self.boundaries()))

    def invertedEnumeratedBoundaries(self):
        """Return dict of {boundary_value: boundary_idx}

        Boundary values are all the start and end points of the paired
            regions in self. Should be unique, otherwise an error is raised.
        Boundary indices are the indices assigned to each start and end point
            during an enumeration of the sorted list. 
        Overall the result is the inverted dictionary of the result of the
            enumeratedBoundaries method.
        """
        eb = self.enumeratedBoundaries()
        result = {}
        for boundary_idx, boundary_value in eb.items():
            if boundary_value in result:
                raise ValueError(\
                    "Boundary value %s is not unique"%(boundary_value))
            result[boundary_value] = boundary_idx
        #result = dict([(v,k) for k,v in eb.items()])
        return result

    def merge(self, other):
        """Merge two PairedRegions objects together

        other -- PairedRegions object

        Duplicate PairedRegion objects are stored only once.
        This methods used the PairedRegion IDs to check for duplications.
        """
        result = PairedRegions()
        seen = {}
        for pr in self+other:
            if pr.Id not in seen:
                result.append(pr)
                seen[pr.Id] = True
        return result

    def conflicting(self, cm=None):
        """Return PairedRegions obj containing regions involved in a conflict

        cm -- ConflictMatrix for this PairedRegions object.

        This method only works if the PairedRegion objects have unique IDs, 
            because a conflict matrix is constructed. This behavior can be 
            changed...
        See PairedRegion.conflicting() for a definition of conflicting
            paired regions.
        """
        if cm is None:
            cm = ConflictMatrix(self)
        id_to_pr = self.byId()

        result = PairedRegions()
        for pr_id in cm.conflicting():
            result.append(id_to_pr[pr_id])
        return result

    def nonConflicting(self, cm=None):
        """Return new PairedRegions object containing non-conflicing regions

        cm -- ConflictMatrix for this PairedRegions object.

        Two PairedRegion objects do not conflict when they are organized
            in a nested fashion.
        """
        if cm is None:
            cm = ConflictMatrix(self)
        id_to_pr = self.byId()

        result = PairedRegions()
        for pr_id in cm.nonConflicting():
            result.append(id_to_pr[pr_id])
        return result

    def conflictCliques(self, cm=None):
        """Return list of PairedRegions objects w/ mutually conflicting regions

        cm -- ConflictMatrix for this PairedRegions object.

        Mutually conflicting regions form a knot-component as defined
            in Rodland 2006.
        Return value is a list of PairedRegions objects. Each PairedRegions
            object contains mutually conflicting regions (knot-components)
        E.g. region A conflicts with B and B conflicts with A and C, and D
            conflicts with E, and F doesn't conflict with any othe region,
            one group would be A, B and C, the other would be D and E.
            F would not be returned in any group, since it isn't conflicting.
        """
        if cm is None:
            cm = ConflictMatrix(self)
        id_to_pr = self.byId()
        cliques = cm.conflictCliques()
        result = []
        for cl in cliques:
            pr = PairedRegions()
            for i in cl:
                pr.append(id_to_pr[i])
            result.append(pr)
        return result

 
def PairedRegionsFromPairs(pairs):
    """Return PairedRegions object from Pairs

    pairs -- Pairs object, no conflicts allowed

    Result is a list of stretches of perfectly nested base pairs, which means
    [[(m,n), (m+1,n-1), (m+2,n-2),...],[(i,j),(i+1,j-1),...]]
    Base pairs are made directed and sorted before the stretches are
        picked out, so the result will be in order.
    IDs of the regions are set as indices (enumeration of all regions). The
        PairedRegion that starts closest to the 5' end will get ID 0, the next
        ID 1, etc.
    """
    p = Pairs(pairs)
    if not p:
        return PairedRegions()
    if p.hasConflicts():
        raise ValueError("Cannot handle base pair conflicts")
    clean_pairs = p.directed()
    clean_pairs.sort()
    regions = []
    curr_region = []
    pr_id = -1 # paired region ID
    for pair in clean_pairs:
        if not curr_region:
            curr_region.append(pair)
        else:
            x,y = curr_region[-1]
            if pair == (x+1, y-1):
                curr_region.append(pair)
            else:
                pr_id += 1
                regions.append(PairedRegionFromPairs(curr_region, Id=pr_id))
                curr_region = [pair]
    if curr_region: # last block
        pr_id += 1
        regions.append(PairedRegionFromPairs(curr_region, Id=pr_id))
    return PairedRegions(regions)

class ConflictMatrix(object):
    """Stores conflict matrix

    A conflict matrix is a matrix that indicates which PairedRegion objects
        are conflicting. Row and column IDs correspond to Region IDs. If
        two regions are conflicting True is stored, otherwise False is stored.
    """

    def __init__(self, data):
        """Initialize new ConflictMatrix object
        
        data -- either a PairedRegions object or a Pairs object,
            or anything that can be made into a Pairs object
            (e.g. a list of tuples)

        This method sets the Matrix attribute to a Dict2D object containing
            conflict information on the PairedRegions. 
        The input data is either a PairedRegions object or it is made into one.
            A ValueError will be raised when the pairs or regions are
            overlapping. Input data that can't be converted to Pairs
            will lead to downstream errors. Pairs doesn't perform
            any validation.
        Row and column IDs are the Identifiers of the PairedRegion objects.
        The RowOrder and ColumnOrder of the Dict2D are the sorted region IDs.
        """
        if isinstance(data, PairedRegions):
            id_to_pr = data.byId()
        elif isinstance(data, Pairs):
            id_to_pr = PairedRegionsFromPairs(data).byId()
        else: # try to convert to Pairs 
            try:
                d = Pairs(data)
                id_to_pr = PairedRegionsFromPairs(d).byId()
            except:
                raise ValueError("Can't convert data to Pairs")

        # handle the rows and columns in order and set RowOrder and ColOrder
        ro = id_to_pr.keys()
        ro.sort()
        co = id_to_pr.keys()
        co.sort()
        conf = {} # dict of conflicts between blocks
        for id1, bl in id_to_pr.items():
            for id2, bl2 in id_to_pr.items():
                if id2 < id1: # minimize number of calculations
                    continue
                if id1 not in conf:
                    conf[id1] = {}
                if id2 not in conf:
                    conf[id2] = {}
                if id1 == id2:
                    conf[id1][id2] = False
                    conf[id2][id1] = False
                    continue
                is_conflicting = bl.conflicting(bl2)
                conf[id1][id2] = is_conflicting
                conf[id2][id1] = is_conflicting
        self.Matrix = Dict2D(conf, RowOrder=ro, ColOrder=co) # create Dict2D

    def conflictsOf(self, pr_id):
        """Return list of region IDs for regions that conflict with pr_id
        
        pr_id -- row ID in the matrix (ID of paired region)

        Input is ID of a particular region, return value are the IDs of
            all regions that conflict with the given region.
        """
        return [k for k,v in self.Matrix[pr_id].items() if v is True]

    def conflicting(self):
        """Return list of region IDs for conflicting regions
        """
        result = []
        cm = self.Matrix
        for pr_id in cm.RowOrder:
            if contains_true(cm[pr_id].values()):
                result.append(pr_id)
        return result

    def nonConflicting(self):
        """Return list of region IDs for non-conflicting regions
        """
        result = []
        cm = self.Matrix
        for pr_id in cm.RowOrder:
            if not contains_true(cm[pr_id].values()):
                result.append(pr_id)
        return result

    def conflictCliques(self):
        """Return list of lists with IDs of mutually conflicting regions

        See documentation on PairedRegions.conflictCliques for more details.
        """
        cm = self.Matrix

        cliques = []
        seen = {}
        
        for pr_id in cm.RowOrder:
            if pr_id in seen:
                continue
            todo = set([pr_id])
            done = set()
            while todo != done:
                collection = [] # collection of conflicts
                for i in todo:
                    if i in done: # no need to do them twice
                        continue
                    conf = [] # conflicting regions
                    for k,v in cm[i].items():
                        if v is True:
                            conf.append(k)
                    collection.extend(conf) # add conflict to collection
                    done.add(i) # register that i is done
                todo.update(collection) # update todo
            if len(done) > 1:
                cliques.append(list(done))
            for i in done:
                seen[i] = True
        return cliques

# =============================================================================
# SCORING FUNCTIONS FOR DYNAMIC PROGRAMMING APPROACH
# =============================================================================

def num_bps(paired_region):
    """Return number of base pairs (=Length) of paired_region
    
    paired_region -- PairedRegion object
    """
    return paired_region.Length

def hydrogen_bonds(seq):
    """Return function to score a PairedRegion by its hydrogen bonds

    seq -- Sequence object or string

    This method counts the number of hydrogen bonds in Watson-Crick and 
        Wobble base pairs. GC pairs score 3, AU and GU pairs score 2.
    """
    HB_SCORE = {('G','C'): 3, ('C','G'): 3,\
                ('A','U'): 2, ('U','A'): 2,\
                ('G','U'): 2, ('U','G'): 2}

    def apply_to(paired_region):
        """Return score of paired_region by its hydrogen bonds
        
        paired_region -- PairedRegion object

        Scores each base pair in the region by giving each GC base pair 3 
            points and each AU or GU base pair 2 points. Other base pairs
            are ignored and don't add anything to the overall score.
        """
        score = 0
        for up,down in paired_region.Pairs:
            seq_pair = (seq[up],seq[down])
            try:
                score += HB_SCORE[seq_pair]
            except KeyError:
                continue
        return score
    return apply_to

# =============================================================================
# HELPER FUNCTIONS FOR DYNAMIC PROGRAMMING APPROACH
# =============================================================================

def contains_true(i):
    """Return True if input contains True
    
    i -- any object that implements __contains__

    Returns True if True is in the input. Both True and 1 count as True.
    Helper function for ConflictMatrix object.
    """
    try:
        if True in i:
            return True
        return False
    except TypeError: # when i is a string
        return False

def empty_matrix(size):
    """Return square matrix as list of lists of specified size.

    size -- int, number of rows and columns in the matrix.

    This function is a helper function of opt_all.
    Each cell is filled with [PairedRegions()]. This is the initialization
        value needed for a dynamic programming matrix that keeps track
        of all optimal solutions. A solution is a single PairedRegions object.
    """
    if size < 1:
        raise ValueError("The size of the matrix should be at least one")
    result = []
    for i in range(size):
        result.append([])
        for j in range(size):
            result[i].append([PairedRegions()])
    return result

def pick_multi_best(candidates, goal='max'):
    """Return list of unique solutions with a maximum/minimum score.

    candidates -- list of PairedRegions objects

    This function returns a list of all PairedRegions objects that
        have an optimal score (maximum or minimum depending on the goal).
        If the list of candidates is empty, a list containing an
        empty PairedRegions object is returned.
    This function is a helper function of dp_matrix_multi.

    NOTE: PairedRegion IDs must be set. They are checked to avoid
        including unsaturated solutions. Maybe implementation should be
        changed, such that (Start, End, Length) tuples are used as IDs?!
    """
    if not candidates:
        return [PairedRegions()]
    result = []
    best_score = None
    seen = {}
    # Candidates have to be processed in order of length
    can_len = [(c.totalLength(), c) for c in candidates]
    can_len.sort()
    can_len.reverse()
    for l, c in can_len:
        c_ids = tuple(c.sortedIds())
        c_ids_set = set(c_ids)
        if not c or c_ids in seen:
            continue
        this_score = c.totalScore()
        if best_score is None:
            best_score = this_score
            result = [c]
            seen[c_ids] = True
        elif this_score == best_score:
            is_sub = False
            for seen_id in seen:
                if len(c_ids_set) == len(c_ids_set & set(seen_id)):
                    is_sub = True
                    break
            if is_sub:
                seen[c_ids] = True
            else:
                result.append(c)
                seen[c_ids] = True
        elif goal == 'max' and this_score < best_score:
            continue
        elif goal == 'min' and this_score > best_score:
            continue
        else:
            result = [c]
            seen[c_ids] = True
            best_score = this_score

    if not result:
        return [PairedRegions()]
    return result


def dp_matrix_multi(paired_regions, goal='max', scoring_function=num_bps): 
    """Return dynamic programming matrix with top-right half filled

    paired_regions -- PairedRegions object
    goal -- str, 'max' or 'min', if the goal is 'max' the routine
        returns the solutions maximizing the score, if the goal
        is 'min' the solutions with a minimum score are returned.
    scoring_function -- function that can be applied to a PairedRegion
        object and that returns a numerical score. 

    This function fills a matrix that calculates the optimal solution for the 
        pseudoknot-removal problem by storing optimal solutions for smaller
        sub-problems. 

    The number of cells in each DP matrix is the number of given paired
        regions times two, because there is one row and column for each
        start and end point of each region. Only the top-right half of 
        the matrix will be filled. 
    A row index is referred to as i (begin_idx in code), a column index
        is referred to as j (end_idx in code). The matrix is initialized
        on the diagonal (where i==j) with a list containing an empty
        solution (an empty PairedRegions object). A list is used because
        we keep track of all possible optimal choices.
    For each cell (i,j) where j>i we collect all the candidate-solutions
        as follows. 
        ** Add all solutions of the cell to the left, which contains the
        best solutions for the area from start point i to end point j-1.
        ** Add all solutions of the cell to the bottom, which contains
        the best solutions for the area from start point i+1 to end point j.
        ** If start point i and end point j are a start and end point of
        the same region, add all possible solutions from the cell to the
        bottom-left plus this region. The cell to the bottom-left contains
        the optimal solutions for the area from start point i+1 to end
        point j-1.
        ** If the lists of solutions at the cells to the left and bottom
        both contained anything different from the empty solution, we need
        to check two more things:
            ** For each combination of a solution in the left cell and a
            solution in the right cell, calculate the highest end point of
            the solution in the left cell and the lowest start point for the
            bottom cell. In a collection of paired regions, every region has
            an end point, and the highest end point is the largest number in
            the list of all end points. The lowest start value is calculated
            in a similar way.
            ** If the highest end point is lower than the lowest start point,
            it means both solutions are disjoint and can be added to form a
            better solution. Thus, merge the two solutions and add them to the
            list of candidate-solutions.
            ** Otherwise, the solutions are not disjoint, but because of the
            pseudoknots sub-solutions of the two solutions might be combined
            to form a better solution. Create a slider k (splitter in code)
            that runs from the lowest start point minus one to the highest
            end point plus one. For each cell (i,k), (k+1,j) merge all 
            possible solutions and add them to the list of candidate-solutions.
        ** Next, store in cell (i,j) all solutions with an optimal score.
        There might be one solution, or there might be multiple solutions.
        ** Finish the calculation when the top-right cell in the matrix is
        filled. This cell contains the optimal solutions for the given set
        of paired regions.
    """
    if goal not in ['max','min']:
        raise ValueError("goal has to be 'min' or 'max', but is '%s'"%(goal))
    prs = paired_regions
    num_cells = len(prs)*2

    # pre-calculate scores 
    for pr in prs:
        pr.score(scoring_function)

    # create and initialize matrix
    result = empty_matrix(num_cells)
    
    # create some lookup dictionaries
    # enumerated start and end points
    enum_boundaries = prs.enumeratedBoundaries()
    # inverted enumerated start/end points {pr.Start/End: position in list}
    inv_enum_boundaries = prs.invertedEnumeratedBoundaries() 
    pos_to_pr = prs.byStartEnd() # {(pr.Start, pr.End): PairedRegion}
    
    # fill the matrix 
    for end_idx in range(num_cells):
        for begin_idx in range(end_idx-1, -1, -1):
            # look up sequence positions that match indices
            begin_pos = enum_boundaries[begin_idx]
            end_pos = enum_boundaries[end_idx]
            # look up solutions in left and bottom cells 
            left_cell = result[begin_idx][end_idx-1]
            bottom_cell = result[begin_idx+1][end_idx]
            # collect candidates
            candidates = []
            
            # add solutions from the left cell
            for sol in bottom_cell:
                candidates.append(sol)
            # add solutions from the bottom cell
            for sol in left_cell:
                candidates.append(sol)
            # if begin_pos and end_pos are paired:
            # ==> add left_bottom + this region
            if (begin_pos, end_pos) in pos_to_pr:
                this_region = pos_to_pr[(begin_pos, end_pos)]
                bottom_left = result[begin_idx+1][end_idx-1]
                for sol in bottom_left:
                    candidates.append(\
                    PairedRegions(sol+PairedRegions([this_region])))
            
            # if we have a solution in the left and in the bottom cell:
            if left_cell !=[[]] and bottom_cell != [[]]:

                # check whether they can be added or iterate
                for sol1 in left_cell:
                    for sol2 in bottom_cell:
                        he_pos = sol1.highestEnd()
                        he_idx = inv_enum_boundaries[he_pos]
                        ls_pos = sol2.lowestStart()
                        ls_idx = inv_enum_boundaries[ls_pos]
                        
                        # If both solutions are disjoint
                        if he_pos < ls_pos:
                            candidates.append(sol1.merge(sol2))
                        else: # not disjoint ==> iterate
                            for splitter in range(ls_idx-1, he_idx+1):
                                cell_to_left = result[begin_idx][splitter]
                                cell_to_bottom = result[splitter+1][end_idx]
                                if cell_to_bottom == [[]]:
                                    break
                                for sub_sol1 in cell_to_left:
                                    for sub_sol2 in cell_to_bottom:
                                        both = sub_sol1.merge(sub_sol2)
                                        candidates.append(both)
            
            # select all the candidates of maximum length
            best_candidate = pick_multi_best(candidates, goal=goal)
            result[begin_idx][end_idx] = best_candidate

    # return the whole matrix 
    return result

def matrix_solutions(paired_regions, goal='max', scoring_function=num_bps):
    """Return the list of solutions in the top-right cell of the DP matrix

    paired_regions -- PairedRegions object

    This methods fills a dynamic programming matrix (by calling
        dp_matrix_multi) and returns the list of solutions in the 
        top-right cell.
    """
    return dp_matrix_multi(paired_regions, goal=goal,\
        scoring_function=scoring_function)[0][-1]

# =============================================================================
# DYNAMIC PROGRAMMING APPROACH
# =============================================================================

# DP function used to remove pseudoknots
def opt_all(pairs, return_removed=False, goal='max',\
    scoring_function=num_bps):
    """Return a list of pseudoknot-free Pairs objects

    pairs -- Pairs object or list of tuples. One base can only interact
        with one other base, otherwise an error will be raised.
    return_removed -- boolean, if True a list of tuples of
        (nested pairs, removed pairs) will be returned. 
        Default is False --> list of nested pairs only is returned.
    goal -- str, 'max' or 'min', if the goal is 'max' the routine
        returns the solutions maximizing the score, if the goal
        is 'min' the solutions with a minimum score are returned.
    scoring_function -- function that can be applied to a PairedRegion
        object and that returns a numerical score.

    OPTIMIZATION, ALL SOLUTIONS (OA) -- PSEUDOKNOT REMOVAL METHOD

    This method will find all nested structures with an optimal score.
        Since there might be multiple optimal solutions, the retun value
        is always a list. If there is only one solution, the list will
        contain a single element.
    The problem is solved by dynamic programming. For each clique of mutually
        conflicting paired regions a matrix is filled out, and the best
        solutions are added to the result. Non-conflicting regions are part
        of the solution, and don't need to be processed.
    The user can specify the goal (maximize or minimize) and the scoring
        function. If one specifies for example goal='max' and 
        scoring_function=num_bps, the routine finds the nested structures
        with the maximum number of base pairs.
    See documentation of dp_matrix_multi for recursion rules used in the 
        approach.
    """
    if not pairs.hasPseudoknots():
        return [pairs]

    prs = PairedRegionsFromPairs(pairs)
    id_to_bl = prs.byId()
    cm = ConflictMatrix(prs)

    nc_regions = prs.nonConflicting(cm=cm)
    cliques = prs.conflictCliques(cm=cm) 
    
    # basis for all nested structures are the non-conflicting regions
    result = [PairedRegions(nc_regions)]
    # resolve conflicts, store survivors and removed
    for cl in cliques: 
        new_result = []
        best = matrix_solutions(cl, goal=goal,\
            scoring_function=scoring_function)
        for best_sol in best:
            for prev_res in result:
                new_result.append(prev_res.merge(best_sol))
        result = new_result

    if return_removed:
        # collect the removed pairs for each solution
        surviving_ids = []
        for sol in result:
            surviving_ids.append(dict.fromkeys([pr.Id for pr in sol]))
        removed = []
        for sol, surv in zip(result, surviving_ids):
            rem = []
            for pr_id in id_to_bl:
                if pr_id not in surv:
                    rem.extend(id_to_bl[pr_id].Pairs)
            rem.sort()
            removed.append(rem)
        nested = [prs.toPairs() for prs in result]
        return zip(nested, removed)
    
    nested = [prs.toPairs() for prs in result] 
    return nested

# =============================================================================
# MAJORITY OF BASE PAIRS -- CONVENIENCE FUNCTIONS
# =============================================================================

def opt_single_random(pairs, return_removed=False, goal='max',\
    scoring_function=num_bps):
    """Return single pseudoknot-free Pairs object with an optimal score

    pairs -- Pairs object or list of tuples. One base can only interact
        with one other base, otherwise an error will be raised.
    return_removed -- boolean, if True a tuple of (nested pairs, removed pairs)
        will be returned. Default is False --> only nested pairs are returned.
    goal -- str, 'max' or 'min', if the goal is 'max' the routine
        returns the solutions maximizing the score, if the goal
        is 'min' the solutions with a minimum score are returned.
    scoring_function -- function that can be applied to a PairedRegion
        object and that returns a numerical score.

    There might be multiple nested structures with an optimal score.
        This method calculates all of them and returns one at random.
    The user can specify the goal (maximize or minimize) and the scoring
        function. If one specifies for example goal='max' and 
        scoring_function=num_bps, the routine finds the nested structures
        with the maximum number of base pairs.
    """
    nested_structs = opt_all(pairs, return_removed, goal,\
        scoring_function)
    return choice(nested_structs)

def opt_single_property(pairs, return_removed=False, goal='max',\
    scoring_function=num_bps):
    """Return single pseudoknot-free Pairs object with max number of bps

    pairs -- Pairs object or list of tuples. One base can only interact
        with one other base, otherwise an error will be raised.
    return_removed -- boolean, if True a tuple of (nested pairs, removed pairs)
        will be returned. Default is False --> only nested pairs are returned.
    goal -- str, 'max' or 'min', if the goal is 'max' the routine
        returns the solutions maximizing the score, if the goal
        is 'min' the solutions with a minimum score are returned.
    scoring_function -- function that can be applied to a PairedRegion
        object and that returns a numerical score.

    There might be multiple nested structures with an optimal score.
        This method calculates all of them and
        returns the best by examining some properties. The first criterion
        is the number of paired regions in the returned structure, the
        second is the average range of the regions, the third is the 
        average start value of the regions. It returns the structure
        with the minimum value for these three properties. If all properties
        are the same for multiple structures, an error is raises. I believe
        this can't happen, but if it does the behavior can be changed.
    The user can specify the goal (maximize or minimize) and the scoring
        function. If one specifies for example goal='max' and 
        scoring_function=num_bps, the routine finds the nested structures
        with the maximum number of base pairs.
    """
    nested_structs = opt_all(pairs, return_removed, goal,\
        scoring_function)
    lookup = {}
    if return_removed:
        for p, p_rem in nested_structs:
            prs = PairedRegionsFromPairs(p)
            num_regions = len(prs)
            avg_range = average([pr.range() for pr in prs])
            avg_start = average([pr.Start for pr in prs])
            three = (num_regions, avg_range, avg_start)
            if three not in lookup:
                lookup[three] = []
            lookup[three].append((p,p_rem))
    else:
        for p in nested_structs:
            prs = PairedRegionsFromPairs(p)
            num_regions = len(prs)
            avg_range = average([pr.range() for pr in prs])
            avg_start = average([pr.Start for pr in prs])
            three = (num_regions, avg_range, avg_start)
            if three not in lookup:
                lookup[three] = []
            lookup[three].append(p)
    min_key = min(lookup.keys())
    min_value = lookup[min_key]
    if len(min_value) == 1:
        return min_value[0]
    else:
        # believe this can never happen, but just to be sure...
        raise ValueError("Multiple solutions found with equal properties")


# =============================================================================
# CONFLICT-ELIMINATION APPROACHES
# =============================================================================

def find_max_conflicts(conflicting_ids, cm, id_to_pr):
    """Return region ID of the region involved in the most conflicts

    conflicting_ids -- list of PairedRegion IDs
    cm -- ConflictMatrix object
    id_to_pr -- dict of {region ID: PairedRegion}. Result of 
        PairedRegions.byId() method.

    This methods returns the region ID (out of conflicting_ids) involved
        in the most conflicts. If there is a single region with the most
        conflicts, return it. Otherwise compare all regions with the 
        max number of conflicts on their gain. Gain is the length of
        the region minus the cumulative length of all of its conflicting
        regions. Return the one with the minimum gain. If both properties
        are equal, return the region that starts closest to the 3' end.
    """
    number_of_conflicts = {}
    for pr_id in conflicting_ids:
        noc = len(cm.conflictsOf(pr_id))
        if noc not in number_of_conflicts:
            number_of_conflicts[noc] = []
        number_of_conflicts[noc].append(pr_id)
    max_noc = max(number_of_conflicts.keys())
    max_ids = number_of_conflicts[max_noc]
    if len(max_ids) == 1:
        return max_ids[0]
    else:
        len_diffs = {}
        for pr_id in max_ids:
            pr_len = id_to_pr[pr_id].Length
            conf_len = sum([id_to_pr[i].Length for i in cm.conflictsOf(pr_id)])
            diff = pr_len - conf_len
            if diff not in len_diffs:
                len_diffs[diff] = []
            len_diffs[diff].append(pr_id)
        min_ld = min(len_diffs.keys())
        min_ids = len_diffs[min_ld]
        if len(min_ids) == 1:
            return min_ids[0]
        else:
            start_vals = {}
            for pr_id in min_ids:
                start = id_to_pr[pr_id].Start
                start_vals[start] = pr_id
            max_start = max(start_vals.keys())
            return start_vals[max_start]

def find_min_gain(conflicting_ids, cm, id_to_pr):
    """Return region ID of the region with the minimum gain

    conflicting_ids -- list of PairedRegion IDs
    cm -- ConflictMatrix object
    id_to_pr -- dict of {region ID: PairedRegion}. Result of 
        PairedRegions.byId() method.

    This methods returns the region ID (out of conflicting_ids) of the region
        that has the minimum gain. Gain is the length of the region minus
        the cumulative length of all of its conflicting regions. It expresses
        how many base pairs are gained if this region is kept and all
        of its conflicts have to be removed. If its gain is positive, it is 
        favorable to keep this region. If its gain is negative, it is better
        to remove this region and keep its conflicts. If there are multiple
        regions with the minimal gain, the one involved in the most conflicts
        is returned. If both properties are equal, the method returns the
        region that starts closest to the 3' end.
    """
    len_diffs = {}
    for pr_id in conflicting_ids:
        pr_len = id_to_pr[pr_id].Length
        conf_len = sum([id_to_pr[i].Length for i in cm.conflictsOf(pr_id)])
        diff = pr_len - conf_len
        if diff not in len_diffs:
            len_diffs[diff] = []
        len_diffs[diff].append(pr_id)
    min_ld = min(len_diffs.keys())
    min_ids = len_diffs[min_ld]
    if len(min_ids) == 1:
        return min_ids[0]
    else:
        number_of_conflicts = {}
        for pr_id in min_ids:
            noc = len(cm.conflictsOf(pr_id))
            if noc not in number_of_conflicts:
                number_of_conflicts[noc] = []
            number_of_conflicts[noc].append(pr_id)
        max_noc = max(number_of_conflicts.keys())
        max_ids = number_of_conflicts[max_noc]
        if len(max_ids) == 1:
            return max_ids[0]
        else:
            start_vals = {}
            for pr_id in min_ids:
                start = id_to_pr[pr_id].Start
                start_vals[start] = pr_id
            max_start = max(start_vals.keys())
            return start_vals[max_start]

def add_back_non_conflicting(paired_regions, removed):
    """Return new PairedRegions object and new dict of removed regions

    paired_regions -- PairedRegions object
    removed -- dict of {region_id: PairedRegion}

    Helper-function for conflict_elimination. 
    Circular removal might occur in conflict-elimination methods. It means
        that a particular region is removed and later in the process all
        of its conflicts are also removed, which result in an eliminated
        region that doens't conflict with any region in the solution anymore.
    This methods adds removed regions back into the solution if 
        they don't conflict with any region in the solution. The order in 
        which regions are tried to add is from 5' to 3' starting point.
    """
    id_to_pr = paired_regions.byId()
    new_removed = removed.copy()
    
    added = True
    # process removed from 5' to 3'
    order = [(pr.Start, pr.Id) for pr in new_removed.values()]
    order.sort() # from low start value to high start value
    while added:
        added = False
        for start, region_id in order:
            pr1 = new_removed[region_id]
            is_conflicting = False
            for pr2 in id_to_pr.values():
                if pr1.conflicting(pr2):
                    is_conflicting = True
                    new_removed[region_id] = pr1
                    break
            if not is_conflicting:
                id_to_pr[region_id] = pr1
                del new_removed[region_id]
                order = [(pr.Start, pr.Id) for pr in new_removed.values()]
                order.sort() # from low start value to high start value
                added = True
                break

    return PairedRegions(id_to_pr.values()), new_removed


# Conflict-elimination heuristic.
def conflict_elimination(pairs, sel_function, add_back=True,\
    return_removed=False):
    """Return pseudoknot-free Pairs object

    pairs -- Pairs object or list of tuples
    sel_function -- function that takes a list of IDs of conflicting regions,
        a conflict matrix and a dict of {region_id: PairedRegion} and returns
        the ID of a paired region that has to be removed.
    add_back -- boolean, if True regions that are removed but not conflicting
        at the end because of circular removal are added back into the 
        solution. If False, regions are only removed. This choice might result
        in too many regions being removed. Default value is True.
    return_removed -- boolean, if True a tuple of (nested pairs, removed pairs)
        will be returned. Default is False --> only nested pairs are returned.
   
    CONFLICT ELIMINATION -- PSEUDOKNOT REMOVAL METHOD
    EC -- sel_function=find_max_conflicts
    EG -- sel_functino=find_min_gain

    This is the general conflict-elimination function that should
        be used to remove pseudoknots from a knotted RNA structure.
    This algorithm removes paired regions one at the time. The order is 
        in which regions are removed is specified by the selection function.
    Different selection functions can be specified. Selection functions should
        take a list of conflicting IDS, a ConflictMatrix and a dict of
        {Region ID: PairedRegion} as input and they should return a 
        single PairedRegion ID.
    Two selection functions are available: find_max_conflicts and 
        find_min_gain. See their documentation for specifications.
    """
    prs = PairedRegionsFromPairs(pairs)
    id_to_pr = prs.byId()
    cm = ConflictMatrix(prs)
    removed = {}

    conf = cm.conflicting()
    while conf:
        to_remove = sel_function(conf, cm, id_to_pr)
        removed[to_remove] = id_to_pr[to_remove]
        prs.remove(id_to_pr[to_remove])
        id_to_pr = prs.byId()
        cm = ConflictMatrix(prs)
        conf = cm.conflicting()
    # potential circular removal: add regions back in
    if add_back:
        # collect IDs of non-conflicting removed regions
        prs, removed = add_back_non_conflicting(prs, removed)
        
    if return_removed:
        rem = PairedRegions(removed.values()).toPairs()
        return prs.toPairs(), rem
    return prs.toPairs()


# =============================================================================
# INCREMENTAL APPROACHES
# =============================================================================

# Incremental in order (IO) method
def inc_order(pairs, reversed=False, return_removed=False):
    """Return pseudoknot-free Pairs object

    pairs -- Pairs object or list of tuples. One base can only interact
        with one other base, otherwise an error will be raised.
    reversed -- boolean, indication whether the algoritm adds paired
        regions from 5' to 3' start values or from 3' to 5' end values. 
        If False, order is 5' to 3', if True, order is 3' to 5'. Default
        is False.
    return_removed -- boolean, if True a tuple of (nested pairs, removed pairs)
        will be returned. Default is False --> only nested pairs are returned.

    INCREMENTAL IN ORDER (IO) -- PSEUDOKNOT REMOVAL METHOD

    This algorithm treats all the paired regions in order, either starting
    at he 5' end (reversed=F) or at the 3' end (reversed=T). It accepts
    all the non-conflicting regions; If a region conflicts with an already
    added region, it is excluded from the solution.
    """
    prs = PairedRegionsFromPairs(pairs)
    id_to_pr = prs.byId()
    cm = ConflictMatrix(prs)

    if reversed:
        by_pos = [(pr.End, pr) for pr in prs]
        by_pos.sort()
        by_pos.reverse()
    else:
        by_pos = [(pr.Start, pr) for pr in prs]
        by_pos.sort()

    excluded = {}
    result = PairedRegions()
    for pos, pr in by_pos:
        if pr.Id in excluded:
            continue
        result.append(pr)
        for k in cm.conflictsOf(pr.Id):
            excluded[k] = True

    if return_removed:
        removed = Pairs([])
        for pr_id in excluded:
            removed.extend(id_to_pr[pr_id].Pairs)
        removed.sort()
        return result.toPairs(), removed
    return result.toPairs()
    
# Incremental by length (IL) method
def inc_length(pairs, reversed=False, return_removed=False):
    """Return pseudoknot-free Pairs object

    pairs -- Pairs object or list of tuples. One base can only interact
        with one other base, otherwise an error will be raised.
    reversed -- boolean. In case of equal lengths, all paired regions
        are processed from 5' to 3' starting position. 
        If reversed is True, regions are processed from 3' to 5' starting
        position.
    return_removed -- boolean, if True a tuple of (nested pairs, removed pairs)
        will be returned. Default is False --> only nested pairs are returned.

    INCREMENTAL BY LENGTH (IL) -- PSEUDOKNOT REMOVAL METHOD

    This algorithm will process the paired regions from the longest
        to the shortest. In case there are multiple regions of the same length,
        the one on the 5' side is added first if reversed=False (3' side is
        preferred if reversed=True). All paired regions that are conflicting
        with an already-added region are excluded.
    """
    prs = PairedRegionsFromPairs(pairs)
    id_to_pr = prs.byId()
    # create conflict matrix to lookup the conflicts
    cm = ConflictMatrix(prs)

    length_pos_data = {} # dict of {region_length: [(pr.start, pr)]}
    for pr in prs:
        if pr.Length not in length_pos_data:
            length_pos_data[pr.Length] = []
        length_pos_data[pr.Length].append((pr.Start, pr))
    for v in length_pos_data.values():
        v.sort()
        if reversed:
            v.reverse()
    
    excluded = {}
    result = PairedRegions()
   
    lengths = length_pos_data.keys()
    lengths.sort()
    lengths.reverse() # longest regions first

    for pr_len in lengths:
        for pr_start, pr in length_pos_data[pr_len]:
            if pr.Id not in excluded:
                result.append(pr)
                # use the conflict matrix to determine which regions to exclude
                for k in cm.conflictsOf(pr.Id):
                    excluded[k] = True

    if return_removed:
        removed = Pairs([])
        for pr_id in excluded:
            removed.extend(id_to_pr[pr_id].Pairs)
        removed.sort()
        return result.toPairs(), removed
    return result.toPairs()

# Incremental by range (IR) method
def inc_range(pairs, reversed=False, return_removed=False):
    """Return pseudoknot-free Pairs object

    pairs -- Pairs object or list of tuples. One base can only interact
        with one other base, otherwise an error will be raised.
    reversed -- boolean. If reversed is True: in case of two regions
        with the same range the region that starts closest to the 5' side
        is added first. If reversed is False: the region that starts closest
        to the 3' side is added first. Default is False (5' region preferred).
    return_removed -- boolean, if True a tuple of (nested pairs, removed pairs)
        will be returned. Default is False --> only nested pairs are returned.
    
    INCREMENTAL BY RANGE (IR) -- PSEUDOKNOT REMOVAL METHOD

    This algorithm will process the paired regions from the one with the 
        shortest range to the one with the longest range. The range of 
        a region is defined as the distance between the highest upstream
        position and the lowest downstream position (-1); in other words, it
        is the number of unpaired bases in the haripin if this region were
        the only paired region in the structure. 
    In case there are multiple regions of the same length,
        the one that starts closest to the 5' side is added first if 
        reversed=False (starting at the 3' side is preferred if reversed=True).
    All paired regions that are conflicting with an already-added region
        are excluded.
    """
    
    prs = PairedRegionsFromPairs(pairs)
    id_to_pr = prs.byId()
    # create conflict matrix to lookup the conflicts
    cm = ConflictMatrix(prs)

    range_pos_data = {} # dict of {region_range: [(pr.start, pr)]}
    for pr in prs:
        rr = pr.range() # region range
        if rr not in range_pos_data:
            range_pos_data[rr] = []
        range_pos_data[rr].append((pr.Start, pr))
    for v in range_pos_data.values():
        v.sort()
        if reversed:
            v.reverse()

    ranges = range_pos_data.keys()
    ranges.sort()

    result = PairedRegions()
    excluded = {}
    for rr in ranges:
        for pr_start, pr in range_pos_data[rr]:
            if pr.Id not in excluded:
                result.append(pr)
                for k in cm.conflictsOf(pr.Id):
                    excluded[k] = True

    if return_removed:
        removed = Pairs([])
        for pr_id in excluded:
            removed.extend(id_to_pr[pr_id].Pairs)
        removed.sort()
        return result.toPairs(), removed
    return result.toPairs()

# =============================================================================
# NUSSINOV RESTRICTED
# =============================================================================

def nussinov_fill(pairs, size):
    """Return filled dynamic programming search matrix with number of base pairs

    pairs -- Pairs object or list of tuples, should be directed (up,down)
    size -- int, number of rows and columns in the matrix (should be
        at least as much as the highest base paired position)

    Applies Nussinov-Jacobson algorithm restricted to input list of pairs.
    This function records the number of base pairs in the optimal
        (sub)solution.
    """
    bp_dict = dict.fromkeys(pairs)
    m = zeros((size, size), int)
    for j in range(size):
        for i in range(j-1,-1,-1):
            m[i,j] = m[i+1,j] # i unpaired
            if m[i,j-1] > m[i,j]: # j unpaired
                m[i,j] = m[i,j-1]
            if (i,j) in bp_dict and m[i+1,j-1]+1 > m[i,j]: # (i,j) pair
                m[i,j] = m[i+1,j-1]+ 1
            for k in range(i+1,j-1): # bifurcation
                if m[i,k] + m[k+1,j] > m[i,j]:
                    m[i,j] = m[i,k]+m[k+1,j]
    return m

def nussinov_traceback(m, i, j, pairs):
    """Return set of base pairs: nested structure with max number of pairs

    m -- filled DP search matrix
    i -- int, row coordinate where traceback should start, normally 0
    j -- int, column coordinate where traceback should start, normally length-1
    pairs -- Pairs object of list of tuples, should be directed (up, down),
        expect the same list as at the fill stage.

    Traceback procedure of the Nussinov-Jacobson algorithm that returns
        a single solution containing the maximum number of base pairs.
    """
    bp_dict = dict.fromkeys(pairs)
    if m[i,j] == 0: #or if i>=j:
        return set()
    if (i,j) in bp_dict and m[i+1,j-1] + 1 == m[i,j]:
        return set([(i,j)]) | nussinov_traceback(m, i+1, j-1, pairs)
    for k in range(i,j):
        if m[i,j] == m[i,k] + m[k+1,j]:
            return nussinov_traceback(m,i,k,pairs) | \
                nussinov_traceback(m,k+1,j,pairs)

def nussinov_restricted(pairs, return_removed=False):
    """Return nested Pairs object containing the maximum number of pairs

    pairs -- Pairs object or list of tuples [(1,10),(2,9),...]. List
        of base pairs has to be conflict-free (one base can only
        pair with one other base)
    return_removed -- boolean, if True a list of tuples of
        (nested pairs, removed pairs) will be returned. 
        Default is False --> list of nested pairs only is returned.

    This function is a modification of the original Nussinov-Jacobson
        algorithm, restricted to the given list of base pairs.
        It calculated a nested structure containing the maximum
        number of base pairs.

    NOTE: This function is very slow. If have have many base pairs in 
        the list, the size of the matrix grows quickly. We recoomend using
        the opt_all (OA) function for larger problems.
    """
    # pairs have to be conflict-free, directed, and sorted
    p = Pairs(pairs)
    if p.hasConflicts():
        raise ValueError("Cannot handle base pair conflicts")
    p = p.directed()
    p.sort()
    if not p.hasPseudoknots():
        nested = p  # if not Pseudoknots: return structure
    else:
        paired_positions = [x for x,y in p] + [y for x,y in p]
        paired_positions.sort()
        if max(paired_positions) > 200:
            # if the 'sequence' is longer than 200, map the numbers
            # to remove the 'unpaired' positions. Avoid waste of space
            mapped_back = dict([(x,y) for x,y in enumerate(paired_positions)])
            mapped = dict([(y,x) for x,y in enumerate(paired_positions)])
            mapped_pairs = []
            for x,y in pairs:
                mapped_pairs.append((mapped[x],mapped[y]))
            max_idx = len(paired_positions)+1
            m = nussinov_fill(mapped_pairs, size=max_idx)
            t = nussinov_traceback(m, 0, max_idx-1, mapped_pairs)
            nested = []
            for x,y in t:
                nested.append((mapped_back[x], mapped_back[y]))
        else:
            # sequence short enough, fill matrix directly
            max_idx = max(filter(None, paired_positions))+1
            m = nussinov_fill(p, size=max_idx)
            nested = nussinov_traceback(m, 0, max_idx-1, p)
            nested = Pairs(list(nested))
        nested.sort()
    if return_removed:
        removed = Pairs([])
        for bp in p:
            if bp not in nested:
                removed.append(bp)
        return nested, removed
    else:
        return nested



if __name__ == "__main__":
    pass
