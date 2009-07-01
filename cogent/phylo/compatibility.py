"""Visualisation of phylogenetic compatibility within an alignment.
Jakobsen & Easteal, CABIOS 12(4), 1996
Jakobsen, Wilson & Easteal, Mol. Biol. Evol. 14(5), 1997 
"""

from __future__ import division
import numpy
import numpy.random
import operator

def order_to_cluster_similar(S, elts=None, start=None):
    """Order so as to keep the most similar parts adjacent to each other
    S is expected to be a square matrix, and the returned list of 
    ordinals is len(S) long."""
    position = {}
    unavailable = set()
    if elts is not None:
        if len(elts) < 2:
            return elts
        elts = set(elts)
        for x in range(len(S)):
            if x is not start and x not in elts:
                unavailable.add(x)
    if start is not None:
        position[start] = (False, [start])
                
    similarity = [numpy.unravel_index(p, S.shape) 
                for p in numpy.argsort(S, axis=None)]
    for (x, y) in similarity[::-1]:
        if x==y or x in unavailable or y in unavailable:
            continue
        (x_end, x_run) = position.get(x, (None, [x]))
        (y_end, y_run) = position.get(y, (None, [y]))
        if x_run is y_run:
            continue
        if x_end is not None:
            if not x_end:
                x_run.reverse()
            unavailable.add(x)
        if y_end is not None:
            if y_end:
                y_run.reverse()
            unavailable.add(y)
        run = x_run + y_run
        position[run[0]] = (False, run)
        position[run[-1]] = (True, run)
    if start is not None:
        if run[-1] == start:
            run.reverse()
        run = run[1:]
    return run

def order_tied_to_cluster_similar(S, scores):
    """Use similarity measure S to make similar elements
    adjacent, but only to break ties in the primary order
    defined by the list of scores"""
    assert S.shape == (len(scores), len(scores))
    (pos,) = numpy.where(numpy.diff(scores) != 0)
    pos = numpy.concatenate(([0],pos+1,[len(scores)]))
    new_order = []
    start = None
    for (a,b) in zip(pos[:-1],pos[1:]):
        sub_order = order_to_cluster_similar(S, range(a,b), start)
        new_order.extend(sub_order)
        start = new_order[-1]
    assert set(new_order) == set(range(len(scores))) 
    return new_order

def bit_encode(x, _bool2num=numpy.array(["0","1"]).take):
    """Convert a boolean array into an integer"""
    return int(_bool2num(x).tostring(), 2)

def binary_partitions(alignment):
    """Returns (sites, columns, partitions) 
    sites[alignment position number] = informative column number
    columns[informative column number] = distinct partition number
    partitions[distinct partition number] = (partition, mask) as ints
    """
    sites = []
    columns = []
    partitions = []
    partition_index = {}
    for column in alignment.Positions:
        column = numpy.array(column)
        (A, T, C, G, R, Y, W, S, U) = [
            bit_encode(column == N) for N in "ATCGRYWSU"]
        T |= U
        for split in [([A, G, R], [C, T, Y]), ([A, T, W], [G, C, S])]:
            halves = []
            for char_group in split:
                X = reduce(operator.or_, char_group)
                if not (X & (X - 1)):
                    break  # fewer than 2 bits set in X
                halves.append(X)
            else:
                (X, Z) = sorted(halves)
                partition = (X,X|Z)
                if partition not in partition_index:
                    partition_index[partition] = len(partitions)
                    partitions.append(partition)
                sites.append(len(columns))
                columns.append(partition_index[partition])
                break  # if R/Y split OK no need to consider W/S split. 
    return (sites, columns, partitions)

def min_edges(columns):
    """Given two boolean arrays each representing an informative alignment 
    position, there are 4 possible combinations for each sequence: 
    TT, TF, FT and FF.
    If N of these 4 possibilities are found then there must be at least 
    N-1 tree edges on which mutations occured"""  
    N = len(columns)
    result = numpy.ones([N, N], int)
    for i in range(0, N-1):
        (a, mask_a) = columns[i]
        for j in range(i+1, N):
            (b, mask_b) = columns[j]
            mask = mask_a & mask_b
            (na, nb) = (~a, ~b)
            combos = [c & mask for c in [a&b, a&nb, na&b, na&nb]]
            combos = [c for c in combos if c]
            result[i,j] = result[j,i] = len(combos) - 1
    return result

def neighbour_similarity_score(matrix):
    left = matrix[:-1]
    right = matrix[1:]
    upper = matrix[:,:-1]
    lower = matrix[:,1:]
    same = (lower == upper).sum() + (left == right).sum()
    neighbours = numpy.product(left.shape)+numpy.product(upper.shape)
    return same / neighbours

def shuffled(matrix):
    assert matrix.shape == (len(matrix), len(matrix)), matrix.shape
    index = numpy.random.permutation(numpy.arange(len(matrix)))
    return matrix[index,:][:,index]
    
def nss_significance(matrix, samples=10000):
    score = neighbour_similarity_score(matrix)
    scores = numpy.empty([samples])
    for i in range(samples):
        s = neighbour_similarity_score(shuffled(matrix))
        scores[i] = s
    scores.sort()
    p = (samples-scores.searchsorted(score)+1) / samples
    return (score, sum(scores)/samples, p)
    
def inter_region_average(a):
    return a.sum()/numpy.product(a.shape)
    
def intra_region_average(a):
    d = numpy.diag(a)    # ignore the diagonal
    return (a.sum()-d.sum())/(numpy.product(a.shape)-len(d))
    
def main(alignment, display=False, samples=0):
    print "%s sequences in %s bp alignment" % (
            alignment.getNumSeqs(), len(alignment))
    (sites, columns, partitions) = binary_partitions(alignment)
    print "%s unique binary partitions from %s informative sites" % (
            len(partitions), len(sites))
    partpart = min_edges(partitions)      # [partition,partition]
    partimatrix = partpart[columns,:]     # [site, partition]
    sitematrix = partimatrix[:,columns]   # [site, site]
    
    # RETICULATE, JE 1996
    
    compatiblity = sitematrix <= 2
    print "Overall compatibility %.6f" % intra_region_average(compatiblity)
    if samples == 0:
        print "Neighbour similarity score = %.6f" % \
                neighbour_similarity_score(compatiblity)
    else:
        print "Neighbour similarity = %.6f, avg random = %.6f, p < %s" % \
                nss_significance(compatiblity, samples=samples)
    if display:
        import pylab
        pylab.matshow(compatiblity, cmap=pylab.cm.gray)
        pylab.show()
        
    # PARTIMATRIX, JWE 1997
    
    # Remove the incomplete partitions with gaps or other ambiguities
    mask = 2**alignment.getNumSeqs()-1
    complete = [i for (i,(x, xz)) in enumerate(partitions) if xz==mask]
    partimatrix = partimatrix[:,complete]
    # For scoring/ordering purposes, also remove the incomplete sequences
    complete_columns = [i for (i,c) in enumerate(columns) if c in complete]
    scoreable_partimatrix = partimatrix[complete_columns, :]
    
    # Order partitions by increasing conflict score
    conflict = (scoreable_partimatrix > 2).sum(axis=0)
    conflict_order = numpy.argsort(conflict)
    partimatrix = partimatrix[:, conflict_order]
    scoreable_partimatrix = partimatrix[complete_columns, :]
    support = (scoreable_partimatrix == 1).sum(axis=0)
    consist = (scoreable_partimatrix <= 2).sum(axis=0)
    conflict = (scoreable_partimatrix > 2).sum(axis=0)
    
    # Similarity measure between partitions
    pp1 = (scoreable_partimatrix <= 2).astype(int).T
    npp = (scoreable_partimatrix > 2).astype(int).T
    O = numpy.inner(pp1, pp1) + numpy.inner(npp, npp)
    s = 1.0*len(complete_columns)
    O = O.astype(float) / s
    p,q = consist/s, conflict/s
    E = numpy.outer(p,p) + numpy.outer(q,q)
    S = (O-E)/numpy.sqrt(E*(1-E)/s)
    
    # Order partitions for better visual grouping
    if "order_by_conflict":
        order = order_tied_to_cluster_similar(S, conflict)
    else:
        order = order_to_cluster_similar(S)
        half = len(order) // 2
        if sum(conflict[order[:half]]) > sum(conflict[order[half:]]):
            order.reverse()
    
    partimatrix = partimatrix[:, order]
            
    if display:
        import pylab
        # 3-colour code: support, compatible, conflict
        partishow = numpy.array([5,4,0]).take(partimatrix-1)
        pylab.matshow(partishow, cmap=pylab.cm.gray)
        pylab.show()
    
if __name__ == '__main__':
    from cogent import LoadSeqs, DNA
    import sys, optparse
    parser = optparse.OptionParser("usage: %prog [-ds] alignment")
    parser.add_option("-d", "--display", action="store_true", default=False, 
            dest="display", help="show matrices via matplotlib")
    parser.add_option("-s", "--samples",
                  dest="samples", default=10000, type="int",
                  help="samples for significance test")
    (options, args) = parser.parse_args()
    alignment = LoadSeqs(args[0], moltype=DNA)
    main(alignment, **vars(options))
    
