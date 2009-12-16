#!usr/bin/env python
"""nonmetric multidimensional scaling (nmds)

see, for example: Jan de Leeuw 2004 (monotone regression), 
Rencher 2002: Methods of multivariate analysis, and the original work: 
Kruskal 1964: Nonmetric multidimensional scaling
"""
from __future__ import division
from numpy import array, multiply, sum, zeros, size, shape, diag, dot, mean,\
    sqrt, transpose, trace, argsort
from numpy.random import seed, normal as random_gauss
from numpy.linalg import norm, svd
from operator import itemgetter
import cogent.maths.scipy_optimize as optimize
from cogent.cluster.metric_scaling import principal_coordinates_analysis

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

class NMDS(object):
    """Generates points using nonmetric scaling
    
    Takes as input an n by n distance/dissimilarity matrix (hereafter called a
    dissimilarity matrix for clarity), and a desired dimension (k).  Using the
    order of these dissimilarities only (the value is considered only to the
    extent that it determines the order), nonmetric_scaling constructs n
    points in k-dimensional space which correspond to the input matrix.
    
    The algorithm attempts to have the order of the pairwise distances between
    points correspond as closely as possible to the order of the dissimilarity
    matrix entries; if dissim[i,j] is the minimum dissimilarity, point i and 
    point j should be the two points closest together.
    
    The algorithm is random in nature, and is not guaranteed to converge to
    the best solution.  Furthermore, in general a dissimilarity matrix cannot
    be represented exactly by n points in k space, for k < n-1.
    
    Currently the convergence test is pretty basic, it just tests for a
    sufficiently small or negative relative improvement in stress, 
    or if a sufficiently tiny value of stress has been reached.
    Alternatively, it will stop after max_iterations 
    
    The basic algorithm is:
        - generate n points in k space
        - compute pairwise distances between these points
        - compare these distances with a set of pseudo-distances (dhats),
        which are constrained to increase in the same order as the input data
        - repeatedly adjust the points and dhats until the point distances
        have nearly the same order as the input data

    Note: increasing MIN_ABS_STRESS causes nans to return from stress fn
    """
    
    def __init__(self, dissimilarity_mtx, initial_pts="pcoa", 
        dimension=2, rand_seed=None, optimization_method=1, verbosity=1,
        max_iterations=50, setup_only=False, min_rel_improvement = 1e-3,
        min_abs_stress = 1e-5):
        """    
        Arguments:
        - dissimilarity_mtx: an n by n numpy float array representing the 
        pairwise dissimilarity of items.  0 on diagonals, symmetric under 
        (i,j) -> (j,i)
        - initial_pts: "random" => random starting points, "pcoa" => 
        pts from pcoa, or a numpy 2d array, ncols = dimension
        - dimension: the desired dimension k of the constructed 
        - rand_seed: used for testing
        - optimization_method: used when points are adjusted to minimize stress:
        0 => justin k's ad hoc method of steepest descent
        1 => cogent's scipy_optimize fmin_bfgs
        """
        self.min_rel_improvement = min_rel_improvement
        self.min_abs_stress = min_abs_stress

        if rand_seed != None:
            seed(rand_seed)
        
        self.verbosity = verbosity
        num_points = len(dissimilarity_mtx)
        point_range = range(num_points)
        self.dimension = dimension
        self.optimization_method = optimization_method
        
        self._calc_dissim_order(dissimilarity_mtx, point_range)
        # sets self.order
        # note that in the rest of the code, only the order matters, the values
        # of the dissimilarity matrix aren't used
        
        if initial_pts == "random":
            self.points = self._get_initial_pts(dimension, point_range)
        elif initial_pts == "pcoa":
            pcoa_pts, pcoa_eigs = principal_coordinates_analysis(\
                dissimilarity_mtx)
            order = argsort(pcoa_eigs)[::-1] # pos to small/neg
            pcoa_pts = pcoa_pts[order].T
            self.points = pcoa_pts[:,:2]
        else:
            self.points = initial_pts
        
        self.points = self._center(self.points)
        self._rescale()
        self._calc_distances() 
        # sets self.dists, ordered according to self.order
        # dists relates to points, not to input data
        
        self._update_dhats()
        # self.dhats are constrained to be monotonic
        
        self._calc_stress()
        # self.stress is calculated from self.dists and self.dhats
        
        self.stresses = [self.stress]
        # stress is the metric of badness of fit used in this code
        # index 0 is the initial stress, with a initial set of 
        # datapoints. index 1 corresponds to iteration 0 of the loop below
        
        if setup_only:
            return
            
        for i in range(max_iterations):
            self._move_points()
            self._calc_distances()
            self._update_dhats()
            self._calc_stress()
            self.stresses.append(self.stress)
            if self.verbosity >= 1:
                print("nonmetric broad iteration, stress: ", i,
                self.stresses[-1])
            
            if (self.stresses[-2]-self.stresses[-1]) / self.stresses[-2] <\
                self.min_rel_improvement:
                if self.verbosity >= 1:
                    print "iteration improvement minimal. converged."
                break
            if (self.stresses[-1] < self.min_abs_stress):
                if self.verbosity >= 1:
                    print "stress below cutoff, done" 
                break
        
        # center and rotate the points, since pos, rotation is arbitrary
        # rotation is to align to principal axes of self.points
        self.points = self._center(self.points)
        u,s,vh = svd(self.points, full_matrices=False)
        S = diag(s)
        self.points = dot(u,S)
        # normalize the scaling, which should not change the stress
        self._rescale()
    
    
    def getPoints(self):
        """Returns (ordered in a list) the n points in k space 
        
        these are the algorithm's attempt at points corresponding to the input
        order of dissimilarities.  Returns a numpy 'd' mtx, points in rows
        """
        return self.points
    
    def getStress(self):
        """Returns a measure of the badness of fit

        not in percent, a typical number for 20 datapoints is .12"""
        return self.stresses[-1]
        
    def getDimension(self):
        """returns the dimensions in which the constructed points lie"""
        return self.dimension
        
    def _center(self, mtx):
        """translate all data (rows of the matrix) to center on the origin
        
        returns a shifted version of the input data.  The new matrix is such 
        that the center of mass of the row vectors is centered at the origin.  
        Returns a numpy float ('d') array
        """
        result = array(mtx, 'd')
        result -= mean(result, 0) 
        # subtract each column's mean from each element in that column
        return result
    
    def _calc_dissim_order(self, dissim_mtx, point_range):
        """calculates the order of the dissim_mtx entries, puts in self.order
        
        First creates a list of dissim elements with structure [i, j, value],
        then sorts that by value and strips the value subelemnt.
        i and j correspond to the row and column of the input dissim matrix 
        """
        
        dissim_list = []
        for i in point_range:
            for j in point_range:
                if j > i:
                    dissim_list.append([i, j, dissim_mtx[i,j]])
        dissim_list.sort(key = itemgetter(2))
        for elem in dissim_list:
            elem.pop()
        self.order = dissim_list

    def _get_initial_pts(self, dimension, pt_range):
        """Generates points randomly with a gaussian distribution (sigma = 1)
        """
        
        # nested list comprehension.  Too dense for good readability?
        points = [[random_gauss(0., 1) for axis in range(dimension)] \
            for pt_idx in pt_range]
        return array(points, 'd')

    def _calc_distances(self):
        """Returns a list of pairwise distances, ordered by self.order
        """
        
        # cProfile indicates this call is the speed bottleneck
        dists = [norm(self.points[i] - self.points[j]) for i, j in self.order]
        
        self.dists = array(dists, 'd')


    def _update_dhats(self):
        """updates self.dhats based on self.dists data"""
        dhats = self.dists.copy()
        dhats = self._do_monotone_regression(dhats)
        self.dhats = array(dhats, 'd')
        
    def _do_monotone_regression(self, dhats):
        """Performs a monotone regression on dhats, returning the result
        
        Assuming the input dhats are the values of the pairwise point 
        distances, this algorithm minimizes the stress while enforcing
        monotonicity of the dhats.
        Jan de Leeuw 2004 (monotone regression) has a rough outline of the
        algorithm.  Basically, as we proceed along the ordered list, 
        if an element is smaller than its preceeding one, the two are averaged
        and grouped together in a block.  The process is repeated until
        the blocks are monotonic, that is block i <= block i+1.
        
        As the list gets long, it is likely that a deviation from monotinicity 
        would not propagate back to the beginning.  
        The current implementation restarts from the beginning
        after fixing each non-monotonicity, potentially needlessly rechecking 
        many blocks.  It could be rewritten to fix the entire list in one pass.
        However, the _move_points step is more likely to be the 
        speed bottleneck here.
        """
        
        initial_len = len(dhats)
        blocklist = [[dhat, 1] for dhat in dhats] 
        # each element is [value, blocksize]
        
        block_idx = 0
        while True:
            if block_idx == (len(blocklist) - 1):
                # we are at the last block => list is monotonic
                are_blocks_monotonic = True
                break
            
            if blocklist[block_idx][0] > blocklist[block_idx+1][0]:
                # fix non-monotonicity and restart from beginning of blocklist
                b1_val = blocklist[block_idx][0]
                b1_size = blocklist[block_idx][1]
                b2_val = blocklist[block_idx+1][0]
                b2_size = blocklist[block_idx+1][1]
                blocklist[block_idx][0] = (b1_val*b1_size + b2_val*b2_size)\
                    /(b1_size + b2_size)
                blocklist[block_idx][1] = b1_size + b2_size
                blocklist.pop(block_idx+1)
                block_idx = 0
            else:
                block_idx += 1
                
        # remake dhats as a flat list, not a blocklist
        result_dhats = []
        for block in blocklist:
            for j in range(block[1]):
                result_dhats.append(block[0])
                
        if len(result_dhats) != initial_len:
            raise RuntimeError("monotone regression changed list size")
        
        return result_dhats
        
    def _calc_stress(self):
        """calculates the stress, or badness of fit from self.dhats, self.dists
        """
        diff = self.dists - self.dhats
        top = sum(multiply(diff, diff))
        self.stress = sqrt(top/sum(multiply(self.dists, self.dists)))
        
    def _rescale(self):
        """ assumes centered, rescales to mean ot-origin dist of 1
        """
    
        factor = array([norm(vec) for vec in self.points]).mean()
        self.points = self.points/factor

    def _move_points(self):
        """ this attempts to move our points in such a manner as to minimize 
        the stress metric, keeping dhats fixed.  If the dists could be chosen 
        without constraints, by assigning each dist[i,j] = dhat[i,j], 
        stress would be zero.
        However, since the distances are computed from points, it is generally
        impossible to change the dists independantly of each other.
        
        a basic algorithm is:
        - move points
        - recompute dists
        - recompute stress
        - if stress decreased, continue in the same manner, otherwise 
        move points in a different manner
        
        self.points often serves as a starting point for optimizaion algorithms
        
        optimization algorithm 0 is justin's hack (steepest descent method)
        """

        
        if self.optimization_method == 0:
            self._steep_descent_move()
        
        elif self.optimization_method == 1:
            numrows, numcols = shape(self.points)
            pts = self.points.ravel().copy()
            optpts = optimize.fmin_bfgs(self._recalc_stress_from_pts, pts,
                disp=self.verbosity, maxiter=100, gtol=1e-3)
            self.points = optpts.reshape((numrows, numcols))
        else:
            raise ValueError
            

    def _steep_descent_move(self,
        rel_step_size=1./100, precision=.00001, max_iters=100):
        """moves self.points. goal: minimize the stress.
        
        Uses steepest descent method.
        
        This is currently an ad-hoc minimization routine, using the method
        of steepest descent.  The default parameters are only shown to work on
        a few simple cases, and aren't optimized.
        
        The gradient is calculated discretely, not via formula.  Each variable
        (there are n points * k dimensions of variables), is adjusted, 
        the stress measured, and the variable returned to its prior value.
        
        If a local minimum is larger than step_size, the algorithm cannot 
        escape.
        """
        num_rows, num_cols = shape(self.points)
        avg_point_dist = sum([norm(point) for point in self.points])/num_rows
        step_size = avg_point_dist*rel_step_size

            
        for iter in range(max_iters):
            
            # initial values
            prestep_stress = self.stress.copy()
            gradient = zeros((num_rows, num_cols))
            
            # get gradient
            for i in range(num_rows):
                for j in range(num_cols):
                    self.points[i,j] += step_size
                    self._calc_distances()
                    self._calc_stress()
                    delta_stress = self.stress - prestep_stress
                    gradient[i,j] = delta_stress/step_size
                    self.points[i,j] -= step_size
            
            grad_mag = norm(gradient)
            
            # step in the direction of the negative gradient
            for i in range(num_rows):
                for j in range(num_cols):
                    self.points[i,j] -= step_size*gradient[i,j]/grad_mag
            self._calc_distances()
            self._calc_stress()
            newstress = self.stress.copy()

            # choose whether to iterate again
            if abs((newstress - prestep_stress)/prestep_stress) < precision:
                if self.verbosity >= 1:
                    print("move pts converged after iteration: ", iter)
                break
            if iter == (max_iters - 1):
                if self.verbosity >= 1:
                    print("move pts didn't converge in ", max_iters)
        
    def _recalc_stress_from_pts(self, pts):
        """returns an updated value for stress based on input pts
        
        a special function for use with external optimization routines.
        pts here is a 1D numpy array"""
        numrows, numcols = shape(self.points)
        self.points = pts.reshape((numrows, numcols))
        self._calc_distances()
        self._calc_stress()
        return self.stress
    

def metaNMDS(iters, *args, **kwargs):
    """ runs NMDS, first with pcoa init, then iters times with random init

    returns NMDS object with lowest stress
    args, kwargs is passed to NMDS(), but must not have initial_pts
    must supply distance matrix
    """
    results = []
    kwargs['initial_pts'] = "pcoa"
    res1 = NMDS(*args,**kwargs)
    results.append(res1)
    kwargs['initial_pts'] = "random"
    for i in range(iters):
        results.append(NMDS(*args, **kwargs))
    stresses = [nmds.getStress() for nmds in results]
    bestidx = stresses.index(min(stresses))
    return results[bestidx]
