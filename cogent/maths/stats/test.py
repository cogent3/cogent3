#!/usr/bin/env python
"""Provides standard statistical tests. Tests produce statistic and P-value.
"""
from __future__ import division
import warnings
from cogent.maths.stats.distribution import chi_high, z_low, z_high, zprob, \
    t_high, t_low, tprob, f_high, f_low, fprob, binomial_high, binomial_low
from cogent.maths.stats.special import lgam, log_one_minus, one_minus_exp,\
    MACHEP
from cogent.maths.stats.ks import psmirnov2x, pkstwo
from cogent.maths.stats.kendall import pkendall, kendalls_tau
from cogent.maths.stats.special import Gamma

from numpy import array, asarray, transpose, ravel, take, nonzero, log, sum,\
        mean, cov, corrcoef, fabs, any, reshape, clip, nan, isnan, isinf, \
        sqrt, exp, median as _median
        #, std - currently incorrect
from numpy.random import permutation
from cogent.maths.stats.util import Numbers
from operator import add
from random import choice

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley", "Rob Knight", "Catherine Lozupone",
                    "Sandra Smit", "Micah Hamady", "Daniel McDonald",
                    "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class IndexOrValueError(IndexError, ValueError): pass

var = cov   #cov will calculate variance if called on a vector

def std_(x, axis=None):
    """Returns standard deviations by axis (similiar to numpy.std)

    The result is unbiased, matching the result from MLab.std
    """
    x = asarray(x)
    
    if axis is None:
        d = x - mean(x)
        return sqrt(sum(d**2)/(len(x)-1))
    elif axis == 0:
        result = []
        for col in range(x.shape[1]):
            vals = x[:,col]
            d = vals - mean(vals)
            result.append(sqrt(sum(d**2)/(len(x)-1)))
        return result
    elif axis == 1:
        result = []
        for row in range(x.shape[0]):
            vals = x[row,:]
            d = vals - mean(vals)
            result.append(sqrt(sum(d**2)/(len(x)-1)))
        return result
    else:
        raise ValueError, "axis out of bounds"
    
# tested only by std
def var(x, axis=None):
    """Returns unbiased standard deviations over given axis.

    Similar with numpy.std, except that it is unbiased. (var = SS/n-1)

    x: a float ndarray or asarray(x) is a float ndarray.
    axis=None: computed for the flattened array by default, or compute along an
     integer axis.

    Implementation Notes:
    Change the SS calculation from:
        SumSq(x-x_bar) to SumSq(x) - SqSum(x)/n
        See p. 37 of Zar (1999) Biostatistical Analysis.
    """
    x = asarray(x)
    #figure out sample size along the axis
    if axis is None:
        n = x.size
    else:
        n = x.shape[axis]
    #compute the sum of squares from the mean(s)
    sample_SS = sum(x**2, axis) - sum(x, axis)**2/n
    return sample_SS / (n-1)

def std(x, axis=None):
    """computed unbiased standard deviations along given axis or flat array.

    Similar with numpy.std, except that it is unbiased. (var = SS/n-1)

    x: a float ndarray or asarray(x) is a float ndarray.
    axis=None: computed for the flattened array by default, or compute along an
      given integer axis.
    """
    try:
        sample_variance = var(x, axis=axis)
    except IndexError, e: #just to avoid breaking the old test code
        raise IndexOrValueError(e)
    return sqrt(sample_variance)

def median(m, axis=None):
    """Returns medians by axis (similiar to numpy.mean)

    numpy.median does not except an axis parameter. Is safe for substition for 
    numpy.median
    """
    median_vals = []
    rows, cols = m.shape

    if axis is None:
        return _median(ravel(m))
    elif axis == 0:
        for col in range(cols):
            median_vals.append(_median(m[:,col]))
    elif axis == 1 or axis == -1:
        for row in range(rows):
            median_vals.append(_median(m[row,:]))
    else:
        raise ValueError, "axis(=%s) out of bounds" % axis

    return array(median_vals)

class ZeroExpectedError(ValueError):
    """Class for handling tests where an expected value was zero."""
    pass
   
def G_2_by_2(a, b, c, d, williams=1, directional=1):
    """G test for independence in a 2 x 2 table.

    Usage: G, prob = G_2_by_2(a, b, c, d, willliams, directional)

    Cells are in the order:
    
        a b
        c d
    
    a, b, c, and d can be int, float, or long.
    williams is a boolean stating whether to do the Williams correction.
    directional is a boolean stating whether the test is 1-tailed.
    
    Briefly, computes sum(f ln f) for cells - sum(f ln f) for
    rows and columns + f ln f for the table.
    
    Always has 1 degree of freedom

    To generalize the test to r x c, use the same protocol:
    2*(cells - rows/cols + table), then with (r-1)(c-1) df.

    Note that G is always positive: to get a directional test,
    the appropriate ratio (e.g. a/b > c/d) must be tested
    as a separate procedure. Find the probability for the
    observed G, and then either halve or halve and subtract from
    one depending on whether the directional prediction was
    upheld. 
    
    The default test is now one-tailed (Rob Knight 4/21/03).

    See Sokal & Rohlf (1995), ch. 17. Specifically, see box 17.6 (p731).
    """
    cells = [a, b, c, d]
    n = sum(cells)
    #return 0 if table was empty
    if not n:
        return (0, 1)
    #raise error if any counts were negative
    if min(cells) < 0:
        raise ValueError, \
        "G_2_by_2 got negative cell counts(s): must all be >= 0."
    
    G = 0
    #Add x ln x for items, adding zero for items whose counts are zero
    for i in filter(None, cells):
        G += i * log(i)
    #Find totals for rows and cols
    ab = a + b
    cd = c + d
    ac = a + c
    bd = b + d
    rows_cols = [ab, cd, ac, bd]
    #exit if we are missing a row or column entirely: result counts as
    #never significant
    if min(rows_cols) == 0:
        return (0, 1)
    #Subtract x ln x for rows and cols
    for i in filter(None, rows_cols):
        G -= i * log(i)
    #Add x ln x for table
    G += n * log(n)
    #Result needs to be multiplied by 2 
    G *= 2

    #apply Williams correction
    if williams:
        q = 1 + ((  ( (n/ab) + (n/cd) ) -1 ) * ( ( (n/ac) + (n/bd) ) -1))/(6*n)
        G /= q

    p = chi_high(max(G,0), 1)
    
    #find which tail we were in if the test was directional
    if directional:
        is_high =  ((b == 0) or (d != 0 and (a/b > c/d)))
        p = tail(p, is_high)
        if not is_high:
            G = -G
    return G, p

def safe_sum_p_log_p(a, base=None):
    """Calculates p * log(p) safely for an array that may contain zeros."""
    flat = ravel(a)
    nz = take(flat, nonzero(flat)[0])
    logs = log(nz)
    if base:
        logs /= log(base)
    return sum(nz * logs, 0)

def G_ind(m, williams=False):
    """Returns G test for independence in an r x c table.
    
    Requires input data as a numpy array. From Sokal and Rohlf p 738.
    """
    f_ln_f_elements = safe_sum_p_log_p(m)
    f_ln_f_rows = safe_sum_p_log_p(sum(m,0))
    f_ln_f_cols = safe_sum_p_log_p(sum(m,1))
    tot = sum(ravel(m))
    f_ln_f_table = tot * log(tot)

    df = (len(m)-1) * (len(m[0])-1)
    G = 2*(f_ln_f_elements-f_ln_f_rows-f_ln_f_cols+f_ln_f_table)
    if williams:
        q = 1+((tot*sum(1.0/sum(m,1))-1)*(tot*sum(1.0/sum(m,0))-1)/ \
            (6*tot*df))
        G = G/q
    return G, chi_high(max(G,0), df)


def calc_contingency_expected(matrix):
        """Calculates expected frequencies from a table of observed frequencies

        The input matrix is a dict2D object and represents a frequency table
        with different variables in the rows and columns. (observed
        frequencies as values)
                                
        The expected value is calculated with the following equation:
            Expected = row_total x column_total / overall_total
                                            
        The returned matrix (dict2D) has lists of the observed and the 
        expected frequency as values
        """
        #transpose matrix for calculating column totals
        t_matrix = matrix.copy()
        t_matrix.transpose()

        overall_total = sum(list(matrix.Items))
        #make new matrix for storing results
        result = matrix.copy()

        #populate result with expected values
        for row in matrix:
            row_sum = sum(matrix[row].values())
            for item in matrix[row]:
                column_sum = sum(t_matrix[item].values())
                #calculate expected frequency
                Expected = (row_sum * column_sum)/overall_total
                result[row][item] = [result[row][item]]
                result[row][item].append(Expected)
        return result

def G_fit(obs, exp, williams=1):
    """G test for fit between two lists of counts.

    Usage: test, prob = G_fit(obs, exp, williams)
    
    obs and exp are two lists of numbers.
    williams is a boolean stating whether to do the Williams correction.
    
    SUM(2 f(obs)ln (f(obs)/f(exp)))
    
    See Sokal and Rohlf chapter 17.
    """
    k = len(obs)
    if k != len(exp):
        raise ValueError, "G_fit requires two lists of equal length."
    G = 0
    n = 0
    
    for o, e in zip(obs, exp):
        if o < 0:
            raise ValueError, \
            "G_fit requires all observed values to be positive."
        if e <= 0:
            raise ZeroExpectedError, \
            "G_fit requires all expected values to be positive."
        if o:   #if o is zero, o * log(o/e) must be zero as well.
            G += o * log(o/e)
            n += o
    
    G *= 2
    if williams:
        q = 1 + (k + 1)/(6*n)
        G /= q

    return G, chi_high(G, k - 1)

def G_fit_from_Dict2D(data):
    """G test for fit on a Dict2D

    data is a dict2D. Values are a list containing the observed
    and expected frequencies (can be created with calc_contingency_expected)
    """
    obs_counts = []
    exp_counts = []
    for item in data.Items:
        if len(item) == 2:
            obs_counts.append(item[0])
            exp_counts.append(item[1])
    g_val, prob = G_fit(obs_counts, exp_counts)
    return g_val, prob

def chi_square_from_Dict2D(data):
    """Chi Square test on a Dict2D

    data is a Dict2D. The values are a list of the observed (O)
    and expected (E) frequencies,(can be created with calc_contingency_expected)

    The chi-square value (test) is the sum of (O-E)^2/E over the items in data

    degrees of freedom are calculated from data as:
    (r-1)*(c-1) if cols and rows are both > 1
    otherwise is just 1 - the # of rows or columns
    (whichever is greater than 1)
    
    """
    test =  sum([((item[0] - item[1]) * (item[0] - item[1]))/item[1] \
                   for item in data.Items])
    num_rows = len(data)
    num_cols = len([col for col in data.Cols])
    if num_rows == 1:
        df = num_cols - 1
    elif num_cols == 1:
        df = num_rows - 1
    elif num_rows == 0 or num_cols == 0:
        raise ValueError, "data matrix must have data"
    else:
        df = (len(data) - 1) * (len([col for col in data.Cols]) - 1)
    
    return test, chi_high(test, df)
    

def likelihoods(d_given_h, priors):
    """Calculate likelihoods through marginalization, given Pr(D|H) and priors.

    Usage: scores = likelihoods(d_given_h, priors)

    d_given_h and priors are equal-length lists of probabilities. Returns
    a list of the same length of numbers (not probabilities).
    """
    #check that the lists of Pr(D|H_i) and priors are equal
    length = len(d_given_h)
    if length != len(priors):
        raise ValueError, "Lists not equal lengths."
    #find weighted sum of Pr(H_i) * Pr(D|H_i)
    wt_sum = 0
    for d, p in zip(d_given_h, priors):
        wt_sum += d * p
    #divide each Pr(D|H_i) by the weighted sum and multiply by its prior
    #to get its likelihood
    return [d/wt_sum for d in d_given_h]

def posteriors(likelihoods, priors):
    """Calculate posterior probabilities given priors and likelihoods.

    Usage: probabilities = posteriors(likelihoods, priors)

    likelihoods is a list of numbers. priors is a list of probabilities.
    Returns a list of probabilities (0-1).
    """
    #Check that there is a prior for each likelihood
    if len(likelihoods) != len(priors):
        raise ValueError, "Lists not equal lengths."
    #Posterior probability is defined as prior * likelihood
    return [l * p for l, p in zip(likelihoods, priors)]

def bayes_updates(ds_given_h, priors = None):
    """Successively apply lists of Pr(D|H) to get Pr(H|D) by marginalization.

    Usage: final_probs = bayes_updates(ds_given_h, [priors])

    ds_given_h is a list (for each form of evidence) of lists of probabilities.
    priors is optionally a list of the prior probabilities.
    Returns a list of posterior probabilities.
    """
    try:
        first_list = ds_given_h[0]
        length = len(first_list)
        #calculate flat prior if none was passed
        if not priors:
            priors = [1/length] * length
        #apply each form of data to the priors to get posterior probabilities
        for index, d in enumerate(ds_given_h):
            #first, ignore the form of data if all the d's are the same
            all_the_same = True
            first_element = d[0]
            for i in d:
                if i != first_element:
                    all_the_same = False
                    break
            if not all_the_same:    #probabilities won't change
                if len(d) != length:
                    raise ValueError, "bayes_updates requires equal-length lists."
                liks = likelihoods(d, priors)
                pr = posteriors(liks, priors)
                priors = pr
        return priors #posteriors after last calculation are 'priors' for next
    #return column of zeroes if anything went wrong, e.g. if the sum of one of
    #the ds_given_h is zero.
    except (ZeroDivisionError, FloatingPointError):
        return [0] * length

def t_paired (a,b, tails=None, exp_diff=0):
    """Returns t and prob for TWO RELATED samples of scores a and b.  
    
    From Sokal and Rohlf (1995), p. 354.
    Calculates the vector of differences and compares it to exp_diff
    using the 1-sample t test.

    Usage:   t, prob = t_paired(a, b, tails, exp_diff)
    
    t is a float; prob is a probability.
    a and b should be equal-length lists of paired observations (numbers).
    tails should be None (default), 'high', or 'low'.
    exp_diff should be the expected difference in means (a-b); 0 by default.
    """
    n = len(a)
    if n != len(b):
        raise ValueError, 'Unequal length lists in ttest_paired.'
    try:
        diffs = array(a) - array(b)
        return t_one_sample(diffs, popmean=exp_diff, tails=tails)
    except (ZeroDivisionError, ValueError, AttributeError, TypeError, \
        FloatingPointError):
        return (None, None)
    

def t_one_sample(a,popmean=0, tails=None):
    """Returns t for ONE group of scores a, given a population mean.

    Usage:   t, prob = t_one_sample(a, popmean, tails)

    t is a float; prob is a probability.
    a should support Mean, StandardDeviation, and Count.
    popmean should be the expected mean; 0 by default.
    tails should be None (default), 'high', or 'low'.
"""
    try:
        n = len(a)
        t = (mean(a) - popmean)/(std(a)/sqrt(n))
    except (ZeroDivisionError, ValueError, AttributeError, TypeError, \
        FloatingPointError):
        return None, None
    if isnan(t) or isinf(t):
        return None, None

    prob = t_tailed_prob(t, n-1, tails)
    return t, prob


def t_two_sample (a, b, tails=None, exp_diff=0):
    """Returns t, prob for two INDEPENDENT samples of scores a, and b.  
    
    From Sokal and Rohlf, p 223.  

    Usage:   t, prob = t_two_sample(a,b, tails, exp_diff)

    t is a float; prob is a probability.
    a and b should be sequences of observations (numbers). Need not be equal
    lengths.
    tails should be None (default), 'high', or 'low'.
    exp_diff should be the expected difference in means (a-b); 0 by default.
 """
    try:
        #see if we need to back off to the single-observation for single-item
        #groups
        n1 = len(a)
        if n1 < 2:
            return t_one_observation(sum(a), b, tails, exp_diff)
        n2 = len(b)
        if n2 < 2:
            return t_one_observation(sum(b), a, reverse_tails(tails), exp_diff)
        #otherwise, calculate things properly
        x1 = mean(a)
        x2 = mean(b)
        df = n1+n2-2
        svar = ((n1-1)*var(a) + (n2-1)*var(b))/df
        t = (x1-x2-exp_diff)/sqrt(svar*(1/n1 + 1/n2))
    except (ZeroDivisionError, ValueError, AttributeError, TypeError, \
        FloatingPointError), e:
        #bail out if the sample sizes are wrong, the values aren't numeric or
        #aren't present, etc.
        return (None, None)
    if isnan(t) or isinf(t):
        return (None, None)

    prob = t_tailed_prob(t, df, tails)
    return t, prob

def t_one_observation(x, sample, tails=None, exp_diff=0):
    """Returns t-test for significance of single observation versus a sample.
    
    Equation for 1-observation t (Sokal and Rohlf 1995 p 228):
    t = obs - mean - exp_diff / (var * sqrt((n+1)/n)) 
    df = n - 1
    """
    try:
        n = len(sample)
        t = (x - mean(sample) - exp_diff)/std(sample)/sqrt((n+1)/n)
    except (ZeroDivisionError, ValueError, AttributeError, TypeError, \
        FloatingPointError):
        return (None, None)
    prob = t_tailed_prob(t, n-1, tails)
    return t, prob

def pearson(x_items, y_items):
    """Returns Pearson correlation coefficient between x and y."""
    x_items, y_items = array(x_items), array(y_items)
    sum_x = sum(x_items)
    sum_y = sum(y_items)
    sum_x_sq = sum(x_items*x_items)
    sum_y_sq = sum(y_items*y_items)
    sum_xy = sum(x_items*y_items)
    n = len(x_items)
    try:
        r = 1.0 * ((n * sum_xy) - (sum_x * sum_y)) / \
           (sqrt((n * sum_x_sq)-(sum_x*sum_x))*sqrt((n*sum_y_sq)-(sum_y*sum_y)))
    except (ZeroDivisionError, ValueError, FloatingPointError): #no variation
        r = 0.0
    #check we didn't get a naughty value for r due to rounding error
    if r > 1.0:
        r = 1.0
    elif r < -1.0:
        r = -1.0
    return r

def correlation(x_items, y_items):
    """Returns Pearson correlation between x and y, and its significance.
    
    WARNING: x_items and y_items must be same length!
    """
    r = pearson(x_items, y_items)
    n = len(x_items)
    if n < 3:
        prob = 1
    else:
        try:
            t = r/sqrt((1 - (r*r))/(n-2))
            prob = tprob(t, n-2)
        except (ZeroDivisionError, FloatingPointError): #r was presumably 1
            prob = 0
    return (r, prob)

def correlation_matrix(series, as_rows=True):
    """Returns pairwise correlations between each pair of series.
    """
    return corrcoef(series, rowvar=as_rows)
    #unused codes below
    if as_rows:
        return corrcoef(transpose(array(series)))
    else:
        return corrcoef(array(series))

def regress(x, y):
    """Returns coefficients to the regression line "y=ax+b" from x[] and y[].  

    Specifically, returns (slope, intercept) as a tuple from the regression of 
    y on x, minimizing the error in y assuming that x is precisely known.

    Basically, it solves 
        Sxx a + Sx b = Sxy
         Sx a +  N b = Sy
    where Sxy = \sum_i x_i y_i, Sx = \sum_i x_i, and Sy = \sum_i y_i.  The
    solution is
        a = (Sxy N - Sy Sx)/det
        b = (Sxx Sy - Sx Sxy)/det
    where det = Sxx N - Sx^2.  In addition,
        Var|a| = s^2 |Sxx Sx|^-1 = s^2 | N  -Sx| / det
           |b|       |Sx  N |          |-Sx Sxx|
        s^2 = {\sum_i (y_i - \hat{y_i})^2 \over N-2}
            = {\sum_i (y_i - ax_i - b)^2 \over N-2}
            = residual / (N-2)
        R^2 = 1 - {\sum_i (y_i - \hat{y_i})^2 \over \sum_i (y_i - \mean{y})^2}
            = 1 - residual/meanerror

    Adapted from the following URL:
    http://www.python.org/topics/scicomp/recipes_in_python.html
    """
    x, y = array(x), array(y)
    N = len(x)
    Sx = sum(x)
    Sy = sum(y)
    Sxx = sum(x*x)
    Syy = sum(y*y)
    Sxy = sum(x*y)
    det = Sxx * N - Sx * Sx
    return (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det

def regress_origin(x,y):
    """Returns coefficients to regression "y=ax+b" passing through origin.

    Requires vectors x and y of same length.
    See p. 351 of Zar (1999) Biostatistical Analysis.

    returns slope, intercept as a tuple.
    """
    x, y = array(x), array(y)
    return sum(x*y)/sum(x*x), 0

def regress_R2(x, y):
    """Returns the R^2 value for the regression of x and y
    
    Used the method explained on pg 334 ofJ.H. Zar, Biostatistical analysis,
    fourth edition. 1999
    """
    slope, intercept = regress(x,y)
    coords = zip(x, y)
    Sx = Sy = Syy = SXY =  0.0
    n = float(len(y))
    for x, y in coords:
        SXY += x*y
        Sx += x
        Sy += y
        Syy += y*y
    Sxy = SXY - (Sx*Sy)/n
    regSS = slope * Sxy
    totSS = Syy - ((Sy*Sy)/n)
    return regSS/totSS

def regress_residuals(x, y):
    """reports the residual (error) for each point from the linear regression"""
    slope, intercept = regress(x, y)
    coords = zip(x, y)
    residuals = []
    for x, y in coords:
        e = y - (slope * x) - intercept
        residuals.append(e)
    return residuals
    
def stdev_from_mean(x):
    """returns num standard deviations from the mean of each val in x[]"""
    x = array(x)
    return (x - mean(x))/std(x)

def regress_major(x, y):
    """Returns major-axis regression line of y on x.

    Use in cases where there is error in both x and y.
    """
    x, y = array(x), array(y)
    N = len(x)
    Sx = sum(x)
    Sy = sum(y)
    Sxx = sum(x*x)
    Syy = sum(y*y)
    Sxy = sum(x*y)
    var_y = (Syy-((Sy*Sy)/N))/(N-1)
    var_x = (Sxx-((Sx*Sx)/N))/(N-1)
    cov = (Sxy-((Sy*Sx)/N))/(N-1)
    mean_y = Sy/N
    mean_x = Sx/N
    D = sqrt((var_y + var_x)*(var_y + var_x) - 4*(var_y*var_x - (cov*cov)))
    eigen_1 = (var_y + var_x + D)/2
    slope = cov/(eigen_1 - var_y)
    intercept = mean_y - (mean_x * slope)
    return (slope, intercept)
    

def z_test(a, popmean=0, popstdev=1, tails=None):
    """Returns z and probability score for a single sample of items.

Calculates the z-score on ONE sample of items with mean x, given a population 
mean and standard deviation (parametric).

Usage:   z, prob = z_test(a, popmean, popstdev, tails)

z is a float; prob is a probability.
a is a sample with Mean and Count.
popmean should be the parametric population mean; 0 by default.
popstdev should be the parametric population standard deviation, 1 by default.
tails should be None (default), 'high', or 'low'.
""" 
    try:
        z = (mean(a) - popmean)/popstdev*sqrt(len(a))
        return z, z_tailed_prob(z, tails)
    except (ValueError, TypeError, ZeroDivisionError, AttributeError, \
        FloatingPointError):
        return None

def z_tailed_prob(z, tails):
    """Returns appropriate p-value for given z, depending on tails."""
    if tails == 'high':
        return z_high(z)
    elif tails == 'low':
        return z_low(z)
    else:
        return zprob(z)

def t_tailed_prob(t, df, tails):
    """Return appropriate p-value for given t and df, depending on tails."""
    if tails == 'high':
        return t_high(t, df)
    elif tails == 'low':
        return t_low(t, df)
    else:
        return tprob(t,df)

def reverse_tails(tails):
    """Swaps high for low or vice versa, leaving other values alone."""
    if tails == 'high':
        return 'low'
    elif tails == 'low':
        return 'high'
    else:
        return tails

def tail(prob, test):
    """If test is true, returns prob/2. Otherwise returns 1-(prob/2).
    """
    prob /= 2
    if test:
        return prob
    else:
        return 1 - prob

def combinations(n, k):
    """Returns the number of ways of choosing k items from n.
    """
    return exp(lgam(n+1) - lgam(k+1) - lgam(n-k+1))

def multiple_comparisons(p, n):
    """Corrects P-value for n multiple comparisons.
    
    Calculates directly if p is large and n is small; resorts to logs
    otherwise to avoid rounding (1-p) to 1
    """
    if p > 1e-6:   #if p is large and n small, calculate directly
        return 1 - (1-p)**n
    else:
        return one_minus_exp(-n * p)

def multiple_inverse(p_final, n):
    """Returns p_initial for desired p_final with n multiple comparisons.
    
    WARNING: multiple_inverse is not very reliable when p_final is very close
    to 1 (say, within 1e-4) since we then take the ratio of two very similar
    numbers.
    """
    return one_minus_exp(log_one_minus(p_final)/n)

def multiple_n(p_initial, p_final):
    """Returns number of comparisons such that p_initial maps to p_final.
    
    WARNING: not very accurate when p_final is very close to 1.
    """
    return log_one_minus(p_final)/log_one_minus(p_initial)

def fisher(probs):
    """Uses Fisher's method to combine multiple tests of a hypothesis.

    -2 * SUM(ln(P)) gives chi-squared distribution with 2n degrees of freedom.
    """
    try:
        return chi_high(-2 * sum(map(log, probs)), 2 * len(probs))
    except OverflowError, e:
        return 0.0 

def f_value(a,b):
    """Returns the num df, the denom df, and the F value.

    a, b: lists of values, must have Variance attribute (recommended to
    make them Numbers objects.

    The F value is always calculated by dividing the variance of a by the 
    variance of b, because R uses the same approach. In f_two_value it's 
    decided what p-value is returned, based on the relative sizes of the
    variances. 
    """
    if not any(a) or not any(b) or len(a) <= 1 or len(b) <= 1:
        raise ValueError, "Vectors should contain more than 1 element"
    F = var(a)/var(b)
    dfn = len(a)-1
    dfd = len(b)-1
    return dfn, dfd, F


def f_two_sample(a, b, tails=None):
    """Returns the dfn, dfd, F-value and probability for two samples a, and b.

    a and b: should be independent samples of scores. Should be lists of
    observations (numbers).

    tails should be None(default, two-sided test), 'high' or 'low'.

    This implementation returns the same results as the F test in R.
    """
    dfn, dfd, F = f_value(a, b)
    if tails == 'low':
        return dfn, dfd, F, f_low(dfn, dfd, F)
    elif tails == 'high':
        return dfn, dfd, F, f_high(dfn, dfd, F)
    else:
        if var(a) >= var(b):
            side='right'
        else:
            side='left'
        return dfn, dfd, F, fprob(dfn, dfd, F, side=side)

def ANOVA_one_way(a):
    """Performs a one way analysis of variance

    a is a list of lists of observed values. Each list is the values
    within a category. The analysis must include 2 or more categories(lists).
    the lists must have a Mean and variance attribute. Recommende to make
    the Numbers objects
    
    An F value is first calculated as the variance of the group means
    divided by the mean of the within-group variances.
    """
    group_means = []
    group_variances = []
    num_cases = 0
    all_vals = []
    for i in a:
        num_cases += len(i)
        group_means.append(i.Mean)
        group_variances.append(i.Variance * (len(i)-1))
        all_vals.extend(i)
    group_means = Numbers(group_means)
    #get within group variances (denominator)
    group_variances = Numbers(group_variances)
    dfd = num_cases - len(group_means)
    within_MS = sum(group_variances)/dfd
    #get between group variances (numerator)
    grand_mean = Numbers(all_vals).Mean
    between_MS = 0
    for i in a:
        diff = i.Mean - grand_mean
        diff_sq = diff * diff
        x = diff_sq * len(i)
        between_MS += x
    dfn = len(group_means) - 1
    between_MS = between_MS/dfn
    F = between_MS/within_MS
    return dfn, dfd, F, between_MS, within_MS, group_means, f_high(dfn, dfd, F)

def MonteCarloP(value, rand_values, tail = 'high'):
    """takes a true value and a list of random values as
        input and returns a p-value

    tail indicates which side of the distribution to look at:
        low = look for smaller values than expected by chance
        high = look for larger values than expected by chance
    """
    pop_size= len(rand_values)
    rand_values.sort()
    if tail == 'high':
        num_better = pop_size
        for i, curr_val in enumerate(rand_values):
            if value <= curr_val:
                num_better = i
                break
        p_val = 1-(num_better / pop_size)
    elif tail == 'low':
        num_better = pop_size
        for i, curr_val in enumerate(rand_values):
            if value < curr_val:
                num_better = i
                break
        p_val = num_better / pop_size
    return p_val

def sign_test(success, trials, alt="two sided"):
    """Returns the probability for the sign test.
    
    Arguments:
        - success: the number of successes
        - trials: the number of trials
        - alt: the alternate hypothesis, one of 'less', 'greater', 'two sided'
          (default).
    """
    lo = ["less", "lo", "lower", "l"]
    hi = ["greater", "hi", "high", "h", "g"]
    two = ["two sided", "2", 2, "two tailed", "two"]
    alt = alt.lower().strip()
    if alt in lo:
        p = binomial_low(success, trials, 0.5)
    elif alt in hi:
        success -= 1
        p = binomial_high(success, trials, 0.5)
    elif alt in two:
        success = min(success, trials-success)
        hi = 1 - binomial_high(success, trials, 0.5)
        lo = binomial_low(success, trials, 0.5)
        p = hi+lo
    else:
        raise RuntimeError("alternate [%s] not in %s" % (lo+hi+two))
    return p

def ks_test(x, y=None, alt="two sided", exact = None, warn_for_ties = True):
    """Returns the statistic and probability from the Kolmogorov-Smirnov test.
    
    Arguments:
        - x, y: vectors of numbers whose distributions are to be compared.
        - alt: the alternative hypothesis, default is 2-sided.
        - exact: whether to compute the exact probability
        - warn_for_ties: warns when values are tied. This should left at True
          unless a monte carlo variant, like ks_boot, is being used.
    
    Note the 1-sample cases are not implemented, although their cdf's are
    implemented in ks.py"""
    # translation from R 2.4
    num_x = len(x)
    num_y = None
    x = [(x[i], 0) for i in range(num_x)]
    lo = ["less", "lo", "lower", "l", "lt"]
    hi = ["greater", "hi", "high", "h", "g", "gt"]
    two = ["two sided", "2", 2, "two tailed", "two", "two.sided"]
    Pval = None
    if y is not None: # in anticipation of actually implementing the 1-sample cases
        num_y = len(y)
        y = [(y[i], 1) for i in range(num_y)]
        n = num_x * num_y / (num_x + num_y)
        combined =  x + y
        combined.sort()
        cumsum = [0.0]
        for index, (val, origin) in enumerate(combined):
            scale = [1/num_x, -1/num_y][origin]
            cumsum.append(scale + cumsum[-1])
        cumsum.pop(0)
        if exact == None:
            exact = num_x * num_y < 1e4
        
        if len(set(combined)) < num_x + num_y:
            ties = True
        else:
            ties = False
        
        if alt in two:
            stat = max(fabs(cumsum))
        elif alt in lo:
            stat = -min(cumsum)
        elif alt in hi:
            stat = max(cumsum)
        else:
            raise RuntimeError, "Unknown alt: %s" % alt
        if exact and alt in two and not ties:
            Pval = 1 - psmirnov2x(stat, num_x, num_y)
    else:
        raise NotImplementedError
    
    if Pval == None:
        if alt in two:
            Pval = 1 - pkstwo(sqrt(n) * stat)
        else:
            Pval = exp(-2 * n * stat**2)
    
    if ties and warn_for_ties:
        warnings.warn("Cannot compute correct KS probability with ties")
    
    try: # if numpy arrays were input, the Pval can be an array of len==1
        Pval = Pval[0]
    except (TypeError, IndexError):
        pass
    return stat, Pval

def ks_boot(x, y, alt = "two sided", num_reps=1000):
    """Monte Carlo (bootstrap) variant of the Kolmogorov-Smirnov test. Useful
    for when there are ties.
    
    Arguments:
        - x, y: vectors of numbers
        - alt: alternate hypothesis, as per ks_test
        - num_reps: number of replicates for the  bootstrap"""
    # based on the ks_boot method in the R Matching package
    # see http://sekhon.berkeley.edu/matching/
    # One important difference is I preserve the original sample sizes
    # instead of making them equal
    tol = MACHEP * 100
    combined = list(x) + list(y)
    observed_stat, _p = ks_test(x, y, exact=False, warn_for_ties=False)
    total_obs = len(combined)
    num_x = len(x)
    num_greater = 0
    for i in range(num_reps):
        # sampling with replacement
        sampled = [choice(combined) for i in range(total_obs)]
        # split into the two populations
        sampled_x = sampled[:num_x]
        sampled_y = sampled[num_x:]
        sample_stat, _p = ks_test(sampled_x, sampled_y, alt=alt, exact=False,
                                    warn_for_ties=False)
        if sample_stat >= (observed_stat - tol):
            num_greater += 1
    return observed_stat, num_greater / num_reps

def permute_2d(m, p):
    """Performs 2D permutation of matrix m according to p."""
    return m[p][:, p]
    #unused below
    m_t = transpose(m)
    r_t = take(m_t, p, axis=0)
    return take(transpose(r_t), p, axis=0)

def mantel(m1, m2, n):
    """Compares two distance matrices. Reports P-value for correlation."""
    m1, m2 = asarray(m1), asarray(m2)
    m1_flat = ravel(m1)
    size = len(m1)
    orig_stat = abs(pearson(m1_flat, ravel(m2)))
    better = 0
    for i in range(n):
        #p2 = m2[permutation(size)][:, permutation(size)]
        p2 = permute_2d(m2, permutation(size))
        r = abs(pearson(m1_flat, ravel(p2)))
        if r >= orig_stat:
            better += 1
    return better/n

def kendall_correlation(x, y, alt="two sided", exact=None, warn=True):
    """returns the statistic (tau) and probability from Kendall's non-parametric
    test of association that tau==0. Uses the large sample approximation when
    len(x) >= 50 or when there are ties, otherwise it computes the probability
    exactly.
    
    Based on the algorithm implemented in R v2.5
    
    Arguments:
        - alt: the alternate hypothesis (greater, less, two sided)
        - exact: when False, forces use of the large sample approximation
          (normal distribution). Not allowed for len(x) >= 50.
        - warn: whether to warn about tied values
    """
    
    assert len(x) == len(y), "data (x, y) not of same length"
    assert len(x) > 2, "not enough observations"
    
    # possible alternate hypotheses arguments
    lo = ["less", "lo", "lower", "l", "lt"]
    hi = ["greater", "hi", "high", "h", "g", "gt"]
    two = ["two sided", "2", 2, "two tailed", "two", "two.sided", "ts"]
    
    ties = False
    num = len(x)
    ties = len(set(x)) != num or len(set(y)) != num
    if ties and warn:
        warnings.warn("Tied values, using normal approximation")
    
    if not ties and num < 50:
        exact = True
    
    if num < 50 and not ties and exact:
        combs = int(num * (num-1) / 2)
        working = []
        for i in range(combs):
            row = [-1 for j in range(combs)]
            working.append(row)
        
        tau = kendalls_tau(x, y, False)
        q = round((tau+1)*num*(num-1) / 4)
        if alt in two:
            if q > num * (num - 1) / 4:
                p = 1 - pkendall(q-1, num, Gamma(num+1), working)
            else:
                p = pkendall(q, num, Gamma(num+1), working)
            p = min(2*p, 1)
        elif alt in hi:
            p = 1 - pkendall(q-1, num, Gamma(num+1), working)
        elif alt in lo:
            p = pkendall(q, num, Gamma(num+1), working)
    else:
        tau, p = kendalls_tau(x, y, True)
        if alt in hi:
            p /= 2
        elif alt in lo:
            p = 1 - p/2
    return tau, p

## Start functions for distance_matrix_permutation_test

def distance_matrix_permutation_test(matrix, cells, cells2=None,\
        f=t_two_sample, tails=None, n=1000, return_scores=False,\
        is_symmetric=True):
    """performs a monte carlo permutation test to determine if the 
    values denoted in cells are significantly different than the rest
    of the values in the matrix

    matrix: a numpy array
    cells: a list of indices of special cells to compare to the rest of the 
        matrix
    cells2: an optional list of indices to compare cells to. If set to None
        (default), compares cells to the rest of the matrix
    f: the statistical test used. Should take a "tails" parameter as input
    tails: can be None(default), 'high', or 'low'. Input into f.
    n: the number of replicates in the Monte Carlo simulations
    is_symmetric: corrects if the matrix is symmetric. Need to only look at
        one half otherwise the degrees of freedom value will be incorrect.
    """
    #if matrix is symmetric convert all indices to lower trangular
    if is_symmetric:
        cells = get_ltm_cells(cells)
        if cells2:
            cells2 = get_ltm_cells(cells2)
    # pull out the special values
    special_values, other_values = \
        get_values_from_matrix(matrix, cells, cells2, is_symmetric)
    # calc the stat and parameteric p-value for real data
    stat, p = f(special_values, other_values, tails)
    #calc for randomized matrices
    count_more_extreme = 0
    stats = []
    indices = range(len(matrix))
    for k in range(n):
        # shuffle the order of indices, and use those to permute the matrix
        permuted_matrix = permute_2d(matrix,permutation(indices))
        special_values, other_values = \
            get_values_from_matrix(permuted_matrix, cells,\
            cells2, is_symmetric)
        # calc the stat and p for a random subset (we don't do anything 
        # with these p-values, we only use the current_stat value)
        current_stat, current_p = f(special_values, other_values, tails)
        stats.append(current_stat)
        if tails == None:
            if abs(current_stat) > abs(stat): count_more_extreme += 1
        elif tails == 'low':
            if current_stat < stat: count_more_extreme += 1
        elif tails == 'high':
            if current_stat > stat: count_more_extreme += 1

    # pack up the parametric stat, parametric p, and empirical p; calc the
    # the latter in the process
    result = [stat, p, count_more_extreme/n]
    # append the scores of the n tests if requested
    if return_scores: result.append(stats)
    return tuple(result)

def get_values_from_matrix(matrix, cells, cells2=None, is_symmetric=True):
    """get values from matrix positions in cells and cells2

        matrix: the numpy array from which values should be taken
        cells: indices of first set of requested values
        cells2: indices of second set of requested values or None
         if they should be randomly selected
        is_symmetric: True if matrix is symmetric 

    """

    # pull cells values
    cells_values = [matrix[i] for i in cells]
    # pull cells2 values
    if cells2:
        cells2_values = [matrix[i] for i in cells2]
    # or generate the indices and grab them if they
    # weren't passed in
    else:
        cells2_values = []
        for i, val_i in enumerate(matrix):
            for j, val in enumerate(val_i):
                if is_symmetric:
                    if (i,j) not in cells and i > j:
                        cells2_values.append(val)
                else:
                    if (i,j) not in cells:
                        cells2_values.append(val)
    return cells_values, cells2_values

def get_ltm_cells(cells):
    """converts matrix indices so all are below the diagonal

        cells: list of indices into a 2D integer-indexable object
         (typically a list or lists of array of arrays)
    
    """
    new_cells = []
    for cell in cells:
        if cell[0] < cell[1]:
            new_cells.append((cell[1], cell[0]))
        elif cell[0] > cell[1]:
            new_cells.append(cell)
    #remove duplicates
    new_cells = set(new_cells)
    return list(new_cells)

## End functions for distance_matrix_permutation_test
