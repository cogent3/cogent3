*****************************
Standard statistical analyses
*****************************

.. authors Tom Elliott, Gavin Huttley, Anuj Pahwa

..
    following is just a list of the filenames that need to be deleted, to be
    appended to after each one is called. Readers don't really need to see
    this housekeeping so I'm 'hiding' this code.

.. doctest::
    :hide:

    >>> filenames_to_delete = []

Random Numbers
==============

Many of the code snippets in this section use random numbers. These can be obtained using functions from Python's ``random`` module, or using ``numpy.random``. To facilitate testing, the examples "seed" the random number generator, which ensures the same results each time the code is run.

.. doctest::

    >>> import random
    >>> random.seed(157)
    >>> random.choice((-1,1))
    1
    >>> random.choice(range(1000))
    224
    >>> random.gauss(mu=50, sigma=3)
    52.7668...
    >>> import numpy as np
    >>> np.random.seed(157)
    >>> np.random.random_integers(-1,1,5)
    array([-1,  1,  1,  1,  0])

For the last example, note that the range includes 0.

.. doctest::

    >>> np.random.normal(loc=50,scale=3,size=2)
    array([ 42.8217253 ,  49.90008293])
    >>> np.random.randn(3)
    array([ 1.26613052,  1.59533412,  0.95612413])

Summary statistics
==================

Population mean and median
--------------------------

PyCogent's functions for statistical analysis operate on ``numpy`` arrays

.. doctest::

    >>> import random
    >>> import numpy as np
    >>> import cogent.maths.stats.test as stats
    >>> random.seed(157)
    >>> nums = [random.gauss(mu=50, sigma=3) for i in range(1000)]
    >>> arr = np.array(nums)
    >>> stats.mean(arr)
    49.9927...

but in some cases they will also accept a simple list of values

.. doctest::

    >>> stats.mean(range(1,8))
    4.0
    >>> stats.var(range(1,8))
    4.66666...

The keyword argument ``axis`` controls whether a function operates by rows (``axis=0``) or by columns (``axis=1``), or on all of the values (``axis=None``)

.. doctest::

    >>> import cogent.maths.stats.test as stats
    >>> import numpy as np
    >>> nums = range(1,6) + [50] + range(10,60,10) + [500]
    >>> arr = np.array(nums)
    >>> arr.shape = (2,6)
    >>> arr
    array([[  1,   2,   3,   4,   5,  50],
           [ 10,  20,  30,  40,  50, 500]])
    >>> stats.mean(arr, axis=0)
    array([   5.5,   11. ,   16.5,   22. ,   27.5,  275. ])
    >>> stats.mean(arr, axis=1)
    array([  10.83333333,  108.33333333])
    >>> stats.mean(arr)
    59.58333...
    >>> stats.median(arr, axis=0)
    array([   5.5,   11. ,   16.5,   22. ,   27.5,  275. ])
    >>> stats.median(arr, axis=1)
    array([  3.5,  35. ])
    >>> stats.median(arr)
    15.0

Population variance and standard deviation
------------------------------------------

.. doctest::

    >>> print stats.var(arr, axis=0)
    [  4.05000000e+01   1.62000000e+02   3.64500000e+02   6.48000000e+02
       1.01250000e+03   1.01250000e+05]
    >>> print stats.std(arr, axis=0)
    [   6.36396103   12.72792206   19.09188309   25.45584412   31.81980515
      318.19805153]
    >>> print stats.var(arr, axis=1)
    [   370.16666667  37016.66666667]
    >>> print stats.std(arr, axis=1)
    [  19.23971587  192.39715868]
    >>> print stats.var(arr, axis=None)
    19586.6287879
    >>> print stats.std(arr, axis=None)
    139.952237524

The variance (and standard deviation) are unbiased

.. doctest::

    >>> import numpy as np
    >>> import cogent.maths.stats.test as stats
    >>> arr = np.array([1,2,3,4,5])
    >>> m = np.mean(arr)
    >>> stats.var(arr)
    2.5
    >>> 1.0 * sum([(n-m)**2 for n in arr]) / (len(arr) - 1)
    2.5

Distributions
=============

Binomial
--------

The binomial distribution can be used for calculating the probability of specific frequencies of states occurring in discrete data. The two alternate states are typically referred to as a success or failure. This distribution is used for sign tests.

.. doctest::

    >>> import cogent.maths.stats.distribution as distr
    >>> distr.binomial_low(successes=5, trials=10, prob=0.5)
    0.623...
    >>> distr.binomial_high(successes=5, trials=10, prob=0.5)
    0.376...
    >>> distr.binomial_exact(successes=5, trials=10, prob=0.5)
    0.246...

Chi-square
----------

A convenience function for computing the probability of a chi-square statistic is provided at the ``stats`` top level.

.. doctest::

    >>> from cogent.maths.stats import chisqprob
    >>> chisqprob(3.84, 1)
    0.05...

which is just a reference to the ``chi_high`` function.

.. doctest::

    >>> from cogent.maths.stats.distribution import chi_high
    >>> chi_high(3.84, 1)
    0.05...

Getting the inverse
^^^^^^^^^^^^^^^^^^^

Given a probability we can determine the corresponding chi-square value for a given degrees-of-freedom.

.. doctest::

    >>> from cogent.maths.stats.distribution import chdtri
    >>> chdtri(1, 0.05)
    3.84...
    >>> chdtri(2, 0.05)
    5.99...

Normal
------

The function ``zprob()`` takes a z-score or standard deviation and computes the fraction of the normal distribution (mean=0, std=1) which lies farther away from the mean than that value.  For example, only about 4.5% of the values are more than 2 standard deviations away from the mean, so that more than 95% of the values are at least that close to the mean.

.. doctest::

    >>> import cogent.maths.stats.distribution as distr
    >>> for z in range(5):
    ...     print '%s %.4f' % (z, distr.zprob(z))
    ...
    0 1.0000
    1 0.3173
    2 0.0455
    3 0.0027
    4 0.0001

Use the functions ``z_low()`` and ``z_high()`` to compute the normal distribution in a directional fashion.  Here we see that a z-score of 1.65 has a value greater than 95% of the values in the distribution, and similarly a z-score of 1.96 has a value greater than 97.5% of the values in the distribution.

.. doctest::

    >>> z = 0
    >>> while distr.z_low(z) < 0.95:
    ...     z += 0.01
    ...
    >>> z
    1.6500...
    >>> z = 0
    >>> while distr.z_low(z) < 0.975:
    ...     z += 0.01
    ...
    >>> z
    1.9600...

Normalizing data (as Z-scores)
==============================

The function ``z_test()`` takes a sample of values as the first argument, and named arguments for the population parameters:  ``popmean`` and ``popstddev`` (with default  values of 0 and 1), and returns the z-score and its probability.

In this example, we grab a sample from a population with ``mean=50`` and ``std=3``, and call ``z_test()`` with the population mean specified as 50 and the ``popstddev`` assuming its default value of 1:

.. Uses the parametric standard deviation

.. doctest::

    >>> import numpy as np
    >>> import cogent.maths.stats.test as stats
    >>> np.random.seed(157)
    >>> arr = np.random.normal(loc=50,scale=3,size=1000)
    >>> round(stats.mean(arr), 1)
    49.8...
    >>> round(stats.std(arr), 1)
    3.1...
    >>> z, prob = stats.z_test(arr, popmean=50.0)
    >>> print z
    -3.08...

.. todo

    TE:  I think the above needs more explanation.  What does this have to do with a Z-score, as in Z = (arr - stats.mean(arr))/stats.std(arr)?

Resampling based statistics
===========================

The Jackknife
-------------

This is a data resampling based approach to estimating the confidence in measures of location (like the mean). The method is based on omission of one member of a sample and recomputing the statistic of interest. This measures the influence of individual observations on the sample and also the confidence in the statistic.

The ``Jackknife`` class relies on our ability to handle a set of indexes for sub-setting our data and re-computing our statistic. The client code must be able to take a indices and generate a new statistic.

We demo using the jackknife the estimate of mean GC% for an alignment. We first write a factory function to compute the confidence in the mean GC% for an alignment by sampling specific columns.

.. doctest::
    
    >>> def CalcGc(aln):
    ...     def calc_gc(indices):
    ...         new = aln.takePositions(indices)
    ...         probs = new.getMotifProbs()
    ...         gc = sum(probs[b] for b in 'CG')
    ...         total = sum(probs[b] for b in 'ACGT')
    ...         return gc / total
    ...     return calc_gc

We then create an instance of this factory function with a specific alignment.

.. doctest::
    
    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs('data/test.paml', moltype=DNA)
    >>> calc_gc = CalcGc(aln)

We now create a ``Jackknife`` instance, passing it the ``calc_gc`` instance we have just made and obtain the sampling statistics. We specify how many elements we're interested in (in this case, the positions in the alignment).

.. doctest::
    
    >>> from cogent.maths.stats.jackknife import JackknifeStats
    >>> jk = JackknifeStats(len(aln), calc_gc)
    >>> print jk.SampleStat
    0.4766...
    >>> print jk.SummaryStats
    Summary Statistics
    ===============================================
    Sample Stat    Jackknife Stat    Standard Error
    -----------------------------------------------
         0.4767            0.4767            0.0584
    -----------------------------------------------

We also display the sub-sample statistics.

.. doctest::
    
    >>> print jk.SubSampleStats
    Subsample Stats
    ============
     i    Stat-i
    ------------
     0    0.4678
     1    0.4678
     2    0.4847
     3    0.4814...

.. note:: You can provide a custom index generation function that omits groups of observations, for instance. This can be assigned to the ``gen_index`` argument of the ``Jackknife`` constructor.

Permutations
============

.. this is really a numpy features

Random
------

.. doctest::

    >>> from numpy.random import permutation as perm
    >>> import numpy as np
    >>> np.random.seed(153)
    >>> arr = np.array(range(5))
    >>> for i in range(3):
    ...     print perm(arr)
    ...
    [2 1 3 0 4]
    [0 3 2 4 1]
    [4 0 1 2 3]

Ordered
-------

*To be written.*

Differences in means
====================

Consider a single sample of 50 value:

.. doctest::

    >>> import numpy as np
    >>> import cogent.maths.stats.test as stats
    >>> np.random.seed(1357)
    >>> nums1 = np.random.normal(loc=45,scale=10,size=50)

Although we don't know the population values for the mean and standard deviation for this sample, we can evaluate the probability that the sample could have been drawn from some population with known values, as shown above in Normalizing data (as Z-scores).

If we have a second sample, whose parent population mean and standard deviation are also unknown:

.. doctest::

    >>> nums2 = np.random.normal(loc=50,scale=10,size=50)

Suppose we believe (before we see any data) that the mean of the first population is different than the second but we don't know in which direction the change lies, we estimate the standard deviation.  We use the standard error of the mean as an estimate for how close the mean of sample 2 is to the mean of its parent population (and vice-versa).

.. doctest::

    >>> mean_nums2 = stats.mean(nums2)
    >>> sd_nums2 = stats.std(nums2)
    >>> se_nums2 = sd_nums2 / np.sqrt(len(nums2))
    >>> se_nums2
    1.1113...
    >>> mean_nums1 = stats.mean(nums1)
    >>> mean_nums1
    46.5727...
    >>> mean_nums2
    50.3825...
    >>> mean_nums1 < mean_nums2 - 1.96 * se_nums2
    True

t-Tests
=======

Small sample sizes can be handled by the use of t-tests.  The function ``t_two_sample()`` is used for two independent samples.

.. doctest::

    >>> subsample1 = nums1[:5]
    >>> [str(round(n,2)) for n in subsample1]
    ['49.25', '38.87', '47.06', '44.49', '43.73']
    >>> stats.mean(subsample1)
    44.67901...
    >>> subsample2 = nums2[:5]
    >>> [str(round(n,2)) for n in subsample2]
    ['51.57', '40.6', '49.62', '46.69', '59.34']
    >>> stats.mean(subsample2)
    49.56494...
    >>> t, prob = stats.t_two_sample(subsample1,subsample2)
    >>> t
    -1.3835...
    >>> prob
    0.20388...

The two sample means are not significantly different.

If there is one small sample and we want to ask whether it is unlikely to have come from a population with a known mean, use the function ``t_one_sample()``

.. doctest::

    >>> import cogent.maths.stats.test as stats
    >>> arr = [52.6, 51.3, 49.8]
    >>> t, prob = stats.t_one_sample(arr, popmean=48, tails='high')
    >>> t
    3.99681...
    >>> prob
    0.02863...

For related samples (pre- and post-treatment), use the function ``t_paired()``

.. doctest::

    >>> import cogent.maths.stats.test as stats
    >>> pre =  [52.6, 51.3, 49.8]
    >>> post = [62.6, 75.0, 65.2]
    >>> t, prob = stats.t_paired(pre, post, tails='low')
    >>> t
    -4.10781...
    >>> prob
    0.02723...

Sign test
=========

This is essentially just a test using the binomial distribution where the probability of success = 0.5.

.. doctest::

    >>> from cogent.maths.stats.test import sign_test
    >>> sign_test(40, 100)
    0.056...

Differences in proportions
==========================

*To be written.*

Association
===========

We create some data for testing for association.

.. doctest::

    >>> import numpy as np
    >>> np.random.seed(13)
    >>> x_nums = range(1,11)
    >>> error = [1.5 * random.random() for i in range(len(x_nums))]
    >>> error = [e * random.choice((-1,1)) for e in error]
    >>> y_nums = [(x * 0.5) + e for x, e in zip(x_nums, error)]
    >>> x_array = np.array(x_nums)
    >>> y_array = np.array(y_nums)

We then compute Kendall's tau and associated probability, which tests the null hypothesis that x and y are not associated.

.. doctest::

    >>> from cogent.maths.stats.test import kendall_correlation
    >>> tau, prob = kendall_correlation(x_array, y_array)
    >>> print tau
    0.688...
    >>> print prob
    0.00468...

Correlation
===========

For this example, we generate y-values as one-half the x-value plus a bit of random error

.. doctest::

    >>> import numpy as np
    >>> np.random.seed(13)
    >>> x_array = np.arange(1,11)
    >>> error = np.random.normal(size=10)
    >>> y_array = x_array * 0.5 + error
    >>> x_array
    array([ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10])
    >>> [str(round(n,2)) for n in y_array]
    ['-0.21', '1.75', '1.46', '2.45', '3.85', '3.53', '4.85', '4.86', '5.98', '3.95']

The function ``correlation()`` returns the Pearson correlation between x and y, as well as its significance

.. doctest::

    >>> import cogent.maths.stats.test as stats
    >>> r, prob = stats.correlation(x_array, y_array)
    >>> r
    0.8907...
    >>> prob
    0.0005...

The function ``regress()`` returns the coefficients to the regression line "y=ax+b"

.. doctest::

    >>> slope, y_intercept = stats.regress(x_array, y_array)
    >>> slope
    0.5514...
    >>> y_intercept
    0.2141...

Calculate the R^2 value for the regression of x and y

.. doctest::

    >>> R_squared = stats.regress_R2(x_array, y_array)
    >>> R_squared
    0.7934...

And finally, the residual error for each point from the linear regression

.. doctest::

    >>> error = stats.regress_residuals(x_array, y_array)
    >>> error = [str(round(e,2)) for e in error]
    >>> error
    ['-0.98', '0.44', '-0.41'...

Differences in variances
========================

*To be written.*

Chi-Squared test
================

.. TODO pick a biological example, perhaps sequence nucleotide composition?  Codon usage for a particular amino acid?

Calculus class data (from Grinstead and Snell, Introduction to Probability).  There seems to be a disparity in the number of 'A' grades awarded when broken down by student gender.  As input to the function ``chi_square_from_Dict2D()`` we need a ``Dict2D`` object containing the observed counts that has been processed by ``calc_contingency_expected()`` to add the expected counts for each element of the table

``Expected = row_total x column_total / overall_total``

.. doctest::

    >>> from cogent.util.dict2d import Dict2D
    >>> import cogent.maths.stats.test as stats
    >>> F_grades = {'A':37,'B':63,'C':47,'F':5}
    >>> M_grades = {'A':56,'B':60,'C':43,'F':8}
    >>> grades = {'F':F_grades,'M':M_grades}
    >>> data = Dict2D(grades)
    >>> data
    {'M': {'A': 56...
    >>> OE_data = stats.calc_contingency_expected(data)
    >>> OE_data
    {'M': {'A': [56, 48.686...
    >>> test, chi_high = stats.chi_square_from_Dict2D(OE_data)
    >>> test
    4.12877...
    >>> chi_high
    0.24789...

Nearly 25% of the time we would expect a Chi-squared statistic as extreme as this one or more (with df = 3), so the result is not significant.

Goodness-of-fit calculation with the same data

.. doctest::

    >>> g_val, prob = stats.G_fit_from_Dict2D(OE_data)
    >>> g_val
    4.1337592429166437
    >>> prob
    0.76424854978813872

Scatterplots
============

In this example, we generate the error as above, but separately from the x-value, and subsequently transform using matrix multiplication

.. doctest::

    >>> import random
    >>> import numpy as np
    >>> import cogent.maths.stats.test as stats
    >>> random.seed(13)
    >>> x_nums = range(1,11)
    >>> error = [1.5 * random.random() for i in range(len(x_nums))]
    >>> error = [e * random.choice((-1,1)) for e in error]
    >>> arr = np.array(x_nums + error)
    >>> arr.shape = (2, len(x_nums))
    >>> arr
    array([[  1.        ,   2.        ,   3.        ,   4.        ,
              5.        ,   6.        ,   7.        ,   8.        ,
              9.        ,  10.        ],
           [  0.38851274,  -1.02788699,  -1.02612288,  -1.27400424,
              0.27858626,   0.34583791,  -0.22073988,  -0.3377444 ,
             -1.1010354 ,   0.19531953]])

We can use a transformation matrix to rotate the points

.. doctest::

    >>> from math import sqrt
    >>> z = 1.0/sqrt(2)
    >>> t = np.array([[z,-z],[z,z]])
    >>> rotated_x, rotated_y = np.dot(t,arr)

The plotting code uses matplotlib_.

.. doctest::

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.scatter(arr[0],arr[1],s=250,color='b',marker='s')
    <matplotlib.collections.RegularPolyCollection object...
    >>> ax.scatter(rotated_x,rotated_y,s=250,color='r',marker='o')
    <matplotlib.collections.CircleCollection object...
    >>> plt.axis('equal')
    (0.0, 12.0, -2.0, 8.0)

Plot the least squares regression lines too

.. doctest::

    >>> slope, y_intercept = stats.regress(rotated_x, rotated_y)
    >>> slope
    0.9547989732316251
    >>> max_x = 10
    >>> ax.plot([0, max_x],[y_intercept, max_x * slope + y_intercept],
    ...     linewidth=4.0, color='k')
    ...
    [<matplotlib.lines.Line2D object...
    >>> slope, y_intercept = stats.regress(arr[0],arr[1])
    >>> ax.plot([0, max_x],[y_intercept, max_x * slope + y_intercept],
    ...     linewidth=4.0, color='0.6')
    ...
    [<matplotlib.lines.Line2D object...
    >>> plt.grid(True)
    >>> plt.savefig('scatter_example.pdf')

(If you want to plot the lines under the points, specify ``zorder=n`` to the plot commands, where ``zorder`` for the lines < ``zorder`` for the points).

..
    Possibly split these out into "Visualizing data"

.. doctest::
    :hide:

    >>> filenames_to_delete.append('scatter_example.pdf')

Histograms
==========

.. doctest::
    :hide:

    >>> plt.clf() # because the plot gets screwed up by operations above

.. doctest::

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> plt.clf()
    >>> mu, sigma = 100, 15
    >>> x = mu + sigma*np.random.randn(10000)
    >>> n, bins, patches = plt.hist(x, 60, normed=1, facecolor='0.75')

add a "best fit" line

.. doctest::

    >>> import matplotlib.mlab as mlab
    >>> y = mlab.normpdf( bins, mu, sigma)
    >>> l = plt.plot(bins, y, 'r--', linewidth=3)
    >>> plt.grid(True)
    >>> plt.savefig('hist_example.png')

.. doctest::
    :hide:

    >>> filenames_to_delete.append('hist_example.png')

Heat Maps
=========

Representing numbers as colors is a powerful data visualization technique.  This example does not actually use any functionality from PyCogent, it simply highlights a convenient matplotlib_ method for constructing a heat map.

.. doctest::

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> data = [i * 0.01 for i in range(100)]
    >>> data = np.array(data)
    >>> data.shape = (10,10)

The plot code

.. doctest::

    >>> fig = plt.figure()
    >>> plt.hot()
    >>> plt.pcolormesh(data)
    <matplotlib.collections.QuadMesh object ...
    >>> plt.colorbar()
    <matplotlib.colorbar.Colorbar instance ...
    >>> plt.savefig('heatmap_example.png')

.. doctest::
    :hide:

    >>> filenames_to_delete.append('heatmap_example.png')
    >>> from cogent.util.misc import remove_files
    >>> remove_files(filenames_to_delete, error_on_missing=False)

.. _matplotlib: http://matplotlib.sourceforge.net/
