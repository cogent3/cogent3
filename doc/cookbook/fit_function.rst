****************************************************
Fitting a function to a giving set of x and y values
****************************************************

.. authors, Antonio Gonzalez Pena, Rob Knight 

Giving a set of values for ``x`` and ``y`` fit a function func that has 
``n_params`` using simplex iterations to minimize the error between the 
model ``func`` to fit and the given values.


Fitting an exponential function
===============================

.. doctest::

    >>> from numpy import array, arange, exp
    >>> from numpy.random import rand
    >>> from cogent.maths.fit_function import fit_function
    >>> #from whatever import fit_function
    >>> 
    >>> # creating x values
    >>> x = arange(-1,1,.01)
    >>> 
    >>> # defining our fitting function
    >>> def f(x,a):
    >>>     return exp(a[0]+x*a[1])
    >>> 
    >>> # getting our real y
    >>> y = f(x,a=[2,5])
    >>> 
    >>> # creating our noisy y
    >>> y_noise = y + rand(len(y))*5
    >>> 
    >>> # fitting our noisy data to the function using 1 iteration
    >>> params = fit_function(x, y_noise, f, 2, 1)
    >>> params
    array([ 2.04676934,  4.95324591])
    >>> 
    >>> # fitting our noisy data to the function using 1 iteration
    >>> params = fit_function(x, y_noise, f, 2, 5)
    >>> params
    array([ 2.04676659,  4.95324502])

