****************
Useful Utilities
****************

Using PyCogent's optimisers for your own functions
==================================================

You have a function that you want to maximise/minimise. The parameters in your function may be bounded (must lie in a specific interval) or not. The cogent optimisers can be applied to these cases. The ``Powell`` (a local optimiser) and ``SimulatedAnnealing`` (a global optimiser) classes in particular have had their interfaces standardised for such use cases. We demonstrate for a very simple function below.

We write a simple factory function that uses a provided value for omega to compute the squared deviation from an estimate.

.. doctest::
    
    >>> import numpy
    >>> def DiffOmega(omega):
    ...     def omega_from_S(S):
    ...         omega_est = S/(1-numpy.e**(-1*S))
    ...         return abs(omega-omega_est)**2
    ...     return omega_from_S

We then import the ``Powell`` optimiser, create our optimisable function with a value of omega to solve for S and instantiate an optimiser instance. Note that we provide lower and upper bounds and an initial guess for our parameter of interest (``S``).

.. doctest::
    
    >>> from cogent.maths.optimisers import Powell
    >>> omega = 0.1
    >>> opt = Powell(DiffOmega(omega),
    ...     xinit=[1.0], # the vector of initial values
    ...     bounds=([-100], [100]), # [(lower),(upper)] bounds for the params
    ...     direction=-1) # -1 is minimise func, 1 is maximise

We then optimise the function, obtaining the fit statistic and the associated estimate of S.

.. doctest::
    
    >>> fit, vec = opt.run(show_progress=False)
    >>> assert 0.0 <= fit < 1e-6
    >>> print 'S=%.4f' % vec[0]
    S=-3.6150

Cartesian products
==================

*To be written.*

.. cogent.util.transform

