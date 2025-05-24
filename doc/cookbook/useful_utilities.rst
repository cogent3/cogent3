****************
Useful Utilities
****************

.. authors, Gavin Huttley

.. include:: union_dict.rst

Using Cogent3's optimisers for your own functions
===================================================

You have a function that you want to maximise/minimise. The parameters in your function may be bounded (must lie in a specific interval) or not. The ``cogent3`` optimisers can be applied to these cases. The ``Powell`` (a local optimiser) and ``SimulatedAnnealing`` (a global optimiser) classes in particular have had their interfaces standardised for such use cases. We demonstrate for a very simple function below.

We write a simple factory function that uses a provided value for omega to compute the squared deviation from an estimate, then use it to create our optimisable function.

.. jupyter-execute::

    import numpy

    def DiffOmega(omega):
        def omega_from_S(S):
            omega_est = S / (1 - numpy.e ** (-1 * S))
            return abs(omega - omega_est) ** 2

        return omega_from_S

    omega = 0.1
    f = DiffOmega(omega)

We then import the minimise function and use it to minimise the function, obtaining the fit statistic and the associated estimate of S. Note that we provide lower and upper bounds (which are optional) and an initial guess for our parameter of interest (``S``).

.. jupyter-execute::

    from cogent3.maths.optimisers import maximise, minimise

    S = minimise(
        f,  # the function
        xinit=1.0,  # the initial value
        bounds=(-100, 100),  # [lower,upper] bounds for the parameter
        local=True,  # just local optimisation, not Simulated Annealing
        show_progress=False,
    )
    assert 0.0 <= f(S) < 1e-6
    print("S=%.4f" % S)

The minimise and maximise functions can also handle multidimensional optimisations, just make xinit (and the bounds) lists rather than scalar values.
