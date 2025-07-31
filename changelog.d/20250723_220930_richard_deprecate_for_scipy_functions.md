### Contributors

- @khiron

### Deprecations

- Modules `cogent3.maths.stats.special.py` and `cogent3.maths.stats.distributions.py`
have been deprecated in favour of `scipy` and `numpy` implementations.
The deprecated modules will be removed in release 2025.10.

### Enhancements

- Moved `theoretical_quantiles` and `probability_points` functions
into `cogent3.maths.stats.test.py`. Now using scipy functions for calculations.
Use can specify additional arguments to the scipy functions
as keyword arguments to these functions, eg. `theoretical_quantiles(data, "t", df=2)`.
