Contributors

- khiron

Deprecations

- `cogent3.math.stats.distribution.binomial_exact` has been deprecated in 
favour of scipy function `scipy.stats.binom.pmf`.  `binomial_exact` will be 
discontinued in cogent3 release 2025.9.  Until then `binomial_exact` will use 
`scipy.stats.binom.pmf` internally and floating point values for `successes` 
and `trials` will be truncated to integers using `math.floor` introducing a 
potential breaking change in behaviour. 
