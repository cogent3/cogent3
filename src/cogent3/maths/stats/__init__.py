#!/usr/bin/env python
"""Provides statistical tests and distributions.

Also provides NumberList and FrequencyDistribution, two classes for
working with statistical data.
"""
from .distribution import chi_high as chisqprob


__all__ = [
    "chisqprob",
    "contingency",
    "distribution",
    "information_criteria",
    "kendall",
    "ks",
    "special",
    "test",
]


__author__ = ""
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Gavin Huttley",
    "Rob Knight",
    "Sandra Smit",
    "Catherine Lozupone",
    "Micah Hamady",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"
