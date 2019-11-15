#!/usr/bin/env python
"""Provides statistical tests and distributions.

Also provides NumberList and FrequencyDistribution, two classes for
working with statistical data.
"""
__all__ = [
    "contingency",
    "distribution",
    "information_criteria",
    "kendall",
    "ks",
    "special",
    "test",
]

from .distribution import chi_high as chisqprob


__author__ = ""
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = [
    "Gavin Huttley",
    "Rob Knight",
    "Sandra Smit",
    "Catherine Lozupone",
    "Micah Hamady",
]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"
