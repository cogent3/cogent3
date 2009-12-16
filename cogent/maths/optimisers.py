#!/usr/bin/env python

"""
wrapper to include various optimiser classes in one include file
"""
from simannealingoptimiser import SimulatedAnnealing
from scipy_optimisers import Powell, DownhillSimplex

__author__ = "Andrew Butterfield"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Andrew Butterfield", "Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"
