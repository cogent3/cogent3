#!/usr/bin/env python

from .align import (
    _align_pairwise,
    classic_align_pairwise,
    global_pairwise,
    local_pairwise,
    make_dna_scoring_dict,
)


__all__ = [
    "align",
    "dp_calculation",
    "indel_model",
    "indel_positions",
    "pairwise",
    "progressive",
    "pycompare",
    "traceback",
]

__author__ = ""
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Peter Maxwell", "Jeremy Widmann", "Gavin Huttley", "Rob Knight"]
__license__ = "BSD-3"
__version__ = "2019.9.13a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"
