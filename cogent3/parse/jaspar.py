import re

from numpy import array

from cogent3.core.profile import MotifCountsArray


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.08.06a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


_brackets = re.compile(r"\[|\]")


def read(filepath):
    """returns matrixid and MotifCountsArray matrix"""
    with open(filepath) as infile:
        matrix = []
        bases = []
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                identifier = line[1:].split()
            elif line:
                line = _brackets.sub("", line)
                line = line.split()
                bases.append(line.pop(0))
                matrix.append([int(i) for i in line])
    matrix = array(matrix, dtype=int).T
    pwm = MotifCountsArray(matrix, bases)
    return identifier, pwm
