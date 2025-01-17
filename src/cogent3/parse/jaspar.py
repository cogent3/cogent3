import re

from numpy import array

from cogent3.core.moltype import get_moltype
from cogent3.core.profile import MotifCountsArray

_brackets = re.compile(r"\[|\]")


def read(filepath):
    """returns matrixid and MotifCountsArray matrix"""
    with open(filepath) as infile:
        matrix = []
        states = []
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                identifier = line[1:].split()
            elif line:
                line = _brackets.sub("", line)
                line = line.split()
                states.append(line.pop(0).upper())
                matrix.append([int(i) for i in line])

    matrix = dict(zip(states, matrix, strict=False))
    name = ("rna" if "U" in states else "dna") if len(states) == 4 else "protein"

    states = list(get_moltype(name))
    matrix = array([matrix[s] for s in states], dtype=int).T

    pwm = MotifCountsArray(matrix, states)
    return identifier, pwm
