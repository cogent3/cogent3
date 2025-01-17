from numpy import array

from cogent3.core.moltype import get_moltype
from cogent3.core.profile import MotifFreqsArray


def read(filepath):
    """returns MotifFreqsArray matrix"""
    try:
        infile = open(filepath)
        data = infile.readlines()
        infile.close()
    except TypeError:
        data = filepath

    data = [l.split() for l in data]
    revised = list(zip(*data, strict=False))
    states = []
    matrix = []
    for row in revised[1:]:
        states.append(row[0])
        matrix.append([float(i) for i in row[1:]])

    matrix = dict(zip(states, matrix, strict=False))
    name = ("rna" if "U" in states else "dna") if len(states) == 4 else "protein"

    states = list(get_moltype(name))
    matrix = [matrix[s] for s in states]
    matrix = array(matrix, dtype=float)

    return MotifFreqsArray(matrix.T, states)
