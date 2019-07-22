import numpy


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.07.10a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


def read(filepath):
    """returns base order, weights matrix"""
    try:
        infile = open(filepath)
        data = infile.readlines()
        infile.close()
    except TypeError:
        data = filepath

    data = [l.split() for l in data]
    revised = list(zip(*data))
    bases = []
    matrix = []
    for row in revised[1:]:
        bases.append(row[0])
        matrix.append([float(i) for i in row[1:]])
    matrix = matrix
    return bases, matrix


if __name__ == "__main__":
    data = """Pos A   C   G   T
1   0.199862209150251   0.244112253499856   0.184496072294034   0.37152946505586
2   0.169704432776104   0.189437336764108   0.0940717974375604  0.546786433022228
3   0.525652161683299   0.067669780378447   0.238564497740711   0.168113560197543
4   0.783953330692583   0.088619284484793   0.0630765767243543  0.0643508080982694
5   0.112432463539848   0.147282530338847   0.0997186496659685  0.640566356455337
6   0.162189909683132   0.283941414680916   0.275549677667935   0.278318997968017
7   0.259011956338907   0.0787969447816471  0.481484588761462   0.180706510117985
8   0.184004977170827   0.285705856273927   0.271771817242226   0.258517349313019
9   0.186869253293474   0.284407819419656   0.181346432217115   0.347376495069755
""".splitlines()
    data = [l.strip() for l in data if l.strip()]
    b, m = read(data)
    print(b[0], m[0])
