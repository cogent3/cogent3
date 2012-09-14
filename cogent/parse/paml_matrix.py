#!/usr/bin/env python

import numpy

Float = numpy.core.numerictypes.sctype2char(float)

from cogent.evolve import substitution_model

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Matthew Wakefield", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production"

three_letter_order = 'ARNDCQEGHILKMFPSTWYV'
aa_order = numpy.array([ord(aa) for aa in three_letter_order])
reorder = numpy.argsort(aa_order)

def numbers_in(f):
    for line in f:
        for word in line.split():
            yield float(word)

def PamlMatrixParser(f):
    """Parses a matrix of amino acid transition probabilities and amino acid
    frequencies in the format used by PAML and returns a symetric array in single
    letter alphabetical order and a dictionary of frequencies for use by
    substitution_model.EmpiricalProteinMatrix"""
    matrix = numpy.zeros([20,20], Float)
    next_number = numbers_in(f).next
    for row in range(1,20):
        for col in range(0, row):
            matrix[row,col] = matrix[col,row] = next_number()
            
    freqs = [next_number() for i in range(20)]
    total = sum(freqs)
    assert abs(total-1) < 0.001, freqs
    freqs = [freq/total for freq in freqs]
    
    matrix = numpy.take(matrix, reorder, 0)
    matrix = numpy.take(matrix, reorder, 1)
    
    assert numpy.alltrue(matrix == numpy.transpose(matrix))
    
    freqs = dict(zip(three_letter_order, freqs))
    
    return (matrix, freqs)
