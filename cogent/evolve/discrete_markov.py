import numpy

from cogent.util.warning import deprecated
from cogent.recalculation.definition import (NonParamDefn, CalcDefn, 
    EvaluatedCell, PartitionDefn, ConstCell, ConstDefn,
    DictArrayTemplate)
    
__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

class PsubMatrixDefn(PartitionDefn):
    "Square 2D array made of 1D partitions"
    
    numeric = False # well, not scalar anyway
    const_by_default = False
    independent_by_default = True
    
    def __init__(self, default=None, name=None, dimensions=None,
            dimension=None, size=None, **kw):
        PartitionDefn.__init__(self, default, name, dimensions,
            dimension, size, **kw)
        
        (dim_name, dim_cats) = self.internal_dimension
        self.internal_dimensions = (dim_name, dim_name+"2")
        self.array_template = DictArrayTemplate(dim_cats, dim_cats)
    
    def _makeDefaultValue(self):
        # Purely flat default doesn't work well so start at approx t=0.5
        flat = numpy.ones([self.size, self.size], float) / self.size
        diag = numpy.identity(self.size, float)
        return (flat + diag) / 2

    def checkValueIsValid(self, value, is_constant):
        if value.shape != (self.size,self.size):
            raise ValueError("Wrong array shape %s for %s, expected (%s,%s)" % 
                    (value.shape, self.name, self.size, self.size))
        for part in value:
            PartitionDefn.checkValueIsValid(self, part, is_constant)
    
    def makeCells(self, input_soup={}, variable=None):
        uniq_cells = []
        all_cells = []
        for (i, v) in enumerate(self.uniq):
            if v is None:
                raise ValueError("input %s not set" % self.name)
            assert hasattr(v, 'getDefaultValue'), v
            value = v.getDefaultValue()
            assert hasattr(value, 'shape'), value
            assert value.shape == (self.size,self.size)
            scope = [key for key in self.assignments
                    if self.assignments[key] is v]
            if v.is_constant or (variable is not None and variable is not v):
                matrix = ConstCell(self.name, value)
            else:
                rows = []
                for part in value:
                    (ratios, partition) = self._makePartitionCell(
                            self.name+'_part', scope, part)
                    all_cells.extend(ratios)
                    rows.append(partition)
                all_cells.extend(rows)
                matrix = EvaluatedCell(self.name, lambda *x:numpy.array(x), 
                        rows)
            all_cells.append(matrix)
            uniq_cells.append(matrix)
        return (all_cells, uniq_cells)

def DiscreteSubstitutionModel(*args, **kw):
    deprecated("class",
        "cogent.evolve.discrete_markov.DiscreteSubstitutionModel",
        "cogent.evolve.substitution_model.DiscreteSubstitutionModel",
        '1.6')
    from cogent.evolve.substitution_model import DiscreteSubstitutionModel
    return DiscreteSubstitutionModel(*args, **kw)

class PartialyDiscretePsubsDefn(object):
    def __init__(self, alphabet, psubs, discrete_edges):
        motifs = tuple(alphabet)
        dpsubs = PsubMatrixDefn(
            name="dpsubs", dimension = ('motif', motifs), default=None, 
            dimensions=('locus', 'edge'))
        self.choices = [psubs, dpsubs]
        self.discrete_edges = discrete_edges
    
    def selectFromDimension(self, dimension, category):
        assert dimension == 'edge', dimension
        special = category in self.discrete_edges
        return self.choices[special].selectFromDimension(dimension, category)
    

