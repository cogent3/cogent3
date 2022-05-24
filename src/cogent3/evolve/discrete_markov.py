import numpy

from cogent3.recalculation.definition import (
    ConstCell,
    DictArrayTemplate,
    EvaluatedCell,
    PartitionDefn,
)


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def _make_array(*x):
    return numpy.array(x)


class PsubMatrixDefn(PartitionDefn):
    "Square 2D array made of 1D partitions"

    numeric = False  # well, not scalar anyway
    const_by_default = False
    independent_by_default = True

    def __init__(
        self, default=None, name=None, dimensions=None, dimension=None, size=None, **kw
    ):
        PartitionDefn.__init__(self, default, name, dimensions, dimension, size, **kw)

        (dim_name, dim_cats) = self.internal_dimension
        self.internal_dimensions = (dim_name, dim_name + "2")
        self.array_template = DictArrayTemplate(dim_cats, dim_cats)

    def _make_default_value(self):
        # Purely flat default doesn't work well so start at approx t=0.5
        flat = numpy.ones([self.size, self.size], float) / self.size
        diag = numpy.identity(self.size, float)
        return (flat + diag) / 2

    def check_value_is_valid(self, value, is_constant):
        if value.shape != (self.size, self.size):
            raise ValueError(
                "Wrong array shape %s for %s, expected (%s,%s)"
                % (value.shape, self.name, self.size, self.size)
            )
        for part in value:
            PartitionDefn.check_value_is_valid(self, part, is_constant)

    def make_cells(self, input_soup=None, variable=None):
        uniq_cells = []
        all_cells = []
        for (i, v) in enumerate(self.uniq):
            if v is None:
                raise ValueError(f"input {self.name} not set")
            assert hasattr(v, "get_default_value"), v
            value = v.get_default_value()
            assert hasattr(value, "shape"), value
            assert value.shape == (self.size, self.size)
            scope = [key for key in self.assignments if self.assignments[key] is v]
            if v.is_constant or (variable is not None and variable is not v):
                matrix = ConstCell(self.name, value)
            else:
                rows = []
                for part in value:
                    (ratios, partition) = self._make_partition_cell(
                        self.name + "_part", scope, part
                    )
                    all_cells.extend(ratios)
                    rows.append(partition)
                all_cells.extend(rows)
                matrix = EvaluatedCell(self.name, _make_array, rows)
            all_cells.append(matrix)
            uniq_cells.append(matrix)
        return (all_cells, uniq_cells)


class PartialyDiscretePsubsDefn(object):
    def __init__(self, alphabet, psubs, discrete_edges):
        motifs = tuple(alphabet)
        dpsubs = PsubMatrixDefn(
            name="dpsubs",
            dimension=("motif", motifs),
            default=None,
            dimensions=("locus", "edge"),
        )
        self.choices = [psubs, dpsubs]
        self.discrete_edges = discrete_edges

    def select_from_dimension(self, dimension, category):
        assert dimension == "edge", dimension
        special = category in self.discrete_edges
        return self.choices[special].select_from_dimension(dimension, category)
