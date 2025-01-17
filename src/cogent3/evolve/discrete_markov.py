import numpy

from cogent3.recalculation.definition import (
    ConstCell,
    DictArrayTemplate,
    EvaluatedCell,
    PartitionDefn,
)


def _make_array(*x):
    return numpy.array(x)


class PsubMatrixDefn(PartitionDefn):
    "Square 2D array made of 1D partitions"

    numeric = False  # well, not scalar anyway
    const_by_default = False
    independent_by_default = True

    def __init__(
        self,
        default=None,
        name=None,
        dimensions=None,
        dimension=None,
        size=None,
        **kw,
    ) -> None:
        PartitionDefn.__init__(self, default, name, dimensions, dimension, size, **kw)

        (dim_name, dim_cats) = self.internal_dimension
        self.internal_dimensions = (dim_name, dim_name + "2")
        self.array_template = DictArrayTemplate(dim_cats, dim_cats)

    def _make_default_value(self):
        # Purely flat default doesn't work well so start at approx t=0.5
        flat = numpy.ones([self.size, self.size], float) / self.size
        diag = numpy.identity(self.size, float)
        return (flat + diag) / 2

    def check_value_is_valid(self, value, is_constant) -> None:
        if value.shape != (self.size, self.size):
            msg = f"Wrong array shape {value.shape} for {self.name}, expected ({self.size},{self.size})"
            raise ValueError(
                msg,
            )
        for part in value:
            PartitionDefn.check_value_is_valid(self, part, is_constant)

    def make_cells(self, input_soup=None, variable=None):
        uniq_cells = []
        all_cells = []
        for _i, v in enumerate(self.uniq):
            if v is None:
                msg = f"input {self.name} not set"
                raise ValueError(msg)
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
                        self.name + "_part",
                        scope,
                        part,
                    )
                    all_cells.extend(ratios)
                    rows.append(partition)
                all_cells.extend(rows)
                matrix = EvaluatedCell(self.name, _make_array, rows)
            all_cells.append(matrix)
            uniq_cells.append(matrix)
        return (all_cells, uniq_cells)


class PartialyDiscretePsubsDefn:
    def __init__(self, alphabet, psubs, discrete_edges) -> None:
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
