"""A recalculation engine, something like a spreadsheet.

Goals:
 - Allow construction of a calculation in a flexible and declarative way.
 - Enable caching at any step in the calculation where it makes sense.

Terms:
 - Definition - defines one cachable step in a complex calculation.
 - ParameterController - Sets parameter scope rules on a DAG of Definitions.
 - Calculator - An instance of an internally caching function.
 - Category - An arbitrary label.
 - Dimension - A named set of categories.
 - Scope - A subset of the categories from each dimension.
 - Setting - A variable (Var) or constant (ConstVal).
 - Assignments - A mapping from Scopes to Settings.
 - Cell - Evaluates one Scope of one Definition.
 - OptPar - A cell with indegree 0.

Structure:
 - A Calculator holds a list of Cells: OptPars and EvaluatedCells.
 - EvaluatedCells take their arguments from other Cells.
 - Each type of cell (motifs, Qs, Psubs) made by a different CalculationDefn.
 - No two cells from the same CalculationDefn have the same inputs, so nothing
   is calculated twice.

Interface:
  1) Define a function for each step in the calculation.
  2) Instantiate a DAG of ParamDefns and CalcDefns, each CalcDefn like
     CalcDefn(f)(*args) where 'f' is one of your functions and the '*args'
     are Defns that correspond to the arguments of 'f'.
  3) With your final CalcDefn called say 'top', PC = ParameterController(top)
     to get a ParameterController.
  4) PC.assign_all(param, value=value, **scope) to define the parameter
     scopes.  'value' can be a constant float or an instance of Var.
  5) calculator = PC.make_calculator() to get a live Calculator.
  6) calculator.optimise() etc.

Caching:
  In addition to the caching provided by the update strategy (not recalculating
  anything that hasn't changed), the calculator keeps a 1-deep cache of the
  previous value for each cell so that it has a 1-deep undo capability.  This
  is ideal for the behaviour of a one-change-at-a-time simanneal optimiser,
  which backtracks when a new value isn't accepted, ie it tries sequences like:
    [0,0,0] [0,0,3] [0,8,0] [7,0,0] [0,0,4] [0,6,0] ...
  when it isn't making progress, and
    [0,0,0] [0,0,3] [0,8,3] [7,8,3] [7,8,9] ...
  when it's having a lucky streak.

  Each cell knows when it's out of date, but doesn't know why (ie: what input
  changed) which limits the undo strategy to all-or-nothing.  An optimiser that
  tried values
    [0,0,0] [0,3,8] [0,3,0] ...
  (ie: the third step is a recombination of the previous two) would not get
  any help from caching.  This does keep things simple and fast though.

Recycling:
  If defn.recycling is True then defn.calc() will be passed the previous
  result as its first argument so it can be reused.  This is to avoid
  having to reallocate memory for say a large numpy array just at the
  very moment that an old one of the same shape is being disposed of.
  To prevent recycling from invalidating the caching system 3 values are
  stored for each cell - current, previous and spare.  The spare value is
  the one to be used next for recycling.
"""


from collections import defaultdict

import numpy

from cogent3.maths.optimisers import ParameterOutOfBoundsError
from cogent3.maths.util import proportions_to_ratios, ratios_to_proportions
from cogent3.util.dict_array import DictArrayTemplate

from .calculation import ConstCell, EvaluatedCell, LogOptPar, OptPar
from .scope import ParameterController, _Defn, _LeafDefn, _NonLeafDefn
from .setting import ConstVal, Var


# In this module we bring together scopes, settings and calculations.
# Most of the classes are 'Defns' with their superclasses in scope.py.
# These supply a make_cells() method which instantiates 'Cell'
# classes from calculation.py


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


DIM_PLURALS = {
    "bin": "bins",
    "edge": "edges",
    "locus": "loci",
    "motif": "motifs",
    "position": "positions",
}


def select_dim_label_value(dim, value):
    """returns singular / plural for dim given value

    Required for construction of param rules"""
    if dim not in DIM_PLURALS:
        return dim

    if type(value) == str:
        return dim

    if len(value) > 1:
        dim = DIM_PLURALS[dim]
    else:
        value = value[0]

    return dim, value


class CalculationDefn(_NonLeafDefn):
    """Defn for a derived value.  In most cases use CalcDefn instead

    The only reason for subclassing this directly would be to override
    .make_calc_function() or setup()."""

    recycling = False

    # positional arguments are inputs to this step of the calculation,
    # keyword arguments are passed on to self.setup(), likely to end up
    # as static attributes of this CalculationDefn, to be used (as self.X)
    # by its 'calc' method.

    def make_likelihood_function(self):
        return ParameterController(self)

    def setup(self):
        pass

    def make_calc_function(self):
        return self.calc

    def make_cell(self, *args):
        calc = self.make_calc_function()
        return EvaluatedCell(
            self.name, calc, args, recycling=self.recycling, default=self.default
        )

    def make_cells(self, input_soup, variable=None):
        # input soups contains all necessary values for calc on self.
        # Going from defns to cells.
        cells = []
        for input_nums in self.uniq:
            args = []
            for (arg, u) in zip(self.args, input_nums):
                arg = input_soup[id(arg)][u]
                args.append(arg)
            cell = self.make_cell(*args)
            cells.append(cell)
        return (cells, cells)


class _FuncDefn(CalculationDefn):
    def __init__(self, calc, *args, **kw):
        self.calc = calc
        CalculationDefn.__init__(self, *args, **kw)


# Use this rather than having to subclass CalculationDefinition
# just to supply the 'calc' method.
class CalcDefn(object):
    """CalcDefn(function)(arg1, arg2)"""

    def __init__(self, calc, name=None, **kw):
        self.calc = calc

        if name is None:
            name = self.calc.__name__
        else:
            assert isinstance(name, str), name
        kw["name"] = name
        self.kw = kw

    def __call__(self, *args):
        return _FuncDefn(self.calc, *args, **self.kw)


class WeightedPartitionDefn(CalculationDefn):
    """Uses a PartitionDefn (ie: N-1 optimiser parameters) to make
    an array of floats with weighted average of 1.0"""

    def __init__(self, weights, name=None):
        N = len(weights.bin_names)
        partition = PartitionDefn(size=N, name=name + "_partition")
        partition.user_param = False
        super(WeightedPartitionDefn, self).__init__(
            weights, partition, name=name + "_distrib"
        )

    def calc(self, weights, values):
        scale = numpy.sum(weights * values)
        return values / scale


class MonotonicDefn(WeightedPartitionDefn):
    """Uses a PartitionDefn (ie: N-1 optimiser parameters) to make
    an ordered array of floats with weighted average of 1.0"""

    def calc(self, weights, increments):
        values = numpy.add.accumulate(increments)
        scale = numpy.sum(weights * values)
        return values / scale


class GammaDefn(MonotonicDefn):
    """Uses 1 optimiser parameter to define a gamma distribution, divides
    the distribution into N equal probability bins and makes an array of
    their medians.  If N > 2 medians are approx means so their average
    is approx 1.0, but not quite so we scale them to make it exactly 1.0"""

    name = "gamma"

    def __init__(
        self, weights, name=None, default_shape=1.0, extra_label=None, dimensions=()
    ):
        name = self.make_name(name, extra_label)
        shape = PositiveParamDefn(
            name + "_shape", default=default_shape, dimensions=dimensions, lower=1e-2
        )
        CalculationDefn.__init__(self, weights, shape, name=name + "_distrib")

    def calc(self, weights, a):
        from cogent3.maths.stats.distribution import gdtri

        weights = weights / numpy.sum(weights)
        percentiles = numpy.add.accumulate(weights) - weights * 0.5
        medians = numpy.array([gdtri(a, a, p) for p in percentiles])
        scale = numpy.sum(medians * weights)
        # assert 0.5 < scale < 2.0, scale # medians as approx. to means.
        return medians / scale


class _InputDefn(_LeafDefn):
    user_param = True

    def __init__(
        self, name=None, default=None, dimensions=None, lower=None, upper=None, **kw
    ):
        _LeafDefn.__init__(self, name=name, dimensions=dimensions, **kw)
        if default is not None:
            if hasattr(default, "__len__"):
                default = numpy.array(default)
            self.default = default
        # these two have no effect on constants
        if lower is not None:
            self.lower = lower
        if upper is not None:
            self.upper = upper

    def make_likelihood_function(self):
        return ParameterController(self)

    def update_from_calculator(self, calc):
        outputs = calc.get_current_cell_values_for_defn(self)
        for (output, setting) in zip(outputs, self.uniq):
            # catch cases where parameters fall outside bounds due to precision
            if setting.is_constant:
                ...  # block trying other conditions
            elif setting.lower and output < setting.lower:
                if not numpy.allclose(output, setting.lower):
                    raise ParameterOutOfBoundsError(
                        f"calculator value {output} for {self.name!r} is < {setting.lower}"
                    )
                output = setting.lower
            elif setting.upper and output > setting.upper:
                if not numpy.allclose(output, setting.upper):
                    raise ParameterOutOfBoundsError(
                        f"calculator value {output} for {self.name!r} is > {setting.upper}"
                    )
                output = setting.upper
            setting.value = output

    def get_num_free_params(self):
        (cells, outputs) = self.make_cells({}, None)
        return len([c for c in cells if isinstance(c, OptPar)])

    def _get_scoped_params(self, keys, dim_indices):
        result = {}
        for index in dim_indices:
            dim_name = self.valid_dimensions[index]
            value = list(sorted(set([k[index] for k in keys])))
            dim_name, value = select_dim_label_value(dim_name, value)
            result[dim_name] = value
        return result

    def get_param_rules(self):
        """returns list of param rule dicts for this parameter"""
        num_valid_dims = len(self.valid_dimensions)
        # todo replace following with self.used_dimensions()
        dimensioned = {
            k: set(v) for k, v in zip(range(num_valid_dims), zip(*self.index))
        }

        discard = [k for k, v in dimensioned.items() if len(v) == 1]
        for k in discard:
            del dimensioned[k]

        dimensioned = list(dimensioned.keys())

        scoped = defaultdict(list)

        for k, v in self.index.items():
            scoped[v].append(k)

        is_global = len(scoped) == 1

        rules = []
        is_probs = isinstance(self, PartitionDefn)
        names = None if not is_probs else self.bin_names
        for index, keys in scoped.items():
            rule = dict(par_name=self.name)
            dimms = self._get_scoped_params(keys, dimensioned)
            rule.update(dimms)
            rule.update(
                self.uniq[index].get_param_rule_dict(names=names, is_probs=is_probs)
            )
            if self.independent_by_default:
                for d in dimms:
                    if type(dimms[d]) == list:
                        rule["is_independent"] = False
                        break
            rules.append(rule)

        if is_global:
            assert len(rules) == 1
            for s, p in DIM_PLURALS.items():
                rules[0].pop(s, None)
                rules[0].pop(p, None)

        return rules


class ParamDefn(_InputDefn):
    """Defn for an optimisable, scalar input to the calculation"""

    numeric = True
    const_by_default = False
    independent_by_default = False
    opt_par_class = OptPar

    # These can be overridden in a subclass or the constructor
    default = 1.0
    lower = -1e10
    upper = +1e10

    def make_default_setting(self):
        return Var(bounds=(self.lower, self.default, self.upper))

    def check_setting_is_valid(self, setting):
        pass

    def make_cells(self, input_soup=None, variable=None):
        input_soup = input_soup or {}
        uniq_cells = []
        for (i, v) in enumerate(self.uniq):
            scope = [key for key in self.assignments if self.assignments[key] is v]
            if v.is_constant or (variable is not None and variable is not v):
                cell = ConstCell(self.name, v.value)
            else:
                cell = self.opt_par_class(self.name, scope, v.get_bounds())
            uniq_cells.append(cell)

        return (uniq_cells, uniq_cells)


# Example / basic ParamDefn subclasses


class PositiveParamDefn(ParamDefn):
    lower = 0.0


class ProbabilityParamDefn(PositiveParamDefn):
    upper = 1.0


class RatioParamDefn(PositiveParamDefn):
    lower = 1e-6
    upper = 1e6
    opt_par_class = LogOptPar


class NonScalarDefn(_InputDefn):
    """Defn for an array or other such object that is an input but
    can not be optimised"""

    user_param = False
    numeric = False
    const_by_default = True
    independent_by_default = False
    default = None

    def make_default_setting(self):
        if self.default is None:
            return None
        else:
            return ConstVal(self.default)

    def check_setting_is_valid(self, setting):
        if not isinstance(setting, ConstVal):
            raise ValueError(f"{self.name} can only be constant")

    def make_cells(self, input_soup=None, variable=None):
        input_soup = input_soup or {}
        if None in self.uniq:
            if [v for v in self.uniq if v is not None]:
                scope = [
                    key for key in self.assignments if self.assignments[key] is None
                ]
                msg = 'Unoptimisable input "%%s" not set for %s' % scope
            else:
                msg = 'Unoptimisable input "%s" not given'
            raise ValueError(msg % self.name)
        uniq_cells = [ConstCell(self.name, v.value) for v in self.uniq]
        return (uniq_cells, uniq_cells)

    def get_num_free_params(self):
        return 0

    def update_from_calculator(self, calc):
        pass


def _ratio_to_proportion(*ratios):
    return numpy.asarray(ratios_to_proportions(1.0, ratios))


class PartitionDefn(_InputDefn):
    """A partition such as mprobs can be const or optimised.  Optimised is
    a bit tricky since it isn't just a scalar."""

    numeric = False  # well, not scalar anyway
    const_by_default = False
    independent_by_default = False

    def __init__(
        self, default=None, name=None, dimensions=None, dimension=None, size=None, **kw
    ):
        assert name
        if size is not None:
            pass
        elif default is not None:
            size = len(default)
        elif dimension is not None:
            size = len(dimension[1])
        self.size = size
        if dimension is not None:
            self.internal_dimension = dimension
            (dim_name, dim_cats) = dimension
            self.bin_names = dim_cats
            self.array_template = DictArrayTemplate(dim_cats)
            self.internal_dimensions = (dim_name,)
        if default is None:
            default = self._make_default_value()
        elif self.array_template is not None:
            default = self.array_template.unwrap(default)
        else:
            default = numpy.asarray(default)
        _InputDefn.__init__(
            self, name=name, default=default, dimensions=dimensions, **kw
        )
        self.check_value_is_valid(default, True)

    def _make_default_value(self):
        return numpy.array([1.0 / self.size] * self.size)

    def make_default_setting(self):
        # return ConstVal(self.default)
        return Var((None, self.default.copy(), None))

    def check_setting_is_valid(self, setting):
        value = setting.get_default_value()
        return self.check_value_is_valid(value, setting.is_constant)

    def check_value_is_valid(self, value, is_constant):
        if value.shape != (self.size,):
            raise ValueError(
                "Wrong array shape %s for %s, expected (%s,)"
                % (value.shape, self.name, self.size)
            )
        for part in value:
            if part < 0:
                raise ValueError(f"Negative probability in {self.name}")
            if part > 1:
                raise ValueError(f"Probability > 1 in {self.name}")
            if not is_constant:
                # 0 or 1 leads to log(0) or log(inf) in optimiser code
                if part == 0:
                    raise ValueError(f"Zeros allowed in {self.name} only when constant")
                if part == 1:
                    raise ValueError(f"Ones allowed in {self.name} only when constant")
        if abs(sum(value) - 1.0) > 0.00001:
            raise ValueError(
                f"Elements of {self.name} must sum to 1.0, not {sum(value)}"
            )

    def _make_partition_cell(self, name, scope, value):
        # This was originally put in its own function so as to provide a
        # closure containing the value of sum(value), which is no longer
        # required since it is now always 1.0.
        assert abs(sum(value) - 1.0) < 0.00001
        ratios = proportions_to_ratios(value)
        ratios = [LogOptPar(name + "_ratio", scope, (1e-6, r, 1e6)) for r in ratios]

        partition = EvaluatedCell(name, _ratio_to_proportion, tuple(ratios))
        return (ratios, partition)

    def make_cells(self, input_soup=None, variable=None):
        input_soup = input_soup or {}
        uniq_cells = []
        all_cells = []
        for (i, v) in enumerate(self.uniq):
            if v is None:
                raise ValueError(f"input {self.name} not set")
            assert hasattr(v, "get_default_value"), v
            value = v.get_default_value()
            assert hasattr(value, "shape"), value
            assert value.shape == (self.size,)
            scope = [key for key in self.assignments if self.assignments[key] is v]
            assert value is not None
            if v.is_constant or (variable is not None and variable is not v):
                partition = ConstCell(self.name, value)
            else:
                (ratios, partition) = self._make_partition_cell(self.name, scope, value)
                all_cells.extend(ratios)
            all_cells.append(partition)
            uniq_cells.append(partition)
        return (all_cells, uniq_cells)


def NonParamDefn(name, dimensions=None, **kw):
    # Just to get 2nd arg as dimensions
    return NonScalarDefn(name=name, dimensions=dimensions, **kw)


class ConstDefn(NonScalarDefn):
    # This isn't really needed - just use NonParamDefn
    name_required = False

    def __init__(self, value, name=None, **kw):
        NonScalarDefn.__init__(self, default=value, name=name, **kw)

    def check_setting_is_valid(self, setting):
        if setting is not None and setting.value is not self.default:
            raise ValueError(f"{self.name} is constant")


class SelectForDimension(_Defn):
    """A special kind of Defn used to bridge from Defns where a particular
    dimension is wrapped up inside an array to later Defns where each
    value has its own Defn, eg: gamma distributed rates"""

    name = "select"
    user_param = True
    numeric = True  # not guarenteed!
    internal_dimensions = ()

    def __init__(self, arg, dimension, name=None):
        assert not arg.activated, arg.name
        if name is not None:
            self.name = name
        _Defn.__init__(self)
        self.args = (arg,)
        self.arg = arg
        self.valid_dimensions = arg.valid_dimensions
        if dimension not in self.valid_dimensions:
            self.valid_dimensions = self.valid_dimensions + (dimension,)
        self.dimension = dimension
        arg.add_client(self)

    def update(self):
        for scope_t in self.assignments:
            scope = dict(list(zip(self.valid_dimensions, scope_t)))
            scope2 = dict(
                (n, v) for (n, v) in list(scope.items()) if n != self.dimension
            )
            input_num = self.arg.output_ordinal_for(scope2)
            pos = self.arg.bin_names.index(scope[self.dimension])
            self.assignments[scope_t] = (input_num, pos)
        self._update_from_assignments()
        self.values = [self.arg.values[i][p] for (i, p) in self.uniq]

    def _select(self, arg, p):
        return arg[p]

    def make_cells(self, input_soup, variable=None):
        cells = []
        distribs = input_soup[id(self.arg)]
        for (input_num, bin_num) in self.uniq:
            cell = EvaluatedCell(
                self.name, (lambda x, p=bin_num: x[p]), (distribs[input_num],)
            )
            cells.append(cell)
        return (cells, cells)


# Some simple CalcDefns

# SumDefn = CalcDefn(lambda *args:sum(args), 'sum')
# ProductDefn = CalcDefn(lambda *args:numpy.product(args), 'product')
# CallDefn = CalcDefn(lambda func,*args:func(*args), 'call')


class SwitchDefn(CalculationDefn):
    name = "switch"

    def calc(self, condition, *args):
        return args[condition]

    def get_shortcut_cell(self, condition, *args):
        if condition.is_constant:
            return self.calc(self, condition.value, *args)


class VectorMatrixInnerDefn(CalculationDefn):
    name = "evolve"

    def calc(self, pi, psub):
        return numpy.inner(pi, psub)

    def get_shortcut_cell(self, pi, psub):
        if psub.is_stationary:
            return pi


class SumDefn(CalculationDefn):
    name = "sum"

    def calc(self, *args):
        return sum(args)


class ProductDefn(CalculationDefn):
    name = "product"

    def calc(self, *args):
        return numpy.product(args)


class CallDefn(CalculationDefn):
    name = "call"

    def calc(self, func, *args):
        return func(*args)


__all__ = [
    "ConstDefn",
    "NonParamDefn",
    "CalcDefn",
    "SumDefn",
    "ProductDefn",
    "CallDefn",
] + [
    n
    for (n, c) in list(vars().items())
    if (isinstance(c, type) and issubclass(c, _Defn) and n[0] != "_")
    or isinstance(c, CalcDefn)
]
