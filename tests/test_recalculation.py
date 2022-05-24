from unittest import TestCase, main

from cogent3.recalculation.definition import CalcDefn, ParamDefn
from cogent3.recalculation.scope import (
    InvalidDimensionError,
    InvalidScopeError,
)


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class RecalculationTest(TestCase):
    def test_recalculation(self):
        def add(*args):
            return sum(args)

        top = CalcDefn(add)(ParamDefn("A"), ParamDefn("B"))
        pc = top.make_likelihood_function()
        f = pc.make_calculator()

        self.assertEqual(f.get_value_array(), [1.0, 1.0])
        self.assertEqual(f([3.0, 4.25]), 7.25)
        self.assertEqual(f.change([(1, 4.5)]), 7.5)
        self.assertEqual(f.get_value_array(), [3.0, 4.5])

        # Now with scopes.  We will set up the calculation
        # result = (Ax+Bx) + (Ay+By) + (Az+Bz)

        # A and B will remain distinct parameters, but x,y and z are merely scopes - ie:
        # it may be the case that Ax = Ay = Az, and that may simplify the calculation, but
        # we will never even notice if Ax = Bx.
        # Each scope dimension (here there is just one, 'category') must be collapsed away
        # at some point towards the end of the calculation if the calculation is to produce
        # a scalar result.  Here this is done with the select_from_dimension method.

        a = ParamDefn("A", dimensions=["category"])
        b = ParamDefn("B", dimensions=["category"])
        mid = CalcDefn(add, name="mid")(a, b)
        top = CalcDefn(add)(
            mid.select_from_dimension("category", "x"),
            mid.select_from_dimension("category", "y"),
            mid.select_from_dimension("category", "z"),
        )

        # or equivalently:
        # top = CalcDefn(add, *mid.acrossDimension('category', ['x', 'y', 'z']))

        pc = top.make_likelihood_function()
        f = pc.make_calculator()

        self.assertEqual(str(f.get_value_array()), "[1.0, 1.0]")

        # There are still only 2 inputs because the default scope
        # is global, ie: Ax == Ay == Az.  If we allow A to be
        # different in the x,y and z categories and set their
        # initial values to 2.0:

        pc.assign_all("A", value=2.0, independent=True)
        f = pc.make_calculator()

        self.assertEqual(str(f.get_value_array()), "[1.0, 2.0, 2.0, 2.0]")

        # Now we have A local and B still global, so the calculation is
        # (Ax+B) + (Ay+B) + (Az+B) with the input parameters being
        # [B, Ax, Ay, Az], so:

        self.assertEqual(f([1.0, 2.0, 2.0, 2.0]), 9.0)
        self.assertEqual(f([0.25, 2.0, 2.0, 2.0]), 6.75)

        # Constants do not appear in the optimisable inputs.
        # Set one of the 3 A values to be a constant and there
        # will be one fewer optimisable parameters:

        pc.assign_all("A", scope_spec={"category": "z"}, const=True)
        f = pc.make_calculator()

        self.assertEqual(str(f.get_value_array()), "[1.0, 2.0, 2.0]")

        # The parameter controller should catch cases where the specified scope
        # does not exist:

        with self.assertRaises(InvalidScopeError):
            pc.assign_all("A", scope_spec={"category": "nosuch"})
        with self.assertRaises(InvalidDimensionError):
            pc.assign_all("A", scope_spec={"nonsuch": "nosuch"})

        # It is complicated guesswork matching the parameters you expect with positions in
        # the value array, let alone remembering whether or not they are presented to the
        # optimiser as logs, so .get_value_array(), .change() and .__call__() should only be
        # used by optimisers.  For other purposes there is an alternative, human friendly
        # interface:

        pc.update_from_calculator(f)
        self.assertEqual(pc.get_param_value("A", category="x"), 2, 0)
        self.assertEqual(pc.get_param_value("B", category=["x", "y"]), 1.0)

        # Despite the name, .get_param_value can get the value from any step in the
        # calculation, so long as it has a unique name.

        self.assertEqual(pc.get_param_value("mid", category="x"), 3.0)

        # For bulk retrieval of parameter values by parameter name and scope name there is
        # the .get_param_value_dict() method:

        vals = pc.get_param_value_dict(["category"])
        self.assertEqual(vals["A"]["x"], 2.0)

        # Here is a function that is more like a likelihood function, in that it has a
        # maximum:

        def curve(x, y):
            return 0 - (x ** 2 + y ** 2)

        top = CalcDefn(curve)(ParamDefn("X"), ParamDefn("Y"))
        pc = top.make_likelihood_function()
        f = pc.make_calculator()

        # Now ask it to find the maximum.  It is a simple function with only one local
        # maximum so local optimisation should be enough:

        f.optimise(local=True, show_progress=False)
        pc.update_from_calculator(f)

        # There were two parameters, X and Y, and at the maximum they should both be 0.0:

        self.assertEqual(pc.get_param_value("Y"), 0.0)
        self.assertEqual(pc.get_param_value("X"), 0.0)

        # Because this function has a maximum it is possible to ask it for a confidence
        # interval around a parameter, ie: how far from 0.0 can we move x before f(x,y)
        # falls bellow f(X,Y)-dropoff:

        self.assertEqual(
            pc.get_param_interval("X", dropoff=4, xtol=0.0), (-2.0, 0.0, 2.0)
        )

        # We test the ability to omit xtol. Due to precision issues we convert the returned value to a string.

        self.assertTrue(
            "-2.0, 0.0, 2.0"
            == "%.1f, %.1f, %.1f" % pc.get_param_interval("X", dropoff=4)
        )

        # And finally intervals can be calculated in bulk by passing a dropoff value to
        # .get_param_value_dict():

        self.assertEqual(
            pc.get_param_value_dict([], dropoff=4, xtol=0.0)["X"], (-2.0, 0.0, 2.0)
        )

        # For likelihood functions it is more convenient to provide 'p' rather than
        # 'dropoff', dropoff = chdtri(1, p) / 2.0.  Also in general you won't need ultra precise answers,
        # so don't use 'xtol=0.0', that's just to make the doctest work.
        pc.graphviz()


if __name__ == "__main__":
    main()
