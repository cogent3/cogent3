#!/usr/bin/env python'
"""Instances of these classes are assigned to different parameter/scopes
by a parameter controller"""

from cogent3.util.misc import adjusted_gt_minprob


class Setting:
    def get_param_rule_dict(self, names=None, is_probs=False):
        """returns component for a param rule dict"""
        value = self.value
        if names is not None:
            assert len(self.value) == len(names)
            if is_probs:
                value = adjusted_gt_minprob(
                    adjusted_gt_minprob(self.value, minprob=1e-6),
                )
            if value.ndim == 1:
                value = dict(zip(names, value, strict=False))
            else:
                result = {}
                for i, n1 in enumerate(names):
                    result[n1] = {}
                    for j, n2 in enumerate(names):
                        result[n1][n2] = value[i, j]
                value = result

        if self.is_constant:
            result = {"value": value, "is_constant": True}
        else:
            result = {"init": value, "lower": self.lower, "upper": self.upper}
        return result


class Var(Setting):
    # placeholder for a single optimiser parameter
    is_constant = False

    def __init__(self, bounds=None) -> None:
        if bounds is None:
            bounds = (None, None, None)
        else:
            assert len(bounds) == 3, bounds
        (self.lower, self.value, self.upper) = bounds

    def get_bounds(self):
        return (self.lower, self.value, self.upper)

    def get_default_value(self):
        return self.value

    def __str__(self) -> str:
        return "Var"  # short as in table

    def __repr__(self) -> str:
        constraints = []
        for template, bound in [
            ("%s<", self.lower),
            ("(%s)", self.value),
            ("<%s", self.upper),
        ]:
            if bound is not None:
                constraints.append(template % bound)
        return f"Var({' '.join(constraints)})"


class ConstVal(Setting):
    # not to be confused with defns.Const.  This is like a Var,
    # assigned to a parameter which may have other Var values
    # for other scopes.
    is_constant = True

    # Val interface
    def __init__(self, value) -> None:
        self.value = value

    def __str__(self) -> str:
        return repr(self.value)  # short as in table

    def __repr__(self) -> str:
        return f"ConstVal({self.value!r})"

    # indep useful sometimes!
    # def __eq__(self, other):
    #    return type(self) is type(other) and other.value == self.value

    def get_default_value(self):
        return self.value

    def get_bounds(self):
        return (None, self.value, None)
