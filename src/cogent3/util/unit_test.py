#!/usr/bin/env python
"""Extension of the built-in unittest framework for floating-point comparisons.

Specific Extensions:

assertFloatEqual, assertFloatEqualAbs, and assertFloatEqualRel give fine-
grained control over how floating point numbers (or lists thereof) are tested
for equality.

assertContains and assertNotContains give more helpful error
messages when testing whether an observed item is present or absent in a set
of possiblities. Ditto assertGreaterThan, assertLessThan, and assertIsProb.

assertSameItems and assertEqualItems test the items in a list
for pairwise identity and equality respectively (i.e. the observed and
expected values must have the same number of each item, though the order can
differ).

assertSimilarMeans and assertSimilarFreqs allow you to test stochastic results
by setting an explicit P-value and checking that the result is not improbable
given the expected P-value. Please use these instead of guessing confidence
intervals! The major advantage is that you can reset the P-value gloabally over
the whole test suite, so that rare failures don't occur every time.

"""
from unittest import TestCase as orig_TestCase
from unittest import TestSuite, findTestCases, main

import numpy

from numpy import (
    array,
    asarray,
    isfinite,
    logical_and,
    logical_or,
    ravel,
    testing,
    zeros,
)

from cogent3.maths.stats.test import G_ind, t_two_sample
from cogent3.util.misc import recursive_flatten


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = [
    "Rob Knight",
    "Peter Maxwell",
    "Sandra Smit",
    "Zongzhi Liu",
    "Micah Hamady",
    "Daniel McDonald",
]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

# SUPPORT2425


class FakeRandom(object):
    """Drop-in substitute for random.random that provides items from list."""

    def __init__(self, data, circular=False):
        """Returns new FakeRandom object, using list of items in data.

        circular: if True (default is False), wraps the list around. Otherwise,
        raises IndexError when we run off the end of the list.

        WARNING: data must always be iterable, even if it's a single item.
        """
        self._data = data
        self._ptr = -1
        self._circular = circular

    def __call__(self, *args, **kwargs):
        """Returns next item from the list in self._data.

        Raises IndexError when we run out of data.
        """
        self._ptr += 1
        # wrap around if circular
        if self._circular:
            if self._ptr >= len(self._data):
                self._ptr = 0
        return self._data[self._ptr]


class TestCase(orig_TestCase):
    """Adds some additional utility methods to unittest.TestCase.

    Notably, adds facilities for dealing with floating point numbers,
    and some common templates for replicated tests.

    BEWARE: Do not start any method with 'test' unless you want it to actually
    run as a test suite in every instance!
    """

    _suite_pvalue = None  # see TestCase._set_suite_pvalue()

    def _get_values_from_matching_dicts(self, d1, d2):
        """Gets corresponding values from matching dicts"""
        if set(d1) != set(d2):
            return None
        # might not be in same order
        return list(d1.values()), [d2[k] for k in d1]

    def errorCheck(self, call, known_errors):
        """Applies function to (data, error) tuples, checking for error
        """
        for (data, error) in known_errors:
            self.assertRaises(error, call, data)

    def valueCheck(self, call, known_values, arg_prefix="", eps=None):
        """Applies function to (data, expected) tuples, treating data as args
        """
        for (data, expected) in known_values:
            observed = eval("call(" + arg_prefix + "data)")
            try:
                allowed_diff = float(eps)
            except TypeError:
                self.assertEqual(observed, expected)
            else:
                self.assertFloatEqual(observed, expected, allowed_diff)

    def assertFloatEqualRel(self, obs, exp, eps=1e-6):
        """Tests whether two floating point numbers/arrays are approx. equal.

        Checks whether the distance is within epsilon relative to the value
        of the sum of observed and expected. Use this method when you expect
        the difference to be small relative to the magnitudes of the observed
        and expected values.

        Note: for arbitrary objects, need to compare the specific attribute
        that's numeric, not the whole object, using this method.
        """
        # do array check first
        # note that we can't use array ops to combine, because we need to check
        # at each element whether the expected is zero to do the test to avoid
        # floating point error.
        # WARNING: numpy iterates over objects that are not regular Python
        # floats/ints, so need to explicitly catch scalar values and prevent
        # cast to array if we want the exact object to print out correctly.
        is_array = False
        if hasattr(obs, "keys") and hasattr(exp, "keys"):  # both dicts?
            result = self._get_values_from_matching_dicts(obs, exp)
            if result:
                obs, exp = result
        else:
            try:
                iter(obs)
                iter(exp)
            except TypeError:
                obs = [obs]
                exp = [exp]
            else:
                try:
                    arr_obs = array(obs)
                    arr_exp = array(exp)
                    arr_diff = arr_obs - arr_exp
                    if arr_obs.shape != arr_exp.shape:
                        self.fail(
                            "Wrong shape: Got %s, but expected %s"
                            % (repr(obs), repr(exp))
                        )
                    obs = arr_obs.ravel()
                    exp = arr_exp.ravel()
                    is_array = True
                except (TypeError, ValueError):
                    pass

        # shape mismatch can still get by...
        # explict cast is to work around bug in certain versions of numpy
        # installed version on osx 10.5
        if asarray(obs, object).shape != asarray(exp, object).shape:
            self.fail("Wrong shape: Got %s, but expected %s" % (obs, exp))

        for observed, expected in zip(obs, exp):
            # try the cheap comparison first
            if observed == expected:
                continue
            try:
                sum = float(observed + expected)
                diff = float(observed - expected)
                if sum == 0:
                    if is_array:
                        self.assertFalse(
                            abs(diff) > abs(eps),
                            "Got %s, but expected %s (diff was %s)"
                            % (repr(arr_obs), repr(arr_exp), repr(arr_diff)),
                        )
                    else:
                        self.assertFalse(
                            abs(diff) > abs(eps),
                            "Got %s, but expected %s (diff was %s)"
                            % (repr(observed), repr(expected), repr(diff)),
                        )

                else:
                    if is_array:
                        self.assertFalse(
                            abs(diff / sum) > abs(eps),
                            "Got %s, but expected %s (diff was %s)"
                            % (repr(arr_obs), repr(arr_exp), repr(arr_diff)),
                        )
                    else:
                        self.assertFalse(
                            abs(diff / sum) > abs(eps),
                            "Got %s, but expected %s (diff was %s)"
                            % (repr(observed), repr(expected), repr(diff)),
                        )
            except (TypeError, ValueError, AttributeError, NotImplementedError):
                self.fail("Got %s, but expected %s" % (repr(observed), repr(expected)))

    def assertFloatEqualAbs(self, obs, exp, eps=1e-6):
        """
        Tests whether two floating point numbers are approximately equal.

        Checks whether the absolute value of (a - b) is within epsilon. Use
        this method when you expect that one of the values should be very
        small, and the other should be zero.
        """
        # do array check first
        # note that we can't use array ops to combine, because we need to check
        # at each element whether the expected is zero to do the test to avoid
        # floating point error.
        if hasattr(obs, "keys") and hasattr(exp, "keys"):  # both dicts?
            result = self._get_values_from_matching_dicts(obs, exp)
            if result:
                obs, exp = result
        else:
            try:
                iter(obs)
                iter(exp)
            except TypeError:
                obs = [obs]
                exp = [exp]
            else:
                try:
                    arr_obs = array(obs)
                    arr_exp = array(exp)
                    if arr_obs.shape != arr_exp.shape:
                        self.fail(
                            "Wrong shape: Got %s, but expected %s"
                            % (repr(obs), repr(exp))
                        )
                    diff = arr_obs - arr_exp
                    self.assertFalse(
                        abs(diff).max() > eps,
                        "Got %s, but expected %s (diff was %s)"
                        % (repr(obs), repr(exp), repr(diff)),
                    )
                    return
                except (TypeError, ValueError):
                    pass
        # only get here if array comparison failed
        for observed, expected in zip(obs, exp):
            # cheap comparison first
            if observed == expected:
                continue
            try:
                diff = observed - expected
                self.assertFalse(
                    abs(diff) > abs(eps),
                    "Got %s, but expected %s (diff was %s)"
                    % (repr(observed), repr(expected), repr(diff)),
                )
            except (TypeError, ValueError, AttributeError, NotImplementedError):
                self.fail("Got %s, but expected %s" % (repr(observed), repr(expected)))

    def assertFloatEqual(self, obs, exp, eps=1e-6, rel_eps=None, abs_eps=None):
        """Tests whether two floating point numbers are approximately equal.

        If one of the arguments is zero, tests the absolute magnitude of the
        difference; otherwise, tests the relative magnitude.

        Use this method as a reasonable default.
        """
        obs = numpy.asarray(obs, dtype="O")
        exp = numpy.asarray(exp, dtype="O")
        obs = numpy.ravel(obs)
        exp = numpy.ravel(exp)

        if obs.shape != exp.shape:
            self.fail("Shape mismatch. Got, %s but expected %s" % (obs, exp))

        for observed, expected in zip(obs, exp):
            if self._is_equal(observed, expected):
                continue
            try:
                rel_eps = rel_eps or eps
                abs_eps = abs_eps or eps
                if (observed == 0) or (expected == 0):
                    self.assertFloatEqualAbs(observed, expected, abs_eps)
                else:
                    self.assertFloatEqualRel(observed, expected, rel_eps)
            except (TypeError, ValueError, AttributeError, NotImplementedError):
                self.fail("Got %s, but expected %s" % (repr(observed), repr(expected)))

    def _is_equal(self, observed, expected):
        """Returns True if observed and expected are equal, False otherwise."""
        # errors to catch: TypeError when obs is None
        tolist_errors = (AttributeError, ValueError, TypeError)

        try:
            obs = observed.tolist()
        except tolist_errors:
            obs = observed
        try:
            exp = expected.tolist()
        except tolist_errors:
            exp = expected
        return obs == exp

    def failUnlessEqual(self, observed, expected, msg=None):
        """Fail if the two objects are unequal as determined by !=

        Overridden to make error message enforce order of observed, expected.
        Use numpy.testing.assert_equal if ValueError, TypeError raised.
        """
        try:
            if not self._is_equal(observed, expected):
                raise self.failureException(
                    msg or "Got %s, but expected %s" % (repr(observed), repr(expected))
                )
        except (ValueError, TypeError) as e:
            # The truth value of an array with more than one element is
            # ambiguous. Use a.any() or a.all()
            # descriptor 'tolist' of 'numpy.generic' object needs an argument
            testing.assert_equal(observed, expected)

    def failIfEqual(self, observed, expected, msg=None):
        """Fail if the two objects are equal as determined by =="""
        try:
            self.assertEqual(observed, expected)
        except self.failureException:
            pass
        else:
            raise self.failureException(
                msg
                or "Observed %s and expected %s: shouldn't test equal"
                % (repr(observed), repr(expected))
            )

        # following needed to get our version instead of unittest's

    assertEqual = assertEquals = failUnlessEqual

    assertNotEqual = assertNotEquals = failIfEqual

    def assertEqualItems(self, observed, expected, msg=None):
        """Fail if the two items contain unequal elements"""
        obs_items = list(observed)
        exp_items = list(expected)
        if len(obs_items) != len(exp_items):
            raise self.failureException(
                msg
                or "Observed and expected are different lengths: %s and %s"
                % (len(obs_items), len(exp_items))
            )

        obs_items.sort()
        exp_items.sort()
        for index, (obs, exp) in enumerate(zip(obs_items, exp_items)):
            if obs != exp:
                raise self.failureException(
                    msg
                    or "Observed %s and expected %s at sorted index %s"
                    % (obs, exp, index)
                )

    def assertSameItems(self, observed, expected, msg=None):
        """Fail if the two items contain non-identical elements"""
        obs_items = list(observed)
        exp_items = list(expected)
        if len(obs_items) != len(exp_items):
            raise self.failureException(
                msg
                or "Observed and expected are different lengths: %s and %s"
                % (len(obs_items), len(exp_items))
            )

        obs_ids = [(id(i), i) for i in obs_items]
        exp_ids = [(id(i), i) for i in exp_items]
        obs_ids.sort()
        exp_ids.sort()
        for index, (obs, exp) in enumerate(zip(obs_ids, exp_ids)):
            o_id, o = obs
            e_id, e = exp
            if o_id != e_id:  # i.e. the ids are different
                raise self.failureException(
                    msg
                    or "Observed %s <%s> and expected %s <%s> at sorted index %s"
                    % (o, o_id, e, e_id, index)
                )

    def assertContains(self, observed, item, msg=None):
        """Fail if item not in observed"""
        try:
            if item in observed:
                return
        except (TypeError, ValueError):
            pass
        raise self.failureException(
            msg or "Item %s not found in %s" % (repr(item), repr(observed))
        )

    def assertNotContains(self, observed, item, msg=None):
        """Fail if item in observed"""
        try:
            if item not in observed:
                return
        except (TypeError, ValueError):
            return
        raise self.failureException(
            msg or "Item %s should not have been in %s" % (repr(item), repr(observed))
        )

    def assertGreaterThan(self, observed, value, msg=None):
        """Fail if observed is <= value"""
        try:
            if value is None or observed is None:
                raise ValueError
            if (asarray(observed) > value).all():
                return
        except:
            pass
        raise self.failureException(
            msg or "Observed %s has elements <= %s" % (repr(observed), repr(value))
        )

    def assertLessThan(self, observed, value, msg=None):
        """Fail if observed is >= value"""
        try:
            if value is None or observed is None:
                raise ValueError
            if (asarray(observed) < value).all():
                return
        except:
            pass
        raise self.failureException(
            msg or "Observed %s has elements >= %s" % (repr(observed), repr(value))
        )

    def assertIsProb(self, observed, msg=None):
        """Fail is observed is not between 0.0 and 1.0"""
        try:
            if observed is None:
                raise ValueError
            if (asarray(observed) >= 0.0).all() and (asarray(observed) <= 1.0).all():
                return
        except:
            pass
        raise self.failureException(
            msg or "Observed %s has elements that are not probs" % (repr(observed))
        )

    def _set_suite_pvalue(self, pvalue):
        """Sets the test suite pvalue to be used in similarity tests

        This value is by default None. The pvalue used in this case is
        specified in the test module itself. The purpose of this method is to
        set the pvalue to be used when running a massive test suite
        """
        self._suite_pvalue = pvalue

    def assertSimilarMeans(self, observed, expected, pvalue=0.01, msg=None):
        """Fail if observed p is lower than pvalue"""
        if self._suite_pvalue:
            pvalue = self._suite_pvalue

        observed, expected = asarray(observed), asarray(expected)

        t, p = t_two_sample(observed, expected)

        # handle case where all elements were the same
        if p is None or not isfinite(p):
            if not observed.shape:
                observed = observed.reshape((1,))
            if not expected.shape:
                expected = expected.reshape((1,))
            if observed[0] == expected[0]:
                return
        elif p > pvalue:
            return
        else:
            raise self.failureException(
                msg or "p-value %s, t-test p %s" % (repr(pvalue), repr(p))
            )

    def assertSimilarFreqs(self, observed, expected, pvalue=0.01, msg=None):
        """Fail if observed p is lower than pvalue"""
        if self._suite_pvalue:
            pvalue = self._suite_pvalue

        obs_ravel = ravel(asarray(observed))
        exp_ravel = ravel(asarray(expected))

        m = zeros((2, len(obs_ravel)))
        m[0, :] = obs_ravel
        m[1, :] = exp_ravel

        G, p = G_ind(m)

        if p > pvalue:
            return
        else:
            raise self.failureException(
                msg or "p-value %s, G-test p %s" % (repr(pvalue), repr(p))
            )

    def assertSameObj(self, observed, expected, msg=None):
        """Fail if 'observed is not expected'"""
        try:
            if observed is expected:
                return
        except:
            pass
        raise self.failureException(
            msg
            or "Observed %s is not the same as expected %s"
            % (repr(observed), repr(expected))
        )

    def assertNotSameObj(self, observed, expected, msg=None):
        """Fail if 'observed is expected'"""
        try:
            if observed is not expected:
                return
        except:
            pass
        raise self.failureException(
            msg
            or "Observed %s is the same as expected %s"
            % (repr(observed), repr(expected))
        )
