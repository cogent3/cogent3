#!/usr/bin/env python
"""Extension of the built-in unittest framework for floating-point comparisons.

Specific Extensions:

assertFloatEqual, assertFloatEqualAbs, and assertFloatEqualRel give fine-
grained control over how floating point numbers (or lists thereof) are tested 
for equality. assertContains and assertNotContains give more helpful error 
messages when testing whether an observed item is present or absent in a set 
of possiblities. assertSameItems and assertEqualItems test the items in a list 
for pairwise identity and equality respectively (i.e. the observed and 
expected values must have the same number of each item, though the order can 
differ); assertNotEqualItems verifies that two lists do not contain the same
set of items.
"""
#from contextlib import contextmanager
import numpy; from numpy import testing, array
from unittest import main, TestCase as orig_TestCase, TestSuite, findTestCases
from cogent.util.misc import recursive_flatten

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Sandra Smit",
                    "Zongzhi Liu", "Micah Hamady"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

## SUPPORT2425
#@contextmanager
#def numpy_err(**kw):
#    """a numpy err context manager.

#    **kw: pass to numpy.seterr(all=None, divide=None, over=None, under=None,
#    invalid=None)

#    Example:
#    with numpy_err(divide='raise'):
#        self.assertRaises(FloatingPointError, log, 0)
#    """
#    ori_err = numpy.geterr()
#    numpy.seterr(**kw)
#    try: yield None
#    finally: numpy.seterr(**ori_err)

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
        #wrap around if circular
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
    def _get_values_from_matching_dicts(self, d1, d2):
        """Gets corresponding values from matching dicts"""
        if set(d1) != set (d2):
            return None
        return d1.values(), [d2[k] for k in d1] #might not be in same order


    def errorCheck(self, call, known_errors):
        """Applies function to (data, error) tuples, checking for error
        """
        for (data, error) in known_errors:
            self.assertRaises(error, call, data)

    def valueCheck(self, call, known_values, arg_prefix='', eps=None):
        """Applies function to (data, expected) tuples, treating data as args
        """
        for (data, expected) in known_values:
                observed = eval('call(' + arg_prefix + 'data)')
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
        #do array check first
        #note that we can't use array ops to combine, because we need to check
        #at each element whether the expected is zero to do the test to avoid
        #floating point error.
        #WARNING: numpy iterates over objects that are not regular Python
        #floats/ints, so need to explicitly catch scalar values and prevent
        #cast to array if we want the exact object to print out correctly.
        is_array = False
        if hasattr(obs, 'keys') and hasattr(exp, 'keys'):   #both dicts?
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
                        self.fail("Wrong shape: Got %s, but expected %s" % \
                            (`obs`, `exp`))
                    obs = arr_obs.ravel()
                    exp = arr_exp.ravel()
                    is_array=True
                except (TypeError, ValueError):
                    pass

        for observed, expected in zip(obs, exp):
            #try the cheap comparison first
            if observed == expected:
                continue
            try:
                sum = float(observed + expected)
                diff = float(observed - expected)
                if (sum == 0):
                    if is_array:
                        self.failIf(abs(diff) > abs(eps), \
                            "Got %s, but expected %s (diff was %s)" % \
                            (`arr_obs`, `arr_exp`, `arr_diff`))
                    else:
                        self.failIf(abs(diff) > abs(eps), \
                            "Got %s, but expected %s (diff was %s)" % \
                            (`observed`, `expected`, `diff`))

                else:
                    if is_array:
                        self.failIf(abs(diff/sum) > abs(eps), \
                            "Got %s, but expected %s (diff was %s)" % \
                            (`arr_obs`, `arr_exp`, `arr_diff`))
                    else:
                        self.failIf(abs(diff/sum) > abs(eps), \
                            "Got %s, but expected %s (diff was %s)" % \
                            (`observed`, `expected`, `diff`))
            except (TypeError, ValueError, AttributeError, NotImplementedError):
                self.fail("Got %s, but expected %s" % \
                    (`observed`, `expected`))
    
    def assertFloatEqualAbs(self, obs, exp, eps=1e-6):
        """
        Tests whether two floating point numbers are approximately equal.

        Checks whether the absolute value of (a - b) is within epsilon. Use
        this method when you expect that one of the values should be very
        small, and the other should be zero.
        """
        #do array check first
        #note that we can't use array ops to combine, because we need to check
        #at each element whether the expected is zero to do the test to avoid
        #floating point error.
        if hasattr(obs, 'keys') and hasattr(exp, 'keys'):   #both dicts?
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
                        self.fail("Wrong shape: Got %s, but expected %s" % \
                            (`obs`, `exp`))
                    diff = arr_obs - arr_exp
                    self.failIf(abs(diff).max() > eps, \
                        "Got %s, but expected %s (diff was %s)" % \
                        (`obs`, `exp`, `diff`))
                    return
                except (TypeError, ValueError):
                    pass
        #only get here if array comparison failed
        for observed, expected in zip(obs, exp):
            #cheap comparison first
            if observed == expected:
                continue
            try:
                diff = observed - expected
                self.failIf(abs(diff) > abs(eps),
                        "Got %s, but expected %s (diff was %s)" % \
                        (`observed`, `expected`, `diff`))
            except (TypeError, ValueError, AttributeError, NotImplementedError):
                self.fail("Got %s, but expected %s" % \
                    (`observed`, `expected`))
    
    def assertFloatEqual(self, obs, exp, eps=1e-6, rel_eps=None, \
                         abs_eps=None):
        """Tests whether two floating point numbers are approximately equal.

        If one of the arguments is zero, tests the absolute magnitude of the
        difference; otherwise, tests the relative magnitude.

        Use this method as a reasonable default.
        """
        obs = numpy.asarray(obs, dtype='O')
        exp = numpy.asarray(exp, dtype='O')
        obs = numpy.ravel(obs)
        exp = numpy.ravel(exp)
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
                self.fail("Got %s, but expected %s" % \
                        (`observed`, `expected`))
                                    
    def _is_equal(self, observed, expected):
        """Returns True if observed and expected are equal, False otherwise."""
        #errors to catch: TypeError when obs is None
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
                raise self.failureException, \
                (msg or 'Got %s, but expected %s' % (`observed`, `expected`))
        except (ValueError, TypeError), e:
            #The truth value of an array with more than one element is
            #ambiguous. Use a.any() or a.all()
            #descriptor 'tolist' of 'numpy.generic' object needs an argument
            testing.assert_equal(observed, expected)

    def failIfEqual(self, observed, expected, msg=None):
        """Fail if the two objects are equal as determined by =="""
        try:
            self.failUnlessEqual(observed, expected)
        except self.failureException:
            pass
        else:
            raise self.failureException, \
            (msg or 'Observed %s and expected %s: shouldn\'t test equal'\
                % (`observed`, `expected`))
        
        #following needed to get our version instead of unittest's
    assertEqual = assertEquals = failUnlessEqual

    assertNotEqual = assertNotEquals = failIfEqual

    def assertEqualItems(self, observed, expected, msg=None):
        """Fail if the two items contain unequal elements"""
        obs_items = list(observed)
        exp_items = list(expected)
        if len(obs_items) != len(exp_items):
            raise self.failureException, \
            (msg or 'Observed and expected are different lengths: %s and %s' \
            % (len(obs_items), len(exp_items)))
            
        obs_items.sort()
        exp_items.sort()
        for index, (obs, exp) in enumerate(zip(obs_items, exp_items)):
            if obs != exp:
                raise self.failureException, \
                (msg or 'Observed %s and expected %s at sorted index %s' \
                % (obs, exp, index))

    def assertSameItems(self, observed, expected, msg=None):
        """Fail if the two items contain non-identical elements"""
        obs_items = list(observed)
        exp_items = list(expected)
        if len(obs_items) != len(exp_items):
            raise self.failureException, \
            (msg or 'Observed and expected are different lengths: %s and %s' \
            % (len(obs_items), len(exp_items)))

        obs_ids = [(id(i), i) for i in obs_items]
        exp_ids = [(id(i), i) for i in exp_items]
        obs_ids.sort()
        exp_ids.sort()
        for index, (obs, exp) in enumerate(zip(obs_ids, exp_ids)):
            o_id, o = obs
            e_id, e = exp
            if o_id != e_id:    #i.e. the ids are different
                raise self.failureException, \
                (msg or \
                'Observed %s <%s> and expected %s <%s> at sorted index %s' \
                % (o, o_id, e, e_id, index))

    def assertNotEqualItems(self, observed, expected, msg=None):
        """Fail if the two items contain only equal elements when sorted"""
        try:
            self.assertEqualItems(observed, expected, msg)
        except:
            pass
        else:
            raise self.failureException, \
            (msg or 'Observed %s has same items as %s'%(`observed`, `expected`))

    def assertContains(self, observed, item, msg=None):
        """Fail if item not in observed"""
        try:
            if item in observed:
                return
        except (TypeError, ValueError):
            pass
        raise self.failureException, \
        (msg or 'Item %s not found in %s' % (`item`, `observed`))

    def assertNotContains(self, observed, item, msg=None):
        """Fail if item in observed"""
        try:
            if item not in observed:
                return
        except (TypeError, ValueError):
            return
        raise self.failureException, \
        (msg or 'Item %s should not have been in %s' % (`item`, `observed`))

