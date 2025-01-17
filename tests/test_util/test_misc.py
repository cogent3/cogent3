"""Unit tests for utility functions and classes."""

from copy import copy, deepcopy
from os import remove, rmdir
from unittest import TestCase

import pytest
from numpy import array
from numpy.testing import assert_allclose

import cogent3
from cogent3.util.misc import (
    ClassChecker,
    ConstrainedContainer,
    ConstrainedDict,
    ConstrainedList,
    ConstraintError,
    Delegator,
    DistanceFromMatrix,
    FunctionWrapper,
    MappedDict,
    MappedList,
    NestedSplitter,
    add_lowercase,
    adjusted_gt_minprob,
    adjusted_within_bounds,
    curry,
    docstring_to_summary_rest,
    extend_docstring_from,
    get_independent_coords,
    get_merged_by_value_coords,
    get_merged_overlapping_coords,
    get_object_provenance,
    get_run_start_indices,
    get_setting_from_environ,
    get_true_spans,
    identity,
    is_char,
    is_char_or_noniterable,
    is_iterable,
    iterable,
    list_flatten,
    negate_condition,
    not_list_tuple,
    recursive_flatten,
)


class UtilsTests(TestCase):
    """Tests of individual functions in utils"""

    def setUp(self):
        """ """
        self.files_to_remove = []
        self.dirs_to_remove = []

    def tearDown(self):
        """ """
        list(map(remove, self.files_to_remove))
        list(map(rmdir, self.dirs_to_remove))

    def test_adjusted_gt_minprob(self):
        """correctly adjust a prob vector so all values > minval"""
        vector = [0.8, 0.2, 0.0, 0.0]
        minprob = 1e-3
        got = adjusted_gt_minprob(vector, minprob=minprob)
        assert got.min() > minprob
        minprob = 1e-6
        got = adjusted_gt_minprob(vector, minprob=minprob)
        assert got.min() > minprob
        minprob = 0
        got = adjusted_gt_minprob(vector, minprob=minprob)
        assert got.min() > minprob

    def test_adjusted_probs_2D(self):
        """correctly adjust a 2D array"""
        data = [
            [1.0000, 0.0000, 0.0000, 0.0000],
            [0.0000, 1.0000, 0.0000, 0.0000],
            [0.1250, 0.1250, 0.6250, 0.1250],
            [0.1250, 0.1250, 0.1250, 0.6250],
        ]
        got = adjusted_gt_minprob(data)
        assert_allclose(got, data, atol=1e-5)

    def test_adjusted_within_bounds(self):
        """values correctly adjusted within specified bounds"""
        l, u = 1e-5, 2
        eps = 1e-6
        got = adjusted_within_bounds(l - eps, l, u, eps=eps)
        assert_allclose(got, l)
        got = adjusted_within_bounds(u + eps, l, u, eps=eps)
        assert_allclose(got, u)
        with pytest.raises(ValueError):
            got = adjusted_within_bounds(u + 4, l, u, eps=eps, action="raise")

        with pytest.raises(ValueError):
            got = adjusted_within_bounds(u - 4, l, u, eps=eps, action="raise")

    def test_identity(self):
        """should return same object"""
        foo = [1, "a", lambda x: x]
        exp = id(foo)
        assert id(identity(foo)) == exp

    def test_iterable(self):
        """iterable(x) should return x or [x], always an iterable result"""
        assert iterable("x") == "x"
        assert iterable("") == ""
        assert iterable(3) == [3]
        assert iterable(None) == [None]
        assert iterable({"a": 1}) == {"a": 1}
        assert iterable(["a", "b", "c"]) == ["a", "b", "c"]

    def test_is_iterable(self):
        """is_iterable should return True for iterables"""
        # test str
        assert is_iterable("aa") is True
        # test list
        assert is_iterable([3, "aa"]) is True
        # test Number, expect False
        assert is_iterable(3) is False

    def test_is_char(self):
        """is_char(obj) should return True when obj is a char"""
        assert is_char("a") is True
        assert is_char("ab") is False
        assert is_char("") is True
        assert is_char([3]) is False
        assert is_char(3) is False

    def test_is_char_or_noniterable(self):
        """is_char_or_noniterable should return True or False"""
        assert is_char_or_noniterable("a") is True
        assert is_char_or_noniterable("ab") is False
        assert is_char_or_noniterable(3) is True
        assert is_char_or_noniterable([3]) is False

    def test_recursive_flatten(self):
        """recursive_flatten should remove all nesting from nested sequences"""
        assert recursive_flatten([1, [2, 3], [[4, [5]]]]) == [1, 2, 3, 4, 5]

        # test default behavior on str unpacking
        assert recursive_flatten(["aa", [8, "cc", "dd"], ["ee", ["ff", "gg"]]]) == [
            "a",
            "a",
            8,
            "c",
            "c",
            "d",
            "d",
            "e",
            "e",
            "f",
            "f",
            "g",
            "g",
        ]

    def test_not_list_tuple(self):
        """not_list_tuple(obj) should return False when obj is list or tuple"""
        assert not_list_tuple([8, 3]) is False
        assert not_list_tuple((8, 3)) is False
        assert not_list_tuple("34") is True

    def test_list_flatten(self):
        """list_flatten should remove all nesting, str is untouched"""
        assert list_flatten(["aa", [8, "cc", "dd"], ["ee", ["ff", "gg"]]]) == [
            "aa",
            8,
            "cc",
            "dd",
            "ee",
            "ff",
            "gg",
        ]

    def test_recursive_flatten_max_depth(self):
        """recursive_flatten should not remover more than max_depth levels"""
        assert recursive_flatten([1, [2, 3], [[4, [5]]]]) == [1, 2, 3, 4, 5]
        assert recursive_flatten([1, [2, 3], [[4, [5]]]], 0) == [1, [2, 3], [[4, [5]]]]
        assert recursive_flatten([1, [2, 3], [[4, [5]]]], 1) == [1, 2, 3, [4, [5]]]
        assert recursive_flatten([1, [2, 3], [[4, [5]]]], 2) == [1, 2, 3, 4, [5]]
        assert recursive_flatten([1, [2, 3], [[4, [5]]]], 3) == [1, 2, 3, 4, 5]
        assert recursive_flatten([1, [2, 3], [[4, [5]]]], 5000) == [1, 2, 3, 4, 5]

    def test_add_lowercase(self):
        """add_lowercase should add lowercase version of each key w/ same val"""
        d = {
            "a": 1,
            "b": "test",
            "A": 5,
            "C": 123,
            "D": [],
            "AbC": "XyZ",
            None: "3",
            "$": "abc",
            145: "5",
        }
        add_lowercase(d)
        assert d["d"] is d["D"]
        d["D"].append(3)
        assert d["D"] == [3]
        assert d["d"] == [3]
        assert d["a"] != d["A"]
        assert d == {
            "a": 1,
            "b": "test",
            "A": 5,
            "C": 123,
            "c": 123,
            "D": [3],
            "d": [3],
            "AbC": "XyZ",
            "abc": "xyz",
            None: "3",
            "$": "abc",
            145: "5",
        }

        # should work with strings
        d = "ABC"
        assert add_lowercase(d) == "ABCabc"
        # should work with tuples
        d = tuple("ABC")
        assert add_lowercase(d) == tuple("ABCabc")
        # should work with lists
        d = list("ABC")
        assert add_lowercase(d) == list("ABCabc")
        # should work with sets
        d = set("ABC")
        assert add_lowercase(d) == set("ABCabc")
        # ...even frozensets
        d = frozenset("ABC")
        assert add_lowercase(d) == frozenset("ABCabc")

    def test_add_lowercase_tuple(self):
        """add_lowercase should deal with tuples correctly"""
        d = {("A", "B"): "C", ("D", "e"): "F", ("b", "c"): "H"}
        add_lowercase(d)
        assert d == {
            ("A", "B"): "C",
            ("a", "b"): "c",
            ("D", "e"): "F",
            ("d", "e"): "f",
            ("b", "c"): "H",
        }

    def test_DistanceFromMatrix(self):
        """DistanceFromMatrix should return correct elements"""
        m = {"a": {"3": 4, 6: 1}, "b": {"3": 5, "6": 2}}
        d = DistanceFromMatrix(m)
        assert d("a", "3") == 4
        assert d("a", 6) == 1
        assert d("b", "3") == 5
        assert d("b", "6") == 2
        self.assertRaises(KeyError, d, "c", 1)
        self.assertRaises(KeyError, d, "b", 3)

    def test_independent_spans(self):
        """get_independent_coords returns truly non-overlapping (decorated) spans"""
        # single span is returned
        data = [(0, 20, "a")]
        got = get_independent_coords(data)
        assert got == data

        # multiple non-overlapping
        data = [(20, 30, "a"), (35, 40, "b"), (65, 75, "c")]
        got = get_independent_coords(data)
        assert got == data

        # over-lapping first/second returns first occurrence by default
        data = [(20, 30, "a"), (25, 40, "b"), (65, 75, "c")]
        got = get_independent_coords(data)
        assert got == [(20, 30, "a"), (65, 75, "c")]
        # but randomly the first or second if random_tie_breaker is chosen
        got = get_independent_coords(data, random_tie_breaker=True)
        assert got in ([(20, 30, "a"), (65, 75, "c")], [(25, 40, "b"), (65, 75, "c")])

        # over-lapping second/last returns first occurrence by default
        data = [(20, 30, "a"), (30, 60, "b"), (50, 75, "c")]
        got = get_independent_coords(data)
        assert got == [(20, 30, "a"), (30, 60, "b")]
        # but randomly the first or second if random_tie_breaker is chosen
        got = get_independent_coords(data, random_tie_breaker=True)
        assert got in ([(20, 30, "a"), (50, 75, "c")], [(20, 30, "a"), (30, 60, "b")])

        # over-lapping middle returns first occurrence by default
        data = [(20, 24, "a"), (25, 40, "b"), (30, 35, "c"), (65, 75, "d")]
        got = get_independent_coords(data)
        assert got == [(20, 24, "a"), (25, 40, "b"), (65, 75, "d")]

        # but randomly the first or second if random_tie_breaker is chosen
        got = get_independent_coords(data, random_tie_breaker=True)
        assert got in (
            [(20, 24, "a"), (25, 40, "b"), (65, 75, "d")],
            [(20, 24, "a"), (30, 35, "c"), (65, 75, "d")],
        )

    def test_get_merged_spans(self):
        """tests merger of overlapping spans"""
        sample = [[0, 10], [12, 15], [13, 16], [18, 25], [19, 20]]
        result = get_merged_overlapping_coords(sample)
        expect = [[0, 10], [12, 16], [18, 25]]
        assert result == expect
        sample = [[0, 10], [5, 9], [12, 16], [18, 20], [19, 25]]
        result = get_merged_overlapping_coords(sample)
        expect = [[0, 10], [12, 16], [18, 25]]
        assert result == expect
        # test with tuples
        sample = tuple(map(tuple, sample))
        result = get_merged_overlapping_coords(sample)
        expect = [[0, 10], [12, 16], [18, 25]]
        assert result == expect

    def test_get_run_start_indices(self):
        """return indices corresponding to start of a run of identical values"""
        #       0  1  2  3  4  5  6  7
        data = [1, 2, 3, 3, 3, 4, 4, 5]
        expect = [[0, 1], [1, 2], [2, 3], [5, 4], [7, 5]]
        got = get_run_start_indices(data)
        assert list(got) == expect

        # raise an exception if try and provide a converter and num digits
        def wrap_gen():  # need to wrap generator so we can actually test this
            gen = get_run_start_indices(data, digits=1, converter_func=lambda x: x)

            def call():
                for _v in gen:
                    pass

            return call

        self.assertRaises(AssertionError, wrap_gen())

    def test_merged_by_value_spans(self):
        """correctly merge adjacent spans with the same value"""
        # initial values same
        data = [[20, 21, 0], [21, 22, 0], [22, 23, 1], [23, 24, 0]]
        assert get_merged_by_value_coords(data) == [
            [20, 22, 0],
            [22, 23, 1],
            [23, 24, 0],
        ]

        # middle values same
        data = [[20, 21, 0], [21, 22, 1], [22, 23, 1], [23, 24, 0]]
        assert get_merged_by_value_coords(data) == [
            [20, 21, 0],
            [21, 23, 1],
            [23, 24, 0],
        ]

        # last values same
        data = [[20, 21, 0], [21, 22, 1], [22, 23, 0], [23, 24, 0]]
        assert get_merged_by_value_coords(data) == [
            [20, 21, 0],
            [21, 22, 1],
            [22, 24, 0],
        ]

        # all unique values
        data = [[20, 21, 0], [21, 22, 1], [22, 23, 2], [23, 24, 0]]
        assert get_merged_by_value_coords(data) == [
            [20, 21, 0],
            [21, 22, 1],
            [22, 23, 2],
            [23, 24, 0],
        ]

        # all values same
        data = [[20, 21, 0], [21, 22, 0], [22, 23, 0], [23, 24, 0]]
        assert get_merged_by_value_coords(data) == [[20, 24, 0]]

        # all unique values to 2nd decimal
        data = [[20, 21, 0.11], [21, 22, 0.12], [22, 23, 0.13], [23, 24, 0.14]]
        assert get_merged_by_value_coords(data) == [
            [20, 21, 0.11],
            [21, 22, 0.12],
            [22, 23, 0.13],
            [23, 24, 0.14],
        ]

        # all values same at 1st decimal
        data = [[20, 21, 0.11], [21, 22, 0.12], [22, 23, 0.13], [23, 24, 0.14]]
        assert get_merged_by_value_coords(data, digits=1) == [[20, 24, 0.1]]

    def test_get_object_provenance(self):
        """correctly deduce object provenance"""
        result = get_object_provenance("abncd")
        assert result == "str"

        DNA = cogent3.get_moltype("dna")
        got = get_object_provenance(DNA)
        assert got == f"{DNA.__module__}.MolType"

        sm = cogent3.get_model("HKY85")
        got = get_object_provenance(sm)
        assert got == "cogent3.evolve.substitution_model.TimeReversibleNucleotide"

        # handle a type
        instance = cogent3.make_unaligned_seqs(
            data={"a": "ACG", "b": "GGG"},
            moltype="dna",
        )
        instance_prov = get_object_provenance(instance)
        assert instance_prov == f"{instance.__module__}.SequenceCollection"
        type_prov = get_object_provenance(type(instance))
        assert instance_prov == type_prov

    def test_get_object_provenance_builtins(self):
        """allow identifying builtins too"""
        from gzip import GzipFile, compress

        obj_prov = get_object_provenance(compress)

        assert obj_prov == "gzip.compress"

        obj_prov = get_object_provenance(GzipFile)
        assert obj_prov == "gzip.GzipFile"

        d = {"a": 23, "b": 1}
        obj_prov = get_object_provenance(d)
        assert obj_prov == "dict"

        obj_prov = get_object_provenance(dict)
        assert obj_prov == "dict"

    def test_NestedSplitter(self):
        """NestedSplitter should make a function which return expected list"""
        # test delimiters, constructor, filter_
        line = "ii=0; oo= 9, 6 5;  ; xx=  8;  "  # noqa
        cmds = [
            "NestedSplitter(';=,')(line)",
            "NestedSplitter([';', '=', ','])(line)",
            "NestedSplitter([(';'), '=', ','], constructor=None)(line)",
            "NestedSplitter([(';'), '=', ','], filter_=None)(line)",
            "NestedSplitter([(';',1), '=', ','])(line)",
            "NestedSplitter([(';',-1), '=', ','])(line)",
        ]
        results = [
            [["ii", "0"], ["oo", ["9", "6 5"]], "", ["xx", "8"], ""],
            [["ii", "0"], ["oo", ["9", "6 5"]], "", ["xx", "8"], ""],
            [["ii", "0"], [" oo", [" 9", " 6 5"]], "  ", [" xx", "  8"], "  "],
            [["ii", "0"], ["oo", ["9", "6 5"]], ["xx", "8"]],
            [["ii", "0"], ["oo", ["9", "6 5;  ; xx"], "8;"]],
            [["ii", "0; oo", ["9", "6 5;  ; xx"], "8"], ""],
        ]
        for cmd, result in zip(cmds, results, strict=False):
            assert eval(cmd) == result

        # test uncontinous level of delimiters
        test = "a; b,c; d,e:f; g:h;"  # g:h should get [[g,h]] instead of [g,h]
        assert NestedSplitter(";,:")(test) == [
            "a",
            ["b", "c"],
            ["d", ["e", "f"]],
            [["g", "h"]],
            "",
        ]

        # test empty
        assert NestedSplitter(";,:")("") == [""]
        assert NestedSplitter(";,:")("  ") == [""]
        assert NestedSplitter(";,:", filter_=None)(" ;, :") == [[[]]]

    def test_curry(self):
        """curry should generate the function with parameters setted"""
        curry_test = curry(lambda x, y: x == y, 5)
        knowns = ((3, False), (9, False), (5, True))
        for arg2, result in knowns:
            assert curry_test(arg2) == result

    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_get_setting_from_environ(self):
        """correctly recovers environment variables"""
        import os

        def make_env_setting(d):
            return ",".join([f"{k}={v}" for k, v in d.items()])

        env_name = "DUMMY_SETTING"
        os.environ.pop(env_name, None)
        setting = {"num_pos": 2, "num_seq": 4, "name": "blah"}
        single_setting = {"num_pos": 2}
        correct_names_types = {"num_pos": int, "num_seq": int, "name": str}
        incorrect_names_types = {"num_pos": int, "num_seq": int, "name": float}

        for stng in (setting, single_setting):
            os.environ[env_name] = make_env_setting(stng)
            got = get_setting_from_environ(env_name, correct_names_types)
            for key in got:
                assert got[key] == setting[key]

        os.environ[env_name] = make_env_setting(setting)
        got = get_setting_from_environ(env_name, incorrect_names_types)
        assert "name" not in got
        for key in got:
            assert got[key] == setting[key]

        # malformed env setting
        os.environ[env_name] = make_env_setting(setting).replace("=", "")
        got = get_setting_from_environ(env_name, correct_names_types)
        assert got == {}

        os.environ.pop(env_name, None)


class _my_dict(dict):
    """Used for testing subclass behavior of ClassChecker"""


class ClassCheckerTests(TestCase):
    """Unit tests for the ClassChecker class."""

    def setUp(self):
        """define a few standard checkers"""
        self.strcheck = ClassChecker(str)
        self.intcheck = ClassChecker(int, int)
        self.numcheck = ClassChecker(float, int, int)
        self.emptycheck = ClassChecker()
        self.dictcheck = ClassChecker(dict)
        self.mydictcheck = ClassChecker(_my_dict)

    def test_init_good(self):
        """ClassChecker should init OK when initialized with classes"""
        assert self.strcheck.Classes == [str]
        assert self.numcheck.Classes == [float, int, int]
        assert self.emptycheck.Classes == []

    def test_init_bad(self):
        """ClassChecker should raise TypeError if initialized with non-class"""
        self.assertRaises(TypeError, ClassChecker, "x")
        self.assertRaises(TypeError, ClassChecker, str, None)

    def test_contains(self):
        """ClassChecker should return True only if given instance of class"""
        assert self.strcheck.__contains__("3") is True
        assert self.strcheck.__contains__("ahsdahisad") is True
        assert self.strcheck.__contains__(3) is False
        assert self.strcheck.__contains__({3: "c"}) is False

        assert self.intcheck.__contains__("ahsdahisad") is False
        assert self.intcheck.__contains__("3") is False
        assert self.intcheck.__contains__(3.0) is False
        assert self.intcheck.__contains__(3) is True
        assert self.intcheck.__contains__(4**60) is True
        assert self.intcheck.__contains__(4**60 * -1) is True

        d = _my_dict()
        assert self.dictcheck.__contains__(d) is True
        assert self.dictcheck.__contains__({"d": 1}) is True
        assert self.mydictcheck.__contains__(d) is True
        assert self.mydictcheck.__contains__({"d": 1}) is False

        assert self.emptycheck.__contains__("d") is False

        assert self.numcheck.__contains__(3) is True
        assert self.numcheck.__contains__(3.0) is True
        assert self.numcheck.__contains__(-3) is True
        assert self.numcheck.__contains__(-3.0) is True
        assert self.numcheck.__contains__(3e-300) is True
        assert self.numcheck.__contains__(0) is True
        assert self.numcheck.__contains__(4**1000) is True
        assert self.numcheck.__contains__("4**1000") is False

    def test_str(self):
        """ClassChecker str should be the same as str(self.Classes)"""
        for c in [
            self.strcheck,
            self.intcheck,
            self.numcheck,
            self.emptycheck,
            self.dictcheck,
            self.mydictcheck,
        ]:
            assert str(c) == str(c.Classes)

    def test_copy(self):
        """copy.copy should work correctly on ClassChecker"""
        c = copy(self.strcheck)
        assert c is not self.strcheck
        assert "3" in c
        assert 3 not in c
        assert c.Classes is self.strcheck.Classes

    def test_deepcopy(self):
        """copy.deepcopy should work correctly on ClassChecker"""
        c = deepcopy(self.strcheck)
        assert c is not self.strcheck
        assert "3" in c
        assert 3 not in c
        assert c.Classes is not self.strcheck.Classes


class modifiable_string(str):
    """Mutable class to allow arbitrary attributes to be set"""


class _list_and_string(list, Delegator):
    """Trivial class to demonstrate Delegator."""

    def __init__(self, items, string):
        Delegator.__init__(self, string)
        self.NormalAttribute = "default"
        self._x = None
        self._constant = "c"
        for i in items:
            self.append(i)

    def _get_rand_property(self):
        return self._x

    def _set_rand_property(self, value):
        self._x = value

    prop = property(_get_rand_property, _set_rand_property)

    def _get_constant_property(self):
        return self._constant

    constant = property(_get_constant_property)


class DelegatorTests(TestCase):
    """Verify that Delegator works with attributes and properties."""

    def test_init(self):
        """Delegator should init OK when data supplied"""
        _list_and_string([1, 2, 3], "abc")
        self.assertRaises(TypeError, _list_and_string, [123])

    def test_getattr(self):
        """Delegator should find attributes in correct places"""
        ls = _list_and_string([1, 2, 3], "abcd")
        # behavior as list
        assert len(ls) == 3
        assert ls[0] == 1
        ls.reverse()
        assert ls == [3, 2, 1]
        # behavior as string
        assert ls.upper() == "ABCD"
        assert len(ls.upper()) == 4
        assert ls.replace("a", "x") == "xbcd"
        # behavior of normal attributes
        assert ls.NormalAttribute == "default"
        # behavior of properties
        assert ls.prop is None
        assert ls.constant == "c"
        # shouldn't be allowed to get unknown properties
        self.assertRaises(AttributeError, getattr, ls, "xyz")
        # if the unknown property can be set in the forwarder, do it there
        flex = modifiable_string("abcd")
        ls_flex = _list_and_string([1, 2, 3], flex)
        ls_flex.blah = "zxc"
        assert ls_flex.blah == "zxc"
        assert flex.blah == "zxc"
        # should get AttributeError if changing a read-only property
        self.assertRaises(AttributeError, setattr, ls, "constant", "xyz")

    def test_setattr(self):
        """Delegator should set attributes in correct places"""
        ls = _list_and_string([1, 2, 3], "abcd")
        # ability to create a new attribute
        ls.xyz = 3
        assert ls.xyz == 3
        # modify a normal attribute
        ls.NormalAttribute = "changed"
        assert ls.NormalAttribute == "changed"
        # modify a read/write property
        ls.prop = "xyz"
        assert ls.prop == "xyz"

    def test_copy(self):
        """copy.copy should work correctly on Delegator"""
        l = ["a"]
        d = Delegator(l)
        c = copy(d)
        assert c is not d
        assert c._handler is d._handler

    def test_deepcopy(self):
        """copy.deepcopy should work correctly on Delegator"""
        l = ["a"]
        d = Delegator(l)
        c = deepcopy(d)
        assert c is not d
        assert c._handler is not d._handler
        assert c._handler == d._handler


class FunctionWrapperTests(TestCase):
    """Tests of the FunctionWrapper class"""

    def test_init(self):
        """FunctionWrapper should initialize with any callable"""
        f = FunctionWrapper(str)
        g = FunctionWrapper(id)
        h = FunctionWrapper(iterable)
        x = 3
        assert f(x) == "3"
        assert g(x) == id(x)
        assert h(x) == [3]

    def test_copy(self):
        """copy should work for FunctionWrapper objects"""
        f = FunctionWrapper(str)
        c = copy(f)
        assert c is not f
        assert c.Function is f.Function

    # NOTE: deepcopy does not work for FunctionWrapper objects because you
    # can't copy a function.


class _simple_container:
    """example of a container to constrain"""

    def __init__(self, data):
        self._data = list(data)

    def __getitem__(self, item):
        return self._data.__getitem__(item)


class _constrained_simple_container(_simple_container, ConstrainedContainer):
    """constrained version of _simple_container"""

    def __init__(self, data):
        _simple_container.__init__(self, data)
        ConstrainedContainer.__init__(self, None)


class ConstrainedContainerTests(TestCase):
    """Tests of the generic ConstrainedContainer interface."""

    def setUp(self):
        """Make a couple of standard containers"""
        self.alphabet = _constrained_simple_container("abc")
        self.numbers = _constrained_simple_container([1, 2, 3])
        self.alphacontainer = "abcdef"
        self.numbercontainer = ClassChecker(int)

    def test_matchesConstraint(self):
        """ConstrainedContainer matchesConstraint should return true if items ok"""
        assert self.alphabet.matches_constraint(self.alphacontainer) is True
        assert self.alphabet.matches_constraint(self.numbercontainer) is False
        assert self.numbers.matches_constraint(self.alphacontainer) is False
        assert self.numbers.matches_constraint(self.numbercontainer) is True

    def test_other_is_valid(self):
        """ConstrainedContainer should use constraint for checking other"""
        assert self.alphabet.other_is_valid("12d8jc") is True
        self.alphabet.constraint = self.alphacontainer
        assert self.alphabet.other_is_valid("12d8jc") is False
        self.alphabet.constraint = list("abcdefghijkl12345678")
        assert self.alphabet.other_is_valid("12d8jc") is True
        assert self.alphabet.other_is_valid("z") is False

    def test_item_is_valid(self):
        """ConstrainedContainer should use constraint for checking item"""
        assert self.alphabet.item_is_valid(3) is True
        self.alphabet.constraint = self.alphacontainer
        assert self.alphabet.item_is_valid(3) is False
        assert self.alphabet.item_is_valid("a") is True

    def test_sequence_is_valid(self):
        """ConstrainedContainer should use constraint for checking sequence"""
        assert self.alphabet.sequence_is_valid("12d8jc") is True
        self.alphabet.constraint = self.alphacontainer
        assert self.alphabet.sequence_is_valid("12d8jc") is False
        self.alphabet.constraint = list("abcdefghijkl12345678")
        assert self.alphabet.sequence_is_valid("12d8jc") is True
        assert self.alphabet.sequence_is_valid("z") is False

    def test_Constraint(self):
        """ConstrainedContainer should only allow valid constraints to be set"""
        try:
            self.alphabet.constraint = self.numbers
        except ConstraintError:
            pass
        else:
            msg = "Failed to raise ConstraintError with invalid constraint."
            raise AssertionError(
                msg,
            )
        self.alphabet.constraint = "abcdefghi"
        self.alphabet.constraint = ["a", "b", "c", 1, 2, 3]
        self.numbers.constraint = list(range(20))
        self.numbers.constraint = range(20)
        self.numbers.constraint = [5, 1, 3, 7, 2]
        self.numbers.constraint = {1: "a", 2: "b", 3: "c"}
        self.assertRaises(ConstraintError, setattr, self.numbers, "constraint", "1")


class ConstrainedListTests(TestCase):
    """Tests that bad data can't sneak into ConstrainedLists."""

    def test_init_good_data(self):
        """ConstrainedList should init OK if list matches constraint"""
        assert ConstrainedList("abc", "abcd") == list("abc")
        assert ConstrainedList("", "abcd") == list("")
        items = [1, 2, 3.2234, (["a"], ["b"]), list("xyz")]
        # should accept anything str() does if no constraint is passed
        assert ConstrainedList(items) == items
        assert ConstrainedList(items, None) == items
        assert ConstrainedList("12345") == list("12345")
        # check that list is formatted correctly and chars are all there
        test_list = list("12345")
        assert ConstrainedList(test_list, "12345") == test_list

    def test_init_bad_data(self):
        """ConstrainedList should fail init with items not in constraint"""
        self.assertRaises(ConstraintError, ConstrainedList, "1234", "123")
        self.assertRaises(ConstraintError, ConstrainedList, [1, 2, 3], ["1", "2", "3"])

    def test_add_prevents_bad_data(self):
        """ConstrainedList should allow addition only of compliant data"""
        a = ConstrainedList("123", "12345")
        b = ConstrainedList("444", "4")
        c = ConstrainedList("45", "12345")
        d = ConstrainedList("x")
        assert a + b == list("123444")
        assert a + c == list("12345")
        self.assertRaises(ConstraintError, b.__add__, c)
        self.assertRaises(ConstraintError, c.__add__, d)
        # should be OK if constraint removed
        b.constraint = None
        assert b + c == list("44445")
        assert b + d == list("444x")
        # should fail if we add the constraint back
        b.constraint = {"4": 1, 5: 2}
        self.assertRaises(ConstraintError, b.__add__, c)

    def test_iadd_prevents_bad_data(self):
        """ConstrainedList should allow in-place addition only of compliant data"""
        a = ConstrainedList("12", "123")
        a += "2"
        assert a == list("122")
        assert a.constraint == "123"
        self.assertRaises(ConstraintError, a.__iadd__, "4")

    def test_imul(self):
        """ConstrainedList imul should preserve constraint"""
        a = ConstrainedList("12", "123")
        a *= 3
        assert a == list("121212")
        assert a.constraint == "123"

    def test_mul(self):
        """ConstrainedList mul should preserve constraint"""
        a = ConstrainedList("12", "123")
        b = a * 3
        assert b == list("121212")
        assert b.constraint == "123"

    def test_rmul(self):
        """ConstrainedList rmul should preserve constraint"""
        a = ConstrainedList("12", "123")
        b = 3 * a
        assert b == list("121212")
        assert b.constraint == "123"

    def test_setitem(self):
        """ConstrainedList setitem should work only if item in constraint"""
        a = ConstrainedList("12", "123")
        a[0] = "3"
        assert a == list("32")
        self.assertRaises(ConstraintError, a.__setitem__, 0, 3)
        a = ConstrainedList("1" * 20, "123")
        self.assertRaises(ConstraintError, a.__setitem__, slice(0, 1, 1), [3])
        self.assertRaises(ConstraintError, a.__setitem__, slice(0, 1, 1), ["111"])
        a[2:9:2] = "3333"
        assert a == list("11313131311111111111")

    def test_append(self):
        """ConstrainedList append should work only if item in constraint"""
        a = ConstrainedList("12", "123")
        a.append("3")
        assert a == list("123")
        self.assertRaises(ConstraintError, a.append, 3)

    def test_extend(self):
        """ConstrainedList extend should work only if all items in constraint"""
        a = ConstrainedList("12", "123")
        a.extend("321")
        assert a == list("12321")
        self.assertRaises(ConstraintError, a.extend, ["1", "2", 3])

    def test_insert(self):
        """ConstrainedList insert should work only if item in constraint"""
        a = ConstrainedList("12", "123")
        a.insert(0, "2")
        assert a == list("212")
        self.assertRaises(ConstraintError, a.insert, 0, [2])

    def test_getslice(self):
        """ConstrainedList getslice should remember constraint"""
        a = ConstrainedList("123333", "12345")
        b = a[2:4]
        assert b == list("33")
        assert b.constraint == "12345"

    def test_setslice(self):
        """ConstrainedList setslice should fail if slice has invalid chars"""
        a = ConstrainedList("123333", "12345")
        a[2:4] = ["2", "2"]
        assert a == list("122233")
        self.assertRaises(ConstraintError, a.__setslice__, 2, 4, [2, 2])
        a[:] = []
        assert a == []
        assert a.constraint == "12345"

    def test_setitem_masks(self):
        """ConstrainedList setitem with masks should transform input"""
        a = ConstrainedList("12333", list(range(5)), lambda x: int(x) + 1)
        assert a == [2, 3, 4, 4, 4]
        self.assertRaises(ConstraintError, a.append, 4)
        b = a[1:3]
        assert b.mask is a.mask
        assert "1" not in a
        assert "2" not in a
        assert 2 in a
        assert "x" not in a


class MappedListTests(TestCase):
    """MappedList should behave like ConstrainedList, but map items."""

    def test_setitem_masks(self):
        """MappedList setitem with masks should transform input"""
        a = MappedList("12333", list(range(5)), lambda x: int(x) + 1)
        assert a == [2, 3, 4, 4, 4]
        self.assertRaises(ConstraintError, a.append, 4)
        b = a[1:3]
        assert b.mask is a.mask
        assert "1" in a
        assert "x" not in a


class ConstrainedDictTests(TestCase):
    """Tests that bad data can't sneak into ConstrainedDicts."""

    def test_init_good_data(self):
        """ConstrainedDict should init OK if list matches constraint"""
        assert ConstrainedDict(dict.fromkeys("abc"), "abcd") == dict.fromkeys("abc")
        assert ConstrainedDict("", "abcd") == dict("")
        items = [1, 2, 3.2234, tuple("xyz")]
        # should accept anything dict() does if no constraint is passed
        assert ConstrainedDict(dict.fromkeys(items)) == dict.fromkeys(items)
        assert ConstrainedDict(dict.fromkeys(items), None) == dict.fromkeys(items)
        assert ConstrainedDict([(x, 1) for x in "12345"]) == dict.fromkeys("12345", 1)
        # check that list is formatted correctly and chars are all there
        test_dict = dict.fromkeys("12345")
        assert ConstrainedDict(test_dict, "12345") == test_dict

    def test_init_sequence(self):
        """ConstrainedDict should init from sequence, unlike normal dict"""
        assert ConstrainedDict("abcda") == {"a": 2, "b": 1, "c": 1, "d": 1}

    def test_init_bad_data(self):
        """ConstrainedDict should fail init with items not in constraint"""
        self.assertRaises(
            ConstraintError,
            ConstrainedDict,
            dict.fromkeys("1234"),
            "123",
        )
        self.assertRaises(
            ConstraintError,
            ConstrainedDict,
            dict.fromkeys([1, 2, 3]),
            ["1", "2", "3"],
        )

    def test_setitem(self):
        """ConstrainedDict setitem should work only if key in constraint"""
        a = ConstrainedDict(dict.fromkeys("12"), "123")
        a["1"] = "3"
        assert a == {"1": "3", "2": None}
        self.assertRaises(ConstraintError, a.__setitem__, 1, "3")

    def test_copy(self):
        """ConstrainedDict copy should retain constraint"""
        a = ConstrainedDict(dict.fromkeys("12"), "123")
        b = a.copy()
        assert a.constraint == b.constraint
        self.assertRaises(ConstraintError, a.__setitem__, 1, "3")
        self.assertRaises(ConstraintError, b.__setitem__, 1, "3")

    def test_fromkeys(self):
        """ConstrainedDict instance fromkeys should retain constraint"""
        a = ConstrainedDict(dict.fromkeys("12"), "123")
        b = a.fromkeys("23")
        assert a.constraint == b.constraint
        self.assertRaises(ConstraintError, a.__setitem__, 1, "3")
        self.assertRaises(ConstraintError, b.__setitem__, 1, "3")
        b["2"] = 5
        assert b == {"2": 5, "3": None}

    def test_setdefault(self):
        """ConstrainedDict setdefault shouldn't allow bad keys"""
        a = ConstrainedDict({"1": None, "2": "xyz"}, "123")
        assert a.setdefault("2", None) == "xyz"
        assert a.setdefault("1", None) is None
        self.assertRaises(ConstraintError, a.setdefault, "x", 3)
        a.setdefault("3", 12345)
        assert a == {"1": None, "2": "xyz", "3": 12345}

    def test_update(self):
        """ConstrainedDict should allow update only of compliant data"""
        a = ConstrainedDict(dict.fromkeys("123"), "12345")
        b = ConstrainedDict(dict.fromkeys("444"), "4")
        c = ConstrainedDict(dict.fromkeys("45"), "12345")
        d = ConstrainedDict([["x", "y"]])
        a.update(b)
        assert a == dict.fromkeys("1234")
        a.update(c)
        assert a == dict.fromkeys("12345")
        self.assertRaises(ConstraintError, b.update, c)
        self.assertRaises(ConstraintError, c.update, d)
        # should be OK if constraint removed
        b.constraint = None
        b.update(c)
        assert b == dict.fromkeys("45")
        b.update(d)
        assert b == {"4": None, "5": None, "x": "y"}
        # should fail if we add the constraint back
        b.constraint = {"4": 1, 5: 2, "5": 1, "x": 1}
        self.assertRaises(ConstraintError, b.update, {4: 1})
        b.update({5: 1})
        assert b == {"4": None, "5": None, "x": "y", 5: 1}

    def test_setitem_masks(self):
        """ConstrainedDict setitem should work only if key in constraint"""
        key_mask = str

        def val_mask(x):
            return int(x) + 3

        d = ConstrainedDict({1: 4, 2: 6}, "123", key_mask, val_mask)
        d[1] = "456"
        assert d == {"1": 459, "2": 9}
        d["1"] = 234
        assert d == {"1": 237, "2": 9}
        self.assertRaises(ConstraintError, d.__setitem__, 4, "3")
        e = d.copy()
        assert e.mask is d.mask
        assert "1" in d
        assert 1 not in d


class MappedDictTests(TestCase):
    """MappedDict should work like ConstrainedDict, but map keys."""

    def test_setitem_masks(self):
        """MappedDict setitem should work only if key in constraint"""
        key_mask = str

        def val_mask(x):
            return int(x) + 3

        d = MappedDict({1: 4, 2: 6}, "123", key_mask, val_mask)
        d[1] = "456"
        assert d == {"1": 459, "2": 9}
        d["1"] = 234
        assert d == {"1": 237, "2": 9}
        self.assertRaises(ConstraintError, d.__setitem__, 4, "3")
        e = d.copy()
        assert e.mask is d.mask
        assert "1" in d
        assert 1 in d
        assert 1 not in list(d.keys())
        assert "x" not in list(d.keys())

    def test_getitem(self):
        """MappedDict getitem should automatically map key."""
        key_mask = str
        d = MappedDict({}, "123", key_mask)
        assert d == {}
        d["1"] = 5
        assert d == {"1": 5}
        assert d[1] == 5

    def test_get(self):
        """MappedDict get should automatically map key."""
        key_mask = str
        d = MappedDict({}, "123", key_mask)
        assert d == {}
        d["1"] = 5
        assert d == {"1": 5}
        assert d.get(1, "x") == 5
        assert d.get(5, "x") == "x"

    def test_has_key(self):
        """MappedDict has_key should automatically map key."""
        key_mask = str
        d = MappedDict({}, "123", key_mask)
        assert d == {}
        d["1"] = 5
        assert "1" in d
        assert 1 in d
        assert "5" not in d


def f():
    """This is a function docstring."""


@extend_docstring_from(f)
def foo_append():
    """I am foo."""


@extend_docstring_from(f)
def foo_mirror():
    pass


@extend_docstring_from(f, pre=True)
def foo_prepend():
    """I am foo."""


class ExtendDocstringTests(TestCase):
    @extend_docstring_from(f)
    def foo_append(self):
        """I am foo."""

    @extend_docstring_from(f)
    def foo_mirror(self):
        pass

    @extend_docstring_from(f, pre=True)
    def foo_prepend(self):
        """I am foo."""

    class TemplateClass:
        """This is a class docstring."""

    @extend_docstring_from(TemplateClass)
    class FooAppend:
        """I am foo."""

    @extend_docstring_from(TemplateClass)
    class FooMirror:
        pass

    @extend_docstring_from(TemplateClass, pre=True)
    class FooPrepend:
        """I am foo."""

    def test_function_append(self):
        assert foo_append.__doc__ == "This is a function docstring.\nI am foo."

    def test_function_mirror(self):
        assert foo_mirror.__doc__ == "This is a function docstring.\n"

    def test_function_prepend(self):
        assert foo_prepend.__doc__ == "I am foo.\nThis is a function docstring."

    def test_method_append(self):
        assert self.foo_append.__doc__ == "This is a function docstring.\nI am foo."

    def test_method_mirror(self):
        assert self.foo_mirror.__doc__ == "This is a function docstring.\n"

    def test_method_prepend(self):
        assert self.foo_prepend.__doc__ == "I am foo.\nThis is a function docstring."

    def test_class_append(self):
        assert self.FooAppend.__doc__ == "This is a class docstring.\nI am foo."

    def test_class_mirror(self):
        assert self.FooMirror.__doc__ == "This is a class docstring.\n"

    def test_class_prepend(self):
        assert self.FooPrepend.__doc__ == "I am foo.\nThis is a class docstring."


def test_not_in_jupyter():
    from cogent3.util.misc import in_jupyter

    assert not in_jupyter()


def test_is_in_jupyter():
    # an ugly hack, the in_jupyter function relies entirely on whether a
    # get_ipython variable exists in the name space
    import cogent3.util.misc as module
    from cogent3.util.misc import in_jupyter

    module.get_ipython = lambda x: x
    assert in_jupyter()
    del module.get_ipython


def foo1():
    """some text"""


def foo2():
    """some text

    Notes
    -----
    body
    """


def foo3():
    """
    Notes
    -----
    body
    """


def foo4(): ...


_sum_expect = "some text"
_body_expect = ["Notes", "-----", "body"]


@pytest.mark.parametrize(
    ("foo", "sum_exp", "body_exp"),
    [
        (foo1, _sum_expect, []),
        (foo2, _sum_expect, _body_expect),
        (foo3, "", _body_expect),
        (foo4, "", []),
    ],
)
def test_docstring_to_summary_rest(foo, sum_exp, body_exp):
    summary, body = docstring_to_summary_rest(foo.__doc__)
    assert summary == sum_exp
    assert body.split() == body_exp


def test_get_true_spans_absolute():
    got = get_true_spans(array([0, 1, 1, 0, 1]))
    assert_allclose(got, array([[1, 2], [4, 1]]))
    got = get_true_spans(array([0, 0]))
    assert not len(got)

    got = get_true_spans(array([0, 0, 1, 1]))
    assert_allclose(got, array([[2, 2]]))
    got = get_true_spans(array([1, 0, 0, 1, 1]))
    assert_allclose(got, array([[0, 1], [3, 2]]))
    got = get_true_spans(array([1, 0, 0, 1, 1, 1, 1, 0, 1, 1]))
    assert_allclose(got, array([[0, 1], [3, 4], [8, 2]]))
    got = get_true_spans(
        # abs  0  1  2  3  4  5  6  7  8  9
        array([1, 0, 0, 1, 1, 1, 1, 0, 1, 1]),
        absolute_pos=False,
        # u       0  1              2
    )
    assert_allclose(got, array([[0, 1], [2, 4], [3, 2]]))

    got = get_true_spans(array([0, 0, 0, 1, 1, 1, 0, 0]))
    assert_allclose(got, array([(3, 3)]))


def test_get_true_spans_not_absolute():
    got = get_true_spans(array([0, 1, 1, 0, 1]), absolute_pos=False)
    assert_allclose(got, array([[1, 2], [2, 1]]))

    got = get_true_spans(array([1, 0, 0, 1, 1, 1, 1, 0, 1, 1]), absolute_pos=False)
    assert_allclose(got, array([[0, 1], [2, 4], [3, 2]]))


def test_negate_condition():
    def greater_than_5(x):
        return x > 5

    negator = negate_condition(greater_than_5)

    result_true = negator(3)
    result_false = negator(8)

    assert result_true is True
    assert result_false is False
    assert greater_than_5(3) != negator(3)
