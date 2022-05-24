#!/usr/bin/env python
"""
Matching motifs: MotifChange("a", "g")
Boolean logic: Any, All, Not (or &, |, ~)

also:
  anypredicate.aliased('shortname')
  UserPredicate(f)
"""

import re
import warnings

import numpy


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


class _CallablePredicate(object):
    # A predicate in the context of a particular model

    def __init__(self, pred, model):
        self.model = model
        self.alphabet = model.get_alphabet()
        self.name = repr(pred)
        self.f = pred.interpret(model)
        self.__doc__ = pred.__doc__

    def __call__(self, x, y):
        return self.f(x, y)

    def __repr__(self):
        return self.name

    def ascii_art(self):
        l = len(self.alphabet.get_gap_motif())
        rows = []
        for i in range(l):
            row = [a[i] for a in list(self.alphabet)]
            rows.append(" " * (l + 1) + " ".join(row))
        for y in self.alphabet:
            row = []
            for x in self.alphabet:
                if not self.model.isinstantanious(x, y):
                    c = " "
                elif self(x, y):
                    c = "*"
                else:
                    c = "-"
                row.append(c)
            rows.append(" ".join([y] + row))
        return "\n".join(rows)


class predicate(object):
    def __and__(self, other):
        return All(self, other)

    def __or__(self, other):
        return Any(self, other)

    def __invert__(self):
        return Not(self)

    def __bool__(self):
        warnings.warn(
            "alphabet predicate used as truth value. Use only binary operators: &, | and ~",
            stacklevel=2,
        )
        return True

    def __eq__(self, other):
        warnings.warn("Warning: alphabet pair predicate used as value. Use parentheses")
        return self is other

    def aliased(self, new_name):
        n = PredicateAlias(new_name, self)
        n.__doc__ = self.__doc__ or repr(self)
        return n

    def make_model_predicate(self, model):
        return _CallablePredicate(self, model)


class PredicateAlias(predicate):
    def __init__(self, name, subpredicate):
        self.name = name
        self.subpredicate = subpredicate

    def __repr__(self):
        return self.name

    def interpret(self, model):
        return self.subpredicate.interpret(model)


class _UnaryPredicate(predicate):
    def __init__(self, subpredicate):
        assert isinstance(subpredicate, predicate), subpredicate
        self.subpredicate = subpredicate
        self.__doc__ = repr(self)

    def __repr__(self):
        if hasattr(self, "_op_repr"):
            return f"{self._op_repr}({self.subpredicate})"
        else:
            return f"{self.__class__.__name__}({self.subpredicate})"


class _GenericPredicate(predicate):
    def __init__(self, *subpredicates):
        for p in subpredicates:
            assert isinstance(p, predicate), p
        self.subpredicates = subpredicates
        self.__doc__ = repr(self)

    def __repr__(self):
        if hasattr(self, "_op_repr"):
            return "(%s)" % (" %s " % self._op_repr).join(
                [repr(p) for p in self.subpredicates]
            )
        else:
            return "%s(%s)" % (
                self.__class__.__name__,
                ",".join(["(%s)" % repr(p) for p in self.subpredicates]),
            )


# Boolean logic on motif pair predicates


class Not(_UnaryPredicate):
    _op_repr = "~"

    def interpret(self, model):
        subpred = self.subpredicate.interpret(model)

        def call(*args):
            return not subpred(*args)

        call.__doc__ = repr(self)
        return call


class All(_GenericPredicate):
    _op_repr = "&"

    def interpret(self, model):
        subpreds = [p.interpret(model) for p in self.subpredicates]

        def call(*args):
            for subpredicate in subpreds:
                if not subpredicate(*args):
                    return False
            return True

        call.__doc__ = repr(self)
        return call


class Any(_GenericPredicate):
    _op_repr = "|"

    def interpret(self, model):
        subpreds = [p.interpret(model) for p in self.subpredicates]

        def call(*args):
            for subpredicate in subpreds:
                if subpredicate(*args):
                    return True
            return False

        call.__doc__ = repr(self)
        return call


class ModelSays(predicate):
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return self.name

    def interpret(self, model):
        return model.get_predefined_predicate(self.name)


class DirectedMotifChange(predicate):
    def __init__(self, from_motif, to_motif, diff_at=None):

        self.from_motif = from_motif
        self.motiflen = len(from_motif)
        self.to_motif = to_motif
        self.diff_at = diff_at

    def __repr__(self):
        if self.diff_at is not None:
            diff = "[%d]" % self.diff_at
        else:
            diff = ""
        return f"{self.from_motif}>{self.to_motif}{diff}"

    def test_motif(self, motifs, query):
        """positions where motif pattern is found in query"""
        positions = set()
        for offset in range(len(query) - self.motiflen + 1):
            for (q, ms) in zip(query[offset : offset + self.motiflen], motifs):
                if q not in ms:
                    break
            else:
                positions.add(offset)
        return positions

    def test_motifs(self, from_motifs, to_motifs, x, y):
        """ "positions where both motifs patterns are found"""
        pre = self.test_motif(from_motifs, x)
        post = self.test_motif(to_motifs, y)
        return pre & post

    def interpret(self, model):
        """Make a callable function which implements this predicate
        specificly for 'alphabet'"""

        # may be looking for a 2nt pattern in a 3nt alphabet, but not
        # 3nt pattern in dinucleotide alphabet.
        alphabet = model.get_alphabet()
        if alphabet.get_motif_len() < self.motiflen:
            raise ValueError(
                "alphabet motifs (%s) too short for %s (%s)"
                % (alphabet.get_motif_len(), repr(self), self.motiflen)
            )

        resolve = model.moltype.ambiguities.__getitem__

        from_motifs = [resolve(m) for m in self.from_motif]
        to_motifs = [resolve(m) for m in self.to_motif]

        def call(x, y):
            diffs = [X != Y for (X, Y) in zip(x, y)]
            matches = []
            for posn in self.test_motifs(from_motifs, to_motifs, x, y):
                diff = list(numpy.nonzero(diffs[posn : posn + self.motiflen])[0])
                if diff and self.diff_at is None or diff == [self.diff_at]:
                    matches.append(posn)
            return len(matches) == 1

        call.__doc__ = repr(self)
        return call


class UndirectedMotifChange(DirectedMotifChange):
    def __repr__(self):
        if self.diff_at is not None:
            diff = "[%d]" % self.diff_at
        else:
            diff = ""
        return f"{self.from_motif}/{self.to_motif}{diff}"

    def test_motifs(self, from_motifs, to_motifs, x, y):
        preF = self.test_motif(from_motifs, x)
        postF = self.test_motif(to_motifs, y)
        preR = self.test_motif(from_motifs, y)
        postR = self.test_motif(to_motifs, x)
        return (preF & postF) | (preR & postR)


def MotifChange(x, y=None, forward_only=False, diff_at=None):
    if y is None:
        y = ""
        for i in range(len(x)):
            if i == diff_at or diff_at is None:
                y += "?"
            else:
                y += x[i]
    if forward_only:
        return DirectedMotifChange(x, y, diff_at=diff_at)
    else:
        return UndirectedMotifChange(x, y, diff_at=diff_at)


class UserPredicate(predicate):
    def __init__(self, f):
        self.f = f

    def __repr__(self):
        return "UserPredicate(%s)" % (getattr(self.f, "__name__", None) or repr(self.f))

    def interpret(self, model):
        return self.f


silent = ModelSays("silent")
replacement = ModelSays("replacement")
omega = ModelSays("omega")


def parse(rule):
    if "|" in rule:
        rules = re.sub("[()]", "", rule)
        preds = [parse(r.strip()) for r in rules.split("|")]
        pred = preds.pop(0)
        for p in preds:
            pred = pred | p
        return pred

    if ":" in rule:
        (label, rule) = rule.split(":")
    else:
        label = None
    if "@" in rule:
        (rule, diff_at) = rule.split("@")
        diff_at = int(diff_at)
    else:
        diff_at = None

    if ">" in rule or "/" in rule:
        forward_only = ">" in rule
        rule = rule.replace(">", "/")
        (x, y) = rule.split("/")
        if not y:
            y = None

        pred = MotifChange(x, y, forward_only, diff_at)
    else:
        pred = ModelSays(rule)

    if label:
        pred = pred.aliased(label)
    return pred
