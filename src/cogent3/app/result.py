import json

from collections import OrderedDict
from collections.abc import MutableMapping
from functools import total_ordering
from pathlib import Path

import numpy

from cogent3.app.data_store import get_data_source
from cogent3.maths.stats import chisqprob
from cogent3.util.misc import extend_docstring_from, get_object_provenance
from cogent3.util.table import Table


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class generic_result(MutableMapping):
    """A dict style container for storing results."""

    _type = "generic_result"
    _item_types = ()

    def __init__(self, source):
        source = get_data_source(source)
        if not isinstance(source, (str, Path)):
            raise ValueError(f"Cannot infer source from type {type(source)}")

        self._store = {}
        self._construction_kwargs = dict(source=source)
        self.source = source

    def __setitem__(self, key, val):
        if isinstance(val, dict):
            type_name = val.get("type", None)
            type_name = type_name or ""
        else:
            type_name = val.__class__.__name__

        for item_type in self._item_types:
            if item_type in type_name:
                break
        else:
            if self._item_types:
                class_name = self.__class__.__name__
                msg = (
                    f"{type_name!r} not a supported value type for {class_name!r}, "
                    f"supported value types are {self._item_types}"
                )
                raise TypeError(msg)

        if not hasattr(val, "to_json"):
            json.dumps(val)

        self._store[key] = val

    def __getitem__(self, key):
        return self._store[key]

    def __delitem__(self, key):
        raise NotImplementedError

    def __len__(self):
        return len(self._store)

    def __iter__(self):
        return iter(self._store)

    def __repr__(self):
        name = self.__class__.__name__
        num = len(self)
        types = [f"{repr(k)}: {self[k].__class__.__name__}" for k in self]
        types = types[:3] + ["..."] if num > 5 else types
        types = ", ".join(types)
        return f"{num}x {name}({types})"

    def __str__(self):
        return repr(self)

    def keys(self):
        return list(self)

    def to_rich_dict(self):
        """returns the rich dict on values"""
        result = {
            "type": get_object_provenance(self),
            "result_construction": self._construction_kwargs,
            "version": __version__,
        }
        items = []
        for key, val in self.items():
            try:
                val = val.to_rich_dict()
            except AttributeError:
                pass
            items.append([key, val])
        result["items"] = items
        return result

    def to_json(self):
        data = self.to_rich_dict()
        return json.dumps(data)

    def deserialised_values(self):
        """deserialises any cogent3 members"""
        from cogent3.util.deserialise import deserialise_object

        for key, value in self.items():
            if isinstance(value, dict):
                type_ = value.get("type", "")
                if "cogent3" in type_:
                    object = deserialise_object(value)
                    self[key] = object
            elif hasattr(value, "deserialised_values"):
                value.deserialised_values()


@total_ordering
class model_result(generic_result):
    """Storage of model results."""

    _type = "model_result"
    _stat_attrs = ("lnL", "nfp", "DLC", "unique_Q")
    _item_types = ("AlignmentLikelihoodFunction",)

    def __init__(
        self,
        name=None,
        stat=sum,
        source=None,
        elapsed_time=None,
        num_evaluations=None,
        evaluation_limit=None,
        lnL=None,
        nfp=None,
        DLC=None,
        unique_Q=None,
    ):
        super(model_result, self).__init__(source)
        if type(stat) == str:
            stat = eval(stat)

        self._construction_kwargs.update(
            dict(
                name=name,
                stat=stat.__name__,
                elapsed_time=elapsed_time,
                num_evaluations=num_evaluations,
            )
        )
        self._store = {}
        self._name = name
        assert stat is sum or stat is max
        self._stat = stat
        self._elapsed_time = elapsed_time
        self._num_evaluations = num_evaluations
        self._evaluation_limit = evaluation_limit
        self._lnL = lnL
        self._nfp = nfp
        self._DLC = DLC
        self._unique_Q = unique_Q

    def _get_repr_data_(self):
        if len(self) == 0:
            return f"{self.__class__.__name__}(name={self.name}, source={self.source})"

        self.deserialised_values()  # making sure we're fully reloaded
        attrs = list(self._stat_attrs)
        header = ["key"] + attrs[:]
        rows = [[repr("")] + [getattr(self, attr) for attr in attrs]]
        if len(self) > 1:
            # we just add keys, lnL and nfp
            for key in self:
                row = [repr(key), self[key].lnL, self[key].nfp, "", ""]
                rows.append(row)
        else:
            rows[0][0] = repr(list(self)[0])

        return Table(header=header, data=rows, title=self.name)

    def _repr_html_(self):
        table = self._get_repr_data_()
        table.set_repr_policy(show_shape=False)
        return table._repr_html_()

    def __repr__(self):
        table = self._get_repr_data_()
        return repr(table)

    def __setitem__(self, key, lf):
        super(self.__class__, self).__setitem__(key, lf)
        self._init_stats()

    def _init_stats(self):
        """reset the values for stat attr to None, triggers recalc in properties"""
        for attr in self._stat_attrs:
            setattr(self, f"_{attr}", None)

    @property
    def num_evaluations(self):
        return self._num_evaluations

    @num_evaluations.setter
    def num_evaluations(self, value):
        value = int(value)
        self._num_evaluations = value
        self._construction_kwargs["num_evaluations"] = value

    @property
    def elapsed_time(self):
        return self._elapsed_time

    @elapsed_time.setter
    def elapsed_time(self, value):
        self._elapsed_time = value
        self._construction_kwargs["elapsed_time"] = value

    @property
    def name(self):
        return self._name

    def simulate_alignment(self):
        self.deserialised_values()
        if len(self) == 1:
            aln = self.lf.simulate_alignment()
            return aln
        # assume we have results from 3 codon positions
        sim = []
        seqnames = None
        for i in sorted(self):
            aln = self[i].simulate_alignment()
            sim.append(aln.to_dict())
            if seqnames is None:
                seqnames = list(sim[-1].keys())

        data = {}
        for n in seqnames:
            seq1, seq2, seq3 = sim[0][n], sim[1][n], sim[2][n]
            seq = "".join(("".join(t) for t in zip(seq1, seq2, seq3)))
            data[n] = seq

        return aln.__class__(data=data)

    def __lt__(self, other):
        self_lnL = self.lnL
        other_lnL = other.lnL
        return self_lnL < other_lnL

    @property
    def lf(self):
        self.deserialised_values()
        self._init_stats()
        if len(self) == 1:
            result = list(self.values())[0]
            result.name = self.name
        else:
            result = OrderedDict()
            for k in sorted(self):
                v = self[k]
                if type(k) == str and k.isdigit():
                    k = int(k)
                result[k] = v
                v.name = f"{self.name} pos-{k}"

        return result

    @property
    def lnL(self):
        if self._lnL is None:
            lnL = 0.0
            for v in self.values():
                l = v.get("lnL") if isinstance(v, dict) else v.lnL
                lnL = self._stat([l, lnL])

            self._lnL = lnL
        return self._lnL

    @property
    def nfp(self):
        if self._nfp is None:
            nfp = 0
            for v in self.values():
                n = v.get("nfp") if isinstance(v, dict) else v.nfp
                nfp = self._stat([n, nfp])

            self._nfp = nfp

        return self._nfp

    @property
    def DLC(self):
        if self._DLC is None:
            DLC = []
            for v in self.values():
                d = v.get("DLC") if isinstance(v, dict) else v.all_psubs_DLC()
                DLC.append(d != False)

            self._DLC = all(DLC)

        return self._DLC

    @property
    def unique_Q(self):
        if self._unique_Q is None:
            unique = []
            for v in self.values():
                if isinstance(v, dict):
                    u = v.get("unique_Q")
                else:
                    try:
                        u = v.all_rate_matrices_unique()
                    except (NotImplementedError, KeyError):
                        # KeyError happens on discrete time model
                        u = None  # non-primary root issue
                unique.append(u != False)

            self._unique_Q = all(unique)

        return self._unique_Q

    def total_length(self, length_as=None):
        """sum of all branch lengths on tree. If split codons, sums across trees

        Parameters
        ----------
        length_as : str or None
            replaces 'length' param with either 'ENS' or 'paralinear'.
            'ENS' is the expected number of substitution, (which will be
            different to standard length if the substitution model is
            non-stationary). 'paralinear' is the measure of Lake 1994.
        """
        if len(self) == 1:
            tree = self.lf.get_annotated_tree(length_as=length_as)
            return tree.total_length()

        total_length = 0
        for lf in self.lf.values():
            tree = lf.get_annotated_tree(length_as=length_as)
            total_length += tree.total_length()

        return total_length

    @property
    def tree(self):
        """an annotated tree with 'ENS' set as the branch length
        Note
        ----
        In the case of a discrete time process, length is 'paralinear'"""
        from cogent3.evolve.ns_substitution_model import (
            DiscreteSubstitutionModel,
        )

        try:
            model = self.lf.model
        except AttributeError:
            model = self.lf[1].model

        if isinstance(model, DiscreteSubstitutionModel):
            length_as = "paralinear"
        else:
            length_as = "ENS"

        if not hasattr(self, "_tree"):
            if len(self) == 1:
                tree = self.lf.get_annotated_tree(length_as=length_as)
            else:
                tree = OrderedDict()
                for k in sorted(self):
                    v = self[k]
                    if type(k) == str and k.isdigit():
                        k = int(k)
                    tree[k] = v.get_annotated_tree(length_as=length_as)

            self._tree = tree

        return self._tree

    @property
    def alignment(self):
        if len(self) == 1:
            result = self.lf.get_param_value("alignment")
        else:
            result = OrderedDict()
            for k in sorted(self):
                v = self[k]
                if type(k) == str and k.isdigit():
                    k = int(k)
                result[k] = v.get_param_value("alignment")
        return result


class model_collection_result(generic_result):
    """Storage of a collection of model_result."""

    _type = "model_collection_result"
    _item_types = ("model_result",)

    def __init__(self, name=None, source=None):
        """
        name : str
            name of this hypothesis
        source : str
            string describing source of the data, e.g. a path
        """
        super(model_collection_result, self).__init__(source)
        self._construction_kwargs.update({"name": name})
        self._name = name

    def _get_repr_data_(self):
        rows = []
        attrs = ["lnL", "nfp", "DLC", "unique_Q"]
        for key, member in self.items():
            member.deserialised_values()  # making sure we're fully reloaded
            row = [repr(key)] + [getattr(member, a) for a in attrs]
            rows.append(row)

        table = Table(header=["key"] + attrs, data=rows, title=self.name)
        table = table.sorted(columns="nfp")
        return table

    def _repr_html_(self):
        table = self._get_repr_data_()
        table.set_repr_policy(show_shape=False)
        return table._repr_html_()

    def __repr__(self):
        if len(self) == 0:
            return f"{self.__class__.__name__}(name={self.name}, source={self.source})"

        table = self._get_repr_data_()
        return str(table._get_repr_())

    @property
    def name(self):
        return self._name

    def select_models(self, stat="aicc", threshold=0.05):
        """returns models satisfying stat threshold.
        Parameters
        ----------
        stat : str
            one of "aicc", "aic" which correspond to
            AIC with correction or AIC.
        threshold : float
            models with exp((minimum stat - model stat) / 2) > threshold are
            considered indistinguishable from the model with minimum stat. Such
            models will be included in the returned result.

        Returns
        -------
        list of models satisfying threshold condition
        """
        self.deserialised_values()
        assert stat in ("aicc", "aic")
        second_order = stat == "aicc"
        results = []
        for m in self.values():
            if isinstance(m.lf, dict):
                # multiple lf's, e.g. split codon position analyses have 3
                val = sum(lf.get_aic(second_order=second_order) for lf in m.lf.values())
            else:
                val = m.lf.get_aic(second_order=second_order)

            results.append((val, m))
        results.sort()
        min_model = results.pop(0)
        min_stat = min_model[0]
        selected = [min_model[1]]
        for v, m in results:
            rel_lik = numpy.exp((min_stat - v) / 2)
            if rel_lik > threshold:
                selected.append(m)

        return selected

    def get_best_model(self, stat="aicc", threshold=0.05):
        """returns model with smallest value of stat
        Parameters
        ----------
        stat : str
            one of "aicc", "aic" which correspond to AIC with correction or AIC
        threshold : float
            models with exp((minimum stat - model stat) / 2) > threshold are
            considered indistinguishable from the model with minimum stat.

        Returns
        -------
        A single model. If multiple models satisfy threshold, the simplest model
        (with the smallest number of free parameters) is returned.
        """
        selected = self.select_models(stat=stat, threshold=threshold)
        if len(selected) != 1:
            selected = list(sorted(self.values(), key=lambda x: x.nfp))
            selected = selected[:1]

        return selected[0]

    def get_hypothesis_result(self, name_null, name_alt):
        """returns a hypothesis result with two models

        Parameters
        ----------
        name_null : str
            name of the null model
        name_alt : str
            name of the alternate model
        """
        result = hypothesis_result(name_of_null=name_null, source=self.source)
        result[name_null] = self[name_null]
        result[name_alt] = self[name_alt]
        return result


class hypothesis_result(model_collection_result):
    """Storage of a collection of model_result instances that are hierarchically
    related."""

    _type = "hypothesis_result"
    _item_types = ("model_result",)

    @extend_docstring_from(model_collection_result.__init__, pre=True)
    def __init__(self, name_of_null, name=None, source=None):
        """
        name_of_null
            key for the null hypothesis
        """
        super(hypothesis_result, self).__init__(name=name, source=source)
        self._construction_kwargs.update(dict(name_of_null=name_of_null))

        self._name_of_null = name_of_null

    def _get_repr_data_(self):
        rows = []
        attrs = ["lnL", "nfp", "DLC", "unique_Q"]
        for key, member in self.items():
            member.deserialised_values()  # making sure we're fully reloaded
            if key == self._name_of_null:
                status_name = ["null", repr(key)]
            else:
                status_name = ["alt", repr(key)]
            row = status_name + [getattr(member, a) for a in attrs]
            rows.append(row)

        table = Table(header=["hypothesis", "key"] + attrs, data=rows, title=self.name)
        table = table.sorted(columns="nfp")
        table.set_repr_policy(show_shape=False)
        stats = [[self.LR, self.df, self.pvalue]]
        col_templates = (
            None
            if self.pvalue is None
            else {
                "pvalue": "%.4f" if self.pvalue > 1e-3 else "%.2e",
            }
        )
        stats = Table(
            header=["LR", "df", "pvalue"],
            data=stats,
            title="Statistics",
            column_templates=col_templates,
        )
        stats.set_repr_policy(show_shape=False)
        return stats, table

    def _repr_html_(self):
        stats, table = self._get_repr_data_()
        result = [t._repr_html_() for t in (stats, table)]
        return "\n".join(result)

    def __repr__(self):
        if len(self) == 0:
            return f"{self.__class__.__name__}(name={self.name}, source={self.source})"

        stats, table = self._get_repr_data_()
        result = []
        for t in (stats, table):
            r, _, _ = t._get_repr_()
            result.append(str(r))

        return "\n".join(result)

    @property
    def null(self):
        return self[self._name_of_null]

    @property
    def alt(self):
        alts = [self[k] for k in self if k != self._name_of_null]
        return max(alts)

    @property
    def LR(self):
        """returns 2 * (alt.lnL - null.lnL)"""
        LR = self.alt.lnL - self.null.lnL
        LR *= 2
        return LR

    @property
    def df(self):
        """returns the degrees-of-freedom (alt.nfp - null.nfp)"""
        return self.alt.nfp - self.null.nfp

    @property
    def pvalue(self):
        """returns p-value from chisqprob(LR, df)

        None if LR < 0"""
        if self.LR == 0:
            pvalue = 1
        elif self.LR > 0:
            pvalue = chisqprob(self.LR, self.df)
        else:
            pvalue = None
        return pvalue


class bootstrap_result(generic_result):
    _type = "bootstrap_result"
    _item_types = ("hypothesis_result", "model_collection_result")

    def __init__(self, source=None):
        super(bootstrap_result, self).__init__(source)

    @property
    def observed(self):
        """the results for the observed data"""
        return self["observed"]

    @observed.setter
    def observed(self, data):
        self.update(dict(observed=data))

    def add_to_null(self, data):
        """add results for a synthetic data set"""
        size = len(self)
        self[size + 1] = data

    @property
    def null_dist(self):
        """returns the LR values corresponding to the synthetic data"""
        return [self[k].LR for k in self if k != "observed"]


class tabular_result(generic_result):
    """stores one or multiple cogent3 Tables, DictArray"""

    _type = "tabular_result"
    _stat_attrs = ("header", "rows")
    _item_types = (
        "Table",
        "DictArray",
        "MotifCounts",
        "MotifFreqs",
        "PSSM",
        "DistanceMatrix",
    )

    def __init__(self, source=None):
        super(tabular_result, self).__init__(source)
