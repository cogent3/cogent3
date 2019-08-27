import json

from collections import OrderedDict
from collections.abc import MutableMapping
from functools import total_ordering

import numpy

from cogent3.maths.stats import chisqprob
from cogent3.util.misc import get_object_provenance


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.8.23a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class generic_result(MutableMapping):
    """a dict style container for storing results. All keys are
    converted to strings to ensure the object can be json serialised"""

    _type = "generic_result"

    def __init__(self, source):
        self._store = dict()
        self._construction_kwargs = dict(source=source)
        self.source = source

    def __setitem__(self, key, val):
        key = str(key)
        self._store[key] = val

    def __getitem__(self, key):
        key = str(key)
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
        types = [f"{repr(k)}: {self[k].__class__.__name__}" for k in self.keys()[:4]]
        types = ", ".join(types)
        result = f"{len(self)}x {name}({types})"
        return result

    def __str__(self):
        return repr(self)

    def keys(self):
        return list(self)

    def to_rich_dict(self):
        """returns the rich dict on values"""
        result = {
            "type": get_object_provenance(self),
            "result_construction": self._construction_kwargs,
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


@total_ordering
class model_result(generic_result):
    _type = "model_result"
    _stat_attrs = ("lnL", "nfp", "DLC", "unique_Q")

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

        self._construction_kwargs = dict(
            name=name,
            stat=stat.__name__,
            source=source,
            elapsed_time=elapsed_time,
            num_evaluations=num_evaluations,
        )
        self._store = dict()
        self._name = name
        assert stat is sum or stat is max
        self._stat = stat
        self._elapsed_time = elapsed_time
        self._num_evaluations = num_evaluations
        self._evaluation_limit = evaluation_limit
        self._lnL = None
        self._nfp = None
        self._DLC = None
        self._unique_Q = None

    def _get_repr_data_(self):
        from cogent3.util.table import Table

        self.lf  # making sure we're fully reloaded
        attrs = ["lnL", "nfp", "DLC", "unique_Q"]
        header = ["key"] + attrs[:]
        rows = [[""] + [getattr(self, attr) for attr in attrs]]
        if len(self) > 1:
            # we just add keys, lnL and nfp
            padd = ["", ""]
            attrs = ["lnL", "nfp"]
            for key in self:
                row = [repr(key), self[key].lnL, self[key].nfp, "", ""]
                rows.append(row)

        table = Table(header=header, rows=rows, title=self.name)
        return table

    def _repr_html_(self):
        table = self._get_repr_data_()
        return table._repr_html_(include_shape=False)

    def __repr__(self):
        table = self._get_repr_data_()
        return repr(table)

    def __setitem__(self, key, lf):
        super(self.__class__, self).__setitem__(key, lf)
        if type(lf) != dict:
            lf.set_name(key)
            lnL = lf.lnL
            nfp = lf.nfp
            DLC = lf.all_psubs_DLC()
            try:
                unique_Q = lf.all_rate_matrices_unique()
            except (NotImplementedError, KeyError):
                # KeyError happens on discrete time model
                unique_Q = None  # non-primary root issue
        else:
            lnL = lf.get("lnL")
            nfp = lf.get("nfp")
            DLC = lf.get("DLC")
            unique_Q = lf.get("unique_Q")

        if self.lnL is not None:
            self.DLC = all([DLC, self.DLC])
            self.unique_Q = all([unique_Q, self.unique_Q])
            self.lnL = self._stat([lnL, self.lnL])
            self.nfp = self._stat([nfp, self.nfp])
        else:
            self.lnL = lnL
            self.nfp = nfp
            self.DLC = DLC
            self.unique_Q = unique_Q

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

        simaln = aln.__class__(data=data)

        return simaln

    def __lt__(self, other):
        self_lnL = self.lnL
        other_lnL = other.lnL
        return self_lnL < other_lnL

    @property
    def lf(self):
        result = list(self.values())
        if type(result[0]) == dict:
            from cogent3.util import deserialise

            # we reset the stat attributes to None
            for attr in self._stat_attrs:
                setattr(self, attr, None)

            for k, v in self.items():
                v = deserialise.deserialise_likelihood_function(v)
                self[k] = v

        if len(self) == 1:
            result = list(self.values())[0]
        else:
            result = OrderedDict()
            for k in sorted(self):
                v = self[k]
                if k.isdigit():
                    k = int(k)
                result[k] = v

        return result

    @property
    def lnL(self):
        return self._lnL

    @lnL.setter
    def lnL(self, value):
        self._lnL = value

    @property
    def nfp(self):
        return self._nfp

    @nfp.setter
    def nfp(self, value):
        self._nfp = value

    @property
    def DLC(self):
        return self._DLC

    @DLC.setter
    def DLC(self, value):
        self._DLC = value

    @property
    def unique_Q(self):
        return self._unique_Q

    @unique_Q.setter
    def unique_Q(self, value):
        self._unique_Q = value

    def total_length(self, length_as=None):
        """sum of all branch lengths on tree

        Parameters
        ----------
        length_as : str or None
            replaces 'length' param with either 'ENS' or 'paralinear'.
            'ENS' is the expected number of substitution, (which will be
            different to standard length if the substitution model is
            non-stationary). 'paralinear' is the measure of Lake 1994.
        """
        tree = self.lf.get_annotated_tree(length_as=length_as)
        return tree.total_length()


class hypothesis_result(generic_result):
    _type = "hypothesis_result"

    def __init__(self, name_of_null, source=None):
        """
        alt
            either a likelihood function instance
        """
        super(hypothesis_result, self).__init__(source)
        self._construction_kwargs = dict(name_of_null=name_of_null, source=source)

        self._name_of_null = name_of_null

    def _get_repr_data_(self):
        from cogent3.util.table import Table

        rows = []
        attrs = ["lnL", "nfp", "DLC", "unique_Q"]
        for key, member in self.items():
            member.lf  # making sure we're fully reloaded
            if key == self._name_of_null:
                status_name = ["null", repr(key)]
            else:
                status_name = ["alt", (repr(key))]
            row = status_name + [getattr(member, a) for a in attrs]
            rows.append(row)

        table = Table(header=["hypothesis", "key"] + attrs, rows=rows)
        table = table.sorted(columns="nfp")
        stats = [[self.LR, self.df, self.pvalue]]
        stats = Table(header=["LR", "df", "pvalue"], rows=stats, title="Statistics")
        return stats, table

    def _repr_html_(self):
        stats, table = self._get_repr_data_()
        result = [t._repr_html_(include_shape=False) for t in (stats, table)]
        return "\n".join(result)

    def __repr__(self):
        stats, table = self._get_repr_data_()
        result = []
        for t in (stats, table):
            r, _ = t._get_repr_()
            result.append(str(r))

        return "\n".join(result)

    @property
    def null(self):
        return self[self._name_of_null]

    @property
    def alt(self):
        alts = [self[k] for k in self if k != self._name_of_null]
        alt = max(alts)
        return alt

    @property
    def LR(self):
        """returns 2 * (alt.lnL - null.lnL)"""
        LR = self.alt.lnL - self.null.lnL
        LR *= 2
        return LR

    @property
    def df(self):
        """returns the degrees-of-freedom (alt.nfp - null.nfp)"""
        df = self.alt.nfp - self.null.nfp
        return df

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
        assert stat in ("aicc", "aic")
        second_order = stat == "aicc"
        results = []
        for k, m in self.items():
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


class bootstrap_result(generic_result):
    _type = "bootstrap_result"

    def __init__(self, source=None):
        super(bootstrap_result, self).__init__(source)
        self._construction_kwargs = dict(source=source)

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
        self.update({size + 1: data.to_rich_dict()})

    @property
    def null_dist(self):
        """returns the LR values corresponding to the synthetic data"""
        result = [self[k].LR for k in self if k != "observed"]
        return result


class tabular_result(generic_result):
    """stores one or multiple tabular data sets, keyed by a title"""

    _type = "tabular_result"
    _stat_attrs = ("header", "rows")

    def __init__(self, source=None):
        super(tabular_result, self).__init__(source)
