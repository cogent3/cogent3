import json
from collections import OrderedDict
from collections.abc import MutableMapping
from functools import total_ordering

from cogent3.util.misc import get_object_provenance
from cogent3.maths.stats import chisqprob


class generic_result(MutableMapping):
    type_ = 'generic_result'

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

    def to_rich_dict(self):
        """returns the rich dict on values"""
        result = {'type': get_object_provenance(self),
                  'result_construction': self._construction_kwargs}
        for key, val in self.items():
            try:
                val = val.to_rich_dict()
            except AttributeError:
                pass
            result[key] = val
        return result

    def to_json(self):
        data = self.to_rich_dict()
        return json.dumps(data)


@total_ordering
class model_result(generic_result):
    type_ = 'model_result'
    _stat_attrs = ('lnL', 'nfp', 'DLC', 'unique_Q')

    def __init__(self, name=None, stat=sum, source=None, elapsed_time=None,
                 num_evaluations=None, evaluation_limit=None,
                 lnL=None, nfp=None, DLC=None, unique_Q=None):
        super(model_result, self).__init__(source)
        if type(stat) == str:
            stat = eval(stat)

        self._construction_kwargs = dict(name=name, stat=stat.__name__,
                                         source=source,
                                         elapsed_time=elapsed_time,
                                         num_evaluations=num_evaluations)
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

    def __setitem__(self, key, lf):
        super(self.__class__, self).__setitem__(key, lf)
        if type(lf) != dict:
            lnL = lf.lnL
            nfp = lf.nfp
            DLC = lf.all_psubs_DLC()
            try:
                unique_Q = lf.all_rate_matrices_unique()
            except NotImplementedError:
                unique_Q = None  # non-primary root issue
        else:
            lnL = lf.get('lnL')
            nfp = lf.get('nfp')
            DLC = lf.get('DLC')
            unique_Q = lf.get('unique_Q')

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
        self._construction_kwargs['num_evaluations'] = value

    @property
    def elapsed_time(self):
        return self._elapsed_time

    @elapsed_time.setter
    def elapsed_time(self, value):
        self._elapsed_time = value
        self._construction_kwargs['elapsed_time'] = value

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
            sim.append(aln.todict())
            if seqnames is None:
                seqnames = list(sim[-1].keys())

        data = {}
        for n in seqnames:
            seq1, seq2, seq3 = sim[0][n], sim[1][n], sim[2][n]
            seq = ''.join((''.join(t) for t in zip(seq1, seq2, seq3)))
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


class hypothesis_result(generic_result):
    type_ = 'hypothesis_result'

    def __init__(self, name_of_null, source=None):
        """
        alt
            either a likelihood function instance
        """
        super(hypothesis_result, self).__init__(source)
        self._construction_kwargs = dict(name_of_null=name_of_null,
                                         source=source)

        self._name_of_null = name_of_null

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


class bootstrap_result(generic_result):
    type_ = 'bootstrap_result'

    def __init__(self, source=None):
        super(bootstrap_result, self).__init__(source)
        self._construction_kwargs = dict(source=source)

    @property
    def observed(self):
        """the results for the observed data"""
        return self['observed']

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
        result = [self[k].LR for k in self if k != 'observed']
        return result
