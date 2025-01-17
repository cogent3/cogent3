import warnings
from typing import Union

import numpy

from cogent3.core.alphabet import Alphabet
from cogent3.evolve.likelihood_tree import make_likelihood_tree_leaf
from cogent3.recalculation.definition import CalcDefn, PartitionDefn


def make_model(mprob_model, tuple_alphabet, mask):
    if mprob_model == "monomers":
        return PosnSpecificMonomerProbModel(tuple_alphabet, mask)
    if mprob_model == "monomer":
        return MonomerProbModel(tuple_alphabet, mask)
    if mprob_model == "conditional":
        return ConditionalMotifProbModel(tuple_alphabet, mask)
    if mprob_model in ["word", "tuple", None]:
        return SimpleMotifProbModel(tuple_alphabet)
    msg = f"Unknown mprob model '{mprob_model!s}'"
    raise ValueError(msg)


class MotifProbModel:
    def __init__(self, *whatever, **kw) -> None:
        raise NotImplementedError

    def calc_word_probs(self, *monomer_probs):
        assert len(monomer_probs) == 1
        return monomer_probs[0]

    def calc_word_weight_matrix(self, *monomer_probs):
        assert len(monomer_probs) == 1
        return monomer_probs[0]

    def make_motif_probs_defn(self):
        """Makes the first part of a parameter controller definition for this
        model, the calculation of motif probabilities"""
        return PartitionDefn(
            name="mprobs",
            default=None,
            dimensions=("locus", "edge"),
            dimension=("motif", tuple(self.get_input_alphabet())),
        )

    def set_param_controller_motif_probs(self, pc, motif_probs, **kw) -> None:
        pc.set_param_rule("mprobs", value=motif_probs, **kw)

    def count_motifs(self, alignment, include_ambiguity=False, recode_gaps=True):
        result = None
        alpha = self.get_counted_alphabet()
        for seq_name in alignment.names:
            sequence = alignment.get_gapped_seq(seq_name, recode_gaps)
            leaf = make_likelihood_tree_leaf(sequence, alpha, seq_name)
            count = leaf.get_motif_counts(include_ambiguity=include_ambiguity)
            if result is None:
                result = count.copy()
            else:
                result += count
        return result

    def adapt_motif_probs(self, motif_probs, **kwargs):
        return adapt_motif_probs(self.get_input_alphabet(), motif_probs)

    def make_equal_motif_probs(self):
        alphabet = self.get_input_alphabet()
        p = 1.0 / len(alphabet)
        return {m: p for m in alphabet}

    def make_sample_motif_probs(self):
        import random

        motif_probs = numpy.array(
            [random.uniform(0.2, 1.0) for m in self.get_counted_alphabet()],
        )
        motif_probs /= sum(motif_probs)
        return motif_probs


class SimpleMotifProbModel(MotifProbModel):
    def __init__(self, alphabet) -> None:
        self.alphabet = alphabet

    def get_input_alphabet(self):
        return self.alphabet

    def get_counted_alphabet(self):
        return self.alphabet

    def make_motif_word_prob_defns(self):
        monomer_probs = self.make_motif_probs_defn()
        return (monomer_probs, monomer_probs, monomer_probs)


class ComplexMotifProbModel(MotifProbModel):
    def __init__(self, tuple_alphabet, mask) -> None:
        """Arguments:
        - tuple_alphabet: series of multi-letter motifs
        - monomers: the monomers from which the motifs are made
        - mask: instantaneous change matrix"""
        self.mask = mask
        self.tuple_alphabet = tuple_alphabet
        self.monomer_alphabet = monomers = tuple_alphabet.moltype.alphabet
        self.word_length = length = tuple_alphabet.motif_len
        size = len(tuple_alphabet)

        # m2w[AC, 1] = C
        # w2m[0, AC, A] = True
        # w2c[ATC, AT*] = 1
        self.m2w = m2w = numpy.zeros([size, length], int)
        self.w2m = w2m = numpy.zeros([length, size, len(monomers)], int)
        contexts = monomers.get_word_alphabet(length - 1)
        self.w2c = w2c = numpy.zeros([size, length * len(contexts)], int)
        for i, word in enumerate(tuple_alphabet):
            for j in range(length):
                monomer = monomers.index(word[j])
                context = contexts.index(word[:j] + word[j + 1 :])
                m2w[i, j] = monomer
                w2m[j, i, monomer] = 1
                w2c[i, context * length + j] = 1

        self.mutated_posn = numpy.zeros(mask.shape, int)
        self.mutant_motif = numpy.zeros(mask.shape, int)
        self.context_indices = numpy.zeros(mask.shape, int)

        for i, _old_word, j, new_word, diff in self._mutations():
            self.mutated_posn[i, j] = diff
            mutant_motif = new_word[diff]
            context = new_word[:diff] + new_word[diff + 1 :]
            self.mutant_motif[i, j] = monomers.index(mutant_motif)
            c = contexts.index(context)
            self.context_indices[i, j] = c * length + diff

    def _mutations(self):
        def diff_pos(x, y):
            return [i for i in range(len(x)) if x[i] != y[i]]

        num_states = len(self.tuple_alphabet)
        for i in range(num_states):
            old_word = self.tuple_alphabet[i]
            for j in range(num_states):
                new_word = self.tuple_alphabet[j]
                if self.mask[i, j]:
                    assert self.mask[i, j] == 1.0
                    diffs = diff_pos(old_word, new_word)
                    assert len(diffs) == 1, (old_word, new_word)
                    diff = diffs[0]
                    yield i, old_word, j, new_word, diff


class MonomerProbModel(ComplexMotifProbModel):
    def get_input_alphabet(self):
        return self.monomer_alphabet

    def get_counted_alphabet(self):
        return self.monomer_alphabet

    def calc_monomer_probs(self, word_probs):
        monomer_probs = numpy.dot(word_probs, self.w2m.sum(axis=0))
        monomer_probs /= monomer_probs.sum()
        return monomer_probs

    def calc_word_probs(self, monomer_probs):
        result = numpy.prod(monomer_probs.take(self.m2w), axis=-1)
        # maybe simpler but slower, works ok:
        # result = numpy.prod(monomer_probs ** (w2m, axis=-1))
        result /= result.sum()
        return result

    def calc_word_weight_matrix(self, monomer_probs):
        return monomer_probs.take(self.mutant_motif) * self.mask

    def make_motif_word_prob_defns(self):
        monomer_probs = self.make_motif_probs_defn()
        word_probs = CalcDefn(self.calc_word_probs, name="wprobs")(monomer_probs)
        mprobs_matrix = CalcDefn(self.calc_word_weight_matrix, name="mprobs_matrix")(
            monomer_probs,
        )
        return (monomer_probs, word_probs, mprobs_matrix)

    def adapt_motif_probs(self, motif_probs, warn=False):
        try:
            motif_probs = adapt_motif_probs(self.monomer_alphabet, motif_probs)
        except ValueError:
            motif_probs = adapt_motif_probs(self.tuple_alphabet, motif_probs)
            if warn:
                warnings.warn("Motif probs over specified", stacklevel=5)
            motif_probs = self.calc_monomer_probs(motif_probs)
        return motif_probs


class PosnSpecificMonomerProbModel(MonomerProbModel):
    def get_counted_alphabet(self):
        return self.tuple_alphabet

    def calc_posn_specific_monomer_probs(self, word_probs):
        monomer_probs = numpy.dot(word_probs, self.w2m)
        monomer_probs /= monomer_probs.sum(axis=1)[..., numpy.newaxis]
        return list(monomer_probs)

    def calc_word_probs(self, monomer_probs):
        positions = list(range(self.word_length))
        assert len(monomer_probs) == self.m2w.shape[1], (
            len(monomer_probs),
            type(monomer_probs),
            self.m2w.shape,
        )
        result = numpy.prod(
            [monomer_probs[i].take(self.m2w[:, i]) for i in positions],
            axis=0,
        )
        result /= result.sum()
        return result

    def calc_word_weight_matrix(self, monomer_probs):
        monomer_probs = numpy.array(monomer_probs)  # so [posn, motif]
        size = monomer_probs.shape[-1]
        # should be constant
        extended_indices = self.mutated_posn * size + self.mutant_motif
        return monomer_probs.take(extended_indices) * self.mask

    def make_motif_word_prob_defns(self):
        monomer_probs = PartitionDefn(
            name="psmprobs",
            default=None,
            dimensions=("locus", "position", "edge"),
            dimension=("motif", tuple(self.get_input_alphabet())),
        )
        monomer_probs3 = monomer_probs.across_dimension(
            "position",
            [str(i) for i in range(self.word_length)],
        )
        monomer_probs3 = CalcDefn(lambda *x: numpy.array(x), name="mprobs")(
            *monomer_probs3,
        )
        word_probs = CalcDefn(self.calc_word_probs, name="wprobs")(monomer_probs3)
        mprobs_matrix = CalcDefn(self.calc_word_weight_matrix, name="mprobs_matrix")(
            monomer_probs3,
        )
        return (monomer_probs, word_probs, mprobs_matrix)

    def set_param_controller_motif_probs(self, pc, motif_probs, **kw) -> None:
        assert len(motif_probs) == self.word_length
        for i, m in enumerate(motif_probs):
            pc.set_param_rule("psmprobs", value=m, position=str(i), **kw)

    def adapt_motif_probs(self, motif_probs, **kwargs):
        try:
            motif_probs = adapt_motif_probs(self.monomer_alphabet, motif_probs)
        except ValueError:
            motif_probs = adapt_motif_probs(self.tuple_alphabet, motif_probs)
            motif_probs = self.calc_posn_specific_monomer_probs(motif_probs)
        else:
            motif_probs = [motif_probs] * self.word_length
        return motif_probs


class ConditionalMotifProbModel(ComplexMotifProbModel):
    def get_input_alphabet(self):
        return self.tuple_alphabet

    def get_counted_alphabet(self):
        return self.tuple_alphabet

    def calc_word_weight_matrix(self, motif_probs):
        context_probs = numpy.dot(motif_probs, self.w2c)
        context_probs[context_probs == 0.0] = numpy.inf
        return motif_probs / context_probs.take(self.context_indices)

    def make_motif_word_prob_defns(self):
        mprobs = self.make_motif_probs_defn()
        mprobs_matrix = CalcDefn(self.calc_word_weight_matrix, name="mprobs_matrix")(
            mprobs,
        )
        return (mprobs, mprobs, mprobs_matrix)


ProbsTypes = Union[dict, numpy.ndarray, list, tuple]
AlphaTypes = Union[Alphabet, tuple, list]


def adapt_motif_probs(alphabet: AlphaTypes, motif_probs: ProbsTypes) -> numpy.ndarray:
    """returns array of motif probs in alphabet order

    Parameters
    ----------
    alphabet
        ordered series of states
    motif_probs
        dict or other series with values present in alphabet

    Raises
    ------
    AssertionError if values do not sum to 1

    Returns
    -------
    array of floats in alphabet order
    """

    if len(alphabet) != len(motif_probs):
        msg = f"Can't match {len(motif_probs)} probs to {len(alphabet)} alphabet"
        raise ValueError(
            msg,
        )

    if hasattr(motif_probs, "keys"):
        # need to wrap the keys method because it can be a DictArray
        if diff := set(motif_probs.keys()) - set(alphabet):
            msg = f"Can't find motif(s) {diff!r} in alphabet"
            raise ValueError(msg)
        return numpy.array([motif_probs[motif] for motif in alphabet])
    motif_probs = numpy.asarray(motif_probs)

    numpy.testing.assert_allclose(
        motif_probs.sum(),
        1,
        err_msg=f"does not summ to 1 {motif_probs}",
    )
    return motif_probs
