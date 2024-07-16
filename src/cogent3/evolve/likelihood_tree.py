"""Leaf and Edge classes that can calculate their likelihoods.
Each leaf holds a sequence.  Used by a likelihood function."""

import numpy

from . import likelihood_tree_numba as likelihood_tree


numpy.seterr(all="ignore")


class _LikelihoodTreeEdge(object):
    def __init__(self, children, edge_name, alignment=None):
        self.edge_name = edge_name
        self.alphabet = children[0].alphabet

        M = children[0].shape[-1]
        for child in children:
            assert child.shape[-1] == M

        # Unique positions are unique combos of input positions
        if alignment is None:
            # The children are pre-aligned gapped sequences
            assignments = [c.index for c in children]
        else:
            self.alignment = alignment
            # The children are ungapped sequences, 'alignment'
            # indicates where gaps need to go.
            assignments = []
            for i, c in enumerate(children):
                a = []
                for align_index in alignment:
                    col = align_index[i]
                    if col is None:
                        u = len(c.uniq) - 1  # gap
                    else:
                        u = c.index[col]
                        assert 0 <= u < len(c.uniq) - 1, (
                            u,
                            len(c.uniq),
                            c.uniq[-1],
                            align_index,
                        )
                    a.append(u)
                assignments.append(a)
        (uniq, counts, self.index) = _indexed(list(zip(*assignments)))

        # extra column for gap
        uniq.append(tuple([len(c.uniq) - 1 for c in children]))
        counts.append(0)

        self.uniq = numpy.asarray(uniq, int)

        # For faster math, a contiguous index array for each child
        self.indexes = numpy.ascontiguousarray(
            [numpy.array(list(ch), int) for ch in numpy.transpose(self.uniq)]
        )

        # If this is the root it will need to weight the total
        # log likelihoods by these counts:
        self.counts = numpy.array(counts, float)

        # For product of child likelihoods
        self._indexed_children = list(zip(self.indexes, children))
        self.shape = [len(self.uniq), M]

        # Derive per-column degree of ambiguity from children's
        ambigs = [child.ambig[index] for (index, child) in self._indexed_children]
        self.ambig = numpy.prod(ambigs, axis=0)

    def get_site_patterns(self, cols):
        # Recursive lookup of Site Patterns aka Alignment Columns
        child_motifs = [
            child.get_site_patterns(index[cols])
            for (index, child) in self._indexed_children
        ]
        return ["".join(child[u] for child in child_motifs) for u in range(len(cols))]

    def restrict_motif(self, input_likelihoods, fixed_motif):
        # for reconstruct_ancestral_seqs
        mask = numpy.zeros([input_likelihoods.shape[-1]], float)
        mask[fixed_motif] = 1.0
        input_likelihoods *= mask

    def select_columns(self, cols):
        children = []
        for index, child in self._indexed_children:
            child = child.select_columns(cols)
            children.append(child)
        return self.__class__(children, self.edge_name)

    def get_full_length_likelihoods(self, likelihoods):
        return likelihoods[self.index]

    def calc_G_statistic(self, likelihoods, return_table=False):
        # A Goodness-of-fit statistic
        from cogent3.util.table import Table

        unambig = (self.ambig == 1.0).nonzero()[0]
        observed = self.counts[unambig].astype(int)
        expected = likelihoods[unambig] * observed.sum()
        # chisq = ((observed-expected)**2 / expected).sum()
        G = 2 * observed.dot(numpy.log(observed / expected))

        if return_table:
            motifs = self.get_site_patterns(unambig)
            rows = list(zip(motifs, observed, expected))
            rows.sort(key=lambda row: (-row[1], row[0]))
            table = Table(
                header=["Pattern", "Observed", "Expected"],
                data=rows,
                index_name="Pattern",
            )
            return (G, table)
        else:
            return G

    def get_edge(self, name):
        if self.edge_name == name:
            return self
        else:
            for i, c in self._indexed_children:
                r = c.get_edge(name)
                if r is not None:
                    return r
        return None

    def make_partial_likelihoods_array(self):
        return numpy.ones(self.shape, float)

    def sum_input_likelihoods(self, *likelihoods):
        result = numpy.ones(self.shape, float)
        self.sum_input_likelihoodsR(result, *likelihoods)
        return result

    def as_leaf(self, likelihoods):
        assert len(likelihoods) == len(self.counts)
        return LikelihoodTreeLeaf(
            likelihoods,
            likelihoods,
            self.counts,
            self.index,
            self.edge_name,
            self.alphabet,
            None,
        )


class LikelihoodTreeEdge(_LikelihoodTreeEdge):
    # Should be a subclass of regular tree edge?

    # For scaling very very small numbers
    BASE = 2.0**100
    LOG_BASE = numpy.log(BASE)

    def sum_input_likelihoodsR(self, result, *likelihoods):
        if not self.indexes.flags["C_CONTIGUOUS"]:
            self.indexes = numpy.ascontiguousarray(self.indexes)
        if not result.flags["C_CONTIGUOUS"]:
            result = numpy.ascontiguousarray(result)
        return likelihood_tree.sum_input_likelihoods(
            self.indexes,
            result,
            likelihoods,
        )

    # For root

    def log_dot_reduce(self, patch_probs, switch_probs, plhs):
        exponent = 0
        state_probs = patch_probs.copy()
        for site in self.index:
            state_probs = numpy.dot(switch_probs, state_probs) * plhs[site]
            while max(state_probs) < 1.0:
                state_probs *= self.BASE
                exponent -= 1
        return numpy.log(sum(state_probs)) + exponent * self.LOG_BASE

    def get_total_log_likelihood(self, input_likelihoods, mprobs):
        lhs = likelihood_tree.inner_product(input_likelihoods, mprobs)
        return self.get_log_sum_across_sites(lhs)

    def get_log_sum_across_sites(self, lhs):
        return likelihood_tree.get_log_sum_across_sites(lhs, self.counts)


def _indexed(values):
    # >>> _indexed(['a', 'b', 'c', 'a', 'a'])
    # (['a', 'b', 'c'], [3, 1, 1], [0, 1, 2, 0, 0])
    index = numpy.zeros([len(values)], int)
    unique = []
    counts = []
    seen = {}
    for c, key in enumerate(values):
        if key in seen:
            i = seen[key]
            counts[i] += 1
        else:
            i = len(unique)
            unique.append(key)
            counts.append(1)
            seen[key] = i
        index[c] = i
    return unique, counts, index


def get_matched_array(alphabet, moltype, motifs, dtype=float) -> numpy.ndarray:
    """Returns an array in which rows are motifs, columns are items in self.

    Result is an array of dtype in which a[i][j] indicates whether the ith
    motif passed in as motifs is a symbol that matches the jth character
    in self. For example, on the DNA alphabet 'TCAG', the degenerate symbol
    'Y' would correspond to the row [1,1,0,0] because Y is a degenerate
    symbol that encompasses T and C but not A or G.
    """
    result = numpy.zeros([len(motifs), len(alphabet)], dtype)
    obj_to_index = alphabet.to_indices
    for u, ambig_motif in enumerate(motifs):
        for motif in moltype.resolve_ambiguity(ambig_motif, alphabet=alphabet):
            result[u, obj_to_index((motif,))] = 1.0
    return result


def make_likelihood_tree_leaf(sequence, alphabet, seq_name):
    motif_len = alphabet.get_motif_len()
    sequence2 = sequence.get_in_motif_size(motif_len)

    # Convert sequence to indexed list of unique motifs
    uniq_motifs, counts, index = _indexed(sequence2)

    # extra column for gap
    uniq_motifs.append("?" * motif_len)
    counts.append(0)

    counts = numpy.array(counts, float)

    # Convert list of unique motifs to array of unique profiles
    try:
        likelihoods = get_matched_array(
            alphabet, alphabet.moltype, uniq_motifs, dtype=float
        )
    except alphabet.AlphabetError as detail:
        motif = str(detail)
        posn = list(sequence2).index(motif) * motif_len
        raise ValueError(
            f"{motif!r} at {seq_name!r}:{posn} not in alphabet"
        ) from detail

    return LikelihoodTreeLeaf(
        uniq_motifs, likelihoods, counts, index, seq_name, alphabet, sequence
    )


class LikelihoodTreeLeaf(object):
    def __init__(self, uniq, likelihoods, counts, index, edge_name, alphabet, sequence):
        if sequence is not None:
            self.sequence = sequence
        self.alphabet = alphabet
        self.name = self.edge_name = edge_name
        self.uniq = uniq
        self.motifs = numpy.asarray(uniq)
        self.input_likelihoods = likelihoods
        self.counts = counts
        self.index = index
        self.shape = likelihoods.shape
        self.ambig = numpy.sum(self.input_likelihoods, axis=-1)

    def backward(self):
        index = numpy.array(self.index[::-1, ...])
        return self.__class__(
            self.uniq,
            self.input_likelihoods,
            self.counts,
            index,
            self.edge_name,
            self.alphabet,
            None,
        )

    def __len__(self):
        return len(self.index)

    def __getitem__(self, index):
        cols = list(range(*index.indices(len(self.index))))
        return self.select_columns(cols)

    def get_motif_counts(self, include_ambiguity=False):
        weights = self.counts / self.ambig
        profile = self.input_likelihoods * weights[..., numpy.newaxis]
        if not include_ambiguity:
            unambig = self.ambig == 1.0
            profile = numpy.compress(unambig, profile, axis=0)
        return numpy.sum(profile, axis=0)

    def get_ambiguous_positions(self):
        ambig = {}
        for i, u in enumerate(self.index):
            if self.ambig[u] != 1.0:
                ambig[i] = self.uniq[u]
        return ambig

    def select_columns(self, cols):
        sub_index = [self.index[i] for i in cols]
        (keep, counts, index) = _indexed(sub_index)
        keep.append(len(self.uniq) - 1)  # extra column for gap
        counts.append(0)
        counts = numpy.array(counts, float)
        uniq = [self.uniq[u] for u in keep]
        likelihoods = self.input_likelihoods[keep]
        return self.__class__(
            uniq, likelihoods, counts, index, self.edge_name, self.alphabet, None
        )

    def get_edge(self, name):
        if self.edge_name == name:
            return self
        else:
            return None

    def get_site_patterns(self, cols):
        return numpy.asarray(self.uniq)[cols]
