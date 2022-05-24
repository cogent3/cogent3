#!/usr/bin/env python
"""Random sequences and random evolution of sequences in a tree"""

import bisect

import numpy


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def argpicks(freqs, random_series):
    partition = numpy.add.accumulate(freqs)
    assert abs(partition[-1] - 1.0) < 1e-6, (freqs, partition)
    while True:
        x = random_series.uniform(0.0, 1.0)
        i = bisect.bisect_left(partition, x)
        yield i


def argpick(freqs, random_series):
    return next(argpicks(freqs, random_series))


def _randomMotifGenerator(random_series, motif_probs):
    motifs = list(motif_probs.keys())
    freqs = [motif_probs[m] for m in motifs]
    for i in argpicks(freqs, random_series):
        yield motifs[i]


def evolve_sequence(
    random_series, motifs, parent_seq, site_cats, psubs, preserved_sites=()
):
    """Evolve a new sequence derived from parent_seq.  Uses psubs[site_cats[i]]
    to pick a new motif derived from parent_seq[i]"""
    seq = []
    randomMotifSources = {}
    for (i, parent_motif) in enumerate(parent_seq):
        if i in preserved_sites:
            edge_motif = preserved_sites[i]
        else:
            if parent_motif not in randomMotifSources:
                mprobs = {}
                parent_motif_index = motifs.index(parent_motif)
                site_cat = site_cats[i]
                psub = psubs[site_cat]
                for (dest_motif_index, dest_motif) in enumerate(motifs):
                    prob = psub[parent_motif_index, dest_motif_index]
                    mprobs[dest_motif] = prob
                randomMotifSources[site_cat, parent_motif] = _randomMotifGenerator(
                    random_series, mprobs
                )
            edge_motif = next(randomMotifSources[site_cat, parent_motif])
        seq.append(edge_motif)
    return seq


def random_sequence(random_series, motif_probs, sequence_length):
    getRootRandomMotif = _randomMotifGenerator(random_series, motif_probs).__next__
    return [getRootRandomMotif() for i in range(sequence_length)]


class AlignmentEvolver(object):
    # Encapsulates settings that are constant throughout the recursive generation
    # of a synthetic alignment.

    def __init__(
        self,
        random_series,
        orig_ambig,
        exclude_internal,
        bin_names,
        site_bins,
        psub_for,
        motifs,
    ):
        self.random_series = random_series
        self.orig_ambig = orig_ambig
        self.exclude_internal = exclude_internal
        self.bin_names = bin_names
        self.site_bins = site_bins
        self.psub_for = psub_for
        self.motifs = motifs

    def __call__(self, tree, root_sequence):
        # probsd = dict(enumerate(self.bin_probs))
        # bprobs = _randomMotifGenerator(self.random_series, probsd)
        # site_bins = [bprobs.next() for c in range(len(root_sequence))]
        return self.generate_simulated_seqs(tree, root_sequence)

    def generate_simulated_seqs(self, parent, parent_seq):
        """recursively generate the descendant sequences by descending the tree
        from root.
        Each child will be set by mutating the parent motif based on the probs
        in the psub matrix of this edge.

        random_series - get a random numer 0-1 by calling random_series.random()
        length - the desired alignment length
        parent - the edge structure.
        parent_seq - the corresponding sequence. This will be mutated for each
        of its children, based on their psub matricies.
        """

        # This depends on parameter names 'mprobs', 'alignment2', 'bprobs' and
        # 'psubs'.  Might be better to integrate it into likelihood_calculation.

        if self.exclude_internal and parent.children:
            simulated_sequences = {}
        else:
            simulated_sequences = {parent.name: "".join(parent_seq)}

        for edge in parent.children:
            # The result for this edge - a list of motifs

            # Keep original ambiguity codes
            if edge.name in self.orig_ambig:
                orig_seq_ambig = self.orig_ambig[edge.name]
            else:
                orig_seq_ambig = {}

            # Matrix of substitution probabilities
            psubs = [self.psub_for(edge.name, bin) for bin in self.bin_names]

            # Make the semi-random sequence for this edge.
            edge_seq = evolve_sequence(
                self.random_series,
                self.motifs,
                parent_seq,
                self.site_bins,
                psubs,
                orig_seq_ambig,
            )

            # Pass this new edge sequence on down the tree
            descendant_sequences = self.generate_simulated_seqs(edge, edge_seq)
            simulated_sequences.update(descendant_sequences)

        return simulated_sequences
