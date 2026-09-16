### Enhancements

- Improve `counts_per_seq()` for sequence collections and alignments using
  array-based monomer counting and fixed-width records for non-overlapping
  motifs, without allocating an additional full k-mer counting alphabet.

### Bug fixes

- Return zero counts for a sequence shorter than `motif_length` in a ragged
  collection, instead of raising an empty-axis error.

### Contributors

- Qinzi27
