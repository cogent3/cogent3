# Per-sequence motif-count benchmark

Run this script using the Python environment with this checkout installed in
editable mode. It compares the current `Alignment.counts_per_seq()` and
`SequenceCollection.counts_per_seq()` methods with the method definitions from a
recorded local Git commit.

```sh
python benchmarks/counts_per_seq.py --repeats 7 --output counts-per-seq-run-01.json
```

The default baseline is `37cda7a777689a13c4e17e71bbca679142d6f15f`. To compare a
different compatible commit, pass `--baseline <commit-or-ref>`. The requested ref
must exist locally. Only use a trusted baseline: the script compiles and executes
its two method definitions in an isolated namespace. It does not modify source
files, replace installed modules, or check out another revision. The other code
and dependencies remain those of the candidate installation, so this measures
the change in these methods rather than comparing two complete installations.

The 48 cases cover aligned and unaligned canonical DNA, ambiguous DNA, protein,
arbitrary bytes, and ragged unaligned DNA. Motif lengths are 1, 2, and 3.
Canonical DNA includes both `exclude_unobserved=True` and the default
`exclude_unobserved=False`, with the default gap and ambiguity filters. Protein
and bytes use observed columns only, avoiding an unnecessarily large Cartesian
alphabet. Ragged inputs include an empty row and lengths from 3 to 30,000; the
known baseline failure for nonempty rows shorter than a motif is covered by the
regression tests instead of being timed.

Data are generated using NumPy's random generator with seed 2215. Before timing,
the script warms both methods and requires identical counts, row and column
labels, shapes, and integer dtypes. It alternates baseline/candidate order on
successive repetitions, records every elapsed time, and reports the ratio of
median latencies (`baseline / candidate`; above 1 means faster). Data generation,
collection construction, correctness checks, and warm-up are excluded.

The JSON includes the resolved baseline commit, candidate commit and source-file
hash, whether that file differs from HEAD, runtime/package versions, OS and CPU
metadata, per-row input lengths, input/output dtypes, parameters, output shape,
and all timing samples. An existing output file is never overwritten. Use
`--repeats 1` only for a smoke run, and choose a new output path for each run.

These are synthetic workloads, not a guarantee of performance for every dataset.
No peak-memory measurement, CPU affinity, or process isolation is performed.
Short cases are sensitive to timer and scheduling noise; inspect the raw samples
and repeat substantial comparisons on a quiet machine. The benchmark requires
Git and the project's normal runtime dependencies, with no extra benchmark
package required.
