"""Benchmark per-sequence counts against methods from a trusted local Git ref.

Run from a development installation with --output <new.json>.
Input construction and the correctness/warm-up calls are excluded from timings.
"""

# This is a standalone developer script, not an importable package.
# ruff: noqa: INP001

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.metadata
import json
import logging
import os
import platform
import shutil
import statistics
import subprocess
import time
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy

from cogent3 import make_aligned_seqs, make_unaligned_seqs
from cogent3.core import alignment

if TYPE_CHECKING:
    from collections.abc import Callable

    from cogent3.core.profile import MotifCountsArray

REPO = Path(__file__).resolve().parents[1]
DEFAULT_BASELINE = "37cda7a777689a13c4e17e71bbca679142d6f15f"
SEED = 2215
LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class Scenario:
    name: str
    moltype: str
    symbols: str | bytes
    lengths: tuple[int, ...]
    include_ambiguity: bool = False
    allow_gap: bool = False
    exclude_unobserved: tuple[bool, ...] = (True,)


SCENARIOS = (
    Scenario("dna-small", "dna", "ACGT", (60,) * 8, exclude_unobserved=(True, False)),
    Scenario("dna-long", "dna", "ACGT", (30000,) * 8, exclude_unobserved=(True, False)),
    Scenario("dna-ambiguous", "dna", "ACGTN?-", (30000,) * 8, True, True),
    Scenario("protein", "protein", "ACDEFGHIKLMNPQRSTVWYBX?-", (6000,) * 8, True, True),
    Scenario(
        "bytes",
        "bytes",
        bytes([0, 1, 45, 63, 127, 128, 254, 255]),
        (6000,) * 8,
        True,
        True,
    ),
    Scenario(
        "dna-ragged",
        "dna",
        "ACGT",
        (0, 3, 60, 600, 3000, 9000, 15000, 30000),
        exclude_unobserved=(True, False),
    ),
)


def git_output(*args: str) -> str:
    """Read repository provenance without changing the checkout or Git config."""
    executable = shutil.which("git")
    if executable is None:
        message = "Git is required to read the baseline methods."
        raise FileNotFoundError(message)
    # No shell is used; the executable is resolved and arguments are explicit.
    return subprocess.check_output(  # noqa: S603
        [executable, "-c", f"safe.directory={REPO.as_posix()}", *args],
        cwd=REPO,
        text=True,
        encoding="utf-8",
    ).strip()


def baseline_methods(
    commit: str,
) -> dict[str, Callable[..., MotifCountsArray | None]]:
    """Compile only counts_per_seq methods from a caller-trusted Git commit."""
    source = git_output("show", f"{commit}:src/cogent3/core/alignment.py")
    tree = ast.parse(source)
    namespace = vars(alignment).copy()
    methods = {}
    for node in tree.body:
        if not isinstance(node, ast.ClassDef) or node.name not in (
            "SequenceCollection",
            "Alignment",
        ):
            continue
        method = next(
            item
            for item in node.body
            if isinstance(item, ast.FunctionDef) and item.name == "counts_per_seq"
        )
        method.name = f"baseline_{node.name}"
        module = ast.Module(body=[method], type_ignores=[])
        # Executing these explicitly requested local definitions is the comparison
        # mechanism. Never pass a ref containing code you do not trust.
        exec(compile(module, f"{commit}/{node.name}", "exec"), namespace)  # noqa: S102
        methods[node.name] = namespace[method.name]
    return methods


def check_equal(old: MotifCountsArray | None, new: MotifCountsArray | None) -> None:
    """Reject timings unless labels, values, dimensions, and dtypes agree."""
    if old is None or new is None:
        if old is not new:
            message = "Baseline and candidate differ in whether they return None."
            raise AssertionError(message)
        return
    numpy.testing.assert_equal(old.to_dict(), new.to_dict())
    numpy.testing.assert_equal(old.template.names, new.template.names)
    numpy.testing.assert_equal(old.motifs, new.motifs)
    numpy.testing.assert_array_equal(old.array, new.array)
    numpy.testing.assert_equal(old.array.dtype, new.array.dtype)


def make_data(
    scenario: Scenario, rng: numpy.random.Generator
) -> dict[str, str | bytes]:
    """Generate fixed-seed inputs without making array encoding part of timing."""
    data: dict[str, str | bytes] = {}
    for index, length in enumerate(scenario.lengths):
        values = rng.choice(list(scenario.symbols), size=length).tolist()
        data[f"seq-{index}"] = (
            bytes(values) if scenario.moltype == "bytes" else "".join(values)
        )
    return data


def time_pair(
    old: Callable[[], MotifCountsArray | None],
    new: Callable[[], MotifCountsArray | None],
    repeats: int,
) -> dict[str, list[float]]:
    """Alternate baseline/candidate execution order to reduce ordering bias."""
    times: dict[str, list[float]] = {"baseline": [], "candidate": []}
    for repeat in range(repeats):
        order = [("baseline", old), ("candidate", new)]
        if repeat % 2:
            order.reverse()
        for name, function in order:
            start = time.perf_counter()
            function()
            times[name].append(time.perf_counter() - start)
    return times


def benchmark_case(
    collection: alignment.Alignment | alignment.SequenceCollection,
    baseline: Callable[..., MotifCountsArray | None],
    scenario: Scenario,
    motif_length: int,
    exclude_unobserved: bool,
    repeats: int,
) -> dict[str, object]:
    parameters = {
        "motif_length": motif_length,
        "include_ambiguity": scenario.include_ambiguity,
        "allow_gap": scenario.allow_gap,
        "exclude_unobserved": exclude_unobserved,
    }
    old = partial(baseline, collection, **parameters)
    new = partial(collection.counts_per_seq, **parameters)
    # Warm both paths and check public results before collecting any timings.
    previous, candidate = old(), new()
    check_equal(previous, candidate)
    times = time_pair(old, new, repeats)
    old_median = statistics.median(times["baseline"])
    new_median = statistics.median(times["candidate"])
    entry = {
        "scenario": scenario.name,
        "kind": type(collection).__name__,
        "moltype": scenario.moltype,
        "rows": len(scenario.lengths),
        "lengths": scenario.lengths,
        "total_characters": sum(scenario.lengths),
        "parameters": parameters,
        "input_dtypes": sorted(
            {str(numpy.array(seq).dtype) for seq in collection.seqs}
        ),
        "baseline_output_dtype": None
        if previous is None
        else str(previous.array.dtype),
        "candidate_output_dtype": None
        if candidate is None
        else str(candidate.array.dtype),
        "output_shape": None if candidate is None else list(candidate.array.shape),
        "times_seconds": times,
        "baseline_median_seconds": old_median,
        "candidate_median_seconds": new_median,
        "speedup": old_median / new_median,
    }
    LOGGER.info(
        "%s %s k=%s exclude_unobserved=%s: %.2fx",
        scenario.name,
        type(collection).__name__,
        motif_length,
        exclude_unobserved,
        entry["speedup"],
    )
    return entry


def run_benchmarks(commit: str, repeats: int) -> list[dict[str, object]]:
    baseline = baseline_methods(commit)
    rng = numpy.random.default_rng(SEED)
    results = []
    for scenario in SCENARIOS:
        data = make_data(scenario, rng)
        makers = [make_unaligned_seqs]
        if len(set(scenario.lengths)) == 1:
            makers.append(make_aligned_seqs)
        for maker in makers:
            collection = maker(data, moltype=scenario.moltype)
            for motif_length in (1, 2, 3):
                results.extend(
                    benchmark_case(
                        collection,
                        baseline[type(collection).__name__],
                        scenario,
                        motif_length,
                        exclude,
                        repeats,
                    )
                    for exclude in scenario.exclude_unobserved
                )
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--baseline", default=DEFAULT_BASELINE)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be at least 1")
    if args.output.exists():
        raise FileExistsError(args.output)
    source_path = Path(cast("str", alignment.__file__)).resolve()
    if source_path != (REPO / "src/cogent3/core/alignment.py").resolve():
        parser.error(
            "Install this checkout in editable mode before running the benchmark."
        )
    commit = git_output(
        "rev-parse", "--verify", "--end-of-options", f"{args.baseline}^{{commit}}"
    )
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    results = run_benchmarks(commit, args.repeats)
    report = {
        "baseline_commit": commit,
        "baseline_requested": args.baseline,
        "candidate_head": git_output("rev-parse", "HEAD"),
        "candidate_alignment_sha256": hashlib.sha256(
            source_path.read_bytes()
        ).hexdigest(),
        "candidate_alignment_modified": bool(
            git_output("diff", "HEAD", "--", str(source_path))
        ),
        "seed": SEED,
        "repeats": args.repeats,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "logical_cpus": os.cpu_count(),
        "python": platform.python_version(),
        "versions": {
            package: importlib.metadata.version(package)
            for package in ("cogent3", "numpy", "scipy", "numba", "scinexus")
        },
        "limitations": [
            "End-to-end method latency excludes input construction and correctness/warm-up calls.",
            "Baseline methods use the candidate checkout's dependencies and remaining code.",
            "Synthetic fixed-seed inputs do not represent all biological workloads.",
            "No peak-memory measurement, CPU affinity, or process isolation is performed.",
            "Short timings and results on a shared machine can be noisy; inspect all repeats.",
        ],
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation also protects against another process choosing this path.
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    LOGGER.info("Saved %s cases to %s", len(results), args.output)


if __name__ == "__main__":
    main()
