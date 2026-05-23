from __future__ import annotations

from collections.abc import Callable, Iterable

import numpy
from scinexus.io_util import PathType, iter_line_blocks

from cogent3.parse.record import RecordError

OutTypes = str | bytes | numpy.ndarray
OptConverterType = Callable[[bytes], OutTypes] | None


def _process_fastq_block(
    block: list[str],
    converter: Callable[[bytes], OutTypes],
    qual_converter: Callable[[bytes], OutTypes],
) -> tuple[str, OutTypes, OutTypes]:
    if len(block) != 4:
        msg = f"fastq record must have 4 lines, got {len(block)}"
        raise RecordError(msg)
    label_line, seq_line, plus_line, qual_line = block
    if not label_line.startswith("@"):
        msg = f"label line must start with '@', got {label_line!r}"
        raise RecordError(msg)
    if not plus_line.startswith("+"):
        msg = f"separator line must start with '+', got {plus_line!r}"
        raise RecordError(msg)
    label = label_line[1:].strip()
    return (
        label,
        converter(seq_line.encode("utf8")),
        qual_converter(qual_line.encode("utf8")),
    )


def iter_fastq_records(
    data: PathType,
    *,
    converter: OptConverterType = None,
    qual_converter: OptConverterType = None,
    chunk_size: int | None = 5_000_000,
) -> Iterable[tuple[str, OutTypes, OutTypes]]:
    """yields (label, sequence, quality) tuples from a fastq source.

    Parameters
    ----------
    data
        path to a fastq file, or any path-like value
    converter
        applied to the sequence line as bytes. If None, the line is decoded
        to str.
    qual_converter
        applied to the quality line as bytes. If None, the line is decoded
        to str.
    chunk_size
        size in bytes of chunks read from the file path.

    Returns
    -------
    iterable of (label, sequence, quality) tuples, with sequence and
    quality typed by their converters.

    Notes
    -----
    Use cogent3.core.alphabet.make_qual_converter to obtain a numerical
    quality converter that maps Phred+33 or Phred+64 ASCII bytes to a
    numpy.uint8 array of quality scores. Sequence converters can be
    built from cogent3.core.alphabet.make_text_to_array_converter.
    """
    seq_conv = converter if converter is not None else bytes.decode
    qual_conv = qual_converter if qual_converter is not None else bytes.decode
    for block in iter_line_blocks(data, num_lines=4, chunk_size=chunk_size):
        yield _process_fastq_block(block, seq_conv, qual_conv)
