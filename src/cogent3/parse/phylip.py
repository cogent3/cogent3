"""Parser for the PHYLIP sequence alignment format.

Supports relaxed PHYLIP (whitespace-delimited labels of arbitrary length) by
default and strict PHYLIP (labels in a fixed 10-character column) via
``strict_mode=True``. Interleaved and sequential layouts are auto-detected from
the line structure when ``interleaved`` is not given.
"""

import typing

import cogent3
from cogent3.parse.fasta import OptConverterType, OutTypes, minimal_converter
from cogent3.parse.record import RecordError

ConverterType = typing.Callable[[bytes], OutTypes]

_LABEL_WIDTH = 10


def _parse_header(line: str) -> tuple[int, int]:
    """Return the number of sequences and the alignment length from the header."""
    parts = line.split()
    return int(parts[0]), int(parts[1])


def _is_labelled(line: str, strict_mode: bool) -> bool:
    """Whether a data line carries a label rather than continuing a sequence."""
    if strict_mode:
        return bool(line[:_LABEL_WIDTH].strip())
    return not line[:1].isspace()


def _split_line(line: str, strict_mode: bool) -> tuple[str | None, str]:
    """Split a data line into its label (or None) and raw sequence fragment."""
    if strict_mode:
        label = line[:_LABEL_WIDTH].strip()
        return (label or None), line[_LABEL_WIDTH:]
    if line[:1].isspace():
        return None, line
    parts = line.split(None, 1)
    return parts[0], parts[1] if len(parts) > 1 else ""


def _detect_layout(lines: list[str], num_seqs: int, strict_mode: bool) -> bool:
    """Whether the alignment is interleaved, inferred from the line structure."""
    if num_seqs == 1:
        return False
    flags = [_is_labelled(line, strict_mode) for line in lines]
    if len(flags) >= 2 and flags[0] and not flags[1]:
        return False
    return all(flags[:num_seqs])


def _emit(
    label: str,
    fragments: list[str],
    seq_len: int,
    converter: ConverterType,
    id_map: dict[str, str] | None,
) -> tuple[str, OutTypes]:
    """Convert assembled fragments to a sequence, validating its length."""
    seq = converter("".join(fragments).encode("utf8"))
    if len(seq) != seq_len:
        msg = (
            f"Length of sequence {label!r} is not the same as in header "
            f"Found {len(seq)}, Expected {seq_len}"
        )
        raise RecordError(msg)
    if id_map is not None:
        label = id_map.get(label, label)
    return label, seq


def _parse_interleaved(
    lines: list[str],
    num_seqs: int,
    seq_len: int,
    strict_mode: bool,
    converter: ConverterType,
    id_map: dict[str, str] | None,
) -> typing.Iterator[tuple[str, OutTypes]]:
    """Yield (label, seq) from interleaved phylip data."""
    labels: list[str] = []
    fragments: list[list[str]] = []
    for index, line in enumerate(lines):
        if index < num_seqs:
            label, fragment = _split_line(line, strict_mode)
            labels.append(label or "")
            fragments.append([fragment])
        else:
            fragments[index % num_seqs].append(line)
    for label, parts in zip(labels, fragments, strict=True):
        yield _emit(label, parts, seq_len, converter, id_map)


def _parse_sequential(
    lines: list[str],
    seq_len: int,
    strict_mode: bool,
    converter: ConverterType,
    id_map: dict[str, str] | None,
) -> typing.Iterator[tuple[str, OutTypes]]:
    """Yield (label, seq) from sequential phylip data."""
    label: str | None = None
    fragments: list[str] = []
    for line in lines:
        this_label, fragment = _split_line(line, strict_mode)
        if this_label is not None:
            if label is not None:
                yield _emit(label, fragments, seq_len, converter, id_map)
            label = this_label
            fragments = [fragment]
        else:
            fragments.append(fragment)
    if label is not None:
        yield _emit(label, fragments, seq_len, converter, id_map)


def MinimalPhylipParser(
    data: typing.Iterable[str],
    id_map: dict[str, str] | None = None,
    interleaved: bool | None = None,
    strict_mode: bool = False,
    converter: OptConverterType = None,
) -> typing.Iterator[tuple[str, OutTypes]]:
    """Yield successive sequences from phylip data as (label, seq) tuples.

    Parameters
    ----------
    data
        an iterable of lines in phylip format
    id_map
        optional mapping from parsed labels to replacement labels
    interleaved
        force interleaved (True) or sequential (False) layout; when None the
        layout is auto-detected from the line structure
    strict_mode
        when True labels occupy a fixed 10-character column; when False labels
        are whitespace-delimited and may be of arbitrary length
    converter
        callable mapping raw sequence bytes to the yielded sequence; when None
        uses a converter that removes whitespace and upper-cases the residues
    """
    lines = [line for line in data if line.strip()]
    if not lines:
        return
    num_seqs, seq_len = _parse_header(lines[0])
    data_lines = lines[1:]
    if not num_seqs or not seq_len or not data_lines:
        return
    if converter is None:
        converter = minimal_converter()
    if interleaved is None:
        interleaved = _detect_layout(data_lines, num_seqs, strict_mode)
    if interleaved:
        yield from _parse_interleaved(
            data_lines, num_seqs, seq_len, strict_mode, converter, id_map
        )
    else:
        yield from _parse_sequential(
            data_lines, seq_len, strict_mode, converter, id_map
        )


def get_align_for_phylip(
    data: typing.Iterable[str],
    id_map: dict[str, str] | None = None,
    strict_mode: bool = False,
) -> "cogent3.core.alignment.Alignment":
    """Return an Alignment object from phylip data.

    Parameters
    ----------
    data
        an iterable of lines in phylip format
    id_map
        optional mapping from parsed labels to replacement labels
    strict_mode
        when True labels occupy a fixed 10-character column
    """
    tuples = list(MinimalPhylipParser(data, id_map, strict_mode=strict_mode))
    return cogent3.make_aligned_seqs(tuples, moltype="text")
