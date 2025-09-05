"""Parsers for FASTA and related formats."""

import io
import os
import pathlib
import re
import string
import typing
from collections.abc import Callable
from functools import singledispatch

import numpy

from cogent3.core.info import Info
from cogent3.parse.record import RecordError
from cogent3.parse.record_finder import LabeledRecordFinder
from cogent3.util.io import is_url, open_

_white_space = re.compile(r"\s+")
strip = str.strip


OutTypes = typing.Union[str, bytes, numpy.ndarray]
OptConverterType = typing.Optional[typing.Callable[[bytes], OutTypes]]
RenamerType = typing.Callable[[str], str]


PathOrIterableType = typing.Union[os.PathLike, list[str], tuple[str]]


@singledispatch
def _prep_data(data) -> typing.Iterable[str]:
    return data


@_prep_data.register
def _(data: str):
    with open_(data) as infile:
        return infile.read().splitlines()


@_prep_data.register
def _(data: os.PathLike):
    return _prep_data(str(data))


def _faster_parser(
    data: typing.Iterable[str],
    label_to_name: RenamerType,
    label_char: str,
) -> typing.Iterable[tuple[str, str]]:
    label = None
    seq = []
    for line in data:
        if not line:
            # ignore empty lines
            continue

        if line[0] in label_char:
            if seq:
                yield label_to_name(label or ""), _white_space.sub("", "".join(seq))

            label = line[1:].strip()
            seq = []
        else:
            seq.append(line.strip())

    if seq:
        yield label_to_name(label or ""), _white_space.sub("", "".join(seq))


def _strict_parser(
    data: typing.Iterable[str],
    label_to_name: RenamerType,
    label_char: str,
) -> typing.Iterable[tuple[str, str]]:
    seq: list[str] = []
    label: str | None = None
    for line in data:
        if not line or line.startswith("#"):
            # ignore empty or comment lines
            continue

        if line[0] in label_char:
            if label is not None:
                if not seq:
                    msg = f"{label} has no data"
                    raise RecordError(msg)
                yield label_to_name(label), _white_space.sub("", "".join(seq))
            elif seq:
                msg = "missing a label"
                raise RecordError(msg)

            label = line[1:].strip()
            seq = []
        else:
            seq.append(line.strip())

    if not seq:
        msg = f"{label} has no data"
        raise RecordError(msg)
    if label is None:
        msg = "missing a label"
        raise RecordError(msg)

    yield label_to_name(label), _white_space.sub("", "".join(seq))


def MinimalFastaParser(
    path: PathOrIterableType,
    strict: bool = True,
    label_to_name: RenamerType = str,
    label_characters: str = ">",
) -> typing.Iterable[tuple[str, str]]:
    """
    Yields successive sequences from infile as (label, seq) tuples.

    Parameters
    ----------
    path
    strict
        raises RecordError when label or seq missing
    label_to_name
        function for converting a label to a name
    label_characters
        character(s) at the start of a label line
    """
    if not path:
        return

    data = _prep_data(path)
    label_char = set(label_characters)
    if strict:
        yield from _strict_parser(data, label_to_name, label_char)
    else:
        yield from _faster_parser(data, label_to_name, label_char)


GdeFinder = LabeledRecordFinder(lambda x: x.startswith(("#", "%")))


def MinimalGdeParser(infile, strict=True, label_to_name=str):
    return MinimalFastaParser(infile, strict, label_to_name, label_characters="%#")


# labeled fields in the NCBI FASTA records
NcbiLabels = {"dbj": "DDBJ", "emb": "EMBL", "gb": "GenBank", "ref": "RefSeq"}


class RichLabel(str):
    """Object for overloaded Fasta labels. Holds an Info object storing keyed
    attributes from the fasta label. The str is created from a provided format
    template that uses the keys from the Info object."""

    def __new__(cls, info, template="%s"):
        """

        Parameters
        ----------
        info
            a cogent3.core.info.info instance
        template
            a string template, using a subset of the keys in info.
            Defaults to just '%s'.

        Example:
            label = RichLabel(Info(name='rat', species='Rattus norvegicus'),
                        '%(name)s')"""
        label = template % info
        new = str.__new__(cls, label)
        new.info = info
        return new


def LabelParser(display_template, field_formatters, split_with=":", DEBUG=False):
    """returns a function for creating a RichLabel's from a string

    Parameters
    ----------
    display_template
        string format template
    field_formatters
        series of
        (field index, field name, coverter function)
    split_with
        characters separating fields in the label.
        The display_template must use at least one of the assigned field
        names.

    """
    indexed = False
    for _index, field, _converter in field_formatters:
        if field in display_template:
            indexed = True
    assert indexed, f"display_template [{display_template}] does not use a field name"
    sep = re.compile(f"[{split_with}]")

    def call(label):
        label = [label, label[1:]][label[0] == ">"]
        label = sep.split(label)
        if DEBUG:
            pass
        info = Info()
        for index, name, converter in field_formatters:
            if isinstance(converter, Callable):
                try:
                    info[name] = converter(label[index])
                except IndexError:
                    msg = f"parsing label {label} failed for property {name} at index {index}"
                    raise IndexError(
                        msg,
                    )
            else:
                info[name] = label[index]
        return RichLabel(info, display_template)

    return call


class minimal_converter:
    """coerces lower case bytes to upper case bytes and removes whitespace"""

    def __init__(self) -> None:
        lc = string.ascii_lowercase.encode("utf8")
        self._translate = b"".maketrans(lc, lc.upper())

    def __call__(self, text: bytes) -> str:
        return text.translate(self._translate, delete=b"\n\r\t ").decode("utf8")


@singledispatch
def iter_fasta_records(
    data,
    converter: OptConverterType = None,
    label_to_name: RenamerType = str,
) -> typing.Iterable[tuple[str, OutTypes]]:
    """generator returns labels and sequences converted from a fasta file

    Parameters
    ----------
    path
        location of the fasta file
    converter
        a callable that converts sequence characters as bytes, deleting
        unwanted characters (newlines, spaces). Whatever type this callable
        returns will be the type of the sequence returned. If None, uses
        minimal_converter() which returns bytes.
    label_to_name
        a callable that takes the sequence label as input and returns a new label.
        Defaults to the label itself.

    Returns
    -------
    the sequence label as a string and the sequence as transformed by converter
    """
    msg = f"iter_fasta_records not implemented for {type(data)}"
    raise TypeError(msg)


@iter_fasta_records.register
def _(
    data: bytes,
    converter: OptConverterType = None,
    label_to_name: RenamerType = str,
) -> typing.Iterable[tuple[str, OutTypes]]:
    if converter is None:
        converter = minimal_converter()

    records = data.split(b">")
    for record in records:
        if not len(record):
            continue

        eol = record.find(b"\n")
        if eol == -1:
            continue
        label = record[:eol].strip().decode("utf8")
        if label_to_name:
            label = label_to_name(label)
        seq = converter(record[eol + 1 :])
        yield label, seq


@iter_fasta_records.register
def _(data: str, converter: OptConverterType = None, label_to_name: RenamerType = str):
    if not is_url(data):
        try:
            os.stat(data)
        except OSError:
            msg = "data is a string but not a file path, directly provided data must be bytes"
            raise TypeError(
                msg,
            )

    with open_(data, mode="rb") as infile:
        data: bytes = infile.read()
    return iter_fasta_records(data, converter=converter, label_to_name=label_to_name)


@iter_fasta_records.register
def _(
    data: pathlib.Path,
    converter: OptConverterType = None,
    label_to_name: RenamerType = str,
):
    with open_(data, mode="rb") as infile:
        data: bytes = infile.read()
    return iter_fasta_records(data, converter=converter, label_to_name=label_to_name)


@iter_fasta_records.register
def _(
    data: io.TextIOWrapper,
    converter: OptConverterType = None,
    label_to_name: RenamerType = str,
):
    read_data: bytes = data.read().encode("utf8")
    return iter_fasta_records(
        read_data, converter=converter, label_to_name=label_to_name
    )


@iter_fasta_records.register
def _(data: list, converter: OptConverterType = None, label_to_name: RenamerType = str):
    return MinimalFastaParser(data, strict=False, label_to_name=label_to_name)
