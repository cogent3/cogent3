import pathlib

import numpy
import pytest

from cogent3.core import alphabet as c3_alphabet
from cogent3.parse.fastq import iter_fastq_records
from cogent3.parse.record import RecordError

EXPECTED_LABELS = [
    "SIM:1:FCX:1:15:6329:1045 1:N:0:ATCCGA",
    "SIM:1:FCX:1:15:6329:1046 1:N:0:ATCCGA",
    "SIM:1:FCX:1:15:6329:1047 1:N:0:ATCCGA",
    "SIM:1:FCX:1:15:6329:1048 1:N:0:ATCCGA",
    "SIM:1:FCX:1:15:6329:1049 1:N:0:ATCCGA",
]
EXPECTED_SEQS = [
    "GCTCAGCATCGAGAAGCTTAGCAACTTGGCAACGT",
    "GAAATATCTGCAAGCCATGTGTGCGTTCCCATATC",
    "AAGTCTCAGGCATACTTGCCTGGCCACAGCAAGAG",
    "CCAGCATCACAACAGTCTGTAGTTCAGGTCTTACG",
    "TCTTGCATAGCCCTGGCATCCGAACTCAGAGCCCG",
]
EXPECTED_QUALS = [
    "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    "9CFFFFFFHHHHHJJJJJJJJIJJJJJJJJJHHHH",
    "CCCFFFFFHHHHHJJJJJJJJJJJJJJIJJJJJJ?",
    "@@CFFDDFGHHHHJJIIJJJJJJJJJJJJJJHJII",
    "?@@DDDDDHHHHFHIIIJJIJJJJIIIIGIJJJJJ",
]


@pytest.mark.parametrize("coerce", [str, lambda p: p])
def test_iter_fastq_records_path(DATA_DIR, coerce):
    records = list(iter_fastq_records(coerce(DATA_DIR / "fastq.txt")))
    assert len(records) == 5
    labels, seqs, quals = zip(*records)
    assert list(labels) == EXPECTED_LABELS
    assert list(seqs) == EXPECTED_SEQS
    assert list(quals) == EXPECTED_QUALS


def test_iter_fastq_records_qual_converter(DATA_DIR):
    path = DATA_DIR / "fastq.txt"
    qual_conv = c3_alphabet.make_qual_converter(c3_alphabet.PhredEncoding.PHRED_33)
    records = list(iter_fastq_records(path, qual_converter=qual_conv))
    label, seq, qual = records[0]
    assert label == EXPECTED_LABELS[0]
    assert seq == EXPECTED_SEQS[0]
    assert isinstance(qual, numpy.ndarray)
    # 'I' is ASCII 73 which under Phred+33 is quality 40
    assert (qual == numpy.full(35, 40, dtype=numpy.uint8)).all()


def test_iter_fastq_records_seq_converter(DATA_DIR):
    path = DATA_DIR / "fastq.txt"
    records = list(iter_fastq_records(path, converter=lambda b: b))
    label, seq, qual = records[0]
    assert seq == EXPECTED_SEQS[0].encode("utf8")
    assert label == EXPECTED_LABELS[0]
    assert qual == EXPECTED_QUALS[0]


def test_iter_fastq_records_explicit_none_converters(DATA_DIR):
    path = DATA_DIR / "fastq.txt"
    records = list(iter_fastq_records(path, converter=None, qual_converter=None))
    assert records[0][1] == EXPECTED_SEQS[0]
    assert records[0][2] == EXPECTED_QUALS[0]


@pytest.mark.parametrize("bad", [42, b"@id\nACGT\n+\nIIII\n"])
def test_iter_fastq_records_unsupported_type(bad):
    with pytest.raises(TypeError):
        list(iter_fastq_records(bad))


def _write_fastq(tmp_path: pathlib.Path, contents: str) -> pathlib.Path:
    path = tmp_path / "input.fastq"
    path.write_text(contents)
    return path


def test_iter_fastq_records_short_record_raises(tmp_path):
    # final record has only 2 lines
    path = _write_fastq(tmp_path, "@a\nACGT\n+\nIIII\n@b\nACGT\n")
    with pytest.raises(RecordError):
        list(iter_fastq_records(path))


def test_iter_fastq_records_bad_label_raises(tmp_path):
    path = _write_fastq(tmp_path, "id\nACGT\n+\nIIII\n")
    with pytest.raises(RecordError):
        list(iter_fastq_records(path))


def test_iter_fastq_records_bad_separator_raises(tmp_path):
    path = _write_fastq(tmp_path, "@id\nACGT\nbad\nIIII\n")
    with pytest.raises(RecordError):
        list(iter_fastq_records(path))


def test_iter_fastq_records_empty_file(tmp_path):
    path = _write_fastq(tmp_path, "")
    assert list(iter_fastq_records(path)) == []
