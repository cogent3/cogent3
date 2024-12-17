import pytest

import cogent3
from cogent3.format.fasta import seqs_to_fasta


@pytest.fixture
def base_data():
    """Basic test data setup."""
    return {
        "strings": ["AAAA", "CCCC", "gggg", "uuuu"],
        "labels": ["1st", "2nd", "3rd", "4th"],
    }


@pytest.fixture
def ascii_sequences(base_data):
    """Create ASCII sequences with labels and names."""
    ASCII = cogent3.get_moltype("text")
    sequences = {
        "with_labels": [ASCII.make_seq(seq=seq) for seq in base_data["strings"]],
        "with_names": [ASCII.make_seq(seq=seq) for seq in base_data["strings"]],
    }

    for label, seq in zip(base_data["labels"], sequences["with_names"], strict=False):
        seq.name = label

    return sequences


@pytest.fixture
def fasta_formats():
    """Different FASTA format strings."""
    return {
        "no_label": ">0\nAAAA\n>1\nCCCC\n>2\ngggg\n>3\nuuuu\n",
        "with_label": ">1st\nAAAA\n>2nd\nCCCC\n>3rd\nGGGG\n>4th\nUUUU\n",
        "with_label_lw2": ">1st\nAA\nAA\n>2nd\nCC\nCC\n>3rd\nGG\nGG\n>4th\nUU\nUU\n",
    }


@pytest.fixture
def alignment_dict():
    """Create alignment dictionary."""
    return {
        "1st": "AAAA",
        "2nd": "CCCC",
        "3rd": "GGGG",
        "4th": "UUUU",
    }


@pytest.fixture
def alignment_object(alignment_dict):
    """Create alignment object."""
    aln = cogent3.make_aligned_seqs(
        alignment_dict,
        moltype="text",
    )
    aln.RowOrder = ["1st", "2nd", "3rd", "4th"]
    return aln


def test_empty_alignment_to_fasta():
    """Test empty alignment conversion to FASTA."""
    assert seqs_to_fasta({}) == ""


def test_alignment_to_fasta(alignment_dict, fasta_formats):
    """Test alignment conversion to FASTA."""
    assert seqs_to_fasta(alignment_dict) == fasta_formats["with_label"]


def test_alignment_to_fasta_with_block_size(alignment_dict, fasta_formats):
    """Test alignment conversion to FASTA with specific block size."""
    assert (
        seqs_to_fasta(alignment_dict, block_size=2) == fasta_formats["with_label_lw2"]
    )
