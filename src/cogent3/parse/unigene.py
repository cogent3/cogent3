#!/usr/bin/env python
"""Parsers for the various files in the UniGene database.
"""
from cogent3.parse.record import (
    ByPairs,
    LineOrientedConstructor,
    MappedRecord,
    equal_pairs,
    int_setter,
    list_adder,
    semi_splitter,
)
from cogent3.parse.record_finder import GbFinder


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

maketrans = str.maketrans
strip = str.strip
rstrip = str.rstrip


def _read_sts(line):
    """Turns an STS line (without label) into a record.

    Infuritatingly, STS lines are not semicolon-delimited, and spaces appear
    in places they shouldn't. This was the case as of 10/9/03: expect this
    'feature' to be unstable!
    """
    filtered = line.replace("=", " ")
    return MappedRecord(list(ByPairs(filtered.split())))


def _read_expression(line):
    """Turns a semicolon-delimited  expression line into list of expressions"""
    return semi_splitter(line)


class UniGeneSeqRecord(MappedRecord):
    Aliases = {
        "ACC": "Accession",
        "CLONE": "CloneId",
        "END": "End",
        "LID": "LibraryId",
        "SEQTYPE": "SequenceType",
        "TRACE": "Trace",
        "EST": "EstId",
        "NID": "NucleotideId",
        "PID": "ProteinId",
    }


class UniGeneProtSimRecord(MappedRecord):
    Aliases = {
        "ORG": "Species",
        "PROTGI": "ProteinGi",
        "ProtId": "ProteinId",
        "PCT": "PercentSimilarity",
        "ALN": "AlignmentScore",
    }


def _read_seq(line):
    """Turns a sequence line into a UniGeneSeqRecord.

    BEWARE: first level delimiter is ';' and second level delimiter is '=', but
    '=' can also appear inside the _value_ of the second level!
    """
    first_level = semi_splitter(line)
    second_level = list(map(equal_pairs, first_level))
    return UniGeneSeqRecord(second_level)


def _read_protsim(line):
    """Turns a protsim line into a UniGeneProtSim record.

    BEWARE: first level delimiter is ';' and second level delimiter is '=', but
    '=' can also appear inside the _value_ of the second level!
    """
    first_level = semi_splitter(line)
    second_level = list(map(equal_pairs, first_level))
    return UniGeneProtSimRecord(second_level)


class UniGene(MappedRecord):
    """Holds data for a UniGene record."""

    Required = {"STS": [], "PROTSIM": [], "SEQUENCE": [], "EXPRESS": []}
    Aliases = {
        "STS": "Sts",
        "PROTSIM": "ProteinSimilarities",
        "SEQUENCE": "SequenceIds",
        "SCOUNT": "SequenceCount",
        "CTYOBAND": "CytoBand",
        "EXPRESS": "ExpressedIn",
        "CHROMOSOME": "Chromosome",
        "ID": "UniGeneId",
        "TITLE": "UniGeneTitle",
        "LOCUSLINK": "LocusLinkId",
    }


def _expressions_setter(obj, field, val):
    """Sets specified field to a list of expressions"""
    setattr(obj, field, semi_splitter(val))


def _sts_adder(obj, field, val):
    """Appends the current STS-type record to specified field"""
    list_adder(obj, field, _read_sts(val))


def _seq_adder(obj, field, val):
    """Appends the current Sequence-type record to specified field"""
    list_adder(obj, field, _read_seq(val))


def _protsim_adder(obj, field, val):
    """Appends the current ProtSim record to specified field"""
    list_adder(obj, field, _read_protsim(val))


LinesToUniGene = LineOrientedConstructor()
LinesToUniGene.Constructor = UniGene
LinesToUniGene.FieldMap = {
    "LOCUSLINK": int_setter,
    "EXPRESS": _expressions_setter,
    "PROTSIM": _protsim_adder,
    "SCOUNT": int_setter,
    "SEQUENCE": _seq_adder,
    "STS": _sts_adder,
}


def UniGeneParser(lines):
    """Treats lines as a stream of unigene records"""
    for record in GbFinder(lines):
        curr = LinesToUniGene(record)
        del curr["//"]  # clean up delimiter
        yield curr


if __name__ == "__main__":
    from sys import argv, stdout

    filename = argv[1]
    count = 0
    for record in UniGeneParser(open(filename)):
        stdout.write(".")
        stdout.flush()
        count += 1
    print(f"read {count} records")
