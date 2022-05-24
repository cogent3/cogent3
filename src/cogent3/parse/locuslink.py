"""Parsers for the LL_tmpl file from LocusLink.

Notes:

The LocusLink format is documented in the README file, but unfortunately this
documentation is mostly lies. Fields that are supposed to be unique are
repeated, fields whose only association with each other is their order are
found out of order, etc.

I suspect that it is impossible to parse the entire file as it was intended,
and writing a parser that conforms to the specification is not useful because
the file does not match the specificiation. Consequently, I chose to break
the assocation between fields that are supposed to form 'sets' within a
subrecord rather than trying to figure out what the sets are from incomplete
data. This means that e.g. products will not be associated with _particular_
RNAs: however, all RNAs and all products produced by a locus will be returned.

The following fields are assumed to be unique (* = required):
    *LOCUSID
     CURRENT_LOCUSID
     LOCUS_CONFIRMED
     LOCUS_TYPE
    *ORGANISM
     STATUS
     OFFICIAL_SYMBOL
     PREFERRED_SYMBOL
     OFFICIAL_GENE_NAME
     PREFERRED_GENE_NAME

All other fields are assumed to be multiple, so will return a 1-item list if
they have a single record rather than returning a single item.

All records will be parsed if possible, typically as MappedRecord objects.
This applies especially to lines with pipe-delimited fields such as GO, CDD,
CONTIG, etc.

It is _likely_, but not necessarily true, that items at corresponding indices
in the lists for grouped fields (e.g. MAP and MAPLINK, PHENOTYPE and
PHENOTYPE_ID,BUTTON and LINK) refer to the same item (i.e. MAP[0] and MAPLINK[0]
are probably a map and its corresponding link).
"""
from cogent3.parse.record import (
    DelimitedSplitter,
    FieldWrapper,
    LineOrientedConstructor,
    MappedRecord,
    int_setter,
    list_adder,
    list_extender,
)
from cogent3.parse.record_finder import LabeledRecordFinder


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


def ll_start(line):
    """Returns True if line looks like the start of a LocusLink record."""
    return line.startswith(">>")


LLFinder = LabeledRecordFinder(ll_start)

pipes = DelimitedSplitter("|", None)
first_pipe = DelimitedSplitter("|")
commas = DelimitedSplitter(",", None)
first_colon = DelimitedSplitter(":", 1)

accession_wrapper = FieldWrapper(["Accession", "Gi", "Strain"], pipes)


def _read_accession(line):
    """Reads accession lines: format is Accession | Gi | Strain."""
    return MappedRecord(accession_wrapper(line))


rell_wrapper = FieldWrapper(["Description", "Id", "IdType", "Printable"], pipes)


def _read_rell(line):
    """Reads RELL lines: format is Description|Id|IdType|Printable"""
    return MappedRecord(rell_wrapper(line))


accnum_wrapper = FieldWrapper(["Accession", "Gi", "Strain", "Start", "End"], pipes)


def _read_accnum(line):
    """Reads ACCNUM lines: format is Accession|Gi|Strain|Start|End."""
    return MappedRecord(accnum_wrapper(line))


map_wrapper = FieldWrapper(["Location", "Source", "Type"], pipes)


def _read_map(line):
    """Reads MAP lines: format is Location|Source|Type."""
    return MappedRecord(map_wrapper(line))


sts_wrapper = FieldWrapper(
    ["Name", "Chromosome", "StsId", "Segment", "SequenceKnown", "Evidence"], pipes
)


def _read_sts(line):
    """Reads STS lines: format is in the full docstring.

    Format:
    Name|Chromosome|StsId|Segment|SequenceKnown|Evidence
    """
    return MappedRecord(sts_wrapper(line))


cdd_wrapper = FieldWrapper(["Name", "Key", "Score", "EValue", "BitScore"], pipes)


def _read_cdd(line):
    """Reads CDD lines: format is Name|Key|Score|EValue|BitScore."""
    return MappedRecord(cdd_wrapper(line))


comp_wrapper = FieldWrapper(
    [
        "TaxonId",
        "Symbol",
        "Chromosome",
        "Position",
        "LocusId",
        "ChromosomeSelf",
        "SymbolSelf",
        "MapName",
    ],
    pipes,
)


def _read_comp(line):
    """Reads COMP lines: format is in the full docstring.

    TaxonId|Symbol|Chromosome|Position|LocusId|ChromosomeSelf|SymbolSelf|MapName
    """
    return MappedRecord(comp_wrapper(line))


grif_wrapper = FieldWrapper(["PubMedId", "Description"], first_pipe)


def _read_grif(line):
    """Reads GRIF lines: format is PubMedId|Description."""
    return MappedRecord(grif_wrapper(line))


def _read_pmid(line):
    """Reads PMID lines: format is comma-delimited list of pubmed IDs."""
    return commas(line)


go_wrapper = FieldWrapper(
    ["Category", "Term", "EvidenceCode", "GoId", "Source", "PubMedId"], pipes
)


def _read_go(line):
    """Reads GO lines. Format: Category|Term|EvidenceCode|GoId|Source|PubMedId"""
    return MappedRecord(go_wrapper(line))


extannot_wrapper = FieldWrapper(
    ["Category", "Term", "EvidenceCode", "Source", "PubMedId"], pipes
)


def _read_extannot(line):
    """Reads EXTANNOT lines. format: Category|Term|EvidenceCode|Source|PubMedId"""
    return MappedRecord(extannot_wrapper(line))


contig_wrapper = FieldWrapper(
    [
        "Accession",
        "Gi",
        "Strain",
        "From",
        "To",
        "Orientation",
        "Chromosome",
        "Assembly",
    ],
    pipes,
)


def _read_contig(line):
    """Reads CONTIG lines. Format described in full docstring.

    Accession|Gi|Strain|From|To|Orientation|Chromosome|Assembly
    """
    return MappedRecord(contig_wrapper(line))


_ll_multi = dict.fromkeys(
    "RELL NG NR NM NC NP PRODUCT TRANSVAR ASSEMBLY CONTIG XG XR EVID XM XP CDD ACCNUM TYPE PROT PREFERRED_PRODUCT ALIAS_SYMBOL ALIAS_PROT PHENOTYPE PHENOTYPE_ID SUMMARY UNIGENE OMIM CHR MAP MAPLINK STS COMP ECNUM BUTTON LINK DB_DESCR DB_LINK PMID GRIF SUMFUNC GO EXTANNOT".split()
)
for i in list(_ll_multi.keys()):
    _ll_multi[i] = []


class LocusLink(MappedRecord):
    """Holds data for a LocusLink record."""

    Required = _ll_multi
    Aliases = {
        "LOCUSID": "LocusLinkId",
        "CURRENT_LOCUS_ID": "CurrentLocusId",
        "LOCUS_CONFIRMED": "LocusConfirmed",
        "LOCUS_TYPE": "LocusType",
        "ORGANISM": "Species",
        "RELL": "RelatedLoci",
        "STATUS": "LocusStatus",
        "PRODUCT": "Products",
        "TRANSVAR": "TranscriptionVariants",
        "ASSEMBLY": "Assemblies",
        "CONTIG": "Contigs",
        "EVID": "ContigEvidenceCodes",
        "ACCNUM": "AccessionNumbers",
        "TYPE": "AccessionTypes",
        "PROT": "ProteinIds",
        "OFFICIAL_SYMBOL": "OfficialSymbol",
        "PREFERRED_SYMBOL": "PreferredSymbol",
        "OFFICIAL_GENE_NAME": "OfficialGeneName",
        "PREFERRED_GENE_NAME": "PreferredGeneName",
        "PREFERRED_PRODUCT": "PreferredProducts",
        "ALIAS_SYMBOL": "SymbolAliases",
        "ALIAS_PROT": "ProteinAliases",
        "PHENOTYPE": "Phenotypes",
        "PHENOTYPE_ID": "PhenotypeIds",
        "SUMMARY": "Summaries",
        "UNIGENE": "UnigeneIds",
        "OMIM": "OmimIds",
        "CHR": "Chromosomes",
        "MAP": "Maps",
        "MAPLINK": "MapLinks",
        "STS": "Sts",
        "COMP": "ComparativeMapLinks",
        "ECNUM": "EcIds",
        "BUTTON": "Buttons",
        "LINK": "Links",
        "DB_DESCR": "DbDescriptions",
        "DB_LINK": "DbLinks",
        "PMID": "PubMedIds",
        "GRIF": "Grifs",
        "SUMFUNC": "FunctionSummaries",
        "GO": "GoIds",
        "EXTANNOT": "ExternalAnnotations",
    }


def _accession_adder(obj, field, line):
    """Adds accessions to relevant field"""
    list_adder(obj, field, _read_accession(line))


def _accnum_adder(obj, field, line):
    """Adds accnum to relevant field"""
    list_adder(obj, field, _read_accnum(line))


def _rell_adder(obj, field, line):
    """Adds rell to relevant field"""
    list_adder(obj, field, _read_rell(line))


def _map_adder(obj, field, line):
    """Adds map to relevant field"""
    list_adder(obj, field, _read_map(line))


def _sts_adder(obj, field, line):
    """Adds sts to relevant field"""
    list_adder(obj, field, _read_sts(line))


def _cdd_adder(obj, field, line):
    """Adds cdd to relevant field"""
    list_adder(obj, field, _read_cdd(line))


def _comp_adder(obj, field, line):
    """Adds comp to relevant field"""
    list_adder(obj, field, _read_comp(line))


def _grif_adder(obj, field, line):
    """Adds grif to relevant field"""
    list_adder(obj, field, _read_grif(line))


def _pmid_adder(obj, field, line):
    """Adds pmid to relevant field"""
    list_extender(obj, field, _read_pmid(line))


def _assembly_adder(obj, field, line):
    """Adds assembly to relevant field"""
    list_adder(obj, field, commas(line))


def _go_adder(obj, field, line):
    """Adds go to relevant field"""
    list_adder(obj, field, _read_go(line))


def _extannot_adder(obj, field, line):
    """Adds pmid to relevant field"""
    list_adder(obj, field, _read_extannot(line))


def _generic_adder(obj, field, line):
    """Adds line to relevant field, unparsed"""
    list_adder(obj, field, line.strip())


def _contig_adder(obj, field, line):
    """Adds contig to relevant field"""
    list_adder(obj, field, _read_contig(line))


_ll_fieldmap = {}
for field in ["LOCUSID", "CURRENT_LOCUSID"]:
    _ll_fieldmap[field] = int_setter
_ll_fieldmap["RELL"] = _rell_adder
_ll_fieldmap["MAP"] = _map_adder
_ll_fieldmap["STS"] = _sts_adder
_ll_fieldmap["COMP"] = _comp_adder
_ll_fieldmap["GRIF"] = _grif_adder
_ll_fieldmap["PMID"] = _pmid_adder
_ll_fieldmap["GO"] = _go_adder
_ll_fieldmap["EXTANNOT"] = _extannot_adder
_ll_fieldmap["MAP"] = _map_adder
_ll_fieldmap["CDD"] = _cdd_adder
_ll_fieldmap["ASSEMBLY"] = _assembly_adder
_ll_fieldmap["CONTIG"] = _contig_adder
for field in "NG ACCNUM".split():
    _ll_fieldmap[field] = _accnum_adder
for field in "NR NM NC NP XG XR XM XP PROT".split():
    _ll_fieldmap[field] = _accession_adder
for field in _ll_multi:
    if field not in _ll_fieldmap:
        _ll_fieldmap[field] = _generic_adder


LinesToLocusLink = LineOrientedConstructor()
LinesToLocusLink.Constructor = LocusLink
LinesToLocusLink.FieldMap = _ll_fieldmap
LinesToLocusLink.LabelSplitter = first_colon


def LocusLinkParser(lines):
    """Treats lines as a stream of unigene records"""
    for record in LLFinder(lines):
        curr = LinesToLocusLink(record)
        yield curr
