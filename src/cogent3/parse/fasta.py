"""Parsers for FASTA and related formats.
"""
import re

from collections.abc import Callable

import cogent3

from cogent3.core.info import Info
from cogent3.core.moltype import ASCII, BYTES
from cogent3.parse.record import RecordError
from cogent3.parse.record_finder import LabeledRecordFinder
from cogent3.util.io import open_


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"


strip = str.strip

Sequence = BYTES.make_seq


def is_fasta_label(x):
    """Checks if x looks like a FASTA label line."""
    return x.startswith(">")


def is_gde_label(x):
    """Checks if x looks like a GDE label line."""
    return x and x[0] in "%#"


def is_blank_or_comment(x):
    """Checks if x is blank or a FASTA comment line."""
    return (not x) or x.startswith("#") or x.isspace()


def is_blank(x):
    """Checks if x is blank."""
    return (not x) or x.isspace()


FastaFinder = LabeledRecordFinder(is_fasta_label, ignore=is_blank_or_comment)


def MinimalFastaParser(
    infile, strict=True, label_to_name=str, finder=FastaFinder, label_characters=">"
):
    """Yields successive sequences from infile as (label, seq) tuples.

    If strict is True (default), raises RecordError when label or seq missing.
    """
    try:
        infile = open_(infile)
        close_at_end = True
    except (TypeError, AttributeError):
        close_at_end = False

    for rec in finder(infile):
        # first line must be a label line
        if not rec[0][0] in label_characters:
            if strict:
                raise RecordError(f"Found Fasta record without label line: {rec}")
            continue
        # record must have at least one sequence
        if len(rec) < 2:
            if strict:
                raise RecordError(f"Found label line without sequences: {rec}")
            else:
                continue

        label = rec[0][1:].strip()
        label = label_to_name(label)
        seq = "".join(rec[1:])

        yield label, seq

    if close_at_end:
        infile.close()


GdeFinder = LabeledRecordFinder(is_gde_label, ignore=is_blank)


def MinimalGdeParser(infile, strict=True, label_to_name=str):
    return MinimalFastaParser(
        infile, strict, label_to_name, finder=GdeFinder, label_characters="%#"
    )


def xmfa_label_to_name(line):
    (loc, strand, contig) = line.split()
    (sp, loc) = loc.split(":")
    (lo, hi) = [int(x) for x in loc.split("-")]
    if strand == "-":
        (lo, hi) = (hi, lo)
    else:
        assert strand == "+"
    return f"{sp}:{contig}:{lo}-{hi}"


def is_xmfa_blank_or_comment(x):
    """Checks if x is blank or an XMFA comment line."""
    return (not x) or x.startswith("=") or x.isspace()


XmfaFinder = LabeledRecordFinder(is_fasta_label, ignore=is_xmfa_blank_or_comment)


def MinimalXmfaParser(infile, strict=True):
    # Fasta-like but with header info like ">1:10-1000 + chr1"
    return MinimalFastaParser(
        infile, strict, label_to_name=xmfa_label_to_name, finder=XmfaFinder
    )


def MinimalInfo(label):
    """Minimal info data maker: returns name, and empty dict for info{}."""
    return label, {}


def NameLabelInfo(label):
    """Returns name as label split on whitespace, and label in Info."""
    return label.split()[0], {"label": label}


def FastaParser(infile, seq_maker=None, info_maker=MinimalInfo, strict=True):
    """Yields successive sequences from infile as (name, sequence) tuples.

    Constructs the sequence using seq_maker(seq, info=Info(info_maker(label))).

    If strict is True (default), raises RecordError when label or seq missing.
    Also raises RecordError if seq_maker fails.

    It is info_maker's responsibility to raise the appropriate RecordError or
    FieldError on failure.

    Result of info_maker need not actually be an info object, but can just be
    a dict or other data that Info can use in its constructor.
    """
    if seq_maker is None:
        seq_maker = Sequence
    for label, seq in MinimalFastaParser(infile, strict=strict):
        if strict:
            # need to do error checking when constructing info and sequence
            try:
                name, info = info_maker(label)  # will raise exception if bad
                yield name, seq_maker(seq, name=name, info=info)
            except Exception:
                raise RecordError(
                    f"Sequence construction failed on record with label {label}"
                )
        else:
            # not strict: just skip any record that raises an exception
            try:
                name, info = info_maker(label)
                yield (name, seq_maker(seq, name=name, info=info))
            except Exception:
                continue


# labeled fields in the NCBI FASTA records
NcbiLabels = {"dbj": "DDBJ", "emb": "EMBL", "gb": "GenBank", "ref": "RefSeq"}


def NcbiFastaLabelParser(line):
    """Creates an Info object and populates it with the line contents.

    As of 11/12/03, all records in genpept.fsa and the human RefSeq fasta
    files were consistent with this format.
    """
    info = Info()
    try:
        ignore, gi, db, db_ref, description = list(map(strip, line.split("|", 4)))
    except ValueError:  # probably got wrong value
        raise RecordError(f"Unable to parse label line {line}")
    info.GI = gi
    info[NcbiLabels[db]] = db_ref
    info.Description = description
    return gi, info


def NcbiFastaParser(infile, seq_maker=None, strict=True):
    return FastaParser(
        infile, seq_maker=seq_maker, info_maker=NcbiFastaLabelParser, strict=strict
    )


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
    for index, field, converter in field_formatters:
        if field in display_template:
            indexed = True
    assert indexed, f"display_template [{display_template}] does not use a field name"
    sep = re.compile(f"[{split_with}]")

    def call(label):
        label = [label, label[1:]][label[0] == ">"]
        label = sep.split(label)
        if DEBUG:
            print(label)
        info = Info()
        for index, name, converter in field_formatters:
            if isinstance(converter, Callable):
                try:
                    info[name] = converter(label[index])
                except IndexError:
                    raise IndexError(
                        "parsing label %s failed for property %s at index %s"
                        % (label, name, index)
                    )
            else:
                info[name] = label[index]
        return RichLabel(info, display_template)

    return call


def GroupFastaParser(
    data,
    label_to_name,
    group_key="Group",
    aligned=False,
    moltype=ASCII,
    done_groups=None,
    DEBUG=False,
):
    """yields related sequences as a separate seq collection

    Parameters
    ----------
    data
        line iterable data source
    label_to_name
        LabelParser callback
    group_key
        name of group key in RichLabel.info object
    aligned
        whether sequences are to be considered aligned
    moltype
        default is ASCII
    done_groups
        series of group keys to be excluded

    """

    done_groups = [[], done_groups][done_groups is not None]
    parser = MinimalFastaParser(data, label_to_name=label_to_name, finder=XmfaFinder)
    group_ids = []
    current_collection = {}
    for label, seq in parser:
        seq = moltype.make_seq(seq, name=label, info=label.info)
        if DEBUG:
            print("str(label) ", str(label), "repr(label)", repr(label))
        if not group_ids or label.info[group_key] in group_ids:
            current_collection[label] = seq
            if not group_ids:
                group_ids.append(label.info[group_key])
        else:
            # we finish off check of current before creating a collection
            if group_ids[-1] not in done_groups:
                info = Info(Group=group_ids[-1])
                if DEBUG:
                    print(
                        "GroupParser collection keys", list(current_collection.keys())
                    )
                seqs = cogent3.make_aligned_seqs(current_collection, moltype=moltype)
                seqs.info = info
                yield seqs
            current_collection = {label: seq}
            group_ids.append(label.info[group_key])
    info = Info(Group=group_ids[-1])
    func = cogent3.make_aligned_seqs if aligned else cogent3.make_unaligned_seqs
    seqs = func(current_collection, moltype=moltype, info=info)
    yield seqs
