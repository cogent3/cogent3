"""Parsers for blast, psi-blast and blat.
"""
from cogent3.parse.record_finder import (
    DelimitedRecordFinder,
    LabeledRecordFinder,
    never_ignore,
)


__author__ = "Micah Hamady"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Micah Hamady", "Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Micah Hamady"
__email__ = "hamady@colorado.edu"
__status__ = "Prototype"

strip = str.strip
upper = str.upper


def iter_finder(line):
    """Split record on rows that start with iteration label."""
    return line.startswith("# Iteration:")


def query_finder(line):
    """Split record on rows that start with query label."""
    return line.startswith("# Query:")


def iteration_set_finder(line):
    """Split record on rows that begin a new iteration."""
    return line.startswith("# Iteration: 1")


def _is_junk(line, t_strs):
    """Ignore empty line, line with blast info, or whitespace line"""
    # empty or white space
    if not line or not line.strip():
        return True
    # blast info line
    for t_str in t_strs:
        if line.startswith(f"# {t_str}"):
            return True
    return False


def is_blast_junk(line):
    """Ignore empty line or lines with blast info"""
    return _is_junk(line, ("BLAST", "TBLAS"))


def is_blat_junk(line):
    """Ignore empty line or lines with blat info"""
    return _is_junk(line, ("BLAT",))


label_constructors = {"ITERATION": int}  # add other label constructors here


def make_label(line):
    """Make key, value for colon-delimited comment lines.

    WARNING: Only maps the data type if the key is in label_constructors above.
    """
    if not line.startswith("#"):
        raise ValueError("Labels must start with a # symbol.")

    if line.find(":") == -1:
        raise ValueError("Labels must contain a : symbol.")

    key, value = list(map(strip, line[1:].split(":", 1)))
    key = key.upper()
    if key in label_constructors:
        value = label_constructors[key](value)
    return key, value


BlatFinder = LabeledRecordFinder(query_finder, constructor=strip, ignore=is_blat_junk)

BlastFinder = LabeledRecordFinder(query_finder, constructor=strip, ignore=is_blast_junk)

PsiBlastFinder = LabeledRecordFinder(
    iter_finder, constructor=strip, ignore=is_blast_junk
)

PsiBlastQueryFinder = LabeledRecordFinder(
    iteration_set_finder, constructor=strip, ignore=is_blast_junk
)


def GenericBlastParser9(lines, finder, make_col_headers=False):
    """Yields successive records from lines (props, data list)

    Infile must in blast9 format

    finder: labeled record finder function

    make_col_header: adds column headers (from fields entry) as first
    row in data output

    props is a dict of {UPPERCASE_KEY:value}.
    data_list is a list of list of strings, optionally with header first.
    """
    for rec in finder(lines):
        props = {}
        data = []
        for line in rec:
            if line.startswith("#"):
                label, value = make_label(line)
                props[label] = value

                # check if need to insert column headers
                if make_col_headers and label == "FIELDS":
                    data.insert(0, list(map(upper, list(map(strip, value.split(","))))))

            else:
                data.append(list(map(strip, line.split("\t"))))
        yield props, data


def TableToValues(table, constructors=None, header=None):
    """Converts table to values according to constructors.

    Returns (table, header).
    Use dict([(val, i) for i, val in enumerate(header)]) to get back
    a dict mapping the fields to indices in each row.
    """
    if header is None:  # assume first row of table
        header = table[0]
        table = table[1:]
    c_list = [constructors.get(k, str) for k in header]
    return [[c(val) for c, val in zip(c_list, row)] for row in table], header


psiblast_constructors = {
    "% identity": float,
    "alignment length": int,
    "mismatches": int,
    "gap openings": int,
    "q. start": int,
    "q. end": int,
    "s. start": int,
    "s. end": int,
    "e-value": float,
    "bit score": float,
}
# make case-insensitive
for key, val in list(psiblast_constructors.items()):
    psiblast_constructors[key.upper()] = val


def PsiBlastTableParser(table):
    return TableToValues(table, psiblast_constructors)


def MinimalBlastParser9(lines, include_column_names=False):
    """Yields succesive records from lines (props, data list).

    lines must be BLAST output format.
    """
    return GenericBlastParser9(lines, BlastFinder, include_column_names)


def MinimalPsiBlastParser9(lines, include_column_names=False):
    """Yields successive records from lines (props, data list)

    lines must be of psi-blast output format
    """
    return GenericBlastParser9(lines, PsiBlastFinder, include_column_names)


def MinimalBlatParser9(lines, include_column_names=True):
    """Yields successive records from lines (props, data list)

    lines must be of blat output (blast9) format
    """
    return GenericBlastParser9(lines, BlatFinder, include_column_names)


def PsiBlastParser9(lines):
    """Returns fully parsed PSI-BLAST result.

    result['query'] gives all the results for specified query sequence.
    result['query'][i] gives result for iteration i (offset by 1: zero-based)
    if x = result['query']['iteration']:
        x[0]['e-value'] gives the e-value of the first result.

    WARNING: designed for ease of use, not efficiency!"""
    result = {}
    for query in PsiBlastQueryFinder(lines):
        first_query = True  # if it's the first, need to make the entry
        for properties, record in MinimalPsiBlastParser9(query, True):
            if first_query:
                curr_resultset = []
                result[properties["QUERY"].split()[0]] = curr_resultset
                first_query = False
            table, header = PsiBlastTableParser(record)
            curr_resultset.append([dict(list(zip(header, row))) for row in table])
    return result


def get_blast_ids(props, data, filter_identity, threshold, keep_values):
    """
    Extract ids from blast output
    """
    fields = list(map(strip, props["FIELDS"].upper().split(",")))

    # get column index of protein ids we want
    p_ix = fields.index("SUBJECT ID")
    # get column index to screen by
    if filter_identity:
        e_ix = fields.index("% IDENTITY")
    else:
        e_ix = fields.index("E-VALUE")
    # no filter, returh all
    if not threshold:
        if keep_values:
            return [(x[p_ix], x[e_ix]) for x in data]
        else:
            return [x[p_ix] for x in data]
    else:
        # will raise exception if invalid threshold passed
        max_val = float(threshold)

        # figure out what we're keeping
        def ok_val(val):
            if threshold:
                return val <= max_val
            return val >= max_val

        if keep_values:
            return [(x[p_ix], x[e_ix]) for x in data if ok_val(float(x[e_ix]))]
        else:
            return [x[p_ix] for x in data if ok_val(float(x[e_ix]))]


def AllProteinIds9(
    lines,
    filter_identity=True,
    threshold=None,
    keep_below_threshold=True,
    output_parser=MinimalPsiBlastParser9,
    keep_values=False,
):
    """Helper to extract just protein ids from each blast search

    lines: output file in output format #9.
    filter_identity: when True, use % identity to filter, else use e-value
    threshold: when None, all results are returned. When not None, used
        as a threshold to filter results.
    keep_below_threshold: when True, keeps any rows below given threshold, else
        keep any rows above threshold
    output_parser: minimal output parser to use (e.g. minimalpsiblast)
    keep_values: if True, returns tuples of (id, value) rather than just ids.

    Note that you can feed it successive output from PsiBlastQueryFinder if
    you have a PSI-BLAST file with multiple input queries.

    Subject ids are stable relative to original order.
    """

    mpbp = output_parser(lines)

    # get last record.
    props = data = None
    out_ids = {}
    out_ct = 1
    for rec in mpbp:
        props, data = rec
        out_ids[out_ct] = get_blast_ids(
            props, data, filter_identity, threshold, keep_values
        )
        out_ct += 1
    return out_ids


def LastProteinIds9(
    lines,
    filter_identity=True,
    threshold=None,
    keep_below_threshold=True,
    output_parser=MinimalPsiBlastParser9,
    keep_values=False,
):
    """Helper to extract just protein ids from last psi-blast iteration.

    lines: output file in output format #9.
    filter_identity: when True, use % identity to filter, else use e-value
    threshold: when None, all results are returned. When not None, used
        as a threshold to filter results.
    keep_below_threshold: when True, keeps any rows below given threshold, else
        keep any rows above threshold
    output_parser: minimal output parser to use (e.g. minimalpsiblast)
    keep_values: if True, returns tuples of (id, value) rather than just ids.

    Note that you can feed it successive output from PsiBlastQueryFinder if
    you have a PSI-BLAST file with multiple input queries.

    Subject ids are stable relative to original order.
    """

    mpbp = output_parser(lines)
    # get last record.
    props = data = None
    for rec in mpbp:
        props, data = rec
    if not (props and data):
        return []
    return get_blast_ids(props, data, filter_identity, threshold, keep_values)


def QMEBlast9(lines):
    """Returns query, match and e-value for each line in Blast-9 output.

    WARNING: Allows duplicates in result.

    WARNING: If you use this on PSI-BLAST output, will not check that you're
    only getting stuff from the last iteration but will give you everything.
    The advantage is that you keep stuff that drops out of the profile. The
    disadvantage is that you keep stuff that drops out of the profile...
    """
    result = []
    for line in lines:
        if line.startswith("#"):
            continue
        try:
            fields = line.split("\t")
            result.append((fields[0], fields[1], float(fields[-2])))
        except (TypeError, ValueError, IndexError):
            pass
    return result


def QMEPsiBlast9(lines):
    """Returns successive query, match, e-value from lines of Psi-Blast run.

    Assumes tabular output. Uses last iteration from each query.

    WARNING: Allows duplicates in result
    """
    result = []
    for query in PsiBlastQueryFinder(lines):
        for iteration in PsiBlastFinder(query):
            pass
        result.extend(QMEBlast9(iteration))
    return result


fastacmd_taxonomy_splitter = DelimitedRecordFinder(delimiter="", ignore=never_ignore)
fasta_field_map = {
    "NCBI sequence id": "seq_id",
    "NCBI taxonomy id": "tax_id",
    "Common name": "common_name",
    "Scientific name": "scientific_name",
}


def FastacmdTaxonomyParser(lines):
    """Yields successive records from the results of fastacmd -T.

    Format is four lines separated by newline:
    NCBI sequence
    NCBI taxonomy
    Common name
    Scientific name

    Result is dict with keys by seq_id, tax_id, common_name, scientific_name.
    """
    for group in fastacmd_taxonomy_splitter(lines):
        result = {}
        for line in group:
            try:
                header, data = line.split(":", 1)
                result[fasta_field_map[header]] = data.strip()
            except (TypeError, ValueError, KeyError):
                continue
        yield result
