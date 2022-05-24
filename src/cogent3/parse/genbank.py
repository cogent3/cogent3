#!/usr/bin/env python
from cogent3.core.annotation import Feature
from cogent3.core.genetic_code import GeneticCodes
from cogent3.core.info import Info
from cogent3.core.moltype import get_moltype
from cogent3.parse.record import FieldWrapper
from cogent3.parse.record_finder import (
    DelimitedRecordFinder,
    LabeledRecordFinder,
)


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Matthew Wakefield", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

maketrans = str.maketrans
strip = str.strip
rstrip = str.rstrip

all_chars = maketrans("", "")
dna_lc = "utacgrywsmkbdhvn"
dna_lc_cmp = "aatgcyrwskmvhdbn"
dna_trans = maketrans(dna_lc + dna_lc.upper(), dna_lc_cmp + dna_lc_cmp.upper())
rna_lc = "utacgrywsmkbdhvn"
rna_lc_cmp = "aaugcyrwskmvhdbn"
rna_trans = maketrans(rna_lc + rna_lc.upper(), rna_lc_cmp + rna_lc_cmp.upper())

locus_fields = [None, "locus", "length", None, "mol_type", "topology", "db", "date"]
_locus_parser = FieldWrapper(locus_fields)

# need to turn off line stripping, because whitespace is significant
GbFinder = DelimitedRecordFinder("//", constructor=rstrip)


class PartialRecordError(Exception):
    pass


def parse_locus(line):
    """Parses a locus line, including conversion of Length to an int.

    WARNING: Gives incorrect results on legacy records that omit the topology.
    All records spot-checked on 8/30/05 had been updated to include the topology
    even when prior versions omitted it.
    """
    result = _locus_parser(line)
    try:
        result["length"] = int(result["length"])
    except KeyError as e:
        raise PartialRecordError(e)

    if None in result:
        del result[None]
    return result


def parse_single_line(line):
    """Generic parser: splits off the label, and return the rest."""
    label, data = line.split(None, 1)
    return data.rstrip()


def indent_splitter(lines):
    """Yields the lines whenever it hits a line with same indent level as first."""
    first_line = True
    curr = []
    for line in lines:
        # skip blank lines
        line = line.rstrip()
        if not line:
            continue
        # need to figure out indent if first line
        if first_line:
            indent = len(line) - len(line.lstrip())
            curr.append(line)
            first_line = False
        elif len(line) > indent and line[indent].isspace():
            curr.append(line)
        else:  # got a line that doesn't match the indent
            yield curr
            curr = [line]
    if curr:
        yield curr


def parse_sequence(lines, constructor="".join):
    """Parses a GenBank sequence block. Doesn't care about ORIGIN line."""
    result = []
    exclude = b"0123456789 \t\n\r/"
    strip_table = dict([(c, None) for c in exclude])

    for i in lines:
        if i.startswith("ORIGIN"):
            continue

        result.append(i.translate(strip_table))
    return constructor(result)


def block_consolidator(lines):
    """Takes block with label and multiline data, and returns (label, [data]).

    [data] will be list of lines of data, including first line w/o label.
    """
    data = []
    first = True
    label = None
    for line in lines:
        if first:  # find label
            line = line.split(None, 1)
            if len(line) == 2:
                label, curr = line
            else:
                label = line[0]
                curr = ""
            data.append(curr)
            first = False
        else:
            data.append(line)
    return label, data


def parse_organism(lines):
    """Takes ORGANISM block. Returns organism, [taxonomy].

    NOTE: Adds species to end of taxonomy if identifiable.
    """
    label, data = block_consolidator(lines)
    # get 'species'
    species = data[0].strip()
    # get rest of taxonomy
    taxonomy = " ".join(data[1:])
    # normalize whitespace, including deleting newlines
    taxonomy = " ".join(taxonomy.split())
    # separate by semicolons
    # get rid of leading/trailing spaces
    taxa = list(map(strip, taxonomy.split(";")))
    # delete trailing period if present
    last = taxa[-1]
    if last.endswith("."):
        taxa[-1] = last[:-1]
    return species, taxa


def is_feature_component_start(line):
    """Checks if a line starts with '/', ignoring whitespace."""
    return line.lstrip().startswith("/")


feature_component_iterator = LabeledRecordFinder(is_feature_component_start)

_join_with_empty = dict.fromkeys(["translation"])
_leave_as_lines = {}


def parse_feature(lines):
    """Parses a feature. Doesn't handle subfeatures.

    Returns dict containing:
    'type': source, gene, CDS, etc.
    'location': unparsed location string
    ...then, key-value pairs for each annotation,
        e.g. '/gene="MNBH"' -> {'gene':['MNBH']} (i.e. quotes stripped)
    All relations are assumed 'to many', and order will be preserved.
    """
    result = {}
    type_, data = block_consolidator(lines)
    result["type"] = type_
    location = []
    found_feature = False
    for curr_line_idx, line in enumerate(data):
        if line.lstrip().startswith("/"):
            found_feature = True
            break
        else:
            location.append(line)
    result["raw_location"] = location
    try:
        result["location"] = parse_location_line(location_line_tokenizer(location))
    except (TypeError, ValueError):
        result["location"] = None
    if not found_feature:
        return result
    fci = feature_component_iterator
    for feature_component in fci(data[curr_line_idx:]):
        first = feature_component[0].lstrip()[1:]  # remove leading space, '/'
        try:
            label, first_line = first.split("=", 1)
        except ValueError:  # sometimes not delimited by =
            label, first_line = first, ""
        # chop off leading quote if appropriate
        if first_line.startswith('"'):
            first_line = first_line[1:]
        remainder = [first_line] + feature_component[1:]
        # chop off trailing quote, if appropriate
        last_line = remainder[-1].rstrip()
        if last_line.endswith('"'):
            remainder[-1] = last_line[:-1]
        if label in _join_with_empty:
            curr_data = "".join(map(strip, remainder))
        elif label in _leave_as_lines:
            curr_data = remainder
        else:
            curr_data = " ".join(map(strip, remainder))
        if label.lower() == "type":
            # some source features have /type=...
            label = "type_field"
        if label not in result:
            result[label.lower()] = []

        result[label.lower()].append(curr_data)
    return result


def location_line_tokenizer(lines):
    """Tokenizes location lines into spans, joins and complements."""
    curr = []
    text = " ".join(map(strip, lines))
    for char in text:
        if char == "(":
            yield "".join(curr).strip() + char
            curr = []
        elif char == ")":
            if curr:
                yield "".join(curr).strip()
            yield char
            curr = []
        elif char == ",":
            if curr:
                yield "".join(curr).strip()
            yield ","
            curr = []
        else:
            curr.append(char)
    if curr:
        yield "".join(curr).strip()


def parse_simple_location_segment(segment):
    """Parses location segment of form a..b or a, incl. '<' and '>'."""
    first_ambiguity, second_ambiguity = None, None
    if ".." in segment:
        first, second = segment.split("..")
        if not first[0].isdigit():
            first_ambiguity = first[0]
            first = int(first[1:])
        else:
            first = int(first)
        if not second[0].isdigit():
            second_ambiguity = second[0]
            second = int(second[1:])
        else:
            second = int(second)

        return Location(
            [
                Location(first, Ambiguity=first_ambiguity),
                Location(second, Ambiguity=second_ambiguity),
            ]
        )
    else:
        if not segment[0].isdigit():
            first_ambiguity = segment[0]
            segment = segment[1:]
        return Location(int(segment), Ambiguity=first_ambiguity)


def parse_location_line(tokens, parser=parse_simple_location_segment):
    """Parses location line tokens into location list."""
    stack = []
    curr = stack
    for t in tokens:
        if t.endswith("("):
            new = [curr, t]
            curr.append(new)
            curr = new
        elif t == ",":  # ignore
            continue
        elif t == ")":
            parent, type_ = curr[:2]
            children = curr[2:]
            if type_ == "complement(":
                children.reverse()
                for c in children:
                    c.Strand *= -1
            curr_index = parent.index(curr)
            del parent[curr_index]
            parent[curr_index:curr_index] = children[:]
            curr = parent
        else:
            curr.append(parser(t))
    return LocationList(stack)


class Location(object):
    """GenBank location object. Integer, or low, high, or 2-base bound.

    data must either be a long, an object that can be coerced to a long, or a
        sequence of two BasePosition objects. It can _not_ be two numbers.
    Ambiguity should be None (the default), '>', or '<'.
    IsBetween should be False (the default), or True.
    IsBounds should be False(the default, indicates range), or True.
    Accession should be an accession, or None (default).
    Db should be a database identifier, or None (default).
    Strand should be 1 (forward, default) or -1 (reverse).

    WARNING: This Location will allow you to do things that can't happen in
    GenBank, such as having a start and stop that aren't from the same
    accession. No validation is performed to prevent this. All reasonable
    cases should work.

    WARNING: Coordinates are based on 1, not 0, as in GenBank format.
    """

    def __init__(
        self,
        data,
        Ambiguity=None,
        IsBetween=False,
        IsBounds=False,
        Accession=None,
        Db=None,
        Strand=1,
    ):
        """Returns new LocalLocation object."""
        try:
            data = int(data)
        except TypeError:
            pass  # assume was two Location objects.
        self._data = data
        self.Ambiguity = Ambiguity
        self.IsBetween = IsBetween
        self.IsBounds = IsBounds
        self.Accession = Accession
        self.Db = Db
        self.Strand = Strand

    def __str__(self):
        """Returns self in string format.

        WARNING: More permissive than GenBank's Backus-Naur form allows. If
        you abuse this object, you'll get results that aren't valid GenBank
        locations.
        """
        if self.IsBetween:  # between two bases
            try:
                first, last = self._data
                curr = f"{first}^{last}"
            except TypeError:  # only one base? must be this or the next
                curr = f"{first}^{first + 1}"
        else:  # not self.IsBetween
            try:
                data = int(self._data)
                # if the above line succeeds, we've got a single item
                if self.Ambiguity:
                    curr = self.Ambiguity + str(data)
                else:
                    curr = str(data)
            except TypeError:
                # if long conversion failed, should have two LocalLocation
                # objects
                first, last = self._data
                if self.IsBounds:
                    curr = f"({first}{'.'}{last})"
                else:
                    curr = f"{first}{'..'}{last}"
        # check if we need to add on the accession and database
        if self.Accession:
            curr = self.Accession + ":" + curr
            # we're only going to add the Db if we got an accession
            if self.Db:
                curr = self.Db + "::" + curr
        # check if it's complemented
        if self.Strand == -1:
            curr = f"complement({curr})"
        return curr

    def first(self):
        """Returns first base self could be."""
        try:
            return int(self._data)
        except TypeError:
            return self._data[0].first()

    def last(self):
        """Returns last base self could be."""
        try:
            return int(self._data)
        except TypeError:
            return self._data[-1].last()


class LocationList(list):
    """List of Location objects.

    WARNING: Coordinates are based on 1, not 0, to match GenBank format.
    """

    BIGNUM = 1e300

    def first(self):
        """Returns first base of self."""
        curr = self.BIGNUM
        for i in self:
            first = i.first()
            if curr > first:
                curr = first
        return curr

    def last(self):
        """Returns last base of self."""
        curr = 0
        for i in self:
            last = i.last()
            if last > curr:
                curr = last
        return curr

    def strand(self):
        """Returns strand of components: 1=forward, -1=reverse, 0=both"""
        curr = {}
        for i in self:
            curr[i.Strand] = 1
        if len(curr) >= 2:  # found stuff on both strands
            return 0
        else:
            return list(curr.keys())[0]

    def __str__(self):
        """Returns (normalized) string representation of self."""
        if len(self) == 0:
            return ""
        elif len(self) == 1:
            return str(self[0])
        else:
            return "join(" + ",".join(map(str, self)) + ")"

    def extract(self, sequence, trans_table=dna_trans):
        """Extracts pieces of self from sequence."""
        result = []
        for i in self:
            first, last = i.first() - 1, i.last()  # inclusive, not exclusive
            # translate to 0-based indices and check if it wraps around
            if first < last:
                curr = sequence[first:last]
            else:
                curr = sequence[first:] + sequence[:last]
            # reverse-complement if necessary
            if i.Strand == -1:
                curr = curr.translate(trans_table)[::-1]
            result.append(curr)
        return "".join(result)


def parse_feature_table(lines):
    """Simple parser for feature table. Assumes starts with FEATURES line."""
    if not lines:
        return []
    if lines[0].startswith("FEATURES"):
        lines = lines[1:]
    return [parse_feature(f) for f in indent_splitter(lines)]


reference_label_marker = " " * 11
reference_field_finder = LabeledRecordFinder(
    lambda x: not x.startswith(reference_label_marker), constructor=None
)


def parse_reference(lines):
    """Simple parser for single reference."""
    result = {}
    for field in reference_field_finder(lines):
        label, data = block_consolidator(field)
        result[label.lower()] = " ".join(map(strip, data))
    return result


def parse_source(lines):
    """Simple parser for source fields."""
    result = {}
    all_lines = list(lines)
    source_field = next(reference_field_finder(all_lines))
    label, data = block_consolidator(source_field)
    result[label.lower()] = " ".join(map(strip, data))
    source_length = len(source_field)
    species, taxonomy = parse_organism(lines[source_length:])
    result["species"] = species
    result["taxonomy"] = taxonomy
    return result


# adaptors to update curr with data from each parser


def locus_adaptor(lines, curr):
    curr.update(parse_locus(lines[0]))


def source_adaptor(lines, curr):
    curr.update(parse_source(lines))


def ref_adaptor(lines, curr):
    if "references" not in curr:
        curr["references"] = []
    curr["references"].append(parse_reference(lines))


def feature_table_adaptor(lines, curr):
    if "features" not in curr:
        curr["features"] = []
    curr["features"].extend(parse_feature_table(lines))


def sequence_adaptor(lines, curr):
    curr["sequence"] = parse_sequence(lines)


def generic_adaptor(lines, curr):
    label, data = block_consolidator(lines)
    curr[label.lower()] = " ".join(map(strip, lines))


handlers = {
    "LOCUS": locus_adaptor,
    "SOURCE": source_adaptor,
    "REFERENCE": ref_adaptor,
    "FEATURES": feature_table_adaptor,
    "ORIGIN": sequence_adaptor,
    "//": lambda lines, curr: None,
    "?": lambda lines, curr: None,
}


def MinimalGenbankParser(lines, handlers=handlers, default_handler=generic_adaptor):
    for rec in GbFinder(lines):
        curr = {}
        bad_record = False
        for field in indent_splitter(rec):
            first_word = field[0].split(None, 1)[0]
            handler = handlers.get(first_word, default_handler)

            try:
                handler(field, curr)
            except:
                bad_record = True
                break

        if not bad_record:
            yield curr


def parse_location_segment(location_segment):
    """Parses a location segment into its component pieces.

    Known possibilities:
    http://www.ebi.ac.uk/embl/Documentation/FT_definitions/feature_table.html

    467             single base
    a..b            range from a to b, including a and b
    <a              strictly before a
    >a              strictly after a
    (a.b)           a single base between a and b, inclusive
    a^b             a site between two adjacent bases between a and b
    accession:a     a occurs in accession, not in the current sequence
    db::accession:a a occurrs in accession in db, not in the current sequence
    """
    s = location_segment  # save some typing...
    lsp = parse_location_segment
    # check if it's a range
    if ".." in s:
        first, second = s.split("..")
        return Location([lsp(first), lsp(second)])
    # check if it's between two adjacent bases
    elif "^" in s:
        first, second = s.split("^")
        return Location([lsp(first), lsp(second)], IsBetween=True)
    # check if it's a single base reference -- but don't be fooled by
    # accessions!
    elif "." in s and s.startswith("(") and s.endswith(")"):
        first, second = s.split(".")
        return Location([lsp(first[1:]), lsp(second[:-1])])


def parse_location_atom(location_atom):
    """Parses a location atom, supposed to be a single-base position."""
    a = location_atom
    if a.startswith("<") or a.startswith(">"):  # fuzzy
        position = int(a[1:])
        return Location(position, Ambiguity=a[0])
    # otherwise, should just be an integer
    return Location(int(a))


wanted_types = dict.fromkeys(["CDS"])


def extract_nt_prot_seqs(rec, wanted=wanted_types):
    """Extracts nucleotide seqs, and, where possible, protein seqs, from recs."""
    rec_seq = rec["sequence"]
    for f in rec["features"]:
        if f["type"] not in wanted:
            continue
        translation = f["translation"][0]
        raw_seq = f["location"].extract(rec_seq)
        print(raw_seq)
        seq = raw_seq[int(f["codon_start"][0]) - 1 :]
        print("dt:", translation)
        print("ct:", GeneticCodes[f.get("transl_table", "1")[0]].translate(seq))
        print("s :", seq)


def RichGenbankParser(
    handle, info_excludes=None, moltype=None, skip_contigs=False, add_annotation=None
):
    """Returns annotated sequences from GenBank formatted file.

    Parameters
    ----------
    info_excludes
        a series of fields to be excluded from the Info object
    moltype
        a MolType instance, such as PROTEIN, DNA. Default is ASCII.
    skip_contigs
        ignores records with no actual sequence data, typically
        a genomic contig.
    add_annotation
        a callback function to create an new annotation from a
        GenBank feature. Function is called with the sequence, a feature dict
        and the feature spans.

    """
    info_excludes = info_excludes or []
    moltype = get_moltype(moltype) if moltype else None
    for rec in MinimalGenbankParser(handle):
        info = Info()
        # populate the info object, excluding the sequence
        for label, value in list(rec.items()):
            if label in info_excludes:
                continue
            info[label] = value

        if moltype is None:
            rec_moltype = rec["mol_type"].lower()
            rec_moltype = (
                rec_moltype if rec_moltype in ("dna", "rna", "protein") else "text"
            )
            rec_moltype = get_moltype(rec_moltype)
        else:
            rec_moltype = moltype

        try:
            seq = rec_moltype.make_seq(
                rec["sequence"].upper(), info=info, name=rec["locus"]
            )
        except KeyError:
            if not skip_contigs:
                if "contig" in rec:
                    yield rec["locus"], rec["contig"]
                elif "WGS" in rec:
                    yield rec["locus"], rec["WGS"]
                else:
                    yield rec["locus"], None
            continue

        for feature in rec["features"]:
            spans = []
            reversed = None
            if feature["location"] is None or feature["type"] in ["source", "organism"]:
                continue
            for location in feature["location"]:
                (lo, hi) = (location.first() - 1, location.last())
                if location.Strand == -1:
                    (lo, hi) = (hi, lo)
                    assert reversed is not False
                    reversed = True
                else:
                    assert reversed is not True
                    reversed = False
                # ensure we don't put in a span that starts beyond the sequence
                if lo > len(seq):
                    continue
                # or that's longer than the sequence
                hi = [hi, len(seq)][hi > len(seq)]
                spans.append((lo, hi))

            if add_annotation:
                add_annotation(seq, feature, spans)
            else:
                for id_field in ["gene", "product", "clone", "note"]:
                    if id_field in feature:
                        name = feature[id_field]
                        if not isinstance(name, str):
                            name = " ".join(name)
                        break
                else:
                    name = None
                seq.add_annotation(Feature, feature["type"], name, spans)

        yield (rec["locus"], seq)


def parse(*args):
    return RichGenbankParser(*args).next()[1]
