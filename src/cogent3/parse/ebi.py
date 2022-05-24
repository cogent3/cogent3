#!/usr/bin/env python
"""Provide a parser for SwissProt EBI format files.
"""
import sys

from pprint import pprint

from cogent3.core.sequence import Sequence
from cogent3.parse.record import FieldError, RecordError
from cogent3.parse.record_finder import (
    DelimitedRecordFinder,
    LabeledRecordFinder,
    TailedRecordFinder,
)
from cogent3.util.misc import NestedSplitter, curry, list_flatten


__author__ = "Zongzhi Liu and Sandra Smit"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Zongzhi Liu",
    "Sandra Smit",
    "Rob Knight",
    "Gavin Huttley",
    "Daniel McDonald",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Zongzhi Liu"
__email__ = "zongzhi.liu@gmail.com"
__status__ = "Development"

maketrans = str.maketrans


def strip(x, chars=None):
    if chars:
        return x.strip(chars)
    else:
        return x.strip()


def rstrip(x, chars=None):
    if chars:
        return x.rstrip(chars)
    else:
        return x.rstrip()


all_chars = bytes(range(256))


def rstrip_(chars=None):
    return curry(rstrip, chars=chars)


EbiFinder = DelimitedRecordFinder("//", constructor=rstrip)

no_indent = lambda s: not s.startswith(" ")
hanging_paragraph_finder = LabeledRecordFinder(no_indent, constructor=None)

endswith_period = lambda x: x.endswith(".")
period_tail_finder = TailedRecordFinder(endswith_period)


#################################
# pairs_to_dict
def pairs_to_dict(
    key_values, dict_mode=None, all_keys=None, handlers=None, default_handler=None
):
    """generate a function which return a dict from a sequence of key_value
    pairs.

    key_values: (key, value) pairs, from any sequence type. Example:
    [('a', 1), ('b', 2), ('b', 3)]

    dict_mode: one of four modes to build a dict from pairs.
    'overwrite_value': default, same as dict(key_values)
    the key_values example get {'a': 1, 'b': 3}
    'no_duplicated_key': raise error when there is duplicated key;
    or a dict with a list of values when there is duplicated key
    'allow_muti_value': a duplicated key will have a list of values,
    the example get {'a': 1, 'b': [2, 3]}
    'always_multi_value': always group value(s) into a list for each key;
    the example get {'a': [1], 'b': [2, 3]}

    all_keys: if a key not found in all_keys, raise error; recommend to use a
    dict for all_keys for efficiency.

    Each value will be converted, if a valid handler can be found in handlers
    or a default_handler is provided.
    handler = handlers.get(key, default_handler)

    When handlers provided but no valid handler is found for a key: raise
    ValueError.
    Always use original value, if no handlers provided, for example,
    pairs_to_dict(adict.items()) will return adict.

    Note: use default_handler=identity is often useful if you want to return
    the original value when no handler found.
    """

    handlers = handlers or {}

    if not dict_mode:
        dict_mode = "overwrite_value"

    # generate add_item for different dict_mode.
    if dict_mode == "always_multi_value":

        def add_item(dictionary, key, value):
            """add key, value to dictionary in place"""
            dictionary.setdefault(key, []).append(value)

    elif dict_mode == "allow_multi_value":
        multiples = {}  # auxillary dict recording the keys with multi_values

        def add_item(dictionary, key, value):
            """add key, value to dictionary in place

            Warning: using outer auxillary dictionary: multiples"""
            if key in dictionary:
                if key not in multiples:
                    multiples[key] = True
                    dictionary[key] = [dictionary[key]]
                dictionary[key].append(value)
            else:
                dictionary[key] = value

    elif dict_mode == "no_duplicated_key":

        def add_item(dictionary, key, value):
            """add key, value to dictionary in place"""
            if key in dictionary:
                raise ValueError("Duplicated Key")
            dictionary[key] = value

    elif dict_mode == "overwrite_value":

        def add_item(dictionary, key, value):
            """add key, value to dictionary in place"""
            dictionary[key] = value

    else:  # unknown dict_mode
        raise ValueError(
            "Unknown dict_mode:%s. \ndict_mode must be one of "
            "overwrite_value, no_duplicated_key, allow_multi_value and "
            "always_multi_value." % dict_mode
        )

    # generate the handle_value function.
    if not handlers and not default_handler:
        handle_value = lambda x, y: (x, y)

    else:  # handlers not empty,

        def handle_value(key, raw_value):
            handler = handlers.get(key, default_handler)
            if handler:
                value = handler(raw_value)
            else:  # no handler found for key
                raise ValueError(f"No handler found for {key}")
            return key, value

    # build the result dict.
    result = {}
    for key, raw_value in key_values:
        if all_keys and key not in all_keys:
            raise ValueError(f"key: {repr(key)} not in all_keys: {all_keys}")
        key, value = handle_value(key, raw_value)
        add_item(result, key, value)
    return result


#################################
# generic parsers


def linecode_maker(line):
    """return the linecode and the line.

    The two-character line code that begins each line is always followed
    by three blanks, so that the actual information begins with the sixth
    character."""
    linecode = line.split("   ", 1)[0]
    return linecode, line


def labeloff(lines, splice_from=5):
    """strip off the first splice_from characters from each line

    Warning: without check!"""
    return [line[splice_from:] for line in lines]


def join_parser(lines, join_str=" ", chars_to_strip=" ;."):
    """return a joined str from a list of lines, strip off chars requested from
    the joined str"""
    # a str will not be joined
    if isinstance(lines, str):
        result = lines
    else:
        result = join_str.join(lines)

    return result.strip(chars_to_strip)


def join_split_parser(
    lines, delimiters=";", item_modifier=strip, same_level=False, **kwargs
):
    """return a nested list from lines, join lines before using NestedSplitter.

    delimiters: delimiters used by NestedSplitter
    item_modifier: passed to NestedSplitter, modify each splitted item.
    kwargs: passed to join_parser

    Examples:
    join_split_parser(['aa; bb;', 'cc.']) -> ['aa', 'bb', 'cc']
    join_split_parser(['aa; bb, bbb;', 'cc.'], delimiters=';,')
    -> ['aa', ['bb','bbb'], 'cc']
    join_split_parser('aa (bb) (cc).', delimiters='(',
    item_modifer=rstrip_(')) -> ['aa','bb','cc']
    """
    result = join_parser(lines, **kwargs)

    return NestedSplitter(delimiters, constructor=item_modifier, same_level=same_level)(
        result
    )


def join_split_dict_parser(
    lines, delimiters=None, dict_mode=None, strict=True, **kwargs
):
    """return a dict from lines, using the splited pairs from
    join_split_parser and pairs_to_dict.

    delimiters, kwargs: pass to join_split_sparser
    strict: when dict() fails -- (a pair not got from the second delimiter).
    return unconstructed list when False or raise error when True (default).

    dict_mode: pass to pairs_to_dict.  if leave as None, will be
    'overwrite_value', which is same as dict(pairs).

    Examples:
    join_split_dict_parser(['aa=1; bb=2,3; cc=4 (if aa=1);'])
    -> {'aa':'1', 'bb': ['2','3'], 'cc': '4 (if aa=1)'}
    """
    delimiters = delimiters or [";", ("=", 1), ","]
    primary_delimiters, value_delimiters = delimiters[:2], delimiters[2:]
    pairs = join_split_parser(
        lines, delimiters=primary_delimiters, same_level=True, **kwargs
    )

    try:
        dict(pairs)  # catch error for any not splitted pair.
    except ValueError:  # dictionary update sequence element #1 has length 1;
        if strict:
            raise ValueError(f"e\nFailed to get a dict from pairs: {pairs}")
        else:
            # return the splitted list without constucting
            return pairs

    if value_delimiters:
        split_value = NestedSplitter(value_delimiters, same_level=False)
        # should raise ValueError here if a pair donot have two elems.
        for i, (k, v) in enumerate(pairs):
            v = split_value(v)
            # modify v only if splitted by the first dilimiter
            if len(v) > 1:
                pairs[i][1] = v

    return pairs_to_dict(pairs, dict_mode)


def mapping_parser(line, fields, delimiters=None, flatten=list_flatten):
    """return a dict of zip(fields, splitted line), None key will be deleted
    from the result dict.

    line: should be a str,  to be splitted.
    fields: field name and optional type constructor for mapping.  example:
        ['EntryName', ('Length', int), 'moltype']
    delimiters: separators used to split the line.
    flatten: a function used to flatten the list from nested splitting.
    """
    delimiters = delimiters or [";", None]
    splits = NestedSplitter(delimiters=delimiters)(line)
    values = flatten(splits)
    result = {}
    for f, v in zip(fields, values):
        if isinstance(f, (tuple, list)):
            name, type = f
            result[name] = type(v)
        else:
            result[f] = v

    if None in result:
        del result[None]

    return result


#################################
# individual parsers
#################################

#################################
# mapping parsers: id, sq
#
def id_parser(lines):
    """return a mapping dict from id lines (only one line).

    The ID (IDentification) line is always the first line of an entry. The general
    form of the ID line is:

    ID   ENTRY_NAME DATA_CLASS; MOLECULE_TYPE; SEQUENCE_LENGTH.
    Example:
    ID   CYC_BOVIN      STANDARD;      PRT;   104 AA.
    """
    lines = labeloff(lines)
    return mapping_parser(
        lines[0],
        delimiters=[";", None],
        fields=("EntryName", "DataClass", "moltype", ("Length", int)),
    )


def sq_parser(lines):
    """return a mapping dict from SQ lines (only one line).

    The SQ (SeQuence header) line marks the beginning of the sequence data and
    gives a quick summary of its content.  The format of the SQ line is:

    SQ   SEQUENCE XXXX AA; XXXXX MW; XXXXXXXXXXXXXXXX CRC64;

    The line contains the length of the sequence in amino acids ('AA')
    followed by the molecular weight ('MW') rounded to the nearest mass unit
    (Dalton) and the sequence 64-bit CRC (Cyclic Redundancy Check) value
    ('CRC64').
    """
    lines = labeloff(lines)
    return mapping_parser(
        lines[0],
        delimiters=[";", None],
        fields=(None, ("Length", int), None, ("MolWeight", int), None, "Crc64"),
    )


def kw_parser(lines):
    """return a list of keywords from KW lines.

    The format of the KW line is:
    KW   Keyword[; Keyword...].
    """
    lines = labeloff(lines)
    return join_split_parser(lines)


def ac_parser(lines):
    """return a list of accession numbers from AC lines.

    The AC (ACcession number) line lists the accession number(s) associated
    with an entry. The format of the AC line is:
    AC   AC_number_1;[ AC_number_2;]...[ AC_number_N;]

    The first accession number is commonly referred to as the 'primary
    accession number'. 'Secondary accession numbers' are sorted
    alphanumerically.
    """
    lines = labeloff(lines)
    return join_split_parser(lines)


def dt_parser(lines):
    """return the origal lines from DT lines.

    Note: not complete parsing

    The DT (DaTe) lines show the date of creation and last modification of the
    database entry.
    The format of the DT line in Swiss-Prot is:
        DT   DD-MMM-YYYY (Rel. XX, Comment)
    The format of the DT line in TrEMBL is:
        DT   DD-MMM-YYYY (TrEMBLrel. XX, Comment)

    There are always three DT lines in each entry, each of them is associated
    with a specific comment:
    * The first DT line indicates when the entry first appeared in the
      database. The comment is 'Created';
    * The second DT line indicates when the sequence data was last modified.
      comment is 'Last sequence update';
    * The third DT line indicates when data (see the note below) other than the
      sequence was last modified. 'Last annotation update'.

    Example of a block of Swiss-Prot DT lines:
        DT   01-AUG-1988 (Rel. 08, Created)
        DT   30-MAY-2000 (Rel. 39, Last sequence update)
        DT   10-MAY-2005 (Rel. 47, Last annotation update)
    """
    lines = labeloff(lines)
    return lines


#################################
# gn_parser


def gn_parser(lines):
    """return a list of dict from GN lines.

    The GN (Gene name) line indicates the name(s) of the gene(s) that code for
    the stored protein sequence. The GN line contains three types of
    information: Gene names, Ordered locus names, ORF names. format:

    GN   name=<name>; Synonyms=<name1>[, <name2>...];
    OrderedLocusNames=<name1>[, <name2>...];
    GN   ORFNames=<name1>[, <name2>...];

    None of the above four tokens are mandatory. But a "Synonyms" token can
    only be present if there is a "name" token.

    If there is more than one gene, GN line blocks for the different genes are
    separated by the following line:
    GN   and
    Example:
    GN   name=Jon99Cii; Synonyms=SER1, SER5, Ser99Da; ORFNames=CG7877;
    GN   and
    GN   name=Jon99Ciii; Synonyms=SER2, SER5, Ser99Db; ORFNames=CG15519;"""
    lines = labeloff(lines)
    return list(map(gn_itemparser, gn_itemfinder(lines)))


gn_itemparser = join_split_dict_parser

gn_itemfinder = DelimitedRecordFinder(
    "and", constructor=None, strict=False, keep_delimiter=False
)


def oc_parser(lines):
    """return a list from OC lines.

    The OC (Organism Classification) lines contain the taxonomic classification
    of the source organism.  The classification is listed top-down as nodes in
    a taxonomic tree in which the most general grouping is given first. format:

    OC   Node[; Node...].
    """
    lines = labeloff(lines)
    return join_split_parser(lines)


def os_parser(lines):
    """return a list from OS lines.

    OS (Organism Species) line specifies the organism which was the source of
    the stored sequence.  The last OS line is terminated by a period.

    The species designation consists, in most cases, of the Latin genus and
    species designation followed by the English name (in parentheses). For
    viruses, only the common English name is given.

    Examples of OS lines are shown here:
        OS   Escherichia coli.
        OS   Solanum melongena (Eggplant) (Aubergine).
        OS   Rous sarcoma virus (strain SchRuppin A) (RSV-SRA) (Avian leukosis
        OS   virus-RSA).
    """
    lines = labeloff(lines)
    return join_split_parser(lines, delimiters="(", item_modifier=rstrip_(") "))


def ox_parser(lines):
    """return a dict from OX lines.

    The OX (Organism taxonomy cross-reference) line is used to indicate the
    identifier of a specific organism in a taxonomic database. The format:
    OX   Taxonomy_database_Qualifier=Taxonomic code;

    Currently the cross-references are made to the taxonomy database of NCBI,
    which is associated with the qualifier 'TaxID' and a one- to six-digit
    taxonomic code."""
    lines = labeloff(lines)
    return join_split_dict_parser(lines)


def og_parser(lines):
    """return a list from OG lines
    The OG (OrGanelle) line indicates if the gene coding for a protein
    originates from the mitochondria, the chloroplast, the cyanelle, the
    nucleomorph or a plasmid.  The format of the OG line is:

    OG   Hydrogenosome.
    OG   Mitochondrion.
    OG   Nucleomorph.
    OG   Plasmid name.
    OG   Plastid.
    OG   Plastid; Apicoplast.
    OG   Plastid; Chloroplast.
    OG   Plastid; Cyanelle.
    OG   Plastid; Non-photosynthetic plastid.
    Where 'name' is the name of the plasmid. example:

    OG   Mitochondrion.
    OG   Plasmid R6-5, Plasmid IncFII R100 (NR1), and
    OG   Plasmid IncFII R1-19 (R1 drd-19)."""
    lines = labeloff(lines)
    result = []
    for item in period_tail_finder(lines):
        item = " ".join(item).rstrip(". ")
        if item.startswith("Plasmid"):
            item = item.replace(" and", "")
            item = list(map(strip, item.split(",")))
        result.append(item)
    return result


#################################
# dr_parser
def dr_parser(lines):
    """return a dict of items from DR lines.

    The DR (Database cross-Reference) lines are used as pointers to information
    related to entries and found in data collections other than Swiss-Prot.
    The format of one of many DR line is:

    DR   DATABASE_IDENTIFIER; PRIMARY_IDENTIFIER; SECONDARY_IDENTIFIER[;
    TERTIARY_IDENTIFIER][; QUATERNARY_IDENTIFIER].
    """
    lines = labeloff(lines)
    keyvalues = list(map(dr_itemparser, period_tail_finder(lines)))
    return pairs_to_dict(keyvalues, "always_multi_value")


def dr_itemparser(lines):
    """return a key, value pair from lines of a DR item."""
    fields = join_split_parser(lines)
    return fields[0], fields[1:]


#################################
# de_parser


def de_parser(lines):
    """return a dict of {OfficalName: str, Synonyms: str, Fragment: bool,
    Contains: [itemdict,],  Includes: [itemdict,]} from DE lines

    The DE (DEscription) lines contain general descriptive information about
    the sequence stored. This information is generally sufficient to identify
    the protein precisely.

    The description always starts with the proposed official name of the
    protein. Synonyms are indicated between brackets. Examples below

    If a protein is known to be cleaved into multiple functional components,
    the description starts with the name of the precursor protein, followed by
    a section delimited by '[Contains: ...]'. All the individual components are
    listed in that section and are separated by semi-colons (';'). Synonyms are
    allowed at the level of the precursor and for each individual component.

    If a protein is known to include multiple functional domains each of which
    is described by a different name, the description starts with the name of
    the overall protein, followed by a section delimited by '[Includes: ]'. All
    the domains are listed in that section and are separated by semi-colons
    (';'). Synonyms are allowed at the level of the protein and for each
    individual domain.

    In rare cases, the functional domains of an enzyme are cleaved, but the
    catalytic activity can only be observed, when the individual chains
    reorganize in a complex. Such proteins are described in the DE line by a
    combination of both '[Includes:...]' and '[Contains:...]', in the order
    given in the following example:

    If the complete sequence is not determined, the last information given on
    the DE lines is '(Fragment)' or '(Fragments)'. Example:

    DE   Dihydrodipicolinate reductase (EC 1.3.1.26) (DHPR) (Fragment).

    DE   Arginine biosynthesis bifunctional protein argJ [Includes: Glutamate
    DE   N-acetyltransferase (EC 2.3.1.35) (Ornithine acetyltransferase)
    DE   (Ornithine transacetylase) (OATase); Amino-acid acetyltransferase
    DE   (EC 2.3.1.1) (N-acetylglutamate synthase) (AGS)] [Contains: Arginine
    DE   biosynthesis bifunctional protein argJ alpha chain; Arginine
    DE   biosynthesis bifunctional protein argJ beta chain] (Fragment).

    Trouble maker:
    DE Amiloride-sensitive amine oxidase [copper-containing] precursor(EC
    DE 1.4.3.6) (Diamine oxidase) (DAO).
    """
    labeloff_lines = labeloff(lines)
    joined = join_parser(labeloff_lines, chars_to_strip="). ")

    keys = ["Includes", "Contains", "Fragment"]
    fragment_label = "(Fragment"
    contains_label = "[Contains:"
    includes_label = "[Includes:"

    # Process Fragment
    fragment = False
    if joined.endswith(fragment_label):
        fragment = True
        joined = joined.rsplit("(", 1)[0]

    # Process Contains
    contains = []
    if contains_label in joined:
        joined, contains_str = joined.split(contains_label)
        contains_str = contains_str.strip(" ]")
        contains = list(map(de_itemparser, contains_str.split("; ")))
    # Process Includes
    includes = []
    if includes_label in joined:
        joined, includes_str = joined.split(includes_label)
        includes_str = includes_str.strip(" ]")
        includes = list(map(de_itemparser, includes_str.split("; ")))

    # Process Primary
    primary = de_itemparser(joined)

    result = dict(list(zip(keys, (includes, contains, fragment))))
    result.update(primary)
    return result


def de_itemparser(line):
    """return a dict of {OfficalName: str, Synonyms: [str,]} from a de_item

    The description item is a str, always starts with the proposed official
    name of the protein. Synonyms are indicated between brackets. Examples
    below

    'Annexin A5 (Annexin V) (Lipocortin V) (Endonexin II)'
    """
    fieldnames = ["OfficalName", "Synonyms"]
    fields = [e.strip(") ") for e in line.split("(")]
    # if no '(', fields[1:] will be []
    return dict(list(zip(fieldnames, [fields[0], fields[1:]])))


def pr_parser(line):
    """Returns a list of [project id, project id, ...]"""
    labeloff(line)
    return [r.strip().split(":")[-1] for r in line.split(";") if r]


#################################
# ft_parser


def ft_parser(lines):
    """return a list of ft items from FT lines.

    The FT (Feature Table) lines lists posttranslational modifications, binding
    sites, enzyme active sites, local secondary structure or other
    characteristics reported in the cited references. Sequence conflicts
    between references are also included in the feature table.

    The FT lines have a fixed format. The column numbers allocated to each of
    the data items within each FT line are shown in the following table (column
    numbers not referred to in the table are always occupied by blanks).

        Columns     Data item
        1-2     FT
        6-13    Key name
        15-20   'From' endpoint
        22-27   'To' endpoint
        35-75   Description

    The key name and the endpoints are always on a single line, but the
    description may require one or more additional lines. The following
    description lines continues from column 35 onwards.  For more information
    about individual ft keys, see
    http://us.expasy.org/sprot/userman.html#FT_keys

    'FROM' and 'TO' endpoints: Numbering start from 1; When a feature is known
    to extend beyond the position that is given in the feature table, the
    endpoint specification will be preceded by '<' for features which continue
    to the left end (N-terminal direction) or by '>' for features which
    continue to the right end (C- terminal direction); Unknown endpoints are
    denoted by '?'. Uncertain endpoints are denoted by a '?' before the
    position, e.g. '?42'.

    Some features (CARBOHYD, CHAIN, PEPTIDE, PROPEP, VARIANT and VARSPLIC) are
    associated with a unique and stable feature identifier (FTId), which allows
    to construct links directly from position-specific annotation in the
    feature table to specialized protein-related databases.  The FTId is always
    the last component of a feature in the description field.

    Examples:
        FT   SIGNAL       <1     10       By similarity.
        FT   MOD_RES      41     41       Arginine amide (G-42 provides amide
        FT                                group) (By similarity).
        FT   CONFLICT    327    327       E -> R (in Ref. 2).
        FT   CONFLICT     77     77       Missing (in Ref. 1).
        FT   CARBOHYD    251    251       N-linked (GlcNAc...).
        FT                                /FTId=CAR_000070.
        FT   PROPEP       25     48
        FT                                /FTId=PRO_0000021449.
        FT   VARIANT     214    214       V -> I.
        FT                                /FTId=VAR_009122.
        FT   VARSPLIC     33     83       TVGRFRRRATP -> PLTSFHPFTSQMPP (in
        FT                                isoform 2).
        FT                                /FTId=VSP_004370.

    Secondary structure (HELIX, STRAND, TURN) - The feature table of sequence
    entries of proteins whose tertiary structure is known experimentally
    contains the secondary structure information extracted from the coordinate
    data sets of the Protein Data Bank (PDB).  Residues not specified in one of
    these classes are in a 'loop' or 'random-coil' structure.
    """
    lines = labeloff(lines)
    fieldnames = "Start End Description".split()
    secondary_structure_keynames = "HELIX STRAND TURN".split()
    result = {}
    for item in hanging_paragraph_finder(lines):
        keyname, start, end, description = ft_basic_itemparser(item)

        # group secondary structures (as a list) into
        # result['SecondaryStructure']
        if keyname in secondary_structure_keynames:
            result.setdefault("SecondaryStructure", []).append((keyname, start, end))
            continue

        # further parser the description for certain keynames
        if keyname in ft_description_parsers:
            description = ft_description_parsers[keyname](description)

        # group current item result (as a dict) into result[keyname]
        curr = dict(list(zip(fieldnames, [start, end, description])))
        result.setdefault(keyname, []).append(curr)
    return result


def ft_basic_itemparser(item_lines):
    """-> (key, start, end, description) from lines of a feature item.

    A feature item (generated by itemfinder) has the same keyname.

    WARNING: not complete, location fields need further work?
    """
    # cut_postions: the postions to split the line into fields
    original_cut_positions = [15, 22, 35]  # see doc of ft_parser
    # keyname will start from 0(instead of 6) after labeloff
    cut_positions = [e - 6 for e in original_cut_positions]

    # unpack the first line to fields
    first_line = item_lines[0]
    keyname, from_point, to_point, description = [
        first_line[i:j].strip()
        for i, j in zip([0] + cut_positions, cut_positions + [None])
    ]

    # extend the description if provided following lines
    if len(item_lines) > 1:
        following_lines = item_lines[1:]
        desc_start = cut_positions[-1]
        following_description = " ".join(
            [e[desc_start:].strip() for e in following_lines]
        )
        description = " ".join((description, following_description))

    # convert start and end points to int, is possible
    from_point, to_point = list(map(try_int, (from_point, to_point)))
    return keyname, from_point, to_point, description.strip(" .")


def try_int(obj):
    """return int(obj), or original obj if failed"""
    try:
        return int(obj)
    except ValueError:  # invalid literal for int()
        return obj


# ft description_parsers below


def ft_id_parser(description):
    """return a dict of {'Description':,'Id':} from raw decription str

    Examples.
    FT   PROPEP       25     48
    FT                                /FTId=PRO_0000021449.
    FT   VARIANT     214    214       V -> I.
    FT                                /FTId=VAR_009122.
    FT   VARSPLIC     33     83       TVGRFRRRATP -> PLTSFHPFTSQMPP (in
    FT                                isoform 2).
    FT                                /FTId=VSP_004370.
    """
    fieldnames = ["Description", "Id"]
    id_sep = "/FTId="
    try:
        desc, id = [i.strip(" .") for i in description.split(id_sep)]
    except:
        desc, id = description, ""

    # replace desc in fields with (desc, id) to get the result
    result = dict(list(zip(fieldnames, [desc, id])))
    return result


def ft_mutation_parser(description, mutation_comment_delimiter="("):
    """return a  dict of {'MutateFrom': , 'MutateTo':,'Comment':} from
    description str

    Warning: if both id and mutation should be parsed, always parse id first.

    Note: will add exceptions later

    Examples.
    FT   VARIANT     214    214       V -> I (in a lung cancer).
    FT                                /FTId=VAR_009122.
    FT   CONFLICT    484    484       Missing (in Ref. 2).
    FT   CONFLICT    802    802       K -> Q (in Ref. 4, 5 and 10).
    """
    fieldnames = "MutateFrom MutateTo Comment".split()

    # split desc into mutation and comment
    desc = description.rstrip(" )")
    try:
        mutation, comment = desc.split(mutation_comment_delimiter, 1)
    except ValueError:  # too many values to unpack
        mutation, comment = desc, ""

    # split mutation into mut_from, mut_to
    # if mut_from/to unknown, the mutation message will be in mut_from
    mutation_delimiter = "->"
    try:
        mut_from, mut_to = list(map(strip, mutation.split(mutation_delimiter, 1)))
    except ValueError:  # too many values to unpack
        mut_from, mut_to = mutation, ""

    # replace desc in fields with mut_from, mut_to and comment to get the
    # result
    result = dict(list(zip(fieldnames, [mut_from, mut_to, comment])))
    return result


def ft_mutagen_parser(description):
    """return a dict from MUTAGEN description

    MUTAGEN - Site which has been experimentally altered.  Examples

    FT   MUTAGEN     119    119       C->R,E,A: Loss of cADPr hydrolase and
    FT                                ADP-ribosyl cyclase activity.
    FT   MUTAGEN     169    177       Missing: Abolishes ATP-binding.
    """
    return ft_mutation_parser(description, mutation_comment_delimiter=":")


def ft_id_mutation_parser(description):
    """return a dict from description str

    Examples.
    FT   VARIANT     214    214       V -> I.
    FT                                /FTId=VAR_009122.
    FT   VARSPLIC     33     83       TVGRFRRRATP -> PLTSFHPFTSQMPP (in
    FT                                isoform 2).
    FT                                /FTId=VSP_004370.
    """
    desc_id_dict = ft_id_parser(description)
    desc = desc_id_dict.pop("Description")
    return dict(desc_id_dict, **ft_mutation_parser(desc))


ft_description_parsers = {
    "VARIANT": ft_id_mutation_parser,
    "VARSPLIC": ft_id_mutation_parser,
    "CARBOHYD": ft_id_parser,
    "CHAIN": ft_id_parser,
    "PEPTIDE": ft_id_parser,
    "PROPEP": ft_id_parser,
    "CONFLICT": ft_mutation_parser,
    "MUTAGEN": ft_mutagen_parser,
}


#################################
# cc_parser
all_cc_topics = dict.fromkeys(
    [
        "ALLERGEN",
        "ALTERNATIVE PRODUCTS",
        "BIOPHYSICOCHEMICAL PROPERTIES",
        "BIOTECHNOLOGY",
        "CATALYTIC ACTIVITY",
        "CAUTION",
        "COFACTOR",
        "DATABASE",
        "DEVELOPMENTAL STAGE",
        "DISEASE",
        "DOMAIN",
        "ENZYME REGULATION",
        "FUNCTION",
        "INDUCTION",
        "INTERACTION",
        "MASS SPECTROMETRY",
        "MISCELLANEOUS",
        "PATHWAY",
        "PHARMACEUTICAL",
        "POLYMORPHISM",
        "PTM",
        "RNA EDITING",
        "SIMILARITY",
        "SUBCELLULAR LOCATION",
        "SUBUNIT",
        "TISSUE SPECIFICITY",
        "TOXIC DOSE",
    ]
)


def cc_parser(lines, strict=False):
    """return a dict of {topic: a list of values} from CC lines.

    some topics have special format and will use specific parsers defined in
    handlers

    The CC lines are free text comments on the entry, and are used to convey
    any useful information. The comments always appear below the last reference
    line and are grouped together in comment blocks; a block is made up of 1 or
    more comment lines. The first line of a block starts with the characters
    '-!-'.  Format:

    CC   -!- TOPIC: First line of a comment block;
    CC       second and subsequent lines of a comment block.

    Examples:
    CC   -!- DISEASE: Defects in PHKA1 are linked to X-linked muscle
    CC       glycogenosis [MIM:311870]. It is a disease characterized by slowly
    CC       progressive, predominantly distal muscle weakness and atrophy.
    CC   -!- DATABASE: NAME=Alzheimer Research Forum; NOTE=APP mutations;
    CC       WWW="http://www.alzforum.org/res/com/mut/app/default.asp".
    CC   -!- INTERACTION:
    CC       Self; NbExp=1; IntAct=EBI-123485, EBI-123485;
    CC       Q9W158:CG4612; NbExp=1; IntAct=EBI-123485, EBI-89895;
    CC       Q9VYI0:fne; NbExp=1; IntAct=EBI-123485, EBI-126770;
    CC   -!- BIOPHYSICOCHEMICAL PROPERTIES:
    CC       Kinetic parameters:
    CC         KM=98 uM for ATP;
    CC         KM=688 uM for pyridoxal;
    CC         Vmax=1.604 mmol/min/mg enzyme;
    CC       pH dependence:
    CC         Optimum pH is 6.0. Active from pH 4.5 to 10.5;
    CC   -!- ALTERNATIVE PRODUCTS:
    CC       Event=Alternative splicing; Named isoforms=3;
    CC         Comment=Additional isoforms seem to exist. Experimental
    CC         confirmation may be lacking for some isoforms;
    CC       name=1; Synonyms=AIRE-1;
    CC         IsoId=O43918-1; Sequence=Displayed;
    CC       name=2; Synonyms=AIRE-2;
    CC         IsoId=O43918-2; Sequence=VSP_004089;
    CC       name=3; Synonyms=AIRE-3;
    CC         IsoId=O43918-3; Sequence=VSP_004089, VSP_004090;
    CC   --------------------------------------------------------------------------
    CC   This SWISS-PROT entry is copyright. It is produced  a collaboration
    CC   removed.
    CC   --------------------------------------------------------------------------
    """
    lines = labeloff(lines)
    # cc_itemfinder yield each topic block
    # cc_basic_itemparser split a topic block into (topic_name,
    # content_as_list)
    topic_contents = list(map(cc_basic_itemparser, cc_itemfinder(lines)))
    # content of a topic further parsed using a content_parser decided by the
    # topic name.  result is grouped into a dict.
    try:
        result = pairs_to_dict(
            topic_contents,
            "always_multi_value",
            handlers=cc_content_parsers,
            default_handler=join_parser,
        )
    except Exception as e:
        pprint(lines)
        raise e

    if strict:
        for topic in result:
            if topic not in all_cc_topics:
                raise FieldError(f"Invalid topic: {topic}")

    return result


def cc_basic_itemparser(topic):
    """return (topic_name, topic_content as a list) from a cc topic block.

    Format of a topic as input of this function: [
    '-!- TOPIC: First line of a comment block;',
    '    second and subsequent lines of a comment block.']
    """
    num_format_leading_spaces = 4  # for topic lines except the first

    # get the keyname and content_head from the first line
    topic_head = topic[0].lstrip(" -!")
    try:
        keyname, content_head = list(map(strip, topic_head.split(":", 1)))
    except ValueError:  # need more than 1 value to unpack
        raise FieldError("Not a valid topic line: %s", topic[0])

    if content_head:
        content = [content_head]
    else:
        content = []

    # the following lines be stripped off the format leading spaces
    if len(topic) > 1:
        content += labeloff(topic[1:], num_format_leading_spaces)

    return keyname, content


def cc_itemfinder(lines):
    """yield each topic/license as a list from CC lines without label and
    leading spaces.

    Warning: hardcoded LICENSE handling"""

    # all the codes except the return line  tries to preprocess the
    # license block

    # two clusters of '-' are used as borders for license, as observed
    license_border = "-" * 74
    license_headstr = "-!- LICENSE:"
    content_start = 4  # the idx where topic content starts

    if license_border in lines:
        # discard the bottom license border
        if lines[-1] == license_border:
            lines.pop()
        else:
            raise FieldError(f"No bottom line for license: {lines}")

        # normalize license lines to the format of topic lines
        license_idx = lines.index(license_border)
        lines[license_idx] = license_headstr
        for i in range(license_idx + 1, len(lines)):
            lines[i] = " " * content_start + lines[i]

    # the return line is all we need, if no license block
    return hanging_paragraph_finder(lines)


# cc_content_parsers here below


def cc_interaction_parser(content_list):
    """return a list of [interactor, {params}] from interaction content.

    Format:
    -!- INTERACTION:
        {{SP_Ac:identifier[ (xeno)]}|Self}; NbExp=n; IntAct=IntAct_Protein_Ac,
        IntAct_Protein_Ac;
    """
    result = []
    for line in content_list:
        interactor, params = line.split(";", 1)
        params = join_split_dict_parser([params])
        result.append((interactor.strip(), params))
    return result


cc_alternative_products_event_finder = LabeledRecordFinder(
    lambda x: x.startswith("Event=")
)
cc_alternative_products_name_finder = LabeledRecordFinder(
    lambda x: x.startswith("name=")
)


def cc_alternative_products_parser(content_list):
    """return a list from AlternativeProucts lines.

    Note: not complete parsing, consider to merge Names to the last Event??
    and make event or name to be the dict key?

    Format:
    CC   -!- ALTERNATIVE PRODUCTS:
    CC       Event=Alternative promoter;
    CC         Comment=Free text;
    CC       Event=Alternative splicing; Named isoforms=n;
    CC         Comment=Optional free text;
    CC       name=Isoform_1; Synonyms=Synonym_1[, Synonym_n];
    CC         IsoId=Isoform_identifier_1[, Isoform_identifier_n]; Sequence=Displayed;
    CC         Note=Free text;
    CC       name=Isoform_n; Synonyms=Synonym_1[, Synonym_n];
    CC         IsoId=Isoform_identifier_1[, Isoform_identifier_n]; Sequence=VSP_identifier_1 [, VSP_identifier_n];
    CC         Note=Free text;
    CC       Event=Alternative initiation;
    CC         Comment=Free text;
    """
    result = []
    for event in cc_alternative_products_event_finder(content_list):
        head_names = list(cc_alternative_products_name_finder(event))
        head, names = head_names[0], head_names[1:]
        event_dict = join_split_dict_parser(head)
        if names:
            event_dict["Names"] = list(map(join_split_dict_parser, names))
        result.append(event_dict)
    return result


def cc_biophysicochemical_properties_parser(content):
    """return a dict from content_list of a ~ topic.

    Format of a ~ topic block:
    CC   -!- BIOPHYSICOCHEMICAL PROPERTIES:
    CC       Absorption:
    CC         Abs(max)=xx nm;
    CC         Note=free_text;
    CC       Kinetic parameters:
    CC         KM=xx unit for substrate [(free_text)];
    CC         Vmax=xx unit enzyme [free_text];
    CC         Note=free_text;
    CC       pH dependence:
    CC         free_text;
    CC       Redox potential:
    CC         free_text;
    CC       Temperature dependence:
    CC         free_text;

    Example of a ~ topic block:
    CC   -!- BIOPHYSICOCHEMICAL PROPERTIES:
    CC       Kinetic parameters:
    CC         KM=98 uM for ATP;
    CC         KM=688 uM for pyridoxal;
    CC         Vmax=1.604 mmol/min/mg enzyme;
    CC       pH dependence:
    CC         Optimum pH is 6.0. Active from pH 4.5 to 10.5;
    """

    def get_sub_key_content(sub_topic):
        """return (sub_key, sub_content as parsed) from lines of a sub_topic"""
        sub_key = sub_topic[0].rstrip(": ")
        # strip the two leading spaces
        sub_content = list(map(strip, sub_topic[1:]))

        # further process the content here
        if sub_key in ["Kinetic parameters", "Absorption"]:
            # group into a dict which allow multiple values.
            subkey_values = join_split_parser(sub_content, delimiters=[";", ("=", 1)])
            sub_content = pairs_to_dict(subkey_values, "allow_multi_value")
        else:
            sub_content = join_parser(sub_content, chars_to_strip="; ")
        return sub_key, sub_content

    sub_key_contents = list(map(get_sub_key_content, hanging_paragraph_finder(content)))
    return pairs_to_dict(sub_key_contents, "no_duplicated_key")


cc_content_parsers = {
    # ? not complete: further group alternative splicing?
    "ALTERNATIVE PRODUCTS": cc_alternative_products_parser,
    "BIOPHYSICOCHEMICAL PROPERTIES": cc_biophysicochemical_properties_parser,
    "INTERACTION": cc_interaction_parser,
    "DATABASE": join_split_dict_parser,
    "MASS SPECTROMETRY": join_split_dict_parser,
}


#################################
# REFs parser
def refs_parser(lines):
    """return a dict of {RN: single_ref_dict}

    These lines comprise the literature citations. The citations indicate the
    sources from which the data has been abstracted.  if several references are
    given, there will be a reference block for each.
    """
    rn_ref_pairs = list(map(single_ref_parser, ref_finder(lines)))
    return pairs_to_dict(rn_ref_pairs)


is_ref_line = lambda x: x.startswith("RN")
ref_finder = LabeledRecordFinder(is_ref_line)

required_ref_labels = "RN RP RL RA/RG RL".split()


def single_ref_parser(lines, strict=False):
    """return rn, ref_dict from lines of a single reference block

    strict: if True (default False), raise RecordError if lacking required
    labels.

    Warning: using global required_ref_labels.

    The reference lines for a given citation occur in a block, and are always
    in the order RN, RP, RC, RX, RG, RA, RT and RL. Within each such reference
    block, the RN line occurs once, the RC, RX and RT lines occur zero or more
    times, and the RP, RG/RA and RL lines occur one or more times.
    """
    # group by linecode
    label_lines = list(map(linecode_maker, lines))
    raw_dict = pairs_to_dict(label_lines, "always_multi_value")

    if strict:
        labels = dict.fromkeys(list(raw_dict.keys()))
        if "RA" in labels or "RG" in labels:
            labels["RA/RG"] = True
        for rlabel in required_ref_labels:
            if rlabel not in labels:
                raise RecordError(f"The reference block lacks required label: {rlabel}")

    # parse each field with relevant parser
    parsed_dict = pairs_to_dict(list(raw_dict.items()), handlers=ref_parsers)
    rn = parsed_dict.pop("RN")

    return rn, parsed_dict


# ref_parsers here below


def rx_parser(lines):
    """return a dict from RX lines.

    The RX (Reference cross-reference) line is an optional line which is used
    to indicate the identifier assigned to a specific reference in a
    bibliographic database. The format:
    RX   Bibliographic_db=IDENTIFIER[; Bibliographic_db=IDENTIFIER...];

    Where the valid bibliographic database names and their associated
    identifiers are:
       MEDLINE    Eight-digit MEDLINE Unique Identifier (UI)
       PubMed    PubMed Unique Identifier (PMID)
       DOI  Digital Object Identifier (DOI), examples:
           DOI=10.2345/S1384107697000225
           DOI=10.4567/0361-9230(1997)42:<OaEoSR>2.0.TX;2-B
           http://www.doi.org/handbook_2000/enumeration.html#2.2
    """
    lines = labeloff(lines)
    return join_split_dict_parser(lines, delimiters=["; ", "="])


def rc_parser(lines):
    """return a dict from RC lines.

    The RC (Reference Comment) lines are optional lines which are used to store
    comments relevant to the reference cited. The format:
    RC   TOKEN1=Text; TOKEN2=Text; ...

    The currently defined tokens and their order in the RC line are:
    STRAIN TISSUE TRANSPOSON PLASMID
    """
    lines = labeloff(lines)
    return join_split_dict_parser(lines)


def rg_parser(lines):
    """return a list of str(group names) from RG lines

    The Reference Group (RG) line lists the consortium name associated with a
    given citation. The RG line is mainly used in submission reference blocks,
    but can also be used in paper references, if the working group is cited as
    an author in the paper. RG line and RA line (Reference Author) can be
    present in the same reference block; at least one RG or RA line is
    mandatory per reference block. example:
    RG   The mouse genome sequencing consortium;
    """
    lines = labeloff(lines)
    return join_split_parser(lines)


def ra_parser(lines):
    """return a list from RA lines.

    The RA (Reference Author) lines list the authors of the paper (or other
    work) cited. RA might be missing in references that cite a reference group
    (see RG line). At least one RG or RA line is mandatory per reference block.

    All of the authors are included, and are listed in the order given in the
    paper. The names are listed surname first followed by a blank, followed by
    initial(s) with periods. The authors' names are separated by commas and
    terminated by a semicolon. Author names are not split between lines. eg:
    RA   Galinier A., Bleicher F., Negre D., Perriere G., Duclos B.;

    All initials of the author names are indicated and hyphens between initials
    are kept.  An author's initials can be followed by an abbreviation such as
    'Jr' (for Junior), 'Sr' (Senior), 'II', 'III' or 'IV'.  Example:
    RA   Nasoff M.S., Baker H.V. II, Wolf R.E. Jr.;
    """
    lines = labeloff(lines)
    return join_split_parser(lines, chars_to_strip=";", delimiters=",")


def rp_parser(lines):
    """return joined str stripped of '.'.

    The RP (Reference Position) lines describe the extent of the work relevant
    to the entry carried out by the authors. format:
    RP   COMMENT.
    """
    lines = labeloff(lines)
    return " ".join(lines).strip(". ")


def rl_parser(lines):
    """return joined str stipped of '.'.

    Note: not complete parsing.

    The RL (Reference Location) lines contain the conventional citation
    information for the reference. In general, the RL lines alone are
    sufficient to find the paper in question.

    a) Journal citations
    RL   Journal_abbrev Volume:First_page-Last_page(YYYY).
    When a reference is made to a paper which is 'in press'
    RL   Int. J. Parasitol. 0:0-0(2005).

    b) Electronic publications
    includes an '(er)' prefix. The format is indicated below:
    RL   (er) Free text.

    c) Book citations
    RL   (In) Editor_1 I.[, Editor_2 I., Editor_X I.] (eds.);
    RL   Book_name, pp.[Volume:]First_page-Last_page, Publisher, City (YYYY).
    Examples:
    RL   (In) Rich D.H., Gross E. (eds.);
    RL   Proceedings symposium, pp.69-72, Pierce
    RL   Chemical Co., Rockford Il. (1981).

    d) Unpublished results, eg:
    RL   Unpublished results, cited by:
    RL   Shelnutt J.A., Rousseau D.L., Dethmers J.K., Margoliash E.;
    RL   Biochemistry 20:6485-6497(1981).

    e) Unpublished observations, format:
    RL   Unpublished observations (MMM-YYYY).

    f) Thesis, format:
    RL   Thesis (Year), Institution_name, Country.

    g) Patent applications, format:
    RL   Patent number Pat_num, DD-MMM-YYYY.

    h) Submissions, format:
    RL   Submitted (MMM-YYYY) to Database_name.
    'Database_name' is one of the following:
    EMBL/GenBank/DDBJ, Swiss-Prot, PDB, PIR.
    """
    lines = labeloff(lines)
    return " ".join(lines).strip(". ")


def rt_parser(lines):
    """return joined line stripped of .";

    The RT (Reference Title) lines give the title of the paper (or other work)
    cited as exactly as possible given the limitations of the computer
    character set. The format of the RT line is:
    RT   "Title.";"""
    lines = labeloff(lines)
    return " ".join(lines).strip('.";')


def rn_parser(lines):
    """return a integer from RN lines (only one line).

    The RN (Reference Number) line gives a sequential number to each reference
    citation in an entry. This number is used to indicate the reference in
    comments and feature table notes. The format of the RN line is:
    RN   [##]
    """
    lines = labeloff(lines)
    return int(lines[0].strip(" []"))


ref_parsers = {
    "RN": rn_parser,
    "RP": rp_parser,
    "RC": rc_parser,
    "RX": rx_parser,
    "RG": rg_parser,
    "RA": ra_parser,
    "RT": rt_parser,
    "RL": rl_parser,
}

required_labels = "ID AC DT DE OS OC OX SQ REF".split() + [""]
#################################
# Minimal Ebi parser


def MinimalEbiParser(lines, strict=True, selected_labels=None):
    """yield each (sequence as a str, a dict of header) from ebi record lines

    if strict (default), raise RecordError if a record lacks required labels.

    Warning: using the global required_labels.

    Line code   Content     Occurrence in an entry
    ID  Identification  Once; starts the entry
    AC  Accession number(s) Once or more
    PR  Project number(s) Once or more
    DT  Date    Three times
    DE  Description Once or more
    GN  Gene name(s)    Optional
    OS  Organism species    Once
    OG  Organelle   Optional
    OC  Organism classification Once or more
    RN  Reference number    Once or more
    RP  Reference position  Once or more
    RC  Reference comment(s)    Optional
    RX  Reference cross-reference(s)    Optional
    RG  Reference group Once or more (Optional if RA line)
    RA  Reference authors   Once or more (Optional if RG line)
    RT  Reference title Optional
    RL  Reference location  Once or more
    CC  Comments or notes   Optional
    DR  Database cross-references   Optional
    KW  Keywords    Optional
    FT  Feature table data  Optional
    SQ  Sequence header Once
    CO  Contig entries, currently skipped
    FH  Emtpy line or uninformative (for human readability)
    XX  Empty line
    AH  TPA and TSA records only, currently skipped
    AS  TPA and TSA records only, currently skipped
    (blanks)    Sequence data   Once or more
    //  Termination line    Once; ends the entry

    The two-character line-type code that begins each line is always followed
    by three blanks, so that the actual information begins with the sixth
    character. Information is not extended beyond character position 75 except
    for one exception: CC lines that contain the 'DATABASE' topic"""
    selected_labels = selected_labels or []
    exclude = b" \t\n\r/"
    strip_table = dict([(c, None) for c in exclude])

    for record in EbiFinder(lines):
        if strict and not record[0].startswith("ID"):
            raise RecordError("Record must begin with ID line")
        del record[-1]  # which must be //, ensured by Finder

        keyvalues = list(map(linecode_merging_maker, record))
        raw_dict = pairs_to_dict(keyvalues, "always_multi_value", all_keys=_parsers)

        if strict:
            for rlabel in required_labels:
                if rlabel not in raw_dict:
                    raise RecordError(f"The record lacks required label: {rlabel}")

        # no sequence found
        if "" not in raw_dict:
            continue

        sequence = raw_dict.pop("")  # which is the linecode for sequence
        sequence = "".join(sequence).translate(strip_table)

        if selected_labels:
            for key in list(raw_dict.keys()):
                if key not in selected_labels:
                    del raw_dict[key]

        header_dict = raw_dict
        yield sequence, header_dict


def linecode_merging_maker(line):
    """return merged linecode and the line.

    All valid reference linecodes merged into REF

    Warning: using global ref_parsers"""
    linecode = linecode_maker(line)[0]
    if linecode in ref_parsers:
        linecode = "REF"
    return linecode, line


#################################
# EbiParser


def parse_header(header_dict, strict=True):
    """Parses a dictionary of header lines"""
    return pairs_to_dict(
        list(header_dict.items()), "no_duplicated_key", handlers=_parsers
    )


_parsers = {
    "ID": id_parser,
    "AC": ac_parser,
    "DE": de_parser,
    "DT": dt_parser,
    "GN": gn_parser,
    "OC": oc_parser,
    "OS": os_parser,
    "OX": ox_parser,
    "OG": og_parser,
    "REF": refs_parser,
    "CC": cc_parser,
    "DR": dr_parser,
    "KW": kw_parser,
    "FT": ft_parser,
    "SQ": sq_parser,
    "PR": pr_parser,
    "XX": None,
    "FH": None,
    "CO": None,
    "AH": None,
    "AS": None,
    "": None,
    "//": None,
}


def EbiParser(
    lines,
    seq_constructor=Sequence,
    header_constructor=parse_header,
    strict=True,
    selected_labels=None,
):
    """Parser for the EBI data format.

    lines: input data (list of lines or file stream)
    seq_constructor: constructor function to construct sequence, 'Sequence' by
        default.
    header_constructor: function to process the header information. Default
        is 'parse_header'
    strict: whether an exception should be raised in case of a problem
        (strict=True) or whether the bad record should be skipped
        (strict=False).
    selected_labels: Labels from the original data format that you want
        returned. All the original header labels are used, except for
        REFERENCES, which is 'REF'.
    """
    selected_labels = selected_labels or []
    for sequence, header_dict in MinimalEbiParser(
        lines, strict=strict, selected_labels=selected_labels
    ):
        if seq_constructor:
            sequence = seq_constructor(sequence)
        try:
            header = header_constructor(header_dict, strict=strict)
        except (RecordError, FieldError, ValueError):
            if strict:
                #!! just raise is better than raise RecordError
                raise  # RecordError, str(e)
            else:
                continue

        yield sequence, header


if __name__ == "__main__":
    from getopt import GetoptError, getopt

    usage = """ Usage: python __.py [options] [source]

Options:
  -h, --help              show this help
  -d                      show debugging information while parsing

Examples:
"""
    try:
        opts, args = getopt(sys.argv[1:], "hd", ["help"])
    except GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(usage)
            sys.exit()
    if args:
        lines = open(args[0])
        print("Parsing the file")
        for i, rec in enumerate(EbiParser(lines, strict=True)):
            print(f"\r {i}: {rec[1]['ID']['EntryName']}", end=" ")
    else:
        lines = """\
ID   Q9U9C5_CAEEL   PRELIMINARY;      PRT;   218 AA.
AC   Q9U9C5;hdfksfsdfs;sdfsfsfs;
DT   01-MAY-2000 (TrEMBLrel. 13, Created)
DT   01-MAY-2000 (TrEMBLrel. 13, Last sequence update)
DT   13-SEP-2005 (TrEMBLrel. 31, Last annotation update)
DE   Basic salivary proline-rich protein 4 allele L (Salivary proline-rich
DE   protein Po) (Parotid o protein) [Contains: Peptide P-D (aa); BB (bb)
DE   (bbb)] (Fragment).
GN   name=nob-1; ORFNames=Y75B8A.2, Y75B8A.2B;
GN   and
GN   name=Jon99Ciii; Synonyms=SER2, SER5, Ser99Db; ORFNames=CG15519;
OS   Caenorhabditis elegans (aa) (bb).
OC   Eukaryota; Metazoa; Nematoda; Chromadorea; Rhabditida; Rhabditoidea;
OC   Rhabditidae; Peloderinae; Caenorhabditis.
OG   Plastid; Apicoplast.
OG   Plasmid R6-5, Plasmid IncFII R100 (NR1), and
OG   Plasmid IncFII R1-19 (R1 drd-19).
OX   NCBI_TaxID=6239;
RN   [1]
RP   NUCLEOTIDE SEQUENCE.
RC   STRAIN=N2;
RX   MEDLINE=20243724; PubMed=10781051; DOI=10.1073/pnas.97.9.4499;
RA   Van Auken K., Weaver D.C., Edgar L.G., Wood W.B.;
RT   "Caenorhabditis elegans embryonic axial patterning requires two
RT   recently discovered posterior-group Hox genes.";
RL   Proc. Natl. Acad. Sci. U.S.A. 97:4499-4503(2000).
RN   [2]
RP   NUCLEOTIDE SEQUENCE.
RC   STRAIN=N2;
RG   The mouse genome sequencing consortium;
RL   Submitted (JUL-1999) to the EMBL/GenBank/DDBJ databases.
CC   -!- SUBCELLULAR LOCATION: Nuclear (By similarity).
CC   -!- DATABASE: NAME=slkdfjAtlas Genet. Cytogenet. Oncol. Haematol.;
CC       WWW="http://www.infobiogen.fr/services/chromcancer/Genes/
CC   -!- DATABASE: NAME=Atlas Genet. Cytogenet. Oncol. Haematol.;
CC       WWW="http://www.infobiogen.fr/services/chromcancer/Genes/
CC       P53ID88.html".
CC   -!- INTERACTION:
CC       P51617:IRAK1; NbExp=1; IntAct=EBI-448466, EBI-358664;
CC       P51617:IRAK1; NbExp=1; IntAct=EBI-448472, EBI-358664;
CC   -!- ALTERNATIVE PRODUCTS:
CC       Event=Alternative splicing; Named isoforms=3;
CC         Comment=Additional isoforms seem to exist. Experimental
CC         confirmation may be lacking for some isoforms;
CC       name=1; Synonyms=AIRE-1;
CC         IsoId=O43918-1; Sequence=Displayed;
CC       name=2; Synonyms=AIRE-2;
CC         IsoId=O43918-2; Sequence=VSP_004089;
CC       name=3; Synonyms=AIRE-3;
CC         IsoId=O43918-3; Sequence=VSP_004089, VSP_004090;
CC   -!- BIOPHYSICOCHEMICAL PROPERTIES:
CC       Kinetic parameters:
CC         KM=98 uM for ATP;
CC         KM=688 uM for pyridoxal;
CC         Vmax=1.604 mmol/min/mg enzyme;
CC       pH dependence:
CC         Optimum pH is 6.0. Active from pH 4.5 to 10.5;
CC   -!- MASS SPECTROMETRY: MW=13822; METHOD=MALDI; RANGE=19-140 (P15522-
CC       2); NOTE=Ref.1.
CC   --------------------------------------------------------------------------
CC   This SWISS-PROT entry is copyright. It is produced through a collaboration
CC   removed.
CC   --------------------------------------------------------------------------
DR   EMBL; AF172090; AAD48874.1; -; mRNA.
DR   EMBL; AL033514; CAC70124.1; -; Genomic_DNA.
DR   HSSP; P02833; 9ANT.
KW   Complete proteome; DNA-binding; Developmental protein; Homeobox;
KW   Hypothetical protein; Nuclear protein.
FT   DNA_BIND    >102    292
FT   REGION        1     44       Transcription activation (acidic).
FT   CHAIN        23    611       Halfway protein.
FT                                /FTId=PRO_0000021413.
FT   VARIANT       1      7       unknown  (in a skin tumor).
FT                                /FTId=VAR_005851.
FT   VARIANT       7      7       D -> H (in a skin tumor).
FT                                /FTId=VAR_005851.
FT   CONFLICT    282    282       R -> Q (in Ref. 18).
FT   STRAND      103    103
FT   NON_TER     80     80        non_ter.
SQ   SEQUENCE   218 AA;  24367 MW;  F24AE5E8A102FAC6 CRC64;
     MISVMQQMIN NDSPEDSKES ITSVQQTPFF WPSAAAAIPS IQGESRSERE SETGSSPQLA
     PSSTGMVMPG TAGMYGFGPS RMPTANEFGM MMNPVYTDFY QNPLASTDIT IPTTAGSSAA
     TTPNAAMHLP WAISHDGKKK RQPYKKDQIS RLEYEYSVNQ YLTNKRRSEL SAQLMLDEKQ
     VKVWFQNRRM KDKKLRQRHS GPFPHGAPVT PCIERLIN
//
ID   Q9U9C5_TEST   PRELIMINARY;      PRT;   218 AA.
DT   ddd.
AC   Q9U9C5;hdfksfsdfs;sdfsfsfs;
SQ   SEQUENCE   218 AA;  24367 MW;  F24AE5E8A102FAC6 CRC64;
     MISVMQQMIN NDSPEDSKES ITSVQQTPFF WPSAAAAIPS IQGESRSERE
//
""".split(
            "\n"
        )
        pprint(list(EbiParser(lines, strict=False, selected_labels=[])))
