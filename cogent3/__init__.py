"""The most commonly used constructors are available from this toplevel module.
The rest are in the subpackages: core, draw, evolve, format, maths, parse
and phylo.
"""

import sys
import re
import pickle
import numpy
import os
import warnings

__author__ = ""
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley", "Rob Knight", "Peter Maxwell",
               "Jeremy Widmann", "Catherine Lozupone", "Matthew Wakefield",
               "Edward Lang", "Greg Caporaso", "Mike Robeson",
               "Micah Hamady", "Sandra Smit", "Zongzhi Liu",
               "Andrew Butterfield", "Amanda Birmingham", "Brett Easton",
               "Hua Ying", "Jason Carnes", "Raymond Sammut",
               "Helen Lindsay", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


if sys.version_info < (3, 6):
    PY_VERSION = ".".join([str(n) for n in sys.version_info])
    raise RuntimeError(
        "Python-3.6 or greater is required, Python-%s used." % PY_VERSION)

NUMPY_VERSION = re.split(r"[^\d]", numpy.__version__)
numpy_version_info = tuple([int(i) for i in NUMPY_VERSION if i.isdigit()])
if numpy_version_info < (1, 3):
    raise RuntimeError("Numpy-1.3 is required, %s found." % NUMPY_VERSION)

version = __version__
version_info = tuple([int(v) for v in version.split(".") if v.isdigit()])


warn_env = 'COGENT3_WARNINGS'

if warn_env in os.environ:
    warnings.simplefilter(os.environ[warn_env])


from cogent3.util.table import Table as _Table
from cogent3.util.misc import get_format_suffixes, open_
from cogent3.parse.table import load_delimited, autogen_reader
from cogent3.core.tree import TreeBuilder, TreeError
from cogent3.parse.tree_xml import parse_string as tree_xml_parse_string
from cogent3.parse.newick import parse_string as newick_parse_string
from cogent3.core.alignment import SequenceCollection
from cogent3.core.alignment import Alignment, ArrayAlignment
from cogent3.parse.sequence import FromFilenameParser
# note that moltype has to be imported last, because it sets the moltype in
# the objects created by the other modules.
from cogent3.core.moltype import (ASCII, DNA, RNA, PROTEIN, STANDARD_CODON,
    CodonAlphabet, get_moltype)


def Sequence(moltype=None, seq=None, name=None, filename=None, format=None):
    if seq is None:
        for (a_name, a_seq) in FromFilenameParser(filename, format):
            if seq is None:
                seq = a_seq
                if name is None:
                    name = a_name
            else:
                raise ValueError("Multiple sequences in '%s'" % filename)
    if moltype is not None:
        moltype = get_moltype(moltype)
        seq = moltype.make_seq(seq)
    elif not hasattr(seq, 'moltype'):
        seq = ASCII.make_seq(seq)
    if name is not None:
        seq.name = name
    return seq


def LoadSeqs(filename=None, format=None, data=None, moltype=None,
             name=None, aligned=True, label_to_name=None, parser_kw=None,
             constructor_kw=None, array_align=True, **kw):
    """Initialize an alignment or collection of sequences.

    Arguments:
    - filename: name of the sequence file
    - format: format of the sequence file
    - data: optional explicit provision of sequences
    - moltype: the moltype, eg DNA, PROTEIN, 'dna', 'protein'
    - aligned: set True if sequences are already aligned and have the same
      length. If False, a SequenceCollection instance is returned instead.
    - array_align: returns alignment as ArrayAlignment
    - label_to_name: function for converting original name into another
      name. Default behavior is to preserve the original FASTA label and
      comment.
      To remove all FASTA label comments, and pass in only the label, pass in:
            label_to_name=lambda x: x.split()[0]
      To look up names in a dict, pass in:
            label_to_name = lambda x: d.get(x, default_name)
      ...where d is a dict that's in scope, and default_name is what you want
      to assign any sequence that isn't in the dict.

    If format is None, will attempt to infer format from the filename
    suffix. If label_to_name is None, will attempt to infer correct
    conversion from the format.
    """
    
    constructor_kw = constructor_kw or {}
    parser_kw = parser_kw or {}

    if filename is None:
        assert data is not None
        assert format is None
        assert not kw, kw
    else:
        assert data is None, (filename, data)
        data = list(FromFilenameParser(filename, format, **parser_kw))

    if aligned:
        klass = ArrayAlignment if array_align else Alignment
    else:  # generic case: return SequenceCollection
        klass = SequenceCollection
    
    return klass(data=data, moltype=moltype, name=name,
                                  label_to_name=label_to_name, **constructor_kw)


def LoadTable(filename=None, sep=None, reader=None, header=None, rows=None,
              row_order=None, digits=4, space=4, title='', missing_data='',
              max_width=1e100, row_ids=None, legend='', column_templates=None,
              dtype=None, static_column_types=False, limit=None,
              data_frame=None, format="simple", **kwargs):
    """
    Arguments:
    - filename: path to file containing a pickled table
    - sep: the delimiting character between columns
    - reader: a parser for reading filename. This approach assumes the first
      row returned by the reader will be the header row.
    - static_column_types: if True, and reader is None, identifies columns
      with a numeric/bool data types from the first non-header row.
      This assumes all subsequent entries in that column are of the same type.
      Default is False.
    - header: column headings
    - rows: a 2D dict, list or tuple. If a dict, it must have column
      headings as top level keys, and common row labels as keys in each
      column.
    - row_order: the order in which rows will be pulled from the twoDdict
    - digits: floating point resolution
    - space: number of spaces between columns or a string
    - title: as implied
    - missing_data: character assigned if a row has no entry for a column
    - max_width: maximum column width for printing
    - row_ids: if True, the 0'th column is used as row identifiers and keys
      for slicing.
    - legend: table legend
    - column_templates: dict of column headings: string format templates
      or a function that will handle the formatting.
    - dtype: optional numpy array typecode.
    - limit: exits after this many lines. Only applied for non pickled data
      file types.
    - data_frame: a pandas DataFrame, supersedes header/rows
    - format: output format when using str(Table)
    """
    sep = sep or kwargs.pop('delimiter', None)
    if filename is not None:
        file_format, compress_format = get_format_suffixes(filename)

    if filename is not None and not (reader or static_column_types):
        if file_format == 'pickle':
            f = open_(filename, mode='rb')
            loaded_table = pickle.load(f)
            f.close()
            return _Table(**loaded_table)
        elif file_format == 'csv':
            sep = sep or ','
        elif file_format == 'tsv':
            sep = sep or '\t'
        
        header, rows, loaded_title, legend = load_delimited(filename,
                                                            delimiter=sep, limit=limit, **kwargs)
        title = title or loaded_title
    elif filename and (reader or static_column_types):
        f = open_(filename, newline=None)
        if not reader:
            if file_format == 'csv':
                sep = sep or ','
            elif file_format == 'tsv':
                sep = sep or '\t'
            elif not sep:
                raise ValueError("static_column_types option requires a value "
                                 "for sep")
            
            reader = autogen_reader(f, sep, limit=limit,
                                    with_title=kwargs.get('with_title', False))

        rows = [row for row in reader(f)]
        f.close()
        header = rows.pop(0)

    table = _Table(header=header, rows=rows, digits=digits, row_order=row_order,
                   title=title,
                   dtype=dtype, column_templates=column_templates, space=space,
                   missing_data=missing_data, max_width=max_width, row_ids=row_ids,
                   legend=legend,
                   data_frame=data_frame, format=format)

    return table


def LoadTree(filename=None, treestring=None, tip_names=None, format=None,
             underscore_unmunge=False):
    """Constructor for tree.

    Arguments, use only one of:
        - filename: a file containing a newick or xml formatted tree.
        - treestring: a newick or xml formatted tree string.
        - tip_names: a list of tip names.

    Note: underscore_unmunging is turned off by default, although it is part
    of the Newick format. Set underscore_unmunge to True to replace underscores
    with spaces in all names read.
    """

    if filename:
        assert not (treestring or tip_names)
        with open_(filename) as tfile:
            treestring = tfile.read()
        if format is None and filename.endswith('.xml'):
            format = "xml"
    if treestring:
        assert not tip_names
        if format is None and treestring.startswith('<'):
            format = "xml"
        if format == "xml":
            parser = tree_xml_parse_string
        else:
            parser = newick_parse_string
        tree_builder = TreeBuilder().create_edge
        # FIXME: More general strategy for underscore_unmunge
        if parser is newick_parse_string:
            tree = parser(treestring, tree_builder,
                          underscore_unmunge=underscore_unmunge)
        else:
            tree = parser(treestring, tree_builder)
        if not tree.name_loaded:
            tree.name = 'root'
    elif tip_names:
        tree_builder = TreeBuilder().create_edge
        tips = [tree_builder([], tip_name, {}) for tip_name in tip_names]
        tree = tree_builder(tips, 'root', {})
    else:
        raise TreeError('filename or treestring not specified')
    return tree
