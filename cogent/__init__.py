"""The most commonly used constructors are available from this toplevel module.
The rest are in the subpackages: phylo, evolve, maths, draw, parse and format."""

import sys, os, re, cPickle
import numpy

__author__ = ""
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Rob Knight", "Peter Maxwell",
                    "Jeremy Widmann", "Catherine Lozupone", "Matthew Wakefield",
                    "Edward Lang", "Greg Caporaso", "Mike Robeson",
                    "Micah Hamady", "Sandra Smit", "Zongzhi Liu",
                    "Andrew Butterfield", "Amanda Birmingham", "Brett Easton",
                    "Hua Ying", "Jason Carnes", "Raymond Sammut",
                    "Helen Lindsay", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

#SUPPORT2425
if sys.version_info < (2, 6):
    py_version = ".".join([str(n) for n in sys.version_info])
    raise RuntimeError("Python-2.6 or greater is required, Python-%s used." % py_version)

numpy_version = re.split("[^\d]", numpy.__version__)
numpy_version_info = tuple([int(i) for i in numpy_version if i.isdigit()])
if numpy_version_info < (1, 3):
    raise RuntimeError("Numpy-1.3 is required, %s found." % numpy_version)

version = __version__
version_info = tuple([int(v) for v in version.split(".") if v.isdigit()])


from cogent.util.table import Table as _Table
from cogent.parse.table import load_delimited, autogen_reader
from cogent.core.tree import TreeBuilder, TreeError
from cogent.parse.tree_xml import parse_string as tree_xml_parse_string
from cogent.parse.newick import parse_string as newick_parse_string
from cogent.core.alignment import SequenceCollection
from cogent.core.alignment import Alignment
from cogent.parse.sequence import FromFilenameParser
from cogent.parse.structure import FromFilenameStructureParser
#note that moltype has to be imported last, because it sets the moltype in
#the objects created by the other modules.
from cogent.core.moltype import ASCII, DNA, RNA, PROTEIN, STANDARD_CODON, \
        CodonAlphabet

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
        seq = moltype.makeSequence(seq)
    elif not hasattr(seq, 'MolType'):
        seq = ASCII.makeSequence(seq)
    if name is not None:
        seq.Name = name
    return seq

def LoadSeqs(filename=None, format=None, data=None, moltype=None,
            name=None, aligned=True, label_to_name=None, parser_kw={},
            constructor_kw={}, **kw):
    """Initialize an alignment or collection of sequences.
    
    Arguments:
    - filename: name of the sequence file
    - format: format of the sequence file
    - data: optional explicit provision of sequences
    - moltype: the MolType, eg DNA, PROTEIN
    - aligned: set True if sequences are already aligned and have the same
      length, results in an Alignment object. If False, a SequenceCollection
      instance is returned instead. If callable, will use as a constructor
      (e.g. can pass in DenseAlignment or CodonAlignment).
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
    
    if filename is None:
        assert data is not None
        assert format is None
        assert not kw, kw
    else:
        assert data is None, (filename, data)
        data = list(FromFilenameParser(filename, format, **parser_kw))

    # the following is a temp hack until we have the load API sorted out.
    if aligned: #if callable, call it -- expect either f(data) or bool
        if hasattr(aligned, '__call__'):
            return aligned(data=data, MolType=moltype, Name=name,
                label_to_name=label_to_name, **constructor_kw)
        else:   #was not callable, but wasn't False
            return Alignment(data=data, MolType=moltype, Name=name,
                label_to_name=label_to_name, **constructor_kw)
    else:   #generic case: return SequenceCollection
        return SequenceCollection(data, MolType=moltype, Name=name,
            label_to_name=label_to_name, **constructor_kw)

def LoadStructure(filename, format=None, parser_kw={}):
    """Initialize a Structure from data contained in filename.
    Arguments:
        - filename: name of the filename to create structure from.
        - format: the optional file format extension.
        - parser_kw: optional keyword arguments for the parser."""
    # currently there is no support for string-input
    assert filename is not None, 'No filename given.'
    return FromFilenameStructureParser(filename, format, **parser_kw)


def LoadTable(filename=None, sep=',', reader=None, header=None, rows=None,
            row_order=None, digits=4, space=4, title='', missing_data='',
            max_width = 1e100, row_ids=False, legend='', column_templates=None,
            dtype=None, static_column_types=False, limit=None, **kwargs):
    """
    Arguments:
    - filename: path to file containing a pickled table
    - sep: the delimiting character between columns
    - reader: a parser for reading filename. This approach assumes the first
      row returned by the reader will be the header row.
    - static_column_types: if True, and reader is None, identifies columns
      with a numeric data type (int, float) from the first non-header row.
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
    """
    # 
    if filename is not None and not (reader or static_column_types):
        if filename[filename.rfind(".")+1:] == 'pickle':
            f = file(filename, 'U')
            loaded_table = cPickle.load(f)
            f.close()
            return _Table(**loaded_table)

        sep = sep or kwargs.pop('delimiter', None)
        header, rows, loaded_title, legend = load_delimited(filename,
                                    delimiter = sep, limit=limit, **kwargs)
        title = title or loaded_title
    elif filename and (reader or static_column_types):
        f = file(filename, "r")
        if not reader:
            reader = autogen_reader(f, sep, limit=limit,
                        with_title=kwargs.get('with_title', False))
        
        rows = [row for row in reader(f)]
        f.close()
        header = rows.pop(0)

    table = _Table(header=header, rows=rows, digits=digits, row_order=row_order,
                title=title,
                dtype=dtype, column_templates=column_templates, space=space,
                missing_data=missing_data, max_width=max_width, row_ids=row_ids,
                legend=legend)

    return table

def LoadTree(filename=None, treestring=None, tip_names=None, format=None, \
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
        treestring = open(filename).read()
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
        tree_builder = TreeBuilder().createEdge
        #FIXME: More general strategy for underscore_unmunge
        if parser is newick_parse_string:
            tree = parser(treestring, tree_builder, \
                    underscore_unmunge=underscore_unmunge)
        else:
            tree = parser(treestring, tree_builder)
        if not tree.NameLoaded:
            tree.Name = 'root'
    elif tip_names:
        tree_builder = TreeBuilder().createEdge
        tips = [tree_builder([], tip_name, {}) for tip_name in tip_names]
        tree = tree_builder(tips, 'root', {})
    else:
        raise TreeError, 'filename or treestring not specified'
    return tree

