"""The most commonly used constructors are available from this toplevel module.
The rest are in the subpackages: phylo, evolve, maths, draw, parse and format."""

import sys, os, logging, re, cPickle
import numpy

__author__ = ""
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Gavin Huttley", "Rob Knight", "Peter Maxwell",
                    "Jeremy Widmann", "Catherine Lozupone", "Matthew Wakefield",
                    "Edward Lang", "Greg Caporaso", "Mike Robeson",
                    "Micah Hamady", "Sandra Smit", "Zongzhi Liu", 
                    "Andrew Butterfield", "Amanda Birmingham", "Brett Easton",
                    "Hua Ying", "Jason Carnes", "Raymond Sammut", 
                    "Helen Lindsay", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

#SUPPORT2425
if sys.version_info < (2, 4):
    py_version = ".".join([str(n) for n in sys.version_info])
    raise RuntimeError("Python-2.4 or 2.5 is required, Python-%s used." % py_version)

numpy_version = re.split("[^\d]", numpy.__version__)
numpy_version_info = tuple([int(i) for i in numpy_version if i.isdigit()])
if numpy_version_info < (1, 0):
    raise RuntimeError("Numpy-1.0 is required, %s found." % numpy_version)

LOG = logging.getLogger('cogent')
LOG.addHandler(logging.StreamHandler())

if 'COGENT_LOG_LEVEL' in os.environ:
    level = os.environ['COGENT_LOG_LEVEL']
    valid_log_levels = ['DEBUG', 'INFO', 'WARNING', 'ERROR']
    assert level in valid_log_levels, valid_log_levels
    logging.basicConfig(level=getattr(logging, level))

if sys.version_info > (2, 5, 3):
    LOG.warning("cogent.align.methods.ACL fails with Python-2.6")

version = __version__
version_info = tuple([int(v) for v in __version__.split(".")])

from cogent.util.table import Table as _Table
from cogent.parse.table import load_delimited
from cogent.core.tree import LoadTree
from cogent.core.alignment import SequenceCollection

from cogent.core.alignment import Alignment
from cogent.parse.sequence import FromFilenameParser
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
            name=None, aligned=True, name_conversion_f=None,
            parser_kw={}, constructor_kw={}, **kw):
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
    - name_conversion_f: function for converting original name into another
      name. Default behavior is to split on first whitespace, i.e. to provide
      FASTA labels. 
      To force names to be left alone, pass in: 
            name_conversion_f=lambda x: x
      To look up names in a dict, pass in:
            name_conversion_f = lambda x: d.get(x, default_name)
      ...where d is a dict that's in scope, and default_name is what you want
      to assign any sequence that isn't in the dict.
    
    If format is None, will attempt to infer format from the filename
    suffix. If name_conversion_f is None, will attempt to infer correct
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
                name_conversion_f=name_conversion_f, **constructor_kw)
        else:   #was not callable, but wasn't False
            return Alignment(data=data, MolType=moltype, Name=name, 
                name_conversion_f=name_conversion_f, **constructor_kw)
    else:   #generic case: return SequenceCollection
        return SequenceCollection(data, MolType=moltype, Name=name, 
            name_conversion_f=name_conversion_f, **constructor_kw)

def LoadTable(filename=None, sep=',', reader=None, header=None, rows=None, 
            row_order=None, digits=4, space=4, title='', missing_data='',
            max_width = 1e100, row_ids=False, legend='', column_templates=None,
            dtype=None,  **kwargs):
    """
    Arguments:
    - filename: path to file containing a pickled table
    - sep: the delimiting character between columns
    - reader: a parser for reading filename. This approach assumes the first
      row returned by the reader will be the header row.
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
    """
    # 
    if filename is not None and reader is None:
        if filename[filename.rfind(".")+1:] == 'pickle':
            f = file(filename, 'U')
            loaded_table = cPickle.load(f)
            f.close()
            return _Table(**loaded_table)
        
        sep = sep or kwargs.pop('delimiter', None)
        header, rows, loaded_title, legend = load_delimited(filename,
                                        delimiter = sep, **kwargs)
        title = title or loaded_title
    elif filename and reader:
        f = file(filename, "r")
        rows = [row for row in reader(f)]
        f.close()
        header = rows.pop(0)
    
    table = _Table(header=header, rows=rows, digits=digits, row_order=row_order,
                title=title,
                dtype=dtype, column_templates=column_templates, space=space,
                missing_data=missing_data, max_width=max_width, row_ids=row_ids,
                legend=legend)
    
    return table
