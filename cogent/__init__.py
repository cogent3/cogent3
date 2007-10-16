"""The most commonly used constructors are available from this toplevel module.
The rest are in the subpackages: phylo, evolve, maths, draw, parse and format."""

import sys, os, logging, re
import numpy

__author__ = ""
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Gavin Huttley", "Rob Knight", "Peter Maxwell",
                    "Jeremy Widmann", "Catherine Lozupone", "Matthew Wakefield",
                    "Edward Lang", "Greg Caporaso", "Mike Robeson",
                    "Micah Hamady", "Sandra Smit", "Zongzhi Liu", 
                    "Andrew Butterfield", "Amanda Birmingham", "Brett Easton",
                    "Hua Ying", "Jason Carnes", "Raymond Sammut", 
                    "Helen Lindsay"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

#SUPPORT2425
if sys.version_info < (2, 4):
    py_version = ".".join([str(n) for n in sys.version_info])
    raise RuntimeError("Python-2.4 is required, Python-%s used." % py_version)

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

version = __version__
version_info = tuple([int(v) for v in __version__.split(".")])

from cogent.util.table import Table
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
