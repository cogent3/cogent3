#!/usr/bin/env python 
"""Application controller for stride."""

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

from cogent.struct.selection import einput
from cogent.parse.stride import stride_parser
from cogent.struct.annotation import xtradata
from cogent.app.util import CommandLineApplication, ResultPath
from cogent.app.parameters import FlagParameter, ValuedParameter
from cogent.format.pdb import PDBWriter


__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"


class Stride(CommandLineApplication):
    """Application controller for stride."""

    _parameters = {
                   # report hydrogen bonds
                   '-h':FlagParameter(Prefix='-', Name='h', Value=False),
                   # calculate contact order
                   '-k':FlagParameter(Prefix='-', Name='k', Value=False),
                   # print out the contact map only
                   '-w':FlagParameter(Prefix='-', Name='w', Value=False),
                   # write output as file
                   '-f':ValuedParameter(Prefix='-', Name='f', Delimiter=''),
                   # write output as molscript file
                   '-m':ValuedParameter(Prefix='-', Name='m', Delimiter=''),
    }

    _command = "stride"
    _input_handler = '_input_as_multiline_string'

    def _input_as_multiline_string(self, data):
        """This allows to feed entities to stride."""
        dummy_file = StringIO()
        PDBWriter(dummy_file, data)
        dummy_file.seek(0)
        data = super(Stride, self)._input_as_multiline_string(dummy_file.read())
        return data


    def _align_out_filename(self):

        if self.Parameters['-f'].isOn():
            aln_filename = self._absolute(str(self.Parameters['-f'].Value))
        else:
            aln_filename = None
        return aln_filename

    def _get_result_paths(self, data):

        result = {}
        if self.Parameters['-f'].isOn():
            out_name = self._align_out_filename()
            result['File'] = ResultPath(Path=out_name, IsWritten=True)
        return result


def stride_xtra(entities, **kwargs):
    """Runs the stride application and maps the result on the residue xtra
    dictionaries.
    """
    structures = einput(entities, 'S')
    if len(structures.values()) > 1:
        raise ValueError('Entities from multiple structures are not supported.')
    stride_app = Stride(**kwargs)
    result = stride_app(entities)['StdOut'].readlines()
    result = stride_parser(result)
    xtradata(result, structures.values()[0][(0,)])
    return result
