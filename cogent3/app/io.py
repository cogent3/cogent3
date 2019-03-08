import os
import zipfile
import json

from cogent3 import LoadSeqs
from cogent3.core.moltype import get_moltype
from cogent3.parse.sequence import PARSERS
from cogent3.format.alignment import WRITERS
from cogent3.core.alignment import ArrayAlignment, SequenceCollection
from cogent3.util.deserialise import deserialise_object
from cogent3.util.misc import open_
from .data_store import (SingleReadDataStore, SKIP, RAISE,
                         OVERWRITE, IGNORE, ReadOnlyZippedDataStore,
                         ReadOnlyDirectoryDataStore,
                         WritableDirectoryDataStore, WritableZippedDataStore, )
from .composable import (ComposableSeq, ComposableAligned, Composable,
                         _checkpointable, )

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def abspath(path):
    """returns absolute path, expanding for user"""
    path = os.path.abspath(os.path.expanduser(path))
    return path


def findall(base_path, suffix='fa', limit=None):
    """returns glob match to suffix

    Parameters
    ----------
    base_path : str
        path to directory or zipped archive
    suffix : str
        suffix of filenames
    limit : int or None
        the number of matches to return
    """
    if not os.path.exists(base_path):
        raise ValueError(f"'{base_path}' does not exist")

    zipped = zipfile.is_zipfile(base_path)
    klass = ReadOnlyZippedDataStore if zipped else ReadOnlyDirectoryDataStore
    data_store = klass(base_path, suffix=suffix, limit=limit)
    return data_store.members


class _seq_loader:
    def __init__(self, data_path, suffix='fasta'):
        """
        Parameters
        ----------
        data_path : str
            path to a directory or zip archive
        """
        if data_path:
            zipped = zipfile.is_zipfile(data_path)
            klass = (
                ReadOnlyZippedDataStore if zipped else ReadOnlyDirectoryDataStore)
            self.data_store = klass(data_path, suffix=suffix)
        else:
            self.data_store = None

        self.func = self.load

    def load(self, path):
        """returns alignment"""
        if self.data_store:
            data_store = self.data_store
        else:
            data_store = SingleReadDataStore(path)
        abs_path = data_store.get_absolute_identifier(path)
        if type(path) == str:
            data = data_store.read(path).splitlines()
            data = dict(record for record in self._parser(data))
            seqs = self.klass(data=data, moltype=self.moltype)
            seqs.info.source = abs_path
        elif not isinstance(path, SequenceCollection):
            seqs = LoadSeqs(data=path, moltype=self.moltype,
                            aligned=self.aligned)
        else:
            seqs = path  # it is a SequenceCollection

        if self._output_type == {'sequences'}:
            seqs = seqs.degap()
            seqs.info.source = abs_path

        return seqs


class load_aligned(_seq_loader, ComposableAligned):
    """loads sequences"""
    klass = ArrayAlignment

    def __init__(self, data_path=None, moltype=None, format='fasta',
                 suffix='fa'):
        """
        Parameters
        ----------
        moltype
            molecular type, string or instance
        format : str
            sequence file format
        data_path : str
            path to a directory or zip archive
        """
        super(ComposableAligned, self).__init__(input_type=None,
                                                output_type='aligned')
        _seq_loader.__init__(self, data_path, suffix=suffix)
        self._formatted_params()
        if moltype:
            moltype = get_moltype(moltype)
        self.moltype = moltype
        self._parser = PARSERS[format.lower()]


class load_unaligned(ComposableSeq, _seq_loader):
    """loads sequences"""
    klass = SequenceCollection

    def __init__(self, data_path=None, moltype=None, format='fasta',
                 suffix='fa'):
        """
        Parameters
        ----------
        data_path : str
            path to a directory or zip archive
        moltype
            molecular type, string or instance
        format : str
            sequence file format
        """
        super(ComposableSeq, self).__init__(input_type=None,
                                            output_type='sequences')
        _seq_loader.__init__(self, data_path, suffix=suffix)
        self._formatted_params()
        if moltype:
            moltype = get_moltype(moltype)
        self.moltype = moltype
        self._parser = PARSERS[format.lower()]


class write_seqs(ComposableSeq):
    def __init__(self, data_path, format='fasta', name_callback=None,
                 create=False, if_exists='skip'):
        """
        Parameters
        ----------
        data_path
            path to write output, if ends with .zip will be a compressed zip
            archive
        format : str
            sequence file format
        name_callback
            function that takes the data object and returns a base
            file name
        create : bool
            whether to create the output directory
        if_exists : str
            behaviour if output exists. Either 'skip', 'raise' (raises an
            exception), 'overwrite'
        """
        super(write_seqs, self).__init__(input_type=('sequences',
                                                     'aligned'),
                                         output_type=('sequences',
                                                      'aligned'))
        self._formatted_params()

        if_exists = if_exists.lower()
        assert if_exists in (SKIP, IGNORE,
                             RAISE, OVERWRITE), 'invalid value for if_exists'
        self._if_exists = if_exists

        zipped = data_path.endswith('.zip')
        klass = WritableZippedDataStore if zipped else WritableDirectoryDataStore
        self.data_store = klass(data_path, suffix=format, create=create,
                                if_exists=if_exists)
        if format != 'fasta':
            # todo refactor alignment formatters
            raise NotImplementedError('must use fasta for now')

        self._format = format
        self._callback = name_callback
        self.func = self.write
        self._formatter = WRITERS[format]
        self._format = format
        self._checkpointable = True

    def _set_checkpoint_loader(self):
        loader = {'sequences': load_unaligned}.get(self._out._type,
                                                   load_aligned)
        loader = loader(format=self._format)
        self._load_checkpoint = loader

    def _make_output_identifier(self, data):
        if self._callback:
            data = self._callback(data)

        identifier = self.data_store.make_absolute_identifier(data)
        return identifier

    def job_done(self, data):
        identifier = self._make_output_identifier(data)
        exists = identifier in self.data_store
        if exists and self._if_exists == 'raise':
            msg = "'%s' already exists" % identifier
            raise RuntimeError(msg)

        if self._if_exists == OVERWRITE:
            exists = False
        return exists

    def write(self, data):
        identifier = self._make_output_identifier(data)
        data.info.stored = self.data_store.write(identifier, data.to_fasta())
        return identifier

