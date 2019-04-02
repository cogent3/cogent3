import os
import zipfile
import json

from cogent3 import LoadSeqs
from cogent3.core.moltype import get_moltype
from cogent3.parse.sequence import PARSERS
from cogent3.format.alignment import WRITERS
from cogent3.core.alignment import ArrayAlignment, SequenceCollection
from cogent3.util.deserialise import deserialise_object
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


def findall(base_path, suffix='fa', limit=None, verbose=False):
    """returns glob match to suffix, path is relative to base_path

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
    data_store = klass(base_path, suffix=suffix, limit=limit, verbose=verbose)
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
            klass = (ReadOnlyZippedDataStore
                     if zipped else ReadOnlyDirectoryDataStore)
            self.data_store = klass(data_path, suffix=suffix)
        else:
            self.data_store = None

        self.func = self.load

    def load(self, path):
        """returns alignment"""
        # if we get a seq object, we try getting abs_path from that now
        try:
            abs_path = path.info.source
        except AttributeError:
            abs_path = ''

        if type(path) == str:
            if self.data_store:
                data_store = self.data_store
            else:
                data_store = SingleReadDataStore(path)

            abs_path = data_store.get_absolute_identifier(path)
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


class write_seqs(_checkpointable):
    def __init__(self, data_path, format='fasta', suffix='fa',
                 name_callback=None, create=False, if_exists=SKIP):
        """
        Parameters
        ----------
        data_path
            path to write output, if ends with .zip will be a compressed zip
            archive
        format : str
            sequence file format
        suffix : str
            filename suffix for output
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
                                                      'aligned',
                                                      'identifier'),
                                         data_path=data_path,
                                         name_callback=name_callback,
                                         create=create, if_exists=if_exists,
                                         suffix=suffix)
        self._formatted_params()
        if format != 'fasta':
            # todo refactor alignment formatters
            raise NotImplementedError('must use fasta for now')

        self._format = format
        self._formatter = WRITERS[format]

    def _set_checkpoint_loader(self):
        loader = {'sequences': load_unaligned}.get(self._out._type,
                                                   load_aligned)
        loader = loader(format=self._format)
        self._load_checkpoint = loader

    def write(self, data):
        identifier = self._make_output_identifier(data)
        data.info.stored = self.data_store.write(identifier, data.to_fasta())
        return identifier


class load_json(Composable):
    _type = 'output'

    def __init__(self, data_path=None, suffix='json'):
        super(load_json, self).__init__(input_type=None,
                                        output_type=('result', 'serialisable'))
        """
        Parameters
        ----------
        data_path : str
            path to a directory or zip archive
        """
        if data_path:
            zipped = zipfile.is_zipfile(data_path)
            klass = (ReadOnlyZippedDataStore
                     if zipped else ReadOnlyDirectoryDataStore)
            self.data_store = klass(data_path, suffix=suffix)
        else:
            self.data_store = None

        self.func = self.read

    def read(self, path):
        """returns alignment"""
        if self.data_store:
            data_store = self.data_store
        else:
            data_store = SingleReadDataStore(path)
        if type(path) == str:
            data = data_store.read(path)
        else:
            raise NotImplementedError

        return deserialise_object(data)


class write_json(_checkpointable):
    _type = 'output'

    def __init__(self, data_path, name_callback=None, create=False,
                 if_exists=SKIP, suffix='json'):
        super(write_json, self).__init__(input_type='serialisable',
                                         output_type=('identifier',
                                                      'serialisable'),
                                         data_path=data_path,
                                         name_callback=name_callback,
                                         create=create, if_exists=if_exists,
                                         suffix=suffix)
        self.func = self.write

    def _set_checkpoint_loader(self):
        self._load_checkpoint = load_json(self.data_store.source,
                                          suffix=self.data_store.suffix)

    def write(self, data):
        identifier = self._make_output_identifier(data)
        out = data.to_json()
        stored = self.data_store.write(identifier, out)
        try:
            data.info.stored = stored
        except AttributeError:
            data.stored = stored
        return identifier
