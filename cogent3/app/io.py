import os
import zipfile
import json
import numpy

from cogent3 import LoadSeqs
from cogent3.core.moltype import get_moltype
from cogent3.parse.sequence import PARSERS
from cogent3.format.alignment import FORMATTERS
from cogent3.core.alignment import ArrayAlignment, SequenceCollection
from cogent3.util.deserialise import deserialise_object
from cogent3.util.table import Table
from .data_store import (SingleReadDataStore, SKIP, RAISE,
                         OVERWRITE, IGNORE, ReadOnlyZippedDataStore,
                         ReadOnlyDirectoryDataStore,
                         WritableDirectoryDataStore, WritableZippedDataStore, )
from .composable import (ComposableSeq, ComposableAligned, Composable,
                         _checkpointable, ComposableTabular,
                         NotCompletedResult, )

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
    data_store = klass(base_path, suffix=suffix,
                       limit=limit, verbose=verbose)
    return data_store.members


def get_data_store(base_path, suffix, limit=None, verbose=False):
    """returns DataStore containing glob matches to suffix in base_path

    Parameters
    ----------
    base_path : str
        path to directory or zipped archive
    suffix : str
        suffix of filenames
    limit : int or None
        the number of matches to return
    Returns
    -------
    ReadOnlyDirectoryDataStore or ReadOnlyZippedDataStore
    """
    if not os.path.exists(base_path):
        raise ValueError(f"'{base_path}' does not exist")
    if not type(suffix) == str:
        raise ValueError(f'{suffix} is not a string')

    zipped = zipfile.is_zipfile(base_path)
    klass = ReadOnlyZippedDataStore if zipped else ReadOnlyDirectoryDataStore
    data_store = klass(base_path, suffix=suffix,
                       limit=limit, verbose=verbose)
    return data_store


class _seq_loader:
    def __init__(self):
        self.func = self.load

    def load(self, path):
        """returns alignment"""
        # if we get a seq object, we try getting abs_path from that now
        try:
            abs_path = path.info.source
        except AttributeError:
            abs_path = str(path)

        if type(path) == str:
            # we use a data store as it's read() handles compression
            path = SingleReadDataStore(path)[0]

        if isinstance(path, str):
            data = path.read().splitlines()
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

    def __init__(self, moltype=None, format='fasta'):
        """
        Parameters
        ----------
        moltype
            molecular type, string or instance
        format : str
            sequence file format
        """
        super(ComposableAligned, self).__init__(input_type=None,
                                                output_type='aligned')
        _seq_loader.__init__(self)
        self._formatted_params()
        if moltype:
            moltype = get_moltype(moltype)
        self.moltype = moltype
        self._parser = PARSERS[format.lower()]


class load_unaligned(ComposableSeq, _seq_loader):
    """loads sequences"""
    klass = SequenceCollection

    def __init__(self, moltype=None, format='fasta'):
        """
        Parameters
        ----------
        moltype
            molecular type, string or instance
        format : str
            sequence file format
        """
        super(ComposableSeq, self).__init__(input_type=None,
                                            output_type='sequences')
        _seq_loader.__init__(self)
        self._formatted_params()
        if moltype:
            moltype = get_moltype(moltype)
        self.moltype = moltype
        self._parser = PARSERS[format.lower()]


class load_tabular(ComposableTabular):
    def __init__(self, with_title=False, with_header=True, limit=None,
                 sep='\t', strict=True):
        """

        Parameters
        ----------
        with_title
            files have a title
        with_header
            files have a header
        limit
            number of records to read
        sep
            field delimiter
        strict
            all rows MUST have the same number of records
        """
        super(ComposableTabular, self).__init__(input_type=None,
                                                output_type='tabular')
        self._formatted_params()
        self._sep = sep
        self._with_title = with_title
        self._with_header = with_header
        self._limit = limit
        self.func = self.load
        self.strict = strict

    def _parse(self, data):
        title = header = None
        sep = self._sep
        strict = self.strict
        read = data.open()
        if self._with_title or self._with_header:
            for line in read:
                line = line.strip()
                if not line:
                    continue
                if self._with_title and title is None:
                    title = line
                elif self._with_header and header is None:
                    line = [e.strip() for e in line.split(sep)]
                    header = line
                    break
        num_records = None if header is None else len(header)
        rows = []
        for i, line in enumerate(read):
            if i == self._limit:
                break
            line = line.strip()
            line = [e.strip() for e in line.split(sep)]
            if num_records is None:
                num_records = len(line)
            if strict and len(line) != num_records:
                msg = (f'Inconsistent number of fields: {len(line)} '
                       '!= {num_records}')
                raise AssertionError(msg)
            rows.append(line)
        data.close()
        records = []
        for record in zip(*rows):
            record = numpy.array(record, dtype='O')
            try:
                record = record.astype(int)
            except ValueError:
                try:
                    record = record.astype(float)
                except ValueError:
                    pass
            records.append(record)
        records = numpy.array(records, dtype='O').T
        table = Table(header, rows=records, title=title)
        return table

    def load(self, path):
        if type(path) == str:
            # we use a data store as it's read() handles compression
            path = SingleReadDataStore(path)[0]

        try:
            result = self._parse(path)
        except Exception as err:
            result = NotCompletedResult('ERROR', self, err.args[0],
                                        source=str(path))

        return result


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
        self._formatter = FORMATTERS[format]

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

    def __init__(self):
        super(load_json, self).__init__(input_type=None,
                                        output_type=('result', 'serialisable'))
        self.func = self.read

    def read(self, path):
        """returns object deserialised from json at path"""
        if type(path) == str:
            path = SingleReadDataStore(path)[0]

        data = path.read()

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
