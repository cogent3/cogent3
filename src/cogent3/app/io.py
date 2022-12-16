import json
import os
import pathlib
import zipfile

from typing import Union

import numpy

from cogent3.app.io_new import _datastore_reader_map, register_datastore_reader
from cogent3.core.alignment import ArrayAlignment, SequenceCollection
from cogent3.core.moltype import get_moltype
from cogent3.core.profile import (
    make_motif_counts_from_tabular,
    make_motif_freqs_from_tabular,
    make_pssm_from_tabular,
)
from cogent3.evolve.fast_distance import DistanceMatrix
from cogent3.format.alignment import FORMATTERS
from cogent3.parse.sequence import PARSERS
from cogent3.util.deserialise import deserialise_object
from cogent3.util.table import Table

from .composable import (
    LOADER,
    WRITER,
    NotCompleted,
    _checkpointable,
    define_app,
)
from .data_store import (
    SKIP,
    SingleReadDataStore,
    WritableTinyDbDataStore,
    get_data_source,
    load_record_from_json,
    make_record_for_json,
)
from .typing import (
    AlignedSeqsType,
    IdentifierType,
    SeqsCollectionType,
    SerialisableType,
    TabularType,
    UnalignedSeqsType,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.10.31a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def findall(base_path, suffix="fa", limit=None, verbose=False):
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
    from cogent3.util.warning import discontinued

    discontinued(
        "function", "findall", "2023.3", "use the new cogent3.open_data_store()"
    )

    if not os.path.exists(base_path):
        raise ValueError(f"'{base_path}' does not exist")

    zipped = zipfile.is_zipfile(base_path)
    klass = _datastore_reader_map.get(".zip" if zipped else None)
    data_store = klass(base_path, suffix=suffix, limit=limit, verbose=verbose)
    return data_store.members


def get_data_store(
    base_path: Union[str, pathlib.Path], suffix=None, limit=None, verbose=False
):  # pragma: no cover
    """DEPRECATED, use top level open_data_store"""
    from cogent3 import open_data_store
    from cogent3.util.warning import deprecated

    deprecated(
        "function",
        "get_data_store",
        "cogent3.open_data_store",
        "2023.3",
        "renamed and now a top-level import",
        stack_level=1,
    )

    return open_data_store(base_path=base_path, suffix=suffix, limit=limit)


def _load_seqs(path, klass, parser, moltype):  # pragma: no cover
    abs_path = str(path)
    if not hasattr(path, "read"):
        # we use a DataStoreMember as it's read() handles zipped compression
        path = SingleReadDataStore(str(path))[0]
    data = path.read().splitlines()
    data = dict(iter(parser(data)))
    seqs = klass(data=data, moltype=moltype)
    seqs.info.source = abs_path
    return seqs


@define_app(app_type=LOADER)
class load_aligned:  # pragma: no cover
    """Loads aligned sequences. Returns an Alignment object."""

    klass = ArrayAlignment

    def __init__(self, moltype=None, format="fasta"):
        """
        Parameters
        ----------
        moltype
            molecular type, string or instance
        format : str
            sequence file format
        """
        from cogent3.util.warning import deprecated

        deprecated(
            "class",
            "load_aligned",
            "use cogent3.get_app('load_aligned')",
            "2023.3",
            "use cogent3.get_app()",
            stack_level=1,
        )

        if moltype:
            moltype = get_moltype(moltype)
        self.moltype = moltype
        self._parser = PARSERS[format.lower()]

    T = Union[SerialisableType, AlignedSeqsType]

    def main(self, path: IdentifierType) -> T:
        """returns alignment"""
        return _load_seqs(path, self.klass, self._parser, self.moltype)


@define_app(app_type=LOADER)
class load_unaligned:  # pragma: no cover
    """Loads unaligned sequences. Returns a SequenceCollection."""

    klass = SequenceCollection

    def __init__(self, moltype=None, format="fasta"):
        """
        Parameters
        ----------
        moltype
            molecular type, string or instance
        format : str
            sequence file format
        """
        from cogent3.util.warning import deprecated

        deprecated(
            "class",
            "load_unaligned",
            "use cogent3.get_app('load_unaligned')",
            "2023.3",
            "use cogent3.get_app()",
            stack_level=1,
        )

        if moltype:
            moltype = get_moltype(moltype)
        self.moltype = moltype
        self._parser = PARSERS[format.lower()]

    T = Union[SerialisableType, UnalignedSeqsType]

    def main(self, path: IdentifierType) -> T:
        """returns sequence collection"""
        seqs = _load_seqs(path, self.klass, self._parser, self.moltype)
        return seqs.degap()


@define_app(app_type=LOADER)
class load_tabular:  # pragma: no coverv
    """Loads delimited data. Returns a Table."""

    def __init__(
        self,
        with_title=False,
        with_header=True,
        limit=None,
        sep="\t",
        strict=True,
        as_type="table",
    ):
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
        from cogent3.util.warning import deprecated

        deprecated(
            "class",
            "load_tabular",
            "use cogent3.get_app('load_tabular')",
            "2023.3",
            "use cogent3.get_app()",
            stack_level=1,
        )

        self._sep = sep
        self._with_title = with_title
        self._with_header = with_header
        self._limit = limit
        self.strict = strict
        self.as_type = as_type

    def _parse(self, data):
        """returns header, records, title"""
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
                read.close()
                raise AssertionError(
                    f"Inconsistent number of fields: {len(line)} != {num_records}"
                )
            rows.append(line)
        data.close()
        records = []
        for record in zip(*rows):
            record = numpy.array(record, dtype="O")
            try:
                record = record.astype(int)
            except ValueError:
                try:
                    record = record.astype(float)
                except ValueError:
                    pass
            records.append(record)
        records = numpy.array(records, dtype="O").T
        return header, records, title

    def main(self, path: IdentifierType) -> TabularType:
        if type(path) == str:
            # we use a data store as it's read() handles compression
            path = SingleReadDataStore(path)[0]

        try:
            header, data, title = self._parse(path)
        except Exception as err:
            return NotCompleted("ERROR", self, err.args[0], source=str(path))

        if self.as_type == "table":
            return Table(header=header, data=data, title=title)

        assert data.shape[1] == 3, "Invalid tabular data"

        if self.as_type == "distances":
            # records is of the form [ [dim-1, dim-2, value] for entries in DistanceMatrix ]
            return DistanceMatrix({(e[0], e[1]): e[2] for e in data})

        if self.as_type == "motif_counts":
            return make_motif_counts_from_tabular(data)
        if self.as_type == "motif_freqs":
            return make_motif_freqs_from_tabular(data)
        if self.as_type == "pssm":
            return make_pssm_from_tabular(data)

        return None


@define_app(app_type=WRITER)
class write_tabular(_checkpointable):  # pragma: no cover
    """writes tabular data"""

    def __init__(
        self, data_path, format="tsv", name_callback=None, create=False, if_exists=SKIP
    ):
        """
        Parameters
        ----------
        data_path
            path to write output, if ends with .zip will be a compressed zip
            archive
        format : str
            one of 'tsv', 'csv', 'tex', 'md'
        name_callback
            function that takes the data object and returns a base
            file name
        create : bool
            whether to create the output directory
        if_exists : str
            behaviour if output exists. Either 'skip', 'raise' (raises an
            exception), 'overwrite'
        """
        from cogent3.util.warning import deprecated

        deprecated(
            "class",
            "write_tabular",
            "use cogent3.get_app('write_tabular')",
            "2023.3",
            "use cogent3.get_app()",
            stack_level=1,
        )

        super().__init__(
            data_path=data_path,
            name_callback=name_callback,
            create=create,
            if_exists=if_exists,
            suffix=format,
        )
        self._format = format

    def main(self, data: TabularType, identifier=None) -> IdentifierType:
        if isinstance(data, NotCompleted):
            return self.data_store.write_incomplete(identifier, data)

        if not data:
            msg = f"{self.__class__.__name__!r} does not support writing {data!r}"
            raise NotImplementedError(msg)

        if identifier is None:
            identifier = self._make_output_identifier(data)

        output = data.to_string(format=self._format)
        return self.data_store.write(identifier, output)


@define_app(app_type=WRITER)
class write_seqs(_checkpointable):  # pragma: no cover
    """Writes sequences to text files in standard format."""

    def __init__(
        self,
        data_path,
        format="fasta",
        suffix="fa",
        name_callback=None,
        create=False,
        if_exists=SKIP,
    ):
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
        from cogent3.util.warning import deprecated

        deprecated(
            "class",
            "write_seqs",
            "use cogent3.get_app('write_seqs')",
            "2023.3",
            "use cogent3.get_app()",
            stack_level=1,
        )

        super().__init__(
            data_path=data_path,
            name_callback=name_callback,
            create=create,
            if_exists=if_exists,
            suffix=suffix,
        )
        self._format = format
        self._formatter = FORMATTERS[format]

    def _set_checkpoint_loader(self):
        loader = {"sequences": load_unaligned}.get(self._out._type, load_aligned)
        loader = loader(format=self._format)
        self._load_checkpoint = loader

    def main(self, data: SeqsCollectionType, identifier=None) -> IdentifierType:
        from cogent3.app.composable import NotCompleted

        if isinstance(data, NotCompleted):
            return self.data_store.write_incomplete(identifier, data)

        if not data:
            msg = f"{self.__class__.__name__!r} does not support writing {data!r}"
            raise NotImplementedError(msg)

        if identifier is None:
            identifier = self._make_output_identifier(data)

        stored = self.data_store.write(identifier, self._formatter(data.to_dict()))
        if hasattr(data, "info"):
            data.info["stored"] = stored
        else:
            try:
                data.stored = stored
            except AttributeError:
                pass

        return stored


@define_app(app_type=LOADER)
class load_json:  # pragma: no cover
    """Loads json serialised cogent3 objects from a json file.
    Returns whatever object type was stored."""

    def __init__(self):
        from cogent3.util.warning import deprecated

        deprecated(
            "class",
            "load_json",
            "use cogent3.get_app('load_json')",
            "2023.3",
            "use cogent3.get_app()",
            stack_level=1,
        )

    def main(self, path: IdentifierType) -> SerialisableType:
        """returns object deserialised from json at path"""
        if type(path) == str:
            path = SingleReadDataStore(path)[0]

        data = path.read()
        identifier, data, completed = load_record_from_json(data)

        result = deserialise_object(data)
        if hasattr(result, "info"):
            result.info["source"] = result.info.get("source", identifier)
        else:
            try:
                identifier = getattr(result, "source", identifier)
                setattr(result, "source", identifier)
            except AttributeError:
                pass
        return result


@define_app(app_type=WRITER)
class write_json(_checkpointable):  # pragma: no cover
    """Writes json serialised objects to individual json files."""

    def __init__(self, data_path, name_callback=None, create=False, if_exists=SKIP):
        """
        Parameters
        ----------
        data_path
            path to write output, if ends with .zip will be a compressed zip
            archive
        name_callback
            function that takes the data object and returns a base
            file name
        create : bool
            whether to create the output directory
        if_exists : str
            behaviour if output exists. Either 'skip', 'raise' (raises an
            exception), 'overwrite'
        """
        from cogent3.util.warning import deprecated

        deprecated(
            "class",
            "write_json",
            "use cogent3.get_app('write_json')",
            "2023.3",
            "use cogent3.get_app()",
            stack_level=1,
        )

        super().__init__(
            data_path=data_path,
            name_callback=name_callback,
            create=create,
            if_exists=if_exists,
            suffix="json",
        )
        self._format = "json"

    def _set_checkpoint_loader(self):
        self._load_checkpoint = self

    def main(self, data: SerialisableType, identifier=None) -> IdentifierType:
        if isinstance(data, NotCompleted):
            return self.data_store.write_incomplete(identifier, data)

        if not data:
            msg = f"{self.__class__.__name__!r} does not support writing {data!r}"
            raise NotImplementedError(msg)

        if identifier is None:
            identifier = self._make_output_identifier(data)

        out = make_record_for_json(os.path.basename(identifier), data, True)
        out = json.dumps(out)
        stored = self.data_store.write(identifier, out)
        # todo is anything actually using this stored attribute? if not, delete this
        #  code and all other cases
        if hasattr(data, "info"):
            data.info["stored"] = stored
        else:
            try:
                data.stored = stored
            except AttributeError:
                pass
        return stored


@define_app(app_type=LOADER)
class load_db:  # pragma: no cover
    """Loads json serialised cogent3 objects from a TinyDB file.
    Returns whatever object type was stored."""

    def __init__(self):
        from cogent3.util.warning import deprecated

        deprecated(
            "class",
            "load_db",
            "use cogent3.get_app('load_db')",
            "2023.3",
            "use cogent3.get_app()",
            stack_level=1,
        )

    def main(self, identifier: IdentifierType) -> SerialisableType:
        """returns object deserialised from a TinyDb"""
        id_ = getattr(identifier, "id", None)
        if id_ is None:
            raise TypeError(
                f"{identifier} not connected to a TinyDB. If a json file path, use io.load_json()"
            )
        data = identifier.read()

        result = deserialise_object(data)
        if hasattr(result, "info"):
            result.info["source"] = result.info.get("source", identifier)
        else:
            try:
                identifier = getattr(result, "source", identifier)
                setattr(result, "source", identifier)
            except AttributeError:
                pass
        return result


@define_app(app_type=WRITER)
class write_db(_checkpointable):  # pragma: no cover
    """Writes json serialised objects to a TinyDB instance."""

    def __init__(self, data_path, name_callback=None, create=False, if_exists=SKIP):
        from cogent3.util.warning import deprecated

        deprecated(
            "class",
            "write_db",
            "use cogent3.get_app('write_db')",
            "2023.3",
            "use cogent3.get_app()",
            stack_level=1,
        )

        super().__init__(
            data_path=data_path,
            name_callback=name_callback,
            create=create,
            if_exists=if_exists,
            suffix="json",
            writer_class=WritableTinyDbDataStore,
        )

    def _set_checkpoint_loader(self):
        self._load_checkpoint = self

    T = Union[SerialisableType, IdentifierType]

    def main(self, data: SerialisableType, identifier=None) -> T:
        """
        Parameters
        ----------
        data
            object that has a `to_json()` method, or can be json serialised
        identifier : str
            if not provided, taken from data.source or data.info.source

        Returns
        -------
        identifier
        """
        identifier = identifier or get_data_source(data)
        identifier = self._make_output_identifier(identifier)

        if isinstance(data, NotCompleted):
            return self.data_store.write_incomplete(identifier, data)

        if not data:
            msg = f"{self.__class__.__name__!r} does not support writing {data!r}"
            raise NotImplementedError(msg)

        identifier = identifier or get_data_source(data)
        identifier = self._make_output_identifier(identifier)
        # todo revisit this when we establish immutability behaviour of database
        try:
            out = data.to_json()
        except AttributeError:
            out = json.dumps(data)
        stored = self.data_store.write(identifier, out)
        # todo is anything actually using this stored attribute? if not, delete this
        #  code and all other cases
        if hasattr(data, "info"):
            data.info["stored"] = stored
        else:
            try:
                data.stored = stored
            except AttributeError:
                pass
        return stored
