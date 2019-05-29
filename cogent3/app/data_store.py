import glob
import os
import re
import shutil
import zipfile
from fnmatch import fnmatch
from pathlib import Path
from pprint import pprint
from warnings import warn
from io import TextIOWrapper

from cogent3.util.misc import open_, get_format_suffixes, atomic_write

# handling archive, member existence
SKIP = 'skip'
OVERWRITE = 'overwrite'
RAISE = 'raise'
IGNORE = 'ignore'


class DataStoreMember(str):
    def __new__(klass, name, parent=None):
        result = str.__new__(klass, name)
        result.name = os.path.basename(name)
        result.parent = parent
        result._file = None
        return result

    def read(self):
        """returns contents"""
        return self.parent.read(self.name)

    def open(self):
        """returns file-like object"""
        if self._file is None:
            self._file = self.parent.open(self.name)
        return self._file

    def close(self):
        """closes file"""
        if self._file is None:
            return
        self._file.close()
        self._file = None


class ReadOnlyDataStoreBase:
    """a read only data store"""

    def __init__(self, source, suffix=None, limit=None, verbose=False):
        """
        Parameters
        ----------
        source
            path to directory / zip file
        suffix
            only members whose name matches the suffix are considered included
        limit
            the maximum number of members to consider
        verbose
            displays files that don't match search (applies only to the Zipped
            variant)
        """
        # assuming delimiter is /

        # todo this approach to caching persistent arguments for reconstruction
        # is fragile. Need an inspect module based approach
        d = locals()
        d.pop('self')
        self._persistent = d

        suffix = suffix or ''
        if suffix != '*':  # wild card search for all
            suffix = re.sub(r'^[\s.*]+', '', suffix)  # tidy the suffix
        source = re.sub(r'/+$', '', source)  # tidy the source

        self.suffix = suffix
        self.source = source
        self.mode = 'r'
        self._members = []
        self.limit = limit
        self._verbose = verbose

    def __getstate__(self):
        data = self._persistent.copy()
        return data

    def __setstate__(self, data):
        new = self.__class__(**data)
        self.__dict__.update(new.__dict__)
        return self

    def __repr__(self):
        if len(self) > 3:
            sample = str(list(self[:3]))
            sample = f'{sample[:-1]}...'
        else:
            sample = list(self)

        num = len(self)
        name = self.__class__.__name__
        txt = f'{num}x member {name}({sample})'
        return txt

    def __str__(self):
        return str(list(self))

    def head(self, n=5):
        """displays top n members"""
        pprint(self[:n])

    def tail(self, n=5):
        """displays last n members"""
        pprint(self[-n:])

    def __iter__(self):
        for member in self.members:
            yield DataStoreMember(self.get_absolute_identifier(member), self)

    def __getitem__(self, index):
        return self.members[index]

    def __len__(self):
        return len(self.members)

    def __contains__(self, identifier):
        """whether relative identifier has been stored"""
        identifier = self.get_relative_identifier(identifier)
        result = False
        for member in self.members:
            if identifier in member:
                result = True
                break
        return result

    def get_member(self, identifier):
        """returns DataStoreMember"""
        identifier = self.get_relative_identifier(identifier)
        for member in self.members:
            if identifier in member:
                return member
        return None

    def get_relative_identifier(self, identifier):
        """returns the identifier relative to store root path
        """
        source = self.source
        identifier = os.path.basename(identifier)
        if source.endswith('.zip'):
            source = source.replace('.zip', '')
            source = os.path.basename(source)
            identifier = f"{source}/{identifier}"
        else:
            if not isinstance(identifier, DataStoreMember):
                identifier = Path(identifier)
            identifier = identifier.name

        return identifier

    def get_absolute_identifier(self, identifier, from_relative=False):
        if not from_relative:
            identifier = self.get_relative_identifier(identifier)
        source = self.source.replace('.zip', '')
        if isinstance(identifier, DataStoreMember):
            identifier = identifier.name
        elif not identifier.startswith(source):
            identifier = f"{source}/{identifier}"
        return identifier

    def read(self, identifier):
        raise NotImplementedError  # override in subclasses

    @property
    def members(self):
        raise NotImplementedError  # override in subclasses

    def open(self, identifier):
        raise NotImplementedError

    def filtered(self, pattern=None, callback=None):
        """returns list of members for which callback returns True"""
        assert any([callback, pattern]
                   ), 'Must provide a pattern or a callback'
        if pattern:
            result = [m for m in self if fnmatch(m, pattern)]
        else:
            result = [m for m in self if callback(m)]
        return result


class ReadOnlyDirectoryDataStore(ReadOnlyDataStoreBase):
    @property
    def members(self):
        if not self._members:
            pattern = '%s/**/*.%s' % (self.source, self.suffix)
            paths = glob.iglob(pattern, recursive=True)
            members = []
            for i, path in enumerate(paths):
                if self.limit and i >= self.limit:
                    break
                member = DataStoreMember(self.get_absolute_identifier(path),
                                         self)
                members.append(member)
            self._members = members
        return self._members

    def read(self, identifier):
        infile = self.open(identifier)
        data = infile.read()
        infile.close()
        return data

    def open(self, identifier):
        identifier = self.get_absolute_identifier(identifier,
                                                  from_relative=False)
        if not os.path.exists(identifier):
            raise ValueError(f"path '{identifier}' does not exist")

        infile = open_(identifier)
        return infile


class SingleReadDataStore(ReadOnlyDirectoryDataStore):
    """simplified for a single file"""

    def __init__(self, source, *args, **kwargs):
        """
        Parameters
        source
            path to one file
        args
            ignored
        kwargs
            ignored
        """
        path = Path(source)
        assert path.exists() and path.is_file()
        super(SingleReadDataStore, self).__init__(str(path.parent),
                                                  suffix=str(path.suffix))
        self._members = [DataStoreMember(path, self)]


class ReadOnlyZippedDataStore(ReadOnlyDataStoreBase):
    @property
    def members(self):
        if os.path.exists(self.source) and not self._members:
            source_path = self.source.replace(Path(self.source).suffix, '')
            pattern = '*.%s' % self.suffix
            members = []
            with zipfile.ZipFile(self.source) as archive:
                names = archive.namelist()
                num_matches = 0
                for name in names:
                    name = os.path.basename(name)
                    if fnmatch(name, pattern):
                        num_matches += 1
                        member = DataStoreMember(
                            os.path.join(source_path, name), self)
                        members.append(member)
                    elif self._verbose:
                        print(f"Did not match {name}")

                    if self.limit and num_matches >= self.limit:
                        break
            self._members = members

        return self._members

    def read(self, identifier):
        record = self.open(identifier)
        data = record.read()
        record.close()
        return data

    def open(self, identifier):
        identifier = self.get_relative_identifier(identifier)
        archive = zipfile.ZipFile(self.source)
        record = archive.open(identifier)
        record = TextIOWrapper(record)
        return record


class WritableDataStoreBase:
    def __init__(self, if_exists=RAISE, create=False):
        """
        Parameters
        ----------
        if_exists : str
             behaviour when the destination already exists. Valid constants are
             defined in this file as OVERWRITE, SKIP, RAISE, IGNORE (they
             correspond to lower case version of the same word)
        create : bool
            if True, the destination is created
        """
        self._members = []
        if_exists = if_exists.lower()
        assert if_exists in (OVERWRITE, SKIP, RAISE, IGNORE)
        if create is False and if_exists == OVERWRITE:
            warn(f"'{OVERWRITE}' reset to '{IGNORE}' and create=True",
                 UserWarning)
            create = True
        self._source_create_delete(if_exists, create)

    def make_relative_identifier(self, data):
        """returns identifier for a new member relative to source"""
        if isinstance(data, DataStoreMember):
            data = data.name
        elif type(data) != str:
            try:
                data = data.info.source
            except AttributeError:
                try:
                    data = data.source
                except AttributeError:
                    raise ValueError('objects for storage require either a '
                                     'source or info.source string attribute')
        basename = os.path.basename(data)
        suffix, comp = get_format_suffixes(basename)
        if suffix and comp:
            pattern = f'.{suffix}.{comp}$'
        elif suffix:
            pattern = f'.{suffix}$'
        elif comp:
            pattern = f'.{comp}*$'
        else:
            raise ValueError(f"unknown name scheme {basename}")

        basename = re.sub(pattern, '', basename)
        basename = f'{basename}.{self.suffix}'
        return basename

    def make_absolute_identifier(self, data):
        """returns a absolute identifier for a new member, includes source"""
        basename = self.make_relative_identifier(data)
        identifier = self.get_absolute_identifier(
            basename, from_relative=True)
        return identifier


class WritableDirectoryDataStore(ReadOnlyDirectoryDataStore,
                                 WritableDataStoreBase):
    def __init__(self, source, suffix, mode='w', if_exists=RAISE, create=False):
        """
        Parameters
        ----------
        source
            path to directory / zip file
        suffix
            only members whose name matches the suffix are considered included
        mode : str
            file opening mode, defaults to write
        if_exists : str
             behaviour when the destination already exists. Valid constants are
             defined in this file as OVERWRITE, SKIP, RAISE, IGNORE (they
             correspond to lower case version of the same word)
        create : bool
            if True, the destination is created
        """
        assert 'w' in mode or 'a' in mode
        ReadOnlyDirectoryDataStore.__init__(
            self, source=source, suffix=suffix)
        WritableDataStoreBase.__init__(
            self, if_exists=if_exists, create=create)
        self.mode = mode

    def _source_create_delete(self, if_exists, create):
        exists = os.path.exists(self.source)
        if exists and if_exists == RAISE:
            raise RuntimeError(f"'{self.source}' exists")
        elif exists and if_exists == OVERWRITE:
            shutil.rmtree(self.source)
        elif not exists and not create:
            raise RuntimeError(f"'{self.source}' does not exist")

        if create:
            os.makedirs(self.source, exist_ok=True)

    def write(self, identifier, data):
        relative_id = self.get_relative_identifier(identifier)
        absolute_id = self.get_absolute_identifier(relative_id,
                                                   from_relative=True)

        with atomic_write(str(absolute_id), in_zip=False) as out:
            out.write(data)

        if relative_id not in self:
            self._members.append(DataStoreMember(relative_id, self))

        return absolute_id


class WritableZippedDataStore(ReadOnlyZippedDataStore, WritableDataStoreBase):
    def __init__(self, source, suffix, mode='a', if_exists=RAISE, create=False):
        """
        Parameters
        ----------
        source
            path to directory / zip file
        suffix
            only members whose name matches the suffix are considered included
        mode : str
            file opening mode, defaults to append
        if_exists : str
             behaviour when the destination already exists. Valid constants are
             defined in this file as OVERWRITE, SKIP, RAISE, IGNORE (they
             correspond to lower case version of the same word)
        create : bool
            if True, the destination is created
        """
        ReadOnlyZippedDataStore.__init__(self, source=source, suffix=suffix)
        WritableDataStoreBase.__init__(
            self, if_exists=if_exists, create=create)
        self.mode = 'a' or mode  # todo does mode 'w' nuke an entire zip?

    def _source_create_delete(self, if_exists, create):
        exists = os.path.exists(self.source)
        dirname = os.path.dirname(self.source)
        if exists and if_exists == RAISE:
            raise RuntimeError(f"'{self.source}' exists")
        elif exists and if_exists == OVERWRITE:
            os.remove(self.source)
        elif dirname and not os.path.exists(dirname) and not create:
            raise RuntimeError(f"'{dirname}' does not exist")

        if create and dirname:
            os.makedirs(dirname, exist_ok=True)

    def write(self, identifier, data):
        relative_id = self.get_relative_identifier(identifier)
        absolute_id = self.get_absolute_identifier(relative_id,
                                                   from_relative=True)

        with atomic_write(str(relative_id), in_zip=self.source) as out:
            out.write(data)

        if relative_id not in self:
            self._members.append(DataStoreMember(relative_id, self))

        return absolute_id
