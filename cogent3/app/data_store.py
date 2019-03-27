import glob
import os
import re
import shutil
import zipfile
from pathlib import Path
from warnings import warn
from fnmatch import fnmatch

from cogent3.util.misc import open_, get_format_suffixes

# handling archive, member existence
SKIP = 'skip'
OVERWRITE = 'overwrite'
RAISE = 'raise'
IGNORE = 'ignore'


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
        suffix = suffix or ''
        suffix = re.sub(r'^[\s.*]+', '', suffix)  # tidy the suffix
        source = re.sub(r'/+$', '', source)  # tidy the source

        self.suffix = suffix
        self.source = source
        self.mode = 'r'
        self._members = []
        self.limit = limit
        self._verbose = verbose

    def __contains__(self, identifier):
        """whether relative identifier has been stored"""
        identifier = self.get_relative_identifier(identifier)
        return identifier in self.members

    def get_relative_identifier(self, identifier):
        """returns the relative identifier"""
        identifier = re.sub(f'{self.source}/', '', identifier)
        return identifier

    def get_absolute_identifier(self, identifier, from_relative=False):
        if not from_relative:
            identifier = self.get_relative_identifier(identifier)

        identifier = f"{self.source}/{identifier}"
        return identifier

    def read(self, identifier):
        raise NotImplementedError  # override in subclasses

    @property
    def members(self):
        raise NotImplementedError  # override in subclasses


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
                members.append(self.get_relative_identifier(path))
            self._members = members
        return self._members

    def read(self, identifier):
        identifier = self.get_absolute_identifier(identifier,
                                                  from_relative=False)
        if not os.path.exists(identifier):
            raise ValueError(f"path '{identifier}' does not exist")

        with open_(identifier) as infile:
            data = infile.read()

        return data


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
        self._members = [str(path.name)]


class ReadOnlyZippedDataStore(ReadOnlyDataStoreBase):
    @property
    def members(self):
        if os.path.exists(self.source) and not self._members:
            pattern = '*.%s' % self.suffix
            members = []
            with zipfile.ZipFile(self.source) as archive:
                names = archive.namelist()
                num_matches = 0
                for name in names:
                    if fnmatch(name, pattern):
                        num_matches += 1
                        members.append(name)
                    elif self._verbose:
                        print(f"Did not match {name}")

                    if self.limit and num_matches >= self.limit:
                        break
            self._members = members

        return self._members

    def read(self, identifier):
        identifier = self.get_relative_identifier(identifier)
        with zipfile.ZipFile(self.source) as archive:
            data = archive.read(identifier)
        data = data.decode('utf8')
        return data



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
        if type(data) != str:
            data = data.info.source
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
        identifier = self.get_absolute_identifier(basename, from_relative=True)
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
        ReadOnlyDirectoryDataStore.__init__(self, source=source, suffix=suffix)
        WritableDataStoreBase.__init__(self, if_exists=if_exists, create=create)
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
        with open_(absolute_id, self.mode) as outfile:
            outfile.write(data)

        if relative_id not in self:
            self._members.append(relative_id)

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
        WritableDataStoreBase.__init__(self, if_exists=if_exists, create=create)
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

        if create:
            os.makedirs(dirname, exist_ok=True)

    def write(self, identifier, data):
        relative_id = self.get_relative_identifier(identifier)
        absolute_id = self.get_absolute_identifier(relative_id,
                                                   from_relative=True)
        with zipfile.ZipFile(self.source, 'a') as out:
            out.writestr(relative_id, data, compress_type=zipfile.ZIP_DEFLATED,
                         compresslevel=9)

        if relative_id not in self:
            self._members.append(relative_id)

        return absolute_id
