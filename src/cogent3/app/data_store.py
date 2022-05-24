import glob
import json
import os
import pathlib
import re
import reprlib
import shutil
import weakref
import zipfile

from collections import defaultdict
from fnmatch import fnmatch, translate
from io import TextIOWrapper
from json import JSONDecodeError
from pathlib import Path
from pprint import pprint
from warnings import warn

from scitrack import get_text_hexdigest
from tinydb import Query, TinyDB
from tinydb.middlewares import CachingMiddleware
from tinydb.storages import JSONStorage

from cogent3.util.deserialise import deserialise_not_completed
from cogent3.util.io import atomic_write, get_format_suffixes, open_
from cogent3.util.misc import extend_docstring_from
from cogent3.util.parallel import is_master_process
from cogent3.util.table import Table
from cogent3.util.union_dict import UnionDict


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


# handling archive, member existence
SKIP = "skip"
OVERWRITE = "overwrite"
RAISE = "raise"
IGNORE = "ignore"


def get_data_source(data) -> str:
    """identifies attribute of data named 'source'

    Notes
    -----
    Alignment objects have a source element in their info dict
    """
    if isinstance(data, (str, pathlib.Path)):
        return str(data)

    if hasattr(data, "source"):
        return str(data.source)

    if hasattr(data, "info"):
        return get_data_source(data.info)

    if isinstance(data, dict):
        value = data.get("source")
        return str(value) if value else None

    return None


def make_record_for_json(identifier, data, completed):
    """returns a dict for storage as json"""
    try:
        data = data.to_rich_dict()
    except AttributeError:
        pass

    data = json.dumps(data)
    return dict(identifier=identifier, data=data, completed=completed)


def load_record_from_json(data):
    """returns identifier, data, completed status from json string"""
    if type(data) == str:
        data = json.loads(data)

    value = data["data"]
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except JSONDecodeError:
            pass

    return data["identifier"], value, data["completed"]


class DataStoreMember(str):
    def __new__(klass, name, parent=None, id=None):
        result = str.__new__(klass, name)
        result.name = os.path.basename(name)
        result.parent = parent
        result._file = None
        result.id = id
        return result

    def read(self):
        """returns contents"""
        return self.parent.read(self)

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

    @property
    def md5(self):
        return self.parent.md5(self, force=True)


class ReadOnlyDataStoreBase:

    store_suffix = None

    def __init__(self, source, suffix=None, limit=None, verbose=False, md5=True):
        """
        Parameters
        ----------
        source
            path to directory / zip file. Forced to end with store_suffix.
        suffix
            only members whose name matches the suffix are considered included
        limit
            the maximum number of members to consider
        verbose
            displays files that don't match search (applies only to the Zipped
            variant)
        md5 : bool
            record md5 hexadecimal checksum of read data when possible
        """
        # assuming delimiter is /

        # todo this approach to caching persistent arguments for reconstruction
        # is fragile. Need an inspect module based approach
        d = locals()
        self._persistent = UnionDict({k: v for k, v in d.items() if k != "self"})

        source = str(source)
        suffix = suffix or ""
        if suffix != "*":  # wild card search for all
            suffix = re.sub(r"^[\s.*]+", "", suffix)  # tidy the suffix
        source = re.sub(r"/+$", "", source)  # tidy the source

        self.suffix = suffix
        if self.store_suffix and not source.endswith(self.store_suffix):
            source = ".".join([source, self.store_suffix])
        self.source = str(pathlib.Path(source).expanduser())
        self.mode = "r"
        self._members = []
        self.limit = limit
        self._verbose = verbose
        self._md5 = md5
        self._checksums = {}

    def __getstate__(self):
        return self._persistent.copy()

    def __setstate__(self, data):
        new = self.__class__(**data)
        self.__dict__.update(new.__dict__)
        return self

    def __repr__(self):
        if len(self) > 3:
            sample = str(list(self[:3]))
            sample = f"{sample[:-1]}..."
        else:
            sample = list(self)

        num = len(self)
        name = self.__class__.__name__
        return f"{num}x member {name}(source='{self.source}', members={sample})"

    def __str__(self):
        return str(list(self))

    def head(self, n=5):
        """displays top n members"""
        pprint(self[:n])

    def tail(self, n=5):
        """displays last n members"""
        pprint(self[-n:])

    def __iter__(self):
        for i, member in enumerate(self.members):
            if not isinstance(member, DataStoreMember):
                member = DataStoreMember(self.get_absolute_identifier(member), self)
                self.members[i] = member
            yield member

    def __getitem__(self, index):
        return self.members[index]

    def __len__(self):
        return len(self.members)

    def __contains__(self, identifier):
        """whether relative identifier has been stored"""
        if isinstance(identifier, DataStoreMember):
            return identifier.parent is self

        if not identifier.endswith(self.suffix):
            suffix = pathlib.Path(identifier).suffix
            # possible an "added" file
            if self.store_suffix == "zip":
                klass = ReadOnlyZippedDataStore
            else:
                klass = ReadOnlyDirectoryDataStore
            new = klass(self.source, suffix=suffix)
            return identifier in new
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
        """returns the identifier relative to store root path"""
        if isinstance(identifier, DataStoreMember) and identifier.parent is self:
            return identifier

        source = self.source
        identifier = os.path.basename(identifier)
        if source.endswith(".zip"):
            # we insert the source path into identifier name
            # for zip members to ensure inflation creates a directory
            # containing them
            source = source.replace(".zip", "")
            source = os.path.basename(source)
            identifier = f"{source}{os.sep}{identifier}"
        else:
            identifier = Path(identifier)
            identifier = identifier.name

        return identifier

    def get_absolute_identifier(self, identifier, from_relative=False):
        """returns the identifier relative to the root path"""
        if not from_relative:
            identifier = self.get_relative_identifier(identifier)
        source = self.source.replace(".zip", "")
        if isinstance(identifier, DataStoreMember):
            identifier = identifier.name
        elif not identifier.startswith(source):
            identifier = f"{source}{os.sep}{identifier}"
        return identifier

    def read(self, identifier):
        """reads data corresponding to identifier"""
        if isinstance(identifier, DataStoreMember) and identifier.parent is self:
            identifier = identifier.name
        source = self.open(identifier)

        data = source.read()
        if self._md5:
            self._checksums[identifier] = get_text_hexdigest(data)
        source.close()
        return data

    @property
    def members(self):
        raise NotImplementedError  # override in subclasses

    def open(self, identifier):
        raise NotImplementedError

    def filtered(self, pattern=None, callback=None):
        """returns list of members for which callback returns True"""
        assert any([callback, pattern]), "Must provide a pattern or a callback"
        if pattern:
            result = [m for m in self if fnmatch(m, pattern)]
        else:
            result = [m for m in self if callback(m)]
        return result

    def md5(self, identifier, force=True):
        """
        Parameters
        ----------
        identifier
            name of data store member
        force : bool
            forces reading of data if not already done

        Returns
        -------
        md5 checksum for the member, if available, None otherwise
        """
        md5_setting = self._md5  # for restoring automatic md5 calc setting
        absoluteid = self.get_absolute_identifier(identifier)
        if force and absoluteid not in self._checksums:
            self._md5 = True
            _ = self.read(absoluteid)

        result = self._checksums.get(absoluteid, None)
        self._md5 = md5_setting
        return result


class ReadOnlyDirectoryDataStore(ReadOnlyDataStoreBase):
    @property
    def members(self):
        if not self._members:
            pattern = f"{self.source}/**/*.{self.suffix}"
            paths = glob.iglob(pattern, recursive=True)
            members = []
            for i, path in enumerate(paths):
                if self.limit and i >= self.limit:
                    break
                member = DataStoreMember(self.get_absolute_identifier(path), self)
                members.append(member)
            self._members = members
        return self._members

    def open(self, identifier):
        identifier = self.get_absolute_identifier(identifier, from_relative=False)
        if not os.path.exists(identifier):
            raise ValueError(f"path '{identifier}' does not exist")

        return open_(identifier)


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
        path = Path(source).expanduser()
        assert path.exists() and path.is_file()
        super(SingleReadDataStore, self).__init__(
            str(path.parent), suffix=str(path.suffix)
        )
        self._members = [DataStoreMember(path, self)]


class ReadOnlyZippedDataStore(ReadOnlyDataStoreBase):
    store_suffix = "zip"

    @property
    def members(self):
        if os.path.exists(self.source) and not self._members:
            source_path = self.source.replace(Path(self.source).suffix, "")
            pattern = f"*.{self.suffix}"
            members = []
            with zipfile.ZipFile(self.source) as archive:
                names = archive.namelist()
                num_matches = 0
                for name in names:
                    name = os.path.basename(name)
                    if fnmatch(name, pattern):
                        num_matches += 1
                        member = DataStoreMember(os.path.join(source_path, name), self)
                        members.append(member)
                    elif self._verbose:
                        print(f"Did not match {name}")

                    if self.limit and num_matches >= self.limit:
                        break
            self._members = members

        return self._members

    def open(self, identifier):
        identifier = self.get_relative_identifier(identifier)
        archive = zipfile.ZipFile(self.source)
        record = archive.open(identifier.replace("\\", "/"))
        record = TextIOWrapper(record, encoding="latin-1")
        return record


class WritableDataStoreBase:
    """a writeable data store"""

    def __init__(self, if_exists=RAISE, create=False):
        """
        if_exists : str
             behaviour when the destination already exists. Valid constants are
             defined in this file as OVERWRITE, RAISE, IGNORE (they
             correspond to lower case version of the same word)
        create : bool
            if True, the destination is created
        """
        d = locals()
        d = UnionDict({k: v for k, v in d.items() if k != "self"})
        if self._persistent:
            self._persistent |= d
        else:
            self._persistent = d

        self._members = []
        if_exists = if_exists.lower()
        assert if_exists in (OVERWRITE, SKIP, RAISE, IGNORE)
        self._source_create_delete(if_exists, create)

    def make_relative_identifier(self, data):
        """returns identifier for a new member relative to source"""

        if isinstance(data, DataStoreMember):
            data = data.name
        elif type(data) != str:
            data = get_data_source(data)
            if data is None:
                raise ValueError(
                    "objects for storage require either a "
                    "source or info.source string attribute"
                )
        basename = os.path.basename(data)
        suffix, comp = get_format_suffixes(basename)
        if suffix and comp:
            pattern = f".{suffix}.{comp}$"
        elif suffix:
            pattern = f".{suffix}$"
        elif comp:
            pattern = f".{comp}*$"
        else:
            pattern = None
        if pattern:
            basename = re.sub(pattern, "", basename)
        basename = f"{basename}.{self.suffix}"
        return basename

    def make_absolute_identifier(self, data):
        """returns a absolute identifier for a new member, includes source"""
        basename = self.make_relative_identifier(data)
        return self.get_absolute_identifier(basename, from_relative=True)

    def add_file(self, path, make_unique=True, keep_suffix=True, cleanup=False):
        """
        Parameters
        ----------
        path : str
            location of file to be added to the data store
        keep_suffix : bool
            new path will retain the suffix of the provided file
        make_unique : bool
            a successive number will be added to the name before the suffix
            until the name is unique
        cleanup : bool
            delete the original
        """
        relativeid = self.make_relative_identifier(path)
        relativeid = Path(relativeid)
        path = Path(path)
        if keep_suffix:
            relativeid = str(relativeid).replace(
                relativeid.suffix, "".join(path.suffixes)
            )
            relativeid = Path(relativeid)

        suffixes = "".join(relativeid.suffixes)
        new = relativeid
        num = 0
        while True:
            if not str(relativeid) in self:
                if num:
                    new = str(relativeid).replace(suffixes, f"-{num}{suffixes}")
                break

            num += 1
        relativeid = new
        data = SingleReadDataStore(path)[0].read()
        self.write(str(relativeid), data)

        if cleanup:
            path.unlink()

        return relativeid

    def write_incomplete(self, identifier, not_completed) -> DataStoreMember:
        """

        Parameters
        ----------
        identifier : str
            identifier for record
        not_completed : NotComplete
            instance that records key details for why incomplete
        Returns
        -------
        None if storage class does not support writing incomplete, otherwise
        a DataStoreMember.
        """
        if self.suffix != "json":
            msg = f"not supported for {self.__class__.__name__}"
            warn(msg, UserWarning)
            return

        record = make_record_for_json(identifier, not_completed, False)
        record = json.dumps(record)
        return self.write(identifier, record)

    def write(self, identifier, data, *args, **kwargs) -> DataStoreMember:
        """
        Parameters
        ----------
        identifier : str
            identifier that data will be saved under. Must have a suffix matching
            self.suffix or ``.log``.
        data
            data to be saved. If a tinydb, must be an object that can be
            converted to json, or has a to_json() method. Otherwise, it must be a string.

        Returns
        -------
        DataStoreMember instance
        """
        if not isinstance(data, str):
            raise TypeError(f"data must be a string type, not {type(data)}")

        id_suffix = identifier.split(".")[-1]
        if id_suffix not in (self.suffix, "log"):
            raise ValueError(
                f"identifier does not end with required suffix {self.suffix}"
            )

    def close(self):
        pass


class WritableDirectoryDataStore(ReadOnlyDirectoryDataStore, WritableDataStoreBase):
    @extend_docstring_from(ReadOnlyDirectoryDataStore.__init__, pre=False)
    @extend_docstring_from(WritableDataStoreBase.__init__, pre=False)
    def __init__(
        self,
        source,
        suffix,
        mode="w",
        if_exists=RAISE,
        create=False,
        md5=True,
        **kwargs,
    ):
        """
        md5 : bool
            record md5 hexadecimal checksum of data when possible
        mode : str
            file opening mode, defaults to write
        """
        assert "w" in mode or "a" in mode
        ReadOnlyDirectoryDataStore.__init__(self, source=source, suffix=suffix, md5=md5)
        WritableDataStoreBase.__init__(self, if_exists=if_exists, create=create)

        d = locals()
        self._persistent = {k: v for k, v in d.items() if k != "self"}
        self.mode = mode

    def _has_other_suffixes(self, path, suffix):
        p = Path(path)
        allowed = {str(suffix).lower(), "log"}
        for f in p.iterdir():
            if get_format_suffixes(str(f))[0] not in allowed:
                return True
        return False

    def _source_create_delete(self, if_exists, create):
        if not is_master_process():
            return

        path = Path(self.source)
        if path.exists() and if_exists == RAISE:
            raise FileExistsError(f"'{self.source}' exists")
        elif path.exists() and if_exists == OVERWRITE:
            if self._has_other_suffixes(self.source, self.suffix):
                raise RuntimeError(
                    f"Unsafe to delete {self.source} as it contains ",
                    f"files other than .{self.suffix} or .log files."
                    " You will need to remove this directly yourself.",
                )
            try:
                shutil.rmtree(self.source)
            except NotADirectoryError:
                os.remove(self.source)
        elif not path.exists() and not create:
            raise FileNotFoundError(f"'{self.source}' does not exist")

        if create:
            path.mkdir(parents=True, exist_ok=True)

    @extend_docstring_from(WritableDataStoreBase.write)
    def write(self, identifier, data):
        if not data:
            return data

        super().write(identifier, data)
        id_suffix = identifier.split(".")[-1]
        if id_suffix not in (self.suffix, "log"):
            raise ValueError(
                f"identifier does not end with required suffix {self.suffix}"
            )

        relative_id = self.get_relative_identifier(identifier)
        absolute_id = self.get_absolute_identifier(relative_id, from_relative=True)

        if self._md5:
            self._checksums[absolute_id] = get_text_hexdigest(data)

        with atomic_write(str(absolute_id), in_zip=False) as out:
            out.write(data)

        member = DataStoreMember(relative_id, self)
        if relative_id not in self and relative_id.endswith(self.suffix):
            self._members.append(member)

        return member


def _db_lockid(path):
    """returns value for pid in LOCK record or None"""
    if not os.path.exists(path):
        return None

    with TinyDB(path) as db:
        query = Query().identifier.matches("LOCK")
        got = db.get(query)
        lockid = None if not got else got["pid"]

    return lockid


class ReadOnlyTinyDbDataStore(ReadOnlyDataStoreBase):
    """A TinyDB based json data store"""

    store_suffix = "tinydb"

    @extend_docstring_from(ReadOnlyDirectoryDataStore.__init__)
    def __init__(self, *args, **kwargs):
        kwargs["suffix"] = "json"
        super(ReadOnlyTinyDbDataStore, self).__init__(*args, **kwargs)
        self._db = None
        self._finish = None

    def __contains__(self, identifier):
        """whether identifier has been stored here"""
        if isinstance(identifier, DataStoreMember):
            return identifier.parent is self

        query = Query().identifier.matches(identifier)
        return self.db.contains(query)

    def __repr__(self):
        txt = super().__repr__()
        query = Query().completed == False
        num = self.db.count(query)
        if num > 0:
            txt = f"{txt}, {num}x incomplete"
        return txt

    @property
    def db(self):
        if self._db is None:
            storage = CachingMiddleware(JSONStorage)
            storage.WRITE_CACHE_SIZE = 50  # todo support for user specifying
            self._db = TinyDB(self.source, storage=storage)
            name = self.__class__.__name__
            if "readonly" in name.lower():
                # remove interface for inserting records making this a read only db
                self._db.insert = None
            else:
                self.lock()

            self._finish = weakref.finalize(self, self._close, self._db)

        return self._db

    def __del__(self):
        self.close()

    @classmethod
    def _close(cls, db):
        try:
            db.storage.flush()
            db.close()
        except ValueError:
            # file probably already closed
            pass

    def close(self):
        """closes the data store"""
        try:
            self.unlock()
            self.db.storage.flush()
        except ValueError:
            # file probably already closed
            pass
        self._finish()
        self._finish.detach()

    def lock(self):
        """if writable, and not locked, locks the database to this pid"""
        if not self.locked:
            self._db.insert(dict(identifier="LOCK", pid=os.getpid()))
            self._db.storage.flush()

    @property
    def locked(self):
        """returns lock pid or None if unlocked or pid matches self"""
        return _db_lockid(self.source) is not None

    def unlock(self, force=False):
        """remove a lock if pid matches. If force, ignores pid."""
        if "readonly" in self.__class__.__name__:
            # not allowed to touch a lock
            return

        query = Query().identifier.matches("LOCK")
        got = self.db.get(query)
        if not got:
            return

        lock_id = got["pid"]
        if lock_id == os.getpid() or force:
            self.db.remove(query)
            self.db.storage.flush()

        return lock_id

    @property
    def incomplete(self):
        """returns database records with completed=False"""
        query = Query().completed == False
        incomplete = []
        for record in self.db.search(query):
            member = DataStoreMember(record["identifier"], self, id=record.doc_id)
            incomplete.append(member)
        return incomplete

    @property
    def summary_incomplete(self):
        """returns a table summarising incomplete results"""
        # detect last exception line
        err_pat = re.compile(r"[A-Z][a-z]+[A-Z][a-z]+\:.+")
        types = defaultdict(list)
        indices = "type", "origin"
        for member in self.incomplete:
            record = member.read()
            record = deserialise_not_completed(record)
            key = tuple(getattr(record, k, None) for k in indices)
            match = err_pat.findall(record.message)
            types[key].append([match[-1] if match else record.message, record.source])

        header = list(indices) + ["message", "num", "source"]
        rows = []
        maxtring = reprlib.aRepr.maxstring
        reprlib.aRepr.maxstring = 45

        for record in types:
            messages, sources = list(zip(*types[record]))
            messages = reprlib.repr(
                ", ".join(m.splitlines()[-1] for m in set(messages))
            )
            sources = reprlib.repr(", ".join(s.splitlines()[-1] for s in sources))
            row = list(record) + [
                messages,
                len(types[record]),
                sources,
            ]
            rows.append(row)

        reprlib.aRepr.maxstring = maxtring  # restoring original val

        return Table(header=header, data=rows, title="incomplete records")

    @property
    def members(self):
        if not self._members:
            pattern = translate(f"*.{self.suffix}") if self.suffix else translate("*")
            members = []
            query = Query()
            query = (query.identifier.matches(pattern)) & (query.completed == True)
            for record in self.db.search(query):
                member = DataStoreMember(record["identifier"], self, id=record.doc_id)
                members.append(member)

                if self.limit and len(members) >= self.limit:
                    break
            self._members = members

        return self._members

    @extend_docstring_from(ReadOnlyDataStoreBase.get_absolute_identifier, pre=True)
    def get_absolute_identifier(self, identifier, from_relative=True):
        """For tinydb, this is the same as the relative identifier"""
        return self.get_relative_identifier(identifier)

    @extend_docstring_from(ReadOnlyDataStoreBase.get_relative_identifier)
    def get_relative_identifier(self, identifier):
        if isinstance(identifier, DataStoreMember) and identifier.parent is self:
            return identifier

        identifier = Path(identifier)
        identifier = identifier.name
        return identifier

    def open(self, identifier):
        if getattr(identifier, "parent", None) is not self:
            member = self.get_member(identifier)
        else:
            member = identifier

        _, record, _ = load_record_from_json(self.db.get(doc_id=member.id))
        return record

    def read(self, identifier):
        data = self.open(identifier)
        if self._md5 and isinstance(data, str):
            self._checksums[identifier] = get_text_hexdigest(data)

        return data

    @extend_docstring_from(ReadOnlyDataStoreBase.md5)
    def md5(self, member, force=True):
        md5_setting = self._md5  # for restoring automatic md5 calc setting
        if not getattr(member, "id", None):
            member = self.filtered(member)[0]

        if force and member not in self._checksums:
            self._md5 = True
            _ = member.read()

        result = self._checksums.get(member, None)
        self._md5 = md5_setting
        return result

    @property
    def logs(self):
        """returns all records with a .log suffix"""
        logfiles = []
        query = Query().identifier.matches(translate("*.log"))
        for record in self.db.search(query):
            member = DataStoreMember(record["identifier"], self, id=record.doc_id)
            logfiles.append(member)
        return logfiles

    @property
    def summary_logs(self):
        """returns a table summarising log files"""
        rows = []
        for record in self.logs:
            data = record.read().splitlines()
            first = data.pop(0).split("\t")
            row = [first[0], record.name]
            key = None
            mapped = {}
            for line in data:
                line = line.split("\t")[-1].split(" : ", maxsplit=1)
                if len(line) == 1:
                    mapped[key] += line[0]
                    continue

                key = line[0]
                mapped[key] = line[1]

            data = mapped
            row.extend(
                [
                    data["python"],
                    data["user"],
                    data["command_string"],
                    data["composable function"],
                ]
            )
            rows.append(row)
        return Table(
            header=["time", "name", "python version", "who", "command", "composable"],
            data=rows,
            title="summary of log files",
        )

    @property
    def describe(self):
        """returns tables describing content types"""
        lock_id = _db_lockid(self.source)
        if lock_id:
            title = (
                f"Locked db store. Locked to pid={lock_id}, current pid={os.getpid()}"
            )
        else:
            title = "Unlocked db store."
        num_incomplete = len(self.incomplete)
        num_complete = len(self.members)
        num_logs = len(self.logs)
        return Table(
            header=["record type", "number"],
            data=[
                ["completed", num_complete],
                ["incomplete", num_incomplete],
                ["logs", num_logs],
            ],
            title=title,
        )


class WritableTinyDbDataStore(ReadOnlyTinyDbDataStore, WritableDataStoreBase):
    @extend_docstring_from(WritableDirectoryDataStore.__init__)
    def __init__(self, *args, **kwargs):
        """

        Notes
        -----
        A TinyDb file can be locked. In which case, ``if_exists=OVERWRITE``
        will be converted to RAISE.
        """
        if_exists = kwargs.pop("if_exists", RAISE)
        create = kwargs.pop("create", True)
        ReadOnlyTinyDbDataStore.__init__(self, *args, **kwargs)
        WritableDataStoreBase.__init__(self, if_exists=if_exists, create=create)

    def _source_create_delete(self, if_exists, create):
        if not is_master_process():
            return

        path = Path(self.source)
        if if_exists == OVERWRITE and path.exists():
            try:
                path.unlink()
            except PermissionError:
                # probably user accidentally created a directory
                shutil.rmtree(path)
            return

        locked_id = _db_lockid(self.source)
        pid = os.getpid()
        if path.exists() and if_exists == RAISE:
            msg = f"'{path}' exists"
            if locked_id is not None:
                msg = (
                    f"{msg}, and is locked by process pid {locked_id}."
                    f" Current pid is {pid}"
                )
            raise FileExistsError(msg)

        if if_exists == IGNORE and locked_id is not None:
            warn(f"'{self.source}' is locked to {locked_id}, current pid is {pid}.")
            return

        if path.parent and not path.parent.exists() and not create:
            raise FileNotFoundError(f"'{path.parent}' does not exist, set create=True")

        path.parent.mkdir(parents=True, exist_ok=True)

    @extend_docstring_from(WritableDataStoreBase.write)
    def write(self, identifier, data):
        # writing into a tinydb has its own logic for conversion to json
        #  so we don't validate data is a string for this case
        from cogent3.app.composable import NotCompleted

        if isinstance(data, NotCompleted):
            return self.write_incomplete(identifier, data)

        super().write(identifier, "")
        id_suffix = identifier.split(".")[-1]
        if id_suffix not in (self.suffix, "log"):
            raise ValueError(
                f"identifier does not end with required suffix {self.suffix}"
            )

        matches = self.filtered(identifier)
        if matches:
            return matches[0]

        relative_id = self.get_relative_identifier(identifier)
        record = make_record_for_json(relative_id, data, True)
        doc_id = self.db.insert(record)

        member = DataStoreMember(relative_id, self, id=doc_id)
        if relative_id.endswith(self.suffix):
            self._members.append(member)

        return member

    def write_incomplete(self, identifier, not_completed):
        """stores an incomplete result object"""

        matches = self.filtered(identifier)
        if matches:
            return matches[0]

        relative_id = self.get_relative_identifier(identifier)
        record = make_record_for_json(relative_id, not_completed, False)
        doc_id = self.db.insert(record)

        return DataStoreMember(relative_id, self, id=doc_id)

    def add_file(self, path, make_unique=True, keep_suffix=True, cleanup=False):
        """
        Parameters
        ----------
        path : str
            location of file to be added to the data store
        keep_suffix : bool
            new path will retain the suffix of the provided file
        make_unique : bool
            a successive number will be added to the name before the suffix
            until the name is unique
        cleanup : bool
            delete the original
        """
        relativeid = self.make_relative_identifier(path)
        relativeid = Path(relativeid)
        path = Path(path)
        if keep_suffix:
            relativeid = str(relativeid).replace(
                relativeid.suffix, "".join(path.suffixes)
            )
            relativeid = Path(relativeid)

        if relativeid.suffix:
            name_wo_suffix = ".".join(relativeid.name.split(".")[:-1])
        else:
            name_wo_suffix = relativeid
        suffixes = "".join(relativeid.suffixes)
        query = Query().identifier.matches(f"{name_wo_suffix}*")
        num = self.db.count(query)
        if num:
            num += 1
            relativeid = str(relativeid).replace(suffixes, f"-{num}{suffixes}")
            relativeid = str(relativeid).replace(suffixes, f"-{num}{suffixes}")

        data = path.read_text()
        m = self.write(str(relativeid), data)

        if cleanup:
            path.unlink()

        return m
