import inspect
import json
import os
import pathlib
import re
import time
import traceback

import scitrack

from cogent3 import make_aligned_seqs, make_unaligned_seqs
from cogent3.core.alignment import SequenceCollection
from cogent3.util import progress_display as UI
from cogent3.util.misc import get_object_provenance, open_

from .data_store import (
    IGNORE,
    OVERWRITE,
    RAISE,
    SKIP,
    WritableDirectoryDataStore,
    WritableZippedDataStore,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.8.30a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def _make_logfile_name(process):
    text = str(process)
    text = re.split(r"\s+\+\s+", text)
    parts = []
    for part in text:
        parts.append(part[: part.find("(")])
    result = "-".join(parts)
    pid = os.getpid()
    result = f"{result}-pid{pid}.log"
    return result


def _get_source(source):
    if isinstance(source, str):
        result = str(source)
    else:
        try:
            result = source.source
        except AttributeError:
            try:
                result = source.info.source
            except AttributeError:
                result = None
    return result


def _get_origin(origin):
    if type(origin) == str:
        result = origin
    else:
        result = origin.__class__.__name__
    return result


class NotCompleted(int):
    """results that failed to complete"""

    def __new__(cls, type, origin, message, source=None):
        """
        Parameters
        ----------
        type : str
            examples are 'ERROR', 'FAIL'
        origin
            where the instance was created, can be an instance
        message : str
            descriptive message, succinct traceback
        source : str or instance with .info.source or .source attributes
            the data operated on that led to this result. Can
        """
        # todo this approach to caching persistent arguments for reconstruction
        # is fragile. Need an inspect module based approach
        origin = _get_origin(origin)
        source = _get_source(source)
        d = locals()
        d = {k: v for k, v in d.items() if k != "cls"}
        result = int.__new__(cls, False)
        args = tuple(d.pop(v) for v in ("type", "origin", "message"))
        result._persistent = args, d

        result.type = type
        result.origin = origin
        result.message = message
        result.source = source
        return result

    def __getnewargs_ex__(self, *args, **kw):
        return self._persistent[0], self._persistent[1]

    def __repr__(self):
        return str(self)

    def __str__(self):
        name = self.__class__.__name__
        source = self.source or "Unknown"
        val = (
            f"{name}(type={self.type}, origin={self.origin}, "
            f'source="{source}", message="{self.message}")'
        )
        return val

    def to_rich_dict(self):
        """returns components for to_json"""
        data = {
            "type": get_object_provenance(self),
            "not_completed_construction": dict(
                args=self._persistent[0], kwargs=self._persistent[1]
            ),
        }
        return data

    def to_json(self):
        """returns json string"""
        return json.dumps(self.to_rich_dict())


class ComposableType:
    _type = None

    def __init__(self, input_types=None, output_types=None, data_types=None):
        """
        Parameters
        ----------
        input_types : str or collection of str
            Allowed input types
        output_types : str or collection of str
            Types of output
        data_types : str or collection of str
            Allowed data types
        """
        input_types = [] if input_types is None else input_types
        output_types = [] if output_types is None else output_types
        data_types = [] if data_types is None else data_types

        if type(input_types) == str:
            input_types = [input_types]
        if type(output_types) == str:
            output_types = [output_types]
        if type(data_types) == str:
            data_types = [data_types]

        self._input_types = frozenset(input_types)
        self._output_types = frozenset(output_types)
        self._data_types = frozenset(data_types)

    def compatible_input(self, other):
        result = other._output_types & self._input_types
        return result != set()

    def _validate_data_type(self, data):
        """checks data class name matches defined compatible types"""
        if not self._data_types:
            # not defined
            return True

        name = data.__class__.__name__
        valid = False
        for type_ in self._data_types:
            if type_ == name:
                valid = True
                break
        if not valid:
            msg = f"invalid data type, '{name}' not in {', '.join(self._data_types)}"
            valid = NotCompleted("ERROR", self, message=msg, source=data)
        return valid


class Composable(ComposableType):
    def __init__(self, **kwargs):
        super(Composable, self).__init__(**kwargs)
        self.func = None  # over-ride in subclass
        self._in = None  # input rules
        self._out = None  # rules receiving output
        # rules operating on result but not part of a chain
        self._checkpointable = False
        self._load_checkpoint = None
        self._formatted = ["type='%s'" % self._type]

    def __str__(self):
        txt = "" if not self.input else str(self.input)
        if txt:
            txt += " + "
        txt += "%s(%s)" % (self.__class__.__name__, ", ".join(self._formatted))
        return txt

    def __repr__(self):
        return str(self)

    def _formatted_params(self):
        stack = inspect.stack()
        stack.reverse()
        for level in stack:
            args = inspect.getargvalues(level.frame).locals
            klass = args.get("__class__", None)
            if klass and isinstance(self, klass):
                break
        args = inspect.getargvalues(level.frame).locals
        params = inspect.signature(args["__class__"].__init__).parameters
        formatted = []
        for p in params:
            if p == "self":
                continue
            try:
                v = args[p]
            except KeyError:
                continue
            try:
                v = v.name
            except AttributeError:
                pass
            try:
                get_ipython()
                if p == "kwargs" and v == {"store_history": True, "silent": False}:
                    continue
            except NameError:
                pass
            formatted.append("%s=%r" % (p, v))
        self._formatted += formatted

    def __add__(self, other):
        if self.output or other.input:
            # can only be part of a single composable function
            self_name = self.__class__.__name__
            other_name = other.__class__.__name__
            if self.output and other.input:
                msg = (
                    f"{self_name} and {other_name} are already part of a"
                    " composed function"
                )
            elif self.output:
                msg = f"{self_name} already part of composed function"
            else:
                msg = f"{other_name} already part of composed function"
            raise AssertionError(f"{msg}, use disconnect() to free them up")

        if not other.compatible_input(self):
            msg = '%s() requires input type "%s", %s() produces "%s"'
            my_name = self.__class__.__name__
            my_output = self._output_types
            their_name = other.__class__.__name__
            their_input = other._input_types
            msg = msg % (their_name, their_input, my_name, my_output)
            raise TypeError(msg)
        self.output = other
        other.input = self
        return other

    def _set_checkpoint_loader(self):
        # over-ride in subclasses that are checkpointable
        self._load_checkpoint = None

    def _make_output_identifier(self, data):
        # over-ride in subclasses that are checkpointable
        pass

    @property
    def checkpointable(self):
        """whether this function is checkpointable"""
        return self._checkpointable

    def job_done(self, *args, **kwargs):
        # over-ride in sub classes
        return False

    def _trapped_call(self, func, val, *args, **kwargs):
        valid = self._validate_data_type(val)
        if not valid:
            return valid
        try:
            val = func(val, *args, **kwargs)
        except Exception:
            val = NotCompleted("ERROR", self, traceback.format_exc(), source=val)
        return val

    def __call__(self, val, *args, **kwargs):
        # initial invocation always transfers call() to first composable
        # element to get input for self
        refobj = kwargs.get("identifier", val)
        if not val:
            return val

        if self.checkpointable:
            job_done = self.job_done(refobj)
            if job_done and self.output:
                result = self._load_checkpoint(refobj)
            elif job_done:
                result = self._make_output_identifier(refobj)

            if job_done:
                return result

        if self.input:
            val = self._in(val, *args, **kwargs)

        if not val:
            return val
        result = self._trapped_call(self.func, val, *args, **kwargs)
        if not result and type(result) != NotCompleted:
            msg = (
                f"The value {result} equates to False. "
                "If unexpected, please post this error message along"
                f" with the code and data '{val}' as an Issue on the"
                " bitbucket project page."
            )
            origin = str(self)
            result = NotCompleted("BUG", origin, msg, source=val)

        return result

    @property
    def input(self):
        return self._in

    @input.setter
    def input(self, other):
        self._in = other

    @property
    def output(self):
        return self._out

    @output.setter
    def output(self, other):
        self._out = other
        self._set_checkpoint_loader()

    def disconnect(self):
        """resets input and output to None

        Breaks all connections among members of a composed function."""
        input = self.input
        if input:
            input.disconnect()

        self._in = None
        self._out = None
        self._load_checkpoint = None

    @UI.display_wrap
    def apply_to(
        self,
        dstore,
        parallel=False,
        mininterval=2,
        par_kw=None,
        logger=True,
        cleanup=False,
        ui=None,
    ):
        """invokes self composable function on the provided data store

        Parameters
        ----------
        dstore
            a path, list of paths, or DataStore to which the process will be
            applied.
        parallel : bool
            run in parallel, according to arguments in par_kwargs. If True,
            the last step of the composable function serves as the master
            process, with earlier steps being executed in parallel for each
            member of dstore.
        par_kw
            dict of values for configuring parallel execution.
        logger
            Argument ignored if not an io.writer. A scitrack logger, a logfile
            name or True. If True, a scitrack logger is created with a name that
            defaults to the composable function names and the process ID,
            e.g. load_unaligned-progressive_align-write_seqs-pid6962.log.
            If string, that name is used as the logfile name. Otherwise the
            logger is used as is.
        cleanup : bool
            after copying of log files into the data store, they are deleted
            from their original location

        Returns
        -------
        Result of the process as a list
        Notes
        -----
        If run in parallel, this instance serves as the master object and
        aggregates results.
        """
        if isinstance(dstore, str):
            dstore = [dstore]

        dstore = [e for e in dstore if e]
        if len(dstore) == 0:
            raise ValueError("dstore is empty")

        start = time.time()
        loggable = hasattr(self, "data_store")
        if not loggable:
            LOGGER = None
        elif type(logger) == scitrack.CachingLogger:
            LOGGER = logger
        elif type(logger) == str:
            LOGGER = scitrack.CachingLogger
            LOGGER.log_file_path = logger
        elif logger == True:
            log_file_path = pathlib.Path(_make_logfile_name(self))
            source = pathlib.Path(self.data_store.source)
            log_file_path = source.parent / log_file_path
            LOGGER = scitrack.CachingLogger()
            LOGGER.log_file_path = str(log_file_path)
        else:
            LOGGER = None

        if LOGGER:
            LOGGER.log_message(str(self), label="composable function")
            LOGGER.log_versions(["cogent3"])
        results = []
        i = 0
        process = self.input if self.input else self
        if self.input:
            # As we will be explicitly calling the input object, we disconnect
            # the two-way interaction between input and self. This means self
            # is not called twice, and self is not unecessarily pickled during
            # parallel execution.
            process.output = None
            self.input = None

        # with a tinydb dstore, this also excludes data that failed to complete
        todo = [m for m in dstore if not self.job_done(m)]

        for result in ui.imap(
            process, todo, parallel=parallel, par_kw=par_kw, mininterval=mininterval
        ):
            outcome = self(result)
            results.append(outcome)
            if LOGGER:
                member = dstore[i]
                LOGGER.log_message(member, label="input")
                LOGGER.log_message(member.md5, label="input md5sum")
                mem_id = self.data_store.make_relative_identifier(member.name)
                if outcome:
                    member = self.data_store.get_member(mem_id)
                    LOGGER.log_message(member, label="output")
                    LOGGER.log_message(member.md5, label="output md5sum")
                else:
                    # we have a NotCompletedResult
                    try:
                        # tinydb supports storage
                        self.data_store.write_incomplete(mem_id, outcome.to_rich_dict())
                    except AttributeError:
                        pass
                    LOGGER.log_message(
                        f"{outcome.origin} : {outcome.message}", label=outcome.type
                    )

            i += 1

        finish = time.time()
        taken = finish - start
        if LOGGER:
            LOGGER.log_message(f"{taken}", label="TIME TAKEN")
            LOGGER.shutdown()
            log_file_path = str(log_file_path)
            self.data_store.add_file(log_file_path, cleanup=cleanup, keep_suffix=True)
            self.data_store.close()

        # now reconnect input
        if process is not self:
            self = process + self

        return results


class ComposableTabular(Composable):
    _type = "tabular"

    def __init__(self, **kwargs):
        super(ComposableTabular, self).__init__(**kwargs)


class ComposableSeq(Composable):

    _type = "sequences"

    def __init__(self, **kwargs):
        super(ComposableSeq, self).__init__(**kwargs)


class ComposableAligned(Composable):
    _type = "aligned"

    def __init__(self, **kwargs):
        super(ComposableAligned, self).__init__(**kwargs)


class ComposableTree(Composable):
    _type = "tree"

    def __init__(self, **kwargs):
        super(ComposableTree, self).__init__(**kwargs)


class ComposableModel(Composable):
    _type = "model"

    def __init__(self, **kwargs):
        super(ComposableModel, self).__init__(**kwargs)


class ComposableHypothesis(Composable):
    _type = "hypothesis"

    def __init__(self, **kwargs):
        super(ComposableHypothesis, self).__init__(**kwargs)


class ComposableDistance(Composable):
    _type = "distance"

    def __init__(self, **kwargs):
        super(ComposableDistance, self).__init__(**kwargs)


class _seq_loader:
    def load(self, data):
        """returns sequences

        Parameters
        ----------
        data
            file path or cogent3 sequence collection / alignment
        """
        if type(data) == str:
            with open_(data) as infile:
                data = dict(record for record in self._parser(infile))
            seqs = self.klass(data=data, moltype=self.moltype)
            seqs.info.path = data
        elif not isinstance(data, SequenceCollection):
            if self.aligned:
                seqs = make_aligned_seqs(data, moltype=self.moltype)
            else:
                seqs = make_unaligned_seqs(data, moltype=self.moltype)

        if self._output_types == "sequence":
            seqs = seqs.degap()

        return seqs


class _checkpointable(Composable):
    def __init__(
        self,
        data_path,
        name_callback=None,
        create=False,
        if_exists=SKIP,
        suffix=None,
        writer_class=None,
        **kwargs,
    ):
        """
        Parameters
        ----------
        data_path
            path to write output
        name_callback
            function that takes the data object and returns a base
            file name
        create : bool
            whether to create the data_path reference
        if_exists : str
            behaviour if output exists. Either 'skip', 'raise' (raises an
            exception), 'overwrite', 'ignore'
        writer_class : type
            constructor for writer
        """
        super(_checkpointable, self).__init__(**kwargs)
        self._formatted_params()

        self._checkpointable = True
        if_exists = if_exists.lower()
        assert if_exists in (
            SKIP,
            IGNORE,
            RAISE,
            OVERWRITE,
        ), "invalid value for if_exists"
        self._if_exists = if_exists

        if writer_class:
            klass = writer_class
        else:
            klass = (
                WritableZippedDataStore
                if data_path.endswith(".zip")
                else WritableDirectoryDataStore
            )
        self.data_store = klass(
            data_path, suffix=suffix, create=create, if_exists=if_exists
        )
        self._callback = name_callback
        self.func = self.write

        # override the following in subclasses
        self._format = None
        self._formatter = None
        self._load_checkpoint = None

    def _make_output_identifier(self, data):
        if self._callback:
            data = self._callback(data)

        identifier = self.data_store.make_absolute_identifier(data)
        return identifier

    def job_done(self, data):
        identifier = self._make_output_identifier(data)
        exists = identifier in self.data_store
        if exists and self._if_exists == RAISE:
            msg = "'%s' already exists" % identifier
            raise RuntimeError(msg)

        if self._if_exists == OVERWRITE:
            exists = False
        return exists

    def write(self, data):
        # over-ride in subclass
        raise NotImplementedError
