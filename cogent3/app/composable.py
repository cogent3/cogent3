import os
import re
import pathlib
import inspect
import time

import scitrack

from cogent3 import LoadSeqs
from cogent3.core.alignment import SequenceCollection
from cogent3.util.misc import open_
from cogent3.util import progress_display as UI
from .data_store import (SKIP, RAISE, OVERWRITE, IGNORE,
                         WritableZippedDataStore,
                         WritableDirectoryDataStore, DataStoreMember, )

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def _make_logfile_name(process):
    text = str(process)
    text = re.split(r'\s+\+\s+', text)
    parts = []
    for part in text:
        parts.append(part[:part.find('(')])
    result = '-'.join(parts)
    pid = os.getpid()
    result = f'{result}-pid{pid}.log'
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


class NotCompletedResult(int):
    """for tracking app results that failed to complete"""

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
        d = locals()
        d.pop('cls')
        result = int.__new__(cls, False)
        args = tuple(d.pop(v) for v in ('type', 'origin', 'message'))
        result._persistent = args, d

        result.type = type
        result.origin = _get_origin(origin)
        result.message = message
        result.source = _get_source(source)
        return result

    def __getnewargs_ex__(self, *args, **kw):
        return self._persistent[0], self._persistent[1]

    def __repr__(self):
        return str(self)

    def __str__(self):
        name = self.__class__.__name__
        source = self.source or 'Unknown'
        val = (f'{name}(type={self.type}, origin={self.origin}, '
               f'source="{source}", message="{self.message}")')
        return val


class ComposableType:
    _type = None

    def __init__(self, input_type=None, output_type=None):
        """
        Parameters
        ----------
        input_type
            Allowed input types
        output_type
            Type of output
        """
        input_type = [] if input_type is None else input_type
        output_type = [] if output_type is None else output_type

        if type(input_type) == str:
            input_type = [input_type]
        if type(output_type) == str:
            output_type = [output_type]

        self._input_type = set(input_type)
        self._output_type = set(output_type)

    def compatible_input(self, other):
        result = other._output_type & self._input_type
        return result != set()


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
        txt = '' if not self.input else str(self.input)
        if txt:
            txt += ' + '
        txt += '%s(%s)' % (self.__class__.__name__,
                           ', '.join(self._formatted))
        return txt

    def __repr__(self):
        return str(self)

    def _formatted_params(self):
        stack = inspect.stack()
        stack.reverse()
        for level in stack:
            args = inspect.getargvalues(level.frame).locals
            klass = args.get('__class__', None)
            if klass and isinstance(self, klass):
                break
        args = inspect.getargvalues(level.frame).locals
        params = inspect.signature(args['__class__'].__init__).parameters
        formatted = []
        for p in params:
            if p == 'self':
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
                if p == 'kwargs' and v == {'store_history': True,
                                           'silent': False}:
                    continue
            except NameError:
                pass
            formatted.append('%s=%r' % (p, v))
        self._formatted += formatted

    def __add__(self, other):
        if self.output or other.input:
            # can only be part of a single composable function
            self_name = self.__class__.__name__
            other_name = other.__class__.__name__
            if self.output and other.input:
                msg = f'{self_name} and {other_name} are already part of a' \
                    ' composed function'
            elif self.output:
                msg = f'{self_name} already part of composed function'
            else:
                msg = f'{other_name} already part of composed function'
            raise AssertionError(f'{msg}, use disconnect() to free them up')

        if not other.compatible_input(self):
            msg = '%s() requires input type "%s", %s() produces "%s"'
            my_name = self.__class__.__name__
            my_output = self._output_type
            their_name = other.__class__.__name__
            their_input = other._input_type
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

    def __call__(self, val, *args, **kwargs):
        # initial invocation always transfers call() to first composable
        # element to get input for self
        if not val:
            return val

        if self.checkpointable:
            job_done = self.job_done(val)
            if job_done and self.output:
                result = self._load_checkpoint(val)
            elif job_done:
                result = self._make_output_identifier(val)

            if job_done:
                return result

        if self.input:
            try:
                val = self._in(val, *args, **kwargs)
            except Exception as err:
                val = NotCompletedResult(
                    'ERROR', self, err.args[0], source=val)
                return val

        if not val:
            return val

        result = self.func(val, *args, **kwargs)
        if not result and type(result) != NotCompletedResult:
            msg = ('This is unexpected! Please post this error message along'
                   ' with the code and data indicated as an Issue on the'
                   ' bitbucket project page.')
            origin = str(self)
            result = NotCompletedResult('BUG', origin, msg, source=val)

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
    def apply_to(self, dstore, parallel=False, mininterval=2, par_kw=None,
                 logger=True, cleanup=False, ui=None):
        """invokes self composable function on the provided data store

        Parameters
        ----------
        dstore : DataStore
            applies composable function to this object
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
        start = time.time()
        loggable = hasattr(self, 'data_store')
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
            LOGGER.log_message(str(self), label='composable function')
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

        for result in ui.imap(process, dstore, parallel=parallel,
                              par_kw=par_kw, mininterval=mininterval):
            outcome = self(result)
            results.append(outcome)
            if LOGGER:
                member = dstore[i]
                LOGGER.log_message(member, label='input')
                LOGGER.log_message(member.md5, label='input md5sum')
                if outcome:
                    mem_id = self.data_store.make_relative_identifier(
                        member.name)
                    member = self.data_store.get_member(mem_id)
                    LOGGER.log_message(member, label='output')
                    LOGGER.log_message(member.md5, label='output md5sum')
                else:
                    # we have a NotCompletedResult
                    LOGGER.log_message(f'{outcome.origin} : {outcome.message}',
                                       label=outcome.type)

            i += 1

        finish = time.time()
        taken = finish - start
        if LOGGER:
            LOGGER.log_message(f'{taken}', label='TIME TAKEN')
            LOGGER.shutdown()
            self.data_store.add_file(str(log_file_path), cleanup=cleanup)

        # now reconnect input
        if process is not self:
            self = process + self

        return results


class ComposableTabular(Composable):
    _type = 'tabular'


class ComposableSeq(Composable):
    _type = 'sequences'


class ComposableAligned(Composable):
    _type = 'aligned'


class ComposableTree(Composable):
    _type = 'tree'


class ComposableModel(Composable):
    _type = 'model'


class ComposableHypothesis(Composable):
    _type = 'hypothesis'


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
            seqs = LoadSeqs(data=data, moltype=self.moltype,
                            aligned=self.aligned)

        if self._output_type == 'sequence':
            seqs = seqs.degap()

        return seqs


class _checkpointable(Composable):
    def __init__(self, input_type, output_type, data_path, name_callback=None,
                 create=False, if_exists=SKIP, suffix=None):
        """
        Parameters
        ----------
        input_type
            compatible input types
        output_type
            output types
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
        """
        super(_checkpointable, self).__init__(input_type=input_type,
                                              output_type=output_type)
        self._formatted_params()

        self._checkpointable = True
        if_exists = if_exists.lower()
        assert if_exists in (SKIP, IGNORE,
                             RAISE, OVERWRITE), 'invalid value for if_exists'
        self._if_exists = if_exists

        klass = (WritableZippedDataStore
                 if data_path.endswith('.zip') else WritableDirectoryDataStore)
        self.data_store = klass(data_path, suffix=suffix, create=create,
                                if_exists=if_exists)
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
