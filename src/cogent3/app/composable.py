import inspect
import json
import os
import pathlib
import re
import textwrap
import time
import traceback
import types

from copy import deepcopy
from enum import Enum
from functools import wraps
from typing import Any, Generator, Tuple
from uuid import uuid4

from scitrack import CachingLogger

from cogent3.app.typing import get_constraint_names, type_tree
from cogent3.util import parallel as PAR
from cogent3.util import progress_display as UI
from cogent3.util.misc import (
    extend_docstring_from,
    get_object_provenance,
    in_jupyter,
)
from cogent3.util.warning import discontinued

from ..util.warning import discontinued
from .data_store import (
    IGNORE,
    OVERWRITE,
    RAISE,
    SKIP,
    DataStoreMember,
    WritableDirectoryDataStore,
)
from .data_store_new import get_data_source


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.10.31a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


_builtin_seqs = list, set, tuple


def _make_logfile_name(process):
    text = str(process)
    text = re.split(r"\s+\+\s+", text)
    parts = [part[: part.find("(")] for part in text]
    result = "-".join(parts)
    pid = os.getpid()
    result = f"{result}-pid{pid}.log"
    return result


def _get_origin(origin):
    return origin if type(origin) == str else origin.__class__.__name__


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
        try:
            source = get_data_source(source)
        except Exception:
            source = None
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
        return f'{name}(type={self.type}, origin={self.origin}, source="{source}", message="{self.message}")'

    def to_rich_dict(self):
        """returns components for to_json"""
        return {
            "type": get_object_provenance(self),
            "not_completed_construction": dict(
                args=self._persistent[0], kwargs=self._persistent[1]
            ),
            "version": __version__,
        }

    def to_json(self):
        """returns json string"""
        return json.dumps(self.to_rich_dict())


class ComposableType:  # pragma: no cover
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
        from cogent3.util.misc import get_object_provenance

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
        prov = get_object_provenance(self)
        if not prov.startswith("cogent3"):
            from warnings import catch_warnings, simplefilter
            from warnings import warn as _warn

            # we have a 3rd party composable
            msg = (
                "Defining composable apps by inheriting from Composable is "
                "deprecated and will be removed in release 2022.11. This "
                "will be a backwards incompatible change! We are moving to a"
                " decorator for defining composable apps. For instructions "
                f" on porting {prov!r} see"
                " https://github.com/cogent3/cogent3/wiki/composable-functions"
            )
            with catch_warnings():
                simplefilter("always")
                _warn(msg, DeprecationWarning, 1)

    def compatible_input(self, other):
        result = other._output_types & self._input_types
        return result != set()

    def _validate_data_type(self, data):
        """checks data class name matches defined compatible types"""
        if not self._data_types:
            # not defined
            return True

        name = data.__class__.__name__
        valid = name in self._data_types
        if not valid:
            msg = f"invalid data type, '{name}' not in {', '.join(self._data_types)}"
            valid = NotCompleted("ERROR", self, message=msg, source=data)
        return valid


class Composable(ComposableType):  # pragma: no cover
    def __init__(self, **kwargs):
        super(Composable, self).__init__(**kwargs)
        self._in = None  # input rules
        # rules operating on result but not part of a chain
        self._checkpointable = False
        self._load_checkpoint = None
        self._formatted = []

    def __str__(self):
        txt = str(self.input) if self.input else ""
        if txt:
            txt += " + "
        txt += f"{self.__class__.__name__}({', '.join(self._formatted)})"
        txt = textwrap.fill(
            txt, width=80, break_long_words=False, break_on_hyphens=False
        )
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

            if in_jupyter():
                if p == "kwargs" and v == {"store_history": True, "silent": False}:
                    continue
            formatted.append(f"{p}={v!r}")
        self._formatted += formatted

    def __add__(self, other):
        if other is self:
            raise ValueError("cannot add an app to itself")

        if other.input:
            # can only be part of a single composable function
            other_name = other.__class__.__name__
            msg = f"{other_name} already part of composed function"
            raise AssertionError(f"{msg}, use disconnect() to free them up")

        if not other.compatible_input(self):
            msg = '%s() requires input type "%s", %s() produces "%s"'
            my_name = self.__class__.__name__
            my_output = self._output_types
            their_name = other.__class__.__name__
            their_input = other._input_types
            msg %= (their_name, their_input, my_name, my_output)
            raise TypeError(msg)
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
            if job_done:
                result = self._make_output_identifier(refobj)
                return result

        if self.input:
            val = self._in(val, *args, **kwargs)

        if not val:
            return val
        result = self._trapped_call(self.main, val, *args, **kwargs)
        if not result and type(result) != NotCompleted:
            msg = (
                f"The value {result} equates to False. "
                "If unexpected, please post this error message along"
                f" with the code and data '{val}' as an Issue on the"
                " github project page."
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
        par_kw=None,
        logger=None,
        cleanup=False,
        ui=None,
        **kwargs,
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
            Argument ignored if not an io.writer. If a scitrack logger not provided,
            one is created with a name that defaults to the composable function names
            and the process ID.
        cleanup : bool
            after copying of log files into the data store, it is deleted
            from the original location

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
        if not dstore:
            raise ValueError("dstore is empty")

        am_writer = hasattr(self, "data_store")
        if am_writer:
            start = time.time()
            if logger is None:
                logger = CachingLogger()
            elif not isinstance(logger, CachingLogger):
                raise TypeError(
                    f"logger must be scitrack.CachingLogger, not {type(logger)}"
                )

            if not logger.log_file_path:
                src = pathlib.Path(self.data_store.source)
                logger.log_file_path = str(src.parent / _make_logfile_name(self))

            log_file_path = str(logger.log_file_path)
            logger.log_message(str(self), label="composable function")
            logger.log_versions(["cogent3"])

        results = []
        process = self.input or self
        if self.input:
            # As we will be explicitly calling the input object, we disconnect
            # the two-way interaction between input and self. This means self
            # is not called twice, and self is not unecessarily pickled during
            # parallel execution.
            process.output = None
            self.input = None

        # with a tinydb dstore, this also excludes data that failed to complete
        # todo update to consider different database backends
        # we want a dict mapping input member names to their md5 sums so these can
        # be logged
        inputs = {}
        for m in dstore:
            input_id = pathlib.Path(
                m if isinstance(m, DataStoreMember) else get_data_source(m)
            )
            suffixes = input_id.suffixes
            input_id = input_id.name.replace("".join(suffixes), "")
            inputs[input_id] = m

        if len(inputs) < len(dstore):
            diff = len(dstore) - len(inputs)
            raise ValueError(
                f"could not construct unique identifiers for {diff} records, "
                "avoid using '.' as a delimiter in names."
            )

        if parallel:
            par_kw = par_kw or {}
            to_do = PAR.as_completed(process, inputs.values(), **par_kw)
        else:
            to_do = map(process, inputs.values())

        for result in ui.series(to_do, count=len(inputs)):
            if process is not self and am_writer:
                # if result is NotCompleted, it will be written as incomplete
                # by data store backend. The outcome is just the
                # associated db identifier for tracking steps below we need to
                # know it's NotCompleted.
                # Note: we directly call .write() so NotCompleted's don't
                # get blocked from being written by __call__()
                outcome = self.main(data=result)
                if result and isinstance(outcome, DataStoreMember):
                    input_id = outcome.name
                else:
                    input_id = get_data_source(result)
                    outcome = result
                input_id = pathlib.Path(pathlib.Path(input_id))
                suffixes = input_id.suffixes
                input_id = input_id.name.replace("".join(suffixes), "")
            elif process is not self:
                outcome = self(result)
            else:
                outcome = result

            results.append(outcome)

            if am_writer:
                # now need to search for the source member
                m = inputs[input_id]
                input_md5 = getattr(m, "md5", None)
                logger.log_message(input_id, label="input")
                if input_md5:
                    logger.log_message(input_md5, label="input md5sum")

                if isinstance(outcome, NotCompleted):
                    # log error/fail details
                    logger.log_message(
                        f"{outcome.origin} : {outcome.message}", label=outcome.type
                    )
                    continue
                elif not outcome:
                    # other cases where outcome is Falsy (e.g. None)
                    logger.log_message(f"unexpected value {outcome!r}", label="FAIL")
                    continue

                logger.log_message(outcome, label="output")
                logger.log_message(outcome.md5, label="output md5sum")

        if am_writer:
            finish = time.time()
            taken = finish - start

            logger.log_message(f"{taken}", label="TIME TAKEN")
            logger.shutdown()
            self.data_store.add_file(log_file_path, cleanup=cleanup, keep_suffix=True)
            self.data_store.close()

        # now reconnect input
        if process is not self:
            self = process + self

        return results


class ComposableTabular(Composable):  # pragma: no cover
    def __init__(self, **kwargs):
        super(ComposableTabular, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableTabular",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class ComposableSeq(Composable):  # pragma: no cover
    def __init__(self, **kwargs):
        super(ComposableSeq, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableSeq",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class ComposableAligned(Composable):  # pragma: no cover
    def __init__(self, **kwargs):
        super(ComposableAligned, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableAligned",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class ComposableTree(Composable):  # pragma: no cover
    def __init__(self, **kwargs):
        super(ComposableTree, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableTree",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class ComposableModel(Composable):  # pragma: no cover
    def __init__(self, **kwargs):
        super(ComposableModel, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableModel",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class ComposableHypothesis(Composable):  # pragma: no cover
    def __init__(self, **kwargs):
        super(ComposableHypothesis, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableHypothesis",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class ComposableDistance(Composable):  # pragma: no cover
    def __init__(self, **kwargs):
        super(ComposableDistance, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableDistance",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class _checkpointable:
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
        data_path = str(data_path)

        if data_path.endswith(".tinydb") and not self.__class__.__name__.endswith("db"):
            raise ValueError("tinydb suffix reserved for write_db")

        if_exists = if_exists.lower()
        assert if_exists in (
            SKIP,
            IGNORE,
            RAISE,
            OVERWRITE,
        ), "invalid value for if_exists"
        self._if_exists = if_exists

        klass = writer_class or WritableDirectoryDataStore
        self.data_store = klass(
            data_path, suffix=suffix, create=create, if_exists=if_exists
        )
        self._callback = name_callback

        # override the following in subclasses
        self._format = None
        self._formatter = None
        self._load_checkpoint = None

    def _make_output_identifier(self, data):
        if self._callback:
            data = self._callback(data)

        return self.data_store.make_absolute_identifier(data)

    def job_done(self, data):
        identifier = self._make_output_identifier(data)
        exists = identifier in self.data_store
        if exists and self._if_exists == RAISE:
            msg = f"'{identifier}' already exists"
            raise RuntimeError(msg)

        if self._if_exists == OVERWRITE:
            exists = False
        return exists

    def main(self, data) -> DataStoreMember:
        # over-ride in subclass
        raise NotImplementedError


class user_function(Composable):  # pragma: no cover
    """wrapper class for user specified function"""

    @extend_docstring_from(ComposableType.__init__, pre=False)
    def __init__(
        self, func, input_types, output_types, *args, data_types=None, **kwargs
    ):
        """
        func : callable
            user specified function
        *args
            positional arguments to append to incoming values prior to calling
            func
        **kwargs
            keyword arguments to include when calling func

        Notes
        -----
        Known types are defined as constants in ``cogent3.app.composable``, e.g.
        ALIGNED_TYPE, SERIALISABLE_TYPE, RESULT_TYPE.

        If you create a function ``foo(arg1, arg2, kwarg1=False)``. You can
        turn this into a user function, e.g.

        >>> ufunc = user_function(foo, in_types, out_types, arg1val, kwarg1=True)

        Then

        >>> ufunc(arg2val) == foo(arg1val, arg2val, kwarg1=True)
        """
        super(user_function, self).__init__(
            input_types=input_types, output_types=output_types, data_types=data_types
        )
        self._user_func = func
        self._args = args
        self._kwargs = kwargs

    def main(self, *args, **kwargs):
        """
        Parameters
        ----------
        args
            self._args + args are passed to the user function
        kwargs
            a deep copy of self._kwargs is updated by kwargs and passed to the
            user function

        Returns
        -------
        the result of the user function
        """
        args = self._args + args
        kwargs_ = deepcopy(self._kwargs)
        kwargs_.update(kwargs)
        # the following enables a decorated user function (via @appify())
        # or directly passed user function
        func = getattr(self._user_func, "__wrapped__", self._user_func)
        return func(*args, **kwargs_)

    def __str__(self):
        txt = str(self.input) if self.input else ""
        if txt:
            txt += " + "
        txt += f"user_function(name='{self._user_func.__name__}', module='{self._user_func.__module__}')"
        txt = textwrap.fill(
            txt, width=80, break_long_words=False, break_on_hyphens=False
        )
        return txt

    def __repr__(self):
        return str(self)


@extend_docstring_from(ComposableType.__init__, pre=True)
def appify(input_types, output_types, data_types=None):  # pragma: no cover
    """function decorator for generating user apps. Simplifies creation of
    user_function() instances, e.g.

    >>> @appify(SEQUENCE_TYPE, SEQUENCE_TYPE, data_types="SequenceCollection")
    ... def omit_seqs(seqs, quantile=None, gap_fraction=1, moltype="dna"):
    ...     return seqs.omit_bad_seqs(quantile=quantile, gap_fraction=gap_fraction, moltype="dna")
    ...

    ``omit_seqs()`` is now an app factory, allowing creating variants of the app.

    >>> omit_bad = omit_seqs(quantile=0.95)

    ``omit_bad`` is now a composable ``user_function`` app. Calling with different
    args/kwargs values returns a variant app, as per the behaviour of builtin
    apps.

    """
    # the 3 nested functions are required to allow setting decorator arguments
    # to allow using functools.wraps so the decorated function has the correct
    # docstring, name etc... And, the final inner one gets to pass the
    # reference to the wrapped function (wrapped_ref) to user_function. This
    # latter is required to enable pickling of the user_function instance.
    def enclosed(func):
        @wraps(func)
        def maker(*args, **kwargs):
            # construct the user_function app
            return user_function(
                wrapped_ref,
                input_types,
                output_types,
                *args,
                data_types=data_types,
                **kwargs,
            )

        wrapped_ref = maker

        return maker

    return enclosed


class AppType(Enum):
    LOADER = "loader"
    WRITER = "writer"
    GENERIC = "generic"
    NON_COMPOSABLE = "non_composable"


# Aliases to use Enum easily
LOADER = AppType.LOADER
WRITER = AppType.WRITER
GENERIC = AppType.GENERIC
NON_COMPOSABLE = AppType.NON_COMPOSABLE


def _get_raw_hints(main_func, min_params):
    _no_value = inspect.Parameter.empty
    params = inspect.signature(main_func)
    if len(params.parameters) < min_params:
        raise ValueError(
            f"{main_func.__name__!r} must have at least {min_params} input parameters"
        )
    # annotation for first parameter other than self, params.parameters is an orderedDict
    first_param_type = [p.annotation for p in params.parameters.values()][
        min_params - 1
    ]
    return_type = params.return_annotation
    if return_type is _no_value:
        raise TypeError("must specify type hint for return type")
    if first_param_type is _no_value:
        raise TypeError("must specify type hint for first parameter")

    if first_param_type is None:
        raise TypeError("NoneType invalid type for first parameter")
    if return_type is None:
        raise TypeError("NoneType invalid type for return value")

    # we disallow type hints with too many levels of testing,
    # e.g set[int] is ok but tuple[set[int]] is not
    msg = (
        "{} type {} nesting level exceeds 2 for {}"
        "we suggest using a custom type, e.g. a dataclass"
    )
    depth, _ = type_tree(first_param_type)
    if depth > 2:
        raise TypeError(msg.format("first_param", first_param_type, depth))

    depth, _ = type_tree(return_type)
    if depth > 2:
        raise TypeError(msg.format("return_type", return_type, depth))

    return first_param_type, return_type


def _get_main_hints(klass) -> Tuple[set, set]:
    """return type hints for main method
    Returns
    -------
    {arg1hints}, {return type hint}
    """
    # Check klass.main exists and is type method
    main_func = getattr(klass, "main", None)
    if (
        main_func is None
        or not inspect.isclass(klass)
        or not inspect.isfunction(main_func)
    ):
        raise ValueError(f"must define a callable main() method in {klass.__name__!r}")

    first_param_type, return_type = _get_raw_hints(main_func, 2)
    first_param_type = get_constraint_names(first_param_type)
    return_type = get_constraint_names(return_type)

    return frozenset(first_param_type), frozenset(return_type)


def _set_hints(main_meth, first_param_type, return_type):
    """adds type hints to main"""
    main_meth.__annotations__["arg"] = first_param_type
    main_meth.__annotations__["return"] = return_type
    return main_meth


# Added new function to decorator, doesn't have function body yet
def _disconnect(self):
    """resets input to None
    Breaks all connections among members of a composed function."""
    if self.app_type is LOADER:
        return
    if self.input:
        self.input.disconnect()

    self.input = None


def _add(self, other):

    if getattr(other, "app_type", None) not in {WRITER, LOADER, GENERIC}:
        raise TypeError(f"{other!r} is not composable")

    if other.input is not None:
        raise ValueError(
            f"{other.__class__.__name__} already part of composed function, use disconnect() to free them up"
        )

    if other is self:
        raise ValueError("cannot add an app to itself")

    # Check order
    if self.app_type is WRITER:
        raise TypeError("Left hand side of add operator must not be of type writer")
    elif other.app_type is LOADER:
        raise TypeError("Right hand side of add operator must not be of type loader")

    if self._return_types & {"SerialisableType", "IdentifierType"}:
        pass
    # validate that self._return_types ia a non-empty set.
    elif not self._return_types:
        raise TypeError(f"return type not defined for {self.__class__.__name__!r}")
    # validate that other._data_types a non-empty set.
    elif not other._data_types:
        raise TypeError(f"input type not defined for {other.__class__.__name__!r}")
    # Check if self._return_types & other._data_types is incompatible.
    elif not (self._return_types & other._data_types):
        raise TypeError(
            f"{self.__class__.__name__!r} return_type {self._return_types} "
            f"incompatible with {other.__class__.__name__!r} input "
            f"type {other._data_types}"
        )
    other.input = self
    return other


def _repr(self):
    val = f"{self.input!r} + " if self.app_type is not LOADER and self.input else ""
    all_args = deepcopy(self._init_vals)
    args_items = all_args.pop("args", None)
    data = ", ".join(f"{v!r}" for v in args_items) if args_items else ""
    kwargs_items = all_args.pop("kwargs", None)
    data += (
        ", ".join(f"{k}={v!r}" for k, v in kwargs_items.items()) if kwargs_items else ""
    )
    data += ", ".join(f"{k}={v!r}" for k, v in all_args.items())
    data = f"{val}{self.__class__.__name__}({data})"
    data = textwrap.fill(data, width=80, break_long_words=False, break_on_hyphens=False)
    return data


def _new(klass, *args, **kwargs):
    obj = object.__new__(klass)

    if hasattr(klass, "_func_sig"):
        # we have a decorated function, the first parameter in the signature
        # is not given to constructor, so we create a new signature excluding that one
        params = klass._func_sig.parameters
        init_sig = inspect.Signature(parameters=list(params.values())[1:])
        bargs = init_sig.bind_partial(*args, **kwargs)
    else:
        init_sig = inspect.signature(klass.__init__)
        bargs = init_sig.bind_partial(klass, *args, **kwargs)
    bargs.apply_defaults()
    init_vals = bargs.arguments
    init_vals.pop("self", None)

    obj._init_vals = init_vals
    return obj


class source_proxy:
    __slots__ = ("_obj", "_src", "_uuid")

    def __init__(self, obj: Any = None) -> None:
        self._obj = obj
        self._src = obj
        self._uuid = uuid4()

    def __hash__(self):
        return hash(self._uuid)

    @property
    def obj(self):
        return self._obj

    def set_obj(self, obj):
        self._obj = obj

    @property
    def source(self):
        """origin of this object, defaults to a random uuid"""
        return self._src

    @source.setter
    def source(self, src: Any):
        # need to check whether src is hashable, how to cope if it isn't?
        # might need to make this instance hashable perhaps using a uuid?
        self._src = src

    def __getattr__(self, name: str) -> Any:
        return getattr(self._obj, name)

    def __setattr__(self, name: str, value: Any) -> None:
        if name.startswith("_"):
            super().__setattr__(name, value)
        else:
            setattr(self._obj, name, value)

    def __bool__(self):
        return bool(self._obj)

    def __repr__(self):
        return self.obj.__repr__()

    def __str__(self):
        return self.obj.__str__()

    def __eq__(self, other):
        return self.obj.__eq__(other)

    def __len__(self):
        return self.obj.__len__()


def _call(self, val, *args, **kwargs):
    source = get_data_source(val) or val
    if val is None:
        val = NotCompleted("ERROR", self, "unexpected input value None", source=source)

    if isinstance(val, NotCompleted):
        return val

    # todo we should get the source information from val here
    if self.app_type is not LOADER and self.input:  # passing to connected app
        val = self.input(val, *args, **kwargs)

    type_checked = self._validate_data_type(val)
    if not type_checked:
        return type_checked

    try:
        result = self.main(val, *args, **kwargs)
    except Exception:
        result = NotCompleted("ERROR", self, traceback.format_exc(), source=source)

    if result is None:
        result = NotCompleted(
            "BUG", self, "unexpected output value None", source=source
        )
    return result


def _validate_data_type(self, data):
    """checks data class name matches defined compatible types"""
    # todo when move to python 3.8 define protocol checks for the two singular types
    if isinstance(data, NotCompleted):
        return data

    if not self._data_types or self._data_types & {
        "SerialisableType",
        "IdentifierType",
    }:
        return True

    if isinstance(data, source_proxy):
        data = data.obj

    if isinstance(data, (list, tuple)) and len(data):
        data = data[0]
    elif isinstance(data, (list, tuple)):
        return NotCompleted("ERROR", self, message=f"empty data", source=data)

    class_name = data.__class__.__name__
    valid = class_name in self._data_types
    if not valid:
        msg = f"invalid data type, '{class_name}' not in {', '.join(list(self._data_types))}"
        valid = NotCompleted("ERROR", self, message=msg, source=data)
    return valid


__app_registry = {}


def _class_from_func(func):
    """make a class based on func

    Notes
    -----
    produces a class consistent with the necessary properties for
    the define_app class decorator.

    func becomes a static method on the class
    """
    # these methods MUST be in function scope so that separate instances are
    # created for each decorated function
    def _init(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs

    def _main(self, arg, *args, **kwargs):
        kw_args = deepcopy(self._kwargs)
        kw_args = {**kw_args, **kwargs}
        args = (arg,) + args + deepcopy(self._args)
        bound = self._func_sig.bind(*args, **kw_args)
        return self._user_func(**bound.arguments)

    module = func.__module__  # to be assigned to the generated class
    sig = inspect.signature(func)
    class_name = func.__name__
    _main = _set_hints(_main, *_get_raw_hints(func, 1))

    _class_dict = {"__init__": _init, "main": _main, "_user_func": staticmethod(func)}

    for method_name, method in _class_dict.items():
        method.__name__ = method_name
        method.__qualname__ = f"{class_name}.{method_name}"

    result = types.new_class(class_name, (), exec_body=lambda x: x.update(_class_dict))
    result.__module__ = module  # necessary for pickle support
    result._func_sig = sig
    return result


def define_app(klass=None, *, app_type: AppType = GENERIC):
    """decorator for building callable apps

    Parameters
    ----------
    klass
        either a class or a function. If a function, it is converted to a
        class with the function bound as a static method.
    app_type
        what type of app, typically you just want GENERIC.

    Notes
    -----

    Instances of ``cogent3`` apps are callable. If an exception occurs,
    the app returns a ``NotCompleted`` instance with logging information.
    Apps defined with app_type LOADER, GENERIC or WRITER can be "composed" (summed together)
    to produce a single callable that sequentially invokes the composed
    apps. For example, the independent usage of app instances
    ``app1`` and ``app2`` as

    >>> app2(app1(data))

    is equivalent to

    >>> combined = app1 + app2
    >>> combined(data)

    The ``app_type`` attribute is used to constrain how apps can be composed.
    ``LOADER`` and ``WRITER`` are special cases. If included, a ``LOADER``
    must always be first, e.g.

    >>> app = a_loader + a_generic

    If included, a ``WRITER`` must always be last, e.g.

    >>> app = a_generic + a_writer

    Changing the order for either of the above will result in a ``TypeError``.

    There are no constraints on ordering of ``GENERIC`` aside from compatability of
    their input and return types (see below).

    In order to be decorated with ``@define_app`` a class **must**

    - implement a method called ``main``
    - type hint the first argument of ``main``
    - type hint the return type for ``main``

    Overlap between the return type hint and first argument hint is required
    for two apps to be composed together.

    An example app definition.

    >>> from cogent3.app.composable import define_app
    >>> from cogent3.app.typing import AlignedSeqsType, SerialisableType
    ...
    ... @define_app
    ... class drop_bad:
    ...     def __init__(self, quantile=None, gap_fraction=1, moltype="dna"):
    ...         self.quantile = quantile
    ...         self.gap_fraction = gap_fraction
    ...         self.moltype = moltype
    ...
    ...     T = Union[AlignedSeqsType, SerialisableType]
    ...
    ...     def main(self, aln: AlignedSeqsType) -> T:
    ...         return aln.omit_bad_seqs(
    ...                                 quantile=self.quantile,
    ...                                 gap_fraction=self.gap_fraction,
    ...                                 moltype=self.moltype
    ...                                 )

    ``drop_bad`` is a composable app with ``app_type=GENERIC``. The input
    data must be a sequence alignment instance. It returns the same type,
    which is also serialisable.  (If invalid input data is provided a
    ``NoteCompleted`` instance is returned.)

    You can also decorate functions. In that case, they will be converted into
    a class with the same name as the original function. The function itself is
    bound to this new class as a ``staticmethod``, e.g.

    >>> from cogent3.app.composable import define_app
    >>> from cogent3.app.typing import AlignedSeqsType, SerialisableType
    >>>
    >>> T = Union[AlignedSeqsType, SerialisableType]
    >>>
    >>> @define_app
    ... def omit_seqs(aln: AlignedSeqsType, quantile=None, gap_fraction=1, moltype="dna") -> T:
    ...         return aln.omit_bad_seqs(
    ...                                 quantile=quantile,
    ...                                 gap_fraction=gap_fraction,
    ...                                 moltype=moltype
    ...                                 )

    ``omit_seqs`` is now an app, allowing creating different variants which
    can be composed as per ones defined via a class

    >>> omit_bad = omit_seqs(quantile=0.95)

    and ``omit_bad`` is an instance of that app.
    """

    if hasattr(klass, "app_type"):
        raise TypeError(
            f"The class {klass.__name__!r} is already decorated, avoid using "
            "inheritance from a decorated class."
        )

    app_type = AppType(app_type)
    composable = app_type is not NON_COMPOSABLE

    def wrapped(klass):
        if inspect.isfunction(klass):
            klass = _class_from_func(klass)

        if not inspect.isclass(klass):
            raise ValueError(f"{klass} is not a class")

        # check if user defined these methods
        # todo derive from __mapping
        method_list = [
            "__call__",
            "__repr__",
            "__str__",
            "__new__",
            "_validate_data_type",
            "as_completed",
            "_source_wrapped",
        ]
        excludes = ["__add__", "disconnect"]
        if app_type is WRITER:
            # only writers get the apply_to method
            excludes.extend(["apply_to", "set_logger"])

        if composable and getattr(klass, "input", None):
            raise TypeError(
                f"remove 'input' attribute in {klass.__name__!r}, this functionality provided by define_app"
            )
        elif composable:
            method_list.extend(excludes)

        for meth in method_list:
            # make sure method not defined by user before adding
            if inspect.isfunction(getattr(klass, meth, None)):
                raise TypeError(
                    f"remove {meth!r} in {klass.__name__!r}, this functionality provided by define_app"
                )
            func = __mapping["__repr__"] if meth == "__str__" else __mapping[meth]
            func.__name__ = meth
            setattr(klass, meth, func)

        # Get type hints of main function in klass
        arg_hints, return_hint = _get_main_hints(klass)
        setattr(klass, "_data_types", arg_hints)
        setattr(klass, "_return_types", return_hint)
        setattr(klass, "app_type", app_type)

        if app_type is not LOADER:
            setattr(klass, "input", None)

        if hasattr(klass, "__slots__"):
            # not supporting this yet
            raise NotImplementedError("slots are not currently supported")

        # register this app
        __app_registry[get_object_provenance(klass)] = composable
        return klass

    return wrapped(klass) if klass else wrapped


def _proxy_input(dstore):
    inputs = []
    for e in dstore:
        if not isinstance(e, source_proxy):
            e = source_proxy(e)
        inputs.append(e)

        # source = get_data_source(e)
        # if source:
        #     s = pathlib.Path(source)
        #     suffixes = "".join(s.suffixes)
        #     source = s.name.replace(suffixes, "")
        #
        # key = f"{source}.{suffix}" if suffix else source
        # proxied = source_proxy(e)
        # inputs[proxied] = proxied
    return inputs


def _source_wrapped(self, value: source_proxy) -> source_proxy:
    # result = result if hasattr(result, "source") else source_proxy(result)
    value.set_obj(self(value.obj))
    return value


def _as_completed(
    self,
    dstore,
    parallel=False,
    par_kw=None,
    # show_progress?
) -> Generator:
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

    Returns
    -------
    Result of the process as a list.

    Notes
    -----
    If run in parallel, this instance serves as the master object and
    aggregates results. If run in serial, results are returned in the
    same order as provided.
    """
    app = self.input if self.app_type is WRITER else self

    if isinstance(dstore, str):
        dstore = [dstore]

    mapped = _proxy_input(dstore)

    if not mapped:
        raise ValueError("dstore is empty")

    if parallel:
        par_kw = par_kw or {}
        to_do = PAR.as_completed(app._source_wrapped, mapped, **par_kw)
    else:
        to_do = map(app._source_wrapped, mapped)

    yield from to_do


def is_composable(obj):
    """checks whether obj has been registered by the composable decorator"""
    return __app_registry.get(get_object_provenance(obj), False)


def _apply_to(
    self,
    dstore,
    get_data_source: callable = get_data_source,
    parallel=False,
    par_kw=None,
    logger=None,
    cleanup=False,
    show_progress=False,
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
        Argument ignored if not an io.writer. If a scitrack logger not provided,
        one is created with a name that defaults to the composable function names
        and the process ID.
    cleanup : bool
        after copying of log files into the data store, it is deleted
        from the original location

    Returns
    -------
    Result of the process as a list

    Notes
    -----
    If run in parallel, this instance serves as the master object and
    aggregates results.
    """
    if not self.input:
        raise RuntimeError(f"{self!r} is not part of a composed function")

    if isinstance(dstore, (str, pathlib.Path)):  # one filename
        dstore = [dstore]

    # todo run check on dstore to make sure we have identifiers for writing output
    # our default function should fail if a source cannot be determined, i.e. does not return None

    # todo this should fail if somebody provides data that cannot produce a unique_id
    # skip records already in data_store
    dstore = [e for e in dstore if get_data_source(e) not in self.data_store]
    if not dstore:  # this should just return datastore, because if all jobs are done!
        raise ValueError("dstore is empty")

    start = time.time()
    self.set_logger(logger)
    logger = self.logger
    log_file_path = str(logger.log_file_path)
    logger.log_message(str(self), label="composable function")
    logger.log_versions(["cogent3"])

    inputs = _proxy_input(dstore)

    if len(inputs) < len(dstore):
        diff = len(dstore) - len(inputs)
        raise ValueError(
            f"could not construct unique identifiers for {diff} records, "
            "avoid using '.' as a delimiter in names."
        )
    app = self.input
    self.input = None
    for result in self.as_completed(inputs, parallel=parallel, par_kw=par_kw):
        member = self(data=result.obj)  # which means writers must return DataMember
        logger.log_message(member, label="output")
        logger.log_message(member.md5, label="output md5sum")

    self.input = app
    taken = time.time() - start

    logger.log_message(f"{taken}", label="TIME TAKEN")
    logger.shutdown()
    log_file_path = pathlib.Path(log_file_path)
    self.data_store.write_log(log_file_path, log_file_path.read_text())
    if cleanup:
        log_file_path.unlink(missing_ok=True)

    return self.data_store


def _set_logger(self, logger=None):
    if logger is None:
        logger = CachingLogger(create_dir=True)
    if not logger.log_file_path:
        src = pathlib.Path(self.data_store.source).parent
        logger.log_file_path = str(src / _make_logfile_name(self))
    self.logger = logger


__mapping = {
    "__new__": _new,
    "__add__": _add,
    "__call__": _call,
    "__repr__": _repr,  # str(obj) calls __repr__ if __str__ missing
    "_validate_data_type": _validate_data_type,
    "disconnect": _disconnect,
    "apply_to": _apply_to,
    "_source_wrapped": _source_wrapped,
    "as_completed": _as_completed,
    "set_logger": _set_logger,
}
