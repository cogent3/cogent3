import inspect
import json
import os
import pathlib
import re
import textwrap
import time
import traceback

from copy import deepcopy
from enum import Enum
from functools import wraps
from typing import Tuple

from scitrack import CachingLogger

from cogent3 import make_aligned_seqs, make_unaligned_seqs
from cogent3.app.typing import get_constraint_names
from cogent3.core.alignment import SequenceCollection
from cogent3.util import parallel as PAR
from cogent3.util import progress_display as UI
from cogent3.util.misc import (
    extend_docstring_from,
    get_object_provenance,
    in_jupyter,
)

from .data_store import (
    IGNORE,
    OVERWRITE,
    RAISE,
    SKIP,
    DataStoreMember,
    WritableDirectoryDataStore,
    get_data_source,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

from ..util.warning import discontinued


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
        source = get_data_source(source)
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


class ComposableType:
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


class Composable(ComposableType):
    def __init__(self, **kwargs):
        super(Composable, self).__init__(**kwargs)
        self._in = None  # input rules
        # rules operating on result but not part of a chain
        self._checkpointable = False
        self._load_checkpoint = None
        self._formatted = []

    def __str__(self):
        txt = "" if not self.input else str(self.input)
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
        if "mininterval" in kwargs:
            discontinued("argument", "mininterval", "2022.10", stack_level=1)

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


class ComposableTabular(Composable):
    def __init__(self, **kwargs):
        super(ComposableTabular, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableTabular",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class ComposableSeq(Composable):
    def __init__(self, **kwargs):
        super(ComposableSeq, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableSeq",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class ComposableAligned(Composable):
    def __init__(self, **kwargs):
        super(ComposableAligned, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableAligned",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class ComposableTree(Composable):
    def __init__(self, **kwargs):
        super(ComposableTree, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableTree",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class ComposableModel(Composable):
    def __init__(self, **kwargs):
        super(ComposableModel, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableModel",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class ComposableHypothesis(Composable):
    def __init__(self, **kwargs):
        super(ComposableHypothesis, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableHypothesis",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


class ComposableDistance(Composable):
    def __init__(self, **kwargs):
        super(ComposableDistance, self).__init__(**kwargs)
        discontinued(
            "class",
            "ComposableDistance",
            "2023.1",
            "see developer docs for new class hierarchy",
        )


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

        data_path = str(data_path)

        if data_path.endswith(".tinydb") and not self.__class__.__name__.endswith("db"):
            raise ValueError("tinydb suffix reserved for write_db")

        self._checkpointable = True
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


class user_function(Composable):
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
        txt = "" if not self.input else str(self.input)
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
def appify(input_types, output_types, data_types=None):
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


### Aliases to use Enum easily
LOADER = AppType.LOADER
WRITER = AppType.WRITER
GENERIC = AppType.GENERIC


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

    _no_value = inspect.Parameter.empty

    params = inspect.signature(main_func)
    if len(params.parameters) < 2:
        raise ValueError("main() method must have at least 1 input parameter")
    # list of annotation of all parameters other than self, params.parameters is orderedDict
    # annotation for first parameter other than self
    first_param_type = [p.annotation for p in params.parameters.values()][1]
    return_type = params.return_annotation
    ## Is TypeError correct?
    if return_type is _no_value:
        raise TypeError("main() method must have valid type hint for return type")
    if first_param_type is _no_value:
        raise TypeError("main() method must have valid type hint for first parameter")

    if first_param_type is None:
        raise TypeError("main() method must not have NoneType first parameter")
    if return_type is None:
        raise TypeError("main() method must not have NoneType return value")

    return frozenset([first_param_type]), frozenset([return_type])


class _ser:
    __slots__ = ("kwargs",)

    def __call__(self):
        return self.kwargs


def _add(self, other):
    ## Check order
    if self.app_type is not LOADER:
        raise TypeError("Left hand side of add operator must be of type loader")
    if other.app_type is LOADER:
        raise TypeError("Left hand side of add operator must not be of type loader")

    ## validate that self._return_types & other._input_types is a non-empty set.
    if not self._return_types:
        raise TypeError(
            "Left hand side of add operator must have non-empty return types"
        )
    if not other._data_types:  ## I am not sure about this
        raise TypeError(
            "right hand side of add operator must have non-empty input types"
        )

    ### Check if self._return_types & other._input_types is incompatible.
    if self._return_types != other._data_types:
        raise TypeError("two apps does not have compatible types")

    other.input = self
    return other


def _repr(self):
    val = f"{self.input!r} + " if self.input else ""
    data = ", ".join(f"{k}={v!r}" for k, v in self._serialisable().items())
    data = f"{self.__class__.__name__}({val}{data})"
    data = textwrap.fill(data, width=80, break_long_words=False, break_on_hyphens=False)
    return data


def _new(klass, *args, **kwargs):
    obj = object.__new__(klass)

    init_sig = inspect.signature(klass.__init__)
    no_default = inspect.Parameter.empty
    arg_order = list(init_sig.parameters)[1:]  # exclude "self"

    kw_args = {
        k: p.default
        for k, p in init_sig.parameters.items()
        if k != "self" and (p.default is not no_default)
    }

    init_vals = {}
    if arg_order:  ### this line is new, for empty args
        for i, v in enumerate(args):
            k, v = arg_order[i], v
            init_vals[k] = v

    init_vals |= kw_args | kwargs
    obj._serialisable.kwargs = init_vals
    return obj


def _call(self, val, *args, **kwargs):
    # logic for trapping call to main

    if val is None:
        ## new message in place of traceback
        return NotCompleted("ERROR", self, "unexpected input value None", source=val)

    if not val:
        return val

    if self.app_type is not LOADER and self.input:  # passing to connected app
        val = self.input(val, *args, **kwargs)
        if not val:  ### moved this if into above if
            return val

    type_checked = self._validate_data_type(val)
    if not type_checked:
        return type_checked

    try:
        result = self.main(val, *args, **kwargs)
    except Exception:
        result = NotCompleted("ERROR", self, traceback.format_exc(), source=val)

    if not result and not isinstance(result, NotCompleted):
        msg = (
            f"The value {result!r} equates to False. "
            "If unexpected, please post this error message along"
            f" with the code and data {val!r} as an Issue on the"
            " github project page."
        )
        result = NotCompleted("BUG", f"{self.__class__.__name__}", msg, source=val)

    return result


def _getstate(self):
    data = self._serialisable()
    return deepcopy(data)


def _setstate(self, data):
    return _new(self.__class__, **data)


def _validate_data_type(self, data):
    """checks data class name matches defined compatible types"""
    # if not self._data_types:
    #     # not defined
    #     return True
    class_name = data.__class__.__name__
    data_types = get_constraint_names(next(iter(self._data_types)))
    valid = class_name in data_types
    if not valid:
        msg = f"invalid data type, '{class_name}' not in {', '.join([p.__name__ for p in self._data_types])}"
        valid = NotCompleted("ERROR", self, message=msg, source=data)
    return valid


class _connected:
    def __init__(self):
        self.storage = {}

    def __get__(self, instance, owner):
        if instance.app_type is LOADER:
            raise TypeError("Looking for input value for Loader app is not permitted.")
        return self.storage.get(id(instance), None)

    def __set__(self, instance, value):
        if instance.app_type is WRITER:
            raise TypeError("Assigning input value for Writer app is not permitted.")
        self.storage[id(instance)] = value


def composable(klass=None, *, app_type: AppType = GENERIC):
    app_type = AppType(app_type)

    def wrapped(klass):
        if not inspect.isclass(klass):
            raise ValueError(f"{klass} is not a class")

        # check if user defined these functions
        method_list = [
            "__call__",
            "__add__",
            "input",
            "__repr__",
            "__str__",
            "__new__",
            "disconnect",
            "apply_to",
        ]
        for meth in method_list:
            if inspect.isfunction(getattr(klass, meth, None)):
                raise TypeError(
                    f"remove {meth!r} method in {klass.__name__!r}, this functionality provided by composable"
                )
        if getattr(klass, "input", None):
            raise TypeError(
                f"remove 'input' attribute in {klass.__name__!r}, this functionality provided by composable"
            )

        _mapping = {
            "__new__": _new,
            "__add__": _add,
            "__call__": _call,
            "__repr__": _repr,  # str(obj) calls __repr__ if __str__ missing
            "__getstate__": _getstate,
            "__setstate__": _setstate,
            "_serialisable": _ser(),
            "_validate_data_type": _validate_data_type,
            "apply_to": _apply_to,
        }

        for meth, func in _mapping.items():
            # check if the developer implements a __getstate__ or
            # __setstate__ don't introduce our own.
            if (
                meth in ["__getstate__", "__setstate__"]
                and hasattr(klass, meth)
                and inspect.isfunction(meth)
            ):
                continue

            if inspect.isfunction(func):
                func.__name__ = meth
            else:
                func.__class__.__name__ = meth
            setattr(klass, meth, func)

        # Get and type hints of main function in klass
        arg_hints, return_hint = _get_main_hints(klass)
        setattr(klass, "_data_types", arg_hints)
        setattr(klass, "_return_types", return_hint)
        setattr(klass, "app_type", app_type)

        slot_attrs = ["_data_types", "_return_types", "input"]
        if app_type == LOADER:
            slot_attrs.remove("input")
        else:
            setattr(klass, "input", _connected())
        # If a developer has defined slots, we should extend them with our own
        # instance variables.
        if hasattr(klass, "__slots__"):
            klass.__slots__ += tuple(slot_attrs)

        return klass

    return wrapped(klass) if klass else wrapped


def _apply_to(
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
    if "mininterval" in kwargs:
        discontinued("argument", "mininterval", "2022.10", stack_level=1)

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
