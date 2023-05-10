import inspect
import json
import re
import textwrap
import time
import traceback
import types

from copy import deepcopy
from enum import Enum
from pathlib import Path
from typing import Any, Generator, Tuple
from uuid import uuid4

from scitrack import CachingLogger

from cogent3._version import __version__
from cogent3.app.typing import get_constraint_names, type_tree
from cogent3.util import parallel as PAR
from cogent3.util import progress_display as UI
from cogent3.util.misc import docstring_to_summary_rest, get_object_provenance

from .data_store_new import DataMember, get_data_source, get_unique_id


_builtin_seqs = list, set, tuple


def _make_logfile_name(process):
    text = re.split(r"\s+\+\s+", str(process))
    parts = []
    for part in text:
        if part.find("(") >= 0:
            part = part[: part.find("(")]
        parts.append(part)
    result = "-".join(parts)
    uid = str(uuid4())
    return f"{result}-{uid[:8]}.log"


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

    # pickling induces infinite recursion on python 3.9/3.10
    # only on Windows, so implementing the following methods explicitly
    def __getstate__(self):
        return self._obj, self._src, self._uuid

    def __setstate__(self, state):
        self._obj, self._src, self._uuid = state


def _call(self, val, *args, **kwargs):
    """
    Parameters
    ----------
    val
        the primary data the app operates on. If the instance type
        does not match tha defined for the first argument of the
        app.main() method, a NotCompleted result is generated.
    args, kwargs
        other positional and keyword arguments of the app.main()
        method.

    Returns
    -------
    The return type of app.main()
    """
    if val is None:
        val = NotCompleted("ERROR", self, "unexpected input value None", source=val)

    if isinstance(val, NotCompleted) and self._skip_not_completed:
        return val

    if self.app_type is not LOADER and self.input:  # passing to connected app
        val = self.input(val, *args, **kwargs)
        if isinstance(val, NotCompleted) and self._skip_not_completed:
            return val

    type_checked = self._validate_data_type(val)
    if not type_checked:
        return type_checked

    try:
        result = self.main(val, *args, **kwargs)
    except Exception:
        result = NotCompleted("ERROR", self, traceback.format_exc(), source=val)

    if result is None:
        result = NotCompleted("BUG", self, "unexpected output value None", source=val)
    return result


def _validate_data_type(self, data):
    """checks data class name matches defined compatible types"""
    # todo when move to python 3.8 define protocol checks for the two singular types
    if isinstance(data, NotCompleted) and self._skip_not_completed:
        return data

    if not self._data_types or self._data_types & {
        "SerialisableType",
        "IdentifierType",
    }:
        return True

    if isinstance(data, source_proxy):
        data = data.obj

    if isinstance(data, _builtin_seqs):
        if len(data):
            data = next(iter(data))
        else:
            return NotCompleted("ERROR", self, message="empty data", source=data)

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
    summary, body = docstring_to_summary_rest(func.__doc__)
    func.__doc__ = None

    _class_dict = {"__init__": _init, "main": _main, "_user_func": staticmethod(func)}

    for method_name, method in _class_dict.items():
        method.__name__ = method_name
        method.__qualname__ = f"{class_name}.{method_name}"

    result = types.new_class(class_name, (), exec_body=lambda x: x.update(_class_dict))
    result.__module__ = module  # necessary for pickle support
    result._func_sig = sig
    result.__doc__ = summary
    result.__init__.__doc__ = body
    return result


def define_app(
    klass=None, *, app_type: AppType = GENERIC, skip_not_completed: bool = True
):
    """decorator for building callable apps

    Parameters
    ----------
    klass
        either a class or a function. If a function, it is converted to a
        class with the function bound as a static method.
    app_type
        what type of app, typically you just want GENERIC.
    skip_not_completed
        if True (default), NotCompleted instances are returned without being
        passed to the app.main() method.

    Notes
    -----

    Instances of ``cogent3`` apps are callable. If an exception occurs,
    the app returns a ``NotCompleted`` instance with logging information.
    Apps defined with app_type ``LOADER``, ``GENERIC`` or ``WRITER`` can be
    "composed" (summed together) to produce a single callable that
    sequentially invokes the composed apps. For example, the independent
    usage of app instances ``app1`` and ``app2`` as

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

    ``define_app`` adds a ``__call__`` method which checks an input value prior
    to passing it to ``app.main()`` as a positional argument. The data checking
    results in ``NotCompleted`` being returned immediately, unless
    ``skip_not_completed==False``. If the input value type is consistent with
    the type hint on the first argument of main it is passed to ``app.main()``.
    If it does not match, a new ``NotCompleted`` instance is returned.

    An example app definition.

    >>> from typing import Union
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
    ``NotCompleted`` instance is returned.)

    You can also decorate functions. In that case, they will be converted into
    a class with the same name as the original function. The function itself is
    bound to this new class as a ``staticmethod``, e.g.

    >>> from typing import Union
    >>> from cogent3.app.composable import define_app
    >>> from cogent3.app.typing import AlignedSeqsType, SerialisableType
    >>>
    >>> T = Union[AlignedSeqsType, SerialisableType]
    >>>
    >>> @define_app
    ... def omit_seqs(aln: AlignedSeqsType, quantile=None, gap_fraction=1, moltype="dna") -> T:
    ...     return aln.omit_bad_seqs(
    ...                             quantile=quantile,
    ...                             gap_fraction=gap_fraction,
    ...                             moltype=moltype
    ...                             )

    ``omit_seqs`` is now an app, allowing creating different variants which
    can be composed as per ones defined via a class.

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

        excludes = []
        if not composable:
            excludes = ["__add__", "disconnect"]
        if app_type is not WRITER:
            excludes.extend(["apply_to", "set_logger"])
        method_list = [item for item in __mapping if item not in excludes] + ["__str__"]
        # check if user defined input for composable
        if composable and getattr(klass, "input", None):
            raise TypeError(
                f"remove 'input' attribute in {klass.__name__!r}, this functionality provided by define_app"
            )

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
        setattr(klass, "_skip_not_completed", skip_not_completed)

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
        if not e:
            continue
        if not isinstance(e, source_proxy):
            e = source_proxy(e)
        inputs.append(e)

    return inputs


def _source_wrapped(self, value: source_proxy) -> source_proxy:
    value.set_obj(self(value.obj))
    return value


@UI.display_wrap
def _as_completed(self, dstore, parallel=False, par_kw=None, **kwargs) -> Generator:
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
    kwargs
        setting a show_progress boolean keyword value here
        affects progress display code, other arguments are passed to
        the cogent3.util.progress_bar.display_wrap decorator

    Notes
    -----
    If run in parallel, this instance serves as the master object and
    aggregates results. If run in serial, results are returned in the
    same order as provided.
    """
    ui = kwargs.pop("ui")
    app = (
        self.input._source_wrapped if self.app_type is WRITER else self._source_wrapped
    )

    if isinstance(dstore, str):
        dstore = [dstore]
    mapped = _proxy_input(dstore)
    if not mapped:
        raise ValueError("dstore is empty")

    if parallel:
        par_kw = par_kw or {}
        to_do = PAR.as_completed(app, mapped, **par_kw)
    else:
        to_do = map(app, mapped)

    return ui.series(to_do, count=len(mapped), **kwargs)


def is_composable(obj):
    """checks whether obj has been registered by the composable decorator"""
    return __app_registry.get(get_object_provenance(obj), False)


def _apply_to(
    self,
    dstore,
    id_from_source: callable = get_unique_id,
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
    id_from_source : callable
        makes the unique identifier from elements of dstore that will be
        used for writing results
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
    show_progress : bool
        controls progress bar display

    Returns
    -------
    The output data store instance

    Notes
    -----
    This is an append only function, meaning that if a member already exists
    in self.data_store for an input, it is skipped.

    If run in parallel, this instance spawns workers and aggregates results.
    """
    if not self.input:
        raise RuntimeError(f"{self!r} is not part of a composed function")

    if isinstance(dstore, (str, Path)):  # one filename
        dstore = [dstore]

    # todo this should fail if somebody provides data that cannot produce a unique_id
    inputs = {}
    for m in dstore:
        input_id = Path(m.unique_id) if isinstance(m, DataMember) else m
        input_id = id_from_source(input_id)
        if input_id in inputs:
            raise ValueError("non-unique identifier detected in data")
        if input_id in self.data_store:  # todo write a test
            continue
        inputs[input_id] = m

    if not dstore:  # this should just return datastore, because if all jobs are done!
        raise ValueError("dstore is empty")

    start = time.time()
    self.set_logger(logger)
    logger = self.logger
    logger.log_message(str(self), label="composable function")
    logger.log_versions(["cogent3"])

    inputs = _proxy_input(inputs.values())
    for result in self.as_completed(
        inputs, parallel=parallel, par_kw=par_kw, show_progress=show_progress
    ):
        member = self.main(data=result.obj, identifier=id_from_source(result.source))
        md5 = member.md5
        logger.log_message(str(member), label="output")
        if md5:
            logger.log_message(md5, label="output md5sum")

    taken = time.time() - start
    logger.log_message(f"{taken}", label="TIME TAKEN")
    log_file_path = Path(logger.log_file_path)
    logger.shutdown()
    self.data_store.write_log(
        unique_id=log_file_path.name, data=log_file_path.read_text()
    )
    if cleanup:
        log_file_path.unlink(missing_ok=True)

    return self.data_store


def _set_logger(self, logger=None):

    if logger is None:
        logger = CachingLogger(create_dir=True)
    if not isinstance(logger, CachingLogger):
        raise TypeError(f"logger must be of type CachingLogger not {type(logger)}")
    if not logger.log_file_path:
        src = Path(self.data_store.source).parent
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
