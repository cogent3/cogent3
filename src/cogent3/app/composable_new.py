import inspect
import pathlib
import time
import traceback
from typing import Any, Generator
from uuid import uuid4

from cogent3.app.composable import (
    GENERIC,
    LOADER,
    NON_COMPOSABLE,
    WRITER,
    AppType,
    NotCompleted,
    _add,
    _class_from_func,
    _disconnect,
    _get_main_hints,
    _make_logfile_name,
    _repr,
)
from cogent3.app.data_store import (
    IGNORE,
    OVERWRITE,
    RAISE,
    SKIP,
    DataStoreMember,
    WritableDirectoryDataStore,
    get_data_source,
)
from cogent3.util import parallel as PAR
from cogent3.util.misc import (
    extend_docstring_from,
    get_object_provenance,
    in_jupyter,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


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


# demo implementation of a proxy class
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
    # todo implement trap for initial identifier type
    # if a result does not have a source attribute, wrap it
    # in a lightweight proxy that captures the source and
    # also forwards all method / attribute / dunder method
    # lookups to the result bound to the proxy
    # (type checking can be made aware of this proxy)
    # or possibly just unwrap the result for all calls
    # except when the app is a WRITER or the final app
    # last is unknowable

    # 1 -- no loader, no writer
    # not being written, so don't care
    # 2 -- no loader, has writer
    # need source attribute
    # 3 -- loader, no writer
    # not being written, so don't care
    # 3 -- loader, writer
    # need source attribute

    # do we just always require / wrap with source

    # logic for trapping call to main

    if val is None:
        # new message in place of traceback
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

    # result = result if hasattr(result, "source") else source_proxy(result, source)
    # result = result if hasattr(result, "source") else source_proxy(result)
    return result


def _validate_data_type(self, data):
    """checks data class name matches defined compatible types"""
    # todo when move to python 3.8 define protocol checks for the two singular types
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


def define_app2(klass=None, *, app_type: AppType = GENERIC):
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
            excludes.append("apply_to")

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

        # Get and type hints of main function in klass
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
    # for result in to_do:
    #     result.source = mapped[result]
    #     yield result


def _apply_to(
    self,
    dstore,
    get_data_source: callable = get_data_source,
    parallel=False,
    par_kw=None,
    logger=None,
    cleanup=False,
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
    for result in self.as_completed(app, inputs, parallel=parallel, par_kw=par_kw):
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
}
