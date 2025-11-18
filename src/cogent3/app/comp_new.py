from __future__ import annotations

import inspect
import json
import pickle
import re
import sys
import textwrap
import time
import traceback
import types
from collections.abc import Generator, Iterable
from enum import Enum
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    ClassVar,
    Concatenate,
    Generic,
    Literal,
    ParamSpec,
    TypeVar,
    cast,
    get_type_hints,
    overload,
)
from uuid import UUID, uuid4

from scitrack import CachingLogger
from typeguard import check_type
from typing_extensions import TypeIs

from cogent3._version import __version__
from cogent3.app import typing as c3_typing
from cogent3.app.data_store import get_data_source
from cogent3.app.source_proxy import propagate_source, proxy_input, source_proxy
from cogent3.util import parallel as PAR
from cogent3.util import progress_display as UI
from cogent3.util.deserialise import register_deserialiser
from cogent3.util.misc import docstring_to_summary_rest, get_object_provenance

from .data_store import DataMemberABC, DataStoreABC, get_unique_id

if TYPE_CHECKING:
    from collections.abc import Callable

T = TypeVar("T")
U = TypeVar("U")

P = ParamSpec("P")
Q = ParamSpec("Q")

R = TypeVar("R")
S = TypeVar("S")


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


class _AppBaseClass(Generic[T, P, R]):
    fn: Callable[Concatenate[T, P], R]
    app_type: ClassVar[AppType]
    original_name: ClassVar[str]
    is_composed: ClassVar[bool] = False
    skip_not_completed: bool = True

    def __init__(self, *args: P.args, **kwargs: P.kwargs) -> None:
        self.args = args
        self.kwargs = kwargs
        self.type_check: bool = True
        self.raise_errors: bool = False
        self.write_data_store = None
        if self.app_type is WRITER:
            for arg in self.args:
                if isinstance(arg, DataStoreABC):
                    self.write_data_store = arg
                    break
            if self.write_data_store is None:
                for kwarg in self.kwargs.values():
                    if isinstance(kwarg, DataStoreABC):
                        self.write_data_store = kwarg
                        break
            if self.write_data_store is None:
                msg = "A writer app must initialise a writable datastore."
                raise ValueError(msg)

    def set_type_check(self, type_check: bool) -> None:
        self.type_check = type_check

    def set_raise_errors(self, raise_errors: bool) -> None:
        self.raise_errors = raise_errors

    def __call__(self, value: T | NotCompleted) -> R | NotCompleted:
        print("Calling")
        print("The app is", self.fn.__name__)
        if isinstance(value, NotCompleted) and self.skip_not_completed:
            return value
        cls = type(self)

        if self.type_check:
            type_hints = get_type_hints(cls.fn)
            if not isinstance(value, NotCompleted):
                check_type(value, type_hints[next(iter(type_hints))])
        else:
            type_hints = {}
        print("Doing function")

        try:
            result = cls.fn(cast("T", value), *self.args, **self.kwargs)
        except Exception as e:
            if isinstance(e, _AppInstantiationError):
                raise e.original.with_traceback(e.tb) from None
            if self.raise_errors:
                raise
            result = NotCompleted("ERROR", self, traceback.format_exc(), source=value)

        print("FINISHED CALL, GOT", type(result))

        if self.type_check and not isinstance(result, NotCompleted):
            print("TYPE IS", type(result), "EXPECTED", type_hints["return"])
            check_type(result, type_hints["return"])
        return result

    def main(self, value: T) -> R:
        raise_errors = self.raise_errors
        self.raise_errors = True
        result = self(value)
        self.raise_errors = raise_errors
        return cast("R", result)

    @UI.display_wrap
    def as_completed(
        self,
        dstore: str
        | Path
        | Iterable[str | Path]
        | Iterable[DataMemberABC]
        | DataStoreABC,
        parallel: bool = False,
        par_kw: dict[str, Any] | None = None,
        id_from_source: Callable[[object], str | None] = get_unique_id,
        **kwargs: Any,
    ) -> Generator[source_proxy[R] | R | NotCompleted, None, None]:
        """invokes self composable function on the provided data store

        Parameters
        ----------
        dstore
            a path, list of paths, or DataStore to which the process will be
            applied.
        parallel : bool
            run in parallel, according to arguments in par_kw. If True,
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
        app = propagate_source(self, id_from_source)

        ui: UI.ProgressContext = kwargs.pop("ui")

        if isinstance(dstore, (str, Path)):
            dstore = [dstore]
        elif isinstance(dstore, DataStoreABC):
            dstore = dstore.completed

        mapped = cast("list[source_proxy[T]]", proxy_input(dstore))
        if not mapped:
            yield from []
            return
        if parallel:
            par_kw = par_kw or {}
            to_do = PAR.as_completed(app, mapped, **par_kw)
        else:
            to_do = map(app, mapped)

        print("TRYING TO PICKLE LOGGER")
        for attr in dir(self):
            try:
                pickle.dumps(getattr(self, attr))
            except Exception as e:
                print("FAILED", attr, e)
            else:
                print("SUCCESS", attr)
        print(self.is_composed)
        print(self.__class__.__name__, self.__class__.__module__)
        print(pickle.dumps(self.app_type))
        print("TRYING TO PICKLE APP")
        print(pickle.dumps(self))
        print("TRYING TO UNPICKLE APP")
        print(pickle.loads(pickle.dumps(self)))
        print("TRYING TO PICKLE MAPPED")
        print(pickle.dumps(mapped))
        print("TRYING TO UNPICKLE MAPPED")
        print(pickle.loads(pickle.dumps(mapped)))
        yield from ui.series(to_do, count=len(mapped), **kwargs)

    def apply_to(
        self,
        dstore: str
        | Path
        | Iterable[str | Path]
        | Iterable[DataMemberABC]
        | DataStoreABC,
        id_from_source: Callable[[object], str | None] = get_unique_id,
        parallel: bool = False,
        par_kw: dict[str, Any] | None = None,
        logger: CachingLogger | None = None,
        cleanup: bool = True,
        show_progress: bool = False,
    ) -> DataStoreABC | None:
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
        print("PICKLING SELF 1")
        print(pickle.dumps(self))
        app = propagate_source(self, id_from_source)

        print("PICKLING SELF 2")
        print(pickle.dumps(self))

        if isinstance(dstore, str | Path):  # one filename
            dstore = [dstore]
        elif isinstance(dstore, DataStoreABC):
            dstore = dstore.completed
        print("PICKLING SELF 3")
        print(pickle.dumps(self))

        # TODO this should fail if somebody provides data that cannot produce a unique_id
        inputs: dict[str, str | Path | DataMemberABC] = {}
        for m in dstore:
            input_id = Path(m.unique_id) if isinstance(m, DataMemberABC) else m
            input_id = id_from_source(input_id)
            if input_id in inputs or not input_id:
                msg = f"non-unique identifier {input_id!r} detected in data"
                raise ValueError(msg)
            # TODO: This used to have `if input_id in self.data_store: continue`... investigate

            inputs[input_id] = m
        print("PICKLING SELF 4")
        print(pickle.dumps(self))

        if (
            not dstore
        ):  # this should just return datastore, because if all jobs are done!
            msg = "dstore is empty"
            raise ValueError(msg)
        print("PICKLING SELF 5", logger)
        print(pickle.dumps(self))

        self.set_logger(logger)
        self.logger = None
        print("PICKLING SELF 6")
        print(pickle.dumps(self))

        start = 0
        if self.logger:
            start = time.time()
            logger = self.logger
            logger.log_message(str(self), label="composable function")
            logger.log_versions(["cogent3"])

        print("PICKLING SELF 7")
        print(pickle.dumps(self))
        for result in self.as_completed(
            inputs,
            parallel=parallel,
            par_kw=par_kw,
            show_progress=show_progress,
        ):
            member = cast(
                "DataMemberABC",
                self.main(
                    cast("source_proxy[T]", result).obj
                    if isinstance(result, source_proxy)
                    else cast("T", result)
                ),
            )
            if self.logger:
                md5 = getattr(member, "md5", None)
                logger.log_message(str(member), label="output")
                if md5:
                    logger.log_message(md5, label="output md5sum")

        if self.logger:
            taken = time.time() - start
            logger.log_message(f"{taken}", label="TIME TAKEN")
            log_file_path = Path(logger.log_file_path)
            logger.shutdown()
            if self.write_data_store:
                self.write_data_store.write_log(
                    unique_id=log_file_path.name,
                    data=log_file_path.read_text(),
                )
            if cleanup:
                log_file_path.unlink(missing_ok=True)

        return self.write_data_store

    def set_logger(self, logger: CachingLogger | Literal[False] | None = None) -> None:
        if logger is False:
            self.logger = None
            return
        if logger is None:
            logger = CachingLogger(create_dir=True)
        if not isinstance(logger, CachingLogger):
            msg = f"logger must be of type CachingLogger not {type(logger)}"
            raise TypeError(msg)
        if not logger.log_file_path:
            if self.write_data_store is None:
                msg = "Without a write_data_store, a directory for the logger must be specified or create_dir set to True."
                raise ValueError(msg)
            src = Path(self.write_data_store.source).parent
            logger.log_file_path = str(src / _make_logfile_name(self))
        self.logger = logger

    def _evaluate_string_annotations(self) -> None:
        cls = type(self)
        resolved_annotations = get_type_hints(cls.fn)
        for ann_name, ann_value in cls.fn.__annotations__.items():
            if isinstance(ann_value, str):
                cls.fn.__annotations__[ann_name] = resolved_annotations[ann_name]

    def __add__(
        self,
        other: _AppBaseClass[R, Q, S],
    ) -> _AppBaseClass[T, P, S]:
        print("COMPOSING", self.fn.__name__, "AND", other.fn.__name__)
        if self.write_data_store is not None and other.write_data_store is not None:
            msg = "Composed app may only contain a singular write_data_store"
            raise ValueError(msg)
        write_data_store = (
            other.write_data_store
            if self.write_data_store is None
            else self.write_data_store
        )

        if self.type_check:
            self._evaluate_string_annotations()
        if other.type_check:
            other._evaluate_string_annotations()

        def new_fn(value: T, *args: P.args, **kwargs: P.kwargs) -> S | NotCompleted:
            intermediate = type(self).fn(value, *args, **kwargs)
            if isinstance(intermediate, NotCompleted) and type(self).skip_not_completed:
                return intermediate
            # TODO: Consider adding intermediate type checking
            return type(other).fn(intermediate, *other.args, **other.kwargs)

        # Copy to avoid modifying original type hints
        new_fn.__annotations__ = type(self).fn.__annotations__.copy()
        new_fn.__annotations__["return"] = type(other).fn.__annotations__["return"]

        new_signature = inspect.signature(type(self).fn)
        new_signature = new_signature.replace(
            return_annotation=inspect.signature(type(other).fn).return_annotation,
        )

        new_fn.__signature__ = new_signature  # type: ignore[attr-defined]

        # Get the original names of the apps to be composed
        left_name = type(self).__name__.removeprefix("composed_")
        if type(self).is_composed:
            left_name = "_".join(left_name.split("_")[:-1])  # Remove uuid
        right_name = type(other).__name__.removeprefix("composed_")
        if type(other).is_composed:
            right_name = "_".join(right_name.split("_")[:-1])
        # print(type(self).is_composed)
        print("LEFT IS", left_name)
        print("RIGHT IS", right_name)
        cls_name = f"composed_{left_name}_{right_name}_{uuid4()}"
        new_fn.__name__ = f"fn_{cls_name}"

        new_original_name = (
            f"{self!r} + " if other.app_type is not LOADER else ""
        ) + repr(other)
        cls = types.new_class(
            cls_name,
            bases=(_AppBaseClass[T, P, S],),
            exec_body=lambda ns: ns.update(
                {
                    "fn": new_fn,
                    "app_type": self.app_type,
                    "original_name": new_original_name,
                    "is_composed": True,
                    "skip_not_completed": self.skip_not_completed
                    and other.skip_not_completed,
                }
            ),
        )

        cls.__module__ = new_fn.__module__
        globals()[cls_name] = cls
        composed_app: _AppBaseClass[T, P, S] = cls(*self.args, **self.kwargs)
        composed_app.set_type_check(self.type_check and other.type_check)
        composed_app.write_data_store = write_data_store
        print("FINISHED COMPOSITION")
        print(self.write_data_store, other.write_data_store, write_data_store)
        return composed_app

    def __radd__(
        self,
        other: _AppBaseClass[U, Q, T],
    ) -> _AppBaseClass[U, Q, R]:
        return other.__add__(self)

    def __repr__(self) -> str:
        if not self.is_composed:
            # Skip the first argument as that is what the app is applied to
            sig = inspect.Signature(
                list(inspect.signature(type(self).fn).parameters.values())[1:]
            )
            bound = sig.bind_partial(*self.args, **self.kwargs)
            bound.apply_defaults()

            parts: list[str] = []
            # Skip the first as that is what the app is applied to
            for name, value in bound.arguments.items():
                # TODO: Consider specially handling *args and **kwargs
                parts.append(f"{name}={value!r}")

            parameter_data = ", ".join(parts)
            data = f"{self.original_name}({parameter_data})"
        else:
            # Parameters already handled at composition time
            data = self.original_name

        return textwrap.fill(
            data,
            width=80,
            break_long_words=False,
            break_on_hyphens=False,
        )


def _make_logfile_name(process: _AppBaseClass[Any, Any, Any]) -> str:
    text = re.split(r"\s+\+\s+", str(process))
    parts: list[str] = []
    for part in text:
        if part.find("(") >= 0:
            part = part[: part.find("(")]
        parts.append(part)
    result = "-".join(parts)
    uid = str(uuid4())
    return f"{result}-{uid[:8]}.log"


def _get_origin(origin: str | object) -> str:
    return origin if isinstance(origin, str) else origin.__class__.__name__


class NotCompleted(int):
    """results that failed to complete"""

    def __new__(
        cls,
        type: str,
        origin: str | object,
        message: str,
        source: str | object | None = None,
    ):
        """
        Parameters
        ----------
        type : str
            examples are 'ERROR', 'FAIL'
        origin
            where the instance was created, can be an instance
        message : str
            descriptive message, succinct traceback
        source : str or instance with .source or .info.source attributes
            the data operated on that led to this result.
        """
        # TODO this approach to caching persistent arguments for reconstruction
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

    def __init__(
        self,
        type: str,
        origin: str | object,
        message: str,
        source: str | object | None = None,
    ) -> None:
        self.type: str
        self.origin: str
        self.message: str
        self.source: str | object | None

    def __getnewargs_ex__(self, *args, **kw) -> tuple[Any, Any]:
        return self._persistent[0], self._persistent[1]

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        name = self.__class__.__name__
        source = self.source or "Unknown"
        return f'{name}(type={self.type}, origin={self.origin}, source="{source}", message="{self.message}")'

    def to_rich_dict(self) -> dict[str, Any]:
        """returns components for to_json"""
        return {
            "type": get_object_provenance(self),
            "not_completed_construction": {
                "args": self._persistent[0],
                "kwargs": self._persistent[1],
            },
            "version": __version__,
        }

    def to_json(self) -> str:
        """returns json string"""
        return json.dumps(self.to_rich_dict())


@register_deserialiser(get_object_provenance(NotCompleted))
def deserialise_not_completed(data: dict[str, Any]) -> NotCompleted:
    """deserialising NotCompletedResult"""
    data.pop("version", None)
    init = data.pop("not_completed_construction")
    args = init.pop("args")
    kwargs = init.pop("kwargs")
    return NotCompleted(*args, **kwargs)


class _AppInstantiationError(Exception):
    def __init__(self, original_exc: Exception, tb: types.TracebackType | None) -> None:
        super().__init__("error instantiating class")
        self.original = original_exc
        self.tb = tb


@overload
def define_app(
    fn: Callable[Concatenate[T, P], R],
    *,
    app_type: AppType,
    skip_not_completed: bool,
) -> type[_AppBaseClass[T, P, R]]: ...
@overload
def define_app(
    fn: None, *, app_type: AppType, skip_not_completed: bool
) -> Callable[[Callable[Concatenate[T, P], R]], type[_AppBaseClass[T, P, R]]]: ...
@overload
def define_app(
    *, app_type: AppType, skip_not_completed: bool
) -> Callable[[Callable[Concatenate[T, P], R]], type[_AppBaseClass[T, P, R]]]: ...
@overload
def define_app(
    *, skip_not_completed: bool
) -> Callable[[Callable[Concatenate[T, P], R]], type[_AppBaseClass[T, P, R]]]: ...


def define_app(
    fn: Callable[Concatenate[T, P], R] | None = None,
    *,
    app_type: AppType = GENERIC,
    skip_not_completed: bool = True,
) -> (
    type[_AppBaseClass[T, P, R]]
    | Callable[[Callable[Concatenate[T, P], R]], type[_AppBaseClass[T, P, R]]]
):
    def make_app(
        fn: Callable[Concatenate[T, P], R], *, app_type: AppType = app_type
    ) -> type[_AppBaseClass[T, P, R]]:
        if is_app(fn):
            msg = f"{fn.__name__} is already an app."
            raise TypeError(msg)

        original_name = fn.__name__
        if isinstance(fn, type):
            tmp_cls = fn

            def new_fn(value: T, *args: P.args, **kwargs: P.kwargs) -> R:
                try:
                    instance: Callable[[T], R] = tmp_cls(*args, **kwargs)
                except Exception as e:
                    tb = sys.exc_info()[2]
                    raise _AppInstantiationError(e, tb) from e
                return instance.main(value)

            # Get the original class globals and attach it to the created function
            new_globals = sys.modules[tmp_cls.__module__].__dict__
            new_globals["_AppBaseClass"] = _AppBaseClass
            new_globals["_AppInstantiationError"] = _AppInstantiationError
            new_globals["sys"] = sys

            new_fn = types.FunctionType(
                new_fn.__code__,
                new_globals,
                name=new_fn.__name__,
                argdefs=new_fn.__defaults__,
                closure=new_fn.__closure__,
            )

            new_fn.__name__ = fn.__name__
            new_fn.__module__ = fn.__module__
            new_fn.__doc__ = fn.__doc__

            globals()[new_fn.__name__] = new_fn

            new_annotations = cast(
                "dict[str, Any]", tmp_cls.main.__annotations__
            ).copy()

            init_annotations = getattr(tmp_cls.__init__, "__annotations__", {})

            # Avoid overriding return annotation for main
            init_annotations.pop("return", None)

            new_annotations.update(init_annotations)
            if "return" in new_annotations:  # Reorder dictionary
                new_annotations["return"] = new_annotations.pop("return")
            new_fn.__annotations__ = new_annotations

            init_sig = inspect.signature(tmp_cls.__init__)
            main_sig = inspect.signature(tmp_cls.main)

            new_fn.__signature__ = inspect.Signature(
                [
                    *tuple(main_sig.parameters.values())[1:],
                    *tuple(init_sig.parameters.values())[1:],
                ],
                return_annotation=main_sig.return_annotation,
            )

            fn = new_fn

        cls_name = fn.__name__
        cls = types.new_class(
            cls_name,
            bases=(_AppBaseClass[T, P, R],),
            exec_body=lambda ns: ns.update(
                {
                    "fn": fn,
                    "app_type": app_type,
                    "original_name": original_name,
                    "skip_not_completed": skip_not_completed,
                }
            ),
        )
        cls.__module__ = fn.__module__
        try:
            mod = sys.modules.get(fn.__module__)
            if mod is not None:
                setattr(mod, cls_name, cls)
                cls.__qualname__ = cls_name
        except Exception:
            pass

        summary, body = docstring_to_summary_rest(fn.__doc__)
        cls.__doc__ = summary
        cls.__init__.__doc__ = body
        return cls

    return make_app if fn is None else make_app(fn)


def is_app_composable(obj: object) -> bool:
    """checks whether obj has been decorated by define_app and it's app_type attribute is not NON_COMPOSABLE"""
    return is_app(obj) and True  # obj.app_type is not NON_COMPOSABLE


def is_app(
    obj: object | type,
) -> TypeIs[_AppBaseClass[Any, Any, Any]] | TypeIs[type[_AppBaseClass[Any, Any, Any]]]:
    """checks whether obj has been decorated by define_app"""
    return isinstance(obj, _AppBaseClass) or (
        isinstance(obj, type) and issubclass(obj, _AppBaseClass)
    )
