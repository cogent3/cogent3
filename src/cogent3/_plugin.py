import dataclasses
import functools
import typing

import stevedore

if typing.TYPE_CHECKING:
    from cogent3.core.new_alignment import (
        AlignedSeqsDataABC,
        Alignment,
        SeqsDataABC,
        SequenceCollection,
    )
    from cogent3.core.tree import PhyloNode
    from cogent3.evolve.fast_distance import DistanceMatrix
    from cogent3.format.sequence import SequenceWriterBase
    from cogent3.parse.sequence import SequenceParserBase

SeqsTypes = typing.Union["SequenceCollection", "Alignment"]

# Entry point for plugins to register themselves as hooks
HOOK_ENTRY_POINT = "cogent3.hook"


def get_quick_tree_hook(
    *,
    name: str | None = None,
) -> typing.Callable[["DistanceMatrix"], "PhyloNode"] | None:
    """returns app instance registered for quick_tree

    Parameters
    ----------
    name
        name of package to get the app from

    Notes
    -----
    The app must take a DistanceMatrix as input and return a PhyloNode
    """
    import cogent3

    mgr = stevedore.hook.HookManager(
        namespace=HOOK_ENTRY_POINT,
        name="quick_tree",
        invoke_on_load=False,
    )
    if name and name != "cogent3":
        for extension in mgr.extensions:
            if name is None or extension.module_name.startswith(name):
                return extension.plugin()

        msg = f"Could not find quick_tree plugin for {name!r}"
        raise ValueError(msg)

    return cogent3.get_app("quick_tree")


# registry for parsing sequence file formats
SEQ_PARSER_ENTRY_POINT = "cogent3.parse.sequence"


@functools.cache
def get_seq_format_parser_plugin(
    *,
    format_name: str | None = None,
    file_suffix: str | None = None,
    unaligned_seqs: bool = True,
) -> "SequenceParserBase":
    """returns sequence format parser plugin

    Parameters
    ----------
    format_name
        name of sequence format
    file_suffix
        suffix of file to parse
    unaligned_seqs
        whether parser is for unaligned sequences

    Notes
    -----
    We default to third-party plugins if they are available, otherwise we
    use a built-in parser.
    """

    mgr = stevedore.ExtensionManager(
        namespace=SEQ_PARSER_ENTRY_POINT,
        invoke_on_load=True,
    )
    built_in = None
    for ext in mgr.extensions:
        plugin = ext.plugin()
        if file_suffix in plugin.supported_suffixes or plugin.name == format_name:
            if ext.module_name.startswith("cogent3."):
                built_in = plugin
                continue

            supports_role = (plugin.supports_unaligned and unaligned_seqs) or (
                plugin.supports_aligned and not unaligned_seqs
            )

            if not supports_role:
                out_type = "unaligned seqs" if unaligned_seqs else "aligned seqs"
                msg = f"{plugin.name} does not support parsing {out_type}"
                raise ValueError(msg)
            return plugin

    if built_in:
        # if we have a built-in plugin, return it
        return built_in

    msg = f"Unknown parser for format {format_name!r} or file suffix {file_suffix!r}"
    raise ValueError(msg)


# registry for writing sequence file formats
SEQ_FORMAT_ENTRY_POINT = "cogent3.format.sequence"


@functools.cache
def get_seq_format_writer_plugin(
    *,
    format_name: str | None = None,
    file_suffix: str | None = None,
    unaligned_seqs: bool = True,
) -> "SequenceWriterBase":
    """returns sequence format writer

    Parameters
    ----------
    format_name
        name of sequence format
    file_suffix
        suffix of file to parse
    unaligned_seqs
        whether format is for unaligned sequences

    Notes
    -----
    We default to third-party plugins if they are available, otherwise we
    use a built-in parser.
    """

    mgr = stevedore.ExtensionManager(
        namespace=SEQ_FORMAT_ENTRY_POINT,
        invoke_on_load=True,
    )
    built_in = None
    for ext in mgr.extensions:
        plugin = ext.plugin()
        if file_suffix in plugin.supported_suffixes or plugin.name == format_name:
            if ext.module_name.startswith("cogent3."):
                built_in = plugin
                continue
            supports_role = (plugin.supports_unaligned and unaligned_seqs) or (
                plugin.supports_aligned and not unaligned_seqs
            )

            if not supports_role:
                out_type = "unaligned seqs" if unaligned_seqs else "aligned seqs"
                msg = f"{plugin.name} does not support writing {out_type}"
                raise ValueError(msg)
            return plugin

    if built_in:
        # if we have a built-in plugin, return it
        return built_in

    msg = f"Unknown writer for format {format_name!r} or file suffix {file_suffix!r}"
    raise ValueError(msg)


# sequence storage drivers
UNALIGNED_SEQ_STORAGE_ENTRY_POINT = "cogent3.storage.unaligned_seqs"


def _get_driver(namespace: str, storage_backend: str) -> stevedore.driver.DriverManager:
    try:
        mgr = stevedore.driver.DriverManager(
            namespace=namespace,
            name=storage_backend,
            invoke_on_load=False,
        )
    except stevedore.exception.NoMatches as err:
        msg = f"Invalid storage backend {storage_backend!r}"
        raise ValueError(msg) from err
    return mgr


def get_unaligned_storage_driver(
    storage_backend: str | None,
) -> typing.Optional["SeqsDataABC"]:
    """returns unaligned sequence storage driver

    Parameters
    ----------
    storage_backend
        name of storage plugin to use
    """
    if not storage_backend:
        return _STORAGE_DEFAULT.unaligned

    mgr = _get_driver(UNALIGNED_SEQ_STORAGE_ENTRY_POINT, storage_backend)
    return mgr.extensions[0].plugin if mgr.extensions else _STORAGE_DEFAULT.unaligned


ALIGNED_SEQ_STORAGE_ENTRY_POINT = "cogent3.storage.aligned_seqs"


def get_aligned_storage_driver(
    storage_backend: str | None,
) -> typing.Optional["AlignedSeqsDataABC"]:
    """returns aligned sequence storage driver

    Parameters
    ----------
    storage_backend
        name of storage plugin to use
    """
    if not storage_backend:
        return _STORAGE_DEFAULT.aligned

    mgr = _get_driver(ALIGNED_SEQ_STORAGE_ENTRY_POINT, storage_backend)
    return mgr.extensions[0].plugin if mgr.extensions else _STORAGE_DEFAULT.aligned


@dataclasses.dataclass
class DefaultStorageDrivers:
    _unaligned: type = dataclasses.field(init=False, default=None)
    _aligned: type = dataclasses.field(init=False, default=None)

    @property
    def unaligned(self) -> "SeqsDataABC":
        if self._unaligned is None:
            from cogent3.core.new_alignment import SeqsData

            self._unaligned = SeqsData

        return self._unaligned

    @unaligned.setter
    def unaligned(self, driver: "SeqsDataABC") -> None:
        self._unaligned = driver

    @property
    def aligned(self) -> "AlignedSeqsDataABC":
        if self._aligned is None:
            from cogent3.core.new_alignment import AlignedSeqsData

            self._aligned = AlignedSeqsData

        return self._aligned

    @aligned.setter
    def aligned(self, driver: "AlignedSeqsDataABC") -> None:
        self._aligned = driver


_STORAGE_DEFAULT = DefaultStorageDrivers()


def set_storage_defaults(
    *,
    unaligned_seqs: str | None = None,
    aligned_seqs: str | None = None,
    reset: bool = False,
) -> None:
    """set default values for storage of unaligned and aligned seqs data

    Parameters
    ----------
    unaligned_seqs
        name of storage backend for unaligned sequences
    aligned_seqs
        name of storage backend for aligned sequences
    reset
        resets defaults to cogent3 objects
    """
    if reset:
        _STORAGE_DEFAULT.unaligned = None
        _STORAGE_DEFAULT.aligned = None
        return

    if not any([unaligned_seqs, aligned_seqs]):
        return

    if unaligned_seqs:
        _STORAGE_DEFAULT.unaligned = get_unaligned_storage_driver(
            storage_backend=unaligned_seqs,
        )

    if aligned_seqs:
        _STORAGE_DEFAULT.aligned = get_aligned_storage_driver(
            storage_backend=aligned_seqs,
        )


# Entry point for plugins to register themselves as apps
APP_ENTRY_POINT = "cogent3.app"

# private global to hold an ExtensionManager instance for apps
# Note that this is relied on for the tests
__apps = None


def get_app_manager() -> stevedore.ExtensionManager:
    """
    Lazy load a stevedore ExtensionManager to collect apps.
    """
    global __apps  # noqa: PLW0603
    if not __apps:
        __apps = stevedore.ExtensionManager(
            namespace=APP_ENTRY_POINT,
            invoke_on_load=False,
        )

    return __apps
