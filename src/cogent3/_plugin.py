import functools
import typing

import stevedore

if typing.TYPE_CHECKING:
    from cogent3.core.new_alignment import Alignment, SequenceCollection
    from cogent3.core.tree import PhyloNode
    from cogent3.evolve.fast_distance import DistanceMatrix
    from cogent3.format.sequence import SequenceWriterBase
    from cogent3.parse.sequence import SequenceParserBase

SeqsTypes = typing.Union["SequenceCollection", "Alignment"]

# Entry point for plugins to register themselves as hooks
HOOK_ENTRY_POINT = "cogent3.hook"


@functools.cache
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
    if name != "cogent3":
        for extension in mgr.extensions:
            if name is None or extension.module_name.startswith(name):
                return extension.plugin()

    return cogent3.get_app("quick_tree")


# registry for parsing sequence file formats
SEQ_PARSER_ENTRY_POINT = "cogent3.parse.sequence"


@functools.cache
def get_seq_format_parser_plugin(
    *,
    format_name: str | None = None,
    file_suffix: str | None = None,
) -> "SequenceParserBase":
    """returns sequence format parser plugin

    Parameters
    ----------
    format_name
        name of sequence format
    file_suffix
        suffix of file to parse

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
) -> "SequenceWriterBase":
    """returns sequence format writer

    Parameters
    ----------
    format_name
        name of sequence format
    file_suffix
        suffix of file to parse

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
            return plugin

    if built_in:
        # if we have a built-in plugin, return it
        return built_in

    msg = f"Unknown writer for format {format_name!r} or file suffix {file_suffix!r}"
    raise ValueError(msg)


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
