import functools
import typing

import stevedore

if typing.TYPE_CHECKING:
    from cogent3.core.tree import PhyloNode
    from cogent3.evolve.fast_distance import DistanceMatrix

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
