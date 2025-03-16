import stevedore

# Entry_point for plugins to register themselves as apps
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
