"""Deprecated compatibility shim for cogent3.app.composable.

This module re-exports the composable-app public surface from
``scinexus.composable``. It exists solely so that external code importing
from ``cogent3.app.composable`` keeps working until the deprecation period
ends in 2026.9. All code within cogent3 itself should import from
``scinexus.composable`` directly.
"""

import warnings

from scinexus.composable import (
    GENERIC,
    LOADER,
    NON_COMPOSABLE,
    WRITER,
    AppBase,
    AppType,
    ComposableApp,
    GetIdFuncType,
    NotCompleted,
    NotCompletedType,
    WriterApp,
    define_app,
    is_app,
    is_app_composable,
    propagate_source,
    source_proxy,
)

warnings.warn(
    "cogent3.app.composable is deprecated and will be removed in 2026.9. "
    "Import from scinexus.composable instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "GENERIC",
    "LOADER",
    "NON_COMPOSABLE",
    "WRITER",
    "AppBase",
    "AppType",
    "ComposableApp",
    "GetIdFuncType",
    "NotCompleted",
    "NotCompletedType",
    "WriterApp",
    "define_app",
    "is_app",
    "is_app_composable",
    "propagate_source",
    "source_proxy",
]
