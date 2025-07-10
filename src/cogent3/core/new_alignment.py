from cogent3.core.alignment import *  # noqa: F403
from cogent3.util.warning import deprecated as _deprecate

_deprecate(
    "module",
    "cogent3.core.new_alignment",
    "cogent3.core.alignment",
    "2025.9",
    reason="has been renamed",
    stack_level=3,
)
