from cogent3.core.sequence import *  # noqa: F403
from cogent3.util.warning import deprecated as _deprecate

_deprecate(
    "module",
    "cogent3.core.new_sequence",
    "cogent3.core.sequence",
    "2025.9",
    reason="has been renamed",
    stack_level=3,
)
