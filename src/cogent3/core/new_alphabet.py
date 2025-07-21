from cogent3.core.alphabet import *  # noqa: F403
from cogent3.util.warning import deprecated as _deprecate

_deprecate(
    "module",
    "cogent3.core.new_alphabet",
    "cogent3.core.alphabet",
    "2025.9",
    reason="has been renamed",
    stack_level=3,
)
