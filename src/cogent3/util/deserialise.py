import warnings
from importlib import import_module

from scinexus.deserialise import deserialise_object, register_deserialiser  # noqa: F401

warnings.warn(
    "cogent3.util.deserialise is discontinued and will be removed in version 2026.9, "
    "use scinexus.deserialise instead",
    DeprecationWarning,
    stacklevel=2,
)


def get_class(provenance: str) -> type:  # pragma: no cover
    index = provenance.rfind(".")
    assert index > 0
    klass = provenance[index + 1 :]
    nc = "NotCompleted"
    klass = nc if nc in klass else klass
    mod = import_module(provenance[:index])
    return getattr(mod, klass)
