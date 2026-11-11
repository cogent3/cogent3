"""Metadata attached to sequence collections and their storage."""

from collections.abc import Mapping
from types import MappingProxyType
from typing import Any

_NO_METADATA: Mapping[str, Any] = MappingProxyType({})
