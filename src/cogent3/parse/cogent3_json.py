"""Supports JSON load/read/write operations on major Cogent3 objects."""

import json
from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

from scinexus.deserialise import deserialise_object
from scinexus.io_util import open_
from scinexus.misc import get_object_provenance

if TYPE_CHECKING:  # pragma: no cover
    from cogent3.util.io import PathType


def load_from_json(filename: "PathType", classes: Iterable[type]) -> Any:
    """Loads objects from json files.

    Parameters
    ----------
    filename : Union[str,Path]
        name of the json file
    classes : Iterable[type]
        A series of the Cogent3 types, for example: (Alignment, SequenceCollection)
    """
    assert all(isinstance(klass, type) for klass in classes), (
        "classes should be a series of Cogent3 types, for example: (Alignment, SequenceCollection)"
    )

    from cogent3.app.data_store import load_record_from_json

    with open_(filename) as f:
        content = json.loads(f.read())
    try:
        _, data, completed = load_record_from_json(content)
        if not completed:
            msg = "json file is a record for type NotCompleted."
            raise TypeError(msg)
    except (KeyError, TypeError):
        data = content

    type_ = data.get("type", None)
    if type_ is None:
        msg = "json does not contain 'type' key"
        raise TypeError(msg)

    valid_types = {get_object_provenance(klass) for klass in classes}
    if type_ not in valid_types:
        msg = f"Invalid data type: {type_} is not one of {valid_types}"
        raise TypeError(msg)

    return deserialise_object(data)
