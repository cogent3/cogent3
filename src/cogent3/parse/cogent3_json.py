"""Supports JSON load/read/write operations on major Cogent3 objects.
"""
import json

from cogent3.app.data_store import load_record_from_json
from cogent3.util.deserialise import deserialise_object
from cogent3.util.io import open_
from cogent3.util.misc import get_object_provenance


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Stephen Ma"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Development"


def load_from_json(filename, classes):
    """Loads objects from json files.

    Parameters
    ----------
    filename : Union[str,Path]
        name of the json file
    classes : Sequence[type]
        A series of the Cogent3 types, for example: (Alignment, ArrayAlignment)
    """
    assert all(
        (isinstance(klass, type) for klass in classes)
    ), "classes should be a series of Cogent3 types, for example: (Alignment, ArrayAlignment)"

    with open_(filename) as f:
        content = json.loads(f.read())
    try:
        _, data, completed = load_record_from_json(content)
        if not completed:
            raise TypeError("json file is a record for type NotCompleted.")
    except (KeyError, TypeError):
        data = content

    type_ = data.get("type", None)
    if type_ is None:
        raise TypeError("json does not contain 'type' key")

    valid_types = {get_object_provenance(klass) for klass in classes}
    if type_ not in valid_types:
        raise TypeError(f"Invalid data type: {type_} is not one of {valid_types}")

    return deserialise_object(data)
