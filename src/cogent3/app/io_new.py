from cogent3.format.alignment import FORMATTERS

from .composable import WRITER, NotCompleted, define_app
from .typing import IdentifierType, SeqsCollectionType


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2022.10.31a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


@define_app(app_type=WRITER)
class write_seqs:
    def __init__(
        self,
        data_store,
        format="fasta",
    ):
        self.data_store = data_store
        self._formatter = FORMATTERS[format]

    def _make_outname(self, obj) -> str:
        try:
            identifier = obj.info.source
        except AttributeError:
            identifier = obj.source

        return identifier

    def main(self, data: SeqsCollectionType, identifier=None) -> IdentifierType:
        identifier = identifier or self._make_outname(data)
        if isinstance(data, NotCompleted):
            return self.data_store.write_not_completed(data.source, data.to_json())

        data = self._formatter(data.to_dict())
        return self.data_store.write(unique_id=identifier, data=data)
