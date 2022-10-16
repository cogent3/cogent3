from pathlib import Path

from .composable import WRITER, NotCompleted
from .composable_new import define_app2
from .data_store_new import SKIP, DataStoreDirectory
from .typing import IdentifierType, SeqsCollectionType


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


@define_app2(app_type=WRITER)
class WriteSeqs:
    def __init__(self, data_path, if_dest_exists=SKIP, if_member_exists=SKIP):
        self.data_store = DataStoreDirectory(
            data_path, if_dest_exists, if_member_exists, suffix="fasta"
        )
        self.data_store = DataStoreDirectory(
            data_path, if_dest_exists, if_member_exists, suffix="fasta"
        )

    def main(self, data: SeqsCollectionType, identifier=None) -> IdentifierType:
        if isinstance(data, NotCompleted):
            return self.data_store.write_not_completed(identifier, data.to_json())

        return self.data_store.write(identifier, data.to_fasta())


@define_app2
def get_bytes(path: IdentifierType) -> bytes:
    path = Path(path)
    return path.read_bytes()
