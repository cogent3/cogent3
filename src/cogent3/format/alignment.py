#!/usr/bin/env python
import os
import re

from cogent3.format.fasta import alignment_to_fasta
from cogent3.format.gde import alignment_to_gde
from cogent3.format.paml import alignment_to_paml
from cogent3.format.phylip import alignment_to_phylip
from cogent3.parse.record import FileFormatError
from cogent3.util.io import atomic_write


__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Thomas La"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

# todo convert formatters so str(formatter) returns correctly formatted
# string, and rename method names, rename base class name (sequences, not
# alignment)

_compression = re.compile(r"\.(gz|bz2)$")


def save_to_filename(alignment, filename, format, **kw):
    """Arguments:
    - alignment: to be written
    - filename: name of the sequence alignment file
    - format: the multiple sequence file format
    """
    if format is None:
        raise FileFormatError("format not known")

    with atomic_write(filename, mode="wt") as f:
        try:
            write_alignment_to_file(f, alignment, format, **kw)
        except Exception:
            try:
                os.unlink(filename)
            except Exception:
                pass
            raise


def write_alignment_to_file(f, alignment, format, **kw):
    format = format.lower()
    if format not in FORMATTERS:
        raise FileFormatError(f"Unsupported file format {format}")
    contents = FORMATTERS[format](alignment, **kw)
    f.write(contents)
    f.close()


FORMATTERS = {
    "phylip": alignment_to_phylip,
    "paml": alignment_to_paml,
    "fasta": alignment_to_fasta,
    "mfa": alignment_to_fasta,
    "fa": alignment_to_fasta,
    "gde": alignment_to_gde,
}
