#!/usr/bin/env python
"""Retrieves records by id from RFAM, the RNA families database."""
from cogent.db.util import UrlGetter

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

rfam_base='http://www.sanger.ac.uk/cgi-bin/Rfam/getalignment.pl?'

rfam_formats = dict.fromkeys('link mul stock fal msf download belvu jalview'.split())
rfam_types = dict.fromkeys(['seed','full'])

class Rfam(UrlGetter):
    """Returns a pdb file."""
    Defaults={'type':'seed','format':'stock', 'acc':None}
    PrintedFields=dict.fromkeys('acc type format'.split())
    BaseUrl = rfam_base

    def __getitem__(self, item):
        """Returns handle to file containing aln of specified rfam id."""
        orig_acc = self.acc
        item = str(item).upper()
        if not item.startswith('RF'):
            item = 'RF'+item
        
        self.acc = item
        result = self.open()
        self.acc = orig_acc
        return result
