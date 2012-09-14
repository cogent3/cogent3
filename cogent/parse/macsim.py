#!/usr/bin/env python

from cogent.core import annotation, moltype

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Raymond Sammut", "Peter Maxwell", "Gavin Huttley",
                    "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

#<?xml version="1.0"?>
#<!DOCTYPE macsim SYSTEM "http://www-bio3d-igbmc.u-strasbg.fr/macsim.dtd">

# As used by BAliBASE

def MacsimParser(doc):
    doc = doc.getElementsByTagName('macsim')[0]
    align = doc.getElementsByTagName('alignment')[0]
    for record in align.getElementsByTagName('sequence'):
        name = record.getElementsByTagName(
                        'seq-name')[0].childNodes[0].nodeValue
        raw_seq = record.getElementsByTagName(
                        'seq-data')[0].childNodes[0].nodeValue
        
        #cast as string to de-unicode
        raw_string = ''.join(str(raw_seq).upper().split())
        name=str(name).strip()
        
        if str(record.getAttribute('seq-type')).lower() == 'protein':
            alphabet = moltype.PROTEIN
        else:
            alphabet = moltype.DNA

        seq = alphabet.makeSequence(raw_string, Name=name)
                
        yield (name, seq)
        
