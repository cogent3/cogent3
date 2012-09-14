#/usr/bin/env python
"""Parsers for CUTG codon and spsum files.

Notes

The CUTG format specifiers are highly inaccurate. For the species sum records,
colons can appear in the organism name (e.g. for some bacterial strains), so 
it's important to split on the _last_ colon.

For the codon records, the 'nt..nt' field now stores the GenBank location line
in all its complexity. The length field precedes the PID field. There is not
really a separator between the title and the descriptions.

Some of the database _names_ contain '/', the field delimiter. For now, these
are just skipped, since writing something sufficiently general to detect 
whether they are quoted is considerably more effort than the reward warrants.
"""
from cogent.parse.record_finder import LabeledRecordFinder
from cogent.parse.record import RecordError, DelimitedSplitter
from cogent.core.info import Info, DbRef
from cogent.util.misc import caps_from_underscores as cfu
from cogent.core.usage import CodonUsage
from string import strip

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

def is_cutg_label(x):
    """Checks if x looks like a CUTG label line."""
    return x.startswith('>')

def is_cutg_species_label(x):
    """Checks if x looks like a CUTG label line."""
    return ':' in x

def is_blank(x):
    """Checks if x is blank."""
    return (not x) or x.isspace()

CutgSpeciesFinder = LabeledRecordFinder(is_cutg_species_label, ignore=is_blank)

CutgFinder = LabeledRecordFinder(is_cutg_label, ignore=is_blank)

codon_order = "CGA CGC CGG CGU AGA AGG CUA CUC CUG CUU UUA UUG UCA UCC UCG UCU AGC AGU ACA ACC ACG ACU CCA CCC CCG CCU GCA GCC GCG GCU GGA GGC GGG GGU GUA GUC GUG GUU AAA AAG AAC AAU CAA CAG CAC CAU GAA GAG GAC GAU UAC UAU UGC UGU UUC UUU AUA AUC AUU AUG UGG UAA UAG UGA".split()

#NOTE: following field order omits Locus/CDS (first field), which needs further
#processing. Use zip(field_order, fields[1:]) and handle first field specially.
field_order = "GenBank Location Length GenPept Species Description".split()

species_label_splitter = DelimitedSplitter(':', -1)

def CutgSpeciesParser(infile, strict=True, constructor=CodonUsage):
    """Yields successive sequences from infile as CodonUsage objects.

    If strict is True (default), raises RecordError when label or seq missing.
    """
    if not strict:  #easier to see logic without detailed error handling
        for rec in CutgSpeciesFinder(infile):
            try:
                label, counts = rec
                if not is_cutg_species_label(label):
                    continue
                species, genes = species_label_splitter(label)
                info = Info({'Species':species, 'NumGenes':int(genes)})
                freqs = constructor(zip(codon_order, map(int, counts.split())),
                    Info=info)
                yield freqs
            except:
                continue
    else:
        for rec in CutgSpeciesFinder(infile):
            try:
                label, counts = rec
            except ValueError:   #can't have got any counts
                raise RecordError, "Found label without sequences: %s" % rec
                
            if not is_cutg_species_label(label):
                raise RecordError, "Found CUTG record without label: %s" % rec
            species, genes = species_label_splitter(label)
            info = Info({'Species':species, 'NumGenes':int(genes)})
            try:
                d = zip(codon_order, map(int, counts.split()))
                freqs = constructor(d, Info=info)
            except:
                raise RecordError, "Unable to convert counts: %s" % counts
            yield freqs

def InfoFromLabel(line):
    """Takes a CUTG codon description line and returns an Info object.

    Raises RecordError if wrong number of fields etc.
    """
    try:
        raw_fields = line.split('\\')
        result = Info(dict(zip(field_order, map(strip, raw_fields[1:]))))
        #extra processing for first field
        first = raw_fields[0]
        if '#' in first:
            locus, cds_num = map(strip, raw_fields[0].split('#'))
        else:
            locus, cds_num = first, '1'
        result['Locus'] = locus[1:]   #remove leading '>'
        result['CdsNumber'] = cds_num
        #additional processing for last field: mostly key="value" pairs
        description = result['Description']
        descrs = description.split('/')
        for d in descrs:
            if '=' in d:    #assume key-value pair
                key, val = map(strip, d.split('=', 1))  #might be '=' in value
                #cut off leading and trailing " if present, but _not_ internal!
                if val.startswith('"'):
                    val = val[1:]
                if val.endswith('"'):
                    val = val[:-1]
                if key == 'db_xref':    #handle cross-refs specially
                    try:
                        key, val = val.split(':')
                    except ValueError:  #missing actual reference?
                        continue        #just skip the bad db records
                    try:
                        if result[key]:
                            result[key].append(val)
                        else:
                            result[key] = [val]
                    except (KeyError, TypeError):   #didn't recognize database
                        result[key] = val
                else:
                    #remember to convert the key to MixedCase naming convention
                    result[cfu(key)] = val
        return result
    except:
        raise RecordError, "Failed to read label line:\n%s" % line

def CutgParser(infile, strict=True, constructor=CodonUsage):
    """Yields successive sequences from infile as CodonUsage objects.

    If strict is True (default), raises RecordError when label or seq missing.
    """
    if not strict:  #not much error checking needed: following makes logic clear
        for rec in CutgFinder(infile):
            try:
                label, counts = rec
                if not is_cutg_label(label):
                    continue
                info = InfoFromLabel(label)
                freqs = constructor(zip(codon_order, map(int, counts.split())),
                    Info=info)
                yield freqs
            except:
                continue
    else:   #need to do more detailed error checking
        count = 0
        for rec in CutgFinder(infile):
            try:
                label, counts = rec
            except ValueError:   #can't have got any counts
                raise RecordError, "Found label without sequences: %s" % rec
            if not is_cutg_label(label):
                raise RecordError, "Found CUTG record without label: %s" % rec
            info = InfoFromLabel(label)
            try:
                freqs = constructor(zip(codon_order, map(int, counts.split())),
                    Info=info)
            except NotImplementedError:
                raise RecordError, "Unable to convert counts: %s" % counts
            yield freqs

