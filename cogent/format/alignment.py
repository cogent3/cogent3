#!/usr/bin/env python
import os
import warnings
from cogent.parse.record import FileFormatError

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def save_to_filename(alignment, filename, format, **kw):
    """Arguments:
            - alignment: to be written
            - filename: name of the sequence alignment file
            - format: the multiple sequence file format
    """
    if format is None:
        raise FileFormatError("format not known")
    
    f = open(filename, 'w')
    try:
        write_alignment_to_file(f, alignment, format, **kw)
    except Exception:
        try:
            os.unlink(filename)
        except Exception:
            pass
        raise
    f.close()

def write_alignment_to_file(f, alignment, format, **kw):
    format = format.lower()
    if format not in WRITERS:
        raise FileFormatError("Unsupported file format %s" % format)
    writer = WRITERS[format](f)
    writer.writealignment(alignment, **kw)

class _AlignmentWriter(file):
    """A virtual class for writing sequence files."""
    
    def __init__(self, f):
        self.file = f
    
    # other utility function
    def setblocksize(self, size):
        """Set the length of the sequence to be printed to each line.
        
        Arguments:
            - size: the sequence length to put on each line."""
        
        self.block_size = size
    
    def setaligninfo(self, alignmentdict, order=[]):
        """Set alignment attributes for writing.
        
        Arguments:
            - alignmentdict: dictionary of seqname -> seqstring
            - order: a list of seqname's in the order for writing
        """
        
        self.number_sequences = len(alignmentdict)
        # supersede the use of alignment length
        self.align_length = len(alignmentdict[alignmentdict.keys()[0]])
        
        if order != [] and len(order) == len(alignmentdict):
            # not testing contents - possibly should.
            self.align_order = order
        else:
            self.align_order = alignmentdict.keys().sort()
    
    def slicestringinblocks(self, seqstring, altblocksize=0):
        """Return a list of string slices of specified length. No line returns.
        
        Arguments:
            - seqstring: the raw sequence string
            - altblocksize: the length of sequence for writing to each
              line, default (0) means default value specified by blocksize
              will be used.
        """
        
        if altblocksize:
            block_size = altblocksize
        else:
            block_size = self.block_size
        
        blocklist = []
        seqlength = len(seqstring)
        for block in range(0, seqlength, block_size):
            if block + block_size < seqlength:
                blocklist.append(seqstring[block: block + block_size])
            else:
                blocklist.append(seqstring[block:])
        
        return blocklist
    
    def wrapstringtoblocksize(self, seqstring, altblocksize=0):
        """Return sequence slices with line returns inserted at the end
        of each slice.
        
        Arguments:
            - seqstring: the raw sequence string
            - altblocksize: the length of sequence for writing to each
              line, default (0) means default value specified by blocksize
              will be used.
        """
        
        if altblocksize:
            self.block_size = altblocksize
        
        strlist = self.slicestringinblocks(seqstring, self.block_size)
        return '\n'.join(strlist) + "\n"
    

class PhylipWriter(_AlignmentWriter):
    def writealignment(self, alignmentdict, block_size=60, order=[]):
        """Write the alignment to a file.
        
        Arguments:
            - alignmentdict: dict of seqname -> seqstring.
            - blocksize: the sequence length to write to each line,
              default is 60
            - order: optional list of sequence names, which order to
              print in.
        (Assumes complete and correct list of names)
        """
        #setup
        if not order:
            order = alignmentdict.keys()
        self.setaligninfo(alignmentdict, order)
        self.setblocksize(block_size)
        
        # header
        self.file.write('%d  %d\n' %(self.number_sequences, self.align_length))
        # sequences (pretty much as writ by Gavin)
        
        for seqname in self.align_order:
            seq = alignmentdict[seqname]
            for block in range(0, self.align_length, self.block_size):
                if not block:
                    # write the otu name
                    if len(seqname) > 9:
                        warnings.warn('Name "%s" too long, truncated to "%s"'
                            % (seqname, seqname[:9]))
                        prefix = '%-10s' % seqname[:9]
                    else: prefix = '%-10s' % seqname
                else: prefix = ' ' * 10
                
                if block + self.block_size > self.align_length:
                        to = self.align_length
                else: to = block + self.block_size
                
                self.file.write('%s%s\n' % (prefix, seq[block:to]))
    

class PamlWriter(_AlignmentWriter):
    def writealignment(self, alignmentdict, block_size=60, order=[]):
        """Write the alignment to a file.
        
        Arguments:
            - alignmentdict: dict of seqname -> seqstring.
            - blocksize: the sequence length to write to each line,
              default is 60
            - order: optional list of sequence names, which order to
              print in.
        (Assumes order is a complete and correct list of names)
        """
        
        #setup
        if not order:
            order = alignmentdict.keys()
        self.setaligninfo(alignmentdict, order)
        self.setblocksize(block_size)
        
        #header
        self.file.write('%d  %d\n' % (self.number_sequences, self.align_length))
        
        #sequences
        for seq in self.align_order:
            self.file.writelines('%s\n%s' % (seq,self.wrapstringtoblocksize(alignmentdict[seq], altblocksize = block_size)))
    

class FastaWriter(_AlignmentWriter):
    def writealignment(self, alignmentdict, block_size=60, order=[]):
        """Write the alignment to a file.
        
        Arguments:
            - alignmentdict: dict of seqname -> seqstring.
            - blocksize: the sequence length to write to each line,
              default is 60
            - order: optional list of sequence names, which order to
              print in.
        
        (Assumes complete and correct list of names)
        """
        #setup
        if not order:
            order = alignmentdict.keys()
        self.setaligninfo(alignmentdict, order)
        self.setblocksize(block_size)
        
        #sequences
        for seq in self.align_order:
            self.file.writelines('>%s\n%s' % (seq,
                        self.wrapstringtoblocksize(alignmentdict[seq],
                                                   altblocksize = block_size)))
    

class GDEWriter(_AlignmentWriter):
    def writealignment(self, alignmentdict, block_size=60, order=[]):
        """Write the alignment to a file.
        
        Arguments:
            - alignmentdict: dict of seqname -> seqstring.
            - blocksize: the sequence length to write to each line,
              default is 60
            - order: optional list of sequence names, which order to
              print in.
        
        (Assumes complete and correct list of names)
        """
        
        #setup
        if not order:
            order = alignmentdict.keys()
        self.setaligninfo(alignmentdict, order)
        self.setblocksize(block_size)
        
        for seq in self.align_order:
            self.file.writelines('%s%s\n%s' % ("%", seq,
                self.wrapstringtoblocksize(alignmentdict[seq],
                                           altblocksize = block_size)))
    

# to add a new file format add it's suffix and class name here
WRITERS  = {
        'phylip': PhylipWriter,
        'paml': PamlWriter,
        'fasta': FastaWriter,
        'mfa': FastaWriter,
        'fa': FastaWriter,
        'gde': GDEWriter,
        }

