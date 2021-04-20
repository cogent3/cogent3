#!/usr/bin/env python

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2021, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2021.04.20a"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def PamlParser(f):
    d = f.readline().split()
    numseqs, seqlen = int(d[0]), int(d[1])
    for i in range(numseqs):
        seqname = f.readline().strip()
        if not seqname:
            raise ValueError("Sequence name missing")
        currseq = []
        length = 0
        while length < seqlen:
            seq_line = f.readline()
            if not seq_line:
                raise ValueError(
                    'Sequence "%s" is short: %s < %s' % (seqname, length, seqlen)
                )
            seq_line = seq_line.strip()
            length += len(seq_line)
            currseq.append(seq_line)

        yield (seqname, "".join(currseq))
