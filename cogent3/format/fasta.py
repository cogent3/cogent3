#!/usr/bin/env python
"""Writer for FASTA sequence format
"""

from cogent3.format.formatter import _AlignmentFormatter

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight", "Gavin Huttley", "Thomas La"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"


class _fake_seq(str):
    """a holder for string sequences that allows provision of a seq.label
    attribute, required by fasta formatting funcs."""
    def __new__(cls, label, seq):
        new = str.__new__(cls, seq)
        new.label = label
        return new

    def __getitem__(self, *args, **kwargs):
        new_seq = str.__getitem__(self, *args, **kwargs)
        return self.__new__(self.__class__, self.label, new_seq)


def fasta_from_sequences(seqs, make_seqlabel=None, line_wrap=None):
    """Returns a FASTA string given a list of sequences. A sequence.label
       attribute takes precedence over sequence.name.

        - seqs can be a list of sequence objects or strings.
        - make_seqlabel: callback function that takes the seq object and returns
          a label str
        - line_wrap: a integer for maximum line width
    """
    fasta_list = []
    for i, seq in enumerate(seqs):
        # Check if it has a label, or one is to be created
        label = str(i)
        if make_seqlabel is not None:
            label = make_seqlabel(seq)
        elif hasattr(seq, 'label') and seq.label:
            label = seq.label
        elif hasattr(seq, 'name') and seq.name:
            label = seq.name

        # wrap sequence lines
        seq_str = str(seq)
        if line_wrap is not None:
            numlines, remainder = divmod(len(seq_str), line_wrap)
            if remainder:
                numlines += 1
            body = ["%s" % seq_str[j * line_wrap:(j + 1) * line_wrap]
                    for j in range(numlines)]
        else:
            body = ["%s" % seq_str]

        fasta_list.append('>' + label)
        fasta_list += body

    return '\n'.join(fasta_list)


def fasta_from_alignment(aln, make_seqlabel=None, line_wrap=None, sorted=True):
    """Returns a FASTA string given an alignment.

        - aln can be an Alignment object or dict.
        - make_seqlabel: callback function that takes the seq object and returns
          a label str
        - line_wrap: a integer for maximum line width
    """
    # get seq output order
    try:
        order = aln.names[:]
    except AttributeError:
        order = list(aln.keys())

    if sorted:
        order.sort()

    try:
        seq_dict = aln.named_seqs
    except AttributeError:
        seq_dict = aln

    ordered_seqs = []
    for label in order:
        seq = seq_dict[label]
        if isinstance(seq, str):
            seq = _fake_seq(label, seq)
        ordered_seqs.append(seq)
    return fasta_from_sequences(ordered_seqs, make_seqlabel=make_seqlabel,
                                line_wrap=line_wrap)


def alignment_to_fasta(alignmentdict, block_size=60, order=[]):
    """Returns a Fasta string given an alignment.
    """
    return FastaFormatter().format(alignmentdict, block_size, order)


class FastaFormatter(_AlignmentFormatter):

    def format(self, alignmentdict, block_size, order):
        """Format the alignment to Fasta.

        Arguments:
            - alignmentdict: dict of seqname -> seqstring.
            - blocksize: the sequence length to write to each line,
              default is 60
            - order: optional list of sequence names, which order to
              print in.

        (Assumes complete and correct list of names)
        """
        # setup
        if not order:
            order = list(alignmentdict.keys())
        self.setaligninfo(alignmentdict, order)
        self.setblocksize(block_size)

        return ''.join(['>%s\n%s' % (seq, self.wrapstringtoblocksize(alignmentdict[seq],
                       altblocksize=block_size)) for seq in self.align_order])
