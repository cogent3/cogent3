import os

_NEW_TYPE = "COGENT3_NEW_TYPE" in os.environ


def nexus_from_alignment(aln, seq_type, wrap=50):
    """returns a nexus formatted string

    Parameters
    ----------
    seq_type
        dna, rna, or protein
    wrap
        the line width
    """
    if aln.is_ragged():
        raise ValueError(
            "Sequences in alignment are not all the same "
            + "length. Cannot generate NEXUS format.",
        )
    num_seq = len(aln.seqs)
    if not aln or not num_seq:
        return ""
    aln_len = aln.seq_len
    nexus_out = ["#NEXUS\n\nbegin data;"]
    nexus_out.append("    dimensions ntax=%d nchar=%d;" % (num_seq, aln_len))
    nexus_out.append(
        f"    format datatype={seq_type} interleave=yes missing=? " + "gap=-;",
    )
    nexus_out.append("    matrix")
    cur_ix = 0
    named_seqs = {a.name: a for a in aln.seqs} if _NEW_TYPE else aln.named_seqs
    names_seqs = sorted(named_seqs.items())
    while cur_ix < aln_len:
        nexus_out.extend(
            [f"    {x}    {y[cur_ix:cur_ix + wrap]}" for x, y in names_seqs],
        )
        nexus_out.append("")
        cur_ix += wrap
    nexus_out.append("    ;\nend;")

    return "\n".join(nexus_out)
