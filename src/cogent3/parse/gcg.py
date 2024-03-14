import contextlib
import warnings


def MsfParser(f):
    """Read sequences from a msf format file"""
    # parse optional header
    # parse optional text information
    # file header and sequence header are seperated by a line ending in '..'
    with contextlib.suppress(AttributeError):
        f = f.read().splitlines()

    f.pop(0)
    for line in f:
        line = line.strip()
        if line.endswith(".."):
            break
    # parse sequence info
    seqinfo = {}
    for line in f:
        line = line.strip()
        if line.startswith("//"):
            break
        line = line.split()
        if line and line[0] == "Name:":
            seqinfo[line[1]] = int(line[3])
    # parse sequences
    sequences = {}
    for line in f:
        if line := line.strip().split():
            if line[0] in sequences:
                sequences[line[0]] += "".join(line[1:])
            elif line[0] in seqinfo:
                sequences[line[0]] = "".join(line[1:])
    # consistency check
    if len(sequences) != len(seqinfo):
        warnings.warn(
            "Number of loaded seqs[%s] not same as "
            "expected[%s]." % (len(sequences), len(seqinfo))
        )
    for name, value_ in sequences.items():
        if len(value_) != seqinfo[name]:
            warnings.warn(
                "Length of loaded seqs [%s] is [%s] not "
                "[%s] as expected." % (name, len(sequences[name]), seqinfo[name])
            )

    # yield sequences
    yield from sequences.items()
