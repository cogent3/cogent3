from collections import defaultdict
from typing import Optional, Union

from cogent3.core.alignment import SequenceCollection
from cogent3.core.genetic_code import GeneticCode, get_code
from cogent3.core.moltype import Alphabet, MolType, get_moltype

from .composable import NotCompleted, define_app
from .typing import SeqsCollectionType, SeqType, SerialisableType

GeneticCodeTypes = Union[str, int, GeneticCode]
MolTypes = Union[str, MolType]


def best_frame(
    seq: SeqsCollectionType,
    gc: GeneticCodeTypes = 1,
    allow_rc: bool = False,
    require_stop: bool = False,
):
    """returns reading frame start that has either no stops or a single
    terminal stop codon

    result will be either 1, 2, 3 (or -1, -2, -3)

    Parameters
    ----------
    gc
        genetic code ID, name or instance
    allow_rc
        If False, forward strand considered only. If True, and
        best frame on rc, it will be negative
    require_stop
        a terminal stop must be present

    Returns
    -------
    int
        1, 2, 3 if the best frame on the +_ strand; -1, -2, -3 if the best
        frame is on the reverse strand

    Raises
    ------
    ValueError
        if the minimum number of stop codons across all frames exceeds 1,
        or the the stop codon is not at the sequence end
    """
    gc = get_code(gc)
    translations = gc.sixframes(seq)
    if not allow_rc:
        translations = translations[:3]

    if not require_stop:
        # don't count stops if they're at the end of the aa sequence
        for i in range(len(translations)):
            if translations[i].endswith("*"):
                translations[i] = translations[i][:-1]

    stops_in_frame = [(tr.count("*"), i) for i, tr in enumerate(translations)]
    stops_in_frame.sort()
    min_stops, frame = stops_in_frame[0]
    # if min_stops > 1, cannot be translated
    if min_stops > 1:
        raise ValueError(f"{seq.name!r} cannot be robustly translated")
    elif min_stops == 0 and require_stop:
        # find seq with 1 stop
        min_stops = 20  # nonsense value
        for n, fr in stops_in_frame:
            if n == 1:
                min_stops, frame = n, fr
                break

    if not 0 <= min_stops <= 1:
        raise ValueError(f"{seq.name!r} cannot be robustly translated")

    if min_stops == 1 and not translations[frame].endswith("*"):
        raise ValueError(f"{seq.name!r} cannot be robustly translated")
    frame += 1
    if allow_rc and frame > 3:
        frame = 3 - frame
    return frame


def translate_frames(
    seq: SeqsCollectionType,
    moltype: Optional[MolTypes] = None,
    gc: GeneticCodeTypes = 1,
    allow_rc: bool = False,
):
    """translates a nucleic acid sequence

    Parameters
    ----------
    moltype
        molecular type, must be either DNA or RNA. Can be string or instance
    gc
        identifer for a genetic code or a genetic code instance
    allow_rc
        includes frames sequence reverse complement

    Returns
    -------
    [(frame, translation), ..]
    Reverse complement frame numbers are negative
    """
    gc = get_code(gc)
    if moltype:
        moltype = get_moltype(moltype)
        seq = moltype.make_seq(seq)

    translations = gc.sixframes(seq)
    if not allow_rc:
        translations = translations[:3]

    return translations


def get_fourfold_degenerate_sets(
    gc: GeneticCodeTypes, alphabet: Optional[Alphabet] = None, as_indices: bool = True
):
    """returns set() of codons that are 4-fold degenerate for genetic code gc

    Parameters
    ----------
    gc
        identifer for a genetic code or a genetic code instance
    alphabet
        nucleic acid Alphabet instance
    as_indices
        codons are represented as indices, rather than strings
    """
    four_fold = set()
    syns = gc.synonyms
    for codons in list(syns.values()):
        if len(codons) < 4:
            continue
        pos12s = defaultdict(list)
        for codon in codons:
            pos12s[codon[:2]].append(codon)

        for groups in list(pos12s.values()):
            if len(groups) == 4:
                four_fold.update([frozenset(groups)])

    if as_indices:
        assert alphabet is not None, "Must provide alphabet to convert to indices"
        ffold = set()
        to_indices = alphabet.to_indices
        for group in four_fold:
            grp = frozenset(tuple(to_indices(element)) for element in group)
            ffold.add(grp)
        four_fold = ffold

    return four_fold


@define_app
class select_translatable:
    """Selects translatable sequences by identifying the most likely reading
    frame. Sequences are truncated to modulo 3. seqs.info has a
    translation_errors entry."""

    def __init__(
        self,
        moltype: MolTypes = "dna",
        gc: GeneticCodeTypes = 1,
        allow_rc: bool = False,
        trim_terminal_stop: bool = True,
        frame: Optional[int] = None,
    ):
        """
        Parameters
        ----------
        moltype
            molecular type, must be either DNA or RNA. Can be string or
            instance
        gc
            identifier for a genetic code or a genetic code instance.
            see https://cogent3.org/doc/cookbook/what_codes.html
        allow_rc
            If False, forward strand considered only. If True, and
            best frame on rc, it will be negative
        trim_terminal_stop
            exclude terminal stop codon from seqs
        frame
            specify the coding frame (as an integer 1, 2, 3). If not specified
            (default) uses best_frame.

        Returns
        -------
        A sequence collection. Sequences that could not be translated
        are excluded.

        Examples
        --------
        Create a sample sequence collection and an app that returns the sequences
        which are translatable.

        >>> from cogent3 import make_unaligned_seqs, get_app
        >>> aln = make_unaligned_seqs(
        ...     {
        ...         "s1": "AATATAAATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTTCATAAAGTCATA",
        ...         "s1_rc": "TATGACTTTATGAAGTAATAAACTGCTGTTCTCATGCTGTAATGAGCTGGCATTTATATT",
        ...     }
        ... )
        >>> app = get_app("select_translatable")
        >>> result = app(aln)
        >>> result.to_dict()
        {'s1': 'AATATAAATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTTCATAAAGTCATA'}

        Use ``allow_rc=True`` to consider the reading frame for reverse strands.

        >>> app_rc = get_app("select_translatable", allow_rc=True)
        >>> result = app_rc(aln)
        >>> result.to_dict()
        {'s1': 'AATATAAATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTTCATAAAGTCATA', 's1_rc': 'AATATAAATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTTCATAAAGTCATA'}

        If you know the sequences are all in a specific frame, you can specify
        that using ``frame=<number>`` where the number is 1, 2 or 3. We make
        two sequences in frame 1 and add an internal stop to one

        >>> aln = make_unaligned_seqs(
        ...     {
        ...         "internal_stop": "AATTAAATGTGA",
        ...         "s2": "TATGACTAA",
        ...     }
        ... )
        >>> app = get_app("select_translatable", frame=1)
        >>> result = app(aln)
        >>> result.to_dict()
        {'s2': 'TATGAC'}
        """
        moltype = get_moltype(moltype)
        assert moltype.label.lower() in ("dna", "rna"), "Invalid moltype"

        if frame is not None:
            assert 1 <= frame <= 3, f"{frame} not 1, 2 or 3"
        self._frame = frame

        self._moltype = moltype
        self._gc = get_code(gc)
        self._allow_rc = allow_rc
        self._trim_terminal_stop = trim_terminal_stop

    T = Union[SerialisableType, SeqsCollectionType]

    def _get_frame(self, seq: SeqType) -> int:
        if self._frame is not None:
            tr = self._gc.translate(str(seq), start=self._frame - 1)
            if "*" in tr[:-1]:
                raise ValueError(f"internal stop in {seq.name}")
            return self._frame
        return best_frame(seq, self._gc, allow_rc=self._allow_rc)

    def main(self, seqs: SeqsCollectionType) -> T:
        """returns the translatable sequences from seqs.

        translation errors are stroed in the info object"""
        seqs = seqs.degap()
        if self._moltype and self._moltype != seqs.moltype:
            seqs = seqs.to_moltype(self._moltype)

        translatable = []
        error_log = []
        for seq in seqs.seqs:
            try:
                frame = self._get_frame(seq)
                if frame < 0:
                    seq = seq.rc()
                    frame *= -1
                frame -= 1  # returned from best frame as 1, 2, 3
                num_codons = (len(seq) - frame) // 3
                seq = seq[frame : frame + (num_codons * 3)]
                if self._trim_terminal_stop:
                    seq = seq.trim_stop_codon(gc=self._gc)
                translatable.append([seq.name, seq])
            except ValueError as msg:
                # TODO handle case where incomplete at end OR beginning
                # plus case where is divisible by 3 but not in frame
                # if not divisible by 3, then calc remainder as len(seq) % 3
                # try translating new[remainder:] and new[:-remainder]
                error_log.append([seq.name, msg.args[0]])

        if translatable:
            translatable = SequenceCollection(
                data=translatable, moltype=self._moltype, info=seqs.info
            )
            translatable.info["translation_errors"] = error_log
        else:
            translatable = NotCompleted("FALSE", self, " ".join(error_log), source=seqs)

        return translatable


@define_app
class translate_seqs:
    """Translates nucleic acid sequences into protein sequences, assumes in
    correct reading frame."""

    def __init__(
        self,
        moltype: MolTypes = "dna",
        gc: GeneticCodeTypes = 1,
        allow_rc: bool = False,
        trim_terminal_stop: bool = True,
    ):
        """
        Parameters
        ----------
        moltype
            molecular type, must be either DNA or RNA. Can be string or
            instance
        gc
            identifier for a genetic code or a genetic code instance.
            see https://cogent3.org/doc/cookbook/what_codes.html
        trim_terminal_stop
            exclude terminal stop codon from seqs

        Returns
        -------
        A sequence collection. Sequences that could not be translated
        are excluded.

        Examples
        --------

        Create a sample sequence collection and an app that translates nucleic
        acid sequences into protein sequences. By default, sequences that end
        with a stop codon are excluded with ``trim_terminal_stop=True``.

        >>> from cogent3 import make_aligned_seqs, get_app
        >>> aln = make_aligned_seqs({"s1": "ATGAGG", "s2": "ATGTAA"})
        >>> app_translate = get_app("translate_seqs", trim_terminal_stop=True)
        >>> result = app_translate(aln)
        >>> print(result.to_pretty())
        s1    MR
        s2    .-
        """
        moltype = get_moltype(moltype)
        assert moltype.label.lower() in ("dna", "rna"), "Invalid moltype"

        self._moltype = moltype
        self._gc = get_code(gc)
        self._trim_terminal_stop = trim_terminal_stop

    T = Union[SerialisableType, SeqsCollectionType]

    def main(self, seqs: SeqsCollectionType) -> T:
        """returns translated sequences"""
        if self._moltype and self._moltype != seqs.moltype:
            seqs = seqs.to_moltype(self._moltype)

        return seqs.get_translation(gc=self._gc, trim_stop=self._trim_terminal_stop)
