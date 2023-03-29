from cogent3.align.pycompare import (
    MatchedSeqPaths,
    SeqKmers,
    _calc_seed_size,
    find_matched_paths,
)
from cogent3.core.moltype import get_moltype
from cogent3.draw.drawable import Drawable
from cogent3.util.union_dict import UnionDict


__author__ = "Rahul Ghangas, Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell", "Rahul Ghangas"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


def suitable_threshold(window, desired_probability):
    """Use cumulative binomial distribution to find the number of identical
    bases which we expect a nucleotide window-mer to have with the desired
    probability"""
    cumulative_p = 0.0
    for matches in range(window, 0, -1):
        mismatches = window - matches
        p = 0.75 ** mismatches
        for i in range(matches, 0, -1):  # n
            p *= i + mismatches
            p /= i
            p *= 0.25
        cumulative_p += p
        if cumulative_p > desired_probability:
            break
    return matches


def len_seq(span):
    """length of a Annotatable map object"""
    return len(span.nongap())


def not_gap(span):
    """whether a span corresponds to a non-gap"""
    return len(span.gaps()) == 0


def _convert_input(seq, moltype):
    """converts to Annotatable map instance and a Sequence object"""
    if hasattr(seq, "map"):
        # already an Aligned instance
        gap_map = seq.map
        seq = seq.data
        return gap_map, seq

    seq = (
        moltype.make_seq(str(seq), seq.name)
        if hasattr(seq, "moltype")
        else moltype.make_seq(seq)
    )

    gap_map, seq = seq.parse_out_gaps()
    return gap_map, seq


def _convert_coords_for_scatter(coords):  # pragma: no cover
    """discontinued, use cogent3.align.compare.MatchedSeqPaths instead"""
    from cogent3.util.warning import discontinued

    discontinued(
        "function",
        "_convert_coords_for_scatter",
        "2023.5",
        "replaced by align.compare.MatchedSeqPaths",
    )
    new = {"x": [], "y": []}
    for (x1, y1), (x2, y2) in coords:
        new["x"].extend([x1, x2, None])
        new["y"].extend([y1, y2, None])

    # drop trailing None
    x = new["x"][:-1]
    y = new["y"][:-1]
    return x, y


def get_dotplot_coords(
    seq1, seq2, window=20, threshold=None, min_gap=0, rc=None, show_progress=False
):  # pragma: no cover
    """discontinued, use cogent3.align.compare.find_matched_paths"""
    from cogent3.align.align import dotplot
    from cogent3.util.warning import discontinued

    discontinued(
        "function",
        "get_dotplot_coords",
        "2023.5",
        "replaced by much faster code in align.compare",
    )

    (len1, len2) = len(seq1), len(seq2)
    if threshold is None:
        universe = (len1 - window) * (len2 - window)
        acceptable_noise = min(len1, len2) / window
        threshold = suitable_threshold(window, acceptable_noise / universe)

    fwd = dotplot(
        str(seq1),
        str(seq2),
        window,
        threshold,
        min_gap,
        None,
        show_progress=show_progress,
    )
    if hasattr(seq1, "reverse_complement") and rc:
        rev = dotplot(
            str(seq1.reverse_complement()),
            str(seq2),
            window,
            threshold,
            min_gap,
            None,
            show_progress=show_progress,
        )
        rev = [((len1 - x1, y1), (len1 - x2, y2)) for ((x1, y1), (x2, y2)) in rev]
        rev = _convert_coords_for_scatter(rev)
    else:
        rev = []

    # convert for Plotly scatter
    fwd = _convert_coords_for_scatter(fwd)
    return fwd, rev


def get_align_coords(map1, map2, aligned=False) -> MatchedSeqPaths:
    """sequence coordinates of aligned segments"""
    from cogent3.align.pycompare import segment

    paths = MatchedSeqPaths("Alignment")
    if not_gap(map1) and not_gap(map2):
        # no gaps
        if aligned:
            assert len(map1) == len(map2), "Aligned sequences inconsistent length"
            paths[0].append((segment(0, len(map1)), segment(0, len(map2))))
        return paths

    assert len(map1) == len(map2), "aligned sequences not equal length"
    # diagonals are places where both sequences are NOT gaps
    # so we record start of a diagonal and when we hit a 'gap'
    # in either sequence, we have the end of the diagonal

    start_x = start_y = None
    for i in range(len(map1)):
        x_not_gap = not_gap(map1[i])
        y_not_gap = not_gap(map2[i])
        if x_not_gap and y_not_gap and start_x is None:
            start_x = len_seq(map1[:i])
            start_y = len_seq(map2[:i])
        elif (not x_not_gap or not y_not_gap) and start_x is not None:
            paths[start_y - start_x].append(
                (
                    segment(start_x, len_seq(map1[:i]) - 1),
                    segment(start_y, len_seq(map2[:i]) - 1),
                )
            )
            start_x = start_y = None

    if start_x is not None:
        paths[start_y - start_x].append(
            (segment(start_x, len_seq(map1) - 1), segment(start_y, len_seq(map2) - 1))
        )

    return paths


class Dotplot(Drawable):
    """calculates matches between sequences and displays as a dotplot"""

    def __init__(
        self,
        seq1,
        seq2,
        is_aligned,
        moltype="text",
        window=20,
        threshold=None,
        k=None,
        min_gap=0,
        rc=False,
        xtitle=None,
        ytitle=None,
        title=None,
        width=500,
        show_progress=False,
    ):
        """
        Parameters
        ----------
        seq1, seq2 : string or sequence object
        moltype : str or MolType instance
            if seq1, seq2 are strings, moltype is used to convert to sequence
            objects
        window : int
            k-mer size for comparison between sequences
        threshold : int
            windows where the sequences are identical >= threshold are a match
        k : int
            size of k-mer to break sequences into. Larger values increase
            speed but reduce resolution. If not specified, is computed as the
            maximum of (window-threshold), (window % k) * k <= threshold.
        min_gap : int
            permitted gap for joining adjacent line segments, default is no gap
            joining
        rc : bool or None
            include dotplot of reverse compliment also. Only applies to Nucleic
            acids moltypes
        xtitle, ytitle
            name of the seq1, seq2. None if included as part of a
            AnnotatedDrawable
        title : str
            title for the plot
        show_progress : bool
            displays progress bar
        """

        # we ensure sequences have gaps parsed and the calculate aspect ratio
        moltype = seq1.moltype if hasattr(seq1, "moltype") else get_moltype(moltype)
        map1, seq1 = _convert_input(seq1, moltype)
        map2, seq2 = _convert_input(seq2, moltype)
        len1, len2 = len(seq1), len(seq2)
        height = width * len2 / len1

        if seq1.name == seq2.name is None:
            seq1.name = "seq1"
            seq2.name = "seq2"

        super(Dotplot, self).__init__(
            visible_axes=True,
            showlegend=True,
            width=width,
            height=height,
            xtitle=xtitle,
            ytitle=ytitle,
        )

        self.seq1 = seq1
        self.seq2 = seq2
        self._aligned_coords = get_align_coords(map1, map2, aligned=is_aligned)

        self.title = title
        self._window = window
        self._min_gap = min_gap
        if threshold is None:
            universe = (len1 - window) * (len2 - window)
            acceptable_noise = min(len1, len2) / window
            threshold = suitable_threshold(window, acceptable_noise / universe)

        self._threshold = threshold

        if k is None:
            k = _calc_seed_size(window, threshold)

        sk = SeqKmers(seq1, k=k, canonical=set(seq1.moltype))
        seq2 = None if seq1 == seq2 else seq2
        fwd = find_matched_paths(
            sk, seq1=seq1, seq2=seq2, window=window, threshold=threshold
        )
        fwd.name = "+ strand"
        if rc:
            if sk.num_seqs == 2:
                sk.drop_seq(seq2.name)

            seq2 = seq1.rc() if seq2 is None else seq2.rc()
            seq2.name = f"{seq2.name}-rc"
            rev = find_matched_paths(
                sk, seq1=seq1, seq2=seq2, window=window, threshold=threshold
            )
            rev.name = "- strand"
        else:
            rev = None

        self._fwd = fwd
        self._rev = rev

    def _build_fig(self, xaxis="x", yaxis="y"):
        # calculate the width based on ratio of seq lengths
        layout = UnionDict()
        if self.xtitle:
            layout |= dict(
                xaxis=dict(
                    title=self.xtitle, mirror=True, showgrid=False, showline=True
                )
            )

        if self.ytitle:
            layout |= dict(
                yaxis=dict(
                    title=self.ytitle, mirror=True, showgrid=False, showline=True
                )
            )

        self.layout |= dict(layout)
        self.layout |= dict(
            yaxis=dict(range=[0, len(self.seq2)]), xaxis=dict(range=[0, len(self.seq1)])
        )

        if self.title is None:
            title = (
                f"Window={self._window}, Matched ≥ {self._threshold}/"
                f"{self._window} & Gap ≤ {self._min_gap}"
            )
        else:
            title = self.title

        self.layout |= dict(title=title)
        trace = UnionDict(
            line=dict(color="blue"),
            xaxis=xaxis,
            yaxis=yaxis,
        )
        trace |= self._fwd.plotly_trace()
        self.add_trace(trace)

        if self._rev:
            trace = UnionDict(
                line=dict(color="red"),
                xaxis=xaxis,
                yaxis=yaxis,
            )
            trace |= self._rev.plotly_trace(rc=True, length=len(self.seq2))
            self.add_trace(trace)

        if self._aligned_coords:
            trace = UnionDict(
                line=dict(color="black", dash="dot"),
                xaxis=xaxis,
                yaxis=yaxis,
            )
            trace |= self._aligned_coords.plotly_trace()
            nudge = 0.2
            for i, v in enumerate(trace.x):
                if v is None:
                    continue
                trace.x[i] = v + nudge
            self.add_trace(trace)
