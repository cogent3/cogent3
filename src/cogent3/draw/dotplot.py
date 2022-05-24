from cogent3.align.align import dotplot
from cogent3.core.moltype import get_moltype
from cogent3.draw.drawable import Drawable
from cogent3.util.union_dict import UnionDict


__author__ = "Rahul Ghangas, Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell", "Rahul Ghangas"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
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

    if not hasattr(seq, "moltype"):
        seq = moltype.make_seq(seq)
    else:
        seq = moltype.make_seq(str(seq), seq.name)
    gap_map, seq = seq.parse_out_gaps()
    return gap_map, seq


def _convert_coords_for_scatter(coords):
    """returns x, y coordinates from [(start_x, start_y), (end_x, end_y), ..]

    Plotly scatter produces disjoint lines if the coordinates are separated
    by None"""
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
):
    """returns coordinates for forward / reverse strand"""
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


def get_align_coords(map1, map2, aligned=False):
    """sequence coordinates of aligned segments"""
    coords = None
    if not_gap(map1) and not_gap(map2):
        # no gaps
        if aligned:
            assert len(map1) == len(map2), "Aligned sequences inconsistent length"
            # but should return an alignment path
            coords = [[(0, 0), (len(map1), len(map2))]]
            coords = _convert_coords_for_scatter(coords)
        return coords

    assert len(map1) == len(map2), "aligned sequences not equal length"
    # diagonals are places where both sequences are NOT gaps
    # so we record start of a diagonal and when we hit a 'gap'
    # in either sequence, we have the end of the diagonal

    start_x = start_y = None
    coords = []
    for i in range(len(map1)):
        x_not_gap = not_gap(map1[i])
        y_not_gap = not_gap(map2[i])
        if x_not_gap and y_not_gap and start_x is None:
            start_x = len_seq(map1[:i])
            start_y = len_seq(map2[:i])
        elif (not x_not_gap or not y_not_gap) and start_x is not None:
            coord = [(start_x, start_y), (len_seq(map1[:i]) - 1, len_seq(map2[:i]) - 1)]
            start_x = start_y = None
            coords.append(coord)

    if start_x is not None:
        coord = [(start_x, start_y), (len_seq(map1) - 1, len_seq(map2) - 1)]
        coords.append(coord)

    return _convert_coords_for_scatter(coords)


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
        if hasattr(seq1, "moltype"):
            moltype = seq1.moltype
        else:
            moltype = get_moltype(moltype)

        map1, seq1 = _convert_input(seq1, moltype)
        map2, seq2 = _convert_input(seq2, moltype)
        len1, len2 = len(seq1), len(seq2)
        height = width * len2 / len1

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

        fwd, rev = get_dotplot_coords(
            self.seq1,
            self.seq2,
            window=window,
            threshold=threshold,
            min_gap=min_gap,
            rc=rc,
            show_progress=show_progress,
        )
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

        fwd, rev = self._fwd, self._rev
        if self.title is None:
            title = (
                f"Window={self._window}, Matched ≥ {self._threshold}/"
                f"{self._window} & Gap ≤ {self._min_gap}"
            )
        else:
            title = self.title

        self.layout |= dict(title=title)
        trace = UnionDict(
            type="scatter",
            x=fwd[0],
            y=fwd[1],
            name="+ strand",
            mode="lines",
            line=dict(color="blue"),
            xaxis=xaxis,
            yaxis=yaxis,
        )
        self.add_trace(trace)

        if rev:
            trace = UnionDict(
                type="scatter",
                x=rev[0],
                y=rev[1],
                name="- strand",
                mode="lines",
                line=dict(color="red"),
                xaxis=xaxis,
                yaxis=yaxis,
            )
            self.add_trace(trace)

        if self._aligned_coords:
            nudge = 0.2
            n, y = self._aligned_coords
            x = []
            for e in n:
                if e is not None:
                    e += nudge
                x.append(e)

            trace = UnionDict(
                type="scatter",
                x=x,
                y=y,
                name="Alignment",
                mode="lines",
                line=dict(color="black", dash="dot"),
                xaxis=xaxis,
                yaxis=yaxis,
            )
            self.add_trace(trace)
