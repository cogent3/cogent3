from cogent3.core.moltype import get_moltype
from cogent3.draw.drawable import Drawable, _iplot
from cogent3.align.align import dotplot

__author__ = "Rahul Ghangas, Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell", "Rahul Ghangas"]
__license__ = "GPL"
__version__ = "3.0a2"
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
            p *= (i + mismatches)
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
    """converts to Annotatable map instyance and a Sequence object"""
    if hasattr(seq, 'map'):
        # already an Aligned instance
        gap_map = seq.map
        seq = seq.data
        return gap_map, seq

    if not hasattr(seq, 'moltype'):
        seq = moltype.make_seq(seq)
    else:
        seq = moltype.make_seq(str(seq), seq.name)
    gap_map, seq = seq.parse_out_gaps()
    return gap_map, seq


def _convert_coords_for_scatter(coords):
    """returns x, y coordinates from [(start_x, start_y), (end_x, end_y), ..]

    Plotly scatter produces disjoint lines if the coordinates are separated
    by None"""
    new = {'x': [], 'y': []}
    for (x1, y1), (x2, y2) in coords:
        new['x'].extend([x1, x2, None])
        new['y'].extend([y1, y2, None])

    # drop trailing None
    x = new['x'][:-1]
    y = new['y'][:-1]
    return x, y


def get_align_coords(map1, map2):
    """sequence coordinates of aligned segments"""
    if not_gap(map1) and not_gap(map2):
        # no gaps
        return None

    assert len(map1) == len(map2), 'aligned sequences not equal length'
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
            coord = [(start_x, start_y),
                     (len_seq(map1[:i]) - 1,
                      len_seq(map2[:i]) - 1)]
            start_x = start_y = None
            coords.append(coord)

    if start_x is not None:
        coord = [(start_x, start_y),
                 (len_seq(map1) - 1,
                  len_seq(map2) - 1)]
        coords.append(coord)

    return _convert_coords_for_scatter(coords)


class Display2D(Drawable):

    def __init__(self, seq1, seq2, moltype='text'):
        if hasattr(seq1, 'moltype'):
            moltype = seq1.moltype
        else:
            moltype = get_moltype(moltype)

        map1, seq1 = _convert_input(seq1, moltype)
        map2, seq2 = _convert_input(seq2, moltype)

        self.seq1 = seq1
        self.seq2 = seq2
        self._aligned_coords = get_align_coords(map1, map2)
        self._cache = {}
        self._trace = None

    def _calc_lines(self, window, threshold, min_gap):
        # Cache dotplot line segment coordinates as they can sometimes
        # be re-used at different resolutions, colours etc.
        (len1, len2) = len(self.seq1), len(self.seq2)
        if threshold is None:
            universe = (len1 - window) * (len2 - window)
            acceptable_noise = min(len1, len2) / window
            threshold = suitable_threshold(
                window, acceptable_noise / universe)

        key = (min_gap, window, threshold)
        if key not in self._cache:
            fwd = dotplot(str(self.seq1), str(self.seq2),
                          window, threshold, min_gap, None)
            if hasattr(self.seq1, "reverse_complement"):
                rev = dotplot(str(self.seq1.reverse_complement()),
                              str(self.seq2), window, threshold, min_gap, None)
                rev = [((len1 - x1, y1), (len1 - x2, y2))
                       for ((x1, y1), (x2, y2)) in rev]
                rev = _convert_coords_for_scatter(rev)
            else:
                rev = []

            # convert for Plotly scatter
            fwd = _convert_coords_for_scatter(fwd)
            self._cache[key] = (fwd, rev)

        return self._cache[key]

    def get_trace(self, window=20, threshold=None, min_gap=0, width=500,
                  title=None):
        # calculate the width based on ratio of seq lengths
        if self._trace is None:
            self._set_initial_layout(window=window, threshold=threshold,
                                     min_gap=min_gap, width=width,
                                     title=title)

        return list(self._trace['data'])

    def set_layout(self, layout_updates):
        self._trace['layout'] = layout_updates

    def update_layout(self, layout_updates):
        self._trace['layout'].update(layout_updates)

    def _set_initial_layout(self, width=500, title=None, window=20,
                            min_gap=0, threshold=None, **kw):
        import plotly.graph_objs as go

        height = width * len(self.seq2) / len(self.seq1)
        super(Display2D, self)._set_initial_layout(
            width, height, visible_axes=True, showlegend=True)

        if title is not None:
            self._trace['layout']['title'] = title

        if self.seq1.name:
            self._trace['layout']['xaxis']['title'] = self.seq1.name

        if self.seq2.name:
            self._trace['layout']['yaxis']['title'] = self.seq2.name

        fwd, rev = self._calc_lines(window, threshold, min_gap)
        trace = go.Scatter(x=fwd[0], y=fwd[1], name='+ strand',
                           mode='lines',
                           line=dict(color='blue'))
        self._trace['data'] = [trace]

        if rev:
            trace = go.Scatter(x=rev[0], y=rev[1], name='- strand',
                               mode='lines',
                               line=dict(color='red'))
            self._trace['data'].append(trace)

        if self._aligned_coords:
            nudge = 0.2
            n, y = self._aligned_coords
            x = []
            for e in n:
                if e:
                    e += nudge
                x.append(e)

            trace = go.Scatter(x=x, y=y, name='Alignment',
                               mode='lines',
                               line=dict(color='black', dash='dot'))
            self._trace['data'].append(trace)
