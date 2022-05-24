from cogent3.draw.drawable import Drawable
from cogent3.draw.letter import letter_stack
from cogent3.util.dict_array import DictArray
from cogent3.util.union_dict import UnionDict


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


def get_mi_char_heights(fa):
    """computes character heights from MI terms
    Parameters
    ----------
    counts : MotifCountsArray

    Returns
    -------
    DictArray
    """
    I = fa.information()
    mit = I * fa.array.T
    return fa.template.wrap(mit.T)


def get_base_logo_layout(axnum, xtick_fontsize, ytick_fontsize):
    """creates default plotly layout for drawing a sequence logo
    Parameters
    ----------
    axnum : int
        axis number for plotly
    xtick_fontsize, ytick_fontsize : int
        font size for tick values

    Returns
    -------
    UnionDict
    """
    layout = UnionDict()

    # set x and y range, set y ticks
    axis_lines = dict(
        mirror=True,
        linewidth=1,
        showgrid=False,
        linecolor="black",
        showline=True,
        visible=True,
        zeroline=False,
    )
    axis = "axis" if axnum == 1 else f"axis{axnum}"
    xanchor = "x" if axnum == 1 else f"x{axnum}"
    yanchor = "y" if axnum == 1 else f"y{axnum}"
    layout[f"x{axis}"] = dict(
        anchor=yanchor,
        tickfont=dict(size=xtick_fontsize),
        ticks="inside",
    )

    layout[f"y{axis}"] = dict(
        tickfont=dict(size=ytick_fontsize),
        title="Bits",
        anchor=xanchor,
        ticks="inside",
    )
    layout[f"x{axis}"] |= axis_lines
    layout[f"y{axis}"] |= axis_lines
    layout.template = "plotly_white"
    return layout


def _char_hts_as_lists(data):
    """returns a [[(base, height), ..], ..]"""
    # data is assumed row-oriented
    result = []
    for d in data:
        try:
            d = d.to_dict()
        except AttributeError:
            # assume it's just a dict
            pass

        if d:
            d = list(d.items())
        else:
            d = None
        result.append(d)

    return result


_dna_colours = dict(A="green", T="red", C="blue", G="orange")


def get_logo(
    char_heights,
    axnum=1,
    height=400,
    width=800,
    ylim=None,
    ydomain=None,
    colours=None,
    layout=None,
):
    """
    Parameters
    ----------
    char_heights
        a seris of [[(base, value), ..], ..] or series of dicts with [{letter1: value, letter2: value}, ...]
        Dicts are coerced to list of lists. If values are < 0, the letter is inverted. Empty elements are ignored.
    axnum : int
        plotly axis number
    height, width: int
        figure dimensions in pixels
    ylim : float
        maximum y-value
    ydomain
        [start, end], specifies vertical positioning for this logo
    colours : dict
        dict mapping characters to colours. Defaults to custom 'dna' colours
        typically used for DNA
    layout : UnionDict
        Customised base layout

    Returns
    -------
    Drawable
    """
    if isinstance(char_heights, DictArray) or isinstance(char_heights[0], dict):
        char_heights = _char_hts_as_lists(char_heights)

    colours = colours or _dna_colours
    if layout is None:
        layout = get_base_logo_layout(axnum, 12, 12)

    stack_data = []
    est_ylim = 0
    for d in char_heights:
        if not d:
            stack_data.append(None)
            continue

        d = sorted(d, key=lambda x: x[1])

        if ylim is None:
            est_ylim = max(est_ylim, max([e[-1] for e in d]))
        stack_data.append(d)

    stacks = []
    for index, stack in enumerate(stack_data):
        if stack is None:
            continue
        middle, stack_shapes = letter_stack(stack, index - 0.5, 1, colours, axnum)
        stacks += stack_shapes

    layout["shapes"] = stacks

    if ylim is None:
        ylim = est_ylim * 1.05

    yaxis = "yaxis" if axnum == 1 else f"yaxis{axnum}"
    layout[yaxis]["range"] = [0, ylim]
    if ydomain:
        layout[yaxis]["domain"] = ydomain

    return Drawable(layout=layout, height=height, width=width)
