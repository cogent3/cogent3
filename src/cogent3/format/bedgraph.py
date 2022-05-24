from cogent3.util.misc import get_merged_by_value_coords


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "alpha"

# following from https://cgwb.nci.nih.gov/goldenPath/help/bedgraph.html
# track type=bedGraph name=track_label description=center_label
#        visibility=display_mode color=r,g,b altColor=r,g,b
#        priority=priority autoScale=on|off alwaysZero=on|off
#        gridDefault=on|off maxHeightPixels=max:default:min
#        graphType=bar|points viewLimits=lower:upper
#        yLineMark=real-value yLineOnOff=on|off
#        windowingFunction=maximum|mean|minimum smoothingWindow=off|2-16
# Data Values
# Bedgraph track data values can be integer or real, positive or negative
# values. Chromosome positions are specified as 0-relative. The first
# chromosome position is 0. The last position in a chromosome of length N
# would be N - 1. Only positions specified have data. Positions not
# specified do not have data and will not be graphed. All positions specified
# in the input data must be in numerical order. The bedGraph format has four
# columns of data:

bedgraph_fields = (
    "name",
    "description",
    "visibility",
    "color",
    "altColor",
    "priority",
    "autoScale",
    "alwaysZero",
    "gridDefault",
    "maxHeightPixels",
    "graphType",
    "viewLimits",
    "yLineMark",
    "yLineOnOff",
    "windowingFunction",
    "smoothingWindow",
)

_booleans = ("autoScale", "alwaysZero", "gridDefault", "yLineOnOff")

valid_values = dict(
    autoScale=["on", "off"],
    graphType=["bar", "points"],
    windowingFunction=["maximum", "mean", "minimum"],
    smoothingWindow=["off"] + list(map(str, list(range(2, 17)))),
)


def raise_invalid_vals(key, val):
    """raises RuntimeError on invalid values for keys"""
    if key not in valid_values:
        return True
    if not str(val) in valid_values[key]:
        raise AssertionError(
            "Invalid bedgraph key/val pair: "
            + f"got {key}={val}; valid values are {valid_values[key]}"
        )


def booleans(key, val):
    """returns ucsc formatted boolean"""
    if val in (1, True, "on", "On", "ON"):
        val = "on"
    else:
        val = "off"
    return val


def get_header(name=None, description=None, color=None, **kwargs):
    """returns header line for bedgraph"""
    min_header = (
        'track type=bedGraph name="%(name)s" '
        + 'description="%(description)s" color=%(color)s'
    )

    assert None not in (name, description, color)
    header = [
        min_header
        % {"name": name, "description": description, "color": ",".join(map(str, color))}
    ]

    if kwargs:
        if not set(kwargs) <= set(bedgraph_fields):
            not_allowed = set(kwargs) - set(bedgraph_fields)
            raise RuntimeError(
                f"incorrect arguments provided to bedgraph {str(list(not_allowed))}"
            )

        if "altColor" in kwargs:
            kwargs["altColor"] = ",".join(map(str, kwargs["altColor"]))

        header_suffix = []
        for key in kwargs:
            if key in _booleans:
                kwargs[key] = booleans(key, kwargs[key])

            raise_invalid_vals(key, kwargs[key])
            header_suffix.append(f"{key}={kwargs[key]}")

        header += header_suffix

    return " ".join(header)


def bedgraph(
    chrom_start_end_val, digits=2, name=None, description=None, color=None, **kwargs
):
    """returns a bed formatted string. Input data must be provided as
    [(chrom, start, end, val), ...]. These will be merged such that adjacent
    records with the same value will be combined.

    Parameters
    ----------
    name
        track name
    description
        track description
    color
        (R,G,B) tuple of ints where max val of int is 255, e.g.
        red is (255, 0, 0)
    **kwargs
        keyword=val, .. valid bedgraph format modifiers
        see https://cgwb.nci.nih.gov/goldenPath/help/bedgraph.html

    """

    header = get_header(name=name, description=description, color=color, **kwargs)

    make_data_row = lambda x: "\t".join(list(map(str, x[:3])) + [f"{x[-1]:.2f}"])
    # get independent spans for each chromosome
    bedgraph_data = []
    data = []
    curr_chrom = None
    for chrom, start, end, val in chrom_start_end_val:
        if curr_chrom is None:
            curr_chrom = chrom

        if curr_chrom != chrom:
            data = get_merged_by_value_coords(data, digits=digits)
            bedgraph_data += [make_data_row([curr_chrom, s, e, v]) for s, e, v in data]
            data = []
            curr_chrom = chrom
        else:
            data.append([start, end, val])

    if data != []:
        data = get_merged_by_value_coords(data, digits=digits)
        bedgraph_data += [make_data_row([curr_chrom, s, e, v]) for s, e, v in data]

    bedgraph_data = [header] + bedgraph_data
    return "\n".join(bedgraph_data)
