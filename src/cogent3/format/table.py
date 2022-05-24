#!/usr/bin/env python
"""
Tool for creating tables and representing them as text, or writing to file for
import into other packages. These classes still under development.

Current formats include restructured text (keyed by 'rest'), latex, html,
columns separated by a provided string, and a simple text format.
"""
import re
import textwrap

from xml.sax.saxutils import escape

import numpy


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell", "Matthew Wakefield", "Jeremy Widmann"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

known_formats = (
    "bedgraph",
    "phylip",
    "rest",
    "rst",
    "markdown",
    "md",
    "latex",
    "tex",
    "html",
    "simple",
    "csv",
    "tsv",
)

css_c3table_template = "\n".join(
    (
        ".c3table table {margin: 10px 0;}",
        ".c3table tr:last-child {border-bottom: 1px solid #000;} ",
        ".c3table tr > th {text-align: left; padding: 0 5px;}",
        ".c3table tr > td {text-align: left; padding: 5px;}",
        ".c3table tr:nth-child(even) {background: #f7f7f7 !important;}",
        ".c3table .ellipsis {background: rgba(0, 0, 0, .01);}",
        ".c3table .index {background: %(colour)s; margin: 10px; font-weight: 600;}",
        ".c3table .head_cell {background: %(head_colour)s; font-weight: bold; text-align: center;}",
        ".c3table caption {color: rgb(250, 250, 250); background: "
        "rgba(30, 140, 200, 1); padding: 3px; white-space: nowrap; "
        "caption-side: top;}",
        ".c3table .cell_title {font-weight: bold;}",
        ".c3col_left { text-align: left !important; display: block;}",
        ".c3col_right { text-align: right !important; display: block;}",
        ".c3col_center { text-align: center !important; display: block;}",
    )
)


def _merged_cell_text_wrap(text, max_line_length, space):
    """left justify wraps text into multiple rows"""
    max_line_width = max_line_length - (2 * space)
    if len(text) < max_line_length:
        return [text]
    buffer = " " * space
    wrapped = textwrap.wrap(
        text, width=max_line_width, initial_indent=buffer, subsequent_indent=buffer
    )
    wrapped = [f"{line.ljust(max_line_width + 2 * space)}" for line in wrapped]
    return wrapped


def _merge_cells(row):
    """merges runs of identical row cells.

    returns a list with structure [((span_start, span_end), cell value),..]"""
    new_row = []
    last = 0
    span = 1  # the minimum
    for i in range(1, len(row), 1):
        if row[i - 1] != row[i]:
            new_row.append(((last, last + span), row[i - 1]))
            last = i
            span = 1
            continue
        span += 1

    new_row.append(((last, last + span), row[i]))
    return new_row


def rich_html(
    rows,
    row_cell_func=None,
    header=None,
    header_cell_func=None,
    element_formatters=None,
    merge_identical=True,
    compact=True,
    caption=None,
):
    """returns just the html Table string

    Parameters
    ----------
    rows
        table rows
    row_cell_func
        callback function that formats the row values. Must
        take the row value and coordinates (row index, column index).
    header
        the table header
    header_cell_func
        callback function that formats the column headings
        must take the header label value and coordinate
    element_formatters
        a dictionary of specific callback funcs for
        formatting individual html table elements.
        e.g. {'table': lambda x: '<table border="1" class="docutils">'}
    merge_identical
        cells within a row are merged to one span.
    caption
        Table title / legend

    Note: header_cell_func and row_cell_func override element_formatters.
    """
    element_formatters = element_formatters or {}
    formatted = element_formatters.get
    data = [formatted("table", "<table>")]
    if caption:
        data.append(
            '<caption style="font-weight: bold;"'
            'background:rgba(30, 140, 200, 1)"; '
            f'align="top">{caption}</caption>'
        )

    if row_cell_func is None:

        def row_cell_func(v, r, c):
            return f"<td>{v}</td>"

    if header_cell_func is None:

        def header_cell_func(v, c):
            return f"<th>{v}</th>"

    if merge_identical:
        row_iterator = _merge_cells
    else:
        row_iterator = enumerate

    if header:
        thead = formatted("thead", '<thead style="font-weight: bold;">')
        row = [header_cell_func(escape(label), i) for i, label in enumerate(header)]
        data += [thead] + row + ["</thead>"]

    formatted_rows = []
    for ridx, row in enumerate(rows):
        new = [formatted("tr", "<tr>")]
        for cidx, cell in row_iterator(row):
            new += [row_cell_func(escape(cell), ridx, cidx)]
        new += ["</tr>"]
        formatted_rows += new

    tbody = formatted("tbody", "<tbody>")
    data += [tbody] + formatted_rows + ["</tbody>"]
    data += ["</table>"]
    if compact:
        data = "".join(data)
    else:
        data = "\n".join(data)
    return data


def latex(
    rows,
    header=None,
    caption=None,
    legend=None,
    justify=None,
    label=None,
    position=None,
):
    """Returns the text a LaTeX table.

    Parameters
    ----------
    rows
        table data in row orientation
    header
        table header
    caption
        title text.
    legend
        If provided, the text is placed in a \\caption*{} command at the
        bottom of the table and the caption is placed at the top.
    justify
        column justification, default is right aligned.
    label
        for cross referencing
    position
        table page position, default is here, top separate page

    Notes
    -----
    The \\caption*{} command is provided with the caption package. See
    https://ctan.org/pkg/caption for more details.
    """

    if not justify:
        numcols = [len(header), len(rows[0])][not header]
        justify = "r" * numcols

    justify = "{ %s }" % " ".join(list(justify))
    if header:
        header = "%s \\\\" % " & ".join([r"\bf{%s}" % head.strip() for head in header])
    rows = [f"{' & '.join(row)} \\\\" for row in rows]
    position = position or "htp!"
    table_format = [
        r"\begin{table}[%s]" % position,
        r"\centering",
        r"\begin{tabular}%s" % justify,
        r"\hline",
        header,
        r"\hline",
        r"\hline",
    ]
    table_format += rows
    table_format.append(r"\hline")
    table_format.append(r"\end{tabular}")

    caption = r"\caption{%s}" % caption if caption else ""
    label = r"\label{%s}" % label if label else ""
    legend = r"\caption*{%s}" % legend if isinstance(legend, str) else None
    if caption and label:
        caption = f"{caption}\n{label}"
    elif caption or label:
        caption = caption or label

    if caption and legend:
        table_format.insert(2, caption)
    elif caption:
        table_format.append(caption)

    if legend is not None:
        table_format.append(legend)
    table_format.append(r"\end{table}")

    return "\n".join(table_format)


def get_continuation_tables(
    header, formatted_table, identifiers=None, space=2, max_width=1e100
):
    """returns series of tables segmented to not exceed max_width"""
    tables = []
    try:
        space = " " * space
    except TypeError:
        pass

    # if we are to split the table, creating sub tables, determine
    # the boundaries
    if len(space.join(header)) < max_width:
        return [(header, formatted_table)]

    # having determined the maximum string lengths we now need to
    # produce subtables of width <= max_width
    col_widths = [len(head) for head in header]
    sep = len(space)
    min_length = col_widths[0]

    if min_length > max_width:
        raise RuntimeError("Maximum width too small for identifiers")

    # if we have an index column, every new table block includes that width
    # in calculating the number of columns; otherwise it's simply the sum
    if identifiers:
        id_width = col_widths[0] + sep
        begin = 1
    else:
        id_width = 0
        begin = 0

    width = id_width
    boundaries = []
    for i in range(begin, len(header)):
        width += col_widths[i] + sep
        if width > max_width:
            boundaries.append((begin, i))
            begin = i
            width = id_width + col_widths[i]

    boundaries.append((begin, len(header)))
    data = {c[0].strip(): c[1:] for c in zip(header, *formatted_table)}
    for start, end in boundaries:
        if identifiers:
            subhead = header[:1] + header[start:end]
        else:
            subhead = header[start:end]
        rows = numpy.array([data[c.strip()] for c in subhead], dtype="<U15")
        if rows.ndim == 1:
            rows = [rows.tolist()]
        else:
            rows = rows.T.tolist()
        tables.append((subhead, rows))

    return tables


def simple_format(
    header,
    formatted_table,
    title=None,
    legend=None,
    max_width=1e100,
    identifiers=None,
    borders=True,
    space=2,
):
    """Returns a table in a simple text format.

    Parameters
    ----------
    header
        series with column headings
    formatted_table
        a two dimensional structure (list/tuple) of strings
        previously formatted to the same width within a column.
    title
        optional table title
    legend
        optional table legend
    max_width
        forces wrapping of table onto successive lines if its'
        width exceeds that specified
    identifiers
        index for the column that uniquely identify rows. Required if table
        width exceeds max_width.
    borders
        whether to display borders.
    space
        minimum number of spaces between columns.

    """
    table = []
    try:
        space = " " * space
    except TypeError:
        pass

    # if we are to split the table, creating sub tables, determine
    # the boundaries
    subtables = get_continuation_tables(
        header, formatted_table, identifiers, space, max_width
    )
    for i, (h, t) in enumerate(subtables):
        st = title if i == 0 else f"continued: {title}"
        if st:
            table.append(st)
        sh = space.join(h)
        length_head = len(sh)
        if borders:
            table.extend(["=" * length_head, sh, "-" * length_head])
        else:
            table.append(sh)
        rows = [space.join(r) for r in t]
        rows = "\n".join(rows)
        if rows:
            table.append(rows)
        if borders:
            table.append("-" * length_head)
        if len(subtables) > 1:
            table.append("")

    # add the legend, wrapped to the table widths
    if legend:
        wrapped = _merged_cell_text_wrap(legend, max_width, 0)
        table += wrapped

    return "\n".join(table)


_pipe = re.compile(r"\|")


def _escape_pipes(formatted_table, header):
    """returns text with | replaced by \\|, adjusting column widths"""
    resized = False
    widths = list(map(len, formatted_table[0]))
    num_rows = len(formatted_table)
    num_cols = len(formatted_table[0])
    for i in range(num_rows):
        for j in range(num_cols):
            cell = formatted_table[i][j]
            if "|" in cell:
                cell = _pipe.sub(r"\|", cell)
                formatted_table[i][j] = cell
                widths[j] = max(len(cell), widths[j])
                resized = True

    if resized:
        for j in range(num_cols):
            header[j] = header[j].center(widths[j])
            for i in range(num_rows):
                cell = formatted_table[i][j]
                formatted_table[i][j] = cell.center(widths[j])

    return formatted_table, header


def markdown(header, formatted_table, space=1, justify=None):
    """Returns a table in Markdown format

    Parameters
    ----------
    header
        series with column headings
    formatted_table
        a two dimensional structure (list/tuple) of strings
        previously formatted to the same width within a column.
    space
        number of spaces surrounding the cell contents, must be >= 1
    justify
        characters indicating alignment of columns
    """
    assert space >= 1, "space must be >= 1"
    if justify is not None:
        assert len(justify) == len(
            header
        ), "column number and justify entries must match"
        justify = [c.lower() for c in justify]

    formatted_table, header = _escape_pipes(formatted_table, header)

    row_template = "| %s |"
    sep = "".join([" " * space, "|", " " * space])
    divider = ["-" * (len(c) + 2 * space) for c in header]
    if justify is not None:
        for i in range(len(divider)):
            d = divider[i]
            if justify[i] == "c":
                d = f":{d[:-2]}:"
            elif justify[i] == "r":
                d = f"{d[:-1]}:"
            elif justify[i] == "l":
                d = f":{d[:-1]}"
            else:
                raise ValueError(f"invalid justfication character '{justify[i]}'")
            divider[i] = d

    divider = f"|{'|'.join(divider)}|"
    rows = [row_template % sep.join(header), divider] + [
        row_template % sep.join(r) for r in formatted_table
    ]
    return "\n".join(rows)


def rst_csv_table(header, formatted_table, title=None, legend=None):
    """Returns a table in restructured text csv-table format

    Parameters
    ----------
    header
        series of strings
    formatted_table
        formatted strings, row based
    title, legend
        combined in this format

    Returns
    -------
    str

    Notes
    -----
    We only support a subset of available attr, see
    https://docutils.sourceforge.io/docs/ref/rst/directives.html#csv-table
    """
    header = ", ".join(f'"{c}"' for c in header)
    header = f"    :header: {header}"
    rows = "\n".join(f"    {', '.join(r)}" for r in formatted_table)

    if title or legend:
        title = f" {title}" if title else ""
        title = f"{title} {legend}" if legend else title
    else:
        title = ""

    table = [f".. csv-table::{title}", header, "", rows]

    return "\n".join(table)


def grid_table_format(header, formatted_table, title=None, legend=None):
    """Returns a table in restructured text grid format.

    Parameters
    ----------
    header
        series with column headings
    formatted_table
        a two dimensional structure (list/tuple) of strings
        previously formatted to the same width within a column.
    title
        optional table title
    legend
        optional table legend

    """
    space = 2
    # make the delineators
    row_delineate = []
    heading_delineate = []
    col_widths = [len(col) for col in header]
    for width in col_widths:
        row_delineate.append("-" * width)
        heading_delineate.append("=" * width)

    row_delineate = "+-" + "-+-".join(row_delineate) + "-+"
    heading_delineate = "+=" + "=+=".join(heading_delineate) + "=+"
    contiguous_delineator = "+" + "-" * (len(row_delineate) - 2) + "+"

    table = []

    # insert the title
    if title:
        table.append(contiguous_delineator)
        if len(title) > len(row_delineate) - 2:
            wrapped = _merged_cell_text_wrap(
                title, len(contiguous_delineator) - 2, space
            )
            for wdex, line in enumerate(wrapped):
                wrapped[wdex] = "|" + line + "|"

            table += wrapped
        else:
            centered = title.center(len(row_delineate) - 2)
            table.append("|" + centered + "|")

    # insert the heading row
    table.append(row_delineate)
    table.append("| " + " | ".join(header) + " |")
    table.append(heading_delineate)

    # concatenate the rows, separating by delineators
    for row in formatted_table:
        table.append("| " + " | ".join(row) + " |")
        table.append(row_delineate)

    if legend:
        if len(legend) > len(row_delineate) - 2:
            wrapped = _merged_cell_text_wrap(
                legend, len(contiguous_delineator) - 2, space
            )
            for wdex, line in enumerate(wrapped):
                wrapped[wdex] = "|" + line + "|"

            table += wrapped
        else:
            ljust = legend.ljust(len(row_delineate) - 3)
            table.append("| " + ljust + "|")

        table.append(contiguous_delineator)

    return "\n".join(table)


def separator_format(header, formatted_table, title=None, legend=None, sep=None):
    """Returns a table with column entries separated by a delimiter. If an entry
    contains the sep character, that entry is put in quotes. Also, title and
    legends (if provided) are forced to a single line and all words forced to
    single spaces.

    Parameters
    ----------
    header
        series with column headings
    formatted_table
        a two dimensional structure (list/tuple) of strings
        previously formatted to the same width within a column.
    sep
        character to separate column entries (eg tab
    title
        optional table title
    legend
        optional table legend

    """
    if sep is None:
        raise RuntimeError("no separator provided")

    if title:
        title = " ".join(" ".join(title.splitlines()).split())

    if legend:
        legend = " ".join(" ".join(legend.splitlines()).split())

    new_table = [sep.join(header)]
    for row in formatted_table:
        for cdex, cell in enumerate(row):
            if sep in cell:
                row[cdex] = f'"{cell}"'

    new_table += [sep.join(row) for row in formatted_table]

    table = "\n".join(new_table)
    # add the title to top of list
    if title:
        table = "\n".join([title, table])
    if legend:
        table = "\n".join([table, legend])

    return table


def format_fields(formats):
    """Formats row fields by index.

    Parameters
    ----------
    formats
        a series consisting of index,formatter callable pairs,
        eg [(0, "'%s'"), (4, '%.4f')]. All non-specified columns are
        formatted as strings.

    """
    index_format = []

    def callable(line, index_format=index_format):
        if not index_format:
            index_format = ["%s" for index in range(len(line))]
            for index, format in formats:
                index_format[index] = format
        formatted = [index_format[i] % line[i] for i in range(len(line))]
        return formatted

    return callable


def separator_formatter(formatter=None, ignore=None, sep=","):
    """Returns a writer for a delimited tabular file. The writer has a
    has_header argument which ignores the formatter for a header line. Default
    format is string. Does not currently handle Titles or Legends.

    Parameters
    ----------
    formatter
        a callable that returns a correctly formatted line.
    ignore
        lines for which ignore returns True are ignored
    sep
        the delimiter deparating fields.

    """
    formatter = formatter or []

    def callable(lines, formatter=formatter, has_header=False):
        if not formatter:
            formatter = format_fields([(i, "%s") for i in range(len(lines[0]))])
        header_done = None
        for line in lines:
            if has_header and not header_done:
                formatted = sep.join([f"{field}" for field in line])
                header_done = True
            else:
                formatted = sep.join(formatter(line))
            yield formatted

    return callable


def formatted_cells(
    rows, header=None, digits=4, column_templates=None, missing_data="", center=False
):
    """Return rows with each columns cells formatted as an equal length
    string.

    Parameters
    ----------
    row
        the series of table rows
    header
        optional header
    digits
        number of decimal places. Can be overridden by following.
    column_templates
        specific format templates for each column.
    missing_data
        default cell value.

    """
    if not header:
        num_col = max(len(row) for row in rows)
        header = [""] * num_col
    else:
        num_col = len(header)

    col_widths = [len(col) for col in header]
    column_templates = column_templates or {}

    float_template = "{0:.%df}" % digits
    # if we have column templates, we use those, otherwise we adaptively
    # apply str/num format
    matrix = []
    for row in rows:
        formatted = []
        for cdex, col_head in enumerate(header):
            try:
                entry = row[cdex]
            except IndexError:
                entry = f"{missing_data}"
            else:
                not_missing = True if isinstance(entry, numpy.ndarray) else entry
                if not not_missing:
                    try:
                        float(entry)  # could numerically be 0, so not missing
                    except (ValueError, TypeError):
                        entry = f"{missing_data}"

            # attempt formatting
            if col_head in column_templates:
                try:  # for functions
                    entry = column_templates[col_head](entry)
                except TypeError:
                    entry = column_templates[col_head] % entry
            elif isinstance(entry, float):
                entry = float_template.format(float(entry))
            else:  # for any other python object
                entry = f"{str(entry)}"

            formatted.append(entry)
            col_widths[cdex] = max(col_widths[cdex], len(entry))
        matrix.append(formatted)

    # now normalise all cell entries to max column widths
    func = {True: lambda x, y: x.center(y)}.get(center, lambda x, y: x.rjust(y))
    new_header = [func(header[i], col_widths[i]) for i in range(num_col)]
    for row in matrix:
        for cdex in range(num_col):
            row[cdex] = func(row[cdex], col_widths[cdex])

    return new_header, matrix


def phylip_matrix(rows, names):
    """Return as a distance matrix in phylip's matrix format."""

    # phylip compatible format is num taxa starting at col 4
    # rows start with taxa names, length 8
    # distances start at 13th col, 2 spaces between each col wrapped
    # at 75th col
    # follow on dists start at col 3
    # outputs a square matrix

    def new_name(names, oldname):
        # the name has to be unique in that number, the best way to ensure that
        # is to determine the number and revise the existing name so it has a
        # int as its end portion
        num = len(names)
        max_num_digits = len(str(num))
        assert max_num_digits < 10, f"can't create a unique name for {oldname}"
        name_base = oldname[: 10 - max_num_digits]
        newname = None
        for i in range(max_num_digits):
            trial_name = f"{name_base}{i}"
            if trial_name not in names:
                newname = trial_name
                break

        if not newname:
            raise RuntimeError(f"Can't create a unique name for {oldname}")
        else:
            print(f"WARN: Seqname {oldname} changed to {newname}")
        return newname

    def append_species(name, formatted_dists, mat_breaks):
        rows = []
        name = name.ljust(12)

        # format the distances first
        for i in range(len(mat_breaks)):
            if i == len(mat_breaks):
                break
            start = mat_breaks[i]
            try:
                end = mat_breaks[i + 1]
            except IndexError:
                end = len(formatted_dists)
            prefix = ["", "  "][i > 0]
            rows.append(f"{prefix}{'  '.join(formatted_dists[start:end])}")
        # mod first row of formatted_dists
        rows[0] = f"{name.ljust(12)}{rows[0]}"
        return rows

    # number of seqs
    numseqs = len(names)

    # determine wrapped table boundaries, if any
    prefix = 13
    mat_breaks = [0]
    line_len = 75  # for the first block
    col_widths = [len(col) for col in rows[0]]
    for i in range(numseqs):
        num_cols = i - mat_breaks[-1]
        if prefix + 2 * num_cols + sum(col_widths[mat_breaks[-1] : i]) > line_len:
            prefix = 3
            line_len = 73
            mat_breaks.append(i)

    # build the formatted distance matrix
    dmat = ["   %d" % numseqs]
    for i in range(numseqs):
        name = names[i].strip()  # we determine white space
        if len(name) > 10:
            name = new_name(names, name)
        dmat += append_species(name, rows[i], mat_breaks)

    return "\n".join(dmat)


def get_continuation_tables_headers(
    cols_widths, index_name=None, space=2, max_width=1e100
):
    """
    returns column headers for continuation tables segmented to not exceed max_width

    Parameters
    ----------
    cols_widths : list
        [[col_name, length of longest string], ...]
    index_name : str
        column name of an index. This column included in all sub table headers.
    space : int
        how much white space between columns
    max_width : int
        maximum width

    Returns
    -------
    list of lists, each inner list is the column names for a subtable
    """
    width_map = dict(cols_widths)
    index_width = 0 if index_name is None else width_map[index_name]
    for name, width in width_map.items():
        if index_width + width > max_width:
            raise ValueError(
                f"{index_name}={index_width} + {name} width={width} > max_width={max_width}"
            )

    if sum(v + space + index_width for _, v in cols_widths) < max_width:
        return [[l for l, _ in cols_widths]]

    headers = []
    curr = [index_name] if index_name is not None else []
    cum_sum = index_width
    for name, width in cols_widths:
        if name == index_name:
            continue

        cum_sum += space + width
        if cum_sum > max_width:
            headers.append(curr)
            curr = [index_name, name] if index_name is not None else [name]
            cum_sum = index_width + space + width
            continue

        curr.append(name)

    headers.append(curr)

    return headers


class _MixedFormatter:
    """handles formatting of mixed data types"""

    def __init__(
        self, alignment, length, precision=4, float_type="f", missing_data=None
    ):
        self.missing_data = missing_data
        self.length = length
        self.alignment = alignment
        self.precision = precision
        self.float_type = float_type

    def __call__(self, val):
        prefix = f"{self.alignment}{self.length}"
        float_spec = f"{prefix}.{self.precision}{self.float_type}"
        int_spec = f"{prefix}d"
        result = str(val)
        if self.missing_data is not None and not result:
            return self.missing_data

        for fspec in (int_spec, float_spec, prefix):
            try:
                result = format(val, fspec)
                break
            except (TypeError, ValueError):
                pass

        return result


def formatted_array(
    series,
    title="",
    precision=4,
    format_spec=None,
    missing_data="",
    pad=True,
    align="r",
):
    """converts elements in a numpy array series to an equal length string.

    Parameters
    ----------
    series
        the series of table rows
    title
        title of series
    precision
        number of decimal places. Can be overridden by following.
    format_spec
        format specification as per the python Format Specification, Mini-Language
        or a callable function.
    missing_data
        default missing data value.
    pad : bool
        Whether to pad all strings to same length. If False, alignment setting is
        ignored.
    align : str
        either 'l', 'c', 'r' for left, center or right alignment, Defaults to 'r'.
        Only applied if pad==True

    Returns
    -------
    list of formatted series, formatted title, maximum string length

    Notes
    -----
    The precedence for formatting is format_spec supersedes pad, precision and
    align values.
    """
    assert isinstance(series, numpy.ndarray), "must be numpy array"
    if pad and align.lower() not in set("lrc"):
        raise ValueError(f"align value '{align}' not in 'l,c,r'")

    if pad:
        align = {"l": "<", "c": "^", "r": ">"}[align]

    if callable(format_spec):
        formatter = format_spec
        format_spec = None
    else:
        formatter = None

    if format_spec and set(format_spec.strip()) <= set("<>^"):
        # format_spec just an alignment character, in which case we assign
        # that to align and reset format_spec as None so other formatting
        # options have an effect
        align = format_spec
        format_spec = None

    if isinstance(format_spec, str):
        format_spec = format_spec.replace("%", "")

    if not any([format_spec, formatter]):
        type_name = series.dtype.name
        if "int" in type_name:
            base_format = "d"
        elif "float" in type_name:
            base_format = f".{precision}f"
        elif "bool" == type_name:
            base_format = ""
        else:
            # handle mixed types with a custom formatter
            formatter = _MixedFormatter(
                align, len(title), precision, missing_data=missing_data
            )
            base_format = ""

        format_spec = base_format

    formatted = []
    max_length = len(title)
    for i, v in enumerate(series):
        if formatter:
            v = formatter(v)
        else:
            try:
                v = format(v, format_spec)
            except (TypeError, ValueError):
                # could be a python object
                v = str(v)

        l = len(v)
        if l > max_length:
            max_length = l

        formatted.append(v)

    if not pad:
        return formatted, title.strip(), max_length

    if format_spec:
        match = re.search("[<>^]", format_spec[:2])
        final_align = align if match is None else match.group()
    else:
        final_align = align

    # now adjust to max_len
    format_spec = f"{final_align}{max_length}s"
    title = format(title, format_spec)
    formatted = [format(v.strip(), format_spec) for v in formatted]
    return formatted, title, max_length


class HtmlElement:
    """wrapper for text to become a HTML element"""

    def __init__(self, text, tag, css_classes=None, newline=False):
        """
        Parameters
        ----------
        text : str
            cell content
        tag : str
            html table cell tag, e.g. 'td', 'th'
        classes : list
            list of custom CSS classes
        newline : bool
            puts the open, close tags on new lines
        """
        self.text = str(text)
        self.tag = tag
        css_classes = [css_classes] if isinstance(css_classes, str) else css_classes
        self.css_classes = css_classes
        self.newline = newline

    def __str__(self):
        txt = self.text
        classes = "" if self.css_classes is None else " ".join(self.css_classes)
        classes = f' class="{classes}"' if classes else ""
        nl = "\n" if self.newline else ""
        return f"{nl}<{self.tag}{classes}>{nl}{txt}{nl}</{self.tag}>"

    def __repr__(self):
        return repr(self.text)


def is_html_markup(text):
    """checks if text contains balanced html markup

    <token ...> body </token>
    """
    pattern = re.compile("(?<=[<])[a-z]+")
    tokens = set(pattern.findall(text))
    if not tokens:
        return False

    for token in tokens:
        num_start = len(re.findall(f"<{token}", text))
        num_end = len(re.findall(f"</{token}", text))
        if num_start != num_end:
            return False

    return True
