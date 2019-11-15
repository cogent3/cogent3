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

from cogent3.util.warning import discontinued


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell", "Matthew Wakefield", "Jeremy Widmann"]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
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


def _merged_cell_text_wrap(text, max_line_length, space):
    """ left justify wraps text into multiple rows"""
    max_line_width = max_line_length - (2 * space)
    if len(text) < max_line_length:
        return [text]
    buffer = " " * space
    wrapped = textwrap.wrap(
        text, width=max_line_width, initial_indent=buffer, subsequent_indent=buffer
    )
    wrapped = ["%s" % line.ljust(max_line_width + 2 * space) for line in wrapped]
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
            return "<td>%s</td>" % v

    if header_cell_func is None:

        def header_cell_func(v, c):
            return "<th>%s</th>" % v

    if merge_identical:
        row_iterator = _merge_cells
    else:
        row_iterator = enumerate

    if header:
        thead = formatted("thead", '<thead style="font-weight: bold;">')
        row = [header_cell_func(escape(label), i) for i, label in enumerate(header)]
        data += [thead] + row + ["</thead>"]

    formatted_rows = []
    td = formatted("td", "<td>")
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


def latex(rows, header=None, caption=None, justify=None, label=None, position=None):
    """Returns the text a LaTeX table.

    Parameters
    ----------
    header
        table header
    position
        table page position, default is here, top separate page
    justify
        column justification, default is right aligned.
    caption
        Table legend
    label
        for cross referencing

    """

    if not justify:
        numcols = [len(header), len(rows[0])][not header]
        justify = "r" * numcols

    justify = "{ %s }" % " ".join(list(justify))
    if header:
        header = "%s \\\\" % " & ".join([r"\bf{%s}" % head.strip() for head in header])
    rows = ["%s \\\\" % " & ".join(row) for row in rows]
    position = position or "htp!"
    table_format = [
        r"\begin{table}[%s]" % position,
        r"\centering",
        r"\begin{tabular}%s" % justify,
    ]
    table_format.append(r"\hline")
    table_format.append(header)
    table_format.append(r"\hline")
    table_format.append(r"\hline")
    table_format += rows
    table_format.append(r"\hline")
    table_format.append(r"\end{tabular}")
    if caption:
        table_format.append(r"\caption{%s}" % caption)
    if label:
        table_format.append(r"\label{%s}" % label)
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

    if not identifiers:
        identifiers = 0
    # having determined the maximum string lengths we now need to
    # produce subtables of width <= max_width
    col_widths = [len(head) for head in header]
    sep = len(space)
    min_length = sep * (identifiers - 1) + sum(col_widths[:identifiers])

    if min_length > max_width:
        raise RuntimeError("Maximum width too small for identifiers")

    begin, width = identifiers, min_length

    boundaries = []
    for i in range(begin, len(header)):
        width += col_widths[i] + sep
        if width > max_width:
            boundaries.append((begin, i, width - col_widths[i] - sep))
            width = min_length + col_widths[i] + sep
            begin = i

    # add the last sub-table
    boundaries.append((begin, len(header), width))
    # generate the table
    for start, end, width in boundaries:
        subhead = header[:identifiers] + header[start:end]
        rows = [row[:identifiers] + row[start:end] for row in formatted_table]
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
        column index for the last column that uniquely identify
        rows. Required if table width exceeds max_width.
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
                d = ":%s:" % d[:-2]
            elif justify[i] == "r":
                d = "%s:" % d[:-1]
            elif justify[i] == "l":
                d = ":%s" % d[:-1]
            else:
                raise ValueError("invalid justfication character '%s'" % justify[i])
            divider[i] = d

    divider = "|%s|" % "|".join(divider)
    rows = [row_template % sep.join(header), divider] + [
        row_template % sep.join(r) for r in formatted_table
    ]
    return "\n".join(rows)


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
                row[cdex] = '"%s"' % cell

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
                formatted = sep.join(["%s" % field for field in line])
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
        num_col = max([len(row) for row in rows])
        header = [""] * num_col
    else:
        num_col = len(header)

    col_widths = [len(col) for col in header]
    num_row = len(rows)
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
                entry = "%s" % missing_data
            else:
                if not entry:
                    try:
                        float(entry)  # could numerically be 0, so not missing
                    except (ValueError, TypeError):
                        entry = "%s" % missing_data

            # attempt formatting
            if col_head in column_templates:
                try:  # for functions
                    entry = column_templates[col_head](entry)
                except TypeError:
                    entry = column_templates[col_head] % entry
            elif isinstance(entry, float):
                entry = float_template.format(float(entry))
            else:  # for any other python object
                entry = "%s" % str(entry)

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
        assert max_num_digits < 10, "can't create a unique name for %s" % oldname
        name_base = oldname[: 10 - max_num_digits]
        newname = None
        for i in range(max_num_digits):
            trial_name = "%s%s" % (name_base, i)
            if trial_name not in names:
                newname = trial_name
                break

        if not newname:
            raise RuntimeError("Can't create a unique name for %s" % oldname)
        else:
            print("WARN: Seqname %s changed to %s" % (oldname, newname))
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
            rows.append("%s%s" % (prefix, "  ".join(formatted_dists[start:end])))
        # mod first row of formatted_dists
        rows[0] = "%s%s" % (name.ljust(12), rows[0])
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
