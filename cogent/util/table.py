#!/usr/bin/env python
"""
A light-weight Table class for manipulating 2D data and representing it as text, or writing to file for import into other packages.

Current output formats include pickle (pythons serialisation format), restructured text (keyed by 'rest'), latex, html, delimited columns, and a simple text format.

Table can read pickled and delimited formats.

"""
from __future__ import division
import cPickle, csv
from gzip import GzipFile

import numpy
from cogent.format import table as table_format, bedgraph

from cogent.util.dict_array import DictArray

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Felix Schill"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

# making reversed characters for use in reverse order sorting
_all_chrs = [chr(i) for i in range(256)]
_all_chrs.reverse()
_reversed_chrs = ''.join(_all_chrs)

def convert2DDict(twoDdict, header = None, row_order = None):
    """Returns a 2 dimensional list.
    
    Arguments:
        - twoDdict: a 2 dimensional dict with top level keys corresponding to
          column headings, lower level keys correspond to row headings but are
          not preserved.
        - header: series with column headings. If not provided, the sorted top
          level dict keys are used.
        - row_order: a specified order to generate the rows.
    """
    if not header:
        header = twoDdict.keys()
        header.sort()
    
    if not row_order: # we assume rows consistent across dict
        row_order = twoDdict[header[0]].keys()
        row_order.sort()
    
    # make twoD list
    table = []
    for row in row_order:
        string_row = []
        for column in header:
            string_row.append(twoDdict[column][row])
        table.append(string_row)
    return table

class _Header(list):
    """a convenience class for storing the Header"""
    def __new__(cls, arg):
        n = list.__new__(cls, list(arg))
        return n
    
    def __setslice__(self, *args):
        """disallowed"""
        raise RuntimeError("Table Header is immutable, use withNewHeader")
    
    def __setitem__(self, *args):
        """disallowed"""
        raise RuntimeError("Table Header is immutable, use withNewHeader")
    

class Table(DictArray):
    def __init__(self, header = None, rows = None, row_order = None, digits = 4,
                space = 4, title = '', missing_data = '', max_width = 1e100,
                row_ids = False, legend = '', column_templates = None,
                dtype = None):
        """
        Arguments:
        - header: column headings
        - rows: a 2D dict, list or tuple. If a dict, it must have column
          headings as top level keys, and common row labels as keys in each
          column.
        - row_order: the order in which rows will be pulled from the twoDdict
        - digits: floating point resolution
        - space: number of spaces between columns or a string
        - title: as implied
        - missing_data: character assigned if a row has no entry for a column
        - max_width: maximum column width for printing
        - row_ids: if True, the 0'th column is used as row identifiers and keys
          for slicing.
        - legend: table legend
        - column_templates: dict of column headings: string format templates
          or a function that will handle the formatting.
        - dtype: optional numpy array typecode.
        """
        try:
            num_cols = len(header)
            assert num_cols > 0
            if type(rows) == numpy.ndarray:
                assert num_cols == rows.shape[1]
            elif type(rows) == dict:
                assert num_cols == len(rows)
            else:
                assert num_cols == len(rows[0])
        except (IndexError, TypeError, AssertionError):
            raise RuntimeError("header and rows must be provided to Table")
        
        header = [str(head) for head in header]
        if isinstance(rows, dict):
            rows = convert2DDict(rows, header = header,
                                row_order = row_order)
        
        # if row_ids, we select that column as the row identifiers
        if row_ids:
            identifiers = [row[0] for row in rows]
        else:
            identifiers = len(rows)
        
        if not dtype:
            dtype = "O"
        DictArray.__init__(self, rows, identifiers, header, dtype = dtype)
        
        # forcing all column headings to be strings
        self._header = _Header([str(head) for head in header])
        self._missing_data = missing_data
        
        self.Title = str(title)
        self.Legend = str(legend)
        try:
            self.Space = ' ' * space
        except TypeError:
            self.Space = space
        self._digits = digits
        self._row_ids = row_ids
        self._max_width = max_width
        
        # some attributes are not preserved in any file format, so always based
        # on args
        self._column_templates = column_templates or {}
    
    def __repr__(self):
        row = []
        try:
            for val in self.array[0][:3]:
                if isinstance(val, float):
                    row.append('%.4f' % val)
                elif isinstance(val, int):
                    row.append('%d' % val)
                else:
                    row.append('%r' % val)
        except IndexError:
            # an empty table
            pass
        
        row_trunc = ', '.join(row)
        header_trunc = ', '.join(map(repr, self.Header[:3]))
        if self.Shape[1] > 3:
            header_trunc = '[%s,..]' % header_trunc
            row_trunc = '[%s,..]' % row_trunc
        else:
            header_trunc = '[%s]' % header_trunc
            row_trunc = '[%s]' % row_trunc
        
        if self.Shape[0] > 1:
            row_trunc = '[%s,..]' % row_trunc
        else:
            row_trunc = '[[%s]]' % row_trunc
        
        if self.Shape[0] == 0:
            row_trunc = str([])
        
        result = 'Table(numrows=%s, numcols=%s, header=%s, rows=%s)' % \
            (self.Shape[0], self.Shape[1], header_trunc, row_trunc)
        return result
    
    def __str__(self):
        return self.tostring()
    
    def __getitem__(self, names):
        (index, remaining) = self.template.interpretIndex(names)
        # if we have two integers, return a single value
        ints = [isinstance(idx, int) for idx in index]
        if len(ints) == 2 and min(ints):
            return self.array[index]
        new_index = list(index)
        for i, idx in enumerate(new_index):
            if isinstance(idx, int):
                new_index[i] = slice(idx, idx+1, None)
            
        index = tuple(new_index)
        rows = self.array[index]
        result = None
        if len(index) > 1:
            header = numpy.asarray(self.Header, dtype="O")[index[1:]]
        else:
            header = self.Header
        if remaining is not None:
            kwargs = self._get_persistent_attrs()
            result = self.__class__(header, rows, **kwargs)
        return result
    
    def __getstate__(self):
        data = self._get_persistent_attrs()
        del(data['column_templates'])
        data.update(dict(header = self.Header, rows = self.getRawData()))
        return data
    
    def __setstate__(self, data):
        limit_ids = data.pop('limit_ids', None)
        if limit_ids is not None:
            data['row_ids'] = limit_ids or False
        new = Table(**data)
        self.__dict__.update(new.__dict__)
        return self
    
    def _get_header(self):
        """returns Header value"""
        return self._header
    
    def _set_header(self, data):
        """disallowed"""
        raise RuntimeError("not allowed to set the Header")
    
    Header = property(_get_header, _set_header)
    
    def withNewHeader(self, old, new, **kwargs):
        """returns a new Table with old header labels replaced by new
        
        Arguments:
            - old: the old column header(s). Can be a string or series of
              them.
            - new: the new column header(s). Can be a string or series of
              them.
        """
        if type(old) == str:
            old = [old]
            new = [new]
        
        assert len(old) == len(new), 'Mismatched number of old/new labels'
        indices = map(self.Header.index, old)
        new_header = list(self.Header)
        for i in range(len(old)):
            new_header[indices[i]] = new[i]
        
        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header = new_header, rows = self.getRawData(), **kw)
    
    def _get_persistent_attrs(self):
        kws = dict(row_ids = self._row_ids, title = self.Title,
                   legend = self.Legend, digits = self._digits,
                   space = self.Space, max_width = self._max_width,
                   missing_data = self._missing_data,
                   column_templates = self._column_templates or None)
        return kws
    
    def setColumnFormat(self, column_head, format_template):
        """Provide a formatting template for a named column.
        
        Arguments:
            - column_head: the column label.
            - format_template: string formatting template or a function that
              will handle the formatting.
        """
        assert column_head in self.Header, \
               "Unknown column heading %s" % column_head
        
        self._column_templates[column_head] = format_template
    
    def tostring(self, borders=True, sep=None, format='', **kwargs):
        """Return the table as a formatted string.
        
        Arguments:
            - format: possible formats are 'rest', 'latex', 'html', 'phylip',
              'bedgraph', or simple text (default).
            - sep: A string separator for delineating columns, e.g. ',' or '\t'.
              Overrides format.
        
        NOTE: If format is bedgraph, assumes that column headers are chrom,
        start, end, value. In that order!
        """
        if format.lower() == 'phylip':
            missing_data = "%.4f" % 0.0
        else:
            missing_data = self._missing_data
        
        # convert self to a 2D list
        formatted_table = self.array.tolist()
        if format != 'bedgraph':
            header, formatted_table = table_format.formattedCells(formatted_table,
                                    self.Header,
                                    digits = self._digits,
                                    column_templates = self._column_templates,
                                    missing_data = missing_data)
            args = (header, formatted_table, self.Title, self.Legend)
        if sep and format != 'bedgraph':
            return table_format.separatorFormat(*args + (sep,))
        elif format == 'rest':
            return table_format.gridTableFormat(*args)
        elif format.endswith('tex'):
            caption = None
            if self.Title or self.Legend:
                caption = " ".join([self.Title or "", self.Legend or ""])
            return table_format.latex(formatted_table, header,
                                caption = caption, **kwargs)
        elif format == 'html':
            rest = table_format.gridTableFormat(*args)
            return table_format.html(rest)
        elif format == 'phylip':
            # need to eliminate row identifiers
            formatted_table = [row[self._row_ids:] for row in formatted_table]
            header = header[self._row_ids:]
            return table_format.phylipMatrix(formatted_table, header)
        elif format == 'bedgraph':
            assert self.Shape[1] == 4, 'bedgraph format is for 4 column tables'
            # assuming that header order is chrom, start, end, val
            formatted_table = bedgraph.bedgraph(self.sorted().array.tolist(),
                **kwargs)
            return formatted_table
        else:
            return table_format.simpleFormat(*args + (self._max_width,
                                self._row_ids, borders, self.Space))
    
    def toRichHtmlTable(self, row_cell_func=None, header_cell_func=None,
                element_formatters={}, merge_identical=True, compact=True):
        """returns just the table html code.
        Arguments:
            - row_cell_func: callback function that formats the row values. Must
              take the row value and coordinates (row index, column index).
            - header_cell_func: callback function that formats the column headings
              must take the header label value and coordinate
            - element_formatters: a dictionary of specific callback funcs for
              formatting individual html table elements.
              e.g. {'table': lambda x: '<table border="1" class="docutils">'}
            - merge_identical: cells within a row are merged to one span.
        """
        formatted_table = self.array.tolist()
        header, formatted_table = table_format.formattedCells(formatted_table,
                                    self.Header,
                                    digits = self._digits,
                                    column_templates = self._column_templates,
                                    missing_data = self._missing_data)
        # but we strip the cell spacing
        header = [v.strip() for v in header]
        rows = [[c.strip() for c in r] for r in formatted_table]
        return table_format.rich_html(rows, row_cell_func=row_cell_func,
                                      header=header,
                                      header_cell_func=header_cell_func,
                                      element_formatters=element_formatters,
                                      compact=compact)
    
    def writeToFile(self, filename, mode = 'w', writer = None, format = None,
                 sep = None, compress=None, **kwargs):
        """Write table to filename in the specified format. If a format is not
        specified, it attempts to use a filename suffix. Note if a sep argument
        is provided, unformatted values are written to file in order to preserve
        numerical accuracy.
        
        Arguments:
            - mode: file opening mode
            - format: Valid formats are those of the tostring method plus
              pickle.
            - writer: a function for formatting the data for output.
            - sep: a character delimiter for fields.
            - compress: if True, gzips the file and appends .gz to the
              filename (if not already added).
        """
        compress = compress or filename.endswith('.gz')
        
        if compress:
            if not filename.endswith('.gz'):
                filename = '%s.gz' % filename
            mode = ['wb', mode][mode == 'w']
            outfile = GzipFile(filename, mode)
        else:
            outfile = file(filename, mode)
        
        if format is None:
            # try guessing from filename suffix
            if compress:
                index = -2
            else:
                index = -1
            suffix = filename.split('.')
            if len(suffix) > 1:
                format = suffix[index]
        
        if writer:
            rows = self.getRawData()
            rows.insert(0, self.Header[:])
            rows = writer(rows, has_header=True)
            outfile.writelines("\n".join(rows))
        elif format == 'pickle':
            data = self.__getstate__()
            cPickle.dump(data, outfile)
        elif sep is not None and format != 'bedgraph':
            writer = csv.writer(outfile, delimiter = sep)
            if self.Title:
                writer.writerow([self.Title])
            writer.writerow(self.Header)
            writer.writerows(self.array)
            if self.Legend:
                writer.writerow([self.Legend])
        else:
            table = self.tostring(format = format, **kwargs)
            outfile.writelines(table + '\n')
        outfile.close()
    
    def toBedgraph(self, chrom_col, start_col, end_col):
        """docstring for toBedgraph"""
        pass
    
    def appended(self, new_column, *tables, **kwargs):
        """Append an arbitrary number of tables to the end of this one.
        Returns a new table object. Optional keyword arguments to the new
        tables constructor may be passed.
        
        Arguments:
            - new_column: provide a heading for the new column, each tables
            title will be placed in it. If value is false, the result is no
            additional column."""
        
        # convert series of tables
        if isinstance(tables[0], tuple) or isinstance(tables[0], list):
            tables = tuple(tables[0])
        # for each table, determine it's number of rows and create an equivalent
        # length vector of its title
        if new_column:
            header = [new_column] + self.Header
        else:
            header = self.Header
        
        big_twoD = ()
        table_series = (self,) + tables
        for table in table_series:
            # check compatible tables
            assert self.Header == table.Header, \
                   "Inconsistent tables -- column headings are not the same."
            new_twoD = []
            for row in table:
                if new_column:
                    new_twoD.append([table.Title] + row.asarray().tolist())
                else:
                    new_twoD.append(row.asarray().tolist())
            new_twoD = tuple(new_twoD)
            big_twoD += new_twoD
        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header, big_twoD, **kw)
    
    def getRawData(self, columns = None):
        """Returns raw data as a 1D or 2D list of rows from columns. If one
        column, its a 1D list.
        
        Arguments:
            - columns: if None, all data are returned"""
        
        if columns is None:
            return self.array.tolist()
        
        if isinstance(columns, str):
            # assumes all column headings are strings.
            columns = (columns,)
        
        column_indices = map(self.Header.index, columns)
        result = self.array.take(column_indices, axis=1)
        
        if len(columns) == 1:
            result = result.flatten()
        
        return result.tolist()
    
    def _callback(self, callback, row, columns=None, num_columns=None):
        if callable(callback):
            row_segment = row.take(columns)
            if num_columns == 1:
                row_segment = row_segment[0]
            return callback(row_segment)
        else:
            return eval(callback, {}, row)
    
    def filtered(self, callback, columns=None, **kwargs):
        """Returns a sub-table of rows for which the provided callback
        function returns True when passed row data from columns. Row data
        is a 1D list if more than one column, raw row[col] value otherwise.
        
        Arguments:
            - columns: the columns whose values determine whether a row is to
            be included.
            - callback: Can be a function, which takes the sub-row delimited
            by columns and returns True/False, or a string representing valid
            python code to be evaluated."""
        
        if isinstance(columns, str):
            columns = (columns,)
        
        if columns:
            num_columns = len(columns)
        else:
            num_columns = None
        
        row_indexes = []
        if not callable(callback):
            data = self
            cols = columns
        else:
            data = self.array
            cols = map(self.Header.index, columns)
        
        for rdex, row in enumerate(data):
            if self._callback(callback, row, cols, num_columns):
                row_indexes.append(rdex)
            
        sub_set = numpy.take(self, row_indexes, 0)
        
        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header = self.Header, rows = sub_set, **kw)
    
    def filteredByColumn(self, callback, **kwargs):
        """Returns a table with columns identified by callback
        
        Arguments:
            - callback: A function which takes the columns delimited
            by columns and returns True/False, or a string representing valid
            python code to be evaluated."""
        data = self.array.transpose()
        column_indices = []
        append = column_indices.append
        for index, row in enumerate(data):
            if callback(row):
                append(index)
        columns = numpy.take(self.Header, column_indices)
        return self.getColumns(columns, **kwargs)
    
    def count(self, callback, columns=None, **kwargs):
        """Returns number of rows for which the provided callback
        function returns True when passed row data from columns. Row data
        is a 1D list if more than one column, raw row[col] value otherwise.
        
        Arguments:
            - columns: the columns whose values determine whether a row is to
            be included.
            - callback: Can be a function, which takes the sub-row delimited
            by columns and returns True/False, or a string representing valid
            python code to be evaluated."""
        
        if isinstance(columns, str):
            columns = (columns,)
        
        if columns:
            num_columns = len(columns)
        else:
            num_columns = None
        
        count = 0
        if not callable(callback):
            data = self
            cols = columns
        else:
            data = self.array
            cols = map(self.Header.index, columns)
        
        for row in data:
            if self._callback(callback, row, cols, num_columns):
                count += 1
        return count
    
    def sorted(self, columns = None, reverse = None, **kwargs):
        """Returns a new table sorted according to columns order.
        
        Arguments:
            - columns: column headings, their order determines the sort order.
            - reverse: column headings, these columns will be reverse sorted.
            
            Either can be provided as just a single string, or a series of
            strings.
        """
        
        if columns is None:
            columns = self.Header
        elif isinstance(columns, str):
            columns = [columns]
        
        indices = [self.Header.index(col) for col in columns]
        
        if not reverse:
            is_reversed = [False] * len(columns)
            reverse_indices = []
        else:
            if type(reverse) == str:
                reverse = [reverse]
            reverse_indices = []
            for index, header_index in enumerate(indices):
                col = self.Header[header_index]
                if col in reverse:
                    reverse_indices += [index]
            
            is_reversed = [col in reverse for col in columns]
        
        reverse_indices = numpy.array(reverse_indices)
        
        dtypes = [(self.Header[i], self.array.dtype) for i in indices]
        
        # applying the decorate-sort-undecorate approach
        aux_list = self.array.take(indices, axis=1)
        
        # we figure out the casting funcs for any reversed elements
        cast = []
        for index in reverse_indices:
            val = aux_list[0, index]
            try:
                val = val.translate(_reversed_chrs)
                func = lambda x: x.translate(_reversed_chrs)
            except AttributeError:
                func = lambda x: x * -1
            func = numpy.vectorize(func)
            aux_list[:, index] = func(aux_list[:, index])
        
        aux_list = numpy.rec.fromarrays(aux_list.copy().T, dtype=dtypes)
        indices = aux_list.argsort()
        new_twoD = self.array.take(indices, axis=0)
        
        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header = self.Header, rows = new_twoD, **kw)
    
    def getColumns(self, columns, **kwargs):
        """Return a slice of columns"""
        # check whether we have integer columns
        
        if isinstance(columns, str):
            columns = [columns]
        
        is_int = min([isinstance(val, int) for val in columns])
        indexes = []
        if is_int:
            indexes = columns
        else:
            indexes = [self.Header.index(head) for head in columns]
        
        if self._row_ids:
            # we disallow reordering of identifiers, and ensure they are only
            # presented once
            for val in range(self._row_ids):
                try:
                    indexes.remove(val)
                except ValueError:
                    pass
            indexes = range(self._row_ids) + indexes
        
        columns = numpy.take(numpy.asarray(self.Header, dtype="O"),
                               indexes)
        new = numpy.take(self.array, indexes, axis=1)
        
        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header = columns, rows = new, **kw)
    
    def getDisjointRows(self, rows, **kwargs):
        """Return the nominated disjoint rows."""
        if isinstance(rows, str):
            rows = [rows]
        
        indexes = []
        for row in rows:
            idx, drop = self.template.interpretIndex(row)
            indexes.append(idx[0])
        
        new = self.array.take(indexes, axis=0)
        
        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header = self.Header, rows = new, **kw)
    
    def withNewColumn(self, new_column, callback, columns = None, **kwargs):
        """Returns a new table with an additional column, computed using
        callback.
        
        Arguments:
            - new_column: new column heading
            - columns: the columns whose values determine whether a row is to
            be included.
            - callback: Can be a function, which takes the sub-row delimited
            by columns and returns True/False, or a string representing valid
            python code to be evaluated."""
        
        if isinstance(columns, str):
            columns = (columns,)
        
        if columns is not None:
            num_columns = len(columns)
        else:
            num_columns = None
        
        if not callable(callback):
            data = self
            cols = columns
        else:
            data = self.array
            cols = map(self.Header.index, columns)
        
        twoD = [list(row) + [self._callback(callback, row, cols,
                num_columns)] for row in data]
        
        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header = self.Header + [new_column], rows = twoD, **kw)
    
    def getDistinctValues(self, column):
        """returns the set of distinct values for the named column(s)"""
        columns = [column, [column]][type(column) == str]
        data = self.getRawData(column)
        
        if len(columns) > 1:
            data = [tuple(row) for row in data]
        
        return set(data)
    
    def joined(self, other_table, columns_self=None, columns_other=None,
                inner_join=True, **kwargs):
        """returns a new table containing the join of this table and
        other_table. Default behaviour is the natural inner join. Checks for
        equality in the specified columns (if provided) or all columns; a
        combined row is included in the output if all indices match exactly. A
        combined row contains first the row of this table, and then columns
        from the other_table that are not key columns (i.e. not specified in
        columns_other). The order (of self, then other)
        is preserved. The column headers of the output are made unique by
        replacing the headers of other_table with
        <other_table.Title>_<other_table.Header>.
        
        Arguments:
            - other_table: A table object which will be joined with this
              table. other_table must have a title.
            - columns_self, columns_other: indices of key columns that will
              be compared in the join operation. Can be either column index,
              or a string matching the column header. The order matters, and
              the dimensions of columns_self and columns_other have to match.
              A row will be included in the output iff
              self[row][columns_self[i]]==other_table[row][columns_other[i]]
              for all i
            - inner_join: if False, the outer join of the two tables is
              returned.
        """
        
        if other_table.Title is None:
            raise RuntimeError, "Cannot join if a other_table.Title is None"
        elif self.Title == other_table.Title:
            raise RuntimeError, "Cannot join if a table.Title's are equal"
        
        columns_self = [columns_self,[columns_self]][type(columns_self)==str]
        columns_other = [columns_other,
                            [columns_other]][type(columns_other)==str]
        if not inner_join:
            assert columns_self is None and columns_other is None, "Cannot "\
                                "specify column indices for an outer join"
            columns_self = []
            columns_other = []
        
        if columns_self is None and columns_other is None:
            # we do the natural inner join
            columns_self=[]
            columns_other=[]
            for col_head in self.Header:
                if col_head in other_table.Header:
                    columns_self.append(self.Header.index(col_head))
                    columns_other.append(other_table.Header.index(col_head))
        elif columns_self is None or columns_other is None:
            # the same column labels will be used for both tables
            columns_self = columns_self or columns_other
            columns_other = columns_self or columns_other
        elif len(columns_self)!=len(columns_other):
            raise RuntimeError("Error during table join: key columns have "\
                  "different dimensions!")
        
        # create new 2d list for the output
        joined_table=[]
        
        #resolve column indices from Header, if necessary
        columns_self_indices=[]
        columns_other_indices=[]
        for col in columns_self:
            if type(col)==int:
                columns_self_indices.append(col)
            else:
                columns_self_indices.append(self.Header.index(col))
        
        for col in columns_other:
            if type(col)==int:
                columns_other_indices.append(col)
            else:
                columns_other_indices.append(other_table.Header.index(col))
        # create a mask of which columns of the other_table will end up in the
        # output
        output_mask_other=[]
        for col in range(0,len(other_table.Header)):
            if not (col in columns_other_indices):
                output_mask_other.append(col)
        # use a dictionary for the key lookup
        # key dictionary for other_table.
        # key is a tuple made from specified columns; data is the row index
        # for lookup...
        key_lookup={}
        row_index=0
        for row in other_table:
            #insert new entry for each row
            key=tuple([row[col] for col in columns_other_indices])
            if key in key_lookup:
                key_lookup[key].append(row_index)
            else:
                key_lookup[key]=[row_index]
            row_index=row_index+1
            
        for this_row in self:
            # assemble key for query of other_table
            key=tuple([this_row[col] for col in columns_self_indices])
            if key in key_lookup:
                for output_row_index in key_lookup[key]:
                    other_row=[other_table[output_row_index,c] \
                                                  for c in output_mask_other]
                    joined_table.append(list(this_row) + other_row)
        
        new_header=self.Header+[other_table.Title+"_"+other_table.Header[c] \
                                                  for c in output_mask_other]
        if not joined_table:
            # YUK, this is to stop dimension check in DictArray causing
            # failures
            joined_table = numpy.empty((0,len(new_header)))
        
        return Table(header=new_header, rows=joined_table, **kwargs)
    
    def summed(self, indices=None, col_sum=True, strict=True, **kwargs):
        """returns the sum of numerical values for column(s)/row(s)
        
        Arguments:
            - indices: column name(s) or indices or row indices
            - col_sum: sums values in the indicated column, the default. If
              False, returns the row sum.
            - strict: if False, ignores cells with non-numeric data in the
              column/row."""
        
        all = indices is None
        
        if type(indices) == str:
            assert col_sum, "Must use row integer indices"
            indices = self.Header.index(indices)
        elif type(indices) == int: # a single index
            indices = [indices]
        elif not all:
            raise RuntimeError("unknown indices type: %s" % indices)
        
        if not all:
            vals = self.array.take([indices], axis=[0,1][col_sum]).flatten()
            if strict:
                result = vals.sum()
            else:
                result = sum(v for v in vals if type(v)!=str)
        else:
            # a multi-rowed result
            if col_sum:
                vals = self.array
            else:
                vals = self.array.transpose() 
            
            if strict:
                result = vals.sum(axis=0).tolist()
            else:
                result = []
                append = result.append
                # we need to iterate over the elements to be summed, so we
                # have to transpose
                for row in vals.transpose():
                    try:
                        num = row.sum()
                    except TypeError:
                        num = sum(r for r in row if type(r) != str)
                    append(num)
            
        
        return result
    
    def normalized(self, by_row=True, denominator_func=None, **kwargs):
        """returns a table with elements expressed as a fraction according
        to the results from func
        
        Arguments:
            - by_row: normalisation done by row
            - denominator_func: a callback function that takes an array and
              returns a value to be used as the denominator. Default is sum."""
        
        if denominator_func:
            data = self.array
            if not by_row:
                data = data.transpose()
            denominators = [denominator_func(row) for row in data]
        else:
            denominators = self.summed(col_sum=not by_row)
        
        if by_row:
            values = self.array
        else:
            values = self.array.transpose()
        
        rows = [values[i]/denom for i, denom in enumerate(denominators)]
        rows = numpy.array(rows)
        
        if not by_row:
            rows = rows.transpose()
        
        return Table(header=self.Header, rows=rows, **kwargs)
    
    def transposed(self, new_column_name, select_as_header=None, **kwargs):
        """returns the transposed table.
        
        Arguments:
            - new_column_name: the existing header will become a column with
              this name
            - select_as_header: current column name containing data to be used
              as the header. Defaults to the first column.
        """
        select_as_header = select_as_header or self.Header[0]
        assert select_as_header in self.Header, \
                    '"%s" not in table Header' % select_as_header
        
        raw_data = self.getRawData()
        raw_data.insert(0, self.Header)
        transposed = numpy.array(raw_data, dtype='O')
        transposed = transposed.transpose()
        
        # indices for the header and non header rows
        header_index = self.Header.index(select_as_header)
        
        data_indices = range(0, header_index)+range(header_index+1,
                                                    len(transposed))
        
        header = list(numpy.take(transposed, [header_index], axis=0)[0])
        header = [new_column_name]+header[1:] # [1:] slice excludes old name
        rows = numpy.take(transposed, data_indices, axis=0)
        return Table(header=header, rows=rows, **kwargs)
    
