#!/usr/bin/env python
"""Dict2D: holds two-dimensional dict, acting as matrix with labeled rows/cols.

The Dict2D is useful for storing arbitrary, sparse data in a way that is easy
to access by labeled rows and columns. It is much slower than a numpy 
array, so only use when the convenience outweighs the performance penalty.
It is especially useful for storing distance matrices between arbitrarily 
labeled objects.
"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class Dict2DError(Exception):
    """All Dict2D-specific errors come from here."""
    pass

class Dict2DInitError(ValueError, Dict2DError):
    """Raised if Dict2D init fails."""
    pass

class Dict2DSparseError(KeyError, Dict2DError):
    """Raised on operations that fail because the Dict2D is sparse."""
    pass


## Following methods based on methods developed by Rob Knight
## These methods are intended to be used by the reflect function of 
## SquareMatrix. Each method takes two values, and returns two values. The 
## idea is to reflect a matrix you take the value from the upper triangle, 
## and the value from the lower triangle perform the intended operation on 
## them, and then return the new values for the upper triangle and the lower
## triangle.
def average(upper, lower):
    """Returns mean of the two values."""
    try:
        val = (upper + lower)/2.0
        return val, val
    except TypeError:
        raise TypeError, "%s or %s invalid types for averaging."\
                % (str(upper), str(lower))

def largest(upper, lower):
    """Returns largest of the two values."""
    val = max(upper, lower)
    return val, val

def smallest(upper, lower):
    """Returns smallest of the two values."""
    val = min(upper, lower)
    return val, val

def swap(upper, lower):
    """Swaps the two values."""
    return lower, upper

def nonzero(upper, lower):
    """Fills both values to whichever evaluates True, or leaves in place."""
    if upper and not lower:
        return upper, upper
    elif lower and not upper:
        return lower, lower
    else:
        return upper, lower

def not_0(upper, lower):
    """Fills both values to whichever is not equal to 0, or leaves in place."""
    # Note that I compare the values to zero istead of using if not lower
    # This is because the values coming in can be anything, including an
    # empty list, which should be considered nonzero
    if upper == 0:
        return lower, lower
    elif lower == 0:
        return upper, upper
    else:
        return upper, lower
     
def upper_to_lower(upper, lower):
    """ return new symm matrix with upper tri copied to lower tri"""
    return upper, upper

def lower_to_upper(upper, lower):
    """ return new symm matrix with upper tri copied to lower tri"""
    return lower, lower
# End methods developed by Rob Knight
    

class Dict2D(dict):
    """Implements dict of dicts with expanded functionality
    
        This class is useful for creating and working with 2D dictionaries.
        
    """
    #override these in subclasses for convenient customization.
    RowOrder = None         #default RowOrder
    ColOrder = None         #default ColOrder
    Default = None          #default Default value when m[r][c] absent
    RowConstructor = dict   #default row constructor
    Pad = False             #default state for automatic col and item padding
    
    #list of attributes that is copied by the copy() method
    _copied_attributes = ['RowOrder', 'ColOrder', 'Default', 'Pad', \
        'RowConstructor']
    
    def __init__(self, data=None, RowOrder=None, ColOrder=None, Default=None,
        Pad=None, RowConstructor=None):
        """Returns new Dict2D with specified parameters.
        
            data : can either be a dict of dicts, or a sequence of 3-item
            sequences giving row, col, value.
            
            RowOrder: list of 'interesting' row keys. If passed in during
            init, all rows in RowOrder will be created. Rows not in RowOrder
            will not be printed or used in most calculations, if they exist.
            Default is None (calculates on the fly from self.keys().
            
            ColOrder: list of 'interesting' column keys. If passed in 
            during init, all columns in ColOrder will be created in all rows.
            Columns not in ColOrder will not be printed or used in most
            calculations, if they exist. Default is None (calculates on the
            fly by examining the keys in each row. This can be expensive!
            
            Default: value returned when m[r][c] doesn't exist.

            Pad: whether or not to pad Cols and Items with the default value 
            instead of raising an exception. Default False.

            RowConstructor: constructor for inner rows. Defaults to class
            value of dict. WARNING: Must be able to accept initialization
            with no parameters (i.e. if you have a class that requires
            parameters to initialize, RowConstructor should be a function
            that supplies all the appropriate defaults.)

            Note that all operations that alter the Dict2D act in place. If 
            you want to operate on a different object you should call the 
            Dict2D copy() to create an identical deep copy of your Dict2D
            and then work on that one, leaving the original untouched.
            See doc string for Dict2D.copy() for usage information.
            
            usage: 
                d = {'a':{'x':1,'y':2}, 'b':{'x':0, 'z':5}}
                m = Dict2D(d)
                m = Dict2D(d,Rows=['a','b'])
                m = Dict2D(d,Cols=['x'],Default=99)
        
        """
        #set the following as instance data if supplied; otherwise, will
        #fall through to class data
        if RowOrder is not None:
            self.RowOrder = RowOrder
        if ColOrder is not None:
            self.ColOrder = ColOrder
        if Default is not None:
            self.Default = Default
        if Pad is not None:
            self.Pad = Pad
        if RowConstructor is not None:
            self.RowConstructor = RowConstructor
       
        #initialize data as an empty dict if data = None
        data = data or {}
       
        init_method = self._guess_input_type(data)
        if not init_method:
            raise Dict2DInitError, \
            "Dict2D init failed (data unknown type, or Row/Col order needed)."
        #if we get here, we got an init method that it's safe to call
        init_method(data)
        #fill in any missing m[r][c] from RowOrder and ColOrder if self.Pad.
        if self.Pad:
            self.pad()

    def _guess_input_type(self, data):
        """Guesses the input type of data, and returns appropriate init method.

        Known init methods are fromDicts, fromIndices, and fromLists.
        Returns None if it can't figure out the data type.
        """
        if isinstance(data, dict):
            #assume dict of dicts
            return self.fromDicts

        else:
            RowOrder = self.RowOrder
            ColOrder = self.ColOrder
            try:
                #guess list of lists or seq of 3-item seqs
                if (RowOrder is not None) and \
                   (ColOrder is not None) and \
                   (len(data) == len(RowOrder)) and \
                   (len(data[0]) == len(ColOrder)):
                       #assume list of lists
                    return self.fromLists
                elif len(data[0]) == 3:
                #assume seq of 3-item seqs
                    return self.fromIndices
            except:
                #if there's any exception, we guessed the wrong type so
                #will return None
                pass

    def fromDicts(self, data):
        """Fills self from dict of dicts."""
        constructor = self.RowConstructor
        try:
            for key, val in data.items():
                self[key] = constructor(val)
        except (TypeError, ValueError, AttributeError):
            raise Dict2DInitError, \
            "Dict2D init from dicts failed."

    def fromIndices(self, data):
        """Fills self from sequence of (row, col, value) sequences."""
        constructor = self.RowConstructor
        try:
            for row, col, val in data:
                curr_row = self.get(row, constructor())
                curr_row[col] = val
                self[row] = curr_row
        except (TypeError, ValueError, AttributeError):
            raise Dict2DInitError, \
            "Dict2D init from indices failed."

    def fromLists(self, data):
        """Fills self from list of lists.

        Note that dimensions of list of lists must match RowOrder x ColOrder."""
        constructor = self.RowConstructor
        if (self.RowOrder is None) or (self.ColOrder is None):
            raise Dict2DInitError, \
            "Must have RowOrder and ColOrder to init Dict2D from list of lists."
        try:
            for key, row in zip(self.RowOrder, data):
                self[key] = dict(zip(self.ColOrder, row))
        except (TypeError):
            raise Dict2DInitError, \
            "Dict2D init from lists failed."
            
    def pad(self, default=None):
        """Ensures self[r][c] exists for r in RowOrder for c in ColOrder.

        default, if not specified, uses self.Default.
        """
        constructor = self.RowConstructor 
        if default is None:
            default = self.Default

        row_order = self.RowOrder or self.rowKeys()
        col_order = self.ColOrder or self.colKeys()
            
        for r in row_order:
            if r not in self:
                self[r] = constructor()
            curr_row = self[r]
            for c in col_order:
                if c not in curr_row:
                    curr_row[c] = default

    def purge(self):
        """Keeps only items self[r][c] if r in RowOrder and c in ColOrder."""
        #first, purge unwanted rows
        if self.RowOrder:
            wanted_keys = dict.fromkeys(self.RowOrder)
            for key in self.keys():
                if not key in wanted_keys:
                    del self[key]
        #then, purge unwanted cols
        if self.ColOrder:
            wanted_keys = dict.fromkeys(self.ColOrder)
            for row in self.values():
                for key in row.keys():
                    if not key in wanted_keys:
                        del row[key]

    def rowKeys(self):
        """Returns list of keys corresponding to all rows.

        Same as list(self).
        """
        return list(self)

    def colKeys(self):
        """Returns list of keys corresponding to all cols."""
        result = {}
        for row in self.values():
            result.update(row)
        return list(result)

    def sharedColKeys(self):
        """Returns list of keys shared by all cols."""
        rows = self.values()
        if not rows:
            return []
        result = rows[0]
        for row in rows:
            for key in result.keys():
                if key not in row:
                    del result[key]
        return list(result)
        

    def square(self, default=None, reset_order=False):
        """Checks RowOrder and ColOrder share keys, and that self[r][c] exists.

        If reset_order is True (default is False), appends additional Cols to 
        RowOrder and sets ColOrder to RowOrder.
        """
        row_order = self.RowOrder or self.rowKeys()
        col_order = self.ColOrder or self.colKeys()
        rows = dict.fromkeys(row_order)
        cols = dict.fromkeys(col_order)
        
        if reset_order:
            if rows != cols:
                for c in cols:
                    if c not in rows:
                        row_order.append(c)
                self.RowOrder = row_order
                #WARNING: we set self.ColOrder to row_order as well, _not_ to
                #col_order, because we want the RowOrder and ColOrder to be
                #the same after squaring.
                self.ColOrder = row_order
        else:
            if rows != cols:
                raise Dict2DError, \
                "Rows and Cols must be the same to square a Dict2D."
        self.pad(default)
            
    def _get_rows(self):
        """Iterates over the rows, using RowOrder/ColOrder.

        Converts the rows to lists of values, so r[i] is the same as
        m[r][m.ColOrder.index(c)] (assuming that the items in ColOrder are
        unique). zip(self.ColOrder, row) will associate the column label
        with each item in the row.

        If you actually want to get a list of the row objects, you probably
        want self.values() or [self[i] for i in self.RowOrder] instead of
        this method.

        If self.Pad is True, will pad rows with self.Default instead of 
        raising an exception.
        """
        row_order = self.RowOrder or self.rowKeys()
        
        if self.Pad:
            col_order = self.ColOrder or self.colKeys()
            
            constructor = self.RowConstructor
            default = self.Default
            for r in row_order:
                curr_row = self.get(r, constructor())
                yield [curr_row.get(c, default) for c in col_order]
        else:
            col_order = self.ColOrder
            if col_order:   #need to get items into the right column order
                try:
                    for r in row_order:
                        curr_row = self[r]
                        yield [curr_row[c] for c in col_order]
                except KeyError:
                    raise Dict2DSparseError, \
                    "Can't iterate over rows of sparse Dict2D."
            else:           #if there's no ColOrder, just return what's there
                for r in row_order:
                    curr_row = self[r]
                    yield curr_row.values()
                
    Rows = property(_get_rows)

    def _get_cols(self):
        """Iterates over the columns, using RowOrder/ColOrder.

        Returns each column as a list of the values in that column, so that
        c[i] = m[m.RowOrder.index(r)][c] (assuming the items in RowOrder are 
        unique).

        zip(self.RowOrder, col) will associate the row label with each item
        in the column.

        If you want to get the column objects as dicts that support named
        lookups, so that c[r] = m[r][c], your best bet is something like:
            cols = self.copy()
            cols.transpose()
            return cols.values()    #or [cols[r] for r in cols.RowOrder]

        Will fail if ColOrder is specified and keys are missing.
        """
        row_order = self.RowOrder or self.rowKeys()
        col_order = self.ColOrder or self.colKeys()
        
        if self.Pad:
            default = self.Default
            constructor = self.RowConstructor
            for c in col_order:
                yield [self.get(r, constructor()).get(c, default) for \
                    r in row_order]
        else:
            try:
                for c in col_order:
                    yield [self[r][c] for r in row_order]
            except KeyError:
                raise Dict2DSparseError, \
                "Can't iterate over cols of sparse Dict2D."
    
    Cols = property(_get_cols)

    def _get_items(self):
        """Iterates over the items, using RowOrder and ColOrder if present.

        Returns a list of the items, by rows rather than by columns.

        self.Pad controls whether to insert the default anywhere a value is
        missing, or to return only the values that exist.
        """
        for row in self.Rows:
            for i in row:
                yield i

    Items = property(_get_items)

    def getRows(self, rows, negate=False):
        """Returns new Dict2D containing only specified rows.
        
        Note that the rows in the new Dict2D will be references to the
        same objects as the rows in the old Dict2D.

        If self.Pad is True, will create new rows rather than raising 
        KeyError.
        """
        result = {}
        if negate:
            #copy everything except the specified rows
            row_lookup = dict.fromkeys(rows)
            for r, row in self.items():
                if r not in row_lookup:
                    result[r] = row
        else:
            #copy only the specified rows
            if self.Pad:
                row_constructor = self.RowConstructor
                for r in rows:
                    result[r] = self.get(r, row_constructor())
            else:
                for r in rows:
                    result[r] = self[r]
        return self.__class__(result)

    def getRowIndices(self, f, negate=False):
        """Returns list of keys of rows where f(row) is True.
        
        List will be in the same order as self.RowOrder, if present.
        Note that the function is applied to the row as given by self.Rows,
        not to the original dict that contains it.
        """
        #negate function if necessary
        if negate:
            new_f = lambda x: not f(x)
        else:
            new_f = f
        #get all the rows where the function is True
        row_order = self.RowOrder or self 
        return [key for key, row in zip(row_order,self.Rows) \
            if new_f(row)]

    def getRowsIf(self, f, negate=False):
        """Returns new Dict2D containing rows where f(row) is True.
        
        Note that the rows in the new Dict2D are the same objects as the
        rows in the old Dict2D, not copies.
        """
        #pass negate to get RowIndices
        return self.getRows(self.getRowIndices(f, negate))

    def getCols(self, cols, negate=False, row_constructor=None):
        """Returns new Dict2D containing only specified cols.
        
        By default, the rows will be dicts, but an alternative constructor
        can be specified.

        Note that getCols should not fail on ragged columns, and will just
        ignore any elements that are not explicitly present in a given row
        whether or not self.Pad is set.
        """
        if row_constructor is None:
            row_constructor = self.RowConstructor
        result = {}
        #if we're negating, pick out all the columns except specified indices
        if negate:
            col_lookup = dict.fromkeys(cols)
            for key, row in self.items():
                result[key] = row_constructor([(i, row[i]) for i in row \
                if (i in row) and (i not in col_lookup)])
        #otherwise, just get the requested indices
        else:
            for key, row in self.items():
                result[key] = row_constructor([(i, row[i]) for i in cols \
                if i in row])
        return self.__class__(result)

    def getColIndices(self, f, negate=False):
        """Returns list of column indices for which f(col) is True."""
        #negate f if necessary
        if negate:
            new_f = lambda x: not f(x)
        else:
            new_f = f
        return [i for i, col in zip(self.ColOrder or self.colKeys(), self.Cols)\
            if new_f(col)]

    def getColsIf(self, f, negate=False, row_constructor=None):
        """Returns new Dict2D containing cols where f(col) is True.

        Note that the rows in the new Dict2D are always new objects. Default
        constructor is list(), but an alternative can be passed in.
        """
        if row_constructor is None:
            row_constructor = self.RowConstructor
        return self.getCols(self.getColIndices(f, negate), \
            row_constructor=row_constructor)

    def getItems(self, items, negate=False):
        """Returns list containing only specified items.
        
        items should be a list of (row_key, col_key) tuples.

        getItems will fail with KeyError if items that don't exist are
        requested, unless self.Pad is True.

        Items will be returned in order (according to self.ColOrder and 
        self.RowOrder) when negate is True; when negate is False, they'll
        be returned in the order in which they were passed in.
        """
        if negate:
            #have to cycle through every item and check that it's not in
            #the list of items to return
            item_lookup = dict.fromkeys(map(tuple, items))
            result = []
            if self.Pad:
                default = self.Default
                row_constructor = self.RowConstructor
                for r in self.RowOrder or self:
                    curr_row = self.get(r, row_constructor())
                    for c in self.ColOrder or curr_row:
                        if (r, c) not in items:
                            result.append(curr_row.get(c, default))
            else:
                for r in self.RowOrder or self:
                    curr_row = self[r]
                    for c in self.ColOrder or curr_row:
                        if c in curr_row and (r, c) not in items:
                            result.append(curr_row[c])
            return result
        #otherwise, just pick the selected items out of the list
        else:
            if self.Pad:
                row_constructor = self.RowConstructor
                default = self.Default
                return [self.get(row, row_constructor()).get(col, default) \
                    for row, col in items]
            else:
                return [self[row][col] for row, col in items]

    def getItemIndices(self, f, negate=False):
        """Returns list of (key,val) tuples where f(self[key][val]) is True."""
        if negate:
            new_f = lambda x: not f(x)
        else:
            new_f = f
        result = [] 
        for row_label in self.RowOrder or self:
            curr_row = self[row_label]
            for col_label in self.ColOrder or curr_row:
                if col_label in curr_row and new_f(curr_row[col_label]):
                    result.append((row_label, col_label))
        return result

    def getItemsIf(self, f, negate=False):
        """Returns list of items where f(self[row][col]) is True."""
        return self.getItems(self.getItemIndices(f, negate))

    def toLists(self, headers=False):
        """Return copy of self as list of lists, in order if specified.

        headers specifies whether to include the row and col headers (default:
        False).

        If the class data specifies RowOrder and ColOrder, can recapture the
        data in a new Dict2D object using self.fromLists() on the result of this
        method.

        Pads with self.Default if self.Pad is True; otherwise, raises exception
        on missing values.

        The result of toLists() can be passed to the array() function of
        numpy to generate a real array object for numerical calculations (if
        headers are False, the default).
        """
        col_order = self.ColOrder or self.colKeys()
        row_order = self.RowOrder or self.rowKeys()
        if self.Pad:
            default = self.Default
            result = []
            for r in row_order:
                result.append([self[r].get(c, default) for c in col_order])
        else:
            try:
                result = []
                for r in row_order:
                    result.append([self[r][c] for c in col_order])
            except KeyError:
                raise Dict2DSparseError, \
                "Unpadded Dict2D can't convert to list of lists if sparse."
        
        if headers:
            for header, row in zip(row_order, result):
                row.insert(0, header)
            result = [['-'] + list(col_order)] + result
        return result
  
    def copy(self):
        """Returns a new deep copy of the data in self (rows are new objects).
            
        NOTE: only copies the attributes in self._copied_attributes. Remember
        to override _copied_attributes in subclasses that need to copy
        additional data.
            
            usage:

                d = {'a':{'a':0}}
                m = Dict2D(d)
                c = m.copy()
                c['a']['a'] = 1
                c['a']['a'] != m['a']['a']
        """

        #create new copy of the data
        data = {}
        for key, row in self.items():
            data[key] = row.copy()
        
        #convert the result to the same class as self
        result = self.__class__(data)

        #copy any of self's attributes that aren't the same as the class values
        for attr in self._copied_attributes:
            curr_value = getattr(self, attr)
            if curr_value is not getattr(self.__class__, attr):
                setattr(result, attr, curr_value)

        return result
    
    def fill(self, val, rows=None, cols=None, set_orders=False):
        """Fills self[r][c] with val for r in rows and c in cols.

        If rows or cols is not specified, fills all the appropriate values
        that are already present in self.

        val: value to fill with. Required.

        rows: row indices to fill. If None, all rows are affected. Every row
        in rows will be created if necessary.

        cols: col indices to fill. If None, all cols are affected. Every col
        in cols will be created if necessary.

        set_orders: if True, sets self.RowOrder to rows and self.ColOrder
        to cols (default False). Otherwise, RowOrder and ColOrder are
        unaffected. 
        
        NOTE: RowOrder and ColOrder are _not_ used as defaults for rows and 
        cols, in order to make it convenient to fill only the elements that
        actually exist.
        """
        if set_orders:
            self.RowOrder = rows
            self.ColOrder = cols
        
        if rows is None:        #affect all rows
            for r in self.values():
                for c in (cols or r):   #if not cols, affect everything in row
                    r[c] = val
        else:                   #affect only specified rows
            constructor = self.RowConstructor   #might need to create new rows
            for r in rows:
                #create row if needed
                if r not in self:
                    self[r] = constructor()
                #bind current row
                curr_row = self[r]
                for c in (cols or curr_row):
                    curr_row[c] = val
            
    def setDiag(self, val):
        """Set the diagonal to val (required).
        
        Note: only affects keys that are rows (i.e. does not create rows for
        keys that are in the cols only). Use self.square() in this situation.
        """
        for k in self:
            self[k][k] = val
                
    def scale(self, f):
        """Applies f(x) to all elements of self."""
        for r in self:
            curr_row = self[r]
            for c, val in curr_row.items():
                curr_row[c] = f(val)

    def transpose(self):
        """Converts self in-place so self[r][c] -> self[c][r].
        
        Also swaps RowOrder and ColOrder.
        """
        t = {}
        for r, row in self.items():
            for c, val in row.items():
                if c not in t:
                    t[c] = {}
                t[c][r] = val
        self.clear()
        self.fromDicts(t)
        self.RowOrder, self.ColOrder = self.ColOrder, self.RowOrder
 
    def reflect(self, method=average):
        """Reflects items across diagonal, as given by RowOrder and ColOrder.

        Fails if RowOrder and ColOrder aren't equal, or if they're unspecified.
        Works fine on triangular matrices (i.e. no need to make the matrix
        square first) -- however, the diagonal won't be set if it wasn't
        present originally.
        """
        row_order = self.RowOrder
        col_order = self.ColOrder
        if (row_order is None) or (col_order is None):
            raise Dict2DError, \
            "Can't reflect a Dict2D without both RowOrder and ColOrder."

        if row_order != col_order:
            raise Dict2DError, \
            "Can only reflect Dict2D if RowOrder and ColOrder are the same."

        constructor = self.RowConstructor
        default = self.Default
        for row_index in range(len(row_order)):
            r = row_order[row_index]
            if r not in self:
                self[r] = constructor()
            curr_row = self[r]
            
            for col_index in range(row_index):
                c = row_order[col_index]
                if c not in self:
                    self[c] = constructor()
                curr_col = self[c]
                # Note that rows and cols need to be transposed, the way the 
                # functions are currently written. I think that changing the
                # functions would break existing code, and they make sense
                # as written.
                curr_row[c], curr_col[r] = \
                    method(curr_col.get(r, default), curr_row.get(c, default))

    def toDelimited(self, headers=True, item_delimiter='\t', \
        row_delimiter = '\n', formatter=str):
        """Printable string of items in self, optionally including headers.

        headers: whether or not to print headers (default True).

        delimiter: inserted between fields (default tab).

        formatter: applied to each element (default str). Note that the 
        formatter is also applied to the headers, so should be a function that
        can handle strings as well as whatever the elements of the matrix are
        (i.e. if it's a number formatter, it should detect strings and return
        them without formatting).

        If RowOrder or ColOrder is present, only prints the relevant rows and
        cols. Always pads missing values with str(self.Default).
        """
        lists = self.toLists(headers)
        return row_delimiter.join([item_delimiter.join(map(formatter, r)) \
            for r in lists])
