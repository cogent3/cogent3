#!/usr/bin/env python
"""Parsers for the AAIndex file format.

AAIndex can be downloaded at: http://www.genome.ad.jp/dbget/aaindex.html

There are two main files: AAIndex1 contains linear measures (one number per
amino acid) of amino acid properties, while AAIndex2 contains pairwise measures
(one number per pair of amino acids, e.g. distance or similarity matrices).
"""
import re
from cogent.parse.record_finder import DelimitedRecordFinder
from string import rstrip
from cogent.maths.matrix.distance import DistanceMatrix

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Greg Caporaso"
__email__ = "caporaso@colorado.edu"
__status__ = "Production"

class AAIndexParser(object):
    """ Abstract class for AAIndex file parsers
        This file is an abstract class for the parsers of the two AAIndex
        files.  The only real difference between the files is that AAIndex1
        has one additional field, labeled in here as Correlating.
    
    """

    def __init__(self):
        """ Initialize the object. """
        
    def __call__(self, infile):
        """ Parse AAIndex file into dict of AAIndex objects with ID as key

            infile = file to parse as file object or list of lines

            Usage:
                aa1p = AAIndex1Parser()
                aaIndex1Objects = aa1p('data/AAIndex1')

                aa2p = AAIndex2Parser()
                aaIndex2Objects = aa2p('data/AAIndex2')
        """
        
        result = {}

        # Break down the file into records delimited by '//' and then
        # parse each record into AAIndexRecord objects which will be stored
        # in a dict keyed by the records unique ID string
        AAIndexRecordFinder = DelimitedRecordFinder('//', constructor=rstrip)
        # parser is a generator of AAIndexRecords from file
        parser = AAIndexRecordFinder(infile)       

        for r in parser:
            new_record = self._parse_record(r)
            if new_record:
                yield new_record

    def _get_field(self, field_identifier, lines):
        """ Returns the field identified as a one line string
        """
        i = 0
        result = ''
        # Concatenate multi-line data with line_split
        line_split = ' '
        # Run through all lines in the current record
        while (i < len(lines)):
            # Check each line to see if it starts with the field
            # identifier we are looking for
            if (lines[i].startswith(field_identifier)):
                # If we find the line we are looking for, include it in
                # the result, unless it's a Data line.
                # Data entries are multi-line, and the first is information
                # that we are not interested in here.
                if (field_identifier != 'I'):
                    result += lines[i]
                    if field_identifier == 'M': result += 'BRK'
                    # Get rid of the line identifier and leading white space
                    result = result[2:]
                # Move to next line
                i += 1
                # and see if it's a continuation from the above line
                while (i < len(lines) and\
                     (lines[i].startswith(' ') or\
                     lines[i].startswith(field_identifier))):
                    # if continuation combine the lines while treating the
                    # spaces nicely, ie, multiple spaces -> one space
                    # this is mostly just important for the
                    # lines that are strings such as title
                    result = result.rstrip() + line_split + lines[i].lstrip()
                    i += 1
                break
            i += 1
        # return the field of interest   
        return result
        
class AAIndex1Parser(AAIndexParser):
    """ Parse AAIndex1 file & return it as dict of AAIndex1 objects"""

    def _parse_record(self, lines):
        """ Parse a single record and return it as a AAIndex1Record Object """
        # init all of the fields each time, this is so that
        # if fields are missing they don't get the value from the last
        # record
        id = None
        description = None
        LITDB = None
        authors = None
        title = None
        citations = None
        comments = None
        correlating = {}
        data = [None] * 20

        id = self._get_field('H', lines)
        description = self._get_field('D', lines)
        LITDB = self._get_field('R', lines)
        authors = self._get_field('A', lines)
        title = self._get_field('T', lines)
        citations = self._get_field('J', lines)
        comments = self._get_field('*', lines)
        correlating = self._parse_correlating(self._get_field('C', lines))
        data = self._parse_data(self._get_field('I', lines))

        return AAIndex1Record(id, description, LITDB, authors,\
                title, citations, comments, correlating, data)
                    

    def _parse_correlating(self, raw):
        """ Parse Correlating entries from the current record """
        keys = []
        values = []
        raw = raw.lstrip()
        # Split by white space
        data = re.split('\s*', raw)

        i=0
        while(i<len(data)):
            # If it's even it's a key
            if((i % 2) == 0):
                keys += [data[i]]
            # if it's not even it's a value
            else:
                # convert values to floats
                try:
                    values += [float(data[i])]
                except ValueError:
                    values += [data[i]]
            i += 1
        result = dict(zip(keys, values))
        return result

    def _parse_data(self, raw):
        """ Parse the data field from current record into a dict
        """
        # init for use in result
        keys = 'ARNDCQEGHILKMFPSTWYV'  
        values = []
        
        # get rid of leading white spaces, it makes../ the reg exp act weird
        raw = raw.lstrip()
        # split by any number/ types of white spaces
        data = re.split('\s*', raw)
        # convert the data to a float while checking for invlaid data,
        # specifically the string 'NA' is present sometimes instead of data
        for i in data:
            try:
                values += [float(i)]
            except ValueError:
                values += i

        result = dict(zip(keys, values))
        # return the dict
        return result


class AAIndex2Parser(AAIndexParser):
    """ Parse AAIndex2 file & return it as dict of AAIndex2 objects"""

    def _parse_record(self, lines):
        """ Parse a single record and return it as a AAIndex2Record Object """
        # Init attributes of each record each run through
       
        id = None
        description = None
        LITDB = None
        authors = None
        title = None
        citations = None
        comments = None
        rowscols = None
        data = []

        # Fill in the values
        id = self._get_field('H', lines)
        description = self._get_field('D', lines)
        LITDB = self._get_field('R', lines)
        authors = self._get_field('A', lines)
        title = self._get_field('T', lines)
        citations = self._get_field('J', lines)
        comments = self._get_field('*', lines)
        raw_data = self._get_field('M', lines)

        rowscols = self._parse_rowscols(raw_data[:raw_data.find('BRK')])
        try:
            data = self._parse_data(raw_data[raw_data.find('BRK')+3:],\
            rowscols[0], rowscols[1])
        except IndexError:
            return None

        return AAIndex2Record(id, description, LITDB, authors,\
            title, citations, comments, data)                       

    def _parse_data(self, raw, rows, cols):
        """ Parse the data field from current record into dict """
        # init result dict
        result = None
        # get rid of leading white spaces, it make the reg exp act weird
        raw = raw.lstrip()
        # split by any number/ types of white spaces
        data = re.split('\s*', raw)


        # If square matrix
        if len(data) == (len(rows)*len(cols)):
            result = dict.fromkeys(rows)
            i = 0
            for r in rows:
                new_row = dict.fromkeys(cols)
                for c in cols:
                    try:
                        new_row[c] = float(data[i])
                    except ValueError:
                        new_row[c] = data[i]
                    i+=1
                result[r] = new_row

        # else if LTM
        elif len(data) == (len(cols)+1) * len(rows)/2 :
            result = dict.fromkeys(rows)
            i = 0
            for r in rows:
                new_row = dict.fromkeys(cols)
                for c in cols:
                    if cols.find(c) <= rows.find(r):
                        try:
                            new_row[c] = float(data[i])
                        except ValueError:
                            new_row[c] = data[i]
                        i += 1
                result[r] = new_row                      
            
        return result

    def _parse_rowscols(self, raw):
        """ Returns two element list, 0: rows info, 1: cols info
        
            This parses the data out of the data description line
            for each record in AAIndex2 so we know what the data is that
            we are looking at.
        """
        p ='[rows|cols]\s=\s([^ \t\n\r\f\v,]*)'
        result = []
        result += re.findall(p, raw)
        return result


class AAIndexRecord(object):
    """ Abstract class, stores records from AAIndex files """

    def __init__(self, id,
                  description, LITDB_entry_num,
                  authors, title,
                  citation, comments, data):
        """ Stores data for individual AAIndex entires """

        self.ID = str(id)
        self.Description = str(description)
        self.LITDBEntryNum = str(LITDB_entry_num)
        self.Authors = str(authors)
        self.Title = str(title)
        self.Citation = str(citation)
        self.Comments = str(comments)
        self.Data = data

    def _toSquareDistanceMatrix(self, include_stops=False):
        """ Converts AAIndex Data to square distance matrix

            This abstract method must be overwritten for each subclass.
            The interface must be identical across subclasses, must
            take self and return new square matrix (for now).
        """
        pass


    def toDistanceMatrix(self, include_stops=False):
        """ Builds a DistanceMatrix object based on self """
        data = self._toSquareDistanceMatrix(include_stops=include_stops)

        # If there is missing or invalid data, data will be None
        # if that's the case return None for easy detection, otherwise
        # return a new DistanceMatrix object
        if data:
            return DistanceMatrix(data=data, info=self)

        return None

class AAIndex1Record(AAIndexRecord):
    """ Stores records from AAIndex1, inherits from AAIndexRecord """

    def __init__(self, id,
                  description, LITDB_entry_num,
                  authors, title,
                  citation, comments,
                  correlating, data):
        """ Stores data for individual AAIndex 1 entires """

        # Call init from super class
        AAIndexRecord.__init__(self, id,
                  description, LITDB_entry_num,
                  authors, title,
                  citation, comments, data)

        self.Correlating = correlating

    def _toSquareDistanceMatrix(self, include_stops=False):
        """ AAIndex1 data to square distance matrix

        """
        keys = self.Data.keys()
        if include_stops : keys += '*'

        # build result dict top layer, start empty
        result = {}
        for r in keys:
            new_row = {}
            for c in keys:
                if (r == '*' or c == '*'):
                    new_row[c] = None
                else:
                    # Build the ditance matrix by subtracting the
                    # value of each aminoacid and then taking the
                    # absolute value.  If the data can not be
                    # turned into a float, it's not a number, so the data
                    # is invalid. Return None for easy detection
                    try:
                        new_row[c] =\
                            abs(float(self.Data[r])
                             - float(self.Data[c]))
                    except ValueError:
                        return None
            result[r] = new_row

        return result


class AAIndex2Record(AAIndexRecord):
    """ Stores records from AAIndex2, inherits from AAIndexRecord  """
    def __init__(self, id,
                  description, LITDB_entry_num,
                  authors, title,
                  citation, comments, data):
        """ Stores data for individual AAIndex 2 entires """

        # Call init from super class
        AAIndexRecord.__init__(self, id,
                  description, LITDB_entry_num,
                  authors, title,
                  citation, comments, data)


    def _toSquareDistanceMatrix(self, include_stops=False):
        """ Returns data as a square matrix

            Note: This method is not currently functional,
            we are awaiting information on how to process data into
            a distance matrix

        """
        # create a new dict based on self.Data so we don't alter self.Data

        result = dict(self.Data)
        # Add in the new row of stop codon data
        if include_stops:
            stop_row = {}
            for i in result:
                stop_row.update({i:None})
            result.update({'*':stop_row})
            for i in result:
                result[i].update({'*':None})

        # Right now we are only dealing with square matrices
        return result

def AAIndexLookup(records):
    """ Build a dict of AAIndexObjects hashed by ID """    
    result = {}
    for r in records:
        result[r.ID] = r

    return result
        
def AAIndex1FromFiles(file):
    """ Taking a file or list of data return a dict of AAIndex1Objects """
    aap = AAIndex1Parser()
    return AAIndexLookup(aap(file))

def AAIndex2FromFiles(file):
    """ Taking a file or list of data return a dict of AAIndex2Objects """
    aap = AAIndex2Parser()
    return AAIndexLookup(aap(file))    

Woese_data = """//
H WOEC730101
D Polar requirement (Woese, 1973)
R PMID:4588588
A Woese, C.R.
T Evolution of genetic code
J Naturwiss. 60, 447-459 (1973)
C GRAR740102    0.960  HOPT810101    0.886  HOPA770101    0.876
  LEVM760101    0.872  PRAM900101    0.871  ROSM880101    0.844
  WOLS870101    0.841  KUHL950101    0.837  OOBM770103    0.835
  VINM940101    0.834  PARJ860101    0.821  FUKS010102    0.820
  FAUJ880110    0.812  OOBM770101    0.804  ROSM880102    0.801
  NADH010102   -0.800  CIDH920105   -0.800  MEIH800103   -0.802
  ISOY800102   -0.803  EISD860103   -0.803  ROSG850102   -0.804
  TANS770103   -0.806  RADA880101   -0.812  BIOV880102   -0.819
  WIMW960101   -0.821  NISK860101   -0.822  PONP800103   -0.823
  CIDH920104   -0.823  RADA880108   -0.825  BIOV880101   -0.829
  PONP800108   -0.831  SWER830101   -0.832  EISD860101   -0.838
  MAXF760102   -0.842  DESM900102   -0.847  FAUJ830101   -0.880
I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
     7.0     9.1    10.0    13.0     5.5     8.6    12.5     7.9     8.4     4.9
     4.9    10.1     5.3     5.0     6.6     7.5     6.6     5.3     5.7     5.6
//
"""

def getWoeseDistanceMatrix():
    """ Return the Woese Polar Requirement Distance Matrix """
    aaindexObjects = AAIndex1FromFiles(Woese_data.split('\n'))
    distance_matrices = {}
    for m in aaindexObjects:
        distance_matrices[m] = aaindexObjects[m].toDistanceMatrix()

    return distance_matrices['WOEC730101']    

