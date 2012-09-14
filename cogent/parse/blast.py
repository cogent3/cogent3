#!/usr/bin/env python
"""Parsers for blast, psi-blast and blat.
"""
from cogent.parse.record_finder import LabeledRecordFinder, \
    DelimitedRecordFinder, never_ignore
from cogent.parse.record import RecordError
from string import strip, upper

__author__ = "Micah Hamady"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Micah Hamady", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Micah Hamady"
__email__ = "hamady@colorado.edu"
__status__ = "Prototype"

def iter_finder(line):
    """Split record on rows that start with iteration label."""
    return line.startswith("# Iteration:")

def query_finder(line):
    """Split record on rows that start with query label."""
    return line.startswith("# Query:")

def iteration_set_finder(line):
    """Split record on rows that begin a new iteration."""
    return line.startswith("# Iteration: 1")

def _is_junk(line, t_strs):
    """Ignore empty line, line with blast info, or whitespace line"""
    # empty or white space
    if not line or not line.strip():
        return True
    # blast info line
    for t_str in t_strs:
        if line.startswith("# %s" % t_str):
            return True
    return False

def is_blast_junk(line):
    """Ignore empty line or lines with blast info"""
    return _is_junk(line, ("BLAST","TBLAS"))

def is_blat_junk(line):
    """Ignore empty line or lines with blat info"""
    return _is_junk(line, ("BLAT",))

label_constructors = {'ITERATION': int} #add other label constructors here

def make_label(line):
    """Make key, value for colon-delimited comment lines.
    
    WARNING: Only maps the data type if the key is in label_constructors above.
    """
    if not line.startswith("#"):
        raise ValueError, "Labels must start with a # symbol."

    if line.find(":") == -1:
        raise ValueError, "Labels must contain a : symbol."

    key, value = map(strip, line[1:].split(":", 1))
    key = key.upper()
    if key in label_constructors:
        value = label_constructors[key](value)
    return key, value

BlatFinder = LabeledRecordFinder(query_finder, constructor=strip, \
    ignore=is_blat_junk)

BlastFinder = LabeledRecordFinder(query_finder, constructor=strip, \
    ignore=is_blast_junk)

PsiBlastFinder = LabeledRecordFinder(iter_finder, constructor=strip, \
    ignore=is_blast_junk)

PsiBlastQueryFinder = LabeledRecordFinder(iteration_set_finder, \
    constructor=strip, ignore=is_blast_junk)

def GenericBlastParser9(lines, finder, make_col_headers=False):
    """Yields successive records from lines (props, data list) 
        
    Infile must in blast9 format
       
    finder: labeled record finder function
    
    make_col_header: adds column headers (from fields entry) as first
    row in data output

    props is a dict of {UPPERCASE_KEY:value}.
    data_list is a list of list of strings, optionally with header first.
    """
    for rec in finder(lines):
        props = {}
        data = []
        for line in rec:
            if line.startswith("#"):
                label, value = make_label(line)
                props[label] = value

                # check if need to insert column headers
                if make_col_headers and label == "FIELDS":
                    data.insert(0, map(upper, map(strip,value.split(","))))

            else:
                data.append(map(strip, line.split("\t")))
        yield props, data

def TableToValues(table, constructors=None, header=None):
    """Converts table to values according to constructors.
    
    Returns (table, header). 
    Use dict([(val, i) for i, val in enumerate(header)]) to get back
    a dict mapping the fields to indices in each row.
    """
    if header is None:  #assume first row of table
        header = table[0]
        table = table[1:]
    c_list = [constructors.get(k, str) for k in header]
    return [[c(val) for c, val in zip(c_list, row)] for row in table], header

psiblast_constructors={'% identity':float, 'alignment length':int, \
    'mismatches':int, 'gap openings':int, 'q. start':int, 'q. end':int, \
    's. start':int, 's. end':int, 'e-value':float, 'bit score':float}
#make case-insensitive
for key, val in psiblast_constructors.items():
    psiblast_constructors[key.upper()] = val

def PsiBlastTableParser(table):
    return TableToValues(table, psiblast_constructors)

def MinimalBlastParser9(lines, include_column_names=False):
    """Yields succesive records from lines (props, data list).

    lines must be BLAST output format.
    """
    return GenericBlastParser9(lines, BlastFinder, include_column_names)

def MinimalPsiBlastParser9(lines, include_column_names=False):
    """Yields successive records from lines (props, data list) 
        
        lines must be of psi-blast output format
    """
    return GenericBlastParser9(lines, PsiBlastFinder, include_column_names)

def MinimalBlatParser9(lines, include_column_names=True):
    """Yields successive records from lines (props, data list) 
        
       lines must be of blat output (blast9) format 
    """
    return GenericBlastParser9(lines, BlatFinder, include_column_names)

def PsiBlastParser9(lines):
    """Returns fully parsed PSI-BLAST result.

    result['query'] gives all the results for specified query sequence.
    result['query'][i] gives result for iteration i (offset by 1: zero-based)
    if x = result['query']['iteration']:
        x[0]['e-value'] gives the e-value of the first result.

    WARNING: designed for ease of use, not efficiency!"""
    result = {}
    for query in PsiBlastQueryFinder(lines):
        first_query = True  #if it's the first, need to make the entry
        for properties, record in MinimalPsiBlastParser9(query, True):
            if first_query:
                curr_resultset = []
                result[properties['QUERY'].split()[0]] = curr_resultset
                first_query = False
            table, header = PsiBlastTableParser(record)
            curr_resultset.append([dict(zip(header, row)) for row in table])
    return result

def get_blast_ids(props, data, filter_identity, threshold, keep_values):
    """
    Extract ids from blast output
    """
    fields = map(strip, props["FIELDS"].upper().split(","))

    # get column index of protein ids we want
    p_ix = fields.index("SUBJECT ID")
    # get column index to screen by
    if filter_identity:
        e_ix = fields.index("% IDENTITY")
    else:
        e_ix = fields.index("E-VALUE")
    # no filter, returh all
    if not threshold:
        if keep_values:
           return [(x[p_ix],x[e_ix]) for x in data]
        else:
            return [x[p_ix] for x in data]
    else:
        # will raise exception if invalid threshold passed 
        max_val = float(threshold)
        
        #figure out what we're keeping
        def ok_val(val):
            if threshold:
                return (val <= max_val)
            return (val >= max_val)
        if keep_values:
            return [(x[p_ix],x[e_ix]) for x in data if ok_val(float(x[e_ix]))]
        else:
            return [x[p_ix] for x in data if ok_val(float(x[e_ix]))]



def AllProteinIds9(lines, filter_identity=True, threshold=None, \
        keep_below_threshold=True, output_parser=MinimalPsiBlastParser9,
        keep_values=False):
    """Helper to extract just protein ids from each blast search 

    lines: output file in output format #9. 
    filter_identity: when True, use % identity to filter, else use e-value
    threshold: when None, all results are returned. When not None, used 
        as a threshold to filter results.
    keep_below_threshold: when True, keeps any rows below given threshold, else
        keep any rows above threshold
    output_parser: minimal output parser to use (e.g. minimalpsiblast)
    keep_values: if True, returns tuples of (id, value) rather than just ids.
    
    Note that you can feed it successive output from PsiBlastQueryFinder if
    you have a PSI-BLAST file with multiple input queries.

    Subject ids are stable relative to original order.
    """

    mpbp = output_parser(lines)

    # get last record. 
    props = data = None
    out_ids = {}
    out_ct = 1
    for rec in mpbp:
        props, data = rec
        out_ids[out_ct] = get_blast_ids(props, data, filter_identity, 
                                        threshold, keep_values)
        out_ct += 1
    return out_ids 

def LastProteinIds9(lines, filter_identity=True, threshold=None, \
        keep_below_threshold=True, output_parser=MinimalPsiBlastParser9,
        keep_values=False):
    """Helper to extract just protein ids from last psi-blast iteration.

    lines: output file in output format #9. 
    filter_identity: when True, use % identity to filter, else use e-value
    threshold: when None, all results are returned. When not None, used 
        as a threshold to filter results.
    keep_below_threshold: when True, keeps any rows below given threshold, else
        keep any rows above threshold
    output_parser: minimal output parser to use (e.g. minimalpsiblast)
    keep_values: if True, returns tuples of (id, value) rather than just ids.
    
    Note that you can feed it successive output from PsiBlastQueryFinder if
    you have a PSI-BLAST file with multiple input queries.

    Subject ids are stable relative to original order.
    """

    mpbp = output_parser(lines)
    # get last record. 
    props = data = None
    for rec in mpbp:
        props, data = rec
    if not (props and data):
        return []
    return get_blast_ids(props, data, filter_identity, threshold, keep_values)

def QMEBlast9(lines):
    """Returns query, match and e-value for each line in Blast-9 output.
  
    WARNING: Allows duplicates in result.
  
    WARNING: If you use this on PSI-BLAST output, will not check that you're
    only getting stuff from the last iteration but will give you everything.
    The advantage is that you keep stuff that drops out of the profile. The
    disadvantage is that you keep stuff that drops out of the profile...
    """
    result = []
    for line in lines:
        if line.startswith('#'):
            continue
        try:
            fields = line.split('\t')
            result.append((fields[0], fields[1], float(fields[-2])))
        except (TypeError, ValueError, IndexError):
            pass
    return result

def QMEPsiBlast9(lines):
    """Returns successive query, match, e-value from lines of Psi-Blast run.

    Assumes tabular output. Uses last iteration from each query.

    WARNING: Allows duplicates in result
    """
    result = []
    for query in PsiBlastQueryFinder(lines):
        for iteration in PsiBlastFinder(query):
            pass
        result.extend(QMEBlast9(iteration))
    return result

class BlastResult(dict):
    """Adds convenience methods to BLAST result dict.
    
    {Query:[[{Field:Value}]]}

    Nesting is:
    query: key/value
    iteration: list
    hit: list
    field: key/value

    For BLAST, there is always exactly one iteration, but PSIBLAST can have
    multiple. Keep interface the same.

    Question: should it be able to construct itself from the result string?
    """
    # FIELD NAMES

    ITERATION = 'ITERATION'
    QUERY_ID = 'QUERY ID'
    SUBJECT_ID = 'SUBJECT ID'
    PERCENT_IDENTITY = '% IDENTITY'
    ALIGNMENT_LENGTH = 'ALIGNMENT LENGTH'
    MISMATCHES = 'MISMATCHES'
    GAP_OPENINGS = 'GAP OPENINGS'
    QUERY_START = 'Q. START'
    QUERY_END = 'Q. END'
    SUBJECT_START = 'S. START'
    SUBJECT_END = 'S. END'
    E_VALUE = 'E-VALUE'
    BIT_SCORE = 'BIT SCORE'

    #standard comparison for each field, e.g.  
    #want long matches, small e-values
    _lt = lambda x, y: cmp(x, y) == -1
    _le = lambda x, y: cmp(x, y) <= 0 
    _gt = lambda x, y: cmp(x, y) == 1 
    _ge = lambda x, y: cmp(x, y) >= 0 
    _eq = lambda x, y: cmp(x, y) == 0 

    FieldComparisonOperators = {
               PERCENT_IDENTITY:(_gt, float),
               ALIGNMENT_LENGTH:(_gt, int),
               MISMATCHES:(_lt, int),
               E_VALUE:(_lt, float),
               BIT_SCORE:(_gt, float)
                }   

    # set up valid blast keys
    HitKeys = set([ ITERATION,
                        QUERY_ID,
                        SUBJECT_ID,
                        PERCENT_IDENTITY,
                        ALIGNMENT_LENGTH,
                        MISMATCHES,
                        GAP_OPENINGS,
                        QUERY_START,
                        QUERY_END,
                        SUBJECT_START,
                        SUBJECT_END,
                        E_VALUE,
                        BIT_SCORE ])

  
    def __init__(self, data, psiblast=False):
        """
        Init using blast results

        data: blast output from the m = 9 output option
        psiblast: if True, will expect psiblast output, else expects 
            blast output

        """
        parser = MinimalBlastParser9
        if psiblast:
            parser = MinimalPsiBlastParser9

        mp = parser(data, True)
        
        
        for props, rec_data in mp:
         
            iteration = 1
            if self.ITERATION in props:
                iteration = int(props[self.ITERATION])

            hits = []
            # check if found any hits
            if len(rec_data) > 1:
                for h in rec_data[1:]:
                    hits.append(dict(zip(rec_data[0], h)))
            else:
                hits.append(dict(zip(rec_data[0], ['' for x in rec_data[0]])))
            
            # get blast version of query id
            query_id = hits[0][self.QUERY_ID]

            if query_id not in self: 
                self[query_id] = [] 
            self[query_id].append(hits)
        
    def iterHitsByQuery(self, iteration=-1):
        """Iterates over set of hits, returning list of hits for each query"""
        for query_id in self:
            yield query_id, self[query_id][iteration] 

    def iterHitsByTarget(self, iteration=-1):
        """Iterates over set of hits, returning list of hits for each target"""
        raise NotImplementedError

    def iterAllHits(self, iteration=-1):
        """Iterates over all hits, one at a time"""
        raise NotImplementedError
        
    def filterByField(self, field='E-value', threshold=0.001):
        """Returns a copy of self containing hits where field better than threshold.
        Uses FieldComparisonOperators to figure out which direction to compare.
        """
        raise NotImplementedError
        
    def filterByFunc(self, f):
        """Returns copy of self containing hits where f(entry) is True."""
        raise NotImplementedError
        
    def bestHitsByQuery(self, iteration=-1,  n=1, field='BIT SCORE', return_self=False):
        """Iterates over all queries and returns best hit for each 
        
        return_self: if False, will not return best hit as itself.

        Uses FieldComparisonOperators to figure out which direction to compare.
        """

        # check that given valid comparison field
        if field not in self.FieldComparisonOperators:
            raise ValueError, "Invalid field: %s. You must specify one of: %s" \
                              % (field, str(self.FieldComparisonOperators))
        cmp_fun, cast_fun = self.FieldComparisonOperators[field]

        # enumerate hits
        for q, hits in self.iterHitsByQuery(iteration=iteration):
            best_hits = [] 
            for hit in hits:
                # check if want to skip self hit 
                if not return_self:
                    if hit[self.SUBJECT_ID] == q:
                        continue
                # check if better hit than ones we have
                if len(best_hits) < n:
                    best_hits.append(hit)
                else:
                    for ix, best_hit in enumerate(best_hits):  
                        new_val = cast_fun(hit[field])
                        old_val = cast_fun(best_hit[field])
                        if cmp_fun(new_val, old_val):
                            best_hits[ix] = hit
                            continue
            yield q, best_hits

    def filterByIteration(self, iteration=-1):
        """Returns copy of self containing only specified iteration.

        Negative indices count backwards."""
    
    #raise error if both field and f passed, uses same dict as filterByField

fastacmd_taxonomy_splitter = DelimitedRecordFinder(delimiter='', \
    ignore=never_ignore)
fasta_field_map = { 'NCBI sequence id':'seq_id',
        'NCBI taxonomy id':'tax_id',
        'Common name':'common_name',
        'Scientific name':'scientific_name'}

def FastacmdTaxonomyParser(lines):
    """Yields successive records from the results of fastacmd -T.

    Format is four lines separated by newline:
    NCBI sequence
    NCBI taxonomy
    Common name
    Scientific name

    Result is dict with keys by seq_id, tax_id, common_name, scientific_name.
    """
    for group in fastacmd_taxonomy_splitter(lines):
        result = {}
        for line in group:
            try:
                header, data = line.split(':', 1)
                result[fasta_field_map[header]] = data.strip()
            except (TypeError, ValueError, KeyError):
                continue
        yield result


