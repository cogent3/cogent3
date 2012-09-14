#!/usr/bin/env python
"""
This module provides a parser for RNAView output.

Authors: 
Greg Caporaso (Gregory.Caporaso@UCHSC.edu)
Sandra Smit (Sandra.Smit@colorado.edu)

RNAView reports annotated base pairs found in a PDB file.
http://ndbserver.rutgers.edu/services/download/index.html

Constants:
RNAVIEW_ACCEPTED -- residues accepted by RNAView
WC_PAIRS -- dict of Watson-Crick pairs
WOBBLE_PAIRS -- dict of Wobble pairs

Important objects:
Base -- holds an RNA residue
BasePair -- holds a base pair
BasePairs -- list of BasePair objects
BaseMultiplet -- holds 3 or more residues
BaseMultiplets -- list of BaseMultiplet objects

Important functions:
RnaviewParser -- parses RNAView output
several selection functions to select a subset of base pairs, for example:
is_canonical -- selects all the pairs annotated as canonical by RNAView

Exceptions:
RnaViewObjectError -- raised when an RNAView object detects a problem
RnaViewParseError -- raised when one of the parser fucntions finds a problem
"""

__author__ = "Greg Caporaso and Sandra Smit"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Greg Caporaso", "Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

RNAVIEW_ACCEPTED = dict.fromkeys('AGUCTIPaguct')
WC_PAIRS = dict.fromkeys(['GC','CG','AU','UA'])
WOBBLE_PAIRS = dict.fromkeys(['UG','GU'])

class RnaViewObjectError(ValueError):
    pass

class RnaViewParseError(SyntaxError):
    pass

# =============================================================================
# RNAVIEW OBJECTS
# =============================================================================

class Base(object):
    """Describes a residue"""

    def __init__(self, ChainId, ResId, ResName, RnaViewSeqPos=None):
        """Initialize Base object.

        ChainId -- Chain identifier, string. If no chain is specified, 
            the ChainId is set to ' ' (for compatibility with PyMol)
        ResId -- Residue identifier, string. Cannot be None
        ResName -- Name of the residue, string. Accepted are AGUCTIPaguct.
            Cannot be None.
        RnaViewSeqPos -- string, position that RnaView assigns to this base.
        """
        if not ResId or not ResName:
            raise RnaViewObjectError(\
            "Found missing attribute: ResId=%s, ResName=%s"\
                %(ResId, ResName))
        self.ChainId = ChainId
        self.ResId = ResId
        self.ResName = ResName
        self.RnaViewSeqPos = RnaViewSeqPos

    def __str__(self):
        """Return string representation of Base object."""
        return ' '.join(map(str,[self.ChainId, self.ResId, self.ResName]))

    def __eq__(self, other):
        """Overwrite the == operator."""
        if self.ChainId != other.ChainId:
            return False
        if self.ResId != other.ResId:
            return False
        if self.ResName != other.ResName:
            return False
        if self.RnaViewSeqPos != other.RnaViewSeqPos:
            return False
        return True

    def __ne__(self,other):
        """Overwrite the != operator."""
        return not self == other

class BasePair(object):
    """Object for storing paired RNA bases.
    
    This object is made to store two Base objects which are involved in a
    base pairing interaction with each other.
    """

    def __init__(self, Up, Down, Edges=None, Orientation=None,\
        Conformation=None, Saenger=None):
        """Initialize the object. 
            
        Up -- the upstream Base object
        Down -- the downstream Base object
        Edges -- string, the pairing edges of the base pair or 'stacked'
        Orientation -- string, the orientation of the bases: 'cis' or 'tran'
        Conformation -- string, 'syn' or 'syn syn'
        Saenger -- string, the Saenger classification of the base pair, either
            roman numeral, 'n/a', or something starting with '!'
        """
        self.Up = Up
        self.Down = Down
        self.Edges = Edges
        self.Orientation = Orientation
        self.Conformation = Conformation
        self.Saenger = Saenger

    def __str__(self):
        """Return string representation of BasePair object."""
        return "Bases: %s -- %s; Annotation: %s -- %s -- %s -- %s;"\
            %(self.Up, self.Down, self.Edges, self.Orientation,\
            self.Conformation,self.Saenger)

    def __eq__(self,other):
        """Overwrite the == operator."""
        if self.Up != other.Up:
            return False
        if self.Down != other.Down:
            return False
        if self.Edges != other.Edges:
            return False
        if self.Orientation != other.Orientation:
            return False
        if self.Conformation != other.Conformation:
            return False
        if self.Saenger != other.Saenger:
            return False

        return True

    def __ne__(self,other):
        """Overwrite the != operator."""
        return not self == other

    def isWC(self):
        """Return True if base pair is a Watson-Crick pair.

        WARNING: this method returns True for GC and AU pairs independent
        of the rest of the annotation. It does not check the specific edges
        in the base pair etc.
        """
        # The complicated looking one-liner just makes an upper-case string
        # out of the two Base identities and checks to see if it exists in
        # WC_PAIRS
        return ''.join([self.Up.ResName, self.Down.ResName]).upper()\
            in WC_PAIRS

    def isWobble(self):
        """Return True if base pair is a wobble pair.

        WARNING: this method returns True for GU pairs independent
        of the rest of the annotation. It does not check the specific edges
        in the base pair etc.
        """
        # The complicated looking one-liner just makes an upper-case string
        # out of the two Base identities and checks to see if it exists in
        # WOBBLE_PAIRS
        return ''.join([self.Up.ResName, self.Down.ResName]).upper()\
            in WOBBLE_PAIRS


class BasePairs(list):
    """A list of BasePair objects."""

    def __init__(self, base_pairs=None):
        """Initialize BasePairs object.

        base_pairs -- list or tuple of BasePair objects.
        """
        if base_pairs is None:
            base_pairs = []
        try:
            self[:] = list(base_pairs)
        except TypeError:
            raise RnaViewObjectError(\
                'base_pairs must be convertable to a list')

    def __str__(self):
        """Return string representation of BasePairs object."""
        header = [\
        "===================================================================",
        "Bases: Up -- Down; Annotation: Edges -- Orient. -- Conf. -- Saenger",
        "==================================================================="]

        return '\n'.join(header + [str(i) for i in self])

    def select(self, selection_function, rest=False):
        """Return a new BasePairs object of pairs that pass the selection.

        selection_function -- should be function that works on a BasePair
            and returns True or False.
        rest -- boolean, if True, in addition to the selection, all the pairs
            that didn't make it into the selection are returned in a list.
            The return value is thus a tuple in this case.
        """
        sel = []
        if rest:
            not_sel = []
        for bp in self:
            if selection_function(bp):
                sel.append(bp)
            else:
                if rest:
                    not_sel.append(bp)
                else:
                    continue
        if rest:
            return BasePairs(sel), BasePairs(not_sel)
        else:
            return BasePairs(sel)

    def _get_present_chains(self):
        """Return the chains present in the BasePairs object (in ORDER)."""
        up_chains = set([i.Up.ChainId for i in self])
        down_chains = set([i.Down.ChainId for i in self])
        result = list(up_chains | down_chains)
        result.sort()
        return result

    PresentChains = property(_get_present_chains)

    def cliques(self):
        """Yield all cliques in the base pairs as new BasePairs objects.

        A clique is a set of interacting chains. All base pairs involved
        in a clique are returned together in a new BasePairs object.
        """
        if len(self.PresentChains) == 1:
            yield self
        else:
            ia = {} #interactions
            for i in self:
                chain_up = i.Up.ChainId
                chain_down = i.Down.ChainId
                if not ia:
                    bp_list = [i]
                    ia[chain_up] = bp_list
                    ia[chain_down] = bp_list
                elif chain_up in ia and not chain_down in ia:
                    ia[chain_up].append(i)
                    ia[chain_down] = ia[chain_up]
                elif chain_down in ia and not chain_up in ia:
                    ia[chain_down].append(i)
                    ia[chain_up] = ia[chain_down]
                elif chain_up not in ia and chain_down not in ia:
                    bp_list = [i]
                    ia[chain_up] = bp_list
                    ia[chain_down] = bp_list
                elif chain_up in ia and chain_down in ia and\
                    ia[chain_up] is ia[chain_down]:
                    ia[chain_up].append(i)
                elif chain_up in ia and chain_down in ia and not\
                    ia[chain_up] is ia[chain_down]:
                    ia[chain_up].append(i)
                    ia[chain_up].extend(ia[chain_down])
                    for k,v in ia.items():
                        if k != chain_down and v is ia[chain_down]:
                            ia[k] = ia[chain_up]
                    ia[chain_down] = ia[chain_up]
                else:
                    continue
            
            uniques = []
            for k,v in ia.items():
                if v not in uniques:
                    uniques.append(v)
            for bps in uniques:
                yield BasePairs(bps)

    def hasConflicts(self, return_conflict=False):
        """Return True if there is a conflict in the list of base pairs.

        return_conflict -- if True, return value is a tuple of a boolean and
            the string representation of the base involved in the conflict.

        A conflict occurs when according to RnaView a base is involved in
        more than one interaction. For example, base A pairs with base B,
        and base A pairs with base C.
        """
        seen = {}
        for bp in self:
            up = bp.Up
            down = bp.Down
            if (up.ChainId, up.ResId, up.ResName) not in seen:
                seen[(up.ChainId, up.ResId, up.ResName)] = None
            else:
                if return_conflict:
                    return True, str(up)
                else:
                    return True
            if (down.ChainId, down.ResId, down.ResName) not in seen:
                seen[(down.ChainId, down.ResId, down.ResName)] = None
            else:
                if return_conflict:
                    return True, str(down)
                else:
                    return True
        if return_conflict:
            return False, None
        else:
            return False         


class BaseMultiplet(list):
    """Hold a base multiplet (3 or more residues) found by RnaView."""
    
    def __init__(self, bases=None, NumberOfBases=None):
        """Initialize a BaseMultiplet object.

        bases -- a list or tuple of Base objects.
        NumberOfBases -- int, the number of bases in this multiplet.
        """
        if bases is None:
            bases = []
        try:
            self[:] = list(bases)
        except TypeError:
            raise RnaViewObjectError(\
                'bases must be convertable to a list')

        self.NumberOfBases = NumberOfBases
                
    def __str__(self):
        """Return string representation of a BaseMultiplet object."""
        return ' -- '.join([str(i) for i in self]) + ';'
       
       
class BaseMultiplets(list):
    """A list of BaseMultiplet objects."""

    def __init__(self, multiplets=None):
        """Initialize a BaseMultiplets object.

        multiplets -- list or tuple of BaseMultiplet objects.
        """
        if multiplets is None:
            multiplets = []
        try:
            self[:] = list(multiplets)
        except TypeError:
            raise RnaViewObjectError(\
                'multiplets must be convertable to a list')

    def __str__(self):
        """Return string representation of BaseMultiplets object."""
        return '\n'.join([str(i) for i in self])

class PairCounts(dict):
    """A dictionary of base pair counts."""
    pass 

# =============================================================================
# SELECTION FUNCTIONS
# =============================================================================

def in_chain(chains):
    """Return selection function that returns True when Up and Down in chains.

    chains -- list of chain IDs, e.g. "AB" or ['A','0']
    """
    def apply_to(bp):
        """Return True if Up and Down both in chains.
        
        bp -- BasePair object
        """
        if bp.Up.ChainId in chains and bp.Down.ChainId in chains:
            return True
        return False
    return apply_to

def is_canonical(bp):
    """Return True if base pair is standard canonical GC, AU, or wobble GU pair.

    bp -- BasePair object

    This functions only looks at the annotation, it doesn't check the base
    identity.
    """
    if bp.Edges == '+/+': #GC pair
        return True
    if bp.Edges == '-/-': #AU pair
        return True
    if bp.Edges == 'W/W' and bp.Orientation == 'cis' and\
        bp.Saenger == 'XXVIII': #GU pair
        return True
    return False

def is_not_canonical(bp):
    """Return True if base pair is not canonical.

    bp -- BasePair object

    See is_canoncical for a definition.
    """
    return not is_canonical(bp)

def is_stacked(bp):
    """Return True if base pair is stacked.
    
    bp -- BasePair object
    """
    if bp.Edges == 'stacked':
        return True
    return False

def is_not_stacked(bp):
    """Return True if bp is not stacked.
    
    bp -- BasePair object
    """
    if bp.Edges == 'stacked':
        return False
    return True

def is_tertiary(bp):
    """Return True if bp is a tertiary interaction.
    
    bp -- BasePair object
    """
    if bp.Saenger and bp.Saenger.startswith('!'):
        return True
    return False

def is_not_stacked_or_tertiary(bp):
    """Return True if the base pair is not stacked or a tertiary interaction.

    bp -- BasePair object
    """
    if bp.Edges == 'stacked':
        return False
    if bp.Saenger.startswith('!'):
        return False
    return True

def is_tertiary_base_base(bp):
    """Return True if base pair is a tertiary interaction between two bases.

    bp -- BasePair object

    Other tertiary interactions can for example be between bases and sugars.
    """
    if bp.Saenger and bp.Saenger.startswith('!1H(b_b)'):
        return True
    return False

#==============================================================================
# RNAVIEW PARSER
#==============================================================================

def is_roman_numeral(s):
    """Return True if s is a roman numeral.
    
    s -- string
    """
    if not s:
        return False
    # there is a comma in the alphabet, because some categories are
    # combined, split by comma
    alphabet = dict.fromkeys("IVXDCM,")
    for i in s:
        if i not in alphabet:
            return False
    return True

def is_edge(s):
    """Return True if s is a description of the base pair edges.
    
    s -- string

    X is added to the list of accepted symbols. X might be used
    for modified bases.
    """
    chars = dict.fromkeys("WHsS+-?.X")
    try:
        first, second = s.split('/')
    except ValueError:
        return False
    if first not in chars or second not in chars:
        return False
    return True


def is_orientation(s):
    """Return True if s is a description of the base pair orientation.
    
    s -- string
    """
    if s == 'cis' or s == 'tran':
        return True
    return False

def parse_annotation(data):
    """Parse the annotation of a base pair, returns tuple of 4 elements.
    
    data -- a single string containing the annotation of a single base pair.
    
    For example: '-/- cis         XX' or 'stacked'
    Edges: could be 'stacked' or two chars separated by '/', e.g. 'W/H'
    Orientation: could be 'cis' or 'tran'
    Conformation: could be 'syn' or 'syn syn'
    Saenger: could be some roman numeral (or two separated by a comma) or
        'n/a' or anything starting with an exclamation mark.
    """
    edges = None
    orientation = None
    conformation = None
    saenger = None
    
    for field in data:
        if field == 'stacked' or is_edge(field):
            edges = field
        elif is_orientation(field):
            orientation = field
        elif is_roman_numeral(field) or field == 'n/a' or\
            field.startswith('!'):
            saenger = field
        elif field == 'syn':
            if not conformation:
                conformation = field
            else:
                conformation += ' '+field
        else:
            raise RnaViewParseError("Unknown annotation field: %s"%(field))

    return edges, orientation, conformation, saenger

def parse_filename(lines):
    """Return PDB filename from the rnaview output file.
    
    lines -- list of strings or anything that behaves like a list of lines
    """
    if len(lines) != 1:
        raise RnaViewParseError(\
        "Parse filename: expected only one line to parse")
    return lines[0].split(':')[1].strip()

def parse_uncommon_residues(lines):
    """Return dictionary of {(chain_id, res_id, res_name): res_assigned}.

    lines -- list of strings or anything that behaves like a list of lines
    
    Parses the uncommon residue lines from an rnaview file.
    The chain_id could be empty (parsed as ' '), all other field should 
    have some value. If any of these is missing, an error will be raised.
    """
    result = {}
    for line in lines:
        res_name = line[17:20].strip()
        res_id = line[20:26].strip()
        chain_id = line[36]
        res_assigned = line.split(':')[1].strip()
        if not res_name or not res_id or not res_assigned:
            raise RnaViewParseError(\
            "Found missing field in uncommon residue line.\n"+\
            "res_id: %s, res_name: %s, res_assigned: %s"\
            %(res_id, res_name, res_assigned))
        result[(chain_id, res_id, res_name)] = res_assigned
    return result

def parse_base_pairs(lines):
    """Return BasePairs object. Parse BASE_PAIR section of rnaview output.

    lines -- list of strings or anything that behaves like a list of lines
   
    An empty chain_id is parsed as ' ' for compatibility with PyMol.
    """
    pairs = []
    for line in lines:
        data = line.split()
        seq_pos_up, seq_pos_down = data[0].strip(", ").split("_")
        up_chain_id = data[1].rstrip(':') or ' '
        down_chain_id = data[5].rstrip(':') or ' '
        up_res_id = data[2].strip()
        down_res_id = data[4].strip()
        up_res_name, down_res_name = data[3].split('-')
        
        #verify found residues are standard in Rnaview
        if up_res_name not in RNAVIEW_ACCEPTED or\
                down_res_name not in RNAVIEW_ACCEPTED:
            raise RnaViewParseError(\
                "Base found that is not generally accepted by Rnaview:"+\
                " %s or %s"%(up_res_name, down_res_name))
        
        #Build upstream and downstream Base objects
        up = Base(up_chain_id, up_res_id, up_res_name, RnaViewSeqPos=seq_pos_up)
        down = Base(down_chain_id, down_res_id, down_res_name,\
            RnaViewSeqPos=seq_pos_down)
        
        e, o, c, s = parse_annotation(data[6:])
        
        pairs.append(BasePair(up, down, Edges=e, Orientation=o,\
            Conformation=c, Saenger=s))
    return BasePairs(pairs)

def parse_base_multiplets(lines):
    """Return BaseMultiplets object. Parses BASE_MULTIPLETS section.

    lines -- list of strings or anything that behaves like a list of lines

    An empty chain_id is parsed as ' ' for compatibility with PyMol.
    """
    residues = []
    for line in lines:
        multi = []
        pos_info, rest = line.strip().split('|',1)
        rnaview_pos = [i.strip() for i in pos_info.split('_') if i != '']
        num_info, bases_info = rest.strip().split(']',1)
        num_bases = int(num_info.split()[1])
        bases = bases_info.split('+')
        if len(rnaview_pos) != len(bases):
            raise RnaViewParseError(\
                "Number of bases (%s) doesn't match number of positions (%s)"\
                %(len(bases),len(rnaview_pos)))
        if num_bases != len(bases):
            raise RnaViewParseError(\
                "Reported number of bases (%s) "%(num_bases)+\
                "doesn't match number of found bases (%s)"%(len(bases)))
        for rv_pos, b in zip(rnaview_pos, bases):
            data = b.split()
            chain_id = data[0].rstrip(':') or ' '
            res_id = data[1].strip()
            res_name = data[2].strip()
            
            #verify found residue is standard in Rnaview
            if res_name not in RNAVIEW_ACCEPTED:
                raise RnaViewParseError(\
                "Base found that is not generally accepted by RNAVIEW: %s"\
                    %(str(res_name)))
                        
            multi.append(Base(chain_id, res_id, res_name, RnaViewSeqPos=rv_pos))
        residues.append(BaseMultiplet(multi, NumberOfBases=num_bases))
    return BaseMultiplets(residues)

def parse_number_of_pairs(lines):
    """Return dict with the number of bases and number of base pairs.
    
    lines -- list of strings or anything that behaves like a list of lines

    The two keys in the dictionary are:
    NUM_PAIRS: contains the number of base pairs (note this number only
        includes normal base pairs, not stacked pairs, or tertiary interactions
    NUM_BASES: contains the number of RNA/DNA bases
    """
    if len(lines) != 1:
        raise RnaViewParseError(\
            "parse_number_of_parse should get only a single line")
    result = {}
    parts = lines[0].split('=')[1].split()
    if len(parts) != 4:
        raise RnaViewParseError("Can't parse 'total base pairs' line")
    result['NUM_PAIRS'] = int(parts[0])
    result['NUM_BASES'] = int(parts[2])
    return result

def parse_pair_counts(lines):
    """Parse summary of base pairs at the end of the rnaview output file.

    lines -- list of strings or anything that behaves like a list of lines
    
    Returns PairCounts object which is a dictionary of {Label: count}
    where label is a string (name from rnaview output), e.g. 'WW--cis',
    and count is an integer representing the base pair count of that
    type of base pair. 
    
    Will raise an error if the number of lines is odd (it expects (multiple)
    label lines followed by the line with the counts.

    If lines is empty it'll return an empty PairCounts object (=empty dict).
    """
    if not len(lines)%2 == 0:
        raise RnaViewParseError("Weird base pair counts format:\n%s"\
            %('\n'.join(lines)))
    res = PairCounts()
    for x in range(0,len(lines),2):
        res.update(zip(lines[x].split(),map(int,lines[x+1].split())))
    return res

def verify_bp_counts(bps, np, pair_counts):
    """Will raise an error on mismatch of reported/actual number of base pairs.
    
    bps -- BasePairs object (the full set from rnaview)
    np -- reported number of pairs (from rnaview output)
    pair_counts -- the dictionary of pair counts

    This function won't return anything. I'll just raise an error
    if the reported and actual numbers don't match.

    The "total base pairs" doesn't have to match the counts
    reported in the dictionary if "X/X" base pairs are present.
    The lines of code dealing with this are removed. Original check was:
    #if np != sum(pair_counts.values()):
    #    raise RnaViewParseError(\
    #    "Reported number of base pairs (%s)"%(np)+\
    #    " doesn't match detailed "+\
    #    "counts (%s)"%(sum(pair_counts.values())))
    """
    subset = bps.select(is_not_stacked_or_tertiary)
    if np != len(subset):
        raise RnaViewParseError(\
        "Reported number of base pairs (%s)"%(np)+\
        " doesn't match number"+\
        " of found base pairs (%s)"%(len(subset)))

def MinimalRnaviewParser(lines):
    """Return line groups for uncommon res, base pairs, base multiples, counts.

    lines -- list of strings or anything that behaves like a list of lines

    This function groups all the lines into particular groups. It recognizes:
    FN: filename, should be a single line
    UC: uncommon residues, lines that start with 'uncommon'
    BP: base pairs, everything between BEGIN_base-pair and END_base_pair
    BM: base multiplets, everything between BEGIN_multiplets and END_multiplets
    NP: number of pairs and bases, single line
    PC: pair counts, lines between '-----------------------'
    """
    result = {'UC':[],'BP':[],'BM':[],'PC':[], 'FN':[], 'NP':[]}
    in_bp = False
    in_bm = False
    in_pc = False
    for line in lines:
        line = line.strip()
        if line:
            if line.startswith('PDB data file name'): # filename
                result['FN'].append(line)
                continue
            if line.startswith('The total base pairs ='): # number of pairs
                result['NP'].append(line)
                continue
            if line.startswith('uncommon'): # uncommon residues
                result['UC'].append(line)
            elif line.startswith('BEGIN_base-pair'): # base pairs
                in_bp = True
            elif line.startswith('END_base-pair'):
                in_bp = False
            elif line.startswith('BEGIN_multiplets'): # base multiplets
                in_bm = True
            elif line.startswith('END_multiplets'):
                in_bm = False
            elif line.startswith('-------'): # pair counts
                in_pc = True
            else:
                if in_bp:
                    result['BP'].append(line)
                elif in_bm:
                    result['BM'].append(line)
                elif in_pc:
                    result['PC'].append(line)
                else:
                    continue
    return result

def RnaviewParser(lines, strict=True):
    """Parse output from the Rnaview program.

    lines -- list of strings or anything that behaves like a list of lines,
        this should be the rnaview output
    strict -- boolean indicating whether an error should be raised or not,
        default=True
    
    The purpose of this parser is to extract every piece of information
    from an Rnaview output file. It parses the filename (FN), uncommon
    residues (UC), base pairs (BP), base multiplets (BM), number of
    base pairs (NP), and the pair counts (PC).
    """
    parsers = {'FN': parse_filename,\
                'UC':parse_uncommon_residues,\
                'BP':parse_base_pairs,\
                'BM':parse_base_multiplets,\
                'PC':parse_pair_counts,
                'NP':parse_number_of_pairs}

    grouped_lines = MinimalRnaviewParser(lines)
    
    # parse the line groups
    result = {}
    for section in parsers.keys():
        try:
            result[section] = parsers[section](grouped_lines[section])
        except RnaViewParseError, e:
            if strict:
                raise RnaViewParseError(str(e))
            else:
                result[section] = None
    
    # verify reported versus actual number of base pairs
    if result['NP'] is not None:
        num_bps_exp = result['NP']['NUM_PAIRS']
        try:
            verify_bp_counts(result['BP'], num_bps_exp, result['PC'])
        except RnaViewParseError, e:
            if strict:
                raise RnaViewParseError(str(e))
            else:
                pass
            
    return result


if __name__ == "__main__":
    pass 
