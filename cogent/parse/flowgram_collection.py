
__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jens Reeder","Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Julia Goodrich"
__email__ = "julia.goodrich@colorado.edu"
__Status__ = "Development"


from copy import copy
from types import GeneratorType

from numpy import transpose
from numpy.random import multinomial

from cogent.util.unit_test import FakeRandom
from cogent.core.sequence import Sequence
from cogent.parse.flowgram_parser import parse_sff
from cogent.parse.flowgram import Flowgram
from cogent.core.alignment import SequenceCollection

default_floworder = "TACG"
default_keyseq = "TCAG"

def assign_sequential_names(ignored, num_seqs, base_name='seq', start_at=0):
    """Returns list of num_seqs sequential, unique names.
    
    First argument is ignored; expect this to be set as a class attribute.
    """
    return ['%s_%s' % (base_name,i) for i in range(start_at,start_at+num_seqs)]

def flows_from_array(a):
    """FlowgramCollection from array of pos x seq: names are integers.
    
    This is an InputHandler for FlowgramCollection. It converts an arbitrary
    array of numbers into Flowgram objects and leaves the flowgrams unlabeled.
    """
    return list(transpose(a)), None, None

def flows_from_flowCollection(flow):
    """FlowgramCollection from array of pos x seq: names are integers.
    
    This is an InputHandler for FlowgramCollection. It converts an arbitrary
    array of numbers into Flowgram objects and leaves the flowgrams unlabeled.
    """
    return flow.flows, flow.Names, [f.header_info for f in flow.flows]

def flows_from_kv_pairs(flows):
    """SequenceCollection from list of (key, val) pairs.

        val can be str or flowgram object
    """
    
    names, flows =  map(list, zip(*flows))
    if isinstance(flows[0], str):
        info = [None]*len(flows)
    else:
        info = [f.header_info for f in flows]
    
    return flows, names, info

def flows_from_empty(obj, *args, **kwargs):
    """SequenceCollection from empty data: raise exception."""
    raise ValueError, "Cannot create empty SequenceCollection."

def flows_from_dict(flows):
    """SequenceCollection from dict of {label:flow_as_str} or {label:flow_obj}.
    """
    names, flows = map(list, zip(*flows.items()))
    if isinstance(flows[0], str):
        info = [None]*len(flows)
    else:
        info = [f.header_info for f in flows]
    return flows, names, info

def flows_from_sff(flows):
    """lines is sff file lines.
    """
    if isinstance(flows, str):
        flows = flows.splitlines()
        
    flows, head = parse_sff(flows)
    return flows_from_generic(flows)

def flows_from_generic(flows):
    """SequenceCollection from generic seq x pos data: seq of seqs of chars.
    
    This is an InputHandler for SequenceCollection. It converts a generic list
    (each item in the list will be mapped onto an object using
    seq_constructor and assigns sequential integers (0-based) as names.
    """
    names = []
    info = []
    for f in flows:
        if hasattr(f, 'Name'):
            names.append(f.Name)
        else:
            names.append(None)
        if hasattr(f, 'header_info'):
            info.append(f.header_info)
        else:
            info.append(None)
            
    return flows, names, info


class FlowgramCollection(object):
    """stores Flowgrams.
    
    - InputHandlers: support for different data types
    - flows: behaves like list of Flowgram objects
    - Names: behaves like list of names for the Flowgram objects
    - NamedFlows: behaves like dict of {name:flow}
    """    
    InputHandlers = {   'array': flows_from_array,
                        'dict': flows_from_dict,
                        'sff' :flows_from_sff,
                        'empty': flows_from_empty,
                        'kv_pairs':flows_from_kv_pairs,
                        'flowcoll':flows_from_flowCollection,
                        'generic':flows_from_generic
                    }
    HeaderInfo = ['Magic Number', 'Version', 'Index Offset', 'Index Length',
                  '# of Reads', 'Header Length', 'Key Length', '# of Flows',
                  'Flowgram Code', 'Flow Chars', 'Key Sequence']

    DefaultNameFunction = assign_sequential_names
    
    def __init__(self, data, Name = None, Names = None,
                 header_info = None, conversion_f = None,
                 name_conversion_f = None, remove_duplicate_names = False):
        """Initialize self with data and optionally Info.
        

        Parameters:
        
        data:           Data to convert into a FlowgramCollection
        
        Name:           Name of the FlowgramCollection.
                
        conversion_f:   Function to convert data into Flowgram.
                        
        name_conversion_f: if present, converts name into f(name).

        header_info: contains info to be printed in the common header of an sff
                    file, it is a dictionary. ex: Key Sequence:"ATCG"
        """

        #read all the data in if we were passed a generator
        if isinstance(data, GeneratorType):
            data = list(data)

        #set the Name
        self.Name = Name

        if header_info is not None:
            self._check_header_info(header_info)
        
            for i in header_info:
                setattr(self,i,header_info[i])
            if 'Key Sequence' in header_info:
                keyseq = header_info['Key Sequence']
            else:
                keyseq = None
            if 'Flow Chars' in header_info:
                floworder = header_info['Flow Chars']
            else:
                floworder = default_floworder
        else:
            keyseq = None
            floworder = default_floworder


        self.header_info = header_info

        per_flow_names, flows, name_order, info = \
                self._names_flows_order(conversion_f, data, Names, \
                name_conversion_f, remove_duplicate_names)
        self.Names = name_order

        #will take only the flows and names that are in name_order
        if per_flow_names != name_order:
            good_indices = []
            for n in name_order:
                good_indices.append(per_flow_names.index(n))
            flows = [flows[i] for i in good_indices]
            info = [info[i] for i in good_indices]
            per_flow_names = name_order

        self.flow_str= flows  
        self.flows = [Flowgram(f,n,keyseq,floworder, i)\
                      for f,n, i in zip(flows,self.Names,info)]

        #create NamedFlows dict for fast lookups
        self.NamedFlows = self._make_named_flows(self.Names, self.flows)
                

    def _strip_duplicates(self, names, flows, info):
        """Internal function to strip duplicates from list of names"""
        if len(set(names)) == len(names):
            return set(), names, flows, info
        #if we got here, there are duplicates
        unique_names = {}
        duplicates = {}
        fixed_names = []
        fixed_flows = []
        fixed_info = []
        for n, f, i in zip(names, flows,info):
            if n in unique_names:
                duplicates[n] = 1
            else:
                unique_names[n] = 1
                fixed_names.append(n)
                fixed_flows.append(f)
                fixed_info.append(i)
                
        return duplicates, fixed_names, fixed_flows, fixed_info


    def _names_flows_order(self, conversion_f, data, Names, \
            name_conversion_f, remove_duplicate_names):
        """Internal function to figure out names, flows, and name_order."""
        #figure out conversion function and whether it's an array
        if not conversion_f:
            input_type = self._guess_input_type(data)
            conversion_f = self.InputHandlers[input_type]
        #set seqs, names, and handler_info as properties
        
        flows, names, info = conversion_f(data)
        if names and name_conversion_f:
            names = map(name_conversion_f, names)

        #if no names were passed in as Names, if we obtained them from
        #the seqs we should use them, but otherwise we should use the
        #default names
        if Names is None:
            if (names is None) or (None in names):
                per_flow_names = name_order = \
                    self.DefaultNameFunction(len(flows))
            else:   #got names from seqs
                per_flow_names = name_order = names
        else:
            #otherwise, names were passed in as Names: use this as the order
            #if we got names from the sequences, but otherwise assign the
            #names to successive sequences in order
            if (names is None) or (None in names):
                per_flow_names = name_order = Names
            else: #got names from seqs, so assume name_order is in Names
                per_flow_names = names
                name_order = Names
        #check for duplicate names
        duplicates, fixed_names, fixed_flows, fixed_info = \
            self._strip_duplicates(per_flow_names, flows, info)
        if duplicates:
            if remove_duplicate_names:
                per_flow_names, flows, info =fixed_names,fixed_flows,fixed_info
                #if name_order doesn't have the same names as per_seq_names,
                #replace it with per_seq_names
                if (set(name_order) != set(per_flow_names)) or\
                    (len(name_order) != len(per_flow_names)):
                    name_order = per_flow_names
            else:
                raise ValueError, \
                "Some names were not unique. Duplicates are:\n" + \
                str(sorted(duplicates.keys()))
        return per_flow_names, flows, name_order, info


    def _check_header_info(self, info):
        for h in info:
            if h not in self.HeaderInfo:
                raise ValueError, "invalid key in header_info"
        
    def _make_named_flows(self, names, flows):
        """Returns NamedFlows: dict of name:flow."""
        name_flow_tuples = zip(names, flows)
        for n, f in name_flow_tuples:
            f.Name = n
        return dict(name_flow_tuples)
    
    def _guess_input_type(self, data):
        """Guesses input type of data; returns result as key of InputHandlers.
               
        Returns 'empty' if check fails, i.e. if it can't recognize the data
        as a specific type. Note that bad data is not guaranteed to
        return 'empty', and may be recognized as another type incorrectly.
        """
        
        if isinstance(data, dict):
            return 'dict'
        if isinstance(data, str):
            return 'sff'
        if isinstance(data, FlowgramCollection):
            return 'flowcoll'
        
        first = None
        try:
            first = data[0]
        except (IndexError, TypeError):
            pass
        try:
            first = iter(data).next()
        except (IndexError, TypeError, StopIteration):
            pass
        if first is None:
            return 'empty'
        try:
            if isinstance(first, Flowgram):     #model sequence base type
                return 'generic'
            if hasattr(first, 'dtype'):    #array object
                return 'array'
            elif isinstance(first, str) and first.startswith('Common'):
                return 'sff'
            else:
                try:
                    dict(data)
                    return 'kv_pairs'
                except (TypeError, ValueError):
                    pass
            return "generic"
        except (IndexError, TypeError), e:
            return 'empty'

    def __len__(self):
        """returns the number of flowgrams in the collection"""
        return len(self.flows)

    def __iter__(self):
        """iterates over the flows in the collection"""
        for f in self.flows:
            yield f

    def __str__(self):
        """returns string like sff file given the flowgrams and header_info"""
        lines = self.createCommonHeader()
        lines.append('')
        lines.extend([f.createFlowHeader() for f in self.flows])
        return '\n'.join(lines)

    def __cmp__(self, other):
        """cmp first tests as dict, then as str."""
        c = cmp(self.NamedFlows, other)
        if not c:
            return 0
        else:
            return cmp(str(self), str(other))

    def keys(self):
        """keys uses self.Names
        
        Note: returns copy, not original.
        """
        return self.Names[:]
    
    def values(self):
        """values returns values corresponding to self.Names."""
        return [self.NamedFlows[n] for n in self.Names]

    def items(self):
        """items returns (name, value) pairs."""
        return [(n, self.NamedFlows[n]) for n in self.Names]
    
    def writeToFile(self, filename=None, **kwargs):
        """Write the flowgrams to a file.
        
        Arguments:
        - filename: name of the output file
        
        """
        
        if filename is None:
            raise DataError('no filename specified')
        
        f = open(filename, 'w')
        f.write(str(self))
        f.close()

    def createCommonHeader(self):
        """header_info dict turned into flowgram common header"""
        lines = ["Common Header:"]
        
        if self.header_info is not None:
            lines.extend(["  %s:\t%s" % (param,self.header_info[param]) \
                      for param in self.header_info])

        return lines

    def toFasta(self, exclude_ambiguous = False,Bases=False,
                make_seqlabel = None):
        """Return flowgram collection in Fasta format
        
        Arguments:
            - make_seqlabel: callback function that takes the seq object and
              returns a label str

        if Bases is True then a fasta string will be made using
            self.Bases instead of translating the flowgram
              
        """
        if exclude_ambiguous:
            flows = flows.omitAmbiguousFlows()
        else:
            flows = self

        seqs = flows.toSequenceCollection(Bases)
            
        return seqs.toFasta(make_seqlabel = make_seqlabel)

    def toPhylip(self, exclude_ambiguous = False,Bases=False,generic_label=True,
                 make_seqlabel=None):
        """
        Return alignment in PHYLIP format and mapping to sequence ids
        
        raises exception if invalid alignment
        
        Arguments:
            - make_seqlabel: callback function that takes the seq object and
              returns a label str

        if Bases is True then a fasta string will be made using
            self.Bases instead of translating the flowgram 
        """
        if exclude_ambiguous:
            flows = flows.omitAmbiguousFlows()
        else:
            flows = self

        seqs = flows.toSequenceCollection( Bases)
        return seqs.toPhylip(make_seqlabel = make_seqlabel,
                             generic_label= generic_label)
    
    def toNexus(self,seq_type, exclude_ambiguous = False, Bases = False,
                interleave_len=50):
        """
        Return alignment in NEXUS format and mapping to sequence ids
        
        **NOTE** Not that every sequence in the alignment MUST come from
            a different species!! (You can concatenate multiple sequences from
            same species together before building tree)
        
        seq_type: dna, rna, or protein
        
        Raises exception if invalid alignment
        if Bases is True then a fasta string will be made using
        self.Bases instead of translating the flowgram
        """

        if exclude_ambiguous:
            flows = flows.omitAmbiguousFlows()
        else:
            flows = self

        seqs = flows.toSequenceCollection(Bases)
        return seqs.toNexus(seq_type,interleave_len = interleave_len)
    
    def toSequenceCollection(self, Bases = False):
        names = self.Names
        flow_dict = self.NamedFlows
        flows = [flow_dict[f].toSeq(Bases = Bases) for f in names]
        return SequenceCollection(flows)

    def addFlows(self, other):
        """Adds flowgrams from other to self. Returns a new object.
        
        other must be of same class as self or coerceable to that class..
        """
        assert not isinstance(other, str), "Must provide a series of flows "+\
                                            "or an flowgramCollection"
        self_flow_class = self.flows[0].__class__
        try:
            combined = self.flows + other.flows
        except AttributeError:
            combined = self.flows + list(other)
        try:
            combined_info =copy(self.header_info)
            combined_info.update(other.header_info)
        except AttributeError:
            combined_info =self.header_info
        
        
        for flow in combined:
            assert flow.__class__ == self_flow_class,\
                "classes different: Expected %s, Got %s" % \
                    (flow.__class__, self_flow_class)
        return self.__class__(data=combined,header_info=combined_info)

    def iterFlows(self, flow_order=None):
        """Iterates over values (sequences) in the alignment, in order.
        
        seq_order: list of keys giving the order in which seqs will be returned.
        Defaults to self.Names. Note that only these sequences will be
        returned, and that KeyError will be raised if there are sequences
        in order that have been deleted from the Alignment. If self.Names
        is None, returns the sequences in the same order as
        self.NamedSeqs.values().
        
        Use map(f, self.seqs()) to apply the constructor f to each seq. f must
        accept a single list as an argument.
        
        Always returns references to the same objects that are values of the
        alignment.
        """
        ns = self.NamedFlows
        get = ns.__getitem__
        for key  in flow_order or self.Names:
            yield get(key)


    def iterItems(self, flow_order=None, pos_order=None):
        """Iterates over elements in the flowgram collection.
        
        seq_order (names) can be used to select a subset of seqs.
        pos_order (positions) can be used to select a subset of positions.
        
        Always iterates along a seq first, then down a position (transposes
        normal order of a[i][j]; possibly, this should change)..
        
        WARNING: FlowgramCollection.iterItems() is not the same as fc.iteritems()
        (which is the built-in dict iteritems that iterates over key-value
        pairs).
        """
        if pos_order:
            for row in self.iterFlows(flow_order):
                for i in pos_order:
                    yield row[i]
        else:
            for row in self.iterFlows(flow_order):
                for i in row:
                    yield i

    Items = property(iterItems)
                    
    def copy(self):
        """Returns deep copy of self."""
        result = self.__class__(self)

    def _take_flows(self): return list(self.iterFlows())
    
    Flows = property(_take_flows)   #access as attribute if using default order.
    
    def takeFlows(self, flows, negate=False, **kwargs):
        """Returns new FlowgramCollection containing only specified flows.
        
        Note that the flows in the new collection will be references to the
        same objects as the seqs in the old collection.
        """
        get = self.NamedFlows.__getitem__
        result = {}
        if negate:
            #copy everything except the specified seqs
            negated_names = []
            row_lookup = dict.fromkeys(flows)
            for r, row in self.NamedFlows.items():
                if r not in row_lookup:
                    result[r] = row
                    negated_names.append(r)
            flows = negated_names #remember to invert the list of names
        
        else:
            #copy only the specified seqs
            for r in flows:
                result[r] = get(r)
        if result:
            return self.__class__(result, Names=flows, **kwargs)
        else:
            return {}   #safe value; can't construct empty collection

    def getFlowIndices(self, f, negate=False, Bases = False):
        """Returns list of keys of flows where f(Seq) is True.
        
        List will be in the same order as self.Names, if present.
        If Bases is True it will use the flowgram Bases attribute otherwise
        it uses the translation to a sequence
        """
        get = self.NamedFlows.__getitem__
        #negate function if necessary
        if negate:
            new_f = lambda x: not f(x)
        else:
            new_f = f
        #get all the seqs where the function is True
        return [key for key in self.Names \
               if new_f(get(key).toSeq(Bases = Bases))]
    
    def takeFlowsIf(self, f, negate=False, Bases = False, **kwargs):
        """Returns new Alignment containing seqs where f(row) is True.
        
        Note that the seqs in the new Alignment are the same objects as the
        seqs in the old Alignment, not copies.

        If Bases is True it will use the flowgram Bases attribute otherwise
        it uses the translation to a sequence
        """
        #pass negate to get SeqIndices
        return self.takeFlows(self.getFlowIndices(f, negate, Bases= Bases),
                              **kwargs)

    def getFlow(self, flowname):
        """Return a flowgram object for the specified seqname.
        """
        return self.NamedFlows[flowname]

    def getFlowNames(self):
        """Return a list of Flowgram names."""
        return self.Names[:]

    def getIntMap(self,prefix='seq_'):
        """Returns a dict with names mapped to enumerates integer names.
            
            - prefix: prefix for sequence label. Default = 'seq_'
            - int_keys is a dict mapping int names to sorted original names.
        """
        get = self.NamedFlows.__getitem__
        int_keys = dict([(prefix+str(i),k) for i,k in \
                enumerate(sorted(self.NamedFlows.keys()))])
        int_map = dict([(k, copy(get(v))) for k,v in int_keys.items()])
        return int_map, int_keys

    def toDict(self):
        """Returns the collection as dict of names -> strings.

        Note: returns strings, NOT Flowgram objects.
        """
        collection_dict = {}
        
        for flow_name in self.Names:
            collection_dict[flow_name] = self.NamedFlows[flow_name]._flowgram
        
        return collection_dict
    
    def omitAmbiguousFlows(self, Bases = False):
        """Returns an object containing only the sequences without N's"""
        is_ambiguous = lambda x: 'N' not in x
        
        return self.takeFlowsIf(is_ambiguous, Bases = Bases)

    def setBases(self):
        """Sets the Bases property for each flowgram using toSeq"""
        for f in self.values():
            f.Bases = f.toSeq()
    

def pick_from_prob_density(pvals,bin_size):
    l = multinomial(1,pvals)
    return (l.nonzero()[0][0]) * bin_size

def seqs_to_flows(seqs, keyseq = default_keyseq, floworder = default_floworder,
                  numflows = None, probs = None, bin_size = 0.01,
                  header_info = {}):
    """ Transfrom a sequence into an ideal flow
        seqs: a list of name sequence object tuples (name,tuple)
        keyseq: the flowgram key Sequence
        floworder: The chars needed to convert seq to flow
        numflows: number of total flows in each flowgram, if it is specified
            the flowgram will be padded to that number
        probs: dictionary defining the probability distribution for each
            homopolymer

        WARNING:each distributions probabilities must add to 1.0
    """
    flows = []
    homopolymer_counter = 1.0
    if probs:
        for p in probs:
            if round(sum(probs[p]),1) != 1.0:
                raise ValueError, 'probs[%s] does not add to 1.0' % p

    for name,seq in seqs:
        flow_seq = FakeRandom(floworder,True)
        flow = []
        seq_len = len(seq)
        for i, nuc in enumerate(seq):
            if i < seq_len-1 and seq[i+1] == nuc:
                homopolymer_counter += 1.0
            else:
                while flow_seq() != nuc:
                    if probs is None:
                        val = 0.0
                    else:
                        val = pick_from_prob_density(probs[0],bin_size)
                    flow.append(val)
                if (probs is None) or (homopolymer_counter > 9):
                    val = homopolymer_counter
                else:
                    val = pick_from_prob_density(probs[int(homopolymer_counter)],bin_size)
                 
                flow.append(val)
                homopolymer_counter = 1.0

        len_flow = len(flow)
        len_order = len(floworder)

        if numflows is not None and numflows % len_order != 0:
            raise ValueError, "numflows must be divisable by the length of floworder"
        if (len_flow % len_order != 0):
            right_missing = len_order - (len_flow % len_order)
            if numflows != (len_flow + right_missing) and numflows is not None:
                right_missing += (numflows - (len_flow+right_missing))
            if probs is None:
                flow.extend([0.0]*right_missing)
            else:
                for i in range(0,right_missing):
                    flow.append(pick_from_prob_density(probs[0],bin_size))

        flows.append((name, Flowgram(flow, id, keyseq, floworder)))
    if keyseq is not None:
        keylen = len(keyseq)
    else:
        keylen = None
    header_info.update({'Key Sequence':keyseq,'Flow Chars':floworder,
                        'Key Length':keylen})
    return FlowgramCollection(flows, header_info = header_info)
