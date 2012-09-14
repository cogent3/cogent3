#!/usr/bin/env python

"""A Flowgram object for 454 sequencing data."""

__author__ = "Jens Reeder, Julia Goodrich"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jens Reeder","Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jens Reeder"
__email__ = "jreeder@colorado.edu"
__status__ = "Development"

from copy import copy

from cogent.util.unit_test import FakeRandom
from cogent.core.sequence import Sequence

DEFAULT_FLOWORDER = "TACG"
DEFAULT_KEYSEQ = "TCAG"

class Flowgram(object):
    """Holds a 454 flowgram object"""

    HeaderInfo = ["Run Prefix", "Region #","XY Location","Run Name",
                  "Analysis Name", "Full Path","Read Header Len","Name Length",
                  "# of Bases","Clip Qual Left","Clip Qual Right",
                  "Clip Adap Left","Clip Adap Right"]
    FlowgramInfo = ['Flow Indexes','Bases','Quality Scores']
    
    def __init__(self, flowgram = '', Name = None, KeySeq = DEFAULT_KEYSEQ,
                 floworder = DEFAULT_FLOWORDER, header_info = None):
        """Initialize a flowgram.
        
        Arguments:
            flowgram: the raw flowgram string, or list, no other type
            guarenteed to work
            as expected, default is ''
            
            Name: the flowgram name

            KeySeq: the 454 key sequence
            
            floworder: flow sequence used to transform flowgram to sequence
        """

        if Name is None and hasattr(flowgram, 'Name'):
            Name = flowgram.Name
        if Name is None and hasattr(flowgram, 'header_info'):
            header_info = flowgram.header_info
                        
        self.Name = Name

        if hasattr(flowgram, '_flowgram'):
            flowgram  = flowgram._flowgram
        if isinstance(flowgram, str):
            self._flowgram = ' '.join(flowgram.split())
        if isinstance(flowgram, list):
            self._flowgram = ' '.join(map(str, flowgram))
        else: 
            self._flowgram = str(flowgram)

        self.flowgram = map(float,self._flowgram.split())

        self.keySeq = KeySeq

        self.floworder = floworder

        if header_info is not None:
            for i in header_info:
                setattr(self,i,header_info[i])

        #Why do we store the information twice, once as attribute and once in header_info?
        self.header_info = header_info

    def __str__(self):
        """__str__ returns self._flowgram unmodified."""
        return '\t'.join(self._flowgram.split())

    def __len__(self):
        """returns the length of the flowgram"""
        return len(self.flowgram)

    def cmpSeqToString(self, other):
        """compares the flowgram's sequence to other which is a string
        will first try to compare by self.Bases, then by self.toSeq"""
        
        if hasattr(self,'Bases') and self.Bases == other:
            return True
        else:
            return self.toSeq() == other

    def cmpByName(self, other):
        """compares based on the name, other must also be a flowgram object"""
        if self is other:
            return 0
        try:
            return cmp(self.Name, other.Name)
        except AttributeError:
            return cmp(type(self), type(other))

    def cmpBySeqs(self, other):
        """compares by the sequences they represent
            other must also be a flowgram object
        """
        if self is other:
            return 0
        try:
            return cmp(self.Bases,other.Bases)
        except AttributeError:
            return cmp(self.toSeq(), other.toSeq())
  
    def hasProperKey(self, keyseq=DEFAULT_KEYSEQ):
        """Checks for the proper key sequence"""

        keylen = len(keyseq)
        keyseq_from_flow = self.toSeq(truncate=False,
                                            Bases=False)[:keylen]
        return (keyseq_from_flow == keyseq)
         
    def __cmp__(self, other):
        """compares flowgram to other which is a string or another flowgram"""

        if isinstance(other, Flowgram):
            other = other._flowgram
        return cmp(self._flowgram, other)
        
    def __iter__(self):
        """yields successive floats in flowgram"""
        for f in self.flowgram:
            yield f

    def __hash__(self):
        """__hash__ behaves like the flowgram string for dict lookup."""
        return hash(self._flowgram)

    def __contains__(self, other):
        """__contains__ checks whether other is in the flowgram string."""
        return other in self._flowgram

    def toSeq(self, Bases=True, truncate=True):
        """Translates flowgram to sequence and returns sequence object
            if Bases is True then a sequence object will be made using
            self.Bases instead of translating the flowgram
        
            truncate: if True strip off lowercase chars (low quality bases)
            """
        if Bases and hasattr(self, "Bases"):
            seq = self.Bases
        else:
            seq = []
            if self.floworder is None:
                raise ValueError, "must have self.floworder set"
            key = FakeRandom(self.floworder,True)

            flows_since_last = 0
            for n in self.flowgram:
                signal = int(round(n))
                seq.extend([key()]* signal)
                if (signal>0):
                    flows_since_last = 0
                else:
                    flows_since_last += 1
                    if(flows_since_last ==4):
                        seq.extend('N')
                        flows_since_last=0
            seq = ''.join(seq)
            #cache the result for next time
            self.Bases = seq

        if(truncate):
            seq = str(seq)
            seq = seq.rstrip("acgtn")
            seq = seq.lstrip("actgn")
        return Sequence(seq, Name = self.Name)

    def toFasta(self, make_seqlabel=None, LineWrap = 80):
        """Return string in FASTA format, no trailing newline

        Will use self.Bases if it is set otherwise it will translate the
            flowgram
        Arguments:
            - make_seqlabel: callback function that takes the seq object and
              returns a label str
        """
        if hasattr(self,'Bases'):
            seq = self.toSeq(Bases = True)
        else:
            seq = self.toSeq()

        seq.LineWrap = LineWrap
        return seq.toFasta(make_seqlabel = make_seqlabel)
     
    def getQualityTrimmedFlowgram(self):
        """Returns trimmed flowgram according to Clip Qual Right"""
        flow_copy = copy(self)
        if (hasattr(self, "Clip Qual Right") and hasattr(self, "Flow Indexes")):
            clip_right   = int(getattr(self, "Clip Qual Right"))
            flow_indices = getattr(self, "Flow Indexes")
            flow_indices = [int(k) for k in flow_indices.split('\t') if k != '']
            
            clip_right_flowgram = flow_indices[clip_right-1]
            #Truncate flowgram
            flow_copy.flowgram = self.flowgram[:clip_right_flowgram]
            flow_copy._flowgram =\
                "\t".join(self._flowgram.split()[:clip_right_flowgram])

            #Update attributes
            if hasattr(flow_copy, "Quality Scores"):
                qual_scores = getattr(flow_copy,"Quality Scores").split('\t')
                setattr(flow_copy, "Quality Scores",
                        "\t".join(qual_scores[:clip_right]))

            if hasattr(flow_copy, "Flow Indexes"):
                setattr(flow_copy, "Flow Indexes", 
                        "\t".join(map(str, flow_indices[:clip_right])))

            if hasattr(flow_copy, "Bases"):
                flow_copy.Bases = self.Bases[:clip_right]

            if hasattr(flow_copy, "# of Bases"):
                setattr(flow_copy, "# of Bases", clip_right)
 
        return flow_copy
                
    def getPrimerTrimmedFlowgram(self, primerseq):
        """Cuts the key and primer sequences of a flowgram.
        
        primerseq: the primer seq to be truncated from flowgram
        """
      
        if(primerseq==""):
            return self
        else:
            flow_copy = copy(self)
            #Key currently not reliable set by FlowgramCollection
            #instead pass key as part of primer
            #key = flow_copy.keySeq or ""
         
            flow_indices = getattr(self, "Flow Indexes")
            flow_indices = [int(k) for k in flow_indices.split('\t') if k != '']

            #position of last primer char in flowgram
            primer_len = len(primerseq)
            pos = flow_indices[primer_len-1]
            signal = flow_copy.flowgram[pos-1] 

            if (signal < 0.5):
                #Flowgram is not consistent with primerseq
                return None

            elif (signal < 1.5):
                pad_num = pos % 4
                #we can simply cut off
                flow_copy.flowgram = flow_copy.flowgram[pos:]
                # and pad flowgram to the left to sync with floworder
                flow_copy.flowgram[:0] = pad_num*[0.00]
                #check that first 4 flows not are all zero

            else:
                pad_num = (pos-1)%4
                # we are cutting within a signal, need to do some flowgram arithmetic
                lastchar = primerseq[-1]
                #get the position in the homopolyemer 
                pos_in_homopoly = len(primerseq) - len(primerseq.rstrip(lastchar))
                                                    
                flow_copy.flowgram = flow_copy.flowgram[pos-1:]
                flow_copy.flowgram[0] =  max(0.00, flow_copy.flowgram[0] - pos_in_homopoly)
                #pad flowgram to the left to sync with floworder
                flow_copy.flowgram[:0] = (pad_num)*[0.00]

            # delete first flow cycle if all <0.5 (otherwise an N would be called)
            if(any([sign>=0.5 for sign in flow_copy.flowgram[:4]])):
                #We are ok
                extra_shift=0
                pass
            else:
                #we truncate the first 4 flows
                flow_copy.flowgram = flow_copy.flowgram[4:]
                extra_shift=4
            #Update "Flow Indexes" attribute
            #shift all flow indices by the deleted amount
            # WARNING: this sets wrong flow indexes, so better set to nothing     
#            setattr(flow_copy, "Flow Indexes", 
#                    "\t".join([ str(a-(pos+extra_shift)+pad_num) for a in\
#                                    flow_indices[primer_len:]]))
            setattr(flow_copy, "Flow Indexes", "")
            #Update flowgram string representation
            flow_copy._flowgram = "\t".join(map(lambda a:"%.2f"%a,
                                                flow_copy.flowgram)) 
            #Update "Quality Scores" attribute
            if hasattr(self, "Quality Scores"):
                qual_scores = getattr(flow_copy,"Quality Scores").split('\t')
                setattr(flow_copy, "Quality Scores",  "\t".join(qual_scores[primer_len:]))
            #Update Bases attribute
            if hasattr(flow_copy, "Bases"):
                flow_copy.Bases = flow_copy.Bases[primer_len:]         
            #Update "# of Bases" attribute
            if hasattr(flow_copy, "# of Bases"):
                setattr(flow_copy, "# of Bases",
                        str(int(getattr(flow_copy, "# of Bases")) - (primer_len)))  
            if hasattr(flow_copy, "Clip Qual Left"):
                setattr(flow_copy, "Clip Qual Left", str(max(0, int(getattr(flow_copy, "Clip Qual Left")) - primer_len)))  
            if hasattr(flow_copy, "Clip Qual Right"):
                setattr(flow_copy, "Clip Qual Right", str(max(0, int(getattr(flow_copy, "Clip Qual Right")) - primer_len)))  

            if hasattr(flow_copy, "Clip Adap Left"):
                setattr(flow_copy, "Clip Adap Left", str(max(0, int(getattr(flow_copy, "Clip Adap Left")) - primer_len)))  
            if hasattr(flow_copy, "Clip Adap Right"):
                setattr(flow_copy, "Clip Adap Right", str(max(0, int(getattr(flow_copy, "Clip Adap Right")) - primer_len))) 

            return flow_copy

    def createFlowHeader(self):
        """header_info dict turned into flowgram header"""
        lines = [">%s\n"%self.Name]
        flow_info = []
        head_info = []
        for i in self.FlowgramInfo:
            if hasattr(self,i):
                flow_info.append('%s:\t%s\n' %
                                 (i, getattr(self, i)))
        for i in self.HeaderInfo:
            if hasattr(self,i):
                head_info.append('  %s:\t%s\n' %
                                 (i, getattr(self, i)))

        lines.extend(head_info)
        lines.extend(flow_info)
        lines.append("Flowgram:\t%s" % str(self))

        return (''.join(lines)+"\n")

def seq_to_flow(seq, id = None, keyseq = None, floworder = DEFAULT_FLOWORDER):
    """ Transform a sequence into an ideal flow.

    seq: sequence to transform to flowgram
    
    id: identifier

    keyseq
    """

    complete_flow = floworder * len(seq) # worst case length
    i = 0  # iterates over seq
    j = 0  # iterates over the flow sequence tcagtcagtcag...
    mask = ""

    while (j < len(complete_flow) and i<len(seq)):
        if(seq[i] == 'N'):
            mask+="0.00 0.00 0.00 0.00"
            i+=1
            j+=4
            continue
        if (complete_flow[j] == seq[i]):
            #check for more than one of this nuc
            c = 1
            i += 1
            while(i < len(seq) and j < len(complete_flow)\
                  and complete_flow[j]==seq[i]):
                i += 1
                c += 1
            mask += "%d" % c
            if (i >= len(seq)):
                break
            j += 1
        else:
            mask += "0"
            j += 1

    # pad mask to finish the last flow to a multiple of the floworder length
    if (len(mask) % len(floworder) != 0):
        right_missing = len(floworder) - (len(mask) % len(floworder))
        mask += "0" * right_missing

    return Flowgram(map(float, mask), id, keyseq, floworder)

def build_averaged_flowgram(flowgrams):
    """Builds an averaged flowgram from a list of raw signals."""

    result=[]
    if(len(flowgrams)==1):
        return flowgrams[0]
    for tuple in map(None, *flowgrams):
        k=0
        sum=0
        for element in tuple: 
            if (element!=None):
                k+=1
                sum +=element
        result.append(round(sum/k,2))
    return result
