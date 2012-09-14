#! /usr/bin/env python
"""A trie and compressed trie data structure."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Prototype"

##################
# A compressed Trie
##################

class _Compressed_Node:
    """ A node in a compressed Trie."""
    
    def __init__(self, key, labels=None):
        """Creates a new Node."""

        self.labels = labels or []
        self.key = key
        self.children = {}
  
    def __nonzero__(self):
        """Checks if Node contains any data."""
        return (self.key!="" or len(self.labels)>0
                or len(self.children.keys())>0)

    def __len__(self):
        """Counts the number of sequences in the Trie"""

        l = 0 
        for n in self.children.values():
            l += len(n)
        l += len(self.labels)
        return l

    def _to_string(self, depth=0):
        """Debugging method to display the Node's content.

        depth: indentation level
        """

        s = ["\n"+depth*'\t'+"key %s"%self.key]
        
        if (self.labels):
            s.append("%s" % str(self.labels))
        s.append("\n")
 
        for n in self.children.values():
            s.append(depth*'\t' +"{%s" %(n._to_string(depth+1)))
            s.append(depth*'\t'+"}\n")
        return "".join(s)

    def size(self):
        """Returns number of nodes."""
        s = 1 
        for n in self.children.values():
            s += n.size()
        return s

    def insert(self, key, label):
        """Insert a key into the Trie.

        key: The key of the entry
        label: The value that should be stored for this key
        """

        node_key_len = len(self.key)
        length = min(node_key_len, len(key))
        #follow the key into the tree
        index=0
        while index < length:
            if(key[index] != self.key[index]):
                #split up node
                new_key_node = _Compressed_Node(key[index:], [label])

                old_key_node = _Compressed_Node(self.key[index:], self.labels)
                old_key_node.children = self.children
                    
                self.children = {key[index]: new_key_node,
                                 self.key[index]: old_key_node}
                self.key = self.key[:index]
                self.labels = []
                return
            index += 1
        #new key matches node key exactly
        if (index == len(self.key) and index == len(key)):
            self.labels.append(label)
            return

        #Key shorter than node key
        if(index < node_key_len):
            lower_node = _Compressed_Node(self.key[index:], self.labels)
            lower_node.children = self.children
            
            self.children = {self.key[index]: lower_node}
            self.key = key
            self.labels = [label]
            return
  
        #new key longer than current node key
        node = self.children.get(key[index])
        if(node):
            # insert into next node
               node.insert(key[index:], label)
        else:
            # create new node
            new_node = _Compressed_Node(key[index:], [label])
            self.children[key[index]] = new_node
        return   

    def find(self, key):
        """Searches for key and returns values stored for the key.

        key: the key whose values should be retuned.
        """

        #key exhausted
        if(len(key) == 0):
            return self.labels
        
        #find matching part of key and node_key
        min_length = min(len(key), len(self.key))
        index = 0   
        while index < min_length:
            if(key[index] != self.key[index]):
                return []
            index+=1

        #key and node_key match exactly
        if (index == len(key)):
            return self.labels

        node = self.children.get(key[index])
        if(node):
            # descend to next node
            return node.find(key[index:])
        else:
            return []

    def prefixMap(self):
        """Builds a prefix map from sequences stored in Trie."""
        
        labels = []
        mapping = {}
        
        # we have a leaf
        if (len(self.children)==0):
            mapping =  {self.labels[0]: self.labels[1:]}
        # we are at an internal node
        else: 
            for child in self.children.values():
                mapping.update(child.prefixMap())
            # get largest group
            n = -1
            key_largest = None
            for (key, value) in mapping.iteritems():
                if (len(value) > n):
                    n = len(value)
                    key_largest = key
            # append this node's values        
            mapping[key_largest].extend(self.labels)

        return mapping
     

class Compressed_Trie:
    """ A compressed Trie"""

    def __init__(self):
        """ Return an empty Trie."""
        self.root = _Compressed_Node("")
        
    def insert(self, key, label):
        """Insert key with label in Trie."""
        self.root.insert(key, label)

    def find(self, key):
        """Returns label for key in Trie."""
        return self.root.find(key)
        
    def __str__(self):
        """Ugly Trie string representation for debugging."""
        return self.root._to_string()
 
    def __nonzero__(self):
        """Checks wheter the trie is empty or not."""
        return (self.root.__nonzero__())
    
    def __len__(self):
        """Returns number of sequences in Trie."""
        return len(self.root)

    def size(self):
        """Returns number of nodes in Trie"""
        return self.root.size()

    def prefixMap(self):
        """builds a prefix map of seqeunces stored in Trie."""
        return self.root.prefixMap()


###################
# An atomic Trie
###################

class _Node:
    """ A node in a Trie."""

    def __init__(self, labels=None):
        """Creates a new node with value label."""
        self.labels = labels or []
        self.children = {}

    def insert(self, key, label):
        """Insert key with value label into Trie.

        key: The key of the entry
        label: The value that should be stored for this key
        """
        curr_node = self
        for ch in key:
            curr_node = curr_node.children.setdefault(ch, _Node())
        
        curr_node.labels.append(label)

    def _insert_unique(self, key, value):
        """Insert key with value if key not already in Trie.

        key: The key of the entry
        label: The value that should be stored for this key
        
        Returns label of key if unique or label of containing key.
        Note: Tries built with this function will only have labels
              on their leaves. 
        """
        
        curr_node = self
        labels = []
        #descend along key and collect all internal node labels
        for ch in key:
            if curr_node.labels != []:
                labels.extend(curr_node.labels)
                curr_node.labels = []
            curr_node = curr_node.children.setdefault(ch, _Node())

        if (curr_node.children=={} and curr_node.labels==[]):
            # we have a new prefix
            curr_node.labels = [value]
            return (curr_node.labels[0], labels)
        else:
            # we are at an internal node or a labeled leaf
            # descend to the "leftmost" leaf 
            while (curr_node.children != {}):
                curr_node = curr_node.children.values()[0]
            return (curr_node.labels[0], [])

    def find(self, key):
        """Retrieves node with key in Trie."""

        next_node = self
        for ch in key:
            try:
                next_node = next_node.children[ch]
            except KeyError:
                return None
        return next_node

class Trie:
    """ A Trie data structure."""

    def __init__(self):
        """Inits an empty Trie."""
        self.root = _Node()
        
    def insert(self, key, label):
        """Insert key with label into Trie.

        key: The key of the entry
        label: The value that should be stored for this key
        """
        return self.root.insert(key, label)
    
    def _insert_unique(self, key, label):
        """Insert unique keys into Trie, skip redundant ones.

        key: The key of the entry
        label: The value that should be stored for this key

        Returns: list of labels removed rom Trie
        """
        return self.root._insert_unique(key, label)
    
    def find(self, key):
        """Retrieves  labels of key in Trie."""
        node = self.root.find(key)

        if (node ==None):
            return []
        return node.labels

def _build_prefix_map(seqs):
    """Builds a prefix map using an atomic Trie.

    seqs: dict like sequence collection
    """
    mapping = {}
    t = Trie()
    for label,seq in seqs:
        (label2, prefix_labels) = t._insert_unique(seq, label)
        if (label2==label):
            #seq got inserted and is new word
            mapping[label] = prefix_labels[:]
            #delete old mappings and transfer old mappings to new label
            for l in prefix_labels:
                mapping[label].extend(mapping[l])
                del(mapping[l])
        else:
            # seq is already in tree and is prefix of label2
            mapping[label2].append(label)

    return mapping

def build_trie(seqs, classname=Compressed_Trie):
    """ build a Trie for a list of (label,seq) pairs.

    seqs: list of (label,seq) pairs

    classname: Constructor to use for tree building
               Compressed_Trie (default) or Trie
     `"""
    if (not(classname==Trie or classname==Compressed_Trie)):
        raise ValueError, "Wrong classname for build_trie."
    t = classname()
    for (label, seq) in seqs:
            t.insert(seq, label)
    return t


def build_prefix_map(seqs, classname=Compressed_Trie):
    """Builds prefix map from seqs.
 
    seqs: list of (label,seq) pairs

    classname: Constructor to use for tree building
               Compressed_Trie (default) or Trie

    This method can be used to filter sequences for identity and for identical prefixes.
    Due to the nature of tries a list of n seqs of length l can be filtered in O(nl) time.
    """
    if (not(classname==Trie or classname==Compressed_Trie)):
        raise ValueError, "Wrong classname for build_trie."
    if (classname == Compressed_Trie):
        t = build_trie(seqs)
        return t.prefixMap()
    else:
        return _build_prefix_map(seqs)
