#!/usr/bin/env python
"""Fast tree class for sequence simulations."""

from operator import add, or_, and_
from random import choice
from numpy import array, zeros, transpose, arange, concatenate, any
from numpy.random import permutation
from cogent.core.tree import PhyloNode

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class RangeNode(PhyloNode):
    """Node object that assigns ids to its leaves and can access leaf blocks.
    
    Note: some of these methods should possibly move to the base class.
    """
    def __init__(self, *args, **kwargs):
        """Returns a new RangeNode object.

        Name: text label
        LeafRange: range of this node's leaves in array. last = index+1
        Id: index of this node in array.
        Children: list of Node objects that are this node's direct 
                  children.
        Parent: Node object that is this node's parent.
        NameLoaded: From cogent.core.tree.TreeNode, undocumented.

        WARNING: Parent/Child relationships are _not_ checked to preserve
        consistency! You must specify both the parent and the children
        explicitly, or the connections will not be made correctly.
        """
        self.LeafRange = kwargs.get('LeafRange', None)
        self.Id = kwargs.get('Id', None)
        super(RangeNode, self).__init__(*args, **kwargs)

    def __int__(self):
        """Returns index of self."""
        return self.Index

    def _find_label(self):
        """Makes up a label for __str__ method."""
        label = None
        if hasattr(self, 'Name'):
            label = self.Name
        if label is None:
            if hasattr(self, 'Id'):
                label = self.Id
        if label is None:
            label = ''
        return label

    def __str__(self):
        """Returns informal Newick-like representation of self."""
        label = self._find_label()
        if self.Children:
            child_string = ','.join(map(str, self.Children))
        else:
            child_string = ''

        if self.Parent is None: #root of tree
            if self.Children:
                return '(%s)%s' % (child_string, label)
            else:
                return '()%s' % label
        else:   #internal node
            if self.Children:
                if hasattr(self, 'Length') and (self.Length!=None):
                    return '(%s)%s:%s' % \
                        (child_string, label, self.Length)
                else:
                    return '(%s)%s' % \
                        (child_string, label)
            else:
                if hasattr(self, 'Length') and (self.Length!=None):
                    return '%s:%s' % (label, self.Length)
                else:
                    return '%s' % (label)

    def traverse(self, self_before=False, self_after=False, include_self=True):
        """Iterates through children of self. Default behavior: leaves only.
        
        self_before: yield self before children (preorder traversal)

        self_after: yield self after children (postorder traversal)

        include_self: if False (default is True), skips self in traversal.
        Primarily useful for skipping the root node.

        If both self_before and self_after are True, the node is returned 
        both before _and_ after all its children are handled. This can be
        useful for certain applications, e.g. in RNA structure.
        """
        return super(RangeNode, self).traverse(self_before, self_after, \
            include_self)

    def indexByAttr(self, attr, multiple=False):
        """Returns dict of node.attr -> node.
        
        WARNING: Assumes all nodes have unique values of attr unless
        multiple is set to True.
        """
        result = {}
        if multiple:
            for n in self.traverse(self_before=True):
                curr = getattr(n, attr)
                if curr not in result:
                    result[curr] = [n]
                else:
                    result[curr].append(n)
        else:
            for n in self.traverse(self_before=True):
                result[getattr(n, attr)] = n
        return result

    def indexByFunc(self, f):
        """Returns dict of f(node) -> [matching nodes]."""
        result = {}
        for n in self.traverse(self_before=True):
            val = f(n)
            if val not in result:
                result[val] = [n]
            else:
                result[val].append(n)
        return result

    def assignIds(self, num_leaves=None):
        """Assigns each node's Id property, based on order in the tree.
        
        WARNING: Will store incorrect data if num_leaves is incorrect.
        """
        if num_leaves is None:
            num_leaves = len(list(self.traverse()))
        last_leaf_id = 0
        last_internal_id = num_leaves
        for node in self.traverse(self_after=True):
            c = node.Children
            if c:
                node.Id = last_internal_id
                last_internal_id += 1
                node.LeafRange = (c[0].LeafRange[0], c[-1].LeafRange[-1])
            else:
                node.Id = last_leaf_id
                node.LeafRange = (last_leaf_id, last_leaf_id + 1)
                last_leaf_id += 1

    def propagateAttr(self, attr, overwrite=False):
        """Propagates self's version of attr to all children without attr.
        
        overwrite: determines whether to overwrite existing attr values.
        """
        curr = getattr(self, attr)
        if overwrite:
            for node in self.traverse(self_after=True):
                setattr(node, attr, curr)
        else:
            for node in self.Children:
                if not hasattr(node, attr):
                    setattr(node, attr, curr)
                node.propagateAttr(attr)

    def delAttr(self, attr):
        """Delets attr in self and all children."""
        for node in self.traverse(self_after=True):
            delattr(node, attr)

    def perturbAttr(self, attr, f, pass_attr=False):
        """Perturbs attr in self and all children according to f() or f(attr).
        
        If pass_attr is False (the default), the branch is assigned f(). If
        pass_attr is True, the branch is assigned f(attr). Make sure that
        f has the correct form or you'll get an error!
        """
        for node in self.traverse(self_after=True):
            if pass_attr:
                setattr(node, attr, f(getattr(node, attr)))
            else:
                setattr(node, attr, f())

    def accumulateAttr(self, attr, towards_leaves=True, f=add):
        """Sets each node's version of attr to f(node, parent|children).
        
        if towards_leaves (the default), node.attr = f(node.attr, parent.attr);
        otherwise, node.attr = f(node.attr, child.attr) for each child.
        """
        if towards_leaves:
            for node in self.traverse(self_before=True):
                parent = node.Parent
                if parent is None:
                    continue
                else:
                    setattr(node, attr, \
                        f(getattr(node,attr), getattr(parent,attr)))
        else:
            for node in self.traverse(self_after=True):
                children = node.Children
                if children:
                    for c in children:
                        setattr(node, attr, \
                            f(getattr(node, attr), getattr(c, attr)))
   
    def accumulateChildAttr(self, attr, f=add):
        """Sets each node's attr based on states in children (only).
        
        Always works from leaves to root. Does not set states in leaves.
        Skips any child where attr is None.
        """
        for node in self.traverse(self_after=True):
            #only reset nodes with children
            if node.Children:
                #get attr from all children that have it
                vals = [getattr(c, attr) for c in node.Children \
                    if hasattr(c, attr)]
                #get rid of None values
                vals = filter(lambda x: x is not None, vals)
                if vals:
                    setattr(node, attr, reduce(f, vals))
                else:
                    setattr(node, attr, None)
                        
    def assignLevelsFromRoot(self):
        """Assigns each node its level reletave to self (self.Level=0)."""
        self.Level = 0
        self.propagateAttr('Level', overwrite=True)
        self.accumulateAttr('Level', towards_leaves=True, f=lambda a,b:b+1)

    def assignLevelsFromLeaves(self, use_min=False):
        """Assigns each node its distance from the leaves.

        use_min: use min distance from leaf instead of max (default:False)
        """
        self.Level = 0
        self.propagateAttr('Level', overwrite=True)
        if use_min:
            self.accumulateAttr('Level', towards_leaves=False, \
                f=lambda a,b: (a and min(a, b+1)) or b+1)
        else:
            self.accumulateAttr('Level', towards_leaves=False, \
                f=lambda a,b:max(a, b+1))
            
    def attrToList(self, attr, default=None, size=None, \
            leaves_only=False):
        """Copies attribute from each node of self into list.

        attr: name of attr to copy.
        size: size of list to copy into (must be >= num nodes).
        leaves_only: only look at leaves, not internal nodes

        WARNING: will fail if the Id attribute of each node has
        not yet been set.
        """
        if leaves_only:
            nodes = list(self.traverse())
        else:
            nodes = list(self.traverse(self_before=True))
        if size is None:
            size = len(nodes)
        result = [default] * size
        for node in nodes:
            result[node.Id] = getattr(node, attr)
        return result

    def attrFromList(self, attr, items, leaves_only=False):
        """Copies items in list into attr of nodes. Must have right # items."""
        for n in self.traverse(self_before = not leaves_only):
            setattr(n, attr, items[n.Id])

    def toBreakpoints(self):
        """Returns list of breakpoints that reconstructs self's topology.

        WARNING: Only works for strictly bifurcating trees.
        """
        result = []
        for node in self.traverse(self_before=True):
            if node.Children:
                result.append(node.Children[0].LeafRange[-1] - 1)
        return result

    def fromBreakpoints(cls, breakpoints):
        """Makes a new RangeNode tree from a sequence of breakpoints.

        Will have one more leaf than breakpoint. Always produces bifurcating
        tree.

        WARNING: will return incorrect results if elements in breakpoints are
        not unique!

        To make a random tree. call fromBreakpoints(permutation(n-1)) where
        n is the number of leaves desired in the tree.
        """
        #return single, leaf node if breakpoints is empty
        if not any(breakpoints):
            return cls(Id=0, LeafRange=(0,1))
        
        num_leaves = len(breakpoints) + 1
        curr_internal_index = num_leaves
        root = cls(Id=curr_internal_index, LeafRange=(0,num_leaves))
        curr_internal_index += 1
        
        #need to walk through the tree for each breakpoint, find the range
        #in which the breakpoint occurs, and make children containing the
        #start (i.e. start:breakpoint+1) and end (i.e. breakpoint+1:end)
        #of the range.
        for b in breakpoints:
            #start at the root
            curr_node = root
            children = curr_node.Children
            #walk down the tree until we find a range without children that
            #the breakpoint is in
            while children:
                middle = children[1].LeafRange[0]
                #SUPPORT2425
                curr_node = children[int(middle <= b)]
                children = curr_node.Children
            #curr_node is now the range that contains the breakpoint
            start, end = curr_node.LeafRange
            #check if left and right nodes are leaves, and assign relevant ids
            #we frequently need the index after the breakpoint, so assign
            #variable after_b to avoid lots of mysterious 'b+1's in the code
            after_b = b+1
            if after_b - start == 1:
                left_id = start
            else:
                left_id = curr_internal_index
                curr_internal_index += 1
            if end - after_b == 1:
                right_id = after_b
            else:
                right_id = curr_internal_index
                curr_internal_index += 1
            #add left and right nodes to current node's children
            left = cls(Parent=curr_node,Id=left_id, LeafRange=(start, after_b))
            right = cls(Parent=curr_node,Id=right_id, LeafRange=(after_b, end))
            curr_node.Children = [left, right]
        return root
            
    fromBreakpoints = classmethod(fromBreakpoints)

    def leafLcaDepths(self, assign_ids=True, assign_levels=True):
        """Returns num_leaves x num_leaves matrix with depth of each LCA.
        
        assign_ids and assign_levels control whether or not to assign
        ids and levels (default: True).

        size: if supplied, sizes the matrix.

        No longer assumes strictly bifurcating tree.
        """
        if assign_ids:
            self.assignIds()
        if assign_levels:
            self.assignLevelsFromLeaves()
        nodes = list(self.traverse(self_before=True))
        #second element of LeafRange should contain largest node index
        #incidentally, will fail if ids not assigned
        num_nodes = self.LeafRange[1]
        result = zeros((num_nodes,num_nodes))
        for node in nodes:
            #skip any nodes that are themselves leaves
            children = node.Children
            if not children:
                continue
            #if node has only one child, can't be anyone's LCA
            if len(children) == 1:
                continue
            if len(children) == 2:
            #if node has two children, is LCA of any descendant of first
            #child w.r.t. any descendant of second child
                curr_level = node.Level
                left, right = children
                for left_index in range(*(left.LeafRange)):
                    for right_index in range(*(right.LeafRange)):
                        result[left_index, right_index] = curr_level
            #otherwise, node is LCA of each child's descendants w.r.t. the
            #descendants of other children
            else:
                curr_level = node.Level
                for first in children:
                    for second in children[1:]:
                        for left_index in range(*(first.LeafRange)):
                            for right_index in range(*(second.LeafRange)):
                                result[left_index, right_index] = curr_level
        result += transpose(result)
        return result

    def randomNode(self):
        """Returns random node from self and children."""
        return choice(list(self.traverse(self_before=True)))
    
    def randomLeaf(self):
        """Returns random leaf descended from self."""
        if self.Children:
            return choice(list(self.traverse()))
        else:
            return self

    def randomNodeWithNLeaves(self, n):
        """Returns random node with exactly the specified number of leaves."""
        try:
            lookup = self.indexByFunc(lambda x: x.LeafRange[1] - x.LeafRange[0])
        except TypeError: #possible that ranges weren't assigned
            self.assignIds()
            lookup = self.indexByFunc(lambda x: x.LeafRange[1] - x.LeafRange[0])
        return choice(lookup[n])

    def randomNodeAtLevel(self, n, from_leaves=True):
        """Returns random node at specified level from root or tips."""
        if from_leaves:
            self.assignLevelsFromLeaves()
        else:
            self.assignLevelsFromRoot()
        lookup = self.indexByAttr('Level', multiple=True)
        return choice(lookup[n])

    def outgroupLast(self, first, second, third, cache=True):
        """Returns tuple of nodes first, second and third, with outgroup last.

        first, second, and third must all be descendants of self, and ids
        must have already been assigned to the trees.

        Sets self._leaf_lca_depths if not already set if cache is True.

        WARNING: if first, second, third are all at the same level of an
        unresolved polytomy, will arbitrarily choose one of the three as
        an outgroup. This choice may be inconsistent between different
        runs of the program.
        """
        #find the leaf lca depths if necessary
        if cache:
            if not hasattr(self, '_leaf_lca_depths'):
                self._leaf_lca_depths = self.leafLcaDepths()
            depths = self._leaf_lca_depths
        else:
            depths = self.leafLcaDepths()
        #get the ids of the nodes
        first_id, second_id, third_id = first.Id, second.Id, third.Id
        lca_12 = depths[first_id, second_id]
        lca_13 = depths[first_id, third_id]
        lca_23 = depths[second_id, third_id]
        #find the shallowest lca and return nodes in appropriate order
        shallowest_lca = min([lca_12, lca_13, lca_23])
        if shallowest_lca == lca_12:
            return (first, second, third)
        elif shallowest_lca == lca_13:
            return first, third, second
        else:
            return second, third, first

    def filter(self, taxa, keep=True):
        """Prunes (inplace) all items in self not leading to a taxon in taxa.

        Taxa must be container of nodes in the tree.

        keep determines whether to keep (if True) or delete (if False) the
        specified taxa.

        Collapses nodes where appropriate (i.e. one-child nodes get deleted).
        Branch lengths are preserved (i.e. if a node is collapsed, its branch
        length is added to the node that collapses onto it).

        WARNING: Root of the tree is always preserved, so you might find that
        all the nodes are in a single-child subtree of the root if all the nodes
        on the other side of the root were deleted. If no nodes are kept, 
        will return an empty root node with no children.
        """
        taxon_ids = dict.fromkeys(map(id, taxa))
        node_ids = self.indexByFunc(id)
        #select specified ids
        for t in taxon_ids:
            if t in node_ids:
                node_ids[t][0]._selected = True
        #unselect root if not specified
        if not id(self) in taxon_ids:
            self._selected = False
        #figure out whether each node is selected: first, accumulate towards
        #tips, then trace back to root
        self.propagateAttr('_selected')
        if keep:
            self.accumulateChildAttr('_selected', f=or_)
        else:
            self.accumulateChildAttr('_selected', f=and_)
        #delete and/or collapse undesired nodes
        for node in list(self.traverse(self_after=True)):
            if node.Parent is None: #back at root
                continue
            #delete if not selected
            if node._selected != keep:
                result = node.Parent.removeNode(node)
                node.Parent = None
            #replace with (already handled) child if single-item
            elif (len(node.Children) == 1) and (node.Parent is not None):
                curr_child = node.Children[0]
                curr_parent = node.Parent
                curr_parent.Children[curr_parent.Children.index(node)] = \
                    curr_child
                curr_child.Parent = curr_parent
                node.Parent = None
                #add branch lengths if present
                if hasattr(node, 'Length'):
                    if hasattr(curr_child, 'Length') and \
                        curr_child.Length is not None:
                        curr_child.Length += node.Length
                    else:
                        curr_child.Length = node.Length
        self.delAttr('_selected')

    def addChildren(self, n):
        """Adds n children to self."""
        constructor = self.__class__
        new_nodes = [constructor(Parent=self) for i in range(n)]

    def makeIdIndex(self):
        """Sets self.IdIndex as index of ids in self."""
        self.assignIds()
        self.IdIndex = self.indexByAttr('Id')
    
    def assignQ(self, q=None, special_qs=None, overwrite=False):
        """Clears and assigns Q matrices.

        q: overall Q for tree.
        special_qs: dict of node id -> q for node and its subtrees.
        overwrite: if True (default is False), overwrites existing Qs
        """
        if q is not None:
            self.Q = q
        if special_qs:
            ids = self.IdIndex
            for k, v in special_qs.items():
                ids[int(k)].Q = v
        if (not hasattr(self, 'Q')) or (not self.Q):
            raise ValueError, "Failed to assign Q matrix to root."
        self.propagateAttr('Q', overwrite)

    def assignP(self):
        """Assigns P-matrices based on current Q-matrices.
        
        Assumes that Length and Q are already set in all nodes.

        WARNING: Assumes that branch lengths represent sequence divergences and
        are at most 0.75 (for DNA, somewhat more for protein), i.e. do not
        exceed saturation. If the branch lengths exceed saturation, will
        probably fail unpredictably. Note that Q.toSimilarProbs generates a
        sequence that matches a specific similarity, not a specific divergence,
        so need to use (1-Length) (with the attendant dangers).

        WARNING: Does not set P-matrix at root (because there's no change
        before the root. Should it instead set the P-matrix at the root to
        the identity matrix of the appropriate size?)
        """
        #don't set P if root
        if self.Parent is not None:
            self.P = self.Q.toSimilarProbs(1-self.Length)
        for c in self.Children:
            c.assignP()

    def assignLength(self, length):
        """Assigns all nodes the specified length."""
        for node in self.traverse(self_before=True):
            node.Length = length

    def evolve(self, seq, field_name='Sequence'):
        """Evolves seq, according to P-matrix on each node.
        
        Assumes that P has been set on all nodes already.

        WARNING: seq is an array, not a Sequence object (at this point).
        """
        #assign seq to root, or evolve from parent's sequence
        if self.Parent is None:
            setattr(self, field_name, seq)
        else:
            setattr(self, field_name, self.P.mutate(seq))
        for c in self.Children:
            c.evolve(getattr(self, field_name), field_name)

    def assignPs(self, rates):
        """Sets many P-matrices from a single Q-matrix, scaled to rates.
        
        Assumes that Branchlength and Q are already set in all nodes.

        rates should be a list of rates at least as long as the seqs to evolve.
        These should all be less than 1 (i.e. the max rate is 1, and other rates
        decline from there).
        """
        if self.Parent is not None:
            self.Ps = [self.Q.toSimilarProbs(1-(self.Length*r)) \
                for r in rates]
        for c in self.Children:
            c.assignPs(rates)

    def evolveSeqs(self, seqs, field_name='Sequences'):
        """Evolves list of seqs according to Q-matrices and rates on each node.

        Assume that Ps has already been set such that P_i -> seq_i at each node,
        e.g. via self.assignPs.

        WARNING: seqs are (currently) assumed to be arrays, not Sequence
        objects.
        """
        if self.Parent is None:
            setattr(self, field_name, seqs)
        else:
            setattr(self, field_name, [p.mutate(seq) for (p, seq) in \
                zip(self.Ps, seqs)])
        for c in self.Children:
            c.evolveSeqs(getattr(self, field_name), field_name)


def balanced_breakpoints(num_leaves):
    """Returns breakpoints for a balanced tree with specified num_leaves.
    
    num_leaves must be at least 1 and must be a power of 2. WARNING: no 
    validation is performed to ensure these conditions are met.

    This algorithm works by figuring the indices of all the nodes that are at
    a particular level, making an array of those indices using arange, and
    then concatenating the arrays in order of level.
    """
    result = []
    curr_step = num_leaves
    curr_start = (curr_step/2) - 1
    while curr_start >= 0:
        result.append(arange(curr_start, num_leaves, curr_step))
        curr_step /= 2
        curr_start = (curr_step/2) -1
    return concatenate(result)

def BalancedTree(num_leaves, node_class=RangeNode):
    """Returns a balanced tree of node_class (num_leaves must be power of 2)."""
    #return node_class.fromBreakpoints(balanced_breakpoints(num_leaves))
    root = node_class()
    curr_children = [root]

    while len(curr_children) < num_leaves:
        tmp = []
        for n in curr_children:
            n.Children[:] = [node_class(Parent=n), node_class(Parent=n)]
            tmp.extend(n.Children)
        curr_children = tmp

    return root

def RandomTree(num_leaves, node_class=RangeNode):
    """Returns a random node_class tree using the breakpoint model."""
    return node_class.fromBreakpoints(permutation(num_leaves-1))

def CombTree(num_leaves, deepest_first=True, node_class=RangeNode):
    """Returns a comb node_class tree."""
    if deepest_first:
        branch_child = 1
    else:
        branch_child = 0

    root = node_class()
    curr = root

    for i in range(num_leaves-1):
        curr.Children[:] = [node_class(Parent=curr),node_class(Parent=curr)]
        curr = curr.Children[branch_child]

    return root

def StarTree(num_leaves, node_class=RangeNode):
    """Returns a star phylogeny, with all leaves equally connected to root."""
    t = node_class()
    t.addChildren(num_leaves)
    return t

def LineTree(depth, node_class=RangeNode):
    """Returns a tree with all nodes arranged in a line."""
    t = node_class()
    curr = t
    for i in range(depth-1):
        new_node = node_class(Parent=curr)
        curr = new_node
    return t
