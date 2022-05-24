"""Classes for storing and manipulating a phylogenetic tree.

These trees can be either strictly binary, or have polytomies
(multiple children to a parent node).

Trees consist of Nodes (or branches) that connect two nodes. The Tree can
be created only from a newick formatted string read either from file or from a
string object. Other formats will be added as time permits.

Tree can:
    -  Deal with either rooted or unrooted tree's and can
       convert between these types.
    -  Return a sub-tree given a list of tip-names
    -  Identify an edge given two tip names. This method facilitates the
       statistical modelling by simplyifying the syntax for specifying
       sub-regions of a tree.
    -  Assess whether two Tree instances represent the same topology.

Definition of relevant terms or abbreviations:
    -  edge: also known as a branch on a tree.
    -  node: the point at which two edges meet
    -  tip: a sequence or species
    -  clade: all and only the nodes (including tips) that descend
       from a node
    -  stem: the edge immediately preceeding a clade
"""
import json
import numbers
import re

from copy import deepcopy
from functools import reduce
from itertools import combinations
from operator import or_
from random import choice, shuffle

from numpy import argsort, ceil, log, zeros

from cogent3.maths.stats.test import correlation
from cogent3.util.io import atomic_write, get_format_suffixes
from cogent3.util.misc import get_object_provenance


__author__ = "Gavin Huttley, Peter Maxwell and Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Gavin Huttley",
    "Peter Maxwell",
    "Rob Knight",
    "Andrew Butterfield",
    "Catherine Lozupone",
    "Micah Hamady",
    "Jeremy Widmann",
    "Zongzhi Liu",
    "Daniel McDonald",
    "Justin Kuczynski",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def distance_from_r_squared(m1, m2):
    """Estimates distance as 1-r^2: no correl = max distance"""
    return 1 - (correlation(m1.flat, m2.flat)[0]) ** 2


def distance_from_r(m1, m2):
    """Estimates distance as (1-r)/2: neg correl = max distance"""
    return (1 - correlation(m1.flat, m2.flat)[0]) / 2


class TreeError(Exception):
    pass


class TreeNode(object):
    """Store information about a tree node. Mutable.

    Parameters:
        name: label for the node, assumed to be unique.
        children: list of the node's children.
        params: dict containing arbitrary parameters for the node.
        name_loaded: ?
    """

    _exclude_from_copy = dict.fromkeys(["_parent", "children"])

    def __init__(
        self,
        name=None,
        children=None,
        parent=None,
        params=None,
        name_loaded=True,
        **kwargs,
    ):
        """Returns new TreeNode object."""
        self.name = name
        self.name_loaded = name_loaded
        if params is None:
            params = {}
        self.params = params
        self.children = []
        if children is not None:
            self.extend(children)
        self._parent = parent
        if (parent is not None) and not (self in parent.children):
            parent.append(self)

    # built-in methods and list interface support
    def __repr__(self):
        """Returns reconstructable string representation of tree.

        WARNING: Does not currently set the class to the right type.
        """
        return f'Tree("{self.get_newick()}")'

    def __str__(self):
        """Returns Newick-format string representation of tree."""
        return self.get_newick()

    # TODO have methods that need to rely on identity of self and
    # other actually do that
    # For now, the following comparison operators are peculiar in that
    # that by omitting eq/ne methods those default to id()
    # whereas sorting should be based on name
    # I think the remove etc .. operations should explicitly
    # used id()
    # def __eq__(self, other):
    # return self.name == other.name

    # def __ne__(self, other):
    # return self.name != other.name

    def __lt__(self, other):
        return self.name < other.name

    def __gt__(self, other):
        return self.name > other.name

    def compare_name(self, other):
        """Compares TreeNode by name"""
        if self is other:
            return True

        return self.name == other.name

    def compare_by_names(self, other):
        """Equality test for trees by name"""
        # if they are the same object then they must be the same tree...
        if self is other:
            return True
        self_names = self.get_node_names()
        other_names = other.get_node_names()
        self_names = [v for v in self_names if v is not None]
        self_names.sort()
        other_names = [v for v in other_names if v is not None]
        other_names.sort()
        return self_names == other_names

    def _to_self_child(self, i):
        """Converts i to self's type, with self as its parent.

        Cleans up refs from i's original parent, but doesn't give self ref to i.
        """
        c = self.__class__
        if isinstance(i, c):
            if i._parent not in (None, self):
                i._parent.children.remove(i)
        else:
            i = c(i)
        i._parent = self
        return i

    def append(self, i):
        """Appends i to self.children, in-place, cleaning up refs."""
        self.children.append(self._to_self_child(i))

    def extend(self, items):
        """Extends self.children by items, in-place, cleaning up refs."""
        self.children.extend(list(map(self._to_self_child, items)))

    def insert(self, index, i):
        """Inserts an item at specified position in self.children."""
        self.children.insert(index, self._to_self_child(i))

    def pop(self, index=-1):
        """Returns and deletes child of self at index (default: -1)"""
        result = self.children.pop(index)
        result._parent = None
        return result

    def remove(self, target):
        """Removes node by name instead of identity.

        Returns True if node was present, False otherwise.
        """
        if isinstance(target, TreeNode):
            target = target.name
        for (i, curr_node) in enumerate(self.children):
            if curr_node.name == target:
                self.remove_node(curr_node)
                return True
        return False

    def __getitem__(self, i):
        """Node delegates slicing to children; faster to access them
        directly."""
        return self.children[i]

    def __setitem__(self, i, val):
        """Node[i] = x sets the corresponding item in children."""
        curr = self.children[i]
        if isinstance(i, slice):
            for c in curr:
                c._parent = None
            coerced_val = list(map(self._to_self_child, val))
            self.children[i] = coerced_val[:]
        else:  # assume we got a single index
            curr._parent = None
            coerced_val = self._to_self_child(val)
            self.children[i] = coerced_val

    def __delitem__(self, i):
        """del node[i] deletes index or slice from self.children."""
        curr = self.children[i]
        if isinstance(i, slice):
            for c in curr:
                c._parent = None
        else:
            curr._parent = None
        del self.children[i]

    def __iter__(self):
        """Node iter iterates over the children."""
        return iter(self.children)

    def __len__(self):
        """Node len returns number of children."""
        return len(self.children)

    # support for copy module
    def copy_recursive(self, memo=None, _nil=None, constructor="ignored"):
        """Returns copy of self's structure, including shallow copy of attrs.

        constructor is ignored; required to support old tree unit tests.
        """
        _nil = _nil or []
        result = self.__class__()
        efc = self._exclude_from_copy
        for k, v in list(self.__dict__.items()):
            if k not in efc:  # avoid infinite recursion
                result.__dict__[k] = deepcopy(self.__dict__[k])
        for c in self:
            result.append(c.copy())
        return result

    def copy(self, memo=None, _nil=None, constructor="ignored"):
        """Returns a copy of self using an iterative approach"""

        _nil = _nil or []

        def __copy_node(n):
            result = n.__class__()
            efc = n._exclude_from_copy
            for k, v in list(n.__dict__.items()):
                if k not in efc:
                    result.__dict__[k] = deepcopy(n.__dict__[k])
            return result

        root = __copy_node(self)
        nodes_stack = [[root, self, len(self.children)]]

        while nodes_stack:
            # check the top node, any children left unvisited?
            top = nodes_stack[-1]
            new_top_node, old_top_node, unvisited_children = top

            if unvisited_children:
                top[2] -= 1
                old_child = old_top_node.children[-unvisited_children]
                new_child = __copy_node(old_child)
                new_top_node.append(new_child)
                nodes_stack.append([new_child, old_child, len(old_child.children)])
            else:  # no unvisited children
                nodes_stack.pop()
        return root

    __deepcopy__ = deepcopy = copy

    def copy_topology(self, constructor=None):
        """Copies only the topology and labels of a tree, not any extra data.

        Useful when you want another copy of the tree with the same structure
        and labels, but want to e.g. assign different branch lengths and
        environments. Does not use deepcopy from the copy module, so _much_
        faster than the copy() method.
        """
        if constructor is None:
            constructor = self.__class__
        children = [c.copy_topology(constructor) for c in self.children]
        return constructor(name=self.name[:], children=children)

    # support for basic tree operations -- finding objects and moving in the
    # tree
    def _get_parent(self):
        """Accessor for parent.

        If using an algorithm that accesses parent a lot, it will be much
        faster to access self._parent directly, but don't do it if mutating
        self._parent! (or, if you must, remember to clean up the refs).
        """
        return self._parent

    def _set_parent(self, parent):
        """Mutator for parent: cleans up refs in old parent."""
        if self._parent is not None:
            self._parent.remove_node(self)
        self._parent = parent
        if (parent is not None) and (self not in parent.children):
            parent.children.append(self)

    parent = property(_get_parent, _set_parent)

    def index_in_parent(self):
        """Returns index of self in parent."""
        return self._parent.children.index(self)

    def is_tip(self):
        """Returns True if the current node is a tip, i.e. has no children."""
        return not self.children

    def is_root(self):
        """Returns True if the current is a root, i.e. has no parent."""
        return self._parent is None

    def traverse(self, self_before=True, self_after=False, include_self=True):
        """Returns iterator over descendants. Iterative: safe for large trees.

        self_before includes each node before its descendants if True.
        self_after includes each node after its descendants if True.
        include_self includes the initial node if True.

        self_before and self_after are independent. If neither is True, only
        terminal nodes will be returned.

        Note that if self is terminal, it will only be included once even if
        self_before and self_after are both True.

        This is a depth-first traversal. Since the trees are not binary,
        preorder and postorder traversals are possible, but inorder traversals
        would depend on the data in the tree and are not handled here.
        """
        if self_before:
            if self_after:
                return self.pre_and_postorder(include_self=include_self)
            else:
                return self.preorder(include_self=include_self)
        else:
            if self_after:
                return self.postorder(include_self=include_self)
            else:
                return self.tips(include_self=include_self)

    def levelorder(self, include_self=True):
        """Performs levelorder iteration over tree"""
        queue = [self]
        while queue:
            curr = queue.pop(0)
            if include_self or (curr is not self):
                yield curr
            if curr.children:
                queue.extend(curr.children)

    def preorder(self, include_self=True):
        """Performs preorder iteration over tree."""
        stack = [self]
        while stack:
            curr = stack.pop()
            if include_self or (curr is not self):
                yield curr
            if curr.children:
                stack.extend(curr.children[::-1])  # 20% faster than reversed

    def postorder(self, include_self=True):
        """performs postorder iteration over tree.

        Notes
        -----

        This is somewhat inelegant compared to saving the node and its index
        on the stack, but is 30% faster in the average case and 3x faster in
        the worst case (for a comb tree).
        """
        child_index_stack = [0]
        curr = self
        curr_children = self.children
        curr_children_len = len(curr_children)
        while 1:
            curr_index = child_index_stack[-1]
            # if there are children left, process them
            if curr_index < curr_children_len:
                curr_child = curr_children[curr_index]
                # if the current child has children, go there
                if curr_child.children:
                    child_index_stack.append(0)
                    curr = curr_child
                    curr_children = curr.children
                    curr_children_len = len(curr_children)
                    curr_index = 0
                # otherwise, yield that child
                else:
                    yield curr_child
                    child_index_stack[-1] += 1
            # if there are no children left, return self, and move to
            # self's parent
            else:
                if include_self or (curr is not self):
                    yield curr
                if curr is self:
                    break
                curr = curr.parent
                curr_children = curr.children
                curr_children_len = len(curr_children)
                child_index_stack.pop()
                child_index_stack[-1] += 1

    def pre_and_postorder(self, include_self=True):
        """Performs iteration over tree, visiting node before and after."""
        # handle simple case first
        if not self.children:
            if include_self:
                yield self
            return
        child_index_stack = [0]
        curr = self
        curr_children = self.children
        while 1:
            curr_index = child_index_stack[-1]
            if not curr_index:
                if include_self or (curr is not self):
                    yield curr
            # if there are children left, process them
            if curr_index < len(curr_children):
                curr_child = curr_children[curr_index]
                # if the current child has children, go there
                if curr_child.children:
                    child_index_stack.append(0)
                    curr = curr_child
                    curr_children = curr.children
                    curr_index = 0
                # otherwise, yield that child
                else:
                    yield curr_child
                    child_index_stack[-1] += 1
            # if there are no children left, return self, and move to
            # self's parent
            else:
                if include_self or (curr is not self):
                    yield curr
                if curr is self:
                    break
                curr = curr.parent
                curr_children = curr.children
                child_index_stack.pop()
                child_index_stack[-1] += 1

    def traverse_recursive(self, self_before=True, self_after=False, include_self=True):
        """Returns iterator over descendants. IMPORTANT: read notes below.

        traverse_recursive is slower than traverse, and can lead to stack
        errors. However, you _must_ use traverse_recursive if you plan to
        modify the tree topology as you walk over it (e.g. in post-order),
        because the iterative methods use their own stack that is not updated
        if you alter the tree.

        self_before includes each node before its descendants if True.
        self_after includes each node after its descendants if True.
        include_self includes the initial node if True.

        self_before and self_after are independent. If neither is True, only
        terminal nodes will be returned.

        Note that if self is terminal, it will only be included once even if
        self_before and self_after are both True.

        This is a depth-first traversal. Since the trees are not binary,
        preorder and postorder traversals are possible, but inorder traversals
        would depend on the data in the tree and are not handled here.
        """
        if self.children:
            if self_before and include_self:
                yield self
            for child in self.children:
                for i in child.traverse_recursive(self_before, self_after):
                    yield i
            if self_after and include_self:
                yield self
        elif include_self:
            yield self

    def ancestors(self):
        """Returns all ancestors back to the root. Dynamically calculated."""
        result = []
        curr = self._parent
        while curr is not None:
            result.append(curr)
            curr = curr._parent
        return result

    def root(self):
        """Returns root of the tree self is in. Dynamically calculated."""
        curr = self
        while curr._parent is not None:
            curr = curr._parent
        return curr

    def isroot(self):
        """Returns True if root of a tree, i.e. no parent."""
        return self.is_root()

    def siblings(self):
        """Returns all nodes that are children of the same parent as self.

        Note: excludes self from the list. Dynamically calculated.
        """
        if self._parent is None:
            return []
        result = self._parent.children[:]
        result.remove(self)
        return result

    def iter_tips(self, include_self=False):
        """Iterates over tips descended from self, [] if self is a tip."""
        # bail out in easy case
        if not self.children:
            if include_self:
                yield self
            return None
        # use stack-based method: robust to large trees
        stack = [self]
        while stack:
            curr = stack.pop()
            if curr.children:
                stack.extend(curr.children[::-1])  # 20% faster than reversed
            else:
                yield curr

    def tips(self, include_self=False):
        """Returns tips descended from self, [] if self is a tip."""
        return list(self.iter_tips(include_self=include_self))

    def iter_nontips(self, include_self=False):
        """Iterates over nontips descended from self, [] if none.

        include_self, if True (default is False), will return the current
        node as part of the list of nontips if it is a nontip."""
        for n in self.traverse(True, False, include_self):
            if n.children:
                yield n

    def nontips(self, include_self=False):
        """Returns nontips descended from self."""
        return list(self.iter_nontips(include_self=include_self))

    def istip(self):
        """Returns True if is tip, i.e. no children."""
        return not self.children

    def tip_children(self):
        """Returns direct children of self that are tips."""
        return [i for i in self.children if not i.children]

    def non_tip_children(self):
        """Returns direct children in self that have descendants."""
        return [i for i in self.children if i.children]

    def child_groups(self):
        """Returns list containing lists of children sharing a state.

        In other words, returns runs of tip and nontip children.
        """
        # bail out in trivial cases of 0 or 1 item
        if not self.children:
            return []
        if len(self.children) == 1:
            return [self.children[0]]
        # otherwise, have to do it properly...
        result = []
        curr = []
        state = None
        for i in self.children:
            curr_state = bool(i.children)
            if curr_state == state:
                curr.append(i)
            else:
                if curr:
                    result.append(curr)
                    curr = []
                curr.append(i)
                state = curr_state
        # handle last group
        result.append(curr)
        return result

    def last_common_ancestor(self, other):
        """Finds last common ancestor of self and other, or None.

        Always tests by identity.
        """
        my_lineage = set([id(node) for node in [self] + self.ancestors()])
        curr = other
        while curr is not None:
            if id(curr) in my_lineage:
                return curr
            curr = curr._parent
        return None

    def lowest_common_ancestor(self, tipnames):
        """Lowest common ancestor for a list of tipnames

        This should be around O(H sqrt(n)), where H is height and n is the
        number of tips passed in.
        """
        if len(tipnames) == 1:
            return self.get_node_matching_name(tipnames[0])

        tipnames = set(tipnames)
        tips = [tip for tip in self.tips() if tip.name in tipnames]

        if len(tips) != len(tipnames):
            missing = tipnames - set(self.get_tip_names())
            raise ValueError(f"tipnames {missing} not present in self")

        # scrub tree
        if hasattr(self, "black"):
            for n in self.traverse(include_self=True):
                if hasattr(n, "black"):
                    delattr(n, "black")

        for t in tips:
            prev = t
            curr = t.parent

            while curr and not hasattr(curr, "black"):
                setattr(curr, "black", [prev])
                prev = curr
                curr = curr.parent

            # increase black count, multiple children lead to here
            if curr:
                curr.black.append(prev)

        curr = self
        while len(curr.black) == 1:
            curr = curr.black[0]

        return curr

    lca = last_common_ancestor  # for convenience

    # support for more advanced tree operations

    def separation(self, other):
        """Returns number of edges separating self and other."""
        # detect trivial case
        if self is other:
            return 0
        # otherwise, check the list of ancestors
        my_ancestors = dict.fromkeys(list(map(id, [self] + self.ancestors())))
        count = 0
        while other is not None:
            if id(other) in my_ancestors:
                # need to figure out how many steps there were back from self
                curr = self
                while not (curr is None or curr is other):
                    count += 1
                    curr = curr._parent
                return count
            else:
                count += 1
                other = other._parent
        return None

    def descendant_array(self, tip_list=None):
        """Returns numpy array with nodes in rows and descendants in columns.

        A value of 1 indicates that the decendant is a descendant of that node/
        A value of 0 indicates that it is not

        Also returns a list of nodes in the same order as they are listed
        in the array.

        tip_list is a list of the names of the tips that will be considered,
        in the order they will appear as columns in the final array. Internal
        nodes will appear as rows in preorder traversal order.
        """

        # get a list of internal nodes
        node_list = [node for node in self.traverse() if node.children]
        node_list.sort()

        # get a list of tip names if one is not supplied
        if not tip_list:
            tip_list = [n.name for n in self.tips()]
            tip_list.sort()
        # make a blank array of the right dimensions to alter
        result = zeros([len(node_list), len(tip_list)])
        # put 1 in the column for each child of each node
        for (i, node) in enumerate(node_list):
            children = [n.name for n in node.tips()]
            for (j, dec) in enumerate(tip_list):
                if dec in children:
                    result[i, j] = 1
        return result, node_list

    def _default_tree_constructor(self):
        return TreeBuilder(constructor=self.__class__).edge_from_edge

    def name_unnamed_nodes(self):
        """sets the Data property of unnamed nodes to an arbitrary value

        Internal nodes are often unnamed and so this function assigns a
        value for referencing."""
        # make a list of the names that are already in the tree
        names_in_use = []
        for node in self.traverse():
            if node.name:
                names_in_use.append(node.name)
        # assign unique names to the Data property of nodes where Data = None
        name_index = 1
        for node in self.traverse():
            if not node.name:
                new_name = "node" + str(name_index)
                # choose a new name if name is already in tree
                while new_name in names_in_use:
                    name_index += 1
                    new_name = "node" + str(name_index)
                node.name = new_name
                names_in_use.append(new_name)
                name_index += 1

    def make_tree_array(self, dec_list=None):
        """Makes an array with nodes in rows and descendants in columns.

        A value of 1 indicates that the decendant is a descendant of that node/
        A value of 0 indicates that it is not

        also returns a list of nodes in the same order as they are listed
        in the array"""
        # get a list of internal nodes
        node_list = [node for node in self.traverse() if node.children]
        node_list.sort()

        # get a list of tips() name if one is not supplied
        if not dec_list:
            dec_list = [dec.name for dec in self.tips()]
            dec_list.sort()
        # make a blank array of the right dimensions to alter
        result = zeros((len(node_list), len(dec_list)))
        # put 1 in the column for each child of each node
        for i, node in enumerate(node_list):
            children = [dec.name for dec in node.tips()]
            for j, dec in enumerate(dec_list):
                if dec in children:
                    result[i, j] = 1
        return result, node_list

    def remove_deleted(self, is_deleted):
        """Removes all nodes where is_deleted tests true.

        Internal nodes that have no children as a result of removing deleted
        are also removed.
        """
        # Traverse tree
        for node in list(self.traverse(self_before=False, self_after=True)):
            # if node is deleted
            if is_deleted(node):
                # Store current parent
                curr_parent = node.parent
                # Set current node's parent to None (this deletes node)
                node.parent = None
                # While there are no chilren at node and not at root
                while (curr_parent is not None) and (not curr_parent.children):
                    # Save old parent
                    old_parent = curr_parent
                    # Get new parent
                    curr_parent = curr_parent.parent
                    # remove old node from tree
                    old_parent.parent = None

    def prune(self):
        """Reconstructs correct topology after nodes have been removed.

        Internal nodes with only one child will be removed and new connections
        will be made to reflect change.
        """
        # traverse tree to decide nodes to be removed.
        nodes_to_remove = []
        for node in self.traverse():
            if (node.parent is not None) and (len(node.children) == 1):
                nodes_to_remove.append(node)
        for node in nodes_to_remove:
            # save current parent
            curr_parent = node.parent
            # save child
            child = node.children[0]
            # remove current node by setting parent to None
            node.parent = None
            # Connect child to current node's parent
            child.parent = curr_parent

    def same_shape(self, other):
        """Ignores lengths and order, so trees should be sorted first"""
        if len(self.children) != len(other.children):
            return False
        if self.children:
            for (self_child, other_child) in zip(self.children, other.children):
                if not self_child.same_shape(other_child):
                    return False
            return True
        else:
            return self.name == other.name

    def to_rich_dict(self):
        """returns {'newick': with node names,
        'edge_attributes': {'tip1': {'length': ...}, ...}}"""
        newick = self.get_newick(
            with_node_names=True, semicolon=False, escape_name=False
        )
        attr = {}
        for edge in self.get_edge_vector(include_root=True):
            attr[edge.name] = edge.params.copy()
        result = dict(
            newick=newick,
            edge_attributes=attr,
            type=get_object_provenance(self),
            version=__version__,
        )
        return result

    def to_json(self):
        """returns json formatted string {'newick': with edges and distances, 'edge_attributes': }"""
        return json.dumps(self.to_rich_dict())

    def get_newick_recursive(
        self,
        with_distances=False,
        semicolon=True,
        escape_name=True,
        with_node_names=False,
    ):
        """Return the newick string for this edge.

        Parameters
        ----------
        with_distances
            whether branch lengths are included.
        semicolon
            end tree string with a semicolon
        escape_name
            if any of these characters []'"(),
            nodes name, wrap the name in single quotes
        with_node_names
            includes internal node names (except 'root')

        """
        newick = []
        subtrees = []
        for child in self.children:
            nwk = child.get_newick_recursive(
                with_distances, semicolon=False, with_node_names=with_node_names
            )
            subtrees.append(nwk)

        if subtrees:
            newick.append(f"({','.join(subtrees)})")

        if self.name_loaded or with_node_names:
            if self.name is None or with_node_names and self.is_root():
                name = ""
            else:
                name = str(self.name)
                if escape_name and not (name.startswith("'") and name.endswith("'")):
                    if re.search("""[]['"(),:;_]""", name):
                        name = "'%s'" % name.replace("'", "''")
                    else:
                        name = name.replace(" ", "_")
            newick.append(name)

        if isinstance(self, PhyloNode):
            if with_distances and self.length is not None:
                newick.append(f":{self.length}")

        if semicolon:
            newick.append(";")

        return "".join(newick)

    def get_newick(
        self,
        with_distances=False,
        semicolon=True,
        escape_name=True,
        with_node_names=False,
    ):
        """Return the newick string for this tree.

        Parameters
        ----------
        with_distances
            whether branch lengths are included.
        semicolon
            end tree string with a semicolon
        escape_name
            if any of these characters []'"(),
            nodes name, wrap the name in single quotes
        with_node_names
            includes internal node names (except 'root')

        NOTE: This method returns the Newick representation of this node
        and its descendents.
        """
        result = ["("]
        nodes_stack = [[self, len(self.children)]]
        node_count = 1

        while nodes_stack:
            node_count += 1
            # check the top node, any children left unvisited?
            top = nodes_stack[-1]
            top_node, num_unvisited_children = top
            if num_unvisited_children:  # has any child unvisited
                top[1] -= 1  # decrease the # of children unvisited
                # - for order
                next_child = top_node.children[-num_unvisited_children]
                # pre-visit
                if next_child.children:
                    result.append("(")
                nodes_stack.append([next_child, len(next_child.children)])
            else:  # no unvisited children
                nodes_stack.pop()
                # post-visit
                if top_node.children:
                    result[-1] = ")"

                if top_node.name_loaded or with_node_names:
                    if top_node.name is None or with_node_names and top_node.is_root():
                        name = ""
                    else:
                        name = str(top_node.name)
                        if escape_name and not (
                            name.startswith("'") and name.endswith("'")
                        ):
                            if re.search("""[]['"(),:;_]""", name):
                                name = "'%s'" % name.replace("'", "''")
                            else:
                                name = name.replace(" ", "_")
                    result.append(name)

                if isinstance(self, PhyloNode):
                    if with_distances and top_node.length is not None:
                        # result.append(":%s" % top_node.length)
                        result[-1] = f"{result[-1]}:{top_node.length}"

                result.append(",")

        len_result = len(result)
        if len_result == 2:  # single node no name
            if semicolon:
                return ";"
            else:
                return ""
        elif len_result == 3:  # single node with name
            if semicolon:
                return f"{result[1]};"
            else:
                return result[1]
        else:
            if semicolon:
                result[-1] = ";"
            else:
                result.pop(-1)
            return "".join(result)

    def remove_node(self, target):
        """Removes node by identity instead of value.

        Returns True if node was present, False otherwise.
        """
        to_delete = None
        for i, curr_node in enumerate(self.children):
            if curr_node is target:
                to_delete = i
                break
        if to_delete is None:
            return False
        else:
            del self[to_delete]
            return True

    def get_edge_names(
        self, tip1name, tip2name, clade=True, stem=False, outgroup_name=None
    ):
        """Return the list of stem and/or sub tree (clade) edge name(s).
        This is done by finding the common intersection, and then getting
        the list of names. If the clade traverses the root, then use the
        outgroup_name argument to ensure valid specification.

        Parameters
        ----------
        tip1/2name
            edge 1/2 names
        stem
            whether the name of the clade stem edge is returned.
        clade
            whether the names of the edges within the clade are
            returned
        outgroup_name
            if provided the calculation is done on a version of
            the tree re-rooted relative to the provided tip.

        Usage:
            The returned list can be used to specify subtrees for special
            parameterisation. For instance, say you want to allow the primates
            to have a different value of a particular parameter. In this case,
            provide the results of this method to the parameter controller
            method `set_param_rule()` along with the parameter name etc..
        """
        # If outgroup specified put it at the top of the tree so that clades are
        # defined by their distance from it.  This makes a temporary tree with
        # a named edge at it's root, but it's only used here then discarded.
        if outgroup_name is not None:
            outgroup = self.get_node_matching_name(outgroup_name)
            if outgroup.children:
                raise TreeError(f"Outgroup ({outgroup_name!r}) must be a tip")
            self = outgroup.unrooted_deepcopy()

        join_edge = self.get_connecting_node(tip1name, tip2name)

        edge_names = []

        if stem:
            if join_edge.isroot():
                raise TreeError(
                    f"LCA({tip1name},{tip2name}) is the root and so has no stem"
                )
            else:
                edge_names.append(join_edge.name)

        if clade:
            # get the list of names contained by join_edge
            for child in join_edge.children:
                branchnames = child.get_node_names(includeself=1)
                edge_names.extend(branchnames)

        return edge_names

    def _getNeighboursExcept(self, parent=None):
        # For walking the tree as if it was unrooted.
        return [
            c
            for c in (tuple(self.children) + (self.parent,))
            if c is not None and c is not parent
        ]

    def _get_distances(self, endpoints=None):
        """Iteratively calcluates all of the root-to-tip and tip-to-tip
        distances, resulting in a tuple of:
            - A list of (name, path length) pairs.
            - A dictionary of (tip1,tip2):distance pairs
        """
        # linearize the tips in postorder.
        # .__start, .__stop compose the slice in tip_order.
        if endpoints is None:
            tip_order = list(self.tips())
        else:
            tip_order = []
            for i, name in enumerate(endpoints):
                node = self.get_node_matching_name(name)
                tip_order.append(node)
        for i, node in enumerate(tip_order):
            node.__start, node.__stop = i, i + 1

        num_tips = len(tip_order)
        result = {}
        # distances from tip to curr node
        tipdistances = zeros((num_tips), float)

        def update_result():
            # set tip_tip distance between tips of different child
            for child1, child2 in combinations(node.children, 2):
                for tip1 in range(child1.__start, child1.__stop):
                    for tip2 in range(child2.__start, child2.__stop):
                        name1 = tip_order[tip1].name
                        name2 = tip_order[tip2].name
                        result[(name1, name2)] = tipdistances[tip1] + tipdistances[tip2]
                        result[(name2, name1)] = tipdistances[tip1] + tipdistances[tip2]

        for node in self.traverse(self_before=False, self_after=True):
            if not node.children:
                continue
            # subtree with solved child wedges
            starts, stops = [], []  # to calc ._start and ._stop for curr node
            for child in node.children:
                if hasattr(child, "length") and child.length is not None:
                    child_len = child.length
                else:
                    child_len = 1  # default length
                tipdistances[child.__start : child.__stop] += child_len
                starts.append(child.__start)
                stops.append(child.__stop)
            node.__start, node.__stop = min(starts), max(stops)
            # update result if nessessary
            if len(node.children) > 1:  # not single child
                update_result()

        from_root = []
        for i, n in enumerate(tip_order):
            from_root.append((n.name, tipdistances[i]))
        return from_root, result

    def get_distances(self, endpoints=None):
        """The distance matrix as a dictionary.

        Usage:
            Grabs the branch lengths (evolutionary distances) as
            a complete matrix (i.e. a,b and b,a).
        """

        (root_dists, endpoint_dists) = self._get_distances(endpoints)
        return endpoint_dists

    def set_max_tip_tip_distance(self):
        """Propagate tip distance information up the tree

        This method was originally implemented by Julia Goodrich with the intent
        of being able to determine max tip to tip distances between nodes on
        large trees efficiently. The code has been modified to track the
        specific tips the distance is between
        """
        for n in self.postorder():
            if n.is_tip():
                n.MaxDistTips = [[0.0, n.name], [0.0, n.name]]
            else:
                if len(n.children) == 1:
                    tip_a, tip_b = n.children[0].MaxDistTips
                    tip_a[0] += n.children[0].length or 0.0
                    tip_b[0] += n.children[0].length or 0.0
                else:
                    tip_info = [(max(c.MaxDistTips), c) for c in n.children]
                    dists = [i[0][0] for i in tip_info]
                    best_idx = argsort(dists)[-2:]
                    tip_a, child_a = tip_info[best_idx[0]]
                    tip_b, child_b = tip_info[best_idx[1]]
                    tip_a[0] += child_a.length or 0.0
                    tip_b[0] += child_b.length or 0.0
                n.MaxDistTips = [tip_a, tip_b]

    def get_max_tip_tip_distance(self):
        """Returns the max tip tip distance between any pair of tips

        Returns (dist, tip_names, internal_node)
        """
        if not hasattr(self, "MaxDistTips"):
            self.set_max_tip_tip_distance()

        longest = 0.0
        names = [None, None]
        best_node = None
        for n in self.nontips(include_self=True):
            tip_a, tip_b = n.MaxDistTips
            dist = tip_a[0] + tip_b[0]

            if dist > longest:
                longest = dist
                best_node = n
                names = [tip_a[1], tip_b[1]]
        return longest, names, best_node

    def max_tip_tip_distance(self):
        """returns the max distance between any pair of tips

        Also returns the tip names  that it is between as a tuple"""
        distmtx, tip_order = self.tip_to_tip_distances()
        idx_max = divmod(distmtx.argmax(), distmtx.shape[1])
        max_pair = (tip_order[idx_max[0]].name, tip_order[idx_max[1]].name)
        return distmtx[idx_max], max_pair

    def _get_sub_tree(
        self, included_names, constructor=None, keep_root=False, tipsonly=False
    ):
        """An equivalent node with possibly fewer children, or None"""

        # Renumber autonamed edges
        if constructor is None:
            constructor = self._default_tree_constructor()

        if self.name in included_names and not tipsonly:
            return self.deepcopy(constructor=constructor)
        elif self.name in included_names and self.istip():
            return self.deepcopy(constructor=constructor)
        else:
            # don't need to pass keep_root to children, though
            # internal nodes will be elminated this way
            children = []
            for child in self.children:
                st = child._get_sub_tree(included_names, constructor, tipsonly=tipsonly)
                children.append(st)

            children = [child for child in children if child is not None]
            if len(children) == 0:
                result = None
            elif len(children) == 1 and not keep_root:
                # Merge parameter dictionaries by adding lengths and making
                # weighted averages of other parameters.  This should probably
                # be moved out of here into a ParameterSet class (Model?) or
                # tree subclass.
                params = {}
                child = children[0]
                if self.length is not None and child.length is not None:
                    shared_params = [
                        n
                        for (n, v) in list(self.params.items())
                        if v is not None
                        and child.params.get(n) is not None
                        and n != "length"
                    ]
                    length = self.length + child.length
                    if length:
                        params = {}
                        for n in shared_params:
                            self_val = self.params[n]
                            child_val = child.params[n]
                            is_scalar = True
                            for i in (self_val, child_val):
                                if not isinstance(i, numbers.Number):
                                    is_scalar = False
                                    break
                            if is_scalar:
                                val = (
                                    self_val * self.length + child_val * child.length
                                ) / length
                            else:
                                val = self_val
                            params[n] = val

                        params["length"] = length
                result = child
                result.params = params
            else:
                result = constructor(self, tuple(children))
        return result

    def get_sub_tree(
        self, name_list, ignore_missing=False, keep_root=False, tipsonly=False
    ):
        """A new instance of a sub tree that contains all the otus that are
        listed in name_list.

        Parameters
        ----------
        ignore_missing
            if False, get_sub_tree will raise a ValueError if
            name_list contains names that aren't nodes in the tree
        keep_root
            if False, the root of the subtree will be the last common
            ancestor of all nodes kept in the subtree. Root to tip distance is
            then (possibly) different from the original tree. If True, the root to
            tip distance remains constant, but root may only have one child node.
        tipsonly
            only tip names matching name_list are allowed

        """
        edge_names = set(self.get_node_names(includeself=1, tipsonly=tipsonly))
        if not ignore_missing:
            # this may take a long time
            for name in name_list:
                if name not in edge_names:
                    raise ValueError(f"edge {name!r} not found in tree")

        new_tree = self._get_sub_tree(name_list, keep_root=keep_root, tipsonly=tipsonly)
        if new_tree is None:
            raise TreeError("no tree created in make sub tree")
        elif new_tree.istip():
            raise TreeError("only a tip was returned from selecting sub tree")
        else:
            new_tree.name = "root"
            # keep unrooted
            if len(self.children) > 2:
                new_tree = new_tree.unrooted()
            return new_tree

    def _edgecount(self, parent, cache):
        """ "The number of edges beyond 'parent' in the direction of 'self',
        unrooted"""
        neighbours = self._getNeighboursExcept(parent)
        key = (id(parent), id(self))
        if key not in cache:
            cache[key] = 1 + sum(
                [child._edgecount(self, cache) for child in neighbours]
            )
        return cache[key]

    def _imbalance(self, parent, cache):
        """The edge count from here, (except via 'parent'), divided into that
        from the heaviest neighbour, and that from the rest of them.  'cache'
        should be a dictionary that can be shared by calls to self.edgecount,
        it stores the edgecount for each node (from self) without having to
        put it on the tree itself."""
        max_weight = 0
        total_weight = 0
        for child in self._getNeighboursExcept(parent):
            weight = child._edgecount(self, cache)
            total_weight += weight
            if weight > max_weight:
                max_weight = weight
                biggest_branch = child
        return (max_weight, total_weight - max_weight, biggest_branch)

    def _sorted(self, sort_order):
        """Score all the edges, sort them, and return minimum score and a
        sorted tree.
        """
        # Only need to duplicate whole tree because of .parent pointers

        constructor = self._default_tree_constructor()

        if not self.children:
            tree = self.deepcopy(constructor)
            score = sort_order.index(self.name)
        else:
            scored_subtrees = [child._sorted(sort_order) for child in self.children]
            scored_subtrees.sort()
            children = tuple(
                [child.deepcopy(constructor) for (score, child) in scored_subtrees]
            )
            tree = constructor(self, children)

            non_null_scores = [
                score for (score, child) in scored_subtrees if score is not None
            ]
            score = (non_null_scores or [None])[0]
        return (score, tree)

    def sorted(self, sort_order=None):
        """An equivalent tree sorted into a standard order. If this is not
        specified then alphabetical order is used.  At each node starting from
        root, the algorithm will try to put the descendant which contains the
        lowest scoring tip on the left.
        """
        sort_order = sort_order or []
        tip_names = self.get_tip_names()
        tip_names.sort()
        full_sort_order = sort_order + tip_names
        (score, tree) = self._sorted(full_sort_order)
        return tree

    def _ascii_art(self, char1="-", show_internal=True, compact=False):
        LEN = 10
        PAD = " " * LEN
        PA = " " * (LEN - 1)
        namestr = self.name or ""  # prevents name of NoneType
        if self.children:
            mids = []
            result = []
            for c in self.children:
                if c is self.children[0]:
                    char2 = "/"
                elif c is self.children[-1]:
                    char2 = "\\"
                else:
                    char2 = "-"
                (clines, mid) = c._ascii_art(char2, show_internal, compact)
                mids.append(mid + len(result))
                result.extend(clines)
                if not compact:
                    result.append("")
            if not compact:
                result.pop()
            (lo, hi, end) = (mids[0], mids[-1], len(result))
            prefixes = (
                [PAD] * (lo + 1) + [PA + "|"] * (hi - lo - 1) + [PAD] * (end - hi)
            )
            mid = (lo + hi) // 2
            prefixes[mid] = char1 + "-" * (LEN - 2) + prefixes[mid][-1]
            result = [p + l for (p, l) in zip(prefixes, result)]
            if show_internal:
                stem = result[mid]
                result[mid] = stem[0] + namestr + stem[len(namestr) + 1 :]
            return (result, mid)
        else:
            return ([char1 + "-" + namestr], 0)

    def ascii_art(self, show_internal=True, compact=False):
        """Returns a string containing an ascii drawing of the tree.

        Parameters
        ----------
        show_internal
            includes internal edge names.
        compact
            use exactly one line per tip.

        """
        (lines, mid) = self._ascii_art(show_internal=show_internal, compact=compact)
        return "\n".join(lines)

    def _getXmlLines(self, indent=0, parent_params=None):
        """Return the xml strings for this edge."""
        params = {}
        if parent_params is not None:
            params.update(parent_params)
        pad = "  " * indent
        xml = [f"{pad}<clade>"]
        if self.name_loaded:
            xml.append(f"{pad}   <name>{self.name}</name>")
        for (n, v) in list(self.params.items()):
            if v == params.get(n, None):
                continue
            xml.append(f"{pad}   <param><name>{n}</name><value>{v}</value></param>")
            params[n] = v
        for child in self.children:
            xml.extend(child._getXmlLines(indent + 1, params))
        xml.append(pad + "</clade>")
        return xml

    def get_xml(self):
        """Return XML formatted tree string."""
        header = ['<?xml version="1.0"?>']  # <!DOCTYPE ...
        return "\n".join(header + self._getXmlLines())

    def write(self, filename, with_distances=True, format=None):
        """Save the tree to filename

        Parameters
        ----------
        filename
            self
        with_distances
            whether branch lengths are included in string.
        format
            default is newick, xml and json are alternate. Argument overrides
            the filename suffix. All attributes are saved in the xml format.
            Value overrides the file name suffix.

        Notes
        -----
        Only the cogent3 json and xml tree formats are supported.

        """
        file_format, _ = get_format_suffixes(filename)
        format = format or file_format
        if format == "json":
            with atomic_write(filename, mode="wt") as f:
                f.write(self.to_json())
            return

        xml = format.lower() == "xml" if format else filename.lower().endswith("xml")
        data = self.get_xml() if xml else self.get_newick(with_distances=with_distances)

        with atomic_write(filename, mode="wt") as outf:
            outf.writelines(data)

    def get_node_names(self, includeself=True, tipsonly=False):
        """Return a list of edges from this edge - may or may not include self.
        This node (or first connection) will be the first, and then they will
        be listed in the natural traverse order.

        Parameters
        ----------
        includeself : bool
            excludes self.name from the result

        tipsonly : bool
            only tips returned
        """
        if tipsonly:
            nodes = self.traverse(self_before=False, self_after=False)
        else:
            nodes = list(self.traverse())
            if not includeself:
                nodes.remove(self)
        return [node.name for node in nodes]

    def get_tip_names(self, includeself=False):
        """return the list of the names of all tips contained by this edge"""
        return self.get_node_names(includeself, tipsonly=True)

    def get_edge_vector(self, include_root=True):
        """Collect the list of edges in postfix order

        Parameters
        ----------
        include_root
            specifies whether root edge included

        """
        if include_root:
            result = [n for n in self.traverse(False, True)]
        else:
            result = [n for n in self.traverse(False, True) if not n.isroot()]
        return result

    def _get_node_matching_name(self, name):
        """
        find the edge with the name, or return None
        """
        for node in self.traverse(self_before=True, self_after=False):
            if node.name == name:
                return node
        return None

    def get_node_matching_name(self, name):
        node = self._get_node_matching_name(name)
        if node is None:
            raise TreeError(f"No node named '{name}' in {self.get_tip_names()}")
        return node

    def get_connecting_node(self, name1, name2):
        """Finds the last common ancestor of the two named edges."""
        edge1 = self.get_node_matching_name(name1)
        edge2 = self.get_node_matching_name(name2)
        lca = edge1.last_common_ancestor(edge2)
        if lca is None:
            raise TreeError(f"No LCA found for {name1} and {name2}")
        return lca

    def get_connecting_edges(self, name1, name2):
        """returns a list of edges connecting two nodes.

        If both are tips, the LCA is excluded from the result."""
        edge1 = self.get_node_matching_name(name1)
        edge2 = self.get_node_matching_name(name2)
        include_parent = False if edge1.istip() and edge2.istip() else True

        LCA = self.get_connecting_node(name1, name2)
        node_path = [edge1]
        node_path.extend(edge1.ancestors())
        # remove nodes deeper than the LCA
        LCA_ind = node_path.index(LCA)
        node_path = node_path[: LCA_ind + 1]
        # remove LCA and deeper nodes from anc list of other
        anc2 = edge2.ancestors()
        LCA_ind = anc2.index(LCA)
        anc2 = anc2[:LCA_ind]
        anc2.reverse()
        node_path.extend(anc2)
        node_path.append(edge2)
        if not include_parent:
            node_path.remove(LCA)
        return node_path

    def get_param_value(self, param, edge):
        """returns the parameter value for named edge"""
        return self.get_node_matching_name(edge).params[param]

    def set_param_value(self, param, edge, value):
        """set's the value for param at named edge"""
        self.get_node_matching_name(edge).params[param] = value

    def reassign_names(self, mapping, nodes=None):
        """Reassigns node names based on a mapping dict

        mapping : dict, old_name -> new_name
        nodes : specific nodes for renaming (such as just tips, etc...)
        """
        if nodes is None:
            nodes = self.traverse()

        for n in nodes:
            if n.name in mapping:
                n.name = mapping[n.name]

    def multifurcating(self, num, eps=None, constructor=None, name_unnamed=False):
        """return a new tree with every node having num or few children

        Parameters
        ----------
        num : int
            the number of children a node can have max
        eps : float
            default branch length to set if self or constructor is of
            PhyloNode type
        constructor
            a TreeNode or subclass constructor. If None, uses self
        name_unnamed : bool
            names unnamed nodes
        """
        if num < 2:
            raise TreeError("Minimum number of children must be >= 2")

        if eps is None:
            eps = 0.0

        if constructor is None:
            constructor = self.__class__

        if hasattr(constructor, "length"):
            set_branchlength = True
        else:
            set_branchlength = False

        new_tree = self.copy()

        for n in new_tree.preorder(include_self=True):
            while len(n.children) > num:
                new_node = constructor(children=n.children[-num:])

                if set_branchlength:
                    new_node.length = eps

                n.append(new_node)

        if name_unnamed:
            alpha = "abcdefghijklmnopqrstuvwxyz"
            alpha += alpha.upper()
            base = "AUTOGENERATED_NAME_%s"

            # scale the random names by tree size
            s = int(ceil(log(len(new_tree.tips()))))

            for n in new_tree.nontips():
                if n.name is None:
                    n.name = base % "".join([choice(alpha) for i in range(s)])

        return new_tree

    def bifurcating(self, eps=None, constructor=None, name_unnamed=False):
        """Wrap multifurcating with a num of 2"""
        return self.multifurcating(2, eps, constructor, name_unnamed)

    def get_nodes_dict(self):
        """Returns a dict keyed by node name, value is node

        Will raise TreeError if non-unique names are encountered
        """
        res = {}

        for n in self.traverse():
            if n.name in res:
                raise TreeError("get_nodes_dict requires unique node names")
            else:
                res[n.name] = n

        return res

    def subset(self):
        """Returns set of names that descend from specified node"""
        return frozenset([i.name for i in self.tips()])

    def subsets(self):
        """Returns all sets of names that come from specified node and its kids"""
        sets = []
        for i in self.traverse(self_before=False, self_after=True, include_self=False):
            if not i.children:
                i.__leaf_set = frozenset([i.name])
            else:
                leaf_set = reduce(or_, [c.__leaf_set for c in i.children])
                if len(leaf_set) > 1:
                    sets.append(leaf_set)
                i.__leaf_set = leaf_set
        return frozenset(sets)

    def compare_by_subsets(self, other, exclude_absent_taxa=False):
        """Returns fraction of overlapping subsets where self and other differ.

        Other is expected to be a tree object compatible with PhyloNode.

        Note: names present in only one of the two trees will count as
        mismatches: if you don't want this behavior, strip out the non-matching
        tips first.
        """
        self_sets, other_sets = self.subsets(), other.subsets()
        if exclude_absent_taxa:
            in_both = self.subset() & other.subset()
            self_sets = [i & in_both for i in self_sets]
            self_sets = frozenset([i for i in self_sets if len(i) > 1])
            other_sets = [i & in_both for i in other_sets]
            other_sets = frozenset([i for i in other_sets if len(i) > 1])
        total_subsets = len(self_sets) + len(other_sets)
        intersection_length = len(self_sets & other_sets)
        if not total_subsets:  # no common subsets after filtering, so max dist
            return 1
        return 1 - 2 * intersection_length / float(total_subsets)

    def tip_to_tip_distances(self, default_length=1):
        """Returns distance matrix between all pairs of tips, and a tip order.

        Warning: .__start and .__stop added to self and its descendants.

        tip_order contains the actual node objects, not their names (may be
        confusing in some cases).
        """
        # linearize the tips in postorder.
        # .__start, .__stop compose the slice in tip_order.
        tip_order = list(self.tips())
        for i, tip in enumerate(tip_order):
            tip.__start, tip.__stop = i, i + 1

        num_tips = len(tip_order)
        result = zeros((num_tips, num_tips), float)  # tip by tip matrix
        # distances from tip to curr node
        tipdistances = zeros((num_tips), float)

        def update_result():
            # set tip_tip distance between tips of different child
            for child1, child2 in combinations(node.children, 2):
                for tip1 in range(child1.__start, child1.__stop):
                    for tip2 in range(child2.__start, child2.__stop):
                        result[tip1, tip2] = tipdistances[tip1] + tipdistances[tip2]

        for node in self.traverse(self_before=False, self_after=True):
            if not node.children:
                continue
            # subtree with solved child wedges
            starts, stops = [], []  # to calc ._start and ._stop for curr node
            for child in node.children:
                if hasattr(child, "length") and child.length is not None:
                    child_len = child.length
                else:
                    child_len = default_length
                tipdistances[child.__start : child.__stop] += child_len
                starts.append(child.__start)
                stops.append(child.__stop)
            node.__start, node.__stop = min(starts), max(stops)
            # update result if nessessary
            if len(node.children) > 1:  # not single child
                update_result()
        return result + result.T, tip_order

    def compare_by_tip_distances(self, other, dist_f=distance_from_r):
        """Compares self to other using tip-to-tip distance matrices.

        Value returned is dist_f(m1, m2) for the two matrices. Default is
        to use the Pearson correlation coefficient, with +1 giving a distance
        of 0 and -1 giving a distance of +1 (the madimum possible value).
        Depending on the application, you might instead want to use
        distance_from_r_squared, which counts correlations of both +1 and -1
        as identical (0 distance).

        Note: automatically strips out the names that don't match (this is
        necessary for this method because the distance between non-matching
        names and matching names is undefined in the tree where they don't
        match, and because we need to reorder the names in the two trees to
        match up the distance matrices).
        """
        self_names = [i.name for i in self.tips()]
        other_names = [i.name for i in other.tips()]
        common_names = frozenset(self_names) & frozenset(other_names)
        if not common_names:
            raise ValueError("No names in common between the two trees." "")
        if len(common_names) <= 2:
            return 1  # the two trees must match by definition in this case
        # figure out correct order of the two name matrices
        self_order = [self_names.index(i) for i in common_names]
        other_order = [other_names.index(i) for i in common_names]
        self_matrix = self.tip_to_tip_distances()[0][self_order][:, self_order]
        other_matrix = other.tip_to_tip_distances()[0][other_order][:, other_order]
        return dist_f(self_matrix, other_matrix)

    def get_figure(self, style="square", **kwargs):
        """
        gets Dendrogram for plotting the phylogeny

        Parameters
        ----------
        style : string
            'square', 'angular', 'radial' or 'circular'
        kwargs
            arguments passed to Dendrogram constructor
        """
        from cogent3.draw.dendrogram import Dendrogram

        style = style.lower()
        types = ("square", "circular", "angular", "radial")
        if style not in types:
            raise ValueError(f"{style} not in supported types {types}")

        return Dendrogram(self, style=style, **kwargs)


class PhyloNode(TreeNode):
    def __init__(self, *args, **kwargs):
        length = kwargs.get("length", None)
        params = kwargs.get("params", {})
        if "length" not in params:
            params["length"] = length
        kwargs["params"] = params
        super(PhyloNode, self).__init__(*args, **kwargs)

    def _set_length(self, value):
        if not hasattr(self, "params"):
            self.params = {}
        self.params["length"] = value

    def _get_length(self):
        return self.params.get("length", None)

    length = property(_get_length, _set_length)

    def __str__(self):
        """Returns string version of self, with names and distances."""
        return self.get_newick(with_distances=True)

    def distance(self, other):
        """Returns branch length between self and other."""
        # never any length between self and other
        if self is other:
            return 0
        # otherwise, find self's ancestors and find the first ancestor of
        # other that is in the list
        self_anc = self.ancestors()
        self_anc_dict = dict([(id(n), n) for n in self_anc])
        self_anc_dict[id(self)] = self

        count = 0
        while other is not None:
            if id(other) in self_anc_dict:
                # found the first shared ancestor -- need to sum other branch
                curr = self
                while curr is not other:
                    if curr.length:
                        count += curr.length
                    curr = curr._parent
                return count
            else:
                if other.length:
                    count += other.length
            other = other._parent
        return None

    def total_descending_branch_length(self):
        """Returns total descending branch length from self"""
        return sum(
            [
                n.length
                for n in self.traverse(include_self=False)
                if n.length is not None
            ]
        )

    def total_length(self):
        """returns the sum of all branch lengths in tree"""
        root = self.root()
        if root is None:
            raise ValueError("no root to this tree!")

        return root.total_descending_branch_length()

    def tips_within_distance(self, distance):
        """Returns tips within specified distance from self

        Branch lengths of None will be interpreted as 0
        """

        def get_distance(d1, d2):
            if d2 is None:
                return d1
            else:
                return d1 + d2

        to_process = [(self, 0.0)]
        tips_to_save = []

        seen = set([id(self)])
        while to_process:
            curr_node, curr_dist = to_process.pop(0)

            # have we've found a tip within distance?
            if curr_node.is_tip() and curr_node != self:
                tips_to_save.append(curr_node)
                continue

            # add the parent node if it is within distance
            parent_dist = get_distance(curr_dist, curr_node.length)
            if (
                curr_node.parent is not None
                and parent_dist <= distance
                and id(curr_node.parent) not in seen
            ):
                to_process.append((curr_node.parent, parent_dist))
                seen.add(id(curr_node.parent))

            # add children if we haven't seen them and if they are in distance
            for child in curr_node.children:
                if id(child) in seen:
                    continue
                seen.add(id(child))

                child_dist = get_distance(curr_dist, child.length)
                if child_dist <= distance:
                    to_process.append((child, child_dist))

        return tips_to_save

    def prune(self):
        """Reconstructs correct tree after nodes have been removed.

        Internal nodes with only one child will be removed and new connections
        and Branch lengths will be made to reflect change.
        """
        # traverse tree to decide nodes to be removed.
        nodes_to_remove = []
        for node in self.traverse():
            if (node.parent is not None) and (len(node.children) == 1):
                nodes_to_remove.append(node)
        for node in nodes_to_remove:
            # save current parent
            curr_parent = node.parent
            # save child
            child = node.children[0]
            # remove current node by setting parent to None
            node.parent = None
            # Connect child to current node's parent
            child.parent = curr_parent
            # Add the length of the removed node to the length of the Child
            if child.length is None or node.length is None:
                child.length = child.length or node.length
            else:
                child.length = child.length + node.length

    def unrooted_deepcopy(self, constructor=None, parent=None):
        # walks the tree unrooted-style, ie: treating self.parent as just
        # another child 'parent' is where we got here from, ie: the neighbour
        # that we don't need to explore.
        if constructor is None:
            constructor = self._default_tree_constructor()

        neighbours = self._getNeighboursExcept(parent)
        children = []
        for child in neighbours:
            children.append(child.unrooted_deepcopy(constructor, parent=self))

        # we might be walking UP the tree, so:
        if parent is None:
            # base edge
            edge = None
        elif parent.parent is self:
            # self's parent is becoming self's child, and edge params are stored
            # by the child
            edge = parent
        else:
            assert parent is self.parent
            edge = self

        result = constructor(edge, tuple(children))
        if parent is None:
            result.name = "root"
        return result

    def balanced(self):
        """Tree 'rooted' here with no neighbour having > 50% of the edges.

        Usage:
            Using a balanced tree can substantially improve performance of
            the likelihood calculations. Note that the resulting tree has a
            different orientation with the effect that specifying clades or
            stems for model parameterisation should be done using the
            'outgroup_name' argument.
        """
        # this should work OK on ordinary 3-way trees, not so sure about
        # other cases.  Given 3 neighbours, if one has > 50% of edges it
        # can only improve things to divide it up, worst case:
        # (51),25,24 -> (50,1),49.
        # If no neighbour has >50% we can't improve on where we are, eg:
        # (49),25,26 -> (20,19),51
        last_edge = None
        edge = self
        known_weight = 0
        cache = {}
        while 1:
            (max_weight, remaining_weight, next_edge) = edge._imbalance(
                last_edge, cache
            )
            known_weight += remaining_weight
            if max_weight <= known_weight + 2:
                break
            last_edge = edge
            edge = next_edge
            known_weight += 1
        return edge.unrooted_deepcopy()

    def same_topology(self, other):
        """Tests whether two trees have the same topology."""
        tip_names = self.get_tip_names()
        root_at = tip_names[0]
        me = self.rooted_with_tip(root_at).sorted(tip_names)
        them = other.rooted_with_tip(root_at).sorted(tip_names)
        return self is other or me.same_shape(them)

    def unrooted(self):
        """A tree with at least 3 children at the root."""
        constructor = self._default_tree_constructor()
        need_to_expand = len(self.children) < 3
        new_children = []
        for oldnode in self.children:
            if oldnode.children and need_to_expand:
                for sib in oldnode.children:
                    sib = sib.deepcopy(constructor)
                    if sib.length is not None and oldnode.length is not None:
                        sib.length += oldnode.length
                    new_children.append(sib)
                need_to_expand = False
            else:
                new_children.append(oldnode.deepcopy(constructor))
        return constructor(self, new_children)

    def rooted_at(self, edge_name):
        """Return a new tree rooted at the provided node.

        Usage:
            This can be useful for drawing unrooted trees with an orientation
            that reflects knowledge of the true root location.
        """
        newroot = self.get_node_matching_name(edge_name)
        if not newroot.children:
            raise TreeError(f"Can't use a tip ({repr(edge_name)}) as the root")
        return newroot.unrooted_deepcopy()

    def rooted_with_tip(self, outgroup_name):
        """A new tree with the named tip as one of the root's children"""
        tip = self.get_node_matching_name(outgroup_name)
        return tip.parent.unrooted_deepcopy()

    def root_at_midpoint(self):
        """return a new tree rooted at midpoint of the two tips farthest apart

        this fn doesn't preserve the internal node naming or structure,
        but does keep tip to tip distances correct.  uses unrooted_deepcopy()
        """
        # max_dist, tip_names = tree.max_tip_tip_distance()
        # this is slow

        max_dist, tip_names = self.max_tip_tip_distance()
        half_max_dist = max_dist / 2.0
        if max_dist == 0.0:  # only pathological cases with no lengths
            return self.unrooted_deepcopy()
        # print tip_names
        tip1 = self.get_node_matching_name(tip_names[0])
        tip2 = self.get_node_matching_name(tip_names[1])
        lca = self.get_connecting_node(tip_names[0], tip_names[1])  # last comm ancestor
        if tip1.distance(lca) > half_max_dist:
            climb_node = tip1
        else:
            climb_node = tip2

        dist_climbed = 0.0
        while dist_climbed + climb_node.length < half_max_dist:
            dist_climbed += climb_node.length
            climb_node = climb_node.parent

        # now midpt is either at on the branch to climb_node's  parent
        # or midpt is at climb_node's parent
        # print dist_climbed, half_max_dist, 'dists cl hamax'
        if dist_climbed + climb_node.length == half_max_dist:
            # climb to midpoint spot
            climb_node = climb_node.parent
            if climb_node.is_tip():
                raise RuntimeError("error trying to root tree at tip")
            else:
                # print climb_node.name, 'clmb node'
                return climb_node.unrooted_deepcopy()

        else:
            # make a new node on climb_node's branch to its parent
            old_br_len = climb_node.length
            new_root = type(self)()
            new_root.parent = climb_node.parent
            climb_node.parent = new_root
            climb_node.length = half_max_dist - dist_climbed
            new_root.length = old_br_len - climb_node.length
            return new_root.unrooted_deepcopy()

    def _find_midpoint_nodes(self, max_dist, tip_pair):
        """returns the nodes surrounding the max_tip_tip_distance midpoint

        WAS used for midpoint rooting.  ORPHANED NOW
        max_dist: The maximum distance between any 2 tips
        tip_pair: Names of the two tips associated with max_dist
        """
        half_max_dist = max_dist / 2.0
        # get a list of the nodes that separate the tip pair
        node_path = self.get_connecting_edges(tip_pair[0], tip_pair[1])
        tip1 = self.get_node_matching_name(tip_pair[0])
        for index, node in enumerate(node_path):
            dist = tip1.distance(node)
            if dist > half_max_dist:
                return node, node_path[index - 1]

    def set_tip_distances(self):
        """Sets distance from each node to the most distant tip."""
        for node in self.traverse(self_before=False, self_after=True):
            if node.children:
                node.TipDistance = max(
                    [c.length + c.TipDistance for c in node.children]
                )
            else:
                node.TipDistance = 0

    def scale_branch_lengths(self, max_length=100, ultrametric=False):
        """Scales BranchLengths in place to integers for ascii output.

        Warning: tree might not be exactly the length you specify.

        Set ultrametric=True if you want all the root-tip distances to end
        up precisely the same.
        """
        self.set_tip_distances()
        orig_max = max([n.TipDistance for n in self.traverse()])
        if not ultrametric:  # easy case -- just scale and round
            for node in self.traverse():
                curr = node.length
                if curr is not None:
                    node.ScaledBranchLength = max(
                        1, int(round(1.0 * curr / orig_max * max_length))
                    )
        else:  # hard case -- need to make sure they all line up at the end
            for node in self.traverse(self_before=False, self_after=True):
                if not node.children:  # easy case: ignore tips
                    node.DistanceUsed = 0
                    continue
                # if we get here, we know the node has children
                # figure out what distance we want to set for this node
                ideal_distance = int(round(node.TipDistance / orig_max * max_length))
                min_distance = max([c.DistanceUsed for c in node.children]) + 1
                distance = max(min_distance, ideal_distance)
                for c in node.children:
                    c.ScaledBranchLength = distance - c.DistanceUsed
                node.DistanceUsed = distance
        # reset the BranchLengths
        for node in self.traverse(self_before=True, self_after=False):
            if node.length is not None:
                node.length = node.ScaledBranchLength
            if hasattr(node, "ScaledBranchLength"):
                del node.ScaledBranchLength
            if hasattr(node, "DistanceUsed"):
                del node.DistanceUsed
            if hasattr(node, "TipDistance"):
                del node.TipDistance

    def _get_distances(self, endpoints=None):
        """Iteratively calcluates all of the root-to-tip and tip-to-tip
        distances, resulting in a tuple of:
            - A list of (name, path length) pairs.
            - A dictionary of (tip1,tip2):distance pairs
        """
        # linearize the tips in postorder.
        # .__start, .__stop compose the slice in tip_order.
        if endpoints is None:
            tip_order = list(self.tips())
        else:
            tip_order = []
            for i, name in enumerate(endpoints):
                node = self.get_node_matching_name(name)
                tip_order.append(node)
        for i, node in enumerate(tip_order):
            node.__start, node.__stop = i, i + 1

        num_tips = len(tip_order)
        result = {}
        # distances from tip to curr node
        tipdistances = zeros((num_tips), float)

        def update_result():
            # set tip_tip distance between tips of different child
            for child1, child2 in combinations(node.children, 2):
                for tip1 in range(child1.__start, child1.__stop):
                    for tip2 in range(child2.__start, child2.__stop):
                        name1 = tip_order[tip1].name
                        name2 = tip_order[tip2].name
                        result[(name1, name2)] = tipdistances[tip1] + tipdistances[tip2]
                        result[(name2, name1)] = tipdistances[tip1] + tipdistances[tip2]

        for node in self.traverse(self_before=False, self_after=True):
            if not node.children:
                continue
            # subtree with solved child wedges
            starts, stops = [], []  # to calc ._start and ._stop for curr node
            for child in node.children:
                if hasattr(child, "length") and child.length is not None:
                    child_len = child.length
                else:
                    child_len = 1  # default length
                tipdistances[child.__start : child.__stop] += child_len
                starts.append(child.__start)
                stops.append(child.__stop)
            node.__start, node.__stop = min(starts), max(stops)
            # update result if nessessary
            if len(node.children) > 1:  # not single child
                update_result()

        from_root = []
        for i, n in enumerate(tip_order):
            from_root.append((n.name, tipdistances[i]))
        return from_root, result

    def get_distances(self, endpoints=None):
        """The distance matrix as a dictionary.

        Usage:
            Grabs the branch lengths (evolutionary distances) as
            a complete matrix (i.e. a,b and b,a)."""

        (root_dists, endpoint_dists) = self._get_distances(endpoints)
        return endpoint_dists

    def tip_to_tip_distances(self, endpoints=None, default_length=1):
        """Returns distance matrix between all pairs of tips, and a tip order.

        Warning: .__start and .__stop added to self and its descendants.

        tip_order contains the actual node objects, not their names (may be
        confusing in some cases).
        """
        all_tips = self.tips()
        if endpoints is None:
            tip_order = list(all_tips)
        else:
            if isinstance(endpoints[0], PhyloNode):
                tip_order = endpoints
            else:
                tip_order = [self.get_node_matching_name(n) for n in endpoints]

        # linearize all tips in postorder
        # .__start, .__stop compose the slice in tip_order.
        for i, node in enumerate(all_tips):
            node.__start, node.__stop = i, i + 1

        # the result map provides index in the result matrix
        result_map = dict([(n.__start, i) for i, n in enumerate(tip_order)])
        num_all_tips = len(all_tips)  # total number of tips
        num_tips = len(tip_order)  # total number of tips in result
        result = zeros((num_tips, num_tips), float)  # tip by tip matrix
        # dist from tip to curr node
        tipdistances = zeros((num_all_tips), float)

        def update_result():
            # set tip_tip distance between tips of different child
            for child1, child2 in combinations(node.children, 2):
                for tip1 in range(child1.__start, child1.__stop):
                    if tip1 not in result_map:
                        continue
                    res_tip1 = result_map[tip1]
                    for tip2 in range(child2.__start, child2.__stop):
                        if tip2 not in result_map:
                            continue
                        result[res_tip1, result_map[tip2]] = (
                            tipdistances[tip1] + tipdistances[tip2]
                        )

        for node in self.traverse(self_before=False, self_after=True):
            if not node.children:
                continue
            # subtree with solved child wedges
            starts, stops = [], []  # to calc ._start and ._stop for curr node
            for child in node.children:
                if hasattr(child, "length") and child.length is not None:
                    child_len = child.length
                else:
                    child_len = default_length
                tipdistances[child.__start : child.__stop] += child_len
                starts.append(child.__start)
                stops.append(child.__stop)
            node.__start, node.__stop = min(starts), max(stops)
            # update result if nessessary
            if len(node.children) > 1:  # not single child
                update_result()
        return result + result.T, tip_order

    def compare_by_tip_distances(
        self, other, sample=None, dist_f=distance_from_r, shuffle_f=shuffle
    ):
        """Compares self to other using tip-to-tip distance matrices.

        Value returned is dist_f(m1, m2) for the two matrices. Default is
        to use the Pearson correlation coefficient, with +1 giving a distance
        of 0 and -1 giving a distance of +1 (the madimum possible value).
        Depending on the application, you might instead want to use
        distance_from_r_squared, which counts correlations of both +1 and -1
        as identical (0 distance).

        Note: automatically strips out the names that don't match (this is
        necessary for this method because the distance between non-matching
        names and matching names is undefined in the tree where they don't
        match, and because we need to reorder the names in the two trees to
        match up the distance matrices).
        """
        self_names = dict([(i.name, i) for i in self.tips()])
        other_names = dict([(i.name, i) for i in other.tips()])
        common_names = frozenset(list(self_names.keys())) & frozenset(
            list(other_names.keys())
        )
        common_names = list(common_names)

        if not common_names:
            raise ValueError("No names in common between the two trees." "")
        if len(common_names) <= 2:
            return 1  # the two trees must match by definition in this case

        if sample is not None:
            shuffle_f(common_names)
            common_names = common_names[:sample]

        self_nodes = [self_names[k] for k in common_names]
        other_nodes = [other_names[k] for k in common_names]

        self_matrix = self.tip_to_tip_distances(endpoints=self_nodes)[0]
        other_matrix = other.tip_to_tip_distances(endpoints=other_nodes)[0]

        return dist_f(self_matrix, other_matrix)


class TreeBuilder(object):
    # Some tree code which isn't needed once the tree is finished.
    # Mostly exists to give edges unique names
    # children must be created before their parents.

    def __init__(self, mutable=False, constructor=PhyloNode):
        self._used_names = {"edge": -1}
        self._known_edges = {}
        self.TreeNodeClass = constructor

    def _unique_name(self, name):
        # Unnamed edges become edge.0, edge.1 edge.2 ...
        # Other duplicates go mouse mouse.2 mouse.3 ...
        if not name:
            name = "edge"
        if name in self._used_names:
            self._used_names[name] += 1
            name += "." + str(self._used_names[name])
            # in case of names like 'edge.1.1'
            name = self._unique_name(name)
        else:
            self._used_names[name] = 1
        return name

    def _params_for_edge(self, edge):
        # default is just to keep it
        return edge.params

    def edge_from_edge(self, edge, children, params=None):
        """Callback for tree-to-tree transforms like get_sub_tree"""
        if edge is None:
            assert not params
            return self.create_edge(children, "root", {}, False)
        else:
            if params is None:
                params = self._params_for_edge(edge)
            return self.create_edge(
                children, edge.name, params, name_loaded=edge.name_loaded
            )

    def create_edge(self, children, name, params, name_loaded=True):
        """Callback for newick parser"""
        if children is None:
            children = []
        node = self.TreeNodeClass(
            children=list(children),
            name=self._unique_name(name),
            name_loaded=name_loaded and (name is not None),
            params=params,
        )
        self._known_edges[id(node)] = node
        return node
