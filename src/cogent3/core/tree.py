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

from __future__ import annotations

import contextlib
import json
import random
import re
import warnings
from copy import deepcopy
from functools import reduce
from itertools import combinations
from operator import or_
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    SupportsIndex,
    TypeVar,
    cast,
    overload,
)

import numpy
import numpy.typing as npt

from cogent3._version import __version__
from cogent3.parse.cogent3_json import load_from_json
from cogent3.parse.newick import parse_string as newick_parse_string
from cogent3.phylo.tree_distance import get_tree_distance_measure
from cogent3.util.deserialise import register_deserialiser
from cogent3.util.io import atomic_write, get_format_suffixes, open_
from cogent3.util.misc import get_object_provenance

if TYPE_CHECKING:  # pragma: no cover
    import os
    import pathlib
    from collections.abc import Callable, Generator, Iterable, Iterator, Sequence

    from typing_extensions import Self

    from cogent3.draw.dendrogram import Dendrogram
    from cogent3.evolve.fast_distance import DistanceMatrix

    PySeq = Sequence
    PySeqStr = PySeq[str]


class TreeError(Exception):
    pass


def _format_node_name(
    node: PhyloNode,
    with_node_names: bool,
    escape_name: bool,
    with_distances: bool,
    with_root_name: bool = False,
) -> str:
    """Helper function to format node name according to parameters"""
    if (node.is_root() and not with_root_name) or (
        not node.is_tip() and not with_node_names
    ):
        node_name = ""
    else:
        node_name = node.name or ""

    if (
        node_name
        and escape_name
        and not (node_name.startswith("'") and node_name.endswith("'"))
    ):
        if re.search("""[]['"(),:;_]""", node_name):
            node_name = "'{}'".format(node_name.replace("'", "''"))
        else:
            node_name = node_name.replace(" ", "_")

    if with_distances and (length := node.length) is not None:
        node_name = f"{node_name}:{length}"

    return node_name


class PhyloNode:
    """Store information about a tree node. Mutable.

    Parameters:
        name: label for the node, assumed to be unique.
        children: list of the node's children.
        parent: parent to this node
        params: dict containing arbitrary parameters for the node.
        name_loaded: ?
    """

    __slots__ = (
        "_parent",
        "children",
        "length",
        "name",
        "name_loaded",
        "params",
        "support",
    )

    _exclude_from_copy = frozenset(["_parent", "children"])

    def __init__(
        self,
        name: str,
        children: Iterable[Self | str] | None = None,
        parent: Self | None = None,
        params: dict[str, Any] | None = None,
        name_loaded: bool = True,
        length: float | None = None,
        support: float | None = None,
    ) -> None:
        """Returns new PhyloNode object."""
        self.name = name
        self.name_loaded = name_loaded
        self.params = params or {}
        self.children: list[Self] = []
        if children:
            self.extend(children)

        self._parent = parent
        if parent is not None and self not in parent.children:
            parent.append(self)

        if "length" in self.params:
            if length is not None and self.params["length"] != length:
                msg = "Got two different lengths in params and as an argument."
                raise ValueError(msg)

            length = self.params.pop("length")

        self.length = length
        self.support = support

    # built-in methods and list interface support
    def __repr__(self) -> str:
        """Returns reconstructable string representation of tree.

        WARNING: Does not currently set the class to the right type.
        """
        return f'Tree("{self.get_newick()}")'

    def __str__(self) -> str:
        """Returns Newick-format string representation of tree."""
        return self.get_newick(with_distances=True)

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

    def __lt__(self, other: PhyloNode) -> bool:
        return self.name < other.name

    def __gt__(self, other: PhyloNode) -> bool:
        return self.name > other.name

    @property
    def source(self) -> str | None:
        return self.params.get("source")

    @source.setter
    def source(self, value: str | None) -> None:
        """Sets the source of the node."""
        if value:
            self.params["source"] = value
        else:
            self.params.pop("source", None)

    def compare_name(self, other: PhyloNode) -> bool:
        """Compares PhyloNode by name"""
        return self is other or self.name == other.name

    def compare_by_names(self, other: PhyloNode) -> bool:
        """Equality test for trees by name"""
        # if they are the same object then they must be the same tree...
        if self is other:
            return True
        self_names = self.get_node_names()
        other_names = other.get_node_names()
        self_names = [v for v in self_names if v]
        self_names.sort()
        other_names = [v for v in other_names if v]
        other_names.sort()
        return self_names == other_names

    def _to_self_child(self, i: Self | str) -> Self:
        """Converts i to self's type, with self as its parent.

        Cleans up refs from i's original parent, but doesn't give self ref to i.
        """
        if isinstance(i, str):
            node = self.__class__(i)
        else:
            if i._parent is not None and i._parent is not self:
                i._parent.children.remove(i)
            node = i
        node._parent = self
        return node

    def append(self, i: Self | str) -> None:
        """Appends i to self.children, in-place, cleaning up refs."""
        self.children.append(self._to_self_child(i))

    def extend(self, items: Iterable[Self | str]) -> None:
        """Extends self.children by items, in-place, cleaning up refs."""
        self.children.extend(self._to_self_child(item) for item in items)

    def insert(self, index: SupportsIndex, i: Self | str) -> None:
        """Inserts an item at specified position in self.children."""
        self.children.insert(index, self._to_self_child(i))

    def pop(self, index: SupportsIndex = -1) -> Self:
        """Returns and deletes child of self at index (default: -1)"""
        result = self.children.pop(index)
        result._parent = None
        return result

    def remove(self, target: Self | str) -> bool:
        """Removes node by name instead of identity.

        Returns True if node was present, False otherwise.
        """
        if isinstance(target, PhyloNode):
            target = target.name
        for _i, curr_node in enumerate(self.children):
            if curr_node.name == target:
                self.remove_node(curr_node)
                return True
        return False

    @overload
    def __getitem__(self, i: SupportsIndex) -> Self: ...
    @overload
    def __getitem__(self, i: slice) -> list[Self]: ...

    def __getitem__(self, i: SupportsIndex | slice) -> Self | list[Self]:
        return self.children[i]

    @overload
    def __setitem__(self, i: SupportsIndex, value: Self | str) -> None: ...
    @overload
    def __setitem__(self, i: slice, value: Iterable[Self | str]) -> None: ...

    def __setitem__(
        self, i: SupportsIndex | slice, value: Self | str | Iterable[Self | str]
    ) -> None:
        """Node[i] = x sets the corresponding item in children."""
        if isinstance(i, slice):
            nodes = cast("Iterable[Self]", value)

            current_children = self.children[i]

            for child in current_children:
                child._parent = None

            self.children[i] = [self._to_self_child(node) for node in nodes]
        else:
            node = cast("Self", value)

            current_child = self.children[i]

            current_child._parent = None
            self.children[i] = self._to_self_child(node)

    def __delitem__(self, i: SupportsIndex | slice) -> None:
        """del node[i] deletes index or slice from self.children."""
        if isinstance(i, slice):
            for child in self.children[i]:
                child._parent = None
        else:
            self.children[i]._parent = None
        del self.children[i]

    def __iter__(self) -> Iterator[Self]:
        """Node iter iterates over the children."""
        return iter(self.children)

    def __len__(self) -> int:
        """Node len returns number of children."""
        return len(self.children)

    @classmethod
    def _copy_node(cls, node: Self, memo: dict[int, Any] | None = None) -> Self:
        result = cls(node.name)
        efc = node._exclude_from_copy
        for k in node.__slots__:
            if k not in efc:
                setattr(result, k, deepcopy(getattr(node, k), memo=memo))
        return result

    def copy(self, memo: dict[int, Any] | None = None) -> Self:
        """Returns a copy of self using an iterative approach"""
        if memo is None:
            memo = {}

        obj_id = id(self)
        if obj_id in memo:
            return memo[obj_id]

        root = self.__class__._copy_node(self)
        nodes_stack = [(root, self, len(self.children))]

        while nodes_stack:
            # check the top node, any children left unvisited?
            top = nodes_stack[-1]
            new_top_node, old_top_node, unvisited_children = top

            if unvisited_children:
                nodes_stack[-1] = (new_top_node, old_top_node, unvisited_children - 1)
                old_child = old_top_node.children[-unvisited_children]
                new_child = self.__class__._copy_node(old_child)
                new_top_node.append(new_child)
                nodes_stack.append((new_child, old_child, len(old_child.children)))
            else:  # no unvisited children
                nodes_stack.pop()
        return root

    __deepcopy__ = deepcopy = copy

    # support for basic tree operations -- finding objects and moving in the
    # tree
    @property
    def parent(self) -> Self | None:
        """parent of this node"""
        return self._parent

    @parent.setter
    def parent(self, parent: Self | None) -> None:
        """parent of this node"""
        if self._parent is not None:
            self._parent.remove_node(self)
        self._parent = parent
        if parent is not None and self not in parent.children:
            parent.children.append(self)

    def index_in_parent(self: Self) -> int:
        """Returns index of self in parent."""
        if self._parent is None:
            msg = "Node has no parent."
            raise TreeError(msg)
        return self._parent.children.index(self)

    def is_tip(self) -> bool:
        """Returns True if the current node is a tip, i.e. has no children."""
        return not self.children

    def is_root(self) -> bool:
        """Returns True if the current is a root, i.e. has no parent."""
        return self._parent is None

    def levelorder(self, include_self: bool = True) -> Generator[Self]:
        """Performs levelorder iteration over tree"""
        queue = [self]
        while queue:
            curr = queue.pop(0)
            if include_self or (curr is not self):
                yield curr
            if curr.children:
                queue.extend(curr.children)

    def preorder(self, include_self: bool = True) -> Generator[Self]:
        """Performs preorder iteration over tree."""
        stack = [self]

        while stack:
            node = stack.pop()
            if include_self or node is not self:
                yield node

            # the stack is last-in-first-out, so we add children
            # in reverse order so they're processed left-to-right
            if node.children:
                stack.extend(node.children[::-1])

    def postorder(self, include_self: bool = True) -> Generator[Self]:
        """performs postorder iteration over tree"""
        stack = [(self, False)]

        while stack:
            node, children_done = stack.pop()
            if children_done:
                if include_self or node is not self:
                    yield node
            else:
                # children still need to be processed
                stack.append((node, True))

                # the stack is last-in-first-out, so we add children
                # in reverse order so they're processed left-to-right
                if node.children:
                    stack.extend((child, False) for child in node.children[::-1])

    def pre_and_postorder(self, include_self: bool = True) -> Generator[Self]:
        """Performs iteration over tree, visiting node before and after."""
        yield from self.preorder(include_self=include_self)
        yield from self.postorder(include_self=include_self)

    def ancestors(self) -> list[Self]:
        """Returns all ancestors back to the root."""
        result: list[Self] = []
        curr = self._parent
        while curr is not None:
            result.append(curr)
            curr = curr._parent
        return result

    def get_root(self) -> Self:
        """Returns root of the tree self is in."""
        curr = self
        while curr._parent is not None:
            curr = curr._parent
        return curr

    def rooted(self, edge_name: str) -> Self:
        """Returns a new tree with split at edge_name

        Parameters
        ----------
        edge_name
            name of the edge to split at. The length of edge_name will be
            halved. The new tree will have two children.
        """
        tree = self.deepcopy()
        if not self.is_root():
            msg = (
                f"cannot apply from non-root node {self.name!r}, "
                "use self.get_root() first"
            )
            raise TreeError(msg)

        if edge_name == "root":
            return tree

        tree.source = None
        node = tree.get_node_matching_name(edge_name)
        is_tip = node.is_tip()
        has_length = any(node.length is not None for node in self.preorder())
        # we put tips on the right
        right_name = edge_name if is_tip else f"{edge_name}-R"
        left_name = f"{edge_name}-root" if is_tip else f"{edge_name}-L"
        length = (node.length or 0) / 2
        parent = cast("Self", node.parent)
        parent.children.remove(node)
        node.parent = None
        left = node.unrooted_deepcopy()
        right = parent.unrooted_deepcopy()
        if is_tip and left.is_tip():
            left.name = right_name
            right.name = left_name
        else:
            left.name = left_name
            right.name = right_name

        if has_length:
            left.length = length
            right.length = length

        result = self.__class__(name="root", children=[left, right])
        result.source = self.source
        result.prune()
        return result

    def isroot(self) -> bool:
        """Returns True if root of a tree, i.e. no parent."""
        return self.is_root()

    def siblings(self) -> list[Self]:
        """Returns all nodes that are children of the same parent as self.

        Note: excludes self from the list. Dynamically calculated.
        """
        if self._parent is None:
            return []
        result = self._parent.children[:]
        result.remove(self)
        return result

    def iter_tips(self, include_self: bool = False) -> Generator[Self]:
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

    def tips(self, include_self: bool = False) -> list[Self]:
        """Returns tips descended from self, [] if self is a tip."""
        return list(self.iter_tips(include_self=include_self))

    def iter_nontips(self, include_self: bool = False) -> Generator[Self]:
        """Iterates over nontips descended from self

        Parameters
        ----------
        include_self
            if True (default is False), will return the current
            node as part of the list of nontips if it is a nontip.
        """
        for n in self.preorder(include_self=include_self):
            if n.children:
                yield n

    def nontips(self, include_self: bool = False) -> list[Self]:
        """Returns nontips descended from self."""
        return list(self.iter_nontips(include_self=include_self))

    def istip(self) -> bool:
        """Returns True if is tip, i.e. no children."""
        return not self.children

    def tip_children(self) -> list[Self]:
        """Returns direct children of self that are tips."""
        return [i for i in self.children if not i.children]

    def non_tip_children(self) -> list[Self]:
        """Returns direct children in self that have descendants."""
        return [i for i in self.children if i.children]

    def last_common_ancestor(self, other: Self) -> Self:
        """Finds last common ancestor of self and other, or None.

        Always tests by identity.
        """
        my_lineage = {id(node) for node in [self, *self.ancestors()]}
        curr: Self | None = other
        while curr is not None:
            if id(curr) in my_lineage:
                return curr
            curr = curr._parent
        msg = "No common ancestor found."
        raise TreeError(msg)

    def lowest_common_ancestor(self, tip_names: list[str]) -> Self:
        """Lowest common ancestor for a list of tipnames

        This should be around O(H sqrt(n)), where H is height and n is the
        number of tips passed in.
        """
        if len(tip_names) == 1:
            return self.get_node_matching_name(tip_names[0])

        tip_names_set: set[str] = set(tip_names)
        tips = [tip for tip in self.tips() if tip.name in tip_names_set]

        if len(tips) != len(tip_names_set):
            missing = tip_names_set - set(self.get_tip_names())
            msg = f"tipnames {missing} not present in self"
            raise ValueError(msg)

        # scrub tree
        for n in self.preorder(include_self=True):
            n.params.pop("black", None)

        for t in tips:
            prev = t
            curr = t.parent

            while curr and "black" not in curr.params:
                curr.params["black"] = [prev]
                prev = curr
                curr = curr.parent

            # increase black count, multiple children lead to here
            if curr:
                curr.params["black"].append(prev)

        curr = self
        while len(curr.params.get("black", [])) == 1:
            curr = curr.params.pop("black")[0]

        return curr

    lca = last_common_ancestor  # for convenience

    # support for more advanced tree operations

    def separation(self, other: Self | None) -> int:
        """Returns number of edges separating self and other."""
        # detect trivial case
        if self is other:
            return 0
        # otherwise, check the list of ancestors
        my_ancestors = dict.fromkeys(list(map(id, [self, *self.ancestors()])))
        count = 0
        while other is not None:
            if id(other) in my_ancestors:
                # need to figure out how many steps there were back from self
                curr: Self | None = self
                while curr is not None and curr is not other:
                    count += 1
                    curr = curr.parent
                return count
            count += 1
            other = other.parent
        msg = "Nodes do not belong to the same tree."
        raise TreeError(msg)

    def descendant_array(
        self, tip_list: list[str] | None = None
    ) -> tuple[npt.NDArray[numpy.bool], list[Self]]:
        """Returns numpy array with nodes in rows and descendants in columns.

        True indicates that the decendant is a descendant of that node
        False indicates that it is not

        Also returns a list of nodes in the same order as they are listed
        in the array.

        tip_list is a list of the names of the tips that will be considered,
        in the order they will appear as columns in the final array. Internal
        nodes will appear as rows in preorder traversal order.
        """

        # get a list of internal nodes
        node_list = [node for node in self.preorder() if node.children]
        node_list.sort()

        # get a list of tip names if one is not supplied
        if not tip_list:
            tip_list = self.get_tip_names()
            tip_list.sort()
        # make a blank array of the right dimensions to alter
        result = numpy.zeros([len(node_list), len(tip_list)], dtype=numpy.bool)
        # put 1 in the column for each child of each node
        for i, node in enumerate(node_list):
            children = [n.name for n in node.tips()]
            for j, dec in enumerate(tip_list):
                if dec in children:
                    result[i, j] = 1
        return result, node_list

    def _default_tree_constructor(
        self,
    ) -> Callable[[Self | None, PySeq[Self], dict[str, Any] | None], Self]:
        return cast(
            "Callable[[Self | None, PySeq[Self], dict[str, Any] | None], Self]",
            TreeBuilder(constructor=self.__class__).edge_from_edge,
        )

    def name_unnamed_nodes(self) -> None:
        """sets the Data property of unnamed nodes to an arbitrary value

        Internal nodes are often unnamed and so this function assigns a
        Internal nodes are often unnamed and so this function assigns a
        value for referencing."""
        # make a list of the names that are already in the tree
        names_in_use = [node.name for node in self.preorder() if node.name]
        # assign unique names to the Data property of nodes where Data = None
        name_index = 1
        for node in self.preorder():
            if not node.name:
                new_name = f"node{name_index!s}"
                # choose a new name if name is already in tree
                while new_name in names_in_use:
                    name_index += 1
                    new_name = f"node{name_index}"
                node.name = new_name
                names_in_use.append(new_name)
                name_index += 1

    def make_tree_array(
        self, dec_list: list[str] | None = None
    ) -> tuple[npt.NDArray[numpy.uint8], list[Self]]:
        """Makes an array with nodes in rows and descendants in columns.

        A value of 1 indicates that the decendant is a descendant of that node/
        A value of 0 indicates that it is not

        also returns a list of nodes in the same order as they are listed
        in the array"""
        # get a list of internal nodes
        node_list = [node for node in self.preorder() if node.children]
        node_list.sort()

        # get a list of tips() name if one is not supplied
        if not dec_list:
            dec_list = self.get_tip_names()
            dec_list.sort()
        # make a blank array of the right dimensions to alter
        result = numpy.zeros((len(node_list), len(dec_list)), dtype=numpy.uint8)
        # put 1 in the column for each child of each node
        for i, node in enumerate(node_list):
            children = [dec.name for dec in node.tips()]
            for j, dec in enumerate(dec_list):
                if dec in children:
                    result[i, j] = 1
        return result, node_list

    def remove_deleted(self, should_delete: Callable[[Self], bool]) -> None:
        """Removes all nodes where should_delete tests true.

        Internal nodes that have no children as a result of removing deleted
        are also removed.
        """
        # Traverse tree
        for node in self.postorder():
            # if node is to be deleted
            if should_delete(node):
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

    def prune(
        self: Self,
        keep_root: bool = False,
        params_merge_callback: Callable[
            [dict[str, Any], dict[str, Any]], dict[str, Any]
        ]
        | None = None,
    ) -> None:
        """removes nodes with one child

        Parameters
        ----------
        keep_root
            If True, a root with a single child is retained.
        params_merge_callback
            How to merge two params dicts when pruning.
            The first argument is the parent node's params,
            the second argument is the child node's params.
            It should return the new params dictionary.

        Notes
        -----
        Mutates the tree in-place. Internal nodes with only one child will be
        merged (except as specified by keep_root).
        """
        while True:
            nodes_to_remove = [
                n
                for n in self.iter_nontips()
                if n.parent is not None and len(n.children) == 1
            ]
            if not nodes_to_remove:
                break

            for node in nodes_to_remove:
                curr_parent = node.parent
                child = node.children[0]
                node.parent = None
                child.parent = curr_parent

                if node.length is not None:
                    child.length = (child.length or 0) + node.length

                # For support, we keep the child's support. Do nothing.

                if params_merge_callback is not None:
                    child.params = params_merge_callback(node.params, child.params)

        # root having one child is edge case
        if not keep_root and len(self.children) == 1:
            child = self.children[0]

            grand_children = list(child.children)

            # Ignore the child's length as a root can't have a length
            # Unclear what should be done about support as well, keep current support.

            if params_merge_callback is not None:
                self.params = params_merge_callback(self.params, child.params)

            self.remove_node(child)
            for grand_child in grand_children:
                grand_child.parent = self

    def same_shape(self, other: Self) -> bool:
        """Ignores lengths and order, so trees should be sorted first"""
        if len(self.children) != len(other.children):
            return False
        if self.children:
            for self_child, other_child in zip(
                self.children,
                other.children,
                strict=False,
            ):
                if not self_child.same_shape(other_child):
                    return False
            return True
        return self.name == other.name

    def to_rich_dict(self) -> dict[str, Any]:
        """returns {'newick': with node names,
        'edge_attributes': {'tip1': {'length': ...}, ...}}"""
        newick = self.get_newick(
            with_node_names=True,
            semicolon=False,
            escape_name=False,
            with_root_name=True,
        )
        attr = {}
        length_and_support = {}
        for edge in self.get_edge_vector(include_root=True):
            attr[edge.name] = edge.params.copy()
            length_and_support[edge.name] = {
                "length": edge.length,
                "support": edge.support,
            }
        return {
            "newick": newick,
            "edge_attributes": attr,
            "length_and_support": length_and_support,
            "type": get_object_provenance(self),
            "version": __version__,
        }

    def to_json(self) -> str:
        """returns json formatted string {'newick': with edges and distances, 'edge_attributes': }"""
        return json.dumps(self.to_rich_dict())

    def get_newick(
        self,
        with_distances: bool = False,
        semicolon: bool = True,
        escape_name: bool = True,
        with_node_names: bool = False,
        with_root_name: bool = False,
    ) -> str:
        """Return the newick string of node and its descendents

        Parameters
        ----------
        with_distances
            include value of node length attribute if present.
        semicolon
            end tree string with a semicolon
        escape_name
            if any of these characters []'"() are within the
            nodes name, wrap the name in single quotes
        with_node_names
            includes internal node names
        with_root_name
            if True and with_node_names, the root node will have
            its name included
        """
        # Stack contains tuples of (tree node, visit flag)
        stack = [(self, False)]
        node_results: dict[int, str] = {}  # results cache

        while stack:
            node, visited = stack.pop()

            if not visited:
                # First visit - push back for processing after children
                stack.append((node, True))
                # add each child to the stack
                stack.extend((child, False) for child in node.children)
            else:
                # children have been seen once
                node_name = _format_node_name(
                    node,
                    with_node_names=with_node_names,
                    escape_name=escape_name,
                    with_distances=with_distances,
                    with_root_name=with_root_name,
                )

                # for tips with parent, the typical case
                if node.is_tip() and node.parent:
                    node_results[id(node)] = node_name
                    continue

                # collecting children
                # Build result for this node
                if children_newick := [
                    node_results[id(child)] for child in node.children
                ]:
                    result = f"({','.join(children_newick)}){node_name}"
                else:
                    result = node_name

                node_results[id(node)] = result

        # final result
        final_result = node_results[id(self)]

        if self.is_root() and semicolon:
            final_result = f"{final_result};"

        return final_result

    def remove_node(self, target: Self) -> bool:
        """Removes node by identity instead of value.

        Returns True if node was present, False otherwise.
        """
        for i, curr_node in enumerate(self.children):
            if curr_node is target:
                del self[i]
                return True
        return False

    def get_edge_names(
        self,
        tip_name_1: str,
        tip_name_2: str,
        clade: bool = True,
        stem: bool = False,
        outgroup_name: str | None = None,
    ) -> list[str]:
        """Return the list of stem and/or sub tree (clade) edge name(s).
        This is done by finding the common intersection, and then getting
        the list of names. If the clade traverses the root, then use the
        outgroup_name argument to ensure valid specification.

        Parameters
        ----------
        tip_name_1/2
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
        root = self
        if outgroup_name is not None:
            outgroup = self.get_node_matching_name(outgroup_name)
            if not outgroup.is_tip():
                msg = f"Outgroup ({outgroup_name!r}) is not a tip"
                raise TreeError(msg)
            root = outgroup.unrooted_deepcopy()

        join_edge = root.get_connecting_node(tip_name_1, tip_name_2)

        edge_names: list[str] = []

        if stem:
            if join_edge.isroot():
                msg = f"LCA({tip_name_1},{tip_name_2}) is the root and so has no stem"
                raise TreeError(
                    msg,
                )
            edge_names.append(join_edge.name)

        if clade:
            # get the list of names contained by join_edge
            for child in join_edge.children:
                branch_names = child.get_node_names(include_self=True)
                edge_names.extend(branch_names)

        return edge_names

    def get_neighbours_except(self, parent: Self | None = None) -> list[Self]:
        # For walking the tree as if it was unrooted.
        return [
            c
            for c in (*self.children, self.parent)
            if c is not None and c is not parent
        ]

    def get_sub_tree(
        self,
        names: Iterable[str],
        ignore_missing: bool = False,
        tips_only: bool = False,
        as_rooted: bool = False,
    ) -> Self:
        """A new instance of a sub tree that contains all the otus that are
        listed in name_list.

        Parameters
        ----------
        ignore_missing
            if False, get_sub_tree will raise a ValueError if
            name_list contains names that aren't nodes in the tree
        tips_only
            only tip names matching name_list are allowed
        as_rooted
            if True, the resulting subtree root will be as resolved. Otherwise,
            the subtree is coerced to have the same number of children as self.
        """
        # find all the selected nodes
        allowed = set(names)
        old_nodes: dict[int, Self] = {}
        found: set[str] = set()
        for old_node in self.preorder(include_self=True):
            if old_node.name not in allowed:
                continue

            found.add(old_node.name)
            old_nodes[id(old_node)] = old_node
            # find all nodes connecting required nodes to root,
            # skipping if already present
            parent = old_node.parent
            while parent is not None and (parent_id := id(parent)) not in old_nodes:
                old_nodes[parent_id] = parent
                parent = parent.parent

            if not tips_only and not old_node.is_tip():
                # add all descendant nodes too
                for n in old_node.preorder():
                    old_nodes[id(n)] = n

        if found != allowed and not ignore_missing:
            msg = f"edges {allowed - found} not found in tree"
            raise ValueError(msg)

        # make new nodes and also map old id's to new id's
        make_node = self.__class__
        self_2_new: dict[int, int] = {}
        new_nodes: dict[int, Self] = {}
        for self_id, old_node in old_nodes.items():
            new_node = make_node(
                old_node.name,
                params=old_node.params.copy(),
                length=old_node.length,
                support=old_node.support,
            )
            new_nodes[id(new_node)] = new_node
            self_2_new[self_id] = id(new_node)

        # connect the nodes
        for self_id, old_node in old_nodes.items():
            if old_node.parent is None:
                continue

            new_node = new_nodes[self_2_new[self_id]]
            new_parent_id = self_2_new[id(old_node.parent)]
            # the following assignment also adds the new_node as
            # a child to parent
            new_node.parent = new_nodes[new_parent_id]

        result_root = new_nodes[self_2_new[id(self)]]
        result_root.prune()
        if as_rooted or len(self.children) == len(result_root.children):
            result_root.name = "root"
            return result_root

        if len(self.children) > 2:
            result_root = result_root.unrooted()
        else:
            # we pick an arbitrary child to root at
            child = result_root.children[0]
            child.name = (
                child.name if child.name else "new-root"
            )  # this is a magic value, which is not good
            result_root = result_root.rooted(child.name)
        result_root.name = "root"
        return result_root

    def _edgecount(self, parent: Self, cache: dict[tuple[int, int], int]) -> int:
        """The number of edges beyond 'parent' in the direction of 'self',
        unrooted"""
        neighbours = self.get_neighbours_except(parent)
        key = (id(parent), id(self))
        if key not in cache:
            cache[key] = 1 + sum(
                [child._edgecount(self, cache) for child in neighbours],
            )
        return cache[key]

    def _imbalance(
        self, parent: Self | None, cache: dict[tuple[int, int], int]
    ) -> tuple[int, int, Self]:
        """The edge count from here, (except via 'parent'), divided into that
        from the heaviest neighbour, and that from the rest of them.  'cache'
        should be a dictionary that can be shared by calls to self.edgecount,
        it stores the edgecount for each node (from self) without having to
        put it on the tree itself."""
        max_weight = 0
        total_weight = 0
        biggest_branch = self
        for child in self.get_neighbours_except(parent):
            weight = child._edgecount(self, cache)
            total_weight += weight
            if weight > max_weight:
                max_weight = weight
                biggest_branch = child
        return max_weight, total_weight - max_weight, biggest_branch

    def sorted(self, sort_order: list[str] | None = None) -> Self:
        """An equivalent tree with tips in sort_order.

        Notes
        -----
        If sort_order is not specified then alphabetical order is used.
        At each node starting from root, the algorithm will try to put
        the descendant which contains the smallest index tip on the left.
        """
        sort_order = sort_order or []
        tip_names = self.get_tip_names()
        tip_names.sort()
        full_sort_order = sort_order + tip_names
        score_map = {name: i for i, name in enumerate(full_sort_order)}

        constructor = self._default_tree_constructor()

        scores: dict[PhyloNode, int | None] = {}
        rebuilt: dict[PhyloNode, PhyloNode] = {}

        infinity = float("inf")
        for node in self.postorder():
            if node.is_tip():
                score = score_map[node.name]
                tree = node.deepcopy()
            else:
                child_info = [(scores[ch], rebuilt[ch]) for ch in node.children]
                # Sort children by score, None is treated as +infinity
                child_info.sort(key=lambda x: (infinity if x[0] is None else x[0]))
                children = tuple(child for _, child in child_info)
                tree = constructor(node, children, None)
                non_null = [s for s, _ in child_info if s is not None]
                score = non_null[0] if non_null else None
            scores[node] = score
            rebuilt[node] = tree

        return rebuilt[self]

    def ladderise(self) -> Self:
        """Return an equivalent tree nodes using a ladderise sort.

        Notes
        -----
        Children are ordered by their number of descendant tips
        with ties broken by alphabetical sort of node names.
        """
        num_tips = {}
        ordered_names_map = {}

        for node in self.postorder():
            if node.is_tip():
                num_tips[node] = 1
                ordered_names_map[node] = [node.name]
            else:
                ordered_kids = sorted(
                    node.children, key=lambda c: (num_tips[c], ordered_names_map[c][0])
                )

                num_tips[node] = sum(num_tips[k] for k in ordered_kids)

                names = []
                for child in ordered_kids:
                    names.extend(ordered_names_map[child])
                ordered_names_map[node] = names

        ordered_names = ordered_names_map[self]
        return self.sorted(sort_order=ordered_names)

    ladderize = ladderise  # a synonym with US spelling

    def _ascii_art(
        self,
        char1: str = "-",
        show_internal: bool = True,
        compact: bool = False,
    ) -> tuple[list[str], int]:
        length = 10
        pad = " " * length
        pa = " " * (length - 1)
        namestr = self.name
        if self.children:
            mids: list[int] = []
            result: list[str] = []
            for c in self.children:
                if c is self.children[0]:
                    char2 = "/"
                elif c is self.children[-1]:
                    char2 = "\\"
                else:
                    char2 = "-"
                clines, mid = c._ascii_art(char2, show_internal, compact)
                mids.append(mid + len(result))
                result.extend(clines)
                if not compact:
                    result.append("")
            if not compact:
                result.pop()
            lo, hi, end = mids[0], mids[-1], len(result)
            prefixes = (
                [pad] * (lo + 1) + [pa + "|"] * (hi - lo - 1) + [pad] * (end - hi)
            )
            mid = (lo + hi) // 2
            prefixes[mid] = char1 + "-" * (length - 2) + prefixes[mid][-1]
            result = [pre + res for pre, res in zip(prefixes, result, strict=False)]
            if show_internal:
                stem = result[mid]
                result[mid] = stem[0] + namestr + stem[len(namestr) + 1 :]
            return result, mid
        return [char1 + "-" + namestr], 0

    def ascii_art(self, show_internal: bool = True, compact: bool = False) -> str:
        """Returns a string containing an ascii drawing of the tree.

        Parameters
        ----------
        show_internal
            includes internal edge names.
        compact
            use exactly one line per tip.

        """
        lines, _ = self._ascii_art(show_internal=show_internal, compact=compact)
        return "\n".join(lines)

    def write(
        self,
        filename: str | os.PathLike[str],
        with_distances: bool = True,
        format_name: str | None = None,
    ) -> None:
        """Save the tree to filename

        Parameters
        ----------
        filename
            path to write the tree to.
        with_distances
            whether branch lengths are included in string.
        format_name
            default is newick, json is alternate. Argument overrides
            the filename suffix. All attributes are saved in the xml format.
            Value overrides the file name suffix.

        Notes
        -----
        Only the cogent3 json and newick tree formats are supported.

        """
        file_format, _ = get_format_suffixes(filename)
        format_name = format_name or file_format
        if format_name == "json":
            with atomic_write(filename, mode="wt") as f:
                f.write(self.to_json())
            return

        data = self.get_newick(with_distances=with_distances)

        with atomic_write(filename, mode="wt") as outf:
            outf.writelines(data)

    def get_node_names(
        self, include_self: bool = True, tips_only: bool = False
    ) -> list[str]:
        """Return a list of edges from this edge - may or may not include self.
        This node (or first connection) will be the first, and then they will
        be listed in the natural traverse order.

        Parameters
        ----------
        include_self : bool
            excludes self.name from the result

        tips_only : bool
            only tips returned
        """
        if tips_only:
            nodes = self.tips(include_self=include_self)
        else:
            nodes = list(self.preorder(include_self=include_self))
        return [node.name for node in nodes]

    def get_tip_names(self, include_self: bool = True) -> list[str]:
        """return the list of the names of all tips contained by this edge"""
        node_names = self.get_node_names(include_self=include_self, tips_only=True)

        if "" in node_names:
            msg = "Tree contains unnamed nodes."
            raise TreeError(msg)

        return node_names

    def get_edge_vector(self, include_root: bool = True) -> list[Self]:
        """Collect the list of edges in postfix order

        Parameters
        ----------
        include_root
            specifies whether root edge included

        """
        return list(self.postorder(include_self=include_root))

    def get_node_matching_name(self, name: str) -> Self:
        """find the edge with the name

        Raises
        ------
        TreeError if no edge with the name is found
        """
        for node in self.preorder(include_self=True):
            if node.name == name:
                break
        else:
            msg = f"No node named '{name}' in {self.get_tip_names()}"
            raise TreeError(msg)
        return node

    def get_connecting_node(self, name1: str, name2: str) -> Self:
        """Finds the last common ancestor of the two named edges."""
        edge1 = self.get_node_matching_name(name1)
        edge2 = self.get_node_matching_name(name2)
        return edge1.last_common_ancestor(edge2)

    def get_connecting_edges(self, name1: str, name2: str) -> list[Self]:
        """returns a list of edges connecting two nodes.

        If both are tips, the LCA is excluded from the result."""
        edge1 = self.get_node_matching_name(name1)
        edge2 = self.get_node_matching_name(name2)
        include_parent = not (edge1.istip() and edge2.istip())

        lca = self.get_connecting_node(name1, name2)
        node_path = [edge1]
        node_path.extend(edge1.ancestors())
        # remove nodes deeper than the LCA
        lca_ind = node_path.index(lca)
        node_path = node_path[: lca_ind + 1]
        # remove LCA and deeper nodes from anc list of other
        anc2 = edge2.ancestors()
        lca_ind = anc2.index(lca)
        anc2 = anc2[:lca_ind]
        anc2.reverse()
        node_path.extend(anc2)
        node_path.append(edge2)
        if not include_parent:
            node_path.remove(lca)
        return node_path

    def get_param_value(self, param: str, edge: str) -> Any:  # noqa: ANN401
        """returns the parameter value for named edge"""
        return self.get_node_matching_name(edge).params[param]

    def set_param_value(self, param: str, edge: str, value: Any) -> None:  # noqa: ANN401
        """set's the value for param at named edge"""
        self.get_node_matching_name(edge).params[param] = value

    def reassign_names(
        self, mapping: dict[str, str], nodes: list[Self] | None = None
    ) -> None:
        """Reassigns node names based on a mapping dict

        mapping : dict, old_name -> new_name
        nodes : specific nodes for renaming (such as just tips, etc...)
        """
        if nodes is None:
            nodes = list(self.preorder())

        for n in nodes:
            if n.name in mapping:
                n.name = mapping[n.name]

    def multifurcating(
        self,
        num: int,
        eps: float | None = None,
        name_unnamed: bool = False,
    ) -> Self:
        """return a new tree with every node having num or few children

        Parameters
        ----------
        num : int
            the number of children a node can have max
        eps : float
            default branch length to set if self or constructor is of
            PhyloNode type
            a PhyloNode or subclass constructor. If None, uses self
        name_unnamed : bool
            names unnamed nodes
        """
        if num < 2:
            msg = "Minimum number of children must be >= 2"
            raise TreeError(msg)

        if eps is None:
            eps = 0.0

        constructor = self.__class__

        new_tree = self.copy()

        for n in new_tree.preorder(include_self=True):
            while len(n.children) > num:
                new_node = constructor("", children=n.children[-num:])

                if new_node[0].length is not None:
                    new_node.length = eps

                n.append(new_node)

        if name_unnamed:
            alpha = "abcdefghijklmnopqrstuvwxyz"
            alpha += alpha.upper()
            base = "AUTOGENERATED_NAME_%s"

            # scale the random names by tree size
            s = int(numpy.ceil(numpy.log(len(new_tree.tips()))))

            for n in new_tree.nontips():
                if not n.name:
                    n.name = base % "".join([random.choice(alpha) for _ in range(s)])

        return new_tree

    def bifurcating(
        self,
        eps: float | None = None,
        name_unnamed: bool = False,
    ) -> Self:
        """Wrap multifurcating with a num of 2"""
        return self.multifurcating(2, eps, name_unnamed)

    def get_nodes_dict(self) -> dict[str, Self]:
        """Returns a dict keyed by node name, value is node

        Will raise TreeError if non-unique names are encountered
        """
        res: dict[str, Self] = {}

        for n in self.preorder():
            if n.name in res:
                msg = "get_nodes_dict requires unique node names"
                raise TreeError(msg)
            res[n.name] = n

        return res

    def subset(self) -> frozenset[str]:
        """Returns set of names that descend from specified node"""
        return frozenset(self.get_tip_names(include_self=False))

    def subsets(self) -> frozenset[frozenset[str]]:
        """Returns all sets of names that come from specified node and its kids"""
        sets: list[frozenset[str]] = []
        for node in self.postorder(include_self=False):
            if not node.children:
                node.params["leaf_set"] = frozenset([node.name])
            else:
                leaf_set: frozenset[str] = reduce(
                    or_, [c.params.pop("leaf_set") for c in node.children]
                )
                if len(leaf_set) > 1:
                    sets.append(leaf_set)
                node.params["leaf_set"] = leaf_set

        # clean up params entry in children of self
        for child in self.children:
            child.params.pop("leaf_set", None)
        return frozenset(sets)

    def compare_by_subsets(
        self, other: Self, exclude_absent_taxa: bool = False
    ) -> float:
        """Returns fraction of overlapping subsets where self and other differ.

        Other is expected to be a tree object compatible with PhyloNode.

        Note: names present in only one of the two trees will count as
        mismatches: if you don't want this behavior, strip out the non-matching
        tips first.
        """
        self_sets, other_sets = self.subsets(), other.subsets()
        if exclude_absent_taxa:
            in_both = self.subset() & other.subset()
            self_sets = frozenset(
                intersection for i in self_sets if len(intersection := i & in_both) > 1
            )
            other_sets = frozenset(
                intersection for i in other_sets if len(intersection := i & in_both) > 1
            )
        total_subsets = len(self_sets) + len(other_sets)
        intersection_length = len(self_sets & other_sets)
        if not total_subsets:  # no common subsets after filtering, so max dist
            return 1
        return 1 - 2 * intersection_length / float(total_subsets)

    def tip_to_tip_distances(
        self, names: PySeqStr | None = None, default_length: float | None = None
    ) -> DistanceMatrix:
        """Returns distance matrix between all pairs of tips, and a tip order"""
        from cogent3.evolve.fast_distance import DistanceMatrix

        if names is not None:
            subtree = self.get_sub_tree(names)
            return subtree.tip_to_tip_distances(
                default_length=default_length,
            )
        default_length = (
            1 if all(node.length is None for node in self.preorder()) else 0
        )
        tips = list(self.tips())

        # For each tip, build path to root with cumulative distances
        paths: dict[
            str, list[tuple[Self, float]]
        ] = {}  # tip name -> list of (node, cumulative distance)
        for tip in tips:
            path: list[tuple[Self, float]] = []
            current: Self | None = tip  # type: ignore[assignment]

            dist = 0.0
            while current is not None:
                path.append((current, dist))
                length = (
                    current.length if current.length is not None else default_length
                )
                dist += length
                current = current.parent
            paths[tip.name] = path  # path from tip to root

        num_tips = len(tips)
        dists = numpy.zeros((num_tips, num_tips), float)
        for i, j in combinations(range(num_tips), 2):
            tip1 = tips[i]
            tip2 = tips[j]
            path1 = {id(node): (node, dist) for node, dist in paths[tip1.name]}
            path2 = {id(node): (node, dist) for node, dist in paths[tip2.name]}
            common = path1.keys() & path2.keys()

            if not common:
                msg = f"No common ancestor for {tip1.name} and {tip2.name}"
                raise ValueError(msg)

            # Find least common ancestor (node with max total depth)
            lca = min(common, key=lambda n: path1[n][1])
            total_dist = path1[lca][1] + path2[lca][1]
            dists[i, j] = dists[j, i] = total_dist

        return DistanceMatrix.from_array_names(dists, self.get_tip_names())

    def get_figure(
        self,
        style: Literal["square", "circular", "angular", "radial"] = "square",
        **kwargs: Any,
    ) -> Dendrogram:
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

        style = cast(
            "Literal['square', 'circular', 'angular', 'radial']", style.lower()
        )
        types = ("square", "circular", "angular", "radial")
        if style not in types:
            msg = f"{style} not in supported types {types}"
            raise ValueError(msg)

        return Dendrogram(self, style=style, **kwargs)

    def balanced(self) -> Self:
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
        cache: dict[tuple[int, int], int] = {}
        while 1:
            max_weight, remaining_weight, next_edge = edge._imbalance(
                last_edge,
                cache,
            )
            known_weight += remaining_weight
            if max_weight <= known_weight + 2:
                break
            last_edge = edge
            edge = next_edge
            known_weight += 1
        return edge.unrooted_deepcopy()

    def same_topology(self, other: Self) -> bool:
        """Tests whether two trees have the same topology."""
        tip_names = self.get_tip_names()
        root_at = tip_names[0]
        me = self.rooted(root_at).sorted(tip_names)
        them = other.rooted(root_at).sorted(tip_names)
        return self is other or me.same_shape(them)

    def unrooted_deepcopy(
        self,
        parent: Self | None = None,
    ) -> Self:
        """
        Returns a deepcopy of the tree using unrooted traversal.

        Each node is treated as connected to its parent and children.
        The resulting tree may contain unary internal nodes, which can
        be cleaned up using `prune()` afterward.
        """
        constructor = self._default_tree_constructor()

        # node_map maps id(original_node) -> new_node
        node_map: dict[int, Self] = {}
        # stack is last in first out
        # stack stores (original_node, parent_we_came_from, state)
        # False state is the first visit, discover neighbors
        # True state is the second visit, construct new node
        stack = [(self, parent, False)]
        while stack:
            node, parent_node, state = stack.pop()

            if not state:
                # put the node, and then it's children on the stack
                stack.append((node, parent_node, True))
                stack.extend(
                    (neigh, node, False)
                    for neigh in node.get_neighbours_except(parent_node)
                )
            else:
                # children are created and in node_map prior to their parents
                # being visited
                children = [
                    node_map[id(neigh)]
                    for neigh in node.get_neighbours_except(parent_node)
                ]

                if parent_node is None:
                    edge = None
                elif parent_node.parent is node:
                    edge = parent_node
                else:
                    edge = node

                new_node = constructor(edge, tuple(children), None)
                node_map[id(node)] = new_node
                if parent_node is None:
                    new_node.name = "root"

        new_root = node_map[id(self)]
        new_root.prune(keep_root=True)
        return new_root

    def unrooted(self) -> Self:
        """A tree with at least 3 children at the root."""
        constructor = self._default_tree_constructor()
        need_to_expand = len(self.children) < 3
        new_children: list[Self] = []
        for oldnode in self.children:
            if oldnode.children and need_to_expand:
                for sib in oldnode.children:
                    new_sib = sib.deepcopy()
                    if new_sib.length is not None and oldnode.length is not None:
                        new_sib.length += oldnode.length
                    new_children.append(new_sib)
                need_to_expand = False
            else:
                new_children.append(oldnode.deepcopy())
        return constructor(self, new_children, None)

    def rooted_at(self, edge_name: str) -> Self:
        """Return a new tree rooted at the provided node.

        Usage:
            This can be useful for drawing unrooted trees with an orientation
            that reflects knowledge of the true root location.
        """
        newroot = self.get_node_matching_name(edge_name)
        if not newroot.children:
            msg = f"Can't use a tip ({edge_name!r}) as the root"
            raise TreeError(msg)
        return newroot.unrooted_deepcopy()

    def rooted_with_tip(self, outgroup_name: str) -> Self:
        """A new tree with the named tip as one of the root's children"""
        tip = self.get_node_matching_name(outgroup_name)
        parent = cast("Self", tip.parent)
        return parent.unrooted_deepcopy()

    def tree_distance(self, other: PhyloNode, method: str | None = None) -> int:
        """Return the specified tree distance between this and another tree.

        Defaults to the Lin-Rajan-Moret distance on unrooted trees.
        Defaults to the Matching Cluster distance on rooted trees.

        Parameters
        ----------
        other: PhyloNode
            The other tree to calculate the distance between.
        method: str | None
            The tree distance metric to use.

            Options are:
            "rooted_robinson_foulds": The Robinson-Foulds distance for rooted trees.
            "unrooted_robinson_foulds": The Robinson-Foulds distance for unrooted trees.
            "matching_cluster": The Matching Cluster distance for rooted trees.
            "lin_rajan_moret": The Lin-Rajan-Moret distance for unrooted trees.
            "rrf": An alias for rooted_robinson_foulds.
            "urf": An alias for unrooted_robinson_foulds.
            "mc": An alias for matching_cluster.
            "lrm": An alias for lin_rajan_moret.
            "rf": The unrooted/rooted Robinson-Foulds distance for unrooted/rooted trees.
            "matching": The Lin-Rajan-Moret/Matching Cluster distance for unrooted/rooted trees.

            Default is "matching".

        Returns
        -------
        int
            the chosen distance between the two trees.

        Notes
        -----
        The Lin-Rajan-Moret distance [2]_ and Matching Cluster distance [1]_
        display superior statistical properties than the Robinson-Foulds
        distance [3]_ on unrooted and rooted trees respectively.

        References
        ----------
        .. [1] Bogdanowicz, D., & Giaro, K. (2013).
           On a matching distance between rooted phylogenetic trees.
           International Journal of Applied Mathematics and Computer Science, 23(3), 669-684.
        .. [2] Lin et al. 2012
           A Metric for Phylogenetic Trees Based on Matching
           IEEE/ACM Transactions on Computational Biology and Bioinformatics
           vol. 9, no. 4, pp. 1014-1022, July-Aug. 2012
        .. [3] Robinson, David F., and Leslie R. Foulds.
           Comparison of phylogenetic trees.
           Mathematical biosciences 53.1-2 (1981): 131-147.
        """

        if method is None:
            method = "matching"

        is_rooted = len(self) == 2
        if (is_rooted and len(other) != 2) or (not is_rooted and len(other) == 2):
            msg = "Both trees must be rooted or both trees must be unrooted."
            raise ValueError(
                msg,
            )

        return get_tree_distance_measure(method, is_rooted)(self, other)

    def lin_rajan_moret(self, tree2: Self) -> int:
        """return the lin-rajan-moret distance between trees

        float
            the Lin-Rajan-Moret distance

        Notes
        -----
        This is a distance measure that exhibits superior statistical
        properties compared to Robinson-Foulds. It can only be applied to
        unrooted trees.

        see: Lin et al. 2012
        A Metric for Phylogenetic Trees Based on Matching
        IEEE/ACM Transactions on Computational Biology and Bioinformatics
        vol. 9, no. 4, pp. 1014-1022, July-Aug. 2012
        """
        from cogent3.phylo.tree_distance import lin_rajan_moret

        return lin_rajan_moret(self, tree2)

    def child_parent_map(self) -> dict[str, str]:
        """return dict of {<child name>: <parent name>, ...}"""
        return {
            e.name: cast("Self", e.parent).name
            for e in self.postorder(include_self=False)
        }

    def distance(self, other: Self) -> float:
        """Returns branch length between self and other."""
        # never any length between self and other
        if self is other:
            return 0
        # otherwise, find self's ancestors and find the first ancestor of
        # other that is in the list
        self_anc = self.ancestors()
        self_anc_dict = {id(n): n for n in self_anc}
        self_anc_dict[id(self)] = self

        count = 0.0

        other_chain: Self | None = other
        while other_chain is not None:
            if id(other_chain) in self_anc_dict:
                # found the first shared ancestor -- need to sum other branch
                curr = self
                while curr is not other_chain:
                    if curr.length:
                        count += curr.length
                    curr = cast("Self", curr._parent)
                return count
            if other_chain.length:
                count += other_chain.length
            other_chain = other_chain._parent
        msg = "The other node is not in the same tree."
        raise TreeError(msg)

    def total_descending_branch_length(self) -> float:
        """Returns total descending branch length from self"""
        return sum(
            n.length for n in self.preorder(include_self=False) if n.length is not None
        )

    def total_length(self) -> float:
        """returns the sum of all branch lengths in tree"""
        root = self.get_root()
        return root.total_descending_branch_length()

    def tips_within_distance(self, distance: float) -> list[Self]:
        """Returns tips within specified distance from self

        Branch lengths of None will be interpreted as 0
        """

        def get_distance(d1: float, d2: float | None) -> float:
            if d2 is None:
                return d1
            return d1 + d2

        to_process = [(self, 0.0)]
        tips_to_save: list[Self] = []

        seen = {id(self)}
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

    def root_at_midpoint(self) -> Self:
        """return a new tree rooted at midpoint of the two tips farthest apart

        this fn doesn't preserve the internal node naming or structure,
        but does keep tip to tip distances correct.  uses unrooted_deepcopy()
        """
        dmat = self.tip_to_tip_distances()
        a, b = dmat.max_pair()
        max_dist = dmat[a, b]
        if max_dist <= 0.0:
            msg = f"{max_dist=} must be > 0"
            raise TreeError(msg)

        mid_point = max_dist / 2.0
        path_nodes = self.get_connecting_edges(a, b)
        cumsum = 0.0
        has_length = any(node.length is not None for node in self.preorder())
        default_length = 0.0 if has_length else 1.0

        node = path_nodes[0]
        for node in path_nodes:
            length = (node.length or default_length) if has_length else default_length
            cumsum += length
            if cumsum >= mid_point:
                break

        parent: Self = cast("Self", node.parent)
        if parent.is_root() and len(parent.children) == 2:
            # already midpoint rooted, but adjust lengths from root
            _adjust_lengths_from_root(tip_name=a, mid_point=mid_point, tree=self)
            return self

        new_tree = self.rooted(node.name)
        _adjust_lengths_from_root(tip_name=a, mid_point=mid_point, tree=new_tree)
        return new_tree

    def get_max_tip_tip_distance(
        self,
    ) -> tuple[float, tuple[str, str], Self]:
        """Returns the max tip-to-tip distance between any pair of tips

        Returns
        -------
        dist, tip_names, internal_node
        """
        dmat = self.tip_to_tip_distances()
        a, b = dmat.max_pair()
        dist = dmat[a, b]
        return float(dist), (a, b), self.get_connecting_node(a, b)

    def max_tip_tip_distance(self) -> tuple[float, tuple[str, str]]:
        """returns the max distance between any pair of tips

        Also returns the tip names  that it is between as a tuple"""
        dist, pair, _ = self.get_max_tip_tip_distance()
        return dist, pair

    @staticmethod
    def parse_token(token: str | None) -> tuple[str | None, float | None]:
        name, support = split_name_and_support(token)
        return name, support

    def tip_to_root_distances(
        self,
        names: list[str] | None = None,
        default_length: float = 1,
        *,
        node_length: bool = False,
    ) -> dict[str, float]:
        """returns the cumulative sum of lengths from each tip to the root

        Parameters
        ----------
        names
            list of tip names to calculate distances for, defaults to all
        default_length
            value to use for edges that no length value
        """
        tips = self.tips()
        if names is not None:
            tips = [t for t in tips if t.name in names]

        if not tips:
            msg = f"No tips matching in {names!r}"
            raise TreeError(msg)

        dists: dict[str, float] = {}
        for tip in tips:
            node = tip
            cum_sum = 0.0
            while node.parent is not None:
                if node_length or node.length is None:
                    cum_sum += default_length
                else:
                    cum_sum += node.length

                node = node.parent
            dists[tip.name] = cum_sum
        return dists


T = TypeVar("T", bound=PhyloNode)


def _adjust_lengths_from_root(
    *, tip_name: str, mid_point: float, tree: PhyloNode
) -> None:
    if len(tree.children) != 2:
        msg = "root node must have 2 children"
        raise TreeError(msg)

    to_tip, other = tree.children
    if tip_name not in to_tip.get_tip_names():
        to_tip, other = other, to_tip

    a_to_root = tree.tip_to_root_distances(names=[tip_name])
    delta = a_to_root[tip_name] - mid_point

    if to_tip.length is not None:
        to_tip.length -= delta
    if other.length is not None:
        other.length += delta


def split_name_and_support(name_field: str | None) -> tuple[str | None, float | None]:
    """Handle cases in the Newick format where an internal node name field
    contains a name or/and support value, like 'edge.98/100'.
    """
    # handle the case where the name field is None or empty string
    if not name_field:
        return None, None

    # if name_field is "24", treat it as support, returns (None, 24.0)
    with contextlib.suppress(ValueError):
        return None, float(name_field)

    # otherwise, split the name field into name and support
    name, *support = name_field.split("/")

    if len(support) == 1:
        try:
            support_value = float(support[0])
        except ValueError as e:
            msg = f"Support value at node: {name!r} should be int or float not {support[0]!r}."
            raise ValueError(
                msg,
            ) from e
    # handle case where mutiple '/' in the name field
    elif len(support) > 1:
        msg = f"Support value at node: {name!r} should be int or float not {'/'.join(support)!r}."
        raise ValueError(
            msg,
        )
    else:
        support_value = None

    return name, support_value


class TreeBuilder:
    # Some tree code which isn't needed once the tree is finished.
    # Mostly exists to give edges unique names
    # children must be created before their parents.

    def __init__(self, constructor: type[PhyloNode] = PhyloNode) -> None:
        self._used_names = {"edge": -1}
        self.PhyloNodeClass = constructor

    def _unique_name(self, name: str | None) -> str:
        # Unnamed edges become edge.0, edge.1 edge.2 ...
        # Other duplicates go mouse mouse.2 mouse.3 ...
        if not name:
            name = "edge"
        if name in self._used_names:
            self._used_names[name] += 1
            name += f".{self._used_names[name]!s}"
            # in case of names like 'edge.1.1'
            name = self._unique_name(name)
        else:
            self._used_names[name] = 1
        return name

    def _params_for_edge(self, edge: PhyloNode) -> dict[str, Any]:
        # default is just to keep it
        return edge.params

    def edge_from_edge(
        self,
        edge: PhyloNode | None,
        children: PySeq[PhyloNode],
        params: dict[str, Any] | None = None,
    ) -> PhyloNode:
        """Callback for tree-to-tree transforms like get_sub_tree"""
        if not isinstance(children, list):
            children = list(children)
        if edge is None:
            if params:
                msg = "No params allowed when edge is None."
                raise ValueError(msg)
            return self.create_edge(
                children,
                "root",
                {},
                None,
                None,
                name_loaded=False,
            )
        if params is None:
            params = self._params_for_edge(edge)
        return self.create_edge(
            children,
            edge.name,
            params,
            edge.length,
            edge.support,
            name_loaded=edge.name_loaded,
        )

    def create_edge(
        self,
        children: PySeq[PhyloNode] | None,
        name: str | None,
        params: dict[str, Any],
        length: float | None,
        support: float | None,
        name_loaded: bool = True,
    ) -> PhyloNode:
        """Callback for newick parser"""
        if children is None:
            children = []
        # split name and support for internal nodes
        elif children != []:
            name, new_support = self.PhyloNodeClass.parse_token(name)

            if new_support is not None:
                if support is not None and new_support != support:
                    msg = f"Got conflicting values for support. In name token '{name}': {new_support}. In constructor {support}."
                    raise ValueError(msg)
                support = new_support

        return self.PhyloNodeClass(
            name=self._unique_name(name),
            children=list(children),
            name_loaded=name_loaded and (name is not None),
            params=params,
            length=length,
            support=support,
        )


def make_tree(
    treestring: str | None = None,
    tip_names: list[str] | None = None,
    format_name: str | None = None,
    underscore_unmunge: bool = False,
    source: str | pathlib.Path | None = None,
) -> PhyloNode:
    """Initialises a tree.

    Parameters
    ----------
    treestring
        a newick or xml formatted tree string
    tip_names
        a list of tip names, returns a "star" topology tree
    format_name
        indicates treestring is either newick or xml formatted, default
        is newick
    underscore_unmunge
        replace underscores with spaces in all names read, i.e. "sp_name"
        becomes "sp name"
    source
        path to file tree came from, string value assigned to tree.source

    Notes
    -----
    Underscore unmunging is turned off by default, although it is part
    of the Newick format.

    Returns
    -------
    PhyloNode
    """

    source = str(source) if source else None
    if tip_names:
        tree_builder = TreeBuilder().create_edge
        tips = [
            tree_builder([], str(tip_name), {}, None, None) for tip_name in tip_names
        ]
        result = tree_builder(tips, "root", {}, None, None)
        result.source = source
        return result

    if not treestring:
        msg = "Must provide either treestring or tip_names."
        raise ValueError(msg)

    if format_name is None and treestring.startswith("<"):
        format_name = "xml"

    tree_builder = TreeBuilder().create_edge
    # FIXME: More general strategy for underscore_unmunge
    tree = newick_parse_string(
        treestring, tree_builder, underscore_unmunge=underscore_unmunge
    )
    if not tree.name_loaded:
        tree.name = "root"

    tree.source = source
    return tree


def load_tree(
    filename: str | pathlib.Path,
    format_name: str | None = None,
    underscore_unmunge: bool = False,
) -> PhyloNode:
    """Constructor for tree.

    Parameters
    ----------
    filename
        a file path containing a newick or xml formatted tree.
    format_name
        either xml or json, all other values default to newick. Overrides
        file name suffix.
    underscore_unmunge
        replace underscores with spaces in all names read, i.e. "sp_name"
        becomes "sp name".

    Notes
    -----
    Underscore unmunging is turned off by default, although it is part
    of the Newick format. Only the cogent3 json and xml tree formats are
    supported.

    filename is assigned to root node tree.source attribute.

    Returns
    -------
    PhyloNode
    """
    fmt, _ = get_format_suffixes(filename)
    format_name = format_name or fmt
    if format_name == "json":
        tree = load_from_json(filename, (PhyloNode,))
        tree.source = str(filename)
        return tree

    with open_(filename) as tfile:
        treestring = tfile.read()

    return make_tree(
        treestring,
        format_name=format_name,
        underscore_unmunge=underscore_unmunge,
        source=filename,
    )


@register_deserialiser("cogent3.core.tree")
def deserialise_tree(
    data: dict[str, str | dict[str, Any]],
) -> PhyloNode:
    """returns a cogent3 PhyloNode instance"""
    # we load tree using make_tree, then populate edge attributes
    edge_attr = cast("dict[str, dict[str, Any]]", data["edge_attributes"])
    length_and_support: dict[str, dict[str, float | None]] | None = cast(
        "dict[str, dict[str, float | None]] | None", data.get("length_and_support")
    )

    if length_and_support is None:
        warnings.warn(
            "Outdated tree json. Please update by regenerating the json on the loaded tree.",
            stacklevel=3,
        )

    tree = make_tree(treestring=data["newick"])
    for edge in tree.preorder():
        params = edge_attr.get(edge.name, {})
        if length_and_support is None:
            edge.length = params.pop("length", None)
            edge.support = params.pop("support", None)
        else:
            edge.length = length_and_support[edge.name]["length"]
            edge.support = length_and_support[edge.name]["support"]
        edge.params.update(params)
    return tree
