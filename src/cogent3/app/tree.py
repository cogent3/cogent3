from cogent3 import make_tree
from cogent3.phylo.nj import gnj

from .composable import ComposableTree

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.8.30a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class scale_branches(ComposableTree):
    """Transforms tree branch lengths from nucleotide to codon, or the converse.
    Returns a Tree."""

    _input_type = frozenset(["tree"])
    _output_type = frozenset(["tree", "serialisable"])
    _data_types = frozenset(["PhyloNode", "TreeNode"])

    def __init__(self, nuc_to_codon=None, codon_to_nuc=None, scalar=1, min_length=1e-6):
        super(scale_branches, self).__init__()
        """returns a new tree with lengths divided by scalar
    
        Parameters
        ----------
        nuc_to_codon : bool
            scalar is 3
        codon_to_nuc : bool
            scalar is 1/3
        scalar
            numerical value. Overridden by either nuc_to_codon or
            codon_to_nuc
        min_length
            set branch length to this value if it's not defined,
            or <= zero
        """
        self._formatted_params()
        assert not all([nuc_to_codon, codon_to_nuc])
        if nuc_to_codon:
            scalar = 3
        elif codon_to_nuc:
            scalar = 1 / 3
        else:
            assert scalar is not None

        self._scalar = scalar
        self._min_length = min_length
        self.func = self._scale_lengths

    def _scale_lengths(self, tree):
        scalar = self._scalar
        min_length = self._min_length
        tree = tree.deepcopy()
        for edge in tree.postorder():
            if edge.name == "root":
                continue  # root has no length
            length = edge.length
            length = length / scalar if length else min_length
            edge.length = abs(length)  # force to be positive

        return tree


class uniformize_tree(ComposableTree):
    """Standardises the orientation of unrooted trees. Returns a Tree."""

    _input_type = frozenset(["tree"])
    _output_type = frozenset(["tree", "serialisable"])
    _data_types = frozenset(["PhyloNode", "TreeNode"])

    def __init__(self, root_at="midpoint", ordered_names=None):
        super(uniformize_tree, self).__init__()
        """returns a new tree with standardised orientation
        
        Parameters
        ----------
        root_at
            edge to root at, defaults to midpoint
        ordered_names
            ordering of names, if not provided, taken from first
            tree
        """
        self._formatted_params()
        self._root_at = root_at
        self._ordered_names = ordered_names
        self.func = self._uniformize

    def _uniformize(self, tree):
        if self._root_at == "midpoint":
            new = tree.root_at_midpoint()
        else:
            new = tree.rooted_with_tip(self._root_at)

        if self._ordered_names is None:
            self._ordered_names = tree.get_tip_names()
        new = new.sorted(self._ordered_names)
        return new


class quick_tree(ComposableTree):
    """Neighbour Joining tree based on pairwise distances."""

    _input_type = frozenset(["pairwise_distances"])
    _output_type = frozenset(["tree", "serialisable"])
    _data_types = frozenset(["ArrayAlignment", "Alignment", "DistanceMatrix"])

    def __init__(self, drop_invalid=False):
        """computes a neighbour joining tree from an alignment"""
        super(quick_tree, self).__init__()
        self._formatted_params()
        self.func = self.quick_tree
        # removes records with missing distances symmetrically.
        self._drop_invalid = drop_invalid

    def quick_tree(self, distance_matrix):
        """estimates a neighbor joining tree"""
        dists = distance_matrix
        size = dists.shape[0]
        dists = dists.drop_invalid() if self._drop_invalid else dists
        if dists is None or dists.shape[0] != size:
            msg = (
                f"some pairwise distances could not be computed with"
                " {self._distance}, pick a different distance"
            )
            raise ValueError(msg)

        # how many species do we have
        species = dists.keys()
        if len(species) == 2:
            dist = list(dists.values())[0] / 2.0
            treestring = "(%s:%.4f,%s:%.4f)" % (species[0], dist, species[1], dist)
            tree = make_tree(treestring=treestring, underscore_unmunge=True)
        else:
            (result,) = gnj(dists.todict(), keep=1, show_progress=False)
            (score, tree) = result

        return tree
