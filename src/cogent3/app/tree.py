from typing import Union

from cogent3 import make_tree
from cogent3.phylo.nj import gnj

from .composable import define_app
from .typing import PairwiseDistanceType, SerialisableType, TreeType


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.10.31a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


@define_app
class scale_branches:
    """Transforms tree branch lengths from nucleotide to codon, or the
    converse. Returns a new tree with lengths divided by scalar"""

    def __init__(self, nuc_to_codon=None, codon_to_nuc=None, scalar=1, min_length=1e-6):
        """
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
        assert not all([nuc_to_codon, codon_to_nuc])
        if nuc_to_codon:
            scalar = 3
        elif codon_to_nuc:
            scalar = 1 / 3
        else:
            assert scalar is not None

        self._scalar = scalar
        self._min_length = min_length

    def main(self, tree: TreeType) -> Union[SerialisableType, TreeType]:
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


@define_app
class uniformize_tree:
    """Standardises the orientation of unrooted trees."""

    def __init__(self, root_at="midpoint", ordered_names=None):
        """
        Parameters
        ----------
        root_at
            edge to root at, defaults to midpoint
        ordered_names
            ordering of names, if not provided, taken from first
            tree
        """
        self._root_at = root_at
        self._ordered_names = ordered_names

    def main(self, tree: TreeType) -> Union[SerialisableType, TreeType]:
        if self._root_at == "midpoint":
            new = tree.root_at_midpoint()
        else:
            new = tree.rooted_with_tip(self._root_at)

        if self._ordered_names is None:
            self._ordered_names = tree.get_tip_names()
        new = new.sorted(self._ordered_names)
        return new


@define_app
class quick_tree:
    """Computes a Neighbour Joining tree from pairwise distances."""

    def __init__(self, drop_invalid=False):
        """
        Parameters
        ----------
        drop_invalid : bool
            drops all rows / columns with an invalid entry. If True, sequences
            for which a distance could not be calculated are excluded and
            the resulting tree will be for the subset of labels with strictly
            valid distances. If False, an ArithmeticError is raised if a
            distance is invalid.
        """
        self._drop_invalid = drop_invalid

    def main(self, dists: PairwiseDistanceType) -> Union[SerialisableType, TreeType]:
        """estimates a neighbor joining tree"""
        size = dists.shape[0]
        dists = dists.drop_invalid() if self._drop_invalid else dists
        if dists is None or (dists.shape[0] != size and not self._drop_invalid):
            raise ValueError("invalid pairwise distances")

        # how many species do we have
        if size == 2:
            dist = dists.array[0, 1] / 2.0
            newick = ",".join(f"{sp}:{dist:.4f}" for sp in dists.names)
            newick = f"({newick})"
            tree = make_tree(treestring=newick, underscore_unmunge=True)
        else:
            (result,) = gnj(dists.to_dict(), keep=1, show_progress=False)
            (score, tree) = result

        return tree
