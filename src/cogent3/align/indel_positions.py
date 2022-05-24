__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def pog_traceback(pogs, aligned_positions):
    upto = [0, 0]
    align_builder = POGBuilder(pogs)
    for posn in aligned_positions:
        assert len(posn) == 2
        for (dim, pos) in enumerate(posn):
            if pos is not None:
                align_builder.add_skipped(dim, upto[dim], pos)
                upto[dim] = pos + 1

        align_builder.add_aligned(posn)

    for dim in [0, 1]:
        align_builder.add_skipped(dim, upto[dim], len(pogs[dim]))

    return align_builder.get_pog()


class POGBuilder(object):
    def __init__(self, children):
        self.children = children
        self.remap = [{} for _ in children]
        self.started = [False, False]
        self.last = [None, None]
        self.result = [[]]
        self.origins = [[]]
        self.aligned_positions = []
        self.states = []

    def add_skipped(self, dim, start, end, old_gap=True):
        for p in range(start, end):
            fp = [None, None]
            fp[dim] = p
            fp = tuple(fp)
            self.add_aligned(fp, old_gap=old_gap)

    def add_aligned(self, posn, old_gap=False):
        pre_merged = set()
        assert len(posn) == 2
        for (dim, pos) in enumerate(posn):
            if pos is None:
                continue
            self.remap[dim][pos] = len(self.aligned_positions)
            self.last[dim] = pos
        self.result.append(pre_merged)
        self.aligned_positions.append(posn)
        if None not in posn:
            state = "m"
        elif posn[0] is None:
            state = "x"
        else:
            state = "y"
        if not old_gap:
            state = state.upper()
        self.states.append(state)

    def get_pog(self):
        jumps = []
        gapmap = {}
        ingap = False
        # Build a list of gaps (ie: segments of X or Y state) in
        # the alignment and a dict which maps from seq posn to the
        # start of the surrounding gap.
        for (i, state) in enumerate(self.states + ["."]):
            gap = state in "XYxy"
            if gap and not ingap:
                start = i
                ingap = True
            elif ingap and not gap:
                jumps.append((start, i))
                ingap = False
            if ingap:
                gapmap[i] = start

        # in case of tail gap
        for (dim, child) in enumerate(self.children):
            pos = len(child)
            self.remap[dim][pos] = len(self.aligned_positions)

        # Keep only those child gaps which sit entirely within a gap
        # in this alignment
        child_jumps = []
        for (dim, pog) in enumerate(self.children):
            r = self.remap[dim]
            for (i, j) in pog.jumps:
                (i, j) = (r[i], r[j])
                if i in gapmap and j in gapmap and gapmap[i] == gapmap[j]:
                    child_jumps.append((i, j))

        pog = POG(len(self.aligned_positions), jumps, child_jumps)
        pog.aligned_positions = self.aligned_positions
        pog.states = "".join(self.states)
        return pog


class POG(object):
    """A representation of the indel positions in a pairwise alignment, ie:
    those segments of the consensus sequence which may be inserts and so absent
    from the common ancestor.  Nearly equivalent to a generic Partial Order
    Graph.

    Indels are represented as tuples of
        (1st posn in indel, 1st posn after indel)

    Two lists of indels are kept, one for indels in the alignment, and one
    for indels in its two children in case they are also alignments.

    This data structure largely inspired by:
    Loytynoja A, Goldman N. 2005. An algorithm for progressive multiple
    alignment of sequences with insertions. PNAS 102:10557-10562
    """

    def __init__(self, length, jumps, child_jumps):
        self.jumps = jumps
        self.child_jumps = child_jumps
        self.all_jumps = self.jumps + self.child_jumps
        self.all_jumps.sort(key=lambda i_j: i_j[1])
        self.length = length
        for (i, j) in self.all_jumps:
            assert i <= j, (length, jumps, child_jumps)
            assert 0 <= i <= length, (length, jumps, child_jumps)
            assert 0 <= j <= length, (length, jumps, child_jumps)

    def traceback(self, other, aligned_positions):
        return pog_traceback([self, other], aligned_positions)

    def as_list_of_pred_lists(self):
        """A representation of the POG as a list of predecessor positions,
        a simple way to represent DAGs eg: [], [0], [1] would be a simple
        sequence of length 3.  Extra start and end positions are added, so
        the length is len(self)+2 and the positions are all offset by 1"""
        result = [[]]
        # First the regular, linear sequence relationships
        for i in range(self.length + 1):
            pre = [i]
            result.append(pre)
        # Then add in the indel jumps.  Given an indel from i to j
        # j could have been ajacent to one of i's predecessors in
        # the ancestral sequence. This depends on all_jumps being sorted
        # by j.
        for (i, j) in self.all_jumps:
            if i == j:
                continue
            assert i < j
            result[j + 1].extend(result[i + 1])
        return result

    def get_aligned_positions(self):
        return self.aligned_positions

    def get_full_aligned_positions(self):
        return self.aligned_positions

    def __len__(self):
        return self.length

    def midlinks(self):
        # for the hirchberg algorithm.
        half = self.length // 2
        jumps = [(i, j) for (i, j) in self.all_jumps if i <= half and j >= half]
        return [(half, half)] + jumps

    def __getitem__(self, index):
        # POGs need to be sliceable for the hirchberg algorithm.
        if index.start is None:
            start = 0
        else:
            start = index.start
        if index.stop is None:
            end = self.length
        else:
            end = index.stop
        assert end >= start, (start, end, index, self.length)

        def moved(i, j):
            i2 = max(min(i, end), start) - start
            j2 = max(min(j, end), start) - start
            return (i2, j2)

        jumps = [moved(i, j) for (i, j) in self.jumps if i < end or j > start]
        cjumps = [moved(i, j) for (i, j) in self.child_jumps if i < end or j > start]
        return POG(end - start, jumps, cjumps)

    def backward(self):
        # Switches predecessors / successors
        # POGs need to be reversable for the hirchberg algorithm.
        length = self.length
        jumps = [(length - j, length - i) for (i, j) in self.jumps]
        cjumps = [(length - j, length - i) for (i, j) in self.child_jumps]
        return POG(length, jumps, cjumps)

    def write_to_dot(self, dot):
        pred_sets = self.as_list_of_pred_lists()
        print("digraph POG {", file=dot)
        for (i, preds) in enumerate(pred_sets):
            # print i, preds
            for pred in preds:
                print("  ", (f"node{pred} -> node{i}"), file=dot)
            if i == 0:
                label = "START"
            elif i == len(pred_sets) - 1:
                label = "END"
            else:
                label = str(i)
            print("  ", (f"node{i}"), f'[label="{label}"]', file=dot)
        print("}", file=dot)
        print("", file=dot)


class LeafPOG(POG):
    """The POG for a known sequence contains no indels."""

    def __init__(self, length):
        self.length = length
        self.all_jumps = []
        self.jumps = []

    def as_list_of_pred_lists(self):
        pog = [[[i]] for i in range(self.length)]
        return [[]] + pog + [[len(pog)]]

    def __len__(self):
        return self.length

    def backward(self):
        return LeafPOG(self.length)


def leaf2pog(leaf):
    return LeafPOG(len(leaf))
