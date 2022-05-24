#!/usr/bin/env python
"""Extracts data from NCBI nodes.dmp and names.dmp files.
"""
from functools import total_ordering

from cogent3.core.tree import TreeNode


__author__ = "Jason Carnes"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Jason Carnes", "Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Jason Carnes"
__email__ = "jason.carnes@sbri.org"
__status__ = "Development"

strip = str.strip


class MissingParentError(Exception):
    pass


# Note: numbers not guaranteed to be consistent if new taxa are invented...
RanksToNumbers = {
    "forma": 1,
    "varietas": 2,
    "subspecies": 3,
    "species": 4,
    "species subgroup": 5,
    "species group": 6,
    "subgenus": 7,
    "genus": 8,
    "subtribe": 9,
    "tribe": 10,
    "subfamily": 11,
    "family": 12,
    "superfamily": 13,
    "parvorder": 14,
    "infraorder": 15,
    "suborder": 16,
    "order": 17,
    "superorder": 18,
    "infraclass": 19,
    "subclass": 20,
    "class": 21,
    "superclass": 22,
    "subphylum": 23,
    "phylum": 24,
    "superphylum": 25,
    "kingdom": 26,
    "superkingdom": 27,
    "no rank": 28,
}


@total_ordering
class NcbiTaxon(object):
    """Extracts taxon information: init from one line of NCBI's nodes.dmp.

    Properties:
        TaxonId     ID of this node
        ParentId    ID of this node's parent
        Rank        Rank of this node: genus, species, etc.
        EmblCode    Locus name prefix; not unique
        DivisionId  From division.dmp
        DivisionInherited  1 or 0; 1 if node inherits division from parent
        TranslTable  ID of this node's genetic code from gencode.dmp
        GCInherit   1 or 0; 1 if node inherits genetic code from parent
        TranslTableMt    ID of this node's mitochondrial code from gencode.dmp
        TranslTableMtInherited 1 or 0; 1 if node inherits mt code from parent
        Hidden      1 or 0; 1 if hidden by default in GenBank's listing
        HiddenSubtreeRoot   1 or 0; 1 if no sequences from this subtree exist
        Comments    free-text comments

        RankId      Arbitrary number corresponding to rank. See RanksToNumbers.
        Name        Name of this node: must get from external source. Thanks
                    so much, NCBI...
                    Expect a string: '' by default.
    """

    Fields = [
        "TaxonId",
        "ParentId",
        "Rank",
        "EmblCode",
        "DivisionId",
        "DivisionInherited",
        "TranslTable",
        "TranslTableInherited",
        "TranslTableMt",
        "TranslTableMtInherited",
        "Hidden",
        "HiddenSubtreeRoot",
        "Comments",
    ]

    def __init__(self, line):
        """Returns new NcbiTaxon from line containing taxonomy data."""
        line_pieces = list(map(strip, line.split("|")))
        for i in [0, 1, 5, 6, 7, 8, 9, 10, 11]:
            line_pieces[i] = int(line_pieces[i])
        # fix trailing delimiter
        last = line_pieces[-1]
        if last.endswith("|"):
            line_pieces[-1] = last[:-1]
        self.__dict__ = dict(list(zip(self.Fields, line_pieces)))
        self.Name = ""  # will get name field from names.dmp; fillNames
        self.RankId = RanksToNumbers.get(self.Rank, None)

    def __str__(self):
        """Writes data out in format we got it."""
        pieces = [str(getattr(self, f)) for f in self.Fields]
        # remember to set the parent of the root to itself
        if pieces[1] == "None":
            pieces[1] = pieces[0]
        return "\t|\t".join(pieces) + "\t|\n"

    def __lt__(self, other):
        """Compare by taxon rank."""
        try:
            return self.RankId < other.RankId
        except AttributeError:
            return True  # always sort ranked nodes above unranked

    def __eq__(self, other):
        """Compare by taxon rank."""
        try:
            return self.RankId == other.RankId
        except AttributeError:
            return True

    def __ne__(self, other):
        """Compare by taxon rank."""
        try:
            return self.RankId != other.RankId
        except AttributeError:
            return True


def NcbiTaxonParser(infile):
    """Returns a sequence of NcbiTaxon objects from sequence of lines."""
    for line in infile:
        if line.strip():
            yield NcbiTaxon(line)


def NcbiTaxonLookup(taxa):
    """Returns dict of TaxonId -> NcbiTaxon object."""
    result = {}
    for t in taxa:
        result[t.TaxonId] = t
    return result


class NcbiName(object):
    """Extracts name information: init from one line of NCBI's names.dmp.

    Properties:
        TaxonId     TaxonId of this node
        Name        Text representation of the name, e.g. Homo sapiens
        UniqueName  The unique variant of this name if Name not unique
        NameClass   Kind of name, e.g. scientific name, synonym, etc.
    """

    Fields = ["TaxonId", "Name", "UniqueName", "NameClass"]

    def __init__(self, line):
        """Returns new NcbiName from line containing name data."""
        line_pieces = list(map(strip, line.split("|")))
        line_pieces[0] = int(line_pieces[0])  # convert taxon_id
        self.__dict__ = dict(list(zip(self.Fields, line_pieces)))

    def __str__(self):
        """Writes data out in similar format as the one we got it from."""
        return "\t|\t".join([str(getattr(self, f)) for f in self.Fields]) + "|\n"


def NcbiNameParser(infile):
    """Returns sequence of NcbiName objects from sequence of lines."""
    for line in infile:
        if line.strip():
            yield NcbiName(line)


def NcbiNameLookup(names):
    """Returns dict mapping taxon id -> NCBI scientific name."""
    result = {}
    for name in names:
        if name.NameClass == "scientific name":
            result[name.TaxonId] = name
    return result


class NcbiTaxonomy(object):
    """Holds root node of a taxonomy tree, plus lookup by id or name."""

    def __init__(self, taxa, names, strict=False):
        """Creates new taxonomy, using data in Taxa and Names.

        taxa should be the product of NcbiTaxonLookup.

        names should be the product of NcbiNameLookup.

        strict, if True, raises an error on finding taxa whose parents don't
        exist. Otherwise, will put them in self.Deadbeats keyed by parent ID.

        Note: because taxa is a dict, nodes will be added in arbitrary order.
        """
        names_to_nodes = {}
        ids_to_nodes = {}
        for t_id, t in taxa.items():
            name_rec = names.get(t_id, None)
            if name_rec:
                name = name_rec.Name
            else:
                name = "Unknown"
            t.Name = name

            node = NcbiTaxonNode(t)
            names_to_nodes[name] = node
            ids_to_nodes[t_id] = node
        self.ByName = names_to_nodes
        self.ById = ids_to_nodes

        deadbeats = {}
        # build the tree by connecting each node to its parent
        for t_id, t in ids_to_nodes.items():
            if t.ParentId == t.TaxonId:
                t.parent = None
            else:
                try:
                    ids_to_nodes[t.ParentId].append(t)
                except KeyError:  # found a child whose parent doesn't exist
                    if strict:
                        raise MissingParentError(
                            f"Node {t_id} has parent {t.ParentId}, which isn't in taxa."
                        )
                    else:
                        deadbeats[t.ParentId] = t
        self.Deadbeats = deadbeats
        self.Root = t.root()

    def __getitem__(self, item):
        """If item is int, returns taxon by id: otherwise, searches by name.

        Returns the relevant NcbiTaxonNode.
        Will raise KeyError if not present.
        """
        try:
            return self.ById[int(item)]
        except ValueError:
            return self.ByName[item]


class NcbiTaxonNode(TreeNode):
    """Provides some additional methods specific to Ncbi taxa."""

    def __init__(self, Data):
        """Returns a new NcbiTaxonNode object; requires NcbiTaxon to initialize."""
        self.Data = Data
        self._parent = None
        self.children = []

    def getRankedDescendants(self, rank):
        """Returns all descendants of self with specified rank as flat list."""
        curr = self.Rank
        if curr == rank:
            result = [self]
        else:
            result = []
        for i in self:
            result.extend(i.getRankedDescendants(rank))
        return result

    def _get_parent_id(self):
        return self.Data.ParentId

    ParentId = property(_get_parent_id)

    def _get_taxon_id(self):
        return self.Data.TaxonId

    TaxonId = property(_get_taxon_id)

    def _get_rank(self):
        return self.Data.Rank

    Rank = property(_get_rank)

    def _get_name(self):
        return self.Data.Name

    Name = property(_get_name)


def NcbiTaxonomyFromFiles(nodes_file, names_file, strict=False):
    """Returns new NcbiTaxonomy fron nodes and names files."""
    taxa = NcbiTaxonLookup(NcbiTaxonParser(nodes_file))
    names = NcbiNameLookup(NcbiNameParser(names_file))
    return NcbiTaxonomy(taxa, names, strict)
