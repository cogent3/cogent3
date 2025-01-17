#!/usr/bin/env python
"""Provides Info, DbRef, DbRefs

Info is a dictionary and is the annotation object of a Sequence object.
"""

from warnings import warn

from cogent3.parse.record import MappedRecord
from cogent3.util.misc import ConstrainedDict, Delegator, FunctionWrapper


class DbRef:
    """Holds a database accession, and optionally other data.

    Accession:      id in the database: str or int
    Db:             database name: str
    name:           short name of the record: str
    Description:    description of the record, possibly lengthy: str
    Data:           any data associated with the record: arbitrary object

    str(DbRef) always returns the accession.
    """

    def __init__(self, Accession, Db="", name="", Description="", Data=None) -> None:
        """Returns new DbRef.

        str(DbRef) always returns the accession as a string.
        """
        self.Accession = Accession
        self.Db = Db
        self.name = name
        self.Description = Description
        self.Data = Data

    def __str__(self) -> str:
        """Returns accession."""
        return str(self.Accession)

    def __int__(self) -> int:
        """Tries to coerce accession to int."""
        return int(self.Accession)

    def __gt__(self, other):
        """Compares by accession: tries numeric first, then alphabetic"""
        try:
            return int(self) > int(other)
        except:
            return str(self) > str(other)

    def __lt__(self, other):
        """Compares by accession: tries numeric first, then alphabetic"""
        try:
            return int(self) < int(other)
        except:
            return str(self) < str(other)

    def __eq__(self, other):
        """Compares by accession: tries numeric first, then alphabetic"""
        try:
            return int(self) == int(other)
        except:
            return str(self) == str(other)

    def __ne__(self, other):
        """Compares by accession: tries numeric first, then alphabetic"""
        try:
            return int(self) != int(other)
        except:
            return str(self) != str(other)


def _make_list(obj):
    """Returns list corresponding to or containing obj, depending on type."""
    if isinstance(obj, list):
        return obj
    if isinstance(obj, tuple):
        return list(obj)
    return [obj]


class DbRefs(MappedRecord, ConstrainedDict):
    """Holds Database -> [Accessions] mapping.

    The accessions for a particular database are always stored as a list.

    DbRefs will ultimately contain methods for actually getting the records
    from known databases.
    """

    value_mask = FunctionWrapper(_make_list)
    DefaultValue = []


KnownDatabases = dict.fromkeys(
    [
        "RefSeq",
        "GenBank",
        "GenNucl",
        "GenPept",
        "GI",
        "SwissProt",
        "PIR",
        "EMBL",
        "DDBJ",
        "NDB",
        "PDB",
        "Taxon",
        "LocusLink",
        "UniGene",
        "OMIM",
        "PubMed",
        "COGS",
        "CDD",
        "Pfam",
        "Rfam",
        "GO",
        "dbEST",
        "IPI",
        "rRNA",
        "EC",
        "HomoloGene",
        "KEGG",
        "BRENDA",
        "EcoCyc",
        "HumanCyc",
        "BLOCKS",
    ],
)


class Info(MappedRecord, Delegator):
    """Dictionary that stores attributes for Sequence objects.

    Delegates to DbRefs for database IDs.
    """

    Required = {"Refs": None}

    def __init__(self, *args, **kwargs) -> None:
        """Returns new Info object. Creates DbRefs if necessary."""
        temp = dict(*args, **kwargs)
        if "Refs" in temp:
            refs = temp["Refs"]
            if not isinstance(refs, DbRefs):
                refs = DbRefs(refs)
        else:
            refs = DbRefs()
        # move keys into refs if they belong there: allows init from flat dict
        for key, val in list(temp.items()):
            if key in KnownDatabases:
                refs[key] = val
                del temp[key]

        Delegator.__init__(self, refs)
        self["Refs"] = refs
        MappedRecord.__init__(self, temp)

    def __getattr__(self, attr):
        """Checks for attr in Refs first."""
        if attr in KnownDatabases:
            return getattr(self.Refs, attr)
        return super().__getattr__(attr)

    def __setattr__(self, attr, val) -> None:
        """Try to set in Refs first."""
        if attr in KnownDatabases:
            return setattr(self.Refs, attr, val)
        return super().__setattr__(attr, val)

    def __getitem__(self, item):
        """Checks for item in Refs first."""
        if item in KnownDatabases:
            return getattr(self.Refs, item)
        return super().__getitem__(item)

    def __setitem__(self, item, val) -> None:
        """Try to set in Refs first."""
        if item in KnownDatabases:
            return setattr(self.Refs, item, val)
        return super().__setitem__(item, val)

    def __contains__(self, item) -> bool:
        """Checks for item in Refs first."""
        if item in KnownDatabases:
            return item in self.Refs
        return super().__contains__(item)

    def update(self, item):
        """updates with another info object and warns when overwriting keys"""
        if overwrites := (set(self) ^ {"Refs"}) & ((set(item)) ^ {"Refs"}):
            warn(
                "Keys overwritten by other sequence: " + "".join(overwrites),
                stacklevel=2,
            )
        return super().update(item)
