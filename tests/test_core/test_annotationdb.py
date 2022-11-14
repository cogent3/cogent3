from cogent3 import load_unaligned_seqs, open_
from cogent3.core.annotation_db import GenbankAnnotationDb, GffAnnotationDb
from cogent3.parse.genbank import MinimalGenbankParser
from cogent3.parse.gff import gff_parser
import os

base_path = os.path.dirname(os.path.dirname(__file__))
data_path = os.path.join(base_path, "data")
gff_path = os.path.join(data_path, "prok_NoLocusTags.gff")
genbank_path = os.path.join(data_path, "NC_000913.3.gb")
fasta_path = os.path.join(data_path, "short.fa")


__author__ = "Kirat Alreja, Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Matthew Wakefield", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Prototype"


def test_make_db():
    """Test if the gff3 database is correctly created
    in terms of the names of columns & the number of entries"""

    db = GffAnnotationDb(gff_path)
    sql_columns = []
    data = db.db.execute("""SELECT * FROM GFF""")
    for column in data.description:
        sql_columns.append(column[0])
    got = sql_columns
    parsed_gff = gff_parser(gff_path)
    expected = list((list(parsed_gff)[0]).keys())[:-1]
    assert got == expected


def test_make_sql_query():
    """Test that the SQLlite queries are correctly formed"""

    db = GffAnnotationDb(gff_path)

    """only the bio_type is provided"""
    expected = ("SELECT * FROM GFF WHERE Type == ?", ["gene"])
    got = db._make_sql_query(
        seq_name=None, bio_type="gene", identifier=None, start=None, end=None
    )
    assert got == expected

    """only the identifier is provided"""
    expected = ("SELECT * FROM GFF WHERE Attributes like ?", ["%RandomIdentifier%"])
    got = db._make_sql_query(
        seq_name=None,
        bio_type=None,
        identifier="RandomIdentifier",
        start=None,
        end=None,
    )
    assert got == expected

    """both identifier and bio_type provided"""
    expected = (
        "SELECT * FROM GFF WHERE Type == ? AND Attributes like ?",
        ["CDS", "%RandomIdentifier%"],
    )
    got = db._make_sql_query(
        seq_name=None,
        bio_type="CDS",
        identifier="RandomIdentifier",
        start=None,
        end=None,
    )
    assert got == expected

    """start provided along with identifier and bio_type"""
    expected = (
        "SELECT * FROM GFF WHERE Type == ? AND Attributes like ? AND Start >= ?",
        ["CDS", "%RandomIdentifier%", 0],
    )
    got = db._make_sql_query(
        seq_name=None, bio_type="CDS", identifier="RandomIdentifier", start=0, end=None
    )
    assert got == expected

    """end provided along with identifier and bio_type"""
    expected = (
        "SELECT * FROM GFF WHERE Type == ? AND Attributes like ? AND End < ?",
        ["CDS", "%RandomIdentifier%", 5000],
    )
    got = db._make_sql_query(
        seq_name=None,
        bio_type="CDS",
        identifier="RandomIdentifier",
        start=None,
        end=5000,
    )
    assert got == expected

    """start and end provided along with identifier and bio_type"""
    expected = (
        "SELECT * FROM GFF WHERE Type == ? AND Attributes like ? AND Start >= ? AND End < ?",
        ["CDS", "%RandomIdentifier%", 0, 5000],
    )
    got = db._make_sql_query(
        seq_name=None, bio_type="CDS", identifier="RandomIdentifier", start=0, end=5000
    )
    assert got == expected

    """all five attributes provided"""
    expected = (
        "SELECT * FROM GFF WHERE SeqID == ? AND Type == ? AND Attributes like ? AND Start >= ? AND End < ?",
        ["1", "CDS", "%RandomIdentifier%", 0, 5000],
    )
    got = db._make_sql_query(
        seq_name="1", bio_type="CDS", identifier="RandomIdentifier", start=0, end=5000
    )
    assert got == expected

    """check exception when both bio_type and identifier are missing"""
    import pytest

    with pytest.raises(ValueError):
        db._make_sql_query(
            seq_name=None, bio_type=None, identifier=None, start=None, end=None
        )

    """check exception when both bio_type and identifier are missing"""

    with pytest.raises(ValueError):
        db._make_sql_query(
            seq_name="1", bio_type=None, identifier=None, start=0, end=1000
        )

    """check exception when only seq_name is provided"""

    with pytest.raises(ValueError):
        db._make_sql_query(
            seq_name="1", bio_type=None, identifier=None, start=None, end=None
        )

    """check exception when only start is provided"""

    with pytest.raises(ValueError):
        db._make_sql_query(
            seq_name=None, bio_type=None, identifier=None, start=0, end=None
        )

    """check exception when only end is provided"""

    with pytest.raises(ValueError):
        db._make_sql_query(
            seq_name=None, bio_type=None, identifier=None, start=None, end=1000
        )


def test_populate_from_file():
    """Test that the database is populated with the correct
    number of columns"""

    db = GffAnnotationDb(gff_path)
    db.db.execute(""" SELECT * FROM GFF """)
    got = len(list(db.db.fetchall()))
    parsed_gff = gff_parser(gff_path)
    expected = len(list(parsed_gff))
    assert got == expected


def test_db_query():

    """Test that the SQL query returns the correct
    number of rows for different combinations of bio_type/identifier"""

    seqs = load_unaligned_seqs(fasta_path, moltype="dna")
    seq = seqs.seqs[0]

    db = GffAnnotationDb(gff_path)

    """multiple hits for the same identifier"""
    got = len(db.db_query(start=0, end=len(seq), identifier="CDS4"))
    expected = 2
    assert got == expected

    """query for an ID and recieve the children to the ID along with the parent"""
    got = len(db.db_query(start=0, end=len(seq), identifier="gene4"))
    expected = 3
    assert got == expected

    """query for an ID, with no children"""
    got = len(db.db_query(start=0, end=len(seq), identifier="trna1"))
    expected = 1
    assert got == expected

    """query for a bio type, with multiple hits"""
    got = len(db.db_query(start=0, end=len(seq), bio_type="gene"))
    expected = 8
    assert got == expected

    """query for an ID & a bio type, with a single hit"""
    got = len(db.db_query(start=0, end=len(seq), identifier="gene0", bio_type="CDS"))
    expected = 1
    assert got == expected

    """query for an ID & a bio type, with multiple hits"""
    got = len(db.db_query(start=0, end=len(seq), identifier="CDS4", bio_type="CDS"))
    expected = 2
    assert got == expected


def test_find_records():

    """Test that the coordinates the grouped correctly, and features
    formed properly"""

    db = GffAnnotationDb(gff_path)
    seqs = load_unaligned_seqs(fasta_path, moltype="dna")
    seq = seqs.seqs[0]

    """combine rows with the same ID"""
    got = len(db.find_records(start=0, end=len(seq), identifier="CDS4"))
    expected = 1
    assert got == expected

    """combine rows with the same ID, when bio_type given"""
    got = len(db.find_records(start=0, end=len(seq), bio_type="CDS"))
    expected = 6
    assert got == expected

    """combine rows with the same ID, when children rows are fetched with the parent"""
    got = len(db.find_records(start=0, end=len(seq), identifier="gene4"))
    expected = 2
    assert got == expected

    """unique ID, single row returned"""
    got = len(db.find_records(start=0, end=len(seq), identifier="id020000"))
    expected = 1
    assert got == expected


def test_describe():

    """Test that the number of unique values returned are correct"""

    db = GffAnnotationDb(gff_path)

    """An empty set when nothing is provided"""
    got = db.describe()
    expected = dict()
    assert got == expected

    """Only the bio_type is provided"""
    got = len(db.describe(bio_type=True)["Type"])
    expected = 9
    assert got == expected

    """Only the seq_name is provided"""
    got = len(db.describe(seq_name=True)["SeqID"])
    expected = 1
    assert got == expected

    """Only the identifier is provided"""
    got = len(db.describe(identifier=True)["identifier"])
    expected = 21
    assert got == expected

    """Check that the set has three keys when all three arguments are True"""
    got = len(db.describe(bio_type=True, seq_name=True, identifier=True))
    expected = 3
    assert got == expected

    """Correct number of Type values when all three arguments are True"""
    got = len(db.describe(bio_type=True, seq_name=True, identifier=True)["Type"])
    expected = 9
    assert got == expected

    """Correct number of identifier values when all three arguments are True"""
    got = len(db.describe(bio_type=True, seq_name=True, identifier=True)["identifier"])
    expected = 21
    assert got == expected

    """Correct number of SeqID values when all three arguments are True"""
    got = len(db.describe(bio_type=True, seq_name=True, identifier=True)["SeqID"])
    expected = 1
    assert got == expected


def test_ordered_values():
    """Check if the dictonary of GFF records is correctly formed into a list for loading into the in-memory database"""
    from cogent3.core.annotation_db import _ordered_values

    """Check the values explicitly"""
    parsed_gff = gff_parser(gff_path)
    example_dict = list(parsed_gff)[0]
    got = _ordered_values(example_dict)
    expected = [
        "sequence001",
        "mine",
        "gene",
        189,
        255,
        ".",
        "+",
        ".",
        {
            "ID": "gene0",
            "Dbxref": "ASAP:ABE-0000006",
            "gene": "thrL",
            "gene_synonym": "ECK0001",
        },
    ]
    assert got == expected

    """Check the data types of the resulting list"""

    def return_type(values):
        types = []
        for v in values:
            types.append(type(v))
        return types

    got = return_type(_ordered_values(example_dict))
    expected = [str, str, str, int, int, str, str, str, dict]
    assert got == expected


def test_populate_from_file_genbank():
    """Test that the GenBank database is populated with the correct
    number of columns"""

    db = GenbankAnnotationDb(genbank_path)

    """test the number of rows populated, after skipping the
    records without 'locus_tag' key"""
    db.db.execute(""" SELECT * FROM Genbank """)
    got = len(list(db.db.fetchall()))

    with open_(genbank_path) as infile:
        data = list(MinimalGenbankParser(infile.readlines()))

    record = data[0]
    features = record["features"]
    expected = 0
    for feature in features:
        if "locus_tag" not in list(feature.keys()):
            continue
        expected += 1

    assert got == expected


def test_db_query_genbank():
    """Test that the SQL query returns the correct
    number of rows for different combinations of bio_type/identifier"""

    import pytest

    db = GenbankAnnotationDb(genbank_path)
    with open_(genbank_path) as infile:
        data = list(MinimalGenbankParser(infile.readlines()))
    record = data[0]

    """only the bio type is provided, along with start and end"""
    got = len(db.db_query(start=0, end=len(record["sequence"]), bio_type="CDS"))
    expected = 4315
    assert got == expected

    """only the bio type is provided, and the end is sliced to half"""
    got = len(db.db_query(end=len(record["sequence"]) / 2, bio_type="CDS"))
    expected = 2199
    assert got == expected

    """only the identifier is provided, along with start and end"""
    got = len(db.db_query(start=0, end=len(record["sequence"]), identifier="b4515"))
    expected = 2
    assert got == expected

    """both the bio type and identifier are provided"""
    got = len(db.db_query(bio_type="CDS", identifier="b4515"))
    expected = 1
    assert got == expected

    """raise value error when bio type or identifer not provided"""

    with pytest.raises(ValueError):
        db.db_query(start=0, end=len(record["sequence"]))
    with pytest.raises(ValueError):
        db.db_query()
    with pytest.raises(ValueError):
        db.db_query(start=0, end=len(record["sequence"]), seq_name="SeqName")
    with pytest.raises(ValueError):
        db.db_query(seq_name="SeqName")


def make_genbank_db():
    """Check if the genbank database created has the correct columns"""
    sql_columns = []
    db = GenbankAnnotationDb(genbank_path)
    data = db.db.execute("""SELECT * FROM GENBANK""")
    for column in data.description:
        sql_columns.append(column[0])
    got = sql_columns
    expected = ["LocusID", "Type", "Spans", "Locus_Tag", "Start", "End", "Strand"]
    assert got == expected


def test_fetch_from_feature():

    """Check if the genbank database created has the correct columns"""
    from cogent3.core.annotation_db import _fetch_from_features

    with open_(genbank_path) as infile:
        data = list(MinimalGenbankParser(infile.readlines()))

    record = data[0]
    features = record["features"]
    got = _fetch_from_features(features[1])
    expected = ["gene", "[[189, 255]]", "b0001", 189, 255, 1]
    assert got == expected

    """Check the datatype of the values returned"""

    def return_type(values):
        types = []
        for v in values:
            types.append(type(v))
        return types

    got = return_type(_fetch_from_features(features[1]))
    expected = [str, str, str, int, int, int]
    assert got == expected


def test_describe_genbank():
    db = GenbankAnnotationDb(genbank_path)

    got = db.describe()
    expected = dict()
    assert got == expected

    """Only when the bio type is provided"""
    got = len(db.describe(bio_type=True)["Type"])
    expected = 6
    assert got == expected

    """Only when the LocusID/seq_name is provided"""
    got = len(db.describe(seq_name=True)["LocusID"])
    expected = 1
    assert got == expected

    """Only when the identifier is provided"""
    got = len(db.describe(identifier=True)["identifier"])
    expected = 4639
    assert got == expected

    """Check the number of keys created"""
    got = len(db.describe(bio_type=True, identifier=True, seq_name=True))
    expected = 3
    assert got == expected

    """Number of unique bio_type values when all three attributes are passed"""
    got = len(db.describe(bio_type=True, identifier=True, seq_name=True)["Type"])
    expected = 6
    assert got == expected

    """Number of unique LocusID values when all three attributes are passed"""
    got = len(db.describe(bio_type=True, identifier=True, seq_name=True)["LocusID"])
    expected = 1
    assert got == expected

    """Number of unique identifer values when all three attributes are passed"""
    got = len(db.describe(bio_type=True, identifier=True, seq_name=True)["identifier"])
    expected = 4639
    assert got == expected


def test_find_records_genbank():
    import pytest

    db = GenbankAnnotationDb(genbank_path)

    """Number of records found when bio type is passed"""
    got = len(db.find_records(seq_name="NC_000913", bio_type="CDS"))
    expected = 4305
    assert got == expected

    "Number of records found when the identifier is passed"
    got = len(db.find_records(seq_name="NC_000913", identifier="b1492"))
    expected = 2
    assert got == expected

    """Number of records found when both the identifier and bio type are passed"""
    got = len(db.find_records(seq_name="NC_000913", identifier="b1492", bio_type="CDS"))
    expected = 1
    assert got == expected

    """Check value error when bio_type and identifier are not passed"""
    with pytest.raises(ValueError):
        db.find_records(seq_name="NC_000913")
    with pytest.raises(ValueError):
        db.find_records(seq_name="NC_000913", start=0, end=5000)
    with pytest.raises(ValueError):
        db.find_records(start=0, end=5000)


def test_make_sql_query_genbank():
    """Test that the SQLlite queries are correctly formed for GenBank files"""

    db = GenbankAnnotationDb(genbank_path)

    """only the bio_type is provided"""
    expected = ("SELECT * FROM GENBANK WHERE Type == ?", ["gene"])
    got = db._make_sql_query(
        seq_name=None, bio_type="gene", identifier=None, start=None, end=None
    )
    assert got == expected

    """only the identifier is provided"""
    expected = ("SELECT * FROM GENBANK WHERE Locus_Tag like ?", ["%RandomIdentifier%"])
    got = db._make_sql_query(
        seq_name=None,
        bio_type=None,
        identifier="RandomIdentifier",
        start=None,
        end=None,
    )
    assert got == expected

    """both identifier and bio_type provided"""
    expected = (
        "SELECT * FROM GENBANK WHERE Type == ? AND Locus_Tag like ?",
        ["CDS", "%RandomIdentifier%"],
    )
    got = db._make_sql_query(
        seq_name=None,
        bio_type="CDS",
        identifier="RandomIdentifier",
        start=None,
        end=None,
    )
    assert got == expected

    """start provided along with identifier and bio_type"""
    expected = (
        "SELECT * FROM GENBANK WHERE Type == ? AND Locus_Tag like ? AND Start >= ?",
        ["CDS", "%RandomIdentifier%", 0],
    )
    got = db._make_sql_query(
        seq_name=None, bio_type="CDS", identifier="RandomIdentifier", start=0, end=None
    )
    assert got == expected

    """end provided along with identifier and bio_type"""
    expected = (
        "SELECT * FROM GENBANK WHERE Type == ? AND Locus_Tag like ? AND End < ?",
        ["CDS", "%RandomIdentifier%", 5000],
    )
    got = db._make_sql_query(
        seq_name=None,
        bio_type="CDS",
        identifier="RandomIdentifier",
        start=None,
        end=5000,
    )
    assert got == expected

    """start and end provided along with identifier and bio_type"""
    expected = (
        "SELECT * FROM GENBANK WHERE Type == ? AND Locus_Tag like ? AND Start >= ? AND End < ?",
        ["CDS", "%RandomIdentifier%", 0, 5000],
    )
    got = db._make_sql_query(
        seq_name=None, bio_type="CDS", identifier="RandomIdentifier", start=0, end=5000
    )
    assert got == expected

    """all five attributes provided"""
    expected = (
        "SELECT * FROM GENBANK WHERE LocusID == ? AND Type == ? AND Locus_Tag like ? AND Start >= ? AND End < ?",
        ["1", "CDS", "%RandomIdentifier%", 0, 5000],
    )
    got = db._make_sql_query(
        seq_name="1", bio_type="CDS", identifier="RandomIdentifier", start=0, end=5000
    )
    assert got == expected

    """check exception when both bio_type and identifier are missing"""
    import pytest

    with pytest.raises(ValueError):
        db._make_sql_query(
            seq_name=None, bio_type=None, identifier=None, start=None, end=None
        )

    """check exception when both bio_type and identifier are missing, even if other attributes"""

    with pytest.raises(ValueError):
        db._make_sql_query(
            seq_name="1", bio_type=None, identifier=None, start=0, end=1000
        )

    """check exception when only seq_name is provided"""

    with pytest.raises(ValueError):
        db._make_sql_query(
            seq_name="1", bio_type=None, identifier=None, start=None, end=None
        )

    """check exception when only start is provided"""

    with pytest.raises(ValueError):
        db._make_sql_query(
            seq_name=None, bio_type=None, identifier=None, start=0, end=None
        )

    """check exception when only end is provided"""

    with pytest.raises(ValueError):
        db._make_sql_query(
            seq_name=None, bio_type=None, identifier=None, start=None, end=1000
        )


if __name__ == "__main__":
    test_db_query()
    test_find_records()
    test_make_sql_query()
    test_populate_from_file()
    test_describe()
    test_populate_from_file_genbank()
    test_db_query_genbank()
    test_fetch_from_feature()
    test_describe_genbank()
    test_find_records_genbank()
