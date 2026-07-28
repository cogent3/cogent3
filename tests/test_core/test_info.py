"""Unit tests for Info class and associated objects (DbRef, DbRefs, etc.)."""

import warnings

import pytest

from cogent3.core.info import DbRef, DbRefs, Info, _make_list


@pytest.fixture
def data():
    return {
        "Accession": "xyz",
        "Db": "abc",
        "name": "qwe",
        "Description": "blah",
        "Data": list(range(20)),
    }


@pytest.fixture
def db(data):
    return DbRef(**data)


def test_init_minimal():
    # DbRef minimal init should fill fields as expected
    d = DbRef("abc")
    assert d.Accession == "abc"
    assert d.Db == ""
    assert d.name == ""
    assert d.Description == ""
    assert d.Data is None
    # empty init not allowed
    with pytest.raises(TypeError):
        DbRef()


@pytest.mark.parametrize(
    ("attr", "val"),
    [
        ("Accession", "xyz"),
        ("Db", "abc"),
        ("name", "qwe"),
        ("Description", "blah"),
        ("Data", list(range(20))),
    ],
)
def test_init(db, attr, val):
    # DbRef init should insert correct data
    assert getattr(db, attr) == val


def test_str(db):
    # DbRef str should be the same as the accession str
    assert str(db) == "xyz"
    db.Accession = 12345
    assert str(db) == "12345"


def test_int(db):
    # DbRef int should be the same as the accession int
    with pytest.raises(ValueError):
        int(db)
    db.Accession = "12345"
    assert int(db) == 12345


def test_cmp():
    # DbRef cmp should first try numeric, then alphabetic, cmp
    assert DbRef("abc") < DbRef("xyz")
    assert DbRef("abc") == DbRef("abc")
    assert DbRef("123") > DbRef("14")
    assert DbRef("123") < DbRef("abc")
    # check that it ignores other attributes
    assert DbRef("x", "y", "z", "a", "b") == DbRef("x")


def test_make_list():
    # _make_list should always return a list
    assert _make_list("abc") == ["abc"]
    assert _make_list([]) == []
    assert _make_list(None) == [None]
    assert _make_list({"x": "y"}) == [{"x": "y"}]
    assert _make_list([1, 2, 3]) == [1, 2, 3]


def test_dbrefs_init_empty():
    # DbRefs empty init should work as expected
    assert DbRefs() == {}


def test_dbrefs_init_data():
    # DbRefs init with data should produce expected results
    d = DbRefs({"GenBank": "ab", "GO": (3, 44), "PDB": ["asdf", "ghjk"]})
    assert d == {"GenBank": ["ab"], "GO": [3, 44], "PDB": ["asdf", "ghjk"]}
    d.GenBank = "xyz"
    assert d["GenBank"] == ["xyz"]


def test_info_init_empty():
    # Info empty init should work as expected
    d = Info()
    assert len(d) == 1
    assert "Refs" in d
    assert d.Refs == DbRefs()
    assert isinstance(d.Refs, DbRefs)


def test_info_init_data():
    # Info init with data should put items in correct places
    # need to check init, setting, and resetting of attributes that belong
    # in the Info object and attributes that belong in Info.Refs. Also need
    # to check __getitem__, __setitem__, and __contains__.
    d = Info({"x": 3, "GO": 12345})
    assert d.x == 3
    assert d.GO == [12345]
    assert d.Refs.GO == [12345]
    with pytest.raises(AttributeError):
        del d.Refs
    d.GenBank = ("qaz", "wsx")
    assert d.GenBank == ["qaz", "wsx"]
    assert "GenBank" in d.Refs
    assert "GenBank" in d
    d.GenBank = "xyz"
    assert d.GenBank == ["xyz"]
    assert d.GenBank is d.Refs.GenBank
    d.GO = "x"
    assert d.GO == ["x"]
    d.GO.append("y")
    assert d.GO == ["x", "y"]
    d.ZZZ = "zzz"
    assert d.ZZZ == "zzz"
    assert "ZZZ" not in d.Refs
    assert "XXX" not in d
    assert d.XXX is None


def test_identity():
    # Info should get its own new Refs when created
    i = Info()
    j = Info()
    assert i is not j
    assert i.Refs is not j.Refs


def test_update():
    # update should warn the user of overlapping keys
    with warnings.catch_warnings(record=True) as w:
        d1 = Info({"key1": "value1", "key2": "value2", "key3": "value3"})
        d2 = Info({"key2": "value2", "key3": "value3", "key4": "value4"})
        d1.update(d2)
        assert len(w) == 1
