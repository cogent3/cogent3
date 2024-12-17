import os
import pathlib

import pytest

import cogent3
from cogent3.app import evo as evo_app
from cogent3.app.data_store import DataMember
from cogent3.app.result import (
    generic_result,
    hypothesis_result,
    model_collection_result,
    model_result,
    tabular_result,
)
from cogent3.util.deserialise import deserialise_object
from cogent3.util.dict_array import DictArray

_NEW_TYPE = "COGENT3_NEW_TYPE" in os.environ


def test_deserialised_values():
    """correctly deserialises values"""
    DNA = cogent3.get_moltype("dna")
    data = {"type": "cogent3.core.moltype.MolType", "moltype": "dna"}
    result = generic_result(source="blah.json")
    result["key"] = data
    result.deserialised_values()
    got = result["key"]
    assert got is DNA
    # if we have a type value without "cogent3", leaves as is
    data = {"type": "core.moltype.MolType", "moltype": "dna"}
    result = generic_result(source="blah.json")
    result["key"] = data
    result.deserialised_values()
    got = result["key"]
    assert got == data
    # or if no "type" entry, leaves as is
    data = {"moltype": "dna"}
    result = generic_result(source="blah.json")
    result["key"] = data
    result.deserialised_values()
    got = result["key"]
    assert got == data


def test_repr_str():
    """it works"""
    data = {"type": "cogent3.core.moltype.MolType", "moltype": "dna"}
    result = generic_result(source="blah.json")
    result["key"] = data
    repr(result)
    str(result)


def test_keys():
    """it works"""
    data = {"type": "cogent3.core.moltype.MolType", "moltype": "dna"}
    result = generic_result(source="blah.json")
    result["key"] = data
    keys = result.keys()
    assert keys == ["key"]


def test_invalid_setitem():
    """generic_result raise TypeError if trying to set invalid item type for json"""
    gr = generic_result("null")
    with pytest.raises(TypeError):
        gr["null"] = {0, 23}


def test_infers_source():
    """flexible handling of data source"""
    # works for string
    source = pathlib.Path("path/blah.fasta")
    aln = cogent3.make_aligned_seqs(
        {"A": "ACGT"},
        info={"source": str(source), "random_key": 1234},
        moltype="dna",
    )
    gr = generic_result(aln)
    assert gr.source == source.name

    # or Path
    aln.info.source = source
    gr = generic_result(aln)
    assert str(gr.source) == source.name

    # or DataMember
    aln.info.source = DataMember(data_store=None, unique_id=source.name)
    gr = generic_result(aln)
    assert str(gr.source) == source.name

    if _NEW_TYPE:
        aln.source = None
    else:
        aln.info = {}
    with pytest.raises(ValueError):
        generic_result(aln)


def test_model_repr():
    """does not fail"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = cogent3.make_aligned_seqs(data=_data, moltype="dna")
    mod = evo_app.model(
        "F81",
        show_progress=False,
        opt_args={"max_evaluations": 1, "limit_action": "ignore"},
    )
    result = mod(aln)
    assert isinstance(repr(result), str)
    # no values set
    assert isinstance(repr(model_result(source="blah")), str)


def test_model_result_alignment():
    """returns alignment from lf"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = cogent3.make_aligned_seqs(data=_data, moltype="dna")
    mod = evo_app.model(
        "F81",
        show_progress=False,
        opt_args={"max_evaluations": 5, "limit_action": "ignore"},
    )
    result = mod(aln)
    got = result.alignment
    assert got.to_dict() == _data


def test_model_name_lf_name():
    """model_result.name is set as lf.name"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = cogent3.make_aligned_seqs(data=_data, moltype="dna")
    mod = evo_app.model(
        "F81",
        name="blah",
        show_progress=False,
        opt_args={"max_evaluations": 5, "limit_action": "ignore"},
    )
    result = mod(aln)
    assert result.name == result.lf.name


def test_model_result_alignment_split_pos_model():
    """returns alignment from lf with split codon positions"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = cogent3.make_aligned_seqs(data=_data, moltype="dna")
    mod = evo_app.model(
        "F81",
        split_codons=True,
        show_progress=False,
        opt_args={"max_evaluations": 5, "limit_action": "ignore"},
    )
    result = mod(aln)
    for i in range(1, 4):
        got = result.alignment[i]
        expect = aln[i - 1 :: 3]
        assert got.to_dict() == expect.to_dict()


def test_model_result_repr_split_pos_model():
    """repr works for model_result of split codon positions"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = cogent3.make_aligned_seqs(data=_data, moltype="dna")
    mod = evo_app.model(
        "F81",
        split_codons=True,
        show_progress=False,
        opt_args={"max_evaluations": 55, "limit_action": "ignore"},
    )
    result = mod(aln)
    repr(result)


def test_model_result_tree_split_pos_model():
    """returns tree from lf with split codon positions"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = cogent3.make_aligned_seqs(data=_data, moltype="dna")
    mod = evo_app.model(
        "F81",
        split_codons=True,
        show_progress=False,
        opt_args={"max_evaluations": 55, "limit_action": "ignore"},
    )
    result = mod(aln)
    assert len(result.tree) == 3
    # check the trees are different by summing lengths
    lengths = {t.total_length() for _, t in result.tree.items()}
    assert len(lengths) > 1


def test_model_result_simulate_alignment():
    """returns tree from lf with split codon positions"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = cogent3.make_aligned_seqs(data=_data, moltype="dna")
    mod = evo_app.model(
        "F81",
        split_codons=True,
        show_progress=False,
        opt_args={"max_evaluations": 55, "limit_action": "ignore"},
    )
    result = mod(aln)
    got = result.simulate_alignment()
    assert len(aln) == len(got)
    assert aln.to_dict() != got.to_dict()


def test_model_result_tree_discrete_time():
    """returns paralinear lengths"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = cogent3.make_aligned_seqs(data=_data, moltype="dna")
    model1 = evo_app.model(
        "BH",
        opt_args={"max_evaluations": 25, "limit_action": "ignore"},
    )
    result = model1(aln)
    got = result.tree
    assert got.children[0].params["length"] == got.children[0].params["paralinear"]


def test_model_result_setitem():
    """TypeError if value a likelihood function, or a dict with correct type"""
    v = {"type": "arbitrary"}
    r = model_result(name="one", source="two")
    with pytest.raises(TypeError):
        r["name"] = v

    with pytest.raises(TypeError):
        r["name"] = 4

    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = cogent3.make_aligned_seqs(data=_data, moltype="dna")
    with pytest.raises(TypeError):
        r["name"] = aln


def test_model_result_repr_str():
    """it works even when no values"""
    mr = model_result(source="blah")
    assert isinstance(repr(mr), str)


def test_model_result_invalid_setitem():
    """model_result raise TypeError if trying to set incorrect item type"""
    mr = model_result(source="blah")
    with pytest.raises(TypeError):
        mr["null"] = 23


@pytest.fixture
def model_results():
    """fixture providing model results for collection tests"""
    if hasattr(model_results, "_cache"):
        return model_results._cache

    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = cogent3.make_aligned_seqs(data=_data, moltype="dna")
    model1 = evo_app.model(
        "F81",
        opt_args={"max_evaluations": 25, "limit_action": "ignore"},
    )
    model2 = evo_app.model(
        "HKY85",
        opt_args={"max_evaluations": 25, "limit_action": "ignore"},
    )
    model3 = evo_app.model(
        "GTR",
        opt_args={"max_evaluations": 25, "limit_action": "ignore"},
    )
    mr1 = model1(aln)
    mr2 = model2(aln)
    mr3 = model3(aln)
    results = {mr1.name: mr1, mr2.name: mr2, mr3.name: mr3}
    model_results._cache = results
    return results


def test_get_best_model(model_results):
    """should correctly identify the best model"""
    coll = model_collection_result(source="blah")
    coll.update(model_results)
    got = coll.get_best_model()
    # we ensure a model_result instance is returned from the possible set
    assert got in model_results.values()


def test_select_model(model_results):
    """correctly select models"""
    # we ensure a series of model_result instances is returned
    coll = model_collection_result(source="blah")
    coll.update(model_results)
    got = coll.select_models()
    assert len(got) > 0
    possible = list(model_results.values())
    for m in got:
        assert m in possible


def test_model_collection_result_repr(model_results):
    """constructed result can do the different repr"""
    result = model_collection_result(source="blah")
    coll = model_collection_result(source="blah")
    coll.update(model_results)
    got = result.__repr__()
    assert isinstance(got, str)
    got = result._repr_html_()
    assert isinstance(got, str)


def test_json_roundtrip(model_results):
    """roundtrip from json correct"""
    coll = model_collection_result(name="blah", source="blah2")
    coll.update(model_results)
    assert coll.name == "blah"
    assert coll.source == "blah2"
    orig = coll.__repr__()
    got = deserialise_object(coll.to_json())
    assert got.__repr__() == orig
    assert isinstance(got, model_collection_result)
    assert got.name == coll.name
    assert got.source == coll.source
    # select_models() should not fail
    got = deserialise_object(coll.to_json())
    m = got.select_models()
    assert isinstance(m[0], model_result)


def test_to_hypothesis(model_results):
    """creates a hypothesis_result from two model results"""
    mr = model_collection_result(source="blah")
    mr.update(model_results)
    hyp = mr.get_hypothesis_result("F81", "HKY85")
    assert isinstance(hyp, hypothesis_result)
    assert hyp.null.name == "F81"


def test_model_collection_repr_str():
    """it works even when no values"""
    mr = model_collection_result(source="blah")
    assert isinstance(repr(mr), str)


def test_model_collection_result_invalid_setitem():
    """model_collection_result raise TypeError if trying to set incorrect item type"""
    mcr = model_collection_result(source="blah")
    with pytest.raises(TypeError):
        mcr["null"] = 23


def test_hypothesis_repr_str():
    """it works even when no values"""
    hr = hypothesis_result(name_of_null="null", source="blah")
    assert isinstance(repr(hr), str)


def test_hypothesis_pvalue():
    """hypothesis test p-value property"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = cogent3.make_aligned_seqs(data=_data, moltype="dna")
    model1 = evo_app.model(
        "F81",
        opt_args={"max_evaluations": 25, "limit_action": "ignore"},
    )
    model2 = evo_app.model(
        "HKY85",
        opt_args={"max_evaluations": 25, "limit_action": "ignore"},
    )
    hyp = evo_app.hypothesis(model1, model2)
    result = hyp(aln)
    assert 0 <= result.pvalue <= 1


def test_hypothesis_invalid_setitem():
    """hypothesis_result raise TypeError if trying to set incorrect item type"""
    hr = hypothesis_result(name_of_null="null", source="blah")
    with pytest.raises(TypeError):
        hr["null"] = {0, 23}


def test_tabular_valid_setitem():
    """tabular_result works when set correct item type"""
    tr = tabular_result("null")
    tr["result"] = cogent3.make_table(data={"A": [0, 1]})
    darr = DictArray({"A": [0, 1]})
    tr["result2"] = darr
    js = tr.to_json()
    assert isinstance(js, str)


def test_tabular_invalid_setitem():
    """tabular_result raise TypeError if trying to set incorrect item type"""
    tr = tabular_result("null")
    with pytest.raises(TypeError):
        tr["null"] = {0, 23}
