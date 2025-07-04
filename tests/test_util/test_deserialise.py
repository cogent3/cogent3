import json

import numpy
import pytest
from numpy.testing import assert_allclose

from cogent3 import (
    get_app,
    get_code,
    get_moltype,
    load_aligned_seqs,
    load_tree,
    make_aligned_seqs,
    make_tree,
    make_unaligned_seqs,
)
from cogent3.app.result import model_collection_result, model_result
from cogent3.evolve.models import get_model
from cogent3.evolve.ns_substitution_model import (
    NonReversibleDinucleotide,
    NonReversibleTrinucleotide,
    _sym_preds,
)
from cogent3.util.deserialise import (
    deserialise_likelihood_function,
    deserialise_object,
)


def test_roundtrip_codon_alphabet():
    """codon alphabet to_json enables roundtrip"""
    standard = get_code(1).get_alphabet()
    data = standard.to_json()
    got = deserialise_object(data)
    assert isinstance(got, type(standard))
    assert list(got) == list(standard)


def test_roundtrip_alphabet():
    """alphabet to_json enables roundtrip"""
    dna = get_moltype("dna")
    data = dna.alphabet.to_json()
    got = deserialise_object(data)
    assert isinstance(got, type(dna.alphabet))
    assert list(got) == list(dna.alphabet)


def test_roundtrip_moltype():
    """moltype to_json enables roundtrip"""
    dna = get_moltype("dna")
    data = dna.to_json()
    got = deserialise_object(data)
    assert isinstance(got, type(dna))
    assert got.name == dna.name


def test_roundtrip_seqcoll():
    """SequenceCollection to_json enables roundtrip"""
    data = {"A": "TTGT", "B": "GGCT"}
    seqcoll = make_unaligned_seqs(data, moltype="dna")
    got = deserialise_object(seqcoll.to_json())
    assert got.rc().to_dict() == seqcoll.rc().to_dict()
    assert got.__class__.__name__ == "SequenceCollection"


def test_roundtrip_align():
    """Alignment to_json enables roundtrip"""
    data = {"A": "TTGTA", "B": "GGCT-"}
    align = make_aligned_seqs(data, moltype="dna")
    got = deserialise_object(align.to_json())
    assert got.rc().to_dict() == align.rc().to_dict()
    assert got.__class__.__name__.endswith("Alignment")


def test_roundtrip_tree():
    """Tree to_json enables roundtrip"""
    tree = make_tree(treestring="(c:01,d:0.3,(a:0.05,b:0.08)xx:0.2)")
    got = deserialise_object(tree.to_json())
    assert_allclose(got.get_node_matching_name("a").length, 0.05)
    assert_allclose(got.get_node_matching_name("xx").length, 0.2)


@pytest.mark.parametrize("model_name", ["HKY85", "GN", "CNFGTR", "GNC"])
def test_roundtrip_submod(model_name):
    """substitution model to_json enables roundtrip"""
    sm = get_model(model_name)
    data = sm.to_json()
    got = deserialise_object(data)
    assert got.get_param_list() == sm.get_param_list()
    assert isinstance(got, type(sm))


def test_roundtrip_non_reversible_submod():
    sm = NonReversibleDinucleotide(_sym_preds)
    data = sm.to_json()
    got = deserialise_object(data)
    assert got.get_param_list() == sm.get_param_list()
    assert isinstance(got, type(sm))

    sm = NonReversibleTrinucleotide(_sym_preds)
    data = sm.to_json()
    got = deserialise_object(data)
    assert got.get_param_list() == sm.get_param_list()
    assert isinstance(got, type(sm))


@pytest.mark.parametrize("motif_length", [1, 2])
def test_roundtrip_discrete_time_submod(motif_length):
    """discrete time substitution models to_json enables roundtrip"""
    sm = get_model("DT", motif_length=motif_length)
    assert sm.word_length == motif_length
    data = sm.to_rich_dict()
    got = deserialise_object(data)
    assert isinstance(got, type(sm))
    assert got.get_param_list() == sm.get_param_list()
    assert got.word_length == motif_length


def test_roundtrip_likelihood_function():
    """likelihood function.to_json enables roundtrip"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = make_aligned_seqs(_data, moltype="dna")
    tree = make_tree(tip_names=aln.names)
    sm = get_model("HKY85")
    lf = sm.make_likelihood_function(tree)
    lf.set_alignment(aln)
    edge_vals = zip(aln.names, (2, 3, 4), strict=False)
    for edge, val in edge_vals:
        lf.set_param_rule("kappa", edge=edge, init=val)
    lnL = lf.get_log_likelihood()
    data = lf.to_json()
    got_obj = deserialise_object(data)
    assert_allclose(got_obj.get_log_likelihood(), lnL)


def test_roundtrip_discrete_time_likelihood_function():
    """discrete time likelihood function.to_json enables roundtrip"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = make_aligned_seqs(_data, moltype="dna")
    tree = make_tree(tip_names=aln.names)
    sm = get_model("BH")
    lf = sm.make_likelihood_function(tree)
    lf.set_alignment(aln)
    lf.optimise(max_evaluations=25, limit_action="ignore", show_progress=False)
    lnL = lf.get_log_likelihood()
    data = lf.to_json()
    got_obj = deserialise_object(data)
    assert_allclose(got_obj.get_log_likelihood(), lnL)


def test_roundtrip_het_lf(DATA_DIR):
    """correctly round trips a site-het model"""
    with open(DATA_DIR / "site-het-param-rules.json") as infile:
        rules = json.load(infile)

    aln = load_aligned_seqs(DATA_DIR / "primates_brca1.fasta", moltype="dna")
    tree = load_tree(DATA_DIR / "primates_brca1.tree")
    rule_lnL = rules.pop("phylohmm-gamma-kappa")
    sm = get_model("HKY85", ordered_param="rate", distribution="gamma")
    lf1 = sm.make_likelihood_function(tree, bins=4, sites_independent=False)
    lf1.set_alignment(aln)
    lf1.apply_param_rules(rule_lnL["rules"])
    data = lf1.to_json()
    got_lf = deserialise_object(data)
    assert_allclose(lf1.lnL, got_lf.lnL)


def test_roundtrip_from_file(tmp_path):
    """correctly roundtrips a likelihood function fro json file"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = make_aligned_seqs(_data, moltype="dna")
    tree = make_tree(tip_names=aln.names)
    sm = get_model("HKY85")
    lf = sm.make_likelihood_function(tree)
    lf.set_alignment(aln)
    edge_vals = zip(aln.names, (2, 3, 4), strict=False)
    for edge, val in edge_vals:
        lf.set_param_rule("kappa", edge=edge, init=val)
    lnL = lf.get_log_likelihood()
    data = lf.to_json()

    outpath = tmp_path / "delme.json"
    with open(outpath, "w") as outfile:
        outfile.write(data)

    got = deserialise_object(outpath)
    assert_allclose(got.get_log_likelihood(), lnL)


def test_roundtrip_model_result():
    """mode_result.to_json enables roundtrip and lazy evaluation"""
    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = make_aligned_seqs(_data, moltype="dna")
    tree = make_tree(tip_names=aln.names)
    sm = get_model("HKY85")
    lf = sm.make_likelihood_function(tree)
    lf.set_alignment(aln)
    edge_vals = zip(aln.names, (2, 3, 4), strict=False)
    for edge, val in edge_vals:
        lf.set_param_rule("kappa", edge=edge, init=val)
    result = model_result(name="test", source="blah")
    result[1] = lf
    assert result[1] is lf
    assert result.nfp == lf.nfp
    assert result.lnL == lf.lnL

    data = result.to_json()
    got_obj = deserialise_object(data)
    # lazy evaluation means initially, the value is a dict
    assert isinstance(got_obj[1], dict)
    # and properties match original
    assert got_obj.lnL == result.lnL
    assert got_obj.nfp == result.nfp
    assert got_obj.DLC == result.DLC
    # when we ask for the lf attribute, it's no longer a dict
    assert not isinstance(got_obj.lf, dict)
    assert got_obj.lf.nfp == got_obj.nfp


def test_roundtrip_model_result2():
    """model_result of split codon correct type after roundtrip"""
    from cogent3.app import evo as evo_app
    from cogent3.evolve.parameter_controller import (
        AlignmentLikelihoodFunction,
    )

    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = make_aligned_seqs(_data, moltype="dna")
    opt_args = {"max_evaluations": 10, "limit_action": "ignore"}
    m1 = evo_app.model("F81", split_codons=True, opt_args=opt_args)
    result = m1(aln)

    data = result.to_json()
    got_obj = deserialise_object(data)
    for i in range(1, 4):
        assert isinstance(got_obj[i], dict)

    # after accessing attribute, should be automatically inflated
    _ = got_obj.lf
    for i in range(1, 4):
        assert isinstance(got_obj[i], AlignmentLikelihoodFunction)

    # or after using the deserialise method
    data = result.to_json()
    got_obj = deserialise_object(data)
    got_obj.deserialised_values()
    for i in range(1, 4):
        assert isinstance(got_obj[i], AlignmentLikelihoodFunction)


def test_model_collection_result():
    """round trip of model collection works"""
    from cogent3.app import evo as evo_app
    from cogent3.evolve.parameter_controller import (
        AlignmentLikelihoodFunction,
    )

    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = make_aligned_seqs(_data, moltype="dna")
    opt_args = {"max_evaluations": 10, "limit_action": "ignore"}
    m1 = evo_app.model("F81", split_codons=True, opt_args=opt_args)
    m2 = evo_app.model("GTR", split_codons=True, opt_args=opt_args)
    models = (m1, m2)
    mc_result = model_collection_result(name="collection", source="blah")
    for model in models:
        mc_result[model.name] = model(aln)

    for model in models:
        for i in range(1, 4):
            assert isinstance(mc_result[model.name][i], AlignmentLikelihoodFunction)

    data = mc_result.to_json()
    got_obj = deserialise_object(data)
    for model in models:
        for i in range(1, 4):
            assert isinstance(got_obj[model.name][i], dict)

    # but after invoking deserialised_values
    got_obj.deserialised_values()
    for model in models:
        for i in range(1, 4):
            assert isinstance(got_obj[model.name][i], AlignmentLikelihoodFunction)


def test_roundtrip_hypothesis_result():
    """nested items retain the correct type after roundtrip"""
    from cogent3.app import evo as evo_app
    from cogent3.evolve.parameter_controller import (
        AlignmentLikelihoodFunction,
    )

    _data = {
        "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = make_aligned_seqs(_data, moltype="dna")
    opt_args = {"max_evaluations": 10, "limit_action": "ignore"}
    m1 = evo_app.model("F81", split_codons=True, opt_args=opt_args)
    m2 = evo_app.model("GTR", split_codons=True, opt_args=opt_args)
    hyp = evo_app.hypothesis(m1, m2)
    result = hyp(aln)
    assert isinstance(result["F81"][1], AlignmentLikelihoodFunction)

    data = result.to_json()
    got_obj = deserialise_object(data)
    for i in range(1, 4):
        for sm in ("F81", "GTR"):
            assert isinstance(got_obj[sm][i], dict)

    # but after invoking  deserialised_values
    got_obj.deserialised_values()
    for i in range(1, 4):
        for sm in ("F81", "GTR"):
            assert isinstance(got_obj[sm][i], AlignmentLikelihoodFunction)


def test_roundtrip_tuple_key():
    """deserialise_result handles tuples as keys"""
    from cogent3.app.result import generic_result

    r = generic_result(source="none")
    r[(1, 2)] = 24
    got = deserialise_object(r.to_json())
    assert got[1, 2] == 24


def test_not_completed_result():
    """correctly reconstructs a NotCompletedResult object"""
    from cogent3.app.composable import NotCompleted

    val = NotCompleted("ERROR", "nothing", "some error", source="here")
    expect = val.to_rich_dict()
    json = val.to_json()
    got = deserialise_object(json)
    assert got.to_rich_dict() == expect


def test_deserialise_tabular_table():
    """correctly deserialises Table"""
    from cogent3 import make_table

    table = make_table(
        header=["id", "foo", "bar"],
        rows=[
            [1, "abc", 11],
            [2, "bca", 22],
            [3, "cab", 33],
            [4, "abc", 44],
            [5, "bca", 55],
        ],
    )
    json = table.to_json()
    got = deserialise_object(json)
    assert got.to_dict() == table.to_dict()


def test_deserialise_tabular_dictarray():
    """correctly deserialises DictArray"""
    from cogent3.util.dict_array import DictArrayTemplate

    template = DictArrayTemplate(5, ["id", "foo", "bar"])
    data = [
        [1, "abc", 11],
        [2, "bca", 22],
        [3, "cab", 33],
        [4, "abc", 44],
        [5, "bca", 55],
    ]
    darr = template.wrap(data)
    json = darr.to_json()
    got = deserialise_object(json)
    assert got.to_dict() == darr.to_dict()


def test_deserialise_tabular_distancematrix():
    """correctly deserialises DistanceMatrix"""
    from cogent3.evolve.fast_distance import DistanceMatrix

    data = {
        ("ABAYE2984", "Atu3667"): None,
        ("ABAYE2984", "Avin_42730"): 0.638,
        ("ABAYE2984", "BAA10469"): None,
        ("Atu3667", "ABAYE2984"): None,
        ("Atu3667", "Avin_42730"): 2.368,
        ("Atu3667", "BAA10469"): None,
        ("Avin_42730", "ABAYE2984"): 0.638,
        ("Avin_42730", "Atu3667"): 2.368,
        ("Avin_42730", "BAA10469"): 1.85,
        ("BAA10469", "ABAYE2984"): None,
        ("BAA10469", "Atu3667"): None,
        ("BAA10469", "Avin_42730"): 1.85,
    }

    dm = DistanceMatrix(data)
    json = dm.to_json()
    got = deserialise_object(json)
    dm_dict = dm.to_dict()
    got_dict = got.to_dict()
    for (a, b), dist in dm_dict.items():
        if dist is None:
            assert numpy.isnan(got_dict[a, b])
        else:
            assert_allclose(dist, got_dict[a, b])


def test_deserialise_python_builtins():
    """any object that does not contain a type key is returned as is"""
    data = {"a": 123, "b": "text"}
    jdata = json.dumps(data)
    got = deserialise_object(jdata)
    assert got == data
    data = range(4)
    got = deserialise_object(data)
    assert got is data


def test_deserialise_likelihood_function1(DATA_DIR):
    """correctly deserialise data into likelihood function"""
    # tests single alignment
    aln = load_aligned_seqs(
        filename=DATA_DIR / "brca1_5.paml",
        moltype="dna",
    )
    tree = make_tree(tip_names=aln.names)
    model = get_model("HKY85")
    lf = model.make_likelihood_function(tree)
    lf.set_alignment(aln)
    lf_rich_dict = lf.to_rich_dict()
    got = deserialise_likelihood_function(lf_rich_dict)
    assert str(lf.defn_for["mprobs"]) == str(got.defn_for["mprobs"])
    assert str(lf.defn_for["alignment"].assignments) == str(
        got.defn_for["alignment"].assignments,
    )


def test_deserialise_likelihood_function_multilocus(DATA_DIR):
    """correctly deserialise data of multilocus likelihood function"""
    # tests multiple alignments
    data = load_aligned_seqs(
        filename=DATA_DIR / "brca1_5.paml",
        moltype="dna",
    )
    half = len(data) // 2
    aln1 = data[:half]
    aln2 = data[half:]
    loci_names = ["1st-half", "2nd-half"]
    loci = [aln1, aln2]
    tree = make_tree(tip_names=data.names)
    model = get_model("HKY85", optimise_motif_probs=True)
    lf = model.make_likelihood_function(tree, loci=loci_names)
    lf.set_alignment(loci)
    lf_rich_dict = lf.to_rich_dict()
    got = deserialise_likelihood_function(lf_rich_dict)
    assert str(lf.defn_for["mprobs"]) == str(got.defn_for["mprobs"])
    assert str(lf.defn_for["alignment"].assignments) == str(
        got.defn_for["alignment"].assignments,
    )
    # now constrain mprobs to be the same
    lf.set_param_rule("mprobs", is_independent=False)
    lf_rich_dict = lf.to_rich_dict()
    got = deserialise_likelihood_function(lf_rich_dict)
    assert str(lf.defn_for["mprobs"]) == str(got.defn_for["mprobs"])
    assert str(lf.defn_for["alignment"].assignments) == str(
        got.defn_for["alignment"].assignments,
    )


def test_custom_deserialiser():
    """correctly registers a function to inflate a custom object"""
    from cogent3.util.deserialise import register_deserialiser

    @register_deserialiser("myfunkydata")
    def astuple(data):
        data.pop("type")
        return tuple(data["data"])

    orig = {"type": "myfunkydata", "data": (1, 2, 3)}
    txt = json.dumps(orig)
    got = deserialise_object(txt)
    assert got == (1, 2, 3)
    assert isinstance(got, tuple)

    with pytest.raises(TypeError):

        @register_deserialiser
        def astupled(data):
            data.pop("type")
            return tuple(data["data"])


def test_convert_annotation_to_annotation_db():
    from cogent3.core.annotation_db import (
        BasicAnnotationDb,
        convert_annotation_to_annotation_db,
    )

    data = {
        "name": "A",
        "data": [
            {
                "annotation_construction": {
                    "type": "gene",
                    "name": "n1",
                    "map": {
                        "spans": [
                            {
                                "start": 2,
                                "end": 5,
                                "value": None,
                                "reverse": False,
                                "type": "cogent3.core.location.Span",
                                "version": "2023.2.12a1",
                            },
                        ],
                        "parent_length": 5,
                        "type": "cogent3.core.location.FeatureMap",
                        "version": "2023.2.12a1",
                    },
                },
                "type": "cogent3.core.annotation.AnnotatableFeature",
                "version": "2023.2.12a1",
            },
        ],
    }
    data = convert_annotation_to_annotation_db(data)
    # which can be turned back into an annotation db
    db = deserialise_object(data)
    assert isinstance(db, BasicAnnotationDb)
    assert db.num_matches() == 1


@pytest.mark.parametrize("rate_matrix_required", [True, False])
def test_roundtrip_TN93_model(rate_matrix_required):
    """model_result of split codon correct type after roundtrip"""
    _data = {
        "a": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "b": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "c": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = make_aligned_seqs(_data, moltype="dna")
    tree = make_tree(tip_names=["a", "b", "c"])
    tn93 = get_model(
        "TN93",
        rate_matrix_required=rate_matrix_required,
    ).make_likelihood_function(tree)
    tn93.set_alignment(aln)

    got = deserialise_object(tn93.to_rich_dict())
    assert_allclose(got.lnL, tn93.lnL)


def test_roundtrip_TN93_model_result():
    _data = {
        "a": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
        "b": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
        "c": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
    }
    aln = make_aligned_seqs(_data, moltype="dna")
    tn93 = get_app("model", "TN93")
    result = tn93(aln)

    got = deserialise_object(result.to_rich_dict())
    assert_allclose(got.lnL, result.lnL)


@pytest.mark.parametrize("mtype", ["dna", "protein"])
def test_roundtrip_seq(mtype):
    """seq to_json enables roundtrip"""
    mtype = get_moltype(mtype)
    seq = mtype.make_seq(seq="ACGGTCGG", name="label", info={"something": 3})
    got = deserialise_object(seq.to_json())
    assert got.info.something == 3
    assert got.name == "label"
    assert got.moltype == seq.moltype
    assert str(got) == str(seq)


@pytest.mark.parametrize("mn", ["BH", "F81", "JC69"])
def test_dser_submodel(mn):
    sm = get_model(mn)
    rd = sm.to_rich_dict(for_pickle=False)
    got = deserialise_object(rd)
    assert got.moltype.alphabet.moltype is got.moltype


def test_deserialise_old_to_new_type_alignment_1():
    from cogent3.util.misc import get_object_provenance

    rd = {
        "seqs": {
            "Human": {
                "name": "Human",
                "seq": {
                    "type": "cogent3.core.sequence.SeqView",
                    "version": "2025.5.8a9",
                    "init_args": {
                        "seq": "TGTGGCACAAATACTCATGC",
                        "seqid": "Human",
                        "step": 1,
                    },
                },
                "moltype": "dna",
                "info": None,
                "type": "cogent3.core.sequence.DnaSequence",
                "version": "2025.5.8a9",
                "annotation_offset": 0,
            },
            "Wombat": {
                "name": "Wombat",
                "seq": {
                    "type": "cogent3.core.sequence.SeqView",
                    "version": "2025.5.8a9",
                    "init_args": {
                        "seq": "---------------NNTGC",
                        "seqid": "Wombat",
                        "step": 1,
                    },
                },
                "moltype": "dna",
                "info": None,
                "type": "cogent3.core.sequence.DnaSequence",
                "version": "2025.5.8a9",
                "annotation_offset": 0,
            },
            "Mouse": {
                "name": "Mouse",
                "seq": {
                    "type": "cogent3.core.sequence.SeqView",
                    "version": "2025.5.8a9",
                    "init_args": {
                        "seq": "TGTGGCACAGATGCTCATGC",
                        "seqid": "Mouse",
                        "step": 1,
                    },
                },
                "moltype": "dna",
                "info": None,
                "type": "cogent3.core.sequence.DnaSequence",
                "version": "2025.5.8a9",
                "annotation_offset": 0,
            },
        },
        "moltype": "dna",
        "info": {"source": "unknown"},
        "type": "cogent3.core.alignment.ArrayAlignment",
        "version": "2025.5.8a9",
    }
    got = deserialise_object(rd)
    got_prov = get_object_provenance(got)
    assert got_prov == "cogent3.core.alignment.Alignment"
    rd = {
        "info": {"seed": "758_443154_73021", "source": "15", "fg_edge": "73021"},
        "moltype": "dna",
        "seqs": {
            "443154": {
                "info": None,
                "moltype": "dna",
                "name": "443154",
                "seq": "CAAGCATAATAA",
                "type": "cogent3.core.sequence.DnaSequence",
                "version": "2022.10.31a1",
            },
            "73021": {
                "info": None,
                "moltype": "dna",
                "name": "73021",
                "seq": "CAAGTATCCTAA",
                "type": "cogent3.core.sequence.DnaSequence",
                "version": "2022.10.31a1",
            },
            "758": {
                "info": None,
                "moltype": "dna",
                "name": "758",
                "seq": "TGAATATCATAA",
                "type": "cogent3.core.sequence.DnaSequence",
                "version": "2022.10.31a1",
            },
        },
        "type": "cogent3.core.alignment.ArrayAlignment",
        "version": "2022.10.31a1",
    }
    got = deserialise_object(rd)
    got_prov = get_object_provenance(got)
    assert got_prov == "cogent3.core.alignment.Alignment"
    assert got.source == "15"
    assert "source" not in got.info


def test_deserialise_old_to_new_type_alignment_2():
    from cogent3.util.misc import get_object_provenance

    rd = {
        "seqs": {
            "Human": {
                "version": "2025.5.8a9",
                "type": "cogent3.core.alignment.Aligned",
                "map_init": {
                    "gap_pos": [],
                    "cum_gap_lengths": [],
                    "termini_unknown": False,
                    "parent_length": 20,
                    "type": "cogent3.core.location.IndelMap",
                    "version": "2025.5.8a9",
                },
                "seq_init": {
                    "name": "Human",
                    "seq": {
                        "type": "cogent3.core.sequence.SeqView",
                        "version": "2025.5.8a9",
                        "init_args": {
                            "seq": "TGTGGCACAAATACTCATGC",
                            "seqid": "Human",
                            "step": 1,
                        },
                    },
                    "moltype": "dna",
                    "info": None,
                    "type": "cogent3.core.sequence.DnaSequence",
                    "version": "2025.5.8a9",
                    "annotation_offset": 0,
                },
            },
            "Wombat": {
                "version": "2025.5.8a9",
                "type": "cogent3.core.alignment.Aligned",
                "map_init": {
                    "gap_pos": [0],
                    "cum_gap_lengths": [15],
                    "termini_unknown": False,
                    "parent_length": 5,
                    "type": "cogent3.core.location.IndelMap",
                    "version": "2025.5.8a9",
                },
                "seq_init": {
                    "name": "Wombat",
                    "seq": {
                        "type": "cogent3.core.sequence.SeqView",
                        "version": "2025.5.8a9",
                        "init_args": {"seq": "NNTGC", "seqid": "Wombat", "step": 1},
                    },
                    "moltype": "dna",
                    "info": None,
                    "type": "cogent3.core.sequence.DnaSequence",
                    "version": "2025.5.8a9",
                    "annotation_offset": 0,
                },
            },
            "Mouse": {
                "version": "2025.5.8a9",
                "type": "cogent3.core.alignment.Aligned",
                "map_init": {
                    "gap_pos": [],
                    "cum_gap_lengths": [],
                    "termini_unknown": False,
                    "parent_length": 20,
                    "type": "cogent3.core.location.IndelMap",
                    "version": "2025.5.8a9",
                },
                "seq_init": {
                    "name": "Mouse",
                    "seq": {
                        "type": "cogent3.core.sequence.SeqView",
                        "version": "2025.5.8a9",
                        "init_args": {
                            "seq": "TGTGGCACAGATGCTCATGC",
                            "seqid": "Mouse",
                            "step": 1,
                        },
                    },
                    "moltype": "dna",
                    "info": None,
                    "type": "cogent3.core.sequence.DnaSequence",
                    "version": "2025.5.8a9",
                    "annotation_offset": 0,
                },
            },
        },
        "moltype": "dna",
        "info": {"source": "unknown"},
        "type": "cogent3.core.alignment.Alignment",
        "version": "2025.5.8a9",
    }
    got = deserialise_object(rd)
    got_prov = get_object_provenance(got)
    assert got_prov == "cogent3.core.alignment.Alignment"
    assert got.source == "unknown"
    assert "source" not in got.info


def test_deserialise_old_to_new_type_seqcoll():
    from cogent3.util.misc import get_object_provenance

    rd = {
        "seqs": {
            "Human": {
                "name": "Human",
                "seq": {
                    "type": "cogent3.core.sequence.SeqView",
                    "version": "2025.5.8a9",
                    "init_args": {
                        "seq": "TGTGGCACAAATACTCATGC",
                        "seqid": "Human",
                        "step": 1,
                    },
                },
                "moltype": "dna",
                "info": None,
                "type": "cogent3.core.sequence.DnaSequence",
                "version": "2025.5.8a9",
                "annotation_offset": 0,
            },
            "Wombat": {
                "name": "Wombat",
                "seq": {
                    "type": "cogent3.core.sequence.SeqView",
                    "version": "2025.5.8a9",
                    "init_args": {
                        "seq": "---------------NNTGC",
                        "seqid": "Wombat",
                        "step": 1,
                    },
                },
                "moltype": "dna",
                "info": None,
                "type": "cogent3.core.sequence.DnaSequence",
                "version": "2025.5.8a9",
                "annotation_offset": 0,
            },
            "Mouse": {
                "name": "Mouse",
                "seq": {
                    "type": "cogent3.core.sequence.SeqView",
                    "version": "2025.5.8a9",
                    "init_args": {
                        "seq": "TGTGGCACAGATGCTCATGC",
                        "seqid": "Mouse",
                        "step": 1,
                    },
                },
                "moltype": "dna",
                "info": None,
                "type": "cogent3.core.sequence.DnaSequence",
                "version": "2025.5.8a9",
                "annotation_offset": 0,
            },
        },
        "moltype": "dna",
        "info": {"source": "unknown"},
        "type": "cogent3.core.alignment.SequenceCollection",
        "version": "2025.5.8a9",
    }
    got = deserialise_object(rd)
    got_prov = get_object_provenance(got)
    assert got_prov == "cogent3.core.alignment.SequenceCollection"
    assert got.source == "unknown"
    assert "source" not in got.info


moltype_old_class_map = {
    "dna": "cogent3.core.sequence.DnaSequence",
    "rna": "cogent3.core.sequence.RnaSequence",
    "text": "cogent3.core.sequence.Sequence",
    "bytes": "cogent3.core.sequence.ByteSequence",
    "protein": "cogent3.core.sequence.ProteinSequence",
    "protein_with_stop": "cogent3.core.sequence.ProteinWithStopSequence",
}


@pytest.mark.parametrize("moltype_name", list(moltype_old_class_map.keys()))
def test_deserialise_old_to_new_type_seqs(moltype_name):
    from cogent3.core import moltype as c3_moltype
    from cogent3.core import sequence

    old_rich_dict = {
        "annotation_offset": 0,
        "info": {"something": 3},
        "moltype": moltype_name,
        "name": "label",
        "seq": {
            "init_args": {"seq": "ACGGTCGG", "seqid": "label", "step": 1},
            "type": "cogent3.core.sequence.SeqView",
            "version": "2025.5.8a9",
        },
        "type": moltype_old_class_map[moltype_name],
        "version": "2025.5.8a9",
    }
    got = deserialise_object(old_rich_dict)
    assert isinstance(got.moltype, c3_moltype.MolType)
    assert isinstance(got, sequence.Sequence)


def test_deserialise_old_to_new_type_moltype():
    from cogent3.core import moltype as c3_moltype

    old_rich_dict = {
        "type": "cogent3.core.moltype.MolType",
        "moltype": "dna",
        "version": "2025.5.8a9",
    }
    got = deserialise_object(old_rich_dict)
    assert isinstance(got, c3_moltype.MolType)


@pytest.mark.parametrize(
    "motifset", ["TCAG", "TCAGNRYWSKMBDHV", "TCAG-", "TCAG-NRYWSKMBDHV?"]
)
def test_deserialise_old_to_new_charalphabet(motifset):
    from cogent3.core import alphabet as c3_alphabet

    old_rich_dict = {
        "motifset": list(motifset),
        "gap": "-",
        "moltype": "dna",
        "type": "cogent3.core.alphabet.CharAlphabet",
        "version": "2025.5.8a9",
    }
    got = deserialise_object(old_rich_dict)
    assert isinstance(got, c3_alphabet.CharAlphabet)


def test_deserialise_old_to_c3_alphabet():
    # this is special to BH
    from cogent3.core import alphabet as c3_alphabet

    old_rich_dict = {
        "motifset": ("T", "C", "A", "G"),
        "gap": "-",
        "moltype": "dna",
        "type": "cogent3.core.alphabet.Alphabet",
        "version": "2025.5.8a9",
    }
    got = deserialise_object(old_rich_dict)
    assert isinstance(got, c3_alphabet.CharAlphabet)


def test_deserialise_old_to_new_type_kmer_alphabet():
    from cogent3.core import alphabet as c3_alphabet

    old_rich_dict = {
        "data": ["TCAG", "TCAG"],
        "gap": None,
        "moltype": "dna",
        "type": "cogent3.core.alphabet.JointEnumeration",
        "version": "2025.5.8a9",
    }
    got = deserialise_object(old_rich_dict)
    assert isinstance(got, c3_alphabet.KmerAlphabet)
    assert len(got) == 16


def test_deserialise_old_to_new_type_codon_alphabet():
    from cogent3.core import alphabet as c3_alphabet

    old_rich_dict = {
        "motifset": (
            "TTT",
            "TTC",
            "TTA",
            "TTG",
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "TAT",
            "TAC",
            "TGT",
            "TGC",
            "TGG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "CCT",
            "CCC",
            "CCA",
            "CCG",
            "CAT",
            "CAC",
            "CAA",
            "CAG",
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "ATT",
            "ATC",
            "ATA",
            "ATG",
            "ACT",
            "ACC",
            "ACA",
            "ACG",
            "AAT",
            "AAC",
            "AAA",
            "AAG",
            "AGT",
            "AGC",
            "AGA",
            "AGG",
            "GTT",
            "GTC",
            "GTA",
            "GTG",
            "GCT",
            "GCC",
            "GCA",
            "GCG",
            "GAT",
            "GAC",
            "GAA",
            "GAG",
            "GGT",
            "GGC",
            "GGA",
            "GGG",
        ),
        "gap": "-",
        "moltype": "dna",
        "type": "cogent3.core.alphabet.Alphabet",
        "version": "2025.5.8a9",
    }
    alpha = deserialise_object(old_rich_dict)
    assert isinstance(alpha, c3_alphabet.SenseCodonAlphabet)


def test_deserialise_old_to_new_type_submodel():
    from cogent3.evolve.ns_substitution_model import DiscreteSubstitutionModel
    from cogent3.evolve.substitution_model import TimeReversibleNucleotide

    old_rich_dict = {
        "alphabet": {
            "motifset": ("T", "C", "A", "G"),
            "gap": "-",
            "moltype": "dna",
            "type": "cogent3.core.alphabet.Alphabet",
            "version": "2025.5.8a9",
        },
        "motif_probs": None,
        "optimise_motif_probs": True,
        "equal_motif_probs": False,
        "motif_probs_from_data": None,
        "motif_probs_alignment": None,
        "mprob_model": "tuple",
        "model_gaps": False,
        "recode_gaps": False,
        "motif_length": 1,
        "name": "BH",
        "motifs": None,
        "type": "cogent3.evolve.ns_substitution_model.DiscreteSubstitutionModel",
        "version": "2025.5.8a9",
    }
    got = deserialise_object(old_rich_dict)
    assert isinstance(got, DiscreteSubstitutionModel)

    old_rich_dict = {
        "alphabet": {
            "motifset": ("T", "C", "A", "G"),
            "gap": "-",
            "moltype": "dna",
            "type": "cogent3.core.alphabet.CharAlphabet",
            "version": "2025.5.8a9",
        },
        "motif_probs": None,
        "optimise_motif_probs": False,
        "equal_motif_probs": False,
        "motif_probs_from_data": None,
        "motif_probs_alignment": None,
        "mprob_model": None,
        "model_gaps": False,
        "recode_gaps": True,
        "motif_length": 1,
        "name": "F81",
        "motifs": None,
        "with_rate": False,
        "ordered_param": None,
        "distribution": None,
        "partitioned_params": None,
        "kw": {},
        "predicates": [],
        "scales": None,
        "type": "cogent3.evolve.substitution_model.TimeReversibleNucleotide",
        "version": "2025.5.8a9",
    }

    got = deserialise_object(old_rich_dict)
    assert isinstance(got, TimeReversibleNucleotide)


def test_deserialise_old_to_new_type_likelihood():
    old_json = {
        "model": {
            "alphabet": {
                "motifset": ["T", "C", "A", "G"],
                "gap": "-",
                "moltype": "dna",
                "type": "cogent3.core.alphabet.Alphabet",
                "version": "2025.5.8a9",
            },
            "motif_probs": None,
            "optimise_motif_probs": True,
            "equal_motif_probs": False,
            "motif_probs_from_data": None,
            "motif_probs_alignment": None,
            "mprob_model": "tuple",
            "model_gaps": False,
            "recode_gaps": False,
            "motif_length": 1,
            "name": "BH",
            "motifs": None,
            "type": "cogent3.evolve.ns_substitution_model.DiscreteSubstitutionModel",
            "version": "2025.5.8a9",
        },
        "tree": {
            "newick": "(Human,Mouse,Opossum)root",
            "edge_attributes": {
                "Human": {"length": None},
                "Mouse": {"length": None},
                "Opossum": {"length": None},
                "root": {"length": None},
            },
            "type": "cogent3.core.tree.PhyloNode",
            "version": "2025.5.8a9",
        },
        "alignment": {
            "seqs": {
                "Human": {
                    "name": "Human",
                    "seq": {
                        "type": "cogent3.core.sequence.SeqView",
                        "version": "2025.5.8a9",
                        "init_args": {
                            "seq": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
                            "seqid": "Human",
                            "step": 1,
                        },
                    },
                    "moltype": "dna",
                    "info": None,
                    "type": "cogent3.core.sequence.DnaSequence",
                    "version": "2025.5.8a9",
                    "annotation_offset": 0,
                },
                "Mouse": {
                    "name": "Mouse",
                    "seq": {
                        "type": "cogent3.core.sequence.SeqView",
                        "version": "2025.5.8a9",
                        "init_args": {
                            "seq": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
                            "seqid": "Mouse",
                            "step": 1,
                        },
                    },
                    "moltype": "dna",
                    "info": None,
                    "type": "cogent3.core.sequence.DnaSequence",
                    "version": "2025.5.8a9",
                    "annotation_offset": 0,
                },
                "Opossum": {
                    "name": "Opossum",
                    "seq": {
                        "type": "cogent3.core.sequence.SeqView",
                        "version": "2025.5.8a9",
                        "init_args": {
                            "seq": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
                            "seqid": "Opossum",
                            "step": 1,
                        },
                    },
                    "moltype": "dna",
                    "info": None,
                    "type": "cogent3.core.sequence.DnaSequence",
                    "version": "2025.5.8a9",
                    "annotation_offset": 0,
                },
            },
            "moltype": "dna",
            "info": {"source": "unknown"},
            "type": "cogent3.core.alignment.ArrayAlignment",
            "version": "2025.5.8a9",
        },
        "likelihood_construction": {
            "bins": 1,
            "loci": 1,
            "optimise_motif_probs": True,
            "motif_probs_from_align": True,
            "default_length": 1.0,
        },
        "param_rules": [
            {
                "par_name": "psubs",
                "edge": "Human",
                "init": {
                    "T": {
                        "T": 0.9971099722856885,
                        "C": 0.002886010264321592,
                        "A": 2.008724995022132e-06,
                        "G": 2.008724995022132e-06,
                    },
                    "C": {"T": 0.125, "C": 0.625, "A": 0.125, "G": 0.125},
                    "A": {"T": 0.125, "C": 0.125, "A": 0.625, "G": 0.12499999999999997},
                    "G": {"T": 0.125, "C": 0.125, "A": 0.125, "G": 0.625},
                },
                "lower": None,
                "upper": None,
            },
            {
                "par_name": "psubs",
                "edge": "Mouse",
                "init": {
                    "T": {"T": 0.625, "C": 0.12499999999999997, "A": 0.125, "G": 0.125},
                    "C": {"T": 0.125, "C": 0.625, "A": 0.125, "G": 0.125},
                    "A": {"T": 0.125, "C": 0.125, "A": 0.625, "G": 0.12499999999999997},
                    "G": {"T": 0.125, "C": 0.125, "A": 0.125, "G": 0.625},
                },
                "lower": None,
                "upper": None,
            },
            {
                "par_name": "psubs",
                "edge": "Opossum",
                "init": {
                    "T": {"T": 0.625, "C": 0.12499999999999997, "A": 0.125, "G": 0.125},
                    "C": {"T": 0.125, "C": 0.625, "A": 0.125, "G": 0.125},
                    "A": {"T": 0.125, "C": 0.125, "A": 0.625, "G": 0.12499999999999997},
                    "G": {"T": 0.125, "C": 0.125, "A": 0.125, "G": 0.625},
                },
                "lower": None,
                "upper": None,
            },
            {
                "par_name": "mprobs",
                "init": {
                    "T": 0.11111111111111109,
                    "C": 0.27777777777777773,
                    "A": 0.15555555555555559,
                    "G": 0.45555555555555555,
                },
                "lower": None,
                "upper": None,
            },
        ],
        "lnL": -98.78628627311876,
        "nfp": 39,
        "motif_probs": {
            "T": 0.11111111111111109,
            "C": 0.27777777777777773,
            "A": 0.15555555555555559,
            "G": 0.45555555555555555,
        },
        "DLC": True,
        "unique_Q": None,
        "type": "cogent3.evolve.parameter_controller.AlignmentLikelihoodFunction",
        "name": "BH",
        "version": "2025.5.8a9",
    }
    lf = deserialise_likelihood_function(old_json)
    assert pytest.approx(lf.lnL) == -98.78628627311876
