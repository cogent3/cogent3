"""Tests for molecular weight."""

from numpy.testing import assert_allclose

from cogent3.data.molecular_weight import ProteinMW, RnaMW


def test_call():
    """WeightCalculator should return correct molecular weight"""
    r = RnaMW
    p = ProteinMW
    assert p("") == 0
    assert r("") == 0
    assert_allclose(p("A"), 89.09)
    assert_allclose(r("A"), 375.17)
    assert_allclose(p("AAA"), 231.27)
    assert_allclose(r("AAA"), 1001.59)
    assert_allclose(r("AAACCCA"), 2182.37)
    big_seq = (
        "MVQQAESLEAESNLPREALDTEEGEFMACSPVALDESDPDWCKTASGHIKRPMNAFMVWS"
        "KIERRKIMEQSPDMHNAEISKRLGKRWKMLKDSEKIPFIREAERLRLKHMADYPDYKYRP"
        "RKKPKMDPSAKPSASQSPEKSAAGGGGGSAGGGAGGAKTSKGSSKKCGKLKAPAAAGAKA"
        "GAGKAAQSGDYGGAGDDYVLGSLRVSGSGGGGAGKTVKCVFLDEDDDDDDDDDELQLQIK"
        "QEPDEEDEEPPHQQLLQPPGQQPSQLLRRYNVAKVPASPTLSSSAESPEGASLYDEVRAG"
        "ATSGAGGGSRLYYSFKNITKQHPPPLAQPALSPASSRSVSTSSSSSSGSSSGSSGEDADD"
        "LMFDLSLNFSQSAHSASEQQLGGGAAAGNLSLSLVDKDLDSFSEGSLGSHFEFPDYCTPE"
        "LSEMIAGDWLEANFSDLVFTY"
    )
    assert_allclose(p(big_seq), 46685.97)
