import xml.dom.minidom
from unittest import TestCase

from cogent3.parse.blast_xml import (
    HIT_XML_FIELDNAMES,
    HSP_XML_FIELDNAMES,
    BlastXMLResult,
    get_tag,
    minimal_blast_parser_7,
    parse_header,
    parse_hit,
    parse_hsp,
    parse_parameters,
)


class GetTagTests(TestCase):
    """Tests for the auxiliary function evaluating the tag objects."""

    def setUp(self):
        self.single_tag = xml.dom.minidom.parseString(
            "<outer>bla <inner>content</inner>bla</outer>",
        )
        self.double_tag = xml.dom.minidom.parseString(
            "<outer><inner>first content</inner><inner>second content</inner></outer>",
        )
        self.empty_tag = xml.dom.minidom.parseString("<outer></outer>")

    def test_get_tag_works(self):
        assert get_tag(self.single_tag, "inner") == "content"
        assert get_tag(self.double_tag, "inner") == "first content"
        assert get_tag(self.empty_tag, "inner") is None
        assert get_tag(self.empty_tag, "inner", "blue elephant") == "blue elephant"
        assert get_tag(self.single_tag, "non-existing tag") is None
        assert (
            get_tag(self.single_tag, "non-existing tag", "pink elephant")
            == "pink elephant"
        )
        assert get_tag(self.single_tag, "inner") == "content"

    def test_get_tag_fail(self):
        """Make sure the tag and name parameters are in the proper types."""
        self.assertRaises(AttributeError, get_tag, None, "h1")
        self.assertRaises(
            AttributeError,
            get_tag,
            "<h1>This is not a XML tag object</h1>",
            "h1",
        )


class MinimalBlastParser7Tests(TestCase):
    """Tests for the functions required by the Blast XML parsers."""

    def setUp(self):
        self.hit1 = xml.dom.minidom.parseString(HIT_WITH_ONE_HSP)
        self.hit2 = xml.dom.minidom.parseString(HIT_WITH_TWO_HSPS)
        self.hsp1 = xml.dom.minidom.parseString(HSP_ONE)
        self.hsp2 = xml.dom.minidom.parseString(HSP_TWO)
        self.hsp_gaps = xml.dom.minidom.parseString(HSP_WITH_GAPS)

        self.param = xml.dom.minidom.parseString(PARAM_XML)
        self.header = xml.dom.minidom.parseString(HEADER_XML)
        self.complete = xml.dom.minidom.parseString(HEADER_COMPLETE)

    def test_parse_header(self):
        """Fields from XML header tag should be available as dict."""
        data = parse_header(self.header)
        assert data.get("application") == "my Grandma"
        assert data.get("version") == "has"
        assert data.get("reference") == "furry"
        assert data.get("query_letters") == 27
        assert data.get("database") == "Cats"

    def test_parse_parameters(self):
        """Fields from XML parameter tag should be available as dict."""
        data = parse_parameters(self.param)
        assert data.get("matrix") == "BLOSUM62"
        assert data.get("expect") == "10"
        assert data.get("gap_open_penalty") == 11.1
        assert data.get("gap_extend_penalty") == 22.2
        assert data.get("filter") == "F"

    def test_parse_header_complete(self):
        """Fields from header+param tag should be available as dict."""
        # try to process header with parameters etc in the XML
        data = parse_header(self.complete)
        assert data.get("database") == "Cats"
        assert data.get("matrix") == "BLOSUM62"

    def test_parse_hit(self):
        """Should return a list with all values for a hit+hsp."""
        data = parse_hit(self.hit1)
        assert len(data) == 1
        d = dict(list(zip(HIT_XML_FIELDNAMES, data[0], strict=False)))
        assert d["SUBJECT_ID"] == "gi|148670104|gb|EDL02051.1|"
        assert (
            d["HIT_DEF"]
            == "insulin-like growth factor 2 receptor, isoform CRA_c [Mus musculus]"
        )
        assert d["HIT_ACCESSION"] == "2001"
        assert int(d["HIT_LENGTH"]) == 707
        # check hit with more HSPs
        data = parse_hit(self.hit2)
        assert len(data) == 2
        assert data[0] != data[1]

    def test_parse_hsp(self):
        """Should return list with all values for a hsp."""
        data = parse_hsp(self.hsp1)
        d = dict(list(zip(HSP_XML_FIELDNAMES, data, strict=False)))
        assert float(d["BIT_SCORE"]) == 1023.46
        assert float(d["SCORE"]) == 2645
        assert float(d["E_VALUE"]) == 0.333
        assert int(d["QUERY_START"]) == 4
        assert int(d["QUERY_END"]) == 18
        assert int(d["SUBJECT_START"]) == 5
        assert int(d["SUBJECT_END"]) == 19
        assert int(d["GAP_OPENINGS"]) == 0
        assert int(d["ALIGNMENT_LENGTH"]) == 14

        assert d["QUERY_ALIGN"] == "ELEPHANTTHISISAHITTIGER"
        assert d["MIDLINE_ALIGN"] == "ORCA-WHALE"
        assert d["SUBJECT_ALIGN"] == "SEALSTHIS---HIT--GER"


class BlastXmlResultTests(TestCase):
    """Tests parsing of output of Blast with output mode 7 (XML)."""

    def setUp(self):
        self.result = BlastXMLResult(COMPLETE_XML, xml=True)

    def test_options(self):
        """Constructor should take parser as an option."""
        result = BlastXMLResult(COMPLETE_XML, parser=minimal_blast_parser_7)
        assert len(list(result.keys())) == 1
        # make sure whether normal Blast parser still works upon code merge!

    def test_parsed_query_sequence(self):
        """The result dict should have one query sequence as a key."""
        # The full query sequence is not given in the XML file.
        # Thus it is not checked explicitly, only whether there is
        # exactly one found.
        assert len(list(self.result.keys())) == 1

    def test_parsed_iterations(self):
        """The result should have the right number of iterations."""
        n_iter = 0
        for _query_id, _hits in self.result.iter_hits_by_query():
            n_iter += 1
        assert n_iter == 1

    def test_parsed_hsps(self):
        """The result should have the right number of hsps."""
        n_hsps = 0
        for _query_id, hsps in self.result.iter_hits_by_query():
            n_hsps += len(hsps)
        assert n_hsps == 3

    def test_parse_hit_details(self):
        """The result should have data from hit fields."""
        for query in self.result:
            first_hsp = self.result[query][0][0]
            assert first_hsp["SUBJECT_ID"] == "gi|148670104|gb|EDL02051.1|"
            assert (
                first_hsp["HIT_DEF"]
                == "insulin-like growth factor 2 receptor, isoform CRA_c [Mus musculus]"
            )
            assert first_hsp["HIT_ACCESSION"] == "2001"
            assert first_hsp["HIT_LENGTH"] == 707

    def test_parse_hsp_details(self):
        """The result should have data from hsp fields."""
        for query in self.result:
            # should check integers in next version.
            first_hsp = self.result[query][0][0]
            assert first_hsp["QUERY ID"] == 1
            assert first_hsp["BIT_SCORE"] == "1023.46"
            assert first_hsp["SCORE"] == "2645"
            assert first_hsp["E_VALUE"] == "0.333"
            assert first_hsp["QUERY_START"] == "4"
            assert first_hsp["QUERY_END"] == "18"
            assert first_hsp["QUERY_ALIGN"] == "ELEPHANTTHISISAHITTIGER"
            assert first_hsp["MIDLINE_ALIGN"] == "ORCA-WHALE"
            assert first_hsp["SUBJECT_ALIGN"] == "SEALSTHIS---HIT--GER"
            assert first_hsp["SUBJECT_START"] == "5"
            assert first_hsp["SUBJECT_END"] == "19"
            assert first_hsp["PERCENT_IDENTITY"] == "55"
            assert first_hsp["POSITIVE"] == "555"
            assert first_hsp["GAP_OPENINGS"] == 0
            assert first_hsp["ALIGNMENT_LENGTH"] == "14"

            gap_hsp = self.result[query][0][1]
            assert gap_hsp["GAP_OPENINGS"] == "33"

    def test_best_hits_by_query(self):
        """Exercising best hits"""
        q, best_hits = next(self.result.best_hits_by_query())
        best_hit = best_hits[0]
        assert best_hit["QUERY ID"] == 1
        assert best_hit["BIT_SCORE"] == "1023.46"
        assert best_hit["SCORE"] == "2645"
        assert best_hit["E_VALUE"] == "0.333"
        assert best_hit["QUERY_START"] == "4"
        assert best_hit["QUERY_END"] == "18"
        assert best_hit["QUERY_ALIGN"] == "ELEPHANTTHISISAHITTIGER"
        assert best_hit["MIDLINE_ALIGN"] == "ORCA-WHALE"
        assert best_hit["SUBJECT_ALIGN"] == "SEALSTHIS---HIT--GER"
        assert best_hit["SUBJECT_START"] == "5"
        assert best_hit["SUBJECT_END"] == "19"
        assert best_hit["PERCENT_IDENTITY"] == "55"
        assert best_hit["POSITIVE"] == "555"
        assert best_hit["GAP_OPENINGS"] == 0
        assert best_hit["ALIGNMENT_LENGTH"] == "14"

    def test_best_hits_unique(self):
        """The result should never contain identical hits"""
        records = next(h for _, h in self.result.best_hits_by_query(n=5))
        assert len(records) == 3
        values = {tuple(h.values()) for h in records}
        assert len(values) == 3


HSP_XML = """
        <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>1023.46</Hsp_bit-score>
              <Hsp_score>2645</Hsp_score>
              <Hsp_evalue>0.333</Hsp_evalue>
              <Hsp_query-from>4</Hsp_query-from>
              <Hsp_query-to>18</Hsp_query-to>
              <Hsp_hit-from>5</Hsp_hit-from>
              <Hsp_hit-to>19</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>55</Hsp_identity>
              %s
              <Hsp_positive>555</Hsp_positive>
              <Hsp_align-len>14</Hsp_align-len>
<Hsp_qseq>ELEPHANTTHISISAHITTIGER</Hsp_qseq>
              <Hsp_hseq>SEALSTHIS---HIT--GER</Hsp_hseq>
              <Hsp_midline>ORCA-WHALE</Hsp_midline>
            </Hsp>
            """
HSP_ONE = HSP_XML % ""
HSP_WITH_GAPS = HSP_XML % "<Hsp_gaps>33</Hsp_gaps>"

HSP_TWO = """
        <Hsp>
              <Hsp_num>2</Hsp_num>
              <Hsp_bit-score>1023.46</Hsp_bit-score>
              <Hsp_score>2645</Hsp_score>
              <Hsp_evalue>0.333</Hsp_evalue>
              <Hsp_query-from>6</Hsp_query-from>
              <Hsp_query-to>22</Hsp_query-to>
              <Hsp_hit-from>5</Hsp_hit-from>
              <Hsp_hit-to>23</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>55</Hsp_identity>
              %s
              <Hsp_positive>555</Hsp_positive>
              <Hsp_align-len>18</Hsp_align-len>
<Hsp_qseq>EPHANT---THISISAHIT-TIGER</Hsp_qseq>
              <Hsp_hseq>ALSWWWTHIS---HITW--GER</Hsp_hseq>
              <Hsp_midline>ORCA-WHALE</Hsp_midline>
            </Hsp>
            """
HIT_XML = """
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gi|148670104|gb|EDL02051.1|</Hit_id>
          <Hit_def>insulin-like growth factor 2 receptor, isoform CRA_c [Mus musculus]</Hit_def>
          <Hit_accession>2001</Hit_accession>
          <Hit_len>707</Hit_len>
          <Hit_hsps>
          %s
          </Hit_hsps>
        </Hit>
"""

HIT_WITH_ONE_HSP = HIT_XML % HSP_ONE
HIT_WITH_TWO_HSPS = HIT_XML % (HSP_WITH_GAPS + HSP_TWO)

PARAM_XML = """
<BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_gap-open>11.1</Parameters_gap-open>
      <Parameters_gap-extend>22.2</Parameters_gap-extend>
      <Parameters_filter>F</Parameters_filter>
    </Parameters>
</BlastOutput_param>
"""

HEADER_XML = """
<BlastOutput>
  <BlastOutput_program>my Grandma</BlastOutput_program>
  <BlastOutput_version>has</BlastOutput_version>
  <BlastOutput_db>Cats</BlastOutput_db>
  <BlastOutput_reference>furry</BlastOutput_reference>
  <BlastOutput_query-len>27</BlastOutput_query-len>

  %s
</BlastOutput>
"""

HIT_PREFIX = """
 <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
"""

HIT_SUFFIX = """
      </Iteration_hits>
    </Iteration>
 </BlastOutput_iterations>
"""

HEADER_COMPLETE = HEADER_XML % (
    PARAM_XML + HIT_PREFIX + HIT_WITH_ONE_HSP + HIT_WITH_TWO_HSPS + HIT_SUFFIX
)
COMPLETE_XML = (
    """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
"""
    + HEADER_COMPLETE
)
