#! /usr/bin/env python
#
# test_blast_xml.py
#

__author__ = "Kristian Rother"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__contributors__ = ["Micah Hamady"]
__credits__ = ["Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Prototype"


import xml.dom.minidom

from unittest import TestCase, main

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
            "<outer>bla <inner>content</inner>bla</outer>"
        )
        self.double_tag = xml.dom.minidom.parseString(
            "<outer><inner>first content</inner><inner>second content</inner></outer>"
        )
        self.empty_tag = xml.dom.minidom.parseString("<outer></outer>")

    def test_get_tag_works(self):
        self.assertEqual(get_tag(self.single_tag, "inner"), "content")
        self.assertEqual(get_tag(self.double_tag, "inner"), "first content")
        self.assertEqual(get_tag(self.empty_tag, "inner"), None)
        self.assertEqual(
            get_tag(self.empty_tag, "inner", "blue elephant"), "blue elephant"
        )
        self.assertEqual(get_tag(self.single_tag, "non-existing tag"), None)
        self.assertEqual(
            get_tag(self.single_tag, "non-existing tag", "pink elephant"),
            "pink elephant",
        )
        self.assertEqual(get_tag(self.single_tag, "inner"), "content")

    def test_get_tag_fail(self):
        """Make sure the tag and name parameters are in the proper types."""
        self.assertRaises(AttributeError, get_tag, None, "h1")
        self.assertRaises(
            AttributeError, get_tag, "<h1>This is not a XML tag object</h1>", "h1"
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
        self.assertEqual(data.get("application"), "my Grandma")
        self.assertEqual(data.get("version"), "has")
        self.assertEqual(data.get("reference"), "furry")
        self.assertEqual(data.get("query_letters"), 27)
        self.assertEqual(data.get("database"), "Cats")

    def test_parse_parameters(self):
        """Fields from XML parameter tag should be available as dict."""
        data = parse_parameters(self.param)
        self.assertEqual(data.get("matrix"), "BLOSUM62")
        self.assertEqual(data.get("expect"), "10")
        self.assertEqual(data.get("gap_open_penalty"), 11.1)
        self.assertEqual(data.get("gap_extend_penalty"), 22.2)
        self.assertEqual(data.get("filter"), "F")

    def test_parse_header_complete(self):
        """Fields from header+param tag should be available as dict."""
        # try to process header with parameters etc in the XML
        data = parse_header(self.complete)
        self.assertEqual(data.get("database"), "Cats")
        self.assertEqual(data.get("matrix"), "BLOSUM62")

    def test_parse_hit(self):
        """Should return a list with all values for a hit+hsp."""
        data = parse_hit(self.hit1)
        self.assertEqual(len(data), 1)
        d = dict(list(zip(HIT_XML_FIELDNAMES, data[0])))
        self.assertEqual(d["SUBJECT_ID"], "gi|148670104|gb|EDL02051.1|")
        self.assertEqual(
            d["HIT_DEF"],
            "insulin-like growth factor 2 receptor, isoform CRA_c [Mus musculus]",
        )
        self.assertEqual(d["HIT_ACCESSION"], "2001")
        self.assertEqual(int(d["HIT_LENGTH"]), 707)
        # check hit with more HSPs
        data = parse_hit(self.hit2)
        self.assertEqual(len(data), 2)
        self.assertNotEqual(data[0], data[1])

    def test_parse_hsp(self):
        """Should return list with all values for a hsp."""
        data = parse_hsp(self.hsp1)
        d = dict(list(zip(HSP_XML_FIELDNAMES, data)))
        self.assertEqual(float(d["BIT_SCORE"]), 1023.46)
        self.assertEqual(float(d["SCORE"]), 2645)
        self.assertEqual(float(d["E_VALUE"]), 0.333)
        self.assertEqual(int(d["QUERY_START"]), 4)
        self.assertEqual(int(d["QUERY_END"]), 18)
        self.assertEqual(int(d["SUBJECT_START"]), 5)
        self.assertEqual(int(d["SUBJECT_END"]), 19)
        self.assertEqual(int(d["GAP_OPENINGS"]), 0)
        self.assertEqual(int(d["ALIGNMENT_LENGTH"]), 14)

        self.assertEqual(d["QUERY_ALIGN"], "ELEPHANTTHISISAHITTIGER")
        self.assertEqual(d["MIDLINE_ALIGN"], "ORCA-WHALE")
        self.assertEqual(d["SUBJECT_ALIGN"], "SEALSTHIS---HIT--GER")


class BlastXmlResultTests(TestCase):
    """Tests parsing of output of Blast with output mode 7 (XML)."""

    def setUp(self):
        self.result = BlastXMLResult(COMPLETE_XML, xml=True)

    def test_options(self):
        """Constructor should take parser as an option."""
        result = BlastXMLResult(COMPLETE_XML, parser=minimal_blast_parser_7)
        self.assertEqual(len(list(result.keys())), 1)
        # make sure whether normal Blast parser still works upon code merge!

    def test_parsed_query_sequence(self):
        """The result dict should have one query sequence as a key."""
        # The full query sequence is not given in the XML file.
        # Thus it is not checked explicitly, only whether there is
        # exactly one found.
        self.assertEqual(len(list(self.result.keys())), 1)

    def test_parsed_iterations(self):
        """The result should have the right number of iterations."""
        n_iter = 0
        for query_id, hits in self.result.iter_hits_by_query():
            n_iter += 1
        self.assertEqual(n_iter, 1)

    def test_parsed_hsps(self):
        """The result should have the right number of hsps."""
        n_hsps = 0
        for query_id, hsps in self.result.iter_hits_by_query():
            n_hsps += len(hsps)
        self.assertEqual(n_hsps, 3)

    def test_parse_hit_details(self):
        """The result should have data from hit fields."""
        for query in self.result:
            first_hsp = self.result[query][0][0]
            self.assertEqual(first_hsp["SUBJECT_ID"], "gi|148670104|gb|EDL02051.1|")
            self.assertEqual(
                first_hsp["HIT_DEF"],
                "insulin-like growth factor 2 receptor, isoform CRA_c [Mus musculus]",
            )
            self.assertEqual(first_hsp["HIT_ACCESSION"], "2001")
            self.assertEqual(first_hsp["HIT_LENGTH"], 707)

    def test_parse_hsp_details(self):
        """The result should have data from hsp fields."""
        for query in self.result:
            # should check integers in next version.
            first_hsp = self.result[query][0][0]
            self.assertEqual(first_hsp["QUERY ID"], 1)
            self.assertEqual(first_hsp["BIT_SCORE"], "1023.46")
            self.assertEqual(first_hsp["SCORE"], "2645")
            self.assertEqual(first_hsp["E_VALUE"], "0.333")
            self.assertEqual(first_hsp["QUERY_START"], "4")
            self.assertEqual(first_hsp["QUERY_END"], "18")
            self.assertEqual(first_hsp["QUERY_ALIGN"], "ELEPHANTTHISISAHITTIGER")
            self.assertEqual(first_hsp["MIDLINE_ALIGN"], "ORCA-WHALE")
            self.assertEqual(first_hsp["SUBJECT_ALIGN"], "SEALSTHIS---HIT--GER")
            self.assertEqual(first_hsp["SUBJECT_START"], "5")
            self.assertEqual(first_hsp["SUBJECT_END"], "19")
            self.assertEqual(first_hsp["PERCENT_IDENTITY"], "55")
            self.assertEqual(first_hsp["POSITIVE"], "555")
            self.assertEqual(first_hsp["GAP_OPENINGS"], 0)
            self.assertEqual(first_hsp["ALIGNMENT_LENGTH"], "14")

            gap_hsp = self.result[query][0][1]
            self.assertEqual(gap_hsp["GAP_OPENINGS"], "33")

    def test_best_hits_by_query(self):
        """Exercising best hits"""
        q, best_hits = next(self.result.best_hits_by_query())
        best_hit = best_hits[0]
        self.assertEqual(best_hit["QUERY ID"], 1)
        self.assertEqual(best_hit["BIT_SCORE"], "1023.46")
        self.assertEqual(best_hit["SCORE"], "2645")
        self.assertEqual(best_hit["E_VALUE"], "0.333")
        self.assertEqual(best_hit["QUERY_START"], "4")
        self.assertEqual(best_hit["QUERY_END"], "18")
        self.assertEqual(best_hit["QUERY_ALIGN"], "ELEPHANTTHISISAHITTIGER")
        self.assertEqual(best_hit["MIDLINE_ALIGN"], "ORCA-WHALE")
        self.assertEqual(best_hit["SUBJECT_ALIGN"], "SEALSTHIS---HIT--GER")
        self.assertEqual(best_hit["SUBJECT_START"], "5")
        self.assertEqual(best_hit["SUBJECT_END"], "19")
        self.assertEqual(best_hit["PERCENT_IDENTITY"], "55")
        self.assertEqual(best_hit["POSITIVE"], "555")
        self.assertEqual(best_hit["GAP_OPENINGS"], 0)
        self.assertEqual(best_hit["ALIGNMENT_LENGTH"], "14")

    def test_best_hits_unique(self):
        """The result should never contain identical hits"""
        records = [h for _, h in self.result.best_hits_by_query(n=5)][0]
        self.assertEqual(len(records), 3)
        values = {tuple(h.values()) for h in records}
        self.assertEqual(len(values), 3)


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


if __name__ == "__main__":
    main()
