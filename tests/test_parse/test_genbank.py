#!/usr/bin/env python
"""Unit tests for the GenBank database parsers.
"""
from unittest import TestCase, main

from cogent3.parse.genbank import (
    Location,
    LocationList,
    RichGenbankParser,
    block_consolidator,
    indent_splitter,
    location_line_tokenizer,
    parse_feature,
    parse_location_line,
    parse_locus,
    parse_organism,
    parse_reference,
    parse_sequence,
    parse_simple_location_segment,
    parse_single_line,
    parse_source,
)


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class GenBankTests(TestCase):
    """Tests of the GenBank main functions."""

    def test_parse_locus(self):
        """parse_locus should give correct results on specimen locus lines"""
        line = "LOCUS       AF108830                5313 bp    mRNA    linear   PRI 19-MAY-1999"
        result = parse_locus(line)
        self.assertEqual(len(result), 6)
        self.assertEqual(result["locus"], "AF108830")
        self.assertEqual(result["length"], 5313)  # note: int, not str
        self.assertEqual(result["mol_type"], "mRNA")
        self.assertEqual(result["topology"], "linear")
        self.assertEqual(result["db"], "PRI")
        self.assertEqual(result["date"], "19-MAY-1999")
        # should work if some of the fields are missing
        line = "LOCUS       AF108830                5313"
        result = parse_locus(line)
        self.assertEqual(len(result), 2)
        self.assertEqual(result["locus"], "AF108830")
        self.assertEqual(result["length"], 5313)  # note: int, not str

    def test_parse_single_line(self):
        """parse_single_line should split off the label and return the rest"""
        line_1 = "VERSION     AF108830.1  GI:4868112\n"
        self.assertEqual(parse_single_line(line_1), "AF108830.1  GI:4868112")
        # should work if leading spaces
        line_2 = "      VERSION     AF108830.1  GI:4868112\n"
        self.assertEqual(parse_single_line(line_2), "AF108830.1  GI:4868112")

    def test_indent_splitter(self):
        """indent_splitter should split lines at correct locations"""
        # if lines have same indent, should not group together
        lines = ["abc    xxx", "def    yyy"]
        self.assertEqual(list(indent_splitter(lines)), [[lines[0]], [lines[1]]])
        # if second line is indented, should group with first
        lines = ["abc    xxx", " def    yyy"]
        self.assertEqual(list(indent_splitter(lines)), [[lines[0], lines[1]]])

        # if both lines indented but second is more, should group with first
        lines = [" abc    xxx", "  def    yyy"]
        self.assertEqual(list(indent_splitter(lines)), [[lines[0], lines[1]]])

        # if both lines indented equally, should not group
        lines = ["   abc    xxx", "   def    yyy"]
        self.assertEqual(list(indent_splitter(lines)), [[lines[0]], [lines[1]]])

        # for more complex situation, should produce correct grouping
        lines = [
            "  xyz",  # 0 -
            "  xxx",  # 1 -
            "   yyy",  # 2
            "   uuu",  # 3
            "   iii",  # 4
            "  qaz",  # 5 -
            "  wsx",  # 6 -
            "   az",  # 7
            "   sx",  # 8
            "        gb",  # 9
            "   bg",  # 10
            "  aaa",  # 11 -
        ]
        self.assertEqual(
            list(indent_splitter(lines)),
            [[lines[0]], lines[1:5], [lines[5]], lines[6:11], [lines[11]]],
        )

        # real example from genbank file
        lines = """LOCUS       NT_016354           92123751 bp    DNA     linear   CON 29-AUG-2006
DEFINITION  Homo sapiens chromosome 4 genomic contig, reference assembly.
ACCESSION   NT_016354 NT_006109 NT_006204 NT_006245 NT_006302 NT_006371
            NT_006397 NT_016393 NT_016589 NT_016599 NT_016606 NT_022752
            NT_022753 NT_022755 NT_022760 NT_022774 NT_022797 NT_022803
            NT_022846 NT_022960 NT_025694 NT_028147 NT_029273 NT_030643
            NT_030646 NT_030662 NT_031780 NT_031781 NT_031791 NT_034703
            NT_034705 NT_037628 NT_037629 NT_079512
VERSION     NT_016354.18  GI:88977422
KEYWORDS    .
SOURCE      Homo sapiens (human)
  ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;
            Catarrhini; Hominidae; Homo.
?
REFERENCE   2  (bases 1 to 92123751)
  AUTHORS   International Human Genome Sequencing Consortium.
  TITLE     Finishing the euchromatic sequence of the human genome""".split(
            "\n"
        )
        self.assertEqual(
            list(indent_splitter(lines)),
            [
                [lines[0]],
                [lines[1]],
                lines[2:8],
                [lines[8]],
                [lines[9]],
                lines[10:15],
                [lines[15]],
                lines[16:],
            ],
        )

    def test_parse_sequence(self):
        """parse_sequence should strip bad chars out of sequence lines"""
        lines = """
ORIGIN
        1 gggagcgcgg cgcgggagcc cgaggctgag actcaccgga ggaagcggcg cgagcgcccc
       61   gccatcgtcc \t\t cggctgaagt 123 \ngcagtg  \n
      121 cctgggctta agcagtcttc45ccacctcagc 
//\n\n\n""".split(
            "\n"
        )
        result = parse_sequence(lines)
        self.assertEqual(
            result,
            "gggagcgcggcgcgggagcccgaggctgagactcaccggaggaagcggcgcgagcgccccgccatcgtcccggctgaagtgcagtgcctgggcttaagcagtcttcccacctcagc",
        )

    def test_block_consolidator(self):
        """block_consolidator should join the block together."""
        lines = """  ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Euarchontoglires; Primates; Catarrhini;
            Hominidae; Homo.""".split(
            "\n"
        )
        label, data = block_consolidator(lines)
        self.assertEqual(label, "ORGANISM")
        self.assertEqual(
            data,
            [
                "Homo sapiens",
                "            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;",
                "            Mammalia; Eutheria; Euarchontoglires; Primates; Catarrhini;",
                "            Hominidae; Homo.",
            ],
        )
        lines = r"""COMMENT
                    Contact: Spindel ER
                    Division of Neuroscience""".splitlines()
        label, data = block_consolidator(lines)
        self.assertEqual(label, "COMMENT")
        self.assertEqual(
            data,
            [
                "",
                "                    Contact: Spindel ER",
                "                    Division of Neuroscience",
            ],
        )

    def test_parse_organism(self):
        """parse_organism should return species, taxonomy (up to genus)"""
        # note: lines modified to include the following:
        # - multiword names
        # - multiword names split over a line break
        # - periods and other punctuation in names
        lines = """  ORGANISM  Homo sapiens
        Eukaryota; Metazoa; Chordata Craniata; Vertebrata; Euteleostomi;
        Mammalia; Eutheria; Euarchontoglires; Primates \t abc.  2.; Catarrhini
        Hominidae; Homo.""".split(
            "\n"
        )
        species, taxonomy = parse_organism(lines)
        self.assertEqual(species, "Homo sapiens")
        self.assertEqual(
            taxonomy,
            [
                "Eukaryota",
                "Metazoa",
                "Chordata Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Euarchontoglires",
                "Primates abc. 2.",
                "Catarrhini Hominidae",
                "Homo",
            ],
        )

    def test_parse_feature(self):
        """parse_feature should return dict containing annotations of feature"""
        example_feature = """     CDS             complement(join(102262..102647,105026..105217,
                     106638..106719,152424..152682,243209..243267))
                     /gene="nad1"
                     /note="Protein sequence is in conflict with the conceptual
                     translation; author given translation (not conceptual
                     translation)
                     start codon is created by C to U RNA editing"
                     /codon_start=1
                     /exception="RNA editing"
                     /product="NADH dehydrogenase subunit 1"
                     /protein_id="NP_064011.1"
                     /db_xref="GI:9838451"
                     /db_xref="IPI:12345"
                     /translation="MYIAVPAEILGIILPLLLGVAFLVLAERKVMAFVQRRKGPDVVG
                     SFGLLQPLADGSKLILKEPISPSSANFSLFRMAPVTTFMLSLVARAVVPFDYGMVLSD
                     PNIGLLYLFAISSLGVYGIIIAGWSSNSKYAFLGALRSAAQMVPYEVSIGLILITVLI
                     CVGPRNSSEIVMAQKQIWSGIPLFPVLVMFFISCLAETNRAPFDLPEAERELVAGYNV
                     EYSSMGSALFFLGEYANMILMSGLCTSLSPGGWPPILDLPISKRIPGSIWFSIKVILF
                     LFLYIWVRAAFPRYRYDQLMGLGRKVFLPLSLARVVAVSGVLVTFQWLP"""
        result = parse_feature(example_feature.split("\n"))
        self.assertEqual(result["type"], "CDS")
        self.assertEqual(
            result["raw_location"],
            [
                "complement(join(102262..102647,105026..105217,",
                "                     106638..106719,152424..152682,243209..243267))",
            ],
        )
        self.assertEqual(result["gene"], ["nad1"])
        self.assertEqual(
            result["note"],
            [
                "Protein sequence is in conflict with the conceptual translation; author given translation (not conceptual translation) start codon is created by C to U RNA editing"
            ],
        )
        self.assertEqual(result["codon_start"], ["1"])
        self.assertEqual(result["exception"], ["RNA editing"])
        self.assertEqual(result["product"], ["NADH dehydrogenase subunit 1"])
        self.assertEqual(result["protein_id"], ["NP_064011.1"])
        self.assertEqual(result["db_xref"], ["GI:9838451", "IPI:12345"])
        self.assertEqual(
            result["translation"],
            [
                "MYIAVPAEILGIILPLLLGVAFLVLAERKVMAFVQRRKGPDVVGSFGLLQPLADGSKLILKEPISPSSANFSLFRMAPVTTFMLSLVARAVVPFDYGMVLSDPNIGLLYLFAISSLGVYGIIIAGWSSNSKYAFLGALRSAAQMVPYEVSIGLILITVLICVGPRNSSEIVMAQKQIWSGIPLFPVLVMFFISCLAETNRAPFDLPEAERELVAGYNVEYSSMGSALFFLGEYANMILMSGLCTSLSPGGWPPILDLPISKRIPGSIWFSIKVILFLFLYIWVRAAFPRYRYDQLMGLGRKVFLPLSLARVVAVSGVLVTFQWLP"
            ],
        )
        self.assertEqual(len(result), 11)

        short_feature = ["D-loop          15418..16866"]
        result = parse_feature(short_feature)
        self.assertEqual(result["type"], "D-loop")
        self.assertEqual(result["raw_location"], ["15418..16866"])
        # can get more than one = in a line
        # from AF260826
        bad_feature = """     tRNA            1173..1238
                     /note="codon recognized: AUC; Cove score = 16.56"
                     /product="tRNA-Ile"
                     /anticodon=(pos:1203..1205,aa:Ile)"""
        result = parse_feature(bad_feature.split("\n"))
        self.assertEqual(result["note"], ["codon recognized: AUC; Cove score = 16.56"])
        # need not always have an = in a line
        # from NC_001807
        bad_feature = '''     mRNA            556
     /partial
     /citation=[6]
     /product="H-strand"'''
        result = parse_feature(bad_feature.split("\n"))
        self.assertEqual(result["partial"], [""])

    def test_location_line_tokenizer(self):
        """location_line_tokenizer should tokenize location lines"""
        llt = location_line_tokenizer
        self.assertEqual(list(llt(["123..456"])), ["123..456"])
        self.assertEqual(
            list(llt(["complement(123..456)"])), ["complement(", "123..456", ")"]
        )
        self.assertEqual(
            list(llt(["join(1..2,3..4)"])), ["join(", "1..2", ",", "3..4", ")"]
        )
        self.assertEqual(
            list(
                llt(
                    [
                        "join(complement(1..2, join(complement( 3..4),",
                        "\n5..6), 7..8\t))",
                    ]
                )
            ),
            [
                "join(",
                "complement(",
                "1..2",
                ",",
                "join(",
                "complement(",
                "3..4",
                ")",
                ",",
                "5..6",
                ")",
                ",",
                "7..8",
                ")",
                ")",
            ],
        )

    def test_parse_simple_location_segment(self):
        """parse_simple_location_segment should parse simple segments"""
        lsp = parse_simple_location_segment
        l = lsp("37")
        self.assertEqual(l._data, 37)
        self.assertEqual(str(l), "37")
        self.assertEqual(l.Strand, 1)
        l = lsp("40..50")
        first, second = l._data
        self.assertEqual(first._data, 40)
        self.assertEqual(second._data, 50)
        self.assertEqual(str(l), "40..50")
        self.assertEqual(l.Strand, 1)
        # should handle ambiguous starts and ends
        l = lsp(">37")
        self.assertEqual(l._data, 37)
        self.assertEqual(str(l), ">37")
        l = lsp("<37")
        self.assertEqual(l._data, 37)
        self.assertEqual(str(l), "<37")
        l = lsp("<37..>42")
        first, second = l._data
        self.assertEqual(first._data, 37)
        self.assertEqual(second._data, 42)
        self.assertEqual(str(first), "<37")
        self.assertEqual(str(second), ">42")
        self.assertEqual(str(l), "<37..>42")

    def test_parse_location_line(self):
        """parse_location_line should give correct list of location objects"""
        llt = location_line_tokenizer
        r = parse_location_line(llt(["123..456"]))
        self.assertEqual(str(r), "123..456")
        r = parse_location_line(llt(["complement(123..456)"]))
        self.assertEqual(str(r), "complement(123..456)")
        r = parse_location_line(llt(["complement(123..456, 345..678)"]))
        self.assertEqual(str(r), "join(complement(345..678),complement(123..456))")
        r = parse_location_line(llt(["complement(join(123..456, 345..678))"]))
        self.assertEqual(str(r), "join(complement(345..678),complement(123..456))")
        r = parse_location_line(
            llt(["join(complement(123..456), complement(345..678))"])
        )
        self.assertEqual(str(r), "join(complement(123..456),complement(345..678))")
        # try some nested joins and complements
        r = parse_location_line(
            llt(
                [
                    "complement(join(1..2,3..4,complement(5..6),",
                    "join(7..8,complement(9..10))))",
                ]
            )
        )
        self.assertEqual(
            str(r),
            "join(9..10,complement(7..8),5..6,complement(3..4),complement(1..2))",
        )

    def test_parse_reference(self):
        """parse_reference should give correct fields"""
        r = """REFERENCE   2  (bases 1 to 2587)
  AUTHORS   Janzen,D.M. and Geballe,A.P.
  TITLE     The effect of eukaryotic release factor depletion on translation
            termination in human cell lines
  JOURNAL   (er) Nucleic Acids Res. 32 (15), 4491-4502 (2004)
   PUBMED   15326224"""
        result = parse_reference(r.split("\n"))
        self.assertEqual(len(result), 5)
        self.assertEqual(result["reference"], "2  (bases 1 to 2587)")
        self.assertEqual(result["authors"], "Janzen,D.M. and Geballe,A.P.")
        self.assertEqual(
            result["title"],
            "The effect of eukaryotic release factor depletion "
            + "on translation termination in human cell lines",
        )
        self.assertEqual(
            result["journal"], "(er) Nucleic Acids Res. 32 (15), 4491-4502 (2004)"
        )
        self.assertEqual(result["pubmed"], "15326224")

    def test_parse_source(self):
        """parse_source should split into source and organism"""
        s = """SOURCE      African elephant.
  ORGANISM  Mitochondrion Loxodonta africana
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Proboscidea; Elephantidae; Loxodonta.""".split(
            "\n"
        )
        r = parse_source(s)
        self.assertEqual(len(r), 3)
        self.assertEqual(r["source"], "African elephant.")
        self.assertEqual(r["species"], "Mitochondrion Loxodonta africana")
        self.assertEqual(
            r["taxonomy"],
            [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Proboscidea",
                "Elephantidae",
                "Loxodonta",
            ],
        )

    def test_rich_parser(self):
        """correctly constructs +/- strand features"""
        # a custom annotation function
        from cogent3.core.annotation import Feature

        def add_annotation(seq, feature, spans):
            if feature["type"] != "CDS":
                return
            name = feature["locus_tag"][0]
            seq.add_annotation(Feature, "CDS", name, spans)

        infile = open("data/annotated_seq.gb")
        parser = RichGenbankParser(infile, add_annotation=add_annotation)

        seq = [s for l, s in parser][0]
        cds = dict([(f.name, f) for f in seq.get_annotations_matching("CDS")])
        expects = {
            "CNA00110": "MAGYDARYGNPLDPMSGGRPSPPETSQQDAYEYSKHGSSSGYLGQLPLGAD"
            "SAQAETASALRTLFGEGADVQALQEPPNQINTLAEGAAVAETGGVLGGDTTRSDNEALAIDPSL"
            "SEQAAPAPKDSTETPDDRSRSPSSGNHHHHHPAVKRKATSRAGMLARGGACEFCKRRKLKCSAEL"
            "PACANCVKSGKECVYAQKKQRSRVKVLEDRLQELEKRLEQGQAGAASASGGDSGAHAASSVYTAP"
            "SLGSGGGSELTVEQTLVHNVDPSLLPPSEYDEAFILHDFDSFADMRKQETQLEPDLMTLADAAAA"
            "DTPAAAAAETNDPWAKMSPEEIVKEIIKVATGGKGEGERIISHLVQTYMNSTVNTWHPLVIPPMD"
            "LVSRVSRTTPDPIHPTLLLSLIPALLPLSPIQSLRHPAIPLLLLPHARAHSVQAITQSDPRVLDT"
            "IIAGVSRAYSFFNEAKNIDGWVDCVAATSLVRAAGLTKQGGVGERFVPEDRVPAERLAKRRREAG"
            "LRALMHKGAIVPPPESWYQFGQRVNLFWTSYICDRAAAIGWGWPSSYNDEDITTPWPKDDYKSVQ"
            "ALLDDTTIHTFLSPLAPAPAPATPDSDLCAQAKSITLLYHAQRLLDSPPELSTPEKTHRLLGLTE"
            "GYMESLEKMRGPRMRAGKLSSVWMILYTTIAVLHSKDGFDKCDPDGADQVSITRVVAAADKVLEL"
            "VSAVQNTGDTHLSSCDVISSVLFLHLARLMIQYTNRLRLRVQDSALVSTLRAKTESFKRALIDQG"
            "ERLVFAQVAAQMLENYHVGAEWKAGEWERADGGDWRGV",
            "CNA00120": "MDFSQFNGAEQAHMSKVIEKKQMQDFMRLYSGLVEKCFNACAQD"
            "FTSKALTTNETTCVQNCTDKFLKHSERVGARFAEHNAGMLSPYGAASLMASQSKCRAP"
            "DSNGLGVFCKWRRIKSTVVLYNHLACIKQMDNRF",
        }

        for locus in cds:
            got = cds[locus].get_slice().trim_stop_codon().get_translation()
            self.assertEqual(str(got), expects[locus])
        infile.close()

    def test_rich_parser_moltype(self):
        """correctly handles moltypes"""
        with open("data/annotated_seq.gb") as infile:
            parser = RichGenbankParser(infile)
            got_1 = [s for _, s in parser][0]

        # name formed from /product value
        got = {f.name for f in got_1.get_annotations_matching("mRNA")}
        self.assertEqual(got, {"conserved hypothetical protein", "chaperone, putative"})

        # the file defines itself as DNA
        self.assertEqual(got_1.moltype.label, "dna")

        # but that is overridden by user setting moltype explicitly
        for moltype in ("dna", "rna", "text"):
            with open("data/annotated_seq.gb") as infile:
                parser = RichGenbankParser(infile, moltype=moltype)
                got_2 = [s for _, s in parser][0]

            self.assertEqual(len(got_1.annotations), len(got_2.annotations))
            self.assertEqual(got_2.moltype.label, moltype)
            got = {f.name for f in got_1.get_annotations_matching("mRNA")}
            self.assertEqual(
                got, {"conserved hypothetical protein", "chaperone, putative"}
            )


class LocationTests(TestCase):
    """Tests of the Location class."""

    def test_init(self):
        """Location should init with 1 or 2 values, plus params."""
        l = Location(37)
        self.assertEqual(str(l), "37")
        l = Location(37, Ambiguity=">")
        self.assertEqual(str(l), ">37")
        l = Location(37, Ambiguity="<")
        self.assertEqual(str(l), "<37")
        l = Location(37, Accession="AB123")
        self.assertEqual(str(l), "AB123:37")
        l = Location(37, Accession="AB123", Db="Kegg")
        self.assertEqual(str(l), "Kegg::AB123:37")

        l1 = Location(37)
        l2 = Location(42)
        l = Location([l1, l2])
        self.assertEqual(str(l), "37..42")
        l3 = Location([l1, l2], IsBounds=True)
        self.assertEqual(str(l3), "(37.42)")
        l4 = Location([l1, l2], IsBetween=True)
        self.assertEqual(str(l4), "37^42")
        l5 = Location([l4, l3])
        self.assertEqual(str(l5), "37^42..(37.42)")
        l5 = Location([l4, l3], Strand=-1)
        self.assertEqual(str(l5), "complement(37^42..(37.42))")


class LocationListTests(TestCase):
    """Tests of the LocationList class."""

    def test_extract(self):
        """LocationList extract should return correct sequence"""
        l = Location(3)
        l2_a = Location(5)
        l2_b = Location(7)
        l2 = Location([l2_a, l2_b], Strand=-1)
        l3_a = Location(10)
        l3_b = Location(12)
        l3 = Location([l3_a, l3_b])
        ll = LocationList([l, l2, l3])
        s = ll.extract("ACGTGCAGTCAGTAGCAT")
        #               123456789012345678
        self.assertEqual(s, "G" + "TGC" + "CAG")
        # check a case where it wraps around
        l5_a = Location(16)
        l5_b = Location(4)
        l5 = Location([l5_a, l5_b])
        ll = LocationList([l5])
        s = ll.extract("ACGTGCAGTCAGTAGCAT")
        self.assertEqual(s, "CATACGT")


if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2 and argv[1] == "x":
        filename = argv[2]
        lines = open(filename)
        for i in indent_splitter(lines):
            print("******")
            print(i[0])
            for j in indent_splitter(i[1:]):
                print("?????")
                for line in j:
                    print(line)
    else:
        main()
