#!/usr/bin/env python
""" Provides tests for EbiParser and related classes and functions.
"""

from unittest import TestCase, main

from cogent3.parse.ebi import (
    EbiFinder,
    EbiParser,
    FieldError,
    MinimalEbiParser,
    RecordError,
    ac_parser,
    cc_alternative_products_parser,
    cc_basic_itemparser,
    cc_biophysicochemical_properties_parser,
    cc_interaction_parser,
    cc_itemfinder,
    cc_parser,
    curry,
    de_itemparser,
    de_parser,
    dr_parser,
    dt_parser,
    ft_basic_itemparser,
    ft_id_parser,
    ft_mutagen_parser,
    ft_mutation_parser,
    ft_parser,
    gn_parser,
    hanging_paragraph_finder,
    id_parser,
    join_parser,
    join_split_dict_parser,
    join_split_parser,
    labeloff,
    linecode_maker,
    linecode_merging_maker,
    mapping_parser,
    oc_parser,
    og_parser,
    os_parser,
    ox_parser,
    pairs_to_dict,
    period_tail_finder,
    pr_parser,
    ra_parser,
    rc_parser,
    required_labels,
    rg_parser,
    rl_parser,
    rn_parser,
    rp_parser,
    rstrip_,
    rt_parser,
    rx_parser,
    single_ref_parser,
    sq_parser,
    try_int,
)


__author__ = "Zongzhi Liu"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Zongzhi Liu", "Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Zongzhi Liu"
__email__ = "zongzhi.liu@gmail.com"
__status__ = "Development"


def item_empty_filter(d):
    """return a dict with only nonempty values"""
    pairs = [(k, v) for (k, v) in d.items() if v]
    return dict(pairs)


class EbiTests(TestCase):
    """Tests ebi parsers and generic parsers and general functions"""

    def setUp(self):
        """Construct some fake data for testing purposes"""
        pass

    def test_item_empty_filter(self):
        """item_empty_filter: known values"""
        inputs = [{1: 0, 2: 1, 3: "", 4: [], 5: False, 6: {}}]
        expects = [{2: 1}]
        self.assertEqual(list(map(item_empty_filter, inputs)), expects)

    def test_rstrip_(self):
        """rstrip_ should generate the expected function"""
        test = " aaa;  "
        self.assertEqual(rstrip_("; ")(test), test.rstrip("; "))
        # test default
        self.assertEqual(rstrip_()(test), test.rstrip())

    def test_hanging_paragraph_finder(self):
        """hanging_paragraph_finder should give expected results"""
        f = hanging_paragraph_finder
        test = ["-aaa:", "  content", "    content", "-bbb:", "  bbb", "c"]  # 3
        self.assertEqual(list(f(test)), [test[0:3], test[3:-1], test[-1:]])

        # test all indent lines
        all_indent = ["  aa", "  bb"]
        self.assertEqual(list(f(all_indent)), [all_indent])

        # test empty lines
        self.assertEqual(list(f(["", " "])), [])

    def test_period_tail_finder(self):
        """period_tail_finder should yield each group of expected lines."""
        test = "a\naa.\nb\nbb.".splitlines()
        self.assertEqual(list(period_tail_finder(test)), [["a", "aa."], ["b", "bb."]])

    def test_EbiFinder(self):
        """EbiFinder should return expected list"""
        test = "a\n//\nb\n//".splitlines()
        self.assertEqual(list(EbiFinder(test)), [["a", "//"], ["b", "//"]])

        test_fail = test + ["c"]
        self.assertRaises(RecordError, list, EbiFinder(test_fail))

    def test_pairs_to_dict(self):
        """pairs_to_dict should return the expected dict"""
        test_dict = {"a": 1, "b": 2, "b": 3}
        sorted_items = [("a", 1), ("b", 2), ("b", 3)]
        add_one = lambda x: x + 1
        double = lambda x: x * 2
        set_zero = lambda x: 0
        handlers = {"a": add_one, "b": double}

        # test default all
        self.assertEqual(pairs_to_dict(sorted_items), {"a": 1, "b": 3})

        # test 'overwrite_value'
        self.assertEqual(
            pairs_to_dict(sorted_items, "overwrite_value"), {"a": 1, "b": 3}
        )

        # test no_duplicated_key, raise
        self.assertRaises(ValueError, pairs_to_dict, sorted_items, "no_duplicated_key")

        # test always_multi_value
        self.assertEqual(
            pairs_to_dict(sorted_items, "always_multi_value"), {"a": [1], "b": [2, 3]}
        )

        # test allow multi_value
        self.assertEqual(
            pairs_to_dict(sorted_items, "allow_multi_value"), {"a": 1, "b": [2, 3]}
        )

        # test raise error when  key not found in all_keys
        f = curry(pairs_to_dict, all_keys=["a", "c"])
        self.assertRaises(ValueError, f, sorted_items)

        # test handlers
        sorted_items.append(("c", 4))
        self.assertEqual(
            pairs_to_dict(sorted_items, handlers=handlers, default_handler=set_zero),
            {"a": 2, "b": 6, "c": 0},
        )

        # test raise error when  no valid handlers were found
        f = curry(pairs_to_dict, handlers=handlers)
        self.assertRaises(ValueError, f, sorted_items)

        # test sanity
        test_dict = dict(sorted_items)
        self.assertEqual(pairs_to_dict(list(test_dict.items())), test_dict)

    def test_linecode_maker(self):
        """linecode_maker: should return expected tuple"""
        tests = ["AA   aa.", "BB    bb.", "CC C   cc.", "DD dd."]
        expected_linecodes = ["AA", "BB", "CC C", "DD dd."]
        # pprint(map(linecode_maker, tests))
        self.assertEqual(
            list(map(linecode_maker, tests)), list(zip(expected_linecodes, tests))
        )

    def test_labeloff(self):
        """labeloff: should return expected lines"""
        tests = ["AA   aa.", "BB    bb.", "CC C   cc.", "DD dd.", "EE", ""]
        # expects = [line[5:] for line in tests]
        expects = ["aa.", " bb.", "  cc.", ".", "", ""]
        self.assertEqual(labeloff(tests), expects)

    def test_join_parser(self):
        """join parser should return expected str."""
        test_list = '"aaa\nbbb \nccc"; \n'.splitlines()
        test_str = "aaa bb. "
        # test default, list input
        self.assertEqual(join_parser(test_list), '"aaa bbb  ccc"')
        # test default, str input
        self.assertEqual(join_parser(test_str), "aaa bb")
        # test no strip
        self.assertEqual(join_parser(test_list, chars_to_strip=""), '"aaa bbb  ccc"; ')
        # test strip
        self.assertEqual(join_parser(test_list, chars_to_strip='"; '), "aaa bbb  ccc")
        # test empty
        self.assertEqual(join_parser([]), "")
        self.assertEqual(join_parser(["", " "]), "")
        self.assertEqual(join_parser(""), "")

    def test_join_split_parser(self):
        """join_split_parser: should return expected"""
        f = join_split_parser
        assertEqual = self.assertEqual

        assertEqual(f(["aa; bb;", "cc."]), ["aa", "bb", "cc"])
        assertEqual(
            f(["aa; bb, bbb;", "cc."], delimiters=";,"), ["aa", ["bb", "bbb"], "cc"]
        )
        # test item_modifer
        got = f("aa (bb) (cc).", "(", item_modifier=rstrip_(") "))
        assertEqual(got, ["aa", "bb", "cc"])
        assertEqual(
            f("aa (bb)xx (cc).", "(", item_modifier=rstrip_(") ")),
            ["aa", "bb)xx", "cc"],
        )

        # test empty
        assertEqual(f([]), [""])
        assertEqual(f(["", " "]), [""])
        assertEqual(f(""), [""])

    def test_join_split_dict_parser(self):
        """join_split_dict_parser: should return expected"""
        f = join_split_dict_parser
        # test default
        self.assertEqual(
            f("aa=1; bb=2,3; cc=4 (if aa=1);"),
            {"aa": "1", "bb": ["2", "3"], "cc": "4 (if aa=1)"},
        )

        self.assertEqual(
            f("aa=1; bb=2,3; cc=4:5", delimiters=";=,:"),
            {"aa": "1", "bb": ["2", "3"], "cc": "4:5"},
        )

        # test strict=False -> splits without dict()
        self.assertEqual(f("aa=1; bb.", strict=False), [["aa", "1"], ["bb"]])

        # test strict -> raise ValueError
        self.assertRaises(ValueError, f, "aa=1; bb.")
        self.assertRaises(ValueError, f, "aa=1; bb=2=3.", ";=")
        self.assertRaises(ValueError, f, "")

    def test_mapping_parser(self):
        """mapping_parser: should return expected dict"""
        fields = [None, "A", "B", ("C", int), ("D", float)]
        line = "blah aa bb; 2 3.1;"

        expect = dict(A="aa", B="bb", C=2, D=3.1)
        self.assertEqual(mapping_parser(line, fields), expect)

        # test more splits -> ignore the last splits
        line_leftover = line + "blah blah"
        self.assertEqual(mapping_parser(line_leftover, fields), expect)

        # test more fields -> truncate the last fields
        fields_leftover = fields + ["E"]
        self.assertEqual(mapping_parser(line, fields_leftover), expect)

        # test empty
        self.assertEqual(mapping_parser("", fields), {})

    def test_linecode_merging_maker(self):
        """linecode_merging_maker:"""
        f = linecode_merging_maker
        lines = ["ID   id.", "RN   rn.", "RR   invalid", "RN  rn."]
        labels = ["ID", "REF", "RR", "RN  rn."]

        self.assertEqual(list(map(f, lines)), list(zip(labels, lines)))

    def test_parse_header(lines):
        pass

    def test_MinimalEbiParser_valid(self):
        """MinimalEbiParser: integrity of output"""
        f = curry(MinimalEbiParser, strict=False)

        # test valid result: sequence, number of records, keys of a header
        valid_result = list(f(fake_records_valid))
        self.assertEqual(len(valid_result), 2)

        sequence, header = valid_result[0]
        self.assertEqual(sequence, "aaccppgghhh")
        # the first fake record use only the required labels, the header is
        # deleted of '', which was assigned to sequence
        self.assertEqual(
            list(sorted(header.keys())), list(sorted(required_labels))[1:]
        )  # [1:] to exclude the ''

        # test selected_labels
        selected_labels = ["ID", "DE"]
        select_result = list(f(fake_records_valid, selected_labels=selected_labels))
        self.assertEqual(
            list(sorted(select_result[0][1].keys())), list(sorted(selected_labels))
        )

        # test bad record - unknown linecode or wrong line format
        self.assertRaises(
            ValueError,
            list,
            f(fake_records_valid + ["ID   id.", "RR   not valid.", "//"]),
        )
        self.assertRaises(
            ValueError,
            list,
            f(fake_records_valid + ["ID   id.", " RN   bad format.", "//"]),
        )
        self.assertRaises(
            ValueError,
            list,
            f(fake_records_valid + ["ID   id.", "RN  bad format.", "//"]),
        )

        # test bad record - not end with '//'
        self.assertRaises(
            RecordError,
            list,
            f(fake_records_valid + ["ID   not end with //", "   seq"]),
        )

        # test strict: lacking required linecodes
        # ?? How to test the error message?  warn message?
        # the first record, [:-5], is valid even when strict=True.
        the_first_valid = list(f(fake_records_valid[:-5], strict=True))[0]
        # [1] get the header_dict
        self.assertEqual(len(the_first_valid[1]), 9)

        self.assertRaises(RecordError, list, f(fake_records_valid, strict=True))

    def test_EbiParser(self):
        """EbiParser:"""
        f = curry(EbiParser, strict=False)
        first_valid = fake_records_valid[:-5]

        # test valid
        self.assertEqual(len(list(f(fake_records_valid))), 2)
        # test skipping bad record which strict=False
        # self.assertEqual(len(list(f(fake_records_valid[:-1] +
        #    ['OX   xx=no equal.', '//']))), 1)
        # test Raise RecordError from parse_head when strict=True
        # self.assertRaises(RecordError, list, f(first_valid[:-1] +
        #    ['OX   xx=no equal.', '//'], strict=True))


class RootParsersKnownValues(TestCase):
    """Test most xx_parsers with known value"""

    def test_id_parser(self):
        """id_parser should return expected dict"""
        id_line = ["ID   CYC_BOVIN      STANDARD;      PRT;   104 AA."]
        self.assertEqual(
            id_parser(id_line),
            {
                "DataClass": "STANDARD",
                "Length": 104,
                "moltype": "PRT",
                "EntryName": "CYC_BOVIN",
            },
        )

    def test_sq_parser(self):
        """sq_parser should return expected dict"""
        lines = ["SQ   SEQUENCE   486 AA;  55639 MW;  D7862E867AD74383 CRC64;"]
        self.assertEqual(
            sq_parser(lines),
            {"Crc64": "D7862E867AD74383", "Length": 486, "MolWeight": 55639},
        )

    def test_ac_parser(self):
        """ac_parser should return expected list"""
        lines = ["AC   Q16653; O00713; O00714;", "AC   Q92892; Q92893;"]
        self.assertEqual(
            ac_parser(lines), ["Q16653", "O00713", "O00714", "Q92892", "Q92893"]
        )

    def test_oc_parser(self):
        """oc_parser should return expected list"""
        lines = ["OC   Eukaryota; Metazoa; Chordata;", "OC   Mammalia;"]
        self.assertEqual(
            oc_parser(lines), ["Eukaryota", "Metazoa", "Chordata", "Mammalia"]
        )

    def test_dt_parser(self):
        """dt_parser should return expected list of list"""
        lines = (
            "DT   01-AUG-1988 (Rel. 08, Created)\n"
            "DT   30-MAY-2000 (Rel. 39, Last sequence update)\n"
            "DT   10-MAY-2005 (Rel. 47, Last annotation update)\n".splitlines()
        )
        self.assertEqual(
            dt_parser(lines),
            [
                "01-AUG-1988 (Rel. 08, Created)",
                "30-MAY-2000 (Rel. 39, Last sequence update)",
                "10-MAY-2005 (Rel. 47, Last annotation update)",
            ],
        )

    def test_de_itemparser(self):
        """de_itemparser: known values"""
        inputs = [" AAA (aa) ", "AAA [xx] (aa)", "AAA", ""]
        expects = [
            {"OfficalName": "AAA", "Synonyms": ["aa"]},
            {"OfficalName": "AAA [xx]", "Synonyms": ["aa"]},
            {"OfficalName": "AAA", "Synonyms": []},
            {"OfficalName": "", "Synonyms": []},
        ]
        # pprint(map(de_itemparser, inputs))
        self.assertEqual(list(map(de_itemparser, inputs)), expects)

    def test_pr_parser(self):
        """pr_parser should return expected list"""
        inpr = "PR   Project:PRJNA38045; Project:PRJNA41539;"
        exp = ["PRJNA38045", "PRJNA41539"]
        obs = pr_parser(inpr)
        self.assertEqual(obs, exp)

    def test_de_parser(self):
        """de_parser should return expected list"""
        inputs = [
            "DE   Annexin [Includes: CCC] [Contains: AAA] (Fragment).",
            "DE   A [Includes: II] (Fragment).",
            "DE   A [Contains: CC].",
            "DE   A (Fragment).",
            "DE   A (aa).",
        ]
        filtered_dicts = [item_empty_filter(de_parser([e])) for e in inputs]
        self.assertEqual(list(map(len, filtered_dicts)), [4, 3, 2, 2, 2])

    def test_os_parser(self):
        """os_parser should return expected list"""
        lines = ["OS   Solanum melongena (Eggplant) (Auber-", "OS   gine)."]
        self.assertEqual(
            os_parser(lines), ["Solanum melongena", "Eggplant", "Auber- gine"]
        )

        lines = """OS   Escherichia coli.""".splitlines()
        self.assertEqual(os_parser(lines), ["Escherichia coli"])

    def test_og_parser(self):
        """og_parser should return expected list"""
        lines = [
            "OG   XXX; xx.",
            "OG   Plasmid R6-5, Plasmid IncFII R100 (NR1), and",
            "OG   Plasmid IncFII R1-19 (R1 drd-19).",
        ]
        self.assertEqual(
            og_parser(lines),
            [
                "XXX; xx",
                [
                    "Plasmid R6-5",
                    "Plasmid IncFII R100 (NR1)",
                    "Plasmid IncFII R1-19 (R1 drd-19)",
                ],
            ],
        )

    def test_ox_parser(self):
        """ox_parser should return expected dict"""
        lines = ["OX   NCBI_TaxID=9606;"]

        self.assertEqual(ox_parser(lines), {"NCBI_TaxID": "9606"})

    def test_gn_parser(self):
        """gn_parser should return expected list of dict"""

        lines = [
            "GN   name=hns; Synonyms=bglY, cur, topS;",
            "GN   OrderedLocusNames=b1237, c1701, ECs1739;",
        ]
        self.assertEqual(
            gn_parser(lines),
            [
                {
                    "Synonyms": ["bglY", "cur", "topS"],
                    "OrderedLocusNames": ["b1237", "c1701", "ECs1739"],
                    "name": "hns",
                }
            ],
        )

        lines = [
            "GN   name=Jon99Cii; Synonyms=SER1, Ser99Da; ORFNames=CG7877;",
            "GN   and",
            "GN   name=Jon99Ciii; Synonyms=SER2, Ser99Db; ORFNames=CG15519;",
        ]
        self.assertEqual(
            gn_parser(lines),
            [
                {
                    "ORFNames": "CG7877",
                    "Synonyms": ["SER1", "Ser99Da"],
                    "name": "Jon99Cii",
                },
                {
                    "ORFNames": "CG15519",
                    "Synonyms": ["SER2", "Ser99Db"],
                    "name": "Jon99Ciii",
                },
            ],
        )

    def test_dr_parser(self):
        """dr_parser should return expected dict"""
        lines = dr_lines
        self.assertEqual(dr_parser(lines), dr_expect)


class FT_Tests(TestCase):
    """Tests for FT parsers."""

    def test_ft_basic_itemparser(self):
        """ft_basic_itemparser: known values"""
        inputs = [
            ["DNA_BIND    >102    292"],
            ["CONFLICT    327    327       E -> R (in Ref. 2)."],
            [
                "PROPEP      ?25     48",
                "                             /FTId=PRO_021449.",
            ],
            [
                "VARIANT     214    214       V -> I.",
                "                             /FTId=VAR_009122.",
            ],
        ]
        expects = [
            ("DNA_BIND", ">102", 292, ""),
            ("CONFLICT", 327, 327, "E -> R (in Ref. 2)"),
            ("PROPEP", "?25", 48, "/FTId=PRO_021449"),
            ("VARIANT", 214, 214, "V -> I. /FTId=VAR_009122"),
        ]
        # pprint(map(ft_basic_itemparser, inputs))
        self.assertEqual(list(map(ft_basic_itemparser, inputs)), expects)

    def test_try_int(self):
        """try_int: known values"""
        inputs = ["9", "0", "-3", "2.3", "<9", ">9", "?", "?35", ""]
        expects = [9, 0, -3, "2.3", "<9", ">9", "?", "?35", ""]
        self.assertEqual(list(map(try_int, inputs)), expects)

    def test_ft_id_parser(self):
        """ft_id_parser: known values"""
        inputs = [
            "",
            "ddd",
            "/FTId=PRO_021449",
            "V -> I. /FTId=VAR_009122",
            "E -> R (tumor). /FTId=VAR_002343",
        ]
        expects = [
            {"Description": "", "Id": ""},
            {"Description": "ddd", "Id": ""},
            {"Description": "", "Id": "PRO_021449"},
            {"Description": "V -> I", "Id": "VAR_009122"},
            {"Description": "E -> R (tumor)", "Id": "VAR_002343"},
        ]

        # pprint(map(ft_id_parser, inputs))
        self.assertEqual(list(map(ft_id_parser, inputs)), expects)

    def test_ft_mutation_parser(self):
        """ft_mutation_parser: known values"""
        inputs = [
            "",
            "ddd",  # should raise error?
            "V -> I. /FTId=xxxxxx",  # should raise error?
            "V -> I",
            "E -> R (tumor)",
            "missing (tumor)",
        ]
        expects = [
            {"MutateFrom": "", "Comment": "", "MutateTo": ""},
            {"MutateFrom": "ddd", "Comment": "", "MutateTo": ""},
            {"MutateFrom": "V", "Comment": "", "MutateTo": "I. /FTId=xxxxxx"},
            {"MutateFrom": "V", "Comment": "", "MutateTo": "I"},
            {"MutateFrom": "E", "Comment": "tumor", "MutateTo": "R"},
            {"MutateFrom": "missing ", "Comment": "tumor", "MutateTo": ""},
        ]

        # pprint(map(ft_mutation_parser, inputs))
        self.assertEqual(list(map(ft_mutation_parser, inputs)), expects)

    def test_ft_mutation_parser_raise(self):
        """ft_mutation_parser: raise ValueError"""
        pass

    def test_ft_mutagen_parser(self):
        """ft_mutagen_parser: known values"""
        inputs = ["C->R,E,A: Loss of cADPr hydrolas", "Missing: Abolishes ATP-binding"]
        expects = [
            {
                "Comment": " Loss of cADPr hydrolas",
                "MutateFrom": "C",
                "MutateTo": "R,E,A",
            },
            {
                "Comment": " Abolishes ATP-binding",
                "MutateFrom": "Missing",
                "MutateTo": "",
            },
        ]

        # pprint(map(ft_mutagen_parser, inputs))
        self.assertEqual(list(map(ft_mutagen_parser, inputs)), expects)

    def test_ft_id_mutation_parser(self):
        """ft_id_mutation_parser: known values"""
        pass

    def test_ft_parser(self):
        """ft_parser should return expected dict"""
        lines = ft_lines
        # pprint(ft_parser(lines))
        self.assertEqual(ft_parser(lines), ft_expect)


class CC_Tests(TestCase):
    """tests for cc_parsers."""

    def test_cc_itemfinder_valid(self):
        """cc_itemfinder: yield each expected block."""
        # pprint(list(cc_itemfinder(labeloff(cc_lines))))
        input_with_license = labeloff(cc_lines)
        self.assertEqual(len(list(cc_itemfinder(input_with_license))), 9)

        input_without_license = labeloff(cc_lines[:-4])
        self.assertEqual(len(list(cc_itemfinder(input_without_license))), 8)

    def test_cc_itemfinder_raise(self):
        """cc_itemfinder: raise RecordError if license block bad."""
        input_with_license_lacking_bottom = labeloff(cc_lines[:-1])
        self.assertRaises(FieldError, cc_itemfinder, input_with_license_lacking_bottom)

    def test_cc_basic_itemparser(self):
        """cc_basic_itemparser: known results or FieldError"""
        valid_topics = [
            ["-!- topic1: first: line", "    second line"],
            ["-!- topic2: ", "    first line", "   second line"],
            [" topic3: not treated invalid topic format"],
        ]
        expects = [
            ("topic1", ["first: line", "second line"]),
            ("topic2", ["first line", "econd line"]),
            ("topic3", ["not treated invalid topic format"]),
        ]

        self.assertEqual(list(map(cc_basic_itemparser, valid_topics)), expects)

        bad_topic = ["-!- bad_topic without colon", "    FieldError"]
        self.assertRaises(FieldError, cc_basic_itemparser, bad_topic)

    def test_cc_interaction_parser(self):
        """cc_interaction_parser: known values"""
        inputs = [
            [
                "Self; NbExp=1; IntAct=EBI-123485, EBI-123485;",
                "Q9W158:CG4612; NbExp=1; IntAct=EBI-123485, EBI-89895;",
                "Q9VYI0:fne; NbExp=1; IntAct=EBI-123485, EBI-126770;",
            ]
        ]
        expects = [
            [
                ("Self", {"NbExp": "1", "IntAct": ["EBI-123485", "EBI-123485"]}),
                (
                    "Q9W158:CG4612",
                    {"NbExp": "1", "IntAct": ["EBI-123485", "EBI-89895"]},
                ),
                ("Q9VYI0:fne", {"NbExp": "1", "IntAct": ["EBI-123485", "EBI-126770"]}),
            ]
        ]
        self.assertEqual(list(map(cc_interaction_parser, inputs)), expects)

    def test_cc_biophysicochemical_properties_parser(self):
        """cc_biophysicochemical_properties_parser: known values"""
        # pprint(cc['BIOPHYSICOCHEMICAL PROPERTIES'])  #topic specific parser
        f = cc_biophysicochemical_properties_parser
        valid_inputs = [
            [
                "Kinetic parameters:",
                "  KM=98 uM for ATP;",
                "  KM=688 uM for pyridoxal;",
                "  Vmax=1.604 mmol/min/mg enzyme;",
                "pH dependence:",
                "  Optimum pH is 6.0. Active pH 4.5 to 10.5;",
            ]
        ]
        expects = [
            {
                "Kinetic parameters": {
                    "KM": ["98 uM for ATP", "688 uM for pyridoxal"],
                    "Vmax": "1.604 mmol/min/mg enzyme",
                },
                "pH dependence": "Optimum pH is 6.0. Active pH 4.5 to 10.5",
            }
        ]
        self.assertEqual(list(map(f, valid_inputs)), expects)

    def test_cc_alternative_products_parser(self):
        """cc_alternative_products_parser: know values"""
        f = cc_alternative_products_parser
        valid_inputs = [
            [
                "Event=Alternative initiation;" "  Comment=Free text;",
                "Event=Alternative splicing; Named isoforms=3;",
                "  Comment=Additional isoforms seem to exist.",
                "  confirmation;",
                "name=1; Synonyms=AIRE-1;",
                "  IsoId=O43918-1; Sequence=Displayed;",
                "name=3; Synonyms=AIRE-3,",
                "ai-2, ai-3;",  # broken the hanging_paragraph_finder
                "  IsoId=O43918-3; Sequence=VSP_004089, VSP_004090;",
            ]
        ]
        expects = [
            [
                {"Comment": "Free text", "Event": "Alternative initiation"},
                {
                    "Comment": "Additional isoforms seem to exist. confirmation",
                    "Event": "Alternative splicing",
                    "Named isoforms": "3",
                    "Names": [
                        {
                            "IsoId": "O43918-1",
                            "name": "1",
                            "Sequence": "Displayed",
                            "Synonyms": "AIRE-1",
                        },
                        {
                            "IsoId": "O43918-3",
                            "name": "3",
                            "Sequence": ["VSP_004089", "VSP_004090"],
                            "Synonyms": ["AIRE-3", "ai-2", "ai-3"],
                        },
                    ],
                },
            ]
        ]

        # pprint(map(f,valid_inputs))
        self.assertEqual(list(map(f, valid_inputs)), expects)

    def test_cc_parser(self):
        """cc_parser: known values and raise when strict"""
        cc = cc_parser(cc_lines)
        # pprint(cc)
        # print cc.keys()
        self.assertEqual(
            list(sorted(cc.keys())),
            [
                "ALLERGEN",
                "ALTERNATIVE PRODUCTS",
                "BIOPHYSICOCHEMICAL PROPERTIES",
                "DATABASE",
                "DISEASE",
                "INTERACTION",
                "LICENSE",
                "MASS SPECTROMETRY",
            ],
        )

        # test Disease topic (default_handler)
        self.assertEqual(
            cc["DISEASE"],
            [
                "Defects in PHKA1 are linked to X-linked muscle glycogenosis "
                "[MIM:311870]",
                "Defects in ABCD1 are the cause of recessive X-linked "
                "adrenoleukodystrophy (X-ALD) [MIM:300100]. X-ALD is a rare "
                "phenotype",
            ],
        )

        # test License (default_handler)
        # pprint(cc['LICENSE'])
        self.assertEqual(
            cc["LICENSE"],
            [
                "This SWISS-PROT entry is copyright. It is produced through a "
                "collaboration removed"
            ],
        )

        # pprint(cc['DATABASE'])  #join_split_dict
        self.assertEqual(
            cc["DATABASE"],
            [
                {
                    "NAME": "CD40Lbase",
                    "NOTE": "European CD40L defect database (mutation db)",
                    "WWW": '"http://www.expasy.org/cd40lbase/"',
                }
            ],
        )

        # test strict
        cc_lines_with_unknown_topic = ["CC   -!- BLAHBLAH: xxxxx"] + cc_lines
        # pprint(cc_parser(cc_lines_with_unknown_topic))
        self.assertEqual(cc_parser(cc_lines_with_unknown_topic)["BLAHBLAH"], ["xxxxx"])
        self.assertRaises(
            FieldError, cc_parser, cc_lines_with_unknown_topic, strict=True
        )


class ReferenceTests(TestCase):
    """Tests for parsers related to reference blocks"""

    def test_ref_finder(self):
        """ref_finder: should return a list of ref blocks"""
        pass

    def test_refs_parser(self):
        """refs_parser: should return a dict of {RN: ref_dict}"""
        pass

    def test_single_ref_parser(self):
        """single_ref_parser: should return the expected dict"""
        fake_ref_block = [
            "RN   [1]",
            "RP   NUCLEOTIDE",
            "RC   STRAIN=Bristol N2;",
            "RX   PubMed=1113;",
            "RA   Zhang L., Wu S.-L., Rubin C.S.;",
            'RT   "A novel ";',
            "RL   J. Biol. Chem. 276:10.",
        ]
        rn, others = single_ref_parser(fake_ref_block)
        self.assertEqual(rn, 1)
        self.assertEqual(len(others), 6)

        # test strict: lacking required labels
        self.assertEqual(
            len(single_ref_parser(fake_ref_block[:-1], strict=False)[1]), 5
        )
        self.assertRaises(RecordError, single_ref_parser, fake_ref_block[:-1], True)

    def test_ra_parser(self):
        """ra_parser should return expected list"""
        lines = (
            "RA   Galinier A., Bleicher F., Negre D.,\n"
            "RA   Cozzone A.J., Cortay J.-C.;\n".splitlines()
        )
        self.assertEqual(
            ra_parser(lines),
            ["Galinier A.", "Bleicher F.", "Negre D.", "Cozzone A.J.", "Cortay J.-C."],
        )

    def test_rx_parser(self):
        """rx_parser should return expected dict"""
        inputs = [
            ["RX   MEDLINE=22709107; PubMed=12788972; DOI=10.1073/pnas.113"],
            [
                "RX   PubMed=14577811; "
                "DOI=10.1597/1545-1569(2003)040<0632:AMMITS>2.0.CO;2;"
            ],
        ]
        expects = [
            {"DOI": "10.1073/pnas.113", "MEDLINE": "22709107", "PubMed": "12788972"},
            {
                "DOI": "10.1597/1545-1569(2003)040<0632:AMMITS>2.0.CO;2",
                "PubMed": "14577811",
            },
        ]

        self.assertEqual(list(map(rx_parser, inputs)), expects)

    def test_rc_parser(self):
        """rc_parser should return expected dict"""
        lines = [
            "RC   PLASMID=R1 (R7268); TRANSPOSON=Tn3;",
            "RC   STRAIN=AL.012, AZ.026;",
        ]

        self.assertEqual(
            rc_parser(lines),
            {
                "TRANSPOSON": "Tn3",
                "PLASMID": "R1 (R7268)",
                "STRAIN": ["AL.012", "AZ.026"],
            },
        )

    def test_rt_parser(self):
        """rt_parser should return expected str"""
        lines = [
            'RT   "New insulin-like proteins',
            'RT   analysis and homology modeling.";',
        ]
        self.assertEqual(
            rt_parser(lines), "New insulin-like proteins analysis and homology modeling"
        )

    def test_rl_parser(self):
        """rl_parser should return expected str"""
        lines = ["RL   J. Mol. Biol. 168:321-331(1983)."]
        self.assertEqual(rl_parser(lines), "J. Mol. Biol. 168:321-331(1983)")

    def test_rn_parser(self):
        """rn_parser should return expected str"""
        lines = ["RN   [8]"]
        self.assertEqual(rn_parser(lines), 8)

    def test_rg_parser(self):
        """rg_parser should return expected str"""
        lines = ["RG   The mouse genome sequencing consortium;"]
        self.assertEqual(rg_parser(lines), ["The mouse genome sequencing consortium"])

    def test_rp_parser(self):
        """rp_parser should return expected str"""
        lines = ["RP   X-RAY CRYSTALLOGRAPHY (1.8 ANGSTROMS)."]
        self.assertEqual(rp_parser(lines), "X-RAY CRYSTALLOGRAPHY (1.8 ANGSTROMS)")


#################################
# global test data
ft_lines = """FT   CHAIN        29    262       Granzyme A.
FT                                /FTId=PRO_0000027394.
FT   ACT_SITE     69     69       Charge relay system.
FT   VARIANT     121    121       T -> M (in dbSNP:3104233).
FT                                /FTId=VAR_024291.
FT   VARIANT       1      7       unknown  (in a skin tumor).
FT                                /FTId=VAR_005851.
FT   CONFLICT    282    282       R -> Q (in Ref. 18).
FT   STRAND       30     30
FT   STRAND       33     34
FT   TURN         37     38""".splitlines()

ft_expect = {
    "ACT_SITE": [{"Start": 69, "End": 69, "Description": "Charge relay system"}],
    "CHAIN": [
        {
            "Description": {"Description": "Granzyme A", "Id": "PRO_0000027394"},
            "End": 262,
            "Start": 29,
        }
    ],
    "CONFLICT": [
        {
            "Description": {
                "Comment": "in Ref. 18",
                "MutateFrom": "R",
                "MutateTo": "Q",
            },
            "End": 282,
            "Start": 282,
        }
    ],
    "SecondaryStructure": [("STRAND", 30, 30), ("STRAND", 33, 34), ("TURN", 37, 38)],
    "VARIANT": [
        {
            "Description": {
                "Comment": "in dbSNP:3104233",
                "Id": "VAR_024291",
                "MutateFrom": "T",
                "MutateTo": "M",
            },
            "End": 121,
            "Start": 121,
        },
        {
            "Description": {
                "Comment": "in a skin tumor",
                "Id": "VAR_005851",
                "MutateFrom": "unknown  ",
                "MutateTo": "",
            },
            "End": 7,
            "Start": 1,
        },
    ],
}

dr_lines = """DR   MIM; 140050; gene.
DR   GO; GO:0001772; C:immunological synapse; TAS.
DR   GO; GO:0005634; C:nucleus; TAS.
DR   GO; GO:0006915; P:apoptosis; TAS.
DR   GO; GO:0006922; P:cleavage of lamin; IDA.
DR   GO; GO:0006955; P:immune response; TAS.
DR   InterPro; IPR001254; Peptidase_S1_S6.
DR   InterPro; IPR001314; Peptidase_S1A.
DR   Pfam; PF00089; Trypsin; 1.""".splitlines()

dr_expect = {
    "GO": [
        ["GO:0001772", "C:immunological synapse", "TAS"],
        ["GO:0005634", "C:nucleus", "TAS"],
        ["GO:0006915", "P:apoptosis", "TAS"],
        ["GO:0006922", "P:cleavage of lamin", "IDA"],
        ["GO:0006955", "P:immune response", "TAS"],
    ],
    "Pfam": [["PF00089", "Trypsin", "1"]],
    "InterPro": [["IPR001254", "Peptidase_S1_S6"], ["IPR001314", "Peptidase_S1A"]],
    "MIM": [["140050", "gene"]],
}

cc_lines = """CC   -!- ALLERGEN: Causes an allergic reaction in human. Binds to IgE.
CC       bovine dander.
CC   -!- ALTERNATIVE PRODUCTS:
CC       Event=Alternative splicing; Named isoforms=3;
CC         Comment=Additional isoforms seem to exist.
CC         confirmation;
CC       name=1; Synonyms=AIRE-1;
CC         IsoId=O43918-1; Sequence=Displayed;
CC       name=2; Synonyms=AIRE-2;
CC         IsoId=O43918-2; Sequence=VSP_004089;
CC       name=3; Synonyms=AIRE-3;
CC         IsoId=O43918-3; Sequence=VSP_004089, VSP_004090;
CC   -!- BIOPHYSICOCHEMICAL PROPERTIES:
CC       Kinetic parameters:
CC         KM=98 uM for ATP;
CC         KM=688 uM for pyridoxal;
CC         Vmax=1.604 mmol/min/mg enzyme;
CC       pH dependence:
CC         Optimum pH is 6.0. Active pH 4.5 to 10.5;
CC   -!- DATABASE: NAME=CD40Lbase;
CC       NOTE=European CD40L defect database (mutation db);
CC       WWW="http://www.expasy.org/cd40lbase/".
CC   -!- DISEASE: Defects in PHKA1 are linked to X-linked muscle
CC       glycogenosis [MIM:311870]. 
CC   -!- DISEASE: Defects in ABCD1 are the cause of recessive X-linked
CC       adrenoleukodystrophy (X-ALD) [MIM:300100]. X-ALD is a rare
CC       phenotype.
CC   -!- INTERACTION:
CC       Self; NbExp=1; IntAct=EBI-123485, EBI-123485;
CC       Q9W158:CG4612; NbExp=1; IntAct=EBI-123485, EBI-89895;
CC       Q9VYI0:fne; NbExp=1; IntAct=EBI-123485, EBI-126770;
CC   -!- MASS SPECTROMETRY: MW=24948; MW_ERR=6; METHOD=MALDI; RANGE=1-228;
CC       NOTE=Ref.2.
CC   --------------------------------------------------------------------------
CC   This SWISS-PROT entry is copyright. It is produced through a collaboration
CC   removed.
CC   --------------------------------------------------------------------------
""".splitlines()

fake_records_valid = """ID   CYC_BOVIN      STANDARD;      PRT;   104 AA.
AC   ac1; ac2;
DT   dt.
DE   de.
OS   os.
OC   oc.
OX   NCBI_TaxID=9606;
RN   [1]
SQ   SEQUENCE   486 AA;  55639 MW;  D7862E867AD74383 CRC64
     aac cpp
     ggh hh
//
ID   idid  std; prt; 104    #-5
OX   NCBI_TaxID=9606;
DE   dede.
     ggaaccpp
//""".splitlines()


# Run tests if called from the command line
if __name__ == "__main__":
    main()
