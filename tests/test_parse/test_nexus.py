"""Unit tests for the Nexus Parser"""

from unittest import TestCase

from cogent3 import load_aligned_seqs
from cogent3.parse.nexus import (
    Partition,
    iter_nexus_align_records,
    find_fields,
    get_BL_table,
    get_tree_info,
    parse_dnd,
    parse_nexus_tree,
    parse_PAUP_log,
    parse_taxa,
    parse_trans_table,
    split_tree_info,
)

# ruff: noqa: SIM905

Nexus_tree = """#NEXUS

Begin trees;  [Treefile saved Wednesday, May 5, 2004  5:02 PM]
[!
>Data file = Grassland_short.nex
>Neighbor-joining search settings:
>  Ties (if encountered) will be broken systematically
>  Distance measure = Jukes-Cantor
>  (Tree is unrooted)
]
	Translate
		1 outgroup25,
		2 AF078391l,
		3 AF078211af,
		4 AF078393l,
		5 AF078187af,
		6 AF078320l,
		7 AF078432l,
		8 AF078290af,
		9 AF078350l,
		10 AF078356l,
		11 AF078306af,
		12 AF078429l,
		13 AF078256af,
		14 AF078443l,
		15 AF078450l,
		16 AF078452l,
		17 AF078258af,
		18 AF078380l,
		19 AF078251af,
		20 AF078179af,
		21 outgroup258
		;
tree PAUP_1 = [&R] (1,(2,(((3,4),(5,(((((6,10),9),(11,18)),((((7,15),19),17),(8,(12,(14,16))))),13))),20)),21);
tree PAUP_2 = [&R] (1,(2,(((3,4),(5,(((((6,10),9),(11,18)),((((7,15),19),17),(8,(12,(14,16))))),13))),20)),21);
End;""".split("\n")

Nexus_tree_2 = """#NEXUS

Begin trees;  [Treefile saved Wednesday, June 14, 2006  11:20 AM]
[!>Neighbor-joining search settings:
>  Ties (if encountered) will be broken systematically
>  Distance measure = uncorrected ("p")
>  (Tree is unrooted)
]
tree nj = [&U] ((((((((((YA10260L1:0.01855,SARAG06_Y:0.00367):0.01965,(((YA270L1G0:0.01095,SARAD10_Y:0.00699):0.01744,YA270L1A0:0.04329):0.00028,((YA165L1C1:0.01241,SARAA02_Y:0.02584):0.00213,((YA165L1H0:0.00092,SARAF10_Y:-0.00092):0.00250,(YA165L1A0:0.00177,SARAH10_Y:0.01226):0.00198):0.00131):0.00700):0.01111):0.11201,(YA160L1F0:0.00348,SARAG01_Y:-0.00122):0.13620):0.01202,((((YRM60L1D0:0.00357,(YRM60L1C0:0.00477,SARAE10_Y:-0.00035):0.00086):0.00092,SARAE03_Y:0.00126):0.00125,SARAC11_Y:0.00318):0.00160,YRM60L1H0:0.00593):0.09975):0.07088,SARAA01_Y:0.02880):0.00190,SARAB04_Y:0.05219):0.00563,YRM60L1E0:0.06099):0.00165,(YRM60L1H0:0.00450,SARAF11_Y:0.01839):0.00288):0.00129,YRM60L1B1:0.00713):0.00194,(YRM60L1G0:0.00990,(YA165L1G0:0.00576,(YA160L1G0:0.01226,SARAA11_Y:0.00389):0.00088):0.00300):0.00614,SARAC06_Y:0.00381);
end;""".split("\n")

Nexus_tree_3 = """#NEXUS

Begin trees;  [Treefile saved Wednesday, May 5, 2004  5:02 PM]
[!
>Data file = Grassland_short.nex
>Neighbor-joining search settings:
>  Ties (if encountered) will be broken systematically
>  Distance measure = Jukes-Cantor
>  (Tree is unrooted)
]
	Translate
		1 outgroup25,
		2 AF078391l,
		3 'AF078211af',
		4 AF078393l,
		5 AF078187af,
		6 AF078320l,
		7 AF078432l,
		8 AF078290af,
		9 AF078350l,
		10 AF078356l,
		11 AF078306af,
		12 AF078429l,
		13 AF078256af,
		14 'AF078443l',
		15 AF078450l,
		16 AF078452l,
		17 AF078258af,
		18 'AF078380l',
		19 AF078251af,
		20 AF078179af,
		21 outgroup258
		;
tree PAUP_1 = [&R] (1,(2,(((3,4),(5,(((((6,10),9),(11,18)),((((7,15),19),17),(8,(12,(14,16))))),13))),20)),21);
tree PAUP_2 = [&R] (1,(2,(((3,4),(5,(((((6,10),9),(11,18)),((((7,15),19),17),(8,(12,(14,16))))),13))),20)),21);
End;""".split("\n")

PAUP_log = """
P A U P *
Version 4.0b10 for Macintosh (PPC/Altivec)
Wednesday, May 5, 2004  5:03 PM

This copy registered to: Scott Dawson
                         UC-Berkeley
                         (serial number = B400784)

      -----------------------------NOTICE-----------------------------
        This is a beta-test version.  Please report any crashes,
        apparent calculation errors, or other anomalous results.
        There are no restrictions on publication of results obtained
        with this version, but you should check the WWW site
        frequently for bug announcements and/or updated versions.
        See the README file on the distribution media for details.
      ----------------------------------------------------------------

Tree description:

  Optimality criterion = parsimony
    Character-status summary:
      Of 500 total characters:
        All characters are of type 'unord'
        All characters have equal weight
        253 characters are constant
        109 variable characters are parsimony-uninformative
        Number of parsimony-informative characters = 138
    Multistate taxa interpreted as uncertainty
    Character-state optimization: Accelerated transformation (ACCTRAN)
                                  AncStates = "standard"

Tree number 1 (rooted using user-specified outgroup)

Branch lengths and linkages for tree #1

                                     Assigned       Minimum       Maximum
                    Connected          branch      possible      possible
   Node              to node           length        length        length
-------------------------------------------------------------------------
    40                root                  0             0             0
outgroup25 (1)*         40                 40            24            52
    39                  40                 57            15            72
AF078391l (2)           39                 56            48            81
    38                  39                 33            17            71
    37                  38                 31            14            48
    22                  37                 20            11            33
AF078211af (3)          22                  4             2             7
AF078393l (4)           22                  1             0             3
    36                  37                 14             5            32
AF078187af (5)          36                 18            10            28
    35                  36                 21            16            45
    34                  35                 10             3            23
    26                  34                  5             3             9
    24                  26                  4             3            13
    23                  24                  0             0             3
AF078320l (6)           23                  1             1             3
AF078356l (10)          23                  2             2             2
AF078350l (9)           24                  5             3             5
    25                  26                  9             2            10
AF078306af (11)         25                  6             4            10
AF078380l (18)          25                  5             3            10
    33                  34                  5             4            15
    29                  33                  3             1             4
    28                  29                  2             2             2
    27                  28                  3             1             3
AF078432l (7)           27                  2             2             2
AF078450l (15)          27                  3             3             4
AF078251af (19)         28                  6             6             7
AF078258af (17)         29                  6             6             6
    32                  33                  4             3            15
AF078290af (8)          32                  9             8            11
    31                  32                  9             6            18
AF078429l (12)          31                  2             1             5
    30                  31                 10             9            18
AF078443l (14)          30                  2             1             6
AF078452l (16)          30                  4             4             5
AF078256af (13)         35                  4             1             6
AF078179af (20)         38                 48            34            79
outgroup258 (21)*       40                 45            27            67
-------------------------------------------------------------------------
Sum                                       509

Tree length = 509
Consistency index (CI) = 0.7151
Homoplasy index (HI) = 0.2849
""".split("\n")

line1 = "    40                root                  0             0           0"
line2 = "outgroup25 (1)*         40                 40            24            52"
line3 = "    39                  40                 57            15            72"
line4 = "AF078391l (2)           39                 56            48            81"


class NexusParserTests(TestCase):
    """Tests of the Nexus Parser functions"""

    def test_parse_nexus_tree(self):
        """parse_nexus_tree returns a dnd string and a translation table list"""
        Trans_table, dnd = parse_nexus_tree(Nexus_tree)

        # check the full dendrogram string is returned
        assert (
            dnd["tree PAUP_1"]
            == "(1,(2,(((3,4),(5,(((((6,10),9),(11,18)),((((7,15),19),17),(8,(12,(14,16))))),13))),20)),21);"
        )

        # check that all taxa are returned in the Trans_table
        assert Trans_table["1"] == "outgroup25"
        assert Trans_table["2"] == "AF078391l"
        assert Trans_table["3"] == "AF078211af"
        assert Trans_table["4"] == "AF078393l"
        assert Trans_table["5"] == "AF078187af"
        assert Trans_table["6"] == "AF078320l"
        assert Trans_table["21"] == "outgroup258"
        assert Trans_table["20"] == "AF078179af"
        assert Trans_table["19"] == "AF078251af"

        # check that Nexus files without translation table work
        Trans_table, dnd = parse_nexus_tree(Nexus_tree_2)
        assert Trans_table is None
        assert (
            dnd["tree nj"]
            == "((((((((((YA10260L1:0.01855,SARAG06_Y:0.00367):0.01965,(((YA270L1G0:0.01095,SARAD10_Y:0.00699):0.01744,YA270L1A0:0.04329):0.00028,((YA165L1C1:0.01241,SARAA02_Y:0.02584):0.00213,((YA165L1H0:0.00092,SARAF10_Y:-0.00092):0.00250,(YA165L1A0:0.00177,SARAH10_Y:0.01226):0.00198):0.00131):0.00700):0.01111):0.11201,(YA160L1F0:0.00348,SARAG01_Y:-0.00122):0.13620):0.01202,((((YRM60L1D0:0.00357,(YRM60L1C0:0.00477,SARAE10_Y:-0.00035):0.00086):0.00092,SARAE03_Y:0.00126):0.00125,SARAC11_Y:0.00318):0.00160,YRM60L1H0:0.00593):0.09975):0.07088,SARAA01_Y:0.02880):0.00190,SARAB04_Y:0.05219):0.00563,YRM60L1E0:0.06099):0.00165,(YRM60L1H0:0.00450,SARAF11_Y:0.01839):0.00288):0.00129,YRM60L1B1:0.00713):0.00194,(YRM60L1G0:0.00990,(YA165L1G0:0.00576,(YA160L1G0:0.01226,SARAA11_Y:0.00389):0.00088):0.00300):0.00614,SARAC06_Y:0.00381);"
        )

    def test_parse_nexus_tree_sq(self):
        """remove single quotes from tree and translate tables"""
        Trans_table, dnd = parse_nexus_tree(Nexus_tree_3)

        # check the full dendrogram string is returned
        assert (
            dnd["tree PAUP_1"]
            == "(1,(2,(((3,4),(5,(((((6,10),9),(11,18)),((((7,15),19),17),(8,(12,(14,16))))),13))),20)),21);"
        )

        # check that all taxa are returned in the Trans_table
        assert Trans_table["1"] == "outgroup25"
        assert Trans_table["2"] == "AF078391l"
        assert Trans_table["3"] == "AF078211af"
        assert Trans_table["4"] == "AF078393l"
        assert Trans_table["5"] == "AF078187af"
        assert Trans_table["6"] == "AF078320l"
        assert Trans_table["21"] == "outgroup258"
        assert Trans_table["20"] == "AF078179af"
        assert Trans_table["19"] == "AF078251af"

    def test_get_tree_info(self):
        """get_tree_info returns the Nexus file section that describes the tree"""
        result = get_tree_info(Nexus_tree)
        assert len(result) == 33
        assert (
            result[0]
            == "Begin trees;  [Treefile saved Wednesday, May 5, 2004  5:02 PM]"
        )
        assert (
            result[31]
            == "tree PAUP_1 = [&R] (1,(2,(((3,4),(5,(((((6,10),9),(11,18)),((((7,15),19),17),(8,(12,(14,16))))),13))),20)),21);"
        )

    def test_split_tree_info(self):
        """split_tree_info splits lines into header, Trans_table, and dnd"""
        tree_info = get_tree_info(Nexus_tree)
        header, trans_table, dnd = split_tree_info(tree_info)

        assert len(header) == 9

        assert len(trans_table) == 22

        assert len(dnd) == 2

        assert (
            header[0]
            == "Begin trees;  [Treefile saved Wednesday, May 5, 2004  5:02 PM]"
        )
        assert header[8] == "\tTranslate"

        assert trans_table[0] == "\t\t1 outgroup25,"
        assert trans_table[21] == "\t\t;"

        assert (
            dnd[0]
            == "tree PAUP_1 = [&R] (1,(2,(((3,4),(5,(((((6,10),9),(11,18)),((((7,15),19),17),(8,(12,(14,16))))),13))),20)),21);"
        )

    def test_parse_trans_table(self):
        """parse_trans_table returns a dict with the taxa names indexed by number"""
        tree_info = get_tree_info(Nexus_tree)
        header, trans_table, dnd = split_tree_info(tree_info)
        Trans_table = parse_trans_table(trans_table)

        assert len(Trans_table) == 21

        # check that taxa are returned in the Trans_table
        assert Trans_table["1"] == "outgroup25"
        assert Trans_table["2"] == "AF078391l"
        assert Trans_table["3"] == "AF078211af"
        assert Trans_table["4"] == "AF078393l"
        assert Trans_table["5"] == "AF078187af"
        assert Trans_table["6"] == "AF078320l"
        assert Trans_table["21"] == "outgroup258"
        assert Trans_table["20"] == "AF078179af"
        assert Trans_table["19"] == "AF078251af"

    def test_parse_dnd(self):
        """parse_dnd returns a dict with dnd indexed by tree name"""
        tree_info = get_tree_info(Nexus_tree)
        header, trans_table, dnd = split_tree_info(tree_info)
        dnd_dict = parse_dnd(dnd)
        assert (
            dnd_dict["tree PAUP_1"]
            == "(1,(2,(((3,4),(5,(((((6,10),9),(11,18)),((((7,15),19),17),(8,(12,(14,16))))),13))),20)),21);"
        )

    # ------------------------------------------------------

    def test_get_BL_table(self):
        """get_BL_table returns the section of the log file w/ the BL table"""
        BL_table = get_BL_table(PAUP_log)
        assert len(BL_table) == 40
        assert (
            BL_table[0]
            == "    40                root                  0             0             0"
        )
        assert (
            BL_table[39]
            == "outgroup258 (21)*       40                 45            27            67"
        )

    def test_find_fields(self):
        """find_fields takes BL table line and returns field names mapped to info"""
        result = find_fields(line1)
        assert result["taxa"] == "40"
        assert result["bl"] == "0"
        assert result["parent"] == "root"

    def test_parse_taxa(self):
        """parse_taxa should return the taxa # from a taxa_field from find_fields"""
        result1 = find_fields(line1)
        result2 = find_fields(line2)
        result3 = find_fields(line3)
        result4 = find_fields(line4)

        assert parse_taxa(result1["taxa"]) == "40"
        assert parse_taxa(result2["taxa"]) == "1"
        assert parse_taxa(result3["taxa"]) == "39"
        assert parse_taxa(result4["taxa"]) == "2"

    def test_parse_PAUP_log(self):
        """parse_PAUP_log extracts branch length info from a PAUP log file"""
        BL_dict = parse_PAUP_log(PAUP_log)
        assert len(BL_dict) == 40
        assert BL_dict["1"] == ("40", 40)
        assert BL_dict["40"] == ("root", 0)
        assert BL_dict["39"] == ("40", 57)
        assert BL_dict["2"] == ("39", 56)
        assert BL_dict["26"] == ("34", 5)
        assert BL_dict["21"] == ("40", 45)

    def test_align_with_comments(self):
        """correctly handle an alignment block containing comments"""
        parser = iter_nexus_align_records("data/nexus_comments.nex")
        got = dict(parser)
        expect = {
            "Ephedra": "TTAAGCCATGCATGTCTAAGTATGAACTAATTCCAAACGGTGA",
            "Gnetum": "TTAAGCCATGCATGTCTATGTACGAACTAATC-AGAACGGTGA",
            "Welwitschia": "TTAAGCCATGCACGTGTAAGTATGAACTAGTC-GAAACGGTGA",
            "Ginkgo": "TTAAGCCATGCATGTGTAAGTATGAACTCTTTACAGACTGTGA",
            "Pinus": "TTAAGCCATGCATGTCTAAGTATGAACTAATTGCAGACTGTGA",
        }
        assert got == expect

    def test_align_with_spaced_seqs(self):
        """correctly handle an alignment block with spaces in seqs"""
        parser = iter_nexus_align_records("data/nexus_dna.nex")
        seqs = dict(parser)
        assert len(seqs) == 10  # 10 taxa
        lengths = {len(seqs[n]) for n in seqs}
        assert lengths == {705}  # all same length and equal 705

    def test_align_from_mixed(self):
        """correctly handle a file with tree and alignment block"""
        parser = iter_nexus_align_records("data/nexus_mixed.nex")
        got = dict(parser)
        expect = {
            "fish": "ACATAGAGGGTACCTCTAAG",
            "frog": "ACATAGAGGGTACCTCTAAG",
            "snake": "ACATAGAGGGTACCTCTAAG",
            "mouse": "ACATAGAGGGTACCTCTAAG",
        }
        assert got == expect

    def test_align_no_blank_columns(self):
        """correctly handle a file with no white space at line starts"""
        parser = iter_nexus_align_records("data/nexus_aa.nxs")
        seqs = dict(parser)
        assert len(seqs) == 10  # 10 taxa
        lengths = {len(seqs[n]) for n in seqs}
        assert lengths == {234}  # all same length and equal 234

    def test_load_seqs_interface(self):
        """load_aligned_seqs correctly loads nexus alignments"""
        aln = load_aligned_seqs("data/nexus_mixed.nex", moltype="text")
        assert aln.num_seqs == 4
        assert len(aln) == 20

        aln = load_aligned_seqs("data/nexus_aa.nxs", moltype="text")
        assert aln.num_seqs == 10
        assert len(aln) == 234


def test_partition_fields():
    # a charset "part1 = 1-100" becomes zero-based slice(0, 100)
    part = Partition(name="part1", ranges=(slice(0, 100),))
    assert part.name == "part1"
    assert part.ranges == (slice(0, 100),)
    assert part.data_type is None
    assert part.model is None
    assert part.align_name is None
    assert part.tree_len is None


def test_range_to_slice():
    from cogent3.parse.nexus import range_to_slice

    # simple range: 1-based inclusive 1-100 -> zero-based half-open [0, 100)
    assert range_to_slice("1-100") == slice(0, 100)
    # a stepped range keeps its step; start shifts by one, stop unchanged
    assert range_to_slice(r"101-249\3") == slice(100, 249, 3)
    # "*" means the whole file -> open-ended from the start
    assert range_to_slice("*") == slice(0, None)
    # "." means "to the end" -> open-ended stop
    assert range_to_slice("1500-.") == slice(1499, None)
    # a bare single site 5 -> just that one column
    assert range_to_slice("5") == slice(4, 5)


def test_parse_charset():
    from cogent3.parse.nexus import parse_charset

    # multiple space-separated ranges, no file, no type
    assert parse_charset("part1 = 1-100 250-384") == Partition(
        name="part1", ranges=(slice(0, 100), slice(249, 384))
    )
    # a stepped codon-position split
    assert parse_charset(r"part2 = 101-249\3 102-249\3") == Partition(
        name="part2", ranges=(slice(100, 249, 3), slice(101, 249, 3))
    )
    # names its own alignment file
    assert parse_charset("dna_gene1 = aln_dna.phy: 1-600") == Partition(
        name="dna_gene1", ranges=(slice(0, 600),), align_name="aln_dna.phy"
    )
    # alignment file plus a data type (comma after the type)
    assert parse_charset("nuclear_cds = aln1.phy:CODON, 1-900") == Partition(
        name="nuclear_cds",
        ranges=(slice(0, 900),),
        data_type="CODON",
        align_name="aln1.phy",
    )
    # data type given without a file (type still followed by a comma)
    assert parse_charset("cds = CODON, 1-900") == Partition(
        name="cds", ranges=(slice(0, 900),), data_type="CODON"
    )
    # "*" whole-file range with a file reference
    assert parse_charset("prot_gene1 = aln_protein.phy: *") == Partition(
        name="prot_gene1", ranges=(slice(0, None),), align_name="aln_protein.phy"
    )


def test_parse_charpartition():
    from cogent3.parse.nexus import parse_charpartition

    # basic: a model assigned to each charset
    assert parse_charpartition("mine = HKY+G:part1, GTR+I+G:part2") == (
        "mine",
        {
            "part1": {"model": "HKY+G", "data_type": None, "tree_len": None},
            "part2": {"model": "GTR+I+G", "data_type": None, "tree_len": None},
        },
    )
    # model-free (file 07): the slot holds a data type, not a model
    assert parse_charpartition("mf = DNA:gene_dna, CODON:cds") == (
        "mf",
        {
            "gene_dna": {"model": None, "data_type": "DNA", "tree_len": None},
            "cds": {"model": None, "data_type": "CODON", "tree_len": None},
        },
    )
    # a {value} branch length attached to the model
    assert parse_charpartition("s = HKY+G{0.5}:part1") == (
        "s",
        {"part1": {"model": "HKY+G", "data_type": None, "tree_len": 0.5}},
    )


def test_partition_set():
    from cogent3.parse.nexus import PartitionSet

    p1 = Partition(name="part1", ranges=(slice(0, 100),), model="HKY+G")
    p2 = Partition(name="part2", ranges=(slice(100, 384),), model="GTR+I+G")
    pset = PartitionSet(name="mine", partitions=(p1, p2))

    assert pset.name == "mine"
    assert len(pset) == 2
    # iterates in order
    assert tuple(pset) == (p1, p2)
    # addressable by partition name
    assert pset["part2"] is p2


def test_parse_nexus_partitions_basic():
    from cogent3.parse.nexus import parse_nexus_partitions

    pset = parse_nexus_partitions("data/01_basic_dna.nex")
    assert pset.name == "mine"
    assert len(pset) == 2
    assert pset["part1"] == Partition(
        name="part1", ranges=(slice(0, 100),), model="HKY+G"
    )
    assert pset["part2"] == Partition(
        name="part2", ranges=(slice(100, 384),), model="GTR+I+G"
    )


def test_parse_nexus_partitions_mixed_files():
    from cogent3.parse.nexus import parse_nexus_partitions

    # file 03: several alignment files, DNA + protein, "*" whole file
    pset = parse_nexus_partitions("data/03_mixed_dna_protein.nex")
    assert pset.name == "mixed"
    assert pset["dna_gene1"] == Partition(
        name="dna_gene1",
        ranges=(slice(0, 600),),
        model="GTR+G",
        align_name="aln_dna.phy",
    )
    assert pset["prot_gene1"] == Partition(
        name="prot_gene1",
        ranges=(slice(0, None),),
        model="WAG+I+G",
        align_name="aln_protein.phy",
    )

    # file 04: data_type comes from the charset (CODON / CODON5)
    pset = parse_nexus_partitions("data/04_codon_models_genetic_code.nex")
    assert pset["nuclear_cds"].data_type == "CODON"
    assert pset["mito_cds"].data_type == "CODON5"
    assert pset["nuclear_cds"].model == "GY"

    # file 07: no models; data_type comes from the charpartition slot
    pset = parse_nexus_partitions("data/07_model_free.nex")
    assert pset["gene_dna"] == Partition(
        name="gene_dna",
        ranges=(slice(0, 1200),),
        data_type="DNA",
        align_name="aln_dna.phy",
    )
    assert pset["cds"].data_type == "CODON"
    assert pset["cds"].model is None


def test_parse_nexus_partitions_data_block_and_comments():
    from cogent3.parse.nexus import parse_nexus_partitions

    # file 08: a DATA block precedes the SETS block and must be skipped
    pset = parse_nexus_partitions("data/08_data_and_partition.nex")
    assert pset.name == "scheme"
    assert pset["part1"] == Partition(
        name="part1", ranges=(slice(0, 100),), model="HKY+G"
    )

    # file 05: [ ... ] comments sit between charset lines; stepped multi-range
    pset = parse_nexus_partitions("data/05_comprehensive_mixed.nex")
    assert pset.name == "full"
    assert len(pset) == 5
    assert pset["cds_pos12"] == Partition(
        name="cds_pos12",
        ranges=(slice(900, 1500, 3), slice(901, 1500, 3)),
        model="GTR+G",
        align_name="genes.phy",
    )
    assert pset["protein"].ranges == (slice(0, None),)

    # a lines iterable (not a path) with a {value} branch length
    lines = [
        "#nexus\n",
        "begin sets;\n",
        "    charset part1 = 1-100;\n",
        "    charpartition s = HKY+G{0.5}:part1;\n",
        "end;\n",
    ]
    pset = parse_nexus_partitions(lines)
    assert pset["part1"].model == "HKY+G"
    assert pset["part1"].tree_len == 0.5


_NEXUS_LOWER = """#nexus
begin data;
  dimensions ntax=2 nchar=4;
  format datatype=dna;
  matrix
    a  acgt
    b  tgca
  ;
end;
"""


def test_iter_nexus_align_records_uppercases():
    got = dict(iter_nexus_align_records(_NEXUS_LOWER.splitlines()))
    assert got == {"a": "ACGT", "b": "TGCA"}


def test_iter_nexus_align_records_custom_converter():
    records = iter_nexus_align_records(_NEXUS_LOWER.splitlines(), converter=bytes.lower)
    got = dict(records)
    assert got == {"a": b"acgt", "b": b"tgca"}
