from unittest import TestCase, main

from cogent3 import ASCII, DNA, make_aligned_seqs
from cogent3.core.annotation import Feature, Variable
# Complete version of manipulating sequence annotations
from cogent3.util.deserialise import deserialise_object


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class FeaturesTest(TestCase):
    """Tests of features in core"""

    def setUp(self):
        # A Sequence with a couple of exons on it.
        self.s = DNA.make_seq(
            "AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAA", name="Orig"
        )
        self.exon1 = self.s.add_annotation(Feature, "exon", "fred", [(10, 15)])
        self.exon2 = self.s.add_annotation(Feature, "exon", "trev", [(30, 40)])
        self.nested_feature = self.exon1.add_feature("repeat", "bob", [(2, 5)])

    def test_exon_extraction(self):
        """exon feature used to slice or directly access sequence"""
        # The corresponding sequence can be extracted either with
        # slice notation or by asking the feature to do it,
        # since the feature knows what sequence it belongs to.

        self.assertEqual(str(self.s[self.exon1]), "CCCCC")
        self.assertEqual(str(self.exon1.get_slice()), "CCCCC")

    def test_get_annotations_matching(self):
        """correctly identifies all features of a given type"""

        # Usually the only way to get a Feature object like exon1
        # is to ask the sequence for it. There is one method for querying
        # annotations by type and optionally by name:

        exons = self.s.get_annotations_matching("exon")
        self.assertEqual(
            str(exons), '[exon "fred" at [10:15]/48, exon "trev" at [30:40]/48]'
        )

    def test_get_nested_annotations_matching(self):
        """correctly identifies all features of a given type when nested annotations"""

        seq = DNA.make_seq("AAAAAAAAA", name="x")
        exon = seq.add_annotation(Feature, "exon", "fred", [(3, 8)])
        nested_exon = exon.add_annotation(Feature, "exon", "fred", [(3, 7)])
        exons = seq.get_annotations_matching("exon", extend_query=True)
        self.assertEqual(len(exons), 2)
        self.assertEqual(str(exons), '[exon "fred" at [3:8]/9, exon "fred" at [3:7]/5]')
        # tests multiple layers of nested annotations
        nested_exon.add_annotation(Feature, "exon", "fred", [(3, 6)])
        exons = seq.get_annotations_matching("exon", extend_query=True)
        self.assertEqual(len(exons), 3)
        self.assertEqual(
            str(exons),
            '[exon "fred" at [3:8]/9, exon "fred" at [3:7]/5, exon "fred" at [3:6]/4]',
        )
        # tests extend_query=False, and only get back the base exon
        exons = seq.get_annotations_matching("exon")
        self.assertEqual(len(exons), 1)
        self.assertEqual(str(exons), '[exon "fred" at [3:8]/9]')

    def test_get_annotations_matching2(self):
        """get_annotations_matching returns empty feature if no matches"""

        # If the sequence does not have a matching feature
        # you get back an empty list, and slicing the sequence
        # with that returns a sequence of length 0.

        dont_exist = self.s.get_annotations_matching("dont_exist")
        self.assertEqual(dont_exist, [])
        self.assertEqual(str(self.s[dont_exist]), "")

    def test_get_region_covering_all(self):
        """combines multiple features into one or their shadow"""

        # To construct a pseudo-feature covering (or excluding)
        # multiple features, use get_region_covering_all:

        exons = self.s.get_annotations_matching("exon")
        self.assertEqual(
            str(self.s.get_region_covering_all(exons)),
            'region "exon" at [10:15, 30:40]/48',
        )
        self.assertEqual(
            str(self.s.get_region_covering_all(exons).get_shadow()),
            'region "not exon" at [0:10, 15:30, 40:48]/48',
        )

        # eg: all the exon sequence:

        self.assertEqual(
            str(self.s.get_region_covering_all(exons).get_slice()), "CCCCCTTTTTAAAAA"
        )

        # or with slice notation:

        self.assertEqual(str(self.s[self.exon1, self.exon2]), "CCCCCTTTTTAAAAA")

    def test_slice_errors_from_merged(self):
        """no overlap in merged features"""

        # Though .get_region_covering_all also guarantees
        # no overlaps within the result, slicing does not:

        exons = self.s.get_annotations_matching("exon")
        self.assertEqual(
            str(self.s.get_region_covering_all(exons + exons)),
            'region "exon" at [10:15, 30:40]/48',
        )
        with self.assertRaises(ValueError):
            self.s[self.exon1, self.exon1, self.exon1, self.exon1, self.exon1]

        # You can use features, maps, slices or integers,
        # but non-monotonic slices are not allowed:

        with self.assertRaises(ValueError):
            self.s[15:20, 5:16]

    def test_feature_slicing(self):
        """features can be sliced"""

        # Features are themselves sliceable:

        self.assertEqual(str(self.exon1[0:3].get_slice()), "CCC")

        # When sequences are concatenated they keep their (non-overlapping) annotations:

        c = self.s[self.exon1[4:]] + self.s
        self.assertEqual(len(c), 49)

        answers = [
            'exon "fred" at [-4-, 0:1]/49',
            'exon "fred" at [11:16]/49',
            'exon "trev" at [31:41]/49',
        ]

        array = []

        for feat in c.annotations:
            array.append(str(feat))

        self.assertEqual(answers, array)

        # Since features know their parents you can't
        # use a feature from one sequence to slice another:

        with self.assertRaises(ValueError):
            c[self.exon1]

    def test_feature_attach_detach(self):
        """correctly associate, disassociate from seq"""

        # Features are generally attached to the thing they annotate,
        # but in those cases where a free-floating feature is created it can
        # ater be attached:

        exons = self.s.get_annotations_matching("exon")
        self.assertEqual(len(self.s.annotations), 2)
        region = self.s.get_region_covering_all(exons)
        self.assertEqual(len(self.s.annotations), 2)
        region.attach()
        self.assertEqual(len(self.s.annotations), 3)
        region.detach()
        self.assertEqual(len(self.s.annotations), 2)

    def test_feature_reverse(self):
        """reverse complement of features"""

        # When dealing with sequences that can be reverse complemented
        # (e.g. DnaSequence) features are **not** reversed.
        # Features are considered to have strand specific meaning
        # (.e.g CDS, exons) and so stay on their original strands.
        # We create a sequence with a CDS that spans multiple exons,
        # and show that after getting the reverse complement we have
        # exactly the same result from getting the CDS annotation.

        plus = DNA.make_seq("AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA", name="plus")
        plus_cds = plus.add_annotation(
            Feature, "CDS", "gene", [(2, 6), (10, 15), (25, 35)]
        )
        self.assertEqual(str(plus_cds.get_slice()), "GGGGCCCCCTTTTTTTTTT")
        minus = plus.rc()
        minus_cds = minus.get_annotations_matching("CDS")[0]
        self.assertEqual(str(minus_cds.get_slice()), "GGGGCCCCCTTTTTTTTTT")

    def test_feature_from_alignment(self):
        """seq features obtained from the alignment"""

        # Sequence features can be accessed via a containing Alignment:

        aln = make_aligned_seqs(
            data=[["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False
        )
        self.assertEqual(str(aln), ">x\n-AAAAAAAAA\n>y\nTTTT--TTTT\n")
        exon = aln.get_seq("x").add_annotation(Feature, "exon", "fred", [(3, 8)])
        aln_exons = aln.get_annotations_from_seq("x", "exon")
        aln_exons = aln.get_annotations_from_any_seq("exon")
        # But these will be returned as **alignment**
        # features with locations in alignment coordinates.

        self.assertEqual(str(exon), 'exon "fred" at [3:8]/9')
        self.assertEqual(str(aln_exons[0]), 'exon "fred" at [4:9]/10')
        self.assertEqual(str(aln_exons[0].get_slice()), ">x\nAAAAA\n>y\n--TTT\n")
        aln_exons[0].attach()
        self.assertEqual(len(aln.annotations), 1)

        # Similarly alignment features can be projected onto the aligned sequences,
        # where they may end up falling across gaps:

        exons = aln.get_projected_annotations("y", "exon")
        self.assertEqual(str(exons), '[exon "fred" at [-2-, 4:7]/8]')
        self.assertEqual(str(aln.get_seq("y")[exons[0].map.without_gaps()]), "TTT")

        # We copy the annotations from another sequence,

        aln = make_aligned_seqs(
            data=[["x", "-AAAAAAAAA"], ["y", "TTTT--CCCC"]], array_align=False
        )
        self.s = DNA.make_seq("AAAAAAAAA", name="x")
        exon = self.s.add_annotation(Feature, "exon", "fred", [(3, 8)])
        exon = aln.get_seq("x").copy_annotations(self.s)
        aln_exons = list(aln.get_annotations_from_seq("x", "exon"))
        self.assertEqual(str(aln_exons), '[exon "fred" at [4:9]/10]')

        # even if the name is different.

        exon = aln.get_seq("y").copy_annotations(self.s)
        aln_exons = list(aln.get_annotations_from_seq("y", "exon"))
        self.assertEqual(str(aln_exons), '[exon "fred" at [3:4, 6:10]/10]')
        self.assertEqual(str(aln[aln_exons]), ">x\nAAAAA\n>y\nTCCCC\n")

        # default for get_annotations_from_any_seq is return all features
        got = aln.get_annotations_from_any_seq()
        self.assertEqual(len(got), 2)

    def test_lost_spans(self):
        """features no longer included in an alignment represented by lost spans"""

        # If the feature lies outside the sequence being copied to, you get a
        # lost span

        aln = make_aligned_seqs(
            data=[["x", "-AAAA"], ["y", "TTTTT"]], array_align=False
        )
        seq = DNA.make_seq("CCCCCCCCCCCCCCCCCCCC", "x")
        seq.add_feature("exon", "A", [(5, 8)])
        aln.get_seq("x").copy_annotations(seq)
        copied = list(aln.get_annotations_from_seq("x", "exon"))
        self.assertEqual(str(copied), '[exon "A" at [5:5, -4-]/5]')
        self.assertEqual(str(copied[0].get_slice()), ">x\n----\n>y\n----\n")

    def test_seq_different_name_with_same_length(self):
        """copying features between sequences"""

        # You can copy to a sequence with a different name,
        # in a different alignment if the feature lies within the length

        aln = make_aligned_seqs(
            data=[["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False
        )
        seq = DNA.make_seq("CCCCCCCCCCCCCCCCCCCC", "x")
        seq.add_feature("exon", "A", [(5, 8)])
        aln.get_seq("y").copy_annotations(seq)
        copied = list(aln.get_annotations_from_seq("y", "exon"))
        self.assertEqual(str(copied), '[exon "A" at [7:10]/10]')

    def test_seq_shorter(self):
        """lost spans on shorter sequences"""

        # If the sequence is shorter, again you get a lost span.

        aln = make_aligned_seqs(
            data=[["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False
        )
        diff_len_seq = DNA.make_seq("CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "x")
        diff_len_seq.add_feature("repeat", "A", [(12, 14)])
        aln.get_seq("y").copy_annotations(diff_len_seq)
        copied = list(aln.get_annotations_from_seq("y", "repeat"))
        self.assertEqual(str(copied), '[repeat "A" at [10:10, -6-]/10]')

    def test_terminal_gaps(self):
        """features in cases of terminal gaps"""

        # We consider cases where there are terminal gaps.

        aln = make_aligned_seqs(
            data=[["x", "-AAAAAAAAA"], ["y", "------TTTT"]], array_align=False
        )
        aln.get_seq("x").add_feature("exon", "fred", [(3, 8)])
        aln_exons = list(aln.get_annotations_from_seq("x", "exon"))
        self.assertEqual(str(aln_exons), '[exon "fred" at [4:9]/10]')
        self.assertEqual(str(aln_exons[0].get_slice()), ">x\nAAAAA\n>y\n--TTT\n")
        aln = make_aligned_seqs(
            data=[["x", "-AAAAAAAAA"], ["y", "TTTT--T---"]], array_align=False
        )
        aln.get_seq("x").add_feature("exon", "fred", [(3, 8)])
        aln_exons = list(aln.get_annotations_from_seq("x", "exon"))
        self.assertEqual(str(aln_exons[0].get_slice()), ">x\nAAAAA\n>y\n--T--\n")

    def test_feature_residue(self):
        """seq features on alignment operate in sequence coordinates"""
        # In this case, only those residues included within the feature are
        # covered - note the omission of the T in y opposite the gap in x.

        aln = make_aligned_seqs(
            data=[["x", "C-CCCAAAAA"], ["y", "-T----TTTT"]],
            moltype=DNA,
            array_align=False,
        )
        self.assertEqual(str(aln), ">x\nC-CCCAAAAA\n>y\n-T----TTTT\n")
        exon = aln.get_seq("x").add_feature("exon", "ex1", [(0, 4)])
        self.assertEqual(str(exon), 'exon "ex1" at [0:4]/9')
        self.assertEqual(str(exon.get_slice()), "CCCC")
        aln_exons = list(aln.get_annotations_from_seq("x", "exon"))
        self.assertEqual(str(aln_exons), '[exon "ex1" at [0:1, 2:5]/10]')
        self.assertEqual(str(aln_exons[0].get_slice()), ">x\nCCCC\n>y\n----\n")

        # Feature.as_one_span(), is applied to the exon that
        # straddles the gap in x. The result is we preserve that feature.

        self.assertEqual(
            str(aln_exons[0].as_one_span().get_slice()), ">x\nC-CCC\n>y\n-T---\n"
        )

        # These properties also are consistently replicated with reverse
        # complemented sequences.

        aln_rc = aln.rc()
        rc_exons = list(aln_rc.get_annotations_from_any_seq("exon"))
        # not using as_one_span, so gap removed from x
        self.assertEqual(str(aln_rc[rc_exons]), ">x\nCCCC\n>y\n----\n")
        self.assertEqual(
            str(aln_rc[rc_exons[0].as_one_span()]), ">x\nC-CCC\n>y\n-T---\n"
        )

        # Features can provide their coordinates, useful for custom analyses.

        all_exons = aln.get_region_covering_all(aln_exons)
        coords = all_exons.get_coordinates()
        assert coords == [(0, 1), (2, 5)]

    def test_annotated_region_masks(self):
        """masking a sequence with specific features"""

        # Annotated regions can be masked (observed sequence characters
        # replaced by another), either through the sequence on which they
        # reside or by projection from the alignment. Note that mask_char must
        # be a valid character for the sequence MolType. Either the features
        # (multiple can be named), or their shadow, can be masked.

        # We create an alignment with a sequence that has two different annotation types.

        aln = make_aligned_seqs(
            data=[["x", "C-CCCAAAAAGGGAA"], ["y", "-T----TTTTG-GTT"]], array_align=False
        )
        self.assertEqual(str(aln), ">x\nC-CCCAAAAAGGGAA\n>y\n-T----TTTTG-GTT\n")
        exon = aln.get_seq("x").add_feature("exon", "norwegian", [(0, 4)])
        self.assertEqual(str(exon.get_slice()), "CCCC")
        repeat = aln.get_seq("x").add_feature("repeat", "blue", [(9, 12)])
        self.assertEqual(str(repeat.get_slice()), "GGG")
        repeat = aln.get_seq("y").add_feature("repeat", "frog", [(5, 7)])
        self.assertEqual(str(repeat.get_slice()), "GG")

        # Each sequence should correctly mask either the single feature,
        # it's shadow, or the multiple features, or shadow.

        self.assertEqual(
            str(aln.get_seq("x").with_masked_annotations("exon", mask_char="?")),
            "????AAAAAGGGAA",
        )
        self.assertEqual(
            str(
                aln.get_seq("x").with_masked_annotations(
                    "exon", mask_char="?", shadow=True
                )
            ),
            "CCCC??????????",
        )
        self.assertEqual(
            str(
                aln.get_seq("x").with_masked_annotations(
                    ["exon", "repeat"], mask_char="?"
                )
            ),
            "????AAAAA???AA",
        )
        self.assertEqual(
            str(
                aln.get_seq("x").with_masked_annotations(
                    ["exon", "repeat"], mask_char="?", shadow=True
                )
            ),
            "CCCC?????GGG??",
        )
        self.assertEqual(
            str(aln.get_seq("y").with_masked_annotations("exon", mask_char="?")),
            "TTTTTGGTT",
        )
        self.assertEqual(
            str(aln.get_seq("y").with_masked_annotations("repeat", mask_char="?")),
            "TTTTT??TT",
        )
        self.assertEqual(
            str(
                aln.get_seq("y").with_masked_annotations(
                    "repeat", mask_char="?", shadow=True
                )
            ),
            "?????GG??",
        )

        # The same methods can be applied to annotated Alignment's.

        self.assertEqual(
            str(aln.with_masked_annotations("exon", mask_char="?")),
            ">x\n?-???AAAAAGGGAA\n>y\n-T----TTTTG-GTT\n",
        )
        self.assertEqual(
            str(aln.with_masked_annotations("exon", mask_char="?", shadow=True)),
            ">x\nC-CCC??????????\n>y\n-?----?????-???\n",
        )
        self.assertEqual(
            str(aln.with_masked_annotations("repeat", mask_char="?")),
            ">x\nC-CCCAAAAA???AA\n>y\n-T----TTTT?-?TT\n",
        )
        self.assertEqual(
            str(aln.with_masked_annotations("repeat", mask_char="?", shadow=True)),
            ">x\n?-????????GGG??\n>y\n-?----????G-G??\n",
        )
        self.assertEqual(
            str(aln.with_masked_annotations(["repeat", "exon"], mask_char="?")),
            ">x\n?-???AAAAA???AA\n>y\n-T----TTTT?-?TT\n",
        )
        self.assertEqual(
            str(aln.with_masked_annotations(["repeat", "exon"], shadow=True)),
            ">x\nC-CCC?????GGG??\n>y\n-?----????G-G??\n",
        )

    def test_projected_to_base(self):
        """tests a given annotation is correctly projected on the base sequence"""

        seq = DNA.make_seq("AAAAAAAAATTTTTTTTT", name="x")
        layer_one = seq.add_feature("repeat", "frog", [(1, 17)])
        layer_two = layer_one.add_feature("repeat", "frog", [(2, 16)])
        got = layer_two._projected_to_base(seq)
        self.assertEqual(got.map.start, 3)
        self.assertEqual(got.map.end, 17)
        self.assertEqual(got.map.parent_length, len(seq))

        layer_three = layer_two.add_feature("repeat", "frog", [(5, 10)])
        got = layer_three._projected_to_base(seq)
        self.assertEqual(got.map.start, 8)
        self.assertEqual(got.map.end, 13)
        self.assertEqual(got.map.parent_length, len(seq))

        layer_four = layer_three.add_feature("repeat", "frog", [(0, 4)])
        layer_five = layer_four.add_feature("repeat", "frog", [(1, 2)])
        got = layer_five._projected_to_base(seq)
        self.assertEqual(got.map.start, 9)
        self.assertEqual(got.map.end, 10)
        self.assertEqual(got.map.parent_length, len(seq))

    def test_nested_annotated_region_masks(self):
        """masking a sequence with specific features when nested annotations"""

        aln = make_aligned_seqs(
            data=[["x", "C-GGCAAAAATTTAA"], ["y", "-T----TTTTG-GTT"]], array_align=False
        )
        gene = aln.get_seq("x").add_feature("gene", "norwegian", [(0, 4)])
        self.assertEqual(str(gene.get_slice()), "CGGC")
        gene.add_feature("repeat", "blue", [(1, 3)])
        # evaluate the sequence directly
        masked = str(
            aln.get_seq("x").with_masked_annotations(
                "repeat", mask_char="?", extend_query=True
            )
        )
        self.assertEqual(masked, "C??CAAAAATTTAA")

        exon = aln.get_seq("y").add_feature("repeat", "frog", [(1, 4)])
        self.assertEqual(str(exon.get_slice()), "TTT")
        # evaluate the sequence directly
        masked = str(
            aln.get_seq("y").with_masked_annotations(
                "repeat", mask_char="?", extend_query=True
            )
        )
        self.assertEqual(masked, "T???TGGTT")

        masked = aln.with_masked_annotations("gene", mask_char="?")
        got = masked.to_dict()
        self.assertEqual(got["x"], "?-???AAAAATTTAA")
        self.assertEqual(got["y"], "-T----TTTTG-GTT")

    def test_annotated_separately_equivalence(self):
        """allow defining features as a series or individually"""

        # It shouldn't matter whether annotated coordinates are entered
        # separately, or as a series.

        data = [["human", "CGAAACGTTT"], ["mouse", "CTAAACGTCG"]]
        as_series = make_aligned_seqs(data=data, array_align=False)
        as_items = make_aligned_seqs(data=data, array_align=False)

        # We add annotations to the sequences as a series.

        self.assertEqual(
            str(
                as_series.get_seq("human").add_feature(
                    "cpgsite", "cpg", [(0, 2), (5, 7)]
                )
            ),
            'cpgsite "cpg" at [0:2, 5:7]/10',
        )
        self.assertEqual(
            str(
                as_series.get_seq("mouse").add_feature(
                    "cpgsite", "cpg", [(5, 7), (8, 10)]
                )
            ),
            'cpgsite "cpg" at [5:7, 8:10]/10',
        )

        # We add the annotations to the sequences one segment at a time.

        self.assertEqual(
            str(as_items.get_seq("human").add_feature("cpgsite", "cpg", [(0, 2)])),
            'cpgsite "cpg" at [0:2]/10',
        )
        self.assertEqual(
            str(as_items.get_seq("human").add_feature("cpgsite", "cpg", [(5, 7)])),
            'cpgsite "cpg" at [5:7]/10',
        )
        self.assertEqual(
            str(as_items.get_seq("mouse").add_feature("cpgsite", "cpg", [(5, 7)])),
            'cpgsite "cpg" at [5:7]/10',
        )
        self.assertEqual(
            str(as_items.get_seq("mouse").add_feature("cpgsite", "cpg", [(8, 10)])),
            'cpgsite "cpg" at [8:10]/10',
        )

    def test_constructor_equivalence(self):
        """ """

        # These different constructions should generate the same output.
        data = [["human", "CGAAACGTTT"], ["mouse", "CTAAACGTCG"]]
        as_series = make_aligned_seqs(data=data, array_align=False)
        as_items = make_aligned_seqs(data=data, array_align=False)

        serial = as_series.with_masked_annotations(["cpgsite"])
        itemwise = as_items.with_masked_annotations(["cpgsite"])
        self.assertEqual(str(serial), str(itemwise))

        # Annotations should be correctly masked,
        # whether the sequence has been reverse complemented or not.
        # We use the plus/minus strand CDS containing sequences created above.
        plus = DNA.make_seq("AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA", name="plus")
        _ = plus.add_annotation(Feature, "CDS", "gene", [(2, 6), (10, 15), (25, 35)])
        minus = plus.rc()
        self.assertEqual(
            str(plus.with_masked_annotations("CDS")),
            "AA????AAAA?????AAAAAAAAAA??????????AAA",
        )
        self.assertEqual(
            str(minus.with_masked_annotations("CDS")),
            "TTT??????????TTTTTTTTTT?????TTTT????TT",
        )

    def test_roundtrip_json(self):
        """features can roundtrip from json"""

        seq = DNA.make_seq("AAAAATATTATTGGGT")
        seq.add_annotation(Feature, "exon", "myname", [(0, 5)])
        got = seq.to_json()
        new = deserialise_object(got)
        feat = new.get_annotations_matching("exon")[0]
        self.assertEqual(str(feat.get_slice()), "AAAAA")

        # now with a list span
        seq = seq[3:]
        feat = seq.get_annotations_matching("exon")[0]
        got = seq.to_json()
        new = deserialise_object(got)
        feat = new.get_annotations_matching("exon")[0]
        self.assertEqual(str(feat.get_slice(complete=False)), "AA")

    def test_roundtripped_alignment(self):
        """Alignment with annotations roundtrips correctly"""
        # annotations just on member sequences
        aln = make_aligned_seqs(
            data=[["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False
        )
        _ = aln.get_seq("x").add_annotation(Feature, "exon", "fred", [(3, 8)])
        seq_exon = list(aln.get_annotations_from_seq("x", "exon"))[0]
        expect = seq_exon.get_slice()

        json = aln.to_json()
        new = deserialise_object(json)
        got_exons = list(new.get_annotations_from_seq("x", "exon"))[0]
        self.assertEqual(got_exons.get_slice().to_dict(), expect.to_dict())

        # annotations just on alignment
        aln = make_aligned_seqs(
            data=[["x", "-AAAAAGGGG"], ["y", "TTTT--CCCC"]], array_align=False
        )
        f = aln.add_annotation(Feature, "generic", "no name", [(1, 4), (6, 10)])
        expect = f.get_slice().to_dict()
        json = aln.to_json()
        new = deserialise_object(json)
        got = list(new.get_annotations_matching("generic"))[0]
        self.assertEqual(got.get_slice().to_dict(), expect)

        # annotations on both alignment and sequence
        aln = make_aligned_seqs(
            data=[["x", "-AAAAAGGGG"], ["y", "TTTT--CCCC"]], array_align=False
        )
        f = aln.add_annotation(Feature, "generic", "no name", [(1, 4), (6, 10)])
        _ = aln.get_seq("x").add_annotation(Feature, "exon", "1", [(3, 8)])
        json = aln.to_json()
        new = deserialise_object(json)
        ## get back the exon
        seq_exon = list(aln.get_annotations_from_seq("x", "exon"))[0]
        expect = seq_exon.get_slice().to_dict()
        got_exons = list(new.get_annotations_from_seq("x", "exon"))[0]
        self.assertEqual(got_exons.get_slice().to_dict(), expect)
        ## get back the generic
        expect = f.get_slice().to_dict()
        got = list(new.get_annotations_matching("generic"))[0]
        self.assertEqual(got.get_slice().to_dict(), expect)

        # check masking of seq features still works
        new = new.with_masked_annotations("exon", mask_char="?")
        self.assertEqual(new[4:9].to_dict(), dict(x="?????", y="--CCC"))

    def test_roundtripped_alignment2(self):
        """Sliced Alignment with annotations roundtrips correctly"""
        # annotations just on member sequences
        aln = make_aligned_seqs(
            data=[["x", "-AAAGGGGGAACCCT"], ["y", "TTTT--TTTTAGGGA"]], array_align=False
        )
        aln.get_seq("x").add_annotation(Feature, "exon", "E1", [(3, 8)])
        aln.get_seq("x").add_annotation(Feature, "exon", "E2", [(10, 13)])
        # at the alignment level
        sub_aln = aln[:-3]
        s = sub_aln.named_seqs["x"]
        s.data.get_annotations_matching("exon", "E2")[0]
        d = s.data[:11]
        json = s.to_json()
        new = deserialise_object(json)
        gf1, gf2 = list(new.data.get_annotations_matching("exon"))
        self.assertEqual(str(gf1.get_slice()), "GGGGG")
        self.assertEqual(str(gf2.get_slice()), "C")
        # the sliced alignment
        json = sub_aln.to_json()
        got = deserialise_object(json)
        x = got.named_seqs["x"]
        self.assertEqual(str(x.data.annotations[0].get_slice()), "GGGGG")
        self.assertEqual(str(x.data.annotations[1].get_slice()), "C")

    def test_roundtrip_rc_annotated_align(self):
        """should work for an alignment that has been reverse complemented"""
        # the key that exposed the bug was a gap in the middle of the sequence
        aln = make_aligned_seqs(
            data=[["x", "-AAAGGGGGAAC-CT"], ["y", "TTTT--TTTTAGGGA"]],
            array_align=False,
            moltype="dna",
        )
        aln.get_seq("x").add_annotation(Feature, "exon", "E1", [(3, 8)])
        aln.get_seq("x").add_annotation(Feature, "exon", "E2", [(10, 13)])

        raln = aln.rc()
        json = raln.to_json()
        got = deserialise_object(json)
        self.assertEqual(got.to_dict(), raln.to_dict())
        orig_annots = {
            a.name: a.get_slice() for a in raln.get_annotations_from_any_seq()
        }
        got_annots = {a.name: a.get_slice() for a in got.get_annotations_from_any_seq()}
        self.assertEqual(got_annots, orig_annots)

    def test_roundtrip_variable(self):
        """should recover the Variable feature type"""
        seq = DNA.make_seq("AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA", name="plus")
        xx_y = [[[2, 6], 2.4], [[10, 15], 5.1], [[25, 35], 1.3]]
        y_valued = seq.add_annotation(Variable, "SNP", "freq", xx_y)
        json = seq.to_json()
        new = deserialise_object(json)
        got = list(new.get_annotations_matching("SNP"))[0]
        # annoyingly, comes back as list of lists
        self.assertEqual(got.xxy_list, [[list(xx), y] for xx, y in y_valued.xxy_list])

    def test_nested_get_slice(self):
        """check the get_slice method works on nested annotations"""
        self.assertEqual(self.nested_feature.get_slice(), "CCC")

    def test_nested_to_rich_dict(self):
        """check the to_rich_dict method works with nested annotations"""
        self.assertEqual(
            self.exon1.to_rich_dict()["annotations"][0],
            self.nested_feature.to_rich_dict(),
        )

    def test_nested_deserialise_annotation(self):
        """nested annotations can be deserialised"""
        got = self.s.to_json()
        new = deserialise_object(got)
        new_exon1 = new.annotations[0]
        new_nested_feature = new_exon1.annotations[0]
        self.assertEqual(
            new_nested_feature.to_rich_dict(), self.nested_feature.to_rich_dict()
        )

    def test_annotate_matches_to(self):
        """annotate_matches_to attaches annotations correctly to a Sequence"""
        seq = DNA.make_seq("TTCCACTTCCGCTT", name="x")
        pattern = "CCRC"
        annot = seq.annotate_matches_to(
            pattern=pattern, annot_type="domain", name="fred", allow_multiple=True
        )
        self.assertEqual([a.get_slice() for a in annot], ["CCAC", "CCGC"])
        annot = seq.annotate_matches_to(
            pattern=pattern, annot_type="domain", name="fred", allow_multiple=False
        )
        self.assertEqual(len(annot), 1)
        fred = annot[0].get_slice()
        self.assertEqual(str(fred), "CCAC")
        # For Sequence objects of a non-IUPAC MolType, annotate_matches_to
        # should return an empty annotation.
        seq = ASCII.make_seq(seq="TTCCACTTCCGCTT")
        annot = seq.annotate_matches_to(
            pattern=pattern, annot_type="domain", name="fred", allow_multiple=False
        )
        self.assertEqual(annot, [])


if __name__ == "__main__":
    main()
