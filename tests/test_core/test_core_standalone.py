#!/usr/bin/env python
import json
import os
import pathlib
import tempfile
import unittest

from tempfile import TemporaryDirectory

from cogent3 import (
    DNA,
    PROTEIN,
    RNA,
    get_format_suffixes,
    load_aligned_seqs,
    load_seq,
    load_unaligned_seqs,
    make_aligned_seqs,
    make_seq,
    make_unaligned_seqs,
)
from cogent3.app.data_store import make_record_for_json
from cogent3.core.alignment import (
    Alignment,
    ArrayAlignment,
    SequenceCollection,
)
from cogent3.core.alphabet import AlphabetError
from cogent3.parse.record import FileFormatError


__author__ = "Peter Maxwell, Gavin Huttley and Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

base_path = os.path.dirname(os.path.dirname(__file__))
data_path = os.path.join(base_path, "data")


class TestConstructorFunctions(unittest.TestCase):
    def test_make_seq(self):
        """test constructor utility function"""
        _seq = "ACGGT"
        seq = make_seq(_seq)
        self.assertEqual(seq.moltype.label, "text")
        seq = make_seq(_seq, moltype="dna")
        self.assertEqual(seq.moltype.label, "dna")
        self.assertEqual(str(seq), _seq)

    def test_make_unaligned_seqs(self):
        """test SequenceCollection constructor utility function"""
        data = {"a": "AGGTT", "b": "AG"}
        got = make_unaligned_seqs(data)
        self.assertIsInstance(got, SequenceCollection)
        self.assertEqual(got.to_dict(), data)
        self.assertEqual(got.info["source"], "unknown")
        # moltype arg works
        got = make_unaligned_seqs(data, moltype="dna")
        self.assertEqual(got.moltype.label, "dna")
        # info works
        got = make_unaligned_seqs(data, info=dict(a=2))
        self.assertEqual(got.info["a"], 2)
        with self.assertRaises(AssertionError):
            _ = make_unaligned_seqs(data, info=2)

        # source works
        for src in ("somewhere", pathlib.Path("somewhere")):
            got = make_unaligned_seqs(data, source=src)
            self.assertEqual(got.info["source"], str(src))

    def test_make_aligned_seqs(self):
        """test Alignment/ArrayAlignment constructor utility function"""
        data = {"a": "AGGTT", "b": "AGAA-"}
        got = make_aligned_seqs(data)
        self.assertIsInstance(got, ArrayAlignment)
        self.assertEqual(got.to_dict(), data)
        self.assertEqual(got.info["source"], "unknown")
        # moltype arg works
        got = make_aligned_seqs(data, moltype="dna")
        self.assertEqual(got.moltype.label, "dna")
        # info works
        got = make_aligned_seqs(data, info=dict(a=2))
        self.assertEqual(got.info["a"], 2)
        with self.assertRaises(AssertionError):
            _ = make_unaligned_seqs(data, info=2)

        # source works
        for src in ("somewhere", pathlib.Path("somewhere")):
            got = make_aligned_seqs(data, source=src)
            self.assertEqual(got.info["source"], str(src))

        # array_align works
        got = make_aligned_seqs(data, array_align=False)
        self.assertIsInstance(got, Alignment)
        self.assertEqual(got.to_dict(), data)
        self.assertEqual(got.info["source"], "unknown")

    def test_load_seq(self):
        """load single sequence"""
        from cogent3 import Sequence

        paths = (
            "c_elegans_WS199_dna_shortened.fasta",
            "annotated_seq.gb",
            "brca1_5.250.paml",
        )
        seq_names = ("I", "AE017341", "NineBande")
        data_dir = pathlib.Path(data_path)
        for i, path in enumerate(paths):
            got = load_seq(data_dir / path)
            assert isinstance(got, Sequence)
            assert got.info.source == str(data_dir / path)
            assert got.name == seq_names[i]

        # try json
        seq = got
        with TemporaryDirectory(dir=".") as dirname:
            outpath = pathlib.Path(dirname) / "seq.json"
            outpath.write_text(seq.to_json())
            got = load_seq(outpath)
            assert str(got) == str(seq)
            assert got.name == seq.name

    def test_load_unaligned_seqs(self):
        """test loading unaligned from file"""
        path = os.path.join(data_path, "brca1_5.paml")
        got = load_unaligned_seqs(path)
        self.assertIsInstance(got, SequenceCollection)
        self.assertTrue("Human" in got.to_dict())
        self.assertEqual(got.info["source"], path)

    def test_load_unaligned_seqs_no_format(self):
        """test loading unaligned from file"""
        with self.assertRaises(ValueError):
            load_unaligned_seqs("somepath")

    def test_load_aligned_seqs(self):
        """test loading aligned from file"""
        path = os.path.join(data_path, "brca1_5.paml")
        got = load_aligned_seqs(path)
        self.assertIsInstance(got, ArrayAlignment)
        self.assertTrue("Human" in got.to_dict())
        self.assertEqual(got.info["source"], path)
        got = load_aligned_seqs(path, moltype="dna")
        self.assertEqual(got.moltype.label, "dna")
        got = load_aligned_seqs(path, moltype="dna", array_align=False)
        self.assertEqual(got.moltype.label, "dna")
        self.assertIsInstance(got, Alignment)

    def test_load_aligned_seqs_no_format(self):
        """test loading unaligned from file"""
        with self.assertRaises(ValueError):
            load_aligned_seqs("somepath")

    def test_load_unaligned_seqs_from_json(self):
        """test loading an unaligned object from json file"""
        with TemporaryDirectory(dir=".") as dirname:
            json_path = os.path.join(dirname, "unaligned.json")
            path = os.path.join(data_path, "brca1_5.paml")
            unaligned = load_unaligned_seqs(path)
            unaligned.write(json_path)

            got = load_unaligned_seqs(json_path)
            self.assertIsInstance(got, SequenceCollection)
            self.assertEqual(got.to_dict(), unaligned.to_dict())
            self.assertEqual(got.info["source"], path)
            # tests json generated by make_record_for_json
            uncompleted_record = make_record_for_json("delme", got, False)
            uncompleted_record_path = os.path.join(dirname, "uncompleted_record.json")
            completed_record = make_record_for_json("delme", got, True)
            completed_record_path = os.path.join(dirname, "completed_record.json")
            with open(uncompleted_record_path, "w") as out:
                out.write(json.dumps(uncompleted_record))
            with open(completed_record_path, "w") as out:
                out.write(json.dumps(completed_record))
            # tests when provided record json file is uncompleted
            with self.assertRaises(TypeError):
                load_unaligned_seqs(uncompleted_record_path)
            # tests when provided record json is completed
            got = load_unaligned_seqs(completed_record_path)
            self.assertIsInstance(got, SequenceCollection)
            self.assertEqual(got.to_dict(), unaligned.to_dict())
            self.assertEqual(got.info["source"], path)

    def test_load_aligned_seqs_from_json(self):
        """tests loading an aligned object from json file"""
        with TemporaryDirectory(dir=".") as dirname:
            path = os.path.join(data_path, "brca1_5.paml")
            alignment = load_aligned_seqs(path, array_align=False, moltype="dna")
            alignment_json_path = os.path.join(dirname, "alignment.json")
            alignment.write(alignment_json_path)
            array_alignment = load_aligned_seqs(path, moltype="dna")
            array_alignment_json_path = os.path.join(dirname, "array_alignment.json")
            array_alignment.write(array_alignment_json_path)
            # tests case Alignment
            got = load_aligned_seqs(alignment_json_path)
            self.assertIsInstance(got, Alignment)
            self.assertEqual(got.moltype.label, "dna")
            self.assertEqual(got.to_dict(), alignment.to_dict())
            self.assertEqual(got.info["source"], path)
            # tests case ArrayAlignment
            got = load_aligned_seqs(array_alignment_json_path)
            self.assertIsInstance(got, ArrayAlignment)
            self.assertEqual(got.moltype.label, "dna")
            self.assertEqual(got.to_dict(), array_alignment.to_dict())
            self.assertEqual(got.info["source"], path)

            # tests json generated by make_record_for_json
            uncompleted_record = make_record_for_json("delme", got, False)
            completed_record = make_record_for_json("delme", got, True)
            uncompleted_record_path = os.path.join(dirname, "uncompleted_record.json")
            completed_record_path = os.path.join(dirname, "completed_record.json")
            with open(uncompleted_record_path, "w") as out:
                out.write(json.dumps(uncompleted_record))
            with open(completed_record_path, "w") as out:
                out.write(json.dumps(completed_record))
            # tests when provided record json file is uncompleted
            with self.assertRaises(TypeError):
                load_unaligned_seqs(uncompleted_record_path)
            # tests when provided record json is completed
            got = load_aligned_seqs(completed_record_path)
            self.assertIsInstance(got, ArrayAlignment)
            self.assertEqual(got.to_dict(), array_alignment.to_dict())
            self.assertEqual(got.info["source"], path)
            # tests wrong input json file
            json_path = os.path.join(dirname, "unaligned.json")
            path = os.path.join(data_path, "brca1_5.paml")
            unaligned = load_unaligned_seqs(path)
            unaligned.write(json_path)
            with self.assertRaises(TypeError):
                load_aligned_seqs(json_path)


class ReadingWritingFileFormats(unittest.TestCase):
    """Testing ability to read file formats."""

    def _loadfromfile(self, filename, test_write=True, **kw):
        filename = os.path.join(data_path, filename)
        aln = load_aligned_seqs(filename, **kw)
        if test_write:
            suffix, cmpr = get_format_suffixes(filename)
            cmpr = "" if not cmpr else f".{cmpr}"
            fn = tempfile.mktemp(suffix="." + suffix + cmpr)
            aln.write(filename=fn)
            os.remove(fn)
            # now use pathlib
            fn = pathlib.Path(fn)
            aln.write(filename=fn)
            fn.unlink()

    def test_write_unknown_raises(self):
        """writing unknown format raises FileFormatError"""
        filename = os.path.join(data_path, "primates_brca1.fasta")
        aln = load_aligned_seqs(filename)
        self.assertRaises(FileFormatError, aln.write, filename="blah")
        self.assertRaises(FileFormatError, aln.write, filename="blah.txt")
        self.assertRaises(
            FileFormatError, aln.write, filename="blah.fasta", format="noway"
        )

    def test_fasta(self):
        self._loadfromfile("formattest.fasta")

    def test_fasta_gzipped(self):
        """correctly load from gzipped"""
        self._loadfromfile("formattest.fasta.gz")

    def test_fasta_bzipped(self):
        """correctly load from bzipped"""
        self._loadfromfile("formattest.fasta.bz2")

    def test_phylipsequential(self):
        self._loadfromfile("formattest.phylip")

    def test_clustal(self):
        self._loadfromfile("formattest.aln", test_write=False)

    def test_phylip_interleaved(self):
        self._loadfromfile(
            "interleaved.phylip", test_write=False, parser_kw=dict(interleaved=True)
        )

    def test_paml(self):
        self._loadfromfile("formattest.paml")

    def test_gde(self):
        self._loadfromfile("formattest.gde")

    def test_msf(self):
        self._loadfromfile("formattest.msf", test_write=False)


class AlignmentTestMethods(unittest.TestCase):
    """Testing Alignment methods"""

    def setUp(self):
        self.alignment = load_aligned_seqs(
            filename=os.path.join(data_path, "brca1_5.paml")
        )

    def test_picklability(self):
        """Pickle an alignment containing an annotated sequence"""
        # This depends on alignments, sequences, features, maps and spans
        # Doesn't test round trip result is correct, which should possibly
        # be done for maps/spans, but seqs/alignments are just simple
        # python classes without __getstate__ etc.
        import pickle as pickle

        seq1 = DNA.make_seq("aagaagaagaccccca")
        seq2 = DNA.make_seq("aagaagaagaccccct")
        seq2.add_feature("exon", "fred", [(10, 15)])
        aln = make_aligned_seqs(data={"a": seq1, "b": seq2})
        # TODO the ability to pickle/unpickle depends on the protocol
        # in Py3 for reasons that are not clear. This needs to be looked
        # more closely
        dmp = pickle.dumps(aln, protocol=1)
        pickle.loads(dmp)

    def test_empty_seq(self):
        """test creation of an alignment from scratch, with one sequence pure gap"""
        new_seqs = {"seq1": "ATGATG", "seq2": "------"}
        align = make_aligned_seqs(moltype=DNA, data=new_seqs)
        assert len(align) == 6, align

    def test_num_seqs(self):
        self.assertEqual(self.alignment.num_seqs, 5)

    def test_numberseqs(self):
        """testing the number of sequences"""
        assert len(self.alignment.names) == 5

    def test_alignlsength(self):
        """testing the alignment length"""
        assert len(self.alignment) == 60

    def test_init_from_strings(self):
        """testing constructing an alignment from a dictionary of strings"""
        new_seqs = {"seq1": "ACGTACGT", "seq2": "ACGTACGT", "seq3": "ACGTACGT"}
        make_aligned_seqs(data=new_seqs)

    def test_get_sub_alignment(self):
        """test slicing otus, and return of new alignment"""
        subset = ["DogFaced", "Human", "HowlerMon", "Mouse"]
        subset.sort()
        sub_align = self.alignment.take_seqs(subset)
        new = sub_align.names
        new.sort()
        assert new == subset, f"included subset didn't work {new}, {subset}"

        # testing exclusion of one
        to_exclude = ["NineBande"]
        sub_align = self.alignment.take_seqs(to_exclude, negate=True)
        new = sub_align.names
        new.sort()
        assert new == subset, f"excluded subset didn't work {new}, {subset}"

        # testing exclusion of two
        subset = ["DogFaced", "HowlerMon", "NineBande"]
        subset.sort()
        to_exclude = ["Human", "Mouse"]
        sub_align = self.alignment.take_seqs(to_exclude, negate=True)
        new = sub_align.names
        new.sort()
        assert new == subset, f"excluded subset didn't work {new}, {subset}"

    def test_slice_align(self):
        """test slicing of sequences"""
        alignment = make_aligned_seqs(
            data={"seq1": "ACGTACGT", "seq2": "ACGTACGT", "seq3": "ACGTACGT"}
        )
        sub_align = alignment[2:5]
        self.assertEqual(len(sub_align), 3)
        self.assertEqual(len(sub_align.names), 3)
        self.assertEqual(
            sub_align.to_dict(), {"seq1": "GTA", "seq2": "GTA", "seq3": "GTA"}
        )

        sub_align = alignment[5:20]
        self.assertEqual(len(sub_align), 3)
        self.assertEqual(len(sub_align.names), 3)
        self.assertEqual(
            sub_align.to_dict(), {"seq1": "CGT", "seq2": "CGT", "seq3": "CGT"}
        )

        sub_align = alignment[2]
        self.assertEqual(len(sub_align), 1)
        self.assertEqual(sub_align.to_dict(), {"seq1": "G", "seq2": "G", "seq3": "G"})

        sub_align = alignment[0]
        self.assertEqual(len(sub_align), 1)
        self.assertEqual(sub_align.to_dict(), {"seq1": "A", "seq2": "A", "seq3": "A"})

        sub_align = alignment[7]
        self.assertEqual(len(sub_align), 1)
        self.assertEqual(sub_align.to_dict(), {"seq1": "T", "seq2": "T", "seq3": "T"})

    def test_sliding_windows(self):
        """test slicing of sequences"""
        alignment = make_aligned_seqs(
            data={"seq1": "ACGTACGT", "seq2": "ACGTACGT", "seq3": "ACGTACGT"}
        )
        result = [bit for bit in alignment.sliding_windows(5, 2)]
        self.assertEqual(
            result[0].to_dict(), {"seq3": "ACGTA", "seq2": "ACGTA", "seq1": "ACGTA"}
        )
        self.assertEqual(
            result[1].to_dict(), {"seq3": "GTACG", "seq2": "GTACG", "seq1": "GTACG"}
        )

        # specify a starting window position
        result = []
        for bit in alignment.sliding_windows(5, 2, start=1):
            result += [bit]
        self.assertEqual(
            result[0].to_dict(), {"seq3": "CGTAC", "seq2": "CGTAC", "seq1": "CGTAC"}
        )
        self.assertEqual(
            result[1].to_dict(), {"seq3": "TACGT", "seq2": "TACGT", "seq1": "TACGT"}
        )

        # specify a ending window position
        result = []
        for bit in alignment.sliding_windows(5, 1, start=1, end=3):
            result += [bit]
        self.assertEqual(
            result[0].to_dict(), {"seq3": "CGTAC", "seq2": "CGTAC", "seq1": "CGTAC"}
        )
        self.assertEqual(
            result[1].to_dict(), {"seq3": "GTACG", "seq2": "GTACG", "seq1": "GTACG"}
        )

        # start conditions < window-size from end don't return a window
        # specify a ending window position
        result = []
        for bit in alignment.sliding_windows(5, 1, start=5):
            result += [bit]
        self.assertEqual(result, [])

        result = []
        for bit in alignment.sliding_windows(5, 1):
            result += [bit]
        self.assertEqual(
            result[0].to_dict(), {"seq3": "ACGTA", "seq2": "ACGTA", "seq1": "ACGTA"}
        )
        self.assertEqual(
            result[1].to_dict(), {"seq3": "CGTAC", "seq2": "CGTAC", "seq1": "CGTAC"}
        )
        self.assertEqual(
            result[2].to_dict(), {"seq3": "GTACG", "seq2": "GTACG", "seq1": "GTACG"}
        )
        self.assertEqual(
            result[3].to_dict(), {"seq3": "TACGT", "seq2": "TACGT", "seq1": "TACGT"}
        )

    def test_omit_gap_pos1(self):
        """test removal of redundant gaps (all entries in alignment column are gaps)"""
        alignment = make_aligned_seqs(
            data={
                "seq1": "--ACGT--GT---",
                "seq2": "--ACGTA-GT---",
                "seq3": "--ACGTA-GT---",
            }
        )
        align_dict = alignment.omit_gap_pos().to_dict()
        self.assertEqual(
            align_dict, {"seq1": "ACGT-GT", "seq2": "ACGTAGT", "seq3": "ACGTAGT"}
        )

    def test_omit_gap_pos2(self):
        """test removal of all gaps (any entries in alignment column are gaps)"""
        alignment = make_aligned_seqs(
            data={
                "seq1": "--ACGT--GT---",
                "seq2": "--ACGTA-GT---",
                "seq3": "--ACGTA-GT---",
            }
        )
        align_dict = alignment.omit_gap_pos(allowed_gap_frac=0).to_dict()
        self.assertEqual(
            align_dict, {"seq1": "ACGTGT", "seq2": "ACGTGT", "seq3": "ACGTGT"}
        )

        alignment = make_aligned_seqs(
            data={"seq1": "ACGT", "seq2": "----", "seq3": "----"}
        )
        result = alignment.omit_gap_pos(allowed_gap_frac=0)
        self.assertEqual(result, None)

    def test_degap(self):
        """test stripping gaps from collections and alignments"""
        aln = make_aligned_seqs(
            data={
                "seq1": "--ACGT--GT---",
                "seq2": "--ACGTA-GT---",
                "seq3": "--ACGTA-GT---",
            }
        )
        observed = aln.degap()
        expect = {"seq1": "ACGTGT", "seq2": "ACGTAGT", "seq3": "ACGTAGT"}
        self.assertEqual(observed.to_dict(), expect)
        collection = make_unaligned_seqs(
            data={
                "seq1": "--ACGT--GT---",
                "seq2": "--ACGTA-GT---",
                "seq3": "--ACGTA-GT---",
            },
            moltype=DNA,
        )
        observed = collection.degap()
        self.assertEqual(observed.to_dict(), expect)
        self.assertEqual(observed.moltype, DNA)

    def test_DnaRna_interconversion(self):
        """test interconversion between Rna and Dna by SequenceCollection and
        Alignment"""
        dna = {
            "seq1": "--ACGT--GT---",
            "seq2": "--ACGTA-GT---",
            "seq3": "--ACGTA-GT---",
        }
        rna = {
            "seq1": "--ACGU--GU---",
            "seq2": "--ACGUA-GU---",
            "seq3": "--ACGUA-GU---",
        }
        collect_Dna = make_unaligned_seqs(data=dna, moltype=DNA)
        collect_Rna = make_unaligned_seqs(data=rna, moltype=RNA)
        self.assertEqual(collect_Rna.to_dna().to_dict(), dna)
        self.assertEqual(collect_Dna.to_rna().to_dict(), rna)

        aln_Dna = make_aligned_seqs(data=dna, moltype=DNA)
        aln_Rna = make_aligned_seqs(data=rna, moltype=RNA)
        rna_from_dna = aln_Dna.to_rna()
        dna_from_rna = aln_Rna.to_dna()
        self.assertEqual(rna_from_dna.to_dict(), rna)
        self.assertEqual(dna_from_rna.to_dict(), dna)

    def test_reverse_complement(self):
        """test reverse complementing of Alignments and SequenceCollection."""
        dna = {
            "seq1": "--ACGT--GT---",
            "seq2": "TTACGTA-GT---",
            "seq3": "--ACGTA-GCC--",
        }
        dna_rc = {
            "seq1": "---AC--ACGT--",
            "seq2": "---AC-TACGTAA",
            "seq3": "--GGC-TACGT--",
        }
        # alignment with gaps
        aln = make_aligned_seqs(data=dna, moltype=DNA)
        aln_rc = aln.rc()
        self.assertEqual(aln_rc.to_dict(), dna_rc)
        # check collection, with gaps
        coll = make_unaligned_seqs(data=dna, moltype=DNA)
        coll_rc = coll.rc()
        self.assertEqual(coll_rc.to_dict(), dna_rc)
        self.assertEqual(coll_rc.to_dict(), coll.reverse_complement().to_dict())
        # collection with no gaps
        dna = {"seq1": "ACGTGT", "seq2": "TTACGTAGT", "seq3": "ACGTAGCC"}
        dna_rc = {"seq1": "ACACGT", "seq2": "ACTACGTAA", "seq3": "GGCTACGT"}
        coll = make_unaligned_seqs(data=dna, moltype=DNA)
        coll_rc = coll.rc()
        self.assertEqual(coll_rc.to_dict(), dna_rc)

    def test_reverse_complement_info(self):
        """reverse_complement should preserve info attribute"""
        dna = {
            "seq1": "--ACGT--GT---",
            "seq2": "TTACGTA-GT---",
            "seq3": "--ACGTA-GCC--",
        }
        # alignment with gaps
        aln = ArrayAlignment(data=dna, moltype=DNA, info={"key": "value"})
        aln_rc = aln.rc()
        self.assertEqual(aln_rc.info["key"], "value")
        # check collection, with gaps
        coll = SequenceCollection(data=dna, moltype=DNA, info={"key": "value"})
        coll_rc = coll.rc()
        self.assertEqual(coll_rc.info["key"], "value")

    def test_reverse_complement_with_ambig(self):
        """correctly reverse complement with ambiguous bases"""
        n = make_aligned_seqs(data=[["x", "?-???AA"], ["y", "-T----T"]], moltype=DNA)
        rc = n.rc()
        self.assertEqual(rc.to_dict(), {"x": "TT???-?", "y": "A----A-"})

    def test_getasdict(self):
        """getting the alignment as a dictionary"""
        seqs = {"seq1": "ACGT--GT", "seq2": "ACGTACGT", "seq3": "ACGTACGT"}
        alignment = make_aligned_seqs(data=seqs)
        align_dict = alignment.to_dict()
        self.assertEqual(align_dict, seqs)

    def test_alignadd(self):
        """testing adding one alignment to another."""
        align1 = make_aligned_seqs(data={"a": "AAAA", "b": "TTTT", "c": "CCCC"})
        align2 = make_aligned_seqs(data={"a": "GGGG", "b": "----", "c": "NNNN"})
        align = align1 + align2
        concatdict = align.to_dict()
        self.assertEqual(
            concatdict, {"a": "AAAAGGGG", "b": "TTTT----", "c": "CCCCNNNN"}
        )

    def test_replace_seqs(self):
        """synchronize gaps between protein seqs and codon seqs"""
        pd = {
            "FlyingFox": "C-TNAH",
            "DogFaced": "CGTNT-",
            "FreeTaile": "-GTDTH",
            "LittleBro": "C-TD-H",
            "TombBat": "C--STH",
        }
        pal = make_aligned_seqs(moltype=PROTEIN, data=pd)

        cu = {
            "TombBat": "TGTAGTACTCAT",
            "FreeTaile": "GGCACAGATACTCAT",
            "FlyingFox": "TGTACAAATGCTCAT",
            "LittleBro": "TGTACAGATCAT",
            "DogFaced": "TGTGGCACAAATACT",
        }

        co = make_unaligned_seqs(moltype=DNA, data=cu)
        cal = pal.replace_seqs(co)
        result = cal.to_dict()
        for taxon, expected_sequence in [
            ("FlyingFox", "TGT---ACAAATGCTCAT"),
            ("DogFaced", "TGTGGCACAAATACT---"),
            ("FreeTaile", "---GGCACAGATACTCAT"),
            ("LittleBro", "TGT---ACAGAT---CAT"),
            ("TombBat", "TGT------AGTACTCAT"),
        ]:
            self.assertEqual(result[taxon], expected_sequence)

    def test_sample(self):
        """Test sample generation"""
        alignment = make_aligned_seqs(
            data={"seq1": "ABCDEFGHIJKLMNOP", "seq2": "ABCDEFGHIJKLMNOP"}
        )
        # effectively permute columns, preserving length
        alignment.sample()
        # ensure length correct
        sample = alignment.sample(10)
        self.assertEqual(len(sample), 10)
        # test columns alignment preserved
        seqs = list(sample.to_dict().values())
        self.assertEqual(seqs[0], seqs[1])
        # ensure each char occurs once as sampling without replacement
        for char in seqs[0]:
            self.assertEqual(seqs[0].count(char), 1)

    def test_sample_with_replacement(self):
        # test with replacement
        alignment = make_aligned_seqs(data={"seq1": "gatc", "seq2": "gatc"})
        sample = alignment.sample(1000, with_replacement=True)

    def test_sample_tuples(self):
        ##### test with motif size != 1 #####
        alignment = make_aligned_seqs(
            data={
                "seq1": "AABBCCDDEEFFGGHHIIJJKKLLMMNNOOPP",
                "seq2": "AABBCCDDEEFFGGHHIIJJKKLLMMNNOOPP",
            }
        )
        shuffled = alignment.sample(motif_length=2)
        # ensure length correct
        sample = alignment.sample(10, motif_length=2)
        self.assertEqual(len(sample), 20)
        # test columns alignment preserved
        seqs = list(sample.to_dict().values())
        self.assertEqual(seqs[0], seqs[1])
        # ensure each char occurs twice as sampling dinucs without replacement
        for char in seqs[0]:
            self.assertEqual(seqs[0].count(char), 2)

    def test_translate(self):
        for seqs in [
            {"seq1": "GATTTT", "seq2": "GATC??"},
            {"seq1": "GAT---", "seq2": "?GATCT"},
        ]:
            alignment = make_aligned_seqs(data=seqs, moltype=DNA)
            self.assertEqual(len(alignment.get_translation()), 2)
            # check for a failure when no moltype specified
            alignment = make_aligned_seqs(data=seqs)
            try:
                alignment.get_translation()
            except AttributeError:
                pass

    def test_seqnames(self):
        s1 = self.alignment.get_seq("Mouse")
        self.assertEqual(s1.get_name(), "Mouse")

    def test_trim_stop_codons(self):
        """test without terminal stop handling"""
        seq_coll = make_unaligned_seqs(
            data={"seq1": "ACGTAA", "seq2": "ACGACG", "seq3": "ACGCGT"}, moltype=DNA
        )
        seq_coll = seq_coll.trim_stop_codons()
        seqs = seq_coll.to_dict()
        self.assertEqual(seqs["seq1"], "ACG")  # note: not 'acg---'
        self.assertEqual(seqs["seq2"], "ACGACG")
        # aligned
        aln = make_aligned_seqs(
            data={"seq1": "ACGTAA", "seq2": "ACGTGA", "seq3": "ACGTAA"}, moltype=DNA
        )
        aln = aln.trim_stop_codons()
        self.assertEqual(
            aln.to_dict(), {"seq1": "ACG", "seq2": "ACG", "seq3": "ACG"}
        )  # note: not 'acg---'
        aln = make_aligned_seqs(
            data={"seq1": "ACGAAA", "seq2": "ACGTGA", "seq3": "ACGTAA"}, moltype=DNA
        )
        aln = aln.trim_stop_codons()
        self.assertEqual(
            aln.to_dict(), {"seq1": "ACGAAA", "seq2": "ACG---", "seq3": "ACG---"}
        )

        # for case where a sequence length is not divisible by 3
        seq_coll = make_unaligned_seqs(
            data={"seq1": "ACGTAA", "seq2": "ACGAC"}, moltype=DNA
        )
        # fail
        self.assertRaises(ValueError, seq_coll.trim_stop_codons)
        # unless explicitly over-ridden with allow_partial
        new_coll = seq_coll.trim_stop_codons(allow_partial=True)
        self.assertEqual(new_coll.to_dict(), dict(seq1="ACG", seq2="ACGAC"))

        # should work for alignments too
        aln = make_aligned_seqs(
            data={"seq1": "ACGTAA---", "seq2": "ACGAC----", "seq3": "ACGCAATTT"},
            moltype=DNA,
        )
        # fail
        self.assertRaises(ValueError, aln.trim_stop_codons)
        # unless explicitly over-ridden with allow_partial
        aln = aln.trim_stop_codons(allow_partial=True)
        self.assertEqual(
            aln.to_dict(),
            {"seq1": "ACG------", "seq2": "ACGAC----", "seq3": "ACGCAATTT"},
        )
        # mixed lengths
        aln = make_aligned_seqs(
            data={"seq1": "ACGTAA---", "seq2": "ACGAC----", "seq3": "ACGCAATGA"},
            moltype=DNA,
        )
        aln = aln.trim_stop_codons(allow_partial=True)
        self.assertEqual(
            aln.to_dict(), {"seq1": "ACG---", "seq2": "ACGAC-", "seq3": "ACGCAA"}
        )
        # longest seq not divisible by 3
        aln = make_aligned_seqs(
            data={"seq1": "ACGTAA--", "seq2": "ACGAC---", "seq3": "ACGC-ATG"},
            moltype=DNA,
        )
        aln = aln.trim_stop_codons(allow_partial=True)
        self.assertEqual(
            aln.to_dict(), {"seq1": "ACG-----", "seq2": "ACGAC---", "seq3": "ACGC-ATG"}
        )

    def test_trim_stop_codons_info(self):
        """trim_stop_codons should preserve info attribute"""
        seq_coll = SequenceCollection(
            data={"seq1": "ACGTAA", "seq2": "ACGACG", "seq3": "ACGCGT"},
            moltype=DNA,
            info={"key": "value"},
        )
        seq_coll = seq_coll.trim_stop_codons()
        self.assertEqual(seq_coll.info["key"], "value")

        # aligned
        aln = ArrayAlignment(
            data={"seq1": "ACGTAA", "seq2": "ACGTGA", "seq3": "ACGTAA"},
            moltype=DNA,
            info={"key": "value"},
        )
        aln = aln.trim_stop_codons()
        self.assertEqual(aln.info["key"], "value")

    def test_has_terminal_stops(self):
        """test truth values for terminal stops"""
        # seq collections
        seq_coll = make_unaligned_seqs(
            data={"seq1": "ACGTAA", "seq2": "ACG", "seq3": "ACGCGT"}, moltype=DNA
        )
        assert seq_coll.has_terminal_stops() == True
        seq_coll = make_unaligned_seqs(
            data={"seq1": "ACGTAC", "seq2": "ACGACG", "seq3": "ACGCGT"}, moltype=DNA
        )
        assert seq_coll.has_terminal_stops() == False
        # alignments
        aln = make_aligned_seqs(
            data={"seq1": "ACGTAA", "seq2": "ACGCAA", "seq3": "ACGCGT"}, moltype=DNA
        )
        assert aln.has_terminal_stops() == True
        aln = make_aligned_seqs(
            data={"seq1": "ACGTAA", "seq2": "ACGTAG", "seq3": "ACGTGA"}, moltype=DNA
        )
        assert aln.has_terminal_stops() == True
        aln = make_aligned_seqs(
            data={"seq1": "ACGCAA", "seq2": "ACGCAA", "seq3": "ACGCGT"}, moltype=DNA
        )
        assert aln.has_terminal_stops() == False

        # ValueError if ragged end
        aln = make_aligned_seqs(
            data={"seq1": "ACGCAA", "seq2": "ACGTAA", "seq3": "ACGCG-"}, moltype=DNA
        )
        self.assertRaises(ValueError, aln.has_terminal_stops)
        self.assertTrue(aln.has_terminal_stops(allow_partial=True))

    def test_slice(self):
        seqs = {"seq1": "ACGTANGT", "seq2": "ACGTACGT", "seq3": "ACGTACGT"}
        alignment = make_aligned_seqs(data=seqs)
        short = {"seq1": "A", "seq2": "A", "seq3": "A"}
        self.assertEqual(alignment[0:1].to_dict(), short)

    def test_get_motifprobs(self):
        """calculation of motif probs"""
        seqs = {"seq1": "ACGTANGT", "seq2": "-CGTACGT", "seq3": "ACGTACGT"}
        aln = make_aligned_seqs(data=seqs, moltype=DNA)
        mprobs = aln.get_motif_probs(allow_gap=False)
        expected = {"A": 5 / 22, "T": 6 / 22, "C": 5 / 22, "G": 6 / 22}
        self.assertEqual(mprobs, expected)
        mprobs = aln.get_motif_probs(allow_gap=True)
        expected = {"A": 5 / 23, "T": 6 / 23, "C": 5 / 23, "G": 6 / 23, "-": 1 / 23}
        self.assertEqual(mprobs, expected)
        mprobs = aln.get_motif_probs(allow_gap=False, include_ambiguity=True)
        expected = {"A": 5.25 / 23, "T": 6.25 / 23, "C": 5.25 / 23, "G": 6.25 / 23}
        self.assertEqual(mprobs, expected)
        mprobs = aln.get_motif_probs(allow_gap=True, include_ambiguity=True)
        expected = {
            "A": 5.25 / 24,
            "T": 6.25 / 24,
            "C": 5.25 / 24,
            "G": 6.25 / 24,
            "-": 1 / 24,
        }
        self.assertEqual(mprobs, expected)
        seqs = {"seq1": "ACGAANGA", "seq2": "-CGAACGA", "seq3": "ACGAACGA"}
        aln = make_aligned_seqs(data=seqs, moltype=DNA)
        mprobs = aln.get_motif_probs(exclude_unobserved=True)
        expected = {"A": 11 / 22, "C": 5 / 22, "G": 6 / 22}
        self.assertEqual(mprobs, expected)


# fileformats doesn't catch an exception when the file has no data!
class SequenceTestMethods(unittest.TestCase):
    """Testing Sequence methods"""

    def setUp(self):
        self.seq = make_seq("ATGACGTTGCGTAGCATAGCTCGA", "dna")

    def test_getlength(self):
        """testing getting length"""
        assert len(self.seq) == 24

    def test_get_in_motif_size(self):
        """test accuracy of chunking various sizes"""
        self.assertEqual(
            self.seq.get_in_motif_size(2),
            ["AT", "GA", "CG", "TT", "GC", "GT", "AG", "CA", "TA", "GC", "TC", "GA"],
        )
        self.assertEqual(
            self.seq.get_in_motif_size(3),
            ["ATG", "ACG", "TTG", "CGT", "AGC", "ATA", "GCT", "CGA"],
        )

    def test_translate(self):
        """test of translating seqs"""
        seq = make_seq("ATGACGTTGCGTAGCATAGCTCGA", moltype=DNA).get_translation()
        self.assertEqual(str(seq), "MTLRSIAR")

    def test_ambig_translate(self):
        """test of translating seqs"""
        seq = make_seq("CGNTGN???---", moltype=DNA).get_translation()
        self.assertEqual(str(seq), "RX?-")

    def test_translate_incomplete(self):
        """test of translating seqs with incomplete codon"""
        seq = make_seq("CGNTGNAC----", moltype=DNA)
        aa = seq.get_translation(incomplete_ok=True)
        self.assertEqual(str(aa), "RX?-")
        with self.assertRaises(AlphabetError):
            _ = seq.get_translation(incomplete_ok=False)

    def test_slidingWindows(self):
        """test sliding window along sequences"""
        result = []
        for bit in self.seq.sliding_windows(5, 2):
            result += [bit]
        self.assertEqual(
            [str(x) for x in result],
            [
                "ATGAC",
                "GACGT",
                "CGTTG",
                "TTGCG",
                "GCGTA",
                "GTAGC",
                "AGCAT",
                "CATAG",
                "TAGCT",
                "GCTCG",
            ],
        )

        result = []
        for bit in self.seq.sliding_windows(5, 1):
            result += [bit]
        self.assertEqual(
            [str(x) for x in result],
            [
                "ATGAC",
                "TGACG",
                "GACGT",
                "ACGTT",
                "CGTTG",
                "GTTGC",
                "TTGCG",
                "TGCGT",
                "GCGTA",
                "CGTAG",
                "GTAGC",
                "TAGCA",
                "AGCAT",
                "GCATA",
                "CATAG",
                "ATAGC",
                "TAGCT",
                "AGCTC",
                "GCTCG",
                "CTCGA",
            ],
        )

        result = []
        for bit in self.seq.sliding_windows(5, 1, start=3, end=6):
            result += [bit]
        self.assertEqual([str(x) for x in result], ["ACGTT", "CGTTG", "GTTGC"])

        # should not get a window when starting conditions don't generate one
        result = []
        for bit in self.seq.sliding_windows(20, 1, start=6):
            result += [bit]
        self.assertEqual(result, [])

    def test_reverse_complement(self):
        """testing reversal and complementing of a sequence"""
        seq = make_seq(moltype=DNA, seq="ACTGTAA")
        rev = seq.reverse_complement()
        self.assertEqual(str(rev), "TTACAGT")
        seq = make_seq(moltype=DNA, seq="ACTG-TAA")
        rev = seq.reverse_complement()
        self.assertEqual(str(rev), "TTA-CAGT")
        # try amigbuities
        seq = make_seq(moltype=DNA, seq="ACHNRTAA")
        rev = seq.reverse_complement()
        self.assertEqual(str(rev), "TTAYNDGT")

    def test_without_terminal_stop_sodon(self):
        """testing deleting terminal stop"""
        # for standard code
        seq = make_seq(moltype=DNA, seq="ACTTAA")
        seq2 = seq.trim_stop_codon()
        self.assertEqual(str(seq2), "ACT")

        # for sequence not divisible by 3
        seq = make_seq(moltype=DNA, seq="ACTTA")
        # fail
        self.assertRaises(ValueError, seq.trim_stop_codon)
        # unless explicitly over-ride length issue using allow_partial
        seq2 = seq.trim_stop_codon(allow_partial=True)

    def test_has_terminal_stop(self):
        """test check for terminal stop codons"""
        seq = make_seq(moltype=DNA, seq="ACTTAA")
        assert seq.has_terminal_stop() == True
        seq = make_seq(moltype=DNA, seq="ACTTAT") == False

        # for sequence not divisible by 3
        seq = make_seq(moltype=DNA, seq="ACTTA")
        # fail
        self.assertRaises(ValueError, seq.has_terminal_stop)
        # unless explicitly over-ride length issue using allow_partial
        # in which case, returns False
        self.assertFalse(seq.has_terminal_stop(allow_partial=True))


if __name__ == "__main__":
    unittest.main()
