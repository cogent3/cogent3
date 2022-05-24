import json
import os
import pathlib
import shutil
import zipfile

from os.path import join
from tempfile import TemporaryDirectory
from unittest import TestCase, main

import numpy

from numpy.testing import assert_allclose

from cogent3 import DNA
from cogent3.app import align as align_app
from cogent3.app import io as io_app
from cogent3.app.composable import NotCompleted
from cogent3.app.data_store import DataStoreMember
from cogent3.app.io import write_db
from cogent3.app.result import generic_result
from cogent3.core.alignment import ArrayAlignment, SequenceCollection
from cogent3.core.profile import PSSM, MotifCountsArray, MotifFreqsArray
from cogent3.evolve.fast_distance import DistanceMatrix
from cogent3.maths.util import safe_log
from cogent3.util.table import Table


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def _get_generic_result(source):
    """creates a generic result with a DNA moltype as the single value"""
    gr = generic_result(source=source)
    gr["dna"] = DNA
    return gr


class TestIo(TestCase):
    basedir = "data"

    def tearDown(self) -> None:
        path = pathlib.Path("delme.tinydb")
        if path.exists():
            path.unlink()

    def test_findall(self):
        """find all files recursively"""
        found = list(io_app.findall(self.basedir, suffix=".fasta"))
        self.assertTrue(len(found) > 1)
        found = list(io_app.findall(self.basedir, suffix=".fasta", limit=2))
        self.assertTrue(len(found) == 2)

        # and with a suffix
        found = list(io_app.findall(self.basedir, suffix=".fasta*"))
        self.assertTrue(len(found) > 2)

    def test_findall_zip(self):
        """find all files recursively in a zip archive"""
        with TemporaryDirectory(dir=".") as dirname:
            zip_path = join(dirname, "new")
            shutil.make_archive(zip_path, "zip", self.basedir)
            zip_path = zip_path + ".zip"  # because shutil adds the suffix
            found = list(io_app.findall(zip_path, suffix="fasta"))
            self.assertTrue(len(found) > 1)
            found = list(io_app.findall(zip_path, suffix="fasta", limit=2))
            self.assertTrue(len(found) == 2)

            # and with a suffix
            found = list(io_app.findall(zip_path, suffix=".fasta*"))
            self.assertTrue(len(found) > 2)

    def test_define_data_store(self):
        """returns an iterable data store"""
        found = io_app.get_data_store(self.basedir, suffix=".fasta")
        self.assertTrue(len(found) > 1)
        found = io_app.get_data_store(self.basedir, suffix=".fasta", limit=2)
        self.assertTrue(len(found) == 2)

        # and with a suffix
        found = list(io_app.get_data_store(self.basedir, suffix=".fasta*"))
        self.assertTrue(len(found) > 2)

        # with a wild-card suffix
        found = list(io_app.get_data_store(self.basedir, suffix="*"))
        self.assertEqual(len(os.listdir(self.basedir)), len(found))

        # raises ValueError if suffix not provided or invalid
        with self.assertRaises(ValueError):
            _ = io_app.get_data_store(self.basedir)

        with self.assertRaises(ValueError):
            _ = io_app.get_data_store(self.basedir, 1)

    def test_load_aligned(self):
        """correctly loads aligned seqs"""

        def validate(paths, loader):
            loaded = list(map(loader, paths))
            for i, aln in enumerate(loaded):
                self.assertTrue(len(aln) > 10)
                self.assertIsInstance(aln, ArrayAlignment)
                self.assertEqual(aln.info.source, paths[i])

        fasta_paths = io_app.get_data_store(self.basedir, suffix=".fasta", limit=2)
        fasta_loader = io_app.load_aligned(format="fasta")
        validate(fasta_paths, fasta_loader)

    def test_load_aligned_nexus(self):
        """should handle nexus too"""
        nexus_paths = io_app.get_data_store(self.basedir, suffix="nex")
        loader = io_app.load_aligned(format="nexus")
        results = [loader(m) for m in nexus_paths]
        for result in results:
            self.assertIsInstance(result, ArrayAlignment)

    def test_load_aligned_paml(self):
        """should handle paml too"""
        paml_paths = io_app.get_data_store(self.basedir, suffix="paml")
        loader = io_app.load_aligned(format="paml")
        results = [loader(m) for m in paml_paths]
        for result in results:
            self.assertIsInstance(result, ArrayAlignment)

    def test_load_aligned_from_zip(self):
        """correctly loads aligned seqs from a zip archive"""

        def validate(paths, loader):
            loaded = list(map(loader, paths))
            for i, aln in enumerate(loaded):
                self.assertTrue(len(aln) > 10)
                self.assertIsInstance(aln, ArrayAlignment)
                # paths is only the basename when workjing with zip archives
                # whereas the inpath will have full path of zip archive
                self.assertEqual(aln.info.source, paths[i])
                self.assertEqual(aln.info.source, paths[i])

        with TemporaryDirectory(dir=".") as dirname:
            zip_path = join(dirname, self.basedir.replace(".zip", ""))
            shutil.make_archive(
                base_name=zip_path, root_dir=".", format="zip", base_dir=self.basedir
            )
            zip_path = zip_path + ".zip"  # because shutil adds the suffix
            fasta_paths = list(io_app.findall(zip_path, suffix=".fasta", limit=2))
            fasta_loader = io_app.load_aligned(format="fasta")
            validate(fasta_paths, fasta_loader)

    def test_load_unaligned(self):
        """load_unaligned returns degapped sequence collections"""
        fasta_paths = io_app.get_data_store(self.basedir, suffix=".fasta", limit=2)
        fasta_loader = io_app.load_unaligned(format="fasta")
        for i, seqs in enumerate(map(fasta_loader, fasta_paths)):
            self.assertIsInstance(seqs, SequenceCollection)
            self.assertTrue("-" not in "".join(seqs.to_dict().values()))
            self.assertEqual(seqs.info.source, fasta_paths[i])

        # returns NotCompleted when it's given an alignment/sequence
        # collection
        got = fasta_loader(seqs)
        self.assertIsInstance(got, NotCompleted)

    def test_write_seqs(self):
        """correctly writes sequences out"""
        fasta_paths = list(io_app.findall(self.basedir, suffix=".fasta", limit=2))
        fasta_loader = io_app.load_aligned(format="fasta")
        alns = list(map(fasta_loader, fasta_paths))
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_seqs(dirname, if_exists="ignore")
            wrote = list(map(writer, alns))
            self.assertIsInstance(wrote[0], DataStoreMember)
            written = list(io_app.findall(dirname, suffix="fasta"))
            for i, wrote in enumerate(written):
                self.assertEqual(alns[i].info.stored, join(dirname, wrote))

    def test_load_json(self):
        """correctly loads an object from json"""
        from cogent3.app.data_store import make_record_for_json

        data = make_record_for_json("delme", DNA, True)
        data = json.dumps(data)
        # straight directory
        with TemporaryDirectory(dir=".") as dirname:
            outpath = join(dirname, "delme.json")
            with open(outpath, "w") as outfile:
                outfile.write(data)
            reader = io_app.load_json()
            got = reader(outpath)
            self.assertIsInstance(got, DNA.__class__)
            self.assertEqual(got, DNA)

        # zipped directory
        with TemporaryDirectory(dir=".") as dirname:
            zip_path = join(dirname, "delme.zip")
            outpath = "delme/delme.json"
            with zipfile.ZipFile(zip_path, "a") as out:
                out.writestr(outpath, data)

            dstore = io_app.get_data_store(zip_path, suffix="json")
            member = dstore.get_member("delme.json")
            reader = io_app.load_json()
            got = reader(member)
            self.assertIsInstance(got, DNA.__class__)
            self.assertEqual(got, DNA)

    def test_write_db_load_db(self):
        """correctly write/load from tinydb"""
        # straight directory
        with TemporaryDirectory(dir=".") as dirname:
            outpath = join(dirname, "delme")
            writer = write_db(outpath, create=True, if_exists="ignore")
            gr = _get_generic_result(join("blah", "delme.json"))
            got = writer(gr)
            self.assertIsInstance(got, DataStoreMember)
            writer.data_store.db.close()
            dstore = io_app.get_data_store(f"{outpath}.tinydb", suffix="json")
            reader = io_app.load_db()
            got = reader(dstore[0])
            dstore.close()
            got.deserialised_values()
            self.assertIsInstance(got["dna"], DNA.__class__)
            self.assertEqual(got["dna"], DNA)

    def test_write_db_load_db2(self):
        """correctly write/load built-in python from tinydb"""
        with TemporaryDirectory(dir=".") as dirname:
            outpath = join(dirname, "delme")
            writer = write_db(outpath, create=True, if_exists="ignore")
            data = dict(a=[1, 2], b="string")
            m = writer(data, identifier=join("blah", "delme.json"))
            writer.data_store.db.close()
            dstore = io_app.get_data_store(f"{outpath}.tinydb", suffix="json")
            reader = io_app.load_db()
            got = reader(dstore[0])
            dstore.close()
            self.assertEqual(got, data)

    def test_load_db_failure_json_file(self):
        """informative load_db error message when given a json file path"""
        # todo this test has a trapped exception about being unable to delete
        # a file
        with TemporaryDirectory(dir=".") as dirname:
            outpath = join(dirname, "delme")
            writer = write_db(outpath, create=True, if_exists="ignore")
            gr = _get_generic_result(join("blah", "delme.json"))
            got = writer(gr)
            writer.data_store.db.close()
            dstore = io_app.get_data_store(f"{outpath}.tinydb", suffix="json")
            reader = io_app.load_db()
            outpath = join(dirname, "dummy.json")
            with open(outpath, mode="w") as outfile:
                outfile.write("\n\n")

            got = reader(outpath)
            self.assertIsInstance(got, NotCompleted)
            self.assertTrue("json" in got.message)

    def test_load_tabular(self):
        """correctly loads tabular data"""
        rows = [[1, 2], [3, 4], [5, 6.5]]
        table = Table(["A", "B"], data=rows)
        load_table = io_app.load_tabular(sep="\t", with_header=True)
        with TemporaryDirectory(dir=".") as dirname:
            outpath = join(dirname, "delme.tsv")
            table.write(outpath)
            new = load_table(outpath)
            self.assertEqual(new.title, "")
            self.assertEqual(type(new[0, "B"]), type(table[0, "B"]))
            self.assertEqual(type(new[0, "A"]), type(table[0, "A"]))
            outpath = join(dirname, "delme2.tsv")
            with open(outpath, "w") as out:
                out.write("\t".join(table.header[:1]) + "\n")
                for row in table.tolist():
                    row = "\t".join(map(str, row))
                    out.write(row + "\n")
            result = load_table(outpath)
            self.assertIsInstance(result, NotCompleted)

    def test_write_tabular_motif_counts_array(self):
        """correctly writes tabular data for MotifCountsArray"""

        data = [[2, 4], [3, 5], [4, 8]]
        mca = MotifCountsArray(data, "AB")
        loader = io_app.load_tabular(sep="\t")
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_tabular(data_path=dirname, format="tsv")
            outpath = join(dirname, "delme.tsv")
            got = writer.write(mca, identifier=outpath)
            self.assertIsInstance(got, DataStoreMember)
            new = loader(outpath)
            # when written to file in tabular form
            # the loaded table will have dim-1 dim-2 as column labels
            # and the key-values pairs listed below; in dict form...
            expected = {
                0: {"dim-1": 0, "dim-2": "A", "value": 2},
                1: {"dim-1": 0, "dim-2": "B", "value": 4},
                2: {"dim-1": 1, "dim-2": "A", "value": 3},
                3: {"dim-1": 1, "dim-2": "B", "value": 5},
                4: {"dim-1": 2, "dim-2": "A", "value": 4},
                5: {"dim-1": 2, "dim-2": "B", "value": 8},
            }
            self.assertEqual(expected, new.to_dict())

    def test_write_tabular_motif_freqs_array(self):
        """correctly writes tabular data for MotifFreqsArray"""

        data = [[0.3333, 0.6667], [0.3750, 0.625], [0.3333, 0.6667]]
        mfa = MotifFreqsArray(data, "AB")
        loader = io_app.load_tabular(sep="\t")
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_tabular(data_path=dirname, format="tsv")
            outpath = join(dirname, "delme.tsv")
            writer.write(mfa, identifier=outpath)
            new = loader(outpath)
            # when written to file in tabular form
            # the loaded table will have dim-1 dim-2 as column labels
            # and the key-values pairs listed below; in dict form...
            expected = {
                0: {"dim-1": 0, "dim-2": "A", "value": 0.3333},
                1: {"dim-1": 0, "dim-2": "B", "value": 0.6667},
                2: {"dim-1": 1, "dim-2": "A", "value": 0.3750},
                3: {"dim-1": 1, "dim-2": "B", "value": 0.6250},
                4: {"dim-1": 2, "dim-2": "A", "value": 0.3333},
                5: {"dim-1": 2, "dim-2": "B", "value": 0.6667},
            }
            self.assertEqual(expected, new.to_dict())

    def test_write_tabular_pssm(self):
        """correctly writes tabular data for PSSM"""

        # data from test_profile
        data = numpy.array(
            [
                [0.1, 0.3, 0.5, 0.1],
                [0.25, 0.25, 0.25, 0.25],
                [0.05, 0.8, 0.05, 0.1],
                [0.7, 0.1, 0.1, 0.1],
                [0.6, 0.15, 0.05, 0.2],
            ]
        )
        pssm = PSSM(data, "ACTG")
        loader = io_app.load_tabular(sep="\t")
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_tabular(data_path=dirname, format="tsv")
            outpath = join(dirname, "delme.tsv")
            writer.write(pssm, identifier=outpath)
            new = loader(outpath)
            expected = safe_log(data) - safe_log(numpy.array([0.25, 0.25, 0.25, 0.25]))
            for i in range(len(expected)):
                j = i // 4
                self.assertTrue(
                    numpy.isclose(new.array[i][2], expected[j][i - j], atol=0.0001)
                )

    def test_write_tabular_distance_matrix(self):
        """correctly writes tabular data for DistanceMatrix"""
        data = {(0, 0): 0, (0, 1): 4, (1, 0): 4, (1, 1): 0}
        matrix = DistanceMatrix(data)
        loader = io_app.load_tabular(sep="\t")
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_tabular(data_path=dirname, format="tsv")
            outpath = join(dirname, "delme.tsv")
            writer.write(matrix, identifier=outpath)
            new = loader(outpath)
            # when written to file in tabular form
            # the loaded table will have dim-1 dim-2 as column labels
            # and the key-values pairs listed below; in dict form...
            expected = {
                0: {"dim-1": 0, "dim-2": 1, "value": 4},
                1: {"dim-1": 1, "dim-2": 0, "value": 4},
            }
            self.assertEqual(expected, new.to_dict())

    def test_write_tabular_table(self):
        """correctly writes tabular data"""
        rows = [[1, 2], [3, 4], [5, 6.5]]
        table = Table(["A", "B"], data=rows)
        loader = io_app.load_tabular(sep="\t")
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_tabular(data_path=dirname, format="tsv")
            outpath = join(dirname, "delme.tsv")
            writer.write(table, identifier=outpath)
            new = loader(outpath)
            self.assertEqual(table.to_dict(), new.to_dict())

    def test_load_tabular_motif_counts_array(self):
        """correctly loads tabular data for MotifCountsArray"""

        data = [[2, 4], [3, 5], [4, 8]]
        mca = MotifCountsArray(data, "AB")
        loader = io_app.load_tabular(sep="\t", as_type="motif_counts")
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_tabular(data_path=dirname, format="tsv")
            outpath = join(dirname, "delme.tsv")
            writer.write(mca, identifier=outpath)
            new = loader(outpath)
            self.assertEqual(mca.to_dict(), new.to_dict())

    def test_load_tabular_motif_freqs_array(self):
        """correctly loads tabular data for MotifFreqsArray"""

        data = [[0.3333, 0.6667], [0.3750, 0.625], [0.3333, 0.6667]]
        mfa = MotifFreqsArray(data, "AB")
        loader = io_app.load_tabular(sep="\t", as_type="motif_freqs")
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_tabular(data_path=dirname, format="tsv")
            outpath = join(dirname, "delme.tsv")
            writer.write(mfa, identifier=outpath)
            new = loader(outpath)
            self.assertEqual(mfa.to_dict(), new.to_dict())

    def test_load_tabular_pssm(self):
        """correctly loads tabular data for PSSM"""

        # data from test_profile
        data = [
            [0.1, 0.3, 0.5, 0.1],
            [0.25, 0.25, 0.25, 0.25],
            [0.05, 0.8, 0.05, 0.1],
            [0.7, 0.1, 0.1, 0.1],
            [0.6, 0.15, 0.05, 0.2],
        ]
        pssm = PSSM(data, "ACTG")
        loader = io_app.load_tabular(sep="\t", as_type="pssm")
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_tabular(data_path=dirname, format="tsv")
            outpath = join(dirname, "delme.tsv")
            writer.write(pssm, identifier=outpath)
            new = loader(outpath)
            assert_allclose(pssm.array, new.array, atol=0.0001)

    def test_load_tabular_distance_matrix(self):
        """correctly loads tabular data for DistanceMatrix"""
        data = {(0, 0): 0, (0, 1): 4, (1, 0): 4, (1, 1): 0}
        matrix = DistanceMatrix(data)
        loader = io_app.load_tabular(sep="\t", as_type="distances")
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_tabular(data_path=dirname, format="tsv")
            outpath = join(dirname, "delme.tsv")
            writer.write(matrix, identifier=outpath)
            new = loader(outpath)
            self.assertEqual(matrix.to_dict(), new.to_dict())

    def test_load_tabular_table(self):
        """correctly loads tabular data"""
        rows = [[1, 2], [3, 4], [5, 6.5]]
        table = Table(["A", "B"], data=rows)
        loader = io_app.load_tabular(sep="\t", as_type="table")
        with TemporaryDirectory(dir=".") as dirname:
            writer = io_app.write_tabular(data_path=dirname, format="tsv")
            outpath = join(dirname, "delme.tsv")
            writer.write(table, identifier=outpath)
            new = loader(outpath)
            self.assertEqual(table.to_dict(), new.to_dict())

    def test_write_json_with_info(self):
        """correctly writes an object with info attribute from json"""
        # create a mock object that pretends like it's been derived from
        # something
        from cogent3.app.result import generic_result

        with TemporaryDirectory(dir=".") as dirname:
            outdir = join(dirname, "delme")

            obj = generic_result(source=join("blah", "delme.json"))
            obj["dna"] = DNA
            writer = io_app.write_json(outdir, create=True)
            got = writer(obj)
            self.assertIsInstance(got, DataStoreMember)
            reader = io_app.load_json()
            got = reader(join(outdir, "delme.json"))
            got.deserialised_values()
            self.assertEqual(got["dna"], DNA)

    def test_write_json_no_info(self):
        """correctly writes an object with out an info attribute from json"""
        # create a mock object that pretends like it's been derived from
        # something
        with TemporaryDirectory(dir=".") as dirname:
            outdir = join(dirname, "delme")
            gr = _get_generic_result(join("blah", "delme.json"))
            writer = io_app.write_json(outdir, create=True)
            _ = writer(gr)
            reader = io_app.load_json()
            got = reader(writer.data_store[0])
            got.deserialised_values()
            self.assertEqual(got["dna"], DNA)

    def test_restricted_usage_of_tinydb_suffix(self):
        """can only use tinydb in a load_db, write_db context"""
        with TemporaryDirectory(dir=".") as dirname:
            outdir = join(dirname, "delme.tinydb")
            for writer_class in (
                io_app.write_seqs,
                io_app.write_json,
                io_app.write_tabular,
            ):
                with self.assertRaises(ValueError):
                    writer_class(outdir, create=True, if_exists="skip")
            # but OK for write_db
            w = io_app.write_db(outdir, create=True, if_exists="skip")
            w.data_store.close()

    def test_write_db_parallel(self):
        """writing with overwrite in parallel should reset db"""
        dstore = io_app.get_data_store(self.basedir, suffix="fasta")
        members = dstore.filtered(callback=lambda x: "brca1.fasta" not in x.split("/"))
        reader = io_app.load_unaligned()
        aligner = align_app.align_to_ref()
        writer = write_db("delme.tinydb", create=True, if_exists="overwrite")
        process = reader + aligner + writer

        r = process.apply_to(members, show_progress=False, parallel=True)

        expect = [str(m) for m in process.data_store]
        process.data_store.close()

        # now get read only and check what's in there
        result = io_app.get_data_store("delme.tinydb")
        got = [str(m) for m in result]
        result.close()

        self.assertEqual(got, expect)


if __name__ == "__main__":
    main()
