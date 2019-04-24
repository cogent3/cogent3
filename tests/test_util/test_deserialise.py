import json
from tempfile import TemporaryDirectory

from cogent3.util.unit_test import main, TestCase
from cogent3.util.deserialise import deserialise_object
from cogent3.evolve.models import get_model
from cogent3.core import moltype, alignment
from cogent3 import LoadSeqs, LoadTree
from cogent3.app.result import model_result

class TestDeserialising(TestCase):
    def test_roundtrip_codon_alphabet(self):
        """codon alphabet to_json enables roundtrip"""
        data = moltype.STANDARD_CODON.to_json()
        got = deserialise_object(data)
        self.assertEqual(type(got), type(moltype.STANDARD_CODON))
        self.assertEqual(list(got), list(moltype.STANDARD_CODON))

    def test_roundtrip_alphabet(self):
        """alphabet to_json enables roundtrip"""
        dna = moltype.get_moltype('dna')
        data = dna.alphabet.to_json()
        got = deserialise_object(data)
        self.assertEqual(type(got), type(dna.alphabet))
        self.assertEqual(list(got), list(dna.alphabet))

    def test_roundtrip_moltype(self):
        """moltype to_json enables roundtrip"""
        dna = moltype.get_moltype('dna')
        data = dna.to_json()
        got = deserialise_object(data)
        self.assertEqual(type(got), type(dna))
        self.assertEqual(list(got), list(dna))
        self.assertEqual(dna, got)

    def test_roundtrip_seq(self):
        """seq to_json enables roundtrip"""
        for mtype in ('dna', 'protein'):
            mtype = moltype.get_moltype(mtype)
            seq = mtype.make_seq('ACGGTCGG', 'label', info={'something': 3})
            got = deserialise_object(seq.to_json())
            self.assertEqual(got.info.something, 3)
            self.assertEqual(got.name, 'label')
            self.assertEqual(got.moltype, seq.moltype)
            self.assertEqual(str(got), str(seq))

    def test_roundtrip_seqcoll(self):
        """SequenceCollection to_json enables roundtrip"""
        data = dict(A='TTGT', B='GGCT')
        seqcoll = LoadSeqs(data=data, moltype='dna', aligned=False)
        got = deserialise_object(seqcoll.to_json())
        self.assertEqual(got.rc().todict(), seqcoll.rc().todict())
        self.assertIsInstance(got, alignment.SequenceCollection)

    def test_roundtrip_arrayalign(self):
        """ArrayAlignment to_json enables roundtrip"""
        data = dict(A='TTGTA', B='GGCT-')
        arrayalign = LoadSeqs(data=data, moltype='dna')
        got = deserialise_object(arrayalign.to_json())
        self.assertEqual(got.rc().todict(), arrayalign.rc().todict())
        self.assertIsInstance(got, alignment.ArrayAlignment)

    def test_roundtrip_align(self):
        """Alignment to_json enables roundtrip"""
        data = dict(A='TTGTA', B='GGCT-')
        align = LoadSeqs(data=data, moltype='dna', array_align=False)
        got = deserialise_object(align.to_json())
        self.assertEqual(got.rc().todict(), align.rc().todict())
        self.assertIsInstance(got, alignment.Alignment)
    
    def test_roundtrip_tree(self):
        """Tree to_json enables roundtrip"""
        tree = LoadTree(treestring='(c:01,d:0.3,(a:0.05,b:0.08)xx:0.2)')
        got = deserialise_object(tree.to_json())
        self.assertFloatEqual(got.get_node_matching_name('a').length, 0.05)
        self.assertFloatEqual(got.get_node_matching_name('xx').length, 0.2)

    def test_roundtrip_submod(self):
        """substitution model to_json enables roundtrip"""
        sm = get_model('HKY85')
        data = sm.to_json()
        got = deserialise_object(data)
        self.assertEqual(got.to_rich_dict(), sm.to_rich_dict())
        sm = get_model('GN')
        data = sm.to_json()
        got = deserialise_object(data)
        self.assertEqual(got.to_rich_dict(), sm.to_rich_dict())
        sm = get_model('CNFGTR')
        data = sm.to_json()
        got = deserialise_object(data)
        self.assertEqual(got.to_rich_dict(), sm.to_rich_dict())

    def test_roundtrip_likelihood_function(self):
        """likelihood function.to_json enables roundtrip"""
        _data = {'Human': 'ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG',
                 'Mouse': 'ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG',
                 'Opossum': 'ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG'}
        aln = LoadSeqs(data=_data, moltype='dna')
        tree = LoadTree(tip_names=aln.names)
        sm = get_model('HKY85')
        lf = sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        edge_vals = zip(aln.names, (2, 3, 4))
        for edge, val in edge_vals:
            lf.set_param_rule('kappa', edge=edge, init=val)
        lnL = lf.get_log_likelihood()
        data = lf.to_json()
        got_obj = deserialise_object(data)
        self.assertFloatEqual(got_obj.get_log_likelihood(), lnL)
        
    def test_roundtrip_from_file(self):
        """correctly roundtrips a likelihood function fro json file"""
        _data = {'Human': 'ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG',
                 'Mouse': 'ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG',
                 'Opossum': 'ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG'}
        aln = LoadSeqs(data=_data, moltype='dna')
        tree = LoadTree(tip_names=aln.names)
        sm = get_model('HKY85')
        lf = sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        edge_vals = zip(aln.names, (2, 3, 4))
        for edge, val in edge_vals:
            lf.set_param_rule('kappa', edge=edge, init=val)
        lnL = lf.get_log_likelihood()
        data = lf.to_json()
        with TemporaryDirectory(dir='.') as dirname:
            outpath = dirname + '/delme.json'
            with open(outpath, 'w') as outfile:
                outfile.write(data)

            got = deserialise_object(outpath)
            self.assertFloatEqual(got.get_log_likelihood(), lnL)

    def test_roundtrip_model_result(self):
        """mode_result.to_json enables roundtrip and lazy evaluation"""
        _data = {'Human': 'ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG',
                 'Mouse': 'ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG',
                 'Opossum': 'ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG'}
        aln = LoadSeqs(data=_data, moltype='dna')
        tree = LoadTree(tip_names=aln.names)
        sm = get_model('HKY85')
        lf = sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        edge_vals = zip(aln.names, (2, 3, 4))
        for edge, val in edge_vals:
            lf.set_param_rule('kappa', edge=edge, init=val)
        result = model_result(name='test')
        result[1] = lf
        self.assertIs(result[1], lf)
        self.assertEqual(result.nfp, lf.nfp)
        self.assertEqual(result.lnL, lf.lnL)

        data = result.to_json()
        got_obj = deserialise_object(data)
        # lazy evaluation means initially, the value is a dict
        self.assertIsInstance(got_obj[1], dict)
        # and properties match original
        self.assertEqual(got_obj.lnL, result.lnL)
        self.assertEqual(got_obj.nfp, result.nfp)
        self.assertEqual(got_obj.DLC, result.DLC)
        # when we ask for the lf attribute, it's no longer a dict
        self.assertNotIsInstance(got_obj.lf, dict)
        self.assertEqual(got_obj.lf.nfp, got_obj.nfp)


if __name__ == '__main__':
    main()
