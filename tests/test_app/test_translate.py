from unittest import TestCase, main

from cogent3 import LoadSeqs, DNA
from cogent3.app.translate import (best_frame, select_translatable,
                                   translate_frames,
                                   get_code,
                                   get_fourfold_degenerate_sets)

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestTranslate(TestCase):
    """testing translation functions"""

    def test_best_frame(self):
        """correctly identify best frame with/without allowing rc"""
        make_seq = DNA.make_seq
        seq = make_seq('ATGCTAACATAAA', name='fake1')
        f = best_frame(seq)
        self.assertEqual(f, 1)
        f = best_frame(seq, require_stop=True)
        self.assertEqual(f, 1)

        # a challenging seq, translatable in 1 and 3 frames, ending on stop in
        # frame 1. Should return frame 1 irrespective of require_stop
        seq = make_seq('ATGTTACGGACGATGCTGAAGTCGAAGATCCACCGCGCCACGGTGACCTGCTGA')
        f = best_frame(seq)
        self.assertEqual(f, 1)

        # a rc seq
        f = best_frame(seq)
        seq = make_seq(
            'AATATAAATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTTCATAAAGTCATA',
            name='fake2')
        f = best_frame(seq, allow_rc=True)
        self.assertEqual(f, 1)
        with self.assertRaises(ValueError):
            f = best_frame(seq, allow_rc=True, require_stop=True)

        rc = seq.rc()
        f = best_frame(rc, allow_rc=True)
        self.assertEqual(f, -1)

    def test_select_translatable(self):
        """correctly get translatable seqs"""
        data = {'a': 'AATATAAATGCCAGCTCATTACAGCATGAGAACA'
                     'GCAGTTTATTACTTCATAAAGTCATA',
                'rc': 'TATGACTTTATGAAGTAATAAACTGCTGTTCTCA'
                      'TGCTGTAATGAGCTGGCATTTATATT'}
        seqs = LoadSeqs(data=data, moltype=DNA, aligned=False)
        trans = select_translatable(allow_rc=False)
        tr = trans(seqs)
        ex = data.copy()
        ex.pop('rc')
        self.assertEqual(tr.todict(), ex)
        trans = select_translatable(allow_rc=True)
        tr = trans(seqs)
        ex = data.copy()
        ex['rc'] = data['a']
        self.assertEqual(tr.todict(), ex)

    def test_translate_frames(self):
        """returns translated sequences"""
        seq = DNA.make_seq('ATGCTGACATAAA', name='fake1')
        tr = translate_frames(seq)
        self.assertEqual(tr, ['MLT*', 'C*HK', 'ADI'])
        # with the bacterial nuclear and plant plastid code
        tr = translate_frames(seq, gc='Euplotid Nuclear')
        self.assertEqual(tr, ['MLT*', 'CCHK', 'ADI'])


class TestFourFoldDegen(TestCase):
    def test_get_fourfold_degenerate_sets(self):
        """correctly identify 4-fold degenerate codons"""
        # using straight characters
        expect = set()
        for di in 'GC', 'GG', 'CT', 'CC', 'TC', 'CG', 'AC', 'GT':
            expect.update([frozenset([di + n for n in 'ACGT'])])

        for i in range(1, 3):
            got = get_fourfold_degenerate_sets(get_code(i), as_indices=False)
            self.assertEqual(got, expect)

        with self.assertRaises(AssertionError):
            # as_indices requires an alphabet
            get_fourfold_degenerate_sets(get_code(1), as_indices=True)

        expect = set()
        for di in 'GC', 'GG', 'CT', 'CC', 'TC', 'CG', 'AC', 'GT':
            codons = list(map(lambda x: tuple(DNA.alphabet.to_indices(x)),
                              [di + n for n in 'ACGT']))
            expect.update([frozenset(codons)])

        for i in range(1, 3):
            got = get_fourfold_degenerate_sets(get_code(i),
                                               alphabet=DNA.alphabet,
                                               as_indices=True)
            self.assertEqual(got, expect)


if __name__ == '__main__':
    main()
