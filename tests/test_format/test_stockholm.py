 #!/usr/bin/env python
"""Tests for Stockholm sequence format writer.
"""
from cogent.util.unit_test import TestCase, main
from cogent.format.stockholm import stockholm_from_alignment
from cogent.core.alignment import Alignment
from cogent.core.sequence import Sequence
from cogent.core.info import Info

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"

class StockholmTests(TestCase):
    """Tests for Stockholm writer.
    """
    def setUp(self):
        """Setup for Stockholm tests."""
        self.unaligned_dict = {'1st':'AAA','2nd':'CCCC','3rd':'GGGG',
            '4th':'UUUU'}
        self.alignment_dict = {'1st':'AAAA','2nd':'CCCC','3rd':'GGGG',
            '4th':'UUUU'}
        #create alignment change order.
        self.alignment_object = Alignment(self.alignment_dict)
        self.alignment_order = ['2nd','4th','3rd','1st']
        self.alignment_object.RowOrder=self.alignment_order
        
        self.gc_annotation = {'SS_cons':'....'}
        self.stockholm_with_label=\
"""# STOCKHOLM 1.0

1st    AAAA
2nd    CCCC
3rd    GGGG
4th    UUUU
//"""
        self.stockholm_with_label_lw2=\
"""# STOCKHOLM 1.0

1st    AA
2nd    CC
3rd    GG
4th    UU

1st    AA
2nd    CC
3rd    GG
4th    UU
//"""

        self.stockholm_with_label_struct=\
"""# STOCKHOLM 1.0

1st             AAAA
2nd             CCCC
3rd             GGGG
4th             UUUU
#=GC SS_cons    ....
//"""
        self.stockholm_with_label_struct_lw2=\
"""# STOCKHOLM 1.0

1st             AA
2nd             CC
3rd             GG
4th             UU
#=GC SS_cons    ..

1st             AA
2nd             CC
3rd             GG
4th             UU
#=GC SS_cons    ..
//"""

        self.stockholm_with_label_reordered=\
"""# STOCKHOLM 1.0

2nd    CCCC
4th    UUUU
3rd    GGGG
1st    AAAA
//"""

        self.stockholm_with_label_lw2_reordered=\
"""# STOCKHOLM 1.0

2nd    CC
4th    UU
3rd    GG
1st    AA

2nd    CC
4th    UU
3rd    GG
1st    AA
//"""

    def test_stockholm_from_alignment_unaligned(self):
        """should raise error with unaligned seqs."""
        self.assertRaises(ValueError,\
            stockholm_from_alignment,self.unaligned_dict)
    
    def test_stockholm_from_alignment(self):
        """should return correct stockholm string."""
        self.assertEqual(stockholm_from_alignment({}),'')
        self.assertEqual(stockholm_from_alignment(self.alignment_dict),\
            self.stockholm_with_label)
        self.assertEqual(stockholm_from_alignment(self.alignment_dict,
                interleave_len=2),self.stockholm_with_label_lw2)
    
    def test_stockholm_from_alignment_struct(self):
        """should return correct stockholm string."""
        self.assertEqual(stockholm_from_alignment({},\
            GC_annotation=self.gc_annotation),'')
        self.assertEqual(stockholm_from_alignment(self.alignment_dict,\
            GC_annotation=self.gc_annotation),\
            self.stockholm_with_label_struct)
        self.assertEqual(stockholm_from_alignment(self.alignment_dict,\
            GC_annotation=self.gc_annotation,\
            interleave_len=2),self.stockholm_with_label_struct_lw2)
    
    def test_stockholm_from_alignment_reordered(self):
        """should return correct stockholm string."""
        self.assertEqual(stockholm_from_alignment(self.alignment_object),\
            self.stockholm_with_label_reordered)
        self.assertEqual(stockholm_from_alignment(self.alignment_object,
                interleave_len=2),self.stockholm_with_label_lw2_reordered)

if __name__ == "__main__":
    main()
