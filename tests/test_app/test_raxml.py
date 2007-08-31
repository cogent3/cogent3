#!/bin/env python

from os import getcwd, remove, rmdir, mkdir
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from cogent.app.raxml import Raxml,raxml_alignment
from cogent.parse.phylip import get_align_for_phylip
from StringIO import StringIO

__author__ = "Micah Hamady"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Micah Hamady", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Micah Hamady"
__email__ = "Micah Hamady"
__status__ = "Development"

class GenericRaxml(TestCase):

    def setUp(self):
        """Setup data for raxml tests"""
        self.seqs1 = ['ACUGCUAGCUAGUAGCGUACGUA','GCUACGUAGCUAC',
            'GCGGCUAUUAGAUCGUA']
        self.labels1 = ['>1','>2','>3']
        self.lines1 = flatten(zip(self.labels1,self.seqs1))
        #self.dnd1 = DND1

        try:
            mkdir('/tmp/ct')
        except OSError: #dir already exists
            pass

        try:
            #create sequence files
            f = open('/tmp/ct/seq1.txt','w')
            f.write('\n'.join(self.lines1))
            f.close()

        except OSError:
            pass

        self.align1 = get_align_for_phylip(StringIO(PHYLIP_FILE))



class RaxmlTests(GenericRaxml):
    """Tests for the Raxml application controller"""

    def test_raxml(self):
        """raxml BaseCommand should return the correct BaseCommand"""
        r = Raxml()
        self.assertEqual(r.BaseCommand, \
            ''.join(['cd ',getcwd(),'/; ','raxmlHPC -e 0.1 -f d -c 50']))
        r.Parameters['-s'].on('seq.nexus')
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd ',getcwd(),'/; ',\
            'raxmlHPC -e 0.1 -f d -c 50 -s seq.nexus']))

   
class raxmlTests(GenericRaxml):
    """Tests for module level functions in raxml.py"""
   
       
    def test_raxml_alignment(self):
        """raxml_alignment should work as expected"""
        phy_node, parsimony_phy_node, log_likelihood, total_exec \
            = raxml_alignment(self.align1)

        #note: we get different results on different runs so can't test this
        #formally, but the above shouldn't raise an error
        #
        #print phy_node
        #print parsimony_phy_node
        #print log_likelihood
        #print total_exec
        
PHYLIP_FILE= """ 7 50
Species001   UGCAUGUCAG UAUAGCUUUA GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
Species002   UGCAUGUCAG UAUAGCUUUA GUGAAACUGC GAAUGGCUNN UUAAAUCAGU
Species003   UGCAUGUCAG UAUAGCAUUA GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
Species004   UGCAUGUCAG UAUAACUUUG GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
Species005   NNNNNNNNNN UAUAUCUUAU GUGAAACUUC GAAUGCCUCA UUAAAUCAGU
Species006   UGCAUGUCAG UAUAGCUUUG GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
Species007   UGCAUGUCAG UAUAACUUUG GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
""" 

if __name__ == '__main__':
    main()
