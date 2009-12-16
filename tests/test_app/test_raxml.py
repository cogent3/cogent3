#!/bin/env python

from os import getcwd, remove, rmdir, mkdir
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from cogent.app.raxml import Raxml,raxml_alignment, build_tree_from_alignment
from cogent.app.util import ApplicationError
from cogent.parse.phylip import get_align_for_phylip
from cogent.core.tree import PhyloNode
from cogent.core.moltype import RNA
from StringIO import StringIO

__author__ = "Micah Hamady"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Micah Hamady", "Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
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

        self.test_model = "GTRCAT"

        self.align1 = get_align_for_phylip(StringIO(PHYLIP_FILE))

        self.test_fn1 = "/tmp/raxml_test1.txt"
        self.test_fn2 = "raxml_test1.txt"
        self.test_fn1_space = "/tmp/raxml test1.txt"

    def writeTmp(self, outname):
        """Write data to temp file"""
        t = open(outname, "w+")
        t.write(PHYLIP_FILE)
        t.close()


class RaxmlTests(GenericRaxml):
    """Tests for the Raxml application controller"""

    def test_raxml(self):
        """raxml BaseCommand should return the correct BaseCommand"""
        r = Raxml()
        self.assertEqual(r.BaseCommand, \
            ''.join(['cd \"',getcwd(),'/\"; ','raxmlHPC -e 0.1 -f d -c 50']))
        r.Parameters['-s'].on('seq.nexus')
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd \"',getcwd(),'/\"; ',\
            'raxmlHPC -e 0.1 -f d -c 50 -s seq.nexus']))


    def test_raxml_params(self):
        """raxml should raise exception if missing required params"""

        r = Raxml(WorkingDir="/tmp")

        r.SuppressStdout = True
        r.SuppressStderr = True
        # raise error by default
        self.assertRaises(ValueError, r)

        # specify output name 
        r.Parameters['-n'].on("test_name")
        self.assertRaises(ApplicationError, r)

        # specify model 
        r.Parameters['-m'].on("GTRCAT")
        self.assertRaises(ApplicationError, r)

        r.Parameters['-s'].on(self.test_fn1)
        self.assertRaises(ApplicationError, r)


        self.writeTmp(self.test_fn1)

        o = r()
        o.cleanUp()

        remove(self.test_fn1)

    
    def test_raxml_from_file(self):
        """raxml should run correctly using filename"""
        r = Raxml(WorkingDir="/tmp")

        r.Parameters['-s'].on(self.test_fn1)
        r.Parameters['-m'].on("GTRCAT")
        r.Parameters['-n'].on("test_me")
       
        # test with abs filename
        cur_out = self.test_fn1
        self.writeTmp(cur_out)
        out = r()
        out.cleanUp()
        remove(cur_out)

        # test with rel + working dir 
        r.Parameters['-s'].on(self.test_fn2)
        r.Parameters['-n'].on("test_me2")
        r.Parameters['-w'].on("/tmp/")
        self.writeTmp(self.test_fn1)
        out = r()
        out.cleanUp()
        remove(self.test_fn1)

        r.Parameters['-s'].on("\"%s\"" % self.test_fn1_space)
        r.Parameters['-n'].on("test_me3")
        r.Parameters['-w'].on("/tmp/")
        #print r.BaseCommand
        self.writeTmp(self.test_fn1_space)
        out = r()
        out.cleanUp()
        remove(self.test_fn1_space)

    def test_raxml_alignment(self):
        """raxml_alignment should work as expected"""
        phy_node, parsimony_phy_node, log_likelihood, total_exec \
            = raxml_alignment(self.align1)

    def test_build_tree_from_alignment(self):
        """Builds a tree from an alignment"""
        tree = build_tree_from_alignment(self.align1, RNA, False)
        self.assertTrue(isinstance(tree, PhyloNode))
        self.assertEqual(len(tree.tips()), 7)
        self.assertRaises(NotImplementedError, build_tree_from_alignment, \
                          self.align1, RNA, True)
   
PHYLIP_FILE= """ 7 50
Species001   UGCAUGUCAG UAUAGCUUUA GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
Species002   UGCAUGUCAG UAUAGCUUUA GUGAAACUGC GAAUGGCUNN UUAAAUCAGU
Species003   UGCAUGUCAG UAUAGCAUUA GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
Species004   UCCAUGUCAG UAUAACUUUG GUGAAACUGC GAAUGGCUCA UUAAAUCAGG
Species005   NNNNNNNNNN UAUAUCUUAU GUGAAACUUC GAAUGCCUCA UUAAAUCAGU
Species006   UGCAUGUCAG UAUAGCUUUG GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
Species007   UGCAUGUCAG UAUAACUUUG GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
""" 

if __name__ == '__main__':
    main()
