import os
from cogent.util.unit_test import TestCase, main

from cogent.db.ensembl.host import HostAccount, get_ensembl_account
from cogent.db.ensembl.assembly import Coordinate, CoordSystem, \
                                    get_coord_conversion
from cogent.db.ensembl.genome import Genome

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

Release = 56

if 'ENSEMBL_ACCOUNT' in os.environ:
    username, password = os.environ['ENSEMBL_ACCOUNT'].split()
    account = HostAccount('127.0.0.1', username, password)
else:
    account = get_ensembl_account(release=Release)


human = Genome(Species = 'human', Release=Release, account=account)
platypus = Genome(Species = 'platypus', Release=Release, account=account)

class TestLocation(TestCase):
    def test_init(self):
        human_loc = Coordinate(CoordName='x', Start=1000, End=10000, Strand=-1,
                        genome = human)
        # TODO: complete test for platpus
        self.assertEqual(human_loc.CoordType, 'chromosome')
        self.assertEqual(human_loc.CoordName, 'x')
        self.assertEqual(human_loc.Start, 1000)
        self.assertEqual(human_loc.End, 10000)
        self.assertEqual(human_loc.Strand, -1)
        self.assertEqual(human_loc.Species, "Homo sapiens")
        self.assertEqual(human_loc.seq_region_id, 27516)
    
    def test_get_coord_conversion(self):
        """should correctly map between different coordinate levels"""
        # not really testing the contig coordinates are correct
        CoordName, Start, End, Strand = '1', 1000, 1000000, 1
        human_loc = Coordinate(CoordName = CoordName, Start = Start, End = End,
                        Strand = Strand, genome = human)
        results = get_coord_conversion(human_loc, 'contig', human.CoreDb)
        for result in results:
            self.assertTrue(result[0].CoordName == CoordName)
            self.assertTrue(result[0].Start >= Start)
            self.assertTrue(result[0].End <= End)
            self.assertTrue(result[0].Strand == Strand)
        
    def test_coord_shift(self):
        """adding coordinates should produce correct results"""
        CoordName, Start, End, Strand = '1', 1000, 1000000, 1
        loc1 = Coordinate(CoordName = CoordName, Start = Start, End = End,
                        Strand = Strand, genome = human)
        for shift in [100, -100]:
            loc2 = loc1.shifted(shift)
            self.assertEqual(loc2.Start, loc1.Start+shift)
            self.assertEqual(loc2.End, loc1.End+shift)
            self.assertEqual(id(loc1.genome), id(loc2.genome))
        self.assertNotEqual(id(loc1), id(loc2))
    
    def test_coord_resize(self):
        """resizing should work"""
        CoordName, Start, End, Strand = '1', 1000, 1000000, 1
        loc1 = Coordinate(CoordName = CoordName, Start = Start, End = End,
                        Strand = Strand, genome = human)
        front_shift = -100
        back_shift = 100
        loc2 = loc1.resized(front_shift, back_shift)
        self.assertEqual(len(loc2), len(loc1)+200)
        self.assertEqual(loc2.Start, loc1.Start+front_shift)
        self.assertEqual(loc2.End, loc1.End+back_shift)
        self.assertEqual(loc1.Strand, loc2.Strand)
    
    def test_adopted(self):
        """coordinate should correctly adopt seq_region_id properties of 
        provided coordinate"""
        CoordName, Start, End, Strand = '1', 1000, 1000000, 1
        c1 = Coordinate(CoordName = CoordName, Start = Start, End = End,
                        Strand = Strand, genome = human)
        CoordName, Start, End, Strand = '2', 2000, 2000000, 1
        c2 = Coordinate(CoordName = CoordName, Start = Start, End = End,
                        Strand = Strand, genome = human)
        c3 = c1.adopted(c2)
        self.assertEqual(c3.CoordName, c2.CoordName)
        self.assertEqual(c3.CoordType, c2.CoordType)
        self.assertEqual(c3.seq_region_id, c2.seq_region_id)
        self.assertEqual(c3.Start, c1.Start)
        self.assertEqual(c3.End, c1.End)
        self.assertEqual(c3.Strand, c1.Strand)
        c3 = c1.adopted(c2, shift = 100)
        self.assertEqual(c3.Start, c1.Start+100)
        self.assertEqual(c3.End, c1.End+100)
    

class TestCoordSystem(TestCase):
    def test_call(self):
        human_chrom = CoordSystem('chromosome', core_db = human.CoreDb,
                        species = 'human')
        human_contig = CoordSystem(1, species = 'human')
        self.assertEqual(human_chrom.coord_system_id, 2)
        self.assertEqual(human_contig.name, 'contig')
        self.assertEqual(human_contig.attr, 'default_version, sequence_level')
    

if __name__ == '__main__':
    main()
