import os

from cogent import DNA
from cogent.util.unit_test import TestCase, main
from cogent.db.ensembl.host import HostAccount, get_ensembl_account
from cogent.db.ensembl.genome import Genome
from cogent.db.ensembl.assembly import CoordSystem, Coordinate, get_coord_conversion
from cogent.db.ensembl.feature_level import FeatureCoordLevels

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

class TestFeatureCoordLevels(TestCase):
    def setUp(self):
        self.chicken = Genome(Species='chicken', Release=Release,
                            account=account)
    
    def test_feature_levels(self):
        ChickenFeatureLevels = FeatureCoordLevels('chicken')
        chicken_feature_levels = ChickenFeatureLevels(
                    feature_types=['gene', 'cpg', 'est'],
                    core_db=self.chicken.CoreDb,
                    otherfeature_db=self.chicken.OtherFeaturesDb)
        self.assertEquals(chicken_feature_levels['repeat'].levels, ['contig'])
        self.assertEquals(set(chicken_feature_levels['cpg'].levels),\
                            set(['contig', 'supercontig', 'chromosome']))
    
    def test_repeat(self):
        # use chicken genome as it need to do conversion
        # chicken coordinate correspondent toRefSeq human IL2A region
        coord = dict(CoordName=9, Start=23817146, End=23818935)
        region = self.chicken.getRegion(**coord)
        # repeat is recorded at contig level, strand is 0
        repeats = region.getFeatures(feature_types = 'repeat')
        expect = [("9", 23817293, 23817321), ("9", 23817803, 23817812),
                  ("9", 23817963, 23817972)]
        obs = []
        for repeat in repeats:
            loc = repeat.Location
            obs.append((loc.CoordName, loc.Start, loc.End))
        self.assertEquals(set(obs), set(expect))
    
    def test_cpg(self):
        # contain 3 CpG island recorded at chromosome level
        coord1 = dict(CoordName=26, Start=110000, End=190000)
        cpgs1 = self.chicken.getFeatures(feature_types = 'cpg', **coord1)
        exp = [("26", 116969, 117955), ("26", 139769, 140694),
               ("26", 184546, 185881)]
        obs = []
        for cpg in cpgs1:
            loc = cpg.Location
            obs.append((loc.CoordName, loc.Start, loc.End))
        self.assertEquals(set(exp), set(obs))
        
        # test cpg features record at supercontig level:
        coord2 = dict(CoordName='Un_random', Start=29434117, End=29439117)
        cpgs2 = self.chicken.getFeatures(feature_types='cpg', **coord2)
        self.assertEquals(len(list(cpgs2)), 1)
    

if __name__ == '__main__':
    main()
