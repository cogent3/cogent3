import os

from cogent import DNA
from cogent.util.unit_test import TestCase, main
from cogent.db.ensembl.host import HostAccount, get_ensembl_account
from cogent.db.ensembl.genome import Genome
from cogent.db.ensembl.assembly import CoordSystem, Coordinate, get_coord_conversion
from cogent.db.ensembl.feature_level import FeatureCoordLevels

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

Release = 76

if 'ENSEMBL_ACCOUNT' in os.environ:
    args = os.environ['ENSEMBL_ACCOUNT'].split()
    host, username, password = args[0:3]
    kwargs = {}
    if len(args) > 3:
        kwargs['port'] = int(args[3])
    account = HostAccount(host, username, password, **kwargs)
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
        self.assertEquals(chicken_feature_levels['repeat'].levels,
                                ['chromosome', 'scaffold'])
        self.assertEquals(set(chicken_feature_levels['cpg'].levels),
                            set(['chromosome', 'scaffold']))
    
    def test_repeat(self):
        # use chicken genome as it need to do conversion
        # chicken coordinate correspondent toRefSeq human IL2A region
        coord = dict(CoordName=9, Start=21727352, End=21729141)
        region = self.chicken.getRegion(**coord)
        # repeat is recorded at contig level, strand is 0
        repeats = region.getFeatures(feature_types = 'repeat')
        expect = [("9", 21727499, 21727527), ("9", 21728009, 21728018),
                  ("9", 21728169, 21728178)]
        obs = []
        for repeat in repeats:
            loc = repeat.Location
            obs.append((str(loc.CoordName), loc.Start, loc.End))
        self.assertEquals(set(obs), set(expect))
    
    def test_cpg(self):
        # contain 3 CpG island recorded at chromosome level
        coord1 = dict(CoordName=26, Start=105184, End=184346)
        cpgs1 = self.chicken.getFeatures(feature_types='cpg', **coord1)
        exp = [("26", 112153, 113139), ("26", 134125, 135050),
               ("26", 178899, 180227)]
        obs = []
        for cpg in cpgs1:
            loc = cpg.Location
            obs.append((str(loc.CoordName), loc.Start, loc.End))
        self.assertEquals(set(obs), set(exp))
        
        # test cpg features record at scaffold level:
        coord2 = dict(CoordName='JH376196.1', Start=1, End=14640)
        cpgs2 = self.chicken.getFeatures(feature_types='cpg', **coord2)
        self.assertEquals(len(list(cpgs2)), 3)
    

if __name__ == '__main__':
    main()
