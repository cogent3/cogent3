#!/usr/bin/env python
#
# test_dihedral.py
#
# Tests the dihedral module.
#
"""Provides tests for functions in the file dihedral.py
"""

__author__ = "Kristian Rother"
__copyright__ = "Copyright 2007-2009 2008, The Cogent Project"
__contributors__ = ["Kristian Rother", "Sandra Smit"]
__credits__ = ["Janusz Bujnicki", "Nils Goldmann"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"

from cogent.util.unit_test import main, TestCase
from cogent.struct.dihedral import dihedral, scalar, angle, \
    DihedralGeometryError, AngleGeometryError
from random import random
from numpy import array
from math import pi, cos, sin

class DihedralTests(TestCase):

    def get_random_array(self):
        """
        Returns a one-dimensional numpy array with three random floats
        in the range between -5.0 and +5.0.
        """
        return array([(random()-0.5)*10,(random()-0.5)*10,(random()-0.5)*10])
        
    def assertAlmostEqualAngle(self, is_value, should_value, digits=7):
        """
        Checks two angles in degrees whether they are the same within
        the given number of digits. This has been implemented to make sure
        that 359.9999991 == 0.0
        """
        maxv = 359.0
        for i in range(digits):
            maxv += 0.9 * 0.1**i
        while is_value < 0.0: is_value += 360
        while is_value > maxv: is_value -= 360
        while should_value < 0.0: should_value += 360
        while should_value > maxv: should_value -= 360
        self.assertAlmostEqual(is_value,should_value, digits)

    def test_scalar(self):
        """Tests the scalar product function for one-dimensional arrays."""
        # test one-dimensional arrays
        self.assertEqual(scalar(array([0]),array([0])),0.0)
        self.assertEqual(scalar(array([2]),array([3])),6.0)
        self.assertEqual(scalar(array([0,0]),array([0,0])),0.0)
        self.assertEqual(scalar(array([-1,-4]),array([1,4])),-17.0)
        self.assertEqual(scalar(array([1,2]),array([3,4])),11.0)
        self.assertEqual(scalar(array([-1,-4]),array([-1,4])),-15.0)                
        self.assertEqual(scalar(array([0,0,0]),array([0,0,0])),0.0)
        self.assertEqual(scalar(array([2.5,0,-1]),array([2.5,0,-1])),7.25)
        self.assertEqual(scalar(array([1,2,3]),array([0,0,0])),0.0)
        self.assertEqual(scalar(array([1,2,3]),array([1,2,3])),14.0)
        self.assertEqual(scalar(array([1,2,3]),array([4,5,6])),32.0)        
        # test two-dimensional arrays (should not be a feature)
        self.assertNotEqual(scalar(array([[0,0],[0,0]]),\
            array([[0,0],[0,0]])),0.0)

    def test_angle_simple(self):
        """Tests the angle function for one- and two-dimensional vectors."""
        # test two-dimensional vectors (not arrays!)
        self.assertEqual(angle(array([0,1]),array([1,0])),0.5*pi)
        self.assertEqual(angle(array([5,0]),array([13,0])),0.0)
        self.assertEqual(angle(array([2,3]),array([26,39])),0.0)
        self.assertEqual(angle(array([2,3]),array([-3,2])),0.5*pi)        
        self.assertEqual(angle(array([-5,0]),array([13,0])),pi)
        # test three-dimensional vectors (not arrays!)        
        self.assertEqual(angle(array([0,0,-1]),array([0,0,1])),pi)
        self.assertEqual(angle(array([0,15,-1]),array([0,-15,1])),pi)
        self.assertEqual(angle(array([0,0,7]),array([14,14,0])),0.5*pi)
        self.assertEqual(angle(array([0,7,7]),array([0,14,14])),0.0)
        self.assertAlmostEqual(angle(array([100000000.0,0,1]),\
            array([1,0,0])),0.0)

    ##    def make_scipy_angles(self):
    ##        """Generates the test data given below. Was commented out 
    ##        for getting rid of the library dependency"""
    ##        from Scientific.Geometry import Vector
    ##        for i in range(20):
    ##            v1 = self.get_random_array()
    ##            v2 = self.get_random_array()
    ##            vec1 = Vector(v1[0],v1[1],v1[2])
    ##            vec2 = Vector(v2[0],v2[1],v2[2])
    ##            scipy_angle = vec1.angle(vec2)
    ##            out = [(vec1[0],vec1[1],vec1[2]),\
    ##                (vec2[0],vec2[1],vec2[2]),scipy_angle]
    ##            print out,","
    
    def test_angle_scipy(self):
        """
        Asserts that dihedral and ScientificPython calculate the same angles.
        """
        for v1,v2,scipy_angle in SCIPY_ANGLES:
            ang = angle(array(v1),array(v2))
            self.assertAlmostEqual(ang, scipy_angle)

    def test_angle_fail(self):
        """The angle function should fail for zero length vectors."""
        # should not work for zero length vectors
        self.assertRaises(AngleGeometryError,angle,\
            array([0,0]),array([0,0]))
        self.assertRaises(AngleGeometryError,angle,\
            array([0,0,0]),array([0,0,0]))

    def test_dihedral_eight_basic_directions(self):
        """Checks dihedrals in all 45 degree intervals."""
        # using vectors with integer positions.
        self.assertAlmostEqualAngle(\
            dihedral([-2,-1,0], [-1,0,0], [1,0,0], [2,-1, 0]),  0.0)
        self.assertAlmostEqualAngle(\
            dihedral([-2,-1,0], [-1,0,0], [1,0,0], [2,-1,-1]), 45.0)
        self.assertAlmostEqualAngle(\
            dihedral([-2,-1,0], [-1,0,0], [1,0,0], [2, 0,-1]), 90.0)
        self.assertAlmostEqualAngle(\
            dihedral([-2,-1,0], [-1,0,0], [1,0,0], [2, 1,-1]),135.0)
        self.assertAlmostEqualAngle(\
            dihedral([-2,-1,0], [-1,0,0], [1,0,0], [2, 1, 0]),180.0)
        self.assertAlmostEqualAngle(\
            dihedral([-2,-1,0], [-1,0,0], [1,0,0], [2, 1, 1]),225.0)
        self.assertAlmostEqualAngle(\
            dihedral([-2,-1,0], [-1,0,0], [1,0,0], [2, 0, 1]),270.0)
        self.assertAlmostEqualAngle(\
            dihedral([-2,-1,0], [-1,0,0], [1,0,0], [2,-1, 1]),315.0)

    def test_dihedral_rotation(self):
        """Checks all angles in 0.2 degree intervals."""
        # constructs vectors using sin/cos and then calculates dihedrals
        precision = 5.0 # the higher the better
        v1 = array([1,0,1])
        v2 = array([0,0,1])
        v3 = array([0,0,2])
        for i in range(int(360*precision)):
            degrees = i/precision
            radians = pi*degrees/180.0
            opp_degrees = 360-degrees
            if opp_degrees == 360.0: opp_degrees = 0.0
            # construct circular motion of vector
            v4 = array([cos(radians), sin(radians),2])
            self.assertAlmostEqualAngle(dihedral(v4,v3,v2,v1), degrees, 5)
            # check rotation in the opposite direction
            self.assertAlmostEqualAngle(dihedral(v1,v2,v3,v4), degrees, 5)


    def test_dihedral_samples(self):
        """Checks values measured manually from atoms in PyMOL."""
        coordinates = [
            [(-1.225,4.621,42.070),(-1.407,4.455,43.516),\
             (-2.495,4.892,44.221),(-3.587,5.523,43.715)],
            [(-2.495,4.892,44.221),(1.513,0.381,40.711),\
             (-3.091,4.715,47.723),(-0.567,3.892,44.433)],
            [(-0.349,5.577,39.446),(-1.559,3.400,41.427),\
             (-4.304,5.563,45.998),(-2.495,4.892,44.221)],
            [(-45.819,84.315,19.372),(-31.124,72.286,14.035),\
             (-27.975,58.688,7.025),(-16.238,78.659,23.731)],
            [(-29.346,66.973,24.152),(-29.977,69.635,24.580),\
             (-30.875,68.788,24.663),(-30.668,67.495,24.449)],
            [(-34.586,84.884,14.064),(-23.351,69.756,11.028),\
             (-40.924,69.442,24.630),(-30.875,68.788,24.663)]
        ]
        angles = [1.201, 304.621, 295.672, 195.184, 358.699, 246.603]
        for i in range(len(coordinates)):
            v1,v2,v3,v4 = coordinates[i]
            self.assertAlmostEqualAngle(dihedral(v1,v2,v3,v4), angles[i],3)

            
    def test_dihedral_linear(self):
        """The dihedral function should fail for collinear vectors."""
        v1 = [1,0,0]
        v2 = [2,0,0]
        v3 = [3,0,0]
        v4 = [4,0,0]
        # print dihedral(v1,v2,v3,v4)
        for i in range(100):
            offset = array([int((random()-0.5)*10),\
                            int((random()-0.5)*10),\
                            int((random()-0.5)*10)])
            v1 = array([int((random()-0.5)*100),\
                        int((random()-0.5)*100),\
                        int((random()-0.5)*100)])
            v2 = v1 * int((random()-0.5)*100) + offset
            v3 = v1 * int((random()-0.5)*100) + offset 
            v4 = v1 * int((random()-0.5)*100) + offset
            v1 += offset
            self.assertRaises(DihedralGeometryError,dihedral,v1,v2,v3,v4)

    def test_dihedral_identical(self):
        """The dihedral function should fail if two vectors are the same."""
        # except for the first and last (the vectors form a triangle),
        # in which case the dihedral angle should be 0.0
        for i in range(100):
            v1 = self.get_random_array()
            v2 = self.get_random_array()
            v3 = self.get_random_array()
            self.assertRaises(DihedralGeometryError,dihedral,v1,v1,v2,v3)
            self.assertRaises(DihedralGeometryError,dihedral,v1,v2,v1,v3)
            self.assertRaises(DihedralGeometryError,dihedral,v1,v2,v2,v3)
            self.assertRaises(DihedralGeometryError,dihedral,v1,v2,v3,v3)
            self.assertRaises(DihedralGeometryError,dihedral,v1,v3,v2,v3)
            # now the triangular case
            # make sure that 359.999998 is equal to 0.0
            torsion = dihedral(v1,v2,v3,v1) + 0.000001
            if torsion > 360.0: torsion -= 360.0
            self.assertAlmostEqualAngle(torsion,0.0,5)  

SCIPY_ANGLES = [
    [(-4.4891521637990852, -1.2310927013330153, -0.96969777583098771), 
     (4.2147455310344171, -3.5069051036633514, 2.2430685816310305),
     2.2088870817461759] ,
    [(0.13959847081794097, 1.7204537912940399, -1.9303780516641089), 
     (0.35412687539602361, -2.9493521724340743, -4.865941405480644), 
     1.2704043143950585] ,
    [(2.3192363837822327, -3.6376441859213848, -2.2337816400479813),
    (-1.0271253661119029, -2.5736009846920425, -4.1470855710278975),
    0.83609068310373857] ,
    [(-0.38347986357358477, -4.1453876196041719, -2.1354583394773785),
    (0.27416747968044608, -2.5747732838982551, -0.68554680652905264),
    0.28348352552436806] ,
    [(-2.4928231204363449, 1.9263125976608209, -0.34275964486924715),
    (0.6721152528064811, 1.5270465172130598, -3.5720133753579564),
    1.3701121510959966] ,
    [(-1.50101403139692, 2.3218275982958292, 1.044582416480222),
    (-3.044743729573085, 2.0655933798619532, 2.9037849925327897),
    0.46218988498537578] ,
    [(4.9648826388927603, -1.7439743270051977, 1.0432135258334796),
    (-3.3694557299188608, 3.7697639370274052, 2.6962018714965055),
    2.3004031013653625] ,
    [(-3.3337033325729051, -0.79660906888508021, -3.4875326261817454),
    (1.4735023133671066, -0.066399047153666846, 0.94171530437632489),
    2.8293957790283595] ,
    [(-2.1249404252000517, 2.7456001658201568, 1.6891202129451799),
    (-0.66412553435299504, 3.371012200444512, -1.1548086037901306),
    0.89846374464990042] ,
    [(3.3993205618602018, 1.2047703532166887, -1.5839949555795063),
    (-4.6759756026580863, -2.8551222890449646, 4.888270217692785),
    2.7825564894291754] ,
    [(-3.966467296275785, 0.75617096138383189, -3.1711352932360248),
    (2.1054362220912326, 4.2867761689586601, -0.65739369331424213),
    1.6933117193742961] ,
    [(0.44413554305522851, 4.6000690382282361, -3.8338383756621819),
    (-2.4947413565865029, 1.8136080147734013, -4.0295344084655405),
    0.73110709489481174] ,
    [(-1.0971777991639065, -3.3166205797568815, 2.7098739534055563),
    (2.1536566381847289, 4.7817155120086055, -0.068554664323454695),
    2.4878958372925202] ,
    [(2.2696914760438136, 4.8841630875833673, 4.9524177412608861),
    (1.2249510822623111, 0.73008672334971658, 1.131607772478449),
    0.45741431674591532] ,
    [(-2.4456899797216938, 4.7894200033986447, -2.839449354837468),
    (0.95035116225980154, 4.2179212878828238, -1.801158217734109),
    0.63144509227951684] ,
    [(-4.6954041297179474, -2.8326266911591391, -1.1804869511610427),
    (3.2585456362256924, -2.2325171051479265, -2.0527260317826901),
    1.8363466110668369] ,
    [(2.1416146613604283, 3.8577375591718677, 3.1463493245087939),
    (-0.32185887240442468, 2.2163051363839505, 2.4704882512058224),
    0.52534998201320993] ,
    [(-2.8493351354335941, -3.8203784990110954, 2.4657357720402273),
    (-2.7799043389229383, 4.358406526726669, -2.8319872383058744),
    2.0906833125217235] ,
    [(2.274223250163784, -3.6086250253596406, 1.7143006579401876),
    (3.2763334328544347, -0.89908959703552171, -4.4068824993431557),
    1.4477009361020545] ,
    [(0.66737672421842809, -3.4628508908383848, 3.9044108358095366),
    (-1.9078974719893915, -0.53231141116878433, 1.3323584972786728),
    1.0932781951137689] ,
    ]

if __name__== '__main__':
    main()

