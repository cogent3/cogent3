#!/usr/bin/env python
"""Unit tests for Mage format writer.
"""
from __future__ import division
from numpy import array
from copy import deepcopy
from cogent.util.unit_test import TestCase, main
from cogent.format.mage import MagePoint, MageList, MageGroup, MageHeader, \
    Kinemage, MagePointFromBaseFreqs
from cogent.core.usage import BaseUsage
from cogent.util.misc import Delegator

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Sandra Smit", "Gavin Huttley", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class MagePointTests(TestCase):
    """Tests of the MagePoint class, holding information about points."""

    def test_init_empty(self):
        """MagePoint should init correctly with no data"""
        m = MagePoint()
        self.assertEqual(str(m), ' '.join(map(str,[0,0,0])))

    def test_init(self):
        """MagePoint should init correctly with normal cases"""
        #label and coords
        m = MagePoint([0.200,0.000,0.800], '0.800')
        self.assertEqual(str(m), '{0.800} ' + \
            ' '.join(map(str, ([0.200,0.000,0.800]))))
        #coords only
        m = MagePoint([0.200,0.000,0.800])
        self.assertEqual(str(m), ' '.join(map(str, ([0.200,0.000,0.800]))))
        #label only
        m = MagePoint(Label='abc')
        self.assertEqual(str(m), '{abc} '+' '.join(map(str,[0,0,0])))
        #all fields occupied
        m = MagePoint(Label='abc', Coordinates=[3, 6, 1.5], Radius=0.5, \
            Width=2, State='P', Color='green')
        self.assertEqual(str(m), \
        '{abc} P green width2 r=0.5 3 6 1.5')

    def test_cmp(self):
        """MagePoint cmp should compare all fields"""
        self.assertEqual(MagePoint([0,0,0]), MagePoint([0,0,0]))
        self.assertNotEqual(MagePoint([0,0,0]), MagePoint([0,0,0], Color='red'))

    def test_get_coord(self):
        """MagePoint _get_coord should return coordinate that is asked for"""
        m = MagePoint([0,1,2])
        self.assertEqual(m.X,m.Coordinates[0])
        self.assertEqual(m.Y,m.Coordinates[1])
        self.assertEqual(m.Z,m.Coordinates[2])
        m = MagePoint()
        self.assertEqual(m.X,m.Coordinates[0])
        self.assertEqual(m.Y,m.Coordinates[1])
        self.assertEqual(m.Z,m.Coordinates[2])

    def test_set_coord(self):
        """MagePoint _get_coord should return coordinate that is asked for"""
        m = MagePoint([0,1,2])
        m.X, m.Y, m.Z = 2,3,4
        self.assertEqual(m.Coordinates,[2,3,4])
        m = MagePoint()
        m.X, m.Y, m.Z = 5,4,3
        self.assertEqual(m.Coordinates,[5,4,3])
        m = MagePoint()
        m.X = 5
        self.assertEqual(m.Coordinates,[5,0,0])

    def test_toCartesian(self):
        """MagePoint toCartesian() should transform coordinates correctly"""
        m = MagePoint([.1,.2,.3])
        self.assertEqual(m.toCartesian().Coordinates,[.6,.7,.5])
        m = MagePoint()
        self.assertEqual(m.toCartesian().Coordinates,[1,1,1])
        m = MagePoint([.25,.25,.25],Color='red',Label='label',State='L')
        self.assertEqual(m.toCartesian().Coordinates,[.5,.5,.5])
        self.assertEqual(m.toCartesian().Color,m.Color)
        self.assertEqual(m.toCartesian().Label,m.Label)
        self.assertEqual(m.toCartesian().State,m.State)
        m = MagePoint([1/3.0,1/3.0,0])
        self.assertFloatEqual(m.toCartesian().Coordinates,
                [2/3.0,1/3.0,2/3.0])
        m = MagePoint([1/3.0,1/3.0,1/3.0])
        self.assertFloatEqual(m.toCartesian().Coordinates,
                [1/3.0,1/3.0,1/3.0])
        m = MagePoint([3,4,5])
        self.assertRaises(ValueError,m.toCartesian)

    def test_fromCartesian(self):
        """MagePoint fromCartesian should transform coordinates correctly"""
        mp = MagePoint([2/3.0,1/3.0,2/3.0])
        self.assertFloatEqual(mp.fromCartesian().Coordinates,[1/3.0,1/3.0,0])
        points = [MagePoint([.1,.2,.3]),MagePoint([.25,.25,.25],Color='red',
                Label='label',State='L'),MagePoint([1/3,1/3,0]),
                MagePoint([0,0,0]),MagePoint([1/7,2/7,3/7])]
        for m in points:
            b = m.toCartesian().fromCartesian()
            self.assertFloatEqual(m.Coordinates,b.Coordinates)
            self.assertEqual(m.Color,b.Color)
            self.assertEqual(m.Label,b.Label)
            self.assertEqual(m.State,b.State)

            #even after multiple iterations?
            mutant = deepcopy(m)
            for x in range(10):
                mutant = mutant.toCartesian().fromCartesian()
            self.assertFloatEqual(m.Coordinates,mutant.Coordinates)

class freqs_label(dict):
    """dict with Label and Id, for testing MagePointFromBaseFreqs"""
    def __init__(self, Label, Id, freqs):
        self.Label = Label
        self.Id = Id
        self.update(freqs)

class freqs_display(dict):
    """dict with display properties, for testing MagePointFromBaseFreqs"""
    def __init__(self, Color, Radius, Id, freqs):
        self.Color = Color
        self.Radius = Radius
        self.Id = Id
        self.update(freqs)
         
class MagePointFromBaseFreqsTests(TestCase):
    """Tests of the MagePointFromBaseFreqs factory function."""
    def setUp(self):
        """Define a few standard frequencies"""
        self.empty = freqs_label(None, None, {})
        self.dna = freqs_label('dna', None, {'A':4, 'T':1, 'G':2, 'C':3})
        self.rna = freqs_label(None, 'rna', {'U':2, 'A':1, 'G':2})
        self.display = freqs_display('green', '0.25', 'xxx', {'A':2})

    def test_MagePointFromBaseFreqs(self):
        """MagePoint should fill itself from base freqs correctly"""
        e = MagePointFromBaseFreqs(self.empty)
        self.assertEqual(str(e), '0.0 0.0 0.0')
        dna = MagePointFromBaseFreqs(self.dna)
        self.assertEqual(str(dna), '{dna} 0.4 0.3 0.2')
        rna = MagePointFromBaseFreqs(self.rna)
        self.assertEqual(str(rna), '{rna} 0.2 0.0 0.4')
        display = MagePointFromBaseFreqs(self.display)
        self.assertEqual(str(display), \
            '{xxx} green r=0.25 1.0 0.0 0.0')

    def test_MagePointFromBaseFreqs_usage(self):
        """MagePoint should init correctly from base freqs"""
        class fake_seq(str, Delegator):
            def __new__(cls, data, *args):
                return str.__new__(cls, data)
            def __init__(self, data, *args):
                Delegator.__init__(self, *args)
                self.__dict__['Info'] = self._handler
                str.__init__(data)

        class has_species(object):
            def __init__(self, sp):
                self.Species = sp
            
        s = fake_seq('AAAAACCCTG', has_species('Homo sapiens'))
        b = BaseUsage(s)
        p = MagePointFromBaseFreqs(b)
        self.assertEqual(str(p), '{Homo sapiens} 0.5 0.3 0.1')

    def test_MagePointFromBaseFreqs_functions(self):
        """MagePointFromBaseFreqs should apply functions correctly"""
        def set_color(x):
            if x.Label == 'dna':
                return 'green'
            else:
                return 'blue'

        def set_radius(x):
            if x.Label == 'dna':
                return 0.25
            else:
                return 0.5

        def set_label(x):
            if x.Id is not None:
                return 'xxx'
            else:
                return 'yyy'

        self.assertEqual(str(MagePointFromBaseFreqs(self.dna, 
            get_label=set_label)),
            '{yyy} 0.4 0.3 0.2')
        self.assertEqual(str(MagePointFromBaseFreqs(self.rna, 
            get_label=set_label)),
            '{xxx} 0.2 0.0 0.4')
        self.assertEqual(str(MagePointFromBaseFreqs(self.dna,
            get_radius=set_radius)),
            '{dna} r=0.25 0.4 0.3 0.2')
        self.assertEqual(str(MagePointFromBaseFreqs(self.rna,
            get_radius=set_radius)),
            '{rna} r=0.5 0.2 0.0 0.4')
        self.assertEqual(str(MagePointFromBaseFreqs(self.dna, 
            get_color=set_color)),
            '{dna} green 0.4 0.3 0.2')
        self.assertEqual(str(MagePointFromBaseFreqs(self.rna, 
            get_color=set_color)),
            '{rna} blue 0.2 0.0 0.4')
        self.assertEqual(str(MagePointFromBaseFreqs(self.rna,
            get_label=set_label, get_radius=set_radius,get_color=set_color)),
            '{xxx} blue r=0.5 0.2 0.0 0.4')

class MageListTests(TestCase):
    """Tests of the MageList class, holding a collection of points."""
    def setUp(self):
        """Define a few standard points and lists of points."""
        self.null = MagePoint([0,0,0])
        self.label = MagePoint([1, 1, 1], 'test')
        self.properties = MagePoint(Width=1, Label='point', State='L',\
            Color='blue', Coordinates=[2.0,4.0,6.0])
        self.radius1 = MagePoint([2,2,2],Radius=.1)
        self.radius2 = MagePoint([3,3,3],Radius=.5)
        self.first_list = [self.null, self.properties]
        self.empty_list = []
        self.minimal_list = [self.null]
        self.single_list = [self.label]
        self.multi_list = [self.properties] * 10
        self.radii = [self.radius1,self.radius2]
        
    def test_init_empty(self):
        """MageList should init correctly with no data"""
        m = MageList()
        self.assertEqual(str(m), "@dotlist")
        m = MageList(self.empty_list)
        self.assertEqual(str(m), "@dotlist")

    def test_init(self):
        """MageList should init correctly with data"""
        m = MageList(self.minimal_list)
        self.assertEqual(str(m), "@dotlist\n" + str(self.null))
        m = MageList(self.multi_list,'x',Off=True,Color='green',NoButton=True)
        self.assertEqual(str(m), "@dotlist {x} off nobutton color=green\n" + \
            '\n'.join(10 * [str(self.properties)]))
        m = MageList(self.first_list,NoButton=True,Color='red', \
            Style='vector', Radius=0.03, Width=3, Label='test')
        self.assertEqual(str(m), "@vectorlist {test} nobutton color=red " + \
        "radius=0.03 width=3\n" + str(self.null) + '\n' + str(self.properties))

    def test_toArray_radii(self):
        """MageList toArray should return the correct array"""
        m = MageList(self.empty_list)
        self.assertEqual(m.toArray(),array(()))
        m = MageList(self.first_list,Radius=.3)
        self.assertEqual(m.toArray(),array([[0,0,0,0.3],[2.0,4.0,6.0,0.3]]))
        m = MageList(self.radii)
        self.assertEqual(m.toArray(), array([[2,2,2,.1],[3,3,3,.5]]))
        m = MageList(self.radii,Radius=.4)
        self.assertEqual(m.toArray(), array([[2,2,2,.1],[3,3,3,.5]]))
        m = MageList(self.single_list) #radius = None
        self.assertRaises(ValueError,m.toArray)
    
    def test_toArray_coords_only(self):
        """MageList toArray should return the correct array"""
        m = MageList(self.empty_list)
        self.assertEqual(m.toArray(include_radius=False),array(()))
        m = MageList(self.first_list,Radius=.3)
        self.assertEqual(m.toArray(include_radius=False),
                array([[0,0,0],[2.0,4.0,6.0]]))
        m = MageList(self.radii)
        self.assertEqual(m.toArray(include_radius=False), 
                array([[2,2,2],[3,3,3]]))
        m = MageList(self.radii,Radius=.4)
        self.assertEqual(m.toArray(include_radius=False),
                array([[2,2,2],[3,3,3]]))
        m = MageList(self.single_list) #radius = None
        self.assertEqual(m.toArray(include_radius=False),array([[1,1,1]]))

    def test_iterPoints(self):
        """MageList iterPoints should yield all points in self"""
        m = MageList(self.single_list)
        for point in m.iterPoints():
            assert isinstance(point,MagePoint)
        self.assertEqual(len(list(m.iterPoints())),1)
        m = MageList(self.multi_list)
        for point in m.iterPoints():
            assert isinstance(point,MagePoint)
        self.assertEqual(len(list(m.iterPoints())),10)

    def test_toCartesian(self):
        """MageList toCartesian should return new list"""
        m = MageList([self.null],Color='green')
        res = m.toCartesian()
        self.assertEqual(len(m), len(res))
        self.assertEqual(m.Color,res.Color)
        self.assertEqual(res[0].Coordinates,[1,1,1])
        m.Color='red'
        self.assertEqual(res.Color,'green')
        m = MageList([self.properties])
        self.assertRaises(ValueError,m.toCartesian)
    
    def test_fromCartesian(self):
        """MageList fromCartesian() should return new list with ACG coordinates
        """
        point = MagePoint([.1,.2,.3])
        m = MageList([point]*5,Color='green')
        res = m.toCartesian().fromCartesian()
        self.assertEqual(str(m),str(res))

class MageGroupTests(TestCase):
    """Test cases for the MageGroup class."""
    def setUp(self):
        """Define some standard lists and groups."""
        self.p1 = MagePoint([0, 1, 0], Color='green', Label='x')
        self.p0 = MagePoint([0,0,0])
        self.min_list = MageList([self.p0]*2,'y')
        self.max_list = MageList([self.p1]*5,'z',Color='blue',Off=True, \
            Style='ball')
        self.min_group = MageGroup([self.min_list], Label="min_group")
        self.max_group = MageGroup([self.min_list, self.max_list], Color='red',
                         Label="max_group", Style='dot')
        self.nested = MageGroup([self.min_group, self.max_group], Label='nest',
            Color='orange', Radius=0.3, Style='vector')
        self.empty = MageGroup(Label='empty',Color='orange', NoButton=True,
                Style='vector',RecessiveOn=False)

    def test_init(self):
        """Nested MageGroups should set subgroup and cascades correctly."""
        exp_lines = [
        '@group {nest} recessiveon',
        '@subgroup {min_group} recessiveon',
        '@vectorlist {y} color=orange radius=0.3',
        str(self.p0),
        str(self.p0),
        '@subgroup {max_group} recessiveon',
        '@dotlist {y} color=red radius=0.3',
        str(self.p0),
        str(self.p0),
        '@balllist {z} off color=blue radius=0.3',
        str(self.p1),
        str(self.p1),
        str(self.p1),
        str(self.p1),
        str(self.p1),
        ]
        s = str(self.nested).split('\n')
        self.assertEqual(str(self.nested), '\n'.join(exp_lines))
        #check that resetting the cascaded values works OK
        nested = self.nested
        str(nested)
        self.assertEqual(nested,self.nested)
        self.assertEqual(nested[0][0].Color,None)
    
    def test_str(self):
        """MageGroup str should print correctly"""
        m = self.empty
        self.assertEqual(str(self.empty),'@group {empty} nobutton')
        m = MageGroup(Label='label',Clone='clone_name',Off=True)
        self.assertEqual(str(m),
                '@group {label} off recessiveon clone={clone_name}')
        m = MageGroup()
        self.assertEqual(str(m),'@group recessiveon')
        
    def test_iterGroups(self):
        """MageGroup iterGroups should behave as expected"""
        groups = list(self.nested.iterGroups())
        self.assertEqual(groups[0],self.min_group)
        self.assertEqual(groups[1],self.max_group)
        self.assertEqual(len(groups),2)

    def test_iterLists(self):
        """MageGroup iterLists should behave as expected"""
        lists = list(self.nested.iterLists())
        self.assertEqual(len(lists),3)
        self.assertEqual(lists[0],self.min_list)
        self.assertEqual(lists[1],self.min_list)
        self.assertEqual(lists[2],self.max_list)

    def test_iterGroupsAndLists(self):
        """MageGroup iterGroupsAndLists should behave as expected"""
        all = list(self.nested.iterGroupsAndLists())
        self.assertEqual(len(all),5)
        self.assertEqual(all[0],self.min_group)
        self.assertEqual(all[4],self.max_list)

    def test_iterPoints(self):
        """MageGroup iterPoints should behave as expected"""
        points = list(self.nested.iterPoints())
        self.assertEqual(len(points),9)
        self.assertEqual(points[1],self.p0)
        self.assertEqual(points[6],self.p1)

    def test_toCartesian(self):
        """MageGroup toCartesian should return a new MageGroup"""
        m = self.nested
        res = m.toCartesian()
        self.assertEqual(len(m),len(res))
        self.assertEqual(m.RecessiveOn,res.RecessiveOn)
        self.assertEqual(m[1][1].Color, res[1][1].Color)
        self.assertEqual(res[1][1][1].Coordinates,[1,0,0])
    
    def test_fromCartesian(self):
        """MageGroup fromCartesian should return a new MageGroup"""
        point = MagePoint([.1,.2,.3])
        l = MageList([point]*5,Color='red')
        m = MageGroup([l],Radius=0.02,Subgroup=True)
        mg = MageGroup([m])
        res = mg.toCartesian().fromCartesian()
        self.assertEqual(str(mg),str(res))

class MageHeaderTests(TestCase):
    """Tests of the MageHeader class.
    
    For now, MageHeader does nothing, so just verify that it gets the string.
    """
    def test_init(self):
        """MageHeader should keep the string it was initialized with."""
        m = MageHeader('@perspective')
        self.assertEqual(str(m), '@perspective')

class KinemageTests(TestCase):
    """Tests of the overall Kinemage class."""
    def setUp(self):
        self.point = MagePoint([0,0,0],'x')
        self.ml = MageList([self.point], Label='y',Color='green')
        self.mg1 = MageGroup([self.ml],Label='z')
        self.mg2 = MageGroup([self.ml,self.ml],Label='b')
        self.kin = Kinemage(1)
        self.kin.Groups = [self.mg1,self.mg2]
    
    def test_init_empty(self):
        """Kinemage empty init should work, but refuse to print"""
        k = Kinemage()
        self.assertEqual(k.Count, None)
        self.assertRaises(ValueError, k.__str__)
                
    def test_init(self):
        """Kinemage should init with any of its usual fields"""
        k = Kinemage(1)
        self.assertEqual(str(k), '@kinemage 1')
        k.Header = '@perspective'
        self.assertEqual(str(k), '@kinemage 1\n@perspective')
        k.Count = 2
        self.assertEqual(str(k), '@kinemage 2\n@perspective')
        k.Header = ''
        k.Caption = 'test caption'
        self.assertEqual(str(k), '@kinemage 2\n@caption\ntest caption')
        k.Caption = None
        k.Text = 'some text'
        self.assertEqual(str(k), '@kinemage 2\n@text\nsome text')
        k.Groups = [self.mg1]
        k.Header = '@test_header'
        k.Caption = 'This is\nThe caption'
        k.Text = 'some text here'
        self.assertEqual(str(k), '@kinemage 2\n@test_header\n@text\n' +\
        'some text here\n' + \
        '@caption\nThis is\nThe caption\n@group {z} recessiveon\n' + \
        '@dotlist {y} color=green\n{x} 0 0 0')

    def test_iterGroups(self):
        """Kinemage iterGroups should behave as expected"""
        k = self.kin
        groups = list(k.iterGroups())
        self.assertEqual(len(groups),2)
        self.assertEqual(groups[0],self.mg1)
        self.assertEqual(groups[1],self.mg2)

    def test_iterLists(self):
        """Kinemage iterLists should behave as expected"""
        k = self.kin
        lists = list(k.iterLists())
        self.assertEqual(len(lists),3)
        self.assertEqual(lists[0],self.ml)
        
    def test_iterPoints(self):
        """Kinemage iterPoints should behave as expected"""
        k = self.kin
        points = list(k.iterPoints())
        self.assertEqual(len(points),3)
        self.assertEqual(points[0],self.point)

    def test_iterGroupAndLists(self):
        """Kinemage iterGroupsAndLists should behave as expected"""
        all = list(self.kin.iterGroupsAndLists())
        self.assertEqual(len(all),5)
        self.assertEqual(all[0],self.mg1)
        self.assertEqual(all[4],self.ml)
    
    def test_toCartesian(self):
        """Kinemage toCartesian should return new Kinemage with UC,UG,UA coords
        """
        k = self.kin
        res = k.toCartesian()
        self.assertEqual(len(k.Groups),len(res.Groups))
        self.assertEqual(k.Text,res.Text)
        self.assertEqual(k.Groups[1].RecessiveOn,res.Groups[1].RecessiveOn)
        self.assertEqual(res.Groups[0][0][0].Coordinates,[1,1,1])
    
    def test_fromCartesian(self):
        """Kinemage fromCartesian should return Kinemage with A,C,G(,U) coords
        """
        point = MagePoint([.1,.2,.3])
        l = MageList([point]*5,Color='red')
        m1 = MageGroup([l],Radius=0.02,Subgroup=True)
        m2 = MageGroup([l],Radius=0.02)
        mg = MageGroup([m1])
        k = Kinemage(Count=1,Groups=[mg,m2])
        res = k.toCartesian().fromCartesian()
        self.assertEqual(str(k),str(res))
    
if __name__ == '__main__':
    main()
