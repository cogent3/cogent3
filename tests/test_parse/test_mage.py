#!/usr/bin/env python
"""Tests of MageParser
"""
from cogent.util.unit_test import TestCase, main
from cogent.parse.mage import MageParser, MageGroupFromString,\
    MageListFromString, MagePointFromString
from cogent.format.mage import MagePoint, MageList, MageGroup, Kinemage

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

class MageGroupFromStringTests(TestCase):
    """Tests for the MageGroupFromString function"""
    
    def test_MageGroupFromString(self):
        """MageGroupFromString should fill itself from string correctly"""
        f = MageGroupFromString
        group = f('@group {GA manifold}   off')
        self.assertEqual(str(group),'@group {GA manifold} off')
        group = f('@group {dna} recessiveon   off dominant    '+\
                'master= {master name} nobutton clone={clone_name}\n')
        self.assertEqual(str(group),\
            '@group {dna} off nobutton recessiveon dominant '+\
            'master={master name} clone={clone_name}')
        group = f(\
            '@subgroup {max_group} recessiveon instance=   {inst \tname} lens')
        self.assertEqual(str(group),
            '@subgroup {max_group} recessiveon lens instance={inst \tname}')
        group = f('@subgroup {Pos 1} on')
        self.assertEqual(str(group), '@subgroup {Pos 1}')
        
    def test_MageGroupFromString_wrong(self):
        """MageGroupFromString should fail on wrong input"""
        f = MageGroupFromString
        # unknown keyword
        self.assertRaises(KeyError,f,'@group {GA manifold}   of')
        # wrong nesting
        self.assertRaises(ValueError,f,'@group {GA manifold} master={blabla')
       
class MageListFromStringTests(TestCase):
    """Tests for the MageListFromString function"""

    def test_MageListFromString(self):
        """MageListFromString should fill itself from string correctly"""
        f = MageListFromString
        l = f('@dotlist {label} color= green radius=0.3 off \t nobutton\n' )
        self.assertEqual(str(l),\
                '@dotlist {label} off nobutton color=green radius=0.3')
        l = f('@vectorlist')
        self.assertEqual(str(l),'@vectorlist')
        l = f('@balllist off angle= 4 width= 2 face=something '+\
                'font=\tother size= 3')
        self.assertEqual(str(l),'@balllist off angle=4 width=2 '+\
                'face=something font=other size=3')
        l = f('@dotlist {} on nobutton color=sky')
        self.assertEqual(str(l),'@dotlist nobutton color=sky')
        
    def test_MageListFromString_wrong(self):
        """MageListFromString should fail on wrong input"""
        f = MageListFromString
        self.assertRaises(ValueError,f,
            '@somelist {label} color= green radius=0.3 off \t nobutton\n')
        self.assertRaises(KeyError,f,
            '@vectorlist {label} colors= green radius=0.3 off \t nobutton\n')
        
class MagePointFromStringTests(TestCase):
    """Tests of the MagePointFromString factory function."""

    def test_MagePointFromString(self):
        """MagePoint should fill itself from string correctly"""
        m = MagePointFromString('{construction}width5  0.000 0.707 -1.225\n')
        self.assertEqual(str(m), \
        '{construction} width5 ' + ' '.join(map(str, [0.0,0.707,-1.225])))
        m = MagePointFromString('3, 4, 5')
        self.assertEqual(str(m), ' '.join(map(str, map(float, [3, 4, 5]))))
        m = MagePointFromString('{b2}P 0.000 0.000 0.000')
        self.assertEqual(str(m), '{b2} P ' + \
            ' '.join(map(str, map(float, [0,0,0]))))
        m = MagePointFromString('P -2650192.000 4309510.000 3872241.000')
        self.assertEqual(str(m), 'P ' + \
            ' '.join(map(str, map(float, [-2650192,4309510,3872241]))))
        m = MagePointFromString('{"}P -2685992.000 5752262.000 535328.000')
        self.assertEqual(str(m), '{"} P ' + \
            ' '.join(map(str, map(float, [-2685992,5752262,535328]))))
        m = MagePointFromString('{ 1, d, 0       } P   1.000, 0.618, 0.000')
        self.assertEqual(str(m), '{ 1, d, 0       } P ' + \
            ' '.join(map(str, map(float, [1.000, 0.618, 0.000]))))
        m = MagePointFromString('{"}width1  -1.022 0.969 -0.131')
        self.assertEqual(str(m), '{"} width1 ' + \
            ' '.join(map(str, map(float, [-1.022,0.969,-0.131]))))
        m = MagePointFromString(\
            'width3 {A label with spaces } A blue r=3.7 5, 6, 7')
        self.assertEqual(m.Width, 3)
        self.assertEqual(m.Label, 'A label with spaces ')
        self.assertFloatEqual(m.Coordinates, [5, 6, 7])
        self.assertFloatEqual(m.Radius, 3.7)
        self.assertEqual(m.Color, 'blue')
        self.assertEqual(m.State, 'A')
        self.assertEqual(str(m),'{A label with spaces } A blue width3 r=3.7 ' +\
            ' '.join(map(str, map(float, [5, 6, 7]))))

class MageParserTests(TestCase):
    """Tests for the MageParser"""
    def test_MageParser(self):
        """MageParser should work on valid input"""
        obs = str(MageParser(EXAMPLE_1.split('\n'))).split('\n')
        exp = EXP_EXAMPLE_1.split('\n')
        assert len(obs) == len(exp)
        #first check per line; easier for debugging
        for x in range(len(obs)):
            self.assertEqual(obs[x],exp[x])
        #double check to see if the whole string is the same
        self.assertEqual(str(MageParser(EXAMPLE_1.split('\n'))),EXP_EXAMPLE_1)

EXAMPLE_1 = """
@text
Kinemage of ribosomal RNA SSU Bacteria
@kinemage1
@caption
SSU Bacteria secondary structure elements
@viewid {oblique}
@zoom 1.05
@zslab 467
@center 0.500 0.289 0.204
@matrix
-0.55836 -0.72046 -0.41133  0.82346 -0.42101 -0.38036  0.10085 -0.55108
0.82833

@2viewid {top}
@2zoom 0.82
@2zslab 470
@2center 0.500 0.289 0.204
@2matrix
-0.38337  0.43731 -0.81351  0.87217 -0.11840 -0.47466 -0.30389 -0.89148
-0.33602

@3viewid {side}
@3zoom 0.82
@3zslab 470
@3center 0.500 0.289 0.204
@3matrix
-0.49808 -0.81559 -0.29450  0.86714 -0.46911 -0.16738 -0.00164 -0.33875
0.94088

@4viewid {End-O-Line}
@4zoom 1.43
@4zslab 469
@4center 0.500 0.289 0.204
@4matrix
 0.00348 -0.99984 -0.01766  0.57533 -0.01244  0.81784 -0.81792 -0.01301
0.57519

@perspective
@fontsizelabel 24
@onewidth
@zclipoff
@localrotation  1 0 0 .5 .866 0 .5 .289 .816

@group {Tetrahedron}
@vectorlist  {Edges}  color=white     nobutton
P   {0 0 0} 0 0 0
     0.5 0 0
  {1 0 0} 1 0 0
     0.5  0.5  0
  {0 1 0} 0 1 0
P    0 0 0
    0  0.5  0
  {0 1 0} 0 1 0
   0  0.5  0.5
  {0 0 1} 0 0 1
P    0 0 0
   0 0 0.5
  {0 0 1} 0 0 1
   0.5  0 0.5
  {1 0 0} 1 0 0

@labellist {labels} color=white    nobutton
  {U}  0  0  0
   {A}  1.1    0  0
   {C}  0    1.05   0
   {G}   0  0   1.08


@group  {Lines}
@vectorlist {A=U&C=G} color= green    off
 P   0   0.5   0.5
     .1  .4  .4
     .25 .25 .25
     .4  .1  .1
 L   0.500, 0.000, 0.000

@vectorlist {A=G&C=U} color= red    off
 P   0.5   0   0.5
     .25 .25 .25
 L   0, 0.500, 0.000

@vectorlist {A=C&G=U} color= red    off
 P   0.5   0.5   0
     .25 .25 .25
 L   0.000, 0.000, 0.500

@group {SSU Bacteria} recessiveon 
@dotlist {Stem} radius=0.03 color= orange
{a} .3 .1 .4
{b} r=.2 .1 .1 .1
@balllist {Junction} radius=.04
{c} red .4 .4 0
{}\t r=.1           green  .3 .2 .1
@group {empty group}
@group {Nested}
@subgroup {First}
@spherelist color=\tpurple 
{e} .1 .1 .1
@subgroup {Second} master=\t {master name}  
@labellist {labels}
{U}  0  0  0
   {A}  1.1    0  0
   {C}  0    1.05   0
   {G}   0  0   1.08
"""

EXP_EXAMPLE_1 =\
"""@kinemage 1
@viewid {oblique}
@zoom 1.05
@zslab 467
@center 0.500 0.289 0.204
@matrix
-0.55836 -0.72046 -0.41133  0.82346 -0.42101 -0.38036  0.10085 -0.55108
0.82833
@2viewid {top}
@2zoom 0.82
@2zslab 470
@2center 0.500 0.289 0.204
@2matrix
-0.38337  0.43731 -0.81351  0.87217 -0.11840 -0.47466 -0.30389 -0.89148
-0.33602
@3viewid {side}
@3zoom 0.82
@3zslab 470
@3center 0.500 0.289 0.204
@3matrix
-0.49808 -0.81559 -0.29450  0.86714 -0.46911 -0.16738 -0.00164 -0.33875
0.94088
@4viewid {End-O-Line}
@4zoom 1.43
@4zslab 469
@4center 0.500 0.289 0.204
@4matrix
0.00348 -0.99984 -0.01766  0.57533 -0.01244  0.81784 -0.81792 -0.01301
0.57519
@perspective
@fontsizelabel 24
@onewidth
@zclipoff
@localrotation  1 0 0 .5 .866 0 .5 .289 .816
@text
Kinemage of ribosomal RNA SSU Bacteria
@caption
SSU Bacteria secondary structure elements
@group {Tetrahedron}
@vectorlist {Edges} nobutton color=white
{0 0 0} P 0.0 0.0 0.0
0.5 0.0 0.0
{1 0 0} 1.0 0.0 0.0
0.5 0.5 0.0
{0 1 0} 0.0 1.0 0.0
P 0.0 0.0 0.0
0.0 0.5 0.0
{0 1 0} 0.0 1.0 0.0
0.0 0.5 0.5
{0 0 1} 0.0 0.0 1.0
P 0.0 0.0 0.0
0.0 0.0 0.5
{0 0 1} 0.0 0.0 1.0
0.5 0.0 0.5
{1 0 0} 1.0 0.0 0.0
@labellist {labels} nobutton color=white
{U} 0.0 0.0 0.0
{A} 1.1 0.0 0.0
{C} 0.0 1.05 0.0
{G} 0.0 0.0 1.08
@group {Lines}
@vectorlist {A=U&C=G} off color=green
P 0.0 0.5 0.5
0.1 0.4 0.4
0.25 0.25 0.25
0.4 0.1 0.1
L 0.5 0.0 0.0
@vectorlist {A=G&C=U} off color=red
P 0.5 0.0 0.5
0.25 0.25 0.25
L 0.0 0.5 0.0
@vectorlist {A=C&G=U} off color=red
P 0.5 0.5 0.0
0.25 0.25 0.25
L 0.0 0.0 0.5
@group {SSU Bacteria} recessiveon
@dotlist {Stem} color=orange radius=0.03
{a} 0.3 0.1 0.4
{b} r=0.2 0.1 0.1 0.1
@balllist {Junction} radius=.04
{c} red 0.4 0.4 0.0
{} green r=0.1 0.3 0.2 0.1
@group {empty group}
@group {Nested}
@subgroup {First}
@spherelist color=purple
{e} 0.1 0.1 0.1
@subgroup {Second} master={master name}
@labellist {labels}
{U} 0.0 0.0 0.0
{A} 1.1 0.0 0.0
{C} 0.0 1.05 0.0
{G} 0.0 0.0 1.08"""

if __name__ == "__main__":
    main()
