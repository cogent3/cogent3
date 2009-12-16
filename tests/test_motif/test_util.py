#!/usr/bin/env python
#file cogent_tests/motif/test_util.py
from __future__ import division
from cogent.util.unit_test import TestCase, main
from cogent.motif.util import Location, ModuleInstance, Module, Motif,\
    MotifResults, MotifFormatter, html_color_to_rgb
from cogent.core.moltype import ASCII

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"

class LocationTests(TestCase):
    """Tests of Location class for holding module location.
    """
    def setUp(self):
        """Setup for Location tests."""
        self.location_no_end = Location('seq1',1)
        self.locations = [
                        Location('seq1',1,5),
                        Location('seq2',3,54),
                        Location('seq1',5,3),
                        Location('seq1',2,3),
                        Location('seq2',54,2),
                        Location('seq0',1,3),
                        ]
        self.locations_sorted = [
                        Location('seq0',1,3),
                        Location('seq1',1,5),
                        Location('seq1',5,3),
                        Location('seq1',2,3),
                        Location('seq2',3,54),
                        Location('seq2',54,2),
                        ]
        
    def test_init_no_end(self):
        """__init__ should properly initialize Location object"""
        self.assertEqual(self.location_no_end.SeqId, 'seq1')
        self.assertEqual(self.location_no_end.Start, 1)
        self.assertEqual(self.location_no_end.End, 2)

    def test_init_complete(self):
        """__init__ should properly initialize Location object"""
        self.assertEqual(self.locations[0].SeqId, 'seq1')
        self.assertEqual(self.locations[0].Start, 1)
        self.assertEqual(self.locations[0].End, 5)

    def test_cmp(self):
        """Location object should sort properly with __cmp__ overwritten."""
        self.locations.sort()
        self.assertEqual(self.locations, self.locations_sorted)

class ModuleInstanceTests(TestCase):
    """Tests for ModuleInstance class."""

    def setUp(self):
        """Setup function for ModuleInstance tests."""
        self.sequences = [
                        'accucua',
                        'caucguu',
                        'accucua',
                        'cgacucg',
                        'cgaucag',
                        'cuguacc',
                        'cgcauca',
                        ]
        self.locations = [
                        Location('seq0',1,3),
                        Location('seq1',2,3),
                        Location('seq1',1,5),
                        Location('seq1',5,3),
                        Location('seq2',3,54),
                        Location('seq2',54,2),
                        Location('seq3',4,0),
                        ]
        self.Pvalues = [
                        .1,
                        .002,
                        .0000000003,
                        .6,
                        .0094,
                        .6,
                        .00201,
                        ]
        self.Evalues = [
                        .006,
                        .02,
                        .9,
                        .0200000001,
                        .09,
                        .0000003,
                        .900001,
                        ]
        self.modules_no_e = []
        for i in xrange(7):
            self.modules_no_e.append(ModuleInstance(self.sequences[i],
                                                    self.locations[i],
                                                    self.Pvalues[i]))
        self.modules_p_and_e = []
        for i in xrange(7):
            self.modules_p_and_e.append(ModuleInstance(self.sequences[i],
                                                       self.locations[i],
                                                       self.Pvalues[i],
                                                       self.Evalues[i]))
        self.modules_no_e_sorted = [
            ModuleInstance(self.sequences[2],self.locations[2],self.Pvalues[2]),
            ModuleInstance(self.sequences[1],self.locations[1],self.Pvalues[1]),
            ModuleInstance(self.sequences[6],self.locations[6],self.Pvalues[6]),
            ModuleInstance(self.sequences[4],self.locations[4],self.Pvalues[4]),
            ModuleInstance(self.sequences[0],self.locations[0],self.Pvalues[0]),
            ModuleInstance(self.sequences[3],self.locations[3],self.Pvalues[3]),
            ModuleInstance(self.sequences[5],self.locations[5],self.Pvalues[5]),
            ]

        self.modules_p_and_e_sorted = [
            ModuleInstance(self.sequences[2],self.locations[2],self.Pvalues[2]),
            ModuleInstance(self.sequences[1],self.locations[1],self.Pvalues[1]),
            ModuleInstance(self.sequences[6],self.locations[6],self.Pvalues[6]),
            ModuleInstance(self.sequences[4],self.locations[4],self.Pvalues[4]),
            ModuleInstance(self.sequences[0],self.locations[0],self.Pvalues[0]),
            ModuleInstance(self.sequences[5],self.locations[5],self.Pvalues[5]),
            ModuleInstance(self.sequences[3],self.locations[3],self.Pvalues[3]),
            ]

    def test_init_no_p_e_values(self):
        """Init should properly initialize ModuleInstance objects."""
        module1 = ModuleInstance(self.sequences[0], self.locations[0])
        module2 = ModuleInstance(self.sequences[1], self.locations[1])
        self.assertEqual(module1.Sequence, 'accucua')
        self.assertEqual(module1.Location.SeqId, 'seq0')
        self.assertEqual(module2.Sequence, 'caucguu')
        self.assertEqual(module2.Location.SeqId, 'seq1')

    def test_init_no_e_values(self):
        """Init should properly initialize ModuleInstance objects."""
        self.modules_no_e.sort()
        self.assertEqual(self.modules_no_e, self.modules_no_e_sorted)

    def test_len(self):
        """len() should return correct length of the ModuleInstance sequence."""
        for module in self.modules_no_e:
            self.assertEqual(len(module), 7)

    def test_str(self):
        """str() should return the correct string for each ModuleInstance."""
        for module, seq in zip(self.modules_no_e, self.sequences):
            self.assertEqual(str(module), seq)

    def test_cmp(self):
        """ModuleInstances should sort properly with __cmp__ overwritten."""
        self.modules_no_e.sort()
        self.modules_p_and_e.sort()
        self.assertEqual(map(str,self.modules_no_e),
                         map(str,self.modules_no_e_sorted))
        self.assertEqual(map(str,self.modules_p_and_e),
                         map(str,self.modules_p_and_e_sorted))

class ModuleTests(TestCase):
    """Tests for Module class."""

    def setUp(self):
        """SetUp for Module class tests."""
        self.sequences = [
                        'accucua',
                        'caucguu',
                        'accucua',
                        'cgacucg',
                        'cgaucag',
                        'cuguacc',
                        'cgcauca',
                        ]
        self.locations = [
                        Location('seq0',1,3),
                        Location('seq1',2,3),
                        Location('seq1',1,5),
                        Location('seq1',5,3),
                        Location('seq2',3,54),
                        Location('seq2',54,2),
                        Location('seq3',4,0),
                        ]
        self.Pvalues = [
                        .1,
                        .002,
                        .0000000003,
                        .6,
                        .0094,
                        .6,
                        .00201,
                        ]
        self.Evalues = [
                        .006,
                        .02,
                        .9,
                        .0200000001,
                        .09,
                        .0000003,
                        .900001,
                        ]
        self.modules_no_e = []
        for i in xrange(7):
            self.modules_no_e.append(ModuleInstance(self.sequences[i],
                                                    self.locations[i],
                                                    self.Pvalues[i]))
        
        self.modules_p_and_e = []
        for i in xrange(7):
            self.modules_p_and_e.append(ModuleInstance(self.sequences[i],
                                                       self.locations[i],
                                                       self.Pvalues[i],
                                                       self.Evalues[i]))
        self.module_no_template = Module(
            {
                (self.modules_no_e[0].Location.SeqId,
                 self.modules_no_e[0].Location.Start):self.modules_no_e[0],
                (self.modules_no_e[1].Location.SeqId,
                 self.modules_no_e[1].Location.Start):self.modules_no_e[1],
                (self.modules_no_e[2].Location.SeqId,
                 self.modules_no_e[2].Location.Start):self.modules_no_e[2],
                (self.modules_no_e[3].Location.SeqId,
                 self.modules_no_e[3].Location.Start):self.modules_no_e[3],
                (self.modules_no_e[4].Location.SeqId,
                 self.modules_no_e[4].Location.Start):self.modules_no_e[4],
                (self.modules_no_e[5].Location.SeqId,
                 self.modules_no_e[5].Location.Start):self.modules_no_e[5],
                (self.modules_no_e[6].Location.SeqId,
                 self.modules_no_e[6].Location.Start):self.modules_no_e[6],
                }
            )

        self.module_with_template = Module(
            {
                (self.modules_no_e[0].Location.SeqId,
                 self.modules_no_e[0].Location.Start):self.modules_no_e[0],
                (self.modules_no_e[1].Location.SeqId,
                 self.modules_no_e[1].Location.Start):self.modules_no_e[1],
                (self.modules_no_e[2].Location.SeqId,
                 self.modules_no_e[2].Location.Start):self.modules_no_e[2],
                (self.modules_no_e[3].Location.SeqId,
                 self.modules_no_e[3].Location.Start):self.modules_no_e[3],
                (self.modules_no_e[4].Location.SeqId,
                 self.modules_no_e[4].Location.Start):self.modules_no_e[4],
                (self.modules_no_e[5].Location.SeqId,
                 self.modules_no_e[5].Location.Start):self.modules_no_e[5],
                (self.modules_no_e[6].Location.SeqId,
                 self.modules_no_e[6].Location.Start):self.modules_no_e[6],
                },
            Template = 'accgucg'
            )
        
    def test_init(self):
        """Init should properly initialize Module object."""
        module = Module(data={(self.modules_no_e[0].Location.SeqId,
                           self.modules_no_e[0].Location.Start): \
                          self.modules_no_e[0]})
        self.assertEqual(module.Template, None)
        self.assertEqual(module.Alphabet, ASCII.Alphabet)
        self.assertEqual(module.Pvalue, None)
        self.assertEqual(module.Evalue, None)
        self.assertEqual(module.keys(),[('seq0',1)])
        self.assertEqual(module.values(),[ModuleInstance(self.sequences[0],
                                                self.locations[0],
                                                self.Pvalues[0])])


    def test_cmp(self):
        """Module objects should sort properly with __cmp__ overwritten."""
        pvals_sorted = [3e-010, 0.002,
                        0.0020100000000000001,
                        0.0094000000000000004,
                        0.10000000000000001,
                        0.59999999999999998,
                        0.59999999999999998]
        evals_sorted = [.9,
                        .02,
                        .900001,
                        .09,
                        .006,
                        .0000003,
                        .0200000001,
                        ]
        modules = []
        for instance, pvalue, evalue in zip(self.modules_no_e,
                                            self.Pvalues,
                                            self.Evalues):
            modules.append(Module({(instance.Location.SeqId,
                                   instance.Location.Start):instance},
                                  Pvalue=pvalue,
                                  Evalue=evalue))
        modules.sort()
        for ans, p, e in zip(modules, pvals_sorted, evals_sorted):
            self.assertEqual(ans.Pvalue, p)
            self.assertEqual(ans.Evalue, e)

    def test_LocationDict(self):
        """LocationDict should return correct dictionary of locations."""
        location_dict_ans = {
                            'seq0':[1],
                            'seq1':[1,2,3],
                            'seq2':[2,3],
                            'seq3':[0],
                            }
        location_dict = self.module_no_template.LocationDict
        self.assertEqual(location_dict,location_dict_ans)


class MotifTests(TestCase):
    """Tests for Motif class."""

    def test_init(self):
        """Init should properly initialize Motif object."""
        module = Module({
                            ('a',3): ModuleInstance('guc', Location('a',3,5)),
                            ('b',3): ModuleInstance('guc', Location('b',3,5)),
                            ('c',8): ModuleInstance('guc', Location('c',8,10)),
                            })
        m = Motif(module)
        self.assertEqual(m.Modules,[module])
        self.assertEqual(m.Info,None)

class MotifResultsTests(TestCase):
    """Tests for MotifResults class."""

    def test_init(self):
        """Init should properly initialize MotifResults object."""
        module = Module({
                            ('a',3): ModuleInstance('guc', Location('a',3,5)),
                            ('b',3): ModuleInstance('guc', Location('b',3,5)),
                            ('c',8): ModuleInstance('guc', Location('c',8,10)),
                            })
        motif = Motif([module])
        results = {'key1':'value1','key2':'value2'}
        parameters = {'parameter1':1,'parameter2':2}
        mr = MotifResults([module],[motif],results,parameters)
        self.assertEqual(mr.Modules,[module])
        self.assertEqual(mr.Motifs,[motif])
        self.assertEqual(mr.Results,results)
        self.assertEqual(mr.parameter1,1)
        self.assertEqual(mr.parameter2,2)

class UtilTests(TestCase):
    """Tests for utility functions."""
    
    def test_html_color_to_rgb(self):
        """Tests for html_to_color_rgb."""
        html_colors = ['#FF0000','#00FF00','#0000FF','545454']
        rgb_colors = [(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0),\
            (0.32941176470588235,0.32941176470588235,0.32941176470588235)]
        for html, rgb in zip(html_colors, rgb_colors):
            self.assertEqual(html_color_to_rgb(html),rgb)
    
    def test_make_remap_dict(self):
        """Tests for make_remap_dict."""
        pass

class MotifFormatterTests(TestCase):
    """Tests for MotifFormatter class."""
    
    def setUp(self):
        """SetUp for MotifFormatter class tests."""
        self.sequences = [
                        'accucua',
                        'caucguu',
                        'accucua',
                        'cgacucg',
                        'cgaucag',
                        'cuguacc',
                        'cgcauca',
                        ]
        self.locations = [
                        Location('seq0',1,3),
                        Location('seq1',2,3),
                        Location('seq1',1,5),
                        Location('seq1',5,3),
                        Location('seq2',3,54),
                        Location('seq2',54,2),
                        Location('seq3',4,0),
                        ]
        self.Pvalues = [
                        .1,
                        .002,
                        .0000000003,
                        .6,
                        .0094,
                        .6,
                        .00201,
                        ]
        self.Evalues = [
                        .006,
                        .02,
                        .9,
                        .0200000001,
                        .09,
                        .0000003,
                        .900001,
                        ]
        self.modules_no_e = []
        for i in xrange(7):
            self.modules_no_e.append(ModuleInstance(self.sequences[i],
                                                    self.locations[i],
                                                    self.Pvalues[i]))

        self.module_with_template = Module(
            {
                (self.modules_no_e[0].Location.SeqId,
                 self.modules_no_e[0].Location.Start):self.modules_no_e[0],
                (self.modules_no_e[1].Location.SeqId,
                 self.modules_no_e[1].Location.Start):self.modules_no_e[1],
                (self.modules_no_e[2].Location.SeqId,
                 self.modules_no_e[2].Location.Start):self.modules_no_e[2],
                (self.modules_no_e[3].Location.SeqId,
                 self.modules_no_e[3].Location.Start):self.modules_no_e[3],
                (self.modules_no_e[4].Location.SeqId,
                 self.modules_no_e[4].Location.Start):self.modules_no_e[4],
                (self.modules_no_e[5].Location.SeqId,
                 self.modules_no_e[5].Location.Start):self.modules_no_e[5],
                (self.modules_no_e[6].Location.SeqId,
                 self.modules_no_e[6].Location.Start):self.modules_no_e[6],
                },
            Template = 'accgucg', ID='1'
            )
        
        self.modules_with_ids =\
                    [Module({
                            ('a',3): ModuleInstance('guc', Location('a',3,5)),
                            ('b',3): ModuleInstance('guc', Location('b',3,5)),
                            ('c',8): ModuleInstance('guc', Location('c',8,10)),
                            },ID='1'),
                    Module({
                            ('a',7): ModuleInstance('cca', Location('a',7,9)),
                            ('b',7): ModuleInstance('cca', Location('b',7,9)),
                            ('c',11): ModuleInstance('cca',Location('c',11,13)),
                            },ID='2'),
                    Module({
                            ('a',10): ModuleInstance('gca',Location('a',10,12)),
                            ('b',10): ModuleInstance('gca',Location('b',10,12)),
                            ('c',14): ModuleInstance('gca',Location('c',14,12)),
                            },ID='3'),
                    Module({
                            ('a',13): ModuleInstance('ggg',Location('a',13,15)),
                            ('b',13): ModuleInstance('ggg',Location('b',13,15)),
                            ('c',18): ModuleInstance('ggg',Location('c',18,20)),
                            },ID='4'),
                    ]
        self.motifs_with_ids = map(Motif,self.modules_with_ids)
        self.motif_results = MotifResults(Modules=self.modules_with_ids,\
            Motifs=self.motifs_with_ids)
        
        self.color_map = {'1':"""background-color: #0000FF; ; font-family: 'Courier New', Courier""",
                          '2':"""background-color: #FFFF00; ; font-family: 'Courier New', Courier""",
                          '3':"""background-color: #00FFFF; ; font-family: 'Courier New', Courier""",
                          '4':"""background-color: #FF00FF; ; font-family: 'Courier New', Courier""",
                          }
        self.color_map_rgb = {
            'color_1':(0.0,0.0,1.0),
            'color_2':(1.0,1.0,0.0),
            'color_3':(0.0,1.0,1.0),
            'color_4':(1.0,0.0,1.0),
            }
        
    def test_getColorMapS0(self):
        """tests for getColorMapS0"""
        mf = MotifFormatter()
        module_ids = ['1','2','3','4']
        self.assertEqual(mf.getColorMapS0(module_ids),self.color_map)
    
    def test_getColorMap(self):
        """tests for getColorMap"""
        mf = MotifFormatter()
        self.assertEqual(mf.getColorMap(self.motif_results),self.color_map)
    
    def test_getColorMapRgb(self):
        """tests for getColorMapRgb"""
        mf = MotifFormatter()
        self.assertEqual(mf.getColorMapRgb(self.motif_results),\
            self.color_map_rgb)


#run if called from command-line
if __name__ == "__main__":
    main()
