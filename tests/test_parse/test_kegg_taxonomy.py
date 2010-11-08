#!/usr/bin/env python
from cogent.util.unit_test import TestCase, main

from cogent.parse.kegg_ko import parse_kegg_taxonomy

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jesse Zaneveld","Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Production"

"""
Test code for kegg_taxonomy.py in cogent.parse.  
"""

class ParseKEGGTaxonomy(TestCase):
    
    def test_parse_kegg_taxonomy(self):
        """parse_kegg_taxonomy should return successive taxonomy entries from
        lines"""
        
        test_lines =\
        ['# Eukaryotes\n',\
         '## Animals\n',\
         '### Vertebrates\n',\
         '#### Mammals\n',\
         'T01001(2000)\thsa\tH.sapiens\tHomo sapiens (human)\n',\
         '#### Birds\n',\
         'T01006(2005)\tgga\tG.gallus\tGallus gallus (chicken)\n',\
         '### Arthropods\n',\
         '#### Insects\n',\
         'T00030(2000)\tdme\tD.melanogaster\tDrosophila melanogaster (fruit fly)\n']
       
        exp =\
        ['Eukaryotes\tAnimals\tVertebrates\tMammals\tT01001(2000)\thsa\tH.sapiens\tHomo sapiens (human)\tHomo\tsapiens\thuman\n',\
         'Eukaryotes\tAnimals\tVertebrates\tBirds\tT01006(2005)\tgga\tG.gallus\tGallus gallus (chicken)\tGallus\tgallus\tchicken\n',\
         'Eukaryotes\tAnimals\tArthropods\tInsects\tT00030(2000)\tdme\tD.melanogaster\tDrosophila melanogaster (fruit fly)\tDrosophila\tmelanogaster\tfruit fly\n']
        obs = parse_kegg_taxonomy(test_lines)
        for i,res in enumerate(obs):
            self.assertEqual(res,exp[i])
                
if __name__=="__main__":
    main()
