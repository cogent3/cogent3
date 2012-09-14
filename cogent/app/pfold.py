#!/usr/bin/env python

"""Application controllers for pfold application package

Run in the same order as the order in this file(instructions from pfold author)
[fasta2col,findphyl,mltree,scfg]
"""
import os

from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"  

#IMPORTANT!!!!
#pfold_path must be set manually to the bin dir in the install dir of pfold
if 'PFOLD_BIN_DIR' not in os.environ:
    raise RuntimeError, \
        "The pfold app controller requires PFOLD_BIN_DIR environment variable"
else:
    pfold_path = os.environ['PFOLD_BIN_DIR']

class fasta2col(CommandLineApplication):
    """Application controller for fasta2col in Pfold package"""

    _command = 'fasta2col'
    _input_handler = '_input_as_string'

    def _input_as_string(self,filename):
        """Overrides  _input_as_string in CommandLineApplication

        Input file need to be modified by sed"""
        sed = '| sed \'s/arbitrary/RNA/g\''
        data = '%s %s' % (filename,sed)
        return data
    
    def _input_as_lines(self,data):
        """Overrides _input_as_lines in CommandLineApplication

        Input file need to be modified by sed""" 
        filename = super(fasta2col,self)._input_as_lines(data)

        sed = '| sed \'s/arbitrary/RNA/g\''
        data = '%s %s' % (filename,sed)
        return data


class findphyl(CommandLineApplication):
    """Application controller for findphyl in Pfold package 

    Find the phylogeny of the sequences using the
    neighbour joining approach"""

    _command = 'findphyl'
    _input_handler = '_input_as_string'

    def _input_as_string(self,filename):
        """Overrides  _input_as_string in CommandLineApplication

        scfg.rate file need to be specified along with the input file""" 
        file = '%s%s' % (pfold_path,'scfg.rate')
        data = '%s %s' % (file,filename)
        return data
    
    def _input_as_lines(self,data):
        """Overrides _input_as_lines in CommandLineApplication

        scfg.rate file need to be specified along with the input file"""
        filename = super(findphyl,self)._input_as_lines(data)

        file = '%s%s' % (pfold_path,'scfg.rate')
        data = '%s %s' % (file,filename)
        return data



class mltree(CommandLineApplication):
    """Application controller for mltree in pfold package

    Performs a maximum likelihood estimate of the branch lengths"""

    _command = 'mltree'
    _input_handler = '_input_as_string'

    def _input_as_string(self,filename):
        """Overrides  _input_as_string in CommandLineApplication

        scfg.rate file need to be specified along with the input file"""
        file = '%s%s' % (pfold_path,'scfg.rate')
        data = '%s %s' % (file,filename)
        return data
    
    def _input_as_lines(self,data):
        """Overrides _input_as_lines in CommandLineApplication

        scfg.rate file need to be specified along with the input file"""
        filename = super(mltree,self)._input_as_lines(data)

        file = '%s%s' % (pfold_path,'scfg.rate')
        data = '%s %s' % (file,filename)
        return data


class scfg(CommandLineApplication):
    """Application controller for scfg in Pfold package

    Performs the analysis

    The file `article.grm' has the grammar and evolutionary model
    that is used for the analysis"""

    _command = 'scfg' 
    _input_handler = '_input_as_string'

    def _input_as_string(self,filename):
        """Overrides  _input_as_string in CommandLineApplication

        Additional input information about tree needed from article.grm file"""
        file = '%s %s%s' % ('--treeinfile',pfold_path,'article.grm')
        data = '%s %s' % (file,filename)
        return data
    
    def _input_as_lines(self,data):
        """Overrides _input_as_lines in CommandLineApplication

        Additional input information about tree needed from article.grm file"""
        filename = super(scfg,self)._input_as_lines(data)

        file = '%s %s%s' % ('--treeinfile',pfold_path,'article.grm')
        data = '%s %s' % (file,filename)
        return data

