#!/usr/bin/env python

from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    MixedParameter,Parameters

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class PknotsRG(CommandLineApplication):
    """Application controller for PknotsRG v1.2 application

    Input: plain seqeunce
    
    pknotsRG is a tool for thermodynamic folding of RNA secondary
    structures, including the class of canonical simple recursive
    pseudoknots.

    Options:
    -m         Use mfe strategy
    -f         Use enf strategy
    -l         Use loc strategy
    -s         Show suboptimals
    -u         no dangling bases (implies -s)
    -o         no suboptimals inside pknots (implies -s -l)
    -e <value> Set energy range for suboptimals (kcal/mole)
    -c <value> Set energy range for suboptimals (%) [10]
    -n <value> Set npp-value [0.3]
    -p <value> Set pkinit-value [9]
    -k <value> Set maximal pknot-length


    """
    _parameters = {
        '-m':FlagParameter(Prefix='-',Name='m'),
        '-f':FlagParameter(Prefix='-',Name='f'),
        '-l':FlagParameter(Prefix='-',Name='l'),
        '-s':FlagParameter(Prefix='-',Name='s'),
        '-u':FlagParameter(Prefix='-',Name='u'),
        '-o':FlagParameter(Prefix='-',Name='o'),
        '-e':ValuedParameter(Prefix='-',Name='e',Delimiter=' '),
        '-c':ValuedParameter(Prefix='-',Name='c',Delimiter=' '),
        '-n':ValuedParameter(Prefix='-',Name='n',Delimiter=' '),
        '-p':ValuedParameter(Prefix='-',Name='p',Delimiter=' '),
        '-k':ValuedParameter(Prefix='-',Name='k',Delimiter=' ')}
        
    _command = 'pknotsRG-1.2-i386-linux-static'
    _input_handler = '_input_as_string'

    def _input_as_string(self,filename):
        """Returns '>filename' to redirect input to stdin"""
        return ''.join(['<',super(PknotsRG,self)._input_as_string(filename)])
    
    def _input_as_lines(self,data):
        """Returns '>temp_filename to redirect input to stdin"""
        return ''.join(['<',super(PknotsRG,self)._input_as_lines(data)])

