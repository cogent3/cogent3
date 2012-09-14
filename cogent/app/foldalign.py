#!/usr/bin/env python

"""Application controller for Foldalign v2.0.3 application

Foldalign takes two sequences as input, these sequences can be in the same
file(2) or separate files.
ex1 Foldalign file(2) 
ex2 Foldalign seq1 seq2
"""

from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    MixedParameter,Parameters, _find_synonym, is_not_None

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class foldalign(CommandLineApplication):
    """Applictation controller for foldalign RNA secondary structure prediction 
    application
    """

    _parameters = {  
        '-max_length':ValuedParameter(Prefix='-',Name='max_length',Delimiter=' '),
        '-max_diff':ValuedParameter(Prefix='-',Name='max_diff',Delimiter=' '),
        '-score_matrix':ValuedParameter(Prefix='-',Name='score_matrix',Delimiter=' '),
        '-format':ValuedParameter(Prefix='-',Name='format',Delimiter=' '),
        '-plot_score':FlagParameter(Prefix='-',Name='plot_score'),
        '-global':FlagParameter(Prefix='-',Name='global'),
        '-summary':FlagParameter(Prefix='-',Name='summary'),}

    _command = 'foldalign'
    _input_handler = '_input_as_string'

    
