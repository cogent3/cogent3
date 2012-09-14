#!/usr/bin/env python 

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

class Contrafold(CommandLineApplication):
    """Application controler for CONTRAfold v1.0"""

    _parameters = {'predict':FlagParameter(Prefix='',Name='predict',Value=True),
                   'train':FlagParameter(Prefix='',Name='train')}

    _command = 'contrafold'
    _input_handler='_input_as_string'
