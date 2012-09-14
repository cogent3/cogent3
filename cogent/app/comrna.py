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

class comRNA(CommandLineApplication):
    """Application controller comRNA v1.80

    for comRNA options type comRNA at the promt
    """

    _parameters = {
        '-L':ValuedParameter(Prefix='',Name='L',Value=None,Delimiter=' '),
        '-E':ValuedParameter(Prefix='',Name='E',Value=None,Delimiter=' '),
        '-S':ValuedParameter(Prefix='',Name='S',Value=None,Delimiter=' '),
        '-Sh':ValuedParameter(Prefix='',Name='Sh',Value=None,Delimiter=' '),
        '-Sl':ValuedParameter(Prefix='',Name='Sl',Value=None,Delimiter=' '),
        '-P':ValuedParameter(Prefix='',Name='P',Value=None,Delimiter=' '),
        '-n':ValuedParameter(Prefix='',Name='n',Value=None,Delimiter=' '),
        '-x':ValuedParameter(Prefix='',Name='x',Value=None,Delimiter=' '),
        '-m':ValuedParameter(Prefix='',Name='m',Value=None,Delimiter=' '),
        '-tp':ValuedParameter(Prefix='',Name='tp',Value=None,Delimiter=' '),
        '-ts':ValuedParameter(Prefix='',Name='ts',Value=None,Delimiter=' '),
        '-a':ValuedParameter(Prefix='',Name='a',Value=None,Delimiter=' '),
        '-o':ValuedParameter(Prefix='',Name='o',Value=None,Delimiter=' '),
        '-c':ValuedParameter(Prefix='',Name='c',Value=None,Delimiter=' '),
        '-j':ValuedParameter(Prefix='',Name='j',Value=None,Delimiter=' '),
        '-r':ValuedParameter(Prefix='',Name='r',Value=None,Delimiter=' '),
        '-f':ValuedParameter(Prefix='',Name='f',Value=None,Delimiter=' '),
        '-v':ValuedParameter(Prefix='',Name='v',Value=None,Delimiter=' '),
        '-g':ValuedParameter(Prefix='',Name='g',Value=None,Delimiter=' '),
        '-d':ValuedParameter(Prefix='',Name='d',Value=None,Delimiter=' '),
        '-wa':ValuedParameter(Prefix='',Name='wa',Value=None,Delimiter=' '),
        '-wb':ValuedParameter(Prefix='',Name='wb',Value=None,Delimiter=' '),
        '-wc':ValuedParameter(Prefix='',Name='wc',Value=None,Delimiter=' '),
        '-wd':ValuedParameter(Prefix='',Name='wd',Value=None,Delimiter=' '),
        '-we':ValuedParameter(Prefix='',Name='we',Value=None,Delimiter=' '),
        '-pk':ValuedParameter(Prefix='',Name='pk',Value=None,Delimiter=' '),
        '-pg':ValuedParameter(Prefix='',Name='pg',Value=None,Delimiter=' '),
        '-pd':ValuedParameter(Prefix='',Name='pd',Value=None,Delimiter=' '),
        '-ps':ValuedParameter(Prefix='',Name='ps',Value=None,Delimiter=' '),
        '-pc':ValuedParameter(Prefix='',Name='pc',Value=None,Delimiter=' ')}

    _command = 'comRNA'
    _input_handler = '_input_as_string'


