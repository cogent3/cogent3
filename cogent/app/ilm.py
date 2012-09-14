#!/usr/bin/env python

from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath
from cogent.app.parameters import Parameter,ValuedParameter,Parameters

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class ILM(CommandLineApplication):
    """Application controller ILM application

    Predict a secondary structure given a score matrix

    Main options:
    -L l:      minimum loop length (default=3)
    -V v:      minimum virtual loop length (default=3)
    -H h:      minimum helix length (default=3)
    -N n:      number of helices selected per iteration (default=1)
    -I i:      number of iterations before termination(default=unlimited)
    """

    _parameters = {
        '-L':ValuedParameter(Prefix='-',Name='L',Delimiter=' '),
        '-V':ValuedParameter(Prefix='-',Name='V',Delimiter=' '),
        '-H':ValuedParameter(Prefix='-',Name='H',Delimiter=' '),
        '-N':ValuedParameter(Prefix='-',Name='N',Delimiter=' '),
        '-I':ValuedParameter(Prefix='-',Name='I',Delimiter=' ')}

    _command = 'ilm'
    _input_handler = '_input_as_string'

class hlxplot(CommandLineApplication):
    """Application controller hlxplot application

    Compute a helix plot score matrix from a sequence alignment

    Options:
    -b B: Set bad pair penalty to B
            (Default = 2)
    -g G: Set good pair score to G
            (Default = 1)
    -h H: Set minimum helix length to H
            (Default = 2)
    -l L: Set minimum loop length to L
            (Default = 3)
    -s S: Set helix length score to S
            (Default = 2.0)
    -t  : Write output in text format
            (Default = Binary format)
    -x X: Set paired gap penalty to X
            (Default = 3)
    """

    _parameters = {
        '-b':ValuedParameter(Prefix='-',Name='b',Delimiter=' '),
        '-g':ValuedParameter(Prefix='-',Name='g',Delimiter=' '),
        '-h':ValuedParameter(Prefix='-',Name='h',Delimiter=' '),
        '-l':ValuedParameter(Prefix='-',Name='l',Delimiter=' '),
        '-s':ValuedParameter(Prefix='-',Name='s',Delimiter=' '),
        '-t':ValuedParameter(Prefix='-',Name='t',Delimiter=' '),
        '-x':ValuedParameter(Prefix='-',Name='x',Delimiter=' ')}

    _command = 'hlxplot'
    _input_handler = '_input_as_string'

class xhlxplot(CommandLineApplication):
    """Application controller xhlxplot application

    Compute an extended helix plot score matrix from a single sequence

    Options:
    -b B: Set bad pair penalty to B
            (Default = 200)
    -h H: Set minimum helix length to H
            (Default = 2)
    -l L: Set minimum loop length to L
            (Default = 3)
    -x X: Set paired gap penalty to X
            (Default = 500)

    -t  : Write output in text format
            (Default = Binary format)
    -c  : No Closing GU
            (Default = allows closing GU)
    """

    _parameters = {
        '-b':ValuedParameter(Prefix='-',Name='b',Delimiter=' '),
        '-h':ValuedParameter(Prefix='-',Name='h',Delimiter=' '),
        '-l':ValuedParameter(Prefix='-',Name='l',Delimiter=' '),
        '-x':ValuedParameter(Prefix='-',Name='x',Delimiter=' '),
        '-t':ValuedParameter(Prefix='-',Name='t',Delimiter=' '),
        '-c':ValuedParameter(Prefix='-',Name='c',Delimiter=' ')}

    _command = 'xhlxplot'
    _input_handler = '_input_as_string'

    
    
