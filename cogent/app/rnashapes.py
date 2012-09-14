#!/usr/bin/env python
#file: RNAshaper.py

"""
Author: Shandy Wikman (ens01svn@cs.umu.se)

Application controller for RNAshapes application

Revision History:
2006 Shandy Wikman created file
"""

from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    MixedParameter,Parameters, _find_synonym, is_not_None

__author__ = "Daniel McDonald and Greg Caporaso"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

class RNAshapes(CommandLineApplication):
    """Application controller for RNAshapes application

    Options:
    -h          Display this information
    -H <option> Display detailed information on <option>
    -v          Show version

    Sequence analysis modes:
    -a          Shape folding (standard mode)
    -s          Complete suboptimal folding
    -p          Shape probabilities
    -q          Shape probabilities (including shreps)
    -P <value>  Shape probabilities for mfe-best shapes
    -i <value>  Sampling with <value> iterations
    -C          Consensus shapes (RNAcast)

    Additional modes (use with any of the above):
    -r          Calculate structure probabilities
    -w <value>  Specify window size
    -W <value>  Specify window position increment (use with -w) [1]
    -m <shape>  Match shape (use with -a, -s, -p, -q or -C)

    Analysis control:
    -e <value>  Set energy range (kcal/mol)
    -c <value>  Set energy range (%) [10]
    -t <value>  Specify shape type (1-5) [5]
    -F <value>  Set probability cutoff filter (use with -p, -q or -P)
    -T <value>  Set probability output filter (use with -p, -q or -P)
    -M <value>  Set maximum loop length [30]  (use -M n for unrestricted)
    -l          Allow lonely base pairs
    -u          Ignore unstable structures (use with -a, -s or -C)

    Input/Output:
    -o <value>  Specify output type (1-4,f) [2]
    -O <string> Specify output format string
    -S <value>  Specify output width for structures
    -# <value>  Print only the first <value> results
    -g <value>  Generate structure graphs for first <value> structures
    -L          Highlight upper case characters in structure graphs
    -N          Do not include additional information in graph output file
    -f <file>   Read input from <file>
    -B          Show progress bar (use with -p, -q or -P)
    -z          Enable colors (in interactive mode: disable colors)
    -Z          Enable colors for dotbracket and shape strings
    -D <string> Convert dotbracket-string to shape (choose type with -t)

    """

    _sequence_analysis = {
         '-a':FlagParameter(Prefix='-',Name='a',Value=True),
         '-s':FlagParameter(Prefix='-',Name='s',Value=False),
         '-p':FlagParameter(Prefix='-',Name='p',Value=False),
         '-q':FlagParameter(Prefix='-',Name='q',Value=False),
         '-P':ValuedParameter(Prefix='-',Name='P',Value=None,Delimiter=' '),
         '-i':ValuedParameter(Prefix='-',Name='i',Value=None,Delimiter=' '),
         '-C':FlagParameter(Prefix='-',Name='C',Value=True)}

    _additional_modes ={
        '-r':FlagParameter(Prefix='-',Name='r',Value=False)}

    _analysis_control = {
        '-e':ValuedParameter(Prefix='-',Name='e',Value=None,Delimiter=' '),
        '-c':ValuedParameter(Prefix='-',Name='c',Value=20,Delimiter=' '),
        '-t':ValuedParameter(Prefix='-',Name='t',Value=None,Delimiter=' '),
        '-F':ValuedParameter(Prefix='-',Name='F',Value=None,Delimiter=' '),
        '-T':ValuedParameter(Prefix='-',Name='T',Value=None,Delimiter=' '),
        '-M':ValuedParameter(Prefix='-',Name='M',Value=None,Delimiter=' '),
        '-l':ValuedParameter(Prefix='-',Name='l',Value=None,Delimiter=' '),
        '-u':ValuedParameter(Prefix='-',Name='u',Value=None,Delimiter=' ')}


    _input_output = {
         '-o':ValuedParameter(Prefix='-',Name='o',Value=None,Delimiter=' '),
         '-S':ValuedParameter(Prefix='-',Name='S',Value=None,Delimiter=' '),
         '-#':ValuedParameter(Prefix='-',Name='#',Value=None,Delimiter=' '),
         '-g':ValuedParameter(Prefix='-',Name='g',Value=None,Delimiter=' '),
         '-L':FlagParameter(Prefix='-',Name='L',Value=False),
         '-N':FlagParameter(Prefix='-',Name='N',Value=False),
         '-f':ValuedParameter(Prefix='-',Name='f',Delimiter=' '),
         '-B':FlagParameter(Prefix='-',Name='B',Value=False),
         '-z':FlagParameter(Prefix='-',Name='z',Value=False),
         '-Z':FlagParameter(Prefix='-',Name='Z',Value=False)}
    
    _parameters = {}
    _parameters.update(_sequence_analysis)
    _parameters.update(_additional_modes)
    _parameters.update(_analysis_control)
    _parameters.update(_input_output)
    
    _command = 'RNAshapes'
    _input_handler = '_input_as_string'

 
    def _input_as_lines(self,data):
        """Makes data the value of a specific parameter"""
        if data:
            self.Parameters['-f']\
                .on(super(RNAshapes,self)._input_as_lines(data))
        return ''

    def _input_as_string(self,data):
        """Makes data the value of a specific parameter
        
        This method returns the empty string. The parameter will be printed
        automatically once set.
        """
        if data:
            self.Parameters['-f'].on(data)
        return ''
        
