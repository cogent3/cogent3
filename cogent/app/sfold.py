#!/usr/bin/env python

"""Application controller for Sfold application

Option -o sets output directory,default ../output/
If working directory is changed during the work flow option -o
need to be reset. -o is currently set to cwd
"""

from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    MixedParameter,Parameters 
from os import getcwd
from sys import stderr

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

#IMPORTANT!!!!!
#must be specified if run outside sfold dir, default is ../param/ 
param_dir = 'your_path/sfold-2.0/param/'


class Sfold(CommandLineApplication):
    """Application controller Sfold application"""

    cwd = getcwd()
    _parameters = {
        '-o':ValuedParameter(Prefix='-',Name='o',Delimiter=' ',Value=cwd),
        '-a':ValuedParameter(Prefix='-',Name='a',Delimiter=' '),
        '-f':ValuedParameter(Prefix='-',Name='f',Delimiter=' '),
        '-l':ValuedParameter(Prefix='-',Name='l',Delimiter=' '),
        '-m':ValuedParameter(Prefix='-',Name='m',Delimiter=' '),
        '-w':ValuedParameter(Prefix='-',Name='w',Delimiter=' '),
        '-p':ValuedParameter(Prefix='-',Name='p',Delimiter=' ',Value=param_dir)
        }
    _command = 'sfold.X86_64.LINUX'
    _input_handler = '_input_as_string'


    def _get_result_paths(self,data):
        """Return a dict of ResultPath objects representing all possible output
        """
        result = {}

        #list of output files
        output_list = ['10structure','10structure_2','Dharmacon_thermo','bp',
                       'bprob','cdf','fe','loopr','oligo','oligo_f','pdf',
                       'sample','sample_1000','sirna','sirna_f','sirna_s',
                       'smfe','sstrand','stability']

        for name in output_list:
            try:
                path = '%s%s.out' % (self.WorkingDir,name)
                f = open(path)
                f.close()
                result[name] = \
                    ResultPath(Path=path)
            except IOError:
                stderr.write('path: %s\n' % self.WorkingDir)
                stderr.write('File: %s.out not written to result dict\n' % name)
                pass

        return result


