#!/usr/bin/env python

"""Application controller for unafold 3.2, hybrid_ss_min application, package

Rest fo unafold package might be added at a later time
"""

from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    Parameters

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"
                                      
class hybrid_ss_min(CommandLineApplication):
    """Application controller for hybrid_ss_min application

    Computes a minimum energy folding for a RNA or DNA sequence
    """
    
    #not all options supported here!
    #left out: -n,-d,-r,-f,-I,-q,-c,-b and all obscure options
    _parameters = {
	't':ValuedParameter(Prefix='-',Name='t',Value=None,Delimiter=' '),
	'T':ValuedParameter(Prefix='-',Name='T',Value=None,Delimiter=' '),
	'i':ValuedParameter(Prefix='-',Name='i',Value=None,Delimiter=' '),
	's':ValuedParameter(Prefix='-',Name='s',Value=None,Delimiter=' '),
	'o':ValuedParameter(Prefix='-',Name='o',Value=None,Delimiter=' '),
	'N':ValuedParameter(Prefix='-',Name='N',Value=None,Delimiter=' '),
	'M':ValuedParameter(Prefix='-',Name='M',Value=None,Delimiter=' '),
	'p':ValuedParameter(Prefix='-',Name='p',Value=None,Delimiter=' '),
	'E':FlagParameter(Prefix='-',Name='E'),
	'F':ValuedParameter(Prefix='-',Name='F',Value=None,Delimiter=' '),
	'm':ValuedParameter(Prefix='-',Name='m',Value=None,Delimiter=' ')}
	


    _command = 'hybrid-ss-min'
    _input_handler = '_input_as_string'

    def _get_result_paths(self,data):
	"""Return a dict of ResultPath objects representing all possible output

    This dictionary will have keys based
    on the name that you'd like to access the file by in the 
    CommandLineAppResult object that will be created, and the values
    which are ResultPath objects."""

	result = {}
        #UNAfold default values
	start_tmp = 37
	step = 1
	end_tmp = 38

	if isinstance(data,list):
	    filename=self._input_filename.split('/')[-1]
	else:
	    filename=data.split('/')[-1]
    
	result['ct']= \
	    ResultPath(Path=(self.WorkingDir+filename+'.ct'))
	result['dG'] = \
	    ResultPath(Path=(self.WorkingDir+filename+'.dG'))
	result['run'] = \
	    ResultPath(Path=(self.WorkingDir+filename+'.run'))

	#if temp interval is not default it will result in more output-files
        #one for every temp
	if self.Parameters['t'].Value is not None:
	    start_tmp = self.Parameters['t'].Value
	if self.Parameters['T'].Value is not None:
	    end_tmp = self.Parameters['T'].Value + 1
	if self.Parameters['i'].Value is not None:
	    step = self.Parameters['i'].Value
	for i in range(start_tmp,end_tmp,step):
	    temp = str(i)
	    result['plot_%d' % i]= \
	    ResultPath(Path=(self.WorkingDir+filename+'.'+temp+'.plot'))
	    result['ext_%d' % i]= \
	    ResultPath(Path=(self.WorkingDir+filename+'.'+temp+'.ext'))

	return result
