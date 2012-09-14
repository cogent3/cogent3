#!/usr/bin/env python

"""Application controller for mfold v3.2 application

!!!To get mfold app.controller to work TmpNameLen should be max 10 
(unknown reason).
"""

from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath
from cogent.app.parameters import Parameter, ValuedParameter,Parameters
from random import choice        

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"      
                        
class Mfold(CommandLineApplication): 
    """Application controller for mfold 3.2 application"""

    #Not all parameters included!
    #skipped: NA_CONC,MG_CONC,LB_FR,ROT_ANG,START,STOP,REUSE
    _parameters = {
        'LC':ValuedParameter(Prefix='',Name='LC=',Value=None,Delimiter=''),
        'T':ValuedParameter(Prefix='',Name='T=',Value=None,Delimiter=''),
        'P':ValuedParameter(Prefix='',Name='P=',Value=None,Delimiter=''),
        'MAXBP':ValuedParameter(Prefix='',Name='MAXBP=',Value=None,Delimiter=''),
        'MAX':ValuedParameter(Prefix='',Name='MAX=',Value=30,Delimiter=''),
        'MAX_LP':ValuedParameter(Prefix='',Name='MAX_LP=',Value=None,Delimiter=''),
        'MAX_AS':ValuedParameter(Prefix='',Name='MAX_AS=',Value=None,Delimiter=''),
        'MODE':ValuedParameter(Prefix='',Name='MODE=',Value=None,Delimiter=''),
        }

    _command = 'mfold'
    _input_handler = '_input_as_string'

    
    def _input_as_string(self,filename):
        """
        mfold dosen't take full paths so a tmp-file is created in the working 
        dir for mfold to read.
        """
        nr = choice(range(150))
        input_file = open(filename).readlines()
        filename = self._input_filename = 'mfold_in%d.txt' % nr
        data_file = open(filename,'w')
        data_to_file = '\n'.join([str(d).strip('\n') for d in input_file])
        data_file.write(data_to_file)
        data_file.close()
        data = '='.join(['SEQ',filename])
        return data
    
    def _input_as_lines(self,data):
        """
        Uses a fixed tmp filename since weird truncation of the generated 
        filename sometimes occured.
        """
        nr = choice(range(150))
        filename = self._input_filename = 'mfold_in%d.txt' % nr
        data_file = open(filename,'w')
        data_to_file = '\n'.join([str(d).strip('\n') for d in data])
        data_file.write(data_to_file)
        data_file.close()
        return '='.join(['SEQ',filename])

    def _get_result_paths(self,data):
        """Return a dict of ResultPath objects representing all possible output
        """
        result = {}
        itr=self.Parameters['MAX'].Value
        if itr == None:
            itr = 30

        filename=self._input_filename.split('/')[-1]
        for i in range(1,itr+1):
            try:
                ct = self.WorkingDir+filename+'_'+str(i)+'.ct'
                f = open(ct)
                f.close()
                result['ct'+str(i)] =\
                    ResultPath(Path=ct)

                pdf = self.WorkingDir+filename+'_'+str(i)+'.pdf'
                f = open(pdf)
                f.close()
                result['pdf'+str(i)] =\
                    ResultPath(Path=pdf)
            except IOError:
                pass
        result['ct_all'] =\
            ResultPath(Path=(self.WorkingDir+filename+'.ct'))

        name = self.WorkingDir+filename
        #output files
        files = ['log',
                 'ann',
                 'h-num',
                 'det',
                 'pnt',
                 'sav',
                 'ss-count',
                 '-local.seq',
                 'rnaml',
                 'out',
                 'plot',
                 'ps',
                 '_1.ps',
                 '_1.ss']
        for f in files:
            if f == '-local.seq':
                file = ''.join([name, f])
            elif f.startswith('_1'):
                file = ''.join([name, f])
            else:
                file = '.'.join([name, f]) 
            result['%s' % f] = ResultPath(Path=file)
 
        return result
