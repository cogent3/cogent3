#!/usr/bin/env python 

from cogent.app.util import CommandLineApplication,\
    CommandLineAppResult, ResultPath
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    MixedParameter,Parameters, _find_synonym

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"
                                      
class Knetfold(CommandLineApplication): 
    """Application controller for Knetfold v1.4.4b application"""

    _parameters = {'-i':ValuedParameter(Prefix='-',Name='i',Delimiter=' '),
                   '-n':ValuedParameter(Prefix='-',Name='n',Delimiter=' '),
                   '-d':ValuedParameter(Prefix='-',Name='d',Delimiter=' '),
                   '-m':ValuedParameter(Prefix='-',Name='m',Delimiter=' '),
                   '-q':ValuedParameter(Prefix='-',Name='q',Delimiter=' '),
                   '-o':ValuedParameter(Prefix='-',Name='o',Delimiter=' '),
                   '-r':ValuedParameter(Prefix='-',Name='r',Delimiter=' '),
                   '-f':ValuedParameter(Prefix='-',Name='f',Delimiter=' '),
                   '-h':ValuedParameter(Prefix='-',Name='h',Delimiter=' '),}
    

    _command = 'knetfold.pl'
    _input_handler = '_input_as_string'

    def _input_as_lines(self,data):
        """
        Infile needs to be set with parameter -i
        """
        filename = self._input_filename = self.getTmpFilename(self.WorkingDir)
        data_file = open(filename,'w')
        data_to_file = '\n'.join([str(d).strip('\n') for d in data])
        data_file.write(data_to_file)
        data_file.write('\n')
        data_file.close()

        #set input flag and give it the input filename
        self.Parameters['-i'].on(filename)
        return ''

    def _input_as_string(self,data):
        """Makes data the value of a specific parameter
        
        This method returns the empty string. The parameter will be printed
        automatically once set.
        """
        if data:
            self.Parameters['-i'].on(data)
        return ''

    def _get_result_paths(self,data):
        """
        Adds output files to the resultpath
        """
        result = {}

        if isinstance(data,list):
            filename=self._input_filename.split('/')[-1]
        else:
            filename=(data.split('/')[-1]).split('.')[0]
        #output files created in extensions list
        extensions = ['ct','coll','sec','fasta','pdf']
        #file = '%s%s%s' % (self.WorkingDir,filename,'_knet')
        file = ''.join([self.WorkingDir, filename, '_knet'])
        for ext in extensions:
            try:
                path = '%s.%s' % (file,ext)
                f = open(path)
                f.close()
                result[ext]=\
                    ResultPath(Path=(path))
            except IOError:
                pass
         
        number = 0
        #Unknown number of mx files, try/except to find all
        #file = '%s%s%s%d' % (self.WorkingDir,filename,'_knet.mx',number)
        file_base = ''.join([self.WorkingDir, filename, '_knet.mx'])
        while(True):
            try:
                file = file_base + str(number)
                f = open(file)
                result['%s%d' % ('mx',number)]= ResultPath(Path=(file))
                f.close()
                number +=1
            except IOError: #No more mx files
                break 
        return result
