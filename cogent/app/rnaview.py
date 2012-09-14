#!/usr/bin/env python
"""Application controller for the RnaView package
"""

from cogent.app.parameters import FlagParameter, ValuedParameter, MixedParameter
from cogent.app.util import CommandLineApplication, ResultPath

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Greg Caporaso", "Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Production"

class RnaView(CommandLineApplication):
    """ The Application controller for the RnaView application 

        There are two known issues with this applications controller
         which are being addressed: 
         (1) The first is that it doesn't handle
         the functionality where -a is passed followed by a filename and
         a float. This will be addressed shortly. 
         (2) The second problem is that under
         some mysterious circumstances, rnaview (the actual application,
         not via cogent) writes extra output files. This seems to be a bug
         in rnaview, and cogent therefore doesn't clean them up. One time 
         when this occurs is when there are spaces in filenames or paths.
         Contact Greg (gregcaporaso@gmail.com) for some example output
         illustrating this issue. Users are warned to look for extra 
         output files either showing up in the working directory or 
         the temp directory (default: /tmp).

    """

    # The functionality necessary for use of the -a parameter is still 
    # under development and is not ready for use
    _parameters = {\
        '-p':FlagParameter(Prefix='-',Name='p'),\
        '-v':FlagParameter(Prefix='-',Name='v'),\
        '-c':ValuedParameter(Prefix='-',Name='c',Delimiter=' '),\
        '-a':FlagParameter(Prefix='-',Name='a'),\
        '-x':FlagParameter(Prefix='-',Name='x')}
    _command = 'rnaview'
    _input_handler = '_input_as_path'

    ### Everything above is necessary for sub-class, code below is for cases 
    ### where files are written (ie. data goes to places other than stdout and
    ### stderr) Complexity increases with the amount of variability in the 
    ### file name. For rnaview the naming of the files is quite complex, hence
    ### the large amount of code necessary.

    def _get_result_paths(self,data):
        # There are two possibilities for what data will be, if -a has been 
        # specified, data will be a file containing a space-delimited list of 
        # pdb files. If -a has not been specified data will be the name of a
        # single pdb file to act on
        result = {}
        if self.Parameters['-a'].isOff():
            # If we have created a temp file containing data we need that
            # temp file name
            if self._input_filename:
                file_prefix = self._get_pdb_filename(self._input_filename)
                out_path = self._get_out_path(self._input_filename)
            # Otherwise we will just be passing a filename as data
            else:
                file_prefix = self._get_pdb_filename(data)
                out_path = self._get_out_path(data)
                
            result.update(\
             self._construct_result_file_set(file_prefix=file_prefix,\
                out_path=out_path))
        else:
            inputs = data.split(' ')
            f = open(inputs[0])
            pdb_files = f.read().split(' ')
            f.close()
            for p in pdb_files:
                file_prefix = self._get_pdb_filename(p)
                key_prefix = ''.join([file_prefix,'_'])
                result.update(\
                 self._construct_result_file_set(key_prefix=key_prefix,\
                 file_prefix=file_prefix,out_path=self._get_out_path(p)))
            
        return result
    
    def _construct_result_file_set(self, key_prefix='', file_prefix='',\
        out_path=''):
        result = {}
        result['bp_stats'] = ResultPath(Path=\
            ''.join([self.WorkingDir,'/base_pair_statistics.out']))
        result[''.join([key_prefix,'base_pairs'])] =\
            ResultPath(Path=''.join([out_path,file_prefix,'.out']))
        result[''.join([key_prefix,'ps'])] =\
            ResultPath(Path=''.join([out_path,file_prefix,'.ps']),\
            IsWritten=self.Parameters['-p'].isOn())
        result[''.join([key_prefix,'vrml'])] =\
            ResultPath(Path=''.join([out_path,file_prefix,'.wrl']),\
            IsWritten=self.Parameters['-v'].isOn())
        result[''.join([key_prefix,'xml'])] =\
            ResultPath(Path=''.join([out_path,file_prefix,'.xml']),\
            IsWritten=self.Parameters['-x'].isOn())
        result['best_pair'] = ResultPath(Path=\
            ''.join([self.WorkingDir,'/best_pair.out']),\
            IsWritten=self.Parameters['-a'].isOn())
        result['pattern_tmp'] = ResultPath(Path=\
            ''.join([self.WorkingDir,'/pattern_tmp.out']),\
            IsWritten=self.Parameters['-a'].isOn())
        result[''.join([key_prefix,'patt'])] =\
            ResultPath(Path=''.join([out_path,file_prefix,'_patt.out']),\
            IsWritten=self.Parameters['-a'].isOn())
        result[''.join([key_prefix,'patt_tmp'])] =\
            ResultPath(Path=''.join([out_path,file_prefix,'_patt_tmp.out']),\
            IsWritten=self.Parameters['-a'].isOn())
        result[''.join([key_prefix,'sort.out'])] =\
            ResultPath(Path=''.join([out_path,file_prefix,'_sort.out']),\
            IsWritten=self.Parameters['-a'].isOn())
    

        return result

    def _accept_exit_status(self,exit_status):
        "Return False if exit_status is not zero """
        if exit_status != 0:
            return False
        return True

    def _get_pdb_filename(self,word):
        """Returns the file prefix of the _input_filename.
        
        If the file is an NMR file, Rnaview creates a new PDB file containing
        only the used model. This file is the original input filename plus
        '_nmr.pdb'. The resulting .out file containing the base pairs uses
        that prefix, e.g. xxx.ent_nmr.pdb.out

        Since NMR and X-RAY files cannot be safely distinguished based on the
        filename, we have to open the file and check the EXPDTA field for 
        'NMR' or detect the word MODEL in there.
        """
        nmr = False
        f = open(word)
        for line in f:
            if line.startswith('EXPDTA') and 'NMR' in line:
                nmr=True
                break
            if line.startswith('MODEL'):
                nmr=True
                break
        f.close()
        start_index = word.rfind('/')
        if nmr:
            return word[start_index+1:].strip() + '_nmr.pdb'
        return word[start_index+1:].strip()

    def _get_out_path(self,word):
        end_index = word.rfind('/')
        if end_index >= 0:
            return word[0:end_index+1].strip()
        return ''
