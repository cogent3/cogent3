#!/usr/bin/env python
"""Provides an application controller for the commandline version of:
DOTUR v1.53
"""
import shutil
from cogent.app.parameters import FlagParameter, ValuedParameter, \
    MixedParameter
from cogent.app.util import CommandLineApplication, ResultPath, \
    get_tmp_filename, FilePath
from cogent.core.alignment import SequenceCollection, Alignment
from cogent.core.moltype import DNA, RNA, PROTEIN
from cogent.format.table import phylipMatrix
from cogent.parse.dotur import OtuListParser

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

class Dotur(CommandLineApplication):
    """Dotur application controller.
    """
    # Options:
    _options = {\
        # -i:		Number of iterations (default = 1000)
        '-i':ValuedParameter('-',Name='i',Delimiter=' '),\
        # -c:		Clustering method - (f) furthest neighbor, (n) nearest
        #           neighbor, (a) average neighbor (default = f)
        '-c':ValuedParameter('-',Name='c',Delimiter=' '),\
        # -p:		Precision of distances for output, increasing can
        #           dramatically lengthen execution times - 10, 100, 1000, 10000 
        #           (default = 100)
        '-p':ValuedParameter('-',Name='p',Delimiter=' '),\
        # -l:		Input file is lower triangular (default = square matrix)
        '-l':FlagParameter('-',Name='l'),\
        # -r:		Calculates rarefaction curves for each parameter, can
        #           dramatically lengthen execution times.  Simple rarefaction
        #           curve always calculated.
        '-r':FlagParameter('-',Name='r'),\
        # -stop:	Stops clustering when cutoff has been reached.
        '-stop':FlagParameter('-',Name='stop'),\
        # -wrep:	Samples with replacement.
        '-wrep':FlagParameter('-',Name='wrep'),\
        # -jumble:	Jumble the order of the distance matrix.
        '-jumble':FlagParameter('-',Name='jumble'),\
        # -sim:		Converts similarity score to distance (D=1-S).
        '-sim':FlagParameter('-',Name='sim'),\
         }
        
    _parameters = {}
    _parameters.update(_options)
    _input_handler = '_input_as_multiline_string'
    _command = 'dotur'

    def getHelp(self):
        """Method that points to the DOTUR documentation."""
        help_str =\
        """
        See DOTUR Documentation page at:
        http://schloss.micro.umass.edu/software/dotur/documentation.html
        """
        return help_str
    
    def _input_as_multiline_string(self, data):
        """Write a multiline string to a temp file and return the filename.

            data: a multiline string to be written to a file.

           * Note: the result will be the filename as a FilePath object 
            (which is a string subclass).

        """
        filename = self._input_filename = \
            FilePath(self.getTmpFilename(self.WorkingDir))
        data_file = open(filename,'w')
        data_file.write(data)
        data_file.close()
        return filename
    
    def _get_cluster_method(self):
        """Returns cluster method as string.
        """
        if self.Parameters['-c'].isOn():
            cluster_method = self._absolute(str(\
                self.Parameters['-c'].Value))+'n'
        else:
            # f (furthest neighbor) is default
            cluster_method = 'fn'
        
        return cluster_method
    
    def _get_result_paths(self,data):
        """Return dict of {key: ResultPath}
        
            - NOTE: Only putting a few files on the results path.  Add more
                here if needed.
        """
        result = {}
        out_name = self._input_filename.split('.txt')[0]
        cluster_method = self._get_cluster_method()
        #only care about Otu, List and Rank, can add others later.
        result['Otu'] = ResultPath(Path=out_name+'.%s.otu'%(cluster_method))
        result['List'] = ResultPath(Path=out_name+'.%s.list'%(cluster_method))
        result['Rank'] = ResultPath(Path=out_name+'.%s.rank'%(cluster_method))
        result['Rarefaction'] = \
            ResultPath(Path=out_name+'.%s.rarefaction'%(cluster_method))
        return result

def remap_seq_names(otu_list, int_map):
    """Returns list with seq names remapped.
        - otu_list: list of lists containing sequence names in an OTU.
        - int_map: mapping between names in otu_list and original names.
    """
    res = []
    for otu in otu_list:
        curr_otu = []
        for seq in otu:
            curr_otu.append(int_map[seq])
        res.append(curr_otu)
    return res

def dotur_from_alignment(aln,moltype,distance_function,params=None):
    """Returns dotur results given an alignment and distance function.
    
        - aln: An Alignment object or something that behaves like one.
            Sequences must be aligned.
        - moltype: cogent.core.moltype object.
        - distance_function: function that can be passed to distanceMatrix()
            method of SequenceCollection.  Must be able to find distance
            between two sequences.
        
        - NOTE:  This function will only return the parsed *.list file, as
            it contains the OTU identities.
            Dotur generates 23 output files, so if this is not the one you
            are looking for, check out the documentation and add the others
            to the result path.
    """
    #construct Alignment object.  This will handle unaligned sequences.
    aln = Alignment(aln, MolType=moltype)
    
    #need to make int map.
    int_map, int_keys = aln.getIntMap()
    #construct Alignment object from int map to use object functionality
    int_map = Alignment(int_map, MolType=moltype)
    order = sorted(int_map.Names)
    
    #Build distance matrix.
    d_matrix_dict = int_map.distanceMatrix(f=distance_function)
    d_matrix_dict.RowOrder=order
    d_matrix_dict.ColOrder=order
    
    #Get distance matrix in list form.
    d_matrix_list = d_matrix_dict.toLists()
    
    #must be strings to use phylipMatrix
    for i,line in enumerate(d_matrix_list):
        d_matrix_list[i]=map(str,line)
    
    #Get phylip formatted string.
    phylip_matrix_string = phylipMatrix(rows=d_matrix_list,names=order)
        
    working_dir = get_tmp_filename(suffix='')
    app = Dotur(InputHandler='_input_as_multiline_string',\
        WorkingDir=working_dir,params=params)
    
    res = app(phylip_matrix_string)
    
    otu_list = OtuListParser(res['List'].readlines())
    
    #remap sequence names
    for i,otu in enumerate(otu_list):
        otu_list[i][2]=remap_seq_names(otu[2], int_keys)
    
    shutil.rmtree(app.WorkingDir)
    
    return otu_list
    

def dotur_from_file(distance_matrix_file_path,params=None):
    """Returns dotur results given a distance matrix file.
    
        - distance_matrix_file_path:  Path to distance matrix file.  This file
             must a PHYLIP formatted square distance matrix.  This format
             is available in cogent.format.table.
             - IMPORANT NOTE:  This distance matrix format allows only 10
                characters for the row labels in the distance matrix.  Also,
                the IDs must be unique and ungapped to be useful when using
                dotur.
        - NOTE:  This function will only return the parsed *.list file, as
            it contains the OTU identities.
            Dotur generates 23 output files, so if this is not the one you
            are looking for, check out the documentation and add the others
            to the result path.
    """
    # Read out the data from the distance_matrix_file_path.
    # This is important so we can run dotur in a temp directory and avoid
    # having to handle all 23 output files.
    d_matrix_string = open(distance_matrix_file_path,'U').read()
    
    working_dir = get_tmp_filename(suffix='')
    app = Dotur(InputHandler='_input_as_multiline_string',\
        WorkingDir=working_dir,params=params)
    
    res = app(d_matrix_string)
    
    otu_list = OtuListParser(res['List'].readlines())
    
    shutil.rmtree(app.WorkingDir)
    
    return otu_list

