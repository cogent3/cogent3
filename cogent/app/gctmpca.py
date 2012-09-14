#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# gctmpca.py

"""Application controller for the Generalized Continuous-Time Markov 
 Process Coevolutionary Algorithm (GCTMPCA). GCTMPCA is presented in:
 
 Detecting coevolution in and among protein domains.
 Yeang CH, Haussler D., PLoS Comput Biol. 2007 Nov;3(11):e211.
 
 Detecting the coevolution of biosequences--an example 
 of RNA interaction prediction. Yeang CH, Darot JF, Noller HF, 
 Haussler D. Mol Biol Evol. 2007 Sep;24(9):2119-31. 
 
This code requires the GCTMPCA package to be installed. As of Nov. 2008,
 that software is available at:
 http://www.sns.ias.edu/~chyeang/coevolution_download.zip

Note that the authors did not name their algorithm or software when they 
 published it. GCTMPCA was suggested as a name by the first author via e-mail.
"""

from __future__ import division
from cogent.app.util import CommandLineApplication, ResultPath,\
    ApplicationError
from cogent.app.parameters import FilePath
from cogent.evolve.models import DSO78_freqs, DSO78_matrix

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Beta"

# Are these values in PyCogent somewhere?
gctmpca_base_order = 'ACGU'
default_gctmpca_rna_priors = {'A':0.2528,'C':0.2372,'G':0.3099,'U':0.2001}
default_gctmpca_rna_sub_matrix = """-1.4150\t0.2372\t0.9777\t0.2001
0.2528\t-1.1940\t0.3099\t0.6313
0.7976\t0.2372\t-1.2349\t0.2001
0.2528\t0.7484\t0.3099\t-1.3111"""

gctmpca_aa_order = 'ARNDCQEGHILKMFPSTWYV'
# By default, the Gctmpca method used the Dayhoff 78 frequencies and rate matrix
default_gctmpca_aa_priors = DSO78_freqs
default_gctmpca_aa_sub_matrix = """-133.941451\t1.104408\t3.962336\t5.624640\t1.205064\t3.404695\t9.806940\t21.266880\t0.773214\t2.397590\t3.499637\t2.092532\t1.062216\t0.715896\t12.670000\t28.456993\t21.719082\t0.000000\t0.717984\t13.461344
2.352429\t-86.970372\t1.293824\t0.000000\t0.769902\t9.410730\t0.049530\t0.797508\t8.068320\t2.360704\t1.280355\t37.343648\t1.327770\t0.556808\t5.220040\t10.714858\t1.522092\t2.109294\t0.239328\t1.553232
8.538446\t1.308928\t-179.776579\t42.419160\t0.000000\t3.940265\t7.330440\t12.317068\t17.985630\t2.840222\t2.902138\t25.593276\t0.014753\t0.556808\t2.128560\t34.440615\t13.406118\t0.241362\t2.842020\t0.970770
10.455240\t0.000000\t36.590960\t-142.144945\t0.000000\t5.126170\t57.108090\t11.076500\t2.891148\t0.885264\t0.000000\t5.714222\t0.000000\t0.000000\t0.658840\t6.609815\t3.863772\t0.000000\t0.000000\t1.164924
3.136572\t0.940792\t0.000000\t0.000000\t-26.760991\t0.000000\t0.000000\t0.974732\t0.941304\t1.622984\t0.000000\t0.000000\t0.000000\t0.000000\t0.962920\t11.201897\t0.936672\t0.000000\t2.871936\t3.171182
7.754303\t10.062384\t4.164496\t6.280848\t0.000000\t-124.487960\t35.463480\t2.481136\t20.372508\t0.663948\t6.231061\t12.313746\t1.681842\t0.000000\t7.754040\t3.896312\t3.102726\t0.000000\t0.000000\t2.265130
17.251146\t0.040904\t5.983936\t54.043416\t0.000000\t27.390580\t-136.769106\t7.177572\t1.445574\t2.250046\t0.938927\t6.680006\t0.442590\t0.000000\t2.584680\t5.496583\t1.990428\t0.000000\t0.658152\t2.394566
20.910480\t0.368136\t5.620048\t5.859000\t0.368214\t1.071140\t4.011930\t-65.418192\t0.336180\t0.000000\t0.597499\t2.173014\t0.250801\t0.596580\t1.723120\t16.281018\t1.756260\t0.000000\t0.000000\t3.494772
2.003921\t9.816960\t21.631120\t4.030992\t0.937272\t23.182530\t2.129790\t0.886120\t-88.051504\t0.258202\t3.755708\t2.092532\t0.000000\t1.909056\t4.763920\t2.435195\t1.287924\t0.283338\t3.799332\t2.847592
5.663255\t2.617856\t3.113264\t1.124928\t1.472856\t0.688590\t3.021330\t0.000000\t0.235326\t-128.487912\t21.936749\t3.702172\t4.957008\t7.795312\t0.608160\t1.669848\t11.240064\t0.000000\t1.106892\t57.534302
3.572207\t0.613560\t1.374688\t0.000000\t0.000000\t2.792615\t0.544830\t0.620284\t1.479192\t9.479702\t-53.327266\t1.448676\t7.774831\t6.244204\t1.621760\t1.182809\t1.931886\t0.482724\t0.837648\t11.325650
2.265302\t18.979456\t12.857376\t3.327912\t0.000000\t5.853015\t4.110990\t2.392524\t0.874068\t1.696756\t1.536426\t-74.828436\t3.584979\t0.000000\t1.672440\t6.679392\t7.961712\t0.000000\t0.388908\t0.647180
6.273144\t3.681360\t0.040432\t0.000000\t0.000000\t4.361070\t1.485900\t1.506404\t0.000000\t12.393696\t44.983139\t19.557126\t-125.902241\t3.659024\t0.861560\t4.313774\t6.088368\t0.000000\t0.000000\t16.697244
1.568286\t0.572656\t0.566048\t0.000000\t0.000000\t0.000000\t0.000000\t1.329180\t1.613664\t7.229656\t13.401049\t0.000000\t1.357276\t-54.612411\t0.557480\t3.200542\t0.761046\t0.797544\t20.881368\t0.776616
21.781750\t4.213112\t1.698144\t0.609336\t0.636006\t5.853015\t2.526030\t3.012808\t3.160092\t0.442632\t2.731424\t2.655906\t0.250801\t0.437492\t-74.727653\t17.046365\t4.566276\t0.000000\t0.000000\t3.106464
35.634943\t6.299216\t20.013840\t4.452840\t5.389314\t2.142280\t3.912870\t20.735208\t1.176630\t0.885264\t1.451069\t7.726272\t0.914686\t1.829512\t12.416600\t-160.924378\t32.198100\t0.787050\t1.017144\t1.941540
32.324117\t1.063504\t9.258928\t3.093552\t0.535584\t2.027515\t1.684020\t2.658360\t0.739596\t7.082112\t2.816781\t10.945552\t1.534312\t0.517036\t3.953040\t38.267350\t-129.918557\t0.000000\t1.256472\t10.160726
0.000000\t8.221704\t0.929936\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\t0.907686\t0.000000\t3.926422\t0.000000\t0.000000\t3.022672\t0.000000\t5.218275\t0.000000\t-24.051571\t1.824876\t0.000000
2.091048\t0.327232\t3.841040\t0.000000\t3.213504\t0.000000\t1.089660\t0.000000\t4.269486\t1.364782\t2.389996\t1.046266\t0.000000\t27.760856\t0.000000\t2.365618\t2.458764\t0.640134\t-54.670490\t1.812104
18.122416\t0.981696\t0.606480\t0.843696\t1.640226\t1.338925\t1.832610\t4.785048\t1.479192\t32.791654\t14.937475\t0.804820\t3.806274\t0.477264\t2.432640\t2.087310\t9.191094\t0.000000\t0.837648\t-98.996468"""

class Gctmpca(CommandLineApplication):
    """ App controller for the GCTMPCA algorithm for detecting sequence coevolution
    
        The Generalized Continuous-Time Markov Process Coevolutionary
         Algorithm (GCTMPCA) is presented in:
         
         Detecting coevolution in and among protein domains.
         Yeang CH, Haussler D., PLoS Comput Biol. 2007 Nov;3(11):e211.
         
         Detecting the coevolution of biosequences--an example 
         of RNA interaction prediction. Yeang CH, Darot JF, Noller HF, 
         Haussler D. Mol Biol Evol. 2007 Sep;24(9):2119-31. 
         
        This code requires the GCTMPCA package to be installed. As of 11/08,
         that software is available at:
          http://www.sns.ias.edu/~chyeang/coevolution_download.zip
    
    """

    _command = 'calculate_likelihood'
    _input_handler = '_gctmpca_cl_input'
    _data = {'mol_type':None,'comparison_type':0,'seqs1':None,\
             'seqs2':'-','tree1':None,'tree2':'-',\
             'seq_names':None,'species_tree':None,\
             'seq_to_species1':None,'seq_to_species2':'-',\
             'char_priors':None,'sub_matrix':None,'epsilon':0.7,\
             'max_gap_threshold':1.0,'max_seq_distance':1.0,\
             'covariation_threshold':0.0,'likelihood_threshold':0.0,\
             'output_path':None,'single_pair_only':0,'family_reps':'-',\
             'pos1':'','pos2':''}
    _parameter_order = ['mol_type','comparison_type','seqs1','seqs2',\
         'tree1','tree2','seq_names','species_tree',\
         'seq_to_species1','seq_to_species2','char_priors',\
         'sub_matrix','epsilon','max_gap_threshold','max_seq_distance',\
         'covariation_threshold','likelihood_threshold','output_path',\
         'single_pair_only','family_reps','pos1','pos2']

    _potential_paths = ['seqs1','tree1','seq_names',\
        'species_tree','seq_to_species1']

    _mol_type_lookup = {'rna':0,'0':0,'protein':1,'1':1}
    _default_priors = {0:default_gctmpca_rna_priors, 1:default_gctmpca_aa_priors}
    _default_sub_matrix = {0:default_gctmpca_rna_sub_matrix, 1:default_gctmpca_aa_sub_matrix}
    _char_order = {0:gctmpca_base_order,1:gctmpca_aa_order}
    _required_parameters = {}.fromkeys(['mol_type','seqs1','tree1',\
     'seq_names','species_tree','seq_to_species1'])

    def _set_command_line_parameters(self,data):
        """ Get the right setting for each command line parameter """
        # This function could be cleaned up.

        # for each command line parameter, set it to the value passed in or
        # the default value.
        for p in self._parameter_order:
            if p not in data:
                if p in self._required_parameters: 
                    raise ApplicationError,\
                     "Required parameter %s missing." % p
                else: data[p] = self._data[p]
            # Write necessary files to disk -- need to modify this so paths
            # to existing files can be passed in.
            if p in self._potential_paths:
                try:
                    data[p] = self._input_as_lines(data[p])
                except TypeError:
                    pass
        if data['single_pair_only'] == 1 and \
           not (data['pos1'] and data['pos2']):
            raise ApplicationError,\
                "Must specify pos1 and pos2 if single_pair_only == 1."

        # Make sure the MolType is in the correct format (i.e., 1 or 0)
        data['mol_type'] = mol_type = \
         self._mol_type_lookup[str(data['mol_type']).lower()]

        char_order = self._char_order[mol_type]
        # If we didn't get several values as parameters, set the defaults. 
        # These are done outside of the above loop b/c they require special 
        # handling.
        if not data['char_priors']: 
            data['char_priors'] = self._default_priors[mol_type]
        data['char_priors'] = \
             self._input_as_lines(\
              self._input_as_gctmpca_char_priors(\
              data['char_priors'],char_order))
        if not data['sub_matrix']: 
            data['sub_matrix'] = \
             self._input_as_multiline_string(\
              self._default_sub_matrix[mol_type])
        else:
            data['sub_matrix'] = \
             self._input_as_lines(\
              self._input_as_gctmpca_rate_matrix(\
              data['sub_matrix'],char_order))
        if not data['output_path']: 
            data['output_path'] = \
             self._input_as_path(self.getTmpFilename())
        return data

    def _gctmpca_cl_input(self,data):
        """ Write the list of 22 command line parameters to a string
        """
        # Get the right setting for each parameter
        data = self._set_command_line_parameters(data)
        # Explicitly disallow intermolecular experiments (I do this here to
        # make sure I'm looking at the final version of data)
        if data['comparison_type'] == 1: 
            raise NotImplementedError,\
                "Intermolecular experiments currently supported only via coevolve_alignments."
        # Create the command line parameter string and return it 
        return ' '.join([str(data[p]) for p in self._parameter_order]).strip()

    def _input_as_gctmpca_char_priors(self,priors,char_order):
        """convert dict of priors to string and write it to tmp file
        """
        # priors t be followed by a newline
        return ['\t'.join([str(priors[c]) for c in char_order]),'']

    def _input_as_gctmpca_rate_matrix(self,matrix,char_order):
        """convert 2D dict rate matrix to string and write it to tmp file
        """
        matrix_rows = []
        for c in char_order:
            matrix_rows.append('\t'.join([str(matrix[c][col_c]) \
                for col_c in char_order]))
        return matrix_rows

    def _get_result_paths(self,data):
        """A single file is written, w/ name specified in command line input
        """
        return {'output':ResultPath(Path=data['output_path'],IsWritten=True)}

