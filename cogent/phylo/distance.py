#!/usr/bin/env python
"""Estimating pairwise distances between sequences.
"""
from cogent.util import parallel, table, warning, progress_display as UI
from cogent.maths.stats.util import Numbers
from cogent import LoadSeqs, LoadTree

from warnings import warn

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell", "Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class EstimateDistances(object):
    """Base class used for estimating pairwise distances between sequences.
    Can also estimate other parameters from pairs."""
    
    def __init__(self, seqs, submodel, threeway=False, motif_probs = None,
                do_pair_align=False, rigorous_align=False, est_params=None,
                modify_lf=None):
        """Arguments:
            - seqs: an Alignment or SeqCollection instance with > 1 sequence
            - submodel: substitution model object Predefined models can
              be imported from cogent.evolve.models
            - threeway: a boolean flag for using threeway comparisons to
              estimate distances. default False. Ignored if do_pair_align is
              True.
            - do_pair_align: if the input sequences are to be pairwise aligned
              first and then the distance will be estimated. A pair HMM based
              on the submodel will be used.
            - rigorous_align: if True the pairwise alignments are actually
              numerically optimised, otherwise the current substitution model
              settings are used. This slows down estimation considerably.
            - est_params: substitution model parameters to save estimates from
              in addition to length (distance)
            - modify_lf: a callback function for that takes a likelihood
              function (with alignment set) and modifies it. Can be used to
              configure local_params, set bounds, optimise using a restriction
              for faster performance.
        
        Note: Unless you know a priori your alignment will be flush ended
        (meaning no sequence has terminal gaps) it is advisable to construct a
        substitution model that recodes gaps. Otherwise the terminal gaps will
        significantly bias the estimation of branch lengths when using
        do_pair_align.
        """
        
        if do_pair_align:
            self.__threeway = False
        else:
            # whether pairwise is to be estimated from 3-way
            self.__threeway = [threeway, False][do_pair_align]
        
        self.__seq_collection = seqs
        self.__seqnames = seqs.getSeqNames()
        self.__motif_probs = motif_probs
        # the following may be pairs or three way combinations
        self.__combination_aligns = None
        self._do_pair_align = do_pair_align
        self._rigorous_align = rigorous_align
        # substitution model stuff
        self.__sm = submodel
        
        self._modify_lf = modify_lf
        # store for the results
        self.__param_ests = {}
        self.__est_params = list(est_params or [])
        
        self.__run = False # a flag indicating whether estimation completed
        # whether we're on the master CPU or not
        self._on_master_cpu = parallel.getCommunicator().Get_rank() == 0
    
    def __str__(self):
        return str(self.getTable())
    
    def __make_pairwise_comparison_sets(self):
        comps = []
        names = self.__seq_collection.getSeqNames()
        n = len(names)
        for i in range(0, n - 1):
            for j in range(i + 1, n):
                comps.append((names[i], names[j]))
        return comps
    
    def __make_threeway_comparison_sets(self):
        comps = []
        names = self.__seq_collection.getSeqNames()
        n = len(names)
        for i in range(0, n - 2):
            for j in range(i + 1, n - 1):
                for k in range(j + 1, n):
                    comps.append((names[i], names[j], names[k]))
        return comps
    
    def __make_pair_alignment(self, seqs, opt_kwargs):
        lf = self.__sm.makeLikelihoodFunction(\
                    LoadTree(tip_names=seqs.getSeqNames()),
                    aligned=False)
        lf.setSequences(seqs.NamedSeqs)
        
        # allow user to modify the lf config
        if self._modify_lf:
            lf = self._modify_lf(lf)
        
        if self._rigorous_align:
            lf.optimise(**opt_kwargs)
        lnL = lf.getLogLikelihood()
        (vtLnL, aln) = lnL.edge.getViterbiScoreAndAlignment()
        return aln
    
    @UI.display_wrap
    def __doset(self, sequence_names, dist_opt_args, aln_opt_args, ui):
        # slice the alignment
        seqs = self.__seq_collection.takeSeqs(sequence_names)
        if self._do_pair_align:
            ui.display('Aligning', progress=0.0, current=.5)
            align = self.__make_pair_alignment(seqs, aln_opt_args)
            ui.display('', progress=.5, current=.5)
            
        else:
            align = seqs
            ui.display('', progress=0.0, current=1.0)
        # note that we may want to consider removing the redundant gaps
        
        # create the tree object
        tree = LoadTree(tip_names = sequence_names)
        
        # make the parameter controller
        lf = self.__sm.makeLikelihoodFunction(tree)
        if not self.__threeway:
            lf.setParamRule('length', is_independent = False)
        
        if self.__motif_probs:
            lf.setMotifProbs(self.__motif_probs)
        
        lf.setAlignment(align)
        
        # allow user modification of lf using the modify_lf
        if self._modify_lf:
            lf = self._modify_lf(lf)
        
        lf.optimise(**dist_opt_args)
                
        # get the statistics
        stats_dict = lf.getParamValueDict(['edge'], 
                params=['length'] + self.__est_params)
        
        # if two-way, grab first distance only
        if not self.__threeway:
            result = {'length': stats_dict['length'].values()[0] * 2.0}
        else:
            result = {'length': stats_dict['length']}
        
        # include any other params requested
        for param in self.__est_params:
            result[param] = stats_dict[param].values()[0]
            
        return result
    
    @UI.display_wrap
    def run(self, dist_opt_args=None, aln_opt_args=None, ui=None, **kwargs):
        """Start estimating the distances between sequences. Distance estimation
        is done using the Powell local optimiser. This can be changed using the
        dist_opt_args and aln_opt_args.
        
        Arguments:
            - show_progress: whether to display progress. More detailed progress
              information from individual optimisation is controlled by the
              ..opt_args.
            - dist_opt_args, aln_opt_args: arguments for the optimise method for
              the distance estimation and alignment estimation respectively."""
        
        if 'local' in kwargs:
              warn("local argument ignored, provide it to dist_opt_args or"\
              " aln_opt_args", DeprecationWarning, stacklevel=2)
        
        ui.display("Distances")
        dist_opt_args = dist_opt_args or {}
        aln_opt_args = aln_opt_args or {}
        # set the optimiser defaults
        dist_opt_args['local'] = dist_opt_args.get('local', True)
        aln_opt_args['local'] = aln_opt_args.get('local', True)
        # generate the list of unique sequence sets (pairs or triples) to be
        # analysed
        if self.__threeway:
            combination_aligns = self.__make_threeway_comparison_sets()
            desc = "triplet "
        else:
            combination_aligns = self.__make_pairwise_comparison_sets()
            desc = "pair "
        labels = [desc + ','.join(names) for names in combination_aligns]
                            
        def _one_alignment(comp):
            result = self.__doset(comp, dist_opt_args, aln_opt_args)
            return (comp, result)
        
        for (comp, value) in ui.imap(_one_alignment, combination_aligns,
                labels=labels):
            self.__param_ests[comp] = value
    
    def getPairwiseParam(self, param, summary_function="mean"):
        """Return the pairwise statistic estimates as a dictionary keyed by
        (seq1, seq2)
        
        Arguments:
            - param: name of a parameter in est_params or 'length'
            - summary_function: a string naming the function used for
              estimating param from threeway distances. Valid values are 'mean'
              (default) and 'median'."""
        summary_func = summary_function.capitalize()
        pairwise_stats = {}
        assert param in self.__est_params + ['length'], \
                "unrecognised param %s" % param
        if self.__threeway and param == 'length':
            pairwise = self.__make_pairwise_comparison_sets()
            # get all the distances involving this pair
            for a, b in pairwise:
                values = Numbers()
                for comp_names, param_vals in self.__param_ests.items():
                    if a in comp_names and b in comp_names:
                        values.append(param_vals[param][a] + \
                                    param_vals[param][b])
                
                pairwise_stats[(a,b)] = getattr(values, summary_func)
        else:
            # no additional processing of the distances is required
            
            for comp_names, param_vals in self.__param_ests.items():
                pairwise_stats[comp_names] = param_vals[param]
            
        return pairwise_stats
    
    def getPairwiseDistances(self,summary_function="mean", **kwargs):
        """Return the pairwise distances as a dictionary keyed by (seq1, seq2).
        Convenience interface to getPairwiseParam.
        
        Arguments:
            - summary_function: a string naming the function used for
              estimating param from threeway distances. Valid values are 'mean'
              (default) and 'median'.
        """
        return self.getPairwiseParam('length',summary_function=summary_function,
                                    **kwargs)
    
    def getParamValues(self, param, **kwargs):
        """Returns a Numbers object with all estimated values of param.
        
        Arguments:
            - param: name of a parameter in est_params or 'length'
            - **kwargs: arguments passed to getPairwiseParam"""
        ests = self.getPairwiseParam(param, **kwargs)
        return Numbers(ests.values())
    
    def getTable(self,summary_function="mean", **kwargs):
        """returns a Table instance of the distance matrix.
        
        Arguments:
            - summary_function: a string naming the function used for
              estimating param from threeway distances. Valid values are 'mean'
              (default) and 'median'."""
        d = \
         self.getPairwiseDistances(summary_function=summary_function,**kwargs)
        if not d:
            d = {}
            for s1 in self.__seqnames:
                for s2 in self.__seqnames:
                    if s1 == s2:
                        continue
                    else:
                        d[(s1,s2)] = 'Not Done'
        twoD = []
        for s1 in self.__seqnames:
            row = [s1]
            for s2 in self.__seqnames:
                if s1 == s2:
                    row.append('')
                    continue
                try:
                    row.append(d[(s1,s2)])
                except KeyError:
                    row.append(d[(s2,s1)])
            twoD.append(row)
        T = table.Table(['Seq1 \ Seq2'] + self.__seqnames, twoD, row_ids = True,
                        missing_data = "*")
        return T
    
    def getNewickTrees(self):
        """Returns a list of Newick format trees for supertree methods."""
        trees = []
        for comp_names, param_vals in self.__param_ests.items():
            tips = []
            for name in comp_names:
                tips.append(repr(name)+":%s" % param_vals[name])
            trees.append("("+",".join(tips)+");")
        
        return trees
    
    def writeToFile(self, filename, summary_function="mean", format='phylip',
            **kwargs):
        """Save the pairwise distances to a file using phylip format. Other
        formats can be obtained by getting to a Table.  If running in parallel,
        the master CPU writes out.
        
        Arguments:
            - filename: where distances will be written, required.
            - summary_function: a string naming the function used for
              estimating param from threeway distances. Valid values are 'mean'
              (default) and 'median'.
            - format: output format of distance matrix
        """
        
        if self._on_master_cpu:
             # only write output from 0th node
             table = self.getTable(summary_function=summary_function, **kwargs)
             table.writeToFile(filename, format=format)
    
