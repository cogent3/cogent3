#!/usr/bin/env python
"""Parses Agilent Microarray output spreadsheet.
"""
from cogent.parse.record_finder import LabeledRecordFinder

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

def MicroarrayParser(lines):
    """Returns tuple: ([ProbeNames],[GeneNames],[LogRatios]) for all dots in
    microarray file.
    """
    probe_names = []
    gene_names = []
    log_ratios = []
    #Make sure lines is not empty
    if lines:
        #Get the block of lines that starts with FEATURES
        features_record = LabeledRecordFinder(\
            lambda x: x.startswith('FEATURES'))
        features_block = list(features_record(lines))
        #Discard first block
        features_block = features_block[1]
        #Get the indices of GeneName and LogRatio from the block
        features_list = features_block[0].split('\t')
        probe_index = features_list.index('ProbeName')
        gene_index = features_list.index('GeneName')
        log_index = features_list.index('LogRatio')
        #Get the lists for GeneName and LogRatio
        for line in features_block[1:]:
            temp = line.split('\t')
            probe_names.append(temp[probe_index].upper())
            gene_names.append(temp[gene_index].upper())
            log_ratios.append(float(temp[log_index]))
    return (probe_names, gene_names, log_ratios)
    
