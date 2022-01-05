"""
Counting gaps per sequence
==========================

We have several different ways of counting sequence gaps, and of visualising the results. By default, the `count_gaps_per_seq()` method returns a matrix of counts without the ability to visualise the results. When setting the argument `unique=True`, the counts are for gaps uniquely induced by each sequence. This can be a useful indicator of highly divergent sequences.
"""

#%%
from cogent3 import load_aligned_seqs


aln = load_aligned_seqs('../../data/brca1.fasta', moltype='dna')

counts = aln.count_gaps_per_seq(unique=True)
counts

#%%
# Plotting counts of unique gaps
# ##############################
# 
# There are three plot types supported. In all cases, placing the mouse pointer over a data point will show hover text with the sequence name.

#%%
# Displaying unique gaps as a bar chart
# *************************************

counts = aln.count_gaps_per_seq(unique=True, drawable='bar')
counts.show(width=500)
#%%
# Displaying unique gaps as a violin plot
# ***************************************

counts = aln.count_gaps_per_seq(unique=True, drawable='violin')
counts.show(width=300, height=500)
#%%
# Displaying unique gaps as a box plot
# ************************************

counts = aln.count_gaps_per_seq(unique=True, drawable='box')
counts.show(width=300, height=500)
