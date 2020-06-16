"""
Evaluating coevolution
======================

A method on the alignment provides an interface to the simpler (and yet robust and fast) methods for estimating coevolution. The default measure is normalised mutual information (NMI).
"""
#%%
# Display coevolution as a heatmap
# ################################

from cogent3 import load_aligned_seqs


aln = load_aligned_seqs("../../data/brca1.fasta", moltype="dna")
aln = aln.no_degenerates(motif_length=3)
aln = aln.get_translation()
aln = aln[:100]  # for compute speed in testing the documentation
coevo = aln.coevolution(show_progress=False, drawable="heatmap")
coevo.show()

#%%
# Display coevolution scores as a Violin plot
# ###########################################

coevo = aln.coevolution(show_progress=False, drawable="violin")
coevo.show(width=300)

#%%
# Display coevolution scores as a Boxplot
# #######################################

coevo = aln.coevolution(show_progress=False, drawable="box")
coevo.show(width=300)
