"""
Display a Phylogenetic Tree Showing Bootstrap Support
=====================================================

We use a tree saved in `json` format from a 100 replicate bootstrap resamplings. The `show_support=True` argument controls whether or not to display support. The `threshold=0.8` argument indicates only nodes with a support level ≤0.8 will have support text displayed.
"""
# %%
from cogent3.app import io


reader = io.load_json()

tree = reader("../../data/tree-with-support.json")
fig = tree.get_figure(show_support=True, threshold=0.8)
fig.scale_bar = None
fig.show(width=500, height=400)

#%%
# Change the placement of support text
# ####################################
# 
# The support text is positioned relative to the `x`, `y` coordinates of the tree node. Control over support text placement is achieved using the `support_xshift` and `support_yshift` attributes. These are expressed in terms of pixels.
# 
# To place the support text internal to the node, we set the yshift=0 (so at the same y-value of the node) and xshift it to the right using a positive integer.

fig.support_xshift = 15
fig.support_yshift = 0
fig.show(width=500, height=400)
