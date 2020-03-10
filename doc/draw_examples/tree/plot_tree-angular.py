"""
Display a Phylogenetic Tree with a Angular Dendrogram Style
===========================================================

This is a left-right style. You'll note that there's overlap of edges at the bottom -- a known issue with this display style.
"""
# %%
from cogent3.app import io


reader = io.load_json()

ens_tree = reader("../../data/GN-tree.json")
fig = ens_tree.get_figure(style="angular", width=600, height=600)
fig.show()

#%%
# With Contemporaneous Tips
# #########################

fig.contemporaneous = True
fig.show()
