"""
Display a Phylogenetic Tree with a Radial Dendrogram Style
==========================================================
"""
# %%
from cogent3.app import io


reader = io.load_json()

ens_tree = reader("../../data/GN-tree.json")
fig = ens_tree.get_figure("radial", width=600, height=600)
fig.show()

#%%
# With Contemporaneous Tips
# #########################

fig.contemporaneous = True
fig.label_pad = 0.23
fig.show()
