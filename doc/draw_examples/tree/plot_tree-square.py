"""
Display a Phylogenetic Tree with a Square Dendrogram Style
==========================================================

We use a tree saved in `json` format from a likelihood function analysis of a non-stationary model. The tree was derived such that the the branch lengths are now "ENS".
"""
# %%
from cogent3.app import io


reader = io.load_json()

ens_tree = reader("../../data/GN-tree.json")
fig = ens_tree.get_figure(width=600, height=600)
fig.show()

#%%
# Changing scale bar placement
# ############################

fig.scale_bar = "top right"
fig.show()

#%%
# Colouring a set of edges
# ########################

fig.style_edges("AfricanEl", tip2="Manatee", legendgroup="Afrotheria", 
                line=dict(color="magenta"))
fig.show()

#%%
# With Contemporaneous Tips
# #########################

fig.contemporaneous = True
fig.show()
