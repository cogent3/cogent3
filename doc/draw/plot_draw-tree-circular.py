"""
Display a Phylogenetic Tree with a Circular Dendrogram Style
============================================================
"""
# %%
from cogent3.app import io

reader = io.load_json()

ens_tree = reader("../data/GN-tree.json")
fig = ens_tree.get_figure("circular", width=600, height=600)
fig.show(renderer="sphinx_gallery")

#%%
# Colouring a set of edges
# ########################

fig.style_edges("AfricanEl", tip2="Manatee", legendgroup="Afrotheria",
                line=dict(color="magenta"))
fig.show(renderer="sphinx_gallery")

#%%
# With Contemporaneous Tips
# #########################

fig.contemporaneous = True
fig.label_pad = 0.23
fig.show(renderer="sphinx_gallery")

