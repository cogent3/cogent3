from cogent3.draw.drawable import Drawable
from cogent3.maths.geometry import SimplexTransform

__author__ = "Rahul Ghangas and Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley", "Rahul Ghangas", "Helmut Simon"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class Simplex(Drawable):
    def __init__(self, probabilties, vertex_labels=None, connect_points=False):
        super(Simplex, self).__init__()
        self.vertex_labels = vertex_labels
        self.connect_points = connect_points
        self.probabilities = probabilties

    def _build_fig(self, width=800, height=800, title=None, **kw):
        import plotly.graph_objs as go
        from numpy import array
        from itertools import combinations
        layout = {k:v for k, v in locals().items() if k != 'self' and v}
        self.layout.update(layout)

        probs = array(self.probabilities)
        simplex = SimplexTransform().q
        coords_3d = list(map(lambda x: x @ simplex, probs))
        coords_3d = array(sorted(coords_3d, key=lambda x: x[0]))

        combinations = array(list(combinations(simplex, 2)))
        vertex_labels = self.vertex_labels if self.vertex_labels else list('ACGT')

        data = [
            go.Scatter3d(
                # Draw the edges of the Tetrahedron
                name='Tetrahedron Outline',
                x=combinations[:, :, 0].ravel(),
                y=combinations[:, :, 1].ravel(),
                z=combinations[:, :, 2].ravel(),
                marker=dict(
                    size=4,
                    color='#1f77b4',
                    colorscale='Viridis',
                ),
                line=dict(
                    color='#1f77b4',
                    width=5
                ),
                mode='lines',
                hoverinfo='skip'
            ),
            go.Scatter3d(
                # Draw the edges of the Tetrahedron
                x=simplex[:, 0],
                y=simplex[:, 1],
                z=simplex[:, 2],
                marker=dict(
                    size=4,
                    color='#1f77b4',
                    colorscale='Viridis',
                ),
                text=vertex_labels,
                mode='markers+text',
                hoverinfo='skip'
            ),
            go.Scatter3d(
                # Draw the 3d point and connect them
                x=coords_3d[:, 0], y=coords_3d[:, 1], z=coords_3d[:, 2],
                marker=dict(
                    size=4,
                    color='#1f77b4',
                    colorscale='Viridis',
                ),
                line=dict(
                    color='#CB2D7C',
                    width=3
                )
            )
        ]

        # Layout attributes for plotly
        layout = go.Layout(scene=dict(
            xaxis=dict(
                title='x',
                visible=False,
                range=[min([x[0] for x in simplex]),
                       max([x[0] for x in simplex])],
            ),
            yaxis=dict(
                title='y',
                visible=False,
                range=[min([x[0] for x in simplex]),
                       max([x[0] for x in simplex])],
            ),
            zaxis=dict(
                title='z',
                visible=False,
                range=[min([x[0] for x in simplex]),
                       max([x[0] for x in simplex])]), ),
            autosize=True,

        )
        self.traces.extend(data)
        self._layout = layout
