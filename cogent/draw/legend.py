#!/usr/bin/env python

from __future__ import division
from cogent.core import moltype, annotation

from matplotlib.collections import PatchCollection
from matplotlib.text import Text
from matplotlib.transforms import Affine2D
from cogent.draw.rlg2mpl import Group, Drawable, figureLayout
from cogent.draw.linear import Display, DisplayPolicy

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class Legend(Drawable):
    """A class for drawing a legend for a display policy
    
    Arguments:
        - policy: a reference to a Display policy class"""
    def __init__(self, policy = DisplayPolicy):
        self.policy = policy
    
    def _makeSampleSequence(self, feature_type):
        seq = moltype.DNA.makeSequence('aaaccggttt' * 7)
        v = seq.addAnnotation(annotation.Feature,
                feature_type, feature_type, [(2,3)])
        v = seq.addAnnotation(annotation.Feature,
                feature_type, feature_type, [(7,18)])
        v = seq.addAnnotation(annotation.Feature,
                feature_type, feature_type, [(20,70)])
        return seq
        
    def populateAxes(self, ax, columns = 3):
        """ Returns the legend as a matplotlib artist
        Arguments:
            - columns: the number of columns of feature / representation
              pairs
        """
        ax.set_xlim(0, 600)
        ax.set_ylim(-800, 50)
        result = []
        x = y = 0
        for track in self.policy()._makeTrackDefns():
            if track.tag is None or track.tag=="Graphs":
                continue
            ax.text(10, y*30, track.tag)
            y -= 1
            for feature in track:
                seq = self._makeSampleSequence(feature)
                display = Display(seq,
                        policy = self.policy,
                        min_feature_height = 10,
                        show_code = False,
                        pad = 0,)
                sample = display.makeArtist()
                #trans = sample.get_transform()
                #offset = Affine2D()
                #offset.translate(x*600+20 / columns, y*30)
                sample.translate(x*600/columns+10, y*30)
                ax.add_artist(sample)
                ax.text(x*600/columns+90, y*30, feature)
                x += 1
                if x % columns == 0:
                    x = 0
                    y -= 1
            if x:
                x = 0
                y -= 1
            ax.axhline((y+.7)*30)
    
    def makeFigure(self, margin=0, default_aspect=1.3, **kw):
        kw['margin'] = margin
        kw['default_aspect'] = default_aspect
        (width, height), posn, kw = figureLayout(leftovers=True, **kw)
        fig = self._makeFigure(width, height)
        ax = fig.add_axes(posn, adjustable="datalim",
            frame_on=False, xticks=[], yticks=[])
        g = self.populateAxes(ax, **kw)
        return fig
        
if __name__ == '__main__':
    Legend().showFigure()
    
