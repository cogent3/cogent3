#!/usr/bin/env python
"""Methods for using trees to generate microarray data.

"""
from cogent.seqsim.tree import RangeNode
from cogent.util.array import mutate_array

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

class MicroarrayNode(RangeNode):

    def __init__(self, Length=0, Array=None, *args, **kwargs):
        """Returns new MicroarrayNode object.

        Length:   float giving the branch length (sd to add to data)
        Array:          array of float giving the expression vector, or None
        Additional args for superclass:
        Name:           usually a text label giving the name of the node
        LeafRange:      range of leaves that the node spans
        Id:             unique numeric identifier from the node
        Children:       list of Node objects specifying the children
        Parent:         Node object specifying the parent
        """
        RangeNode.__init__(self, *args, **kwargs)
        self.Length = Length
        self.Array = Array

    def setExpression(self, vec):
        """Sets expression in self and all children.
        
        WARNING: Will overwrite existing array with new array passed in.
        Expects vec to be a floating-point Numeric array.
        """
        #if it's the root or zero branch length, just set the vector
        if not self.Length:
            self.Array = vec.copy()
        else:
            self.Array = mutate_array(vec, self.Length)

        curr_children = self.Children
        curr_arrays = [self.Array] * len(curr_children)

        while len(curr_children):
            new_children = []
            new_arrays = []

            for c,a in zip(curr_children, curr_arrays):
                if not c.Length:
                    c.Array = a.copy()
                else:
                    c.Array = mutate_array(a, c.Length)

                new_children.extend(c.Children)
                new_arrays.extend([c.Array] * len(c.Children))

            curr_children = new_children
            curr_arrays = new_arrays

