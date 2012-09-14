#!/usr/bin/env python
"""Provides Histogram, which bins arbitrary objects.
"""
from cogent.maths.stats.util import Freqs

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"

class Histogram(object):
    """Stores a list of _bins and a list of corresponding objects.

    It contains always a similar number of bins and value lists.
    A value list can contain multiple objects.
    _bins should implement __contains__.
    """
    Multi = False
    Mapping = None
    
    def __init__(self, data='', bins=None, Mapping=None, Multi=None):
        """Returns a new Histogram object.

        Data is any sequence of data
        bins is a list of objects that implement __contains__
        Mapping is a function to be applied on each object in data. Gives 
            you the option to order the objects according to some value
            e.g. make a distribution of sequence lengths.
            def seq_length(s): return len(s)
            h = histogram(Data=[sequences],_bins=[Spans],Mapping=seq_length)
        Multi determins whether an object might end up in multiple bins.
            False by default.
        
        All parameters other than Data are None by default. This gives you 
        the opportunity to subclass Histogram and store _bins, Multi, and 
        Mapping as class data.
        """
        if bins is not None:
            self._bins = bins
        if Mapping is not None:
            self.Mapping = Mapping
        if Multi is not None:
            self.Multi = Multi
        self.clear()
        self(data)


    def __iter__(self):
        """Iterates through Bin, Values pairs"""
        return iter(zip(self._bins, self._values))

    def __call__(self, data):
        """Puts all objects in data in the bins.

        Keeps old data that was already in histogram, so updates the data.
        """
        function = self.Mapping
        if function:
            transformed = map(function, data)
        else:
            transformed = data
            
        for t, d in zip(transformed, data):
            in_bin = False
            for bin, values in self:
                if t in bin:
                    values.append(d)
                    in_bin = True
                    if not self.Multi:
                        break
                else:
                    continue
            if not in_bin:
                self.Other.append(d)
    
    def clear(self):
        """Erases all data in the bins.
        
        Note: creates new, empty lists, so if you have references to the
        old data elsewhere in your code they won't be updated.
        """
        self._values = [[] for i in self._bins]
        self.Other = []
    
    def __str__(self):
        """Returns string representation of the histogram"""
        result = []
        for bin, values in self:
            result.append(str(bin)+'\t'+str(values))
        return '\n'.join(result)
       
    def toFreqs(self):
        """Returns a Freqs object based on the histogram.

        Labels of Freqs will be _bins converted into strings
        Values of Freqs will be the number of objects in a Bin
        """
        result = Freqs()
        for bin,values in self:
            result[str(bin)] = len(values)
        return result
