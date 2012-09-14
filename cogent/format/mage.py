#!/usr/bin/env python
"""
MAGE format writer, especially useful for writing the RNA/DNA simplex.

The MAGE format is documented here:
    ftp://kinemage.biochem.duke.edu/pub/kinfiles/docsDemos/KinFmt6.19.txt

Implementation notes

Mage seems to have a dynamic range of 'only' 2-3 orders of magnitude for balls
and spheres. Consequently, although we have to truncate small values when
writing in the RNA simplex, no additional output accuracy is gained by scaling
everything by a large constant factor (e.g. by 10,000). Consequently, this code
preserves the property of writing the simplex in the interval (0,1) since the
numbers reflect the fractions of A, C, G and U directly.
"""

from __future__ import division
from numpy import array, fabs

from copy import deepcopy
from cogent.util.misc import extract_delimited
from cogent.maths.stats.util import Freqs
from cogent.maths.stats.special import fix_rounding_error

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley", "Sandra Smit", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

def is_essentially(ref_value, value, round_error=1e-14):
    """Returns True if value is within round_error distance from ref_value.

    ref_value -- int or float, number
    value -- int or float, number
    round_error -- limit on deviation from ref_value, default=1e-14

    This function differs from fix_rounding_error:
    - it returns a boolean, not 1 or 0. 
    - it checks the lower bound as well, not only if 
        ref_value < value < ref_value+round_error
    """
    if ref_value-round_error < value < ref_value+round_error:
        return True
    return False

#Note: in addition to the following 24 colors, there's also deadwhite and
#deadblack (force white or black always). The color names cannot be changed.
#MAGE uses an internal color table to look up 4 'steps' in each color 
#corresponding to distance from the camera, and varies the colors depending
#on whether the background is black or white. A @fullrgbpalette declaration
#can be used to reassign colors by number, but it seems to be tricky to
#get it right.
MageColors = [ 'red', 'orange', 'gold', 'yellow', 'lime', 'green', 'sea',
               'cyan', 'sky', 'blue', 'purple', 'magenta', 'hotpink', 'pink',
               'lilac', 'peach', 'peachtint', 'yellowtint', 'greentint',
               'bluetint', 'lilactint', 'pinktint', 'white', 'gray', 'brown']

#use these for 7-point scales
RainbowColors = ['red', 'orange', 'yellow', 'green', 'cyan', 'blue', 'purple']

#use these for matched series of up to 3 items. The first 5 series work pretty
#well, but the last series doesn't really match so only use it if desperate.
MatchedColors = [['red', 'green', 'blue'],
                 ['pinktint', 'greentint', 'bluetint'],
                 ['orange', 'lime', 'sky'],
                 ['magenta', 'yellow', 'cyan'],
                 ['peach', 'sea', 'purple'],
                 ['white', 'gray', 'brown']]
    
class MagePoint(object):
    """Holds information for a Mage point: label, attributes, coordinates."""
    MinRadius = 0.001
    DefaultRadius = MinRadius

    def __init__(self, Coordinates=None, Label=None, Color=None, State=None, 
                 Width=None, Radius=None):
        """Returns a new MagePoint with specified coordinates.

        Usage: m = MagePoint(Coordinates=None, Label=None, Color=None,
                   Attributes=None, Width=None, Radius=None)

        Coordinates is a list of floats: x, y, z
        Label is a text label (warning: does not check length).
        Color is the name of the point's color
        State can be one of:
            P   beginning of new point in vectorlist
            U   unpickable
            L   line to this point
            B   ball at this point in a vectorlist
            S   sphere at this point in a vectorlist
            R   ring at this point in a vectorlist
            Q   square at this point in a vectorlist
            A   unknown (arrow?): something to do with vectorlist
            T   interpret as part of triangle in vectorlist
        Width is an integer that defines pen width
        Radius is a float that defines the point radius for balls, spheres, etc.
        
        If Coordinates is None, the Coordinates will be set to [0,0,0]
        """
        if Coordinates is None:
            self.Coordinates = [0,0,0]
        else:
            self.Coordinates = Coordinates
        self.Label = Label
        self.Color = Color
        self.State = State
        self.Width = Width
        self.Radius = Radius
        
        
    def __str__(self):
        """Prints out current point as string."""
        coords = self.Coordinates
        if len(coords) != 3:
            raise ValueError, "Must have exactly 3 coordinates: found %s" \
                % coords
        pieces = []
        #add the pointID
        if self.Label is not None:
            pieces.append('{%s}' % self.Label)
        #collect attributes and add them if necessary
        attr = []
        if self.State:
            attr.append(self.State)
        if self.Color:
            attr.append(self.Color)
        if self.Width:
            attr.append('width%s' % self.Width)
        r = self.Radius
        #NOTE: as of version 6.25, MAGE could not handle radius < 0.001
        if r is not None:
            if r > self.MinRadius:
                attr.append('r=%s' % r)
            else:
                attr.append('r=%s' % self.DefaultRadius)
        if attr:
            pieces += attr
        #add the coordinates
        pieces.extend(map(str, coords))
        #return the result
        return ' '.join(pieces)

    def __cmp__(self, other):
        """Checks equality/inequality on all fields.
        
        Order: Coordinates, Label, Color, State, Width, Radius"""
        try:
            return cmp(self.__dict__, other.__dict__)
        except AttributeError:  #some objects don't have __dict__...
            return cmp(self.__dict__, other)

    def _get_coord(coord_idx):
        """Gets the coordinate asked for; X is first, Y second, Z third coord
        """
        def get_it(self):
            return self.Coordinates[coord_idx]
        return get_it
    
    def _set_coord(coord_idx):
        """Sets the given coordinate; X is first, Y second, Z third coord"""
        def set_it(self,val):
            self.Coordinates[coord_idx] = val
        return set_it

    _lookup = {'x':(_get_coord(0),_set_coord(0)),
               'y':(_get_coord(1),_set_coord(1)),
               'z':(_get_coord(2),_set_coord(2))}

    X = property(*_lookup['x'])
    Y = property(*_lookup['y'])
    Z = property(*_lookup['z'])
       
    def toCartesian(self, round_error=1e-14):
        """Returns new MagePoint with UC,UG,UA coordinates from A,C,G(,U).

        round_error -- float, accepted rounding error

        x=u+c, y=u+g, z=u+a

        This will only work for coordinates in a simplex, where all of them 
        (inlcuding the implicit one) add up to one.
        """
        # all values have to be between 0 and 1
        a = fix_rounding_error(self.X)
        c = fix_rounding_error(self.Y)
        g = fix_rounding_error(self.Z)
        #u = fix_rounding_error(1-a-c-g)
        if is_essentially(1, a+c+g, round_error=round_error):
            u = 0.0
        else:
            u = fix_rounding_error(1-a-c-g)
        
        for coord in [a,c,g,u]:
            if not 0 <= coord <= 1:
                raise ValueError,\
                "%s is not in unit simplex (between 0 and 1)"%(coord)
        cart_x, cart_y, cart_z = u+c, u+g, u+a
        result = deepcopy(self)
        result.Coordinates = [cart_x, cart_y, cart_z]
        return result
    
    def fromCartesian(self):
        """Returns new MagePoint with A,C,G(,U) coordinates from UC,UG,UA.

        From UC,UG,UA to A,C,G(,U).

        This will only work when the original coordinates come from a simplex,
        where U+C+A+G=1
        """
        # x=U+C, y=U+G, z=U+A, U+C+A+G=1
        # U=(1-x-y-z)/-2
        x,y,z = self.X, self.Y, self.Z
        u = fix_rounding_error((1-x-y-z)/-2)
        a, c, g = map(fix_rounding_error,[z-u, x-u, y-u])
        result = deepcopy(self)
        result.Coordinates = [a, c, g]
        return result
        
def MagePointFromBaseFreqs(freqs, get_label=None, get_color=None, \
    get_radius=None):
    """Returns a MagePoint from an object with counts for the bases.
    
    get_label should be a function that calculates a label from the freqs.
    If get_label is not supplied, checks freqs.Label, freqs.Species, freqs.Id,
    freqs.Accession, and freqs.Name in that order. If get_label fails or none
    of the attributes is found, no label is written.

    get_color should be a function that calculates a color from the freqs. 
    Default is no color (i.e. the point has the color for the series), which
    will also happen if get_color fails.

    get_radius is similar to get_color.
    """
    label = None
    if get_label:
        try:
            label = get_label(freqs)
        except:
            pass    #label will be assigned None below
    else:
        for attr in ['Label', 'Species', 'Id', 'Accession', 'Name']:
            if hasattr(freqs, attr):
                label = getattr(freqs, attr)
                #keep going if the label is empty
                if label is not None and label != '':
                    break
    if not label and label != 0:
        label = None
    if get_color:
        try:
            color = get_color(freqs)
        except:
            color=None
    else:
        if hasattr(freqs, 'Color'):
            color = freqs.Color
        else:
            color = None
            
    if get_radius:
        try:
            radius = get_radius(freqs)
        except:
            radius=None
    else:
        if hasattr(freqs, 'Radius'):
            try:
                radius = float(freqs.Radius)
            except:
                radius = None
        else:
            radius = None
            
    relevant = Freqs({'A':freqs.get('A',0), 'C':freqs.get('C',0), 
        'G':freqs.get('G',0), 'U':freqs.get('U',0) or  freqs.get('T',0)})
    relevant.normalize()
    return MagePoint((relevant['A'],relevant['C'],relevant['G']), Label=label,\
        Color=color, Radius=radius)

class MageList(list):
    """Holds information about a list of Mage points.
    """
    KnownStyles = ['vector', 'dot', 'label', 'word', 'ball', 'sphere', \
                   'triangle', 'ribbon', 'arrow', 'mark', 'ring', 'fan']
    BooleanOptions = ['Off', 'NoButton']
    KeywordOptions = ['Color','Radius','Angle','Width','Face','Font','Size']
    def __init__(self, Data=None, Label=None, Style=None, Color=None, \
        Off=False, NoButton=False, Radius=None, Angle=None, Width=None, \
        Face=None, Font=None, Size=None):
        """Returns new MageList object. Must initialize with list of MagePoints.
        """
        if Data is None:
            super(MageList,self).__init__([])
        else:
            super(MageList, self).__init__(Data)
        self.Label = Label
        self.Style = Style
        self.Color = Color
        self.Off = Off
        self.NoButton = NoButton
        self.Radius = Radius
        self.Angle = Angle
        self.Width = Width
        self.Face = Face
        self.Font = Font
        self.Size = Size

    def __str__(self):
        """Returns data as a mage-format string, one point to a line."""
        pieces = []
        curr_style = self.Style or 'dot'    #dotlists by default
        pieces.append('@%slist' % curr_style.lower())
        if self.Label:
            pieces.append('{%s}' % self.Label)
        for opt in self.BooleanOptions:
            if getattr(self, opt):
                pieces.append(opt.lower())
        for opt in self.KeywordOptions:
            curr = getattr(self, opt)
            if curr:
                pieces.append('%s=%s' % (opt.lower(), curr))
        lines = [' '.join(pieces)] + [str(d) for d in self]
        return '\n'.join(lines)

    def toArray(self, include_radius=True):
        """Returns an array of the three coordinates and the radius for all MP
        """
        elem = []
        if include_radius:
            for point in self:
                coords = point.Coordinates
                radius = point.Radius or self.Radius
                if radius is None:
                    raise ValueError,\
                        "Radius is not set for point or list, %s"%(str(point))
                else:
                    elem.append(coords + [radius])
        else:
            for point in self:
                coords = point.Coordinates
                elem.append(coords)
        return array(elem)
    
    def iterPoints(self):
        """Iterates over all points in the list"""
        for point in self:
            yield point

    def toCartesian(self, round_error=1e-14):
        """Returns new MageList where points are in Cartesian coordinates"""
        #create new list
        result = MageList()
        #copy all attributes from the original
        result.__dict__ = self.__dict__.copy()
        #fill with new points
        for point in self.iterPoints():
            result.append(point.toCartesian(round_error=round_error))
        return result

    def fromCartesian(self):
        """Returns new MageList where points are in ACG coordinates"""
        #create new list
        result = MageList()
        #copy all attributes from the original
        result.__dict__ = self.__dict__.copy()
        #fill with new points
        for point in self.iterPoints():
            result.append(point.fromCartesian())
        return result

class MageGroup(list):
    """Holds information about a MAGE Group, which has a list of lists/groups.

    Specificially, any nested MageGroups will be treated as MAGE @subgroups, 
    while any dotlists, vectorlists, etc. will be treated directly.
    """
    BooleanOptions = ['Off', 'NoButton', 'RecessiveOn', 'Dominant', 'Lens']
    KeywordOptions = ['Master', 'Instance', 'Clone']
    Cascaders = ['Style', 'Color','Radius','Angle','Width','Face','Font','Size']
    
    def __init__(self, Data=None, Label=None, Subgroup=False, Off=False, \
        NoButton=False, RecessiveOn=True, Dominant=False, Master=None,   \
        Instance=None, Clone=None, Lens=False, Style=None, Color=None,   \
        Radius=None, Angle=None, Width=None, Face=None, Font=None, Size=None):
        """Returns new MageList object. Must initialize with list of MagePoints.
        """
        if Data is None:
            super(MageGroup,self).__init__([])
        else:
            super(MageGroup,self).__init__(Data)
        self.Label = Label
        self.Subgroup = Subgroup
        self.Off = Off
        self.NoButton = NoButton
        self.RecessiveOn = RecessiveOn
        self.Dominant = Dominant
        self.Master = Master
        self.Instance = Instance
        self.Clone = Clone
        self.Lens = Lens
        #properties of lists. Group settings will be expressed when Kinemage
        #object is printed
        self.Style = Style
        self.Color = Color
        self.Radius = Radius
        self.Angle = Angle
        self.Width = Width
        self.Face = Face
        self.Font = Font
        self.Size = Size
 
    def __str__(self):
        """Returns data as a mage-format string, one point to a line.
        
        Children are temporarily mutated to output the correct string. Changes
        are undone afterwards, so the original doesn't change.
        """
        pieces = []
        if self.Subgroup:
            pieces.append('@subgroup')
        else:
            pieces.append('@group')
        if self.Label:
            pieces.append('{%s}' % self.Label)
        for opt in self.BooleanOptions:
            if getattr(self, opt):
                pieces.append(opt.lower())
        for opt in self.KeywordOptions:
            curr = getattr(self, opt)
            if curr:
                pieces.append('%s={%s}' % (opt.lower(), curr))
        for d in self:
            if isinstance(d, MageGroup):
                d.Subgroup = True
        # cascade values 1 level down (to either child lists or groups)
        changed = {}
        for attr in self.Cascaders:
            val = getattr(self,attr)
            if val is not None:
                changed[attr] = []
                for item in self: #either group or list
                    if getattr(item,attr) is None:
                        setattr(item,attr,val)
                        changed[attr].append(item)
        #gather all lines
        lines = [' '.join(pieces)] + [str(d) for d in self]
        #reset cascaded values
        for k,v in changed.items():
            for l in v:
                setattr(l,k,None)
        return '\n'.join(lines)
    
    def iterGroups(self):
        """Iterates over all groups in this group"""
        for item in self:
            if isinstance(item,MageGroup):
                yield item
                for j in item.iterGroups():
                    yield j
                
    def iterLists(self):
        """Iterates over all lists in this group"""
        for item in self:
            if isinstance(item,MageList):
                yield item
            else:
                for j in item.iterLists():
                    yield j
    
    def iterGroupsAndLists(self):
        """Iterates over all groups and lists in this group
        
        Groups and Lists have multiple elements in common. You might want to 
        change a property for both groups and lists. This iterator doesn't 
        distinguish between them, but returns them all. 
        """
        for i in self:
            if isinstance(i,MageGroup):
                yield i
                for j in i.iterGroupsAndLists():
                    yield j
            elif isinstance(i,MageList):
                yield i
        
    def iterPoints(self):
        """Iterates over all points in this group"""
        for item in self.iterLists():
            for p in item.iterPoints():
                yield p

    def toCartesian(self, round_error=1e-14):
        """Returns a new MageGroup where all points are Cartesian coordinates
        """
        #create an empty group
        result = MageGroup()
        #copy the attributes of the original
        result.__dict__ = self.__dict__.copy()
        #fill with new groups and lists
        for item in self:
            result.append(item.toCartesian(round_error=round_error))
        return result
   
    def fromCartesian(self):
        """Returns a new MageGroup where all points are ACG coordinates
        """
        #create an empty group
        result = MageGroup()
        #copy the attributes of the original
        result.__dict__ = self.__dict__.copy()
        #fill with new groups and lists
        for item in self:
            result.append(item.fromCartesian())
        return result

class MageHeader(object):
    """Defines the header for a kinemage.

    For now, just returns the text string it was initialized with.
    """
    def __init__(self, data):
        """Returns a new MageHeader object."""
        self.Data = data

    def __str__(self):
        """Writes the header information out as a string."""
        return str(self.Data)

class SimplexHeader(MageHeader):
    """Defines the header for the RNA simplex.

    May have a bunch of scaling and coloring options later.
    """
    def __init__(self, *args, **kwargs):
        """Returns a new SimplexHeader object.

        May have options added later, but none for now.
        """
        self.Data = \
"""
@viewid {oblique}
@zoom 1.05
@zslab 467
@center 0.500 0.289 0.204
@matrix
-0.55836 -0.72046 -0.41133  0.82346 -0.42101 -0.38036  0.10085 -0.55108
0.82833

@2viewid {top}
@2zoom 0.82
@2zslab 470
@2center 0.500 0.289 0.204
@2matrix
-0.38337  0.43731 -0.81351  0.87217 -0.11840 -0.47466 -0.30389 -0.89148
-0.33602

@3viewid {side}
@3zoom 0.82
@3zslab 470
@3center 0.500 0.289 0.204
@3matrix
-0.49808 -0.81559 -0.29450  0.86714 -0.46911 -0.16738 -0.00164 -0.33875
0.94088

@4viewid {End-O-Line}
@4zoom 1.43
@4zslab 469
@4center 0.500 0.289 0.204
@4matrix
 0.00348 -0.99984 -0.01766  0.57533 -0.01244  0.81784 -0.81792 -0.01301
0.57519

@perspective
@fontsizelabel 24
@onewidth
@zclipoff
@localrotation  1 0 0 .5 .866 0 .5 .289 .816

@group {Tetrahedron}
@vectorlist  {Edges}  color=white     nobutton
P   {0 0 0} 0 0 0
     0.5 0 0
  {1 0 0} 1 0 0
     0.5  0.5  0
  {0 1 0} 0 1 0
P    0 0 0
    0  0.5  0
  {0 1 0} 0 1 0
   0  0.5  0.5
  {0 0 1} 0 0 1
P    0 0 0
   0 0 0.5
  {0 0 1} 0 0 1
   0.5  0 0.5
  {1 0 0} 1 0 0

@labellist {labels} color=white    nobutton
  {U}  0  0  0
   {A}  1.1    0  0
   {C}  0    1.05   0
   {G}   0  0   1.08


@group  {Lines}
@vectorlist {A=U&C=G} color= green    off
 P   0   0.5   0.5
     .1  .4  .4
     .25 .25 .25
     .4  .1  .1
 L   0.500, 0.000, 0.000

@vectorlist {A=G&C=U} color= red    off
 P   0.5   0   0.5
     .25 .25 .25
 L   0, 0.500, 0.000

@vectorlist {A=C&G=U} color= red    off
 P   0.5   0.5   0
     .25 .25 .25
 L   0.000, 0.000, 0.500
"""

class CartesianSimplexHeader(MageHeader):
    """Defines the header for the RNA simplex: coordinates in Cartesian system.

    May have a bunch of scaling and coloring options later.
    """
    def __init__(self, *args, **kwargs):
        """Returns a new CartesianSimplexHeader object.

        May have options added later, but none for now.
        """
        self.Data = \
"""
@1viewid {Oblique}
@1span 2.2
@1zslab 200.0
@1center 0.5 0.5 0.5
@1matrix 0.2264 0.406 0.8854 0.6132 0.6468 -0.4535 -0.7567 0.6456 -0.1025

@2viewid {Down UA axis}
@2span 2.0
@2zslab 200.0
@2center 0.5 0.5 0.5
@2matrix 1 0 0 0 1 0 0 0 1

@3viewid {Down UC axis}
@3span 2.0
@3zslab 200.0
@3center 0.5 0.5 0.5
@3matrix 0 0 1 0 1 0 -1 0 0

@4viewid {Down UG axis}
@4span 2.0
@4zslab 200.0
@4center 0.5 0.5 0.5
@4matrix 0 -1 0 0 0 1 -1 0 0

@5viewid {Cube}
@5span 2.0
@5zslab 200.0
@5center .5 0.5 0.5
@5matrix 0.956 -0.0237 0.2923 0.2933 0.0739 -0.9532 0.001 0.997 0.0776

@group {Tetrahedron}
@vectorlist {Edges} nobutton color= white
{Edges}P 1 1 1
1 0 0
0 0 1
0 1 0
1 0 0
1 1 1
0 0 1
1 1 1
0 1 0

@labellist {labels} nobutton color= white
{U}1.05 1.05 1.05
{A}-0.05 -0.05 1.05
{C}1.05 -0.05 -0.05
{G}-0.05 1.05 -0.05

@group {Lines}
@vectorlist {A=U&C=G} color= green
{A=U&C=G}P 0.5 0.5 1
0.5 0.5 0
@vectorlist {A=G&C=U} color= red
{A=G&C=U}P 0 0.5 0.5
1 0.5 0.5
@vectorlist {A=C&G=U} color= red
{A=C&G=U}P 0.5 0 0.5
0.5 1 0.5

@group {Cube} off
@vectorlist nobutton color=yellow 
{"}P 0 0 0
{"}1 0 0
{"}1 1 0
{"}0 1 0
{"}0 0 0
{"}P 0 0 1
{"}1 0 1
{"}1 1 1
{"}0 1 1
{"}0 0 1
{"}P 0 0 0
{"}0 0 1
{"}P 1 0 0
{"}1 0 1
{"}P 1 1 0
{"}1 1 1
{"}P 0 1 0
{"}0 1 1

@labellist {Cube labels} color= red off
{0,0,0}-0.1 -0.1 -0.1
{1,1,1}1.1 1.1 1.1
{0,0,1} -.1 -.1 1.1
{0,1,1} -.1 1.1 1.1
{0,1,0} -.1 0.9 0.05
{1,0,1} 1.05 -.1 0.95
{1,1,0} 1.1 1.1 -.1
{1,0,0} 1.1 -.1 -.1
{x++}1.03 0.5 0.5
{y++}0.5 1.1 0.475
{z++}0.5 0.5 1.05
{x--}-0.075 0.5 0.5
{y--}0.475 -0.15 0.5
{z--}0.5 0.5 -0.08

"""

class Kinemage(object):
    """Stores information associated with a kinemage: header, caption, groups.

    A single file can have multiple kinemages.
    """
    def __init__(self,Count=None,Header=None,Groups=None,Caption=None,\
        Text=None):
        """Returns a new Kinemage object.

        Usage: k = Kinemage(count, header, groups, caption=None, text=None)
        """
        self.Count = Count  #integer, required for MAGE
        self.Header = Header
        self.Groups = Groups or []
        self.Caption = Caption
        self.Text = Text

    def __str__(self):
        """String representation suitable for writing to file."""
        if not self.Count:
            raise ValueError, "Must set a count to display a kinemage."
        pieces = ['@kinemage %s' % self.Count]
        if self.Header:
            pieces.append(str(self.Header))
        if self.Text:
            pieces.append('@text')
            pieces.append(str(self.Text))
        if self.Caption:
            pieces.append('@caption')
            pieces.append(str(self.Caption))
        for g in self.Groups:
            pieces.append(str(g))
        return '\n'.join(pieces)

    def iterGroups(self):
        """Iterates over all groups in Kinemage object"""
        for i in self.Groups:
            yield i
            for j in i.iterGroups():
                yield j

    def iterLists(self):
        """Iterates over all lists in Kinemage object"""
        for gr in self.Groups:
            for l in gr.iterLists():
                yield l

    def iterPoints(self):
        """Iterates over all points in Kinemage object"""
        for gr in self.Groups:
            for p in gr.iterPoints():
                yield p
   
    def iterGroupsAndLists(self):
        """Iterates over all groups and lists in Kinemage object"""
        for gr in self.Groups:
           yield gr
           for j in gr.iterGroupsAndLists():
               yield j

    def toCartesian(self, round_error=1e-14):
        """Returns a new Kinemage where all coordinates are UC,UG,UA"""
        result = deepcopy(self)
        result.Groups = []
        for item in self.Groups:
            result.Groups.append(item.toCartesian(round_error=round_error))
        return result

    def fromCartesian(self):
        """Returns a new Kinemage where all coordinates are ACG again"""
        result = deepcopy(self)
        result.Groups = []
        for item in self.Groups:
            result.Groups.append(item.fromCartesian())
        return result
