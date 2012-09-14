#!/usr/bin/env python
"""Provides a parser to create a Kinemage object from a kinemage file
"""
import re
from cogent.format.mage import Kinemage,MageGroup,MageList,MagePoint,\
    SimplexHeader
from cogent.util.misc import extract_delimited
from cogent.parse.record_finder import LabeledRecordFinder

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"

def MageGroupFromString(line):
    """Returns a new MageGroup, created from a string representation"""
    result = MageGroup([],RecessiveOn=False) 
    trans = {'off':('Off',True),\
            'on':('Off',False),\
            'recessiveon':('RecessiveOn',True),\
            'color':'Color',\
            'radius':'Radius',\
            'nobutton':('NoButton',True),\
            'dominant':('Dominant',True),\
            'lens':('Lens',True),
            'master':'Master',\
            'instance':'Instance',\
            'clone':'Clone'}
    
    #extract all delimited fields: label & KeyWordOptions (master, etc)
    delimited_fields = []
    while 1:
        part = extract_delimited(line,'{','}')
        if part is not None:
            delimited_fields.append(part)
            line = line.replace('{'+part+'}','')
        else:
            break
    #the first one is always the label
    label = delimited_fields[0] 
    #the later ones (starting with 1) are keyword options
    field_idx = 1
    #gather all left-over pieces
    pieces = line.split()
    
    if 'sub' in pieces[0]: #@(sub)group
        result.Subgroup = True
    result.Label = label

    #process all optional pieces
    for piece in pieces[1:]:
        try:
            #here we're finding the key. The value will be '', because it
            #is stored in delimited_fields (accesible by field_idx)
            key,value = piece.split('=')
            setattr(result,trans[key],delimited_fields[field_idx])
            field_idx += 1
        except ValueError:
            setattr(result,trans[piece][0],trans[piece][1])
    return result

def MageListFromString(line):
    """Returns a new MageList, created from a string representation"""
    result = MageList()
    trans = {'off':('Off',True),\
            'on':('Off',False),\
            'color':'Color',\
            'radius':'Radius',\
            'nobutton':('NoButton',True),\
            'angle':'Angle',\
            'width':'Width',\
            'face':'Face',\
            'font':'Font',\
            'size':'Size'}

    line = re.sub('=\s*','=',line)
    label = extract_delimited(line, '{', '}')
    if label is not None:
        pieces = line.replace('{'+label+'}', '').split()
    else:
        pieces = line.split()
   
    style = pieces[0].strip('@')[:-4] #take off the 'list' part
    if style in MageList.KnownStyles:
        result.Style = style
    else:
        raise ValueError,"Unknown style: %s"%(style)
    result.Label = label 
    
    #process all optional pieces
    for piece in pieces[1:]:
        try:
            key,value = [item.strip() for item in piece.split('=')]
            key = trans[key]
        except ValueError: #unpack list of wrong size
            key,value = trans[piece][0],trans[piece][1]
        #KeyError will be raised in case of an unkown key
        setattr(result,key,value)
    return result


def MagePointFromString(line):
    """Constructs a new MagePoint from a one-line string."""
    #handle the label if there is one: note that it might contain spaces
    line = line.strip()
    result = MagePoint()
    label = extract_delimited(line, '{', '}')
    if label:
        pieces = line.replace('{'+label+'}', '').split()
    else:
        pieces = line.split()
    fields = []
    #also have to take into account the possibility of comma-delimited parts
    for p in pieces:
        fields.extend(filter(None, p.split(',')))
    pieces = fields
    result.Label = label
    #get the coordinates and remove them from the list of items
    result.Coordinates = map(float, pieces[-3:])
    pieces = pieces[:-3]
    #parse the remaining attributes in more detail
    result.State = None
    result.Width = None
    result.Radius = None
    result.Color = None
    for attr in pieces:
        #handle radius
        if attr.startswith('r='):  #radius: note case sensitivity
            result.Radius = float(attr[2:])
        #handle single-character attributes
        elif len(attr) == 1:
            result.State = attr
        #handle line width
        elif attr.startswith('width'):
            result.Width = int(attr[5:])
        else:
            #otherwise assume it's a color label
            result.Color = attr
    return result

def _is_keyword(line):
    if line.startswith('@'):
        return True
    return False

KeywordFinder = LabeledRecordFinder(_is_keyword)

def MageParser(infile):
    """MageParser returns a new kinemage object, created from a string repr.

    infile: should be an iterable file object
    
    The MageParser works only on ONE kinemage object, so files containing more
    than one kinemage should be split beforehand. This can easily be adjusted
    if it would be useful in the future.

    The MageParser handles only certain keywords (@kinemage, @text, @caption,
    @___group, @____list) and MagePoints at this point in time. All unkown 
    keywords are assumed to be part of the header, so you can find them in 
    the header information. 
    The lists that are part of the Simplex header are treated as normal lists.
    The 'text' and 'caption' are printed between the header and the first group
    (see cogent.format.mage). All text found after @text keywords is grouped
    and outputted as Kinemage.Text.
    
    WARNING: MageParser should work on all .kin files generated by the code in
    cogent.format.mage. There are no guarantees using it on other kinemage
    files, so always check your output!!!
    """
    text = []
    caption = []
    header = []
    group_pat = re.compile('@(sub)?group')
    list_pat = re.compile('@\w{3,8}list')
    last_group = None
    
    for rec in KeywordFinder(infile):
        first = rec[0]
        other=None
        if len(rec)>1:
            other = rec[1:]
        if first.startswith('@kinemage'): #kinemage
            k = Kinemage()
            count = int(first.replace('@kinemage','').strip())
            k.Count = count
        elif first.startswith('@text'): #text
            if other:
                text.extend(other)
        elif first.startswith('@caption'): #caption
            if other:
                caption.extend(other)
        elif group_pat.match(first): #group
            m = MageGroupFromString(first)
            if m.Subgroup: #subgroup, should be appended to some other group
                if last_group is None:
                    raise ValueError,"Subgroup found before first group"
                last_group.append(m)
            else: #normal group
                k.Groups.append(m)
            last_group = m
        elif list_pat.match(first): #list
            l = MageListFromString(first)
            if other:
                points = [MagePointFromString(p) for p in other]
                l.extend(points)
            last_group.append(l)
        else: #something else
            header.append(first)
            if other:
                header.extend(other)
    
    if text:
        k.Text = '\n'.join(text)
    if caption:
        k.Caption = '\n'.join(caption)
    if header:
        k.Header = '\n'.join(header)
    return k
