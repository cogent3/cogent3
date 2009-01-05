#!/usr/bin/env python

from cogent import DNA
from cogent.core import alignment, alphabet, annotation
from cogent.draw import *

from reportlab.platypus import Paragraph, SimpleDocTemplate, KeepTogether, Spacer
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4
from xml.sax.saxutils import escape

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
                    "Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.3.0.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

WIDTH = A4[0] - 40

def makeSampleSequence():
    seq = DNA.makeSequence('aaaccggttt' * 10)
    v = seq.addAnnotation(annotation.Feature, 'exon', 'exon', [(20,35)])
    v = seq.addAnnotation(annotation.Feature, 'repeat_unit', 'repeat_unit', [(39,49)])
    v = seq.addAnnotation(annotation.Feature, 'repeat_unit', 'rep2', [(49,60)])
    return seq

def makeSampleAlignment():
    # must be an esier way to make an alignment of annotated sequences!
    from cogent.align.align import global_pairwise, make_dna_scoring_dict
    DNA = make_dna_scoring_dict(10, -8, -8)
    seq1 = makeSampleSequence()[:-2]
    seq2 = makeSampleSequence()[2:]
    seq1.Name = 'FAKE01'
    seq2.Name = 'FAKE02'
    names = (seq1.getName(), seq2.getName())
    align = global_pairwise(seq1, seq2, DNA, 2, 1)
    align.addAnnotation(annotation.Variable, 'redline', 'align', [((0,15),1),((15,30),2),((30,45),3)])
    align.addAnnotation(annotation.Variable, 'blueline', 'align', [((0,15),1.5),((15,30),2.5),((30,45),3.5)])
    return align
    
seq = makeSampleSequence()
a = seq.addAnnotation(annotation.Variable, 'blueline', 'seq', [((0,15),1),((15,30),2),((30,45),3)])
v = seq.addAnnotation(annotation.Feature, 'gene', 'gene', [(0,15),(20,35),(40,55)])
b = v.addAnnotation(annotation.Variable, 'redline', 'feat', [((0,15),1.5),((15,30),2.5),((30,45),3.5)])

align = makeSampleAlignment()

title_style = ParagraphStyle('normal')
story = []

def fig(msg, seq_display, width=WIDTH, **kw):
    kw['total_width'] = width
    kw['margin'] = 10
    msg = '<para>' + escape(msg + ' ' + repr(kw)) + '</para>'
    story.append(KeepTogether([
        seq_display.asDrawing(**kw),
        Paragraph(msg, ParagraphStyle(title_style, leftIndent=10)),
    ]))
    story.append(Spacer(30, 30))

doc = SimpleDocTemplate("draw.test.pdf",
        leftMargin=10, rightMargin=10, pagesize=A4)

for obj in [seq, align]:
    disp = Display(obj, min_feature_height=10)
    fig('small', disp, width=200)
    fig('normal', disp)
    fig('colour', Display(obj, colour_sequences=True))
    fig('wrapped', disp, rowlen=50)
    fig('wrappped slice of display', disp[10:50], rowlen=20)
    fig('slice of display', disp[10:50])
    fig('slice of %s' % type(obj).__name__, Display(obj[10:50]))

# LEGEND
story.append(Spacer(20, 20))
story.append(KeepTogether([
            Legend().asDrawing(WIDTH),
            Paragraph('<para>Legend with default DisplayPolicy</para>',
                     ParagraphStyle(title_style, leftIndent=10)),
          ]))
story.append(Spacer(20, 20))

# DOTPLOT

from cogent.draw.dotplot import Display2D
fig('2d', Display2D(seq, seq))

#fig('reversed', Display(seq[50:10]), 500)
# no good because seqs slice like lists: ie len==0
# complement() probably doesn't know about annotations either.

# TREES
def tree_fig(msg, seq_display, **kw):
    msg = '<para>' + escape(msg + ' ' + repr(kw)) + '</para>'
    story.append(KeepTogether([
            seq_display.asDrawing(width=WIDTH, height=200, **kw),
            Paragraph(msg, ParagraphStyle(title_style)),
        ]))
    story.append(Spacer(20, 20))

from cogent import LoadTree
t = LoadTree(treestring="((A:1,B:2)ab:3,((C:4,D:5)cd:5.5,E:6)cde:7,F:2)")
from cogent.draw.dendrogram import *
for klass in [UnrootedDendrogram, SquareDendrogram, ContemporaneousDendrogram, ShelvedDendrogram, StraightDendrogram,
    ContemporaneousStraightDendrogram]:
    tree_fig(klass.__name__ + ' Color', klass(t), shade_param="length")
    tree_fig(klass.__name__ + ' Label', klass(t), show_params=["length"])

def callback(edge):
    return ["blue", "red"][edge.Name=="A"]
tree_fig("Highlight edge A", UnrootedDendrogram(t), edge_color_callback=callback)

doc.build(story)

#top = ComparisonDisplay(d, d).shape(300)
