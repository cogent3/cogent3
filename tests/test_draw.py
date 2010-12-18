#!/usr/bin/env python

from cogent import DNA
from cogent.core import alignment, alphabet, annotation
from cogent.draw.linear import *

from cogent import LoadTree
from cogent.draw.dendrogram import *

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
                    "Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

interactive = True
if not interactive:
    import matplotlib
    matplotlib.use('PDF') # or Agg
    file_ext = "pdf" # or png

def makeSampleSequence(mut=False):
    repeat = 'aaaccggtwt'
    if mut: repeat = repeat[:-1]
    seq = DNA.makeSequence(repeat * 10)
    v = seq.addAnnotation(annotation.Feature, 'exon', 'exon', [(20,35)])
    v = seq.addAnnotation(annotation.Feature, 'repeat_unit', 'repeat_unit', [(39,49)])
    v = seq.addAnnotation(annotation.Feature, 'repeat_unit', 'rep2', [(49,60)])
    return seq

def makeSampleAlignment():
    # must be an esier way to make an alignment of annotated sequences!
    from cogent.align.align import global_pairwise, make_dna_scoring_dict
    DNA = make_dna_scoring_dict(10, -8, -8)
    seq1 = makeSampleSequence()[:-2]
    seq2 = makeSampleSequence(mut=True)[2:]
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

def fig(msg, seq_display, **kw):
    if interactive:
        seq_display.showFigure(title=msg, top=.5, **kw)
    else:
        fname = 'draw_test_%s.%s' % (msg.replace(' ', '_'), file_ext)
        seq_display.makeFigure(**kw).savefig(fname)

def green_cg(seq):
    seq = str(seq)
    posn = 0
    result = []
    while True:
        last = posn
        posn = seq.find('CG', posn)
        if posn < 0: break
        result.append('k' * (posn-last)+'gg')
        posn += 2
    result.append('k' * (len(seq)-last))
    return list(''.join(result))

def test_seqs():
    seqd = Display(seq)
    fig('sequence wrapped at 50', 
        seqd, rowlen=50)
    small = FontProperties(size=7, stretch='extra-condensed')
    fig('squashed sequence', 
        seqd.copy(seq_font=small, colour_sequences=True))
    fig('seq display slice from 5 to 45 starts "GGT"', 
        seqd[5:45])

def test_alns():
    alignd = Display(align, colour_sequences=True, min_feature_height=10)
    fig('coloured text alignment', 
        alignd)
    fig('coloured alignment no text', 
        alignd.copy(show_text=False))
    fig('no text and no colour', 
        alignd.copy(show_text=False, colour_sequences=False))
    fig('no shapes', 
        alignd.copy(show_text=False, draw_bases=False))            
    fig('no text or colour or shapes', 
        alignd.copy(show_text=False, colour_sequences=False, draw_bases=False))
    fig('green seqs', 
        alignd.copy(seq_color_callback=green_cg))

# LEGEND
def test_legend():
    from cogent.draw.legend import Legend
    fig('Feature Legend', Legend())

# DOTPLOT

def test_dotplot():
    from cogent.draw.dotplot import Display2D
    fig('2d', Display2D(seq, seq[:40], show_text=False, draw_bases=False))

#fig('reversed', Display(seq[50:10]), 500)
# no good because seqs slice like lists: ie len==0
# complement() probably doesn't know about annotations either.

# TREES

def test_trees():
    treestring = "((A:.1,B:.22)ab:.3,((C:.4,D:.5)cd:.55,E:.6)cde:.7,F:.2)"
    for edge in 'ABCDEF':
        treestring = treestring.replace(edge, edge+edge.lower()*10)
    t = LoadTree(treestring=treestring)
    for klass in [
            UnrootedDendrogram, 
            SquareDendrogram, 
            ContemporaneousDendrogram, 
            ShelvedDendrogram, 
    #        StraightDendrogram, 
    #        ContemporaneousStraightDendrogram
            ]:
        dendro = klass(t)
        dendro.getConnectingNode('Ccccccccccc', 'Eeeeeeeeeee').setCollapsed(
            color="green", label="C, D and E")
        fig(klass.__name__, dendro, shade_param="length", 
            show_params=["length"])
 
    def callback(edge):
        return ["blue", "red"][edge.Name.startswith("A")]
    fig("Highlight edge A", UnrootedDendrogram(t), edge_color_callback=callback)
    
if __name__ == '__main__':
    test_seqs()
    test_alns()
    test_legend()
    test_dotplot()
    test_trees()
    
