#!/usr/bin/env python

"""This script/module can do any of 5 different actions with the figures it makes.
When run as a script, one of the actions must be specified:

    - exercise
      The default when used as a module.  Drawing code is used to make a 
      matplotlib figure object but that is all.
    - record
      Save PNG images of the figures in draw_results/baseline/
      To be used before making changes in drawing code.
    - check
      Compare figures with saved baseline figures.  Fail if they don't match.
      Failed figures are saved in draw_results/current.  Also makes an HTML 
      page comparing them with the baseline images.
    - compare
      Save ALL differing figures in draw_results/current and make HTML page 
      comparing them with the baseline images.
    - view
      Save all differing figures in draw_results/current and make HTML page 
      comparing them with the baseline images, along with all the matching
      figures too.
    
"""
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib',
                        lineno=881)

import matplotlib
matplotlib.use('Agg')

import unittest
import sys, os, cStringIO

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


def file_for_test(msg, baseline=False, prefixed=True):
    file_ext = "png"
    dirname = 'baseline' if baseline else 'current'
    if prefixed:
        dirname = os.path.join('draw_results', dirname)
    fname = msg.replace(' ', '_') + '.' + file_ext
    return os.path.join(dirname, fname)

def display2png(display, kw):
    f = cStringIO.StringIO()
    display.makeFigure(**kw).savefig(f, format='png')
    return f.getvalue()

def writefile(fname, content):
    dirname = os.path.dirname(fname)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    f = open(fname, 'wb')
    f.write(content)
    f.close()

def exercise(msg, display, **kw):
    display.makeFigure(**kw)
    
def record(msg, display, **kw):
    png = display2png(display, kw)
    fname = file_for_test(msg, True)
    writefile(fname, png)

class CheckOutput(object):
    def __init__(self, failOnDifference=True, showAll=False):
        if not os.path.exists('draw_results/baseline'):
            raise RuntimeError(
                'No baseline found.  Run "test_draw.py record" first')
        self.results = []
        self.failOnDifference = failOnDifference
        self.showAll = showAll
        self.anyFailures = False
    
    def __call__(self, msg, display, **kw):
        fname = file_for_test(msg, True)
        observed = display2png(display, kw)
        if os.path.exists(fname):
            expected = open(fname, 'rb').read()
            different = observed != expected
            self.results.append((msg, different))
            if different:
                self.anyFailures = True
                writefile(file_for_test(msg, False), observed)
                if self.failOnDifference:
                    raise AssertionError('See draw_results/comparison.html')
                else:
                    print 'difference from', fname
        else:
            raise RuntimeError('No baseline image at %s' % fname)
    
    def writeHTML(self):
        html = ['<html><head><title>Drawing Test Output</title></head>',
                '<body>']
        html.append('<p>%s figures of which %s differ from baseline' % (
                len(self.results), sum(d for (m,d) in self.results)))
        for (msg, different) in self.results:
            fn1 = file_for_test(msg, True, False)
            fn2 = file_for_test(msg, False, False)
            if different:
                html.append('<h2>%s</h2>' % msg)
                html.append('<h3>Old</h3>')
                html.append('<img src="%s"/>' % fn1)
                html.append('<h3>New</h3>')
                html.append('<img src="%s"/>' % fn2)
            elif self.showAll:
                html.append('<h2>%s</h2>' % msg)
                html.append('<img src="%s"/>' % fn1)
            else:
                html.append('<p>%s</p>' % msg)                
            html.append('<hr/>')
        html.append('</body></html>')
        html = '\n'.join(html)
        f = open('draw_results/comparison.html', 'w')
        f.write(html)
        f.close()

    def report(self):
        self.writeHTML()
        if self.anyFailures or self.showAll:
            if sys.platform == 'darwin':
                import subprocess
                subprocess.call(['open', 'draw_results/comparison.html'])
            else:
                print "See draw_results/comparison.html"


def makeSampleSequence():
    seq = 'tgccnwsrygagcgtgttaaacaatggccaactctctaccttcctatgttaaacaagtgagatcgcaggcgcgccaaggc'
    seq = DNA.makeSequence(seq)
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
    seq2 = seq2[:30] + seq2[50:]
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

class DrawingTests(unittest.TestCase):

    def test_seqs(self):
        seqd = Display(seq)
        fig('sequence wrapped at 50', 
            seqd, rowlen=50)
        small = FontProperties(size=7, stretch='extra-condensed')
        fig('squashed sequence', 
            seqd.copy(seq_font=small, colour_sequences=True))
        fig('seq display slice from 5 to 45 starts %s' % seq[5:8], 
            seqd[5:45])

    def test_alns(self):
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

    def test_legend(self):
        from cogent.draw.legend import Legend
        fig('Feature Legend', Legend())

    def test_dotplot(self):
        from cogent.draw.dotplot import Display2D
        fig('2d', Display2D(seq, seq[:40], show_text=False, draw_bases=False))

    def test_trees(self):
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


                
if __name__ != "__main__":
    fig = exercise
else:
    myargs = []
    for arg in ['exercise', 'record', 'check', 'compare', 'view']:
        if arg in sys.argv:
            sys.argv.remove(arg)
            myargs.append(arg)
    if len(myargs) != 1:
        print 'Need one action, got', myargs
        print __doc__
        sys.exit(1)
    action = myargs[0]
    if action == 'record':
        fig = record
    elif action == 'check':
        fig = CheckOutput(True)
    elif action == 'compare':
        fig = CheckOutput(False)
    elif action == 'view':
        fig = CheckOutput(False, True)
    elif action == 'exercise':
        fig = exercise        
    else:
        raise RuntimeError('Unknown action %s' % action)
        
    try:
        unittest.main()
    finally:
        if hasattr(fig, 'report'):
            fig.report()
       
    
