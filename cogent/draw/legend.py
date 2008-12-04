#!/usr/bin/env python
from cogent.core import moltype, annotation

from cogent.draw.linear import Display, DisplayPolicy

from reportlab.lib import colors
from reportlab.platypus import Paragraph, SimpleDocTemplate, Table
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class Legend(object):
    """A class for drawing a legend for a display policy
    
    Arguments:
        - policy: a reference to a Display policy class"""
    def __init__(self, policy = DisplayPolicy):
        self.policy = policy
        self.table_style = [
                ('GRID',(0,0),(-1,-1),0.5,colors.grey),
                ('VALIGN',(0,0),(-1,-1),'MIDDLE'),]
    
    def _reformat_table(self, table_data, columns):
        new=[]
        padding_req = (columns) - (len(table_data) % (columns))
        if padding_req != columns:
            table_data += padding_req*[["",""]]
        for i in range(columns-1,len(table_data),columns):
            row = []
            for j in range(i-(columns-1),i+1,1):
                row += table_data[j]
            new.append(row)
        return new
    
    def _makeSampleSequence(self, feature_type):
        seq = moltype.DNA.makeSequence('aaaccggttt' * 7)
        v = seq.addAnnotation(annotation.Feature,
                feature_type, feature_type, [(2,3)])
        v = seq.addAnnotation(annotation.Feature,
                feature_type, feature_type, [(7,18)])
        v = seq.addAnnotation(annotation.Feature,
                feature_type, feature_type, [(20,70)])
        return seq
    
    def _makeSubTable(self, track, columns, column_width):
        diagram_width = column_width * 0.8
        if track.tag is None:
            return None
        temp_table_data = []
        for feature in track:
            if feature == 'blueline' or feature == 'redline':
                continue
            seq = self._makeSampleSequence(feature)
            label = Paragraph(feature, ParagraphStyle('normal'))
            display = Display(seq,
                    policy = self.policy,
                    label_width = 0,
                    min_feature_height = 10,
                    show_scale = False,
                    show_code = False,
                    pad = 0,)
            sample = display.asDrawing(
                    total_width=diagram_width,
                    margin=0,
                    withTrackLabelColumn=False,
                    border=False,),
            temp_table_data.append([label, sample])
        table_data = [
            [Paragraph("<b>"+track.tag+"</b>", ParagraphStyle('normal'))] +
            (2 * columns - 1) * [""]]
        return table_data + self._reformat_table(temp_table_data, columns)
    
    def asDrawing(self, total_width, columns = 3):
        """ Returns the legend as a reportlab.platypus table
        Arguments:
            - total_width: the width of the table in points
            - columns: the number of columns of feature / representation
              pairs
        """
        column_width = total_width / (columns *2)
        table_data = []
        for track in self.policy()._makeTrackDefns():
            sub_table = self._makeSubTable(track, columns, column_width)
            if sub_table:
                table_data += sub_table
        
        return Table(table_data,
                columns*[column_width,column_width],
                style=self.table_style)
    
    def drawToPDF(self, filename, pagesize=A4, *args, **kw):
        """ Writes the Legend to a file
        Arguments:
            - filename: the name of the file
            - pagesize: a tuple of the page dimentions (in points)
              Default is A4
            - columns: the number of columns of feature / representation
              pairs"""
        doc = SimpleDocTemplate(filename,leftMargin=10, rightMargin=10,
                pagesize=pagesize)
        doc.build([self.asDrawing(pagesize[0]*0.8, *args, **kw)])
    

