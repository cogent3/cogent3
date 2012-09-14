#!/usr/bin/env python
"""Format classes for MotifResults objects."""

from __future__ import division
from matplotlib import use
use('Agg')  #suppress graphical rendering

from cogent.motif.util import MotifFormatter
from cogent.format.pdb_color import get_matching_chains,\
    align_subject_to_pdb, PYMOL_FUNCTION_STRING, MAIN_FUNCTION_STRING
from cogent.format.rna_struct import color_on_structure, draw_structure
from cogent.format.fasta import fasta_from_alignment
from cogent.core.moltype import PROTEIN, RNA, DNA
from cogent.core.alignment import Alignment,SequenceCollection
from cogent.util.dict2d import Dict2D
from numpy import zeros, nonzero
from cogent.align.weights.util import AlnToProfile
from zipfile import ZipFile
from cogent.app.util import get_tmp_filename
from gzip import GzipFile
from pylab import savefig,clf

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Micah Hamady"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Prototype"


def _format_number(number):
    if number is None:
        return "None"
    return "%.2e" % number

def avg(l):
    """Returns the average of a list of numbers."""

    if not l:
        return None
    return sum(l)/len(l)

class MotifStatsBySequence(MotifFormatter):
    """Generates HTML table with Motifs organized by sequence.

        - Each sequence is listed as follows:

        Seq-ID  Combined-P-Value    Motif-ID    Motif-P-Value   Motif-Seq
    """

    def __init__(self,MotifResults=None):
        """init function for MotifStatsBySequence class"""
        MotifFormatter.__init__(self)
        self.MotifResults = MotifResults
        self.ColorMap = self.getColorMap(MotifResults)
        
    def __call__(self,order=None, wrap_html=False, title_class="ntitle",
                 normal_class="normal", cons_thresh=.9):
        """Call method for MotifStatsBySequence class.

        wrap_html: if True, wrap in html + body, else just return table
        table_class: css class to use to format table
        title_class: css class to use to format titles
        normal_class: css class to use to format normal text
        cons_thresh: conservation threshold
        """
        self.ConservationThresh = cons_thresh

        html_list = []
        #if MotifResults is not None
        if self.MotifResults:
            #Find out if the alignment has a combined P value from results
            if 'CombinedP' in self.MotifResults.Results:
                combined_p_string = 'Combined P-Value'
                self.combinedP = True
            else:
                combined_p_string = '&nbsp;'
                self.combinedP = False
            
            #Start HTML string with table and table headers
            html_list = ["""<table cellpadding=2 cellspacing=2 border=0>
                             <tr class="%s">
                                <td>Sequence ID</td>
                                <td>%s</td>
                                <td>Motif ID</td>
                                <td>Motif P-Value</td>
                                <td>Motif Sequence</td>
                             </tr>""" % ( title_class, combined_p_string)]
 
            #For each sequence in alignment get HTML for that sequence
            if not order:
                order = sorted(self.MotifResults.Alignment.keys())

            for seqID in order:
                html_list.append(self.seqLines(seqID, title_class, normal_class))

            html_list.append("</table>")

            if wrap_html:
                return """<html><head><title>Motif Finder Results</title></head><body>%s</table></body></html>""" % ''.join(html_list)
            return  ''.join(html_list) 
            
        return "" 

    def _get_location_dict(self):
        """Builds dict of all locations.
            {module:{seqID:[indices]}}
        """
        location_dict = {}  #Dict with locations of every motif keyed by module
        #Build dict of all the locations:
        # {module:{seqID:[indices]}}
        if self.MotifResults:
            for motif in self.MotifResults.Motifs:
                for module in motif.Modules:
                    location_dict[module]=module.LocationDict
        return location_dict

    Locations = property(_get_location_dict)

    def seqLines(self, seqID, title_class="ntitle", normal_class="normal"):
        """Returns HTML string for single sequence in alignment.

         - Must call for each sequence in the alignment.
         normal_class: css class to use to format rows
         title_class: css class to use to format title cells
        """
        #Variable which signifies if given sequence in alignment contains motifs
        contains_motifs=False
        #Generate first row in table
        html_list = ["""<tr bgcolor="eeeeee" class="%s"><td class="%s">%s</td><td colspan=4>&nbsp;</td></tr>"""%(normal_class, title_class, seqID)]
        #If there is a combined P for the sequences, put it in first row
        if self.combinedP:
            try:
                html_list = ["""<tr bgcolor="#eeeeee" class="%s">
                                <td class="%s">%s</td>
                                <td>%s</td>
                                <td>&nbsp;</td>
                                <td>&nbsp;</td>
                                <td>&nbsp;</td>
                                </tr>""" % ( normal_class, title_class, seqID,
                                         _format_number(float(
                        self.MotifResults.Results['CombinedP'][seqID])))]
            except KeyError:
                pass

        #For each module
        for module in self.Locations.keys():

            cons_seq, cons_con_seq =  self._make_conservation_consensus(module)

            #Check to see if it appeared in the sequence
            if seqID in self.Locations[module]:
                contains_motifs=True
                for index in self.Locations[module][seqID]:
                    cur_seq =  str(module.NamedSeqs[(seqID,index)])
                    html_list.append("""<tr class="%s">
            <td colspan=2>&nbsp;</td>
            <td>%s</td>
            <td>%s</td>
            <td><span style="%s">%s</span><br>%s</td>
            </tr>""" %(normal_class, module.ID,
                       _format_number(module.Pvalue),
                       self.ColorMap[module.ID],
                       cur_seq,
                       """<font face="Courier New, Courier, monospace">%s</font>""" % ''.join( 
                       self._flag_conserved_consensus(cons_con_seq, cons_seq, cur_seq))
                       ))

        if not contains_motifs:
            html_list=[]
        return ''.join(html_list)
                    

class MotifLocationsBySequence(MotifFormatter):
    """Generates HTML table with Motifs organized by sequence.

        - Each sequence is listed as follows:

        Seq-ID  #_bases-module_sequence-#_bases-module_sequence-#_bases-etc 
    """

    def __init__(self,MotifResults=None):
        """init function for MotifLocationsBySequence class"""
        MotifFormatter.__init__(self)
        self.MotifResults = MotifResults
        self.ColorMap = self.getColorMap(MotifResults)

    def _get_location_dict(self):
        """Builds dict of all locations.
            {seqID:{index:module}}
        """
        #Dict with locations of every module keyed by seqID
        locations_list = []
        module_map = {}
        #If MotifResults object exists
        if self.MotifResults:
            #Build dict of all the locations:
            # {seqID:{index:module}}
            #For each motif in MotifResults object
            for motif in self.MotifResults.Motifs:
                #For each module in the Motif
                for module in motif.Modules:
                    #Get the location dict for that module
                    location_dict = module.LocationDict
                    #For each sequence the module is in
                    for seqID in location_dict.keys():
                        #For each module instance in the sequence
                        for index in location_dict[seqID]:
                            #Add module to dict
                            locations_list.append((seqID,index,module))
                            module_map[(seqID,index)] = module
        locations = Dict2D()
        locations.fromIndices(locations_list)
        self.ModuleMap = module_map
        return locations

    Locations = property(_get_location_dict)

    def formatLine(self, seqID, max_index_len):
        """Returns motif line """
        #Variable which signifies if given sequence in alignment contains motifs
        contains_motifs=False
        #List of strings for sequence line
        seq_line_list = []
        #Get all indices for the sequence
        indices = []
        if seqID in self.Locations:
            indices = self.Locations[seqID].keys()
            contains_motifs=True
        indices.sort()
        #Current position in sequence
        pos=0
        #For each index in sequence
        fmt_str = "%0" + max_index_len + "d"
        for index in indices:
            #Add distance between motifs or ends of sequence to list if > 0

            seq_line_list.append(fmt_str % (index-pos))
            #Add module instance sequence to list

            mod_id = self.ModuleMap[(seqID, index)].ID
            seq_line_list.append( """<span style="%s" onmouseover="return overlib('Motif ID: %s');" onmouseout="return nd();">%s</span>""" % \
            (self.ColorMap[mod_id], mod_id,\
            self.Locations[seqID][index].NamedSeqs[(seqID,index)].Sequence))
            #Find new position in sequence
            pos = \
            self.Locations[seqID][index].NamedSeqs[(seqID,index)].Location.End
        #Add distance from end of last module to end of sequence to list if > 0
        if len(self.MotifResults.Alignment.NamedSeqs[seqID])-pos > 0:
            seq_line_list.append(\
            str(len(self.MotifResults.Alignment.NamedSeqs[seqID])\
                                     -pos))

        return '-'.join(seq_line_list), contains_motifs
       
    def seqLines(self, seqID, max_index_len, title_class="ntitle", normal_class="normal"):
        """Returns HTML string for single sequence in alignment.

         - Must call for each sequence in the alignment.

         title_class: css class to use to format title
         normal_class: css class to use to format text
        """
        seq_line_list, contains_motifs = self.formatLine(seqID, max_index_len)
        if contains_motifs:
            #Return sequence string
            return """<tr class="%s">
                        <td class="%s">%s</td>
                        <td>%s</td></tr>"""%( normal_class, title_class,
                            seqID, seq_line_list)
        else:
            return ''
        

    def __call__(self, order=None, wrap_html=False, title_class="ntitle",
                normal_class="normal"):
        """Call method for MotifLocationsBySequence class.
        
            - must pass in an alignment order
        """
        html_list = []

        # need to calculate this acrosss all

        all_indicies = []
        for locs in  self.Locations.values():
            all_indicies.extend(locs.keys())

        max_index_len = str(len(str(max(all_indicies))))

        #If MotifResults is not None
        if self.MotifResults:
            #For each sequence in alignment get HTML for that sequence
            html_list.append("""<tr class="%s"><td>Sequence ID</td><td>Motif Locations</td></tr>""" % title_class)
            if not order:
                order=sorted(self.Locations.keys())
            for seqID in order:
                html_list.append(self.seqLines(seqID, max_index_len, title_class, normal_class))
           
            out_str = "<table cellpadding=2 cellspacing=2 border=0>%s</table>" % ''.join(html_list)
             
            if wrap_html: 
                return """<html><head><title>Motif Finder Results</title></head><body>%s</body></html>""" % out_str 
            return out_str

        return "" 

        
class SequenceByMotif(MotifFormatter):
    """Generates HTML table with sequences organized by Motif.  """

    def __init__(self,MotifResults=None):
        """init function for MotifLocationsBySequence class"""
        MotifFormatter.__init__(self)
        self.MotifResults = MotifResults
        self.ColorMap = self.getColorMap(MotifResults)

    def _get_location_dict(self):
        """Build dict of all the locations:
            {module:{seqID:[indices]}}
        """
        location_dict = {}  #Dict with locations of every motif keyed by module
        if self.MotifResults:
            for motif in self.MotifResults.Motifs:
                for module in motif.Modules:
                    location_dict[module]=module.LocationDict
        return location_dict

    Locations = property(_get_location_dict)

    def __call__(self, wrap_html=False, title_class="ntitle", 
                 normal_class="normal", cons_thresh=.9):
        """Call method for SequenceByMotif class.  """

        #Start HTML string with table and table headers
        html_list = []

        #Get modules
        modules = self.Locations.keys()
        modules.sort()
        self.ConservationThresh = cons_thresh

        #For each module
        for module in modules:
            html_list.append(self.moduleLines(module, title_class, normal_class))

        out_str = """<table cellpadding=2 cellspacing=2 border=0>
               <tr class="%s"><td>Motif ID</td><td>Combined P-Value</td>
                    <td>Sequence ID</td>
                    <td>Motif Sequence</td>
               </tr>
               %s
               </table>""" % (title_class, ''.join(html_list))


        if html_list:
            if wrap_html:
                return"""<html><head><title>Motif Finder Results</title></head><body>%s</body></html>""" % out_str
            else:
                return out_str
        return ""

    def _highlightConsensus(self, con_seq, cons_con_seq, cur_seq, cur_color):
        """
        Hightlight positions identical to consensus
        """
        grey_style = """background-color: #dddddd; font-family: 'Courier New', Courier"""

        span_fmt = """<span style="%s">%s</span>"""
        h_str = []
        for ix in range(len(cur_seq)):
            cur_c = cur_seq[ix]
            if cur_c == cons_con_seq[ix]:
                h_str.append(span_fmt % (cur_color,cur_c))
            elif cur_c == con_seq[ix]:
                h_str.append(span_fmt % (grey_style,cur_c))
            else:
                h_str.append(cur_c)
        return ''.join(h_str)

    def moduleLines(self, module, title_class="ntitle", normal_class="normal"):
        """Returns HTML string for single module.

         - Must call for each module.
        """
        cons_seq, cons_con_seq = self._make_conservation_consensus(module)

        cur_color =  self.ColorMap[module.ID]

        #Generate first row in table
        html_list = ['<tr bgcolor="#eeeeee" class="%s"><td class="%s">%s</td><td>%s</td><td>&nbsp;</td><td><span style="%s">%s</span></td></tr>'%\
                     (normal_class, title_class, 
                      module.ID,
                      _format_number(module.Pvalue),
                      cur_color, 
                      cons_seq, 
                      )]
        sequences = self.Locations[module].keys()
        sequences.sort()
        #For each sequence
        for seq in sequences:
            cur_seq = str(module.NamedSeqs[(seq,self.Locations[module][seq][0])] )
            
            html_list.append("""<tr class="%s">
                                    <td>&nbsp;</td>
                                    <td>&nbsp;</td>
                                    <td class="%s">%s</td>
                                    <td style="font-family: 'Courier New', Courier, monospace">%s</td>
                               </tr>""" % (normal_class, title_class,
               seq, 
               #_format_number(module.Pvalue), 
               self._highlightConsensus(cons_seq, cons_con_seq, cur_seq, cur_color)
               ))
        return ''.join(html_list)


class HighlightMotifs(MotifFormatter):
    """Generates HTML table with sequences highlighted """

    def makeModuleMap(self, motif_results):
        """
        Need to extract this b/c can't pickle motif_results... grr.

        motif_results: MotifResults object
        keep_module_ids: list of module ids to keep
        """
        module_map = {}  #Dict with locations of every motif keyed by module
        if motif_results:
            for motif in motif_results.Motifs:
                for module in motif.Modules:
                    mod_len = len(module)
                    mod_id = str(module.ID)
                    for skey, indexes in module.LocationDict.items():
                        if skey not in module_map:
                            module_map[skey] = []
                        for ix in indexes:
                            module_map[skey].append((ix, mod_id, mod_len))
        return module_map



    def __init__(self, MotifResults, NodeOrder=None, KeepIds=None,\
        KeepAll=False, MolType=PROTEIN):
        """Set up color map and motif results

        ModuleMap: flattened map (b/c of pickle problem.)
                generate using make_module_map() function
        Alignment: SequenceCollection or Alignment object
        KeepIds: list of module ids to keep
        KeepAll: When True, ignores KeepIds and highlights all motifs
        """
        MotifFormatter.__init__(self)
        ModuleMap = self.makeModuleMap(MotifResults)
        module_ids = set([])
        for skey, slist in ModuleMap.items():
            for stup in slist:
                module_ids.add(stup[1])
        self.ColorMap = self.getColorMapS0(sorted(list(module_ids)))
        self.ModuleMap = ModuleMap 
        self.Alignment = MotifResults.Alignment
        if KeepIds is None:
            KeepIds = []
        self.KeepIds = set(KeepIds)
        self.KeepAll = KeepAll
        if not NodeOrder:
            NodeOrder=self.Alignment.Names
        self.NodeOrder = NodeOrder
        self.MolType = MolType
        self.GapMap = self.getGapMap()
        self.HighlightMap = {}

    def __call__(self, title_class="ntitle", normal_class="normal", row_class="highlight", table_style=''):
        """Call method for HightlightMotifs class.  """

        #Start HTML string with table and table headers
        html_list = []

        #For each module
        for seq_id in self.NodeOrder:
            html_list.append(self.highlightSeq(seq_id, title_class, row_class))

        out_str = """<table cellpadding=2 cellspacing=2 border=0 %s>
               <tr class="ntitle">
                    <td colspan=2><p>Selected Motifs Highlighted on Sequences:</p> </td>
               </tr>
         
               <tr class="%s"><td nowrap>Sequence ID</td>
                    <td>Sequence</td>
               </tr>
               %s
               </table>""" % (table_style, title_class, ''.join(html_list))


        if html_list:
                return out_str
        return ""
    
    def getGapMap(self):
        """Returns dict mapping gapped_coord to ungapped_coord in self.Alignment
        
            - {seq_id:{gapped_coord:ungapped_coord}}
        """
        gap_map = {}
        for k,v in self.Alignment.items():
            gapped, ungapped = self.MolType.gapMaps(v)
            gap_map[k] = gapped
        return gap_map

    def highlightSeq(self, seq_id, title_class="ntitle", row_class="highlight"):
        """Returns HTML string for single sequence.

        seq_id: sequence_id to highlight 
        """
        #Generate first row in table
        row_tmpl = """<tr class="%s">
                          <td class="%s">%s</td>
                          <td>%s</td>
                      </tr>"""
  
        mo_span_tmpl = """<span style="%s" onmouseover="return overlib('%s');" onmouseout="return nd();">%s</span>"""
        seq_list = list(self.Alignment.NamedSeqs[seq_id])
        seq_len = len(seq_list)
        seq_mask = zeros(seq_len)
        mod_id_map = {}

        if seq_id in self.ModuleMap:
            for mod_tup in self.ModuleMap[seq_id]:
                ix, mod_id, mod_len = mod_tup
                
                # skip modules we con't care about
                if not self.KeepAll and mod_id not in self.KeepIds:
                    continue

                mod_mask = zeros(seq_len)

                # mask motif region
                for i in range(ix,ix+mod_len):
                    gapped_ix = self.GapMap[seq_id][i]
                    mod_mask[gapped_ix] = 1

                # add to sequence map
                seq_mask += mod_mask

                # map module ids to indexes
                for jx in range(ix,ix+mod_len):
                    gapped_jx = self.GapMap[seq_id][jx]
                    if gapped_jx not in mod_id_map:
                        mod_id_map[gapped_jx] = [] 
                    mod_id_map[gapped_jx].append(mod_id)

            # get module regions
            
            #need to take [0] element of nonzero() since numpy returns a tuple
            # where Numeric did not
            for jx in nonzero(seq_mask)[0]:

                # if overlapping use red background, otherwise display color
                if seq_mask[jx] > 1:
                    style = "background-color: red"
                else:
                    #style = self.ColorMap[mod_id_map[jx][0]].replace("font-family: 'Courier New', Courier, monospace", "")
                    style = self.ColorMap[mod_id_map[jx][0]]

                seq_list[jx] = mo_span_tmpl % (style,
                                                "<br>".join(['Motif ID: %s' % \
                                                x for x in mod_id_map[jx]]),
                                                seq_list[jx])

        # cache data
        self.HighlightMap[seq_id] = ''.join(seq_list)
        # return row output 
        return row_tmpl % (row_class, title_class, seq_id, ''.join(seq_list))



class HighlightMotifsForm(MotifFormatter):
    """Generates HTML form to submit module ids """

    def __init__(self,
                 MotifResults,
                 FormAction="/cgi-bin/motifcluster/highlight.py",
                 FormTarget="_blank"):
        """Set up color map and motif results

        MotifResults: MotifResults object
        ModuleIds: List of module ids to highlight 
        FormAction: Form action cgi-script
        FormTarget: Form target 
        """
        MotifFormatter.__init__(self)
        self.MotifResults = MotifResults
        self.ColorMap = self.getColorMap(MotifResults)
        self.Modules= self._get_modules()
        self.FormAction = FormAction
        self.FormTarget = FormTarget

    def _get_modules(self):
        """Build map of modules
            {seq_id:[(ix, module_id, module_len)]
        """
        modules = [] 
        if self.MotifResults:
            for motif in self.MotifResults.Motifs:
                for module in motif.Modules:
                    modules.append(module)
        return modules


    def __call__(self, title_class="ntitle", normal_class="normal", highlight_class="highlight"):
        """Call method for HightlightMotifs class.  """

        #Start HTML string with table and table headers
        cells = []

        # format cells
        for module in self.Modules:
            cur_tup = (module.Pvalue, len(module.LocationDict), 
                        self.moduleRow(module))
            cells.append(cur_tup)

        # sort by p value, then frequency
        cells = [x[-1] for x in sorted(cells)]
       
        cur_cells = []
        header_tmpl = """<td bgcolor=ffffff>&nbsp;</td><td>ID</td> <td>Motif</td> <td>Fequency</td> <td>P-Value</td>"""

        header_cells = []
        num_headers = len(cells)
        if num_headers > 3:
            num_headers = 3
        for i in range(num_headers):
            header_cells.append(header_tmpl)
        header_row = """<tr bgcolor=eeeeee class="%s">%s</tr>""" % (title_class, ''.join(header_cells))

        html_out = []
        html_out.append("""<tr class="%s">""" % highlight_class)

        for ix, cell in enumerate(cells):
            if ix % 3 == 0:
                if cur_cells:
                    html_out.append("""%s</tr><tr class="%s">""" % (''.join(cur_cells), highlight_class))
                cur_cells = []
                cur_cells.append(cell)
            else:
                cur_cells.append(cell)

        html_out.append("%s</tr>" % ''.join(cur_cells))
        out_str = """
               <form action="%s" method="POST" target="%s">
               <table cellpadding=2 cellspacing=2 border=0>
               %s
               %s
               </table>
               <p>
               <input type="Submit" value="Highlight Selected Motifs" />
               </p>
               </form> 
               """ % (self.FormAction, self.FormTarget,
                      header_row, ''.join(html_out))

        if cells:
            return out_str
        return "" 

    def moduleRow(self, module):
        """Returns HTML string for single module.

        module: module to generate 
        """
        #Generate first row in table
        cells_tmpl = """<td bgcolor=dddddd><input type="checkbox" name="module_ids" value="%s" checked /></td><td>%s</td><td ><span style="%s">%s</span></td><td>%d</td><td>%s</td>"""
  
        # return row output 
        return cells_tmpl % (module.ID, 
                             module.ID, 
                             self.ColorMap[module.ID], 
                             str(module),
                             #module.ConsensusSequence,
                             len(module.LocationDict),
                             _format_number(module.Pvalue))

   
class HighlightOnCrystal(MotifFormatter):
    """Generates pymol script to highlight motifs on crystal structure.
    """
    
    def __init__(self, MotifResults=None, cons_thresh=0.9, MolType=PROTEIN):
        """init function for HighlightOnAlignment class.
       
        MotifResults: motif results object
        """
        MotifFormatter.__init__(self)
        self.ConservationThresh = float(cons_thresh)
        ModuleMap, ModuleConsMap = self.makeModuleMap(MotifResults)
        module_ids = set([])
        for skey, slist in ModuleMap.items():
            for stup in slist:
                module_ids.add(stup[1])
        self.ColorMapHex = self.getColorMapS0(sorted(list(module_ids)))
        self.ModuleMap = ModuleMap 
        self.ModuleConsMap = ModuleConsMap 
        self.MotifResults = MotifResults
        self.ColorMap = self.getColorMapRgb(MotifResults)
        self.GapMap = {}
        self.HighlightMap = {}
        self.MolType=MolType

        self.RunScriptString = \
'''
from pymol import cmd
cmd.load("%s")
cmd.do("run %s")
'''
        self.ColorFunctionString = \
'''
color_map = %s
color_command_list = %s
sticks_command_list = %s
#Set color list using color_map
set_color_list(list(color_map.items()))
#Set seq colors
for color_cmd in color_command_list:
    colors,indices,chain_id = color_cmd
    set_seq_colors(colors,indices,chain_id)
#Set sticks
for sticks_cmd in sticks_command_list:
    indices,chain_id = sticks_cmd
    set_show_shapes(indices,chain_id)

'''
    
    def __call__(self,seq_id,pdb_id,\
                sequence_type='Protein',\
                zipfile_dir='.',
                pdb_dir='/quicksand2/hamady/data/cron_sync/pdb/'):
        """call method for HighlightOnCrystal class.
        Generates pymol script for highlighting motifs on crystal structure and
        creates .zip archive with pdb file and pymol script.
        """
        #Get PDB file
        curr_pdb = \
            [x.rstrip("\n") for x in self.getPdb(pdb_id,pdb_dir).readlines()]
        

        #Get subject sequence
        subject_seq = self.MotifResults.Alignment.NamedSeqs[seq_id]
        #Get PDB chains
        pdb_matching, ungapped_to_pdb = \
            get_matching_chains(subject_seq,curr_pdb,sequence_type)
        
        pdb_aligned = align_subject_to_pdb(subject_seq,pdb_matching)
        #get color command list
        color_command_list, found_seq_motifs, missed_seq_motifs = \
            self.makeColorCommandLists(seq_id,pdb_aligned,ungapped_to_pdb)
        #get sticks command list
        sticks_command_list = \
            self.makeSticksCommandsConservedPositions(seq_id, pdb_aligned,\
                                                    ungapped_to_pdb)
        
        #Generate pdb file
        pdb_out = pdb_id+'.pdb'
        
        #Generate pymol script
        pymol_script_list = [PYMOL_FUNCTION_STRING,MAIN_FUNCTION_STRING]

        pymol_script_list.append(self.ColorFunctionString % (self.ColorMap,\
                                color_command_list,sticks_command_list))
        pymol_script_string = ''.join(pymol_script_list)
        pymol_script_name = '%s_motif_coloring.pml' % (pdb_id)
        
        pymol_execute_name = '%s_double_click_me.pml' % (pdb_id)
        pymol_execute_string = \
            self.RunScriptString % (pdb_out,pymol_script_name)
        
        
        
        #Generate zip file

        output_pre = get_tmp_filename(zipfile_dir, prefix="pdb_%s_" % pdb_id) 
        if output_pre.endswith(".txt"):
            output_pre = output_pre[:-4]
        zip_dir = output_pre.split("/")[-1]
        output_filename = output_pre + ".zip"
        web_name = output_filename.split("/")[-1]

        curr_zip = ZipFile(output_filename,'w')
        curr_zip.writestr(zip_dir + "/" + pymol_script_name,pymol_script_string)

        curr_zip.writestr(zip_dir + "/" + pdb_out,'\n'.join(curr_pdb))
        
        curr_zip.writestr(zip_dir + "/" + pymol_execute_name,pymol_execute_string)
        curr_zip.close()
        
        alignment_html = {}
        for k,v in pdb_aligned.items():
            alignment_html[k]=self.highlightSeq(seq_id,v[0],pdb_id,v[1])

        #set up return dictionary
        return_dir = { "output_filename":output_filename, 
                       "web_name":web_name,
                       "colored_alignment":alignment_html,
                       #"found_seq_motifs":found_seq_motifs,
                       "found_seq_motifs":[(module.ID, _format_number(module.Pvalue)) for module in self.MotifResults.Modules if module.ID in found_seq_motifs],
                       "missed_seq_motifs":missed_seq_motifs,
                       "all_motifs":[(module.ID, _format_number(module.Pvalue)) for module in self.MotifResults.Modules],
                       "all_motif_colors":self.ColorMapHex}
        
        return return_dir

    def getConservedPositions(self,min_conservation=1.0):
        """Returns dict mapping motif id to list of conserved positions.
        """
        conserved_positions = {}
        for motif in self.MotifResults.Motifs:
            for module in motif.Modules:
                curr_id = module.ID
                conserved_positions[curr_id]=[]
                curr_profile = AlnToProfile(module,self.MotifResults.MolType)
                for ix,pos in enumerate(curr_profile.rowMax()):
                    if pos >= min_conservation:
                        conserved_positions[curr_id].append(ix)
        return conserved_positions

    def makeColorCommandLists(self,seq_id, pdb_aligned, ungapped_to_pdb):
        """Returns lists of (colors, indices, and chain_id) for coloring.
        
            - each chain is a separate tuple in the list.
            - colors are named by motif id
        """
        #list of motifs found in pdb sequence
        found_seq_motifs = []
        #list of motifs not in pdb sequence
        missed_seq_motifs = []
        
        color_command_list = []
        #Get locations by sequence
        locations = \
            MotifLocationsBySequence(self.MotifResults).Locations[seq_id]
        for chain, aligned in pdb_aligned.items():
            #Get subject gap map
            subject_gapped,subject_ungapped = \
                self.MotifResults.MolType.gapMaps(aligned[0])
            #Get pdb gap map
            pdb_gapped,pdb_ungapped = \
                self.MotifResults.MolType.gapMaps(aligned[1])
            

            for curr_ix, curr_module in locations.items():
                curr_module_len = len(str(curr_module))
                curr_module_id = curr_module.ID
                #curr_color = self.ColorMap[curr_module_id]
                curr_color = "color_" + str(curr_module_id)
                #get index list
                ix_list = []
                #for each position in motif
                for i in range(curr_ix,curr_ix+curr_module_len):
                    #Get the gapped index of the motif in the subject seq
                    try:
                        sub_gap = subject_gapped[i]
                    except KeyError:
                        continue

                    #Get the ungapped index of the motif in the pdb seq
                    try:
                        pdb_ungap = pdb_ungapped[sub_gap]
                        
                        #Get the index of the position in pdb coordinates

                        pdb_ix = ungapped_to_pdb[chain][pdb_ungap]
                        ix_list.append(pdb_ix)
                    except KeyError:
                        continue

                    
                if ix_list:
                    ix_string = '+'.join(map(str,ix_list))
                    color_command_list.append(([curr_color],[ix_string],chain))
                    found_seq_motifs.append(curr_module_id)
                else:
                    missed_seq_motifs.append(curr_module_id)
        

        return color_command_list, found_seq_motifs, missed_seq_motifs
                
    
    def makeSticksCommandsConservedPositions(self,seq_id, pdb_aligned,\
        ungapped_to_pdb, min_conservation=1.0):
        """Returns list of (indices, chain_id) to show sticks at indices.
        
            - each chain is a separate tuple in the list.
        """
        conserved_positions = self.getConservedPositions()
        show_command_list = []
        #Get locations by sequence
        locations = \
            MotifLocationsBySequence(self.MotifResults).Locations[seq_id]
        for chain, aligned in pdb_aligned.items():
            #Get subject gap map
            subject_gapped,subject_ungapped = \
                self.MotifResults.MolType.gapMaps(aligned[0])
            #Get pdb gap map
            pdb_gapped,pdb_ungapped = \
                self.MotifResults.MolType.gapMaps(aligned[1])
            
            for curr_ix, curr_module in locations.items():
                curr_module_id = curr_module.ID
                curr_module_conserved = conserved_positions[curr_module_id]
                #get index list
                ix_list = []
                #for each position in motif
                for i in curr_module_conserved:
                    #Get the gapped index of the motif in the subject seq
                    try:
                        sub_gap = subject_ungapped[i]
                    except KeyError:
                        continue
                    #Get the ungapped index of the motif in the pdb seq
                    try:
                        pdb_ungap = pdb_gapped[sub_gap]
                        #Get the index of the position in pdb coordinates
                        
                        pdb_ix = ungapped_to_pdb[pdb_ungap]
                        ix_list.append(pdb_ix)
                    except KeyError:
                        continue
                
                if ix_list:
                    show_command_list.append((ix_list,chain))
        
        return show_command_list
    
    
    def getPdb(self,pdb_id, pdb_dir):
        """Returns open pdb file.
        
            - currently gets pdb file from pdb website.
        """
        pdb_file = pdb_dir + "pdb%s.ent" % pdb_id.lower()
        of = None
        try:
            of = open(pdb_file)
        except Exception, e:
            of = GzipFile(pdb_file + ".gz")
        return of 
    
    def getGapMap(self,seq_id,gapped_seq):
        """Returns dict mapping gapped_coord to ungapped_coord in self.Alignment
        
            - {seq_id:{gapped_coord:ungapped_coord}}
        """
        gap_map = {}
        gapped, ungapped = self.MolType.gapMaps(gapped_seq)
        gap_map[seq_id] = gapped
        return gap_map

    def highlightSeq(self, seq_id, seq_aligned, pdb_id, pdb_aligned,\
        title_class="ntitle", row_class="highlight"):
        """Returns HTML string for single sequence.

        seq_id: sequence_id to highlight 
        """
        #Generate first row in table
        seq_row_tmpl = """<tr class="%s">
                          <td class="%s">%s</td>
                          <td>%s</td>
                      </tr>"""
  
        mo_span_tmpl = """<span style="%s" onmouseover="return overlib('%s');" onmouseout="return nd();">%s</span>"""
        seq_list = list(seq_aligned)
        pdb_list = list(pdb_aligned)
        seq_len = len(seq_list)
        seq_mask = zeros(seq_len)
        #pdb sequence mask
        pdb_mask = zeros(seq_len)
        mod_id_map = {}
        
        self.GapMap = self.getGapMap(seq_id,seq_aligned)

        cons_str = list(" " * len(seq_list)) 
        if seq_id in self.ModuleMap:
            for mod_tup in self.ModuleMap[seq_id]:
                ix, mod_id, mod_len = mod_tup
                

                cons_seq, cons_con_seq =  self.ModuleConsMap[mod_id]
                cur_seq =  ''.join(seq_list[ix:ix+mod_len])
                cc_str = self._flag_conserved_consensus(cons_con_seq, cons_seq, cur_seq)
                cons_str[ix:ix+mod_len] = cc_str

                mod_mask = zeros(seq_len)
                #pdb module mask
                pdb_mod_mask = zeros(seq_len)

                # mask motif region
                for i in range(ix,ix+mod_len):
                    gapped_ix = self.GapMap[seq_id][i]
                    mod_mask[gapped_ix] = 1
                    #only allow coloring of ungapped pdb sequence
                    if pdb_list[gapped_ix] != '-':
                        pdb_mod_mask[gapped_ix]=1

                # add to sequence map
                seq_mask += mod_mask
                pdb_mask += pdb_mod_mask

                # map module ids to indexes
                for jx in range(ix,ix+mod_len):
                    gapped_jx = self.GapMap[seq_id][jx]
                    if gapped_jx not in mod_id_map:
                        mod_id_map[gapped_jx] = [] 
                    mod_id_map[gapped_jx].append(mod_id)

            mm_str = []
            for ix in range(len(seq_list)):
                if seq_list[ix] == pdb_list[ix]:
                    mm_str.append("|")
                else:
                    mm_str.append("*")

            # get module regions
            for jx in nonzero(seq_mask)[0]:

                # if overlapping use red background, otherwise display color
                if seq_mask[jx] > 1:
                    style = "background-color: red"
                else:
                    style = self.ColorMapHex[mod_id_map[jx][0]].replace("font-family: 'Courier New', Courier, monospace", "")

                seq_list[jx] = mo_span_tmpl % (style,
                                                "<br>".join(['Motif ID: %s' % \
                                                x for x in mod_id_map[jx]]),
                                                seq_list[jx])
                pdb_list[jx] = mo_span_tmpl % (style,
                                                "<br>".join(['Motif ID: %s' % \
                                                x for x in mod_id_map[jx]]),
                                                pdb_list[jx])
        
        clean_cons_str = []
        for item in cons_str:
            if item == " ":
                clean_cons_str.append("&nbsp;")
            else:
                clean_cons_str.append(item)

        cons_str = ''.join(clean_cons_str)

        # cache data
        self.HighlightMap[seq_id] = ''.join(seq_list)
        # return row output 
        seq_row = seq_row_tmpl % (row_class, title_class, seq_id, ''.join(seq_list))
        pdb_row = seq_row_tmpl % (row_class, title_class, pdb_id, ''.join(pdb_list))
        high_row = seq_row_tmpl % (row_class, title_class, "Cons", ''.join(cons_str))
        mm_row = seq_row_tmpl % (row_class, title_class, "Mismatch", ''.join(mm_str))

        return ''.join(['<table>',high_row,seq_row,mm_row, pdb_row,'</table>'])


    def makeModuleMap(self, motif_results):
        """
        Need to extract this b/c can't pickle motif_results... grr.

        motif_results: MotifResults object
        keep_module_ids: list of module ids to keep
        """
        module_map = {}  #Dict with locations of every motif keyed by module
        module_cons_map = {}
        if motif_results:
            for motif in motif_results.Motifs:
                for module in motif.Modules:
                    mod_id = str(module.ID)
                    mod_len = len(str(module))
                    
                    if mod_id not in module_cons_map:
                        module_cons_map[mod_id] = self._make_conservation_consensus(module)
                    for skey, indexes in module.LocationDict.items():
                        if skey not in module_map:
                            module_map[skey] = []
                        for ix in indexes:
                            module_map[skey].append((ix, mod_id, mod_len))

        return module_map, module_cons_map

class ColorSecondaryStructurePostscript(MotifFormatter):
    """Generates postscript file with motifs highlighted on 2D structure  """

    def makeModuleMap(self, motif_results):
        """
        Need to extract this b/c can't pickle motif_results... grr.

        motif_results: MotifResults object
        keep_module_ids: list of module ids to keep
        """
        module_map = {}  #Dict with locations of every motif keyed by module
        if motif_results:
            for motif in motif_results.Motifs:
                for module in motif.Modules:
                    mod_len = len(module)
                    mod_id = str(module.ID)
                    for skey, indexes in module.LocationDict.items():
                        if skey not in module_map:
                            module_map[skey] = []
                        for ix in indexes:
                            module_map[skey].append((ix, mod_id, mod_len))
        return module_map


    def __init__(self, MotifResults, KeepIds=None,\
        KeepAll=True, MolType=RNA, strict=True,circle_motif_id=None,\
        SkipIds=None):
        """Set up color map and motif results

        ModuleMap: flattened map (b/c of pickle problem.)
                generate using make_module_map() function
        Alignment: SequenceCollection or Alignment object
        KeepIds: list of module ids to keep
        KeepAll: When True, ignores KeepIds and highlights all motifs
        """
        MotifFormatter.__init__(self)
        self.ModuleMap = self.makeModuleMap(MotifResults)
        module_ids = set([])
        for skey, slist in self.ModuleMap.items():
            for stup in slist:
                module_ids.add(stup[1])
        self.ColorMap = self.getColorMapRgb(MotifResults)
        overlap_color = {'overlap_color':(1.0,0.0,0.0)}
        self.ColorMap.update(overlap_color)
        
        self.Alignment = MotifResults.Alignment
        if KeepIds is None:
            KeepIds = []
        self.KeepIds = set(KeepIds)
        self.KeepAll = KeepAll
        self.MolType = MolType
        self.GapMap = self.getGapMap()
        self.HighlightMap = {}
        self.Strict=strict
        self.CircleId = circle_motif_id
        self.SkipIds=SkipIds

    def __call__(self,seq_id,sequence,struct,write_dir='.'):
        """Call method for ColorSecondaryStructure class.  """

        indices, colors = self.getColorIndices(seq_id,sequence)
        circle_indices = []
        if self.CircleId:
            circle_indices = \
                self.getCircleIndices(seq_id,sequence,self.CircleId)
        
        structure_postscript = color_on_structure(\
            sequence=sequence,\
            struct=struct,\
            color_map = self.ColorMap,\
            indices=indices,\
            colors=colors,\
            circle_indices=circle_indices)
        file_id = seq_id.split()[0]
        ps_out_path = \
            write_dir+'/'+file_id+'_secondary_struct.ps'
        ps_out = open(ps_out_path,'w')
        ps_out.write(structure_postscript)
        ps_out.close()
        return ps_out_path
            
    def getGapMap(self):
        """Returns dict mapping gapped_coord to ungapped_coord in self.Alignment
        
            - {seq_id:{gapped_coord:ungapped_coord}}
        """
        gap_map = {}
        for k,v in self.Alignment.items():
            gapped, ungapped = self.MolType.gapMaps(v)
            gap_map[k] = gapped
        return gap_map

    def getColorIndices(self, seq_id, sequence):
        """Returns list of indices and colors for a given sequence.

        seq_id: sequence ID to highlight on structure
        """
        #seq_list = list(self.Alignment.NamedSeqs[seq_id])
        seq_list = list(sequence)
        seq_len = len(seq_list)
        seq_mask = zeros(seq_len)
        mod_id_map = {}
        indices = []
        colors = []
        gapped,ungapped = self.MolType.gapMaps(sequence)
        self.GapMap[seq_id]= gapped
        if self.Strict:
            if seq_id not in self.ModuleMap:
                raise IndexError, 'seq_id %s not in ModuleMap'%(seq_id)
        else:
            if seq_id not in self.ModuleMap:
                return [],[]

        for mod_tup in self.ModuleMap[seq_id]:
            ix, mod_id, mod_len = mod_tup
            
            # skip modules we con't care about
            if not self.KeepAll and mod_id not in self.KeepIds:
                continue
            elif mod_id in self.SkipIds:
                continue

            mod_mask = zeros(seq_len)

            # mask motif region
            for i in range(ix,ix+mod_len):
                gapped_ix = self.GapMap[seq_id][i]
                mod_mask[gapped_ix] = 1
            # add to sequence map
            seq_mask += mod_mask

            # map module ids to indexes
            for jx in range(ix,ix+mod_len):
                gapped_jx = self.GapMap[seq_id][jx]
                if gapped_jx not in mod_id_map:
                    mod_id_map[gapped_jx] = [] 
                mod_id_map[gapped_jx].append(mod_id)


        # get module regions
        #need to take [0] element of nonzero() since numpy returns a tuple
        # where Numeric did not
        for kx in nonzero(seq_mask)[0]:
            # if overlapping use red background, otherwise display color
            if seq_mask[kx] > 1:
                curr_color = 'overlap_color'
            else:
                mod_id = mod_id_map[kx][0]
                curr_color = 'color_'+mod_id
            #append indices.  Must start at 1, not 0 for RNAplot to work
            indices.append(kx+1)
            colors.append(curr_color)
            
        return indices, colors

    def getCircleIndices(self, seq_id, sequence, motif_id):
        """Returns list of indices to be circled for a given sequence.

            seq_id: sequence ID to highlight on structure
            sequence: sequence string
            motif_id: ID of motif to circle.
        """
        seq_list = list(sequence)
        seq_len = len(seq_list)
        seq_mask = zeros(seq_len)
        mod_id_map = {}
        indices = []

        gapped,ungapped = self.MolType.gapMaps(sequence)
        self.GapMap[seq_id]= gapped
        if self.Strict:
            if seq_id not in self.ModuleMap:
                raise IndexError, 'seq_id %s not in ModuleMap'%(seq_id)
        else:
            if seq_id not in self.ModuleMap:
                return [],[]

        for mod_tup in self.ModuleMap[seq_id]:
            ix, mod_id, mod_len = mod_tup
            if mod_id != motif_id:
                continue
            
            # skip modules we con't care about
            if not self.KeepAll and mod_id not in self.KeepIds:
                continue

            mod_mask = zeros(seq_len)

            # mask motif region
            for i in range(ix,ix+mod_len):
                gapped_ix = self.GapMap[seq_id][i]
                mod_mask[gapped_ix] = 1
            # add to sequence map
            seq_mask += mod_mask

            # map module ids to indexes
            for jx in range(ix,ix+mod_len):
                gapped_jx = self.GapMap[seq_id][jx]
                if gapped_jx not in mod_id_map:
                    mod_id_map[gapped_jx] = [] 
                mod_id_map[gapped_jx].append(mod_id)


        # get module regions
        #need to take [0] element of nonzero() since numpy returns a tuple
        # where Numeric did not
        for kx in nonzero(seq_mask)[0]:
            #append indices.  Must start at 1, not 0 for RNAplot to work
            indices.append(kx+1)
            
        return indices

class ColorSecondaryStructureMatplotlib(MotifFormatter):
    """Generates png file with motifs highlighted on 2D structure  """

    def makeModuleMap(self, motif_results):
        """
        Need to extract this b/c can't pickle motif_results... grr.

        motif_results: MotifResults object
        keep_module_ids: list of module ids to keep
        """
        module_map = {}  #Dict with locations of every motif keyed by module
        if motif_results:
            for motif in motif_results.Motifs:
                for module in motif.Modules:
                    mod_len = len(module)
                    mod_id = str(module.ID)
                    for skey, indexes in module.LocationDict.items():
                        if skey not in module_map:
                            module_map[skey] = []
                        for ix in indexes:
                            module_map[skey].append((ix, mod_id, mod_len))
        return module_map


    def __init__(self, MotifResults, KeepIds=None,\
        KeepAll=True, MolType=RNA, strict=True,circle_motif_id=None,\
        SkipIds=None,square_motif_id=None,square_label=None):
        """Set up color map and motif results

        ModuleMap: flattened map (b/c of pickle problem.)
                generate using make_module_map() function
        Alignment: SequenceCollection or Alignment object
        KeepIds: list of module ids to keep
        KeepAll: When True, ignores KeepIds and highlights all motifs
        """
        MotifFormatter.__init__(self)
        self.ModuleMap = self.makeModuleMap(MotifResults)
        module_ids = set([])
        for skey, slist in self.ModuleMap.items():
            for stup in slist:
                module_ids.add(stup[1])
        self.ColorMap = self.getColorMapRgb(MotifResults)
        overlap_color = {'overlap_color':(1.0,0.0,0.0)}
        self.ColorMap.update(overlap_color)
        
        self.Alignment = MotifResults.Alignment
        if KeepIds is None:
            KeepIds = []
        self.KeepIds = set(KeepIds)
        self.KeepAll = KeepAll
        self.MolType = MolType
        self.GapMap = self.getGapMap()
        self.HighlightMap = {}
        self.Strict=strict
        self.CircleId = circle_motif_id
        self.SquareId = square_motif_id
        self.SquareLabel = square_label
        self.SkipIds=SkipIds

    def __call__(self,seq_id,sequence,struct,write_dir='.'):
        """Call method for ColorSecondaryStructure class.  """

        indices, colors = self.getColorIndices(seq_id,sequence)
        circle_indices = []
        if self.CircleId:
            circle_indices = \
                self.getMarkedIndices(seq_id,sequence,self.CircleId)
        square_indices = []
        if self.SquareId:
            square_indices = \
                self.getMarkedIndices(seq_id,sequence,self.SquareId)
        if self.SquareLabel is None:
            self.SquareLabel = ''
        
        draw_structure(sequence, struct, indices=indices, colors=colors,\
            circle_indices=circle_indices, square_indices=square_indices)
        file_id = seq_id.split()[0]
        struct_out_path = '%s/%s_%s_secondary_struct.png'%(write_dir,\
            file_id,self.SquareLabel)
        savefig(struct_out_path,format='png',dpi=150)
        clf()
        return struct_out_path
            
    def getGapMap(self):
        """Returns dict mapping gapped_coord to ungapped_coord in self.Alignment
        
            - {seq_id:{gapped_coord:ungapped_coord}}
        """
        gap_map = {}
        for k,v in self.Alignment.items():
            gapped, ungapped = self.MolType.gapMaps(v)
            gap_map[k] = gapped
        return gap_map

    def getColorIndices(self, seq_id, sequence):
        """Returns list of indices and colors for a given sequence.

        seq_id: sequence ID to highlight on structure
        """
        #seq_list = list(self.Alignment.NamedSeqs[seq_id])
        seq_list = list(sequence)
        seq_len = len(seq_list)
        seq_mask = zeros(seq_len)
        mod_id_map = {}
        indices = []
        colors = []
        gapped,ungapped = self.MolType.gapMaps(sequence)
        self.GapMap[seq_id]= gapped
        if self.Strict:
            if seq_id not in self.ModuleMap:
                raise IndexError, 'seq_id %s not in ModuleMap'%(seq_id)
        else:
            if seq_id not in self.ModuleMap:
                return [],[]

        for mod_tup in self.ModuleMap[seq_id]:
            ix, mod_id, mod_len = mod_tup
            
            # skip modules we con't care about
            if not self.KeepAll and mod_id not in self.KeepIds:
                continue
            elif mod_id in self.SkipIds:
                continue

            mod_mask = zeros(seq_len)

            # mask motif region
            for i in range(ix,ix+mod_len):
                gapped_ix = self.GapMap[seq_id][i]
                mod_mask[gapped_ix] = 1
            # add to sequence map
            seq_mask += mod_mask

            # map module ids to indexes
            for jx in range(ix,ix+mod_len):
                gapped_jx = self.GapMap[seq_id][jx]
                if gapped_jx not in mod_id_map:
                    mod_id_map[gapped_jx] = [] 
                mod_id_map[gapped_jx].append(mod_id)


        # get module regions
        #need to take [0] element of nonzero() since numpy returns a tuple
        # where Numeric did not
        for kx in nonzero(seq_mask)[0]:
            # if overlapping use red background, otherwise display color
            if seq_mask[kx] > 1:
                curr_color = 'overlap_color'
            else:
                mod_id = mod_id_map[kx][0]
                curr_color = 'color_'+mod_id
            #append indices.
            indices.append(kx)
            colors.append(self.ColorMap[curr_color])
            
        return indices, colors

    def getMarkedIndices(self, seq_id, sequence, motif_id):
        """Returns list of indices to be marked for a given sequence.

            seq_id: sequence ID to highlight on structure
            sequence: sequence string
            motif_id: ID of motif to circle.
        """
        seq_list = list(sequence)
        seq_len = len(seq_list)
        seq_mask = zeros(seq_len)
        mod_id_map = {}
        indices = []

        gapped,ungapped = self.MolType.gapMaps(sequence)
        self.GapMap[seq_id]= gapped
        if self.Strict:
            if seq_id not in self.ModuleMap:
                raise IndexError, 'seq_id %s not in ModuleMap'%(seq_id)
        else:
            if seq_id not in self.ModuleMap:
                return [],[]

        for mod_tup in self.ModuleMap[seq_id]:
            ix, mod_id, mod_len = mod_tup
            if mod_id != motif_id:
                continue
            
            # skip modules we con't care about
            if not self.KeepAll and mod_id not in self.KeepIds:
                continue

            mod_mask = zeros(seq_len)

            # mask motif region
            for i in range(ix,ix+mod_len):
                gapped_ix = self.GapMap[seq_id][i]
                mod_mask[gapped_ix] = 1
            # add to sequence map
            seq_mask += mod_mask

            # map module ids to indexes
            for jx in range(ix,ix+mod_len):
                gapped_jx = self.GapMap[seq_id][jx]
                if gapped_jx not in mod_id_map:
                    mod_id_map[gapped_jx] = [] 
                mod_id_map[gapped_jx].append(mod_id)


        # get module regions
        #need to take [0] element of nonzero() since numpy returns a tuple
        # where Numeric did not
        for kx in nonzero(seq_mask)[0]:
            #append indices.
            indices.append(kx)
            
        return indices

class MotifsUpperCase(MotifFormatter):
    """Generates postscript file with motifs highlighted on 2D structure  """

    def makeModuleMap(self, motif_results):
        """
        Need to extract this b/c can't pickle motif_results... grr.

        motif_results: MotifResults object
        keep_module_ids: list of module ids to keep
        """
        module_map = {}  #Dict with locations of every motif keyed by module
        if motif_results:
            for motif in motif_results.Motifs:
                for module in motif.Modules:
                    mod_len = len(module)
                    mod_id = str(module.ID)
                    for skey, indexes in module.LocationDict.items():
                        if skey not in module_map:
                            module_map[skey] = []
                        for ix in indexes:
                            module_map[skey].append((ix, mod_id, mod_len))
        return module_map


    def __init__(self, MotifResults, KeepIds=None,\
        KeepAll=True, MolType=RNA):
        """Set up color map and motif results

        ModuleMap: flattened map (b/c of pickle problem.)
                generate using make_module_map() function
        Alignment: SequenceCollection or Alignment object
        KeepIds: list of module ids to keep
        KeepAll: When True, ignores KeepIds and highlights all motifs
        """
        MotifFormatter.__init__(self)
        self.ModuleMap = self.makeModuleMap(MotifResults)
        
        self.Alignment = MotifResults.Alignment
        if KeepIds is None:
            KeepIds = []
        self.KeepIds = set(KeepIds)
        self.KeepAll = KeepAll
        self.MolType = MolType
        self.GapMap = self.getGapMap()
        self.HighlightMap = {}

    def __call__(self,aln,outfile_prefix='',write_dir='.'):
        """Call method for ColorSecondaryStructure class.  """
        aln = SequenceCollection(aln)
        new_aln = {}
        for k,v in aln.NamedSeqs.items():
            new_aln[k]=self.getUpperCaseSequence(k,v)
        
        
        out_path = \
            write_dir+'/'+outfile_prefix+'_binding_site_alignment.fasta'
        out_file = open(out_path,'w')
        out_file.write(fasta_from_alignment(new_aln))
        out_file.close()
        print out_path
        return new_aln
            
    def getGapMap(self):
        """Returns dict mapping gapped_coord to ungapped_coord in self.Alignment
        
            - {seq_id:{gapped_coord:ungapped_coord}}
        """
        gap_map = {}
        for k,v in self.Alignment.items():
            gapped, ungapped = self.MolType.gapMaps(v)
            gap_map[k] = gapped
        return gap_map

    def getUpperCaseSequence(self, seq_id, sequence):
        """Returns list of indices and colors for a given sequence.

        seq_id: sequence ID to highlight on structure
        """
        sequence=str(sequence)
        seq_len = len(sequence)
        seq_mask = zeros(seq_len)
        mod_id_map = {}
        
        gapped,ungapped = self.MolType.gapMaps(sequence)
        self.GapMap[seq_id] = gapped
        
        #if there are no motifs, return lower case sequence.
        if seq_id not in self.ModuleMap:
            return sequence.lower()

        for mod_tup in self.ModuleMap[seq_id]:
            ix, mod_id, mod_len = mod_tup
            
            # skip modules we con't care about
            if not self.KeepAll and mod_id not in self.KeepIds:
                continue

            mod_mask = zeros(seq_len)

            # mask motif region
            for i in range(ix,ix+mod_len):
                gapped_ix = self.GapMap[seq_id][i]
                mod_mask[gapped_ix] = 1
            # add to sequence map
            seq_mask += mod_mask
        
        # get upper case in module regions
        new_seq = []
        for kx, lc, uc in zip(seq_mask,sequence.lower(), sequence.upper()):
            # if not motif region, use lower
            if kx < 1:
                new_seq.append(lc)
            else:
                new_seq.append(uc)
            
        return ''.join(new_seq)
    
class MotifSequenceConstraints(MotifFormatter):
    """Generates postscript file with motifs highlighted on 2D structure  """

    def makeModuleMap(self, motif_results):
        """
        Need to extract this b/c can't pickle motif_results... grr.

        motif_results: MotifResults object
        keep_module_ids: list of module ids to keep
        """
        module_map = {}  #Dict with locations of every motif keyed by module
        if motif_results:
            for motif in motif_results.Motifs:
                for module in motif.Modules:
                    mod_len = len(module)
                    mod_id = str(module.ID)
                    for skey, indexes in module.LocationDict.items():
                        if skey not in module_map:
                            module_map[skey] = []
                        for ix in indexes:
                            module_map[skey].append((ix, mod_id, mod_len))
        return module_map


    def __init__(self, MotifResults, KeepIds=None,\
        KeepAll=True, MolType=RNA, strict=True):
        """Set up color map and motif results

        ModuleMap: flattened map (b/c of pickle problem.)
                generate using make_module_map() function
        Alignment: SequenceCollection or Alignment object
        KeepIds: list of module ids to keep
        KeepAll: When True, ignores KeepIds and highlights all motifs
        """
        MotifFormatter.__init__(self)
        self.ModuleMap = self.makeModuleMap(MotifResults)
        
        self.Alignment = MotifResults.Alignment
        if KeepIds is None:
            KeepIds = []
        self.KeepIds = set(KeepIds)
        self.KeepAll = KeepAll
        self.MolType = MolType
        self.GapMap = self.getGapMap()
        self.Strict=strict
        self.MotifCharacter = \
        list('0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')

    def __call__(self,seq_id,sequence,struct,write_dir='.'):
        """Call method for ColorSecondaryStructure class.  """

        #indices, colors = self.getColorIndices(seq_id,sequence)
        indices, colors = self.getColorIndices(seq_id,sequence)
        circle_indices = []
        if self.CircleId:
            circle_indices = \
                self.getCircleIndices(seq_id,sequence,self.CircleId)
        
        structure_postscript = color_on_structure(\
            sequence=sequence,\
            struct=struct,\
            color_map = self.ColorMap,\
            indices=indices,\
            colors=colors,\
            circle_indices=circle_indices)
        file_id = seq_id.split()[0]
        ps_out_path = \
            write_dir+'/'+file_id+'_secondary_struct.ps'
        ps_out = open(ps_out_path,'w')
        ps_out.write(structure_postscript)
        ps_out.close()
        return ps_out_path
            
    def getGapMap(self):
        """Returns dict mapping gapped_coord to ungapped_coord in self.Alignment
        
            - {seq_id:{gapped_coord:ungapped_coord}}
        """
        gap_map = {}
        for k,v in self.Alignment.items():
            gapped, ungapped = self.MolType.gapMaps(v)
            gap_map[k] = gapped
        return gap_map

    def getSeqMask(self, seq_id, sequence):
        """Returns vector where motifs are present in sequence.

            - seq_id: sequence ID to make seq mask.
            - sequence: sequence to make seq mask with.
        """
        #seq_list = list(self.Alignment.NamedSeqs[seq_id])
        seq_list = list(sequence)
        seq_len = len(seq_list)
        seq_mask = zeros(seq_len)
        mod_id_map = {}

        gapped,ungapped = self.MolType.gapMaps(sequence)
        self.GapMap[seq_id]= gapped
        if self.Strict:
            if seq_id not in self.ModuleMap:
                raise IndexError, 'seq_id %s not in ModuleMap'%(seq_id)
        else:
            if seq_id not in self.ModuleMap:
                return '',''

        for mod_tup in self.ModuleMap[seq_id]:
            ix, mod_id, mod_len = mod_tup
            
            # skip modules we con't care about
            if not self.KeepAll and mod_id not in self.KeepIds:
                continue
            elif mod_id in self.SkipIds:
                continue

            mod_mask = zeros(seq_len)

            # mask motif region
            for i in range(ix,ix+mod_len):
                gapped_ix = self.GapMap[seq_id][i]
                mod_mask[gapped_ix] = 1
            # add to sequence map
            seq_mask += mod_mask

        return seq_mask
        
    def getOverlapDicts(self,alignment):
        """Returns dicts of motifs that overlap in start or end positions.
        """
        #start overlap dict
        start_overlap = {}
        #end overlap dict
        end_overlap = {}
        
        #Get seq_mask_dict.  Calling this will construct self.GapMap for given
        # sequence.
        for seq_id,seq in alignment.items():
            curr_seq_mask = self.getSeqMask(seq_id,seq)
                    
            #for each module
            for mod_tup in self.ModuleMap[seq_id]:
                ix, mod_id, mod_len = mod_tup
                gapped_start = self.GapMap[seq_id][ix]
                gapped_end = self.GapMap[seq_id][ix+mod_len]
                if curr_seq_mask[gapped_start]>1:
                    start_overlap[mod_id]=seq_id
                if curr_seq_mask[gapped_end]>1:
                    end_overlap[mod_id]=seq_id
        
        return start_overlap, end_overlap
    
    def getConstraintStrings(self,alignment):
        """Returns dict of constraint strings for each sequence in alignment.
        
            - {seq_id:{1:constraint_string_1,2:constraint_string_2}}
        """
        start_overlap, end_overlap = self.getOverlapDicts(alignment)
        
        
