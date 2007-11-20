#!/usr/bin/env python
"""Provides an application controller for the commandline version of:
CLUSTALW
"""
from cogent.app.parameters import FlagParameter, ValuedParameter, \
    MixedParameter
from cogent.app.util import CommandLineApplication, ResultPath

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Sandra Smit", "Micah Hamady", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"

class Clustalw(CommandLineApplication):
    """ clustalw application controller 
   
    The parameters are organized by function to give some idea of how the 
    program works. However, no restrictions are put on any combinations 
    of parameters. Misuse of parameters can lead to errors or otherwise
    strange results.

    You are supposed to choose one action for the program to perform. (align, 
    profile, sequences, tree, or bootstrap). If you choose multiple, only the
    dominant action (see order above) will be executed. By DEFAULT, the -align
    parameter is turned on. If you decide to turn another one on, you should 
    turn '-align' off IN ADDITION!
    
    Some references to help pages are available in the 'getHelp' method.
    Some might be useful to you.
    """
    _actions = {\
        '-align':FlagParameter('-','align',Value=True),
        '-profile':FlagParameter('-','profile'),
        '-sequences':FlagParameter('-','sequences'),
        '-tree':FlagParameter('-','tree'),
        '-bootstrap':MixedParameter('-','bootstrap',Delimiter='=')}

    #sequence file for alignment, or alignment file for bootstrap and tree
    #actions
    _input = {'-infile':ValuedParameter('-','infile',Delimiter='=',IsPath=True)}
   
    # matrix and dnamatrix can be filenames as well, but not always.
    # They won't be treated as filenames and thus not quoted.
    # Therefore filepaths containing spaces might result in errors.
    _multiple_alignment={\
        '-quicktree':FlagParameter('-','quicktree'),
        '-type':ValuedParameter('-','type',Delimiter='='),
        '-matrix':ValuedParameter('-','matrix',Delimiter='='),
        '-dnamatrix':ValuedParameter('-','dnamatrix',Delimiter='='),
        '-gapopen':ValuedParameter('-','gapopen',Delimiter='='),
        '-gapext':ValuedParameter('-','gapext',Delimiter='='),
        '-endgaps':FlagParameter('-','endgaps'),
        '-gapdist':ValuedParameter('-',Name='gapdist',Delimiter='='),
        '-nopgap':FlagParameter('-','nopgap'),
        '-nohgap':FlagParameter('-','nohgap'),
        '-hgapresidues':ValuedParameter('-','hgapresidues',Delimiter='='),
        '-maxdiv':ValuedParameter('-',Name='maxdiv',Delimiter='='),
        '-negative':FlagParameter('-','negative'),
        '-transweight':ValuedParameter('-',Name='transweight',Delimiter='='),
        '-newtree':ValuedParameter('-','newtree',Delimiter='=',IsPath=True),
        '-usetree':ValuedParameter('-','usetree',Delimiter='=',IsPath=True)}

    _fast_pairwise={\
        '-ktuple':ValuedParameter('-',Name='ktuple',Delimiter='='),
        '-topdiags':ValuedParameter('-',Name='topdiags',Delimiter='='),
        '-window':ValuedParameter('-',Name='window',Delimiter='='),
        '-pairgap':ValuedParameter('-',Name='pairgap',Delimiter='='),
        '-score':ValuedParameter('-',Name='score',Delimiter='=')}
        
    # pwmatrix and pwdnamatrix can be filenames as well, but not always.
    # They won't be treated as filenames and thus not quoted.
    # Therefore filepaths containing spaces might result in errors.
    _slow_pairwise={\
        '-pwmatrix':ValuedParameter('-',Name='pwmatrix',Delimiter='='),
        '-pwdnamatrix':ValuedParameter('-',Name='pwdnamatrix',Delimiter='='),
        '-pwgapopen':ValuedParameter('-',Name='pwgapopen',Delimiter='='),
        '-pwgapext':ValuedParameter('-',Name='pwgapext',Delimiter='=')}

    #plus -bootstrap
    _tree={\
        '-kimura':FlagParameter('-',Name='kimura'),
        '-tossgaps':FlagParameter('-',Name='tossgaps'),
        '-bootlabels':ValuedParameter('-',Name='bootlabels',Delimiter='='),
        '-seed':ValuedParameter('-',Name='seed',Delimiter='='),
        '-outputtree':ValuedParameter('-',Name='outputtree',Delimiter='=')}

    _output={\
        '-outfile':ValuedParameter('-',Name='outfile',Delimiter='=',\
            IsPath=True),
        '-output':ValuedParameter('-',Name='output',Delimiter='='),
        '-case':ValuedParameter('-',Name='case',Delimiter='='),
        '-outorder':ValuedParameter('-',Name='outorder',Delimiter='='),
        '-seqnos':ValuedParameter('-',Name='seqnos',Delimiter='=')}

    _profile_alignment={\
        '-profile1':ValuedParameter('-','profile1',Delimiter='=',IsPath=True),
        '-profile2':ValuedParameter('-','profile2',Delimiter='=',IsPath=True),
        '-usetree1':ValuedParameter('-','usetree1',Delimiter='=',IsPath=True),
        '-usetree2':ValuedParameter('-','usetree2',Delimiter='=',IsPath=True),
        '-newtree1':ValuedParameter('-','newtree1',Delimiter='=',IsPath=True),
        '-newtree2':ValuedParameter('-','newtree2',Delimiter='=',IsPath=True)}
    
    _structure_alignment={\
        '-nosecstr1':FlagParameter('-',Name='nosecstr1'),
        '-nosecstr2':FlagParameter('-',Name='nosecstr2'),
        '-helixgap':ValuedParameter('-',Name='helixgap',Delimiter='='),
        '-strandgap':ValuedParameter('-',Name='strandgap',Delimiter='='),
        '-loopgap':ValuedParameter('-',Name='loopgap',Delimiter='='),
        '-terminalgap':ValuedParameter('-',Name='terminalgap',Delimiter='='),
        '-helixendin':ValuedParameter('-',Name='helixendin',Delimiter='='),
        '-helixendout':ValuedParameter('-',Name='helixendout',Delimiter='='),
        '-strandendin':ValuedParameter('-',Name='strandendin',Delimiter='='),
        '-strandendout':ValuedParameter('-',Name='strandendout',Delimiter='='),
        '-secstrout':ValuedParameter('-',Name='secstrout',Delimiter='=')}
    
        #NOT SUPPORTED
        #'-help':FlagParameter('-','help'),
        #'-check':FlagParameter('-','check'),
        #'-options':FlagParameter('-','options'),
        #'-convert':FlagParameter('-','convert'),
        #'-batch':FlagParameter('-','batch'),
        #'-noweights':FlagParameter('-','noweights'),
        #'-novgap':FlagParameter('-','novgap'),
        #'-debug':ValuedParameter('-',Name='debug',Delimiter='='),

    _parameters = {}
    _parameters.update(_actions)
    _parameters.update(_input)
    _parameters.update(_multiple_alignment)
    _parameters.update(_fast_pairwise)
    _parameters.update(_slow_pairwise)
    _parameters.update(_tree)
    _parameters.update(_output)
    _parameters.update(_profile_alignment)
    _parameters.update(_structure_alignment)
 
    _command = 'clustalw'
   
    def getHelp(self):
        """Methods that points to the documentation"""
        help_str =\
        """
        There are several help pages available online. For example:
        http://searchlauncher.bcm.tmc.edu/multi-align/Help/
            clustalw_help_1.8.html
        http://hypernig.nig.ac.jp/homology/clustalw-e_help.html
        http://www.genebee.msu.su/clustal/help.html
        
        A page that give reasonable insight in use of the parameters:
        http://bioweb.pasteur.fr/seqanal/interfaces/clustalw.html
        """
        return help_str
   
    def _input_as_lines(self,data):
        """Writes data to tempfile and sets -infile parameter

        data -- list of lines, ready to be written to file
        """
        if data:
            self.Parameters['-infile']\
                .on(super(Clustalw,self)._input_as_lines(data))
        return ''

    def _input_as_seqs(self,data):
        """writes sequences to tempfile and sets -infile parameter

        data -- list of sequences

        Adds numbering to the sequences: >1, >2, etc.
        """
        lines = []
        for i,s in enumerate(data):
            #will number the sequences 1,2,3,etc.
            lines.append(''.join(['>',str(i+1)]))
            lines.append(s)
        return self._input_as_lines(lines)

    def _input_as_string(self,data):
        """Makes data the value of a specific parameter
    
        This method returns the empty string. The parameter will be printed
        automatically once set.
        """
        if data:
            self.Parameters['-infile'].on(data)
        return ''

    def _suffix(self):
        """Return appropriate suffix for alignment file"""
        _output_formats={'GCG':'.msf',
                        'GDE':'.gde',
                        'PHYLIP':'.phy',
                        'PIR':'.pir',
                        'NEXUS':'.nxs'}

        if self.Parameters['-output'].isOn():
            return _output_formats[self.Parameters['-output'].Value]
        else:
            return '.aln'
    
    def _aln_filename(self,prefix):
        """Return name of file containing the alignment
        
        prefix -- str, prefix of alignment file.
        """
        if self.Parameters['-outfile'].isOn():
            aln_filename = self._absolute(self.Parameters['-outfile'].Value)
        else:
            aln_filename = prefix + self._suffix()
        return aln_filename

    def _get_result_paths(self,data):
        """Return dict of {key: ResultPath}
        """
        
        #clustalw .aln is used when no or unkown output type specified
        _treeinfo_formats = {'nj':'.nj',
                            'dist':'.dst',
                            'nexus':'.tre'}

        result = {}
        par = self.Parameters
        abs = self._absolute
        
        if par['-align'].isOn():
            prefix = par['-infile'].Value.rsplit('.', 1)[0]
            #prefix = par['-infile'].Value.split('.')[0]
            aln_filename = self._aln_filename(prefix)
            if par['-newtree'].isOn():
                dnd_filename = abs(par['-newtree'].Value)
            elif par['-usetree'].isOn():
                dnd_filename = abs(par['-usetree'].Value)
            else:
                dnd_filename = abs(prefix + '.dnd')
            result['Align'] = ResultPath(Path=aln_filename,IsWritten=True)
            result['Dendro'] = ResultPath(Path=dnd_filename,IsWritten=True)
        elif par['-profile'].isOn():
            prefix1 = par['-profile1'].Value.rsplit('.', 1)[0]
            prefix2 = par['-profile2'].Value.rsplit('.', 1)[0]
            #prefix1 = par['-profile1'].Value.split('.')[0]
            #prefix2 = par['-profile2'].Value.split('.')[0]
            aln_filename = ''; aln_written = True
            dnd1_filename = ''; tree1_written = True
            dnd2_filename = ''; tree2_written = True
            aln_filename = self._aln_filename(prefix1)
            #usetree1
            if par['-usetree1'].isOn():
                tree1_written = False
            #usetree2
            if par['-usetree2'].isOn():
                tree2_written = False
            if par['-newtree1'].isOn():
                dnd1_filename = abs(par['-newtree1'].Value)
                aln_written=False
            else:
                dnd1_filename = abs(prefix1 + '.dnd')
            if par['-newtree2'].isOn():
                dnd2_filename = abs(par['-newtree2'].Value)
                aln_written=False
            else:
                dnd2_filename = abs(prefix2 + '.dnd')
            result['Align'] = ResultPath(Path=aln_filename,
                IsWritten=aln_written)
            result['Dendro1'] = ResultPath(Path=dnd1_filename,
                IsWritten=tree1_written)
            result['Dendro2'] = ResultPath(Path=dnd2_filename,
                IsWritten=tree2_written)
        elif par['-sequences'].isOn():
            prefix1 = par['-profile1'].Value.rsplit('.', 1)[0]
            prefix2 = par['-profile2'].Value.rsplit('.', 1)[0]
            #prefix1 = par['-profile1'].Value.split('.')[0] #alignment
            #prefix2 = par['-profile2'].Value.split('.')[0] #sequences
            aln_filename = ''; aln_written = True
            dnd_filename = ''; dnd_written = True
            
            aln_filename = self._aln_filename(prefix2)
            if par['-usetree'].isOn():
                dnd_written = False
            elif par['-newtree'].isOn():
                aln_written = False
                dnd_filename = abs(par['-newtree'].Value)
            else:
                dnd_filename = prefix2 + '.dnd'  
            result['Align'] = ResultPath(Path=aln_filename,\
                IsWritten=aln_written)
            result['Dendro'] = ResultPath(Path=dnd_filename,\
                IsWritten=dnd_written)
        elif par['-tree'].isOn():
            prefix = par['-infile'].Value.rsplit('.', 1)[0]
            #prefix = par['-infile'].Value.split('.')[0]
            tree_filename = ''; tree_written = True
            treeinfo_filename = ''; treeinfo_written = False
            tree_filename = prefix + '.ph'
            if par['-outputtree'].isOn() and\
                par['-outputtree'].Value != 'phylip':
                treeinfo_filename = prefix +\
                    _treeinfo_formats[par['-outputtree'].Value]
                treeinfo_written = True
            result['Tree'] = ResultPath(Path=tree_filename,\
                IsWritten=tree_written)
            result['TreeInfo'] = ResultPath(Path=treeinfo_filename,\
                IsWritten=treeinfo_written)
            
        elif par['-bootstrap'].isOn():
            prefix = par['-infile'].Value.rsplit('.', 1)[0]
            #prefix = par['-infile'].Value.split('.')[0]   
            boottree_filename = prefix + '.phb'
            result['Tree'] = ResultPath(Path=boottree_filename,IsWritten=True)
        
        return result

        
#SOME FUNCTIONS TO EXECUTE THE MOST COMMON TASKS
def alignUnalignedSeqs(seqs,add_seq_names=True,WorkingDir=None,\
    SuppressStderr=None,SuppressStdout=None):
    """Aligns unaligned sequences

    seqs: either list of sequence objects or list of strings
    add_seq_names: boolean. if True, sequence names are inserted in the list
        of sequences. if False, it assumes seqs is a list of lines of some
        proper format that the program can handle
    """
    if add_seq_names:
        app = Clustalw(InputHandler='_input_as_seqs',\
            WorkingDir=WorkingDir,SuppressStderr=SuppressStderr,\
            SuppressStdout=SuppressStdout)
    else:
        app = Clustalw(InputHandler='_input_as_lines',\
            WorkingDir=WorkingDir,SuppressStderr=SuppressStderr,\
            SuppressStdout=SuppressStdout)
    return app(seqs)

def alignUnalignedSeqsFromFile(filename,WorkingDir=None,SuppressStderr=None,\
    SuppressStdout=None):
    """Aligns unaligned sequences from some file (file should be right format)

    filename: string, the filename of the file containing the sequences
        to be aligned in a valid format.
    """
    app = Clustalw(WorkingDir=WorkingDir,SuppressStderr=SuppressStderr,\
        SuppressStdout=SuppressStdout)
    return app(filename)

def alignTwoAlignments(aln1,aln2,outfile,WorkingDir=None,SuppressStderr=None,\
    SuppressStdout=None):
    """Aligns two alignments. Individual sequences are not realigned
    
    aln1: string, name of file containing the first alignment
    aln2: string, name of file containing the second alignment
    outfile: you're forced to specify an outfile name, because if you don't 
        aln1 will be overwritten. So, if you want aln1 to be overwritten, you 
        should specify the same filename.
    WARNING: a .dnd file is created with the same prefix as aln1. So an 
    existing dendrogram might get overwritten.
    """
    app = Clustalw({'-profile':None,'-profile1':aln1,\
        '-profile2':aln2,'-outfile':outfile},SuppressStderr=\
        SuppressStderr,WorkingDir=WorkingDir,SuppressStdout=SuppressStdout)
    app.Parameters['-align'].off()
    return app()

def addSeqsToAlignment(aln1,seqs,outfile,WorkingDir=None,SuppressStderr=None,\
        SuppressStdout=None):
    """Aligns sequences from second profile against first profile
    
    aln1: string, name of file containing the alignment
    seqs: string, name of file containing the sequences that should be added
        to the alignment.
    outfile: string, name of the output file (the new alignment)
    """
    app = Clustalw({'-sequences':None,'-profile1':aln1,\
        '-profile2':seqs,'-outfile':outfile},SuppressStderr=\
        SuppressStderr,WorkingDir=WorkingDir, SuppressStdout=SuppressStdout)
        
    app.Parameters['-align'].off()
    return app()

def buildTreeFromAlignment(filename,WorkingDir=None,SuppressStderr=None):
    """Builds a new tree from an existing alignment
    
    filename: string, name of file containing the seqs or alignment
    """
    app = Clustalw({'-tree':None,'-infile':filename},SuppressStderr=\
        SuppressStderr,WorkingDir=WorkingDir)
    app.Parameters['-align'].off()
    return app()

def align_and_build_tree(seqs, best_tree=False, params=None):
    """Returns an alignment and a tree from Sequences object seqs.
    
    seqs: an xxx.Sequences object, or data that can be used to build one.
    
    best_tree: if True (default:False), uses a slower but more accurate
    algorithm to build the tree.

    params: dict of parameters to pass in to the Clustal app controller.

    The result will be a tuple containing a xxx.Alignment and an xxx.Tree
    object (or None for the alignment and/or tree if either fails).
    """
    raise NotImplementedError
    
def build_tree_from_alignment(aln, best_tree=False, params=None):
    """Returns a tree from Alignment object aln.

    aln: an xxx.Alignment object, or data that can be used to build one.

    best_tree: if True (default:False), uses a slower but more accurate
    algorithm to build the tree.

    params: dict of parameters to pass in to the Clustal app controller.

    The result will be an xxx.Alignment object, or None if tree fails.
    """
    raise NotImplementedError
    
def add_seqs_to_alignment(seqs, aln, params=None):
    """Returns an Alignment object from seqs and existing Alignment.

    seqs: an xxx.Sequences object, or data that can be used to build one.

    aln: an xxx.Alignment object, or data that can be used to build one

    params: dict of parameters to pass in to the Clustal app controller.
    """
    raise NotImplementedError

def align_two_alignments(aln1, aln2, params=None):
    """Returns an Alignment object from two existing Alignments.

    aln1, aln2: xxx.Alignment objects, or data that can be used to build them.

    params: dict of parameters to pass in to the Clustal app controller.
    """
    raise NotImplementedError
