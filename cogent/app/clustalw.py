#!/usr/bin/env python
"""Provides an application controller for the commandline version of:
CLUSTALW v1.83
"""
from cogent.app.parameters import FlagParameter, ValuedParameter, \
    MixedParameter, FilePath
from cogent.app.util import CommandLineApplication, ResultPath, remove
from cogent.core.alignment import SequenceCollection, Alignment
from cogent.parse.tree import DndParser 
from cogent.parse.clustal import ClustalParser
from cogent.core.tree import PhyloNode
from cogent.core.moltype import RNA, DNA, PROTEIN
from numpy.random import randint

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Sandra Smit", "Micah Hamady", "Rob Knight", "Jeremy Widmann",
                "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
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
   
    def _input_as_multiline_string(self, data):
        """Writes data to tempfile and sets -infile parameter

        data -- list of lines
        """
        if data:
            self.Parameters['-infile']\
                .on(super(Clustalw,self)._input_as_multiline_string(data))
        return ''

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
    
    def _tempfile_as_multiline_string(self, data):
        """Write a multiline string to a temp file and return the filename.

            data: a multiline string to be written to a file.

           * Note: the result will be the filename as a FilePath object 
            (which is a string subclass).

        """
        filename = FilePath(self.getTmpFilename(self.TmpDir))
        data_file = open(filename,'w')
        data_file.write(data)
        data_file.close()
        return filename

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
    opoutfile: string, name of the output file (the new alignment)
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

def align_and_build_tree(seqs, moltype, best_tree=False, params=None):
    """Returns an alignment and a tree from Sequences object seqs.
    
    seqs: an cogent.core.alignment.SequenceCollection object, or data that can
    be used to build one.
    
    moltype: cogent.core.moltype.MolType object

    best_tree: if True (default:False), uses a slower but more accurate
    algorithm to build the tree.

    params: dict of parameters to pass in to the Clustal app controller.

    The result will be a tuple containing a cogent.core.alignment.Alignment
    and a cogent.core.tree.PhyloNode
    object (or None for the alignment and/or tree if either fails).
    """
    aln = align_unaligned_seqs(seqs, moltype=moltype, params=params)
    tree = build_tree_from_alignment(aln, moltype, best_tree, params)
    return {'Align':aln,'Tree':tree}
    
def build_tree_from_alignment(aln, moltype, best_tree=False, params=None):
    """Returns a tree from Alignment object aln.

    aln: an cogent.core.alignment.Alignment object, or data that can be used
    to build one.

    moltype: cogent.core.moltype.MolType object

    best_tree: if True (default:False), uses a slower but more accurate
    algorithm to build the tree.

    params: dict of parameters to pass in to the Clustal app controller.

    The result will be an cogent.core.tree.PhyloNode object, or None if tree
    fails.
    """
    # Create instance of app controller, enable tree, disable alignment
    app = Clustalw(InputHandler='_input_as_multiline_string', params=params, \
                   WorkingDir='/tmp')
    app.Parameters['-align'].off()
    
    #Set params to empty dict if None.
    if params is None:
        params={}

    if moltype == DNA or moltype == RNA:
        params['-type'] = 'd'
    elif moltype == PROTEIN:
        params['-type'] = 'p'
    else:
        raise ValueError, "moltype must be DNA, RNA, or PROTEIN"

    # best_tree -> bootstrap
    if best_tree:
        if '-bootstrap' not in params:
            app.Parameters['-bootstrap'].on(1000)
        if '-seed' not in params:
            app.Parameters['-seed'].on(randint(0,1000))
        if '-bootlabels' not in params:
            app.Parameters['-bootlabels'].on('nodes')
    else:
        app.Parameters['-tree'].on()

    # Setup mapping. Clustalw clips identifiers. We will need to remap them.
    seq_collection = SequenceCollection(aln)
    int_map, int_keys = seq_collection.getIntMap()
    int_map = SequenceCollection(int_map)
    
    # Collect result
    result = app(int_map.toFasta())

    # Build tree
    tree = DndParser(result['Tree'].read(), constructor=PhyloNode)
    for node in tree.tips():
        node.Name = int_keys[node.Name]

    # Clean up
    result.cleanUp()
    del(seq_collection, app, result, int_map, int_keys)

    return tree
    
def bootstrap_tree_from_alignment(aln, seed=None, num_trees=None, params=None):
    """Returns a tree from Alignment object aln with bootstrap support values.

    aln: an cogent.core.alignment.Alignment object, or data that can be used
    to build one.

    seed: an interger, seed value to use
    
    num_trees: an integer, number of trees to bootstrap against

    params: dict of parameters to pass in to the Clustal app controller.

    The result will be an cogent.core.tree.PhyloNode object, or None if tree
    fails.

    If seed is not specifed in params, a random integer between 0-1000 is used.
    """
    # Create instance of controllor, enable bootstrap, disable alignment,tree
    app = Clustalw(InputHandler='_input_as_multiline_string', params=params, \
                   WorkingDir='/tmp')
    app.Parameters['-align'].off()
    app.Parameters['-tree'].off()

    if app.Parameters['-bootstrap'].isOff():
        if num_trees is None:
            num_trees = 1000

        app.Parameters['-bootstrap'].on(num_trees)

    if app.Parameters['-seed'].isOff():
        if seed is None:
            seed = randint(0,1000)

        app.Parameters['-seed'].on(seed)

    if app.Parameters['-bootlabels'].isOff():
        app.Parameters['-bootlabels'].on("node")

    # Setup mapping. Clustalw clips identifiers. We will need to remap them.
    seq_collection = SequenceCollection(aln)
    int_map, int_keys = seq_collection.getIntMap()
    int_map = SequenceCollection(int_map)
    
    # Collect result
    result = app(int_map.toFasta())

    # Build tree
    tree = DndParser(result['Tree'].read(), constructor=PhyloNode)
    for node in tree.tips():
        node.Name = int_keys[node.Name]

    # Clean up
    result.cleanUp()
    del(seq_collection, app, result, int_map, int_keys)

    return tree

def align_unaligned_seqs(seqs, moltype, params=None):
    """Returns an Alignment object from seqs.

    seqs: cogent.core.alignment.SequenceCollection object, or data that can be
    used to build one.
    
    moltype: a MolType object.  DNA, RNA, or PROTEIN.

    params: dict of parameters to pass in to the Clustal app controller.
    
    Result will be a cogent.core.alignment.Alignment object.
    """
    #create SequenceCollection object from seqs
    seq_collection = SequenceCollection(seqs,MolType=moltype)
    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seq_collection.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)
    #Create Clustalw app.
    app = Clustalw(InputHandler='_input_as_multiline_string',params=params)
    #Get results using int_map as input to app
    res = app(int_map.toFasta())
    #Get alignment as dict out of results
    alignment = dict(ClustalParser(res['Align'].readlines()))
    #Make new dict mapping original IDs
    new_alignment = {}
    for k,v in alignment.items():
        new_alignment[int_keys[k]]=v
    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment,MolType=moltype)
    #Clean up
    res.cleanUp()
    del(seq_collection,int_map,int_keys,app,res,alignment)

    return new_alignment

def add_seqs_to_alignment(seqs, aln, moltype, params=None):
    """Returns an Alignment object from seqs and existing Alignment.

    seqs: a cogent.core.alignment.SequenceCollection object, or data that can
    be used to build one.

    aln: a cogent.core.alignment.Alignment object, or data that can be used to
    build one

    params: dict of parameters to pass in to the Clustal app controller.
    """
    #create SequenceCollection object from seqs
    seq_collection = SequenceCollection(seqs,MolType=moltype)
    #Create mapping between abbreviated IDs and full IDs
    seq_int_map, seq_int_keys = seq_collection.getIntMap()
    #Create SequenceCollection from int_map.
    seq_int_map = SequenceCollection(seq_int_map,MolType=moltype)
    
    #create Alignment object from aln
    aln = Alignment(aln,MolType=moltype)
    #Create mapping between abbreviated IDs and full IDs
    aln_int_map, aln_int_keys = aln.getIntMap(prefix='seqn_')
    #Create SequenceCollection from int_map.
    aln_int_map = Alignment(aln_int_map,MolType=moltype)
    
    #Update seq_int_keys with aln_int_keys
    seq_int_keys.update(aln_int_keys)
    
    #Create Mafft app.
    app = Clustalw(InputHandler='_input_as_multiline_string',\
        params=params,
        SuppressStderr=True)
    app.Parameters['-align'].off()
    app.Parameters['-infile'].off()
    app.Parameters['-sequences'].on()
    
    #Add aln_int_map as profile1
    app.Parameters['-profile1'].on(\
        app._tempfile_as_multiline_string(aln_int_map.toFasta()))
    
    #Add seq_int_map as profile2
    app.Parameters['-profile2'].on(\
        app._tempfile_as_multiline_string(seq_int_map.toFasta()))
    #Get results using int_map as input to app
    res = app()
    
    #Get alignment as dict out of results
    alignment = dict(ClustalParser(res['Align'].readlines()))
    
    #Make new dict mapping original IDs
    new_alignment = {}
    for k,v in alignment.items():
        new_alignment[seq_int_keys[k]]=v
    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment,MolType=moltype)
    #Clean up
    res.cleanUp()
    remove(app.Parameters['-profile1'].Value)
    remove(app.Parameters['-profile2'].Value)
    del(seq_collection,seq_int_map,seq_int_keys,\
        aln,aln_int_map,aln_int_keys,app,res,alignment)

    return new_alignment

def align_two_alignments(aln1, aln2, moltype, params=None):
    """Returns an Alignment object from two existing Alignments.

    aln1, aln2: cogent.core.alignment.Alignment objects, or data that can be
    used to build them.

    params: dict of parameters to pass in to the Clustal app controller.
    """
    #create SequenceCollection object from seqs
    aln1 = Alignment(aln1,MolType=moltype)
    #Create mapping between abbreviated IDs and full IDs
    aln1_int_map, aln1_int_keys = aln1.getIntMap()
    #Create SequenceCollection from int_map.
    aln1_int_map = Alignment(aln1_int_map,MolType=moltype)
    
    #create Alignment object from aln
    aln2 = Alignment(aln2,MolType=moltype)
    #Create mapping between abbreviated IDs and full IDs
    aln2_int_map, aln2_int_keys = aln2.getIntMap(prefix='seqn_')
    #Create SequenceCollection from int_map.
    aln2_int_map = Alignment(aln2_int_map,MolType=moltype)
    
    #Update aln1_int_keys with aln2_int_keys
    aln1_int_keys.update(aln2_int_keys)
    
    #Create Mafft app.
    app = Clustalw(InputHandler='_input_as_multiline_string',\
        params=params,
        SuppressStderr=True)
    app.Parameters['-align'].off()
    app.Parameters['-infile'].off()
    app.Parameters['-profile'].on()
    
    #Add aln_int_map as profile1
    app.Parameters['-profile1'].on(\
        app._tempfile_as_multiline_string(aln1_int_map.toFasta()))
    
    #Add seq_int_map as profile2
    app.Parameters['-profile2'].on(\
        app._tempfile_as_multiline_string(aln2_int_map.toFasta()))
    #Get results using int_map as input to app
    res = app()
    
    #Get alignment as dict out of results
    alignment = dict(ClustalParser(res['Align'].readlines()))
    
    #Make new dict mapping original IDs
    new_alignment = {}
    for k,v in alignment.items():
        new_alignment[aln1_int_keys[k]]=v
    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment,MolType=moltype)
    #Clean up
    res.cleanUp()
    remove(app.Parameters['-profile1'].Value)
    remove(app.Parameters['-profile2'].Value)
    del(aln1,aln1_int_map,aln1_int_keys,\
        aln2,aln2_int_map,aln2_int_keys,app,res,alignment)

    return new_alignment
    
