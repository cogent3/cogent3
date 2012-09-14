#!/usr/bin/env python

"""Parses MEME output file and creates Module objects.

    Supports MEME version 3.0 - 4.8.1
"""
from cogent.parse.record_finder import LabeledRecordFinder
from cogent.parse.record import DelimitedSplitter
from cogent.motif.util import Location, ModuleInstance, Module, Motif,\
     MotifResults, make_remap_dict
from cogent.core.moltype import DNA, RNA, PROTEIN

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

def getDataBlock(lines):
    """Returns main block of data as list.
    """
    #Get main data block: All lines following "COMMAND LINE SUMMARY"
    meme_command = LabeledRecordFinder(lambda x: x.startswith('COMMAND'))
    main_block = list(meme_command(lines))
    alphabet = getMolType(main_block[0])
    return main_block[1], alphabet

def getMolType(lines):
    """Returns alphabet type that sequences belong to.
    """
    for line in lines:
        if 'ALPHABET' in line:
            alphabet_line = line
    #Split on equal sign
    alphabet_line = alphabet_line.strip().split('ALPHABET= ')[1].split()[0]
    #Get Set of alphabet letters
    alphabet = set(alphabet_line)
    #get Protein Set
    protein_order = set(PROTEIN.Alphabet)
    #get RNA Set
    rna_order = set(RNA.Alphabet)
    #Find out which alphabet is used
    
    if len(alphabet) >= 20:
        return PROTEIN
    elif alphabet == rna_order:
        return RNA
    else:
        return DNA
    

def getCommandModuleBlocks(main_block):
    """Returns command line summary block and list of module blocks.
    """
    #Get Command line summary and all module information
    meme_module = LabeledRecordFinder(lambda x: x.startswith('MOTIF'))
    main_block = list(meme_module(main_block))
    command_block = main_block[0]
    module_blocks = []
    if len(main_block) > 1:
        module_blocks = main_block[1:]
    return command_block, module_blocks

def getSummaryBlock(module_blocks):
    """Returns summary of motifs block.
    """
    meme_summary = LabeledRecordFinder(lambda x: x.startswith('SUMMARY'),\
        constructor=None,ignore=lambda x: x.startswith(' '))
    summary_block = list(meme_summary(module_blocks))
    return summary_block[1]

def dictFromList(data_list):
    """Returns a dict given a list.

        - Dict created from a list where list contains alternating key, value
        pairs.
        - ex: [key1, value1, key2, value2] returns: {key1:value1, key2:value2}
    """
    data_dict = {}
    for i in range(0,len(data_list)-1,2):
        #If there is already a value for the given key
        if data_list[i] in data_dict:
            #Add the rest of data to the value string
            data_dict[data_list[i]] = data_dict[data_list[i]] + ' ' + \
                                     data_list[i+1]
        else:
            #Otherwise add value to given key
            data_dict[data_list[i]] = data_list[i+1]
    return data_dict

def extractCommandLineData(command_block):
    """Returns a dict of all command line data from MEME output.
    """
    data_dict = {}
    #Get only necessary Command Line Summary data
    ignore = lambda x: x.startswith('*')
    meme_model = LabeledRecordFinder(lambda x: 'model:' in x,
                                    ignore=ignore)
    cmd_data = list(meme_model(command_block))
    cmd_data = cmd_data[1]
    cmd_data = cmd_data[:-4]

    #Just return list of strings rather than parse data
    """
    cmd_data = '^'.join(cmd_data)
    cmd_data = cmd_data.split()
    cmd_data = ' '.join(cmd_data)
    cmd_data = cmd_data.split(': ')
    lastkarat = DelimitedSplitter('^',-1)
    cmd_data_temp = []
    for line in cmd_data:
        cmd_data_temp.extend(lastkarat(line))
    cmd_data = '>'.join(cmd_data_temp)
    cmd_data = cmd_data.replace('= ','=')
    cmd_data = cmd_data.replace('^',' ')
    cmd_data = cmd_data.split('>')
    """

    return cmd_data

def getModuleDataBlocks(module_blocks):
    """Returns list data blocks for each module.
    """
    #Get blocks of module information for each module
    meme_module_data = LabeledRecordFinder(lambda x: x.startswith('Motif'))
    module_data_blocks = []
    for module in module_blocks:
        module_data_blocks.append(list(meme_module_data(module)))
    return module_data_blocks

def extractModuleData(module_data, alphabet, remap_dict):
    """Creates Module object given module_data list.

        - Only works on 1 module at a time: only pass in data from one module.

    """
    #Create Module object
    meme_module = {}
    
    #Only keep first 3 elements of the list
    module_data = module_data[:3]

    #Get Module general information: module_data[0]
    #Only need to keep first line
    general_dict = getModuleGeneralInfo(module_data[0][0])
    module_length = int(general_dict['width'])

    #Get ModuleInstances: module_data[2]
    instance_data = module_data[2][4:-2]
    for i in xrange(len(instance_data)):
        instance_data[i] = instance_data[i].split()
    #Create a ModuleInstance object and add it to Module for each instance
    for instance in instance_data:
        seqId = remap_dict[instance[0]]
        start = int(instance[1])-1
        Pvalue = float(instance[2])
        sequence = instance[4]
        #Create Location object for ModuleInstance
        location = Location(seqId, start, start + module_length)
        #Create ModuleInstance
        mod_instance = ModuleInstance(sequence,location,Pvalue)

        #Add ModuleInstance to Module
        meme_module[(seqId,start)] = mod_instance
    
    meme_module = Module(meme_module, MolType=alphabet)
    #Get Multilevel Consensus Sequence
    meme_module.ConsensusSequence = getConsensusSequence(module_data[1])
    #Pull out desired values from dict    
    meme_module.Llr = int(general_dict['llr'])
    meme_module.Evalue = float(general_dict['E-value'])
    meme_module.ID = general_dict['MOTIF']

    return meme_module

def getConsensusSequence(first_block):
    """Returns multilevel consensus sequences string.
    """
    for line in first_block:
        if line.upper().startswith('MULTILEVEL'):
            return line.split()[1]

def getModuleGeneralInfo(module_general):
    """Returns dict with Module general information.

        - Module general information includes:
            - width, sites, llr, E-value
    """
    module_id = module_general[:8]
    module_general = module_general[8:].strip().replace(' =','')
    module_general = module_general.split()
    #Get dict of Module general info from list
    general_dict = dictFromList(module_general)
    general_dict['MOTIF']=module_id[5:].strip()
    return general_dict

def extractSummaryData(summary_block):
    """Returns dict of sequences and combined P values.

        - {'CombinedP':{
            'seqId1': Pvalue1,
            'seqId2': Pvalue2,}
            }
    """
    #Get slice of necessary data from summary_block
    summary = summary_block[7:]
    #print summary
    summary_dict = {}
    #Split on whitespace
    for i in xrange(len(summary)):
        summary[i] = summary[i].split()
    #Add necesary data to dict
    for seq in summary:
        #Stop when '--------------------' is found: end of data
        if seq[0].startswith('--------------------'):
            break
        if len(seq) < 3:
            continue
        summary_dict[seq[0]] = float(seq[1])
    return {'CombinedP':summary_dict}

def MemeParser(lines, allowed_ids=[]):
    """Returns a MotifResults object given a MEME results file.
    """
    warnings = []
    #Create MotifResults object
    meme_motif_results = MotifResults()
    #Get main block and alphabet
    main_block, alphabet = getDataBlock(lines)
    #Add alphabet to MotifResults object
    meme_motif_results.MolType = alphabet
    #Get command line summary block and module blocks
    command_block, module_blocks = getCommandModuleBlocks(main_block)
    if command_block:
        #Extract command line data and put in dict
        parameters_list = extractCommandLineData(command_block)
        #Add parameters dict to MotifResults object parameters
        meme_motif_results.Parameters = parameters_list
    #make sure modules were found
    if len(module_blocks) > 0:
        #Get Summary of motifs block
        summary_block = getSummaryBlock(module_blocks[-1])
        #Extract summary data and get summary_dict
        summary_dict = extractSummaryData(summary_block)
        seq_names = summary_dict['CombinedP'].keys()
        if allowed_ids:
            remap_dict,warning = make_remap_dict(seq_names,allowed_ids)
            if warning:
                warnings.append(warning)
            sd = {}
            for k,v in summary_dict['CombinedP'].items():
                sd[remap_dict[k]]=v
            summary_dict['CombinedP']=sd
        else:
            remap_dict = dict(zip(seq_names,seq_names))
        
        #Add summary dict to MotifResults object
        meme_motif_results.Results = summary_dict
        
        #Add warnings to MotifResults object
        meme_motif_results.Results['Warnings']=warnings
        
        #Get blocks for each module
        module_blocks = getModuleDataBlocks(module_blocks)
        #Extract modules and put in MotifResults.Modules list
        for module in module_blocks:
            meme_motif_results.Modules.append(extractModuleData(module,\
                alphabet,remap_dict))
        for module in meme_motif_results.Modules:
            meme_motif_results.Motifs.append(Motif(module))
    return meme_motif_results

