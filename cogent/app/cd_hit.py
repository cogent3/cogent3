#!/usr/bin/env python
"""Application controller for CD-HIT v3.1.1"""

import shutil
from os import remove
from cogent.app.parameters import ValuedParameter
from cogent.app.util import CommandLineApplication, ResultPath,\
        get_tmp_filename
from cogent.core.moltype import RNA, DNA, PROTEIN
from cogent.core.alignment import SequenceCollection
from cogent.parse.fasta import MinimalFastaParser

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class CD_HIT(CommandLineApplication):
    """cd-hit Application Controller

    Use this version of CD-HIT if your MolType is PROTEIN
    """

    _command = 'cd-hit'
    _input_handler = '_input_as_multiline_string'
    _parameters = {
        # input input filename in fasta format, required
        '-i':ValuedParameter('-',Name='i',Delimiter=' ',IsPath=True),

        # output filename, required
        '-o':ValuedParameter('-',Name='o',Delimiter=' ',IsPath=True),

        # sequence identity threshold, default 0.9
        # this is the default cd-hit's "global sequence identity" calc'd as :
        # number of identical amino acids in alignment
        # divided by the full length of the shorter sequence
        '-c':ValuedParameter('-',Name='c',Delimiter=' '),

        # use global sequence identity, default 1
        # if set to 0, then use local sequence identity, calculated as :
        # number of identical amino acids in alignment
        # divided by the length of the alignment
        # NOTE!!! don't use -G 0 unless you use alignment coverage controls
        # see options -aL, -AL, -aS, -AS
        '-g':ValuedParameter('-',Name='g',Delimiter=' '),

        # band_width of alignment, default 20
        '-b':ValuedParameter('-',Name='b',Delimiter=' '),

        # max available memory (Mbyte), default 400
        '-M':ValuedParameter('-',Name='M',Delimiter=' '),

        # word_length, default 8, see user's guide for choosing it
        '-n':ValuedParameter('-',Name='n',Delimiter=' '),

        # length of throw_away_sequences, default 10
        '-l':ValuedParameter('-',Name='l',Delimiter=' '),

        # tolerance for redundance, default 2
        '-t':ValuedParameter('-',Name='t',Delimiter=' '),

        # length of description in .clstr file, default 20
        # if set to 0, it takes the fasta defline and stops at first space
        '-d':ValuedParameter('-',Name='d',Delimiter=' '),

        # length difference cutoff, default 0.0
        # if set to 0.9, the shorter sequences need to be
        # at least 90% length of the representative of the cluster
        '-s':ValuedParameter('-',Name='s',Delimiter=' '),

        # length difference cutoff in amino acid, default 999999
        # f set to 60, the length difference between the shorter sequences
        # and the representative of the cluster can not be bigger than 60
        '-S':ValuedParameter('-',Name='S',Delimiter=' '),

        # alignment coverage for the longer sequence, default 0.0
        # if set to 0.9, the alignment must covers 90% of the sequence
        '-aL':ValuedParameter('-',Name='aL',Delimiter=' '),

        # alignment coverage control for the longer sequence, default 99999999
        # if set to 60, and the length of the sequence is 400,
        # then the alignment must be >= 340 (400-60) residues
        '-AL':ValuedParameter('-',Name='AL',Delimiter=' '),

        # alignment coverage for the shorter sequence, default 0.0
        # if set to 0.9, the alignment must covers 90% of the sequence
        '-aS':ValuedParameter('-',Name='aS',Delimiter=' '),

        # alignment coverage control for the shorter sequence, default 99999999
        # if set to 60, and the length of the sequence is 400,
        # then the alignment must be >= 340 (400-60) residues
        '-AS':ValuedParameter('-',Name='AS',Delimiter=' '),

        # 1 or 0, default 0, by default, sequences are stored in RAM
        # if set to 1, sequence are stored on hard drive
        # it is recommended to use -B 1 for huge databases
        '-B':ValuedParameter('-',Name='B',Delimiter=' '),

        # 1 or 0, default 0
        # if set to 1, print alignment overlap in .clstr file
        '-p':ValuedParameter('-',Name='p',Delimiter=' '),

        # 1 or 0, default 0
        # by cd-hit's default algorithm, a sequence is clustered to the first 
        # cluster that meet the threshold (fast cluster). If set to 1, the program
        # will cluster it into the most similar cluster that meet the threshold
        # (accurate but slow mode)
        # but either 1 or 0 won't change the representatives of final clusters
        '-g':ValuedParameter('-',Name='g',Delimiter=' '),

        # print this help
        '-h':ValuedParameter('-',Name='h',Delimiter=' ')
    }
    _synonyms = {'Similarity':'-c'}
 
    def getHelp(self):
        """Method that points to documentation"""
        help_str =\
        """
        CD-HIT is hosted as an open source project at:
        http://www.bioinformatics.org/cd-hit/

        The following papers should be cited if this resource is used:

        Clustering of highly homologous sequences to reduce thesize of large 
        protein database", Weizhong Li, Lukasz Jaroszewski & Adam Godzik
        Bioinformatics, (2001) 17:282-283

        Tolerating some redundancy significantly speeds up clustering of large
        protein databases", Weizhong Li, Lukasz Jaroszewski & Adam Godzik 
        Bioinformatics, (2002) 18:77-82
        """
        return help_str

    def _input_as_multiline_string(self, data):
        """Writes data to tempfile and sets -i parameter

        data -- list of lines
        """
        if data:
            self.Parameters['-i']\
                    .on(super(CD_HIT,self)._input_as_multiline_string(data))
        return ''

    def _input_as_lines(self, data):
        """Writes data to tempfile and sets -i parameter

        data -- list of lines, ready to be written to file
        """
        if data:
            self.Parameters['-i']\
                    .on(super(CD_HIT,self)._input_as_lines(data))
        return ''

    def _input_as_seqs(self, data):
        """Creates a list of seqs to pass to _input_as_lines

        data -- list like object of sequences
        """
        lines = []
        for i,s in enumerate(data):
            # will number the sequences 1,2,3, etc...
            lines.append(''.join(['>',str(i+1)]))
            lines.append(s)
        return self._input_as_lines(lines)

    def _input_as_string(self, data):
        """Makes data the value of a specific parameter"""
        if data:
            self.Parameters['-i'].on(str(data))
        return ''

    def _get_seqs_outfile(self):
        """Returns the absolute path to the seqs outfile"""
        if self.Parameters['-o'].isOn():
            return self.Parameters['-o'].Value
        else:
            raise ValueError, "No output file specified"

    def _get_clstr_outfile(self):
        """Returns the absolute path to the clstr outfile"""
        if self.Parameters['-o'].isOn():
            return ''.join([self.Parameters['-o'].Value, '.clstr'])
        else:
            raise ValueError, "No output file specified"

    def _get_result_paths(self, data):
        """Return dict of {key: ResultPath}"""
        result = {}
        result['FASTA'] = ResultPath(Path=self._get_seqs_outfile())
        result['CLSTR'] = ResultPath(Path=self._get_clstr_outfile())
        return result

class CD_HIT_EST(CD_HIT):
    """cd-hit Application Controller

    Use this version of CD-HIT if your MolType is PROTEIN
    """

    _command = 'cd-hit-est'
    _input_handler = '_input_as_multiline_string'
    _parameters = CD_HIT._parameters
    _parameters.update({\
        # 1 or 0, default 0, by default only +/+ strand alignment
        # if set to 1, do both +/+ & +/- alignments
        '-r':ValuedParameter('-',Name='r',Delimiter=' ')
        })

def cdhit_clusters_from_seqs(seqs, moltype, params=None):
    """Returns the CD-HIT clusters given seqs

    seqs        : dict like collection of sequences
    moltype     : cogent.core.moltype object
    params      : cd-hit parameters

    NOTE: This method will call CD_HIT if moltype is PROTIEN,
        CD_HIT_EST if moltype is RNA/DNA, and raise if any other
        moltype is passed.
    """
    # keys are not remapped. Tested against seq_ids of 100char length
    seqs = SequenceCollection(seqs, MolType=moltype)
    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seqs.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)
    
    # setup params and make sure the output argument is set
    if params is None:
        params = {}
    if '-o' not in params:
        params['-o'] = get_tmp_filename()

    # call the correct version of cd-hit base on moltype
    working_dir = get_tmp_filename()
    if moltype is PROTEIN:
        app = CD_HIT(WorkingDir=working_dir, params=params)
    elif moltype is RNA:
        app = CD_HIT_EST(WorkingDir=working_dir, params=params)
    elif moltype is DNA:
        app = CD_HIT_EST(WorkingDir=working_dir, params=params)
    else:
        raise ValueError, "Moltype must be either PROTEIN, RNA, or DNA"

    # grab result
    res = app(int_map.toFasta())
    clusters = parse_cdhit_clstr_file(res['CLSTR'].readlines())

    remapped_clusters = []
    for c in clusters:
        curr = [int_keys[i] for i in c]
        remapped_clusters.append(curr)

    # perform cleanup
    res.cleanUp()
    shutil.rmtree(working_dir)
    remove(params['-o'] + '.bak.clstr')

    return remapped_clusters

def cdhit_from_seqs(seqs, moltype, params=None):
    """Returns the CD-HIT results given seqs

    seqs    : dict like collection of sequences
    moltype : cogent.core.moltype object
    params  : cd-hit parameters

    NOTE: This method will call CD_HIT if moltype is PROTIEN,
        CD_HIT_EST if moltype is RNA/DNA, and raise if any other
        moltype is passed.
    """
    # keys are not remapped. Tested against seq_ids of 100char length
    seqs = SequenceCollection(seqs, MolType=moltype)

    # setup params and make sure the output argument is set
    if params is None:
        params = {}
    if '-o' not in params:
        params['-o'] = get_tmp_filename()

    # call the correct version of cd-hit base on moltype
    working_dir = get_tmp_filename()
    if moltype is PROTEIN:
        app = CD_HIT(WorkingDir=working_dir, params=params)
    elif moltype is RNA:
        app = CD_HIT_EST(WorkingDir=working_dir, params=params)
    elif moltype is DNA:
        app = CD_HIT_EST(WorkingDir=working_dir, params=params)
    else:
        raise ValueError, "Moltype must be either PROTEIN, RNA, or DNA"

    # grab result
    res = app(seqs.toFasta())
    new_seqs = dict(MinimalFastaParser(res['FASTA'].readlines()))

    # perform cleanup
    res.cleanUp()
    shutil.rmtree(working_dir)
    remove(params['-o'] + '.bak.clstr')

    return SequenceCollection(new_seqs, MolType=moltype)

def clean_cluster_seq_id(id):
    """Returns a cleaned cd-hit sequence id

    The cluster file has sequence ids in the form of:
    >some_id...
    """
    return id[1:-3]

def parse_cdhit_clstr_file(lines):
    """Returns a list of list of sequence ids representing clusters"""
    clusters = []
    curr_cluster = []

    for l in lines:
        if l.startswith('>Cluster'):
            if not curr_cluster:
                continue
            clusters.append(curr_cluster)
            curr_cluster = []
        else:
            curr_cluster.append(clean_cluster_seq_id(l.split()[2]))

    if curr_cluster:
        clusters.append(curr_cluster)

    return clusters

