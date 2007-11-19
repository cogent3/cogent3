"""
Application controller for lagan20 applications.

Presently restricted to mlagan
"""
import os, re

from cogent.app.parameters import FlagParameter, ValuedParameter, FilePath
from cogent.app.util import CommandLineApplication, ResultPath, \
                        get_tmp_filename, ApplicationError, guess_input_handler
from cogent.core.moltype import ASCII
from cogent.parse.fasta import MinimalFastaParser

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

if 'LAGAN_DIR' not in os.environ:
    raise RuntimeError, \
        'mlagan cannot run if the LAGAN_DIR environment variable is not set.'

suffix = re.compile("[.][^.]+$")

class Mlagan(CommandLineApplication):
    """Mlagan application controller."""
    
    _options ={

        # -tree "string" 
        # You need to specify a phylogenetic tree for the sequences. This must be a pairwise tree,
        # with parenthesis specifying nodes. Here are a few examples:
        # "(human (mouse rat))"
        # "((human mouse)(fugu zebrafish))"
        # The name of each sequence must be specified somewhere on the fasta line of the input sequence:
        # >g324325|Homo sapiens    human
        # ACTGG....
        # Either "Homo" or "sapiens" or "human" are valid names to call the sequence.
        '-tree':ValuedParameter('-',Name='tree',Delimiter=' '),
        # -translate [default off] 
        # Use translated anchoring (homology done on the amino acid level). This is useful
        # for distant (human/chicken, human/fish, and the like) comparisons.
        '-translate':FlagParameter(Prefix='-',Name='translate'),
        # -fastreject [default off]
        # Abandon the alignment if the homology looks weak. Currently tuned for 
        # human/mouse distance, or closer. Please contact the authors for more 
        # details on this option.
        '-fastreject':FlagParameter(Prefix='-',Name='fastreject'),
        # -out filename [default standard out]
        # Output the alignment to filename, rather than standard out.
        '-out':ValuedParameter('-',Name='out',Delimiter=' ', IsPath=True)
    
    }
    
    _parameters = {}
    _parameters.update(_options)
    _command = "mlagan"
    _input_handler = "_input_as_aln"
    
    def __init__(self, **kwargs):
        CommandLineApplication.__init__(self, **kwargs)
        self._input_fasta_names = []
    
    def _input_as_seqs(self, data):
        # need to generate temp files for each sequence in aln,
        # these will be cleaned up after run
        for seq in data:
            fn = self.getTmpFilename(self.WorkingDir)
            self._input_fasta_names.append(FilePath(fn))
            f = open(fn, "w")
            f.writelines(seq.toFasta())
            f.close()
        return ""
    
    def _input_as_lines(self, lines):
        parser = MinimalFastaParser(lines)
        make_seq = ASCII.makeSequence
        data=[make_seq(seq, label) for label, seq in parser]
        return self._input_as_seqs(data)
    
    def _align_out_filename(self):
        if self.Parameters['-out'].isOn():
            aln_filename = self._absolute(str(self.Parameters['-out'].Value))
        else:
            raise ValueError, "No output file specified."
        return aln_filename
    
    def _input_as_aln(self, aln):
        """input as a cogent alignment instance, writes individual sequences to
        separate tmp files"""
        self._input_as_seqs(aln.Seqs)
        return ""
    
    def _get_Base_Command(self):
        bc = self._get_base_command()
        bc = bc.split(";")
        cd = bc[0]
        ml = bc[1].strip()
        
        ml = ml.split(self._command_delimiter)
        for fn in self._input_fasta_names:
            ml.insert(1, fn)
        ml = self._command_delimiter.join(ml)
        command = "; ".join([cd, ml])
        return command
    
    BaseCommand = property(_get_Base_Command)
    
    def _get_result_paths(self, data):
        """mlagan requires input all sequences in separate files, writes
        pairwise .anchors files, etc.."""
        results = {}
        
        if self.Parameters['-out'].isOn():
            out_name = self._align_out_filename()
            result['Align'] = ResultPath(Path=out_name,IsWritten=True)
        
        file_paths = self._input_fasta_names[:]
        for index, fn in enumerate(file_paths):
            results["SeqIn-%d" % index] = ResultPath(fn)
        
        for i in range(len(file_paths)-1):
            fn_i = suffix.sub("", os.path.split(file_paths[i])[1])
            for j in range(i+1, len(file_paths)):
                # strip trailing suffix
                fn_j = suffix.sub("", os.path.split(file_paths[j])[1])
                fn = "%s%s.anchors" % (fn_j, fn_i)
                fn = os.path.join(self.WorkingDir, fn)
                results["Anchors-%d-%d" % (i,j)] = ResultPath(fn,IsWritten=True)
        
        return results
        
