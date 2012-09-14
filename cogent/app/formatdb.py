#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# formatdb.py

""" Description
File created on 16 Sep 2009.

"""
from __future__ import division
from optparse import OptionParser
from os.path import split, splitext
from os import remove
from glob import glob
from cogent.app.util import CommandLineApplication, ResultPath, get_tmp_filename
from cogent.app.parameters import ValuedParameter, FilePath

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Production"

class FormatDb(CommandLineApplication):
    """ ApplicationController for formatting blast databases
    
        Currently contains a minimal parameter set.
    """

    _command = 'formatdb'
    _parameters = {\
     '-i':ValuedParameter(Prefix='-',Name='i',Delimiter=' ',IsPath=True),\
     '-l':ValuedParameter(Prefix='-',Name='l',Delimiter=' ',IsPath=True),\
     '-o':ValuedParameter(Prefix='-',Name='o',Delimiter=' ',Value='T'),\
     '-p':ValuedParameter(Prefix='-',Name='p',Delimiter=' ',Value='F'),\
     '-n':ValuedParameter(Prefix='-',Name='n',Delimiter=' ')
     }
    _input_handler = '_input_as_parameter'
    _suppress_stdout = True
    _suppress_stderr = True

    def _input_as_parameter(self,data):
        """ Set the input path and log path based on data (a fasta filepath)
        """
        self.Parameters['-i'].on(data)
        # access data through self.Parameters so we know it's been cast
        # to a FilePath
        input_filepath = self.Parameters['-i'].Value
        input_file_dir, input_filename = split(input_filepath)
        input_file_base, input_file_ext = splitext(input_filename)
        # FIXME: the following all other options
        # formatdb ignores the working directory if not name is passed.
        self.Parameters['-l'].on(FilePath('%s.log') % input_filename)
        self.Parameters['-n'].on(FilePath(input_filename))
        return ''

    def _get_result_paths(self,data):
        """ Build the dict of result filepaths
        """
        # access data through self.Parameters so we know it's been cast
        # to a FilePath
        wd = self.WorkingDir
        db_name = self.Parameters['-n'].Value
        log_name = self.Parameters['-l'].Value
        result = {}
        result['log'] = ResultPath(Path=wd + log_name, IsWritten=True)
        if self.Parameters['-p'].Value == 'F':
            extensions = ['nhr','nin','nsq','nsd','nsi']
        else:
            extensions = ['phr','pin','psq','psd','psi']
        for extension in extensions:
            for file_path in glob(wd + (db_name + '*' + extension)):
                # this will match e.g. nr.01.psd and nr.psd
                key = file_path.split(db_name + '.')[1]
                result_path = ResultPath(Path=file_path, IsWritten=True)
                result[key] = result_path
        return result

    def _accept_exit_status(self,exit_status):
        """ Return True when the exit status was 0
        """
        return exit_status == 0
        
def build_blast_db_from_fasta_path(fasta_path,is_protein=False,\
    output_dir=None,HALT_EXEC=False):
    """Build blast db from fasta_path; return db name and list of files created
    
        **If using to create temporary blast databases, you can call
        cogent.util.misc.remove_files(db_filepaths) to clean up all the
        files created by formatdb when you're done with the database.
    
        fasta_path: path to fasta file of sequences to build database from
        is_protein: True if working on protein seqs (default: False)
        output_dir: directory where output should be written
         (default: directory containing fasta_path)
        HALT_EXEC: halt just before running the formatdb command and
         print the command -- useful for debugging
    """
    fasta_dir, fasta_filename = split(fasta_path)
    if not output_dir:
        output_dir = fasta_dir or '.'
        # Will cd to this directory, so just pass the filename
        # so the app is not confused by relative paths
        fasta_path = fasta_filename
        
    if not output_dir.endswith('/'):
        db_name = output_dir + '/' + fasta_filename
    else:
        db_name = output_dir + fasta_filename

    # instantiate the object
    fdb = FormatDb(WorkingDir=output_dir,HALT_EXEC=HALT_EXEC)
    if is_protein:
        fdb.Parameters['-p'].on('T')
    else:
        fdb.Parameters['-p'].on('F')
    app_result = fdb(fasta_path)
    db_filepaths = []
    for v in app_result.values():
        try:
            db_filepaths.append(v.name)
        except AttributeError:
            # not a file object, so no path to return
            pass
    return db_name, db_filepaths
    
def build_blast_db_from_fasta_file(fasta_file,is_protein=False,\
    output_dir=None,HALT_EXEC=False):
    """Build blast db from fasta_path; return db name and list of files created
    
        **If using to create temporary blast databases, you can call
        cogent.util.misc.remove_files(db_filepaths) to clean up all the
        files created by formatdb when you're done with the database.
    
        fasta_path: path to fasta file of sequences to build database from
        is_protein: True if working on protein seqs (default: False)
        output_dir: directory where output should be written
         (default: directory containing fasta_path)
        HALT_EXEC: halt just before running the formatdb command and
         print the command -- useful for debugging
    """
    output_dir = output_dir or '.'
    fasta_path = get_tmp_filename(\
     tmp_dir=output_dir, prefix="BLAST_temp_db_", suffix=".fasta")
    
    fasta_f = open(fasta_path,'w')
    for line in fasta_file:
        fasta_f.write('%s\n' % line.strip())
    fasta_f.close()
    
    blast_db, db_filepaths = build_blast_db_from_fasta_path(\
     fasta_path, is_protein=is_protein, output_dir=None, HALT_EXEC=HALT_EXEC)
     
    db_filepaths.append(fasta_path)
    
    return blast_db, db_filepaths
    
def build_blast_db_from_seqs(seqs,is_protein=False,\
    output_dir='./',HALT_EXEC=False):
    """Build blast db from seqs; return db name and list of files created
    
        **If using to create temporary blast databases, you can call
        cogent.util.misc.remove_files(db_filepaths) to clean up all the
        files created by formatdb when you're done with the database.
    
        seqs: sequence collection or alignment object
        is_protein: True if working on protein seqs (default: False)
        output_dir: directory where output should be written
         (default: current directory)
        HALT_EXEC: halt just before running the formatdb command and
         print the command -- useful for debugging
    """
    
    # Build a temp filepath
    tmp_fasta_filepath = get_tmp_filename(\
     prefix='Blast_tmp_db',suffix='.fasta')
    # open the temp file
    tmp_fasta_file = open(tmp_fasta_filepath,'w')
    # write the sequence collection to file
    tmp_fasta_file.write(seqs.toFasta())
    tmp_fasta_file.close()
    
    # build the bast database
    db_name, db_filepaths = build_blast_db_from_fasta_path(\
     tmp_fasta_filepath,is_protein=is_protein,\
     output_dir=output_dir,HALT_EXEC=HALT_EXEC)
     
    # clean-up the temporary file
    remove(tmp_fasta_filepath)
    
    # return the results
    return db_name, db_filepaths


def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = 'usage: %prog [options] fasta_filepath' 
    version = 'Version: %prog 0.1'
    parser = OptionParser(usage=usage, version=version)

    # A binary 'verbose' flag
    parser.add_option('-p','--is_protein',action='store_true',\
        dest='is_protein',default=False,\
        help='Pass if building db of protein sequences '+\
        '[default: False, nucleotide db]')

    parser.add_option('-o','--output_dir',action='store',\
          type='string',dest='output_dir',default=None,
          help='the output directory '+\
          '[default: directory containing input fasta_filepath]')

    opts,args = parser.parse_args()
    num_args = 1
    if len(args) != num_args:
       parser.error('Must provide single filepath to build database from.')

    return opts,args


if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    
    fasta_filepath = args[0]
    is_protein = opts.is_protein
    output_dir = opts.output_dir
    
    db_name, db_filepaths = build_blast_db_from_fasta_path(\
     fasta_filepath,is_protein=is_protein,output_dir=output_dir)
    
    
