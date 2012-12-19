#!/usr/bin/env python

"""Application controller for cdbfasta

Code obtained from http://sourceforge.net/projects/cdbfasta/

cdbfasta Version 0.99, dated 07-22-10 on download
cdbyank Version 0.981, dated 07-22-10 on download
"""

from os import remove, path
from cogent.app.parameters import FlagParameter, ValuedParameter
from cogent.app.util import CommandLineApplication, ResultPath, \
    get_tmp_filename, guess_input_handler

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Prototype"

class cdbfasta(CommandLineApplication):
    """cdbfasta application controller"""

    _options ={
        # -o the index file will be named <index_file>; if not given,
        # the index filename is database name plus the suffix '.cidx'
        '-o':ValuedParameter('-',Name='anchorspacing',Delimiter=' '),
        
        # -r <record_delimiter> a string of characters at the beginning of line
        # marking the start of a record (default: '>')
        '-clw':FlagParameter(Prefix='-',Name='clw'),

        # -Q treat input as fastq format, i.e. with '@' as record delimiter
        # and with records expected to have at least 4 lines
        '-Q':FlagParameter(Prefix='-',Name='Q'),

        # -z database is compressed into the file <compressed_db>
        # before indexing (<fastafile> can be "-" or "stdin" 
        # in order to get the input records from stdin)
        '-z':ValuedParameter('-',Name='z',Delimiter=' '),

        # -s strip extraneous characters from *around* the space delimited
        # tokens, for the multikey options below (-m,-n,-f);
        # Default <stripendchars> set is: '",`.(){}/[]!:;~|><+-
        '-s':ValuedParameter('-',Name='s',Delimiter=' '),

        # -m ("multi-key" option) create hash entries pointing to 
        # the same record for all tokens found in
        # the defline
        '-m':FlagParameter('-',Name='m'),

        # -n <numkeys> same as -m, but only takes the first <numkeys>
        # tokens from the defline
        '-n':ValuedParameter('-',Name='n',Delimiter=' '),

        # -f indexes *space* delimited tokens (fields) in the defline as given
        # by LIST of fields or fields ranges (the same syntax as UNIX 'cut')
        '-f':ValuedParameter('-',Name='f',Delimiter=''),

        # -w <stopwordslist> exclude from indexing all the words found
        # in the file <stopwordslist> (for options -m, -n and -k)
        '-w':ValuedParameter('-',Name='w',Delimiter=' '),

        # -i do case insensitive indexing (i.e. create additional keys for 
        # all-lowercase tokens used for indexing from the defline 
        '-i':FlagParameter('-',Name='i'),

        # -c for deflines in the format: db1|accession1|db2|accession2|...,
        # only the first db-accession pair ('db1|accession1') is taken as key
        '-c':FlagParameter('-',Name='c'),

        # -C like -c, but also subsequent db|accession constructs are indexed,
        # along with the full (default) token; additionally,
        # all nrdb concatenated accessions found in the defline 
        # are parsed and stored (assuming 0x01 or '^|^' as separators)
        '-C':FlagParameter('-', Name='C'),

        # -a accession mode: like -C option, but indexes the 'accession'
        # part for all 'db|accession' constructs found
        '-a':FlagParameter('-', Name='a'),

        # -A like -a and -C together (both accessions and 'db|accession'
        # constructs are used as keys
        '-A':FlagParameter('-', Name='A'),

        # -v show program version and exit
        '-v':FlagParameter('-', Name='v')
        } 

    _parameters = {}
    _parameters.update(_options)
    _command = "cdbfasta"
    _input_file = ""

    def _input_as_string(self, data):
        """Index a single file"""
        if not data:
            raise ValueError, "Expected a file!"
        if not path.exists(data):
            raise ValueError, "File to index doesn't exist: %s" % data

        self._input_file = data
        return ""

    def _get_result_paths(self,data):
        if self.Parameters['-v'].isOn():
            return {}

        output = {}
        if self.Parameters['-o'].isOn():
            output['cidx'] = ResultPath(self.Parameters['-o'].Value)
        else:
            output['cidx'] = ResultPath(self._input_file + '.cidx')
        return output
    
    def _get_base_command(self):
        """Yay for positional arguments..."""
        command_parts = []
        cd_command = ''.join(['cd ',str(self.WorkingDir),';'])
        if self._command is None:
            raise ApplicationError, '_command has not been set.'
        command = self._command
        parameters = sorted([str(x) for x in self.Parameters.values() 
                            if str(x)])

        synonyms = self._synonyms

        command_parts.append(cd_command)
        command_parts.append(command)
        command_parts.append(self._input_file) # Positional argument
        command_parts += parameters

        return self._command_delimiter.join(filter(None,command_parts)).strip()

    BaseCommand = property(_get_base_command)

class cdbyank(CommandLineApplication):
    """cdbyank application controller"""

    _options ={
        # -a <key> the sequence name (accession) for a fasta record to be
        # retrieved; if not given, a list of accessions is expected
        # at stdin
        '-a':ValuedParameter('-',Name='a', Delimiter=' '),

        # -d <fasta_file> is the fasta file to pull records from; 
        # if not specified, cdbyank will look in the same directory
        # where <index_file> resides, for a file with the same name
        # but without the ".cidx" suffix
        '-d':ValuedParameter('-', Name='d', Delimiter=' '),

        # -o the records found are written to file <outfile> instead of stdout
        '-o':ValuedParameter('-', Name='o', Delimiter=' '),

        # -x allows retrieval of multiple records per key, if the indexed 
        # database had records with the same key (non-unique keys);
        # (without -x only one record for a given key is retrieved)
        '-x':FlagParameter('-', Name='x'),

        # -i case insensitive query (expects the <index_file> to have been 
        # created with cdbfasta -i option)
        '-i':FlagParameter('-', Name='i'),

        # -Q output the query key surrounded by character '%' before the
        # corresponding record
        '-Q':FlagParameter('-', Name='Q'),

        # -q same as -Q but use character <char> instead of '%'
        '-q':ValuedParameter('-', Name='q', Delimiter=' '),

        # -w enable warnings (sent to stderr) when a key is not found
        '-w':FlagParameter('-', Name='w'),

        # -F pulls only the defline for each record (discard the sequence)
        '-F':FlagParameter('-', Name='F'),

        # -P only displays the position(s) (file offset) within the 
        # database file, for the requested record(s)
        '-P':FlagParameter('-', Name='P'),

        # -R sequence range extraction: expects the input <key(s)> to have 
        # the format: '<seq_name> <start> <end>'
        # and pulls only the specified sequence range
        '-R':ValuedParameter('-', Name='R', Delimiter=' '),

        # -z decompress the entire file <dbfasta.cdbz>
        # (assumes it was built using cdbfasta with '-z' option)
        '-z':ValuedParameter('-', Name='z', Delimiter=' '),

        # -v show version number and exit
        '-v':FlagParameter('-', Name='v'),

        ###
        # Index file statistics (no database file needed):
        # -n display the number of records indexed
        '-n':FlagParameter('-', Name='n'),

        # -l list all keys stored in <index_file>
        '-l':FlagParameter('-', Name='l'),

        # -s display indexing summary info
        '-s':FlagParameter('-', Name='s')
    }

    _parameters = {}
    _parameters.update(_options)
    _command = "cdbyank"
    _input_file = ""
    _queries = []

    def _input_as_string(self, data):
        """File path for an index"""
        if not data:
            raise ValueError, "Expected a file!"
        if not path.exists(data):
            raise ValueError, "Index doesn't exist: %s" % data

        self._input_file = data
        return ""
       
    def _get_result_paths(self, data):
        if self.Parameters['-v'].isOn():
            return {}

        output = {}
        if self.Parameters['-o'].isOn():
            output['seqs'] = ResultPath(self.Parameters['-o'].Value)
        return output
   
    def _get_base_command(self):
        """Yay for positional arguments..."""
        command_parts = []
        cd_command = ''.join(['cd ',str(self.WorkingDir),';'])
        if self._command is None:
            raise ApplicationError, '_command has not been set.'
        command = self._command
        parameters = sorted([str(x) for x in self.Parameters.values() 
                            if str(x)])

        synonyms = self._synonyms

        if self._queries:
            bulk_query = 'echo "%s" | ' % " ".join(self._queries)
        else:
            bulk_query = ""

        command_parts.append(cd_command)
        command_parts.append(bulk_query)
        command_parts.append(command)
        command_parts.append(self._input_file) # Positional argument
        command_parts += parameters

        return self._command_delimiter.join(filter(None,command_parts)).strip()

    BaseCommand = property(_get_base_command)

    def setQueries(self, queries):
        """Sets queries"""
        self._queries = queries

def index_fasta(filename, params=None):
    """Index a fasta file, returns the absolute path for the index"""
    if params is None:
        params = {}

    app = cdbfasta(params)

    return app(filename)['cidx'].name

def index_fastq(filename, params=None):
    """Index a fastq file, returns the absolute path for the index"""
    if params is None:
        params = {'-Q':True}

    app = cdbfasta(params)

    return app(filename)['cidx'].name

def query_indexed_seqs(queries, path_to_index, params=None):
    """Returns sequences from the indexed file
    
    works with both fasta and fastq, returns in the file type of the indexed
    """
    if params is None:
        params = {}

    app = cdbyank(params)
    app.setQueries(queries)

    res = app(path_to_index)

    if 'seqs' not in res:
        res['seqs'] = res['StdOut']

    return res['seqs'].read()
