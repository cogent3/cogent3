import sqlalchemy as sql

from cogent import DNA
from cogent.core.location import Map

from cogent.db.ensembl.species import Species
from cogent.db.ensembl.util import NoItemError, asserted_one
from cogent.db.ensembl.assembly import CoordSystem, Coordinate, \
                                    get_coord_conversion

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

# local reference to the sqlalchemy substring function
substr = sql.sql.func.substr

def _assemble_seq(frags, start, end, frag_positions):
    """returns a single string in which missing sequence is replaced by 'N'"""
    prev_end = start
    assert len(frag_positions) == len(frags), "Mismatched number of "\
                                                    "fragments and positions"
    assembled = []
    for index, (frag_start, frag_end) in enumerate(frag_positions):
        diff = frag_start - prev_end
        assert diff >= 0, 'fragment position start < previous end: %s, %s' %\
                                                (frag_start, prev_end)
        assembled += ['N'*diff, frags[index]]
        prev_end = frag_end
    diff = end - frag_end
    assert diff >= 0, 'end[%s] < previous frag_end[%s]' % (end, frag_end)
    assembled += ['N' * diff]
    return DNA.makeSequence(''.join(assembled))

def _make_coord(genome, coord_name, start, end, strand):
    """returns a Coordinate"""
    return Coordinate(CoordName=coord_name, Start=start, End=end,
                Strand=strand, genome = genome)

def get_sequence(coord=None, genome=None, coord_name=None, start=None,
                                        end=None, strand = 1, DEBUG=False):
    # TODO clean up use of a coord
    if coord is None:
        coord = _make_coord(genome, coord_name, start, end, 1)
    else:
        coord = coord.copy()
    
    start, end, strand = coord.Start, coord.End, coord.Strand
    coord_name = coord.CoordName
    
    genome = coord.genome
    # no matter what strand user provide, we get the + sequence first
    coord.Strand = 1
    species = genome.Species
    coord_type = CoordSystem(species=species,core_db=genome.CoreDb,
                             seq_level=True)
    if DEBUG:
        print 'Created Coordinate:',coord,coord.EnsemblStart,coord.EnsemblEnd
        print coord.CoordType, coord_type
    assemblys = get_coord_conversion(coord, coord_type, genome.CoreDb)
    
    if not assemblys:
        raise NoItemError, 'no assembly for %s' % coord
    
    dna = genome.CoreDb.getTable('dna')
    seqs, positions = [], []
    for q_loc, t_loc in assemblys:
        assert q_loc.Strand == 1
        length = len(t_loc)
        # get MySQL to do the string slicing via substr function
        query = sql.select([substr(dna.c.sequence,
                                  t_loc.EnsemblStart,
                                  length).label('sequence')],
                            dna.c.seq_region_id == t_loc.seq_region_id)
        record = asserted_one(query.execute().fetchall())
        seq = record['sequence']
        seq = DNA.makeSequence(seq)
        if t_loc.Strand == -1:
            seq = seq.rc()
        seqs.append(str(seq))
        positions.append((q_loc.Start, q_loc.End))
    sequence = _assemble_seq(seqs, coord.Start, coord.End, positions)
    
    if strand == -1:
        sequence = sequence.rc()
    return sequence

