import sqlalchemy as sql

from cogent import DNA
from cogent.core.location import Map

from cogent.db.ensembl.species import Species
from cogent.db.ensembl.util import NoItemError, asserted_one
from cogent.db.ensembl.assembly import CoordSystem, Coordinate, \
                                    get_coord_conversion

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
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
                Strand=strand, genome=genome)


def get_lower_coord_conversion(coord, species, core_db):
    coord_system = CoordSystem(species=species, core_db=core_db)
    seq_level_coord_type = CoordSystem(species=species,core_db=core_db,
                             seq_level=True)
    query_rank = coord_system[coord.CoordType].rank
    seq_level_rank = coord_system[seq_level_coord_type].rank
    assemblies = None
    for rank in range(query_rank+1, seq_level_rank):
        coord_type = None
        for key in coord_system.keys():
            if coord_system[key].rank == rank:
                coord_type = coord_system[key].name
                break
        
        if coord_type is None:
            continue
        
        assemblies = get_coord_conversion(coord, coord_type, core_db)
        
        if assemblies: 
            break
        
    
    return assemblies

def _get_sequence_from_direct_assembly(coord=None, DEBUG=False):
    # TODO clean up use of a coord
    genome = coord.genome
    # no matter what strand user provide, we get the + sequence first
    coord.Strand = 1
    species = genome.Species
    coord_type = CoordSystem(species=species,core_db=genome.CoreDb,
                             seq_level=True)
    
    if DEBUG:
        print 'Created Coordinate:',coord,coord.EnsemblStart,coord.EnsemblEnd
        print coord.CoordType, coord_type
    
    assemblies = get_coord_conversion(coord, coord_type, genome.CoreDb)
    
    if not assemblies:
        raise NoItemError, 'no assembly for %s' % coord
    
    dna = genome.CoreDb.getTable('dna')
    seqs, positions = [], []
    for q_loc, t_loc in assemblies:
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
    return sequence

def _get_sequence_from_lower_assembly(coord, DEBUG):
    assemblies = get_lower_coord_conversion(coord, coord.genome.Species,
                                            coord.genome.CoreDb)
    if not assemblies:
        raise NoItemError, 'no assembly for %s' % coord
    
    if DEBUG:
        print '\nMedium_level_assemblies = ', assemblies
    
    seqs, positions = [], []
    for q_loc, t_loc in assemblies:
        t_strand = t_loc.Strand
        temp_seq = _get_sequence_from_direct_assembly(t_loc, DEBUG)
        if t_strand == -1:
            temp_seq = temp_seq.rc()
        
        if DEBUG:
            print q_loc
            print t_loc
            print 'temp_seq = ', temp_seq[:10], '\n'
        
        seqs.append(str(temp_seq))
        positions.append((q_loc.Start, q_loc.End))
    
    sequence = _assemble_seq(seqs, coord.Start, coord.End, positions)
    return sequence

def get_sequence(coord=None, genome=None, coord_name=None, start=None, end=None, strand=1, DEBUG=False):
    if coord is None:
        coord = _make_coord(genome, coord_name, start, end, 1)
    else:
        coord = coord.copy()
    
    strand = coord.Strand
    
    try: 
        sequence = _get_sequence_from_direct_assembly(coord, DEBUG)
    except NoItemError:
        ## means there is no assembly, so we do a thorough assembly by converting according to the "rank"
        sequence = _get_sequence_from_lower_assembly(coord, DEBUG)
    
    if strand == -1:
        sequence = sequence.rc()
    return sequence

