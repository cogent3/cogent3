import sqlalchemy as sql

from cogent.core.location import Map
from cogent.db.ensembl.species import Species as _Species
from cogent.db.ensembl.util import asserted_one, convert_strand, DisplayString
from cogent.db.ensembl.host import DbConnection

__author__ = "Hua Ying"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Hua Ying"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Hua Ying"
__email__ = "Hua.Ying@anu.edu.au"
__status__ = "alpha"

def location_query(table, query_start, query_end,
    start_col = 'seq_region_start', end_col = 'seq_region_end', query = None,
    where = 'overlap'):
    # TODO should we allow for spans, overlaps, within?
    # the union result is a complex query and has be appended to any other queries
    # in which it's being employed
    # should we be setting default values here regarding the columns that start/end
    # are pulled from, or explicitly state which columns
    if query is None:
        query = sql.select([table])
    
    if where == 'within':
        query.append_whereclause(sql.and_(table.c[start_col] < query_start,
                                         table.c[end_col] > query_end))
    else:
        query.append_whereclause(
            sql.or_(sql.and_(table.c[start_col] < query_start,
                                        table.c[end_col] > query_end),
                                sql.and_(table.c[start_col] >= query_start,
                                        table.c[start_col] <= query_end),
                                sql.and_(table.c[end_col] >= query_start,
                                        table.c[end_col] <= query_end)))
    # the union is only being used here to order the results
    # that usage imposes the limitation this function must be appended to
    # other queries components being built into a fuller SQL query
    # makes me think it shouldn't be here?
    query = query.order_by(table.c[start_col])
    return query

def _get_coord_type_and_seq_region_id(coord_name, core_db):
    seq_region_table = core_db.getTable('seq_region')
    rows = sql.select([seq_region_table]).\
        where(seq_region_table.c.name == coord_name).execute().fetchall()
    species_coord_sys = CoordSystem(species=core_db.db_name.Species,
                                    core_db = core_db)
    try:
        selected_row = asserted_one(rows)
    except ValueError:
        selected_row = None
        for row in rows:
            # not a default_version
            if not row['coord_system_id'] in species_coord_sys:
                continue
            elif not selected_row:
                selected_row = row
                break
        if selected_row is None:
            raise ValueError("Ambigous coordinate name: %s" % coord_name)
    coord_type = species_coord_sys[selected_row['coord_system_id']].name
    return selected_row, coord_type

class Coordinate(object):
    def __init__(self, genome, CoordName, Start, End, Strand = 1,
            CoordType = None, seq_region_id = None, ensembl_coord=False):
        if not CoordType or not (seq_region_id or Start or End):
            seq_region_data, CoordType = \
                _get_coord_type_and_seq_region_id(CoordName, genome.CoreDb)
            seq_region_id = seq_region_data['seq_region_id']
            Start = Start or 0
            End = End or seq_region_data['length']
        # TODO allow creation with just seq_region_id
        self.Species = genome.Species
        self.CoordType = DisplayString(CoordType, repr_length=4,
                                       with_quotes=False)
        self.CoordName = DisplayString(CoordName, repr_length=4,
                                       with_quotes=False)
        # if Start == End, we +1 to End, unless these are ensembl_coord's
        if ensembl_coord:
            Start -= 1
        elif Start == End:
            End += 1
        
        if Start > End:
            assert Strand == -1,\
                    "strand incorrect for start[%s] > end[%s]" % (Start, End)
            Start, End = End, Start
        
        self.Start = Start
        self.End = End
        self.Strand = convert_strand(Strand)
        self.seq_region_id = seq_region_id
        self.genome = genome
    
    def __len__(self):
        return self.End - self.Start
    
    def __cmp__(self, other):
        return cmp((self.CoordName,self.Start), (other.CoordName,other.Start))
    
    def _get_ensembl_start(self):
        # ensembl counting starts from 1
        return self.Start + 1
    
    EnsemblStart = property(_get_ensembl_start)
    
    def _get_ensembl_end(self):
        return self.End
    
    EnsemblEnd = property(_get_ensembl_end)
    
    def __str__(self):
        return '%s:%s:%s:%d-%d:%d' % (self.Species, self.CoordType,
                    self.CoordName, self.Start, self.End, self.Strand)
    
    def __repr__(self):
        my_type = self.__class__.__name__
        name = _Species.getCommonName(self.Species)
        coord_type = self.CoordType
        c = '%s(%r,%r,%r,%d-%d,%d)'%(my_type, name, coord_type,
                    self.CoordName, self.Start, self.End, self.Strand)
        return c.replace("'", "")
    
    def adopted(self, other, shift=False):
        """adopts the seq_region_id (including CoordName and CoordType) of
        another coordinate.
        
        Arguments:
            - shift: an int or True/False. If int, it's added to Start/End.
              If bool, other.Start is added to Start/End"""
        if type(shift) == bool:
            shift = [0, other.Start][shift]
        return self.__class__(other.genome, CoordName=other.CoordName,
                            Start=self.Start+shift, End=self.End+shift,
                            Strand=other.Strand,
                            seq_region_id=other.seq_region_id)
    
    def shifted(self, value):
        """adds value to Start/End coords, returning a new instance."""
        new = self.copy()
        new.Start += value
        new.End += value
        assert len(new) > 0, 'shift generated a negative length'
        return new
    
    def copy(self):
        """returns a copy"""
        return self.__class__(genome=self.genome, CoordName=self.CoordName,
            Start=self.Start, End=self.End, Strand = self.Strand,
            CoordType = self.CoordType, seq_region_id = self.seq_region_id)
    
    def resized(self, from_start, from_end):
        """returns a new resized Coordinate with the
        Start=self.Start+from_start and End = self.End+from_end.
        
        If you want to shift Start upstream, add a -ve number"""
        new = self.copy()
        new.Start += from_start
        new.End += from_end
        try:
            assert len(new) >= 0, 'resized generated a negative length: %s' % new
        except (ValueError, AssertionError):
            raise ValueError
        return new
    
    def makeRelativeTo(self, other, make_relative=True):
        """returns a new coordinate with attributes adopted from other, and
        positioned relative to other."""
        
        if other.Strand != self.Strand:
            Start = other.End-self.End
        elif make_relative:
            Start = self.Start-other.Start
        else:
            Start = self.Start+other.Start
        
        End = Start+len(self)
        
        return self.__class__(other.genome, CoordName=other.CoordName,
                            Start=Start, End=End, Strand=other.Strand,
                            seq_region_id=other.seq_region_id)


class _CoordRecord(object):
    """store one record of the coord"""
    def __init__(self, attrib, rank, name = None, coord_system_id=None):
        self.coord_system_id = coord_system_id
        self.name = name
        self.rank = rank
        self.attr = attrib
    
    def __str__(self):
        return "coord_system_id = %d; name = %s; rank = %d; attr = %s "\
                % (self.coord_system_id, self.name, self.rank, self.attr)


class CoordSystemCache(object):
    """store coord_system table from core database.
    (only read default_version as stated in attrib column)
    There are two ways to get information about coordinate system:
    (1) use coord_type (e.g contig) which is at coord_system.c.name,
        and which are keys of _species_coord_systems[species]
    (2) use coord_system_id (e.g 17 refers to chromosome) which are also keys
        of _species_coord_systems[species]
    (3) to get which level of system is used for storing dna table, check
        'attrib' column of coord_system as default_version, sequence_level.
    """
    # Problem: multiple species (for compara) --> organized as {species: coordsystem}
    # TODO: simplify _species_coord_systems?
    # we place each species coord-system in _species_coord_systems, once, so
    # this attribute is a very _public_ attribute, and serves as a cache to
    # reduce unecessary lookups
    _species_coord_systems = {}
    columns = ['coord_system_id', 'name', 'rank', 'attrib'] # columns needed from coord_system table
    # the attrib property has sequence_level, which means this the coordinate system employed for sequence
    def _set_species_system(self, core_db, species):
        if species in self._species_coord_systems:
            return
        self._species_coord_systems[species] = {}
        coord_table = core_db.getTable('coord_system')
        records = sql.select([coord_table]).where(coord_table.c.attrib.like('default%')).\
                        execute().fetchall()    # only select default version
        for record in records:
            attr = self._species_coord_systems[species]
            for key in ['coord_system_id', 'name']:
                key_val = record[key]
                vals = {}
                for column in self.columns:
                    val = record[column]
                    if isinstance(val, set): # join items in set to one string
                        try:
                            val = ", ".join(val)
                        except TypeError:
                            pass
                    vals[column] = val
                attr[key_val] = _CoordRecord(**vals)
    
    def _get_seq_level_system(self, species):
        """returns the sequence level system for species"""
        sp_sys = self._species_coord_systems[species]
        for key, val in sp_sys.items():
            if 'sequence_level' in val.attr:
                return val.name
        
        raise RuntimeError, 'no coord system for %s' % species
    
    def __call__(self, coord_type = None, core_db = None, species = None,
                 seq_level=False):
        """coord_type can be coord_type or coord_system_id"""
        # TODO should only pass in core_db here, not that and Species, or just
        # the genome - what if someone wants to compare different ensembl
        # releases? keying by species is then a bad idea! better to key by
        # id(object)
        # change identifier to coord_system, handle either string val or int
        # (see MySQL table) as is this shouldn't be a __call__, see line 168
        # for reason why we should have a method to set data: setSpeciesCoord
        # call then then just returns the coords for the named species
        species = _Species.getSpeciesName(species or core_db.db_name.Species)
        self._set_species_system(core_db, species)
        if seq_level:
            result = self._get_seq_level_system(species)
        elif coord_type:
            result = self._species_coord_systems[species][coord_type]
        else:
            result = self._species_coord_systems[species]
        return result

CoordSystem = CoordSystemCache()

def _rank_checking(query_coord_type, target_coord_type, core_db, species):
    # assiting in constructingthe query language for assembly
    
    # in order to convert between coordinate systems, we need to establish the
    # ranking for coordinate types
    # rank defines the order of conversion between coord system 'levels'
    # chromosome has rank 1
    # super contig has rank 2
    # contig has rank 4
    # clone has rank 3
    # converting requires changing columns between 'asm' and 'cmp'
    # converting from clone -> contig, use 'asm' column
    # converting from contig -> clone, use 'cmp' column
    query_rank = CoordSystem(core_db = core_db, species = species,
                             coord_type=query_coord_type).rank
    target_rank = CoordSystem(core_db = core_db, species = species,
                             coord_type=target_coord_type).rank
    
    if query_rank < target_rank:
        query_prefix, target_prefix = 'asm', 'cmp'
    elif query_rank > target_rank:
        query_prefix, target_prefix = 'cmp', 'asm'
    else:
        query_prefix, target_prefix = '', ''
    return query_prefix, target_prefix

def _get_equivalent_coords(query_coord, assembly_row, query_prefix,
                target_prefix, target_coord_type):
    # TODO better function name
    start = query_coord.EnsemblStart
    end = query_coord.EnsemblEnd
    strand = query_coord.Strand
    
    ori = assembly_row['ori']
    q_strand, t_strand = strand, strand * ori
    if 'seq_region' not in query_prefix:
        q_seq_region_id = assembly_row['%s_seq_region_id' % query_prefix]
        t_seq_region_id = assembly_row['%s_seq_region_id' % target_prefix]
    else:
        q_seq_region_id = assembly_row['_'.join([query_prefix, 'id'])]
        t_seq_region_id = assembly_row['_'.join([target_prefix, 'id'])]
    
    # d -- distance
    d_start = max(0, start - int(assembly_row['%s_start' % query_prefix]))
    d_end = max(0, int(assembly_row['%s_end' % query_prefix]) - end)
    # q -- query (to differ from the origin query block)
    q_start = int(assembly_row['%s_start'%query_prefix]) + d_start
    q_end = int(assembly_row['%s_end'%query_prefix]) - d_end
    
    if int(assembly_row['ori']) == -1:
        d_start, d_end = d_end, d_start
    # t -- target
    t_start = int(assembly_row['%s_start' % target_prefix]) + d_start
    t_end = int(assembly_row['%s_end' % target_prefix]) - d_end
    
    q_location = Coordinate(CoordName=query_coord.CoordName, Start=q_start,
                    End=q_end, Strand=q_strand,
                    CoordType=query_coord.CoordType,
                    seq_region_id=q_seq_region_id,
                    genome = query_coord.genome, ensembl_coord=True)
    t_location = Coordinate(CoordName=assembly_row['name'], Start=t_start,
                    End=t_end, Strand=t_strand, CoordType=target_coord_type,
                    seq_region_id=t_seq_region_id,
                    genome = query_coord.genome,
                    ensembl_coord=True)
    return [q_location, t_location]

def assembly_exception_coordinate(loc):
    """returns a coordinate conversion for one with an assembly exception"""
    genome = loc.genome
    assemb_except_table = genome.CoreDb.getTable('assembly_exception')
    seq_region_table = genome.CoreDb.getTable('seq_region')
    
    query = sql.select([assemb_except_table, seq_region_table.c.name],
                sql.and_(
                assemb_except_table.c.seq_region_id == \
                                            loc.seq_region_id,
                assemb_except_table.c.exc_seq_region_id == \
                                            seq_region_table.c.seq_region_id))
    query = location_query(assemb_except_table,
                    loc.Start, loc.End, query = query)
    record = asserted_one(query.execute().fetchall())
    s, conv_loc = _get_equivalent_coords(loc, record, "seq_region",
                    "exc_seq_region", loc.CoordType)
    return conv_loc

def get_coord_conversion(query_location,target_coord_type,core_db,where=None):
    """returns the ???"""
    where = where or 'overlap'
    # TODO better function name
    species = core_db.db_name.Species
    assert query_location.Species == species
    assembly = core_db.getTable('assembly')
    seq_region = core_db.getTable('seq_region')
    target_coord_system_id = CoordSystem(target_coord_type, core_db=core_db,
                                        species=species).coord_system_id
    
    query_prefix, target_prefix = _rank_checking(query_location.CoordType,
                                        target_coord_type, core_db, species)
    if query_prefix == target_prefix:
        return [[query_location, query_location]]
    # TODO: deal with query_prefix == target_prefix == '' --> could happen
    # when query features.
    query = sql.select([assembly, seq_region.c.name], sql.and_(assembly.c\
            ['%s_seq_region_id' % target_prefix] == seq_region.c.seq_region_id,
            seq_region.c.coord_system_id == target_coord_system_id,
            assembly.c['%s_seq_region_id' % query_prefix] ==\
             query_location.seq_region_id))
    query = location_query(assembly, query_location.EnsemblStart,
                           query_location.EnsemblEnd,
                           start_col = "%s_start" % query_prefix,
                           end_col = "%s_end" % query_prefix, query = query,
                           where=where)
    assembly_rows = query.execute().fetchall()
    results = []
    for assembly_row in assembly_rows:
        results.append(_get_equivalent_coords(query_location, assembly_row,
            query_prefix, target_prefix, target_coord_type))
    return results

