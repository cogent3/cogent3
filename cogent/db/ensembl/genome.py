import re
import sqlalchemy as sql

from cogent.db.ensembl.species import Species as _Species
from cogent.db.ensembl.util import LazyRecord, asserted_one,\
            convert_strand, DisplayString
from cogent.db.ensembl.host import get_ensembl_account, get_latest_release
from cogent.db.ensembl.database import Database
from cogent.db.ensembl.assembly import CoordSystem, Coordinate, \
                                        get_coord_conversion, location_query
from cogent.db.ensembl.region import Gene, Variation, GenericRegion, \
                                    CpGisland, Repeat, Est
from cogent.db.ensembl.feature_level import FeatureCoordLevels
from cogent.util.misc import flatten


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

class FeatureTypeCache(LazyRecord):
    """stores critical indices for different feature types"""
    def __init__(self, genome):
        super(FeatureTypeCache, self).__init__()
        self.genome = genome
        self._type_func_map = dict(CpGisland=self._get_cpg_island_analysis_id,
                                    Repeat=self._get_repeat_id)
    
    def _get_cpg_island_analysis_id(self):
        analysis_description_table = \
                             self.genome.CoreDb.getTable('analysis_description')
        query = sql.select([analysis_description_table.c.analysis_id],
                    analysis_description_table.c.display_label.like('%CpG%'))
        record = asserted_one(query.execute())
        self._table_rows['analysis_description'] = record
        quoted_limited = lambda x : DisplayString(x, with_quotes=True, 
                                                  num_words=2)
        self._populate_cache_from_record(
                    [('CpGisland','analysis_id',quoted_limited)],
                    'analysis_description')
    
    def _get_cpg_island_id(self):
        return self._get_cached_value('CpGisland',
                                      self._get_cpg_island_analysis_id)
    
    CpGisland = property(_get_cpg_island_id)
    
    def _get_repeat_id(self):
        raise NotImplementedError
    
    Repeat = property(_get_repeat_id)
    
    def get(self, feature_type):
        """returns the analysis_id for feature_type"""
        try:
            func = self._type_func_map[feature_type]
        except KeyError:
            raise RuntimeError,"Unknown feature type: %s" % feature_type
        return self._get_cached_value(feature_type, func)
    

class Genome(object):
    """An Ensembl Genome"""
    def __init__(self, Species, Release, account=None, pool_recycle=None):
        super(Genome, self).__init__()
        
        assert Release, 'invalid release specified'
        if account is None:
            account = get_ensembl_account(release=Release)
        
        self._account = account
        self._pool_recycle = pool_recycle
        
        # TODO: check Release may not be necessary because: assert Release above
        if Release is None:
            Release = get_latest_release(account=account)
        
        self._gen_release = None
        
        # TODO make name and release immutable properties
        self.Species = _Species.getSpeciesName(Species)
        self.Release = str(Release)
        
        # the db connections
        self._core_db = None
        self._var_db = None
        self._other_db = None
        self._feature_type_ids = FeatureTypeCache(self)
        self._feature_coord_levels = FeatureCoordLevels(self.Species)
    
    def __str__(self):
        my_type = self.__class__.__name__
        return "%s(Species='%s'; Release='%s')" % (my_type, self.Species,
                                                       self.Release)
    
    def __repr__(self):
        return self.__str__()
    
    def __cmp__(self, other):
        return cmp(self.CoreDb, other.CoreDb)
    
    def _connect_db(self, db_type):
        connection = dict(account=self._account, release=self.Release,
                        species=self.Species, pool_recycle=self._pool_recycle)
        if self._core_db is None and db_type == 'core':
            self._core_db = Database(db_type='core', **connection)
            gen_rel = self.CoreDb.db_name.GeneralRelease
            gen_rel = int(re.findall(r'^\d+', str(gen_rel))[0])
            self._gen_release = gen_rel
        elif self._var_db is None and db_type == 'variation':
            self._var_db = Database(db_type='variation', **connection)
        elif self._other_db is None and db_type == 'otherfeatures':
            self._other_db = Database(db_type='otherfeatures', **connection)
    
    def _get_core_db(self):
        self._connect_db('core')
        return self._core_db
    
    CoreDb = property(_get_core_db)
    
    def _get_var_db(self):
        self._connect_db('variation')
        return self._var_db
    
    VarDb = property(_get_var_db)
    
    def _get_other_db(self):
        self._connect_db('otherfeatures')
        return self._other_db
    
    OtherFeaturesDb = property(_get_other_db)
    
    @property
    def GeneralRelease(self):
        """returns True if the general Ensembl release is >= 65"""
        # General release is used here as to support Ensembl genomes
        if self._gen_release is None:
            self.CoreDb
        
        return self._gen_release
    
    def _get_biotype_description_condition(self, gene_table, Description=None, BioType=None, like=True):
        assert Description or BioType, "no valid argument provided"
        btype, descr = None, None
        
        if BioType:
            if like:
                btype = gene_table.c.biotype.like('%'+BioType+'%')
            else:
                btype = gene_table.c.biotype==BioType
        if Description:
            if like:
                descr = gene_table.c.description.like('%'+Description+'%')
            else:
                descr = gene_table.c.description.op('regexp')(
                                    '[[:<:]]%s[[:>:]]' % Description)
        
        if btype is not None and descr is not None:
            condition = sql.and_(btype, descr)
        elif btype is not None:
            condition = btype
        elif descr is not None:
            condition = descr
        
        return condition
    
    def _build_gene_query(self, db, condition, gene_table, gene_id_table, xref_table=None):
        if gene_id_table is None: # Ensembl releases later than >= 65
            join_obj = gene_table
            select_obj = [gene_table]
        else:
            join_obj = gene_id_table.join(gene_table, gene_id_table.c.gene_id==gene_table.c.gene_id)
            select_obj = [gene_id_table.c.stable_id, gene_table]
        
        if db.Type == 'core':
            join_obj = join_obj.outerjoin(xref_table, gene_table.c.display_xref_id==xref_table.c.xref_id)
            select_obj.append(xref_table.c.display_label)
        query = sql.select(select_obj, from_obj=[join_obj], whereclause=condition)
        return query
    
    def _get_symbol_from_synonym(self, db, synonym):
        """returns the gene symbol for a synonym"""
        synonym_table = db.getTable('external_synonym')
        xref_table = db.getTable('xref')
        joinclause = xref_table.join(synonym_table,
                        xref_table.c.xref_id==synonym_table.c.xref_id)
        whereclause = synonym_table.c.synonym==synonym
        query = sql.select([xref_table.c.display_label], from_obj=[joinclause],
            whereclause=whereclause).distinct()
        result = query.execute().fetchall()
        if result:
            try:
                symbol = flatten(result)[0]
            except IndexError:
                symbol = None
        else:
            symbol = None
        return symbol
    
    def _get_gene_query(self, db, Symbol=None, Description=None, StableId=None,
                         BioType=None, synonym=None, like=True):
        xref_table = [None, db.getTable('xref')][db.Type == 'core']
        gene_table = db.getTable('gene')
        
        # after release 65, the gene_id_table is removed. The following is to maintain
        # support for earlier releases
        release_ge_65 = self.GeneralRelease >= 65
        if release_ge_65:
            gene_id_table = None
        else:
            gene_id_table = db.getTable('gene_stable_id')
        
        assert Symbol or Description or StableId or BioType, "no valid argument provided"
        if Symbol:
            condition = xref_table.c.display_label==Symbol
        elif StableId and release_ge_65:
            condition = gene_table.c.stable_id==StableId
        elif StableId:
            condition = gene_id_table.c.stable_id==StableId
        else:
            condition = self._get_biotype_description_condition(gene_table, Description, BioType, like)
        
        query = self._build_gene_query(db, condition, gene_table, gene_id_table, xref_table)
        
        return query
    
    def makeLocation(self, CoordName, Start=None, End=None, Strand=1,
                ensembl_coord=False):
        """returns a location in the genome"""
        return Coordinate(self, CoordName=CoordName, Start=Start, End=End,
                            Strand=Strand, ensembl_coord=ensembl_coord)
    
    def getGeneByStableId(self, StableId):
        """returns the gene matching StableId, or None if no record found"""
        query = self._get_gene_query(self.CoreDb, StableId=StableId)
        try:
            record = list(query.execute())[0]
            gene = Gene(self, self.CoreDb, data=record)
        except IndexError:
            gene = None
        return gene
    
    def getGenesMatching(self, Symbol=None, Description=None, StableId=None,
                         BioType=None, like=True):
        """Symbol: HGC gene symbol, case doesn't matter
        description: a functional description
        StableId: the ensebl identifier
        BioType: the biological encoding type"""
        # TODO additional arguments to satisfy: external_ref, go_terms
        if Symbol is not None:
            Symbol = Symbol.lower()
        # biotype -> gene
        # description -> gene
        # Symbols -> xref
        # StableId -> gene_stable_id
        # XREF table calls
        # for gene symbols, these need to be matched against the display_label
        # attribute of core.xref table
        # for description, these need to be matched against the description 
        # field of the xref table
        
        # TODO catch conditions where user passes in both a symbol and a
        # biotype
        args = dict(Symbol=Symbol, Description=Description, 
                     StableId=StableId, BioType=BioType, like=like)
        query = self._get_gene_query(self.CoreDb, **args)
        records = query.execute()
        if records.rowcount == 0 and Symbol is not None:
            # see if the symbol has a synonym
            Symbol = self._get_symbol_from_synonym(self.CoreDb, Symbol)
            if Symbol is not None:
                args['Symbol'] = Symbol
                records = self._get_gene_query(self.CoreDb, **args).execute()
            else:
                records = []
        
        for record in records:
            gene = Gene(self, self.CoreDb, data=record)
            yield gene
    
    def getEstMatching(self, StableId):
        """returns an Est object from the otherfeatures db with the StableId"""
        query = self._get_gene_query(self.OtherFeaturesDb, StableId=StableId)
        records = query.execute()
        for record in records:
            yield Est(self,self.OtherFeaturesDb,StableId=StableId,data=record)
    
    def _get_seq_region_id(self, CoordName):
        """returns the seq_region_id for the provided CoordName"""
        seq_region_table = self.CoreDb.getTable('seq_region')
        coord_systems = CoordSystem(core_db=self.CoreDb)
        coord_system_ids = [k for k in coord_systems if not isinstance(k, str)]
        record = sql.select([seq_region_table.c.seq_region_id],
                    sql.and_(seq_region_table.c.name == CoordName,
                seq_region_table.c.coord_system_id.in_(coord_system_ids)))
        record = asserted_one(record.execute().fetchall())
        return record['seq_region_id']
    
    def _get_simple_features(self, db, klass, target_coord, query_coord,
                             where_feature):
        """returns feature_type records for the query_coord from the
        simple_feature table. The returned coord is referenced to
        target_coord. At present, only CpG islands being queried."""
        simple_feature_table = db.getTable('simple_feature')
        feature_types = ['CpGisland']
        feature_type_ids=[self._feature_type_ids.get(f) for f in feature_types]
        # fix the following
        query = sql.select([simple_feature_table],
            sql.and_(simple_feature_table.c.analysis_id.in_(feature_type_ids),
            simple_feature_table.c.seq_region_id == query_coord.seq_region_id))
        query = location_query(simple_feature_table,query_coord.EnsemblStart,
                        query_coord.EnsemblEnd, query=query,
                        where=where_feature)
        records = query.execute()
        for record in records:
            coord = Coordinate(self, CoordName=query_coord.CoordName,
                            Start=record['seq_region_start'],
                            End = record['seq_region_end'],
                            seq_region_id=record['seq_region_id'],
                            Strand = record['seq_region_strand'],
                            ensembl_coord=True)
            if query_coord.CoordName != target_coord.CoordName:
                coord = asserted_one(get_coord_conversion(coord, target_coord.CoordType, self.CoreDb))[1]
                
            # coord = coord.makeRelativeTo(query_coord) #TODO: fix here if query_coord and target_coord have different coordName
            # coord = coord.makeRelativeTo(target_coord, False)
            yield klass(self, db, Location=coord, Score=record['score'])
    
    def _get_repeat_features(self, db, klass, target_coord, query_coord,
                             where_feature):
        """returns Repeat region instances"""
        # we build repeats using coordinates from repeat_feature table
        # the repeat_consensus_id is required to get the repeat name, class
        # and type
        repeat_feature_table = db.getTable('repeat_feature')
        query = sql.select([repeat_feature_table],
            repeat_feature_table.c.seq_region_id == query_coord.seq_region_id)
        query = location_query(repeat_feature_table, query_coord.EnsemblStart,
                    query_coord.EnsemblEnd, query=query, where=where_feature)
        for record in query.execute():
            coord = Coordinate(self, CoordName=query_coord.CoordName,
                            Start=record['seq_region_start'],
                            End = record['seq_region_end'],
                            seq_region_id=record['seq_region_id'],
                            Strand = record['seq_region_strand'],
                            ensembl_coord=True)
            if query_coord.CoordName != target_coord.CoordName:
                coord = asserted_one(get_coord_conversion(coord, target_coord.CoordType, self.CoreDb))[1]
            # coord = coord.makeRelativeTo(query_coord) #TODO: fix here if query_coord and target_coord have different coordName
            # coord = coord.makeRelativeTo(target_coord, False)
            yield klass(self, db, Location=coord, Score=record['score'],
                        data=record)
    
    def _get_gene_features(self, db, klass, target_coord, query_coord,
                           where_feature):
        """returns all genes"""
        xref_table = [None, db.getTable('xref')][db.Type == 'core']
        gene_table = db.getTable('gene')
        
        # after release 65, the gene_id_table is removed. The following is to maintain
        # support for earlier releases.
        if self.GeneralRelease >= 65:
            gene_id_table = None
        else:
            gene_id_table = db.getTable('gene_stable_id')
        
        # note gene records are at chromosome, not contig, level
        condition = gene_table.c.seq_region_id == query_coord.seq_region_id
        query = self._build_gene_query(db, condition, gene_table, gene_id_table, xref_table)
        query = location_query(gene_table, query_coord.EnsemblStart,
                    query_coord.EnsemblEnd, query=query, where=where_feature)
        
        for record in query.execute():
            new = Coordinate(self, CoordName=query_coord.CoordName,
                            Start=record['seq_region_start'],
                            End = record['seq_region_end'],
                            Strand = record['seq_region_strand'], 
                            seq_region_id=record['seq_region_id'],
                            ensembl_coord=True)
            if query_coord.CoordName != target_coord.CoordName:
                coord = asserted_one(get_coord_conversion(coord, target_coord.CoordType, self.CoreDb))[1]
            
            # TODO: check coord, used 'new' here. where is coord (above line) used? 
            gene = klass(self, db, Location=new, data=record)
            yield gene
        
    
    def _get_variation_features(self, db, klass, target_coord, query_coord,
                        where_feature):
        """returns variation instances within the specified region"""
        # variation features at supercontig level
        var_feature_table = self.VarDb.getTable('variation_feature')
        # note gene records are at chromosome, not contig, level
        query = sql.select([var_feature_table],
            var_feature_table.c.seq_region_id == query_coord.seq_region_id)
        query = location_query(var_feature_table, query_coord.EnsemblStart,
                    query_coord.EnsemblEnd, query=query, where=where_feature)
        for record in query.execute():
            yield klass(self, self.CoreDb, Symbol=record['variation_name'],
                            data=record)
    
    def _get_feature_coord_levels(self, feature_types):
        dbs = dict(core_db = self.CoreDb)
        if 'variation' in feature_types:
            dbs["var_db"] = self.VarDb
        if 'est' in feature_types:
            dbs["otherfeature_db"] = self.OtherFeaturesDb
        feature_coord_levels = self._feature_coord_levels(self.Species, 
                                feature_types = feature_types,**dbs)
        return feature_coord_levels
    
    def _feature_coord_levels(self):
        if str(self._feature_coord_levels):
            return self._feature_coord_levels
        feature_types = ['gene', 'est', 'variation', 'cpg', 'repeat']
        feature_coord_levels = self._get_feature_coord_levels(feature_types)
        return self._feature_coord_levels
    
    FeatureCoordLevels = property(_feature_coord_levels)
    
    def getFeatures(self, region=None, feature_types=None, where_feature=None,
                    CoordName=None, Start=None, End=None, Strand=None,
                    ensembl_coord=False):
        """returns Region instances for the specified location"""
        if isinstance(feature_types, str):
            feature_types = [feature_types]
        feature_types = [ft.lower() for ft in feature_types]
        feature_coord_levels = self._get_feature_coord_levels(feature_types)
        
        if region is None:
            seq_region_id = self._get_seq_region_id(CoordName)
            region = Coordinate(self,CoordName=CoordName, Start=Start,
                        End=End,
                        Strand = convert_strand(Strand),
                        seq_region_id=seq_region_id,
                        ensembl_coord=ensembl_coord)
        elif hasattr(region, 'Location'):
            region = region.Location
        
        coord = region
        # the coordinate system at which locations are to be referenced, and
        # the processing function
        target_coords_funcs = \
            dict(cpg = (self._get_simple_features, CpGisland),
                 repeat = (self._get_repeat_features, Repeat),
                 gene = (self._get_gene_features, Gene),
                 est = (self._get_gene_features, Est),
                 variation = (self._get_variation_features, Variation))
        
        known_types = set(target_coords_funcs.keys())
        if not set(feature_types) <= known_types:
            raise RuntimeError, 'Unknown feature[%s], valid feature_types \
                are: %s' % (set(feature_types)^known_types, known_types)
        
        for feature_type in feature_types:
            target_func, target_class = target_coords_funcs[feature_type]
            db = self.CoreDb
            if feature_type == 'est':
                db = self.OtherFeaturesDb
            
            feature_coords = feature_coord_levels[feature_type].levels
            for feature_coord in feature_coords:
                chrom_other_coords = get_coord_conversion(coord, feature_coord,
                                            db, where=where_feature)
                for chrom_coord, other_coord in chrom_other_coords:
                    for region in target_func(db, target_class, chrom_coord,
                                            other_coord, where_feature):
                        yield region
    
    def getVariation(self, Effect=None, Symbol=None, like=True,
                     validated=False):
        """returns a generator of Variation instances
        
        Arguments:
            - Effect: the coding impact, eg. nonsynonymous
            - like: Effect is exactly matched against records like that
              provided
            - Symbol: the external or ensembl identifier - returns the exact
              match
            - validated: variant has validation_status != None"""
        var_feature_table = self.VarDb.getTable('variation_feature')
        
        assert Effect or Symbol, "No arguments provided"
        #  if we don't have Symbol, then we deal with Effect
        if Effect is not None:
            if like:
                query = \
                    var_feature_table.c.consequence_type.like('%'+Effect+'%')
            else:
                query = var_feature_table.c.consequence_type == Effect
        else:
            query = var_feature_table.c.variation_name == Symbol
        
        if validated:
            # in Release 65, the default validated status is now ''
            # why?? thanks Ensembl!
            null = None
            if int(self.Release) >= 65:
                null = ''
                
            query = sql.and_(query,var_feature_table.c.validation_status!=null)
        
        query = sql.select([var_feature_table],
                    query).order_by(var_feature_table.c.seq_region_start)
        for record in query.execute():
            yield Variation(self, self.CoreDb, Effect = Effect, Symbol=Symbol,
                            data=record)
        
    
    def getRegion(self, region=None, CoordName=None, Start=None, End=None,
                  Strand=None, ensembl_coord=False):
        """returns a single generic region for the specified coordinates
        Arguments:
            - region: a genomic region or a Coordinate instance
            - ensembl_coords: if True, follows indexing system of Ensembl
              where indexing starts at 1"""
        if region is None:
            seq_region_id = self._get_seq_region_id(CoordName)
            region = Coordinate(self,CoordName=CoordName, Start=Start,
                        End=End,
                        Strand = convert_strand(Strand),
                        seq_region_id=seq_region_id,
                        ensembl_coord=ensembl_coord)
        elif hasattr(region, 'Location'):
            region = region.Location
        
        return GenericRegion(self, self.CoreDb, CoordName=CoordName,
                             Start=Start, End=End, Strand=Strand,
                            Location=region, ensembl_coord=ensembl_coord)
    
    def getDistinct(self, property_type):
        """returns the Ensembl data-bases distinct values for the named
        property_type.
        
        Arguments:
            - property_type: valid values are biotype, status, effect"""
        property_type = property_type.lower()
        if property_type == 'effect':
            db = self.VarDb
        else:
            db = self.CoreDb
        
        property_map = {'effect': ('variation_feature', 'consequence_type'),
                        'biotype': ('gene', 'biotype'),
                        'status': ('gene', 'status')}
        
        if property_type not in property_map:
            raise RuntimeError,\
                "ERROR: Unknown property type: %s" % property_type
        
        table_name, column = property_map[property_type]
        return list(db.getDistinct(table_name, column))

