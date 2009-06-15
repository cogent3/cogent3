#!/usr/bin/env python
import sqlalchemy as sql

from cogent.util.table import Table
from cogent.db.ensembl.assembly import CoordSystem
from cogent.db.ensembl.species import Species as _Species

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

class _FeatureLevelRecord(object):
    def __init__(self, feature_type, coord_system_names, default=None):
        self.feature_type = feature_type
        self.levels = ', '.join(coord_system_names)
        if default is not None:
            self.level_in_use = default   # should use rank information from coord_system
        else:
            self.level_in_use = None
    
    def __str__(self):
        return 'feature = %s; Levels = %s; Level in use = %s'%(self.feature_type, self.levels, self.level_in_use)
    

class FeatureCoordLevelsCache(object):
    _species_feature_levels = {}
    _species_feature_dbs = {}
    def __init__(self, species):
        self.Species = _Species.getSpeciesName(species)
    
    def __repr__(self):
        """print table format"""
        header = ['Type', 'Levels', 'Level_in_use']
        result = []
        for species in self._species_feature_levels.keys():
            feature_levels = self._species_feature_levels[species]
            collate = []
            for feature in feature_levels.keys():
                collate.append([feature, feature_levels[feature].levels, feature_levels[feature].level_in_use])
            t = Table(header, collate, title=species)
            result.append(str(t))
        result = '\n'.join(result)
        return result
    
    def _get_meta_coord_records(self, db):
        meta_coord = db.getTable('meta_coord')
        if 'core' in str(db.db_name):
            query = sql.select([meta_coord]).where(meta_coord.c.table_name.\
                        in_(['gene', 'simple_feature', 'repeat_feature']))
            query = query.order_by(meta_coord.c.table_name)
        elif 'variation' in str(db.db_name): 
            query = sql.select([meta_coord]).where(meta_coord.c.table_name == 'variation_feature')
        else:
            assert 'otherfeature' in str(db.db_name)
            query = sql.select([meta_coord]).where(meta_coord.c.table_name == 'gene')
        records = query.execute().fetchall()
        return records
    
    def _add_species_feature_levels(self, species, records, db_type, coord_system):
        if db_type == 'core':
            features = ['cpg', 'repeat', 'gene', 'est']
            tables = ['simple_feature', 'repeat_feature', 'gene', 'gene']
        elif db_type == 'var':
            features, tables = ['variation'], ['variation_feature']
        else:
            assert db_type == 'otherfeature'
            features, tables = ['est'], ['gene']
            
        for feature, table_name in zip(features, tables):
            feature_coord_ids = [r['coord_system_id'] for r in records if r['table_name'] == table_name]
            feature_coord_systems = [coord_system[coord_id] for coord_id in feature_coord_ids]
            feature_coord_systems = [(s.rank, s.name) for s in feature_coord_systems]
            feature_coord_systems.sort()
            levels = [name for rank, name in feature_coord_systems]
            self._species_feature_levels[species][feature] = _FeatureLevelRecord(feature, levels, levels[0])
    
    def _set_species_feature_levels(self, species, core_db, feature_types, var_db, otherfeature_db):
        if species not in self._species_feature_levels:
            self._species_feature_levels[species] = {}
            self._species_feature_dbs[species] = []
        coord_system = CoordSystem(core_db = core_db)
        if set(feature_types).intersection(set(['cpg', 'repeat', 'gene'])):
            if 'core_db' not in self._species_feature_dbs[species]:
                self._species_feature_dbs[species].append('core_db')
                records = self._get_meta_coord_records(core_db)
                self._add_species_feature_levels(species, records, 'core', coord_system)
        if 'variation' in feature_types:
            if 'var_db' not in self._species_feature_dbs[species]:
                self._species_feature_dbs[species].append('var_db')
                assert var_db is not None
                records = self._get_meta_coord_records(var_db)
                self._add_species_feature_levels(species, records, 'var', coord_system)
        if 'est' in feature_types:
            if 'otherfeature_db' not in self._species_feature_dbs[species]:
                self._species_feature_dbs[species].append('otherfeature_db')
                assert otherfeature_db is not None
                records = self._get_meta_coord_records(otherfeature_db)
                self._add_species_feature_levels(species, records, 'otherfeature', coord_system)
    
    def __call__(self, species = None, core_db=None, feature_types=None, var_db=None, otherfeature_db=None):
        if 'variation' in feature_types:
            assert var_db is not None
        species = _Species.getSpeciesName(core_db.db_name.Species or species)
        self._set_species_feature_levels(species, core_db, feature_types, var_db, otherfeature_db)
        return self._species_feature_levels[species]
    


class FeatureCoordLevels(FeatureCoordLevelsCache):
    def __init__(self, species):
        self.Species = _Species.getSpeciesName(species)
    
    def __repr__(self):
        """print table format"""
        header = ['Type', 'Levels', 'Levels_in_use']
        if self.Species not in self._species_feature_levels:
            result = ''
        else:
            collate = []
            feature_levels = self._species_feature_levels[self.Species]
            for feature in feature_levels.keys():
                record = feature_levels[feature]
                collate.append([feature, record.levels, record.level_in_use])
            result = str(Table(header, collate, title=self.Species))
        return result
    

