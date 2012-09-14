import sqlalchemy as sql
from numpy import empty

from cogent.util.table import Table
from cogent.db.ensembl.species import Species as _Species
from cogent.db.ensembl.util import NoItemError, asserted_one
from cogent.db.ensembl.host import get_ensembl_account, get_latest_release
from cogent.db.ensembl.database import Database
from cogent.db.ensembl.assembly import Coordinate, location_query
from cogent.db.ensembl.genome import Genome
from cogent.db.ensembl.related_region import RelatedGenes, SyntenicRegions

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Hua Ying", "Jason Merkin"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

class Compara(object):
    """comaparison among genomes"""
    def __init__(self, species, Release, account=None, pool_recycle=None,
            division=None):
        assert Release, 'invalid release specified'
        self.Release = str(Release)
        if account is None:
            account = get_ensembl_account(release=Release)
        self._account = account
        self._pool_recycle = pool_recycle
        self._compara_db = None
        sp = sorted([_Species.getSpeciesName(sp) for sp in set(species)])
        self.Species = tuple(sp)
        self._genomes = {}
        self._attach_genomes()
        
        self._species_id_map = None
        self._species_db_map = None
        self._species_set = None
        self._method_species_link = None
        self.division = division
    
    def _attach_genomes(self):
        for species in self.Species:
            attr_name = _Species.getComparaName(species)
            genome = Genome(Species=species, Release=self.Release,
                            account=self._account)
            self._genomes[species] = genome
            setattr(self, attr_name, genome)
    
    def __str__(self):
        my_type = self.__class__.__name__
        return "%s(Species=%s; Release=%s; connected=%s)" % \
                    (my_type, self.Species, self.Release,
                    self.ComparaDb is not None)
    
    def _connect_db(self):
        # TODO can the connection be all done in init?
        connection = dict(account=self._account, release=self.Release,
                        pool_recycle=self._pool_recycle)
        if self._compara_db is None:
            self._compara_db = Database(db_type='compara',
                                        division=self.division, **connection)
    
    def _get_compara_db(self):
        self._connect_db()
        return self._compara_db
    
    ComparaDb = property(_get_compara_db)
    
    def _make_species_id_map(self):
        """caches the taxon id's for the self.Species"""
        if self._species_id_map is not None:
            return self._species_id_map
        ncbi_table = self.ComparaDb.getTable('ncbi_taxa_name')
        conditon = sql.select([ncbi_table.c.taxon_id, ncbi_table.c.name],
                    ncbi_table.c.name.in_([sp for sp in self.Species]))
        # TODO this should make the dict values the actual Genome instances
        id_genome = []
        for r in conditon.execute():
            id_genome += [(r['taxon_id'], self._genomes[r['name']])]
        self._species_id_map = dict(id_genome)
        return self._species_id_map
    
    taxon_id_species = property(_make_species_id_map)
    
    def _get_genome_db_ids(self):
        if self._species_db_map is not None:
            return self._species_db_map
        genome_db_table = self.ComparaDb.getTable('genome_db')
        query = sql.select([genome_db_table.c.genome_db_id,
                           genome_db_table.c.taxon_id],
                 genome_db_table.c.taxon_id.in_(self.taxon_id_species.keys()))
        records = query.execute()
        self._species_db_map = \
                    dict([(r['genome_db_id'],r['taxon_id']) for r in records])
        return self._species_db_map
    
    genome_taxon = property(_get_genome_db_ids)
    
    def _get_species_set(self):
        if self._species_set is not None:
            return self._species_set
        # we make sure the species set contains all species
        species_set_table = self.ComparaDb.getTable('species_set')
        query = sql.select([species_set_table],
               species_set_table.c.genome_db_id.in_(self.genome_taxon.keys()))
        species_sets = {}
        for record in query.execute():
            gen_id = record['genome_db_id']
            sp_set_id = record['species_set_id']
            if sp_set_id in species_sets:
                species_sets[sp_set_id].update([gen_id])
            else:
                species_sets[sp_set_id] = set([gen_id])
        
        expected = set(self.genome_taxon.keys())
        species_set_ids = []
        for sp_set, gen_id in species_sets.items():
            if expected <= gen_id:
                species_set_ids.append(sp_set)
        self._species_set = species_set_ids
        return self._species_set
    
    species_set = property(_get_species_set)
    
    def _get_method_link_species_set(self):
        if self._method_species_link is not None:
            return self._method_species_link
        
        method_link_table = self.ComparaDb.getTable('method_link')
        query = sql.select([method_link_table],
            method_link_table.c['class'].like('%'+'alignment'+'%'))
        methods = query.execute().fetchall()
        method_link_ids = dict([(r['method_link_id'], r) for r in methods])
        method_link_species_table = \
                            self.ComparaDb.getTable('method_link_species_set')
        query = sql.select([method_link_species_table],
            sql.and_(
            method_link_species_table.c.species_set_id.in_(self.species_set),
            method_link_species_table.c.method_link_id.in_(
                                                    method_link_ids.keys())))
        records = query.execute().fetchall()
        # store method_link_id, type, species_set_id,
        # method_link_species_set.name, class
        header = ['method_link_species_set_id', 'method_link_id',
                'species_set_id', 'align_method', 'align_clade']
        rows = []
        for record in records:
            ml_id = record['method_link_id']
            sp_set_id = record['species_set_id']
            ml_sp_set_id = record['method_link_species_set_id']
            clade_name = record['name']
            aln_name = method_link_ids[ml_id]['type']
            rows += [[ml_sp_set_id, ml_id, sp_set_id, aln_name, clade_name]]
        
        if rows == []:
            rows = empty((0,len(header)))
        
        t = Table(header=header, rows=rows, space=2, row_ids=True,
                  title='Align Methods/Clades')
        self._method_species_link = t
        return t
    
    method_species_links = property(_get_method_link_species_set)
    
    def getRelatedGenes(self, gene_region=None, StableId=None,
                                Relationship=None, DEBUG=False):
        """returns a RelatedGenes instance.
        
        Arguments:
            - gene_region: a Gene instance
            - StableId: ensembl stable_id identifier
            - Relationship: the types of related genes sought"""
        assert gene_region is not None or StableId is not None,\
            "No identifier provided"
        assert Relationship is not None, "No Relationship specified"
        
        # TODO understand why this has become necessary to suppress warnings
        # in SQLAlchemy 0.6
        Relationship = u'%s' % Relationship
        
        StableId = StableId or gene_region.StableId
        
        member_table = self.ComparaDb.getTable('member')
        homology_member_table = self.ComparaDb.getTable('homology_member')
        homology_table = self.ComparaDb.getTable('homology')
        
        member_ids = sql.select([member_table.c.member_id],
            member_table.c.stable_id == StableId)
        
        member_ids = [r['member_id'] for r in member_ids.execute()]
        if not member_ids:
            return None
        
        if DEBUG: print "member_ids", member_ids
        
        homology_ids = sql.select([homology_member_table.c.homology_id,
                          homology_member_table.c.member_id],
                          homology_member_table.c.member_id.in_(member_ids))
        homology_ids = [r['homology_id'] for r in homology_ids.execute()]
        if not homology_ids:
            return None
        
        if DEBUG: print "1 - homology_ids", homology_ids
        
        homology_records = \
                sql.select([homology_table.c.homology_id,
                           homology_table.c.description,
                           homology_table.c.method_link_species_set_id],
                    sql.and_(homology_table.c.homology_id.in_(homology_ids),
                    homology_table.c.description == Relationship))
        
        homology_ids = []
        for r in homology_records.execute():
            homology_ids.append((r["homology_id"],
                        (r["description"], r["method_link_species_set_id"])))
        homology_ids = dict(homology_ids)
        
        if DEBUG: print "2 - homology_ids", homology_ids
        if not homology_ids:
            return None
        
        ortholog_ids = sql.select([homology_member_table.c.member_id,
                                homology_member_table.c.homology_id],
                homology_member_table.c.homology_id.in_(homology_ids.keys()))
        
        ortholog_ids = dict([(r['member_id'], r['homology_id']) \
                                      for r in ortholog_ids.execute()])
        
        if DEBUG: print "ortholog_ids", ortholog_ids
        if not ortholog_ids:
            return None
        
        # could we have more than one here?
        relationships = set()
        for memid, homid in ortholog_ids.items():
            relationships.update([homology_ids[homid][0]])
        relationships = tuple(relationships)
        
        gene_set = sql.select([member_table],
                sql.and_(member_table.c.member_id.in_(ortholog_ids.keys()),
                  member_table.c.taxon_id.in_(self.taxon_id_species.keys())))
        data = []
        for record in gene_set.execute():
            genome = self.taxon_id_species[record['taxon_id']]
            StableId = record['stable_id']
            gene = list(genome.getGenesMatching(StableId=StableId))
            assert len(gene) == 1, "Error in selecting genes: %s" % gene
            gene = gene[0]
            gene.Location.Strand = record['chr_strand']
            data += [gene]
        
        if not data:
            return None
        return RelatedGenes(self, data, Relationships=relationships)
    
    def _get_dnafrag_id_for_coord(self, coord):
        """returns the dnafrag_id for the coordnate"""
        dnafrag_table = self.ComparaDb.getTable('dnafrag')
        genome_db_table = self.ComparaDb.getTable('genome_db')
        
        # column renamed between versions
        prefix = coord.genome.Species.lower()
        if int(self.Release) > 58:
            prefix = _Species.getEnsemblDbPrefix(prefix)
        
        query = sql.select([dnafrag_table.c.dnafrag_id,
                           dnafrag_table.c.coord_system_name],
                  sql.and_(dnafrag_table.c.genome_db_id ==\
                                            genome_db_table.c.genome_db_id,
                                genome_db_table.c.name == prefix,
                                dnafrag_table.c.name == coord.CoordName))
        try:
            record = asserted_one(query.execute().fetchall())
            dnafrag_id = record['dnafrag_id']
        except NoItemError:
            raise RuntimeError, 'No DNA fragment identified'
        return dnafrag_id
    
    def _get_genomic_align_blocks_for_dna_frag_id(self, method_clade_id,
                                                  dnafrag_id, coord):
        genomic_align_table = self.ComparaDb.getTable('genomic_align')
        query = sql.select([genomic_align_table.c.genomic_align_id,
                           genomic_align_table.c.genomic_align_block_id],
                sql.and_(genomic_align_table.c.method_link_species_set_id ==\
                                                    method_clade_id,
                        genomic_align_table.c.dnafrag_id == dnafrag_id))
        query = location_query(genomic_align_table,
                               coord.EnsemblStart,
                               coord.EnsemblEnd,
                               start_col = 'dnafrag_start',
                               end_col = 'dnafrag_end',
                               query = query)
        
        return query.execute().fetchall()
    
    def _get_joint_genomic_align_dnafrag(self, genomic_align_block_id):
        genomic_align_table = self.ComparaDb.getTable('genomic_align')
        dnafrag_table = self.ComparaDb.getTable('dnafrag')
        query = sql.select([genomic_align_table.c.genomic_align_id,
                           genomic_align_table.c.genomic_align_block_id,
                           genomic_align_table.c.dnafrag_start,
                           genomic_align_table.c.dnafrag_end,
                           genomic_align_table.c.dnafrag_strand,
                           dnafrag_table],
            sql.and_(genomic_align_table.c.genomic_align_block_id == \
                                                        genomic_align_block_id,
                genomic_align_table.c.dnafrag_id == dnafrag_table.c.dnafrag_id,
                dnafrag_table.c.genome_db_id.in_(self.genome_taxon.keys())))
        return query.execute().fetchall()
    
    def getSyntenicRegions(self, Species=None, CoordName=None, Start=None,
            End=None, Strand=1, ensembl_coord=False, region=None,
            align_method=None, align_clade=None, method_clade_id=None):
        """returns a SyntenicRegions instance
        
        Arguments:
            - Species: the species name
            - CoordName, Start, End, Strand: the coordinates for the region
            - ensembl_coord: whether the coordinates are in Ensembl form
            - region: a region instance or a location, in which case the
              CoordName etc .. arguments are ignored
            - align_method, align_clade: the alignment method and clade to use
              Note: the options for this instance can be found by printing
              the method_species_links attribute of this object.
            - method_clade_id: over-rides align_method/align_clade. The entry
              in method_species_links under method_link_species_set_id
              """
        assert (align_method and align_clade) or method_clade_id, \
                'Must specify (align_method & align_clade) or method_clade_id'
        if method_clade_id is None:
            for row in self.method_species_links:
                if align_method.lower() in row['align_method'].lower() and\
                        align_clade.lower() in row['align_clade'].lower():
                    method_clade_id = row['method_link_species_set_id']
        
        if method_clade_id is None:
            raise RuntimeError, "Invalid align_method[%s] or align_clade "\
                                "specified[%s]" % (align_method, align_clade)
        
        if region is None:
            ref_genome = self._genomes[_Species.getSpeciesName(Species)]
            region = ref_genome.makeLocation(CoordName=CoordName,
                                Start=Start, End=End, Strand=Strand,
                                ensembl_coord=ensembl_coord)
        elif hasattr(region, 'Location'):
            region = region.Location
        
        # make sure the genome instances match
        ref_genome = self._genomes[region.genome.Species]
        if ref_genome is not region.genome:
            # recreate region from our instance
            region = ref_genome.makeLocation(CoordName=region.CoordName,
                                Start=region.Start, End=region.End,
                                Strand=region.Strand)
        
        ref_dnafrag_id = self._get_dnafrag_id_for_coord(region)
        blocks=self._get_genomic_align_blocks_for_dna_frag_id(method_clade_id,
                                                    ref_dnafrag_id, region)
        for block in blocks:
            genomic_align_block_id = block['genomic_align_block_id']
            # we get joint records for these identifiers from
            records = self._get_joint_genomic_align_dnafrag(
                                                genomic_align_block_id)
            members = []
            ref_location = None
            for record in records:
                taxon_id = self.genome_taxon[record.genome_db_id]
                genome = self.taxon_id_species[taxon_id]
                # we have a case where we getback different coordinate system
                # results for the ref genome. We keep only those that match
                # the CoordName of region
                
                if genome is region.genome and \
                        record.name == region.CoordName:
                    # this is the ref species and we adjust the ref_location
                    # for this block
                    diff_start = record.dnafrag_start-region.EnsemblStart
                    shift_start = [0, diff_start][diff_start > 0]
                    diff_end = record.dnafrag_end-region.EnsemblEnd
                    shift_end = [diff_end, 0][diff_end > 0]
                    try:
                        ref_location = region.resized(shift_start, shift_end)
                    except ValueError:
                        # we've hit some ref genome fragment that matches
                        # but whose coordinates aren't right
                        continue
                elif genome is region.genome:
                    continue
                members += [(genome, record)]
            assert ref_location is not None, "Failed to make the reference"\
                                                            " location"
            yield SyntenicRegions(self, members, ref_location=ref_location)
        
    
    def getDistinct(self, property_type):
        """returns the Ensembl data-bases distinct values for the named
        property_type.
        
        Arguments:
            - property_type: valid values are relationship"""
        property_type = property_type.lower()
        db = self.ComparaDb
        property_map = {'relationship': ('homology', 'description'),
                        'clade': ('method_link_species_set', 'name')}
        if property_type not in property_map:
            raise RuntimeError, "ERROR: Unknown property type: %s"%property_type
        table_name, column = property_map[property_type]
        return list(db.getDistinct(table_name, column))

