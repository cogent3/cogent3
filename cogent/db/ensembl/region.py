import sys
import sqlalchemy as sql

from cogent import DNA
from cogent.core.annotation import Feature
from cogent.core.location import Map
from cogent.util.table import Table
from cogent.db.ensembl.util import LazyRecord, asserted_one, DisplayString, \
                                    NoItemError
from cogent.db.ensembl.assembly import Coordinate, CoordSystem, \
                                location_query, assembly_exception_coordinate
from cogent.db.ensembl.sequence import get_sequence

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

DEFAULT_PARENT_LENGTH = 2 **30

# some common string display formatting
_quoted = lambda x : DisplayString(x, with_quotes=True)
_limit_words = lambda x : DisplayString(x, with_quotes=True, num_words=3)

class _Region(LazyRecord):
    """a simple genomic region object"""
    Type = None
    
    def __init__(self):
        super(_Region, self).__init__()
        self._attr_ensembl_table_map = None
        self._location_column_prefix = 'seq_region_'
    
    def __len__(self):
        return len(self.Location)
    
    def __cmp__(self, other):
        return cmp(self.Location, other.Location)
    
    def _make_location(self):
        row = self._table_rows[self._attr_ensembl_table_map['Location']]
        if row is None:
            return
        seq_region_id = row['%sid' % self._location_column_prefix]
        start = row['%sstart' % self._location_column_prefix]
        end = row['%send' % self._location_column_prefix]
        strand = row['%sstrand' % self._location_column_prefix]
        seq_region_table = self.db.getTable('seq_region')
        query = sql.select([seq_region_table.c.name],
                         seq_region_table.c.seq_region_id == seq_region_id)
        result = asserted_one(query.execute().fetchall())
        coord_name = result['name']
        
        coord = Coordinate(genome = self.genome, CoordName=coord_name,
                    Start=start, End=end, Strand=strand,
                    seq_region_id=seq_region_id,
                    ensembl_coord=True)
        self._cached['Location'] = coord
    
    def _get_location_record(self):
        """makes the Location data"""
        if not self._attr_ensembl_table_map['Location'] in self._table_rows:
            # we use a bit of magic to figure out what method will be required
            # this magic assumes the method for obtaining a record from a table
            # are named _get_tablename_record
            dep_record_func = getattr(self, '_get_%s_record' % \
                        self._attr_ensembl_table_map['Location'])
            dep_record_func()
        self._make_location()
    
    def _get_location(self):
        return self._get_cached_value('Location', self._get_location_record)
    
    Location = property(_get_location)
    
    def _get_sequence(self):
        if 'Seq' not in self._cached:
            try:
                seq = get_sequence(self.Location)
            except NoItemError:
                try:
                    alt_loc = assembly_exception_coordinate(self.Location)
                    seq = get_sequence(alt_loc)
                except NoItemError:
                    seq = DNA.makeSequence("N"*len(self))
            seq.Name = str(self.Location)
            self._cached['Seq'] = seq
        return self._cached['Seq']
    
    Seq = property(_get_sequence)
    
    def _get_symbol(self):
        # override in subclasses
        return None
    
    Symbol = property(_get_symbol)
    
    def getFeatures(self, feature_types, where_feature=None):
        """queries the parent genome for feature types corresponding to this
        region
        where_feature: the returned region can either lie 'within' this region,
        'overlap' this region, or 'span' this region"""
        return self.genome.getFeatures(self.Location,
                                        feature_types=feature_types,
                                        where_feature=where_feature)
    
    def _get_variants(self):
        """constructs the variants attribute"""
        if 'Variants' not in self._cached:
            variants = self.genome.getFeatures(feature_types='variation', region=self)
            self._cached['Variants'] = tuple(variants)
        
        return self._cached['Variants']
    
    Variants = property(_get_variants)
    
    def featureData(self, parent_map):
        symbol = self.Symbol or getattr(self, 'StableId', '')
        assert not parent_map.Reverse
        feat_map = parent_map[self.Location.Start:self.Location.End]
        if feat_map.useful:
            if self.Location.Strand == -1:
                # this map is relative to + strand
                feat_map = feat_map.reversed()
            data = (self.Type, str(symbol), feat_map)
        else:
            data = None
        return data
    
    def getAnnotatedSeq(self, feature_types=None, where_feature=None):
        regions = list(self.getFeatures(feature_types = feature_types, 
                        where_feature = where_feature))
        # seq_map is on the + strand, regardless the actual strand of sequence
        seq_map = Map(locations = [(self.Location.Start, self.Location.End)], 
                    parent_length = DEFAULT_PARENT_LENGTH)
        seq_map = seq_map.inverse()
        
        for region in regions:
            data = region.featureData(seq_map)
            if data is None:
                continue
            # this will consider the strand information of actual sequence
            feature_map = [data[-1],
                    data[-1].nucleicReversed()][self.Location.Strand == -1]
            self.Seq.addAnnotation(Feature, data[0], data[1], feature_map)
            
            if region.Type == 'gene':  # TODO: SHOULD be much simplified
                sub_data = region.subFeatureData(seq_map)
                for feature_type, feature_name, feature_map in sub_data:
                    if self.Location.Strand == -1:
                        # again, change feature map to -1 strand sequence if
                        # needed.
                        feature_map = feature_map.nucleicReversed()
                    self.Seq.addAnnotation(Feature, feature_type, 
                                            feature_name, feature_map)
        
        return self.Seq
    


class GenericRegion(_Region):
    """a generic genomic region"""
    
    Type = 'generic_region'
    
    def __init__(self, genome, db, Location=None, CoordName=None, Start=None,
                 End=None, Strand=1, ensembl_coord=False):
        super(GenericRegion, self).__init__()
        self.genome = genome
        self.db = db
        
        if Location is None and CoordName:
            self._get_seq_region_record(str(CoordName))
            if End is not None:
                assert self._table_rows['seq_region']['length'] > End, \
                        'Requested End[%s] too large' % End
            seq_region_id = self._table_rows['seq_region']['seq_region_id']
            Location = Coordinate(genome=genome, CoordName=str(CoordName),
                            Start=Start, End=End, Strand=Strand,
                            seq_region_id=seq_region_id,
                            ensembl_coord=ensembl_coord)
        
        if Location is not None:
            self._cached['Location'] = Location
    
    def __str__(self):
        my_type = self.__class__.__name__
        return "%s(Species='%s'; CoordName='%s'; Start=%s; End=%s;"\
                " length=%s; Strand='%s')" %  (my_type,
                self.genome.Species,
                self.Location.CoordName, self.Location.Start,
                self.Location.End, len(self), '-+'[self.Location.Strand>0])
    
    def _get_seq_region_record(self, CoordName):
        # override the _Region class method, since, we take the provided Start
        # etc .. attributes
        # CoordName comes from seq_region_table.c.name
        # matched, by coord_system_id, to default coord system
        seq_region_table = self.genome.db.getTable('seq_region')
        coord_systems = CoordSystem(core_db=self.genome.CoreDb)
        coord_system_ids = [k for k in coord_systems if not isinstance(k, str)]
        record = sql.select([seq_region_table],
                    sql.and_(seq_region_table.c.name == CoordName,
                seq_region_table.c.coord_system_id.in_(coord_system_ids)))
        record = asserted_one(record.execute().fetchall())
        self._table_rows['seq_region'] = record
    


class _StableRegion(GenericRegion):
    """region with a stable_id"""
    
    _member_types = None
    
    def __init__(self, genome, db, **kwargs):
        super(_StableRegion, self).__init__(genome, db, **kwargs)
    
    def __repr__(self):
        my_type = self.__class__.__name__
        return '%s(%s; %s)' % (my_type, self.genome.Species, self.StableId)
    
    def _get_record_for_stable_id(self):
        # subclasses need to provide a function for loading the correct
        # record for obtaining a stable_id
        table_name = self._attr_ensembl_table_map['StableId']
        
        if self.genome.GeneralRelease >= 65:
            func_name = '_get_%s_record' % (table_name + '_stable_id')
        else:
            func_name = '_get_%s_record' % table_name
        
        func = getattr(self, func_name)
        func()
        attr_column_map = [('StableId', 'stable_id', _quoted)]
        
        self._populate_cache_from_record(attr_column_map, table_name)
    
    def _get_stable_id(self):
        return self._get_cached_value('StableId',
                                       self._get_record_for_stable_id)
    
    StableId = property(_get_stable_id)
    
    def getMember(self, StableId, member_types=None):
        """returns the associated member with matching StableId or None if not
        found.
        
        Arguments:
            - member_types: the property to be searched, depends on self.Type.
              Transcripts for genes, Exons/TranslatedExons for Transcripts."""
        
        member_types = member_types or self._member_types
        if type(member_types) == str:
            member_types = [member_types]
        
        for member_type in member_types:
            member = getattr(self, member_type, None)
            if member is None:
                raise AttributeError,\
                    "%s doesn't have property %s" % (self.Type, member_type)
            for element in member:
                if element.StableId == StableId:
                    return element
        return None

class Gene(_StableRegion):
    """a gene region"""
    Type = 'gene'
    _member_types = ['Transcripts']
    
    def __init__(self, genome, db, StableId=None, Symbol=None, Location=None, data=None):
        """constructed by a genome instance"""
        super(Gene, self).__init__(genome, db, Location=Location)
        
        self._attr_ensembl_table_map = dict(StableId=['gene_stable_id',
                                                      'gene'][genome.GeneralRelease >= 65],
                                            Symbol='xref',
                                            Description='gene', BioType='gene', Location='gene',
                                            CanonicalTranscript='gene',
                                            Transcripts='transcript',
                                            Exons='transcript')
        
        if data is None:
            args = [dict(StableId=StableId), dict(Symbol=Symbol)][StableId is None]
            assert args
            data = asserted_one(list(self.genome._get_gene_query(db, **args).execute()))
        
        for name, func in \
                    [('StableId',self._get_gene_stable_id_record),
                     ('BioType', self._get_gene_record),
                     ('Description', self._get_gene_record),
                     ('Symbol', self._get_xref_record),
                     ('Location', self._get_gene_record)]:
            if name == 'Symbol' and 'display_label' not in data.keys(): # For EST
                continue
            self._table_rows[self._attr_ensembl_table_map[name]] = data
            func() # this populates the attributes
    
    def __str__(self):
        my_type = self.__class__.__name__
        vals = ['%s=%r' % (key, val) for key, val in self._cached.items() \
                                if val is not None]
        vals.sort()
        vals.insert(0, "Species='%s'" % self.genome.Species)
        return '%s(%s)' % (my_type, '; '.join(vals))
    
    def __repr__(self):
        my_type = self.__class__.__name__
        vals = ['%s=%r' % (key, val) for key, val in self._cached.items() \
                                if val is not None]
        vals.sort()
        vals.insert(0, 'Species=%r' % self.genome.Species)
        return '%s(%s)' % (my_type, '; '.join(vals))
    
    def _get_gene_record(self):
        """adds the gene data to self._table_rows"""
        attr_column_map = [('BioType', 'biotype', _quoted),
                           ('Status', 'status', _quoted),
                           ('Description', 'description', _limit_words)]
        # we set all the attributes that derive from this
        self._populate_cache_from_record(attr_column_map, 'gene')
        return
    
    def _get_gene_stable_id_record(self):
        """adds the gene_stable_id data to self._table_rows"""
        attr_column_map = [('StableId', 'stable_id', _quoted)]
        self._populate_cache_from_record(attr_column_map,
                                        self._attr_ensembl_table_map['StableId'])
        return
    
    def _get_xref_record(self):
        attr_column_map = [('Symbol','display_label', _quoted)]
        self._populate_cache_from_record(attr_column_map, 'xref')
        return
    
    def _get_biotype(self):
        return self._get_cached_value('BioType', self._get_gene_record)
    
    BioType = property(_get_biotype)
    
    def _get_symbol(self):
        if 'xref' in self._table_rows:
            return self._get_cached_value('Symbol', self._get_xref_record)
        self._set_null_values(['Symbol'])
        return self._cached['Symbol']
    
    Symbol = property(_get_symbol)
    
    def _get_description(self):
        return self._get_cached_value('Description', self._get_gene_record)
    
    Description = property(_get_description)
    
    def _get_status(self):
        return self._get_cached_value('Status', self._get_gene_record)
    
    Status = property(_get_status)
    
    def _make_canonical_transcript(self):
        if not 'gene' in self._table_rows:
            self._get_gene_record()
        canonical_id = self._table_rows['gene']['canonical_transcript_id']
        transcript_table = self.db.getTable('transcript')
        query = sql.select([transcript_table],
                 transcript_table.c.transcript_id==canonical_id)
        records = query.execute().fetchall()
        assert len(records) == 1,\
                "wrong number of records from CanonicalTranscript"
        record = records[0]
        transcript = Transcript(self.genome, self.db, canonical_id,
                                        data=record)
        self._cached['CanonicalTranscript'] = transcript
    
    def _get_canonical_transcript(self):
        return self._get_cached_value('CanonicalTranscript',
                                    self._make_canonical_transcript)
    
    CanonicalTranscript = property(_get_canonical_transcript)
    
    def _make_transcripts(self):
        if not 'gene' in self._table_rows:
            self._get_gene_record()
        gene_id = self._table_rows['gene']['gene_id']
        transcript_table = self.db.getTable('transcript')
        query = sql.select([transcript_table],
                           transcript_table.c.gene_id==gene_id)
        records = query.execute().fetchall()
        if not records:
            self._set_null_values(['Transcripts'], 'transcript')
            return
        transcripts = []
        for record in records:
            transcript_id = record['transcript_id']
            transcripts.append(Transcript(self.genome, self.db, transcript_id,
                                            data=record))
        self._cached['Transcripts'] = tuple(transcripts)
    
    def _get_transcripts(self):
        return self._get_cached_value('Transcripts', self._make_transcripts)
    
    Transcripts = property(_get_transcripts)
    
    def subFeatureData(self, parent_map):
        """returns data for making a cogent Feature. These can be
        automatically applied to the Seq by the getAnnotatedSeq method.
        Returns None if self lies outside parent's span.
        """
        features = []
        for transcript in self.Transcripts:
            transcript_data = transcript.featureData(parent_map)
            if transcript_data:
                features.append(transcript_data)
                data = transcript.subFeatureData(parent_map)
                features.extend(data)
        return features
    
    def getCdsLengths(self):
        """returns the Cds lengths from transcripts with the same biotype.
        returns None if no transcripts."""
        if self.Transcripts is self.NULL_VALUE:
            return None
        l = [ts.getCdsLength() for ts in self.Transcripts
                                            if ts.BioType == self.BioType]
        return l
    
    def getLongestCdsTranscript(self):
        """returns the Transcript with the longest Cds and the same biotype"""
        result = sorted([(ts.getCdsLength(), ts) for ts in self.Transcripts
                                            if ts.BioType == self.BioType])
        
        if result: # last result is longest
            result = result[-1][1]
        
        return result
    

class Transcript(_StableRegion):
    Type = 'transcript'
    _member_types = ['Exons', 'TranslatedExons']
    
    def __init__(self, genome, db, transcript_id, data, Location=None):
        """created by Gene"""
        super(Transcript, self).__init__(genome, db, Location=Location)

        self._attr_ensembl_table_map = dict(StableId=['transcript_stable_id',
                                                'transcript'][genome.GeneralRelease >= 65],
                                           Location='transcript',
                                           Status = 'transcript',
                                           TranslatedExons='translation')

        self._am_prot_coding = None
        self.transcript_id = transcript_id
        self._table_rows['transcript'] = data
        self._set_transcript_record()
    
    def _set_transcript_record(self):
        attr_column_map = [('BioType', 'biotype', _quoted),
                           ('Status', 'status', _quoted)]
        self._populate_cache_from_record(attr_column_map, 'transcript')
        self._am_prot_coding=self._cached['BioType'].lower()=='protein_coding'
    
    def _get_status(self):
        return self._cached['Status']
    
    Status = property(_get_status)
    
    def _get_biotype(self):
        return self._cached['BioType']
    
    BioType = property(_get_biotype)
    
    def _get_transcript_stable_id_record(self):
        table_name = self._attr_ensembl_table_map['StableId']
        if table_name in self._table_rows:
            return
        transcript_id = self.transcript_id
        table = self.db.getTable(table_name)
        query = sql.select([table], table.c.transcript_id == transcript_id)
        record = asserted_one(query.execute())
        self._table_rows[table_name] = record
        
    def _get_exon_transcript_records(self):
        transcript_id = self.transcript_id
        exon_transcript_table = self.db.getTable('exon_transcript')
        query = sql.select([exon_transcript_table],
                       exon_transcript_table.c.transcript_id == transcript_id)
        records = query.execute()
        exons = []
        for record in records:
            exons.append(Exon(self.genome, self.db, record['exon_id'],
                                record['rank']))
        exons.sort()
        self._cached['Exons'] = tuple(exons)
    
    def _get_exons(self):
        return self._get_cached_value('Exons',
                        self._get_exon_transcript_records)
    
    Exons = property(_get_exons)
    
    def _get_intron_transcript_records(self):
        if len(self.Exons) < 2:
            self._set_null_values(["Introns"])
            return
        
        exon_positions = [(exon.Location.Start, exon.Location.End) 
                            for exon in self.Exons]
        exon_positions.sort()
        end = exon_positions[-1][-1]
        exon_map = Map(locations=exon_positions, parent_length=end)
        intron_map = exon_map.shadow()
        
        intron_positions = [(span.Start, span.End) 
                            for span in intron_map.spans if span.Start != 0]
        
        chrom = self.Location.CoordName
        strand = self.Location.Strand
        introns = []
        rank = 1
        if strand == -1:
            intron_positions.reverse()
        for s, e in intron_positions:
            coord = self.genome.makeLocation(CoordName=chrom, Start=s, End=e, 
                                            Strand=strand, ensembl_coord=False)
            introns.append(Intron(self.genome, self.db, rank, self.StableId,
                                    coord))
            rank += 1
        
        self._cached['Introns'] = tuple(introns)
    
    def _get_introns(self):
        return self._get_cached_value('Introns', 
                        self._get_intron_transcript_records)
    
    Introns = property(_get_introns)
    
    def _get_translation_record(self):
        transcript_id = self.transcript_id
        translation_table = self.db.getTable('translation')
        query = sql.select([translation_table],
                        translation_table.c.transcript_id == transcript_id)
        try:
            record = asserted_one(query.execute())
        except NoItemError:
            self._set_null_values(['TranslatedExons'], 'translation')
            return
        
        self._table_rows['translation'] = record
    
    def _get_transcript(self):
        self._get_translation_record()
        record = self._table_rows['translation']
        if record == self.NULL_VALUE:
            return
        
        start_exon_id = record['start_exon_id']
        end_exon_id = record['end_exon_id']
        # because this will be used to shift a coord. Note: these are relative
        # to the exon start but ignore strand, so we have to decide whether
        # the coord shifts need to be flipped
        seq_start = record['seq_start'] - 1
        seq_end = record['seq_end']
        flip_coords = self.Exons[0].Location.Strand == -1
        
        start_index = None
        end_index = None
        for index, exon in enumerate(self.Exons):
            if exon.exon_id == start_exon_id:
                start_index = index
            if exon.exon_id == end_exon_id:
                end_index = index
        assert None not in (start_index, end_index), \
                                    'Error in matching transcript and exons'
        
        start_exon = self.Exons[start_index]
        
        if start_index == end_index:
            shift_start=[seq_start,len(start_exon)-seq_end][flip_coords]
            shift_end = [seq_end-len(start_exon), -1*seq_start][flip_coords]
        else:
            shift_start = [seq_start, 0][flip_coords]
            shift_end = [0, -1*seq_start][flip_coords]
        
        coord = start_exon.Location.resized(shift_start, shift_end)
        
        DEBUG = False
        if DEBUG:
            out = ['\nseq_start=%d; seq_end=%d' %(seq_start, seq_end),
                'shift_start=%d; shift_end=%d'%(shift_start, shift_end),
                'len=%s'% len(coord)]
            sys.stderr.write('\n'.join(map(str,out))+'\n')
        
        new_start_exon = Exon(self.genome, self.db, start_exon.exon_id,
                                start_exon.Rank, Location=coord)
        translated_exons=(new_start_exon,)+\
                        self.Exons[start_index+1:end_index]
        if start_index != end_index:
            end_exon = self.Exons[end_index]
            shift_start = [0, len(end_exon)-seq_end][flip_coords]
            shift_end = [seq_end-len(end_exon), 0][flip_coords]
            coord=end_exon.Location.resized(shift_start, shift_end)
            new_end_exon = Exon(self.genome, self.db, end_exon.exon_id,
                                end_exon.Rank, Location=coord)
            translated_exons += (new_end_exon,)
        self._cached['TranslatedExons'] = translated_exons
    
    def _get_translated_exons(self):
        return self._get_cached_value('TranslatedExons', self._get_transcript)
    
    TranslatedExons = property(_get_translated_exons)
    
    def _calculate_Utr_exons(self):
        # TODO clean up this code
        exons = self.Exons
        translated_exons = self.TranslatedExons
        num_exons = len(self.Exons)
        if not translated_exons:
            self._set_null_values(["UntranslatedExons5", "UntranslatedExons3"])
            return
        untranslated_5exons, untranslated_3exons = [], []
        start_exon, end_exon = translated_exons[0], translated_exons[-1]
        flip_coords = start_exon.Location.Strand == -1
        
        for exon in exons[0:start_exon.Rank]:   # get 5'UTR
            coord = exon.Location.copy()
            if exon.StableId == start_exon.StableId:
                coord.Start = [coord.Start, start_exon.Location.End][flip_coords]
                coord.End = [start_exon.Location.Start, coord.End][flip_coords]
            if len(coord) != 0:
                untranslated_5exons.append(Exon(self.genome, self.db,
                                exon.exon_id, exon.Rank, Location = coord))
        for exon in exons[end_exon.Rank -1: num_exons]: # get 3'UTR
            coord = exon.Location.copy()
            if exon.StableId == end_exon.StableId:
                coord.Start = [end_exon.Location.End, coord.Start][flip_coords]
                coord.End = [coord.End, end_exon.Location.Start][flip_coords]
            if len(coord) != 0:
                untranslated_3exons.append(Exon(self.genome, self.db,
                                exon.exon_id, exon.Rank,  Location = coord))
        
        self._cached["UntranslatedExons5"] = tuple(untranslated_5exons)
        self._cached["UntranslatedExons3"] = tuple(untranslated_3exons)
    
    def _get_5prime_untranslated_exons(self):
        return self._get_cached_value("UntranslatedExons5",
                                        self._calculate_Utr_exons)
    
    UntranslatedExons5 = property(_get_5prime_untranslated_exons)
    
    def _get_3prime_untranslated_exons(self):
        return self._get_cached_value("UntranslatedExons3",
                                        self._calculate_Utr_exons)
    
    UntranslatedExons3 = property(_get_3prime_untranslated_exons)
    
    def _make_utr_seq(self):
        if self.UntranslatedExons5 is None and self.UntranslatedExons3 is None:
            self._cached["Utr5"] = self.NULL_VALUE
            self._cached["Utr3"] = self.NULL_VALUE
            return
        Utr5_seq, Utr3_seq = DNA.makeSequence(""), DNA.makeSequence("")
        for exon in self.UntranslatedExons5:
            Utr5_seq += exon.Seq
        for exon in self.UntranslatedExons3:
            Utr3_seq += exon.Seq
        self._cached["Utr5"] = Utr5_seq
        self._cached["Utr3"] = Utr3_seq
    
    def _get_utr5_seq(self):
        return self._get_cached_value("Utr5", self._make_utr_seq)
    
    Utr5 = property(_get_utr5_seq)
    
    def _get_utr3_seq(self):
        return self._get_cached_value("Utr3", self._make_utr_seq)
    
    Utr3 = property(_get_utr3_seq)
    
    def _make_cds_seq(self):
        if self.Exons is self.NULL_VALUE:
            self._cached['Cds'] = self.NULL_VALUE
            return
        
        exons = [self.Exons, self.TranslatedExons][self._am_prot_coding]
        full_seq = None
        for exon in exons:
            if full_seq is None:
                full_seq = exon.Seq
                continue
            full_seq += exon.Seq
        
        # check first exon PhaseStart is 0 and last exon PhaseEnd
        if exons[0].PhaseStart > 0:
            fill = DNA.makeSequence('N' * exons[0].PhaseStart, Name=full_seq.Name)
            full_seq = fill + full_seq
        
        if exons[-1].PhaseEnd > 0:
            fill = DNA.makeSequence('N' * exons[-1].PhaseEnd, Name=full_seq.Name)
            full_seq += fill
        
        self._cached['Cds'] = full_seq
    
    def _get_cds(self):
        return self._get_cached_value('Cds', self._make_cds_seq)
    
    Cds = property(_get_cds)
    
    def getCdsLength(self):
        """returns the length of the Cds. If this property is not available,
        returns None."""
        if self.Cds is self.NULL_VALUE:
            return None
        exons = [self.Exons, self.TranslatedExons][self._am_prot_coding]
        return sum(map(len, exons))
    
    def _make_protein_seq(self):
        if not self._am_prot_coding or self.Cds is self.NULL_VALUE:
            self._cached['ProteinSeq'] = self.NULL_VALUE
            return
        
        DEBUG = False
        # enforce multiple of 3
        cds = self.Cds
        length = len(cds)
        cds = cds[: length - (length % 3)]
        try:
            cds = cds.withoutTerminalStopCodon()
        except AssertionError:
            if not DEBUG:
                raise
            out = ['\n****\nFAILED=%s' % self.StableId]
            for exon in self.TranslatedExons:
                out += ['TranslatedExon[rank=%d]\n' % exon.Rank, exon,
                        exon.Location, 
                        '%s ... %s'%(exon.Seq[:20],exon.Seq[-20:])]
                sys.stderr.write('\n'.join(map(str, out))+'\n')
            raise
        
        self._cached['ProteinSeq'] = cds.getTranslation()
    
    def _get_protein_seq(self):
        return self._get_cached_value('ProteinSeq', self._make_protein_seq)
    
    ProteinSeq = property(_get_protein_seq)
    
    def _get_exon_feature_data(self, parent_map):
        """returns the exon feature data"""
        features = []
        if self.Exons is self.NULL_VALUE:
            return features
        for exon in self.Exons:
            feature_data = exon.featureData(parent_map)
            if feature_data is None:
                continue
            features.append(feature_data)
        return features
    
    def _get_intron_feature_data(self, parent_map):
        """return the intron feature data"""
        features = []
        if self.Introns is self.NULL_VALUE:
            return features
        for intron in self.Introns:
            feature_data = intron.featureData(parent_map)
            if feature_data is None:
                continue
            features.append(feature_data)
        return features
    
    def _get_translated_exon_feature_data(self, parent_map):
        """returns featureD data for translated exons"""
        features = []
        if self.TranslatedExons is self.NULL_VALUE:
            return features
        cds_spans = []
        for exon in self.TranslatedExons:
            feature_data = exon.featureData(parent_map)
            if feature_data is None:
                continue
            cds_spans.extend(feature_data[-1].spans)
        if cds_spans:
            # TODO: check strand
            cds_map = Map(spans = cds_spans, parent_length = parent_map.parent_length)
            features.append(('CDS', str(self.StableId), cds_map))
        return features
    
    def _get_Utr_feature_data(self, parent_map):
        # TODO: Simplify this part
        features = []
        utr5_spans, utr3_spans = [], []
        for exon in self.UntranslatedExons5:
            feature_data = exon.featureData(parent_map)
            if feature_data is None:
                continue
            utr5_spans.extend(feature_data[-1].spans)
        for exon in self.UntranslatedExons3:
            feature_data = exon.featureData(parent_map)
            if feature_data is None:
                continue
            utr3_spans.extend(feature_data[-1].spans)
        if utr5_spans:
            utr5_map = Map(spans = utr5_spans,
                                parent_length = parent_map.parent_length)
            features.append(("5'UTR", str(self.StableId), utr5_map))
        if utr3_spans:
            utr3_map = Map(spans = utr3_spans,
                                parent_length = parent_map.parent_length)
            features.append(("3'UTR", str(self.StableId), utr3_map))
        return features
    
    def subFeatureData(self, parent_map):
        """returns data for making a cogent Feature. This can be automatically
        applied to the Seq by the getAnnotatedSeq method. Returns None if
        self lies outside parent's span.
        """
        features = self._get_exon_feature_data(parent_map)
        features += self._get_intron_feature_data(parent_map)
        features += self._get_translated_exon_feature_data(parent_map)
        if self.TranslatedExons:
            features += self._get_Utr_feature_data(parent_map)
        return features
    

class Exon(_StableRegion):
    Type = 'exon'
    
    def __init__(self, genome, db, exon_id, Rank, Location=None):
        """created by a Gene"""
        _StableRegion.__init__(self, genome, db, Location=Location)

        self._attr_ensembl_table_map = dict(StableId=['exon_stable_id',
                                                      'exon'][genome.GeneralRelease >= 65],
                                            Location='exon')

        self.exon_id = exon_id
        self.Rank = Rank
    
    def __str__(self):
        return self.__repr__()
    
    def __repr__(self):
        my_type = self.__class__.__name__
        return '%s(StableId=%s, Rank=%s)' % (my_type, self.StableId, self.Rank)
    
    def __cmp__(self, other):
        return cmp(self.Rank, other.Rank)
    
    def _get_exon_stable_id_record(self):
        if self.genome.GeneralRelease >= 65:
            # release >= 65, data is just in the exon table
            self._get_exon_record()
            return
        
        table_name = self._attr_ensembl_table_map['StableId']
        exon_stable_id_table = self.db.getTable(table_name)
        query = sql.select([exon_stable_id_table.c.stable_id],
                           exon_stable_id_table.c.exon_id == self.exon_id)
        records = query.execute()
        record = asserted_one(records.fetchall())
        self._table_rows[table_name] = record
    
    def _get_exon_record(self):
        # this will be called by _Region parent class to make the location
        exon_table = self.db.getTable('exon')
        query = sql.select([exon_table], exon_table.c.exon_id == self.exon_id)
        records = query.execute()
        record = asserted_one(records.fetchall())
        self._table_rows['exon'] = record
    
    def _make_symbol(self):
        self._cached['Symbol'] = '%s-%s' % (self.StableId, self.Rank)
    
    def _get_symbol(self):
        return self._get_cached_value('Symbol', self._make_symbol)
    
    Symbol = property(_get_symbol)
    
    def _make_phase(self):
        """creates the exon phase attributes"""
        if 'exon' not in self._table_rows:
            self._get_exon_record()
        
        exon = self._table_rows['exon']
        self._cached['PhaseStart'] = exon['phase']
        self._cached['PhaseEnd'] = exon['end_phase']
    
    @property
    def PhaseStart(self):
        """reading frame start for this exon"""
        return self._get_cached_value('PhaseStart', self._make_phase)
    
    @property
    def PhaseEnd(self):
        """reading frame end for this exon"""
        return self._get_cached_value('PhaseEnd', self._make_phase)
    
    


class Intron(GenericRegion):
    Type = 'intron'
    def __init__(self, genome, db, rank, transcript_stable_id, Location=None):
        GenericRegion.__init__(self, genome, db, Location=Location)
        self.TranscriptStableId = transcript_stable_id
        self.Rank = rank
    
    def __str__(self):
        return self.__repr__()
    
    def __repr__(self):
        my_type = self.__class__.__name__
        return '%s(TranscriptId=%s, Rank=%s)' % (my_type,
                            self.TranscriptStableId, self.Rank)
    
    def _make_symbol(self):
        self._cached['Symbol'] = '%s-%s'%(self.TranscriptStableId, self.Rank)
    
    def _get_symbol(self):
        return self._get_cached_value('Symbol', self._make_symbol)
    
    Symbol = property(_get_symbol)


class Est(Gene):
    """an EST region"""
    Type = 'est'


def _set_to_string(val):
    if type(val) in (str, type(None)):
        return val
    val = list(val)
    while len(val) == 1 and type(val) in (tuple, list):
        val = val[0]
    return val

class Variation(_Region):
    """genomic variation"""
    Type = 'variation'
    
    def __init__(self, genome, db=None, Effect=None, Symbol=None, data=None):
        self.genome = genome
        
        get_table = genome.VarDb.getTable
        self.variation_feature_table = get_table('variation_feature')
        self.transcript_variation_table = get_table('transcript_variation')
        self.flanking_sequence_table = get_table('flanking_sequence')
        self.allele_table = get_table('allele')
        try:
            self.allele_code_table = get_table('allele_code')
        except sql.exceptions.ProgrammingError:
            self.allele_code_table = None
        
        super(Variation, self).__init__()
        
        self._attr_ensembl_table_map = dict(Effect='variation_feature',
                                            Symbol='variation_feature',
                                            Validation='variation_feature',
                                            MapWeight='variation_feature',
                                            FlankingSeq='flanking_sequence',
                                            PeptideAlleles='transcript_variation',
                                            TranslationLocation='transcript_variation',
                                            Location='variation_feature',
                                            AlleleFreqs='allele')
        
        assert data is not None, 'Variation record created in an unusual way'
        for name, value, func in \
                [('Effect',Effect,self._get_variation_table_record),
                 ('Symbol',Symbol,self._get_variation_table_record)]:
            if value is not None:
                self._table_rows[self._attr_ensembl_table_map[name]] = data
                if func is not None:
                    func() # this populates the attributes
        
        self.db = db or self.genome.CoreDb
    
    def __len__(self):
        """return the length of the longest allelic variant"""
        return max(map(len, self._split_alleles()))
    
    def __str__(self):
        my_type = self.__class__.__name__
        
        return "%s(Symbol=%r; Effect=%r; Alleles=%r)" % \
            (my_type, self.Symbol, self.Effect, self.Alleles)
    
    def _get_variation_table_record(self):
        consequence_type = 'consequence_type'
        if self.genome.GeneralRelease > 67:
            consequence_type += 's' # change to plural column name
        
        attr_name_map = [('Effect', consequence_type, _set_to_string),
                         ('Alleles', 'allele_string', _quoted),
                         ('Symbol', 'variation_name', _quoted),
                         ('Validation', 'validation_status', _set_to_string),
                         ('MapWeight', 'map_weight', int)]
        self._populate_cache_from_record(attr_name_map, 'variation_feature')
        # TODO handle obtaining the variation_feature if we were created in
        # any way other than through the Symbol or Effect
    
    def _get_seq_region_record(self, seq_region_id):
        # should this be on a parent class? or a generic function in assembly?
        seq_region_table = self.db.getTable('seq_region')
        query = sql.select([seq_region_table],
                        seq_region_table.c.seq_region_id == seq_region_id)
        record = asserted_one(query.execute())
        return record
    
    def _get_flanking_seq_data(self):
        # maps to flanking_sequence through variation_feature_id
        # if this fails, we grab from genomic sequence
        variation_id = self._table_rows['variation_feature']['variation_id']
        flanking_seq_table = self.flanking_sequence_table
        query = sql.select([flanking_seq_table],
                    flanking_seq_table.c.variation_id == variation_id)
        record = asserted_one(query.execute())
        self._table_rows['flanking_sequence'] = record
        up_seq = record['up_seq']
        down_seq = record['down_seq']
        # the following two lines are because -- wait for it -- someone has
        # entered the string 'NULL' instead of NULL in the MySQL tables!!!
        up_seq = [up_seq, None][up_seq == 'NULL']
        down_seq = [down_seq, None][down_seq == 'NULL']
        seqs = dict(up=up_seq, down=down_seq)
        for name, seq in seqs.items():
            if seq is not None:
                seq = DNA.makeSequence(seq)
            else:
                resized = [(-301, -1), (1, 301)][name == 'down']
                if self.Location.Strand == -1:
                    resized = [(1, 301), (-301, -1)][name == 'down']
                flank = self.Location.resized(*resized)
                flanking = self.genome.getRegion(region=flank)
                seq = flanking.Seq
            seqs[name] = seq
        
        self._cached[('FlankingSeq')] = (seqs['up'][-300:],seqs['down'][:300])
    
    def _get_flanking_seq(self):
        return self._get_cached_value('FlankingSeq',
                                self._get_flanking_seq_data)
    
    FlankingSeq = property(_get_flanking_seq)
    
    def _get_effect(self):
        return self._get_cached_value('Effect',
                        self._get_variation_table_record)
    
    Effect = property(_get_effect)
    
    def _get_alleles(self):
        return self._get_cached_value('Alleles',
                        self._get_variation_table_record)
    
    Alleles = property(_get_alleles)
    
    def _get_allele_table_record(self):
        variation_id = self._table_rows['variation_feature']['variation_id']
        allele_table = self.allele_table
        query = sql.select([allele_table],
                allele_table.c.variation_id == variation_id)
        records = [r for r in query.execute()]
        
        if len(records) == 0:
            self._cached[('AlleleFreqs')] = self.NULL_VALUE
            return
        
        # property change from >= 65, allele ids need to be looked up in
        # the allele_code table
        allele_code = self.allele_code_table
        
        self._table_rows['allele_table'] = records
        data = []
        for rec in records:
            if not rec['sample_id']:
                continue
            
            if allele_code is None:
                allele = rec['allele']
            else:
                allele_query = sql.select([allele_code.c.allele],
                                allele_code.c.allele_code_id == rec['allele_code_id'])
                allele = list(allele_query.execute())[0][0]
            
            data.append((allele, rec['frequency'], rec['sample_id']))
        
        if not data:
            self._cached[('AlleleFreqs')] = self.NULL_VALUE
            return
        
        table = Table(header='allele freq sample_id'.split(), rows=data)
        self._cached[('AlleleFreqs')] = table.sorted(['sample_id', 'allele'])
    
    def _get_allele_freqs(self):
        return self._get_cached_value('AlleleFreqs',
                        self._get_allele_table_record)
    
    AlleleFreqs = property(_get_allele_freqs)
    
    def _get_symbol(self):
        return self._get_cached_value('Symbol',
                        self._get_variation_table_record)
    
    Symbol = property(_get_symbol)
    
    def _get_validation(self):
        return self._get_cached_value('Validation',
                        self._get_variation_table_record)
    
    Validation = property(_get_validation)
    
    def _get_map_weight(self):
        return self._get_cached_value('MapWeight',
                        self._get_variation_table_record)
    
    MapWeight = property(_get_map_weight)
    
    def _get_transcript_record(self):
        if not 'variation_feature' in self._table_rows:
            raise NotImplementedError
        
        try:
            effects = [self.Effect.lower()]
        except AttributeError:
            effects = [v.lower() for v in self.Effect]
        
        effects = set(effects)
        nsyn = set(('non_synonymous_coding', 'non_synonymous_codon'))
        if not effects & nsyn:
            self._cached['PeptideAlleles'] = self.NULL_VALUE
            self._cached['TranslationLocation'] = self.NULL_VALUE
            return
        
        table_name = self._attr_ensembl_table_map['PeptideAlleles']
        loc = lambda x: int(x)-1
        
        # column name changed between releases, so we check to see which
        # one is being used for this instance and set the column strings
        
        # TODO can we modify the table on loading? This would give better
        # performance.
        
        if self.genome.VarDb.tableHasColumn(table_name, 'pep_allele_string'):
            pep_allele_string = 'pep_allele_string'
            consequence_type = 'consequence_types'
        else:
            pep_allele_string = 'peptide_allele_string'
            consequence_type = 'consequence_type'
        
        attr_column_map = [
                ('PeptideAlleles', pep_allele_string, _quoted),
                ('TranslationLocation', 'translation_start', loc)]
        
        if table_name in self._table_rows:
            self._populate_cache_from_record(attr_column_map, table_name)
            return
        
        var_feature_record = self._table_rows['variation_feature']
        var_feature_id = var_feature_record['variation_feature_id']
        table = self.genome.VarDb.getTable(table_name)
        self_effect = set([self.Effect,[self.Effect]][type(self.Effect)==str])
        query = sql.select([table.c.variation_feature_id,
                           table.columns[pep_allele_string],
                           table.c.translation_start,
                           table.columns[consequence_type]],
                    sql.and_(table.c.variation_feature_id == var_feature_id,
                            table.columns[pep_allele_string] != None))
        records = query.execute().fetchall()
        pep_alleles = []
        translation_location = []
        for record in records:
            if not record[consequence_type] & self_effect:
                continue
            pep_alleles += [record[pep_allele_string]]
            translation_location += [record['translation_start']]
        
        if not pep_alleles:
            print 'Expected at least a single record'
            raise RuntimeError
        
        # we only want unique allele strings
        allele_location = dict(zip(pep_alleles, translation_location))
        pep_alleles = list(set(pep_alleles))
        pep_alleles = [pep_alleles, pep_alleles[0]][len(pep_alleles)==1]
        if type(pep_alleles) != str:
            for pep_allele in pep_alleles:
                translation_location = allele_location[pep_allele]
        else:
            translation_location = allele_location[pep_alleles]
        
        self._table_rows[table_name] = dict(pep_allele_string=pep_alleles,
                                            translation_start=translation_location)
        self._populate_cache_from_record(attr_column_map, table_name)
    
    def _get_peptide_variation(self):
        return self._get_cached_value('PeptideAlleles',
                                self._get_transcript_record)
    
    PeptideAlleles = property(_get_peptide_variation)
    
    def _get_translation_location(self):
        return self._get_cached_value('TranslationLocation',
                                self._get_transcript_record)
    
    TranslationLocation = property(_get_translation_location)
    
    def _split_alleles(self):
        return self.Alleles.split('/')
    
    def _get_number_alleles(self):
        result = self._split_alleles()
        return len(result)
    
    NumAlleles = property(_get_number_alleles)

class CpGisland(GenericRegion):
    Type = 'CpGisland'
    
    def __init__(self, genome, db, Location, Score):
        super(CpGisland, self).__init__(genome=genome, db=db, 
                                        Location=Location)
        self.Score = Score
    
    def __str__(self):
        my_type = self.__class__.__name__
        
        return "%s(CoordName='%s'; Start=%s; End=%s; length=%s;"\
            " Strand='%s', Score=%.1f)" % (my_type,
            self.Location.CoordName,
            self.Location.Start,
            self.Location.End,
            len(self),
            '-+'[self.Location.Strand>0], self.Score)
    

class Repeat(GenericRegion):
    Type = 'repeat'
    
    def __init__(self, genome, db, Location, Score, data):
        super(Repeat, self).__init__(genome=genome, db=db, Location=Location)
        self._attr_ensembl_table_map = dict(Symbol='repeat_consensus',
                                            RepeatType='repeat_consensus',
                                            RepeatClass='repeat_consensus',
                                            Consensus='repeat_consensus')

        self.Score = Score
        # assume always created from repeat_feature table
        self._table_rows['repeat_feature']= data
    
    def __str__(self):
        my_type = self.__class__.__name__
        
        return "%s(CoordName='%s'; Start=%s; End=%s; length=%s;"\
        " Strand='%s', Score=%.1f)" %  (my_type,
            self.Location.CoordName,
            self.Location.Start, self.Location.End, len(self),
            '-+'[self.Location.Strand>0], self.Score)
    
    def _get_repeat_consensus_record(self):
        repeat_consensus_table = self.db.getTable('repeat_consensus')
        repeat_consensus_id = self._table_rows['repeat_feature']['repeat_consensus_id']
        record = sql.select([repeat_consensus_table],
        repeat_consensus_table.c.repeat_consensus_id == repeat_consensus_id)
        record = asserted_one(record.execute().fetchall())
        self._table_rows['repeat_consensus'] = record
        limit_length = lambda x : DisplayString(x, repr_length=10)
        attr_column_map = [('Symbol', 'repeat_name', _quoted),
                           ('RepeatClass', 'repeat_class', _quoted),
                           ('RepeatType', 'repeat_type', _quoted),
                           ('Consensus', 'repeat_consensus', limit_length)]
        self._populate_cache_from_record(attr_column_map, 'repeat_consensus')
    
    def _get_symbol(self):
        return self._get_cached_value('Symbol',
            self._get_repeat_consensus_record)
    
    Symbol = property(_get_symbol)
    
    def _get_repeat_class(self):
        return self._get_cached_value('RepeatClass',
            self._get_repeat_consensus_record)
    
    RepeatClass = property(_get_repeat_class)
    
    def _get_repeat_type(self):
        return self._get_cached_value('RepeatType',
            self._get_repeat_consensus_record)
    
    RepeatType = property(_get_repeat_type)
    
    def _get_consensus(self):
        return self._get_cached_value('Consensus',
            self._get_repeat_consensus_record)
    
    Consensus = property(_get_consensus)
    
