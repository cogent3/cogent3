from pprint import pprint
import sqlalchemy as sql

from cogent import DNA
from cogent.core.alignment import SequenceCollection, Alignment, Aligned
from cogent.parse import cigar 

from cogent.db.ensembl.util import LazyRecord, asserted_one, NoItemError
from cogent.db.ensembl.assembly import location_query
from cogent.db.ensembl.species import Species

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

class _RelatedRegions(LazyRecord):
    # a basic related region, capable of providing the sequences
    # obtaining SyntenicRegions -- for getting aligned blocks -- is delegated
    # to compara
    Type = None
    def __init__(self):
        super(_RelatedRegions, self).__init__()
    
    def __str__(self):
        # temporary string method, just to demo correct assembly
        # TODO StableID, Species and Description
        my_type = self.__class__.__name__
        
        data = map(repr, self.Members)
        data.insert(0, '%s(' % my_type)
        data.append(')')
        return "\n\t".join(data)
    
    def getSeqCollection(self, feature_types=None, where_feature=None):
        """returns a SequenceCollection instance of the unaligned sequences"""
        seqs = []
        for member in self.Members:
            if feature_types:
                seq = member.getAnnotatedSeq(feature_types,where_feature)
            else:
                seq = member.Seq
            if seq is None:
                continue
            seqs.append((seq.Name, seq))
        return SequenceCollection(data=seqs, MolType=DNA)
    
    def getSeqLengths(self):
        """returns a vector of lengths"""
        return [len(member) for member in self.Members]
    
    def getSpeciesSet(self):
        """returns the latin names of self.Member species as a set"""
        return set([m.Location.Species for m in self.Members])
    

class RelatedGenes(_RelatedRegions):
    Type = 'related_genes'
    def __init__(self, compara, Members, Relationships):
        super(RelatedGenes, self).__init__()
        self.compara = compara
        self.Members = tuple(Members)
        self.Relationships = Relationships
    
    def __str__(self):
        my_type = self.__class__.__name__
        
        display = ['%s:' % my_type,
                   ' Relationships=%s' % self.Relationships]
        display += ['  %s' % m for m in self.Members]
        return '\n'.join(display)
    
    def __repr__(self):
        return self.__str__()
    
    def getMaxCdsLengths(self):
        """returns the vector of maximum Cds lengths from member transcripts"""
        return [max(member.getCdsLengths()) for member in self.Members]


class SyntenicRegion(LazyRecord):
    """a class that takes the genome, compara instances and is used to build
    Aligned sequences for Ensembl multiple alignments"""
    def __init__(self, parent, genome, identifiers_values, am_ref_member,
                    Location=None):
        # create with method_link_species_set_id, at least, in
        # identifiers_values
        super(SyntenicRegion, self).__init__()
        self.parent = parent
        self.compara = parent.compara
        self.genome = genome
        self.am_ref_member = am_ref_member
        self.aln_map = None
        self.aln_loc = None
        self._make_map_func = [self._make_map_from_ref,
                                  self._make_ref_map][am_ref_member]
        
        if Location is not None:
            if hasattr(Location, 'Location'): # likely to be a feature region
                region = Location
            else:
                region = genome.getRegion(region=Location)
            self._cached['Region'] = region
        
        for identifier, value in dict(identifiers_values).items():
            self._cached[identifier] = value
        
    
    def __len__(self):
        return len(self._get_cached_value('Region', self._make_map_func))
    
    def _get_location(self):
        region = self._get_cached_value('Region', self._make_map_func)
        return region.Location
    
    Location = property(_get_location)
    
    def _get_region(self):
        region = self._get_cached_value('Region', self._make_map_func)
        return region
    
    Region = property(_get_region)
    
    def _get_cigar_record(self):
        genomic_align_table = \
                self.parent.compara.ComparaDb.getTable('genomic_align')
        query = sql.select([genomic_align_table.c.cigar_line],
                    genomic_align_table.c.genomic_align_id == \
                                        self._cached['genomic_align_id'])
        record = asserted_one(query.execute())
        self._cached['cigar_line'] = record['cigar_line']
        return record
    
    def _get_cigar_line(self):
        return self._get_cached_value('cigar_line', self._get_cigar_record)
    
    cigar_line = property(_get_cigar_line)
    
    def _make_ref_map(self):
        if self.aln_map and self.aln_loc is not None:
            return
        
        ref_record = self._cached
        record_start = ref_record['dnafrag_start']
        record_end = ref_record['dnafrag_end']
        record_strand = ref_record['dnafrag_strand']
        
        block_loc = self.genome.makeLocation(CoordName=ref_record['name'],
                                       Start=record_start,
                                       End=record_end,
                                       Strand=record_strand,
                                       ensembl_coord=True)
        
        ref_location = self.parent.ref_location
        relative_start = ref_location.Start-block_loc.Start
        relative_end = relative_start + len(ref_location)
        if block_loc.Strand != 1:
            relative_start = len(block_loc) - relative_end
            relative_end = relative_start + len(ref_location)
        
        aln_map, aln_loc = cigar.slice_cigar(self.cigar_line, relative_start,
                                         relative_end, by_align = False)
        
        self.aln_map = aln_map
        self.aln_loc = aln_loc
        region_loc = ref_location.copy()
        region_loc.Strand = block_loc.Strand
        region = self.genome.getRegion(region=region_loc)
        self._cached['Region'] = region
    
    def _make_map_from_ref(self):
        # this is the 'other' species
        if self.aln_loc and self.aln_map is not None:
            return
        record = self._cached
        try:
            aln_map, aln_loc = cigar.slice_cigar(self.cigar_line,
                                            self.parent.CigarStart,
                                            self.parent.CigarEnd,
                                            by_align=True)
            self.aln_map = aln_map
            self.aln_loc = aln_loc # probably unnecesary to store??
            
            # we make a loc for the aligned region
            block_loc = self.genome.makeLocation(CoordName=record['name'],
                                             Start=record['dnafrag_start'],
                                             End = record['dnafrag_end'],
                                             Strand=record['dnafrag_strand'],
                                             ensembl_coord=True)
            relative_start = aln_loc[0]
            relative_end = aln_loc[1]
            # new location with correct length
            loc = block_loc.copy()
            loc.End = loc.Start+(relative_end-relative_start)
            
            if block_loc.Strand != 1:
                shift = len(block_loc) - relative_end
            else:
                shift = relative_start
            loc = loc.shifted(shift)
            region = self.genome.getRegion(region=loc)
        except IndexError: # TODO ask Hua where these index errors occur
            region = None
        self._cached['Region'] = region
    
    def _make_aligned(self, feature_types = None, where_feature=None):
        if self.aln_loc is None or self.aln_map is None: # is this required?
            self._make_map_func()
        region = self._cached['Region']
        if region is None:
            self._cached['AlignedSeq'] = None
            return
        if feature_types:
            seq = region.getAnnotatedSeq(feature_types, where_feature)
        else:
            seq = region.Seq
        
        # we get the Seq objects to allow for copying of their annotations
        gapped_seq = Aligned(self.aln_map, seq)
        
        self._cached['AlignedSeq'] = gapped_seq
    
    def _get_aligned_seq(self):
        aligned = self._get_cached_value('AlignedSeq', self._make_aligned)
        return aligned
    
    AlignedSeq = property(_get_aligned_seq)
    
    def getAnnotatedAligned(self, feature_types, where_feature=None):
        """returns aligned seq annotated for the specified feature types"""
        region = self._get_cached_value('Region', self._make_map_func)
        if region is None:
            return None
        self._make_aligned(feature_types=feature_types,
                                    where_feature=where_feature)
        return self.AlignedSeq
    

class SyntenicRegions(_RelatedRegions):
    Type = 'syntenic_regions'
    def __init__(self, compara, Members, ref_location):
        super(SyntenicRegions, self).__init__()
        self.compara = compara
        members = []
        ref_member = None
        self.ref_location = ref_location
        for genome, data in Members:
            if genome is ref_location.genome:
                ref_member = SyntenicRegion(self, genome, dict(data),
                                    am_ref_member=True, Location=ref_location)
            else:
                members += [SyntenicRegion(self, genome, dict(data),
                                                    am_ref_member=False)]
        
        assert ref_member is not None, "Can't match a member to ref_location"
        self.ref_member = ref_member
        self.Members = tuple([ref_member] + members)
        self.NumMembers = len(self.Members)
        self.aln_loc = None
        self._do_rc = None
    
    def __str__(self):
        my_type = self.__class__.__name__
        
        display = ['%s:' % my_type]
        display += ['  %r' % m.Location for m in self.Members \
                                                if m.Region is not None]
        return '\n'.join(display)
    
    def __repr__(self):
        return self.__str__()
    
    def _populate_ref(self):
        """near (don't actually get the sequence) completes construction of
        ref sequence"""
        self.ref_member._make_map_func()
        self._cached['CigarStart'] = self.ref_member.aln_loc[0]
        self._cached['CigarEnd'] = self.ref_member.aln_loc[1]
    
    def _get_rc_state(self):
        """determines whether the ref_member strand is the same as that from
        the align block, if they diff we will rc the alignment, seqs,
        seq_names"""
        if self._do_rc is not None:
            return self._do_rc
        self._populate_ref()
        inferred = self.ref_member._cached['Region'].Location.Strand
        self._do_rc = self.ref_location.Strand != inferred
        return self._do_rc
    
    _rc = property(fget=_get_rc_state)
    
    def __len__(self):
        return self.CigarEnd - self.CigarStart
    
    def _get_ref_start(self):
        return self._get_cached_value('CigarStart', self._populate_ref)
    
    CigarStart = property(_get_ref_start)
    
    def _get_ref_end(self):
        return self._get_cached_value('CigarEnd', self._populate_ref)
    
    CigarEnd = property(_get_ref_end)
    
    def getAlignment(self, feature_types=None, where_feature=None,
                        omit_redundant=True):
        """Arguments:
            - feature_types: annotations to be applied to the returned
              sequences
            - omit_redundant: exclude redundant gap positions"""
        seqs = []
        annotations = {}
        
        for member in self.Members:
            if feature_types:
                seq = member.getAnnotatedAligned(feature_types, where_feature)
            else:
                seq = member.AlignedSeq
            if seq is None:
                continue
            name = seq.Name
            
            if self._rc: # names should reflect change to strand
                loc = member.Location.copy()
                loc.Strand *= -1
                name = str(loc)
            
            annotations[name] = seq.data.annotations
            seq.Name = seq.data.Name = name
            seqs += [(name, seq)]
        
        if seqs is None:
            return None
        
        aln = Alignment(data=seqs, MolType=DNA)
        
        if self._rc:
            aln = aln.rc()
        
        if omit_redundant:
            aln = aln.filtered(lambda x: set(x) != set('-'))
        
        return aln
    

