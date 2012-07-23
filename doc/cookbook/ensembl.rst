Note that much more extensive documentation is available in :ref:`query-ensembl`.

Connecting
----------

.. Gavin Huttley

`Ensembl <http://www.ensembl.org>`_ provides access to their MySQL databases directly or users can download and run those databases on a local machine. To use the Ensembl's UK servers for running queries, nothing special needs to be done as this is the default setting for PyCogent's ``ensembl`` module. To use a different Ensembl installation, you create an account instance:

.. doctest::

    >>> from cogent.db.ensembl import HostAccount
    >>> account = HostAccount('fastcomputer.topuni.edu', 'username',
    ...                       'canthackthis')

To specify a specific port to connect to MySQL on:

.. doctest::

    >>> from cogent.db.ensembl import HostAccount
    >>> account = HostAccount('anensembl.server.edu', 'someuser',
    ...                       'somepass', port=3306)

.. we create valid account now to work on my local machines here at ANU

.. doctest::
    :hide:

    >>> import os
    >>> hotsname, uname, passwd = os.environ['ENSEMBL_ACCOUNT'].split()
    >>> account = HostAccount(hotsname, uname, passwd)

Species to be queried
---------------------

To see what existing species are available

.. doctest::

    >>> from cogent.db.ensembl import Species
    >>> print Species
    ================================================================================
           Common Name                   Species Name              Ensembl Db Prefix
    --------------------------------------------------------------------------------
             A.aegypti                  Aedes aegypti                  aedes_aegypti
            A.clavatus           Aspergillus clavatus           aspergillus_clavatus...

If Ensembl has added a new species which is not yet included in ``Species``, you can add it yourself.

.. doctest::

    >>> Species.amendSpecies('A latinname', 'a common name')

You can get the common name for a species

.. doctest::

    >>> Species.getCommonName('Procavia capensis')
    'Rock hyrax'

and the Ensembl database name prefix which will be used for all databases for this species.

.. doctest::

    >>> Species.getEnsemblDbPrefix('Procavia capensis')
    'procavia_capensis'

Species common names are used to construct attributes on PyCogent ``Compara`` instances). You can get the name that will be using the ``getComparaName`` method. For species with a real common name

.. doctest::
    
    >>> Species.getComparaName('Procavia capensis')
    'RockHyrax'

or with a shortened species name

.. doctest::
    
    >>> Species.getComparaName('Caenorhabditis remanei')
    'Cremanei'

Get genomic features
--------------------

Find a gene by gene symbol
^^^^^^^^^^^^^^^^^^^^^^^^^^

We query for the *BRCA2* gene for humans.

.. doctest::

    >>> from cogent.db.ensembl import Genome
    >>> human = Genome('human', Release=67, account=account)
    >>> print human
    Genome(Species='Homo sapiens'; Release='67')
    >>> genes = human.getGenesMatching(Symbol='BRCA2')
    >>> for gene in genes:
    ...     if gene.Symbol == 'BRCA2':
    ...         print gene
    ...         break
    Gene(Species='Homo sapiens'; BioType='protein_coding'; Description='breast cancer 2,...'; StableId='ENSG00000139618'; Status='KNOWN'; Symbol='BRCA2')

Find a gene by Ensembl Stable ID
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use the stable ID for *BRCA2*.

.. doctest::

    >>> from cogent.db.ensembl import Genome
    >>> human = Genome('human', Release=67, account=account)
    >>> gene = human.getGeneByStableId(StableId='ENSG00000139618')
    >>> print gene
    Gene(Species='Homo sapiens'; BioType='protein_coding'; Description='breast cancer 2,...'; StableId='ENSG00000139618'; Status='KNOWN'; Symbol='BRCA2')

Find genes matching a description
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We look for breast cancer related genes that are estrogen induced.

.. doctest::

    >>> from cogent.db.ensembl import Genome
    >>> human = Genome('human', Release=67, account=account)
    >>> genes = human.getGenesMatching(Description='breast cancer anti-estrogen')
    >>> for gene in genes:
    ...     print gene
    Gene(Species='Homo sapiens'; BioType='lincRNA'; Description='breast cancer anti-estrogen...'; StableId='ENSG00000262117'; Status='NOVEL'; Symbol='BCAR4')...

We can also require that an exact (case insensitive) match to the word(s) occurs within the description by setting ``like=False``.

.. doctest::
    
    >>> genes = human.getGenesMatching(Description='breast cancer anti-estrogen',
    ...                                  like=False)
    >>> for gene in genes:
    ...     print gene
    Gene(Species='Homo sapiens'; BioType='lincRNA'; Description='breast cancer anti-estrogen...'; StableId='ENSG00000262117'; Status='NOVEL'; Symbol='BCAR4')...

Get canonical transcript for a gene
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We get the canonical transcripts for *BRCA2*.

.. doctest::

    >>> from cogent.db.ensembl import Genome
    >>> human = Genome('human', Release=67, account=account)
    >>> brca2 = human.getGeneByStableId(StableId='ENSG00000139618')
    >>> transcript = brca2.CanonicalTranscript
    >>> print transcript
    Transcript(Species='Homo sapiens'; CoordName='13'; Start=32889610; End=32973347; length=83737; Strand='+')

Get the CDS for a transcript
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.db.ensembl import Genome
    >>> human = Genome('human', Release=67, account=account)
    >>> brca2 = human.getGeneByStableId(StableId='ENSG00000139618')
    >>> transcript = brca2.CanonicalTranscript
    >>> cds = transcript.Cds
    >>> print type(cds)
    <class 'cogent.core.sequence.DnaSequence'>
    >>> print cds
    ATGCCTATTGGATCCAAAGAGAGGCCA...

Look at all transcripts for a gene
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.db.ensembl import Genome
    >>> human = Genome('human', Release=67, account=account)
    >>> brca2 = human.getGeneByStableId(StableId='ENSG00000139618')
    >>> for transcript in brca2.Transcripts:
    ...     print transcript
    Transcript(Species='Homo sapiens'; CoordName='13'; Start=32889610; End=32973347; length=83737; Strand='+')
    Transcript(Species='Homo sapiens'; CoordName='13'; Start=32889641; End=32907428; length=17787; Strand='+')...

Get the first exon for a transcript
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We show just for the canonical transcript.

.. doctest::

    >>> from cogent.db.ensembl import Genome
    >>> human = Genome('human', Release=67, account=account)
    >>> brca2 = human.getGeneByStableId(StableId='ENSG00000139618')
    >>> print brca2.CanonicalTranscript.Exons[0]
    Exon(StableId=ENSE00001184784, Rank=1)

Get the introns for a transcript
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We show just for the canonical transcript.

.. doctest::

    >>> from cogent.db.ensembl import Genome
    >>> human = Genome('human', Release=67, account=account)
    >>> brca2 = human.getGeneByStableId(StableId='ENSG00000139618')
    >>> for intron in brca2.CanonicalTranscript.Introns:
    ...     print intron
    Intron(TranscriptId=ENST00000380152, Rank=1)
    Intron(TranscriptId=ENST00000380152, Rank=2)
    Intron(TranscriptId=ENST00000380152, Rank=3)...


Inspect the genomic coordinate for a feature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.db.ensembl import Genome
    >>> human = Genome('human', Release=67, account=account)
    >>> brca2 = human.getGeneByStableId(StableId='ENSG00000139618')
    >>> print brca2.Location.CoordName
    13
    >>> print brca2.Location.Start
    32889610
    >>> print brca2.Location.Strand
    1

Get repeat elements in a genomic interval
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We query the genome for repeats within a specific coordinate range on chromosome 13.

.. doctest::

    >>> from cogent.db.ensembl import Genome
    >>> human = Genome('human', Release=67, account=account)
    >>> repeats = human.getFeatures(CoordName='13', Start=32879610, End=32889610, feature_types='repeat')
    >>> for repeat in repeats:
    ...     print repeat.RepeatClass
    ...     print repeat
    ...     break
    SINE/Alu
    Repeat(CoordName='13'; Start=32879362; End=32879662; length=300; Strand='-', Score=2479.0)

Get CpG island elements in a genomic interval
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We query the genome for CpG islands within a specific coordinate range on chromosome 11.

.. doctest::

    >>> from cogent.db.ensembl import Genome
    >>> human = Genome('human', Release=67, account=account)
    >>> islands = human.getFeatures(CoordName='11', Start=2150341, End=2170833, feature_types='cpg')
    >>> for island in islands:
    ...     print island
    ...     break
    CpGisland(CoordName='11'; Start=2158951; End=2162484; length=3533; Strand='-', Score=3254.0)

Get SNPs
--------

For a gene
^^^^^^^^^^

We find the genetic variants for the canonical transcript of *BRCA2*.

.. note:: The output is significantly truncated!

.. doctest::

    >>> from cogent.db.ensembl import Genome
    >>> human = Genome('human', Release=67, account=account)
    >>> brca2 = human.getGeneByStableId(StableId='ENSG00000139618')
    >>> transcript = brca2.CanonicalTranscript
    >>> print transcript.Variants
    (<cogent.db.ensembl.region.Variation object at ...
    >>> for variant in transcript.Variants:
    ...     print variant
    ...     break
    Variation(Symbol='rs55880202'; Effect=['2KB_upstream_variant', '5_prime_UTR_variant', '5KB_upstream_variant']; Alleles='C/T')...

Get a single SNP
^^^^^^^^^^^^^^^^

We get a single SNP and print it's allele frequencies.

.. doctest::
    
    >>> snp = list(human.getVariation(Symbol='rs34213141'))[0]
    >>> print snp.AlleleFreqs
    =============================
    allele      freq    sample_id
    -----------------------------
         A    0.0303          933
         G    0.9697          933
         G    1.0000        11208
         G    1.0000        11519
         A    0.0110       113559
         G    0.9889       113559...

What alignment types available
------------------------------

We create a ``Compara`` instance for human, chimpanzee and macaque.

.. doctest::

    >>> from cogent.db.ensembl import Compara
    >>> compara = Compara(['human', 'chimp', 'macaque'], Release=67,
    ...                  account=account)
    >>> print compara.method_species_links
    Align Methods/Clades
    ===================================================================================================================
    method_link_species_set_id  method_link_id  species_set_id      align_method                            align_clade
    -------------------------------------------------------------------------------------------------------------------
                           580              10           34468             PECAN           19 amniota vertebrates Pecan
                           548              13           34115               EPO                         6 primates EPO
                           578              13           34466               EPO               12 eutherian mammals EPO
                           582              14           34697  EPO_LOW_COVERAGE  35 eutherian mammals EPO_LOW_COVERAGE
    -------------------------------------------------------------------------------------------------------------------

Get genomic alignment for a gene region
---------------------------------------

We first get the syntenic region corresponding to human gene *BRCA2*.

.. doctest::

    >>> from cogent.db.ensembl import Compara
    >>> compara = Compara(['human', 'chimp', 'macaque'], Release=67,
    ...                  account=account)
    >>> human_brca2 = compara.Human.getGeneByStableId(StableId='ENSG00000139618')
    >>> regions = compara.getSyntenicRegions(region=human_brca2, align_method='EPO', align_clade='primates')
    >>> for region in regions:
    ...     print region
    SyntenicRegions:
      Coordinate(Human,chro...,13,32889610-32973805,1)
      Coordinate(Macaque,chro...,17,11686607-11778803,1)
      Coordinate(Chimp,chro...,13,31957346-32040817,1)...

We then get a cogent ``Alignment`` object, requesting that sequences be annotated for gene spans.

.. doctest::

    >>> aln = region.getAlignment(feature_types='gene')
    >>> print repr(aln)
    3 x 98805 dna alignment: Homo sapiens:chromosome:13:3288...

Parsing syntenic regions
------------------------

Not all regions in a given genome have a syntenic alignment, and some have more than one alignment.
To illustrate these cases, we can consider an alignment between mouse and human, using the ``PECAN`` 
alignment method in the vertebrates clade:

.. doctest::

    >>> species = ["mouse", "human"]
    >>> compara = Compara(species, Release=67, account=account)
    >>> clade = "vertebrates"
    >>> chrom, start, end, strand = "X", 165754928, 165755079, "-"
    >>> regions = compara.getSyntenicRegions(Species="mouse", CoordName=chrom, 
    ...                                      Start=start, End=end, align_method="PECAN", 
    ...                                      align_clade=clade, Strand=strand)     
    >>> aligned_pairs = [r for r in regions]
    >>> alignment = aligned_pairs[0]                                                            
    >>> aligned_regions = [m for m in alignment.Members
    ...                    if m.Region is not None]
    >>> source_region, target_region = aligned_regions
    >>> print source_region.Location.CoordName, source_region.Location.Start, source_region.Location.End
    X 165754928 165755079
    >>> print target_region.Location.CoordName, target_region.Location.Start, target_region.Location.End
    X 11132954 11133105

.. note:: We took the aligned regions from the ``regions`` generator and put them in a list for convenience.

If there are no regions returned (i.e. ``num_pairs`` is zero), then no alignment could be found. In the case of 
the above region, an exon in the *Hccs* gene, there is only one alignment. We then accessed the coordinates of the 
alignment using the ``Members`` attribute of the region. Each element of ``aligned_regions`` is a ``SyntenicRegion``
instance, whose coordinates can be pulled from the ``Location`` attribute.

This example shows that mouse region ``X:165754928-165755079`` aligns only to human region ``X:11132954-11133105``.

.. note:: Sometimes, the genomic coordinates given to ``getSyntenicRegions`` will contain multiple alignments between the pair of genomes, in which case two or more regions will be returned in ``aligned_pairs``.

Getting related genes
---------------------

What gene relationships are available
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.db.ensembl import Compara
    >>> compara = Compara(['human', 'chimp', 'macaque'], Release=67,
    ...                  account=account)
    >>> print compara.getDistinct('relationship')
    [u'ortholog_one2many', u'contiguous_gene_split', u'ortholog_one2one',...

Get one-to-one orthologs
^^^^^^^^^^^^^^^^^^^^^^^^

We get the one-to-one orthologs for *BRCA2*.

.. doctest::

    >>> from cogent.db.ensembl import Compara
    >>> compara = Compara(['human', 'chimp', 'macaque'], Release=67,
    ...                  account=account)
    >>> orthologs = compara.getRelatedGenes(StableId='ENSG00000139618',
    ...                  Relationship='ortholog_one2one')
    >>> print orthologs
    RelatedGenes:
     Relationships=ortholog_one2one
      Gene(Species='Macaca mulatta'; BioType='protein_coding'; Description=...

We iterate over the related members.

.. doctest::
    
    >>> for ortholog in orthologs.Members:
    ...     print ortholog
    Gene(Species='Macaca mulatta'; BioType='protein_coding'; Description=...

We get statistics on the ortholog CDS lengths.

.. doctest::
    
    >>> print orthologs.getMaxCdsLengths()
    [10008, 10257, 10257]

We get the sequences as a sequence collection, with annotations for gene.

.. doctest::
    
    >>> seqs = orthologs.getSeqCollection(feature_types='gene')

Get CDS for all one-to-one orthologs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We sample all one-to-one orthologs for a group of species, generating a FASTA formatted string that can be written to file. We check all species have an ortholog and that all are translatable.

.. doctest::
    
    >>> from cogent.core.alphabet import AlphabetError
    >>> common_names = ["mouse", "rat", "human", "opossum"]
    >>> latin_names = set([Species.getSpeciesName(n) for n in common_names])
    >>> latin_to_common = dict(zip(latin_names, common_names))
    >>> compara = Compara(common_names, Release=67, account=account)
    >>> for gene in compara.Human.getGenesMatching(BioType='protein_coding'):
    ...     orthologs = compara.getRelatedGenes(gene,
    ...                                  Relationship='ortholog_one2one')
    ...     # make sure all species represented
    ...     if orthologs is None or orthologs.getSpeciesSet() != latin_names:
    ...         continue
    ...     seqs = []
    ...     for m in orthologs.Members:
    ...         try: # if sequence can't be translated, we ignore it
    ...             # get the CDS without the ending stop
    ...             seq = m.CanonicalTranscript.Cds.withoutTerminalStopCodon()
    ...             # make the sequence name
    ...             seq.Name = '%s:%s:%s' % \
    ...         (latin_to_common[m.genome.Species], m.StableId, m.Location)
    ...             aa = seq.getTranslation()
    ...             seqs += [seq]
    ...         except (AlphabetError, AssertionError):
    ...             seqs = [] # exclude this gene
    ...             break
    ...     if len(seqs) == len(common_names):
    ...         fasta = '\n'.join(s.toFasta() for s in seqs)
    ...         break

Get within species paralogs
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> paralogs = compara.getRelatedGenes(StableId='ENSG00000164032',
    ...             Relationship='within_species_paralog')
    >>> print paralogs
    RelatedGenes:
     Relationships=within_species_paralog
      Gene(Species='Homo sapiens'; BioType='protein_coding'; Description='H2A histone...

