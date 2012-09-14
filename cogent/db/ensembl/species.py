from cogent.util.table import Table
from cogent.db.ensembl.util import CaseInsensitiveString

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Jason Merkin"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

_species_common_map = [['Acyrthosiphon pisum', 'A.pisum'],
                       ['Aedes aegypti', 'A.aegypti'],
                       ['Anolis carolinensis', 'Anole lizard'],
                       ['Anopheles gambiae', 'A.gambiae'],
                       ['Apis mellifera', 'Honeybee'],
                       ['Arabidopsis lyrata', 'A.lyrata'],
                       ['Arabidopsis thaliana', 'A.thaliana'],
                       ['Aspergillus clavatus', 'A.clavatus'],
                       ['Aspergillus flavus', 'A.flavus'],
                       ['Aspergillus fumigatus', 'A.fumigatus'],
                       ['Aspergillus fumigatusa1163', 'A.fumigatusa1163'],
                       ['Aspergillus nidulans', 'A.nidulans'],
                       ['Aspergillus niger', 'A.niger'],
                       ['Aspergillus oryzae', 'A.oryzae'],
                       ['Aspergillus terreus', 'A.terreus'],
                       ['Bacillus collection', 'Bacillus'],
                       ['Borrelia collection', 'Borrelia'],
                       ['Bos taurus', 'Cow'],
                       ['Brachypodium distachyon', 'B.distachyon'],
                       ['Buchnera collection', 'Buchnera'],
                       ['Caenorhabditis brenneri', 'C.brenneri'],
                       ['Caenorhabditis briggsae', 'C.briggsae'],
                       ['Caenorhabditis elegans', 'C.elegans'],
                       ['Caenorhabditis japonica', 'C.japonica'],
                       ['Caenorhabditis remanei', 'C.remanei'],
                       ['Callithrix jacchus', 'Marmoset'],
                       ['Canis familiaris', 'Dog'],
                       ['Cavia porcellus', 'Guinea Pig'],
                       ['Choloepus hoffmanni', 'Sloth'],
                       ['Ciona intestinalis', 'C.intestinalis'],
                       ['Ciona savignyi', 'C.savignyi'],
                       ['Culex quinquefasciatus', 'C.quinquefasciatus'],
                       ['Danio rerio', 'Zebrafish'],
                       ['Dasypus novemcinctus', 'Armadillo'],
                       ['Dictyostelium discoideum', 'D.discoideum'],
                       ['Dipodomys ordii', 'Kangaroo rat'],
                       ['Drosophila ananassae', 'D.ananassae'],
                       ['Drosophila erecta', 'D.erecta'],
                       ['Drosophila grimshawi', 'D.grimshawi'],
                       ['Drosophila melanogaster', 'D.melanogaster'],
                       ['Drosophila mojavensis', 'D.mojavensis'],
                       ['Drosophila persimilis', 'D.persimilis'],
                       ['Drosophila pseudoobscura', 'D.pseudoobscura'],
                       ['Drosophila sechellia', 'D.sechellia'],
                       ['Drosophila simulans', 'D.simulans'],
                       ['Drosophila virilis', 'D.virilis'],
                       ['Drosophila willistoni', 'D.willistoni'],
                       ['Drosophila yakuba', 'D.yakuba'],
                       ['Echinops telfairi', 'Tenrec'],
                       ['Equus caballus', 'Horse'],
                       ['Erinaceus europaeus', 'Hedgehog'],
                       ['Escherichia shigella', 'E.shigella'],
                       ['Felis catus', 'Cat'],
                       ['Gallus gallus', 'Chicken'],
                       ['Gasterosteus aculeatus', 'Stickleback'],
                       ['Gorilla gorilla', 'Gorilla'],
                       ['Homo sapiens', 'Human'],
                       ['Ixodes scapularis', 'I.scapularis'],
                       ['Loxodonta africana', 'Elephant'],
                       ['Macaca mulatta', 'Macaque'],
                       ['Macropus eugenii', 'Wallaby'],
                       ['Microcebus murinus', 'Mouse lemur'],
                       ['Monodelphis domestica', 'Opossum'],
                       ['Mus musculus', 'Mouse'],
                       ['Mycobacterium collection', 'Mycobacterium'],
                       ['Myotis lucifugus', 'Microbat'],
                       ['Neisseria collection', 'Neisseria'],
                       ['Nematostella vectensis', 'N.vectensis'],
                       ['Neosartorya fischeri', 'N.fischeri'],
                       ['Neurospora crassa', 'N.crassa'],
                       ['Ochotona princeps', 'Pika'],
                       ['Ornithorhynchus anatinus', 'Platypus'],
                       ['Oryctolagus cuniculus', 'Rabbit'],
                       ['Oryza indica', 'O.indica'],
                       ['Oryza sativa', 'O.sativa'],
                       ['Oryzias latipes', 'Medaka'],
                       ['Otolemur garnettii', 'Bushbaby'],
                       ['Pan troglodytes', 'Chimp'],
                       ['Pediculus humanus', 'P.humanus'],
                       ['Petromyzon marinus', 'Lamprey'],
                       ['Phaeodactylum tricornutum', 'P.tricornutum'],
                       ['Physcomitrella patens', 'P.patens'],
                       ['Phytophthora infestans', 'P.infestans'],
                       ['Phytophthora ramorum', 'P.ramorum'],
                       ['Phytophthora sojae', 'P.sojae'],
                       ['Plasmodium berghei', 'P.berghei'],
                       ['Plasmodium chabaudi', 'P.chabaudi'],
                       ['Plasmodium falciparum', 'P.falciparum'],
                       ['Plasmodium knowlesi', 'P.knowlesi'],
                       ['Plasmodium vivax', 'P.vivax'],
                       ['Pongo pygmaeus', 'Orangutan'],
                       ['Populus trichocarpa', 'P.trichocarpa'],
                       ['Pristionchus pacificus', 'P.pacificus'],
                       ['Procavia capensis', 'Rock hyrax'],
                       ['Pteropus vampyrus', 'Flying fox'],
                       ['Puccinia graministritici', 'P.graministritici'],
                       ['Pyrococcus collection', 'Pyrococcus'],
                       ['Rattus norvegicus', 'Rat'],
                       ['Saccharomyces cerevisiae', 'S.cerevisiae'],
                       ['Schistosoma mansoni', 'S.mansoni'],
                       ['Schizosaccharomyces pombe', 'S.pombe'],
                       ['Sorex araneus', 'Shrew'],
                       ['Sorghum bicolor', 'S.bicolor'],
                       ['Spermophilus tridecemlineatus', 'Ground Squirrel'],
                       ['Staphylococcus collection', 'Staphylococcus'],
                       ['Streptococcus collection', 'Streptococcus'],
                       ['Strongylocentrotus purpuratus', 'S.purpuratus'],
                       ['Sus scrofa', 'Pig'],
                       ['Taeniopygia guttata', 'Zebra finch'],
                       ['Takifugu rubripes', 'Fugu'],
                       ['Tarsius syrichta', 'Tarsier'],
                       ['Tetraodon nigroviridis', 'Tetraodon'],
                       ['Thalassiosira pseudonana', 'T.pseudonana'],
                       ['Trichoplax adhaerens', 'T.adhaerens'],
                       ['Tupaia belangeri', 'Tree Shrew'],
                       ['Tursiops truncatus', 'Bottlenose dolphin'],
                       ['Vicugna pacos', 'Alpaca'],
                       ['Vitis vinifera', 'V.vinifera'],
                       ['Wolbachia collection', 'Wolbachia'],
                       ['Xenopus tropicalis', 'Xenopus'],
                       ['Zea mays', 'Z.mays']]

class SpeciesNameMap(dict):
    """mapping between common names and latin names"""
    def __init__(self, species_common = _species_common_map):
        """provides latin name:common name mappings"""
        self._species_common = {}
        self._common_species = {}
        self._species_ensembl = {}
        self._ensembl_species = {}
        for species_name, common_name in species_common:
            self.amendSpecies(CaseInsensitiveString(species_name),
                              CaseInsensitiveString(common_name))
    
    def __str__(self):
        rows = []
        for common in self._common_species:
            species = self._common_species[common]
            ensembl = self._species_ensembl[species]
            rows += [[common, species, ensembl]]
        return str(Table(['Common Name', 'Species Name', 'Ensembl Db Prefix'],
                    rows=rows, space=2).sorted())
    
    def __repr__(self):
        return 'Available species: %s' % ("'"+\
                "'; '".join(self._common_species.keys())+"'")
    
    def getCommonName(self, name):
        """returns the common name for the given name (which can be either a
        species name or the ensembl version)"""
        name = CaseInsensitiveString(name)
        if name in self._ensembl_species:
            name = self._ensembl_species[name]
        
        if name in self._species_common:
            common_name = self._species_common[name]
        elif name in self._common_species:
            common_name = name
        else:
            raise RuntimeError("Unknown species: %s" % name)
        return str(common_name)
    
    def getSpeciesName(self, name, level='ignore'):
        """returns the species name for the given common name"""
        name = CaseInsensitiveString(name)
        if name in self._species_common:
            return str(name)
        species_name = None
        level = level.lower().strip()
        name = name
        for data in [self._common_species, self._ensembl_species]:
            if name in data:
                species_name = data[name]
        if species_name is None:
            msg = "Unknown common name: %s" % name
            if level == 'raise':
                raise RuntimeError(msg)
            elif level == 'warn':
                print "WARN: %s" % msg
        return str(species_name)
    
    def getSpeciesNames(self):
        """returns the list of species names"""
        names = self._species_common.keys()
        names.sort()
        return [str(n) for n in names]
    
    def getEnsemblDbPrefix(self, name):
        """returns a string of the species name in the format used by
        ensembl"""
        name = CaseInsensitiveString(name)
        if name in self._common_species:
            name = self._common_species[name]
        try:
            species_name = self.getSpeciesName(name, level='raise')
        except RuntimeError:
            if name not in self._species_common:
                raise RuntimeError("Unknown name %s" % name)
            species_name = name
        
        return str(species_name.lower().replace(" ","_"))
    
    def getComparaName(self, name):
        """returns string matching a compara instance attribute name for a
        species"""
        name = self.getCommonName(name)
        if '.' in name:
            name = name.replace('.', '')
        else:
            name = name.title()
        
        name = name.split()
        return ''.join(name)
    
    def _purge_species(self, species_name):
        """removes a species record"""
        species_name = CaseInsensitiveString(species_name)
        if not species_name in self._species_common:
            return
        common_name = self._species_common.pop(species_name)
        ensembl_name= self._species_ensembl.pop(species_name)
        self._ensembl_species.pop(ensembl_name)
        self._common_species.pop(common_name)
    
    def amendSpecies(self, species_name, common_name):
        """add a new species, and common name"""
        species_name = CaseInsensitiveString(species_name)
        common_name = CaseInsensitiveString(common_name)
        assert "_" not in species_name,\
                        "'_' in species_name, not a Latin name?"
        self._purge_species(species_name) # remove if existing
        self._species_common[species_name] = common_name
        self._common_species[common_name] = species_name
        ensembl_name = species_name.lower().replace(" ","_")
        self._species_ensembl[species_name] = ensembl_name
        self._ensembl_species[ensembl_name] = species_name
        return

Species = SpeciesNameMap()
