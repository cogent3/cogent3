#!/usr/bin/env python
"""Unit tests for the pdb parser.
"""
from cogent.util.unit_test import TestCase, main
from cogent.parse.pdb import dict2pdb, dict2ter, pdb2dict, get_symmetry, \
                             get_coords_offset, get_trailer_offset, \
                             parse_header, parse_coords, parse_trailer, \
                             PDBParser
from cogent.core.entity import Structure
from cogent.core.entity import StructureBuilder
from numpy import array, allclose

__author__ = "Marcin Knight"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Production"


class pdbTests(TestCase):
    """Tests of cogent.parse.pdb functions."""

    def test_PDBParser(self):
        """tests the UI parsing function.
        """
        fh = open('data/2E12.pdb')
        structure = PDBParser(fh, 'JUNK')
        assert type(structure) is Structure
        assert len(structure) == 1
        assert (0,) in structure
        assert structure.header['space_group'] == 'P 21 21 21'
        assert structure.header['experiment_type'] == 'X-RAY DIFFRACTION'
        assert structure.header['r_free'] == '0.280'
        assert structure.header['dbref_acc'] == 'Q8P4R5'
        assert structure.header['cryst1'] == '49.942   51.699   82.120  90.00  90.00  90.00'
        assert structure.header['matthews'] == '2.29'
        model = structure[(0,)]
        assert len(model) == 2
        assert structure.raw_header
        assert structure.raw_trailer
        assert structure.header
        assert structure.trailer ==  {}
        assert structure.getId() == ('JUNK', )
        
    def test_parse_trailer(self):
        """testing trailer parsing dummy."""
        d = parse_trailer(None)
        assert isinstance(d, dict)
        
    def test_parse_coords(self):
        """testing minimal structure building and coords parsing.
        """
        builder = StructureBuilder()
        builder.initStructure('JUNK')
        atom = 'ATOM     10  CA  PRO A   2      51.588  38.262  31.417  1.00  6.58           C  \n'
        hetatm = 'HETATM 1633  O   HOH B 164      17.979  35.529  38.171  1.00  1.02           O  \n'
        lines = ['MODEL ', atom, hetatm]
        z = parse_coords(builder, lines)
        assert len(z[(0,)]) == 2
        assert len(z[(0,)][('A',)]) == 1
        assert len(z[(0,)][('B',)]) == 1
        z.setTable()
        
        atom1 = z.table['A'][('JUNK', 0, 'A', ('PRO', 2, ' '), ('CA', ' '))]
        hetatm1 = z.table['A'][('JUNK', 0, 'B', ('H_HOH', 164, ' '), ('O', ' '))]
        
        self.assertAlmostEqual([51.588 , 38.262 , 31.417][2], list(atom1.coords)[2])
        self.assertAlmostEqual([17.979 ,  35.529 , 38.171][2], list(hetatm1.coords)[2])

    def test_parse_header(self):
        """testing header parsing.
        """
        header = ['HEADER    TRANSLATION                             17-OCT-06   2E12              \n',
 'TITLE     THE CRYSTAL STRUCTURE OF XC5848 FROM XANTHOMONAS CAMPESTRIS           \n',
 'TITLE    2 ADOPTING A NOVEL VARIANT OF SM-LIKE MOTIF                            \n',
 'COMPND    MOL_ID: 1;                                                            \n',
 'COMPND   2 MOLECULE: HYPOTHETICAL PROTEIN XCC3642;                              \n',
 'COMPND   3 CHAIN: A, B;                                                         \n',
 'COMPND   4 SYNONYM: SM-LIKE MOTIF;                                              \n',
 'COMPND   5 ENGINEERED: YES                                                      \n',
 'SOURCE    MOL_ID: 1;                                                            \n',
 'SOURCE   2 ORGANISM_SCIENTIFIC: XANTHOMONAS CAMPESTRIS PV. CAMPESTRIS;          \n',
 'SOURCE   3 ORGANISM_TAXID: 340;                                                 \n',
 'SOURCE   4 STRAIN: PV. CAMPESTRIS;                                              \n',
 'SOURCE   5 EXPRESSION_SYSTEM: ESCHERICHIA COLI;                                 \n',
 'SOURCE   6 EXPRESSION_SYSTEM_TAXID: 562                                         \n',
 'KEYWDS    NOVEL SM-LIKE MOTIF, LSM MOTIF, XANTHOMONAS CAMPESTRIS, X-            \n',
 'KEYWDS   2 RAY CRYSTALLOGRAPHY, TRANSLATION                                     \n',
 'EXPDTA    X-RAY DIFFRACTION                                                     \n',
 'AUTHOR    K.-H.CHIN,S.-K.RUAN,A.H.-J.WANG,S.-H.CHOU                             \n',
 'REVDAT   2   24-FEB-09 2E12    1       VERSN                                    \n',
 'REVDAT   1   30-OCT-07 2E12    0                                                \n',
 'JRNL        AUTH   K.-H.CHIN,S.-K.RUAN,A.H.-J.WANG,S.-H.CHOU                    \n',
 'JRNL        TITL   XC5848, AN ORFAN PROTEIN FROM XANTHOMONAS                    \n',
 'JRNL        TITL 2 CAMPESTRIS, ADOPTS A NOVEL VARIANT OF SM-LIKE MOTIF          \n',
 'JRNL        REF    PROTEINS                      V.  68  1006 2007              \n',
 'JRNL        REFN                   ISSN 0887-3585                               \n',
 'JRNL        PMID   17546661                                                     \n',
 'JRNL        DOI    10.1002/PROT.21375                                           \n',
 'REMARK   1                                                                      \n',
 'REMARK   2                                                                      \n',
 'REMARK   2 RESOLUTION.    1.70 ANGSTROMS.                                       \n',
 'REMARK   3                                                                      \n',
 'REMARK   3 REFINEMENT.                                                          \n',
 'REMARK   3   PROGRAM     : CNS                                                  \n',
 'REMARK   3   AUTHORS     : BRUNGER,ADAMS,CLORE,DELANO,GROS,GROSSE-              \n',
 'REMARK   3               : KUNSTLEVE,JIANG,KUSZEWSKI,NILGES, PANNU,             \n',
 'REMARK   3               : READ,RICE,SIMONSON,WARREN                            \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  REFINEMENT TARGET : NULL                                            \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  DATA USED IN REFINEMENT.                                            \n',
 'REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 1.70                           \n',
 'REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 30.00                          \n',
 'REMARK   3   DATA CUTOFF            (SIGMA(F)) : 5.000                          \n',
 'REMARK   3   DATA CUTOFF HIGH         (ABS(F)) : NULL                           \n',
 'REMARK   3   DATA CUTOFF LOW          (ABS(F)) : NULL                           \n',
 'REMARK   3   COMPLETENESS (WORKING+TEST)   (%) : 99.1                           \n',
 'REMARK   3   NUMBER OF REFLECTIONS             : 6937                           \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  FIT TO DATA USED IN REFINEMENT.                                     \n',
 'REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT                      \n',
 'REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM                          \n',
 'REMARK   3   R VALUE            (WORKING SET) : 0.220                           \n',
 'REMARK   3   FREE R VALUE                     : 0.280                           \n',
 'REMARK   3   FREE R VALUE TEST SET SIZE   (%) : NULL                            \n',
 'REMARK   3   FREE R VALUE TEST SET COUNT      : NULL                            \n',
 'REMARK   3   ESTIMATED ERROR OF FREE R VALUE  : NULL                            \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  FIT IN THE HIGHEST RESOLUTION BIN.                                  \n',
 'REMARK   3   TOTAL NUMBER OF BINS USED           : NULL                         \n',
 'REMARK   3   BIN RESOLUTION RANGE HIGH       (A) : 1.70                         \n',
 'REMARK   3   BIN RESOLUTION RANGE LOW        (A) : 1.75                         \n',
 'REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%) : 97.00                        \n',
 'REMARK   3   REFLECTIONS IN BIN    (WORKING SET) : NULL                         \n',
 'REMARK   3   BIN R VALUE           (WORKING SET) : 0.2400                       \n',
 'REMARK   3   BIN FREE R VALUE                    : 0.2200                       \n',
 'REMARK   3   BIN FREE R VALUE TEST SET SIZE  (%) : NULL                         \n',
 'REMARK   3   BIN FREE R VALUE TEST SET COUNT     : NULL                         \n',
 'REMARK   3   ESTIMATED ERROR OF BIN FREE R VALUE : 0.012                        \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.                    \n',
 'REMARK   3   PROTEIN ATOMS            : 1512                                    \n',
 'REMARK   3   NUCLEIC ACID ATOMS       : 0                                       \n',
 'REMARK   3   HETEROGEN ATOMS          : 0                                       \n',
 'REMARK   3   SOLVENT ATOMS            : 122                                     \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  B VALUES.                                                           \n',
 'REMARK   3   FROM WILSON PLOT           (A**2) : 24.00                          \n',
 'REMARK   3   MEAN B VALUE      (OVERALL, A**2) : NULL                           \n',
 'REMARK   3   OVERALL ANISOTROPIC B VALUE.                                       \n',
 'REMARK   3    B11 (A**2) : NULL                                                 \n',
 'REMARK   3    B22 (A**2) : NULL                                                 \n',
 'REMARK   3    B33 (A**2) : NULL                                                 \n',
 'REMARK   3    B12 (A**2) : NULL                                                 \n',
 'REMARK   3    B13 (A**2) : NULL                                                 \n',
 'REMARK   3    B23 (A**2) : NULL                                                 \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  ESTIMATED COORDINATE ERROR.                                         \n',
 'REMARK   3   ESD FROM LUZZATI PLOT        (A) : NULL                            \n',
 'REMARK   3   ESD FROM SIGMAA              (A) : NULL                            \n',
 'REMARK   3   LOW RESOLUTION CUTOFF        (A) : NULL                            \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  CROSS-VALIDATED ESTIMATED COORDINATE ERROR.                         \n',
 'REMARK   3   ESD FROM C-V LUZZATI PLOT    (A) : NULL                            \n',
 'REMARK   3   ESD FROM C-V SIGMAA          (A) : NULL                            \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  RMS DEVIATIONS FROM IDEAL VALUES.                                   \n',
 'REMARK   3   BOND LENGTHS                 (A) : 0.007                           \n',
 'REMARK   3   BOND ANGLES            (DEGREES) : 1.32                            \n',
 'REMARK   3   DIHEDRAL ANGLES        (DEGREES) : NULL                            \n',
 'REMARK   3   IMPROPER ANGLES        (DEGREES) : NULL                            \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  ISOTROPIC THERMAL MODEL : NULL                                      \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  ISOTROPIC THERMAL FACTOR RESTRAINTS.    RMS    SIGMA                \n',
 'REMARK   3   MAIN-CHAIN BOND              (A**2) : NULL  ; NULL                 \n',
 'REMARK   3   MAIN-CHAIN ANGLE             (A**2) : NULL  ; NULL                 \n',
 'REMARK   3   SIDE-CHAIN BOND              (A**2) : NULL  ; NULL                 \n',
 'REMARK   3   SIDE-CHAIN ANGLE             (A**2) : NULL  ; NULL                 \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  BULK SOLVENT MODELING.                                              \n',
 'REMARK   3   METHOD USED : NULL                                                 \n',
 'REMARK   3   KSOL        : NULL                                                 \n',
 'REMARK   3   BSOL        : NULL                                                 \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  NCS MODEL : NULL                                                    \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  NCS RESTRAINTS.                         RMS   SIGMA/WEIGHT          \n',
 'REMARK   3   GROUP  1  POSITIONAL            (A) : NULL  ; NULL                 \n',
 'REMARK   3   GROUP  1  B-FACTOR           (A**2) : NULL  ; NULL                 \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  PARAMETER FILE  1  : NULL                                           \n',
 'REMARK   3  TOPOLOGY FILE  1   : NULL                                           \n',
 'REMARK   3                                                                      \n',
 'REMARK   3  OTHER REFINEMENT REMARKS: NULL                                      \n',
 'REMARK   4                                                                      \n',
 'REMARK   4 2E12 COMPLIES WITH FORMAT V. 3.15, 01-DEC-08                         \n',
 'REMARK 100                                                                      \n',
 'REMARK 100 THIS ENTRY HAS BEEN PROCESSED BY PDBJ ON 19-OCT-06.                  \n',
 'REMARK 100 THE RCSB ID CODE IS RCSB026092.                                      \n',
 'REMARK 200                                                                      \n',
 'REMARK 200 EXPERIMENTAL DETAILS                                                 \n',
 'REMARK 200  EXPERIMENT TYPE                : X-RAY DIFFRACTION                  \n',
 'REMARK 200  DATE OF DATA COLLECTION        : 28-JUL-06                          \n',
 'REMARK 200  TEMPERATURE           (KELVIN) : 100                                \n',
 'REMARK 200  PH                             : 8.0                                \n',
 'REMARK 200  NUMBER OF CRYSTALS USED        : 10                                 \n',
 'REMARK 200                                                                      \n',
 'REMARK 200  SYNCHROTRON              (Y/N) : Y                                  \n',
 'REMARK 200  RADIATION SOURCE               : NSRRC                              \n',
 'REMARK 200  BEAMLINE                       : BL13B1                             \n',
 'REMARK 200  X-RAY GENERATOR MODEL          : NULL                               \n',
 'REMARK 200  MONOCHROMATIC OR LAUE    (M/L) : M                                  \n',
 'REMARK 200  WAVELENGTH OR RANGE        (A) : 0.96437, 0.97983                   \n',
 'REMARK 200  MONOCHROMATOR                  : NULL                               \n',
 'REMARK 200  OPTICS                         : NULL                               \n',
 'REMARK 200                                                                      \n',
 'REMARK 200  DETECTOR TYPE                  : CCD                                \n',
 'REMARK 200  DETECTOR MANUFACTURER          : ADSC QUANTUM 315                   \n',
 'REMARK 200  INTENSITY-INTEGRATION SOFTWARE : DENZO                              \n',
 'REMARK 200  DATA SCALING SOFTWARE          : HKL-2000                           \n',
 'REMARK 200                                                                      \n',
 'REMARK 200  NUMBER OF UNIQUE REFLECTIONS   : 6937                               \n',
 'REMARK 200  RESOLUTION RANGE HIGH      (A) : 1.700                              \n',
 'REMARK 200  RESOLUTION RANGE LOW       (A) : 30.000                             \n',
 'REMARK 200  REJECTION CRITERIA  (SIGMA(I)) : 2.000                              \n',
 'REMARK 200                                                                      \n',
 'REMARK 200 OVERALL.                                                             \n',
 'REMARK 200  COMPLETENESS FOR RANGE     (%) : 99.7                               \n',
 'REMARK 200  DATA REDUNDANCY                : 4.500                              \n',
 'REMARK 200  R MERGE                    (I) : 0.24000                            \n',
 'REMARK 200  R SYM                      (I) : 0.06000                            \n',
 'REMARK 200  <I/SIGMA(I)> FOR THE DATA SET  : 8.0000                             \n',
 'REMARK 200                                                                      \n',
 'REMARK 200 IN THE HIGHEST RESOLUTION SHELL.                                     \n',
 'REMARK 200  HIGHEST RESOLUTION SHELL, RANGE HIGH (A) : 1.70                     \n',
 'REMARK 200  HIGHEST RESOLUTION SHELL, RANGE LOW  (A) : NULL                     \n',
 'REMARK 200  COMPLETENESS FOR SHELL     (%) : 97.5                               \n',
 'REMARK 200  DATA REDUNDANCY IN SHELL       : 4.50                               \n',
 'REMARK 200  R MERGE FOR SHELL          (I) : 0.06000                            \n',
 'REMARK 200  R SYM FOR SHELL            (I) : 0.24000                            \n',
 'REMARK 200  <I/SIGMA(I)> FOR SHELL         : 7.900                              \n',
 'REMARK 200                                                                      \n',
 'REMARK 200 DIFFRACTION PROTOCOL: MAD                                            \n',
 'REMARK 200 METHOD USED TO DETERMINE THE STRUCTURE: MAD                          \n',
 'REMARK 200 SOFTWARE USED: AMORE                                                 \n',
 'REMARK 200 STARTING MODEL: NULL                                                 \n',
 'REMARK 200                                                                      \n',
 'REMARK 200 REMARK: NULL                                                         \n',
 'REMARK 280                                                                      \n',
 'REMARK 280 CRYSTAL                                                              \n',
 'REMARK 280 SOLVENT CONTENT, VS   (%): 46.26                                     \n',
 'REMARK 280 MATTHEWS COEFFICIENT, VM (ANGSTROMS**3/DA): 2.29                     \n',
 'REMARK 280                                                                      \n',
 'REMARK 280 CRYSTALLIZATION CONDITIONS: PH 8.0, VAPOR DIFFUSION, SITTING         \n',
 'REMARK 280  DROP, TEMPERATURE 298K                                              \n',
 'REMARK 290 REMARK: NULL                                                         \n',
 'REMARK 300                                                                      \n',
 'REMARK 300 BIOMOLECULE: 1                                                       \n',
 'REMARK 300 SEE REMARK 350 FOR THE AUTHOR PROVIDED AND/OR PROGRAM                \n',
 'REMARK 300 GENERATED ASSEMBLY INFORMATION FOR THE STRUCTURE IN                  \n',
 'REMARK 300 THIS ENTRY. THE REMARK MAY ALSO PROVIDE INFORMATION ON               \n',
 'REMARK 300 BURIED SURFACE AREA.                                                 \n',
 'REMARK 465                                                                      \n',
 'REMARK 465 MISSING RESIDUES                                                     \n',
 'REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE                       \n',
 'REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN               \n',
 'REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)                \n',
 'REMARK 465                                                                      \n',
 'REMARK 465   M RES C SSSEQI                                                     \n',
 'REMARK 465     LEU A    94                                                      \n',
 'REMARK 465     GLY A    95                                                      \n',
 'REMARK 465     ALA A    96                                                      \n',
 'REMARK 465     PRO A    97                                                      \n',
 'REMARK 465     GLN A    98                                                      \n',
 'REMARK 465     VAL A    99                                                      \n',
 'REMARK 465     MET A   100                                                      \n',
 'REMARK 465     PRO A   101                                                      \n',
 'REMARK 465     LEU B    94                                                      \n',
 'REMARK 465     GLY B    95                                                      \n',
 'REMARK 465     ALA B    96                                                      \n',
 'REMARK 465     PRO B    97                                                      \n',
 'REMARK 465     GLN B    98                                                      \n',
 'REMARK 465     VAL B    99                                                      \n',
 'REMARK 465     MET B   100                                                      \n',
 'REMARK 465     PRO B   101                                                      \n',
 'REMARK 500                                                                      \n',
 'REMARK 500 GEOMETRY AND STEREOCHEMISTRY                                         \n',
 'REMARK 500 SUBTOPIC: CLOSE CONTACTS IN SAME ASYMMETRIC UNIT                     \n',
 'REMARK 500                                                                      \n',
 'REMARK 500 THE FOLLOWING ATOMS ARE IN CLOSE CONTACT.                            \n',
 'REMARK 500                                                                      \n',
 'REMARK 500  ATM1  RES C  SSEQI   ATM2  RES C  SSEQI           DISTANCE          \n',
 'REMARK 500   O    HOH A   127     O    HOH A   149              2.05            \n',
 'REMARK 500                                                                      \n',
 'REMARK 500 REMARK: NULL                                                         \n',
 'REMARK 500                                                                      \n',
 'REMARK 500 GEOMETRY AND STEREOCHEMISTRY                                         \n',
 'REMARK 500 SUBTOPIC: TORSION ANGLES                                             \n',
 'REMARK 500                                                                      \n',
 'REMARK 500 TORSION ANGLES OUTSIDE THE EXPECTED RAMACHANDRAN REGIONS:            \n',
 'REMARK 500 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN IDENTIFIER;               \n',
 'REMARK 500 SSEQ=SEQUENCE NUMBER; I=INSERTION CODE).                             \n',
 'REMARK 500                                                                      \n',
 'REMARK 500 STANDARD TABLE:                                                      \n',
 'REMARK 500 FORMAT:(10X,I3,1X,A3,1X,A1,I4,A1,4X,F7.2,3X,F7.2)                    \n',
 'REMARK 500                                                                      \n',
 'REMARK 500 EXPECTED VALUES: GJ KLEYWEGT AND TA JONES (1996). PHI/PSI-           \n',
 'REMARK 500 CHOLOGY: RAMACHANDRAN REVISITED. STRUCTURE 4, 1395 - 1400            \n',
 'REMARK 500                                                                      \n',
 'REMARK 500  M RES CSSEQI        PSI       PHI                                   \n',
 'REMARK 500    ASN A  64     -175.84   -178.56                                   \n',
 'REMARK 500    HIS A  71     -156.72   -164.33                                   \n',
 'REMARK 500    LEU A  72      -70.52   -135.73                                   \n',
 'REMARK 500    ALA A  74      -75.47    -29.45                                   \n',
 'REMARK 500    SER A  75       -5.54   -145.62                                   \n',
 'REMARK 500    GLN A  76     -178.36     65.32                                   \n',
 'REMARK 500    GLU A  77      115.33     61.52                                   \n',
 'REMARK 500    MET A  92      -36.89     93.21                                   \n',
 'REMARK 500    LEU B  25       37.31    -77.89                                   \n',
 'REMARK 500    GLN B  28       37.09     32.85                                   \n',
 'REMARK 500    ARG B  30      132.20    -36.84                                   \n',
 'REMARK 500    ASN B  64     -172.78   -175.93                                   \n',
 'REMARK 500    GLN B  76       67.64     34.46                                   \n',
 'REMARK 500    PRO B  91     -156.99    -48.61                                   \n',
 'REMARK 500    MET B  92      -37.52   -160.44                                   \n',
 'REMARK 500                                                                      \n',
 'REMARK 500 REMARK: NULL                                                         \n',
 'REMARK 525                                                                      \n',
 'REMARK 525 SOLVENT                                                              \n',
 'REMARK 525                                                                      \n',
 'REMARK 525 THE SOLVENT MOLECULES HAVE CHAIN IDENTIFIERS THAT                    \n',
 'REMARK 525 INDICATE THE POLYMER CHAIN WITH WHICH THEY ARE MOST                  \n',
 'REMARK 525 CLOSELY ASSOCIATED. THE REMARK LISTS ALL THE SOLVENT                 \n',
 'REMARK 525 MOLECULES WHICH ARE MORE THAN 5A AWAY FROM THE                       \n',
 'REMARK 525 NEAREST POLYMER CHAIN (M = MODEL NUMBER;                             \n',
 'REMARK 525 RES=RESIDUE NAME; C=CHAIN IDENTIFIER; SSEQ=SEQUENCE                  \n',
 'REMARK 525 NUMBER; I=INSERTION CODE):                                           \n',
 'REMARK 525                                                                      \n',
 'REMARK 525  M RES CSSEQI                                                        \n',
 'REMARK 525    HOH B 115        DISTANCE =  6.82 ANGSTROMS                       \n',
 'REMARK 525    HOH A 116        DISTANCE =  6.52 ANGSTROMS                       \n',
 'REMARK 525    HOH B 119        DISTANCE =  5.12 ANGSTROMS                       \n',
 'REMARK 525    HOH B 121        DISTANCE =  5.21 ANGSTROMS                       \n',
 'REMARK 525    HOH B 123        DISTANCE =  5.18 ANGSTROMS                       \n',
 'REMARK 525    HOH A 124        DISTANCE =  6.99 ANGSTROMS                       \n',
 'REMARK 525    HOH B 124        DISTANCE =  5.13 ANGSTROMS                       \n',
 'REMARK 525    HOH B 134        DISTANCE =  7.25 ANGSTROMS                       \n',
 'REMARK 525    HOH B 140        DISTANCE =  5.54 ANGSTROMS                       \n',
 'REMARK 525    HOH B 141        DISTANCE =  5.94 ANGSTROMS                       \n',
 'REMARK 525    HOH B 142        DISTANCE =  6.60 ANGSTROMS                       \n',
 'REMARK 525    HOH B 143        DISTANCE =  7.39 ANGSTROMS                       \n',
 'REMARK 525    HOH A 145        DISTANCE =  9.25 ANGSTROMS                       \n',
 'REMARK 525    HOH A 150        DISTANCE =  6.01 ANGSTROMS                       \n',
 'REMARK 525    HOH B 152        DISTANCE =  5.46 ANGSTROMS                       \n',
 'REMARK 525    HOH B 153        DISTANCE =  9.74 ANGSTROMS                       \n',
 'REMARK 525    HOH B 154        DISTANCE =  9.32 ANGSTROMS                       \n',
 'REMARK 525    HOH B 155        DISTANCE =  5.41 ANGSTROMS                       \n',
 'REMARK 525    HOH B 163        DISTANCE =  5.16 ANGSTROMS                       \n',
 'DBREF  2E12 A    1   101  UNP    Q8P4R5   Q8P4R5_XANCP     1    101             \n',
 'DBREF  2E12 B    1   101  UNP    Q8P4R5   Q8P4R5_XANCP     1    101             \n',
 'SEQRES   1 A  101  MET PRO LYS TYR ALA PRO HIS VAL TYR THR GLU GLN ALA          \n',
 'SEQRES   2 A  101  GLN ILE ALA THR LEU GLU HIS TRP VAL LYS LEU LEU ASP          \n',
 'SEQRES   3 A  101  GLY GLN GLU ARG VAL ARG ILE GLU LEU ASP ASP GLY SER          \n',
 'SEQRES   4 A  101  MET ILE ALA GLY THR VAL ALA VAL ARG PRO THR ILE GLN          \n',
 'SEQRES   5 A  101  THR TYR ARG ASP GLU GLN GLU ARG GLU GLY SER ASN GLY          \n',
 'SEQRES   6 A  101  GLN LEU ARG ILE ASP HIS LEU ASP ALA SER GLN GLU PRO          \n',
 'SEQRES   7 A  101  GLN TRP ILE TRP MET ASP ARG ILE VAL ALA VAL HIS PRO          \n',
 'SEQRES   8 A  101  MET PRO LEU GLY ALA PRO GLN VAL MET PRO                      \n',
 'SEQRES   1 B  101  MET PRO LYS TYR ALA PRO HIS VAL TYR THR GLU GLN ALA          \n',
 'SEQRES   2 B  101  GLN ILE ALA THR LEU GLU HIS TRP VAL LYS LEU LEU ASP          \n',
 'SEQRES   3 B  101  GLY GLN GLU ARG VAL ARG ILE GLU LEU ASP ASP GLY SER          \n',
 'SEQRES   4 B  101  MET ILE ALA GLY THR VAL ALA VAL ARG PRO THR ILE GLN          \n',
 'SEQRES   5 B  101  THR TYR ARG ASP GLU GLN GLU ARG GLU GLY SER ASN GLY          \n',
 'SEQRES   6 B  101  GLN LEU ARG ILE ASP HIS LEU ASP ALA SER GLN GLU PRO          \n',
 'SEQRES   7 B  101  GLN TRP ILE TRP MET ASP ARG ILE VAL ALA VAL HIS PRO          \n',
 'SEQRES   8 B  101  MET PRO LEU GLY ALA PRO GLN VAL MET PRO                      \n',
 'FORMUL   3  HOH   *122(H2 O)                                                    \n',
 'HELIX    1   1 GLU A   11  LEU A   24  1                                  14    \n',
 'HELIX    2   2 GLU B   11  LEU B   25  1                                  15    \n',
 'SHEET    1   A 3 ILE A  51  ARG A  55  0                                        \n',
 'SHEET    2   A 3 GLU A  61  ASP A  70 -1  O  ASN A  64   N  GLN A  52           \n',
 'SHEET    3   A 3 GLN A  79  TRP A  82 -1  O  ILE A  81   N  LEU A  67           \n',
 'SHEET    1   B 5 ILE A  51  ARG A  55  0                                        \n',
 'SHEET    2   B 5 GLU A  61  ASP A  70 -1  O  ASN A  64   N  GLN A  52           \n',
 'SHEET    3   B 5 MET A  40  VAL A  45 -1  N  THR A  44   O  ASP A  70           \n',
 'SHEET    4   B 5 ARG A  30  LEU A  35 -1  N  ILE A  33   O  ILE A  41           \n',
 'SHEET    5   B 5 ILE A  86  PRO A  91 -1  O  VAL A  87   N  GLU A  34           \n',
 'SHEET    1   C 5 PRO B  78  TRP B  82  0                                        \n',
 'SHEET    2   C 5 GLN B  66  ASP B  70 -1  N  ILE B  69   O  GLN B  79           \n',
 'SHEET    3   C 5 MET B  40  VAL B  47 -1  N  ALA B  46   O  ARG B  68           \n',
 'SHEET    4   C 5 VAL B  31  LEU B  35 -1  N  ILE B  33   O  ILE B  41           \n',
 'SHEET    5   C 5 ILE B  86  HIS B  90 -1  O  VAL B  87   N  GLU B  34           \n',
 'SHEET    1   D 2 GLN B  52  ARG B  55  0                                        \n',
 'SHEET    2   D 2 GLU B  61  ASN B  64 -1  O  ASN B  64   N  GLN B  52           \n',
 'CRYST1   49.942   51.699   82.120  90.00  90.00  90.00 P 21 21 21    8          \n']
        
        correct_header = {
       'bio_cmx': [[[('A',), ('B',)], 1]], 
        'uc_mxs': array([[[  1.    ,   0.    ,   0.    ,   0.    ],\
        [  0.    ,   1.    ,   0.    ,   0.    ],\
        [  0.    ,   0.    ,   1.    ,   0.    ],\
        [  0.    ,   0.    ,   0.    ,   1.    ]],\

       [[ -1.    ,   0.    ,   0.    ,  24.971 ],\
        [  0.    ,  -1.    ,   0.    ,   0.    ],\
        [  0.    ,   0.    ,   1.    ,  41.06  ],\
        [  0.    ,   0.    ,   0.    ,   1.    ]],\

       [[ -1.    ,   0.    ,   0.    ,   0.    ],\
        [  0.    ,   1.    ,   0.    ,  25.8495],\
        [  0.    ,   0.    ,  -1.    ,  41.06  ],\
        [  0.    ,   0.    ,   0.    ,   1.    ]],\

       [[  1.    ,   0.    ,   0.    ,  24.971 ],\
        [  0.    ,  -1.    ,   0.    ,  25.8495],\
        [  0.    ,   0.    ,  -1.    ,   0.    ],\
        [  0.    ,   0.    ,   0.    ,   1.    ]]]), \
        'dbref_acc_full': 'Q8P4R5_XANCP', \
        'name': 'TRANSLATION', \
        'solvent_content': '46.26', \
        'expdta': 'X-RAY', \
        'bio_mxs': array([[[ 1.,  0.,  0.,  0.],\
        [ 0.,  1.,  0.,  0.],\
        [ 0.,  0.,  1.,  0.],\
        [ 0.,  0.,  0.,  1.]]]), 
        'uc_omx': array([[ 49.94256605,   0.        ,   0.        ],\
       [  0.        ,  51.69828879,   0.        ],\
       [  0.        ,   0.        ,  82.12203334]]), \
       'space_group': 'P 21 21 21', 'r_free': '0.280', \
       'cryst1': '49.942   51.699   82.120  90.00  90.00  90.00', \
       'experiment_type': 'X-RAY DIFFRACTION', \
       'uc_fmx': array([[ 0.020023,  0.      ,  0.      ],\
       [ 0.      ,  0.019343,  0.      ],\
       [ 0.      ,  0.      ,  0.012177]]),\
        'date': '17-OCT-06', \
        'matthews': '2.29', \
        'resolution': '1.70', \
        'id': '2E12', \
        'dbref_acc': 'Q8P4R5'}
    
        parsed_header = parse_header(header)
        for key, val in parsed_header.items():
            assert val == correct_header[key]
    
    def test_get_trailer_offset(self):
        lines = ['ATOM','CONNECT']
        assert get_trailer_offset(lines) == 1
        
    def test_get_coords_offset(self):
        lines = ['dummy','ATOM','CONNECT']
        assert get_coords_offset(lines) == 1
            
    def test_get_symmetry(self):
        """testing parsing of symmetry operators
        """
        lines = [
 'REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000            \n',
 'REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000            \n',
 'REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000            \n',
 'REMARK 290   SMTRY1   2 -1.000000  0.000000  0.000000       24.97100            \n',
 'REMARK 290   SMTRY2   2  0.000000 -1.000000  0.000000        0.00000            \n',
 'REMARK 290   SMTRY3   2  0.000000  0.000000  1.000000       41.06000            \n',
 'REMARK 290   SMTRY1   3 -1.000000  0.000000  0.000000        0.00000            \n',
 'REMARK 290   SMTRY2   3  0.000000  1.000000  0.000000       25.84950            \n',
 'REMARK 290   SMTRY3   3  0.000000  0.000000 -1.000000       41.06000            \n',
 'REMARK 290   SMTRY1   4  1.000000  0.000000  0.000000       24.97100            \n',
 'REMARK 290   SMTRY2   4  0.000000 -1.000000  0.000000       25.84950            \n',
 'REMARK 290   SMTRY3   4  0.000000  0.000000 -1.000000        0.00000            \n',
 'REMARK 290                                                                      \n',
 'REMARK 290 REMARK: NULL                                                         \n',
 'REMARK 350                                                                      \n',
 'REMARK 350 COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN           \n',
 'REMARK 350 BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE                \n',
 'REMARK 350 MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS          \n',
 'REMARK 350 GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND                          \n',
 'REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.                               \n',
 'REMARK 350                                                                      \n',
 'REMARK 350 BIOMOLECULE: 1                                                       \n',
 'REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC                           \n',
 'REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B                                  \n',
 'REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000            \n',
 'REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000            \n',
 'REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000            \n',
 'CRYST1   49.942   51.699   82.120  90.00  90.00  90.00 P 21 21 21    8          \n',
 'ORIGX1      1.000000  0.000000  0.000000        0.00000                         \n',
 'ORIGX2      0.000000  1.000000  0.000000        0.00000                         \n',
 'ORIGX3      0.000000  0.000000  1.000000        0.00000                         \n',
 'SCALE1      0.020023  0.000000  0.000000        0.00000                         \n',
 'SCALE2      0.000000  0.019343  0.000000        0.00000                         \n',
 'SCALE3      0.000000  0.000000  0.012177        0.00000                         \n']
        sym = get_symmetry(lines)
        correct_sym = {
        'bio_cmx': [[[('A',), ('B',)], 1]], 
        'uc_mxs': array([[[  1.    ,   0.    ,   0.    ,   0.    ],\
        [  0.    ,   1.    ,   0.    ,   0.    ],\
        [  0.    ,   0.    ,   1.    ,   0.    ],\
        [  0.    ,   0.    ,   0.    ,   1.    ]],\

        [[ -1.    ,   0.    ,   0.    ,  24.971 ],\
        [  0.    ,  -1.    ,   0.    ,   0.    ],\
        [  0.    ,   0.    ,   1.    ,  41.06  ],\
        [  0.    ,   0.    ,   0.    ,   1.    ]],\

        [[ -1.    ,   0.    ,   0.    ,   0.    ],\
        [  0.    ,   1.    ,   0.    ,  25.8495],\
        [  0.    ,   0.    ,  -1.    ,  41.06  ],\
        [  0.    ,   0.    ,   0.    ,   1.    ]],\

        [[  1.    ,   0.    ,   0.    ,  24.971 ],\
        [  0.    ,  -1.    ,   0.    ,  25.8495],\
        [  0.    ,   0.    ,  -1.    ,   0.    ],\
        [  0.    ,   0.    ,   0.    ,   1.    ]]]), \
        'bio_mxs': array([[[ 1.,  0.,  0.,  0.],\
        [ 0.,  1.,  0.,  0.],\
        [ 0.,  0.,  1.,  0.],\
        [ 0.,  0.,  0.,  1.]]]), 
        'uc_omx': array([[ 49.94256605,   0.        ,   0.        ],\
        [  0.        ,  51.69828879,   0.        ],\
        [  0.        ,   0.        ,  82.12203334]]), \
        'uc_fmx': array([[ 0.020023,  0.      ,  0.      ],\
        [ 0.      ,  0.019343,  0.      ],\
        [ 0.      ,  0.      ,  0.012177]]),}   
        
        for key in sym:
            try:
                assert sym[key] == correct_sym[key]
            except ValueError:
                assert allclose(sym[key], correct_sym[key])
            
    def test_dict2pdb(self):
        """testing pdb dict round-trip.
        """
        line = 'ATOM      1  N   MET A   1      53.045  42.225  33.724  1.00  2.75           N\n'
        d = pdb2dict(line)
        line2 = dict2pdb(d)
        assert line == line2
        d.pop('coords')
        assert d == {'ser_num': 1, 'res_long_id': ('MET', 1, ' '), 
                     'h_flag': ' ', 
                     'at_name': ' N  ', 
                     'at_long_id': ('N', ' '), 
                     'bfactor': 2.75, 'chain_id': 'A', 
                     'occupancy': 1.0, 'element': ' N', 
                     'res_name': 'MET',
                     'seg_id': '    ', 'at_id': 'N', 
                     'alt_loc': ' ', 
                     'res_ic': ' ', 
                     'res_id': 1, 
                     'at_type': 'ATOM  '}

    def test_dict2ter(self):
        d = {'ser_num': 1, 'chain_id': 'A', 'res_name': 'MET', 'res_ic': ' ', \
              'res_id': 1,}
        assert dict2ter(d) == 'TER       2      MET A   1 \n'
if __name__ == '__main__':
    main()
