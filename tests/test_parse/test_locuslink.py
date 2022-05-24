#!/usr/bin/env python
"""Unit tests for locuslink-specific classes
"""
from unittest import TestCase, main

from cogent3.parse.locuslink import (
    LinesToLocusLink,
    LLFinder,
    _read_accession,
    _read_accnum,
    _read_cdd,
    _read_comp,
    _read_contig,
    _read_extannot,
    _read_go,
    _read_grif,
    _read_map,
    _read_pmid,
    _read_rell,
    _read_sts,
)


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class locuslinkTests(TestCase):
    """Tests toplevel functions."""

    def test_read_accession(self):
        """_read_accession should perform correct conversions"""
        self.assertEqual(
            _read_accession("NP_035835|6755985|na\n"),
            {"Accession": "NP_035835", "Gi": "6755985", "Strain": "na"},
        )
        # check that it ignores additional fields
        self.assertEqual(
            _read_accession("NG_002740|30172554|na|1|1315\n"),
            {"Accession": "NG_002740", "Gi": "30172554", "Strain": "na"},
        )

    def test_read_rell(self):
        """_read_rell should perform correct conversions"""
        self.assertEqual(
            _read_rell("related mRNA|AK090391|n|NM_153775--AK090391\n"),
            {
                "Description": "related mRNA",
                "Id": "AK090391",
                "IdType": "n",
                "Printable": "NM_153775--AK090391",
            },
        )

    def test_read_accnum(self):
        """_read_accnum should perform correct conversions"""
        self.assertEqual(
            _read_accnum("NG_002740|30172554|na|1|1315\n"),
            {
                "Accession": "NG_002740",
                "Gi": "30172554",
                "Strain": "na",
                "Start": "1",
                "End": "1315",
            },
        )

    def test_read_map(self):
        """_read_map should perform correct conversions"""
        self.assertEqual(
            _read_map("10 C1|RefSeq|C|\n"),
            {"Location": "10 C1", "Source": "RefSeq", "Type": "C"},
        )

    def test_read_sts(self):
        """_read_sts should perform correct conversions"""
        self.assertEqual(
            _read_sts("RH35858|2|37920|na|seq_map|epcr\n"),
            {
                "Name": "RH35858",
                "Chromosome": "2",
                "StsId": "37920",
                "Segment": "na",
                "SequenceKnown": "seq_map",
                "Evidence": "epcr",
            },
        )

    def test_read_cdd(self):
        """_read_cdd should perform correct conversions"""
        self.assertEqual(
            _read_cdd("Immunoglobulin C-2 Type|smart00408|103|na|4.388540e+01\n"),
            {
                "Name": "Immunoglobulin C-2 Type",
                "Key": "smart00408",
                "Score": "103",
                "EValue": "na",
                "BitScore": "4.388540e+01",
            },
        )

    def test_read_comp(self):
        """_read_comp should perform correct conversions"""
        self.assertEqual(
            _read_comp("10090|Map2k6|11|11  cM|26399|17|MAP2K6|ncbi_mgd\n"),
            {
                "TaxonId": "10090",
                "Symbol": "Map2k6",
                "Chromosome": "11",
                "Position": "11  cM",
                "LocusId": "26399",
                "ChromosomeSelf": "17",
                "SymbolSelf": "MAP2K6",
                "MapName": "ncbi_mgd",
            },
        )

    def test_read_grif(self):
        """_read_grif should perform correct conversions"""
        self.assertEqual(
            _read_grif("12037672|interaction with pRb\n"),
            {"PubMedId": "12037672", "Description": "interaction with pRb"},
        )

    def test_read_pmid(self):
        """_read_pmid should perform correct conversions"""
        self.assertEqual(
            _read_pmid("12875969,12817023,12743034\n"),
            ["12875969", "12817023", "12743034"],
        )

    def test_read_go(self):
        """_read_go should perform correct conversions"""
        self.assertEqual(
            _read_go("molecular function|zinc ion binding|IEA|GO:0008270|GOA|na\n"),
            {
                "Category": "molecular function",
                "Term": "zinc ion binding",
                "EvidenceCode": "IEA",
                "GoId": "GO:0008270",
                "Source": "GOA",
                "PubMedId": "na",
            },
        )

    def test_read_extannot(self):
        """_read_extannot should perform correct conversions"""
        self.assertEqual(
            _read_extannot("cellular role|Pol II transcription|NR|Proteome|8760285\n"),
            {
                "Category": "cellular role",
                "Term": "Pol II transcription",
                "EvidenceCode": "NR",
                "Source": "Proteome",
                "PubMedId": "8760285",
            },
        )

    def test_read_contig(self):
        """_read_contig should perform correct conversions"""
        self.assertEqual(
            _read_contig("NT_011109.15|29800594|na|31124734|31133047|-|19|reference\n"),
            {
                "Accession": "NT_011109.15",
                "Gi": "29800594",
                "Strain": "na",
                "From": "31124734",
                "To": "31133047",
                "Orientation": "-",
                "Chromosome": "19",
                "Assembly": "reference",
            },
        )

    def test_LinesToLocusLink(self):
        """LinesToLocusLink should give expected results on sample data"""
        fake_file = """>>1
LOCUSID: 1
LOCUS_CONFIRMED: yes
LOCUS_TYPE: gene with protein product, function known or inferred
ORGANISM: Homo sapiens
STATUS: REVIEWED
NM: NM_130786|21071029|na
NP: NP_570602|21071030
CDD: Immunoglobulin C-2 Type|smart00408|103|na|4.388540e+01
PRODUCT: alpha 1B-glycoprotein
ASSEMBLY: AF414429,AK055885,AK056201
CONTIG: NT_011109.15|29800594|na|31124734|31133047|-|19|reference
EVID: supported by alignment with mRNA
XM: NM_130786|21071029|na
XP: NP_570602|21071030|na
ACCNUM: AC010642|9929687|na|43581|41119
TYPE: g
ACCNUM: AF414429|15778555|na|na|na
TYPE: m
PROT: AAL07469|15778556
ACCNUM: AK055885|16550723|na|na|na
TYPE: m
ACCNUM: AK056201|16551539|na|na|na
TYPE: m
ACCNUM: BC035719|23273475|na|na|na
TYPE: m
PROT: AAH35719|23273476
ACCNUM: none|na|na|na|na
TYPE: p
PROT: P04217|23503038
OFFICIAL_SYMBOL: A1BG
OFFICIAL_GENE_NAME: alpha-1-B glycoprotein
ALIAS_SYMBOL: A1B
ALIAS_SYMBOL: ABG
ALIAS_SYMBOL: GAB
PREFERRED_PRODUCT: alpha 1B-glycoprotein
SUMMARY: Summary: The protein encoded by this gene is a plasma glycoprotein of unknown function. The protein shows sequence similarity to the variable regions of some immunoglobulin supergene family member proteins.
CHR: 19
STS: RH65092|-|10673|na|na|epcr
STS: WI-16009|-|52209|na|na|epcr
STS: G59506|-|136670|na|na|epcr
COMP: 10090|A1bg|na|na|117586|19|A1BG|ncbi_mgd
COMP: 10090|A1bg|7|7  cM|117586|19|A1BG|ncbi_mgd
BUTTON: unigene.gif
LINK: http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Hs&CID=390608
UNIGENE: Hs.390608
OMIM: 138670
MAP: 19q13.4|RefSeq|C|
MAPLINK: default_human_gene|A1BG
BUTTON: snp.gif
LINK: http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?locusId=1
BUTTON: homol.gif
LINK: http://www.ncbi.nlm.nih.gov/HomoloGene/homolquery.cgi?TEXT=1[loc]&TAXID=9606
BUTTON: ensembl.gif
LINK: http://www.ensembl.org/Homo_sapiens/contigview?geneid=NM_130786
BUTTON: ucsc.gif
LINK: http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&position=NM_130786
BUTTON: mgc.gif
LINK: http://mgc.nci.nih.gov/Genes/GeneInfo?ORG=Hs&CID=390608
PMID: 12477932,8889549,3458201,2591067
GO: molecular function|molecular_function unknown|ND|GO:0005554|GOA|3458201
GO: biological process|biological_process unknown|ND|GO:0000004|GOA|na
GO: cellular component|extracellular|IDA|GO:0005576|GOA|3458201
>>386590
LOCUSID: 386590
LOCUS_CONFIRMED: yes
LOCUS_TYPE: gene with protein product, function known or inferred
ORGANISM: Danio rerio
ACCNUM: AF510108|31323727|na|na|na
TYPE: m
PROT: AAP47138|31323728
OFFICIAL_SYMBOL: tra1
OFFICIAL_GENE_NAME: tumor rejection antigen (gp96) 1
BUTTON: zfin.gif
LINK: http://zfin.org/cgi-bin/ZFIN_jump?record=ZDB-GENE-031002-1
PMID: 14499652"""

        records = list(LLFinder(fake_file.split("\n")))
        self.assertEqual(len(records), 2)
        first, second = list(map(LinesToLocusLink, records))

        # test the second one first, since it's shorter
        self.assertEqual(second.LOCUSID, 386590)
        self.assertEqual(second.LOCUS_CONFIRMED, "yes")
        self.assertEqual(
            second.LOCUS_TYPE, "gene with protein product, function known or inferred"
        )
        self.assertEqual(second.ORGANISM, "Danio rerio")
        self.assertEqual(
            second.ACCNUM,
            [
                {
                    "Accession": "AF510108",
                    "Gi": "31323727",
                    "Strain": "na",
                    "Start": "na",
                    "End": "na",
                }
            ],
        )
        self.assertEqual(second.TYPE, ["m"])
        self.assertEqual(second.PROT, [{"Accession": "AAP47138", "Gi": "31323728"}])
        self.assertEqual(second.OFFICIAL_SYMBOL, "tra1")
        self.assertEqual(second.OFFICIAL_GENE_NAME, "tumor rejection antigen (gp96) 1")
        self.assertEqual(second.BUTTON, ["zfin.gif"])
        self.assertEqual(
            second.LINK, ["http://zfin.org/cgi-bin/ZFIN_jump?record=ZDB-GENE-031002-1"]
        )
        self.assertEqual(second.PMID, ["14499652"])

        # now for the annoying test on the longer record
        self.assertEqual(first.LOCUSID, 1)
        self.assertEqual(first.LOCUS_CONFIRMED, "yes")
        self.assertEqual(first.ORGANISM, "Homo sapiens")
        self.assertEqual(
            first.LOCUS_TYPE, "gene with protein product, function known or inferred"
        )
        self.assertEqual(first.STATUS, "REVIEWED")
        self.assertEqual(
            first.NM, [{"Accession": "NM_130786", "Gi": "21071029", "Strain": "na"}]
        )
        self.assertEqual(first.NP, [{"Accession": "NP_570602", "Gi": "21071030"}])
        self.assertEqual(
            first.CDD,
            [
                {
                    "Name": "Immunoglobulin C-2 Type",
                    "Key": "smart00408",
                    "Score": "103",
                    "EValue": "na",
                    "BitScore": "4.388540e+01",
                }
            ],
        )
        self.assertEqual(first.PRODUCT, ["alpha 1B-glycoprotein"])
        self.assertEqual(first.ASSEMBLY, [["AF414429", "AK055885", "AK056201"]])
        self.assertEqual(
            first.CONTIG,
            [
                {
                    "Accession": "NT_011109.15",
                    "Gi": "29800594",
                    "Strain": "na",
                    "From": "31124734",
                    "To": "31133047",
                    "Orientation": "-",
                    "Chromosome": "19",
                    "Assembly": "reference",
                }
            ],
        )
        self.assertEqual(first.EVID, ["supported by alignment with mRNA"])
        self.assertEqual(
            first.XM, [{"Accession": "NM_130786", "Gi": "21071029", "Strain": "na"}]
        )
        self.assertEqual(
            first.XP, [{"Accession": "NP_570602", "Gi": "21071030", "Strain": "na"}]
        )
        self.assertEqual(
            first.ACCNUM,
            [
                {
                    "Accession": "AC010642",
                    "Gi": "9929687",
                    "Strain": "na",
                    "Start": "43581",
                    "End": "41119",
                },
                {
                    "Accession": "AF414429",
                    "Gi": "15778555",
                    "Strain": "na",
                    "Start": "na",
                    "End": "na",
                },
                {
                    "Accession": "AK055885",
                    "Gi": "16550723",
                    "Strain": "na",
                    "Start": "na",
                    "End": "na",
                },
                {
                    "Accession": "AK056201",
                    "Gi": "16551539",
                    "Strain": "na",
                    "Start": "na",
                    "End": "na",
                },
                {
                    "Accession": "BC035719",
                    "Gi": "23273475",
                    "Strain": "na",
                    "Start": "na",
                    "End": "na",
                },
                {
                    "Accession": "none",
                    "Gi": "na",
                    "Strain": "na",
                    "Start": "na",
                    "End": "na",
                },
            ],
        )
        self.assertEqual(first.TYPE, ["g", "m", "m", "m", "m", "p"])
        self.assertEqual(
            first.PROT,
            [
                {"Accession": "AAL07469", "Gi": "15778556"},
                {"Accession": "AAH35719", "Gi": "23273476"},
                {"Accession": "P04217", "Gi": "23503038"},
            ],
        )
        self.assertEqual(first.OFFICIAL_SYMBOL, "A1BG")
        self.assertEqual(first.OFFICIAL_GENE_NAME, "alpha-1-B glycoprotein")
        self.assertEqual(first.ALIAS_SYMBOL, ["A1B", "ABG", "GAB"])
        self.assertEqual(first.PREFERRED_PRODUCT, ["alpha 1B-glycoprotein"])
        self.assertEqual(
            first.SUMMARY,
            [
                """Summary: The protein encoded by this gene is a plasma glycoprotein of unknown function. The protein shows sequence similarity to the variable regions of some immunoglobulin supergene family member proteins."""
            ],
        )
        self.assertEqual(first.CHR, ["19"])
        self.assertEqual(
            first.STS,
            [
                {
                    "Name": "RH65092",
                    "Chromosome": "-",
                    "StsId": "10673",
                    "Segment": "na",
                    "SequenceKnown": "na",
                    "Evidence": "epcr",
                },
                {
                    "Name": "WI-16009",
                    "Chromosome": "-",
                    "StsId": "52209",
                    "Segment": "na",
                    "SequenceKnown": "na",
                    "Evidence": "epcr",
                },
                {
                    "Name": "G59506",
                    "Chromosome": "-",
                    "StsId": "136670",
                    "Segment": "na",
                    "SequenceKnown": "na",
                    "Evidence": "epcr",
                },
            ],
        )
        self.assertEqual(
            first.COMP,
            [
                {
                    "TaxonId": "10090",
                    "Symbol": "A1bg",
                    "Chromosome": "na",
                    "Position": "na",
                    "LocusId": "117586",
                    "ChromosomeSelf": "19",
                    "SymbolSelf": "A1BG",
                    "MapName": "ncbi_mgd",
                },
                {
                    "TaxonId": "10090",
                    "Symbol": "A1bg",
                    "Chromosome": "7",
                    "Position": "7  cM",
                    "LocusId": "117586",
                    "ChromosomeSelf": "19",
                    "SymbolSelf": "A1BG",
                    "MapName": "ncbi_mgd",
                },
            ],
        )
        self.assertEqual(
            first.BUTTON,
            [
                "unigene.gif",
                "snp.gif",
                "homol.gif",
                "ensembl.gif",
                "ucsc.gif",
                "mgc.gif",
            ],
        )
        self.assertEqual(
            first.LINK,
            [
                "http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Hs&CID=390608",
                "http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?locusId=1",
                "http://www.ncbi.nlm.nih.gov/HomoloGene/homolquery.cgi?TEXT=1[loc]&TAXID=9606",
                "http://www.ensembl.org/Homo_sapiens/contigview?geneid=NM_130786",
                "http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&position=NM_130786",
                "http://mgc.nci.nih.gov/Genes/GeneInfo?ORG=Hs&CID=390608",
            ],
        )
        self.assertEqual(first.UNIGENE, ["Hs.390608"])
        self.assertEqual(first.OMIM, ["138670"])
        self.assertEqual(
            first.MAP, [{"Location": "19q13.4", "Source": "RefSeq", "Type": "C"}]
        )
        self.assertEqual(first.MAPLINK, ["default_human_gene|A1BG"])
        self.assertEqual(first.PMID, ["12477932", "8889549", "3458201", "2591067"])
        self.assertEqual(
            first.GO,
            [
                {
                    "Category": "molecular function",
                    "Term": "molecular_function unknown",
                    "EvidenceCode": "ND",
                    "GoId": "GO:0005554",
                    "Source": "GOA",
                    "PubMedId": "3458201",
                },
                {
                    "Category": "biological process",
                    "Term": "biological_process unknown",
                    "EvidenceCode": "ND",
                    "GoId": "GO:0000004",
                    "Source": "GOA",
                    "PubMedId": "na",
                },
                {
                    "Category": "cellular component",
                    "Term": "extracellular",
                    "EvidenceCode": "IDA",
                    "GoId": "GO:0005576",
                    "Source": "GOA",
                    "PubMedId": "3458201",
                },
            ],
        )


if __name__ == "__main__":
    main()
