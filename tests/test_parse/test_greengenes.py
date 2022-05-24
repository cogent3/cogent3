#!/usr/bin/env python

from unittest import TestCase, main

from cogent3.parse.greengenes import (
    MinimalGreengenesParser,
    SpecificGreengenesParser,
    make_ignore_f,
)


__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2022, The Cogent Project"  # consider project name
# remember to add yourself if you make changes
__credits__ = ["Daniel McDonald"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Prototype"


class ParseGreengenesRecordsTests(TestCase):
    def setUp(self):
        pass

    def test_MinimalGreengenesParser_mock(self):
        """Test MinimalGreengenesParser against mock data"""
        res = MinimalGreengenesParser(
            mock_data.splitlines(), RecStart="my_starting", RecEnd="my_ending"
        )

        records = list(res)

        exp = [
            {"a": "1", "b": "2", "c": "3", "d": "", "e": "5"},
            {"q": "asdasd", "c": "taco"},
        ]

        self.assertEqual(records, exp)

    def test_MinimalGreengenesParser_real(self):
        """Test MinimalGreengenesParser against real data"""
        res = MinimalGreengenesParser(real_data.splitlines())
        record1, record2 = list(res)

        self.assertEqual(record1["G2_chip_tax_string"], "Unclassified")
        self.assertEqual(
            record1["authors"],
            "Hernanandez-Eugenio,G., Silva-Rojas,H.V., Zelaya-Molina,L.X.",
        )
        self.assertEqual(record1["bel3_div_ratio"], "")
        self.assertEqual(len(record1), 72)

        self.assertEqual(record2["ncbi_acc_w_ver"], "FJ832719.1")
        self.assertEqual(record2["timestamp"], "2010-03-23 14:08:27")
        self.assertEqual(
            record2["title"],
            "Developmental Microbial Ecology of the Crop of the Folivorous Hoatzin",
        )

    def test_SpecificGreengenesParser_real(self):
        """Test SpecificGreengenesParser against real data"""
        fields = ["prokMSA_id", "journal"]
        res = SpecificGreengenesParser(real_data.splitlines(), fields)
        records = list(res)
        exp = [("604868", ""), ("604867", "ISME J (2010) In press")]
        self.assertEqual(records, exp)

        ids = ["604867", "12312312323"]
        res = SpecificGreengenesParser(real_data.splitlines(), fields, ids)
        records = list(res)
        exp = [("604867", "ISME J (2010) In press")]
        self.assertEqual(records, exp)

    def test_make_ignore_f(self):
        """Properly ignore empty records and the start line"""
        f = make_ignore_f("testing")
        self.assertFalse(f(["asasdasd", ""]))
        self.assertFalse(f(["test", ""]))
        self.assertFalse(f(["testing2", ""]))
        self.assertFalse(f(["testing", "asd"]))
        self.assertTrue(f(["", ""]))
        self.assertTrue(f(None))
        self.assertTrue(f(["", ""]))
        self.assertTrue(f(["testing", ""]))


mock_data = """my_starting
a=1
b=2
c=3
d=
e=5
my_ending

my_starting
q=asdasd
c=taco
my_ending
"""

real_data = """BEGIN
G2_chip_tax_string=Unclassified
G2_chip_tax_string_format_2=Unclassified
HOMD_tax_string=
HOMD_tax_string_format_2=
Hugenholtz_tax_string=Unclassified
Hugenholtz_tax_string_format_2=Unclassified
Ludwig_tax_string=Unclassified
Ludwig_tax_string_format_2=Unclassified
Pace_tax_string=Unclassified
Pace_tax_string_format_2=Unclassified
RDP_tax_string=Unclassified
RDP_tax_string_format_2=Unclassified
Silva_tax_string=Unclassified
Silva_tax_string_format_2=Unclassified
authors=Hernanandez-Eugenio,G., Silva-Rojas,H.V., Zelaya-Molina,L.X.
bel3_div_ratio=
bellerophon=
blast_perc_ident_to_template=
clone=51a 
contact_info=Irrigacion, Universidad Autonoma Chapingo, Carretera Mexico-Texcoco Km 37.5, Texcoco, Mexico 56230, Mexico
core_set_member=
core_set_member2=
country=Mexico: Mexico City 
create_date=21-NOV-2009
db_name=
decision=clone
description=Uncultured bacterium clone 51a 16S ribosomal RNA gene, partial sequence
email=
gold_id=
img_oid=
isolate=
isolation_source=mesophilic anaerobic reactor fed with effluent from the chemical industry 
journal=
longest_insertion=
medline_ids=
ncbi_acc=
ncbi_acc_w_ver=FJ461956.1
ncbi_gi=213390944
ncbi_seq_length=1512
ncbi_tax_id=77133
ncbi_tax_string=Bacteria; environmental samples
ncbi_tax_string_format_2=Unclassified
non_ACGT_count=
non_ACGT_percent=
note=
organism=uncultured bacterium
perc_ident_to_invariant_core=
prokMSA_id=604868
prokMSAname=Microbial ecology industrial digestor mesophilic anaerobic reactor fed effluent chemical industry clone 51a
pubmed_ids=
remark=
replaced_by=
single_nt_runs_over_7=
small_gap_intrusions=
source=uncultured bacterium
span_aligned=1..2
specific_host=
status=0
strain=
study_id=38002
sub_species=
submit_date=24-OCT-2008
template=
timestamp=2010-03-23 14:08:27
title=Microbial ecology of industrial anaerobic digestor
unaligned_length=
update_date=21-NOV-2009
warning=
wigeon95=
wigeon99=
wigeon_std_dev=
aligned_seq=unaligned
END

BEGIN
G2_chip_tax_string=Unclassified
G2_chip_tax_string_format_2=Unclassified
HOMD_tax_string=
HOMD_tax_string_format_2=
Hugenholtz_tax_string=Unclassified
Hugenholtz_tax_string_format_2=Unclassified
Ludwig_tax_string=Unclassified
Ludwig_tax_string_format_2=Unclassified
Pace_tax_string=Unclassified
Pace_tax_string_format_2=Unclassified
RDP_tax_string=Unclassified
RDP_tax_string_format_2=Unclassified
Silva_tax_string=Unclassified
Silva_tax_string_format_2=Unclassified
authors=Brodie,E.L., Dominguez-Bello,M.G., Garcia-Amado,M.A., Godoy-Vitorino,F., Goldfarb,K.C., Michelangeli,F.
bel3_div_ratio=
bellerophon=
blast_perc_ident_to_template=
clone=J3Q101_11C02 
contact_info=Biology, University of Puerto Rico, Rio Piedras Campus, PO Box 23360, San Juan, PR 00931-3360, USA
core_set_member=
core_set_member2=
country=Venezuela 
create_date=10-DEC-2009
db_name=
decision=clone
description=Uncultured bacterium clone J3Q101_11C02 16S ribosomal RNA gene, partial sequence
email=
gold_id=
img_oid=
isolate=
isolation_source=crop contents 
journal=ISME J (2010) In press
longest_insertion=
medline_ids=
ncbi_acc=
ncbi_acc_w_ver=FJ832719.1
ncbi_gi=226447371
ncbi_seq_length=1326
ncbi_tax_id=77133
ncbi_tax_string=Bacteria; environmental samples
ncbi_tax_string_format_2=Unclassified
non_ACGT_count=
non_ACGT_percent=
note=
organism=uncultured bacterium
perc_ident_to_invariant_core=
prokMSA_id=604867
prokMSAname=Microbial Ecology Crop Folivorous Hoatzin crop contents clone J3Q101_11C02
pubmed_ids=
remark=
replaced_by=
single_nt_runs_over_7=
small_gap_intrusions=
source=uncultured bacterium
span_aligned=1..2
specific_host=
status=0
strain=
study_id=37901
sub_species=
submit_date=16-MAR-2009
template=
timestamp=2010-03-23 14:08:27
title=Developmental Microbial Ecology of the Crop of the Folivorous Hoatzin
unaligned_length=
update_date=10-DEC-2009
warning=
wigeon95=
wigeon99=
wigeon_std_dev=
aligned_seq=unaligned
END
"""

if __name__ == "__main__":
    main()
