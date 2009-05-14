Querying NCBI for VWF
=====================

This example is taken from the PyCogent paper (Knight et al. Genome Biol, 8(8):R171, 2007). **Note:** Due to changes by NCBI in the structure of there records, it is no longer easy (but not impossible) to use the Sequences ``Info`` attribute to restrict sequences to Swissprot entries as we had done previously.

We query the NCBI protein data base for von Willebrand Factor (VWF), a 2813 amino acid glycoprotein required for platelet adhesion in blood coagulation. Missense mutations in this molecule have been associated with von Willebrand disease, a heterogeneous disorder characterized by prolonged bleeding.

We import ``EUtils`` for querying NCBI and search the protein data-base, restricting our search to mammal sequences.

.. pycode::
    
    >>> from cogent.db.ncbi import EUtils
    >>> db = EUtils(db="protein", rettype="gp")
    >>> query = '"VWf"[gene] AND Mammalia[orgn]'
    >>> records = db[query].readlines()

We have requested the GenBank record format. We use the ``RichGenbankParser`` to grab features from the feature annotations of this format. We illustrate grabbing additional content from the feature tables by extracting the locations of SNPs and whether those SNPs have a disease association. An inspection of the feature tables for human entries reveals that the Swissprot entry had the most complete and accessible content regarding the SNP data: region features with ``region_name="Variant"`` and a note field that contains the indicated amino acid difference along with an indication of whether the SNP was associated with von Willebrand disease (denoted by the symbol VWD). Finally, we seek to extract the protein domain locations for presentation purposes.

We simply limit our attention to large, relatively complete sequences. From the human record we extract annotation data. We use regular expressions to assist with extracting data regarding the amino acid change and the domain names. We also name sequences based on their [G]enus [spe]cies and accession.

The selected species are accumulated in a ``seqs`` dictionary, keyed by their name. The feature data are accumulated in a list.

.. pycode::
    
    >>> import re
    >>> from cogent.parse.genbank import RichGenbankParser
    >>> parser = RichGenbankParser(records)
    >>> seqs = {}
    >>> rows = []
    >>> for accession, seq in parser:
    ...     if len(seq) < 2800:
    ...         continue
    ...     # we extract annotation data only from the human record
    ...     if "Homo" in seq.Info.species:
    ...         for feature in seq.Info.features:
    ...             if "region_name" not in feature:
    ...                 continue # ignore this one, go to the next feature
    ...             if "Variant" in feature["region_name"]:
    ...                 note = feature["note"][0]
    ...                 variant = " ".join(re.findall("[A-Z] *-> *[A-Z]",
    ...                                               note))
    ...                 disease = "VWD" in note
    ...                 lo = feature["location"].first() - 1
    ...                 hi = feature["location"].last()
    ...                 rows.append(["SNP", lo, hi, variant, disease])
    ...             else:
    ...                 region_name = feature["region_name"][0]
    ...                 if region_name == "Domain":
    ...                     note = [field.strip() \
    ...         for field in re.split("[.;]", feature["note"][0]) if field]
    ...                     if len(note) == 1:
    ...                         note += [""]
    ...                     lo, hi = feature["location"].first() - 1, feature["location"].last()
    ...                     rows.append(["Domain", lo, hi, ' '.join(note), None])
    ...     species = seq.Info.species.split()
    ...     seq_name = "%s.%s" % (species[0][0] + species[1][:3], accession)
    ...     seqs[seq_name] = seq

We convert the sequences to a ``SequenceCollection`` using ``LoadSeqs`` and then save to a fasta formatted file.

.. pycode::
    
    >>> from cogent import LoadSeqs
    >>> seqs = LoadSeqs(data=seqs, aligned=False)
    >>> print seqs.NamedSeqs['Clup.Q28295'].toFasta()
    >Clup.Q28295
    MSPTRLVRVLLALALI...

We convert the features into a PyCogent ``Table`` object, which requires we specify column headings. This can be saved to file if desired, but we don't do that here.
    
.. pycode::
    
    >>> from cogent import LoadTable
    >>> feature_table = LoadTable(header=["Type", "Start", "Stop", "Note",
    ...                    "Disease"], rows=rows)
    >>> print feature_table
    ===========================================================================================
      Type    Start    Stop                                                     Note    Disease
    -------------------------------------------------------------------------------------------
    Domain       33     240                                                  VWFD 1            
       SNP      272     273                                                   R -> W       True
    Domain      294     348                                                   TIL 1            
       SNP      376     377                                                   W -> C       True
    Domain      386     598                                                  VWFD 2            
       SNP      483     484                                                   H -> R      False
       SNP      527     528                                                   N -> S       True
       SNP      549     550                                                   G -> R       True
    Domain      651     707                                                   TIL 2            
    Domain      775     827                                                   TIL 3            
       SNP      787     788                                                   C -> Y       True
       SNP      788     789                                                   T -> A      False
       SNP      790     791                                                   T -> M      False
       SNP      815     816                                                   R -> W      False
       SNP      851     852                                                   R -> Q      False
       SNP      853     854                                                   R -> Q      False
       SNP      856     857                                                   N -> D      False
    Domain      865    1074                                                  VWFD 3            
       SNP     1059    1060                                                   C -> R       True
    Domain     1145    1196                                                   TIL 4            
       SNP     1265    1266                                                   P -> L       True
       SNP     1267    1268                                                   H -> D       True
       SNP     1271    1272                                                   C -> R       True
    Domain     1276    1453         VWFA 1 binding site for platelet glycoprotein Ib           
       SNP     1305    1306                                                   R -> W       True
       SNP     1307    1308                                                   R -> C       True
       SNP     1312    1313                                                   W -> C       True
       SNP     1313    1314                                                   V -> L       True
       SNP     1315    1316                                                   V -> M       True
       SNP     1317    1318                                                   V -> L       True
       SNP     1323    1324                                                   G -> S       True
       SNP     1340    1341                                                   R -> Q       True
       SNP     1373    1374                                                   R -> C       True
       SNP     1373    1374                                                   R -> H       True
       SNP     1380    1381                                                   A -> T      False
       SNP     1398    1399                                                   R -> H      False
       SNP     1459    1460                                                   L -> V       True
       SNP     1460    1461                                                   A -> V       True
       SNP     1471    1472                                                   H -> D      False
    Domain     1497    1665                                                  VWFA 2            
       SNP     1513    1514                                                   F -> C       True
       SNP     1539    1540                                                   L -> P       True
       SNP     1564    1565                                                   V -> L      False
       SNP     1569    1570                                                   Y -> C      False
       SNP     1583    1584                                                   Y -> C      False
       SNP     1596    1597                                                   R -> G       True
       SNP     1596    1597                                                   R -> Q       True
       SNP     1596    1597                                                   R -> W       True
       SNP     1606    1607                                                   V -> D       True
       SNP     1608    1609                                                   G -> R       True
       SNP     1612    1613                                                   S -> P       True
       SNP     1627    1628                                                   I -> T       True
       SNP     1637    1638                                                   E -> K       True
       SNP     1647    1648                                                   P -> S       True
       SNP     1664    1665                                                   V -> E       True
    Domain     1690    1871    VWFA 3 main binding site for collagens type I and III           
    Domain     1948    2153                                                  VWFD 4            
       SNP     2062    2063                                                   P -> S       True
    Domain     2254    2328                                                  VWFC 1            
       SNP     2361    2362                                                   C -> F       True
    Domain     2428    2495                                                  VWFC 2            
       SNP     2545    2546                                                   N -> Y       True
    Domain     2579    2645                                                  VWFC 3            
    Domain     2723    2812                                                    CTCK            
       SNP     2738    2739                                                   C -> Y       True
       SNP     2772    2773                                                   C -> R       True
    -------------------------------------------------------------------------------------------


