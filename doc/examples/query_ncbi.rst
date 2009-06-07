Querying NCBI for VWF
=====================

.. sectionauthor:: Gavin Huttley

This example is taken from the PyCogent paper (Knight et al. Genome Biol, 8(8):R171, 2007).

.. note:: Due to changes by NCBI in the structure of there records, it is no longer easy (but not impossible) to use the Sequences ``Info`` attribute to restrict sequences to Swissprot entries as we had done previously.

We query the NCBI protein data base for von Willebrand Factor (VWF), a 2813 amino acid glycoprotein required for platelet adhesion in blood coagulation. Missense mutations in this molecule have been associated with von Willebrand disease, a heterogeneous disorder characterized by prolonged bleeding.

We import ``EUtils`` for querying NCBI and search the protein data-base, restricting our search to mammal sequences.

.. doctest::
    
    >>> from cogent.db.ncbi import EUtils
    >>> db = EUtils(db="protein", rettype="gp")
    >>> query = '"VWf"[gene] AND Mammalia[orgn]'
    >>> records = db[query].readlines()

We have requested the GenBank record format. We use the ``RichGenbankParser`` to grab features from the feature annotations of this format. We illustrate grabbing additional content from the feature tables by extracting the locations of SNPs and whether those SNPs have a disease association. An inspection of the feature tables for human entries reveals that the Swissprot entry had the most complete and accessible content regarding the SNP data: region features with ``region_name="Variant"`` and a note field that contains the indicated amino acid difference along with an indication of whether the SNP was associated with von Willebrand disease (denoted by the symbol VWD). Finally, we seek to extract the protein domain locations for presentation purposes.

We simply limit our attention to large, relatively complete sequences. From the human record we extract annotation data. We use regular expressions to assist with extracting data regarding the amino acid change and the domain names. We also name sequences based on their [G]enus [spe]cies and accession.

The selected species are accumulated in a ``seqs`` dictionary, keyed by their name. The feature data are accumulated in a list.

.. doctest::
    
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
    ...                 rows.append(["SNP", lo, hi, variant.strip(), disease])
    ...             else:
    ...                 region_name = feature["region_name"][0]
    ...                 if region_name == "Domain":
    ...                     note = [field.strip() \
    ...         for field in re.split("[.;]", feature["note"][0]) if field]
    ...                     if len(note) == 1:
    ...                         note += [""]
    ...                     lo = feature["location"].first() - 1
    ...                     hi = feature["location"].last()
    ...                     rows.append(["Domain", lo, hi, ' '.join(note).strip(), None])
    ...     species = seq.Info.species.split()
    ...     seq_name = "%s.%s" % (species[0][0] + species[1][:3], accession)
    ...     seqs[seq_name] = seq

We convert the sequences to a ``SequenceCollection`` using ``LoadSeqs`` and then save to a fasta formatted file.

.. doctest::
    
    >>> from cogent import LoadSeqs
    >>> seqs = LoadSeqs(data=seqs, aligned=False)
    >>> print seqs.NamedSeqs['Clup.Q28295'].toFasta()
    >Clup.Q28295
    MSPTRLVRVLLALALI...

We convert the features into a PyCogent ``Table`` object, which requires we specify column headings. This can be saved to file if desired, but we don't do that here. For display purposes, we just print the first 10 records.

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> from cogent import LoadTable
    >>> feature_table = LoadTable(header=["Type", "Start", "Stop", "Note",
    ...                    "Disease"], rows=rows)

Printing ``feature_table[:10]`` should result in something like:

.. code-block:: python
    
    ============================================
      Type    Start    Stop      Note    Disease
    --------------------------------------------
    Domain       33     240    VWFD 1           
       SNP      272     273    R -> W       True
    Domain      294     348     TIL 1           
       SNP      317     318    N -> K      False
       SNP      376     377    W -> C       True
    Domain      386     598    VWFD 2           
       SNP      483     484    H -> R      False
       SNP      527     528    N -> S       True
       SNP      549     550    G -> R       True
    Domain      651     707     TIL 2           
    --------------------------------------------
