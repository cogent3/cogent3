Translate DNA sequences
-----------------------

.. jupyter-execute::

    from cogent3 import get_code

    standard_code = get_code(1)
    standard_code.translate("TTTGCAAAC")

Conversion to a ``ProteinSequence`` from a ``DnaSequence`` is shown in :ref:`translation`.

Translate all six frames
------------------------

.. jupyter-execute::

    from cogent3 import get_code, make_seq

    standard_code = get_code(1)
    seq = make_seq("ATGCTAACATAAA", moltype="dna")
    translations = standard_code.sixframes(seq)
    print(translations)

Find out how many stops in a frame
----------------------------------

.. jupyter-execute::

    from cogent3 import get_code, make_seq

    standard_code = get_code(1)
    seq = make_seq("ATGCTAACATAAA", moltype="dna")
    stops_frame1 = standard_code.get_stop_indices(seq, start=0)
    stops_frame1
    stop_index = stops_frame1[0]
    seq[stop_index : stop_index + 3]

Translate a codon
-----------------

.. jupyter-execute::

    from cogent3 import get_code, make_seq

    standard_code = get_code(1)
    standard_code["TTT"]

or get the codons for a single amino acid

.. jupyter-execute::

    standard_code["A"]

Look up the amino acid corresponding to a single codon
------------------------------------------------------

.. jupyter-execute::

    from cogent3 import get_code

    standard_code = get_code(1)
    standard_code["TTT"]

Get all the codons for one amino acid
-------------------------------------

.. jupyter-execute::

    from cogent3 import get_code

    standard_code = get_code(1)
    standard_code["A"]

Get all the codons for a group of amino acids
---------------------------------------------

.. jupyter-execute::

    targets = ["A", "C"]
    codons = [standard_code[aa] for aa in targets]
    codons
    flat_list = sum(codons, [])
    flat_list

Converting the ``CodonAlphabet`` to codon series
------------------------------------------------

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("AGTACACTGGTT", moltype="dna")
    sorted(my_seq.codon_alphabet())
    len(my_seq.codon_alphabet())

Obtaining the codons from a ``DnaSequence`` object
--------------------------------------------------

Use the method ``get_in_motif_size``

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("ATGCACTGGTAA", name="my_gene", moltype="dna")
    codons = my_seq.get_in_motif_size(3)
    print(codons)

Translating a DNA sequence with a terminating stop codon
--------------------------------------------------------

You can't translate a sequence that contains a stop codon.

.. jupyter-execute::
    :raises: AlphabetError

    pep = my_seq.get_translation()

By removing the trailing stop codon first
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("ATGCACTGGTAA", name="my_gene", moltype="dna")
    seq = my_seq.trim_stop_codon()
    pep = seq.get_translation()
    print(pep.to_fasta())
    print(type(pep))

By slicing the ``DnaSequence`` first
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("CAAATGTATTAA", name="my_gene", moltype="dna")
    pep = my_seq[:-3].get_translation()
    print(pep.to_fasta())
