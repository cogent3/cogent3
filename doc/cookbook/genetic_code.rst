.. _genetic-codes:

Using genetic codes
^^^^^^^^^^^^^^^^^^^

.. note:: **Alpha Release of the New GeneticCode API**

   We are pleased to announce an alpha release of our new ``GeneticCode`` API. This version can be accessed by specifying the argument ``new_type=True`` in the ``get_code()`` function. 
   
   Please be aware that this alpha release has not been fully integrated with the entire library. Users are encouraged to explore its capabilities but should proceed with caution!

Selecting codes in methods that support them
""""""""""""""""""""""""""""""""""""""""""""

In cases where a ``cogent3`` object method has a ``gc`` argument, you can just use the number under "Code ID" column.

For example, I've created a partial codon in ``"s1"``

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    data = {
        "s1": "GCTCATGCCAGCTCTTTACAGCATGAGAACA--AGT",
        "s2": "ACTCATGCCAACTCATTACAGCATGAGAACAGCAGT",
        "s3": "ACTCATGCCAGCTCATTACAGCATGAGAACAGCAGT",
        "s4": "ACTCATGCCAGCTCATTACAGCATGAGAACAGCAGT",
        "s5": "ACTCATGCCAGCTCAGTACAGCATGAGAACAGCAGT",
    }

    nt_seqs = make_aligned_seqs(data=data, moltype="dna")
    nt_seqs

We specify the genetic code, and we allow incomplete codons. In this case, if a codon contains a gap, they are converted to ``?`` in the translation.

.. jupyter-execute::

    nt_seqs.get_translation(gc=1, incomplete_ok=True)


Translate DNA sequences
"""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import get_code

    standard_code = get_code(1)
    standard_code.translate("TTTGCAAAC")

Conversion to a ``ProteinSequence`` from a ``DnaSequence`` is shown in :ref:`translation`.

Translate all six frames
""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import get_code, make_seq

    standard_code = get_code(1)
    seq = make_seq("ATGCTAACATAAA", moltype="dna")
    translations = standard_code.sixframes(seq)
    print(translations)

Find out how many stops in a frame
""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import get_code, make_seq

    standard_code = get_code(1)
    seq = make_seq("ATGCTAACATAAA", moltype="dna")
    stops_frame1 = standard_code.get_stop_indices(seq, start=0)
    stops_frame1

.. jupyter-execute::

    stop_index = stops_frame1[0]
    seq[stop_index : stop_index + 3]

Translate a codon
"""""""""""""""""

.. jupyter-execute::

    from cogent3 import get_code, make_seq

    standard_code = get_code(1)
    standard_code["TTT"]

or get the codons for a single amino acid

.. jupyter-execute::

    standard_code["A"]

Look up the amino acid corresponding to a single codon
""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import get_code

    standard_code = get_code(1)
    standard_code["TTT"]

Get all the codons for one amino acid
"""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import get_code

    standard_code = get_code(1)
    standard_code["A"]

Get all the codons for a group of amino acids
"""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    targets = ["A", "C"]
    codons = [standard_code[aa] for aa in targets]
    codons

.. jupyter-execute::

    flat_list = sum(codons, [])
    flat_list

Converting the ``CodonAlphabet`` to codon series
""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import get_code

    gc = get_code(1)
    alphabet = gc.get_alphabet()
    print(alphabet)

Obtaining the codons from a ``DnaSequence`` object
""""""""""""""""""""""""""""""""""""""""""""""""""

Use the method ``get_in_motif_size``

.. jupyter-execute::

    from cogent3 import make_seq

    my_seq = make_seq("ATGCACTGGTAA", name="my_gene", moltype="dna")
    codons = my_seq.get_in_motif_size(3)
    codons

Translating a DNA sequence
""""""""""""""""""""""""""

The defaults for ``get_translation()`` include using the standard genetic code and trimming a terminating stop if it exists.

.. jupyter-execute::

    pep = my_seq.get_translation()
    pep

Translating a DNA sequence containing stop codons
"""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::
    :hide-code:

    from cogent3.core.alphabet import AlphabetError

Making a sequence that contains both internal and terminating stop codons.

.. jupyter-execute::
    :raises:

    from cogent3 import make_seq

    seq = make_seq("ATGTGATGGTAA", name="s1", moltype="dna")

Translating this will fail with default settings.

.. jupyter-execute::
    :raises: AlphabetError

    pep = seq.get_translation()

Unless you explicitly allow stop codons

.. jupyter-execute::

    pep = seq.get_translation(include_stop=True)
    pep

