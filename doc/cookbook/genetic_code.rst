.. _genetic-codes:

Using genetic codes
^^^^^^^^^^^^^^^^^^^

.. note:: These docs now use the ``new_type`` core objects via the following setting.

    .. jupyter-execute::

        import os

        # using new types without requiring an explicit argument
        os.environ["COGENT3_NEW_TYPE"] = "1"

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

From a string

.. jupyter-execute::

    from cogent3 import get_code

    standard_code = get_code(1)
    standard_code.translate("TTTGCAAAC")

This can also be applied to a numpy array.

.. jupyter-execute::

    import numpy
    from cogent3 import get_code

    standard_code = get_code(1)

    standard_code.translate(numpy.array([0, 0, 0, 3, 1, 2, 2, 2, 1], dtype=numpy.uint8))

Conversion to a ``ProteinSequence`` from a ``DnaSequence`` is shown in :ref:`translation`.

Translate all six frames
""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import get_code, make_seq

    standard_code = get_code(1)
    seq = make_seq("ATGCTAACATAAA", moltype="dna")
    translations = standard_code.sixframes(seq)
    print(translations)

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

Getting the alphabet for the genetic code
"""""""""""""""""""""""""""""""""""""""""

The default for the ``get_alphabet()`` method is to return an alphabet representing just the sense codons (a ``SenseCodonAlphabet`` instance).

.. jupyter-execute::

    from cogent3 import get_code

    gc = get_code(1)
    alphabet = gc.get_alphabet()
    len(alphabet)

Setting ``include_stop=True`` returns all codons.

.. jupyter-execute::

    from cogent3 import get_code

    gc = get_code(1)
    alphabet = gc.get_alphabet(include_stop=True)
    type(alphabet)

You can also include "gap state" (i.e. ``"---"``) or "missing state" (``"???"``) codons with the arguments ``include_gap`` and ``include_missing`` respectively.


