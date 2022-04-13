.. jupyter-execute::
    :hide-code:

    import set_working_directory

********************************************
Tracking records that could not be processed
********************************************

.. _not_completed:

The ``NotCompleted`` object
===========================

This is an object returned when a composable function did not complete a task. ``NotCompleted`` is a special result type that can be produced by a composable app. These objects evaluate to ``False``.

An app can return a ``NotCompleted`` result for one of 2 reasons. The object contains information regarding the input data, where the issue arose and whatever message was provided by the code (like an exception traceback).

``NotCompleted`` FALSE type
===========================

The results when a condition was not met. For example, below I create an app that will return alignments that with 2 specific sequences but I'm including one that does not exist ("Mouse"). So this will fail.

.. jupyter-execute::

    from cogent3.app import io, sample

    reader = io.load_aligned(format="fasta")
    select_seqs = sample.take_named_seqs("Mouse", "Human")
    aln = reader("data/primate_brca1.fasta")
    result = select_seqs(aln)
    result

.. jupyter-execute::

    result == False
    result.type
    result.message

``NotCompleted`` ERROR type
===========================

An ``ERROR`` type is returned if an exception is raised during the calculation. We trigger it in this case by trying to open a non-existent file.

.. jupyter-execute::
    :raises:

    result = reader("primate_brca1.fasta")
    result

Composed functions propagate ``NotCompleted`` results
=====================================================

.. jupyter-execute::

    process = reader + select_seqs
    result = process("data/primate_brca1.fasta")
    result

and

.. jupyter-execute::
    :raises:

    result = process("primate_brca1.fasta")
    result
