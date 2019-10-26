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

.. doctest::

    >>> from cogent3.app import io, sample
    >>> reader = io.load_aligned(format="fasta")
    >>> select_seqs = sample.take_named_seqs("Mouse", "Human")
    >>> aln = reader("data/primate_brca1.fasta")
    >>> result = select_seqs(aln)
    >>> result
    NotCompleted(type=FALSE, origin=take_named_seqs, source="data/primate_brca1.fasta", message="named seq(s) {'Mouse'} not in ['Chimpanzee', 'Galago', 'Gorilla', 'HowlerMon', 'Human', 'Orangutan', 'Rhesus']")


.. doctest::

    >>> result == False
    True
    >>> result.type
    'FALSE'
    >>> result.message
    "named seq(s) {'Mouse'} not in ['Chimpanzee', 'Galago', 'Gorilla', 'HowlerMon', 'Human', 'Orangutan', 'Rhesus']"

``NotCompleted`` ERROR type
===========================

An ``ERROR`` type is returned if an exception is raised during the calculation. We trigger it in this case by trying to open a non-existent file.

.. doctest::

    >>> result = reader("primate_brca1.fasta")
    >>> result
    NotCompleted(type=ERROR, origin=load_aligned, source="primate_brca1.fasta", message="Traceback (most...

Composed functions propagate ``NotCompleted`` results
=====================================================

.. doctest::

    >>> process = reader + select_seqs
    >>> result = process("data/primate_brca1.fasta")
    >>> result
    NotCompleted(type=FALSE, origin=take_named_seqs, source="data/primate_brca1.fasta", message="named seq(s) {'Mouse'} not in ['Chimpanzee', 'Galago', 'Gorilla', 'HowlerMon', 'Human', 'Orangutan', 'Rhesus']")

and

.. doctest::

    >>> result = process("primate_brca1.fasta")
    >>> result
    NotCompleted(type=ERROR, origin=load_aligned, source="primate_brca1.fasta", message="Traceback (most...

If you write results into a ``tinydb`` data store, ``NotCompleted`` objects are saved. After an analysis, you can get a summary of those using methods on the ``WritableTinyDbDataStore`` instance.