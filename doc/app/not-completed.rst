.. jupyter-execute::
    :hide-code:

    import set_working_directory

Tracking records that could not be processed
============================================

.. _not_completed:

The ``NotCompleted`` object
---------------------------

``NotCompleted`` is a special result type that can be produced by a composable app. These objects evaluate to ``False``.

An app can return a ``NotCompleted`` result for one of 2 reasons. The object contains information regarding the input data, where the issue arose and whatever message was provided by the code (like an exception traceback).

``NotCompleted`` FALSE type
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The results when a condition was not met. For example, below I create an app that will return alignments with 2 specific sequences. I'm applying this to a data set where a "Mouse" sequence does not exist. This will produce a FALSE type.

.. jupyter-execute::

    from cogent3 import get_app

    reader = get_app("load_aligned", format="fasta")
    select_seqs = get_app("take_named_seqs", "Mouse", "Human")
    aln = reader("data/primate_brca1.fasta")
    result = select_seqs(aln)
    assert result == False
    result

The ``NotCompleted`` instance has attributes identifying what data failed,

.. jupyter-execute::

    result.source

where the failure occurred

.. jupyter-execute::

    result.origin

and the reason for the failure

.. jupyter-execute::

    result.message


``NotCompleted`` ERROR type
^^^^^^^^^^^^^^^^^^^^^^^^^^^

An ``ERROR`` type is returned if an unexpected condition occurs. This can be an exception raised during the calculation. In our example, we illustrate this by trying to open a file with an incorrect path.

.. jupyter-execute::
    :raises:

    result = reader("primate_brca1.fasta")
    result

Composed functions propagate ``NotCompleted`` results
-----------------------------------------------------

If you have a composed function, with multiple steps and an error occurs then the resulting ``NotCompleted`` result is returned without any of the other steps being applied to it. For example, we make a composed app from both of the above apps.

.. jupyter-execute::

    app = reader + select_seqs
    result = app("data/primate_brca1.fasta")
    result

and

.. jupyter-execute::
    :raises:

    result = app("primate_brca1.fasta")
    result
