.. _not_completed:

Tracking records that could not be processed
=============================================

``NotCompleted`` is a special result type provided by the `scinexus <https://pypi.org/project/scinexus/>`_ package. Apps return a ``NotCompleted`` result when a condition is not met (``FALSE`` type) or an unexpected error occurs (``ERROR`` type). These objects evaluate to ``False``, carry information about the failure, and propagate through composed pipelines so that subsequent steps are skipped.

.. todo::

    Link to scinexus NotCompleted documentation once available.

..
    Placeholder: <scinexus NotCompleted docs URL>

A cogent3 example
-----------------

.. jupyter-execute::
    :hide-code:

    import set_working_directory

Below, an app that selects specific sequences returns a ``NotCompleted`` because the requested sequence ``"Mouse"`` does not exist in the primate data set.

.. jupyter-execute::

    from cogent3 import get_app

    reader = get_app("load_aligned", format_name="fasta")
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

Composed functions propagate ``NotCompleted`` results
-----------------------------------------------------

If you have a composed function and an error occurs, the ``NotCompleted`` result is returned without any subsequent steps being applied.

.. jupyter-execute::

    app = reader + select_seqs
    result = app("data/primate_brca1.fasta")
    result
