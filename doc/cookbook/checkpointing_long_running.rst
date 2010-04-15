Checkpointing optimisation runs
===============================

.. sectionauthor Gavin Huttley

A common problem in HPC systems is to make sure a long running process is capable of restarting after interruptions by restoring to the last check pointed state. The optimiser class code has this capability, for example and we'll illustrate that here. We first construct a likelihood function object.

.. doctest::
    
    >>> from cogent import LoadSeqs, LoadTree
    >>> from cogent.evolve.models import F81
    >>> aln = LoadSeqs('data/primate_brca1.fasta')
    >>> tree = LoadTree('data/primate_brca1.tree')
    >>> sub_model = F81()
    >>> lf = sub_model.makeLikelihoodFunction(tree)
    >>> lf.setAlignment(aln)

We then start an optimisation, providing a filename for checkpointing and specifying a time-interval in (which we make very short here to ensure something get's written, for longer running functions the default ``interval`` setting is fine). Calling ``optimise`` then results in the notice that checkpoint's are being written.

.. doctest::
    
    >>> checkpoint_fn = 'checkpoint_this.txt'
    >>> lf.optimise(filename=checkpoint_fn, interval=100, show_progress = False)
    CHECKPOINTING to file 'checkpoint_this.txt'...

Recovering from a real run that was interrupted generates an additional notification: ``RESUMING from file ..``. For the purpose of this snippet we just show that the checkpoint file exists.

.. doctest::
    
    >>> import cPickle
    >>> data = cPickle.load(open(checkpoint_fn))
    >>> print data
    <cogent.maths.simannealingoptimiser.AnnealingRun object...

.. following cleans up files

.. doctest::
    :hide:
    
    >>> from cogent.util.misc import remove_files
    >>> remove_files([checkpoint_fn], error_on_missing=False)

