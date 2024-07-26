Why write apps?
===============

Cogent3 apps solve several problems. Here are just a few.

.. dropdown:: **An analysis that requires multiple apps can be composed into a single process.**

    Defining and using composed apps is simple, and the code is easier to read and check. Compare using a "composed" app, made up of three separate apps

    .. code-block:: python

        # making a composed
        my_steps = step_1 + step_2 + step_3
        step_3_output = my_steps(step_1_input)

    to using each app in order.

    .. code-block:: python

        step_1_output = step_1(step_1_input)
        step_2_output = step_2(step_1_output)
        step_3_output = step_3(step_2_input)

.. dropdown:: **Applying a process to multiple data records can be done without looping.**

    Cogent3 :ref:`data stores <data_stores>` simplify selecting data for analysis and can be applied for batch execution by an app.

    .. code-block:: python

        step_1_inputs = open_data_store("path/to/seqs", suffix="fasta")
        results = list(my_steps.as_completed(step_1_inputs))

.. dropdown:: **Distributing your app as a Cogent3 plugin reduces the learning curve for users.**

    Cogent3 provides mechanisms for discovering your app, getting help on it with an example of how to use it. This mechanism is the same for all apps, which means users do not need to know the structure of your package (see the :ref:`app overview <app_start>`).

.. dropdown:: **The app infrastructure comes with a freebies!**

    When used as part of a composed function, they automatically provide logging. Composed apps can record run failures (everyone encounters bad data ☹️) and make it easier to see how much of an analysis was affected. They also greatly simplify running analyses in parallel. 

.. dropdown:: **Apps are seriously easy to write!**

    Read on to see just how easy they are. In summary, just two lines of code extra plus two extra lines in your projects' `pyproject.toml`, and you have everything you need to distribute your app.

    Use the `cogent3 app cookiecutter template <https://github.com/cogent3/app_template>`_ to get a head start.

Writing your own apps
=====================

Consider the case where you have a directory of files with the same format. You want to apply a consistent procedure to each file separately and produce a corresponding output. If one of the files is incompatible with the calculations, you want to record that and continue processing the rest. This is the use case for which Cogent3 apps are designed. By building your own apps, you can easily incorporate them as a part of a more substantial algorithm. You also get invaluable capabilities such as automated logging and simplifying their execution in parallel.

To write a Cogent3 app, you must make some decisions. What type of processing will your app do? What data type(s) will it accept as input? What data type will it produce as output? Do you need to deal with data that doesn't satisfy a condition?

Defining apps is achieved with the ``define_app``  decorator. Below, we describe the different configuration options and give three examples of writing apps.

.. dropdown:: Specifying the type of app

    The :ref:`supported app types <app_types>` are indicated by the ``AppType`` enum

    .. jupyter-execute::

        from cogent3.app.composable import AppType

        list(AppType)

    The decorator has a default value of ``"generic"``. This means that the app does data transformation and does not, for example, load data from disk or write data to disk (those are the ``loader`` and ``writer`` types).

    If your application is not intended to be applied sequentially, before or after other cogent3 apps, to a series of independent data records of the same type, then you set ``define_app(app_type=AppType.NON_COMPOSABLE)`` (or equivalently ``define_app(app_type="non_composable")``).
    
    .. note:: Non-composable apps cannot be added (or composed) with other apps.

.. dropdown:: Handling not completed values

    Cogent3 apps are designed to handle conditions under which data cannot be completely processed. For example, any errors during execution are stored in a :ref:`not completed <not_completed>` object (``NotCompleted``). All subsequent steps are skipped whenever an app returns a ``NotCompleted`` object.

    These objects can also be used to stop the processing of a particular data record. For example, say you're writing an application that requires sequences to have some minimum length. If an input sequence is shorter than this, then you can create a :ref:`not completed <not_completed>` object and return it. This prevents any future processing, but the reasons for the failure (which you get to specify) can be saved for future reference.

    You may want access to ``NotCompleted`` instances if, for example, you are developing an ``AppType.WRITER``. You can allow your code to "see" them with ``define_app(skip_not_completed=False)``.

.. dropdown:: Supported cogent3 types

    You can use the existing type hints if your function takes or returns ``cogent3`` types. To see these, use the ``defined_types()`` function.

    .. jupyter-execute::

        from cogent3.app.typing import defined_types

        defined_types()

    .. note:: You don't have to use cogent3 types. You can also use standard python types.

Defining a cogent3 app from a function
--------------------------------------

We will write a function that takes a Cogent3 alignment and returns the first *n* positions where the user defines *n*.

.. jupyter-execute::

    from cogent3.app.composable import define_app
    from cogent3.app.typing import AlignedSeqsType

    @define_app
    def n_positions(val: AlignedSeqsType, n=2) -> AlignedSeqsType:
        return val[:n]

The critical elements of a function being defined as an app are:

1. The ``define_app`` decorator is used.
2. Type hints are specified for the function's first argument and its return type.

.. note: Currently, your function can only have one required argument. It can have any number of optional arguments.

.. dropdown:: Using the custom app

    We create an app instance for a specific value of ``n``.

    .. jupyter-execute::

        first4 = n_positions(n=4)
        first4

    The instance's ``repr()`` indicates the wrapped function and the argument values. You use ``first4()`` like all composable apps, e.g.

    .. jupyter-execute::

        from cogent3 import make_aligned_seqs

        aln = make_aligned_seqs(
            data=dict(a="GCAAGCGTTTAT", b="GCTTTTGTCAAT"), moltype="dna"
        )
        result = first4(aln)
        result

Defining a cogent3 app from a class
-----------------------------------

.. jupyter-execute::

    from cogent3.app.composable import define_app
    from cogent3.app.typing import AlignedSeqsType

    @define_app
    class n_positions:
        def __init__(self, n=2):
            self.n = n

        def main(self, val: AlignedSeqsType) -> AlignedSeqsType:
            return val[:self.n]

The critical elements of a class being defined as an app are:

1. The ``define_app`` decorator is used.
2. The class has a ``main()`` method.
3. Type hints are specified for the ``main()`` methods first argument and its return type.

.. dropdown:: Using the custom app

    This is identical to what we did above.

    .. jupyter-execute::

        first4 = n_positions(n=4)

        # we use the alignment defined above

        result = first4(aln)
        result

Custom apps for standard python types
-------------------------------------

In this example, we have two apps that process pure Python types only.

.. jupyter-execute::

    from cogent3.app.composable import define_app

    @define_app
    def lower(arg: str | bytes) -> str | bytes:
        return arg.lower()

    @define_app
    def space(arg: str | bytes) -> str | bytes:
        sep = " " if isinstance(arg, str) else b" "
        arg = arg.split()
        return sep.join(arg)

The only difference to the above examples is we use standard python types for the type hints.

.. dropdown:: Using the two apps

    We create a composed app and apply it.
    
    .. jupyter-execute::

        app = lower() + space()
        app

    .. jupyter-execute::

        app(b"HELLO   there")


App naming conventions
----------------------

Use words in lower case separated by underscores (e.g. ``lower_case``) to name your apps. Apps are callable, just like functions, and the `PEP8 guidelines <https://peps.python.org/pep-0008/#function-and-variable-names>`_ specify this naming style.

If you will make your app available on the Python package index, we recommend prefixing each app with your package name. For example, the `piqtree2 <https://pypi.org/project/piqtree2>`_ library distributes apps with names such as ``piqtree_phylo``.
