.. jupyter-execute::
    :hide-code:

    import set_working_directory

Specifying data for analysis
============================

We introduce the concept of a “data store”. This represents the data record(s) that you want to analyse. It can be a single file, a directory of files, a zipped directory of files or a single ``tinydb`` file containing multiple data records.

We represent this concept by a ``DataStore`` class. There are different flavours of these:

-  directory based
-  zip archive based
-  TinyDB based (this is a NoSQL json based data base)

These can be read only or writable. All of these types support being indexed, iterated over, filtered, etc.. The ``tinydb`` variants do have some unique abilities (discussed below).

A read only data store
----------------------

To create one of these, you provide a ``path`` AND a ``suffix`` of the files within the directory / zip that you will be analysing. (If the path ends with ``.tinydb``, no file suffix is required.)

.. jupyter-execute::

    from cogent3.app.io import get_data_store

    dstore = get_data_store("data/raw.zip", suffix="fa*", limit=5)
    dstore

Data store “members”
--------------------

These are able to read their own raw data.

.. jupyter-execute::

    m = dstore[0]
    m

.. jupyter-execute::

    m.read()[:20]  # truncating

Showing the last few members
----------------------------

Use the ``head()`` method to see the first few.

.. jupyter-execute::

    dstore.tail()

Filtering a data store for specific members
-------------------------------------------

.. jupyter-execute::

    dstore.filtered("*ENSG00000067704*")

Looping over a data store
-------------------------

.. jupyter-execute::

    for m in dstore:
        print(m)

Making a writeable data store
-----------------------------

The creation of a writeable data store is handled for you by the different writers we provide under ``cogent3.app.io``.

.. warning:: The ``WritableZippedDataStore`` is deprecated.

TinyDB data stores are special
------------------------------

When you specify a TinyDB data store as your output (by using ``io.write_db()``), you get additional features that are useful for dissecting the results of an analysis.

One important issue to note is the process which creates a TinyDB “locks” the file. If that process exits unnaturally (e.g. the run that was producing it was interrupted) then the file may remain in a locked state. If the db is in this state, ``cogent3`` will not modify it unless you explicitly unlock it.

This is represented in the display as shown below.

.. jupyter-execute::

    dstore = get_data_store("data/demo-locked.tinydb")
    dstore.describe

To unlock, you execute the following:

.. jupyter-execute::

    dstore.unlock(force=True)

Interrogating run logs
~~~~~~~~~~~~~~~~~~~~~~

If you use the ``apply_to(logger=true)`` method, a ``scitrack`` logfile will be included in the data store. This includes useful information regarding the run conditions that produced the contents of the data store.

.. jupyter-execute::

    dstore.summary_logs

Log files can be accessed vial a special attribute.

.. jupyter-execute::

    dstore.logs

Each element in that list is a ``DataStoreMember`` which you can use to get the data contents.

.. jupyter-execute::

    print(dstore.logs[0].read()[:225])  # truncated for clarity
