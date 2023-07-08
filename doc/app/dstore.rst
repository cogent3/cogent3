.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _data_stores:

Data stores -- collections of data records
==========================================

If you download :download:`raw.zip <data/raw.zip>` and unzip it, you will see it contains 1,035 files ending with a ``.fa`` filename suffix. (It also contains a tab delimited file and a log file, which we ignore for now.) The directory ``raw`` is a "data store" and the ``.fa`` files are "members" of it. In summary, a :index:`data store` is a collection of members of the same "type". This means we can apply the same application to every member.

Types of data store
-------------------

.. csv-table::
    :header: "Class Name", "Supported Operations", "Supported Data Types", "Identifying Suffix"

    ``DataStoreDirectory``, "read / write / append", text, None
    ``ReadOnlyDataStoreZipped``, read, text, ``.zip``
    ``DataStoreSqlite``, "read, write, append", "text or bytes", ``.sqlitedb``

.. note:: The ``ReadOnlyDataStoreZipped`` is just a compressed ``DataStoreDirectory``.

The structure of data stores
----------------------------

If a directory was not created by ``cogent3`` as a ``DataStoreDirectory`` then it has only the structure that existed previously.

If a data store was created by ``cogent3``, either as a directory or as a ``sqlitedb``, then it contains four types of data: completed records, *not* completed records, log files and md5 files. In a ``DataStoreDirectory``, these are organised using the file system. The :index:`completed members` are valid data records (as distinct from :ref:`not completed <not_completed>`) and are at the top level. The remaining types are in subdirectories.

::

    demo_dstore
    ├── logs
    ├── md5
    ├── not_completed
    └── ... <the completed members>

``logs/`` stores scitrack_ log files produced by ``cogent3.app`` writer apps. ``md5/`` stores plain text files with the md5 sum of a corresponding data member which are used to check the integrity of the data store.

The ``DataStoreSqlite`` stores the same information, just in SQL tables.

.. _scitrack: https://github.com/HuttleyLab/scitrack

Supported operations on a data store
------------------------------------

All data store classes can be iterated over, indexed, checked for membership. These operations return a ``DataMember`` object. In addition to providing access to members, the data store classes have convenience methods for describing their contents and providing summaries of log files that are included and of the ``NotCompleted`` members (see :ref:`not_completed`).

Opening a data store
--------------------

Use the :index:`open_data_store()` function, illustrated below. Use the mode argument to identify whether to open as read only (``mode="r"``), write (``mode=w``) or append(``mode="a"``).

Opening a read only data store
------------------------------

We open the zipped directory described above, defining the filenames ending in ``.fa`` as the data store members. All files within the directory become members of the data store (unless we use the ``limit`` argument).

.. jupyter-execute::

    from cogent3 import open_data_store

    dstore = open_data_store("data/raw.zip", suffix="fa", mode="r")
    print(dstore)

Summarising the data store
--------------------------

The ``.describe`` property demonstrates that there are only completed members.

.. jupyter-execute::

    dstore.describe

Looping over a data store
-------------------------

.. jupyter-execute::

    for m in dstore[:5]:
        print(m)

.. _data_member:

Data store “members”
--------------------

These are able to read their own raw data.

.. jupyter-execute::

    m = dstore[0]
    m

.. jupyter-execute::

    m.read()[:20]  # truncating

.. note:: For a ``DataStoreSqlite`` member, the default data storage format is as bytes. So reading the content of an individual record is best done using the ``load_db`` app.

Making a writeable data store
-----------------------------

The creation of a writeable data store is specified with ``mode="w"``, or (to append) ``mode="a"``. In the former case, any existing records are overwritten. In the latter case, existing records are ignored.

``DataStoreSqlite`` stores serialised data
------------------------------------------

When you specify a Sqlitedb data store as your output (by using ``open_data_store()``) you write multiple records into a single file making distribution easier.

One important issue to note is the process which creates a Sqlitedb “locks” the file. If that process exits unnaturally (e.g. the run that was producing it was interrupted) then the file may remain in a locked state. If the db is in this state, ``cogent3`` will not modify it unless you explicitly unlock it.

This is represented in the display as shown below.

.. jupyter-execute::

    dstore = open_data_store("data/demo-locked.sqlitedb")
    dstore.describe

To unlock, you execute the following:

.. jupyter-execute::

    dstore.unlock(force=True)

Interrogating run logs
----------------------

If you use the ``apply_to()`` method, a scitrack_ logfile will be stored in the data store. This includes useful information regarding the run conditions that produced the contents of the data store.

.. jupyter-execute::

    dstore.summary_logs

Log files can be accessed vial a special attribute.

.. jupyter-execute::

    dstore.logs

Each element in that list is a ``DataMember`` which you can use to get the data contents.

.. jupyter-execute::

    print(dstore.logs[0].read()[:225])  # truncated for clarity

Pulling it all together
=======================

We will translate the DNA sequences in ``raw.zip`` into amino acid and store them as sqlite database. We will interrogate the generated data store to gtet a synopsis of the results.

Defining the data stores for analysis
-------------------------------------

Loading our input data

.. jupyter-execute::

    from cogent3 import open_data_store

    in_dstore = open_data_store("data/raw.zip", suffix="fa")

Creating our output ``DataStoreSqlite``

.. jupyter-execute::

    out_dstore = open_data_store("translated.sqlitedb", mode="w")

Create an app and apply it
--------------------------

We need apps to load the data, translate it and then to write the translated sequences out. We define those and compose into a single app.

.. jupyter-execute::

    from cogent3 import get_app

    load = get_app("load_unaligned", moltype="dna")
    translate = get_app("translate_seqs")
    write = get_app("write_db", data_store=out_dstore)
    app = load + translate + write
    app

We apply the app to all members of ``in_dstore``. The results will be written to ``out_dstore``.

.. jupyter-execute::

    out_dstore = app.apply_to(in_dstore)

Inspecting the outcome
----------------------

The ``.describe`` method gives us an analysis level summary.

.. jupyter-execute::

    out_dstore.describe

We confirm the data store integrity

.. jupyter-execute::

    out_dstore.validate()

We can examine why some input data could not be processed by looking at the summary of the not completed records.

.. jupyter-execute::

    out_dstore.summary_not_completed

We see they all came from the ``translate_seqs`` step. Some had a terminal stop codon while others had a length that was not divisible by 3.

.. note::
    
    The ``.completed`` and ``.not_completed`` attributes give access to the different types of members while the ``.members`` attribute gives them all. For example,

    .. jupyter-execute::

        len(out_dstore.not_completed)

    is the same as in the ``describe`` output and each element is a ``DataMember``.

    .. jupyter-execute::

        out_dstore.not_completed[:2]

.. jupyter-execute::
    :hide-code:

    import pathlib

    fn = pathlib.Path("translated.sqlitedb")
    fn.unlink()
