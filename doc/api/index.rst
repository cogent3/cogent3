###
API
###

*********************
Utility Functions For
*********************

Loading data from file
======================

These are all top level imports. For example,

.. code-block:: python

    from cogent3 import load_unaligned_seqs

.. toctree::
    :maxdepth: 1

    __init__/cogent3.__init__.load_seq
    __init__/cogent3.__init__.load_aligned_seqs
    __init__/cogent3.__init__.load_unaligned_seqs
    __init__/cogent3.__init__.load_delimited
    __init__/cogent3.__init__.load_table
    __init__/cogent3.__init__.load_tree
    __init__/cogent3.__init__.load_annotations
    __init__/cogent3.__init__.open_data_store
    __init__/cogent3.__init__.open_
    
    

Making cogent3 types from standard Python types
===============================================

These are all top level imports. For example,

.. code-block:: python

    from cogent3 import make_unaligned_seqs

.. toctree::
    :maxdepth: 1

    __init__/cogent3.__init__.make_seq
    __init__/cogent3.__init__.make_aligned_seqs
    __init__/cogent3.__init__.make_unaligned_seqs
    __init__/cogent3.__init__.make_table
    __init__/cogent3.__init__.make_tree

Getting commonly used cogent3 types
===================================

These are all top level imports. For example,

.. code-block:: python

    from cogent3 import get_code

.. toctree::
    :maxdepth: 1

    __init__/cogent3.__init__.get_code
    __init__/cogent3.__init__.get_moltype
    __init__/cogent3.__init__.get_model

Displaying cogent3 builtins
===========================

These are all top level imports. For example,

.. code-block:: python

    from cogent3 import get_code

.. toctree::
    :maxdepth: 1

    __init__/cogent3.__init__.available_codes
    __init__/cogent3.__init__.available_moltypes
    __init__/cogent3.__init__.available_models
    __init__/cogent3.__init__.available_apps

****************************
The Major cogent3 Data Types
****************************

.. toctree::
    :maxdepth: 1

    alignment/alignment
    annotation_db/annotation_db
    sequence/sequence
    genetic_code/genetic_code
    moltype/moltype
    alphabet/alphabet
    table/table
    tree/tree
    cogent3.core.annotation.Feature

************************
Defining Composable Apps
************************

.. toctree::
    :maxdepth: 1

    apps/composable.define_app
    apps/datastore/DataMember
    apps/datastore/DataStoreDirectory
    apps/datastore/ReadOnlyDataStoreZipped
    apps/datastore/DataStoreSqlite
