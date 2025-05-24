.. jupyter-execute::
    :hide-code:

    import set_working_directory

***********************
Available genetic codes
***********************

.. note:: These docs now use the ``new_type`` core objects via the following setting.

    .. jupyter-execute::

        import os

        # using new types without requiring an explicit argument
        os.environ["COGENT3_NEW_TYPE"] = "1"

.. jupyter-execute::

    from cogent3 import available_codes

    available_codes()


Getting a genetic code with ``get_code()``
==========================================

This function can be used directly to get a genetic code. We will get the code with ID 4.

.. jupyter-execute::

    from cogent3 import get_code

    gc = get_code(4)
    gc

