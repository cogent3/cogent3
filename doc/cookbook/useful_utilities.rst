****************
Useful Utilities
****************

.. authors, Daniel McDonald, Gavin Huttley, Antonio Gonzalez Pena, Rob Knight

.. include:: union_dict.rst

Using Cogent3's optimisers for your own functions
===================================================

You have a function that you want to maximise/minimise. The parameters in your function may be bounded (must lie in a specific interval) or not. The ``cogent3`` optimisers can be applied to these cases. The ``Powell`` (a local optimiser) and ``SimulatedAnnealing`` (a global optimiser) classes in particular have had their interfaces standardised for such use cases. We demonstrate for a very simple function below.

We write a simple factory function that uses a provided value for omega to compute the squared deviation from an estimate, then use it to create our optimisable function.

.. doctest::

    >>> import numpy
    >>> def DiffOmega(omega):
    ...     def omega_from_S(S):
    ...         omega_est = S/(1-numpy.e**(-1*S))
    ...         return abs(omega-omega_est)**2
    ...     return omega_from_S
    >>> omega = 0.1
    >>> f = DiffOmega(omega)

We then import the minimise function and use it to minimise the function, obtaining the fit statistic and the associated estimate of S. Note that we provide lower and upper bounds (which are optional) and an initial guess for our parameter of interest (``S``).

.. doctest::

    >>> from cogent3.maths.optimisers import minimise, maximise
    >>> S = minimise(f,  # the function
    ...     xinit=1.0, # the initial value
    ...     bounds=(-100, 100), # [lower,upper] bounds for the parameter
    ...     local=True, # just local optimisation, not Simulated Annealing
    ...     show_progress=False)
    >>> assert 0.0 <= f(S) < 1e-6
    >>> print('S=%.4f' % S)
    S=-3.6150

The minimise and maximise functions can also handle multidimensional optimisations, just make xinit (and the bounds) lists rather than scalar values.

Miscellaneous functions
=======================

.. index:: cogent3.util.misc

Identity testing
----------------

Basic ``identity`` function to avoid having to test explicitly for None

.. doctest::

    >>> from cogent3.util.misc import identity
    >>> my_var = None
    >>> if identity(my_var):
    ...   print("foo")
    ... else:
    ...   print("bar")
    ...
    bar

Force a variable to be iterable
-------------------------------

This support method will force a variable to be an iterable, allowing you to guarantee that the variable will be safe for use in, say, a ``for`` loop.

.. doctest::

    >>> from cogent3.util.misc import iterable
    >>> my_var = 10
    >>> for i in my_var:
    ...   print("will not work")
    ...
    Traceback (most recent call last):
    TypeError: 'int' object is not iterable
    >>> for i in iterable(my_var):
    ...   print(i)
    ...
    10

Curry a function
----------------

curry(f,x)(y) = f(x,y) or = lambda y: f(x,y). This was modified from the Python Cookbook. Docstrings are also carried over.

.. doctest::

    >>> from cogent3.util.misc import curry
    >>> def foo(x,y):
    ...   """Some function"""
    ...   return x + y
    ...
    >>> bar = curry(foo, 5)
    >>> print(bar.__doc__)
     curry(foo,5)
    == curried from foo ==
     Some function
    >>> bar(10)
    15

Test to see if an object is iterable
------------------------------------

Perform a simple test to see if an object supports iteration

.. doctest::

    >>> from cogent3.util.misc import is_iterable
    >>> can_iter = [1,2,3,4]
    >>> cannot_iter = 1.234
    >>> is_iterable(can_iter)
    True
    >>> is_iterable(cannot_iter)
    False

Test to see if an object is a single char
-----------------------------------------

Perform a simple test to see if an object is a single character

.. doctest::

    >>> from cogent3.util.misc import is_char
    >>> class foo:
    ...   pass
    ...
    >>> is_char('a')
    True
    >>> is_char('ab')
    False
    >>> is_char(foo())
    False

Flatten a deeply nested iterable
--------------------------------

To flatten a deeply nested iterable, use ``recursive_flatten``. This method supports multiple levels of nesting, and multiple iterable types

.. doctest::

    >>> from cogent3.util.misc import recursive_flatten
    >>> l = [[[[1,2], 'abcde'], [5,6]], [7,8], [9,10]]
    >>> recursive_flatten(l)
    [1, 2, 'a', 'b', 'c', 'd', 'e', 5, 6, 7, 8, 9, 10]

Test to determine if ``list`` of ``tuple``
------------------------------------------

Perform a simple check to see if an object is not a list or a tuple

.. doctest::

    >>> from cogent3.util.misc import not_list_tuple
    >>> not_list_tuple(1)
    True
    >>> not_list_tuple([1])
    False
    >>> not_list_tuple('ab')
    True

Create a case-insensitive iterable
----------------------------------

Create a case-insensitive object, for instance, if you want the key 'a' and 'A' to point to the same item in a dict

.. doctest::

    >>> from cogent3.util.misc import add_lowercase
    >>> d = {'A':5,'B':6,'C':7,'foo':8,42:'life'}
    >>> add_lowercase(d)  # doctest: +SKIP
    {'A': 5, 'a': 5, 'C': 7, 'B': 6, 42: 'life', 'c': 7, 'b': 6, 'foo': 8}

Construct a distance matrix lookup function
-------------------------------------------

Automatically construct a distance matrix lookup function. This is useful for maintaining flexibility about whether a function is being computed or if a lookup is being used

.. doctest::

    >>> from cogent3.util.misc import DistanceFromMatrix
    >>> from numpy import array
    >>> m = array([[1,2,3],[4,5,6],[7,8,9]])
    >>> f = DistanceFromMatrix(m)
    >>> f(0,0)
    1
    >>> f(1,2)
    6

Check class types
-----------------

Check an object against base classes or derived classes to see if it is acceptable

.. doctest::

    >>> from cogent3.util.misc import ClassChecker
    >>> class not_okay(object):
    ...   pass
    ...
    >>> no = not_okay()
    >>> class okay(object):
    ...   pass
    ...
    >>> o = okay()
    >>> class my_dict(dict):
    ...   pass
    ...
    >>> md = my_dict()
    >>> cc = ClassChecker(str, okay, dict)
    >>> o in cc
    True
    >>> no in cc
    False
    >>> 5 in cc
    False
    >>> {'a':5} in cc
    True
    >>> 'asasas' in cc
    True
    >>> md in cc
    True

Delegate to a separate object
-----------------------------

Delegate object method calls, properties and variables to the appropriate object. Useful to combine multiple objects together while assuring that the calls will go to the correct object.

.. doctest::

    >>> from cogent3.util.misc import Delegator
    >>> class ListAndString(list, Delegator):
    ...   def __init__(self, items, string):
    ...     Delegator.__init__(self, string)
    ...     for i in items:
    ...       self.append(i)
    ...
    >>> ls = ListAndString([1,2,3], 'ab_cd')
    >>> len(ls)
    3
    >>> ls[0]
    1
    >>> ls.upper()
    'AB_CD'
    >>> ls.split('_')
    ['ab', 'cd']

Wrap a function to hide from a class
------------------------------------

Wrap a function to hide it from a class so that it isn't a method.

.. doctest::

    >>> from cogent3.util.misc import FunctionWrapper
    >>> f = FunctionWrapper(str)
    >>> f
    <cogent3.util.misc.FunctionWrapper object at ...
    >>> f(123)
    '123'

Construct a constrained container
---------------------------------

Wrap a container with a constraint. This is useful for enforcing that the data contained is valid within a defined context. Cogent3 provides a base ``ConstrainedContainer`` which can be used to construct user-defined constrained objects. Cogent3 also provides ``ConstrainedString``, ``ConstrainedList``, and ``ConstrainedDict``. These provided types fully cover the builtin types while staying integrated with the ``ConstrainedContainer``.

Here is a light example of the ``ConstrainedDict``

.. doctest::

    >>> from cogent3.util.misc import ConstrainedDict
    >>> d = ConstrainedDict({'a':1,'b':2,'c':3}, constraint='abc')
    >>> d  # doctest: +SKIP
    {'a': 1, 'c': 3, 'b': 2}
    >>> d['d'] = 5
    Traceback (most recent call last):
    cogent3.util.misc.ConstraintError: Item 'd' not in constraint 'abc'

