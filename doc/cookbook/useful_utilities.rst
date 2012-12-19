****************
Useful Utilities
****************

.. authors, Daniel McDonald, Gavin Huttley, Antonio Gonzalez Pena, Rob Knight

Using PyCogent's optimisers for your own functions
==================================================

You have a function that you want to maximise/minimise. The parameters in your function may be bounded (must lie in a specific interval) or not. The cogent optimisers can be applied to these cases. The ``Powell`` (a local optimiser) and ``SimulatedAnnealing`` (a global optimiser) classes in particular have had their interfaces standardised for such use cases. We demonstrate for a very simple function below.

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
    
    >>> from cogent.maths.optimisers import minimise, maximise
    >>> S = minimise(f,  # the function
    ...     xinit=1.0, # the initial value
    ...     bounds=(-100, 100), # [lower,upper] bounds for the parameter
    ...     local=True) # just local optimisation, not Simulated Annealing
    >>> assert 0.0 <= f(S) < 1e-6
    >>> print 'S=%.4f' % S
    S=-3.6150

The minimise and maximise functions can also handle multidimensional optimisations, just make xinit (and the bounds) lists rather than scalar values.

Fitting a function to a giving set of x and y values
====================================================

Giving a set of values for ``x`` and ``y`` fit a function ``func`` that has ``n_params`` using simplex iterations to minimize the error between the model ``func`` to fit and the given values. Here we fitting an exponential function.

.. doctest::
    :hide:
    
    >>> from numpy.random import seed
    >>> seed(42) # so the results are not volatile


.. doctest::

    >>> from numpy import array, arange, exp
    >>> from numpy.random import rand, seed
    >>> from cogent.maths.fit_function import fit_function
    >>> # creating x values
    >>> x = arange(-1,1,.01)
    >>>
    >>> # defining our fitting function
    >>> def f(x,a):
    ...     return exp(a[0]+x*a[1])
    ...
    >>> # getting our real y
    >>> y = f(x,a=[2,5])
    >>>
    >>> # creating our noisy y
    >>> y_noise = y + rand(len(y))*5
    >>>
    >>> # fitting our noisy data to the function using 1 iteration
    >>> params = fit_function(x, y_noise, f, 2, 1)
    >>> params
    array([ 2.0399908 ,  4.96109191])
    >>>
    >>> # fitting our noisy data to the function using 1 iteration
    >>> params = fit_function(x, y_noise, f, 2, 5)
    >>> params
    array([ 2.0399641 ,  4.96112469])

Cartesian products
==================

*To be written.*

.. cogent.util.transform

Miscellaneous functions
=======================

.. index:: cogent.util.misc

Identity testing
^^^^^^^^^^^^^^^^

Basic ``identity`` function to avoid having to test explicitly for None

.. doctest::

    >>> from cogent.util.misc import identity
    >>> my_var = None
    >>> if identity(my_var):
    ...   print "foo"
    ... else:
    ...   print "bar"
    ... 
    bar

One-line if/else statement
^^^^^^^^^^^^^^^^^^^^^^^^^^

Convenience function for performing one-line if/else statements. This is similar to the C-style ternary operator:

.. doctest::

    >>> from cogent.util.misc import if_
    >>> result = if_(4 > 5, "Expression is True", "Expression is False")
    >>> result
    'Expression is False'

However, the value returned is evaluated, but not called. For instance:

.. doctest::

    >>> from cogent.util.misc import if_
    >>> def foo():
    ...   print "in foo"
    ... 
    >>> def bar():
    ...   print "in bar"
    ...
    >>> if_(4 > 5, foo, bar)
    <function bar at...

Force a variable to be iterable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This support method will force a variable to be an iterable, allowing you to guarantee that the variable will be safe for use in, say, a ``for`` loop.

.. doctest::

    >>> from cogent.util.misc import iterable
    >>> my_var = 10
    >>> for i in my_var:
    ...   print "will not work"
    ... 
    Traceback (most recent call last):
    TypeError: 'int' object is not iterable
    >>> for i in iterable(my_var):
    ...   print i
    ... 
    10

Obtain the index of the largest item
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To determine the index of the largest item in any iterable container, use ``max_index``:

.. doctest::

    >>> from cogent.util.misc import max_index
    >>> l = [5,4,2,2,6,8,0,10,0,5]
    >>> max_index(l)
    7

.. note:: Will return the lowest index of duplicate max values

Obtain the index of the smallest item
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To determine the index of the smallest item in any iterable container, use ``min_index``:

.. doctest::

    >>> from cogent.util.misc import min_index
    >>> l = [5,4,2,2,6,8,0,10,0,5]
    >>> min_index(l)
    6

.. note:: Will return the lowest index of duplicate min values

Remove a nesting level
^^^^^^^^^^^^^^^^^^^^^^

To flatten a 2-dimensional list, you can use ``flatten``:

.. doctest::

    >>> from cogent.util.misc import flatten
    >>> l = ['abcd','efgh','ijkl']
    >>> flatten(l)
    ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']

Convert a nested tuple into a list
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Conversion of a nested ``tuple`` into a ``list`` can be performed using ``deep_list``:

.. doctest::

    >>> from cogent.util.misc import deep_list
    >>> t = ((1,2),(3,4),(5,6))
    >>> deep_list(t)
    [[1, 2], [3, 4], [5, 6]]

Simply calling ``list`` will not convert the nested items:

.. doctest::

    >>> list(t)
    [(1, 2), (3, 4), (5, 6)]

Convert a nested list into a tuple
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Conversion of a nested ``list`` into a ``tuple`` can be performed using ``deep_list``:

.. doctest::

    >>> from cogent.util.misc import deep_tuple
    >>> l = [[1,2],[3,4],[5,6]]
    >>> deep_tuple(l)
    ((1, 2), (3, 4), (5, 6))

Simply calling ``tuple`` will not convert the nested items:

.. doctest::

    >>> tuple(l)
    ([1, 2], [3, 4], [5, 6])

Testing if an item is between two values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Same as: min <= number <= max, although it is quickly readable within code

.. doctest::

    >>> from cogent.util.misc import between
    >>> between((3,5),4)
    True
    >>> between((3,5),6)
    False

Return combinations of items
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Combinate`` returns all k-combinations of items. For instance:

.. doctest::

    >>> from cogent.util.misc import combinate
    >>> list(combinate([1,2,3],0))
    [[]]
    >>> list(combinate([1,2,3],1))
    [[1], [2], [3]]
    >>> list(combinate([1,2,3],2))
    [[1, 2], [1, 3], [2, 3]]
    >>> list(combinate([1,2,3],3))
    [[1, 2, 3]]

Save and load gzip'd files
^^^^^^^^^^^^^^^^^^^^^^^^^^

These handy methods will ``cPickle`` an object and automagically gzip the file. You can also then reload the object at a later date.

.. doctest::

    >>> from cogent.util.misc import gzip_dump, gzip_load
    >>> class foo(object):
    ...   some_var = 5
    ... 
    >>> bar = foo()
    >>> bar.some_var = 10
    >>> # gzip_dump(bar, 'test_file')
    >>> # new_bar = gzip_load('test_file')
    >>> # isinstance(new_bar, foo)

.. note:: The above code does work, but cPickle won't write out within doctest

Curry a function
^^^^^^^^^^^^^^^^

curry(f,x)(y) = f(x,y) or = lambda y: f(x,y). This was modified from the Python Cookbook. Docstrings are also carried over.

.. doctest::

    >>> from cogent.util.misc import curry
    >>> def foo(x,y):
    ...   """Some function"""
    ...   return x + y
    ... 
    >>> bar = curry(foo, 5)
    >>> print bar.__doc__
     curry(foo,5)
    == curried from foo ==
     Some function
    >>> bar(10)
    15

Test to see if an object is iterable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Perform a simple test to see if an object supports iteration

.. doctest::

    >>> from cogent.util.misc import is_iterable
    >>> can_iter = [1,2,3,4]
    >>> cannot_iter = 1.234
    >>> is_iterable(can_iter)
    True
    >>> is_iterable(cannot_iter)
    False

Test to see if an object is a single char
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Perform a simple test to see if an object is a single character

.. doctest::

    >>> from cogent.util.misc import is_char
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To flatten a deeply nested iterable, use ``recursive_flatten``. This method supports multiple levels of nesting, and multiple iterable types

.. doctest::

    >>> from cogent.util.misc import recursive_flatten
    >>> l = [[[[1,2], 'abcde'], [5,6]], [7,8], [9,10]]
    >>> recursive_flatten(l)
    [1, 2, 'a', 'b', 'c', 'd', 'e', 5, 6, 7, 8, 9, 10]

Test to determine if ``list`` of ``tuple``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Perform a simple check to see if an object is not a list or a tuple

.. doctest::

    >>> from cogent.util.misc import not_list_tuple
    >>> not_list_tuple(1)
    True
    >>> not_list_tuple([1])
    False
    >>> not_list_tuple('ab')
    True

Unflatten items to row-width
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unflatten an iterable of items to a specified row-width. This does reverse the effect of ``zip`` as the lists produced are not interleaved.

.. doctest::

    >>> from cogent.util.misc import unflatten
    >>> l = [1,2,3,4,5,6,7,8]
    >>> unflatten(l,1)
    [[1], [2], [3], [4], [5], [6], [7], [8]]
    >>> unflatten(l,2)
    [[1, 2], [3, 4], [5, 6], [7, 8]]
    >>> unflatten(l,3)
    [[1, 2, 3], [4, 5, 6]]
    >>> unflatten(l,4)
    [[1, 2, 3, 4], [5, 6, 7, 8]]

Unzip items
^^^^^^^^^^^

Reverse the effects of a ``zip`` method, i.e. produces separate lists from tuples

.. doctest::

    >>> from cogent.util.misc import unzip
    >>> l = ((1,2),(3,4),(5,6))
    >>> unzip(l)
    [[1, 3, 5], [2, 4, 6]]

Select items in order
^^^^^^^^^^^^^^^^^^^^^

Select items in a specified order

.. doctest::

    >>> from cogent.util.misc import select
    >>> select('ea', {'a':1,'b':5,'c':2,'d':4,'e':6})
    [6, 1]
    >>> select([0,4,8], 'abcdefghijklm')
    ['a', 'e', 'i']

Obtain the index sort order
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Obtain the indices for items in sort order. This is similar to numpy.argsort, but will work on any iterable that implements the necessary ``cmp`` methods

.. doctest::

    >>> from cogent.util.misc import sort_order
    >>> sort_order([4,2,3,5,7,8])
    [1, 2, 0, 3, 4, 5]
    >>> sort_order('dcba')
    [3, 2, 1, 0]

Find overlapping pattern occurrences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Find all of the overlapping occurrences of a pattern within a text

.. doctest::

    >>> from cogent.util.misc import find_all
    >>> text = 'aaaaaaa'
    >>> pattern = 'aa'
    >>> find_all(text, pattern)
    [0, 1, 2, 3, 4, 5]
    >>> text = 'abababab'
    >>> pattern = 'aba'
    >>> find_all(text, pattern)
    [0, 2, 4]

Find multiple pattern occurrences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Find all of the overlapping occurrences of multiple patterns within a text. Returned indices are sorted, each index is the start position of one of the patterns

.. doctest::

    >>> from cogent.util.misc import find_many
    >>> text = 'abababcabab'
    >>> patterns = ['ab','abc']
    >>> find_many(text, patterns)
    [0, 2, 4, 4, 7, 9]

Safely remove a trailing underscore
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

'Unreserve' a mutation of Python reserved words

.. doctest::

    >>> from cogent.util.misc import unreserve
    >>> unreserve('class_')
    'class'
    >>> unreserve('class')
    'class'

Create a case-insensitive iterable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a case-insensitive object, for instance, if you want the key 'a' and 'A' to point to the same item in a dict

.. doctest::

    >>> from cogent.util.misc import add_lowercase
    >>> d = {'A':5,'B':6,'C':7,'foo':8,42:'life'}
    >>> add_lowercase(d)
    {'A': 5, 'a': 5, 'C': 7, 'B': 6, 42: 'life', 'c': 7, 'b': 6, 'foo': 8}

Extract data delimited by differing left and right delimiters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Extract data from a line that is surrounded by different right/left delimiters

.. doctest::

    >>> from cogent.util.misc import extract_delimited
    >>> line = "abc[def]ghi"
    >>> extract_delimited(line,'[',']')
    'def'

Invert a dictionary
^^^^^^^^^^^^^^^^^^^

Get a dictionary with the values set as keys and the keys set as values

.. doctest::

    >>> from cogent.util.misc import InverseDict
    >>> d = {'some_key':1,'some_key_2':2}
    >>> InverseDict(d)
    {1: 'some_key', 2: 'some_key_2'}

.. note:: An arbitrary key will be set if there are multiple keys with the same value

Invert a dictionary with multiple keys having the same value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get a dictionary with the values set as keys and the keys set as values. Can handle the case where multiple keys point to the same values

.. doctest::

    >>> from cogent.util.misc import InverseDictMulti
    >>> d = {'some_key':1,'some_key_2':1}
    >>> InverseDictMulti(d)
    {1: ['some_key_2', 'some_key']}
    >>> 

Get mapping from sequence item to all positions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``DictFromPos`` returns the positions of all items seen within a sequence. This is useful for obtaining, for instance, nucleotide counts and positions

.. doctest::

    >>> from cogent.util.misc import DictFromPos
    >>> seq = 'aattggttggaaggccgccgttagacg'
    >>> DictFromPos(seq)
    {'a': [0, 1, 10, 11, 22, 24], 'c': [14, 15, 17, 18, 25], 't': [2, 3, 6, 7, 20, 21], 'g': [4, 5, 8, 9, 12, 13, 16, 19, 23, 26]}

Get the first index of occurrence for each item in a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``DictFromFirst`` will return the first location of each item in a sequence

.. doctest::
    
    >>> from cogent.util.misc import DictFromFirst
    >>> seq = 'aattggttggaaggccgccgttagacg'
    >>> DictFromFirst(seq)
    {'a': 0, 'c': 14, 't': 2, 'g': 4}

Get the last index of occurrence for each item in a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``DictFromLast`` will return the last location of each item in a sequence

.. doctest::

    >>> from cogent.util.misc import DictFromLast
    >>> seq = 'aattggttggaaggccgccgttagacg'
    >>> DictFromLast(seq)
    {'a': 24, 'c': 25, 't': 21, 'g': 26}

Construct a distance matrix lookup function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Automatically construct a distance matrix lookup function. This is useful for maintaining flexibility about whether a function is being computed or if a lookup is being used

.. doctest::

    >>> from cogent.util.misc import DistanceFromMatrix
    >>> from numpy import array
    >>> m = array([[1,2,3],[4,5,6],[7,8,9]])
    >>> f = DistanceFromMatrix(m)
    >>> f(0,0)
    1
    >>> f(1,2)
    6

Get all pairs from groups
^^^^^^^^^^^^^^^^^^^^^^^^^

Get all of the pairs of items present in a list of groups. A key will be created (i,j) iff i and j share a group

.. doctest::

    >>> from cogent.util.misc import PairsFromGroups
    >>> groups = ['ab','xyz']
    >>> PairsFromGroups(groups)
    {('a', 'a'): None, ('b', 'b'): None, ('b', 'a'): None, ('x', 'y'): None, ('z', 'x'): None, ('y', 'y'): None, ('x', 'x'): None, ('y', 'x'): None, ('z', 'y'): None, ('x', 'z'): None, ('a', 'b'): None, ('y', 'z'): None, ('z', 'z'): None}

Check class types
^^^^^^^^^^^^^^^^^

Check an object against base classes or derived classes to see if it is acceptable

.. doctest::

    >>> from cogent.util.misc import ClassChecker
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Delegate object method calls, properties and variables to the appropriate object. Useful to combine multiple objects together while assuring that the calls will go to the correct object.

.. doctest::

    >>> from cogent.util.misc import Delegator
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Wrap a function to hide it from a class so that it isn't a method. 

.. doctest::

    >>> from cogent.util.misc import FunctionWrapper
    >>> f = FunctionWrapper(str)
    >>> f
    <cogent.util.misc.FunctionWrapper object at ...
    >>> f(123)
    '123'

Construct a constrained container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Wrap a container with a constraint. This is useful for enforcing that the data contained is valid within a defined context. PyCogent provides a base ``ConstrainedContainer`` which can be used to construct user-defined constrained objects. PyCogent also provides ``ConstrainedString``, ``ConstrainedList``, and ``ConstrainedDict``. These provided types fully cover the builtin types while staying integrated with the ``ConstrainedContainer``.

Here is a light example of the ``ConstrainedDict``

.. doctest::

    >>> from cogent.util.misc import ConstrainedDict
    >>> d = ConstrainedDict({'a':1,'b':2,'c':3}, Constraint='abc')
    >>> d
    {'a': 1, 'c': 3, 'b': 2}
    >>> d['d'] = 5
    Traceback (most recent call last):
    ConstraintError: Item 'd' not in constraint 'abc'

PyCogent also provides mapped constrained containers for each of the default types provided, ``MappedString``, ``MappedList``, and ``MappedDict``. These behave the same, except that they map a mask onto ``__contains__`` and ``__getitem__``

.. doctest::

    >>> def mask(x):
    ...   return str(int(x) + 3)
    ... 
    >>> from cogent.util.misc import MappedString
    >>> s = MappedString('12345', Constraint='45678', Mask=mask)
    >>> s
    '45678'
    >>> s + '123'
    '45678456'
    >>> s + '9'
    Traceback (most recent call last):
    ConstraintError: Sequence '9' doesn't meet constraint

Check the location of an application
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Determine if an application is available on a system

.. doctest::

    >>> from cogent.util.misc import app_path
    >>> app_path('ls')
    '/bin/ls'
    >>> app_path('does_not_exist')
    False
