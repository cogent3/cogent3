.. _coding-guidelines:

Coding guidelines
=================

As project size increases, consistency increases in importance. Unit testing and a consistent style are critical to having trusted code to integrate. Also, guesses about names and interfaces will be correct more often.

What should I call my variables?
--------------------------------

- *Choose the name that people will most likely guess.* Make it descriptive, but not too long: ``curr_record`` is better than ``c``, or ``curr``, or ``current_genbank_record_from_database``.

- *Good names are hard to find.* Don't be afraid to change names except when they are part of interfaces that other people are also using. It may take some time working with the code to come up with reasonable names for everything: if you have unit tests, it's easy to change them, especially with global search and replace.

- *Use singular names for individual things, plural names for collections.* For example, you'd expect ``self.name`` to hold something like a single string, but ``self.names`` to hold something that you could loop through like a list or dict. Sometimes the decision can be tricky: is ``self.index`` an int holding a positon, or a dict holding records keyed by name for easy lookup? If you find yourself wondering these things, the name should probably be changed to avoid the problem: try ``self.position`` or ``self.look_up``.

- *Don't make the type part of the name.* You might want to change the implementation later. Use ``Records`` rather than ``RecordDict`` or ``RecordList``, etc. Don't use Hungarian Notation either (i.e. where you prefix the name with the type).

- *Make the name as precise as possible.* If the variable is the name of the input file, call it ``infile_name``, not ``input`` or ``file`` (which you shouldn't use anyway, since they're keywords), and not ``infile`` (because that looks like it should be a file object, not just its name).

- *Use* ``result`` *to store the value that will be returned from a method or function.* Use ``data`` for input in cases where the function or method acts on arbitrary data (e.g. sequence data, or a list of numbers, etc.) unless a more descriptive name is appropriate.

- *One-letter variable names should only occur in math functions or as loop iterators with limited scope.* Limited scope covers things like ``for k in keys: print k``, where ``k`` survives only a line or two. Loop iterators should refer to the variable that they're looping through: ``for k in keys, i in items``, or ``for key in keys, item in items``. If the loop is long or there are several 1-letter variables active in the same scope, rename them.

- *Limit your use of abbreviations.* A few well-known abbreviations are OK, but you don't want to come back to your code in 6 months and have to figure out what ``sptxck2`` is. It's worth it to spend the extra time typing ``species_taxon_check_2``, but that's still a horrible name: what's check number 1? Far better to go with something like ``taxon_is_species_rank`` that needs no explanation, especially if the variable is only used once or twice.

Acceptable abbreviations
^^^^^^^^^^^^^^^^^^^^^^^^

The following list of abbreviations can be considered well-known and used with impunity within mixed name variables, but some should not be used by themselves as they would conflict with common functions, python built-in's, or raise an exception. Do not use the following by themselves as variable names: ``dir``,  ``exp`` (a common ``math`` module function), ``in``, ``max``, and ``min``. They can, however, be used as part of a name, eg ``matrix_exp``.

+--------------------+--------------+
|        Full        |  Abbreviated |
+====================+==============+
|          alignment |          aln |
+--------------------+--------------+
|           archaeal |         arch |
+--------------------+--------------+
|          auxillary |          aux |
+--------------------+--------------+
|          bacterial |         bact |
+--------------------+--------------+
|           citation |         cite |
+--------------------+--------------+
|            current |         curr |
+--------------------+--------------+
|           database |           db |
+--------------------+--------------+
|         dictionary |         dict |
+--------------------+--------------+
|          directory |          dir |
+--------------------+--------------+
|        end of file |          eof |
+--------------------+--------------+
|         eukaryotic |          euk |
+--------------------+--------------+
|          frequency |         freq |
+--------------------+--------------+
|           expected |          exp |
+--------------------+--------------+
|              index |          idx |
+--------------------+--------------+
|              input |           in |
+--------------------+--------------+
|            maximum |          max |
+--------------------+--------------+
|            minimum |          min |
+--------------------+--------------+
|      mitochondrial |           mt |
+--------------------+--------------+
|             number |          num |
+--------------------+--------------+
|           observed |          obs |
+--------------------+--------------+
|           original |         orig |
+--------------------+--------------+
|             output |          out |
+--------------------+--------------+
|          parameter |        param |
+--------------------+--------------+
|          phylogeny |        phylo |
+--------------------+--------------+
|           previous |         prev |
+--------------------+--------------+
|        probability |         prob |
+--------------------+--------------+
|            protein |         prot |
+--------------------+--------------+
|             record |          rec |
+--------------------+--------------+
|          reference |          ref |
+--------------------+--------------+
|           sequence |          seq |
+--------------------+--------------+
| standard deviation |        stdev |
+--------------------+--------------+
|         statistics |        stats |
+--------------------+--------------+
|             string |          str |
+--------------------+--------------+
|          structure |       struct |
+--------------------+--------------+
|          temporary |         temp |
+--------------------+--------------+
|          taxonomic |          tax |
+--------------------+--------------+
|           variance |          var |
+--------------------+--------------+

What are the naming conventions?
--------------------------------

We aim to adhere, to a large extent, to `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_!

How should I write comments?
----------------------------

.. todo::

    refer to the numpy way of documenting
    
- *Always update the docstring when the code changes.* Like outdated comments, outdated docstrings can waste a lot of time. "Correct examples are priceless, but incorrect examples are worse than worthless." `Jim Fulton`_.

How should I format my code?
----------------------------

- *Use 4 spaces for indentation.* Do not use tabs (set your editor to convert tabs to spaces). The behaviour of tabs is not predictable across platforms, and will cause syntax errors. If we all use the same indentation, collaboration is much easier.

- *Lines should not be longer than 79 characters.* Long lines are inconvenient in some editors. Use \\ for line continuation. Note that there cannot be whitespace after the \\.

- *Blank lines should be used to highlight class and method definitions.* Separate class definitions by two blank lines. Separate methods by one blank line.

How should I test my code ?
---------------------------

.. TODO update to refer to more recent discussions on testing in python

Tests are an opportunity to invent the interface(s) you want. Write the test for a method before you write the method: often, this helps you figure out what you would want to call it and what parameters it should take. It's OK to write the tests a few methods at a time, and to change them as your ideas about the interface change. However, you shouldn't change them once you've told other people what the interface is.

Never treat prototypes as production code. It's fine to write prototype code without tests to try things out, but when you've figured out the algorithm and interfaces you must rewrite it *with tests* to consider it finished. Often, this helps you decide what interfaces and functionality you actually need and what you can get rid of.

"Code a little test a little". For production code, write a couple of tests, then a couple of methods, then a couple more tests, then a couple more methods, then maybe change some of the names or generalize some of the functionality. If you have a huge amount of code where 'all you have to do is write the tests', you're probably closer to 30% done than 90%. Testing vastly reduces the time spent debugging, since whatever went wrong has to be in the code you wrote since the last test suite. And remember to use python's interactive interpreter for quick checks of syntax and ideas.

Run the test suite when you change `anything`. Even if a change seems trivial, it will only take a couple of seconds to run the tests and then you'll be sure. This can eliminate long and frustrating debugging sessions where the change turned out to have been made long ago, but didn't seem significant at the time.

Some ``unittest`` pointers
^^^^^^^^^^^^^^^^^^^^^^^^^^

- *Use the* ``unittest`` *framework with tests in a separate file for each module.* Name the test file ``test_module_name.py``. Keeping the tests separate from the code reduces the temptation to change the tests when the code doesn't work, and makes it easy to verify that a completely new implementation presents the same interface (behaves the same) as the old.

- *Use* ``evo.unit_test`` *if you are doing anything with floating point numbers or permutations* (use ``assertFloatEqual``). Do *not* try to compare floating point numbers using ``assertEqual`` if you value your sanity. ``assertFloatEqualAbs`` and ``assertFloatEqualRel`` can specifically test for absolute and relative differences if the default behavior is not giving you what you want. Similarly, ``assertEqualItems``, ``assertSameItems``, etc. can be useful when testing permutations.

- *Test the interface of each class in your code by defining at least one* ``TestCase`` *with the name* ``ClassNameTests``. This should contain tests for everything in the public interface.

- *If the class is complicated, you may want to define additional tests with names* ``ClassNameTests_test_type``. These might subclass ``ClassNameTests`` in order to share ``setUp`` methods, etc.

- *Tests of private methods should be in a separate* ``TestCase`` *called* ``ClassNameTests_private``. Private methods may change if you change the implementation. It is not required that test cases for private methods pass when you change things (that's why they're private, after all), though it is often useful to have these tests for debugging.

- *Test `all` the methods in your class.* You should assume that any method you haven't tested has bugs. The convention for naming tests is ``test_method_name``. Any leading and trailing underscores on the method name can be ignored for the purposes of the test; however, *all tests must start with the literal substring* ``test`` *for* ``unittest`` *to find them.* If the method is particularly complex, or has several discretely different cases you need to check, use ``test_method_name_suffix``, e.g. ``test_init_empty``, ``test_init_single``, ``test_init_wrong_type``, etc. for testing ``__init__``.

- *Write good docstrings for all your test methods.* When you run the test with the ``-v`` command-line switch for verbose output, the docstring for each test will be printed along with ``...OK`` or ``...FAILED`` on a single line. It is thus important that your docstring is short and descriptive, and makes sense in this context.

    **Good docstrings:** ::

        NumberList.var should raise ValueError on empty or 1-item list
        NumberList.var should match values from R if list has >2 items
        NumberList.__init__ should raise error on values that fail float()
        FrequencyDistribution.var should match corresponding NumberList var

    **Bad docstrings:** ::

        var should calculate variance           # lacks class name, not descriptive
        Check initialization of a NumberList    # doesn't say what's expected
        Tests of the NumberList initialization. # ditto

- *Module-level functions should be tested in their own* ``TestCase``\ *, called* ``modulenameTests``. Even if these functions are simple, it's important to check that they work as advertised.

- *It is much more important to test several small cases that you can check by hand than a single large case that requires a calculator.* Don't trust spreadsheets for numerical calculations -- use R instead!

- *Make sure you test all the edge cases: what happens when the input is None, or '', or 0, or negative?* What happens at values that cause a conditional to go one way or the other? Does incorrect input raise the right exceptions? Can your code accept subclasses or superclasses of the types it expects? What happens with very large input?

- *To test permutations, check that the original and shuffled version are different, but that the sorted original and sorted shuffled version are the same.* Make sure that you get *different* permutations on repeated runs and when starting from different points.

- *To test random choices, figure out how many of each choice you expect in a large sample (say, 1000 or a million) using the binomial distribution or its normal approximation.* Run the test several times and check that you're within, say, 3 standard deviations of the mean.

.. _`Jim Fulton`: http://www.python.org/pycon/dc2004/papers/4/PyCon2004DocTestUnit.pdf