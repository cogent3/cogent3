Generating application commandlines
===================================

.. sectionauthor:: Daniel McDonald

In this example we will generate application command lines. This tool is useful for creating jobs scripts in a cluster or supercomputer environment, as well as when varying parameters. First, we must create a ``ParameterIterBase`` subclass instance. For this example, we will use the ``ParameterCombinations`` object and will vary parameters from the ``Clustalw`` application controller.

.. doctest::

    >>> from cogent.app.clustalw import Clustalw
    >>> from cogent.app.util import ParameterCombinations

Lets go ahead and vary the parameters ``-gapdist`` and ``-kimura``, specifying that we always want a value for the ``-gapdist`` parameter. ``-gapdist`` is a ``ValuedParameter``, so we must specify a list of values we wish to use. ``-kimura`` is a ``FlagParameter``, so we only need to say that we would like it on (``True``). If the parameter is not specified in ``AlwaysOn`` then the off state will be added to the range of possible values for the parameter.

.. doctest::

    >>> params_to_vary = {'-gapdist':[1,2,3], '-kimura':True}
    >>> always_on = ['-gapdist']
    >>> param_comb = ParameterCombinations(Clustalw, params_to_vary, always_on)

The returned instance is a generator that will yield parameter dictionaries that can be passed to the application controller. For this example, we will instead use the generator to construct command line strings.

.. doctest::
    
    >>> from cogent.app.util import cmdline_generator

To generate command lines, you must specify a ``ParameterIterBase`` subclass instance, the full path to the application, an optional initial binary (such as Python), how to handle inputs and outputs, specify the actual inputs, how to handle ``stdin/stdout/stderr``, and if you would like unique outputs to be created. This sounds like a lot, but not all applications support by PyCogent work the same. This generator is designed to handle every application supported by PyCogent. In this example, we are not specifying how to handle ``stderr`` and ``stdout``. They are by default thrown to ``/dev/null``.

.. note:: Output from printing ``cmd`` is truncated for document formatting

.. doctest::
    
    >>> path_to_bin = '' # we do not need an initial binary
    >>> path_to_cmd = '/usr/bin/clustalw'
    >>> paths_to_inputs = ['/home/user/input1','/home/user/input2']
    >>> path_to_output = '/home/user/output'
    >>> unique_outputs = True
    >>> input_param = '-infile'
    >>> output_param = '-outfile'
    >>> cmd_gen = cmdline_generator(param_comb, PathToBin=path_to_bin, \
    ... PathToCmd=path_to_cmd, PathsToInputs=paths_to_inputs, \
    ... PathToOutput=path_to_output, UniqueOutputs=unique_outputs,\
    ... InputParam=input_param, OutputParam=output_param)
    >>> for cmd in cmd_gen:
    ...   print cmd
    ... 
     /usr/bin/clustalw -align -gapdist=1 -kimura -infile="/home/user/input1" -outfile="/home/...
     /usr/bin/clustalw -align -gapdist=1 -kimura -infile="/home/user/input2" -outfile="/home/...
     /usr/bin/clustalw -align -gapdist=1 -infile="/home/user/input1" -outfile="/home/user/out...
     /usr/bin/clustalw -align -gapdist=1 -infile="/home/user/input2" -outfile="/home/user/out...
     /usr/bin/clustalw -align -gapdist=2 -kimura -infile="/home/user/input1" -outfile="/home/...
     /usr/bin/clustalw -align -gapdist=2 -kimura -infile="/home/user/input2" -outfile="/home/...
     /usr/bin/clustalw -align -gapdist=2 -infile="/home/user/input1" -outfile="/home/user/out...
     /usr/bin/clustalw -align -gapdist=2 -infile="/home/user/input2" -outfile="/home/user/out...
     /usr/bin/clustalw -align -gapdist=3 -kimura -infile="/home/user/input1" -outfile="/home/...
     /usr/bin/clustalw -align -gapdist=3 -kimura -infile="/home/user/input2" -outfile="/home/...
     /usr/bin/clustalw -align -gapdist=3 -infile="/home/user/input1" -outfile="/home/user/out...
     /usr/bin/clustalw -align -gapdist=3 -infile="/home/user/input2" -outfile="/home/user/out...

