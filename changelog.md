
<a id='changelog-2025.12.10a2'></a>
# Changes in release "2025.12.10a2"

Release to address development dependency affecting dependent projects
on Python 3.14.

## Contributors

- @GavinHuttley
- @rmcar17

<a id='changelog-2025.12.10a1'></a>
# Changes in release "2025.12.10a1"

Releasing this to avoid poor user experience due to the
numpy 2.4 pre-release being incompatible with numba.

## Documentation

- numba released for python 3.14, remove instructions to use the --prerelease option
  as this causes problems by allowing numpy 2.4 to be installed, which breaks numba

<a id='changelog-2025.9.8a11'></a>
# Changes in release "2025.9.8a11"

Updated to support python 3.14 plus minor bug fix.

## Contributors

- @GavinHuttley

## Enhancements

- We now support Python 3.14. Specify `--prerelease=allow` when doing
  `pip install` (the pre-release flag is required for `numba` and
  `llvmlite`, which are both at release candidate stage).

## Bug fixes

- The counts_per_seq() methods on both Alignment and SequenceCollection
  are now consistent in honouring the `exclude_unobserved` argument.

## Deprecations

- Functions to convert `tinydb` data stores marked for removal by version 2026.1.

<a id='changelog-2025.9.8a10'></a>
# Changes in release "2025.9.8a10"

A performance enhancement release.

## Contributors

- @GavinHuttley, import speed improvements
- @rmcar17, code review

## Enhancements

- Improve top level import speed. Done by implementing a
  custom lazy import mechanism for the top level cogent3.
  Additionally, move imports inside functions / methods
  for some modules.

<a id='changelog-2025.9.8a9'></a>
# Changes in release "2025.9.8a9"

Minor tweaks to support piqtree developers.

## API

- Modify args `make_tree(tip_names)` and `Alignment.iter_seqs(seq_order)` to accept any `Iterable[str]`, rather than strictly `list[str]`.

## Contributors

- @rmcar17 generalised type hints on `make_tree(tip_names)` and `Alignment.iter_seqs(seq_order)`

<a id='changelog-2025.9.8a8'></a>
# Changes in release "2025.9.8a8"

This release contains minor bug fix and enhancements.

## Contributors

- @GavinHuttley

## Enhancements

- `cogent3.available_apps()` now displays the software licenses for
  apps 🎉

## Bug fixes

- Removed obsolete xmfa parser registration. This was causing errors in
  some cases.

## Deprecations

- Renamed `<sequence collection>.rename_seqs()` to
  `<sequence collection>.renamed_seqs()` because former
  name implies modifying the instance, rather than returning
  a new instance where the names are new.

<a id='changelog-2025.9.8a7'></a>
# Changes in release "2025.9.8a7"

Version bump due to interrupted upload to pypi.

<a id='changelog-2025.9.8a6'></a>
# Changes in release "2025.9.8a6"

This is a minor enhancement release.

## Contributors

- @rmcar17, for more typing related commits!
- @gavinHuttley, performance improvements for feature querying

## Enhancements

- Modify `SequenceCollection.get_features()` to improve performance
  when the collection has large numbers of sequences.

<a id='changelog-2025.9.8a5'></a>
# Changes in release "2025.9.8a5"

Minor enhancement, major advances to type hinting allowing static
type checking for core classes. This greatly improves the IDE experience,
thanks @rmcar17!

## Contributors

- @GavinHuttley
- @rmcar17, major refactoring of code and type hinting core
  classes and functions🎉🚀

## API

- `Alignment` and `SequenceCollection` now inherit shared methods from
  `CollectionBase`. `Alignment` no longer inherits from `SequenceCollection`,
  this enables static type-checking across these classes.
- `MolType` is no longer a dataclass so we can have a constructor argument
  share the same name as a property.

## Enhancements

- `PhyloNode.renamed_nodes()` takes a dict mapping the current names to a
  new name and returns a new `PhyloNode` instance with the new names. If a
  node name is not present in the dict, it's left as is.

<a id='changelog-2025.9.8a4'></a>
# Changes in release "2025.9.8a4"

This is a minor feature enhancement release.

## Contributors

- @GavinHuttley

## Enhancements

- Add a plugin hook for counting kmers. There are no packages doing this
  ... yet. Unless a package name is specified, the first app for the hook
  is returned.
- Added a `SequenceCollection.count_kmers()` method and support plugins for
  both this and `Sequence.count_kmers()` via `use_hook` argument. Additional
  arguments to the plugin are passed via optional kwargs.

<a id='changelog-2025.9.8a3'></a>
# Changes in release "2025.9.8a3"

## Contributors

- First-time contribution from @iosonofabio 🎉🚀
- @GavinHuttley

## Enhancements

- `PhyloNode.sorted()` now uses an iterative algorithm so it will be
  more efficient on large trees.

## Documentation

- @iosonofabio added a gallery example demonstrating integration of cogent3 trees with
  the third-party tree visualisation package iplotx.

<a id='changelog-2025.9.8a2'></a>
# Changes in release "2025.9.8a2"

Minor feature enhancement release.

## API

- `SeqsDataABC` subclasses now need to implement a `.get_hash()` method that
  takes the seqid and returns a string representing the hash. We make no
  assumptions about the hash algorithm used. This hash is used to identify
  duplicated sequences. In the case of the `cogent` `AlignedSeqsData` class,
  the hash is done off the gapped numpy array.
- `AlignedSeqsDataABC.get_positions()` renamed to `AlignedSeqsDataABC.get_pos_range()`

## Contributors

- @GavinHuttley

## Enhancements

- Added property for `MolType.gapped_missing_alphabet`.
- `AlignedSeqsDataABC.get_positions()` now takes a series of integers
  representing (potentially) disjoint positions
- `AlignedSeqsDataABC.variable_positions()` returns the indices for
  positions that are variable (no filtering by type of variation)
- Added `cogent3.core.sequence.count_kmers()` function. Returns an array
  of counts of k-mers consisting of only canonical moltype characters.
  The function is `numba.jit` decorated and very fast. k-mers containing
  a non-canonical state (e.g. an ambiguity code) are excluded.
- Added `Sequence.count_kmers()`, which provides a convenient interface
  to the `count_kmers()` function.

<a id='changelog-2025.9.8a1'></a>
# Changes in release "2025.9.8a1"

A minor feature enhancement release.

## Contributors

- @GavinHuttley

## Enhancements

- Added `PhyloNode.ladderise()` and the US spelling version
  `PhyloNode.ladderize()`. This feature was requested by @yudalang3. The method
  sorts nodes in ascending order based on their number of descendants. We break
  ties using alphabetical sorting of tip names.
- Added `validate` arguments to the following methods:
  `<alphabet classes>.to_indices()`, `<kmer alphabet>.to_index()`,
  `<moltype>.has_ambiguity()`, `<moltype>.rc()`, `<moltype>.is_degenerate()`,
  `<moltype>.is_gapped()`, `<moltype>.get_degenerate_positions()`,
  `<moltype>.count_degenerate()`. When true, input is checked against the
  `<alphabet>.is_valid()` method and raises an `AlphabetError` if that returns
  `False`.

## Deprecations

- Deleted modules, functions / methods and arguments marked for removal by this
  release 🎉.

<a id='changelog-2025.7.10a10'></a>
# Changes in release "2025.7.10a10"

A minor feature enhancement to support the cogent3-h5seqs project.

## Contributors

- @GavinHuttley

## Enhancements

- SeqCollection.modified attribute. Returns whether the collection is a modification
  of the underlying storage. This is useful for third-party writers. Conditions
  checked are subset of names, renamed sequences, reverse complemented, or sliced
  (Alignment only).

<a id='changelog-2025.7.10a9'></a>
# Changes in release "2025.7.10a9"

Minor bug fix release with some advances to docs.

## Contributors

- @rmcar17, bug reports and more type hinting
- @GavinHuttley, documentation, more type hinting and other

## Bug fixes

- Sequence.is_valid() method was not being tested and, unsurprisingly,
  was broken. Now fixed.

## Documentation

- Added documentation for using the cogent3 parsers for basic sequence formats
  (fasta and GenBank). These return Python primitives.

<a id='changelog-2025.7.10a8'></a>
# Changes in release "2025.7.10a8"

This is a maintenance release with deprecation warnings for statistical
tests that are available in other packages and improvements to the type
hinting in core modules.

## Contributors

- @khiron deprecating statistical tests
- @rmcar17 type hinting more cogent3.core modules: alphabet.py,
  genetic_code.py, location.py, moltype.py

## Deprecations

- The following cogent3 functions in /maths/stats/test.py have been deprecated
in favor of their SciPy, numpy, or statsmodel counterparts:

  - `var`
  - `std`
  - `likelihoods`
  - `posteriors`
  - `bayes_updates`
  - `t_paired`
  - `t_one_sample`
  - `t_two_sample`
  - `_t_test_no_variance`
  - `mc_t_two_sample`
  - `_permute_observations`
  - `t_one_observation`
  - `spearman`
  - `pearson_correlation`
  - `correlation`
  - `correlation_test`
  - `correlation_matrix`
  - `regress`
  - `regress_origin`
  - `regress_R2`
  - `regress_residuals`
  - `stdev_from_mean`
  - `regress_major`
  - `z_test`
  - `z_tailed_prob`
  - `t_tailed_prob`
  - `reverse_tails`
  - `tail`
  - `multiple_comparisons`
  - `multiple_inverse`
  - `multiple_n`
  - `fisher`
  - `f_value`
  - `f_two_sample`
  - `ANOVA_one_way`
  - `MonteCarloP`
  - `sign_test`
  - `ks_test`

<a id='changelog-2025.7.10a7'></a>
# Changes in release "2025.7.10a7"

This release has major changes to tree classes.

## Contributors

- @rmcar17 refactored the cogent3 tree classes.

## Enhancements

- `TreeNode` and `PhyloNode` have been merged into a single class.
- `PhyloNode` has been fully type hinted.
- Partial type hinting for some modules related to `PhyloNode`.
- `PhyloNode` now uses `__slots__`.
- `PhyloNode` now has explicit length and support attributes instead of the
  values being stored in the params dictionary.
- `PhyloNode.descendant_array` now returns a numpy array of bools.
- `PhyloNode.prune` now has a callback parameter for merging params of
  different nodes. This replaces the behaviour where values in the params
  dict were all added.
- `PhyloNode.last_common_ancestor` now throws an error if the nodes given
  do not belong to the same tree.
- Rich dict and json representations of `PhyloNode` have been enhanced to
  allow proper capturing of length and support attributed with params dict.
- `PhyloNode.get_neighbours_except` has been promoted to a public method.
  This is useful for manual tree traversal.
- `PhyloNode.get_tip_names` now defaults to including itself if it is a tip.

## Bug fixes

- `PhyloNode.lowest_common_ancestor` now behaves correctly after
  successive calls.
- `PhyloNode.get_node_names` now propogates `include_self` when `tips_only`
  is True.
- `Alignment.pretty_print()` no longer fails if zero-length alignment and
  wrap argument is not None.

## Deprecations

- `PhyloNode.get_distances` will be removed in 2025.9. Use
  `PhyloNode.tip_to_tip_distances` instead for same behaviour.
- `to_image` on `Drawable` objects will have `format` parameter replaced by
  `format_name` in 2025.9.

## Discontinued

- `cogent3.core.tree.distance_from_r` will be removed in 2025.9
- `cogent3.core.tree.TreeNode` will be removed in 2025.9
- `PhyloNode` copy functions will have parameters `_nil` and `constructor`
  removed in 2025.9 (were unused).
- `PhyloNode.copy_topology` will be removed in 2025.9. Use the existing
  `.copy` method.
- `PhyloNode.child_groups` will be removed in 2025.9.
- `constructor` parameter will be removed from `PhyloNode.unrooted_deepcopy`,
  `PhyloNode.multifurcating` and `PhyloNode.bifurcating` in 2025.9.
- Support for XML representation of trees will be removed in 2025.9. Use
  newick or json instead.

<a id='changelog-2025.7.10a6'></a>
# Changes in release "2025.7.10a6"

## Contributors

- @khiron
- @GavinHuttley

## Bug fixes

- Counting states on sequences and sequence collections now robust to zero-length
  sequences. This fixes the string formatting method failures (e.g. `to_pretty()`).

## Enhancements

- Moved `theoretical_quantiles` and `probability_points` functions
into `cogent3.maths.stats.test.py`. Now using scipy functions for calculations.
Use can specify additional arguments to the scipy functions
as keyword arguments to these functions, eg. `theoretical_quantiles(data, "t", df=2)`.

## Deprecations

- Modules `cogent3.maths.stats.special.py` and `cogent3.maths.stats.distributions.py`
have been deprecated in favour of `scipy` and `numpy` implementations.
The deprecated modules will be removed in release 2025.10.

<a id='changelog-2025.7.10a5'></a>
# Changes in release "2025.7.10a5"

Minor enhancement to support new version of [ensembl-tui](https://github.com/cogent3/ensembl_tui)
which implements a custom annotation db.

## Enhancements

- Sequence.get_features() now supports additional keyword arguments to be
  passed to the bound annotation db. This allows using third-party annotation
  databases with additional features.

<a id='changelog-2025.7.10a4'></a>
# Changes in release "2025.7.10a4"

This is a bugfix release.

## Contributors

- @GavinHuttley

## Bug fixes

- Fixed another edge case affecting feature querying on sequences. If a sequence,
  or a sequence collection, is created with an annotation offset we now consistently
  apply that offset.

<a id='changelog-2025.7.10a3'></a>
# Changes in release "2025.7.10a3"

A minor bug fix release.

## Contributors

- GavinHuttley

## Bug fixes

- Feature querying on alignments now respects user-provided annotation offsets.

<a id='changelog-2025.7.10a2'></a>
# Changes in release "2025.7.10a2"

This is a minor enhancement release to improve the experience of users
who have relied on the old type sequence classes to have a default moltype.

## Enhancements

- No longer fail if a user does not specify the moltype. Instead, we
  provide a UserWarning and default to ASCII.

<a id='changelog-2025.7.10a1'></a>
# Changes in release "2025.7.10a1"

Major changes in this release! The old style core objects have now been removed.
Please report any problems as a [issue ](https://github.com/cogent3/cogent3/issues).

## Contributors

- @GavinHuttley

## Enhancements

- `new_type` migration has now been completed 🚀🎉
  For all modules whose name began with `new_`, we have
  dropped this prefix, replacing the original classes.
  Finally, fully integrated, plugin-backed, modernised
  sequence collections and accessory objects!

## Discontinued

- The `cogent3.core.new_...` modules are now all marked to
  be discontinued as of 2025.9.
- All uses of `format` as a parameter are now deprecated in favour
  of `format_name` to avoid overriding the builtin function of the same
  name. They will be removed by 2025.9.

<a id='changelog-2025.5.8a10'></a>
# Changes in release "2025.5.8a10"

## Contributors

- @GavinHuttley, misc
- @khiron

## Enhancements

- Added support deserialising old type sequence collections and the old type
  Alignment and ArrayAlignment into the new type classes. We do not bring
  over all attributes (e.g. annotation data is ignored). These transformations
  are transparent to the user.

<a id='changelog-2025.5.8a9'></a>
# Changes in release "2025.5.8a9"

This is a bug fix, minor enhancements and improvements under the hood.

## Contributors

- GavinHuttley, assorted

## Enhancements

- Refactoring some methods on the tree classes to use iterative, rather
  than recursive, algorithms. Also made some algorithms more general.
- `tree.get_sub_tree(as_rooted: bool)` argument when True preserves the
  number of children of the resolved sub-tree. The default is to coerce
  the number of children to match the original tree.
- `tree.prune(keep_root: bool)` preserves nodes from the root to the
  first node with >= 2 children.

## Bug fixes

- Improved robutsness of the tree rooting methods. Previous code
  could return trees with nodes that had a single child. This has
  been fixed.

## Deprecations

- The Tree `root()` method is being replaced by `get_root()` as the
  original name is misleading. This will be finalised in release 2025.6.
- Many methods on tree classes that are not being used anywhere else,
  or are not well suited to being on a tree class, are marked for deletion.
  These include `compare_by_tip_distances()`, `scale_branch_lengths()`,
  `set_tip_distances()`. These wll be finalised in release 2025.6.
- The `tree.traverse()` method is being discontinued as its functionality is
  readily replaced with existing methods. To traverse all nodes on a tree
  use `tree.preorder()` or `tree.postorder()`. If you want to walk up and
  down the tree use `tree.pre_and_postorder()`.
- `tree.get_sub_tree(keep_root)` argument is being discontinued and no longer
  has an effect.This is a special case. If you want to preserve the original
  tip-to-root distances, call the new `tree.tip_to_root_distances()` method
  with the names of interest.
- `cogent3.maths.stats.distribution.binomial_exact` has been deprecated in
  favour of scipy function `scipy.stats.binom.pmf`.  `binomial_exact` will be
  discontinued in cogent3 release 2025.9.  Until then `binomial_exact` will use
  `scipy.stats.binom.pmf` internally and floating point values for `successes`
  and `trials` will be truncated to integers using `math.floor` introducing a
  potential breaking change in behaviour.

<a id='changelog-2025.5.8a8'></a>
# Changes in release "2025.5.8a8"

A minor bug fix release.

## Bug fixes

- The newick string for some trees created using the recently introduced
  `rooted(edge_name)` method were not correct if `edge_name` was a tip.

<a id='changelog-2025.5.8a7'></a>
# Changes in release "2025.5.8a7"

This is a minor release with some bugfixes and new features.

## Contributors

- rmcar17, fixed display of support in dendrograms
- GavinHuttley, assorted

## Enhancements

- `PhyloNode.rooted()` and `TreeNode.rooted()` now return trees that bifurcate
  at a named edge. Addresses an issue raised by @xonq -- thanks!
- `Dendrogram(support_is_percent=True)` argument means that thresholds and support
  statistics on dendrograms are treated as percentages. These are rounded to the
  nearest integer for display. This is now the default.

## Bug fixes

- `cogent3.util.parallel.is_master_process()` now works correctly
  when running using MPI. This was affecting the ability to create
  data stores when running with MPI.
- fix for displaying support statistics on dendrograms

## Deprecations

- Deprecate `TreeNode.root()` and `PhyloNode.root()` in favour of
  `<class>.get_root()`.

<a id='changelog-2025.5.8a6'></a>
# Changes in release "2025.5.8a6"

This is a minor enhancements release with most improvements being to the documentation.

## Contributors

- @GavinHuttley

## Enhancements

- Now support loading an annotation database that has been serialised using json and
  has a registered deserialiser.
- `Feature.as_one_span()` and `Feature.shadow()` methods now allow providing a name argument
  to overide the default naming. The latter can result in very long names.

## Documentation

- All the documentation has now been switched to using the `new_type` only.
- Added examples of using third-party storage for plugins and for the `quick_tree()` hook.

<a id='changelog-2025.5.8a5'></a>
# Changes in release "2025.5.8a5"

This is a minor bug fix release with some improvements to consistency under the hood.

## Contributors

- @GavinHuttley

## Enhancements

- Dotplots now show the alignment path whenever Alignment.dotplot() is invoked.
  Previously, the path was only displayed if there gaps, which could be confusing
  since that's behaviour of SequenceCollection.dotplot().
- Added Enum for Strand information. Developers writing their own annotation db
  classes should employ this Enum to ensure that the strand information is consistent.

## Bug fixes

- Alignment paths were incorrectly displayed on dotplots when invoked from the
  dotplot method. This has been fixed.

<a id='changelog-2025.5.8a4'></a>
# Changes in release "2025.5.8a4"

This is a minor bug fix release.

## Bug fixes

- some methods on new type sequence collections with renamed sequences
  could give errors due to failure to propagate the name_map attribute
  correctly. This has been fixed.

<a id='changelog-2025.5.8a3'></a>
# Changes in release "2025.5.8a3"

This is a minor release, with updated documentation in the README.

## Contributors

- @GavinHuttley

## Enhancements

- The aligned sequence view defaults to a view of the gapped sequence.

<a id='changelog-2025.5.8a2'></a>
# Changes in release "2025.5.8a2"

## Contributors

- @GavinHuttley

## Enhancements

- To support third-party plugins for improved format parsers and writer we have
  added boolean properties to the base classes, indicating whether the parser/writer
  supports unaligned or aligned sequences. The default is both. Plugin authors should
  over-ride these as required.

<a id='changelog-2025.5.8a1'></a>
# Changes in release "2025.5.8a1"

## Major announcement

We are changing the migration strategy from old type to new type classes.
At present we have old type and new type implementations for sequences,
sequence collections, alignments, molecular types, alphabets and genetic
codes. At present, one can select the new types using `new_type=True`
or by using the `COGENT3_NEW_TYPE` environment variable. We have established
that it is not viable to support both simultaneously. Therefore,
the first release after July 1st 2025 will remove all of the old type
classes. Arguments specific to the old type classes will be deprecated
at that point. While this is a major change, we have been using these
ourselves consistently and feel confident that the disruption to users
should be small at worst. However, we strongly advise all users to migrate
now and report any errors. To do this, add the following statement to the top
of your scripts.

```python
import os

os.environ["COGENT3_NEW_TYPE"] = "1"
```

## Contributors

- @GavinHuttley

## Enhancements

- Implement plugin architecture for `cogent3` sequence format parsers and
  writers. For developers, third party parsers and writers can be added to
  `cogent3` by creating a plugin. Parser plugin's must inherit from
  `cogent3.parse.sequence.SequenceParserBase`. Writer plugins
  must inherit from `cogent3.format.sequence.SequenceWriterBase`. Look at
  the `cogent3` `pyproject.toml` for the plugin entry points.
  We default to using third party plugins for the sequence file formats,
  falling back to `cogent3`'s internal plugins if no others are available.
- New type `SequenceCollection` and `Alignment` classes now have
  a public `name_map` property. This records the mapping between
  the sequence names in the collection and the names in the
  underlying storage. This is required for using annotations.
  `name_map` is an immutable dictionary.
- Added public `.storage` attribute to new type `SequenceCollection` and
  `Alignment` classes. The attribute is immutable.
- Added functions for describing and loading data sets that are now included
  with `cogent3`! The main `brca1.fasta` alignment and associated Murphy
  tree are now included with the package. Use
  `cogent3.available_datasets()` to display what's available and
  `cogent3.get_dataset()` to create the instance you want.
- Added `duplicated_seqs()` method to new type `SequenceCollection` and
  `Alignments`. As the name implies, it identifies the names of sequences
  that are equal. On `Alignments`, the gapped sequence is used. Added
  companion method `drop_duplicated_seqs()` which does as the name implies.
  The first name in each group of duplicates is retained.
- Added new type `Sequence.sample()` method for generating random sequences.

## Discontinued

- We have changed the API for the `SupportsFeatures` protocol. It no longer
  inherits from `SerialisableType`. Thus, subclasses of `AnnotationDbABC` no
  longer have to implement `to_rich_dict()`, `from_dict()` and
  `to_json()` methods.
- We have changed the API for `AnnotationDbABC`. The following methods
  have been converted from abstract methods to concrete methods:
  - `add_feature()`
  - `add_records()`
  - `subset()`
  - `update()`
  - `union()`
  - These concrete methods raise a `NotImplementedError`. Thus
    developers need not implement them in their subclasses.

### Backwards incompatible changes

- `cogent3.format.alignment` has been renamed `cogent3.format.sequence`.
- `cogent3.util.table` has been moved to `cogent3.core.table`. Usage of the
  top level functions for making and loading tables is unaffected.

<a id='changelog-2025.3.22a4'></a>
# Changes in release "2025.3.22a4"

Minor bug-fix and small enhancement for the ensembl-tui project.

## Enhancements

- new type Alignment.with_masked_annotations() masking is only applied to the
  sequences that have matching annotations.
- new type Alignment.with_masked_annotations(seqid) allows applying a mask to
  only one sequence in an alignment.

## Bug fixes

- Masking annotations on sequences with shadow=True would mask all positions
  if the feature ws missing from the sequence.

<a id='changelog-2025.3.22a3'></a>
# Changes in release "2025.3.22a3"

This is a minor bug-fix release.

## Bug fixes

- New type Alignment methods could fail when on a instance that
  had sequences renamed. We were passing the incorrect values for
  constructing the new AlignSeqsData instance.

<a id='changelog-2025.3.22a2'></a>
# Changes in release "2025.3.22a2"

## Contributors

- First time contributions from @c-jonesy, @Yixuan567, @Qinzi27
  and @je634 👏🚀🎉
- @GavinHuttley bug fixes and maintenance

## Enhancements

- The cogent3.core.annotation.Feature class now gets an xattr property.
  This allows developers to pack additional data into the Feature object
  without compromising the core API.
- cogent3 now supports hook-style plugins🎉. The first hook is for computing
  quick_tree() trees. This method is on our alignment class and on our distance
  matrix class. It currently defaults to the cogent3 quick_tree() app, which
  applies the NJ algorithm. For developers, you can see what hooks are supported
  by looking in cogent3._plugin at functions named `get_<something>_hook()`. All
  third-party hooks must adhere to the cogent3 app interface.

  Developers, register your hook by adding a section `[project.entry-points."cogent3.hook"]`
  section to your `pyproject.toml` file. Assign to the hook name (e.g. `quick_tree`) the string
  reference to the app supporting the hook. For example,
  `quick_tree = "my_package.my_module.my_quick_tree_app"`.

  Users, specify which package you want to use by setting `<obj>.quick_tree(use_hook="<package name>")`.
  For example, `use_hook="my_package"` would invoke the one defined above. Set `use_hook="cogent3"`
  to force using the cogent3 default.

- IndelMap.merge_maps(other: "IndelMap", alignd_indices: bool) handles case where
  `other.gap_pos` are in alignment coordinates. This is essential for the case
  of progressive alignments where the alignment coordinates are sequentially
  adjusted as more profiles are aligned. This is the key change that fixes the
  bug in the progressive aligner.

## Bug fixes

- Dotplot on new_type Alignments now show the alignment path 🎉
- Fixed a MAJOR bug affecting the cogent3 progressive aligner. Gaps
  were being incorrectly merged during the progressive steps. Users
  are advised to re-run any analyses that used the progressive aligner.
  This DOES NOT affect creating alignment objects from existing data,
  e.g. alignments produced by other algorithms.

## Discontinued

- new_alignment.Alignment.quick_tree() no longer supports doing bootstraps.
  This is because bootstrapping is not quick! If you want this, post a feature
  request and we will implement a separate app for this.

<a id='changelog-2024.12.19a2'></a>
# Changes in release "2024.12.19a2"

This is a minor bug-fix release with a minor new feature.

## Enhancements

- Slicing a new type sequence with a feature defaults to applying
  the feature.name attribute to the resulting seq.name. This can be
  controlled by the new apply_name argument to the
  Feature.get_slice() method.

## Bug fixes

- new_alphabet.KmerAlphabet.to_indices() consistently applies alphabet
  dtype to the resulting array.

<a id='changelog-2024.12.19a1'></a>
# Changes in release "2024.12.19a1"

## Contributors

- @GavinHuttley

## Enhancements

- Improved load balancing in parallel execution of `new_alignment.Alignment.distance_matrix()`
  calculations. This can deliver pretty substantial speedups for large numbers of sequences.

- Major effort in improving compatibility of `new_type` core objects with the
  rest of the codebase. There remain some missing capabilities, and some edge-case
  test failures. However, the core functionality is now working.

## Deprecations

- The nominated release for setting `new_type=True` as the default for new
  core objects has been pushed back to 2025.1. This will allow time for
  the remaining issues to be resolved.

- All `annotate_from_gff()` methods are deprecated, they will be removed in 2025.6.
  Users should assign the result of `cogent3.load_annotations()` to an object's
  `annotation_db` attribute, or provide the `annotation_path` argument to the standard
  `load_(un)aligned_seqs()` functions.

<a id='changelog-2024.11.29a1'></a>
# Changes in release "2024.11.29a1"

## Contributors

- Many thousand lines of code --- @KatherineCaley 🚀
- @YapengLang contributed bug reports and bug fixes 👏

## Enhancements

- `cogent3.open_()` now automatically handles .xz and .lzma compressed files.
- sequence / alignment loading functions recognise .phy as phylip format.
- Formatting of fasta records is now much quicker.
- The grand rewrite of alignment classes is ready for use! This is a major effort,
  unifying the capabilities of the `Alignment` and `AlignmentArray` classes into a
  single class. The logic is cleaner, the performance is better, and the API is
  largely (but not entirely) backwards compatible. The new approach gives us the
  foundation for major performance improvements in the future. As with the
  moltype, alphabet, genetic code and SequenceCollection, you can select the new
  class via `make_aligned_seqs()` or `load_aligned_seqs()` by specifying `new_type=True`.
  Please post any bugs or issues to the issue tracker.
- Refactored a subset of the pairwise distance calculation demonstrating the performance
  benefits of the new alignment classes. These are accessible via the Alignment.distance_matrix()
  method, and also include support for parallel execution. They are wicked fast 🚀🚀🚀!

## Discontinued

- We no longer support python 3.9.


<a id='changelog-2024.7.19a9'></a>
# Changes in release "2024.7.19a9"

## Enhancements

- `cogent3.open_()` now automatically handles `.xz` and `.lzma` compressed files.
- Formatting of fasta records is now much quicker.
- Standardise arrow head sizes on annotation feature displays.

<a id='changelog-2024.7.19a8'></a>
# Changes in release "2024.7.19a8"

A single bugfix release.

## Contributors

- Bug identified and fixed by @YapengLang 🚀🤩

## Bug fixes

- Fixed a regression in Trees arising from new name / support split functionality.

<a id='changelog-2024.7.19a7'></a>
# Changes in release "2024.7.19a7"

## Contributors

- @GavinHuttley multiple commits and maintenance
- A first time contribution from @petergoodman 🚀🎉!
- @YapengLang added the key function for parsing support from newick format
  to tree object 💪

## Enhancements

- Provide convenience class methods on both DictArray and DistanceMatrix to simplify
  creation when you have a numpy array and the names corresponding to the dimension
  element labels. These methods are mainly of use for developers.
- Major rewrite of the mutual information based coevolution statistics. These are
  most easily accessed from the `Alignment.coevolution()` method. The rewrite delivers
  several orders of magnitude performance improvement by using `numba.jit` compiled
  functions based on numpy arrays, plus caching of intermediate calculations.
  The speed of the resampled Mutual Information calculation is now near identical to
  that for the simpler MI and NMI statistics.
- The coevolution statistics can now be run in parallel on a single machine.
  This is most easily accessed by setting `<alignment>.coevolution(parallel=True)`.
  The parallelisation is at the level of chunks of position pairs.
- cogent3 now supports parsing from newick format with node support
  statistics. The support measure is stored in `<node>.params["support"]`.

## Bug fixes

- Dotplots of sequences that had gaps in common is now correctly handled. Previously,
  the displayed aligned path was very funky because aligned segments could be
  interrupted.

## Documentation

- A new tutorial on using a nonstationary amino-acid model from @petergoodman.
  This will appear in the examples documentation in the next release. Nice
  showcase of cogent3's ability to disentangle clade specific substitution
  processes. Plus a first time contribution!

## Deprecations

- All the old coevolution functions and statistics are marked for removal
  by the last major release of 2024. Their implementations were pre-2010!,
  and so not well suited to current data sizes. The mutual information based
  statistics are retained, but rewritten to be much faster (see the Enhancement
  section).
- The old coevolution example document is removed as its out-of-date.

<a id='changelog-2024.7.19a6'></a>
# Changes in release "2024.7.19a6"

Mainly a bug fix release. But one important new feature is support for parallelisation
in jupyter notebooks!

## Contributors

- @GavinHuttley

## Enhancements

- Use loky for parallel processing instead of multiprocessing alone. loky is more robust
  and supports parallel execution in jupyter notebooks .. yay!
- Using a feature to slice a sequence or an alignment should now preserve the annotation
  when possible (a feature has a single contiguous span), or removes the annotation_db
  when not querying it could be misleading.

## Bug fixes

- Fixed issue where read only zipped data stores would include files
  starting with a "._"
- Fixed regression affecting app.as_completed()

<a id='changelog-2024.7.19a5'></a>
# Changes in release "2024.7.19a5"

There are some new features and minor bug fixes in this release. 

## Contributors

- @GavinHuttley

## Enhancements

- genbank.iter_genbank_records() is a revised genbank parser that returns
  primitive python types. It defaults to returning the genbank record LOCUS
  values as the name. If the convert_features argument is None, then
  genbank record metadata is left as a string. This can boost performance
  10x over the old MinimalGenbankParser. The default is to convert features
  into a dict as before.
- genbank.minimal_parser() is a new function that provides a minimal
  parser for Genbank files, replacing the MinimalGenbankParser. It uses
  the genbank.iter_genbank_records() function.
- genbank.rich_parser() replaces the RichGenbankParser. It uses
  genbank.minimal_parser() but has the same behaviour as the older version.
- the <app writer>.apply_to() method now accepts logger=False. This turns off
  logging completely.

## Bug fixes

- Fixed a regression affecting setting tip font sizes on dendrograms. You can now
  assign either an integer to `<dendrogram>.tip_font` or a dict,
  e.g. `{"size": 10, "family": "Inconsolata, monospace"}`.

## Deprecations

- genbank.MinimalGenbankParser is being discontinued in favor of
  genbank.minimal_parser (see above).
- genbank.RichGenbankParser is being discontinued in favor of
  genbank.rich_parser (see above).

<a id='changelog-2024.7.19a4'></a>
# Changes in release "2024.7.19a4"

This is a minor bug fix release.

## Contributors

- Contribution from @rmcar17 on fixing an encoding error when loading small newick trees.

## Bug fixes

- Fixed an issue with `load_tree` incorrectly detecting the encoding of very small newick files.

<a id='changelog-2024.7.19a3'></a>
# Changes in release "2024.7.19a3"

This is a minor bug fix release.

## Contributors

- Contributions from @GavinHuttley and @rmcar17 to Cogent3 repo
  maintenance, developer tooling

## Enhancements

- The cogent3 select_translatable() app gets a `frame` argument. This allows
  specifying one frame of reference for the translation for all sequences.

## Bug fixes

- Fixed a regression affecting loading sequences when the path had the ~
  character in it.

<a id='changelog-2024.7.19a1'></a>
# Changes in release "2024.7.19a1"

## Contributors

- Contributions from @YapengLang, @KatherineCaley, @rmcar17, @fredjaya,
  @GavinHuttley to cogent3 issues, code AND code reviews!
- Contributions from @GavinHuttley and @khiron to Cogent3 repo
  maintenance, developer tooling, issues etc.. related to the SMBE2024
  workshop!

## Enhancements

- Defined an abstract base class for Gff records and support users providing
  their own class via the cogent3.parse.gff.gff_parser(make_record) argument.
  This enables users to devise a custom implementation tuned to the features
  of different provider gff data.
- Contributions from @KatherineCaley, @fredjaya and @GavinHuttley on massive
  refactor of core objects concerning genetic codes, molecular types, 
  alphabets, sequences and sequence collections ... huge efforts 💪 🚀 👷🏼‍♀️!
- Introduced alpha versions of new implementations of core objects: `Sequence`,
  `SequenceCollection`, `MolType`, `GeneticCode`, and `Alphabet`. These
  "new-style" objects enhance performance by supporting the access of the
  underlying data in various formats (i.e. numpy arrays, bytes or strings). The
  "new-style" objects can be accessed by setting the `new_type=True` argument
  in top-level functions (`make_seq`, `load_seq`, `make_unaligned_seqs`,
  `get_moltype`, `get_code`). These are not yet the default and are not fully
  integrated into the existing code, however, we encourage experimentation in
  cases where integration with old objects is NOT required and look forward to
  any feedback!

## Bug fixes

- @YapengLang identified and fixed a bug in `cogent3.app.evo.model` where
  settings for upper / lower bounds where being ignored.

<a id='changelog-2024.5.7a1'></a>
# Changes in release "2024.5.7a1"

## Contributors

- @KatherineCaley for reviewing PRs and code contributions
  related to SeqView taking on a more prominent role in recording
  history of operations on sequences.
- First-time contributions from @Yutong-Shao, @Zongjing-Han, @yantonglu,
  @ShunyiYang, @KyleSu12, @Ruizhe-wang, @Firimp 🎉
- Adopted suggested documentation example proposed by Lizhaozhe123,
  with some amendments, for the trim_stop_codons app.
- Multiple contributions from @GavinHuttley
- @rmcar17 identified and fixed a bug 🚀!
- @rmcar17 addressed deprecation of using `end` as a column name
- Implementation of plug-in support :building_construction: by @khiron
- We had many contributors to project maintenance, making
  the contributor experience better 🎉!!
  @rmcar17, @GavinHuttley @fredjaya reviewed PR's 🛠️🚀
  @khiron, @KatherineCaley, @YapengLang, @fredjaya, @rmcar17
  and @GavinHuttley all commented on Issues.
  @khiron, @fredjaya and @rmcar17 contributed to general
  maintenance tasks, including developer docs and troubleshooting
  infrastructure issues.

## Enhancements

- Annotation databases are now preserved after renaming sequences.
  This is made possible because SeqView has the seqid attribute which
  is preserved and is independent of the name of the enclosing Sequence.
- Sequence.parent_coordinates() now returns information related to its
  parent. For example, we have a sequence with name "1" which we then
  sliced and rename the slice to "gene". Calling gene.parent_coordinates()
  would return "1", slice start, slice end, strand.
- methods on likelihood function objects that work with continuous-time
  Markov models no longer fail if the model also has discrete-time edges.
  These include lf.get_lengths_as_ens(), lf.get_annotated_tree(),
  lf.get_paralinear_metric().
- Added new lf.get_ens_tree(). This returns trees with the expected number
  of substitutions as the branch length. On discrete-time edges the branch
  length is set to None. Think of this tree as the true "evolutionary tree".
  Thanks to Von Bing Yap for suggesting this!
- Plugin support that will allow 3rd-party developers to add their custom
  functionality to the cogent3 pipeline as first class apps
- The new IndelMap class uses numpy arrays to store information about gap
  locations for aligned sequences. This will greatly reduce the memory
  overhead for aligned sequences. The class also provides explicit methods
  for inter-converting between sequence and alignment coordinates. An important
  difference to the original Map implementation is that IndelMap is memoryless,
  meaning the history of slice operations is now fully delegated to the SeqView
  class.
- The new FeatureMap class is a slightly modified version of the original Map
  class. It is used solely for handling sequence feature mappings to sequences
  and alignments. Like IndelMap, it is memoryless but it does record its
  orientation with respect to the parent sequence.
- Note that both IndelMap and FeatureMap replace the spans attribute with a
  generator.

## Bug fixes

- The saving of numpy arrays in annotation db's was not cross-platform
  compatible due to differences in default types between OS's. Fixed by
  using numpy.save and numpy.load. A warning is raised if an old style
  format is detected and a function provided for updating the format.
- Alignment.quick_tree() default settings referred to a distance calculator
  that has been removed. Now updated to the pdist calculator.

## Documentation

- added cookbook docs for io apps

## Deprecations

- `end` is a SQL keyword and so should be avoided. It is now replaced
  by `stop`. Annotation db's with the old column name are dynamically
  updated.

## Discontinued

- The cogent3.core.location.Map class is now marked for deprecation. It is
  being replaced by two classes, IndelMap and FeatureMap. The latter has
  largely the same functionality of Map.

<a id='changelog-2024.2.5a1'></a>
# Changes in release "2024.2.5a1"

## Contributors

- @fredjaya, documentation, bug fixes and tests. Thanks @fredjaya 🚀!!
- @KatherineCaley and @khiron, maintenance 🛠️. Thanks @KatherineCaley and @khiron 🏎️!
- @GavinHuttley, miscellaneous 🍥

## Enhancements

- Now support python3.12 🚀
- Pin numpy version to < 2 until we can test against released version.
- AnnotationDb.subset() method. Returns a new instance matching provided
  conditions.
- SeqView.parent_start, SeqView.parent_stop properties report position
  on original sequence. These help keep track of the segment on original
  after slicing operations. They always **return** plus strand orientation.
- SeqView offset argument added to constructor
- SeqView.copy() method. Supports slicing the original data (and
  recording the offset).
- Sequence.parent_coordinates() returns parent_start, parent_stop, strand
  values from underlying SeqView. strand is +/- 1, which is relative to
  the original sequence.
- AnnotationDb methods now expect the strand argument to have value "+"/"-"
  or None.
- make_seq() returns a seq as is if it's already the correct molecular
  type. This preserves the AnnotationDb attribute.

## Bug fixes

- Fixed ArrayAlignment.get_degapped_relative_to(). It was returning the
  transpose of the alignment.

## Documentation

- Improvements to docstrings for some cogent3.apps. We have added code
  snippets that can be copy/paste into a python session.
- Major updates to the developer docs to guide new contributors 🎉! Check
  them out [at the c3dev wiki](https://github.com/cogent3/c3dev/wiki).

## Deprecations

- We drop support for python3.8
- Assorted arguments marked for deprecation

<a id='changelog-2023.12.15a1'></a>
# Changes in release "2023.12.15a1"

## Contributors

- @wjjmjh (aka Stephen Ma) fixed an error in the seq features gallery drawing example. Turns out, if you want a thumbnail, you actually have to write it out 🤦‍♂️!!
- Multiple PRs from first-time contributor Fred Jaya. Thanks @fredjaya 🚀!!
- Kath Caley for the sweet new cogent3 logo 🥳!
- Kath Caley for numerous other commits (bug fixes, enhancements, etc..). Thanks @KatherineCaley 🚀!
- First time contribution from @cosmastech 🎉
- Important new tree metrics from Robert McArthur. Thanks @rmcar17 🚀!

## Enhancements

- Added from cogent3.app.typing import defined_types function. This displays
  the cogent3 defined type hints and the objects they represent.

- added available_apps(name_filter) which allows filtering apps by
  their name (thanks to @cosmastech)

- The new AnnotationDb.subset() method creates a subset of the instance
  with records matching the provided conditions.

- app.as_completed() and app.apply_to() no longer raise a
  ValueError if there's no work to be done.

- Added the ever popular Robinson-Fould tree topology measure, but
  you shouldn't use it except for comparison to the far superior
  (statistically at least) matching distance measures. These are all
  available on the tree objects via a new PhyloNode.tree_distance()
  method (thanks to @rmcar17).

- Added methods `min_pair()` and `max_pair()` to `DistanceMatrix`, to return the names of the sequences with minimum and maximum pairwise distance, respectively (thanks to @KatherineCaley).

## Bug fixes

- The sequence returned from Alignment.get_seq() now reflects
  slice operations applied to the Alignment.

## Documentation

- The sequence features drawing example now has a gallery thumbnail.

- Doc showing how the dotplot can be used to highlight alignment errors.

- Added acknowledgement that this project has received funding support from
  the Australian National University. We are grateful!

## Deprecations

- Table.tolist() is being replaced by Table.to_list()

- Reverse slicing of Alignment and ArrayAlignment are now consistent
  with reverse slicing of a string. Previously reverse slicing an
  Alignment instance, e.g. `aln[5:1]` would reverse complement a
  nucleic acid object, but fail if it was any other molecular type.
  This behaviour was different to ArrayAlignment. For both objects, use
  a normal slice and reverse complement, e.g. `aln[1:5].rc()`.

## Discontinued

- Began a major refactor of the sequence collection classes. The major change
  is they no longer accept fasta formatted string as input. As we have simplified
  the data conversion functions, all previously used public functions have now been
  marked as being discontinued 2024.3.

<a id='changelog-2023.9.22a1'></a>
# Changes in release "2023.9.22a1"

## Contributors

- YapengLang for reporting the bug in computing the ENS for
  non-stationary processes.

## Enhancements

- MinimalFastaParser now removes internal spaces
  in sequence blocks

- Added app.evo.model(tree_func: Optional[Callable]) argument. A callback
  function assigned to tree_func is passed an alignment and returns a tree.
  This could be done via loading a tree from disk for that alignment source,
  or estimation of the tree using the alignment. The argument overrides the
  other tree related arguments (tree, unique_trees).

## Bug fixes

- The calculation of the expected number of substitutions on a tree branch
  was incorrect. It was not using the motif probs from the parent node.

- Degapped sequence and collections now retain features

- Fixed issue #993. We provide a new default_length argument to
  LF.make_likelihood_function() to be applied when a provided tree
  has zero (or no) lengths. This is set to be 1.

<a id='changelog-2023.7.18a1'></a>
# Changes in release "2023.7.18a1"

## Contributers

- Richard Morris
- AoboHill made their first contribution!
- pithirat-horvichien made their first contribution!
- First contributions from Robert McArthur and Vijini Mallawaarachchi!
- Kath Caley, massive contributions to the sequence annotation refactor And the documentation and ...!
- Thanks to Dr Minh Bui, author of IQ-tree, for a sample extended newick output
  from IQ-tree!

## Enhancements

- Added automated changelog management to the project
- serializable deprecation of function and method arguments using decorator
- Added `SequenceCollection.distance_matrix()` method. This method provides
  a mechanism for k-mer based approximation of genetic distances
  between sequence pairs. It is applicable only to DNA or RNA moltypes.
  Sequences are converted into 10-mers and the Jaccard distance is computed
  between them. This distance is converted into an estimate of a proportion
  distance using a 10-th degree polynomial. (That polynomial was derived from
  regression to distances from 116 mammal alignments.) The final step is applying
  JC69 to these approximated proportion distances.

- Robert and Vijini added a function for computing the matching distance
  statistic between tree toplogies.

  `from cogent3.phylo.tree_distance.lin_rajan_moret`

  This is a better behaved statistic than Robinson-Foulds. The original
  author was Dr Yu Lin who tragically passed away in 2022. He was our
  dear friend, colleague and mentor.

  Lin et al. 2012 "A Metric for Phylogenetic Trees Based on Matching"
  IEEE/ACM Transactions on Computational Biology and Bioinformatics
  vol. 9, no. 4, pp. 1014-1022, July-Aug. 2012

- Major rewrite of annotation handling! In short, we now use an in-memory SQlite
  database to store all annotations from data sources such as GFF or GenBank. New
  classes provide an interface to this database to support adding, querying records
  that match certain conditions. The new database is added to `Sequence` or `Alignment`
  instances on a `.annotation_db` attribute. When sequences are part of a collection
  (`SequenceCollection` or `Alignment`) they share the same data base. Features are now
  created on demand via the `Sequence` or `Alignment` instances and behave much as the
  original `_Annotatable` subclasses did. There are notable exception to this, as
  outlined in the deprecated and discontinued sections.
  This approach brings a massive performance boost in terms of both speed and memory
  A microbial genome sequence and all it's annotates can be loaded in less than a second.
- A new `cogent3.load_annotations()` function allows loading an annotation
  db instance from one, or more, flatfiles. If you provide an existing annotation
  db instance, the records are added to that db.

- Capture extended newick formatted node data. This is stored in
  `TreeNode.params["other"]` as a raw string.

- The `tree_align()` function now uses new approximation method for faster
  estimation of distances for a obtaining guide tree. This is controlled by
  the `approx_dists` argument. The additional argument `iters` can be used to
  do multiple iterations, using genetic distances computed from the alignment
  produced by the preceding iteration to generate the new guide tree.

  If `approx_dists` is `False`, or the moltype of chosen model is not a nucleic acid
  compatible type, distances are computed by the slower method of performing
  all pairwise alignments first to estimate the distances.

- Added new alignment quality measures as apps, and the ability to invoke them
  from the `Alignment.alignment_quality()` method. The new apps are the
  Information content measure of Hertz and Stormo (denoted `ic_score`), a
  variant on the the sum of pairs measure of Carillo and Lipman
  (denoted `sp_score`), and the log-liklelihood produced by the cogent3
  progressive-HMM aligner (denoted `cogent3_score`). If these apps cannot
  compute a score (e.g. the alignment has only 1 sequence), the return a
  `NotCompleted` instance. Instances of that class evaluates to `False`.

- Added optional argument `lower` to `app.model()`. This provides a global
  mechanism for setting the lower bound on rate and length parameters.

- `load_unaligned_seqs()` now handles glob patterns. If the filename is a glob
  pattern, assumes a listing of files containing a single sequence. The `load_seq()`
  function is applied to each file and a single `SequenceCollection` is returned. To
  see progress, set `show_progress=True`.

- `Table.joined(col_prefix)` argument allows users to specify the prefix of
  columns corresponding to the second table, i.e.
  `result = table.inner_join(table2, col_prefix="")`
  ensures result.header is the sum of table.header and table2.header
  (minus the index column).

- Added `trim_stop` argument to `get_translation()` methods. This means
  translating DNA to protein can be done with one method call, instead of
  two.

## Documentation

- Thanks to Katherine Caley for awesome new docs describing the
  new feature and annotation DB capabilities!

## Deprecations

- Removed the original Composable app classes and related decorators for
  user based apps. `user_function` and `appify` are replaced by the
  `define_app` decorator.

- The function TreeAlign is to be deprecated in 2023.9 and replaced with tree_align

- Every method that has "annotation" in it is now deprecated with a replacement
  indicated by their deprecation warnings. Typically, there's a new method with the
  name "feature" in it.

- `<collection>.has_terminal_stops()` is being deprecated for
  `<collection>.has_terminal_stop()`, because it returns True if a single
  sequence has a terminal stop.

## Discontinued

- Removed methods on `TreeNode` that are a recursion variant of an
  existing methods. `TreeNode.copy_recursive()`, `TreeNode.traverse_recursive()`
  `TreeNode.get_newick_recursive()` all have standard implementations that can
  be used instead.
- `PhyloNode` inherits from `TreeNode` and is distinguished from it only by
  have a length attribute on nodes. All methods that rely on length
  have now been moved to `PhyloNode`. These methods are `PhyloNode.get_distances()`,
  `PhyloNode.set_max_tip_tip_distance()`, `PhyloNode.get_max_tip_tip_distance()`,
  `PhyloNode.max_tip_tip_distance()`, `PhyloNode.compare_by_tip_distances()`.
  One exception is `TreeNode.get_newick()`. When `with_distance=True`, this
  method grabs the "length" attribute from nodes.
- All methods that do not depend on the length attribute are moved to `TreeNode`.
  These methods are `TreeNode.balanced()`, `TreeNode.same_topology()`,
  `TreeNode.unrooted_deepcopy()`, `TreeNode.unrooted()`, `TreeNode.rooted_at()`,
  `TreeNode.rooted_with_tip()`.

- The `SequenceCollection.annotate_from_gff()` method now accept file
  paths only.

- Renaming a sequence in a sequence collection is not applied
  to annotations. Users need to modify names prior to binding
  annotations.

- Dropping support for attaching / detaching individual annotation
  instances from an alignment.

- Backwards incompatible change! `Sequence` and `Alignment` no longer inherit from
  `_Annotatable`, so the methods and attributes from that mixin class are no longer
  available. (As there was no migration strategy, please let us know if it broke
  your code and need help in updating it.)

  Major differences include: the `.annotations` attribute is gone; individual
  annotations can no longer be copied; annotations are not updated on sequence
  operations (you need to re-query).

<a id='changelog-2023.2.12a1'></a>

# Changes in release "2023.2.12a1"

## Contributors

- Gavin Huttley
- Katherine Caley
- Nick Shahmaras
- Richard Morris

Thanks to dgslos who raised the issue regarding IUPAC consensus. Thanks to users active on the GitHub Discussions!

## Enhancements

- get_object_provenance() now allows builtins
- jaccard() distance measure added and older approach deprecated

### Composable apps

- app_help() and get_app() available as top-level imports.
   - app_help() takes the name of the app as a string and displays its summary, parameters and how to create one using get_app().
   - get_app() creates an app instance given its name as a string and constructor arguments.
- added skip_not_completed keyword parameter to define_app decorator.
   - Some apps need to process NotCompleted instances. The current `app.__call__` method returns these instances immediately without passing through to the apps `.main()` method. The change introduces a semi-private attribute `_skip_not_completed` to the class. If it's False, the instance will be passed to `main()`.
- composable data validation now allows NotCompleted
   - if <app>.input returned a NotCompleted, it was being treated as an invalid data type rather than preserving the original cause for failure. The data validation method now immediately returns a provided NotCompleted instance
- add argument id_from_source to all writer apps for a naming callback
  - It should be a callable that generates a unique ID from input data
  - Defaults to new get_unique_id() function, which extracts base name by calling get_data_source() and processing the result, removing file suffixes identified by get_format_suffixes().
  - this means filename suffixes are dropped more cleanly
- new app to_primitive(). This uses a to_rich_dict() method if available otherwise it just returns the original object.
- new app from_primitive(). This takes a dict and deserialises the object using the standard cogent3 functions.
- new app pickle_it(). Does as the name implies.
- new app unpickle_it(). Does as the name implies.
- new app compress(). Compresses using a provided compress function. Defaults to gzip compress.
- new app decompress(). Deompresses using a provided decompress function. Defaults to gzip decompress.
- new app to_json(). Converts result of to_primitive() to json string.
- new app from_json(). Converts json string to python primitives suitable for from_primitive().
- added DEFAULT_SERIALISER and a corresponding DEFAULT_DESERIALISER app instances
  - these are to_primitive() + pickle_it() (and the reverse)
- app.typing.get_constraint_names() now supports all standard python Sequence built-ins (list, tuple, set).
- add type resolver for nested types
  - function resolves the type tree of nested types and also returns the depth of that type tree
  - ensure custom apps don't have excessive nested types. The motivation for this check is it is difficult to efficiently resolve, so we advise the developer (via a TypeError message) to define a custom class for such complex types. They can then choose to validate construction of those class attributes themselves.

### DataStores

These have been completely rewritten and have different behaviour from the original versions. Deprecation warnings are raised when the old ones are employed.

- Loading and creating data stores should now be done using open_data_store(), a top-level import.
  - It replaces the (now deprecated) get_data_store() function.
  - It adds a mode argument, "r" is read only, "w" write, and "a" append. This function should now be used for all creation of new data store instances.
  - Supports opening in-memory sqlitedb for writing, just use ":memory:" as the data_path. If mode is read only, raises a NotImplementedError.
- added new DataStoreSqlite for a more flexibile data store backed
  - supports all python types via pickling, including compression of that data
  - is part of the standard library
  - uses the new DEFAULT_SERIALISER for serialisation. The corresponding DEFAULT_DESERIALISER can be used for reversing that.
  - specified using the suffix ".sqlitedb" or using ":memory:" for an in memory sqlitedb
  - record_type property is the type of completed records
- DataStore's have completed and not_completed properties
  - Iteration on data stores is across *both* of those
  - Iterate over the completed property for subsequent analyses
- DataStore's have drop_not_completed method
- All data stores record NotCompleted and md5 data
- DataStore's have a .validate() method, which checks all records match their recorded md5.
- DataStore's provide separate methods for writing different types
  - write, write_not_completed, write_log
  - all require keyword style arguments
  - all return a DataMember
- DataStoreABC.validate() now records missing md5
- DataStores now have a summary_not_completed property
- repr(DataStore) now displays the construction statement. The str(DataStore) returns the output previously displayed by repr().

### Alignments

- iupac_consensus() method now allows ignoring gaps using the allow_gaps argument.

## Deprecations

- get_data_store -> open_data_store
- all previous data store, data member, writer, loader classes
- Data store summary_incomplete property is renamed summary_not_completed

## Discontinued

- We are discontinuing support for tinydb.
  - added convert_tinydb_to_sqlite() function for converting old tinydb to sqlitedb
  - adds a log recording the conversion
- All previous data store types are discontinued, use open_data_store() function for getting a data store instead of a direct import path.

## Bug fixes

- progress display in notebooks now works again

# Changes since release 2022.8.24a1

## Contributors

Thanks to our contributors!

### Accepted PRs from

- Gavin Huttley
- KatherineCaley
- Nick Shahmaras
- Xingjian Leng

### Identified a Bug

- StephenRogers1

## Enhancements

Significant refactor of composable apps. This is the backwards incompatible change we warned of in the last release! We now use a decorator `define_app` (on classes or functions) instead of class inheritance. Please see [the c3dev wiki](https://github.com/cogent3/cogent3/wiki/composable-functions) for examples on how to port from old-style to new-style composable apps.

We updated to the latest NCBI versions of genetic codes. Note, the name of genetic code 1 has changed from "Standard Nuclear" to "Standard".

## Bug fixes

- Fix progressive alignment bug when a guide-tree with zero edge lengths was encountered.
- Non-stationary independent tuple models can now be serialised.

## Deprecations

- We have removed support for python 3.7.
- We have made scipy a dependency and begun deprecating statistical functions that are available in scipy. All deprecated functions have a warning that indicates the scipy replacement. The deprecated functions are: combinations, chi_high, chdtri, z_high, z_low function, chi_low, binomial_high, binomial_low, f_high, f_low, t_low and t_high.

# Changes since release 2022.5.25a1

## Contributors

- Gavin Huttley
- Nick Shahmaras

## Notice of upcoming major changes

The definition of composable apps by inheritance from the `Composable` app base class will no longer be supported from 2022.11. This is to be replaced by a decorator which greatly simplifies constructing new composable apps. Guidance on how to port existing code will be [posted to a dedicated cogent3 wiki page](https://github.com/cogent3/cogent3/wiki/composable-functions). The `cogent3` composable apps will continue to be available and to work as per usual.

We will drop support for python 3.7 by release 2022.10.

## API

- The following method names are marked for deprecation. `Sequence.gettype` for `Sequence.get_type`, `Sequence.resolveambiguities` for `Sequence.resolved_ambiguities`.

## Enhancements

- Refactor dotplot, fixes #1060. This is a major performance improvement in dotplot code, making it much faster for long DNA sequences. I've implemented a seed-and-extend algorithm which identifies matching k-mers between sequences and only works on extending those. The key interfaces are `find_matched_paths()` which takes a `SeqKmers()` instance, the two sequences to be compared and returns a `MatchedSeqPaths()` instance.
- Added Sequence methods to compute all k-mers, fixes #1012
- `Sequence.iter_kmers()`, generator yielding all overlapping k-mers
- Enhanced placement of dendrogram scale bars
- `cogent3.open_` now handles urls too! This enables files to be downloaded from the web using the convenience functions, e.g. `load_aligned_seqs("https://path/to/somefile.fasta")` will now work. Many thanks to Nick Shahmaras for assistance in doing this!

# Changes since release 2022.4.20a1

## Contributors

- Gavin Huttley

## Enhancements

- new `cogent3.util.parallel.as_completed()` generator function
  - `as_completed()` wraps MPI or `concurrent.futures` executors and delivers results as they are completed. In contrast, `parallel.imap()` / `parallel.map()` deliver results in the same order as the input series. The advantage of `as_completed()` is the interval of result arrival at the parent process is better distributed.
- new function `cogent3.load_seq()` loads a single sequence from a file
- convert substitution model `__str__` to `__repr__`; more useful since `__repr__` is called also by str().

## Bug fixes

- fixes to `annotation_from_gff()` method on annotatable sequence / alignment objects
  - method would break if GFF records had no ID. This situation is quite common in some Ensembl gff3 files. We generate a "no-id-#" identifier in those cases.
  - we now add features are added to their parent feature.
- improve consistency in setting motif_probs on likelihood function
  - only apply a pseudocount if optimising motif probs and at least one state has zero frequency, default pseudocount is 0.5. Thanks to StephenRogers1 for finding this issue!

## Documentation

- document the API of the new `load_seq()` function

# Changes since release 2022.4.15a1

## Deprecations

- added warning that we will drop support for python 3.7 by 2022.10. This means the developer version will switch to python 3.8 from 2022.6.
- discontinued delimiter argument from parse.table.load_delimited

# Changes since release 2021.10.12a1

## Contributors

- Gavin Huttley
- u6675275

## API

- moved all io related functions classes from util.misc to util.io, indicating their removal after version 2022.4
- app.result objects require source instance of str or pathlib.Path
- fail if users set motif prob optimisation via sm_args in app.evo.model as value is over ridden by the explicit argument, need to block this as effect is major

## Bug fixes

- RichGenbankParser moltype argument now overrides file spec, if provided, this defines the moltype of the returned sequence, otherwise the moltype is determined from the genbank file meta-data
- fix initialise from nested params for codon models
- load_tree now handles pathlib.Path's as input, fixes #991
- writer composable apps apply_to now handles provided logger
- fixed serialisation of multi-locus likelihood functions with constrained motif probs
- support multiple calls of to_rich_dict()
- solve case where optimiser gets an invalid starting vector
- solved case where optimised parameter values are outside bounds

## Deprecations

- removed deprecated function for median, use numpy.median instead
- removed deprecated index argument from table constructors, use index_name instead
- cogent3.math periodicity classes method names pep8, old names retained, with deprecation warnings

## Enhancements

- Drawable.plotly_figure property returns plotly graph object Figure instance
- refactor of cogent3.app.composable.appify so decorated functions can be pickled
- app.evo.model handles sequential fitting of models with mixed process. Sequential fitting now works if lf_args includes specifying edges for using a discrete-time Markov process
- add optimise_motif_probs argument to app.evo.model
- add upper argument to app.evo.model
- now support python 3.10
- added register_model decorator to cogent3.evolve.models. Used for simplifying discovery of canned substitution models. Users can now use this mechanism too for adding their own custom models. Doing this smoothes usage of custom models with cogent3.app.evo.model. A further benefit is the inclusion of a model to the appropriate module attributes is now done automatically.
- generalise Jensen-Shannon calculations to > 2 distributions
- the register_deserialiser class takes a series of strings that serve to uniquely identify the "type" value in a dict to be reconstituted using the decorated function. This enables support for user defined custom json storage.
- add type hint for input paths to most commonly used loaders
- time-heterogeneity support mixed discrete and continuous-time models
- more compact representation of datastore summary_incomplete
- more refinements on summary_logs
- cogent3.app.io.register_datastore_reader enables development of third party readers / loaders to be developed. Registering a reader class requires decorating it with the filename suffix that will distinguish that content type. Still limited to reading from files only.
- improve general stationary model numerical precision tolerance

# Since release 2021.5.7a1

## Contributors

- GavinHuttley

## DEV

- added missing `dev` requires-extras to pyproject.toml for installing all packages required for development

# Since release 2021.5.7a

## Contributors

- GavinHuttley
- jamesmartini
- KatherineCaley

## API

- ValueError if any tips missing in TreeNode.lowest_common_ancestor()
- added index_name argument to Table.to_categorical(), allows specifying the category column and getting the categorical table in one statement.

## Bug fixes

- DataStore.write() requires identifiers end with indicated suffix
- cogent3.app.tree.quicktree() now works for 2 sequences
- Alignment.degap() now preserves sequence names
- cogent3.app.io.load_aligned() handles paml format
- fast_slow_dist results can now be saved by write_tabular, a DistanceMatrix.source attribute is created on-the-fly by the fast_slow_dist calculator, enabling it be written
- Alignment.variable_positions(), always report a position as variable if > 1 non-gap characters are present
- SequenceCollection.dotplot() method defaults handle single sequence
- app.evo.boostrap() can now be composed,

## DEV

- change to using flit for package management. This change requires you `python -m pip install flit`. If you clone this repository and want to do a developer install, you should first remove your existing one

    ````bash
    $ python -m pip uninstall cogent3
    ````

    then

    ```bash
    $ flit install -s --python `which python`
    ```

## Deprecations

- removed WritableZippedDataStore, the zip archive format is inefficient for incremental inclusion of files. Use a tinydb instead.
- replaced interleave_len argument with wrap in sequence format writers
- removed Table.to_rich_html() method, use Table.to_html() instead

## Enhancements

- More robust alignment to reference algorithm. Builds a multiple sequence alignment from a series of pairwise alignments to a reference sequence. cogent3.app.align.align_to_ref() now retains gaps in the reference. This will be modestly slower than previously, but avoids losing information if the choice of reference sequence is a bad one.
- cogent3.app.composable.appify decorator class, simplifies converting a user defined function into a cogent3 composable app
- JSD calculation now uses more accurate math.fsum()
