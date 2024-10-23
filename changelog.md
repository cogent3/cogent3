
<a id='changelog-2024.7.19a7'></a>
# Changes in release "2024.7.19a7"

## Contributors

- @GavinHuttley multiple commits and maintenance
- A first time contribution from @petergoodman 🚀🎉!
- @YapengLang added the key function for parsing support from newick format
  to tree object 💪

## ENH

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

## BUG

- Dotplots of sequences that had gaps in common is now correctly handled. Previously,
  the displayed aligned path was very funky because aligned segments could be
  interrupted.

## DOC

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

## ENH

- Use loky for parallel processing instead of multiprocessing alone. loky is more robust
  and supports parallel execution in jupyter notebooks .. yay!
- Using a feature to slice a sequence or an alignment should now preserve the annotation
  when possible (a feature has a single contiguous span), or removes the annotation_db
  when not querying it could be misleading.

## BUG

- Fixed issue where read only zipped data stores would include files
  starting with a "._"
- Fixed regression affecting app.as_completed()

<a id='changelog-2024.7.19a5'></a>
# Changes in release "2024.7.19a5"

There are some new features and minor bug fixes in this release. 

## Contributors

- @GavinHuttley

## ENH

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

## BUG

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

## BUG

- Fixed an issue with `load_tree` incorrectly detecting the encoding of very small newick files.

<a id='changelog-2024.7.19a3'></a>
# Changes in release "2024.7.19a3"

This is a minor bug fix release.

## Contributors

- Contributions from @GavinHuttley and @rmcar17 to Cogent3 repo
  maintenance, developer tooling

## ENH

- The cogent3 select_translatable() app gets a `frame` argument. This allows
  specifying one frame of reference for the translation for all sequences.

## BUG

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

## ENH

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

## BUG

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

## ENH

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

## BUG

- The saving of numpy arrays in annotation db's was not cross-platform
  compatible due to differences in default types between OS's. Fixed by
  using numpy.save and numpy.load. A warning is raised if an old style
  format is detected and a function provided for updating the format.
- Alignment.quick_tree() default settings referred to a distance calculator
  that has been removed. Now updated to the pdist calculator.

## DOC

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

## ENH

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

## BUG

- Fixed ArrayAlignment.get_degapped_relative_to(). It was returning the
  transpose of the alignment.

## DOC

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

## ENH

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

## BUG

- The sequence returned from Alignment.get_seq() now reflects
  slice operations applied to the Alignment.

## DOC

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

## ENH

- MinimalFastaParser now removes internal spaces
  in sequence blocks

- Added app.evo.model(tree_func: Optional[Callable]) argument. A callback
  function assigned to tree_func is passed an alignment and returns a tree.
  This could be done via loading a tree from disk for that alignment source,
  or estimation of the tree using the alignment. The argument overrides the
  other tree related arguments (tree, unique_trees).

## BUG

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

## ENH

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

## DOC

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

## BUG

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

## ENH

Significant refactor of composable apps. This is the backwards incompatible change we warned of in the last release! We now use a decorator `define_app` (on classes or functions) instead of class inheritance. Please see [the c3dev wiki](https://github.com/cogent3/cogent3/wiki/composable-functions) for examples on how to port from old-style to new-style composable apps.

We updated to the latest NCBI versions of genetic codes. Note, the name of genetic code 1 has changed from "Standard Nuclear" to "Standard".

## BUG

- Fix progressive alignment bug when a guide-tree with zero edge lengths was encountered.
- Non-stationary independent tuple models can now be serialised.

## DEP

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

## ENH

- Refactor dotplot, fixes #1060. This is a major performance improvement in dotplot code, making it much faster for long DNA sequences. I've implemented a seed-and-extend algorithm which identifies matching k-mers between sequences and only works on extending those. The key interfaces are `find_matched_paths()` which takes a `SeqKmers()` instance, the two sequences to be compared and returns a `MatchedSeqPaths()` instance.
- Added Sequence methods to compute all k-mers, fixes #1012
- `Sequence.iter_kmers()`, generator yielding all overlapping k-mers
- Enhanced placement of dendrogram scale bars
- `cogent3.open_` now handles urls too! This enables files to be downloaded from the web using the convenience functions, e.g. `load_aligned_seqs("https://path/to/somefile.fasta")` will now work. Many thanks to Nick Shahmaras for assistance in doing this!

# Changes since release 2022.4.20a1

## Contributors

- Gavin Huttley

## ENH

- new `cogent3.util.parallel.as_completed()` generator function
  - `as_completed()` wraps MPI or `concurrent.futures` executors and delivers results as they are completed. In contrast, `parallel.imap()` / `parallel.map()` deliver results in the same order as the input series. The advantage of `as_completed()` is the interval of result arrival at the parent process is better distributed.
- new function `cogent3.load_seq()` loads a single sequence from a file
- convert substitution model `__str__` to `__repr__`; more useful since `__repr__` is called also by str().

## BUG

- fixes to `annotation_from_gff()` method on annotatable sequence / alignment objects
  - method would break if GFF records had no ID. This situation is quite common in some Ensembl gff3 files. We generate a "no-id-#" identifier in those cases.
  - we now add features are added to their parent feature.
- improve consistency in setting motif_probs on likelihood function
  - only apply a pseudocount if optimising motif probs and at least one state has zero frequency, default pseudocount is 0.5. Thanks to StephenRogers1 for finding this issue!

## DOC

- document the API of the new `load_seq()` function

# Changes since release 2022.4.15a1

## DEP

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

## BUG

- RichGenbankParser moltype argument now overrides file spec, if provided, this defines the moltype of the returned sequence, otherwise the moltype is determined from the genbank file meta-data
- fix initialise from nested params for codon models
- load_tree now handles pathlib.Path's as input, fixes #991
- writer composable apps apply_to now handles provided logger
- fixed serialisation of multi-locus likelihood functions with constrained motif probs
- support multiple calls of to_rich_dict()
- solve case where optimiser gets an invalid starting vector
- solved case where optimised parameter values are outside bounds

## DEP

- removed deprecated function for median, use numpy.median instead
- removed deprecated index argument from table constructors, use index_name instead
- cogent3.math periodicity classes method names pep8, old names retained, with deprecation warnings

## ENH

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

## BUG

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

## DEP

- removed WritableZippedDataStore, the zip archive format is inefficient for incremental inclusion of files. Use a tinydb instead.
- replaced interleave_len argument with wrap in sequence format writers
- removed Table.to_rich_html() method, use Table.to_html() instead

## ENH

- More robust alignment to reference algorithm. Builds a multiple sequence alignment from a series of pairwise alignments to a reference sequence. cogent3.app.align.align_to_ref() now retains gaps in the reference. This will be modestly slower than previously, but avoids losing information if the choice of reference sequence is a bad one.
- cogent3.app.composable.appify decorator class, simplifies converting a user defined function into a cogent3 composable app
- JSD calculation now uses more accurate math.fsum()
