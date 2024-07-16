<!--
A new scriv changelog fragment.

Uncomment the section that is right (remove the HTML comment wrapper).
-->


### Contributors

- katherine caley (@KatherineCaley) 
- fred jaya (@fredjaya) 
- gavin huttley (@GavinHuttley)


### ENH

- Introduced alpha versions of the new implementations of core objects: `Sequence`, `SequenceCollection`, `MolType`, `GeneticCode`, and `Alphabet`. These "new-style" objects enhance performance by supporting the access of the underlying data in various formats (i.e., numpy arrays, bytes or strings). The "new-style" objects can be accessed by setting the `new_type=True` argument in top-level functions (`make_seq`, `load_seq`, `make_unaligned_seqs`, `get_moltype`, `get_code`). These are not yet the default as they are not fully integrated into the existing code, however, we encourage experimentation in cases where integration with old objects is NOT required and look forward to any feedback!

<!--
### BUG

- A bullet item for the BUG category.

-->
<!--
### DOC

- A bullet item for the DOC category.

-->
<!--
### Deprecations

- A bullet item for the Deprecations category.

-->
<!--
### Discontinued

- A bullet item for the Discontinued category.

-->
