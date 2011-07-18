***********************
Building a tree of life
***********************

.. authors, Greg Caporaso

Building a tree of life with PyCogent
======================================

This cookbook example runs through how to construct construct a tree of life from 16S rRNA sequences to test whether the three domains of life are visible as three separate clusters in a phylogenetic tree. This example covers compiling sequences, building a multiple sequence alignment, building a phylogenetic tree from that sequence alignment, and visualizing the tree. 

Step 0. Set up your python environment
--------------------------------------

For this tutorial you'll need cogent, muscle, and FastTree installed on your system.

Start an interactive python session by entering the following into a command terminal::

	python

You should now see the python command prompt::

	>>>

Step 1: Download sequences from NCBI
------------------------------------

Here we'll work with archaeal, bacteria, and eukaryotic sequences obtained from NCBI using the PyCogent EUtils wrappers. Run the following commands to obtain these sequences::

	from cogent.db.ncbi import EUtils
	e = EUtils()
	arc16s = list(MinimalFastaParser(e['"small subunit rRNA"[ti] AND archaea[orgn]']))
	bac16s = list(MinimalFastaParser(e['"small subunit rRNA"[ti] AND bacteria[orgn]']))
	euk16s = list(MinimalFastaParser(e['"small subunit rRNA"[ti] AND eukarya[orgn]']))

You can check how many sequences you obtained for each query by running::

	len(arc16s)
	len(bac16s)
	len(euk16s)



.. note:: In this example you'll notice that you have relatively few sequences for each query. You'd obtain many more if you replaced the ``rRNA`` in the query with ``ribosomal RNA``, but the runtime would also be significantly longer. For the purpose of these tutorial we'll therefore stick with this command that returns fewer sequences.

Step 2: Load the sequences
--------------------------

We'll begin by loading the sequences that have been downloaded, applying a filter to retain only those that we consider to be of good quality, and renaming the sequences in the process. Sequences fewer than 750 bases or sequences containing one or more ``N``  characters will be ignored (these represent ambiguous bases).

First, define a function to load and filter the sequences::

	from cogent.parse.fasta import MinimalFastaParser
	
	def load_and_filter_seqs(seqs):
	    filtered_seqs = []
	    id_lookup = {}
	    for seq_id, seq in seqs:
	        # Extract the accession number from the sequence.
	        # Use that as the new sequence identifier.
	        new_id = seq_id.split('|')[1]
	        # Create a mapping from new to old identifiers.
	        id_lookup[new_id] = seq_id
	        if len(seq) > 750 and seq.count('N') < 1:
	            filtered_seqs.append((new_id,seq))
	    return filtered_seqs, id_lookup

Next, load and filter the three sequence sets::

	arc16s_filtered, arc_ids = load_and_filter_seqs(arc16s)
	bac16s_filtered, bac_ids = load_and_filter_seqs(bac16s)
	euk16s_filtered, euk_ids = load_and_filter_seqs(euk16s)
	
	len(arc16s_filtered)
	len(bac16s_filtered)
	len(euk16s_filtered)


Step 3: Select a random subset of the sequences
-----------------------------------------------

Import shuffle from the random module to extract a random collection of sequences::

	from random import shuffle
	shuffle(arc16s_filtered)
	shuffle(bac16s_filtered)
	shuffle(euk16s_filtered)

Select some random sequences from each domain. Note that only a few sequences are chosen to facilitate a quick analysis::

	combined16s = arc16s_filtered[:3] + bac16s_filtered[:10] + euk16s_filtered[:6]
	len(combined16s)

Step 4: Load the sequences into a SequenceCollection object
-----------------------------------------------------------

Use ``LoadSeqs`` to load the unaligned sequences into a ``SequenceCollection`` object::

	from cogent import LoadSeqs, DNA
	seqs = LoadSeqs(data=combined16s,moltype=DNA,aligned=False)

You can explore some properties of this sequence collection. For example, you can count how many sequences are in the sequence collection object::

	seqs.getNumSeqs()

Step 5: Align the sequences using muscle
----------------------------------------

Load an aligner function, and align the sequences. Here we'll align with muscle via the muscle application controller. The sequences will be loaded into an ``Alignment`` object called ``aln``.
::

	from cogent.app.muscle import align_unaligned_seqs
	aln = align_unaligned_seqs(seqs,DNA)

Step 6: Build a tree from the alignment using FastTree
------------------------------------------------------

Load a tree-building function, and build a tree from the alignment. Here we'll use FastTree. The tree will be stored in a ``PhyloNode`` object called ``tree``.
::

	from cogent.app.fasttree import build_tree_from_alignment
	tree = build_tree_from_alignment(aln,DNA)

Step 7: Visualize the tree
------------------------------------------

Load a drawing function to generate a prettier picture of the tree::

	from cogent.draw.dendrogram import UnrootedDendrogram 
	dendrogram = UnrootedDendrogram(tree)

Have a quick look at the unrooted dendrogram::

	dendrogram.showFigure()

You should see something like this:

	.. image:: ../images/tol_not_gap_filtered.png
	   :width: 700

Figure 1: A tree of life build from 16S rRNA sequences. The tree domains are not clearly distinct.


Step 8: Filter highly gapped positions from the alignment
---------------------------------------------------------

To try to improve the quality of the alignment and therefore the tree, it's often a good idea to removed positions that contain a high proportion of gap characters from the alignment. These generally represent non-homologous regions of the sequence of interest, and therefore contribute little to our understanding of the evolutionary history of the sequence.

To remove positions that are greater than 10% gap characters from the alignment, run the following command::

	gap_filtered_aln = aln.omitGapPositions(allowed_gap_frac=0.10)

If you count the positions in both the full and reduced alignments you'll see that your alignment is now a lot shorter::

	len(aln)
	len(gap_filtered_aln)

Step 9: Rebuild the tree and visualize the result
-------------------------------------------------

Rebuild the tree and visualize the result as before::

	gap_filtered_tree = build_tree_from_alignment(gap_filtered_aln,DNA)
	gap_filtered_dendrogram = UnrootedDendrogram(gap_filtered_tree)
	gap_filtered_dendrogram.showFigure()

You should now see something that much more clearly looks like a tree with three distinct groups. For example:

	.. image:: ../images/tol_gap_filtered.png
	   :width: 700

Figure 2: A tree of life build from 16S rRNA sequences. There now appear to be three distinct domains.

Step 10: Save the tree as a PDF
-------------------------------

Finally, you can save this tree as a PDF for sharing or later viewing::

	gap_filtered_dendrogram.drawToPDF('./tol.pdf')

