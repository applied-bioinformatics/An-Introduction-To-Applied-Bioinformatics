
# Phylogenetic reconstruction <link src='adb2e7'/>

In this chapter we'll begin to explore the goals, approaches, and challenges for creating phylogenetic trees, or phylogenies. Phylogenies, such as the two presented in Figure 1, represent hypotheses about the evolutionary history of a group of individuals, who are represented by the *tips* in the tree. You can explore an interactive version of the three-domain tree presented in Figure 1b online, through the [Interactive Tree of Life project](http://itol.embl.de/itol.cgi#).

** Figure 1: Add (a) the Darwin OoS tree and (b) the Pace three domain tree. **

## Why build phylogenies? <link src="Q4cFRp"/>

Reconstructing the phylogeny of a group of individuals is useful for many reasons. Probably the most obvious of these is understanding the evolutionary relationship between a group of organisms. For example, over the past half-century we've gained great insight into the evolution of our species, *Homo sapiens*, by studying features of our closest relatives, both extant (still existing) organisms, such as the *Pan* (chimpanzees) and *Gorilla* genera, and extinct (no longer living) species, including *Homo neanderthalensis*, *Homo erectus*, *Homo habilus*, and many species in the *Australopithecus* genus. (The [Smithsonian Museum's Human Origins Initiative](http://humanorigins.si.edu/) is an excellent resource for learning more about the fascinating subject of human evolution.) Phylogenetic trees such as these are also useful for understanding evolution itself. In fact, they're so useful that the single image that Charles Darwin found important enough to include in *On the Origin of Species* was the phylogenetic tree presented in Figure 1a.

The *individuals* represented at the tips of our trees don't necessarily have to be organisms though. In another important application of phylogenetic trees, we can study the evolution of genes, which can help us gain a deeper understanding of gene function. Specifically, we can learn about families of related genes. A classic example of this is the globin family, which includes the proteins hemoglobin and myoglobin, molecules that can reversibly bind oxygen (meaning they can bind to it, and then let go of it). You've probably heard of hemoglobin (if not globins in general), as this molecule binds to oxygen where it is present in high concentration (such as in your lung) and releases it where it is present in low concentration (such as in the bicep, where it is ultimately used to power your arm). Hemoglobin and myoglobin are paralogs, meaning that they are related by a gene duplication and subsequent divergence. If you were to compare an unknown globin sequence to either of these you could detect homology, but a tree such as the one present in Figure 2, would help you understand the type of homologous relationship (i.e., whether it was orthology or paralogy).

** Figure 2: Add a globin tree with organism and gene names. **

Phylogenetic trees are used for many other diverse applications in bioinformatics, so it's therefore important that a bioinformatican have an understanding of they are built and how they should be interpreted. An additional application that we'll cover in this text is comparing the composition of communities of organisms, but we'll come back to that [later](alias://2bb2cf).

## How phylogenies is reconstructed <link src="nluhSw"/>

Phylogenies are reconstructed using a variety of different algorithms, some of which we'll cover in this chapter. These algorithms all work by comparing a set of *features* of extant organisms, and inferring the evolutionary distance between those organisms based on the similarity of their features. The features that are compared can be nearly anything, for example ranging from


<figure>
    <img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/alignment.png">
    <figcaption>Figure 3: A spider phylogeny. This image is an unmodified version of [Figure 1](https://doi.org/10.7717/peerj.1719/fig-1), presented by [Garrison et al., 2016](https://peerj.com/articles/1719/).</figcaption>
</figure>
<p>


## Simulating evolution <link src="bR7jKb"/>

```python
>>> # import some functions from python's random module - these will
... # be used in the modeling process
... from random import choice, randint
>>> # import some math functions from the numpy library (note that this
... # isn't part of the python standard library)
... from numpy import log10, average
>>> # import argv from the sys module to support basic command line
... # control of this script
... from sys import argv
...
>>> #####
... # Start of function definitions
... #####
...
... def count_differences(sequence1,sequence2):
...      """Count the number of differences between two sequences of the same length
...      """
...      # confirm that the two sequences are the same length and throw
...      # an error if they aren't
...      assert len(sequence1) == len(sequence2), "Sequences differ in length"
...      # initiate a counter for the number of differences
...      result = 0
...      # iterate over the two sequences and count the number of
...      # positions which are not identical
...      for base1,base2 in zip(sequence1,sequence2):
...          if base1 != base2:
...              # this is a commonly used shortcut for incrementing a count
...              # and is equivalent to the following statement
...              # result = result + 1
...              result += 1
...      return result
...
>>> def evolve_seq(sequence,
...                substitution_probability=0.01,
...                mutation_choices=['A','C','G','T']):
...     """Return two child sequences simulating point mutations
...
...        An error occurs with probability substitution_probability
...         independently at each position of each child sequence.
...     """
...     # Generate two lists for storing the resulting sequences
...     r1 = []
...     r2 = []
...
...     range_length = 10 ** (-1 * log10(substitution_probability))
...
...     for base in sequence:
...         if randint(0,range_length) == 0:
...             # a point mutation will occur at this position
...             # what's wrong with the following statement?
...             r1.append(choice(mutation_choices))
...         else:
...             # no point mutation at this position
...             r1.append(base)
...         if randint(0,range_length) == 0:
...             # a point mutation will occur at this position
...             # what's wrong with the following statement?
...             r2.append(choice(mutation_choices))
...         else:
...             # no point mutation at this position
...             r2.append(base)
...     # convert the lists to strings and return them
...     return ''.join(r1), ''.join(r2)
...
>>> def main(root_sequence,generations,verbose=False):
...     # initial some values and perform some basic error checking
...     assert generations > 0, "Must simulate one or more generations."
...     # can you simplify the following test?
...     for base in root_sequence:
...         assert base != 'A' or base != 'C' or base != 'G' or base != 'T',\
...          "Invalid base identified: %s. Only A, C, G, or T are allowed." % base
...     # initialize a list of the previous generations sequences - this gets used
...     # in the for loop below. since we'll start with the first generation of
...     # children, root_sequence is the previous generation's sequence
...     previous_generation_sequences = [root_sequence]
...
...     # iterate over each generation (why do we add one to generations?)
...     for i in range(1,generations+1):
...         # print the generation number and the current number of sequences
...         print("Generation: %d (Number of child sequences: %d)" % (i,2**i))
...         # create a list to store the current generation of sequences
...         current_generation_sequences = []
...         # create a list to store the differences in each current generation
...         # sequence from the root sequence
...         difference_counts = []
...         # iterate over the sequences of the previous generation
...         for parent_sequence in previous_generation_sequences:
...             # evolve two child sequences
...             r1, r2 = evolve_seq(parent_sequence)
...             # count the differences in the first sequence (from root_sequence)
...             r1_diffs = count_differences(root_sequence,r1)
...             # append the count of differences to the list of difference counts
...             difference_counts.append(r1_diffs)
...             # add the new sequence to the list of this generation's sequences
...             current_generation_sequences.append(r1)
...             # count the differences in the second sequence (from root_sequence)
...             r2_diffs = count_differences(root_sequence,r2)
...             # append the count of differences to the list of difference counts
...             difference_counts.append(r2_diffs)
...             # add the new sequence to the list of this generation's sequences
...             current_generation_sequences.append(r2)
...             if verbose:
...                 # if the caller specified verbose output, print the actual sequences
...                 print("  %s %d" % (r1, r1_diffs))
...                 print("  %s %d" % (r2, r2_diffs))
...         # print summary information: the average number of differences in the current
...         # generation from root_sequence
...         print("Mean differences %1.3f\n" % average(difference_counts))
...         # current_generation_sequences becomes the next generation's
...         # previous_generation_sequences
...         previous_generation_sequences = current_generation_sequences
...
...     # upon completion of all generations, return the last generation's sequences
...     return previous_generation_sequences
...
>>> main('ACCGGGGGGAACCCATTTTACACACACACAC', 5)
Generation: 1 (Number of child sequences: 2)
Mean differences 1.000

Generation: 2 (Number of child sequences: 4)
Mean differences 1.250

Generation: 3 (Number of child sequences: 8)
Mean differences 1.500

Generation: 4 (Number of child sequences: 16)
Mean differences 1.562

Generation: 5 (Number of child sequences: 32)
Mean differences 1.750
['ACCGGGGGGAACCCATTTTACTCACACACAC',
 'ACCGGGGGGAACCCATTTTACTCACACACAC',
 'ACCGGGGGGAACCCATTTTACTCACACACAC',
 'ACCGGGGGGAACCCATTTTACTCACACACAC',
 'ACCGGGGGGAACCCATTTTACTCACACACAC',
 'GCCGGGGGGAACCCATTTTACTCACACACAC',
 'ACCGGGGGGAACCCATTTTACTCACACACAC',
 'ACCGGGGGGAACCCATTTTACTCACACACAC',
 'ACCGGTGGGAACCCATTTTACTCACACACAC',
 'ACCGGTGGGAACCCATTTTACTCACACACAC',
 'ACCGGTGGGAACCCAATTTACTCACACACAC',
 'ACCGGTGGGAACCCATTTTACTCACACACAC',
 'ACCGGGGGGAACCCATTTTACTCACCCACAC',
 'ACCGGGGGGAACCCATTTTACTCACACACAC',
 'ACCGGGGGGAACCCATTTTACTTACACACAC',
 'ACCGGGGGGAACCCATTTTACTTACACACAC',
 'ACCGGGGGGTACCCATTTTACCCACACACAC',
 'ACCGGGGGGTACCCATTTTACCCACACACAC',
 'ACCGGGGGGTACCCATTTTACCCACACACAC',
 'ACCGGGGGGTACCCATTTTACCCACACACAC',
 'ACCGGGGGGAACCCATTTTACCCACACACAC',
 'ACCGGGGGGAACCCATTTTACCCACACACAC',
 'ACCGGAGGGAACCCATTTTACCCACACACAC',
 'ACCGGGGGGAACCCATTTTACCCACACACAC',
 'ACCGGGGGGAACCCATTTTATCCACACACAC',
 'ACCGGGGGGAACCCATTTTATCCACACACAC',
 'ACCGGGGGGAACCCATTTTATCCACACACAC',
 'ACCGGGGGGAACCCATTTTATCCACACACAC',
 'ACCGGGGGGAACCCATTTTATCCACACACAC',
 'ACCGGGGGGAACCCATTTTATCCACACACAC',
 'ACCGAGGGGAACCCATTTTATCCCCACACAC',
 'ACCGGGGGGAACCCATTTTATCCACACACAC']
```


** Pick up here **

Old notes:

By comparing traits of extant organisms. In our case, traits are columns in a multiple sequence alignment. Many algorithms and tools exist for achieving this, and they vary widely in runtime and quality of results. We're going to begin by learning about one of the oldest and simplest methods for doing this: *Unweighted Pair Group Method with Arithmetic Mean* or UPGMA. (Don't be scared by the name - it's actually fairly simple.)

UPGMA is a heirarchical clustering algorithm. It is widely used, though it's application in phylogenetics is usually restricted to building preliminary trees to "guide" the process of multiple sequence alignment, as it makes some assumptions that don't work well for inferring relationships between organsims. We're going to start with it here however for a few reasons. First, the underlying math is very basic, so we don't need to assume anything about your background. Second, there are some other applications of UPGMA that we'll explore later, including grouping samples based on their species compositions. In general, my strategy with teaching this material is to start by giving you a basic introdution into how a process works so you can visualize and do it. From there, we can get more complex.



## Some terminology <link src='7bde92'/>

Here's the goal (well, *a* goal, but this is the basic idea):

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/basic-rooted-tree1.jpg" width=600>

Each *leaf* (or *tip*, or *terminal node*) in this tree represents a sequence, and the length of the horizonal branches between them indicate their dissimilarity to one another. This is a *rooted tree*, which means that it includes an assumption about the last common ancestor of all sequences represented in the tree.

Note that the vertical lines in this tree are used for layout purposes only - they do not represent dissimilarity between sequences.

An **unrooted tree**, like the following, doesn't include an assumption about the last common ancestor of all sequences:

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/basic-unrooted-tree1.jpg" width=600>

**Terminal nodes, tips or leaves** extant organisms, frequently called operational taxonomic units or OTUs. OTUs are families of related organisms.

**Internal nodes** hypothetical ancestors - we postulate their existence but often don't have direct evidence.

**Clade** a node and all nodes "below" it (i.e., toward the tips)

**Root** the internal node defining the clade which contains all nodes in the tree

**Branches** representative of the distance between the nodes.

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/tree-schematic1.png">

**Monophyletic group** the last common ancestor was a member of the group (e.g., multicellular organisms)

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/tree-monophyly.png">

**Polyphyletic group** the last common ancestor was not a member of the group (e.g., flying animals)

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/tree-polyphyly.png">

We assume that as time progresses, sequences will diverge such that more similar sequences have diverged more recently. The problem of phylogenetic reconstruction however is that we only have the tips. We don't have sequences for the internal nodes, and so we use modern sequences to develop a hypothesis about the evolutionary history of a sequence (and hopefully of the organisms who encode those sequences in their genomes).

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/sequence-evo-tree.png">

How many (rooted) trees are there for `n` sequences? This topic is discussed in detail in Chapter 3 of [Inferring Phylogenies](http://www.amazon.com/Inferring-Phylogenies-Joseph-Felsenstein/dp/0878931775/ref=sr_1_1?s=books&ie=UTF8&qid=1393288952&sr=1-1&keywords=inferring+phylogenies), the definitive text on this topic, but the basic answer is **a lot**.

Because of the massive number of possible trees for any reasonable number of sequences, we can't just search all trees to figure out which one best matches our data. Instead, we must take a heurtistic approach to tree building. Some design features of heurtistic methods that are used in practice are that they:

1. Look at a subset of the possible trees, and don't guarantee to find the best tree.
2. Scale to trees for many OTUs (how well they scale depends on the method, and there is a lot of variability)
3. Often provide a single tree, so do not include information on how likely other tree topologies are (we'll talk about methods, such as bootstrapping, to address this).

## Distances and distance matrices <link src='c11a4e'/>

Computing a UPGMA tree for a group of sequences relies on first computing *distances* between each pair of those sequences. A *distance*, in this sense, is a technical term. It's a measure of dissimilarity between two items, `x` and `y`, which meets a few criteria:

1. `d(x,y) >= 0` (non-negativity)
2. `d(x,y) = 0` [`iff`](http://en.wikipedia.org/wiki/If_and_only_if) `x = y` (identity of indiscernibles)
3. `d(x,y) = d(y,x)` (symmetry)
4. `d(x,z) <= d(x,y) + d(y,z)` (triangle inequality)

Let's start with something more similar than sequences. Let's compute the distances between points in a Cartesian plane, to explore what each of these mean.

```python
>>> %pylab inline
>>> from __future__ import print_function
...
>>> from IPython.core import page
>>> import matplotlib.pyplot as plt
...
>>> page.page = print
...
>>> x_vals = [1, 4, 2]
>>> y_vals = [3, 3, 6]
>>> coordinates = list(zip(x_vals, y_vals))
...
>>> fig, ax = plt.subplots()
>>> ax.scatter(x_vals, y_vals)
>>> ax.set_xlim(0)
>>> ax.set_ylim(0)
```

```python
>>> from scipy.spatial.distance import euclidean
>>> %psource euclidean
```

Does Euclidean distance appear to meet the above four requirements of a distance metric?

We can now compute the distances between all of these pairs of coordinates, and we can store these in what's known as a *distance matrix*.

```python
>>> dm = []
...
>>> for c1 in coordinates:
...     row = []
...     for c2 in coordinates:
...         row.append(euclidean(c1, c2))
...     dm.append(row)
...
>>> print(dm)
```

The [scikit-bio](https://github.com/biocore/scikit-bio) project defines defines a ``DistanceMatrix`` object that we can use to work with this.

```python
>>> from skbio.stats.distance import DistanceMatrix
...
>>> dm = DistanceMatrix(dm, map(str,range(3)))
```

One feature that this gets us is nicer printing.

```python
>>> print(dm)
```

```python
>>> dm
```

```python
>>> import scipy.spatial.distance
>>> from skbio import DistanceMatrix
>>> dm = DistanceMatrix(scipy.spatial.distance.pdist(coordinates, "euclidean"))
>>> print(dm)
```

The conditions of a distance matrix listed above lead to a few specific features: the distance matrix is symmetric (if you flip the upper triangle over the diagonal, the values are the same as those in the lower triangle), the distance matrix is *hollow* (the diagonal is all zeros), and there are no negative values.

## Computing the distance between pairs of sequences. <link src='9d757a'/>

Most often, distances between pairs of sequences are derived from a multiple sequence alignment. These differ from the pairwise alignments that we've looked at thus far, but use the same underlying algorithms (and we'll be coming back to this in the next chapter).

Let's load up some aligned sequences, and compute a distance matrix. For now, we'll compute distances between the sequences using the ``hamming`` function that we worked with in the pairwise alignment chapter.

```python
>>> from skbio import TabularMSA, DNA
>>> import scipy.spatial.distance
...
>>> aln = TabularMSA([DNA('ACCGTGAAGCCAATAC', {'id': 's1'}),
...                   DNA('A-CGTGCAACCATTAC', {'id': 's2'}),
...                   DNA('AGCGTGCAGCCAATAC', {'id': 's3'}),
...                   DNA('AGGGTGCCGC-AATAC', {'id': 's4'}),
...                   DNA('AGGGTGCCAC-AATAC', {'id': 's5'})])
```

```python
>>> def hamming(seq1, seq2):
...     return float(scipy.spatial.distance.hamming(seq1.values, seq2.values))
...
>>> master_dm = DistanceMatrix.from_iterable(aln, metric=hamming, key='id')
>>> print(master_dm)
```

Once we have these distances, we can cluster the sequences based on their similiaries/dissimilarities. This is the first process that we'll explore for tree building.

**NOTE:** The example below assumes that each value in this distance matrix is multiplied by the sequence length, so we'll do that here and work work with the resulting distance matrix.

```python
>>> master_dm = DistanceMatrix(master_dm.data*16, master_dm.ids)
>>> print(master_dm)
```

## Hierarchical clustering with UPGMA <link src='73d028'/>

Unweighted Pair-Group Method with Arithmetic mean

Unweighted: all tip-to-tip distances contribute equally
Pair-group: all branch points lead to exactly two clades
Arithmetic mean: distances to each clade are the mean of distances to all members of that clade

Steps
-----

1. Identify the smallest distance in the matrix and define a clade containing only those members. Draw that clade, and set the *total* branch length to the distance between the tips.
2. Create a new distance matrix with an entry representing the clade created in step 1.
3. Calculate the distance matrix entries for the clade as the mean distance from each of the tips of the new clade to all other tips in the distance matrix.
4. If there is only one distance (below or above the diagonal) in the distance matrix, use it to connect the remaing clades, and stop. Otherwise repeat step 1.

Let's start, working from the above distance matrix.

Iteration 1
------------

Step 1.1: The smallest distance in the above matrix is `1.00`, between `s4` and `s5`. So, we'll draw that clade and set each branch length to `0.5`.

Step 1.2: Next, we'll create a new, smaller distance matrix where the sequences `s4` and `s5` are now represented by a single clade, `(s4, s5)`.

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/upgma-tree-iter1.png">

```python
>>> iter1_ids = ['s1', 's2', 's3', '(s4, s5)']
>>> iter1_dm = [[0.0,   4.0,  2.0, None],
...             [4.0,   0.0,  3.0, None],
...             [2.0,   3.0,  0.0, None],
...             [None, None, None, None]]
```

Step 1.3: We'll now fill in the values from the new clade to each of the existing sequences (or clades). The distance will be the mean between each pre-existing clade, and each of the sequences in the new clade. For example, the distance between `s1` and `(s4, s5)` is the mean of the distance between `s1` and `s4` and `s1` and `s5`:

```python
>>> import numpy as np
...
>>> np.mean([master_dm[0][3], master_dm[0][4]])
```

Step 1.3 (continued): Similarly, the distance between `s2` and `(s4, s5)` is the mean of the distance between `s2` and `s4` and `s2` and `s5`:

```python
>>> np.mean([master_dm[1][3], master_dm[1][4]])
```

Step 1.3 (continued): And finally, the distance between `s3` and `(s4, s5)` is the mean of the distance between `s3` and `s4` and the distance between `s3` and `s5`:

```python
>>> np.mean([master_dm[2][3], master_dm[2][4]])
```

Step 1.3 (continued): If we fill these values in (note that they will be the same above and below the diagonal) the post-iteration 1 distance matrix looks like the following:

```python
>>> iter1_dm = [[0.0, 4.0, 2.0, 5.5],
...       [4.0, 0.0, 3.0, 5.5],
...       [2.0, 3.0, 0.0, 3.5],
...       [5.5, 5.5, 3.5, 0.0]]
...
>>> iter1_dm = DistanceMatrix(iter1_dm, iter1_ids)
>>> print(iter1_dm)
```

Step 1.4: There is still more than one value below the diagonal, so we start a new iteration.

Iteration 2
------------

Step 2.1: The smallest distance in the above matrix is `2.00`, between `s1` and `s3`. So, we'll draw that clade and set each branch length to `1.0`.

Step 2.2: Next, we'll create a new, smaller distance matrix where the sequences `s1` and `s3` are now represented by a single clade, `(s1, s3)`.

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/upgma-tree-iter2.png">

```python
>>> iter2_ids = ['(s1, s3)', 's2', '(s4, s5)']
...
>>> iter2_dm = [[None, None, None],
...       [None,  0.0, 5.5],
...       [None,  5.5, 0.0]]
```

Step 2.3: We'll now fill in the values from the new clade to each of the existing sequences (or clades). The distance will be the mean between each pre-existing clade, and each of the sequences in the new clade. For example, the distance between `s2` and `(s1, s3)` is the mean of the distance between `s2` and `s1` and `s2` and `s3`:

```python
>>> np.mean([master_dm[1][0], master_dm[1][2]])
```

Step 2.3 (continued): Next, we need to find the distance between `(s1, s3)` and `(s4, s5)`. This is computed as the mean of the distance between `s1` and `s4`, the distance between `s1` and `s5`, the distance between `s3` and `s4`, and the distance between `s3` and `s5`. Note that we are going back to our master distance matrix here for these distances, **not** our iteration 1 distance matrix.

```python
>>> np.mean([master_dm[0][3], master_dm[0][4], master_dm[2][3], master_dm[2][4]])
```

Step 2.3 (continued): We can now fill in all of the distances in our iteration 2 distance matrix.

```python
>>> iter2_dm = [[0.0, 3.5, 4.5],
...             [3.5, 0.0, 5.5],
...             [4.5, 5.5, 0.0]]
...
>>> iter2_dm = DistanceMatrix(iter2_dm, iter2_ids)
>>> print(iter2_dm)
```

Step 2.4: There is still more than one value below the diagonal, so we start a new iteration.

Iteration 3
------------

Step 3.1: The smallest distance in the above matrix is `3.50`, between `(s1, s3)` and `s2`. So, we'll draw that clade and set each branch length to `1.75`.

Step 3.2: Next, we'll create a new, smaller distance matrix where the clade `(s1, s3)` and the sequence `s2` are now represented by a single clade, `((s1, s3), s2)`.

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/upgma-tree-iter3.png">

```python
>>> iter3_ids = ['((s1, s3), s2)', '(s4, s5)']
...
>>> iter3_dm = [[None, None],
...             [None,  0.0]]
```

Step 3.3: We'll now fill in the values from the new clade to each of the existing sequences (or clades). This is computed as the mean of the distance between `s1` and `s4`, the distance between `s1` and `s5`, the distance between `s3` and `s4`, the distance between `s3` and `s5`, the distance between `s2` and `s4`, and the distance between `s2` and `s5`. Again, note that we are going back to our master distance matrix here for these distances, **not** our iteration 1 or iteration 2 distance matrix.

```python
>>> np.mean([master_dm[0][3], master_dm[0][4], master_dm[2][3], master_dm[2][4], master_dm[1][3], master_dm[1][4]])
```

Step 3.3 (continued): We can now fill in all of the distances in our iteration 3 distance matrix.

```python
>>> iter3_dm = [[0.0, 4.8],
...             [4.8, 0.0]]
...
>>> iter3_dm = DistanceMatrix(iter3_dm, iter3_ids)
>>> print(iter3_dm)
```

Step 3.4: At this stage, there is only one distance below the diagonal in our distance matrix. So, we can use that distance to draw the final branch in our tree, setting the total branch length to 4.8.

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/upgma-tree-final.png">

[SciPy](http://www.scipy.org/) contains support for running UPGMA and generating *dendrograms* (or basic tree visualizations). We can apply this to our distance matrix as follows. You can explore other options for hierarchical clustering in SciPy [here](http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html) (see the *routines for agglomerative clustering*).

```python
>>> # scipy.cluster.hierarchy.average is an implementation of UPGMA
... from scipy.cluster.hierarchy import average, dendrogram
>>> lm = average(master_dm.condensed_form())
>>> d = dendrogram(lm, labels=master_dm.ids, orientation='right',
...                link_color_func=lambda x: 'black')
```

## Acknowledgements <link src='99ad11'/>

The material in this section was compiled while consulting the following sources:

1. The Phylogenetic Handbook (Lemey, Salemi, Vandamme)
2. Inferring Phylogeny (Felsenstein)
3. [Richard Edwards's teaching website](http://www.southampton.ac.uk/~re1u06/teaching/upgma/)
