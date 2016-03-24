
# Phylogenetic reconstruction <link src='adb2e7'/>

In this chapter we'll begin to explore the goals, approaches, and challenges for creating phylogenetic trees, or phylogenies. Phylogenies, such as the two presented in Figure 1, represent hypotheses about the evolutionary history of a group of individuals, who are represented by the *tips* in the tree. You can explore an interactive version of the three-domain tree presented in Figure 1b online, through the [Interactive Tree of Life project](http://itol.embl.de/itol.cgi#).

<table border=0>
<tr>
<td>
<figure>
    <img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/Darwins_tree_of_life_1859.png">
</figure>
</td>
<td>
<figure>
    <img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/Pace_Big_Tree.png" width=50%>
</figure>
</td>
</tr>
<tr >
<td colspan=2>
<b>Figure 1a</b>: Evolutionary tree presented by Charles Darwin in <i>On the Origin of Species</i>. <b>Figure 1b</b>: A hypothesis of evolutionary relationships between the three domains of life. This image was created by the <a href="http://pacelab.colorado.edu/PI_NormPace.html">Norman Pace Laboratory.</a>
</td>
</tr>
</table>

## Why build phylogenies? <link src="Q4cFRp"/>

Reconstructing the phylogeny of a group of individuals is useful for many reasons. Probably the most obvious of these is understanding the evolutionary relationship between a group of organisms. For example, over the past half-century we've gained great insight into the evolution of our species, *Homo sapiens*, by studying features of our closest relatives, both extant (still existing) organisms, such as the *Pan* (chimpanzees) and *Gorilla* genera, and extinct (no longer living) species, including *Homo neanderthalensis*, *Homo erectus*, *Homo habilus*, and many species in the *Australopithecus* genus. (The [Smithsonian Museum's Human Origins Initiative](http://humanorigins.si.edu/) is an excellent resource for learning more about the fascinating subject of human evolution.) In this time, we've also improved our understanding of the deeper branches in the tree of life. For example, [Woese and Fox (1977)](http://www.pnas.org/content/74/11/5088.full.pdf) used phylogenetic reconstruction to first illustrate that the "prokaryotes" really represented two ancient lineages which they called the eubacteria and the archaebacteria, which ultimately led to the proposal of a "three domain" tree of life composed of the three deep branching linages, the archaea, the bacteria, and the eucarya ([Woese, Kandler and Wheelis (1990)](http://www.pnas.org/content/87/12/4576.full.pdf)).

Phylogenetic trees such as these are also useful for understanding evolution itself. In fact, they're so useful that the single image that Charles Darwin found important enough to include in *On the Origin of Species* was the phylogenetic tree presented in Figure 1a.

The *individuals* represented at the tips of our trees don't necessarily have to be organisms though. In another important application of phylogenetic trees, we can study the evolution of genes, which can help us gain a deeper understanding of gene function. Specifically, we can learn about families of related genes. A classic example of this is the globin family, which includes the proteins hemoglobin and myoglobin, molecules that can reversibly bind oxygen (meaning they can bind to it, and then let go of it). You've probably heard of hemoglobin (if not globins in general), as this molecule binds to oxygen where it is present in high concentration (such as in your lung) and releases it where it is present in low concentration (such as in the bicep, where it is ultimately used to power your arm). Hemoglobin and myoglobin are paralogs, meaning that they are related by a gene duplication and subsequent divergence. If you were to compare an unknown globin sequence to either of these you could detect homology, but a tree such as the one present in Figure 2, would help you understand the type of homologous relationship (i.e., whether it was orthology or paralogy).

<figure>
    <img src="https://static-content.springer.com/image/art%3A10.1186%2F1471-2148-6-31/MediaObjects/12862_2005_Article_216_Fig5_HTML.jpg" width=50%>
    <figcaption><b>Figure 2</b>: A tree representing members of the globin gene family from diverse taxa. This image is an unmodified reproduction of Figure 5 from <a href="http://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-6-31"><i>A phylogenomic profile of globins</i></a> by Vinogradov et al (2006).

<p>

Phylogenetic trees are used for many other diverse applications in bioinformatics, so it's therefore important that a bioinformatican have an understanding of they are built and how they should be interpreted. An additional application that we'll cover in this text is comparing the composition of communities of organisms, but we'll come back to that [later](alias://2bb2cf).

## How phylogenies are reconstructed <link src="nluhSw"/>

Phylogenies are reconstructed using a variety of different algorithms, some of which we'll cover in this chapter. These algorithms all work by comparing a set of *features* of organisms, and inferring the evolutionary distance between those organisms based on the similarity of their features. The features that are compared can be nearly anything that is observable, either from extant organisms or fossilized representatives of extinct organisms.

As an example, let's consider the reconstruction of the phylogeny of spiders (the order Araneae), a hypothesis of which is presented in Figure 3. Of the extant spiders, some are orb-weavers (meaning they spin circular, sticky webs), and others are not. Entomologists have debated whether orb-weaving is a monophyletic trait (meaning that it evolved one time), or whether it is polyphyletic (meaning that it evolved multiple times, such as flight, which has evolved independently in birds, flying dinosaurs, insects, and mammals). If orb-weaving is monophyletic, it would mean that over the course of evolution, extant spiders which don't weave orb webs have lost that ability. Some researchers doubt this as it's a very effective means of catching prey, and losing that ability would likely constitute an evolutionary disadvantage. If orb-weaving is polyphyletic, it would means that in at least two different spider lineages, this trait arose independently, which other researchers consider to be very unlikely due to the complexity of engineering these webs. Examples of the evolution of monophyletic and polyphyletic traits are presented in Figures 3 and 4, respectively.


<figure>
    <img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/tree-monophyly.png">
    <figcaption><b>Figure 3</b>: Example phylogeny illustrating a monophyletic trait shared by a group of organisms. In a monophyletic group, the last common ancestor was also member of the group (e.g., multicellular organisms).</figcaption>
 </figure>

<figure>
    <img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/tree-polyphyly.png">
    <figcaption><b>Figure 4</b>: Example phylogeny illustrating a polyphyletic trait shared by a group of organisms. In a polyphyletic group the last common ancestor was not a member of the group (e.g., flying animals).</figcaption>
</figure>

<p>

Earlier work on understanding the relations between the spider lineages focused on comparing traits that entomologists would observe, for example by watching spiders in action or by dissecting them. For example, in 1986 through 1991, [Johnathan Coddington](http://entomology.si.edu/StaffPages/coddington.html) published several studies that tabulated and compared about sixty features of 32 spider taxa (Coddington J. 1986. The monophyletic origin of the orb web. In: Shear W, ed. Spiders: webs, behavior, and evolution. Stanford, California: Stanford University Press. 319-363; Coddington JA. 1991. Cladistics and spider classification: araneomorph phylogeny and the monophyly of orbweavers (Araneae: Araneomorphae; Orbiculariae) Acta Zoologica Fennica 190:75-87). Features included whether a spider wraps its prey when it attacks (a behavioral trait), and how branched the spider's trachea is (a morphological trait). By determining which spiders were more similar and different across these traits, Dr. Coddington provided early evidence for the hypothesis that orb-weaving is an ancient monophyletic trait.

More recently, several research teams have used features of spider genomes to reconstruct the spider phylogeny (<a href="http://www.cell.com/current-biology/abstract/S0960-9822(14)00750-7">Bond et al., 2014</a>, [Garrison et al., 2016](https://peerj.com/articles/1719/)). Using this approach, the features become the nucleotides observed at particular positions in the genome, which are observed first by sequencing specific genes that the researchers target that are present in all members of the group, and then aligning those sequences with multiple sequence alignment. This has several advantages over feature matrices derived from morphological and behavioral traits, including that many more features can be observed. For example, ([Garrison et al., 2016](https://peerj.com/articles/1719/)), compared approximately 700,000 amino acid positions from nearly 4000 loci around the genomes of 70 spider taxa. Compare the number of features here to the number mentioned in the previous paragraph. These *phylogenomic* studies have further supported the idea that orb-weaving is an ancient monophyletic trait, and have provided much finer scale information on the evolution of spiders. Supported by these data, researchers hypothesize that the loss of orb-weaving might not be that surprising. While it does provide an effective means of catching flying insects, many insects which are potential prey for spiders don't fly. Further, orb webs may attract predators of spiders, as they are easily observable signals of where a spider can be found.

<figure>
    <img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/spider-tree.png">
    <figcaption><b>Figure 5</b>: A spider phylogeny. Numbers at internal nodes correspond to the taxonomic groups described in <a href="https://peerj.com/articles/1719/#table-1">Table 1 of Garrison et al., 2016</a>. This image is an unmodified version of <a href="https://doi.org/10.7717/peerj.1719/fig-1">Figure 1</a> of <a href="https://peerj.com/articles/1719/">Garrison et al., 2016</a>.</figcaption>
</figure>
<p>

For the remainder of this chapter, we'll consider methods for phylogenetic reconstruction that use genome sequence data as features.

## Some terminology <link src='7bde92'/>

Next, let's cover a few terms using the tree diagram in Figure 6.

<figure>
    <img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/tree-schematic1.png">
    <figcaption><b>Figure 6</b>: A schematic of a phylogenetic tree illustrating important terms.
</figure>

<p>

*Terminal nodes or tips* typically represent extant organisms, also frequently called operational taxonomic units or OTUs. OTU is a generic way of referring to a grouping of organisms (such as a species, a genus, or a phylum), without specifically identifying what that grouping is.

*Internal nodes* in a phylogenetic tree represent hypothetical ancestors. We postulate their existence but often don't have direct evidence. The *root node* is the internal node from which all other nodes in the tree descend. This is often referred to as the *last common ancestor (LCA)* of the OTUs represented in the tree. In a universal tree of life, the LCA is often referred to as *LUCA*, the *last universal common ancestor*. All nodes in the tree can be referred to as OTUs.

*Branches* connect the nodes in the tree, and generally represent time or some amount of evolutionary change between the OTUs. The specific meaning of the branches will be dependent on the method that was used to build the phylogenetic tree.

A *clade* in a tree refers to some node (either internal or terminal) and all nodes "below" it (i.e., toward the tips).

We'll use all of these terms below as we begin to explore phylogenetic trees.

## Simulating evolution <link src="bR7jKb"/>

Before we jump into how to reconstruct a phylogeny from DNA sequence data, we're going to perform a simulation of the process of evolution of a DNA sequence. In this simulation, we're going to model sequence evolution with a Python function, and then we're going to run that function to simulate multiple generations of evolution.

Bioinformatics developers often use simulations to understand how their algorithms work, as they uniquely provide an opportunity to know what the correct answer is. This provides a way to compare algorithms to each other to figure out which performs best under which circumstances. In our simulation we're going to have control over the starting sequence, and the probability of incurring a substitution mutation or an insertion/deletion mutation at each position of the sequence in each generation. This would, for example, let us understand whether different algorithms for phylogenetic reconstruction are better or worse for more closely or distantly related sequences.

Our simulation will works as follows. We'll have one function that we primarily interact with called ``evolve_generations``. This will take a starting sequence and the number of generations that we want to simulate. It will also take the probability that we want a substitution mutation to occur at each position in the starting sequence, and the probability that we want either an insertion or deletion mutation to occur at each position. In each generation, every sequence will spawn two new sequences, randomly incurring mutations at the prescribed rates. This effectively simulates a clonal process of reproduction, a form of asexual reproduction common in single cellular organisms, where a parent cell divides into two cells, each containing a copy of the parent's genome with some errors introduced.

Let's inspect this code and then run our simulation beginning with a random sequence.

```python
>>> %pylab inline
...
>>> from IPython.core import page
>>> page.page = print
```

First we'll look at the function used to simulate the evolution of a single sequence. This is where most of the important evolutionary modeling of sequence evolution happens.

```python
>>> from iab.algorithms import evolve_sequence
>>> %psource evolve_sequence
```

Next, take a look at the function that models a single generation of a single sequence. This is where the clonal reproduction (i.e., one parent sequence becoming two child sequences) occurs.

```python
>>> from iab.algorithms import evolve_generation
>>> %psource evolve_generation
```

Finally, take a look at our entry point function. This is where we provide the parameters of our simulation, including the starting sequence, the number of generations, and the mutation probabilities. Notice how each of these functions builds on the prior functions.

```python
>>> from iab.algorithms import evolve_generations
>>> %psource evolve_generations
```

Now we'll run our simulation. We'll start with a random DNA sequence, and then evolve three generations. Before running this, can you predict how many child sequences we'll end up with after three generations?

When we call ``evolve_generations``, we'll pass the parameter ``verbose=True``. This will tell the function to print out some information throughout the process. This will let us inspect our evolutionary process: something that is impossible to do with real sequences.

```python
>>> from iab.algorithms import random_sequence
>>> import skbio
>>> sequence = random_sequence(skbio.DNA, 50)
```

```python
>>> sequences = evolve_generations(sequence, 3, 0.1, 0.05, verbose=True)
```

We now have a new variable, sequences, which contains the child sequences from the last generation. Take a minute to look at the ids of the parent and child sequences above, and the ids of a couple of the final generation sequences. These ids are constructed so that each sequence contains the identifiers of its ancestral sequences, and then either ``1`` or ``2``. Notice that all sequence identifiers start with ``0``, the identifier of the last common ancestor (or our starting sequence) of all of the sequences. These identifiers will help us interpret whether the phylogenies that we reconstruct accurately represent the evolutionary relationships between the sequences.

Also, notice that at this point we only have the sequences from the last generation. We no longer have the ancestral sequences (which would correspond to the internal nodes in the tree). This models the real world, were we only have sequences from extant organisms, but not their ancestors.

Take a minute to compare the two sequences below. What types of mutations happened over the course of their evolution?

```python
>>> print(len(sequences))
```

```python
>>> sequences[0]
```

```python
>>> sequences[-1]
```

In our simulation, each sequence is directly derived from exactly one sequence from the previous generation, and the evolution of all of the sequences traces back to starting sequence that we provided. This means that our final sequences are all homologous. And because we have modeled this process, we know where each sequence fits in relation to all of the other sequences in the phylogeny. Our goal with the algorithms we'll study for the rest of this chapter is to reconstruct that phylogeny given only the last generation of sequences. We'll use the fact that we know the true phylogeny to help us evaluate the relative performance of the different methods.

<figure>
    <img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/sequence-evo-tree.png">
    <figcaption><b>Figure 7</b>: Schematic of a simulated evolutionary process. Bases in red indicate mutation since the last common ancestor. The bottom panel illustrates the real-world equivalent of our final product, where we wouldn't know the true phylogeny (indicated by the dashed branches), the sequence of the last common ancestor, or what positions have changed since the last common ancestor.</figcaption>
</figure>

```python
>>> %matplotlib inline
...
>>> import numpy as np
>>> sequences = evolve_generations(sequence, 10, 0.01, 0.005, verbose=False)
```

```python
>>> sequences = np.random.choice(sequences, 10, replace=False)
```

### A cautionary word about simulations <link src="Wbhke4"/>

While simulations are extremely powerful for comparing algorithms, they can also be misleading. This is because when we model evolution we simplify the evolutionary process. For example, in the simulation above, we assume that the rate of substitution mutations doesn't change in different parts of our phylogeny. Imagine in the real-world that the environment changed drastically for some descendants (for example, if a geological event created new thermal vents in a lake they inhabited, resulting in an increase in mean water temperature), but not for others. The descendants who experience the environmental change might have an increased rate of substitutions as their genomes adapt to the new environment. The increased substitution rate may be temporary or permanent.

If we use our simulation code to evaluate phylogeny reconstruction algorithms, it will tell us nothing about which algorithms better handle different evolutionary rates in different branches of the tree. This is one limitation of our simulation that we know about, but because we don't have a perfect understanding of sequence evolution, there are limitations that we don't know about. For this reason, you always want to understand what assumptions a simulation is making, and consider those when determining how confident you are in the results of an evaluation based on simulation. One assumption that our simulation is making is that the evolutionary rate is constant across all branches of the tree. What are some other assumptions that are being made?

On the opposite end of the spectrum from simulations for algorithm comparison is comparisons based on real data. The trade-off however is that with real data we don't know what the right answer is (in our case, the correct phylogeny) so it's harder to determine which algorithms are doing better or worse. The take-away message here is that neither approach is perfect, and often researchers will use a combination of simulated and real data to evaluate algorithms.

**NOTE: The text below here gets very rough as these sections are currently being written or re-written.**

## Parsimony-based approaches to phylogenetic reconstruction <link src="7WKikt"/>

This section is currently a placeholder. You can track progress on this section through [issue #119](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/issues/119). In the meantime, I recommend Chapter 8 of *[The Phylogenetic Handbook](http://www.amazon.com/gp/product/0521730716/ref=as_li_tl?ie=UTF8&camp=1789&creative=9325&creativeASIN=0521730716&linkCode=as2&tag=anintrotoappl-20&linkId=YLNAKVFX7BV4W5TW")*, by Lemey, Salemi, and Vandamme for discussion of this topic.

### How many possible phylogenies are there for a given collection of sequences? <link src="DguiCU"/>

This section is currently a placeholder. You can track progress on this section through [issue #119](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/issues/119). In the meantime, I recommend Chapter 3 of [Inferring Phylogenies](http://www.amazon.com/Inferring-Phylogenies-Joseph-Felsenstein/dp/0878931775/ref=sr_1_1?s=books&ie=UTF8&qid=1393288952&sr=1-1&keywords=inferring+phylogenies), the definitive text on this topic.

## Distance-based approaches to phylogenetic reconstruction <link src="XrznYP"/>

The first approaches we'll take for phylogenetic reconstruction rely on computing distances between sequences. **This section is currently being written. What currently follows is an outline derived from an earlier version of this chapter.** You can track progress on this section through [issue #119](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/issues/119). In the meantime, an additional resource that you may find useful is the UPGMA section of [Richard Edwards's teaching website](http://www.southampton.ac.uk/~re1u06/teaching/upgma/).

### Distances and distance matrices <link src='c11a4e'/>

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

### Alignment-free distances between sequences <link src='9d757a'/>

```python
>>> from iab.algorithms import guide_tree_from_sequences
...
>>> guide_tree_from_sequences(sequences, display_tree=True)
```

### Alignment-based distances between sequences <link src="xAIbfm"/>

```python
>>> from skbio.alignment import global_pairwise_align_nucleotide
>>> from iab.algorithms import progressive_msa_and_tree
>>> from functools import partial
>>> gpa = partial(global_pairwise_align_nucleotide, penalize_terminal_gaps=True)
>>> msa, tree = progressive_msa_and_tree(sequences, pairwise_aligner=gpa,
...                                      display_tree=True, display_aln=True)
```

### Phylogenetic reconstruction with UPGMA <link src='73d028'/>

**OLD TEXT**

UPGMA is a heirarchical clustering algorithm. It is widely used, though it's application in phylogenetics is usually restricted to building preliminary trees to "guide" the process of multiple sequence alignment, as it makes some assumptions that don't work well for inferring relationships between organsims. We're going to start with it here however for a few reasons. First, the underlying math is very basic, so we don't need to assume anything about your background. Second, there are some other applications of UPGMA that we'll explore later, including grouping samples based on their species compositions. In general, my strategy with teaching this material is to start by giving you a basic introdution into how a process works so you can visualize and do it. From there, we can get more complex.


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

### Phylogenetic reconstruction with neighbor-joining <link src="JlqeYq"/>

This section is currently a placeholder. You can track progress on this section through [issue #119](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/issues/119). In the meantime, I recommend Chapter 5.2.2 of *[The Phylogenetic Handbook](http://www.amazon.com/gp/product/0521730716/ref=as_li_tl?ie=UTF8&camp=1789&creative=9325&creativeASIN=0521730716&linkCode=as2&tag=anintrotoappl-20&linkId=YLNAKVFX7BV4W5TW")*, by Lemey, Salemi, and Vandamme for discussion of this topic. You can also refer to the [scikit-bio implementation of Neighbor Joining](http://scikit-bio.org/docs/latest/generated/skbio.tree.nj.html), which will be used here (the source code is linked from that page).

### Limitations of distance-based approaches <link src="4hDWma"/>

This section is currently a placeholder. You can track progress on this section through [issue #119](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/issues/119).

## Statistical approaches to phylogenetic reconstruction <link src="VLXTEL"/>

### Bayesian methods <link src="pTnwaZ"/>

This section is currently a placeholder. You can track progress on this section through [issue #119](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/issues/119). In the meantime, I recommend Chapter 7 of *[The Phylogenetic Handbook](http://www.amazon.com/gp/product/0521730716/ref=as_li_tl?ie=UTF8&camp=1789&creative=9325&creativeASIN=0521730716&linkCode=as2&tag=anintrotoappl-20&linkId=YLNAKVFX7BV4W5TW")*, by Lemey, Salemi, and Vandamme for discussion of this topic.

### Maximum likelihood methods <link src="pTnwaZ"/>

This section is currently a placeholder. You can track progress on this section through [issue #119](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/issues/119). In the meantime, I recommend Chapter 6 of *[The Phylogenetic Handbook](http://www.amazon.com/gp/product/0521730716/ref=as_li_tl?ie=UTF8&camp=1789&creative=9325&creativeASIN=0521730716&linkCode=as2&tag=anintrotoappl-20&linkId=YLNAKVFX7BV4W5TW")*, by Lemey, Salemi, and Vandamme for discussion of this topic.

## Rooted versus unrooted trees <link src="HCGmey"/>

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/basic-rooted-tree1.jpg" width=600>

Each *leaf* (or *tip*, or *terminal node*) in this tree represents a sequence, and the length of the horizontal branches between them indicate their dissimilarity to one another. This is a *rooted tree*, which means that it includes an assumption about the last common ancestor of all sequences represented in the tree.

Note that the vertical lines in this tree are used for layout purposes only - they do not represent dissimilarity between sequences.

An **unrooted tree**, like the following, doesn't include an assumption about the last common ancestor of all sequences:

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/basic-unrooted-tree1.jpg" width=600>

## Acknowledgements <link src='99ad11'/>

The material in this section was compiled while consulting the following sources:

1. The Phylogenetic Handbook (Lemey, Salemi, Vandamme)
2. Inferring Phylogeny (Felsenstein)
3. [Richard Edwards's teaching website](http://www.southampton.ac.uk/~re1u06/teaching/upgma/)
