
# Pairwise sequence alignment <link src='a76822'/>

One of the most fundamental problems in bioinformatics is determining how "similar" a pair of biological sequences are. There are many applications for this, including inferring the function or source organism of an unknown gene sequence, developing hypotheses about the relatedness of organisms, or grouping sequences from closely related organisms. On the surface this seems like a pretty straight-forward problem, not one that would have been at the center of decades of research and the subject of [one of the most cited papers](http://scholar.google.com/citations?view_op=view_citation&hl=en&user=VRccPlQAAAAJ&citation_for_view=VRccPlQAAAAJ:u-x6o8ySG0sC) in modern biology. In this chapter we'll explore why determining sequence similarity is harder than it might initially seem, and learn about *pairwise sequence alignment*, the standard approach for determining sequence similarity.

Imagine you have three sequences - call them ``r1``and ``r2`` (*r* is for *reference*) and ``q1`` (*q* is for *query*) - and you want to know whether ``q1`` is more similar to ``r1`` or ``r2``. On the surface, it seems like you could just count the number of positions where they differ (i.e., compute the [Hamming distance](http://en.wikipedia.org/wiki/Hamming_distance) between them) to figure this out. Here's what this would look like.

```python
>>> %pylab inline
>>> from __future__ import division, print_function
...
>>> import numpy as np
>>> from IPython.core.display import HTML
>>> from IPython.core import page
>>> page.page = print
```

```python
>>> from scipy.spatial.distance import hamming
>>> from skbio import DNA
...
>>> r1 = DNA("ACCCAGGTTAACGGTGACCAGGTACCAGAAGGGTACCAGGTAGGACACACGGGGATTAA")
>>> r2 = DNA("ACCGAGGTTAACGGTGACCAGGTACCAGAAGGGTACCAGGTAGGAGACACGGCGATTAA")
>>> q1 = DNA("TTCCAGGTAAACGGTGACCAGGTACCAGTTGCGTTTGTTGTAGGAGACACGGGGACCCA")
...
>>> %psource hamming
```

```python
>>> print(hamming(r1, q1))
>>> print(hamming(r2, q1))
```

In this case, ``q1`` has a smaller distance to ``r1`` than it does to ``r2``, so ``q1`` is more similar to ``r1`` than ``r2``. But it's not always that simple.

Here we've assumed that only *substitution events* have occurred, meaning one DNA base was substituted for with another. Let's define ``q2``, which is the same as ``q1`` except that a single base has been deleted at the beginning of the sequence, and a single base has been inserted at the end of the sequence.

```python
>>> q2 = DNA("TCCAGGTAAACGGTGACCAGGTACCAGTTGCGTTTGTTGTAGGAGACACGGGGACCCAT")
>>> print(hamming(r1, q2))
```

This change had a big effect on the distance between the two sequences. In this case, the deletion event at the beginning of ``q2`` has shifted that sequence relative to ``r1``, which resulted in many of the bases "downstream" of the deleted base being different. However the sequences do still seem fairly similar, so perhaps this relatively large distance isn't biologically justified.

What we'd really want to do is have a way to indicate that a deletion seems to have occurred in ``q2``. Let's define ``q3``, where we use a ``-`` character to indicate a deletion with respect to ``r1``. This results in what seems like a more reasonable distance between the two sequences:

```python
>>> q3 = DNA("-TCCAGGTAAACGGTGACCAGGTACCAGTTGCGTTTGTTGTAGGAGACACGGGGACCCA")
>>> print(hamming(r1, q3))
```

What we've done here is create a pairwise alignment of ``r1`` and ``q3``. In other words, we've **aligned** the positions that minimize the similarity of the two sequences, using the ``-`` to fill in spaces where one character is missing with respect to the other sequence. We refer to ``-`` characters in aligned sequences as **gap characters**, or gaps.

The *alignment* is of these two sequences is clear if we print them  out, one on top of the other:

```python
>>> print(r1)
>>> print(q3)
```

Scanning through these two sequences, we can see that they are largely identical, with the exception of one ``-`` character, and about 25% *substitutions* of one base for another.

## What is a sequence alignment? <link src='e63a4f'/>

Let's take a minute to think about sequence evolution and what a biological sequence alignment actually is. Over the course of evolution, the sequence will likely change, most frequently due to random errors in replication (or the copying of a DNA sequence), or **mutations**. Some of the types of mutation events that can occur are:

* **substitutions**, where one base (or amino acid, in protein sequences) is replaced with another;
* **insertions**, where one or more contiguous bases are inserted into a sequence;
* and **deletions**, where one or more contiguous bases are deleted from a sequence.

(Other types of mutation events can occur, but we're going to focus on these for now.)

Figure 1 illustrates how one ancestral DNA sequence (Figure 1a), over time, might evolve into two derived sequences (Figure 1b). When two or more sequences are derived from a single ancestral sequence, as is the case in this example, those sequences are said to be **homologs** of one another, or homologous sequences. On a piece of paper, make a hypothesis about which of these types of mutation events occurred where over our hypothetical evolution of these sequences.

<figure>
    <img src="images/alignment.png">
    <figcaption>Figure 1: Sequence evolution and pairwise sequence alignment.</figcaption>
</figure>
<p>

**The goal of pairwise sequence alignment is, given two sequences, to generate a hypothesis about which sequence positions derived from a common ancestral sequence position.** In practice, we develop this hypothesis by aligning the sequences to one another inserting gaps as necessary, in a way that maximizes their similarity. This is a **maximum parsimony** approach, where we assume that the simplest explanation (the one involving the fewest or least extreme mutation events) is the most likely.

In nearly all cases, the only sequences we have to work with are the modern (derived) sequences, as illustrated in Figure 1c. The ancestral sequence is not something we have access to (for example, because the organism whose genome it was present in went extinct 100 million years ago).

Figure 1d illustrates one possible alignment of these two sequences. Just as the notes you made about which types of mutation events may have happened where represents your *hypothesis* about the evolutionary events that took place, a sequence alignment that you might get from a computer program such as BLAST is also only a hypothesis. You can think of an alignment as a table, where the rows are sequences and the columns are positions in those sequences. When you have two or more aligned sequences, there will, by definition, always be the same number of columns in each row. Each column in your alignment represents a hypothesis about the evolutionary events that occurred at that position since the last ancestor of the aligned sequences (the sequence in Figure 1a in our example).

One thing that's worth pointing out at this point is that because we don't know what the ancestral sequence was, when we encounter a gap in a pairwise alignment, we generally won't know whether a deletion occurred in one sequence, or an insertion occurred in the other. For that reason, you will often see the term **indel** used to refer to these either an insertion or deletion events.

In the next section we'll work through our first bioinformatics algorithm, in this case a very simple (and also simplistic) method for aligning a pair of sequences. As you work through this exercise, think about why it might be too simple given what you know about biological sequences.

## A simple procedure for aligning a pair of sequences <link src='86c6b7'/>

Lets define two sequences, ``seq1`` and ``seq2``, and develop an approach for aligning them.

```python
>>> seq1 = "ACCGGTGGAACCGGTAACACCCAC"
>>> seq2 = "ACCGGTAACCGGTTAACACCCAC"
```

I'm going to use a function in the following cells called ``show_table`` to display a table that we're going to use to develop our alignment. Once a function has been imported, you can view the source code for that function. This will be useful as we begin to explore some of the algorithms that are in use throughout these notebooks, and you should spend time reading the source code until you're sure that you understand what's happening.

Here's how you'd import a function and then view its source code:

```python
>>> from iab.algorithms import show_table
```

```python
>>> %psource show_table
```

Now let's look at how to align these sequences.

**Step 1.** Create a matrix (or *array*), where the columns represent the positions in ``seq1`` and the rows represent the positions in ``seq2``. We'll initialize this matrix with zeros.

```python
>>> num_rows = len(seq2)
>>> num_cols = len(seq1)
>>> data = np.zeros(shape=(num_rows, num_cols), dtype=np.int)
...
>>> HTML(show_table(seq1, seq2, data))
```

**Step 2.** Score the cells so if the characters at the corresponding row and column are the same the value is changed from zero to one. We can then revew the resulting matrix. For clarity, we'll have ``show_table`` hide the zero values.

```python
>>> for row_number, row_character in enumerate(seq2):
...     for col_number, col_character in enumerate(seq1):
...         if row_character == col_character:
...             data[row_number, col_number] = 1
...
>>> HTML(show_table(seq1, seq2, data, hide_zeros=True))
```

**Step 3**: Identify the longest diagonal stretches of non-zero characters (we'll call these *diagonals*). Diagonals indicate segments of the two sequences that are identical and uninterrupted by mismatched characters (substitution events) or indel events.

We can identify the longest diagonals as follows:

```python
>>> # create a copy of our data matrix to work with, so we
... # leave the original untouched.
... summed_data = data.copy()
>>> # iterate over the cells in our data matrix, starting in
... # the second row and second column
... for i in range(1, summed_data.shape[0]):
...     for j in range(1, summed_data.shape[1]):
...         # if the value in the current cell is greater than zero
...         # (i.e., the characters at the corresponding pair of
...         # sequence positions are the same), add the value from the
...         # cell that is diagonally up and to the left.
...         if summed_data[i, j] > 0:
...             summed_data[i, j] += summed_data[i-1][j-1]
...
>>> # Identify the longest diagonal
... print("The longest diagonal is %d characters long." % summed_data.max())
>>> HTML(show_table(seq1, seq2, summed_data, hide_zeros=True))
```

**Step 4**: Next, we'd want to transcribe some of the possible alignments that arise from this process.

We're going to gloss over how to do this algorithmically for the moment, as we'll come back to that in a lot of detail later in this chapter. Briefly, what we want to do is start with the longest diagonal and work bacwards to transcribe the alignment by writing down the matching characters. When we encounter a break in the diagonal, we find the next highest score diagonal that starts in a cell that is up and/or to the left of the cell when the previous diagonal you were following ends. For every cell that you move upwards, you'd insert a gap in the sequence on the horizontal axis of your matrix. For every cell that you move leftwards, you'd insert a gap in the sequence on the vertical axis of your matrix.

If this is confusing, don't worry about it for the moment. We'll be back to this in a lot more detail soon.

Here are two possible alignments:

Alignment 1:
```
ACCGGTGGAACCGG-TAACACCCAC
ACCGGT--AACCGGTTAACACCCAC
```

Alignment 2:
```
ACCGGTGGAACCGGTAACACCCAC
ACCGGT--------TAACACCCAC
```

**Remember that an alignment represents a hypothesis about the evolutionary history of a sequence.  Which of these hypotheses do you think is more likely to be true based on what you know about sequence evolution?** Why might the first alignment be the more biologically relevant one? Why might the second?

**As an exercise**, go back to where we defined `seq1` and `seq2` and re-define one or both of those as other sequences. Execute the code through here and see how the matrices change.

### Why this simple procedure is too simplistic <link src="jzshiO"/>

I suggested above that you keep a list of assumptions that are made by this approach. Here are a couple of the very problematic ones.

1. We're scoring all matches as 1 and all mismatches as 0. This suggests that all matches are equally likely, and all mismatches are equally unlikely. What's a more biologically meaningful way to do this (think about protein sequences here)?
2. Similarly, every gap that is introduced results in the same penalty being incurred. Based on what we know about how insertion/deletion events occur, it likely makes more sense to score *opening a new gap* differently from *extending an existing gap*.

All scoring schemes have limitations, and you should consider alignments that you get (e.g., from systems such as [BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi)) as hypotheses. You'll need to determine if you agree with the result that a computational system gives you. Algorithms such as the one we just explored are there to help you do your work, but their answers are based on models (for example, how we model matches, mismatches and gaps here) and as we're learning here, the models are not perfect. Be skeptical!

Another important consideration as we think about algorithms for aligning pairs of sequences is how long this will take to run (or in technical terminology, the [computational complexity](http://bigocheatsheet.com/) of the algorithm). When searching a novel sequence against a database, you may have billions of bases to search against (which would correspond to billions of columns in these matrices). How can this be done efficiently?

Over the next several sections we'll explore ways of addressing the two issues noted above. We'll introduce the problem of the computational complexity at the end of this chapter, and explore approaches for addressing that (i.e., making database searching faster) in the next chapter.

## Substitution matrices <link src='9f5e71'/>

**TODO: pick up here! I started editing this section before realizing that I needed to update some stuff above.**

The first of the limitations we identified above was that all matches and mismatches were scored equally when aligning a pair of sequences. To understand why this is a problem, let's think about the meaning of match and mismatches.

When we align a pair of

 though we know that that isn't the most biologically meaningful way to score an alignment. This is particularly the case for protein sequences, where each amino acid residue has different chemical properties, that impact the structure and function of the protein.

We'll next explore a more general approach to the problem of *global sequence alignment* for protein sequences, or aligning a pair of protein sequences from beginning to end. We'll start by defining a **substitution matrix which defines the score associated with substitution of one amino acid for another**.

Early work on defining protein substitution matrices was performed by Dayhoff in the 1970s and by Henikoff and Henikoff in the early 1990s. We'll start by working with a substitution matrix known as the blosum 50 matrix, which was [presented in PNAS in 1992](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC50453/). There are many different criteria that can be used to define how likely the substitution of Briefly, these matrices are defined empirically, by aligning sequences manually or through automated systems, and counting how frequent certain substitutions are. There is a good [wikipedia article on this topic](http://en.wikipedia.org/wiki/BLOSUM).

We can import this from the ``iab`` support module, and then look up scores for substitutions of amino acids. Based on these scores and the biochemistry of the amino acids ([see the molecular structures on Wikipedia](http://en.wikipedia.org/wiki/Amino_acid)), does a positive score represent a more or less favorable substitution?

```python
>>> from iab.algorithms import blosum50
...
>>> print(blosum50['A']['G'])
>>> print(blosum50['G']['A'])
...
>>> print(blosum50['W']['K'])
...
>>> print(blosum50['A']['A'])
>>> print(blosum50['W']['W'])
```

Here's a global view of the matrix.

```python
>>> aas = list(blosum50.keys())
>>> aas.sort()
>>> data = []
>>> for aa1 in aas:
...     row = []
...     for aa2 in aas:
...         row.append(blosum50[aa1][aa2])
...     data.append(row)
...
>>> HTML(show_table(aas, aas, data))
```

## Needleman-Wunsch global pairwise sequence alignment <link src='15efc2'/>

Now let's get started on using this to align a pair of sequences.

**Step 1.** Create a matrix, where the columns represent the positions in ``seq1`` and the rows represent the positions in ``seq2``.

```python
>>> ## Example adapted from Biological Sequence Analysis Chapter 2.
...
... seq1 = "HEAGAWGHEE"
>>> seq2 = "PAWHEAE"
```

```python
>>> data = []
>>> for p in seq2:
...     data.append(['-']*len(seq1))
...
>>> HTML(show_table(seq1, seq2, data))
```

**Step 2**:  Using a substitution matrix, score each cell in the matrix.

```python
>>> from iab.algorithms import generate_score_matrix
...
>>> %psource generate_score_matrix
```

```python
>>> score_matrix = generate_score_matrix(seq1,seq2,blosum50)
...
>>> HTML(show_table(seq1,
...                     seq2,
...                     score_matrix))
```

**Step 3**: Generate the dynamic programming and traceback matrices.

In the next step we determine the best alignment given the sequences and scoring scheme in what we'll call the **dynamic programming matrix**, and then define programmatically how to transcribe the alignment in what we'll call the **traceback matrix** to yield a pair of aligned sequences.

For the convenience of coding this algorithm, it helps to define the dynamic programming matrix with one extra row and one extra column relative to the score matrix, and make these the first column and row of the matrix. These then represent the beginning of the alignment position $(0, 0)$. The score $F$ for cell $(i, j)$, where $i$ represents the row number and $j$ represents the column number, is defined for the first row and column as follows:

$$
\begin{align}
& F(0, 0) = 0\\
& F(i, 0) = F(i-1, 0) - d\\
& F(0, j) = F(0, j-1) - d\\
\end{align}
$$

This matrix, pre-initialization, would look like the following. As an exercise, try computing the values for the cells in the first four rows in column zero, and the first four columns in row zero. As you fill in the value for a cell, for all cells with a score based on another score in the matrix (i.e., everything except for $F(0, 0)$), draw an arrow from that cell to the cell whose score it depends on.

For the sake of this exercise, define the gap penalty, $d$, as $d=8$.

```python
>>> data = []
>>> # This is a hack: to pad the matrix with an
... # extra row and column at the beginning I'm just prepending a
... # space to each sequence. Need to improve handling of that.
... padded_seq1 = " " + seq1
>>> padded_seq2 = " " + seq2
...
>>> for p in padded_seq2:
...     data.append(['-']*len(padded_seq1))
...
>>> HTML(show_table(padded_seq1, padded_seq2, data))
```

Initializing this would result in the following.

```python
>>> # We'll define the gap penalty as 8.
... d = 8
...
>>> data[0][0] = 0
>>> for i in range(1,len(padded_seq2)):
...     data[i][0] = data[i-1][0] - d
...
>>> for j in range(1,len(padded_seq1)):
...     data[0][j] = data[0][j-1] - d
...
>>> HTML(show_table(padded_seq1, padded_seq2, data, cell_width=4))
```

Next, we'll compute the scores for all of the other cells in the matrix, starting at position $(1, 1)$.

In a Needleman-Wunsch alignment, the score $F$ for cell $(i, j)$ (where $i$ is the row number and $j$ is the column number, and $i > 0$ and $j > 0$) is computed as the maximum of three possible values.

$$
F(i, j) = max \left(\begin{align}
& F(i-1, j-1) + s(c_i, c_j)\\
& F(i-1, j) - d\\
& F(i, j-1) - d
\end{align}\right)
$$

In this notation, $s$ refers to the substitution matrix, $c_i$ and $c_j$ refers to characters in `seq1` and `seq2`, and $d$ again is the gap penalty. Describing the scoring function in English, we score a cell with the maximum of three values: either the value of the cell up and to the left plus the score for the substitution taking place in the current cell (which you find by looking up the substitution in the substitution matrix); the value of the cell above minus the gap penalty; or the value of the cell to the left minus the gap penalty. In this way, you're determining whether the best (highest) score is obtained by inserting a gap in sequence 1 (corresponding to $F(i-1, j) - d$), inserting a gap in sequence 2 (corresponding to $F(i, j-1) - d$), or aligning the characters in sequence 1 and sequence 2 (corresponding to $F(i-1, j-1) + s(c_i, c_j)$).

As an exercise, fill in the values of cells $(1, 1)$, $(1, 2)$, and $(2, 1)$ by hand. Remember to insert arrows indicating which cell each score was derived from as you fill in the matrix, and notice the situation that you encounter when computing the value for $(2, 1)$. Which arrow do you draw there? Keep this issue in mind, and think about how it might affect your final result.

The function in the next cell generates the dynamic programming and traceback matrices for us. You should review this code to understand exactly how it's working.

```python
>>> from iab.algorithms import format_dynamic_programming_matrix, format_traceback_matrix
>>> from skbio.alignment._pairwise import _compute_score_and_traceback_matrices
...
>>> %psource _compute_score_and_traceback_matrices
```

You can now apply this function to `seq1` and `seq2` to compute the dynamic programming and traceback matrices. Based on the arrows in your traceback matrix, what do you think the four different values used in this traceback matrix represent?

```python
>>> from skbio.sequence import Protein
>>> from skbio.alignment import TabularMSA
...
>>> seq1 = TabularMSA([Protein("HEAGAWGHEE")])
>>> seq2 = TabularMSA([Protein("PAWHEAE")])
...
>>> nw_matrix, traceback_matrix = _compute_score_and_traceback_matrices(
...     seq1, seq2, 8, 8, blosum50)
...
>>> print(format_dynamic_programming_matrix(seq1, seq2, nw_matrix))
```

```python
>>> print(format_traceback_matrix(seq1, seq2, traceback_matrix))
```

**Step 4**: Transcribe the alignment.

We can now read the dynamic programming and traceback matrices to transcribe and score the alignment of sequences 1 and 2. To do this, we start at the bottom right of the matrices and traceback the arrows to cell $(0, 0)$.

* Every time we hit a vertical arrow (represented by `|`), we consume a character from sequence 2 (the vertical sequence) and add a gap to sequence 1.
* Every time we hit a horizontal arrow (represented by `-`), we consume a character from sequence 1 (the horizontal sequence) and add a gap to sequence 2.
* Every time we hit a diagonal arrow (represented by `\`), we consume a character from sequence 1 and sequence 2.

As you transcribe the alignment, write sequence 1 on top of sequence 2, and work from right to left (since you are working backwards through the matrix).

The score in the bottom right cell of the matrix is the score for the alignment.

Work through this process on paper, and then review the function in the next cell to see how this looks in python.

```python
>>> from skbio.alignment._pairwise import _traceback
>>> %psource _traceback
```

You can then execute this as follows, and print out the resulting alignment.

```python
>>> aln1, aln2, score, _, _ = _traceback(traceback_matrix,nw_matrix,seq1,seq2, nw_matrix.shape[0]-1, nw_matrix.shape[1]-1)
...
>>> print(aln1[0])
>>> print(aln2[0])
>>> print(score)
```

**Next steps**: All of those steps are a bit ugly, so as a developer you'd want to make this functionality generally accessible to users. To do that, you'd want to define a function that takes all of the necessary input and provides the aligned sequences and the score as output, without requiring the user to make several function calls. What are the required inputs? What steps would this function need to perform?

```python
>>> from skbio.alignment import global_pairwise_align
>>> %psource global_pairwise_align
```

```python
>>> aln, score, _ = global_pairwise_align(Protein("HEAGAWGHEE"), Protein("PAWHEAE"), 8, 8, blosum50, penalize_terminal_gaps=True)
...
>>> print(aln)
>>> print(score)
```

## Global versus local alignment <link src='c80f21'/>

The alignment we just constructed is what's known as a *global alignment*, meaning we align both sequences from the beginning through the end. This has some important specific applications: for example, if we have two full-length protein sequences, and we have a crystal structure for one of them, we may want to have a direct mapping between all positions in both sequences. Probably more common however, is local alignment: where we have a pair of sequences that we suspect may partially overlap each other, and we want to know what the best possible alignment of all or part of one sequence is with all or part of the other sequences. Perhaps the most widely used application of this is in the BLAST search algorithm, where we have a query sequence and we want to find the closest match (or matches) in a reference database containing many different gene sequences. In this case, the whole reference database could be represented as a single sequence, as we could perform a local alignment against it to find the region that contains the highest scoring match.

## Smith-Waterman local sequence alignment <link src='c9656e'/>

The Smith-Waterman algorithm is used for performing pairwise local alignment. It is nearly identical to Needleman-Wunsch, with three small important differences.

First, initialization is easier:
$$
\begin{align}
& F(0, 0) = 0\\
& F(i, 0) = 0\\
& F(0, j) = 0
\end{align}
$$

```python
>>> data = []
>>> # This is a hack: to pad the matrix with an
... # extra row and column at the beginning I'm just prepending a
... # space to each sequence. Need to improve handling of that.
... padded_seq1 = " " + str(seq1[0])
>>> padded_seq2 = " " + str(seq2[0])
...
>>> for p in padded_seq2:
...     data.append(['-']*len(padded_seq1))
...
>>> HTML(show_table(padded_seq1, padded_seq2, data))
```

```python
>>> data[0][0] = 0
>>> for i in range(1,len(padded_seq2)):
...     data[i][0] = 0
...
>>> for j in range(1,len(padded_seq1)):
...     data[0][j] = 0
...
>>> HTML(show_table(padded_seq1, padded_seq2, data))
```

Next, there is one additional term in the scoring function:

$$
F(i, j) = max \left(\begin{align}
& 0\\
& F(i-1, j-1) + s(c_i, c_j)\\
& F(i-1, j) - d\\
& F(i, j-1) - d)
\end{align}\right)
$$

```python
>>> %psource _compute_score_and_traceback_matrices
```

```python
>>> from skbio.alignment._pairwise import _init_matrices_sw
...
>>> sw_matrix, traceback_matrix = \
...     _compute_score_and_traceback_matrices(seq1, seq2, 8, 8, blosum50, new_alignment_score=0.0, init_matrices_f=_init_matrices_sw)
...
>>> print(format_dynamic_programming_matrix(seq1, seq2, sw_matrix, 7))
```

```python
>>> print(format_traceback_matrix(seq1, seq2, traceback_matrix, cell_width=5))
```

And finally, during the traceback step, you begin in the cell with the highest value, rather than the bottom right cell of the matrix.

```python
>>> %psource _traceback
```

```python
>>> max_value = 0.0
>>> max_i = 0
>>> max_j = 0
>>> for i in range(sw_matrix.shape[0]):
...     for j in range(sw_matrix.shape[1]):
...         if sw_matrix[i, j] > max_value:
...             max_i, max_j = i, j
...             max_value = sw_matrix[i, j]
...
>>> aln1, aln2, score, start_a1, start_a2 = _traceback(traceback_matrix, sw_matrix, seq1, seq2, max_i, max_j)
>>> print(aln1[0])
>>> print(aln2[0])
>>> print(score)
```

Again, we can define a *convenience function*, which will allow us to provide the required input and just get our aligned sequences back.

```python
>>> from skbio.alignment import local_pairwise_align
...
>>> %psource local_pairwise_align
```

And we can take the *convenience function* one step further, and wrap `sw_align` and `nw_align` up in a more general `align` function, which takes a boolean parameter (i.e., `True` or `False`) indicating where we want a local or global alignment.

```python
>>> def align(sequence1, sequence2, gap_penalty, substitution_matrix, local):
...     if local:
...         return local_pairwise_align(sequence1, sequence2, gap_penalty, gap_penalty, substitution_matrix)
...     else:
...         return global_pairwise_align(sequence1, sequence2, gap_penalty, gap_penalty, substitution_matrix)
```

```python
>>> aln, score, _ = align(Protein('HEAGAWGHEE'), Protein('PAWHEAE'), 8, blosum50, True)
...
>>> print(aln)
>>> print(score)
```

```python
>>> aln, score, _ = align(Protein('HEAGAWGHEE'), Protein('PAWHEAE'), 8, blosum50, False)
...
>>> print(aln)
>>> print(score)
```

So there you have it: the basics of pairwise sequence alignment, which is easily the most fundamental algorithm in bioinformatics.

### Smith-Waterman local alignment with affine gap scoring <link src='976169'/>

The second limitation of the our simple alignment algorithm, and one that is also present in our version of Smith-Waterman as implemented above, is that all gaps are scored equally whether they represent the opening of a new insertion/deletion, or the extension of an existing insertion/deletion. This isn't ideal based on what we know about how insertion/deletion events occur (see [this discussion of replication slippage](http://www.ncbi.nlm.nih.gov/books/NBK21114/)). Instead, **we might want to incur a large penalty for opening a gap, but a smaller penalty for extending a gap**. To do this, **we need to make two small changes to our scoring scheme**. When we compute the score for a gap, we should incur a *gap open penalty* if the previous max score was derived from inserting a gap character in the same sequence. If we represent our traceback matrix as $T$, our gap open penalty as $d^0$, and our gap extend penalty as $d^e$, our scoring scheme would look like the following:

$$
F(i, j) = max \left(\begin{align}
& 0\\
& F(i-1, j-1) + s(c_i, c_j)\\
& \left\{\begin{array}{l l} F(i-1, j) - d^e \quad \text{if $T(i-1, j)$ is gap}\\ F(i-1, j) - d^o \quad \text{if $T(i-1, j)$ is not gap} \end{array}  \right\} \\
& \left\{\begin{array}{l l} F(i, j-1) - d^e \quad \text{if $T(i, j-1)$ is gap}\\ F(i, j-1) - d^o \quad \text{if $T(i, j-1)$ is not gap} \end{array}  \right\}
 \end{align}\right)
$$

Notice how we only use the gap extend penalty if the previous max score resulted from a gap in the same sequence (which we know by looking in the traceback matrix) because it represents the continuation of an existing gap in that sequence. This is why we check for a specific type of gap in $T$, rather than checking whether $T$ `!= '\'`.

Here is our ``_compute_score_and_traceback_matrices`` function again for reference.

```python
>>> %psource _compute_score_and_traceback_matrices
```

Take a look at how the scores differ with these additions.

```python
>>> seq1 = TabularMSA([Protein("HEAGAWGHEE")])
>>> seq2 = TabularMSA([Protein("PAWHEAE")])
...
>>> sw_matrix, traceback_matrix = _compute_score_and_traceback_matrices(seq1, seq2, 8, 1, blosum50)
...
>>> print(format_dynamic_programming_matrix(seq1, seq2, sw_matrix))
```

```python
>>> print(format_traceback_matrix(seq1, seq2, traceback_matrix))
```

The convenience functions we worked with above all take ``gap_open_penalty`` and ``gap_extend_penalty``, so we can use those to explore sequence alignment with affine gap scoring. Here I define `seq1` to be slightly different than what I have above. Notice how we get different alignments when we use affine gap penalties (i.e., ``gap_extend_penalty`` is not equal to ``gap_open_penalty``) versus equal gap open and gap extend penalties.

```python
>>> help(local_pairwise_align)
```

```python
>>> seq1 = TabularMSA([Protein("HEAGAWGFHEE")])
>>> seq2 = TabularMSA([Protein("PAWHEAE")])
```

```python
>>> aln, score, _ = global_pairwise_align(seq1, seq2, 8, 8, blosum50)
>>> print(aln)
>>> print(score)
```

```python
>>> aln, score, _ = global_pairwise_align(seq1, seq2, 8, 1, blosum50)
>>> print(aln)
>>> print(score)
```

## How long does pairwise sequence alignment take? <link src='ac446d'/>

The focus of this book is *applied* bioinformatics, and **some of the practical considerations we need to think about when developing applications is their runtime and memory requirements**. The third issue we mentioned above is general to the problem of sequence alignment: runtime can be problematic. Over the next few cells we'll explore the runtime of sequence alignment.

We just worked through a few algorithms for pairwise sequence alignment, and used some toy examples with short sequences. What if we wanted to scale this up to align much longer sequences, or to align relatively short sequences against a large database?

To explore runtime, let's use the IPython [magic function](http://ipython.org/ipython-doc/dev/interactive/tutorial.html#magic-functions) called ``timeit``, which runs a given command many times and reports the average time it takes to fun. We'll use this to see how long global alignment takes to run. Note that we don't care about getting the actual alignment back anymore, we just want the runtime in seconds.

```python
>>> %timeit global_pairwise_align(seq1, seq2, 8, 1, blosum50)
```

Next, let's apply this to pairs of sequences where we vary the length. We don't really care what the sequences are here, so we'll use python's ``random`` module to get random pairs of sequences. Let's play with that first to see how it can be applied to generate random sequences, as that's generally useful functionality.

```python
>>> from random import choice
...
>>> def random_sequence(moltype, length):
...     result = []
...     alphabet = list(moltype.nondegenerate_chars)
...     for e in range(length):
...         result.append(choice(alphabet))
...     return moltype(''.join(result))
```

```python
>>> print(random_sequence(DNA, 10))
>>> print(random_sequence(DNA, 10))
>>> print(random_sequence(DNA, 25))
>>> print(random_sequence(DNA, 50))
```

Next, let's define a loop where we align, randomly, pairs of sequences of increasing length, and compile the time it took to align the sequences. Here we're going to use a faster version of pairwise alignment that's implemented in scikit-bio, to facilitate testing with more alignments.

```python
>>> import timeit
>>> from skbio.alignment import local_pairwise_align_ssw
...
>>> times = []
>>> seq_lengths = range(50,100000,20000)
...
>>> def get_time_function(seq_length):
...     def f():
...         seq1 = random_sequence(DNA, seq_length)
...         seq2 = random_sequence(DNA, seq_length)
...         local_pairwise_align_ssw(seq1, seq2)
...     return f
...
>>> for seq_length in seq_lengths:
...     times.append(min(timeit.Timer(get_time_function(seq_length)).repeat(repeat=3, number=3)))
```

If we look at the run times, we can see that they are increasing with increasing sequence lengths.

```python
>>> for seq_length, t in zip(seq_lengths, times):
...     print("%d\t%1.4f sec" % (seq_length, t))
```

That's expected, but what we care about is how they're increasing. Can we use this information to project how well this alignment would work if our sequences were much longer? This is where plotting becomes useful.

```python
>>> import matplotlib.pyplot as plt
...
>>> plt.plot(seq_lengths, times)
>>> plt.xlabel('Sequence Length')
>>> plt.ylabel('Runtime (s)')
```

**One good question is whether developing a version of this algorithm which can run in parallel would be an effective way to make it scale to larger data sets.** In the next cell, we look and how the plot would change if we could run the alignment process over four processors. This would effectively make each alignment run four times as fast (so each runtime would be divided by four) but it doesn't solve our scalability problem.

```python
>>> # if we could split this process over more processors (four, for example)
... # that would effectively reduce the runtime by 1/4
... times = [t / 4 for t in times]
...
...
>>> plt.plot(seq_lengths, times)
>>> plt.xlabel('Sequence Length')
>>> plt.ylabel('Runtime (s)')
```

**Notice that the runtimes in the plot (the y-axis) are smaller, but shape of the curve is the same.** This tells us that we won't be in trouble as soon (we can run bigger alignments in a reasonable amount of time), but we'll still be in trouble eventually. While parallelization can help with this class of computational problem -- one that scales [quadratically](http://en.wikipedia.org/wiki/Quadratic_time) -- it doesn't resolve the problem completely.

How an algorithm scales with input size is referred to as its computational complexity. You can explore the computational complexity of different types of algorithms in the [Big-O cheatsheet](http://bigocheatsheet.com/), though it's a fairly advanced introduction to the topic (and one that's usually covered in the second or third year for Computer Science majors).

In the next chapter we'll begin exploring ways to address this scalability issue by approximating solutions to the problem.
