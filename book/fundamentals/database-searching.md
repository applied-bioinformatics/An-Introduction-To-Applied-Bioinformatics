
# Sequence homology searching <link src='d22e6b'/>

In this chapter we'll talk about using pairwise alignment to search databases of biological sequences with the goal of identifying sequence homology. We [previously defined homology](alias://e63a4f) between a pair of sequences to mean that those sequences are derived from a common ancestral sequence. Homology searching is an essential part of making inferences about where a biological sequence came from, and/or what it does. In most cases, if you have an unannotated biological sequence, such as the following protein sequence, it's very hard (really, impossible) to know what it is without more information.

```
> mystery-sequence1
PQITLWQRPLVTIRIGGQLKEALLDTGADDTVLEEMNLPGKWKPKMIGGIGGFIKVRQYDQIPVEIAHKAIGTVLVGPTPVNIIGRNLLTQIGATLNF
```

What a researcher will often do is search this sequence, their *query*, against some *reference database* of annotated sequences to learn what function the sequence performs (if the reference database contains functional annotation of sequences) and/or what organisms are likely to encode this sequence in their genome (if the reference database contains taxonomic annotation of sequences).

Whose genome is the above sequence encoded in? What is its function? Take a minute now to answer these questions using the [Protein BLAST homology search tool on the NCBI website](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp).

In the context of database searching, a query sequence and a reference sequence that we hypothesize to be homologous can be identical to one another, or they can differ as a result of mutation events. When sequences differ, we're often then interested in how much their differ, or their pairwise similarity, which can help us identify the most closely related of several homologs in the reference database. There is an important distinction in the terms *homology* and *similarity*: homology is a discrete variable, and similarity is a continuous variable. A pair of biological sequences either *are* or *are not* derived from a common ancestor, but they can be more or less similar to each other. Saying that two sequences are 80% homologous doesn't make sense. What people generally mean when they say this is that two sequences are 80% similar, and as a result they are hypothesizing homology between the sequences.

## Defining the problem

As mentioned above, if we want to perform a homology search we'll have one or more *query sequences*, and for each we want to know which sequence(s) in a reference database it is most similar to.

Sequence homology searching can be implemented in a few ways. In this chapter, we'll use the local alignment function that we worked with in [the Pairwise Alignment chapter](alias://a76822), ``local_pairwise_align_ssw``, run it many times to search one *query* sequence against many *reference* sequences, and investigate the highest scoring alignment(s) to identify the best database match. Remember that you can always get help with a function by passing it as an argument to ``help``:

```python
>>> from skbio.alignment import local_pairwise_align_ssw
>>> help(local_pairwise_align_ssw)
```

When our reference database starts getting hundreds of millions of bases long (as would be the case if we were searching against 97% OTUs from the [Greengenes small-subunit ribosomal RNA (SSU rRNA) reference database](http://www.ncbi.nlm.nih.gov/pubmed/22134646)), billions of bases long (as would be the case if we were searching against [the human genome](https://genome.ucsc.edu/cgi-bin/hgGateway)) or trillions of bases long (as would be the case if we were searching against the [NCBI non-redundant nucleotide database](http://www.ncbi.nlm.nih.gov/refseq/)), runtime becomes an important consideration. For that reason, learning about *heuristic algorithms* is an essential part of learning about sequence homology searching. Heuristic algorithms apply some rules (i.e., heuristics) to approximate the correct solution to a problem in a fraction of the runtime that would be required if we wanted to be guaranteed to find the correct solution. Heuristic algorithms are very common in bioinformatics, and we'll use them in several other places in this book.

While we'll be aligning nucleotide sequences in this chapter, the same concepts apply to protein homology searching.

## Loading annotated sequences <link src="gAKBxE"/>

The first thing we'll do as we learn about sequence homology searching is load some annotated sequences. The sequences that we're going to work with are derived from the [Greengenes](http://greengenes.secondgenome.com/) [13_8](ftp://greengenes.microbio.me/greengenes_release/gg_13_5/) database, and we're accessing them using the [QIIME default reference project](https://github.com/biocore/qiime-default-reference) project. Greengenes is a database of 16S rRNA gene sequences, a component of the archaeal and bacterial [ribosome](http://www.nature.com/scitable/definition/ribosome-194) (the molecular machine that drives translation of mRNA to proteins). This gene is of a lot of interest to biologists because it's one of about 200 genes that are encoded in the genomes of all known cellular organisms. We'll come back to this gene a few times in the book, notably in [Studying Biological Diversity](alias://2bb2cf). The sequences in Greengenes are taxonomically annotated, meaning that we'll have a collection of gene sequences and the taxonomic identify of the organism whose genome the sequence is found in. If we search an unannotated 16S rRNA query sequence against this database, we can make inferences about what organism our query sequence is from. 

First, let's load Greengenes into a list of ``skbio.DNA`` sequence objects, and associate the taxonomy of each sequence as sequence metadata.

```python
>>> %pylab inline
...
>>> from IPython.core import page
>>> page.page = print
```

```python
>>> import qiime_default_reference as qdr
>>> import skbio
>>> import locale
>>> locale.setlocale(locale.LC_ALL, 'en_US')
...
>>> # Load the taxonomic data
... reference_taxonomy = {}
>>> for e in open(qdr.get_reference_taxonomy()):
...     seq_id, seq_tax = e.strip().split('\t')
...     reference_taxonomy[seq_id] = [t.strip() for t in seq_tax.split(';')]
...
>>> # Load the reference sequences, and associate the taxonmic annotation with
... # each as metadata
... reference_db = []
>>> for e in skbio.io.read(qdr.get_reference_sequences(), format='fasta', constructor=skbio.DNA):
...     seq_tax = reference_taxonomy[e.metadata['id']]
...     e.metadata['domain'] = seq_tax[0][3:] or 'unknown'
...     e.metadata['phylum'] = seq_tax[1][3:] or 'unknown'
...     e.metadata['class'] = seq_tax[2][3:] or 'unknown'
...     e.metadata['order'] = seq_tax[3][3:] or 'unknown'
...     e.metadata['family'] = seq_tax[4][3:] or 'unknown'
...     e.metadata['genus'] = seq_tax[5][3:] or 'unknown'
...     e.metadata['species'] = seq_tax[6][3:] or 'unknown'
...     reference_db.append(e)
...
>>> print("%s sequences were loaded from the reference database." % locale.format("%d", len(reference_db), grouping=True))
```

Next, we'll just inspect a couple of the sequences we loaded. Notice how the specificity of our taxonomic annotations (i.e., how many taxonomic levels are annotated and unknown) differs for different sequences.

```python
>>> reference_db[0]
```

```python
>>> reference_db[-1]
```

For the sake of runtime, we're going to work through this chapter using a random sample of 10,000 sequences from this database. Here we'll use Python's [random module](https://docs.python.org/3/library/random.html) to select sequences at random.

```python
>>> import random
...
>>> reference_db = random.sample(reference_db, k=10000)
>>> print("%s are present in the subsampled database." % locale.format("%d", len(reference_db), grouping=True))
```

We'll also extract some sequences from Greengenes to use as query sequences in our database searches. This time we won't annotate them (to simulate no knowing what organisms they're from). We'll also trim these sequences so they're shorter than the full length references. This will simulate obtaining a partial gene sequence, as is most common with the current sequencing technologies, but will also help to make the examples run faster.

Note that some of our query sequences may also be in our subsampled reference database and some won't. This is realistic: sometimes we're working with sequences that are exact matches to known sequences, and sometimes we're working with sequences that don't match any known sequences (or at least any in the reference database that we're working with).

```python
>>> queries = []
>>> for e in skbio.io.read(qdr.get_reference_sequences(), format='fasta', constructor=skbio.DNA):
...     e = e[100:300]
...     queries.append(e)
>>> queries = random.sample(queries, k=100)
```

Let's also inspect a couple of the query sequences that we'll work with.

```python
>>> queries[0]
```

```python
>>> queries[-1]
```

## Defining a homology search function <link src="0C9FCS"/>

**PICK UP HERE**

Next, we'll define our homology search function. This function will take as input a query sequence and the list of sequences that will serve as our reference database. The results of this will be (Our search function also optionally takes a function to use for performing the alignments, in case for example we wanted to perform global alignments instead of local alignments during our search.)

```python
>>> from iab.algorithms import local_alignment_search
...
>>> %psource local_alignment_search
```

And now we can now perform some database searches. Experiment with different sequences to see how they align by modifying which query sequence you're searching against the database in the next cell and then executing it (remember that you'll need to execute all of the above cells before executing this one).

Also, think about the runtime here. How many sequences are we searching, and how long are they? Does this strategy seem scalable?

```python
>>> import time
...
>>> query_sequence = queries[0]
>>> start_time = time.time()
>>> results = local_alignment_search(query_sequence, reference_db)
>>> stop_time = time.time()
>>> print("Runtime: %1.4f sec" % (stop_time - start_time))
```

```python
>>> for alignment, score, reference_id in results:
...     print('Percent similarity between query and reference: %1.2f%%' % (100 * (1. - alignment[0].distance(alignment[1]))))
...     print('Length of the alignment: %d' % alignment.shape[1])
...     print('Alignment score: %d' % score)
...     print('Reference taxonomy:\n %s' % '; '.join(reference_taxonomy[reference_id]))
...     print('The actual taxonomy of your query is:\n %s' % '; '.join(reference_taxonomy[query_sequence.metadata['id']]))
...     print()
```

In the next cell, I took a shorter exact match from `query1`. What is the effect on our database base search here? How can this happen?

```python
>>> query2 = query1[25:35]
...
>>> start_time = time()
>>> a1, a2, score, ref_id = local_alignment_search(query2, reference_db)
>>> stop_time = time()
...
>>> alignment_length = len(a1)
>>> percent_id = 1 - (hamming(a1, a2)/alignment_length)
...
>>> print(a1)
>>> print(a2)
>>> print(score)
>>> print(ref_id)
>>> print(alignment_length)
>>> print(percent_id)
>>> print("Runtime: %1.4f sec" % (stop_time - start_time))
```

## Using heuristics to reduce runtime for database searches <link src='0f9232'/>

**As illustrated above, runtimes for performing pairwise alignments can be prohibitive for database searching.** The Smith-Waterman implementation we're using here, [SSW](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082138), is very fast. As we covered in the previous chapter, if we were to create an even faster implementation, that would provide some initial help (it'd reduce our runtime by some factor $f$) but ultimately, we're still going to have a problem regardless of how big $f$ is, because **the runtime of the algorithm scales quadratically with sequence lengths**. Experiment with different values of $f$ to see how it changes the curve below.

$f$ represents the fold-reduction in runtime (e.g., $f=2$ represents a two-fold reduction, or halving, of runtime, and $f=10$ equals a ten-fold reduction in runtime).

```python
>>> seq_lengths = range(25)
>>> times = [t * t for t in range(25)]
>>> f = 10000 # no scaling factor
>>> times = [t / f for t in times]
...
>>> plt.plot(range(25), times)
>>> plt.xlabel('Sequence Length')
>>> plt.ylabel('Runtime (s)')
```

Database searching is a slightly different problem however. There are a few different scenarios here:

1. we may have a database that is growing in size (for example, over months and years as more sequences are discovered);
2. we may have a fixed database, but increasingly be obtaining larger numbers of sequences that we want to search against that database;
3. or, the situation that we find ourselves in as of this writing: both.

For the purposes of an exercise, think of the database as one long sequence that we want to align against and a collection of query sequence as another long sequence that we want to search. What do you expect a curve for each of the above to look like?

**The core issue here is that it just takes too long to search each query sequence against the whole database. If we can search against a subset of the database by quickly deciding on certain sequences that are unlikely to match, that may help us to reduce runtime.** A heurtistic in this case would be a rule that we apply to determine which sequecnces we're going to align and which sequences we're not going to align. If we decided to not align against a given reference sequence, it becomes possible that we might miss the best alignment, so we want our heurtistics to be good to make that unlikely. So, when thinking about heurtistics, there are some important considerations:

1. How often do I fail to find the best alignment?
2. Is my runtime reduced enough that I can tolerate not getting the best alignment this often?

Let's look at a few heuristics, starting with a [straw man](https://en.wikipedia.org/wiki/Straw_man). We'll select a random `p` percent of database to align our query against. We'll start by defining `p` as 10%.

```python
>>> from iab.algorithms import approximated_local_alignment_search_random
>>> %psource heuristic_local_alignment_search_random
```

Let's pass our initial `query1` (again, an exact match to a database sequence) and see if we get the right answer, and how much runtime is reduced.

```python
>>> start_time = time()
>>> a1, a2, score, ref_id = heuristic_local_alignment_search_random(query1, reference_db)
>>> stop_time = time()
...
>>> print(a1)
>>> print(a2)
>>> print(score)
>>> print(ref_id)
>>> print("Runtime: %1.4f sec" % (stop_time - start_time))
```

In this case, we know what the right answer is, so we can run this a bunch of times and figure out how often we'll get that right answer. In this case, we'll iterate over the first 100 sequences in the reference database and search 50 bases of each against the full database.

```python
>>> results = []
...
>>> for query_seq in reference_db[:100]:
...     query_seq = query_seq[50:101]
...     a1, a2, score, ref_id = approximated_local_alignment_search_random(query_seq, reference_db)
...     results.append(ref_id == query_seq.metadata['id'])
>>> fraction_correct = results.count(True) / len(results)
...
>>> print("We get the right answer %.2f%% of the time." % (fraction_correct * 100))
```

How much was the run time reduced here? How often were we right? What do you think: good heurtistic?

Let's go with something a little smarter. **We can hypothesize that if the overall composition of a query sequence is different than the overall composition of a reference sequence, it's unlikely that the best alignment will result from that pairwise alignment.** One metric of sequence composition is GC content, so let's use that and only align against sequences whose GC content is within `p` percent of the query sequence.

```python
>>> from iab.algorithms import gc_content, approximated_local_alignment_search_gc
>>> %psource gc_content
```

```python
>>> gc_contents = {}
>>> for seq in reference_db:
...     gc_contents[seq.metadata['id']] = gc_content(seq)
...
>>> sns.set(style="white", palette="muted")
>>> ax = sns.distplot(list(gc_contents.values()))
```

```python
>>> %psource approximated_local_alignment_search_gc
```

If we run our `query1` again, do we get the right answer? How much did we reduce runtime? Do you think this is a better or worse heurtistic?

```python
>>> start_time = time()
>>> a1, a2, score, ref_id = approximated_local_alignment_search_gc(query1, reference_db, gc_contents)
>>> stop_time = time()
...
>>> print(a1)
>>> print(a2)
>>> print(score)
>>> print(ref_id)
>>> print("Runtime: %1.4f sec" % (stop_time - start_time))
```

Let's make the alignment a litte harder by deleting some bases from our query. Now what happens?

```python
>>> query2 = query1[[i < 75 or i > 85 for i in range(len(query1))]]
>>> start_time = time()
>>> a1, a2, score, ref_id = approximated_local_alignment_search_gc(query2, reference_db, gc_contents)
>>> stop_time = time()
...
>>> print(a1)
>>> print(a2)
>>> print(score)
>>> print(ref_id)
>>> print("Runtime: %1.4f sec" % (stop_time - start_time))
```

Again, let's look at how many times we'd be right if we ran this on a subsequence of each reference sequence in the database.

```python
>>> results = []
...
>>> for query_seq in reference_db[:20]:
...     query_seq = query_seq[50:151]
...     a1, a2, score, ref_id = approximated_local_alignment_search_gc(query_seq, reference_db, gc_contents)
...     results.append(ref_id == query_seq.metadata['id'])
>>> fraction_correct = results.count(True) / len(results)
...
>>> print("We get the right answer %.2f%% of the time." % (fraction_correct * 100))
```

So this is looking a little better. It doesn't get our runtime quite as low, but we seem to get the right answer a lot more often. Our evaluation here was pretty limited though: our query sequences are always perfect subsequences of sequences in the reference database. If you were evaluating heuristics for a real-world application, what are some of the other things you might want to test?

**TODO**: Port k-mer composition comparison code from [multiple sequence alignment notebook](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/algorithms/4-multiple-sequence-alignment.ipynb) to look at a third heuristic.

## Is my alignment "good"? Determining whether an alignment is statistically significant. <link src='87c92f'/>

You may have noticed that the score you get back for an alignment isn't extremely informative. It's dependent on the query and reference sequence lengths (and possibly composition, depending on your substitution matrix). An important question then is: **is my alignment score good?**

Remember that an alignment of a pair of sequences represents a hypothesis about homology between those sequences. So when we are asking whether an alignment is good, what we really want to know is: **what fraction of the time would I obtain a score at least this good if my sequences are not homologous?** If this fraction is high, then our alignment is not good. If it's low, then our alignment is good. What is defined as high and low in this context is dependent on **how often you are willing to be wrong**.

If being wrong 5% of the time is acceptable (i.e., you can tolerate a false positive, or calling a pair of sequences homologous when they are actually not, one in twenty times) then you'd set your *cut-off fraction* as 0.05. This fraction is usually call our **alpha**. You have to balance this with how often you can accept false negatives, or deciding that a pair of sequences are not homologous when they actually are. If alpha is high, then you will err on the side of false positives. If alpha is low, then you will err on the side of false negatives. **There is not a hard-and-fast rule for whether false positives or false negatives are better**. It's application specific, so you need to understand your application when making this decision.

In general, when might you prefer to have false positives? When might you prefer to have false negatives?

**In this section, we are going to empiricially determine if a pairwise alignment is better than we would expect by chance.** For each pair of sequences, we're going to align them to determine the score of the alignment, and then we're going to align pairs of sequences that are similar to the query and reference, but that we know are not homologous. We'll do this by *shuffling* or randomizing the order of the bases in the query sequences, and performing another pairwise alignment.

```python
>>> from random import shuffle, choice
>>> from collections import Counter
```

We're going to use python's `random.shuffle` and `random.choice` functions for this. `random.choice` randomly selects an element from a sequence, so we can use it to contruct a random sequence of length `n` as follows:

```python
>>> n = 10
>>> seq = [choice('ACGT') for e in range(n)]
>>> print("".join(seq), Counter(seq))
```

`random.shuffle` randomly re-orders the order of the elements in a sequence, but keeps the composition and length of the sequence the same. Run this next cell a few times to see the sequences that are generated.

```python
>>> shuffle(seq)
>>> print("".join(seq), Counter(seq))
>>> shuffle(seq)
>>> print("".join(seq), Counter(seq))
```

Let's generate a random query sequence. Then we'll generate 99 random variants of that sequence with ``shuffle`` and compute the pairwise alignment for each of those variants against the query sequence. We'll then look at the distribution of those scores.

```python
>>> alphabet = list(skbio.DNA.nondegenerate_chars)
>>> query_seq = skbio.DNA(''.join([choice(alphabet) for e in range(40)]))
```

```python
>>> from iab.algorithms import generate_random_score_distribution
>>> %psource generate_random_score_distribution
```

```python
>>> random_scores = generate_random_score_distribution(query_seq, query_seq)
>>> print(random_scores)
```

```python
>>> import matplotlib.pyplot as plt
...
>>> n, bins, patches = plt.hist(random_scores, facecolor='green', alpha=0.5, bins=range(0,100,1), normed=1)
>>> plt.xlabel('Score')
>>> plt.ylabel('Frequency')
```

Next, we'll compute the score for aligning the query sequence against itself. How does the actual score compare to the random distribution of scores? What does that suggest about our alignment?

```python
>>> _, score, _ = local_pairwise_align_ssw(query_seq, query_seq)
>>> print(score)
...
>>> # plot the distribution of random scores, but add in the actual score
... n, bins, patches = plt.hist(random_scores + [score], facecolor='green', alpha=0.5, bins=range(0,100,1), normed=1)
>>> plt.xlabel('Score')
>>> plt.ylabel('Frequency')
```

**Let's do this experiment again, but this time quanitfy the result by computing the fraction of the random alignments that achieve equal or better scores than the random sequences.**

```python
>>> from iab.algorithms import fraction_better_or_equivalent_alignments
>>> %psource fraction_better_or_equivalent_alignments
```

```python
>>> print(fraction_better_or_equivalent_alignments(query_seq, query_seq))
```

What does this tell us about the quality of our alignment?

Let's now try this for some harder cases, where the query and subject sequences are not identical.

First, let's generate a longer subject sequence at random. Then, we'll create a random query sequence and compare it. Since we're doing this in two random steps, we know that these sequences are not homologous. Does the resulting fraction reflect that?

```python
>>> def random_sequence(moltype, length):
...     result = []
...     alphabet = list(moltype.nondegenerate_chars)
...     for e in range(length):
...         result.append(choice(alphabet))
...     return moltype(''.join(result))
...
>>> subject = random_sequence(skbio.DNA, 250)
>>> query = random_sequence(skbio.DNA, 250)
...
>>> print(query)
>>> print(subject)
```

```python
>>> print(fraction_better_or_equivalent_alignments(query,subject))
```

**We've now looked at two extremes: where sequences are obviously homologous, and where sequences are obviously not homologous. Next, we'll explore the region between these.** We'll obscure the homology of a pair of sequences by randomly introducing some number of substitutions to make them approximately ``percent_id`` equal. By doing this, we can explore how this strategy works for increasingly more distantly related pairs of sequences.

```python
>>> def query_at_percent_id(percent_id, subject):
...     result = []
...     for b in subject:
...         if random() < percent_id:
...             result.append(str(b))
...         else:
...             # choose a base at random that is not the current base
...             # i.e., simulate a substitution event
...             result.append(choice([c for c in subject.nondegenerate_chars if c != b]))
...     return type(subject)(''.join(result))
```

```python
>>> q = query_at_percent_id(0.95,subject)
>>> print(q)
>>> print(subject)
>>> print(fraction_better_or_equivalent_alignments(q,subject))
```

```python
>>> q = query_at_percent_id(0.25,subject)
>>> print(q)
>>> print(subject)
>>> print(fraction_better_or_equivalent_alignments(q,subject))
```

In this case we know that our input sequences are "homologous" because `query` is derived from `subject`. Our method detected that homology when `query` was roughly 95% identical to `subject` (because we got a low fraction) but did not detect that homology when `query` was roughly 25% identicial to `subject`. This gives us an idea of the limit of detection of this method, and is a **real-world problem that biologists face: as sequences are more divergent from one another, detecting homology becomes increasingly difficult.**

**If we want to gain some insight into our limit of detection, we can run a simulation.** If we simulate alignment of different pairs of sequences in steps of different percent identities, we can see where we start failing to observe homology. The following cell illustrates a very simplistic simulation, though this still takes a few minutes to run.

What does this tell us about our limit of detection for homology? What are some things that we might want to do more robustly if we weren't as concerned about runtime?

```python
>>> percent_ids = np.arange(0.0, 1.0, 0.05)
>>> num_trials = 10
>>> results = []
...
>>> for percent_id in percent_ids:
...     p_values = []
...     for i in range(num_trials):
...         subject = random_sequence(skbio.DNA, 250)
...         q = query_at_percent_id(percent_id,subject)
...         p = fraction_better_or_equivalent_alignments(q,subject)
...         p_values.append(p)
...     results.append((percent_id, np.median(p_values), np.mean(p_values)))
>>> pd.DataFrame(results, columns=["Percent id between query and subject",
...                                "Median p-value", "Mean p-value"])
```

```python

```
