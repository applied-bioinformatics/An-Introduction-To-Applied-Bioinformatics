
# Sequence homology searching <link src='d22e6b'/>

In this chapter we'll talk about using pairwise alignment to search databases of biological sequences with the goal of identifying sequence homology. We [previously defined homology](alias://e63a4f) between a pair of sequences to mean that those sequences are derived from a common ancestral sequence. Homology searching is an essential part of making inferences about where a biological sequence came from, and/or what it does. In most cases, if you have an unannotated biological sequence, such as the following protein sequence, it's very hard (really, impossible) to know what it is without more information.

```
>mystery-sequence1
PQITLWQRPLVTIRIGGQLKEALLDTGADDTVLEEMNLPGKWKPKMIGGIGGFIKVRQYDQIPVEIAHKAIGTVLVGPTPVNIIGRNLLTQIGATLNF
```

What a researcher will often do is search this sequence, their *query*, against some *reference database* of annotated sequences to learn what function the sequence performs (if the reference database contains functional annotation of sequences) and/or what organisms are likely to encode this sequence in their genome (if the reference database contains taxonomic annotation of sequences).

Whose genome is the above sequence encoded in? What is its function? Take a minute now to answer these questions using the [Protein BLAST homology search tool on the NCBI website](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome).

In the context of database searching, a query sequence and a reference sequence that we hypothesize to be homologous can be identical to one another, or they can differ as a result of mutation events. When sequences differ, we're often then interested in how much they differ, or their pairwise similarity, which can help us identify the most closely related of several homologs in the reference database. There is an important distinction in the terms *homology* and *similarity*: homology is a discrete variable, and similarity is a continuous variable. A pair of biological sequences either *are* or *are not* derived from a common ancestor, but they can be more or less similar to each other. Saying that two sequences are 80% homologous doesn't make sense. What people generally mean when they say this is that two sequences are 80% similar, and as a result they are hypothesizing homology between the sequences.

## Defining the problem <link src="P7kHdy"/>

As mentioned above, if we want to perform a homology search we'll have one or more *query sequences*, and for each we want to know which sequence(s) in a reference database it is most similar to.

Sequence homology searching can be implemented in a few ways. In this chapter, we'll use the local alignment function that we worked with in [the Pairwise Alignment chapter](alias://a76822), ``local_pairwise_align_ssw``, run it many times to search one *query* sequence against many *reference* sequences, and investigate the highest scoring alignment(s) to identify the best database match. Remember that you can always get help with a function by passing it as an argument to ``help``:

```python
>>> from skbio.alignment import local_pairwise_align_ssw
>>> help(local_pairwise_align_ssw)
```

When our reference database starts getting hundreds of millions of bases long (as would be the case if we were searching against 97% OTUs from the [Greengenes small-subunit ribosomal RNA (SSU rRNA) reference database](http://www.ncbi.nlm.nih.gov/pubmed/22134646)), billions of bases long (as would be the case if we were searching against [the human genome](https://genome.ucsc.edu/cgi-bin/hgGateway)) or trillions of bases long (as would be the case if we were searching against the [NCBI non-redundant nucleotide database](http://www.ncbi.nlm.nih.gov/refseq/)), runtime becomes an important consideration. For that reason, learning about *heuristic algorithms* is an essential part of learning about sequence homology searching. Heuristic algorithms apply some rules (i.e., heuristics) to approximate the correct solution to a problem in a fraction of the runtime that would be required if we wanted to be guaranteed to find the correct solution. Heuristic algorithms are very common in bioinformatics, and we'll use them in several other places in this book.

While we'll be aligning nucleotide sequences in this chapter, the same concepts apply to protein homology searching.

## Loading annotated sequences <link src="gAKBxE"/>

The first thing we'll do as we learn about sequence homology searching is load some annotated sequences. The sequences that we're going to work with are derived from the [Greengenes](http://greengenes.secondgenome.com/) database, and we're accessing them using the [QIIME default reference project](https://github.com/biocore/qiime-default-reference). Greengenes is a database of 16S rRNA gene sequences, a component of the archaeal and bacterial [ribosome](http://www.nature.com/scitable/definition/ribosome-194) (the molecular machine that drives translation of mRNA to proteins). This gene is of a lot of interest to biologists because it's one of about 200 genes that are encoded in the genomes of all known cellular organisms. We'll come back to this gene a few times in the book, notably in [Studying Biological Diversity](alias://2bb2cf). The sequences in Greengenes are taxonomically annotated, meaning that we'll have a collection of gene sequences and the taxonomic identity of the organism whose genome the sequence is found in. If we search an unannotated 16S rRNA query sequence against this database, we can make inferences about what organism our query sequence is from.

First, let's load Greengenes into a list of ``skbio.DNA`` sequence objects, and associate the taxonomy of each sequence as sequence metadata.

```python
>>> %pylab inline
...
>>> from IPython.core import page
>>> page.page = print
```

```python
>>> from iab.algorithms import load_taxonomy_reference_database
...
>>> %psource load_taxonomy_reference_database
```

```python
>>> reference_taxonomy, reference_db = load_taxonomy_reference_database()
```

Next, we'll just inspect a couple of the sequences we loaded. Notice how the specificity of our taxonomic annotations (i.e., how many taxonomic levels are annotated and unknown) differs for different sequences.

```python
>>> reference_db[0]
```

```python
>>> reference_db[-1]
```

For the sake of runtime, we're going to work through this chapter using a random sample of sequences from this database. Here we'll use Python's [random module](https://docs.python.org/3/library/random.html) to select sequences at random.

```python
>>> import random
...
>>> reference_db = random.sample(reference_db, k=5000)
>>> print("%s sequences are present in the subsampled database." % len(reference_db))
```

We'll also extract some sequences from Greengenes to use as query sequences in our database searches. This time we won't annotate them (to simulate not knowing what organisms they're from). We'll also trim these sequences so they're shorter than the full length references. This will simulate obtaining a partial gene sequence, as is most common with the current sequencing technologies (as of this writing), but will also help to make the examples run faster.

Note that some of our query sequences may also be in our subsampled reference database and some won't. This is realistic: sometimes we're working with sequences that are exact matches to known sequences, and sometimes we're working with sequences that don't match any known sequences (or at least any in the reference database that we're working with).

```python
>>> from iab.algorithms import load_taxonomy_query_sequences
...
>>> %psource load_taxonomy_query_sequences
```

```python
>>> queries = load_taxonomy_query_sequences()
>>> queries = random.sample(queries, k=50)
```

Let's inspect a couple of the query sequences that we'll work with.

```python
>>> queries[0]
```

```python
>>> queries[-1]
```

## Defining the problem <link src="SZ9u3S"/>

The problem that we are going to address here is as follows. We now have a query sequence ($q_i$) which is not taxonomically annotated (meaning we don't know the taxonomy of the organism whose genome it is found in), and a reference database ($R$) of taxonomically annotated sequences ($r_1, r_2, r_3, ... r_n$). We want to infer a taxonomic annotation for $q_i$. We'll do this by identifying the most similar sequence(s) in $R$ and associating their taxonomy with $q_i$. Because we actually do know the taxonomy of $q_i$ (to the extent that we trust the annotations in $R$), we can evaluate how well this approach works.

There are a few realistic features of the situation that we've set up here that I want you to be aware of.

1. All of the query and reference sequences are homologous. In this case, they are all sequences of the 16S rRNA gene from archaea and bacteria. This may or may not be the case in real-world applications. Sometimes you'll work with gene-specific databases such as Greengenes, and sometimes you'll work with non-specific databases such as the NCBI nucleotide database (nr). Regardless, the search process is similar.
2. The distance between each query sequence and its most closely related sequences in $R$ will vary widely. Sometimes $q$ will be an exact match to a reference sequence $r_i$, and sometimes we may have as little as $50\%$ similarity.

As we work through the next sections, imagine that we're exploring scaling this system up, so that instead of searching just one or a few query sequences against the reference database, we ultimately want to apply this to search millions of sequences against the database. This would be the real-world problem we faced if we had collected 16S rRNA sequences from the environment (which would of course be unannotated) using high-throughput DNA sequencing.

## A complete homology search function <link src="0C9FCS"/>

Let's define a homology search function that aligns each provided query sequences $q_i$ with each of our reference database sequences ($r_1, r_2, r_3, ... r_n$). This function will take as input one or more query sequences, and the reference database. We'll call the top scoring alignments for each $q_i$ the *best hits*, and we'll specifically request some number (`n`) of best hits for each $q_i$. The output of this function will be a summary of the `n` best hits for each query sequence, including some technical information about the alignment and the taxonomy associated with the corresponding reference sequence. We'll then review the taxonomy annotations for our best hits, and from those make an inference about the taxonomy annotation for $q_i$.

Spend a minute looking at this function and try to understand what it's doing.

```python
>>> from iab.algorithms import local_alignment_search
...
>>> %psource local_alignment_search
```

Now let's perform some database searches. You can run the remaining code cells in this section a few times to experiment with searching different query sequences against the same reference database.

This next cell, which is the one that actually performs the database searches, will take a little bit of time to run (maybe up to a minute or two). There is some code in this cell that will track the runtime. As it's running, think about how many query sequences we're searching against how many reference sequences, and refer back to the number of sequences in the full reference database. Does this strategy seem scalable to millions of sequences, which as mentioned above might be our ultimate goal? When you know the per-sequence runtime of this search, estimate how long it would take to do this in seconds for one million sequences. Convert the time in seconds to a unit that will be more meaningful to you.

```python
>>> import time
...
>>> start_time = time.time()
>>> current_queries = random.sample(queries, k=4)
>>> results = local_alignment_search(current_queries, reference_db)
>>> stop_time = time.time()
>>> print("Runtime: %1.4f sec per query" % ((stop_time - start_time) / len(current_queries)))
>>> results
```

Now, let's try to answer our initial question: what is the most likely taxonomic annotation for each of our query sequences? Spend a few minutes reviewing this information, and write down what you think the most likely taxonomic annotation is for each of the query sequences. Here are some hints to help you out:

 * The ``k``, ``p``, ``c``, ``o``, ``f``, ``g``, and ``s`` refer to *kingdom*, *phylum*, *class*, *order*, *family*, *genus*, and *species*, respectively. If you see an annotation for a reference sequence that looks like ``g__``, that means that the genus is unknown for that sequence.
 * Just as the reference taxonomy annotations don't always go down to the species level, your taxonomic annotations don't have to either. Not assigning at a given level implies that you're uncertain about what the annotation should be at that level, and it's usually better just to indicate that you're uncertain rather than make a bad guess. If you're uncertain of what the species is, assign the query ``s__`` and try to decide what the most likely genus is. If you're uncertain of the genus, assign ``g__``, and try to decide what the most likely family is...
 * As you look at each of the reference taxonomy annotations below, refer back to the table above to look at the percent similarity between each query and reference, and maybe the length of the alignments and their scores. These values give you an idea of how confident you should be in each of your taxonomic annotations.

```python
>>> for q in current_queries:
...     q_id = q.metadata['id']
...     print('Closest taxonomies for query %s (in order):' % q_id)
...     for e in results['reference taxonomy'][q_id]:
...         print(' ', e)
...     print()
```

Because we have taxonomic annotations for all of the Greengenes sequences (though as you probably have noticed by now, they differ in their specificity), we can next look at taxonomy associated with each of our queries in Greengenes. How do your annotations compare to those from Greengenes, which we'll print out in the next cell?

```python
>>> for q in current_queries:
...     q_id = q.metadata['id']
...     print('Known taxonomy for query %s:\n %s' % (q_id, reference_taxonomy[q_id]))
...     print()
```

## Reducing the runtime for database searches <link src='0f9232'/>

In the examples above, it's taking on the order of 5-15 seconds to search a single sequence against our subset of Greengenes. This makes sense when you think about the computations that are being performed. For every sequence in our reference database (5000, if you haven't modified the database subsampling step) it is computing the $F$ and $T$ matrices described in [the Pairwise Alignment chapter](alias://a76822), and then tracing back the matrix to compute the aligned sequences. Given all of that, the fact that computation only takes 5-15 seconds is pretty incredible. However, that doesn't change the fact that this doesn't scale to real-world applications because we'd have to wait way too long for results. Performing all pairwise alignments is prohibitively expensive for database searching.

As we discussed in the previous chapter, the run time of pairwise alignment scales quadratically with sequence length. Database searching, at least in the example we're exploring in this chapter, is a bit of a different problem however. Our sequence lengths aren't changing, but rather it takes a long time because we're performing a computationally expensive step, pairwise alignment, many times. Our database is fixed in that the number of sequences in it doesn't change and the sequences themselves don't change. Our query sequences are all exactly the same length in this example (remember that we set that above, when we sliced a single region from reference database sequences to create our query sequences). Let's explore how the runtime of this database search scales under these constraints.

```python
>>> import pandas as pd
>>> import itertools
...
>>> def tabulate_local_alignment_search_runtime(queries, reference_db, n_query_sequences,
...                                             n_reference_sequences, search_function):
...     data = []
...     # we'll iterate over the pairs of number of query sequences
...     # and number of reference sequences, and compute the runtime
...     # of the database search three times for each pair (so we
...     # have some idea of the variance in the runtimes). this is
...     # achieved here with a nested for loop (i.e., a for loop
...     # within a for loop).
...     for nq, nr in itertools.product(n_query_sequences, n_reference_sequences):
...         for i in range(3):
...             # select nq query sequences at random
...             current_queries = random.sample(queries, k=nq)
...             # select nr reference sequences at random
...             temp_reference_db = random.sample(reference_db, k=nr)
...             # run the search and store its runtime
...             start_time = time.time()
...             _ = search_function(current_queries, temp_reference_db)
...             stop_time = time.time()
...             median_query_sequence_len = np.median([len(q) for q in current_queries])
...             median_reference_sequence_len = np.median([len(r) for r in temp_reference_db])
...             data.append((nq, nr, median_query_sequence_len, median_reference_sequence_len,
...                          stop_time - start_time))
...     runtimes = pd.DataFrame(data=np.asarray(data),
...                             columns=["Number of query seqs", "Number of reference seqs",
...                                      "Median query seq length", "Median reference seq length",
...                                      "Runtime (s)"] )
...     return runtimes
...
>>> # we'll temporarily work with a smaller reference database
... # so this will run a lot faster. this will be of fixed size.
... n_reference_sequences = [100]
>>> # since our database is smaller, we can work with some slightly
... # larger numbers of sequences.
... n_query_sequences = [1, 5, 10, 15]
...
>>> local_alignment_search_runtimes = tabulate_local_alignment_search_runtime(queries, reference_db,
...                                                                           n_query_sequences, n_reference_sequences,
...                                                                           local_alignment_search)
>>> local_alignment_search_runtimes
```

This table shows that we've tried a few variations on number of query sequences but kept the number of reference sequences constant. There is no variance in the query sequence length, and there is a relatively small amount of variance in reference sequence length (they're all of the same order of magnitude). There is also relatively little variance in runtime for fixed numbers of query and reference sequences.

This table clearly shows that there is an increase in runtime with an increasing number of query sequences, which we'd of course expect. What we care about is how runtime is increasing as a function of number of query sequences. Let's plot runtime versus the number of query sequences to help us understand that relationship.

```python
>>> import seaborn as sns
>>> ax = sns.regplot(x="Number of query seqs", y="Runtime (s)", data=local_alignment_search_runtimes)
>>> ax.set_xlim(0)
>>> ax.set_ylim(0)
>>> ax
```

What we see here is pretty clearly a linear relationship: $runtime \approx constant \times number\ of\ query\ sequences$. This is because as we increase the number of query sequences, we're increasing the number of pairwise alignments that we need to perform. If we have 5 queries and 10 reference sequences, we compute $5 \times 10 = 50$ pairwise alignments. If we have 10 queries and 100 reference sequences, we compute $10 \times 100 = 1000$ pairwise alignments. There are a few practical ways to reduce the runtime of a process like this.

The first seems obvious, and even silly at first: perform fewer alignments. This could be achieved in a few ways. You could reduce the number of query sequences, though this might be something a researcher is resistant to: they have some collection of unknown sequences, and they want to know what they all are. You could alternatively reduce the number of reference sequences, but you might run into the same issues there: we wouldn't want to exclude reference sequences that might provide us with useful information about our query sequences. Finally, we might be able to figure out some ways to perform fewer alignments by not searching all of the query sequences against all of the reference sequences. If we could come up with some procedure to approximate which pairwise alignments were likely to be good (i.e., high scoring) and which were likely to be bad (i.e., low scoring) that is faster than performing the pairwise alignments, we could apply that procedure and only align a pair of sequences when we expect to get a high score. That could potentially allow us to reduce the number of alignments we need to perform, and therefore the runtime of the algorithm.

Another approach to reducing the runtime of this process would be to create a faster implementation of the algorithm (though at some point that won't be possible anymore), use a faster computer, or run the process in parallel on multiple processors. All of these would be ways to reduce the runtime of the search by some factor $f$, where $new\ runtime \approx \frac{runtime}{f}$.

In practice, for a production-scale sequence database search application like BLAST, we'd combine these approaches. In the next section we'll explore ways to reduce the runtime of database searching for a fixed number of query sequences and a fixed number of reference sequences by reducing the number of pairwise alignments that the search function will perform.

## Heuristic algorithms <link src="mUArdw"/>

As mentioned above, it just takes too long to search individual query sequences against a large database. This problem also isn't going away anytime soon. While computers are getting faster (or cheaper), the size of our sequences collections are getting bigger because sequencing is getting cheaper. In fact, many people think that obtaining DNA sequences is getting cheaper faster than computers are getting cheaper. As our number of query sequences increases because we are able to obtain more for the same amount of money, and the size of our reference databases increases (because we're continuously obtaining more sequence data) this will increasingly become a bigger problem. Figures 1 and 2, respectively, illustrate that these are both real-world issues. Notice that the axes are on a log scale in both cases.

```python
>>> import IPython.display
>>> IPython.display.IFrame(width="600", height="394", src="https://docs.google.com/spreadsheets/d/1vUkUuZsRlLW5U05rXXUn8B2sDYwShkClRMGa8Wiu6bc/pubchart?oid=1844125885&amp;format=interactive")
```

Figure 1: Genome sequencing costs.

```python
>>> import IPython.display
>>> IPython.display.IFrame(width="763", height="371", src="https://docs.google.com/spreadsheets/d/1vUkUuZsRlLW5U05rXXUn8B2sDYwShkClRMGa8Wiu6bc/pubchart?oid=2103353397&amp;format=interactive")
```

Figure 2: Size of GenBank.

One way that we can deal with this problem is by recognizing that most of the alignments that are performed in a database search are unlikely to be very good alignments. An algorithm developer could therefore improve runtime by defining a heuristic (or a rule) that is applied to determine which reference sequences are likely to result in good alignments, and only aligning the query against those. For it to be useful, making the decision to align or not (i.e., applying the heuristic) must be *much faster* than actually performing the pairwise alignment. The heuristic also needs to make *good* choices about which reference sequences to align the query against. If the algorithm chooses to not align against a specific reference, that reference is ruled out as a possible result of the database search. A good heuristic for sequence homology searching would therefore be very unlikely to exclude the best alignment(s). When thinking about heuristic algorithms in general, there are some important considerations:

1. How often does the heuristic algorithm fail to get the right answer (in our case, does it make good choices about which reference sequences to align against)?
2. How much faster is the heuristic than the "complete" approach, and is that reduction in runtime enough to justify not being guaranteed to get the best answer?

We'll now look at a few heuristics in the context of these questions.

### Random reference sequence selection <link src="bEQxHf"/>

Our first heuristic will be a [straw man](https://en.wikipedia.org/wiki/Straw_man) that we use as a baseline. We'll select a random $p\%$ of the reference sequences to align our query against. This will clearly result in a large decrease in the number of sequence alignments that we need to perform because we'll go from performing $R_s$ (the reference database size) sequence alignments to $p \times R_s$ sequence alignments for each query sequence $q_i$.

Here's the source code for this. You can see that we're just wrapping our ``local_alignment_search`` function in a function that samples down to $p\%$ of the reference sequences.

```python
>>> from iab.algorithms import heuristic_local_alignment_search_random
>>> %psource heuristic_local_alignment_search_random
```

Let's select some new queries and see how the results compare to our known taxonomies.

```python
>>> current_queries = random.sample(queries, k=10)
```

```python
>>> results = heuristic_local_alignment_search_random(current_queries, reference_db, p=0.10)
...
>>> for q in current_queries:
...     q_id = q.metadata['id']
...     print('Closest taxonomies for query %s (in order):' % q_id)
...     for e in results['reference taxonomy'][q_id]:
...         print(' ', e)
...     print()
```

```python
>>> for q in current_queries:
...     q_id = q.metadata['id']
...     print('Known taxonomy for query %s:\n %s' % (q_id, reference_taxonomy[q_id]))
```

What we need now is a way to know how often we get the "right answer", and how long this heuristic algorithm takes relative to the complete algorithm. We therefore first need to define what the "right answer" is. How about this: if the most common taxonomy assignment resulting from the database search at `taxonomy_levels` levels of taxonomy (i.e., how deep or specific our assignment is) matches the known taxonomy, then our algorithm has achieved the right answer. We can vary `taxonomy_levels` to see how the different heuristics perform at different levels.

Here's what this would look like:

```python
>>> import collections
...
>>> def evaluate_search(queries, reference_db, reference_taxonomy, search_function, taxonomy_levels, n=5, aligner=local_pairwise_align_ssw):
...     start_time = time.time()
...     search_results = search_function(current_queries, reference_db, n=n, aligner=aligner)
...     stop_time = time.time()
...     runtime = stop_time - start_time
...     per_query_runtime = runtime/len(queries)
...     data = []
...     indices = []
...     for q in queries:
...         q_id = q.metadata['id']
...         indices.append(q_id)
...         q_known_taxonomy = tuple(reference_taxonomy[q_id].split('; ')[:taxonomy_levels])
...         q_observed_taxonomies = collections.Counter()
...         for e in search_results['reference taxonomy'][q_id]:
...             q_observed_taxonomies[tuple(e.split('; ')[:taxonomy_levels])] += 1
...         q_observed_taxonomy = q_observed_taxonomies.most_common()[0][0]
...         data.append((q_known_taxonomy, q_observed_taxonomy))
...     index = pd.Index(indices, name='Query ID')
...     data = pd.DataFrame(data, index=index, columns=['Known taxonomy', 'Observed taxonomy'])
...     number_correct = np.sum(data['Known taxonomy'] == data['Observed taxonomy'])
...     fraction_correct = number_correct / data.shape[0]
...     return per_query_runtime, fraction_correct, data
```

First let's see how this works for our full database search algorithm. What's the runtime, and how often do we get the correct answer? We'll start with five levels of taxonomy (which corresponds to the family level). **This step will take a couple of minutes to run, because it's doing the full database search.**

```python
>>> taxonomy_levels = 5
```

```python
>>> runtime, fraction_correct, data = evaluate_search(current_queries, reference_db, reference_taxonomy,
...                                                   local_alignment_search, taxonomy_levels=taxonomy_levels)
>>> print('%1.2f seconds per query sequence' % runtime)
>>> print('%1.2f%% correct answers' % (fraction_correct * 100.0))
>>> print('Result details:')
>>> for q_id in data.index:
...     print(q_id)
...     print(' ', data['Known taxonomy'][q_id])
...     print(' ', data['Observed taxonomy'][q_id])
...     print()
```

Next let's see how this compares to our random heuristic search algorithm. Try running this a few times, as you might get different answers due to different random selections of the database.

```python
>>> import functools
...
>>> heuristic_local_alignment_search_random_10 = functools.partial(heuristic_local_alignment_search_random, p=0.10)
...
>>> runtime, fraction_correct, data = evaluate_search(current_queries, reference_db, reference_taxonomy,
...                                                   heuristic_local_alignment_search_random_10, taxonomy_levels=taxonomy_levels)
...
>>> print('%1.2f seconds per query sequence' % runtime)
>>> print('%1.2f%% correct answers' % (fraction_correct * 100.0))
>>> print('Result details:')
>>> for q_id in data.index:
...     print(q_id)
...     print(' ', data['Known taxonomy'][q_id])
...     print(' ', data['Observed taxonomy'][q_id])
...     print()
```

Again, what's the runtime, and how often do we get the correct answer? Based on comparison to the full search, what do you think: is this a good heuristic?

After performing many trials of the above searches, I get the correct genus-level assignment about half as often with the random reference database heuristic relative to the full database search. Your results might differ from that due to differences in the random selection of query and reference sequences. Try running all the cells in this section a few times.

Go back to the beginning of this section and try running this check based on fewer levels of taxonomy (i.e., decreased taxonomic specificity, such as the phylum) and on more levels of taxonomy (i.e., increased taxonomic specificity, such as the species level). How does that impact how often we get the right answer?

### Composition-based reference sequence collection <link src="P4vQ4b"/>

While the random selection of database sequences can vastly reduce the runtime for database searching, we don't get the right answer very often. Let's try some heuristics that are a bit smarter. How about this: if the overall nucleotide composition of a query sequence is very different than the overall nucleotide composition of a reference sequence, it's unlikely that the best alignment will result from that pairwise alignment, so don't align the query to that reference sequence. Given that, how do we define "overall nucleotide composition" in a useful way?

#### GC content <link src="8yiKFO"/>

One metric of sequence composition that we can compute quickly (because remember, this has to be a lot faster than computing the alignment for it to be worth it) is GC content. Let's define a heuristic that only performs a pairwise alignment for the reference sequences that have the most similar GC content to the query sequence. The number of alignments that we'll perform will be defined as ``database_subset_size``.

```python
>>> database_subset_size = 500
```

```python
>>> from iab.algorithms import heuristic_local_alignment_search_gc
...
>>> %psource heuristic_local_alignment_search_gc
```

If we run our queries again, how often do we get the right answer? How much did we reduce runtime? Do you think this is a better or worse heuristic than what we implemented above?

```python
>>> heuristic_local_alignment_search_gc_2 = functools.partial(heuristic_local_alignment_search_gc, database_subset_size=database_subset_size)
...
>>> runtime, fraction_correct, data = evaluate_search(current_queries, reference_db, reference_taxonomy,
...                                                   heuristic_local_alignment_search_gc_2, taxonomy_levels=taxonomy_levels)
...
>>> print('%1.2f seconds per query sequence' % runtime)
>>> print('%1.2f%% correct answers' % (fraction_correct * 100.0))
>>> print('Result details:')
>>> for q_id in data.index:
...     print(q_id)
...     print(' ', data['Known taxonomy'][q_id])
...     print(' ', data['Observed taxonomy'][q_id])
...     print()
```

Try increasing and decreasing the number of sequences we'll align by increasing or decreasing ``database_subset_size``. How does this impact the runtime and fraction of time that we get the correct answer?

#### kmer content <link src="QblTRV"/>

Another metric of sequence composition is *kmer composition*. A [kmer](alias://C7hMX5) is simply a word (or list of adjacent characters) of length *k* found within a sequence. Here are the kmer frequencies in a short DNA sequence. The ``overlap=True`` parameter here means that our kmers can overlap one another.

```python
>>> import skbio
...
>>> skbio.DNA('ACCGTGACCAGTTACCAGTTTGACCAA').kmer_frequencies(k=5, overlap=True)
```

In our next heuristic, we'll only align our query to the reference sequences with the largest fraction of the kmers that are observed in the query sequence are also present in the reference sequence. This makes a lot of sense to use as an alignment heuristic: we're only aligning sequences when it looks like they'll have multiple length-``k`` stretches of nucleotides that are not interrupted by substitutions or insertion/deletion mutations.

In our next heuristic, we'll only align our query to the reference sequences with the largest fraction of the kmers that are observed in the query sequence. This makes a lot of sense to use as an alignment heuristic: we're only aligning sequences when it looks like they'll have multiple length-``k`` stretches of nucleotides that are not interrupted by substitutions or insertion/deletion mutations.


Here's the source code:

```python
>>> from iab.algorithms import heuristic_local_alignment_search_kmers
...
>>> %psource heuristic_local_alignment_search_kmers
```

```python
>>> k = 7
```

Let's apply this and see how it does. How does the runtime and fraction of correct assignments compare to our GC content-based search and our full database search?

```python
>>> heuristic_local_alignment_search_kmers_50 = \
>>> functools.partial(heuristic_local_alignment_search_kmers, k=k, database_subset_size=database_subset_size)
...
>>> runtime, fraction_correct, data = evaluate_search(current_queries, reference_db, reference_taxonomy,
...                                                   heuristic_local_alignment_search_kmers_50,
...                                                   taxonomy_levels=taxonomy_levels)
...
>>> print('%1.2f seconds per query sequence' % runtime)
>>> print('%1.2f%% correct answers' % (fraction_correct * 100.0))
>>> print('Result details:')
>>> for q_id in data.index:
...     print(q_id)
...     print(' ', data['Known taxonomy'][q_id])
...     print(' ', data['Observed taxonomy'][q_id])
...     print()
```

#### Further optimizing composition-based approaches by pre-computing reference database information <link src="HQmZgF"/>

One important feature of composition-based approaches is that, because the reference database doesn't change very often, we can pre-compute features of the reference sequences and re-use them. This can help us to vastly decrease the runtime of our heuristic searches. For example, the computation of all of the reference database kmer frequencies is a lot of work. If we can compute that outside of our database search, we can avoid doing that step for every database search, and therefore remove that computationally expensive (i.e., slow) step of the process.

Here we'll compute all of the reference database kmer frequencies. Notice that this step takes about a minute to complete. This is a minute of compute time that we can save on every database search!

```python
>>> reference_db_kmer_frequencies = {r.metadata['id']: r.kmer_frequencies(k=k, overlap=True) for r in reference_db}
```

We'll now pass our pre-computed kmer frequencies into our search function. How does the runtime and accuracy of this search compare to the searches above? This last database search that we've implemented here is very similar to how BLAST works.

```python
>>> heuristic_local_alignment_search_kmers_50 = \
>>>  functools.partial(heuristic_local_alignment_search_kmers, reference_db_kmer_frequencies=reference_db_kmer_frequencies,
...                    k=k, database_subset_size=database_subset_size)
...
>>> runtime, fraction_correct, data = evaluate_search(current_queries, reference_db, reference_taxonomy,
...                                                   heuristic_local_alignment_search_kmers_50,
...                                                   taxonomy_levels=taxonomy_levels)
...
>>> print('%1.2f seconds per query sequence' % runtime)
>>> print('%1.2f%% correct answers' % (fraction_correct * 100.0))
>>> print('Result details:')
>>> for q_id in data.index:
...     print(q_id)
...     print(' ', data['Known taxonomy'][q_id])
...     print(' ', data['Observed taxonomy'][q_id])
...     print()
```

## Determining the statistical significance of a pairwise alignment <link src='87c92f'/>

One thing you may have noticed is that the score you get back for a pairwise alignment is hard to interpret. It's dependent on the query and reference sequence lengths (and possibly their composition, depending on your substitution matrix). So an important question is how to determine *how good* a given pairwise alignment is. Here we'll learn about a statistical approach for answering that.

### Metrics of alignment quality <link src="YdqCls"/>

In the examples above, we compared features such as how long the alignment is (relevant for local but not global alignment), the pairwise similarity between the aligned query and reference, and the score. If you've used a system like BLAST, you'll know that there are other values that are often reported about an alignment, like the number of substitutions, or the number of insertion/deletion (or gap) positions. None of these metrics are useful on their own. Let's look at an example to see why.

Imagine we're aligning these two sequences:

```
GAAGCAGCAC
GAACAGAAC
```

If we tell our search algorithm that we're interested in the alignment with the fewest number of substitutions, the following alignment would get us zero substitutions, but there are a lot of bases that look homologous which are not aligned.

```
GAAGCAGCAC-----
GAA------CAGAAC
```

On the other hand, if we want to find the alignment with the fewest number of gaps, this one would get us that result, but we now have a lot of substitution events, and some regions that clearly look misaligned (such as the ``CAG`` sequence in the middle of both).

```
GAAGCAGCAC
GAACAGA-AC
```

The alignment score that has been reported by our pairwise aligners helps us to balance these different features, and we can adjust the scoring scheme to weight things differently (e.g., so that gaps are penalized more or less than certain substitutions). The problem is that the scores are hard to interpret, particularly when we have only one or a few of them.

### False positives, false negatives, p-values, and alpha  <link src="TjqwVY"/>

Remember that an alignment of a pair of sequences represents a hypothesis about homology between those sequences. One way that we think about determining if an alignment is good or not is to ask: *what fraction of the time would I obtain a score at least this good if my sequences are not homologous?* This fraction is usually referred to as our *p-value*, and this is computed in many different ways. If our p-value is high (e.g., 25%), then our alignment is probably not very good since it means that many non-homologous pairs of sequences would achieve a score at least that high. If our p-value is low (say 0.001%), then our alignment is probably good since scores that high are achieved only infrequently.

Our threshold for defining what we consider to be a high versus low p-value is dependent on how often we're willing to be wrong. We would set this value, which is usually referred to as $\alpha$, to some fraction, and if our p-value is less than $\alpha$, we say that the alignment is statistically significant. If our p-value is greater than $\alpha$, we say that our alignment is not statistically significant.

There are a couple of ways that we could be wrong when we do sequence homology searching, and we need to consider these when we determine what value we want to define as $\alpha$. First, we could say a pair of sequences are homologous when they're not, which would be a *false positive* or a *type 1 error*. Or, we could say that a pair of sequences are not homologous when they are, which would be a *false negative*, or a *type 2 error*.

If incurring a false positive about 5% of the time is acceptable (i.e., you're ok with calling a pair of sequences homologous when they actually are not about one in twenty times) then you'd set your $\alpha$ to 0.05. Setting $\alpha$ to a value this high likely means that the method will err on the side of false positives, and only infrequently will it say that a pair of sequences are not homologous when they actually are (i.e., achieve a false negative). If $\alpha$ were set to be very low on the other hand (say, $1 \times 10^{-50}$), then you will err on the side of false negatives. Only infrequently will you say that a pair of non-homologous sequences are homologous, but you might call many pairs of homologous sequences non-homologous. You should think of $\alpha$ as a dial. If you turn the dial toward higher values, you'll increase your false positive rate and decrease your false negative rate. If you turn the dial toward lower values, you'll decrease your false positive rate and increase your false negative rate.

There is not a hard-and-fast rule for whether false positives or false negatives are better, which makes choosing $\alpha$ hard. It's application specific, so you need to understand the biological question your asking when making this decision, and the ramifications of false positives versus false negatives. In general, when might you prefer to have false positives? When might you prefer to have false negatives?

### Interpreting alignment scores in context <link src="a0nqBH"/>

In this section, we are going to learn about how to interpret alignment scores by empirically determining if a pairwise alignment that we obtain is better than we would expect if the pair of sequences we're working with were definitely not homologous. For a given pair of sequences that we want to align, we're first going to align them and compute the score of the alignment. We're then going to align many pairs of sequences that are similar to the query and reference, but that we know are not homologous. We'll do this by shuffling or randomizing the order of the bases in the query sequences, and performing another pairwise alignment.

First, we'll define a function that can generate random sequences for us. This will take a scikit-bio sequence object (either ``skbio.DNA``, ``skbio.RNA``, or ``skbio.Protein``) and a length, and it will randomly generate a sequence of that type and length for us.

```python
>>> import random
>>> def random_sequence(moltype, length):
...     result = []
...     alphabet = list(moltype.nondegenerate_chars)
...     for e in range(length):
...         result.append(random.choice(alphabet))
...     return moltype(''.join(result))
```

We can now run this a few times to generate some random sequences:

```python
>>> random_sequence(skbio.DNA, 50)
```

```python
>>> random_sequence(skbio.DNA, 50)
```

Next, we need a function that will shuffle the characters in a sequence, and give us a new sequence back. We'll use this to generate a sequence that is similar (in length and composition) to our input sequence, but which we know is not homologous. We'll use Pythons `random.shuffle` function, which randomly re-orders the order of the elements in a sequence, but keeps the composition and length of the sequence the same.

```python
>>> from iab.algorithms import shuffle_sequence
>>> %psource shuffle_sequence
```

Now we can define a random sequence and shuffle it. Notice how the sequences are different (in their order), but their compositions (e.g., length and GC content) are the same. Shuffling will change the order of the bases, but it won't change the frequency at which each base is present - it's exactly analogous to shuffling a deck of cards.

```python
>>> seq = random_sequence(skbio.DNA, 50)
>>> seq
```

```python
>>> shuffle_sequence(seq)
```

Let's generate a random query sequence and align it against itself to see what that score would be.

```python
>>> query_seq = random_sequence(skbio.DNA, 50)
>>> _, actual_score, _ = local_pairwise_align_ssw(query_seq, query_seq)
>>> print("Score: %1.2f" % actual_score)
```

Next let's generate 99 random variants of that sequence with ``shuffle_sequence`` and compute the pairwise alignment for each of those variants against the query sequence. We'll then look at the distribution of those scores.

```python
>>> from iab.algorithms import generate_random_score_distribution
>>> %psource generate_random_score_distribution
```

```python
>>> random_scores = generate_random_score_distribution(query_seq, query_seq, 99)
>>> print(random_scores)
```

How does the actual score of aligning the sequence to itself compare to the score of aligning it to many similar but non-homologous sequences? Let's plot these to get a better idea.

```python
>>> import seaborn as sns
...
>>> def plot_score_distribution(actual_score, random_scores):
...     ax = sns.distplot(random_scores, kde=False, label="Random scores", color="b")
...     ax.plot([actual_score, actual_score], ax.get_ylim(), '--', label="Actual score")
...     # set the range of the x axis to be zero through 110% of the actual score
...     ax.set_xlim(0, actual_score + actual_score * 0.1)
...     ax.legend(loc=9, fontsize='large')
...     return ax
```

```python
>>> plot_score_distribution(actual_score, random_scores)
```

What does this tell us about our alignment score and therefore about our alignment? Is it good or bad?

We finally have information that we can use to evaluate an alignment score, and therefore to evaluate the quality of an alignment. Let's use this information to quantify the quality of the alignment by computing a p-value. As we described above, this is simply the probability that we would obtain an alignment score at least this good if the sequences being aligned are not homologous. Since we have a lot of scores now from sequences that are similar but not homologous, if we just count how many are at least as high as our actual score and divide by the number of scores we compute, that is an empirical (data-driven) way of determining our p-value.

To determine if our alignment is statistically significant, we need to define $\alpha$ before computing the p-value so the p-value does not impact our choice of $\alpha$. Let's define $\alpha$ as 0.05. This choice means if we obtain a p-value less than 0.05 we will consider the alignment statistically significant and accept the hypothesis that the sequences are homologous.

Here's what all of this looks like:

```python
>>> from iab.algorithms import fraction_better_or_equivalent_alignments
>>> %psource fraction_better_or_equivalent_alignments
```

```python
>>> print("Fraction of alignment scores at least as good as the alignment score: %r" %
...       fraction_better_or_equivalent_alignments(query_seq, query_seq, 99))
```

The fraction that we get back here is ``0.01``, which is lower than $\alpha$, so we would accept the hypothesis that our sequences are homologous.

A few notes on these empirically defined p-values. First, here's what the formula for computing this looks like:

$p\ value = \frac{number\ of\ computed\ aligned\ scores\ greater\ than\ or\ equal\ to\ the\ actual\ alignment\ score}{number\ of\ alignment\ scores\ computed}$

The numerator and the denominator both include the actual alignment score, so the lowest p-value that can be achieved is $\frac{1}{99 + 1}$, where the $1$ in the numerator corresponds to our actual alignment score (which is of course equal to itself), where the $99$ in the denominator is the number of permutations, and the $1$ in the denominator is a constant which corresponds the computation of the actual score. If we increase the number of permutations, say to 999, we could achieve greater precision (more significant digits) in our p-value.

```python
>>> print("Fraction of alignment scores at least as good as the alignment score: %r" %
...       fraction_better_or_equivalent_alignments(query_seq, query_seq, 999))
```

When we achieve the lowest possible value for a given test, as is the case here, we report the p-value as being less than that value, since we've yet to observe a random alignment score at least that high. For example, here we would report something like:

*The alignment of our query and reference sequence was statistically significant, as determined by comparing our actual alignment score to random variants ($p < 0.001$).*

Let's now try this for some harder cases, where the query and subject sequences are not identical. First, let's generate a longer subject sequence at random. Then, we'll create a random query sequence and compare it. Since we're doing this in two random steps, we know that these sequences are not homologous. Does the resulting p-value reflect that?

```python
>>> sequence1 = random_sequence(skbio.DNA, 250)
>>> sequence1
```

```python
>>> sequence2 = random_sequence(skbio.DNA, 250)
>>> sequence2
```

```python
>>> print("Fraction of alignment scores at least as good as the alignment score: %r" %
...       fraction_better_or_equivalent_alignments(sequence1,sequence2))
```

We've now looked at two extremes: where sequences are obviously homologous (because they were the same), and where sequences are obviously not homologous (because they were both independently randomly generated). Next, we'll explore the region between these, where this gets interesting. We'll now create a partially randomized sequence to create a pair of sequences where the homology is more obscure. We'll do this again using the Python ``random`` module, but this time we'll introduce mutations only at some positions to create a pair of sequences that are approximately ``percent_id`` identical.

Let's define a function to do this, and then compute a sequence that is 95% identical to our ``sequence1``.

```python
>>> def partially_randomize_sequence(percent_id, sequence):
...     result = []
...     for c in sequence:
...         if random.random() < percent_id:
...             result.append(str(c))
...         else:
...             # choose a base at random that is not the current base
...             # i.e., simulate a substitution event
...             result.append(choice([r for r in sequence.nondegenerate_chars if r != c]))
...     return sequence.__class__(''.join(result))
```

```python
>>> sequence1_95 = partially_randomize_sequence(0.95, sequence1)
```

```python
>>> sequence1
```

```python
>>> sequence1_95
```

Notice how these sequences are almost identical, but have some differences. Let's apply our approach to determine if it would identify these sequences as being homologous based on $\alpha = 0.05$.

```python
>>> print("Fraction of alignment scores at least as good as the alignment score: %r" %
...       fraction_better_or_equivalent_alignments(sequence1, sequence1_95))
```

You likely got a significant p-value there, telling you that the sequences are homologous.

Now let's simulate much more distantly related sequences by introducing substitutions at many more sites.

```python
>>> sequence1_25 = partially_randomize_sequence(0.25, sequence1)
```

```python
>>> sequence1
```

```python
>>> sequence1_25
```

```python
>>> print("Fraction of alignment scores at least as good as the alignment score: %r" %
...       fraction_better_or_equivalent_alignments(sequence1, sequence1_25))
```

### Exploring the limit of detection of sequence homology searches <link src="n8OGGW"/>

In the example above, we know that our input sequences are "homologous" because `sequence1_25` and `sequence1_95` are both derived from `sequence1`. Our method detected that homology for `sequence1_95`, when we simulated very closely related sequences, but not for ``sequence1_25``, when we simulated much more distantly related sequences. This gives us an idea of the limit of detection of this method, and is a real-world problem that biologists face: as sequences are more divergent from one another, detecting homology becomes increasingly difficult.

Lets run a simulation to gain some more insight into the limit of detection of this method. We'll run this approach for pairs of sequences where we vary the ``percent_id`` parameter, and identify when our approach stops identifying sequence pairs as being homologous. This is important to know as a bioinformatician, because it tells us around what pairwise similarity we will no longer be able to identify homology using this approach.

```python
>>> # First, let's define the range of percent identities that we'll test
... percent_ids = np.arange(0.0, 1.0, 0.05)
>>> # Then, we'll define the number of random sequences we'll test at each percent identity
... num_trials = 20
>>> # Then, we'll define the sequence length that we want to work with, and num_trials random sequences
... sequence_length = 150
>>> random_sequences = [random_sequence(skbio.DNA, sequence_length) for i in range(num_trials)]
...
>>> results = []
...
>>> for percent_id in percent_ids:
...     # at each percent_id, we'll track the p-values for each trial (random sequence)
...     p_values = []
...     for sequence in random_sequences:
...         # partially randomize the sequence, compute its p-value, and record that p-value
...         sequence_at_percent_id = partially_randomize_sequence(percent_id, sequence)
...         p = fraction_better_or_equivalent_alignments(sequence, sequence_at_percent_id)
...         p_values.append(p)
...     results.append((percent_id, np.median(p_values), np.mean(p_values)))
>>> pd.DataFrame(results, columns=["Percent id between query and subject",
...                                "Median p-value", "Mean p-value"])
```

What does this simulation tell us about our limit of detection for homology (i.e., how similar must a pair of sequences be for us to reliably be able to identify homology between them)? Is this higher or lower than you expected?

With respect to our simulation, I took a few shortcuts here to keep the runtime low. What are some things that could be improved to make this simulation more robust, if we weren't as concerned about runtime?
