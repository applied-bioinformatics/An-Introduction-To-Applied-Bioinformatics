
# Sequence homology searching <link src='d22e6b'/>

In this chapter we'll talk about using pairwise alignment to search databases of biological sequences with the goal of identifying sequence homology. We [previously defined homology](alias://e63a4f) between a pair of sequences to mean that those sequences are derived from a common ancestral sequence. Homology searching is an essential part of making inferences about where a biological sequence came from, and/or what it does. In most cases, if you have an unannotated biological sequence, such as the following protein sequence, it's very hard (really, impossible) to know what it is without more information.

```
> mystery-sequence1
PQITLWQRPLVTIRIGGQLKEALLDTGADDTVLEEMNLPGKWKPKMIGGIGGFIKVRQYDQIPVEIAHKAIGTVLVGPTPVNIIGRNLLTQIGATLNF
```

What a researcher will often do is search this sequence, their *query*, against some *reference database* of annotated sequences to learn what function the sequence performs (if the reference database contains functional annotation of sequences) and/or what organisms are likely to encode this sequence in their genome (if the reference database contains taxonomic annotation of sequences).

Whose genome is the above sequence encoded in? What is its function? Take a minute now to answer these questions using the [Protein BLAST homology search tool on the NCBI website](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp).

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
>>> locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
...
>>> # Load the taxonomic data
... reference_taxonomy = {}
>>> for e in open(qdr.get_reference_taxonomy()):
...     seq_id, seq_tax = e.strip().split('\t')
...     reference_taxonomy[seq_id] = seq_tax
...
>>> # Load the reference sequences, and associate the taxonmic annotation with
... # each as metadata
... reference_db = []
>>> for e in skbio.io.read(qdr.get_reference_sequences(), format='fasta', constructor=skbio.DNA):
...     seq_tax = reference_taxonomy[e.metadata['id']]
...     e.metadata['taxonomy'] = seq_tax
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

For the sake of runtime, we're going to work through this chapter using a random sample of sequences from this database. Here we'll use Python's [random module](https://docs.python.org/3/library/random.html) to select sequences at random.

```python
>>> import random
...
>>> reference_db = random.sample(reference_db, k=5000)
>>> print("%s sequences are present in the subsampled database." % locale.format("%d", len(reference_db), grouping=True))
```

We'll also extract some sequences from Greengenes to use as query sequences in our database searches. This time we won't annotate them (to simulate not knowing what organisms they're from). We'll also trim these sequences so they're shorter than the full length references. This will simulate obtaining a partial gene sequence, as is most common with the current sequencing technologies (as of this writing), but will also help to make the examples run faster.

Note that some of our query sequences may also be in our subsampled reference database and some won't. This is realistic: sometimes we're working with sequences that are exact matches to known sequences, and sometimes we're working with sequences that don't match any known sequences (or at least any in the reference database that we're working with).

```python
>>> queries = []
>>> for e in skbio.io.read(qdr.get_reference_sequences(), format='fasta', constructor=skbio.DNA):
...     e = e[100:300]
...     queries.append(e)
>>> queries = random.sample(queries, k=500)
```

Let's inspect a couple of the query sequences that we'll work with.

```python
>>> queries[0]
```

```python
>>> queries[-1]
```

## Defining the problem <link src="SZ9u3S"/>

The problem that we are going to address here is as follows. We now have a query sequence ($q_i$) which is not taxonomically annotated (meaning we don't know the taxonomy of the organism whose genome it is found in), and a reference database ($R$) of taxonomically annotated sequences ($r_1, r_2, r_3, ... r_n$). We want to infer a taxonomic annotation for $q_i$. We'll do this by identifying the most closely evolutionary related sequences in $R$ and associating their taxonomy with $q_i$. Because we actually do know the taxonomy of $q_i$ (to the extent that we trust the annotations in $R$), we can evaluate how well this approach works.

There are a few realistic features of the situtation that we've set up here that I want you to be aware of.

1. All of the query and reference sequences are homologous. In this case, they are all sequences of the 16S rRNA gene from archaea and bacteria. This may or may not be the case in real-world applications. Sometimes you'll work with gene-specific databases such as Greengenes, and sometimes you'll work with non-specific databases such as the NCBI nucleotide database (nr). Regardless, the search process is similar.
2. The evolutionary distance between each query sequence and its most closely related sequences in $R$ will vary widely. Sometimes $q$ will be an exact match to a reference sequence $r_i$, and sometimes we may have as little as $50\%$ similarity.

As we work through the next sections, imagine that we're exploring scaling this system up, so that instead of search just one or a few query sequences against the reference database, we ultimately want to apply this to search millions of sequences against the database. This would be the real-world problem we faced if we had collected 16S rRNA sequences from the environment (which would of course be unannotated) using high-throughput DNA sequencing.

## A complete homology search function <link src="0C9FCS"/>

Let's define a homology search function that aligns each provided query sequences $q_i$ with each of our reference database sequences ($r_1, r_2, r_3, ... r_n$). This function will take as input one or more query sequences, and the reference database. We'll call the top scoring alignments for each $q_i$ the *best hits*, and we'll specifically request some number (`n`) of best hits for each $q_i$. The output of this function will be a summary of the `n` best hits for each query sequence, including some technical information about the alignment and the taxonomy associated with the corresponding reference sequence. We'll then review the taxonomy annotations for our best hits, and from those make an inference about the taxonomy annotation for $q$.

Spend a minute looking at this function to understand what it's doing.

```python
>>> from iab.algorithms import local_alignment_search
...
>>> %psource local_alignment_search
```

Now let's perform some database searches. You can run the remaining code cells in this section a few times to experiment with searching different query sequences against the same reference database.

This next cell, which is the one that actually performs the database searches, will take a little bit of time to run (maybe up to a minute or two). There is some code in this cell that will track the runtime. As it's running, think about how many query sequences we're searching against how many reference sequences, and refer back to the number of sequences in the full reference database. Does this strategy seem scalable to millions or sequences, which as mentioned above might be our ultimate goal? When you know the per-sequence runtime of this search, estimate how long it would take to do this in seconds for one million sequences. Convert the time in seconds to a unit that will be more meaningful to you.

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

As we discussed in the previous chapter, the run time of pairwise alignment scales quadratically with sequence length. Database searching, at least in the example we're exploring in this chapter, is a bit of a different problem however. Our sequence lengths aren't changing, but rather it takes a long time because we're performing a computationally expensive step, pairwise alignment, many times. Our database is fixed in that the number of sequences in it doesn't change and the sequences themselves don't change. Our query sequences are all exactly the same length in this example (remember that we set that above, when we sliced a single region from reference database sequences to create our query sequences). Let's explore how the runtime of this database search scales under these constriants.

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

This table shows that we've tried a few variations on number of query sequences but kept the number of reference sequences constant. The is no variance in the query sequence length, and there is a relatively small amount of variance in reference sequence length (they're all of the same order of magnitude). There is also relatively little variance in runtime for fixed numbers of query and reference sequences.

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

Another approach to reducing the runtime of this process would be to create a faster implemention of the algorithm (though at some point that won't be possible anymore), use a faster computer, or run the process in parallel on multiple processors. All of these would be ways to reduce the runtime of the search by some factor $f$, where $new\ runtime \approx \frac{runtime}{f}$.

In practice, for a production scale sequence database search application like BLAST, we'd combine these approaches. In the next section we'll explore ways to reduce the runtime of database searching for a fixed number of query sequences and a fixed number of reference sequences by reducing the number of pairwise alignments that the search function will perform.

## Heuristic algorithms <link src="mUArdw"/>

As mentioned above, it just takes too long to search individual query sequences against a large database. This problem also isn't going away anytime soon. While computers are getting faster (or cheaper), the size of our sequences collections are getting bigger because sequencing is getting cheaper. In fact, many people think that obtaining DNA sequences is getting cheaper faster than computers are getting cheaper. As our number of query sequences increases because we are able to obtain more for the same amount of money, or the size of our reference databases increases (because we're continuously obtaining more sequence data) this will increasing become a problem. Figures 1 and 2, respectively, illustrate that these are both real-world issues.

```python
>>> import IPython.display
>>> IPython.display.IFrame(width="600", height="394", src="https://docs.google.com/spreadsheets/d/1vUkUuZsRlLW5U05rXXUn8B2sDYwShkClRMGa8Wiu6bc/pubchart?oid=1844125885&amp;format=interactive")
```

Figure 1: Genome sequencing costs.

```python
>>> import IPython.display
>>> IPython.display.IFrame(width="763", height="371", src="https://docs.google.com/spreadsheets/d/1vUkUuZsRlLW5U05rXXUn8B2sDYwShkClRMGa8Wiu6bc/pubchart?oid=2103353397&amp;format=interactive")
```

Figure 2: Size of GenBank.</figcaption>

One way that we can deal with this problem is by recognizing that most of the alignments that are performed in a database search are unlikely to be the best alignment. An algorithm developer could therefore improve runtime by defining a heuristic (or a rule) that is applied to determine which sequences we're going to align and which sequences we're not going to align. For it to be useful, deciding whether to align a pair of sequences (i.e., applying the heurtistic) must be much faster than performing the pairwise alignment. It is important to note that if we decide to not align a query against a given reference sequence, we exclude that reference sequence as a possible result of our search. A good heuristic makes it very unlikely to exclude good alignments. When thinking about heuristic algorithms, there are some important considerations:

1. How often do I fail to get the right answer?
2. Is my runtime reduced enough that I'm willing to tolerate not getting the best alignment this often?

We'll now look at a few heuristics in the context of these questions.

### Random reference sequence selection <link src="bEQxHf"/>

Our first heuristic will be a [straw man](https://en.wikipedia.org/wiki/Straw_man) that we use as a baseline. We'll select a random $p\%$ of the reference sequences to align our query against. This will clearly result in a large decrease in the number of sequence alignments that we need to perform because we'll go from performing $R_s$ (the reference database size) sequence alignments to $p \times R_s$ sequence alignments for each query sequence $q_i$.

Here's the source code for this. You can see that we're just wrapping our ``local_alignment_search`` function in a function that samples down to $p\%$ of the reference sequences.

```python
>>> from iab.algorithms import heuristic_local_alignment_search_random
>>> %psource heuristic_local_alignment_search_random
```

Let's pass our initial `queries` to see how the results compare to our known taxonomies.

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

First let's see how this works for our full database search algorithm. What's the runtime, and how often do we get the correct answer? We'll start with five levels of taxonomy (which corresponds to the family level).

```python
>>> runtime, fraction_correct, data = evaluate_search(current_queries, reference_db, reference_taxonomy,
...                                                   local_alignment_search, taxonomy_levels=5)
>>> print('%1.2f seconds per query sequence' % runtime)
>>> print('%1.2f%% correct answers' % (fraction_correct * 100.0))
>>> print('Result details:')
>>> for q_id in data.index:
...     print(q_id)
...     print(' ', data['Known taxonomy'][q_id])
...     print(' ', data['Observed taxonomy'][q_id])
...     print()
```

Next let's see how this compare to our random heuristic search algorithm. Try running this a few times, as you might get different answers due to different random selections of the database.

```python
>>> import functools
...
>>> heuristic_local_alignment_search_random_10 = functools.partial(heuristic_local_alignment_search_random, p=0.10)
...
>>> runtime, fraction_correct, data = evaluate_search(current_queries, reference_db, reference_taxonomy,
...                                                   heuristic_local_alignment_search_random_10, taxonomy_levels=5)
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

### Composition-based reference sequence collection <link src="P4vQ4b"/>

While the random selection of database sequences can vastly reduce the runtime for database searching, we don't get the right answer very often. Let's try a heuristic that's a bit smarter. How about this: if the overall nucleotide composition of a query sequence is different than the overall nucleotide composition of a reference sequence, it's unlikely that the best alignment will result from that pairwise alignment. One metric of sequence composition that we can compute quickly (because remember, this has to be a lot faster than computing the alignment for it to be worth it) is GC content. Let's define a heuristic that only performs a pairwise alignment if the GC content of the reference sequence is within $p%%$ of the GC content of the query sequence.

First, let's get an idea of what the range of GC contents is for all of our database sequences, as that will help us decide on a reasonable value for ``p``.

```python
>>> sns.set(style="white", palette="muted")
>>> reference_db_gc_contents = {r.metadata['id'] : r.gc_content() for r in reference_db}
>>> ax = sns.distplot(list(reference_db_gc_contents.values()))
>>> _ = ax.set_xlim((0, 1))
```

The distribution of GC content values is fairly narrow, so we'll pick a small value for ``p``, maybe $2.5%%$, to keep the number of sequences that we search small. As you might image, the smaller we make ``p``, the faster our search will run, but we'll be more likely to obtain an incorrect taxonomy assignment.

```python
>>> from iab.algorithms import heuristic_local_alignment_search_gc
...
>>> %psource heuristic_local_alignment_search_gc
```

If we run our queries again, how often do we get the right answer? How much did we reduce runtime? Do you think this is a better or worse heurtistic?

```python
>>> heuristic_local_alignment_search_gc_2 = functools.partial(heuristic_local_alignment_search_gc, p=0.02)
...
>>> runtime, fraction_correct, data = evaluate_search(current_queries, reference_db, reference_taxonomy,
...                                                   heuristic_local_alignment_search_gc_2, taxonomy_levels=5)
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

** Revision of this chapter is in progress. Pick up here!**

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
...         if random.random() < percent_id:
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
