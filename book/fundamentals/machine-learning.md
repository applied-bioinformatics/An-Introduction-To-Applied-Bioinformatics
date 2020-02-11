# Machine learning in bioinformatics (work-in-progress) <link src="7SFnTr"/>

**This chapter is currently a work-in-progress, and is incomplete.**

Machine learning algorithms are commonly used in bioinformatics for a variety of tasks. Typically, the common thread in these tasks is that the user would like the algorithm to assist in the identification of patterns in a complex data set. In this chapter we'll implement a few machine learning algorithms so we can gain an in-depth understanding of how they work. In practice though, there are many mature machine learning libraries that you'd want to use. [scikit-learn](http://scikit-learn.org/) is a popular and well-documented Python library for machine learning which many bioinformatics researchers and software developers use in their work.

These algorithms generally work beginning with a collection of samples and some user-defined features of those samples. These data are typically represented in a matrix, where samples are the rows and features are the columns. There are a few different high-level tasks that are common in machine learning, including classification, regression, and dimensionality reduction. In a classification task, a user provides examples of data that fall into certain discrete classes (for example, _healthy_ and _disease_), and tries to have the computer develop a model that can differentiate those classes based on the defined features. If successful, the resulting model could be applied to data where the class isn't known ahead of time, in attempt to predict the class from the features. A regression task is similar, except that a continuous value will be predicted rather than a discrete value. Dimensionality reduction tasks, on the other hand, generally don't have classes or labels assigned ahead of time, and the user is hoping to identify which samples are most similar to each other based on new features that are defined by the algorithm. The goal here might be to reduce the number of features from thousands or more to around two or three that explain most of the variation in the data. This allows the user to explore the samples visually, for example in a scatter plot, which would not be feasible if there were thousands of features.

In this chapter we'll explore two classification algorithms and one dimensionality reduction task in the context of some real-world examples.

## Defining a classification problem <link src="2R6CTy"/>

We'll explore machine learning classifiers in the context of a familiar topic: taxonomic classification of 16S rRNA sequences. We previously explored this problem in [Sequence Homology Searching](alias://d22e6b), so it is likely worth spending a few minutes skimming that chapter if it's not fresh in your mind.

Briefly, the problem that we are going to address here is as follows. We have a query sequence ($q_i$) which is not taxonomically annotated (meaning we don't know the taxonomy of the organism whose genome it is found in), and a reference database ($R$) of taxonomically annotated sequences ($r_1, r_2, r_3, ... r_n$). We want to infer a taxonomic annotation for $q_i$. We'll again work with the [Greengenes](http://greengenes.secondgenome.com/) database, which we'll access using [QIIME default reference project](https://github.com/biocore/qiime-default-reference). Greengenes is a database of 16S rRNA gene sequences. (This should all sound very familiar - if not, I again suggest that you review [Sequence Homology Searching](alias://d22e6b).)

This time, instead of using sequence alignment to identify the most likely taxonomic origin of a sequence, we'll train classifiers by building [kmer](alias://C7hMX5)-based models of the 16S sequences of taxa in our reference database. We'll then run our query sequences through those models to identify the most likely taxonomic origin of each query sequence. Since we know the taxonomic origin of our query sequences in this case, we can evaluate the accuracy of our classifiers by seeing how often they return the known taxonomy assignment. If our training and testing approaches are well-designed, the performance on our tests will inform us of how accurate we can expect our classifier to be on data where the actual taxonomic origin is unknown.

Let's jump in...

### Naive Bayes classifiers <link src="H8vYPu"/>

The first classifier we'll explore is the popular and relatively simple Naive Bayes classifier. This classifier uses Bayes Theorem to determine the most likely label for an unknown input based on a probabilistic model it has constructed from training data. (_The preceding text needs work._) The model that is constructed is based on user-defined features of the sequences. The most commonly used features for sequence classification tasks such as this is overlapping [kmers](alias://C7hMX5).

We'll begin by importing some libraries that we'll use in this chapter, and then [preparing our reference database and query sequences as we did previously](alias://gAKBxE).

```python
>>> %pylab inline
...
>>> from IPython.core import page
>>> page.page = print
...
>>> import pandas as pd
>>> import skbio
>>> import numpy as np
>>> import itertools
>>> import collections
```

```python
>>> from iab.algorithms import load_taxonomy_reference_database
...
>>> %psource load_taxonomy_reference_database
```

```python
>>> reference_taxonomy, reference_db = load_taxonomy_reference_database()
```

```python
>>> reference_db[0]
```

```python
>>> reference_db[-1]
```

We'll select a random subset of the reference database to work with here.

```python
>>> reference_db = np.random.choice(reference_db, 500, replace=False)
>>> print("%s sequences are present in the subsampled database." % len(reference_db))
```

The first thing our Naive Bayes classifier will need is the set of all possible words of length ``k``. This will be dependent on the value of ``k`` and the characters in our alphabet (i.e., the characters that we should expect to find in the reference database). This set is referred to as ``W``, and can be computed as follows. Given the following alphabet, how many kmers of length 2 are there (i.e., 2-mers)? How many 7-mers are there? How many 7-mers are there if there are twenty characters in our alphabet (as would be the case if we were working with protein sequences instead of DNA sequences)?

```python
>>> alphabet = skbio.DNA.nondegenerate_chars
>>> k = 2
...
>>> def compute_W(alphabet, k):
>>>  return set(map(''.join, itertools.product(alphabet, repeat=k)))
...
>>> W = compute_W(alphabet, k)
>>> print('Alphabet contains the characters: %s' % ', '.join(alphabet))
>>> print('For an alphabet size of %d, W contains %d length-%d kmers.' % (len(alphabet), len(W), k))
```

scikit-bio provides methods for identifying all kmers in a ``skbio.DNA`` sequence object, and for computing the kmer frequencies. This information can be obtained for one of our reference sequences as follows:

```python
>>> kmers = reference_db[0].iter_kmers(k=k)
>>> for kmer in kmers:
...     print(kmer, end=' ')
```

```python
>>> print(reference_db[0].kmer_frequencies(k=k))
```

This information can be convenient to store in a pandas ``Series`` object:

```python
>>> pd.Series(reference_db[0].kmer_frequencies(k=k), name=reference_db[0].metadata['id'])
```

To train our taxonomic classifier, we next need to define a few things. First, at what level of taxonomic specificity do we want to classify our sequences? We should expect to achieve higher accuracy at less specific taxonomic levels such as phylum or class, but these are likely to be less informative biologically than more specific levels such as genus or species. Let's start classifying at the phylum level to keep our task simple, since we're working with a small subset of the reference database here. In Greengenes, phylum is the second level of the taxonomy.

Next, how long should our kmers be? We don't have a good idea of this to start with. The longer our kmers, the more likely they are to be specific to certain taxa, which is good because that will help with classification. However, if they get too long it becomes less likely that we'll observe those kmers in sequences that aren't represented in our database because the longer the sequence is the more likely we are to see variation across other organisms that are assigned to the same taxonomy. Based on some of my own work in this area, I'll start us out with 7-mers (i.e., kmers of length 7).

Finally, we'll need to know the value of `W`, defined above as the set of all possible kmers given our alphabet and the value of `k`.

As an exercise, I recommend exploring the impact of the value of `k` and `taxonomic_level` on the accuracy of our classifier after reading this chapter.

```python
>>> taxonomic_level = 2
>>> k = 7
>>> alphabet = skbio.DNA.nondegenerate_chars
```

Next, we'll compute a table of the per-sequence kmer counts for all kmers in `W` for all sequences in our reference database. We'll also store the taxonomic label of each of our reference sequences at our specified taxonomic level. We can store this information in a pandas `DataFrame`, and then view the first 25 rows of that table.

```python
>>> def get_taxon_at_level(taxon, level):
...     taxon = [l.strip() for l in taxon.split(';')]
...     return '; '.join(taxon[:level])
...
>>> W = compute_W(alphabet, k)
...
>>> per_sequence_kmer_counts = []
>>> for reference_sequence in reference_db:
...     taxon = get_taxon_at_level(reference_sequence.metadata['taxonomy'], taxonomic_level)
...     kmer_counts = dict.fromkeys(W, 0)
...     kmer_counts.update(reference_sequence.kmer_frequencies(k=k))
...     per_sequence_kmer_counts.append(pd.Series(kmer_counts, name=taxon))
...
>>> per_sequence_kmer_counts = pd.DataFrame(data=per_sequence_kmer_counts).fillna(0).T
>>> per_sequence_kmer_counts[:25]
```

With this information, we'll next compute our "kmer probability table" (EXISTING NAME FOR THIS?). The content of this table will be the probability of observing each kmer in W given a taxon. This is computed based on a few values:

$N$ : The total number of sequences in the training set.

$n(w_i)$ : The number of total sequences containing kmer _i_.

$P_i$ : The probability of observing kmer _i_. Initially it might seem as though this would be computed as $n(w_i) / N$, but this neglects the possibility of that a kmer observed in a query sequence might not be represented in our reference database, so a small pseudocount is added to the numerator and denomenator.

$P(w_i | taxon)$ : The probability of observing a kmer given a taxon. Again, it would seem that this would be computed as the proportion of sequences in the taxon containing the kmer, but this would neglect that we'll likely observe kmers in our query sequences that are not represented in our reference database. As pseudocount is therefore added again to the numerator and denominator. This time the pseudocount in the numerator is scaled by how frequent the kmer is in the reference database as a whole: specifically, it is $P_i$.

Our "kmer probability table" is $P(w_i | taxon)$ computed for all kmers in W and all taxa represented in our reference database. We'll compute that and again look at the first 25 rows.

```python
>>> def compute_kmer_probability_table(per_sequence_kmer_counts):
...     N = len(per_sequence_kmer_counts) # number of training sequences
...
...     # number of sequences containing kmer wi
...     n_wi = per_sequence_kmer_counts.astype(bool).sum(axis=1)
...     n_wi.name = 'n(w_i)'
...
...     # probabilities of observing each kmer
...     Pi = (n_wi + 0.5) / (N + 1)
...     Pi.name = 'P_i'
...
...     # number of times each taxon appears in training set
...     taxon_counts = collections.Counter(per_sequence_kmer_counts.columns)
...     n_taxon_members_containing_kmer = per_sequence_kmer_counts.astype(bool).groupby(level=0, axis=1).sum()
...
...     # probabilities of observing each kmer in each taxon
...     p_wi_t = []
...     for taxon, count in taxon_counts.items():
...         p_wi_t.append(pd.Series((n_taxon_members_containing_kmer[taxon] + Pi) / (count + 1), name=taxon))
...
...     return pd.DataFrame(p_wi_t).T
```

```python
>>> kmer_probability_table = compute_kmer_probability_table(per_sequence_kmer_counts)
```

```python
>>> kmer_probability_table[:25]
```

With our kmer probability table we are now ready to classify unknown sequences. We'll begin by defining some query sequences. We'll pull these at random from our reference sequences, which means that some of the query sequences will be represented in our reference database and some won't be. This is the sitatuation that is typically encountered in practice. To simulate real-world 16S taxonomy classification tasks, we'll also trim out 200 bases of our reference sequences since (as of this writing) we typically don't obtain full-length 16S sequences from a DNA sequencing instrument.

```python
>>> from iab.algorithms import load_taxonomy_query_sequences
...
>>> %psource load_taxonomy_query_sequences
```

```python
>>> import random
...
>>> queries = load_taxonomy_query_sequences()
>>> queries = random.sample(queries, k=50)
```

```python
>>> queries[0]
```

For a given query sequence, its taxonomy will be classified as follows. First, the set of all kmers will be extracted from the sequence. This is referred to as $V$. Then, for all taxa in the kmer probability table, the probability of observing the query sequence will be computed given that taxon: $P(query | taxon)$. This is computed as the product of all its kmer probabilities for the given taxon. (It should be clear based on this formula why it was necessary to add pseudocounts when computing our kmer probability table - if not, kmer probabilities of zero would result in a zero probability of the sequence being derived from that taxon at this step.)

After computing $P(query | taxon)$ for all taxa, the taxonomy assignment return is simply the one achieving the maximum probability. Here we'll classify a sequence and look at the resulting taxonomy assignment.

```python
>>> def classify_V(V, kmer_probability_table):
...     P_S_t = [] # probability of the sequence given the taxon
...     for taxon in kmer_probability_table:
...         kmer_probabilities = kmer_probability_table[taxon]
...         probability = 1.0
...         for v_i in V:
...             probability *= kmer_probabilities[v_i]
...         P_S_t.append((probability, taxon))
...     return max(P_S_t)[1], V
...
>>> def classify_sequence(query_sequence, kmer_probability_table, k):
...     V = list(map(str, query_sequence.iter_kmers(k=k)))
...     return classify_V(V, kmer_probability_table)
```

```python
>>> taxon_assignment, V = classify_sequence(queries[0], kmer_probability_table, k)
>>> print(taxon_assignment)
```

Since we know the actual taxonomy assignment for this sequence, we can look that up in our reference database. Was your assignment correct? Try this with a few query sequences and keep track of how many times the classifier achieved the correct assignment.

```python
>>> get_taxon_at_level(reference_taxonomy[queries[0].metadata['id']], taxonomic_level)
```

Because the query and reference sequences that were working with were randomly selected from the full reference database, each time you run this notebook you should observe different results. Chances are however that if you run the above steps multiple times you'll get the wrong taxonomy assignment at least some of the time. Up to this point, we've left out an important piece of information: how confident should we be in our assignment, or in other words, how dependent is our taxonomy assignment on our specific query? If there were slight differences in our query (e.g., because we observed a very closely related organism, such as one of the same species but a different strain, or because we sequenced a different region of the 16S sequence) would we obtain the same taxonomy assignment? If so, we should have higher confidence in our assignment. If not, we should have lower confidence in our assignment.

We can quantify confidence using an approach called bootstrapping. With a bootstrap approach, we'll get our taxonomy assignment as we did above, but then for some user-specified number of times, we'll create random subsets of V sampled with replacement (DEFINE THIS). We'll then assign taxonomy to each random subset of V, and count the number of times the resulting taxonomy assignment is the same that we achieved when assigning taxonomy to V. The count divided by the number of iterations we've chosen to run will be our confidence value. If the assignments are often the same we'll have a high confidence value. If the assignments are often different, we'll have a low confidence value.

Let's now assign taxonomy and compute a confidence for that assignment.

```python
>>> def classify_sequence_with_confidence(sequence, kmer_probability_table, k,
...                                       confidence_iterations=100):
...     taxon, V = classify_sequence(sequence, kmer_probability_table, k)
...
...     count_same_taxon = 0
...     subsample_size = int(len(V) * 0.1)
...     for i in range(confidence_iterations):
...         subsample_V = np.random.choice(V, subsample_size, replace=True)
...         subsample_taxon, _ = classify_V(subsample_V, kmer_probability_table)
...         if taxon == subsample_taxon:
...             count_same_taxon += 1
...     confidence = count_same_taxon / confidence_iterations
...
...     return (taxon, confidence)
```

```python
>>> taxon_assignment, confidence = classify_sequence_with_confidence(queries[0], kmer_probability_table, k)
>>> print(taxon_assignment)
>>> print(confidence)
```

How did the computed confidence compare to the accuracy taxonomy assignment?

We don't have an _a priori_ idea for what good versus bad confidence scores are, but we can use our reference database to explore that. We might want this information so we can come up with a confidence threshold, above which we would accept a taxonomy assignment and below which we might reject it. To explore this, let's compute taxonomy assignments and confidence for all of our query sequences and then see what the distributions of confidence scores look like for correct assignments and incorrect assignments.

```python
>>> correct_assignment_confidences = []
>>> incorrect_assignment_confidences = []
>>> summary = []
...
>>> for query in queries:
...     predicted_taxonomy, confidence = classify_sequence_with_confidence(query, kmer_probability_table, k)
...     actual_taxonomy = get_taxon_at_level(reference_taxonomy[query.metadata['id']], taxonomic_level)
...     if actual_taxonomy == predicted_taxonomy:
...         correct_assignment_confidences.append(confidence)
...     else:
...         incorrect_assignment_confidences.append(confidence)
...
...     summary.append([predicted_taxonomy, actual_taxonomy, confidence])
>>> summary = pd.DataFrame(summary, columns=['Predicted taxonomy', 'Actual taxonomy', 'Confidence'])
```

```python
>>> import seaborn as sns
...
>>> ax = sns.boxplot(data=[correct_assignment_confidences, incorrect_assignment_confidences])
>>> ax = sns.swarmplot(data=[correct_assignment_confidences, incorrect_assignment_confidences], color="black")
>>> _ = ax.set_xticklabels(['Correct assignments', 'Incorrect assignments'])
>>> _ = ax.set_ylabel('Confidence')
```

What does this plot tell you about how well setting a confidence threshold is likely to work? If you never wanted to reject a correct assignment, how often would you accept an incorrect assignment? If you never wanted to accept an incorrect assignment, how often would you reject a correct assignment?

```python
>>> summary # maybe explore whether certain taxa are more frequently wrong than others...
```

### Random Forest classifiers <link src="N7CyaN"/>

Coming soon...

## Defining a dimensionality reduction problem <link src="Y0TtkW"/>

[This content](alias://b1cdbe) will be adapted and ported here.
