# Machine learning in bioinformatics <link src="7SFnTr"/>

## Supervised v unsupervised classification <link src="b05Cya"/>

## Training data, test data, and cross validation <link src="QIewra"/>

## scikit-learn <link src="ifjF4O"/>

In this chapter we'll implement several machine learning classifiers so we can gain an in-depth understanding of how they work. In practice though, there are many mature machine learning libraries that you'd want to use. [scikit-learn](http://scikit-learn.org/) is a popular and well-documented Python library for machine learning which many bioinformatics researchers and software developers use in their work.

## Defining the problem <link src="2R6CTy"/>

We'll explore machine learning classifiers in the context of a familiar topic: taxonomic classification of 16S rRNA sequences. We previously explored this problem in [Sequence Homology Searching](alias://d22e6b), so it is likely worth spending a few minutes skimming that chapter if it's not fresh in your mind.

Briefly, the problem that we are going to address here is as follows. We have a query sequence ($q_i$) which is not taxonomically annotated (meaning we don't know the taxonomy of the organism whose genome it is found in), and a reference database ($R$) of taxonomically annotated sequences ($r_1, r_2, r_3, ... r_n$). We want to infer a taxonomic annotation for $q_i$. We'll again work with the [Greengenes](http://greengenes.secondgenome.com/) database, which we'll access using [QIIME default reference project](https://github.com/biocore/qiime-default-reference). Greengenes is a database of 16S rRNA gene sequences. (This should all sound very familiar - if not, I again suggest that you review [Sequence Homology Searching](alias://d22e6b).)

This time, instead of using sequence alignment to identify the most likely taxonomic origin of a sequence, we'll train classifiers (i.e., build kmer-based models) of the 16S sequences of taxa in our reference database. We'll then run our query sequences through those models to identify the most likely taxonomic origin of each query sequence. Since we know the taxonomic origin of our query sequences in this case, we can test our classifiers by seeing how often they return the known taxonomy assignment. If our training and testing approaches are well-designed, the performance on our tests will inform us of how accurate we can expect our classifier to be on data where the actual taxonomic origin is unknown.

Let's jump in...

## Naive Bayes classifiers <link src="H8vYPu"/>

The first classifier we'll explore is the popular and relatively simple Naive Bayes classifier. This classifier uses Bayes Theorem to determine the most likely label for an unknown input based on a probabilistic model is has constructed from training data. (_The preceding text needs work._) The model that is constructed is based on user-defined features of the sequences. The most commonly used features for sequence classification tasks such as this is overlapping [kmers](alias://C7hMX5).

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
>>> import qiime_default_reference as qdr
...
>>> # Load the taxonomic data
... reference_taxonomy = {}
>>> for e in open(qdr.get_reference_taxonomy()):
...     seq_id, seq_tax = e.strip().split('\t')
...     reference_taxonomy[seq_id] = seq_tax
...
>>> # Load the reference sequences, and associate the taxonomic annotation with
... # each as metadata
... reference_db = []
>>> for e in skbio.io.read(qdr.get_reference_sequences(), format='fasta', constructor=skbio.DNA):
...     seq_tax = reference_taxonomy[e.metadata['id']]
...     e.metadata['taxonomy'] = seq_tax
...     reference_db.append(e)
...
>>> print("%s sequences were loaded from the reference database." % len(reference_db))
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
>>> alphabet = ['A', 'C', 'G', 'T']
>>> k = 2
...
>>> def compute_W(alphabet, k):
>>>  return set(map(''.join, itertools.product(alphabet, repeat=k)))
...
>>> W = compute_W(alphabet, k)
>>> len(W)
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

```python
>>> def get_taxon_at_level(taxon, level):
...     taxon = [l.strip() for l in taxon.split(';')]
...     return '; '.join(taxon[:level])
```

```python
>>> import pandas as pd
...
>>> taxonomic_level = 2
>>> k = 7
>>> alphabet = ['A', 'C', 'G', 'T', 'N']
...
>>> W = compute_W(alphabet, k)
...
>>> per_sequence_kmer_counts = []
>>> for reference_sequnence in reference_db:
...     taxon = get_taxon_at_level(reference_sequnence.metadata['taxonomy'], taxonomic_level)
...     kmer_counts = dict.fromkeys(W, 0)
...     kmer_counts.update(reference_sequnence.kmer_frequencies(k=k))
...     per_sequence_kmer_counts.append(pd.Series(kmer_counts, name=taxon))
>>> # if W is passed as the index can I drop the fromkeys step above?
... per_sequence_kmer_counts = pd.DataFrame(data=per_sequence_kmer_counts).fillna(0).T
>>> per_sequence_kmer_counts
```

```python
>>> def compute_training_table(per_sequence_kmer_counts):
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
...     # probabilities of observing each word in each taxon
...     p_wi_t = []
...     for taxon, count in taxon_counts.items():
...         p_wi_t.append(pd.Series((n_taxon_members_containing_kmer[taxon] + Pi) / (count + 1), name=taxon))
...
...     return pd.DataFrame(p_wi_t).T
```

```python
>>> training_table = compute_training_table(per_sequence_kmer_counts)
```

```python
>>> training_table
```

```python
>>> def classify_kmers(V, training_table):
...     P_S_t = [] # probability of the sequence given the taxon
...     for taxon in training_table:
...         kmer_probabilities = training_table[taxon]
...         probability = 1.0
...         for v_i in V:
...             probability *= kmer_probabilities[v_i]
...         P_S_t.append((probability, taxon))
...     return(max(P_S_t)[1])
...
>>> def classify_sequence(query, training_table, k,
...                       compute_confidence=True, confidence_iterations=100):
...     V = list(map(str, query.iter_kmers(k=k)))
...     taxon = classify_kmers(V, training_table)
...
...     if compute_confidence:
...         count_same_taxon = 0
...         subsample_size = int(len(V) * 0.1)
...         for i in range(confidence_iterations):
...             subsample_V = np.random.choice(V, subsample_size, replace=True)
...             subsample_taxon = classify_kmers(subsample_V, training_table)
...             if taxon == subsample_taxon:
...                 count_same_taxon += 1
...         confidence = count_same_taxon / confidence_iterations
...
...     else:
...         confidence = -1.
...
...     return(taxon, confidence)
```

```python
>>> queries = []
>>> for e in skbio.io.read(qdr.get_reference_sequences(), format='fasta', constructor=skbio.DNA):
...     e = e[100:300]
...     queries.append(e)
>>> # can't figure out why np.random.choice isn't working here...
... np.random.shuffle(queries)
>>> queries = queries[:50]
```

```python
>>> queries[0]
```

```python
>>> correct_assignment_confidences = []
>>> incorrect_assignment_confidences = []
>>> summary = []
...
>>> for query in queries:
...     predicted_taxonomy, confidence = classify_sequence(query, training_table, k)
...     actual_taxonomy = get_taxon_at_level(reference_taxonomy[query.metadata['id']], taxonomic_level)
...     if actual_taxonomy == predicted_taxonomy:
...         correct_assignment_confidences.append(confidence)
...     else:
...         incorrect_assignment_confidences.append(confidence)
...
...     summary.append([predicted_taxonomy, actual_taxonomy, confidence])
>>> summary = pd.DataFrame(summary, columns=['Predicted taxonomy', 'Actual taxonomy', 'Confidence'])
>>> print(np.median(correct_assignment_confidences))
>>> print(np.median(incorrect_assignment_confidences))
```

```python
>>> summary
```

## Random Forest classifiers <link src="N7CyaN"/>

## Neural networks and "deep learning" <link src="DgnnyS"/>
