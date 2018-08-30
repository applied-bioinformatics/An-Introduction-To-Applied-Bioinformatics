# Machine learning in bioinformatics <link src="7SFnTr"/>

## Supervised v unsupervised classification

## Training data, test data, and cross validation

## Naive Bayes classifiers

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

```python
>>> reference_db = np.random.choice(reference_db, 500, replace=False)
>>> print("%s sequences are present in the subsampled database." % len(reference_db))
```

```python
>>> alphabet = ['A', 'C', 'G', 'T', 'N']
>>> word_length = 2
...
>>> def compute_W(alphabet, word_length):
>>>  return set(map(''.join, itertools.product(alphabet, repeat=word_length)))
...
>>> W = compute_W(alphabet, word_length)
>>> len(W)
```

```python
>>> pd.Series(reference_db[0].kmer_frequencies(k=word_length), name=reference_db[0].metadata['id'])
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
>>> word_length = 7
>>> alphabet = ['A', 'C', 'G', 'T', 'N']
...
>>> W = compute_W(alphabet, word_length)
...
>>> per_sequence_word_counts = []
>>> for reference_sequnence in reference_db:
...     taxon = get_taxon_at_level(reference_sequnence.metadata['taxonomy'], taxonomic_level)
...     word_counts = dict.fromkeys(W, 0)
...     word_counts.update(reference_sequnence.kmer_frequencies(k=word_length))
...     per_sequence_word_counts.append(pd.Series(word_counts, name=taxon))
>>> # if W is passed as the index can I drop the fromkeys step above?
... per_sequence_word_counts = pd.DataFrame(data=per_sequence_word_counts).fillna(0).T
>>> per_sequence_word_counts
```

```python
>>> def compute_training_table(per_sequence_word_counts):
...     N = len(per_sequence_word_counts) # number of training sequences
...
...     # number of sequences containing word wi
...     n_wi = per_sequence_word_counts.astype(bool).sum(axis=1)
...     n_wi.name = 'n(w_i)'
...
...     # probabilities of observing each word
...     Pi = (n_wi + 0.5) / (N + 1)
...     Pi.name = 'P_i'
...
...     # number of times each taxon appears in training set
...     taxon_counts = collections.Counter(per_sequence_word_counts.columns)
...     n_taxon_members_containing_word = per_sequence_word_counts.astype(bool).groupby(level=0, axis=1).sum()
...
...     # probabilities of observing each word in each taxon
...     p_wi_t = []
...     for taxon, count in taxon_counts.items():
...         p_wi_t.append(pd.Series((n_taxon_members_containing_word[taxon] + Pi) / (count + 1), name=taxon))
...
...     return pd.DataFrame(p_wi_t).T
```

```python
>>> training_table = compute_training_table(per_sequence_word_counts)
```

```python
>>> training_table
```

```python
>>> def classify_words(V, training_table):
...     P_S_t = [] # probability of the sequence given the taxon
...     for taxon in training_table:
...         word_probabilities = training_table[taxon]
...         probability = 1.0
...         for v_i in V:
...             probability *= word_probabilities[v_i]
...         P_S_t.append((probability, taxon))
...     return(max(P_S_t)[1])
...
>>> def classify_sequence(query, training_table, word_length,
...                       compute_confidence=True, confidence_iterations=100):
...     V = list(map(str, query.iter_kmers(k=word_length)))
...     taxon = classify_words(V, training_table)
...
...     if compute_confidence:
...         count_same_taxon = 0
...         subsample_size = int(len(V) * 0.1)
...         for i in range(confidence_iterations):
...             subsample_V = np.random.choice(V, subsample_size, replace=True)
...             subsample_taxon = classify_words(subsample_V, training_table)
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
...     predicted_taxonomy, confidence = classify_sequence(query, training_table, word_length)
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

## Random Forest classifiers

## Neural networks and "deep learning"
