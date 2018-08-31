# Glossary <link src="wDGwYI"/>

## Pairwise alignment (noun) <link src="bSQQ5s"/>

A hypothesis about which bases or amino acids in two biological sequences are derived from a common ancestral base or amino acid. By definition, the *aligned sequences* will be of equal length with gaps (usually denoted with ``-``, or ``.`` for terminal gaps) indicating hypothesized insertion deletion events. A pairwise alignment may be represented as follows:

```
ACC---GTAC
CCCATCGTAG
```

## kmer (noun) <link src="C7hMX5"/>

A kmer is simply a word (or list of adjacent characters) in a sequence of length k. For example, the overlapping kmers in the sequence ``ACCGTGACCAGTTACCAGTTTGACCAA`` are as follows:

```python
>>> import skbio
>>> skbio.DNA('ACCGTGACCAGTTACCAGTTTGACCAA').kmer_frequencies(k=5, overlap=True)
```
