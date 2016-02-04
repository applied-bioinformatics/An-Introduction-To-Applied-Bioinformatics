
# Sequence homology searching <link src='d22e6b'/>

In this chapter we'll talk about using pairwise alignment to search databases of biological sequences with the goal of identifying sequence homology. We [previously defined homology](alias://e63a4f) between a pair of sequences to mean that those sequences are derived from a common ancestral sequence. Homology searching is an essential part of making inferences about where a biological sequence came from, and/or what it does. In most cases, if you have an unannotated biological sequence, such as the following protein sequence, it's very hard (really, impossible) to know what it is without more information.

```
> mystery-sequence1
PQITLWQRPLVTIRIGGQLKEALLDTGADDTVLEEMNLPGKWKPKMIGGIGGFIKVRQYDQIPVEIAHKAIGTVLVGPTPVNIIGRNLLTQIGATLNF
```

What a researcher will often do is search this sequence, their *query*, against some *reference database* of annotated sequences to learn what function the sequence performs (if the reference database contains functional annotation of sequences) and/or what organisms are likely to encode this sequence in their genome (if the reference database contains taxonomic annotation of sequences).

Whose genome is the above sequence encoded in? What is its function? Take a minute now to answer these questions using the [Protein BLAST homology search tool on the NCBI website](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp).

In the context of database searching, a query sequence and a reference sequence that we hypothesize to be homologous can be identical to one another, or they can differ as a result of mutation events. When sequences differ, we're often then interested in how much they differ, or their pairwise similarity, which can help us identify the most closely related of several homologs in the reference database. There is an important distinction in the terms *homology* and *similarity*: homology is a discrete variable, and similarity is a continuous variable. A pair of biological sequences either *are* or *are not* derived from a common ancestor, but they can be more or less similar to each other. Saying that two sequences are 80% homologous doesn't make sense. What people generally mean when they say this is that two sequences are 80% similar, and as a result they are hypothesizing homology between the sequences.

## Defining the problem

As mentioned above, if we want to perform a homology search we'll have one or more *query sequences*, and for each we want to know which sequence(s) in a reference database it is most similar to.

Sequence homology searching can be implemented in a few ways. In this chapter, we'll use the local alignment function that we worked with in [the Pairwise Alignment chapter](alias://a76822), ``local_pairwise_align_ssw``, run it many times to search one *query* sequence against many *reference* sequences, and investigate the highest scoring alignment(s) to identify the best database match. Remember that you can always get help with a function by passing it as an argument to ``help``:

```python
>>> from skbio.alignment import local_pairwise_align_ssw
>>> help(local_pairwise_align_ssw)
Help on function local_pairwise_align_ssw in module skbio.alignment._pairwise:

local_pairwise_align_ssw(sequence1, sequence2, **kwargs)
    Align query and target sequences with Striped Smith-Waterman.
    
    State: Experimental as of 0.4.0.
    
    Parameters
    ----------
    sequence1 : DNA, RNA, or Protein
        The first unaligned sequence
    sequence2 : DNA, RNA, or Protein
        The second unaligned sequence
    
    Returns
    -------
    tuple
        ``TabularMSA`` object containing the aligned sequences, alignment score
        (float), and start/end positions of each input sequence (iterable
        of two-item tuples). Note that start/end positions are indexes into the
        unaligned sequences.
    
    Notes
    -----
    This is a wrapper for the SSW package [1]_.
    
    For a complete list of optional keyword-arguments that can be provided,
    see ``skbio.alignment.StripedSmithWaterman``.
    
    The following kwargs will not have any effect: `suppress_sequences`,
    `zero_index`, and `protein`
    
    If an alignment does not meet a provided filter, `None` will be returned.
    
    References
    ----------
    .. [1] Zhao, Mengyao, Wan-Ping Lee, Erik P. Garrison, & Gabor T.
       Marth. "SSW Library: An SIMD Smith-Waterman C/C++ Library for
       Applications". PLOS ONE (2013). Web. 11 July 2014.
       http://www.plosone.org/article/info:doi/10.1371/journal.pone.0082138
    
    See Also
    --------
    skbio.alignment.StripedSmithWaterman
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
Populating the interactive namespace from numpy and matplotlib
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
99,322 sequences were loaded from the reference database.
```

Next, we'll just inspect a couple of the sequences we loaded. Notice how the specificity of our taxonomic annotations (i.e., how many taxonomic levels are annotated and unknown) differs for different sequences.

```python
>>> reference_db[0]
DNA
-----------------------------------------------------------------------
Metadata:
    'description': ''
    'id': '1111883'
    'taxonomy': 'k__Bacteria; p__Gemmatimonadetes; c__Gemm-1; o__; f__;
                 g__; s__'
Stats:
    length: 1428
    has gaps: False
    has degenerates: False
    has non-degenerates: True
    GC-content: 61.90%
-----------------------------------------------------------------------
0    GCTGGCGGCG TGCCTAACAC ATGTAAGTCG AACGGGACTG GGGGCAACTC CAGTTCAGTG
60   GCAGACGGGT GCGTAACACG TGAGCAACTT GTCCGACGGC GGGGGATAGC CGGCCCAACG
...
1320 GCCGCGGTGA ATACGTTCCC GGGCCTTGTA CACACCGCCC GTCACGCCAT GGAAGCCGGA
1380 GGGACCCGAA ACCGGTGGGC CAACCGCAAG GGGGCAGCCG TCTAAGGT
```

```python
>>> reference_db[-1]
DNA
----------------------------------------------------------------------
Metadata:
    'description': ''
    'id': '4483258'
    'taxonomy': 'k__Archaea; p__Crenarchaeota; c__Thermoprotei;
                 o__Thermoproteales; f__Thermoproteaceae; g__; s__'
Stats:
    length: 2123
    has gaps: False
    has degenerates: False
    has non-degenerates: True
    GC-content: 58.36%
----------------------------------------------------------------------
0    CTGGTTGATC CTGCCGGACC CGACCGCTAT CGGGGTGGGG CTTAGCCATG CGAGTCAAGC
60   GCCCCAGGGA CCCGCTGGGG TGCGGCGCAC GGCTCAGTAA CACGTGGCCA ACCTACCCTC
...
2040 ATAATCTCCT TATTGTCTGA TCCTTATGCA TTTTCCTTTG GCCCATCCCG TGAATACGCG
2100 CGGTGAATAC GTCCCTGCCC CTT
```

For the sake of runtime, we're going to work through this chapter using a random sample of sequences from this database. Here we'll use Python's [random module](https://docs.python.org/3/library/random.html) to select sequences at random.

```python
>>> import random
...
>>> reference_db = random.sample(reference_db, k=5000)
>>> print("%s sequences are present in the subsampled database." % locale.format("%d", len(reference_db), grouping=True))
5,000 sequences are present in the subsampled database.
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
DNA
---------------------------------------------------------------------
Metadata:
    'description': ''
    'id': '1129362'
Stats:
    length: 200
    has gaps: False
    has degenerates: False
    has non-degenerates: True
    GC-content: 55.00%
---------------------------------------------------------------------
0   TGACACGTGG GCAACCTGCC TGTAAGATCG GGATAACTCC GGGAAACCGG GGCTAATACC
60  GGATAATACT TTTTCTCGCG TGGGGAGAAG TTGAAAGATG GCTTCGGCTA TCGCTTACAG
120 ATGGGCCCGC GGCGCATTAG CTAGTTGGTA GGGTAACGGC CTACCAAGGC AACGATGCGT
180 AGCCGACCTG AGAGGGTGAT
```

```python
>>> queries[-1]
DNA
---------------------------------------------------------------------
Metadata:
    'description': ''
    'id': '241186'
Stats:
    length: 200
    has gaps: False
    has degenerates: False
    has non-degenerates: True
    GC-content: 58.00%
---------------------------------------------------------------------
0   TAACGGCCTA CCAAGGCGAC GATCCGTAGC TGGTTTGAGA GGACGACCAG CCACACTGGG
60  ACTGAGACAC GGCCCAGACT CCTACGGGAG GCAGCAGTGG GGAATTTTGG ACAATGGGGG
120 AAACCCTGAT CCAGCCATCC CGCGTGTGCG ATGAAGGCCT TCGGGTTGTA AAGCACTTTT
180 GGCAGGAAAG AAACGTCGCG
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
[0;32mdef[0m [0mlocal_alignment_search[0m[0;34m([0m[0mqueries[0m[0;34m,[0m [0mreference_db[0m[0;34m,[0m [0mn[0m[0;34m=[0m[0;36m5[0m[0;34m,[0m[0;34m[0m
[0;34m[0m                           [0maligner[0m[0;34m=[0m[0mlocal_pairwise_align_ssw[0m[0;34m)[0m[0;34m:[0m[0;34m[0m
[0;34m[0m    [0mresults[0m [0;34m=[0m [0;34m[[0m[0;34m][0m[0;34m[0m
[0;34m[0m    [0mindices[0m [0;34m=[0m [0;34m[[0m[0;34m][0m[0;34m[0m
[0;34m[0m    [0;32mfor[0m [0mq[0m [0;32min[0m [0mqueries[0m[0;34m:[0m[0;34m[0m
[0;34m[0m        [0;31m# first we'll compute all of the alignments and their associated scores[0m[0;34m[0m
[0;34m[0m        [0mhits[0m [0;34m=[0m [0;34m[[0m[0;34m][0m[0;34m[0m
[0;34m[0m        [0;32mfor[0m [0mr[0m [0;32min[0m [0mreference_db[0m[0;34m:[0m[0;34m[0m
[0;34m[0m            [0maln[0m[0;34m,[0m [0mscore[0m[0;34m,[0m [0m_[0m [0;34m=[0m [0maligner[0m[0;34m([0m[0mq[0m[0;34m,[0m [0mr[0m[0;34m)[0m[0;34m[0m
[0;34m[0m            [0mhits[0m[0;34m.[0m[0mappend[0m[0;34m([0m[0;34m[[0m[0mr[0m[0;34m.[0m[0mmetadata[0m[0;34m[[0m[0;34m'id'[0m[0;34m][0m[0;34m,[0m [0mscore[0m[0;34m,[0m [0maln[0m[0;34m,[0m[0;34m[0m
[0;34m[0m                         [0mr[0m[0;34m.[0m[0mmetadata[0m[0;34m[[0m[0;34m'taxonomy'[0m[0;34m][0m[0;34m][0m[0;34m)[0m[0;34m[0m
[0;34m[0m        [0;31m# then we reverse-sort them by score, and return the n highest[0m[0;34m[0m
[0;34m[0m        [0;31m# scoring alignments (this needs to be updated so we only[0m[0;34m[0m
[0;34m[0m        [0;31m# ever keep track of the n highest scoring alignments)[0m[0;34m[0m
[0;34m[0m        [0mbest_hits[0m [0;34m=[0m [0msorted[0m[0;34m([0m[0mhits[0m[0;34m,[0m [0mkey[0m[0;34m=[0m[0;32mlambda[0m [0me[0m[0;34m:[0m [0me[0m[0;34m[[0m[0;36m1[0m[0;34m][0m[0;34m,[0m [0mreverse[0m[0;34m=[0m[0;32mTrue[0m[0;34m)[0m[0;34m[[0m[0;34m:[0m[0mn[0m[0;34m][0m[0;34m[0m
[0;34m[0m        [0;31m# then we compile and track some information about the n best hits[0m[0;34m[0m
[0;34m[0m        [0;32mfor[0m [0mr_id[0m[0;34m,[0m [0mscore[0m[0;34m,[0m [0maln[0m[0;34m,[0m [0mr_tax[0m [0;32min[0m [0mbest_hits[0m[0;34m:[0m[0;34m[0m
[0;34m[0m            [0mpercent_similarity[0m [0;34m=[0m [0;34m([0m[0;36m100[0m [0;34m*[0m [0;34m([0m[0;36m1.[0m [0;34m-[0m [0maln[0m[0;34m[[0m[0;36m0[0m[0;34m][0m[0;34m.[0m[0mdistance[0m[0;34m([0m[0maln[0m[0;34m[[0m[0;36m1[0m[0;34m][0m[0;34m)[0m[0;34m)[0m[0;34m)[0m[0;34m[0m
[0;34m[0m            [0maln_length[0m [0;34m=[0m [0maln[0m[0;34m.[0m[0mshape[0m[0;34m[[0m[0;36m1[0m[0;34m][0m[0;34m[0m
[0;34m[0m            [0mindices[0m[0;34m.[0m[0mappend[0m[0;34m([0m[0;34m([0m[0mq[0m[0;34m.[0m[0mmetadata[0m[0;34m[[0m[0;34m'id'[0m[0;34m][0m[0;34m,[0m [0mr_id[0m[0;34m)[0m[0;34m)[0m[0;34m[0m
[0;34m[0m            [0mresults[0m[0;34m.[0m[0mappend[0m[0;34m([0m[0;34m([0m[0mr_tax[0m[0;34m,[0m [0mpercent_similarity[0m[0;34m,[0m[0;34m[0m
[0;34m[0m                            [0maln_length[0m[0;34m,[0m [0mscore[0m[0;34m)[0m[0;34m)[0m[0;34m[0m
[0;34m[0m    [0mindex[0m [0;34m=[0m [0mpd[0m[0;34m.[0m[0mMultiIndex[0m[0;34m.[0m[0mfrom_tuples[0m[0;34m([0m[0mindices[0m[0;34m,[0m [0mnames[0m[0;34m=[0m[0;34m[[0m[0;34m'query'[0m[0;34m,[0m [0;34m'reference'[0m[0;34m][0m[0;34m)[0m[0;34m[0m
[0;34m[0m    [0mcolumns[0m [0;34m=[0m [0;34m[[0m[0;34m'reference taxonomy'[0m[0;34m,[0m [0;34m'percent similarity'[0m[0;34m,[0m[0;34m[0m
[0;34m[0m               [0;34m'alignment length'[0m[0;34m,[0m [0;34m'score'[0m[0;34m][0m[0;34m[0m
[0;34m[0m    [0mresults[0m [0;34m=[0m [0mpd[0m[0;34m.[0m[0mDataFrame[0m[0;34m([0m[0mresults[0m[0;34m,[0m [0mindex[0m[0;34m=[0m[0mindex[0m[0;34m,[0m [0mcolumns[0m[0;34m=[0m[0mcolumns[0m[0;34m)[0m[0;34m[0m
[0;34m[0m    [0;32mreturn[0m [0mresults[0m[0;34m[0m[0m
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
Runtime: 8.7565 sec per query
                                                  reference taxonomy  \
query   reference                                                      
316642  326804     k__Bacteria; p__Firmicutes; c__Clostridia; o__...   
        311096     k__Bacteria; p__Firmicutes; c__Clostridia; o__...   
        323155     k__Bacteria; p__Firmicutes; c__Clostridia; o__...   
        251395     k__Bacteria; p__Firmicutes; c__Clostridia; o__...   
        333675     k__Bacteria; p__Firmicutes; c__Clostridia; o__...   
3900488 4384293    k__Bacteria; p__Firmicutes; c__Clostridia; o__...   
        1689987    k__Bacteria; p__Firmicutes; c__Clostridia; o__...   
        4311536    k__Bacteria; p__Firmicutes; c__Clostridia; o__...   
        152014     k__Bacteria; p__Firmicutes; c__Clostridia; o__...   
        333471     k__Bacteria; p__Firmicutes; c__Clostridia; o__...   
4478817 719879     k__Bacteria; p__Bacteroidetes; c__[Saprospirae...   
        1104589    k__Bacteria; p__Bacteroidetes; c__[Saprospirae...   
        4384311    k__Bacteria; p__Bacteroidetes; c__[Saprospirae...   
        808539     k__Bacteria; p__Bacteroidetes; c__[Saprospirae...   
        208595     k__Bacteria; p__Bacteroidetes; c__[Saprospirae...   
254871  575373     k__Bacteria; p__Bacteroidetes; c__Sphingobacte...   
        103997     k__Bacteria; p__Bacteroidetes; c__Cytophagia; ...   
        4302951    k__Bacteria; p__Bacteroidetes; c__Cytophagia; ...   
        4342029    k__Bacteria; p__Bacteroidetes; c__Cytophagia; ...   
        826061     k__Bacteria; p__Bacteroidetes; c__Cytophagia; ...   

                   percent similarity  alignment length  score  
query   reference                                               
316642  326804              90.594059               202    330  
        311096              89.447236               199    294  
        323155              88.944724               199    289  
        251395              88.888889               198    286  
        333675              87.939698               199    279  
3900488 4384293            100.000000               200    400  
        1689987             92.574257               202    336  
        4311536             81.407035               199    236  
        152014              81.909548               199    233  
        333471              78.500000               200    230  
4478817 719879              72.857143               210    217  
        1104589             79.899497               199    210  
        4384311             77.832512               203    198  
        808539              76.470588               204    193  
        208595              77.889447               199    179  
254871  575373              82.587065               201    225  
        103997              82.000000               200    220  
        4302951             80.295567               203    220  
        4342029             78.325123               203    210  
        826061              79.207921               202    202
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
Closest taxonomies for query 316642 (in order):
  k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__
  k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__
  k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__
  k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__
  k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__

Closest taxonomies for query 3900488 (in order):
  k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__
  k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__
  k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__
  k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__
  k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__

Closest taxonomies for query 4478817 (in order):
  k__Bacteria; p__Bacteroidetes; c__[Saprospirae]; o__[Saprospirales]; f__Saprospiraceae; g__; s__
  k__Bacteria; p__Bacteroidetes; c__[Saprospirae]; o__[Saprospirales]; f__Saprospiraceae; g__; s__
  k__Bacteria; p__Bacteroidetes; c__[Saprospirae]; o__[Saprospirales]; f__Saprospiraceae; g__; s__
  k__Bacteria; p__Bacteroidetes; c__[Saprospirae]; o__[Saprospirales]; f__Saprospiraceae; g__; s__
  k__Bacteria; p__Bacteroidetes; c__[Saprospirae]; o__[Saprospirales]; f__Saprospiraceae; g__; s__

Closest taxonomies for query 254871 (in order):
  k__Bacteria; p__Bacteroidetes; c__Sphingobacteriia; o__Sphingobacteriales; f__; g__; s__
  k__Bacteria; p__Bacteroidetes; c__Cytophagia; o__Cytophagales; f__Flammeovirgaceae; g__; s__
  k__Bacteria; p__Bacteroidetes; c__Cytophagia; o__Cytophagales; f__Cyclobacteriaceae; g__Indibacter; s__alkaliphilus
  k__Bacteria; p__Bacteroidetes; c__Cytophagia; o__Cytophagales; f__Cyclobacteriaceae; g__; s__
  k__Bacteria; p__Bacteroidetes; c__Cytophagia; o__Cytophagales; f__Flammeovirgaceae; g__; s__
```

Because we have taxonomic annotations for all of the Greengenes sequences (though as you probably have noticed by now, they differ in their specificity), we can next look at taxonomy associated with each of our queries in Greengenes. How do your annotations compare to those from Greengenes, which we'll print out in the next cell?

```python
>>> for q in current_queries:
...     q_id = q.metadata['id']
...     print('Known taxonomy for query %s:\n %s' % (q_id, reference_taxonomy[q_id]))
...     print()
Known taxonomy for query 316642:
 k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__; s__

Known taxonomy for query 3900488:
 k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__

Known taxonomy for query 4478817:
 k__Bacteria; p__Bacteroidetes; c__[Saprospirae]; o__[Saprospirales]; f__Saprospiraceae; g__; s__

Known taxonomy for query 254871:
 k__Bacteria; p__Bacteroidetes; c__Cytophagia; o__Cytophagales; f__Flammeovirgaceae; g__; s__
```

## Reducing the runtime for database searches <link src='0f9232'/>

In the examples above, it's taking on the order of 5-15 seconds to search a single sequence against our subset of Greengenes. This makes sense when you think about the computations that are being performed. For every sequence in our reference database (5000, if you haven't modified the database subsampling step) it is computing the $F$ and $T$ matrices described in [the Pairwise Alignment chapter](alias://a76822), and then tracing back the matrix to compute the aligned sequences. Given all of that, the fact that computation only takes 5-15 seconds is pretty incredible. However, that doesn't change the fact that this very likely won't scale to modern-sized sequence data sets. In other words, we'd have to wait way too long for results. Performing all pairwise alignments is prohibitively expensive for database searching.

As we discussed in the previous chapter, the run time of pairwise alignment scales quadratically with sequence length. Database searching, at least in the example we're exploring in this chapter, is a bit of a different problem however. Our sequence lengths aren't changing, but rather it takes a long time because we're performing a computationally expensive step, pairwise alignment, many times. Our database is fixed in that the number of sequences in it doesn't change, the sequences themselves don't change, and they're roughly the same length. Our query sequences are also exactly the same length in this example (remember that we set that above, when we sliced a single region from reference database sequences to create our query sequences). Let's explore how the runtime of this database search scales under these constriants.

```python
>>> import pandas as pd
...
>>> def tabulate_local_alignment_search_runtime(queries, reference_db, n_query_sequences, n_reference_sequences,
...                                             search_function):
...     data = []
...     # we'll iterate over the pairs of number of query sequences
...     # and number of reference sequences, and compute the runtime
...     # of the database search three times for each pair (so we
...     # have some idea of the variance in the runtimes). this is
...     # achieved here with a nested for loop (i.e., a for loop
...     # within a for loop).
...     for nq, nr in zip(n_query_sequences, n_reference_sequences):
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
... n_reference_sequences = [100, 100, 100]
>>> # since our database is smaller, we can work with some slightly
... # larger numbers of sequences.
... n_query_sequences = [1, 5, 10, 15]
...
>>> local_alignment_search_runtimes = tabulate_local_alignment_search_runtime(queries, reference_db,
...                                                                           n_query_sequences, n_reference_sequences,
...                                                                           local_alignment_search)
>>> local_alignment_search_runtimes
   Number of query seqs  Number of reference seqs  Median query seq length  \
0                     1                       100                      200   
1                     1                       100                      200   
2                     1                       100                      200   
3                     5                       100                      200   
4                     5                       100                      200   
5                     5                       100                      200   
6                    10                       100                      200   
7                    10                       100                      200   
8                    10                       100                      200   

   Median reference seq length  Runtime (s)  
0                       1424.5     0.182923  
1                       1450.0     0.170758  
2                       1414.0     0.172627  
3                       1437.0     0.853667  
4                       1446.0     0.867526  
5                       1437.5     0.887163  
6                       1408.0     1.696071  
7                       1428.0     1.729055  
8                       1422.0     1.708480
```

This table shows that we've tried a few variations on number of query sequences and kept the number of reference sequences constant. The is no variance in the query sequence length, and there is a relatively small amount of variance in reference sequence length (they're all of the same order of magnitude). There is also relatively little variance in runtime for fixed numbers of query and reference sequences.

The interesting aspect of this table is that there is an increase in runtime with an increasing number of query sequences. We expect that, but what we care about is how runtime is increasing as a function of number of query sequences. Let's plot runtime versus the number of query sequences to help us understand the relationship.

```python
>>> import seaborn as sns
>>> ax = sns.regplot(x="Number of query seqs", y="Runtime (s)", data=local_alignment_search_runtimes)
>>> ax.set_xlim(0)
>>> ax.set_ylim(0)
>>> ax
```

What we see here is pretty clearly a linear relationship: $runtime \approx constant \times number\ of\ query\ sequences$. This is because as we increase the number of query sequences, we're increasing the number of pairwise alignments that we need to perform. If we have 5 queries and 10 reference sequences, we compute $5 \times 10 = 50$ pairwise alignments. If we have 10 queries and 100 reference sequences, we compute $10 \times 100 = 1000$ pairwise alignments. There are a few practical ways to reduce the runtime of a process like this. 

The first seems obvious, and even silly at first: perform fewer alignments. This could be achieved in a few ways. You could reduce the number of query sequences, though this might be something a researcher is resistant to: they have some collection of unknown sequences, and they want to know what they all are. You could alternatively reduce the number of reference sequences, but you might run into the same issues there: we wouldn't want to exclude reference sequences that might provide us with useful information about our query sequences. Finally, we might be able to figure out some ways to perform fewer alignments by not searching all of the query sequences against all of the reference sequences. 

Another approach to reducing the runtime of this process would be to create a faster implemention of the algorithm (though at some point that won't be possible anymore), use a faster computer, or run the process in parallel on multiple processors. All of these would be ways to reduce the runtime of the search by some factor $f$, where $new\ runtime \approx \frac{runtime}{f}$. 

In practice, for a production scale sequence database search application like BLAST, we'd combine these approaches. In the next section we'll explore wa

## Heurisitic algorithms <link src="mUArdw"/>

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
