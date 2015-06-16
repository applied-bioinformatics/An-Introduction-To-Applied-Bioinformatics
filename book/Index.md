
# An Introduction To Applied Bioinformatics

## Table of Contents

<div style="float: right; margin-left: 30px; width: 200px"><img title="Logo by @gregcaporaso." style="float: right;margin-left: 30px;" src="https://raw.github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/images/logo.png" align=right height=250/></div>

1. Getting started
  1. [Reading *An Introduction To Applied Bioinformatics*](getting-started/reading-iab.ipynb)
2. Fundamentals
  1. [Pairwise alignment](fundamentals/pairwise-alignment.ipynb) ([exercise](fundamentals/pairwise-alignment-exercises.ipynb))
  2. [Database searching and determining the statistical significance of an alignment](fundamentals/database-searching.ipynb)
  3. [Phylogeny reconstruction: distances, distances matrices and hierarchical clustering with UPGMA](fundamentals/phylogeny-reconstruction.ipynb)
  4. [Multiple sequence alignment](fundamentals/multiple-sequence-alignment.ipynb) ([exercise](fundamentals/msa-assignment.ipynb))
  5. [Read mapping and clustering](fundamentals/sequence-mapping-and-clustering.ipynb?create=1)
3. Applications
  1. [Studying biological diversity](applications/biological-diversity.ipynb)
4. Wrapping up

## Disclaimer

**This project is in very early development stage.** It's not ready for prime-time by any means, but I fall firmly into the "publish early, publish often" mindset, hence its public availability. I am very interested in feedback. The best way to get feedback to me is through [the IAB issue tracker](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/issues), and the best way to contribute is through [pull requests](https://help.github.com/articles/using-pull-requests).

The code in the iab module is **not sufficiently tested, documented, or optimized for production use**. As code reaches those quality standards it will be ported to [scikit-bio](http://www.scikit-bio.org). I do not recommend using the code in the iab module outside of these notebooks. In other words, don't `import iab` outside of the notebooks - if you want access to the functionality in your own code, you should `import skbio`.

Currently, the **best example of where I'm hoping to go with these materials** is the [multiple sequence alignment](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/algorithms/4-multiple-sequence-alignment.ipynb) chapter.