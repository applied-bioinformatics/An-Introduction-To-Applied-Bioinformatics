[![Build Status](https://travis-ci.org/gregcaporaso/An-Introduction-To-Applied-Bioinformatics.png?branch=master)](https://travis-ci.org/gregcaporaso/An-Introduction-To-Applied-Bioinformatics) [![DOI](https://zenodo.org/badge/11339/gregcaporaso/An-Introduction-To-Applied-Bioinformatics.svg)](http://dx.doi.org/10.5281/zenodo.16419)

An Introduction To Applied Bioinformatics
=========================================

<div style="float: right; margin-left: 30px;"><img title="Logo by @gregcaporaso." style="float: right;margin-left: 30px;" src="https://raw.github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/images/logo.png" align=right height=300/></div>

Bioinformatics, as I see it, is the application of the tools of computer science (things like programming languages, algorithms, and databases) to address biological problems (for example, inferring the evolutionary relationship between a group of organisms based on fragments of their genomes, or understanding if or how the community of microorganisms that live in my gut changes if I modify my diet). Bioinformatics is a rapidly growing field, largely in response to the vast increase in the quantity of data that biologists now grapple with. Students from varied disciplines (e.g., biology, computer science, statistics, and biochemistry) and stages of their educational careers (undergraduate, graduate, or postdoctoral) are becoming interested in bioinformatics.

I teach bioinformatics at the undergraduate and graduate levels at Northern Arizona University. This repository contains some of the materials that I've developed in these courses, and represents an initial attempt to organize these materials in a standalone way. If you'd like to read a little more about the project, see my [blog post](http://microbe.net/2014/05/01/teaching-bioinformatics-using-ipython-notebooks/) on [microbe.net](http://microbe.net).

Disclaimer
----------

**This project is in very early development stage.** It's not ready for prime-time by any means, but I fall firmly into the "publish early, publish often" mindset, hence its public availability. I am very interested in feedback. The best way to get feedback to me is through [the IAB issue tracker](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/issues), and the best way to get me contributions is through [pull requests](https://help.github.com/articles/using-pull-requests).

The code in the iab module is **not sufficiently tested, documented, or optimized for production use**. As code reaches those quality standards it will be ported to [scikit-bio](http://www.scikit-bio.org). I do not recommend using the code in the iab module outside of these notebooks. In other words, don't `import iab` outside of the notebooks - if you want access to the functionality in your own code, you should `import skbio`.

Currently, the **best example of where I'm hoping to go with these materials** is the [multiple sequence alignment](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/algorithms/4-multiple-sequence-alignment.ipynb) chapter.

Outline
-------

To browse the book, [start here](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/Index.ipynb).

0. Getting started
1. Fundamentals
  1. Pairwise alignment (contains an exercise)
  2. Database searching and determining the statistical significance of an alignment
  3. Phylogeny reconstruction: distances, distances matrices and hierarchical clustering with UPGMA
  4. Multiple sequence alignment (contains an exercise)
  5. Read mapping and clustering
2. Applications
  1. Studying biological diversity

How to read the book
--------------------

There are two ways to read *An Introduction To Applied Bioinformatics*:

* The *recommended* way to read the book is to download and run the IPython notebooks interactively. You can do this by cloning the GitHub repository, installing the package and its dependencies, and running the notebooks interactively. Instructions for doing this are provided below in the **Installation** section.

* The *easiest* way to read the book is to view the static notebooks online using [nbviewer](http://nbviewer.ipython.org/). You should [start here](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/Index.ipynb).

If you're new to using IPython or the IPython Notebook, you can find more information at the [IPython website](http://www.ipython.org/), [IPython Notebook website](http://ipython.org/notebook), and the [IPython Notebook example gallery](https://github.com/ipython/ipython/wiki/A-gallery-of-interesting-IPython-Notebooks).

Installation
------------

If you're going to read the book interactively (recommended), you'll need to clone this repository, install some dependencies, and launch the IPython Notebook. For example, the following commands should work for Linux and Mac OS X users:

    pip install numpy
    pip install https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/archive/master.zip

Finally, launch the IPython Notebook to get started (be sure that you're in the ``An-Introduction-To-Applied-Bioinformatics`` directory when you run this command):

    ipython notebook --pylab inline Index.ipynb

That's it!

More information
----------------

These materials are primarily being developed by [Greg Caporaso](http://caporasolab.us/people/greg-caporaso/) (GitHub: [@gregcaporaso](https://github.com/gregcaporaso)) in the [Caporaso Lab](http://www.caporasolab.us) at [Northern Arizona University](http://www.nau.edu). You can find information on the courses I teach on [my teaching website](http://www.caporasolab.us/teaching) and information on my research and lab on [my lab website](http://www.caporasolab.us).

See the repository's [contributors page](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/graphs/contributors) for information on who has contributed to the project.

Acknowledgements
----------------
Development of *An Introduction to Applied Bioinformatics* was supported in part by [Arizona's Technology and Research Initiative Fund](http://nau.edu/Research/Funding/Technology-Research-Initiative-Fund/). The style of the project was inspired by [Bayesian Methods for Hackers](http://camdavidsonpilon.github.io/Probabilistic-Programming-and-Bayesian-Methods-for-Hackers/).

I want to thank the [IPython Developers](https://github.com/ipython/ipython/graphs/contributors) for all of their work on the [IPython Notebook](http://www.ipython.org/notebook), as well as the [QIIME developers](https://github.com/biocore/qiime/graphs/contributors) and [scikit-bio developers](https://github.com/biocore/scikit-bio/graphs/contributors) for the countless discussions over the years that helped me develop my understanding of the material presented here. This project wouldn't be possible without all of you, and I look forward to many more years of productive, fun and exciting work together!

License
-------

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" href="http://purl.org/dc/dcmitype/InteractiveResource" property="dct:title" rel="dct:type">An Introduction to Applied Bioinformatics</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="http://www.caporasolab.us" property="cc:attributionName" rel="cc:attributionURL">The Caporaso Laboratory</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics" rel="dct:source">https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics</a>.
