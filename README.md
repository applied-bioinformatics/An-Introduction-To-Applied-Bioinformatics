[![Build Status](https://travis-ci.org/gregcaporaso/An-Introduction-To-Applied-Bioinformatics.png?branch=master)](https://travis-ci.org/gregcaporaso/An-Introduction-To-Applied-Bioinformatics)

An Introduction To Applied Bioinformatics
=========================================

<div style="float: right; margin-left: 30px;"><img title="Logo by @gregcaporaso." style="float: right;margin-left: 30px;" src="https://raw.github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/images/logo.png" align=right height=300/></div>

Bioinformatics, as I see it, is the application of the tools of computer science (things like programming languages, algorithms, and databases) to address biological problems (for example, inferring the evolutionary relationship between a group of organisms based on fragments of their genomes, or understanding if or how the community of microorganisms that live in my gut changes if I modify my diet). Bioinformatics is a rapidly growing field, largely in response to the vast increase in the quantity of data that biologists now grapple with. Students from varied disciplines (e.g., biology, computer science, statistics, and biochemistry) and stages of their educational careers (undergraduate, graduate, or postdoctoral) are becoming interested in bioinformatics.

I teach bioinformatics at the undergraduate and graduate levels at Northern Arizona University. This repository contains some of the materials that I've developed in these courses, and represents an initial attempt to organize these materials in a standalone way. In some cases, I'm just linking out to other materials for now. 

Disclaimer
----------

This project is in very early development stage. It's not ready for prime-time by any means, but I fall firmly into the "publish early, publish often" mindset, hence its public availability. I am very interested in feedback in the form of email (gregcaporaso@gmail.com) or [pull requests](https://help.github.com/articles/using-pull-requests).

Currently, the **best example of where I'm hoping to go with these materials** is the [multiple sequence alignment](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/algorithms/4-multiple-sequence-alignment.ipynb) chapter.

Outline
-------

0. [Getting started](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/getting-started/0-overview.ipynb)
1. [Fundamentals](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/algorithms/0-overview.ipynb)
  1. [Pairwise alignment](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/algorithms/1-pairwise-alignment.ipynb) ([exercise](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/algorithms/pairwise-alignment-exercises.ipynb))
  2. [Database searching and determining the statistical significance of an alignment](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/algorithms/2-database-searching.ipynb)
  3. [Phylogeny reconstruction: distances, distances matrices and hierarchical clustering with UPGMA](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/algorithms/3-phylogeny-reconstruction.ipynb)
  4. [Multiple sequence alignment](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/algorithms/4-multiple-sequence-alignment.ipynb) ([exercise](http://nbviewer.ipython.org/github/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/algorithms/msa-assignment.ipynb))
  5. [Read mapping and clustering](http://nbviewer.ipython.org/urls/raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/algorithms/5-sequence-mapping-and-clustering.ipynb?create=1)
2. Applications
  1. [Studying biological diversity](http://nbviewer.ipython.org/urls/raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/applications/1-biological-diversity.ipynb)
  2. Interrogating genomes
3. Wrapping up

How to read the book
--------------------

There are two ways to read *An Introduction To Applied Bioinformatics*:

* The *recommended* way to read the book is to download and run the IPython notebooks interactively. You can do this by cloning the GitHub repository, installing the package and its dependencies, and running the notebooks interactively. Instructions for doing this are provided below in the **Installation** section.

* The *easiest* way to read the book is to view the static notebooks online using [nbviewer](http://nbviewer.ipython.org/). You can access the notebooks on nbviewer from the links in the **Outline** section above.

If you're new to using IPython or the IPython Notebook, you can find more information at the [IPython website](http://www.ipython.org/), [IPython Notebook website](http://ipython.org/notebook), and the [IPython Notebook example gallery](https://github.com/ipython/ipython/wiki/A-gallery-of-interesting-IPython-Notebooks).

Installation
------------

If you're going to read the book interactively (recommended), you'll need to clone this repository, install some dependencies, and launch the IPython Notebook. For example, the following commands should work for Linux and Mac OS X users:

    git clone https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics.git
    cd An-Introduction-To-Applied-Bioinformatics
    pip install .

If you're using a newer version of pip and scikit-bio fails to install, you may need to modify the pip command to be:

    pip install --allow-all-external --allow-unverified scikit-bio --process-dependency-links .

Finally, launch the IPython Notebook and browse to one of the notebook directories to read a chapter:

    ipython notebook

That's it!

If you'd like to install the book's dependencies manually (or some other way than using pip), here's what you'll need:

- [Python](http://www.python.org/) 2.7
- [numpy](http://www.numpy.org/) >= 1.7
- [scipy](http://www.scipy.org/) >= 0.13.0
- [matplotlib](http://www.matplotlib.org/) >= 1.1.0
- [IPython](http://www.ipython.org/)
- [tornado](http://www.tornadoweb.org/en/stable/)
- [pyzmq](http://zeromq.github.io/pyzmq/)
- [jinja2](http://jinja.pocoo.org/)
- [scikit-bio](http://scikit-bio.org/) (latest GitHub version)
- [biom-format](http://www.biom-format.org)
- [pyqi](http://biocore.github.io/pyqi/doc/index.html)

More information
----------------

These materials are primarily being developed by [Greg Caporaso](http://caporasolab.us/people/greg-caporaso/) (GitHub: [@gregcaporaso](https://github.com/gregcaporaso)) at [Northern Arizona University](http://www.nau.edu). You can find information on the courses I teach on [my teaching website](http://www.caporasolab.us/teaching) and information on my research and lab on [my lab website](http://www.caporasolab.us).

See the repository's [contributors page](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/graphs/contributors) for information on who has contributed to the project.

License
-------

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" href="http://purl.org/dc/dcmitype/InteractiveResource" property="dct:title" rel="dct:type">An Introduction to Applied Bioinformatics</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="http://www.caporasolab.us" property="cc:attributionName" rel="cc:attributionURL">The Caporaso Laboratory</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics" rel="dct:source">https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics</a>.
