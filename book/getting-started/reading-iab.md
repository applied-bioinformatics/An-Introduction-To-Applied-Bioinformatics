# Reading An Introduction to Applied Bioinformatics <link src='74cadb'/>

**Bioinformatics, as I see it, is the application of the tools of computer science (things like programming languages, algorithms, and databases) to address biological problems (for example, inferring the evolutionary relationship between a group of organisms based on fragments of their genomes, or understanding if or how the community of microorganisms that live in my gut changes if I modify my diet).** Bioinformatics is a rapidly growing field, largely in response to the vast increase in the quantity of data that biologists now grapple with. Students from varied disciplines (e.g., biology, computer science, statistics, and biochemistry) and stages of their educational careers (undergraduate, graduate, or postdoctoral) are becoming interested in bioinformatics.

*An **I**ntroduction to **A**pplied **B**ioinformatics*, or **IAB**, is an open source, interactive bioinformatics text. **It introduces readers to the core concepts of bioinformatics in the context of their implementation and application to real-world problems and data.** IAB is closely tied to the [scikit-bio](http://www.scikit-bio.org) python package, which provides production-ready implementations of core bioinformatics algorithms and data structures. Readers therefore learn the concepts in the context of tools they can use to develop their own bioinformatics software and pipelines, enabling them to rapidly get started on their own projects. While some theory is discussed, the focus of IAB is on what readers need to know to be effective, practicing bioinformaticians.

IAB is interactive, being **based on IPython Notebooks** which can be installed on a reader's computer or viewed statically online. As readers are learning a concept, for example, pairwise sequence alignment, they are presented with its scikit-bio implementation directly in the text. scikit-bio code is well annotated (adhering to the [pep8](https://www.python.org/dev/peps/pep-0008/) and [numpydoc](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt) conventions), so readers can use it to assist with their understanding of the concept. And, because IAB is presented as an IPython Notebook, readers can execute the code directly in the text. For example, when learning pairwise alignment, users can align sequences provided in IAB (or their own sequences) and modify parameters (or even the algorithm itself) to see how changes affect the resulting alignments.

IAB is **completely open access**, with all software being BSD-licensed, and all text being licensed under Creative Commons Attribution Only (i.e., CC BY-NC-SA 4.0). All development and publication is coordinated under [public revision control on GitHub](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics).

IAB is also an **electronic-only resource**. There are currently no plans to commercialize it or to create a print version. This means that, unlike printed bioinformatics texts which are generally out of date before the ink dries, IAB can be updated as the field changes.

**The life cycle of IAB is more like a software package than a book.** There will be development and release versions of IAB, where the release versions are more polished but won't always contain the latest content, and the development versions will contain all of the latest materials, but won't necessarily be copy-edited and polished.

We are in the process of developing a **project status page** that will detail the plans for IAB. This will include the full table of contents, and what stage you can expect chapters to be at at different times. You can track progress of this on [IAB #97](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/issues/97).

[My](alias://fedd13) goal for IAB is for it to make bioinformatics as accessible as possible to students from varied backgrounds, and to get more people into this hugely exciting field. I'm very interested in hearing from readers and instructors who are using IAB, so get in touch if you have corrections, suggestions for how to improve the content, or any other thoughts or comments on the text. In the spirit of openness, I'd prefer to be contacted via the [IAB issue tracker](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/issues/). I'll respond to direct e-mail as well, but I'm always backlogged (just ask my students), so e-mail responses are likely to be slower.

I hope you find IAB useful, and that you enjoy reading it!

## Who should read IAB? <link src='976b11'/>

IAB is written for scientists, software developers, and students interested in understanding and applying bioinformatics methods, and ultimately in developing their own bioinformatics analysis pipelines or software.

IAB was initially developed for an undergraduate course cross-listed in computer science and biology with no pre-requisites. It therefore assumes little background in biology or computer science, however some basic background is very helpful. For example, an understanding of the roles of and relationship between DNA and protein in a cell, and the ability to read and follow well-annotated python code, are both helpful (but not necessary) to get started.

In the *Getting started with [Biology](alias://cf88ac) and [Computer Science](alias://6ad7e1)* sections below I provide some suggestions for other texts that will help you to get started.

## How to read IAB <link src='de1bc2'/>

There are two ways to read *An Introduction To Applied Bioinformatics*:

* The *recommended* way is to install it and work with it interactively. See the [install instructions below](alias://96d24f).
* The *easiest* way is to view the static notebooks online at  [readIAB.org](http://readIAB.org/). You should:
 * [start here to view the latest release version](http://readIAB.org/book/latest/), or
 * [start here to view the development version](http://readIAB.org/book/dev/) (which will have the latest content, but will be less polished and possibly buggy).

IAB is split into four different sections: *Getting started*, *Fundamentals*, *Applications*, and *Wrapping up*. You should start reading IAB by working through the *Getting started* and *Fundamentals* chapters in order. You should then read the *Applications* chapters and *Wrapping up* in any order, based on your own interest.

## Installation <link src='96d24f'/>

If you're going to read *An Introduction to Applied Bioinformatics* interactively (recommended), you'll need to install it. The easiest way to install IAB is using [``conda``](http://conda.pydata.org/docs/intro.html), a free package and environment manager for Python. You'll first need to install Miniconda following the instructions [here](http://conda.pydata.org/docs/install/quick.html#miniconda-quick-install-requirements).

After installing Miniconda, change to the directory where you'd like to store the IAB notebooks (e.g., ``cd $HOME/Documents``). Then, run:

```bash
conda create --name iab -c http://conda.binstar.org/gregcaporaso python=3.4 iab
```

This will create a new ``iab`` environment, install IAB and its dependencies in it, and download the IAB notebooks. Once that completes, you should be able to run:

```bash
source activate iab
cd $HOME/Documents/IAB-notebooks
ipython notebook index.ipynb
```

to get started. If you installed in a directory other than ``$HOME/Documents``, you'll need to adapt the second command to refer to the directory containing ``IAB-notebooks``. Whenever you want to start reading and working with IAB, you'll run these three commands.

## Using the IPython Notebook <link src='fdf5dd'/>

IAB is built using the IPython Notebook, an interactive HTML-based computing environment. The main source for information about the IPython Notebook is the [IPython Notebook website](http://ipython.org/notebook). You can find information there on how to use the IPython Notebook, and also on how to set up and run and IPython Notebook server (for example, if you'd like to make one available to your students).

Most of the code that is used in IAB comes from [scikit-bio](http://scikit-bio.org) package, or other python scientific computing tools. You can access these in the same way that you would in a python script. For example:

```python
>>> import skbio
>>> from __future__ import print_function
>>> from IPython.core import page
>>> page.page = print
```

We can then access functions, variables, and classes from these modules.

```python
>>> print(skbio.title)
>>> print(skbio.art)
```

We'll inspect a lot of source code in IAB as we explore bioinformatics algorithms. If you're ever interested in seeing the source code for some functionality that we're using, you can do that using IPython's ``psource`` magic.

```python
>>> from skbio.alignment import TabularMSA
>>> %psource TabularMSA.conservation
```

The documentation for scikit-bio is also very extensive (though the package itself is still in early development). You can view the documentation for the `TabularMSA` object, for example, [here](http://scikit-bio.org/docs/latest/generated/skbio.alignment.TabularMSA.html). These documents will be invaluable for learning how to use the objects.

## Reading list

### Getting started with Biology <link src='cf88ac'/>

If you're new to biology, these are some books and resources that will help you get started.

* [The Processes of Life](http://amzn.to/1P0dc2E) by Lawrence Hunter. *An introduction to biology for computer scientists.*

* The [NIH Bookshelf](http://www.ncbi.nlm.nih.gov/books/) A lot of free biology texts, some obviously better than others.

* [Molecular Biology of the Cell](http://amzn.to/1UkwYpL) by Bruce Alberts, Alexander Johnson, Julian Lewis, Martin Raff, Keith Roberts, Peter Walter. *One of the best texts on molecular biology. This is fairly advanced (it's generally used in upper division molecular biology courses) so it may not be the best place to start. You'll find it invaluable though if you plan to go on in Bioinformatics. This book is available via the NIH Bookshelf (for example, from Chapter 1: [The Universal Features of Cells on Earth](http://www.ncbi.nlm.nih.gov/books/NBK26864/) and [The Diversity of Genomes and the Tree of Life](http://www.ncbi.nlm.nih.gov/books/NBK26866/).*

* [Brock Biology of Microorganisms](http://amzn.to/1P0derp) by Michael T. Madigan, John M. Martinko, David Stahl, David P. Clark. *One of the best textbooks on microbiology. This is also fairly advanced, but if you're interested in microbial ecology or other aspects of microbiology it will likely be extremely useful.*

### Getting started with Computer Science and programming <link src='6ad7e1'/>

If you're new to Computer Science and programming, these are some books and resources that will help you get started.

* [Software Carpentry](www.software-carpentry.org) *Online resources for learning scientific computing skills, and regular in-person workshops all over the world. Taking a Software Carpentry workshop **will** pay off for biology students interested in a career in research.*

* [Practical Computing for Biologists](http://amzn.to/1Ukx5S6) by Steven Haddock and Casey Dunn. *A great introduction to many computational skills that are required of modern biologists. I *highly* recommend this book to all Biology undergraduate and graduate students.*

* [Practical Programming: A Introduction to Computer Science Using Python](http://amzn.to/1P0dmqM) by Jennifer Campbell, Paul Gries, Jason Montojo, Greg Wilson. *An introduction to the python programming language and basic computer science. This is a great first programming book for people interested in bioinformatics or scientific computing in general.*

* [The Pragmatic Programmer](http://amzn.to/1P0dl6i) by Andrew Hunt. *A more advanced book on becoming a better programmer. This book is excellent, and I highly recommend it for anyone developing bioinformatics software. You should know how to program and have done some software development before jumping into this.*

### Philosophy of biology and popular science books

These are some books that I've enjoyed, and which helped get me thinking about biological systems. These are generally written for a more popular audience, so should be accessible to any readers of *An Introduction to Applied Bioinformatics*.

* [The Selfish Gene](http://amzn.to/1UkyQ1R) by Richard Dawkins.

* [Ever Since Darwin](http://amzn.to/1Ukzdt7) by Stephen Jay Gould. *This is the first book in a series of collections of short essays.*

* [The Demon Haunted World](http://amzn.to/1UkyIzi) by Carl Sagan.

* [Sex and Death](http://amzn.to/1UkySXg) by Kim Sterelny.

* [Gödel, Escher, Bach](http://amzn.to/1UkzxYL) by Douglas Hofstadter.

## Need help? <link src='da2930'/>

If you're having issues getting *An Introduction to Applied Bioinformatics* running on your computer, or you have corrections or suggestions on the content, you should get in touch through the [IAB issue tracker](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/issues). This will generally be much faster than e-mailing the author directly, as there are multiple people who monitor the issue tracker. It also helps us manage our technical support load if we can consolidate all requests and responses in one place.

## Contributing to IAB <link src='b1c06a'/>

If you're interested in contributing content or features to IAB, you should start by reviewing [CONTRIBUTING.md](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/CONTRIBUTING.md) which provides guidelines on how to get involved.

## About the author <link src='fedd13'/>

My name is Greg Caporaso. I'm the primary author of *An Introduction to Applied Bioinformatics*, but there are [other contributors](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/graphs/contributors) and I hope that list will grow.

<div style="float: left; margin-left: 20px; margin-right: 20px; margin-top: 25px; margin-bottom: 15px; width: 200px">
 <img title="The author in Telluride, CO (April, 2015)." style="float: right;" src="https://raw.github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/getting-started/images/greg-in-telluride.png" height=250/>
</div>

I have degrees in Computer Science (B.S., University of Colorado, 2001) and Biochemistry (B.A., University of Colorado, 2004; Ph.D., University of Colorado 2009). Following my formal training, I joined the [Rob Knight Laboratory](http://knightlab.ucsd.edu), then at the University of Colorado, for approximately 2 years as a post-doctoral scholar. In 2011, I joined the faculty at [Northern Arizona University (NAU)](www.nau.edu) where I'm an Assistant Professor in the Biological Sciences department. I [teach](http://www.caporasolab.us/teaching/) one course per year in bioinformatics for graduate and undergraduate students of Biology and Computer Science. I also run a [research lab](http://www.caporasolab.us/) in the [Center for Microbial Genetics and Genomics](http://www.mggen.nau.edu/), which is focused on developing bioinformatics software and studying the human microbiome.

I'm not the world expert on the topics that I present in IAB, but I have a passion for bioinformatics, open source software, writing, and education. When I'm learning a new bioinformatics concept, for example an algorithm like pairwise alignment or a statistical technique like Monte Carlo simulation, implementing it is usually the best way for me to wrap my head around it. This led me to start developing IAB, as I found that my implementations helped my students learn the concepts too. I think that one of my strongest skills is the ability to break complex ideas into accessible components. I do this well for bioinformatics because I remember (and still regularly experience) the challenges of learning it, so can relate to newcomers in the field.

I'm very active in open source bioinformatics software development. I am most widely known for my involvement in the development of the [QIIME software package](http://www.qiime.org), and more recently for leading the development of [scikit-bio](http://www.scikit-bio.org). I am also involved in many other bioinformatics software projects (see my [GitHub page](http://github.com/gregcaporaso)). IAB is one of the projects that I'm currently most excited about, so I truly hope that it's as useful for you as it is fun for me.

For updates on IAB and various other things, you should [follow me on Twitter](https://twitter.com/gregcaporaso).

## Acknowledgements <link src="CLeix6"/>

<div style="float: right; margin-left: 30px;"><a href="http://www.sloan.org"><img title="Logo by @gregcaporaso." style="float: right;margin-left: 30px;" src="http://www.sloan.org/fileadmin/media/images/logos/SloanLogo-primary-black.png" align=right></a></div>

*An Introduction to Applied Bioinformatics* is funded in part by the [Alfred P. Sloan Foundation](www.sloan.org). Initial prototyping was funded by [Arizona's Technology and Research Initiative Fund](http://nau.edu/Research/Funding/Technology-Research-Initiative-Fund/). The style of the project was inspired by [Bayesian Methods for Hackers](http://camdavidsonpilon.github.io/Probabilistic-Programming-and-Bayesian-Methods-for-Hackers/).

See the repository's [contributors page](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/graphs/contributors) for information on who has contributed to the project.
