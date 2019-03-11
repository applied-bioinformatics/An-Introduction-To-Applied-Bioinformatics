# Reading An Introduction to Applied Bioinformatics <link src='74cadb'/>

**Bioinformatics, as I see it, is the application of the tools of computer science (things like programming languages, algorithms, and databases) to address biological problems (for example, inferring the evolutionary relationship between a group of organisms based on fragments of their genomes, or understanding if or how the community of microorganisms that live in my gut changes if I modify my diet).** Bioinformatics is a rapidly growing field, largely in response to the vast increase in the quantity of data that biologists now grapple with. Students from varied disciplines (e.g., biology, computer science, statistics, and biochemistry) and stages of their educational careers (undergraduate, graduate, or postdoctoral) are becoming interested in bioinformatics.

*An **I**ntroduction to **A**pplied **B**ioinformatics*, or **IAB**, is a free, open access bioinformatics text available at http://readIAB.org. **It introduces readers to the core concepts of bioinformatics in the context of their implementation and application to real-world problems and data.** IAB makes extensive use of the [scikit-bio](http://www.scikit-bio.org) Python package, which provides production-ready implementations of core bioinformatics algorithms and data structures. As readers are learning a concept, for example, pairwise sequence alignment, they are presented with its scikit-bio implementation directly in the text. scikit-bio code is well annotated (adhering to the [pep8](https://www.python.org/dev/peps/pep-0008/) and [numpydoc](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt) conventions), so readers can use it to assist with their understanding of the concept. Readers of IAB also therefore learn the concepts in the context of tools they can use to develop their own bioinformatics software and pipelines, enabling them to rapidly get started on their own projects. While some theory is discussed, the focus of IAB is on what readers need to know to be effective, practicing bioinformaticians.

IAB is **completely open access**, with all software being BSD-licensed, and all text being licensed under Creative Commons Attribution Only (i.e., CC BY-NC-SA 4.0). All development and publication is coordinated under [public revision control](https://github.com/applied-bioinformatics/An-Introduction-To-Applied-Bioinformatics).

[My](alias://fedd13) goal for IAB is for it to make bioinformatics as accessible as possible to students from varied backgrounds, and to get more and diverse people into this hugely exciting field. I'm very interested in hearing from readers and instructors who are using IAB, so get in touch if you have corrections, suggestions for how to improve the content, or any other thoughts or comments on the text. In the spirit of openness, I'd prefer to be contacted via the [IAB issue tracker](https://github.com/applied-bioinformatics/An-Introduction-To-Applied-Bioinformatics/issues/). I'll respond to direct e-mail as well, but I'm often backlogged on e-mail (just ask my students), so e-mail responses are likely to be slower.

I hope you find IAB useful, and that you enjoy reading it!

## Who should read IAB? <link src='976b11'/>

IAB is written for scientists, software developers, and students interested in understanding and applying bioinformatics methods, and ultimately in developing their own bioinformatics analysis pipelines or software.

IAB was initially developed for an undergraduate course cross-listed in computer science and biology with no pre-requisites. It therefore assumes little background in biology or computer science, however some basic background is very helpful. For example, an understanding of the roles of and relationship between DNA and protein in a cell, and the ability to read and follow well-annotated python code, are both helpful (but not necessary) to get started.

In the *Getting started with [Biology](alias://cf88ac) and [Computer Science](alias://6ad7e1)* sections below I provide some suggestions for other texts that will help you to get started.

## How to read IAB <link src='de1bc2'/>

IAB can be read interactively as a series of [Jupyter Notebooks](https://jupyter.org) or read statically. Due to popular demand, a print version may ultimately be available for a fee, but the full and most recent version of IAB will always be available for free on the [project website](http://readIAB.org). The recommended way to read IAB is interactively as this allows readers to execute code directly in the text. For example, when learning pairwise alignment, users can align sequences provided in IAB (or their own sequences) and modify parameters (or even the algorithm itself) to see how changes affect the resulting alignments.

IAB is constantly being updated. As I teach with it, I will often update text or add new chapters in an effort to keep up with advances in the field. The project website contains the most up-to-date recommendations on how to read IAB or teach with IAB, including strategies for dealing with changing content. (For example, if you're teaching with IAB, you can fork the IAB repository and only pull updates into your fork when you're ready for them. If _forking repositories_ and _pulling updates_ are terms that don't mean anything to you right now, you can safely ignore this!)

IAB is split into four different sections: *Getting started*, *Fundamentals*, *Applications*, and *Wrapping up*. You should start reading IAB by working through the *Getting started* and *Fundamentals* chapters in order. You should then read the *Applications* chapters and *Wrapping up* in any order, based on your own interest.

## Using Jupyter Notebooks to read IAB interactively <link src='fdf5dd'/>

IAB can be read interactively as a series of Jupyter Notebooks. The main source for information about Jupyter Notebooks is the [Jupyter website](https://jupyter.org). You can find information there on how to use Jupyter Notebooks as well as setting up and running a Jupyter Notebook server (for example, if you'd like to make one available to your students).

Most of the code that is used in IAB comes from [scikit-bio](http://scikit-bio.org) package, or other Python scientific computing tools. You can access these in the same way that you would in a Python script. For example:

```python
>>> import skbio
>>> from IPython.core import page
>>> page.page = print
```

We can then access functions, variables, and classes from these modules.

```python
>>> print(skbio.title)
>>> print(skbio.art)
```

We'll inspect a lot of source code in IAB as we explore bioinformatics algorithms. If you're ever interested in seeing the source code for some functionality that we're using, you can do that using Jupyter's ``psource`` magic.

```python
>>> from skbio.alignment import TabularMSA
>>> %psource TabularMSA.conservation
```

The documentation for scikit-bio is also very extensive. You can view the documentation for the `TabularMSA` object, for example, [here](http://scikit-bio.org/docs/latest/generated/skbio.alignment.TabularMSA.html). These documents will be invaluable for learning how to use the objects.

## Reading list <link src="22AJOO"/>

### Getting started with Biology <link src='cf88ac'/>

If you're new to biology, these are some books and resources that will help you get started.

* [The Processes of Life](http://amzn.to/1P0dc2E) by Lawrence Hunter. *An introduction to biology for computer scientists.*

* The [NIH Bookshelf](http://www.ncbi.nlm.nih.gov/books/) A lot of free biology texts, some obviously better than others.

### Getting started with Computer Science and programming <link src='6ad7e1'/>

If you're new to Computer Science and programming, these are some books and resources that will help you get started.

* [Software Carpentry](http://www.software-carpentry.org) *Online resources for learning scientific computing skills, and regular in-person workshops all over the world. Taking a Software Carpentry workshop **will** pay off for biology students interested in a career in research.*

* [Practical Computing for Biologists](http://amzn.to/1Ukx5S6) by Steven Haddock and Casey Dunn. *A great introduction to many computational skills that are required of modern biologists. I *highly* recommend this book to all Biology undergraduate and graduate students.*

* [Practical Programming: A Introduction to Computer Science Using Python](http://amzn.to/1P0dmqM) by Jennifer Campbell, Paul Gries, Jason Montojo, Greg Wilson. *An introduction to the python programming language and basic computer science. This is a great first programming book for people interested in bioinformatics or scientific computing in general.*

* [The Pragmatic Programmer](http://amzn.to/1P0dl6i) by Andrew Hunt. *A more advanced book on becoming a better programmer. This book is excellent, and I highly recommend it for anyone developing bioinformatics software. You should know how to program and have done some software development before jumping into this.*

### Philosophy of biology and popular science books <link src="bwMyRz"/>

These are some books that I've enjoyed, that have also helped me think about biological systems. These are generally written for a more popular audience, so should be accessible to any readers of *An Introduction to Applied Bioinformatics*.

* [The Selfish Gene](http://amzn.to/1UkyQ1R) by Richard Dawkins.

* [Ever Since Darwin](http://amzn.to/1Ukzdt7) by Stephen Jay Gould. *This is the first book in a series of collections of short essays.*

* [The Demon Haunted World](http://amzn.to/1UkyIzi) by Carl Sagan.

* [Sex and Death](http://amzn.to/1UkySXg) by Kim Sterelny.

* [Gödel, Escher, Bach](http://amzn.to/1UkzxYL) by Douglas Hofstadter.

## Need help? <link src='da2930'/>

If you're having issues getting *An Introduction to Applied Bioinformatics* running on your computer, or you have corrections or suggestions on the content, you should get in touch through the [IAB issue tracker](https://github.com/applied-bioinformatics/An-Introduction-To-Applied-Bioinformatics/issues). This will generally be much faster than e-mailing the author directly, as there are multiple people who monitor the issue tracker. It also helps us manage our technical support load if we can consolidate all requests and responses in one place.

## Contributing and Code of Conduct <link src='b1c06a'/>

If you're interested in contributing content or features to IAB, you should start by reviewing the project's [Code of Conduct](https://github.com/applied-bioinformatics/An-Introduction-To-Applied-Bioinformatics/blob/master/CODE-OF-CONDUCT.md) and [Contributing Guide](https://github.com/applied-bioinformatics/An-Introduction-To-Applied-Bioinformatics/blob/master/CONTRIBUTING.md).

## Acknowledgements <link src="CLeix6"/>

*An Introduction to Applied Bioinformatics* was funded in part by the [Alfred P. Sloan Foundation](www.sloan.org). Initial prototyping was funded by [Arizona's Technology and Research Initiative Fund](http://nau.edu/Research/Funding/Technology-Research-Initiative-Fund/). The style of the project was inspired by [Bayesian Methods for Hackers](http://camdavidsonpilon.github.io/Probabilistic-Programming-and-Bayesian-Methods-for-Hackers/).

See the repository's [contributors page](https://github.com/applied-bioinformatics/An-Introduction-To-Applied-Bioinformatics/graphs/contributors) for information on who has contributed to the project.
