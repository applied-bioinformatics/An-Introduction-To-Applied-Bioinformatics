# Contributing to *An Introduction to Applied Bioinformatics*

*An Introduction to Applied Bioinformatics* (IAB), is an open source project and we welcome community contributions. Contributions should generally be made in the form of GitHub pull requests. We've simplified the process of submitting these for changes to the IAB content, as described below.

You will need a (free) GitHub account to submit a pull request.

Before considering a contribution to IAB, please read the project [Code of Conduct](https://github.com/applied-bioinformatics/An-Introduction-to-Applied-Bioinformatics/blob/master/CODE-OF-CONDUCT.md).

## How IAB is developed

IAB is written in [markdown](http://commonmark.org/). [build-iab](https://github.com/caporaso-lab/build-iab) is then used to convert that markdown to html (for static online viewing) and ipynb (for interactive use). If you're submitting changes to content, you'll be submitting changes to markdown files. These are much more manageable that submitting changes to IPython notebooks, as it's much easier to diff the content.

## Small contributions

Contributions such as typo fixes, wording changes, and small code updates are relatively easy to submit. Every unit, chapter, section, and subsection heading in IAB has an *Edit on GitHub* link below it. You should click the link corresponding to the section that you want to edit, which will take you to the GitHub online editor. You can make your changes and submit a [pull request](https://help.github.com/articles/using-pull-requests/) from that page. This will initiate the IAB tests, which will ensure that your change hasn't broken anything. After the tests pass, we will review your changes and either merge them, request modifications before merging, or let you know why we won't integrate your proposed changes.

Watch a five-minute YouTube video illustrating this process [here](https://www.youtube.com/watch?v=s9-aZrX5CY8).

## Large contributions

If you're interested in making contributions involving code refactoring, new chapters or sections, restructuring of content, etc, you should first comment on existing issues (or create new issues) indicating what you want to work on, and wait for us to discuss the changes with you before you get started. This will ensure that the changes you want to make are in line with the direction of the project, compatible with future plans, and that there is no one else already working on related contributions. This will avoid wasted time if your change involves something that we wouldn't ultimately be merged. If in doubt, bring it up on the [issue tracker](https://github.com/caporaso-lab/An-Introduction-To-Applied-Bioinformatics/issues/) first.

## Who's contributing?

See the list of [IAB contributors](https://github.com/caporaso-lab/An-Introduction-To-Applied-Bioinformatics/graphs/contributors) to find out who is involved with the project. If you submit a pull request that is merged, your GitHub account will be automatically listed on that page.

## Technical points

### Building IAB locally

If you're interested in building the IAB html and/or IPython Notebooks locally, you'll need to install IAB and build-iab. You can do this as follows:

```
pip install numpy
pip install https://github.com/caporaso-lab/build-iab/archive/master.zip
wget https://github.com/caporaso-lab/An-Introduction-To-Applied-Bioinformatics/archive/master.zip
unzip master.zip
cd An-Introduction-To-Applied-Bioinformatics-master/
pip install .
```

Then, to build the IPython Notebooks, you can run:

```
biab notebook -i book -o ipynb
```

or to build the HTML version, you can run:

```
biab html -i book -o html
```

### Linking to other sections of the text

All section headings must have ids associated with them. Should be generated as follows:

```bash
$ biab idgen
<link src="9mM4Bb"/>
```

When you define a section heading, you'd end it with the tag returned from the above command. For example:

```markdown
## Some section <link src="9mM4Bb"/>
```

If you then wanted to link to that section from somewhere else in the text, you could do that with a markdown link as follows:

```markdown
This concept is discussed in further detail [above](alias://9mM4Bb).
```

You should always link using these ids, and never statically link to other sections of the text with URLs (because a section name might change, but its id won't).

## License and license changes

The IAB license is available [here](https://github.com/caporaso-lab/An-Introduction-To-Applied-Bioinformatics/blob/master/LICENSE). This license may change over time, but the online version of IAB will always be the most current version, and will be available free of charge.

By contributing to IAB, you are agreeing that Greg Caporaso has sole discretion over the license and any future changes to the license. If a paid (e.g., printed) copy of IAB is ever created, contributors are not entitled to payments or royalties of any kind. Your contribution of content represents your agreement with these terms.
