# Publishing a release of *An Introduction to Applied Bioinformatics*

These notes describe how to build and publish a release of IAB, and how to upload the content to [readIAB.org](http://www.readIAB.org).

## Part 1: Prepare the release.

1. Ensure that the Travis build is passing on the master branch.

1. Update version strings. As of this writing, these need to be updated in ``setup.py``, ``iab/__init__.py``, ``book/getting-started/reading-iab.md``, and ``conda-recipe/meta.yaml``, but it's a good idea to do a global find on the code base for the last release version number to ensure that you're not missing anything. (In the future, the version number will only be specified in one place to simplify this process.)

1. Update [``CHANGELOG.md``](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/CHANGELOG.md) to include descriptions of the changes since the last release.

1. Submit a pull request with these changes and let Travis run.

## Part 2: Build and test the html and IPython Notebook content.

1. Install [build-iab](https://github.com/gregcaporaso/build-iab) if you don't already have it installed.

1. From the top-level ``An-Introduction-to-Applied-Bioinformatics`` directory, build the html and IPython Notebooks with the following commands. Poke through these files to ensure that they look right, and that the *Edit* links work correctly.

   ```bash
   biab html -i book -o html
   biab notebook -i book -o ipynb
   ```

## Part 3: Tag the release.

From the [An Introduction to Applied Bioinformatics GitHub page](https://github.com/gregcaporaso/An-Introduction-to-Applied-Bioinformatics), click on the ``release`` link and draft a new release. Use the version number for the tag name and create the tag against master. Fill in a release title that is consistent with previous release titles and add a summary of the release (linking to ``CHANGELOG.md`` is a good idea). This release summary will be the primary information that we point users to when we announce the release.

## Part 4: Build the conda packages

*These notes are in draft stage right now as they were developed between releases. Will test and finalize during the next release.*

1. Configure your conda build environment by installing ``Miniconda`` and then running ``conda install conda-build``.

1. Download the release that was just created on GitHub, and unzip that. Change into the resulting directory.

1. Build the conda package for your OS and upload it (since I'm running OS X, this will build the OS X package).

   ```bash
   conda build conda-recipe --python 3.4
   binstar upload $HOME/.miniconda/conda-bld/osx-64/iab-0.1.1-np19_0.tar.bz2
   ```

1. Build conda packages for other OSes and upload them (since I'm running OS X, I'll build the Linix packages).

   ```bash
   conda convert --platform linux-32 $HOME/.miniconda/conda-bld/osx-64/iab-0.1.1-np19_0.tar.bz2
   binstar upload linux-32/iab-0.1.1-np19_0.tar.bz2
   conda convert --platform linux-64 $HOME/.miniconda/conda-bld/osx-64/iab-0.1.1-np19_0.tar.bz2
   binstar upload linux-64/iab-0.1.1-np19_0.tar.bz2
   ```

## Part 5: Update the website so the new version of IAB is live.

1. Download the release that was just created on GitHub, and unzip that. Change into the resulting directory.

1. Perform the s3 build and upload, which will required that you have the ``AWS_ACCESS_KEY_ID`` and ``AWS_SECRET_ACCESS_KEY`` environment variables set.

   ```bash
   biab s3build -i book -o s3built
   biab s3upload -i s3built -v latest
   biab s3upload -i s3built -v <version>
   ```
   Note that ``<version>`` should be replaced with the current release version.

1. Optional: If there are changes to the README.md (which is the content that lives at readIAB.org/index.html), convert that to html:

   ```bash
   ipython -c "import biab.util; biab.util.md_to_html('README.md', 'index.html')"
   ```
   Upload the resulting HTML file to the readiab.org S3 bucket.

These steps are pretty clunky now, but will be updated in the near future to be a single step.

## Part 6:

1. Update the version strings to indicate the development release number. This should be ``<version>-dev``, where ``<version>`` is the release that you just published.

1. Update CHANGELOG.md to create a section for changes since the last release.

1. Submit a pull request with these changes. Let the Travis build run, and if the tests pass merge the pull request.

1. Tweet about the release pointing to the release notes on GitHub.

1. :beers:

## Acknowledgements

These notes were developed based on [scikit-bio's RELEASE.md](https://github.com/biocore/scikit-bio/blob/master/RELEASE.md).
