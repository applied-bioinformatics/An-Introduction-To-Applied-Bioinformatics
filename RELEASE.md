# Published a release of *An Introduction to Applied Bioinformatics*

## Part 1: Prepare the release.

1. Ensure that the Travis build is passing on the master branch.

2. Update version strings. As of this writing, these need to be updated in ``setup.py``, ``iab/__init__.py`` and ``book/getting-started/reading-iab.md``, but it's a good idea to do a global find on the code base for the last release version number to ensure that you're not missing anything. (In the future, the version number will only be specified in one place to simplify this process.)

3. Update [``CHANGELOG.md``](https://github.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/blob/master/CHANGELOG.md) to include descriptions of the changes since the last release.

4. Submit a pull request with these changes and let Travis run.

## Part 2: Build and test the html and IPython Notebook content.

1. Install [build-iab](https://github.com/gregcaporaso/build-iab) if you don't already have it installed.

1. Build the html and IPython Notebooks with the following commands. Poke through these files to ensure that they look right, and that the *Edit* links work correctly. (You may need to adjust the input directory being passed to ``biab``.)

    ```bash
    biab html -i An-Introduction-To-Applied-Bioinformatics/book -o html
biab notebook -i An-Introduction-To-Applied-Bioinformatics/book -o ipynb
    ```

## Part 3: Tag the release.

From the [An Introduction to Applied Bioinformatics GitHub page](https://github.com/gregcaporaso/An-Introduction-to-Applied-Bioinformatics), click on the releases tab and draft a new release. Use the version number for the tag name and create the tag against master. Fill in a release title that is consistent with previous release titles and add a summary of the release (linking to ``CHANGELOG.md`` is a good idea). This release summary will be the primary information that we point users to when we announce the release.


## Part 4: Update the website so the new version of IAB is live.

1. Download the release that was just created on GitHub, and unzip that. Change into the resulting directory.

1. Perform the s3 build, which will required that you have the ``AWS_ACCESS_KEY_ID`` and ``AWS_SECRET_ACCESS_KEY`` environment variables set:

    ```bash
    biab s3 -i book -v <version>
    ```
    Note that ``<version>`` should be replaced with the current release version.

## Part 5:

1. Update the version strings to indicate the development release.

1. Update CHANGELOG.md to create a section for changes since the last release.

1. Submit a pull request with these changes. Let the Travis build run, and if the tests pass merge the pull request.

1. Tweet about the release pointing to the release notes on GitHub.

1. :beers:

## Acknowledgements

These notes were developed while consulting the [scikit-bio release notes](https://github.com/biocore/scikit-bio/blob/master/RELEASE.md).
