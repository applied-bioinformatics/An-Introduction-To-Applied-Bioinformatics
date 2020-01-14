#!/usr/bin/env python

# -----------------------------------------------------------------------------
# This work is licensed under the Creative Commons
# Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
# copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/.
# -----------------------------------------------------------------------------

__version__ = '0.1.4.dev0'

from setuptools import find_packages, setup
import sys

# Check that Python version is 3, since build will
# currently complete in Python 2, but the notebooks
# won't work.
python_version = sys.version_info.major
if python_version != 3:
  sys.exit("IAB can only be used with Python 3. "
           "You are currently running Python %d." % python_version)

classes = """
    Development Status :: 1 - Planning
    Framework :: IPython
    Intended Audience :: Developers
    Intended Audience :: Education
    Intended Audience :: Science/Research
    Natural Language :: English
    Operating System :: MacOS :: MacOS X
    Operating System :: POSIX
    Operating System :: Unix
    Programming Language :: Python
    Programming Language :: Python :: 3.5
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

description = ("An Introduction To Applied Bioinformatics (IAB): "
               "Interactive lessions in bioinformatics.")

setup(name='An-Introduction-To-Applied-Bioinformatics',
      version=__version__,
      license='CC BY-NC-SA 4.0',
      description=description,
      long_description=description,
      author='Greg Caporaso',
      author_email='gregcaporaso@gmail.com',
      maintainer='Greg Caporaso',
      maintainer_email='gregcaporaso@gmail.com',
      url='http://readIAB.org',
      packages=find_packages(),
      install_requires=['scikit-bio >= 0.5.5, < 0.6.0',
                        'jupyter', 'seaborn',
                        'qiime-default-reference >= 0.1.3, < 0.2.0',
                        'pandas',
                        'markdown2 >= 2.3.0',
                        'tabulate',
                        'networkx == 2.3.0',
                        'ete3',
                        'ipymd >= 0.1.2'],
      classifiers=classifiers)
