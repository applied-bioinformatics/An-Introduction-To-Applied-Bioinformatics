#!/usr/bin/env python

# -----------------------------------------------------------------------------
# This work is licensed under the Creative Commons
# Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
# copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/.
# -----------------------------------------------------------------------------

__version__ = '0.0.0-dev'

from setuptools import find_packages, setup

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
    Programming Language :: Python :: 2.7
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

description = "An Introduction To Applied Bioinformatics"

setup(name='An-Introduction-To-Applied-Bioinformatics',
      version=__version__,
      license='CC BY-NC-SA 4.0',
      description=description,
      long_description=description,
      author='Greg Caporaso',
      author_email='gregcaporaso@gmail.com',
      maintainer='Greg Caporaso',
      maintainer_email='gregcaporaso@gmail.com',
      url='https://github.com/gregcaporaso/'
          'An-Introduction-To-Applied-Bioinformatics',
      packages=find_packages(),
      install_requires=['scikit-bio == 0.0.0-dev', 'ipython', 'tornado',
                        'pyzmq', 'jinja2', 'nose', 'runipy', 'biom-format',
                        'pyqi', 'pandas'],
      dependency_links=[
          'https://github.com/biocore/scikit-bio/archive/master.zip#egg='
          'scikit-bio-0.0.0-dev'
      ],
      classifiers=classifiers)
