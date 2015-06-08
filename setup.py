#! /usr/bin/env python

##############################################################################
##  DendroBites
##
##  Copyright 2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Package setup and installation.
"""

import sys
import os

###############################################################################
# Identification

from dendrobites import __version__, revision_description, description
sys.stderr.write("-setup.py: {}\n".format(description()))

###############################################################################
# setuptools/distutils/etc. import and configuration

try:
    import ez_setup
    try:
        ez_setup_path = " ('" + os.path.abspath(ez_setup.__file__) + "')"
    except OSError:
        ez_setup_path = ""
    sys.stderr.write("-setup.py: using ez_setup{}\n".format(ez_setup_path))
    ez_setup.use_setuptools()
    import setuptools
    try:
        setuptools_path = " ('" +  os.path.abspath(setuptools.__file__) + "')"
    except OSError:
        setuptools_path = ""
    sys.stderr.write("-setup.py: using setuptools{}\n".format(setuptools_path))
    from setuptools import setup, find_packages
except ImportError as e:
    sys.stderr.write("-setup.py: using distutils\n")
    from distutils.core import setup
    sys.stderr.write("-setup.py: using canned package list\n")
    PACKAGES = [
            "dendrobites",
            ]
else:
    sys.stderr.write("-setup.py: searching for packages\n")
    PACKAGES = find_packages()
EXTRA_KWARGS = dict(
    install_requires = ['setuptools'],
    include_package_data = True,
    #test_suite = "dendropy.test",
    zip_safe = True,
    )

PACKAGE_DIRS = [p.replace(".", os.path.sep) for p in PACKAGES]
PACKAGE_INFO = [("{p[0]:>40} : {p[1]}".format(p=p)) for p in zip(PACKAGES, PACKAGE_DIRS)]
sys.stderr.write("-setup.py: packages identified:\n{}\n".format("\n".join(PACKAGE_INFO)))
ENTRY_POINTS = {}

###############################################################################
# Script paths

SCRIPT_SUBPATHS = [
    #['applications', 'sumtrees', 'sumtrees.py'],
]
SCRIPTS = [os.path.join(*i) for i in SCRIPT_SUBPATHS]
sys.stderr.write("\n-setup.py: scripts identified: {}\n".format(", ".join(SCRIPTS)))

###############################################################################
# setuptools/distuils command extensions

try:
    from setuptools import Command
except ImportError:
    sys.stderr.write("-setup.py: setuptools. Command could not be imported: setuptools extensions not available\n")
else:
    sys.stderr.write("-setup.py: setuptools command extensions are available\n")
    command_hook = "distutils.commands"
    ENTRY_POINTS[command_hook] = []


###############################################################################
# Main setup

### compose long description ###
long_description = open('README.md').read()
long_description = long_description.replace("DendroBites-0.x.x", "DendroBites-{}".format(__version__))

revision_text = revision_description()
long_description = long_description + ("""\

Current Release
===============

The current release of DendroBites is version {}{}.

""".format(__version__, revision_text))

setup(name='DendroBites',
      version=__version__,
      author='Jeet Sukumaran and Mark T. Holder',
      author_email='jeetsukumaran@gmail.com and mtholder@ku.edu',
      #url='http://packages.python.org/DendroPy/',
      description="A Python library demonstrating how to use the DendroPy library for phylogenetics.",
      license='BSD',
      packages=PACKAGES,
      package_dir=dict(zip(PACKAGES, PACKAGE_DIRS)),
      # not needed?
      # package_data={
      #     # "dendropy.utility" : ["libexec/*"],
      #     },
      scripts = SCRIPTS,
      long_description=long_description,
      entry_points = ENTRY_POINTS,
      classifiers = [
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.1",
            "Programming Language :: Python :: 3.2",
            "Programming Language :: Python :: 3.3",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            ],
      keywords='phylogenetics phylogeny phylogenies phylogeography evolution evolutionary biology systematics coalescent population genetics phyloinformatics bioinformatics',
      **EXTRA_KWARGS
      )
