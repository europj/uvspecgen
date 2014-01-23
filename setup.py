#!/usr/bin/env python
#
# Generate UV-Vis spectra from electronic structure TDHF/TDDFT output files. 
# Copyright (C) 2014 Li Research Group (University of Washington) 
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Installation script for the ``uvspec`` package and ``uvspecgen`` script.

This installer also looks for the necessary dependencies including the
``cclib`` computational chemistry log file parsing library and the 
``argparse`` and ``configparse`` modules.  If these modules are not found, it
will perform these installations as well.

"""
from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

from uvspec.version import get_version 


version = get_version()

with open('README.rst') as f:
    long_description = f.read()

setup(
    name='uvspec',
    version = version,
    packages = find_packages(),

    # automatic generation of the ``uvspecgen`` script
    entry_points = {'console_scripts': ['uvspecgen = uvspec.generate:main']},

    # package dependencies
    install_requires = ['argparse', 'cclib>=1.0', 'configparser'],

    # package metadata
    description = 'Gaussian UV-Vis spectrum generation',
    maintainer = 'Li Research Group',
    maintainer_email = 'liresearchgroupuw@gmail.com',
    url = 'http://github.com/gaussiantoolkit/uvspecgen/',
    license = 'http://www.gnu.org/licenses/gpl.html',
    long_description = long_description,
)
