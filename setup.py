#!/usr/bin/env python
#
# Generate UV-Vis spectra from Gaussian09 TDHF/TDDFT log files. 
# Copyright (C) 2013  Gaussian Toolkit
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

"""Installation script for the uvspec package and uvspecgen script.

This installer also determines if the user has the cclib computational
chemistry log file parsing library installed.  If not, it will perform that
installation as well.

"""
from setuptools import setup
from uvspec import get_version 


version = get_version()

with open('README.md') as f:
    long_description = f.read()


setup(name='uvspec',
      version=version,
      packages=['uvspec', 'uvspec.config', 'uvspec.lib'],
      scripts=['uvspec/scripts/uvspecgen'],
      package_data={'uvspec': ['tests/*']},
      install_requires=['cclib>=1.0'],
      description='Gaussian UV-Vis spectrum generation',
      maintainer='Li Research Group',
      maintainer_email='gaussiantoolkit@gmail.com',
      url='http://github.com/gaussiantoolkit/uvspecgen/',
      license='http://www.gnu.org/licenses/gpl.html',
      long_description=long_description,
     )
