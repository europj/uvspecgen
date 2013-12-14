#!/usr/bin/env python

from distutils.core import setup
from uvspec import get_version 

version = get_version()

with open('README.md') as f:
    long_description = f.read()

setup(name='uvspecgen',
      version=version,
      packages=['uvspec', 'uvspec.config', 'uvspec.lib'],
      scripts=['uvspec/scripts/uvspecgen'],
      package_data={'uvspec': ['tests/*']},
      description='Gaussian UV-Vis spectrum generation',
      maintainer='Li Research Group',
      maintainer_email='gaussiantoolkit@gmail.com',
      url='http://github.com/gaussiantoolkit/uvspecgen/',
      license='http://www.gnu.org/licenses/gpl.html',
      long_description=long_description,
     )
