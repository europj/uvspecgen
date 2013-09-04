#!/usr/bin/env python

from distutils.core import setup
from uvspec import version as vrsn

version = vrsn.get_version()

with open('README.md') as f:
    long_description = f.read()

setup(name='uvspecgen',
      version=version,
      py_modules=['lib.argparse'],
      packages=['uvspec'],
      scripts=['scripts/uvspecgen'],
      description='Gaussian UV-Vis spectrum generation',
      maintainer='Li Research Group',
      maintainer_email='gaussiantoolkit@gmail.com',
      url='http://github.com/gaussiantoolkit/uvspecgen/',
      license='http://www.gnu.org/licenses/gpl.html',
      long_description=long_description,
     )
