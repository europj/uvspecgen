uvspecgen
=========

Generate UV-Vis absorption spectra from TDHF/TDDFT data in ADF, GAMESS,
Gaussian, and Jaguar output files.


Overview
--------
The discrete spectrum obtained from an electronic structure TDHF/TDDFT
calculation is not the most intuitive method for visualizing a simulated
UV-Vis absorption spectrum.  Experimentally, such spectra have broad peaks
and are described by their line shape.  This program extracts the excited
state energies and oscillator strengths from the output files of TDHF/TDDFT 
calculations performed by several widely-used electronic structure programs
and generates a line shape function by summing together Gaussian functions
fit to each peak.

This package includes the `uvspecgen` script for quickly processing electronic
structure TDHF/TDDFT output files to produce a file containing the discrete
spectrum and line shape function for quick plotting.  This script uses the
`uvspec` module, which can be imported into your own Python scripts for use of
the `AbsorptionSpectrum` class, which parses a TDHF/TDDFT output file and
stores the excited state energies, oscillator strengths, and line shape
function as attributes of the class.


Supported Programs
------------------
This program uses the `cclib`[1] Python library for parsing and interpreting
the results of computational chemistry packages.  It currently supports parsing
the results of TDHF/TDDFT calculations for the following electronic structure
programs:
 * ADF
 * GAMESS
 * Gaussian03
 * Gaussian09
 * Jaguar


Installation
------------
Installations are best performed using the `setuptools` Python package via
the included `setup.py` file. To learn more about custom installations, visit
the Installing Python Modules documentation. Standard installations are
described below.

#### Install as Root ####
If you have root permissions on the target system, run the following commands
at the command prompt:

tar xzvf uvspecgen-$VERSION.tar.gz
cd uvspecgen-$VERSION
python setup.py build
sudo python setup.py install

#### Install a Local Copy (UNIX only) ####
For users on UNIX systems that do not have root privileges, a local
installation can be performed.  By default, the program and its modules will
be installed in your home directory at ~/bin/uvspecgen/.  This behavior can
be changed by editing the setup.cfg file located in the config directory.
The installation should then be performed as follows:

tar xvzf uvspecgen-$VERSION.tar.gz
cd uvspecgen-$VERSION
cp config/setup.cfg .
python setup.py install

The source files uvspecgen-$VERSION/ and uvspecgen-$VERSION.tar.gz can be
deleted after installation.


References
----------
 [1] N. M. O'Boyle, A. L. Tenderholt, K. M. Langner, cclib: a library for
     package-independent computational chemistry algorithms, J. Comp. Chem.
     29 (5), pp. 839-845, 2008.
