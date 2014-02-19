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

This package includes the ``uvspecgen`` script for quickly processing
electronic structure TDHF/TDDFT output files to produce a file containing
the discrete spectrum and line shape function for plotting.  This script uses
the included ``uvspec`` module, which can be imported into your own Python
scripts for use of the ``AbsorptionSpectrum`` class, which parses a TDHF/TDDFT
output file and stores the excited state energies, oscillator strengths, and
line shape function as attributes of the class.


Supported Programs
------------------
This program uses the ``cclib`` Python library (ref. 1) for parsing and
interpreting the results of computational chemistry packages.  It currently
supports parsing the results of TDHF/TDDFT calculations for the following
electronic structure programs

* ADF
* GAMESS
* Gaussian03
* Gaussian09
* Jaguar


Dependencies
------------
The ``uvspecgen`` script requires the ``argparse`` and ``configparser`` modules,
which are installed using the ``setup.py`` script as described below.  For
plotting spectra, the ``matplotlib`` package is required.  For information on
how to install ``matplotlib``, visit http://www.matplotlib.org.

The ``uvspec`` module uses the ``cclib`` parsing library to extract the
TDHF/TDDFT excited state energies and oscillator strengths, which is installed
when using the ``setup.py`` script as described below.  The ``cclib`` requires
the ``numpy`` package, which is **not** installed using the ``setup.py`` script.
For information on how to install ``numpy``, visit http://www.numpy.org.


Installation
------------
Installations are best performed using the ``setuptools`` Python package via
the included ``setup.py`` file. To learn more about custom installations, visit
the Installing Python Modules documentation. Standard installations are
described below.

**Install as Root**

If you have root permissions on the target system, run the following commands
at the command prompt::

  tar xzvf uvspec-$VERSION.tar.gz
  cd uvspec-$VERSION
  sudo python setup.py install

**Install a Local Copy**

For users that do not have root privileges, a local installation can be
performed.  The installation location is system-specific, but can quickly
determined by running::

  ./setup.py install --help

and look for the text accompanying the ``--user`` option.  More information on
local installations can be found at
http://docs.python.org/2/install/#alternate-installation-the-user-scheme.

The installation should then be performed as follows::

  tar xzf uvspec-$VERSION.tar.gz
  cd uvspec-$VERSION
  python setup.py install --user

The source files ``uvspec-$VERSION/`` and ``uvspec-$VERSION.tar.gz`` can be
deleted after installation.


Using the AbsorptionSpectrum Class
----------------------------------
The ``AbsorptionSpectrum`` class and methods can be used in your own scripts.
An example of its usage is provided below::

    from uvspec.spectrum import AbsorptionSpectrum

    AS = AbsorptionSpectrum()

    AS.excited_state_energy = [1,2,3]
    AS.oscillator_strength = [0.1,0.2,0.3]

    AS.generate()
    AS.write('api-test')
    AS.plot()

Rather than manually assigning the attributes ``excited_state_energy`` and
``oscillator_strength``, you can call the ``extract(logfile)`` method with a
logfile name provided as a string.  This will extract the excited state energies
and oscillator strengths into these attributes automatically.  The ``plot()``
method only needs to be called if you want to generate a ``matplotlib`` plot of
the line shape function.


References
----------
1. N. M. O'Boyle, A. L. Tenderholt, K. M. Langner, cclib: a library for
   package-independent computational chemistry algorithms, J. Comp. Chem.
   29 (5), pp. 839-845, 2008.
