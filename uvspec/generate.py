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

"""Generate UV-Vis spectra from electronic structure TDHF/TDDFT output files.

The top-level and main functions for the ``uvspecgen`` script are stored
here.  Functions for updating and resetting the Gaussian fit parameters are
defined here.

"""
import sys

from uvspec.config import settings
from uvspec.config.settings import ConfigFile
from uvspec.spectrum import AbsorptionSpectrum


def _update_fit_parameters():
    # Update the fit parameters in the configuration file with the values
    # provdied by the user at run-time.  If no values are specified, the
    # current default values are maintained.
    config = ConfigFile()
    for parameter, value in settings.parameters.iteritems():
        config.update(parameter, value)
    print ' Fit paramters have been updated'


def _reset_fit_parameters():
    # Reset the fit parameters to their originally installed default values.
    ConfigFile().reset()
    print ' Fit parameters have been reset to their original default values'


def _generate_spectrum():
    # Spectrum generation using the ``AbsorptionSpectrum`` class and methods.
    # If multiple output files are to be joined, an instance of the
    # ``AbsorptionSpectrum`` class is created using the first output
    # filename.  The remianing output filenames are passed to the ``join()``
    # method and processed there.  Creation of the output file and plotting
    # of the spectrum are also processed here.
    spectrum = AbsorptionSpectrum(settings.logfile[0],
                                  settings.parameters,
                                  settings.outfile)
    if settings.join:
        spectrum.join(settings.logfile[1:])
    else:
        spectrum.generate()
    
    spectrum.write(settings.output, settings.nometa)
    
    if settings.plot:
        spectrum.plot()


def main():
    """The core function that drives the ``uvspecgen`` program.

    First handle updates/resetting of the configuration file containing the
    Gaussian fit parameters.  The ``uvspecgen`` program can be run without
    an electronic structure output filename specified solely for the purposes
    of updating/resetting the Gaussian fit parameters in the configuration
    file.

    If the ``--save`` or ``--reset`` flag is not specified, and a logfile is
    not given, the program will terminate with an error message and usage
    instructions.

    """
    if settings.save:
        _update_fit_parameters()
    elif settings.reset:
        _reset_fit_parameters() 

    if settings.logfile:
        _generate_spectrum()
    elif not settings.save and not settings.reset:
        settings.parser.error('Must specify at least one logfile name')
    else:
        sys.exit()
