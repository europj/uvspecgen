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

"""UV-Vis spectrum generation from Gaussian09 TDHF/TDDFT log files.

This program uses the uvspec Python module for generating UV-Vis spectra
from Gaussian09 TDHF/TDDFT log files.  Alternatively, the uvspec module can
be imported into your own programs for use of the AbsorptionSpectrum class.

"""
import sys

from uvspec.config import settings
from uvspec.config.settings import ConfigFile
from uvspec.spectrum import AbsorptionSpectrum


def config_fit_parameters():
    config = ConfigFile()
    for parameter,value in settings.parameters.iteritems():
        config.update(parameter, value)
    print ' Fit paramters have been updated'


def reset_fit_parameters():
    config = ConfigFile()
    config.reset()
    print ' Fit parameters have been reset to their original default values'


def generate_spectrum():
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
    if settings.save:
        config_fit_parameters()
    elif settings.reset:
        reset_fit_parameters() 

    if settings.logfile:
        generate_spectrum()
    elif not settings.save and not settings.reset:
        print ' [ERROR] Must specify at least one logfile name'
        sys.exit()
    else:
        sys.exit()
