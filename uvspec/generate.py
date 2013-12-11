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
from uvspec.spectrum import AbsorptionSpectrum, generate_outfile_name, \
        join_spectra
from uvspec.ui import CommandLineInput 


def main():
    options = CommandLineInput(defaults)
        
    logfile = options.logfile[0]
    outfile_name = generate_outfile_name(options.outfile, logfile)
 
    if options.join:
        spectrum = join_spectra(options.logfile, options.parameters)
    else:
        spectrum = AbsorptionSpectrum(logfile, options.parameters)
    
    spectrum.create_outfile(outfile_name, options.output, options.nometa)
    
    if options.plot:
        spectrum.plot_spectrum()
    
    print ' Spectrum generation complete: output written to %s' % outfile_name
    return

