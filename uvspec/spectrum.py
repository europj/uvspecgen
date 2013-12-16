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

"""Gaussian09 TDHF/TDDFT UV-Vis absorption spectrum class.

This module contains the AbsorptionSpectrum class for parsing excited state
energies and oscillator strengths from Gaussian09 TDHF/TDDFT log files.  It
then fits the extracted 'stick' spectrum with a Gaussian line shape function
to generate a UV-Vis spectrum.

This module can be imported into your own Python programs for access to the
AbsorptionSpectrum class, or may be used as a stand-alone program for
generating plottable output files using the uvspecgen script.

"""
import datetime
import math
import os
import os.path
import sys

from uvspec.logfile import Logfile


class AbsorptionSpectrum(object):
    """Gaussian UV-Vis spectrum object with 'stick' spectrum and line shape.

    The object requires a logfile name and fit parameters for the line
    shape function to be passed in.  Methods for extracting the 'stick'
    spectrum from the logfile and for generating the line shape function
    are provided.  An output method prints the results (with optional
    user control) to a file with the extension '.spec.txt'.

    A dictionary of parameters containing the following entries must be
    passed in: 
    
        grid
            grid spacing for the energy axis of the line shape function
        range
            applied to the lowest and highest excited state energies to set
            the start and end points of the grid generated for the line shape
            function
        sigma
            broadening parameter for the Gaussian line shape function;
            applies to each Gaussian fit to each excited state
        shift
            shift the energy scale of the line shape function 

    This class provides the following attributes that define the spectrum:
    
        excited_state_energy
            list of excited state energies in eV extracted from a Gaussian09
            log file
        excited_state_wavelength
            list of excited state wavelengths in nm extracted from a
            Gaussian09 log file
        oscillator_strength
            list of oscillator strengths in a.u. for each excitation extracted
            from a Gaussian09 log file
        energy
            list of equally spaced energy grid points in eV for plotting the
            line shape function
        wavelength
            list of equally spaced wavelength grid points in nm for plotting
            the line shape function
        absorbance
            list of absorbance values in a.u. for plotting the line shape
            function

    This class provides the following methods for generating the spectrum by
    extracting the excited state data from the logfile and fitting the 'stick'
    spectrum to generate the line shape function.  Methods are also provided
    for writing the absorption spectrum data to a file and for plotting the
    absorption spectrum:

        generate()
            Runs the extractor method that parses the excited state energies
            and oscillator strengths from a TDHF/TDDFT log file.  The energies
            and oscillator strengths are stored in lists as attributes of the
            class.  Also runs the method for creating the line shape function
            by forming the sum of Gaussian functions fit to each 'stick'.  The
            fitting parameters can be modified by passing a dictionary of fit
            parameters to the instance call of the class.  The energy scale and
            intensities are stored in lists as attributes of the class.
    
        write()
            The method for generating the ouput file.  By default, the output
            file contains the line shape function, the 'stick' spectrum, and
            the fit metadata.  The output can be modified using the
            command-line input options.

        plot()
            This method uses the matplotlib to visualize the line shape
            function.  It catches the ImportError in the event the matplotlib
            plotting module is not available.

    
    """
    def __init__(self, logfile, params, outfile=None):
        """Initialize the AbsorptionSpectrum object.
        
        The object is initialized with a logfile name, parameters for the
        Gaussian line shape function, and metadata to be printed in the output
        file.  Lists are initialized to store the 'stick' and line shape data,
        but are not populated until the appropriate methods are called on the
        object.

        """
        self.logfile = logfile
        # If outfile is given, use it, otherwise, generate it from the logfile
        outfile_name = outfile if outfile else logfile
        self.outfile = self._get_outfile_name(outfile_name)
        self.time = get_time() 

        try:
            self.grid = params['grid']
            self.sigma = params['sigma']
            self.shift = params['shift']
            self.plot_range = params['range']
        except KeyError:
            print (' [ERROR] incomplete parameter definitions in',
                   ' AbsorptionSpectrum')

        # Metadata to be printed in the ouput file
        self.metadata_tags = ['Logfile:',
                              'Sigma:',
                              'Grid:',
                              'Shift:',
                              'Range:',
                              'Created:']
        
        self.metadata = [self.logfile,
                         self.sigma,
                         self.grid,
                         self.shift,
                         self.plot_range,
                         self.time]
        
        # List attributes to store the extracted excited state
        # and line shape data
        self.excited_state_energy = [] 
        self.excited_state_wavelength = []
        self.oscillator_strength = []
        self.absorbance = []                                    
        self.wavelength = []
        self.energy = []                                 

    def generate(self):
        """Generate the absorption spectrum from the logfile.
        
        The excited state energies, wavelengths, and oscillator strengths are
        extracted from the logfile and a Gaussian line shape function is
        generated from the stick spectrum.

        """
        logfile = Logfile(self.logfile)
        self._get_excited_state_data(logfile)
        self._generate_spectrum()

    def join(self, logfiles):
        """Join the excited states from logfiles into the current specturm.

        For each logfile input, a Logfile object is created.  The excited
        state energies and oscillator strengths are compared and duplicates
        are removed.  The current AbsorptionSpectrum object is updated with
        the joined exicted states from all given logfiles.

        """
        logfiles.insert(0, self.logfile)
        energy, strength = self._remove_duplicates(logfiles)
        # Need a Logfile obejct to assign energy and strength for use in the
        # _get_excited_state_data() and _generate_spectrum() methods.
        _logfile = Logfile(self.logfile)
        _logfile.excited_state_energy = energy
        _logfile.oscillator_strength = strength
        self._get_excited_state_data(_logfile)
        self._generate_spectrum()
        # Update the metadata
        self.metadata_tags.append('Joined Files:')
        self.metadata.append(', '.join(logfiles))

    def write(self, output='all', nometa=False):
        """Write the output file.
    
        Write the UV-Vis stick spectrum and line shape data to the output file. 
        The optional input parameters 'output' and 'nometa' control the level
        of output written to the file.  Output can be 'all', 'curve', or
        'sticks' and nometa can be either True or False.
    
        """
        # Setup print control flags; by default, everything is printed
        printout = dict(curve=True, sticks=True, meta=True) 
        if output == 'curve':
            printout['sticks'] = False
        if output == 'sticks':
            printout['curve'] = False 
        if nometa == True:
            printout['meta'] = False 
        
        # Determine the number of lines to print
        # The curve is assumed to always be the longest item to print.  The
        # metadata and stick spectrum are compared to determine which is
        # longest if the curve is not being printed.
        if printout['curve']:
            lines_to_print = len(self.energy)
        elif printout['meta']:
            lines_to_print = max(len(self.excited_state_energy),
                                 len(self.metadata_tags))
        else:
            lines_to_print = len(self.excited_state_energy)

        # Create the output file
        outfile = open(self.outfile, 'w')

        for i in range(lines_to_print):
            # The header line and current line are dynamically generated
            # based on what informaiton was requested for printing 
            header_line = []
            current_line = []
            if printout['curve']: 
                if i == 0:
                    header = '%(energy)15s %(wavelength)17s %(intensity)19s' % \
                              {'energy': 'Energy (eV)',
                               'wavelength': 'Wavelength (nm)',
                               'intensity': 'Intensity (au)'}
                    header_line.append(header)
                
                line = '%(energy)15.3F %(wavelength)17.3F %(intensity)19.5F' % \
                        {'energy': self.energy[i],
                         'wavelength': self.wavelength[i],
                         'intensity': self.absorbance[i]}
                current_line.append(line)

            if printout['sticks']: 
                if i == 0:
                    header = '%(state)8s %(energy)15s %(wavelength)17s %(intensity)19s' % \
                              {'state': 'State', 'energy': 'Energy (eV)',
                               'wavelength': 'Wavelength (nm)',
                               'intensity': 'Intensity (au)'}
                    header_line.append(header)
               
                if i < len(self.excited_state_energy):
                    line = '%(state)8i %(energy)15.3F %(wavelength)17.3F %(intensity)19.5F' % \
                            {'state': i+1,
                             'energy': self.excited_state_energy[i],
                             'wavelength': self.excited_state_wavelength[i],
                             'intensity': self.oscillator_strength[i]}
                elif i < len(self.metadata_tags):
                    line = '%(blank)62s' % {'blank': ' '}
                else:
                    line = ''
                current_line.append(line)
            
            if printout['meta'] and i < len(self.metadata_tags):
                line = '%(tag)15s %(value)20s' % \
                        {'tag': self.metadata_tags[i],
                         'value': self.metadata[i]}
                current_line.append(line)

            if i == 0:
                header_line.append('\n')
                outfile.write(' '.join(header_line))
            current_line.append('\n')
            outfile.write(' '.join(current_line))
        
        outfile.close()
        print (' Spectrum generation complete: output written to %s'
               % self.outfile)
    
    def plot(self):
        """Visualize a plot of the line shape function.
        
        Plot the line shape function as absorbance versus energy with
        appropriate axis labels.  The matplotlib module is not distributed
        as part of the Python Standard Library, so perform a check to
        determine if the module is available.
    
        """
        try:
            import matplotlib.pyplot as plot 
            plot.xlabel('Energy (eV)')
            plot.ylabel('Oscillator strength')
            plot.plot(self.energy, self.absorbance, 'k')
            plot.show()
        except ImportError:
            print ' [ERROR] matplotlib is required to plot the spectrum'

    def _get_outfile_name(self, outfile):
        """Generate the output filename.
        
        The output filename is formed using the log filename prefix and
        replacing the .log extension with the .spec.txt extension.
        
        """
        ext = '.spec.txt'
        if outfile.endswith('.log'):
            name = outfile[:-4]
        else:
            name = outfile
        return str(name + ext) 

    def _get_excited_state_data(self, logfile):
        """Given Logfile object, set ex. state attributes with correct units.
        
        The excited_state_energy attribute is set with units of eV, the
        oscillator_strength attribute is set with arbitrary units, and the
        excited_state_wavelengths is set with units of nm.
    
        """
        es_energy = logfile.excited_state_energy
        # Convert from cm-1 to eV
        self.excited_state_energy = convert_units(es_energy, 'cm-1', 'eV')
        self.oscillator_strength = logfile.oscillator_strength
        self.excited_state_wavelength = \
                convert_units(self.excited_state_energy, 'eV', 'nm')

    def _generate_spectrum(self):
        """Fit the sticks with Gaussians to generate the line shape function.
    
        Generate a grid of energy data points separated by 'grid' within the
        range 'plot_range' above and below the largest and smallest excited
        state energies.  Fit each 'stick' with a Gaussian function of width
        'sigma' shifted by 'shift'.  Sum each Gaussian function to give the
        line shape function.  This method can only be executed after
        get_excited_state_data() has ran.
        
        """
        # Range of graph/output is between 'plot_range' less than the smallest
        # excited state energy and 'plot_range' greater than the largest
        # excited state energy 
        max_energy = self.plot_range + max(self.excited_state_energy) 
        min_energy = min(self.excited_state_energy) - self.plot_range
        point = min_energy
        
        # Generate grid of energy points for the absorbance spectrum
        while point <= max_energy:
            self.energy.append(point)
            point += self.grid
    
        # Generate grid of wavelength points for output of absorbance spectrum
        self.wavelength = convert_units(self.energy, 'eV', 'nm')

        # For every grid point in the energy list, compute the absorbance
        # according to the following equation:
        #   Abs(X) = SUM_{S} Int_{S} * EXP[-0.5 * [(X + SFT - ES_{S}) /
        #              SIG ]**2]
        # where the absorbance, Abs, is a sum of Gaussian functions fit to
        # each S excited state and is a function of the energy, X.  Int is
        # the oscillator strength, SFT is the shift, ES is the excited state,
        # and SIG is the sigma broadening constant.
        for point in self.energy:
            gau_fit = 0.0
            for state in range(len(self.excited_state_energy)):
                gau_fit += self.oscillator_strength[state] * \
                            math.e**(-0.5 * ((point + self.shift - \
                            self.excited_state_energy[state])**2)/
                            (self.sigma**2))
            self.absorbance.append(gau_fit)

    def _remove_duplicates(self, logfiles):
        # Return lists of excited state energies and oscillator strengths with
        # all duplicates removed given a list of logfile names.
        _esdata = []
        for logfile in logfiles:
            _logfile = Logfile(logfile)
            _esdata.extend(zip(_logfile.excited_state_energy,
                               _logfile.oscillator_strength))
        _esdata = list(set(_esdata))
        _esdata.sort()
        energy, strength = zip(*_esdata)
        return energy, strength


def convert_units(data, fromunits, tounits):
    """Return list of floats with units converted from fromunits to tounits."""
    _converter = {'eV_to_nm': lambda x: 1239.84/x,
                  'cm-1_to_eV': lambda x: x/8065.6}
    try:
        converted_data = []
        for value in data:
            if abs(value) <= 1E-8:
                converted_value = 0.0
            else:
                converted_value = _converter['%s_to_%s'
                                             % (fromunits, tounits)](value)
            converted_data.append(converted_value)
        return converted_data
    except KeyError:
        print ' [ERROR] Unsupported conversion in convert_units'
        sys.exit(1)


def get_time():
    """Return the date and time of program execution as MM-DD-YYYY @ HH:MM."""
    now = datetime.datetime.now()
    return now.strftime('%m-%d-%Y @ %H:%M') 
