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

"""UVSpecGen

This module parses Gaussian09 TDHF/TDDFT log files for excited states and
oscillator strengths.  It then fits this extracted 'stick' spectrum with a
Gaussian line shape function to generate a UV-Vis spectrum.

The module contains the AbsorptionSpectrum class with the following methods 
for extracting the 'stick' spectrum and for generating the line shape:

AbsorptionSpectrum(logfile, params)

    - get_excited_states() -- The extractor method that acts on an instance of
        the AbsorptionSpectrum class to parse the excited state energies and
        oscillator strengths from a Gaussian TDHF/TDDFT log file.  The energies
        and oscillator strengths are stored in lists as attributes of the
        class.
    
    - generate_spectrum() -- The method for creating the line shape function by
        forming the sum of Gaussian functions fit to each 'stick'.  The fitting
        parameters can be modified by passing a dictionary of fit parameters to
        the instance call of the class.  The energy scale and intensities are
        stored in lists as attributes of the class.  This method should be run
        only after get_excited_states() has been executed.
    
    - create_outfile(outfile_name) -- The method for generating the ouput file.
        By default, the output file contains the line shape function, the
        'stick' spectrum, and the fit metadata.  The output can be modified
        using the command-line input options.

    - plot_spectrum() -- This method uses the matplotlib to visualize the line
        shape function.  It catches the ImportError in the event the matplotlib
        plotting module is not available.

The first two methods listed above are executed for each instance of the class
making immediately available the following attributes of the spectrum:

    excited_state_energy
        list of excited state energies in eV extracted from a Gaussian09
        log file
    oscillator_strength
        list of oscillator strengths in a.u. for each excitation extracted
        from a Gaussian09 log file
    energy
        list of equally spaced energy grid points in eV for plotting the line
        shape function
    wavelength
        list of equally spaced wavelength grid points in nm for plotting the
        line shape function
    absorbance
        list of absorbance values in a.u. for plotting the line shape function

A dictionary of parameters must be passed to the AbsorptionSpectrum class and
must include the following entries:

    grid
        grid spacing for the energy axis of the line shape function
    range
        applied to the lowest and highest excited state energies to set the
        start and end points of the grid generated for the line shape function
    sigma
        broadening parameter for the Gaussian line shape function; applies to
        each Gaussian fit to each excited state
    shift
        shift the energy scale of the line shape function 

This module can be imported into your own Python programs for access to the
AbsorptionSpectrum class, or may be used as a stand-alone program for
generating plottable output files using the uvspecgen script.

"""
from math import e
from datetime import datetime

# Conversion factor for eV to nm
EV2NM = 1239.84 

class AbsorptionSpectrum:
    """Gaussian UV-Vis spectrum object with 'stick' spectrum and line shape.

    The object is instantiated with a logfile name and fit parameters for the
    line shape function.  Methods for extracting the 'stick' spectrum from the
    logfile and for generating the line shape function are provided.  An output
    method prints the results (with optional user control) to a file with the
    .spec.txt extension.

    """
    def __init__(self, logfile, params):
        """Initialize the AbsorptionSpectrum object.
        
        The object is initialized with a logfile name, parameters for the
        Gaussian line shape function, and metadata to be printed in the output
        file.  Lists are initialized to store the 'stick' and line shape data,
        but are not populated until the appropriate methods are called on the
        object.

        """
        self.logfile_name = logfile
        self.grid = params['grid']
        self.sigma = params['sigma']
        self.shift = params['shift']
        self.plot_range = params['range']
        self.time = get_time() 
        
        # Metadata to be printed in the ouput file
        self.metadata_tags = ['Logfile:', 'Sigma:', 'Grid:', 'Shift:',
                               'Range:', 'Created:']
        self.metadata = [self.logfile_name, self.sigma, self.grid, self.shift,
                          self.plot_range, self.time]
        
        # List attributes to store the extracted excited state and line shape
        # data
        self.excited_state_energy = [] 
        self.excited_state_wavelength = []
        self.oscillator_strength = []
        self.absorbance = []                                    
        self.wavelength = []
        self.energy = []                                 

        # Extract the excited states and generate the line shape when the class
        # is instantiated
        self.get_excited_states()
        self.generate_spectrum()

    def get_excited_states(self):
        """Read in the excited states and oscillator strengths from log file.
        
        Read each line of the logfile looking for the 'Excited State'
        keyword at the beginning of the line.  Extract the excited state
        energy and oscillator strength into the appropriate list.  Convert eV
        to nm for the excited state wavelengths.
    
        """
        with open(self.logfile_name) as logfile:
            for line in logfile:                                
                if line.startswith(' Excited State '):      
                    words = line.split(' ')
                    words = delspaces(words)
                    self.excited_state_energy += [float(words[4])]
                    self.oscillator_strength += [float(words[8][2:])]
        self.excited_state_wavelength = \
                convert_units(self.excited_state_energy, EV2NM, True)
    
    def generate_spectrum(self):
        """Fit the sticks with Gaussians to generate the line shape function.
    
        Generate a grid of energy data points separated by 'grid' within the
        range 'plot_range' above and below the largest and smallest excited
        state energies.  Fit each 'stick' with a Gaussian function of width
        'sigma' shifted by 'shift'.  Sum each Gaussian function to give the
        line shape function.  This method can only be executed after
        get_excited_states() has ran.
        
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
        self.wavelength = convert_units(self.energy, EV2NM, True)

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
                gau_fit += self.oscillator_strength[state]*e**(-0.5*((point +
                            self.shift - self.excited_state_energy[state])**2)/
                            (self.sigma**2))
            self.absorbance.append(gau_fit)
            
    def create_outfile(self, outfile_name, output, nometa):
        """Write the output file.
    
        Write the UV-Vis stick spectrum and line shape data to the output file. 
    
        """
        # Setup print control flags; by default, everything is printed
        printout = dict(curve=True, sticks=True, meta=True) 
        if output == "curve":
            printout['sticks'] = False
        if output == "sticks":
            printout['curve'] = False 
        if nometa == True:
            printout['meta'] = False 
        
        # Determine number of lines to print
        if printout['curve']:
            lines_to_print = len(self.energy)
        elif printout['meta']:
            lines_to_print = max(len(self.excited_state_energy),
                                 len(self.metadata_tags))
        else:
            lines_to_print = len(self.excited_state_energy)

        # Create the output file
        outfile = open(outfile_name, 'w')

        for i in range(lines_to_print):
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
    
    def plot_spectrum(self):
        """Visualize a plot of the line shape function.
        
        Plot the line shape function as absorbance versus energy with
        appropriate axis labels.  The matplotlib module is not distributed
        as part of the Python Standard Library, so perform a check to
        determine if the module is available.
    
        """
        try:
            import matplotlib.pyplot as plt 
            plt.xlabel("Energy (eV)")
            plt.ylabel("Oscillator strength")
            plt.plot(self.energy, self.absorbance, 'k')
            plt.show()
        except ImportError:
            print ' [ERROR] matplotlib.pyplot is required to plot the spectrum'

def generate_outfile_name(outfile, logfile):
    """Generate the output filename.
    
    The output filename is formed using the log filename prefix and replacing
    the .log extension with the .spec.txt extension.
    
    """
    if outfile == "<logfile>.spec.txt":
        outfile_name = logfile[:-3]+'spec.txt'
    else:
        outfile_name = outfile + ".spec.txt"
    return outfile_name

def combine_spectra(logfiles, params):
    """Merge the stick spectra from multiple log files.

    Given multiple Gaussian TDDFT log files, this routine will extract the
    stick spectrum from each file, combine the excited states into a single
    array removing any duplicate states, and fit the resulting spectrum with a
    Gaussian line shape.

    """
    stick_spectra = []

    for logfile in logfiles:
        stick_spectra.append(AbsorptionSpectrum(logfile, params))
    print stick_spectra

def get_time():
    """Return the date and time of program execution as MM-DD-YYYY @ HH:MM."""
    dt = datetime.now()
    mth = dt.month
    day = dt.day
    yr = dt.year
    hr = dt.hour
    mnt = dt.minute
    time_stamp = '%(month)02i-%(day)02i-%(year)4i @ %(hour)02i:%(minute)02i' %\
            {'month': mth, 'day': day, 'year': yr, 'hour': hr, 'minute': mnt}
    return time_stamp 

def delspaces(L):
    """Remove all items in a list containing only empty strings."""
    if len(L)==0:
        return L
    elif L[0]=='':
        return delspaces(L[1:])
    elif L[0]!='':
        return [L[0]] + delspaces(L[1:])

def convert_units(data, cfactor, inverse=False):
    """Convert the units for the data stored in a list.
    
    The inverse control is used if the conversion involves data that
    is inversely proportional.

    """
    converted_data = []
    for value in data:
        if inverse == False:
            converted_value = value*cfactor
        else:
            converted_value = cfactor/value
        converted_data.append(converted_value)
    return converted_data
