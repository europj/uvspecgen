# Generate UV-Vis spectra from electronic structure TDHF/TDDFT log files. 
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

"""UV-Vis absorption spectrum class for TDHF/TDDFT log files.

This module contains the ``AbsorptionSpectrum`` class for parsing excited
state energies and oscillator strengths from TDHF/TDDFT electronic structure
log files.  It then fits the extracted discrete spectrum with a Gaussian
line shape function to generate a UV-Vis spectrum.

This module can be imported into your own Python programs for access to the
``AbsorptionSpectrum`` class, or may be used as a stand-alone program for
generating plottable output files using the ``uvspecgen`` script.

"""
import datetime
import math

from uvspec.config.settings import error, DEFAULT_FIT_PARAMETERS
from uvspec.logfile import Logfile


class AbsorptionSpectrum(object):
    """UV-Vis absorption spectrum object with discrete and broadened spectra.

    The object has no arguments.  Methods for extracting the discrete
    spectrum from the TDHF/TDDFT log file and for generating the line shape
    function are provided.  An output method prints the results (with optional
    user control) to a file with the extension `.spec.txt`.  If ``matplotlib``
    is available, an optional plotting method can be called to plot the line
    shape function as oscillator strength versus energy.

    This class provides the following attributes that define the spectrum:
    
        excited_state_energy
            list of excited state energies for the discrete spectrum in eV
        excited_state_wavelength
            list of excited state wavelengths for the discrete spectrum in nm
        oscillator_strength
            list of oscillator strengths in a.u. for each excited state in
            the discrete spectrum
        energy
            list of equally spaced energy grid points in eV for plotting
            the line shape function
        wavelength
            list of equally spaced wavelength grid points in nm for plotting
            the line shape function
        absorbance
            list of absorbance values in a.u. for plotting the line shape
            function

    This class provides the following methods for generating the spectrum by
    extracting the excited state data from the TDHF/TDDFT log file and
    fitting the discrete spectrum to generate the line shape function.
    Methods are also provided for writing the absorption spectrum data to a
    file and for plotting the absorption spectrum:

        extract(logfile)
            Runs the extractor method that parses the excited state energies
            and oscillator strengths from a TDHF/TDDFT log file.  The energies
            and oscillator strengths are stored in lists as attributes of the
            class.  A string containing the logfile name must be passed in.

        join(logfiles)
            This method takes a list of logfile names, parses the excited
            state energies and oscillator strengths and combines them with
            the current values of the object.  Any duplicate energy-strength
            pairs are removed.  When finished, the current object contains
            a single, combined UV-Vis spectrum of the logfiles.  Note that the
            ``generate()`` method must still be run to create the line shape.

        generate([params])
            Runs the method for creating the line shape function by forming
            the sum of Gaussian functions fit to each excited state.  The
            fitting parameters can be modified by passing a dictionary of fit
            parameters to the ``params`` argument.  The energy scale and
            intensities are stored in lists as attributes of the class.
    
        write([outfile[, output[, nometa]]])
            The method for generating the ouput file.  By default, the output
            file contains the line shape function, the discrete spectrum, and
            the fit metadata.  The output can be modified by passing the
            keyword arguments ``output`` (takes one of the following strings:
            `all`, `sticks`, or `curve`) and ``nometa`` (takes a boolean).

        plot()
            This method uses ``matplotlib`` to visualize the line shape
            function.  It catches the ``ImportError`` in the event the
            ``matplotlib`` plotting module is not available.
 
    """
    def __init__(self):
        """Initialize the ``AbsorptionSpectrum`` object.
 
        The object is initialized with lists to store the discrete and line
        shape data and for fit metadata, but are not populated until the
        appropriate methods are called on the object.

        """
        self.logfile = None 
        self.outfile = None
        self.joined = None

        self.grid = None
        self.range = None
        self.sigma = None
        self.shift = None

        self.time = get_time() 

        # Metadata to be printed in the ouput file
        self.metadata_tags = ['logfile',
                              'sigma',
                              'grid',
                              'shift',
                              'range',
                              'created']
        
        # List attributes to store the extracted excited state
        # and line shape data
        #
        # Discrete spectrum
        self.excited_state_energy = [] 
        self.excited_state_wavelength = []
        self.oscillator_strength = []
        #
        # Line shape function
        self.absorbance = []                                    
        self.wavelength = []
        self.energy = []                                 

    def extract(self, logfile):
        """Parse excited state energies and oscillator strengths from logfile.
        
        A string containing the logfile name must be passed in.  The excited
        state energies, wavelengths, and oscillator strengths are extracted
        from the logfile as lists and stored in the ``excited_state_energy``,
        ``excited_state_wavelength``, and ``oscillator_strength`` attributes,
        respectively.

        """
        setattr(self, 'logfile', logfile)
        _logfile = Logfile(logfile)
        _logfile.parse()
        self._get_excited_state_data(_logfile)

    def generate(self, params=DEFAULT_FIT_PARAMETERS):
        """Generate the absorption spectrum from the TDHF/TDDFT data.
        
        A Gaussian line shape function is generated from the discrete spectrum
        stored in the ``excited_state_energy`` and ``oscillator_strength``
        attributes.  The optional dictionary of Gaussian fit parameters that
        can be passed in must contain the following keys, where each value
        must be a float:
     
            grid
                spacing of the gridpoints along the energy axis of the line
                shape function
            range
                applied to the lowest and highest excited state energies to
                set the start and end points of the grid generated for the
                line shape function
            sigma
                broadening parameter for the Gaussian line shape function;
                applies to each Gaussian fit to each excited state
            shift
                apply a shift to the energy scale of the line shape function 

        """
        try:
            setattr(self, 'grid', params['grid'])
            setattr(self, 'range', params['range'])
            setattr(self, 'sigma', params['sigma'])
            setattr(self, 'shift', params['shift'])
        except KeyError:
            error('Incomplete parameter definitions in AbsorptionSpectrum')

        if any(self.excited_state_energy) and any(self.oscillator_strength):
            self._generate_spectrum()
        else:
            error('The attributes ``excited_state_energy`` and/or'
                  ' ``oscillator_strength`` are empty')

    def join(self, logfiles):
        """Merge excited states from ``logfiles`` into current specturm.

        For each given logfile, a ``Logfile`` object is created.  The excited
        state energies and oscillator strengths are compared and duplicates
        are removed.  The current ``AbsorptionSpectrum`` object is updated
        with the joined excited states from all given logfiles.

        """
        # The outfile will be named after the first logfile joined, unless
        # specified otherwise in the ``write()`` method
        setattr(self, 'logfile', logfiles[0])

        energy, strength = self._remove_duplicates(logfiles)

        # Need a ``Logfile`` object to assign energy and strength for use
        # in the ``_get_excited_state_data()`` and ``_generate_spectrum()``
        # methods.
        _logfile = Logfile()
        _logfile.excited_state_energy = energy
        _logfile.oscillator_strength = strength

        self._get_excited_state_data(_logfile)

        # Update the metadata
        self.metadata_tags.append('joined files')
        setattr(self, 'joined', ', '.join(logfiles))

    def write(self, outfile=None, output='all', nometa=False):
        """Write the output file.
    
        Write the UV-Vis discrete spectrum, line shape function, and fit
        metadata to the output file.  The ``outfile`` argument is optional so
        long as the ``extract()`` method was run, otherwise, the ``outfile``
        argument is required.  If the ``extract()`` method was run and the
        ``outfile`` argument is not passed in, the logfile name will be used.
        The optional arguments ``output`` and ``nometa`` control the level of
        output written to the file.
        
        The ``output`` argument accepts one of the following strings:
        `all`, `sticks`, or `both` with `all` being the default value.
        
        The ``nometa`` argument accepts a boolean with `False` being the
        default value.

        """
        # Use ``outfile`` if given, otherwise, generate it from ``logfile``
        # permitting the ``extract()`` method was called
        if outfile:
            outfile_name = outfile
        elif self.logfile:
            outfile_name = self.logfile
        else:
            error('An ``outfile`` name must be specified')
        self.outfile = self._get_outfile_name(outfile_name)

        # Setup print control flags; by default, everything is printed
        printout = dict(curve=True, sticks=True, meta=True) 
        if output not in printout and output != 'all':
            print (' [ERROR] Invalid ``output`` parameter in'
                   ' AbsorptionSpectrum.write(), using `all`')
            output = 'all'
        printout['sticks'] = False if output == 'curve' else True
        printout['curve'] = False if output == 'sticks' else True
        printout['meta'] = False if nometa else True
        
        # Determine the number of lines to print
        # The curve is assumed to always be the longest item to print.  The
        # metadata and discrete spectrum are compared to determine which is
        # longest if the curve is not being printed.
        if printout['curve']:
            lines_to_print = len(self.energy)
        elif printout['meta']:
            lines_to_print = max(len(self.excited_state_energy),
                                 len(self.metadata_tags))
        else:
            lines_to_print = len(self.excited_state_energy)

        # Assemble the metadata
        metadata = [self.logfile,
                    self.sigma,
                    self.grid,
                    self.shift,
                    self.range,
                    self.time,
                    self.joined]

        # Create the output file
        outfile = open(self.outfile, 'w')

        for i in range(lines_to_print):
            # The header line and current line are dynamically generated
            # based on what informaiton was requested for printing 
            header_line = []
            current_line = []
            if printout['curve']: 
                if i == 0:
                    header = ('%(energy)15s %(wavelength)17s %(intensity)19s'
                              % {'energy': 'Energy (eV)',
                                 'wavelength': 'Wavelength (nm)',
                                 'intensity': 'Intensity (au)'})
                    header_line.append(header)
                
                line = ('%(energy)15.3F %(wavelength)17.3F %(intensity)19.5F'
                        % {'energy': self.energy[i],
                           'wavelength': self.wavelength[i],
                           'intensity': self.absorbance[i]})
                current_line.append(line)

            if printout['sticks']: 
                if i == 0:
                    header = ('%(state)8s %(energy)15s %(wavelength)17s %(intensity)19s'
                              % {'state': 'State', 'energy': 'Energy (eV)',
                                 'wavelength': 'Wavelength (nm)',
                                 'intensity': 'Intensity (au)'})
                    header_line.append(header)
               
                if i < len(self.excited_state_energy):
                    line = ('%(state)8i %(energy)15.3F %(wavelength)17.3F %(intensity)19.5F'
                            % {'state': i+1,
                               'energy': self.excited_state_energy[i],
                               'wavelength': self.excited_state_wavelength[i],
                               'intensity': self.oscillator_strength[i]})
                elif i < len(self.metadata_tags):
                    line = '%(blank)62s' % {'blank': ' '}
                else:
                    line = ''
                current_line.append(line)
            
            if printout['meta'] and i < len(self.metadata_tags):
                line = '%(tag)15s: %(value)20s' % \
                        {'tag': self.metadata_tags[i].title(),
                         'value': metadata[i]}
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
        appropriate axis labels using the ``matplotlib`` package. 
    
        """
        try:
            import matplotlib.pyplot as plot 
            plot.xlabel('Energy (eV)')
            plot.ylabel('Oscillator strength')
            plot.plot(self.energy, self.absorbance, 'k')
            plot.show()
        except ImportError:
            print (' [ERROR] The ``matplotlib`` package is required'
                   ' to plot the spectrum\n'
                   '         ``matplotlib`` is free to download at'
                   ' http://www.matplotlib.org')

    def _get_outfile_name(self, outfile):
        # Generate the output filename.
        #
        # The output filename is formed using the log filename prefix and
        # replacing the `.log` extension with the `.spec.txt` extension.
        #
        ext = '.spec.txt'
        name = outfile[:-4] if outfile.endswith('.log') else outfile
        return str(name + ext) 

    def _get_excited_state_data(self, logfile):
        # Given Logfile object, set ex. state attributes with correct units.
        #
        # The excited_state_energy attribute is set with units of eV and the
        # oscillator_strength attribute is set with arbitrary units.
        #
        es_energy = logfile.excited_state_energy
        # Convert from cm-1 to eV
        self.excited_state_energy = convert_units(es_energy, 'cm-1', 'eV')
        self.oscillator_strength = logfile.oscillator_strength

    def _generate_spectrum(self):
        # Fit the sticks with Gaussians to generate the line shape function.
        # 
        # Generate a grid of energy data points separated by ``grid`` within
        # the range ``range`` above and below the largest and smallest
        # excited state energies.  Fit each `stick` with a Gaussian function
        # of width ``sigma`` shifted by ``shift``.  Sum each Gaussian function
        # to give the line shape function.  This method can only be executed
        # after ``_get_excited_state_data()`` has run.
        # 
        # Range of graph/output is between 'range' less than the smallest
        # excited state energy and 'range' greater than the largest excited
        # state energy 
        max_energy = self.range + max(self.excited_state_energy) 
        min_energy = min(self.excited_state_energy) - self.range
        point = min_energy
        
        # Generate grid of energy points for the absorbance spectrum
        while point <= max_energy:
            self.energy.append(point)
            point += self.grid
    
        # Generate grid of wavelength points for output of absorbance spectrum
        self.wavelength = convert_units(self.energy, 'eV', 'nm')
        self.excited_state_wavelength = \
                convert_units(self.excited_state_energy, 'eV', 'nm')

        # For every grid point in the energy list, compute the absorbance
        # according to the following equation:
        #
        # Abs(X) = SUM_{S} Int_{S} * EXP[-0.5 * [(X + SFT - ES_{S}) / SIG ]**2]
        #
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
        # Given a list of logfile names, return lists of excited state
        # energies and oscillator strengths with all duplicates removed. 
        _esdata = []
        for logfile in logfiles:
            _logfile = Logfile(logfile)
            _logfile.parse()
            _esdata.extend(zip(_logfile.excited_state_energy,
                               _logfile.oscillator_strength))
        # The Python built-in function ``set()`` removes all duplicate
        # energy-strength pairs from the ``_esdata`` list.
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
        error('Unsupported conversion in ``convert_units``')


def get_time():
    """Return string of the current date and time as `MM-DD-YYYY @ HH:MM`."""
    now = datetime.datetime.now()
    return now.strftime('%m-%d-%Y @ %H:%M') 
