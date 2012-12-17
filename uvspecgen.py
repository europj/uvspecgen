#!/Library/Frameworks/Python.framework/Versions/Current/bin/python

'''
#!/usr/bin/python
'''

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
#
#                                   UVSPECGEN
#
# Copyright (C) 2012; Li Research Group; University of Washington; Seattle, WA
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

'''
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*                               FILE DESCRIPTION
*
* Filename      : uvspecgen.py
* Version       : 1.0.0
* Programmer(s) : JWM
* Updated on    : 12.17.12
* Created on    : 07.11.12
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
'''

__version__ = '1.0.0'

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
#
# PROGRAM DEFAULTS 
#
# The fit parameters can be changed using the command-line options.  However,
# if a set of values is to be used multiple times, users will save time by
# directly modifying the system defaults below.
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
defaults = dict(
    grid  = 0.02,
    range = 2.50,
    sigma = 0.12,
    shift = 0.00)


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
#
# MODULES
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

import argparse as arp
from math import e


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
#
# CLASSES
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

class CommandLineInput():
    """CommandLineInput()

    An instance of this class parses command-line input and provides program
    help documentation.  All needed user-input is handled with this class,
    including input/output file names and fit parameters.

    """
    def __init__(self, version=__version__, default=defaults):
        self.header = '%(prog)s' + ' ' + version
        self.input = self._parse_command_line_input(self.header, default)
        self.logfile = self.input.logfile
        self.outfile = self.input.outfile
        self.plot = self.input.plot
        self.fit_parameters = dict(
            grid = self.input.grid,
            range = self.input.range,
            sigma = self.input.sigma,
            shift = self.input.shift)


    def _parse_command_line_input(self, header, default):
        """_parse_command_line_input(header : string, default : dictionary) ->
        [list]

        Create an ArgumentParser object and define program options.  The
        default dictionary contains the default values for several parameters
        required to perform the Gaussian fit.  The dictionary is defined at the
        top of this file so that users may easily set the default values
        themself.  Additionally, the dictionary allows a single object to be
        passed into the CommandLineInput class and the herein described
        function.

        """
        parser = arp.ArgumentParser(
            formatter_class = arp.ArgumentDefaultsHelpFormatter,
            description = 'Generate UV-vis spectrom from TDHF/TDDFT data',
            epilog = 'Report bugs to <gaussiantoolkit@gmail.com>')
        
        # Required positional argument, the input (.log) file
        parser.add_argument(
            'logfile',
            type = arp.FileType('r'))
        
        # Optional positional argument, the output file name
        parser.add_argument(
            'outfile',
            nargs = '?')
        
        # Optional arguments, to control the parameters of the Gaussian fit
        parser.add_argument(
            '-g',
            '--grid',
            default = default['grid'],
            type = float,
            help = 'set grid spacing value for the Gaussian curve')
        parser.add_argument(
            '-r',
            '--range',
            default = default['range'],
            type = float,
            help = 'set the spacing below and above the smallest and largest\
            excited state energy for plotting')
        parser.add_argument(
            '-s',
            '--sigma',
            default = default['sigma'],
            type = float,
            help = 'set value for Gaussian broadening constant sigma')
        parser.add_argument(
            '--shift',
            default = default['shift'],
            type = float,
            help = 'set shift for the starting point of the Gaussian curve')

        # Optional flag to plot the spectrum using matplotlib, if available
        parser.add_argument(
            '-p',
            '--plot',
            action = 'store_true',
            help = 'display a plot of the absorbance spectrum if matplotlib is\
                available')
        
        # Print program version
        parser.add_argument(
            '--version',
            action = 'version',
            version = header)
        args = parser.parse_args()
        return args

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
#
# LOCAL FUNCTIONS
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def generate_outfile(logfile):
    """generate_outfile(logfile : string) -> [string]
    
    Removes the log extension from the input file (of the form filename.log)
    and appends the spec extension for the output file

    """
    outfile = logfile[:-3]+'spec'
    return outfile


def get_excited_states(logfile, excited_state, intensity):
    """get_excited_states(logfile : string, excited_state : list, intensity :
        list)
    
    Reads each line of the logfile looking for the 'Excited State' keyword at
    the beginning of the line.  Extracts the excited state energy and
    oscillator strength and stores the values in the appropriate list.

    """
    for line in logfile:                                
        if line.startswith(' Excited State '):      
            words = line.split(' ')
            words = delspaces(words)
            excited_state += [float(words[4])]
            intensity += [float(words[8][2:])]
            

def generate_spectrum(plot_range, grid_spacing, excited_state, intensity,
                      energy, absorbance):
    """generate_spectrum(plot_range : string, grid_spacing : string,
    excited_state : list, intensity : list, energy : list, absorbance : list)

    Given the lists excited_state and intensity, this function generates a
    grid of energy data points separated by grid_spacing between the range
    plot_range above and below the largest and smallest energy excited states.
    The excited state peaks are then fit with Gaussian functions, which are
    summed together to give the final absorbance.  The lists energy and
    absorbance are populated with the data from the fit.
    
    """
    # Range of graph/output is between plot_range less than the smallest
    # excited state energy and plot_range greater than the largest excited
    # state energy 
    max_energy = plot_range + max(excited_state) 
    min_energy = min(excited_state) - plot_range
    point = min_energy
    
    # Generate grid of energy points for the absorbance spectrum
    while point <= max_energy:
        energy.append(point)
        point += grid_spacing

    # For every grid point in the energy list, compute the absorbance according
    # to the following equation:
    #   Abs(X) = SUM_{S} Int_{S} * EXP[-0.5 * [(X + SFT - ES_{S}) / SIG ]**2]
    # where the absorbance, Abs, is a sum of Gaussian functions fit to each S
    # excited state and is a function of the energy, X.  Int is the intensity,
    # SFT is the shift, ES is the excited state, and SIG is the sigma
    # broadening constant.
    for point in energy:
        absorb_fit = 0.0
        for state in range(len(excited_state)):
            absorb_fit += intensity[state]*e**(-0.5*((point + shift - 
                                             excited_state[state])**2)/(sigma**2))
        absorbance.append(absorb_fit)
        
    
def create_outfile(outfile_name, grid_spacing, sigma, shift, excited_state,
                   intensity, energy, absorbance):
    """create_outfile(outfile_name : string, grid_spacing : string, sigma :
        string, shift : string, excited_state : list, intensity : list, 
        energy : list, absorbance : list)

    Writes the contents of the excited_state, intensity, energy, and absorbance
    lists to the output file named outfile_name.

    """
    outfile = open(outfile_name, 'w')
   
    # Details of the fit, including the grid spacing, sigma broadening
    # parameter, and shift magnitude are printed at the top of the ouput file.
    outfile.write('Grid Spacing = '.rjust(17))
    outfile.write(str(grid_spacing).rjust(5))
    outfile.write('\n')
    outfile.write('Sigma = '.rjust(17))
    outfile.write(str(sigma).rjust(5))
    outfile.write('\n')
    outfile.write('Shift = '.rjust(17))
    outfile.write(str(shift).rjust(5))
    outfile.write('\n\n')

    # The data from the logfile (excited state number, energy, and oscillator
    # strength are printed at the top of the file.
    outfile.write('State'.rjust(7))
    outfile.write('Energy (eV)'.rjust(16))
    outfile.write('Intensity (au)\n'.rjust(21))
    for state in range(len(excited_state)):
        outfile.write(str(state+1).rjust(7))
        outfile.write(str(excited_state[state]).rjust(16))
        outfile.write(str(intensity[state]).rjust(20))
        outfile.write('\n')
    outfile.write('\n')
    
    # The fitted Gaussian curve is printed in the second half of the outfile.
    outfile.write('Gaussian Fit\n'.rjust(15))
    outfile.write('Energy (eV)'.rjust(25))
    outfile.write('Intensity (au)\n'.rjust(41))
    for point in range(len(energy)):
        outfile.write(repr(energy[point]).rjust(25))
        outfile.write(repr(absorbance[point]).rjust(40))
        outfile.write('\n')
    outfile.close()


def plot_spectrum(energy, absorbance):
    """plot_spectrum(energy : list, absorbance : list)
    
    Plot the absorbance spectrum as A versus Evec with appropriate axis labels.
    This is done using the matplotlib module, which is not distributed as part
    of the Python Standard Library.  Thus, a check is performed to determine if
    this module is available, only if the plot option is selected; otherwise,
    an error message is returned.

    """
    try:
        import matplotlib.pyplot as plt 
        plt.xlabel("Energy (eV)")
        plt.ylabel("Oscillator strength")
        plt.plot(energy, absorbance, 'k')
        plt.show()
    except ImportError:
        print '  [ERROR] matplotlib.pylot is required to plot the spectrum'


def delspaces(L):
    '''delspaces(L : list) -> [list]
    
    Input a list and recursively remove all items containing only empty
    strings and return the list without empty strings.

    '''
    if len(L)==0:
        return L
    elif L[0]=='':
        return delspaces(L[1:])
    elif L[0]!='':
        return [L[0]] + delspaces(L[1:])
    

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
#
# MAIN PROGRAM
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

# options is an instance of the CommandLineInput class, which provides the
# command-line user interface and help documentation.  All input parameters are
# attributes of this class, including: logfile, outfile, the plot flag and the
# fit parameters stored in the fit_parameters dictionary with the keys grid,
# range, sigma, and shift.
options = CommandLineInput()
logfile = options.logfile

# User may choose to specify name of output file, otherwise, the default name
# will be the logfile name with the .spec extension.
if options.outfile == None:
    outfile_name = generate_outfile(logfile.name)
else:
    outfile_name = options.outfile

# Assign plot flag and fit parameters (attributes of options) to variables
plot = options.plot
grid = options.fit_parameters['grid']
plot_range = options.fit_parameters['range']
sigma = options.fit_parameters['sigma']
shift = options.fit_parameters['shift']

# Create lists for the excited states and oscillator strengths extracted from
# the logfile, as well as for the energies and absorbance values generated from
# the fit.
excited_state = []                                 
intensity = []
energy = [] 
absorbance = []                                    

# Extract the excited state energies and intensities from the logfile
get_excited_states(logfile, excited_state, intensity)

# Fit the extracted excited states with Gaussians and sum to obtain the total
# absorbance spectrum
generate_spectrum(plot_range, grid, excited_state, intensity, energy,
                  absorbance)

# Write the extracted excited states and oscillator strengths, as well as the
# Gaussian-fitted absorbance spectrum to the outfile
create_outfile(outfile_name, grid, sigma, shift, excited_state,
               intensity, energy, absorbance)

# Plot the spectrum, if requested
if plot:
    plot_spectrum(energy, absorbance)

# Indicate successful termination by printing the name of the outfile
print '  Spectrum generation complete: output written to %s' % outfile_name
